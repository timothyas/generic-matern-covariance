
import numpy as np
import xarray as xr
import zarr

from scipy.special import gamma, kv

from pych import read_jacobi_iters
from matern import MaternField
from new_smooth_store import open_smoothdataset
from llcutils import get_pacific
from timer import Timer

class CorrelationCalculator():
    """Read in samples, compute correlation for subset of domain, save"""

    n_range             = None
    log10tol            = None
    horizontal_factor   = 1
    n_applications      = 1
    n_samples           = 1000
    isoxy               = False
    baro_rad            = False

    # set dimensions and subregion
    n_shift             = 50
    dimlist             = ("ix", "iy", "k")
    outer               = {"ix": slice( 40, 140), "iy": slice(230, 40)}
    inner               = {"ix": slice( 60, 120), "iy": slice(210, 60)}
    select_outer_subregion = True

    # dask/zarr stuff
    drop_coords         = True
    persist             = True
    load_samples        = False
    work_chunks         = {"sample" : 100,
                           "k"      : 50,
                           "iy"     : -1,
                           "ix"     : -1}
    save_chunks         = {"shifty" : 1,
                           "k"      : 1,
                           "iy"     : -1,
                           "ix"     : -1}

    @property
    def mean_differentiability(self):
        return 2*self.n_applications - 3/2


    @property
    def main_run_dir(self):
        if self.isoxy:
            return f"/scratch2/tsmith/generic-matern-covariance/sampling/llc90/matern-isoxy"
        elif self.baro_rad:
            return f"/scratch2/tsmith/generic-matern-covariance/sampling/llc90/matern-barotropic-radius"

        else:
            return f"/scratch2/tsmith/generic-matern-covariance/sampling/llc90/matern-{self.n_applications:02d}apps"

    @property
    def run_dir(self):
        return f"{self.main_run_dir}/log10tol{self.log10tol:03d}-3D-C/run.{self.n_range:02d}dx.{self.horizontal_factor:02}xi"


    @property
    def diag_dir(self):
        return self.run_dir + "/smooth-output"

    @property
    def main_zstore_path(self):
        if self.isoxy:
            return f"/scratch2/tsmith/generic-matern-covariance/sampling/llc90/zstores/matern-corr-isoxy"
        elif self.baro_rad:
            return f"/scratch2/tsmith/generic-matern-covariance/sampling/llc90/zstores/matern-corr-barotropic-radius"
        else:
            return f"/scratch2/tsmith/generic-matern-covariance/sampling/llc90/zstores/matern-corr-{self.n_applications:02d}apps"

    @property
    def zstore_path(self):
        return f"{self.main_zstore_path}/log10tol{self.log10tol:03d}.{self.n_range:02d}dx.{self.horizontal_factor:02}xi.{self.n_samples:04d}samples"


    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            setattr(self, key, val)

        if self.load_samples and self.persist:
            warnings.warn("CorrelationCalculator.__init__: load_samples and persist set, ignoring persist and will load ginv_norm field into memory")
            self.persist = False


    def __call__(self):

        # Get samples
        ds = self.open_dataset()

        # Compute correlation dataset
        cds = self.calc_correlation(ds["ginv_norm"], ds["ginv_norm_mean"])
        cds["avg_sor_iters"] = ds["sor_iters"].mean("sample")

        cds = self.expand_dims(cds)
        self.save_results(cds)


    def open_dataset(self):
        ds = open_smoothdataset(data_dir=self.diag_dir,
                                grid_dir=self.run_dir,
                                geometry="llc",
                                k_chunksize=self.work_chunks["k"],
                                iter_stop=self.n_samples+1)



        # Drop unnecessary stuff, keep maskC
        ds = self._drop_coords(ds) if self.drop_coords else ds

        # Isolate subregion
        ds = get_pacific(ds)
        if self.select_outer_subregion:
            ds = ds.sel(self.outer)

        # rechunk
        ds = ds.chunk(self.work_chunks)

        # Compute some things
        ds = self.calc_ginv_norm(ds)

        # Get number of iterations for solve
        ds = self.get_sor_iters(ds)

        ds = self.persist_or_load(ds)
        return ds


    def calc_correlation(self, xda, xda_mean):

        # make container with ideal/expected correlation
        cds = xr.Dataset()
        shifty = np.arange(-self.n_shift, self.n_shift+1)
        cds["shifty"] = xr.DataArray(shifty, {"shifty": shifty}, ("shifty",))
        cds["rho_hat"] = np.abs(cds['shifty'])
        cds["ideal_corr"] = self.ideal_correlation(distance=cds["rho_hat"])

        # do this part once
        x_deviation = xda - xda_mean
        selection   = x_deviation.sel(self.inner)
        x_ssr_inv   = (1 / np.sqrt( (selection**2).sum("sample") ) )
        if self.persist:
            x_deviation = x_deviation.persist()
            selection   = selection.persist()
            x_ssr_inv   = x_ssr_inv.persist()

        # shift through each dimension separately
        for dim in self.dimlist:

            corrfld = []
            for s in shifty:
                y_deviation = x_deviation.shift({dim: s}).sel(self.inner)
                numerator = (selection * y_deviation).sum("sample")
                y_ssr = np.sqrt( (y_deviation**2).sum("sample") )

                corr = numerator * x_ssr_inv / y_ssr
                corrfld.append(corr)

            cds[f"corr_{dim}"] = xr.concat(corrfld, dim="shifty")

        cds.attrs = {"description": f"1D correlations computed from {len(xda.sample)} samples"}
        return cds


    def ideal_correlation(self, distance):
        """Compute correlation in ideal space"""

        factor = 2 ** (1-self.mean_differentiability) / gamma(self.mean_differentiability)
        arg = np.sqrt(8*self.mean_differentiability)/self.n_range * distance

        left = arg ** self.mean_differentiability
        where = xr.where if isinstance(distance, xr.DataArray) else np.where
        right = where(distance == 0, 1, kv(self.mean_differentiability, arg))
        result = factor * left * right
        return where(distance == 0 , 1, result)


    def expand_dims(self, ds):
        return ds.expand_dims({'log10tol': [self.log10tol], 'n_range': [self.n_range]})


    def save_results(self, ds):
        store = zarr.NestedDirectoryStore(path=self.zstore_path)
        ds = ds.chunk(self.save_chunks)
        ds.to_zarr(store=store)

    @staticmethod
    def _drop_coords(ds):
        ds = ds.reset_coords("maskC")
        ds = ds.reset_coords(drop=True)
        keep = ['maskC', 'smooth3Dnorm001', 'smooth3Dfld001', 'smooth3Dmean001']
        remove = [v for v in list(ds.data_vars) if v not in keep]
        ds = ds.drop(remove)
        return ds


    @staticmethod
    def calc_ginv_norm(ds):
        with xr.set_options(keep_attrs=True):
            ds['smooth3Dnorm001'] = ds['smooth3Dnorm001'].where(ds['maskC'])
        ds['ginv_norm'] = ds['smooth3Dfld001'] * ds['smooth3Dnorm001']
        ds['ginv_norm_mean'] = ds['smooth3Dmean001'] * ds['smooth3Dnorm001']
        return ds


    def get_sor_iters(self, ds):
        myiters = read_jacobi_iters(f"{self.run_dir}/STDOUT.0000", which_jacobi="3D")

        myiters = np.array([ np.sum(myiters[i:i+self.n_applications]) for i in np.arange(0, self.n_samples*self.n_applications, self.n_applications)])
        myiters = xr.DataArray(myiters[:len(ds.sample)], ds['sample'].coords, ds['sample'].dims,
                               attrs={'description': 'Number of iterations to converge, given tolerance and range'})
        ds['sor_iters'] = myiters
        return ds

    def persist_or_load(self, ds):
        if self.persist:
            for key in ["maskC", "smooth3Dnorm001", "smooth3Dmean001", "ginv_norm", "ginv_norm_mean"]:
                ds[key] = ds[key].persist()

        elif self.load_samples:
            ds["ginv_norm_mean"].load();
            ds["ginv_norm"].load();

        return ds

if __name__ == "__main__":

    # Note to self, this is faster to load in and read, rather than parallelize
    #from dask_jobqueue import SLURMCluster
    #from dask.distributed import Client

    #cluster = SLURMCluster(log_directory="/scratch2/tsmith/dask-jobqueue-space")
    #cluster.adapt(minimum=0, maximum=10)
    #client = Client(cluster)

    stdout = "stdout.correlation-timing.1000samples.log"
    localtime = Timer(filename=stdout)
    walltime = Timer(filename=stdout)

    walltime.start("Starting job")

    n_applications = 1
    for n_range in [5, 10, 15, 20]:
        for log10tol in [-1, -2, -3, -4, -7, -11, -15]:

            localtime.start(f"n_range = {n_range}, log10tol = {log10tol}, n_applications = {n_applications}")
            cc = CorrelationCalculator(n_range=n_range,
                                       log10tol=log10tol,
                                       n_applications=n_applications,
                                       n_samples=1000,
                                       isoxy=False,
                                       load_samples=True,
                                       persist=False)
            cc()
            localtime.stop()

    log10tol = -3
    for n_range in [5, 10, 15, 20]:
        for n_applications in [1, 2, 4, 8]:

            localtime.start(f"n_range = {n_range}, log10tol = {log10tol}, n_applications = {n_applications}")
            cc = CorrelationCalculator(n_range=n_range,
                                       log10tol=log10tol,
                                       n_applications=n_applications,
                                       n_samples=1000,
                                       isoxy=False,
                                       load_samples=True,
                                       persist=False)
            cc()
            localtime.stop()

    walltime.stop("Total Walltime")
