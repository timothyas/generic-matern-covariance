
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
    n_samples           = 1000
    mean_differentiability = 1/2

    # set dimensions and subregion
    n_shift             = 50
    dimlist             = ("ix", "iy", "k")
    outer               = {"ix": slice( 40, 120), "iy": slice(200,  45)}
    inner               = {"ix": slice( 60, 100), "iy": slice(160,  80)}

    # dask/zarr stuff
    persist             = True
    work_chunks         = {"sample" : 100,
                           "k"      : 50,
                           "iy"     : -1,
                           "ix"     : -1}
    save_chunks         = {"shifty" : 1,
                           "k"      : 1,
                           "iy"     : -1,
                           "ix"     : -1}

    @property
    def run_dir(self):
        return f"/scratch2/tsmith/generic-matern-covariance/sampling/llc90/matern-sample-log10tol{self.log10tol:03d}-3D-C/run.{self.n_range:02d}dx.{self.horizontal_factor:02}xi"


    @property
    def diag_dir(self):
        return self.run_dir + "/smooth-output"


    @property
    def zstore_path(self):
        return f"/scratch2/tsmith/generic-matern-covariance/sampling/llc90/zstores/matern-correlation.log10tol{self.log10tol:03d}.{self.n_range:02d}dx.{self.horizontal_factor:02}xi"


    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            setattr(self, key, val)


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
        ds = ds.reset_coords("maskC")
        ds = ds.reset_coords(drop=True)
        keep = ['maskC', 'smooth3Dnorm001', 'smooth3Dfld001', 'smooth3Dmean001']
        remove = [v for v in list(ds.data_vars) if v not in keep]
        ds = ds.drop(remove)

        # Isolate subregion
        ds = get_pacific(ds)
        ds = ds.sel(self.outer)

        # rechunk
        ds = ds.chunk(self.work_chunks)

        # Compute some things
        with xr.set_options(keep_attrs=True):
            ds['smooth3Dnorm001'] = ds['smooth3Dnorm001'].where(ds['maskC'])
        ds['ginv_norm'] = ds['smooth3Dfld001'] * ds['smooth3Dnorm001']
        ds['ginv_norm_mean'] = ds['smooth3Dmean001'] * ds['smooth3Dnorm001']

        # Get number of iterations for solve
        myiters = read_jacobi_iters(f"{self.run_dir}/STDOUT.0000", which_jacobi="3D")
        myiters = xr.DataArray(myiters[:len(ds.sample)], ds['sample'].coords, ds['sample'].dims,
                               attrs={'description': 'Number of iterations to converge, given tolerance and range'})
        ds['sor_iters'] = myiters

        if self.persist:
            for key in ["maskC", "smooth3Dnorm001", "smooth3Dmean001", "ginv_norm", "ginv_norm_mean"]:
                ds[key] = ds[key].persist()

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


if __name__ == "__main__":

    from dask_jobqueue import SLURMCluster
    from dask.distributed import Client

    cluster = SLURMCluster(log_directory="/scratch2/tsmith/dask-jobqueue-space")
    cluster.adapt(minimum=0, maximum=20)
    client = Client(cluster)

    stdout = "stdout.correlation-timing.log"
    localtime = Timer(filename=stdout)
    walltime = Timer(filename=stdout)

    walltime.start("Starting job")

    for n_range in [5, 10, 15, 20]:
        for log10tol in [-1, -2, -4, -6, -8, -10, -12, -14]:
            localtime.start(f"n_range = {n_range}, log10tol = {log10tol}")
            cc = CorrelationCalculator(n_range=n_range, log10tol=log10tol)
            cc()
            localtime.stop()

    walltime.stop("Total Walltime")
