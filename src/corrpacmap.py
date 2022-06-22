
import numpy as np
import xarray as xr
import zarr

from scipy.special import gamma, kv

from pych import read_jacobi_iters
from matern import MaternField
from new_smooth_store import open_smoothdataset
from llcutils import get_pacific
from timer import Timer

from corrcalc import CorrelationCalculator

class CorrelationPacMap(CorrelationCalculator):
    """Slightly different than correlation calculator,
        -> don't automatically subset outer region
        -> compute correlation at single point in 3D
    """

    spot        = {"ix": 90, "iy": 165, "k":25}

    @property
    def spot_xy(self):
        return {k:spot[k] for k in ["ix","iy"]}

    @property
    def spot_xz(self):
        return {k:spot[k] for k in ["ix","iy"]}

    @property
    def zstore_path(self):
        return super().zstore_path.replace("matern-correlation", "matern-correlation-pacmap")

    def __call__(self):

        # Get samples
        ds = self.open_dataset()

        # Compute correlation dataset
        cds = self.calc_correlation(ds["ginv_norm"], ds["ginv_norm_mean"])

        cds = self.expand_dims(cds)
        self.save_results(cds)


    def open_dataset(self):
        ds = open_smoothdataset(data_dir=self.diag_dir,
                                grid_dir=self.run_dir,
                                geometry="llc",
                                k_chunksize=self.work_chunks["k"],
                                iter_stop=self.n_samples+1)

        # Drop unnecessary stuff, keep maskC
        if self.drop_coords:
            ds = ds.reset_coords("maskC")
            ds = ds.reset_coords(drop=True)
            keep = ['maskC', 'smooth3Dnorm001', 'smooth3Dfld001', 'smooth3Dmean001']
            remove = [v for v in list(ds.data_vars) if v not in keep]
            ds = ds.drop(remove)

        # Isolate subregion
        ds = get_pacific(ds)

        # rechunk
        ds = ds.chunk(self.work_chunks)

        # Compute some things
        with xr.set_options(keep_attrs=True):
            ds['smooth3Dnorm001'] = ds['smooth3Dnorm001'].where(ds['maskC'])
        ds['ginv_norm'] = ds['smooth3Dfld001'] * ds['smooth3Dnorm001']
        ds['ginv_norm_mean'] = ds['smooth3Dmean001'] * ds['smooth3Dnorm001']

        if self.persist:
            for key in ["maskC", "smooth3Dnorm001", "smooth3Dmean001", "ginv_norm", "ginv_norm_mean"]:
                ds[key] = ds[key].persist()

        return ds


    def calc_correlation(self, xda, xda_mean, client=None):

        # make container with ideal/expected correlation
        cds = xr.Dataset()

        nz_shift = min(24, self.n_shift)
        shift_z  = np.arange(-nz_shift, nz_shift+1)

        # do this part once
        x_deviation = xda - xda_mean
        spot_deviation   = x_deviation.sel(self.spot).drop(["ix", "iy", "k"]).compute()
        spot_ssr_inv   = (1 / np.sqrt( (spot_deviation**2).sum("sample") ) ).compute()
        if self.persist:
            x_deviation = x_deviation.persist()

        # for each k shift, compute xy correlations
        if client is None:

            for ks in shift_z:

                corrlist = []
                this_k = self.spot["k"] + ks
                tmp = self.calc_xyshifts(xdev=x_deviation.sel(k=this_k),
                                         spot_deviation=spot_deviation,
                                         spot_ssr_inv=spot_ssr_inv)
                corrlist.append(tmp.expand_dims({"k": [this_k]}))
        else:


            xdevs = [x_deviation.sel(k=self.spot["k"]+ks) for ks in shift_z]
            futures = client.map(self.calc_xyshifts, xdevs,
                                 spot_deviation=spot_deviation,
                                 spot_ssr_inv=spot_ssr_inv)

            corrlist = client.gather(futures)

        corrfld = xr.concat(corrlist, dim="k")
        corrfld.attrs = {"description": f"3D correlation field computed about spot",
                         "spot": self.spot}
        return corrfld


    def calc_xyshifts(self, xdev, spot_deviation, spot_ssr_inv):
        """
        Args:
            xdev: 2D field at shifted k
            spot_deviation, spot_ssr_inv: deviation and ssr inverse at single point
            zs

        """

        corr_xy = xr.zeros_like(xdev.sel(sample=0).drop("sample"))

        shift_xy = np.arange(-self.n_shift, self.n_shift+1)
        for yshift in shift_xy:
            for xshift in shift_xy.copy():
                shift = {"ix": self.spot["ix"] + xshift,
                         "iy": self.spot["iy"] + yshift}

                y_deviation = xdev.sel(shift).drop(["ix", "iy", "k"])
                numerator = (spot_deviation * y_deviation).sum("sample")
                y_ssr = np.sqrt( (y_deviation**2).sum("sample") )

                corr = numerator * spot_ssr_inv / y_ssr
                corr_xy.loc[shift] = corr

        return corr_xy
