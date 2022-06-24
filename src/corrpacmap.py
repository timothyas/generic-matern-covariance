
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

    name = "equator"
    spot = {"ix": 90, "iy": 165, "k":25}

    drop_coords = False
    select_outer_subregion  = False

    save_chunks = {"k": 1, "iy": -1, "ix": -1}

#    @property
#    def spot_xy(self):
#        return {k:spot[k] for k in ["ix","iy"]}
#
#    @property
#    def spot_xz(self):
#        return {k:spot[k] for k in ["ix","iy"]}

    @property
    def main_zstore_path(self):
        return super().main_zstore_path.replace("matern-correlation", "matern-pacmap")

    @property
    def zstore_path(self):
        return f"{self.main_zstore_path}/{name}.{self.n_range:02d}dx.{self.horizontal_factor:02}xi"


    def __call__(self, client=None):

        # Get samples
        ds = self.open_dataset()

        # Compute correlation dataset
        cds = self.calc_correlation(ds["ginv_norm"], ds["ginv_norm_mean"], client=client)

        cds = self.expand_dims(cds)
        self.save_results(cds)


    def calc_correlation(self, xda, xda_mean, client=None):

        # make container with ideal/expected correlation
        cds = xr.Dataset()

        # do this part once
        x_deviation = xda - xda_mean
        spot_deviation   = x_deviation.sel(self.spot).drop(["ix", "iy", "k"]).compute()
        spot_ssr_inv   = (1 / np.sqrt( (spot_deviation**2).sum("sample") ) ).compute()
        if self.persist:
            x_deviation = x_deviation.persist()

        # for each k shift, compute xy correlations
        all_k  = self.get_shifted_indices("k", xda["k"])
        if client is None:

            for this_k in all_k:

                corrlist = []
                tmp = self.calc_xyshifts(xdev=x_deviation.sel(k=this_k),
                                         spot_deviation=spot_deviation,
                                         spot_ssr_inv=spot_ssr_inv)
                corrlist.append(tmp.expand_dims({"k": [this_k]}))
        else:


            xdevs = [x_deviation.sel(k=this_k) for this_k in all_k]
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

        all_y = self.get_shifted_indices("iy", xdev["iy"])
        all_x = self.get_shifted_indices("ix", xdev["ix"])
        for this_y in all_y:
            for this_x in all_x:
                this_xy = {"ix": this_x, "iy": this_y}

                y_deviation = xdev.sel(this_xy).drop(["ix", "iy", "k"])
                numerator = (spot_deviation * y_deviation).sum("sample")
                y_ssr = np.sqrt( (y_deviation**2).sum("sample") )

                corr = numerator * spot_ssr_inv / y_ssr
                corr_xy.loc[this_xy] = corr

        return corr_xy


    def get_shifted_indices(self, dim, arr):
        """get lower and upper bounds for shifting, make an array"""

        lo = max(int(arr.min()), self.spot[dim] - self.n_shift)
        hi = min(int(arr.max()), self.spot[dim] + self.n_shift)
        return np.arange(lo, hi)


if __name__ == "__main__":

    from dask.distributed import Client

    stdout = "stdout.pacmap-timing.1000samples.log"
    localtime = Timer(filename=stdout)
    walltime = Timer(filename=stdout)

    walltime.start("Starting job")

    equator = {"ix":  90, "iy": 165, "k": 25}
    coast   = {"ix": 130, "iy": 180, "k":  0}

    for name, spot in zip(["equator", "coast"], [equator, coast]):
        cc = CorrelationPacMap(name=name,
                               spot=spot,
                               n_range=20,
                               log10tol=-15,
                               n_samples=1000,
                               drop_coords=False,
                               load_samples=True,
                               persist=False)

        client = Client()
        cc(client=client)
        client.close()
        localtime.stop()

    walltime.stop("Total Walltime")
