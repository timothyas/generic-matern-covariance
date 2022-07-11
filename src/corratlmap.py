
import numpy as np
import xarray as xr

from new_smooth_store import open_smoothdataset
from llcutils import get_atlantic
from timer import Timer

from corrpacmap import CorrelationPacMap

class CorrelationAtlMap(CorrelationPacMap):
    """Just grab a different region
    """

    name = "shelf"
    spot = {"ix": 65, "iy": 110, "k":0}

    @property
    def main_zstore_path(self):
        return super().main_zstore_path.replace("/matern-pacmap", "/matern-atlmap")

    def open_dataset(self):
        ds = open_smoothdataset(data_dir=self.diag_dir,
                                grid_dir=self.run_dir,
                                geometry='llc',
                                k_chunksize=self.work_chunks["k"],
                                iter_stop=self.n_samples+1)

        # Drop unnecessary stuff, keep maskC
        ds = self._drop_coords(ds) if self.drop_coords else ds

        # Isolate subregion
        ds = get_atlantic(ds)
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


if __name__ == "__main__":

    from dask.distributed import Client

    stdout = "stdout.atlmap-timing.1000samples.log"
    localtime = Timer(filename=stdout)
    walltime = Timer(filename=stdout)

    walltime.start("Starting job")

    localtime.start("shelf")
    cc = CorrelationAtlMap(n_range=1,
                           horizontal_factor=0.01,
                           log10tol=-4,
                           n_samples=1000,
                           baro_rad=True,
                           drop_coords=False,
                           load_samples=True,
                           persist=False)

    ds = cc.open_dataset()
    client = Client()
    cds = cc.calc_correlation(ds.ginv_norm, ds.ginv_norm_mean, client=client)
    client.close()

    cc.save_results(cds.to_dataset(name=f"corr_{cc.name}"))

    localtime.stop()

    walltime.stop("Total Walltime")
