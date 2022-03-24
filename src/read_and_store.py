"""
Read smooth3Dfld* eagerly because this is nuts... store to zarr
"""

import os
import numpy as np
import xarray as xr
import zarr

import ecco_v4_py
from xmitgcm import open_mdsdataset

from smooth_store import open_smoothdataset

def eager_smooth_reader(run_dir, fname, n_samples=1000, chunksize=200, less_output=True):
    smoothfld = []
    for skip in range(0,n_samples,chunksize):
        smoothfld.append(ecco_v4_py.read_llc_to_tiles(run_dir, fname,
                                                      nk=50,
                                                      nl=chunksize,
                                                      skip=skip*50,
                                                      filetype='>f4',
                                                      less_output=less_output))
    return np.concatenate(smoothfld, axis=0)

def get_dataset(run_dir):

    ds = open_mdsdataset(run_dir, iters=None, geometry="llc")
    ds = open_smoothdataset(run_dir, xdalike=ds.hFacC, iters=None, geometry="llc", chunks={"face":13})
    ds["smooth3Dfld001"].values = eager_smooth_reader(run_dir, fname="smooth3Dfld001.data")
    ds["smooth3Dfld001"] = ds["smooth3Dfld001"].chunk({'sample':1, 'k':1, 'face':None, 'j':None, 'i':None})

    del ds["ginv_norm"]
    return ds


if __name__=="__main__":

    main_dir = "/scratch2/tsmith/generic-matern-covariance/sampling/llc90/sample-3D-C"

    for n_range in [5, 10, 15, 20]:
        for xi in [.5, 1, 2]:
            run_dir = f"{main_dir}/run.{n_range:02d}dx.{xi:02}xi"
            ds = get_dataset(run_dir)

            path = f"{main_dir}/zstores/matern.{n_range:02d}dx.{xi:02}xi.zarr"
            store = zarr.NestedDirectoryStore(path=path)
            ds.to_zarr(store=store, consolidated=True)
