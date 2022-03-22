
import numpy as np
import xarray as xr
from xmitgcm import open_mdsdataset
import pych.pigmachine as pm


if __name__ == "__main__":

    main_run='/scratch2/tsmith/generic-matern-covariance/sampling/llc90'
    ds = open_mdsdataset(f"{main_run}/grid", iters=None, geometry='llc')

    # --- directories and dicts
    dirs = {'main_run'      : main_run,
            'netcdf'        : main_run+'/ncfiles'}

    dsim = {'machine'       : 'sverdrup',
            'n_procs'       : 96,
            'exe_file'      : main_run+'/build/mitgcmuv',
            'binary_dir'    : main_run+'/bin',
            'namelist_dir'  : main_run+'/input',
            'time'          : '72:00:00'}

    slurm = {'be_nice':True,'max_job_submissions':9,'dependency':'afterany'}

    # --- Launch
    driver = pm.SampleDriver(f'sample-3D')
    driver.start(dirs=dirs,
                 dsim=dsim,
                 mymodel=ds['maskC'],
                 ctrl_ds=ds,
                 NxList=[5, 10, 15, 20],
                 xiList=[.5, 1, 2],
                 sorDict={.5:1.6, 1:1.3, 2:1.06},
                 slurm=slurm,
                 n_samples=1000,
                 smooth2DDims=None)
