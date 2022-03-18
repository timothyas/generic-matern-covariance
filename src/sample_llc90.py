
import numpy as np
import xarray as xr
from xmitgcm import open_mdsdataset
import pych.pigmachine as pm

if __name__ == "__main__":

    ds = open_mdsdataset('/scratch/tsmith/grids/llc90', iters=None)
    ds = ds.sortby(['Z','YC','XC'])

    # --- directories and dicts
    main_run='/scratch2/tsmith/generic-matern-covariance/sampling/llc90'
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
    driver = pm.SampleDriver(f'llc90-3D')
    driver.start(dirs=dirs,
                 dsim=dsim,
                 mymodel=ds['maskC'],
                 ctrl_ds=ds,
                 NxList=[15],
                 xiList=[2],
                 slurm=slurm,
                 n_samples=1000,
                 smooth2DDims=None,
                 )
