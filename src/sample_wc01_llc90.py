"""Run 3D LLC sampling"""

from xmitgcm import open_mdsdataset

from wcsampledriver import WCSampleDriver


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
    for n_time_steps in [100, 150, 200, 250, 300, 400, 500, 1000, 1500, 2000]:
        driver = WCSampleDriver(f'wc01-3D-nbt{n_time_steps:03d}')
        driver.start(dirs=dirs,
                     dsim=dsim,
                     mymodel=ds['maskC'],
                     ctrl_ds=ds,
                     NxList=[10],
                     xiList=[1],
                     n_time_steps=n_time_steps,
                     slurm=slurm,
                     n_samples=1000)
