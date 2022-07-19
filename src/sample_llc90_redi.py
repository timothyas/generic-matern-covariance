"""Run 3D LLC sampling"""

from xmitgcm import open_mdsdataset

from sampledriver import SampleDriver


if __name__ == "__main__":

    main_run='/scratch2/tsmith/generic-matern-covariance/sampling/llc90'
    ds = open_mdsdataset(f"{main_run}/grid", iters=None, geometry='llc')

    # --- directories and dicts
    dirs = {'main_run'      : main_run,
            'netcdf'        : main_run+'/ncfiles'}

    dsim = {'machine'       : 'sverdrup',
            'n_procs'       : 96,
            'exe_file'      : main_run+'/build.redi/mitgcmuv',
            'binary_dir'    : main_run+'/bin',
            'namelist_dir'  : main_run+'/input',
            'time'          : '72:00:00'}

    slurm = {'be_nice':True,'max_job_submissions':4,'dependency':'afterany'}

    # --- Launch
    log10tol = -3
    for n_apps in [2]:
        driver = SampleDriver(f'matern-redi-{n_apps:02d}apps/log10tol{log10tol:03d}-3D')
        driver.start(dirs=dirs,
                     dsim=dsim,
                     mymodel=ds['maskC'],
                     ctrl_ds=ds,
                     NxList=[5, 10, 20],
                     xiList=[1],
                     sorDict={.5:1., 1:1., 2:1.},
                     slurm=slurm,
                     n_samples=1000,
                     smooth2DDims=None,
                     elliptic_tol=10**log10tol,
                     n_applications=n_apps,
                     redi_rotation=True,
                     isoxy=True,
                     density_path='/scratch2/tsmith/generic-matern-covariance/sampling/llc90/zstores/rhoanoma.2014-05.zarr'
                     )

