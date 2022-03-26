"""Compute correlations, for now in k

TODO:
    - Replace this with the notebook version that does it in dask
    - add i/j dimensions when ready
"""


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from dask.distributed import Client, performance_report

from xmitgcm import open_mdsdataset
from smooth_store import open_smoothdataset
from matern import MaternField

def calc_correlation_field(xda, mask,
                           dimlist=['Z','YC'],
                           n_shift=15,
                           sample_mean=None):
    """calculate the correlation field for each shifted distance

    Parameters
    ----------
    xda : xarray.DataArray
        The field to compute correlations on, over the 'sample' dimension
    mask : xarra.DataArray
        True/False inside/outside of domain
    dimlist : list of str
        denoting dimensions to compute shifted correlations
    n_shift : int
        number of shifts to do
    """

    xds = xr.Dataset()
    shifty = np.arange(-n_shift,n_shift+1)
    shifty = xr.DataArray(shifty,coords={'shifty':shifty},dims=('shifty',))
    xds['shifty'] = shifty

    # Do this part once
    sample_mean = xda.mean('sample') if sample_mean is None else sample_mean
    template = xda.isel(sample=0).drop('sample')
    x_deviation = (xda - sample_mean).where(mask)
    x_ssr = np.sqrt( (x_deviation**2).sum('sample') )
    for dim in dimlist:
        corrfld = f'corr_{dim.lower()}'
        xds[corrfld] = xr.zeros_like(shifty*template)

        for s in shifty.values:
            y_deviation = x_deviation.shift({dim:int(s)})
            numerator = (x_deviation*y_deviation).sum('sample')

            y_ssr = np.sqrt( (y_deviation**2).sum('sample') )
            denominator = x_ssr*y_ssr

            xds[corrfld].loc[{'shifty':s}] = numerator / denominator
    return xds


def plot_correlation(cdsd, dimlist=['YC','Z']):

    plt.style.use("thesis")
    ncols = len(dimlist)
    nrows = 1

    fig,axs = plt.subplots(nrows, ncols,
                           figsize=(18,6*nrows),
                           sharey=True)

    axs = [axs] if ncols*nrows == 1 else axs

    for col, (dim, ax) in enumerate(zip(['i'],axs)):
        for curve, (n_range, xds) in enumerate(cdsd.items()):

            # Plot ideal correlation curve
            xaxis  = xds[f"delta_{dim}hat"]
            plotme = xds[f"ideal_{dim.lower()}"]
            label = r"$r\,(\hat{\rho},||\hat{x}_1-\hat{x}_2||)$" if col*curve==0 else None

            ax.plot(xaxis, plotme, label=label, color='black')


            # Plot empirical correlation curve
            fld = f'corr_{dim.lower()}'
            label=r'$\hat{\rho}$ = %d' % n_range if col*curve == 0  else None

            ax.fill_between(xaxis,
                            xds[f"avg_{fld}"]-xds[f"std_{fld}"],
                            xds[f"avg_{fld}"]+xds[f"std_{fld}"],
                            alpha=.3,
                            label=label)

        ax.set(xlabel=r'$\delta\hat{%s}$'%dim[0].lower(),ylabel='',title='')
        if col == 0:
            ax.set_ylabel('Correlation')

    fig.subplots_adjust(wspace=.05)
    fig.legend(ncol=len(dimlist)+1,
               loc='center',
               bbox_to_anchor=(.5,-0.075),
               frameon=False)
    return fig, axs


if __name__ == "__main__":

    client = Client()

    n_range = 10
    horizontal_factor = 1
    main_dir = '/scratch2/tsmith/generic-matern-covariance/sampling/llc90/sample-3D-C'
    run_dir = f'{main_dir}/run.{n_range:02d}dx.{horizontal_factor:02}xi'
    kw = {'iters': None, 'geometry': 'llc'}
    chunks = {'k': None, 'face': None}
    filename = "performance/correlation_k.html"

    with performance_report(filename=filename):

        # Open and rechunk
        ds = open_mdsdataset(run_dir, **kw)
        ds = open_smoothdataset(run_dir, xdalike=ds.maskC, **kw)

        for key in ['smooth3Dmean001', 'smooth3Dfld001', 'ginv_norm', 'smooth3Dnorm001', 'maskC']:
            ds[key] = ds[key].chunk(chunks)

        # Compute
        cds = calc_correlation_field(ds['ginv_norm'],
                                     mask=ds['maskC'],
                                     dimlist=['k'],
                                     sample_mean=ds['smooth3Dmean001'],
                                     n_shift=horizontal_factor*2*n_range)

        # These operations should get wrapped into a function
        mf = MaternField(ds['maskC'], n_range=n_range, horizontal_factor=horizontal_factor)
        cds['dist'] = np.abs(cds['shifty'])
        cds['delta_khat'] = cds['dist'].copy(deep=True)
        cds['ideal_k'] = mf.ideal_correlation(distance=cds['delta_khat'])

        # compute emirical correlations, be persistent
        cds['avg_corr_k'] = cds['corr_k'].mean(['k','face','j','i']).compute()
        cds['std_corr_k'] = cds['corr_k'].std(['k','face','j','i']).compute()

        fig, _ = plot_correlations({n_range: cds}, dimlist=['k'])
        fig.savefig("../figures/nondimensional_correlation.pdf", bbox_inches="tight")

        cds.to_zarr("/scratch2/tsmith/generic-matern-covariance/sample/llc90/ncfiles/correlation_k.zarr")
