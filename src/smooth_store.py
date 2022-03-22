
import os
import numpy as np
import xarray as xr
import xmitgcm


def open_smoothdataset(run_dir,
                       xdalike,
                       read_filternorm=False,
                       n_samples=1000,
                       **kwargs):
    """
    TODO:
        - don't pass xdalike ... because would need one for every operator type
        - ideally... write smooth package to use diagnostics to write out the smoothed field
        - read_filternorm necessary for OIDriver still?
        - split up different operators, remove number from variable name (add as dimension)
        - add 2D variables... etc
    """

    # get all smooth*.data prefixes, then their operator numbers
    prefixes = [x for x in os.listdir(run_dir) if 'smooth' in x and '.data' in x]
    numbers = [x.split('.')[0][-3:] for x in prefixes]

    # Right now operators will be a list of strings, since going to want string not int anyway...
    operators = []
    for x in numbers:
        if x not in operators:
            operators.append(x)

    extras = [get_smooth_vars(op) for op in operators]
    custom_grid_variables = {}
    for fld in extras:
        custom_grid_variables.update(fld)

    xds = xmitgcm.open_mdsdataset(run_dir, custom_grid_variables=custom_grid_variables, **kwargs)
    smoothvars = [vname for vname in xds.coords if 'smooth' in vname]
    xds = xds.reset_coords(smoothvars)

    # rechunk to have full 2D field
    # Now get the field
    # HACK: assumes LLC, assumes 3D operator 1, number of samples ... ideally don't need to input
    sample = np.arange(n_samples)
    xds['sample'] = xr.DataArray(sample,
                          coords={'sample':sample},
                          dims=('sample',))
    smoothfld = xmitgcm.utils.read_3d_llc_data(run_dir+'/smooth3Dfld001.data', nx=len(xds.i), nz=len(xds.k),
                                               nrecs=n_samples,
                                               dtype=xds.Depth.dtype,
                                               memmap=False)
    template = xdalike.broadcast_like(xds.sample)
    xds['smooth3Dfld001'] = xr.DataArray(smoothfld, coords=template.coords, dims=template.dims)
    for vname in xds.data_vars:
        xds[vname] = xds[vname].chunk({'face':None})

    # Need to deal with these too...
    xds['ginv_norm'] = xds['smooth3Dfld001'] * xds['smooth3Dnorm001']
    xds['ginv_norm'].attrs = {'label': r'$XC\mathbf{z}$'}
    return xds


def get_smooth_vars(smoothOpNb):

    extra_variables = {
        f'smooth3DKux{smoothOpNb}': {
            'dims': ['k', 'j', 'i_g'],
            'attrs': {
                'standard_name': 'kux',
                'long_name': 'k11'
            }
        },
        f'smooth3DKvy{smoothOpNb}': {
            'dims': ['k', 'j', 'i_g'],
            'attrs': {
                'standard_name': 'kvy',
                'long_name': 'k22'
            }
        },
        f'smooth3DKwz{smoothOpNb}': {
            'dims': ['k_l', 'j', 'i'],
            'attrs': {
                'standard_name': 'kwz',
                'long_name': 'k33'
            }
        },
        f'smooth3DDelta{smoothOpNb}': {
            'dims': ['k', 'j', 'i'],
            'attrs': {
                'standard_name': 'delta',
                'long_name': 'delta'
            }
        },
        f'smooth3Dnorm{smoothOpNb}': {
            'dims': ['k', 'j', 'i'],
            'attrs': {
                'standard_name': 'norm',
                'long_name': 'norm',
            }
        },
        f'smooth3Dmean{smoothOpNb}': {
            'dims': ['k', 'j', 'i'],
            'attrs': {
                'standard_name': 'sample_mean',
                'long_name': 'sample_mean',
            }
        },
        f'smooth3DRandNorm{smoothOpNb}': {
            'dims': ['k', 'j', 'i'],
            'attrs': {
                'standard_name': 'rand_norm',
                'long_name': 'rhs_factor',
            }
        },
    }
    return extra_variables
