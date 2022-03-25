
import os
import numpy as np
import xarray as xr
import xmitgcm

def open_smoothdataset(run_dir,
                       read_filternorm=False,
                       **kwargs):
    """
    Note:
        This is a really clean version of open_smoothdataset, but it relies on having
        the smooth*Dfld written out with the 10 digit iteration corresponding to the sample
        number. It turns out that this still gives xmitgcm.open_mdsdataset a really hard time.
        That's why I cooked up read_and_store.py... to do this the eager and not
        generalized way - it worked.


    TODO:
        - read_filternorm necessary for OIDriver still?
        - split up different operators, remove number from variable name (add as dimension)
        - add 2D variables... etc
    """

    # get all smooth*.data prefixes, then their operator numbers
    prefixes = [x for x in os.listdir(run_dir) if 'smooth' in x and '.data' in x]
    numbers = [x.split('.')[0][-3:] for x in prefixes]

    # Right now operators will be a list of strings,
    # since going to want string not int anyway...
    operators = []
    for x in numbers:
        if x not in operators:
            operators.append(x)

    # Get the variable names
    extra_variables = {}
    custom_grid_variables = {}
    for op in operators:
        extra_variables.update(_extra_variables(op))
        custom_grid_variables.update(_custom_grid_variables(op))

    xds = xmitgcm.open_mdsdataset(run_dir,
                                  extra_variables=extra_variables,
                                  custom_grid_variables=custom_grid_variables,
                                  **kwargs)

    # Move all smooth variables to data_vars status
    smoothvars = [vname for vname in xds.coords if 'smooth' in vname]
    xds = xds.reset_coords(smoothvars)

    # These get read in with a "time" dimension, make this "sample"
    xds['sample'] = xr.DataArray(np.arange(len(xds['time'])),
                                 coords=xds['time'].coords,
                                 dims=xds['time'].dims)
    xds = xds.set_coords('sample').swap_dims({'time':'sample'})
    xds = xds.drop(['time','iter'])

    # Need to deal with these too...
    xds['ginv_norm'] = xds['smooth3Dfld001'] * xds['smooth3Dnorm001']
    xds['ginv_norm'].attrs = {'label': r'$XC\mathbf{z}$'}
    return xds


def _custom_grid_variables(smoothOpNb):

    custom_grid_variables = {
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
    return custom_grid_variables

def _extra_variables(smoothOpNb):

    extra_variables = {
        f'smooth3Dfld{smoothOpNb}': {
            'dims': ['k', 'j', 'i'],
            'attrs': {
                'standard_name': 'smooth_fld',
                'long_name': r'$C\mathbf{z}$',
            }
        },
    }
    return extra_variables
