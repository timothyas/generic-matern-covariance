
import os
import numpy as np
import xarray as xr
import xmitgcm
from xmitgcm import llcreader
from fsspec.implementations.local import LocalFileSystem

def open_smoothdataset(data_dir, grid_dir=None,
                       geometry='sphericalpolar',
                       read_filternorm=False,
                       **kwargs):
    """Create an xarray Dataset object with pkg/smooth outputs, using
    xmitgcm.open_mdsdataset if not 'llc' geometry, and xmitgcm.llcreader.get_dataset if 'llc'

    TODO:
        - read_filternorm necessary for OIDriver still?

    Parameters
    ----------
    data_dir : str
        Path to directory with smooth*.meta/data files
    grid_dir : str, optional
        Path to smooth*.meta/data files that don't have iteration numbers
    geometry : {'sphericalpolar', 'cartesian', 'llc', 'curvilinear', 'cs'}
        MITgcm grid geometry specifier
    read_filternorm : bool, optional
        TBD
    kwargs : optional
        additional arguments are passed to xmitgcm.open_mdsdataset or
        xmitgcm.llcreader.get_dataset, as relevant

    Returns
    -------
    xds : xarray.Dataset
        with all smooth fields
    """

    # Get all relevant smooth*D files, from data_dir and grid_dir
    extra_variables, custom_grid_variables = get_vars_from_run_dir(data_dir)
    if grid_dir is not None:
        ev2, cgv2 = get_vars_from_run_dir(grid_dir)
        extra_variables.update(ev2)
        custom_grid_variables.update(cgv2)

    if geometry == 'llc':

        fs = LocalFileSystem()
        store = llcreader.BaseStore(fs, base_path=data_dir, grid_path=grid_dir)
        model = SmoothLLC90Model(store)

        # Add grid variables to this store
        model.grid_varnames = list(model.grid_varnames) + list(custom_grid_variables.keys())

        for key in ['iter_start', 'iter_stop', 'iter_step']:
            if key in kwargs:
                setattr(model, key, kwargs.pop(key))

        read_grid = False if grid_dir is None else True
        xds = model.get_dataset(varnames=list(extra_variables.keys()),
                                read_grid=read_grid,
                                extra_variables={**extra_variables,**custom_grid_variables},
                                **kwargs)

    else:
        xds = xmitgcm.open_mdsdataset(data_dir,
                                      grid_dir,
                                      extra_variables=extra_variables,
                                      custom_grid_variables=custom_grid_variables,
                                      geometry=geometry,
                                      **kwargs)

    # Add mask if not there
    for key in ['maskC', 'maskW', 'maskS']:
        if key not in xds:
            xds[key] = xds[key.replace('mask','hFac')] > 0
            xds = xds.set_coords(key)

    # Move all smooth variables to data_vars status
    smoothvars = [vname for vname in xds.coords if 'smooth' in vname]
    xds = xds.reset_coords(smoothvars)

    # These get read in with a "time" dimension, make this "sample"
    xds['sample'] = xr.DataArray(np.arange(len(xds['time'])),
                                 coords=xds['time'].coords,
                                 dims=xds['time'].dims)
    xds = xds.set_coords('sample').swap_dims({'time':'sample'})
    for key in ['time','iter','niter']:
        if key in xds:
            xds = xds.drop(key)

    # Need to deal with these too...
    #xds['ginv_norm'] = xds['smooth3Dfld001'] * xds['smooth3Dnorm001']
    #xds['ginv_norm'].attrs = {'label': r'$XC\mathbf{z}$'}
    #TODO:
    # - ginvnorm thing
    return xds


class SmoothLLC90Model(llcreader.LLC90Model):
    grid_varnames = list(('AngleCS', 'AngleSN', 'DRC', 'DRF',
                          'DXC', 'DXG', 'DXF', 'DXV',
                          'DYC', 'DYG', 'DYF', 'DYU',
                          'Depth', 'PHrefC', 'PHrefF',
                          'RAC', 'RAS', 'RAW', 'RAZ', 'RC', 'RF',
                          'RhoRef', 'XC', 'XG', 'YC', 'YG',
                          'hFacC', 'hFacS', 'hFacW',
                     ))
    iter_start = 1
    iter_stop = 1001
    iter_step = 1


def get_vars_from_run_dir(run_dir):
    """Get all potential smooth*D output from the run directory
    """
    # get all smooth*.data prefixes, then their operator numbers
    prefixes = [x.split('.')[0] for x in os.listdir(run_dir) if 'smooth' in x and '.data' in x]
    numbers = [x[-3:] for x in prefixes]

    # Right now operators will be a list of strings,
    # since going to want string not int anyway...
    # Make this a dictionary, key is dimensionality, values are operator numbers
    operators = {2:[], 3:[]}
    for prefix in prefixes:
        number = prefix.split('.')[0][-3:]
        for d in operators.keys():
            if f'smooth{d}D' in prefix and number not in operators[d]:
                operators[d].append(number)

    # Get possible variable names, and remove ones not in the data directory
    extra_variables = {}
    custom_grid_variables = {}
    for dim, opnumbers in operators.items():
        for op in opnumbers:
            ev  = {key: val for key, val in _extra_variables(dim, op).items() if key in prefixes}
            cgv = {key: val for key, val in _custom_grid_variables(dim, op).items() if key in prefixes}

            extra_variables.update(ev)
            custom_grid_variables.update(cgv)

    return extra_variables, custom_grid_variables


def _custom_grid_variables(n_dim, smoothOpNb):

    custom_grid_variables = {
        f'smooth{n_dim}DKux{smoothOpNb}': {
            'dims': _dims(n_dim, 'W'),
            'attrs': {
                'standard_name': 'kux',
                'long_name': 'k11'
            }
        },
        f'smooth{n_dim}DKvy{smoothOpNb}': {
            'dims': _dims(n_dim,'S'),
            'attrs': {
                'standard_name': 'kvy',
                'long_name': 'k22'
            }
        },
        f'smooth{n_dim}DKwz{smoothOpNb}': {
            'dims': _dims(n_dim,'L'),
            'attrs': {
                'standard_name': 'kwz',
                'long_name': 'k33'
            }
        },
        f'smooth{n_dim}DDelta{smoothOpNb}': {
            'dims': _dims(n_dim, 'C'),
            'attrs': {
                'standard_name': 'delta',
                'long_name': 'delta'
            }
        },
        f'smooth{n_dim}Dnorm{smoothOpNb}': {
            'dims': _dims(n_dim, 'C'),
            'attrs': {
                'standard_name': 'norm',
                'long_name': 'norm',
            }
        },
        f'smooth{n_dim}Dmean{smoothOpNb}': {
            'dims': _dims(n_dim, 'C'),
            'attrs': {
                'standard_name': 'sample_mean',
                'long_name': 'sample_mean',
            }
        },
        f'smooth{n_dim}DRandNorm{smoothOpNb}': {
            'dims': _dims(n_dim, 'C'),
            'attrs': {
                'standard_name': 'rand_norm',
                'long_name': 'rhs_factor',
            }
        },
    }
    return custom_grid_variables

def _extra_variables(n_dim, smoothOpNb):

    extra_variables = {
        f'smooth{n_dim}Dfld{smoothOpNb}': {
            'dims': _dims(n_dim, 'C'),
            'attrs': {
                'standard_name': 'smooth_fld',
                'long_name': r'$C\mathbf{z}$',
            }
        },
    }
    return extra_variables

def _dims(n_dim, loc):
    if loc == 'C':
        return ['k','j','i'] if n_dim ==3 else ['j','i']
    elif loc == 'W':
        return ['k','j','i_g'] if n_dim ==3 else ['j','i_g']
    elif loc == 'S':
        return ['k','j_g','i'] if n_dim ==3 else ['j_g','i']
    elif loc == 'L':
        return ['k_l','j','i'] if n_dim ==3 else ['j','i']
