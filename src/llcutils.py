import numpy as np
import xarray as xr

def get_pacific(xds):
    xds = xds.isel(face=slice(7,None))
    mlist = []

    for face in [7,8,9]:
        maskc = xds.sel(face=face)
        i = 9-face
        inew = np.arange(89,-1,-1)+i*90
        maskc['i_new'] = xr.DataArray(inew, xds.i.coords, xds.i.dims)
        maskc = maskc.swap_dims({'i':'i_new'}).drop('i').rename({'i_new':'iy'})
        mlist.append(maskc)

    pac1 = xr.concat(mlist, dim='iy').rename({'j':'ix'})

    mlist = []
    for face in [10,11,12]:
        maskc = xds.sel(face=face)
        i = 12-face
        inew = np.arange(89,-1,-1)+i*90
        maskc['i_new'] = xr.DataArray(inew, xds.i.coords, xds.i.dims)
        maskc = maskc.swap_dims({'i':'i_new'}).drop('i').rename({'i_new':'iy'})
        jnew = np.arange(90) + 90
        maskc['j_new'] = xr.DataArray(jnew, xds.j.coords, xds.j.dims)
        maskc = maskc.swap_dims({'j':'j_new'}).drop('j').rename({'j_new':'ix'})
        mlist.append(maskc)

    pac2 = xr.concat(mlist, dim='iy')
    pac = xr.concat([pac1,pac2], dim='ix')
    #pac = pac.sortby(['ix','iy'])
    pac['ix'].attrs = {'long_name':'Real x direction',
                       'standard_name': 'x_index'
                      }
    pac['iy'].attrs = {'long_name':'Real y direction',
                       'standard_name': 'y_index'
                      }
    return pac