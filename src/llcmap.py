"""
Map class for creating lon/lat projections of LLC grid
This was originally copied form the xgcm documentation:

    https://xgcm.readthedocs.io/en/latest/example_eccov4.html#A-Pretty-Map

Because wow, it's well written
"""

from copy import copy
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import pyresample as pr

class PacificMap:

    def __init__(self, ds, dx=0.25, dy=0.25):
        # Extract LLC 2D coordinates
        #lons = xr.where(ds.XC>0, ds.XC-360, ds.XC)
        lons_1d = ds.XC.values.ravel()
        lats_1d = ds.YC.values.ravel()

        # Define original grid
        self.orig_grid = pr.geometry.SwathDefinition(lons=lons_1d, lats=lats_1d)

        # Longitudes latitudes to which we will we interpolate
        lon_tmp = np.arange(-180, 180, dx) + dx/2
        lat_tmp = np.arange(-90, 90, dy) + dy/2

        # Define the lat lon points of the two parts.
        self.new_grid_lon, self.new_grid_lat = np.meshgrid(lon_tmp, lat_tmp)
        self.new_grid  = pr.geometry.GridDefinition(lons=self.new_grid_lon,
                                                    lats=self.new_grid_lat)

    def __call__(self, da, ax=None, projection=ccrs.Robinson(central_longitude=-120),
                 lon_0=-120,
                 lon_bds=[-180, -60],
                 lat_bds=[-30, 30],
                 show_cbar=True,
                 cbar_label='',
                 **plt_kwargs):

        if ax is None:
            _, ax = plt.subplots(subplot_kw={'projection':projection})

        field = self.regrid(da)

        vmax = plt_kwargs.pop('vmax', np.nanmax(field))
        vmin = plt_kwargs.pop('vmin', np.nanmin(field))
        if vmax*vmin < 0:
            vmax = np.nanmax(np.abs([vmax,vmin]))
            vmin = -vmax

        # Handle colorbar and NaN color
        default_cmap = 'RdBu_r' if vmax*vmin < 0 else 'viridis'
        cmap = plt_kwargs.pop('cmap', default_cmap)
        if type(cmap)==str:
            cmap = copy(plt.cm.get_cmap(cmap))
        #cmap.set_bad(color='gray',alpha=.6)

        if lon_bds is not None and lat_bds is not None:
            ax.set_extent([lon_bds[0],lon_bds[1],lat_bds[0],lat_bds[1]], ccrs.PlateCarree())
        x,y = self.new_grid_lon, self.new_grid_lat

        # Plot each separately
        p = ax.pcolormesh(x, y, field,
                           vmax=vmax, vmin=vmin, cmap=cmap,
                           transform=ccrs.PlateCarree(), zorder=1, **plt_kwargs)

        # Add land and coastlines
        ax.add_feature(cf.LAND.with_scale('50m'), zorder=3, color=[.8, .8, .8])
        ax.add_feature(cf.COASTLINE.with_scale('50m'), zorder=3)

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(),
                          xlocs=np.arange(lon_bds[0], lon_bds[1]+1, 30),
                          ylocs=np.arange(lat_bds[0], lat_bds[1]+1, 15),
                          draw_labels=True,
                          linewidth=1, color='gray', alpha=0.2, linestyle='-')
        gkw = {}
        spec = ax.get_subplotspec()
        gkw["left_labels"] = spec.is_first_col()
        gkw["right_labels"] = spec.is_last_col()
        gkw["top_labels"] = spec.is_first_row()
        gkw["bottom_labels"] = spec.is_last_row()
        for key, val in gkw.items():
            setattr(gl, key, val)

        # Colorbar...
        if show_cbar:
            cb=plt.colorbar(p,ax=ax,shrink=.8,label=cbar_label,
                            orientation='horizontal',pad=0.05)

        return p, gl


    def regrid(self,xda):
        """regrid xda based on llcmap grid"""
        return pr.kd_tree.resample_nearest(self.orig_grid, xda.values,
                                           self.new_grid,
                                           radius_of_influence=100000,
                                           fill_value=None)
