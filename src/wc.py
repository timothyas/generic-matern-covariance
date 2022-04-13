import os
import numpy as np
import xarray as xr
import xmitgcm

from MITgcmutils import wrmds
import ecco_v4_py
import scipy.special

# Some issues
# - with LLC, ordering by k is not the same as ordering by Z ... so this would not work in the OIDriver ...
# - Here... not worrying about any of this

from matern import _get_dims,  _horizontal_length_scale, _vertical_length_scale

class WCField():

    def __init__(self, xdalike, n_range, horizontal_factor,
                 isotropic=False):

        self.xdalike = xdalike
        self.n_range = n_range
        self.horizontal_factor = horizontal_factor
        self.isotropic = isotropic

        # Set some properties that don't change
        self.dims = xdalike.dims
        self.llc = 'face' in xdalike.dims or 'tile' in xdalike.dims
        self.n_dims = xdalike.ndim - 1 if self.llc else xdalike.ndim
        self.xyz = _get_dims(self.dims)

        try:
            assert self.n_dims == 2 or self.n_dims == 3
        except:
            raise TypeError(f"Only 2D or 3D data, got n_dims={self.n_dims} with data {self.xdalike}")

        if self.llc and self.n_dims == 3 and self.xyz['z'] != 'k':
            raise NotImplementedError(f"vertical dimension needs to be 'k' due to line 1709 in xmitgcm.utils...")

        if self.llc and self.n_dims == 2 and self.xyz['z'] is not None:
            raise NotImplementedError(f"Due to writing routines (and probably other problems), can't do X-Z or Y-Z slice in LLC grid ... it's a global grid anyway!")

        self.horizontal_length = _horizontal_length_scale(xdalike)
        self.vertical_length = _vertical_length_scale(xdalike)
        self.cell_volume = self.get_cell_volume()
        self.L = self.get_appropriate_lengths()


    def get_appropriate_lengths(self):

        if self.n_dims == 2:

            if self.xyz['z'] is None:
                L = {'x': self.horizontal_length,
                     'y': self.horizontal_length}

            elif self.xyz['y'] is None:
                L = {'x': self.horizontal_factor * self.horizontal_length,
                     'z': self.vertical_length}

            elif self.xyz['x'] is None:
                L = {'y': self.horizontal_factor * self.horizontal_length,
                     'z': self.vertical_length}

            else:
                raise TypeError("get_appropriate_lengths dims problem 2d")

        else:

            L = {'x': self.horizontal_factor * self.horizontal_length,
                 'y': self.horizontal_factor * self.horizontal_length,
                 'z': self.vertical_length}

        for key in L.keys():
            L[key] = L[key].broadcast_like(self.xdalike)

        return L


    def get_cell_volume(self):

        if self.n_dims == 2:
            if self.xyz['z'] is None:
                V = self.horizontal_length**2

            else:
                V = self.vertical_length * self.horizontal_length

        else:
            V = self.vertical_length * self.horizontal_length**2

        return V


    def write_binaries(self, write_dir, smoothOpNb, dtype=None):

        if not os.path.isdir(write_dir):
            os.makedirs(write_dir)

        dataprec = str(self.horizontal_length.dtype) if dtype is None else dtype
        dataprec = "float32" if "f4" in dataprec else dataprec
        dataprec = "float64" if "f8" in dataprec else dataprec

        # Write length scales
        for key in self.L.keys():
            writeme = self.prepare_to_write(self.L[key])

            fname = os.path.join(write_dir, f"smooth{self.n_dims}DL{key}{smoothOpNb:03}")
            wrmds(fname, writeme, dataprec=dataprec)


    def prepare_to_write(self, xda):

        # If vertical dimension is present, we want to sort it from surface to depth
        # this is "descending" sort for Z (or Zl, Zu, Zp1...) coordinate, and
        # and "ascending" for k (or k_l, etc..)
        writeme = xda.where(~np.isnan(xda), 0.)
        if self.xyz['z'] is not None:
            ascending = False if 'Z' in self.xyz['z'] else True
            writeme = writeme.sortby(self.xyz['z'], ascending=ascending)

        # Deal with LLC
        if self.llc:

            writeme = ecco_v4_py.llc_tiles_to_compact(writeme.values)

        else:
            writeme = writeme.values

        return writeme
