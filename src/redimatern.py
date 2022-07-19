
import numpy as np
import xarray as xr

import ecco_v4_py

from matern import MaternField

class RediMaternField(MaternField):

    def __init__(self, xdalike, n_range, horizontal_factor,
                 density_path,
                 isotropic=False,
                 isoxy=True,
                 n_applications=1):

        try:
            assert isoxy
        except:
            raise NotImplementedError("Only isoxy=True is implemented ...")

        try:
            assert not isotropic
        except:
            raise NotImplementedError("isotropic version not implemented")

        self.density_path = density_path
        super().__init__(xdalike=xdalike,
                         n_range=n_range,
                         horizontal_factor=horizontal_factor,
                         isotropic=isotropic,
                         isoxy=isoxy,
                         n_applications=n_applications)

        try:
            assert self.n_dims == 3
        except:
            raise NotImplementedError("Redi version only implemented for 3D fields")




    def get_operator_elements(self):
        """Note that det(Phi) remains unchanged since rotation is unitary,
        so just need to modify laplacian tensor

        Simplified tensor looks like

                    1  0  Sx
        K = Lx^2 *  0  1  Sy
                    Sx Sy Kwz

        Kwz = Sx^2 + Sy^2 + (Lz^2 / Lx^2)
        Sx and Sy are slopes
            $ Sx = - \partial_x\sigma / \partial_z\sigma $
            $ Sy = - \partial_y\sigma / \partial_z\sigma $
        """

        super().get_operator_elements()

        open_dataset = xr.open_dataset if 'zarr' not in self.density_path else xr.open_zarr
        ds = open_dataset(self.density_path).squeeze()
        grid = ecco_v4_py.get_llc_grid(ds)

        ds['drCl'] = xr.DataArray(ds.drC[:-1].values, coords=ds.k_l.coords, dims=('k_l',))

        rhox = grid.diff(ds.RHOAnoma.where(ds.maskC), 'X', boundary='fill', fill_value=np.nan) / ds.dxC
        rhoy = grid.diff(ds.RHOAnoma.where(ds.maskC), 'Y', boundary='fill', fill_value=np.nan) / ds.dyC
        # need a negative due to flipped z axis with MITgcm
        rhoz = -grid.diff(ds.RHOAnoma.where(ds.maskC), 'Z', boundary='fill', fill_value=np.nan) /ds.drCl

        Sx = xr.DataArray(-rhox.values / rhoz.values, coords=ds.maskC.coords, dims=ds.maskC.dims)
        Sy = xr.DataArray(-rhoy.values / rhoz.values, coords=ds.maskC.coords, dims=ds.maskC.dims)

        # Mask it
        Sx = Sx.where(~np.isnan(Sx), 0.)
        Sy = Sy.where(~np.isnan(Sy), 0.)

        # Rename tile->face
        Sx = Sx.rename({'tile':'face'})
        Sy = Sy.rename({'tile':'face'})
        ds = ds.rename({'tile':'face'})

        # Kux and Kvy remain the same, modify others
        LzInv = 1/self.Lz
        LxInv = xr.where(np.isnan(self.Lx) | (self.Lx==0.), 0., 1/self.Lx)
        self.K['wz'] = LzInv * (Sx**2 + Sy**2 + (self.Lz*LxInv)**2)
        self.K['uz'] = LzInv * Sx
        self.K['vz'] = LzInv * Sy
        self.K['wx'] = LzInv * Sx
        self.K['wy'] = LzInv * Sy
