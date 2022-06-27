import numpy as np
import xarray as xr

from matern import MaternField

class BarotropicRadiusMaternField(MaternField):
    """Matern field with horizontal length scales defined by
    the Barotropic Rossby radius of deformation

    Uses "regularized" form of Rossby radius following definiton by
    Hallberg, 2013.
    """

    def get_horizontal_length_scale(self, xds):
        """Get horizontal length scale from data array or dataset"""

        xds = xds.to_dataset(name='tmp') if isinstance(xds, xr.DataArray) else xds

        assert 'YC' in xds and 'Depth' in xds, "Must have YC and Depth in dataset"

        g       = 9.81 # m/s^2
        omega   = 7.2921e-5
        earthrad= 6_371_000
        radlat  = np.deg2rad(xds["YC"])
        f       = 2 * omega * np.sin(radlat)
        beta    = 2 * omega * np.cos(radlat) / earthrad

        Cg      = g*xds["Depth"]
        denom   = f**2 + 2*beta*Cg

        Lx = np.sqrt(Cg**2 / denom)
        Lx = Lx.astype(xds["YC"].dtype)
        Ly = Lx.copy()

        Lx.name = 'Lx'
        Ly.name = 'Ly'
        return Lx, Ly
