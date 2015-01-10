from omi  import OMI
from gfio import GFIO

if __name__ == "__main__":

    """Simple unit testing..."""

    from mie import getMie
   
#   Read OMI data
#   -------------
    nymd = 20080629
    nhms = 0
    omi = OMI('ods/orig/omi.aero_tc8.obs.20080629.ods',nymd,nhms)

    aerosols = 'data/a5_arctas_02.inst3d_aer_v.20080629_0000z.nc'

    g = GFIO(aerosols)

    ps = g.interp('ps',omi.lon,omi.lat)


