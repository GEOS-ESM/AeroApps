from pyobs import LIDAR_L2, NPZ
from numpy import savez, load, tile

from mie    import getTopo, getAlbedo
from mieobs import *    

MISSING = -9999
def vlidort_scalar (c, channels, ps, albedo, tau, ssa, g, pe, ze, te, verbose=0 ):
        """

        *** THIS NEED TO BE FINISHED AND TESTED ****
        
        Uses VLIDORT to compute radiances. On Input,

           tau   --- aerosol optical depth
           ssa   --- aerosol single scattering albedo
           g     --- aerosol asymmetry factor
           pe    --- pressure at edges [Pa]
           ze    --- height at edges [m]
           te    --- temperature at edges [K]

        """

        nch = len(channels)
        N = len(c.lon)
        
        radiance, ai, reflectivity, rc = OMI_.scalar(nymd,nhms,channels,
                                                     tau, ssa, g, pe, he, te, 
                                                     c.lon, c.lat, ps, albedo,
                                                     solar_zenith, relat_azymuth,
                                                     sensor_zenith, MISSING, radiance, 
                                                     ai, verbose)
        if rc != 0:
            raise ValueError, "on return from OMI_.scalar/vector, rc = "+str(rc)

        return (radiance, ai, reflectivity)


#-------------------------------------------------------------------------------
if __name__ == "__main__":

    channels = [532.,]

    topo_fn  = '/nobackup/1/VLIDORT/topo/topography.1152x721.nc'
    albedo_fn = '/nobackup/MODIS/Level3/ALBEDO/albedo_clim.ctl'
    
    c_dn = '/nobackup/LIDAR/dR_Fortuna-2-4-b4/' # Pete's Calipso interp
    a_dn = '/home/adasilva/GAAS/LIDAR/'         # aer_Nv + interp

    c_fn = 'dR_Fortuna-2-4-b4.calipso_532nm.20090715.nc' # Pete's interp
    a_fn = 'aer_obs/dR_Fortuna-2-4-b4.calipso_aer.20090715.npz' # aer collocation

    # Read relevant data
    # ------------------
    c = LIDAR_L2(c_dn+c_fn)   # Pete's interp with CALPSO coordinates
    a = NPZ(a_dn+a_fn)        # aer_v interpolated to obs location

    km = a.delp.shape[1]
    
#   ---------------------------------------------------------------------------    
#   NOTE:
#
#   For creating the npz file above with the aer_v data interpolated to obs ,
#   location you do this:
#      c.sampleFile(aer_v,npzFile=npzFile)   
#   where aer_v is the full, gridded, netcdf file.
#   In principle, you do not need to write the npzFile as the result of the
#   sampling is also available as attribute *sample*. So, instead of reading
#   the npzFile you could do
#     a = c.sample
#   But sampling takes time, so saving a "sampling  file" is convenient.
#   BTW, npzFile is too python specific. Soon, I'll implement a NetCDF option,
#   so you will be able to say
#     c.sample(aer_v,ncFile=NetCDF_filename)
#   See bottom of GMAO_pyobs/pyobs.lidar_l2.py file for an example.
#   ---------------------------------------------------------------------------    

    # Sample albedo at lat lon
    # ------------------------
    albedo = getAlbedo(albedo_fn,c.tyme[0],c.lon,c.lat)
    savez('albedo.npz',albedo=albedo.data)

    # once you save the albedo, comment the 2 lines above and
    # read the file instead.
    #albedo = load('albedo.npz')['albedo']

    # I need to look into this...
    # ---------------------------
    albedo[albedo<0] = 0.07 # make something up for now
    
    # Get the VLIDORT inputs
    # ----------------------
    zs = getTopo(topo_fn,c.lon,c.lat)
    pe, ze, te = getEdgeVars(a)
    tau, ssa, g = getAOPscalar(a)
    he = ze + tile(zs,(km+1,1))

    # call vlidort
