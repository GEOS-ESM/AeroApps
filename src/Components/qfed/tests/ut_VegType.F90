   program TestVeg

     use VegType_Mod

     type(VegType) :: veg
     integer, parameter :: nobs = 21
     real :: lons(nobs), lats(nobs)
     integer :: iVeg(nobs)

     lons = (/ &
157.1995, 142.1126, 142.1061, 146.063,  143.3694, 143.364,  143.393, 142.0272, &
148.1296, 147.0529, 147.1054, 142.6649, 147.2151, 147.2297, 147.213, 147.2276, &
147.1822, 147.211,  147.2089, 147.2235, 142.6765 /)

     lats = (/ &
 -8.013471, 18.0738,  18.07808, 19.25864, 19.04274, 19.04079, 19.04668, &
-19.00074,  20.15362, 20.024,   20.0794,  19.31606, 20.19115, 20.19373, &
 20.20216, -20.20471, 20.20818, 20.21319, 20.22419, 20.22671, 19.41586 /)

     call VegType_Initialize(veg,'/output/Emissions/Vegetation/GL_OGE_INPE/OGE')
     call VegType_GetSimple(veg, lons, lats, nobs, iVeg)
     do  i = 1, nobs
        print *, lons(i), lats(i), iveg(i)
     end do
     call VegType_Finalize(veg)

   end program TestVeg
