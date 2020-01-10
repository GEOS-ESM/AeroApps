Module LidarAngles


contains

! NAME:
!	JULDAY
!
! PURPOSE:
!	Calculate the Julian Day Number for a given month, day, and year.
!	This is the inverse of the library function CALDAT.
!	See also caldat, the inverse of this function.
!
! INPUTS:
!	MONTH:	Number of the desired month (1 = January, ..., 12 = December).
!
!	DAY:	Number of day of the month.
!
!	YEAR:	Number of the desired year. Year parameters must be valid
!              values from the civil calendar. Years B.C.E. are represented
!              as negative integers. Years in the common era are represented
!              as positive integers. In particular, note that there is no
!              year 0 in the civil calendar. 1 B.C.E. (-1) is followed by
!              1 C.E. (1).
!
!	HOUR:	Number of the hour of the day.
!
!	MINUTE:	Number of the minute of the hour.
!
!	SECOND:	Number of the second of the minute (double precision number).
!
! OUTPUTS:
!	JULDAY returns the Julian Day Number (which begins at noon) of the
!	specified calendar date.
!
! RESTRICTIONS:
!	Accuracy using IEEE double precision numbers is approximately
!      1/10000th of a second, with higher accuracy for smaller (earlier)
!      Julian dates.
!
! MODIFICATION HISTORY:
!	Translated from "Numerical Recipies in C", by William H. Press,
!	Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
!	Cambridge University Press, 1988 (second printing).
!
!	AB, September, 1988
!	DMS, April, 1995, Added time of day.
!
!      Translated into Fortran90 by Vijay Natraj, JPL, February 24 2012

function JULDAY(MONTH, DAY, YEAR, Hour, Minute, Second)

implicit none

!  Inputs

integer(kind=4), intent(in)  :: MONTH
integer(kind=4), intent(in)  :: DAY
integer(kind=4), intent(in)  :: YEAR
integer(kind=4), intent(in)  :: HOUR
integer(kind=4), intent(in)  :: MINUTE
real(kind=8),    intent(in)  :: SECOND

!  Outputs

real(kind=8)                 :: JULDAY

!  Local variables

integer(kind=4)              :: GREG
integer(kind=4)              :: min_calendar
integer(kind=4)              :: max_calendar
integer(kind=4)              :: bc
integer(kind=4)              :: L_YEAR
integer(kind=4)              :: inJanFeb
integer(kind=4)              :: JY
integer(kind=4)              :: JM
integer(kind=4)              :: JA
integer(kind=4)              :: JUL
real(kind=8)                 :: eps

! Gregorian Calendar was adopted on Oct. 15, 1582
! skipping from Oct. 4, 1582 to Oct. 15, 1582

GREG = 2299171  ! incorrect Julian day for Oct. 25, 1582

! check if date is within allowed range

min_calendar = -4716
max_calendar = 5000000
IF ((YEAR .LT. min_calendar) .OR. (YEAR .GT. max_calendar)) &
   write(0,*) 'Value of Julian date is out of allowed range'

IF (YEAR .LT. 0) THEN
   bc = 1
ELSE
   bc = 0
ENDIF
L_YEAR = YEAR + bc

IF (MONTH .LE. 2) THEN
   inJanFeb = 1
ELSE
   inJanFeb = 0
ENDIF

JY = YEAR - inJanFeb
JM = MONTH + 1 + 12*inJanFeb

JUL = FLOOR(365.25d0 * JY) + FLOOR(30.6001d0 * JM) + DAY + 1720995

! Test whether to change to Gregorian Calendar.

IF (JUL .GE. GREG) THEN
   JA = FLOOR(0.01d0 * JY)
   JUL = JUL + 2 - JA + FLOOR(0.25d0 * JA)
ENDIF

! Add a small offset so we get the hours, minutes, & seconds back correctly
! if we convert the Julian dates back. This offset is proportional to the
! Julian date, so small dates (a long, long time ago) will be "more" accurate.

! eps = (MACHAR(/DOUBLE)).eps

eps = 2.2204460d-16 ! For Ganesha, calculated from IDL output of above statement

IF (ABS(JUL) .GE. 1) THEN
   eps = eps*ABS(JUL)
ENDIF

! For Hours, divide by 24, then subtract 0.5, in case we have unsigned integers.

JULDAY = real(JUL,kind=8) + real(Hour,kind=8)/24.d0 - &
         0.5d0 + real(Minute,kind=8)/1440.d0 + &
         Second/86400.d0 + eps

RETURN
END FUNCTION JULDAY

! NAME:
!	CALDAT
!
! PURPOSE:
!	Return the calendar date and time given julian date.
!	This is the inverse of the function JULDAY.
!
! INPUTS:
!	JULIAN contains the Julian Day Number (which begins at noon) of the
!	specified calendar date.  It should be a long integer.
!
! OUTPUTS:
!	(Trailing parameters may be omitted if not required.)
!	MONTH:	Number of the desired month (1 = January, ..., 12 = December).
!
!	DAY:	Number of day of the month.
!
!	YEAR:	Number of the desired year.
!
!	HOUR:	Hour of the day
!
!	Minute: Minute of the day
!
!	Second: Second (and fractions) of the day.
!
! RESTRICTIONS:
!	Accuracy using IEEE double precision numbers is approximately
!	1/10000th of a second.
!
! MODIFICATION HISTORY:
!	Translated from "Numerical Recipies in C", by William H. Press,
!	Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
!	Cambridge University Press, 1988 (second printing).
!
!	DMS, July 1992.
!	DMS, April 1996, Added HOUR, MINUTE and SECOND keyword
!	AB, 7 December 1997, Generalized to handle array input.
!	AB, 3 January 2000, Make seconds output as DOUBLE in array output.
!      CT, Nov 2006: For Hour/Min/Sec, tweak the input to make sure hours
!      and minutes are correct. Restrict hours to 0-23 & min to 0-59.
!
!      Translated into Fortran90 by Vijay Natraj, JPL, February 24 2012

SUBROUTINE CALDAT(julian, month, day, year, hour, minute, second)

implicit none

!  Inputs
  
real(kind=8), intent(in)      :: julian

!  Outputs

integer(kind=4), intent(out)  :: month
integer(kind=4), intent(out)  :: day
integer(kind=4), intent(out)  :: year
integer(kind=4), intent(out)  :: hour
integer(kind=4), intent(out)  :: minute
integer(kind=4), intent(out)  :: second
   
!  Local variables

real(kind=8)                  :: min_julian
real(kind=8)                  :: max_julian
integer(kind=4)               :: igreg
integer(kind=4)               :: julLong
integer(kind=4)               :: jalpha
integer(kind=4)               :: ja
integer(kind=4)               :: jb
integer(kind=4)               :: jc
integer(kind=4)               :: jd
integer(kind=4)               :: je
real(kind=8)                  :: fraction
real(kind=8)                  :: eps

min_julian = -1095.d0
max_julian = 1827933925.d0
IF ((julian .LT. min_julian) .OR. (julian .GT. max_julian)) THEN
   write(0,*) 'Value of Julian date is out of allowed range.'
ENDIF

igreg = 2299161                 ! Beginning of Gregorian calendar
julLong = FLOOR(julian + 0.5d0) ! Better be long

IF (julLong .GE. igreg) THEN    ! Gregorian
   jalpha = FLOOR(((julLong - 1867216) - 0.25d0) / 36524.25d0)
   ja = julLong + 1 + jalpha - FLOOR(0.25d0 * jalpha)
ELSE IF (julLong .LT. 0) THEN
   ja = julLong + 36525 * (1-julLong/36525)
ELSE
   ja = julLong
ENDIF

jb = ja + 1524
jc = FLOOR(6680.d0 + ((jb-2439870)-122.1d0)/365.25d0)
jd = FLOOR(365.d0 * jc + (0.25d0 * jc))
je = FLOOR((jb - jd) / 30.6001d0)

day = jb - jd - FLOOR(30.6001d0 * je)
month = je - 1
month = MOD(month-1,12) + 1
year = jc - 4715
IF (month .GT. 2) year = year - 1
IF (year .LE. 0) year = year - 1
IF (julLong .LT. 0) year = year - 100 * (1-julLong/36525)

!  hours, minutes, seconds

fraction = julian + 0.5d0 - julLong
IF (ABS(julLong) .GE. 1) THEN
  eps = 1.d-12*ABS(julLong)
ELSE
  eps = 1.d-12
ENDIF
hour = FLOOR(fraction * 24.d0 + eps)
fraction = fraction - hour/24.d0
minute = FLOOR(fraction*1440.d0 + eps)
second = (fraction - minute/1440.d0)*86400.d0

RETURN
END SUBROUTINE CALDAT

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    solar_angles
!
! PURPOSE
!   Calculate Elevation and Azimuth Angles of Sun at given location and time
!
! INPUT
!    yr  = 4-digit year
!    mon = month of year, 1-12
!    day = day of month
!    hr  = hour, 0-23
!    min = minute, 0-59
!    sec = second, [0.d0,60.d0)
!    tz  = local time zone (0.d0 = UTC, -8.d0 = PST)
!    lat = Latitude of location in degrees, positive North, negative South
!    lon = Longitude of location in degrees, positive East, negative West
!    
! OUTPUT
!    solar_angles is returned array of length 2.
!    solar_angles(1) is Solar Azimuth Angle in Degrees, 0-360.
!    solar_angles(2) is Solar Elevation Angle in Degrees, 0-90. Negative is Sun below Horizon.
!
! REFERENCES
!    NOAA Solar Calculator at http://www.esrl.noaa.gov/gmd/grad/solcalc/
!      (Decl, Eqtime, Az, and EL agree with values returned by above.)
!    General info about calculation details is at 
!      http://www.srrb.noaa.gov/highlights/sunrise/calcdetails.html
!         "The calculations in the NOAA Sunrise/Sunset and Solar Position Calculators
!          are based on equations from Astronomical Algorithms, by Jean Meeus."
!    Code used by NOAA Solar Calculator to get Declination and Equation of Time is at
!       http://www.srrb.noaa.gov/highlights/sunrise/program.txt
!
! HISTORY
!    written February 2010 by Joyce Wolf for Annmarie Eldering and Susan Kulawik (JPL)
!    translated into Fortran90 24 Feb 2012 by Vijay Natraj (JPL)
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    

FUNCTION solar_angles(yr,mo,day,hr,min,sec,tz,lat,lon)

implicit none

!  Inputs
  
integer(kind=4), intent(in)  :: yr
integer(kind=4), intent(in)  :: mo
integer(kind=4), intent(in)  :: day
integer(kind=4), intent(in)  :: hr
integer(kind=4), intent(in)  :: min
real(kind=8),    intent(in)  :: sec
real(kind=8),    intent(in)  :: tz
real(kind=8),    intent(in)  :: lat
real(kind=8),    intent(in)  :: lon

!  Outputs
   
real(kind=8)                :: solar_angles(2)
   
!  Local variables

real(kind=8)                 :: d2r, r2d
real(kind=8)                 :: julianday
real(kind=8)                 :: t, meanlon, mean_anom
real(kind=8)                 :: mrad, sinm, sin2m, sin3m, c
real(kind=8)                 :: truelon, true_anom
real(kind=8)                 :: omega, app_lon, e, seconds, eps0
real(kind=8)                 :: epsilon, sint, decl
real(kind=8)                 :: m, lz, y, sin2lz, cos2lz, sin4lz
real(kind=8)                 :: eqtime, utc, ha
real(kind=8)                 :: declr, latr, har
real(kind=8)                 :: elr, yaz, xaz, azr, az, el

d2r = acos(-1.d0)/180.d0
r2d = 180.d0/acos(-1.d0)

! t = Julian centuries past J2000.0
! Adjust if time is local

julianday = JULDAY(mo,day,yr,hr,min,sec)
! print *, julianday
t = (julianday-2451545.d0)/36525.d0 - tz/24.d0

! Geometric Mean Longitude of Sun in degrees

meanlon = 280.46646d0 + t * (36000.76983d0 + 0.0003032d0 * t)
do while (meanlon .GT. 360.d0)
   meanlon = meanlon-360.d0
end do
do while (meanlon .LT. 0.d0)
   meanlon = meanlon + 360.d0
end do

! Geometric Mean Anomaly of Sun in degrees

mean_anom = 357.52911d0 + t * (35999.05029d0 - 0.0001537d0 * t)

! Sun Equation of Center

mrad = d2r*mean_anom
sinm = sin(mrad)
sin2m = sin(2.d0*mrad)
sin3m = sin(3.d0*mrad)
c = sinm * (1.914602d0 - t * (0.004817d0 + 0.000014d0 * t)) + &
    sin2m * (0.019993d0 - 0.000101d0 * t) + sin3m * 0.000289d0

! Sun's True Longitude in degrees

truelon = meanlon + c

! Sun's True Anomaly in degrees

true_anom = mean_anom + c

! Sun's Apparent Longitude in degrees

omega = 125.04d0 - 1934.136d0 * t
app_lon = truelon - 0.00569d0 - 0.00478d0 * sin(d2r*omega)

! Eccentricity of Earth Orbit

e = 0.016708634d0 - t * (0.000042037d0 + 0.0000001267d0 * t)

! Mean Obliquity of Ecliptic in degrees

seconds = 21.448d0 - t*(46.8150d0 + t*(0.00059d0 - t*(0.001813d0)))
eps0 = 23.d0 + (26.d0 + (seconds/60.d0)) / 60.d0

! Corrected Obliquity

epsilon = eps0 + 0.00256d0 * cos(d2r*omega)

! Sun's Declination in degrees

sint = sin(d2r*epsilon) * sin(d2r*app_lon)
decl = r2d*asin(sint)

! Equation of Time in minutes

m = mean_anom
lz = meanlon
y = (tan(d2r*epsilon/2.d0))**2
sin2lz = sin(2.d0 * d2r*lz)
sinm   = sin(d2r*m)
cos2lz = cos(2.d0 * d2r*lz)
sin4lz = sin(4.d0 * d2r*lz)
sin2m  = sin(2.d0 * d2r*m)
eqtime = y * sin2lz - 2.d0 * e * sinm &
          + 4.d0 * e * y * sinm * cos2lz &
          - 0.5d0 * y * y * sin4lz - 1.25d0 * e * e * sin2m
eqtime = r2d * eqtime * 4.d0

! Now use DECL and EQTIME to calculate Solar Elevation and Azimuth at (Lat,Lon) 
! for time specified by yr,mo,day,hr,min,sec

! UTC time in Decimal Hours

utc = real(hr,kind=8) + (real(min,kind=8)+sec/60.d0) / 60.d0 - tz

! Hour Angle in degrees, negative is before noon.

ha = 180.d0 - 15.d0*utc - lon - eqtime/4.d0

! Convert to Radians

declr = d2r*decl
latr = d2r*lat
har = d2r*ha

! Solar Elevation Angle

elr = asin(sin(latr)*sin(declr)+cos(latr)*cos(declr)*cos(har))

! Solar Azimuthal Angle
! atan2(y,x) = arctan(y,x) is in range (-pi,pi].

yaz = -sin(har)*cos(declr)
xaz = cos(har)*cos(declr)*sin(latr) - sin(declr)*cos(latr)

! Above formulas are for az measured from due South
! So rotate by 180 to to get az measured from due North

azr = atan2(-yaz,-xaz)
az = r2d* azr
if (az .LT. 0.d0) az = az + 360.d0
el = r2d*elr

solar_angles(1) = az
solar_angles(2) = el

RETURN
end FUNCTION solar_angles

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    view_angles_from_target
!
! PURPOSE
!    Calculate viewing angles of satellite at a given target location P, 
!    subsatellite point SSP, and satellite altitude h.
!
! INPUT
!    latp  = latitude of target in degrees (-90 to 90)
!    lonp  = longitude of target in degrees (-180 to 180)
!    latss = latitude of subsatellite point in degrees (-90 to 90)
!    lonss = longitude of subsatellite point in degrees (-180 to 180)
!    h     = altitude above Earth in km
!
! OUTPUT
!    view_angles is returned array of length 2.
!    view_angles(1) is Satellite Azimuth Angle in Degrees, 0-360, clockwise from North at Target.
!    view_angles(2) is Satellite Elevation Angle in Degrees, 0-90, measured at target between
!    satellite and local horizontal (Angle EPSILON in Reference).
!    Example in Reference: Given LATP,LONP = 22, -160; LATSSP,LONSSP = 10, -175; h=1000 km;
!    view_angles = [232.5, 14.42] NOTE: We need view angles for satellite as viewed from target, 
!    hence results are for this scenario and not that in the reference.
!    
! REFERENCE
!    Space Mission Analysis and Design, 
!        by James R. Wertz, Wiley J. Larson. 
!        Chapter 5, Space Mission Geometry, pp 112-114. 
!        (These pages are available online in
!         http://astrobooks.com/files/SMAD3Err3rd.pdf)
!    See also http://www.aoe.vt.edu/~cdhall/courses/aoe4140/missa.pdf
!    and http://en.wikipedia.org/wiki/Spherical_law_of_cosines
!
! HISTORY
!    written 1 April 2010 by Joyce Wolf for Annemarie Eldering and Susan Kulawik (JPL)
!    updated 30 April 2010 to add warning if target is outside 
!    region visible to spacecraft (spacecraft is below horizon at target).
!    translated into Fortran90 24 Feb 2012 by Vijay Natraj
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    

FUNCTION view_angles_from_target(latp, lonp, latss, lonss, h)

implicit none

!  Inputs

real(kind=8),    intent(in)  :: latp
real(kind=8),    intent(in)  :: lonp
real(kind=8),    intent(in)  :: latss
real(kind=8),    intent(in)  :: lonss
real(kind=8),    intent(in)  :: h

!  Outputs
   
real(kind=8)                 :: view_angles_from_target(2)

!  Local variables

real(kind=8)                 :: d2r, re, srho
real(kind=8)                 :: deltaL, cdeltaL, clatss, slatss, clatp, slatp
real(kind=8)                 :: clambda, slambda, cphiv, phiv
real(kind=8)                 :: taneta, eta, ceps, eps, r2d
real(kind=8)                 :: lambda, lambda0

d2r = acos(-1.d0)/180.d0

! Earth radius in km

re = 6378.d0

! sin(angular radius of Earth at height h)

srho = re/(re+h)

deltaL = abs(lonss-lonp)
cdeltaL = cos(d2r*deltaL)
clatss = cos(d2r*latss)
slatss = sin(d2r*latss)
clatp = cos(d2r*latp)
slatp = sin(d2r*latp)

! Find lambda, central angle of great circle arc connecting P and SSP.
! use Law of Cosines for spherical triangle with vertices P, SSP, and North Pole.
! sides are central angles (great circle arcs) 90-LatP, 90-LatSS, and Lambda.
! Law of Cosines: cos(c) = cos(a) cos(b) +sin(a) sin(b) cos(C),
! where a, b, c are the sides and C is the corner angle opposite side c.

! cos(lambda)

clambda = slatss*slatp + clatss*clatp*cdeltaL

! sin(lambda) 

slambda = sin(acos(clambda))

! cos phiv (Phi_V is azimuth of satellite measured from North at target P).
! Use Law of Cosines on Spherical Triangle formed by P, North Pole, SSP.

cphiv = (slatss - slatp*clambda) / ( clatp*slambda)
if (cphiv .GT. 1.d0) cphiv = 1.d0
if (cphiv .LT. -1.d0) cphiv = -1.d0
phiv = acos(cphiv)

! tan eta

taneta = srho*slambda / (1.d0-srho*clambda)
eta = atan(taneta)

! cos epsilon

ceps = sin(eta)/srho
eps = acos(ceps)

r2d = 180.d0/acos(-1.d0)
view_angles_from_target(1) = r2d*phiv
if (lonp-lonss .GT. 0.d0) view_angles_from_target(1) = 360.d0 - view_angles_from_target(1)
view_angles_from_target(2) = r2d*eps

! Check for spacecraft below horizon at target

lambda = r2d*acos(clambda)
lambda0 = r2d*acos(srho)
if (lambda .GT. lambda0) then
   write(0,*) 'WARNING: SPACECRAFT BELOW HORIZON AT TARGET'
   write(0,*) 'Lambda  (Central Angle between Target and SSP) = ', lambda
   write(0,*) 'Lambda0 (Central Angle Visible to Spacecraft)  = ', lambda0
   view_angles_from_target(2) = -view_angles_from_target(2)
endif

return
end function view_angles_from_target

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    view_angles
!
! PURPOSE
!    Calculate viewing angles of satellite at a given target location P, 
!    subsatellite point SSP, and satellite altitude h.
!
! INPUT
!    latp  = latitude of target in degrees (-90 to 90)
!    lonp  = longitude of target in degrees (-180 to 180)
!    latss = latitude of subsatellite point in degrees (-90 to 90)
!    lonss = longitude of subsatellite point in degrees (-180 to 180)
!    h     = altitude above Earth in km
!
! OUTPUT
!    view_angles is returned array of length 2.
!    view_angles(1) is Satellite Azimuth Angle in Degrees, 0-360, clockwise from North at Target.
!    view_angles(2) is Satellite Azimuth Angle in Degrees, 0-90, measured at satellite between
!    satellite and target (Angle ETA in Reference).
!    Example in Reference: Given LATP,LONP = 22, -200; LATSSP,LONSSP = 10, -185; h=1000 km;
!    view_angles = [48.3, 56.8] 
!    
! REFERENCE
!    Space Mission Analysis and Design, 
!        by James R. Wertz, Wiley J. Larson. 
!        Chapter 5, Space Mission Geometry, pp 112-114. 
!
! HISTORY
!    written 1 April 2010 by Joyce Wolf for Annemarie Eldering and Susan Kulawik (JPL)
!    updated 30 April 2010 to add warning if target is outside 
!    region visible to spacecraft (spacecraft is below horizon at target).
!    translated into Fortran90 24 Feb 2012 by Vijay Natraj
!    adapted for view from spacecraft by P. Castellanos Dec 2019
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    

FUNCTION view_angles_from_space(latp, lonp, latss, lonss, h)

implicit none

!  Inputs

real(kind=8),    intent(in)  :: latp
real(kind=8),    intent(in)  :: lonp
real(kind=8),    intent(in)  :: latss
real(kind=8),    intent(in)  :: lonss
real(kind=8),    intent(in)  :: h

!  Outputs
   
real(kind=8)                 :: view_angles_from_space(2)

!  Local variables

real(kind=8)                 :: d2r, re, srho
real(kind=8)                 :: deltaL, cdeltaL, clatss, slatss, clatp, slatp
real(kind=8)                 :: clambda, slambda, cphiv, phiv
real(kind=8)                 :: taneta, eta, ceps, eps, r2d
real(kind=8)                 :: lambda, lambda0

d2r = acos(-1.d0)/180.d0

! Earth radius in km

re = 6378.d0

! sin(angular radius of Earth at height h)

srho = re/(re+h)

deltaL = abs(lonss-lonp)
cdeltaL = cos(d2r*deltaL)
clatss = cos(d2r*latss)
slatss = sin(d2r*latss)
clatp = cos(d2r*latp)
slatp = sin(d2r*latp)

! Find lambda, central angle of great circle arc connecting P and SSP.
! use Law of Cosines for spherical triangle with vertices P, SSP, and North Pole.
! sides are central angles (great circle arcs) 90-LatP, 90-LatSS, and Lambda.
! Law of Cosines: cos(c) = cos(a) cos(b) +sin(a) sin(b) cos(C),
! where a, b, c are the sides and C is the corner angle opposite side c.

! cos(lambda)

clambda = slatss*slatp + clatss*clatp*cdeltaL

! sin(lambda) 

slambda = sin(acos(clambda))

! cos phiv (Phi_V is azimuth of satellite measured from North at target P).
! Use Law of Cosines on Spherical Triangle formed by P, North Pole, SSP.

cphiv = (slatp - clambda*slatss) / (slambda*clatss)

if (cphiv .GT. 1.d0) cphiv = 1.d0
if (cphiv .LT. -1.d0) cphiv = -1.d0
phiv = acos(cphiv)

! tan eta

taneta = srho*slambda / (1.d0-srho*clambda)
eta = atan(taneta)

! cos epsilon

ceps = sin(eta)/srho
eps = acos(ceps)

r2d = 180.d0/acos(-1.d0)

! azimuth angle
view_angles_from_space(1) = r2d*phiv
! if target is west of satellite, phiv should be negative
! subtract from 360 to get azimuth angles in 0-360
if (lonp-lonss .LT. 0.d0) view_angles_from_space(1) = 360.d0 - view_angles_from_space(1)

! zenith angle
view_angles_from_space(2) = r2d*eta

! Check for spacecraft below horizon at target

lambda = r2d*acos(clambda)
lambda0 = r2d*acos(srho)
if (lambda .GT. lambda0) then
   write(0,*) 'WARNING: SPACECRAFT BELOW HORIZON AT TARGET'
   write(0,*) 'Lambda  (Central Angle between Target and SSP) = ', lambda
   write(0,*) 'Lambda0 (Central Angle Visible to Spacecraft)  = ', lambda0
   view_angles_from_space(2) = -view_angles_from_space(2)
endif

return
end function view_angles_from_space



!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    satellite_angles
!
! PURPOSE
!    For a satellite, at a given target location,
!    calculate spacecraft Azimuth Viewing Angle VAA
!    and spacecraft Zenith Viewing Angle VZA;
!    also given UTC time calculate
!    Solar Zenith Angle (SZA) and Solar Azimuth Angle(SAA).
!    (User may then calculate Relative Azimuth Angle RAA = VAA-SAA)
!
! INPUTS
!    Yr,Mo,Day,Hr,Min,Sec = UTC Time
!    tz    = local time zone (0.d0 = UTC, -8.d0 = PST)
!    Latp  = latitude of target in degrees (-90 to 90)
!    Lonp  = longitude of target in degrees (-180 to 180)
!    Latss = latitude of satellite in degrees (-90 to 90)
!    Lonss = longitude of satellite in degrees (-180 to 180)
!    hgtss = altitude of satellite in km
!
! OUTPUTS
!   satellite_angles    = returned array of length 4
!   satellite_angles(1) = Spacecraft Viewing Azimuth Angle in deg (0-360)
!   satellite_angles(2) = Spacecraft Viewing Zenith Angle in deg (0-90)
!   satellite_angles(3) = Solar Azimuth Angle in degrees (0-360)
!   satellite_angles(4) = Solar Zenith Angles in degrees (0-180, >90 is below horizon)
!
! REQUIRED
!   solar_angles.f90
!   view_angles_from_space.f90
!   julday.f90
!
! HISTORY
!   written 29 April 2010 by Joyce Wolf for Annmarie Eldering & S. Kulawik
!   translated into Fortran90 24 Feb 2012 by Vijay Natraj
!   fixed code to give angles from space not target Dec 2019 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION satellite_angles(yr,mo,day,hr,min,sec,tz,latp,lonp,latss,lonss,hgtss)

implicit none

!  Inputs

integer(kind=4), intent(in) :: yr
integer(kind=4), intent(in) :: mo
integer(kind=4), intent(in) :: day
integer(kind=4), intent(in) :: hr
integer(kind=4), intent(in) :: min
real(kind=8),    intent(in) :: sec
real(kind=8),    intent(in) :: tz
real(kind=8),    intent(in) :: latp  ! latitude of target in degrees (-90 to 90)
real(kind=8),    intent(in) :: lonp  ! longitude of target in degrees (-180 to 180)
real(kind=8),    intent(in) :: latss ! latitude of satellite in degrees (-90 to 90)
real(kind=8),    intent(in) :: lonss ! longitude of satellite in degrees (-180 to 180)
real(kind=8),    intent(in) :: hgtss ! altitude of satellite in km

!  Outputs

real(kind=8) :: satellite_angles(4)
                             ! 1: Spacecraft Viewing Azimuth Angle in deg (0-360)
                             ! 2: Spacecraft Viewing Zenith Angle in deg (0-90)
                             ! 3: Solar Azimuth Angle in degrees (0-360)
                             ! 4: Solar Zenith Angles in degrees (0-180), 
                             !    >90 is below horizon)

!  Local variables

real(kind=8)                :: vangles(2), phiv, thetav
real(kind=8)                :: azel(2), saa, sza

vangles = view_angles_from_space(latp,lonp,latss,lonss,hgtss)
phiv    = vangles(1)
thetav  = vangles(2)

azel = solar_angles(yr,mo,day,hr,min,sec,tz,latp,lonp)
saa = azel(1)
sza = 90.d0 - azel(2)

satellite_angles(1) = phiv
satellite_angles(2) = thetav
satellite_angles(3) = saa
satellite_angles(4) = sza

return
end function satellite_angles

end Module LidarAngles

! -------------------------- Python Interface -------------------

subroutine satAngles(n,yr,mo,day,hr,min,sec,latp,lonp,latss,lonss,tz,hgtss, &
                     view_azimuth, view_zenith,  &
                     solar_azimuth, solar_zenith )

use LidarAngles

implicit none

!  Inputs

integer(kind=4), intent(in) :: n        ! number of points
integer(kind=4), intent(in) :: yr(n)    ! UTC year
integer(kind=4), intent(in) :: mo(n)    ! UTC month
integer(kind=4), intent(in) :: day(n)   ! UTC day
integer(kind=4), intent(in) :: hr(n)    ! UTC hour
integer(kind=4), intent(in) :: min(n)   ! UTC min
real(kind=8),    intent(in) :: sec(n)   ! UTC second
real(kind=8),    intent(in) :: lonp(n)  ! longitude of target in degrees (-180 to 180)
real(kind=8),    intent(in) :: latp(n)  ! latitude of target in degrees (-90 to 90)
real(kind=8),    intent(in) :: lonss(n) ! longitude of satellite in degrees (-180 to 180)
real(kind=8),    intent(in) :: latss(n) ! latitude of satellite in degrees (-90 to 90)
real(kind=8),    intent(in) :: tz       ! local time zone (0.d0 = UTC, -8.d0 = PST)
real(kind=8),    intent(in) :: hgtss    ! altitude of satellite in km



!  Outputs

real(kind=8), intent(out) :: view_azimuth(n)
real(kind=8), intent(out) :: view_zenith(n)
real(kind=8), intent(out) :: solar_azimuth(n)
real(kind=8), intent(out) :: solar_zenith(n)

!                ---

  integer(kind=4) :: i

  real(kind=8) :: sat_angles(4)

  do i = 1, n

     sat_angles = satellite_angles(yr(i),mo(i),day(i),hr(i),min(i),sec(i),tz, &
                                   latp(i),lonp(i),latss(i),lonss(i),hgtss)

     view_azimuth(i)  = sat_angles(1)
     view_zenith(i)   = sat_angles(2)
     solar_azimuth(i) = sat_angles(3)
     solar_zenith(i)  = sat_angles(4)

  end do

end subroutine satAngles

subroutine solarAngles(n,yr,mo,day,hr,min,sec,latp,lonp,tz, &
                     solar_azimuth, solar_zenith )

use LidarAngles

implicit none

!  Inputs

integer(kind=4), intent(in) :: n        ! number of points
integer(kind=4), intent(in) :: yr(n)    ! UTC year
integer(kind=4), intent(in) :: mo(n)    ! UTC month
integer(kind=4), intent(in) :: day(n)   ! UTC day
integer(kind=4), intent(in) :: hr(n)    ! UTC hour
integer(kind=4), intent(in) :: min(n)   ! UTC min
real(kind=8),    intent(in) :: sec(n)   ! UTC second
real(kind=8),    intent(in) :: lonp(n)  ! longitude of target in degrees (-180 to 180)
real(kind=8),    intent(in) :: latp(n)  ! latitude of target in degrees (-90 to 90)
real(kind=8),    intent(in) :: tz       ! local time zone (0.d0 = UTC, -8.d0 = PST)



!  Outputs

real(kind=8), intent(out) :: solar_azimuth(n)
real(kind=8), intent(out) :: solar_zenith(n)

!                ---

  integer(kind=4) :: i

  real(kind=8) :: azel(2), sza, saa

  do i = 1, n

     azel = solar_angles(yr(i),mo(i),day(i),hr(i),min(i),sec(i),tz,latp(i),lonp(i))
     saa = azel(1)
     sza = 90.d0 - azel(2)


     solar_azimuth(i) = saa
     solar_zenith(i)  = sza

  end do

end subroutine solarAngles

