!
!  Python interface to ods
!

#define INDEF -99999
#define UNDEF 1.0E+15

subroutine pyods_open(filename,fid,rc)
  implicit NONE
  character(len=*), intent(in) :: filename
  integer, intent(out) :: fid,rc
  integer ier
  call ODS_Open ( fid, filename, 'r', rc )
end subroutine pyods_open

subroutine pyods_close(fid,rc)
  implicit NONE
  integer, intent(in) :: fid
  integer, intent(out) :: rc
!---
  call ODS_Close ( fid, ' ', rc )
end subroutine pyods_close

subroutine pyods_nget(fid,nymd,nhms,nobs,rc)
  implicit NONE
  integer, intent(in) :: fid
  integer, intent(in) :: nymd, nhms
  integer, intent(out) :: nobs
  integer, intent(out) :: rc
!---
  integer :: jday, syn_hour
  integer :: ods_julian, ods_nget
  jday = ODS_Julian ( nymd )
  syn_hour = nhms / 10000
  nobs = ODS_NGet ( fid, jday, syn_hour, rc )

end subroutine pyods_nget

subroutine pyods_getInt(fid,name,nymd,nhms,nobs,v,rc)
  implicit NONE
  integer, intent(in) :: fid
  character(len=*) :: name
  integer, intent(in) :: nymd, nhms, nobs
  integer, intent(out) :: v(nobs)
  integer, intent(out) :: rc
!---
  integer :: ier, jday, syn_hour
  integer :: ods_julian, ods_nget
  jday = ODS_Julian ( nymd )
  syn_hour = nhms / 10000
  call ODS_GetI(fid,name,jday,syn_hour,nobs,v,rc)
end subroutine pyods_getInt

subroutine pyods_getFloat(fid,name,nymd,nhms,nobs,v,rc)
  implicit NONE
  integer, intent(in) :: fid
  character(len=*) :: name
  integer, intent(in) :: nymd, nhms, nobs
  real, intent(out) :: v(nobs)
  integer, intent(out) :: rc
!---
  integer :: ier, jday, syn_hour
  integer :: ods_julian, ods_nget
  real :: v8(nobs)
  jday = ODS_Julian ( nymd )
  syn_hour = nhms / 10000
  call ODS_GetR(fid,name,jday,syn_hour,nobs,v8,rc) ! assumes ODS is r8
  v(1:nobs) = v8(1:nobs)
end subroutine pyods_getFloat

subroutine pyods_getAll(filename,nymd,nhms,nobs, &
                        kid, lat, lon, lev, &
                        kx, kt, ks, xm, time, &
                        obs, omf, oma, xvec, qcexcl, qchist, &
                        rc )
  use m_ODS
  implicit NONE
  integer, parameter :: NOBS_MAX = 5000 * 1000
  character(len=*), intent(in) :: filename
  integer, intent(in) :: nymd, nhms

  integer, intent(out) :: nobs             ! number of observations

  integer, intent(out) :: kid(NOBS_MAX)    ! Obs identification index

  real,    intent(out) :: lat(NOBS_MAX)    ! latitute     of obs (degrees)
  real,    intent(out) :: lon(NOBS_MAX)    ! longitude    of obs (degrees)
  real,    intent(out) :: lev(NOBS_MAX)    ! level        of obs (hPa)

  integer, intent(out) :: kx(NOBS_MAX)     ! data source index
  integer, intent(out) :: kt(NOBS_MAX)     ! data type   index
  integer, intent(out) :: ks(NOBS_MAX)     ! sounding    index
  real,    intent(out) :: xm(NOBS_MAX)     ! atomic metadata (depends on kx)
  integer, intent(out) :: time(NOBS_MAX)   ! time (relative to the input/output
                                    !   synoptic date and time)
  real,    intent(out) :: obs(NOBS_MAX)    ! observation value (units depend on kt)
  real,    intent(out) :: OmF(NOBS_MAX)    ! obs minus forecast (O-F) innovations
  real,    intent(out) :: OmA(NOBS_MAX)    ! obs minus analysis (O-A) residuals
  real,    intent(out) :: Xvec(NOBS_MAX)   ! PSAS CG solution vector
  integer, intent(out) :: qcexcl(NOBS_MAX) ! On-line QC exclusion flag
  integer, intent(out) :: qchist(NOBS_MAX) ! On-line QC history flag

  integer, intent(out) :: rc

!                             ----

  character(len=256) :: ftype
  type(ods_vect) :: ods
  
  call ods_get (filename, nymd, nhms, ftype, ods, rc)
  if ( rc .ne. 0 ) return

  if ( nobs > NOBS_MAX ) then
     rc = -1
     return
  end if

  nobs = ods%data%nobs
  kid(1:nobs) = ods%data%kid(1:nobs)
  lat(1:nobs) = ods%data%lat(1:nobs)
  lon(1:nobs) = ods%data%lon(1:nobs)
  lev(1:nobs) = ods%data%lev(1:nobs)
  kx(1:nobs) = ods%data%kx(1:nobs)
  kt(1:nobs) = ods%data%kt(1:nobs)
  ks(1:nobs) = ods%data%ks(1:nobs)
  xm(1:nobs) = ods%data%xm(1:nobs)
  time(1:nobs) = ods%data%time(1:nobs)
  obs(1:nobs) = ods%data%obs(1:nobs)
  omf(1:nobs) = ods%data%omf(1:nobs)
  oma(1:nobs) = ods%data%oma(1:nobs)
  xvec(1:nobs) = ods%data%xvec(1:nobs)
  qcexcl(1:nobs) = ods%data%qcexcl(1:nobs)
  qchist(1:nobs) = ods%data%qchist(1:nobs)

  call ods_clean(ods,rc)

end subroutine pyods_getAll

subroutine pyods_putAll(filename, ftype, nymd,nhms, nsyn, nobs, &
                        kid, lat, lon, lev, &
                        kx, kt, ks, xm, time, &
                        obs, omf, oma, xvec, qcexcl, qchist, &
                        rc )
  use m_ODS
  implicit NONE
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: ftype  ! 'pre_anal' or 'post_anal'
  integer, intent(in)  :: nymd, nhms
  integer, intent(in)  :: nsyn

  integer, intent(in)  :: nobs        ! number of observations

  integer, intent(in) :: kid(nobs)    ! Obs identification index

  real,    intent(in) :: lat(nobs)    ! latitute     of obs (degrees)
  real,    intent(in) :: lon(nobs)    ! longitude    of obs (degrees)
  real,    intent(in) :: lev(nobs)    ! level        of obs (hPa)

  integer, intent(in) :: kx(nobs)     ! data source index
  integer, intent(in) :: kt(nobs)     ! data type   index
  integer, intent(in) :: ks(nobs)     ! sounding    index
  real,    intent(in) :: xm(nobs)     ! atomic metadata (depends on kx)
  integer, intent(in) :: time(nobs)   ! time (relative to the input/output
                                      !   synoptic date and time)
  real,    intent(in) :: obs(nobs)    ! observation value (units depend on kt)
  real,    intent(in) :: OmF(nobs)    ! obs minus forecast (O-F) innovations
  real,    intent(in) :: OmA(nobs)    ! obs minus analysis (O-A) residuals
  real,    intent(in) :: Xvec(nobs)   ! PSAS CG solution vector
  integer, intent(in) :: qcexcl(nobs) ! On-line QC exclusion flag
  integer, intent(in) :: qchist(nobs) ! On-line QC history flag

  integer, intent(out) :: rc

! Note: Metadata not filled.

!                             ----

  integer :: nobs_
  type(ods_vect) :: ods
  
  rc = 0
  nobs_ = nobs

  if ( nobs == 0 ) return

  call ods_init ( ods, nobs_, rc )
  if ( rc /= 0 ) return 

  ods%meta%nsyn = nsyn

  ods%data%kid(1:nobs)    =   kid(1:nobs) 
  ods%data%lat(1:nobs)    =   lat(1:nobs) 
  ods%data%lon(1:nobs)    =   lon(1:nobs) 
  ods%data%lev(1:nobs)    =   lev(1:nobs) 
  ods%data%kx(1:nobs)     =   kx(1:nobs) 
  ods%data%kt(1:nobs)     =   kt(1:nobs) 
  ods%data%ks(1:nobs)     =   ks(1:nobs) 
  ods%data%xm(1:nobs)     =   xm(1:nobs) 
  ods%data%time(1:nobs)   =   time(1:nobs) 
  ods%data%obs(1:nobs)    =   obs(1:nobs) 
  ods%data%omf(1:nobs)    =   omf(1:nobs) 
  ods%data%oma(1:nobs)    =   oma(1:nobs) 
  ods%data%xvec(1:nobs)   =   xvec(1:nobs) 
  ods%data%qcexcl(1:nobs) =   qcexcl(1:nobs) 
  ods%data%qchist(1:nobs) =   qchist(1:nobs) 

!
! Note: kt/kx lists are hardwired for now. Better solution
!       is to have this code fragmented created from kx/kt lists
!       on the fly during build.
!

  ods%meta%kt_names(1) = 'Surface (10m) zonal wind'
  ods%meta%kt_names(2) = 'Surface (10m) meridional wind'
  ods%meta%kt_names(3) = 'Sea level pressure'
  ods%meta%kt_names(4) = 'Upper-air zonal wind'
  ods%meta%kt_names(5) = 'Upper-air meridional wind'
  ods%meta%kt_names(6) = 'Upper-air geopotential height'
  ods%meta%kt_names(7) = 'Upper-air water vapor mixing'
  ods%meta%kt_names(8) = 'Upper-air temperature'
  ods%meta%kt_names(9) = 'Upper-air dew-point temperature'
  ods%meta%kt_names(10) = 'Upper-air relative humidity'
  ods%meta%kt_names(11) = 'Upper-air specific humidity'
  ods%meta%kt_names(12) = 'Surface (10m) wind speed'
  ods%meta%kt_names(13) = 'Surface (10m) temperature Kelvin'
  ods%meta%kt_names(14) = 'Surface (10m) dew-point temperature'
  ods%meta%kt_names(15) = 'Surface (10m) relative humidity'
  ods%meta%kt_names(16) = 'Surface (10m) specific humidity'
  ods%meta%kt_names(17) = 'Precipitation rate'
  ods%meta%kt_names(18) = 'Total precipitable water'
  ods%meta%kt_names(19) = 'Total cloud liquid water'
  ods%meta%kt_names(20) = 'Fractional cloud cover'
  ods%meta%kt_names(21) = 'Column ozone'
  ods%meta%kt_names(22) = 'Ozone layer'
  ods%meta%kt_names(23) = 'Upper-air thickness'
  ods%meta%kt_names(24) = 'Surface (2m) wind speed'
  ods%meta%kt_names(25) = 'Surface (2m) maximum wind'
  ods%meta%kt_names(26) = 'Surface (2m) maximum wind'
  ods%meta%kt_names(27) = 'Surface (2m) temperature'
  ods%meta%kt_names(28) = 'Surface (2m) maximum temperature'
  ods%meta%kt_names(29) = 'Surface (2m) minimum temperature'
  ods%meta%kt_names(30) = 'Surface (2m) dew-point temperature'
  ods%meta%kt_names(31) = 'Surface (2m) relative humidity'
  ods%meta%kt_names(32) = 'Surface (2m) specific humidity'
  ods%meta%kt_names(33) = 'Surface (2m) pressure'
  ods%meta%kt_names(34) = 'Visibility'
  ods%meta%kt_names(35) = 'Snow depth'
  ods%meta%kt_names(36) = 'Weather condition'
  ods%meta%kt_names(37) = 'Upper-air potential temperature (ps=1000hPa)'
  ods%meta%kt_names(38) = 'Surface skin temperature'
  ods%meta%kt_names(39) = 'Sea Surface temperature'
  ods%meta%kt_names(40) = 'Brightness temperature'
  ods%meta%kt_names(41) = 'Layer-mean virtual temperature'
  ods%meta%kt_names(42) = 'Line-of-Sight Wind (LOS)'
  ods%meta%kt_names(43) = 'Log-transformed AOD'
  ods%meta%kt_names(44) = 'Upper-air virtual temperature'
  ods%meta%kt_names(45) = 'Aerosol Optical Depth'
  ods%meta%kt_names(46) = 'Angstrom Exponent'
  ods%meta%kt_names(47) = 'Surface Albedo'
  ods%meta%kt_names(48) = 'Reflectance'
  ods%meta%kt_names(49) = 'Aerosol Absorption Optical Depth'
  ods%meta%kt_names(50) = 'Single Scattering Albedo'
  ods%meta%kt_names(51) = 'CO mixing ratio'
  ods%meta%kt_names(52) = 'NO2 mixing ratio'
  ods%meta%kt_names(53) = 'unknown'
  ods%meta%kt_names(54) = 'Solar Zenith Angle'
  ods%meta%kt_names(55) = 'Solar Azimuth Angle'
  ods%meta%kt_names(56) = 'Sensor Zenith Angle'
  ods%meta%kt_names(57) = 'Sensor Azimuth Angle'
  ods%meta%kt_names(58) = 'Aerosol Scattering Angle'
  ods%meta%kt_names(59) = 'unknown'
  ods%meta%kt_names(60) = 'Aerosol Type'
  ods%meta%kt_names(61) = 'Fine model Optical thickness'
  ods%meta%kt_names(62) = 'Fine model Angstrom exponent'
  ods%meta%kt_names(63) = 'unknown'
  ods%meta%kt_names(64) = 'unknown'
  ods%meta%kt_names(65) = 'unknown'
  ods%meta%kt_names(66) = 'unknown'
  ods%meta%kt_names(67) = 'unknown'
  ods%meta%kt_names(68) = 'unknown'
  ods%meta%kt_names(69) = 'unknown'
  ods%meta%kt_names(70) = 'Aerosol Optical Depth Ratio'
  ods%meta%kt_names(71) = 'unknown'
  ods%meta%kt_names(72) = 'unknown'
  ods%meta%kt_names(73) = 'unknown'
  ods%meta%kt_names(74) = 'Quality Assurance Flag'
  ods%meta%kt_names(75) = 'unknown'
  ods%meta%kt_names(76) = 'Aerosol Refractive index (Real)'
  ods%meta%kt_names(77) = 'Aerosol Refractive index (Imaginary)'
  ods%meta%kt_names(78) = 'Reflectivity'
  ods%meta%kt_names(79) = 'Height of aerosol layer'
  ods%meta%kt_names(80) = 'Radiance'
  ods%meta%kt_names(81) = 'Aerosol Index'
  ods%meta%kt_names(82) = 'Fire Radiative Power Megawatts'
  ods%meta%kt_names(83) = 'Aerosol Asymmetry Factor'
  ods%meta%kt_names(84) = 'Aerosol Mass Concentration'
  ods%meta%kt_names(85) = 'Aerosol Effective Radius'
  ods%meta%kt_names(86) = 'GSI precip rate ln(1+rain_rate)'
  ods%meta%kt_names(87) = 'Ozone mixing ratio'

  ods%meta%kx_names(1) = 'Surface Land Obs - 1  '
  ods%meta%kx_names(2) = 'Surface Land Obs - 2  '
  ods%meta%kx_names(3) = 'Surface Ship Obs - 1  '
  ods%meta%kx_names(4) = 'Surface Ship Obs - 2  '
  ods%meta%kx_names(5) = 'Environment Buoy     '
  ods%meta%kx_names(6) = 'Drifting Buoy     '
  ods%meta%kx_names(7) = 'Rawinsonde      '
  ods%meta%kx_names(8) = 'Pilot Wind     '
  ods%meta%kx_names(9) = 'Ship Released Rawinsonde    '
  ods%meta%kx_names(10) = 'Dropwinsonde      '
  ods%meta%kx_names(11) = 'Radar-tracked Rawinsonde     '
  ods%meta%kx_names(12) = 'Rocketsonde      '
  ods%meta%kx_names(13) = 'Balloon      '
  ods%meta%kx_names(14) = 'Aircraft - Air/Sat Relay   '
  ods%meta%kx_names(15) = 'Aircraft - Int. Data Sys  '
  ods%meta%kx_names(16) = 'Aircraft Report     '
  ods%meta%kx_names(17) = 'Aircraft Coded Report    '
  ods%meta%kx_names(18) = 'Aircraft - ALPX    '
  ods%meta%kx_names(19) = 'Cloud Track Wind - Wisc E1 '
  ods%meta%kx_names(20) = 'Cloud Track Wind - Wisc E2 '
  ods%meta%kx_names(21) = 'Cloud Track Wind - Wisc W '
  ods%meta%kx_names(22) = 'Cloud Track Wind - Wisc Ocean '
  ods%meta%kx_names(23) = 'Cloud Track Wind - Repr. Jap. '
  ods%meta%kx_names(24) = 'Cloud Track Wind - US1 Picture triplet'
  ods%meta%kx_names(25) = 'Cloud Track Wind - US2 Picture triplet'
  ods%meta%kx_names(26) = 'Cloud Track Wind - European IR automated'
  ods%meta%kx_names(27) = 'Cloud Track Wind - Japanese IR automated'
  ods%meta%kx_names(28) = 'SCATTEROMETER Research Mode    '
  ods%meta%kx_names(29) = 'WSAT 55     '
  ods%meta%kx_names(30) = 'WSAT 57     '
  ods%meta%kx_names(31) = 'Limb Infrared Monitor - Strat  '
  ods%meta%kx_names(32) = 'Aircraft EU experiment    '
  ods%meta%kx_names(33) = 'NESDIS NH Land AM type A '
  ods%meta%kx_names(34) = 'NESDIS SH Land AM type A '
  ods%meta%kx_names(35) = 'NESDIS NH Land AM type B /'
  ods%meta%kx_names(36) = 'NESDIS SH Land AM type B /'
  ods%meta%kx_names(37) = 'NESDIS NH Land AM type C '
  ods%meta%kx_names(38) = 'NESDIS SH Land AM type C '
  ods%meta%kx_names(39) = 'NESDIS NH Ocean AM type A '
  ods%meta%kx_names(40) = 'NESDIS SH Ocean AM type A '
  ods%meta%kx_names(41) = 'NESDIS NH Ocean AM type B /'
  ods%meta%kx_names(42) = 'NESDIS SH Ocean AM type B /'
  ods%meta%kx_names(43) = 'NESDIS NH Ocean AM type C '
  ods%meta%kx_names(44) = 'NESDIS SH Ocean AM type C '
  ods%meta%kx_names(45) = 'NESDIS NH Land PM type A '
  ods%meta%kx_names(46) = 'NESDIS SH Land PM type A '
  ods%meta%kx_names(47) = 'NESDIS NH Land PM type B /'
  ods%meta%kx_names(48) = 'NESDIS SH Land PM type B /'
  ods%meta%kx_names(49) = 'NESDIS NH Land PM type C '
  ods%meta%kx_names(50) = 'NESDIS SH Land PM type C '
  ods%meta%kx_names(51) = 'NESDIS NH Ocean PM type A '
  ods%meta%kx_names(52) = 'NESDIS SH Ocean PM type A '
  ods%meta%kx_names(53) = 'NESDIS NH Ocean PM type B /'
  ods%meta%kx_names(54) = 'NESDIS SH Ocean PM type B /'
  ods%meta%kx_names(55) = 'NESDIS NH Ocean PM type C '
  ods%meta%kx_names(56) = 'NESDIS SH Ocean PM type C '
  ods%meta%kx_names(57) = 'Special Sat NH - Ocn A '
  ods%meta%kx_names(58) = 'Special Sat SH - Ocn A '
  ods%meta%kx_names(59) = 'Special Sat NH - Ocn B '
  ods%meta%kx_names(60) = 'Special Sat SH - Ocn B '
  ods%meta%kx_names(61) = 'Special Sat NH - Ocn C '
  ods%meta%kx_names(62) = 'Special Sat SH - Ocn C '
  ods%meta%kx_names(63) = 'VAS NH Land - type A '
  ods%meta%kx_names(64) = 'VAS SH Land - type A '
  ods%meta%kx_names(65) = 'VAS NH Land - type B '
  ods%meta%kx_names(66) = 'VAS SH Land - type B '
  ods%meta%kx_names(67) = 'VAS NH Ocean - type A '
  ods%meta%kx_names(68) = 'VAS SH Ocean - type A '
  ods%meta%kx_names(69) = 'VAS NH Ocean - type B '
  ods%meta%kx_names(70) = 'VAS SH Ocean - type B '
  ods%meta%kx_names(71) = 'NASA-GLA NH Land type A  '
  ods%meta%kx_names(72) = 'NASA-GLA SH Land type A  '
  ods%meta%kx_names(73) = 'NASA-GLA NH Land type B  '
  ods%meta%kx_names(74) = 'NASA-GLA SH Land type B  '
  ods%meta%kx_names(75) = 'NASA-GLA NH Land type C  '
  ods%meta%kx_names(76) = 'NASA-GLA SH Land type C  '
  ods%meta%kx_names(77) = 'NASA-GLA NH Land type D  '
  ods%meta%kx_names(78) = 'NASA-GLA SH Land type D  '
  ods%meta%kx_names(79) = 'NASA-GLA NH Oceantype A   '
  ods%meta%kx_names(80) = 'NASA-GLA SH Ocean type A  '
  ods%meta%kx_names(81) = 'NASA-GLA NH Ocean type B  '
  ods%meta%kx_names(82) = 'NASA-GLA SH Ocean type B  '
  ods%meta%kx_names(83) = 'NASA-GLA NH Ocean type C  '
  ods%meta%kx_names(84) = 'NASA-GLA SH Ocean type C  '
  ods%meta%kx_names(85) = 'NASA-GLA NH Ocean type D  '
  ods%meta%kx_names(86) = 'NASA-GLA SH Ocean type D  '
  ods%meta%kx_names(87) = 'Pseudo-1000 mb Heights    '
  ods%meta%kx_names(88) = 'ER-2 Aircraft / MMS Data  '
  ods%meta%kx_names(89) = 'Aircraft reports (ACARS)    '
  ods%meta%kx_names(90) = 'Surface METAR     '
  ods%meta%kx_names(91) = 'UARS / MLS 1   '
  ods%meta%kx_names(92) = 'UARS / MLS 2   '
  ods%meta%kx_names(93) = 'DAOTOVS Land AM Type: Clr HIRS/MSU/SSU '
  ods%meta%kx_names(94) = 'DAOTOVS Land AM Type: Clr HIRS/MSU '
  ods%meta%kx_names(95) = 'DAOTOVS Land AM Type: MSU/SSU  '
  ods%meta%kx_names(96) = 'DAOTOVS Land AM Type: MSU  '
  ods%meta%kx_names(97) = 'DAOTOVS Land AM Type: SSU  '
  ods%meta%kx_names(98) = 'DAOTOVS Ocean AM Type: Clr HIRS/MSU/SSU '
  ods%meta%kx_names(99) = 'DAOTOVS Ocean AM Type: Clr HIRS/MSU '
  ods%meta%kx_names(100) = 'DAOTOVS Ocean AM Type: MSU/SSU  '
  ods%meta%kx_names(101) = 'DAOTOVS Ocean AM Type: MSU  '
  ods%meta%kx_names(102) = 'DAOTOVS Ocean AM Type: SSU  '
  ods%meta%kx_names(103) = 'DAOTOVS Ice AM Type: Clr HIRS/MSU/SSU '
  ods%meta%kx_names(104) = 'DAOTOVS Ice AM Type: Clr HIRS/MSU '
  ods%meta%kx_names(105) = 'DAOTOVS Ice AM Type: MSU/SSU  '
  ods%meta%kx_names(106) = 'DAOTOVS Ice AM Type: MSU  '
  ods%meta%kx_names(107) = 'DAOTOVS Ice AM Type: SSU  '
  ods%meta%kx_names(108) = 'DAOTOVS Land AM Type: Cld Clr HIRS/MSU/SSU'
  ods%meta%kx_names(109) = 'DAOTOVS Land AM Type: Cld Clr HIRS/MSU'
  ods%meta%kx_names(110) = 'DAOTOVS Ocean AM Type: Cld Clr HIRS/MSU/SSU'
  ods%meta%kx_names(111) = 'DAOTOVS Ocean AM Type: Cld Clr HIRS/MSU'
  ods%meta%kx_names(112) = 'DAOTOVS Ice AM Type: Cld Clr HIRS/MSU/SSU'
  ods%meta%kx_names(113) = 'DAOTOVS Ice AM Type: Cld Clr HIRS/MSU'
  ods%meta%kx_names(114) = 'Cloud Track Wind - INDIAN OCEAN '
  ods%meta%kx_names(115) = 'Cloud Track Wind - GOES 9 (04/04/97-03/11/98)'
  ods%meta%kx_names(116) = 'SSM/I NESDIS Precipitation    '
  ods%meta%kx_names(117) = 'TOMS      '
  ods%meta%kx_names(118) = 'SBUV/2      '
  ods%meta%kx_names(119) = 'Cloud Track Wind - US1 IR automated'
  ods%meta%kx_names(120) = 'Cloud Track Wind - US1 Water Vapor,Deep-Layer'
  ods%meta%kx_names(121) = 'Cloud Track Wind - US1 Water Vapor,Cloud-Top'
  ods%meta%kx_names(122) = 'Cloud Track Wind - US2 IR automated'
  ods%meta%kx_names(123) = 'Cloud Track Wind - US2 Water Vapor,Deep-Layer'
  ods%meta%kx_names(124) = 'Cloud Track Wind - US2 Water Vapor,Cloud-Top'
  ods%meta%kx_names(125) = 'DAOTOVS Land PM Type: Clr HIRS/MSU/SSU '
  ods%meta%kx_names(126) = 'DAOTOVS Land PM Type: Clr HIRS/MSU '
  ods%meta%kx_names(127) = 'DAOTOVS Land PM Type: MSU/SSU  '
  ods%meta%kx_names(128) = 'DAOTOVS Land PM Type: MSU  '
  ods%meta%kx_names(129) = 'DAOTOVS Land PM Type: SSU  '
  ods%meta%kx_names(130) = 'DAOTOVS Ocean PM Type: Clr HIRS/MSU/SSU '
  ods%meta%kx_names(131) = 'DAOTOVS Ocean PM Type: Clr HIRS/MSU '
  ods%meta%kx_names(132) = 'DAOTOVS Ocean PM Type: MSU/SSU  '
  ods%meta%kx_names(133) = 'DAOTOVS Ocean PM Type: MSU  '
  ods%meta%kx_names(134) = 'DAOTOVS Ocean PM Type: SSU  '
  ods%meta%kx_names(135) = 'DAOTOVS Ice PM Type: Clr HIRS/MSU/SSU '
  ods%meta%kx_names(136) = 'DAOTOVS Ice PM Type: Clr HIRS/MSU '
  ods%meta%kx_names(137) = 'DAOTOVS Ice PM Type: MSU/SSU  '
  ods%meta%kx_names(138) = 'DAOTOVS Ice PM Type: MSU  '
  ods%meta%kx_names(139) = 'DAOTOVS Ice PM Type: SSU  '
  ods%meta%kx_names(140) = 'DAOTOVS Land PM Type: Cld Clr HIRS/MSU/SSU'
  ods%meta%kx_names(141) = 'DAOTOVS Land PM Type: Cld Clr HIRS/MSU'
  ods%meta%kx_names(142) = 'DAOTOVS Ocean PM Type: Cld Clr HIRS/MSU/SSU'
  ods%meta%kx_names(143) = 'DAOTOVS Ocean PM Type: Cld Clr HIRS/MSU'
  ods%meta%kx_names(144) = 'DAOTOVS Ice PM Type: Cld Clr HIRS/MSU/SSU'
  ods%meta%kx_names(145) = 'DAOTOVS Ice PM Type: Cld Clr HIRS/MSU'
  ods%meta%kx_names(146) = 'Cloud Track Wind - US1 visible automated'
  ods%meta%kx_names(147) = 'Cloud Track Wind - US2 visible automated'
  ods%meta%kx_names(148) = 'SSM/I WENTZ Speed only   '
  ods%meta%kx_names(149) = 'SSM/I NESDIS Speed only   '
  ods%meta%kx_names(150) = 'SSM/I WENTZ F8/F14    '
  ods%meta%kx_names(151) = 'SCATTEROMETER ERS1     '
  ods%meta%kx_names(152) = 'SCATTEROMETER ERS2     '
  ods%meta%kx_names(153) = 'SCATTEROMETER NSCAT     '
  ods%meta%kx_names(154) = 'SCATTEROMETER QuikSCAT     '
  ods%meta%kx_names(155) = 'SCATTEROMETER SEAWINDS     '
  ods%meta%kx_names(156) = 'TRMM      '
  ods%meta%kx_names(157) = 'TMI      '
  ods%meta%kx_names(158) = 'GADS(BRITISH AIRWAYS)     '
  ods%meta%kx_names(159) = 'SAT am, WATER, CLFRC<=.10   '
  ods%meta%kx_names(160) = 'SAT am, WATER, CLFRC<=.40   '
  ods%meta%kx_names(161) = 'SAT am, WATER, CLFRC> .40  '
  ods%meta%kx_names(162) = 'SAT am, Other, CLFRC<=.10   '
  ods%meta%kx_names(163) = 'SAT am, Other, CLFRC<=.40   '
  ods%meta%kx_names(164) = 'SAT am, Other, CLFRC> .40  '
  ods%meta%kx_names(165) = 'SAT pm, WATER, CLFRC<=.10   '
  ods%meta%kx_names(166) = 'SAT pm, WATER, CLFRC<=.40   '
  ods%meta%kx_names(167) = 'SAT pm, WATER, CLFRC> .40  '
  ods%meta%kx_names(168) = 'SAT pm, Other, CLFRC<=.10   '
  ods%meta%kx_names(169) = 'SAT pm, Other, CLFRC<=.40   '
  ods%meta%kx_names(170) = 'SAT pm, Other, CLFRC> .40  '
  ods%meta%kx_names(171) = 'GPS DAO Refractive Retrieval   '
  ods%meta%kx_names(172) = 'GPS DAO Bending Angle Retrieval  '
  ods%meta%kx_names(173) = 'GPS DAO Observation Retrieval   '
  ods%meta%kx_names(174) = 'GPS JPL Refractive Retrieval   '
  ods%meta%kx_names(175) = 'GPS JPL Bending Angle Retrieval  '
  ods%meta%kx_names(176) = 'GPS JPL Observation Retrieval   '
  ods%meta%kx_names(177) = 'GPS NCAR Refractive Retrieval   '
  ods%meta%kx_names(178) = 'GPS NCAR Bending Angle Retrieval  '
  ods%meta%kx_names(179) = 'GPS NCAR Observation Retrieval   '
  ods%meta%kx_names(180) = 'GPS INST Refractive Retrieval   '
  ods%meta%kx_names(181) = 'GPS INST Bending Angle Retrieval  '
  ods%meta%kx_names(182) = 'GPS INST Observation Retrieval   '
  ods%meta%kx_names(183) = 'SSM/I WENTZ F10/F15    '
  ods%meta%kx_names(184) = 'SSM/I WENTZ F11    '
  ods%meta%kx_names(185) = 'SSM/I WENTZ F13    '
  ods%meta%kx_names(186) = 'DAOTOVS Land AM Type: Clr HIRS/AMSUA/AMSUB '
  ods%meta%kx_names(187) = 'DAOTOVS Land AM Type: Clr HIRS/AMSUA '
  ods%meta%kx_names(188) = 'DAOTOVS Land AM Type: AMSUA/AMSUB  '
  ods%meta%kx_names(189) = 'DAOTOVS Land AM Type: AMSUA  '
  ods%meta%kx_names(190) = 'DAOTOVS Land AM Type: AMSUB  '
  ods%meta%kx_names(191) = 'DAOTOVS Ocean AM Type: Clr HIRS/AMSUA/AMSUB '
  ods%meta%kx_names(192) = 'DAOTOVS Ocean AM Type: Clr HIRS/AMSUA '
  ods%meta%kx_names(193) = 'DAOTOVS Ocean AM Type: AMSUA/AMSUB  '
  ods%meta%kx_names(194) = 'DAOTOVS Ocean AM Type: AMSUA  '
  ods%meta%kx_names(195) = 'DAOTOVS Ocean AM Type: AMSUB  '
  ods%meta%kx_names(196) = 'DAOTOVS Ice AM Type: Clr HIRS/AMSUA/AMSUB '
  ods%meta%kx_names(197) = 'DAOTOVS Ice AM Type: Clr HIRS/AMSUA '
  ods%meta%kx_names(198) = 'DAOTOVS Ice AM Type: AMSUA/AMSUB  '
  ods%meta%kx_names(199) = 'DAOTOVS Ice AM Type: AMSUA  '
  ods%meta%kx_names(200) = 'DAOTOVS Ice AM Type: AMSUB  '
  ods%meta%kx_names(201) = 'DAOTOVS Land AM Type: Cld Clr HIRS/AMSUA/AMSUB'
  ods%meta%kx_names(202) = 'DAOTOVS Land AM Type: Cld Clr HIRS/AMSUA'
  ods%meta%kx_names(203) = 'DAOTOVS Ocean AM Type: Cld Clr HIRS/AMSUA/AMSUB'
  ods%meta%kx_names(204) = 'DAOTOVS Ocean AM Type: Cld Clr HIRS/AMSUA'
  ods%meta%kx_names(205) = 'DAOTOVS Ice AM Type: Cld Clr HIRS/AMSUA/AMSUB'
  ods%meta%kx_names(206) = 'DAOTOVS Ice AM Type: Cld Clr HIRS/AMSUA'
  ods%meta%kx_names(207) = 'DAOTOVS Land PM Type: Clr HIRS/AMSUA/AMSUB '
  ods%meta%kx_names(208) = 'DAOTOVS Land PM Type: Clr HIRS/AMSUA '
  ods%meta%kx_names(209) = 'DAOTOVS Land PM Type: AMSUA/AMSUB  '
  ods%meta%kx_names(210) = 'DAOTOVS Land PM Type: AMSUA  '
  ods%meta%kx_names(211) = 'DAOTOVS Land PM Type: AMSUB  '
  ods%meta%kx_names(212) = 'DAOTOVS Ocean PM Type: Clr HIRS/AMSUA/AMSUB '
  ods%meta%kx_names(213) = 'DAOTOVS Ocean PM Type: Clr HIRS/AMSUA '
  ods%meta%kx_names(214) = 'DAOTOVS Ocean PM Type: AMSUA/AMSUB  '
  ods%meta%kx_names(215) = 'DAOTOVS Ocean PM Type: AMSUA  '
  ods%meta%kx_names(216) = 'DAOTOVS Ocean PM Type: AMSUB  '
  ods%meta%kx_names(217) = 'DAOTOVS Ice PM Type: Clr HIRS/AMSUA/AMSUB '
  ods%meta%kx_names(218) = 'DAOTOVS Ice PM Type: Clr HIRS/AMSUA '
  ods%meta%kx_names(219) = 'DAOTOVS Ice PM Type: AMSUA/AMSUB  '
  ods%meta%kx_names(220) = 'DAOTOVS Ice PM Type: AMSUA  '
  ods%meta%kx_names(221) = 'DAOTOVS Ice PM Type: AMSUB  '
  ods%meta%kx_names(222) = 'DAOTOVS Land PM Type: Cld Clr HIRS/AMSUA/AMSUB'
  ods%meta%kx_names(223) = 'DAOTOVS Land PM Type: Cld Clr HIRS/AMSUA'
  ods%meta%kx_names(224) = 'DAOTOVS Ocean PM Type: Cld Clr HIRS/AMSUA/AMSUB'
  ods%meta%kx_names(225) = 'DAOTOVS Ocean PM Type: Cld Clr HIRS/AMSUA'
  ods%meta%kx_names(226) = 'DAOTOVS Ice PM Type: Cld Clr HIRS/AMSUA/AMSUB'
  ods%meta%kx_names(227) = 'DAOTOVS Ice PM Type: Cld Clr HIRS/AMSUA'
  ods%meta%kx_names(228) = 'NESTAT AIRS Land PM type A '
  ods%meta%kx_names(229) = 'NESTAT AIRS Land PM type B '
  ods%meta%kx_names(230) = 'NESTAT AIRS Land PM type C '
  ods%meta%kx_names(231) = 'NESTAT AIRS Ocean PM type A '
  ods%meta%kx_names(232) = 'NESTAT AIRS Ocean PM type B '
  ods%meta%kx_names(233) = 'NESTAT AIRS Ocean PM type C '
  ods%meta%kx_names(234) = 'NESPHY AIRS Land PM type A '
  ods%meta%kx_names(235) = 'NESPHY AIRS Land PM type B '
  ods%meta%kx_names(236) = 'NESPHY AIRS Land PM type C '
  ods%meta%kx_names(237) = 'NESPHY AIRS Ocean PM type A '
  ods%meta%kx_names(238) = 'NESPHY AIRS Ocean PM type B '
  ods%meta%kx_names(239) = 'NESPHY AIRS Ocean PM type C '
  ods%meta%kx_names(240) = 'DAOAIRS Land Type: Clr AIRS/AMSU/HSB  '
  ods%meta%kx_names(241) = 'DAOAIRS Land Type: Clr AIRS/AMSU  '
  ods%meta%kx_names(242) = 'DAOAIRS Land Type: Clr AIRS/HSB  '
  ods%meta%kx_names(243) = 'DAOAIRS Land Type: AMSU/HSB   '
  ods%meta%kx_names(244) = 'DAOAIRS Land Type: AMSU   '
  ods%meta%kx_names(245) = 'DAOAIRS Land Type: HSB   '
  ods%meta%kx_names(246) = 'DAOAIRS Ocean Type: Clr AIRS/AMSU/HSB  '
  ods%meta%kx_names(247) = 'DAOAIRS Ocean Type: Clr AIRS/AMSU  '
  ods%meta%kx_names(248) = 'DAOAIRS Ocean Type: Clr AIRS/HSB  '
  ods%meta%kx_names(249) = 'DAOAIRS Ocean Type: AMSU/HSB   '
  ods%meta%kx_names(250) = 'DAOAIRS Ocean Type: AMSU   '
  ods%meta%kx_names(251) = 'DAOAIRS Ocean Type: HSB   '
  ods%meta%kx_names(252) = 'DAOAIRS Ice Type: Clr AIRS/AMSU/HSB  '
  ods%meta%kx_names(253) = 'DAOAIRS Ice Type: Clr AIRS/AMSU  '
  ods%meta%kx_names(254) = 'DAOAIRS Ice Type: Clr AIRS/HSB  '
  ods%meta%kx_names(255) = 'DAOAIRS Ice Type: AMSU/HSB   '
  ods%meta%kx_names(256) = 'DAOAIRS Ice Type: AMSU   '
  ods%meta%kx_names(257) = 'DAOAIRS Ice Type: HSB   '
  ods%meta%kx_names(258) = 'DAOAIRS Land Type: Cld Clr AIRS/AMSU/HSB '
  ods%meta%kx_names(259) = 'DAOAIRS Land Type: Cld Clr AIRS/AMSU '
  ods%meta%kx_names(260) = 'DAOAIRS Land Type: Cld Clr AIRS/HSB '
  ods%meta%kx_names(261) = 'DAOAIRS Ocean Type: Cld Clr AIRS/AMSU/HSB '
  ods%meta%kx_names(262) = 'DAOAIRS Ocean Type: Cld Clr AIRS/AMSU '
  ods%meta%kx_names(263) = 'DAOAIRS Ocean Type: Cld Clr AIRS/HSB '
  ods%meta%kx_names(264) = 'DAOAIRS Ice Type: Cld Clr AIRS/AMSU/HSB '
  ods%meta%kx_names(265) = 'DAOAIRS Ice Type: Cld Clr AIRS/AMSU '
  ods%meta%kx_names(266) = 'DAOAIRS Ice Type: Cld Clr AIRS/HSB '
  ods%meta%kx_names(267) = 'Cloud Track Wind - European Water Vapor'
  ods%meta%kx_names(268) = 'Cloud Track Wind - European visible automated'
  ods%meta%kx_names(269) = 'Cloud Track Wind - Japanese Water Vapor'
  ods%meta%kx_names(270) = 'Cloud Track Wind - Japanese visible automated'
  ods%meta%kx_names(271) = 'Microwave Limb Sounder (MLS) at 205 GHz'
  ods%meta%kx_names(272) = 'Microwave Limb Sounder (MLS) at 183 GHz'
  ods%meta%kx_names(273) = 'Global Ozone Monitoring Experiment (GOME)  '
  ods%meta%kx_names(274) = 'Cloud Track Wind - MISR  '
  ods%meta%kx_names(275) = 'Cloud Track Wind - European ELW***) IR'
  ods%meta%kx_names(276) = 'Cloud Track Wind - European ELW***) visible'
  ods%meta%kx_names(277) = 'Cloud Track Wind - European ELW***) Water'
  ods%meta%kx_names(278) = 'Cloud Track Wind - European High-resolution Visible'
  ods%meta%kx_names(279) = 'Cloud Track Wind - European Clear Sky'
  ods%meta%kx_names(280) = 'TIDE GAUGE STATION    '
  ods%meta%kx_names(281) = 'Cloud Track Wind - TERRA Modis Water'
  ods%meta%kx_names(282) = 'Cloud Track Wind - TERRA Modis IR-CIMSS'
  ods%meta%kx_names(283) = 'Super-pressure balloon     '
  ods%meta%kx_names(284) = 'VAD (NEXRAD) WINDS    '
  ods%meta%kx_names(285) = 'WIND PROFILER     '
  ods%meta%kx_names(286) = 'Cloud Track Wind-European High-resolution Water Vapor Wind(HWW)'
  ods%meta%kx_names(287) = 'GOES IRET - US1 Land Type: Clr'
  ods%meta%kx_names(288) = 'GOES IRET - US1 Ocean Type: Clr'
  ods%meta%kx_names(289) = 'GOES IRET - US2 Land Type: Clr'
  ods%meta%kx_names(290) = 'GOES IRET - US2 Ocean Type: Clr'
  ods%meta%kx_names(291) = 'GOES IRET - US3 Land Type: Clr'
  ods%meta%kx_names(292) = 'GOES IRET - US3 Ocean Type: Clr'
  ods%meta%kx_names(293) = 'HRDI - High Resolution Doppler Imager '
  ods%meta%kx_names(294) = 'SABER - Sounding of the Atmosphere Broadband'
  ods%meta%kx_names(295) = 'Cloud Track Wind - AQUA Modis Water'
  ods%meta%kx_names(296) = 'Cloud Track Wind - AQUA Modis IR-CIMSS'
  ods%meta%kx_names(297) = 'Cloud Track Wind - TERRA Modis Water'
  ods%meta%kx_names(298) = 'Cloud Track Wind - TERRA Modis IR'
  ods%meta%kx_names(299) = 'Cloud Track Wind - AQUA Modis Water'
  ods%meta%kx_names(300) = 'Cloud Track Wind - AQUA Modis IR'
  ods%meta%kx_names(301) = 'TERRA MODIS Aerosol (Ocean Algorithm)  '
  ods%meta%kx_names(302) = 'TERRA MODIS Aerosol (Dark Target Land Algorithm)  '
  ods%meta%kx_names(303) = 'EPTOMS Aerosol (Ocean Algorithm)   '
  ods%meta%kx_names(304) = 'EPTOMS Aerosol (Land Algorithm)   '
  ods%meta%kx_names(305) = 'MSG Cloud track wind - Infrared channel'
  ods%meta%kx_names(306) = 'MSG Cloud track wind - Water vapor'
  ods%meta%kx_names(307) = 'MSG Cloud track wind - Visible channel'
  ods%meta%kx_names(308) = 'MSG Clear sky Water Vapor channel wind'
  ods%meta%kx_names(309) = 'TERRA MODIS Aerosol (Deep Blue Ocean Algorithm)'
  ods%meta%kx_names(310) = 'TERRA MODIS Aerosol (Deep Blue Land Algorithm)'
  ods%meta%kx_names(311) = 'AQUA MODIS Aerosol (Ocean Algorithm)  '
  ods%meta%kx_names(312) = 'AQUA MODIS Aerosol (Dark Target Land Algorithm)  '
  ods%meta%kx_names(313) = 'MISR (Multi-angle Imaging SpectroRadiometer)   '
  ods%meta%kx_names(314) = 'OMI (Ozone Monitoring Instrument)   '
  ods%meta%kx_names(315) = 'Mozaic Aircraft Data    '
  ods%meta%kx_names(316) = 'Parasol (Ocean Algorithm)    '
  ods%meta%kx_names(317) = 'Parasol (Land Algorithm)    '
  ods%meta%kx_names(318) = 'MOPITT (Measurements Of Pollution In The Troposphere)'
  ods%meta%kx_names(319) = 'AQUA MODIS Aerosol (Deep Blue Ocean Algorithm)'
  ods%meta%kx_names(320) = 'AQUA MODIS Blue Aerosol (Deep Blue Land Algorithm)'
  ods%meta%kx_names(321) = 'TERRA MODIS pixFire    '
  ods%meta%kx_names(322) = 'AQUA MODIS Pixfire    '
  ods%meta%kx_names(323) = 'AERONET'
  ods%meta%kx_names(324) = 'AVHRR PATMOSX Aerosol Retrievals'
  ods%meta%kx_names(325) = 'OMSO2 (OMI Sulfer Dioxide)'
  ods%meta%kx_names(326) = 'GOCI Aerosol Retrievals (Yonsei University)'

  call ods_put (filename, ftype, nymd, nhms, ods, rc)
  if ( rc .ne. 0 ) return

  call ods_clean(ods,rc)

end subroutine pyods_putAll
