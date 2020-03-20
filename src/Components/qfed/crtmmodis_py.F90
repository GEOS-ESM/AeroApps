!
!  Assumed mapping of CRTM to MODIS channel
!
!  CRTM Channel         MODIS Channel
!                        #       nm
! -------------      ------------------
!      1                 20     3750
!      2*                21     3959  
!      3                 22     3959
!      4                 23     4050
!      5                 24     4465.5
!      6                 25     4525.5
!      7                 27     6715
!      8                 28     7325
!      9                 29     8550
!     10                 30     9730
!     11*                31    11030
!     12                 32    12020
!     13                 33    13335
!     14                 34    13635
!     15                 35    13935
!     16                 36    14235

!      
! (*) MODIS fire channels
!

#define COEFFnccs '/gpfsm/dnb52/projects/p10/gmao_ops/fvInput_4dvar/gsi/etc/fix_ncep20130901/REL-2.1.3/CRTM_Coeffs/Little_Endian/'
#define COEFFcalc '/Users/adasilva/data/CRTM/REL-2.1.3/CRTM_Coeffs/Little_Endian/'

!---
subroutine getSfcTrans ( km, nobs,                      &
                          U, V, T, Qv, O3, DELP,        &
                          Ts, u10m, v10m,               &
                          sensor_zenith_angle,          & ! local zenith angle
                          source_zenith_angle,          & ! solar zenith angle
                          sensor_scan_angle,            & ! scan angle
                          nch, channels, coefPath,      &
                          sfcTrans, rc )

  use CRTM_Module
  use CRTM_Atmosphere_Define, only: H2O_ID, O3_ID, VOLUME_MIXING_RATIO_UNITS, MASS_MIXING_RATIO_UNITS
  implicit NONE
    
  integer, intent(in) :: km
  integer, intent(in) :: nobs
  integer, intent(in) :: nch
  integer, intent(in) :: channels(nch) ! desired channels, CRTM ordering, not MODIS

  real,    intent(in) :: U(km,nobs)
  real,    intent(in) :: V(km,nobs)
  real,    intent(in) :: T(km,nobs)
  real,    intent(in) :: Qv(km,nobs)
  real,    intent(in) :: O3(km,nobs)
  real,    intent(in) :: DELP(km,nobs)
  real,    intent(in) :: Ts(nobs), u10m(nobs), v10m(nobs)

  real,    intent(in) :: sensor_zenith_angle(nobs)  ! local zenith angle
  real,    intent(in) :: source_zenith_angle(nobs)  ! solar zenith angle
  real,    intent(in) :: sensor_scan_angle(nobs)    ! scan angle

!f2py character(len=255), intent(in), optional :: coefPath='/Users/adasilva/data/CRTM/REL-2.1.3/CRTM_Coeffs/Little_Endian/'
  character(len=255), intent(in)  :: coefPath

  real,    intent(out) :: sfcTrans(nch,nobs)

  integer, intent(out) :: rc
  
!                -----------

    type (crtm_channelinfo_type) :: channelInfo(1)
    character(len=32)            :: sensorList(1)

   type (crtm_atmosphere_type)   :: atmosphere(1)
   type (crtm_surface_type)      :: surface(1)
   type (crtm_geometry_type)     :: geometry(1)

   type (crtm_rtsolution_type), allocatable :: rtsolution(:,:)

   integer, parameter ::  n_absorbers = 2
   integer, parameter ::  n_clouds = 0
   integer, parameter ::  n_aerosols = 0

   integer :: n_channels, n, k, i, ic
   real    :: alpha
   real    :: pm(km), pe(km+1)
   real(kind=8) :: optical_depth
   real, parameter    :: ptop = 1.

!  Hardwire these for now
!  ----------------------
   rc = 0
   sfcTrans = 1.
   sensorList(1) = 'modis_terra'
   
!  Initialize CRTM
!  ---------------
   rc = CRTM_Init(sensorList,channelInfo,File_Path=coefPath)
   if ( rc /= 0 ) return

   n_channels = channelInfo(1)%n_channels
   
   allocate(RTSolution(n_channels,1))
   
   call CRTM_Atmosphere_Create(atmosphere,km,n_absorbers,n_clouds,n_aerosols) 
   call CRTM_Surface_Create(surface,n_channels)
   call CRTM_rtSolution_Create(rtSolution,km)
   
! Profile independent settings
! ----------------------------
  atmosphere(1)%n_layers = km
  atmosphere(1)%absorber_id(1) = H2O_ID
  atmosphere(1)%absorber_id(2) = O3_ID
  atmosphere(1)%absorber_units(1) = MASS_MIXING_RATIO_UNITS
  atmosphere(1)%absorber_units(2) = VOLUME_MIXING_RATIO_UNITS
  atmosphere(1)%level_pressure(0) = TOA_PRESSURE

  surface(1)%sensordata%n_channels = channelinfo(1)%n_channels
  surface(1)%sensordata%sensor_id  = channelinfo(1)%sensor_id
  surface(1)%sensordata%wmo_sensor_id = channelinfo(1)%wmo_sensor_id
  surface(1)%sensordata%wmo_satellite_id = channelinfo(1)%wmo_satellite_id
  surface(1)%sensordata%sensor_channel  = channelinfo(1)%sensor_channel
  
  ! Loop over observations (profiles)
  ! ---------------------------------
  do n = 1, nobs

     ! Layer and level pressures
     ! -------------------------
     pe(1) = ptop
     do k = 1, km
        pe(k+1) = pe(k) + delp(k,n)
     end do
     do k = 1, km
        pm(k) = (pe(k)+pe(k+1))/2.
     end do

     ! Set geometry
     ! ------------
     geometry(1)%sensor_zenith_angle = sensor_zenith_angle(n)
     geometry(1)%source_zenith_angle = source_zenith_angle(n)
     geometry(1)%sensor_scan_angle   = sensor_scan_angle(n)
     geometry(1)%ifov                = 0

     ! Set surface
     ! -----------
     surface(1)%sensordata%tb(:)  = ts(n)
     surface(1)%wind_speed        = sqrt(u10m(n)**2+v10m(n)**2)
     surface(1)%water_coverage    = 1.0 ! check if this is OK
     surface(1)%ice_coverage      = 0.0
     surface(1)%land_coverage     = 0.0
     surface(1)%snow_coverage     = 0.0
     surface(1)%water_temperature = max(ts(n),273.)
     surface(1)%land_temperature  = ts(n)
     surface(1)%ice_temperature   = ts(n)
     surface(1)%snow_temperature  = ts(n)
     surface(1)%soil_temperature  = ts(n)

     ! Set Atmosphere
     ! --------------
     atmosphere(1)%level_pressure(1:)    = pe(2:) / 100.
     atmosphere(1)%pressure(:)          = pm / 100.   
     atmosphere(1)%temperature(:)       = T(:,n)
     atmosphere(1)%absorber(:,1)        = Qv(:,n)/(1-Qv(:,n))
     atmosphere(1)%absorber(:,2)        = O3(:,n)
     
     ! Go CRTM, Go!
     ! ------------
     rc = CRTM_Forward ( Atmosphere, Surface, Geometry, ChannelInfo,  &
          RTSolution )
    
     if ( rc /= 0 ) return

     ! Compute surface transmittance
     ! -----------------------------
     do i = 1, nch 
        ic = channels(i) 
        optical_depth = 0.0
        do  k = 1,km 
           optical_depth  = optical_depth + RTSolution(ic,1)%layer_optical_depth(k) 
        end do
        sfcTrans(ic,n) = exp(-min(limit_exp,optical_depth)) ! secant term removed 
     end do

     
  end do

  ! Finalize
  ! --------
  !rc = CRTM_Destroy(channelInfo)
  CALL CRTM_RTSolution_Destroy(RTSolution)
  CALL CRTM_Atmosphere_Destroy(atmosphere)
  CALL CRTM_Surface_Destroy(surface)
  
  deallocate(RTSolution)

end subroutine getSfcTrans
