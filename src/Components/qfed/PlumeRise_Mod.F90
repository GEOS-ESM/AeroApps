!
! This is the the main Plume Rise code base without any ESMF code. It is
! essentially the same code I got from Saulo, except that it has been
! made more object oriented by eliminating static data segments and
! making it thread safe.
!
! Arlindo da Silva (October, 2006)
!..........................................................................

Module PlumeRise_Mod

implicit none

Private
Public PlumeRise
Public PlumeRise_Initialize
Public PlumeRise_Run
Public PlumeRise_Finalize

integer, public, parameter :: N_BIOME = 4
integer, public, parameter :: TROPICAL_FOREST = 1
integer, public, parameter :: XTROPICAL_FOREST = 2
integer, public, parameter :: SAVANNA = 3
integer, public, parameter :: GRASSLAND = 4

character(len=15), parameter :: veg_name(N_BIOME) = (/ &
                               'Tropical Forest', &
                               'Boreal Forest  ', &
                               'Savanna        ', &
                               'Grassland      ' /)

! Hardwired problem size for now
! ------------------------------
  integer, public, parameter :: nkp  = 200, ntime = 200

! The PlumeRise Object
! --------------------
  type PlumeRise

   private

   real,dimension(nkp) ::  w,t,qv,qc,qh,qi,sc,  &  ! blob
                           vth,vti,rho,txs,     &
                           est,qsat,qpas,qtotal
   
   real,dimension(nkp) ::  wc,wt,tt,qvt,qct,qht,qit,sct
   real,dimension(nkp) ::  dzm,dzt,zm,zt,vctr1,vctr2,    &
                           vt3dc,vt3df,vt3dk,vt3dg,scr1

!
                                                    ! environment at plume grid
   real,dimension(nkp) ::  pke,the,thve,thee,pe,te,qvenv,rhe,dne,sce 
                                                    ! environment at RAMS  grid

!-ams: required by srf mods of 04may2009
   real,dimension(nkp) :: vel_p, rad_p, rad_t, vel_t, upe, vpe, vel_e

                          

   real,dimension(nkp) ::  ucon,vcon,wcon,thtcon ,rvcon,picon,tmpcon, &
                           dncon,prcon,zcon,zzcon,scon 

   real :: DZ,DQSDZ,VISC(nkp),VISCOSITY,TSTPF   
   integer :: N,NM1,L
   
   real :: ADVW,ADVT,ADVV,ADVC,ADVH,ADVI,CVH(nkp),CVI(nkp),ADIABAT, &
           WBAR,ALAST(10),VHREL,VIREL  ! advection

   real :: ZSURF,ZBASE,ZTOP
   integer :: LBASE

   real :: AREA,RSURF,ALPHA,RADIUS(nkp)  ! entrain

   real :: HEATING(ntime),FMOIST,BLOAD   ! heating

   real :: DT,TIME,TDUR
   integer :: MINTIME,MDUR,MAXTIME

   real :: ztop_(ntime)

end type PlumeRise

interface PlumeRise_Run
   module procedure PlumeRise_Run_1d 
   module procedure PlumeRise_Run_HF
   module procedure PlumeRise_Run_VMD
end interface



CONTAINS

!..........................................................................

!                        --------------
!                        Public Methods
!                        --------------

subroutine PlumeRise_Initialize ( This )
   use rconstants
   implicit none
#  include "PlumeRise.h"   
   call zero_plumegen_coms(this)
   call set_grid(this) ! define vertical grid of plume model
                       ! zt(k) =  thermo and water levels
                       ! zm(k) =  dynamical levels 
 end subroutine PlumeRise_Initialize

subroutine PlumeRise_Finalize ( This )
   use rconstants
   implicit none
#  include "PlumeRise.h"   
 end subroutine PlumeRise_Finalize


subroutine PlumeRise_Run_1d ( this, km,     &
                              u_, v_, T_, q_, delp_, ptop, areas,   &
                              z1_plume, z2_plume, dz_plume, &
                              p1_plume, p2_plume, dp_plume, &
                              k1_plume, k2_plume)


   use rconstants
   implicit none

   integer, intent(in) :: km                            ! number of vertical layers 
   real,    intent(in) :: ptop                          ! pressure edge top [Pa]
   real,    intent(in) :: u_(km)                        ! Zonal wind [m/s]
   real,    intent(in) :: v_(km)                        ! Meridional wind [m/s]
   real,    intent(in) :: T_(km)                        ! Dry temperature [K]
   real,    intent(in) :: q_(km)                        ! Specific humidity [kg/kg]
   real,    intent(in) :: delp_(km)                     ! layer pressure thickness [Pa]
   real,    intent(in) :: areas(N_BIOME)                ! fire areas [m^2]

   real,    OPTIONAL, intent(out) :: z1_plume(N_BIOME)  ! lower plume height
   real,    OPTIONAL, intent(out) :: z2_plume(N_BIOME)  ! upper plume height
   real,    OPTIONAL, intent(out) :: dz_plume(N_BIOME)  ! z2-z1; z2>z1
   real,    OPTIONAL, intent(out) :: p1_plume(N_BIOME)  ! upper plume pressure
   real,    OPTIONAL, intent(out) :: p2_plume(N_BIOME)  ! lower plume pressure
   real,    OPTIONAL, intent(out) :: dp_plume(N_BIOME)  ! p2 - p1; p2>p1
   integer, OPTIONAL, intent(out) :: k1_plume(N_BIOME)  ! lower vertical index
   integer, OPTIONAL, intent(out) :: k2_plume(N_BIOME)  ! upper vertical index

!                                    -------


!  Derived quantities
!  ------------------
   real  :: theta(km)           ! potential temperature
   real  :: airdens(km)         ! Dry air density
   real  :: mixr(km)            ! mixing ratio
   real  :: pe_(km+1)            ! pressure at edges
   real  :: pm(km)              ! pressure at mid-layers
   real  :: he(km+1)            ! height (m) of edges
   real  :: hm(km)              ! height (m) of mid-layers


  integer :: m1,m2,m3,ibcon,mynum,i,j,k,i_biome,imm,k1,k2,ixx

  integer :: kmt, kkmax

  real                    :: dz_flam, rhodzi, fire_area
  real,    dimension(2)   :: ztopmax, ptopmax
  real,    dimension(2)   :: g5_ztopmax, g5_ptopmax

#include "PlumeRise.h"   
                         ! only : ucon,vcon,wcon,thtcon,rvcon,picon,tmpcon&
                         ! ,dncon,prcon,zcon,zzcon,scon,zero_plumegen_coms

!----------------------------------------------------------------------

! Compute derived quantities
! --------------------------
  call Derived_( theta, airdens, mixr, pe_, pm, he, hm, &
                 km, T_, q_, delp_, ptop )


! Copy input sounding into internal data structure.
! Vertical levels n GEOS-5 are from top to bottom, while the PlumeRise
! model expect them bottom to top: flip levels here
! --------------------------------------------------------------------
  m1 = km + 1
  thtcon(1:km) = theta(km:1:-1)
  dncon(1:km) = airdens(km:1:-1)
  prcon(1:km) = pm(km:1:-1)
  picon(1:km) = cp * (pm(km:1:-1)/p00)**rocp
  rvcon(1:km) = mixr(km:1:-1)
  zcon(1:km) = hm(km:1:-1) - he(km+1)
  zzcon(1:m1) = he(km+1:1:-1) - he(km+1)
  ucon(1:km) = u_(km:1:-1)
  vcon(1:km) = v_(km:1:-1)

! Get envinronmental state (temp, water vapor mix ratio, ...)
! -----------------------------------------------------------
  call get_env_condition(this,m1-1,kmt)
 	    
! Loop over major biomes with fires; see enumeration above
! TROPICAL_FOREST, XTROPICAL_FOREST, SAVANNA, GRASSLAND
! -------------------------------------------------------
  do i_biome=1,N_BIOME
 
     fire_area = areas(i_biome)

!    No fire, nothing to do
!    ----------------------
     if ( fire_area == 0.0 ) then
        if (present(z1_plume))  z1_plume(i_biome) = 0.0
        if (present(z2_plume))  z2_plume(i_biome) = 0.0
        if (present(dz_plume))  dz_plume(i_biome) = 0.0
        if (present(p1_plume))  p1_plume(i_biome) = 0.0
        if (present(p2_plume))  p2_plume(i_biome) = 0.0
        if (present(dp_plume))  dp_plume(i_biome) = 0.0
        if (present(k1_plume))  k1_plume(i_biome) = 0
        if (present(k2_plume))  k2_plume(i_biome) = 0
        cycle
     end if

!    Loop at the maximum and minimum heat flux for this biome
!    --------------------------------------------------------
     do imm=1,2
	      
!       Get fire properties (burned area, plume radius, heating rates ...)
!       ------------------------------------------------------------------
        call get_fire_properties(this,imm,i_biome,fire_area)
     
!       Special handling for Grassland - ask Saulo: WTF is this?
!       --------------------------------------------------------
        if(i_biome == GRASSLAND  .and. imm == 2) then 
           ztopmax(2)=ztopmax(1)
           ztopmax(1)=zzcon(1)
           cycle
        endif
	     
!       ----------------------
        ixx = i_biome*10 + imm
!       ----------------------           

!       Calculate the plume
!       -------------------
        call makeplume (this,kmt,ztopmax(imm),ixx,kkmax,imm) 

!       ZTOPMAX if the final height (above local surface) level
!       of the plume for each heat flux  		 
	   
!!!!    ptopmax(imm) = pe_(km-kkmax+1) ! recall, pe_ is top to bottom
        
     end do ! loop over min/max heat flux (biome is fixed)
	   
!    Define the vertical domain where the flaming emission will be released 
!    for this biome (forest, savanna, grassland)
!    NOTE: k1>k2 are top->bottom GEOS-5 *layer* indices
!    ----------------------------------------------------------------------
     call set_geos5_vert(this,ztopmax,he,km,k1,k2)

!    Bracketing pressure edges based on these (k1,k2)
!    At this point we should have k1 > k2
!    ------------------------------------------------
     g5_ptopmax(1) = pe_(k1)
     g5_ptopmax(2) = pe_(k2+1)
     g5_ztopmax(1) = he(k1)
     g5_ztopmax(2) = he(k2+1)

!    Fill in output
!    --------------
     if (present(z1_plume)) z1_plume(i_biome) = g5_ztopmax(1)
     if (present(z2_plume)) z2_plume(i_biome) = g5_ztopmax(2)
     if (present(dz_plume)) dz_plume(i_biome) = (g5_ztopmax(2)-g5_ztopmax(1))
     if (present(p1_plume)) p1_plume(i_biome) = g5_ptopmax(1)
     if (present(p2_plume)) p2_plume(i_biome) = g5_ptopmax(2) 
     if (present(dp_plume)) dp_plume(i_biome) = (g5_ptopmax(1)-g5_ptopmax(2)) 
     if (present(k1_plume)) k1_plume(i_biome) = k1 
     if (present(k2_plume)) k2_plume(i_biome) = k2 
     
  enddo ! loop over biome

end subroutine PlumeRise_Run_1d

!..........................................................................

subroutine PlumeRise_Run_HF ( this, km,     &
                              u_, v_, T_, q_, delp_, ptop, &
                              fire_area, heat_flux_kW,   &
                              z_plume, p_plume, k_plume )

   use rconstants
   implicit none

   integer, intent(in) :: km                            ! number of vertical layers 
   real,    intent(in) :: ptop                          ! pressure edge top [Pa]
   real,    intent(in) :: u_(km)                        ! Zonal wind [m/s]
   real,    intent(in) :: v_(km)                        ! Meridional wind [m/s]
   real,    intent(in) :: T_(km)                        ! Dry temperature [K]
   real,    intent(in) :: q_(km)                        ! Specific humidity [kg/kg]
   real,    intent(in) :: delp_(km)                     ! layer pressure thickness [Pa]

   real,    intent(in) :: fire_area                     ! fire areas [m^2]
   real,    intent(in) :: heat_flux_kW                  ! Heat Flux [kW/m^2]

   real,    OPTIONAL, intent(out) :: z_plume           ! top plume height
   real,    OPTIONAL, intent(out) :: p_plume           ! top plume pressure
   integer, OPTIONAL, intent(out) :: k_plume           ! top vertical index

!                                    -------

!  Derived quantities
!  ------------------
   real  :: theta(km)           ! potential temperature
   real  :: airdens(km)         ! Dry air density
   real  :: mixr(km)            ! mixing ratio
   real  :: pe_(km+1)            ! pressure at edges
   real  :: pm(km)              ! pressure at mid-layers
   real  :: he(km+1)            ! height (m) of edges
   real  :: hm(km)              ! height (m) of mid-layers


  integer :: m1,m2,m3,ibcon,mynum,i,j,k,i_biome,imm,k1,k2,ixx

  integer :: kmt, kkmax

  real                    :: dz_flam, rhodzi
  real,    dimension(2)   :: ztopmax, ptopmax
  real,    dimension(2)   :: g5_ztopmax, g5_ptopmax

#include "PlumeRise.h"   
                         ! only : ucon,vcon,wcon,thtcon,rvcon,picon,tmpcon&
                         ! ,dncon,prcon,zcon,zzcon,scon,zero_plumegen_coms

!----------------------------------------------------------------------

! Compute derived quantities
! --------------------------
  call Derived_( theta, airdens, mixr, pe_, pm, he, hm, &
                 km, T_, q_, delp_, ptop )


! Copy input sounding into internal data structure.
! Vertical levels n GEOS-5 are from top to bottom, while the PlumeRise
! model expect them bottom to top: flip levels here
! --------------------------------------------------------------------
  m1 = km + 1
  thtcon(1:km) = theta(km:1:-1)
  dncon(1:km) = airdens(km:1:-1)
  prcon(1:km) = pm(km:1:-1)
  picon(1:km) = cp * (pm(km:1:-1)/p00)**rocp
  rvcon(1:km) = mixr(km:1:-1)
  zcon(1:km) = hm(km:1:-1) - he(km+1)
  zzcon(1:m1) = he(km+1:1:-1) - he(km+1)
  ucon(1:km) = u_(km:1:-1)
  vcon(1:km) = v_(km:1:-1)

! Get envinronmental state (temp, water vapor mix ratio, ...)
! -----------------------------------------------------------
  call get_env_condition(this,m1-1,kmt)
 	    
! No fire, nothing to do
! ----------------------
  if ( fire_area == 0.0 ) then
     if (present(z_plume))  z_plume = 0.0
     if (present(p_plume))  p_plume = 0.0
     if (present(k_plume))  k_plume = 0
     return
  end if
	      
! Get fire properties (burned area, plume radius, heating rates ...)
! ------------------------------------------------------------------
  call get_fire_properties_HF(this,heat_flux_kW,fire_area)
     
! ----------------------
  i_biome = 1 ! not relevant here
  imm = 1     ! not relevant here
  ixx = i_biome*10 + imm
! ----------------------           

! Calculate the plume
! -------------------
  call makeplume (this,kmt,ztopmax(imm),ixx,kkmax,imm) 

! Define the vertical domain where the flaming emission will be released 
! for this biome (forest, savanna, grassland)
! NOTE: k1>k2 are top->bottom GEOS-5 *layer* indices
! ----------------------------------------------------------------------
  ztopmax(2)=ztopmax(1)
  call set_geos5_vert(this,ztopmax,he,km,k1,k2)

! Bracketing pressure edges based on these (k1,k2)
! At this point we should have k1 > k2
! ------------------------------------------------
  g5_ptopmax(1) = pe_(k1)
  g5_ptopmax(2) = pe_(k2+1)
  g5_ztopmax(1) = he(k1)
  g5_ztopmax(2) = he(k2+1)

! Fill in output
! --------------
  if (present(z_plume)) z_plume = ztopmax(1)
!!!  if (present(z_plume)) z_plume = g5_ztopmax(1)
  if (present(p_plume)) p_plume = g5_ptopmax(1)
  if (present(k_plume)) k_plume = k1 
     
end subroutine PlumeRise_Run_HF

!..........................................................................

subroutine PlumeRise_Run_VMD ( this, km,     &
                              u_, v_, T_, q_, delp_, ptop, &
                              fire_area, heat_flux_kW,   &
                              z_i, z_d, z_a, z_f, &
                              z_plume, w_plume )

   use rconstants
   implicit none

   integer, intent(in) :: km                            ! number of vertical layers 
   real,    intent(in) :: ptop                          ! pressure edge top [Pa]
   real,    intent(in) :: u_(km)                        ! Zonal wind [m/s]
   real,    intent(in) :: v_(km)                        ! Meridional wind [m/s]
   real,    intent(in) :: T_(km)                        ! Dry temperature [K]
   real,    intent(in) :: q_(km)                        ! Specific humidity [kg/kg]
   real,    intent(in) :: delp_(km)                     ! layer pressure thickness [Pa]

   real,    intent(in) :: fire_area                     ! fire areas [m^2]
   real,    intent(in) :: heat_flux_kW                  ! Heat Flux [kW/m^2]

   real, intent(out) :: z_i    ! level of maximum W
   real, intent(out) :: z_d    ! level of maximum detrainment
   real, intent(out) :: z_a    ! average height in (z_i,z_f), weighted by -dw/dz
   real, intent(out) :: z_f    ! level when w<1.

   real, intent(out) :: z_plume(nkp) ! vertical levels
   real, intent(out) :: w_plume(nkp) ! vertical velocity


!                                    -------

!  Derived quantities
!  ------------------
   real  :: theta(km)           ! potential temperature
   real  :: airdens(km)         ! Dry air density
   real  :: mixr(km)            ! mixing ratio
   real  :: pe_(km+1)            ! pressure at edges
   real  :: pm(km)              ! pressure at mid-layers
   real  :: he(km+1)            ! height (m) of edges
   real  :: hm(km)              ! height (m) of mid-layers

  integer :: m1,m2,m3,ibcon,mynum,i,j,k,i_biome,imm,k1,k2,ixx
  integer :: kmt, kkmax
  real,  dimension(2) :: ztopmax


#include "PlumeRise.h"   
                         ! only : ucon,vcon,wcon,thtcon,rvcon,picon,tmpcon&
                         ! ,dncon,prcon,zcon,zzcon,scon,zero_plumegen_coms

!----------------------------------------------------------------------

! Compute derived quantities
! --------------------------
  call Derived_( theta, airdens, mixr, pe_, pm, he, hm, &
                 km, T_, q_, delp_, ptop )


! Copy input sounding into internal data structure.
! Vertical levels n GEOS-5 are from top to bottom, while the PlumeRise
! model expect them bottom to top: flip levels here
! --------------------------------------------------------------------
  m1 = km + 1
  thtcon(1:km) = theta(km:1:-1)
  dncon(1:km) = airdens(km:1:-1)
  prcon(1:km) = pm(km:1:-1)
  picon(1:km) = cp * (pm(km:1:-1)/p00)**rocp
  rvcon(1:km) = mixr(km:1:-1)
  zcon(1:km) = hm(km:1:-1) - he(km+1)
  zzcon(1:m1) = he(km+1:1:-1) - he(km+1)
  ucon(1:km) = u_(km:1:-1)
  vcon(1:km) = v_(km:1:-1)

! Get envinronmental state (temp, water vapor mix ratio, ...)
! -----------------------------------------------------------
  call get_env_condition(this,m1-1,kmt)
 	    
! No fire, nothing to do
! ----------------------
  if ( fire_area == 0.0 ) then
     z_i = 0.0
     z_d = 0.0
     z_a = 0.0
     z_f = 0.0
     return
  end if
	      
! Get fire properties (burned area, plume radius, heating rates ...)
! ------------------------------------------------------------------
  call get_fire_properties_HF(this,heat_flux_kW,fire_area)
     
! ----------------------
  i_biome = 1 ! not relevant here
  imm = 1     ! not relevant here
  ixx = i_biome*10 + imm
! ----------------------           

! Calculate the plume
! -------------------
  call makeplume (this,kmt,ztopmax(imm),ixx,kkmax,imm) 

! Get plume extent parameters
! ---------------------------
  call get_VMD(z_i,z_d,z_a,z_f,this%w,this%zm,nkp)

! Return plume profile as well
! ----------------------------
  z_plume = this%zm
  w_plume = this%w

end subroutine PlumeRise_Run_VMD

!...................................................................................

  subroutine Derived_( theta, airdens, mixr, pe, pm, he, hm, &
                       km, T, q, delp_, ptop )

   use rconstants
   implicit NONE

   integer, intent(in) :: km                 ! No. vertical layers
   real,    intent(in) :: T(km)              ! Dry temperature [K]
   real,    intent(in) :: q(km)              ! Specific humidity [kg/kg]
   real,    intent(in) :: delp_(km)          ! layer pressure thicjness [Pa]
   real,    intent(in) :: ptop               ! top (edge) pressure [Pa]

   real, intent(out)  :: theta(km)           ! potential temperature
   real, intent(out)  :: airdens(km)         ! Dry air density
   real, intent(out)  :: mixr(km)            ! mixing ratio

   real, intent(out)  :: pe(km+1)            ! pressure at edges
   real, intent(out)  :: pm(km)              ! pressure at mid-layers

   real, intent(out)  :: he(km+1)            ! height at edges
   real, intent(out)  :: hm(km)              ! height at mid-layers

!
! ------------------------------- pe(k),   he(k)
!
! ............................... pm(k),   hm(k), delp(k)
!
! ------------------------------- pe(k+1), he(k+1)
!

   integer :: k
   real :: Tv, kappa = rocp
   real :: delh(km), delhm(km), pek(km+1)

!  Mixing ratio from specific humidity
!  -----------------------------------
   mixr = q / ( 1.0 - q )

!  Construct edge pressures
!  ------------------------
   pe(1) = ptop
   do k = 1, km
      pe(k+1) = pe(k) + delp_(k)
   end do

!  Construct mid-layer pressures and layer thickness
!  -------------------------------------------------
   do k = 1, km
      Tv = T(k) * ( 1 + 0.61 * mixr(k) )
      pm(k)    = 0.5 * ( pe(k) + pe(k+1) )
      delh(k)  = Rgas * Tv * log(pe(k+1)/pe(k)) / g
      delhm(k) = Rgas * Tv * log(pe(k+1)/pm(k)) / g
   end do

!  Compute Geo-potential height at edges and mid-layer
!  ---------------------------------------------------
   he(km+1) = 0.0  ! set the surface to zero height
   do k = km, 1, -1
      he(k) = he(k+1) + delh(k)
      hm(k) = he(k+1) + delhm(k)
   end do

   airdens = delp_ / ( g * delh )

!  Compute theta
!  -------------
   pek = (pe/p00) ** kappa
   do k = 1, km
      theta(k) =  (g/cp) * delh(k) / ( pek(k+1) - pek(k) ) 
   end do

 end subroutine Derived_

!..........................................................................

!                        ---------------
!                        Private Methods
!                        ---------------

real FUNCTION ESAT_PR (TEM)  
!
! ******* Vapor Pressure  A.L. Buck JAM V.20 p.1527. (1981) ***********
!
real, PARAMETER :: CI1 = 6.1115, CI2 = 22.542, CI3 = 273.48
real, PARAMETER :: CW1 = 6.1121, CW2 = 18.729, CW3 = 257.87, CW4 = 227.3
real, PARAMETER :: TMELT = 273.3

real temc , tem,esatm
!
!     formulae from Buck, A.L., JAM 20,1527-1532
!     custom takes esat wrt water always. formula for h2o only
!     good to -40C so:
!
!
TEMC = TEM - TMELT  
IF (TEMC.GT. - 40.0) GOTO 230  
ESATM = CI1 * EXP (CI2 * TEMC / (TEMC + CI3) )  !ice, millibars  
ESAT_PR = ESATM / 10.	!kPa			  

RETURN  
!
230 ESATM = CW1 * EXP ( ( (CW2 - (TEMC / CW4) ) * TEMC) / (TEMC + CW3))
                          
ESAT_PR = ESATM / 10.	!kPa			  
RETURN  
END function ESAT_PR
!---------------------------------------------------------------------------

  subroutine zero_plumegen_coms ( this )
  
# include "PlumeRise.h"  

w=0.0;t=0.0;qv=0.0;qc=0.0;qh=0.0;qi=0.0;sc=0.0
vth=0.0;vti=0.0;rho=0.0;txs=0.0
est=0.0;qsat=0.0;qpas=0.0;qtotal=0.0
wc=0.0;wt=0.0;tt=0.0;qvt=0.0;qct=0.0;qht=0.0;qit=0.0;sct=0.0
dzm=0.0;dzt=0.0;zm=0.0;zt=0.0;vctr1=0.0;vctr2=0.0
vt3dc=0.0;vt3df=0.0;vt3dk=0.0;vt3dg=0.0;scr1=0.0
pke=0.0;the=0.0;thve=0.0;thee=0.0;pe=0.0;te=0.0;qvenv=0.0;rhe=0.0;dne=0.0;sce=0.0 
ucon=0.0;vcon=0.0;wcon=0.0;thtcon =0.0;rvcon=0.0;picon=0.0;tmpcon=0.0;dncon=0.0;prcon=0.0 
zcon=0.0;zzcon=0.0;scon=0.0 
dz=0.0;dqsdz=0.0;visc=0.0;viscosity=0.0;tstpf=0.0
advw=0.0;advt=0.0;advv=0.0;advc=0.0;advh=0.0;advi=0.0;cvh=0.0;cvi=0.0;adiabat=0.0
wbar=0.0;alast=0.0;vhrel=0.0;virel=0.0  
zsurf=0.0;zbase=0.0;ztop=0.0;area=0.0;rsurf=0.0;alpha=0.0;radius=0.0;heating=0.0
fmoist=0.0;bload=0.0;dt=0.0;time=0.0;tdur=0.0
ztop_=0.0

n=0;nm1=0;l=0;lbase=0;mintime=0;mdur=0;maxtime=0

!-srf-04may2009-AWE
vel_p=0.;rad_p=0.;rad_t=0.;vel_t=0.
upe=0.;vpe=0.;vel_e=0.


end subroutine zero_plumegen_coms
!...........................................................................

subroutine get_env_condition(this,k2,kmt)

use rconstants
implicit none
integer :: k1,k2,k,kcon,klcl,kmt,nk,nkmid,i
real :: znz,themax,tlll,plll,rlll,zlll,dzdd,dzlll,tlcl,plcl,dzlcl,dummy
!ams integer :: n_setgrid = 0 
#include "PlumeRise.h"

#ifdef DEBUG
!  Sample output over amazon
!  -------------------------
   print *
   print *, '          Inside PlumeRise'
   print *, '    H(m)  p (hPa) THETA(K)  rho    q(g/kg)'
   do k = 1, k2
      print 10, zcon(k), prcon(k)/100., thtcon(k), dncon(k),& 
                1000*rvcon(k)
   end do
10 format(1x,5f8.2)
#endif

k1=2 !first level of thermodynamic grid point

#if 0
!ams -------------------------------------------------------------
 The call set_grid() has been moved to PlumeRise_Initialize where
 it belongs.

if( n_setgrid == 0) then
  n_setgrid = 1
  call set_grid(this) ! define vertical grid of plume model
                      ! zt(k) =  thermo and water levels
                      ! zm(k) =  dynamical levels 
!ams -------------------------------------------------------------

#endif

znz=zcon(k2)
kmt = -1
do k=nkp,1,-1
  if(zt(k).lt.znz) then
     kmt=k
     exit
  endif
enddo

!ams: Saulo, before you had kmt=k outside the loop --- this
!     does not work in general because the do-loop counter k
!     is usually undefined after a loop.

if ( kmt == -1 ) then
   print *, 'get_env_condition: could not find kmt, aborting...'
   call exit(1)
endif

nk=k2-k1+1
!call htint(nk, wcon,zzcon(k1),kmt,wpe,zt)

!-srf-04may2009-AWE
 call htint(nk,  ucon,zcon(k1:),kmt,upe,zt)
 call htint(nk,  vcon,zcon(k1:),kmt,vpe,zt)
!-srf-end
 
 call htint(nk,thtcon,zcon(k1:),kmt,the  ,zt)
 call htint(nk, rvcon,zcon(k1:),kmt,qvenv,zt)
do k=1,kmt
  qvenv(k)=max(qvenv(k),1e-8)
enddo

pke(1)=picon(1)
do k=1,kmt
  thve(k)=the(k)*(1.+.61*qvenv(k)) ! virtual pot temperature
enddo
do k=2,kmt
  pke(k)=pke(k-1)-g*2.*(zt(k)-zt(k-1))  & ! exner function
        /(thve(k)+thve(k-1))
enddo
do k=1,kmt
  te(k)  = the(k)*pke(k)/cp         ! temperature (K) 
  pe(k)  = (pke(k)/cp)**cpor*p00    ! pressure (Pa)
  dne(k)= pe(k)/(rgas*te(k)*(1.+.61*qvenv(k))) !  dry air density (kg/m3)

!-srf-04may2009-AWE
  vel_e(k) = sqrt(upe(k)**2+vpe(k)**2)
!-srf-end

enddo
!-srf-04may2009-AWE
! to turn off the env. wind effect - do this at compile time
!!!#ifndef WIND_EFF
!!!vel_e(:) = 0.
!!!#endif
!-srf-end


!--------- converte press de Pa para kPa para uso modelo de plumerise
do k=1,kmt
 pe(k) = pe(k)*1.e-3
enddo 

return 

#ifdef TESTING

! only for testing ---------------------------

themax=0.
nkmid=20
kcon=2
do k=2,nkmid
  if(thee(k).gt.themax)then
    themax=thee(k)
    kcon=k
  endif
enddo

kcon=max(2,kcon)
tlll=(te(kcon)+te(kcon+1)+te(kcon-1))/3.
plll=pe(kcon)*1.e+3 ! converte de volta para Pa
rlll=(qvenv(kcon)+qvenv(kcon+1)+qvenv(kcon-1))/3.
zlll=zt(kcon)

call lcl(tlll,plll,rlll,tlcl,plcl,dzlcl)

!  FIND THE CLOSEST LEVEL ON THE CONVECTIVE GRID TO THE LCL

dzlll=1e20
do k=1,kmt
  dzdd=abs(zt(k)-(zlll+dzlcl))
  if(dzdd.lt.dzlll)then
    dzlll=dzdd
    klcl=k
  endif
enddo
zbase=zt(klcl)

!print*,' ---------- LCL --------------'
!print*,kmt,klcl,zt(klcl)

#endif

end subroutine get_env_condition
!-------------------------------------------------------------------------

subroutine set_grid(this)
implicit none
integer :: k,mzp
#include "PlumeRise.h"  

dz=100. ! set constant grid spacing of plume grid model(meters)

mzp=nkp
zt(1) = zsurf
zm(1) = zsurf
zt(2) = zt(1) + 0.5*dz
zm(2) = zm(1) + dz
do k=3,mzp
 zt(k) = zt(k-1) + dz ! thermo and water levels
 zm(k) = zm(k-1) + dz ! dynamical levels	
enddo
!print*,zsurf
!Print*,zt(:)
do k = 1,mzp-1
   dzm(k) = 1. / (zt(k+1) - zt(k))
enddo 
dzm(mzp)=dzm(mzp-1)

do k = 2,mzp
   dzt(k) = 1. / (zm(k) - zm(k-1))
enddo
dzt(1) = dzt(2) * dzt(2) / dzt(3)
   
!   dzm(1) = 0.5/dz
!   dzm(2:mzp) = 1./dz
return
end subroutine set_grid
!-------------------------------------------------------------------------
 
subroutine set_flam_vert(this,ztopmax,k1,k2)
  implicit none
  integer imm,k,k1,k2,kk
  real,    dimension(2)  :: ztopmax
  integer, dimension(2)  :: k_lim
                        ! , only : nkp,zzcon
# include "PlumeRise.h"   

  do imm=1,2
! checar 
!    do k=1,m1-1
    do k=1,nkp-1
       kk = k
              if(zzcon(k) > ztopmax(imm) ) exit
    enddo
    k_lim(imm) = kk
  enddo           
  k1=max(3,k_lim(1))
  k2=max(3,k_lim(2))
  
  if(k2 < k1) then
    !print*,'1: ztopmax k=',ztopmax(1), k1
    !print*,'2: ztopmax k=',ztopmax(2), k2
    k2=k1
    !stop 1234
  endif
    
end subroutine set_flam_vert

!-------------------------------------------------------------------------

subroutine set_geos5_vert(this,ztopmax,he,km,k1,k2)
  implicit none
  integer, intent(in)  :: km         ! number of GEOS-5 layers
  real, intent(in)     :: ztopmax(2) ! bracking heights from PR model
  real, intent(in)     :: he(km+1)   ! GEOS-5 height above sfc
                                     ! This is atop to bottom array
  integer, intent(out) :: k1, k2     ! GEOS-5 layer indices
!
! Given vertical range ztopmax (im meters), find the
! corresponding layer index range in GEOS-5.
!

  integer :: k, k_bot, k_top

# include "PlumeRise.h"   

! Bottom
! ------
  k1 = -1
  k2 = -1
  do k = 1, km
     k_top = k
     k_bot = k+1
     if ( (ztopmax(1)>=he(k_bot)) .AND. (ztopmax(1)<he(k_top)) ) then
           k1 = k
           exit
     end if
  end do

! Top
! ---
  do k = 1, km
     k_top = k
     k_bot = k+1
     if ( (ztopmax(2)>=he(k_bot)) .AND. (ztopmax(2)<he(k_top)) ) then
           k2 = k
           exit
     end if
  end do

! Make sure k1<=k2
! ----------------
  if ( k1>k2 ) then
     k = k1
     k1 = k2
     k2 = k
  end if

! To be sure
! ----------
  k1 = max(k1,1)
  k2 = min(k2,km)

end subroutine set_geos5_vert

function vmd(zm,nkp,z_c,delta) result(v)

!
! Computes normalized vertical mass distribution given z_d and delta.
!
   implicit NONE
   integer, intent(in) :: nkp
   real, intent(in)  :: zm(nkp)
   real, intent(in)  :: z_c   ! level of maximum detrainment
   real, intent(in)  :: delta ! width of vertical mass distribution
   real :: v(nkp)
!                   ---

!  Analytic function
!  -----------------
   where ( abs(zm-z_c) <= delta )
         v = 3.*(delta**2 - (zm-z_c)**2)/(4.*delta**3) 
   elsewhere
         v = 0.0
   end where

!  Ensure normalization on grid
!  ----------------------------
   v = v / sum(v)

 end function vmd

!-------------------------------------------------------------------------

subroutine get_VMD(z_i,z_d,z_a,z_f,w,zm,nkp)

   implicit NONE
   integer, intent(in) :: nkp
   real, intent(in)  :: w(nkp) ! vertical velocity
   real, intent(in)  :: zm(nkp) ! vertical coordinate
   real, intent(out) :: z_i    ! level of maximum W
   real, intent(out) :: z_d    ! level of maximum detrainment
   real, intent(out) :: z_a    ! average height in (z_i,z_f), weighted by -dw/dz
   real, intent(out) :: z_f    ! level when w<1.

!
!  Find bottom, mid and top plume height. If using the parabolic VMD function above,
!  there are several options:
!
!  a) Like Saulo:
!                     z_m   = (z_f+z_i)/2
!                     delta = (z_f-z_i)/2
!
!  b) Preserve bottom half:
!                     z_m   = z_d
!                     delta = z_d - z_i
!
!  c) Preserve upper half:
!                     z_m   = z_d
!                     delta = z_f - z_d
!
!                       ---

   integer :: k, k_i, k_f, k_d
   !integer :: k_thresh=10
   integer :: k_thresh=10
   real    :: w_max, div_max, div, dz, z_h, w_thresh=0.2
   real*8  :: acc, wgt

   dz = zm(2) - zm(1)

!  Find level of maximum vertical velocity
!  ---------------------------------------
   w_max = 0.
   k_i = 0
   do k = nkp-10, k_thresh, -1
      if ( w(k) .gt. w_max ) then
         w_max = w(k)
         k_i = k
      end if
   end do

!  Could not find w_max, stop here
!  -------------------------------
   !if ( (k_i==0) .or. (k_i==k_thresh) ) then
   if ( (k_i==0) ) then
        z_i = -1.
        z_d = -1.
        z_a = -1.
        z_f = -1.
#if 0
        do k = nkp, 1, -1
           print *, zm(k), w(k)
        end do
#endif
        return
   else
        z_i = zm(k_i)
   end if

!  Find top plume level
!  --------------------
   k_f = 0
   do k = k_i,nkp-10
      if ( w(k) < w_thresh ) then
           k_f = k
           exit
      end if
   end do

   if ( k_f==0 ) then
        z_a = -2.
        z_d = -2.
        z_f = -2.
        return
   else
        if ( w(k_f) >= 0.0 ) then ! make sure top w is upward
           z_f = zm(k_f)
        else
           z_f = zm(k_f+1)
        end if
   end if

!  Find level of maximum detraining above w_max
!  --------------------------------------------
   div_max = 0.0
   k_d = 0
   do k = k_i, k_f
      div = -(w(k+1)-w(k))/dz
      if ( div .gt. div_max ) then
         div_max = div
         k_d = k
      end if
   end do

   if ( k_d==0 ) then
        z_a = -3.
        z_d = -3.
        return
   else
        z_d = zm(k_d)
   end if


!  Find div-weighted average height
!  --------------------------------------------
   wgt = 0.0
   acc = 0.0
   do k = k_i, k_f 
      div = -(w(k+1)-w(k))/dz
      z_h = (zm(k+1)+zm(k))/2.
      acc = acc + z_h * div
      wgt = wgt +       div
   end do
   z_a = acc / wgt

 end subroutine get_VMD

!-------------------------------------------------------------------------

 SUBROUTINE set_VMD_as_Saulo(this,nkp,VMD)

    !- version 2
    INTEGER,  INTENT(IN)  :: nkp                ! number of (edge) levels
    REAL    , INTENT(OUT), OPTIONAL :: VMD(nkp) ! vertical mass distribution function

    real   w_thresold,xxx
    integer k_i,k_f,ko,kk4,kl

# include "PlumeRise.h"   
    
    !- version 2    
    !- vertical mass distribution
    w_thresold = 1.
    
    vmd(1:nkp)= 0.
    k_i = 0
    k_f   = 0
    
    !- define range of the upper detrainemnt layer
    do ko=nkp-10,2,-1
        if(w(ko) < w_thresold) cycle
        if(k_f==0) k_f=ko
        if(w(ko)-1. > w(ko-1)) then
          k_i=ko
          exit
        endif
     enddo

     !- if there is a non zero depth layer, make the mass vertical distribution 
     if(k_f > 0 .and. k_i > 0) then 
       
        k_i=int((k_f+k_i)*0.5)
       
        !- parabolic vertical distribution between k_i and k_f
        kk4 = k_f-k_i+2
        do ko=1,kk4-1
           kl=ko+k_i-1
           VMD(kl) = 6.* float(ko)/float(kk4)**2 * (1. - float(ko)/float(kk4))
        enddo
        if(sum(VMD(1:NKP)) .ne. 1.) then
           xxx = ( 1.- sum(VMD(1:NKP)) )/float(k_f-k_i+1)
           do ko=k_i,k_f
              VMD(ko) = VMD(ko) + xxx !- values between 0 and 1.
           end do
        endif

     end if !k_f > 0 .and. k_i > 

  END SUBROUTINE set_VMD_as_Saulo


!-------------------------------------------------------------------------

subroutine get_fire_properties(this,imm,iveg_ag,fire_area)
implicit none
integer ::  moist,  i,  icount,imm,iveg_ag
real::   bfract,  effload,  heat,  hinc ,fire_area,heat_fluxW
real,    dimension(2,N_BIOME) :: heat_flux


data heat_flux/  &
!---------------------------------------------------------------------
!  heat flux      !IGBP Land Cover	    ! 
! min  ! max      !Legend and		    ! reference
!    kW/m^2       !description  	    ! 
!--------------------------------------------------------------------
 30.0,	 80.0,   &! Tropical Forest         ! igbp 2 & 4
 30.0,   80.0,   &! Boreal forest           ! igbp 1 & 3
  4.4,	 23.0,   &! cerrado/woody savanna   | igbp  5 thru 9
  3.3,	  3.3    /! Grassland/cropland      ! igbp 10 thru 17
!--------------------------------------------------------------------

#include "PlumeRise.h"  


!-- fire at the surface
!
!area = 20.e+4   ! area of burn, m^2
area = fire_area ! area of burn, m^2

!fluxo de calor para o bioma
heat_fluxW = heat_flux(imm,iveg_ag) * 1000. ! converte para W/m^2

mdur = 53        ! duration of burn, minutes
bload = 10.      ! total loading, kg/m**2 
!!!moist = 10       ! fuel moisture, %. average fuel moisture,percent dry
moist = 40       !ams per Val Martin (2013), Trentmann et al. (2006)
!srf-24jan2007
!maxtime =mdur+2  ! model time, min
 maxtime =mdur-1  ! model time, min
!srf-24jan2007 - end
!heat = 21.e6    !- joules per kg of fuel consumed                   
!heat = 15.5e6   !joules/kg - cerrado
heat = 19.3e6    !joules/kg - floresta em alta floresta (mt)
!srf - 04may2009
!alpha = 0.1      !- entrainment constant
alpha = 0.05      !- entrainment constant
!srf-end
!-------------------- printout ----------------------------------------

!!WRITE ( * ,  * ) ' SURFACE =', ZSURF, 'M', '  LCL =', ZBASE, 'M'  
!
!PRINT*,'======================================================='
!print * , ' FIRE BOUNDARY CONDITION   :'  
!print * , ' DURATION OF BURN, MINUTES =',MDUR  
!print * , ' AREA OF BURN, HA	      =',AREA*1.e-4
!print * , ' HEAT FLUX, kW/m^2	      =',heat_fluxW*1.e-3
!print * , ' TOTAL LOADING, KG/M**2    =',BLOAD  
!print * , ' FUEL MOISTURE, %	      =',MOIST !average fuel moisture,percent dry
!print * , ' MODEL TIME, MIN.	      =',MAXTIME  
!
!
!
! ******************** fix up inputs *********************************
!
                                             
!IF (MOD (MAXTIME, 2) .NE.0) MAXTIME = MAXTIME+1  !make maxtime even
                                                  
MAXTIME = MAXTIME * 60  ! and put in seconds
!
RSURF = SQRT (AREA / 3.14159) !- entrainment surface radius (m)

FMOIST   = MOIST / 100.       !- fuel moisture fraction
!
!
! calculate the energy flux and water content at lboundary.
! fills heating() on a minute basis. could ask for a file at this po
! in the program. whatever is input has to be adjusted to a one
! minute timescale.
!
                        
  DO I = 1, ntime         !- make sure of energy release
    HEATING (I) = 0.0001  !- avoid possible divide by 0
  enddo  
!                                  
  TDUR = MDUR * 60.       !- number of seconds in the burn

  bfract = 1.             !- combustion factor

  EFFLOAD = BLOAD * BFRACT  !- patchy burning
  
!     spread the burning evenly over the interval
!     except for the first few minutes for stability
  ICOUNT = 1  
!
  if(MDUR > NTIME) STOP 'Increase time duration (ntime) in min - see file "plumerise_mod.f90"'

  DO WHILE (ICOUNT.LE.MDUR)                             
!  HEATING (ICOUNT) = HEAT * EFFLOAD / TDUR  ! W/m**2 
!  HEATING (ICOUNT) = 80000.  * 0.55         ! W/m**2 

   HEATING (ICOUNT) = heat_fluxW  * 0.55     ! W/m**2 (0.55 converte para energia convectiva)
   ICOUNT = ICOUNT + 1  
  ENDDO  
!srf-24jan2007
!     ramp for 5 minutes
!  HINC = HEATING (1) / 4.  
!  HEATING (1) = 0.1  
!  HEATING (2) = HINC  
!  HEATING (3) = 2. * HINC  
!  HEATING (4) = 3. * HINC  
!
 if(imm==1) then
  HINC = HEATING (1) / 4.  
  HEATING (1) = 0.1  
  HEATING (2) = HINC  
  HEATING (3) = 2. * HINC  
  HEATING (4) = 3. * HINC 
 else 
  HINC = (HEATING (1) - heat_flux(imm-1,iveg_ag) * 1000. *0.55)/ 4.
  HEATING (1) = heat_flux(imm-1,iveg_ag) * 1000. *0.55 + 0.1  
  HEATING (2) = HEATING (1)+ HINC  
  HEATING (3) = HEATING (2)+ HINC  
  HEATING (4) = HEATING (3)+ HINC 
 endif
!srf-24jan2007 - end

return
end subroutine get_fire_properties
!-------------------------------------------------------------------------------
!

! AMS: For the next, the heat fux is given on input, not hardwired.

subroutine get_fire_properties_HF(this,hflux_kW,fire_area)
implicit none
real, intent(in) :: hflux_kW,fire_area
integer ::  moist,  i,  icount,imm,iveg_ag
real::   bfract,  effload,  heat,  hinc ,heat_fluxW

#include "PlumeRise.h"  


!-- fire at the surface
area = fire_area ! area of burn, m^2

!fluxo de calor 
heat_fluxW = hflux_kW * 1000. ! converte para W/m^2

mdur = 53        ! duration of burn, minutes
bload = 10.      ! total loading, kg/m**2 
!!!moist = 10    ! fuel moisture, %. average fuel moisture,percent dry
moist = 40       !ams per Val Martin (2013), Trentmann et al. (2006)
maxtime = mdur-1 ! model time, min
heat = 19.3e6    !joules/kg - floresta em alta floresta (mt)
alpha = 0.05      !- entrainment constant
                                                  
MAXTIME = MAXTIME * 60  ! and put in seconds

RSURF = SQRT (AREA / 3.14159) !- entrainment surface radius (m)

FMOIST   = MOIST / 100.       !- fuel moisture fraction

!
!
! calculate the energy flux and water content at lboundary.
! fills heating() on a minute basis. could ask for a file at this po
! in the program. whatever is input has to be adjusted to a one
! minute timescale.
!
                        
  DO I = 1, ntime         !- make sure of energy release
    HEATING (I) = 0.0001  !- avoid possible divide by 0
  enddo  
!                                  
  TDUR = MDUR * 60.       !- number of seconds in the burn

  bfract = 1.             !- combustion factor

  EFFLOAD = BLOAD * BFRACT  !- patchy burning
  
! spread the burning evenly over the interval
! except for the first few minutes for stability
  ICOUNT = 1  
!
  if(MDUR > NTIME) STOP 'Increase time duration (ntime) in min - see file "plumerise_mod.f90"'

  DO WHILE (ICOUNT.LE.MDUR)                             
   HEATING (ICOUNT) = heat_fluxW  * 0.55     ! W/m**2 (0.55 converte para energia convectiva)
   ICOUNT = ICOUNT + 1  
  ENDDO  

! ramp for 5 minutes

  HINC = HEATING (1) / 4.  
  HEATING (1) = 0.1  
  HEATING (2) = HINC  
  HEATING (3) = 2. * HINC  
  HEATING (4) = 3. * HINC 

return
end subroutine get_fire_properties_HF
!-------------------------------------------------------------------------------
!

!srf-24jan2007
!SUBROUTINE MAKEPLUME (this,kmt,ztopmax,ixx,kkmax)  
SUBROUTINE MAKEPLUME (this,kmt,ztopmax,ixx,kkmax,imm)  
!srf-24jan2007 - end
!
! *********************************************************************
!
!    EQUATION SOURCE--Kessler Met.Monograph No. 32 V.10 (K)
!    Alan Weinstein, JAS V.27 pp 246-255. (W),
!    Ogura and Takahashi, Monthly Weather Review V.99,pp895-911 (OT)
!    Roger Pielke,Mesoscale Meteorological Modeling,Academic Press,1984
!
!
! ************************ VARIABLE ID ********************************
!
!     DT=COMPUTING TIME INCREMENT (SEC)
!     DZ=VERTICAL INCREMENT (M)
!     LBASE=LEVEL ,CLOUD BASE
!
!     CONSTANTS:
!       G = GRAVITATIONAL ACCELERATION 9.80796 (M/SEC/SEC).
!       R = DRY AIR GAS CONSTANT (287.04E6 JOULE/KG/DEG K)
!       CP = SPECIFIC HT. (1004 JOULE/KG/DEG K)
!       HEATCOND = HEAT OF CONDENSATION (2.5E6 JOULE/KG)
!       HEATFUS = HEAT OF FUSION (3.336E5 JOULE/KG)
!       HEATSUBL = HEAT OF SUBLIMATION (2.83396E6 JOULE/KG)
!       EPS = RATIO OF MOL.WT. OF WATER VAPOR TO THAT OF DRY AIR (0.622)
!       DES = DIFFERENCE BETWEEN VAPOR PRESSURE OVER WATER AND ICE (MB)
!       TFREEZE = FREEZING TEMPERATURE (K)
!
!
!     PARCEL VALUES:
!       T = TEMPERATURE (K)
!       TXS = TEMPERATURE EXCESS (K)
!       QH = HYDROMETEOR WATER CONTENT (G/G DRY AIR)
!       QHI = HYDROMETEOR ICE CONTENT (G/G DRY AIR)
!       QC = WATER CONTENT (G/G DRY AIR)
!       QVAP = WATER VAPOR MIXING RATIO (G/G DRY AIR)
!       QSAT = SATURATION MIXING RATIO (G/G DRY AIR)
!       RHO = DRY AIR DENSITY (G/M**3) MASSES = RHO*Q'S IN G/M**3
!       ES = SATURATION VAPOR PRESSURE (kPa)
!
!     ENVIRONMENT VALUES:
!       TE = TEMPERATURE (K)
!       PE = PRESSURE (kPa)
!       QVENV = WATER VAPOR (G/G)
!       RHE = RELATIVE HUMIDITY FRACTION (e/esat)
!       DNE = dry air density (kg/m^3)
!
!     HEAT VALUES:
!       HEATING = HEAT OUTPUT OF FIRE (WATTS/M**2)
!       MDUR = DURATION OF BURN, MINUTES
!
!       W = VERTICAL VELOCITY (M/S)
!       RADIUS=ENTRAINMENT RADIUS (FCN OF Z)
!	RSURF = ENTRAINMENT RADIUS AT GROUND (SIMPLE PLUME, TURNER)
!	ALPHA = ENTRAINMENT CONSTANT
!       MAXTIME = TERMINATION TIME (MIN)
!
!
!**********************************************************************
!**********************************************************************               
implicit none 
!logical :: endspace  
character (len=10) :: varn
integer ::  izprint, iconv,  itime, k, kk, kkmax, deltak,ilastprint,kmt &
           ,ixx,nrectotal,i_micro,n_sub_step
real ::  vc, g,  r,  cp,  eps,  &
         tmelt,  heatsubl,  heatfus,  heatcond, tfreeze, &
         ztopmax, wmax, rmaxtime, es, esat, heat,dt_save
character (len=2) :: cixx	 
!srf-24jan2007
integer :: imm

real :: z_i, z_d, z_a, z_f

!srf-24jan2007-end
!
!
! ******************* SOME CONSTANTS **********************************
!
!      XNO=10.0E06 median volume diameter raindrop (K table 4)
!      VC = 38.3/(XNO**.125) mean volume fallspeed eqn. (K)
!
parameter (vc = 5.107387)  
parameter (g = 9.80796, r = 287.04, cp = 1004., eps = 0.622,  tmelt = 273.3)
parameter (heatsubl = 2.834e6, heatfus = 3.34e5, heatcond = 2.501e6)
parameter (tfreeze = 269.3)  
#include "PlumeRise.h" 
!
!
tstpf = 2.0  	!- timestep factor
viscosity = 500.!- viscosity constant (original value: 0.001)

nrectotal=150
!
!*************** PROBLEM SETUP AND INITIAL CONDITIONS *****************
mintime = 1  
ztopmax = 0. 
ztop    = 0. 
   time = 0.  
     dt = 1.
   wmax = 1. 
kkmax   = 10
deltaK  = 20
ilastprint=0
L       = 1   ! L initialization

!--- initialization
!srf-24jan2007 
! CALL INITIAL(this,kmt)  
if(imm==1) CALL INITIAL(this,kmt)  
!srf-24jan2007 - end

!--- initial print fields:
izprint  = 0          ! if = 0 => no printout
if (izprint.ne.0) then
 write(cixx(1:2),'(i2.2)') ixx
 open(2, file = 'debug.'//cixx//'.dat')  
 open(19,file='plumegen9.'//cixx//'.gra',         &
     form='unformatted',access='direct',status='unknown',  &
     recl=8*nrectotal)  !PC   
!     recl=1*nrectotal) !sx6 e tupay
 call printout (this,izprint,nrectotal)
 ilastprint=2
endif     

! ******************* model evolution ******************************
rmaxtime = float(maxtime)
!
 DO WHILE (TIME.LE.RMAXTIME)  !beginning of time loop

   !do itime=1,120

!-- set model top integration
    nm1 = min(kmt, kkmax + deltak)
                                    
!-- set timestep
    !dt = (zm(2)-zm(1)) / (tstpf * wmax)  
    !dt = min(5.,(zm(2)-zm(1)) / (tstpf * wmax))
    dt = min(1.,(zm(2)-zm(1)) / (tstpf * wmax))
                                
!-- elapsed time, sec
    time = time+dt 
!-- elapsed time, minutes                                      
    mintime = 1 + int (time) / 60     
    wmax = 1.  !no zeroes allowed.
!************************** BEGIN SPACE LOOP **************************

!-- zerout all model tendencies
    call tend0_plumerise(this)

!-- bounday conditions (k=1)
    L=1
    call lbound(this)

!-- dynamics for the level k>1 
!-- W advection 
!   call vel_advectc_plumerise(NM1,WC,WT,DNE,DZM)
    call vel_advectc_plumerise(NM1,WC,WT,RHO,DZM)
  
!-- scalars advection 1
    call scl_advectc_plumerise(this,'SC',NM1)

!-- scalars entrainment, adiabatic
    call scl_misc(this,NM1)
!srf-04may2009-AWE
!-- scalars dinamic entrainment
    call  scl_dyn_entrain(this,NM1)
!-srf-end


!-- gravity wave damping using Rayleigh friction layer fot T
    call damp_grav_wave(1,nm1,deltak,dt,zt,zm,w,t,tt,qv,qh,qi,qc,te,pe,qvenv&
                       ,vel_p,vel_t,vel_e)!srf-04may2009-AWE

!-- microphysics
!   goto 101 ! bypass microphysics
    dt_save=dt
    n_sub_step=3
    dt=dt/float(n_sub_step)

    do i_micro=1,n_sub_step
!-- sedim ?
     call fallpart(this,NM1)
!-- microphysics
     do L=2,nm1-1
        WBAR    = 0.5*(W(L)+W(L-1))
        !ES      = 0.1*ESAT (T(L))   !BLOB SATURATION VAPOR PRESSURE, EM KPA
                                   ! rotina do plumegen calcula em kPa
        ES       = ESAT_PR (T(L))  ! BLOB SATURATION VAPOR PRESSURE, EM KPA
        QSAT(L) = (EPS * ES) / (PE(L) - ES)  !BLOB SATURATION LWC G/G DRY AIR
        EST (L) = ES  
        RHO (L) = 3483.8 * PE (L) / T (L) ! AIR PARCEL DENSITY , G/M**3
	IF (W(L) .ge. 0.) then 
	   DQSDZ = (QSAT(L+1) - QSAT(L-1)) / (ZT(L+1 )-ZT(L-1))
	ELSE
	   DQSDZ = (QSAT(L+1) - QSAT(L-1)) / (ZT(L+1) -ZT(L-1))
	ENDIF 
	
	call waterbal(this)
     enddo
    enddo
    dt=dt_save
!
    101 continue
!
!-- W-viscosity for stability 
    call visc_W(this,nm1,deltak,kmt)

!-- update scalars
    call update_plumerise(this,nm1,'S')
    
    call hadvance_plumerise(1,nm1,dt,WC,WT,W,mintime) 

!-- Buoyancy
    call buoyancy_plumerise(NM1, T, TE, QV, QVENV, QH, QI, QC, WT, SCR1)
 
!-- Entrainment 
    call entrainment(NM1,W,WT,RADIUS,ALPHA,vel_p,vel_e)

!-- update W
    call update_plumerise(this,nm1,'W')

    call hadvance_plumerise(2,nm1,dt,WC,WT,W,mintime) 


!-- misc
    do k=2,nm1
!    pe esta em kpa  - esat do rams esta em mbar = 100 Pa = 0.1 kpa
!     es       = 0.1*esat (t(k)) !blob saturation vapor pressure, em kPa
!    rotina do plumegen calcula em kPa
     es       = esat_pr (t(k))  !blob saturation vapor pressure, em kPa
     qsat(k) = (eps * es) / (pe(k) - es)  !blob saturation lwc g/g dry air
     est (k) = es  
     txs (k) = t(k) - te(k)
     rho (k) = 3483.8 * pe (k) / t (k) ! air parcel density , g/m**3
                                       ! no pressure diff with radius
				       
     if((abs(wc(k))).gt.wmax) wmax = abs(wc(k)) ! keep wmax largest w
    enddo  

! Gravity wave damping using Rayleigh friction layer for W
    call damp_grav_wave(2,nm1,deltak,dt,zt,zm,w,t,tt,qv,qh,qi,qc,te,pe,qvenv&
                       ,vel_p,vel_t,vel_e)!srf-04may2009-AWE
!---
!-srf-04may2009-AWE
!update radius 
    do k=2,nm1
     radius(k) = rad_p(k)
    enddo
!-srf-end

!-- try to find the plume top (from the top)
    kk = 1
    do while (w (kk) .gt. 1.)  
     kk = kk + 1  
     ztop =  zm(kk) 
!     print*,'W=',w (kk)
    enddo  

!
    ztop_(mintime) = ztop
    ztopmax = max (ztop, ztopmax) 
    kkmax   = max (kk  , kkmax  ) 

!
!srf-27082005
! if the solution is going to a stationary phase, exit
!srf-24jan2007 - end
!   if(mintime > 10) then
   if(mintime > 20) then
!srf-24jan2007 - end
   
    if( abs(ztop_(mintime)-ztop_(mintime-10)) < DZ ) then
        if ( isNAN(w(nm1/2)) ) ztopmax = w(nm1/2) ! NaN
        exit
    end if

   endif

    if(ilastprint == mintime) then
      call printout (this,izprint,nrectotal)  
      ilastprint = mintime+1
    endif      
                               
END DO   !do next timestep

!ams-- redefine ztopmax as the level of maximum detrainment

    if ( .NOT. isnan(ztopmax) ) then

       call get_VMD(z_i,z_d,z_a,z_f,w,zm,nkp)

#if 0
    print *
    print *, '>>> z_i,z_d,z_a,z_f: ', int(z_i),int(z_d),nint(z_a),int(z_f)
    do kk = nm1, 1, -1
       if ( z_f .eq. zm(kk) ) then
           print *,'f W =', zm(kk), w(kk), t(kk)-273.25
       else if ( z_i .eq. zm(kk) ) then
           print *,'i W =', zm(kk), w(kk), t(kk)-273.25
       else if ( z_d .eq. zm(kk) ) then
           print *,'d W =', zm(kk), w(kk), t(kk)-273.25
       else if ( ztopmax .eq. zm(kk) ) then
           print *,'+ W =', zm(kk), w(kk), t(kk)-273.25
       else
           print *,'  W =', zm(kk), w(kk), t(kk)-273.25
       end if
    end do

    print * ,'ztopmax=', mintime,'mn ',ztop_(mintime), ztopmax, z_a
#endif

    ztopmax = z_a

   end if

!ams

!print * ,' ztopmax=',ztopmax,'m',mintime,'mn '
!print*,'======================================================='
!

!the last printout
if (izprint.ne.0) then
 call printout (this,izprint,nrectotal)  
 close (2)            
 close (19)            
endif

RETURN  
END SUBROUTINE MAKEPLUME
!-------------------------------------------------------------------------------
!
SUBROUTINE BURN(this,EFLUX, WATER)  
!	
!- calculates the energy flux and water content at lboundary
!real, parameter :: HEAT = 21.E6 !Joules/kg
!real, parameter :: HEAT = 15.5E6 !Joules/kg - cerrado
real, parameter :: HEAT = 19.3E6 !Joules/kg - floresta em Alta Floresta (MT)
real :: eflux,water
#include "PlumeRise.h"                               
!
! The emission factor for water is 0.5. The water produced, in kg,
! is then  fuel mass*0.5 + (moist/100)*mass per square meter.
! The fire burns for DT out of TDUR seconds, the total amount of
! fuel burned is AREA*BLOAD*(DT/TDUR) kg. this amount of fuel is
! considered to be spread over area AREA and so the mass burned per
! unit area is BLOAD*(DT/TDUR), and the rate is BLOAD/TDUR.
!        
IF (TIME.GT.TDUR) THEN !is the burn over?   
   EFLUX = 0.000001    !prevent a potential divide by zero
   WATER = 0.  
   RETURN  
ELSE  
!                                                   
   EFLUX = HEATING (MINTIME)                          ! Watts/m**2                                                   
!  WATER = EFLUX * (DT / HEAT) * (0.5 + FMOIST)       ! kg/m**2 
   WATER = EFLUX * (DT / HEAT) * (0.5 + FMOIST) /0.55 ! kg/m**2 
   WATER = WATER * 1000.                              ! g/m**2
!
!        print*,'BURN:',time,EFLUX/1.e+9
ENDIF  
!
RETURN  
END SUBROUTINE BURN
!-------------------------------------------------------------------------------
!
SUBROUTINE LBOUND (this)  
!
! ********** BOUNDARY CONDITIONS AT ZSURF FOR PLUME AND CLOUD ********
!
! source of equations: J.S. Turner Buoyancy Effects in Fluids
!                      Cambridge U.P. 1973 p.172,
!                      G.A. Briggs Plume Rise, USAtomic Energy Commissio
!                      TID-25075, 1969, P.28
!
! fundamentally a point source below ground. at surface, this produces
! a velocity w(1) and temperature T(1) which vary with time. There is
! also a water load which will first saturate, then remainder go into
! QC(1).
! EFLUX = energy flux at ground,watt/m**2 for the last DT
!
implicit none
real, parameter :: g = 9.80796, r = 287.04, cp = 1004.6, eps = 0.622,tmelt = 273.3
real, parameter :: tfreeze = 269.3, pi = 3.14159, e1 = 1./3., e2 = 5./3.
real :: es,  esat, eflux, water,  pres, c1,  c2, f, zv,  denscor, xwater ! ,ESAT_PR
!            
#include "PlumeRise.h"  
QH (1) = QH (2)   !soak up hydrometeors
QI (1) = QI (2)              
QC (1) = 0.       !no cloud here
!
!
   CALL BURN (this,EFLUX, WATER)  
!
!  calculate parameters at boundary from a virtual buoyancy point source
!
   PRES = PE (1) * 1000.   !need pressure in N/m**2
                              
   C1 = 5. / (6. * ALPHA)  !alpha is entrainment constant

   C2 = 0.9 * ALPHA  

   F = EFLUX / (PRES * CP * PI)  
                             
   F = G * R * F * AREA  !buoyancy flux
                 
   ZV = C1 * RSURF  !virtual boundary height
                                   
   W (1) = C1 * ( (C2 * F) **E1) / ZV**E1  !boundary velocity
                                         
   DENSCOR = C1 * F / G / (C2 * F) **E1 / ZV**E2   !density correction

   T (1) = TE (1) / (1. - DENSCOR)    !temperature of virtual plume at zsurf
   
!
   WC(1) = W(1)

!srf-04may2009-AWE
   VEL_P(1) = 0.
   rad_p(1) = rsurf
!-srf-end

   !SC(1) = SCE(1)+F/1000.*dt  ! gas/particle (g/g)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     match dw/dz,dt/dz at the boundary. F is conserved.
!
   !WBAR = W (1) * (1. - 1. / (6. * ZV) )  
   !ADVW = WBAR * W (1) / (3. * ZV)  
   !ADVT = WBAR * (5. / (3. * ZV) ) * (DENSCOR / (1. - DENSCOR) )  
   !ADVC = 0.  
   !ADVH = 0.  
   !ADVI = 0.  
   !ADIABAT = - WBAR * G / CP  
   VTH (1) = - 4.  
   VTI (1) = - 3.  
   TXS (1) = T (1) - TE (1)  

   VISC (1) = VISCOSITY  

   RHO (1) = 3483.8 * PE (1) / T (1)   !air density at level 1, g/m**3

   XWATER = WATER / (W (1) * DT * RHO (1) )   !firewater mixing ratio
                                            
   QV (1) = XWATER + QVENV (1)  !plus what's already there 


!  PE esta em kPa  - ESAT do RAMS esta em mbar = 100 Pa = 0.1 kPa
!   ES       = 0.1*ESAT (T(1)) !blob saturation vapor pressure, em kPa
!  rotina do plumegen ja calcula em kPa
   ES       = ESAT_PR (T(1))  !blob saturation vapor pressure, em kPa

   EST  (1)  = ES                                  
   QSAT (1) = (EPS * ES) / (PE (1) - ES)   !blob saturation lwc g/g dry air
  
   IF (QV (1) .gt. QSAT (1) ) THEN  
       QC (1) = QV   (1) - QSAT (1) + QC (1)  !remainder goes into cloud drops
       QV (1) = QSAT (1)  
   ENDIF  
!
   CALL WATERBAL(this)  
!
RETURN  
END SUBROUTINE LBOUND
!-------------------------------------------------------------------------------
!
SUBROUTINE INITIAL (this, kmt)  
!
! ************* SETS UP INITIAL CONDITIONS FOR THE PROBLEM ************
implicit none 
real, parameter :: tfreeze = 269.3
integer ::  isub,  k,  n1,  n2,  n3,  lbuoy,  itmp,  isubm1 ,kmt
real ::     xn1,  xi,  es,  esat !,ESAT_PR
#include "PlumeRise.h" 
!
N=kmt
! initialize temperature structure,to the end of equal spaced sounding,
  do k = 1, N			  
  TXS (k) = 0.0  
    W (k) = 0.0             
    T (k) = TE(k)       !blob set to environment		  
    WC(k) = 0.0
    WT(k) = 0.0
    QV(k) = QVENV (k)   !blob set to environment             
   VTH(k) = 0.		!initial rain velocity = 0	                     
   VTI(k) = 0.		!initial ice  velocity = 0	                     
    QH(k) = 0.		!no rain			     
    QI(k) = 0.		!no ice 			     
    QC(k) = 0.		!no cloud drops	                     
!  PE esta em kPa  - ESAT do RAMS esta em mbar = 100 Pa = 0.1 kPa
!   ES       = 0.1*ESAT (T(k)) !blob saturation vapor pressure, em kPa
!  rotina do plumegen calcula em kPa
   ES       = ESAT_PR (T(k))  !blob saturation vapor pressure, em kPa
   EST  (k) = ES  
   QSAT (k) = (.622 * ES) / (PE (k) - ES) !saturation lwc g/g
   RHO  (k) = 3483.8 * PE (k) / T (k) 	!dry air density g/m**3    
!srf-04may2009
   VEL_P(k) = 0.
   rad_p(k) = 0.
  enddo  

! Initialize the entrainment radius, Turner-style plume
  radius(1) = rsurf
  do k=2,N
     radius(k) = radius(k-1)+(6./5.)*alpha*(zt(k)-zt(k-1))
     !-srf-04may2009
     rad_p(k)  = radius(k)
  enddo
  !-srf-04may2009
  rad_p(1) = rsurf
    
!  Initialize the viscosity
   VISC (1) = VISCOSITY
   do k=2,N
     VISC (k) = VISCOSITY!max(1.e-3,visc(k-1) - 1.* VISCOSITY/float(nkp))
   enddo
!--   Initialize gas/concentration
  !DO k =10,20
  !   SC(k) = 20.
  !ENDDO
  !stop 333

   CALL LBOUND(this)

RETURN  
END SUBROUTINE INITIAL
!-------------------------------------------------------------------------------
!
subroutine damp_grav_wave(ifrom,nm1,deltak,dt,zt,zm,w,t,tt,qv,qh,qi,qc,te,pe,qvenv&
                         ,vel_p,vel_t,vel_e)
implicit none
integer nm1,ifrom,deltak
real dt
real, dimension(nm1) :: w,t,tt,qv,qh,qi,qc,te,pe,qvenv,dummy,zt,zm&
                       ,vel_p,vel_t,vel_e

if(ifrom==1) then
 call friction(ifrom,nm1,deltak,dt,zt,zm,t,tt    ,te)
! call friction(ifrom,nm1,deltak,dt,zt,zm,vel_p,vel_t,vel_e)
!call friction(ifrom,nm1,dt,zt,zm,qv,qvt,qvenv)
 return
endif 

dummy(:) = 0.
if(ifrom==2) call friction(ifrom,nm1,deltak,dt,zt,zm,w,dummy ,dummy)
!call friction(ifrom,nm1,dt,zt,zm,qi,qit ,dummy)
!call friction(ifrom,nm1,dt,zt,zm,qh,qht ,dummy)
!call friction(ifrom,nm1,dt,zt,zm,qc,qct ,dummy)
return
end subroutine damp_grav_wave
!-------------------------------------------------------------------------------
!
subroutine friction(ifrom,nm1,deltak,dt,zt,zm,var1,vart,var2)
implicit none
real, dimension(nm1) :: var1,var2,vart,zt,zm
integer k,nfpt,kf,nm1,ifrom,deltak
real zmkf,ztop,distim,c1,c2,dt

!nfpt=50
!kf = nm1 - nfpt
kf = nm1 - int(deltak)

zmkf = zm(kf) !old: float(kf )*dz
ztop = zm(nm1)
!distim = min(4.*dt,200.)
!distim = 60.
distim = min(3.*dt,60.)

c1 = 1. / (distim * (ztop - zmkf))
c2 = dt * c1

if(ifrom == 1) then  
  do k = nm1,2,-1
   if (zt(k) .le. zmkf) cycle
   vart(k) = vart(k)   + c1 * (zt(k) - zmkf)*(var2(k) - var1(k))
  enddo
elseif(ifrom == 2) then
  do k = nm1,2,-1
   if (zt(k) .le. zmkf) cycle
   var1(k) =  var1(k) + c2 * (zt(k) - zmkf)*(var2(k) - var1(k))
  enddo
endif
return
end subroutine friction
!-------------------------------------------------------------------------------
!
subroutine vel_advectc_plumerise(m1,wc,wt,rho,dzm)

implicit none
integer :: k,m1
real, dimension(m1) :: wc,wt,flxw,dzm,rho
real, dimension(m1) :: dn0 ! var local
real :: c1z

!dzm(:)= 1./dz

dn0(1:m1)=rho(1:m1)*1.e-3 ! converte de cgs para mks

flxw(1) = wc(1) * dn0(1) 

do k = 2,m1-1
   flxw(k) = wc(k) * .5 * (dn0(k) + dn0(k+1))
enddo

! Compute advection contribution to W tendency

c1z = .5 

do k = 2,m1-2

   wt(k) = wt(k)  &
      + c1z * dzm(k) / (dn0(k) + dn0(k+1)) *     (   &
	(flxw(k) + flxw(k-1))  * (wc(k) + wc(k-1))   &
      - (flxw(k) + flxw(k+1))  * (wc(k) + wc(k+1))   &
      + (flxw(k+1) - flxw(k-1)) * 2.* wc(k)       )

enddo

return
end subroutine vel_advectc_plumerise
!-------------------------------------------------------------------------------
!
subroutine hadvance_plumerise(iac,m1,dt,wc,wt,wp,mintime)

implicit none
integer :: k,iac
integer :: m1,mintime
real, dimension(m1) :: dummy, wc,wt,wp
real eps,dt
!     It is here that the Asselin filter is applied.  For the velocities
!     and pressure, this must be done in two stages, the first when
!     IAC=1 and the second when IAC=2.


eps = .2
if(mintime == 1) eps=0.5

!     For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P.
!
call predict_plumerise(m1,wc,wp,wt,dummy,iac,2.*dt,eps)
!print*,'mintime',mintime,eps
!do k=1,m1
!   print*,'W-HAD',k,wc(k),wp(k),wt(k)
!enddo
return
end subroutine hadvance_plumerise
!-------------------------------------------------------------------------------
!
subroutine predict_plumerise(npts,ac,ap,fa,af,iac,dtlp,epsu)
implicit none
integer :: npts,iac,m
real :: epsu,dtlp
real, dimension(*) :: ac,ap,fa,af

!     For IAC=3, this routine moves the arrays AC and AP forward by
!     1 time level by adding in the prescribed tendency. It also
!     applies the Asselin filter given by:

!              {AC} = AC + EPS * (AP - 2 * AC + AF)

!     where AP,AC,AF are the past, current and future time levels of A.
!     All IAC=1 does is to perform the {AC} calculation without the AF
!     term present.  IAC=2 completes the calculation of {AC} by adding
!     the AF term only, and advances AC by filling it with input AP
!     values which were already updated in ACOUSTC.
!

if (iac .eq. 1) then
   do m = 1,npts
      ac(m) = ac(m) + epsu * (ap(m) - 2. * ac(m))
   enddo
   return
elseif (iac .eq. 2) then
   do m = 1,npts
      af(m) = ap(m)
      ap(m) = ac(m) + epsu * af(m)
   enddo
!elseif (iac .eq. 3) then
!   do m = 1,npts
!      af(m) = ap(m) + dtlp * fa(m)
!   enddo
!   if (ngrid .eq. 1 .and. ipara .eq. 0) call cyclic(nzp,nxp,nyp,af,'T')
!   do m = 1,npts
!      ap(m) = ac(m) + epsu * (ap(m) - 2. * ac(m) + af(m))
!   enddo
endif

do m = 1,npts
  ac(m) = af(m)
enddo
return
end subroutine predict_plumerise
!-------------------------------------------------------------------------------
!
subroutine  buoyancy_plumerise(m1, T, TE, QV, QVENV, QH, QI, QC, WT, scr1)
implicit none
integer :: k,m1
real, parameter :: g = 9.8, eps = 0.622, gama = 0.5 ! mass virtual coeff.
real, dimension(m1) :: T, TE, QV, QVENV, QH, QI, QC, WT, scr1
real :: TV,TVE,QWTOTL,umgamai
real, parameter :: mu = 0.15 

!- orig
umgamai = 1./(1.+gama) ! compensa a falta do termo de aceleracao associado `as
                       ! das pertubacoes nao-hidrostaticas no campo de pressao

!- new                 ! Siesbema et al, 2004
!umgamai = 1./(1.-2.*mu)

do k = 2,m1-1

    TV =   T(k) * (1. + (QV(k)   /EPS))/(1. + QV(k)   )  !blob virtual temp.                                        	   
    TVE = TE(k) * (1. + (QVENV(k)/EPS))/(1. + QVENV(k))  !and environment

    QWTOTL = QH(k) + QI(k) + QC(k)                       ! QWTOTL*G is drag
!- orig
   !scr1(k)= G*( umgamai*(  TV - TVE) / TVE   - QWTOTL) 
    scr1(k)= G*  umgamai*( (TV - TVE) / TVE   - QWTOTL) 

    !if(k .lt. 10)print*,'BT',k,TV,TVE,TVE,QWTOTL
enddo

do k = 2,m1-2
    wt(k) = wt(k)+0.5*(scr1(k)+scr1(k+1))
!   print*,'W-BUO',k,wt(k),scr1(k),scr1(k+1)
enddo

end subroutine  buoyancy_plumerise
!-------------------------------------------------------------------------------
!
subroutine ENTRAINMENT(m1,w,wt,radius,ALPHA,vel_p,vel_e)
implicit none
integer :: k,m1
real, dimension(m1) :: w,wt,radius,vel_p,vel_e
real DMDTM,ALPHA,WBAR,RADIUS_BAR,umgamai,DYN_ENTR
real, parameter :: mu = 0.15 ,gama = 0.5 ! mass virtual coeff.

!- new - Siesbema et al, 2004
!umgamai = 1./(1.-2.*mu)

!- orig
!umgamai = 1
umgamai = 1./(1.+gama) ! compensa a falta do termo de aceleracao associado `as
                       !  pertubacoes nao-hidrostaticas no campo de pressao

!
!-- ALPHA/RADIUS(L) = (1/M)DM/DZ  (W 14a)
  do k=2,m1-1

!-- for W: WBAR is only W(k)
!     WBAR=0.5*(W(k)+W(k-1))           
      WBAR=W(k)          
      RADIUS_BAR = 0.5*(RADIUS(k) + RADIUS(k-1))
! orig
     !DMDTM =           2. * ALPHA * ABS (WBAR) / RADIUS_BAR  != (1/M)DM/DT
      DMDTM = umgamai * 2. * ALPHA * ABS (WBAR) / RADIUS_BAR  != (1/M)DM/DT

!--  DMDTM*W(L) entrainment,
      wt(k) = wt(k)  - DMDTM*ABS (WBAR)


!-   dynamic entrainment
     DYN_ENTR =  (2./3.1416)*0.5*ABS (VEL_P(k)-VEL_E(k)+VEL_P(k-1)-VEL_E(k-1)) /RADIUS_BAR

     wt(k) = wt(k)  - DYN_ENTR*ABS (WBAR)

      !print*,'W-ENTR=',k,w(k),- DMDTM*ABS (WBAR)
  enddo
end subroutine  ENTRAINMENT
!-------------------------------------------------------------------------------
!
subroutine scl_advectc_plumerise(this,varn,mzp)
implicit none
integer :: mzp
character(len=*) :: varn
real :: dtlto2
integer :: k
#include "PlumeRise.h"

!  wp => w
!- Advect  scalars
   dtlto2   = .5 * dt
!  vt3dc(1) =      (w(1) + wc(1)) * dtlto2 * dne(1)
   vt3dc(1) =      (w(1) + wc(1)) * dtlto2 * rho(1)*1.e-3!converte de CGS p/ MKS
   vt3df(1) = .5 * (w(1) + wc(1)) * dtlto2 * dzm(1)

   do k = 2,mzp
!     vt3dc(k) =  (w(k) + wc(k)) * dtlto2 *.5 * (dne(k) + dne(k+1))
      vt3dc(k) =  (w(k) + wc(k)) * dtlto2 *.5 * (rho(k) + rho(k+1))*1.e-3
      vt3df(k) =  (w(k) + wc(k)) * dtlto2 *.5 *  dzm(k)
     !print*,'vt3df-vt3dc',k,vt3dc(k),vt3df(k)
   enddo

 
!-srf-24082005
!  do k = 1,mzp-1
  do k = 1,mzp
     vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
     vctr2(k) = (zm(k)   - zt(k)) * dzm(k)
!    vt3dk(k) = dzt(k) / dne(k)
     vt3dk(k) = dzt(k) /(rho(k)*1.e-3)
     !print*,'VT3dk',k,dzt(k) , dne(k)
  enddo

!      scalarp => scalar_tab(n,ngrid)%var_p
!      scalart => scalar_tab(n,ngrid)%var_t

!- temp advection tendency (TT)
   scr1=T
   call fa_zc_plumerise(mzp                   &
             	       ,T	  ,scr1  (1:)  &
             	       ,vt3dc (1:) ,vt3df (1:)  &
             	       ,vt3dg (1:) ,vt3dk (1:)  &
             	       ,vctr1,vctr2	      )

   call advtndc_plumerise(mzp,T,scr1(1:),TT,dt)

!- water vapor advection tendency (QVT)
   scr1=QV
   call fa_zc_plumerise(mzp                  &
             	       ,QV	  ,scr1  (1:)  &
             	       ,vt3dc (1:) ,vt3df (1:)  &
             	       ,vt3dg (1:) ,vt3dk (1:)  &
             	       ,vctr1,vctr2	     )

   call advtndc_plumerise(mzp,QV,scr1(1:),QVT,dt)

!- liquid advection tendency (QCT)
   scr1=QC
   call fa_zc_plumerise(mzp                  &
             	       ,QC	  ,scr1  (1:)  &
             	       ,vt3dc (1:) ,vt3df (1:)  &
             	       ,vt3dg (1:) ,vt3dk (1:)  &
             	       ,vctr1,vctr2	     )

   call advtndc_plumerise(mzp,QC,scr1(1:),QCT,dt)

!- ice advection tendency (QIT)
   scr1=QI
   call fa_zc_plumerise(mzp                  &
             	       ,QI	  ,scr1  (1:)  &
             	       ,vt3dc (1:) ,vt3df (1:)  &
             	       ,vt3dg (1:) ,vt3dk (1:)  &
             	       ,vctr1,vctr2	     )

   call advtndc_plumerise(mzp,QI,scr1(1:),QIT,dt)

!- hail/rain advection tendency (QHT)
!   if(ak1 > 0. .or. ak2 > 0.) then

      scr1=QH
      call fa_zc_plumerise(mzp                  &
             	          ,QH	    ,scr1  (1:)  &
             	          ,vt3dc (1:) ,vt3df (1:)  &
             	          ,vt3dg (1:) ,vt3dk (1:)  &
             	          ,vctr1,vctr2	       )

      call advtndc_plumerise(mzp,QH,scr1(1:),QHT,dt)
!   endif

!-srf-04may2009-AWE
!- horizontal wind advection tendency (VEL_T)

      scr1=VEL_P
      call fa_zc_plumerise(mzp                  &
             	          ,VEL_P     ,scr1  (1:)  &
             	          ,vt3dc (1:) ,vt3df (1:)  &
             	          ,vt3dg (1:) ,vt3dk (1:)  &
             	          ,vctr1,vctr2	       )

      call advtndc_plumerise(mzp,VEL_P,scr1(1:),VEL_T,dt)

!- vertical radius transport

      scr1=rad_p
      call fa_zc_plumerise(mzp                  &
             	          ,rad_p     ,scr1  (1:)  &
             	          ,vt3dc (1:) ,vt3df (1:)  &
             	          ,vt3dg (1:) ,vt3dk (1:)  &
             	          ,vctr1,vctr2	       )

      call advtndc_plumerise(mzp,rad_p,scr1(1:),rad_t,dt)
!srf-end
   return

!- gas/particle advection tendency (SCT)
!    if(varn == 'SC')return
   scr1=SC
   call fa_zc_plumerise(mzp		    &
   	     	       ,SC	 ,scr1  (1:)  &
   	     	       ,vt3dc (1:) ,vt3df (1:)  &
   	     	       ,vt3dg (1:) ,vt3dk (1:)  &
   	     	       ,vctr1,vctr2	     )
   
   call advtndc_plumerise(mzp,SC,scr1(1:),SCT,dt)


return
end subroutine scl_advectc_plumerise
!-------------------------------------------------------------------------------
!
subroutine fa_zc_plumerise(m1,scp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2)

implicit none
integer :: m1,k
real :: dfact
real, dimension(m1) :: scp,scr1,vt3dc,vt3df,vt3dg,vt3dk
real, dimension(m1) :: vctr1,vctr2

dfact = .5

! Compute scalar flux VT3DG
      do k = 1,m1-1
         vt3dg(k) = vt3dc(k)                   &
                  * (vctr1(k) * scr1(k)        &
                  +  vctr2(k) * scr1(k+1)      &
                  +  vt3df(k) * (scr1(k) - scr1(k+1)))
      enddo
      
! Modify fluxes to retain positive-definiteness on scalar quantities.
!    If a flux will remove 1/2 quantity during a timestep,
!    reduce to first order flux. This will remain positive-definite
!    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
!    both fluxes are evacuating the box.

do k = 1,m1-1
 if (vt3dc(k) .gt. 0.) then
   if (vt3dg(k) * vt3dk(k)    .gt. dfact * scr1(k)) then
	 vt3dg(k) = vt3dc(k) * scr1(k)
   endif
 elseif (vt3dc(k) .lt. 0.) then
   if (-vt3dg(k) * vt3dk(k+1) .gt. dfact * scr1(k+1)) then
	 vt3dg(k) = vt3dc(k) * scr1(k+1)
   endif
 endif

enddo

! Compute flux divergence

do k = 2,m1-1
    scr1(k) = scr1(k)  &
            + vt3dk(k) * ( vt3dg(k-1) - vt3dg(k) &
            + scp  (k) * ( vt3dc(k)   - vt3dc(k-1)))
enddo
return
end subroutine fa_zc_plumerise
!-------------------------------------------------------------------------------
!
subroutine advtndc_plumerise(m1,scp,sca,sct,dtl)
implicit none
integer :: m1,k
real :: dtl,dtli
real, dimension(m1) :: scp,sca,sct

dtli = 1. / dtl
do k = 2,m1-1
   sct(k) = sct(k) + (sca(k)-scp(k)) * dtli
enddo
return
end subroutine advtndc_plumerise
!-------------------------------------------------------------------------------
!
subroutine tend0_plumerise(this)
                        !!!, only: nm1,wt,tt,qvt,qct,qht,qit,sct
#include "PlumeRise.h"  
 wt(1:nm1)  = 0.
 tt(1:nm1)  = 0.
qvt(1:nm1)  = 0.
qct(1:nm1)  = 0.
qht(1:nm1)  = 0.
qit(1:nm1)  = 0.
vel_t(1:nm1)  = 0.
rad_t(1:nm1)  = 0.
!sct(1:nm1)  = 0.
end subroutine tend0_plumerise

!     ****************************************************************

subroutine scl_misc(this,m1)
implicit none
real, parameter :: g = 9.81, cp=1004.
integer m1,k
real dmdtm
#include "PlumeRise.h"

 do k=2,m1-1
      WBAR    = 0.5*(W(k)+W(k-1))  
!-- dry adiabat
      ADIABAT = - WBAR * G / CP 
!      
!-- entrainment     
      DMDTM = 2. * ALPHA * ABS (WBAR) / RADIUS (k)  != (1/M)DM/DT
      
!-- tendency temperature = adv + adiab + entrainment
      TT(k) = TT(K) + ADIABAT - DMDTM * ( T  (k) -    TE (k) ) 

!-- tendency water vapor = adv  + entrainment
      QVT(K) = QVT(K)         - DMDTM * ( QV (k) - QVENV (k) )

      QCT(K) = QCT(K)	      - DMDTM * ( QC (k)  )
      QHT(K) = QHT(K)	      - DMDTM * ( QH (k)  )
      QIT(K) = QIT(K)	      - DMDTM * ( QI (k)  )

!srf-04may2009
!-- tendency horizontal speed = adv  + entrainment
      VEL_T(K) = VEL_T(K)     - DMDTM * ( VEL_P (k) - VEL_E (k) )
!-- tendency horizontal speed = adv  + entrainment
      rad_t(K) = rad_t(K)     + 0.5*DMDTM*(6./5.)*RADIUS (k)
!srf-end

!-- tendency gas/particle = adv  + entrainment
!      SCT(K) = SCT(K)         - DMDTM * ( SC (k) -   SCE (k) )

enddo
end subroutine scl_misc

!     ****************************************************************
!srf-04may2009-AWE
subroutine scl_dyn_entrain(this,m1)
!use plumegen_coms
implicit none
real, parameter :: g = 9.81, cp=1004., pi=3.1416!, tunp =3.0
integer m1,k
real dmdtm
#include "PlumeRise.h"

 do k=2,m1-1
!      
!-- tendency horizontal radius from dyn entrainment
      rad_t(K) = rad_t(K)   + ABS((vel_e(k)-vel_p(k)))/pi

!-- entrainment     
      DMDTM = (2./3.1416)  *  ABS(VEL_E (k) - VEL_P (k)) / RADIUS (k)  

!-- tendency horizontal speed  from dyn entrainment
      VEL_T(K) = VEL_T(K)     - DMDTM * ( VEL_P (k) - VEL_E (k) )

!     if(VEL_P (k) - VEL_E (k) > 0.) cycle

!-- tendency temperature  from dyn entrainment
      TT(k) = TT(K)           - DMDTM * ( T (k) - TE  (k) ) 

!-- tendency water vapor  from dyn entrainment
      QVT(K) = QVT(K)         - DMDTM * ( QV (k) - QVENV (k) )

      QCT(K) = QCT(K)	      - DMDTM * ( QC (k)  )
      QHT(K) = QHT(K)	      - DMDTM * ( QH (k)  )
      QIT(K) = QIT(K)	      - DMDTM * ( QI (k)  )

!-- tendency gas/particle  from dyn entrainment
      !SCT(K) = SCT(K)         - DMDTM * ( SC (k) - SCE (k) )

enddo
end subroutine scl_dyn_entrain
!-srf-end
!     ****************************************************************

subroutine visc_W(this,m1,deltak,kmt)
implicit none
integer m1,k,deltak,kmt,m2
real dz1t,dz1m,dz2t,dz2m,d2wdz,d2tdz  ,d2qvdz ,d2qhdz ,d2qcdz ,d2qidz ,d2scdz&
    ,d2vel_pdz,d2rad_dz
#include "PlumeRise.h"

!srf--- 17/08/2005
!m2=min(m1+deltak,kmt)
m2=min(m1,kmt)

!do k=2,m1-1
do k=2,m2-1
 DZ1T   = 0.5*(ZT(K+1)-ZT(K-1))
 DZ2T   = VISC (k) / (DZ1T * DZ1T)  
 DZ1M   = 0.5*(ZM(K+1)-ZM(K-1))
 DZ2M   = VISC (k) / (DZ1M * DZ1M)  
 D2WDZ  = (W  (k + 1) - 2 * W  (k) + W  (k - 1) ) * DZ2M  
 D2TDZ  = (T  (k + 1) - 2 * T  (k) + T  (k - 1) ) * DZ2T  
 D2QVDZ = (QV (k + 1) - 2 * QV (k) + QV (k - 1) ) * DZ2T  
 D2QHDZ = (QH (k + 1) - 2 * QH (k) + QH (k - 1) ) * DZ2T 
 D2QCDZ = (QC (k + 1) - 2 * QC (k) + QC (k - 1) ) * DZ2T  
 D2QIDZ = (QI (k + 1) - 2 * QI (k) + QI (k - 1) ) * DZ2T  
 !D2SCDZ = (SC (k + 1) - 2 * SC (k) + SC (k - 1) ) * DZ2T 
 d2vel_pdz=(vel_P  (k + 1) - 2 * vel_P  (k) + vel_P  (k - 1) ) * DZ2T
 d2rad_dz =(rad_p  (k + 1) - 2 * rad_p  (k) + rad_p  (k - 1) ) * DZ2T
 
  WT(k) =   WT(k) + D2WDZ 
  TT(k) =   TT(k) + D2TDZ                          
 QVT(k) =  QVT(k) + D2QVDZ 
 QCT(k) =  QCT(k) + D2QCDZ 
 QHT(k) =  QHT(k) + D2QHDZ 
 QIT(k) =  QIT(k) + D2QIDZ     

 VEL_T(k) =   VEL_T(k) + d2vel_pdz

 rad_t(k) =   rad_t(k) + d2rad_dz
 !SCT(k) =  SCT(k) + D2SCDZ
 !print*,'W-VISC=',k,D2WDZ
enddo  

end subroutine visc_W

!     ****************************************************************
subroutine update_plumerise(this,m1,varn)
integer m1,k
character(len=*) :: varn
#include "PlumeRise.h"

if(varn == 'W') then

 do k=2,m1-1
   W(k) =  W(k) +  WT(k) * DT  
 enddo
 return

else 
do k=2,m1-1
   T(k) =  T(k) +  TT(k) * DT  

  QV(k) = QV(k) + QVT(k) * DT  

  QC(k) = QC(k) + QCT(k) * DT !cloud drops travel with air 
  QH(k) = QH(k) + QHT(k) * DT  
  QI(k) = QI(k) + QIT(k) * DT 
! SC(k) = SC(k) + SCT(k) * DT 

!srf---18jun2005  
  QV(k) = max(0., QV(k))
  QC(k) = max(0., QC(k))
  QH(k) = max(0., QH(k))
  QI(k) = max(0., QI(k))
  
!srf-04may2009
  VEL_P(k) =  VEL_P(k) + VEL_T(k) * DT  

  rad_p(k) =  rad_p(k) + rad_t(k) * DT  
! SC(k) = max(0., SC(k))

 enddo
endif
end subroutine update_plumerise
!-------------------------------------------------------------------------------
!
subroutine fallpart(this,m1)
integer m1,k
real vtc, dfhz,dfiz,dz1
!srf==================================
!   verificar se o gradiente esta correto 
!  
!srf==================================
!
!     XNO=1.E7  [m**-4] median volume diameter raindrop,Kessler
!     VC = 38.3/(XNO**.125), median volume fallspeed eqn., Kessler
!     for ice, see (OT18), use F0=0.75 per argument there. rho*q
!     values are in g/m**3, velocities in m/s

real, PARAMETER :: VCONST = 5.107387, EPS = 0.622, F0 = 0.75  
real, PARAMETER :: G = 9.81, CP = 1004.
#include "PlumeRise.h"

!
do k=2,m1-1

   VTC = VCONST * RHO (k) **.125   ! median volume fallspeed (KTable4)
                                
!  hydrometeor assembly velocity calculations (K Table4)
!  VTH(k)=-VTC*QH(k)**.125  !median volume fallspeed, water            
   VTH (k) = - 4.	    !small variation with qh
   
   VHREL = W (k) + VTH (k)  !relative to surrounding cloud
 
!  rain ventilation coefficient for evaporation
   CVH(k) = 1.6 + 0.57E-3 * (ABS (VHREL) ) **1.5  
!
!  VTI(k)=-VTC*F0*QI(k)**.125    !median volume fallspeed,ice             
   VTI (k) = - 3.                !small variation with qi

   VIREL = W (k) + VTI (k)       !relative to surrounding cloud
!
!  ice ventilation coefficient for sublimation
   CVI(k) = 1.6 + 0.57E-3 * (ABS (VIREL) ) **1.5 / F0  
!
!
   IF (VHREL.GE.0.0) THEN  
    DFHZ=QH(k)*(RHO(k  )*VTH(k  )-RHO(k-1)*VTH(k-1))/RHO(k-1)
   ELSE  
    DFHZ=QH(k)*(RHO(k+1)*VTH(k+1)-RHO(k  )*VTH(k  ))/RHO(k)
   ENDIF  
   !
   !
   IF (VIREL.GE.0.0) THEN  
    DFIZ=QI(k)*(RHO(k  )*VTI(k  )-RHO(k-1)*VTI(k-1))/RHO(k-1)
   ELSE  
    DFIZ=QI(k)*(RHO(k+1)*VTI(k+1)-RHO(k  )*VTI(k  ))/RHO(k)
   ENDIF
   
   DZ1=ZM(K)-ZM(K-1)
   
   qht(k) = qht(k) - DFHZ / DZ1 !hydrometeors don't
   		  
   qit(k) = qit(k) - DFIZ / DZ1  !nor does ice? hail, what about

enddo
end subroutine fallpart
!-------------------------------------------------------------------------------
!
subroutine printout (this,izprint,nrectotal)  
real, parameter :: tmelt = 273.3
integer, save :: nrec
data nrec/0/
integer :: ko,izprint,interval,nrectotal
real :: pea, btmp,etmp,vap1,vap2,gpkc,gpkh,gpki,deficit

real mass(nkp),xx(nkp),xxx, w_thresold

!ams real dummy(nxp,nkp)

integer :: k_initial, k_final, kk4, kl ! ams

#include "PlumeRise.h"


interval = 1              !debug time interval,min

!-srf-04-may-2009- AWE
!- define range of the upper detrainemnt layer
do ko=nkp-10,1,-1
 
 if(w(ko) < w_thresold) cycle
 
 if(k_final==0) k_final=ko
 
 if(w(ko)-1. > w(ko-1)) then
   k_initial=ko
   exit
  endif
  xxx=xxx+mass(ko)
enddo
!print*,'ki kf=',k_initial,k_final


!- if there is a non zero depth layer, make the mass vertical distribution 
if(k_final > 0 .and. k_initial > 0) then 
    
     k_initial=int((k_final+k_initial)*0.5)
    
    
     !- get the normalized mass distribution
     do ko=nkp-10,1,-1
        if(w(ko) < w_thresold) cycle
      
        if(w(ko)-1.0 > w(ko-1)) exit
        xx(ko) = mass(ko)/xxx
     enddo
   
    !v-2.5
    !- parabolic vertical distribution between k_initial and k_final
    xx=0.
    KK4 = k_final-k_initial+2
    do ko=1,kk4-1
      kl=ko+k_initial-1
      xx(kl) = 6.* float(ko)/float(kk4)**2 * (1. - float(ko)/float(kk4))
    enddo
    mass=xx
    !print*,'ki kf=',int(k_initial+k_final)/2,k_final,sum(mass)*100.
    
    !- check if the integral is 1 (normalized)
    if(sum(mass) .ne. 1.) then
 	xxx= ( 1.- sum(mass) )/float(k_final-k_initial+1)
 	do ko=k_initial,k_final
 	  mass(ko) = mass(ko)+ xxx !- values between 0 and 1.
 	enddo
      ! print*,'new mass=',sum(mass)*100.,xxx
      !pause
    endif
    
endif !k_final > 0 .and. k_initial > 


!- vertical mass distribution (horizontally integrated : kg/m)
!mass=3.14*rad_p**2*dne*sc*100.
!mass=mass/sum(mass)

!srf-040may2009-end

IF (IZPRINT.EQ.0) RETURN  


IF(MINTIME == 1) nrec = 0
!
WRITE (2, 430) MINTIME, DT, TIME  
WRITE (2, 431) ZTOP  
WRITE (2, 380)  
!
! do the print
!
 DO 390 KO = 1, nrectotal, interval  
                             
   PEA = PE (KO) * 10.       !pressure is stored in decibars(kPa),print in mb;
   BTMP = T (KO) - TMELT     !temps in Celsius
   ETMP = T (KO) - TE (KO)   !temperature excess
   VAP1 = QV (KO)   * 1000.  !printout in g/kg for all water,
   VAP2 = QSAT (KO) * 1000.  !vapor (internal storage is in g/g)
   GPKC = QC (KO)   * 1000.  !cloud water
   GPKH = QH (KO)   * 1000.  !raindrops
   GPKI = QI (KO)   * 1000.  !ice particles 
   DEFICIT = VAP2 - VAP1     !vapor deficit
!
   WRITE (2, 400) zt(KO)/1000., PEA, W (KO), BTMP, ETMP, VAP1, &
    VAP2, GPKC, GPKH, GPKI, VTH (KO), SC(KO)
!
!
!                                    !end of printout
   
  390 CONTINUE		  

   nrec=nrec+1
   write (19,rec=nrec) (W (KO), KO=1,nrectotal)
   nrec=nrec+1
   write (19,rec=nrec) (T (KO), KO=1,nrectotal)
   nrec=nrec+1
   write (19,rec=nrec) (TE(KO), KO=1,nrectotal)
   nrec=nrec+1
   write (19,rec=nrec) (QV(KO)*1000., KO=1,nrectotal)
   nrec=nrec+1
   write (19,rec=nrec) (QC(KO)*1000., KO=1,nrectotal)
   nrec=nrec+1
   write (19,rec=nrec) (QH(KO)*1000., KO=1,nrectotal)
   nrec=nrec+1
   write (19,rec=nrec) (QI(KO)*1000., KO=1,nrectotal)
   nrec=nrec+1
!   write (19,rec=nrec) (SC(KO), KO=1,nrectotal)
   write (19,rec=nrec) (QSAT(KO)*1000., KO=1,nrectotal)
   nrec=nrec+1
   write (19,rec=nrec) (QVENV(KO)*1000., KO=1,nrectotal)


 
!print*,'ntimes=',nrec/(9)
!
RETURN  
!
! ************** FORMATS *********************************************
!
  380 FORMAT(/,' Z(KM) P(MB) W(MPS) T(C)  T-TE   VAP   SAT   QC    QH' &
'     QI    VTH(MPS) SCAL'/)
!
  400 FORMAT(1H , F4.1,F7.2,F7.2,F6.1,6F6.2,F7.2,1X,F6.2)  
!
  430 FORMAT(1H ,//I5,' MINUTES       DT= ',F6.2,' SECONDS   TIME= ' &
        ,F8.2,' SECONDS')
  431 FORMAT(' ZTOP= ',F10.2)  
!
end subroutine printout
!
! *********************************************************************
SUBROUTINE WATERBAL(this)
#include "PlumeRise.h"  
!
                                        
IF (QC (L) .LE.1.0E-10) QC (L) = 0.  !DEFEAT UNDERFLOW PROBLEM
IF (QH (L) .LE.1.0E-10) QH (L) = 0.  
IF (QI (L) .LE.1.0E-10) QI (L) = 0.  
!
CALL EVAPORATE(this)    !vapor to cloud,cloud to vapor  
!                             
CALL SUBLIMATE(this)    !vapor to ice  
!                            
CALL GLACIATE(this)     !rain to ice 
                           
CALL MELT(this)         !ice to rain
!         
!if(ak1 > 0. .or. ak2 > 0.) &
CALL CONVERT (this) !(auto)conversion and accretion 
!CALL CONVERT2 () !(auto)conversion and accretion 
!

RETURN  
END SUBROUTINE WATERBAL
! *********************************************************************
SUBROUTINE EVAPORATE(this) 
!
!- evaporates cloud,rain and ice to saturation
!
implicit none
!
!     XNO=10.0E06
!     HERC = 1.93*1.E-6*XN035        !evaporation constant
!
real, PARAMETER :: HERC = 5.44E-4, CP = 1.004, HEATCOND = 2.5E3  
real, PARAMETER :: HEATSUBL = 2834., TMELT = 273., TFREEZE = 269.3

real, PARAMETER :: FRC = HEATCOND / CP, SRC = HEATSUBL / CP

real :: evhdt, evidt, evrate, evap, sd,	quant, dividend, divisor, devidt
#include "PlumeRise.h"  

!
!
SD = QSAT (L) - QV (L)  !vapor deficit
IF (SD.EQ.0.0)  RETURN  
!IF (abs(SD).lt.1.e-7)  RETURN  


EVHDT = 0.  
EVIDT = 0.  
!evrate =0.; evap=0.; sd=0.0; quant=0.0; dividend=0.0; divisor=0.0; devidt=0.0
                                 
EVRATE = ABS (WBAR * DQSDZ)   !evaporation rate (Kessler 8.32)
EVAP = EVRATE * DT            !what we can get in DT
                                  

IF (SD.LE.0.0) THEN  !     condense. SD is negative

   IF (EVAP.GE.ABS (SD) ) THEN    !we get it all
                                  
      QC (L) = QC  (L) - SD  !deficit,remember?
      QV (L) = QSAT(L)       !set the vapor to saturation  
      T  (L) = T   (L) - SD * FRC  !heat gained through condensation
                                !per gram of dry air
      RETURN  

   ELSE  
                                 
      QC (L) = QC (L) + EVAP         !get what we can in DT 
      QV (L) = QV (L) - EVAP         !remove it from the vapor
      T  (L) = T  (L) + EVAP * FRC   !get some heat

      RETURN  

   ENDIF  
!
ELSE                                !SD is positive, need some water
!
! not saturated. saturate if possible. use everything in order
! cloud, rain, ice. SD is positive
                                         
   IF (EVAP.LE.QC (L) ) THEN        !enough cloud to last DT  
!
                                         
      IF (SD.LE.EVAP) THEN          !enough time to saturate
                                         
         QC (L) = QC (L) - SD       !remove cloud                                          
         QV (L) = QSAT (L)          !saturate
         T (L) = T (L) - SD * FRC   !cool the parcel                                          
         RETURN  !done
!
                                         
      ELSE   !not enough time
                                        
         SD = SD-EVAP               !use what there is
         QV (L) = QV (L) + EVAP     !add vapor
         T (L) = T (L) - EVAP * FRC !lose heat
         QC (L) = QC (L) - EVAP     !lose cloud
	                            !go on to rain.                                      
      ENDIF     
!
   ELSE                !not enough cloud to last DT
!      
      IF (SD.LE.QC (L) ) THEN   !but there is enough to sat
                                          
         QV (L) = QSAT (L)  !use it
         QC (L) = QC (L) - SD  
         T  (L) = T (L) - SD * FRC  
         RETURN  
	                              
      ELSE            !not enough to sat
         SD = SD-QC (L)  
         QV (L) = QV (L) + QC (L)  
         T  (L) = T (L) - QC (L) * FRC         
         QC (L) = 0.0  !all gone
                                          
      ENDIF       !on to rain                           
   ENDIF          !finished with cloud
!
!  but still not saturated, so try to use some rain
!  this is tricky, because we only have time DT to evaporate. if there
!  is enough rain, we can evaporate it for dt. ice can also sublimate
!  at the same time. there is a compromise here.....use rain first, then
!  ice. saturation may not be possible in one DT time.
!  rain evaporation rate (W12),(OT25),(K Table 4). evaporate rain first
!  sd is still positive or we wouldn't be here.


   IF (QH (L) .LE.1.E-10) GOTO 33                                  

!srf-25082005
!  QUANT = ( QC (L)  + QV (L) - QSAT (L) ) * RHO (L)   !g/m**3
   QUANT = ( QSAT (L)- QC (L) - QV (L)   ) * RHO (L)   !g/m**3
!
   EVHDT = (DT * HERC * (QUANT) * (QH (L) * RHO (L) ) **.65) / RHO (L)
!             rain evaporation in time DT
                                         
   IF (EVHDT.LE.QH (L) ) THEN           !enough rain to last DT

      IF (SD.LE.EVHDT) THEN  		!enough time to saturate	  
         QH (L) = QH (L) - SD   	!remove rain	  
         QV (L) = QSAT (L)  		!saturate	  
         T (L) = T (L) - SD * FRC  	!cool the parcel		  
	 
	 RETURN  			!done
!                       
      ELSE                               !not enough time
         SD = SD-EVHDT  		 !use what there is
         QV (L) = QV (L) + EVHDT  	 !add vapor
         T (L) = T (L) - EVHDT * FRC  	 !lose heat
         QH (L) = QH (L) - EVHDT  	 !lose rain

      ENDIF  				  !go on to ice.
!                                    
   ELSE  !not enough rain to last DT
!
      IF (SD.LE.QH (L) ) THEN             !but there is enough to sat
         QV (L) = QSAT (L)                !use it
         QH (L) = QH (L) - SD  
         T (L) = T (L) - SD * FRC  
         RETURN  
!                            
      ELSE                              !not enough to sat
         SD = SD-QH (L)  
         QV (L) = QV (L) + QH (L)  
         T (L) = T (L) - QH (L) * FRC    
         QH (L) = 0.0                   !all gone
                                          
      ENDIF                             !on to ice
!
                                          
   ENDIF                                !finished with rain
!
!
!  now for ice
!  equation from (OT); correction factors for units applied
!
   33    continue
   IF (QI (L) .LE.1.E-10) RETURN            !no ice there
!
   DIVIDEND = ( (1.E6 / RHO (L) ) **0.475) * (SD / QSAT (L) &
            - 1) * (QI (L) **0.525) * 1.13
   DIVISOR = 7.E5 + 4.1E6 / (10. * EST (L) )  
                                                 
   DEVIDT = - CVI(L) * DIVIDEND / DIVISOR   !rate of change
                                                  
   EVIDT = DEVIDT * DT                      !what we could get
!
! logic here is identical to rain. could get fancy and make subroutine
! but duplication of code is easier. God bless the screen editor.
!
                                         
   IF (EVIDT.LE.QI (L) ) THEN             !enough ice to last DT
!
                                         
      IF (SD.LE.EVIDT) THEN  		  !enough time to saturate
         QI (L) = QI (L) - SD   	  !remove ice
         QV (L) = QSAT (L)  		  !saturate
         T (L) = T (L) - SD * SRC  	  !cool the parcel
	 
         RETURN  			  !done
!
                                          
      ELSE                                !not enough time
                                          
         SD = SD-EVIDT  		  !use what there is
         QV (L) = QV (L) + EVIDT  	  !add vapor
          T (L) =  T (L) - EVIDT * SRC  	  !lose heat
         QI (L) = QI (L) - EVIDT  	  !lose ice
                                          
      ENDIF  				  !go on,unsatisfied
!                                          
   ELSE                                   !not enough ice to last DT
!                                         
      IF (SD.LE.QI (L) ) THEN             !but there is enough to sat
                                          
         QV (L) = QSAT (L)                !use it
         QI (L) = QI   (L) - SD  
          T (L) =  T   (L) - SD * SRC  
	  
         RETURN  
!
      ELSE                                 !not enough to sat
         SD = SD-QI (L)  
         QV (L) = QV (L) + QI (L)  
         T (L) = T (L) - QI (L) * SRC             
         QI (L) = 0.0                      !all gone

      ENDIF                                !on to better things
                                           !finished with ice
   ENDIF  
!                                 
ENDIF                                      !finished with the SD decision
!
RETURN  
!
END SUBROUTINE EVAPORATE
!
! *********************************************************************
SUBROUTINE CONVERT (this)  
!
!- ACCRETION AND AUTOCONVERSION
!
!
real,      PARAMETER ::  AK1 = 0.001    !conversion rate constant
real,      PARAMETER ::  AK2 = 0.0052   !collection (accretion) rate
real,      PARAMETER ::  TH  = 0.5      !Kessler threshold
integer,   PARAMETER ::iconv = 1        !- Kessler conversion (=0)
                                     
!real, parameter :: ANBASE =  50.!*1.e+6 !Berry-number at cloud base #/m^3(maritime)
 real, parameter :: ANBASE =100000.!*1.e+6 !Berry-number at cloud base #/m^3(continental)
!real, parameter :: BDISP = 0.366       !Berry--size dispersion (maritime)
 real, parameter :: BDISP = 0.146       !Berry--size dispersion (continental)
real, parameter :: TFREEZE = 269.3  !ice formation temperature
!
real ::   accrete, con, q, h, bc1,   bc2,  total
#include "PlumeRise.h"  


IF (T (L)  .LE. TFREEZE) RETURN  !process not allowed above ice
!
IF (QC (L) .EQ. 0.     ) RETURN  

ACCRETE = 0.  
CON = 0.  
Q = RHO (L) * QC (L)  
H = RHO (L) * QH (L)  
!
!     selection rules
!                         
!            
IF (QH (L) .GT. 0.     ) ACCRETE = AK2 * Q * (H**.875)  !accretion, Kessler
!
IF (ICONV.NE.0) THEN   !select Berry or Kessler
!
!old   BC1 = 120.  
!old   BC2 = .0266 * ANBASE * 60.  
!old   CON = BDISP * Q * Q * Q / (BC1 * Q * BDISP + BC2) 	  

   CON = Q*Q*Q*BDISP/(60.*(5.*Q*BDISP+0.0366*ANBASE))
!
ELSE  
!                             
!   CON = AK1 * (Q - TH)   !Kessler autoconversion rate
!      
!   IF (CON.LT.0.0) CON = 0.0   !havent reached threshold
 
   CON = max(0.,AK1 * (Q - TH)) ! versao otimizada
!
ENDIF  
!
!
TOTAL = (CON + ACCRETE) * DT / RHO (L)  

!
IF (TOTAL.LT.QC (L) ) THEN  
!
   QC (L) = QC (L) - TOTAL  
   QH (L) = QH (L) + TOTAL    !no phase change involved
   RETURN  
!
ELSE  
!              
   QH (L) = QH (L) + QC (L)    !uses all there is
   QC (L) = 0.0  
!
ENDIF  
!
RETURN  
!
END SUBROUTINE CONVERT
!
!**********************************************************************
!
SUBROUTINE CONVERT2 (this)  
implicit none
LOGICAL  AEROSOL
parameter(AEROSOL=.true.)
!
real, parameter :: TNULL=273.16, LAT=2.5008E6 &
                  ,EPSI=0.622 ,DB=1. ,NB=1500. !ALPHA=0.2 
real :: KA,KEINS,KZWEI,KDREI,VT	
real :: A,B,C,D, CON,ACCRETE,total 
      
real Y(6),ROH
#include "PlumeRise.h"  
      
A=0.
B=0.
Y(1) = T(L)
Y(4) = W(L)
y(2) = QC(L)
y(3) = QH(L)
Y(5) = RADIUS(L)
ROH =  RHO(L)*1.e-3 ! dens (MKS) ??


! autoconversion

KA = 0.0005 
IF( Y(1) .LT. 258.15 )THEN
!   KEINS=0.00075
    KEINS=0.0009 
    KZWEI=0.0052
    KDREI=15.39
ELSE
    KEINS=0.0015
    KZWEI=0.00696
    KDREI=11.58
ENDIF
      
!   ROH=PE/RD/TE
VT=-KDREI* (Y(3)/ROH)**0.125

 
IF (Y(4).GT.0.0 ) THEN
 IF (AEROSOL) THEN
   A = 1/y(4)  *  y(2)*y(2)*1000./( 60. *( 5. + 0.0366*NB/(y(2)*1000.*DB) )  )
 ELSE
   IF (y(2).GT.(KA*ROH)) THEN
   !print*,'1',y(2),KA*ROH
   A = KEINS/y(4) *(y(2) - KA*ROH )
   ENDIF
 ENDIF
ELSE
   A = 0.0
ENDIF

! accretion

IF(y(4).GT.0.0) THEN
   B = KZWEI/(y(4) - VT) * MAX(0.,y(2)) *   &
       MAX(0.001,ROH)**(-0.875)*(MAX(0.,y(3)))**(0.875)
ELSE
   B = 0.0
ENDIF
   
   
      !PSATW=610.7*EXP( 17.25 *( Y(1) - TNULL )/( Y(1)-36. ) )
      !PSATE=610.7*EXP( 22.33 *( Y(1) - TNULL )/( Y(1)- 2. ) )

      !QSATW=EPSI*PSATW/( PE-(1.-EPSI)*PSATW )
      !QSATE=EPSI*PSATE/( PE-(1.-EPSI)*PSATE )
      
      !MU=2.*ALPHA/Y(5)

      !C = MU*( ROH*QSATW - ROH*QVE + y(2) )
      !D = ROH*LAT*QSATW*EPSI/Y1/Y1/RD *DYDX1

      
      !DYDX(2) = - A - B - C - D  ! d rc/dz
      !DYDX(3) = A + B            ! d rh/dz
 
 
      ! rc=rc+dydx(2)*dz
      ! rh=rh+dydx(3)*dz
 
CON      = A
ACCRETE  = B
 
TOTAL = (CON + ACCRETE) *(1/DZM(L)) /ROH     ! DT / RHO (L)  

!print*,'L=',L,total,QC(L),dzm(l)

!
IF (TOTAL.LT.QC (L) ) THEN  
!
   QC (L) = QC (L) - TOTAL  
   QH (L) = QH (L) + TOTAL    !no phase change involved
   RETURN  
!
ELSE  
!              
   QH (L) = QH (L) + QC (L)    !uses all there is
   QC (L) = 0.0  
!
ENDIF  
!
RETURN  
!
END SUBROUTINE CONVERT2
! ice - effect on temperature
!      TTD = 0.0 
!      TTE = 0.0  
!       CALL ICE(QSATW,QSATE,Y(1),Y(2),Y(3), &
!               TTA,TTB,TTC,DZ,ROH,D,C,TTD,TTE)
!       DYDX(1) = DYDX(1) + TTD  + TTE ! DT/DZ on Temp
!
!**********************************************************************
!
SUBROUTINE SUBLIMATE(this) 
!
! ********************* VAPOR TO ICE (USE EQUATION OT22)***************
!
real, PARAMETER :: EPS = 0.622, HEATFUS = 334., HEATSUBL = 2834., CP = 1.004
real, PARAMETER :: SRC = HEATSUBL / CP, FRC = HEATFUS / CP, TMELT = 273.3
real, PARAMETER :: TFREEZE = 269.3

real ::dtsubh,  dividend,divisor, subl
#include "PlumeRise.h"  
!
DTSUBH = 0.  
!
!selection criteria for sublimation
IF (T (L)  .GT. TFREEZE  ) RETURN  
IF (QV (L) .LE. QSAT (L) ) RETURN  
!
!     from (OT); correction factors for units applied
!
 DIVIDEND = ( (1.E6 / RHO (L) ) **0.475) * (QV (L) / QSAT (L) &
            - 1) * (QI (L) **0.525) * 1.13
 DIVISOR = 7.E5 + 4.1E6 / (10. * EST (L) )  
!
                                         
 DTSUBH = ABS (DIVIDEND / DIVISOR)   !sublimation rate
 SUBL = DTSUBH * DT                  !and amount possible
!
!     again check the possibilities
!
IF (SUBL.LT.QV (L) ) THEN  
!
   QV (L) = QV (L) - SUBL             !lose vapor
   QI (L) = QI (L) + SUBL  	      !gain ice
   T (L) = T (L) + SUBL * SRC         !energy change, warms air

	 !print*,'5',l,qi(l),SUBL

   RETURN  
!
ELSE  
!                                     
   QI (L) = QV (L)                    !use what there is
   T  (L) = T (L) + QV (L) * SRC      !warm the air
   QV (L) = 0.0  
	 !print*,'6',l,qi(l)
!
ENDIF  
!
RETURN  
END SUBROUTINE SUBLIMATE
!
! *********************************************************************
!
SUBROUTINE GLACIATE(this)
!
! *********************** CONVERSION OF RAIN TO ICE *******************
!     uses equation OT 16, simplest. correction from W not applied, but
!     vapor pressure differences are supplied.
!
!
real, PARAMETER :: HEATFUS = 334., CP = 1.004, EPS = 0.622, HEATSUBL = 2834.
real, PARAMETER :: FRC = HEATFUS / CP, FRS = HEATSUBL / CP, TFREEZE =  269.3
real, PARAMETER :: GLCONST = 0.025   !glaciation time constant, 1/sec
real dfrzh
#include "PlumeRise.h"  
!
                                    
 DFRZH = 0.    !rate of mass gain in ice
!
!selection rules for glaciation
IF (QH (L) .LE. 0.       ) RETURN  
IF (QV (L) .LT. QSAT (L) ) RETURN                                        
IF (T  (L) .GT. TFREEZE  ) RETURN  
!
!      NT=TMELT-T(L)
!      IF (NT.GT.50) NT=50
!
                                    
 DFRZH = DT * GLCONST * QH (L)    ! from OT(16)
!
IF (DFRZH.LT.QH (L) ) THEN  
!
   QI (L) = QI (L) + DFRZH  
   QH (L) = QH (L) - DFRZH  
   T (L) = T (L) + FRC * DFRZH  !warms air
   
   	 !print*,'7',l,qi(l),DFRZH

   
   RETURN  
!
ELSE  
!
   QI (L) = QI (L) + QH (L)  
   T  (L) = T  (L) + FRC * QH (L)  
   QH (L) = 0.0  

 !print*,'8',l,qi(l), QH (L)  
!
ENDIF  
!
RETURN  
!
END SUBROUTINE GLACIATE
!
!
! *********************************************************************
SUBROUTINE MELT(this)  
!
! ******************* MAKES WATER OUT OF ICE **************************
!                                              
real, PARAMETER :: FRC = 332.27, TMELT = 273., F0 = 0.75   !ice velocity factor
real DTMELT
#include "PlumeRise.h"  
!                                    
 DTMELT = 0.   !conversion,ice to rain
!
!selection rules
IF (QI (L) .LE. 0.0  ) RETURN  
IF (T (L)  .LT. TMELT) RETURN  
!
                                                      !OT(23,24)
 DTMELT = DT * (2.27 / RHO (L) ) * CVI(L) * (T (L) - TMELT) * ( (RHO(L)  &
         * QI (L) * 1.E-6) **0.525) * (F0** ( - 0.42) )
                                                      !after Mason,1956
!
!     check the possibilities
!
IF (DTMELT.LT.QI (L) ) THEN  
!
   QH (L) = QH (L) + DTMELT  
   QI (L) = QI (L) - DTMELT  
   T  (L) = T (L) - FRC * DTMELT     !cools air
   	 !print*,'9',l,qi(l),DTMELT

   
   RETURN  
!
ELSE  
!
   QH (L) = QH (L) + QI (L)   !get all there is to get
   T  (L) = T (L) - FRC * QI (L)  
   QI (L) = 0.0  
   	 !print*,'10',l,qi(l)
!
ENDIF  
!
RETURN  
!
!-----------------------------------------------------------------------------
END SUBROUTINE MELT
SUBROUTINE htint (nzz1, vctra, eleva, nzz2, vctrb, elevb)
  IMPLICIT NONE
  INTEGER, INTENT(IN ) :: nzz1
  INTEGER, INTENT(IN ) :: nzz2
  REAL,    INTENT(IN ) :: vctra(nzz1)
  REAL,    INTENT(OUT) :: vctrb(nzz2)
  REAL,    INTENT(IN ) :: eleva(nzz1)
  REAL,    INTENT(IN ) :: elevb(nzz2)

  INTEGER :: l
  INTEGER :: k
  INTEGER :: kk
  REAL    :: wt

  l=1

  DO k=1,nzz2
     DO
        IF ( (elevb(k) <  eleva(1)) .OR. &
             ((elevb(k) >= eleva(l)) .AND. (elevb(k) <= eleva(l+1))) ) THEN
           wt       = (elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
           vctrb(k) = vctra(l)+(vctra(l+1)-vctra(l))*wt
           EXIT
        ELSE IF ( elevb(k) >  eleva(nzz1))  THEN
           wt       = (elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
           vctrb(k) = vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
           EXIT
        END IF

        l=l+1
        IF(l == nzz1) THEN
           !PRINT *,'htint:nzz1',nzz1
           DO kk=1,l
              !PRINT*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
           END DO
           STOP 'htint'
        END IF
     END DO
     !print*, k, vctrb(k), wt, eleva(k), elevb(k), vctra(k)
  END DO
END SUBROUTINE htint

End Module PlumeRise_Mod
