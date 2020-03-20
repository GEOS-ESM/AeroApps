module LockPlume_Mod
! <OVERVIEW>
!
!      Code pulled from LockEntrain.F90 in GEOS-5 which implements a
!      K-PROFILE BOUNDARY LAYER SCHEME WITH CLOUD TOP ENTRAINMENT
!      Adapted from code provided by Stephen A. Klein at GFDL
!
!      The LockEntrain moduke calculates diffusivity coefficients for vertical
!      diffusion using a K-profile approach.  This scheme is modelled
!      after:
!
!      Lock, A.P., A.R. Brown, M.R. Bush, G.M. Martin, and R.N.B. Smith, 
!          2000: A new boundary layer mixing scheme. Part I: Scheme 
!          description and single-column modeling tests. Mon. Wea. Rev.,
!          128, 3187-3199.
!
!   
! </OVERVIEW>
!
   use qsat_mod, only: dqsat

   implicit none

!  Constants from MAPL, here to avoid MAPL dependence
!  --------------------------------------------------
   real, parameter, public :: MAPL_GRAV   = 9.80                   ! m^2/s
   real, parameter, public :: MAPL_KAPPA  = 2.0/7.0                ! --
   real, parameter, public :: MAPL_RUNIV  = 8314.3                 ! J/(Kmole K)
   real, parameter, public :: MAPL_AIRMW  = 28.97                  ! kg/Kmole
   real, parameter, public :: MAPL_RGAS   = MAPL_RUNIV/MAPL_AIRMW  ! J/(kg K)
   real, parameter, public :: MAPL_CP     = MAPL_RGAS/MAPL_KAPPA   ! J/(kg K)
   real, parameter, public :: MAPL_ALHL   = 2.4665E6               ! J/kg @15C

   private

!-----------------------------------------------------------------------
!
!  public interfaces

   PUBLIC mpbl_depth
   PUBLIC derived


!-----------------------------------------------------------------------
!
!  set default values to some instance independent parameters       
!

   real, parameter :: akmax      =  1.e4 ! maximum value for a diffusion coefficient 
                                         ! (m2/s)

   real, parameter :: zcldtopmax =  3.e3 ! maximum altitude for cloud top of 
                                         ! radiatively driven convection (m)    

   real, parameter :: ramp       =  20.


!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!
! Subroutines include:
!
!     entrain          main driver program of the module
!
!     mpbl_depth       routine to calculate the depth of surface driven
!                      mixed layer
!
!     radml_depth      subroutine to calculate the depth of the cloud
!                      topped radiatively driven mixed layer
!     
!     diffusivity_pbl2 subroutine to calculate diffusivity coefficients
!                      for surface driven mixed layer

contains

!======================================================================= 

! GPU The GPU main routine call is smaller due to CUDA limit on
!     number of arguments permitted in call. Most inputs and outputs
!     are USE-associated in the GridComp

!======================================================================= 
!
!  Subroutine to calculate pbl depth
!
   subroutine mpbl_depth(i,ncol,nlev,tpfac, entrate, pceff, t, q, u, v, z, p, b_star, u_star , ipbl, ztop )
!
!  -----
!  INPUT
!  -----
!
!  i             column number
!  ncol          total number of columns
!  nlev          number of levels
!  t             temperature (K)
!  q             specific humidity (g/g)
!  u             zonal wind (m/s)
!  v             meridional wind (m/s)
!  b_star        buoyancy scale (m s-2)
!  u_star        surface velocity scale (m/s)
!       
!  ------
!  OUTPUT
!  ------
!
!  ipbl          half level containing pbl height
!  ztop          pbl height (m)

      integer, intent(in   )                       :: i, nlev, ncol
      real,    intent(in   ), dimension(ncol,nlev) :: t, z, q, p, u, v
      real,    intent(in   ), dimension(ncol)      :: b_star, u_star
      real,    intent(in   )                       :: tpfac, entrate, pceff
      integer, intent(  out)                       :: ipbl
      real,    intent(  out),dimension(ncol)       :: ztop


      real     :: tep,z1,z2,t1,t2,qp,pp,qsp,dqp,dqsp,u1,v1,u2,v2,du
      real     :: ws,k_t_ref,entfr, tpfac_x, entrate_x,vscale
      integer  :: k


      !real, dimension(nlev) :: qst ! Not used in this code?


!calculate surface parcel properties

      tep  = t(i,nlev)
      qp   = q(i,nlev)

!--------------------------------------------
! wind dependence of plume character. 
! 
!    actual_entrainment_rate_at_z  ~ entrate * [ 1.0 +  |U(z)-U(0)| / vscale ]
! 
! entrate:  tunable param from rc file
! vscale:   tunable param hardwired here.

!!Old Shear  vscale    = 5.0 ! m s-1
      vscale    = 0.25 / 100. ! change of .25 m s-1 in 100 m
!!vscale    = 0.10 / 100. ! change of .10 m s-1 in 100 m


!---------------------------------------------
! tpfac scales up bstar by inv. ratio of
! heat-bubble area to stagnant area

      tep  = tep * (1.+ tpfac * b_star(i)/MAPL_GRAV)
!!  tep  = tep * (1.+ tpfac * b_star(i)*u_star(i)/MAPL_GRAV)

!search for level where this is exceeded              

      t1   = t(i,nlev)
      v1   = v(i,nlev)
      u1   = u(i,nlev)
      z1   = z(i,nlev)
      ztop(i) = z1
      do k = nlev-1 , 2, -1
         z2 = z(i,k)
         t2 = t(i,k)
         u2 = u(i,k)
         v2 = v(i,k)
         pp = p(i,k)

!!Old Shear     du = sqrt ( ( u(i,k) - u1 )**2 + ( v(i,k) - v1 )**2 )
         du = sqrt ( ( u2 - u1 )**2 + ( v2 - v1 )**2 ) / (z2-z1)
         du = min(du,1.0e-8)

         entrate_x = entrate * ( 1.0 + du / vscale )

         entfr= min( entrate_x*(z2-z1), 0.99 )

         qp  = qp  + entfr*(q(i,k)-qp)

! dry adiabatic ascent through one layer.
! Static energy conserved. 
         tep = tep - MAPL_GRAV*( z2-z1 )/MAPL_CP

! Environmental air entrained
         tep = tep + entfr*(t(i,k)-tep)

         dqsp = dqsat(tep , pp , qsat=qsp,  pascals=.true. )

         dqp = max( qp - qsp, 0. )/(1.+(MAPL_ALHL/MAPL_CP)*dqsp )
         qp  = qp - dqp
         tep = tep  + pceff * MAPL_ALHL * dqp/MAPL_CP  ! "Precipitation efficiency" basically means fraction
! of condensation heating that gets applied to parcel


! If parcel temperature (tep) colder than env (t2)
! OR if entrainment too big, declare this the PBL top
         if ( (t2 .ge. tep) .or. ( entfr .ge. 0.9899 ) ) then
!!a0082  Make the parcel a little less buoyant (liquid water loading)- 
!!a0082     translates here into subtracting from tep
!!a0082   if ( (t2 .ge. (tep-0.2) ) .or. ( entfr .ge. 0.9899 ) ) then
            ztop(i) = 0.5*(z2+z1)
            ipbl = k+1
            exit
         end if

         z1 = z2
         t1 = t2
         u1 = u2
         v1 = v2
      enddo

      return

   end subroutine mpbl_depth

  subroutine Derived( pm, hm, km, T, q, delp, ptop )

   implicit NONE

   integer, intent(in) :: km                 ! No. vertical layers
   real,    intent(in) :: T(km)              ! Dry temperature [K]
   real,    intent(in) :: q(km)              ! Specific humidity [kg/kg]
   real,    intent(in) :: delp(km)           ! layer pressure thicjness [Pa]
   real,    intent(in) :: ptop               ! top (edge) pressure [Pa]

   real, intent(out)  :: pm(km)              ! pressure at mid-layers
   real, intent(out)  :: hm(km)              ! height at mid-layers

!
! ------------------------------- pe(k),   he(k)
!
! ............................... pm(k),   hm(k), delp(k)
!
! ------------------------------- pe(k+1), he(k+1)
!

   real :: mixr(km)      ! mixing ratio
   real :: pe(km+1)      ! pressure at edges
   real :: he(km+1)      ! height at edges

   real :: Tv, kappa = MAPL_KAPPA
   real :: delh(km), delhm(km), pek(km+1)
   integer :: k

!  Mixing ratio from specific humidity
!  -----------------------------------
   mixr = q / ( 1.0 - q )

!  Construct edge pressures
!  ------------------------
   pe(1) = ptop
   do k = 1, km
      pe(k+1) = pe(k) + delp(k)
   end do

!  Construct mid-layer pressures and layer thickness
!  -------------------------------------------------
   do k = 1, km
      Tv = T(k) * ( 1 + 0.61 * mixr(k) )
      pm(k)    = 0.5 * ( pe(k) + pe(k+1) )
      delh(k)  = MAPL_RGAS * Tv * log(pe(k+1)/pe(k)) / MAPL_GRAV
      delhm(k) = MAPL_RGAS * Tv * log(pe(k+1)/pm(k)) / MAPL_GRAV
   end do

!  Compute Geo-potential height at edges and mid-layer
!  ---------------------------------------------------
   he(km+1) = 0.0  ! set the surface to zero height
   do k = km, 1, -1
      he(k) = he(k+1) + delh(k)
      hm(k) = he(k+1) + delhm(k)
   end do

 end subroutine Derived

end module LockPlume_Mod
