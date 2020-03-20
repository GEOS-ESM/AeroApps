! **********************************************
!
!  Simple f77 wrapper for the Python interface.
!
! **********************************************

subroutine plume ( km, u, v, T, q, delp, ptop, hflux_kW, area, & 
                   p, z, k, rc)

   use PlumeRise_Mod

   implicit none

!  !ARGUMENTS:

   integer, intent(in)  :: km             ! number of vertical layers

   real,    intent(in)  :: u(km)          ! zonal wind (m/s)
   real,    intent(in)  :: v(km)          ! meridional wind (m/s)
   real,    intent(in)  :: T(km)          ! potential temperature [K]
   real,    intent(in)  :: q(km)          ! specific humidity [kg/kg]
   real,    intent(in)  :: delp(km)       ! pressure thickness [Pa]
   real,    intent(in)  :: ptop           ! top edge pressure [Pa]
   real,    intent(in)  :: area           ! fire area [m^2]
   real,    intent(in)  :: hflux_kW       ! fire heat flux [kW/m2]
   
   real,    intent(out) :: p              ! upper plume pressure
   real,    intent(out) :: z              ! upper plume height
   integer, intent(out) :: k              ! upper plume vertical index
   integer, intent(out) :: rc             ! error code

 !                       ----

!  Local variables
!  ---------------
   type(PlumeRise) :: pr
   
!                                ----

!  Initialize
!  ----------
   call PlumeRise_Initialize(pr)

!  Run
!  ---
   call PlumeRise_Run (pr, km,     &
                       u, v, T, q, delp, ptop, &
                       area, hflux_kW,   &
                       z, p, k )

!  Finalize
!  ---------
   call PlumeRise_Finalize(pr)

   rc = 0

 end subroutine plume

subroutine plumeVMD ( km, u, v, T, q, delp, ptop, hflux_kW, area, & 
                      z_i, z_d, z_a, z_f, &
                      z_plume, w_plume, rc)

   use PlumeRise_Mod

   implicit none

!  !ARGUMENTS:

   integer, intent(in)  :: km             ! number of vertical layers

   real,    intent(in)  :: u(km)          ! zonal wind (m/s)
   real,    intent(in)  :: v(km)          ! meridional wind (m/s)
   real,    intent(in)  :: T(km)          ! potential temperature [K]
   real,    intent(in)  :: q(km)          ! specific humidity [kg/kg]
   real,    intent(in)  :: delp(km)       ! pressure thickness [Pa]
   real,    intent(in)  :: ptop           ! top edge pressure [Pa]
   real,    intent(in)  :: area           ! fire area [m^2]
   real,    intent(in)  :: hflux_kW       ! fire heat flux [kW/m2]
   
   real,    intent(out) :: z_i  ! height of maximum W (bottom of plume)
   real,    intent(out) :: z_d  ! height of maximum detrainment
   real,    intent(out) :: z_a  ! average height in (z_i,z_f), weighted by -dw/dz
   real,    intent(out) :: z_f  ! height where w<1 (top of plume)

   integer, parameter :: kp = 200
   real, intent(out) :: z_plume(kp) ! vertical levels
   real, intent(out) :: w_plume(kp) ! vertical velocity

   integer, intent(out) :: rc

 !                       ----

!  Local variables
!  ---------------
   type(PlumeRise) :: pr
   
!                                ----


!  Consistency check
!  -----------------
   if ( nkp .ne. kp ) then
        rc = 1
        return
   end if

!  Initialize
!  ----------
   call PlumeRise_Initialize(pr)

!  Run
!  ---
   call PlumeRise_Run (pr, km,     &
                       u, v, T, q, delp, ptop, &
                       area, hflux_kW,   &
                       z_i, z_d, z_a, z_f, &
                       z_plume, w_plume )

!
!  Find bottom, mid and top plume height. If using the parabolic VMD as in getVMD() below,
!  there are several options:
!
!  a) Like Saulo:
!                     z_c   = (z_f+z_i)/2
!                     delta = (z_f-z_i)/2
!
!  b) Preserve bottom half:
!                     z_c   = z_d
!                     delta = z_d - z_i
!
!  c) Preserve upper half:
!                     z_c   = z_d
!                     delta = z_f - z_d
!
!                       ---


!  Finalize
!  ---------
   call PlumeRise_Finalize(pr)

   rc = 0

 end subroutine plumeVMD

!..................................................................................

subroutine getVMD(km,z,z_c,delta,v) 

!
! Computes normalized vertical mass distribution (VMD) given z_c and delta. Assumes a gaussian.
!
   implicit NONE
   integer, intent(in)  :: km             ! number of vertical layers
   real,    intent(in)  :: z(km) ! height above surface
   real,    intent(in)  :: z_c   ! level of maximum detrainment (center of plume)
   real,    intent(in)  :: delta ! width of vertical mass distribution
   real,    intent(out) :: v(km) ! normalized vertical mass distribution
!                   ---

!  Analytic function
!  -----------------
   where ( abs(z-z_c) <= delta )
         v = 3.*(delta**2 - (z-z_c)**2)/(4.*delta**3) 
   elsewhere
         v = 0.0
   end where

!  Ensure normalization on grid
!  ----------------------------
   v = v / sum(v)

 end subroutine getVMD

!..................................................................................

subroutine biome ( km, u, v, T, q, delp, ptop, area, ibiome, & 
                   p1, p2, z1, z2, k1, k2, rc)

   use PlumeRise_Mod

   implicit none

!  !ARGUMENTS:

   integer, intent(in)  :: km             ! number of vertical layers

   real,    intent(in)  :: u(km)          ! zonal wind (m/s)
   real,    intent(in)  :: v(km)          ! meridional wind (m/s)
   real,    intent(in)  :: T(km)          ! potential temperature [K]
   real,    intent(in)  :: q(km)          ! specific humidity [kg/kg]
   real,    intent(in)  :: delp(km)       ! pressure thickness [Pa]
   real,    intent(in)  :: ptop           ! top edge pressure [Pa]
   real,    intent(in)  :: area           ! fire area [m^2]
   integer, intent(in)  :: ibiome         ! biome index
   
   real,    intent(out) :: p1             ! lower pressure
   real,    intent(out) :: p2             ! upper pressure
   real,    intent(out) :: z1             ! lower plume height
   real,    intent(out) :: z2             ! upper plume height
   integer, intent(out) :: k1             ! upper plume vertical index
   integer, intent(out) :: k2             ! lower plume vertical index
   integer, intent(out) :: rc             ! error code

 !                       ----

!  Local variables
!  ---------------
   type(PlumeRise) :: pr
   
   real    :: p1_p(N_BIOME),p2_p(N_BIOME) 
   real    :: z1_p(N_BIOME),z2_p(N_BIOME) 
   integer :: k1s(N_BIOME), k2s(N_BIOME)
   real    :: areas(N_BIOME)

!                                ----

!  Biome check
!  -----------
   if ( ibiome < 1 .OR. ibiome > N_BIOME ) then
      rc = 1
      return
   else
      areas = 0.0
      areas(ibiome) = area
   end if

!   print *, '----'
!   print *, 'f2py: km, ibiome, ptop, areas = ', km, ibiome, ptop, area 
!   print *, 'f2py:      T min/max   = ', minval(T), maxval(T)
!   print *, 'f2py:      q min/max   = ', minval(q), maxval(q)
!   print *, 'f2py:   delp min/max   = ', minval(delp), maxval(delp)

!  Initialize
!  ----------
   call PlumeRise_Initialize(pr)

!  Run
!  ---
   call PlumeRise_Run (pr, km, u, v, T, q, delp, ptop, areas, &
                       z1_plume=z1_p, z2_plume=z2_p,    &
                       p1_plume=p1_p ,p2_plume=p2_p,    &
                       k1_plume=k1s, k2_plume=k2s)
                           
!  Finalize
!  ---------
   call PlumeRise_Finalize(pr)

   p1 = p1_p(ibiome)  ! bottom pressure
   p2 = p2_p(ibiome)  ! top pressure
   z1 = z1_p(ibiome)
   z2 = z2_p(ibiome)
   k1 = k1s(ibiome) 
   k2 = k2s(ibiome)  
   rc = 0

!  write(*,'(a,2F10.2,2I4,2F10.2)') 'f2py: p, k, z = ', p1, p2, k1, k2, z1, z2 

 end subroutine Biome

!..................................................................................

  subroutine kPBL( k_pbl, km, T, q, delp, ptop, pblh )

   use rconstants
   implicit NONE

   integer, intent(in) :: km                 ! No. vertical layers
   real,    intent(in) :: T(km)              ! Dry temperature [K]
   real,    intent(in) :: q(km)              ! Specific humidity [kg/kg]
   real,    intent(in) :: delp(km)           ! layer pressure thickness [Pa]
   real,    intent(in) :: ptop               ! top (edge) pressure [Pa]
   real,    intent(in) :: pblh               ! PBL height [m]
   
   integer, intent(out) :: k_pbl             ! vertical layer index of PBL

!
! Returns the vertical layer index of the PBL.
!

!
! ------------------------------- pe(k),   he(k)
!
! ............................... pm(k),   hm(k), delp(k)
!
! ------------------------------- pe(k+1), he(k+1)
!

   integer :: k
   real :: Tv, kappa = rocp
   real :: delh(km), pe(km+1), he(km+1), mixr(km)

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
      delh(k)  = Rgas * Tv * log(pe(k+1)/pe(k)) / g
    end do

!  Compute Geo-potential height at edges
!  -------------------------------------
   he(km+1) = 0.0  ! set the surface to zero height
   do k = km, 1, -1
      he(k) = he(k+1) + delh(k)
   end do

!  Now find the PBL layer
!  ----------------------
   k_pbl = -1  ! so that we can check later for abnormal behavior 
   do k = km, 1, -1
      if ( he(k) >= pblh ) then
           k_pbl = k
           exit
      end if
   end do

 end subroutine kPBL
