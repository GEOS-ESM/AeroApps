
 subroutine getEdgeCoords ( pe, ze, ns, nf, km, T, q, delp, zs, ptop )

   use mcsRegrid_Mod
   implicit NONE

   integer, intent(in) :: ns, nf             ! No. scans, frames
   integer, intent(in) :: km                 ! No. vertical layers
   real,    intent(in) :: T(ns,nf,km)        ! Dry temperature [K]
   real,    intent(in) :: q(ns,nf,km)        ! Specific humidity [kg/kg]
   real,    intent(in) :: delp(ns,nf,km)     ! layer pressure thicjness [Pa]
   real,    intent(in) :: zs(ns,nf)          ! model surface height [m]
   real,    intent(in) :: ptop               ! top (edge) pressure [Pa]

   real, intent(out)   :: pe(ns,nf,km+1)     ! pressure at edges
   real, intent(out)   :: ze(ns,nf,km+1)     ! height (above sea-level) at edges

!                   -------
   integer :: i, j

   do j = 1, nf
      do i = 1, ns
         call getEdgeCoords1 ( pe(i,j,:), ze(i,j,:), km, T(i,j,:), q(i,j,:), delp(i,j,:), ptop )
         ze(i,j,:) = zs(i,j) + ze(i,j,:) ! now ze is above sea-level
      end do
   end do

 end subroutine getEdgeCoords

!---

  subroutine interpPressure ( di_pe, ps, ns, nf, nk, ks, &
                              di_ze, zs, g5_pe, g5_ze, km, rc ) 

     use mcsRegrid_Mod
     implicit NONE

     integer, intent(in)  :: ns, nf              ! No. scans, frames
     integer, intent(in)  :: nk                  ! number of DISORT layers
     integer, intent(in)  :: km                  ! number of GEOS-5 layers

     real, intent(in)     :: zs(ns,nf)           ! Terrain height
     real, intent(in)     :: di_ze(nk+1)         ! DISORT edge height
     real, intent(in)     :: g5_pe(ns,nf,km+1)   ! GEOS-5 edge pressure
     real, intent(in)     :: g5_ze(ns,nf,km+1)   ! GEOS-5 edge height
 
     real, intent(out)    :: di_pe(ns,nf,nk+1)   ! DISORT edge pressure
     real, intent(out)    :: ps(ns,nf)           ! Surface pressure adjusted for Terrain height
     integer, intent(out) :: ks(ns,nf)    ! surface level 
     integer, intent(out) :: rc

!                          ----

   integer :: i, j

   do j = 1, nf
      do i = 1, ns
         call interpPressure1 ( di_pe(i,j,:), nk, ps(i,j), ks(i,j),    &
                                di_ze(:), zs(i,j), g5_pe(i,j,:), &
                                g5_ze(i,j,:), km, rc )
      end do
   end do

   end subroutine interpPressure

!---
   subroutine edgeT ( di_Te, ns, nf, nk, ks, di_pe, g5_pe, g5_Tm, g5_Tskin, km, rc ) 

        use mcsRegrid_Mod, only: edgeT1
        use mod_consts, only: UNDEF
        implicit NONE

        integer, intent(in)  :: ns         ! No. scans, frames
        integer, intent(in)  :: nf         ! No. scans, frames
        integer, intent(in)  :: nk         ! number of DISORT layers
        integer, intent(in)  :: km         ! number of GEOS-5 layers

        integer, intent(in)  :: ks(ns,nf)  ! DISORT surface level
  
        real, intent(in)     :: di_pe(ns,nf,nk+1)
        real, intent(in)     :: g5_pe(ns,nf,km+1)
        real, intent(in)     :: g5_Tm(ns,nf,km)
        real, intent(in)     :: g5_Tskin(ns,nf)

        real,    intent(out) :: di_Te(ns,nf,nk+1)
        integer, intent(out) :: rc

!                            ---
     integer :: i, j, nk_

     di_Te = UNDEF
     do j = 1, nf
        do i = 1, ns
           nk_ = ks(i,j) - 1
           call edgeT1 ( di_Te(i,j,:), nk, ks(i,j), di_pe(i,j,:), g5_pe(i,j,:), &
                         g5_Tm(i,j,:), g5_Tskin(i,j), km, rc ) 
           if ( rc /= 0 ) return
        end do
     end do

   end subroutine edgeT

!---
   subroutine edgeQ ( di_Qe, ns, nf, nk, ks, di_pe, g5_pe, g5_Qm, km, rc) 

        use mcsRegrid_Mod, only: edgeQ1
        use mod_consts, only: UNDEF
        implicit NONE

        integer, intent(in)  :: ns         ! No. scans, frames
        integer, intent(in)  :: nf         ! No. scans, frames
        integer, intent(in)  :: nk         ! number of DISORT layers
        integer, intent(in)  :: km         ! number of GEOS-5 layers

        integer, intent(in)  :: ks(ns,nf)  ! DISORT surface level
  
        real, intent(in)     :: di_pe(ns,nf,nk+1)
        real, intent(in)     :: g5_pe(ns,nf,km+1)
        real, intent(in)     :: g5_Qm(ns,nf,km)

        real,    intent(out) :: di_Qe(ns,nf,nk+1)
        integer, intent(out) :: rc

!                            ---
     integer :: i, j, nk_

     rc = 0 ! for symmetry, always zero here
     di_Qe = UNDEF
     do j = 1, nf
        do i = 1, ns
           nk_ = ks(i,j) - 1
           call edgeQ1 ( di_Qe(i,j,:), nk, ks(i,j), di_pe(i,j,:), g5_pe(i,j,:), &
                         g5_Qm(i,j,:), km ) 
        end do
     end do

   end subroutine edgeQ

!---

   subroutine regridTracer ( di_q, di_pe, ns, nf, nk, ks, &
                             g5_q, g5_pe, km )

     use mcsRegrid_Mod, only: map_ppm
     use mod_consts, only: UNDEF

     implicit NONE

     integer, intent(in)  :: ns, nf              ! No. scans, frames
     integer, intent(in)  :: nk                  ! number of DISORT layers
     integer, intent(in)  :: km                  ! number of GEOS-5 layers
     integer, intent(in)  :: ks(ns,nf)           ! lowest DISORT edge level 

     real, intent(in)     :: di_pe(ns,nf,nk+1)   ! DISORT edge pressure
     real, intent(in)     :: g5_q(ns,nf,km)      ! GEOS-5 tracer
     real, intent(in)     :: g5_pe(ns,nf,km+1)   ! GEOS-5 edge pressure
  
     real, intent(out)    :: di_q(ns,nf,nk)      ! interpolated tracer on DISORT grid  

!                            ---
     integer :: iv, kord, i, j, nk_

! pt = 1, 3, pk
! u  = 1, 3, pe
! q  = 0, 3, pe

     iv = 0
     kord = 3
     di_q = UNDEF
     do j = 1, nf
        do i = 1, ns
           nk_ = ks(i,j) - 1
           call map_ppm ( di_q(i,j,:), nk_, di_pe(i,j,:), &
                          g5_pe(i,j,:), g5_q(i,j,:), km, iv, kord )
        end do
     end do

   end subroutine regridTracer

!---

   subroutine regridMet ( di_q, di_pe, ns, nf, nk, ks, &
                          g5_q, g5_pe, km )

     use mcsRegrid_Mod, only: map_ppm
     use mod_consts, only: UNDEF

     implicit NONE

     integer, intent(in)  :: ns, nf              ! No. scans, frames
     integer, intent(in)  :: nk                  ! number of DISORT layers
     integer, intent(in)  :: km                  ! number of GEOS-5 layers
     integer, intent(in)  :: ks(ns,nf)           ! lowest DISORT edge level 

     real, intent(in)     :: di_pe(ns,nf,nk+1)   ! DISORT edge pressure (or proxy)
     real, intent(in)     :: g5_q(ns,nf,km)      ! GEOS-5 tracer
     real, intent(in)     :: g5_pe(ns,nf,km+1)   ! GEOS-5 edge pressure (or proxy)
  
     real, intent(out)    :: di_q(ns,nf,nk)      ! interpolated tracer on DISORT grid  

!                            ---
     integer :: iv, kord, i, j, nk_

! pt = 1, 3, pk
! u  = 1, 3, pe
! q  = 0, 3, pe

     iv = 1
     kord = 3
     do j = 1, nf
        do i = 1, ns
           nk_ = ks(i,j) - 1
           call map_ppm ( di_q(i,j,:), nk_, di_pe(i,j,:), &
                          g5_pe(i,j,:), g5_q(i,j,:), km, iv, kord )
        end do
     end do

   end subroutine regridMet
