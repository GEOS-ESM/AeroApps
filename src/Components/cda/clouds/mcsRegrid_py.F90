! pre-2016       original MCS version
! 020816    pmn  horizontal dims collapsed to 1D for SCS


 subroutine getEdgeCoords ( pe, ze, np, km, T, q, delp, zs, ptop )

   use mcsRegrid_Mod
   implicit NONE

   integer, intent(in) :: np              ! No. of pixels
   integer, intent(in) :: km              ! No. vertical layers
   real,    intent(in) :: T(np,km)        ! Dry temperature [K]
   real,    intent(in) :: q(np,km)        ! Specific humidity [kg/kg]
   real,    intent(in) :: delp(np,km)     ! layer pressure thicjness [Pa]
   real,    intent(in) :: zs(np)          ! model surface height [m]
   real,    intent(in) :: ptop            ! top (edge) pressure [Pa]

   real, intent(out)   :: pe(np,km+1)     ! pressure at edges
   real, intent(out)   :: ze(np,km+1)     ! height (above sea-level) at edges

!                   -------
   integer :: n

   do n = 1, np
     call getEdgeCoords1 ( pe(n,:), ze(n,:), km, T(n,:), q(n,:), delp(n,:), ptop )
     ze(n,:) = zs(n) + ze(n,:) ! now ze is above sea-level
   end do

 end subroutine getEdgeCoords

!---

  subroutine interpPressure ( di_pe, ps, np, nk, ks, &
                              di_ze, zs, g5_pe, g5_ze, km, rc ) 

     use mcsRegrid_Mod
     implicit NONE

     integer, intent(in)  :: np               ! No. of pixels
     integer, intent(in)  :: nk               ! number of DISORT layers
     integer, intent(in)  :: km               ! number of GEOS-5 layers

     real, intent(in)     :: zs(np)           ! Terrain height
     real, intent(in)     :: di_ze(nk+1)      ! DISORT edge height
     real, intent(in)     :: g5_pe(np,km+1)   ! GEOS-5 edge pressure
     real, intent(in)     :: g5_ze(np,km+1)   ! GEOS-5 edge height
 
     real, intent(out)    :: di_pe(np,nk+1)   ! DISORT edge pressure
     real, intent(out)    :: ps(np)           ! Surface pressure adjusted for Terrain height
     integer, intent(out) :: ks(np)           ! surface level 
     integer, intent(out) :: rc

!                          ----

   integer :: n

   do n = 1, np
     call interpPressure1 ( di_pe(n,:), nk, ps(n), ks(n),    &
                            di_ze(:), zs(n), g5_pe(n,:), &
                            g5_ze(n,:), km, rc )
   end do

   end subroutine interpPressure

!---
   subroutine edgeT ( di_Te, np, nk, ks, di_pe, g5_pe, g5_Tm, g5_Tskin, km, rc ) 

        use mcsRegrid_Mod, only: edgeT1
        use mod_consts, only: UNDEF
        implicit NONE

        integer, intent(in)  :: np      ! No. of pixels
        integer, intent(in)  :: nk      ! number of DISORT layers
        integer, intent(in)  :: km      ! number of GEOS-5 layers

        integer, intent(in)  :: ks(np)  ! DISORT surface level
  
        real, intent(in)     :: di_pe(np,nk+1)
        real, intent(in)     :: g5_pe(np,km+1)
        real, intent(in)     :: g5_Tm(np,km)
        real, intent(in)     :: g5_Tskin(np)

        real,    intent(out) :: di_Te(np,nk+1)
        integer, intent(out) :: rc

!                            ---
     integer :: n

     di_Te = UNDEF
     do n = 1, np
       call edgeT1 ( di_Te(n,:), nk, ks(n), di_pe(n,:), g5_pe(n,:), &
                     g5_Tm(n,:), g5_Tskin(n), km, rc ) 
       if ( rc /= 0 ) return
     end do

   end subroutine edgeT

!---
   subroutine edgeQ ( di_Qe, np, nk, ks, di_pe, g5_pe, g5_Qm, km, rc) 

        use mcsRegrid_Mod, only: edgeQ1
        use mod_consts, only: UNDEF
        implicit NONE

        integer, intent(in)  :: np      ! No. of pixels
        integer, intent(in)  :: nk      ! number of DISORT layers
        integer, intent(in)  :: km      ! number of GEOS-5 layers

        integer, intent(in)  :: ks(np)  ! DISORT surface level
  
        real, intent(in)     :: di_pe(np,nk+1)
        real, intent(in)     :: g5_pe(np,km+1)
        real, intent(in)     :: g5_Qm(np,km)

        real,    intent(out) :: di_Qe(np,nk+1)
        integer, intent(out) :: rc

!                            ---
     integer :: n

     rc = 0 ! for symmetry, always zero here
     di_Qe = UNDEF
     do n = 1, np
       call edgeQ1 ( di_Qe(n,:), nk, ks(n), di_pe(n,:), g5_pe(n,:), &
                     g5_Qm(n,:), km ) 
     end do

   end subroutine edgeQ

!---

   subroutine regridTracer ( di_q, di_pe, np, nk, ks, &
                             g5_q, g5_pe, km )

     use mcsRegrid_Mod, only: map_ppm
     use mod_consts, only: UNDEF

     implicit NONE

     integer, intent(in)  :: np               ! No. of pixels
     integer, intent(in)  :: nk               ! number of DISORT layers
     integer, intent(in)  :: km               ! number of GEOS-5 layers
     integer, intent(in)  :: ks(np)           ! lowest DISORT edge level 

     real, intent(in)     :: di_pe(np,nk+1)   ! DISORT edge pressure
     real, intent(in)     :: g5_q(np,km)      ! GEOS-5 tracer
     real, intent(in)     :: g5_pe(np,km+1)   ! GEOS-5 edge pressure
  
     real, intent(out)    :: di_q(np,nk)      ! interpolated tracer on DISORT grid  

!                            ---
     integer :: iv, kord, n, nk_

! pt = 1, 3, pk
! u  = 1, 3, pe
! q  = 0, 3, pe

     iv = 0
     kord = 3
     di_q = UNDEF
     do n = 1, np
       nk_ = ks(n) - 1
       call map_ppm ( di_q(n,:), nk_, di_pe(n,:), &
                      g5_pe(n,:), g5_q(n,:), km, iv, kord )
     end do

   end subroutine regridTracer

!---

   subroutine regridMet ( di_q, di_pe, np, nk, ks, &
                          g5_q, g5_pe, km )

     use mcsRegrid_Mod, only: map_ppm
     use mod_consts, only: UNDEF

     implicit NONE

     integer, intent(in)  :: np               ! No. of pixels
     integer, intent(in)  :: nk               ! number of DISORT layers
     integer, intent(in)  :: km               ! number of GEOS-5 layers
     integer, intent(in)  :: ks(np)           ! lowest DISORT edge level 

     real, intent(in)     :: di_pe(np,nk+1)   ! DISORT edge pressure (or proxy)
     real, intent(in)     :: g5_q(np,km)      ! GEOS-5 tracer
     real, intent(in)     :: g5_pe(np,km+1)   ! GEOS-5 edge pressure (or proxy)
  
     real, intent(out)    :: di_q(np,nk)      ! interpolated tracer on DISORT grid  

!                            ---
     integer :: iv, kord, n, nk_

! pt = 1, 3, pk
! u  = 1, 3, pe
! q  = 0, 3, pe

     iv = 1
     kord = 3
     do n = 1, np
       nk_ = ks(n) - 1
       call map_ppm ( di_q(n,:), nk_, di_pe(n,:), &
                      g5_pe(n,:), g5_q(n,:), km, iv, kord )
     end do

   end subroutine regridMet
