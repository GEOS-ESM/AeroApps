! file: mod_reff.f90
! Effective Radii as used by GEOS-5 (Dec 2012)

module mod_reff

   implicit none
 
   private
   public reff, get_tempor

contains

   ! tempor variable needed by reff:
   ! Defaults to 0.
   ! Set to one if zonal wind speed is stronger than 4 m/s eastward
   ! in any of layers containing reference pressure 925 hPa and below

   subroutine get_tempor (lm, pref, u, tempor)

      implicit none

      ! pref = reference air pressure (vert, edges) [Pa]
      ! u = zonal (easward) wind speed (vert, centers) [m/s]

      integer, intent(in) :: lm
      real, intent(in) :: pref(lm+1)
      real, intent(in) :: u(lm)
      real, intent(out) :: tempor
      integer :: l
     
      tempor = 0.
      do l = max(1,count(pref < 92500.)), lm
         if (u(l) .gt. 4.) then
           tempor = 1.
           exit
         endif
      end do

   end subroutine get_tempor

   subroutine REFF( &
         TE,        &  ! [K]
         PL,        &  ! [hPa]
         QCiLS,     &  ! [kg/kg]
         QCiAN,     &  ! [kg/kg]
         TEMPOR,    &  ! see get_tempor()
         RAD_RL,    &  ! [m]
         RAD_RI)       ! [m]

      use mod_consts, only : MIN_RI, MAX_RI, RI_ANV
      implicit none

      real, intent(in ) :: TE, PL, QCiAN, QCiLS, TEMPOR
      real, intent(out) :: RAD_RL, RAD_RI

      real :: RAD_RI_AN

      !!!!!!!!!!!
      !!! ICE !!!
      !!!!!!!!!!!

      ! basic pressure form
      if (PL < 150.) then
         RAD_RI = MAX_RI
      else
         RAD_RI = MAX_RI*150./PL
      end if

      ! weighted anvil reff
      RAD_RI_AN = RAD_RI
      if (QCiLS + QCiAN > 0.0) then
         RAD_RI_AN =(QCiLS + QCiAN) / (QCiLS/RAD_RI + QCiAN/RI_ANV)
      end if

      ! choose the smallest
      RAD_RI = MIN( RAD_RI, RAD_RI_AN )
      RAD_RI = MAX( RAD_RI, MIN_RI )

      !!!!!!!!!!!!!!
      !!! LIQUID !!!
      !!!!!!!!!!!!!!

      ! basic pressure form
      if (PL < 300.) then
         RAD_RL = 21.e-6
      else
         RAD_RL = 21.e-6*300./PL
      end if
      RAD_RL = MAX( RAD_RL, 10.e-6 )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Thicken low high lat clouds !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( PL .GE. 775.  .AND. TE .LE.  275. .AND. (tempor.eq.1.) ) then
         RAD_RL = max(min(-0.1 * PL + 87.5, 10.),5.)*1.e-6
      end if
      if ( PL .GE. 825.  .AND. TE .LE.  282. .AND. (tempor.eq.1.) ) then
         RAD_RL = max(0.71 * TE - 190.25, 5.)*1.e-6
      end if
      if ( PL .GE. 775.  .AND. PL .LT. 825. .AND. TE .LE.  282. .AND. TE .GT. 275. .AND. (tempor.eq.1.) ) then
         RAD_RL = min(-0.1*PL + 0.71 * TE - 107.75, 10.)*1.e-6
      end if
      if ( PL .GE. 825.  .AND. TE .LE.  275. .AND. (tempor.eq.1.) ) then
         RAD_RL = 5.*1.e-6
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Thin low tropical clouds !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
         RAD_RL = min(2.2 * TE - 617., 21.)*1.e-6
      end if
      if ( PL .GE. 925.  .AND. TE .GE.  290. ) then
         RAD_RL = min(0.44 * PL - 397., 21.)*1.e-6
      end if
      if ( PL .GE. 925.  .AND. PL .LT. 950. .AND. TE .GT.  285. .AND. TE .LT. 290.) then
         RAD_RL = max(min(0.44*PL + 2.2 * TE - 1035., 21.),10.)*1.e-6
      end if
      if ( PL .GE. 950.  .AND. TE .GE.  290. ) then
         RAD_RL = 21.*1.e-6
      end if

   end subroutine REFF

end module mod_reff
