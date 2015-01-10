  subroutine getAeroPBL ( QABL, QPBL, &
                          Q, PE, HE, PBLH, km, jm, im)

    implicit NONE
    integer, intent(in) :: km, jm,im
    real,    intent(in) :: Q(km,jm,im)
    real,    intent(in) :: PE(km+1,jm,im)
    real,    intent(in) :: HE(km+1,jm,im)
    real,    intent(in) :: PBLH(jm,im)

    real,   intent(out) :: QABL(jm,im)
    real,   intent(out) :: QPBL(jm,im)

!      Returns average concentration above and below PBL.


!                           ----

    real*8  :: abl, pbl
    real    :: DZ(km), RHO(km)
    real    :: MAPL_GRAV = 9.8
    integer :: i, j, k

    i_loop: do i = 1, im
     j_loop: do j = 1, jm
       QABL(j,i) = 0.0
       QPBL(j,i) = 0.0
       abl = 0.0
       pbl = 0.0
       do k = 1, km
          DZ(k)  = (HE(k,j,i)-HE(k+1,j,i))
          RHO(k) = (PE(k+1,j,i)-PE(k,j,i))/(MAPL_GRAV*DZ(k))
       end do
       k_loop: do k = 1, km
          if ( (HE(k,j,i)-HE(km,j,i)) >= PBLH(j,i) ) then
               QABL(j,i) = QABL(j,i) + Q(k,j,i) * RHO(k) * DZ(k)
               abl = abl + DZ(k)
          else
               QPBL(j,i) = QPBL(j,i) + Q(k,j,i) * RHO(k) * DZ(k)
               pbl = pbl + DZ(k)
          end if
       end do k_loop
       if ( abl>0. ) QABL(j,i) = QABL(j,i) / abl ! averare concentration
       if ( pbl>0. ) QPBL(j,i) = QPBL(j,i) / pbl ! average concentration
     end do j_loop
    end do i_loop

    QABL = 1.E9 * QABL
    QPBL = 1.E9 * QPBL

  end subroutine getAeroPBL

  subroutine getAeroLev ( QLEV, Q, PE, HE, CTP, km, jm, im)

    use MAPL_ConstantsMod, only: MAPL_GRAV
    implicit NONE
    integer, intent(in) :: km, jm,im
    real,    intent(in) :: Q(km,jm,im)
    real,    intent(in) :: PE(km+1,jm,im)
    real,    intent(in) :: HE(km+1,jm,im)
    real,    intent(in) :: CTP(jm,im)

    real,   intent(out) :: QLEV(jm,im)

!      Returns concentration at cloud top
!      Note: On input, CTP assumed to be in Pa.

!                           ----

    integer :: i, j, k
    real    :: RHO

    i_loop: do i = 1, im
     j_loop: do j = 1, jm

       QLEV(j,i) = 0.0
       k_loop: do k = 1, km
          if ( (CTP(j,i)>PE(k,j,i)) .AND. (CTP(j,i)<=PE(k+1,j,i) ) ) then 
               RHO  = -(PE(k+1,j,i)-PE(k,j,i))/(HE(k+1,j,i)-HE(k,j,i))
               QLEV(j,i) = Q(k,j,i) * RHO / MAPL_GRAV
               exit k_loop
          end if
       end do k_loop

     end do j_loop
    end do i_loop

    QLEV = 1.E9 * QLEV

  end subroutine getAeroLev
  

  subroutine reviseCTP ( rCTP, CTP, CTT, &
                         PE, HE, PBLH, km, jm, im)

    use MAPL_ConstantsMod, only: MAPL_GRAV, MAPL_RGAS, MAPL_CP, MAPL_KAPPA
    implicit NONE
    integer, intent(in) :: km, jm,im
    real,    intent(in) :: CTP(jm,im)
    real,    intent(in) :: CTT(jm,im)
    real,    intent(in) :: PE(km+1,jm,im)
    real,    intent(in) :: HE(km+1,jm,im)
    real,    intent(in) :: PBLH(jm,im)

    real,   intent(out) :: rCTP(jm,im)

!   Revises cloud top pressure (CTP) in Collection 5.
!   The idea here is to retain good CO2 slicing retrievals
!   (above 550 hPA) but rederive the IR based CTP.
!
!   Note: On input, CTP assumed to be in Pa.

!                           ----

    real*8  :: pblTHETA
    real    :: DZ, DPK, THETA
    real    :: gcp = MAPL_GRAV/MAPL_CP
    integer :: i, j, k, n

    i_loop: do i = 1, im
     j_loop: do j = 1, jm

       rCTP(j,i) = CTP(j,i)
       if ( CTP(j,i) < 55000. ) then
            cycle j_loop
       end if

       n = 0
       pblTHETA = 0.0
       k_loop: do k = km, 1, -1

          DZ  =  HE(k,j,i) - HE(k+1,j,i)
          DPK = (PE(k+1,j,i)**MAPL_KAPPA) - (PE(k,j,i)**MAPL_KAPPA)
          THETA = gcp * DZ/DPK

          if ( (HE(k,j,i)-HE(km,j,i)) >= PBLH(j,i) ) then
               if ( k==km ) then
                  pblTHETA = THETA
                  n = 1
               end if
               exit k_loop
          end if
          pblTHETA = pblTHETA + THETA
          n = n + 1
       end do k_loop

       pblTHETA = pblTHETA / n

       rCTP(j,i) = (CTT(j,i)/pblTHETA)**(1./MAPL_KAPPA)

       rCTP(j,i) = min(rCTP(j,i),PE(km+1,j,i)) ! not below the surface 
            
     end do j_loop
    end do i_loop

  end subroutine reviseCTP

  subroutine reviseCTP2 ( rCTP, CTP, CTT, &
                          PE, HE, km, jm, im)

    use MAPL_ConstantsMod, only: MAPL_GRAV, MAPL_RGAS
    implicit NONE
    integer, intent(in) :: km, jm,im
    real,    intent(in) :: CTP(jm,im)
    real,    intent(in) :: CTT(jm,im)
    real,    intent(in) :: PE(km+1,jm,im)
    real,    intent(in) :: HE(km+1,jm,im)

    real,   intent(out) :: rCTP(jm,im)

!   Revises cloud top pressure (CTP) in Collection 5.
!   The idea here is to retain good CO2 slicing retrievals
!   (above 550 hPA) but rederive the IR based CTP.
!
!   Temperature match from below.
!
!   Note: On input, CTP assumed to be in Pa.

!                           ----

    real    :: DZ, DLP, TV
    real    :: factor = MAPL_GRAV/MAPL_RGAS
    integer :: i, j, k

    i_loop: do i = 1, im
     j_loop: do j = 1, jm

       rCTP(j,i) = CTP(j,i)
       if ( CTP(j,i) < 55000. ) then
            cycle j_loop
       end if

       k_loop: do k = km, 1, -1

          DZ  =  HE(k,j,i) - HE(k+1,j,i)
          DLP =  log(PE(k+1,j,i)) - log(PE(k,j,i))
          TV = factor * DZ/DLP

          if ( CTT(j,i) >= TV  ) then
             rCTP(j,i) = 0.5 * (PE(k,j,i)+PE(k+1,j,i))
             exit k_loop
          end if

       end do k_loop

     end do j_loop
    end do i_loop

  end subroutine reviseCTP2

  subroutine getRH(RH,QV,T,P,jm,im)
    use MAPL_SatVaporMod
    implicit NONE

    integer, intent(in) :: jm, im
    real, intent(in)    :: QV(jm,im) 
    real, intent(in)    :: T(jm,im) 
    real, intent(in)    :: P(jm,im) 
    
    real, intent(out)   :: RH(jm,im)

!             -----
    integer :: i, j

    i_loop: do i = 1, im
       j_loop: do j = 1, jm

          RH(j,i) = QV(j,i) / MAPL_EQsat(T(j,i),P(j,i))
    
       end do j_loop
    end do i_loop

  end subroutine getRH
