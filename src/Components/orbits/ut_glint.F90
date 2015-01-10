      PROGRAM driver_glint 
          USE glintangle_mod

        IMPLICIT  NONE

        INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(15,307)
        REAL(dp), DIMENSION(1:3)  :: pos           ! satellite position vec (3x1) ECEF
        REAL(dp)  :: GLINT_ANGLE, sunzenith1, sunazimuth1 


!      ***********************************************************************************
!       !! Followings are inputs
!      ***********************************************************************************

! ! !       integer, PARAMETER :: nymd = 20081231 !date
! ! !       integer, PARAMETER :: nhms = 015057   !time
! ! 
! ! !        integer, PARAMETER :: nymd = 20030817 !date
! ! !        integer, PARAMETER :: nhms = 083030   !time

! 1- i=1 
       integer, PARAMETER :: nymd = 20091222 !date
       integer, PARAMETER :: nhms = 120057   !time
!2- j=7 i = 34
!        integer, PARAMETER :: nymd = 20091222 !date
!        integer, PARAMETER :: nhms = 120657   !time

! ! !        REAL(dp), PARAMETER  :: long = -105.1786   ! degrees pixel location
! ! !        REAL(dp), PARAMETER  :: lat  =  70.7425    ! degrees pixel location
! ! !        REAL(dp), PARAMETER  :: alt  =  1830.14    ! meters above mean sea level

! 1- i=1 
       REAL(dp), PARAMETER  :: longpix = 74.661635571182742   ! degrees pixel location
       REAL(dp), PARAMETER  :: latpix  =  68.067108483911937    ! degrees pixel location
       REAL(dp), PARAMETER  :: altpix  =  0    ! meters above mean sea level

!      POint on the earth (pixel) that will be used in the glint angle calculation
!2- j=7 i=34
!        REAL(dp), PARAMETER  :: longpix = 70.617638901981366   ! degrees pixel location
!        REAL(dp), PARAMETER  :: latpix  = -7.532203739000860   ! degrees pixel location
!        REAL(dp), PARAMETER  :: altpix  =  0    ! meters above mean sea level glint 100.34



! Satellite location
! 1- i=1
       REAL(dp), PARAMETER  :: tlon =   88.581539038744083 ! -0.046911860813809 * 180/pi       ! sat track lon degrees
       REAL(dp), PARAMETER  :: tlat  =  81.187756667752197 ! 1.221955472080043 * 180/pi        ! sat track lat degrees
       REAL(dp), PARAMETER  :: talt  =  722.9549514316926     ! km above mean sea level

! ! 2  i=34 j=7
!        REAL(dp), PARAMETER  :: tlon =  -2.687851633736419 ! -0.046911860813809 * 180/pi       ! sat track lon degrees
!        REAL(dp), PARAMETER  :: tlat  =  70.012891303102535 ! 1.221955472080043 * 180/pi        ! sat track lat degrees
!        REAL(dp), PARAMETER  :: talt  =  721.0658035813328     ! km above mean sea level



! Satellite location ___ ECEF coordinate system
       !DATA pos/99.7573, 1820.1013, -6842.42033/ ! ECEF
! 1- i=1
       DATA pos/0027.008796023224, 1090.741232272660, 6995.639314167769/ !ECEF
! 2  i=34 j=7
!        DATA pos/2430.365981369558,  -0114.096701326629,   6649.167424806347/ !ECEF  In km ... 

!*******************************************************************
!*******************************************************************


!      Others
       INTEGER :: i,j

       CALL calc_glint(GLINT_ANGLE, nymd, nhms, longpix, latpix, altpix, pos)
       print*, "1-", GLINT_ANGLE

       CALL calc_glint(GLINT_ANGLE, nymd, nhms, longpix, latpix, altpix, tlon, tlat, talt)
       print*, "2-", GLINT_ANGLE
!      CALL CALC_sun_zen_az(sunzenith1, sunazimuth1, nymd, nhms, long, lat, alt)
!      print*, "out ---> ", sunzenith1, sunazimuth1 

!   ***************************************************





! ! 
! !       !OPEN (UNIT=1, FILE='mask_mat') ! 74by99
! !       ! DO i = 1, size(mask,1)
! !       !    WRITE (1, "(99I2,1x)") (mask(i,j),j=1,size(mask,2)) 
! !       ! END DO
! !       !CLOSE(1,STATUS='KEEP') 
! ! 
!        OPEN (UNIT=2, FILE='m_lat') ! 74by99
!         DO i = 1, size(tlat,1)
!            WRITE (2, "(F30.10,1x)") (tlat(i)) 
!         END DO
!        CLOSE(2,STATUS='KEEP') 
!  
!        OPEN (UNIT=3, FILE='m_lon') ! 74by99
!         DO i = 1, size(tlon,1)
!            WRITE (3, "(F30.10,1x)") (tlon(i)) 
!         END DO
!        CLOSE(3,STATUS='KEEP') 
! !  
! ! 
! !       OPEN (UNIT=4, FILE='slon1') ! 74by99
! !        DO i = 1, size(slon,2)
! !           WRITE (4, "(F30.10,1x)") (slon(1,i)) 
! !        END DO
! !       CLOSE(4,STATUS='KEEP') 
! ! 
! !       OPEN (UNIT=5, FILE='slon2') ! 74by99
! !        DO i = 1, size(slon,2)
! !           WRITE (5, "(F30.10,1x)") (slon(2,i)) 
! !        END DO
! !       CLOSE(5,STATUS='KEEP') 
! ! 
! !       OPEN (UNIT=6, FILE='slon3') ! 74by99
! !        DO i = 1, size(slon,2)
! !           WRITE (6, "(F30.10,1x)") (slon(3,i)) 
! !        END DO
! !       CLOSE(6,STATUS='KEEP') 
! ! 
! ! 
! !       OPEN (UNIT=7, FILE='slat1') ! 74by99
! !        DO i = 1, size(slat,2)
! !           WRITE (7, "(F30.10,1x)") (slat(1,i)) 
! !        END DO
! !       CLOSE(7,STATUS='KEEP') 
! ! 
! !       OPEN (UNIT=8, FILE='slat2') ! 74by99
! !        DO i = 1, size(slat,2)
! !           WRITE (8, "(F30.10,1x)") (slat(2,i)) 
! !        END DO
! !       CLOSE(8,STATUS='KEEP') 
! ! 
! !       OPEN (UNIT=9, FILE='slat3') ! 74by99
! !        DO i = 1, size(slat,2)
! !           WRITE (9, "(F30.10,1x)") (slat(3,i)) 
! !        END DO
! !       CLOSE(9,STATUS='KEEP') 
! 



END PROGRAM driver_glint 

