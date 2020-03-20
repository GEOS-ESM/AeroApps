      PROGRAM ut_sgp4
        
        USE sgp4_mod
        USE TLE_mod 

        IMPLICIT  NONE

!      Test dates 1-
       integer, dimension(1:2), PARAMETER :: nymd = (/20090323, 20090430/) ! start/end date
       integer, dimension(1:2), PARAMETER :: nhms = (/081500, 095500/)     ! start/end time

       integer, parameter :: deltat = 60 
       real(r8), dimension(2), parameter :: swath_size = (/100, 100/)  ! swath in km [left right]
                                           ! on the going direction left right 
       LOGICAL                           :: wrapon 

       INTEGER   :: rc  ! error code
       REAL(r8), DIMENSION(:), POINTER :: tlat, tlon, latlinear, lonlinear, obstimelinear 
       REAL(r8), DIMENSION(:,:), POINTER :: slat, slon ! swath coordinates

       INTEGER :: i,j

       integer :: Sat_num

       type (TLE) :: p

       logical :: flagr = .TRUE.
                
        
!      Input Sat name or number 
!      ------------------------
#if 0
      Character(len=12) :: Sat_name = 'SGP4-VER.TLE'
      Character(len=12) :: Sat_name = 'AQUA.TLE'      !Sat-1 
      Character(len=12) :: Sat_name = 'CALIPSO.TLE'   !Sat-2
      Character(len=12) :: Sat_name = 'CLOUDSAT.TLE'  !Sat-3
      Character(len=12) :: Sat_name = 'AURA.TLE'      !Sat-4
      Character(len=12) :: Sat_name = 'TERRA.TLE'     !Sat-5
#endif

      Character(len=*), parameter :: Sat_name = 'TLE/QUIKSCAT.TLE'   !Sat-6

#if 0
        !character(len=20), PARAMETER :: Sat_name = "Aqua"
        !character(len=20), PARAMETER :: Sat_name = "Calipso"
        !character(len=20), PARAMETER :: Sat_name = "CloudSat"
        !character(len=20), PARAMETER :: Sat_name = "Aura"
        !character(len=20), PARAMETER :: Sat_name = "Terra"
#endif

        CALL SGP4_Track(tlon, tlat, Sat_name, nymd, nhms, deltat, rc )
        CALL SGP4_Swath(slon, slat, Sat_name, nymd, nhms, deltat, swath_size, rc, wrapon)

        OPEN (UNIT=2, FILE='sgp4_track.txt') !
         DO i = 1, size(tlat,1)
            WRITE (2,*) tlon(i), tlat(i) 
         END DO
        CLOSE(2,STATUS='KEEP') 

    END PROGRAM ut_sgp4












