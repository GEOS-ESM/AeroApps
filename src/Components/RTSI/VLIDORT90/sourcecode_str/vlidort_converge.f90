! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6                                  #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #              VLIDORT_CONVERGE                               #
! #                                                             #
! ###############################################################


      MODULE vlidort_converge_module

      PRIVATE
      PUBLIC :: VLIDORT_CONVERGE

      CONTAINS

      SUBROUTINE VLIDORT_CONVERGE ( &
        AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM, &
        DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, &
        DO_SS_EXTERNAL, & ! New 15 March 2012
        DO_SSFULL, DO_DOUBLE_CONVTEST, &
        DO_RAYLEIGH_ONLY, DO_UPWELLING, &
        NSTOKES, NSTREAMS, &
        NLAYERS, VLIDORT_ACCURACY, &
        N_USER_RELAZMS, N_USER_LEVELS, &
        DO_TOA_CONTRIBS, &
        DO_ALL_FOURIER, DO_DBCORRECTION, &
        DO_NO_AZIMUTH, N_CONVTESTS, &
        LOCAL_UM_START, N_DIRECTIONS, &
        WHICH_DIRECTIONS, &
        N_OUT_STREAMS, VZA_OFFSETS, &
        STOKES_SS, STOKES_DB, SS_CONTRIBS, &
        STOKES_F, MS_CONTRIBS_F, IB, &
        STOKES, FOURIER_SAVED, CONTRIBS, &
        TESTCONV, LOCAL_ITERATION )

!  convergence test on the Stokes vector intensity

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::           LOCAL_N_USERAZM
      DOUBLE PRECISION, INTENT (IN) ::  AZMFAC &
          ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )
      LOGICAL, INTENT (IN) ::           DO_SSCORR_NADIR
      LOGICAL, INTENT (IN) ::           DO_SSCORR_OUTGOING
!  New 15 March 2012
      LOGICAL, INTENT (IN) ::           DO_SS_EXTERNAL
      LOGICAL, INTENT (IN) ::           DO_SSFULL
      LOGICAL, INTENT (IN) ::           DO_DOUBLE_CONVTEST
      LOGICAL, INTENT (IN) ::           DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      DOUBLE PRECISION, INTENT (IN) ::  VLIDORT_ACCURACY
      INTEGER, INTENT (IN) ::           N_USER_RELAZMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS
      LOGICAL, INTENT (IN) ::           DO_ALL_FOURIER
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH
      INTEGER, INTENT (IN) ::           N_CONVTESTS
      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      INTEGER, INTENT (IN) ::           N_DIRECTIONS
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::           N_OUT_STREAMS
      INTEGER, INTENT (IN) ::           VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::  STOKES_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  SS_CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::  MS_CONTRIBS_F &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
      INTEGER, INTENT (IN) ::           IB

      DOUBLE PRECISION, INTENT (INOUT) :: STOKES &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      INTEGER, INTENT (INOUT) ::          FOURIER_SAVED ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (INOUT) :: CONTRIBS &
          ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
      INTEGER, INTENT (INOUT) ::          TESTCONV
      LOGICAL, INTENT (INOUT) ::          LOCAL_ITERATION

!  local variables
!  ---------------

      INTEGER ::          COUNT, COUNT_A, V, N
      INTEGER ::          I, IDIR, UT, UA, W, O
      DOUBLE PRECISION :: TNEW ( MAXSTOKES ), ACCUR, &
                          TOLD ( MAXSTOKES ), TAZM ( MAXSTOKES )

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on STOKES = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on STOKES = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED, do not get calculated)

!  Full single scatter calculation is initialized to zero here (Version

        IF ( .not. DO_SSFULL ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = LOCAL_UM_START, N_OUT_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = VZA_OFFSETS(IB,I) + UA
                  DO O = 1, NSTOKES
                    STOKES(UT,V,O,W) = STOKES_F(UT,I,IB,O,W)
                  ENDDO
                ENDDO
               ENDDO
            ENDDO
          ENDDO
        ELSE
          FOURIER_SAVED(IB) = FOURIER_COMPONENT
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = LOCAL_UM_START, N_OUT_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = VZA_OFFSETS(IB,I) + UA
                  DO O = 1, NSTOKES
                    STOKES(UT,V,O,W) = ZERO
                  ENDDO
                ENDDO
               ENDDO
            ENDDO
          ENDDO
        ENDIF

!  TOA contribution functions (only if flagged)

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( .not. DO_SSFULL ) THEN
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IB,I) + UA
                DO O = 1, NSTOKES
                  DO N = 1, NLAYERS
                    CONTRIBS(V,O,N) = MS_CONTRIBS_F(I,IB,O,N)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IB,I) + UA
                DO O = 1, NSTOKES
                  DO N = 1, NLAYERS
                    CONTRIBS(V,O,N) = ZERO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!    STOKES. Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        IF ( DO_SSFULL .OR. &
             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN
         DO IDIR = 1, N_DIRECTIONS
          W = WHICH_DIRECTIONS(IDIR)
          DO UT = 1, N_USER_LEVELS
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS (IB,I) + UA
                DO O = 1, NSTOKES
!         if (v.eq.1)write(*,*)UT,STOKES(UT,V,O,W),STOKES_SS(UT,V,1,W)
                  STOKES(UT,V,O,W) = &
                    STOKES(UT,V,O,W) + STOKES_SS(UT,V,O,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
         ENDDO
        ENDIF
!        pause

!    CONTRIBS. Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS (IB,I) + UA
                DO O = 1, NSTOKES
                  DO N = 1, NLAYERS
                    CONTRIBS(V,O,N) = &
                      CONTRIBS(V,O,N) + SS_CONTRIBS(V,O,N)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
!        pause

!    STOKES. Add the direct beam component if flagged (upwelling only)
!       Convergence on STOKES = STOKES + DBEXACT
!    Full single scatter option added Version 2.3

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
        IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) &
             .AND.DO_UPWELLING ) THEN
          DO UT = 1, N_USER_LEVELS
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS(IB,I) + UA
                DO O = 1, NSTOKES
                  STOKES(UT,V,O,UPIDX) = &
                    STOKES(UT,V,O,UPIDX) + STOKES_DB(UT,V,O)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!  For each azimuth, add Fourier component

          DO UA = 1, LOCAL_N_USERAZM

!     - for direction, user optical depth, out stream

            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = LOCAL_UM_START, N_OUT_STREAMS
                  V = VZA_OFFSETS(IB,I) + UA
                  DO O = 1, NSTOKES
                    TOLD(O) = STOKES(UT,V,O,W)
                    TAZM(O) = AZMFAC(I,IB,UA,O)*STOKES_F(UT,I,IB,O,W)
                    STOKES(UT,V,O,W) = TOLD(O) + TAZM(O)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

          ENDDO

!  Examine convergence on intensity only
!  -------------------------------------

!  convergence test applied to ALL directions AND
!                              ALL stream values (except near zenith) AN
!                              ALL azimuths taken together
!                              ALL user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy lev

          COUNT = 0
          DO UA = 1, N_USER_RELAZMS
            COUNT_A = 0
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = LOCAL_UM_START, N_OUT_STREAMS
                  V = VZA_OFFSETS(IB,I) + UA
                  DO O = 1, NSTOKES
                    TOLD(O) = STOKES(UT,V,O,W)
                    TAZM(O) = AZMFAC(I,IB,UA,O)*STOKES_F(UT,I,IB,O,W)
                    TNEW(O) = TOLD(O) + TAZM(O)
                  ENDDO
                  IF ( TAZM(1) .NE. ZERO ) THEN
                    ACCUR     = DABS(TAZM(1)/TNEW(1))
                    IF ( ACCUR .LT. VLIDORT_ACCURACY ) THEN
                      COUNT   = COUNT + 1
                      COUNT_A = COUNT_A + 1
                    ENDIF
                  ELSE
                    COUNT   = COUNT + 1
                    COUNT_A = COUNT_A + 1
                  ENDIF
                  DO O = 1, NSTOKES
                    STOKES(UT,V,O,W)     = TNEW(O)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                  LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
                LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IB) = FOURIER_COMPONENT
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IB) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF

!  TOA_CONTRIBS: For each azimuth, add Fourier component

        IF ( DO_TOA_CONTRIBS ) THEN
          DO UA = 1, LOCAL_N_USERAZM
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              V = VZA_OFFSETS(IB,I) + UA
              DO N = 1, NLAYERS
                DO O = 1, NSTOKES
                  TOLD(O) = CONTRIBS(V,O,N)
                  TAZM(O) = AZMFAC(I,IB,UA,O)*MS_CONTRIBS_F(I,IB,O,N)
                  CONTRIBS(V,O,N) = TOLD(O) + TAZM(O)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER_COMPONENT .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IB) = FOURIER_COMPONENT
          ENDIF
        ENDIF

!  For all Fourier, keep saveing the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IB) = FOURIER_COMPONENT
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CONVERGE

      END MODULE vlidort_converge_module

