C ###############################################################
C #                                                             #
C #                    THE VECTOR LIDORT MODEL                  #
C #                                                             #
C #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
C #   -      --         -        -        -         -           #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, inc.                          #
C #            9 Channing Street                                #
C #             Cambridge, MA 02138, USA                        #
C #            Tel: (617) 492 1183                              #
C #                                                             #
C #  Email :      rtsolutions@verizon.net                       #
C #                                                             #
C #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT           #
C #  Release Date :   December 2005  (2.0)                      #
C #  Release Date :   March 2007     (2.2)                      #
C #  Release Date :   October 2007   (2.3)                      #
C #  Release Date :   December 2008  (2.4)                      #
C #  Release Date :   April/May 2009 (2.4R)                     #
C #  Release Date :   July 2009      (2.4RT)                    #
C #                                                             #
C #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
C #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
C #       NEW: Thermal Emission Treatment     (2.4RT)           #
C #                                                             #
C ###############################################################

C    #####################################################
C    #                                                   #
C    #   This Version of VLIDORT comes with a GNU-style  #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #              VLIDORT_CONVERGE                               #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_CONVERGE
     I      ( AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM,
     O        TESTCONV, IB, LOCAL_ITERATION )

C  convergence test on the Stokes vector intensity

C  Include files
C  -------------

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables 

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include files of already calculated results

      INCLUDE '../includes/VLIDORT_SINGSCAT.VARS'

C  Include file of final results

      INCLUDE '../includes/VLIDORT_RESULTS.VARS'

C  input variables

      INTEGER          FOURIER_COMPONENT, LOCAL_N_USERAZM
      DOUBLE PRECISION AZMFAC
     &     (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS,MAXSTOKES)

C  modified/output variables

      LOGICAL          LOCAL_ITERATION
      INTEGER          TESTCONV, IB

C  local variables
C  ---------------

C  local variables

      INTEGER          COUNT, COUNT_A, V
      INTEGER          I, IDIR, UT, UA, W, O
      DOUBLE PRECISION TNEW(MAXSTOKES), ACCUR,
     &                 TOLD(MAXSTOKES), TAZM(MAXSTOKES)

C  ###################
C  Fourier 0 component
C  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

C  Copy DIFFUSE Fourier component at all output angles and optical depths
C    If no SSCORR and no DBCORR, then two options apply:
C     (a) Convergence on STOKES = DIFFUSE + SSTRUNCATED + DBTRUNCATED
C              (full radiance, no SS correction, no DB correction)
C     (b) Convergence on STOKES = DIFFUSE alone (MS only mode)
C              (SSTRUNCATED + DBTRUNCATED, do not get calculated)

C  Full single scatter calculation is initialized to zero here (Version 2.3)

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

C    Add the single scatter component if flagged
C     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
C     Version 2.1.   Added outgoing correction flag to this.....
C     Version 2.3    Added Full single scatter flag

        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         DO IDIR = 1, N_DIRECTIONS
          W = WHICH_DIRECTIONS(IDIR)
          DO UT = 1, N_USER_LEVELS
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS (IB,I) + UA
                DO O = 1, NSTOKES
c         if (v.eq.1)write(45,*)UT,STOKES(UT,V,O,W),STOKES_SS(UT,V,1,W)
                  STOKES(UT,V,O,W) =
     &              STOKES(UT,V,O,W) + STOKES_SS(UT,V,O,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
         ENDDO
        ENDIF
c        pause

C   (New code, 6 May 2005)
C    Add the direct beam component if flagged (upwelling only)
C       Convergence on STOKES = STOKES + DBEXACT
C    Full single scatter option added Version 2.3

        IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
          DO UT = 1, N_USER_LEVELS
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS(IB,I) + UA
                DO O = 1, NSTOKES
                  STOKES(UT,V,O,UPIDX) =
     &              STOKES(UT,V,O,UPIDX) + STOKES_DB(UT,V,O)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

C  ######################
C  Fourier components > 0
C  ######################

      ELSE

C  No examination of convergence
C  -----------------------------

C  For Rayleigh atmosphere or if All Fourier components are required,
C     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY. OR. DO_ALL_FOURIER ) THEN

C  For each azimuth, add Fourier component

          DO UA = 1, LOCAL_N_USERAZM

C     - for direction, user optical depth, out stream

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

C  Examine convergence on intensity only 
C  -------------------------------------

C  convergence test applied to ALL directions AND
C                              ALL stream values (except near zenith) AND
C                              ALL azimuths taken together
C                              ALL user optical depths

        ELSE

C  Count number of occasions Fourier term addition is below accuracy level

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

C  set convergence counter TESTCONV

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

C  end convergence clause

        ENDIF

C  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER_COMPONENT .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IB) = FOURIER_COMPONENT
          ENDIF
        ENDIF

C  For all Fourier, keep saveing the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IB) = FOURIER_COMPONENT
        ENDIF

C  Finish iteration loop
   
      ENDIF

C  Finish

      RETURN
      END
