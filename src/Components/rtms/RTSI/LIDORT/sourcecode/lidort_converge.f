C ###########################################################
C #                                                         #
C #                    THE LIDORT FAMILY                    #
C #                                                         #
C #      (LInearized Discrete Ordinate Radiative Transfer)  #
C #       --         -        -        -         -          #
C #                                                         #
C ###########################################################

C ###########################################################
C #                                                         #
C #  Author :      Robert. J. D. Spurr                      #
C #                                                         #
C #  Address :     RT Solutions, Inc.                       #
C #                9 Channing Street                        #
C #                Cambridge, MA 02138, USA                 #
C #                                                         #
C #  Tel:          (617) 492 1183                           #
C #  Email :        rtsolutions@verizon.net                 #
C #                                                         #
C #  This Version :   3.3                                   #
C #  Release Date :   September 2007                        #
C #                                                         #
C #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
C #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
C #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
C #                                                         #
C ###########################################################

C    #####################################################
C    #                                                   #
C    #   This Version of LIDORT comes with a GNU-style   #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ###############################################################
C #                                                             #
C # Subroutines in this Module                                  #
C #                                                             #
C #            LIDORT_CONVERGE (master)                         #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_CONVERGE
     I      ( AZMFAC, FOURIER_COMPONENT, LOCAL_N_USERAZM,
     O        TESTCONV, IBEAM, LOCAL_ITERATION )

C  convergence testing on the Radiance intensity

C  Include files
C  -------------

      INCLUDE '../includes/LIDORT.PARS'

C  Include file of input variables
C  Include file of bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include files of already calculated results

      INCLUDE '../includes/LIDORT_SETUPS.VARS'
      INCLUDE '../includes/LIDORT_SINGSCAT.VARS'

C  Include file of final results

      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  input variables

      INTEGER          FOURIER_COMPONENT, LOCAL_N_USERAZM
      DOUBLE PRECISION AZMFAC
     &     (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)


C  modified/output variables

      LOGICAL          LOCAL_ITERATION
      INTEGER          TESTCONV, IBEAM

C  local variables
C  ---------------

C  local variables

      INTEGER          COUNT, COUNT_A
      INTEGER          I, IDIR, UT, UA, W, V
      DOUBLE PRECISION TNEW, ACCUR, TOLD, TAZM

C  ###################
C  Fourier 0 component
C  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

C  Copy DIFFUSE Fourier component at all output angles and optical depths
C    If no SSCORR and no DBCORR, then two options apply:
C     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
C              (full radiance, no SS correction, no DB correction)
C     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
C              (SSTRUNCATED + DBTRUNCATED do not get calculated)

        IF ( .not. DO_SSFULL ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_OUT_USERTAUS
              DO I = LOCAL_UM_START, N_OUT_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = UMOFF(IBEAM,I) + UA
                  INTENSITY(UT,V,W) = INTENSITY_F(UT,I,IBEAM,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_OUT_USERTAUS
              DO I = LOCAL_UM_START, N_OUT_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = UMOFF(IBEAM,I) + UA
                  INTENSITY(UT,V,W) = ZERO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C    Add the single scatter component if flagged
C     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
C     Version 3.2.   Added outgoing correction flag to this.....
C     Version 3.3    Added Full single scatter flag

        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
          DO IDIR = 1, N_DIRECTIONS
           W = WHICH_DIRECTIONS(IDIR)
           DO UT = 1, N_OUT_USERTAUS
             DO I = LOCAL_UM_START, N_OUT_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                 V = UMOFF(IBEAM,I) + UA
                 INTENSITY(UT,V,W) =
     &              INTENSITY(UT,V,W) + INTENSITY_SS(UT,V,W)
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ENDIF

C   (New code, 6 May 2005)
C    Add the  component if flagged (upwelling only)
C       Convergence on Radiance = Radiance_sofar + DBEXACT
C    Full single scatter option added Version 3.3

        IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
          DO UT = 1, N_OUT_USERTAUS
            DO I = LOCAL_UM_START, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V = UMOFF(IBEAM,I) + UA
                INTENSITY(UT,V,UPIDX) =
     &            INTENSITY(UT,V,UPIDX) + INTENSITY_DB(UT,V)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C------------------------------------------------- SECTION REMOVED 04/02/07
C  Multiple scatter source terms
C  -----------------------------
c        IF ( SAVE_LAYER_MSST ) THEN
C  downwelling output
C     Layer terms only, copy Fourier component 0
c          IF ( DO_DNWELLING ) THEN
c            DO N = 1, N_LAYERSOURCE_DN
c              DO I = LOCAL_UM_START, N_USER_STREAMS
c                DO UA = 1, LOCAL_N_USERAZM
c                  V = UMOFF(IBEAM,I) + UA
c                  INTENSITY_MSST_LAYER(N,V,DNIDX) = 
c     &                    INTENSITY_MSST_LAYER_F(N,I,IBEAM,DNIDX)
c                ENDDO
c              ENDDO
c            ENDDO
c          ENDIF
C  Upwelling output
c          IF ( DO_UPWELLING ) THEN
C  Direct beam output at BOA
C    Use exact value if flagged, otherwise copy Fourier m = 0
c            IF ( DO_DBCORRECTION ) THEN
c              DO I = LOCAL_UM_START, N_USER_STREAMS
c                DO UA = 1, LOCAL_N_USERAZM
c                  V = UMOFF(IBEAM,I) + UA
c                  INTENSITY_DBST_BOA(V) = EXACTDB_SOURCE(V)
c                ENDDO
c              ENDDO
c            ELSE
c              DO I = LOCAL_UM_START, N_USER_STREAMS
c                DO UA = 1, LOCAL_N_USERAZM
c                  V = UMOFF(IBEAM,I) + UA
c                  INTENSITY_DBST_BOA(V) =
c     &                  INTENSITY_DBST_BOA_F(I,IBEAM)
c                ENDDO
c              ENDDO
c            ENDIF
C  Diffuse source term at BOA (copy Fourier m = 0 )
c            DO I = LOCAL_UM_START, N_USER_STREAMS
c              DO UA = 1, LOCAL_N_USERAZM
c                V = UMOFF(IBEAM,I) + UA
c                INTENSITY_MSST_BOA(V) =
c     &                    INTENSITY_MSST_BOA_F(I,IBEAM)
c              ENDDO
c            ENDDO
C  Layer terms, copy Fourier component 0
c            DO I = LOCAL_UM_START, N_USER_STREAMS
c              DO UA = 1, LOCAL_N_USERAZM
c                V = UMOFF(IBEAM,I) + UA
c                DO N = NLAYERS, N_LAYERSOURCE_UP, -1
c                  INTENSITY_MSST_LAYER(N,V,UPIDX) = 
c     &                  INTENSITY_MSST_LAYER_F(N,I,IBEAM,UPIDX)
c                ENDDO
c              ENDDO
c            ENDDO
C  End upwelling
c          ENDIF
C  end save MS output
c        ENDIF
C------------------------------------------------- SECTION REMOVED 04/02/07

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
              DO UT = 1, N_OUT_USERTAUS
                DO I = LOCAL_UM_START, N_OUT_STREAMS
                  V = UMOFF(IBEAM,I) + UA
                  TOLD = INTENSITY(UT,V,W)
                  TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F(UT,I,IBEAM,W)
                  INTENSITY(UT,V,W) = TOLD + TAZM
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
              DO UT = 1, N_OUT_USERTAUS
                DO I = LOCAL_UM_START, N_OUT_STREAMS
                  V = UMOFF(IBEAM,I) + UA
                  TOLD = INTENSITY(UT,V,W)
                  TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F(UT,I,IBEAM,W)
                  TNEW = TOLD + TAZM
                  IF ( TAZM .NE. ZERO ) THEN
                    ACCUR     = DABS(TAZM/TNEW)
                    IF ( ACCUR .LT. LIDORT_ACCURACY ) THEN
                      COUNT   = COUNT + 1
                      COUNT_A = COUNT_A + 1
                    ENDIF
                  ELSE
                    COUNT   = COUNT + 1
                    COUNT_A = COUNT_A + 1
                  ENDIF
                  INTENSITY(UT,V,W)     = TNEW
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
              FOURIER_SAVED(IBEAM) = FOURIER_COMPONENT
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IBEAM) = 2*NSTREAMS - 1
          ENDIF

C  end convergence clause

        ENDIF

C  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER_COMPONENT .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IBEAM) = FOURIER_COMPONENT
          ENDIF
        ENDIF

C  For all Fourier, keep saveing the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IBEAM) = FOURIER_COMPONENT
        ENDIF

C------------------------------------------------- SECTION REMOVED 04/02/07
C  Update series for Multiple scatter source terms (if flagged)
C  ------------------------------------------------------------
c        IF ( SAVE_LAYER_MSST ) THEN
c         DO UA = 1, LOCAL_N_USERAZM
C  (Downwelling) - just the layer source terms
c          IF ( DO_DNWELLING ) THEN
c           W = DNIDX
c           DO N = 1, N_LAYERSOURCE_DN
c             DO I = LOCAL_UM_START, N_USER_STREAMS
c               V = UMOFF(IBEAM,I) + UA
c               TOLD = INTENSITY_MSST_LAYER(N,V,W)
c               TAZM= AZMFAC(UA)*INTENSITY_MSST_LAYER_F(N,I,IBEAM,W)
c               INTENSITY_MSST_LAYER(N,V,W) = TOLD + TAZM
c             ENDDO
c           ENDDO
c          ENDIF
C  (Upwelling)
c          IF ( DO_UPWELLING ) THEN
c           W = UPIDX
C  ..... Upwelling, do BOA source terms
C        (contributions are zero if no surface reflection)
C        (only do direct beam if no correction)
c           DO I = LOCAL_UM_START, N_USER_STREAMS
c            V = UMOFF(IBEAM,I) + UA
c            TOLD = INTENSITY_MSST_BOA(V)
c            TAZM = AZMFAC(UA)*INTENSITY_MSST_BOA_F(I,IBEAM)
c            INTENSITY_MSST_BOA(V) = TOLD + TAZM
c            IF ( .NOT.DO_DBCORRECTION ) THEN
c              TOLD = INTENSITY_DBST_BOA(V)
c              TAZM = AZMFAC(UA)*INTENSITY_DBST_BOA_F(I,IBEAM)
c              INTENSITY_DBST_BOA(V) = TOLD + TAZM
c            ENDIF
c           ENDDO
C  ..... Upwelling, do layer source terms
c           DO I = LOCAL_UM_START, N_USER_STREAMS
c            V = UMOFF(IBEAM,I) + UA
c            DO N = NLAYERS, N_LAYERSOURCE_UP, -1
c              TOLD = INTENSITY_MSST_LAYER(N,V,W)
c              TAZM = AZMFAC(UA)*INTENSITY_MSST_LAYER_F(N,I,IBEAM,W)
c              INTENSITY_MSST_LAYER(N,V,W) = TOLD + TAZM
c            ENDDO
c           ENDDO
c          ENDIF
C  End azimuth loop
c         ENDDO
C  End save layer clause
c        ENDIF
C------------------------------------------------- SECTION REMOVED 04/02/07

C  Finish iteration loop
   
      ENDIF

C  Finish

      RETURN
      END
