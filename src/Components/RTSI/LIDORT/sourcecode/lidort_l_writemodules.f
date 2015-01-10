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
C #              LIDORT_L_WRITEFOURIER                          #
C #              LIDORT_L_WRITERESULTS                          #
C #              LIDORT_L_WRITEINPUT                            #
C #              LIDORT_L_WRITESCEN                             #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_L_WRITEFOURIER
     &     ( DO_INCLUDE_SURFACE, FUNIT, FOURIER_COMPONENT )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of linearization input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include file with results for writing

      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  input variables
C  ---------------

      INTEGER          FOURIER_COMPONENT, FUNIT
      LOGICAL          DO_INCLUDE_SURFACE
      
C  local variables
C  ---------------

      INTEGER          I, IDIR, UT, WDIR, Z, Q, K, IB

C  Profile Atmospheric Weighting function output
C  #############################################

      IF ( DO_PROFILE_LINEARIZATION ) THEN

C  overall header

        WRITE(FUNIT,'(a)')' '
        write(FUNIT,'(/a,i3/a/)')
     & 'Atmospheric Profile Weighting function: Fourier component',
     &         FOURIER_COMPONENT,
     & '----------------------------------------------------------'

        WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
        WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS
        WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of solar zenith angles = ',NBEAMS

C  detailed output

        DO IDIR = 1, N_DIRECTIONS
          WDIR = WHICH_DIRECTIONS(IDIR)

C  linearization control

          DO K = 1, NLAYERS
           IF ( LAYER_VARY_FLAG(K) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(K)

C  direction headers

             IF (WDIR .EQ. UPIDX ) THEN
              WRITE(FUNIT,'(/A/A,I2/A,A31)')
     &     '  --> Upwelling Atmospheric Weighting functions',
     &     '        * for variations in layer = ',K,
     &     '        * with respect to : ',PROFILEWF_NAMES(Q)
             ELSE IF (WDIR .EQ. DNIDX ) THEN
              WRITE(FUNIT,'(/A/A,I2/A,A31)')
     &     '  --> Downwelling Atmospheric Weighting functions',
     &     '        * for variations in layer = ',K,
     &     '        * with respect to : ',PROFILEWF_NAMES(Q)
             ENDIF

C  output loops

             DO IB = 1, NBEAMS
              WRITE(FUNIT,'(/a,2x,f10.4/)')'SZA (degs)',BEAM_SZAS(IB)
              write(FUNIT,'(a/)')
     &             'output angles | Optical depth values-----> '
              WRITE(FUNIT,'(14X,50(1pe15.5))')
     *          ( USER_TAUS_INPUT(UT),UT = 1, N_OUT_USERTAUS)
              DO I = 1, N_OUT_STREAMS
               WRITE(FUNIT,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
     &           (ATMOSWF_F(Q,K,UT,I,IB,WDIR),UT = 1, N_OUT_USERTAUS)
              ENDDO
             ENDDO

C  close control and loops

            ENDDO
           ENDIF
          ENDDO
        ENDDO

C  End atmospheric weighting function stuff

      ENDIF

C  Column Atmospheric Weighting function output
C  ############################################

      IF ( DO_COLUMN_LINEARIZATION ) THEN

C  overall header

        WRITE(FUNIT,'(a)')' '
        write(FUNIT,'(/a,i3/a/)')
     & 'Atmospheric Column Weighting function: Fourier component',
     &         FOURIER_COMPONENT,
     & '----------------------------------------------------------'

        WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
        WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS
        WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of solar zenith angles = ',NBEAMS

C  detailed output

        DO IDIR = 1, N_DIRECTIONS
          WDIR = WHICH_DIRECTIONS(IDIR)
          DO Q = 1, N_TOTALCOLUMN_WFS

C  direction headers

            IF (WDIR .EQ. UPIDX ) THEN
              WRITE(FUNIT,'(/A/A,A31)')
     &     '  --> Upwelling Column Atmospheric Weighting functions',
     &     '        * with respect to : ',COLUMNWF_NAMES(Q)
            ELSE IF (WDIR .EQ. DNIDX ) THEN
              WRITE(FUNIT,'(/A/A,A31)')
     &     '  --> Downwelling Column Atmospheric Weighting functions',
     &     '        * with respect to : ',COLUMNWF_NAMES(Q)
            ENDIF

C  output loops

            DO IB = 1, NBEAMS
              WRITE(FUNIT,'(/a,2x,f10.4/)')'SZA (degs)',BEAM_SZAS(IB)
              write(FUNIT,'(a/)')
     &             'output angles | Optical depth values-----> '
              WRITE(FUNIT,'(14X,50(1pe15.5))')
     *          ( USER_TAUS_INPUT(UT),UT = 1, N_OUT_USERTAUS)
              DO I = 1, N_OUT_STREAMS
               WRITE(FUNIT,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
     &           (ATMOSWF_F(Q,0,UT,I,IB,WDIR),UT = 1, N_OUT_USERTAUS)
              ENDDO
            ENDDO

C  close control and loops

          ENDDO
        ENDDO

C  End atmospheric bulk/column weighting function stuff

      ENDIF

C  Surface Weighting function output
C  #################################

      IF ( DO_SURFACE_LINEARIZATION.AND.DO_INCLUDE_SURFACE ) THEN

C  only required for non-Lambertian cases

        IF ( DO_LAMBERTIAN_SURFACE ) RETURN

C  overall header

        WRITE(FUNIT,'(a)')' '
        write(FUNIT,'(/a,i3/a/)')
     & 'Surface Weighting function output for Fourier component',
     &         FOURIER_COMPONENT,
     & '----------------------------------------------------------'

        WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of optical depths  = ',N_OUT_USERTAUS
        WRITE(FUNIT,FMT_INTEGER) 
     &      'Total number of output angles   = ',N_OUT_STREAMS
        WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of Surface wt.fns. = ',N_TOTALBRDF_WFS

C  detailed output

        DO IDIR = 1, N_DIRECTIONS
          WDIR = WHICH_DIRECTIONS(IDIR)

C  direction headers

          IF (WDIR .EQ. UPIDX ) THEN
            WRITE(FUNIT,'(/A)')
     &     '  --> Upwelling Surface Weighting functions (Fourier)'
          ELSE IF (WDIR .EQ. DNIDX ) THEN
            WRITE(FUNIT,'(/A)')
     &     '  --> Downwelling Surface Weighting functions (Fourier)'
          ENDIF

C  for each surface weighting function

          DO Z = 1, N_TOTALBRDF_WFS

            WRITE(FUNIT,'(/A,I2)')
     &     '  -->  Surface Weighting functions, #: ',Z

C  output loops

            DO IB = 1, NBEAMS
              WRITE(FUNIT,'(/a,2x,f10.4/)')'SZA (degs)',BEAM_SZAS(IB)
              write(FUNIT,'(a/)')
     &             'output angles | Optical depth values-----> '
              WRITE(FUNIT,'(14X,50(1pe15.5))')
     *          ( USER_TAUS_INPUT(UT),UT = 1, N_OUT_USERTAUS)
              DO I = 1, N_OUT_STREAMS
               WRITE(FUNIT,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
     &           (SURFACEWF_F(Z,UT,I,IB,WDIR),UT = 1, N_OUT_USERTAUS)
              ENDDO
            ENDDO

C  end loops

          ENDDO
        ENDDO

C  End surface weighting function stuff

      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_WRITERESULTS (RUN,LOCAL_NUSERAZMS)

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of linearization input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  include files with results for writing

      INCLUDE '../includes/LIDORT_RESULTS.VARS'
      INCLUDE '../includes/LIDORT_L_RESULTS.VARS'

C  input variables
C  ---------------

      INTEGER      RUN
      INTEGER      LOCAL_NUSERAZMS

C  local variables
C  ---------------

      INTEGER      I, UA, IDIR, UT, WDIR, Z
      INTEGER      Q, K, IB, V, FMAX, F
      INTEGER      VINDEX(MAXBEAMS,MAX_OUT_STREAMS,MAX_USER_RELAZMS)

C  Beam Attenuation summary
C  ========================

      write(RUN,'(/a/a/)')'Results Output summary',
     &                     '----------------------'

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        write(RUN,'(a)')
     &       'Classical solution of beam particular integral'
      ELSE
        write(RUN,'(a)')
     &       'Green function solution of beam particular integral'
      ENDIF

      IF ( DO_PLANE_PARALLEL ) THEN
        write(RUN,'(a)')
     &       'Plane parallel approximation to beam attenuation'
      ELSE
        write(RUN,'(a)')
     &       'Pseudo-spherical approximation to beam attenuation'
        write(RUN,'(a)')
     &      ' * Average secant approximation was used'
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          write(RUN,'(a)')
     &      ' * Refractive geometry was used (Snells Law)'
        ELSE
          write(RUN,'(a)')
     &      ' * Straight line geometry was used (no refraction)'
        ENDIF
      ENDIF

      IF ( DO_FULLRAD_MODE ) THEN
        write(RUN,'(a)')
     &       'Full radiance calculation has been performed'
        IF ( DO_SSCORR_NADIR ) THEN
          write(RUN,'(a)')
     &  '  --> Nakajima-Tanaka TMS single scatter correction (nadir)'
        ENDIF
        IF ( DO_SSCORR_OUTGOING ) THEN
          write(RUN,'(a)')
     &  '  --> Nakajima-Tanaka TMS single scatter correction (outgoing)'
        ENDIF
        IF ( .NOT.DO_SSCORR_NADIR .AND. .NOT.DO_SSCORR_OUTGOING ) THEN
          write(RUN,'(a)')
     &  '  --> No single scatter correction has been applied'
        ENDIF
      ELSE
        write(RUN,'(a)')
     &  'ONLY Multiple-scatter radiance calculation has been performed'
c        IF ( SAVE_LAYER_MSST ) THEN
c          write(RUN,'(a)')
c     &      '--> Layer multiple-scatter source terms were output'
c        ENDIF
      ENDIF

      IF ( .NOT.DO_NO_AZIMUTH ) THEN
        IF ( DO_DOUBLE_CONVTEST ) THEN
          write(RUN,'(/a)')
     &   'Double convergence test was used for Fourier Azimuth Series'
        ELSE
          write(RUN,'(/a)')
     &   'Single convergence test was used for Fourier Azimuth Series'
        ENDIF
        write(RUN,'(a,F10.7/)')
     &   ' --> Accuracy level was pre-set at : ',LIDORT_ACCURACY

        write(RUN,'(a,I5)')' - Number Solar zenith angles : ',NBEAMS
        write(run,'(/a/)')
     &         '   SZA  | # Fourier | Fourier breakdown -->'
        fmax = -100
        DO ib = 1, nbeams
          fmax = MAX(fmax,fourier_saved(ib))
        END DO
        DO ib = 1, nbeams
           WRITE(run,'(1x,f7.2,4x,i3,6x,50(L1))')
     &            beam_szas(ib),fourier_saved(ib),
     &            (do_multibeam(ib,f),f=0,fmax)
        END DO
      ELSE
        write(RUN,'(/a)')
     &   'Azimuth independent output only (Fourier = 0)'
      ENDIF

C  Atmospheric Profile Weighting function output
C  #############################################

      IF ( DO_PROFILE_LINEARIZATION ) THEN

C  control point for avoiding weighting function output
         
        IF ( DO_MVOUT_ONLY ) GO TO 401

C  local number of azimuths

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_NUSERAZMS = 1
        ELSE
          LOCAL_NUSERAZMS = N_USER_RELAZMS
        ENDIF

C  Indexing
C   This loop moved to its current localtion. 27 November 2006.

        DO IB = 1, NBEAMS
         DO I = 1, N_OUT_STREAMS
          DO UA = 1, LOCAL_NUSERAZMS
           VINDEX(IB,I,UA) = UMOFF(IB,I) + UA
          ENDDO
         ENDDO
        ENDDO

C  overall header

        write(RUN,'(/a/a/)')
     $     'Atmospheric Profile Weighting function output',
     &     '---------------------------------------------'

C  start beam loop

        DO IB = 1, NBEAMS
         write(RUN,FMT_REAL)
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',BEAM_SZAS(IB)
         write(RUN,'(a)')' '

C  start azimuth loop

         DO UA = 1, LOCAL_NUSERAZMS

C  azimuth angle header

          IF ( DO_NO_AZIMUTH ) THEN
            write(RUN,FMT_CHAR)
     *       '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
          ELSE
            write(RUN,FMT_REAL)
     * '** RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
          ENDIF

          WRITE(RUN,'(a)')' '
          WRITE(RUN,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
          WRITE(RUN,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS

C  detailed output

          DO IDIR = 1, N_DIRECTIONS
           WDIR = WHICH_DIRECTIONS(IDIR)

C  linearization control

           DO K = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(K) ) THEN
             DO Q = 1, LAYER_VARY_NUMBER(K)

C  direction headers

              IF (WDIR .EQ. UPIDX ) THEN
               WRITE(RUN,'(/A/A,I2/A,A31)')
     &     '  --> Upwelling Atmospheric Weighting functions',
     &     '        * for variations in layer = ',K,
     &     '        * with respect to : ',PROFILEWF_NAMES(Q)
              ELSE IF (WDIR .EQ. DNIDX ) THEN
               WRITE(RUN,'(/A/A,I2/A,A31)')
     &     '  --> Downwelling Atmospheric Weighting functions',
     &     '        * for variations in layer = ',K,
     &     '        * with respect to : ',PROFILEWF_NAMES(Q)
              ENDIF

C  output loop

              write(RUN,'(a)') 'output angles | optical depths ---> '
              WRITE(RUN,'(14x,50(2x,1pe11.4,2x))')
     &           (USER_TAUS_INPUT(UT), UT = 1, N_OUT_USERTAUS)
              DO I = 1, N_OUT_STREAMS
               V = UMOFF(IB,I) + UA
               WRITE(RUN,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
     &         (PROFILEWF(Q,K,UT,V,WDIR), UT = 1, N_OUT_USERTAUS)
              ENDDO

C  close control and loops

             ENDDO
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO

C  ------------------------------------------ REMOVED, 02 April 2007
C  Mult-scat source term output
C  ----------------------------
c        IF ( .NOT. SAVE_LAYER_MSST ) GO TO 401
C  overall header
c        write(RUN,'(/a/a/)')
c     &  'Mult-scat source term atmospheric weighting function output',
c     &  '-----------------------------------------------------------'
C  start azimuth, SZA loops
c        DO IB = 1, NBEAMS
c         write(RUN,FMT_REAL)
c     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',BEAM_SZAS(IB)
c         write(RUN,'(a)')' '
c         DO UA = 1, LOCAL_NUSERAZMS
C  azimuth angle header
c          IF ( DO_NO_AZIMUTH ) THEN
c            write(RUN,FMT_CHAR)
c     *       '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
c          ELSE
c            write(RUN,FMT_REAL)
c     * '* * RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
c          ENDIF
C  User angle header
c          DO I = 1, N_USER_STREAMS
c           V = UMOFF(IB,I) + UA
c           write(RUN,FMT_REAL)
c     * '* * RESULTS FOR USER ANGLE (degs)=',USER_ANGLES_INPUT(I)
C  detailed output
c           DO IDIR = 1, N_DIRECTIONS
c            WDIR = WHICH_DIRECTIONS(IDIR)
c            IF (WDIR .EQ. UPIDX ) THEN
c             DO K = 1, NLAYERS
c              IF ( LAYER_VARY_FLAG(K) ) THEN
c               DO Q = 1, LAYER_VARY_NUMBER(K)
c                WRITE(RUN,'(/A/A,I2/A,I2)')
c     &   '  --> Upwelling Mult-scat Source Term Weighting functions',
c     &     '        * W.R.T parameters varying in layer = ',K,
c     &     '        * Parameter number = ',Q
c                DO N = N_LAYERSOURCE_UP, NLAYERS
c                 WRITE(RUN,'(1X,I3,3X,10(1PE15.5))') N,
c     &              ATMOSWF_MSCATSTERM(Q,K,N,V,WDIR)
c                ENDDO
c                WRITE(RUN,'(/A)')
c     &     '  --> Upwelling BOA multiple source terms'
c                WRITE(RUN,'(A7,10(1PE15.5))')'diffuse',
c     &               ATMOSWF_MSCATBOA_STERM(Q,K,V)
c                WRITE(RUN,'(A7,10(1PE15.5))')'direct ',
c     &               ATMOSWF_DIRECTBOA_STERM(Q,K,V)
c               ENDDO
c              ENDIF
c             ENDDO
c            ELSE IF (WDIR .EQ. DNIDX ) THEN
c             DO K = 1, NLAYERS
c              IF ( LAYER_VARY_FLAG(K) ) THEN
c               DO Q = 1, LAYER_VARY_NUMBER(K)
c                WRITE(RUN,'(/A/A,I2/A,I2)')
c     &   '  --> Downwelling Mult-scat Source Term Weighting functions',
c     &     '        * W.R.T parameters varying in layer = ',K,
c     &     '        * Parameter number = ',Q
c                DO N = 1, N_LAYERSOURCE_DN
c                 WRITE(RUN,'(1X,I3,3X,10(1PE15.5))') N,
c     &              ATMOSWF_MSCATSTERM(Q,K,N,V,WDIR)
c                ENDDO
c               ENDDO
c              ENDIF
c             ENDDO
c            ENDIF
C  Finish loops
c           ENDDO
c          ENDDO
c         ENDDO
c        ENDDO
C  ------------------------------------------ REMOVED, 02 April 2007

C  Mean-value output
C  -----------------

401     CONTINUE
        IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

          write(RUN,'(/a/a)')
     &          'Mean value atmospheric weighting function output',
     &          '------------------------------------------------'

C  start beam loop

          DO IB = 1, NBEAMS

           write(RUN,FMT_REAL)
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',BEAM_SZAS(IB)
           write(RUN,'(a)')' '

C  detailed output

           DO IDIR = 1, N_DIRECTIONS
            WDIR = WHICH_DIRECTIONS(IDIR)

C  linearization control

            DO K = 1, NLAYERS
             IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)

C  direction headers

               IF (WDIR .EQ. UPIDX ) THEN
                WRITE(RUN,'(/A/A,I2/A,I2/)')
     & '  --> Upwelling Mean intensity and Flux Weighting functions',
     &     '        * W.R.T parameters varying in layer = ',K,
     &     '        * Parameter number = ',Q
               ELSE IF (WDIR .EQ. DNIDX ) THEN
                WRITE(RUN,'(/A/A,I2/A,I2/)')
     & '  --> Downwelling Mean intensity and Flux Weighting functions',
     &     '        * W.R.T parameters varying in layer = ',K,
     &     '        * Parameter number = ',Q
               ENDIF

C  optical depth header

               WRITE(RUN,'(a/)')
     &      'optical depth     mean int. WF      flux WF'

C  optical depth loop

               DO UT = 1, N_OUT_USERTAUS
                WRITE(RUN,'(2x,F9.5,3x,2E17.7)')
     &                    USER_TAUS_INPUT(UT),
     &                    MINT_ATMOSWF(Q,K,UT,IB,WDIR),
     &                    FLUX_ATMOSWF(Q,K,UT,IB,WDIR)
               ENDDO

C  close control and loops

              ENDDO
             ENDIF
            ENDDO
           ENDDO
          ENDDO
        ENDIF

C  End atmospheric pprofile weighting function stuff

      ENDIF

C  Atmospheric Column Weighting function output
C  ############################################

      IF ( DO_COLUMN_LINEARIZATION ) THEN

C  control point for avoiding weighting function output
         
        IF ( DO_MVOUT_ONLY ) GO TO 4401

C  local number of azimuths

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_NUSERAZMS = 1
        ELSE
          LOCAL_NUSERAZMS = N_USER_RELAZMS
        ENDIF

C  overall header

        write(RUN,'(/a/a/)')
     $     'Atmospheric Column Weighting function output',
     &     '--------------------------------------------'

C  start beam loop

        DO IB = 1, NBEAMS
         write(RUN,FMT_REAL)
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',BEAM_SZAS(IB)
         write(RUN,'(a)')' '

C  start azimuth loop

         DO UA = 1, LOCAL_NUSERAZMS

C  azimuth angle header

          IF ( DO_NO_AZIMUTH ) THEN
            write(RUN,FMT_CHAR)
     *       '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
          ELSE
            write(RUN,FMT_REAL)
     * '** RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
          ENDIF

          WRITE(RUN,'(a)')' '
          WRITE(RUN,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
          WRITE(RUN,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS

C  detailed output

          DO IDIR = 1, N_DIRECTIONS
            WDIR = WHICH_DIRECTIONS(IDIR)
 
C  linearization control

            DO Q = 1, N_TOTALCOLUMN_WFS

C  direction headers

              IF (WDIR .EQ. UPIDX ) THEN
               WRITE(RUN,'(/A/A,A31)')
     &     '  --> Upwelling Atmospheric Column Weighting functions',
     &     '        * with respect to : ',COLUMNWF_NAMES(Q)
              ELSE IF (WDIR .EQ. DNIDX ) THEN
               WRITE(RUN,'(/A/A,A31)')
     &     '  --> Downwelling Atmospheric Column Weighting functions',
     &     '        * with respect to : ',COLUMNWF_NAMES(Q)
              ENDIF

C  output loop

              write(RUN,'(a)') 'output angles | optical depths ---> '
              WRITE(RUN,'(14x,50(2x,1pe11.4,2x))')
     &           (USER_TAUS_INPUT(UT), UT = 1, N_OUT_USERTAUS)
              DO I = 1, N_OUT_STREAMS
               V = UMOFF(IB,I) + UA
               WRITE(RUN,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
     &         (COLUMNWF(Q,UT,V,WDIR), UT = 1, N_OUT_USERTAUS)
              ENDDO

C  close control and loops

             ENDDO
           ENDDO
         ENDDO
        ENDDO

C Commented out-----------------------------------------------------
C  Mult-scat source term output
C  ----------------------------
c
c        IF ( .NOT. SAVE_LAYER_MSST ) GO TO 401
C  overall header
c        write(RUN,'(/a/a/)')
c     &  'Mult-scat source term column weighting function output',
c     &  '------------------------------------------------------'
C  start azimuth, SZA loops
c        DO IB = 1, NBEAMS
c         write(RUN,FMT_REAL)
c     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',BEAM_SZAS(IB)
c         write(RUN,'(a)')' '
c         DO UA = 1, LOCAL_NUSERAZMS
C  azimuth angle header
c          IF ( DO_NO_AZIMUTH ) THEN
c            write(RUN,FMT_CHAR)
c     *       '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
c          ELSE
c            write(RUN,FMT_REAL)
c     * '* * RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
c          ENDIF
C  User angle header
c          DO I = 1, N_USER_STREAMS
c           V = UMOFF(IB,I) + UA
c           write(RUN,FMT_REAL)
c     * '* * RESULTS FOR USER ANGLE (degs)=',USER_ANGLES_INPUT(I)
C  detailed output
c           DO IDIR = 1, N_DIRECTIONS
c             WDIR = WHICH_DIRECTIONS(IDIR)
c             IF (WDIR .EQ. UPIDX ) THEN
c               DO Q = 1, N_TOTALCOLUMN_WFS
c                WRITE(RUN,'(/A/A,I2)')
c     &   '  --> Upwelling MS sourceterm column Weighting functions',
c     &     '        * Parameter number = ',Q
c                DO N = N_LAYERSOURCE_UP, NLAYERS
c                 WRITE(RUN,'(1X,I3,3X,10(1PE15.5))') N,
c     &              COLUMNWF_MSCATSTERM(Q,N,V,WDIR)
c                ENDDO
c                WRITE(RUN,'(/A)')
c     &     '  --> Upwelling BOA multiple source terms'
c                WRITE(RUN,'(A7,10(1PE15.5))')'diffuse',
c     &               COLUMNWF_MSCATBOA_STERM(Q,V)
c                WRITE(RUN,'(A7,10(1PE15.5))')'direct ',
c     &               COLUMNWF_DIRECTBOA_STERM(Q,V)
c               ENDDO
c             ELSE IF (WDIR .EQ. DNIDX ) THEN
c               DO Q = 1, N_TOTALCOLUMN_WFS
c                WRITE(RUN,'(/A/A,I2)')
c     &   '  --> Downwelling MS sourceterm column Weighting functions',
c     &     '        * Parameter number = ',Q
c                DO N = 1, N_LAYERSOURCE_DN
c                 WRITE(RUN,'(1X,I3,3X,10(1PE15.5))') N,
c     &              COLUMNWF_MSCATSTERM(Q,N,V,WDIR)
c                ENDDO
c               ENDDO
c             ENDIF
C  Finish loops
c           ENDDO
c          ENDDO
c         ENDDO
c        ENDDO
C---------------------------------------Commented out

C  Mean-value output
C  -----------------

4401     CONTINUE
        IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

          write(RUN,'(/a/a)')
     &  'Mean value atmospheric Column weighting function output',
     &  '-------------------------------------------------------'

C  start beam loop

          DO IB = 1, NBEAMS

           write(RUN,FMT_REAL)
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',BEAM_SZAS(IB)
           write(RUN,'(a)')' '

C  detailed output

           DO IDIR = 1, N_DIRECTIONS
             WDIR = WHICH_DIRECTIONS(IDIR)

C  linearization control

             DO Q = 1, N_TOTALCOLUMN_WFS

C  direction headers

               IF (WDIR .EQ. UPIDX ) THEN
                WRITE(RUN,'(/A/A,I2/)')
     & '  --> Upwelling Mean-I & Flux: Column Weighting functions',
     &     '        * Parameter number = ',Q
               ELSE IF (WDIR .EQ. DNIDX ) THEN
                WRITE(RUN,'(/A/A,I2/)')
     & '  --> Downwelling Mean_I & Flux: Column Weighting functions',
     &     '        * Parameter number = ',Q
               ENDIF

C  optical depth header

               WRITE(RUN,'(a/)')
     &      'optical depth     mean int. WF      flux WF'

C  optical depth loop

               DO UT = 1, N_OUT_USERTAUS
                WRITE(RUN,'(2x,F9.5,3x,2E17.7)')
     &                    USER_TAUS_INPUT(UT),
     &                    MINT_ATMOSWF(Q,0,UT,IB,WDIR),
     &                    FLUX_ATMOSWF(Q,0,UT,IB,WDIR)
               ENDDO

C  close control and loops

             ENDDO
           ENDDO
          ENDDO
        ENDIF

C  End atmospheric column weighting function stuff

      ENDIF

C  Surface Weighting function output
C  #################################

      IF ( DO_SURFACE_LINEARIZATION ) THEN

C  control point for avoiding weighting function output
         
        IF ( DO_MVOUT_ONLY ) GO TO 402

C  local number of azimuths

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_NUSERAZMS = 1
        ELSE
          IF ( DO_LAMBERTIAN_SURFACE ) THEN
            LOCAL_NUSERAZMS = 1
          ELSE
            LOCAL_NUSERAZMS = N_USER_RELAZMS
          ENDIF
        ENDIF

C  overall header

        write(RUN,'(/a/a/)')'Surface Weighting function output',
     &                      '---------------------------------'

C  start beam loop

        DO IB = 1, NBEAMS

         write(RUN,FMT_REAL)
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',BEAM_SZAS(IB)
         write(RUN,'(a)')' '

C  start azimuth loop

         DO UA = 1, LOCAL_NUSERAZMS

C  azimuth angle header

          IF ( DO_NO_AZIMUTH ) THEN
           write(RUN,FMT_CHAR)
     *       '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
          ELSE
           IF ( DO_LAMBERTIAN_SURFACE ) THEN
            write(RUN,FMT_CHAR)
     *   '** Lambertian Surface --> Results AZIMUTH-INDEPENDENT **'
           ELSE
            write(RUN,FMT_REAL)
     * '** RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
           ENDIF
          ENDIF

          WRITE(RUN,'(a)')' '
          WRITE(RUN,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
          WRITE(RUN,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS

C  detailed output

          DO IDIR = 1, N_DIRECTIONS
           WDIR = WHICH_DIRECTIONS(IDIR)

C  direction headers

           IF (WDIR .EQ. UPIDX ) THEN
            WRITE(RUN,'(/A)')
     &     '  --> Upwelling Surface Weighting functions'
           ELSE IF (WDIR .EQ. DNIDX ) THEN
            WRITE(RUN,'(/A)')
     &     '  --> Downwelling Surface Weighting functions'
           ENDIF

C  for each Surface weighting function

           DO Z = 1, N_TOTALBRDF_WFS

            WRITE(RUN,'(/A,I2)')
     &     '  -->  Surface Weighting functions, #: ',Z

C  output loop

            write(RUN,'(a)') 'output angles | optical depths ---> '
            WRITE(RUN,'(14x,50(2x,1pe11.4,2x))')
     &         (USER_TAUS_INPUT(UT), UT = 1, N_OUT_USERTAUS)
            DO I = 1, N_OUT_STREAMS
             V = UMOFF(IB,I) + UA
             WRITE(RUN,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
     &         (SURFACEWF(Z,UT,V,WDIR), UT = 1, N_OUT_USERTAUS)
            ENDDO

C  close control and loops

           ENDDO
          ENDDO
         ENDDO
        ENDDO

C  ------------------------------------------ REMOVED, 02 April 2007
C  Mult-scat source term output
C  ----------------------------
C        IF ( .NOT. SAVE_LAYER_MSST ) GO TO 402
C  overall header
c        write(RUN,'(/a/a/)')
c     &     'Mult-scat source term surface weighting function output',
c     &     '-------------------------------------------------------'
C  start beam loop
c        DO IB = 1, NBEAMS
c         write(RUN,FMT_REAL)
c     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',BEAM_SZAS(IB)
c         write(RUN,'(a)')' '
C  start azimuth loop
c         DO UA = 1, LOCAL_NUSERAZMS
C  azimuth angle header
c          IF ( DO_NO_AZIMUTH ) THEN
c           write(RUN,FMT_CHAR)
c     *       '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
c          ELSE
c           write(RUN,FMT_REAL)
c     * '** RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
c          ENDIF
C  detailed output
c          DO IDIR = 1, N_DIRECTIONS
c           WDIR = WHICH_DIRECTIONS(IDIR)
c           IF (WDIR .EQ. UPIDX ) THEN
c            DO Z = 1, N_TOTALBRDF_WFS
c             WRITE(RUN,'(/A/A,I2)')
c     &     '  --> Upwelling Mult-scat Source Term Weighting functions',
c     &     '  -->   For BRDF Kernel: ',Z
c             WRITE(RUN,'(/A,10(2x,F10.5,1X))')
c     &     'angle->', (USER_ANGLES_INPUT(I),I=1,N_USER_STREAMS)
c             WRITE(RUN,'(A)')' layer'      
c             DO N = N_LAYERSOURCE_UP, NLAYERS
c              WRITE(RUN,'(1X,I3,3X,10(1PE13.5))')
c     &         N,(SURFACEWF_MSCATSTERM(Z,N,VINDEX(IB,I,UA),WDIR),
c     &               I=1,N_USER_STREAMS)
c             ENDDO
c             WRITE(RUN,'(/A/)')
c     & '  --> Upwelling Mult-scat BOA Source Term Weighting function'
c             WRITE(RUN,'(7x,10(1PE13.5))')
c     &           (SURFACEWF_MSCATBOA_STERM(Z,VINDEX(IB,I,UA)),
c     &             I=1,N_USER_STREAMS)
c            ENDDO
c           ELSE IF (WDIR .EQ. DNIDX ) THEN
c            DO Z = 1, N_TOTALBRDF_WFS
c             WRITE(RUN,'(/A/A,I2)')
c     &   '  --> Downwelling Mult-scat Source Term Weighting functions',
c     &   '  -->   For BRDF Kernel: ',Z
c             WRITE(RUN,'(/A,10(2x,F10.5,1X))')
c     &     'angle->', (USER_ANGLES_INPUT(I),I=1,N_USER_STREAMS)
c             WRITE(RUN,'(A)')' layer'      
c             DO N = 1, N_LAYERSOURCE_DN
c              WRITE(RUN,'(1X,I3,3X,10(1PE13.5))')
c     &        N,(SURFACEWF_MSCATSTERM(Z,N,VINDEX(IB,I,UA),WDIR),
c     &                I=1,N_USER_STREAMS)
c             ENDDO
c            ENDDO
c           ENDIF
c          ENDDO
c         ENDDO
c        ENDDO
C  ------------------------------------------ REMOVED, 02 April 2007

C  Mean-value output
C  -----------------

402     CONTINUE
        IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

         write(RUN,'(/a/a)')
     &          'Mean value surface weighting function output',
     &          '--------------------------------------------'

C  start beam loop

         DO IB = 1, NBEAMS

         write(RUN,FMT_REAL)
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',BEAM_SZAS(IB)
         write(RUN,'(a)')' '

C  detailed output

          DO IDIR = 1, N_DIRECTIONS
            WDIR = WHICH_DIRECTIONS(IDIR)

C  direction headers

            IF (WDIR .EQ. UPIDX ) THEN
              WRITE(RUN,'(/A/A/)')
     & '  --> Upwelling Mean intensity and Flux Weighting functions',
     &     '      * W.R.T surface variation '
            ELSE IF (WDIR .EQ. DNIDX ) THEN
              WRITE(RUN,'(/A/A/)')
     & '  --> Downwelling Mean intensity and Flux Weighting functions',
     &     '      * W.R.T surface variation '
            ENDIF

C  Loop over kernels

            DO Z = 1, N_TOTALBRDF_WFS

C  optical depth header

              WRITE(RUN,'(a,I2/a/)')
     &      '  --> For Surface WF number :',Z,
     &      'optical depth     mean int. WF      flux WF'

C  optical depth loop

              DO UT = 1, N_OUT_USERTAUS
                WRITE(RUN,'(2x,F9.5,3x,2E17.7)')
     &            USER_TAUS_INPUT(UT),
     &            MINT_SURFACEWF(Z,UT,IB,WDIR),
     &            FLUX_SURFACEWF(Z,UT,IB,WDIR)
              ENDDO

C  close control and loops

            ENDDO
          ENDDO
         ENDDO
        ENDIF

C  End surface weighting function stuff

      ENDIF

C  Surfbb Weighting function output
C  ################################

C      IF ( DO_SURFBB_LINEARIZATION ) THEN

C  PLACEHOLDER

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_WRITEINPUT ( IUNIT )

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  include files with input variables to be written out

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  argument

      INTEGER            IUNIT

C  local variable

      INTEGER            I

C  Linearization control

      WRITE(IUNIT, FMT_HEADING) 'Linearization control'

      IF ( DO_SIMULATION_ONLY ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' Output will be intensity only (no weight. fns.)'
      ELSE
        IF ( DO_PROFILE_LINEARIZATION ) THEN
          WRITE(IUNIT, FMT_CHAR)
     & ' Atmospheric profile weighting functions will be output w.r.t:'
          DO I = 1, N_TOTALATMOS_WFS
            WRITE(IUNIT, FMT_CHAR)PROFILEWF_NAMES(I)
          ENDDO
        ENDIF
        WRITE(IUNIT,'(A)')' '
        IF ( DO_COLUMN_LINEARIZATION ) THEN
          WRITE(IUNIT, FMT_CHAR)
     & 'Atmospheric column weighting functions will be output w.r.t:'
          DO I = 1, N_TOTALCOLUMN_WFS
            WRITE(IUNIT, FMT_CHAR)COLUMNWF_NAMES(I)
          ENDDO
        ENDIF
        WRITE(IUNIT,'(A)')' '
        IF ( DO_SURFACE_LINEARIZATION ) THEN
          WRITE(IUNIT, FMT_CHAR)
     &      ' Surface weighting functions will be output'
        ENDIF
        IF ( DO_SURFBB_LINEARIZATION ) THEN
          WRITE(IUNIT, FMT_CHAR)
     &      ' Surface blackbody weighting functions will be output'
        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_L_WRITESCEN ( SUNIT )

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  Include file of linearization input variables

      INCLUDE '../includes/LIDORT_L_INPUTS.VARS'

C  Module input

      INTEGER             SUNIT

C  local variables

      INTEGER             N, V

C  set up

C   WF input
C   --------

C  layer variational quantities

      IF ( DO_ATMOS_LINEARIZATION ) THEN

        WRITE(SUNIT,FMT_SECTION)
     &     'Atmospheric input for Weighting function calculation'

        WRITE(SUNIT, FMT_HEADING)
     &        'Single scattering albedo variations'
        WRITE(SUNIT,'(a/)')
     &     'Layer varying | parameter numbers-->'
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
             WRITE(SUNIT,'(I3,14x,5(1pe12.4))')
     &        N,(L_OMEGA_TOTAL_INPUT(V,N),V=1,LAYER_VARY_NUMBER(N))
          ENDIF
        ENDDO

        WRITE(SUNIT, FMT_HEADING) 'Optical depth variations'
        WRITE(SUNIT,'(a/)')
     &     'Layer varying | parameter numbers-->'
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            WRITE(SUNIT,'(I3,14x,5(1pe12.4))')
     &         N,(L_DELTAU_VERT_INPUT(V,N),V=1,LAYER_VARY_NUMBER(N))
          ENDIF
        ENDDO

      ENDIF

C  Finish

      RETURN
      END

