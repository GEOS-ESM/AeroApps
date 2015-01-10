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
C #              VLIDORT_L_WRITEFOURIER                         #
C #              VLIDORT_L_WRITERESULTS                         #
C #              VLIDORT_L_WRITEINPUT                           #
C #              VLIDORT_L_WRITESCEN                            #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_L_WRITEFOURIER
     &     ( DO_INCLUDE_SURFACE, FUNIT, FOURIER_COMPONENT )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of linearization input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  include file with results for writing

      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  input variables
C  ---------------

      INTEGER          FOURIER_COMPONENT, FUNIT
      LOGICAL          DO_INCLUDE_SURFACE
      
C  local variables
C  ---------------

      INTEGER          I, IDIR, UT, WDIR, Z, Q, K, IB, O1

C  Profile Atmospheric Weighting function output
C  #############################################

      IF ( DO_PROFILE_LINEARIZATION ) THEN

C  overall header

        WRITE(FUNIT,'(a)')' '
        write(FUNIT,'(/a,i3/a/)')
     & 'Atmospheric Profile Weighting function: Fourier component',
     &         FOURIER_COMPONENT,
     & '--------------------------------------------------------------'

        WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of output levels = ',N_USER_LEVELS
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
              WRITE(FUNIT,'(/a,2x,f10.4/)')'SZA (degs)',SZANGLES(IB)
              write(FUNIT,'(a/)')
     &             'output angles | output levels values-----> '
              WRITE(FUNIT,'(14X,50(1pe15.5))')
     *          ( USER_LEVELS(UT),UT = 1, N_USER_LEVELS)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_OUT_STREAMS
                  WRITE(FUNIT,'(3x,F9.5,2x,I2,50(1PE15.5))')
     *                 OUT_ANGLES(I),UT,
     &           (ATMOSWF_F(Q,K,UT,I,IB,O1,WDIR), O1 = 1, NSTOKES)
                ENDDO
              ENDDO
             ENDDO

C  close control and loops

            ENDDO
           ENDIF
          ENDDO
        ENDDO

C  End profile weighting functions

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
     &      'Total number of output levels = ',N_USER_LEVELS
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
              WRITE(FUNIT,'(/a,2x,f10.4/)')'SZA (degs)',SZANGLES(IB)
              write(FUNIT,'(a/)')
     &             'output angles | Output level values-----> '
              WRITE(FUNIT,'(14X,50(1pe15.5))')
     *          ( USER_LEVELS(UT),UT = 1, N_USER_LEVELS)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_OUT_STREAMS
                  WRITE(FUNIT,'(3x,F9.5,2x,I2,50(1PE15.5))')
     *                 OUT_ANGLES(I),UT,
     &           (ATMOSWF_F(Q,0,UT,I,IB,O1,WDIR), O1 = 1, NSTOKES)
                ENDDO
              ENDDO
            ENDDO

C  close control and loops

          ENDDO
        ENDDO

C  End atmospheric  bulk/column weighting function stuff

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
     &      'Total number of output levels  = ',N_USER_LEVELS
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
              WRITE(FUNIT,'(/a,2x,f10.4/)')'SZA (degs)',SZANGLES(IB)
              write(FUNIT,'(a/)')
     &             'output angles | Output level values-----> '
              WRITE(FUNIT,'(14X,50(1pe15.5))')
     *          ( USER_LEVELS(UT),UT = 1, N_USER_LEVELS)
              DO I = 1, N_OUT_STREAMS
c               WRITE(FUNIT,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
c     &           (SURFACEWF_F(Z,UT,I,IB,WDIR),UT = 1, N_USER_LEVELS)
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

      SUBROUTINE VLIDORT_L_WRITERESULTS (RUN)

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of linearization input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

C  include files with results for writing

      INCLUDE '../includes/VLIDORT_RESULTS.VARS'
      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  input variables
C  ---------------

      INTEGER      RUN

C  local variables
C  ---------------

      INTEGER      I, UA, LOCAL_NUSERAZMS, IDIR, UT, WDIR, Z, O1
      INTEGER      Q, K, IB, V, FMAX, F
      INTEGER      VINDEX(MAXBEAMS,MAX_USER_STREAMS,MAX_USER_RELAZMS)
      CHARACTER*11 STOKESNAM(4)
      DATA STOKESNAM / 'I-component','Q-component',
     &                 'U-component','V-component'/

C  Beam Attenuation summary
C  ========================

      write(RUN,'(/a/a/)')'Results Output summary',
     &                     '----------------------'

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
     &       'Full Stokes calculation has been performed'
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
     &   ' --> Accuracy level was pre-set at : ',VLIDORT_ACCURACY

        write(RUN,'(a,I5)')' - Number Solar zenith angles : ',NBEAMS
        write(run,'(/a/)')
     &         '   SZA  | # Fourier | Fourier breakdown -->'
        fmax = -100
        DO ib = 1, nbeams
          fmax = MAX(fmax,fourier_saved(ib))
        END DO
        DO ib = 1, nbeams
           WRITE(run,'(1x,f7.2,4x,i3,6x,50(L1))')
     &            szangles(ib),fourier_saved(ib),
     &            (do_multibeam(ib,f),f=0,fmax)
        END DO
      ELSE
        write(RUN,'(/a)')
     &   'Azimuth independent output only (Fourier = 0)'
      ENDIF

C  local number of azimuths

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_NUSERAZMS = 1
      ELSE
        LOCAL_NUSERAZMS = N_USER_RELAZMS
      ENDIF

C  Indexing

      DO IB = 1, NBEAMS
       DO I = 1, N_OUT_STREAMS
        DO UA = 1, LOCAL_NUSERAZMS
         VINDEX(IB,I,UA) = VZA_OFFSETS(IB,I) + UA
        ENDDO
       ENDDO
      ENDDO

C  Atmospheric Profile Weighting function output
C  #############################################

      IF ( DO_PROFILE_LINEARIZATION ) THEN

C  control point for avoiding weighting function output
         
        IF ( DO_MVOUT_ONLY ) GO TO 401

C  overall header

        write(RUN,'(/a/a/)')
     &       'Atmospheric Profile Weighting function output',
     &       '---------------------------------------------'

C  start beam loop

        DO IB = 1, NBEAMS
         write(RUN,FMT_REAL)
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
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
     &      'Total number of output levels = ',N_USER_LEVELS
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

              DO O1 = 1, NSTOKES
               write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
               write(RUN,'(a)') 'output angles | output levels ---> '
               WRITE(RUN,'(14x,50(2x,1pe11.4,2x))')
     &           (USER_LEVELS(UT), UT = 1, N_USER_LEVELS)
               DO I = 1, N_OUT_STREAMS
                V = VZA_OFFSETS(IB,I) + UA
                WRITE(RUN,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
     &         (PROFILEWF(Q,K,UT,V,O1,WDIR), UT = 1, N_USER_LEVELS)
               ENDDO
              ENDDO

C  close control and loops

             ENDDO
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO

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
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
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

C  output

               DO O1 = 1, NSTOKES
                write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
                WRITE(RUN,'(a/)')
     &      'output level     mean int. WF      flux WF'
                DO UT = 1, N_USER_LEVELS
                WRITE(RUN,'(2x,F9.5,3x,2E17.7)')
     &                    USER_LEVELS(UT),
     &                    MINT_ATMOSWF(Q,K,UT,IB,O1,WDIR),
     &                    FLUX_ATMOSWF(Q,K,UT,IB,O1, WDIR)
                ENDDO
               ENDDO

C  close control and loops

              ENDDO
             ENDIF
            ENDDO
           ENDDO
          ENDDO
        ENDIF

C  End profile atmospheric weighting function stuff

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
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
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
     &      'Total number of output levels = ',N_USER_LEVELS
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

              DO O1 = 1, NSTOKES
               write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
               write(RUN,'(a)') 'output angles | output levels ---> '
               WRITE(RUN,'(14x,50(2x,1pe11.4,2x))')
     &           (USER_LEVELS(UT), UT = 1, N_USER_LEVELS)
               DO I = 1, N_OUT_STREAMS
                V = VZA_OFFSETS(IB,I) + UA
                WRITE(RUN,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
     &         (COLUMNWF(Q,UT,V,O1,WDIR), UT = 1, N_USER_LEVELS)
               ENDDO
              ENDDO

C  close control and loops

             ENDDO
           ENDDO
         ENDDO
        ENDDO

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
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
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

C  Main output

               DO O1 = 1, NSTOKES
                write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
                 WRITE(RUN,'(a/)')
     &      'output level     mean int. WF      flux WF'
                 DO UT = 1, N_USER_LEVELS
                  WRITE(RUN,'(2x,F9.5,3x,2E17.7)')
     &                    USER_LEVELS(UT),
     &                    MINT_ATMOSWF(Q,0,UT,IB,O1,WDIR),
     &                    FLUX_ATMOSWF(Q,0,UT,IB,O1,WDIR)
                 ENDDO
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
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
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
     &      'Total number of output levels = ',N_USER_LEVELS
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

            DO O1 = 1, NSTOKES
             write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
             write(RUN,'(a)') 'output angles | output levels ---> '
             WRITE(RUN,'(14x,50(2x,1pe11.4,2x))')
     &         (USER_LEVELS(UT), UT = 1, N_USER_LEVELS)
             DO I = 1, N_OUT_STREAMS
              V = VZA_OFFSETS(IB,I) + UA
              WRITE(RUN,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
     &         (SURFACEWF(Z,UT,V,O1,WDIR), UT = 1, N_USER_LEVELS)
             ENDDO
            ENDDO

C  close control and loops

           ENDDO
          ENDDO
         ENDDO
        ENDDO

C  Mean-value output
C  -----------------

 402     continue

        IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

         write(RUN,'(/a/a)')
     &          'Mean value surface weighting function output',
     &          '--------------------------------------------'

C  start beam loop

         DO IB = 1, NBEAMS

         write(RUN,FMT_REAL)
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
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

C  output level header

              WRITE(RUN,'(a,I2/a/)')
     &      '  --> For Surface WF number :',Z,
     &      'output level     mean int. WF      flux WF'

C  output level loop

              DO O1 = 1, NSTOKES
               write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
               DO UT = 1, N_USER_LEVELS
                WRITE(RUN,'(2x,F9.5,3x,2E17.7)')
     &            USER_LEVELS(UT),
     &            MINT_SURFACEWF(Z,UT,IB,O1,WDIR),
     &            FLUX_SURFACEWF(Z,UT,IB,O1,WDIR)
               ENDDO
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

      SUBROUTINE VLIDORT_L_WRITEINPUT ( IUNIT )

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files with input variables to be written out

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

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

      SUBROUTINE VLIDORT_L_WRITESCEN ( SUNIT )

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  Include file of linearization input variables

      INCLUDE '../includes/VLIDORT_L_INPUTS.VARS'

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

        WRITE(SUNIT, FMT_HEADING) 'Output level variations'
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

