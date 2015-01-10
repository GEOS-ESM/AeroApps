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
C #              VLIDORT_WRITEFOURIER                           #
C #              VLIDORT_WRITERESULTS                           #
C #              VLIDORT_WRITEINPUT                             #
C #              VLIDORT_WRITESCEN                              #
C #                                                             #
C ###############################################################

      SUBROUTINE VLIDORT_WRITEFOURIER
     &     ( FUNIT, FOURIER_COMPONENT )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setups

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Include file of result variables

      INCLUDE '../includes/VLIDORT_RESULTS.VARS'

C  input variables
C  ---------------

      INTEGER             FOURIER_COMPONENT, FUNIT

C  local variables
C  ---------------

      INTEGER          I, S, UT, UTA, NT, IB, IDIR, WDIR

      DOUBLE PRECISION DT, USER_HEIGHTS(MAX_USER_LEVELS)
      DOUBLE PRECISION USER_OPDEPS (MAX_USER_LEVELS)

      CHARACTER*15       STOKES_LABEL(MAXSTOKES)
      DATA STOKES_LABEL /
     &          ' I component   ',' Q component   ',
     &          ' U component   ',' V component   '/

C  write header

      WRITE(FUNIT,'(a)')' '
      write(FUNIT,'(/a,i3/a/)')
     &  'Stokes Vector output for Fourier component',
     &   FOURIER_COMPONENT,
     &  '------------------------------------------'

      WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of output levels = ',N_USER_LEVELS
      WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS
      WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of solar zenith angles = ',N_SZANGLES

C  Fix output

      DO UTA = 1, N_USER_LEVELS
        NT = INT(USER_LEVELS(UTA))
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          DT = USER_LEVELS(UTA) - DBLE(NT)
          USER_HEIGHTS(UTA) = HEIGHT_GRID(NT-1) * (ONE-DT)
     &                      + HEIGHT_GRID(NT)   * DT
          USER_OPDEPS(UTA) = TAUGRID_INPUT(NT-1)
     &                       + DELTAU_VERT_INPUT(NT)*DT
        ELSE
          USER_OPDEPS(UTA)  = TAUGRID_INPUT(NT)
          USER_HEIGHTS(UTA) = HEIGHT_GRID(NT)
        ENDIF
      ENDDO

C  detailed output

      DO IDIR = 1, N_DIRECTIONS

        WDIR = WHICH_DIRECTIONS(IDIR)

C  direction header

        IF (WDIR .EQ. UPIDX ) THEN
          WRITE(FUNIT,'(/A)')
     & '  --> Upwelling Stokes Vector, all optical depths & angles'
        ELSE IF (WDIR .EQ. DNIDX ) THEN
            WRITE(FUNIT,'(/A)')
     & '  --> Downwelling Stokes Vector, all optical depths & angles'
        ENDIF

C  output loops

        DO IB = 1, N_SZANGLES
          WRITE(FUNIT,'(/a,2x,f10.4/)') 'SZA (degs)',SZANGLES(IB)
          DO UT = 1, N_USER_LEVELS
            WRITE(FUNIT,'(/a,3(f10.5,2x))')
     *        'Layer #, Height and Optical depth',
     *       USER_LEVELS(UT),USER_HEIGHTS(UT), USER_OPDEPS(UT) 
            write(FUNIT,'(a,4(a15)/)')
     &       'output angle |  ',(STOKES_LABEL(S),S=1,NSTOKES)
            DO I = 1, N_OUT_STREAMS
              WRITE(FUNIT,'(2x,F9.5,2x,5(1PE15.5))')OUT_ANGLES(I),
     &              (STOKES_F(UT,I,IB,S,WDIR),S=1,NSTOKES)
            ENDDO
          ENDDO
        ENDDO

      ENDDO

C  Finish

      END

C

      SUBROUTINE VLIDORT_WRITERESULTS (RUN)

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setups

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  Include file of result variables

      INCLUDE '../includes/VLIDORT_RESULTS.VARS'

C  argument
C  --------

      INTEGER             RUN

C  local variables
C  ---------------

      INTEGER          I, UA, IB, V, UT, S, F, FMAX, NT, UTA
      INTEGER          LOCAL_NUSERAZMS, IDIR, WDIR

      DOUBLE PRECISION DT, USER_HEIGHTS(MAX_USER_LEVELS)
      DOUBLE PRECISION USER_OPDEPS (MAX_USER_LEVELS)

      CHARACTER*15        STOKES_LABEL(MAXSTOKES)
      DATA STOKES_LABEL /
     &          ' I component   ',' Q component   ',
     &          ' U component   ',' V component   '/
      
C  Beam Attenuation summary

      write(RUN,'(/a/a/)')'Results Output summary',
     &                      '----------------------'

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

C  Fix output

      DO UTA = 1, N_USER_LEVELS
        NT = INT(USER_LEVELS(UTA))
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          DT = USER_LEVELS(UTA) - DBLE(NT)
          USER_HEIGHTS(UTA) = HEIGHT_GRID(NT-1) * (ONE-DT)
     &                      + HEIGHT_GRID(NT)   * DT
          USER_OPDEPS(UTA) = TAUGRID_INPUT(NT-1)
     &                       + DELTAU_VERT_INPUT(NT)*DT
        ELSE
          USER_OPDEPS(UTA)  = TAUGRID_INPUT(NT)
          USER_HEIGHTS(UTA) = HEIGHT_GRID(NT)
        ENDIF
      ENDDO

C  control point for avoiding intensity output
         
      IF ( DO_MVOUT_ONLY ) GO TO 400

C  Stokes vector output
C  --------------------

C  local number of azimuths

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_NUSERAZMS = 1
      ELSE
        LOCAL_NUSERAZMS = N_USER_RELAZMS
      ENDIF

C  overall header

      write(RUN,'(/a/a/)')'Stokes vector output',
     &                      '--------------------'

C  start beam loop

      DO IB = 1, NBEAMS

       write(RUN,FMT_REAL)
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)

C  start azimuth loop

       DO UA = 1, LOCAL_NUSERAZMS

C  azimuth angle header

        IF ( DO_NO_AZIMUTH ) THEN
          write(RUN,FMT_CHAR)
     *       '* * Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
        ELSE
          write(RUN,FMT_REAL)
     * '* * RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
        ENDIF

        WRITE(RUN,'(a)')' '
        WRITE(RUN,FMT_INTEGER)
     &      'Total number of output levels = ',N_USER_LEVELS
        WRITE(RUN,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS

C  detailed output

        DO IDIR = 1, N_DIRECTIONS

          WDIR = WHICH_DIRECTIONS(IDIR)

C  direction header

          IF (WDIR .EQ. UPIDX ) THEN
            WRITE(RUN,'(/A)')
     &     '  --> Upwelling intensities all output levels and angles'
          ELSE IF (WDIR .EQ. DNIDX ) THEN
            WRITE(RUN,'(/A)')
     &     '  --> Downwelling intensities all output levels and angles'
          ENDIF

C  output loop

          DO UT = 1, N_USER_LEVELS
            WRITE(RUN,'(/a,3(f10.5,2x))')
     *        'Layer #, Height and Optical depth',
     *       USER_LEVELS(UT),USER_HEIGHTS(UT), USER_OPDEPS(UT)
            write(RUN,'(a,4(a15))')
     &         'output angle |  ',(STOKES_LABEL(S),S=1,NSTOKES)
            DO I = 1, N_OUT_STREAMS
              V = VZA_OFFSETS(IB,I) + UA
              WRITE(RUN,'(2x,F9.5,2x,5(1PE15.5))')OUT_ANGLES(I),
     &              (STOKES(UT,V,S,WDIR),S=1,NSTOKES)
            ENDDO
          ENDDO

C  direction and solar/azimuth loops - end

        ENDDO
       ENDDO
      ENDDO

C  integrated value output
C  -----------------------

400   CONTINUE
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

       write(RUN,'(/a/a)')'integrated value output',
     &                       '-----------------------'

       DO IB = 1, NBEAMS

        write(RUN,FMT_REAL)
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
        write(RUN,'(a)')' '

C  detailed output

        DO IDIR = 1, N_DIRECTIONS

          WDIR = WHICH_DIRECTIONS(IDIR)

C  direction header

          IF (WDIR .EQ. UPIDX ) THEN
            WRITE(RUN,'(/A/)')
     &'  --> Upwelling mean values & fluxes, all output levels'
          ELSE IF (WDIR .EQ. DNIDX ) THEN
            WRITE(RUN,'(/A/)')
     &'  --> Downwelling mean values & fluxes, all output levels'
          ENDIF

C  Mean values

          write(RUN,'(a/)')
     &      '      ** Mean-values for Stokes components  ----> '
          write(RUN,'(a13,3x,4(a15)/)')'output levels',
     &              (STOKES_LABEL(S),S=1,NSTOKES)
          DO UT = 1, N_USER_LEVELS
            WRITE(RUN,'(2x,F9.5,3x,1p4e15.5)')USER_LEVELS(UT),
     &                 (MEAN_STOKES(UT,IB,S,WDIR),S=1,NSTOKES)
          ENDDO

C  fluxes

          write(RUN,'(/a/)')
     &      '      ** Fluxes for Stokes components  ----> '
          write(RUN,'(a13,3x,4(a15)/)')'output levels',
     &              (STOKES_LABEL(S),S=1,NSTOKES)
          DO UT = 1, N_USER_LEVELS
            WRITE(RUN,'(2x,F9.5,3x,1p4e15.5)')USER_LEVELS(UT),
     &                 (FLUX_STOKES(UT,IB,S,WDIR),S=1,NSTOKES)
          ENDDO

C  end direction loop

        ENDDO

       ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_WRITEINPUT ( IUNIT )

C  Include files
C  -------------

C  Write to file of all control input

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  argument
C  --------

      INTEGER             IUNIT

C  local variables

      INTEGER             BRDF

C  heading and version number

      WRITE(IUNIT, FMT_SECTION)
     &     ' Control input variables for run of VLIDORT'
      WRITE(IUNIT, FMT_CHAR)
     &     ' VLIDORT Version number = ',VLIDORT_VERSION_NUMBER

C  general control

      WRITE(IUNIT, FMT_HEADING) ' control integers'

      WRITE(IUNIT, FMT_INTEGER)
     &   ' Double Gaussian quadrature over [-1,0] & [0,1]'
      WRITE(IUNIT, FMT_INTEGER)
     &     '  ** Number of half space streams = ', NSTREAMS
      WRITE(IUNIT, FMT_INTEGER)
     &     '  ** Number of atmospheric layers = ', NLAYERS
      WRITE(IUNIT, FMT_INTEGER)
     &   '  ** Number of Greek moments (input) = ',
     &             NGREEK_MOMENTS_INPUT

C      IF ( DO_THERMAL_EMISSION ) THEN
C        WRITE(IUNIT, FMT_INTEGER)
C     &  '  ** Number of thermal emission coefficients = ',
C     &       N_THERMAL_COEFFS
C      ENDIF

      WRITE(IUNIT, FMT_HEADING) ' flux/accuracy control'
      WRITE(IUNIT, FMT_REAL)
     &     ' Flux constant = ', FLUX_FACTOR
      IF ( .NOT.DO_NO_AZIMUTH ) THEN
        WRITE(IUNIT, FMT_REAL)
     & ' accuracy criterion (Fourier convergence) = ',VLIDORT_ACCURACY
      ELSE
        WRITE(IUNIT, FMT_CHAR)
     &  '  ** No Fourier series -- Azimuth-independent term only'
      ENDIF

C  RTE control

      WRITE(IUNIT, FMT_HEADING) ' RTE solution control'

      IF ( DO_DIRECT_BEAM ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' Direct beam will be included in solution'
      ELSE
        WRITE(IUNIT, FMT_CHAR)
     &      ' Direct beam will NOT be included in solution'
      ENDIF

C      IF ( DO_THERMAL_EMISSION ) THEN
C        WRITE(IUNIT, FMT_CHAR)
C     &      ' Thermal Emission will be included in solution'
C      ELSE
C        WRITE(IUNIT, FMT_CHAR)
C     &      ' NO Thermal emission in the solution'
C      ENDIF

      IF ( DO_SURFACE_EMISSION ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' Surface Thermal Emission will be included in solution'
      ELSE
        WRITE(IUNIT, FMT_CHAR)
     &      ' NO Surface Thermal emission in the solution'
      ENDIF

C  surface input write

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' Surface will be treated as Lambertian'
      ELSE
        WRITE(IUNIT, FMT_CHAR)
     &    ' Bidirectional surface input: Selected BRDF kernels are:'
        DO BRDF = 1, N_BRDF_KERNELS
          WRITE(IUNIT, '(20X,A4,A10)')' ** ',BRDF_NAMES(BRDF)
        ENDDO
      ENDIF

C  Other writes. Isotropic-only option now removed. 01/17/06.

c      IF ( DO_ISOTROPIC_ONLY ) THEN
c        WRITE(IUNIT, FMT_CHAR)' Medium is isotropic'
c      ENDIF

      IF ( DO_RAYLEIGH_ONLY ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' Medium has Rayleigh scattering only'
      ELSE
        WRITE(IUNIT, FMT_CHAR)
     &      ' Medium has general scattering phase function'
      ENDIF

C      IF ( DO_TRANSMITTANCE_ONLY ) THEN
C        WRITE(IUNIT, FMT_CHAR)
C     &      ' RTE solution is Transmission-only (no scattering)'
C      ELSE
C        WRITE(IUNIT, FMT_CHAR)
C     &      ' RTE solution is full multiple-scattering'
C      ENDIF

C  Beam particular integral control

      WRITE(IUNIT, FMT_HEADING) 'Beam solution control'

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &  ' Beam solution determined by classical (Chandrasekhar) method'
      ELSE
        WRITE(IUNIT, FMT_CHAR)
     &     ' Beam solution determined by Greens function method'
      ENDIF

      IF ( DO_PLANE_PARALLEL ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' Beam solution in Plane-parallel approximation'
        ELSE
        WRITE(IUNIT, FMT_CHAR)
     &      ' Beam solution in Pseudo-spherical approximation'
        WRITE(IUNIT, FMT_CHAR)
     &  '  ** attenuation will use average secant estimation'
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          WRITE(IUNIT, FMT_CHAR)
     &  '  ** Refractive geometry'
        ELSE
          WRITE(IUNIT, FMT_CHAR)
     &  '  ** Straight line geometry'
        ENDIF
      ENDIF

C  output control

      WRITE(IUNIT, FMT_HEADING) 'Output control'

      IF ( DO_ADDITIONAL_MVOUT ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' Output of intensities AND fluxes & mean intensities'
      ELSE
        IF ( DO_MVOUT_ONLY ) THEN
          WRITE(IUNIT, FMT_CHAR)
     &      ' Output for fluxes & mean intensities ONLY'
        ELSE
          WRITE(IUNIT, FMT_CHAR)
     &      ' Output for intensities ONLY'
        ENDIF
      ENDIF

      WRITE(IUNIT, FMT_INTEGER)
     &     ' Number of Solar Zenith angles = ', NBEAMS
      
      IF ( .NOT.DO_NO_AZIMUTH ) THEN
        WRITE(IUNIT, FMT_INTEGER)
     &     ' Number of user-defined azimuth angles = ', N_USER_RELAZMS
      ENDIF

      WRITE(IUNIT, FMT_INTEGER)
     &       ' Total number of output levels = ',
     &        N_USER_LEVELS

      IF ( DO_USER_STREAMS ) THEN
        IF ( DO_QUAD_OUTPUT ) THEN
          WRITE(IUNIT, FMT_CHAR)
     &       ' Stream angle output will include quadratures'
          WRITE(IUNIT, FMT_INTEGER)
     &       '  ** Total Number of output stream angles = ',
     &                    NSTREAMS + N_USER_STREAMS
        ELSE
          WRITE(IUNIT, FMT_CHAR)
     &       ' Stream angle output for user-defined angles only'
          WRITE(IUNIT, FMT_INTEGER)
     &       '  ** Total Number of output stream angles = ',
     &                    N_USER_STREAMS
        ENDIF
      ELSE
        IF ( .NOT. DO_MVOUT_ONLY ) THEN
          WRITE(IUNIT, FMT_CHAR)
     &      ' Stream output at Quadrature angles only'
          WRITE(IUNIT, FMT_INTEGER)
     &       '  ** Total Number of output stream angles = ',NSTREAMS
        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE VLIDORT_WRITESCEN ( SUNIT )

C  Include files
C  -------------

C  write to file of all geophysical VLIDORT input

      INCLUDE '../includes/VLIDORT.PARS'

C  Include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of setups

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'

C  argument

      INTEGER             SUNIT

C  local variables

      INTEGER             S, N, L, BRDF

C  heading and version number

      WRITE(SUNIT,FMT_SECTION)
     &     ' Geophysical scenario variables for run of VLIDORT'
      WRITE(SUNIT,FMT_CHAR)
     &     ' VLIDORT Version number = ',VLIDORT_VERSION_NUMBER

C   basic input
C   -----------

      WRITE(SUNIT,FMT_SECTION)
     &     ' Basic atmospheric input for stokes vector calculation'

C  Layer optical depths and single scatter albedos

      IF ( DO_DELTAM_SCALING ) THEN
        WRITE(SUNIT, FMT_HEADING)
     &      '  scaled and unscaled inputs with Delta-M turned ON'
        WRITE(SUNIT, FMT_HEADING)
     &      'Layer optical depths and single-scatter-albedos '
        WRITE(SUNIT,'(a,a25/a,a25)')
     &       '      (unscaled) (scaled) ','  (unscaled)   (scaled)  ',
     &       'Layer  Op-depth  Op-depth ','  s.s albedo  s.s albedo '
        WRITE(SUNIT,'(a)')' '
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(I3,T6,2(f9.5,1x),2x,2(f9.5,3x))')
     &            N,TAUGRID_INPUT(N),TAUGRID(N),
     &            OMEGA_TOTAL_INPUT(N),OMEGA_TOTAL(N)
        ENDDO
      ELSE
        WRITE(SUNIT, FMT_HEADING)
     &      '  Unscaled inputs only: Delta-M turned OFF'
        WRITE(SUNIT, FMT_HEADING)
     &    'Layer optical depths and single-scatter-albedos '
        WRITE(SUNIT,'(a,a13/a,a13)')
     &        '      (unscaled) ','  (unscaled) ',
     &        'Layer  Op-depth  ','  s.s albedo '
        WRITE(SUNIT,'(a)')' '
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(I3,T6,f9.5,4x,f9.5)')
     &            N,TAUGRID_INPUT(N),OMEGA_TOTAL_INPUT(N)
        ENDDO
      ENDIF

C   phase matrix elements

      IF ( DO_DELTAM_SCALING ) THEN
        WRITE(SUNIT, FMT_HEADING)
     &       'Phase matrix elements (Greek constants), SCALED'
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(a, I3/)')
     &              'Matrix elements (A-F) for layer ',N
          DO L = 0, NMOMENTS
            WRITE(SUNIT,'(I2,1x,1p6e12.4)')
     &     L,(GREEKMAT_TOTAL(L,N,GREEKMAT_INDEX(S)),S = 1, 6)
          ENDDO
        ENDDO
      ELSE
        WRITE(SUNIT, FMT_HEADING)
     &       'Phase matrix elements (Greek constants), UNSCALED'
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(a, I3/)')
     &              'Matrix elements (A-F) for layer ',N
          DO L = 0, NGREEK_MOMENTS_INPUT
            WRITE(SUNIT,'(I2,1x,1p6e12.4)')
     &    L,(GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(S)),S = 1, 6)
          ENDDO
        ENDDO
      ENDIF

C  Commented out thermal expansion coefficients
C      IF ( DO_THERMAL_EMISSION ) THEN
C        WRITE(SUNIT,FMT_HEADING)'thermal emission coefficients'
C        WRITE(SUNIT,'(a/)')
C     &       'Layer | thermal emission expansion coeffs-->'
C        DO N = 1, NLAYERS
C          WRITE(SUNIT,'(I3,4x,10(f10.5))')
C     &          N,(THERMAL_COEFFS(N,S),S=1,N_THERMAL_COEFFS)
C        ENDDO
C      ENDIF

C  surface property

      WRITE(SUNIT,FMT_HEADING)'Surface reflecting property'

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        WRITE(SUNIT,FMT_REAL)
     &       '(Lambertian) surface albedo is',LAMBERTIAN_ALBEDO
        IF ( DO_SURFACE_EMISSION ) THEN
          WRITE(SUNIT,FMT_REAL)
     &      '(Lambertian) emissivity = ',ONE-LAMBERTIAN_ALBEDO
          WRITE(SUNIT,FMT_REAL)
     &      'Surface blackbody function is', SURFBB
        ENDIF
      ELSE
        WRITE(SUNIT,FMT_CHAR)
     &         ' BRDF kernels and parameters are:'
        WRITE(SUNIT,'(A)')
     &  ' Kernels    Amplitudes    Parameters --->'
        DO BRDF = 1, N_BRDF_KERNELS
          WRITE(SUNIT, '(A10,2x,f10.5,2x,1p3e12.4)')
     &     BRDF_NAMES(WHICH_BRDF(BRDF)),  BRDF_FACTORS(BRDF),
     &     (BRDF_PARAMETERS(BRDF,N),N=1,N_BRDF_PARAMETERS(BRDF))
        ENDDO
      ENDIF

C   Geometry input
C   --------------

      WRITE(SUNIT,FMT_SECTION)' Geometry input'

      WRITE(SUNIT,FMT_HEADING)
     &          'Solar zenith angles'
      WRITE(SUNIT,'(a/)')'Number |   Angle'
      DO N = 1, N_SZANGLES
        WRITE(SUNIT,'(I3,5x,1x,f9.5)') N,SZANGLES(N)
      ENDDO

      WRITE(SUNIT,FMT_HEADING)
     &         'Computational (quadrature) angles in the half space'
      WRITE(SUNIT,'(a/)')
     &           'Stream  |  Angle    |   Cosine   |   Weight'
      DO N = 1, NSTREAMS
        WRITE(SUNIT,'(I3,5x,3(1x,f9.5,3x))')N,QUAD_ANGLES(N),
     &        QUAD_STREAMS(N),QUAD_WEIGHTS(N)
      ENDDO

      WRITE(SUNIT,FMT_HEADING) 'User-defined viewing zenith angles'
      WRITE(SUNIT,'(a/)')'Number  |  Angle    |   Cosine'
      DO N = 1, N_USER_VZANGLES
        WRITE(SUNIT,'(I3,5x,2(1x,f9.5,3x))')
     *  N,USER_VZANGLES(N),DCOS(USER_VZANGLES(N)*DEG_TO_RAD)
      ENDDO

      IF ( DO_NO_AZIMUTH ) THEN
        WRITE(SUNIT,FMT_HEADING) 'No azimuth angles'
      ELSE
        WRITE(SUNIT,FMT_HEADING)
     &          'User-defined relative azimuth angles'
        WRITE(SUNIT,'(a/)')'Number |   Angle'
        DO N = 1, N_USER_RELAZMS
          WRITE(SUNIT,'(I3,5x,1x,f9.5)') N,USER_RELAZMS(N)
        ENDDO
      ENDIF

C  Output levels
C  -------------

      WRITE(SUNIT,FMT_SECTION)
     &       ' User-defined levels for output'
      WRITE(SUNIT,'(a/)')' # | Level/Layer of occurrence'
      DO N = 1, N_USER_LEVELS
        WRITE(SUNIT,'(i3,4x,f9.5)')N, USER_LEVELS(N)
      ENDDO

C  Finish

      RETURN
      END

