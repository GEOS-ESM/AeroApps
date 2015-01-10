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
C #              LIDORT_WRITEFOURIER                            #
C #              LIDORT_WRITERESULTS                            #
C #              LIDORT_WRITEINPUT                              #
C #              LIDORT_WRITESCEN                               #
C #                                                             #
C ###############################################################

      SUBROUTINE LIDORT_WRITEFOURIER
     &     ( FUNIT, FOURIER_COMPONENT )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input/output variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'
      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  input variables
C  ---------------

      INTEGER             FOURIER_COMPONENT, FUNIT

C  local variables
C  ---------------

      INTEGER            I, UT, IB, IDIR, WDIR

C  write header

      WRITE(FUNIT,'(a)')' '
      write(FUNIT,'(/a,i3/a/)')
     &  'Intensity output for Fourier component',FOURIER_COMPONENT,
     &  '--------------------------------------'

      WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
      WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS
      WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of solar zenith angles = ',NBEAMS

C  detailed output

      DO IDIR = 1, N_DIRECTIONS

        WDIR = WHICH_DIRECTIONS(IDIR)

C  direction header

        IF (WDIR .EQ. UPIDX ) THEN
          WRITE(FUNIT,'(/A)')
     & '  --> Upwelling Radiances, all optical depths & angles'
        ELSE IF (WDIR .EQ. DNIDX ) THEN
            WRITE(FUNIT,'(/A)')
     & '  --> Downwelling Radiances, all optical depths & angles'
        ENDIF

C  output loops

        DO IB = 1, NBEAMS
          WRITE(FUNIT,'(/a,2x,f10.4/)') 'SZA (degs)',BEAM_SZAS(IB)
          write(FUNIT,'(a/)')
     &             'output angles | Optical depth values-----> '
          WRITE(FUNIT,'(14X,50(1pe15.5))')
     *      ( USER_TAUS_INPUT(UT),UT = 1, N_OUT_USERTAUS)
          DO I = 1, N_OUT_STREAMS
            WRITE(FUNIT,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
     &           (INTENSITY_F(UT,I,IB,WDIR),UT = 1, N_OUT_USERTAUS)
          ENDDO
        ENDDO

      ENDDO

C  Finish

      END

C

      SUBROUTINE LIDORT_WRITERESULTS (RUN)

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input/output variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'
      INCLUDE '../includes/LIDORT_RESULTS.VARS'

C  argument
C  --------

      INTEGER             RUN

C  local variables
C  ---------------

      INTEGER             I, UA, IB, UT, F, FMAX, V
      INTEGER             LOCAL_NUSERAZMS, IDIR, WDIR

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

C  control point for avoiding intensity output
         
      IF ( DO_MVOUT_ONLY ) GO TO 400

C  Intensity output
C  ----------------

C  local number of azimuths

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_NUSERAZMS = 1
      ELSE
        LOCAL_NUSERAZMS = N_USER_RELAZMS
      ENDIF

C  overall header

      write(RUN,'(/a/a/)')'Radiance output',
     &                    '---------------'

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
     *       '* * Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
        ELSE
          write(RUN,FMT_REAL)
     * '* * RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
        ENDIF

        WRITE(RUN,'(a)')' '
        WRITE(RUN,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
        WRITE(RUN,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS

C  detailed output

        DO IDIR = 1, N_DIRECTIONS

          WDIR = WHICH_DIRECTIONS(IDIR)

C  direction header

          IF (WDIR .EQ. UPIDX ) THEN
            WRITE(RUN,'(/A)')
     &     '  --> Upwelling intensities all optical depths and angles'
          ELSE IF (WDIR .EQ. DNIDX ) THEN
            WRITE(RUN,'(/A)')
     &     '  --> Downwelling intensities all optical depths and angles'
          ENDIF

C  output loop

          write(RUN,'(a)') 'output angles | optical depths ---> '
          WRITE(RUN,'(14x,50(2x,1pe11.4,2x))')
     &           (USER_TAUS_INPUT(UT), UT = 1, N_OUT_USERTAUS)
          DO I = 1, N_OUT_STREAMS
            V = UMOFF(IB,I) + UA
            WRITE(RUN,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I),
     &         (INTENSITY(UT,V,WDIR), UT = 1, N_OUT_USERTAUS)
          ENDDO

C  direction and solar/azimuth loops - end

        ENDDO
       ENDDO
      ENDDO

C ------------------------------------- Section removed, 2 April 2007
C  Multiple scatter layer source term output
C  -----------------------------------------
c      IF ( .NOT. SAVE_LAYER_MSST ) GO TO 400
C  overall header
c      write(RUN,'(/a/a/)')'Multiple scatter source term output',
c     &                    '-----------------------------------'
C  start azimuth, SZA loops
c      DO IB = 1, NBEAMS
c       write(RUN,FMT_REAL)
c     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',BEAM_SZAS(IB)
c       write(RUN,'(a)')' '
c       DO UA = 1, LOCAL_NUSERAZMS
C  azimuth angle header
c        IF ( DO_NO_AZIMUTH ) THEN
c          write(RUN,FMT_CHAR)
c     *       '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
c        ELSE
c          write(RUN,FMT_REAL)
c     * '* * RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
c        ENDIF
C  User angle header
c        DO I = 1, N_USER_STREAMS
c         V = UMOFF(IB,I) + UA
c         write(RUN,FMT_REAL)
c     * '* * RESULTS FOR USER ANGLE (degs)=',USER_ANGLES_INPUT(I)
C  detailed output
c         DO IDIR = 1, N_DIRECTIONS
c          WDIR = WHICH_DIRECTIONS(IDIR)
c          IF (WDIR .EQ. UPIDX ) THEN
c            WRITE(RUN,'(/A)')
c     &  '  --> Upwelling mult-scat layer source terms'
c            DO N = N_LAYERSOURCE_UP, NLAYERS
c              WRITE(RUN,'(1X,I3,3X,10(1PE15.5))') N,
c     &              INTENSITY_MSST_LAYER(N,V,WDIR)
c            ENDDO
c            WRITE(RUN,'(/A)')
c     &     '  --> Upwelling BOA multiple source terms'
c            WRITE(RUN,'(A7,10(1PE15.5))')'diffuse',
c     &               INTENSITY_MSST_BOA(V)
c            WRITE(RUN,'(A7,10(1PE15.5))')'direct ',
c     &               INTENSITY_DBST_BOA(V)
c          ELSE IF (WDIR .EQ. DNIDX ) THEN
c            WRITE(RUN,'(/A)')
c     &  '  --> Downwelling mult-scat layer source terms'
c            DO N = 1, N_LAYERSOURCE_DN
c              WRITE(RUN,'(1X,I3,3X,10(1PE15.5))') N,
c     &              INTENSITY_MSST_LAYER(N,V,WDIR)
c            ENDDO
c          ENDIF
C  Loops closed
c         ENDDO
c        ENDDO
c       ENDDO
c      ENDDO
C 400  CONTINUE
C ------------------------------------- Section removed, 2 April 2007

C  integrated value output
C  -----------------------

 400  CONTINUE
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

       write(RUN,'(/a/a)')'integrated value output',
     &                       '-----------------------'

       DO IB = 1, NBEAMS

        write(RUN,FMT_HEADING)
     * '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',BEAM_SZAS(IB)

C  detailed output

        DO IDIR = 1, N_DIRECTIONS

          WDIR = WHICH_DIRECTIONS(IDIR)

C  direction header

          IF (WDIR .EQ. UPIDX ) THEN
            WRITE(RUN,'(/A/)')
     &'  --> Upwelling mean values & fluxes, all optical depths'
          ELSE IF (WDIR .EQ. DNIDX ) THEN
            WRITE(RUN,'(/A/)')
     &'  --> Downwelling mean values & fluxes, all optical depths'
          ENDIF

C  Mean intensity and Flux integrals

          write(RUN,'(a/)')' ** Mean intensities and Flux integrals'
          write(RUN,'(a)')'optical depth |   Mean Value     Intensity'
          DO UT = 1, N_OUT_USERTAUS
            WRITE(RUN,'(2x,F9.5,3x,1p2e15.5)')USER_TAUS_INPUT(UT),
     &     MEAN_INTENSITY(UT,IB,WDIR),FLUX_INTEGRAL(UT,IB,WDIR)
          ENDDO

C  end direction loop

        ENDDO

       ENDDO
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LIDORT_WRITEINPUT ( IUNIT )

C  Read, check and write to file of all control input

C  Include files
C  -------------

C  include file of parameters

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'

C  argument
C  --------

      INTEGER             IUNIT

C  local variables

      INTEGER             BRDF

C  heading and version number

      WRITE(IUNIT, FMT_SECTION)
     &     ' Control input variables for run of LIDORT'
      WRITE(IUNIT, FMT_CHAR)
     &     ' LIDORT Version number = ',LIDORT_VERSION_NUMBER

C  general control

      WRITE(IUNIT, FMT_HEADING) ' control integers'

      WRITE(IUNIT, FMT_INTEGER)
     &   ' Double Gaussian quadrature over [-1,0] & [0,1]'
      WRITE(IUNIT, FMT_INTEGER)
     &     '  ** Number of half space streams = ', NSTREAMS
      WRITE(IUNIT, FMT_INTEGER)
     &     '  ** Number of atmospheric layers = ', NLAYERS
      WRITE(IUNIT, FMT_INTEGER)
     &   '  ** Number of Phase function moments (input) = ',
     &             NMOMENTS_INPUT

      WRITE(IUNIT, FMT_HEADING) ' flux/accuracy control'
      WRITE(IUNIT, FMT_REAL)
     &     ' Flux constant = ', FLUX_FACTOR
      IF ( .NOT.DO_NO_AZIMUTH ) THEN
        WRITE(IUNIT, FMT_REAL)
     & ' accuracy criterion (Fourier convergence) = ',LIDORT_ACCURACY
      ELSE
        WRITE(IUNIT, FMT_CHAR)
     &  '  ** No Fourier series -- Azimuth-independent term only'
      ENDIF

C  Phase function control

      WRITE(IUNIT, FMT_HEADING) 'Phase function control'

      IF ( DO_ISOTROPIC_ONLY ) THEN
        WRITE(IUNIT, FMT_CHAR)' Medium is isotropic'
      ENDIF

      IF ( DO_RAYLEIGH_ONLY ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' Medium has Rayleigh scattering only'
      ENDIF

      IF ( .NOT.DO_RAYLEIGH_ONLY.AND..NOT.DO_ISOTROPIC_ONLY ) THEN
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

      IF ( DO_SOLAR_SOURCES ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' solar sources will be included in solution'
      ELSE
        WRITE(IUNIT, FMT_CHAR)
     &      ' Solar sources will NOT be included in solution'
      ENDIF

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_CLASSICAL_SOLUTION ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &  ' Beam solutions determined by classical (Chandrasekhar) method'
       ELSE
        WRITE(IUNIT, FMT_CHAR)
     &     ' Beam solutions determined by Greens function method'
       ENDIF
       IF ( DO_PLANE_PARALLEL ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' Beam solutions in Plane-parallel approximation'
       ELSE
        WRITE(IUNIT, FMT_CHAR)
     &      ' Beam solutions in Pseudo-spherical approximation'
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
      ENDIF

C  surface input control

      WRITE(IUNIT, FMT_HEADING) 'surface input control'

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

C  Thermal  input control

      WRITE(IUNIT, FMT_HEADING) 'thermal input control'

      IF ( DO_THERMAL_EMISSION ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' Thermal Emission will be included in solution'
        WRITE(IUNIT, FMT_INTEGER)
     &  '  ** Number of thermal emission coefficients = ',
     &       N_THERMAL_COEFFS
      ELSE
        WRITE(IUNIT, FMT_CHAR)
     &      ' NO Thermal emission in the solution'
      ENDIF

      IF ( DO_SURFACE_EMISSION ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &      ' Surface Thermal Emission will be included in solution'
      ELSE
        WRITE(IUNIT, FMT_CHAR)
     &      ' NO Surface Thermal emission in the solution'
      ENDIF

C  output control

      WRITE(IUNIT, FMT_HEADING) 'Output control'

      IF ( DO_ADDITIONAL_MVOUT ) THEN
        WRITE(IUNIT, FMT_CHAR)
     &  ' Output of intensities + flux integerals & mean intensities'
      ELSE
        IF ( DO_MVOUT_ONLY ) THEN
          WRITE(IUNIT, FMT_CHAR)
     &      ' Output for flux integrals & mean intensities ONLY'
        ELSE
          WRITE(IUNIT, FMT_CHAR)
     &      ' Output for intensities ONLY'
        ENDIF
      ENDIF
c      IF ( SAVE_LAYER_MSST ) THEN
c        WRITE(IUNIT, FMT_CHAR)
c     &   ' Output of multiple-scatter layer source terms also given'
c      ENDIF

      WRITE(IUNIT, FMT_INTEGER)
     &     ' Number of Solar Zenith angles = ', NBEAMS
      
      IF ( .NOT.DO_NO_AZIMUTH ) THEN
        WRITE(IUNIT, FMT_INTEGER)
     &     ' Number of user-defined azimuth angles = ', N_USER_RELAZMS
      ENDIF

      IF ( DO_USER_TAUS ) THEN
        WRITE(IUNIT, FMT_INTEGER)
     &       ' Total number of optical depth output levels = ',
     &        N_OUT_USERTAUS
      ENDIF

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

      SUBROUTINE LIDORT_WRITESCEN ( SUNIT )

C  write to file of all geophysical LIDORT input

C  Include files
C  -------------

C  include file of parameters

      INCLUDE '../includes/LIDORT.PARS'

C  Include files of input variables

      INCLUDE '../includes/LIDORT_INPUTS.VARS'
      INCLUDE '../includes/LIDORT_BOOKKEEP.VARS'
      INCLUDE '../includes/LIDORT_SETUPS.VARS'

C  argument

      INTEGER             SUNIT

C  local variables

      INTEGER             N, L, UT, BRDF

C  heading and version number

      WRITE(SUNIT,FMT_SECTION)
     &     ' Geophysical scenario variables for run of LIDORT'
      WRITE(SUNIT,FMT_CHAR)
     &     ' LIDORT Version number = ',LIDORT_VERSION_NUMBER

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

C  First 5 and final phase function moments

      IF ( .NOT.DO_RAYLEIGH_ONLY .AND. .NOT.DO_ISOTROPIC_ONLY ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          WRITE(SUNIT, FMT_HEADING)
     &         '1-5 and final SCALED TOTAL phase fn. moments'
          WRITE(SUNIT,'(a,T53,a6,T60,a6,I3/)')
     &       'Layer | Moments 0 through 4-->',' ... ','Moment',NMOMENTS
          DO N = 1, NLAYERS
            WRITE(SUNIT,'(I3,5x,5f9.5,a5,f9.5)')
     &             N,(PHASMOMS_TOTAL(L,N),L=0,4),' ... ',
     &             PHASMOMS_TOTAL(NMOMENTS,N)
          ENDDO
        ELSE
          WRITE(SUNIT, FMT_HEADING)
     &         '1-5 and final unscaled TOTAL phase fn. moments'
          WRITE(SUNIT,'(a,T53,a6,T60,a6,I3/)')
     &       'Layer | Moments 0 through 4-->',' ... ','Moment',NMOMENTS
          DO N = 1, NLAYERS
            WRITE(SUNIT,'(I3,5x,5f9.5,a5,f9.5)')
     &             N,(PHASMOMS_TOTAL_INPUT(L,N),L=0,4),' ... ',
     &             PHASMOMS_TOTAL_INPUT(NMOMENTS,N)
          ENDDO
        ENDIF
      ENDIF

C Thermal emission inputs

      IF ( DO_THERMAL_EMISSION.OR.DO_SURFACE_EMISSION ) THEN
        WRITE(SUNIT,FMT_HEADING)'thermal/surface emission inputs'
      ENDIF

      IF ( DO_THERMAL_EMISSION ) THEN
        WRITE(SUNIT,FMT_HEADING)'thermal atmospheric input'
        WRITE(SUNIT,'(a/)')
     &       'Level | optical depth | thermal black body function '
        DO N = 0, NLAYERS
          WRITE(SUNIT,'(I3,7x,F9.5,7X,1PE15.6)')
     &      N,TAUGRID_INPUT(N),THERMAL_BB_INPUT(N)
        ENDDO
      ENDIF

      IF ( DO_SURFACE_EMISSION ) THEN
        WRITE(SUNIT,FMT_HEADING)'thermal surface input'
        WRITE(SUNIT,FMT_REAL)'Surface black body function',SURFBB
      ENDIF

C  surface property

      WRITE(SUNIT,FMT_HEADING)'Surface reflecting property'

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        WRITE(SUNIT,FMT_REAL)
     &       '(Lambertian) surface albedo is',LAMBERTIAN_ALBEDO
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
      DO N = 1, NBEAMS
        WRITE(SUNIT,'(I3,5x,1x,f9.5)') N,BEAM_SZAS(N)
      ENDDO

      WRITE(SUNIT,FMT_HEADING)
     &         'Computational (quadrature) angles in the half space'
      WRITE(SUNIT,'(a/)')
     &           'Stream  |  Angle    |   Cosine   |   Weight'
      DO N = 1, NSTREAMS
        WRITE(SUNIT,'(I3,5x,3(1x,f9.5,3x))')N,XANG(N),X(N),A(N)
      ENDDO

      WRITE(SUNIT,FMT_HEADING) 'User-defined stream angles'
      WRITE(SUNIT,'(a/)')'Stream  |  Angle    |   Cosine'
      DO N = 1, N_USER_STREAMS
        WRITE(SUNIT,'(I3,5x,2(1x,f9.5,3x))')
     *              N,USER_ANGLES_INPUT(N),USER_STREAMS(N)
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

C  Optical depth input
C  -------------------

      WRITE(SUNIT,FMT_SECTION)
     &       ' User-defined optical depths for output'
      WRITE(SUNIT,'(a/)')'Optical depth | Level/Layer of occurrence'
      UT = 0
      DO N = 1, N_OUT_USERTAUS
        IF ( OFFGRID_UTAU_OUTFLAG(N) ) THEN
          UT = UT + 1
          WRITE(SUNIT,'(1x,f9.5,2x,a9,i3)')
     &         USER_TAUS_INPUT(N),'    layer',OFFGRID_UTAU_LAYERIDX(UT)
        ELSE
          WRITE(SUNIT,'(1x,f9.5,2x,a9,i3)')
     &         USER_TAUS_INPUT(N),'    level',UTAU_LEVEL_MASK_UP(N)
        ENDIF
      ENDDO

C  Finish

      RETURN
      END

