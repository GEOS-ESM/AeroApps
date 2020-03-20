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
! #                   2.5, 2.6, 2.7                             #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
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
! #              VLIDORT_WRITERESULTS                           #
! #              VLIDORT_WRITEINPUT                             #
! #              VLIDORT_WRITESCEN                              #
! #              VLIDORT_WRITE_STD_INPUT                        #
! #              VLIDORT_WRITE_SUP_BRDF_INPUT                   #
! #              VLIDORT_WRITE_SUP_SS_INPUT                     #
! #              VLIDORT_WRITE_SUP_SLEAVE_INPUT                 #
! #                                                             #
! ###############################################################


      MODULE vlidort_writemodules

      PRIVATE
      PUBLIC :: VLIDORT_WRITERESULTS, &
                VLIDORT_WRITEINPUT, &
                VLIDORT_WRITESCEN, &
                VLIDORT_WRITE_STD_INPUT, &
                VLIDORT_WRITE_SUP_BRDF_INPUT, &
                VLIDORT_WRITE_SUP_SS_INPUT, &
                VLIDORT_WRITE_SUP_SLEAVE_INPUT

      CONTAINS

      SUBROUTINE VLIDORT_WRITERESULTS ( &
        RUN, &
        DO_FULLRAD_MODE, DO_SSCORR_NADIR, &
        DO_SSCORR_OUTGOING, DO_DOUBLE_CONVTEST, &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, &
        NSTOKES, VLIDORT_ACCURACY, &
        SZANGLES, N_USER_RELAZMS, &
        USER_RELAZMS, N_USER_LEVELS, &
        USER_LEVELS, HEIGHT_GRID, &
        DELTAU_VERT_INPUT, &
        DO_CLASSICAL_SOLUTION, DO_NO_AZIMUTH, &
        NBEAMS, N_DIRECTIONS, &
        WHICH_DIRECTIONS, N_OUT_STREAMS, &
        OUT_ANGLES, PARTLAYERS_OUTFLAG, &
        VZA_OFFSETS, &
        TAUGRID_INPUT, DO_MULTIBEAM, &
        STOKES, MEAN_STOKES, &
        FLUX_STOKES, FOURIER_SAVED )

!  include file of dimensions and numbers

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          RUN
      LOGICAL, INTENT (IN) ::          DO_FULLRAD_MODE
      LOGICAL, INTENT (IN) ::          DO_SSCORR_NADIR
      LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT (IN) ::          DO_DOUBLE_CONVTEST
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (IN) ::          DO_MVOUT_ONLY
      INTEGER, INTENT (IN) ::          NSTOKES
      DOUBLE PRECISION, INTENT (IN) :: VLIDORT_ACCURACY
      DOUBLE PRECISION, INTENT (IN) :: SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      DOUBLE PRECISION, INTENT (IN)::  USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: USER_LEVELS ( MAX_USER_LEVELS )
      DOUBLE PRECISION, INTENT (IN) :: HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::          DO_NO_AZIMUTH
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::          N_OUT_STREAMS
      DOUBLE PRECISION, INTENT (IN) :: OUT_ANGLES ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          VZA_OFFSETS &
          ( MAX_SZANGLES, MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: TAUGRID_INPUT ( 0:MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_MULTIBEAM ( MAXBEAMS, 0:MAXFOURIER )
      DOUBLE PRECISION, INTENT (IN) :: STOKES &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: MEAN_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: FLUX_STOKES &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::          FOURIER_SAVED ( MAX_SZANGLES )

!  local variables
!  ---------------

      INTEGER ::          I, UA, IB, V, UT, S, F, FMAX, NT, UTA
      INTEGER ::          LOCAL_NUSERAZMS, IDIR, WDIR

      DOUBLE PRECISION :: DT, USER_HEIGHTS(MAX_USER_LEVELS)
      DOUBLE PRECISION :: USER_OPDEPS (MAX_USER_LEVELS)

      CHARACTER (LEN=15), DIMENSION (MAXSTOKES) :: &
        STOKES_LABEL = (/ ' I component   ',' Q component   ', &
                          ' U component   ',' V component   ' /)

!  Beam Attenuation summary

      write(RUN,'(/a/a/)')'Results Output summary', &
                            '----------------------'

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        write(RUN,'(a)') &
             'Classical solution of beam particular integral'
      ELSE
        write(RUN,'(a)') &
             'Green function solution of beam particular integral'
      ENDIF

      IF ( DO_PLANE_PARALLEL ) THEN
        write(RUN,'(a)') &
             'Plane parallel approximation to beam attenuation'
      ELSE
        write(RUN,'(a)') &
             'Pseudo-spherical approximation to beam attenuation'
        write(RUN,'(a)') &
            ' * Average secant approximation was used'
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          write(RUN,'(a)') &
            ' * Refractive geometry was used (Snells Law)'
        ELSE
          write(RUN,'(a)') &
            ' * Straight line geometry was used (no refraction)'
        ENDIF
      ENDIF

      IF ( DO_FULLRAD_MODE ) THEN
        write(RUN,'(a)') &
             'Full Stokes calculation has been performed'
        IF ( DO_SSCORR_NADIR ) THEN
          write(RUN,'(a)') &
        '  --> Nakajima-Tanaka TMS single scatter correction (nadir)'
        ENDIF
        IF ( DO_SSCORR_OUTGOING ) THEN
          write(RUN,'(a)') &
        '  --> Nakajima-Tanaka TMS single scatter correction (outgoing)'
        ENDIF
        IF ( .NOT.DO_SSCORR_NADIR .AND. .NOT.DO_SSCORR_OUTGOING ) THEN
          write(RUN,'(a)') &
        '  --> No single scatter correction has been applied'
        ENDIF
      ELSE
        write(RUN,'(a)') &
        'ONLY Multiple-scatter radiance calculation has been performed'
      ENDIF

      IF ( .NOT.DO_NO_AZIMUTH ) THEN
        IF ( DO_DOUBLE_CONVTEST ) THEN
          write(RUN,'(/a)') &
         'Double convergence test was used for Fourier Azimuth Series'
        ELSE
          write(RUN,'(/a)') &
         'Single convergence test was used for Fourier Azimuth Series'
        ENDIF
        write(RUN,'(a,F10.7/)') &
         ' --> Accuracy level was pre-set at : ',VLIDORT_ACCURACY

        write(RUN,'(a,I5)')' - Number Solar zenith angles : ',NBEAMS
        write(run,'(/a/)') &
               '   SZA  | # Fourier | Fourier breakdown -->'
        FMAX = -100
        DO IB = 1, NBEAMS
          FMAX = MAX(FMAX,FOURIER_SAVED(IB))
        END DO
        DO IB = 1, NBEAMS
           write(RUN,'(1X,F7.2,4X,I3,6X,50(L1))') &
                  SZANGLES(IB),FOURIER_SAVED(IB), &
                  (DO_MULTIBEAM(IB,F),F=0,FMAX)
        END DO
      ELSE
        write(RUN,'(/A)') &
         'Azimuth independent output only (Fourier = 0)'
      ENDIF

!  Fix output

      DO UTA = 1, N_USER_LEVELS
        NT = INT(USER_LEVELS(UTA))
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          DT = USER_LEVELS(UTA) - DBLE(NT)
          USER_HEIGHTS(UTA) = HEIGHT_GRID(NT-1) * (ONE-DT) &
                            + HEIGHT_GRID(NT)   * DT
          USER_OPDEPS(UTA) = TAUGRID_INPUT(NT-1) &
                             + DELTAU_VERT_INPUT(NT)*DT
        ELSE
          USER_OPDEPS(UTA)  = TAUGRID_INPUT(NT)
          USER_HEIGHTS(UTA) = HEIGHT_GRID(NT)
        ENDIF
      ENDDO

!  control point for avoiding intensity output

      IF ( DO_MVOUT_ONLY ) GO TO 400

!  Stokes vector output
!  --------------------

!  local number of azimuths

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_NUSERAZMS = 1
      ELSE
        LOCAL_NUSERAZMS = N_USER_RELAZMS
      ENDIF

!  overall header

      write(RUN,'(/a/a/)')'Stokes vector output', &
                            '--------------------'

!  start beam loop

      DO IB = 1, NBEAMS

       write(RUN,FMT_REAL) &
       '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)

!  start azimuth loop

       DO UA = 1, LOCAL_NUSERAZMS

!  azimuth angle header

        IF ( DO_NO_AZIMUTH ) THEN
          write(RUN,FMT_CHAR) &
             '* * Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
        ELSE
          write(RUN,FMT_REAL) &
       '* * RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
        ENDIF

        WRITE(RUN,'(a)')' '
        WRITE(RUN,FMT_INTEGER) &
            'Total number of output levels = ',N_USER_LEVELS
        WRITE(RUN,FMT_INTEGER) &
            'Total number of output angles = ',N_OUT_STREAMS

!  detailed output

        DO IDIR = 1, N_DIRECTIONS

          WDIR = WHICH_DIRECTIONS(IDIR)

!  direction header

          IF (WDIR .EQ. UPIDX ) THEN
            WRITE(RUN,'(/A)') &
           '  --> Upwelling intensities all output levels and angles'
          ELSE IF (WDIR .EQ. DNIDX ) THEN
            WRITE(RUN,'(/A)') &
           '  --> Downwelling intensities all output levels and angles'
          ENDIF

!  output loop

          DO UT = 1, N_USER_LEVELS
            WRITE(RUN,'(/a,3(f10.5,2x))') &
              'Layer #, Height and Optical depth', &
             USER_LEVELS(UT),USER_HEIGHTS(UT), USER_OPDEPS(UT)
            write(RUN,'(a,4(a15))') &
               'output angle |  ',(STOKES_LABEL(S),S=1,NSTOKES)
            DO I = 1, N_OUT_STREAMS
              V = VZA_OFFSETS(IB,I) + UA
              WRITE(RUN,'(2x,F9.5,2x,5(1PE15.5))')OUT_ANGLES(I), &
                    (STOKES(UT,V,S,WDIR),S=1,NSTOKES)
            ENDDO
          ENDDO

!  direction and solar/azimuth loops - end

        ENDDO
       ENDDO
      ENDDO

!  integrated value output
!  -----------------------

400   CONTINUE
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

       write(RUN,'(/a/a)')'integrated value output', &
                             '-----------------------'

       DO IB = 1, NBEAMS

        write(RUN,FMT_REAL) &
       '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
        write(RUN,'(a)')' '

!  detailed output

        DO IDIR = 1, N_DIRECTIONS

          WDIR = WHICH_DIRECTIONS(IDIR)

!  direction header

          IF (WDIR .EQ. UPIDX ) THEN
            WRITE(RUN,'(/A/)') &
      '  --> Upwelling mean values & fluxes, all output levels'
          ELSE IF (WDIR .EQ. DNIDX ) THEN
            WRITE(RUN,'(/A/)') &
      '  --> Downwelling mean values & fluxes, all output levels'
          ENDIF

!  Mean values

          write(RUN,'(a/)') &
            '      ** Mean-values for Stokes components  ----> '
          write(RUN,'(a13,3x,4(a15)/)')'output levels', &
                    (STOKES_LABEL(S),S=1,NSTOKES)
          DO UT = 1, N_USER_LEVELS
            WRITE(RUN,'(2x,F9.5,3x,1p4e15.5)')USER_LEVELS(UT), &
                       (MEAN_STOKES(UT,IB,S,WDIR),S=1,NSTOKES)
          ENDDO

!  fluxes

          write(RUN,'(/a/)') &
            '      ** Fluxes for Stokes components  ----> '
          write(RUN,'(a13,3x,4(a15)/)')'output levels', &
                    (STOKES_LABEL(S),S=1,NSTOKES)
          DO UT = 1, N_USER_LEVELS
            WRITE(RUN,'(2x,F9.5,3x,1p4e15.5)')USER_LEVELS(UT), &
                       (FLUX_STOKES(UT,IB,S,WDIR),S=1,NSTOKES)
          ENDDO

!  end direction loop

        ENDDO

       ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_WRITERESULTS

!

      SUBROUTINE VLIDORT_WRITEINPUT ( &
        IUNIT, &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        DO_RAYLEIGH_ONLY, DO_QUAD_OUTPUT, &
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, &
        NSTREAMS, NLAYERS, &
        NGREEK_MOMENTS_INPUT, VLIDORT_ACCURACY, &
        FLUX_FACTOR, N_USER_RELAZMS, &
        N_USER_LEVELS, &
        DO_LAMBERTIAN_SURFACE, &
        DO_THERMAL_EMISSION, &
        DO_SURFACE_EMISSION, &
        DO_DIRECT_BEAM, DO_CLASSICAL_SOLUTION, &
        DO_NO_AZIMUTH, NBEAMS, &
        N_USER_STREAMS, DO_USER_STREAMS )

!  Write to file of all control input

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          IUNIT
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT
      LOGICAL, INTENT (IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (IN) ::          DO_MVOUT_ONLY
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
      DOUBLE PRECISION, INTENT (IN) :: VLIDORT_ACCURACY
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_THERMAL_EMISSION
      LOGICAL, INTENT (IN) ::          DO_SURFACE_EMISSION
      LOGICAL, INTENT (IN) ::          DO_DIRECT_BEAM
      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION
      LOGICAL, INTENT (IN) ::          DO_NO_AZIMUTH
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS

!  heading and version number

      WRITE(IUNIT, FMT_SECTION) &
           ' Control input variables for run of VLIDORT'
      WRITE(IUNIT, FMT_CHAR) &
           ' VLIDORT Version number = ',VLIDORT_VERSION_NUMBER

!  general control

      WRITE(IUNIT, FMT_HEADING) ' control integers'

      WRITE(IUNIT, FMT_INTEGER) &
         ' Double Gaussian quadrature over [-1,0] & [0,1]'
      WRITE(IUNIT, FMT_INTEGER) &
           '  ** Number of half space streams = ', NSTREAMS
      WRITE(IUNIT, FMT_INTEGER) &
           '  ** Number of atmospheric layers = ', NLAYERS
      WRITE(IUNIT, FMT_INTEGER) &
         '  ** Number of Greek moments (input) = ', &
                   NGREEK_MOMENTS_INPUT

!      IF ( DO_THERMAL_EMISSION ) THEN
!        WRITE(IUNIT, FMT_INTEGER)
!     &  '  ** Number of thermal emission coefficients = ',
!     &       N_THERMAL_COEFFS
!      ENDIF

      WRITE(IUNIT, FMT_HEADING) ' flux/accuracy control'
      WRITE(IUNIT, FMT_REAL) &
           ' Flux constant = ', FLUX_FACTOR

!mick thermal fix - DO_NO_AZIMUTH not defined at this point
!      IF ( .NOT.DO_NO_AZIMUTH ) THEN
!        WRITE(IUNIT, FMT_REAL) &
!       ' accuracy criterion (Fourier convergence) = ',VLIDORT_ACCURACY
!      ELSE
!        WRITE(IUNIT, FMT_CHAR) &
!        '  ** No Fourier series -- Azimuth-independent term only'
!      ENDIF

!  RTE control

      WRITE(IUNIT, FMT_HEADING) ' RTE solution control'

      IF ( DO_DIRECT_BEAM ) THEN
        WRITE(IUNIT, FMT_CHAR) &
            ' Direct beam will be included in solution'
      ELSE
        WRITE(IUNIT, FMT_CHAR) &
            ' Direct beam will NOT be included in solution'
      ENDIF

!      IF ( DO_THERMAL_EMISSION ) THEN
!        WRITE(IUNIT, FMT_CHAR)
!     &      ' Thermal Emission will be included in solution'
!      ELSE
!        WRITE(IUNIT, FMT_CHAR)
!     &      ' NO Thermal emission in the solution'
!      ENDIF

      IF ( DO_SURFACE_EMISSION ) THEN
        WRITE(IUNIT, FMT_CHAR) &
            ' Surface Thermal Emission will be included in solution'
      ELSE
        WRITE(IUNIT, FMT_CHAR) &
            ' NO Surface Thermal emission in the solution'
      ENDIF

!  surface input write

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        WRITE(IUNIT, FMT_CHAR) &
            ' Surface will be treated as Lambertian'
      ELSE
        WRITE(IUNIT, FMT_CHAR) &
          ' Supplementary Bidirectional surface input will be used'
      ENDIF

!  Other writes. Isotropic-only option now removed. 01/17/06.

!      IF ( DO_ISOTROPIC_ONLY ) THEN
!        WRITE(IUNIT, FMT_CHAR)' Medium is isotropic'
!      ENDIF

      IF ( DO_RAYLEIGH_ONLY ) THEN
        WRITE(IUNIT, FMT_CHAR) &
            ' Medium has Rayleigh scattering only'
      ELSE
        WRITE(IUNIT, FMT_CHAR) &
            ' Medium has general scattering phase function'
      ENDIF

!      IF ( DO_TRANSMITTANCE_ONLY ) THEN
!        WRITE(IUNIT, FMT_CHAR)
!     &      ' RTE solution is Transmission-only (no scattering)'
!      ELSE
!        WRITE(IUNIT, FMT_CHAR)
!     &      ' RTE solution is full multiple-scattering'
!      ENDIF

!  Beam particular integral control

      WRITE(IUNIT, FMT_HEADING) 'Beam solution control'

      IF ( DO_CLASSICAL_SOLUTION ) THEN
        WRITE(IUNIT, FMT_CHAR) &
        ' Beam solution determined by classical (Chandrasekhar) method'
      ELSE
        WRITE(IUNIT, FMT_CHAR) &
           ' Beam solution determined by Greens function method'
      ENDIF

      IF ( DO_PLANE_PARALLEL ) THEN
        WRITE(IUNIT, FMT_CHAR) &
            ' Beam solution in Plane-parallel approximation'
        ELSE
        WRITE(IUNIT, FMT_CHAR) &
            ' Beam solution in Pseudo-spherical approximation'
        WRITE(IUNIT, FMT_CHAR) &
        '  ** attenuation will use average secant estimation'
        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          WRITE(IUNIT, FMT_CHAR) &
        '  ** Refractive geometry'
        ELSE
          WRITE(IUNIT, FMT_CHAR) &
        '  ** Straight line geometry'
        ENDIF
      ENDIF

!  output control

      WRITE(IUNIT, FMT_HEADING) 'Output control'

      IF ( DO_ADDITIONAL_MVOUT ) THEN
        WRITE(IUNIT, FMT_CHAR) &
            ' Output of intensities AND fluxes & mean intensities'
      ELSE
        IF ( DO_MVOUT_ONLY ) THEN
          WRITE(IUNIT, FMT_CHAR) &
            ' Output for fluxes & mean intensities ONLY'
        ELSE
          WRITE(IUNIT, FMT_CHAR) &
            ' Output for intensities ONLY'
        ENDIF
      ENDIF

      WRITE(IUNIT, FMT_INTEGER) &
           ' Number of Solar Zenith angles = ', NBEAMS

!mick thermal fix - DO_NO_AZIMUTH not defined at this point
!      IF ( .NOT.DO_NO_AZIMUTH ) THEN
!        WRITE(IUNIT, FMT_INTEGER) &
!           ' Number of user-defined azimuth angles = ', N_USER_RELAZMS
!      ENDIF

      WRITE(IUNIT, FMT_INTEGER) &
             ' Total number of output levels = ', &
              N_USER_LEVELS

      IF ( DO_USER_STREAMS ) THEN
        IF ( DO_QUAD_OUTPUT ) THEN
          WRITE(IUNIT, FMT_CHAR) &
             ' Stream angle output will include quadratures'
          WRITE(IUNIT, FMT_INTEGER) &
             '  ** Total Number of output stream angles = ', &
                          NSTREAMS + N_USER_STREAMS
        ELSE
          WRITE(IUNIT, FMT_CHAR) &
             ' Stream angle output for user-defined angles only'
          WRITE(IUNIT, FMT_INTEGER) &
             '  ** Total Number of output stream angles = ', &
                          N_USER_STREAMS
        ENDIF
      ELSE
        IF ( .NOT. DO_MVOUT_ONLY ) THEN
          WRITE(IUNIT, FMT_CHAR) &
            ' Stream output at Quadrature angles only'
          WRITE(IUNIT, FMT_INTEGER) &
             '  ** Total Number of output stream angles = ',NSTREAMS
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_WRITEINPUT

!

      SUBROUTINE VLIDORT_WRITESCEN ( &
        SUNIT, &
        DO_DELTAM_SCALING, NSTREAMS, &
        NLAYERS, NGREEK_MOMENTS_INPUT, &
        N_SZANGLES, SZANGLES, &
        N_USER_RELAZMS, USER_RELAZMS, &
        N_USER_VZANGLES, USER_VZANGLES, &
        N_USER_LEVELS, USER_LEVELS, &
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        DO_THERMAL_EMISSION, N_THERMAL_COEFFS, &
        DO_SURFACE_EMISSION, SURFBB, &
        QUAD_STREAMS, QUAD_WEIGHTS, &
        QUAD_ANGLES, DO_NO_AZIMUTH, &
        NMOMENTS, GREEKMAT_INDEX, &
        TAUGRID_INPUT, OMEGA_TOTAL, &
        GREEKMAT_TOTAL, TAUGRID )

!  write to file of all geophysical VLIDORT input

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          SUNIT
      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
      INTEGER, INTENT (IN) ::          N_SZANGLES
      DOUBLE PRECISION, INTENT (IN) :: SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      DOUBLE PRECISION, INTENT (IN) :: USER_VZANGLES ( MAX_USER_VZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: USER_LEVELS ( MAX_USER_LEVELS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      LOGICAL, INTENT (IN) ::          DO_THERMAL_EMISSION
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS
      LOGICAL, INTENT (IN) ::          DO_SURFACE_EMISSION
      DOUBLE PRECISION, INTENT (IN) :: SURFBB
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_ANGLES  ( MAXSTREAMS )
      LOGICAL, INTENT (IN) ::          DO_NO_AZIMUTH
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          GREEKMAT_INDEX ( 6 )
      DOUBLE PRECISION, INTENT (IN) :: TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) :: TAUGRID ( 0:MAXLAYERS )

!  local variables

      INTEGER ::            S, N, L

!  heading and version number

      WRITE(SUNIT,FMT_SECTION) &
           ' Geophysical scenario variables for run of VLIDORT'
      WRITE(SUNIT,FMT_CHAR) &
           ' VLIDORT Version number = ',VLIDORT_VERSION_NUMBER

!   basic input
!   -----------

      WRITE(SUNIT,FMT_SECTION) &
           ' Basic atmospheric input for stokes vector calculation'

!  Layer optical depths and single scatter albedos

      IF ( DO_DELTAM_SCALING ) THEN
        WRITE(SUNIT, FMT_HEADING) &
            '  scaled and unscaled inputs with Delta-M turned ON'
        WRITE(SUNIT, FMT_HEADING) &
            'Layer optical depths and single-scatter-albedos '
        WRITE(SUNIT,'(a,a25/a,a25)') &
             '      (unscaled) (scaled) ','  (unscaled)   (scaled)  ', &
             'Layer  Op-depth  Op-depth ','  s.s albedo  s.s albedo '
        WRITE(SUNIT,'(a)')' '
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(I3,T6,2(f9.5,1x),2x,2(f9.5,3x))') &
                  N,TAUGRID_INPUT(N),TAUGRID(N), &
                  OMEGA_TOTAL_INPUT(N),OMEGA_TOTAL(N)
        ENDDO
      ELSE
        WRITE(SUNIT, FMT_HEADING) &
            '  Unscaled inputs only: Delta-M turned OFF'
        WRITE(SUNIT, FMT_HEADING) &
          'Layer optical depths and single-scatter-albedos '
        WRITE(SUNIT,'(a,a13/a,a13)') &
              '      (unscaled) ','  (unscaled) ', &
              'Layer  Op-depth  ','  s.s albedo '
        WRITE(SUNIT,'(a)')' '
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(I3,T6,f9.5,4x,f9.5)') &
                  N,TAUGRID_INPUT(N),OMEGA_TOTAL_INPUT(N)
        ENDDO
      ENDIF

!   phase matrix elements

      IF ( DO_DELTAM_SCALING ) THEN
        WRITE(SUNIT, FMT_HEADING) &
             'Phase matrix elements (Greek constants), SCALED'
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(a, I3/)') &
                    'Matrix elements (A-F) for layer ',N
          DO L = 0, NMOMENTS
            WRITE(SUNIT,'(I2,1x,1p6e12.4)') &
           L,(GREEKMAT_TOTAL(L,N,GREEKMAT_INDEX(S)),S = 1, 6)
          ENDDO
        ENDDO
      ELSE
        WRITE(SUNIT, FMT_HEADING) &
             'Phase matrix elements (Greek constants), UNSCALED'
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(a, I3/)') &
                    'Matrix elements (A-F) for layer ',N
          DO L = 0, NGREEK_MOMENTS_INPUT
            WRITE(SUNIT,'(I2,1x,1p6e12.4)') &
          L,(GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(S)),S = 1, 6)
          ENDDO
        ENDDO
      ENDIF

!  Commented out thermal expansion coefficients
!      IF ( DO_THERMAL_EMISSION ) THEN
!        WRITE(SUNIT,FMT_HEADING)'thermal emission coefficients'
!        WRITE(SUNIT,'(a/)')
!     &       'Layer | thermal emission expansion coeffs-->'
!        DO N = 1, NLAYERS
!          WRITE(SUNIT,'(I3,4x,10(f10.5))')
!     &          N,(THERMAL_COEFFS(N,S),S=1,N_THERMAL_COEFFS)
!        ENDDO
!      ENDIF

!  surface property

      WRITE(SUNIT,FMT_HEADING)'Surface reflecting property'

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        WRITE(SUNIT,FMT_REAL) &
             '(Lambertian) surface albedo is',LAMBERTIAN_ALBEDO
        IF ( DO_SURFACE_EMISSION ) THEN
          WRITE(SUNIT,FMT_REAL) &
            '(Lambertian) emissivity = ',ONE-LAMBERTIAN_ALBEDO
          WRITE(SUNIT,FMT_REAL) &
            'Surface blackbody function is', SURFBB
        ENDIF
      ELSE
        WRITE(SUNIT,FMT_CHAR) &
               ' BRDF Supplementary input'
      ENDIF

!   Geometry input
!   --------------

      WRITE(SUNIT,FMT_SECTION)' Geometry input'

      WRITE(SUNIT,FMT_HEADING) &
                'Solar zenith angles'
      WRITE(SUNIT,'(a/)')'Number |   Angle'
      DO N = 1, N_SZANGLES
        WRITE(SUNIT,'(I3,5x,1x,f9.5)') N,SZANGLES(N)
      ENDDO

      WRITE(SUNIT,FMT_HEADING) &
               'Computational (quadrature) angles in the half space'
      WRITE(SUNIT,'(a/)') &
                 'Stream  |  Angle    |   Cosine   |   Weight'
      DO N = 1, NSTREAMS
        WRITE(SUNIT,'(I3,5x,3(1x,f9.5,3x))')N,QUAD_ANGLES(N), &
              QUAD_STREAMS(N),QUAD_WEIGHTS(N)
      ENDDO

      WRITE(SUNIT,FMT_HEADING) 'User-defined viewing zenith angles'
      WRITE(SUNIT,'(a/)')'Number  |  Angle    |   Cosine'
      DO N = 1, N_USER_VZANGLES
        WRITE(SUNIT,'(I3,5x,2(1x,f9.5,3x))') &
        N,USER_VZANGLES(N),DCOS(USER_VZANGLES(N)*DEG_TO_RAD)
      ENDDO

      IF ( DO_NO_AZIMUTH ) THEN
        WRITE(SUNIT,FMT_HEADING) 'No azimuth angles'
      ELSE
        WRITE(SUNIT,FMT_HEADING) &
                'User-defined relative azimuth angles'
        WRITE(SUNIT,'(a/)')'Number |   Angle'
        DO N = 1, N_USER_RELAZMS
          WRITE(SUNIT,'(I3,5x,1x,f9.5)') N,USER_RELAZMS(N)
        ENDDO
      ENDIF

!  Output levels
!  -------------

      WRITE(SUNIT,FMT_SECTION) &
             ' User-defined levels for output'
      WRITE(SUNIT,'(a/)')' # | Level/Layer of occurrence'
      DO N = 1, N_USER_LEVELS
        WRITE(SUNIT,'(i3,4x,f9.5)')N, USER_LEVELS(N)
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_WRITESCEN

!

      SUBROUTINE VLIDORT_WRITE_STD_INPUT ( &
        DO_FULLRAD_MODE,DO_SSCORR_TRUNCATION,DO_SS_EXTERNAL,DO_SSFULL,&
        DO_THERMAL_EMISSION,DO_SURFACE_EMISSION,DO_PLANE_PARALLEL,&
        DO_UPWELLING,DO_DNWELLING,DO_QUAD_OUTPUT,&
        DO_TOA_CONTRIBS,DO_LAMBERTIAN_SURFACE,&
        DO_SPECIALIST_OPTION_1,DO_SPECIALIST_OPTION_2,DO_SPECIALIST_OPTION_3,&
        DO_SURFACE_LEAVING,DO_SL_ISOTROPIC,&
        NSTOKES,NSTREAMS,NLAYERS,&
        NFINELAYERS,N_THERMAL_COEFFS,VLIDORT_ACCURACY,&
        NLAYERS_NOMS,NLAYERS_CUTOFF,FLUX_FACTOR,N_USER_LEVELS,&
        HEIGHT_GRID,PRESSURE_GRID,TEMPERATURE_GRID,&
        FINEGRID,RFINDEX_PARAMETER,&
        DELTAU_VERT_INPUT,GREEKMAT_TOTAL_INPUT,THERMAL_BB_INPUT,&
        LAMBERTIAN_ALBEDO,SURFBB,&  ! ********
        DO_DEBUG_WRITE,DO_WRITE_INPUT,DO_WRITE_SCENARIO,&
        DO_WRITE_FOURIER,DO_WRITE_RESULTS,&
        INPUT_WRITE_FILENAME,SCENARIO_WRITE_FILENAME,&
        FOURIER_WRITE_FILENAME,RESULTS_WRITE_FILENAME,&
        DO_SSCORR_NADIR,DO_SSCORR_OUTGOING,DO_FO_CALC,DO_DOUBLE_CONVTEST,&
        DO_SOLAR_SOURCES,DO_REFRACTIVE_GEOMETRY,DO_CHAPMAN_FUNCTION,&
        DO_RAYLEIGH_ONLY,DO_DELTAM_SCALING,DO_SOLUTION_SAVING,&
        DO_BVP_TELESCOPING,DO_USER_VZANGLES,DO_ADDITIONAL_MVOUT,&
        DO_MVOUT_ONLY,DO_THERMAL_TRANSONLY,DO_OBSERVATION_GEOMETRY,&
        NGREEK_MOMENTS_INPUT,N_SZANGLES,SZANGLES,&
        N_USER_RELAZMS,USER_RELAZMS,N_USER_VZANGLES,USER_VZANGLES,&
        USER_LEVELS,GEOMETRY_SPECHEIGHT,N_USER_OBSGEOMS,USER_OBSGEOMS,&
        EARTH_RADIUS,OMEGA_TOTAL_INPUT)

!  3/28/14. Changes for Version 2.7. remove LTE linearization references

      USE VLIDORT_PARS

      IMPLICIT NONE

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

      LOGICAL, INTENT(IN) ::            DO_FULLRAD_MODE
      LOGICAL, INTENT(IN) ::            DO_SSCORR_TRUNCATION
      LOGICAL, INTENT(IN) ::            DO_SS_EXTERNAL
      LOGICAL, INTENT(IN) ::            DO_SSFULL
      LOGICAL, INTENT(IN) ::            DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN) ::            DO_SURFACE_EMISSION
      LOGICAL, INTENT(IN) ::            DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN) ::            DO_UPWELLING
      LOGICAL, INTENT(IN) ::            DO_DNWELLING
      LOGICAL, INTENT(IN) ::            DO_QUAD_OUTPUT
      LOGICAL, INTENT(IN) ::            DO_TOA_CONTRIBS
      LOGICAL, INTENT(IN) ::            DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT(IN) ::            DO_SPECIALIST_OPTION_1
      LOGICAL, INTENT(IN) ::            DO_SPECIALIST_OPTION_2
      LOGICAL, INTENT(IN) ::            DO_SPECIALIST_OPTION_3
      LOGICAL, INTENT(IN) ::            DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN) ::            DO_SL_ISOTROPIC

      INTEGER, INTENT(IN) ::            NSTOKES
      INTEGER, INTENT(IN) ::            NSTREAMS
      INTEGER, INTENT(IN) ::            NLAYERS
      INTEGER, INTENT(IN) ::            NFINELAYERS
      INTEGER, INTENT(IN) ::            N_THERMAL_COEFFS
      DOUBLE PRECISION, INTENT(IN) ::   VLIDORT_ACCURACY
      INTEGER, INTENT(IN) ::            NLAYERS_NOMS
      INTEGER, INTENT(IN) ::            NLAYERS_CUTOFF

      DOUBLE PRECISION, INTENT(IN) ::   FLUX_FACTOR

      INTEGER, INTENT(IN) ::            N_USER_LEVELS

      DOUBLE PRECISION, INTENT(IN) ::   HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   PRESSURE_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   TEMPERATURE_GRID ( 0:MAXLAYERS )
      INTEGER, INTENT(IN) ::            FINEGRID ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   RFINDEX_PARAMETER

      DOUBLE PRECISION, INTENT(IN) ::   DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   GREEKMAT_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT(IN) ::   THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT(IN) ::   SURFBB

!  Code superseded, Version 2.7. Remove this
!      DOUBLE PRECISION, INTENT(IN) ::   LTE_DELTAU_VERT_INPUT ( 2, MAXLAYERS )
!      DOUBLE PRECISION, INTENT(IN) ::   LTE_THERMAL_BB_INPUT ( 0:MAXLAYERS )

      LOGICAL, INTENT(IN) ::            DO_DEBUG_WRITE
      LOGICAL, INTENT(IN) ::            DO_WRITE_INPUT
      LOGICAL, INTENT(IN) ::            DO_WRITE_SCENARIO
      LOGICAL, INTENT(IN) ::            DO_WRITE_FOURIER
      LOGICAL, INTENT(IN) ::            DO_WRITE_RESULTS
      CHARACTER (LEN=60), INTENT(IN) :: INPUT_WRITE_FILENAME
      CHARACTER (LEN=60), INTENT(IN) :: SCENARIO_WRITE_FILENAME
      CHARACTER (LEN=60), INTENT(IN) :: FOURIER_WRITE_FILENAME
      CHARACTER (LEN=60), INTENT(IN) :: RESULTS_WRITE_FILENAME

!  --------------------------
!  Standard Inputs - Variable
!  --------------------------

      LOGICAL, INTENT(IN) ::          DO_SSCORR_NADIR
      LOGICAL, INTENT(IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT(IN) ::          DO_FO_CALC
      LOGICAL, INTENT(IN) ::          DO_DOUBLE_CONVTEST
      LOGICAL, INTENT(IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN) ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT(IN) ::          DO_CHAPMAN_FUNCTION
      LOGICAL, INTENT(IN) ::          DO_RAYLEIGH_ONLY
      LOGICAL, INTENT(IN) ::          DO_DELTAM_SCALING
      LOGICAL, INTENT(IN) ::          DO_SOLUTION_SAVING
      LOGICAL, INTENT(IN) ::          DO_BVP_TELESCOPING
      LOGICAL, INTENT(IN) ::          DO_USER_VZANGLES
      LOGICAL, INTENT(IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT(IN) ::          DO_MVOUT_ONLY
      LOGICAL, INTENT(IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT(IN) ::          DO_OBSERVATION_GEOMETRY

      INTEGER, INTENT(IN) ::          NGREEK_MOMENTS_INPUT

      INTEGER, INTENT(IN) ::          N_SZANGLES
      DOUBLE PRECISION, INTENT(IN) :: SZANGLES ( MAX_SZANGLES )

      INTEGER, INTENT(IN) ::          N_USER_RELAZMS
      DOUBLE PRECISION, INTENT(IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER, INTENT(IN) ::          N_USER_VZANGLES
      DOUBLE PRECISION, INTENT(IN) :: USER_VZANGLES ( MAX_USER_VZANGLES )
      DOUBLE PRECISION, INTENT(IN) :: USER_LEVELS ( MAX_USER_LEVELS )
      DOUBLE PRECISION, INTENT(IN) :: GEOMETRY_SPECHEIGHT
      INTEGER, INTENT(IN) ::          N_USER_OBSGEOMS
      DOUBLE PRECISION, INTENT(IN) :: USER_OBSGEOMS ( MAX_USER_OBSGEOMS, 3 )

      DOUBLE PRECISION, INTENT(IN) :: EARTH_RADIUS

      DOUBLE PRECISION, INTENT(IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SZA,ULEV,LAY,SS,LAYI,LAYJ,MOM,NGREEK,K,STRM,S,USTRM,I,IB,URA,&
                 STRMI,STRMJ,UVA,UOG, GMASK(8)
      data GMASK / 1, 2, 5, 6, 11, 12, 15, 16 /

!  Open output file

      OUTUNIT = 101
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_STD_INPUT.dbg',&
            status = 'replace')

!  Define local variable.
!  Changed 16 oct 14, only output non-zero entries of GreekMat

      NGREEK = 1
      IF ( NSTOKES.eq.3 ) NGREEK = 5
      IF ( NSTOKES.eq.4 ) NGREEK = 8

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Fixed'
      WRITE(OUTUNIT,'(A)') '-----------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_FULLRAD_MODE        = ',DO_FULLRAD_MODE
      WRITE(OUTUNIT,DWFL)  'DO_SSCORR_TRUNCATION   = ',DO_SSCORR_TRUNCATION
      WRITE(OUTUNIT,DWFL)  'DO_SS_EXTERNAL         = ',DO_SS_EXTERNAL
      WRITE(OUTUNIT,DWFL)  'DO_SSFULL              = ',DO_SSFULL
      WRITE(OUTUNIT,DWFL)  'DO_THERMAL_EMISSION    = ',DO_THERMAL_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_EMISSION    = ',DO_SURFACE_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_PLANE_PARALLEL      = ',DO_PLANE_PARALLEL
      WRITE(OUTUNIT,DWFL)  'DO_UPWELLING           = ',DO_UPWELLING
      WRITE(OUTUNIT,DWFL)  'DO_DNWELLING           = ',DO_DNWELLING
      WRITE(OUTUNIT,DWFL)  'DO_QUAD_OUTPUT         = ',DO_QUAD_OUTPUT
      WRITE(OUTUNIT,DWFL)  'DO_TOA_CONTRIBS        = ',DO_TOA_CONTRIBS
      WRITE(OUTUNIT,DWFL)  'DO_LAMBERTIAN_SURFACE  = ',DO_LAMBERTIAN_SURFACE
      WRITE(OUTUNIT,DWFL)  'DO_SPECIALIST_OPTION_1 = ',DO_SPECIALIST_OPTION_1
      WRITE(OUTUNIT,DWFL)  'DO_SPECIALIST_OPTION_2 = ',DO_SPECIALIST_OPTION_2
      WRITE(OUTUNIT,DWFL)  'DO_SPECIALIST_OPTION_3 = ',DO_SPECIALIST_OPTION_3
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_LEAVING     = ',DO_SURFACE_LEAVING
      WRITE(OUTUNIT,DWFL)  'DO_SL_ISOTROPIC        = ',DO_SL_ISOTROPIC

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NSTOKES          = ',NSTOKES
      WRITE(OUTUNIT,DWFI)  'NSTREAMS         = ',NSTREAMS
      WRITE(OUTUNIT,DWFI)  'NLAYERS          = ',NLAYERS
      WRITE(OUTUNIT,DWFI)  'NFINELAYERS      = ',NFINELAYERS
      WRITE(OUTUNIT,DWFI)  'N_THERMAL_COEFFS = ',N_THERMAL_COEFFS
      WRITE(OUTUNIT,DWFR)  'VLIDORT_ACCURACY = ',VLIDORT_ACCURACY
      WRITE(OUTUNIT,DWFI)  'NLAYERS_NOMS     = ',NLAYERS_NOMS
      WRITE(OUTUNIT,DWFI)  'NLAYERS_CUTOFF   = ',NLAYERS_CUTOFF

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'FLUX_FACTOR      = ',FLUX_FACTOR

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_LEVELS    = ',N_USER_LEVELS

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' HEIGHT_GRID(LAY)      = ',HEIGHT_GRID(LAY)
      END DO
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' PRESSURE_GRID(LAY)    = ',PRESSURE_GRID(LAY)
      END DO
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' TEMPERATURE_GRID(LAY) = ',TEMPERATURE_GRID(LAY)
      END DO
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFI1)  'LAY = ',LAY,&
          ' FINEGRID (LAY)        = ',FINEGRID(LAY)
      END DO
      WRITE(OUTUNIT,DWFR)  'RFINDEX_PARAMETER = ',RFINDEX_PARAMETER

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' DELTAU_VERT_INPUT(LAY) = ',DELTAU_VERT_INPUT(LAY)
      END DO
      DO K=1,NGREEK
        SS = GMASK(K)
        DO LAY=1,NLAYERS
          DO MOM=0,NGREEK_MOMENTS_INPUT
            WRITE(OUTUNIT,DWFR3)  'SS = ',SS,' LAY = ',LAY,' MOM = ',MOM,&
              ' GREEKMAT_TOTAL_INPUT(MOM,LAY,SS) = ',&
                GREEKMAT_TOTAL_INPUT(MOM,LAY,SS)
          END DO
        END DO
      END DO
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' THERMAL_BB_INPUT(LAY) = ',THERMAL_BB_INPUT(LAY)
      END DO
      WRITE(OUTUNIT,DWFR)  'LAMBERTIAN_ALBEDO = ',LAMBERTIAN_ALBEDO
      WRITE(OUTUNIT,DWFR)  'SURFBB            = ',SURFBB

!  Next two loops removed, Version 2.7
 !     DO LAY=1,NLAYERS
 !       DO I=1,2
 !         WRITE(OUTUNIT,DWFR2)  'LAY = ',LAY,' I = ',I,&
 !           ' LTE_DELTAU_VERT_INPUT(I,LAY) = ',LTE_DELTAU_VERT_INPUT(I,LAY)
 !       END DO
 !     END DO
 !     DO LAY=0,NLAYERS
 !       WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
 !         ' LTE_THERMAL_BB_INPUT(LAY) = ',LTE_THERMAL_BB_INPUT(LAY)
 !     END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_DEBUG_WRITE          = ',DO_DEBUG_WRITE
      WRITE(OUTUNIT,DWFL)  'DO_WRITE_INPUT          = ',DO_WRITE_INPUT
      WRITE(OUTUNIT,DWFL)  'DO_WRITE_SCENARIO       = ',DO_WRITE_SCENARIO
      WRITE(OUTUNIT,DWFL)  'DO_WRITE_FOURIER        = ',DO_WRITE_FOURIER
      WRITE(OUTUNIT,DWFL)  'DO_WRITE_RESULTS        = ',DO_WRITE_RESULTS
      WRITE(OUTUNIT,DWFC)  'INPUT_WRITE_FILENAME    = ',INPUT_WRITE_FILENAME
      WRITE(OUTUNIT,DWFC)  'SCENARIO_WRITE_FILENAME = ',SCENARIO_WRITE_FILENAME
      WRITE(OUTUNIT,DWFC)  'FOURIER_WRITE_FILENAME  = ',FOURIER_WRITE_FILENAME
      WRITE(OUTUNIT,DWFC)  'RESULTS_WRITE_FILENAME  = ',RESULTS_WRITE_FILENAME

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '--------------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Variable'
      WRITE(OUTUNIT,'(A)') '--------------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_SSCORR_NADIR         = ',DO_SSCORR_NADIR
      WRITE(OUTUNIT,DWFL)  'DO_SSCORR_OUTGOING      = ',DO_SSCORR_OUTGOING
      WRITE(OUTUNIT,DWFL)  'DO_FO_CALC              = ',DO_FO_CALC
      WRITE(OUTUNIT,DWFL)  'DO_DOUBLE_CONVTEST      = ',DO_DOUBLE_CONVTEST
      WRITE(OUTUNIT,DWFL)  'DO_SOLAR_SOURCES        = ',DO_SOLAR_SOURCES
      WRITE(OUTUNIT,DWFL)  'DO_REFRACTIVE_GEOMETRY  = ',DO_REFRACTIVE_GEOMETRY
      WRITE(OUTUNIT,DWFL)  'DO_CHAPMAN_FUNCTION     = ',DO_CHAPMAN_FUNCTION
      WRITE(OUTUNIT,DWFL)  'DO_RAYLEIGH_ONLY        = ',DO_RAYLEIGH_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_DELTAM_SCALING       = ',DO_DELTAM_SCALING
      WRITE(OUTUNIT,DWFL)  'DO_SOLUTION_SAVING      = ',DO_SOLUTION_SAVING
      WRITE(OUTUNIT,DWFL)  'DO_BVP_TELESCOPING      = ',DO_BVP_TELESCOPING
      WRITE(OUTUNIT,DWFL)  'DO_USER_VZANGLES        = ',DO_USER_VZANGLES
      WRITE(OUTUNIT,DWFL)  'DO_ADDITIONAL_MVOUT     = ',DO_ADDITIONAL_MVOUT
      WRITE(OUTUNIT,DWFL)  'DO_MVOUT_ONLY           = ',DO_MVOUT_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_THERMAL_TRANSONLY    = ',DO_THERMAL_TRANSONLY
      WRITE(OUTUNIT,DWFL)  'DO_OBSERVATION_GEOMETRY = ',DO_OBSERVATION_GEOMETRY

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NGREEK_MOMENTS_INPUT    = ',NGREEK_MOMENTS_INPUT

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_SZANGLES     = ',N_SZANGLES
      DO SZA=1,N_SZANGLES
        WRITE(OUTUNIT,DWFR1)  'SZA = ',SZA,&
          ' SZANGLES(SZA) = ',SZANGLES(SZA)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_RELAZMS = ',N_USER_RELAZMS
      DO URA=1,N_USER_RELAZMS
        WRITE(OUTUNIT,DWFR1)  'URA = ',URA,&
          ' USER_RELAZMS(URA) = ',USER_RELAZMS(URA)
      END DO
      WRITE(OUTUNIT,DWFI)  'N_USER_VZANGLES = ',N_USER_VZANGLES
      DO UVA=1,N_USER_VZANGLES
        WRITE(OUTUNIT,DWFR1)  'UVA = ',UVA,&
          ' USER_VZANGLES(UVA) = ',USER_VZANGLES(UVA)
      END DO
      DO ULEV=1,N_USER_LEVELS
        WRITE(OUTUNIT,DWFR1)  'ULEV = ',ULEV,&
          ' USER_LEVELS(ULEV)  = ',USER_LEVELS(ULEV)
      END DO
      WRITE(OUTUNIT,DWFR)  'GEOMETRY_SPECHEIGHT = ',GEOMETRY_SPECHEIGHT
      WRITE(OUTUNIT,DWFI)  'N_USER_OBSGEOMS     = ',N_USER_OBSGEOMS
      DO UOG=1,N_USER_OBSGEOMS
        WRITE(OUTUNIT,DWFR1_3)  'UOG = ',UOG,&
          ' USER_OBSGEOMS(UOG,1:3) = ',USER_OBSGEOMS(UOG,1),&
                                   ',',USER_OBSGEOMS(UOG,2),&
                                   ',',USER_OBSGEOMS(UOG,3)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'EARTH_RADIUS        = ',EARTH_RADIUS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' OMEGA_TOTAL_INPUT(LAY) = ',OMEGA_TOTAL_INPUT(LAY)
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_STD_INPUT

!

      SUBROUTINE VLIDORT_WRITE_SUP_BRDF_INPUT ( &
        NSTOKES,NSTREAMS,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,&
        EXACTDB_BRDFUNC,BRDF_F_0,BRDF_F,USER_BRDF_F_0,USER_BRDF_F,&
        EMISSIVITY,USER_EMISSIVITY)

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   N_SZANGLES
      INTEGER, INTENT(IN) ::   N_USER_VZANGLES
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS

      DOUBLE PRECISION, INTENT(IN) ::   EXACTDB_BRDFUNC &
          ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT(IN) ::   USER_BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   USER_BRDF_F &
          ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT(IN) ::   EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT(IN) ::   USER_EMISSIVITY &
          ( MAXSTOKES, MAX_USER_STREAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SZA,ULEV,LAY,SS,MOM,NREFL,K,STRM,S,USTRM,I,IB,URA,&
                 STRMI,STRMJ,UVA, RMASK3(9), RMASK4(16), RMASK(16)
      INTEGER :: NMOMENTS
      data RMASK3 / 1, 2, 3, 5, 6, 7, 9, 10, 11 /
      data RMASK4 / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 /

!  Open output file

      OUTUNIT = 102
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_SUP_BRDF_INPUT.dbg',&
            status = 'replace')

!  Define local variables
!   Changed 16 October 2014, according to correct entries in 4x4 BRDF matrices

      NREFL = NSTOKES * NSTOKES ; RMASK(1) = 1
      IF ( NSTOKES.eq.3 ) RMASK(1:NREFL) = RMASK3(1:NREFL)
      IF ( NSTOKES.eq.4 ) RMASK(1:NREFL) = RMASK4(1:NREFL)
      NMOMENTS = 2*NSTREAMS

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '----------------------'
      WRITE(OUTUNIT,'(A)') 'BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '----------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,N_SZANGLES
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_VZANGLES
            DO K=1,NREFL
              SS = RMASK(K)
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,' SS = ',SS,&
                ' EXACTDB_BRDFUNC(SS,USTRM,URA,IB) = ',&
                  EXACTDB_BRDFUNC(SS,USTRM,URA,IB)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO IB=1,N_SZANGLES
        DO STRM=1,NSTREAMS
          DO K=1,NREFL
            SS = RMASK(K)
            DO MOM=0,NMOMENTS
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' STRM = ',STRM,' SS = ',SS,' MOM = ',MOM,&
                ' BRDF_F_0(MOM,SS,STRM,IB) = ',BRDF_F_0(MOM,SS,STRM,IB)
            END DO
          END DO
        END DO
      END DO
      DO STRMJ=1,NSTREAMS
        DO STRMI=1,NSTREAMS
          DO K=1,NREFL
            SS = RMASK(K)
            DO MOM=0,NMOMENTS
              WRITE(OUTUNIT,DWFR4) &
                'STRMJ = ',STRMJ,' STRMI = ',STRMI,' SS = ',SS,' MOM = ',MOM,&
                ' BRDF_F(MOM,SS,STRMI,STRMJ) = ',BRDF_F(MOM,SS,STRMI,STRMJ)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO IB=1,N_SZANGLES
        DO USTRM=1,N_USER_VZANGLES
          DO K=1,NREFL
            SS = RMASK(K)
            DO MOM=0,NMOMENTS
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' USTRM = ',USTRM,' SS = ',SS,' MOM = ',MOM,&
                ' USER_BRDF_F_0(MOM,SS,USTRM,IB) = ',&
                  USER_BRDF_F_0(MOM,SS,USTRM,IB)
            END DO
          END DO
        END DO
      END DO
      DO STRM=1,NSTREAMS
        DO USTRM=1,N_USER_VZANGLES
          DO K=1,NREFL
            SS = RMASK(K)
            DO MOM=0,NMOMENTS
              WRITE(OUTUNIT,DWFR4) &
                'STRM = ',STRM,' USTRM = ',USTRM,' SS = ',SS,' MOM = ',MOM,&
                ' USER_BRDF_F(MOM,SS,USTRM,STRM) = ',&
                  USER_BRDF_F(MOM,SS,USTRM,STRM)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO STRM=1,NSTREAMS
        DO S=1,NSTOKES
          WRITE(OUTUNIT,DWFR2)  'STRM = ',STRM,' S = ',S,&
            ' EMISSIVITY (S,STRM) = ',EMISSIVITY (S,STRM)
        END DO
      END DO
      DO USTRM=1,N_USER_VZANGLES
        DO S=1,NSTOKES
          WRITE(OUTUNIT,DWFR2)  'USTRM = ',USTRM,' S = ',S,&
            ' USER_EMISSIVITY (S,USTRM) = ',USER_EMISSIVITY(S,USTRM)
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_SUP_BRDF_INPUT

!

      SUBROUTINE VLIDORT_WRITE_SUP_SS_INPUT ( &
        NSTOKES,N_USER_LEVELS,&
        STOKES_SS,STOKES_DB)

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   N_USER_LEVELS

      DOUBLE PRECISION, INTENT(IN) ::  STOKES_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) ::  STOKES_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: DIR,GEO,S,ULEV

!  Open output file

      OUTUNIT = 103
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_SUP_SS_INPUT.dbg',&
            status = 'replace')

!  Write all single-scatter input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '--------------------'
      WRITE(OUTUNIT,'(A)') 'SS Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '--------------------'

      WRITE(OUTUNIT,*)
      DO DIR=1,MAX_DIRECTIONS
        DO S=1,NSTOKES
          DO GEO=1,MAX_GEOMETRIES
            DO ULEV=1,N_USER_LEVELS
            WRITE(OUTUNIT,DWFR4) &
              'DIR = ',DIR,' S = ',S,' GEO = ',GEO,' ULEV = ',ULEV,&
              ' STOKES_SS(ULEV,GEO,S,DIR) = ',STOKES_SS(ULEV,GEO,S,DIR)
            END DO
          END DO
        END DO
      END DO
      DO S=1,NSTOKES
        DO GEO=1,MAX_GEOMETRIES
          DO ULEV=1,N_USER_LEVELS
          WRITE(OUTUNIT,DWFR3) &
            'S = ',S,' GEO = ',GEO,' ULEV = ',ULEV,&
            ' STOKES_DB(ULEV,GEO,S) = ',STOKES_DB(ULEV,GEO,S)
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_SUP_SS_INPUT

!

      SUBROUTINE VLIDORT_WRITE_SUP_SLEAVE_INPUT ( &
        NSTOKES,NSTREAMS,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,&
        SLTERM_ISOTROPIC,SLTERM_USERANGLES,SLTERM_F_0,USER_SLTERM_F_0)

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   N_SZANGLES
      INTEGER, INTENT(IN) ::   N_USER_VZANGLES
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS

      DOUBLE PRECISION, INTENT(IN) ::   SLTERM_ISOTROPIC &
          ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   SLTERM_USERANGLES &
          ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

      DOUBLE PRECISION, INTENT(IN) ::   SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   USER_SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SZA,ULEV,LAY,SS,MOM,STRM,S,USTRM,I,IB,URA,&
                 STRMI,STRMJ,UVA
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_SUP_SLEAVE_INPUT.dbg',&
            status = 'replace')

!  Define local variable

      NMOMENTS = 2*NSTREAMS

!  Write all surface-leaving input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '------------------------'
      WRITE(OUTUNIT,'(A)') 'SLEAVE Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '------------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,N_SZANGLES
        DO S=1,NSTOKES
          WRITE(OUTUNIT,DWFR2) &
            'IB = ',IB,' S = ',S,&
            ' SLTERM_ISOTROPIC(S,IB) = ',SLTERM_ISOTROPIC(S,IB)
        END DO
      END DO
      DO IB=1,N_SZANGLES
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_VZANGLES
            DO S=1,NSTOKES
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,' S = ',S,&
                ' SLTERM_USERANGLES(S,USTRM,URA,IB) = ',&
                  SLTERM_USERANGLES(S,USTRM,URA,IB)
            END DO
          END DO
        END DO
      END DO
      DO IB=1,N_SZANGLES
        DO STRM=1,NSTREAMS
          DO S=1,NSTOKES
            DO MOM=0,NMOMENTS
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' STRM = ',STRM,' S = ',S,' MOM = ',MOM,&
                ' SLTERM_F_0(MOM,S,STRM,IB) = ',SLTERM_F_0(MOM,S,STRM,IB)
            END DO
          END DO
        END DO
      END DO
      DO IB=1,N_SZANGLES
        DO USTRM=1,N_USER_VZANGLES
          DO S=1,NSTOKES
            DO MOM=0,NMOMENTS
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' USTRM = ',USTRM,' S = ',S,' MOM = ',MOM,&
                ' USER_SLTERM_F_0(MOM,S,USTRM,IB) = ',&
                  USER_SLTERM_F_0(MOM,S,USTRM,IB)
            END DO
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_SUP_SLEAVE_INPUT

      END MODULE vlidort_writemodules

