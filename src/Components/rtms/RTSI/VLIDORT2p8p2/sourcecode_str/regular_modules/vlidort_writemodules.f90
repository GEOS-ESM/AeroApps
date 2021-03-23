
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p2                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p2                             #
! #  Release Date :   15 April 2020                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released August 2014                        #
! #      2.8   F90, released May    2017                        #
! #      2.8.1 F90, released August 2019                        # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.2, released 15 April 2020.                            #
! #     ==> Geometry (FO/MS), check/derive separation           #
! #     ==> New setup_master for Geometry/Check/Derive          #
! #     ==> Reduction of zeroing, some dynamic memory           #
! #     ==> Use of F-matrixes only in FO code                   #
! #     ==> Use I/O type structures directly                    #
! #     ==> Doublet geometry post-processing option             #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.2 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2020.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p2 ( Version 2.8.2 )            #
! #                                                                 #
! # VLIDORT_2p8p2 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p2 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p2  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #              VLIDORT_WRITERESULTS                           #
! #              VLIDORT_WRITESCEN                              #
! #              VLIDORT_WRITE_STD_INPUT                        #
! #              VLIDORT_WRITE_SUP_BRDF_INPUT                   #
! #              VLIDORT_WRITE_SUP_SS_INPUT                     #
! #              VLIDORT_WRITE_SUP_SLEAVE_INPUT                 #
! #                                                             #
! ###############################################################

!  Upgrade, Version 2.8.1, September 2019. Write new STD input variables
!    -- (TOA/BOA illumination, Planetary, WLEAVING

!  4/15/20. Version 2.8.2. NMOMENTS instead of NGREEK_MOMENTS_INPUT, NMOMS, Remove it.
!  4/15/20. Version 2.8.2. Remove DO_SSCORR_TRUNCATION, DO_SSCORR_USEFMAT
!  4/15/20. Version 2.8.2. Remove VLIDORT_WRITE_INPUT routine (superceded)
!  4/15/20. Version 2.8.2. Add Doublet_geometry flag to VLIDORT_WRITE_STD_INPUT

      MODULE vlidort_writemodules_m

      PRIVATE
      PUBLIC :: VLIDORT_WRITERESULTS,         &
                VLIDORT_WRITESCEN,            &
                VLIDORT_WRITE_STD_INPUT,      &
                VLIDORT_WRITE_SUP_BRDF_INPUT, &
                VLIDORT_WRITE_SUP_SS_INPUT,   &
                VLIDORT_WRITE_SUP_SLEAVE_INPUT

      CONTAINS

      SUBROUTINE VLIDORT_WRITERESULTS ( &
        RUN, DO_FULLRAD_MODE, DO_FOCORR_NADIR,     &
        DO_FOCORR_OUTGOING, DO_DOUBLE_CONVTEST,    &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, & 
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,        &
        NSTOKES, VLIDORT_ACCURACY,                 &
        SZANGLES, N_USER_RELAZMS,                  &
        USER_RELAZMS, N_USER_LEVELS, &
        USER_LEVELS, HEIGHT_GRID,    &
        DELTAU_VERT_INPUT,           &
        DO_NO_AZIMUTH,               &
        NBEAMS, N_DIRECTIONS,        &
        WHICH_DIRECTIONS, N_OUT_STREAMS, &
        OUT_ANGLES, PARTLAYERS_OUTFLAG,  &
        VZA_OFFSETS,                     &
        TAUGRID_INPUT, DO_MULTIBEAM,     &
        STOKES, MEAN_STOKES,             &
        FLUX_STOKES, FOURIER_SAVED )

!  include file of dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS, MAXLAYERS, MAX_SZANGLES,           &
                                 MAX_USER_RELAZMS, MAX_USER_VZANGLES, MAX_USER_STREAMS, MAX_USER_LEVELS, &
                                 MAX_USER_OBSGEOMS, MAXSTOKES_SQ, MAXFOURIER, MAX_GEOMETRIES,            &
                                 FMT_SECTION, FMT_HEADING, FMT_CHAR, FMT_REAL, FMT_INTEGER, ONE, UPIDX, DNIDX

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          RUN
      LOGICAL, INTENT (IN) ::          DO_FULLRAD_MODE
      LOGICAL, INTENT (IN) ::          DO_FOCORR_NADIR
      LOGICAL, INTENT (IN) ::          DO_FOCORR_OUTGOING
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

      LOGICAL, INTENT (IN) ::          DO_NO_AZIMUTH
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::          N_OUT_STREAMS

      DOUBLE PRECISION, INTENT (IN) :: OUT_ANGLES ( MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          VZA_OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )

      DOUBLE PRECISION, INTENT (IN) :: TAUGRID_INPUT ( 0:MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_MULTIBEAM  ( MAXBEAMS, 0:MAXFOURIER )

      DOUBLE PRECISION, INTENT (IN) :: STOKES        ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: MEAN_STOKES   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: FLUX_STOKES   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )
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

!  Removed for Version 2.8 7/6/16
  !    IF ( DO_CLASSICAL_SOLUTION ) THEN
  !      write(RUN,'(a)') &
  !           'Classical solution of beam particular integral'
  !    ELSE
  !      write(RUN,'(a)') &
  !           'Green function solution of beam particular integral'
  !    ENDIF

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
        IF ( DO_FOCORR_NADIR ) THEN
          write(RUN,'(a)') &
        '  --> Nakajima-Tanaka TMS single scatter correction (nadir)'
        ENDIF
        IF ( DO_FOCORR_OUTGOING ) THEN
          write(RUN,'(a)') &
        '  --> Nakajima-Tanaka TMS single scatter correction (outgoing)'
        ENDIF
        IF ( .NOT.DO_FOCORR_NADIR .AND. .NOT.DO_FOCORR_OUTGOING ) THEN
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

      SUBROUTINE VLIDORT_WRITESCEN ( &
        SUNIT, DO_DELTAM_SCALING, NSTREAMS, &
        NLAYERS, N_SZANGLES, SZANGLES,      &
        N_USER_RELAZMS, USER_RELAZMS,       &
        N_USER_VZANGLES, USER_VZANGLES,     &
        N_USER_LEVELS, USER_LEVELS,         &
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,  &
        DO_LAMBERTIAN_SURFACE, LAMBERTIAN_ALBEDO, &
        DO_SURFACE_EMISSION, SURFBB, &
        QUAD_STREAMS, QUAD_WEIGHTS,  &
        QUAD_ANGLES, DO_NO_AZIMUTH,  &
        NMOMENTS, GREEKMAT_INDEX,    &
        TAUGRID_INPUT, OMEGA_TOTAL,  &
        GREEKMAT_TOTAL, TAUGRID )

!  4/15/20. Version 2.8.2. 
!     -- Remove NGREEK_MOMENTS_INPUT from list, Remove MAXMOMENTS_INPUT from parameter list
!     -- Greek moments dimension, now MAXMOMENTS
!     -- Remove thermal coefficients

!  write to file of all geophysical VLIDORT input

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTREAMS, MAXLAYERS, MAX_SZANGLES, MAXSTOKES_SQ, ONE, DEG_TO_RAD, &
                                 MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAX_USER_LEVELS, MAX_USER_OBSGEOMS,        &
                                 VLIDORT_VERSION_NUMBER, FMT_SECTION, FMT_HEADING, FMT_CHAR, FMT_REAL 

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          SUNIT
      LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

      INTEGER, INTENT (IN) ::          N_SZANGLES
      DOUBLE PRECISION, INTENT (IN) :: SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER, INTENT (IN) ::          N_USER_VZANGLES
      DOUBLE PRECISION, INTENT (IN) :: USER_VZANGLES ( MAX_USER_VZANGLES )
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) :: USER_LEVELS ( MAX_USER_LEVELS )
      DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
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
      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) :: TAUGRID ( 0:MAXLAYERS )

!  4/15/20. Version 2.8.2. These removed
!      INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
!      DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

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
!  4/15/20. Version 2.8.2. NGREEK_MOMENTS_INPUT ==> NMOMENTS 

      IF ( DO_DELTAM_SCALING ) THEN
        WRITE(SUNIT, FMT_HEADING) 'Phase matrix elements (Greek constants), SCALED'
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(a, I3/)') 'Matrix elements (A-F) for layer ',N
          DO L = 0, NMOMENTS
            WRITE(SUNIT,'(I2,1x,1p6e12.4)') &
           L,(GREEKMAT_TOTAL(L,N,GREEKMAT_INDEX(S)),S = 1, 6)
          ENDDO
        ENDDO
      ELSE
        WRITE(SUNIT, FMT_HEADING) 'Phase matrix elements (Greek constants), UNSCALED'
        DO N = 1, NLAYERS
          WRITE(SUNIT,'(a, I3/)')  'Matrix elements (A-F) for layer ',N
          DO L = 0, 2*NSTREAMS
            WRITE(SUNIT,'(I2,1x,1p6e12.4)') L,(GREEKMAT_TOTAL_INPUT(L,N,GREEKMAT_INDEX(S)),S = 1, 6)
          ENDDO
        ENDDO
      ENDIF

!  Commented out thermal expansion coefficients

!      IF ( DO_THERMAL_EMISSION ) THEN
!        WRITE(SUNIT,FMT_HEADING)'thermal emission coefficients'
!        WRITE(SUNIT,'(a/)') 'Layer | thermal emission expansion coeffs-->'
!        DO N = 1, NLAYERS
!          WRITE(SUNIT,'(I3,4x,10(f10.5))')  N,(THERMAL_COEFFS(N,S),S=1,N_THERMAL_COEFFS)
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

      SUBROUTINE VLIDORT_WRITE_STD_INPUT ( DO_DEBUG_WRITE, &
        DO_SOLAR_SOURCES,   DO_THERMAL_EMISSION, DO_THERMAL_TRANSONLY, DO_SURFACE_EMISSION,       & ! Sources
        DO_FULLRAD_MODE,    DO_FOCORR, DO_FOCORR_EXTERNAL, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING,   & ! FOCORR
        DO_RAYLEIGH_ONLY,   DO_DELTAM_SCALING,   DO_DOUBLE_CONVTEST,   DO_SOLUTION_SAVING,   DO_BVP_TELESCOPING, & ! Performance
        DO_UPWELLING,       DO_DNWELLING,        DO_PLANE_PARALLEL,    DO_CHAPMAN_FUNCTION,       & ! RT Model
        DO_USER_VZANGLES,   DO_OBSERV_GEOMETRY,  DO_DOUBLET_GEOMETRY,  DO_ADDITIONAL_MVOUT,  DO_MVOUT_ONLY,      & ! RT model
        DO_REFRAC_GEOMETRY, DO_ALBTRN_MEDIA,     DO_PLANETARY_PROBLEM, DO_TOAFLUX, DO_BOAFLUX,    & ! Media/Planetary/Flux
        DO_LAMBERT_SURFACE, DO_SURFACE_LEAVING,  DO_SL_ISOTROPIC,      DO_EXTERNAL_WLEAVE,        & ! Surface
        DO_WATER_LEAVING,   DO_FLUORESCENCE,     DO_WLADJUSTED_OUTPUT, DO_TF_ITERATION,           & ! Surface/Media-props
        TAYLOR_ORDER, TF_MAXITER, NSTOKES, NSTREAMS, NLAYERS, NFINELAYERS, N_THERMAL_COEFFS,      & ! Main numbers
        DO_TOA_CONTRIBS,  DO_SPECIAL_1, DO_SPECIAL_2, DO_SPECIAL_3, NLAYERS_NOMS, NLAYERS_CUTOFF, & ! Specialist Options
        N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES, N_USER_OBSGEOMS, N_USER_LEVELS,              & ! Geometry and level numbers
        SZANGLES, USER_RELAZMS, USER_VZANGLES, USER_OBSGEOMS, USER_LEVELS,                        & ! Geometry and level control
        VLIDORT_ACCURACY, FLUX_FACTOR, EARTH_RADIUS, RFINDEX_PARAMETER, TF_CRITERION,             & ! Flux/Acc/Radius
        HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, FINEGRID, GEOMETRY_SPECHEIGHT,              & ! Grids
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                               & ! Optical
        FMATRIX_UP, FMATRIX_DN, LAMBERTIAN_ALBEDO, TOAFLUX, BOAFLUX,                              & ! Optical & albedo & Fluxes
        THERMAL_BB_INPUT, SURFACE_BB_INPUT, ATMOS_WAVELENGTH,                                     & ! Thermal & spectral
        DO_WRITE_INPUT,    DO_WRITE_SCENARIO,   DO_WRITE_FOURIER  , DO_WRITE_RESULTS,             & ! Debug control
        INPUT_WRITE_FILE,  SCENARIO_WRITE_FILE, FOURIER_WRITE_FILE, RESULTS_WRITE_FILE )            ! Debug control

!  3/28/14. Changes for Version 2.7. Remove LTE linearization references
!  3/1/17 . Changes for Version 2.8. Complete reorganization of argument lists
!mick fix 9/19/2017 - added TAYLOR_ORDER to argument list

!  Additional Control for SLEAVE (DO_WLADJUSTED_OUTPUT,DO_EXTERNAL_WLEAVE).
!     Introduced 3/18/19 for Version 2.8.1

!   Version 2.8.1, Controls for TOA/BOA isotropic illumination added, 3/23/19
!  --- Introduce flags and values for including TOA/BOA isotropic illumination

!  4/15/20. Version 2.8.2. Argument list rearranged/renamed, follows LIDORT Version 3.8.2.
!  4/15/20. Version 2.8.2. MAXMOMENTS instead of MAXMOMENTS_INPUT, NGREEK_MOMENTS_INPUT dropped

!  4/15/20. Add DO_DOUBLET_GEOMETRY flag to argument list
        
      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAX_SZANGLES, MAX_USER_RELAZMS, &
                                 MAX_USER_VZANGLES, MAX_USER_LEVELS, MAX_USER_OBSGEOMS, &
                                 MAX_GEOMETRIES, MAXSTOKES_SQ,                          &
                                 DWFC, DWFL, DWFI, DWFI1, DWFR, DWFR1, DWFR2, DWFR3, DWFR4, DWFR1_3

      IMPLICIT NONE

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

      LOGICAL, INTENT(IN) ::          DO_FULLRAD_MODE

      LOGICAL, INTENT(IN) ::          DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN) ::          DO_SURFACE_EMISSION
      LOGICAL, INTENT(IN) ::          DO_PLANE_PARALLEL

      LOGICAL, INTENT(IN) ::          DO_UPWELLING
      LOGICAL, INTENT(IN) ::          DO_DNWELLING

      LOGICAL, INTENT(IN) ::          DO_LAMBERT_SURFACE
      LOGICAL, INTENT(IN) ::          DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN) ::          DO_SL_ISOTROPIC
      LOGICAL, INTENT(IN) ::          DO_WATER_LEAVING
      LOGICAL, INTENT(IN) ::          DO_FLUORESCENCE
      LOGICAL, INTENT(IN) ::          DO_TF_ITERATION

      LOGICAL, INTENT(IN) ::          DO_TOA_CONTRIBS
      LOGICAL, INTENT(IN) ::          DO_SPECIAL_1
      LOGICAL, INTENT(IN) ::          DO_SPECIAL_2
      LOGICAL, INTENT(IN) ::          DO_SPECIAL_3

      LOGICAL, INTENT(IN) ::          DO_TOAFLUX ! New 3/23/19
      LOGICAL, INTENT(IN) ::          DO_BOAFLUX ! New 3/23/19

!  integers

      INTEGER, INTENT(IN) ::          TAYLOR_ORDER
      INTEGER, INTENT(IN) ::          TF_MAXITER

      INTEGER, INTENT(IN) ::          NSTOKES
      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          NLAYERS
      INTEGER, INTENT(IN) ::          NFINELAYERS
      INTEGER, INTENT(IN) ::          N_THERMAL_COEFFS
      INTEGER, INTENT(IN) ::          N_USER_LEVELS

      INTEGER, INTENT(IN) ::          NLAYERS_NOMS
      INTEGER, INTENT(IN) ::          NLAYERS_CUTOFF

!  Reals

      DOUBLE PRECISION, INTENT(IN) :: VLIDORT_ACCURACY
      DOUBLE PRECISION, INTENT(IN) :: FLUX_FACTOR

      DOUBLE PRECISION, INTENT(IN) :: HEIGHT_GRID      ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: PRESSURE_GRID    ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: TEMPERATURE_GRID ( 0:MAXLAYERS )
      INTEGER, INTENT(IN) ::          FINEGRID ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: RFINDEX_PARAMETER
      DOUBLE PRECISION, INTENT(IN) :: TF_CRITERION
      DOUBLE PRECISION, INTENT(IN) :: EARTH_RADIUS

!  4/15/20. Version 2.8.2. MAXMOMENTS instead of MAXMOMENTS_INPUT

      DOUBLE PRECISION, INTENT(IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT(IN) :: FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 )
      DOUBLE PRECISION, INTENT(IN) :: FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 )

      DOUBLE PRECISION, INTENT(IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT(IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) :: SURFACE_BB_INPUT
      DOUBLE PRECISION, INTENT(IN) :: ATMOS_WAVELENGTH

!  TOA/BOA Fluxes (new 3/23/19). Version 2.8.1

      DOUBLE PRECISION, INTENT (IN) :: TOAFLUX
      DOUBLE PRECISION, INTENT (IN) :: BOAFLUX

!  These two added for Version 2.8.1, 4/26/19, 4/28/19

      LOGICAL, INTENT(IN) ::          DO_ALBTRN_MEDIA(2)
      LOGICAL, INTENT(IN) ::          DO_PLANETARY_PROBLEM
      
!  New flags for Version 2.8.1. Introduced 3/18/19.
      
      LOGICAL, INTENT(IN) ::          DO_WLADJUSTED_OUTPUT
      LOGICAL, INTENT(IN) ::          DO_EXTERNAL_WLEAVE

!  Code superseded, Version 2.7. Remove this
      !DOUBLE PRECISION, INTENT(IN) ::   LTE_DELTAU_VERT_INPUT ( 2, MAXLAYERS )
      !DOUBLE PRECISION, INTENT(IN) ::   LTE_THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  debug

      LOGICAL, INTENT(IN) ::          DO_DEBUG_WRITE
      LOGICAL, INTENT(IN) ::          DO_WRITE_INPUT
      LOGICAL, INTENT(IN) ::          DO_WRITE_SCENARIO
      LOGICAL, INTENT(IN) ::          DO_WRITE_FOURIER
      LOGICAL, INTENT(IN) ::          DO_WRITE_RESULTS

      CHARACTER (LEN=60), INTENT(IN) :: INPUT_WRITE_FILE
      CHARACTER (LEN=60), INTENT(IN) :: SCENARIO_WRITE_FILE
      CHARACTER (LEN=60), INTENT(IN) :: FOURIER_WRITE_FILE
      CHARACTER (LEN=60), INTENT(IN) :: RESULTS_WRITE_FILE

!  --------------------------
!  Standard Inputs - Variable
!  --------------------------
!mick mod 9/19/2017 - DO_FOCORR_ALONE now defined internally

      LOGICAL, INTENT(IN) ::          DO_FOCORR
      LOGICAL, INTENT(IN) ::          DO_FOCORR_EXTERNAL
      !LOGICAL, INTENT(IN) ::          DO_FOCORR_ALONE
      LOGICAL, INTENT(IN) ::          DO_FOCORR_NADIR
      LOGICAL, INTENT(IN) ::          DO_FOCORR_OUTGOING
      
!  4/15/20. Version 2.8.2. These two variables no longer required.
!      LOGICAL, INTENT(IN) ::          DO_SSCORR_TRUNCATION
!      LOGICAL, INTENT(IN) ::          DO_SSCORR_USEFMAT

!  4/15/20. Add DO_DOUBLET_GEOMETRY flag to argument list
      
      LOGICAL, INTENT(IN) :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN) :: DO_REFRAC_GEOMETRY
      LOGICAL, INTENT(IN) :: DO_CHAPMAN_FUNCTION

      LOGICAL, INTENT(IN) :: DO_RAYLEIGH_ONLY

      LOGICAL, INTENT(IN) :: DO_DOUBLE_CONVTEST
      LOGICAL, INTENT(IN) :: DO_DELTAM_SCALING
      LOGICAL, INTENT(IN) :: DO_SOLUTION_SAVING
      LOGICAL, INTENT(IN) :: DO_BVP_TELESCOPING

      LOGICAL, INTENT(IN) :: DO_USER_VZANGLES
      LOGICAL, INTENT(IN) :: DO_OBSERV_GEOMETRY
      LOGICAL, INTENT(IN) :: DO_DOUBLET_GEOMETRY
      
      LOGICAL, INTENT(IN) :: DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT(IN) :: DO_MVOUT_ONLY

      LOGICAL, INTENT(IN) :: DO_THERMAL_TRANSONLY

!  Geometry

      INTEGER, INTENT(IN) ::          N_SZANGLES
      DOUBLE PRECISION, INTENT(IN) :: SZANGLES ( MAX_SZANGLES )

      INTEGER, INTENT(IN) ::          N_USER_RELAZMS
      DOUBLE PRECISION, INTENT(IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )

      INTEGER, INTENT(IN) ::          N_USER_VZANGLES
      DOUBLE PRECISION, INTENT(IN) :: USER_VZANGLES ( MAX_USER_VZANGLES )

      DOUBLE PRECISION, INTENT(IN) :: USER_LEVELS   ( MAX_USER_LEVELS )
      DOUBLE PRECISION, INTENT(IN) :: GEOMETRY_SPECHEIGHT

      INTEGER, INTENT(IN) ::          N_USER_OBSGEOMS
      DOUBLE PRECISION, INTENT(IN) :: USER_OBSGEOMS ( MAX_USER_OBSGEOMS, 3 )

!  Local variables

      INTEGER :: OUTUNIT, NGEOMS
      INTEGER :: GEO,K,LAY,MOM,NGREEK,SS,SZA,ULEV,URA,UVA,UOG
      INTEGER :: GMASK(8)
      DATA GMASK / 1, 2, 5, 6, 11, 12, 15, 16 /

!  Define local variables

      NGEOMS = N_SZANGLES
      IF ( .NOT.DO_OBSERV_GEOMETRY ) NGEOMS = N_SZANGLES * N_USER_VZANGLES * N_USER_RELAZMS

!  Open output file

      OUTUNIT = 101
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_STD_INPUT.dbg',status = 'replace')

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
      WRITE(OUTUNIT,DWFL)  'DO_THERMAL_EMISSION    = ',DO_THERMAL_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_EMISSION    = ',DO_SURFACE_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_PLANE_PARALLEL      = ',DO_PLANE_PARALLEL
      WRITE(OUTUNIT,DWFL)  'DO_UPWELLING           = ',DO_UPWELLING
      WRITE(OUTUNIT,DWFL)  'DO_DNWELLING           = ',DO_DNWELLING

      !WRITE(OUTUNIT,DWFL)  'DO_QUAD_OUTPUT         = ',DO_QUAD_OUTPUT  ! removed Version 2.8
      WRITE(OUTUNIT,DWFL)  'DO_TOA_CONTRIBS        = ',DO_TOA_CONTRIBS
      WRITE(OUTUNIT,DWFL)  'DO_LAMBERTIAN_SURFACE  = ',DO_LAMBERT_SURFACE
      WRITE(OUTUNIT,DWFL)  'DO_SPECIALIST_OPTION_1 = ',DO_SPECIAL_1
      WRITE(OUTUNIT,DWFL)  'DO_SPECIALIST_OPTION_2 = ',DO_SPECIAL_2
      WRITE(OUTUNIT,DWFL)  'DO_SPECIALIST_OPTION_3 = ',DO_SPECIAL_3

      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_LEAVING     = ',DO_SURFACE_LEAVING
      WRITE(OUTUNIT,DWFL)  'DO_SL_ISOTROPIC        = ',DO_SL_ISOTROPIC
      WRITE(OUTUNIT,DWFL)  'DO_WATER_LEAVING       = ',DO_WATER_LEAVING
      WRITE(OUTUNIT,DWFL)  'DO_FLUORESCENCE        = ',DO_FLUORESCENCE
      WRITE(OUTUNIT,DWFL)  'DO_TF_ITERATION        = ',DO_TF_ITERATION
      WRITE(OUTUNIT,DWFL)  'DO_WLADJUSTED_OUTPUT   = ',DO_WLADJUSTED_OUTPUT   !    Introduced 3/18/19 for Version 2.8.1
      WRITE(OUTUNIT,DWFL)  'DO_TOA_ILLUMINATION    = ',DO_TOAFLUX             !    Introduced 3/23/19 for Version 2.8.1
      WRITE(OUTUNIT,DWFL)  'DO_BOA_ILLUMINATION    = ',DO_BOAFLUX             !    Introduced 3/23/19 for Version 2.8.1

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'TAYLOR_ORDER     = ',TAYLOR_ORDER
      WRITE(OUTUNIT,DWFI)  'NSTOKES          = ',NSTOKES
      WRITE(OUTUNIT,DWFI)  'NSTREAMS         = ',NSTREAMS
      WRITE(OUTUNIT,DWFI)  'NLAYERS          = ',NLAYERS
      WRITE(OUTUNIT,DWFI)  'NFINELAYERS      = ',NFINELAYERS
      WRITE(OUTUNIT,DWFI)  'N_THERMAL_COEFFS = ',N_THERMAL_COEFFS

      WRITE(OUTUNIT,DWFR)  'VLIDORT_ACCURACY = ',VLIDORT_ACCURACY

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NLAYERS_NOMS     = ',NLAYERS_NOMS
      WRITE(OUTUNIT,DWFI)  'NLAYERS_CUTOFF   = ',NLAYERS_CUTOFF

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'TF_MAXITER       = ',TF_MAXITER
      WRITE(OUTUNIT,DWFR)  'TF_CRITERION     = ',TF_CRITERION

!  4/26/19. Record the Media-problem and Planetary-problem inputs

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_ALBTRN_MEDIA(1)   = ',DO_ALBTRN_MEDIA(1)
      WRITE(OUTUNIT,DWFL)  'DO_ALBTRN_MEDIA(2)   = ',DO_ALBTRN_MEDIA(2)
      WRITE(OUTUNIT,DWFL)  'DO_PLANETARY_PROBLEM = ',DO_PLANETARY_PROBLEM

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'FLUX_FACTOR      = ',FLUX_FACTOR
      
      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'TOA_ILLUMINATION = ',TOAFLUX          !    Introduced 3/23/19 for Version 2.8.1
      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'BOA_ILLUMINATION = ',BOAFLUX          !    Introduced 3/23/19 for Version 2.8.1

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_LEVELS    = ',N_USER_LEVELS

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' HEIGHT_GRID(LAY)      = ',HEIGHT_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' PRESSURE_GRID(LAY)    = ',PRESSURE_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' TEMPERATURE_GRID(LAY) = ',TEMPERATURE_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFI1)  'LAY = ',LAY,' FINEGRID (LAY)        = ',FINEGRID(LAY)
      END DO

 !  4/15/20. Version 2.8.2. EARTH_RADIUS moved here

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'RFINDEX_PARAMETER = ',RFINDEX_PARAMETER
      WRITE(OUTUNIT,DWFR)  'EARTH_RADIUS        = ',EARTH_RADIUS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' DELTAU_VERT_INPUT(LAY) = ',DELTAU_VERT_INPUT(LAY)
      END DO

 !  4/15/20. Version 2.8.2. OMEGA_TOTAL_INPUT moved here

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' OMEGA_TOTAL_INPUT(LAY) = ',OMEGA_TOTAL_INPUT(LAY)
      END DO

 !  4/15/20. Version 2.8.2. Only have 2*NSTREAMS as the number of usable moments

      WRITE(OUTUNIT,*)
      DO K=1,NGREEK
        SS = GMASK(K)
        DO LAY=1,NLAYERS
        DO MOM=0,2*NSTREAMS
            WRITE(OUTUNIT,DWFR3)  'SS = ',SS,' LAY = ',LAY,' MOM = ',MOM,&
              ' GREEKMAT_TOTAL_INPUT(MOM,LAY,SS) = ',GREEKMAT_TOTAL_INPUT(MOM,LAY,SS)
          END DO
        END DO
      END DO

!  New section for 2.8 code. F-Matrices

      WRITE(OUTUNIT,*)
      DO K=1,6
        DO GEO=1,NGEOMS
          DO LAY=1,NLAYERS
            WRITE(OUTUNIT,DWFR3)  'K = ',K,' GEO = ',GEO,' LAY = ',LAY,&
              ' FMATRIX_UP(LAY,GEO,K) = ',FMATRIX_UP(LAY,GEO,K)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO K=1,6
        DO GEO=1,NGEOMS
          DO LAY=1,NLAYERS
            WRITE(OUTUNIT,DWFR3)  'K = ',K,' GEO = ',GEO,' LAY = ',LAY,&
              ' FMATRIX_DN(LAY,GEO,K) = ',FMATRIX_DN(LAY,GEO,K)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'LAMBERTIAN_ALBEDO = ',LAMBERTIAN_ALBEDO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' THERMAL_BB_INPUT(LAY) = ',THERMAL_BB_INPUT(LAY)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'SURFACE_BB_INPUT  = ',SURFACE_BB_INPUT

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'ATMOS_WAVELENGTH  = ',ATMOS_WAVELENGTH

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

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFC)  'INPUT_WRITE_FILENAME    = ',INPUT_WRITE_FILE
      WRITE(OUTUNIT,DWFC)  'SCENARIO_WRITE_FILENAME = ',SCENARIO_WRITE_FILE
      WRITE(OUTUNIT,DWFC)  'FOURIER_WRITE_FILENAME  = ',FOURIER_WRITE_FILE
      WRITE(OUTUNIT,DWFC)  'RESULTS_WRITE_FILENAME  = ',RESULTS_WRITE_FILE

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '--------------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Variable'
      WRITE(OUTUNIT,'(A)') '--------------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR               = ',DO_FOCORR
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR_EXTERNAL      = ',DO_FOCORR_EXTERNAL
      !WRITE(OUTUNIT,DWFL)  'DO_FOCORR_ALONE         = ',DO_FOCORR_ALONE
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR_NADIR         = ',DO_FOCORR_NADIR
      WRITE(OUTUNIT,DWFL)  'DO_FOCORR_OUTGOING      = ',DO_FOCORR_OUTGOING
!      WRITE(OUTUNIT,DWFL)  'DO_SSCORR_TRUNCATION    = ',DO_SSCORR_TRUNCATION
!      WRITE(OUTUNIT,DWFL)  'DO_SSCORR_USEFMAT       = ',DO_SSCORR_USEFMAT

      WRITE(OUTUNIT,DWFL)  'DO_EXTERNAL_WLEAVE      = ',DO_EXTERNAL_WLEAVE   !    Introduced 3/18/19 for Version 2.8.1
      WRITE(OUTUNIT,DWFL)  'DO_DOUBLE_CONVTEST      = ',DO_DOUBLE_CONVTEST
      WRITE(OUTUNIT,DWFL)  'DO_SOLAR_SOURCES        = ',DO_SOLAR_SOURCES
      WRITE(OUTUNIT,DWFL)  'DO_REFRACTIVE_GEOMETRY  = ',DO_REFRAC_GEOMETRY
      WRITE(OUTUNIT,DWFL)  'DO_CHAPMAN_FUNCTION     = ',DO_CHAPMAN_FUNCTION
      WRITE(OUTUNIT,DWFL)  'DO_RAYLEIGH_ONLY        = ',DO_RAYLEIGH_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_DELTAM_SCALING       = ',DO_DELTAM_SCALING
      WRITE(OUTUNIT,DWFL)  'DO_SOLUTION_SAVING      = ',DO_SOLUTION_SAVING
      WRITE(OUTUNIT,DWFL)  'DO_BVP_TELESCOPING      = ',DO_BVP_TELESCOPING
      WRITE(OUTUNIT,DWFL)  'DO_USER_VZANGLES        = ',DO_USER_VZANGLES
      WRITE(OUTUNIT,DWFL)  'DO_ADDITIONAL_MVOUT     = ',DO_ADDITIONAL_MVOUT
      WRITE(OUTUNIT,DWFL)  'DO_MVOUT_ONLY           = ',DO_MVOUT_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_THERMAL_TRANSONLY    = ',DO_THERMAL_TRANSONLY
      WRITE(OUTUNIT,DWFL)  'DO_OBSERVATION_GEOMETRY = ',DO_OBSERV_GEOMETRY
      
!  4/15/20. Add DO_DOUBLET_GEOMETRY flag to output list

      WRITE(OUTUNIT,DWFL)  'DO_DOUBLET_GEOMETRY     = ',DO_DOUBLET_GEOMETRY
      
      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NMOMENTS = (2*NSTREAMS) ',2*NSTREAMS

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_SZANGLES     = ',N_SZANGLES
      DO SZA=1,N_SZANGLES
        WRITE(OUTUNIT,DWFR1)  'SZA = ',SZA,' SZANGLES(SZA) = ',SZANGLES(SZA)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_RELAZMS = ',N_USER_RELAZMS
      DO URA=1,N_USER_RELAZMS
        WRITE(OUTUNIT,DWFR1)  'URA = ',URA,' USER_RELAZMS(URA) = ',USER_RELAZMS(URA)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_VZANGLES = ',N_USER_VZANGLES
      DO UVA=1,N_USER_VZANGLES
        WRITE(OUTUNIT,DWFR1)  'UVA = ',UVA,' USER_VZANGLES(UVA) = ',USER_VZANGLES(UVA)
      END DO

      WRITE(OUTUNIT,*)
      DO ULEV=1,N_USER_LEVELS
        WRITE(OUTUNIT,DWFR1)  'ULEV = ',ULEV,' USER_LEVELS(ULEV)  = ',USER_LEVELS(ULEV)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'GEOMETRY_SPECHEIGHT = ',GEOMETRY_SPECHEIGHT
      WRITE(OUTUNIT,DWFI)  'N_USER_OBSGEOMS     = ',N_USER_OBSGEOMS

      IF (N_USER_OBSGEOMS > 0) WRITE(OUTUNIT,*)
      DO UOG=1,N_USER_OBSGEOMS
        WRITE(OUTUNIT,DWFR1_3)  'UOG = ',UOG,&
          ' USER_OBSGEOMS(UOG,1:3) = ',USER_OBSGEOMS(UOG,1),&
                                   ',',USER_OBSGEOMS(UOG,2),&
                                   ',',USER_OBSGEOMS(UOG,3)
      END DO

!  5/22/20. Verison 2.8.2 Upgrade. This code moved up the routie a bit.

!      WRITE(OUTUNIT,*)
!      WRITE(OUTUNIT,DWFR)  'EARTH_RADIUS        = ',EARTH_RADIUS
!      WRITE(OUTUNIT,*)
!      DO LAY=1,NLAYERS
!        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' OMEGA_TOTAL_INPUT(LAY) = ',OMEGA_TOTAL_INPUT(LAY)
!      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_STD_INPUT

!

      SUBROUTINE VLIDORT_WRITE_SUP_BRDF_INPUT ( &
        DO_USER_STREAMS, DO_SURFACE_EMISSION, &
        NSTOKES,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,NSTREAMS, &
        EXACTDB_BRDFUNC,BRDF_F_0,BRDF_F,USER_BRDF_F_0,USER_BRDF_F,  &
        EMISSIVITY,USER_EMISSIVITY)

!mick mod 9/19/2017 - added DO_USER_STREAMS & DO_SURFACE_EMISSION for better output control
!                   - added NMOMENTS for VLIDORT internal consistency

!  4/15/20. Version 2.8.2. Remove NMOMENTS and MAXMOMENTS, Only printing First Fourier component.

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS, & 
                                 MAXSTOKES_SQ, DWFR2, DWFR4

      IMPLICIT NONE

      LOGICAL, INTENT(IN) ::   DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::   DO_SURFACE_EMISSION

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   N_SZANGLES
      INTEGER, INTENT(IN) ::   N_USER_VZANGLES
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   NSTREAMS
!      INTEGER, INTENT(IN) ::   NMOMENTS

      DOUBLE PRECISION, INTENT(IN) ::   EXACTDB_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  4/15/20. Version 2.8.2.  Only have the M = 0 Component to write now.
!      DOUBLE PRECISION, INTENT(IN) ::   BRDF_F_0      ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
!      DOUBLE PRECISION, INTENT(IN) ::   BRDF_F        ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
!      DOUBLE PRECISION, INTENT(IN) ::   USER_BRDF_F_0 ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
!      DOUBLE PRECISION, INTENT(IN) ::   USER_BRDF_F   ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

      DOUBLE PRECISION, INTENT(IN) ::   BRDF_F_0      ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   BRDF_F        ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT(IN) ::   USER_BRDF_F_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   USER_BRDF_F   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

      DOUBLE PRECISION, INTENT(IN) ::   EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT(IN) ::   USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SS,MOM,NREFL,K,STRM,S,USTRM,IB,URA,&
                 STRMI,STRMJ, RMASK3(9), RMASK4(16), RMASK(16)

      data RMASK3 / 1, 2, 3, 5, 6, 7, 9, 10, 11 /
      data RMASK4 / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 /

!  Open output file

      OUTUNIT = 102
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_SUP_BRDF_INPUT.dbg',status = 'replace')

!  Define local variables
!   Changed 16 October 2014, according to correct entries in 4x4 BRDF matrices

      NREFL = NSTOKES * NSTOKES ; RMASK(1) = 1
      IF ( NSTOKES.eq.3 ) RMASK(1:NREFL) = RMASK3(1:NREFL)
      IF ( NSTOKES.eq.4 ) RMASK(1:NREFL) = RMASK4(1:NREFL)

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '----------------------'
      WRITE(OUTUNIT,'(A)') 'BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '----------------------'

      WRITE(OUTUNIT,*)
      DO K=1,NREFL
        SS = RMASK(K)
        DO IB=1,N_SZANGLES
          DO URA=1,N_USER_RELAZMS
            DO USTRM=1,N_USER_VZANGLES
              WRITE(OUTUNIT,DWFR4) &
                'SS = ',SS,' IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,&
                ' EXACTDB_BRDFUNC(SS,USTRM,URA,IB) = ',&
                  EXACTDB_BRDFUNC(SS,USTRM,URA,IB)
            END DO
          END DO
        END DO
      END DO

!  4/15/20. Version 2.8.2. Remove MOM loops
!  4/15/20. Version 2.8.2. BRDFs only defined locally, without the Fourier index dimension
!mick mod 3/2/2020 - added "MOM = 0" notes

      MOM = 0

      WRITE(OUTUNIT,*)
      IF ( .NOT.DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING TWO BRDF INPUT SETS:'
      ELSE
        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING FOUR BRDF INPUT SETS:'
      ENDIF

      WRITE(OUTUNIT,*)
      DO K=1,NREFL
        SS = RMASK(K)
        DO IB=1,N_SZANGLES ; DO STRM=1,NSTREAMS
           WRITE(OUTUNIT,DWFR4) &
                'SS = ',SS,' IB = ',IB,' STRM = ',STRM,' MOM = ',MOM,&
                ' BRDF_F_0(0,SS,STRM,IB) = ',BRDF_F_0(SS,STRM,IB)
        END DO ; END DO
      END DO

      WRITE(OUTUNIT,*)
      DO K=1,NREFL
        SS = RMASK(K)
        DO STRMJ=1,NSTREAMS ; DO STRMI=1,NSTREAMS
          WRITE(OUTUNIT,DWFR4) &
                'SS = ',SS,' STRMJ = ',STRMJ,' STRMI = ',STRMI,' MOM = ',MOM,&
                ' BRDF_F(0,SS,STRMI,STRMJ) = ',BRDF_F(SS,STRMI,STRMJ)
        END DO ; END DO
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO K=1,NREFL
          SS = RMASK(K)
          DO IB=1,N_SZANGLES ; DO USTRM=1,N_USER_VZANGLES
            WRITE(OUTUNIT,DWFR4) &
                  'SS = ',SS,' IB = ',IB,' USTRM = ',USTRM,' MOM = ',MOM,&
                  ' USER_BRDF_F_0(0,SS,USTRM,IB) = ',USER_BRDF_F_0(SS,USTRM,IB)
          END DO ; END DO
        END DO
        WRITE(OUTUNIT,*)
        DO K=1,NREFL
          SS = RMASK(K)
          DO STRM=1,NSTREAMS ; DO USTRM=1,N_USER_VZANGLES
            WRITE(OUTUNIT,DWFR4) &
              'SS = ',SS,' STRM = ',STRM,' USTRM = ',USTRM,' MOM = ',MOM,&
              ' USER_BRDF_F(0,SS,USTRM,STRM) = ',USER_BRDF_F(SS,USTRM,STRM)
          END DO ; END DO
        END DO
      END IF

      IF ( DO_SURFACE_EMISSION ) THEN
        WRITE(OUTUNIT,*)
        DO S=1,NSTOKES
          DO STRM=1,NSTREAMS
            WRITE(OUTUNIT,DWFR2)  'S = ',S,' STRM = ',STRM,&
              ' EMISSIVITY (S,STRM) = ',EMISSIVITY (S,STRM)
          END DO
        END DO

        WRITE(OUTUNIT,*)
        DO S=1,NSTOKES
          DO USTRM=1,N_USER_VZANGLES
            WRITE(OUTUNIT,DWFR2)  'S = ',S,' USTRM = ',USTRM,&
              ' USER_EMISSIVITY (S,USTRM) = ',USER_EMISSIVITY(S,USTRM)
          END DO
        END DO
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_SUP_BRDF_INPUT

!

      SUBROUTINE VLIDORT_WRITE_SUP_SS_INPUT ( &
        NSTOKES,N_USER_LEVELS,STOKES_SS,STOKES_DB)

      USE VLIDORT_PARS_m, Only : MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS, DWFR3, DWFR4

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   N_USER_LEVELS

      DOUBLE PRECISION, INTENT(IN) ::  STOKES_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) ::  STOKES_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: DIR,GEO,S,ULEV

!  Open output file

      OUTUNIT = 103
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_SUP_SS_INPUT.dbg',status = 'replace')

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

      WRITE(OUTUNIT,*)
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
        DO_USER_STREAMS,NSTOKES,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,NSTREAMS,&
        SLTERM_ISOTROPIC,SLTERM_USERANGLES,SLTERM_F_0,USER_SLTERM_F_0)

!mick mod 9/19/2017 - added DO_USER_STREAMS for better output control
!                   - added NMOMENTS for VLIDORT internal consistency

!  4/15/20. Version 2.8.2. Remove NMOMS and MAXMOMENTS, Only printing First Fourier componnet.

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAX_USER_STREAMS, &
                                 MAX_USER_RELAZMS, MAXBEAMS, DWFR2, DWFR4

      IMPLICIT NONE

      LOGICAL, INTENT(IN) ::   DO_USER_STREAMS

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   N_SZANGLES
      INTEGER, INTENT(IN) ::   N_USER_VZANGLES
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   NSTREAMS

      DOUBLE PRECISION, INTENT(IN) :: SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) :: SLTERM_USERANGLES ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

!      DOUBLE PRECISION, INTENT(IN) :: SLTERM_F_0 ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
!      DOUBLE PRECISION, INTENT(IN) :: USER_SLTERM_F_0  ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

      DOUBLE PRECISION, INTENT(IN) :: SLTERM_F_0      ( MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) :: USER_SLTERM_F_0 ( MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: S,USTRM,IB,URA,MOM,STRM

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_SUP_SLEAVE_INPUT.dbg',status = 'replace')

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

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO IB=1,N_SZANGLES
          DO URA=1,N_USER_RELAZMS
            DO USTRM=1,N_USER_VZANGLES
              DO S=1,NSTOKES
                WRITE(OUTUNIT,DWFR4) &
                  'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,' S = ',S,&
                  ' SLTERM_USERANGLES(S,USTRM,URA,IB) = ',SLTERM_USERANGLES(S,USTRM,URA,IB)
              END DO
            END DO
          END DO
        END DO
      END IF

!  4/15/20. Version 2.8.2. Remove MOM loop. Only require first moment
!mick mod 3/2/2020 - added "MOM = 0" notes

      MOM = 0

      WRITE(OUTUNIT,*)
      IF ( .NOT.DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING SLEAVE INPUT SET:'
      ELSE
        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING TWO SLEAVE INPUT SETS:'
      ENDIF

      WRITE(OUTUNIT,*)
      DO IB=1,N_SZANGLES
        DO STRM=1,NSTREAMS ; DO S=1,NSTOKES
          WRITE(OUTUNIT,DWFR4) &
            'IB = ',IB,' STRM = ',STRM,' S = ',S,' MOM = ',MOM,&
            ' SLTERM_F_0(0,S,STRM,IB) = ',SLTERM_F_0(S,STRM,IB)
        END DO ; END DO
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO IB=1,N_SZANGLES ; DO USTRM=1,N_USER_VZANGLES
          DO S=1,NSTOKES
             WRITE(OUTUNIT,DWFR4) &
               'IB = ',IB,' USTRM = ',USTRM,' S = ',S,' MOM = ',MOM,&
                  ' USER_SLTERM_F_0(0,S,USTRM,IB) = ',USER_SLTERM_F_0(S,USTRM,IB)
          END DO
        END DO ; END DO
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_SUP_SLEAVE_INPUT

      END MODULE vlidort_writemodules_m

