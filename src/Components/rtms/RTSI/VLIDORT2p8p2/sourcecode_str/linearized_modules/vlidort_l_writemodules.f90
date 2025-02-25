
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
! #              VLIDORT_L_WRITERESULTS                         #
! #              VLIDORT_L_WRITESCEN                            #
! #              VLIDORT_WRITE_LIN_INPUT                        #
! #              VLIDORT_WRITE_LIN_SUP_BRDF_INPUT               #
! #              VLIDORT_WRITE_LIN_SUP_SS_INPUT                 #
! #              VLIDORT_WRITE_LIN_SUP_SLEAVE_INPUT             #
! #                                                             #
! ###############################################################

!  Version 2.8 Changes, 7/20/16
!     1. the flag DO_SURFBB_LINEARIZATION is gone.

!  4/15/20. Version 2.8.2. NMOMENTS instead of NGREEK_MOMENTS_INPUT, NMOMS, Remove it.
!  4/15/20. Version 2.8.2. MAXMOMENTS instead of MAXMOMENTS_INPUT.
!  4/15/20. Version 2.8.2. Remove VLIDORT_L_WRITEINPUT routine (superceded)

      MODULE vlidort_l_writemodules_m

      PRIVATE
      PUBLIC :: VLIDORT_L_WRITERESULTS,           &
                VLIDORT_L_WRITESCEN,              &
                VLIDORT_WRITE_LIN_INPUT,          &
                VLIDORT_WRITE_LIN_SUP_BRDF_INPUT, &
                VLIDORT_WRITE_LIN_SUP_SS_INPUT,   &
                VLIDORT_WRITE_LIN_SUP_SLEAVE_INPUT

      CONTAINS

      SUBROUTINE VLIDORT_L_WRITERESULTS ( &
        RUN, DO_FULLRAD_MODE, DO_SSCORR_NADIR,     &
        DO_SSCORR_OUTGOING, DO_DOUBLE_CONVTEST,    &
        DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, &
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,        &
        NSTOKES, NLAYERS, VLIDORT_ACCURACY,        &
        SZANGLES, N_USER_RELAZMS, USER_RELAZMS,    &
        N_USER_LEVELS, USER_LEVELS,                &
        DO_LAMBERTIAN_SURFACE,                     &
        DO_NO_AZIMUTH, NBEAMS,                     &
        N_DIRECTIONS, WHICH_DIRECTIONS,            &
        N_OUT_STREAMS,                             &
        OUT_ANGLES, VZA_OFFSETS,                   &
        DO_MULTIBEAM, FOURIER_SAVED,               &
        DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION,     &
        N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        DO_SURFACE_LINEARIZATION, N_SURFACE_WFS, &
        PROFILEWF_NAMES, COLUMNWF_NAMES, SURFACEWF, PROFILEWF, COLUMNWF,  &
        MINT_SURFACEWF, MINT_ATMOSWF, FLUX_SURFACEWF, FLUX_ATMOSWF )

      USE VLIDORT_PARS_m, Only : MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS, MAXLAYERS, MAX_SZANGLES,           &
                                 MAX_USER_RELAZMS, MAX_USER_VZANGLES, MAX_USER_STREAMS, MAX_USER_LEVELS, &
                                 MAX_USER_OBSGEOMS, MAXSTOKES_SQ, MAXFOURIER, MAX_GEOMETRIES,            &
                                 MAX_ATMOSWFS, MAX_SURFACEWFS, ONE, UPIDX, DNIDX,                        &
                                 FMT_SECTION, FMT_HEADING, FMT_CHAR, FMT_REAL, FMT_INTEGER


      IMPLICIT NONE

      INTEGER, INTENT (IN) ::            RUN
      LOGICAL, INTENT (IN) ::            DO_FULLRAD_MODE
      LOGICAL, INTENT (IN) ::            DO_SSCORR_NADIR
      LOGICAL, INTENT (IN) ::            DO_SSCORR_OUTGOING
      LOGICAL, INTENT (IN) ::            DO_DOUBLE_CONVTEST
      LOGICAL, INTENT (IN) ::            DO_PLANE_PARALLEL
      LOGICAL, INTENT (IN) ::            DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT (IN) ::            DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (IN) ::            DO_MVOUT_ONLY
      INTEGER, INTENT (IN) ::            NSTOKES
      INTEGER, INTENT (IN) ::            NLAYERS
      DOUBLE PRECISION, INTENT (IN) ::   VLIDORT_ACCURACY
      DOUBLE PRECISION, INTENT (IN) ::   SZANGLES ( MAX_SZANGLES )
      INTEGER, INTENT (IN) ::            N_USER_RELAZMS
      DOUBLE PRECISION, INTENT (IN) ::   USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER, INTENT (IN) ::            N_USER_LEVELS
      DOUBLE PRECISION, INTENT (IN) ::   USER_LEVELS ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::            DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::            DO_NO_AZIMUTH
      INTEGER, INTENT (IN) ::            NBEAMS
      INTEGER, INTENT (IN) ::            N_DIRECTIONS
      INTEGER, INTENT (IN) ::            WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::            N_OUT_STREAMS
      DOUBLE PRECISION, INTENT (IN) ::   OUT_ANGLES ( MAX_USER_STREAMS )
      INTEGER, INTENT (IN) ::            VZA_OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )

      LOGICAL, INTENT (IN) ::            DO_MULTIBEAM ( MAXBEAMS, 0:MAXFOURIER )
      INTEGER, INTENT (IN) ::            FOURIER_SAVED ( MAX_SZANGLES )
      LOGICAL, INTENT (IN) ::            DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::            DO_COLUMN_LINEARIZATION
      INTEGER, INTENT (IN) ::            N_TOTALCOLUMN_WFS
      LOGICAL, INTENT (IN) ::            LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::            LAYER_VARY_NUMBER ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::            DO_SURFACE_LINEARIZATION
      INTEGER, INTENT (IN) ::            N_SURFACE_WFS
      CHARACTER (LEN=31), INTENT (IN) :: PROFILEWF_NAMES ( MAX_ATMOSWFS )
      CHARACTER (LEN=31), INTENT (IN) :: COLUMNWF_NAMES  ( MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) ::   SURFACEWF &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::   PROFILEWF &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::   COLUMNWF &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT (IN) ::   MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::   MINT_ATMOSWF &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT (IN) ::   FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::   FLUX_ATMOSWF &
          ( MAX_ATMOSWFS, 0:MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER ::            I, UA, LOCAL_NUSERAZMS, IDIR, UT, WDIR, Z, O1
      INTEGER ::            Q, K, IB, V, FMAX, F
      INTEGER ::            VINDEX ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS )
      CHARACTER (LEN=11) :: STOKESNAM(4)

      STOKESNAM = (/ 'I-component','Q-component', &
                     'U-component','V-component'/)

!  Beam Attenuation summary
!  ========================

      write(RUN,'(/a/a/)')'Results Output summary', &
                           '----------------------'

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
           WRITE(RUN,'(1X,F7.2,4X,I3,6X,50(L1))') &
                  SZANGLES(IB),FOURIER_SAVED(IB), &
                  (DO_MULTIBEAM(IB,F),F=0,FMAX)
        END DO
      ELSE
        write(RUN,'(/a)') &
         'Azimuth independent output only (Fourier = 0)'
      ENDIF

!  local number of azimuths

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_NUSERAZMS = 1
      ELSE
        LOCAL_NUSERAZMS = N_USER_RELAZMS
      ENDIF

!  Indexing

      DO IB = 1, NBEAMS
       DO I = 1, N_OUT_STREAMS
        DO UA = 1, LOCAL_NUSERAZMS
         VINDEX(IB,I,UA) = VZA_OFFSETS(IB,I) + UA
        ENDDO
       ENDDO
      ENDDO

!  Atmospheric Profile Weighting function output
!  #############################################

      IF ( DO_PROFILE_LINEARIZATION ) THEN

!  control point for avoiding weighting function output

        IF ( DO_MVOUT_ONLY ) GO TO 401

!  overall header

        write(RUN,'(/a/a/)') &
             'Atmospheric Profile Weighting function output', &
             '---------------------------------------------'

!  start beam loop

        DO IB = 1, NBEAMS
         write(RUN,FMT_REAL) &
       '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
         write(RUN,'(a)')' '

!  start azimuth loop

         DO UA = 1, LOCAL_NUSERAZMS

!  azimuth angle header

          IF ( DO_NO_AZIMUTH ) THEN
            write(RUN,FMT_CHAR) &
             '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
          ELSE
            write(RUN,FMT_REAL) &
       '** RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
          ENDIF

          WRITE(RUN,'(a)')' '
          WRITE(RUN,FMT_INTEGER) &
            'Total number of output levels = ',N_USER_LEVELS
          WRITE(RUN,FMT_INTEGER) &
            'Total number of output angles = ',N_OUT_STREAMS

!  detailed output

          DO IDIR = 1, N_DIRECTIONS
           WDIR = WHICH_DIRECTIONS(IDIR)

!  linearization control

           DO K = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(K) ) THEN
             DO Q = 1, LAYER_VARY_NUMBER(K)

!  direction headers

              IF (WDIR .EQ. UPIDX ) THEN
               WRITE(RUN,'(/A/A,I2/A,A31)') &
           '  --> Upwelling Atmospheric Weighting functions', &
           '        * for variations in layer = ',K, &
           '        * with respect to : ',PROFILEWF_NAMES(Q)
              ELSE IF (WDIR .EQ. DNIDX ) THEN
               WRITE(RUN,'(/A/A,I2/A,A31)') &
           '  --> Downwelling Atmospheric Weighting functions', &
           '        * for variations in layer = ',K, &
           '        * with respect to : ',PROFILEWF_NAMES(Q)
              ENDIF

!  output loop

              DO O1 = 1, NSTOKES
               write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
               write(RUN,'(a)') 'output angles | output levels ---> '
               WRITE(RUN,'(14x,50(2x,1pe11.4,2x))') &
                 (USER_LEVELS(UT), UT = 1, N_USER_LEVELS)
               DO I = 1, N_OUT_STREAMS
                V = VZA_OFFSETS(IB,I) + UA
                WRITE(RUN,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I), &
               (PROFILEWF(Q,K,UT,V,O1,WDIR), UT = 1, N_USER_LEVELS)
               ENDDO
              ENDDO

!  close control and loops

             ENDDO
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO

!  Mean-value output
!  -----------------

401     CONTINUE
        IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

          write(RUN,'(/a/a)') &
                'Mean value atmospheric weighting function output', &
                '------------------------------------------------'

!  start beam loop

          DO IB = 1, NBEAMS

           write(RUN,FMT_REAL) &
       '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
           write(RUN,'(a)')' '

!  detailed output

           DO IDIR = 1, N_DIRECTIONS
            WDIR = WHICH_DIRECTIONS(IDIR)

!  linearization control

            DO K = 1, NLAYERS
             IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)

!  direction headers

               IF (WDIR .EQ. UPIDX ) THEN
                WRITE(RUN,'(/A/A,I2/A,I2/)') &
       '  --> Upwelling Mean intensity and Flux Weighting functions', &
           '        * W.R.T parameters varying in layer = ',K, &
           '        * Parameter number = ',Q
               ELSE IF (WDIR .EQ. DNIDX ) THEN
                WRITE(RUN,'(/A/A,I2/A,I2/)') &
       '  --> Downwelling Mean intensity and Flux Weighting functions', &
           '        * W.R.T parameters varying in layer = ',K, &
           '        * Parameter number = ',Q
               ENDIF

!  output

               DO O1 = 1, NSTOKES
                write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
                WRITE(RUN,'(a/)') &
            'output level     mean int. WF      flux WF'
                DO UT = 1, N_USER_LEVELS
                WRITE(RUN,'(2x,F9.5,3x,2E17.7)') &
                          USER_LEVELS(UT), &
                          MINT_ATMOSWF(Q,K,UT,IB,O1,WDIR), &
                          FLUX_ATMOSWF(Q,K,UT,IB,O1, WDIR)
                ENDDO
               ENDDO

!  close control and loops

              ENDDO
             ENDIF
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  End profile atmospheric weighting function stuff

      ENDIF

!  Atmospheric Column Weighting function output
!  ############################################

      IF ( DO_COLUMN_LINEARIZATION ) THEN

!  control point for avoiding weighting function output

        IF ( DO_MVOUT_ONLY ) GO TO 4401

!  local number of azimuths

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_NUSERAZMS = 1
        ELSE
          LOCAL_NUSERAZMS = N_USER_RELAZMS
        ENDIF

!  overall header

        write(RUN,'(/a/a/)') &
           'Atmospheric Column Weighting function output', &
           '--------------------------------------------'

!  start beam loop

        DO IB = 1, NBEAMS
         write(RUN,FMT_REAL) &
       '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
         write(RUN,'(a)')' '

!  start azimuth loop

         DO UA = 1, LOCAL_NUSERAZMS

!  azimuth angle header

          IF ( DO_NO_AZIMUTH ) THEN
            write(RUN,FMT_CHAR) &
             '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
          ELSE
            write(RUN,FMT_REAL) &
       '** RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
          ENDIF

          WRITE(RUN,'(a)')' '
          WRITE(RUN,FMT_INTEGER) &
            'Total number of output levels = ',N_USER_LEVELS
          WRITE(RUN,FMT_INTEGER) &
            'Total number of output angles = ',N_OUT_STREAMS

!  detailed output

          DO IDIR = 1, N_DIRECTIONS
            WDIR = WHICH_DIRECTIONS(IDIR)

!  linearization control

            DO Q = 1, N_TOTALCOLUMN_WFS

!  direction headers

              IF (WDIR .EQ. UPIDX ) THEN
               WRITE(RUN,'(/A/A,A31)') &
           '  --> Upwelling Atmospheric Column Weighting functions', &
           '        * with respect to : ',COLUMNWF_NAMES(Q)
              ELSE IF (WDIR .EQ. DNIDX ) THEN
               WRITE(RUN,'(/A/A,A31)') &
           '  --> Downwelling Atmospheric Column Weighting functions', &
           '        * with respect to : ',COLUMNWF_NAMES(Q)
              ENDIF

!  output loop

              DO O1 = 1, NSTOKES
               write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
               write(RUN,'(a)') 'output angles | output levels ---> '
               WRITE(RUN,'(14x,50(2x,1pe11.4,2x))') &
                 (USER_LEVELS(UT), UT = 1, N_USER_LEVELS)
               DO I = 1, N_OUT_STREAMS
                V = VZA_OFFSETS(IB,I) + UA
                WRITE(RUN,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I), &
               (COLUMNWF(Q,UT,V,O1,WDIR), UT = 1, N_USER_LEVELS)
               ENDDO
              ENDDO

!  close control and loops

             ENDDO
           ENDDO
         ENDDO
        ENDDO

!  Mean-value output
!  -----------------

4401     CONTINUE
        IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

          write(RUN,'(/a/a)') &
        'Mean value atmospheric Column weighting function output', &
        '-------------------------------------------------------'

!  start beam loop

          DO IB = 1, NBEAMS

           write(RUN,FMT_REAL) &
       '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
           write(RUN,'(a)')' '

!  detailed output

           DO IDIR = 1, N_DIRECTIONS
             WDIR = WHICH_DIRECTIONS(IDIR)

!  linearization control

             DO Q = 1, N_TOTALCOLUMN_WFS

!  direction headers

               IF (WDIR .EQ. UPIDX ) THEN
                WRITE(RUN,'(/A/A,I2/)') &
       '  --> Upwelling Mean-I & Flux: Column Weighting functions', &
           '        * Parameter number = ',Q
               ELSE IF (WDIR .EQ. DNIDX ) THEN
                WRITE(RUN,'(/A/A,I2/)') &
       '  --> Downwelling Mean_I & Flux: Column Weighting functions', &
           '        * Parameter number = ',Q
               ENDIF

!  Main output

               DO O1 = 1, NSTOKES
                write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
                 WRITE(RUN,'(a/)') &
            'output level     mean int. WF      flux WF'
                 DO UT = 1, N_USER_LEVELS
                  WRITE(RUN,'(2x,F9.5,3x,2E17.7)') &
                          USER_LEVELS(UT), &
                          MINT_ATMOSWF(Q,0,UT,IB,O1,WDIR), &
                          FLUX_ATMOSWF(Q,0,UT,IB,O1,WDIR)
                 ENDDO
               ENDDO
!  close control and loops

             ENDDO
           ENDDO
          ENDDO
        ENDIF

!  End atmospheric column weighting function stuff

      ENDIF

!  Surface Weighting function output
!  #################################

      IF ( DO_SURFACE_LINEARIZATION ) THEN

!  control point for avoiding weighting function output

        IF ( DO_MVOUT_ONLY ) GO TO 402

!  local number of azimuths

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_NUSERAZMS = 1
        ELSE
          IF ( DO_LAMBERTIAN_SURFACE ) THEN
            LOCAL_NUSERAZMS = 1
          ELSE
            LOCAL_NUSERAZMS = N_USER_RELAZMS
          ENDIF
        ENDIF

!  overall header

        write(RUN,'(/a/a/)')'Surface Weighting function output', &
                            '---------------------------------'

!  start beam loop

        DO IB = 1, NBEAMS

         write(RUN,FMT_REAL) &
       '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
         write(RUN,'(a)')' '

!  start azimuth loop

         DO UA = 1, LOCAL_NUSERAZMS

!  azimuth angle header

          IF ( DO_NO_AZIMUTH ) THEN
           write(RUN,FMT_CHAR) &
             '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
          ELSE
           IF ( DO_LAMBERTIAN_SURFACE ) THEN
            write(RUN,FMT_CHAR) &
         '** Lambertian Surface --> Results AZIMUTH-INDEPENDENT **'
           ELSE
            write(RUN,FMT_REAL) &
       '** RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
           ENDIF
          ENDIF

          WRITE(RUN,'(a)')' '
          WRITE(RUN,FMT_INTEGER) &
            'Total number of output levels = ',N_USER_LEVELS
          WRITE(RUN,FMT_INTEGER) &
            'Total number of output angles = ',N_OUT_STREAMS

!  detailed output

          DO IDIR = 1, N_DIRECTIONS
           WDIR = WHICH_DIRECTIONS(IDIR)

!  direction headers

           IF (WDIR .EQ. UPIDX ) THEN
            WRITE(RUN,'(/A)') &
           '  --> Upwelling Surface Weighting functions'
           ELSE IF (WDIR .EQ. DNIDX ) THEN
            WRITE(RUN,'(/A)') &
           '  --> Downwelling Surface Weighting functions'
           ENDIF

!  for each Surface weighting function

           DO Z = 1, N_SURFACE_WFS

            WRITE(RUN,'(/A,I2)') &
           '  -->  Surface Weighting functions, #: ',Z

!  output loop

            DO O1 = 1, NSTOKES
             write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
             write(RUN,'(a)') 'output angles | output levels ---> '
             WRITE(RUN,'(14x,50(2x,1pe11.4,2x))') &
               (USER_LEVELS(UT), UT = 1, N_USER_LEVELS)
             DO I = 1, N_OUT_STREAMS
              V = VZA_OFFSETS(IB,I) + UA
              WRITE(RUN,'(3x,F9.5,2x,50(1PE15.5))')OUT_ANGLES(I), &
               (SURFACEWF(Z,UT,V,O1,WDIR), UT = 1, N_USER_LEVELS)
             ENDDO
            ENDDO

!  close control and loops

           ENDDO
          ENDDO
         ENDDO
        ENDDO

!  Mean-value output
!  -----------------

 402     continue

        IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

         write(RUN,'(/a/a)') &
                'Mean value surface weighting function output', &
                '--------------------------------------------'

!  start beam loop

         DO IB = 1, NBEAMS

         write(RUN,FMT_REAL) &
       '* * RESULTS FOR SOLAR ZENITH ANGLE (degs)=',SZANGLES(IB)
         write(RUN,'(a)')' '

!  detailed output

          DO IDIR = 1, N_DIRECTIONS
            WDIR = WHICH_DIRECTIONS(IDIR)

!  direction headers

            IF (WDIR .EQ. UPIDX ) THEN
              WRITE(RUN,'(/A/A/)') &
       '  --> Upwelling Mean intensity and Flux Weighting functions' // &
           '      * W.R.T surface variation '
            ELSE IF (WDIR .EQ. DNIDX ) THEN
              WRITE(RUN,'(/A/A/)') &
       '  --> Downwelling Mean intensity and Flux Weighting functions' // &
           '      * W.R.T surface variation '
            ENDIF

!  Loop over kernels

            DO Z = 1, N_SURFACE_WFS

!  output level header

              WRITE(RUN,'(a,I2/a/)') &
            '  --> For Surface WF number :',Z, &
            'output level     mean int. WF      flux WF'

!  output level loop

              DO O1 = 1, NSTOKES
               write(RUN,'(/a,a)') ' %%%%%%% STOKES ',STOKESNAM(O1)
               DO UT = 1, N_USER_LEVELS
                WRITE(RUN,'(2x,F9.5,3x,2E17.7)') &
                  USER_LEVELS(UT), &
                  MINT_SURFACEWF(Z,UT,IB,O1,WDIR), &
                  FLUX_SURFACEWF(Z,UT,IB,O1,WDIR)
               ENDDO
              ENDDO

!  close control and loops

            ENDDO
          ENDDO
         ENDDO
        ENDIF

!  End surface weighting function stuff

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_WRITERESULTS


!

      SUBROUTINE VLIDORT_L_WRITESCEN ( &
        SUNIT, NLAYERS, DO_ATMOS_LINEARIZATION, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER,     &
        L_OMEGA_TOTAL_INPUT, L_DELTAU_VERT_INPUT )

      USE VLIDORT_PARS_m, Only : MAX_ATMOSWFS, MAXLAYERS, FMT_SECTION, FMT_HEADING

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::          SUNIT
      INTEGER, INTENT (IN) ::          NLAYERS
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )

!  local variables

      INTEGER :: N, V

!  set up

!   WF input
!   --------

!  layer variational quantities

      IF ( DO_ATMOS_LINEARIZATION ) THEN

        WRITE(SUNIT,FMT_SECTION) &
           'Atmospheric input for Weighting function calculation'

        WRITE(SUNIT, FMT_HEADING) &
              'Single scattering albedo variations'
        WRITE(SUNIT,'(a/)') &
           'Layer varying | parameter numbers-->'
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
             WRITE(SUNIT,'(I3,14x,5(1pe12.4))') &
              N,(L_OMEGA_TOTAL_INPUT(V,N),V=1,LAYER_VARY_NUMBER(N))
          ENDIF
        ENDDO

        WRITE(SUNIT, FMT_HEADING) 'Output level variations'
        WRITE(SUNIT,'(a/)') &
           'Layer varying | parameter numbers-->'
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            WRITE(SUNIT,'(I3,14x,5(1pe12.4))') &
               N,(L_DELTAU_VERT_INPUT(V,N),V=1,LAYER_VARY_NUMBER(N))
          ENDIF
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_WRITESCEN

!

      SUBROUTINE VLIDORT_WRITE_LIN_INPUT ( &
        NSTOKES, NLAYERS, NSTREAMS, N_SZANGLES, N_USER_RELAZMS, &
        N_USER_VZANGLES, DO_OBSERVATION_GEOMETRY, DO_SIMULATION_ONLY,       &
        N_TOTALCOLUMN_WFS, N_TOTALPROFILE_WFS, N_SURFACE_WFS, N_SLEAVE_WFS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                                 &
        COLUMNWF_NAMES, PROFILEWF_NAMES, DO_ATMOS_LBBF, DO_SURFACE_LBBF,    &
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,   &
        L_FMATRIX_UP, L_FMATRIX_DN,                                         &
        DO_COLUMN_LINEARIZATION, DO_PROFILE_LINEARIZATION, DO_ATMOS_LINEARIZATION,&
        DO_SURFACE_LINEARIZATION, DO_LINEARIZATION, DO_SLEAVE_WFS )

!  3/28/14. Changes for Version 2.7. remove LTE linearization references. Add LBBF

!  4/15/20. Version 2.8.2. Argument list rearranged/renamed, follows LIDORT Version 3.8.2
!  4/15/20. Version 2.8.2. MAXMOMENTS instead of MAXMOMENTS_INPUT, NGREEK_MOMENTS_INPUT dropped

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAX_GEOMETRIES, MAXSTOKES_SQ, MAX_ATMOSWFS, &
                                 DWFL, DWFL1, DWFI, DWFI1, DWFR, DWFR2, DWFR4

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   NLAYERS
      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   N_SZANGLES
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   N_USER_VZANGLES

      LOGICAL, INTENT(IN) ::   DO_OBSERVATION_GEOMETRY

!  -------------------------
!  Linearized Inputs - Fixed
!  -------------------------

      LOGICAL, INTENT(IN) ::            DO_SIMULATION_ONLY

      INTEGER, INTENT(IN) ::            N_TOTALCOLUMN_WFS
      INTEGER, INTENT(IN) ::            N_TOTALPROFILE_WFS
      INTEGER, INTENT(IN) ::            N_SURFACE_WFS
      INTEGER, INTENT(IN) ::            N_SLEAVE_WFS

      LOGICAL, INTENT(IN) ::            LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT(IN) ::            LAYER_VARY_NUMBER ( MAXLAYERS )

      CHARACTER (LEN=31), INTENT(IN) :: COLUMNWF_NAMES  ( MAX_ATMOSWFS )
      CHARACTER (LEN=31), INTENT(IN) :: PROFILEWF_NAMES ( MAX_ATMOSWFS )

!  New Version 2.7. LTE and SURFBB linearization flags replaced

!      LOGICAL, INTENT(IN) ::            DO_LTE_LINEARIZATION, DO_SURFBB_LINEARIZATION
      LOGICAL, INTENT(IN) ::            DO_ATMOS_LBBF,DO_SURFACE_LBBF

!  4/15/20. Version 2.8.2. MAXMOMENTS instead of MAXMOMENTS_INPUT

      DOUBLE PRECISION, INTENT(IN) ::   L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT(IN) ::   L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT(IN) ::   L_FMATRIX_UP ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES,  6 )
      DOUBLE PRECISION, INTENT(IN) ::   L_FMATRIX_DN ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES,  6 )

!  ----------------------------
!  Linearized Inputs - Variable
!  ----------------------------

      LOGICAL, INTENT(IN) ::         DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT(IN) ::         DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT(IN) ::         DO_ATMOS_LINEARIZATION

      LOGICAL, INTENT(IN) ::         DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT(IN) ::         DO_LINEARIZATION

      LOGICAL, INTENT(IN) ::         DO_SLEAVE_WFS

!  Local variables

      INTEGER :: OUTUNIT, NGEOMS
      INTEGER :: GEO,K,LAY,MOM,NGREEK,SS,WF,GMASK(8)
      data GMASK / 1, 2, 5, 6, 11, 12, 15, 16 /

      INTEGER :: N_WFS
      CHARACTER (LEN=9) :: DWFC1S

!  Define local variables

      NGEOMS = N_SZANGLES
      IF ( .NOT.DO_OBSERVATION_GEOMETRY ) NGEOMS = N_SZANGLES * N_USER_VZANGLES * N_USER_RELAZMS

!  Open output file

      OUTUNIT = 111
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_LIN_INPUT.dbg',status = 'replace')

!  Define local variables
!  Changed 16 oct 14, only output non-zero entries of GreekMat

      NGREEK = 1
      IF ( NSTOKES.eq.3 ) NGREEK = 5
      IF ( NSTOKES.eq.4 ) NGREEK = 8

      N_WFS = 0
      IF (DO_COLUMN_LINEARIZATION) THEN
        N_WFS = N_TOTALCOLUMN_WFS
      ELSE IF (DO_PROFILE_LINEARIZATION) THEN
        N_WFS = N_TOTALPROFILE_WFS
      END IF

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized Inputs - Fixed'
      WRITE(OUTUNIT,'(A)') '-------------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_SIMULATION_ONLY      = ',DO_SIMULATION_ONLY
!      WRITE(OUTUNIT,DWFL)  'DO_SURFBB_LINEARIZATION = ',DO_SURFBB_LINEARIZATION

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_TOTALCOLUMN_WFS  = ',N_TOTALCOLUMN_WFS
      WRITE(OUTUNIT,DWFI)  'N_TOTALPROFILE_WFS = ',N_TOTALPROFILE_WFS
      WRITE(OUTUNIT,DWFI)  'N_SURFACE_WFS      = ',N_SURFACE_WFS
      WRITE(OUTUNIT,DWFI)  'N_SLEAVE_WFS       = ',N_SLEAVE_WFS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFL1)  'LAY = ',LAY,&
          ' LAYER_VARY_FLAG(LAY)   = ',LAYER_VARY_FLAG(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFI1)  'LAY = ',LAY,&
          ' LAYER_VARY_NUMBER(LAY) = ',LAYER_VARY_NUMBER(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DWFC1S = '(A,I3,3A)'
      IF (DO_COLUMN_LINEARIZATION) THEN
        DO WF=1,N_TOTALCOLUMN_WFS
          WRITE(OUTUNIT,DWFC1S)  'WF = ',WF,&
            ' COLUMNWF_NAMES(WF) = |',COLUMNWF_NAMES(WF),'|'
        END DO
      ELSE IF (DO_PROFILE_LINEARIZATION) THEN
        DO WF=1,N_TOTALPROFILE_WFS
          WRITE(OUTUNIT,DWFC1S)  'WF = ',WF,&
            ' PROFILEWF_NAMES(WF) = |',PROFILEWF_NAMES(WF),'|'
        END DO
      END IF

!  replaced
!      WRITE(OUTUNIT,*)
!      WRITE(OUTUNIT,DWFL)  'DO_LTE_LINEARIZATION    = ',DO_LTE_LINEARIZATION

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  'LAY = ',LAY,' WF = ',WF,&
            ' L_DELTAU_VERT_INPUT(WF,LAY) = ',L_DELTAU_VERT_INPUT(WF,LAY)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  'LAY = ',LAY,' WF = ',WF,&
            ' L_OMEGA_TOTAL_INPUT(WF,LAY) = ',L_OMEGA_TOTAL_INPUT(WF,LAY)
        END DO
      END DO

 !  4/15/20. Version 2.8.2. Only have 2*NSTREAMS as the number of usable moments

      WRITE(OUTUNIT,*)
      DO K=1,NGREEK
        SS = GMASK(K)
        DO LAY=1,NLAYERS
          DO MOM=0,2*NSTREAMS
            DO WF=1,N_WFS
              WRITE(OUTUNIT,DWFR4) &
                'SS = ',SS,' LAY = ',LAY,' MOM = ',MOM,' WF = ',WF,&
                ' L_GREEKMAT_TOTAL_INPUT(WF,MOM,LAY,SS) = ',&
                  L_GREEKMAT_TOTAL_INPUT(WF,MOM,LAY,SS)
            END DO
          END DO
        END DO
      END DO

!  New section for 2.8 code. F-matrices

      WRITE(OUTUNIT,*)
      DO K=1,6
        DO GEO=1,NGEOMS
          DO LAY=1,NLAYERS
            DO WF=1,N_WFS
              WRITE(OUTUNIT,DWFR4) 'K = ',K,' GEO = ',GEO,' LAY = ',LAY,' WF = ',WF, &
                ' L_FMATRIX_UP(WF,LAY,GEO,K) = ',L_FMATRIX_UP(WF,LAY,GEO,K)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO K=1,6
        DO GEO=1,NGEOMS
          DO LAY=1,NLAYERS
            DO WF=1,N_WFS
              WRITE(OUTUNIT,DWFR4) 'K = ',K,' GEO = ',GEO,' LAY = ',LAY,' WF = ',WF, &
                ' L_FMATRIX_DN(WF,LAY,GEO,K) = ',L_FMATRIX_DN(WF,LAY,GEO,K)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '----------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized Inputs - Variable'
      WRITE(OUTUNIT,'(A)') '----------------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_COLUMN_LINEARIZATION  = ',DO_COLUMN_LINEARIZATION
      WRITE(OUTUNIT,DWFL)  'DO_PROFILE_LINEARIZATION = ',DO_PROFILE_LINEARIZATION
      WRITE(OUTUNIT,DWFL)  'DO_ATMOS_LINEARIZATION   = ',DO_ATMOS_LINEARIZATION
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_LINEARIZATION = ',DO_SURFACE_LINEARIZATION
      WRITE(OUTUNIT,DWFL)  'DO_LINEARIZATION         = ',DO_LINEARIZATION
      WRITE(OUTUNIT,DWFL)  'DO_SLEAVE_WFS            = ',DO_SLEAVE_WFS

!  3/28/14. Changes for Version 2.7. remove LTE linearization. Add LBBF

      WRITE(OUTUNIT,DWFL)  'DO_ATMOS_LBBF            = ',DO_ATMOS_LBBF
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_LBBF          = ',DO_SURFACE_LBBF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_LIN_INPUT

!

      SUBROUTINE VLIDORT_WRITE_LIN_SUP_BRDF_INPUT ( &
        DO_USER_STREAMS, DO_SURFACE_EMISSION, &
        NSTOKES,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,&
        NSTREAMS,N_SURFACE_WFS,&
        LS_EXACTDB_BRDFUNC,LS_BRDF_F_0,LS_BRDF_F,&
        LS_USER_BRDF_F_0,LS_USER_BRDF_F,&
        LS_EMISSIVITY,LS_USER_EMISSIVITY)

!mick mod 9/19/2017 - added DO_USER_STREAMS & DO_SURFACE_EMISSION for better output control
!                   - added NMOMENTS for VLIDORT internal consistency

!  4/15/20. Version 2.8.2. Remove NMOMENTS and MAXMOMENTS, Only printing First Fourier component.

      USE VLIDORT_PARS_m, Only : MAX_SURFACEWFS, MAXSTOKES, MAXSTOKES_SQ, MAXSTREAMS, &
                                 MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS, DWFR3, DWFR5

      IMPLICIT NONE

!  Input

      LOGICAL, INTENT(IN) ::   DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::   DO_SURFACE_EMISSION

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   N_SZANGLES
      INTEGER, INTENT(IN) ::   N_USER_VZANGLES
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   N_SURFACE_WFS

      DOUBLE PRECISION, INTENT(IN) ::   LS_EXACTDB_BRDFUNC &
                    ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  4/15/20. Version 2.8.2.  Only have the M = 0 Component to write now.

      DOUBLE PRECISION, INTENT(IN) ::   LS_BRDF_F_0      ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   LS_BRDF_F        ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT(IN) ::   LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   LS_USER_BRDF_F   ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

      DOUBLE PRECISION, INTENT(IN) ::   LS_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT(IN) ::   LS_USER_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTOKES, MAX_USER_STREAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SS,MOM,NREFL,K,STRM,S,USTRM,IB,URA,&
                 STRMI,STRMJ,SWF, RMASK3(9), RMASK4(16), RMASK(16)

      data RMASK3 / 1, 2, 3, 5, 6, 7, 9, 10, 11 /
      data RMASK4 / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 /

!  Open output file

      OUTUNIT = 112
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_LIN_SUP_BRDF_INPUT.dbg',&
            status = 'replace')

!  Define local variables
!   Changed 16 October 2014, according to correct entries in 4x4 BRDF matrices

      NREFL = NSTOKES * NSTOKES ; RMASK(1) = 1
      IF ( NSTOKES.eq.3 ) RMASK(1:NREFL) = RMASK3(1:NREFL)
      IF ( NSTOKES.eq.4 ) RMASK(1:NREFL) = RMASK4(1:NREFL)

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '---------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '---------------------------------'

      WRITE(OUTUNIT,*)
      DO K=1,NREFL
        SS = RMASK(K)
        DO IB=1,N_SZANGLES
          DO URA=1,N_USER_RELAZMS
            DO USTRM=1,N_USER_VZANGLES
              DO SWF=1,N_SURFACE_WFS
                WRITE(OUTUNIT,DWFR5) &
                  'SS = ',SS,' IB = ',IB,' URA = ',URA,&
                  ' USTRM = ',USTRM,' SWF = ',SWF,&
                  ' LS_EXACTDB_BRDFUNC(SWF,SS,USTRM,URA,IB) = ',&
                    LS_EXACTDB_BRDFUNC(SWF,SS,USTRM,URA,IB)
              END DO
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
        DO IB=1,N_SZANGLES ; DO STRM=1,NSTREAMS ; DO SWF=1,N_SURFACE_WFS
           WRITE(OUTUNIT,DWFR5) &
                  'SS = ',SS,' IB = ',IB,' STRM = ',STRM,' MOM = ',MOM,' SWF = ',SWF,&
                  ' LS_BRDF_F_0(SWF,MOM,SS,STRM,IB) = ',LS_BRDF_F_0(SWF,SS,STRM,IB)
        END DO ; END DO ; END DO
      END DO

      WRITE(OUTUNIT,*)
      DO K=1,NREFL
        SS = RMASK(K)
        DO STRMJ=1,NSTREAMS ; DO STRMI=1,NSTREAMS ; DO SWF=1,N_SURFACE_WFS
           WRITE(OUTUNIT,DWFR5) &
                  'SS = ',SS,' STRMJ = ',STRMJ,' STRMI = ',STRMI,' MOM = ',MOM,' SWF = ',SWF,&
                  ' LS_BRDF_F(SWF,MOM,SS,STRMI,STRMJ) = ',LS_BRDF_F(SWF,SS,STRMI,STRMJ)
        END DO ; END DO ; END DO
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO K=1,NREFL
          SS = RMASK(K)
          DO IB=1,N_SZANGLES ; DO USTRM=1,N_USER_VZANGLES ; DO SWF=1,N_SURFACE_WFS
             WRITE(OUTUNIT,DWFR5) &
                    'SS = ',SS,' IB = ',IB,' USTRM = ',USTRM,' MOM = ',MOM,' SWF = ',SWF,&
                    ' LS_USER_BRDF_F_0(SWF,MOM,SS,USTRM,IB) = ',LS_USER_BRDF_F_0(SWF,SS,USTRM,IB)
          END DO ; END DO ; END DO
        END DO

        WRITE(OUTUNIT,*)
        DO K=1,NREFL
          SS = RMASK(K)
          DO STRM=1,NSTREAMS ; DO USTRM=1,N_USER_VZANGLES ; DO SWF=1,N_SURFACE_WFS
             WRITE(OUTUNIT,DWFR5) &
                    'SS = ',SS,' STRM = ',STRM,' USTRM = ',USTRM,' MOM = ',MOM,' SWF = ',SWF,&
                    ' LS_USER_BRDF_F(SWF,MOM,SS,USTRM,STRM) = ',LS_USER_BRDF_F(SWF,SS,USTRM,STRM)
          END DO ; END DO ; END DO
        END DO
      END IF

      IF ( DO_SURFACE_EMISSION ) THEN
        WRITE(OUTUNIT,*)
        DO S=1,NSTOKES
          DO STRM=1,NSTREAMS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR3)  'S = ',S,' STRM = ',STRM,' SWF = ',SWF,&
                ' LS_EMISSIVITY(SWF,S,STRM) = ',LS_EMISSIVITY(SWF,S,STRM)
            END DO
          END DO
        END DO

        IF ( DO_USER_STREAMS ) THEN
          WRITE(OUTUNIT,*)
          DO S=1,NSTOKES
            DO USTRM=1,N_USER_VZANGLES
              DO SWF=1,N_SURFACE_WFS
                WRITE(OUTUNIT,DWFR3)  'S = ',S,' USTRM = ',USTRM,' SWF = ',SWF,&
                  ' LS_USER_EMISSIVITY(SWF,S,USTRM) = ',&
                    LS_USER_EMISSIVITY(SWF,S,USTRM)
              END DO
            END DO
          END DO
        END IF
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_LIN_SUP_BRDF_INPUT

!

      SUBROUTINE VLIDORT_WRITE_LIN_SUP_SS_INPUT ( &
        DO_COLUMN_LINEARIZATION,DO_PROFILE_LINEARIZATION,&
        DO_SURFACE_LINEARIZATION,&
        NSTOKES,NLAYERS,N_USER_LEVELS,&
        N_TOTALCOLUMN_WFS,N_TOTALPROFILE_WFS,N_TOTALSURFACE_WFS,&
        COLUMNWF_SS,COLUMNWF_DB,PROFILEWF_SS,PROFILEWF_DB,SURFACEWF_DB)

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAX_ATMOSWFS, MAX_SURFACEWFS, MAXLAYERS, MAX_USER_LEVELS,  &
                                 MAX_GEOMETRIES, MAX_DIRECTIONS, DWFR4, DWFR5, DWFR6

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT(IN) :: DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT(IN) :: DO_SURFACE_LINEARIZATION

      INTEGER, INTENT(IN) :: NSTOKES
      INTEGER, INTENT(IN) :: NLAYERS
      INTEGER, INTENT(IN) :: N_USER_LEVELS
      INTEGER, INTENT(IN) :: N_TOTALCOLUMN_WFS
      INTEGER, INTENT(IN) :: N_TOTALPROFILE_WFS
      INTEGER, INTENT(IN) :: N_TOTALSURFACE_WFS

      DOUBLE PRECISION, INTENT(IN) :: COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: COLUMNWF_DB &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION, INTENT(IN) :: PROFILEWF_SS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: PROFILEWF_DB &
          ( MAX_ATMOSWFS, MAXLAYERS, &
            MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION, INTENT(IN) :: SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: DIR,GEO,S,ULEV,LAY,WF,SWF

!  Open output file

      OUTUNIT = 113
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_LIN_SUP_SS_INPUT.dbg',&
            status = 'replace')

!  Write all single-scatter input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized SS Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '-------------------------------'

      IF (DO_COLUMN_LINEARIZATION) THEN
        WRITE(OUTUNIT,*)
        DO DIR=1,MAX_DIRECTIONS
          DO S=1,NSTOKES
            DO GEO=1,MAX_GEOMETRIES
              DO ULEV=1,N_USER_LEVELS
                DO WF=1,N_TOTALCOLUMN_WFS
                  WRITE(OUTUNIT,DWFR5) &
                    'DIR = ',DIR,' S = ',S,' GEO = ',GEO,&
                    ' ULEV = ',ULEV,' WF = ',WF,&
                    ' COLUMNWF_SS(WF,ULEV,GEO,S,DIR) = ',&
                      COLUMNWF_SS(WF,ULEV,GEO,S,DIR)
                END DO
              END DO
            END DO
          END DO
        END DO

        WRITE(OUTUNIT,*)
        DO S=1,NSTOKES
          DO GEO=1,MAX_GEOMETRIES
            DO ULEV=1,N_USER_LEVELS
              DO WF=1,N_TOTALCOLUMN_WFS
                WRITE(OUTUNIT,DWFR4) &
                  'S = ',S,' GEO = ',GEO,' ULEV = ',ULEV,' WF = ',WF,&
                  ' COLUMNWF_DB(WF,ULEV,GEO,S) = ',&
                    COLUMNWF_DB(WF,ULEV,GEO,S)
              END DO
            END DO
          END DO
        END DO
      END IF

      IF (DO_PROFILE_LINEARIZATION) THEN
        WRITE(OUTUNIT,*)
        DO DIR=1,MAX_DIRECTIONS
          DO S=1,NSTOKES
            DO GEO=1,MAX_GEOMETRIES
              DO ULEV=1,N_USER_LEVELS
                DO LAY=1,NLAYERS
                  DO WF=1,N_TOTALPROFILE_WFS
                    WRITE(OUTUNIT,DWFR6) &
                      'DIR = ',DIR,' S = ',S,' GEO = ',GEO,&
                      ' ULEV = ',ULEV,' LAY = ',LAY,' WF = ',WF,&
                      ' PROFILEWF_SS(WF,LAY,ULEV,GEO,S,DIR) = ',&
                        PROFILEWF_SS(WF,LAY,ULEV,GEO,S,DIR)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO

        WRITE(OUTUNIT,*)
        DO S=1,NSTOKES
          DO GEO=1,MAX_GEOMETRIES
            DO ULEV=1,N_USER_LEVELS
              DO LAY=1,NLAYERS
                DO WF=1,N_TOTALPROFILE_WFS
                  WRITE(OUTUNIT,DWFR5) &
                    'S = ',S,' GEO = ',GEO,' ULEV = ',ULEV,&
                    ' LAY = ',LAY,' WF = ',WF,&
                    ' PROFILEWF_DB(WF,LAY,ULEV,GEO,S) = ',&
                      PROFILEWF_DB(WF,LAY,ULEV,GEO,S)
                END DO
              END DO
            END DO
          END DO
        END DO
      END IF

      IF (DO_SURFACE_LINEARIZATION) THEN
        WRITE(OUTUNIT,*)
        DO S=1,NSTOKES
          DO GEO=1,MAX_GEOMETRIES
            DO ULEV=1,N_USER_LEVELS
              DO SWF=1,N_TOTALSURFACE_WFS
                WRITE(OUTUNIT,DWFR4) &
                  'S = ',S,' GEO = ',GEO,' ULEV = ',ULEV,' SWF = ',SWF,&
                  ' SURFACEWF_DB(SWF,ULEV,GEO,S) = ',&
                    SURFACEWF_DB(SWF,ULEV,GEO,S)
              END DO
            END DO
          END DO
        END DO
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_LIN_SUP_SS_INPUT

!

      SUBROUTINE VLIDORT_WRITE_LIN_SUP_SLEAVE_INPUT ( &
        DO_USER_STREAMS,NSTOKES,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,NSTREAMS,N_SLEAVE_WFS,&
        LSSL_SLTERM_ISOTROPIC,LSSL_SLTERM_USERANGLES,LSSL_SLTERM_F_0,LSSL_USER_SLTERM_F_0)

!mick mod 9/19/2017 - added DO_USER_STREAMS for better output control
!                   - added NMOMENTS for VLIDORT internal consistency

!  4/15/20. Version 2.8.2. Remove NMOMS and MAXMOMENTS, Only printing First Fourier componnet.

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAX_SLEAVEWFS, MAX_USER_STREAMS,  &
                                 MAXSTREAMS, MAXBEAMS, MAX_USER_RELAZMS, DWFR3, DWFR5

      IMPLICIT NONE

      LOGICAL, INTENT(IN) ::   DO_USER_STREAMS

      INTEGER, INTENT(IN) ::   NSTOKES
      INTEGER, INTENT(IN) ::   N_SZANGLES
      INTEGER, INTENT(IN) ::   N_USER_VZANGLES
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   N_SLEAVE_WFS

!  4/15/20. Version 2.8.2. Drop MAXMOMENTS dimension from Fourier terms.

      DOUBLE PRECISION, INTENT(IN) ::   LSSL_SLTERM_ISOTROPIC ( MAX_SLEAVEWFS, MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   LSSL_SLTERM_USERANGLES &
                  ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )
      DOUBLE PRECISION, INTENT(IN) ::   LSSL_SLTERM_F_0    ( MAX_SLEAVEWFS,  MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   LSSL_USER_SLTERM_F_0 ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SWF,MOM,STRM,S,USTRM,IB,URA

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'VLIDORT_WRITE_LIN_SUP_SLEAVE_INPUT.dbg',status = 'replace')

!  Write all surface-leaving input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized SLEAVE Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '-----------------------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,N_SZANGLES
        DO S=1,NSTOKES
          DO SWF=1,N_SLEAVE_WFS
            WRITE(OUTUNIT,DWFR3) &
              'IB = ',IB,' S = ',S,' SWF = ',SWF,&
              ' LSSL_SLTERM_ISOTROPIC(SWF,S,IB) = ',&
                LSSL_SLTERM_ISOTROPIC(SWF,S,IB)
          END DO
        END DO
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO IB=1,N_SZANGLES
          DO URA=1,N_USER_RELAZMS
            DO USTRM=1,N_USER_VZANGLES
              DO S=1,NSTOKES
                DO SWF=1,N_SLEAVE_WFS
                  WRITE(OUTUNIT,DWFR5) &
                    'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,&
                    ' S = ',S,' SWF = ',SWF,&
                    ' LSSL_SLTERM_USERANGLES(SWF,S,USTRM,URA,IB) = ',&
                      LSSL_SLTERM_USERANGLES(SWF,S,USTRM,URA,IB)
                END DO
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
        DO STRM=1,NSTREAMS ; DO S=1,NSTOKES ; DO SWF=1,N_SLEAVE_WFS
           WRITE(OUTUNIT,DWFR5) &
                  'IB = ',IB,' STRM = ',STRM,' S = ',S,' MOM = ',MOM,' SWF = ',SWF,&
                  ' LSSL_SLTERM_F_0(SWF,MOM,S,STRM,IB) = ',LSSL_SLTERM_F_0(SWF,S,STRM,IB)
        END DO ; END DO ; END DO
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO IB=1,N_SZANGLES
          DO USTRM=1,N_USER_VZANGLES ; DO S=1,NSTOKES ; DO SWF=1,N_SLEAVE_WFS
             WRITE(OUTUNIT,DWFR5) &
                    'IB = ',IB,' USTRM = ',USTRM,' S = ',S,' MOM = ',MOM,' SWF = ',SWF,&
                    ' LSSL_USER_SLTERM_F_0(SWF,MOM,S,USTRM,IB) = ',LSSL_USER_SLTERM_F_0(SWF,S,USTRM,IB)
          END DO ; END DO ; END DO
        END DO
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VLIDORT_WRITE_LIN_SUP_SLEAVE_INPUT

      END MODULE vlidort_l_writemodules_m
