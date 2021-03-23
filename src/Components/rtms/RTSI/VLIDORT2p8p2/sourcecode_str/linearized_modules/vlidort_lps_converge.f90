
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

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #       VLIDORT_LPS_CONVERGE                             #
! #       VLIDORT_LPS_CONVERGE_OBSGEO                      #
! #       VLIDORT_LPS_CONVERGE_DOUBLET                     #
! #                                                        #
! ##########################################################

!  4/15/20. Version 2.8.2. Separate module for the converge routines
!          ==> Uses type structures directly.
!          ==> Addition of new VLIDORT_LPS_CONVERGE_DOUBLET subroutine


      MODULE vlidort_lps_converge_m

!  Dependencies

      USE VLIDORT_LinOutputs_def_m
      USE VLIDORT_LinSup_SS_def_m

!  Everything public

      PUBLIC  :: VLIDORT_LPS_CONVERGE, &
                 VLIDORT_LPS_CONVERGE_OBSGEO, &
                 VLIDORT_LPS_CONVERGE_DOUBLET


      CONTAINS

      SUBROUTINE VLIDORT_LPS_CONVERGE ( &
        DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH, DO_DBCORRECTION,        & ! Input flags
        DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_LAMBERTIAN_SURFACE,       & ! Input flags
        NSTOKES, NLAYERS, N_USER_STREAMS, LOCAL_UM_START, N_OUT_STREAMS, N_USER_LEVELS,  & ! Input control numbers
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,                 & ! Input linearization control
        IBEAM, FOURIER, VZA_OFFSETS, N_DIRECTIONS, WHICH_DIRECTIONS,                     & ! Input bookkeeping
        LOCAL_N_USERAZM, AZMFAC, PROFILEWF_F, SURFACEWF_F, VLIDORT_LinSS, VLIDORT_LinOut ) ! Input Azm, fields

!  Just upgrades the weighting function Fourier cosine series
!   Version 2.8, 3/1/17. Logic for FOCORR variables changed

!  4/15/20. Version 2.8.2
!   -- Replaced SURFACE_DB/PROFILEWF_SS/DB (FO inputs), PROFILEWF/SURFACEWF (final outputs) with Type structure variables
!   -- Rearranged argument list, removed DO_FOCORR_EXTERNAL

!  module, dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAX_GEOMETRIES, MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS, &
                                 MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS,  &
                                 MAX_USER_LEVELS, MAX_ATMOSWFS, MAX_SURFACEWFS, ZERO, UPIDX

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_FOCORR
      LOGICAL, INTENT (IN) ::           DO_FOCORR_ALONE
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH

!  More flags

      LOGICAL, INTENT (IN) ::           DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE

!  Control integers

      INTEGER, INTENT (IN) ::           NSTOKES, NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      INTEGER, INTENT (IN) ::           N_OUT_STREAMS

!  Linearization control

      LOGICAL, INTENT (IN) ::           LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           N_SURFACE_WFS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

!  Fourier component and beam index

      INTEGER, INTENT (IN) ::           FOURIER, IBEAM

!  Bookkeeping: Offsets for geometry indexing

      INTEGER, INTENT (IN) ::           VZA_OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )

!  Directional control

      INTEGER, INTENT (IN) ::           N_DIRECTIONS
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Local number of azimuths and azimuth factors

      INTEGER, INTENT (IN) ::           LOCAL_N_USERAZM
      DOUBLE PRECISION, INTENT (IN) ::  AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )

!  Fourier-component Profile/Surface weighting functions at user angles

      DOUBLE PRECISION, INTENT (IN) :: SURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
                                                     MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                     MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Type structure for Single scatter and Direct-beam weighting functions; We want Prof/Surf substructures

      TYPE(VLIDORT_LinSup_SS), intent(in)     :: VLIDORT_LinSS

!  modified/output variables
!  -------------------------

!  Full Jacobian results, type definition; We want Prof/Surf substructures

      TYPE(VLIDORT_LinOutputs), intent(inout)      :: VLIDORT_LinOut

!  local variables
!  ---------------

      INTEGER :: I, IDIR, UT, UA, Q, N, W, Z, V, O1

!   For Single scatter corrections, the quadrature output flag is
!    turned off, so that N_OUT_STREAMS = N_USER_STREAMS, and the
!    offsetting is consistent here with the usage elsewhere.

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  single scatter calculation is initialized to zero here (Version 2.7)

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_OUT_STREAMS
                        DO UA = 1, LOCAL_N_USERAZM
                          V = VZA_OFFSETS(IBEAM,I) + UA
                          DO O1 = 1, NSTOKES
                            VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) = PROFILEWF_F(Q,N,UT,I,IBEAM,O1,W)
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ELSE
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_OUT_STREAMS
                        DO UA = 1, LOCAL_N_USERAZM
                          V = VZA_OFFSETS(IBEAM,I) + UA
                          DO O1 = 1, NSTOKES
                            VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) = ZERO
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!    Version 2.8 some renaming of variables
          !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!          IF ( DO_SSFULL .OR. &
!               ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

          IF ( DO_FOCORR ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = LOCAL_UM_START, N_OUT_STREAMS
                        DO UA = 1, LOCAL_N_USERAZM
                          V = VZA_OFFSETS(IBEAM,I) + UA
                          DO O1 = 1, NSTOKES
                            !CALL TP23E (FOURIER,Q,N,UT,V,O1,W,PROFILEWF,PROFILEWF_SS)
                            VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) = &
                              VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) + VLIDORT_LinSS%Prof%TS_PROFILEWF_SS(Q,N,UT,V,O1,W)
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Add the direct bounce to the upwelling if flagged
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!       Version 2.8    Changed Logic for SS terms

          !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
!          IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN
          IF ( (DO_FOCORR.OR.DO_DBCORRECTION) .AND.DO_UPWELLING ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO UT = 1, N_USER_LEVELS
                    DO I = LOCAL_UM_START, N_OUT_STREAMS
                      DO UA = 1, LOCAL_N_USERAZM
                        V = VZA_OFFSETS(IBEAM,I) + UA
                        DO O1 = 1, NSTOKES
                          !CALL TP23F (FOURIER,Q,N,UT,V,O1,PROFILEWF,PROFILEWF_DB)
                          VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,UPIDX) = &
                  VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,UPIDX) + VLIDORT_LinSS%Prof%TS_PROFILEWF_DB(Q,N,UT,V,O1)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  end atmospheric WF clause

        ENDIF

!  Reflectance weighting functions
!  -------------------------------

!  Diffuse field at all output angles
!  Alternative - zero the output (single scatter case)

        IF ( DO_SURFACE_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
!mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
            DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = LOCAL_UM_START, N_USER_STREAMS
                    DO UA = 1, LOCAL_N_USERAZM
                      V = VZA_OFFSETS(IBEAM,I) + UA
                      DO O1 = 1, NSTOKES
                        VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,W) = SURFACEWF_F(Z,UT,I,IBEAM,O1,W)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
!mick fix 9/12/2012 - added N_SLEAVE_WFS for Z
            DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = LOCAL_UM_START, N_USER_STREAMS
                    DO UA = 1, LOCAL_N_USERAZM
                      V = VZA_OFFSETS(IBEAM,I) + UA
                      DO O1 = 1, NSTOKES
                        VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,W) = ZERO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 2.3
!     Version 2.8    Changed Logic for SS terms

          !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
          !IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL)  .AND.DO_UPWELLING ) THEN

          IF ( (DO_FOCORR.OR.DO_DBCORRECTION) .AND.DO_UPWELLING ) THEN
!mick fix 9/12/2012 - added N_SLEAVE_WFS for Z
            DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO UT = 1, N_USER_LEVELS
                DO I = LOCAL_UM_START, N_OUT_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = VZA_OFFSETS(IBEAM,I) + UA
                    DO O1 = 1, NSTOKES
                      !CALL TP26F (FOURIER,Z,UT,V,O1,SURFACEWF,SURFACEWF_DB)
                      VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,UPIDX) = &
                  VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,UPIDX) + VLIDORT_LinSS%Surf%TS_SURFACEWF_DB(Z,UT,V,O1)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  End surface linearization WF clause

        ENDIF

!  If no_azimuth, then exit
!  ------------------------

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Add next Fourier component to output

          DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                DO UA = 1, LOCAL_N_USERAZM
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_OUT_STREAMS
                        V = VZA_OFFSETS(IBEAM,I) + UA
                        DO O1 = 1, NSTOKES
                          !CALL TP23G (FOURIER,Q,N,UT,I,IBEAM,V,O1,W,PROFILEWF,PROFILEWF_F)
                          VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) = &
                     VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) + PROFILEWF_F(Q,N,UT,I,IBEAM,O1,W)*AZMFAC(I,IBEAM,UA,O1)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO

!  End profile atmospheric WF clause

        ENDIF

!  albedo weighting functions (non-Lambertian only)
!  ------------------------------------------------

        IF ( DO_SURFACE_LINEARIZATION ) THEN
          IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN

!  Copy full output

!mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
            DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO UA = 1, LOCAL_N_USERAZM
                DO IDIR = 1, N_DIRECTIONS
                  W = WHICH_DIRECTIONS(IDIR)
                  DO UT = 1, N_USER_LEVELS
                    DO I = 1, N_OUT_STREAMS
                      V = VZA_OFFSETS(IBEAM,I) + UA
                      DO O1 = 1, NSTOKES
                        !CALL TP26G (FOURIER,Z,UT,I,IBEAM,V,O1,W,SURFACEWF,SURFACEWF_F)
                        VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,W) = &
                 VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,W) + SURFACEWF_F(Z,UT,I,IBEAM,O1,W)*AZMFAC(I,IBEAM,UA,O1)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

!  End surface linearization clause

          ENDIF
        ENDIF

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LPS_CONVERGE

!

      SUBROUTINE VLIDORT_LPS_CONVERGE_DOUBLET ( &
        DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH, DO_DBCORRECTION,      & ! Input flags
        DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_LAMBERTIAN_SURFACE,     & ! Input flags
        NSTOKES, NLAYERS, N_USER_STREAMS, N_OUT_STREAMS, N_USER_LEVELS,                & ! Input control numbers
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,               & ! Input linearization control
        IBEAM, FOURIER, SZD_OFFSETS, N_DIRECTIONS, WHICH_DIRECTIONS,                   & ! Input bookkeeping
        AZMFAC, PROFILEWF_F, SURFACEWF_F, VLIDORT_LinSS, VLIDORT_LinOut )                ! Input Azm, fields

!  Just upgrades the weighting function Fourier cosine series
!   Version 2.8, 3/1/17. Logic for FOCORR variables changed

!  4/15/20. Version 2.8.2. Brand New subroutine for the doublet convergence option.

      USE VLIDORT_PARS_m, Only : MAX_GEOMETRIES, MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS, &
                                 MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS,  &
                                 MAX_USER_LEVELS, MAX_ATMOSWFS, MAX_SURFACEWFS, ZERO, UPIDX

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_FOCORR
      LOGICAL, INTENT (IN) ::           DO_FOCORR_ALONE
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_SURFACE_LINEARIZATION

      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE

!  Numbers, basic

      INTEGER, INTENT (IN) ::           FOURIER, IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS

      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           N_OUT_STREAMS

!  Linearization control

      LOGICAL, INTENT (IN) ::           LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           N_SURFACE_WFS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

!  Numbers, derived

      INTEGER, INTENT (IN) ::           N_DIRECTIONS
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Azimuth factors

      DOUBLE PRECISION, INTENT (IN) ::  AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )
      INTEGER         , INTENT (IN) ::  SZD_OFFSETS ( MAX_SZANGLES )

!  Fourier components

      DOUBLE PRECISION, INTENT (IN) :: SURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
                                                     MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                     MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Type structure for Single scatter and Direct-beam weighting functions; We want Col/Surf substructures

      TYPE(VLIDORT_LinSup_SS), intent(in)     :: VLIDORT_LinSS

!  modified/output variables
!  -------------------------

!  Full Jacobian results, type definition; We want Col/Surf substructures

      TYPE(VLIDORT_LinOutputs), intent(inout)      :: VLIDORT_LinOut

!  local variables

      INTEGER :: LUA = 1
      INTEGER :: I, IDIR, UT, Q, N, W, Z, V, O1

!   For Single scatter corrections, the quadrature output flag is
!    turned off, so that N_OUT_STREAMS = N_USER_STREAMS, and the
!    offsetting is consistent here with the usage elsewhere.

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  single scatter calculation is initialized to zero here (Version 2.7)

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_OUT_STREAMS
                        V = SZD_OFFSETS(IBEAM) + I
                        DO O1 = 1, NSTOKES
                          VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) = PROFILEWF_F(Q,N,UT,I,IBEAM,O1,W)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ELSE
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_OUT_STREAMS
                        V = SZD_OFFSETS(IBEAM) + I
                        DO O1 = 1, NSTOKES
                          VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) = ZERO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!    Version 2.8 some renaming of variables
          !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!          IF ( DO_SSFULL .OR. &
!               ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

          IF ( DO_FOCORR ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_OUT_STREAMS
                        V = SZD_OFFSETS(IBEAM) + I
                        DO O1 = 1, NSTOKES
                          VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) = &
                            VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) + VLIDORT_LinSS%Prof%TS_PROFILEWF_SS(Q,N,UT,V,O1,W)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Add the direct bounce to the upwelling if flagged
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!       Version 2.8    Changed Logic for SS terms

          !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
!          IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN
          IF ( (DO_FOCORR.OR.DO_DBCORRECTION) .AND.DO_UPWELLING ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO UT = 1, N_USER_LEVELS
                    DO I = 1, N_OUT_STREAMS
                      V = SZD_OFFSETS(IBEAM) + I
                      DO O1 = 1, NSTOKES
                        !CALL TP23F (FOURIER,Q,N,UT,V,O1,PROFILEWF,PROFILEWF_DB)
                        VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,UPIDX) = &
                          VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,UPIDX) + VLIDORT_LinSS%Prof%TS_PROFILEWF_DB(Q,N,UT,V,O1)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  end atmospheric WF clause

        ENDIF

!  Reflectance weighting functions
!  -------------------------------

!  Diffuse field at all output angles
!  Alternative - zero the output (single scatter case)

        IF ( DO_SURFACE_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
!mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
            DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    V = SZD_OFFSETS(IBEAM) + I
                    DO O1 = 1, NSTOKES
                      VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,W) = SURFACEWF_F(Z,UT,I,IBEAM,O1,W)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
!mick fix 9/12/2012 - added N_SLEAVE_WFS for Z
            DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    V = SZD_OFFSETS(IBEAM) + I
                    DO O1 = 1, NSTOKES
                      VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,W) = ZERO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 2.3
!     Version 2.8    Changed Logic for SS terms

          !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
          !IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL)  .AND.DO_UPWELLING ) THEN

          IF ( (DO_FOCORR.OR.DO_DBCORRECTION) .AND.DO_UPWELLING ) THEN
!mick fix 9/12/2012 - added N_SLEAVE_WFS for Z
            DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_OUT_STREAMS
                  V = SZD_OFFSETS(IBEAM) + I
                  DO O1 = 1, NSTOKES
                    !CALL TP26F (FOURIER,Z,UT,V,O1,SURFACEWF,SURFACEWF_DB)
                    VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,UPIDX) = &
                      VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,UPIDX) + VLIDORT_LinSS%Surf%TS_SURFACEWF_DB(Z,UT,V,O1)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  End surface linearization WF clause

        ENDIF

!  If no_azimuth, then exit
!  ------------------------

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Add next Fourier component to output

          DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                DO IDIR = 1, N_DIRECTIONS
                  W = WHICH_DIRECTIONS(IDIR)
                  DO UT = 1, N_USER_LEVELS
                    DO I = 1, N_OUT_STREAMS
                      V = SZD_OFFSETS(IBEAM) + I
                      DO O1 = 1, NSTOKES
                        VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) = &
                        VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,V,O1,W) + PROFILEWF_F(Q,N,UT,I,IBEAM,O1,W)*AZMFAC(I,IBEAM,LUA,O1)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO

!  End profile atmospheric WF clause

        ENDIF

!  albedo weighting functions (non-Lambertian only)
!  ------------------------------------------------

        IF ( DO_SURFACE_LINEARIZATION ) THEN
          IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN

!  Copy full output

!mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
            DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_OUT_STREAMS
                    V = SZD_OFFSETS(IBEAM) + I
                    DO O1 = 1, NSTOKES
                      !CALL TP26G (FOURIER,Z,UT,I,IBEAM,V,O1,W,SURFACEWF,SURFACEWF_F)
                      VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,W) = &
                        VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,O1,W) + SURFACEWF_F(Z,UT,I,IBEAM,O1,W)*AZMFAC(I,IBEAM,LUA,O1)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

!  End surface linearization clause

          ENDIF
        ENDIF

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LPS_CONVERGE_DOUBLET

!

      SUBROUTINE VLIDORT_LPS_CONVERGE_OBSGEO ( &
        DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH, DO_DBCORRECTION,           & ! Input flags
        DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_LAMBERTIAN_SURFACE,          & ! Input flags
        NSTOKES, NLAYERS, N_USER_LEVELS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, & ! Input linearization control
        N_SLEAVE_WFS, IBEAM, FOURIER, N_DIRECTIONS, WHICH_DIRECTIONS,                       & ! Input bookkeeping
        AZMFAC, PROFILEWF_F, SURFACEWF_F, VLIDORT_LinSS, VLIDORT_LinOut )                     ! Input Azm, fields

!  Just upgrades the weighting function Fourier cosine series
!   Version 2.8, 3/1/17. Logic for FOCORR variables changed

!  4/15/20. Version 2.8.2
!   -- Replaced SURFACE_DB/PROFILEWF_SS/DB (FO inputs), PROFILEWF/SURFACEWF (final outputs) with Type structure variables
!   -- Rearranged argument list, removed DO_FOCORR_EXTERNAL

!  module, dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAX_GEOMETRIES, MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS, &
                                 MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS,  &
                                 MAX_USER_LEVELS, MAX_ATMOSWFS, MAX_SURFACEWFS, ZERO, UPIDX

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_FOCORR
      LOGICAL, INTENT (IN) ::           DO_FOCORR_ALONE
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH

!  More flags

      LOGICAL, INTENT (IN) ::           DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE

!  Control integers

      INTEGER, INTENT (IN) ::           NSTOKES, NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS

!  Linearization control

      LOGICAL, INTENT (IN) ::           LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           N_SURFACE_WFS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

!  Fourier component and beam index

      INTEGER, INTENT (IN) ::           FOURIER, IBEAM

!  Directional control

      INTEGER, INTENT (IN) ::           N_DIRECTIONS
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  azimuth factors

      DOUBLE PRECISION, INTENT (IN) ::  AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )

!  Fourier-component Profile/Surface weighting functions at user angles

      DOUBLE PRECISION, INTENT (IN) :: SURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
                                                     MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                     MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Type structure for Single scatter and Direct-beam weighting functions; We want Prof/Surf substructures

      TYPE(VLIDORT_LinSup_SS), intent(in)     :: VLIDORT_LinSS

!  modified/output variables
!  -------------------------

!  Full Jacobian results, type definition; We want Prof/Surf substructures

      TYPE(VLIDORT_LinOutputs), intent(inout)      :: VLIDORT_LinOut

!  local variables
!  ---------------

      INTEGER :: IDIR, UT, Q, N, W, Z, O1, IB, LUM, LUA

!  Local user indices

      IB  = IBEAM
      LUM = 1
      LUA = 1

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  single scatter calculation is initialized to zero here (Version 2.7)

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO O1 = 1, NSTOKES
                        VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,IB,O1,W) = PROFILEWF_F(Q,N,UT,LUM,IB,O1,W)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ELSE
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO O1 = 1, NSTOKES
                        VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,IB,O1,W) = ZERO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!    Version 2.8 some renaming of variables

          !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!          IF ( DO_SSFULL .OR. &
!             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

          IF ( DO_FOCORR ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO O1 = 1, NSTOKES
                        !CALL TP23E2 (FOURIER,Q,N,UT,IB,O1,W,PROFILEWF,PROFILEWF_SS)
                        VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,IB,O1,W) = &
                 VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,IB,O1,W) + VLIDORT_LinSS%Prof%TS_PROFILEWF_SS(Q,N,UT,IB,O1,W)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 2.3
!     Version 2.8    Changed Logic for SS terms

          !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
          !IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN

          IF ( (DO_FOCORR.OR.DO_DBCORRECTION) .AND.DO_UPWELLING ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO UT = 1, N_USER_LEVELS
                     DO O1 = 1, NSTOKES
                       !CALL TP23F2 (FOURIER,Q,N,UT,IB,O1,PROFILEWF,PROFILEWF_DB)
                       VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,IB,O1,UPIDX) = &
                 VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,IB,O1,UPIDX) + VLIDORT_LinSS%Prof%TS_PROFILEWF_DB(Q,N,UT,IB,O1)
                     ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  end atmospheric WF clause

        ENDIF

!  Reflectance weighting functions
!  -------------------------------

!  Diffuse field at all output angles
!  Alternative - zero the output (single scatter case)

        IF ( DO_SURFACE_LINEARIZATION ) THEN

          IF ( .not. DO_FOCORR_ALONE ) THEN
!mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
            DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO O1 = 1, NSTOKES
                    VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,O1,W) = SURFACEWF_F(Z,UT,LUM,IB,O1,W)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
!mick fix 9/12/2012 - added N_SLEAVE_WFS for Z
            DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO O1 = 1, NSTOKES
                    VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,O1,W) = ZERO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 2.3
!     Version 2.8    Changed Logic for SS terms

          !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
          !IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN

          IF ( (DO_FOCORR.OR.DO_DBCORRECTION) .AND.DO_UPWELLING ) THEN
!mick fix 9/12/2012 - added N_SLEAVE_WFS for Z
            DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO UT = 1, N_USER_LEVELS
                DO O1 = 1, NSTOKES
                  !CALL TP26F2 (FOURIER,Z,UT,IB,O1,SURFACEWF,SURFACEWF_DB)
                  VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,O1,UPIDX) = &
         VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,O1,UPIDX) + VLIDORT_LinSS%Surf%TS_SURFACEWF_DB(Z,UT,IB,O1)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  End surface linearization clause

        ENDIF

!  If no_azimuth, then exit
!  ------------------------

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Add next Fourier component to output

          DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                DO IDIR = 1, N_DIRECTIONS
                  W = WHICH_DIRECTIONS(IDIR)
                  DO UT = 1, N_USER_LEVELS
                    DO O1 = 1, NSTOKES
                      !CALL TP23G2 (FOURIER,Q,N,UT,LUM,IB,O1,W,PROFILEWF,PROFILEWF_F)
                      VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,IB,O1,W) = &
          VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,UT,IB,O1,W) + PROFILEWF_F(Q,N,UT,LUM,IB,O1,W)*AZMFAC(LUM,IB,LUA,O1)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO

!  End profile atmospheric WF clause

       ENDIF

!  albedo weighting functions (non-Lambertian only)
!  ------------------------------------------------

       IF ( DO_SURFACE_LINEARIZATION ) THEN
         IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN

!  Copy full output

!mick fix 9/6/2012 - added N_SLEAVE_WFS for Z

           DO Z = 1, N_SURFACE_WFS + N_SLEAVE_WFS
             DO IDIR = 1, N_DIRECTIONS
               W = WHICH_DIRECTIONS(IDIR)
               DO UT = 1, N_USER_LEVELS
                 DO O1 = 1, NSTOKES
                   !CALL TP26G2 (FOURIER,Z,UT,LUM,IB,O1,W,SURFACEWF,SURFACEWF_F)
                   VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,O1,W) = &
             VLIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,O1,W) + SURFACEWF_F(Z,UT,LUM,IB,O1,W)*AZMFAC(LUM,IB,LUA,O1)
                 ENDDO
               ENDDO
             ENDDO
           ENDDO

!  End surface linearization clause

         ENDIF
       ENDIF

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LPS_CONVERGE_OBSGEO

!  End Module

      END MODULE vlidort_lps_converge_m

