
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
! #  For linearizations with regular BVP                        #
! #                                                             #
! #            LC_BVP_SOLUTION_MASTER                           #
! #            LC_BVP_COLUMN_SETUP  (new for version 2.4)       #
! #            LC_BEAMSOLUTION_NEQK (new for version 2.4)       #
! #                                                             #
! #  For linearizations with telescoped BVP                     #
! #                                                             #
! #            LC_BVPTEL_SOLUTION_MASTER                        #
! #            LC_BVPTEL_COLUMN_SETUP                           #
! #            LC_BVPTEL_SURFACE_SETUP  (New, Version 2.8)      #
! #                                                             #
! ###############################################################

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.
!     * Replaced GOTO constructions with IF blocks
!     * Implementation of BRDFs in the Linearized Telescoped problem

!  4/15/20. Version 2.8.2. BRDF and SLEAVE arrays are defined locally, each Fourier.

      MODULE vlidort_lc_bvproblem_m

      private :: LC_BVP_COLUMN_SETUP,    LC_BEAMSOLUTION_NEQK, &
                 LC_BVPTEL_COLUMN_SETUP, LC_BVPTEL_SURFACE_SETUP
      public  :: LC_BVP_SOLUTION_MASTER, LC_BVPTEL_SOLUTION_MASTER

      CONTAINS

      SUBROUTINE LC_BVP_SOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTBEAM,       & ! Flags
        DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, FOURIER,     & ! Flags, FOURIER
        DO_LAYER_SCATTERING, IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_COLUMNWFS,    & ! Numbers (Basic)
        NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, & ! Numbers (derived)
        MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR, QUAD_STRMWTS,         & ! Bookkeeping
        DELTAU_SLANT, L_DELTAU_VERT, RF_DIRECT_BEAM, SLTERM,    & ! Optical and direct beam
        BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, BVEC,         & ! Beam and PI Solutions
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON,         & ! Homogeneous Solutions
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT, ALBEDO, BRDF_F,       & ! BVP matrices, surface
        L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG, L_T_WUPPER, L_T_WLOWER, & ! Linearized Homog/Thermal
        LC_TRANS_ATMOS, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_BVEC,       & ! Linearized Beam/PI solutions
        NCON, PCON, L_WLOWER, L_WUPPER,                         & ! output solutions
        STATUS, MESSAGE, TRACE )                                  ! Exception handling

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ,   &
                                 MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS

      USE VLIDORT_LPC_BVPROBLEM_m, Only : L_BVP_BACKSUB

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM

      LOGICAL, INTENT (IN) ::          DO_WATER_LEAVING     
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS

!  Basic control numbers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_COLUMNWFS

!  other numbers


      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG

!  Layer scattering

      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  4/9/19. RF_DIRECT_BEAM is the reflected beam (excludes SLTERM)

      DOUBLE PRECISION, INTENT (IN) :: RF_DIRECT_BEAM   ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SLTERM           ( MAXSTREAMS, MAXSTOKES )

!  optical

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  homogeneous RTE solutions and inegration constants

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Beam and particular integral vector

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  BVP matrices

      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOT   ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2   ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTRMSTKS_2 )

!  Albedo, Fourier components of BRDF incident quadrature streams ==> reflected quadrature streams
!  4/15/20. Version 2.8.2. BRDF array defined locally, each Fourier, remove MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

!  Linearized homogeneous solutions

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized thermal particular integerals

      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)

!  Linearized PI vector, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!  4/9/19 Add linearization of TRANS_ATMOS for the waterleaving contribution

      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  :: LC_TRANS_ATMOS   ( MAXBEAMS, MAX_ATMOSWFS )

!  output linearized integration constants

      DOUBLE PRECISION, INTENT (OUT) ::  NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  output particular integrals

      DOUBLE PRECISION, INTENT (OUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  local help variables

      INTEGER :: STATUS_SUB, VAR_INDEX

!  Arrays for linearized column vectors

      DOUBLE PRECISION :: COL2_WF  ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  module status and message initialization

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Variation index  = 0 for the Column weighting function case

      VAR_INDEX = 0

!  Linearization of the regular BVP case
!  =====================================

!  Set up the column vectors for Bulk linearizations
!  -------------------------------------------------

!  Bulk: Compute the main column B' where AX = B'
! 4/9/19. Add water-leaving control, reflected direct beam, surface leaving linearization contribution
!  4/15/20. Version 2.8.2. BRDF arrays defined locally, each Fourier, remove MAXMOMENTS dimension

      CALL LC_BVP_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTBEAM,     & ! Input, Flags
        DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, FOURIER,   & ! Input, Flags/FOURIER
        DO_LAYER_SCATTERING, IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_COLUMNWFS,  & ! Input, Basic numbers
        NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2,                     & ! Input, Other numbers
        DELTAU_SLANT, L_DELTAU_VERT, MUELLER_INDEX, K_REAL, K_COMPLEX,        & ! Input, Optical/Bookkeeping
        SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SLTERM, & ! Input, Surface Stuff
        BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, BVEC,   & ! Beam/solarPI solutions
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON,   & ! Homogeneous  solutions
        L_SOLA_XPOS, L_SOLB_XNEG, L_T_DELT_EIGEN, L_T_WUPPER, L_T_WLOWER, & ! Linearized Homog/thermal solutions
        LC_TRANS_ATMOS, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_BVEC,       & ! Linearized Beam/SolarPI  solutions
        L_WLOWER, L_WUPPER, COL2_WF, SCOL2_WF )             ! PI Solutions and Column vectors

!  Back-substitution

      CALL L_BVP_BACKSUB ( &
        VAR_INDEX, N_COLUMNWFS, NLAYERS, NTOTAL,     & ! Numbers
        N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS,          & ! Numbers
        NSTKS_NSTRMS_2, K_REAL, K_COMPLEX,           & ! Numbers
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT,            & ! Input BVP matrices
        COL2_WF, SCOL2_WF,                           & ! Input BVP vectors
        NCON, PCON, STATUS_SUB, MESSAGE, TRACE )       ! Output and excpetion handling

!  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS ; RETURN
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LC_BVP_SOLUTION_MASTER

!

      SUBROUTINE LC_BVP_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTBEAM,     & ! Input, Flags
        DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, FOURIER,   & ! Input, Flags/FOURIER
        DO_LAYER_SCATTERING, IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_COLUMNWFS,  & ! Input, Basic numbers
        NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2,                     & ! Input, Other numbers
        DELTAU_SLANT, L_DELTAU_VERT, MUELLER_INDEX, K_REAL, K_COMPLEX,        & ! Input, Optical/Bookkeeping
        SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SLTERM, & ! Input, Surface Stuff
        BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, BVEC,   & ! Beam/solarPI solutions
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON,   & ! Homogeneous  solutions
        L_SOLA_XPOS, L_SOLB_XNEG, L_T_DELT_EIGEN, L_T_WUPPER, L_T_WLOWER, & ! Linearized Homog/thermal solutions
        LC_TRANS_ATMOS, LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_BVEC,       & ! Linearized Beam/SolarPI  solutions
        L_WLOWER, L_WUPPER, COL2_WF, SCOL2_WF )             ! PI Solutions and Column vectors

!  Linearized column vector setup (bulk property linearization)
! 4/9/19. Add water-leaving control, reflected direct beam, surface leaving linearization contribution
  
!  4/15/20. Version 2.8.2. BRDF array defined locally, each Fourier.

!  module, dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ,   &
                                 MAXBANDTOTAL, MAXTOTAL, ZERO

      USE VLIDORT_LPC_BVPROBLEM_m, Only : L_BVP_SURFACE_SETUP

      IMPLICIT NONE

!  Subroutine arguments
!  ====================

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_WATER_LEAVING
   
      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Basic control numbers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_COLUMNWFS

!  other numbers


      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  optical 

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  Beam and particular integral vector

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  homogeneous RTE solutions and integration constants

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Surface inputs
!    Fourier components of BRDF incident quadrature streams, reflected quadrature streams
!    4/15/20. Version 2.8.2. BRDF arrays defined locally, each Fourier, remove MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  4/9/19. RF_DIRECT_BEAM is the reflected beam (excludes SLTERM)

      DOUBLE PRECISION, INTENT (IN) :: RF_DIRECT_BEAM   ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SLTERM           ( MAXSTREAMS, MAXSTOKES )

!  Linearized
!  ----------

!  Linearized PI Vector, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!  4/9/19 Add linearization of TRANS_ATMOS for the waterleaving contribution

      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_TRANS_ATMOS    ( MAXBEAMS,  MAX_ATMOSWFS )

!  Linearized homogeneous solutions

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized thermal particular integerals

      DOUBLE PRECISION, INTENT (IN) :: L_T_WUPPER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_WLOWER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)

!  Output arguments
!  ----------------

!  output particular integrals

      DOUBLE PRECISION, INTENT (OUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  output Column vectors

      DOUBLE PRECISION, INTENT (OUT) :: COL2_WF ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER ::          Q, N, N1, I, I1, IR, IROW, CM, C0, O1
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: CPOS, CNEG, L_HOM_R, L_HOM_CR, L_BEAM, FAC
      DOUBLE PRECISION :: T1, T2, T1R, T1I, T2R, T2I, LCTERM
      DOUBLE PRECISION :: L_HOM_U_R  , L_HOM_D_R
      DOUBLE PRECISION :: L_HOM_U_CR , L_HOM_D_CR
      LOGICAL ::          MODIFIED_BOUNDARY

!  reflected arrays

      DOUBLE PRECISION :: R2_L_BEAM ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

!  initialize
!  ----------

!  zero the results vectors

!  Enhancement # 1, 6/27/16
      COL2_WF(1:NTOTAL,1:MAX_ATMOSWFS) = ZERO

!  Copy already existing thermal linearizations
!    This is a very important zeroing.................!!!!!

!  Enhancement # 2, 6/27/16
      IF ( DO_INCLUDE_THERMEMISS ) THEN
        L_WUPPER(1:NSTREAMS_2, 1, 1:NLAYERS, 1:N_COLUMNWFS) = L_T_WUPPER(1:NSTREAMS_2, 1:NLAYERS, 1:N_COLUMNWFS)
        L_WLOWER(1:NSTREAMS_2, 1, 1:NLAYERS, 1:N_COLUMNWFS) = L_T_WLOWER(1:NSTREAMS_2, 1:NLAYERS, 1:N_COLUMNWFS)
        L_WUPPER(1:NSTREAMS_2, 2:NSTOKES, 1:NLAYERS, 1:N_COLUMNWFS) = ZERO
        L_WLOWER(1:NSTREAMS_2, 2:NSTOKES, 1:NLAYERS, 1:N_COLUMNWFS) = ZERO
      ELSE
        L_WUPPER(1:NSTREAMS_2, 1:NSTOKES, 1:NLAYERS, 1:N_COLUMNWFS) = ZERO
        L_WLOWER(1:NSTREAMS_2, 1:NSTOKES, 1:NLAYERS, 1:N_COLUMNWFS) = ZERO
      ENDIF

!  Top of first layer (TOA), UPPER boundary condition
!  --------------------------------------------------

      N = 1

!  Get the linearized beam solution for the first layer

      IF ( DO_SOLAR_SOURCES ) THEN
        CALL LC_BEAMSOLUTION_NEQK ( &
          FOURIER, IBEAM, N, N_COLUMNWFS, NSTOKES,      & ! basic control
          NSTREAMS_2, DO_LAYER_SCATTERING, BEAM_CUTOFF, & ! Bookkeeping
          T_DELT_MUBAR, INITIAL_TRANS, BVEC,            & ! Beam/SolarPI nputs
          LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_BVEC,   & ! Linearized Beam/SolarPI 
          L_WUPPER, L_WLOWER )                            ! Output solutions
      ENDIF

!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          DO Q = 1, N_COLUMNWFS

!  beam solution linearization at top of layer

            L_BEAM = - L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 3, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!            L_HOM_R = sum ( ( LCON(1:K_REAL(N),N) * L_SOLA_XPOS(I,O1,1:K_REAL(N),N,Q) ) + &
!               ( MCON(1:K_REAL(N),N) * ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLB_XNEG(I,O1,1:K_REAL(N),N,Q) + &
!                 L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLB_XNEG(I,O1,1:K_REAL(N),N) ) ) )
            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
              CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                   L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
              T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
              T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contribution

            COL2_WF(IROW,Q) = L_BEAM - L_HOM_R - L_HOM_CR
          ENDDO
        ENDDO
      ENDDO

!  Intermediate boundary conditions
!  --------------------------------

      DO N = 1, NLAYERS - 1

!  N1 is the layer below, C0 is the offset

        N1 = N + 1
        C0 = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  Get the linearized beam solution for the next layer

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL LC_BEAMSOLUTION_NEQK ( &
            FOURIER, IBEAM, N1, N_COLUMNWFS, NSTOKES,     & ! basic control
            NSTREAMS_2, DO_LAYER_SCATTERING, BEAM_CUTOFF, & ! Bookkeeping
            T_DELT_MUBAR, INITIAL_TRANS, BVEC,            & ! Beam/SolarPI nputs
            LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_BVEC,   & ! Linearized Beam/SolarPI 
            L_WUPPER, L_WLOWER )                            ! Output solutions
        ENDIF

!  .. 2 contributions to L_BEAM, from variations L_WUPPER L_WLOWER
!  .. 2 contributions to L_HOM,  from variations above and below

        DO I = 1, NSTREAMS_2
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_COLUMNWFS

!  Beam contributions

            L_BEAM  = + L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions above

!  Enhancement # 4, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!            L_HOM_R =sum( ( LCON(1:K_REAL(N),N) *  L_SOLA_XPOS(I,O1,1:K_REAL(N),N,Q) ) + &
!                ( MCON(1:K_REAL(N),N) * ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLB_XNEG(I,O1,1:K_REAL(N),N,Q) + &
!                 L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLB_XNEG(I,O1,1:K_REAL(N),N) ) ) )
            L_HOM_U_R = ZERO
            DO K = 1, K_REAL(N1)
              CPOS = L_SOLA_XPOS(I,O1,K,N1,Q)
              CNEG = T_DELT_EIGEN(K,N1)   * L_SOLB_XNEG(I,O1,K,N1,Q) + &
                   L_T_DELT_EIGEN(K,N1,Q) *   SOLB_XNEG(I,O1,K,N1)
              T1 = LCON(K,N1) * CPOS
              T2 = MCON(K,N1) * CNEG
              L_HOM_U_R = L_HOM_U_R + T1 + T2
            ENDDO

!  Linearized Real homogeneous solution contributions below

            L_HOM_D_R = ZERO
            DO K = 1, K_REAL(N)
              CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
              CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + &
                   L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_D_R = L_HOM_D_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions above

            L_HOM_U_CR  = ZERO
            KO1 = K_REAL(N1) + 1
            DO K = 1, K_COMPLEX(N1)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1 = L_SOLA_XPOS(I,O1,K1,N1,Q) * LCON(K1,N1) - &
                   L_SOLA_XPOS(I,O1,K2,N1,Q) * LCON(K2,N1)
              T2R =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                   - T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                   + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K1,N1) &
                   - L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
              T2I =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                   + T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                   + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K2,N1) &
                   + L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
              T2 =  T2R * MCON(K1,N1) - T2I * MCON(K2,N1)
              L_HOM_U_CR = L_HOM_U_CR + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions below

            L_HOM_D_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) - &
                   L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
              T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                   - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                   - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
              T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                   + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                   + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                   + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              L_HOM_D_CR = L_HOM_D_CR + T1 + T2
            ENDDO

!  Final contribution

            L_HOM_R  =  L_HOM_U_R  - L_HOM_D_R
            L_HOM_CR =  L_HOM_U_CR - L_HOM_D_CR
            COL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

          ENDDO
         ENDDO
        ENDDO

!  End layer

      ENDDO

!  FINAL layer
!  -----------

      N = NLAYERS
      MODIFIED_BOUNDARY = .TRUE.

!  get the linearized downward-reflected term
!  4/15/20. Version 2.8.2. BRDF array defined locally, remove Fourier index.

      CALL L_BVP_SURFACE_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, MODIFIED_BOUNDARY, & ! Flags
        IBEAM, FOURIER, N_COLUMNWFS, NSTOKES, NSTREAMS, NLAYERS,      & ! Numbers
        SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F,                 & ! Surface input
        NSTKS_NSTRMS, MUELLER_INDEX, K_REAL, K_COMPLEX,               & ! bookkeeping
        SOLA_XPOS, T_DELT_EIGEN, L_T_DELT_EIGEN,                      & ! RT Solutions
        L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER,                           & ! Linearized RT Solutions
        R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )                               ! Output reflected solutions

!  offsets

      C0  = (N-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

!  start loops

      DO I = 1, NSTREAMS
        IR = NSTOKES*(I-1)
        I1 = I + NSTREAMS
        DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM = C0 + IROW
          DO Q = 1, N_COLUMNWFS

!  Beam contributions

            L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 5, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!            L_HOM_R = sum( ( LCON(1:K_REAL(N),N) * ( ( T_DELT_EIGEN(1:K_REAL(N),N)   * L_SOLA_XPOS(I1,O1,1:K_REAL(N),N,Q) &
!                + L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLA_XPOS(I1,O1,1:K_REAL(N),N) ) - R2_L_HOMP(I,O1,1:K_REAL(N),Q) ) ) + &
!                ( MCON(1:K_REAL(N),N) * ( L_SOLB_XNEG(I1,O1,1:K_REAL(N),N,Q) - R2_L_HOMM(I,O1,1:K_REAL(N),Q) ) ) )
            L_HOM_R  = ZERO
            DO K = 1, K_REAL(N)
              CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                  + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
              CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
              CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
              CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
              T1 = LCON(K,N) * CPOS
              T2 = MCON(K,N) * CNEG
              L_HOM_R = L_HOM_R + T1 + T2
            ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

            L_HOM_CR  = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                  - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
              T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                  + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                  + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                  + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
              T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
              T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
              T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
              T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
              T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
              T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
              L_HOM_CR = L_HOM_CR + T1 + T2
            ENDDO

!  Final contributions

            COL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

          ENDDO
        ENDDO
      ENDDO

!  Add direct beam variation to Final boundary
!  -------------------------------------------

!  Must use the reflected beam here
      
      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
         DO Q = 1, N_COLUMNWFS
            DO I = 1, NSTREAMS
               IR = NSTOKES*(I-1) ; I1 = I + NSTREAMS
               DO O1 = 1, NSTOKES
                  IROW = IR + O1 ;  CM = C0 + IROW
                  FAC = - RF_DIRECT_BEAM(I,IBEAM,O1)
                  L_BEAM = DOT_PRODUCT(L_DELTAU_VERT(Q,1:NLAYERS),DELTAU_SLANT(N,1:NLAYERS,IBEAM))
!                  COL2_WF(CM,Q) = COL2_WF(CM,Q) - FAC * L_BEAM   Bug 8/12/19
                  COL2_WF(CM,Q) = COL2_WF(CM,Q) + FAC * L_BEAM
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  4/9/19 Add the linearization due to Adjusted waterleaving term

      IF ( DO_WATER_LEAVING .and. FOURIER .eq. 0 ) THEN
         DO Q = 1, N_COLUMNWFS
            LCTERM = LC_TRANS_ATMOS(IBEAM,Q) ; O1 = 1
            DO I = 1, NSTREAMS
               IR = NSTOKES*(I-1) ; I1 = I + NSTREAMS
               IROW = IR + O1 ;  CM = C0 + IROW
               COL2_WF(CM,Q) = COL2_WF(CM,Q) + LCTERM * SLTERM(I,O1)
            ENDDO
         ENDDO
      ENDIF
        
!  copy the single layer vector

      IF ( NLAYERS .EQ. 1 ) THEN
!  Enhancement # 6, 6/27/16
        SCOL2_WF(1:NTOTAL,N_COLUMNWFS) = COL2_WF(1:NTOTAL,N_COLUMNWFS)
      ENDIF

!  debug

!      if ( do_debug_write ) then
!        DO N = 1, NTOTAL
!          write(95,'(3i4,1p4e17.9)')
!     &      FOURIER,IBEAM,N,
!     &                 COL2_WF(N,1),COL2_WF(N,2)
!        ENDDO
!      ENDIF
!       pause'f95'

!    if ( fourier_component.eq.1 ) then
!      do n = 1, ntotal
!         write(*,*)n,Col2_wf(n,1),col2_wf(n,2)
!      enddo
!    endif


!  finish

      RETURN
      END SUBROUTINE LC_BVP_COLUMN_SETUP

!

      SUBROUTINE LC_BEAMSOLUTION_NEQK ( &
        FOURIER, IB, N, N_COLUMNWFS, NSTOKES,         & ! basic control
        NSTREAMS_2, DO_LAYER_SCATTERING, BEAM_CUTOFF, & ! Bookkeeping
        T_DELT_MUBAR, INITIAL_TRANS, BVEC,            & ! Beam/SolarPI nputs
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_BVEC,   & ! Linearized Beam/SolarPI 
        L_WUPPER, L_WLOWER )                            ! Output solutions

!  Linearization of beam particular integral in layer N, due to column.
!   This is the bulk property linearization

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, MAXSTREAMS_2
                               

      IMPLICIT NONE

!  Basic control

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IB
      INTEGER, INTENT (IN) ::          N
      INTEGER, INTENT (IN) ::          N_COLUMNWFS
      INTEGER, INTENT (IN) ::          NSTOKES

!  bookkeeping

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Beam/SolarPI input

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  Linearized Beam/SolarPI input

      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  output PI 

      DOUBLE PRECISION, INTENT (INOUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER ::          I, O1, Q
      DOUBLE PRECISION :: CONST, WDEL, VAR1, VAR2, VAR_U, TRANS2
      DOUBLE PRECISION :: LBSOL(MAXSTREAMS_2,MAXSTOKES,MAX_ATMOSWFS,2)

!  No linearized particular solution beyond the cutoff layer. ALSO--
!  Nothing if the solution saving mode is on and layer is inactive
!   ALSO - Nothing if layer has no scattering (Do not need solution savi

      IF (.NOT.DO_LAYER_SCATTERING(FOURIER,N) &
             .OR. (N .GT.BEAM_CUTOFF(IB))) THEN
        RETURN
      ENDIF

!  Classical solution
!  ==================

!  Very simple, same code for all situations

      CONST   = INITIAL_TRANS(N,IB)
      WDEL    = T_DELT_MUBAR(N,IB)
      TRANS2  = CONST * WDEL
      DO Q = 1, N_COLUMNWFS
        VAR1 = LC_T_DELT_MUBAR (N,IB,Q) * CONST
        VAR2 = LC_INITIAL_TRANS(N,IB,Q)
        DO I = 1, NSTREAMS_2
          DO O1 = 1, NSTOKES
          VAR_U = VAR2 * BVEC(I,O1,N) + LC_BVEC(I,O1,N,Q)
            LBSOL(I,O1,Q,1) = CONST  * VAR_U
            LBSOL(I,O1,Q,2) = TRANS2 * VAR_U + VAR1 * BVEC(I,O1,N)
          ENDDO
        ENDDO
      ENDDO

! Add to existing solution

!  Enhancement # 7, 6/27/16
      L_WUPPER(1:NSTREAMS_2,1:NSTOKES,N,1:N_COLUMNWFS) = &
        L_WUPPER(1:NSTREAMS_2,1:NSTOKES,N,1:N_COLUMNWFS) + LBSOL(1:NSTREAMS_2,1:NSTOKES,1:N_COLUMNWFS,1)
      L_WLOWER(1:NSTREAMS_2,1:NSTOKES,N,1:N_COLUMNWFS) = &
        L_WLOWER(1:NSTREAMS_2,1:NSTOKES,N,1:N_COLUMNWFS) + LBSOL(1:NSTREAMS_2,1:NSTOKES,1:N_COLUMNWFS,2)

!  Finish

      RETURN
      END SUBROUTINE LC_BEAMSOLUTION_NEQK

!

      SUBROUTINE LC_BVPTEL_SOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTBEAM, & ! Flags
        FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_COLUMNWFS,          & ! Basic Control Numbers
        NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, DO_LAYER_SCATTERING,    & ! Numbers and scattering
        N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,              & ! BVPTel Control
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,                & ! BVPTel Control
        MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR,            & ! Bookkeeping
        DELTAU_SLANT, L_DELTAU_VERT, DIRECT_BEAM,               & ! Optical and direct beam
        T_DELT_DISORDS, L_T_DELT_DISORDS,                       & ! Discrete ordinate transmittances
        BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, BVEC, WLOWER, & ! Beam and PI Solutions
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON,         & ! Homogeneous Solutions
        QUAD_STRMWTS, LAMBERTIAN_ALBEDO, BRDF_F,                & ! Surface inputs
        L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG,               & ! Linearized Homog Solutions
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_BVEC,             & ! Linearized Beam/PI solutions
        BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT,                 & ! BVPTel matrices
        L_WLOWER, L_WUPPER, NCON, PCON,                         & ! output solutions
        STATUS, MESSAGE, TRACE )                                  ! Exception handling

!  Linearization of the Telescoped Boundary Problem Solution
!   Version 2.8. Major extension to BRDFs with TELESCOPING. Implemented May 2016.

!  4/15/20. Version 2.8.2. BRDF and SLEAVE arrays are defined locally, each Fourier.

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ,   &
                                 MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO

      USE LAPACK_TOOLS_m, Only : DGBTRS, DGETRS

      USE VLIDORT_LPC_BVPROBLEM_m, Only : L_BVP_BACKSUB

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM

!  Basic control numbers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_COLUMNWFS

!  other numbers


      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  Layer scattering

      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Telescoped BVP Control

      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (IN) ::          NLAYERS_TEL
      INTEGER, INTENT (IN) ::          ACTIVE_LAYERS ( MAXLAYERS )

!  optical and direct beam

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM   ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

!  discrete ordinate stuff

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR

!  homogeneous RTE solutions and inegration constants

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Beam and particular integral vector

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  Surface inputs
!    Fourier components of BRDF incident quadrature streams, reflected quadrature streams
!  4/15/20. Version 2.8.2. BRDF array defined locally, each Fourier, remove MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Linearized homogeneous solutions

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Beam and PI vector

      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  BVP matrices

      DOUBLE PRECISION, INTENT (IN) :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOTTEL   ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2   ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTRMSTKS_2 )

!  output linearized integration constants

      DOUBLE PRECISION, INTENT (OUT) ::  NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) ::  PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!  output particular integrals

      DOUBLE PRECISION, INTENT (OUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  error tracing variables

      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI

!  Other local help variables

      INTEGER ::          I, I1, N, N1, NAL, NAL1, O1, IC, ICOW, Q
      INTEGER ::          K, KO1, K0, K1, K2, C0, NS
      INTEGER ::          IR, IROW, IROW1, IROW_S, IROW1_S
      DOUBLE PRECISION :: SPAR, SHOM, L_HOM1, L_HOM2, SHOM_R
      DOUBLE PRECISION :: SHOM_CR, L_HOM1CR, L_HOM2CR
      DOUBLE PRECISION :: LXR, MXR, NXR, PXR, LLXR, MLXR
      DOUBLE PRECISION :: LXR1, MXR1, NXR1, PXR1, LLXR1, MLXR1
      DOUBLE PRECISION :: LXR2, MXR2, NXR2, PXR2, LLXR2, MLXR2

!  Arrays for linearized column vectors

      DOUBLE PRECISION :: COLTEL2_WF ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION :: SCOL2_WF   ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  Initialize

      STATUS = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  What is this ?

      Q = 0

!  Bulk: Compute the main column B' where AX = B'
!  4/15/20. Version 2.8.2. BRDF array defined locally, each Fourier, remove MAXMOMENTS dimension

      CALL LC_BVPTEL_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTBEAM, & ! Flags
        FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_COLUMNWFS,          & ! Basic Control Numbers
        NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, DO_LAYER_SCATTERING,    & ! Numbers and scattering
        N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,         & ! BVPTel Control
        MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR,       & ! Bookkeeping
        DELTAU_SLANT, L_DELTAU_VERT, DIRECT_BEAM,               & ! Optical and direct beam
        T_DELT_DISORDS, L_T_DELT_DISORDS,                       & ! Discrete ordinate transmittances
        BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, BVEC, WLOWER, & ! Beam and PI Solutions
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON,         & ! Homogeneous Solutions
        QUAD_STRMWTS, LAMBERTIAN_ALBEDO, BRDF_F,                & ! Surface inputs
        L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG,               & ! Linearized Homog Solutions
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_BVEC,             & ! Linearized Beam/PI solutions
        L_WLOWER, L_WUPPER, COLTEL2_WF, SCOL2_WF )                ! output solutions and Column vector

!  Solve linearized BVP: Several Active layers
!  ===========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  BV solution for linearized integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUBDIAG, &
              N_BVTELMATRIX_SUPDIAG, N_COLUMNWFS, &
              BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, &
              COLTEL2_WF, MAXTOTAL, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (NLAYERS>1)in LC_BVPTEL_SOLUTION_MASTER'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set linearized integration constants, active layers

        C0 = - NSTKS_NSTRMS_2
        DO NS = 1, NLAYERS_TEL
         N = ACTIVE_LAYERS(NS)
         C0 = C0 + NSTKS_NSTRMS_2

!  set real constants from the solution vector

!  Enhancement # 8, 6/27/16
         NCON(1:K_REAL(N),N,1:N_COLUMNWFS) = COLTEL2_WF(C0+1:C0+K_REAL(N), 1:N_COLUMNWFS)
         PCON(1:K_REAL(N),N,1:N_COLUMNWFS) = COLTEL2_WF(C0+NSTKS_NSTRMS+1:C0+NSTKS_NSTRMS+K_REAL(N), 1:N_COLUMNWFS)

!  set complex constants from the solution vector

         KO1 = K_REAL(N) + 1
         DO K = 1, K_COMPLEX(N)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(N)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = K + K_REAL(N) + K_COMPLEX(N)
          IROW1_S = IROW_S + NSTKS_NSTRMS
!  Enhancement # 9, 6/27/16
          NCON(K1,N,1:N_COLUMNWFS) = COLTEL2_WF(C0+IROW,    1:N_COLUMNWFS)
          NCON(K2,N,1:N_COLUMNWFS) = COLTEL2_WF(C0+IROW_S,  1:N_COLUMNWFS)
          PCON(K1,N,1:N_COLUMNWFS) = COLTEL2_WF(C0+IROW1,   1:N_COLUMNWFS)
          PCON(K2,N,1:N_COLUMNWFS) = COLTEL2_WF(C0+IROW1_S, 1:N_COLUMNWFS)
         ENDDO

!  End number of telescoped layers

        ENDDO

!  Solve linearized BVP: Single Layer only
!  =======================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

        NAL = ACTIVE_LAYERS(1)

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS &
           ( 'N', NSTKS_NSTRMS_2, N_COLUMNWFS, &
              SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
              SCOL2_WF, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call (NLAYERS=1)in LC_BVPTEL_SOLUTION_MASTER'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  set real constants from the solution vector

!  Enhancement # 10, 6/27/16
        NCON(1:K_REAL(NAL),NAL,1:N_COLUMNWFS) = SCOL2_WF(1:K_REAL(NAL),1:N_COLUMNWFS)
        PCON(1:K_REAL(NAL),NAL,1:N_COLUMNWFS) = SCOL2_WF(NSTKS_NSTRMS+1:NSTKS_NSTRMS+K_REAL(NAL),1:N_COLUMNWFS)

!  set complex constants from the solution vector

        KO1 = K_REAL(NAL) + 1
        DO K = 1, K_COMPLEX(NAL)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(NAL)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = IROW + K_COMPLEX(NAL)
          IROW1_S = IROW_S + NSTKS_NSTRMS
!  Enhancement # 11, 6/27/16
          NCON(K1,NAL,1:N_COLUMNWFS) = SCOL2_WF(IROW,    1:N_COLUMNWFS)
          NCON(K2,NAL,1:N_COLUMNWFS) = SCOL2_WF(IROW_S,  1:N_COLUMNWFS)
          PCON(K1,NAL,1:N_COLUMNWFS) = SCOL2_WF(IROW1,   1:N_COLUMNWFS)
          PCON(K2,NAL,1:N_COLUMNWFS) = SCOL2_WF(IROW1_S, 1:N_COLUMNWFS)
        ENDDO

      ENDIF

!  Set linearized integration constants for non-active layers
!  ==========================================================

!  Now we propagate the results upwards and downwards through the
!  appropriate non-active layers where there is no scattering.

!  Transmittance layers ABOVE active layer(s)
!  -----------------------------------------

!   --NCON values are zero (no downwelling radiation)
!   --PCON values propagated upwards from top of first active layer

!  layer immediately above first active layer
!   --- Require linearized solutions at top of first active layer
!   --- Additional linearizations required, first layer is always active

      NAL = ACTIVE_LAYERS(1)
      IF ( NAL .GT. 1 ) THEN
        N1 = NAL - 1

!  start stream, stokes and parameter loops

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = ( I1 - 1 ) * NSTOKES
          IC = ( I - 1  ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            ICOW = IC + O1
            DO Q = 1, N_COLUMNWFS

!  real homogeneous solutions
!    Needs checking------------------------ 22 March 2007

!  Enhancement # 11, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!             SHOM_R = sum( ( NCON(1:K_REAL(NAL),NAL,Q) * SOLA_XPOS(I1,O1,1:K_REAL(NAL),NAL) + & 
!                       LCON(1:K_REAL(NAL),NAL)   * L_SOLA_XPOS(I1,O1,1:K_REAL(NAL),NAL,Q) ) + &
!                     ( T_DELT_EIGEN(1:K_REAL(NAL),NAL) * & 
!                       ( PCON(1:K_REAL(NAL),NAL,Q) * SOLB_XNEG(I1,O1,1:K_REAL(NAL),NAL) + & 
!                       MCON(1:K_REAL(NAL),NAL) * L_SOLB_XNEG(I1,O1,1:K_REAL(NAL),NAL,Q) ) + &
!                       L_T_DELT_EIGEN(1:K_REAL(NAL),NAL,Q) * MCON(1:K_REAL(NAL),NAL) * & 
!                       SOLB_XNEG(I1,O1,1:K_REAL(NAL),NAL) ) )
              SHOM_R = ZERO
              DO K = 1, K_REAL(NAL)
                NXR  = NCON(K,NAL,Q) *   SOLA_XPOS(I1,O1,K,NAL)
                PXR  = PCON(K,NAL,Q) *   SOLB_XNEG(I1,O1,K,NAL)
                MXR  = MCON(K,NAL)   *   SOLB_XNEG(I1,O1,K,NAL)
                LLXR = LCON(K,NAL)   * L_SOLA_XPOS(I1,O1,K,NAL,Q)
                MLXR = MCON(K,NAL)   * L_SOLB_XNEG(I1,O1,K,NAL,Q)
                L_HOM1 = NXR + LLXR
                L_HOM2 =   T_DELT_EIGEN(K,NAL) * ( PXR + MLXR ) + L_T_DELT_EIGEN(K,NAL,Q) * MXR
                SHOM_R = SHOM_R + L_HOM1 + L_HOM2
              ENDDO

!  complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(NAL) + 1
              DO K = 1, K_COMPLEX(NAL)
                K0 = 2*K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  = NCON(K1,NAL,Q) *   SOLA_XPOS(I1,O1,K1,NAL) - NCON(K2,NAL,Q) *   SOLA_XPOS(I1,O1,K2,NAL)
                PXR1  = PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL) - PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL)
                PXR2  = PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL) + PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL)
                MXR1  = MCON(K1,NAL) *   SOLB_XNEG(I1,O1,K1,NAL) - MCON(K2,NAL) *   SOLB_XNEG(I1,O1,K2,NAL)
                MXR2  = MCON(K1,NAL) *   SOLB_XNEG(I1,O1,K2,NAL) + MCON(K2,NAL) *   SOLB_XNEG(I1,O1,K1,NAL)
                LLXR1 = LCON(K1,NAL) * L_SOLA_XPOS(I1,O1,K1,NAL,Q) - LCON(K2,NAL) * L_SOLA_XPOS(I1,O1,K2,NAL,Q)
                MLXR1 = MCON(K1,NAL) * L_SOLB_XNEG(I1,O1,K1,NAL,Q) - MCON(K2,NAL) * L_SOLB_XNEG(I1,O1,K2,NAL,Q)
                MLXR2 = MCON(K1,NAL) * L_SOLB_XNEG(I1,O1,K2,NAL,Q) + MCON(K2,NAL) * L_SOLB_XNEG(I1,O1,K1,NAL,Q)
                L_HOM1CR = NXR1 + LLXR1
                L_HOM2CR =   T_DELT_EIGEN(K1,NAL)     * ( PXR1 + MLXR1 ) -   T_DELT_EIGEN(K2,NAL)   * ( PXR2 + MLXR2 ) &
                           + L_T_DELT_EIGEN(K1,NAL,Q) *   MXR1           - L_T_DELT_EIGEN(K2,NAL,Q) * MXR2
                SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
              ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

              SHOM = SHOM_R + SHOM_CR
              SPAR = L_WUPPER(I1,O1,NAL,Q)
              PCON(ICOW,N1,Q) = SPAR + SHOM
              NCON(ICOW,N1,Q) = ZERO

!  End loops

            ENDDO
          ENDDO
        ENDDO

!  End active layer varying

      ENDIF

!  For remaining non-active atmospheric layers to TOA, propagate upwards
!   Additional linearizations if you are passing through the varying lay

      DO N = NAL - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
!  Enhancement # 13, 6/27/16
            NCON(IROW,N,1:N_COLUMNWFS) = ZERO
            PCON(IROW,N,1:N_COLUMNWFS) = T_DELT_DISORDS(I,N1) * PCON(IROW,N,1:N_COLUMNWFS) &
                                     + L_T_DELT_DISORDS(I,N1,1:N_COLUMNWFS) * MCON(IROW,N1)
          ENDDO
        ENDDO
      ENDDO

!  Transmittance layers below active layer(s)
!  -----------------------------------------

!  This section substantially revised for Version 2.8
!   - General surface treatment for the Telescpied BVP

!       ** Only do this if active scattering is above (not adjacent to) the surface layer

!   -- NCON values are  propagated downwards from bottom of last active layer
!   -- PCON values also propagated downwards, BUT only present if surface condition
!  1.   Require linearized solutions at bottom of last active layer
!  2.   Set values for layer immediately below last active layer
!  3.   Remaining layers to bottom, just propagate using discrete-ordinate transmittances

      NAL = ACTIVE_LAYERS(NLAYERS_TEL) ; NAL1 = NAL + 1
      IF ( NAL .LT. NLAYERS ) THEN

!  N-constants, always required

!  start stream, stokes, parameter loops

        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_COLUMNWFS

!  real homogeneous solutions
!   Needs checking------------------------- 22 March 2007

              SHOM_R = ZERO
              DO K = 1, K_REAL(NAL)
                NXR  = NCON(K,NAL,Q) *   SOLA_XPOS(I,O1,K,NAL)
                PXR  = PCON(K,NAL,Q) *   SOLB_XNEG(I,O1,K,NAL)
                LXR  = LCON(K,NAL)   *   SOLA_XPOS(I,O1,K,NAL)
                LLXR = LCON(K,NAL)   * L_SOLA_XPOS(I,O1,K,NAL,Q)
                MLXR = MCON(K,NAL)   * L_SOLB_XNEG(I,O1,K,NAL,Q)
                L_HOM2 = PXR + MLXR
                L_HOM1 =   T_DELT_EIGEN(K,NAL) * ( NXR + LLXR ) +  L_T_DELT_EIGEN(K,NAL,Q) * LXR
                SHOM_R = SHOM_R + L_HOM1 + L_HOM2
              ENDDO

!  complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(NAL) + 1
              DO K = 1, K_COMPLEX(NAL)
                K0 = 2*K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1  =   NCON(K1,NAL,Q) *   SOLA_XPOS(I,O1,K1,NAL) - NCON(K2,NAL,Q) *   SOLA_XPOS(I,O1,K2,NAL)
                NXR2  =   NCON(K1,NAL,Q) *   SOLA_XPOS(I,O1,K2,NAL) + NCON(K2,NAL,Q) *   SOLA_XPOS(I,O1,K1,NAL)
                PXR1  =   PCON(K1,NAL,Q) *   SOLB_XNEG(I,O1,K1,NAL) - PCON(K2,NAL,Q) *   SOLB_XNEG(I,O1,K2,NAL)
                LXR1  =   LCON(K1,NAL) *   SOLA_XPOS(I,O1,K1,NAL) - LCON(K2,NAL) *   SOLA_XPOS(I,O1,K2,NAL)
                LXR2  =   LCON(K1,NAL) *   SOLA_XPOS(I,O1,K2,NAL) + LCON(K2,NAL) *   SOLA_XPOS(I,O1,K1,NAL)
                LLXR1  =   LCON(K1,NAL) * L_SOLA_XPOS(I,O1,K1,NAL,Q) - LCON(K2,NAL) * L_SOLA_XPOS(I,O1,K2,NAL,Q)
                LLXR2  =   LCON(K1,NAL) * L_SOLA_XPOS(I,O1,K2,NAL,Q) + LCON(K2,NAL) * L_SOLA_XPOS(I,O1,K1,NAL,Q)
                MLXR1  =   MCON(K1,NAL) * L_SOLB_XNEG(I,O1,K1,NAL,Q) - MCON(K2,NAL) * L_SOLB_XNEG(I,O1,K2,NAL,Q)
                L_HOM2CR = PXR1 + MLXR1
                L_HOM1CR = +   T_DELT_EIGEN(K1,NAL)   * ( NXR1 + LLXR1 ) -   T_DELT_EIGEN(K2,NAL)   * ( NXR2 + LLXR2 ) &
                           + L_T_DELT_EIGEN(K1,NAL,Q) *   LXR1           - L_T_DELT_EIGEN(K2,NAL,Q) *   LXR2
                SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
              ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

              SHOM = SHOM_R + SHOM_CR
              SPAR = L_WLOWER(I,O1,NAL,Q)
              NCON(IROW,NAL1,Q) = SPAR + SHOM

!  End loops

            ENDDO 
          ENDDO
        ENDDO

!  other layers to bottom of medium: propagate downwards.
!   Additional variation, since you are passing through varying layers.

        DO N = NAL + 2, NLAYERS
          N1 = N - 1
          DO I = 1, NSTREAMS
            IR = ( I - 1 ) * NSTOKES
            DO O1 = 1, NSTOKES
              IROW = IR + O1
!  Enhancement # 14, 6/27/16
              NCON(IROW,N,1:N_COLUMNWFS) =   T_DELT_DISORDS(I,N1) * NCON(IROW,N,1:N_COLUMNWFS) &
                                         + L_T_DELT_DISORDS(I,N1,1:N_COLUMNWFS) * LCON(IROW,N1)
            ENDDO
          ENDDO
        ENDDO

!  P-Constants need to be determined if there is a surface condition. Otherwise zero.

        IF ( DO_INCLUDE_SURFACE ) THEN

!  start stream, stokes and parameter loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IC = ( I - 1  ) * NSTOKES
            DO O1 = 1, NSTOKES
              ICOW = IC + O1
              DO Q = 1, N_COLUMNWFS

!  real homogeneous solutions

                SHOM_R = ZERO
                DO K = 1, K_REAL(NAL)
                  NXR  = NCON(K,NAL,Q) *   SOLA_XPOS(I1,O1,K,NAL)
                  PXR  = PCON(K,NAL,Q) *   SOLB_XNEG(I1,O1,K,NAL)
                  MXR  = MCON(K,NAL)   *   SOLB_XNEG(I1,O1,K,NAL)
                  LLXR = LCON(K,NAL)   * L_SOLA_XPOS(I1,O1,K,NAL,Q)
                  MLXR = MCON(K,NAL)   * L_SOLB_XNEG(I1,O1,K,NAL,Q)
                  L_HOM1 = NXR + LLXR
                  L_HOM2 =   T_DELT_EIGEN(K,NAL) * ( PXR + MLXR ) + L_T_DELT_EIGEN(K,NAL,Q) * MXR
                  SHOM_R = SHOM_R + L_HOM1 + L_HOM2
                ENDDO

!  complex homogeneous solutions

                SHOM_CR = ZERO
                KO1 = K_REAL(NAL) + 1
                DO K = 1, K_COMPLEX(NAL)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  NXR1  = NCON(K1,NAL,Q) *   SOLA_XPOS(I1,O1,K1,NAL) - NCON(K2,NAL,Q) *   SOLA_XPOS(I1,O1,K2,NAL)
                  PXR1  = PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL) - PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL)
                  PXR2  = PCON(K1,NAL,Q) *   SOLB_XNEG(I1,O1,K2,NAL) + PCON(K2,NAL,Q) *   SOLB_XNEG(I1,O1,K1,NAL)
                  MXR1  = MCON(K1,NAL) *   SOLB_XNEG(I1,O1,K1,NAL) - MCON(K2,NAL) *   SOLB_XNEG(I1,O1,K2,NAL)
                  MXR2  = MCON(K1,NAL) *   SOLB_XNEG(I1,O1,K2,NAL) + MCON(K2,NAL) *   SOLB_XNEG(I1,O1,K1,NAL)
                  LLXR1 = LCON(K1,NAL) * L_SOLA_XPOS(I1,O1,K1,NAL,Q) - LCON(K2,NAL) * L_SOLA_XPOS(I1,O1,K2,NAL,Q)
                  MLXR1 = MCON(K1,NAL) * L_SOLB_XNEG(I1,O1,K1,NAL,Q) - MCON(K2,NAL) * L_SOLB_XNEG(I1,O1,K2,NAL,Q)
                  MLXR2 = MCON(K1,NAL) * L_SOLB_XNEG(I1,O1,K2,NAL,Q) + MCON(K2,NAL) * L_SOLB_XNEG(I1,O1,K1,NAL,Q)
                  L_HOM1CR = NXR1 + LLXR1
                  L_HOM2CR =   T_DELT_EIGEN(K1,NAL)     * ( PXR1 + MLXR1 ) -   T_DELT_EIGEN(K2,NAL)   * ( PXR2 + MLXR2 ) &
                             + L_T_DELT_EIGEN(K1,NAL,Q) *   MXR1           - L_T_DELT_EIGEN(K2,NAL,Q) * MXR2
                  SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
                ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

                SHOM = SHOM_R + SHOM_CR
                SPAR = L_WUPPER(I1,O1,NAL,Q)
                PCON(ICOW,NAL1,Q) = SPAR + SHOM

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  other constants propagated

          DO N = NAL + 2, NLAYERS
            N1 = N - 1
            DO I = 1, NSTREAMS
              I1  = I + NSTREAMS
              IC = ( I - 1  ) * NSTOKES
              DO O1 = 1, NSTOKES
                ICOW = IC + O1
                DO Q = 1, N_COLUMNWFS
                  PCON(ICOW,N,Q) = ( PCON(ICOW,N1,Q) - L_T_DELT_DISORDS(I,N,Q) * MCON(ICOW,N) ) / T_DELT_DISORDS(I,N) 
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  Otherwise all P-constants are zero
        
        ELSE

          DO N = NAL1, NLAYERS
            DO I = 1, NSTREAMS
              I1  = I + NSTREAMS
              IC = ( I - 1  ) * NSTOKES
              DO O1 = 1, NSTOKES
                ICOW = IC + O1
                DO Q = 1, N_COLUMNWFS
                  PCON(ICOW,N,Q)   = ZERO
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  End general surface treatment

        ENDIF

!  End clause for non-active layers below telescoped problem

      ENDIF

!  finish

      RETURN
      END SUBROUTINE LC_BVPTEL_SOLUTION_MASTER

!

      SUBROUTINE LC_BVPTEL_SURFACE_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, MODIFIED_BOUNDARY, & ! Flags
        IBEAM, FOURIER, NSTOKES, NSTREAMS, N, NLAYERS, N_COLUMNWFS,   & ! Numbers
        SURFACE_FACTOR, QUAD_STRMWTS, LAMBERTIAN_ALBEDO, BRDF_F,      & ! Surface input
        MUELLER_INDEX, K_REAL, K_COMPLEX,                             & ! bookkeeping
        T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input transmittances
        SOLA_XPOS,  SOLB_XNEG, WLOWER, L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER, & ! RT Solutions
        R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output reflected solutions

!  Linearized surface reflectance terms, Telescoping
!    Version 2.8. Major extension to BRDFs with TELESCOPING. Implemented May 2016.

!  4/15/20. Version 2.8.2. BRDF and SLEAVE arrays are defined locally, each Fourier.

!  module, dimensions and numbers
!  4/15/20. Version 2.8.2. Remove MAXMOMENTS dimension.


      USE VLIDORT_PARS_m, Only : MAXLAYERS, MAXSTOKES, MAXSTREAMS, MAX_ATMOSWFS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTOKES_SQ, ZERO, ONE

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          MODIFIED_BOUNDARY

!  Numbers

      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS, N
      INTEGER, INTENT (IN) ::          N_COLUMNWFS

!  Surface inputs
!  Fourier components of BRDF, Incident--> reflected streams
!  4/15/20. Version 2.8.2. BRDF array defined locally, each Fourier, remove MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

!  Bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  Solutions

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  Linearized solutions

      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Discrete ordinate transmittances

      DOUBLE PRECISION, intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output - linearized surface-reflected solutions

      DOUBLE PRECISION, INTENT (OUT) :: R2_L_BEAM ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: R2_L_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: R2_L_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

!  cumulative tranmsittance (and linearization) from bottom of lowest active layer to surface
!     Calculated here, but  could be done earlier and passed in

      DOUBLE PRECISION, intent (out) :: CUMTRANS(MAXSTREAMS)
      DOUBLE PRECISION, intent (out) :: L_CUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      DOUBLE PRECISION  :: QCUMTRANS(MAXSTREAMS)
      DOUBLE PRECISION  :: L_QCUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)

!  help arrays

      DOUBLE PRECISION :: PV_W ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: HV_P ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: HV_M ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: PS_W ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: HS_P ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: HS_M ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

!  other variables

      INTEGER ::          I, J, O1, IB, Q, N1, M
      INTEGER ::          K, KO1, K0, K1, K2, MO2(MAXSTOKES)

      DOUBLE PRECISION :: H1R, H1I, H2R, H2I, KMULT
      DOUBLE PRECISION :: H1, H2, HP, HM, L_H1, L_H2, L_H1R, L_H2R, L_H1I, L_H2I
      DOUBLE PRECISION :: BEAM, L_BEAM, HOMP, HOMM, L_HOMP, L_HOMM
      DOUBLE PRECISION :: HOMPR, HOMMR, L_HOMPR, L_HOMMR, HOMPI, HOMMI, L_HOMPI, L_HOMMI
      DOUBLE PRECISION :: H1_S_CR, H2_S_CR, H1_S_CI, H2_S_CI

!  Initial section
!  ---------------

!  Always zero the result to start

      R2_L_BEAM = ZERO
      R2_L_HOMP = ZERO
      R2_L_HOMM = ZERO

      CUMTRANS   = ONE
      L_CUMTRANS = ZERO

!  Beam index and Fourier component

      IB = IBEAM
      M = FOURIER

!  Return if no surface contributions

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Return if Fourier component > 0 (Lambertian)

      IF ( DO_LAMBERTIAN_SURFACE .and. M.gt.0 ) RETURN

!  Cumulative transmittance and its linearization

      CUMTRANS(1:NSTREAMS) = ONE ; L_CUMTRANS(1:NSTREAMS,:) = ZERO
      DO N1 = NLAYERS, N+1, -1
         DO Q = 1, N_COLUMNWFS
           L_CUMTRANS(1:NSTREAMS,Q) = L_CUMTRANS(1:NSTREAMS,Q) * T_DELT_DISORDS  (1:NSTREAMS,N1) &
                                      + CUMTRANS(1:NSTREAMS)   * L_T_DELT_DISORDS(1:NSTREAMS,N1,Q)
         ENDDO
         CUMTRANS(1:NSTREAMS)  = CUMTRANS(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
      ENDDO

!  stored variables

      QCUMTRANS(1:NSTREAMS) = CUMTRANS(1:NSTREAMS) * QUAD_STRMWTS(1:NSTREAMS)
      DO Q = 1, N_COLUMNWFS
         L_QCUMTRANS(1:NSTREAMS,Q) = L_CUMTRANS(1:NSTREAMS,Q) * QUAD_STRMWTS(1:NSTREAMS)
      ENDDO

!  Set up Auxiliary arrays
!  -----------------------

!  Particular integral

      DO J = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          PS_W(J,O1) = WLOWER(J,O1,N) * QCUMTRANS(J)
          DO Q = 1, N_COLUMNWFS
            PV_W(J,O1,Q) = L_WLOWER(J,O1,N,Q) * QCUMTRANS(J) + WLOWER(J,O1,N) * L_QCUMTRANS(J,Q)
          ENDDO
        ENDDO
      ENDDO

!    Modified boundary condition: homogeneous parts

      IF ( MODIFIED_BOUNDARY ) THEN

!  start loops

        DO J = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  real homogeneous solution contributions

            DO K = 1, K_REAL(N)
              H1 = SOLA_XPOS(J,O1,K,N) * T_DELT_EIGEN(K,N) 
              H2 = SOLB_XNEG(J,O1,K,N)
              HS_P(J,O1,K) = QCUMTRANS(J)*H1
              HS_M(J,O1,K) = QCUMTRANS(J)*H2
              DO Q = 1, N_COLUMNWFS
                L_H1 = L_SOLA_XPOS(J,O1,K,N,Q) *   T_DELT_EIGEN(K,N) + &
                         SOLA_XPOS(J,O1,K,N)   * L_T_DELT_EIGEN(K,N,Q)
                L_H2 = L_SOLB_XNEG(J,O1,K,N,Q)
                HV_P(J,O1,K,Q) = QCUMTRANS(J) * L_H1 + L_QCUMTRANS(J,Q) * H1
                HV_M(J,O1,K,Q) = QCUMTRANS(J) * L_H2 + L_QCUMTRANS(J,Q) * H2
              ENDDO
            ENDDO

!  Complex homogeneous solution contributions

            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              H1R = SOLA_XPOS(J,O1,K1,N) * T_DELT_EIGEN(K1,N) - SOLA_XPOS(J,O1,K2,N) * T_DELT_EIGEN(K2,N)
              H1I = SOLA_XPOS(J,O1,K1,N) * T_DELT_EIGEN(K2,N) + SOLA_XPOS(J,O1,K2,N) * T_DELT_EIGEN(K1,N)
              H2R = SOLB_XNEG(J,O1,K1,N)
              H2I = SOLB_XNEG(J,O1,K2,N)
              HS_P(J,O1,K1) = QCUMTRANS(J) * H1R
              HS_P(J,O1,K2) = QCUMTRANS(J) * H1I
              HS_M(J,O1,K1) = QCUMTRANS(J) * H2R
              HS_M(J,O1,K2) = QCUMTRANS(J) * H1R
              DO Q = 1, N_COLUMNWFS
                L_H1R = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K1,N)   &
                      - L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K2,N)   &
                      +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K1,N,Q) &
                      -   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K2,N,Q)
                L_H1I = L_SOLA_XPOS(J,O1,K1,N,Q) *   T_DELT_EIGEN(K2,N)   &
                      + L_SOLA_XPOS(J,O1,K2,N,Q) *   T_DELT_EIGEN(K1,N)   &
                      +   SOLA_XPOS(J,O1,K1,N)   * L_T_DELT_EIGEN(K2,N,Q) &
                      +   SOLA_XPOS(J,O1,K2,N)   * L_T_DELT_EIGEN(K1,N,Q)
                L_H2R = L_SOLB_XNEG(J,O1,K1,N,Q)
                L_H2I = L_SOLB_XNEG(J,O1,K2,N,Q)
                HV_P(J,O1,K1,Q) = QCUMTRANS(J) * L_H1R + L_QCUMTRANS(J,Q) * H1R
                HV_P(J,O1,K2,Q) = QCUMTRANS(J) * L_H1I + L_QCUMTRANS(J,Q) * H1I
                HV_M(J,O1,K1,Q) = QCUMTRANS(J) * L_H2R + L_QCUMTRANS(J,Q) * H2R
                HV_M(J,O1,K2,Q) = QCUMTRANS(J) * L_H2I + L_QCUMTRANS(J,Q) * H2I
              ENDDO
            ENDDO

!  End loops

          ENDDO
        ENDDO

!  End modified boundary condition

      ENDIF

!  Integrated Downward reflection (Lambertian case)
!  ================================================

      if ( DO_LAMBERTIAN_SURFACE ) THEN

!  reflection

        KMULT = SURFACE_FACTOR * LAMBERTIAN_ALBEDO

!  only 1 component

        O1 = 1

!  Particular solution (only for the first Stokes component)

        BEAM = SUM (PS_W(1:NSTREAMS,O1) )
        DO Q = 1, N_COLUMNWFS
          L_BEAM = SUM(PV_W(1:NSTREAMS,O1,Q))
          R2_L_BEAM(1:NSTREAMS,O1,Q) = KMULT * ( L_BEAM * CUMTRANS(1:NSTREAMS) + BEAM * L_CUMTRANS(1:NSTREAMS,Q) )
        ENDDO

!  Homogeneous solutions for the modified condition

        IF ( MODIFIED_BOUNDARY ) THEN

!  Homogeneous real solutions

          DO K = 1, K_REAL(N)
            HOMP = SUM(HS_P(1:NSTREAMS,O1,K))
            HOMM = SUM(HS_M(1:NSTREAMS,O1,K))
            DO Q = 1, N_COLUMNWFS
              L_HOMP = SUM(HV_P(1:NSTREAMS,O1,K,Q))
              L_HOMM = SUM(HV_M(1:NSTREAMS,O1,K,Q))
              R2_L_HOMP(1:NSTREAMS,O1,K,Q) = KMULT * ( L_HOMP * CUMTRANS(1:NSTREAMS) + HOMP * L_CUMTRANS(1:NSTREAMS,Q) )
              R2_L_HOMM(1:NSTREAMS,O1,K,Q) = KMULT * ( L_HOMM * CUMTRANS(1:NSTREAMS) + HOMM * L_CUMTRANS(1:NSTREAMS,Q) )
            ENDDO
          ENDDO

!  Homogeneous Complex solutions

          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            HOMPR = SUM(HV_P(1:NSTREAMS,O1,K1,Q))
            HOMPI = SUM(HV_P(1:NSTREAMS,O1,K2,Q))
            HOMMR = SUM(HV_M(1:NSTREAMS,O1,K1,Q))
            HOMMI = SUM(HV_M(1:NSTREAMS,O1,K2,Q))
            DO Q = 1, N_COLUMNWFS
              L_HOMPR = SUM(HV_P(1:NSTREAMS,O1,K1,Q))
              L_HOMPI = SUM(HV_P(1:NSTREAMS,O1,K2,Q))
              L_HOMMR = SUM(HV_M(1:NSTREAMS,O1,K1,Q))
              L_HOMMI = SUM(HV_M(1:NSTREAMS,O1,K2,Q))
              R2_L_HOMP(1:NSTREAMS,O1,K1,Q) = KMULT * ( L_HOMPR * CUMTRANS(1:NSTREAMS) + HOMPR * L_CUMTRANS(1:NSTREAMS,Q) )
              R2_L_HOMP(1:NSTREAMS,O1,K2,Q) = KMULT * ( L_HOMPI * CUMTRANS(1:NSTREAMS) + HOMPI * L_CUMTRANS(1:NSTREAMS,Q) )
              R2_L_HOMM(1:NSTREAMS,O1,K1,Q) = KMULT * ( L_HOMMR * CUMTRANS(1:NSTREAMS) + HOMMR * L_CUMTRANS(1:NSTREAMS,Q) )
              R2_L_HOMM(1:NSTREAMS,O1,K2,Q) = KMULT * ( L_HOMMI * CUMTRANS(1:NSTREAMS) + HOMMI * L_CUMTRANS(1:NSTREAMS,Q) )
            ENDDO
          ENDDO

!  End modified boundary condition clause

        ENDIF

!  BRDF surface condition
!  ======================

!  Integrated Downward reflection (Calculation, Bidirectional case)
!     homogeneous and particular solutions.
!  4/15/20. Version 2.8.2. BRDF array defined locally, remove M=Fourier index

      ELSE

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            MO2(1:NSTOKES) = MUELLER_INDEX(O1,1:NSTOKES)

!  Particular solution

            BEAM = ZERO
            DO J = 1, NSTREAMS
              H1 = DOT_PRODUCT(PS_W(J,1:NSTOKES),BRDF_F(MO2(1:NSTOKES),J,I))
              BEAM = BEAM + H1
            ENDDO
            DO Q = 1, N_COLUMNWFS
              L_BEAM = ZERO
              DO J = 1, NSTREAMS
                L_H1 = DOT_PRODUCT(PV_W(J,1:NSTOKES,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                L_BEAM = L_BEAM + L_H1
              ENDDO
              R2_L_BEAM(I,O1,Q) = SURFACE_FACTOR * ( L_BEAM * CUMTRANS(I) + BEAM * L_CUMTRANS(I,Q) )
            ENDDO

!  Homogeneous solutions for the modified condition

            IF ( MODIFIED_BOUNDARY ) THEN

!  Homogeneous real solutions

              DO K = 1, K_REAL(N)
                HOMP = ZERO ; HOMM = ZERO
                DO J = 1, NSTREAMS
                  HP = DOT_PRODUCT(HS_P(J,1:NSTOKES,K),BRDF_F(MO2(1:NSTOKES),J,I))
                  HM = DOT_PRODUCT(HS_M(J,1:NSTOKES,K),BRDF_F(MO2(1:NSTOKES),J,I))
                  HOMP = HOMP + HP ; HOMM = HOMM + HM
                ENDDO
                DO Q = 1, N_COLUMNWFS
                  L_HOMP = ZERO ; L_HOMM = ZERO
                  DO J = 1, NSTREAMS
                    HP = DOT_PRODUCT(HV_P(J,1:NSTOKES,K,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                    HM = DOT_PRODUCT(HV_M(J,1:NSTOKES,K,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                    L_HOMP = L_HOMP + HP ; L_HOMM = L_HOMM + HM
                  ENDDO
                  R2_L_HOMP(I,O1,K,Q) = SURFACE_FACTOR * ( L_HOMP * CUMTRANS(I) + HOMP * L_CUMTRANS(I,Q) )
                  R2_L_HOMM(I,O1,K,Q) = SURFACE_FACTOR * ( L_HOMM * CUMTRANS(I) + HOMM * L_CUMTRANS(I,Q) )
                ENDDO
              ENDDO

!  homogeneous complex solutions

              KO1 = K_REAL(NLAYERS) + 1
              DO K = 1, K_COMPLEX(NLAYERS)
                K0 = 2*K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                HOMPR = ZERO ; HOMMR = ZERO ; HOMPI = ZERO ; HOMMI = ZERO
                DO J = 1, NSTREAMS
                  H1_S_CR = DOT_PRODUCT(HS_P(J,1:NSTOKES,K1),BRDF_F(MO2(1:NSTOKES),J,I))
                  H2_S_CR = DOT_PRODUCT(HS_M(J,1:NSTOKES,K1),BRDF_F(MO2(1:NSTOKES),J,I))
                  H1_S_CI = DOT_PRODUCT(HS_P(J,1:NSTOKES,K2),BRDF_F(MO2(1:NSTOKES),J,I))
                  H2_S_CI = DOT_PRODUCT(HS_M(J,1:NSTOKES,K2),BRDF_F(MO2(1:NSTOKES),J,I))
                  HOMPR = HOMPR + H1_S_CR
                  HOMMR = HOMMR + H2_S_CR
                  HOMPI = HOMPI + H1_S_CI
                  HOMMI = HOMMI + H2_S_CI
                ENDDO
                DO Q = 1, N_COLUMNWFS
                  L_HOMPR = ZERO ; L_HOMMR = ZERO ; L_HOMPI = ZERO ; L_HOMMI = ZERO
                  H1_S_CR = DOT_PRODUCT(HV_P(J,1:NSTOKES,K1,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                  H2_S_CR = DOT_PRODUCT(HV_M(J,1:NSTOKES,K1,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                  H1_S_CI = DOT_PRODUCT(HV_P(J,1:NSTOKES,K2,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                  H2_S_CI = DOT_PRODUCT(HV_M(J,1:NSTOKES,K2,Q),BRDF_F(MO2(1:NSTOKES),J,I))
                  L_HOMPR = L_HOMPR + H1_S_CR
                  L_HOMMR = L_HOMMR + H2_S_CR
                  L_HOMPI = L_HOMPI + H1_S_CI
                  L_HOMMI = L_HOMMI + H2_S_CI
                  R2_L_HOMP(I,O1,K1,Q) = SURFACE_FACTOR * ( L_HOMPR * CUMTRANS(I) + HOMPR * L_CUMTRANS(I,Q) )
                  R2_L_HOMM(I,O1,K1,Q) = SURFACE_FACTOR * ( L_HOMMR * CUMTRANS(I) + HOMMR * L_CUMTRANS(I,Q) )
                  R2_L_HOMP(I,O1,K2,Q) = SURFACE_FACTOR * ( L_HOMPI * CUMTRANS(I) + HOMPI * L_CUMTRANS(I,Q) )
                  R2_L_HOMM(I,O1,K2,Q) = SURFACE_FACTOR * ( L_HOMMI * CUMTRANS(I) + HOMMI * L_CUMTRANS(I,Q) )
                ENDDO
              ENDDO

!  End modified condition

            ENDIF

!  end stream and Stokes loops

          ENDDO
        ENDDO

!  End clause BRDF vs. Lambertian

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LC_BVPTEL_SURFACE_SETUP

!

      SUBROUTINE LC_BVPTEL_COLUMN_SETUP ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTBEAM, & ! Flags
        FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_COLUMNWFS,          & ! Basic Control Numbers
        NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, DO_LAYER_SCATTERING,    & ! Numbers and scattering
        N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,         & ! BVPTel Control
        MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR,       & ! Bookkeeping
        DELTAU_SLANT, L_DELTAU_VERT, DIRECT_BEAM,               & ! Optical and direct beam
        T_DELT_DISORDS, L_T_DELT_DISORDS,                       & ! Discrete ordinate transmittances
        BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, BVEC, WLOWER, & ! Beam and PI Solutions
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON,         & ! Homogeneous Solutions
        QUAD_STRMWTS, LAMBERTIAN_ALBEDO, BRDF_F,                & ! Surface inputs
        L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG,               & ! Linearized Homog Solutions
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_BVEC,             & ! Linearized Beam/PI solutions
        L_WLOWER, L_WUPPER, COLTEL2_WF, SCOL2_WF )                ! output solutions and Column vector

!  Column setup for the linearized telescoped BVP

!  4/15/20. Version 2.8.2. BRDF arrays are defined locally, each Fourier.

!  module, dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXLAYERS, MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ, MAXTOTAL, ZERO


      IMPLICIT NONE

!  Flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM

!  Numbers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_COLUMNWFS

!  More numbers

      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  Layer scattering

      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  BVP Tel control

      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::          NLAYERS_TEL
      INTEGER, INTENT (IN) ::          ACTIVE_LAYERS ( MAXLAYERS )

!  optical and direct beam

      DOUBLE PRECISION, INTENT (IN) :: DELTAU_SLANT  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DIRECT_BEAM   ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

!  discrete ordinate stuff

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR

!  homogeneous RTE solutions and inegration constants

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Beam and particular integral vector

      INTEGER, INTENT (IN) ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  Surface inputs
!    Fourier components of BRDF, incident quadrature streams, reflected quadrature streams
!  4/15/20. Version 2.8.2. BRDF array defined locally, each Fourier, remove MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Linearized homogeneous solutions

      DOUBLE PRECISION, INTENT (IN) :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Beam and PI vector

      DOUBLE PRECISION, INTENT (IN) :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (IN) :: LC_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  output particular integrals

      DOUBLE PRECISION, INTENT (OUT) :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  output Column vectos

      DOUBLE PRECISION, INTENT (OUT) :: COLTEL2_WF ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

!  local variables
!  ---------------

!  help variables

      LOGICAL ::          MODIFIED_BOUNDARY
      INTEGER ::          Q, N, I, I1, IR, CM, C0, IROW, O1, N1, NS
      INTEGER ::          K, KO1, K0, K1, K2, STATUS
      DOUBLE PRECISION :: CPOS, CNEG, L_HOM_R, L_HOM_CR, L_BEAM, BEAM
      DOUBLE PRECISION :: L_HOM_U_R, L_HOM_U_CR, L_HOM_D_R, L_HOM_D_CR
      DOUBLE PRECISION :: T1,T2,T1R,T1I,T2R,T2I,FAC3

!  Output arguments from the Surface setup (reflectances and cumulative transmittances)
!  cumulative tranmsittance (and linearization) from bottom of lowest active layer to surface

      DOUBLE PRECISION :: R2_L_BEAM ( MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: R2_L_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: CUMTRANS(MAXSTREAMS)
      DOUBLE PRECISION :: L_CUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)

!  status

      status = 0

!  Try this safety-first zeroing, copied from the scalar LIDORT code.
!    Very important to do this.   Bug solved, July 14th 2009

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
!  Enhancement # 15, 6/27/16
        L_WUPPER(1:NSTREAMS_2,1:NSTOKES,N,1:N_COLUMNWFS) = ZERO
        L_WLOWER(1:NSTREAMS_2,1:NSTOKES,N,1:N_COLUMNWFS) = ZERO
      ENDDO

!  Get the linearized solutions for the layer that is varying
!    Always need this, regardless of number of active layers

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        CALL LC_BEAMSOLUTION_NEQK ( &
          FOURIER, IBEAM, N, N_COLUMNWFS, NSTOKES,      & ! basic control
          NSTREAMS_2, DO_LAYER_SCATTERING, BEAM_CUTOFF, & ! Bookkeeping
          T_DELT_MUBAR, INITIAL_TRANS, BVEC,            & ! Beam/SolarPI nputs
          LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_BVEC,   & ! Linearized Beam/SolarPI 
          L_WUPPER, L_WLOWER )                            ! Output solutions
      ENDDO

!  General case. NLAYERS_TEL > 1
!  =============================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  zero column vector
!   Enhancement # 16, 6/27/16
        COLTEL2_WF(1:N_BVTELMATRIX_SIZE, 1:MAX_ATMOSWFS) = ZERO

!  top of first active layer, first boundary condition
!  ---------------------------------------------------

        NS = 1
        N = ACTIVE_LAYERS(NS)
        C0 = 0

!  require homogeneous and beam solution linearizations

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            DO Q = 1, N_COLUMNWFS

!  Beam contribution

              L_BEAM  = - L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

!  Enhancement # 17, 6/27/16. NOT IMPLEMENTED, LACKS CLARITY
!            L_HOM_R = sum( LCON(1:K_REAL(N),N) * L_SOLA_XPOS(I,O1,1:K_REAL(N),N,Q) + &
!                      MCON(1:K_REAL(N),N) * ( T_DELT_EIGEN(1:K_REAL(N),N) * L_SOLB_XNEG(I,O1,1:K_REAL(N),N,Q) + &
!                      L_T_DELT_EIGEN(1:K_REAL(N),N,Q) *   SOLB_XNEG(I,O1,1:K_REAL(N),N) )  )
              L_HOM_R  = ZERO
              DO K = 1, K_REAL(N)
                CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
                CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                     L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                T1 = LCON(K,N) * CPOS
                T2 = MCON(K,N) * CNEG
                L_HOM_R = L_HOM_R + T1 + T2
              ENDDO

!  Linearized Complex homogeneous solution contributions

              L_HOM_CR  = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2*K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                     L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
                T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                     - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                     + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                     - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                     + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                     + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                     + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
                T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                L_HOM_CR = L_HOM_CR + T1 + T2
              ENDDO

!  Final contribution

              COLTEL2_WF(CM,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  end loops

            ENDDO
          ENDDO
        ENDDO

!  Intermediate boundaries between active layers
!  ---------------------------------------------

        DO NS = 1, NLAYERS_TEL - 1

!  offsets

          N  = ACTIVE_LAYERS(NS)
          N1 = N + 1
          C0 = NS*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  Get the linearized beam solution for the next layer N1

          DO I = 1, NSTREAMS_2
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              DO Q = 1, N_COLUMNWFS

!  Beam contributions

                L_BEAM = L_WUPPER(I,O1,N1,Q) - L_WLOWER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions above

                L_HOM_U_R = ZERO
                DO K = 1, K_REAL(N1)
                  CPOS = L_SOLA_XPOS(I,O1,K,N1,Q)
                  CNEG = T_DELT_EIGEN(K,N1)   * L_SOLB_XNEG(I,O1,K,N1,Q) + &
                       L_T_DELT_EIGEN(K,N1,Q) *   SOLB_XNEG(I,O1,K,N1)
                  T1 = LCON(K,N1) * CPOS
                  T2 = MCON(K,N1) * CNEG
                  L_HOM_U_R = L_HOM_U_R + T1 + T2
                ENDDO

!  Linearized Real homogeneous solution contributions below

                L_HOM_D_R = ZERO
                DO K = 1, K_REAL(N)
                  CNEG = L_SOLB_XNEG(I,O1,K,N,Q)
                  CPOS = T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q) + &
                       L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                  T1 = LCON(K,N) * CPOS
                  T2 = MCON(K,N) * CNEG
                  L_HOM_D_R = L_HOM_D_R + T1 + T2
                ENDDO

!  Linearized Complex homogeneous solution contributions above

                L_HOM_U_CR  = ZERO
                KO1 = K_REAL(N1) + 1
                DO K = 1, K_COMPLEX(N1)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  T1 = L_SOLA_XPOS(I,O1,K1,N1,Q) * LCON(K1,N1) - &
                       L_SOLA_XPOS(I,O1,K2,N1,Q) * LCON(K2,N1)
                  T2R =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                       - T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                       + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K1,N1) &
                       - L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K2,N1)
                  T2I =  T_DELT_EIGEN(K1,N1)   * L_SOLB_XNEG(I,O1,K2,N1,Q) &
                       + T_DELT_EIGEN(K2,N1)   * L_SOLB_XNEG(I,O1,K1,N1,Q) &
                       + L_T_DELT_EIGEN(K1,N1,Q) *   SOLB_XNEG(I,O1,K2,N1) &
                       + L_T_DELT_EIGEN(K2,N1,Q) *   SOLB_XNEG(I,O1,K1,N1)
                  T2 =  T2R * MCON(K1,N1) - T2I * MCON(K2,N1)
                  L_HOM_U_CR = L_HOM_U_CR + T1 + T2
                ENDDO
    
!  Linearized Complex homogeneous solution contributions below

                L_HOM_D_CR  = ZERO
                KO1 = K_REAL(N) + 1
                DO K = 1, K_COMPLEX(N)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  T2 = L_SOLB_XNEG(I,O1,K1,N,Q) * MCON(K1,N) - &
                       L_SOLB_XNEG(I,O1,K2,N,Q) * MCON(K2,N)
                  T1R =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                       - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                       + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) &
                       - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                  T1I =  T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I,O1,K2,N,Q) &
                       + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I,O1,K1,N,Q) &
                       + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) &
                       + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                  T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                  L_HOM_D_CR = L_HOM_D_CR + T1 + T2
                ENDDO

!  Final contribution

                L_HOM_R  =  L_HOM_U_R  - L_HOM_D_R
                L_HOM_CR =  L_HOM_U_CR - L_HOM_D_CR
                COLTEL2_WF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  End loop over intermediate active layer boundaries

        ENDDO

!  Final boundary, bottom of lowest active layer
!  ---------------------------------------------

        NS = NLAYERS_TEL
        N  = ACTIVE_LAYERS(NS)
        C0 = (NS-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS
        MODIFIED_BOUNDARY = .true.

!  Old code. Condition is now completely general.
!  If this is the surface and Specialist option #2 is in place
!       if ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 &
!             .AND.  N.EQ.NLAYERS ) THEN

!  get the linearized downward-reflected term. New 6/29/16 Rob Fix
!  4/15/20. Version 2.8.2. BRDF and SLEAVE arrays are defined locally, each Fourier.

        CALL LC_BVPTEL_SURFACE_SETUP                                           &
          ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, MODIFIED_BOUNDARY,      & ! Flags
            IBEAM, FOURIER, NSTOKES, NSTREAMS, N, NLAYERS, N_COLUMNWFS,        & ! Numbers
            SURFACE_FACTOR, QUAD_STRMWTS, LAMBERTIAN_ALBEDO, BRDF_F,           & ! Surface input
            MUELLER_INDEX, K_REAL, K_COMPLEX,                                  & ! bookkeeping
            T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input transmittances
            SOLA_XPOS, SOLB_XNEG, WLOWER, L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER,  & ! RT Solutions
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output reflected solutions

!  last active layer, With surface
!  ###############################

        IF ( DO_INCLUDE_SURFACE ) THEN

!  Staert loops

          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              DO Q = 1, N_COLUMNWFS

!  Beam contributions

                L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

!  Linearized Real homogeneous solution contributions

                L_HOM_R  = ZERO
                DO K = 1, K_REAL(N)
                  CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                      + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                  CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
                  CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
                  CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
                  L_HOM_R = L_HOM_R + LCON(K,N) * CPOS + MCON(K,N) * CNEG
                ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

                L_HOM_CR  = ZERO
                KO1 = K_REAL(N) + 1
                DO K = 1, K_COMPLEX(N)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  T1R =   T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &  
                      -   T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)   & 
                      - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                  T1I =   T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) & 
                      +   T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)   &
                      + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                  T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
                  T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
                  T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                  T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
                  T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
                  T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                  L_HOM_CR = L_HOM_CR + T1 + T2
                ENDDO

!  Final contributions

                COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  last active layer, No surface
!  #############################

        ELSE

!  Start loops

          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              DO Q = 1, N_COLUMNWFS

!  Beam contributions

                L_BEAM = L_WLOWER(I1,O1,N,Q)

!  Linearized Real homogeneous solution contributions

                L_HOM_R  = ZERO
                DO K = 1, K_REAL(N)
                  CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                      + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                  CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
                  L_HOM_R = L_HOM_R + LCON(K,N) * CPOS + MCON(K,N) * CNEG
                ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

                L_HOM_CR  = ZERO
                KO1 = K_REAL(N) + 1
                DO K = 1, K_COMPLEX(N)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  T1R =   T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &  
                      -   T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N)   & 
                      - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                  T1I =   T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) & 
                      +   T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N)   &
                      + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                  T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                  T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) 
                  T2I = L_SOLB_XNEG(I1,O1,K2,N,Q)
                  T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                  L_HOM_CR = L_HOM_CR + T1 + T2
                ENDDO

!  Final contributions

                COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  End last active layer computation

        ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------

!  Add direct beam solution. This is new code, R. Spurr 06/29/16
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.
!    --- Formerly: Only for Specialist option 2, Fourier = 0 Lambertian, lowest active layer at surface

        IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              BEAM = DIRECT_BEAM(I,IBEAM,O1)
              DO Q = 1, N_COLUMNWFS
                L_BEAM = ZERO
                DO K = 1, NLAYERS
                  FAC3 = - BEAM * DELTAU_SLANT(NLAYERS,K,IBEAM)
                  L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K) * FAC3
                ENDDO
                COLTEL2_WF(CM,Q) = COLTEL2_WF(CM,Q) + L_BEAM * CUMTRANS(I) + L_CUMTRANS(I,Q) * BEAM
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  single-active-layer case
!  ========================

      ELSE

!  zero column vector. Extremely important.
!  Enhancement # 18, 6/27/16
        SCOL2_WF(1:NSTKS_NSTRMS_2, 1:MAX_ATMOSWFS) = ZERO

!  top of active layer
!  -------------------

        NS = 1
        N = ACTIVE_LAYERS(NS)

!  layer that is varying,
!       then require homogeneous and beam solution linearizations

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_COLUMNWFS

!  beam solution linearization at top of layer

              L_BEAM = - L_WUPPER(I,O1,N,Q)

!  Linearized Real homogeneous solution contributions

              L_HOM_R  = ZERO
              DO K = 1, K_REAL(N)
                CPOS = L_SOLA_XPOS(I,O1,K,N,Q)
                CNEG = T_DELT_EIGEN(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q) + &
                     L_T_DELT_EIGEN(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                L_HOM_R = L_HOM_R + LCON(K,N) * CPOS + MCON(K,N) * CNEG
              ENDDO

!  Linearized Complex homogeneous solution contributions

              L_HOM_CR  = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2*K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                T1 = L_SOLA_XPOS(I,O1,K1,N,Q) * LCON(K1,N) - &
                     L_SOLA_XPOS(I,O1,K2,N,Q) * LCON(K2,N)
                T2R =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                     - T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                     + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) &
                     - L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                T2I =  T_DELT_EIGEN(K1,N)   * L_SOLB_XNEG(I,O1,K2,N,Q) &
                     + T_DELT_EIGEN(K2,N)   * L_SOLB_XNEG(I,O1,K1,N,Q) &
                     + L_T_DELT_EIGEN(K1,N,Q) *   SOLB_XNEG(I,O1,K2,N) &
                     + L_T_DELT_EIGEN(K2,N,Q) *   SOLB_XNEG(I,O1,K1,N)
                T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                L_HOM_CR = L_HOM_CR + T1 + T2
              ENDDO

!  Final contribution

              SCOL2_WF(IROW,Q) = L_BEAM - L_HOM_R - L_HOM_CR

!  end loops

            ENDDO
          ENDDO
        ENDDO

!  Bottom of active layer
!  ----------------------

        MODIFIED_BOUNDARY = .true.
        C0 = NSTKS_NSTRMS

!  Old code. Condition is nowcompletely general.
!  If this is the surface and Specialist option #2 is in place
!       if ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 &
!             .AND.  N.EQ.NLAYERS ) THEN

!  get the linearized downward-reflected term. New 6/29/16 Rob Fix
!  4/15/20. Version 2.8.2. BRDF and SLEAVE arrays are defined locally, each Fourier.

        CALL LC_BVPTEL_SURFACE_SETUP                                           &
          ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, MODIFIED_BOUNDARY,      & ! Flags
            IBEAM, FOURIER, NSTOKES, NSTREAMS, N, NLAYERS, N_COLUMNWFS,        & ! Numbers
            SURFACE_FACTOR, QUAD_STRMWTS, LAMBERTIAN_ALBEDO, BRDF_F,           & ! Surface input
            MUELLER_INDEX, K_REAL, K_COMPLEX,                                  & ! bookkeeping
            T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input transmittances
            SOLA_XPOS, SOLB_XNEG, WLOWER, L_SOLA_XPOS, L_SOLB_XNEG, L_WLOWER,  & ! RT Solutions
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output reflected solutions

!  Last Active layer with surface
!  ##############################

        IF ( DO_INCLUDE_SURFACE ) THEN

!  Start loops

          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              DO Q = 1, N_COLUMNWFS

!  Beam contributions

                L_BEAM = L_WLOWER(I1,O1,N,Q) - R2_L_BEAM(I,O1,Q)

!  Linearized Real homogeneous solution contributions

                L_HOM_R  = ZERO
                DO K = 1, K_REAL(N)
                  CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                      + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                  CPOS = CPOS - R2_L_HOMP(I,O1,K,Q)
                  CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
                  CNEG = CNEG - R2_L_HOMM(I,O1,K,Q)
                  L_HOM_R = L_HOM_R + LCON(K,N) * CPOS + MCON(K,N) * CNEG
                ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

                L_HOM_CR  = ZERO
                KO1 = K_REAL(N) + 1
                DO K = 1, K_COMPLEX(N)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                      - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                      - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                  T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                      + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                      + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                  T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
                  T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
                  T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                  T2R = L_SOLB_XNEG(I1,O1,K1,N,Q) - R2_L_HOMM(I,O1,K1,Q)
                  T2I = L_SOLB_XNEG(I1,O1,K2,N,Q) - R2_L_HOMM(I,O1,K2,Q)
                  T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                  L_HOM_CR = L_HOM_CR + T1 + T2
                ENDDO

!  Final contributions

                SCOL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  Last Active layer, No surface
!  #############################

        ELSE

!  Start loops

          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              DO Q = 1, N_COLUMNWFS

!  Beam contributions

                L_BEAM = L_WLOWER(I1,O1,N,Q) 

!  Linearized Real homogeneous solution contributions

                L_HOM_R  = ZERO
                DO K = 1, K_REAL(N)
                  CPOS =  T_DELT_EIGEN(K,N)   * L_SOLA_XPOS(I1,O1,K,N,Q) &
                      + L_T_DELT_EIGEN(K,N,Q) *   SOLA_XPOS(I1,O1,K,N)
                  CNEG = L_SOLB_XNEG(I1,O1,K,N,Q)
                  L_HOM_R = L_HOM_R + LCON(K,N) * CPOS + MCON(K,N) * CNEG
                ENDDO

!  Linearized Complex homogeneous solution contributions
!    Bug Fixed 16 December 2005.

                L_HOM_CR  = ZERO
                KO1 = K_REAL(N) + 1
                DO K = 1, K_COMPLEX(N)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  T1R = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                      - T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K1,N) &
                      - L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K2,N)
                  T1I = T_DELT_EIGEN(K1,N)   * L_SOLA_XPOS(I1,O1,K2,N,Q) &
                      + T_DELT_EIGEN(K2,N)   * L_SOLA_XPOS(I1,O1,K1,N,Q) &
                      + L_T_DELT_EIGEN(K1,N,Q) *   SOLA_XPOS(I1,O1,K2,N) &
                      + L_T_DELT_EIGEN(K2,N,Q) *   SOLA_XPOS(I1,O1,K1,N)
                  T1R = T1R - R2_L_HOMP(I,O1,K1,Q)
                  T1I = T1I - R2_L_HOMP(I,O1,K2,Q)
                  T1 =  T1R * LCON(K1,N) - T1I * LCON(K2,N)
                  T2R = L_SOLB_XNEG(I1,O1,K1,N,Q)
                  T2I = L_SOLB_XNEG(I1,O1,K2,N,Q)
                  T2 =  T2R * MCON(K1,N) - T2I * MCON(K2,N)
                  L_HOM_CR = L_HOM_CR + T1 + T2
                ENDDO

!  Final contributions

                SCOL2_WF(CM,Q) = - L_BEAM - L_HOM_R - L_HOM_CR

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  End surface inclusion cause

        ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------

!  Add direct beam solution. This is new code, R. Spurr 06/29/16
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.
!    --- Formerly: Only for Specialist option 2, Fourier = 0 Lambertian, lowest active layer at surface

        IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              BEAM = DIRECT_BEAM(I,IBEAM,O1)
              DO Q = 1, N_COLUMNWFS
                L_BEAM = ZERO
                DO K = 1, NLAYERS
                  FAC3 = - BEAM * DELTAU_SLANT(NLAYERS,K,IBEAM)
                  L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K) * FAC3
                ENDDO
                SCOL2_WF(CM,Q) = SCOL2_WF(CM,Q) + L_BEAM * CUMTRANS(I) + L_CUMTRANS(I,Q) * BEAM
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  End layer clause

      ENDIF

!  finish

      RETURN
      END SUBROUTINE LC_BVPTEL_COLUMN_SETUP

!  End Module

      END MODULE vlidort_lc_bvproblem_m

