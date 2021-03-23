
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
! #            VLIDORT_SETUP_MASTER (top-level master)          #
! #                                                             #
! #            VLIDORT_Bookkeep_Init                            #
! #            VLIDORT_MSGeom_Init                              #
! #            VLIDORT_FOGeom_Init                              #
! #                                                             #
! ###############################################################

!  VLIDORT HISTORY NOTES
!  ====================

!  Notes for Version 3.5
!  ---------------------

!  Construction August 26-31, 2010
!  R. Spurr, RT SOLUTIONS Inc., 9 Channing Street, Cambridge MA 02138

!  Notes for Version 2.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module VLIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter

!  Notes for Version 2.8
!  ---------------------

!       - New inputs for Phase matrices
!       - New interface for FO code, old SS corrections removed
!       - BVP Telescoping and reflecting surfaces.
!       - Updates on Water-leaving and BRDF kernels
!       - Bookkeeping overhaul

!  Notes for Version v.8.1
!  -----------------------

!       - Revised Water-leaving ocean-optics
!       - Complete new treatment of water-leaving flux-adjustment
!       - Media-property module introduced.
!       - Planetary problem output (S and T) generated.
!       - TOA/BOA illumination control.

MODULE VLIDORT_Setup_Master_m

!  Module for Doing setup operations for VLIDORT 2.8.2
!    - Checking, Bookkeeping and Geometry
!    - mick mod: added subroutines VLIDORT_Bookkeep_Init, VLIDORT_MSGeom_Init, VLIDORT_FOGeom_Init

!  Parameter types

      USE VLIDORT_pars_m

!  I/O Type structure dependencies

      USE VLIDORT_Inputs_def_m
      USE VLIDORT_Outputs_def_m
      USE VLIDORT_Setups_def_m

!  VLIDORT module dependencies
!    "Inputs" for the Checking and Bookkeeping, and MS geometry
!    "FOGeometry" is a separate module....

      USE VLIDORT_inputs_m,     only : VLIDORT_CHECK_INPUT, VLIDORT_DERIVE_INPUT, &
                                       VLIDORT_CHAPMAN_INPUT

      USE VFO_Geometry_Master_m

!  Public

      PUBLIC  :: VLIDORT_Setup_Master, &
                 VLIDORT_Bookkeep_Init, &
                 VLIDORT_MSGeom_Init, &
                 VLIDORT_FOGeom_Init

      CONTAINS

SUBROUTINE VLIDORT_Setup_Master ( &
         VLIDORT_FixIn, VLIDORT_ModIn, & ! INPUTS (possibly modified)
         VLIDORT_Bookkeep, VLIDORT_FOGeom, VLIDORT_MSGeom, VLIDORT_Out )

!  Implicit none

      IMPLICIT NONE

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs)   , INTENT (IN)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (INOUT) :: VLIDORT_ModIn

!  VLIDORT Type structures, Bookkeeping and FO-Geometry outputs

      TYPE(VLIDORT_Bookkeeping), INTENT(INOUT)      :: VLIDORT_Bookkeep

!  Type structures, FO and MS Geometry outputs
  
      TYPE(VLIDORT_Geometry_FO), INTENT(INOUT)      :: VLIDORT_FOGeom
      TYPE(VLIDORT_Geometry_MS), INTENT(INOUT)      :: VLIDORT_MSGeom

!  Exception handling

      TYPE(VLIDORT_Outputs), INTENT (INOUT)         :: VLIDORT_Out

!  Exception handling
!  ==================

!  Exception handling for Input Checking. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER              :: STATUS_INPUTCHECK
      INTEGER              :: NCHECKMESSAGES
      CHARACTER(Len=120)   :: CHECKMESSAGES(0:MAX_MESSAGES)
      CHARACTER(Len=120)   :: ACTIONS (0:MAX_MESSAGES)

!  Exception handling for Model Calculation. New code, 18 May 2010

      LOGICAL              :: FAIL
      INTEGER              :: STATUS_CALCULATION
      CHARACTER(Len=120)   :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

!  Checking
!  --------    

!  %% Major revision of I/O output list, 25 October 2012.
!     ---- Observational Geometry control, New, 25 October 2012
!     ---- Automatic setting of NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, DO_USER_STREAMS, DO_NO_AZIMUTH

!  %% DO_FULLRAD_MODE argument added (Second line). R. Spurr, 05 March 2013
!  %%   Needed to ensure MS-only output in all cases when flagged

!  %% Revision, Version 2.7, 10/10/13. Include TAYLOR_ORDER parameter, revise argument list
!  %% Revision, Version 2.8, 03/10/17. Include TF variables, revise argument list
!mick fix 3/22/2017 - replaced USER_ANGLES with USER_ANGLES_INPUT in the argument list
!  Add Flag DO_WLADJUSTED_OUTPUT (Water-leaving output). 4/22/19 for Version 3.8.1. 

!  Initialize

      CHECKMESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS      (1:MAX_MESSAGES) = ' '
      NCHECKMESSAGES  = 0
      CHECKMESSAGES(0) = 'Successful Check of VLIDORT Basic Input'
      ACTIONS(0)       = 'No Action required for this Task'

!  Check

      CALL VLIDORT_CHECK_INPUT &
        ( VLIDORT_FixIn, VLIDORT_ModIn,                             & ! InOut
          STATUS_INPUTCHECK, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS ) ! Output

!  Exception handling

      IF ( STATUS_INPUTCHECK .EQ. VLIDORT_SERIOUS ) THEN
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VLIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        VLIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        VLIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ELSE IF ( STATUS_INPUTCHECK .EQ. VLIDORT_WARNING ) THEN
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VLIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        VLIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        VLIDORT_Out%Status%TS_ACTIONS           = ACTIONS
      ENDIF

!  Initialize VLIDORT derived inputs
!  ---------------------------------

      CALL VLIDORT_Bookkeep_Init ( VLIDORT_Bookkeep )
      CALL VLIDORT_MSGeom_Init   ( VLIDORT_MSGeom )
      IF ( VLIDORT_ModIn%MBool%TS_DO_FOCORR ) CALL VLIDORT_FOGeom_Init ( VLIDORT_FOGeom )

!  Get derived inputs
!  ------------------

      CALL VLIDORT_DERIVE_INPUT &
      ( VLIDORT_FixIn, VLIDORT_ModIn,      & ! Input type structures
        VLIDORT_Bookkeep, VLIDORT_MSGeom )   ! Output Type structures

!  Chapman function calculation, MS Only
!  -------------------------------------

      IF ( VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION .and. VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES ) THEN

         CALL VLIDORT_CHAPMAN_INPUT &
            ( VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Bookkeep, & ! Input type structures
              VLIDORT_MSGeom, FAIL, MESSAGE, TRACE_1 )        ! Output TS and exception handling

         IF (FAIL) THEN
            TRACE_2 = 'Failure from VLIDORT_CHAPMAN_INPUT, Called in VLIDORT_SETUP_MASTER'
            STATUS_CALCULATION = VLIDORT_SERIOUS
            VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
            VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            VLIDORT_Out%Status%TS_TRACE_3 = TRACE_2
            RETURN
         ENDIF

      ENDIF

!  First Order Geometry Calculation
!  --------------------------------

      IF ( VLIDORT_ModIn%MBool%TS_DO_FOCORR ) THEN

         CALL VFO_GEOMETRY_MASTER &
           ( VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Bookkeep,   & ! VLIDORT Input and Bookkeeping
             VLIDORT_FOGeom, FAIL, MESSAGE, TRACE_1, TRACE_2 ) ! Exception-Handling

         IF (FAIL) THEN
            TRACE_3 = 'Failure from VFO_GEOMETRY_MASTER, Called in VLIDORT_SETUP_MASTER'
            STATUS_CALCULATION = VLIDORT_SERIOUS
            VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
            VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
         ENDIF

      ENDIF

!  End of routine

END SUBROUTINE VLIDORT_Setup_Master

!

SUBROUTINE VLIDORT_Bookkeep_Init ( VLIDORT_Bookkeep )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_Setups_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT derived input structure

      TYPE(VLIDORT_Bookkeeping), INTENT(INOUT) :: VLIDORT_Bookkeep

!  Initialize VLIDORT bookkeep inputs
!  =================================

!  Mode of operation

      VLIDORT_Bookkeep%DO_MSMODE_VLIDORT = .FALSE.
      VLIDORT_Bookkeep%DO_MSMODE_THERMAL = .FALSE.
      VLIDORT_Bookkeep%DO_FOCORR_ALONE   = .FALSE.

!  Actual number of moments used in calculations

      VLIDORT_Bookkeep%NMOMENTS = 0

!  Total number of streams and layers
!  Number of super and sub diagonals in Band Matrix storage

      VLIDORT_Bookkeep%NSTREAMS_2     = 0
      VLIDORT_Bookkeep%NSTKS_NSTRMS   = 0
      VLIDORT_Bookkeep%NSTKS_NSTRMS_2 = 0

      VLIDORT_Bookkeep%NTOTAL    = 0
      VLIDORT_Bookkeep%N_SUBDIAG = 0
      VLIDORT_Bookkeep%N_SUPDIAG = 0

!  Number of directions (1 or 2) and directional array

      VLIDORT_Bookkeep%N_DIRECTIONS     = 0
      VLIDORT_Bookkeep%WHICH_DIRECTIONS = 0

!  Output optical depth masks and indices

      VLIDORT_Bookkeep%PARTLAYERS_OUTFLAG  = .FALSE.
      VLIDORT_Bookkeep%PARTLAYERS_OUTINDEX = 0
      VLIDORT_Bookkeep%UTAU_LEVEL_MASK_UP  = 0
      VLIDORT_Bookkeep%UTAU_LEVEL_MASK_DN  = 0

!  Off-grid optical depths (values, masks, indices)

      VLIDORT_Bookkeep%DO_PARTLAYERS       = .FALSE.
      VLIDORT_Bookkeep%N_PARTLAYERS        = 0
      VLIDORT_Bookkeep%PARTLAYERS_LAYERIDX = 0
      VLIDORT_Bookkeep%PARTLAYERS_VALUES   = ZERO

!  Partial heights

      VLIDORT_Bookkeep%PARTLAYERS_HEIGHTS = ZERO

!  Number of convergences

      VLIDORT_Bookkeep%N_OUT_STREAMS = 0
      VLIDORT_Bookkeep%N_CONVTESTS   = 0

!  Layer masks and limits for doing integrated source terms

      VLIDORT_Bookkeep%STERM_LAYERMASK_UP = .FALSE.
      VLIDORT_Bookkeep%STERM_LAYERMASK_DN = .FALSE.

      VLIDORT_Bookkeep%N_ALLLAYERS_UP = 0
      VLIDORT_Bookkeep%N_ALLLAYERS_DN = 0

!  Offsets for geometry indexing

      VLIDORT_Bookkeep%LOCAL_N_USERAZM = 0
      VLIDORT_Bookkeep%N_GEOMETRIES    = 0
      VLIDORT_Bookkeep%SZA_OFFSETS     = 0
      VLIDORT_Bookkeep%VZA_OFFSETS     = 0

!  Post-processing masks

      VLIDORT_Bookkeep%N_PPSTREAMS   = 0
      VLIDORT_Bookkeep%PPSTREAM_MASK = 0

!  Mueller index

      VLIDORT_Bookkeep%MUELLER_INDEX = 0

!  Greekmat indices, DMAT, FLUXVEC and related quantities (VLIDORT only)

      VLIDORT_Bookkeep%GREEKMAT_INDEX = 0
      VLIDORT_Bookkeep%FLUXVEC        = ZERO
      VLIDORT_Bookkeep%DFLUX          = ZERO
      VLIDORT_Bookkeep%DMAT           = ZERO

!  Number of particular solutions (not enabled yet)

      VLIDORT_Bookkeep%NPARTICSOLS        = 0
      VLIDORT_Bookkeep%DMAT_PSOLS         = ZERO
      VLIDORT_Bookkeep%DMAT_PSOLS_FLUXVEC = ZERO

!  Misc control flags

      VLIDORT_Bookkeep%DO_ALL_FOURIER  = .FALSE.
      VLIDORT_Bookkeep%DO_DBCORRECTION = .FALSE.

!  Finish

END SUBROUTINE VLIDORT_Bookkeep_Init

!

SUBROUTINE VLIDORT_MSGeom_Init ( VLIDORT_MSGeom )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_Setups_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT derived input structure

      TYPE(VLIDORT_Geometry_MS), INTENT(INOUT) :: VLIDORT_MSGeom

!  Initialize VLIDORT MS geometry inputs
!  ====================================

!  Local solar zenith angles Cosines (regular case)

      VLIDORT_MSGeom%COS_SZANGLES = ZERO
      VLIDORT_MSGeom%SIN_SZANGLES = ZERO

!  Quadrature weights and abscissae, and product

      VLIDORT_MSGeom%QUAD_STREAMS = ZERO
      VLIDORT_MSGeom%QUAD_HALFWTS = ZERO
      VLIDORT_MSGeom%QUAD_WEIGHTS = ZERO
      VLIDORT_MSGeom%QUAD_STRMWTS = ZERO
      VLIDORT_MSGeom%QUAD_SINES   = ZERO
      VLIDORT_MSGeom%QUAD_ANGLES  = ZERO

!  User stream cosines, secants, angles

      VLIDORT_MSGeom%USER_STREAMS = ZERO
      VLIDORT_MSGeom%USER_SECANTS = ZERO
      VLIDORT_MSGeom%USER_SINES   = ZERO
      VLIDORT_MSGeom%USER_ANGLES  = ZERO

!  Chapman Factors (Level and partials),
!  Local solar zenith angles at nadir, average cosines

      VLIDORT_MSGeom%CHAPMAN_FACTORS  = ZERO
      VLIDORT_MSGeom%PARTIAL_CHAPFACS = ZERO

      VLIDORT_MSGeom%SZA_LEVEL_OUTPUT = ZERO
      VLIDORT_MSGeom%SUNLAYER_COSINES = ZERO

!  Finish

END SUBROUTINE VLIDORT_MSGeom_Init

!

SUBROUTINE VLIDORT_FOGeom_Init ( VLIDORT_FOGeom )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_Setups_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT derived input structure

      TYPE(VLIDORT_Geometry_FO), INTENT(INOUT) :: VLIDORT_FOGeom

!  Initialize VLIDORT FO geometry inputs
!  ====================================

!  SOLAR Variables
!  ===============

!  Flag for the Nadir case

      VLIDORT_FOGeom%doNadir = .FALSE.

!  Ray constants

      VLIDORT_FOGeom%Raycon_up = ZERO
      VLIDORT_FOGeom%Raycon_dn = ZERO

!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      VLIDORT_FOGeom%Mu0_up = ZERO
      VLIDORT_FOGeom%Mu1_up = ZERO
      VLIDORT_FOGeom%Mu0_dn = ZERO
      VLIDORT_FOGeom%Mu1_dn = ZERO

!  Cosine scattering angles
   
      VLIDORT_FOGeom%COSSCAT_UP = ZERO
      VLIDORT_FOGeom%COSSCAT_DN = ZERO

!  Legendre function for Phase function calculation

      VLIDORT_FOGeom%GENSPHER_UP  = ZERO
      VLIDORT_FOGeom%GENSPHER_DN  = ZERO
      VLIDORT_FOGeom%ROTATIONS_UP = ZERO
      VLIDORT_FOGeom%ROTATIONS_DN = ZERO

!  LOS Quadratures for Enhanced PS

      VLIDORT_FOGeom%nfinedivs      = 0
      VLIDORT_FOGeom%xfine          = ZERO
      VLIDORT_FOGeom%wfine          = ZERO

      VLIDORT_FOGeom%nfinedivs_p_up = 0
      VLIDORT_FOGeom%xfine_p_up     = ZERO
      VLIDORT_FOGeom%wfine_p_up     = ZERO

      VLIDORT_FOGeom%nfinedivs_p_dn = 0
      VLIDORT_FOGeom%xfine_p_dn     = ZERO
      VLIDORT_FOGeom%wfine_p_dn     = ZERO

!  LOS path lengths

      VLIDORT_FOGeom%LosW_paths = ZERO
      VLIDORT_FOGeom%LosP_paths = ZERO

!  Chapman factor outputs

      VLIDORT_FOGeom%chapfacs   = ZERO
      VLIDORT_FOGeom%chapfacs_p = ZERO

!  Level-boundary solar paths

      VLIDORT_FOGeom%ntraverse_up     = 0
      VLIDORT_FOGeom%sunpaths_up      = ZERO

      VLIDORT_FOGeom%ntraversefine_up = 0
      VLIDORT_FOGeom%sunpathsfine_up  = ZERO

      VLIDORT_FOGeom%ntraverse_dn     = 0
      VLIDORT_FOGeom%sunpaths_dn      = ZERO

      VLIDORT_FOGeom%ntraversefine_dn = 0
      VLIDORT_FOGeom%sunpathsfine_dn  = ZERO

!  Partial-level outputs (Sunpaths, fine-sunpaths)

      VLIDORT_FOGeom%ntraverse_p_up     = 0
      VLIDORT_FOGeom%sunpaths_p_up      = ZERO

      VLIDORT_FOGeom%ntraversefine_p_up = 0
      VLIDORT_FOGeom%sunpathsfine_p_up  = ZERO

      VLIDORT_FOGeom%ntraverse_p_dn     = 0
      VLIDORT_FOGeom%sunpaths_p_dn      = ZERO

      VLIDORT_FOGeom%ntraversefine_p_dn = 0
      VLIDORT_FOGeom%sunpathsfine_p_dn  = ZERO

!  THERMAL GEOMETRY variables
!  ==========================

!  Geometry
 
      VLIDORT_FOGeom%Mu1_LOS        = ZERO
      VLIDORT_FOGeom%LosW_paths_LOS = ZERO
      VLIDORT_FOGeom%LosP_paths_LOS = ZERO

!  LOS Quadratures for Enhanced PS

      VLIDORT_FOGeom%nfinedivs_LOS      = 0
      VLIDORT_FOGeom%xfine_LOS          = ZERO
      VLIDORT_FOGeom%wfine_LOS          = ZERO

      VLIDORT_FOGeom%nfinedivs_p_LOS_up = 0
      VLIDORT_FOGeom%xfine_p_LOS_up     = ZERO
      VLIDORT_FOGeom%wfine_p_LOS_up     = ZERO

      VLIDORT_FOGeom%nfinedivs_p_LOS_dn = 0
      VLIDORT_FOGeom%xfine_p_LOS_dn     = ZERO
      VLIDORT_FOGeom%wfine_p_LOS_dn     = ZERO

!  Finish

END SUBROUTINE VLIDORT_FOGeom_Init

!  End of Module
     
end Module VLIDORT_Setup_Master_m


