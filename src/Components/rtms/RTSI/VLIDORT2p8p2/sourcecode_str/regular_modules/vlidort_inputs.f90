
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
! #            VLIDORT_INPUT_MASTER (master), calls             #
! #             VLIDORT_INIT_INPUTS                             #
! #             VLIDORT_READ_INPUTS                             #
! #                                                             #
! #            VLIDORT_Sup_Init                                 #
! #            VLIDORT_BRDF_Sup_Init                            #
! #            VLIDORT_SLEAVE_Sup_Init                          #
! #            VLIDORT_SS_Sup_Init                              #
! #                                                             #
! #    These routines are called by the VLIDORT SETUP MASTER    #
! #                                                             #
! #            VLIDORT_CHECK_INPUT                              #
! #            VLIDORT_DERIVE_INPUT                             #
! #            VLIDORT_CHAPMAN_INPUT                            #
! #                                                             #
! #    These routines are called by the Main VLIDORT masters    #
! #                                                             #
! #            VLIDORT_CHECK_INPUT_DIMS                         #
! #            VLIDORT_CHECK_INPUT_OPTICAL                      #
! #                                                             #
! ###############################################################

!  Upgrade to Version 2.8.1, September 2019
!    --- New inputs to control WLeaving and Planetary problems
!    --- New inputs to control TOA/BOA isotropic illumination added, 3/23/19

!  4/15/20. Version 2.8.2. Changes as follows
!     1. New Module VLIDORT_Setups_def_m, with Type Structures Bookkeep, Geometry_MS
!     2. DERIVE_INPUT subroutine now fills these structures.
!     3. Check input routine uses input type structures directly.
!     4. VLIDORT_CHAPMAN_INPUT subroutine introduced (formerly in vlidort_geometry.f90)
!     5. Doublet geometry post-processing operation.

      MODULE vlidort_inputs_m

!  Parameter types

      USE VLIDORT_pars_m

!  I/O type definitions

      USE VLIDORT_Inputs_def_m
      USE VLIDORT_Outputs_def_m
      USE VLIDORT_Setups_def_m

!  Auxiliary subroutine

      USE vlidort_aux_m, only : GFINDPAR, FINDPAR_ERROR, LEN_STRING, &
                                GETQUAD2, RSSORT, RQSORT_IDX

      !PRIVATE
      PUBLIC

      CONTAINS

      SUBROUTINE VLIDORT_INPUT_MASTER ( &
        FILNAM,             & ! INPUT
        VLIDORT_FixIn,      & ! OUTPUTS
        VLIDORT_ModIn,      & ! OUTPUTS
        VLIDORT_InputStatus ) ! OUTPUTS
        
!   Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination
!  --- Also include value of this illumination

!  Parameter types

      USE VLIDORT_PARS_m, Only : MAX_SZANGLES, MAX_USER_RELAZMS, MAX_USER_VZANGLES, &
                                 MAX_USER_LEVELS, MAX_USER_OBSGEOMS, MAXLAYERS, MAX_MESSAGES, &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, VLIDORT_INUNIT

      IMPLICIT NONE

!  Inputs
!  ------

      CHARACTER (LEN=*), intent(in) :: FILNAM

!  Outputs
!  -------

      TYPE(VLIDORT_Fixed_Inputs)            , INTENT (OUT) :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs)         , INTENT (OUT) :: VLIDORT_ModIn
      TYPE(VLIDORT_Input_Exception_Handling), INTENT (OUT) :: VLIDORT_InputStatus

!  Local variables
!  ===============

      INTEGER       :: FILUNIT 
      INTEGER       :: ISTAT, STATUS_SUB 

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER       :: STATUS
      INTEGER       :: NMESSAGES
      CHARACTER*120 :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*120 :: ACTIONS (0:MAX_MESSAGES)

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Initialize variables
!     Version 2.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 2.6: New: Observation Geometry variables, 25 October 2012
!     Version 2.8: DO_WATER_LEAVING and DO_FLUORESCENCE flags (added (1/31/16)
!     Version 2.8: Transflux iteration control for the WATER_LEAVING case. 3 "TF" variables added 7/6/16.
!     Version 2.8: FO/SS flags upgraded. 9/19/16, 3/1/17. Argument list rearranged.
!                  VLIDORT input type structures now used for argument passing - 9/19/2017
!                  All VLIDORT fixed and modified inputs now initialized - 9/19/2017

      CALL VLIDORT_INIT_INPUTS ( &
         VLIDORT_FixIn, VLIDORT_ModIn ) !Outputs

!  Open VLIDORT config file

      FILUNIT = VLIDORT_INUNIT
      OPEN(VLIDORT_INUNIT,FILE=FILNAM,IOSTAT=ISTAT,STATUS='OLD')

!  Open file error

      IF (ISTAT .GT. 0) THEN

!  Define the exception

         STATUS = VLIDORT_SERIOUS
         NMESSAGES = NMESSAGES + 1
         MESSAGES(NMESSAGES) = 'Openfile failure for ' // Trim(FILNAM)
         ACTIONS(NMESSAGES)  = 'Find the Right File!!'

!  Copy to exception handling

         VLIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
         VLIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
         VLIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
         VLIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS

         RETURN
      ENDIF

!  Read standard inputs

      CALL VLIDORT_READ_INPUTS ( &
         VLIDORT_FixIn, VLIDORT_ModIn,& !InOut
         STATUS_SUB,                  & !Output
         NMESSAGES, MESSAGES, ACTIONS ) !InOut

!  Read file error

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) STATUS = VLIDORT_SERIOUS

!  Copy to exception handling

      VLIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
      VLIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
      VLIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
      VLIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS
      CLOSE(FILUNIT)

!  Return

      RETURN

!  Finish

      END SUBROUTINE VLIDORT_input_master

!

      SUBROUTINE VLIDORT_INIT_INPUTS &
      ( VLIDORT_FixIn, VLIDORT_ModIn ) !Outputs

!  Initialize variables

!     Version 2.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 2.6: New: Observation Geometry variables, 25 October 2012
!     Version 2.8: DO_WATER_LEAVING and DO_FLUORESCENCE flags (added 1/31/16).
!     Version 2.8: Transflux iteration control for the WATER_LEAVING case. 3 "TF" variables added 7/6/16.
!     Version 2.8: FO/SS flags upgraded. 9/19/16, 3/1/17. Argument list rearranged.
!                  VLIDORT input type structures now used for argument passing - 9/19/2017
!                  All VLIDORT fixed and modified inputs now initialized - 9/19/2017

!   Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination
!  --- Also include value of this illumination
   
!  Initialises all inputs for VLIDORT
!  ---------------------------------

!  Module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS, &
                                 MAX_USER_LEVELS, MAXBEAMS, ZERO

!  Implicit none

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT standard input structures

      TYPE(VLIDORT_Fixed_Inputs), INTENT(OUT)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT(OUT) :: VLIDORT_ModIn

!  Initialize inputs
!  =================

!  VLIDORT Fixed Boolean
!  ---------------------

      VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE        = .FALSE.

      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION    = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION    = .FALSE.

      VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL      = .FALSE.

!  Surface control (New, 23 March 2010).  removed for Version 2.8, 7/6/16
      !VLIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE        = .FALSE.

      VLIDORT_FixIn%Bool%TS_DO_UPWELLING           = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_DNWELLING           = .FALSE.

!  Stream angle flag. removed for Version 2.8, 7/6/16

      !VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT         = .FALSE.

!  Contributions (RT Solutions Use Only)

      VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS        = .FALSE.

      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE  = .FALSE.

!  Special Options (RT Solutions Use Only)

      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_1 = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_2 = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_3 = .FALSE.

!  Surface leaving Control. New 17 May 2012

      VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING     = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC        = .FALSE.

!  Water leaving flag added 28 October 2015 (2.7a)
!  Fluorescence  flag added 31 January 2016  (2.8)

      VLIDORT_FixIn%Bool%TS_DO_WATER_LEAVING       = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_FLUORESCENCE        = .FALSE.

!  Water-leaving transmittance iteration flag added 7/6/16 (2.8)

      VLIDORT_FixIn%Bool%TS_DO_TF_ITERATION        = .FALSE.

!  Water-leaving output flag. 3/18/19 for Version 2.8.1
      
      VLIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT   = .FALSE.

!  Planetary problem, and ALBTRN_MEDIA
!   4/26/19 Added control for the media problem. Version 2.8.1

      VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA        = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM   = .FALSE.

!  TOA/BOA Illumination flags. 3/23/19 for Version 2.8.1

      VLIDORT_FixIn%Bool%TS_DO_TOA_ILLUMINATION    = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_BOA_ILLUMINATION    = .FALSE.

!  VLIDORT Modified Boolean
!  ------------------------

!  Flags for FO and SS Corr. Newly re-arranged, Version 3.8, 3/3/17

!  FO (First-Order) choices (Now all modified Booleans)
!  ----------------------------------------------------

!  SINGLE-SCATTER and DIRECT-BOUNCE for Solar Sources
!  DIRECT_PLANCKF and DIRECT-SURFBB for Thermal Sources

!  FO choices completely rewritten, Version 2.8. 9/18/16, 3/1/17, 9/19/17.
!mick mod 9/19/2017 - reordered FO variables in a manner similar to LIDORT3.8

!  Flag for Computing the corrected FO solution, using FO code Version 1.5
!     New 5 Jul 2013, when it was originally named DO_FO_CALC.
!     - if not set, then VLIDORT will perform a truncated pseudo-spherical SS calculation
!     - If not set, then all other SS choices are turned off

!mick mod 9/19/2017 - DO_FOCORR_ALONE now defined internally

      VLIDORT_ModIn%MBool%TS_DO_FOCORR               = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL      = .FALSE.
      !VLIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE         = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR         = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING      = .FALSE.

!  4/15/20. Version 2.8.2. These two flags are removed
!      VLIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION    = .FALSE.
!      VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT       = .FALSE.

!  Criticality flag is new for this version
!  4/15/20. Version 2.8.2. New

      VLIDORT_ModIn%MBool%TS_DO_FOCORR_DOCRIT        = .FALSE.

!  Other flags

      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST      = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES        = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY  = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION     = .FALSE.

!     Isotropic, no_azimuth and All-Fourier flags disabled.

      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY        = .FALSE.
      !VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY       = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH           = .FALSE.
      !VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER          = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING       = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING      = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING      = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES        = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT     = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY           = .FALSE.

      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY    = .FALSE.

!  Observation-Geometry input control. New 25 October 2012

      VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = .FALSE.

!  4/15/20. Doublet-Geometry input control.

      VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY     = .FALSE.

!  Additional Control for Externalized water-leaving inputs. 3/18/19 for Version 2.8.1

      VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE      = .FALSE.

!  VLIDORT Fixed Control
!  ---------------------

!  Taylor ordering parameter. Should be set to 2 or 3
!     Added, 2/19/14 for Taylor-series expansions

      VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER     = 0

      VLIDORT_FixIn%Cont%TS_NSTOKES          = 0
      VLIDORT_FixIn%Cont%TS_NSTREAMS         = 0
      VLIDORT_FixIn%Cont%TS_NLAYERS          = 0
      VLIDORT_FixIn%Cont%TS_NFINELAYERS      = 0
      VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = 0

      VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY = ZERO

!  Special Options (RT Solutions Use Only)

      VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS     = 0
      VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF   = 0

!  Water-leaving: Control for iterative calculation of transmittanaces
!    Variables added for Version 2.8, 7/6/16

      VLIDORT_FixIn%Cont%TS_TF_MAXITER       = 0
      VLIDORT_FixIn%Cont%TS_TF_CRITERION     = ZERO

!  TOA/BOA Illumination. 3/23/19 for Version 2.8.1
!    Must be solar-flux normalized
      
      VLIDORT_FixIn%Cont%TS_TOA_ILLUMINATION = ZERO
      VLIDORT_FixIn%Cont%TS_BOA_ILLUMINATION = ZERO

!  VLIDORT Modified Control
!  ------------------------

!  4/15/20. Version 2.8.2.
!    Two  quantities (SPECHEIGHT, EARTHRADIUS), moved from elsewhere
!    ACRIT is new for this version
      
      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 0
      VLIDORT_ModIn%MCont%TS_GEOMETRY_SPECHEIGHT  = ZERO
      VLIDORT_ModIn%MCont%TS_EARTH_RADIUS         = ZERO
      VLIDORT_ModIn%MCont%TS_DO_FOCORR_ACRIT      = ZERO

!  VLIDORT Fixed Sunrays
!  ---------------------

      VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR   = ZERO

!  VLIDORT Modified Sunrays
!  ------------------------

      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES   = 0
      VLIDORT_ModIn%MSunrays%TS_SZANGLES     = ZERO

!  VLIDORT Fixed UserValues
!  ------------------------

      VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS = 0

!  VLIDORT Modified UserValues
!  ---------------------------

      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS       = 0
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS         = ZERO

      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES      = 0
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT  = ZERO

!  User-defined vertical level output. From Top-of-atmosphere. Example for 4 outputs:
!     USER_LEVELS(1) = 0.0           --> Top-of-atmosphere
!     USER_LEVELS(2) = 1.1           --> One tenth of the way down into Layer 2
!     USER_LEVELS(3) = 3.5           --> One half  of the way down into Layer 4
!     USER_LEVELS(4) = dble(NLAYERS) --> Bottom of atmosphere

      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS          = ZERO

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS      = 0
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT  = ZERO

!  VLIDORT Fixed Chapman
!  ---------------------

      VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID       = ZERO

      VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID     = ZERO
      VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID  = ZERO

      VLIDORT_FixIn%Chapman%TS_FINEGRID          = 0

      VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = ZERO

!  4/15/20. Version 2.8.2, VLIDORT Modified Chapman/Optical, now retired.

!  VLIDORT Fixed Optical
!  ---------------------

      VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT     = ZERO
      VLIDORT_FixIn%Optical%TS_OMEGA_TOTAL_INPUT     = ZERO
      VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT  = ZERO

!  F-matrix input as an alternative to the use of expansion coefficients in
!   the single-scatter (SSCORR and FO) codes.
!  Introduced for Version 2.8 by R. Spurr, 02/08/16
!mick mod 9/19/2017 - changed names from "FMATRIX_INPUT_UP/DN" to "FMATRIX_UP/DN"
!                     for consistency with the companion linearized input variables

      VLIDORT_FixIn%Optical%TS_FMATRIX_UP            = ZERO
      VLIDORT_FixIn%Optical%TS_FMATRIX_DN            = ZERO

      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO     = ZERO

      VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT      = ZERO
      VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT      = ZERO

!  Special LTE variables (RT Solutions Use Only).
!   This has been superseded in Version 2.7. No longer required

      !VLIDORT_FixIn%Optical%TS_LTE_DELTAU_VERT_INPUT = ZERO
      !VLIDORT_FixIn%Optical%TS_LTE_THERMAL_BB_INPUT  = ZERO

!  Rob Fix 3/18/15. Add Wavelength (Microns) as a Diagnostic

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH      = ZERO

!  VLIDORT Fixed Write
!  -------------------

      VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE          = .FALSE.

      VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT          = .FALSE.
      VLIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME    = ' '

      VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO       = .FALSE.
      VLIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME = ' '

      VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER        = .FALSE.
      VLIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME  = ' '

      VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS        = .FALSE.
      VLIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME  = ' '

! Finish

      RETURN
      END SUBROUTINE VLIDORT_INIT_INPUTS

!

      SUBROUTINE VLIDORT_Sup_Init ( VLIDORT_Sup )

      USE VLIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT supplement input structure

      TYPE(VLIDORT_Sup_InOut), INTENT(INOUT) :: VLIDORT_Sup

!  Initialize VLIDORT supplement inputs
!  ====================================

      CALL VLIDORT_BRDF_Sup_Init   ( VLIDORT_Sup )
      CALL VLIDORT_SLEAVE_Sup_Init ( VLIDORT_Sup )
      CALL VLIDORT_SS_Sup_Init     ( VLIDORT_Sup )

!  Finish

      END SUBROUTINE VLIDORT_Sup_Init

!

      SUBROUTINE VLIDORT_BRDF_Sup_Init ( VLIDORT_Sup )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT supplement input structure

      TYPE(VLIDORT_Sup_InOut), INTENT(INOUT) :: VLIDORT_Sup

!  Initialize VLIDORT brdf supplement inputs
!  =========================================

      VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC = ZERO
      VLIDORT_Sup%BRDF%TS_BRDF_F_0        = ZERO
      VLIDORT_Sup%BRDF%TS_BRDF_F          = ZERO
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0   = ZERO
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F     = ZERO

      VLIDORT_Sup%BRDF%TS_EMISSIVITY      = ZERO
      VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY = ZERO

!  Finish

      END SUBROUTINE VLIDORT_BRDF_Sup_Init

!

      SUBROUTINE VLIDORT_SLEAVE_Sup_Init ( VLIDORT_Sup )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT supplement input structure

      TYPE(VLIDORT_Sup_InOut), INTENT(INOUT) :: VLIDORT_Sup

!  Initialize VLIDORT sleave supplement inputs
!  ===========================================

      VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC  = ZERO
      VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES = ZERO
      VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0        = ZERO
      VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0   = ZERO

!  Finish

      END SUBROUTINE VLIDORT_SLEAVE_Sup_Init

!

      SUBROUTINE VLIDORT_SS_Sup_Init ( VLIDORT_Sup )

      USE VLIDORT_PARS_m, Only : ZERO
      USE VLIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  VLIDORT supplement input structure

      TYPE(VLIDORT_Sup_InOut), INTENT(INOUT) :: VLIDORT_Sup

!  Initialize VLIDORT single-scatter supplement inputs
!  ===================================================

      VLIDORT_Sup%SS%TS_STOKES_SS = ZERO
      VLIDORT_Sup%SS%TS_STOKES_DB = ZERO

!  Finish

      END SUBROUTINE VLIDORT_SS_Sup_Init

!

      SUBROUTINE VLIDORT_READ_INPUTS ( &
      VLIDORT_FixIn, VLIDORT_ModIn,    & !InOut
      STATUS,                          & !Output
      NMESSAGES, MESSAGES, ACTIONS )     !InOut

!  Read all control inputs for VLIDORT
!  -----------------------------------

!     Version 2.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 2.8: DO_WATER_LEAVING and DO_FLUORESCENCE flags (added (1/31/16)
!     Version 2.8: Transflux iteration control for the WATER_LEAVING case. 3 "TF" variables added 7/6/16.
!     Version 2.8: FO/SS flags upgraded. 9/19/16, 3/1/17. Argument list rearranged.

!   Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination
!  --- Also include value of this illumination

!  Module, dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAXMOMENTS_INPUT, MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAX_SZANGLES,   &
                                 MAX_MESSAGES, MAX_USER_RELAZMS, MAX_USER_VZANGLES, MAX_USER_LEVELS, &
                                 MAX_USER_OBSGEOMS, MAX_THERMAL_COEFFS, MAXFINELAYERS,               &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, VLIDORT_INUNIT, ONE

      IMPLICIT NONE

!  InOut
!  -----

      TYPE(VLIDORT_Fixed_Inputs), INTENT(INOUT)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT(INOUT) :: VLIDORT_ModIn

!  Outputs
!  -------

!  Exception handling
!     Message Length should be at least 120 Characters

      INTEGER, INTENT(OUT) ::                STATUS
      INTEGER, INTENT(INOUT) ::              NMESSAGES
      CHARACTER (LEN=*), INTENT(INOUT) ::    MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=*), INTENT(INOUT) ::    ACTIONS  ( 0:MAX_MESSAGES )

!  Local variables
!  ---------------

!  Flag for Full Radiance  calculation
!    If not set, just produce diffuse (Multiple scatter) field

      LOGICAL ::   DO_FULLRAD_MODE

!  FO (First-Order) choices (Now all modified Booleans)
!  ----------------------------------------------------

!  SINGLE-SCATTER and DIRECT-BOUNCE for Solar Sources
!  DIRECT_PLANCKF and DIRECT-SURFBB for Thermal Sources

!  FO choices completely rewritten, Version 2.8, 3/1/17

!  Flag for Computing the corrected FO solution, using FO code Version 1.5
!     New 5 Jul 2013, when it was originally named DO_FO_CALC.
!     - if not set, then VLIDORT will perform a truncated pseudo-spherical SS calculation
!     - If not set, then all other SS choices are turned off

      LOGICAL   :: DO_FOCORR

!  Flag for Use of Externally-derived FO results
!     - DO_FOCORR must be set first.

      LOGICAL   :: DO_FOCORR_EXTERNAL

!  Flag for Doing FO calculation alone (no Multiple scatter)
!     - Formerly called DO_SSFULL (confusingly!)
!     - DO_FOCORR must be set first.
!mick mod 9/19/2017 - DO_FOCORR_ALONE now defined in VLIDORT internally

      !LOGICAL   :: DO_FOCORR_ALONE

!   sphericity options.
!     - Solar scattering : FOCORR_NADIR and FOCORR_OUTGOING are mutually exclusive. This is checked.
!     - Solar scattering : If both FOCORR_NADIR and FOCORR_OUTGOING, then PLANE-PARALLEL
!     - Direct Planck    : FOCORR_OUTGOING or PLANE-PARALLEL
!     - DO_FOCORR must be set first.

      LOGICAL   :: DO_FOCORR_NADIR
      LOGICAL   :: DO_FOCORR_OUTGOING

!  Outgoing sphericity - criticality flag and attenuation criterion
!  4/15/20. New for Version 2.8.2 

      LOGICAL   :: DO_FOCORR_DOCRIT
      DOUBLE PRECISION  :: DO_FOCORR_ACRIT

!  Additional SSCORR flags
!  -----------------------

!  Flag for performing SSCORR truncation. Single-scattering
!     - Kept here, but disabled for Version 2.8.
!   4/15/20. Version 2.8.2. Flag removed.
!      LOGICAL   :: DO_SSCORR_TRUNCATION

!  Flag for using Phase Matrices in Single-scatter calculations (instead of GSF Expansion Coefficients)
!     - Introduced for Version 2.8, 3/3/17.  R. Spurr
!     - DO_FOCORR must be set first.
!     - Does not apply to thermal case.
!   4/15/20. Version 2.8.2. Flag retired.
!      LOGICAL   :: DO_SSCORR_USEFMAT

!  Local SS variables (Version 2.7 and earlier)
!      LOGICAL   :: DO_SSCORR_NADIR
!      LOGICAL   :: DO_SSCORR_OUTGOING
!      LOGICAL   :: DO_SSCORR_TRUNCATION
!      LOGICAL   :: DO_SS_EXTERNAL
!      LOGICAL   :: DO_SSFULL 

!  Other flags
!  -----------

!  Basic top-level control

      LOGICAL :: DO_SOLAR_SOURCES
      LOGICAL :: DO_THERMAL_EMISSION

!  directional control

      LOGICAL :: DO_UPWELLING
      LOGICAL :: DO_DNWELLING

!  stream angle flag. Normally required for post-processing solutions
!    ( Exception : DO_MVOUT_ONLY is set, then only want Flux output)

      LOGICAL :: DO_USER_VZANGLES

!  Observation-Geometry input control. 10/25/12

      LOGICAL :: DO_OBSERVATION_GEOMETRY

!  4/15/20. Add DO_DOUBLET_GEOMETRY flag

      LOGICAL :: DO_DOUBLET_GEOMETRY

!  Beam particular solution, plane parallel flag
!    - Not normally required; pseudo-spherical if not set

      LOGICAL :: DO_PLANE_PARALLEL

!  Transmittance only for thermal mode.

      LOGICAL :: DO_THERMAL_TRANSONLY

!  Flag for use of Lambertian surface
!    - If not set, default to BRDF surface

      LOGICAL :: DO_LAMBERTIAN_SURFACE

!  Surface emission flag

      LOGICAL :: DO_SURFACE_EMISSION

!  mean value control (1). If set --> Flux output AS WELL AS Intensities

      LOGICAL :: DO_ADDITIONAL_MVOUT

!  mean value control (2). If set --> only Flux output (No Intensities)
!    - DO_USER_STREAMS should be turned off

      LOGICAL :: DO_MVOUT_ONLY

!  Beam particular solution: Flag for calculating solar beam paths
!    ( Chapman factors = slant/vertical path-length ratios)
!     - This should normally be set. 

      LOGICAL :: DO_CHAPMAN_FUNCTION

!  Beam particular solution: Flag for using refraction in solar paths
!     - This should NOT normally be set. 

      LOGICAL :: DO_REFRACTIVE_GEOMETRY

!  Flag for Use of Delta-M scaling
!    - Should normally be set
!    - Not required for DO_RAYLEIGH_ONLY or DO_ISOTROPIC_ONLY

      LOGICAL :: DO_DELTAM_SCALING

!  double convergence test flag

      LOGICAL :: DO_DOUBLE_CONVTEST

!  Performance flags
!    -- SOLUTION_SAVING gets rid of unneeded RTE computations
!    -- BVP_TELESCOPING creates reduced Boundary value problems
!    -- These flags should be used with CAUTION
!    -- Best, Rayleigh atmospheres with few contiguous cloud/aerosol layers

      LOGICAL :: DO_SOLUTION_SAVING
      LOGICAL :: DO_BVP_TELESCOPING

!  scatterers and phase function control
!    - Rayleigh only, if set, make sure that scattering Law is Rayleigh!

      LOGICAL :: DO_RAYLEIGH_ONLY

!  Surface leaving control
!     Version 2.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 2.8: DO_WATER_LEAVING (added 10/28/15) and DO_FLUORESCENCE (added 1/31/16) flags.

      LOGICAL :: DO_SURFACE_LEAVING
      LOGICAL :: DO_SL_ISOTROPIC

!  Version 2.8: DO_WATER_LEAVING (added 10/28/15) and DO_FLUORESCENCE (added 1/31/16) flags.

      LOGICAL :: DO_WATER_LEAVING
      LOGICAL :: DO_FLUORESCENCE

!  Version 2.8: Transflux iteration control for the WATER_LEAVING case. 3 "TF" variables added 7/6/16.

      LOGICAL  ::  DO_TF_ITERATION
      INTEGER  ::  TF_MAXITER

      DOUBLE PRECISION  ::   TF_CRITERION
      
!  Water-leaving output flag. 3/18/19 for Version 2.8.1

      LOGICAL  ::  DO_WLADJUSTED_OUTPUT

!  4/29/19. Planetary problem and media properties flags,  Version 2.8.1

      LOGICAL  ::  DO_ALBTRN_MEDIA(2)
      LOGICAL  ::  DO_PLANETARY_PROBLEM

!  TOA/BOA illumination control, New Version 2.8.1, 3/23/19

      LOGICAL  ::  DO_TOA_ILLUMINATION
      LOGICAL  ::  DO_BOA_ILLUMINATION

      DOUBLE PRECISION :: TOA_ILLUMINATION
      DOUBLE PRECISION :: BOA_ILLUMINATION

!  Additional Control for Externalized water-leaving inputs. 3/18/19 for Version 2.8.1

      LOGICAL  ::  DO_EXTERNAL_WLEAVE

!  Following removed from Version 2.8, 7/6/16
!      LOGICAL ::             DO_QUAD_OUTPUT

!  debug output control

      LOGICAL ::   DO_DEBUG_WRITE
      LOGICAL ::   DO_WRITE_INPUT
      LOGICAL ::   DO_WRITE_SCENARIO
      LOGICAL ::   DO_WRITE_FOURIER
      LOGICAL ::   DO_WRITE_RESULTS

      CHARACTER (LEN=60) ::  INPUT_WRITE_FILENAME
      CHARACTER (LEN=60) ::  SCENARIO_WRITE_FILENAME
      CHARACTER (LEN=60) ::  FOURIER_WRITE_FILENAME
      CHARACTER (LEN=60) ::  RESULTS_WRITE_FILENAME

!  2. CONTROL INTEGERS
!  ===================

!  Order of Taylor series (including terms up to EPS^n). Introduced 10/10/13 for Version 3.7
      
      INTEGER :: TAYLOR_ORDER

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  number of computational layers

      INTEGER :: NLAYERS

!  number of Stokes vector elements

      INTEGER :: NSTOKES

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing single scattering correction )

      INTEGER :: NFINELAYERS

!  number of GSF F-matrix expansion moments (External to the model)

      INTEGER :: NGREEK_MOMENTS_INPUT

!  number of solar beams to be processed

      INTEGER :: N_SZANGLES

!  Number of user-defined relative azimuths

      INTEGER :: N_USER_RELAZMS

!  Number of User-defined viewing zenith angles (0 to 90 degrees)

      INTEGER :: N_USER_VZANGLES

!  Number of User-defined vertical levels for  output

      INTEGER :: N_USER_LEVELS

!  Number of thermal coefficients (2 should be the default)

      INTEGER :: N_THERMAL_COEFFS

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      INTEGER :: N_USER_OBSGEOMS

!  3. CONTROL NUMBERS
!  ==================

!  Flux factor ( should be 1 or pi ). Same for all beams.

      DOUBLE PRECISION :: FLUX_FACTOR

!  accuracy for convergence of Fourier series

      DOUBLE PRECISION :: VLIDORT_ACCURACY

!  Zenith tolerance (nearness of output zenith cosine to 1.0 )
!    removed 02 June 2010
!      DOUBLE PRECISION :: ZENITH_TOLERANCE

!  Atmospheric wavelength, new Version 2.8, bookkeeping

      DOUBLE PRECISION :: ATMOS_WAVELENGTH

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT

      DOUBLE PRECISION :: EARTH_RADIUS

!  Refractive index parameter
!  ( Only required for refractive geometry attenuation of the solar beam)

      DOUBLE PRECISION :: RFINDEX_PARAMETER

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007
!    -- See special note below

      DOUBLE PRECISION :: GEOMETRY_SPECHEIGHT

!  Lambertian Albedo

      DOUBLE PRECISION :: LAMBERTIAN_ALBEDO

!  BOA solar zenith angles (degrees)

      DOUBLE PRECISION :: SZANGLES ( MAX_SZANGLES )

!  user-defined relative azimuths (degrees) (mandatory for Fourier > 0)

      DOUBLE PRECISION :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  User-defined viewing zenith angles input (degrees) 

      DOUBLE PRECISION :: USER_VZANGLES ( MAX_USER_VZANGLES )

!  User-defined vertical levels for output
!    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
!    E.g. For 4.5, this means half way down the 5th layer
!    E.g. For 0.0, this is output at TOA

      DOUBLE PRECISION :: USER_LEVELS  (MAX_USER_LEVELS)

!  User-defined Observation Geometry angle input
!   New variable, 25 OCtober 2012, for Observational Geometry input

      DOUBLE PRECISION :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Other variables
!  ===============

      CHARACTER (LEN=9), PARAMETER :: PREFIX = 'VLIDORT -'
      LOGICAL            :: ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER            :: I, FILUNIT, NM

!  Initialize status

      STATUS = VLIDORT_SUCCESS
      ERROR  = .FALSE.
      NM     = 0

!  These are already initialized in calling routine
!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
!      ACTIONS(0)      = 'No Action required for this Task'

!  File unit

      FILUNIT = VLIDORT_INUNIT

!  1. READ ALL CONTROL VARIABLES (BOOLEAN INPUTS)
!  ==============================================

!  Operation modes
!  ---------------

!  Full Stokes vector calculation, SS + MS fields
!    if False, then VLIDORT only does a multiple scatter calculation.

      PAR_STR = 'Do full Stokes vector calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FULLRAD_MODE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  FO (First-Order) choices (Now all modified Booleans)
!  ----------------------------------------------------

!  SINGLE-SCATTER and DIRECT-BOUNCE for Solar Sources
!  DIRECT_PLANCKF and DIRECT-SURFBB for Thermal Sources

!  FO choices completely rewritten, Version 2.8. 9/18/16, 3/1/17.

!  Flag for Computing the corrected FO solution, using FO code Version 1.5
!     New 5 Jul 2013, when it was originally named DO_FO_CALC.
!     - if not set, then VLIDORT will perform a truncated pseudo-spherical SS calculation
!     - If not set, then all other SS choices are turned off

      PAR_STR = 'Do First-Order (FO) correction?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All other options depend on this flag.

      IF ( DO_FOCORR ) THEN

!  Flag for Use of Externally-derived FO results
!     - DO_FOCORR must be set first. New 15 March 2012.

        PAR_STR = 'Do external First-Order correction?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_EXTERNAL
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Flag for Doing FO calculation alone (no Multiple scatter)
!     - Formerly called DO_SSFULL (confusingly!)
!     - DO_FOCORR must be set first.
!mick mod 9/19/2017 - Turned off.  Now defined in VLIDORT internally

        !PAR_STR = 'Do First-Order correction alone?'
        !IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_ALONE
        !CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Control for executing the FO code Version 1.5.
!    -  Only if the external field is turned off.

        IF ( .not. DO_FOCORR_EXTERNAL ) THEN

!   sphericity options.
!     - Solar scattering : FOCORR_NADIR and FOCORR_OUTGOING are mutually exclusive. This is checked.
!     - Solar scattering : If both FOCORR_NADIR and FOCORR_OUTGOING, then PLANE-PARALLEL
!     - Direct Planck    : FOCORR_OUTGOING or PLANE-PARALLEL
!     - DO_FOCORR must be set first. DO_FOCORR_EXTERNAL should be off.

          PAR_STR = 'Do nadir First-Order correction?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_NADIR
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'Do outgoing First-Order correction?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_OUTGOING
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  outgoing sphericity --> Flag for criticality, and criticality attenuation.
!  4/15/20. Version 2.8.2

          IF ( DO_FOCORR_OUTGOING ) THEN
             PAR_STR = 'Do Criticality for First-Order correction?'
             IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_DOCRIT
             CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

             PAR_STR = 'Criticality Attenuation for First-Order correction'
             IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_ACRIT
             CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          ENDIF

!  Flag for using the F-matrix in the SS calculations (instead of Greekmat Coefficients)
!     - Introduced for Version 2.8, 7/7/16.  R. Spurr
!   4/15/20. Version 2.8.2. Flag retired.
!          PAR_STR = 'Do Fmatrix usage in single scatter correction?'
!          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SSCORR_USEFMAT
!          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Additional truncation scaling for single-scatter corrections
!   --- Disabled for Version 2.8, do we require it again ???
!        PAR_STR = 'Do truncation scaling on single scatter corrections?'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SSCORR_TRUNCATION
!        CALL FINDPAR_ERROR  ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  End First-Order options

        ENDIF

      ENDIF

!  Direct beam correction (BRDF options only)
!      PAR_STR = 'Do direct beam correction?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DBCORRECTION
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Multiple scatter source function output control. Removed. 30 March 2007
!      PAR_STR = 'Output multiple scatter layer source functions?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) SAVE_LAYER_MSST
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  4/29/19 for Version 2.8.1. AlbTrn Media control
      
      PAR_STR = 'Do Media properties calculation with TOA illumination?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ALBTRN_MEDIA(1)
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      PAR_STR = 'Do Media properties calculation with BOA illumination?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ALBTRN_MEDIA(2)
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!   4/29/19 for Version 2.8.1. Planetary problem input

      PAR_STR = 'Do planetary problem calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_PLANETARY_PROBLEM
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Solar beam control
!  ------------------

!  Basic control

      PAR_STR = 'Use solar sources?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SOLAR_SOURCES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Other control options for the solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

!  Pseudo-spherical control

        PAR_STR = 'Do plane-parallel treatment of direct beam?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_PLANE_PARALLEL
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Internal Chapman function calculation

        PAR_STR = 'Do internal Chapman function calculation?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_CHAPMAN_FUNCTION
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Refractive atmosphere

        PAR_STR = 'Do refractive geometry?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_REFRACTIVE_GEOMETRY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  End control

      ENDIF

!  Direct-beam control
!      PAR_STR = 'Include direct beam?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DIRECT_BEAM
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!  Solution method (Classical only). Variable is set in DERIVE_INPUTS.
!      DO_CLASSICAL_SOLUTION = .TRUE.
!      PAR_STR = 'Use classical beam solution?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_CLASSICAL_SOLUTION
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Thermal Control
!  ---------------

!  Atmospheric thermal emission, Basic control

      PAR_STR = 'Do thermal emission?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_THERMAL_EMISSION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_THERMAL_EMISSION ) THEN

!  Thermal sources, transmittance only

        PAR_STR = 'Do thermal emission, transmittance only?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_THERMAL_TRANSONLY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Number of coefficients (includes a dimensioning check)

        PAR_STR = 'Number of thermal coefficients'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_THERMAL_COEFFS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_THERMAL_COEFFS .GT. MAX_THERMAL_COEFFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = &
             'Entry under "Number of thermal coefficients" >'// &
             ' allowed Maximum dimension'
            ACTIONS(NM)  = &
             'Re-set input value or increase MAX_THERMAL_COEFFS dimension '// &
             'in VLIDORT_PARS'
          STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
        ENDIF
      ENDIF

!  Surface emission control

      PAR_STR = 'Do surface emission?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SURFACE_EMISSION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Surface control
!  ---------------

!  Lambertian surface

      PAR_STR = 'Do Lambertian surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_LAMBERTIAN_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!    START Surface Leaving control (New 17 May 2012)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

!  Basic control

      PAR_STR = 'Do surface-leaving term?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SURFACE_LEAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_SURFACE_LEAVING ) THEN
         
!  Additional Control for Externalized water-leaving inputs. Introduced 3/18/19 for Version 2.8.1
!    -- This is a general flag, you can still adjust the water-leaving (if flagged)         
         
         PAR_STR = 'Do external water-leaving production?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_EXTERNAL_WLEAVE
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Isotropic control

         PAR_STR = 'Do isotropic surface-leaving term?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SL_ISOTROPIC
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Water-leaving control. New 28 October 2015

         PAR_STR = 'Do Water-leaving option?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WATER_LEAVING
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  This section added 7/6/16 for Version 2.8.
!   The 3 TF inputs will control the water-leaving Transmittance calculation 
!      Also included now is the Water-leaving output flag. 3/18/19 for Version 2.8.1

         IF ( DO_WATER_LEAVING ) THEN
           PAR_STR = 'Do iterative calculation of Water-leaving transmittance?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_TF_ITERATION
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
 
           PAR_STR = 'Flag for output of transmittance-adjusted water-leaving radiances'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WLADJUSTED_OUTPUT
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

           IF ( DO_TF_ITERATION ) THEN

             PAR_STR = 'Maximum number of iterations in calculation of Water-leaving transmittance'
             IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) TF_MAXITER
             CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

             PAR_STR = 'Convergence criterion for iterative calculation of Water-leaving transmittance'
             IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) TF_CRITERION
             CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

           ENDIF
         ENDIF

!  Fluorescence control. New 31 January 2016
!   ( This is mutually exclusive from Water-leaving; condition will be checked )

         PAR_STR = 'Do Fluorescence option?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FLUORESCENCE
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!     END Surface Leaving control
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!  TOA/BOA Illumination. 3/23/19 for Version 2.8.1  NEW SECTION
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

      PAR_STR = 'Do TOA Illumination for Airglow?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_TOA_ILLUMINATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_TOA_ILLUMINATION ) THEN
         PAR_STR = ' TOA Illumination Flux (sun-normalized)'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) TOA_ILLUMINATION
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

      PAR_STR = 'Do BOA Illumination for nighttime?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_BOA_ILLUMINATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_BOA_ILLUMINATION ) THEN
         PAR_STR = ' BOA Illumination Flux (sun-normalized)'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) BOA_ILLUMINATION
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Performance control
!  -------------------

!  Delta-M scaling. Now set regardless of source
!    Surely a mistake (Versions 2.7 and earlier) --> Should only be set for solar beam sources

!      IF ( DO_SOLAR_SOURCES ) THEN
      PAR_STR = 'Do delta-M scaling?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DELTAM_SCALING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      ENDIF

!  Double convergence test

      PAR_STR = 'Do double convergence test?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DOUBLE_CONVTEST
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Solution saving mode.
!    New code, RJDS, RT Solutions, Inc. 4/11/05.

      PAR_STR = 'Do solution saving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SOLUTION_SAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Boundary value problem (BVP) telescope mode.
!    New code, RJDS, RT Solutions, Inc. 4/11/05.

      PAR_STR = 'Do boundary-value telescoping?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_BVP_TELESCOPING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  New flag, Observational Geometry
!    Added, 25 October 12
!mick mod 9/19/2017 - added ELSE section
!  4/15/20. Add DO_DOUBLET_GEOMETRY flag. Turn off if no solar sources

      IF ( DO_SOLAR_SOURCES ) THEN
         PAR_STR = 'Do Observation Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_OBSERVATION_GEOMETRY
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            DO_DOUBLET_GEOMETRY = .false.
         ELSE
            PAR_STR = 'Do Doublet Geometry?'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DOUBLET_GEOMETRY
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
         ENDIF 
      ELSE
         DO_OBSERVATION_GEOMETRY = .FALSE.
         DO_DOUBLET_GEOMETRY     = .FALSE.
      ENDIF

!  User-defined output control
!  ---------------------------

!  Directional output control

      PAR_STR = 'Do upwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_UPWELLING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do downwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DNWELLING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Stream angle and optical depth output control
!  ---------------------------------------------

!  User-defined viewing zenith angle

!  Removed.
!      PAR_STR = 'Include quadrature angles in output?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR) READ (FILUNIT,*,ERR=998) DO_QUAD_OUTPUT
!      CALL FINDPAR_ERROR( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_USER_VZANGLES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Old code, version 2.3 and earlier............
!      PAR_STR = 'User-defined optical depths?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_USER_TAUS
!      CALL FINDPAR_ERROR( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'Layer boundary optical depths?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_LBOUND_TAUS
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Mean-value output control
!  -------------------------

      PAR_STR = 'Do mean-value output additionally?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ADDITIONAL_MVOUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do only mean-value output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_MVOUT_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Numerical control (azimuth series)
!  ----------------------------------

!  Scatterers and phase function control

      PAR_STR='Do Rayleigh atmosphere only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_RAYLEIGH_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Removed isotropic-only option, 17 January 2006.
!      PAR_STR='Isotropic atmosphere only?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ISOTROPIC_ONLY
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  No azimuth dependence (TRUE means Fourier m = 0 only )
!    Removed. Now this is a bookkeeping variable.
!    17 January 2006.
!      PAR_STR = 'No azimuth dependence in the solution?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_NO_AZIMUTH
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All possible Fourier components (2N-1). Debug only
!      PAR_STR = 'Compute all Fourier components?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ALL_FOURIER
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Write control
!  -------------

!  Output write flags

      PAR_STR = 'Do debug write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DEBUG_WRITE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do input control write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WRITE_INPUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do input scenario write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WRITE_SCENARIO
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do Fourier component output write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WRITE_FOURIER
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'DO results write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WRITE_RESULTS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Output filenames

      PAR_STR = 'filename for input write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,'(a)',ERR=998) INPUT_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'filename for scenario write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,'(a)',ERR=998) SCENARIO_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'filename for Fourier output write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,'(a)',ERR=998) FOURIER_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'filename for main output'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,'(a)',ERR=998) RESULTS_WRITE_FILENAME
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  2. READ ALL THE CONTROL NUMBERS (INTEGER INPUTS)
!  ================================================

!  Streams/layers/finelayers/moments (INTEGER input)
!  -------------------------------------------------

!  Taylor order parameter added 2/19/14, Version 2p7

      PAR_STR = 'Number of small-number terms in Taylor series expansions'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) TAYLOR_ORDER
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Number of Stokes parameters

      PAR_STR = 'Number of Stokes vector components'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) NSTOKES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( NSTOKES .GT. MAXSTOKES ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under' &
              //' "Number of Stokes parameters" >'//' allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value to 4 or less'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Number of computational streams

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
               'Entry under "Number of half-space streams" >'// ' allowed Maximum dimension'
        ACTIONS(NM)  = &
             'Re-set input value or increase MAXSTREAMS dimension '// 'in VLIDORT_PARS'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Number of atmospheric layers

      PAR_STR = 'Number of atmospheric layers'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) NLAYERS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
               'Entry under "Number of atmospheric layers" >'// ' allowed Maximum dimension'
        ACTIONS(NM)  = &
             'Re-set input value or increase MAXLAYERS dimension '// 'in VLIDORT_PARS'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Change for Version 2.8. 3/1/17
!mick fix 9/19/2017 - now added DO_FOCORR & DO_FOCORR_EXTERNAL if conditions
!                     to DO_FOCORR_OUTGOING if condition

      !IF ( DO_SSCORR_OUTGOING .OR. DO_SSFULL ) THEN
      IF ( DO_FOCORR ) THEN
        IF ( .NOT.DO_FOCORR_EXTERNAL ) THEN
          IF ( DO_FOCORR_OUTGOING ) THEN

!  Number of fine layers

            PAR_STR = 'Number of fine layers (outgoing sphericity option only)'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) NFINELAYERS
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

            IF ( NFINELAYERS .GT. MAXFINELAYERS ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Entry under "Number of fine layers..." >'//' allowed Maximum dimension'
              ACTIONS(NM)  = 'Re-set input value or increase MAXFINELAYERS dimension in VLIDORT_PARS'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              RETURN
            ENDIF

!  Input Geometry specification height

            PAR_STR = 'Input geometry specification height (km)'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
               READ (FILUNIT,*,ERR=998) GEOMETRY_SPECHEIGHT
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          ENDIF
        ENDIF
      ENDIF

!  Number of scattering matrix expansion coefficients

      PAR_STR = 'Number of scattering matrix expansion coefficients'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) NGREEK_MOMENTS_INPUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( NGREEK_MOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under "Number of input Scattering Matrix expansion coefficients" >'//' allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXMOMENTS_INPUT dimension '//'in VLIDORT_PARS'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Geometry inputs
!  ---------------

!  Observational Geometry control.  New, 25 October 2012

      IF ( DO_OBSERVATION_GEOMETRY ) THEN

!  Number of Observational Geometry inputs

        PAR_STR = 'Number of Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_USER_OBSGEOMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of Observation Geometry inputs" > allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_OBSGEOMS dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
          NMESSAGES    = NM
          RETURN
        ENDIF

!  Automatic setting of N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS and DO_USER_VZANGLES

        N_SZANGLES       = N_USER_OBSGEOMS
        N_USER_VZANGLES  = N_USER_OBSGEOMS
        N_USER_RELAZMS   = N_USER_OBSGEOMS
        DO_USER_VZANGLES = .true.

!  Observational Geometry inputs

        PAR_STR = 'Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           DO I = 1, N_USER_OBSGEOMS
             READ (FILUNIT,*,ERR=998) USER_OBSGEOMS(I,1), USER_OBSGEOMS(I,2), USER_OBSGEOMS(I,3)
           ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!   Automatic setting of SZANGLES, USER_VZANGLES, and USER_RELAZMS

        SZANGLES     (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
        USER_VZANGLES(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        USER_RELAZMS (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)

!  Doublet Geometry Control
!  4/15/20. Add DO_DOUBLET_GEOMETRY inputs

      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN

!  1. Number of Solar zenith angles
!     ---- check not exceeding dimensioned number

        PAR_STR = 'Number of solar zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_SZANGLES
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_SZANGLES .GT. MAX_SZANGLES ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of solar zenith angles" >'//' allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_SZANGLES dimension '//'in VLIDORT_PARS'
          STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
        ENDIF

!  BOA solar zenith angle inputs

        PAR_STR = 'Solar zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_SZANGLES
            READ (FILUNIT,*,ERR=998) SZANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  3. User defined viewing zenith angles (should be positive)
!     ---- check not exceeding dimensioned number

        IF ( DO_USER_VZANGLES ) THEN
          PAR_STR = 'Number of user-defined viewing zenith angles'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_VZANGLES
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          IF ( N_USER_VZANGLES .GT. MAX_USER_VZANGLES ) THEN
            NM = NM + 1
            MESSAGES(NM) = &
               'Entry under "Number of .....viewing zenith angles" >'// &
               ' allowed Maximum dimension'
            ACTIONS(NM)  = &
               'Re-set input value or increase MAX_USER_VZANGLES dimension '// &
               'in VLIDORT_PARS'
            STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
          ENDIF

          PAR_STR = 'User-defined viewing zenith angles (degrees)'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            DO I = 1, N_USER_VZANGLES
              READ (FILUNIT,*,ERR=998) USER_VZANGLES(I)
            ENDDO
          ENDIF
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        ENDIF

!  3. Number of azimuths. MUST EQUAL N_USER_VZANGLES
!     ---- check not exceeding dimensioned number

        N_USER_RELAZMS = N_USER_VZANGLES
        IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of ..... azimuth angles" >'//' allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension '//'in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS ; NMESSAGES    = NM ; RETURN
        ENDIF

!  Azimuth angles

        PAR_STR = 'User-defined relative azimuth angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_RELAZMS
            READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Lattice Geometry control
!  4/15/20. Add DO_DOUBLET_GEOMETRY flag to this option

      ELSE IF ( .not. DO_OBSERVATION_GEOMETRY .and..not. DO_DOUBLET_GEOMETRY ) THEN

!  1. Number of Solar zenith angles
!     ---- check not exceeding dimensioned number

        PAR_STR = 'Number of solar zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_SZANGLES
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_SZANGLES .GT. MAX_SZANGLES ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of solar zenith angles" >'//' allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_SZANGLES dimension '//'in VLIDORT_PARS'
          STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
        ENDIF

!  BOA solar zenith angle inputs

        PAR_STR = 'Solar zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_SZANGLES
            READ (FILUNIT,*,ERR=998) SZANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  2. Number of azimuths
!     ---- check not exceeding dimensioned number

        PAR_STR = 'Number of user-defined relative azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of ..... azimuth angles" >'//' allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension '//'in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS ; NMESSAGES    = NM ; RETURN
        ENDIF

!  Azimuth angles

        PAR_STR = 'User-defined relative azimuth angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_RELAZMS
            READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  3. User defined viewing zenith angles (should be positive)
!     ---- check not exceeding dimensioned number

        IF ( DO_USER_VZANGLES ) THEN
          PAR_STR = 'Number of user-defined viewing zenith angles'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_VZANGLES
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          IF ( N_USER_VZANGLES .GT. MAX_USER_VZANGLES ) THEN
            NM = NM + 1
            MESSAGES(NM) = &
               'Entry under "Number of .....viewing zenith angles" >'// &
               ' allowed Maximum dimension'
            ACTIONS(NM)  = &
               'Re-set input value or increase MAX_USER_VZANGLES dimension '// &
               'in VLIDORT_PARS'
            STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
          ENDIF

          PAR_STR = 'User-defined viewing zenith angles (degrees)'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            DO I = 1, N_USER_VZANGLES
              READ (FILUNIT,*,ERR=998) USER_VZANGLES(I)
            ENDDO
          ENDIF
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        ENDIF

!  End lattice geometry clause

      ENDIF

!  Continuation point for skipping lattice geometry input-reads. 
!5665 CONTINUE   !!! Goto statement removed, version 2.8, replaced by If-block.

!  4. Number of output levels
!     ---- check not exceeding dimensioned number

!  This is designed to ensure that output is not related to optical
!  depth (which is wavelength-dependent); we use a height-based system.

!  User defined boundary (whole layer) and off-boundary (partial layer)
!  output choices are specified as follows.
!       USER_LEVELS(1) = 0.0    Top of the first layer
!       USER_LEVELS(2) = 1.0    Bottom of the first layer
!       USER_LEVELS(3) = 17.49  Output is in Layer 18, at a distance of
!                               0.49 of the way down from top (in height

      PAR_STR = 'Number of user-defined output levels'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) N_USER_LEVELS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( N_USER_LEVELS .GT. MAX_USER_LEVELS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
             'Entry under "Number of ..... output levels" > allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_LEVELS dimension in VLIDORT_PARS'
        STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
      ENDIF

!  Vertical output levels

      PAR_STR = 'User-defined output levels'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_LEVELS
          READ (FILUNIT,*,ERR=998) USER_LEVELS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  3. READ ALL THE FLOATING-POINT INPUTS
!  =====================================

!  Flux constant. Should be set to 1 if no solar sources.
!   Formerly a Book-keeping variable set to 1 in "derive inputs"
!   Now (July 2009) allowed to vary because of thermal emission
!   Must take physical values if using solar + thermal.

      IF ( DO_SOLAR_SOURCES ) THEN
        PAR_STR = 'Solar flux constant'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) FLUX_FACTOR
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ELSE
        FLUX_FACTOR = ONE
      ENDIF

!  Note the following possibility (not yet allowed for)
!      PAR_STR = 'TOA flux vector'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) (FLUXVEC(I),I=1,MAXSTOKES)
!      CALL FINDPAR_ERROR( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Accuracy criterion

      PAR_STR = 'Fourier series convergence'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) VLIDORT_ACCURACY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Zenith tolerance level. Now removed. 17 January 2006.

!  Atmospheric Wavelength [Microns].
!  Configuration file input added for Version 2.8, 1/31/16.
!  Just a diagnostic, only used for comparison with BRDF/SLEAVE wavelengths

      PAR_STR = 'Atmospheric Wavelength [Microns]'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) ATMOS_WAVELENGTH
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Earth radius and RF parameter
!  ---- only for Chapman function calculation
!mick mod 9/19/2017 - added DO_SOLAR_SOURCES if condition

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_CHAPMAN_FUNCTION ) THEN
          PAR_STR = 'Earth radius (km)'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
               READ (FILUNIT,*,ERR=998) EARTH_RADIUS
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          IF ( DO_REFRACTIVE_GEOMETRY ) THEN
            PAR_STR = 'Refractive index parameter'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
                 READ (FILUNIT,*,ERR=998) RFINDEX_PARAMETER
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          ENDIF
        ENDIF
      ENDIF

!  Lambertian surface

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        PAR_STR = 'Lambertian albedo'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) LAMBERTIAN_ALBEDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Define VLIDORT std inputs
!  =========================

!mick mod 9/19/2017 - initializing of all std VLIDORT type structure input variables is done
!                     in subroutine VLIDORT_INIT_INPUTS; only those read in from the VLIDORT 
!                     config file are possibly modified here.  IF conditions are applied where
!                     appropriate.
      
!  TOA/BOA Illumination. 3/23/19 for Version 2.8.1

!  VLIDORT Fixed Boolean Inputs

      VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE       = DO_FULLRAD_MODE
      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION   = DO_THERMAL_EMISSION
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION   = DO_SURFACE_EMISSION
      IF ( DO_SOLAR_SOURCES ) &
         VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL  = DO_PLANE_PARALLEL 
      VLIDORT_FixIn%Bool%TS_DO_UPWELLING          = DO_UPWELLING
      VLIDORT_FixIn%Bool%TS_DO_DNWELLING          = DO_DNWELLING
      !VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT        = DO_QUAD_OUTPUT
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = DO_LAMBERTIAN_SURFACE
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING    = DO_SURFACE_LEAVING
      
      IF ( DO_SURFACE_LEAVING ) THEN
         VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC    = DO_SL_ISOTROPIC
         VLIDORT_FixIn%Bool%TS_DO_WATER_LEAVING   = DO_WATER_LEAVING
         IF ( DO_WATER_LEAVING ) THEN
            VLIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT = DO_WLADJUSTED_OUTPUT   ! New 3/18/19 Version 2.8.1
            VLIDORT_FixIn%Bool%TS_DO_TF_ITERATION = DO_TF_ITERATION
            IF ( DO_TF_ITERATION ) THEN
               VLIDORT_FixIn%Cont%TS_TF_MAXITER   = TF_MAXITER
               VLIDORT_FixIn%Cont%TS_TF_CRITERION = TF_CRITERION
            ENDIF
         ENDIF
         VLIDORT_FixIn%Bool%TS_DO_FLUORESCENCE    = DO_FLUORESCENCE
      ENDIF

!  TOA/BOA Illumination. New code for Version 2.8.1, 3/18/19 
      
      VLIDORT_FixIn%Bool%TS_DO_TOA_ILLUMINATION   = DO_TOA_ILLUMINATION
      if ( DO_TOA_ILLUMINATION ) THEN
         VLIDORT_FixIn%Cont%TS_TOA_ILLUMINATION   = TOA_ILLUMINATION
      endif
      VLIDORT_FixIn%Bool%TS_DO_BOA_ILLUMINATION   = DO_BOA_ILLUMINATION
      if ( DO_BOA_ILLUMINATION ) THEN
         VLIDORT_FixIn%Cont%TS_BOA_ILLUMINATION   = BOA_ILLUMINATION
      endif

!  Flag for the Planetary problem calculation. Version 2.8.1, 4/29/19 

      VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM   = DO_PLANETARY_PROBLEM

!  Flags for the media properties calculation. Version 2.8.1, 4/29/19 

      VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA   = DO_ALBTRN_MEDIA

!  VLIDORT Modified Boolean Inputs
!mick mod 9/19/2017 - DO_FOCORR_ALONE now defined internally
!mick fix 6/4/2019 - added "DO_SURFACE_LEAVING" IF block for DO_EXTERNAL_WLEAVE
!  4/15/20. Version 2.8.2, USEFMAT gone, SSCORR_TRUNCATION gone
!  4/15/20. Version 2.8.2, EARTH_RADIUS, GEOMETRY_SPECHEIGHT now assigned different Types.
!  4/15/20. Version 2.8.2, DOCRIT and ACRIT are new

      VLIDORT_ModIn%MBool%TS_DO_FOCORR                   = DO_FOCORR
      IF ( DO_FOCORR ) THEN
         VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL       = DO_FOCORR_EXTERNAL
         !VLIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE          = DO_FOCORR_ALONE
         IF ( .NOT.DO_FOCORR_EXTERNAL ) THEN
            VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR       = DO_FOCORR_NADIR 
            VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING    = DO_FOCORR_OUTGOING
            IF ( DO_FOCORR_OUTGOING ) THEN
              VLIDORT_FixIn%Cont%TS_NFINELAYERS          = NFINELAYERS
              VLIDORT_ModIn%MCont%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT
              VLIDORT_ModIn%MBool%TS_DO_FOCORR_DOCRIT    = DO_FOCORR_DOCRIT
              VLIDORT_ModIn%MCont%TS_DO_FOCORR_ACRIT     = DO_FOCORR_ACRIT
            ENDIF
            !VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT     = DO_SSCORR_USEFMAT
            !VLIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION  = DO_SSCORR_TRUNCATION
         ENDIF
      ENDIF
      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST          = DO_DOUBLE_CONVTEST
      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES            = DO_SOLAR_SOURCES

!  4/15/20. Add DO_DOUBLET_GEOMETRY flag

      IF ( DO_SOLAR_SOURCES ) THEN
         VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION      = DO_CHAPMAN_FUNCTION
         IF ( DO_CHAPMAN_FUNCTION ) &
            VLIDORT_ModIn%MCont%TS_EARTH_RADIUS       = EARTH_RADIUS

         VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY   = DO_REFRACTIVE_GEOMETRY
         IF ( DO_CHAPMAN_FUNCTION .AND. DO_REFRACTIVE_GEOMETRY ) &
            VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = RFINDEX_PARAMETER

         VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY  = DO_OBSERVATION_GEOMETRY
         VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY      = DO_DOUBLET_GEOMETRY
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS    = N_USER_OBSGEOMS
            VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,1:3) = &
                                     USER_OBSGEOMS(1:N_USER_OBSGEOMS,1:3)
         ENDIF
      ENDIF
      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY            = DO_RAYLEIGH_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY           = DO_ISOTROPIC_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH               = DO_NO_AZIMUTH
      !VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER              = DO_ALL_FOURIER
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING           = DO_DELTAM_SCALING
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING          = DO_SOLUTION_SAVING
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING          = DO_BVP_TELESCOPING
      VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES            = DO_USER_VZANGLES
      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT         = DO_ADDITIONAL_MVOUT
      VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY               = DO_MVOUT_ONLY
      IF ( DO_THERMAL_EMISSION ) &
         VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY     = DO_THERMAL_TRANSONLY

      ! New 4/22/19 Version 2.8.1.
      IF ( DO_SURFACE_LEAVING ) THEN
         VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE       = DO_EXTERNAL_WLEAVE
      ELSE
         VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE       = .false.  
      ENDIF

!  VLIDORT Fixed Control Inputs

      VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER    = TAYLOR_ORDER
      VLIDORT_FixIn%Cont%TS_NSTOKES         = NSTOKES
      VLIDORT_FixIn%Cont%TS_NSTREAMS        = NSTREAMS
      VLIDORT_FixIn%Cont%TS_NLAYERS         = NLAYERS
      IF ( DO_THERMAL_EMISSION ) &
         VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = N_THERMAL_COEFFS
      VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY = VLIDORT_ACCURACY

!  VLIDORT Modified Control Inputs

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = NGREEK_MOMENTS_INPUT

!  VLIDORT Fixed Sunrays Inputs:

      VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR  = FLUX_FACTOR

!  VLIDORT Modified Sunrays Inputs

      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES = N_SZANGLES
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:N_SZANGLES) = SZANGLES(1:N_SZANGLES)

!  VLIDORT Fixed UserValues Inputs

      VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS = N_USER_LEVELS

!  VLIDORT Modified UserValues Inputs

      !IF ( .NOT. DO_NO_AZIMUTH ) THEN
         VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS = N_USER_RELAZMS
         VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS) = &
                                   USER_RELAZMS(1:N_USER_RELAZMS) 
      !ENDIF
      IF ( DO_USER_VZANGLES ) THEN
         VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = N_USER_VZANGLES
         VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:N_USER_VZANGLES) = &
                                   USER_VZANGLES(1:N_USER_VZANGLES)
      ENDIF
      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS) = &
                                USER_LEVELS(1:N_USER_LEVELS)

!  VLIDORT Fixed Optical Inputs

      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = LAMBERTIAN_ALBEDO
      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH  = ATMOS_WAVELENGTH

!  VLIDORT Fixed Write Inputs

      VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE          = DO_DEBUG_WRITE
      VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT          = DO_WRITE_INPUT
      VLIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME    = INPUT_WRITE_FILENAME
      VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO       = DO_WRITE_SCENARIO
      VLIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME = SCENARIO_WRITE_FILENAME
      VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER        = DO_WRITE_FOURIER
      VLIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME  = FOURIER_WRITE_FILENAME
      VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS        = DO_WRITE_RESULTS
      VLIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME  = RESULTS_WRITE_FILENAME

!  Normal return

!mick fix
      NMESSAGES = NM

      RETURN

!  Line read error - abort immediately

998   CONTINUE
      NM = NM + 1
      STATUS       = VLIDORT_SERIOUS
      MESSAGES(NM) = 'Read failure for entry below String: ' //Trim(Adjustl(PAR_STR))
      ACTIONS(NM)  = 'Re-set value: Entry wrongly formatted in Input file'
      NMESSAGES    = NM

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_READ_INPUTS

!

      SUBROUTINE VLIDORT_CHECK_INPUT_DIMS &
      ( VLIDORT_FixIn, VLIDORT_ModIn, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

!  module, dimensions and numbers
!   Rob Fix 3/18/15. Added MAX_GEOMETRIES to this list..

      USE VLIDORT_pars_m, only : MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAX_USER_LEVELS, &
                                 MAX_SZANGLES, MAXSTREAMS, MAXLAYERS, MAX_GEOMETRIES,  &
                                 MAXMOMENTS_INPUT, MAX_THERMAL_COEFFS, MAXSTOKES,      &
                                 MAX_USER_OBSGEOMS, MAXFINELAYERS, MAX_MESSAGES,       &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS

      USE VLIDORT_Inputs_def_m

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs)   , INTENT (IN) :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (IN) :: VLIDORT_ModIn

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS
      NM = NMESSAGES

!  Check VLIDORT input dimensions against maximum dimensions
!  =========================================================

!  1a. Basic dimensions - always checked

      IF ( VLIDORT_FixIn%Cont%TS_NSTOKES .GT. MAXSTOKES ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number Stokes parameters NSTOKES > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value to 4 or less'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

      IF ( VLIDORT_FixIn%Cont%TS_NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams NSTREAMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

      IF ( VLIDORT_FixIn%Cont%TS_NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of layers NLAYERS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXLAYERS dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

      IF ( VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS .GT. MAX_USER_LEVELS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of user vertical output levels N_USER_LEVELS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_LEVELS dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

      IF ( VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Expansion coefficients NGREEK_MOMENTS_INPUT > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXMOMENTS_INPUT dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

!  1b. Basic dimensions - conditionally checked

      IF ( VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING  ) THEN
        IF ( VLIDORT_FixIn%Cont%TS_NFINELAYERS .GT. MAXFINELAYERS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of fine layers NFINELAYERS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAXFINELAYERS dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION ) THEN
        IF ( VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS .GT. MAX_THERMAL_COEFFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of thermal coefficients" > allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_THERMCOEFFS dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  2a. Geometry dimensions - always checked

      IF ( VLIDORT_ModIn%MSunrays%TS_N_SZANGLES .GT. MAX_SZANGLES) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar zenith angles N_SZANGLES > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_SZANGLES dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

      IF ( VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of relative azimuths N_USER_RELAZMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension in VLIDORT_PARS'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

!  2b. Geometry dimensions - conditionally checked

      IF ( VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES ) THEN
        IF ( VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES .GT. MAX_USER_VZANGLES ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of user streams N_USER_VZANGLES > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_VZANGLES dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) THEN
        IF ( VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Observation Geometries > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_OBSGEOMS dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Rob Fix 3/18/15. Geometry dimension to be checked for the Lattice case
!  4/15/20. Add DO_DOUBLET_GEOMETRY flag

      IF ( .not. VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY &
           .and. .not.VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY ) THEN
         IF ( MAX_GEOMETRIES.ne.MAX_SZANGLES*MAX_USER_VZANGLES*MAX_USER_RELAZMS ) then
          NM = NM + 1
          MESSAGES(NM) = 'Lattice Case: MAX_GEOMETRIES not equal to MAX_SZANGLES*MAX_USER_VZANGLES*MAX_USER_RELAZMS'
          ACTIONS(NM)  = 'Re-set  MAX_GEOMETRIES dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Rob Fix 4/15/20. Geometry dimension to be checked for the Doublet case

      IF ( VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY ) THEN
         IF ( MAX_GEOMETRIES.ne.MAX_SZANGLES*MAX_USER_VZANGLES ) then
          NM = NM + 1
          MESSAGES(NM) = 'Doublet Case: MAX_GEOMETRIES not equal to MAX_SZANGLES*MAX_USER_VZANGLES'
          ACTIONS(NM)  = 'Re-set  MAX_GEOMETRIES dimension in VLIDORT_PARS'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CHECK_INPUT_DIMS

!
      
      SUBROUTINE VLIDORT_CHECK_INPUT &
        ( VLIDORT_FixIn, VLIDORT_ModIn,           & ! InOut
          STATUS, NMESSAGES, MESSAGES, ACTIONS )    ! Output

!  Check the inputs

!  4/15/20. Version 2.8.2
!  ===========================

!  Use type structures directly for the checking

!  PREVIOUS VERSIONS - NOTES
!  =========================

!  %% DO_FULLRAD_MODE argument added (First line). R. Spurr, 05 March 2013
!  %%   Needed to ensure MS-only output in all cases when flagged
      
!  Order of Taylor series (including terms up to EPS^n). Introduced 2/19/14 for Version 2.7
!  %% Revision, 3/3/17.   Complete overhaul of FO and SSCORR flags (Version 2.8)

!  Add Flag DO_WLADJUSTED_OUTPUT (Water-leaving output). 3/18/19 for Version 2.8.1. 
!  Add Flag DO_EXTERNAL_WLEAVE   (Water-leaving output). 4/22/19 for Version 2.8.1. 


      USE VLIDORT_PARS_m, Only : MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ, MAX_SZANGLES, MAX_USER_RELAZMS, &
                                 MAX_USER_VZANGLES, MAX_USER_LEVELS, MAX_USER_OBSGEOMS, MAX_USER_STREAMS,   &
                                 MAX_ALLSTRMS_P1, MAX_TAYLOR_TERMS, MAX_MESSAGES,                           &
                                 VLIDORT_WARNING, VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO, ONE, OMEGA_SMALLNUM

      IMPLICIT NONE

!  Type structure inputs for checking
!  ----------------------------------

      TYPE(VLIDORT_Fixed_Inputs), INTENT(IN)       :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT(INOUT) :: VLIDORT_ModIn

!  Exception handling
!  ------------------

!   Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER, INTENT (OUT) ::              STATUS
      INTEGER, INTENT (INOUT) ::            NMESSAGES
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGES( 0:MAX_MESSAGES )
      CHARACTER (LEN=*), INTENT (INOUT) ::  ACTIONS ( 0:MAX_MESSAGES )

!  Local variables

      INTEGER ::           I, N, UTA, NSTART, NALLSTREAMS, NM, N_OUT_STREAMS
      INTEGER ::           NSTREAMS, NUSERG, NLAYERS
      INTEGER ::           INDEX_ANGLES   ( MAX_USER_STREAMS )
      DOUBLE PRECISION  :: XT, ALL_ANGLES ( MAX_USER_STREAMS )
      DOUBLE PRECISION  :: OUT_ANGLES     ( MAX_USER_STREAMS )

      CHARACTER (LEN=2) :: C2
      LOGICAL           :: LOOP, DO_FOCORR_ALONE

!  Initialize output status

      STATUS = VLIDORT_SUCCESS

!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Check of VLIDORT Basic Input'
!      ACTIONS(0)      = 'No Action required for this Task'

      NM = NMESSAGES

!  Define DO_FOCORR_ALONE flag before input checks

      DO_FOCORR_ALONE = ( .NOT.VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE .AND. VLIDORT_ModIn%MBool%TS_DO_FOCORR )

!  FOCORR_ALONE flag cancels some other flags. Not the DB correction !!!!
!    (Revision, version 2.8, new name for FOCORR_ALONE, used to be DO_SSFULL)

      IF ( DO_FOCORR_ALONE ) THEN
        VLidort_ModIn%MBool%TS_DO_DOUBLE_CONVTEST  = .FALSE.
        VLidort_ModIn%MBool%TS_DO_SOLUTION_SAVING  = .FALSE.
        VLidort_ModIn%MBool%TS_DO_BVP_TELESCOPING  = .FALSE.
        VLidort_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT = .FALSE.
        VLidort_ModIn%MBool%TS_DO_MVOUT_ONLY       = .FALSE.
      ENDIF

!  some proxies

      NLAYERS        = VLidort_FixIn%Cont%TS_NLAYERS
      NSTREAMS       = VLidort_FixIn%Cont%TS_NSTREAMS
      NUSERG         = VLidort_ModIn%MUserVal%TS_N_USER_OBSGEOMS

!  Automatic input
!    Flux factor set to unity. (21 December 2005)

!      DO_ALL_FOURIER        = .FALSE.
!      DO_CLASSICAL_SOLUTION = .TRUE.  ! Removed for Version 2.8
!      DO_DIRECT_BEAM        = .TRUE.

!  Check top level Solar/Thermal options, set warnings
!  ===================================================

!  Check thermal or Solar sources present

      IF ( .NOT.VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES.AND..NOT.VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: DO_SOLAR_SOURCES,DO_THERMAL_EMISSION both not set'
        ACTIONS(NM)  = 'Abort: must set DO_SOLAR_SOURCES and/or DO_THERMAL_EMISSION'
        STATUS = VLIDORT_SERIOUS
      ENDIF

!  Switch off several flags with thermal-only option
!    Set default, regardless of whether solar sources are on. Revised, 2.8
!mick fix 9/19/2017 - turned off

      !IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
      ! IF ( DO_FOCORR .and. DO_FOCORR_NADIR ) THEN
      !   NM = NM + 1
      !   MESSAGES(NM) = 'Switch off FO Nadir correction, not needed for thermal-only'
      !   ACTIONS(NM)  = 'Warning: FOCORR_NADIR correction flag turned off internally'
      !   STATUS = VLIDORT_WARNING
      !   DO_FOCORR_NADIR = .FALSE.
      !  ENDIF
      !ENDIF

!  Set number of solar angles N_SZANGLES to 1 for the thermal-only default

      IF ( .NOT.VLidort_ModIn%MBool%TS_DO_SOLAR_SOURCES.AND.VLidort_FixIn%Bool%TS_DO_THERMAL_EMISSION ) THEN
        IF ( VLidort_ModIn%MSunRays%TS_N_SZANGLES .NE. 1 ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: N_SZANGLES not set to 1 for thermal-only'
         ACTIONS(NM) = 'Warning: N_SZANGLES is set to 1 internally'
         STATUS = VLIDORT_WARNING
         VLidort_ModIn%MSunRays%TS_N_SZANGLES = 1
        ENDIF
      ENDIF

!  Set number of azimuths N_USER_RELAZMS to 1 for the thermal-only default

      IF ( .NOT.VLidort_ModIn%MBool%TS_DO_SOLAR_SOURCES.AND.VLidort_FixIn%Bool%TS_DO_THERMAL_EMISSION ) THEN
        IF ( VLidort_ModIn%MUserVal%TS_N_USER_RELAZMS .NE. 1 ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: N_USER_RELAZMS not set to 1 for thermal-only'
         ACTIONS(NM)  = 'Warning: N_USER_RELAZMS is set to 1 internally'
         STATUS = VLIDORT_WARNING
         VLidort_ModIn%MUserVal%TS_N_USER_RELAZMS = 1
        ENDIF
      ENDIF

!  New code for the Observational Geometry
!  ---------------------------------------

      IF ( VLidort_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) THEN

        IF ( .NOT.VLidort_ModIn%MBool%TS_DO_SOLAR_SOURCES ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: DO_SOLAR_SOURCES not set for Observation Geometry option'
          ACTIONS(NM)  = 'Abort: must set DO_SOLAR_SOURCE'
          STATUS = VLIDORT_SERIOUS
        ENDIF

!%%% No reason you cannot have Observation Geometry with Cross-Over
!%%% Clause commented out, 05 March 2013
!        IF ( DO_THERMAL_EMISSION ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = &
!            'Bad input: DO_THERMAL_EMISSION should not be set for Observation Geometry option'
!          ACTIONS(NM)  = &
!            'Abort: must turn off DO_THERMAL EMISSION'
!          STATUS = VLIDORT_SERIOUS
!        ENDIF
!%%% End Clause commented out, 05 March 2013

!  Observational Geometry control
!   New, 25 October 2012
!     ---- Automatic setting of N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS,
!          and DO_USER_VZANGLES

         NUSERG = VLidort_ModIn%MUserVal%TS_N_USER_OBSGEOMS

         VLidort_ModIn%MSunrays%TS_N_SZANGLES      = NUSERG
         VLidort_ModIn%MUserVal%TS_N_USER_VZANGLES = NUSERG
         VLidort_ModIn%MUserVal%TS_N_USER_RELAZMS  = NUSERG
         VLidort_ModIn%MBool%TS_DO_USER_VZANGLES   = .true.

!     ---- Automatic setting of BEAM_SZAS, USER_ANGLES_INPUT, and USER_RELAZMS

         VLidort_ModIn%MSunrays%TS_SZANGLES           (1:NUSERG) = VLidort_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:NUSERG,1)
         VLidort_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:NUSERG) = VLidort_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:NUSERG,2)
         VLidort_ModIn%MUserVal%TS_USER_RELAZMS       (1:NUSERG) = VLidort_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:NUSERG,3)

      ENDIF

!  Check Taylor series parameter is not out of range.
!   This code was added 2/19/14 for Version 2.7

      IF ( (VLidort_FixIn%Cont%TS_TAYLOR_ORDER .GT. MAX_TAYLOR_TERMS-2) &
               .OR. (VLidort_FixIn%Cont%TS_TAYLOR_ORDER .LT. 0) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Taylor series Order parameter out of range'
        ACTIONS(NM)  = 'Re-set input value, should be in range 0-4'
        STATUS       = VLIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

!  Check inputs (both file-read and derived)
!  -----------------------------------------

!  Check Chapman function options

      IF ( VLidort_ModIn%MBool%TS_DO_SOLAR_SOURCES ) THEN
       IF ( .NOT. VLidort_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION ) THEN
        IF ( VLidort_FixIn%Bool%TS_DO_PLANE_PARALLEL ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Chapman Function not set, plane parallel'
          ACTIONS(NM)  = 'Warning: Chapman function set internally'
          STATUS       = VLIDORT_WARNING
          VLidort_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION = .TRUE.
        ELSE
          NM = NM + 1
          MESSAGES(NM) = 'Chapman Function not set, pseudo-spherical treatment'
          ACTIONS(NM)  = 'Have you set the CHAPMAN_FACTORS values?'
          STATUS = VLIDORT_SERIOUS
        ENDIF
       ENDIF
      ENDIF

!  Check sphericity corrections for Solar....Cannot both be turned on
!    --------- New code 31 January 2007
!mick fix 9/19/2017 - turned off "DO_SOLAR_SOURCES" IF condition
!                   - included additional DO_FOCORR_NADIR/DO_FOCORR_OUTGOING check 

!  6/17/20. Bug. Must allow FOCORR to operate with PLANE_PARALLEL possibility

      !IF ( DO_SOLAR_SOURCES ) THEN
        IF ( VLidort_ModIn%MBool%TS_DO_FOCORR ) THEN

          IF ( VLidort_ModIn%MBool%TS_DO_FOCORR_NADIR .AND. VLidort_ModIn%MBool%TS_DO_FOCORR_OUTGOING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Cannot have both FO corrections on for Solar single-scatter'
            ACTIONS(NM)  = 'Turn off DO_FOCORR_NADIR and/or DO_FOCORR_OUTGOING'
            STATUS = VLIDORT_SERIOUS
          ENDIF

! 6/17/20 comment out this clause

!          IF (      ( VLidort_ModIn%MBool%TS_DO_FOCORR_NADIR .AND. VLidort_ModIn%MBool%TS_DO_FOCORR_OUTGOING ) .OR. &
!               .NOT.( VLidort_ModIn%MBool%TS_DO_FOCORR_NADIR .OR.  VLidort_ModIn%MBool%TS_DO_FOCORR_OUTGOING )         )  THEN
!            NM = NM + 1
!            MESSAGES(NM) = 'Cannot have both FO correction types on or off when FO correction is in use'
!            ACTIONS(NM)  = 'Set one of DO_FOCORR_NADIR or DO_FOCORR_OUTGOING'
!            STATUS = VLIDORT_SERIOUS
!          ENDIF

        ENDIF
      !ENDIF

!  New for Version 2.8, Check surface leaving inputs
!  -------------------------------------------------

      IF ( VLidort_FixIn%Bool%TS_DO_SURFACE_LEAVING) THEN

!  Check compatibility. 1/31/16

        IF ( VLidort_FixIn%Bool%TS_DO_WATER_LEAVING .and. VLidort_FixIn%Bool%TS_DO_FLUORESCENCE ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Cannot have both Water-leaving and Fluorescence turned on'
          ACTIONS(NM)    = 'Turn off DO_WATER_LEAVING or DO_FLUORESCENCE'
          STATUS = VLIDORT_SERIOUS
        ENDIF

        IF ( .not.VLidort_FixIn%Bool%TS_DO_WATER_LEAVING .and. .not.VLidort_FixIn%Bool%TS_DO_FLUORESCENCE ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Surface-leaving is ON, but both Water-leaving and Fluorescence are OFF!'
          ACTIONS(NM)    = 'Turn on either DO_WATER_LEAVING or DO_FLUORESCENCE'
          STATUS = VLIDORT_SERIOUS
        ENDIF

!  Check Water-leaving output flag. 3/18/19 for Version 2.8.1
       
        IF ( VLidort_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT .and..not.VLidort_FixIn%Bool%TS_DO_WATER_LEAVING  ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Water-leaving output is ON, but main Water-leaving Flag is OFF!'
          ACTIONS(NM)    = 'Turn on DO_WATER_LEAVING if you want to use DO_WLADJUSTED_OUTPUT'
          STATUS = VLIDORT_SERIOUS
        ENDIF

!  Check Use of external water-leaving source. 4/22/19 for Version 2.8.1
        
        IF ( VLidort_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE .and..not.VLidort_FixIn%Bool%TS_DO_WATER_LEAVING  ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'External Water-leaving source flag is ON, but main Water-leaving Flag is OFF!'
          ACTIONS(NM)    = 'Turn on DO_WATER_LEAVING if you want to use DO_EXTERNAL_WLEAVE'
          STATUS = VLIDORT_SERIOUS
        ENDIF
        IF ( VLidort_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE .and. &
               ( VLidort_FixIn%Bool%TS_DO_WATER_LEAVING .and. VLidort_FixIn%Bool%TS_DO_TF_ITERATION ) ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Both External Water-leaving source flag and TF_ITERATION flag are on'
          ACTIONS(NM)    = 'Turn off DO_TF_ITERATION if you want to use DO_EXTERNAL_WLEAVE'
          STATUS = VLIDORT_SERIOUS
        ENDIF

!  Check TF conditions. 7/8/16

        IF ( VLidort_FixIn%Bool%TS_DO_WATER_LEAVING .and. VLidort_FixIn%Bool%TS_DO_TF_ITERATION ) THEN
          IF ( VLidort_FixIn%Cont%TS_TF_MAXITER .gt. 10 ) then
            NM = NM + 1
            MESSAGES(NM)   = 'Water-leaving Transmittance calculation, cannot have more than 10 iterations'
            ACTIONS(NM)    = 'Set input TF_MAXITER to 10 or less'
            STATUS = VLIDORT_SERIOUS
          ENDIF
          IF ( VLidort_FixIn%Cont%TS_TF_CRITERION .gt. 0.01d0 .or. VLidort_FixIn%Cont%TS_TF_CRITERION.lt. 1.0d-07 ) then
            NM = NM + 1
            MESSAGES(NM)   = 'Water-leaving Transmittance calculation, criterion not in range [0.0000001 to 0.01]'
            ACTIONS(NM)    = 'Set input TF_CRITERION in range [0.0000001,0.01]'
            STATUS = VLIDORT_SERIOUS
          ENDIF
        ENDIF

      ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  05 March 2013  Start %%%%%%%%%%%%%%%%%

!   New check to make sure output is MS-only when DO_FULLRAD_MODE is off
!     -- Do not want Internal FO Corrections in this case.
!     -- Turn off Internal FO Corrections, with Warning.

!  Old code, pre 2.8
!      IF ( .not.DO_SS_EXTERNAL .and. .not.DO_FULLRAD_MODE ) then
!        IF ( DO_SSCORR_NADIR .or. DO_SSCORR_OUTGOING ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Internal FO calculation must be Off when MS-only is desired'
!          ACTIONS(NM)  = 'DO_SSCORR_NADIR/DO_SSCORR_OUTGOING flags turned off Internally'
!          STATUS = VLIDORT_WARNING
!          DO_SSCORR_NADIR    = .false.
!          DO_SSCORR_OUTGOING = .false.
!        ENDIF
!      ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  05 March 2013 Finish %%%%%%%%%%%%%%%%

!   R. Spurr, RT SOLUTIONS Inc.  Revised for Version 2.8
!mick mod 9/19/2017 - turned off now to allow tightening of some VLIDORT input
!                     definitions related to type of rad solution VLIDORT returns in
!                     different cases

      !IF ( DO_FOCORR .and. .not.DO_FULLRAD_MODE ) then
      !   NM = NM + 1
      !   MESSAGES(NM) = 'Internal FO calculation must be off when MS-only is desired'
      !   ACTIONS(NM)  = 'DO_FOCORR/DO_FOCORR_NADIR/DO_FOCORR_OUTGOING flags turned off Internally'
      !   STATUS = VLIDORT_WARNING
      !   DO_FOCORR          = .false.
      !   DO_FOCORR_NADIR    = .false.
      !   DO_FOCORR_OUTGOING = .false.
      !ENDIF

!  New 15 March 2012
!    Turn off FOCORR flags if the external FO calculation applies

      IF ( VLidort_ModIn%MBool%TS_DO_FOCORR_EXTERNAL ) then
        IF ( VLidort_ModIn%MBool%TS_DO_FOCORR_NADIR ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'External FO calculation: Cannot have Nadir single scatter correction'
          ACTIONS(NM)  = 'Turn off DO_FOCORR_NADIR flag'
          STATUS = VLIDORT_SERIOUS
        ENDIF
        IF ( VLidort_ModIn%MBool%TS_DO_FOCORR_OUTGOING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'External FO calculation: Cannot have Outgoing single scatter correction'
          ACTIONS(NM)  = 'Turn off DO_FOCORR_OUTGOING flag'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  New 02 Jul 2013. Check removed, 19 March 2015, with advent of Lattice option
!  Check dedicated FO code / ObsGeo mode compatibility

!      IF (DO_FO_CALC) THEN
!        IF ( .NOT. DO_OBSERVATION_GEOMETRY ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Internal SS calculation: Must set Observation Geometry option'
!          ACTIONS(NM)  = 'Turn on DO_OBSERVATION_GEOMETRY flag'
!          STATUS = VLIDORT_SERIOUS
!        ENDIF
!      ENDIF

!  Check beam mode operation

      IF ( VLidort_ModIn%MBool%TS_DO_SOLAR_SOURCES ) THEN
       IF ( VLidort_FixIn%Bool%TS_DO_PLANE_PARALLEL ) THEN
        IF ( VLidort_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: plane-parallel and refractive flags both set'
         ACTIONS(NM)  = 'Warning: turn off Refraction internally'
         STATUS       = VLIDORT_WARNING
         VLidort_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = .FALSE.
        ENDIF
       ENDIF
      ENDIF

!  Check consistency of mean value input control
!  ---------------------------------------------

      IF ( VLidort_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT ) THEN
        IF ( VLidort_ModIn%MBool%TS_DO_MVOUT_ONLY ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Cannot have both mean-value flags set'
          ACTIONS(NM)  = 'Warning: disable DO_MVOUT_ONLY flag internally'
          STATUS = VLIDORT_WARNING
          VLidort_ModIn%MBool%TS_DO_MVOUT_ONLY = .FALSE.
        ENDIF
      ENDIF

!  Removed DO_NO_AZIMUTH check. 17 January 2006

      IF ( .NOT.VLidort_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT ) THEN
        IF ( VLidort_ModIn%MBool%TS_DO_MVOUT_ONLY ) THEN
          IF ( VLidort_ModIn%MBool%TS_DO_USER_VZANGLES ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Mean-value option needs quadratures only'
            ACTIONS(NM)  = 'Warning: DO_USER_VZANGLES flag disabled internally'
            STATUS = VLIDORT_WARNING
            VLidort_ModIn%MBool%TS_DO_USER_VZANGLES = .FALSE.
          ENDIF
        ENDIF
      ENDIF

!  Check consistency of multiple scatter source term output control
!   ---Specialist options. Removed 30 March 2007.
!      IF ( SAVE_LAYER_MSST ) THEN
!        IF ( .NOT. DO_USER_VZANGLES ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Bad input: MSCAT. source term - USER_VZANGLES flag not set'
!          ACTIONS(NM)  = 'Check DO_USER_VZANGLES and SAVE_LAYER_MSST flags'
!          STATUS = VLIDORT_SERIOUS
!        ENDIF
!      ENDIF

!  Check consistency of BVP_TELESCOPING and SOLUTION_SAVING flags
!  ---Warning. Set solution-saving internally

      IF (VLidort_ModIn%MBool%TS_DO_BVP_TELESCOPING .AND. .NOT.VLidort_ModIn%MBool%TS_DO_SOLUTION_SAVING) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: BVP telescoping -> solution saving must be set'
        ACTIONS(NM)  = 'Warning:  Solution saving was set internally'
        STATUS = VLIDORT_WARNING
        VLidort_ModIn%MBool%TS_DO_SOLUTION_SAVING = .TRUE.
      ENDIF

!  Check consistency of Rayleigh-only and Isotropic-only cases
!   Removed 17 January 2006, Isotropic option removed.

!      IF ( VLidort_ModIn%MBool%TS_DO_RAYLEIGH_ONLY .AND. VLidort_ModIn%MBool%TS_DO_ISOTROPIC_ONLY ) THEN
!        NM = NM + 1
!        MESSAGES(NM) = 'Bad input: Isotropic_only & Rayleigh-only flags both set'
!        ACTIONS(NM)  = 'Check DO_RAYLEIGH_ONLY and DO_ISOTROPIC_ONLY flags'
!        STATUS = VLIDORT_SERIOUS
!      ENDIF

!  No Delta-M scaling with Rayleigh only
!   ---Warning. Turn off delta-M scaling.

      IF ( VLidort_ModIn%MBool%TS_DO_RAYLEIGH_ONLY ) THEN
        IF ( VLidort_ModIn%MBool%TS_DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: No delta-M scaling with Rayleigh-only'
          ACTIONS(NM)  = 'Warning: DO_DELTAM_SCALING turned off internally'
          STATUS = VLIDORT_WARNING
          VLidort_ModIn%MBool%TS_DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  -----------Note the following in the scalar code -------------------
!  No Delta-M scaling with Isotropic only
!   ---Warning. Turn off delta-M scaling.
!   Removed 17 January 2006, Isotropic option removed.

!      IF ( VLidort_ModIn%MBool%TS_DO_ISOTROPIC_ONLY ) THEN
!        IF ( VLidort_ModIn%MBool%TS_DO_DELTAM_SCALING ) THEN
!         NM = NM + 1
!          MESSAGES(NM) = 'Bad input: No delta-M scaling with Isotropic-only'
!          ACTIONS(NM)  = 'Warning: DO_DELTAM_SCALING turned off internally'
!          STATUS = VLIDORT_WARNING
!          VLidort_ModIn%MBool%TS_DO_DELTAM_SCALING = .FALSE.
!        ENDIF
!      ENDIF

!  Check directional input

      IF ( .NOT.VLidort_FixIn%Bool%TS_DO_UPWELLING .AND. .NOT. VLidort_FixIn%Bool%TS_DO_DNWELLING ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: no directional input is set'
        ACTIONS(NM)  = 'Check DO_UPWELLING & DO_DNWELLING: one must be set!'
        STATUS = VLIDORT_SERIOUS
      ENDIF

!  Check number of input GSF expansion coefficient moments (non-Rayleigh)
!  ----------------------------------------------------------------------

!  Checks for general scattering case
!    Isotropic part Removed 17 January 2006, Isotropic option removed.

      IF ( .NOT.VLidort_ModIn%MBool%TS_DO_RAYLEIGH_ONLY .AND. .NOT. DO_FOCORR_ALONE ) THEN
        IF ( VLidort_ModIn%MBool%TS_DO_DELTAM_SCALING ) THEN
          IF ( VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT.LT.2*NSTREAMS ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Fewer than 2N expansion moments with delta-M'
            ACTIONS(NM)  = 'Warning: Re-set NGREEK_MOMENTS_INPUT to 2N internally'
            STATUS = VLIDORT_WARNING
            VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 2*NSTREAMS
          ENDIF
        ELSE
!          IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN   !  Changed, Version 2.8
          IF ( VLidort_ModIn%MBool%TS_DO_FOCORR ) THEN
            IF ( VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT.LT.2*NSTREAMS-1 ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: Fewer than 2N-1 expansion moments without delta-M'
              ACTIONS(NM)  = 'Warning: Re-set NGREEK_MOMENTS_INPUT to 2N-1 internally'
              STATUS = VLIDORT_WARNING
              VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 2*NSTREAMS - 1
            ENDIF
          ENDIF
        ENDIF

      ELSE

!  Checks for Rayleigh-only option
!   All warnings.

        IF ( VLidort_ModIn%MBool%TS_DO_RAYLEIGH_ONLY ) THEN

          IF ( VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT.NE.2 ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Rayleigh-only, expansion momemts NOT = 2'
            ACTIONS(NM)  = 'Warning: Set NGREEK_MOMENTS_INPUT = 2 internally'
            STATUS       = VLIDORT_WARNING
            VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 2
          ENDIF

          IF ( VLidort_ModIn%MBool%TS_DO_BVP_TELESCOPING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Bvp telescoping not possible, Rayleigh only'
            ACTIONS(NM)  = 'Warning: Turn off BVP_TELESCOPING internally'
            STATUS       = VLIDORT_WARNING
            VLidort_ModIn%MBool%TS_DO_BVP_TELESCOPING = .FALSE.
          ENDIF

          IF ( VLidort_ModIn%MBool%TS_DO_SOLUTION_SAVING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Solution saving not possible, Rayleigh only'
            ACTIONS(NM)  = 'Warning: Turn off SOLUTION_SAVING internally'
            STATUS       = VLIDORT_WARNING
            VLidort_ModIn%MBool%TS_DO_SOLUTION_SAVING = .FALSE.
          ENDIF

        ENDIF

!  Checks for Isotropic only option. Removed, 17 January 2006. (Still present in Scalar LIDORT code)

      ENDIF

!  Reset solution saving and BVP telescoping flags
!  Do not need the isotropic-only options
!      DO_SOLUTION_SAVING =  ( DO_SOLUTION_SAVING .AND. ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )
!      DO_BVP_TELESCOPING =  ( DO_BVP_TELESCOPING .AND. ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )

      VLidort_ModIn%MBool%TS_DO_SOLUTION_SAVING = &
           ( VLidort_ModIn%MBool%TS_DO_SOLUTION_SAVING .AND. .NOT.VLidort_ModIn%MBool%TS_DO_RAYLEIGH_ONLY )
      VLidort_ModIn%MBool%TS_DO_BVP_TELESCOPING = &
           ( VLidort_ModIn%MBool%TS_DO_BVP_TELESCOPING .AND. .NOT.VLidort_ModIn%MBool%TS_DO_RAYLEIGH_ONLY )

!  BVP telescoping doesn't work with non-Lambertian surfaces
!   Condition relaxed, Version 2.8, now possible.
!      IF (  DO_BVP_TELESCOPING ) THEN
!        IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'BVP telescoping disabled for non-Lambertian'
!          ACTIONS(NM)  = 'Turn off DO_BVP_TELESCOPING flag, done internally'
!          STATUS = VLIDORT_WARNING
!          DO_BVP_TELESCOPING = .FALSE.
!        ENDIF
!      ENDIF

!  Check azimuth-only conditions
!  -----------------------------

!  Check no-Azimuth flag. Now set internally
!    ---WARNING. Do-no-Azimuth Flag turned on
!      IF ( .NOT.DO_NO_AZIMUTH ) THEN
!        IF ( DO_USER_VZANGLES. AND. N_USER_VZANGLES.EQ.1 ) THEN
!          IF ( USER_ANGLES_INPUT(1) .EQ. ZERO ) THEN
!            NM = NM + 1
!            MESSAGES(NM) = 'Bad input: zenith-sky output requires no azimuth'
!            ACTIONS(NM)  = 'Warning: DO_NO_AZIMUTH flag set true internally'
!            STATUS = VLIDORT_WARNING
!          ENDIF
!        ENDIF
!      ENDIF

!  Check: OLD single scattering correction and Do Rayleigh
!    ---WARNING. SS Flag turned off
!  Check only required for the diffuse field calculations (Version 2.3)

!  @@@ Rob Fix 5/20/13. NO turn-off DO_SSCORR_NADIR for Rayleigh-BRDF situation
!     if DO_RAYLEIGH_ONLY.and..NOT.DO_LAMBERTIAN_SURFACE and.not.DO_SS_EXTERNAL,
!     then one of the SSCORR flags should be set internally
!  Old Code -------
!      IF ( DO_SSCORR_NADIR ) THEN
!        IF ( DO_RAYLEIGH_ONLY .AND. .NOT. DO_SSFULL ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Bad input: No SS correction for Rayleigh only'
!          ACTIONS(NM)  = 'Warning: DO_SSCORR_NADIR turned off internally'
!          STATUS = VLIDORT_WARNING
!          DO_SSCORR_NADIR = .FALSE.
!        ENDIF
!      ENDIF
!  New code --------  ! Rob Fix 8/18/15
!      IF ( DO_RAYLEIGH_ONLY .AND. .NOT.DO_LAMBERTIAN_SURFACE ) THEN
!        IF ( .NOT.DO_SSCORR_OUTGOING .OR. .NOT.DO_SSCORR_NADIR ) THEN
!          IF ( .NOT.DO_SS_EXTERNAL ) THEN
!            if ( .not.DO_SSCORR_OUTGOING ) then
!              NM = NM + 1
!              MESSAGES(NM) = 'Bad input: Rayleigh-only+BRDF: set an SSCORR flag!'
!              ACTIONS(NM)  = 'Warning: DO_SSCORR_NADIR turned on internally'
!              STATUS = VLIDORT_WARNING
!              DO_SSCORR_NADIR = .true.   ! Rob Fix 8/18/15
!            endif
!          ENDIF
!        ENDIF
!      ENDIF
!  @@@ Rob Fix 5/20/13. End of Fix ===================================

!  Version 2.8, Fix 9/19/16. 3/1/17.

      IF ( VLidort_ModIn%MBool%TS_DO_RAYLEIGH_ONLY .AND. .NOT.VLidort_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE ) THEN
        IF ( .NOT.VLidort_ModIn%MBool%TS_DO_FOCORR .AND. .NOT.VLidort_ModIn%MBool%TS_DO_FOCORR_EXTERNAL ) THEN
           NM = NM + 1
           MESSAGES(NM) = 'Bad input: Rayleigh-only+BRDF: Need to set FOCORR flags!'
           ACTIONS(NM)  = 'Warning: DO_FOCORR and DO_FOCORR_NADIR turned on internally'
           STATUS = VLIDORT_WARNING
           VLidort_ModIn%MBool%TS_DO_FOCORR       = .true. 
           VLidort_ModIn%MBool%TS_DO_FOCORR_NADIR = .true. 
        ENDIF
      ENDIF

!  Just the single scatter, enabled 25 September 2007. Name changed, Version 2.8
!   Single scatter corrections must be turned on
!   No longer required, Version 2.8
!      IF ( DO_SSCORR_ALONE ) THEN
!        IF ( .NOT.DO_SSCORR_NADIR .AND. .NOT.DO_SSCORR_OUTGOING ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Bad input: Full SS, must have one SSCORR flag set'
!          ACTIONS(NM)  = 'Full SS: default to use outgoing SS correction'
!          STATUS = VLIDORT_WARNING
!          DO_SSCORR_NADIR    = .FALSE.
!          DO_SSCORR_OUTGOING = .TRUE.
!        ENDIF
!      ENDIF

!  Just the FO correction alone, enabled 25 September 2007. Name changed, Version 2.8
!   Diffuse-field Delta-M scaling must be turned off

      IF ( DO_FOCORR_ALONE ) THEN
        IF ( VLidort_ModIn%MBool%TS_DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: First-Order ALONE --> diffuse-field delta-M on'
          ACTIONS(NM)  = 'First-Order ALONE: internal default to deltam_scaling = false'
          STATUS = VLIDORT_WARNING
          VLidort_ModIn%MBool%TS_DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  Check thermal inputs
!  --------------------

!  If thermal transmittance only, check thermal flag

      IF ( VLidort_ModIn%MBool%TS_DO_SOLAR_SOURCES .or. .NOT. VLidort_FixIn%Bool%TS_DO_THERMAL_EMISSION ) THEN
       IF ( VLidort_ModIn%MBool%TS_DO_THERMAL_TRANSONLY ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: No thermal, must turn off transonly flag'
         ACTIONS(NM)  = 'Warning: DO_THERMAL_TRANSONLY turned off internally'
         STATUS = VLIDORT_WARNING
         VLidort_ModIn%MBool%TS_DO_THERMAL_TRANSONLY = .FALSE.
       ENDIF
      ENDIF

!  Switch off a bunch of flags. Revised for 2.8.2
!    Solution saving is required though!

      IF ( VLidort_FixIn%Bool%TS_DO_THERMAL_EMISSION ) THEN
        IF ( VLidort_ModIn%MBool%TS_DO_THERMAL_TRANSONLY ) THEN
          VLidort_ModIn%MBool%TS_DO_RAYLEIGH_ONLY  = .FALSE.
          VLidort_ModIn%MBool%TS_DO_DELTAM_SCALING  = .FALSE.
          VLidort_ModIn%MBool%TS_DO_SOLUTION_SAVING = .TRUE.
        ENDIF
      ENDIF

!  No solar sources for thermal transmittance

      IF ( VLidort_ModIn%MBool%TS_DO_THERMAL_TRANSONLY ) THEN
        IF ( VLidort_ModIn%MBool%TS_DO_SOLAR_SOURCES ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: thermal tranmsittance, must turn off solar'
         ACTIONS(NM)  = 'Warning: DO_SOLAR_SOURCES turned off internally'
         STATUS = VLIDORT_WARNING
         VLidort_ModIn%MBool%TS_DO_SOLAR_SOURCES = .FALSE.
       ENDIF
      ENDIF

!  Check viewing geometry input
!  ----------------------------

!  Check Earth radius (Chapman function only)
!    ---WARNING. Default value of 6371.0 will be set

      IF ( VLidort_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION ) THEN
        IF ( .NOT. VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL ) THEN
          IF ( VLIDORT_ModIn%MCont%TS_EARTH_RADIUS.LT.6320.0D0 .OR. &
               VLIDORT_ModIn%MCont%TS_EARTH_RADIUS.GT.6420.0D0 ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Earth radius outside of [6320-6420]'
            ACTIONS(NM)  = 'Warning: default value of 6371.0 was set'
            STATUS = VLIDORT_WARNING
            VLIDORT_ModIn%MCont%TS_EARTH_RADIUS = 6371.0D0
          ENDIF
        ENDIF
      ENDIF

!  Check dimensioning on Legendre numbers (refractive geometry only)

      IF ( VLidort_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY ) THEN
        NALLSTREAMS = VLidort_ModIn%MSunrays%TS_N_SZANGLES * NLAYERS + NSTREAMS + VLidort_ModIn%MUserVal%TS_N_USER_VZANGLES
        IF ( NALLSTREAMS .GT. MAX_ALLSTRMS_P1 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Dimensioning error for refractive beam angles'
          ACTIONS(NM)  = 'Increase dimension MAX_ALLSTRMS_P1'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  Check GEOMETRY_SPECHEIGHT (only for outgoing sphericity correction)
!    GEOMETRY_SPECHEIGHT cannot be greater than HEIGHT_GRID(NLAYERS)

      IF ( VLidort_ModIn%MBool%TS_DO_FOCORR .and. VLidort_ModIn%MBool%TS_DO_FOCORR_OUTGOING ) THEN
        IF ( VLIDORT_ModIn%MCont%TS_GEOMETRY_SPECHEIGHT .GT. VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(NLAYERS) ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'GEOMETRY_SPECHEIGHT must be =< Input BOA-HEIGHT '
          ACTIONS(NM)  = 'Warning: Internal Re-set of GEOMETRY_SPECHEIGHT '
          STATUS = VLIDORT_WARNING
          VLIDORT_ModIn%MCont%TS_GEOMETRY_SPECHEIGHT  = VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(NLAYERS)
        ENDIF
      ENDIF

!  Check solar zenith angle input

!mick hold - 9/26/2012
      !IF (DO_SOLAR_SOURCES ) THEN
        DO I = 1, VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
          IF ( VLIDORT_ModIn%MSunrays%TS_SZANGLES(I) .LT. ZERO .OR. &
               VLIDORT_ModIn%MSunrays%TS_SZANGLES(I) .GE. 90.0D0 ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: out-of-range solar angle, no. '//C2
            ACTIONS(NM)  = 'Look at SZANGLES input, should be < 90 & > 0'
            STATUS = VLIDORT_SERIOUS
          ENDIF
        ENDDO
      !ELSE IF ( .NOT.DO_SOLAR_SOURCES .AND. DO_THERMAL_EMISSION ) THEN
      !  SZANGLES(1:N_SZANGLES) = ZERO
      !ENDIF

!  Check relative azimuths

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS)
        I = I + 1
        IF ( VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(I) .GT. 360.0D0   .OR. &
             VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(I) .LT. ZERO ) THEN
          WRITE(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: out-of-range azimuth angle, no. '//C2
          ACTIONS(NM)  = 'Look at azimuth angle input, should be in [0,360]'
          LOOP = .FALSE.
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDDO

!  Limits on user-defined options

!      IF ( .NOT. DO_USER_VZANGLES .AND..NOT.DO_QUAD_OUTPUT ) THEN
!        NM = NM + 1
!        MESSAGES(NM) = 'Bad input: No angular stream output is specified'
!        ACTIONS(NM)  = 'Check DO_USER_VZANGLES and DO_QUAD_OUTPUT flags'
!        STATUS = VLIDORT_SERIOUS
!      ENDIF

!  Check user-defined stream angles (should always be [0,90])

      IF ( VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES)
          I = I + 1
          IF ( VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(I) .GT. 90.0D0   .OR. &
               VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(I) .LT. ZERO ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: out-of-range user stream, no. '//C2
            ACTIONS(NM)  = 'Look at user viewing zenith angle input'
            LOOP         = .FALSE.
            STATUS       = VLIDORT_SERIOUS
          ENDIF
        ENDDO
      ENDIF

!  Re-order the input angles
!  -------------------------

!  This section has been truncated. 28 March 2007

!mick fix - moved "N_OUT_STREAMS = N_USER_VZANGLES" outside if block
!           & initialized OUT_ANGLES to ZERO in case DO_USER_VZANGLES = .FALSE.

      N_OUT_STREAMS  = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
      OUT_ANGLES     = ZERO
      IF ( VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES ) THEN
        IF ( N_OUT_STREAMS .EQ. 1 ) THEN
          OUT_ANGLES(1) =  VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1)
        ELSE
!%%%%%%%%%%%%%%%% ROB clause 05 March 2013 %%%%%%%%%%%%%%%%%
!%%%%%% Re-ordering only for Lattice computations %%%%%%%%%%

!  revised Call, 3/17/17, with RQSORT_IDX

           IF ( .NOT.VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) THEN
              DO I = 1, VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
                ALL_ANGLES(I)   = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(I)
                INDEX_ANGLES(I) = I                 ! 3/17/17  Need to initialize
              ENDDO
              CALL RQSORT_IDX(N_OUT_STREAMS,ALL_ANGLES,N_OUT_STREAMS,INDEX_ANGLES)
              DO I = 1, N_OUT_STREAMS
                OUT_ANGLES(I) = ALL_ANGLES(INDEX_ANGLES(I))
                VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(I) = OUT_ANGLES(I)
              ENDDO
           ELSE
              DO I = 1, N_OUT_STREAMS
                OUT_ANGLES(I) = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(I)
              ENDDO
           ENDIF

!%%%%%%%%%%%%%%%%%%%%%%% End ROB clause 05 March 2013 %%%%%%%%%%%%%
        ENDIF
      ENDIF

!  Check height grid input (Chapman function only)

      IF ( VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.NLAYERS)
          I = I + 1
          IF ( VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(I-1).LE.VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(I) ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Height-grid not monotonically decreasing; Layer '//C2
            ACTIONS(NM)  = 'Look at Height-grid input'
            LOOP = .FALSE.
            STATUS = VLIDORT_SERIOUS
          ENDIF
        ENDDO
      ENDIF

!  Check vertical outputs
!  ----------------------

!  Check vertical output levels (should always be within atmosphere!)

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS)
        I = I + 1
        IF ( VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(I) .GT. DBLE(NLAYERS) .OR. &
             VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(I) .LT. ZERO )  THEN
          WRITE(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Out of range for level choice # '//C2
          ACTIONS(NM)  = 'Re-set level output '
          LOOP         = .FALSE.
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDDO

!  Check repetition of vertical output choices

      UTA = 0
      LOOP = .TRUE.
      DO WHILE ( LOOP .AND. UTA .LT. VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS )
        UTA = UTA + 1
        XT = VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(UTA)
        NSTART = 0
        DO N = 1, VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS
          IF ( XT .EQ. VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(N)) NSTART = NSTART + 1
        ENDDO
        IF ( NSTART .NE. 1 ) THEN
          LOOP = .FALSE.
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: repetition of vertical output choice'
          ACTIONS(NM)  = 'Re-set level output '
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDDO

!  Check Planetary problem and media flags

!  4/26/19, 9/24/19. Type structure variable Status_inputcheck was not set, now done

!  5/22/20. Version 2.8.2 Upgrades 
!   ==> Relax condition on Rayleigh only, and on FOCORR_NADIR
!              ( Planetary problem works for Aerosols and FOCORR_NADIR )
!   ==> Now, only fails for (dark-surface, Lambertian case, no thermal)

!  Here is the older code..........(Pre 5/22/20)
!      IF ( VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM .or. &
!           VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA(1)   .or. &
!           VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA(2) ) then
!         IF ( .not.VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE.or.VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING .or. &
!              VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION .or. VLIDORT_ModIn%MBool%TS_DO_FOCORR  ) then
!            NM = NM + 1 ; STATUS = VLIDORT_SERIOUS
!            MESSAGES(NM) = 'Media_problem input check not valid'
!            ACTIONS(NM)  = 'Check thermal/BRDF/FOCorr/Sleave flags for this option'           
!         ENDIF
!      ENDIF

!  Here is the New Code..........(5/22/20 Upgrade)

      IF ( VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM .or. &
           VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA(1)   .or. &
           VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA(2) ) then
         IF ( (.not.VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE) .or. VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING .or. &
                    VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION    .or. &
                    (VLIDORT_ModIn%MBool%TS_DO_FOCORR.and.VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING) ) then
            NM = NM + 1 ; STATUS = VLIDORT_SERIOUS
            MESSAGES(NM) = 'Media/Planetary problems: input check not valid'
            ACTIONS(NM)  = 'Turn off ThermalEmission/FOCorr/FOCorr_Outgoing/Sleave flags, Turn on Lambertian flag'           
         ENDIF
      ENDIF

!  Check geophysical scattering inputs
!  -----------------------------------

!  This section has been removed from Input checking, which is now confied to non-optical inputs
!  4/15/20. VLIDORT 2.8.2

!  Check single scatter albedos
!      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
!       DO L = 1, NLAYERS
!        IF ( OMEGA_TOTAL_INPUT(L).GT.ONE-OMEGA_SMALLNUM ) THEN
!          NM = NM + 1 ; WRITE(C3,'(I3)')L ; STATUS = VLIDORT_SERIOUS
!          MESSAGES(NM) = 'Bad input: SS-albedo too close to 1, layer '//C3
!          ACTIONS(NM)  = 'Check SS-albedo input'
!        ELSE IF ( OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
!          NM = NM + 1 ; WRITE(C3,'(I3)')L ; STATUS = VLIDORT_SERIOUS
!          MESSAGES(NM) = 'Bad input: SS-albedo too close to 0, layer '//C3
!          ACTIONS(NM)  = 'Check SS-albedo input'
!        ENDIF
!       ENDDO
!      ENDIF
!  Check first phase function moments
!    Version 2.8, Upgrade, 3/1/17. if not using the Greekmat coefficients in FOCORR
!      IF ( ( .NOT.DO_THERMAL_TRANSONLY ) .OR. .NOT.( DO_FOCORR.AND.DO_SSCORR_USEFMAT ) ) THEN
!       DO L = 1, NLAYERS
!        IF ( GREEKMAT_TOTAL_INPUT(0,L,1).NE.ONE ) THEN
!          NM = NM + 1 ; WRITE(C3,'(I3)')L ; STATUS = VLIDORT_SERIOUS
!          MESSAGES(NM) = 'First phase moment (GREEK_11) not 1 for layer '//C3
!          ACTIONS(NM)  = 'Check First phase function moment'
!        ENDIF
!       ENDDO
!      ENDIF
!  Additional check on non-negativity of (1,1) moments.
!    Furnished 27 September 2007, as a result of input-checking.
!      DO N = 1, NLAYERS
!        DO L = 1, NGREEK_MOMENTS_INPUT
!          IF ( GREEKMAT_TOTAL_INPUT(L,N,1).LT.ZERO ) THEN
!          NM = NM + 1 ; WRITE(C3,'(I3)')N ; STATUS = VLIDORT_SERIOUS
!            MESSAGES(NM) = 'Some moments (GREEK_11) are NEGATIVE for layer '//C3
!            ACTIONS(NM)  = 'Check Greek moments input - some bad values!'
!          ENDIF
!        ENDDO
!      ENDDO

!  Specialist options, Check the given input
!  =========================================

      IF ( VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_2 ) THEN
        IF ( VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS .GT. NLAYERS-1 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input specialist option 2: NLAYERS_NOMS'
          ACTIONS(NM)  = 'Check NLAYERS_NOMS must be less than NLAYERS'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_3 ) THEN
        IF ( VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF .lt. 2 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input specialist option 3: NLAYERS_CUTOFF'
          ACTIONS(NM)  = 'Check NLAYERS_CUTOFF must be greater than 1'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  TOA contributions
!    TOA Upwelling must be set, if you are using this flag

      IF ( VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS ) THEN
        IF ( .not. VLIDORT_FixIn%Bool%TS_DO_UPWELLING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input TOA contributions, Upwelling Not set'
          ACTIONS(NM)  = 'Must set the DO_UPWELLING flag'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS ) THEN
        N = 0
        DO UTA = 1, VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS
          if ( VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(UTA) .ne. 0.0d0 ) N = N + 1
        ENDDO
        IF ( N .EQ. VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input TOA contributions, Level TOA not set'
          ACTIONS(NM)  = 'Must set one level output to be TOA (0.0)'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!  New section, Checking on Doublet Geometry
!  4/15/20. Add DO_DOUBLET_GEOMETRY flag

      IF ( VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY ) THEN
        IF ( VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input Doublet Geometry: TOA contributions not allowed '
          ACTIONS(NM)  = 'Turn off TOA_CONTRIBS flag'
          STATUS = VLIDORT_SERIOUS
        ENDIF
        IF ( VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input Doublet Geometry: No thermal emission with Doublet Geometry'
          ACTIONS(NM)  = 'Turn off DO_THERMAL_EMISSION flags'
          STATUS = VLIDORT_SERIOUS
        ENDIF
      ENDIF

!mick fix - define NMESSAGES
!  Number of messages

      NMESSAGES = NM

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CHECK_INPUT

!

      SUBROUTINE VLIDORT_CHECK_INPUT_OPTICAL &
          ( NLAYERS, DO_THERMAL_TRANSONLY, DO_LOCAL_MEDIAPROBLEM, LAMBERTIAN_ALBEDO, & ! Input
            DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,              & ! Input
            STATUS, NMESSAGES, MESSAGES, ACTIONS )            ! output

!  Check the optical property inputs

!  4/15/20. Version 2.8.2. MAXMOMENTS_INPUT ==> MAXMOMENTS

!  Module, dimensions and numbers

      USE VLIDORT_pars_m, only : fpk, MAXLAYERS, MAXMOMENTS, MAX_MESSAGES, &
                                 MAXSTOKES_SQ, VLIDORT_SUCCESS, VLIDORT_SERIOUS, &
                                 ZERO, ONE, OMEGA_SMALLNUM

      IMPLICIT NONE

!  Module input
!  ------------

      INTEGER  , intent(in)  :: NLAYERS

!  Thermal Transonly and Local Media problem flags

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY, DO_LOCAL_MEDIAPROBLEM

!  Multilayer optical property (bulk) inputs

      DOUBLE PRECISION, intent(inout) :: OMEGA_TOTAL_INPUT  ( MAXLAYERS )
      DOUBLE PRECISION, intent(in)    :: DELTAU_VERT_INPUT  ( MAXLAYERS )

!  F-matrix GSF expansion coefficients
!   Versions 2.8.1 and Earlier  : Include all that you require for exact single scatter calculations
!   4/15/20, Version 2.8.2 : Only need up to MAXMOMENTS

      DOUBLE PRECISION, intent(in)  :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  Lambertian Surface control

      DOUBLE PRECISION, intent(in)  :: LAMBERTIAN_ALBEDO

!  Module output
!  -------------

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER          :: L, NM
      CHARACTER*3      :: C3

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS

!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Check of VLIDORT Optical Input'
!      ACTIONS(0)      = 'No Action required for this Task'

      NM = NMESSAGES

!  Check Thread-dependent optical property inputs
!  ----------------------------------------------

!  Special clause for transmittance-only calculation
!   4/15/20. Version 2.8.2. Moved to here from CHECK_INPUT

      IF ( DO_THERMAL_TRANSONLY ) OMEGA_TOTAL_INPUT(1:NLAYERS) = ZERO

!  Special clause for Planetary problem. Albedo must be zero
!   4/15/20. Version 2.8.2. Moved to here from MAIN Master

      IF ( LAMBERTIAN_ALBEDO .NE. ZERO .and. DO_LOCAL_MEDIAPROBLEM ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: Lambertian albedo is NON-ZERO for Planetary/Media problem'
        ACTIONS(NM)  = 'Check albedo input, set to zero'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

!  Make sure the Lambertian surface is in range

      IF ( LAMBERTIAN_ALBEDO.LT.ZERO .OR. LAMBERTIAN_ALBEDO.GT.ONE ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: Lambertian albedo not in range [0,1]'
        ACTIONS(NM)  = 'Check albedo input'
        STATUS       = VLIDORT_SERIOUS
      ENDIF

!  Check non-negative optical thickness values

      DO L = 1, NLAYERS
        IF ( DELTAU_VERT_INPUT(L).LE.ZERO ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: optical thickness <= 0, layer '//C3
          ACTIONS(NM)  = 'Check optical thickness input'
          STATUS       = VLIDORT_SERIOUS
        ENDIF
      ENDDO

!  Check single scatter albedos Not Thermal_Transonly

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
        DO L = 1, NLAYERS
          IF ( OMEGA_TOTAL_INPUT(L).GT.ONE-OMEGA_SMALLNUM ) THEN
            WRITE(C3,'(I3)')L ; NM = NM + 1 ; STATUS = VLIDORT_SERIOUS
            MESSAGES(NM) = 'Bad input: SS-albedo too close to 1, layer '//C3
            ACTIONS(NM)  = 'Check SS-albedo input'         
          ENDIF
        ENDDO
      ENDIF

!  solar beam, Omega cannot be too small. Not Thermal_Transonly

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
        DO L = 1, NLAYERS
          IF ( OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
            WRITE(C3,'(I3)')L ; NM = NM + 1 ; STATUS = VLIDORT_SERIOUS
            MESSAGES(NM) = 'Bad input: SS-albedo too close to 0, layer '//C3
            ACTIONS(NM)  = 'Check SS-albedo input'
          ENDIF
        ENDDO
      ENDIF

!  Check first phase function moments

      DO L = 1, NLAYERS
        IF ( GREEKMAT_TOTAL_INPUT(0,L,1).NE.ONE ) THEN
          WRITE(C3,'(I3)')L ; NM = NM + 1 ; STATUS = VLIDORT_SERIOUS
          MESSAGES(NM) = 'First Greek matrix moment not 1 for layer '//C3
          ACTIONS(NM)  = 'Check First Greek matrix moment input'
        ENDIF
      ENDDO

!  Number of messages

      NMESSAGES = NM

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CHECK_INPUT_OPTICAL

!

      SUBROUTINE VLIDORT_DERIVE_INPUT &
      ( VLIDORT_FixIn, VLIDORT_ModIn,      & ! Input type structures
        VLIDORT_Bookkeep, VLIDORT_MSGeom )   ! Output Type structures

!  4/15/20. Version 2.8.2
!  ======================

!   Use Input type structures, and output Bookkeeping structure.

!  PREVIOUS VERSIONS, NOTES
!  ========================

!  This is the basic bookkeeping routine

!  Revised Call, 3/17/17
!mick fix 9/19/2017 - added DO_FOCORR_NADIR & DO_FOCORR_OUTGOING to input
!                   - removed SZA_LOCAL_INPUT & SUN_SZA_COSINES from I/O

!  Miscellaneous and layer input.
!   ( 6/21/10)  Important point: Must use ADJUSTED VZangles as input here.
!   ( 11/19/14) DISAGREE WITH THIS POINT------------------!!!!!!!!!!!!!!!!!
!   revised list, Version 3.8, 3/3/17
!mick fix 3/22/2017 - added DO_FOCORR_NADIR & DO_FOCORR_OUTGOING to input

!  module, dimensions and numbers

      USE VLIDORT_PARS_m, Only : MAXMOMENTS_INPUT, MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXLAYERS,      &
                                 MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_LEVELS, MAX_USER_STREAMS,  &
                                 MAX_PARTLAYERS, MAX_DIRECTIONS, MAX_PSOLS, MAXSTOKES_SQ,             &
                                 ZERO, ONE, HALF, THREE, DEG_TO_RAD, UPIDX, DNIDX

      USE VLIDORT_AUX_m, Only : GETQUAD2, RSSORT

      IMPLICIT NONE

!  Type structure inputs for Derivations
!  -------------------------------------

      TYPE(VLIDORT_Fixed_Inputs), INTENT(IN)       :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT(INOUT) :: VLIDORT_ModIn

!  Type structure bookkeeping and MS Geometry
!  ------------------------------------------

      TYPE(VLIDORT_Bookkeeping), INTENT(INOUT)    :: VLIDORT_Bookkeep   ! mostly output
      TYPE(VLIDORT_Geometry_MS), INTENT(INOUT)    :: VLIDORT_MSGeom     ! some outputs

!  Local variables
!  ---------------

      DOUBLE PRECISION :: DT, RT, MU, H1, DF
      INTEGER ::          I, UT, N, UTA, NSTART
      INTEGER ::          O1, O2, O1_S, UM, IBEAM
      INTEGER ::          NSTREAMS, NSTOKES, NBEAMS, NLAYERS
      INTEGER ::          N_USER_VZANGLES, N_USER_LEVELS, N_USER_RELAZMS, N_USER_OBSGEOMS
      INTEGER ::          N_ALLLAYERS_UP, N_ALLLAYERS_DN, N_CONVTESTS, N_CONV_STREAMS, N_OUT_STREAMS
      INTEGER ::          N_GEOMETRIES, NSOURCES, NVIEWING, NSTOKES_SQ, NSTREAMS_2
      LOGICAL ::          LOCAL_NADIR_ONLY
      INTEGER ::          CONV_START_STREAMS

!  Set derived input
!  =================

!  Proxies

      NSTOKES         = VLidort_FixIn%Cont%TS_NSTOKES
      NBEAMS          = VLidort_ModIn%MSunrays%TS_N_SZANGLES
      NLAYERS         = VLidort_FixIn%Cont%TS_NLAYERS
      NSTREAMS        = VLidort_FixIn%Cont%TS_NSTREAMS
      N_USER_LEVELS   = VLidort_FixIn%UserVal%TS_N_USER_LEVELS
      N_USER_VZANGLES = VLidort_ModIn%MUserVal%TS_N_USER_VZANGLES
      N_USER_RELAZMS  = VLidort_ModIn%MUserVal%TS_N_USER_RELAZMS
      N_USER_OBSGEOMS = VLidort_ModIn%MUserVal%TS_N_USER_OBSGEOMS

!  Automatic input
!mick fix 3/2/2020 - reactived DO_ALL_FOURIER for now

      VLIDORT_Bookkeep%DO_ALL_FOURIER = .FALSE.
!      DO_CLASSICAL_SOLUTION = .TRUE.   ! Removed from the Code, Version 2.8, 7/6/16
!      DO_DIRECT_BEAM        = .TRUE.

!  Define DO_FOCORR_ALONE flag
!  FOCORR_ALONE flag cancels some other flags. Not the DB correction !!!!
!    (Revision, version 2.8, new name for FOCORR_ALONE, used to be DO_SSFULL)

      VLIDORT_Bookkeep%DO_FOCORR_ALONE = &
         ( .NOT.VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE .AND. VLIDORT_ModIn%MBool%TS_DO_FOCORR )

!  Set DB correction flag (this section 11 October 2010)

      VLIDORT_Bookkeep%DO_DBCORRECTION = .false.
      IF ( VLIDORT_ModIn%MBool%TS_DO_FOCORR ) VLIDORT_Bookkeep%DO_DBCORRECTION = .true.

!  Mode of operation
!   SS outgoing sphericity option, added 31 January 2007
!   SS full calculation option added 25 September 2007.
!  Version 2.8. MSMODE VLIDORT is never True now.
!mick mod 9/19/2017 - modified IF structure to accomodate some modification of some
!                     VLIDORT input definitions related to the type of rad solution VLIDORT
!                     returns in different cases.
 !      IF ( DO_FULLRAD_MODE ) THEN
!!        IF ( .NOT. DO_FOCORR_ALONE ) THEN
!!          IF ( DO_FOCORR_NADIR .OR. DO_FOCORR_OUTGOING ) DO_MSMODE_LIDORT = .TRUE.
!!        ENDIF
!      ELSE
!        IF ( .NOT. DO_FOCORR_ALONE ) DO_MSMODE_LIDORT = .TRUE.
!      ENDIF

      VLIDORT_Bookkeep%DO_MSMODE_VLIDORT = .FALSE.
      IF ( VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE ) THEN
        IF ( VLIDORT_ModIn%MBool%TS_DO_FOCORR ) THEN
          !Case: MS trunc + FO corr
          VLIDORT_Bookkeep%DO_MSMODE_VLIDORT = .TRUE.
        ENDIF
      ELSE
        IF ( .NOT.VLIDORT_ModIn%MBool%TS_DO_FOCORR ) THEN
          !Case: MS trunc only
          VLIDORT_Bookkeep%DO_MSMODE_VLIDORT = .TRUE.
        ENDIF
      ENDIF

!  Set thermal MS flag. 1/13/12
!mick fix 3/22/2017 - changed def of DO_MSMODE_THERMAL flag based on implementation
!                     of new FO code.  When DO_FOCORR set, both solar AND THERMAL
!                     direct now come from the FO code.

      !DO_MSMODE_THERMAL = (.not.DO_FULLRAD_MODE) .and. ( DO_SURFACE_EMISSION .and. DO_THERMAL_EMISSION )
      VLIDORT_Bookkeep%DO_MSMODE_THERMAL = VLIDORT_Bookkeep%DO_MSMODE_VLIDORT .and. &
             ( VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION .and. VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION )

!  Directional indices

      IF ( VLIDORT_FixIn%Bool%TS_DO_UPWELLING .AND. VLIDORT_FixIn%Bool%TS_DO_DNWELLING ) THEN
        VLIDORT_Bookkeep%N_DIRECTIONS = 2
        VLIDORT_Bookkeep%WHICH_DIRECTIONS(1) = UPIDX
        VLIDORT_Bookkeep%WHICH_DIRECTIONS(2) = DNIDX
      ELSE
        VLIDORT_Bookkeep%N_DIRECTIONS = 1
        VLIDORT_Bookkeep%WHICH_DIRECTIONS(2) = 0
        IF ( VLidort_FixIn%Bool%TS_DO_UPWELLING ) THEN
          VLIDORT_Bookkeep%WHICH_DIRECTIONS(1) = UPIDX
        ELSE IF ( VLidort_FixIn%Bool%TS_DO_DNWELLING) THEN
          VLIDORT_Bookkeep%WHICH_DIRECTIONS(1) = DNIDX
        ENDIF
      ENDIF

!  New section. Setting the DO_NO_AZIMUTH flag
!    Rt Solutions. 17 January 2006. R. Spurr and V. Natraj.
!  DO_NO_AZIMUTH should be set internally only when:
!     (a) NSTOKES = 1 and nadir view only, and no SS correction
!     (b) DO_MVOUT_ONLY flag is true.
!     (c) SSFULL calculation flag is set

      LOCAL_NADIR_ONLY = .FALSE.
!mick fix 3/30/2015 - added outer if block
      IF ( VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES ) THEN
        IF ( N_USER_VZANGLES.EQ.1 .AND. VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1).EQ.ZERO ) THEN
          LOCAL_NADIR_ONLY = .TRUE.
        ENDIF
      ELSE
        LOCAL_NADIR_ONLY = .TRUE.
      END IF

!  Check for no-azimuth conditions
!mick fix 3/2/2020 - reactivated and made IF condition similar to VLIDORT 2.8.1
      IF ( ( VLIDORT_FixIn%Cont%TS_NSTOKES.EQ.1 .AND. LOCAL_NADIR_ONLY .AND. .NOT.VLIDORT_ModIn%MBool%TS_DO_FOCORR ) .OR. &
             VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY ) THEN
        VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH = .TRUE.
      ENDIF

!  Total quadratures (up and down)

      NSTREAMS_2 = 2*NSTREAMS
      VLIDORT_Bookkeep%NSTREAMS_2 = 2*NSTREAMS

!  Additional quantities (Stokes)

      NSTOKES_SQ     = NSTOKES * NSTOKES
      VLIDORT_Bookkeep%NSTKS_NSTRMS   = NSTOKES * NSTREAMS
      VLIDORT_Bookkeep%NSTKS_NSTRMS_2 = NSTOKES * NSTREAMS_2

!  Size of boundary value problem matrices and vectors

      VLIDORT_Bookkeep%NTOTAL = NLAYERS * VLIDORT_Bookkeep%NSTKS_NSTRMS_2

!  Number of sub and super diagonals in band matrix (boundary value problem)
!mick fix 3/2/2020 - replaced NSTKS_NSTRMS_2 with NSTKS_NSTRMS when defining
!                    N_SUBDIAG & N_SUPDIAG in the NLAYERS > 1 case
 
      IF ( NLAYERS .EQ. 1 ) THEN
        VLIDORT_Bookkeep%N_SUBDIAG = VLIDORT_Bookkeep%NSTKS_NSTRMS_2 - 1
        VLIDORT_Bookkeep%N_SUPDIAG = VLIDORT_Bookkeep%NSTKS_NSTRMS_2 - 1
      ELSE
        VLIDORT_Bookkeep%N_SUBDIAG = 3*VLIDORT_Bookkeep%NSTKS_NSTRMS - 1
        VLIDORT_Bookkeep%N_SUPDIAG = 3*VLIDORT_Bookkeep%NSTKS_NSTRMS - 1
      ENDIF

!  Number of moments. Isotropic option removed 17 January 2006.

      IF ( VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY ) THEN
        VLIDORT_Bookkeep%NMOMENTS = 2
      ELSE
        VLIDORT_Bookkeep%NMOMENTS = MIN ( NSTREAMS_2 - 1, VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT )
      ENDIF

!  Mueller index

      DO O1 = 1, MAXSTOKES
        O1_S = MAXSTOKES*(O1 - 1)
        DO O2 = 1, MAXSTOKES
          VLIDORT_Bookkeep%MUELLER_INDEX(O1,O2) = O1_S + O2
        ENDDO
      ENDDO

!  Set Quadrature abscissae and weights
!    - REVISED for Version 2.8, 3/17/17

      IF ( NSTREAMS .EQ. 1 ) THEN
        VLIDORT_MSGeom%QUAD_STREAMS(1) = SQRT ( ONE / THREE )
        VLIDORT_MSGeom%QUAD_WEIGHTS(1) = ONE
      ELSE
        CALL GETQUAD2 ( ZERO, ONE, NSTREAMS, &
           VLIDORT_MSGeom%QUAD_STREAMS(1:NSTREAMS),VLIDORT_MSGeom%QUAD_WEIGHTS(1:NSTREAMS))
      ENDIF

!  Set auxiliary quantities

      DO I = 1, NSTREAMS
        MU = VLIDORT_MSGeom%QUAD_STREAMS(I)
        VLIDORT_MSGeom%QUAD_STRMWTS(I) = VLIDORT_MSGeom%QUAD_STREAMS(I)*VLIDORT_MSGeom%QUAD_WEIGHTS(I)
        VLIDORT_MSGeom%QUAD_HALFWTS(I) = HALF * VLIDORT_MSGeom%QUAD_WEIGHTS(I)
        VLIDORT_MSGeom%QUAD_ANGLES(I)  = ACOS(MU)/DEG_TO_RAD
        VLIDORT_MSGeom%QUAD_SINES(I)   = SQRT(ONE - MU*MU)
      ENDDO

!  Flux vector set unity input.

      VLIDORT_Bookkeep%FLUXVEC(1) = ONE ; VLIDORT_Bookkeep%FLUXVEC(2:4) = ZERO

!  Convert Surface emission input.............. NOT USED....
!  Input values should be Watts/ sq m
!      IF ( DO_SURFACE_EMISSION ) THEN
!        FP_SURFBB = PI4 * SURFBB
!      ENDIF

!  Greek matrix index

      VLIDORT_Bookkeep%GREEKMAT_INDEX(1) = 1
      VLIDORT_Bookkeep%GREEKMAT_INDEX(2) = 6
      VLIDORT_Bookkeep%GREEKMAT_INDEX(3) = 2
      VLIDORT_Bookkeep%GREEKMAT_INDEX(4) = 11
      VLIDORT_Bookkeep%GREEKMAT_INDEX(5) = 12
      VLIDORT_Bookkeep%GREEKMAT_INDEX(6) = 16

!  Check number of particular solution modes
!   Current default is NPARTICSOLS = 1

      IF ( VLIDORT_Bookkeep%FLUXVEC(3).EQ.ZERO.AND.VLIDORT_Bookkeep%FLUXVEC(4).EQ.ZERO ) THEN
        VLIDORT_Bookkeep%NPARTICSOLS = 1
      ELSE
        VLIDORT_Bookkeep%NPARTICSOLS = 2
      ENDIF

!  D matrices, and multiply by Fluxvector

      VLIDORT_Bookkeep%DFLUX = zero ; VLIDORT_Bookkeep%DMAT = zero ; VLIDORT_Bookkeep%DMAT_PSOLS = zero
      DO O1 = 1, 2
        O2 = O1 + 2
        VLIDORT_Bookkeep%DMAT(O1,O1) = ONE
        VLIDORT_Bookkeep%DMAT(O2,O2) = -ONE
        VLIDORT_Bookkeep%DMAT_PSOLS(O1,O1,1) = ONE
        VLIDORT_Bookkeep%DMAT_PSOLS_FLUXVEC(O1,1) = VLIDORT_Bookkeep%FLUXVEC(O1)
        VLIDORT_Bookkeep%DFLUX(O1) = VLIDORT_Bookkeep%FLUXVEC(O1)
      ENDDO
      IF ( VLIDORT_Bookkeep%NPARTICSOLS .EQ. 2 ) THEN
        DO O2 = 3, 4
          VLIDORT_Bookkeep%DMAT_PSOLS(O2,O2,2) = ONE
          VLIDORT_Bookkeep%DMAT_PSOLS_FLUXVEC(O2,1) = VLIDORT_Bookkeep%FLUXVEC(O2)
          VLIDORT_Bookkeep%DFLUX(O2) = VLIDORT_Bookkeep%FLUXVEC(O2)
        ENDDO
      ENDIF

!  Solar zenith angle cosines/sines

      DO I = 1, NBEAMS
         VLIDORT_MSGeom%COS_SZANGLES(I) = COS( VLIDORT_ModIn%MSunrays%TS_SZANGLES(I) * DEG_TO_RAD )
         VLIDORT_MSGeom%SIN_SZANGLES(I) = SIN( VLIDORT_ModIn%MSunrays%TS_SZANGLES(I) * DEG_TO_RAD )
      ENDDO

!mick fix 9/19/2017 - moved calculation of SUN_SZA_COSINES from here
!                     to LEVELS_GEOMETRY_PREPARE to resolve an I/O contradiction
!  Set average cosines in the refractive geometry case
!      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        DO I = 1, N_SZANGLES
!          MU1 = COS(SZA_LOCAL_INPUT(0,I)*DEG_TO_RAD)
!          DO N = 1, NLAYERS
!            MU2 = COS(SZA_LOCAL_INPUT(N,I)*DEG_TO_RAD)
!            SUN_SZA_COSINES(N,I) = HALF * ( MU1 + MU2 )
!            MU1 = MU2
!          ENDDO
!        ENDDO
!      ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Set performance flags
!  4/15/20. Version 2.8.2
!     This whole section has been moved to vlidort_miscsetups.f90, where it apperas
!     as the subroutine VLIDORT_PERFORMANCE_SETUP. This is in line with new LIDORT practice.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  User stream cosines and secants
!   Adjusted values removed, Version 2.8
!mick fix 3/2/2020 - added defining of VLIDORT_MSGeom%USER_ANGLES

      IF ( VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES ) THEN
        DO I = 1, VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
          VLIDORT_MSGeom%USER_ANGLES(I)  = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(I)
          VLIDORT_MSGeom%USER_STREAMS(I) = COS(DEG_TO_RAD*VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(I))
!          USER_STREAMS(I) = COS(DEG_TO_RAD*USER_VZANGLES_ADJUST(I))
          VLIDORT_MSGeom%USER_SECANTS(I) = ONE / VLIDORT_MSGeom%USER_STREAMS(I)
          VLIDORT_MSGeom%USER_SINES(I)   = SQRT(ONE-VLIDORT_MSGeom%USER_STREAMS(I)**2)
        ENDDO
      ENDIF

!  Number of tests to be applied for convergence
!  Number of streams for convergence
!    - Zenith tolerance no longer applies for the vector code
!    - Set CONV_START_STREAMS = 1. Do not examine for "Zenith_tolerance"
!mick fix 3/2/2020  - defined VLIDORT_Bookkeep%N_OUT_STREAMS also

      N_OUT_STREAMS  = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
      CONV_START_STREAMS = 1
      N_CONV_STREAMS = N_OUT_STREAMS - CONV_START_STREAMS + 1

!  4/15/20. Add Doublet geometry option

      IF ( VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) THEN
        N_CONVTESTS = VLIDORT_Bookkeep%N_DIRECTIONS * N_USER_LEVELS
      ELSE IF ( VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY ) THEN
        N_CONVTESTS = N_CONV_STREAMS * VLIDORT_Bookkeep%N_DIRECTIONS
        N_CONVTESTS = N_CONVTESTS * N_USER_LEVELS
      ELSE
        N_CONVTESTS = N_USER_RELAZMS * N_CONV_STREAMS * VLIDORT_Bookkeep%N_DIRECTIONS
        N_CONVTESTS = N_CONVTESTS * N_USER_LEVELS
      ENDIF

!  Old code, pre Version 2.8.2
!      N_CONVTESTS = VLIDORT_Bookkeep%N_DIRECTIONS * N_USER_LEVELS
!      IF ( .NOT.VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) THEN
!        N_CONVTESTS = N_USER_RELAZMS * N_CONV_STREAMS * N_CONVTESTS
!      ENDIF

      VLIDORT_Bookkeep%N_OUT_STREAMS = N_OUT_STREAMS
      VLIDORT_Bookkeep%N_CONVTESTS   = N_CONVTESTS

!  Number of user azimuth angles
!mick fix 3/2/2020 - defined this new bookkeep variable

      VLIDORT_Bookkeep%LOCAL_N_USERAZM = N_USER_RELAZMS

!  Sort out User vertical level outputs
!  ------------------------------------

!  Sort in ascending order. 
!    - Revised Call, 3/17/17

      IF ( N_USER_LEVELS .GT. 1 ) THEN
        CALL RSSORT(N_USER_LEVELS,VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS))
      ENDIF

!  Mark all output levels not equal to layer boundary values

      NSTART = 0
      UT = 0
      DO UTA = 1, N_USER_LEVELS
        DT = VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(UTA)
        RT = DT - DBLE(INT(DT))
        N = INT(DT) + 1
        IF ( RT.GT.ZERO ) THEN
          UT = UT + 1
          VLIDORT_Bookkeep%PARTLAYERS_OUTFLAG(UTA)  = .TRUE.
          VLIDORT_Bookkeep%PARTLAYERS_OUTINDEX(UTA) = UT
          VLIDORT_Bookkeep%PARTLAYERS_LAYERIDX(UT)  = N
          VLIDORT_Bookkeep%UTAU_LEVEL_MASK_UP(UTA)  = N
          VLIDORT_Bookkeep%UTAU_LEVEL_MASK_DN(UTA)  = N - 1
          VLIDORT_Bookkeep%PARTLAYERS_VALUES(UT)    = RT
        ELSE
          VLIDORT_Bookkeep%PARTLAYERS_OUTFLAG(UTA)  = .FALSE.
          VLIDORT_Bookkeep%PARTLAYERS_OUTINDEX(UTA) = 0
          VLIDORT_Bookkeep%UTAU_LEVEL_MASK_UP(UTA)  = N - 1
          VLIDORT_Bookkeep%UTAU_LEVEL_MASK_DN(UTA)  = N - 1
        ENDIF
      ENDDO
      VLIDORT_Bookkeep%N_PARTLAYERS = UT
      VLIDORT_Bookkeep%DO_PARTLAYERS =  ( VLIDORT_Bookkeep%N_PARTLAYERS .gt. 0 )

!  4/15/20. Version 2.8.2. Partial heights

      do UT = 1, VLIDORT_Bookkeep%N_PARTLAYERS
         N  = VLIDORT_Bookkeep%PARTLAYERS_LAYERIDX(UT)
         RT = VLIDORT_Bookkeep%PARTLAYERS_VALUES(UT)
         H1 = VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(N-1) ; DF = H1 - VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(N) 
         VLIDORT_Bookkeep%PARTLAYERS_HEIGHTS(UT) = H1 - RT * DF
      enddo
      
!  Set masking and number of layer source terms
!  --------------------------------------------
!mick fix 6/25/2012 - initialize STERM_LAYERMASK_UP & STERM_LAYERMASK_DN outside IF blocks
!mick fix 3/2/2020  - make N_ALLLAYERS_UP & N_ALLLAYERS_DN elements of VLIDORT_Bookkeep and define

!   .. for upwelling

      DO N = 1, NLAYERS
        VLIDORT_Bookkeep%STERM_LAYERMASK_UP(N) = .FALSE.
      ENDDO
      IF ( VLIDORT_FixIn%Bool%TS_DO_UPWELLING ) THEN
        UTA = 1 ; UT  = 1
        IF ( .NOT. VLIDORT_Bookkeep%PARTLAYERS_OUTFLAG(UTA) ) THEN
          N_ALLLAYERS_UP = VLIDORT_Bookkeep%UTAU_LEVEL_MASK_UP(UTA) + 1
        ELSE
          N_ALLLAYERS_UP = VLIDORT_Bookkeep%PARTLAYERS_LAYERIDX(UT) + 1
        ENDIF
        DO N = NLAYERS, N_ALLLAYERS_UP, -1
          VLIDORT_Bookkeep%STERM_LAYERMASK_UP(N) = .TRUE.
        ENDDO
        VLIDORT_Bookkeep%N_ALLLAYERS_UP = N_ALLLAYERS_UP
      ENDIF

!   .. for downwelling

      DO N = 1, NLAYERS
        VLIDORT_Bookkeep%STERM_LAYERMASK_DN(N) = .FALSE.
      ENDDO
      IF ( VLIDORT_FixIn%Bool%TS_DO_DNWELLING ) THEN
        UTA = N_USER_LEVELS ; UT  = VLIDORT_Bookkeep%N_PARTLAYERS
        IF ( .NOT. VLIDORT_Bookkeep%PARTLAYERS_OUTFLAG(UTA) ) THEN
          N_ALLLAYERS_DN = VLIDORT_Bookkeep%UTAU_LEVEL_MASK_DN(UTA)
        ELSE
          N_ALLLAYERS_DN = VLIDORT_Bookkeep%PARTLAYERS_LAYERIDX(UT)
        ENDIF
        DO N = 1, N_ALLLAYERS_DN
          VLIDORT_Bookkeep%STERM_LAYERMASK_DN(N) = .TRUE.
        ENDDO
        VLIDORT_Bookkeep%N_ALLLAYERS_DN = N_ALLLAYERS_DN
      ENDIF

!  Number of sources.
!  Save some offsets for indexing geometries (lattice only)
!   --This section revised for the Observational Geometry option
!   ==> 4/15/20. Version 2.8.2. Add Do_Doublet offset.

      NSOURCES = 1 ; IF ( VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES ) NSOURCES = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
      IF ( VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) THEN
        NVIEWING = N_USER_OBSGEOMS
        N_GEOMETRIES = N_USER_OBSGEOMS
      ELSE IF ( VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY ) THEN
        NVIEWING     = N_USER_VZANGLES
        N_GEOMETRIES = NSOURCES * NVIEWING
        DO IBEAM = 1, NBEAMS
          VLIDORT_Bookkeep%SZD_OFFSETS(IBEAM) = NVIEWING * ( IBEAM - 1 )
        END DO
      ELSE
        NVIEWING = N_USER_VZANGLES * VLIDORT_Bookkeep%LOCAL_N_USERAZM
        N_GEOMETRIES = NSOURCES * NVIEWING
        DO IBEAM = 1, NBEAMS
          VLIDORT_Bookkeep%SZA_OFFSETS(IBEAM) = NVIEWING * ( IBEAM - 1 )
          DO UM = 1, N_USER_VZANGLES
            VLIDORT_Bookkeep%VZA_OFFSETS(IBEAM,UM) = &
                VLIDORT_Bookkeep%SZA_OFFSETS(IBEAM) + VLIDORT_Bookkeep%LOCAL_N_USERAZM * (UM - 1)
          END DO
        END DO
      ENDIF
      VLIDORT_Bookkeep%N_GEOMETRIES = N_GEOMETRIES

!  Local post-processing control. Added, 4/9/19

      VLIDORT_Bookkeep%PPSTREAM_MASK = 0 ; VLIDORT_Bookkeep%N_PPSTREAMS = N_USER_VZANGLES
      IF ( VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) VLIDORT_Bookkeep%N_PPSTREAMS = 1
      DO IBEAM = 1, NBEAMS
        IF ( VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) THEN
          VLIDORT_Bookkeep%PPSTREAM_MASK(1,IBEAM) = IBEAM
        else
          do UM = 1, VLIDORT_Bookkeep%N_PPSTREAMS
            VLIDORT_Bookkeep%PPSTREAM_MASK(UM,IBEAM) = UM
          enddo
        endif
     enddo

!  Define TauGrid_Input. Moved from here, 4/15/20 Version 2.8.2
!  --------------------

!      TAUGRID_INPUT(0) = ZERO
!      DO N = 1, NLAYERS
!        TAUGRID_INPUT(N) = TAUGRID_INPUT(N-1) + DELTAU_VERT_INPUT(N)
!      END DO

!  Finish

      RETURN
END SUBROUTINE VLIDORT_DERIVE_INPUT

!

SUBROUTINE VLIDORT_CHAPMAN_INPUT &
      ( VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Bookkeep, & ! Input type structures
        VLIDORT_MSGeom, FAIL, MESSAGE, TRACE )            ! Output TS and exception handling

!  4/15/20. Version 2.8.2
!  ======================

!   Use Input type structures, and output Bookkeeping structure.
!   Formerly "vlidort_geometry.f90" module, restructured and condensed to 1 subroutine

!  The following options apply:

!   1. If the plane-parallel flag is on, no further inputs are required
!   2. If the plane-parallel flag is off, then Pseudo-spherical:
!       (a) Straight line geometry, must specify
!               Earth_radius, height grid
!       (b) Refractive geometry, must specify
!               Earth_radius, height grid
!               pressure grid, temperature grid

!  The logic will be checked before the module is called.

!  Newly programmed by R. Spurr, RT SOLUTIONS Inc. 5/5/05.

!    Based round a call to a pure geometry module which returns slant
!    path distances which was adapted for use in the Radiant model by
!    R. Spurr during an OCO L2 intensive April 24-29, 2005.

!  Major upgrade to include partial Chapman factors. 1/9/18 for Version 2.8
!    -- Cleaned up code
!    -- Introduced New subroutine for partials.

!  module, dimensions and numbers

      USE VLIDORT_PARS_m, only : MAX_SZANGLES, MAXBEAMS, MAXFINELAYERS, MAXLAYERS, MAX_PARTLAYERS, &
                                 ZERO, HALF, ONE, DEG_TO_RAD

!  implicit none

      IMPLICIT NONE

!  Type structure inputs for Derivations
!  -------------------------------------

      TYPE(VLIDORT_Fixed_Inputs), INTENT(IN)       :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT(INOUT) :: VLIDORT_ModIn

!  Type structure bookkeeping and MS Geometry
!  ------------------------------------------

      TYPE(VLIDORT_Bookkeeping), INTENT(INOUT)    :: VLIDORT_Bookkeep   ! inputs
      TYPE(VLIDORT_Geometry_MS), INTENT(INOUT)    :: VLIDORT_MSGeom     ! some outputs

!  Exception handling. 

      LOGICAL      , intent(out)  :: FAIL
      CHARACTER*(*), intent(out)  :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  Proxies

      INTEGER      :: N_PARTLAYERS, NLAYERS, NBEAMS
      LOGICAL      :: DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY

!  fine layer gridding for refraction

      REAL(fpk)    :: ZRFINE(MAXLAYERS,0:MAXFINELAYERS)
      REAL(fpk)    :: PRFINE(MAXLAYERS,0:MAXFINELAYERS)
      REAL(fpk)    :: TRFINE(MAXLAYERS,0:MAXFINELAYERS)

!  number of iterations (refractive case only)
!      This is debug output

      INTEGER      :: ITERSAVE(MAXLAYERS)
!      INTEGER      :: PITERSAVE(MAX_PARTLAYERS)   ! NOT YET ENABLED

!  local height arrays

      REAL(fpk)    :: RADII(0:MAXLAYERS)
      REAL(fpk)    :: DELZ (MAXLAYERS)

!  Standard temperature (K) and pressure (mbar).
!  Loschmidt's number (particles/cm2/km).

      REAL(fpk), PARAMETER  :: T_STANDARD = 273.16D0
      REAL(fpk), PARAMETER  :: P_STANDARD = 1013.25D0
      REAL(fpk), PARAMETER  :: STP_RATIO = T_STANDARD / P_STANDARD
      REAL(fpk), PARAMETER  :: RHO_STANDARD = 2.68675D+24

!  help variables

      INTEGER      :: N, J, NRFINE, K, ITER, IB, UT
      LOGICAL      :: LOOP
      REAL(fpk)    :: MU_NEXT, Z1, Z0, Z, T1, T0, T, P1, P0, Q1, Q0, Q, FU, FL, HP, RP
      REAL(fpk)    :: SZA_GEOM_TRUE, TH_TOA, MU_TOA, SN_TOA, SECMU0
      REAL(fpk)    :: LAYER_DIST, MU_PREV, STH1, SINTH1, STH2, SINTH2, LOCAL_SUBTHICK, &
                      PHI, PHI_0, PHI_CUM, SINPHI, DELPHI, REFRAC, RATIO,              &
                      RE_LOWER, RE_UPPER, DIST, STH2D, SINTH2D, SNELL, MU1, MU2
      CHARACTER*3  :: C3

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  INITIAL SECTION, RADII, REGRIDDING, INITIALIZATION
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Local Proxies

      NLAYERS = VLIDORT_FixIn%Cont%TS_NLAYERS
      NBEAMS  = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
      DO_PLANE_PARALLEL      = VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL
      DO_REFRACTIVE_GEOMETRY = VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY
      N_PARTLAYERS           = VLIDORT_Bookkeep%N_PARTLAYERS

!  earth radii and heights differences

      RADII(0) = VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(0) + VLIDORT_ModIn%MCont%TS_EARTH_RADIUS
      DO N = 1, NLAYERS
        RADII(N) = VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(N) + VLIDORT_ModIn%MCont%TS_EARTH_RADIUS
        DELZ(N)  = RADII(N-1) - RADII(N)
      ENDDO

!  initialise exception handling

      FAIL     = .FALSE.
      MESSAGE  = ' '
      ITERSAVE = 0

!  derive the fine values

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        Z0 = VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(0)      ; P0 = VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID(0) 
        T0 = VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(0) ; Q0 = LOG(P0)
        DO N = 1, NLAYERS
          NRFINE = VLIDORT_FixIn%Chapman%TS_FINEGRID(N)
          LOCAL_SUBTHICK = DELZ(N) / DBLE(NRFINE)
          P1 = VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID(N)    ; Z1 = VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(N) 
          T1 = VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(N) ; Q1 = LOG(P1)
          ZRFINE(N,0) = Z0 ; PRFINE(N,0) = P0 ; TRFINE(N,0) = T0
          DO J = 1, NRFINE - 1
            Z  = Z0 - DBLE(J)*LOCAL_SUBTHICK
            FL = ( Z0 - Z ) / DELZ(N) ; FU = ONE - FL
            Q  = FL * Q1 + FU * Q0    ; T  = FL * T0 + FU * T1
            PRFINE(N,J) = EXP (Q) ; TRFINE(N,J) = T ; ZRFINE(N,J) = Z
          ENDDO
          PRFINE(N,NRFINE) = P1 ; TRFINE(N,NRFINE) = T1 ; ZRFINE(N,NRFINE) = Z1
          Z0 = Z1 ; P0 = P1 ; T0 = T1 ; Q0 = Q1
        ENDDO
      ENDIF

!  @@@@@@@@@@@@@@@@@@@@@@@@@
!  MAIN SECTION, RAY-TRACING
!  @@@@@@@@@@@@@@@@@@@@@@@@@

!  start beam loop

      DO IB = 1, NBEAMS

!  BOA and nadir-TOA values

        SZA_GEOM_TRUE = VLIDORT_ModIn%MSunrays%TS_SZANGLES(IB)
        TH_TOA  = SZA_GEOM_TRUE * DEG_TO_RAD
        MU_TOA = COS(TH_TOA)
        SN_TOA = SQRT ( ONE - MU_TOA * MU_TOA )

!  PART 1. WHOLE-LAYER CHAPMAN FACTORS
!  ===================================

!  Initializing as you need

        VLIDORT_MSGeom%SUNLAYER_COSINES(1:NLAYERS,IB) = zero
        VLIDORT_MSGeom%SZA_LEVEL_OUTPUT(0:NLAYERS,IB) = zero
        DO N = 1, NLAYERS
          VLIDORT_MSGeom%CHAPMAN_FACTORS(N,1:N,IB) = Zero
        ENDDO

!  TOA value

        VLIDORT_MSGeom%SZA_LEVEL_OUTPUT(0,IB) = SZA_GEOM_TRUE

!  plane-parallel case. 4/28/19. Add the sunlayers calculation
      
        IF ( DO_PLANE_PARALLEL ) THEN
          SECMU0 = ONE / MU_TOA
          DO N = 1, NLAYERS
            VLIDORT_MSGeom%SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE
            VLIDORT_MSGeom%CHAPMAN_FACTORS(N,1:N,IB) = SECMU0
          ENDDO
          VLIDORT_MSGeom%SUNLAYER_COSINES(1:NLAYERS,IB) = COS(SZA_GEOM_TRUE*DEG_TO_RAD)
        ENDIF
      
!  Refractive Geometry case
!  ------------------------

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN

!  Rob fix 9/9/14. Sun zero case

          IF ( SN_TOA .eq. zero ) then
            DO N = 1, NLAYERS
              VLIDORT_MSGeom%CHAPMAN_FACTORS(N,1:N,IB) = one
            ENDDO
          ELSE

!  Start value of SZA cosine

            MU_PREV = MU_TOA ; STH2D  = Zero
            DO N = 1, NLAYERS

!  start values

              SINTH1 = SN_TOA * RADII(N) / RADII(0)
              STH1   = ASIN(SINTH1)
              PHI_0  = TH_TOA - STH1
              NRFINE = VLIDORT_FixIn%Chapman%TS_FINEGRID(N)

!  iteration loop

              ITER = 0 ;  LOOP = .TRUE.
              DO WHILE (LOOP.AND.ITER.LT.100)
                ITER = ITER + 1
                PHI_CUM = Zero
                RE_UPPER = ZRFINE(1,0) + RADII(NLAYERS)
                RATIO  = PRFINE(1,0) * STP_RATIO / TRFINE(1,0)
                REFRAC = One + VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER * RATIO
                SNELL = REFRAC * RE_UPPER * SINTH1
                DO K = 1, N
                  LAYER_DIST = Zero
                  LOCAL_SUBTHICK = DELZ(K) / DBLE(NRFINE)
                  DO J = 0, NRFINE - 1
                    RATIO  = PRFINE(K,J) * STP_RATIO / TRFINE(K,J)
                    REFRAC = One + VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER * RATIO
                    RE_LOWER = RE_UPPER - LOCAL_SUBTHICK
                    SINTH2 = SNELL/ (REFRAC * RE_UPPER )
                    IF ( SINTH2.GT.One ) SINTH2 = One
                    STH2 = ASIN(SINTH2)
                    SINTH2D = RE_UPPER * SINTH2 / RE_LOWER
                    IF ( SINTH2D .GT. One) THEN
                      C3 = '000' ; write(C3,'(I3)')IB
                      MESSAGE = 'refraction yields angles > 90 some levels'
                      TRACE = 'Level Chapman calculation failed for Beam # '//C3
                      FAIL = .TRUE. ; RETURN
                    ENDIF
                    STH2D = ASIN(SINTH2D)
                    PHI = STH2D - STH2
                    SINPHI = SIN(PHI)
                    PHI_CUM = PHI_CUM + PHI
                    DIST = RE_UPPER * SINPHI / SINTH2D
                    LAYER_DIST = LAYER_DIST +  DIST
                    RE_UPPER = RE_LOWER
                  ENDDO
                  VLIDORT_MSGeom%CHAPMAN_FACTORS(N,K,IB) = LAYER_DIST / DELZ(K)
                ENDDO

!  examine convergence

                DELPHI = PHI_0 - PHI_CUM
                LOOP = (ABS(DELPHI/PHI_CUM).GT.0.0001d0)

!  Fudge factors to speed up the iteration

                IF ( SZA_GEOM_TRUE .GT. 88.7d0 ) THEN
                  STH1 = STH1 + 0.1d0 * DELPHI ; PHI_0 = TH_TOA - STH1
                ELSE IF ( SZA_GEOM_TRUE .LT. 80.0d0 ) THEN
                  PHI_0 = PHI_CUM ; STH1 = TH_TOA - PHI_0
                ELSE
                  STH1 = STH1 + 0.3d0 * DELPHI ; PHI_0 = TH_TOA - STH1
                ENDIF
                SINTH1 = SIN(STH1)

!  End Iteration loop

              ENDDO

!  failure

              IF ( LOOP ) THEN
                C3 = '000' ; write(C3,'(I3)')IB
                MESSAGE = 'refractive iteration not converged'
                TRACE = 'Level Chapman calculation failed for Beam # '//C3
                FAIL = .TRUE. ; RETURN
              ENDIF

!  Update and save angle output

              MU_NEXT = COS(STH2D)
              MU_PREV = MU_NEXT
              VLIDORT_MSGeom%SZA_LEVEL_OUTPUT(N,IB) = ACOS(MU_NEXT) / DEG_TO_RAD
              ITERSAVE(N) = ITER

!  End layer loop

            ENDDO

!mick fix 1/19/2018 - moved calculation of SUNLAYER_COSINES from LIDORT_DERIVE_INPUT
!                     to here to resolve an I/O contradiction
!  Set average cosines in the refractive geometry case

            MU1 = COS(VLIDORT_MSGeom%SZA_LEVEL_OUTPUT(0,IB)*DEG_TO_RAD)
            DO N = 1, NLAYERS
              MU2 = COS(VLIDORT_MSGeom%SZA_LEVEL_OUTPUT(N,IB)*DEG_TO_RAD)
              VLIDORT_MSGeom%SUNLAYER_COSINES(N,IB) = HALF * ( MU1 + MU2 )
              MU1 = MU2
            ENDDO

!  End general case SZA

          ENDIF

!  Straight line geometry
!  ----------------------

        ELSE
 
!  Rob fix 9/9/14. Sun zero case

          IF ( SN_TOA .eq. zero ) then
            DO N = 1, NLAYERS
              VLIDORT_MSGeom%CHAPMAN_FACTORS(N,1:N,IB) = one
            ENDDO
          ELSE

!  Non-zero case....(resume code)

            DO N = 1, NLAYERS

!  start values

              STH2D  = Zero
              SINTH1 = SN_TOA * RADII(N) / RADII(0)
              STH1   = ASIN(SINTH1)
              RE_UPPER = RADII(0)

!  solar zenith angles are all the same = input value

              VLIDORT_MSGeom%SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE

! loop over layers K from 1 to layer N
!  sine-rule; PHI = earth-centered angle

              DO K = 1, N
                RE_LOWER = RE_UPPER - DELZ(K)
                SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
                STH2   = ASIN(SINTH2)
                PHI    = STH2 - STH1
                SINPHI = SIN(PHI)
                DIST = RE_UPPER * SINPHI / SINTH2
                VLIDORT_MSGeom%CHAPMAN_FACTORS(N,K,IB) = DIST / DELZ(K)
                RE_UPPER = RE_LOWER
                SINTH1 = SINTH2
                STH1   = STH2
              ENDDO

!  finish main layer loop and general SZA clause

            ENDDO
          ENDIF

!mick fix 1/19/2018 - moved calculation of SUNLAYER_COSINES from LIDORT_DERIVE_INPUT
!                     to here to resolve an I/O contradiction
!  Set cosines in the non-refractive geometry case (same all layers)

          VLIDORT_MSGeom%SUNLAYER_COSINES(1:NLAYERS,IB) = COS(SZA_GEOM_TRUE*DEG_TO_RAD)

!  Finish refractive vs. Straightline

        ENDIF

!  PART 2. PARTIAL-LAYER CHAPMAN FACTORS
!  =====================================

!  New, 1/9/18 for Version 2.8, Partial-layer Chapman Factors
!  Add this line when you get the refractive geometry working
!            RFINDEX_PARAMETER, FINEGRID, HEIGHTS, PRESSURES, TEMPERATURES,                    & ! Input

        IF ( N_PARTLAYERS .gt. 0 ) then

!  initialise output

          DO UT = 1, N_PARTLAYERS
            VLIDORT_MSGeom%PARTIAL_CHAPFACS(UT,1:NLAYERS,IB) = Zero
          ENDDO

!  Temporary Fix - No refractive Geometry

          IF ( DO_REFRACTIVE_GEOMETRY ) THEN
            C3 = '000' ; write(C3,'(I3)')IB
            MESSAGE = 'Partial Chapman Factors not enabled yet for refractive geometry'
            TRACE   = 'Derive Chapman subroutine failed for Beam # '//C3
            FAIL = .TRUE. ; RETURN
          ENDIF

!  plane-parallel case

          IF ( DO_PLANE_PARALLEL ) THEN
            SECMU0 = ONE / MU_TOA
            DO UT = 1, N_PARTLAYERS
              N = VLIDORT_Bookkeep%PARTLAYERS_LAYERIDX(UT)
              VLIDORT_MSGeom%PARTIAL_CHAPFACS(UT,1:N,IB) = SECMU0
            ENDDO
          ELSE IF ( .not. DO_REFRACTIVE_GEOMETRY ) then
      
!  Straight line geometry (refractive not enabled - see above)

! Sun zero case

            IF ( SN_TOA .eq. zero ) then
              DO UT = 1, N_PARTLAYERS
                N = VLIDORT_Bookkeep%PARTLAYERS_LAYERIDX(UT)
                VLIDORT_MSGeom%PARTIAL_CHAPFACS(UT,1:N,IB) = ONE
              ENDDO
            ELSE

              DO UT = 1, N_PARTLAYERS
                N = VLIDORT_Bookkeep%PARTLAYERS_LAYERIDX(UT)
                HP = DELZ(N) * VLIDORT_Bookkeep%PARTLAYERS_VALUES(UT)
                RP = RADII(N-1) - HP

!  start values

                STH2D  = ZERO
                SINTH1 = SN_TOA * RP / RADII(0)
                STH1   = ASIN(SINTH1)
                RE_UPPER = RADII(0)

! loop over layers K from 1 to layer N - 1
!  sine-rule; PHI = earth-centered angle

                DO K = 1, N-1
                  RE_LOWER = RE_UPPER - DELZ(K)
                  SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
                  STH2   = ASIN(SINTH2)
                  PHI    = STH2 - STH1
                  SINPHI = SIN(PHI)
                  DIST = RE_UPPER * SINPHI / SINTH2
                  VLIDORT_MSGeom%PARTIAL_CHAPFACS(UT,K,IB) = DIST / DELZ(K)
                  RE_UPPER = RE_LOWER
                  SINTH1 = SINTH2
                  STH1   = STH2
                ENDDO

!  Partial layer

                RE_LOWER = RP
                SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
                STH2   = ASIN(SINTH2)
                PHI    = STH2 - STH1
                SINPHI = SIN(PHI)
                DIST = RE_UPPER * SINPHI / SINTH2
                VLIDORT_MSGeom%PARTIAL_CHAPFACS(UT,N,IB) = DIST / HP

!  Check on the values, use sine-rule distancing.
!          SINTH1 = SN_TOA * RP / RADII(0)
!          STH1   = ASIN(SINTH1)
!          DIST = DOT_PRODUCT(PARTIAL_CHAPFACS(UT,1:N-1,IB),DELZ(1:N-1)) + HP*PARTIAL_CHAPFACS(UT,N,IB)
!          write(*,*)'Check ',UT,DIST,SIN(TH_TOA-STH1)*RADII(0)/SN_TOA

!  finish main partial layer loop

              ENDDO

!  End general sun case, PP.vs.spherical, partial layer clauses (in that order)

            ENDIF
          ENDIF
        ENDIF

!  End Beam Loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE VLIDORT_CHAPMAN_INPUT

!  finish

END MODULE vlidort_inputs_m

