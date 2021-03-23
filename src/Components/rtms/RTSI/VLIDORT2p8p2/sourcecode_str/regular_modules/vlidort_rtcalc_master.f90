
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
! #            VLIDORT_RTCalc_Master (top-level master)         #
! #            VLIDORT_FOURIER                                  #
! #                                                             #
! ###############################################################

!  Notes for Version 2.8.1
!  -----------------------

!       - Revised Water-leaving ocean-optics
!       - Complete new treatment of water-leaving flux-adjustment
!       - Media-property module introduced.
!       - Planetary problem output (S and T) generated.
!       - TOA/BOA illumination control.

!  Notes for Version 2.8.2
!  -----------------------

!       - Geometry function separated
!       - New FOGeometry type structure input
!       - Use of type structures directly for output.
!       - Limitations on zeroing and interna copying.
!       - Use of dynamic memory allocation (started).
!       - Option to use doublet-geometry post-processing

MODULE VLIDORT_RTCalc_Master_m

!  Module for calling VLIDORT, radiances only

!  Parameter types

   USE VLIDORT_pars_m

!  everything public, except the Fourier Routine

   PUBLIC  :: VLIDORT_RTCalc_Master
   PRIVATE :: VLIDORT_FOURIER

   CONTAINS

   SUBROUTINE VLIDORT_RTCalc_Master ( &
      do_debug_input, nstreams_dm, nmoments_dm, nlayers_dm,  & ! debug write and Proxies for dynamic memory (DM)
      VLIDORT_FixIn,    VLIDORT_ModIn,  VLIDORT_Sup,            & ! INPUTS (possibly modified)
      VLIDORT_Bookkeep, VLIDORT_FOGeom, VLIDORT_MSGeom,         & ! SETUPS (possibly modified)
      VLIDORT_Out )                                               ! OUTPUTS

!  Parameter type

      USE VLIDORT_pars_m

!  I/O Type structures

      USE VLIDORT_Inputs_def_m
      USE VLIDORT_Sup_InOut_def_m
      USE VLIDORT_Outputs_def_m

!  Setups Structures for VLIDORT (Bookkeeping, Geometries)
!    4/15/20. Version 2.8.2. New structure

      USE VLIDORT_Setups_def_m

!  internal Type structure

      USE VLIDORT_Work_def_m

!  Version 2.8: VLIDORT module dependencies
!  ----------------------------------------

!  4/15/20. Version 2.8.2. Inputs Check/Derive and Geometry gone now, done in setup_master
!mick fix 3/2/2020 - re-added use of VLIDORT_CHECK_INPUT_DIMS

      USE vlidort_inputs_m      , only : VLIDORT_CHECK_INPUT_DIMS, VLIDORT_CHECK_INPUT_OPTICAL

!  Misc setups and multipliers

      USE vlidort_miscsetups_m  , only : VLIDORT_MISCSETUPS
      USE VLIDORT_MULTIPLIERS_m , only : EMULT_MASTER, EMULT_MASTER_OBSGEO

!  Thermal

      USE VLIDORT_THERMALSUP_m  , only : THERMAL_SETUP

!  4/15/20. Version 2.8.2. Converge module (2 subroutines) now separated

      USE VLIDORT_CONVERGE_m

!  Pack routines

      USE VLIDORT_PACK_m

!  write routines

      USE VLIDORT_WRITEMODULES_m

!  7/7/16, RT Solutions. Version 2.8, only require the FO interface. CORRECTIONS removed
!      USE VLIDORT_CORRECTIONS         !  removed.

      USE VLIDORT_VFO_INTERFACE_m

!  4/9/19. TRANSFLUX superceded, replaced by internal "Adjusted_Backsub" routine
!  VLIDORT 2.8, 9/25/15. RT Solutions, new "transflux" module (Mark1, Mark2)
!  VLIDORT 2.8, 2/3/16 . RT Solutions, new "transflux" module (Mark 3)
!    USE Vlidort_transflux_Master_m
!      USE vlidort_transflux_m

!   Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination, and flux
!  --- THIS IS SEPARATE from the MEDIA_PROPERTIES installation.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      IMPLICIT NONE

!  4/15/20. Version 2.8.2. Arguments

      logical  , intent(in) :: do_debug_input
      integer  , intent(in) :: nstreams_dm, nmoments_dm, nlayers_dm

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs), INTENT (IN)       :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (INOUT) :: VLIDORT_ModIn

!  VLIDORT supplements structure

      TYPE(VLIDORT_Sup_InOut), INTENT (INOUT)       :: VLIDORT_Sup

!  VLIDORT Setups, Type structures, Bookkeeping, FO and MS Geometry outputs
!  4/15/20. Version 2.8.2. New

      TYPE(VLIDORT_Bookkeeping), INTENT (INOUT)      :: VLIDORT_Bookkeep
      TYPE(VLIDORT_Geometry_FO), INTENT (INOUT)      :: VLIDORT_FOGeom
      TYPE(VLIDORT_Geometry_MS), INTENT (INOUT)      :: VLIDORT_MSGeom

!  VLIDORT output structure
!  4/15/20. Version 2.8.2. Return to Intent(InOut). Output must be pre-zeroed.

      TYPE(VLIDORT_Outputs), INTENT (INOUT)           :: VLIDORT_Out

!  VLIDORT local variables
!  +++++++++++++++++++++++

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

!  Boolean
!  =======

!  Flag for Full Radiance  calculation
!    If not set, just produce diffuse (Multiple scatter) field

      LOGICAL ::    DO_FULLRAD_MODE

!  Surface and thermal emission flags

      LOGICAL ::    DO_THERMAL_EMISSION
      LOGICAL ::    DO_SURFACE_EMISSION

!  Beam particular solution, plane parallel flag
!    - Not normally required; pseudo-spherical if not set

      LOGICAL ::    DO_PLANE_PARALLEL

!  Flag for use of LAMBERTIAN surface
!    - If not set, default to BRDF surface

      LOGICAL ::    DO_LAMBERTIAN_SURFACE

!  Directional control

      LOGICAL ::    DO_UPWELLING
      LOGICAL ::    DO_DNWELLING

!  Specialist choices

      LOGICAL ::            DO_TOA_CONTRIBS
      LOGICAL ::            DO_SPECIALIST_OPTION_1
      LOGICAL ::            DO_SPECIALIST_OPTION_2
      LOGICAL ::            DO_SPECIALIST_OPTION_3

!  New 17 May 2012, surface leaving flags

      LOGICAL ::            DO_SURFACE_LEAVING
      LOGICAL ::            DO_SL_ISOTROPIC

!   Version 2.8: DO_WATER_LEAVING and DO_FLUORESCENCE (added 3/3/17) flags.
!   Version 2.8: Transflux iteration control for WATER_LEAVING case.

      LOGICAL ::            DO_WATER_LEAVING       ! New 10/28/15
      LOGICAL ::            DO_FLUORESCENCE        ! New 10/28/15
      LOGICAL ::            DO_TF_ITERATION        ! New 7/7/16

!  Water-leaving output flag. 3/18/19 for Version 2.8.1
      
      LOGICAL ::            DO_WLADJUSTED_OUTPUT 

!   4/26/19 Added control for the media problem. Version 2.8.1
!     Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!       1 = Isotropic illumination from TOA, 2 = Isotropic illumination from BOA
      
      LOGICAL ::            DO_ALBTRN_MEDIA(2)

!   4/28/19 Added control for the planetary problem. Version 2.8.1
      !     Linked to the Media-problem, requires some flux and tranmsittance output (BOA unit illumination)

      LOGICAL ::            DO_PLANETARY_PROBLEM

!  TOA and BOA Illumination flags. 3/23/19 for Version 2.8.1
!   Airglow and Nighttime-viewing scenarios
      
      LOGICAL ::            DO_TOAFLUX
      LOGICAL ::            DO_BOAFLUX

!  Control
!  =======

!  Order of Taylor series (including terms up to EPS^n). Introduced 2/19/14 for Version 2.7
      
      INTEGER ::            TAYLOR_ORDER

!  Number of Stokes components

      INTEGER :: NSTOKES

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  Number of computational layers

      INTEGER :: NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing single scattering correction )

      INTEGER :: NFINELAYERS

!  Number of thermal coefficients (2 should be the default)

      INTEGER :: N_THERMAL_COEFFS

!  Accuracy for convergence of Fourier series

      DOUBLE PRECISION :: VLIDORT_ACCURACY
      
!  TOA Illumination, Flux value. 3/23/19 for Version 2.8.1
!    Must be solar-flux normalized. Designed for Airglow Studies
      
      DOUBLE PRECISION ::   TOAFLUX
      
!  BOA Illumination, Flux value. 3/23/19 for Version 2.8.1
!    Must be solar-flux normalized. Designed for Nighttime Studies.
      
      DOUBLE PRECISION ::   BOAFLUX

!  Special inputs (RT Solutions use only)

      INTEGER ::            NLAYERS_NOMS
      INTEGER ::            NLAYERS_CUTOFF

!  Water-leaving: Control for iterative calculation of transmittanaces
!    Variables added for Version 2.8, 7/6/16

      INTEGER          ::   TF_MAXITER
      DOUBLE PRECISION ::   TF_CRITERION

!  Flux factor

      DOUBLE PRECISION ::   FLUX_FACTOR

!  Chapman
!  =======

!  height grid

      DOUBLE PRECISION :: HEIGHT_GRID ( 0:MAXLAYERS )

!  4/15/20. VErsion 2.8.2. Only required if you want them printed out for debug

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT
!  Multilayer atmospheric inputs (P and T), fine layer gradations 
!   (only required for the Chapman function calculations, w/wo refractive geometry)
!  Refractive index parameter ( Only required for refractive geometry attenuation of the solar beam)

      DOUBLE PRECISION :: PRESSURE_GRID    ( 0:MAXLAYERS )
      DOUBLE PRECISION :: TEMPERATURE_GRID ( 0:MAXLAYERS )
      INTEGER          :: FINEGRID ( MAXLAYERS )
      DOUBLE PRECISION :: EARTH_RADIUS, RFINDEX_PARAMETER

!  Optical
!  =======

!  Multilayer optical property (bulk) inputs

      DOUBLE PRECISION ::   DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION ::   OMEGA_TOTAL_INPUT  ( MAXLAYERS )

!  Phase matrices GSF expansion coefficients
!        Version 2.8.1 and earlier. Include all that you require for exact single scatter calculations
!  4/15/20. Version 2.8.2. Replace 0:MAXMOMENTS_INPUT by 0:MAXMOMENTS. Only needed for MS calculations

      DOUBLE PRECISION ::   GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  F-matrix optical properties (New for Version 2.8, 7/7/16)

      DOUBLE PRECISION ::   FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 )
      DOUBLE PRECISION ::   FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 )

!  Thermal Black Body functions

      DOUBLE PRECISION ::   THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  Lambertian Surface control

      DOUBLE PRECISION :: ALBEDO

!  Surface Black body input

      DOUBLE PRECISION :: SURFACE_BB_INPUT

!  Atmos wavelength (microns) as a diagnostic

      DOUBLE PRECISION :: ATMOS_WAVELENGTH

!  Write
!  =====

      LOGICAL ::            DO_DEBUG_WRITE
      LOGICAL ::            DO_WRITE_INPUT
      LOGICAL ::            DO_WRITE_SCENARIO
      LOGICAL ::            DO_WRITE_FOURIER
      LOGICAL ::            DO_WRITE_RESULTS

      CHARACTER (LEN=60) :: INPUT_WRITE_FILENAME
      CHARACTER (LEN=60) :: SCENARIO_WRITE_FILENAME
      CHARACTER (LEN=60) :: FOURIER_WRITE_FILENAME
      CHARACTER (LEN=60) :: RESULTS_WRITE_FILENAME

!  --------------------------
!  Standard Inputs - Modified
!  --------------------------

!  Boolean
!  =======

!  FO (First-Order) choices (Now all modified Booleans)
!  ----------------------------------------------------

!  SINGLE-SCATTER and DIRECT-BOUNCE for Solar Sources
!  DIRECT_PLANCKF and DIRECT-SURFBB for Thermal Sources

!  FO choices completely rewritten, Version 2.8, 3/1/17.

!  Flag for Computing the corrected FO solution, using FO code Version 1.5
!     New 2 Jul 2013, when it was originally named DO_FO_CALC.
!     - if not set, then VLIDORT will perform a truncated pseudo-spherical SS calculation
!     - If not set, then all other SS choices are turned off

      LOGICAL     :: DO_FOCORR

!  Flag for Use of Externally-derived FO results
!     - DO_FOCORR must be set first.  renamed 2.8

      LOGICAL     :: DO_FOCORR_EXTERNAL

!   sphericity options.
!     - Solar scattering : FOCORR_NADIR and FOCORR_OUTGOING are mutually exclusive. This is checked.
!     - Solar scattering : If both FOCORR_NADIR and FOCORR_OUTGOING, then PLANE-PARALLEL
!     - Direct Planck    : FOCORR_OUTGOING or PLANE-PARALLEL
!     - DO_FOCORR must be set first.

      LOGICAL    :: DO_FOCORR_NADIR
      LOGICAL    :: DO_FOCORR_OUTGOING

!  Additional SSCORR flags
!  -----------------------

!  Flag for performing SSCORR truncation. Single-scattering
!     - Kept here, but disabled for Version 2.8.
!      LOGICAL     :: DO_SSCORR_TRUNCATION

!  Flag for using Phase function in Single-scatter calculations (instead of Legendre Coefficients)
!     - Introduced for Version 2.8, 3/3/17.  R. Spurr
!     - DO_FOCORR must be set first.
!     - Does not apply to thermal case.
!      LOGICAL     :: DO_SSCORR_USEFMAT 

!  Double convergence test flag

      LOGICAL ::    DO_DOUBLE_CONVTEST      ! May be re-set after Checking

!  Basic top-level solar beam control

      LOGICAL ::    DO_SOLAR_SOURCES

!  Beam particular solution: Flag for using refraction in solar paths
!     - This should NOT normally be set.

      LOGICAL ::    DO_REFRACTIVE_GEOMETRY  ! May be re-set after Checking

!  Beam particular solution: Flag for calculating solar beam paths
!    ( Chapman factors = slant/vertical path-length ratios)
!     - This should normally be set.

      LOGICAL ::    DO_CHAPMAN_FUNCTION     ! May be re-set after Checking

!    - Rayleigh only, if set, make sure that scattering Law is Rayleigh!

      LOGICAL ::    DO_RAYLEIGH_ONLY        ! May be re-set after Checking

!  Flag for Use of Delta-M scaling
!    - Should normally be set
!    - Not required for DO_RAYLEIGH_ONLY 

      LOGICAL ::    DO_DELTAM_SCALING       ! May be re-set after Checking

!  Performance flags
!    -- SOLUTION_SAVING gets rid of unneeded RTE computations
!    -- BVP_TELESCOPING creates reduced Boundary value problems
!    -- These flags should be used with CAUTION
!    -- Best, Rayleigh atmospheres with few contiguous cloud/aerosol layers

      LOGICAL ::          DO_SOLUTION_SAVING
      LOGICAL ::          DO_BVP_TELESCOPING

!  Stream angle flag. Normally required for post-processing solutions
!    ( Exception : DO_MVOUT_ONLY is set, then only want Flux output)

      LOGICAL ::          DO_USER_VZANGLES

!  Mean value control (1). If set --> Flux output AS WELL AS Intensities

      LOGICAL ::          DO_ADDITIONAL_MVOUT     ! May be re-set after Checking

!  Mean value control (2). If set --> only Flux output (No Intensities)
!    - DO_USER_STREAMS should be turned off

      LOGICAL ::          DO_MVOUT_ONLY           ! May be re-set after Checking

!  Flag for use of BRDF surface
!  Transmittance only for thermal mode.

      LOGICAL ::          DO_THERMAL_TRANSONLY

!  Observation-Geometry input control. 10/25/12, Version 3.6

      LOGICAL ::          DO_OBSERVATION_GEOMETRY

!  4/15/20. Version 2.8.2. Add Doublet geometry flag

      LOGICAL ::          DO_DOUBLET_GEOMETRY

!  Additional Control for Externalized input (SLEAVE). Introduced 3/18/19 for Version 2.8.1
      
      LOGICAL ::          DO_EXTERNAL_WLEAVE

!  Control
!  4/15/20. Version 2.8.2. No Longer required
!  Number of GSF expansion moments
!      INTEGER ::    NGREEK_MOMENTS_INPUT

!  Chapman
!  4/15/20. Version 2.8.2. No Longer required. Now Supplied by Geometry module
!  Chapman factors (from pseudo-spherical geometry) !mick fix 1/19/2018 - added PARTIAL_CHAPFACS
!      REAL(fpk) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
!      REAL(fpk) :: PARTIAL_CHAPFACS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  UserValues
!  ==========

!  Number and values of User-defined viewing zenith angles (0 to 90 degrees)

      INTEGER          :: N_USER_VZANGLES
      DOUBLE PRECISION :: USER_VZANGLES ( MAX_USER_VZANGLES )

!  Number and values of of user-defined relative azimuths (degrees) (mandatory for Fourier > 0)

      INTEGER   :: N_USER_RELAZMS
      REAL(fpk) :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  Number and values of User-defined vertical levels for output
!    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
!    E.g. For 4.5, this means half way down the 5th layer
!    E.g. For 0.0, this is output at TOA

      INTEGER          :: N_USER_LEVELS
      DOUBLE PRECISION :: USER_LEVELS ( MAX_USER_LEVELS )

!  User-defined Observation Geometry angle input
!   New variables, 25 October 2012, for Observational Geometry input

!  Observation-Geometry input control. 10/25/12

      INTEGER          :: N_USER_OBSGEOMS
      DOUBLE PRECISION :: USER_OBSGEOMS ( MAX_USER_OBSGEOMS, 3 )

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007
!    -- See special note below

      REAL(fpk) :: GEOMETRY_SPECHEIGHT

!  Number and values of solar zenith angles (degrees) to be processed

      INTEGER          :: N_SZANGLES
      DOUBLE PRECISION :: SZANGLES ( MAX_SZANGLES )

!  -----------------------
!  Standard Supplement I/O
!  -----------------------

!  BRDF Inputs
!  -----------

!  Exact (direct bounce) BRDF

      DOUBLE PRECISION ::   EXACTDB_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of BRDF, in the following order

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

!  4/15/20. Version 2.8.2.  Local Fourier-component dimension (0:MAXMOMENTS) has been dropped
!                                 Now using only the terms you want.

      DOUBLE PRECISION ::   BRDF_F_0 ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   BRDF_F   ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION ::   USER_BRDF_F_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_BRDF_F   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

!  Emissivity

      DOUBLE PRECISION ::   EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION ::   USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  SS and DB I/O
!  -------------

!  Single scatter (SS) and Direct bounce (DB) Intensity Results at all angles and optical depths
!   4/15/20. Version 2.8.2. Not used except for debug output write.

      DOUBLE PRECISION ::  STOKES_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION ::  STOKES_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Contribution function (TOA Upwelling only)
!  SS component of Diffuse Field
!      DOUBLE PRECISION ::  SS_CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  New Surface-Leaving Inputs, 17 May 12
!  -------------------------------------

!  Isotropic Surface leaving term (if flag set)

      DOUBLE PRECISION ::   SLTERM_ISOTROPIC  ( MAXSTOKES, MAXBEAMS )

!  Exact Surface-Leaving term

      DOUBLE PRECISION ::   SLTERM_USERANGLES ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams
!  4/15/20. Version 2.8.2. Remove the Fourier component dimension (0:MAXMOMENTS)

      DOUBLE PRECISION ::   SLTERM_F_0      ( MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_SLTERM_F_0 ( MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  ----------------
!  Standard Outputs
!  ----------------

!  Main
!  ====

!  Fourier-summed Stokes vectors at all angles and optical depths
!  4/15/20. Version 2.8.2. Use type structure variable directly
!      DOUBLE PRECISION :: STOKES ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  Results for mean-value output
!mick mod 9/19/2017 - renamed fluxes for separating diffuse and direct

      DOUBLE PRECISION :: MEANST_DIFFUSE ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FLUX_DIFFUSE   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )

!  Direct Fluxes only

      DOUBLE PRECISION :: DNMEANST_DIRECT ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION :: DNFLUX_DIRECT   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

!  4/26-28/19. Special Media-property output. -- Introduced by R. Spurr.
!     ** Output for User-angles and fluxes, Transmittances for solar beam

      DOUBLE PRECISION :: ALBMED_USER ( MAXSTOKES, MAX_USER_STREAMS ), ALBMED_FLUXES(MAXSTOKES,2)    !  TOA illumination
      DOUBLE PRECISION :: TRNMED_USER ( MAXSTOKES, MAX_USER_STREAMS ), TRNMED_FLUXES(MAXSTOKES,2)    !  BOA illumination
      DOUBLE PRECISION :: TRANSBEAM   ( MAXSTOKES, MAXBEAMS )                                        !  Planetary problem

!  Transmittance/Flux Output   for the Waterleaving case
!    - Complete/Direct Actinic and Regular Fluxes
!    - RT Solutions, 9/25/15. removed 2/3/16 (superceded by Mark 3 output)
!      DOUBLE PRECISION :: MEANST_DIFFUSE_TF  ( MAX_SZANGLES, MAXSTOKES )
!      DOUBLE PRECISION :: FLUX_DIFFUSE_TF    ( MAX_SZANGLES, MAXSTOKES )
!      DOUBLE PRECISION :: DNMEANST_DIRECT_TF ( MAX_SZANGLES, MAXSTOKES )
!      DOUBLE PRECISION :: DNFLUX_DIRECT_TF   ( MAX_SZANGLES, MAXSTOKES )

!  Contribution functions (TOA Upwelling only)
!  Fourier component of Diffuse Field
!      DOUBLE PRECISION ::  MS_CONTRIBS_F ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
!  Fourier-summed values
!      DOUBLE PRECISION ::  CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  Ancillary Output
!  ================

!  Fourier numbers used

      INTEGER      :: FOURIER_SAVED ( MAX_SZANGLES )

!  Exception handling
!  ==================

!  Exception handling for Input Checking. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER              :: STATUS_INPUTCHECK
      INTEGER              :: NCHECKMESSAGES
      CHARACTER(Len=120)   :: CHECKMESSAGES(0:MAX_MESSAGES)
      CHARACTER(Len=120)   :: ACTIONS (0:MAX_MESSAGES)

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER              :: STATUS_CALCULATION
      CHARACTER(Len=120)   :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

!  ============================================================
!  ============================================================

!  SPECIAL NOTE on variable GEOMETRY_SPECHEIGHT

!    This is only required when the outgoing sphericity correction is
!    in operation. Otherwise, the regular pseudo-spherical correction
!    (wiht or without an exact single-scatter correction) uses the same
!    set of angles all the up the nadir from the bottom of atmosphere.

!     This height is normally set equal to the height at the lowest
!     level of the atmosphere: GEOMETRY_SPECHEIGHT = HEIGHT_GRID(NLAYERS)
!     In this case, no adjustment to the geometrical inputs is needed
!     for the outgoing sphericity correction.

!     If there is a situation GEOMETRY_SPECHEIGHT < HEIGHT_GRID(NLAYERS),
!     then an adjustment to the geometrical inputs is needed for the
!     outgoing sphericity correction. This adjustment is internal and
!     the model output will still be given at the geometrical angles
!     as specified by the user, even though these angles may not be the
!     ones at which the calculations were done. This situation will occur
!     when we are given a BOA geometry but we want to make a calculation
!     for a reflecting surface (such as a cloud-top) which is above the
!     BOA level. In this case, GEOMETRY_SPECHEIGHT = 0.0, and the lowest
!     height HEIGHT_GRID(NLAYERS) = cloud-top.

!     This height cannot be greater than HEIGHT_GRID(NLAYERS). If this is
!     the case, this height will be set equal to HEIGHT_GRID(NLAYERS), and
!     the calculation will go through without the adjustment. A warning
!     about this incorrect input choice will be sent to LOGFILE.

!  ============================================================
!  ============================================================

! ######################################################################
!                       Local Arguments
! ######################################################################

!  Local bookkeeping
!  -----------------

!               Intent(In) To the Fourier  routine
!               Intent(In) To the Converge routine

!  Mode flags (Robfix, 13 January 2012. Add MSMODE_THERMAL flag)

      LOGICAL         :: DO_MSMODE_VLIDORT
      LOGICAL         :: DO_MSMODE_THERMAL

!  Flag for Doing FO calculation alone (no Multiple scatter)
!     - Formerly called DO_SSFULL (confusingly!) - internal variable now (Version 2.8)
!     - DO_FOCORR must be set first.

      LOGICAL         :: DO_FOCORR_ALONE

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER         :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS, NSTKS_NSTRMS = NSTOKES * NSTREAMS
!  total number of layers and streams NTOTAL = NSTKS_NSTRMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER         :: NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NSTOKES_SQ
      INTEGER         :: NTOTAL
      INTEGER         :: N_SUBDIAG, N_SUPDIAG

!  Output optical depth masks and indices

      LOGICAL          :: PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER          :: PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER          :: UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER          :: UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!  Off-grid optical depths (values, masks, indices)

      LOGICAL          :: DO_PARTLAYERS
      INTEGER          :: N_PARTLAYERS
      INTEGER          :: PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION :: PARTLAYERS_VALUES   ( MAX_PARTLAYERS )

!  Layer masks for doing integrated source terms

      LOGICAL         :: STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL         :: STERM_LAYERMASK_DN ( MAXLAYERS )
      INTEGER         :: N_ALLLAYERS_UP, N_ALLLAYERS_DN

!  Post-processing masks. Introduced 4/9/19.

      INTEGER         :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Number of directions (1 or 2) and directional array

      INTEGER         :: N_DIRECTIONS, WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Number of convergence tests, Azimuths Geometry

      INTEGER         :: N_CONVTESTS, LOCAL_N_USERAZM, N_GEOMETRIES

!  Solar beam flags (always internal)

      LOGICAL         :: DO_MULTIBEAM ( MAXBEAMS, 0:MAXFOURIER )

!  Real Eigensolution mask

      LOGICAL         :: DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )

!  Offsets for geometry indexing
!    4/15/20. Version 2.8.2. Add Doublet geometry offsets (SZD_OFFSETS) 

      INTEGER         :: VZA_OFFSETS(MAXBEAMS,MAX_USER_STREAMS)
      INTEGER         :: SZD_OFFSETS(MAXBEAMS)

!  Mueller index

      INTEGER         :: MUELLER_INDEX(MAXSTOKES,MAXSTOKES)

!  Greekmat indices, DMAT, FLUXVEC and related quantities (VLIDORT only)

      INTEGER          :: GREEKMAT_INDEX(6)
      DOUBLE PRECISION :: FLUXVEC(MAXSTOKES), DFLUX(MAXSTOKES)
      DOUBLE PRECISION :: DMAT(MAXSTOKES,MAXSTOKES)

!  Number of particular solutions (not enabled yet)

      INTEGER          :: NPARTICSOLS
      DOUBLE PRECISION :: DMAT_PSOLS(MAXSTOKES,MAXSTOKES,2)
      DOUBLE PRECISION :: DMAT_PSOLS_FLUXVEC(MAXSTOKES,2)

!  Local MS-Geometry. 4/15/20. Version 2.8.2. Code reorganized
!  ----------------------------------------------------------------

!  Local solar zenith angles Cosines (regular case)

      DOUBLE PRECISION :: COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION :: SIN_SZANGLES ( MAX_SZANGLES )

!  Quadrature weights and abscissae, and product

      DOUBLE PRECISION :: QUAD_STREAMS (MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_WEIGHTS (MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_STRMWTS (MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_HALFWTS (MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_ANGLES  (MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_SINES   (MAXSTREAMS)

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      DOUBLE PRECISION :: USER_ANGLES   (MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_STREAMS  (MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SECANTS  (MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES   ( MAX_USER_STREAMS )

!  Local input solar zenith angles Cosines
!  ( Only required for refractive geometry attenuation of the solar beam)

      DOUBLE PRECISION :: SUN_SZA_COSINES  ( MAXLAYERS  , MAX_SZANGLES )
!      DOUBLE PRECISION :: SZA_LEVEL_OUTPUT ( 0:MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION :: SUNLAYER_COSINES ( MAXLAYERS  , MAX_SZANGLES )  ! Average cosines for refractive geometry case

!  Chapman Factors (Level and (1/9/18) partials, Local solar zenith angles at nadir, average cosines
!    These will be calculated either from FO Geometry, or stand-alone (Better)

      DOUBLE PRECISION :: CHAPMAN_FACTORS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: PARTIAL_CHAPFACS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Arrays for setups and Corrections
!  ---------------------------------

!               Intent(In) To the Fourier routine

!  Local flags for the solution saving option
!      INTEGER          :: LAYER_MAXMOMENTS ( MAXLAYERS )

!  Initial transmittances * (secants)

      DOUBLE PRECISION :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Saved arrays for truncation factor and Delta-M scaling

      DOUBLE PRECISION :: TRUNC_FACTOR ( MAXLAYERS )
      DOUBLE PRECISION :: FAC1 ( MAXLAYERS )

!  Derived Solar-beam Transmittance at all levels
!mick fix 9/19/2017 - added LEVELS_SOLARTRANS & PARTIALS_SOLARTRANS
!                     to facilitate correction of direct flux

      DOUBLE PRECISION :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )

!  Derived Solar-beam Transmittance at BOA (Diagnostic output). 
!  Rob Fix 10/24/14, 11/17/14.
!    SOLARBEAM_BOATRANS is computed with Unscaled optical depths
!    TRANS_SOLAR_BEAM with scaled ODs. [ Same if  no Deltam-scaling]

      DOUBLE PRECISION :: SOLARBEAM_BOATRANS ( MAXBEAMS )

!  Derived Slant optical thickness inputs
!mick fix 9/19/2017 - added DELTAU_SLANT_UNSCALED & PARTAU_SLANT_UNSCALED
!                     to facilitate correction of direct flux

      DOUBLE PRECISION :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: DELTAU_SLANT_UNSCALED ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: PARTAU_SLANT_UNSCALED ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Scaled SSAs and phase function moments

      DOUBLE PRECISION :: OMEGA_TOTAL    ( MAXLAYERS )
      DOUBLE PRECISION :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  L'Hopital's rule logical variables

      LOGICAL          :: EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Coefficient functions for user-defined angles

      DOUBLE PRECISION :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Fourier component output
!  ------------------------

!               Intent(Out) from the Fourier routine
!               Intent(in)  to   the Converge routine

!  Fourier components for user-defined solutions

      DOUBLE PRECISION :: STOKES_F ( MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Help arrays from the SS/DB correction routines
!  ==============================================

!  THIS SECTION HAS BEEN REMOVED for VERSION 2.8

!  Saved Legendre polynomials
!      !? DOUBLE PRECISION :: SS_PLEG_UP(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
!      !? DOUBLE PRECISION :: SS_PLEG_DN(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
!  Saved TMS (Nakajima-Tanaka) factor
!      DOUBLE PRECISION :: TMS ( MAXLAYERS )
!  Local truncation factors for additional DELTAM scaling
!      DOUBLE PRECISION :: SSFDEL ( MAXLAYERS )
!  Exact Phase function calculations
      !DOUBLE PRECISION :: EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS)
      !DOUBLE PRECISION :: EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS)
!      DOUBLE PRECISION :: ZMAT_UP ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
!      DOUBLE PRECISION :: ZMAT_DN ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
!  Cumulative single scatter source terms
!      DOUBLE PRECISION :: SS_CUMSOURCE_UP ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
!      DOUBLE PRECISION :: SS_CUMSOURCE_DN ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
!  Atmospheric attenuation before reflection
!      DOUBLE PRECISION :: ATTN_DB_SAVE ( MAX_GEOMETRIES )
!  Exact direct beam source terms
!      DOUBLE PRECISION :: EXACTDB_SOURCE ( MAX_GEOMETRIES, MAXSTOKES )
!  Cumulative direct bounce source terms
!      DOUBLE PRECISION :: DB_CUMSOURCE ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
!  Solar beam attenuation to BOA (required for exact DB calculation)
!      DOUBLE PRECISION :: BOA_ATTN ( MAX_GEOMETRIES )
!  Outgoing sphericity stuff - Whole and part-layer LOS transmittance factors
!      DOUBLE PRECISION :: UP_LOSTRANS    ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_LOSTRANS    ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: UP_LOSTRANS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_LOSTRANS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )
!  Whole and part-layer multipliers
!      DOUBLE PRECISION :: UP_MULTIPLIERS    ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_MULTIPLIERS    ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: UP_MULTIPLIERS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_MULTIPLIERS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )

!  Arrays required at the Top level
!  ================================

!               Intent(In) To the Fourier routine

!  Input optical properties after delta-M scaling

      DOUBLE PRECISION :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Local input solar zenith angles by levels
!  ( Only required for refractive geometry attenuation of the solar beam)
!  These will be set internally if the refraction flag is set.

      DOUBLE PRECISION :: SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER          :: BEAM_CUTOFF ( MAXBEAMS )

!  Average-secant and initial tramsittance factors for solar beams.

      DOUBLE PRECISION :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: LOCAL_CSZA     ( 0:MAXLAYERS, MAXBEAMS )

!  Solar beam attenuation

      DOUBLE PRECISION :: TRANS_SOLAR_BEAM ( MAXBEAMS )

!  Reflectance flags

      LOGICAL         :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Local flags for the solution saving option

      LOGICAL         :: DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Local flags,  BVP telescoping enhancement

      LOGICAL         :: BVP_REGULAR_FLAG ( 0:MAXMOMENTS )

!  Masking for regular case. Required again for linearization
      !? INTEGER         :: LCONMASK ( MAXSTREAMS, MAXLAYERS )
      !? INTEGER         :: MCONMASK ( MAXSTREAMS, MAXLAYERS )

!  Telescoping initial flag (modified argument)

      LOGICAL         :: DO_BVTEL_INITIAL
      INTEGER         :: BVTEL_FOURIER

!  Number of telescoped layers, active layers,  Size of BVP matrix 
      !INTEGER         :: NLAYERS_TEL, ACTIVE_LAYERS ( MAXLAYERS ), N_BVTELMATRIX_SIZE

!  Set up for band matrix compression
      !? INTEGER         :: BMAT_ROWMASK    ( MAXTOTAL, MAXTOTAL )
      !? INTEGER         :: BTELMAT_ROWMASK ( MAXTOTAL, MAXTOTAL )

!  Transmittance setups
!  --------------------

!               Intent(In) To the Fourier routine

!  Discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      DOUBLE PRECISION :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_UTUP_DISORDS ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UTDN_DISORDS ( MAXSTREAMS, MAX_PARTLAYERS )

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage

      DOUBLE PRECISION :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      DOUBLE PRECISION :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage 

      DOUBLE PRECISION :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Cumulative transmittance

      DOUBLE PRECISION :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )

!  Multiplier arrays
!  -----------------

!               Intent(In) To the Fourier routine

!  Forcing term multipliers (saved for whole atmosphere)

      DOUBLE PRECISION :: EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )

!  Partial layer multipliers

      DOUBLE PRECISION :: UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  Thermal Setup outputs
!  ---------------------

!  Optical depth powers

      DOUBLE PRECISION :: DELTAU_POWER ( MAXLAYERS,      MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )

!  Thermal coefficients, bookkeeping

      DOUBLE PRECISION :: THERMCOEFFS ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: TCOM1       ( MAXLAYERS, MAX_THERMAL_COEFFS )

!  Tranmsittance solutions

      DOUBLE PRECISION :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION :: T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Vector FO output
!  ----------------
!mick mod 9/19/2017 - added FO_STOKES_ATMOS & FO_STOKES_SURF (new output from vector FO code)

!  from VFO_MASTER_INTERFACE
!    SS/DB (solar), DTA/DTS (thermal) and total

      DOUBLE PRECISION :: FO_STOKES_SS    ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_DB    ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_STOKES_DTA   ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_DTS   ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION :: FO_STOKES_ATMOS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_SURF  ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_STOKES       ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  Transmitted Fluxes for Adjusted Water-leaving
!  ---------------------------------------------

!  4/9/19. Additional output from the SFO Interface, for the sleave correction

      DOUBLE PRECISION :: FO_CUMTRANS ( max_user_levels, MAX_GEOMETRIES )

!  4/9/19. Surface leaving FO assignation (un-adjusted)

      DOUBLE PRECISION :: FO_SLTERM ( MAXSTOKES, MAX_GEOMETRIES )

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  First Introduced 3/22/17 for LIDORT, based on VLIDORT code 

      DOUBLE PRECISION :: TRANS_ATMOS_FINAL ( MAXBEAMS )

!  Other Local variables
!  =====================

!  Local quantities for the water-leaving configuration

      DOUBLE PRECISION :: CUMSOURCE_DB, SLTERM_LOCAL, TRANS, TFACTOR

!  Local flags for media problem and Planetary problem control
      
      LOGICAL   ::        LOCAL_DO_ALBTRN_MEDIA(2)
      LOGICAL   ::        LOCAL_DO_PLANETARY_PROBLEM
      LOGICAL   ::        LOCAL_DO_MEDIA_PROBLEM

      LOGICAL ::          LOCAL_DO_NO_AZIMUTH
      LOGICAL ::          SAVE_DO_NO_AZIMUTH
      LOGICAL ::          LOCAL_ITERATION

!  inclusion flags

      LOGICAL ::          DO_INCLUDE_SURFACE
      LOGICAL ::          DO_INCLUDE_SURFEMISS
      LOGICAL ::          DO_INCLUDE_THERMEMISS

!  Local fourier index

      INTEGER ::          FOURIER, N_FOURIERS

!  Misc. help (4/15/20, Version 2.8.2, some additions, removed IUNIT)

      INTEGER ::          OFF, O1, UA, UM, IB, G, G1, G2, K, L, V, UT, UTA, LUM, LUA, NMOMS, N
      INTEGER ::          TESTCONV, STATUS_SUB
      INTEGER ::          SUNIT, RUNIT ! IUNIT

!  For Convergence routines

      DOUBLE PRECISION :: AZM_ARGUMENT, DFC
      DOUBLE PRECISION :: AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )

      INTEGER ::          IBEAM_COUNT, IBEAM, NSOURCES
      LOGICAL ::          BEAM_ITERATION ( MAXBEAMS )
      INTEGER ::          BEAM_TESTCONV  ( MAXBEAMS )

!  Adjusted geometries. New, 2007.
!  -------------------------------

!  No longer required for Version 2.8, Have removed the Correction routines
!               Intent(Out) from the Adjust-geometry routine
!               Intent(In)  to   the Correction      routine
!      DOUBLE PRECISION :: USER_VZANGLES_ADJUST ( MAX_USER_VZANGLES )
!      DOUBLE PRECISION :: SZANGLES_ADJUST      ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
!      DOUBLE PRECISION :: USER_RELAZMS_ADJUST  ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
!      LOGICAL ::          ADJUST_SURFACE
!      DOUBLE PRECISION :: MODIFIED_ERADIUS

!  Fluxes for Water-leaving, Master routine (Mark 3)
!  -------------------------------------------------
!    - Introduced 1/12/16, Validated 2/3/16. DISABLED Version 2.8.1.
!      DOUBLE PRECISION :: FLUX_DIFFUSE_FINAL  ( MAX_SZANGLES )

!  Helper variables
!  ----------------

!  flags

      LOGICAL ::          DO_NO_AZIMUTH
      LOGICAL ::          DO_ALL_FOURIER
!      LOGICAL ::          DO_CLASSICAL_SOLUTION  !!! removed for Version 2.8
      LOGICAL ::          DO_DBCORRECTION

!  Single scatter flux multipier

      DOUBLE PRECISION :: SS_FLUX_MULTIPLIER

!  output control

      INTEGER ::          LOCAL_UM_START
      INTEGER ::          N_OUT_STREAMS
      DOUBLE PRECISION :: OUT_ANGLES ( MAX_USER_STREAMS )

! #################### INNOVATIONS 5/5/20 ##############
!  4/29/20. SPECIALIST OUTPUT. Rayleigh Fourier output (TOA Upwelling only)
!      LOGICAL ::          SPECIAL_RAYF_OUTPUT
! #################### INNOVATIONS 5/5/20 ##############

!  Local optical depths

      DOUBLE PRECISION :: TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION :: TAUGRID ( 0:MAXLAYERS )

!  Contribution functions (TOA Upwelling only)
!  -------------------------------------------

!  Fourier component of Diffuse Field

      DOUBLE PRECISION ::  MS_CONTRIBS_F ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )

!  4/15/20/ Version 2.8.2. Use type structure  variables directly.
!  Fourier-summed values
!      DOUBLE PRECISION ::  CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
!  Single scatter
!      DOUBLE PRECISION :: SS_CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  Work type structures
!  --------------------

      TYPE(VLIDORT_Work_Miscellanous) :: Misc
      TYPE(VLIDORT_Work_Thermal)      :: Therm
      TYPE(VLIDORT_Work_Multiplier)   :: Mult
      !TYPE(VLIDORT_Work_Corrections)  :: Corr
      !TYPE(VLIDORT_Work_FirstOrder)   :: Fo

!  Intermediate quantities

      DOUBLE PRECISION :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION :: PIMM_KM ( MAX_ALLSTRMS_P1 )

!  Local error handling
!      LOGICAL ::          FAIL
!  Test variables
!      LOGICAL ::          DO_FDTEST =.FALSE.

!  This is the flag for a complete write-up of the inputs
!  4/15/20/ Version 2.8.2. Now an input argument.
!      LOGICAL ::          DO_DEBUG_INPUT=.FALSE.
!      LOGICAL ::          DO_DEBUG_INPUT=.TRUE.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!                S T A R T   T H E    C O D E
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialize some variables
!  -------------------------

!  Main status

      STATUS_CALCULATION = VLIDORT_SUCCESS
      STATUS_INPUTCHECK  = VLIDORT_SUCCESS

!mick fix 6/29/11 - initialize "Input checks"
!  Input checks

      NCHECKMESSAGES = 0
      CHECKMESSAGES  = ' '
      ACTIONS        = ' '

!  Model calculation

      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      TRACE_3 = ' '

!  Local user indices

      LUM = 1
      LUA = 1

!  Check input dimensions. New for Version 2.7
!  -------------------------------------------

      CALL VLIDORT_CHECK_INPUT_DIMS &
      ( VLIDORT_FixIn, VLIDORT_ModIn, &
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VLIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        VLIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        VLIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  Fixed Boolean inputs
!  --------------------

      DO_FULLRAD_MODE        = VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE

      DO_THERMAL_EMISSION    = VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION
      DO_SURFACE_EMISSION    = VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION
      DO_PLANE_PARALLEL      = VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL

      DO_UPWELLING           = VLIDORT_FixIn%Bool%TS_DO_UPWELLING
      DO_DNWELLING           = VLIDORT_FixIn%Bool%TS_DO_DNWELLING

!      DO_QUAD_OUTPUT         = VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT  ! removed 7/7/16
      DO_TOA_CONTRIBS        = VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS

      DO_LAMBERTIAN_SURFACE  = VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE
      DO_SPECIALIST_OPTION_1 = VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_1
      DO_SPECIALIST_OPTION_2 = VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_2
      DO_SPECIALIST_OPTION_3 = VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_3

!  New 17 May 2012. Surface leaving

      DO_SURFACE_LEAVING     = VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC        = VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC

!  New Version 2.8. Water leaving flags.

      DO_WATER_LEAVING       = VLIDORT_FixIn%Bool%TS_DO_WATER_LEAVING   ! introduced 10/28/15
      DO_FLUORESCENCE        = VLIDORT_FixIn%Bool%TS_DO_FLUORESCENCE    ! introduced 10/28/15
      DO_TF_ITERATION        = VLIDORT_FixIn%Bool%TS_DO_TF_ITERATION    ! introduced 07/07/16
      
!  3/18/19. New for Version 2.8.1. Water leaving output flag.
      
      DO_WLADJUSTED_OUTPUT   = VLIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT

!  4/26/19. new for Version 2.8.1, introduced by R. Spurr
!    -- Control for Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA

      DO_ALBTRN_MEDIA        = VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA
      
!  4/28/19. new for Version 2.8.1, introduced by R. Spurr
!    -- Control for Planetary problem.

      DO_PLANETARY_PROBLEM = VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM
      
!  TOA/BOA Illumination flags. 3/23/19 for Version 2.8.1
      
      DO_TOAFLUX    = VLIDORT_FixIn%Bool%TS_DO_TOA_ILLUMINATION
      DO_BOAFLUX    = VLIDORT_FixIn%Bool%TS_DO_BOA_ILLUMINATION

!  Fixed Control inputs
!  --------------------

!  Taylor parameter new, Version 2p7

      TAYLOR_ORDER     = VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER 

!  Numbers

      NSTOKES          = VLIDORT_FixIn%Cont%TS_NSTOKES
      NSTREAMS         = VLIDORT_FixIn%Cont%TS_NSTREAMS
      NMOMS            = 2*NSTREAMS - 1
      NLAYERS          = VLIDORT_FixIn%Cont%TS_NLAYERS
      NFINELAYERS      = VLIDORT_FixIn%Cont%TS_NFINELAYERS
      N_THERMAL_COEFFS = VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS

      VLIDORT_ACCURACY = VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY

      NLAYERS_NOMS     = VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS
      NLAYERS_CUTOFF   = VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF

!  New for Version 2.8, Water-leaving iteration control. 7/7/16

      TF_MAXITER       = VLIDORT_FixIn%Cont%TS_TF_MAXITER
      TF_CRITERION     = VLIDORT_FixIn%Cont%TS_TF_CRITERION

!  TOA/BOA Illumination. 3/23/19 for Version 2.8.1
      
      TOAFLUX       = VLIDORT_FixIn%Cont%TS_TOA_ILLUMINATION
      BOAFLUX       = VLIDORT_FixIn%Cont%TS_BOA_ILLUMINATION

!  Fixed Beam inputs
!  -----------------

      FLUX_FACTOR      = VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR

!  Fixed User Value inputs
!  -----------------------

      N_USER_LEVELS    = VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS

!  Fixed Chapman Function inputs
!  -----------------------------

!  height grid

      HEIGHT_GRID(0:NLAYERS)      = VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(0:NLAYERS)

!  4/15/20. Version 2.8.2. None of these required here, but may be required for debug input
!      PRESSURE_GRID(0:NLAYERS)    = VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID(0:NLAYERS)
!      TEMPERATURE_GRID(0:NLAYERS) = VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(0:NLAYERS)
!      FINEGRID(1:NLAYERS)         = VLIDORT_FixIn%Chapman%TS_FINEGRID(1:NLAYERS)
!      RFINDEX_PARAMETER           = VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER

!  Fixed Optical inputs
!  --------------------

!  4/15/20. Version 2.8.2. No longer need NMOMENTS_INPUT, replace by NSTREAMS_2 = 2 * NSTREAMS
!  4/15/20. Version 2.8.2. OMEGA_TOTAL_INPUT now Fixed-Optical, move here.

      DELTAU_VERT_INPUT(1:NLAYERS) = VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(1:NLAYERS)
      OMEGA_TOTAL_INPUT(1:NLAYERS) = VLIDORT_FixIn%Optical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS)
      THERMAL_BB_INPUT(0:NLAYERS)  = VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT(0:NLAYERS)

      ALBEDO           = VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO
      SURFACE_BB_INPUT = VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT
      ATMOS_WAVELENGTH = VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH

!  New version 2.8, F-matrix input. 7/7/16
!  4/15/20. Version 2.8.2. Local copying of F-matrices and coefficients.
!mick fix 9/19/2017 - swapped layer & geo indices. Local is same as input now,

      NSTREAMS_2   = VLIDORT_Bookkeep%NSTREAMS_2
      N_GEOMETRIES = VLIDORT_Bookkeep%N_GEOMETRIES
      IF ( VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES ) THEN
         DO N = 1, NLAYERS
            FMATRIX_UP(N,1:N_GEOMETRIES,:) = VLIDORT_FixIn%Optical%TS_FMATRIX_UP(N,1:N_GEOMETRIES,:)
            FMATRIX_DN(N,1:N_GEOMETRIES,:) = VLIDORT_FixIn%Optical%TS_FMATRIX_DN(N,1:N_GEOMETRIES,:)
         ENDDO
      ENDIF
      IF ( .not. VLIDORT_Bookkeep%DO_FOCORR_ALONE ) THEN
        DO N = 1, NLAYERS
          GREEKMAT_TOTAL_INPUT(0:NSTREAMS_2,N,:) = VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:NSTREAMS_2,N,:)
        ENDDO
      ENDIF

!  Fixed Write inputs
!  ------------------

      DO_DEBUG_WRITE          = VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE

      DO_WRITE_INPUT          = VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT
      INPUT_WRITE_FILENAME    = VLIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME

      DO_WRITE_SCENARIO       = VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO
      SCENARIO_WRITE_FILENAME = VLIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME

      DO_WRITE_FOURIER        = VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER
      FOURIER_WRITE_FILENAME  = VLIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME

      DO_WRITE_RESULTS        = VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS
      RESULTS_WRITE_FILENAME  = VLIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME

!  Modified Boolean inputs
!  -----------------------

!  FOCORR and SSCORR Booleans. Completely reorganized for Version 2.8
!mick mod 9/19/2017 - DO_FOCORR_ALONE now defined internally

      DO_FOCORR               = VLIDORT_ModIn%MBool%TS_DO_FOCORR            !New 02 Jul 2013
      DO_FOCORR_EXTERNAL      = VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL
      !DO_FOCORR_ALONE         = VLIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE
      DO_FOCORR_NADIR         = VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR
      DO_FOCORR_OUTGOING      = VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING

!  4/15/20. Version 2.8.2. Flags removed
!      DO_SSCORR_TRUNCATION    = VLIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION
!      DO_SSCORR_USEFMAT       = VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT

!  Other Booleans
!mick fix 3/2/2020 - added defining of DO_NO_AZIMUTH since it is now done outside of VLIDORT

      DO_SOLAR_SOURCES        = VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES
      DO_REFRACTIVE_GEOMETRY  = VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY
      DO_CHAPMAN_FUNCTION     = VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION

      DO_RAYLEIGH_ONLY        = VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY
      DO_DELTAM_SCALING       = VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING
      DO_DOUBLE_CONVTEST      = VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST
      DO_SOLUTION_SAVING      = VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING
      DO_BVP_TELESCOPING      = VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING

      DO_NO_AZIMUTH           = VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH

      DO_USER_VZANGLES        = VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES
      DO_ADDITIONAL_MVOUT     = VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT
      DO_MVOUT_ONLY           = VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

      DO_THERMAL_TRANSONLY    = VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY
      DO_OBSERVATION_GEOMETRY = VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY

!  4/15/20. Version 2.8.2. Add Doublet geometry flag

      DO_DOUBLET_GEOMETRY     = VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY

!  Additional Control for Externalized Water-leaving input. Introduced 3/18/19 for Version 2.8.1

      DO_EXTERNAL_WLEAVE     = VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE

!  Modified Beam inputs
!  --------------------

      N_SZANGLES             = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
      SZANGLES(1:N_SZANGLES) = VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:N_SZANGLES)

!  Modified User Value inputs
!  --------------------------

      N_USER_RELAZMS                 = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      USER_RELAZMS(1:N_USER_RELAZMS) = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS)

      N_USER_VZANGLES                  = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
      USER_VZANGLES(1:N_USER_VZANGLES) = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:N_USER_VZANGLES)

      USER_LEVELS(1:N_USER_LEVELS)     = VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS)

!  4/15/20.  Version 2.8.2. Type structure change for GEOMETRY_SPECHEIGHT

      GEOMETRY_SPECHEIGHT              = VLIDORT_ModIn%MCont%TS_GEOMETRY_SPECHEIGHT

      N_USER_OBSGEOMS                      = VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
      USER_OBSGEOMS(1:N_USER_OBSGEOMS,1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,1:3)

!  4/15/20. Version 2.8.2. Modified Control/Chapman/Optical inputs no longer required here
!      NGREEK_MOMENTS_INPUT = VLIDORT_ModIn%MCont%TS_NMOMENTS_INPUT                              ! Not needed now
!      CHAPMAN_FACTORS      = VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS                          ! Output from geometry
!      EARTH_RADIUS         = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS                             ! Input  to   geometry
!      OMEGA_TOTAL_INPUT(1:NLAYERS)  = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS)    ! Moved to FixIn%Optical

!  BRDF inputs
!  -----------

!mick mod 9/19/2017 - added IF conditions
!                   - note: emissivities left out of BRDF SURFACE block due to possibility
!                           of thermal being used in the Lambertian case

!  4/15/20. Version 2.8.2. BRDF copying only direct BRDF, Fourier BRDF copying moved into Fourier loop.

      IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
        IF ( DO_USER_VZANGLES ) THEN
          EXACTDB_BRDFUNC (:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
            VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC(:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
        ENDIF
      ENDIF

!  @@@ Rob fix 1/31/11, Emissivities from earlier structure
!        Lambertian case was OK, as used internal definitions

      IF ( DO_SURFACE_EMISSION ) THEN
        EMISSIVITY(1:NSTOKES,1:NSTREAMS) = VLIDORT_Sup%BRDF%TS_EMISSIVITY(1:NSTOKES,1:NSTREAMS)
        IF ( DO_USER_VZANGLES ) THEN
          USER_EMISSIVITY(1:NSTOKES,1:N_USER_VZANGLES) = VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY(1:NSTOKES,1:N_USER_VZANGLES)
        ENDIF
      ENDIF

!  SLEAVE inputs  (This code introduced 17 May 2012)
!  -------------

!mick mod 9/19/2017 - added IF conditions

!  4/15/20. Version 2.8.2. SLEAVE copying only for Isotropic and direct TERM.
!                                Fourier SLEAVE copying now moved into Fourier loop.

      IF ( DO_SURFACE_LEAVING ) THEN
        SLTERM_ISOTROPIC(1:NSTOKES,1:N_SZANGLES) = &
          VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(1:NSTOKES,1:N_SZANGLES)
        IF ( DO_USER_VZANGLES ) THEN
          SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
            VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
        ENDIF
      ENDIF

!  SS inputs
!  ---------

!  New 12 March 2012 --> IF SS results already available copy them. Modified flagging, Version 2.8
!  4/15/20. Version 2.8.2. No Longer requires this code, now using type structure inputs directly.
!      IF ( DO_FOCORR_EXTERNAL ) THEN
!        STOKES_SS(1:N_USER_LEVELS,:,1:NSTOKES,:) = &
!          VLIDORT_Sup%SS%TS_STOKES_SS(1:N_USER_LEVELS,:,1:NSTOKES,:)
!        STOKES_DB(1:N_USER_LEVELS,:,1:NSTOKES) = &
!          VLIDORT_Sup%SS%TS_STOKES_DB(1:N_USER_LEVELS,:,1:NSTOKES)
!      ENDIF
     
! #################### INNOVATIONS 5/5/20 ##############
!  4/29/20. Set Flag: Special Fourier-component output for Rayleigh + Planetary-problem TOA Upwelling situations
!      SPECIAL_RAYF_OUTPUT = &
!       ( DO_PLANETARY_PROBLEM .and.DO_UPWELLING .and..not.DO_FOCORR .and.DO_RAYLEIGH_ONLY .and. (USER_LEVELS(1).eq.zero) )
! #################### INNOVATIONS 5/5/20 ##############

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  ====================================================
!  BEGIN COPY BOOKKEEPING VARIABLES TO LOCAL Quantities
!  ====================================================

!  MSMode of operation, Add MSMODE_THERMAL flag (1/13/12. Version 3.6)
!mick mod 3/22/2017 - DO_FOCORR_ALONE now defined internally, as bookkeeping variable

      DO_MSMODE_VLIDORT = VLIDORT_Bookkeep%DO_MSMODE_VLIDORT
      DO_MSMODE_THERMAL = VLIDORT_Bookkeep%DO_MSMODE_THERMAL
      DO_FOCORR_ALONE   = VLIDORT_Bookkeep%DO_FOCORR_ALONE

!  NSTOKES etc...

      NSTOKES_SQ     = NSTOKES * NSTOKES
      NSTKS_NSTRMS   = VLIDORT_Bookkeep%NSTKS_NSTRMS
      NSTKS_NSTRMS_2 = VLIDORT_Bookkeep%NSTKS_NSTRMS_2

!  Actual number of moments =  2 x NSTREAMS - 1 ; NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      NMOMENTS   = VLIDORT_Bookkeep%NMOMENTS
      NTOTAL     = VLIDORT_Bookkeep%NTOTAL
      N_SUBDIAG  = VLIDORT_Bookkeep%N_SUBDIAG
      N_SUPDIAG  = VLIDORT_Bookkeep%N_SUPDIAG

!  Output optical depth masks and indices

      UTAU_LEVEL_MASK_UP(1:N_USER_LEVELS) = VLIDORT_Bookkeep%UTAU_LEVEL_MASK_UP(1:N_USER_LEVELS)
      UTAU_LEVEL_MASK_DN(1:N_USER_LEVELS) = VLIDORT_Bookkeep%UTAU_LEVEL_MASK_DN(1:N_USER_LEVELS)

!  Off-grid optical depths (values, masks, indices)
!   --- Need to zero-out the non-active partial layer values. 12/15/19.

      DO_PARTLAYERS =  VLIDORT_Bookkeep%DO_PARTLAYERS
      IF ( DO_PARTLAYERS ) THEN
         N_PARTLAYERS = VLIDORT_Bookkeep%N_PARTLAYERS
         PARTLAYERS_LAYERIDX (1:N_PARTLAYERS)  = VLIDORT_Bookkeep%PARTLAYERS_LAYERIDX (1:N_PARTLAYERS)
         PARTLAYERS_VALUES   (1:N_PARTLAYERS)  = VLIDORT_Bookkeep%PARTLAYERS_VALUES   (1:N_PARTLAYERS)
         PARTLAYERS_OUTFLAG  (1:N_USER_LEVELS) = VLIDORT_Bookkeep%PARTLAYERS_OUTFLAG  (1:N_USER_LEVELS) 
         PARTLAYERS_OUTINDEX (1:N_USER_LEVELS) = VLIDORT_Bookkeep%PARTLAYERS_OUTINDEX (1:N_USER_LEVELS)
      else
         N_PARTLAYERS = 0
         PARTLAYERS_LAYERIDX = 0 
         PARTLAYERS_VALUES   = ZERO
         PARTLAYERS_OUTFLAG  = .FALSE.
         PARTLAYERS_OUTINDEX = 0
      endif

!  Layer masks and limits for doing integrated source terms
!mick fix 3/2/2020 - added N_ALLLAYERS_UP & N_ALLLAYERS_DN

      STERM_LAYERMASK_UP(1:NLAYERS) = VLIDORT_Bookkeep%STERM_LAYERMASK_UP(1:NLAYERS)
      STERM_LAYERMASK_DN(1:NLAYERS) = VLIDORT_Bookkeep%STERM_LAYERMASK_UP(1:NLAYERS)

      N_ALLLAYERS_UP = VLIDORT_Bookkeep%N_ALLLAYERS_UP
      N_ALLLAYERS_DN = VLIDORT_Bookkeep%N_ALLLAYERS_DN

!  Post-processing masks. Introduced 4/9/19. Offsets for geometry indexing
!mick fix 3/2/2020 - trimmed dimensions on bookkeeping arrays in the next three subsections here
!                  - added IF condition for VZA_OFFSETS
!    4/15/20. Version 2.8.2. Add Logic for Doublet geometry offset

      N_PPSTREAMS = VLIDORT_Bookkeep%N_PPSTREAMS
      PPSTREAM_MASK(1:N_USER_VZANGLES,1:N_SZANGLES) = VLIDORT_Bookkeep%PPSTREAM_MASK(1:N_USER_VZANGLES,1:N_SZANGLES)

      SZD_OFFSETS(1:N_SZANGLES)   = 0
      VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES) = 0

      IF ( DO_DOUBLET_GEOMETRY ) THEN
        SZD_OFFSETS(1:N_SZANGLES) = VLIDORT_Bookkeep%SZD_OFFSETS(1:N_SZANGLES)
      ELSE IF ( .not.DO_OBSERVATION_GEOMETRY ) THEN
        VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES) = VLIDORT_Bookkeep%VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES)
      ENDIF

!  Number of directions (1 or 2) and directional array
!  Number of convergence tests, Azimuths, Geometries
!mick fix 3/2/2020 - defined N_OUT_STREAMS

      N_DIRECTIONS     = VLIDORT_Bookkeep%N_DIRECTIONS
      WHICH_DIRECTIONS(1:N_DIRECTIONS) = VLIDORT_Bookkeep%WHICH_DIRECTIONS(1:N_DIRECTIONS)

      N_OUT_STREAMS    = VLIDORT_Bookkeep%N_OUT_STREAMS
      N_CONVTESTS      = VLIDORT_Bookkeep%N_CONVTESTS
      LOCAL_N_USERAZM  = VLIDORT_Bookkeep%LOCAL_N_USERAZM

!  Polarization quantities
!mick note 3/2/2020 - particular solution variables not enabled --> turned off

      MUELLER_INDEX(1:NSTOKES,1:NSTOKES) = VLIDORT_Bookkeep%MUELLER_INDEX(1:NSTOKES,1:NSTOKES)
      GREEKMAT_INDEX                     = VLIDORT_Bookkeep%GREEKMAT_INDEX
      DMAT(1:NSTOKES,1:NSTOKES)          = VLIDORT_Bookkeep%DMAT(1:NSTOKES,1:NSTOKES)
      DFLUX(1:NSTOKES)                   = VLIDORT_Bookkeep%DFLUX(1:NSTOKES)
      FLUXVEC(1:NSTOKES)                 = VLIDORT_Bookkeep%FLUXVEC(1:NSTOKES)
      !NPARTICSOLS                        = VLIDORT_Bookkeep%NPARTICSOLS
      !DMAT_PSOLS(1:NSTOKES,1:NSTOKES,2)  = VLIDORT_Bookkeep%DMAT_PSOLS(1:NSTOKES,1:NSTOKES,2)
      !DMAT_PSOLS_FLUXVEC(1:NSTOKES,2)    = VLIDORT_Bookkeep%DMAT_PSOLS_FLUXVEC(1:NSTOKES,2)

!  Misc control flags
!mick fix 3/2/2020 - defined DO_ALL_FOURIER & DO_DBCORRECTION

      DO_ALL_FOURIER  = VLIDORT_Bookkeep%DO_ALL_FOURIER
      DO_DBCORRECTION = VLIDORT_Bookkeep%DO_DBCORRECTION

!  ====================================================
!  END COPY BOOKKEEPING VARIABLES TO LOCAL Quantities
!  ====================================================

!  ====================================================
!  BEGIN COPY MS-GEOMETRY VARIABLES TO LOCAL Quantities
!  ====================================================

!  Local input solar zenith angles Cosines. New usage for direct-beam radiance, 4/9/19
!    (Also required for refractive geometry attenuation of the solar beam)
!mick fix 3/2/2020 - added DO_SOLAR_SOURCES & DO_CHAPMAN_FUNCTION IF conditions

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_CHAPMAN_FUNCTION ) THEN
          SUNLAYER_COSINES(1:NLAYERS,1:N_SZANGLES)= VLIDORT_MSGeom%SUNLAYER_COSINES(1:NLAYERS,1:N_SZANGLES)
        ENDIF
      ENDIF

!  Local solar zenith angles Cosines (regular case)
!mick fix 3/2/2020 - added DO_SOLAR_SOURCES IF condition

      !IF ( DO_SOLAR_SOURCES ) THEN
        COS_SZANGLES(1:N_SZANGLES) = VLIDORT_MSGeom%COS_SZANGLES(1:N_SZANGLES)
        SIN_SZANGLES(1:N_SZANGLES) = VLIDORT_MSGeom%SIN_SZANGLES(1:N_SZANGLES)
      !ENDIF

!  Quadrature weights and abscissae, and product

      QUAD_STREAMS(1:NSTREAMS) = VLIDORT_MSGeom%QUAD_STREAMS(1:NSTREAMS)
      QUAD_WEIGHTS(1:NSTREAMS) = VLIDORT_MSGeom%QUAD_WEIGHTS(1:NSTREAMS)
      QUAD_STRMWTS(1:NSTREAMS) = VLIDORT_MSGeom%QUAD_STRMWTS(1:NSTREAMS)
      QUAD_HALFWTS(1:NSTREAMS) = VLIDORT_MSGeom%QUAD_HALFWTS(1:NSTREAMS)
      QUAD_SINES  (1:NSTREAMS) = VLIDORT_MSGeom%QUAD_SINES  (1:NSTREAMS)
      QUAD_ANGLES (1:NSTREAMS) = VLIDORT_MSGeom%QUAD_ANGLES (1:NSTREAMS)

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles
!mick fix 3/2/2020 - add IF condition

      IF ( DO_USER_VZANGLES ) THEN
        USER_ANGLES  (1:N_USER_VZANGLES) = VLIDORT_MSGeom%USER_ANGLES  (1:N_USER_VZANGLES)
        USER_STREAMS (1:N_USER_VZANGLES) = VLIDORT_MSGeom%USER_STREAMS (1:N_USER_VZANGLES)
        USER_SINES   (1:N_USER_VZANGLES) = VLIDORT_MSGeom%USER_SINES   (1:N_USER_VZANGLES)
        USER_SECANTS (1:N_USER_VZANGLES) = VLIDORT_MSGeom%USER_SECANTS (1:N_USER_VZANGLES)
      ENDIF

!  Chapman factors
!mick fix 3/2/2020 - trim dimensions on CHAPMAN_FACTORS
!                  - added DO_SOLAR_SOURCES IF condition

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_CHAPMAN_FUNCTION ) THEN
          DO IB = 1, N_SZANGLES
            !CHAPMAN_FACTORS(1:NLAYERS,1:NLAYERS,IB) = VLIDORT_MSGeom%CHAPMAN_FACTORS(1:NLAYERS,1:NLAYERS,IB)
            DO N = 1, NLAYERS
              CHAPMAN_FACTORS(N,1:N,IB) = VLIDORT_MSGeom%CHAPMAN_FACTORS(N,1:N,IB)
            ENDDO
            IF ( DO_PARTLAYERS ) THEN
              !PARTIAL_CHAPFACS(1:N_PARTLAYERS,1:NLAYERS,IB) = VLIDORT_MSGeom%PARTIAL_CHAPFACS(1:N_PARTLAYERS,1:NLAYERS,IB)
              DO UT = 1, N_PARTLAYERS
                N = PARTLAYERS_LAYERIDX(UT)
                DO K = 1, N
                  PARTIAL_CHAPFACS(UT,K,IB) = VLIDORT_MSGeom%PARTIAL_CHAPFACS(UT,K,IB)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDIF

!  ====================================================
!  END COPY MS-GEOMETRY VARIABLES TO LOCAL Quantities
!  ====================================================

!  VLIDORT input debug
!  ===================

!  4/15/20. Version 2.8.2. Need to set some local variables, when debugging
!     Some geometry, BRDF, and SLEAVE inputs and SS/DB results
!mick fix 3/2/2020 - added BRDF and SLEAVE blocks here to still enable a spotcheck
!                    on the Fourier=0 moments of these arrays; modified above comments

      IF (DO_DEBUG_INPUT) THEN

!  Set a few inputs from the geometry setup routines

         PRESSURE_GRID   (0:NLAYERS) = VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID(0:NLAYERS)
         TEMPERATURE_GRID(0:NLAYERS) = VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(0:NLAYERS)
         FINEGRID(1:NLAYERS)         = VLIDORT_FixIn%Chapman%TS_FINEGRID(1:NLAYERS)
         RFINDEX_PARAMETER           = VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER
         EARTH_RADIUS                = VLIDORT_ModIn%MCont%TS_EARTH_RADIUS

!  Set local BRDF Fourier arrays for Fourier=0

         IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
           FOURIER = 0
           BRDF_F_0(1:MAXSTOKES_SQ,1:NSTREAMS,1:N_SZANGLES) = &
             VLIDORT_Sup%BRDF%TS_BRDF_F_0(FOURIER,1:MAXSTOKES_SQ,1:NSTREAMS,1:N_SZANGLES)
           BRDF_F  (1:MAXSTOKES_SQ,1:NSTREAMS,1:NSTREAMS)   = &
             VLIDORT_Sup%BRDF%TS_BRDF_F  (FOURIER,1:MAXSTOKES_SQ,1:NSTREAMS,1:NSTREAMS)
           IF ( DO_USER_VZANGLES ) THEN
             USER_BRDF_F_0 (1:MAXSTOKES_SQ,1:N_USER_VZANGLES,1:N_SZANGLES) = &
               VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0(FOURIER,1:MAXSTOKES_SQ,1:N_USER_VZANGLES,1:N_SZANGLES)
             USER_BRDF_F   (1:MAXSTOKES_SQ,1:N_USER_VZANGLES,1:NSTREAMS)   = &
               VLIDORT_Sup%BRDF%TS_USER_BRDF_F  (FOURIER,1:MAXSTOKES_SQ,1:N_USER_VZANGLES,1:NSTREAMS)
           ENDIF
         ENDIF

!  New 12 March 2012 --> IF SS results already available, copy them otherwise initialize

         IF ( DO_FOCORR_EXTERNAL ) THEN
            STOKES_SS = VLIDORT_Sup%SS%TS_STOKES_SS ; STOKES_DB = VLIDORT_Sup%SS%TS_STOKES_DB
         ELSE
            STOKES_SS = ZERO ; STOKES_DB = ZERO
         ENDIF

!  Set local SLEAVE Fourier arrays for Fourier=0

         IF ( DO_SURFACE_LEAVING ) THEN
           FOURIER = 0
           SLTERM_F_0(1:MAXSTOKES,1:NSTREAMS,1:N_SZANGLES) = &
             VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0(FOURIER,1:MAXSTOKES,1:NSTREAMS,1:N_SZANGLES)
           IF ( DO_USER_VZANGLES ) THEN
             USER_SLTERM_F_0(1:MAXSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES) = &
               VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(FOURIER,1:MAXSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES)
           ENDIF
         ENDIF

!  Debug Call

         CALL VLIDORT_DEBUG_INPUT_MASTER()

      END IF

!  Pre-Zeroing (initialize outputs)
!  ================================

!  Main outputs (Radiances and fluxes)

!  4/15/20. Version 2.8.2. No Longer required, now using type structure inputs directly.
!      STOKES(1:N_USER_LEVELS,:,1:NSTOKES,:) = ZERO

      MEANST_DIFFUSE (1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = ZERO
      FLUX_DIFFUSE   (1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = ZERO
      DNMEANST_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)   = ZERO
      DNFLUX_DIRECT  (1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)   = ZERO

!  Special Media-property output. -- Introduced 4/26/19 R. Spurr. Pre-initialized here.
!     ** Output for User-angle streams, also fluxes. TRANSBEAM for the planetary problem.

      ALBMED_USER = zero ; ALBMED_FLUXES = zero
      TRNMED_USER = zero ; TRNMED_FLUXES = zero
      TRANSBEAM   = zero

!  4/28/19. Initialize the planetary problem outputs
!  4/15/20. Version 2.8.2. No Longer required, Should be pre-zeroed.
      !VLIDORT_Out%Main%TS_PLANETARY_SBTERM    = ZERO
      !VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM = ZERO

!  Checking of Inputs
!  ------------------

!  4/15/20. Version 2.8.2. Only the wavelength dependent optical checks survive.
!                                All other input checking is now done in SETUP_MASTER

!  Check input optical values (IOPs, albedo)

      LOCAL_DO_MEDIA_PROBLEM = DO_PLANETARY_PROBLEM .or. DO_ALBTRN_MEDIA(1) .or. DO_ALBTRN_MEDIA(2) 
      CALL VLIDORT_CHECK_INPUT_OPTICAL &
         ( NLAYERS, DO_THERMAL_TRANSONLY, LOCAL_DO_MEDIA_PROBLEM, ALBEDO,  & ! Input
           DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,     & ! Input
           STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS)               ! Input/Output

!  Exception handling

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VLIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        VLIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        VLIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  4/15/20 Version 2.8.2. Define TAUGRID_INPUT. Added here, formerly in DERIVE_INPUT

      TAUGRID_INPUT(0) = ZERO
      DO N = 1, NLAYERS
        TAUGRID_INPUT(N) = TAUGRID_INPUT(N-1) + DELTAU_VERT_INPUT(N)
      END DO

!  Preparation for Fourier call
!  ----------------------------

!  if there's no azimuth dependence, just do one value in azimuth loop
!  4/15/20. Version 2.8.2. Included in new DERIVE_INPUT subroutine 

!  Local post-processing control. Added, 4/9/19
!  4/15/20. Version 2.8.2. Included in new DERIVE_INPUT subroutine 

!  4/15/20. Version 2.8.2. Included in new DERIVE_INPUT subroutine 
!  Number of sources --> Save some offsets for indexing geometries (lattice only)
!   This section revised for the Observational Geometry option

!  4/15/20. Version 2.8.2. Moved in new DERIVE_INPUT subroutine 
!  Get derived inputs ==> Miscellaneous and layer input.
!   ( 6/21/10)  Important point: Must use ADJUSTED VZangles as input here.
!   ( 11/19/14) DISAGREE WITH THIS POINT------------------!!!!!!!!!!!!!!!!!
!    Geometry adjustment, Coding removed for Version 2.8
!   revised list, Version 2.8, 3/3/17 !mick fix 3/22/2017 - added DO_FOCORR_NADIR & DO_FOCORR_OUTGOING to input

!  4/15/20. Version 2.8.2. Chapman function calculation now done in SETUP_MASTER
!mick fix 1/19/2018 - moved the call to LIDORT_CHAPMAN from before LIDORT_DERIVE_INPUT
!                     to here to resolve an I/O contradiction

!  #################
!  Set up operations
!  #################

!  Setups for Fourier = 0
!  ======================

!  Each call is followed by a packing routine

!  1.  VLIDORT_MISCSETUPS  : Performance, Delta-M, average-secant formulation, transmittances
!  2.  THERMAL_SETUP       : Coefficients, direct multipliers
!  3a. EMULT_MASTER /      : Beam source function multipliers.  Not required for the
!  3b. EMULT_MASTER_OBSGEO     Full SS calculation in outgoing mode

!  1. VLIDORT MISCSETUPS
!  ---------------------

!  1/30/08. Telescoping setup is done in DERIVE_INPUTS (unlike LIDORT scalar code)
!  7/8/16 . Revision of I/O lists for Version 2.8.
!mick fix 9/19/2017 - the following arguments added to facilitate correction of direct flux:
!                     to input  - DO_SOLAR_SOURCES, PARTIAL_CHAPFACS
!                     to output - DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED, 
!                                 LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS

!  4/15/20. Version 2.8.2. Now added additional performance stuff.

      CALL VLIDORT_MISCSETUPS ( &
        DO_SOLAR_SOURCES, DO_SOLUTION_SAVING,   DO_BVP_TELESCOPING, DO_DELTAM_SCALING,      & ! Input flags
        DO_RAYLEIGH_ONLY, DO_THERMAL_TRANSONLY, DO_PLANE_PARALLEL,  DO_REFRACTIVE_GEOMETRY, & ! Input flags
        DO_PARTLAYERS, DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY,        & ! Input flags
        DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3, DO_TOA_CONTRIBS, & ! Input flags
        NSTOKES, NLAYERS, NSTREAMS, N_USER_VZANGLES,                     & ! Input Numbers
        N_SZANGLES, NMOMENTS, NLAYERS_CUTOFF, NLAYERS_NOMS,              & ! Input Numbers
        MUELLER_INDEX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,         & ! Input level/Stokes indices
        N_PARTLAYERS,  PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,         & ! Input Partlayers
        QUAD_STREAMS, USER_SECANTS, COS_SZANGLES, SUN_SZA_COSINES,     & ! Input streams, SZA cosines
        OMEGA_TOTAL_INPUT, DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT,    & ! Input Optical
        TAUGRID_INPUT, CHAPMAN_FACTORS, PARTIAL_CHAPFACS,              & ! Input Chapman
        DO_LAYER_SCATTERING, BVP_REGULAR_FLAG, DO_REAL_EIGENSOLVER,    & ! Output from PERFORMANCE
        DELTAU_VERT, PARTAU_VERT, OMEGA_TOTAL, GREEKMAT_TOTAL,         & ! Output from DELTAMSCALE
        TAUGRID, DELTAU_SLANT, DELTAU_SLANT_UNSCALED,                  & ! Output from DELTAMSCALE
        PARTAU_SLANT_UNSCALED, LEVELS_SOLARTRANS,                      & ! Output from DELTAMSCALE
        PARTIALS_SOLARTRANS, SOLARBEAM_BOATRANS, TRUNC_FACTOR, FAC1,   & ! Output from DELTAMSCALE
        OMEGA_GREEK,                                                   & ! Output from SSALBINIT
        DO_REFLECTED_DIRECTBEAM, BEAM_CUTOFF, TRANS_SOLAR_BEAM,        & ! Output QSPREP
        INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA,                     & ! output QSPREP
        T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                & ! Output PREPTRANS (Discrete Ords.)
        T_DELT_MUBAR, T_UTUP_MUBAR, T_UTDN_MUBAR,                      & ! Output PREPTRANS (Solar beams)
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                      & ! Output PREPTRANS (User-vza)
        CUMTRANS, ITRANS_USERM,                                        & ! Output PREPTRANS (auxiliary)       
        STATUS_SUB, MESSAGE, TRACE_1 )                                   ! Output PERFORMANCE
      
!  4/15/20. Version 2.8.2. Now added exception handling (for performance output).

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS) THEN
         TRACE_2 = 'PERFORMANCE SETUP failed, called by VLIDORT_MISCSETUPS'
         TRACE_3 = 'VLIDORT_MISCSETUPS failed, called by VLIDORT_RTCALC_MASTER'
         STATUS_CALCULATION = VLIDORT_SERIOUS
         VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
         VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
         VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
         VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
         VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
         RETURN
       ENDIF

!  Pack selected MISCSETUPS results

      CALL VLIDORT_PACK_MISC ( &
        NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS, N_USER_VZANGLES, N_SZANGLES, NMOMENTS, & ! Input Numbers
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT,                                          & ! Input Optical
        LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, OMEGA_GREEK,                             & ! Input Optical
        BEAM_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,                          & ! Input Solar
        T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                                  & ! Input Trans. D.O.
        T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,            & ! Input Trans solar/User
        CUMTRANS, INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA,                             & ! Input Av. secant.
        Misc )                                                                             ! Output structure to be packed

!  2. THERMAL SETUPS
!  -----------------

!  7/8/16 . Revision of I/O lists for Version 2.8.

      IF ( DO_THERMAL_EMISSION ) THEN

        CALL THERMAL_SETUP ( &
          DO_USER_VZANGLES, DO_UPWELLING, DO_DNWELLING,                & ! Flags
          DO_PARTLAYERS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,      & ! Flags
          NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES,    & ! Numbers basic
          THERMAL_BB_INPUT, USER_STREAMS,                              & ! thermal input, streams
          PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Level control
          OMEGA_TOTAL, DELTAU_VERT, PARTAU_VERT,                       & ! Input optical
          T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                    & ! Input transmittances
          THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,                & ! output thermal setups
          T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN )     ! output thermal direct solutions

        CALL VLIDORT_PACK_THERM ( &
          NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, & ! Input Numbers
          THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Input Thermal setup
          T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN, & ! Input Thermal solutions
          Therm )                                                     ! Output

      END IF

!  3a/b. EMULT_MASTERS.
!  --------------------

!  7/8/16 . Revision of I/O lists for Version 2.8.

!mick fix 3/30/2015 - modified if condition
      !IF ( DO_SOLAR_SOURCES ) THEN
!Rob fix 3/1/2017 - modified if condition
      !IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN

      IF ( DO_SOLAR_SOURCES .AND. DO_USER_VZANGLES) THEN
        IF (.NOT.DO_FOCORR_ALONE) THEN

          IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
            CALL EMULT_MASTER ( &
                DO_UPWELLING, DO_DNWELLING,                                       & ! Input flags
                NLAYERS, N_SZANGLES, N_USER_VZANGLES, TAYLOR_ORDER, N_PARTLAYERS, & ! Input numbers
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input Level control
                USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                           & ! Input optical, streams   
                BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,          & ! Input solar beam
                T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, ITRANS_USERM,           & ! Input User-stream Trans.
                EMULT_HOPRULE, SIGMA_M, SIGMA_P,                                  & ! Output Multipliers
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )                      ! Output Multipliers
          ELSE
            CALL EMULT_MASTER_OBSGEO ( &
                DO_UPWELLING, DO_DNWELLING,                                  & ! Input flags
                NLAYERS, N_SZANGLES, TAYLOR_ORDER, N_PARTLAYERS,             & ! Input numbers
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Input Level control
                USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                      & ! Input optical, streams   
                BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,     & ! Input solar beam
                T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, ITRANS_USERM,      & ! Input User-stream Trans.
                EMULT_HOPRULE, SIGMA_M, SIGMA_P,                             & ! Output Multipliers
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )                 ! Output Multipliers
          ENDIF

          CALL VLIDORT_PACK_MULT ( &
            NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES, & ! Input
            EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,       & ! Input
            Mult )                                                ! Output

        ENDIF
      ENDIF

!  WATER-LEAVING TRANSMITTANCE CALCULATION
!  =======================================
      
!  - 04/09/19. Radical new departure. Now moved to the BVPROBLEM as an adjusted Backsub calculation
!              Scaling of SLTERMS by TRANS_ATMOS_FINAL is now done either inside the Fourier routine (MS)
!              or below for the direct term contribution (After Fourier call for Fourier 0)
!              The  LIDORT_TRANSFLUX_MASTER routine has been removed.


!  New section for Version 2.7a and 2.8
!  - 09/25/15. First programmed by R. Spurr for Version 2.7a, RT Solutions Inc.
!  - 12/24/15. Drop SL_Isotropic Constraint (if non-Isotropy, additional SL terms need transmittance scaling)
!  - 02/03/16. Mark 1 and Mark2 codes (Dark Surface Result). COMMENTED OUT IN FAVOR of Mark 3
!  - 07/08/16. Mark 3 code given iteration control (3 new inputs).

!  HERE IS THE OLD CODE.....
!      IF ( DO_SURFACE_LEAVING .and. DO_WATER_LEAVING .and. DO_SL_ISOTROPIC ) THEN
!      IF ( DO_SURFACE_LEAVING .and. DO_WATER_LEAVING ) THEN
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         Mark 1, Mark2 code, Dark Surface Result. COMMENTED OUT 2/3/16
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         CALL VLIDORT_TRANSFLUX & 
!           ( NSTOKES, NSTREAMS, NLAYERS, N_SZANGLES, NMOMENTS, NSTREAMS_2,           & ! Input
!             NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,         & ! Input
!             FLUX_FACTOR, FLUXVEC, COS_SZANGLES, MUELLER_INDEX, DMAT, DFLUX,     & ! Input
!             QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS, LOCAL_CSZA, & ! Input
!             DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, DELTAU_VERT, OMEGA_GREEK, & ! Input
!             BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT,      & ! Input
!             MEANST_DIFFUSE_TF, FLUX_DIFFUSE_TF, DNMEANST_DIRECT_TF, DNFLUX_DIRECT_TF,     & ! InOut
!             STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                               ! Output
!       write(*,*)'MK2 Result (Dark)   = ',0,FLUX_DIFFUSE_TF(1,1),FLUX_DIFFUSE_TF(1,1) / flux_factor / local_csza(nlayers,1)
!  debug output. 17 december 2015.
!     It is the last entry that is correct.....
!         write(779,'(A,A)')'   SZA    COSSZA   F0      Fdirect      Tdirect      Fdiffuse   ',&
!                                                    '  TDiffuse     Fcomplete    Tcomplete  [T* = F*/F0/cossza]'
!         do ibeam = 1, N_SZANGLES
!!            write(*,'(a,i2,2(1p3e15.6,2x))')'Transflux',ibeam,MEANST_DIFFUSE_TF(ibeam,1:3),FLUX_DIFFUSE_TF(ibeam,1:3)
!             write(779,'(F8.3,F9.5,F6.2,2x,1p6e13.5)')szangles(ibeam),local_csza(nlayers,ibeam),flux_factor,&
!              DNFLUX_DIRECT_TF(ibeam,1),DNFLUX_DIRECT_TF(ibeam,1)/flux_factor/local_csza(nlayers,ibeam),&
!              FLUX_DIFFUSE_TF(ibeam,1)-DNFLUX_DIRECT_TF(ibeam,1),&
!              (FLUX_DIFFUSE_TF(ibeam,1)-DNFLUX_DIRECT_TF(ibeam,1))/flux_factor/local_csza(nlayers,ibeam),&
!              FLUX_DIFFUSE_TF(ibeam,1),FLUX_DIFFUSE_TF(ibeam,1)/flux_factor/local_csza(nlayers,ibeam)
!         enddo
!  error handling
!         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
!            TRACE_3 = ' Called by VLIDORT_MASTER ' ; STATUS_CALCULATION = VLIDORT_SERIOUS
!            VLIDORT_Out%Status%TS_MESSAGE = MESSAGE ; VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!            VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2 ; VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!            VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION ; RETURN
!         ENDIF
!  Scale = Flux/F0/mu0 = Atmospheric Transmittance (Dark Result)
!         TRANS_ATMOS(1:N_SZANGLES) = FLUX_DIFFUSE_TF(1:N_SZANGLES,1) / flux_factor / local_csza(nlayers,1:N_SZANGLES)
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         Mark 3 code, Self-Consistent Result, by iteration.   2/3/16
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Very important to set the surface flag
!         DO_INCLUDE_SURFACE = .true.
!   Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination
!  --- Also include value of this illumination
!  --- Also need the discrete ordinate transmittances  
!         CALL VLIDORT_TRANSFLUX_MASTER &
!        ( DO_TF_ITERATION, TF_MAXITER, TF_CRITERION, DO_TOAFLUX, TOAFLUX,     & ! Input TF control, TOAISO (3/23/19 new)
!          NSTOKES, NSTREAMS, NLAYERS, N_SZANGLES, NMOMENTS, NSTREAMS_2,                 & ! Input Numbers
!          NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,               & ! Input Numbers
!          FLUX_FACTOR, FLUXVEC, COS_SZANGLES, MUELLER_INDEX, DMAT, DFLUX,           & ! Input Bookkeeping
!          DO_INCLUDE_SURFACE, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,            & ! Input surface 1/8/16 
!          DO_LAMBERTIAN_SURFACE, ALBEDO, BRDF_F, BRDF_F_0,               & ! Input surface 1/8/16
!          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_F_0,        & ! Input surface 1/8/16
!          QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS, LOCAL_CSZA,       & ! Input Quadrature
!          DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, DELTAU_VERT, OMEGA_GREEK,       & ! Input Optical
!          T_DELT_DISORDS, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT, & ! Input solar beam, 3/23/19 new
!          TRANS_ATMOS_FINAL, FLUX_DIFFUSE_FINAL,                                    & ! Output
!          STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                                     ! Output
!  Exception handling
!            IF ( STATUS_SUB == VLIDORT_SERIOUS ) THEN
!               TRACE_3 = 'VLIDORT_TRANSFLUX_MASTER failed, VLIDORT_MASTER'
!               STATUS_CALCULATION = VLIDORT_SERIOUS
!               VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
!               VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
!               VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!               VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
!               VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!               RETURN
!            ENDIF
!  Scale the Isotropic term
!         do ibeam = 1, N_SZANGLES
!            SLTERM_ISOTROPIC(1,IBEAM) = SLTERM_ISOTROPIC(1,IBEAM) * TRANS_ATMOS_FINAL(IBEAM)
!         enddo
!  Scale the Non-isotropic terms, 24 December 2015, if flagged
!         if ( .not. DO_SL_ISOTROPIC ) THEN
!           do ibeam = 1, N_SZANGLES
!             TFACTOR = TRANS_ATMOS_FINAL(ibeam)
!             SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IBEAM) = &
!             SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IBEAM) * TFACTOR
!             SLTERM_F_0(0:NMOMS,1:NSTOKES,1:NSTREAMS,IBEAM) = &
!             SLTERM_F_0(0:NMOMS,1:NSTOKES,1:NSTREAMS,IBEAM) * TFACTOR
!             USER_SLTERM_F_0(0:NMOMS,1:NSTOKES,1:N_USER_VZANGLES,IBEAM) = &
!             USER_SLTERM_F_0(0:NMOMS,1:NSTOKES,1:N_USER_VZANGLES,IBEAM) * TFACTOR
!           enddo
!         endif
!  End of Transmittance calculation
!      ENDIF

!  SINGLE-SCATTER & DIRECT-BOUNCE CALCULATIONS
!  ===========================================

!  Version 2.8, Major revision, 7/8/16

!    - SSCORR and DBCORRECTION routines completely removed
!    - Calculations only done using the FO code, Version 1.5 (replaces FO Version 1.4 added 7/2/13)
!    - Rob Fix 3/17/15. Exception handling updated for serious error.
!    - Rob Fix 7/08/16. Master interface updated for FO 1.5, includes FMATRIX inputs.

!  Not required if no solar sources, and no user-angles
!mick fix 9/19/2017 - deactivated the DO_SOLAR_SOURCES IF condition due to the
!                     enhanced nature of the new internal FO code
!mick mod 9/19/2017 - added FO_STOKES_ATMOS & FO_STOKES_SURF (new output from vector FO code)
!                   - using them instead of FO_STOKES_SS & FO_STOKES_DB now (as in LIDORT)

!  4/15/20. Version 2.8.2. Single scatter correction: flux multiplier. Now always F / 4pi

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

      !IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_USER_VZANGLES ) THEN
          IF ( DO_FOCORR ) THEN

!  Call to interface. Updated, 17 September 2016.
!   -- 4/9/19. Added output FO Surface-leaving assignation + cumulative transmittances.
!              Added input, water-leaving control

!  4/15/20. Version 2.8.2.
!      -- Argument list considerably reduced in size
!      -- No geometrical variables (all removed), now precalculated in FOGeometry (input)
!      -- No greekmoms_total_input setting, no DO_FMATRIX and no ngreek_moments_input
!      -- ssflux and na_offset are input directly (not re-set here)
!      -- Exception handling no longer required
!      -- Add DO_DOUBLET_GEOMETRY flag and Offset array (SZD_OFFSET)

            CALL VFO_MASTER_INTERFACE ( &
                DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                       & ! Input Sources flags
                DO_PLANE_PARALLEL, DO_FOCORR_OUTGOING, DO_DELTAM_SCALING,                         & ! Input Model Flags
                DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,                              & ! Input Model Flags
                DO_DOUBLET_GEOMETRY, DO_PARTLAYERS, DO_LAMBERTIAN_SURFACE,                        & ! Input Model Flags
                DO_SURFACE_LEAVING, DO_WATER_LEAVING, DO_SL_ISOTROPIC,                            & ! Input Sources flags
                NSTOKES, NLAYERS, N_GEOMETRIES, N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS,      & ! Input numbers
                VZA_OFFSETS, SZD_OFFSETS, N_USER_LEVELS, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,  & ! Input Offsets/Levels
                N_PARTLAYERS, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,       & ! Input partials
                VLIDORT_FOGeom, HEIGHT_GRID, SS_FLUX_MULTIPLIER, FLUXVEC,                         & ! Input Geometry/Solar
                DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, DELTAU_VERT, FMATRIX_UP, FMATRIX_DN,        & ! Inputs (Optical)
                TRUNC_FACTOR, ALBEDO, EXACTDB_BRDFUNC, SLTERM_ISOTROPIC, SLTERM_USERANGLES,       & ! Inputs (Surface)
                THERMAL_BB_INPUT, SURFACE_BB_INPUT, USER_EMISSIVITY,                              & ! Inputs (Thermal)
                FO_STOKES_SS, FO_STOKES_DB, FO_STOKES_DTA, FO_STOKES_DTS,                         & ! Output Stokes vectors
                FO_STOKES_ATMOS, FO_STOKES_SURF, FO_STOKES, FO_CUMTRANS, FO_SLTERM )                ! Output Stokes vectors

!  Exception handling (no longer present, 4/15/20, Version 2.8.2)

!  4/15/20. Version 2.8.2. Copy FO results directly to Type structure arrays

!mick note 9/19/2017 - Important! STOKES_SS & STOKES_DB contain BOTH solar AND thermal
!                      direct terms when computing in the crossover region!

            if ( do_upwelling ) then
               do v = 1, N_GEOMETRIES
                  VLIDORT_Sup%SS%TS_STOKES_SS(1:N_USER_LEVELS,v,1:NSTOKES,UPIDX) = &
                              FO_STOKES_ATMOS(1:N_USER_LEVELS,v,1:NSTOKES,UPIDX)
                  VLIDORT_Sup%SS%TS_STOKES_DB(1:N_USER_LEVELS,v,1:NSTOKES)       = &
                              FO_STOKES_SURF (1:N_USER_LEVELS,v,1:NSTOKES)
               enddo
            endif
            if ( do_dnwelling ) then
               do v = 1, N_GEOMETRIES
                  VLIDORT_Sup%SS%TS_STOKES_SS(1:N_USER_LEVELS,v,1:NSTOKES,DNIDX) = &
                              FO_STOKES_ATMOS(1:N_USER_LEVELS,v,1:NSTOKES,DNIDX)
               enddo
            endif

          ENDIF

!  End user-angle and solar-source if blocks

        ENDIF
      !ENDIF

!  ####################
!   MAIN FOURIER LOOP
!  ####################

!  Initialise Fourier loop
!  =======================

!  Zeroing the Fourier-Saved Diagnostic. 4/9/19

      FOURIER_SAVED(1:N_SZANGLES) = 0

!  Set Number of Fourier terms (NMOMENTS = Maximum).
!    ( Starting from 0 = Fundamental )

      SAVE_DO_NO_AZIMUTH  = DO_NO_AZIMUTH
      LOCAL_DO_NO_AZIMUTH = DO_NO_AZIMUTH

!  No azimuth dependency for following cases
!    Thermal-only treatment is isotropic, no azimuth

!      IF ( DO_TRANSMITTANCE_ONLY ) THEN
!        LOCAL_DO_NO_AZIMUTH = .TRUE.
!      ENDIF

!mick mod 3/30/2015 - consolidated IF conditions
      IF ( .NOT.DO_SOLAR_SOURCES .OR. DO_MVOUT_ONLY ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

!  set Fourier number (2 for Rayleigh only)

      IF ( LOCAL_DO_NO_AZIMUTH .OR. DO_FOCORR_ALONE ) THEN
        N_FOURIERS = 0
      ELSE
        IF ( DO_RAYLEIGH_ONLY  ) THEN
          N_FOURIERS = 2
        ELSE
          N_FOURIERS = NMOMENTS
        ENDIF
      ENDIF

!  re-set no-azimuth flag

      DO_NO_AZIMUTH = LOCAL_DO_NO_AZIMUTH

!  Initialise BVP telescoping. Important.

      DO_BVTEL_INITIAL = DO_BVP_TELESCOPING

!  Fourier loop
!  ============

!  Initialize

      LOCAL_ITERATION = .TRUE.
      FOURIER         = -1
      TESTCONV        = 0

!  set up solar beam flags. Required in all cases.
!   ---Even required for the thermal.....

      NSOURCES = 1 ; IF ( DO_SOLAR_SOURCES ) NSOURCES = N_SZANGLES
      DO IBEAM = 1, NSOURCES
        BEAM_TESTCONV  ( IBEAM ) = 0
        BEAM_ITERATION ( IBEAM ) = .TRUE.
        DO_MULTIBEAM ( IBEAM,0:MAXFOURIER ) = .TRUE.   !mick eff 3/22/2017
      ENDDO

!  Start Fourier loop
!  ------------------

      DO WHILE ( LOCAL_ITERATION .AND. FOURIER.LT.N_FOURIERS )

!  Fourier counter

         FOURIER = FOURIER + 1

!  Local start of user-defined streams. Should always be 1.
!    No zenith tolerance now.

         LOCAL_UM_START = 1

!  4/28/19. Local media problem and planetary-problem flags
!    -- if the Planetary problem is set, then must have LOCAL_ALBTRN for BOA Unit illumination #2
!       (regardless of the input values of DO_ALBTRN_MEDIA)        
        
         LOCAL_DO_ALBTRN_MEDIA = .false.
         IF ( DO_ALBTRN_MEDIA(1) ) LOCAL_DO_ALBTRN_MEDIA(1) = ( FOURIER == 0 )
         IF ( DO_ALBTRN_MEDIA(2) ) LOCAL_DO_ALBTRN_MEDIA(2) = ( FOURIER == 0 )

         LOCAL_DO_PLANETARY_PROBLEM = .false.
         IF ( DO_PLANETARY_PROBLEM  ) then
            LOCAL_DO_PLANETARY_PROBLEM = ( FOURIER == 0 )
            LOCAL_DO_ALBTRN_MEDIA(2)   = ( FOURIER == 0 )
         ENDIF

!  Main call to VLidort Fourier module
!  -----------------------------------

!  7/8/16. Version 2.8 revision. Cleanup I/O listings, remove "Classical_Solution, Quad_Output"
        
!   Version 2.8.1, Control for TOA/BOA isotropic illumination added, 3/23/19
        
!  4/9/19. Added Inputs,  Water-leaving control, SOLARBEAM_BOATRANS
!          Added Output, TRANS_ATMOS_FINAL (for water-leaving self-consistency)

!  4/28/19 Module for Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!    -- introduced by R. Spurr 4/26/19. Controlled by flags LOCAL_DO_ALBTRN_MEDIA, LOCAL_DO_PLANETARY_PROBLEM
!    -- Associated outputs are TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES   

!  4/15/20. Version 2.8.2. Copy Local BRDF Fourier-component Input (only what you need)

         IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
            BRDF_F_0(:,1:NSTREAMS,1:N_SZANGLES) = &
               VLIDORT_Sup%BRDF%TS_BRDF_F_0(FOURIER,:,1:NSTREAMS,1:N_SZANGLES)
            BRDF_F  (:,1:NSTREAMS,1:NSTREAMS)   = &
               VLIDORT_Sup%BRDF%TS_BRDF_F  (FOURIER,:,1:NSTREAMS,1:NSTREAMS)
            IF ( DO_USER_VZANGLES ) THEN
               USER_BRDF_F_0(:,1:N_USER_VZANGLES,1:N_SZANGLES) = &
                  VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0(FOURIER,:,1:N_USER_VZANGLES,1:N_SZANGLES)
               USER_BRDF_F(:,1:N_USER_VZANGLES,1:NSTREAMS)      = &
                  VLIDORT_Sup%BRDF%TS_USER_BRDF_F(FOURIER,:,1:N_USER_VZANGLES,1:NSTREAMS)
            ENDIF
         ENDIF

!  4/15/20. Version 2.8.2. Copy Local SLEAVE Fourier-component Input (only what you need)

         IF ( DO_SURFACE_LEAVING ) THEN
            SLTERM_F_0(1:NSTOKES,1:NSTREAMS,1:N_SZANGLES) = &
               VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0(FOURIER,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES)
            IF ( DO_USER_VZANGLES ) THEN
               USER_SLTERM_F_0(1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES) = &
                  VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(FOURIER,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES)
            ENDIF
         ENDIF

!  Call to Fourier routine

         CALL VLIDORT_FOURIER ( FOURIER, &
            DO_UPWELLING, DO_DNWELLING, DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY,                  & !Input flags (RT operation)
            DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,           & !Input flags (RT operation)
            DO_LAYER_SCATTERING, DO_SOLUTION_SAVING, DO_BVTEL_INITIAL,                              & !Input flags (performance)
            DO_MSMODE_VLIDORT, DO_MULTIBEAM, LOCAL_DO_ALBTRN_MEDIA, LOCAL_DO_PLANETARY_PROBLEM,     & !Input flags (Beam/Planetary)
            DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_REFLECTED_DIRECTBEAM,    & !Input flags (Surface)
            DO_MSMODE_THERMAL, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_THERMAL_TRANSONLY,      & !Input flags (thermal)
            DO_FOCORR_ALONE, DO_DBCORRECTION, DO_TOA_CONTRIBS, DO_PARTLAYERS, DO_REAL_EIGENSOLVER,  & !Input Bookkeeping
            DO_TOAFLUX, TOAFLUX, DO_BOAFLUX, BOAFLUX,                                               & !Input TOA/BOA Illumination
            DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION, TF_MAXITER, TF_CRITERION,        & !Input Water-leaving control
            NSTOKES, NSTREAMS, NLAYERS, N_SZANGLES, N_USER_VZANGLES, N_USER_LEVELS, N_THERMAL_COEFFS, NMOMENTS, NSTREAMS_2, & !Input
            NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, N_PARTLAYERS, N_DIRECTIONS, TAYLOR_ORDER, & !Input
            N_ALLLAYERS_UP, N_ALLLAYERS_DN, FLUX_FACTOR, FLUXVEC, COS_SZANGLES, SZA_LOCAL_INPUT, SUN_SZA_COSINES, & !Input
            USER_STREAMS, USER_SECANTS, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS,                   & !Input
            MUELLER_INDEX, DMAT, BVP_REGULAR_FLAG, LOCAL_UM_START, WHICH_DIRECTIONS, SOLARBEAM_BOATRANS,          & !Input
            UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & !Input
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DFLUX, SURFACE_BB_INPUT, EMISSIVITY, USER_EMISSIVITY,         & !Input
            Misc, Therm, Mult,                                                                                    & !Input
            ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,  & !Input
            PIMM_11, PIMM_KM, BVTEL_FOURIER, DO_INCLUDE_THERMEMISS, DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS,     & !Output
            STOKES_F, MEANST_DIFFUSE, FLUX_DIFFUSE, DNMEANST_DIRECT, DNFLUX_DIRECT, MS_CONTRIBS_F,                & !Output MAIN
            TRANS_ATMOS_FINAL, TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,                 & !Output 4/26/19
            STATUS_SUB, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )                                                        !Output Status

!  Useful debug
!        if ( fourier.lt.3)write(*,*)FOURIER,TRANS_ATMOS_FINAL(1),STOKES_F(1:2,1,1,1,1),STOKES_F(2,1,1,1,2)
!        if ( FOURIER.eq.0) write(*,*)FOURIER,STOKES_F(1:5,1,1,1,2)
!        if ( FOURIER.eq.0) write(*,*)FO_STOKES_ATMOS(1:5,1,1,2)
!        if ( FOURIER.eq.0) write(*,*)FO_STOKES_SURF(1:5,1,1)

!  error handling

         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            if ( Len(TRACE_3) == 0 ) TRACE_3 = ' ** Called in VLIDORT_RTCalc_Master'
            STATUS_CALCULATION = VLIDORT_SERIOUS
            VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
            VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
         ENDIF

!  4/9/19. FO-DB adjustment for Fourier 0, using self-consistent water-leaving term
!                          Intensity only, Surface-leaving not polarized.         
!  4/15/20. Version 2.8.2. Use Type structure variables directly. ==> Replace array STOKES_DB
!  4/15/20. Version 2.8.2. Add Doublet geometry option

         if ( FOURIER .eq. 0 .and. DO_FOCORR .and. DO_WATER_LEAVING ) then
            O1 = 1
            if ( do_OBSERVATION_GEOMETRY ) THEN
               DO IB = 1, N_SZANGLES
                  SLTERM_LOCAL = FO_SLTERM(O1,ib) * TRANS_ATMOS_FINAL(ib)
                  do uta = 1, n_user_levels
                     CUMSOURCE_DB = FO_CUMTRANS(UTA,IB) * SLTERM_LOCAL
                     VLIDORT_Sup%SS%TS_STOKES_DB(UTA,IB,O1) = VLIDORT_Sup%SS%TS_STOKES_DB(UTA,IB,O1) + CUMSOURCE_DB
                  enddo
               enddo
            else if ( DO_DOUBLET_GEOMETRY ) THEN
               DO IB = 1, N_SZANGLES
                  g1 = SZD_OFFSETS(IB) + 1 ; g2 = SZD_OFFSETS(IB) + N_USER_VZANGLES
                  SLTERM_LOCAL = FO_SLTERM(O1,g1) * TRANS_ATMOS_FINAL(ib)
                  do uta = 1, n_user_levels
                     do g = g1, g2
                        CUMSOURCE_DB = FO_CUMTRANS(UTA,g) * SLTERM_LOCAL
                        VLIDORT_Sup%SS%TS_STOKES_DB(UTA,g,O1) = VLIDORT_Sup%SS%TS_STOKES_DB(UTA,g,O1) + CUMSOURCE_DB
                     enddo
                  enddo
               enddo
            else
               DO IB = 1, N_SZANGLES ; DO UM = 1, N_USER_VZANGLES
                  g1 = VZA_OFFSETS(IB,UM) + 1 ; g2 = VZA_OFFSETS(IB,UM) + n_user_relazms
                  SLTERM_LOCAL = FO_SLTERM(O1,g1) * TRANS_ATMOS_FINAL(ib)
                  do uta = 1, n_user_levels
                     do g = g1, g2
                        CUMSOURCE_DB = FO_CUMTRANS(UTA,g) * SLTERM_LOCAL
                        VLIDORT_Sup%SS%TS_STOKES_DB(UTA,g,O1) = VLIDORT_Sup%SS%TS_STOKES_DB(UTA,g,O1) + CUMSOURCE_DB
                     enddo
                  enddo
               enddo; enddo
            Endif
         endif

!  Output for WLADJUSTED Water-Leaving . Introduced 4/22/19 for Version 2.8.1
      ! Sleave Results need to be modified from their original inputs.
      ! Debug results - USE WITH CAUTION. Note the preliminary zeroing to avoid unassigned arrays.
!  4/15/20. Version 2.8.2. SLTERM_F_0 arrays defined locally, remove M=Fourier index

         if ( FOURIER .eq. 0 .and. DO_WATER_LEAVING .and. DO_WLADJUSTED_OUTPUT ) then
            VLIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC = zero
            VLIDORT_Out%WLOut%TS_WLADJUSTED_DIRECT    = zero
            VLIDORT_Out%WLOut%TS_WLADJUSTED_F_Ords_0  = zero
            VLIDORT_Out%WLOut%TS_WLADJUSTED_F_User_0  = zero
            DO IB = 1, N_SZANGLES
               TFACTOR = TRANS_ATMOS_FINAL(ib) ; O1 = 1
               VLIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC(O1,IB) = SLTERM_ISOTROPIC(O1,IB) * TFACTOR
               VLIDORT_Out%WLOut%TS_WLADJUSTED_F_Ords_0(0,O1,1:NSTREAMS,IB) = TFACTOR * SLTERM_F_0(O1,1:NSTREAMS,IB)
               IF ( DO_USER_VZANGLES ) THEN
                  VLIDORT_Out%WLOut%TS_WLADJUSTED_DIRECT(O1,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IB) = &
                                          TFACTOR * SLTERM_USERANGLES(O1,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IB)
                  VLIDORT_Out%WLOut%TS_WLADJUSTED_F_User_0(0,O1,1:N_USER_VZANGLES,IB) = &
                                          TFACTOR * USER_SLTERM_F_0(O1,1:N_USER_VZANGLES,IB)
               ENDIF
            ENDDO
         ENDIF

! #################### INNOVATIONS 5/22/20 ##############
!  4/29/20. Special Fourier-component output for Rayleigh + Planetary-problem TOA Upwelling situations
!         IF ( SPECIAL_RAYF_OUTPUT ) then
!           UTA = 1
!           IF ( DO_OBSERVATION_GEOMETRY ) THEN
!             DO IB =1, NBEAMS ; DO O1 = 1, NSTOKES ; LUM = 1
!               VLIDORT_Out%Main%TS_TOAUP_RAYSTOKES_FOURIER(LUM,IB,O1,FOURIER) =  STOKES_F(UTA,LUM,IB,O1,UPIDX)
!             ENDDO ; ENDDO
!           ELSE
!             DO IB =1, NBEAMS ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
!               VLIDORT_Out%Main%TS_TOAUP_RAYSTOKES_FOURIER(UM,IB,O1,FOURIER) =  STOKES_F(UTA,UM,IB,O1,UPIDX)
!             ENDDO ; ENDDO ; ENDDO
!           ENDIF
!         ENDIF
! #################### INNOVATIONS 5/22/20 ##############

!  Fourier summation and Convergence examination
!  ---------------------------------------------

!mick fix 3/30/2015 - added if condition and moved azimuth block from
!                     before call to VLIDORT_FOURIER to here

!  Begin convergence if block

        IF ( .NOT.DO_MVOUT_ONLY ) THEN

!  azimuth cosine/sine factors
!    - Use of adjusted geometries retired for Version 2.8
!mick eff 3/22/2017 - both IF and ELSE sections
!    4/15/20. Version 2.8.2. Add Doublet geometry option

          IF ( FOURIER .GT. 0 ) THEN
            DFC = DBLE(FOURIER)
            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UA = 1, LOCAL_N_USERAZM
                AZM_ARGUMENT = USER_RELAZMS(UA) * DFC * DEG_TO_RAD
                AZMFAC(1:N_USER_VZANGLES,1:N_SZANGLES,UA,1)  = COS(AZM_ARGUMENT)
                AZMFAC(1:N_USER_VZANGLES,1:N_SZANGLES,UA,3)  = SIN(AZM_ARGUMENT)
                AZMFAC(1:N_USER_VZANGLES,1:N_SZANGLES,UA,2)  = AZMFAC(1:N_USER_VZANGLES,1:N_SZANGLES,UA,1)
                AZMFAC(1:N_USER_VZANGLES,1:N_SZANGLES,UA,4)  = AZMFAC(1:N_USER_VZANGLES,1:N_SZANGLES,UA,3)
              ENDDO
            ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
              DO UM = 1, N_USER_VZANGLES
                AZM_ARGUMENT = USER_RELAZMS(UM) * DFC * DEG_TO_RAD
                AZMFAC(UM,1:N_SZANGLES,LUA,1)  = COS(AZM_ARGUMENT)
                AZMFAC(UM,1:N_SZANGLES,LUA,3)  = SIN(AZM_ARGUMENT)
                AZMFAC(UM,1:N_SZANGLES,LUA,2)  = AZMFAC(UM,1:N_SZANGLES,LUA,1)
                AZMFAC(UM,1:N_SZANGLES,LUA,4)  = AZMFAC(UM,1:N_SZANGLES,LUA,3)
              ENDDO
            ELSE
              DO IB = 1, NSOURCES
                AZM_ARGUMENT = USER_RELAZMS(IB) * DFC * DEG_TO_RAD
                AZMFAC(LUM,IB,LUA,1)  = COS(AZM_ARGUMENT)
                AZMFAC(LUM,IB,LUA,3)  = SIN(AZM_ARGUMENT)
                AZMFAC(LUM,IB,LUA,2)  = AZMFAC(LUM,IB,LUA,1)
                AZMFAC(LUM,IB,LUA,4)  = AZMFAC(LUM,IB,LUA,3)
              ENDDO
            ENDIF
          ENDIF

!   -- only done for beams which are still not converged
!      This is controlled by flag DO_MULTIBEAM

!   -- new criterion, SS is added for Fourier = 0, as this means that
!      higher-order terms will be relatively smaller, which implies
!      faster convergence in some circumstances (generally not often).

          IBEAM_COUNT = 0
          DO IBEAM = 1, NSOURCES
            IF ( DO_MULTIBEAM ( IBEAM, FOURIER ) ) THEN

!  Convergence and radiance summation. Version 2.8, Clean up I/O listings, 7/8/16
!  4/15/20. Version 2.8.2. Use Type structure variables directly. 
!                                ==> Replace arrays STOKES, STOKES_SS, STOKES_DB
!  4/15/20. Version 2.8.2. Add Doublet geometry option (new convergence routine)

              IF ( DO_OBSERVATION_GEOMETRY ) THEN
                CALL VLIDORT_CONVERGE_OBSGEO ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_DBCORRECTION, DO_UPWELLING, DO_NO_AZIMUTH, & ! Input flags
                  DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST, DO_TOA_CONTRIBS,    & ! Input flags
                  NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, IBEAM, FOURIER,                & ! Input numbers
                  N_CONVTESTS, VLIDORT_ACCURACY, AZMFAC, N_DIRECTIONS, WHICH_DIRECTIONS,    & ! Input Bookkeep, Conv.
                  STOKES_F, MS_CONTRIBS_F, VLIDORT_Sup%SS, VLIDORT_Out%Main,                & ! Input/Output fields
                  FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )                ! Output diagnostics
              ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
                CALL VLIDORT_CONVERGE_DOUBLET ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_DBCORRECTION, DO_UPWELLING,                          & ! Input flags
                  DO_NO_AZIMUTH, DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST,                & ! Input flags
                  NSTOKES, NSTREAMS, NLAYERS, N_OUT_STREAMS, N_USER_LEVELS, IBEAM, FOURIER,           & ! Input numbers
                  N_CONVTESTS, VLIDORT_ACCURACY, SZD_OFFSETS, AZMFAC, N_DIRECTIONS, WHICH_DIRECTIONS, & ! Input Bookkeep, Conv.
                  STOKES_F, VLIDORT_Sup%SS, VLIDORT_Out%Main,                                         & ! Input/Output fields
                  FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM)  )                         ! Output Convergence
              ELSE
                CALL VLIDORT_CONVERGE ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_DBCORRECTION, DO_UPWELLING, DO_NO_AZIMUTH, & ! Input flags
                  DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST, DO_TOA_CONTRIBS,    & ! Input flags
                  NSTOKES, NSTREAMS, NLAYERS, N_OUT_STREAMS, N_USER_RELAZMS, N_USER_LEVELS, & ! Input control numbers
                  IBEAM, FOURIER, N_CONVTESTS, VLIDORT_ACCURACY, VZA_OFFSETS, AZMFAC,       & ! Input numbers, convergence 
                  N_DIRECTIONS, WHICH_DIRECTIONS, LOCAL_UM_START, LOCAL_N_USERAZM,          & ! Input bookkeeping
                  STOKES_F, MS_CONTRIBS_F, VLIDORT_Sup%SS, VLIDORT_Out%Main,                & ! Input and output fields
                  FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )                ! Output diagnostics
              ENDIF

!  Check number of beams already converged

              IF ( BEAM_ITERATION(IBEAM) ) THEN
                IBEAM_COUNT = IBEAM_COUNT + 1
              ELSE
                DO L = FOURIER+1,MAXFOURIER
                  DO_MULTIBEAM (IBEAM,L) = .FALSE.
                ENDDO
              ENDIF

!  end beam count loop

            ENDIF
          END DO

!  If all beams have converged, stop iteration

          IF ( IBEAM_COUNT .EQ. 0 ) LOCAL_ITERATION = .FALSE.

!  End convergence if block

        END IF

!  Fourier output
!  --------------

!  Disabled, Version 2.8
!  Open file if Fourier = 0
!  Write Standard Fourier output
!  Close file if iteration has finished
!  New comment:
!    If the SS correction is set, Fourier=0 will include SS field
!        IF ( DO_WRITE_FOURIER ) THEN
!          FUNIT = VLIDORT_FUNIT
!          IF ( FOURIER .EQ. 0 ) THEN
!            OPEN(FUNIT,FILE=FOURIER_WRITE_FILENAME,STATUS='UNKNOWN')
!          ENDIF
!          CALL VLIDORT_WRITEFOURIER ( FUNIT, FOURIER )
!          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
!        ENDIF

!  4/15/20. Version 2.8.2. 
!  Copying for Sleave Fourier components must be done here
!  Additional Control for Externalized input (SLEAVE). Introduced 4/22/19 for Version 2.8.1
!  Sleave Results May have been modified from their original inputs.
!    Allows you to come out with modified SLEAVE, so you can use it!!!!!!!!!!!!

         IF ( DO_EXTERNAL_WLEAVE ) THEN
            O1 = 1
            IF ( FOURIER .eq. 0 ) THEN
               VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(O1,1:N_SZANGLES) = SLTERM_ISOTROPIC(O1,1:N_SZANGLES)
               IF ( DO_USER_VZANGLES ) THEN
                  VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(O1,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
                                        SLTERM_USERANGLES(O1,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
               ENDIF
            ENDIF
            VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0(FOURIER,O1,1:NSTREAMS,1:N_SZANGLES) = SLTERM_F_0(O1,1:NSTREAMS,1:N_SZANGLES)
            VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(FOURIER,O1,1:N_USER_VZANGLES,1:N_SZANGLES) = &
                                  USER_SLTERM_F_0(O1,1:N_USER_VZANGLES,1:N_SZANGLES)
         ENDIF

!  End Fourier loop

      ENDDO

!  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

!  Major result output
!  ===================

!  Revised Version 2.8, 7/8/16. Cleaned up I/O listings. 

!  Standard output
!  ---------------

      IF ( DO_WRITE_RESULTS ) THEN
        RUNIT = VLIDORT_RESUNIT
        OPEN(RUNIT,FILE=RESULTS_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL VLIDORT_WRITERESULTS ( &
          RUNIT, DO_FULLRAD_MODE, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_DOUBLE_CONVTEST,    &
          DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,      &
          NSTOKES, VLIDORT_ACCURACY, SZANGLES, N_USER_RELAZMS, USER_RELAZMS,                  &
          N_USER_LEVELS, USER_LEVELS, HEIGHT_GRID, DELTAU_VERT_INPUT, DO_NO_AZIMUTH,          &
          N_SZANGLES, N_DIRECTIONS, WHICH_DIRECTIONS, N_OUT_STREAMS, OUT_ANGLES,              &
          PARTLAYERS_OUTFLAG, VZA_OFFSETS, TAUGRID_INPUT, DO_MULTIBEAM,                       &
          VLIDORT_Out%Main%TS_STOKES, MEANST_DIFFUSE, FLUX_DIFFUSE, FOURIER_SAVED )
        CLOSE(RUNIT)
      ENDIF

!  Geophysical input (scenario) write
!  ----------------------------------

!  4/15/20. Version 2.8.2. NGREEK_MOMENTS_INPUT Dropped, thermal emission also

      IF ( DO_WRITE_SCENARIO ) THEN
        SUNIT = VLIDORT_SCENUNIT
        OPEN(SUNIT,FILE=SCENARIO_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL VLIDORT_WRITESCEN ( &
          SUNIT, DO_DELTAM_SCALING, NSTREAMS, NLAYERS, &
          N_SZANGLES, SZANGLES, N_USER_RELAZMS, USER_RELAZMS,                     &
          N_USER_VZANGLES, USER_VZANGLES, N_USER_LEVELS, USER_LEVELS,             &
          OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, DO_LAMBERTIAN_SURFACE,         &
          ALBEDO, DO_SURFACE_EMISSION, SURFACE_BB_INPUT,                          &
          QUAD_STREAMS, QUAD_WEIGHTS, QUAD_ANGLES, DO_NO_AZIMUTH, NMOMENTS,       &
          GREEKMAT_INDEX, TAUGRID_INPUT, OMEGA_TOTAL, GREEKMAT_TOTAL, TAUGRID )
        CLOSE(SUNIT)
      ENDIF

!  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

!  ========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (IN/OUT variables)
!  ========================================================

!  FOCORR Booleans reorganized for Version 2.8
!mick mod 9/19/2017 - reordered FO variables to conform to newly modified input type structure
!                   - DO_FOCORR_ALONE now defined internally
      !VLIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE        = DO_FOCORR_ALONE

      VLIDORT_ModIn%MBool%TS_DO_FOCORR              = DO_FOCORR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL     = DO_FOCORR_EXTERNAL
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR        = DO_FOCORR_NADIR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING     = DO_FOCORR_OUTGOING

!  4/15/20. Version 2.8.2. Remove DO_SSCORR_TRUNCATION altogether.
!  4/15/20. Version 2.8.2. Remove DO_SSCORR_USEFMAT altogether.
!      VLIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION   = DO_SSCORR_TRUNCATION
!      VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT      = DO_SSCORR_USEFMAT

!  Additional Control for Externalized input (SLEAVE). Introduced 3/18/19 for Version 2.8.1

      VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE     = DO_EXTERNAL_WLEAVE
      
!  Solar control

      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = DO_SOLAR_SOURCES
      VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = DO_REFRACTIVE_GEOMETRY

!  4/15/20. Version 2.8.2. DO_CHAPMAN function flag no longer required.
!     - Chapman functions now calculcated in Geometry routine
!      VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = DO_CHAPMAN_FUNCTION

!  Performance control

      !VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY      = DO_ISOTROPIC_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH          = DO_NO_AZIMUTH
      !VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER         = DO_ALL_FOURIER

      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY       = DO_RAYLEIGH_ONLY
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING      = DO_DELTAM_SCALING
      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST     = DO_DOUBLE_CONVTEST
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING     = DO_SOLUTION_SAVING
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING     = DO_BVP_TELESCOPING

!  RT Model control
!  4/15/20. Version 2.8.2. Copy back DO_DOUBLET_GEOMETRY

      VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES       = DO_USER_VZANGLES
      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = DO_ADDITIONAL_MVOUT
      VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY          = DO_MVOUT_ONLY
      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY   = DO_THERMAL_TRANSONLY
      VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY= DO_OBSERVATION_GEOMETRY
      VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY    = DO_DOUBLET_GEOMETRY

!  Modified control inputs
!  4/15/20. Version 2.8.2. MAXMOMENTS_INPUT not required in LIDORT_Main.
!      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT   = NGREEK_MOMENTS_INPUT

!  Modified beam inputs

      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES             = N_SZANGLES
      VLIDORT_ModIn%MSunRays%TS_SZANGLES(1:N_SZANGLES) = SZANGLES(1:N_SZANGLES)

!  Modified user value inputs

      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS      = N_USER_RELAZMS
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS)           = USER_RELAZMS(1:N_USER_RELAZMS)

!mick fix 9/19/2017 - added IF condition
      IF ( DO_USER_VZANGLES ) THEN
        VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES      = N_USER_VZANGLES
        VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:N_USER_VZANGLES) = USER_VZANGLES(1:N_USER_VZANGLES)
      ENDIF

!  4/15/20. Version 2.8.2. GEOMETRY_SPECHEIGHT moved to MCont type-structure (from MUserVal)

      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS) = USER_LEVELS(1:N_USER_LEVELS)
      VLIDORT_ModIn%MCont%TS_GEOMETRY_SPECHEIGHT             = GEOMETRY_SPECHEIGHT

      IF ( DO_OBSERVATION_GEOMETRY ) THEN
         VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS     = N_USER_OBSGEOMS
         VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,:) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,:)
      ENDIF

!  Modified Chapman function inputs
!  4/15/20. Version 2.8.2. No Longer required, as these quantities are calculated in Geometry module.
!    VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS = CHAPMAN_FACTORS
!      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS    = EARTH_RADIUS

!  Modified optical inputs, 
!  4/15/20. Version 2.8.2. OMEGA_TOTAL_INPUT now Intent(in), so this statement noLonger required.
!      VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS) = OMEGA_TOTAL_INPUT(1:NLAYERS)

!  ==========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (pure OUT variables)
!  ==========================================================

!  Stokes vectors and fluxes
!    Output direct Mean/Flux added 17 May 2012

!  4/15/20. Version 2.8.2. Stokes vectors - Type structure filled directly
!            VLIDORT_Out%Main%TS_STOKES(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)    = &
!                                STOKES(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)

!  4/15/20. Version 2.8.2. Integrated values - Only  copy if necessary.

      IF ( DO_MVOUT_ONLY .or. DO_ADDITIONAL_MVOUT ) THEN
         IF ( DO_UPWELLING ) THEN
            VLIDORT_Out%Main%TS_MEANST_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,UPIDX) = &
                                MEANST_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,UPIDX)
            VLIDORT_Out%Main%TS_FLUX_DIFFUSE  (1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,UPIDX) = &
                                FLUX_DIFFUSE  (1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,UPIDX)
         ENDIF
         IF ( DO_DNWELLING ) THEN
            VLIDORT_Out%Main%TS_MEANST_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,DNIDX) = &
                                MEANST_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,DNIDX)
            VLIDORT_Out%Main%TS_FLUX_DIFFUSE  (1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,DNIDX) = &
                                FLUX_DIFFUSE  (1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,DNIDX)

            VLIDORT_Out%Main%TS_DNMEANST_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = &
                                DNMEANST_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)
            VLIDORT_Out%Main%TS_DNFLUX_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = &
                                DNFLUX_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)
         ENDIF
      ENDIF

!  Verify TF calculation. Debug 2/3/16. Worked out.
!      do ibeam = 1, n_szangles
!         write(*,*)'MK2 Result = ',FLUX_DIFFUSE(2,1,1,2)
!         write(778,'(a,i2,1p3e15.6,2x)')'Original', ibeam,MEANST_DIFFUSE(2,ibeam,1,2),FLUX_DIFFUSE(2,ibeam,1,2),&
!         DNFLUX_DIRECT(2,ibeam,1)
!      enddo

!  4/26/19. Media properties output.
!   4/15/20. Version 2.8.2. Only fill out what you need

      IF ( DO_ALBTRN_MEDIA(1) ) THEN
         VLIDORT_Out%Main%TS_ALBMED_USER(1:NSTOKES,1:N_USER_VZANGLES) = ALBMED_USER(1:NSTOKES,1:N_USER_VZANGLES)
         VLIDORT_Out%Main%TS_ALBMED_FLUXES(1:NSTOKES,1:2) = ALBMED_FLUXES(1:NSTOKES,1:2)
      ENDIF
      IF ( DO_ALBTRN_MEDIA(2) ) THEN
         VLIDORT_Out%Main%TS_TRNMED_USER(1:NSTOKES,1:N_USER_VZANGLES) = TRNMED_USER(1:NSTOKES,1:N_USER_VZANGLES)
         VLIDORT_Out%Main%TS_TRNMED_FLUXES(1:NSTOKES,1:2) = TRNMED_FLUXES(1:NSTOKES,1:2)
      ENDIF

!  4/28/19. Planetary problem output
!  ---------------------------------

!  5/22/20. Version 2.8.2 Upgrades 
!    ==> Revised code to use only the first index of TRANSBEAM ==> makes the Q-problem valid

!  Here is the Old code (pre 5/22/20 Upgrade) 
!      IF ( DO_PLANETARY_PROBLEM ) THEN
!         VLIDORT_Out%Main%TS_PLANETARY_SBTERM = TRNMED_FLUXES(1,2)
!         IF ( DO_OBSERVATION_GEOMETRY ) THEN
!            DO IB =1, N_GEOMETRIES ; DO O1 = 1, NSTOKES
!               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,IB) = TRANSBEAM(O1,IB) * TRNMED_USER(O1,IB) / PIE
!            ENDDO ; ENDDO
!         ELSE
!            DO IB =1, N_SZANGLES ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
!               OFF = VZA_OFFSETS(IB,UM) ; TRANS = TRANSBEAM(O1,IB) * TRNMED_USER(O1,UM) / PIE
!               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,OFF+1:OFF+N_USER_RELAZMS) = TRANS
!            ENDDO ; ENDDO ; ENDDO
!         ENDIF   
!      ENDIF
      
!  Here is the New code (5/22/20 Upgrade)  

      IF ( DO_PLANETARY_PROBLEM ) THEN
         VLIDORT_Out%Main%TS_PLANETARY_SBTERM = TRNMED_FLUXES(1,2)
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            DO IB =1, N_GEOMETRIES ; DO O1 = 1, NSTOKES
               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,IB) = TRANSBEAM(1,IB) * TRNMED_USER(O1,IB) / PIE
            ENDDO ; ENDDO
         ELSE
            DO IB =1, N_SZANGLES ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
               OFF = VZA_OFFSETS(IB,UM) ; TRANS = TRANSBEAM(1,IB) * TRNMED_USER(O1,UM) / PIE
               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,OFF+1:OFF+N_USER_RELAZMS) = TRANS
            ENDDO ; ENDDO ; ENDDO
         ENDIF   
      ENDIF

!  new 12 March 2012
!   IF SS results already available, no need to copy them !
!   4/15/20. Version 2.8.2. No Longer required, as Type structure filled directly/
!      IF ( .NOT. DO_FOCORR_EXTERNAL ) THEN
!         VLIDORT_Sup%SS%TS_STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:) = &
!           STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)
!         VLIDORT_Sup%SS%TS_STOKES_DB(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES) = &
!           STOKES_DB(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
!      ENDIF
      
!  Bookkeeping

      VLIDORT_Out%Main%TS_FOURIER_SAVED(1:N_SZANGLES) = FOURIER_SAVED(1:N_SZANGLES)
      VLIDORT_Out%Main%TS_N_GEOMETRIES                = N_GEOMETRIES

!  Solar Beam Transmittance to BOA
!  rob fix 11/17/2014, for diagnostic use only

      VLIDORT_Out%Main%TS_SOLARBEAM_BOATRANS (1:N_SZANGLES) = SOLARBEAM_BOATRANS (1:N_SZANGLES)

!  Exception handling

      VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  = STATUS_INPUTCHECK
      VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION

      VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
      VLIDORT_Out%Status%TS_CHECKMESSAGES(0:NCHECKMESSAGES)  = CHECKMESSAGES(0:NCHECKMESSAGES)
      VLIDORT_Out%Status%TS_ACTIONS(0:NCHECKMESSAGES)        = ACTIONS(0:NCHECKMESSAGES)

      VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
      VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
      VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
      VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3

!  ===================================
!  END COPY LOCAL VARIABLES TO OUTPUTS
!  ===================================

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE VLIDORT_DEBUG_INPUT_MASTER()

!  3/28/14. Changes for Version 2.7. remove LTE linearization references
!  7/8/16. Removed DO_QUAD_OUTPUT for Version 2.8. I/O list needs a clean-up
!  3/1/17 . Changes for Version 2.8. Complete reorganization of argument lists
!mick mod 9/19/2017 - added DO_FLUORESCENCE, DO_TF_ITERATION, TAYLOR_ORDER, FMATRIX_UP, FMATRIX_DN,
!                           ATMOS_WAVELENGTH, TF_MAXITER, TF_CRITERION to argument list

!  Additional Control for SLEAVE (DO_WLADJUSTED_OUTPUT,DO_EXTERNAL_WLEAVE).
!     Introduced 3/18/19 for Version 2.8.1
         
!  4/26/19. Record the Media-problem inputs     (DO_ALBTRN_MEDIA)
!  4/28/19. Record the Planetary-problem inputs (DO_PLANETARY_PROBLEM)
!  3/23/19, Version 2.8.1, Record Control for TOA/BOA illumination
        
!  4/15/20. Version 2.8.2. NMOMS, NMOMENTS_INPUT --> NMOMENTS
!  4/15/20. Version 2.8.2. Remove DO_SSCORR_TRUNCATION, DO_SSCORR_USEFMAT
!  4/15/20. Version 2.8.2. Add DO_DOUBLET flag to argument list

      CALL VLIDORT_WRITE_STD_INPUT ( DO_DEBUG_WRITE, &
        DO_SOLAR_SOURCES,      DO_THERMAL_EMISSION,  DO_THERMAL_TRANSONLY, DO_SURFACE_EMISSION,      & ! Sources
        DO_FULLRAD_MODE,       DO_FOCORR,            DO_FOCORR_EXTERNAL,   DO_FOCORR_NADIR,          & ! Model
        DO_FOCORR_OUTGOING,    DO_RAYLEIGH_ONLY,     DO_DELTAM_SCALING,    DO_DOUBLE_CONVTEST,       & ! Model
        DO_SOLUTION_SAVING,    DO_BVP_TELESCOPING,   DO_UPWELLING,         DO_DNWELLING,             & ! Model
        DO_PLANE_PARALLEL,     DO_CHAPMAN_FUNCTION,  DO_USER_VZANGLES,     DO_OBSERVATION_GEOMETRY,  & ! Model
        DO_DOUBLET_GEOMETRY,   DO_ADDITIONAL_MVOUT,  DO_MVOUT_ONLY,        DO_REFRACTIVE_GEOMETRY,   & ! RT model
        DO_ALBTRN_MEDIA,       DO_PLANETARY_PROBLEM, DO_TOAFLUX,           DO_BOAFLUX,               & ! Media/Planetary/Flux
        DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING,   DO_SL_ISOTROPIC,      DO_EXTERNAL_WLEAVE,       & ! Surface
        DO_WATER_LEAVING,      DO_FLUORESCENCE,      DO_WLADJUSTED_OUTPUT, DO_TF_ITERATION,          & ! Surface/Media-props
        TAYLOR_ORDER, TF_MAXITER, NSTOKES, NSTREAMS, NLAYERS, NFINELAYERS, N_THERMAL_COEFFS,      & ! Main numbers
        DO_TOA_CONTRIBS,  DO_SPECIALIST_OPTION_1,  DO_SPECIALIST_OPTION_2,                        & ! Specialist Options
        DO_SPECIALIST_OPTION_3, NLAYERS_NOMS, NLAYERS_CUTOFF,                                     & ! Specialist Options
        N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES, N_USER_OBSGEOMS, N_USER_LEVELS,              & ! Geometry and level numbers
        SZANGLES, USER_RELAZMS, USER_VZANGLES, USER_OBSGEOMS, USER_LEVELS,                        & ! Geometry and level control
        VLIDORT_ACCURACY, FLUX_FACTOR, EARTH_RADIUS, RFINDEX_PARAMETER, TF_CRITERION,             & ! Flux/Acc/Radius
        HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, FINEGRID, GEOMETRY_SPECHEIGHT,              & ! Grids
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                               & ! Optical
        FMATRIX_UP, FMATRIX_DN, ALBEDO, TOAFLUX, BOAFLUX,                                         & ! Optical & albedo & Fluxes
        THERMAL_BB_INPUT, SURFACE_BB_INPUT, ATMOS_WAVELENGTH,                                     & ! Thermal & spectral
        DO_WRITE_INPUT, DO_WRITE_SCENARIO, DO_WRITE_FOURIER, DO_WRITE_RESULTS,                         & ! Debug control
        INPUT_WRITE_FILENAME, SCENARIO_WRITE_FILENAME, FOURIER_WRITE_FILENAME, RESULTS_WRITE_FILENAME )  ! Debug control

!  4/15/20. Version 2.8.2. Remove NMOMENTS, Only printing First Fourier component.

      IF (.NOT. DO_LAMBERTIAN_SURFACE) THEN
!mick mod 9/19/2017 - added NMOMENTS to input
        CALL VLIDORT_WRITE_SUP_BRDF_INPUT ( &
          DO_USER_VZANGLES, DO_SURFACE_EMISSION, &
          NSTOKES,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,NSTREAMS,&
          EXACTDB_BRDFUNC,BRDF_F_0,BRDF_F,USER_BRDF_F_0,USER_BRDF_F,&
          EMISSIVITY,USER_EMISSIVITY)
      END IF

      IF (DO_FOCORR_EXTERNAL) THEN
        CALL VLIDORT_WRITE_SUP_SS_INPUT ( &
          NSTOKES,N_USER_LEVELS,STOKES_SS,STOKES_DB)
      END IF

!  4/15/20. Version 2.8.2. Remove NMOMENTS, Only printing First Fourier component.

      IF (DO_SURFACE_LEAVING) THEN
        CALL VLIDORT_WRITE_SUP_SLEAVE_INPUT ( &
          DO_USER_VZANGLES, &
          NSTOKES,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,NSTREAMS,&
          SLTERM_ISOTROPIC,SLTERM_USERANGLES,SLTERM_F_0,USER_SLTERM_F_0)
      END IF

      END SUBROUTINE VLIDORT_DEBUG_INPUT_MASTER

      END SUBROUTINE VLIDORT_RTCalc_Master

!

      SUBROUTINE VLIDORT_FOURIER ( FOURIER, &        
            DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY,                   & !Input flags (RT operation)
            DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,           & !Input flags (RT operation)
            DO_LAYER_SCATTERING, DO_SOLUTION_SAVING, DO_BVTEL_INITIAL,                              & !Input flags (performance)
            DO_MSMODE_VLIDORT, DO_MULTIBEAM, DO_ALBTRN_MEDIA, DO_PLANETARY_PROBLEM,                 & !Input flags (Beam/Planetary)
            DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_REFLECTED_DIRECTBEAM,    & !Input flags (Surface)
            DO_MSMODE_THERMAL, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_THERMAL_TRANSONLY,      & !Input flags (thermal)
            DO_FOCORR_ALONE, DO_DBCORRECTION, DO_TOA_CONTRIBS, DO_PARTLAYERS, DO_REAL_EIGENSOLVER,  & !Input Bookkeeping
            DO_TOAFLUX, TOAFLUX, DO_BOAFLUX, BOAFLUX,                                               & !Input TOA/BOA illumination
            DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION, TF_MAXITER, TF_CRITERION,        & !Input Water-leaving control
            NSTOKES, NSTREAMS, NLAYERS, N_SZANGLES, N_USER_STREAMS, N_USER_LEVELS, N_THERMAL_COEFFS, NMOMENTS, NSTREAMS_2, & !Input
            NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, N_PARTLAYERS, N_DIRECTIONS, TAYLOR_ORDER, & !Input
            N_ALLLAYERS_UP, N_ALLLAYERS_DN, FLUX_FACTOR, FLUXVEC, COS_SZANGLES, SZA_LOCAL_INPUT, SUN_SZA_COSINES, & !Input
            USER_STREAMS, USER_SECANTS, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS,                   & !Input
            MUELLER_INDEX, DMAT, BVP_REGULAR_FLAG, LOCAL_UM_START, WHICH_DIRECTIONS, SOLARBEAM_BOATRANS,          & !Input
            UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & !Input
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DFLUX, SURFBB, EMISSIVITY, USER_EMISSIVITY,                   & !Input
            Misc, Therm, Mult,                                                                                    & !Input
            ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,  & !Input
            PIMM_11, PIMM_KM, BVTEL_FOURIER, DO_INCLUDE_THERMEMISS, DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS,     & !Output
            STOKES_F, MEANST_DIFFUSE, FLUX_DIFFUSE, DNMEANST_DIRECT, DNFLUX_DIRECT, MS_CONTRIBS_F,                & !output MAIN
            TRANS_ATMOS_FINAL, TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,                 & !Output 4/26/19
            STATUS, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )                                                            !Output Status

!  This line was removed, 7/8/16.
!        DO_SPECIALIST_OPTION_2, DO_DEBUG_WRITE, DO_FDTEST,                                         & !Input

!   Version 2.8.1, Control for TOA/BOA isotropic illumination added, 3/23/19
!      --- NOT THE SAME AS THE ALBTRN_MEDIA condition.         
        
!  Argument list revised for Version 2.8.1, 4/9/19, 3/23/19
!    1.    Water-leaving adjustment control added (DO_WATER_LEAVING, DO_TF_ITERATION, TF_MAXITER, TF_CRITERION)
!    2.    TOA/BOA isotropic illumination control, DO_TOAFLUX, TOAFLUX, DO_BOAFLUX, BOAFLUX    
!    3.    TRANS_ATMOS_FINAL, added to the output list (Self-adjusting water-leaving formulation)
!    4.    SOLARBEAM_BOATRANS, TRANS_SOLAR_BEAM distinguished.
      
!  Argument list revised for Version 2.8.1, 4/26/19, 4/28/19
!  4/26/19 Module for Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!    -- introduced by R. Spurr 4/26/19. Controlled by flags DO_ALBTRN_MEDIA, DO_PLANETARY_PROBLEM
     
!  Complete Fourier component calculation for the Standard Code.
!    Consolidated V2p4RTC: Lambertian or BRDF Surface

      USE VLIDORT_PARS_m

      USE VLIDORT_Work_def_m

      USE VLIDORT_MISCSETUPS_m , Only : VLIDORT_DIRECTRADIANCE, VLIDORT_PIMATRIX_SETUP_OMP
      USE VLIDORT_THERMALSUP_m , Only : THERMAL_CLSOLUTION, THERMAL_STERMS_UP, THERMAL_STERMS_DN
      USE VLIDORT_MULTIPLIERS_m, Only : HMULT_MASTER
      USE VLIDORT_SOLUTIONS_m
      USE VLIDORT_BVPROBLEM_m  , Only : BVP_MATRIXSETUP_MASTER,    BVP_SOLUTION_MASTER, &
                                        BVPTEL_MATRIXSETUP_MASTER, BVPTEL_SOLUTION_MASTER
      USE VLIDORT_INTENSITY_m  , Only : VLIDORT_UPUSER_INTENSITY, VLIDORT_DNUSER_INTENSITY, VLIDORT_INTEGRATED_OUTPUT
      USE VLIDORT_UNPACK_m

!  4/26/19. Media-properties routine, Mark II, 4/28/19.
   
      USE VLIDORT_mediaprops_m

!  Implicit none

      IMPLICIT NONE

!  INPUT, INTENT(IN) ONLY
!  ======================

!  Input Fourier component number

      INTEGER, INTENT (IN) ::          FOURIER

!  Flags
!  -----

!  SS control

!      LOGICAL, INTENT (IN) ::          DO_SSCORR_NADIR  ! Dropped 2.8
      LOGICAL, INTENT (IN) ::          DO_FOCORR_ALONE   ! Renamed 2.8

!  Solar control

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY

!  RT control

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (IN) ::          DO_MVOUT_ONLY
      LOGICAL, INTENT (IN) ::          DO_SOLUTION_SAVING
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT

!  Convergence tracker

      LOGICAL, INTENT (IN) ::          DO_MULTIBEAM ( MAXBEAMS, 0:MAXFOURIER )

!  4/26/19. Control flags for the isotropic-illuminated Media calculations

      LOGICAL, intent(in)  ::          DO_ALBTRN_MEDIA(2)
      
!  4/28/19  Added control for the planetary problem      

      LOGICAL, intent(in)  ::          DO_PLANETARY_PROBLEM
      
!  4/9/19. Additional Water-leaving control

      LOGICAL, intent(in)  ::          DO_WATER_LEAVING
      LOGICAL, intent(in)  ::          DO_EXTERNAL_WLEAVE
      LOGICAL, intent(in)  ::          DO_TF_ITERATION
      INTEGER, INTENT (IN) ::          TF_MAXITER
      DOUBLE PRECISION, INTENT (IN) :: TF_CRITERION

!  surface

      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
!      LOGICAL, INTENT (IN) ::          DO_DIRECT_BEAM
      LOGICAL, INTENT (IN) ::          DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::          DO_SURFACE_EMISSION

!  surface leaving

      LOGICAL, INTENT (IN) ::          DO_SURFACE_LEAVING
      LOGICAL, INTENT (IN) ::          DO_SL_ISOTROPIC

!  Output control

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS

!  thermal

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_THERMAL
      LOGICAL, INTENT (IN) ::          DO_THERMAL_EMISSION

!  Bookkeeping

      LOGICAL, INTENT (IN) ::          DO_REAL_EIGENSOLVER  ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING  ( 0:MAXMOMENTS, MAXLAYERS )

!  Other flags

      LOGICAL, INTENT (IN) ::          DO_TOA_CONTRIBS

!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19

      LOGICAL         , INTENT(IN) :: DO_TOAFLUX, DO_BOAFLUX 
      DOUBLE PRECISION, INTENT(IN) :: TOAFLUX, BOAFLUX

!  Disabled flags, removed 7/8/16 for Version 2.8

      !LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
      !LOGICAL, INTENT (IN) ::          DO_FO_CALC !New 02 Jul 2013
      !LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      !LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      !LOGICAL, INTENT (IN) ::          DO_BVP_TELESCOPING
!      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL
!      LOGICAL, INTENT (IN) ::          DO_FDTEST
!      LOGICAL, INTENT (IN) ::          DO_DEBUG_WRITE
!      LOGICAL, INTENT (IN) ::          DO_SPECIALIST_OPTION_2
!      LOGICAL, INTENT (IN) ::          DO_SPECIALIST_OPTION_3
!      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT             ! removed, Version 2.8 7/8/16
!      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION      ! removed, Version 2.8 7/8/16

!  Integers
!  --------

!  Basic

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

      INTEGER, INTENT (IN) ::          N_SZANGLES
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS

!  derived

      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL

      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_DIRECTIONS

      INTEGER, INTENT (IN) ::          N_ALLLAYERS_UP
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_DN

!  Version 2p7 input, 2/19/14
      INTEGER, INTENT (IN) ::          TAYLOR_ORDER

!  Disabled

      !INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
      !INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      !INTEGER, INTENT (IN) ::          NFINELAYERS
      !INTEGER, INTENT (IN) ::          N_GEOMETRIES
      !INTEGER, INTENT (IN) ::          NLAYERS_CUTOFF

!  Floating point and Misc.
!  ------------------------

!  FLux control

      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )

!  SZA control

      DOUBLE PRECISION, INTENT (IN) :: COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SUN_SZA_COSINES ( MAXLAYERS, MAX_SZANGLES )

!  User-angle control

      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )

!  Quadrature (discrete ordinates)

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )

!  solution bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DMAT ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DFLUX ( MAXSTOKES )
      LOGICAL, INTENT (IN) ::          BVP_REGULAR_FLAG ( 0:MAXMOMENTS )

      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      
!  Rob fix 11/27/14. Proxy for new output.
!     NOTE. SOLARBEAM_BOATRANS is computed with Unscaled optical depths, TRANS_SOLAR_BEAM with scaled ODs.

      DOUBLE PRECISION, INTENT (IN) :: SOLARBEAM_BOATRANS ( MAXBEAMS )

!  Level output control

      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )

!  Disabled

      !DOUBLE PRECISION, INTENT (IN) :: SS_FLUX_MULTIPLIER
      !DOUBLE PRECISION, INTENT (IN) :: SZANGLES ( MAX_SZANGLES )
      !DOUBLE PRECISION, INTENT (IN) :: SZANGLES_ADJUST ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      !DOUBLE PRECISION, INTENT (IN) :: SIN_SZANGLES ( MAX_SZANGLES )
      !DOUBLE PRECISION, INTENT (IN) :: USER_VZANGLES ( MAX_USER_VZANGLES )
      !DOUBLE PRECISION, INTENT (IN) :: USER_VZANGLES_ADJUST ( MAX_USER_VZANGLES )
      !INTEGER, INTENT (IN) ::          VZA_OFFSETS  ( MAX_SZANGLES, MAX_USER_VZANGLES )
      !DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      !DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS_ADJUST ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      !DOUBLE PRECISION, INTENT (IN) :: USER_LEVELS   ( MAX_USER_LEVELS )
      !INTEGER, INTENT (IN) ::          LAYER_MAXMOMENTS ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: PARTLAYERS_VALUES ( MAX_PARTLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: CHAPMAN_FACTORS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: EARTH_RADIUS
      !DOUBLE PRECISION, INTENT (IN) :: HEIGHT_GRID ( 0:MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      !DOUBLE PRECISION, INTENT (IN) :: TAUGRID_INPUT ( 0:MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  surface

      DOUBLE PRECISION, INTENT (IN) :: SURFBB

!  Work type structures

      TYPE(VLIDORT_Work_Miscellanous), INTENT (IN) :: Misc
      TYPE(VLIDORT_Work_Thermal),      INTENT (IN) :: Therm
      TYPE(VLIDORT_Work_Multiplier),   INTENT (IN) :: Mult

!  From BRDF supplement

      DOUBLE PRECISION, INTENT (IN) :: ALBEDO

!  4/15/20. Version 2.8.2.  Local Fourier-component dimension (0:MAXMOMENTS) has been dropped
!                                 Now using only the terms you want.

      DOUBLE PRECISION, INTENT (IN) :: BRDF_F        ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F_0      ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )

      DOUBLE PRECISION, INTENT (IN) :: EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  New surface-leaving stuff 17 May 2012

      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )

!  4/15/20. Version 2.8.2.  Local Fourier-component dimension (0:MAXMOMENTS) has been dropped
!                                 Now using only the terms you want.

      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_F_0 ( MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   USER_SLTERM_F_0 ( MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  OUTPUTS, Intent(INOUT)
!  ======================

!  Flags

      LOGICAL, INTENT (INOUT) ::          DO_BVTEL_INITIAL
      INTEGER, INTENT (INOUT) ::          BVTEL_FOURIER

!  Help arrays for the PI MAtrix Setup to make if thread-safe

      DOUBLE PRECISION, INTENT (INOUT) :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, INTENT (INOUT) :: PIMM_KM ( MAX_ALLSTRMS_P1 )

!  Mean-value results will be done for Fourier = 0, then input again.
 
      DOUBLE PRECISION, INTENT (INOUT) :: MEANST_DIFFUSE ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_DIFFUSE   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: DNMEANST_DIRECT( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: DNFLUX_DIRECT  ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  First Introduced 3/22/17 for LIDORT, based on VLIDORT code 

      DOUBLE PRECISION, INTENT (INOUT) :: TRANS_ATMOS_FINAL  ( MAX_SZANGLES )

!  4/26/19. Special Media-property output. -- Introduced by R. Spurr.
!     ** Output for User-angles and fluxes. Output of Beam transmittance for Planetary problem.

      REAL(fpk), INTENT (INOUT) :: ALBMED_USER ( MAXSTOKES, MAX_USER_STREAMS ), ALBMED_FLUXES( MAXSTOKES, 2 ) !  TOA illumination
      REAL(fpk), INTENT (INOUT) :: TRNMED_USER ( MAXSTOKES, MAX_USER_STREAMS ), TRNMED_FLUXES( MAXSTOKES, 2 ) !  BOA illumination

!  4/28/19. Special Output of Beam transmittance for Planetary problem.

      real(fpk), intent(inout)  :: TRANSBEAM   ( MAXSTOKES, MAXBEAMS )

!  Disabled

      !DOUBLE PRECISION, INTENT (INOUT) :: TAUGRID ( 0:MAXLAYERS )
      !DOUBLE PRECISION, INTENT (INOUT) :: OMEGA_TOTAL ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (INOUT) :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

      !DOUBLE PRECISION, INTENT (INOUT) :: STOKES_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      !DOUBLE PRECISION, INTENT (INOUT) :: STOKES_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      !DOUBLE PRECISION, INTENT (INOUT) :: SS_CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  OUTPUT, Intent(OUT)
!  ===================

!  Inclusion flags

      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_THERMEMISS

!  Fourier components

      DOUBLE PRECISION, INTENT (OUT) ::  MS_CONTRIBS_F ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  STOKES_F ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
                                                    MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

!  Local variables
!  ================

!  From VLIDORT_MISCSETUPS
!  -----------------------

!  deltam-scaled optical
!mick fix 9/19/2017 - added LEVELS_SOLARTRANS & PARTIALS_SOLARTRANS to facilitate correction of direct flux

      DOUBLE PRECISION :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Transmittances

      DOUBLE PRECISION :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_UTUP_DISORDS ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UTDN_DISORDS ( MAXSTREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

      DOUBLE PRECISION :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )

!  Beam pseudo-spherical

      INTEGER ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

      DOUBLE PRECISION :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )

!  Disabled

      !DOUBLE PRECISION :: TRUNC_FACTOR ( MAXLAYERS )
      !DOUBLE PRECISION :: FAC1 ( MAXLAYERS )
      !DOUBLE PRECISION :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      !DOUBLE PRECISION :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  From THERMAL_SETUP
!  ------------------

!  Bookkeeping

      DOUBLE PRECISION :: THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: TCOM1        ( MAXLAYERS, MAX_THERMAL_COEFFS )

!  Driect thermal solutions

      DOUBLE PRECISION :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_UT_DIRECT_UP  ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UT_DIRECT_DN  ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  From EMULT_MASTER
!  -----------------

!  Multipliers

      DOUBLE PRECISION :: EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  Disabled

      !LOGICAL ::          EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      !DOUBLE PRECISION :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      !DOUBLE PRECISION :: SIGMA_P  ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Other local variables
!  ---------------------
!mick mod 9/19/2017 - added and activated the variable DO_LOCALBEAM as in LIDORT

!  Local direct beam reflectance

      LOGICAL ::          DO_LOCALBEAM ( MAXBEAMS )

!  Indices and help for the Planetary problem

      INTEGER ::          I, IBEAM, IPARTIC, LAYER, N, O1, K, K0, K1, K2, KO1, LUI
      DOUBLE PRECISION :: TRANSQUAD(MAXSTREAMS,MAXSTOKES), TRANSDIRECT
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, LXR, MXR, LXR_CR, LXR_CI, HOM1CR,  MXR_CR

!  Local inclusion flags
!   Version 2.8.1, Control for TOA/BOA  illumination added, 3/23/19

      LOGICAL ::          DO_INCLUDE_MVOUTPUT
      LOGICAL ::          DO_INCLUDE_DIRECTRF   ! 4/9/19 renamed
      LOGICAL ::          DO_INCLUDE_DIRECTSL   ! 4/9/19 New
!      LOGICAL ::          DO_INCLUDE_DIRECTBEAM  ! Replaced 4/9/19
      LOGICAL ::          DO_INCLUDE_TOAFLUX
      LOGICAL ::          DO_INCLUDE_BOAFLUX

!  Local variables for the iterated water-leaving

      DOUBLE PRECISION :: TFACTOR, SL
      
!  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION ::  FLUX_MULTIPLIER
      DOUBLE PRECISION ::  DELTA_FACTOR
      DOUBLE PRECISION ::  SURFACE_FACTOR

!  Error tracing

!      LOGICAL ::           FAIL
      INTEGER ::           STATUS_SUB
      CHARACTER (LEN=2) :: CF

!  SOLUTION VARIABLES
!  ==================

!  Pi matrix setups
!  ----------------

!  At quadrature angles

      DOUBLE PRECISION :: PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQM ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION :: PI_XQM_POST ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQM_PRE  ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQP_PRE  ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  At user angles

      DOUBLE PRECISION :: PI_XUP ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XUM ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION :: PI_XUM_POST ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XUP_PRE  ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  At solar angles

      DOUBLE PRECISION :: PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Thermal solutions of the RTE
!  ----------------------------

!  Discrete ordinate solutions

      DOUBLE PRECISION :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS )

!  User-solutions

      DOUBLE PRECISION :: U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS )

!  Post processed solutions

      DOUBLE PRECISION :: LAYER_TSUP_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: LAYER_TSUP_DN ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Homogeneous RTE
!  ---------------

!  Eigenmatrices

      DOUBLE PRECISION :: SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )

      DOUBLE PRECISION :: EIGENMAT_SAVE ( MAXEVALUES, MAXEVALUES, MAXLAYERS )

!  Eigensolutions

      DOUBLE PRECISION :: REAL_KSQ   ( MAXSTRMSTKS )
      DOUBLE PRECISION :: IMAG_KSQ   ( MAXSTRMSTKS )
      DOUBLE PRECISION :: KEIGEN_CSQ ( MAXEVALUES )
      DOUBLE PRECISION :: LEFT_EVEC  ( MAXSTRMSTKS, MAXSTRMSTKS )
      DOUBLE PRECISION :: RITE_EVEC  ( MAXSTRMSTKS, MAXSTRMSTKS )

      LOGICAL ::          EIGENDEGEN  ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER ::          EIGENMASK_R ( MAXEVALUES )
      INTEGER ::          EIGENMASK_C ( MAXEVALUES )

!  auxiliary solutions

      DOUBLE PRECISION :: FWD_SUMVEC  ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: FWD_DIFVEC  ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

!  Main solution variables

      INTEGER ::          K_REAL     ( MAXLAYERS )
      INTEGER ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION :: KEIGEN     ( MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Eigenstream transmittances

      DOUBLE PRECISION :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  User homogeneous solutions

      DOUBLE PRECISION :: HELPSTOKES ( 0:MAXMOMENTS, MAXEVALUES, MAXSTOKES )
      DOUBLE PRECISION :: ZETA_M ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: ZETA_P ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Homogeneous solution multipliers

      LOGICAL ::          HSINGO  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION :: UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Solar beam particular integral
!  ------------------------------

!  associated arrays for solving 

      DOUBLE PRECISION :: QSUMVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: QDIFVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: QVEC_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION :: QDIF_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION :: QMAT_SAVE ( MAXSTRMSTKS, MAXSTRMSTKS )
      INTEGER ::          QPIVOT ( MAXSTRMSTKS )

!  Quadrature Solutions

      DOUBLE PRECISION :: BVEC   ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  User solutions

      DOUBLE PRECISION :: HELPSTOKES_BEAM ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  BVProblem arrays
!  ----------------

!  Regular BVP

      DOUBLE PRECISION :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER ::          IPIVOT   ( MAXTOTAL )

!  One-layer BVP

      DOUBLE PRECISION :: SMAT2    ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER ::          SIPIVOT  ( MAXSTRMSTKS_2 )

!  Column vectors [mick fix 7/29/2014 - added to make VLIDORT threadsafe]

      DOUBLE PRECISION :: COL2    ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION :: SCOL2   ( MAXSTRMSTKS_2, MAXBEAMS )

!  Integration constants (BVP results)

      DOUBLE PRECISION :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Telescoped BVP setup

      INTEGER ::          ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER ::          N_BVTELMATRIX_SIZE
      INTEGER ::          N_BVTELMATRIX_SUPDIAG
      INTEGER ::          N_BVTELMATRIX_SUBDIAG
      INTEGER ::          NLAYERS_TEL

!  Telescoped BVP

      DOUBLE PRECISION :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER ::          IPIVOTTEL   ( MAXTOTAL )
      DOUBLE PRECISION :: COLTEL2     ( MAXTOTAL, MAXBEAMS ) ! [mick fix 7/29/2014 - added to make VLIDORT threadsafe]

!  Surface-reflected solutions
!  ---------------------------

!  Thermal BOA term

      DOUBLE PRECISION :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  Direct beam solutions, 4/9/19 renamed

      DOUBLE PRECISION :: ATMOS_ATTN ( MAXBEAMS )
      DOUBLE PRECISION :: RF_DIRECT_BEAM      ( MAXSTREAMS,       MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION :: RF_USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
!      DOUBLE PRECISION :: DIRECT_BEAM      ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
!      DOUBLE PRECISION :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  4/9/19 Surface-leaving contributions, added
      
      DOUBLE PRECISION :: SL_QUADTERM ( MAXSTREAMS,       MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION :: SL_USERTERM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      
!  BRDF auxiliary array

      DOUBLE PRECISION :: AXBID_F  ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!  reflected downwelling solutions

      DOUBLE PRECISION :: STOKES_DOWNSURF ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: R2_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: R2_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: R2_BEAM ( MAXSTREAMS, MAXSTOKES )

!  Cumulative discrete-ordinate transmittances to surface (Telescoped problem only)

      DOUBLE PRECISION :: CUMTRANSDOM(MAXSTREAMS)
      DOUBLE PRECISION :: CUMQUADDOM (MAXSTREAMS)

!  Cumulative source terms
!  -----------------------

!  Post-processed

      DOUBLE PRECISION :: CUMSOURCE_UP ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION :: CUMSOURCE_DN ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )

!  ##############
!  initialization
!  ##############

!  module status and message initialization

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      TRACE_3 = ' '

!  Exponent
      
      LUI = 1

!  Set local flags
!  ---------------

!  inclusion of thermal surface emission term, only for Fourier = 0

      DO_INCLUDE_SURFEMISS = .FALSE.
      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_SURFEMISS = .TRUE.
        ENDIF
      ENDIF

!  inclusion of thermal emission term, only for Fourier = 0

      DO_INCLUDE_THERMEMISS = .FALSE.
      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_THERMEMISS = .TRUE.
        ENDIF
      ENDIF

!  inclusion of TOA  illumination, only for Fourier = 0
!     Version 2.8.1, Control for TOA  illumination added, 3/23/19

      DO_INCLUDE_TOAFLUX = .FALSE.
      IF ( DO_TOAFLUX ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_TOAFLUX = .TRUE.
        ENDIF
      ENDIF

!  inclusion of BOA  illumination, only for Fourier = 0
!     Version 2.8.1, Control for BOA  illumination added, 3/23/19

      DO_INCLUDE_BOAFLUX = .FALSE.
      IF ( DO_BOAFLUX ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_BOAFLUX = .TRUE.
        ENDIF
      ENDIF

!  Surface flag (for inclusion of some kind of reflecting boundary)
!    Zero albedo case relaxed, if there is surface leaving. 12/17/15

      DO_INCLUDE_SURFACE = .TRUE.
      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        IF ( FOURIER .NE. 0 ) THEN
            DO_INCLUDE_SURFACE = .FALSE.
        ELSE
          IF ( ALBEDO .EQ. ZERO .and..not.DO_SURFACE_LEAVING ) THEN
            DO_INCLUDE_SURFACE = .FALSE.
          ENDIF
        ENDIF
      ENDIF

!  surface reflectance factors

      IF ( FOURIER .EQ. 0 ) THEN
        SURFACE_FACTOR = TWO
        DELTA_FACTOR   = ONE
      ELSE
        SURFACE_FACTOR = ONE
        DELTA_FACTOR   = TWO
      ENDIF

!  Flux multiplier

!       = 1 / 4.pi with beam sources,  = 1 for Thermal alone.    ! OLD
!       = 1                                                      ! NEW

!      FLUX_MULTIPLIER   = SS_FLUX_MULTIPLIER * DELTA_FACTOR     ! OLD
      FLUX_MULTIPLIER = DELTA_FACTOR                             ! NEW

      IF ( DO_INCLUDE_THERMEMISS.AND..NOT.DO_SOLAR_SOURCES) THEN
        FLUX_MULTIPLIER = DELTA_FACTOR
      ENDIF

!  inclusion of mean value output

      DO_INCLUDE_MVOUTPUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_MVOUTPUT = .TRUE.
        ENDIF
      ENDIF

!  ##################
!  mick chg 7/11/2014 - moved misc setup, thermal setup, emult setup, sscor codes (2),
!                       db code, and fo code up to VLIDORT_MASTER.  There, the needed
!                       vars are packed into type structures to bypass f90 continuation
!                       line limits and passed down.  Here, they are unpacked.
!                     - (Later) Only the misc-setups, thermal-setups, emult setups are UNPACKED HERE
!  (Rev. 7/8/16, 2.8) - Cleaned up

      CALL VLIDORT_UNPACK_MISC ( Misc,                                          & ! Input structure to be unpacked
        NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS, N_USER_STREAMS, N_SZANGLES, NMOMENTS, & ! Input Numbers
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT,                                 & ! unpacked Optical
        LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, OMEGA_GREEK,                    & ! unpacked Optical
        BEAM_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,                 & ! unpacked Solar
        T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                         & ! unpacked Trans. D.O.
        T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,   & ! unpacked Trans solar/User
        CUMTRANS, INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA )                     ! unpacked Av. secant.

      IF ( DO_THERMAL_EMISSION ) THEN
        CALL VLIDORT_UNPACK_THERM ( Therm,                          & !Input Sstructure to be unpacked
          NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,  & ! Input Numbers
          THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Output
          T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN )  ! Output
      END IF

!  Flag changed 3/1/17.
!        IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN

      IF ( DO_SOLAR_SOURCES .AND. DO_USER_STREAMS ) THEN
        IF ( .NOT. DO_FOCORR_ALONE ) THEN
          CALL VLIDORT_UNPACK_MULT ( Mult,                  & ! Input structure to be unpacked
            NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_STREAMS,  & ! Input numbers
            EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )    ! Output
        ENDIF
      ENDIF

!  ##################

!  Direct beam flag (only if above albedo flag has been set)
!mick fix 7/11/2014 - this IF block moved from above and modified to accomodate move of
!                     "Fourier = 0" subroutines
!mick mod 9/19/2017 - commented out this IF block (DO_REFLECTED_DIRECTBEAM already defined
!                     in VLIDORT_MASTER & passed down VIA the "Misc" type structure)
!                   - activated DO_LOCALBEAM as in LIDORT

      !IF ( FOURIER .GT. 0 ) THEN
      !  IF ( DO_DIRECT_BEAM .AND. DO_SOLAR_SOURCES ) THEN
      !    IF ( DO_INCLUDE_SURFACE ) THEN
      !      DO IBEAM = 1, N_SZANGLES
      !        DO_REFLECTED_DIRECTBEAM(IBEAM) = .TRUE.
      !      ENDDO
      !    ELSE
      !      DO IBEAM = 1, N_SZANGLES
      !        DO_REFLECTED_DIRECTBEAM(IBEAM) = .FALSE.
      !      ENDDO
      !    ENDIF
      !  ENDIF
      !ENDIF

      IF ( DO_SOLAR_SOURCES .and. DO_INCLUDE_SURFACE ) THEN
        DO_LOCALBEAM(1:N_SZANGLES) = DO_REFLECTED_DIRECTBEAM(1:N_SZANGLES)
      ELSE
        DO_LOCALBEAM(1:N_SZANGLES) = .FALSE.
      ENDIF

!  Surface direct beam (Not required if no solar sources)
!    - New Surface-Leaving arguments added, 17 May 2012
!    - Revised I/O listings for Version 2.8, 7/8/16
!mick mod 9/19/2017 - replaced DO_REFLECTED_DIRECTBEAM with DO_LOCALBEAM

!  4/9/19. Argument list refined to include Water-leaving control
!  4/9/19. Here, SL output ONLY for non water-leaving, or water-leaving external, otherwise zeroed

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( .not. DO_FOCORR_ALONE ) THEN
          CALL VLIDORT_DIRECTRADIANCE ( &
                DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_REFRACTIVE_GEOMETRY,          & ! Input flags (general)
                DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING,             & ! Input flags (surface)
                DO_SL_ISOTROPIC, DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE,                     & ! input
                FOURIER, NSTOKES, NSTREAMS, NLAYERS, N_SZANGLES, N_USER_STREAMS,               & ! Input numbers
                MUELLER_INDEX, FLUX_FACTOR, FLUXVEC, DELTA_FACTOR,                         & ! Input flux-factors, indices
                SZA_LOCAL_INPUT, COS_SZANGLES, TRANS_SOLAR_BEAM, DO_LOCALBEAM,             & ! Input solar beam
                ALBEDO, BRDF_F_0, USER_BRDF_F_0,                                           & ! Input surface reflectance
                SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                             & ! Input surface leaving
                ATMOS_ATTN, RF_DIRECT_BEAM, RF_USER_DIRECT_BEAM, SL_QUADTERM, SL_USERTERM )  ! Output reflected direct beams.
        ENDIF
     ENDIF

!  Useful debug
!      do IBEAM = 1, N_SZANGLES
!         if ( FOURIER.eq.0 ) write(*,*)IBEAM,ATMOS_ATTN(IBEAM),RF_DIRECT_BEAM(1,IBEAM,1), RF_USER_DIRECT_BEAM(1,IBEAM,1)
!         if ( FOURIER.eq.0 ) write(*,*)IBEAM,SL_QUADTERM(1,IBEAM,1), SL_USERTERM(1,IBEAM,1)
!      enddo

!  ###################
!  Spherical functions
!  ###################

!  Get Pi Matrices for this Fourier componentIncluding all beam angles !!
!   - Version 2.7. Use "OMP" routine for thread-safety. Mick fix.
!   - Version 2.8. Clean up arguments 7/8/16
!mick fix 3/2/2020 - added DO_SOLAR_SOURCES to VLIDORT_PIMATRIX_SETUP_OMP inputs

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
        CALL VLIDORT_PIMATRIX_SETUP_OMP ( &
          DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_USER_STREAMS, FOURIER, & ! Input flags and Fourier
          NSTOKES, NSTREAMS, NLAYERS, N_SZANGLES, N_USER_STREAMS,             & ! Input numbers
          QUAD_STREAMS, USER_STREAMS, COS_SZANGLES, SUN_SZA_COSINES,          & ! Input streams and angles
          NMOMENTS, MUELLER_INDEX, DMAT,                                      & ! auxiliary inputs
          PIMM_11, PIMM_KM, PI_XQP, PI_XQM, PI_XUP, PI_XUM, PI_X0P,           & ! Output Pi-Matrix stuff
          PI_XQM_POST, PI_XQM_PRE, PI_XQP_PRE, PI_XUM_POST, PI_XUP_PRE )        ! Output Pi-Matrix stuff
      ENDIF

!  ########################################
!  RT differential equation Eigensolutions
!  ########################################

!  Version 2.8, remove GOTO statements
!  Go to continuation point for thermal transmittance only
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 8899

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Start layer loop

        DO LAYER = 1, NLAYERS

!  Get Discrete ordinate solutions for this layer
!   Version 2.8, Cleaned-up I/O presentation. 7/8/16.

          CALL VLIDORT_QHOM_SOLUTION ( &
            DO_SOLUTION_SAVING, LAYER, FOURIER,                                   & ! Input Flag, Fourier/Layer
            NSTOKES, NSTREAMS, N_USER_LEVELS, NMOMENTS, NSTKS_NSTRMS,             & ! Input numbers
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,         & ! Input partial control
            QUAD_STREAMS, QUAD_HALFWTS, DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, & ! Quadrature and bookkeeping
            DELTAU_VERT, PARTAU_VERT, OMEGA_GREEK, PI_XQP, PI_XQM, PI_XQM_PRE,    & ! Input optical and PI matrices
            T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                       & ! Input transmittances (discrete Ords.)
            SAB, DAB, EIGENMAT_SAVE, REAL_KSQ, IMAG_KSQ,                          & ! Output Eigenproblem
            LEFT_EVEC, RITE_EVEC, EIGENDEGEN, EIGENMASK_R, EIGENMASK_C,           & ! Output Eigenproblem
            FWD_SUMVEC, FWD_DIFVEC, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,     & ! Output Homogeneous solutions
            K_REAL, K_COMPLEX, KEIGEN, KEIGEN_CSQ, SOLA_XPOS, SOLB_XNEG,          & ! Output Homogeneous solutions
            STATUS_SUB, MESSAGE, TRACE_1 )                                          ! Exception handling

!  .. error tracing

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER
            TRACE_2 = 'Called in VLIDORT_FOURIER, Fourier component '//CF
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

!  Get Post-processing ("user") solutions for this layer
!   Version 2.8, Cleaned-up I/O presentation. 7/8/16.

          IF  ( STERM_LAYERMASK_UP(LAYER) .OR. STERM_LAYERMASK_DN(LAYER) ) THEN
            IF ( DO_USER_STREAMS ) THEN
              CALL VLIDORT_UHOM_SOLUTION ( &
                DO_UPWELLING, DO_DNWELLING,                                  & ! Input flags
                LAYER, FOURIER, NSTOKES, NSTREAMS, N_USER_STREAMS,           & ! Input numbers
                NMOMENTS, LOCAL_UM_START, DO_LAYER_SCATTERING,               & ! Input bookkeeping
                QUAD_HALFWTS, USER_SECANTS, OMEGA_GREEK,                     & ! Input quadratures + optical
                PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE, PI_XUM_POST, PI_XUP_PRE, & ! Input PI Matrices
                K_REAL, K_COMPLEX, KEIGEN, KEIGEN_CSQ, SOLA_XPOS, SOLB_XNEG, & ! Input Homog. RTE solutions
                UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                  & ! Output user Homog
                HSINGO, HELPSTOKES, ZETA_M, ZETA_P )                           ! Output user Homog
            ENDIF
          ENDIF

!  end layer loop

        ENDDO

!  Prepare homogeneous solution multipliers
!    - mick fix 3/30/2015 - added if DO_USER_STREAMS condition
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.

        IF ( DO_USER_STREAMS ) THEN
          CALL HMULT_MASTER ( &
            DO_UPWELLING, DO_DNWELLING,                                   & ! flags
            NLAYERS, N_USER_STREAMS, N_USER_LEVELS, TAYLOR_ORDER,         & ! Basic numbers
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                       & ! whole-layer control
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! partial-layer control
            LOCAL_UM_START, USER_SECANTS, DELTAU_VERT, PARTAU_VERT,       & ! secants, optical
            T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                     & ! User-stream transmittances
            T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                     & ! eigenstreamTransmittances
            K_REAL, K_COMPLEX, HSINGO, ZETA_M, ZETA_P,                    & ! RTE Eigen solution
            HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD ) ! OUTPUT Multipliers
        END IF

!  ############################################
!   boundary value problem - MATRIX PREPARATION
!  ############################################

!  standard case using compression of band matrices, etc..
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.

        IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN

          CALL BVP_MATRIXSETUP_MASTER ( &
            DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,       & ! Input Flags
            FOURIER, NSTOKES, NSTREAMS, NLAYERS,             & ! Input Numbers
            NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,        & ! Input Numbers
            NTOTAL, N_SUBDIAG, N_SUPDIAG, SURFACE_FACTOR,    & ! Input Numbers
            QUAD_STRMWTS, ALBEDO, BRDF_F,                    & ! Input Surface stuff
            MUELLER_INDEX, K_REAL, K_COMPLEX,                & ! Input Bookkeeping
            SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,              & ! Input RTE stuff
            R2_HOMP, R2_HOMM, AXBID_F,                       & ! Output Surface reflection
            BANDMAT2, IPIVOT, SMAT2, SIPIVOT,                & ! Output BVP Matrices
            STATUS_SUB, MESSAGE, TRACE_1 )                     ! Exception handling

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER
            TRACE_2 = 'Error from BVP_MATRIXSETUP_MASTER, '// 'Called in VLIDORT_FOURIER, Fourier # '//CF
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

!  Telescoped case

        ELSE

          CALL BVPTEL_MATRIXSETUP_MASTER ( &
            DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,          & ! Input Flags
            DO_LAYER_SCATTERING, FOURIER, NSTOKES, NSTREAMS,    & ! Input Numbers
            NLAYERS, NSTREAMS_2, NSTKS_NSTRMS_2, NSTKS_NSTRMS,  & ! Input Numbers
            DO_BVTEL_INITIAL, BVTEL_FOURIER,                    & ! BVP Tel Control
            N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,     & ! BVP Tel Control
            N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,       & ! BVP Tel Control
            SURFACE_FACTOR, MUELLER_INDEX, K_REAL, K_COMPLEX,   & ! Input Bookkeeping
            QUAD_STRMWTS, ALBEDO, BRDF_F,            & ! Input Surface inputs
            SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, T_DELT_DISORDS, & ! Input RTE stuff
            R2_HOMP, R2_HOMM, AXBID_F, CUMTRANSDOM, CUMQUADDOM, & ! Output Surface reflection
            BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT,             & ! Output BVP Matrices
            STATUS_SUB, MESSAGE, TRACE_1 )                        ! Exception handling

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER
            TRACE_2 = 'Error from BVPTEL_MATRIXSETUP_MASTER, '//'Called in VLIDORT_FOURIER, Fourier # '//CF
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

        ENDIF

!  4/26/19. Add Call to Media properties subroutine.
!   -- Stand-alone output, but this is a necessary call for the planetary problem

        IF ( DO_ALBTRN_MEDIA(1) .or. DO_ALBTRN_MEDIA(2) .or. DO_PLANETARY_PROBLEM ) THEN
          CALL VLIDORT_MediaProps &
          ( DO_USER_STREAMS, DO_ALBTRN_MEDIA, NSTOKES,                       & ! Input
            NLAYERS, NSTREAMS, N_USER_STREAMS, NSTKS_NSTRMS, NSTKS_NSTRMS_2, & ! Input
            NTOTAL, N_SUBDIAG, N_SUPDIAG, QUAD_STRMWTS, DELTAU_VERT,         & ! Input
            K_REAL, K_COMPLEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG,           & ! Input Homog solutions
            USER_STREAMS, T_DELT_USERM, UHOM_UPDN, UHOM_UPUP,                & ! Input User solutions
            HMULT_1, HMULT_2, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,              & ! Input Multipliers, BVP
            ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,          & ! output
            STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                            ! Output
          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
             TRACE_3 = 'Error from VLIDORT_MediaProps, '//'Called in VLIDORT_FOURIER'
             STATUS = VLIDORT_SERIOUS ; RETURN
          ENDIF
        ENDIF
      
!  End scattering calculation (replaces GOTO 8899)

      ENDIF

!  Continuation point for avoiding the scattering calculations
! 8899 continue. Removed, Version 2.8

!  ################
!  Thermal Solution
!  ################

!  1. Find the Particular solution (also for transmittance only)
!  2. Compute thermal layer source terms. (Upwelling and Downwelling)
!     These will be scaled up by factor 4.pi if solar beams as well

!    ** REMARK. No Green's function treatment here
!    ** Version 2.8, Cleaned-up I/O presentations, all thermal routines. 7/8/16.

      IF ( DO_INCLUDE_THERMEMISS ) THEN

!  discrete ordinate particular integral for thermal sources.
!    - Rob fix 1/31/11, This line was wrong (see below)     PI_XQP, PI_XQM, PI_XUP, PI_XQM_PRE, &
!mick fix 9/19/2017 - added DO_USER_STREAMS to input

        CALL THERMAL_CLSOLUTION ( &
          DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                 & ! input flags
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_THERMAL_TRANSONLY,    & ! Input flags
          NSTREAMS, NLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,         & ! Input basic numbers
          NMOMENTS, NSTREAMS_2, LOCAL_UM_START, N_PARTLAYERS,          & ! Input numbers
          PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Input level control
          QUAD_STREAMS, QUAD_HALFWTS, OMEGA_GREEK, SAB, DAB,           & ! Input optical and SAB/DAB
          T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,              & ! Input Discrete Ord. Trans.
          PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE,                          & ! Input PI matrices, Corrected 1/31/11
          THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,                & ! Input thermal setups
          T_WUPPER, T_WLOWER, UT_T_PARTIC,                             & ! Output thermal solutions
          U_TPOS1, U_TNEG1, U_TPOS2, U_TNEG2,                          & ! Output User thermal solutions
          STATUS_SUB, MESSAGE, TRACE_1 )                                 ! Exception handling

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
          TRACE_2 = 'Called in VLIDORT_FOURIER, Fourier 0'
          STATUS = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  User solutions and post processing.
!    - mick fix 3/30/2015 - modified if conditions for up & dn thermal source terms

        IF ( DO_UPWELLING .AND. DO_USER_STREAMS) THEN
          CALL THERMAL_STERMS_UP ( &
            DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, DO_PARTLAYERS,    & ! Input flags
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,  & ! Input numbers
            N_ALLLAYERS_UP, PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP,  & ! Input level control
            LOCAL_UM_START, USER_STREAMS, T_DELT_USERM, T_UTUP_USERM, & ! Input User streams and transmittances
            DELTAU_POWER, XTAU_POWER, U_TPOS1, U_TPOS2,               & ! Input Thermal setups/solutions
            T_DIRECT_UP, T_UT_DIRECT_UP,                              & ! Input thermal direct solutions
            LAYER_TSUP_UP, LAYER_TSUP_UTUP )                            ! Output user RTE thermal
        ENDIF

        IF ( DO_DNWELLING .AND. DO_USER_STREAMS ) THEN
          CALL THERMAL_STERMS_DN ( &
            DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, DO_PARTLAYERS,    & ! Input flags
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,  & ! Input numbers
            N_ALLLAYERS_DN, PARTLAYERS_LAYERIDX, STERM_LAYERMASK_DN,  & ! Input level control
            LOCAL_UM_START, USER_STREAMS, T_DELT_USERM, T_UTDN_USERM, & ! Input User streams and transmittances
            DELTAU_POWER, XTAU_POWER, U_TNEG1, U_TNEG2,               & ! Input Thermal setups/solutions
            T_DIRECT_DN, T_UT_DIRECT_DN,                              & ! Input thermal direct solutions
            LAYER_TSUP_DN, LAYER_TSUP_UTDN )                            ! Output user RTE thermal
        ENDIF

!  End include thermal emission

      ENDIF

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Skip the thermal-only section if there are solar sources
!      IF ( DO_SOLAR_SOURCES ) GO TO 455

!  Version 2.8, Get rid of the GO TO 455 statement

      IF ( .not. DO_SOLAR_SOURCES ) THEN

!  Only one solution, local direct_beam flag NOT set
!mick note 9/19/2017 - DO_INCLUDE_DIRECTBEAM in the thermal-only case here
!  covers the roles of both DO_LOCALBEAM(IBEAM) and DO_INCLUDE_DIRECTBEAM
!  in the solar case later

        IPARTIC = 1
!        DO_INCLUDE_DIRECTBEAM = .FALSE.       ! replaced 4/28/19
        DO_INCLUDE_DIRECTRF = .FALSE.
        DO_INCLUDE_DIRECTSL = .FALSE.

!  Version 2.8, Get rid of the GOTO 566 statement
!  Avoid the scattering solutions if not flagged
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 566

        IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Thermal-only. Find the BVP solution and intensity field
!  -------------------------------------------------------

!  set the BVP PI solution at the lower/upper boundaries
!          O1 = 1
!          DO LAYER = 1, NLAYERS
!            DO I = 1, NSTREAMS_2
!              WUPPER(I,O1,LAYER) = T_WUPPER(I,LAYER)
!              WLOWER(I,O1,LAYER) = T_WLOWER(I,LAYER)
!            ENDDO
!          ENDDO
!mick fix 2/17/11 - include remaining stokes components as needed
!Rob Fix  April 2016. Original code was correct
!          IF (NSTOKES .EQ. 1) THEN
!            DO LAYER = 1, NLAYERS
!              DO I = 1, NSTREAMS_2
!                WUPPER(I,1,LAYER) = T_WUPPER(I,LAYER) ; WLOWER(I,1,LAYER) = T_WLOWER(I,LAYER)
!              ENDDO
!            ENDDO
!          ELSE
!            DO LAYER = 1, NLAYERS
!              DO O1 = 1, NSTOKES
!                DO I = 1, NSTREAMS_2
!                  WUPPER(I,O1,LAYER) = T_WUPPER(I,LAYER) ; WLOWER(I,O1,LAYER) = T_WLOWER(I,LAYER)
!                ENDDO
!              ENDDO
!            ENDDO
!          END IF

!mick fix 9/19/2017 - still bug with original code;
!                     new fix: initialize all elements of WUPPER and WLOWER

          WUPPER = ZERO
          WLOWER = ZERO
          O1 = 1
          WUPPER(1:NSTREAMS_2,O1,1:NLAYERS) = T_WUPPER(1:NSTREAMS_2,1:NLAYERS)
          WLOWER(1:NSTREAMS_2,O1,1:NLAYERS) = T_WLOWER(1:NSTREAMS_2,1:NLAYERS)

!  Solve the boundary value problem. No telescoping here
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
      
!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19
!      -- Rob fix  4/9/2019  - Major overhaul for adjusted BVP solution

          CALL BVP_SOLUTION_MASTER ( & 
             DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_SURFEMISS,          & ! Input Surface Flags
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_INCLUDE_DIRECTRF,                 & ! Input Surface Flags
             DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX, FOURIER, IPARTIC,                 & ! Input illumination flags, indices
             DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,                    & ! Input Water-leaving flags
             NSTOKES, NSTREAMS, NLAYERS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,     & ! Numbers
             NTOTAL, N_SUBDIAG, N_SUPDIAG, K_REAL, K_COMPLEX, MUELLER_INDEX,           & ! Numbers, bookkeeping
             TF_MAXITER, TF_CRITERION, FLUX_FACTOR, FLUXVEC, SURFACE_FACTOR,           & ! Input factors/WL control
             QUAD_STRMWTS, ALBEDO, AXBID_F, SURFBB, EMISSIVITY,                        & ! Input surface terms
             SLTERM_ISOTROPIC, SLTERM_F_0, TOAFLUX, BOAFLUX,                           & ! Input Sleave/illumination
             SOLARBEAM_BOATRANS, LOCAL_CSZA, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, & ! Input Direct-flux
             RF_DIRECT_BEAM, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, WUPPER, WLOWER,       & ! Input RTE solutions
             BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                                         & ! Input BVP matrices/pivots             
             COL2, SCOL2, TRANS_ATMOS_FINAL, SL_QUADTERM, R2_BEAM, LCON, MCON,         & ! Modified input/output
             STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                                     ! Exception handling

!  Exception handling

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER
            TRACE_3 = 'Error return from BVP_SOLUTION_MASTER (thermal-only), ' // &
                      'Called in VLIDORT_FOURIER, Fourier # ' // CF
            STATUS = VLIDORT_SERIOUS ; RETURN
          ENDIF

!  Continuation point for avoiding thermal scattering
! 566  CONTINUE. removed, Versio 2.8

!  End thermal scattering clause

        ENDIF

!  Post-processing - upwelling thermal-only field
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!    - Version 2.8.1, Control for BOA illumination added, 3/23/19
!    - Version 2.8.1, 4/9/19. DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL added
!    - Version 2.8.1, 4/9/19. RF_USER_DIRECT_BEAM, SL_USERTERM have been added
        
        IF ( DO_UPWELLING ) THEN
          CALL VLIDORT_UPUSER_INTENSITY ( &
            DO_USER_STREAMS,    DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,              & ! Input flags (RT mode)
            DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,             & ! Input flags (sources)
            DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,   DO_INCLUDE_SURFEMISS,             & ! Input flags (Surface)
            DO_DBCORRECTION, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL, DO_LAYER_SCATTERING,& ! Input flags (Beam/scattering)
            DO_MSMODE_VLIDORT,  DO_MSMODE_THERMAL, DO_TOA_CONTRIBS, DO_INCLUDE_BOAFLUX,    & ! Input flags (RT mode)
            FOURIER, IPARTIC, NSTOKES, NSTREAMS, NLAYERS, N_USER_STREAMS,                  & ! Input numbers (basic)
            LOCAL_UM_START, N_USER_LEVELS, UTAU_LEVEL_MASK_UP, MUELLER_INDEX,              & ! Input bookkeeping + levels
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                  & ! Input partial-layer control
            FLUX_MULTIPLIER,  BOAFLUX, SURFACE_FACTOR, QUAD_WEIGHTS, QUAD_STRMWTS,         & ! Input Flux and quadrature
            T_DELT_DISORDS,  T_DELT_USERM, T_UTUP_USERM, CUMTRANS,                         & ! Input Transmittances
            ALBEDO, BRDF_F, USER_BRDF_F, SURFBB, EMISSIVITY, USER_EMISSIVITY,              & ! Input Surface BRDF/Emiss.
            K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,                         & ! Input Homog. RTE Soln.
            WLOWER, LCON, MCON, T_WLOWER, LAYER_TSUP_UP, LAYER_TSUP_UTUP,                  & ! Input RTE PI and thermal
            RF_USER_DIRECT_BEAM, SL_USERTERM, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2,  & ! Input User solutions
            HMULT_1, HMULT_2, EMULT_UP,  UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,            & ! Input multipliers
            STOKES_DOWNSURF, CUMSOURCE_UP, STOKES_F, MS_CONTRIBS_F, BOA_THTONLY_SOURCE )     ! OUTPUT
        ENDIF

!  Post-processing - Downwelling thermal-only field
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!    - Version 2.8.1, Control for TOA illumination added, 3/23/19
        
        IF ( DO_DNWELLING ) THEN
          CALL VLIDORT_DNUSER_INTENSITY ( &
            DO_USER_STREAMS,     DO_OBSERVATION_GEOMETRY, DO_INCLUDE_TOAFLUX,    & ! Input flags (RT mode)
            DO_SOLAR_SOURCES,    DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,  & ! Input flags (sources)
            DO_LAYER_SCATTERING, DO_MSMODE_VLIDORT,                              & ! Input flags (RT mode)
            FOURIER, IPARTIC, NSTOKES, N_USER_STREAMS,                           & ! Input numbers (basic)
            LOCAL_UM_START, N_USER_LEVELS, UTAU_LEVEL_MASK_DN,                   & ! Input bookkeeping + levels
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,        & ! Input partial-layer control
            FLUX_MULTIPLIER, TOAFLUX, T_DELT_USERM, T_UTDN_USERM,                & ! Input Transmittances, Flux
            K_REAL, K_COMPLEX, LCON, MCON, LAYER_TSUP_DN, LAYER_TSUP_UTDN,       & ! Input RTE Sol + thermal
            UHOM_DNDN, UHOM_DNUP, UPAR_DN_1, UPAR_DN_2,                          & ! Input User solutions
            HMULT_1, HMULT_2, EMULT_DN,  UT_HMULT_DD, UT_HMULT_DU, UT_EMULT_DN,  & ! Input multipliers
            CUMSOURCE_DN, STOKES_F )                                               ! OUTPUT
        ENDIF

!  mean value (integrated) output for thermal-only field
!    - Can also use this to get debug Quadrature output
!    - @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16. Removed QUAD_OUTPUT flag.
!mick mod 9/19/2017 - added LEVELS_SOLARTRANS & PARTIALS_SOLARTRANS to input arguments

!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19
        
        IF ( DO_INCLUDE_MVOUTPUT ) THEN
          CALL VLIDORT_INTEGRATED_OUTPUT ( &
            DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTRF, DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX, & ! Input flags
            DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,     & ! Input flags
            IPARTIC, NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, N_DIRECTIONS,  & ! Input numbers
            UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, WHICH_DIRECTIONS,          & ! Input output levels
            PARTLAYERS_LAYERIDX, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,      & ! Input partial layer control
            FLUX_MULTIPLIER, TOAFLUX, BOAFLUX, FLUX_FACTOR,                    & ! Flux inputs
            FLUXVEC, BEAM_CUTOFF, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,    & ! Quadrature
            LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,                            & ! Solar beam transmittances
            T_DELT_MUBAR, T_UTDN_MUBAR, INITIAL_TRANS, LOCAL_CSZA,             & ! Solar beam parameterization
            T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                    & ! Discrete ordinate transmittances
            K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, LCON, MCON,               & ! Input RTE
            T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, WLOWER, WUPPER,          & ! Input RTE
            T_WLOWER, T_WUPPER, UT_T_PARTIC, BOA_THTONLY_SOURCE,               & ! Input Thermal
            MEANST_DIFFUSE, FLUX_DIFFUSE, DNMEANST_DIRECT, DNFLUX_DIRECT )       ! Output
        ENDIF

!  Finish Thermal only.

        RETURN

!  Continuation point. removed for Version 2.8
! 455  CONTINUE

!  End thermal only clause (replaces GOTO)

      ENDIF

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

!   (classical flag removed, 7/8/16 for Version 2.8)
!      *** "IF ( DO_CLASSICAL_SOLUTION ) THEN" [after "IF ( DO_MULTIBEAM(IBEAM,FOURIER) ) THEN"]
!      *** Solar beam Particular solutions (Green's function) Clause removed......

!  Start loop over desired beam solutions

      DO IBEAM = 1, N_SZANGLES

        IF ( DO_MULTIBEAM(IBEAM,FOURIER) ) THEN

!  Solar beam Particular solutions 
!  -------------------------------

          DO LAYER = 1, NLAYERS

!  get the solution
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.

            CALL VLIDORT_QBEAM_SOLUTION ( &
              LAYER, FOURIER, IBEAM,                                  & ! Input indices
              NSTOKES, NSTREAMS, NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS,  & ! Input numbers
              FLUX_FACTOR, DFLUX, QUAD_STREAMS,                       & ! Input Flux and quadrature
              DO_LAYER_SCATTERING, OMEGA_GREEK, BEAM_CUTOFF,          & ! Input bookkeeping + optical
              T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT,            & ! Input beam attenuation
              PI_XQP, PI_XQM, PI_X0P, PI_XQM_POST,                    & ! Input PI matrices
              SAB, DAB, EIGENMAT_SAVE,                                & ! Input matrices from RTE
              QSUMVEC_SAVE, QDIFVEC_SAVE, QVEC_SAVE, QDIF_SAVE,       & ! Output beam auxiliary vectors
              QMAT_SAVE, QPIVOT, BVEC, WUPPER, WLOWER,                & ! Output matrix and Beam solutions
              STATUS_SUB, MESSAGE, TRACE_1 )                            ! Exception handling

!  error handling

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
              TRACE_2 = 'Error return from VLIDORT_QBEAM_SOLUTION, '//'Called in VLIDORT_FOURIER, Fourier # '//CF
              STATUS = VLIDORT_SERIOUS
              RETURN
            ENDIF

!  user solutions
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.

            IF ( STERM_LAYERMASK_UP(LAYER) .OR. STERM_LAYERMASK_DN(LAYER) ) THEN
              IF ( DO_USER_STREAMS ) THEN
                CALL VLIDORT_UBEAM_SOLUTION ( &
                  DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,                & ! Input flags
                  LAYER, FOURIER, IBEAM, NSTOKES, NSTREAMS, N_USER_STREAMS,           & ! Input numbers
                  NMOMENTS, LOCAL_UM_START, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,   & ! Input bookkeeping
                  DO_LAYER_SCATTERING, FLUX_FACTOR, DFLUX, BEAM_CUTOFF, QUAD_HALFWTS, & ! Input flux/quadrature/Bookkeeping
                  OMEGA_GREEK, PI_XQP, PI_XUP, PI_XUM, PI_X0P, PI_XQM_PRE, BVEC,      & ! Input optical/PI/Beam-PI
                  HELPSTOKES_BEAM, UPAR_DN_1, UPAR_DN_2, UPAR_UP_1, UPAR_UP_2 )         ! Output user solutions
              END IF
            END IF

!  end layer loop

          END DO

!  Add thermal solutions if flagged
!  ---------------------------------

!    NO modulus on the thermal contribution (Bug fixed 26 January 2010)

          IF ( DO_INCLUDE_THERMEMISS ) THEN
            O1 = 1
            DO N = 1, NLAYERS
             DO I = 1, NSTREAMS_2
! 1/26/10       WUPPER(I,O1,N) = WUPPER(I,O1,N) + pi4*T_WUPPER(I,N)
! 1/26/10       WLOWER(I,O1,N) = WLOWER(I,O1,N) + pi4*T_WLOWER(I,N)
                WUPPER(I,O1,N) = WUPPER(I,O1,N) + T_WUPPER(I,N)
                WLOWER(I,O1,N) = WLOWER(I,O1,N) + T_WLOWER(I,N)
              ENDDO
            ENDDO
          ENDIF

!  Solve boundary value problem
!  ----------------------------

!  standard case using compressed-band matrices, etc..
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.

          IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN

!  Get the BVP solution (regular case) for this beam component
!mick mod 9/19/2017 - replaced DO_REFLECTED_DIRECTBEAM(IBEAM) with DO_LOCALBEAM(IBEAM)

!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19
!      -- Rob fix  4/9/2019  - Major overhaul for adjusted BVP solution

            CALL BVP_SOLUTION_MASTER ( & 
               DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_SURFEMISS,          & ! Input Surface Flags
               DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_LOCALBEAM(IBEAM),                 & ! Input Surface Flags
               DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX, FOURIER, IBEAM,                   & ! Input illumination flags, indices
               DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,                    & ! Input Water-leaving flags
               NSTOKES, NSTREAMS, NLAYERS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,     & ! Numbers
               NTOTAL, N_SUBDIAG, N_SUPDIAG, K_REAL, K_COMPLEX, MUELLER_INDEX,           & ! Numbers, bookkeeping
               TF_MAXITER, TF_CRITERION, FLUX_FACTOR, FLUXVEC, SURFACE_FACTOR,           & ! Input factors/WL control
               QUAD_STRMWTS, ALBEDO, AXBID_F, SURFBB, EMISSIVITY,                        & ! Input surface terms
               SLTERM_ISOTROPIC, SLTERM_F_0, TOAFLUX, BOAFLUX,                           & ! Input Sleave/illumination
               SOLARBEAM_BOATRANS, LOCAL_CSZA, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, & ! Input Direct-flux
               RF_DIRECT_BEAM, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, WUPPER, WLOWER,       & ! Input RTE solutions
               BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                                         & ! Input BVP matrices/pivots
               COL2, SCOL2, TRANS_ATMOS_FINAL, SL_QUADTERM, R2_BEAM, LCON, MCON,         & ! Modified input/output
               STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                                     ! Exception handling

!  Exception handling

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
              TRACE_3 = 'Error return from BVP_SOLUTION_MASTER (Beam solution), '//&
                         'Called in VLIDORT_FOURIER, Fourier # '//CF
              STATUS = VLIDORT_SERIOUS ; RETURN
            ENDIF

!  Telescoped case
!mick mod 9/19/2017 - replaced DO_REFLECTED_DIRECTBEAM(IBEAM) with DO_LOCALBEAM(IBEAM)
!    - Version 2.8, Substantial upgrade to this routine (BRDF with Telescoping)
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.

          ELSE

            CALL BVPTEL_SOLUTION_MASTER ( &
              DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,                   & ! Flags
              DO_LOCALBEAM(IBEAM), FOURIER, IBEAM, NSTOKES,                & ! Numbers
              NSTREAMS, NLAYERS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, & ! Numbers
              N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,              & ! BVP Tel Control
              N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,                & ! BVP Tel Control
              MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR,            & ! Input Bookkeeping
              ALBEDO, AXBID_F, CUMTRANSDOM, CUMQUADDOM,                    & ! Input Surface inputs
              WLOWER, WUPPER, RF_DIRECT_BEAM,                              & ! Input RTE stuff
              SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, T_DELT_DISORDS,          & ! Input RTE stuff
              BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT, COLTEL2, SCOL2,      & ! Input BVProblem
              R2_BEAM, LCON, MCON, STATUS_SUB, MESSAGE, TRACE_1 )            ! Output and Status

!  Exception handling

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
              TRACE_2 = 'Error return from BVPTEL_SOLUTION_MASTER, '//'Called in VLIDORT_FOURIER, Fourier # '//CF
              STATUS = VLIDORT_SERIOUS
              RETURN
            ENDIF

          ENDIF

!  New 4/28/19. Diffuse and Direct solar flux at the bottom of the atmosphere
!  here is where you save the solutions you need for the planetary problem
!     SBTERM is already available, TRANSTERM(IB,UM) = TRANSBEAM(IB) * TRNMED_USER(UM) / PIE

          IF ( DO_PLANETARY_PROBLEM ) THEN
             N = NLAYERS ;  KO1 = K_REAL(N) + 1
             DO I = 1, NSTREAMS
                DO O1 = 1, NSTOKES
                   SHOM_R = ZERO ; SHOM_CR = ZERO
                   DO K = 1, K_REAL(N)
                      LXR  = LCON(K,N)*SOLA_XPOS(I,O1,K,N) * T_DELT_EIGEN(K,N)
                      MXR  = MCON(K,N)*SOLB_XNEG(I,O1,K,N)
                      SHOM_R = SHOM_R + LXR + MXR
                   ENDDO
                   DO K = 1, K_COMPLEX(N)
                      K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                      LXR_CR =  LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                      LXR_CI =  LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                      MXR_CR =  MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                      HOM1CR = LXR_CR*T_DELT_EIGEN(K1,N) - LXR_CI*T_DELT_EIGEN(K2,N)
                      SHOM_CR = SHOM_CR + HOM1CR + MXR_CR
                   ENDDO
                   TRANSQUAD(I,O1) = SHOM_R + SHOM_CR + WLOWER(I,O1,N)
                ENDDO
             ENDDO
             TRANSDIRECT = LOCAL_CSZA(NLAYERS,IBEAM) * FLUX_FACTOR * TRANS_SOLAR_BEAM(IBEAM)      ! Direct
             DO O1 = 1, NSTOKES
                TRANSBEAM(O1,IBEAM) = PI2 * DOT_PRODUCT(TRANSQUAD(1:NSTREAMS,O1),QUAD_STRMWTS(1:NSTREAMS)) ! Diffuse
                TRANSBEAM(O1,IBEAM) = TRANSBEAM(O1,IBEAM) + TRANSDIRECT * FLUXVEC(O1)
             ENDDO
          ENDIF            
           
!  Stokes Vector Post Processing
!  -----------------------------

!  upwelling
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.

          IF ( DO_UPWELLING ) THEN

!  VLIDORT COMMENTS:  Direct beam inclusion flag:
!   This now has the DBCORRECTION option: if the DBCORRECTION Flag
!   is set, then we will be doing exact calculations of the reflected
!   directbeam, so we do not need to include it in the Post-processing.
!   However, the direct beam will need to be included in the basic RT
!   solution (the BVP), and this is controlled separately by the
!   DO_REFLECTED_DIRECTBEAM(IBEAM) flags.
!     R. Spurr, RT Solutions, Inc., 19 August 2005.
             
!   Version 2.8.1,  3/23/2019. Introduce Control for including BOA illumination

!Rob fix 4/12/12 - added DO_MSMODE_VLIDORT
!            DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND. &
!                 (DO_REFLECTED_DIRECTBEAM(IBEAM) .AND. .NOT.DO_DBCORRECTION) ) .AND. .NOT.DO_MSMODE_VLIDORT

!Rob Fix 5/27/19. DO_INCLUDE_DIRECTRF same as before (DO_INCLUDE_DIRECTBEAM)
!                 DO_INCLUDE_DIRECTSL similarly constructed, using DO_DBCORRECTION = DO_FOCORR         

            DO_INCLUDE_DIRECTRF = &
             ( DO_UPWELLING .AND. ( DO_LOCALBEAM(IBEAM).AND..NOT.DO_DBCORRECTION) ) .and. .not. DO_MSMODE_VLIDORT
            DO_INCLUDE_DIRECTSL = &
             ( FOURIER.eq.0 .and. ( DO_SURFACE_LEAVING .AND..NOT.DO_DBCORRECTION) ) .and. .not. DO_MSMODE_VLIDORT

!  4/28/19. Version 2.8.1. Need to set the User-defined surface leaving, if not set.
!     -- Get adjusted User-term surface-leaving contribution
!  4/15/20. Version 2.8.2. SLEAVE arrays defined locally, remove Fourier index

            IF (  DO_INCLUDE_DIRECTSL .and. (DO_WATER_LEAVING .and..not.DO_EXTERNAL_WLEAVE) ) then
               O1 = 1 ; TFACTOR = TRANS_ATMOS_FINAL(IBEAM) * FLUX_FACTOR / DELTA_FACTOR
               IF ( DO_SL_ISOTROPIC ) THEN
                  SL = SLTERM_ISOTROPIC(O1,IBEAM) * TFACTOR
                  IF ( DO_USER_STREAMS ) THEN
                     IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
                        SL_USERTERM(1:N_USER_STREAMS,IBEAM,O1) =  SL
                     ELSE
                        SL_USERTERM(LUI,IBEAM,O1) = SL
                     ENDIF
                  ENDIF
               ELSE
                  IF ( DO_USER_STREAMS ) THEN
                     IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
                        SL_USERTERM(1:N_USER_STREAMS,IBEAM,O1) = USER_SLTERM_F_0(O1,1:N_USER_STREAMS,IBEAM) * TFACTOR
                     ELSE
                        SL_USERTERM(LUI,IBEAM,O1) = USER_SLTERM_F_0(O1,LUI,IBEAM) * TFACTOR
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF

!            write(*,*)'SLUSERTERM',IBEAM,DO_INCLUDE_DIRECTSL,FOURIER,SL_USERTERM(1,IBEAM,1)
            
!  Now, call the post-processing routine
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!    - Version 2.8.1, Control for BOA illumination added, 3/23/19
!    - Version 2.8.1, 4/9/19. DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL added
!    - Version 2.8.1, 4/9/19. RF_USER_DIRECT_BEAM, SL_USERTERM have been added
          
            CALL VLIDORT_UPUSER_INTENSITY ( &
              DO_USER_STREAMS,    DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,            & ! Input flags (RT mode)
              DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,           & ! Input flags (sources)
              DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,   DO_INCLUDE_SURFEMISS,           & ! Input flags (Surface)
              DO_DBCORRECTION, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL, DO_LAYER_SCATTERING,& ! Input flags (Beam/scattering)
              DO_MSMODE_VLIDORT, DO_MSMODE_THERMAL, DO_TOA_CONTRIBS, DO_INCLUDE_BOAFLUX,   & ! Input flags (RT mode)
              FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_USER_STREAMS,                  & ! Input numbers (basic)
              LOCAL_UM_START, N_USER_LEVELS, UTAU_LEVEL_MASK_UP, MUELLER_INDEX,            & ! Input bookkeeping + levels
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                & ! Input partial-layer control
              FLUX_MULTIPLIER, BOAFLUX, SURFACE_FACTOR, QUAD_WEIGHTS, QUAD_STRMWTS,        & ! Input Flux and quadrature
              T_DELT_DISORDS,  T_DELT_USERM, T_UTUP_USERM, CUMTRANS,                       & ! Input Transmittances
              ALBEDO, BRDF_F, USER_BRDF_F, SURFBB, EMISSIVITY, USER_EMISSIVITY,            & ! Input Surface BRDF/Emiss.
              K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,                       & ! Input Homog. RTE Soln.
              WLOWER, LCON, MCON, T_WLOWER, LAYER_TSUP_UP, LAYER_TSUP_UTUP,                & ! Input RTE PI and thermal
              RF_USER_DIRECT_BEAM, SL_USERTERM, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2,& ! Input User solutions
              HMULT_1, HMULT_2, EMULT_UP,  UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,          & ! Input multipliers
              STOKES_DOWNSURF, CUMSOURCE_UP, STOKES_F, MS_CONTRIBS_F, BOA_THTONLY_SOURCE )   ! OUTPUT
          ENDIF

!  Downwelling
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!    - Version 2.8.1,  3/23/19. Introduce Control for including TOA illumination
        
          IF ( DO_DNWELLING ) THEN
            CALL VLIDORT_DNUSER_INTENSITY ( &
              DO_USER_STREAMS,     DO_OBSERVATION_GEOMETRY, DO_INCLUDE_TOAFLUX,    & ! Input flags (RT mode)
              DO_SOLAR_SOURCES,    DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,  & ! Input flags (sources)
              DO_LAYER_SCATTERING, DO_MSMODE_VLIDORT,                              & ! Input flags (RT mode)
              FOURIER, IBEAM, NSTOKES, N_USER_STREAMS,                             & ! Input numbers (basic)
              LOCAL_UM_START, N_USER_LEVELS, UTAU_LEVEL_MASK_DN,                   & ! Input bookkeeping + levels
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,        & ! Input partial-layer control
              FLUX_MULTIPLIER, TOAFLUX, T_DELT_USERM, T_UTDN_USERM,                & ! Input Transmittances, Flux
              K_REAL, K_COMPLEX, LCON, MCON, LAYER_TSUP_DN, LAYER_TSUP_UTDN,       & ! Input RTE Sol + thermal
              UHOM_DNDN, UHOM_DNUP, UPAR_DN_1, UPAR_DN_2,                          & ! Input User solutions
              HMULT_1, HMULT_2, EMULT_DN,  UT_HMULT_DD, UT_HMULT_DU, UT_EMULT_DN,  & ! Input multipliers
              CUMSOURCE_DN, STOKES_F )                                               ! OUTPUT
          ENDIF

!  mean value (integrated) output
!    - Can also use this to get debug Quadrature output
!    - @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16. Removed QUAD_OUTPUT flag.
!mick mod 9/19/2017 - replaced DO_INCLUDE_DIRECTBEAM with DO_LOCALBEAM(IBEAM)
!                   - added LEVELS_SOLARTRANS & PARTIALS_SOLARTRANS to input arguments

!   Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination
!  --- Also include value of this illumination
!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL VLIDORT_INTEGRATED_OUTPUT ( &
              DO_INCLUDE_MVOUTPUT, DO_LOCALBEAM(IBEAM), DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX, & ! Input flags
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,    & ! Input flags
              IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, N_DIRECTIONS,   & ! Input numbers
              UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, WHICH_DIRECTIONS,         & ! Input output levels
              PARTLAYERS_LAYERIDX, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,     & ! Input partial layer control
              FLUX_MULTIPLIER, TOAFLUX, BOAFLUX, FLUX_FACTOR, FLUXVEC, BEAM_CUTOFF,   & ! Flux inputs
              QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                         & ! Quadrature
              LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,                           & ! Solar beam transmittances
              T_DELT_MUBAR, T_UTDN_MUBAR, INITIAL_TRANS, LOCAL_CSZA,            & ! Solar beam parameterization
              T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                   & ! Discrete ordinate transmittances
              K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, LCON, MCON,              & ! Input RTE
              T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, WLOWER, WUPPER,         & ! Input RTE
              T_WLOWER, T_WUPPER, UT_T_PARTIC, BOA_THTONLY_SOURCE,              & ! Input Thermal
              MEANST_DIFFUSE, FLUX_DIFFUSE, DNMEANST_DIRECT, DNFLUX_DIRECT )      ! Output
          ENDIF

!  End loop over beam solutions

        END IF
      END DO

!  Debug write
!      IF ( DO_DEBUG_WRITE ) THEN
!        write(*,'(a,I3,a,100I3)')'Fourier ',FOURIER, ' : # complex evalues by layer: ', (K_COMPLEX(LAYER),LAYER=1,NLAYERS)
!      ENDIF

!  ######
!  finish
!  ######

      RETURN
      END SUBROUTINE VLIDORT_FOURIER

!  End Module

      END MODULE VLIDORT_RTCalc_Master_m

