
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p1                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr                          #
! #                Matt Christi                                 #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p1                             #
! #  Release Date :   31 August 2019                            #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released August 2014                        #
! #      2.8   F90, released May    2017                        #
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
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.1 of the VLIDORT software library.          #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2019.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! #                                                                 #
! # This file is part of VLIDORT_2p8p1 ( Version 2.8.1 )            #
! #                                                                 #
! # VLIDORT_2p8p1 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p1 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p1  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            VLIDORT_MASTER  (master, public)                 #
! #            VLIDORT_FOURIER (called by master routine)       #
! #                                                             #
! ###############################################################


      MODULE vlidort_masters_m

      PUBLIC  :: VLIDORT_MASTER
      PRIVATE :: VLIDORT_FOURIER

      CONTAINS

      SUBROUTINE VLIDORT_MASTER ( &
        VLIDORT_FixIn, &
        VLIDORT_ModIn, &
        VLIDORT_Sup,   &
        VLIDORT_Out )

!  parameter file

      USE VLIDORT_PARS_m

!  I/O Type structures

      USE VLIDORT_Inputs_def_m
      USE VLIDORT_Sup_InOut_def_m
      USE VLIDORT_Outputs_def_m

!  internal Type structure

      USE VLIDORT_Work_def_m

!  Modules. Version 2.8, only use selections

      USE VLIDORT_INPUTS_m      , only : VLIDORT_CHECK_INPUT_DIMS, VLIDORT_CHECK_INPUT, VLIDORT_DERIVE_INPUT
      USE VLIDORT_GEOMETRY_m    , only : VLIDORT_CHAPMAN
      USE VLIDORT_MISCSETUPS_m  , only : VLIDORT_MISCSETUPS
      USE VLIDORT_THERMALSUP_m  , only : THERMAL_SETUP
      USE VLIDORT_MULTIPLIERS_m , only : EMULT_MASTER, EMULT_MASTER_OBSGEO
      USE VLIDORT_INTENSITY_m   , only : VLIDORT_CONVERGE, VLIDORT_CONVERGE_OBSGEO

      USE VLIDORT_PACK_m
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

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs), INTENT (IN)       :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (INOUT) :: VLIDORT_ModIn

!  VLIDORT supplements structure

      TYPE(VLIDORT_Sup_InOut), INTENT (INOUT)       :: VLIDORT_Sup

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs), INTENT (OUT)           :: VLIDORT_Out

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

      LOGICAL ::            DO_FULLRAD_MODE

      LOGICAL ::            DO_THERMAL_EMISSION
      LOGICAL ::            DO_SURFACE_EMISSION
      LOGICAL ::            DO_PLANE_PARALLEL

      LOGICAL ::            DO_UPWELLING
      LOGICAL ::            DO_DNWELLING

!      LOGICAL ::            DO_QUAD_OUTPUT       ! removed, 7/7/16
      LOGICAL ::            DO_TOA_CONTRIBS

      LOGICAL ::            DO_LAMBERTIAN_SURFACE

      LOGICAL ::            DO_SPECIALIST_OPTION_1
      LOGICAL ::            DO_SPECIALIST_OPTION_2
      LOGICAL ::            DO_SPECIALIST_OPTION_3

!  New 17 May 2012, surface leaving flags
      LOGICAL ::            DO_SURFACE_LEAVING
      LOGICAL ::            DO_SL_ISOTROPIC
      LOGICAL ::            DO_WATER_LEAVING       ! New 10/28/15
      LOGICAL ::            DO_FLUORESCENCE        ! New 10/28/15
      LOGICAL ::            DO_TF_ITERATION        ! New 7/7/16

!  Water-leaving output flag. 3/18/19 for Version 2.8.1
      
      LOGICAL ::            DO_WLADJUSTED_OUTPUT 

!   4/26/19 Added control for the media problem. Version 2.8.1
!     Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!       1 = Isotropic illumination from TOA, 2 = Isotropic illumination from BOA
      
      LOGICAL ::            DO_ALBTRN_MEDIA(2)

!   4/28/19 Added control for the planetary problem. Version 3.8.1
      !     Linked to the Media-problem, requires some flux and tranmsittance output (BOA unit illumination)

      LOGICAL ::            DO_PLANETARY_PROBLEM

!  TOA and BOA Illumination flags. 3/23/19 for Version 2.8.1
!   Airglow and Nighttime-viewing scenarios
      
      LOGICAL ::            DO_TOAFLUX
      LOGICAL ::            DO_BOAFLUX

!  Order of Taylor series (including terms up to EPS^n).
!        Introduced 2/19/14 for Version 2.7
      
      INTEGER ::            TAYLOR_ORDER

!  Basic control numbers

      INTEGER ::            NSTOKES
      INTEGER ::            NSTREAMS
      INTEGER ::            NLAYERS
      INTEGER ::            NFINELAYERS
      INTEGER ::            N_THERMAL_COEFFS

!  Accuracy

      DOUBLE PRECISION ::   VLIDORT_ACCURACY

!  Special inputs (RT Solutions use only)

      INTEGER ::            NLAYERS_NOMS
      INTEGER ::            NLAYERS_CUTOFF

!  Water-leaving: Control for iterative calculation of transmittanaces
!    Variables added for Version 2.8, 7/6/16

      INTEGER          ::   TF_MAXITER
      DOUBLE PRECISION ::   TF_CRITERION

!  TOA Illumination, Flux value. 3/23/19 for Version 2.8.1
!    Must be solar-flux normalized. Designed for Airglow Studies
      
      DOUBLE PRECISION ::   TOAFLUX
      
!  BOA Illumination, Flux value. 3/23/19 for Version 2.8.1
!    Must be solar-flux normalized. Designed for Nighttime Studies.
      
      DOUBLE PRECISION ::   BOAFLUX

!  Flux factor

      DOUBLE PRECISION ::   FLUX_FACTOR

!  User levels

      INTEGER ::            N_USER_LEVELS

!  PTH and fine grid

      DOUBLE PRECISION ::   HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION ::   PRESSURE_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION ::   TEMPERATURE_GRID ( 0:MAXLAYERS )
      INTEGER ::            FINEGRID ( MAXLAYERS )
      DOUBLE PRECISION ::   RFINDEX_PARAMETER

!  Optical properties

      DOUBLE PRECISION ::   DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION ::   GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  F-matrix optical properties (New for Version 2.8, 7/7/16)
      DOUBLE PRECISION ::   FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 )
      DOUBLE PRECISION ::   FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 )

      DOUBLE PRECISION ::   ALBEDO
      DOUBLE PRECISION ::   THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION ::   SURFBB

!  These "LTE" inputs were in Version 2.6.
!      DOUBLE PRECISION ::   LTE_DELTAU_VERT_INPUT ( 2, MAXLAYERS )
!      DOUBLE PRECISION ::   LTE_THERMAL_BB_INPUT ( 0:MAXLAYERS )

      DOUBLE PRECISION ::   ATMOS_WAVELENGTH

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

!  FOCorr and SSCorr Booleans. Completely revised for Version 2.8, 3/1/17.

      LOGICAL ::          DO_FOCORR                ! New 02 Jul 2013, used to be DO_FO_CALC
      LOGICAL ::          DO_FOCORR_EXTERNAL       ! renamed 2.8
      LOGICAL ::          DO_FOCORR_ALONE          ! Used to be DO_SSFULL - internal variable now (Version 2.8)
      LOGICAL ::          DO_FOCORR_NADIR
      LOGICAL ::          DO_FOCORR_OUTGOING

      LOGICAL ::          DO_SSCORR_TRUNCATION
      LOGICAL ::          DO_SSCORR_USEFMAT        ! New, 2.8

!  Additional Control for Externalized input (SLEAVE). Introduced 3/18/19 for Version 2.8.1
      
      LOGICAL ::          DO_EXTERNAL_WLEAVE

!  Solar

      LOGICAL ::          DO_SOLAR_SOURCES
      LOGICAL ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL ::          DO_CHAPMAN_FUNCTION

!  Performance flags

      LOGICAL ::          DO_DOUBLE_CONVTEST
      LOGICAL ::          DO_RAYLEIGH_ONLY
      LOGICAL ::          DO_DELTAM_SCALING
      LOGICAL ::          DO_SOLUTION_SAVING
      LOGICAL ::          DO_BVP_TELESCOPING

!  Other RT model flags

      LOGICAL ::          DO_USER_VZANGLES
      LOGICAL ::          DO_ADDITIONAL_MVOUT
      LOGICAL ::          DO_MVOUT_ONLY
      LOGICAL ::          DO_THERMAL_TRANSONLY
      LOGICAL ::          DO_OBSERVATION_GEOMETRY

!  Integers

      INTEGER ::          NGREEK_MOMENTS_INPUT

      INTEGER ::          N_SZANGLES
      INTEGER ::          N_USER_RELAZMS
      INTEGER ::          N_USER_VZANGLES
      INTEGER ::          N_USER_OBSGEOMS

!  Angles and levels

      DOUBLE PRECISION :: SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      DOUBLE PRECISION :: USER_VZANGLES ( MAX_USER_VZANGLES )
      DOUBLE PRECISION :: USER_LEVELS ( MAX_USER_LEVELS )
      DOUBLE PRECISION :: USER_OBSGEOMS ( MAX_USER_OBSGEOMS, 3 )

!  Spec height and earth radius

      DOUBLE PRECISION :: GEOMETRY_SPECHEIGHT
      DOUBLE PRECISION :: EARTH_RADIUS

!  Single scattering albedo

      DOUBLE PRECISION :: OMEGA_TOTAL_INPUT ( MAXLAYERS )

!  -----------------------
!  Standard Supplement I/O
!  -----------------------

!  BRDF Inputs
!  -----------

      DOUBLE PRECISION ::   EXACTDB_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION ::   BRDF_F_0 ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   BRDF_F   ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION ::   USER_BRDF_F_0 ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_BRDF_F   ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION ::   EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION ::   USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  SS and DB I/O
!  -------------

!  SS Intensity Results at all angles and optical depths

      DOUBLE PRECISION ::  STOKES_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  DB Intensity Results at all angles and optical depths

      DOUBLE PRECISION ::  STOKES_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Contribution function (TOA Upwelling only)
!  SS component of Diffuse Field
!      DOUBLE PRECISION ::  SS_CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  New Surface-Leaving Inputs, 17 May 12
!  -------------------------------------

      DOUBLE PRECISION ::   SLTERM_ISOTROPIC  ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION ::   SLTERM_USERANGLES ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

      DOUBLE PRECISION ::   SLTERM_F_0      ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_SLTERM_F_0 ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  ----------------
!  Standard Outputs
!  ----------------

!  Fourier-summed values

      DOUBLE PRECISION :: STOKES ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  Results for mean-value output
!  -----------------------------
!mick mod 9/19/2017 - renamed fluxes for separating diffuse and direct

!  Complete Actinic and Regular Fluxes (including Direct terms)

      DOUBLE PRECISION :: MEANST_DIFFUSE ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FLUX_DIFFUSE   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )

!  Direct Fluxes only

      DOUBLE PRECISION :: DNMEANST_DIRECT ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION :: DNFLUX_DIRECT   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

!  MEDIA-PROPERTY Output
!  ---------------------
      
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
!  -------------------------------------------

!  Fourier component of Diffuse Field
!      DOUBLE PRECISION ::  MS_CONTRIBS_F ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
!  Fourier-summed values
!      DOUBLE PRECISION ::  CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  Ancillary Output
!  ----------------

!  Fourier numbers used

      INTEGER ::          FOURIER_SAVED ( MAX_SZANGLES )

!  Number of geometries computed

      INTEGER ::          N_GEOMETRIES

!  Offsets for indexing geometries

      INTEGER ::          SZA_OFFSETS ( MAX_SZANGLES )
      INTEGER ::          VZA_OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )

!  Error handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NCHECKMESSAGES
      CHARACTER (LEN=120) :: CHECKMESSAGES (0:MAX_MESSAGES)
      CHARACTER (LEN=120) :: ACTIONS       (0:MAX_MESSAGES)
      INTEGER ::             STATUS_CALCULATION
      CHARACTER (LEN=120) :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

! ######################################################################

!                       Local Arguments
!                       +++++++++++++++

! ######################################################################

!  Local bookkeeping
!  ----------------

!               Intent(In) To the Fourier  routine
!               Intent(In) To the Converge routine

!  Mode flags (Robfix, 13 January 2012. Add MSMODE_THERMAL flag)

      LOGICAL         :: DO_MSMODE_VLIDORT
      LOGICAL         :: DO_MSMODE_THERMAL

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER         :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER         :: NSTREAMS_2
      INTEGER         :: NTOTAL
      INTEGER         :: N_SUBDIAG, N_SUPDIAG

!  Quadrature weights and abscissae, and product

      DOUBLE PRECISION :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_HALFWTS ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_SINES   ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_ANGLES  ( MAXSTREAMS )

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      !DOUBLE PRECISION :: USER_VZANGLES ( MAX_USER_STREAMS )
      DOUBLE PRECISION :: USER_SINES   ( MAX_USER_STREAMS )
      DOUBLE PRECISION :: USER_STREAMS ( MAX_USER_STREAMS )
      DOUBLE PRECISION :: USER_SECANTS ( MAX_USER_STREAMS )

!  Output optical depth masks and indices

      LOGICAL          :: PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER          :: PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER          :: UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER          :: UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!  Off-grid optical depths (values, masks, indices)

      INTEGER          :: N_PARTLAYERS
      INTEGER          :: PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      DOUBLE PRECISION :: PARTLAYERS_VALUES   ( MAX_PARTLAYERS )

!  Layer masks for doing integrated source terms

      LOGICAL         :: STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL         :: STERM_LAYERMASK_DN ( MAXLAYERS )

!  Indexing numbers

      INTEGER         :: N_VIEWING

!  Offsets for geometry indexing

      !INTEGER         :: IBOFF(MAXBEAMS)
      !INTEGER         :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Local input solar zenith angles Cosines
!  ( Only required for refractive geometry attenuation of the solar beam)

      DOUBLE PRECISION :: SUN_SZA_COSINES ( MAXLAYERS, MAX_SZANGLES )

!  Local solar zenith angles Cosines (regular case)

      DOUBLE PRECISION :: COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION :: SIN_SZANGLES ( MAX_SZANGLES )

!  Solar beam flags (always internal)

      LOGICAL         :: DO_MULTIBEAM ( MAXBEAMS, 0:MAXFOURIER )

!  Number of directions (1 or 2) and directional array

      INTEGER         :: N_DIRECTIONS
      INTEGER         :: WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Number of convergence tests

      INTEGER         :: N_CONVTESTS
      INTEGER         :: N_CONV_STREAMS

!  Arrays for setups and Corrections
!  ---------------------------------

!               Intent(In) To the Fourier routine

!  Local flags for the solution saving option

      INTEGER          :: LAYER_MAXMOMENTS ( MAXLAYERS )

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

      DOUBLE PRECISION :: STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

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

!  Local variables
!  ---------------

!  Local quantities for the water-leaving configuration

      DOUBLE PRECISION :: CUMSOURCE_DB, SLTERM_LOCAL, TRANS, TFACTOR

!  Local flags for media problem and Planetary problem control
      
      LOGICAL ::          LOCAL_DO_ALBTRN_MEDIA(2)
      LOGICAL ::          LOCAL_DO_PLANETARY_PROBLEM
      
!  Local flags

      LOGICAL ::          LOCAL_DO_NO_AZIMUTH
      LOGICAL ::          SAVE_DO_NO_AZIMUTH
      LOGICAL ::          LOCAL_ITERATION

!  inclusion flags

      LOGICAL ::          DO_INCLUDE_SURFACE
      LOGICAL ::          DO_INCLUDE_SURFEMISS
      LOGICAL ::          DO_INCLUDE_THERMEMISS

!  Fourier

      INTEGER ::          FOURIER
      INTEGER ::          N_FOURIERS

!  Misc. help

      INTEGER ::          OFF, O1, UA, UM, IB, G, G1, G2, L, UTA, LUM, LUA, NMOMS
      INTEGER ::          TESTCONV, LOCAL_N_USERAZM, STATUS_SUB
      INTEGER ::          IUNIT, SUNIT, RUNIT

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

!  Transmittances and Fluxes for Water-leaving, Master routine (Mark 3)
!  --------------------------------------------------------------------
!    - Introduced 1/12/16, Validated 2/3/16. DISABLED Version 2.8.1.
!      DOUBLE PRECISION :: TFACTOR
!      DOUBLE PRECISION :: TRANS_ATMOS_FINAL  ( MAX_SZANGLES )
!      DOUBLE PRECISION :: FLUX_DIFFUSE_FINAL  ( MAX_SZANGLES )

!  Helper variables
!  ----------------

!  flags

      LOGICAL ::          DO_NO_AZIMUTH
      LOGICAL ::          DO_ALL_FOURIER
      LOGICAL ::          DO_DIRECT_BEAM
!      LOGICAL ::          DO_CLASSICAL_SOLUTION  !!! removed for Version 2.8
      LOGICAL ::          DO_DBCORRECTION

!  Numbers

      INTEGER ::          NSTOKES_SQ
      INTEGER ::          NSTKS_NSTRMS
      INTEGER ::          NSTKS_NSTRMS_2
      INTEGER ::          NBEAMS
      INTEGER ::          NPARTICSOLS

!  Single scatter flux multipier

      DOUBLE PRECISION :: SS_FLUX_MULTIPLIER

!  Solution bookkeeping

      INTEGER ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: DMAT ( MAXSTOKES, MAXSTOKES )
      INTEGER ::          GREEKMAT_INDEX ( 6 )
      LOGICAL ::          DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION :: DFLUX ( MAXSTOKES )

!  output control

      INTEGER ::          LOCAL_UM_START
      INTEGER ::          N_OUT_STREAMS
      DOUBLE PRECISION :: OUT_ANGLES ( MAX_USER_STREAMS )

! #################### INNOVATIONS 5/5/20 ##############
!  4/29/20. SPECIALIST OUTPUT. Rayleigh Fourier output (TOA Upwelling only)
!      LOGICAL ::          SPECIAL_RAYF_OUTPUT
! #################### INNOVATIONS 5/5/20 ##############

!  level and partial control

      LOGICAL ::          DO_PARTLAYERS
      INTEGER ::          N_LAYERSOURCE_UP
      INTEGER ::          N_LAYERSOURCE_DN
      INTEGER ::          N_ALLLAYERS_UP
      INTEGER ::          N_ALLLAYERS_DN

!  Chapman factors
!mick fix 9/19/2017 - added PARTIAL_CHAPFACS

      DOUBLE PRECISION :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: PARTIAL_CHAPFACS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Local optical depths

      DOUBLE PRECISION :: TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION :: TAUGRID ( 0:MAXLAYERS )

!  Contribution functions (TOA Upwelling only)
!  -------------------------------------------

!  Fourier component of Diffuse Field

      DOUBLE PRECISION ::  MS_CONTRIBS_F ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )

!  Fourier-summed values

      DOUBLE PRECISION ::  CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  Single scatter

      DOUBLE PRECISION :: SS_CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

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

      LOGICAL ::          FAIL

!  Test variables
!      LOGICAL ::          DO_FDTEST =.FALSE.

!  This is the flag for a complete write-up of the inputs

      LOGICAL ::          DO_DEBUG_INPUT=.FALSE.
!      LOGICAL ::          DO_DEBUG_INPUT=.TRUE.

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

      HEIGHT_GRID(0:NLAYERS)      = VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(0:NLAYERS)
      PRESSURE_GRID(0:NLAYERS)    = VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID(0:NLAYERS)
      TEMPERATURE_GRID(0:NLAYERS) = VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(0:NLAYERS)
      FINEGRID(1:NLAYERS)         = VLIDORT_FixIn%Chapman%TS_FINEGRID(1:NLAYERS)
      RFINDEX_PARAMETER           = VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER

!  Fixed Optical inputs
!  --------------------

      DELTAU_VERT_INPUT(1:NLAYERS) = VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(1:NLAYERS)
      GREEKMAT_TOTAL_INPUT(0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,1:NLAYERS,:) = &
               VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,1:NLAYERS,:)

!  New version 2.8, F-matrix input. 7/7/16
!mick fix 9/19/2017 - swapped layer & geo indices

      !FMATRIX_UP(:,1:NLAYERS,:) = VLIDORT_FixIn%Optical%TS_FMATRIX_UP(:,1:NLAYERS,:)
      !FMATRIX_DN(:,1:NLAYERS,:) = VLIDORT_FixIn%Optical%TS_FMATRIX_DN(:,1:NLAYERS,:)
      FMATRIX_UP(1:NLAYERS,:,:) = VLIDORT_FixIn%Optical%TS_FMATRIX_UP(1:NLAYERS,:,:)
      FMATRIX_DN(1:NLAYERS,:,:) = VLIDORT_FixIn%Optical%TS_FMATRIX_DN(1:NLAYERS,:,:)

      ALBEDO     = VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO

      THERMAL_BB_INPUT(0:NLAYERS)  = VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT(0:NLAYERS)
      SURFBB                = VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT

!  @@@ Rob Fix 1/31/11, TS_EMISSIVITY, TS_USER_EMISSIVITY are defined in
!                       in Structure VLIDORT_Sup%BRDF, so do not need to be copied here.
!      EMISSIVITY            = VLIDORT_FixIn%Optical%TS_EMISSIVITY
!      USER_EMISSIVITY       = VLIDORT_FixIn%Optical%TS_USER_EMISSIVITY
!  Replaced by LBBF facility. Version 2.7
!      LTE_DELTAU_VERT_INPUT(1:2,1:NLAYERS) = &
!        VLIDORT_FixIn%Optical%TS_LTE_DELTAU_VERT_INPUT(1:2,1:NLAYERS)
!      LTE_THERMAL_BB_INPUT(0:NLAYERS)      = &
!        VLIDORT_FixIn%Optical%TS_LTE_THERMAL_BB_INPUT(0:NLAYERS)

      ATMOS_WAVELENGTH        = VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH

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

      DO_SSCORR_TRUNCATION    = VLIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION
      DO_SSCORR_USEFMAT       = VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT

!  Other Booleans

      DO_DOUBLE_CONVTEST      = VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST

      DO_SOLAR_SOURCES        = VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES
      DO_REFRACTIVE_GEOMETRY  = VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY
      DO_CHAPMAN_FUNCTION     = VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION

      DO_RAYLEIGH_ONLY        = VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY
      DO_DELTAM_SCALING       = VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING

      DO_SOLUTION_SAVING      = VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING
      DO_BVP_TELESCOPING      = VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING

      DO_USER_VZANGLES        = VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES
      DO_ADDITIONAL_MVOUT     = VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT
      DO_MVOUT_ONLY           = VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY
      DO_THERMAL_TRANSONLY    = VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY
      DO_OBSERVATION_GEOMETRY = VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY

!  Additional Control for Externalized Water-leaving input. Introduced 3/18/19 for Version 2.8.1

      DO_EXTERNAL_WLEAVE     = VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE

!  Modified Control inputs
!  -----------------------

      NGREEK_MOMENTS_INPUT = VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT

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

      GEOMETRY_SPECHEIGHT              = VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT

      N_USER_OBSGEOMS                      = VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
      USER_OBSGEOMS(1:N_USER_OBSGEOMS,1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,1:3)

!  Modified Chapman Function inputs
!  --------------------------------

      !CHAPMAN_FACTORS     = VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS
      EARTH_RADIUS        = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS

!  Modified Optical inputs
!  -----------------------

      OMEGA_TOTAL_INPUT(1:NLAYERS) = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS)

!  BRDF inputs
!  -----------
!mick mod 9/19/2017 - added IF conditions
!                   - note: emissivities left out of BRDF SURFACE block due to possibility
!                           of thermal being used in the Lambertian case

      IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
        BRDF_F_0(0:NMOMS,:,1:NSTREAMS,1:N_SZANGLES) = &
          VLIDORT_Sup%BRDF%TS_BRDF_F_0(0:NMOMS,:,1:NSTREAMS,1:N_SZANGLES)
        BRDF_F(0:NMOMS,:,1:NSTREAMS,1:NSTREAMS)     = &
          VLIDORT_Sup%BRDF%TS_BRDF_F  (0:NMOMS,:,1:NSTREAMS,1:NSTREAMS)
        IF ( DO_USER_VZANGLES ) THEN
          EXACTDB_BRDFUNC (:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
            VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC(:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
          USER_BRDF_F_0(0:NMOMS,:,1:N_USER_VZANGLES,1:N_SZANGLES) = &
            VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0(0:NMOMS,:,1:N_USER_VZANGLES,1:N_SZANGLES)
          USER_BRDF_F(0:NMOMS,:,1:N_USER_VZANGLES,1:NSTREAMS)      = &
            VLIDORT_Sup%BRDF%TS_USER_BRDF_F(0:NMOMS,:,1:N_USER_VZANGLES,1:NSTREAMS)
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

      IF ( DO_SURFACE_LEAVING ) THEN
        SLTERM_ISOTROPIC(1:NSTOKES,1:N_SZANGLES) = &
          VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(1:NSTOKES,1:N_SZANGLES)
        SLTERM_F_0(0:NMOMS,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES) = &
          VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0(0:NMOMS,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES)
        IF ( DO_USER_VZANGLES ) THEN
          SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
            VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
          USER_SLTERM_F_0(0:NMOMS,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES) = &
            VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(0:NMOMS,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES)
        ENDIF
      ENDIF

!  SS inputs
!  ---------
!  New 12 March 2012 --> IF SS results already available copy them !
!     Modified flagging, Version 2.8

      IF ( DO_FOCORR_EXTERNAL ) THEN
        STOKES_SS(1:N_USER_LEVELS,:,1:NSTOKES,:) = &
          VLIDORT_Sup%SS%TS_STOKES_SS(1:N_USER_LEVELS,:,1:NSTOKES,:)
        STOKES_DB(1:N_USER_LEVELS,:,1:NSTOKES) = &
          VLIDORT_Sup%SS%TS_STOKES_DB(1:N_USER_LEVELS,:,1:NSTOKES)
      ENDIF
     
! #################### INNOVATIONS 5/5/20 ##############
!  4/29/20. Set Flag: Special Fourier-component output for Rayleigh + Planetary-problem TOA Upwelling situations
!      SPECIAL_RAYF_OUTPUT = &
!       ( DO_PLANETARY_PROBLEM .and.DO_UPWELLING .and..not.DO_FOCORR .and.DO_RAYLEIGH_ONLY .and. (USER_LEVELS(1).eq.zero) )
! #################### INNOVATIONS 5/5/20 ##############

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  VLIDORT input debug

      IF (DO_DEBUG_INPUT) THEN
        CALL VLIDORT_DEBUG_INPUT_MASTER()
      END IF

!  Initialize outputs
!  ------------------

!  Main outputs (Radiances and fluxes)

      STOKES(1:N_USER_LEVELS,:,1:NSTOKES,:) = ZERO
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

      VLIDORT_Out%Main%TS_PLANETARY_SBTERM    = ZERO
      VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM = ZERO

!  TF outputs, Rob fix, RT Solutions, 9/25/15, Superceded 2/3/16
!      MEANST_DIFFUSE_TF(1:N_SZANGLES,1:NSTOKES) = ZERO
!      FLUX_DIFFUSE_TF(1:N_SZANGLES,1:NSTOKES) = ZERO
!      DNMEANST_DIRECT_TF(1:N_SZANGLES,1:NSTOKES) = ZERO
!      DNFLUX_DIRECT_TF(1:N_SZANGLES,1:NSTOKES) = ZERO

!  SS inputs
!  ---------

!  New 12 March 2012 --> IF SS results already available copy them !
!     Modified flagging, Version 2.8

      IF ( .NOT. DO_FOCORR_EXTERNAL ) THEN
         STOKES_SS(1:N_USER_LEVELS,:,1:NSTOKES,:) = ZERO
         STOKES_DB(1:N_USER_LEVELS,:,1:NSTOKES)   = ZERO
      ENDIF

!  Bookkeeping
!  -----------
      
!  Zeroing offsets 4/9/19

      FOURIER_SAVED(1:N_SZANGLES) = 0
      N_GEOMETRIES                = 0
      SZA_OFFSETS(1:N_SZANGLES)   = 0
      VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES) = 0

!  Single scatter correction: flux multiplier

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

!  Define DO_FOCORR_ALONE flag before input checks

      DO_FOCORR_ALONE = ( .NOT.DO_FULLRAD_MODE .AND. DO_FOCORR )

!  Check input
!  -----------

!    Major revision of I/O output list, 25 October 2012.
!     ---- Observational Geometry control, New, 25 October 2012
!     ---- Automatic setting of NBEAMS, N_USER_VZANGLES, N_USER_RELAZMS, DO_USER_VZANGLES

!  %% DO_FULLRAD_MODE argument added (First line). R. Spurr, 05 March 2013
!  %%   Needed to ensure MS-only output in all cases when flagged
!  Add Flag DO_WLADJUSTED_OUTPUT (Water-leaving output). 3/18/19 for Version 2.8.1. 

!    Major revision of I/O output list, 07 july 2016, for Version 2.8.
!      Inputs are arranged by type.

      CALL VLIDORT_CHECK_INPUT ( &
           DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,                & ! Input
           DO_THERMAL_TRANSONLY, DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_FLUORESCENCE, & ! Input
           DO_WATER_LEAVING, DO_TF_ITERATION, DO_WLADJUSTED_OUTPUT, DO_EXTERNAL_WLEAVE,      & ! Boolean 
           DO_FULLRAD_MODE,  DO_RAYLEIGH_ONLY, DO_SSCORR_USEFMAT, DO_DIRECT_BEAM,            & ! Boolean
           DO_FOCORR,  DO_OBSERVATION_GEOMETRY,  DO_USER_VZANGLES,  DO_DELTAM_SCALING,       & ! Boolean
           DO_FOCORR_EXTERNAL, DO_FOCORR_ALONE, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING,         & ! Boolean
           DO_REFRACTIVE_GEOMETRY, DO_CHAPMAN_FUNCTION, DO_PLANE_PARALLEL,                   & ! Boolean
           DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,       & ! Boolean
           DO_TOA_CONTRIBS, DO_ALL_FOURIER, DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3,  & ! Boolean
           TAYLOR_ORDER, NSTREAMS, NLAYERS, N_USER_LEVELS, TF_MAXITER,        & ! Integer
           N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS, N_USER_OBSGEOMS,      & ! Integer
           N_OUT_STREAMS, NLAYERS_NOMS, NLAYERS_CUTOFF, NGREEK_MOMENTS_INPUT, & ! Integer
           SZANGLES, USER_VZANGLES, USER_RELAZMS, USER_OBSGEOMS, USER_LEVELS, & ! Floating point
           OUT_ANGLES, TF_CRITERION, EARTH_RADIUS, GEOMETRY_SPECHEIGHT,       & ! Floating point
           HEIGHT_GRID, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,              & ! Floating point
           STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )               ! Exception handling

!  Exception handling

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
        VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
        VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
        RETURN
      ELSE IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = VLIDORT_WARNING
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
!  ########## Rob Change, 3/28/2011 ################
!  Program will execute - these outputs are set at the end.
!          VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
!          VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
!          VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
!  ########## Rob Change, 3/28/2011 ################
     ENDIF
     
!  4/26/19. Additional checks on the Isotropic-illumination input.
!  ---------------------------------------------------------------

!  9/25/19. Bug: VLIDORT_Out%Status%TS_STATUS_INPUTCHECK was not set

!  5/5/20. Version 2.8.1 Upgrades 
!   ==> Relax condition on Rayleigh only, and on FOCORR_NADIR
!              ( Planetary problem works for Aerosols and FOCORR_NADIR )
!   ==> Now, only fails for (dark-surface, Lambertian case, no thermal)

!  Here is the older code..........(Pre 5/5/20)
!      IF (DO_PLANETARY_PROBLEM .or. DO_ALBTRN_MEDIA(1) .or. DO_ALBTRN_MEDIA(2) ) then
!         IF ( .not.DO_LAMBERTIAN_SURFACE .or. DO_SURFACE_LEAVING .or. DO_THERMAL_EMISSION &
!              .or. DO_FOCORR .or. (ALBEDO .ne.zero) ) then
!            STATUS_INPUTCHECK = VLIDORT_SERIOUS
!            VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK ! 9/25/19 this was not set
!            VLIDORT_Out%Status%TS_NCHECKMESSAGES = 1
!            VLIDORT_Out%Status%TS_CHECKMESSAGES(1) = 'Media_problem input check not valid'
!            VLIDORT_Out%Status%TS_ACTIONS(1)       = 'Check thermal/Rayleigh/Lambertian flags for this option'
!            RETURN
!         ENDIF
!      ENDIF

!  Here is the New Code..........(5/5/20 Upgrade)

      IF ( DO_PLANETARY_PROBLEM .or. DO_ALBTRN_MEDIA(1) .or. DO_ALBTRN_MEDIA(2) ) then
         IF ( (.not. DO_LAMBERTIAN_SURFACE) .or. DO_SURFACE_LEAVING .or. DO_THERMAL_EMISSION &
              .or. (ALBEDO .ne.zero) .or. (DO_FOCORR.and.DO_FOCORR_OUTGOING) ) then
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
            VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK ! 9/25/19 this was not set
            VLIDORT_Out%Status%TS_NCHECKMESSAGES   = 2
            VLIDORT_Out%Status%TS_CHECKMESSAGES(1:2) = 'Media/Planetary problem: input check not valid'
            VLIDORT_Out%Status%TS_ACTIONS(1)       = 'Either: Turn off Thermal_Emission/FOCORR_OUTGOING/Surface_Leaving flags'
            VLIDORT_Out%Status%TS_ACTIONS(2)       = 'Or    : Make sure Lambertian surface flag and set albedo to zero'
            RETURN
         ENDIF
      ENDIF

!  Geometry adjustment
!  -------------------

!  Coding removed for Version 2.8

!  Adjust surface condition

!      ADJUST_SURFACE = .FALSE.
!      IF ( DO_SSCORR_OUTGOING ) THEN
!        IF (HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) THEN
!         ADJUST_SURFACE = .TRUE.
!        ENDIF
!      ENDIF

!  Perform adjustment

!      MODIFIED_ERADIUS = EARTH_RADIUS + GEOMETRY_SPECHEIGHT
!mick hold - 9/26/2012
!      IF ( DO_SOLAR_SOURCES ) THEN

!      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
!        CALL MULTI_OUTGOING_ADJUSTGEOM                                &
!         ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS,         & ! Input
!           N_USER_VZANGLES,   N_SZANGLES,   N_USER_RELAZMS,           & ! Input
!           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE,    & ! Input
!           USER_VZANGLES, SZANGLES, USER_RELAZMS,                     & ! Input
!           USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST,& ! Output
!           FAIL, MESSAGE, TRACE_1 )                                     ! Output
!      ELSE
!        CALL OBSGEOM_OUTGOING_ADJUSTGEOM                               &
!         ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS,          & ! Input
!           N_USER_VZANGLES,   N_SZANGLES,   N_USER_RELAZMS,            & ! Input
!           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE,     & ! Input
!           USER_VZANGLES, SZANGLES, USER_RELAZMS,                      & ! Input
!           USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST, & ! Output
!           FAIL, MESSAGE, TRACE_1 )                                      ! Output
!      ENDIF

!      ELSE
!        CALL LOSONLY_OUTGOING_ADJUSTGEOM                           &
!         ( MAX_USER_VZANGLES, N_USER_VZANGLES,                     & ! Input
!           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE, & ! Input
!           USER_VZANGLES,                                          & ! Input
!           USER_VZANGLES_ADJUST,                                   & ! Output
!           FAIL, MESSAGE, TRACE_1 )                                  ! Output
!      ENDIF

!  Update exception handling. October 2010 2p4RTC

!      IF ( FAIL ) THEN
!        TRACE_2 = ' Failure in multi_outgoing_adjustgeom'
!        TRACE_3 = ' ** VLIDORT_MASTER'
!        STATUS_CALCULATION = VLIDORT_SERIOUS
!        VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
!        VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
!        VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!        VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
!        VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!        RETURN
!      ENDIF

!  write input variables
!  ---------------------

!  Version 2.8 revision 7/8/16, remove QUAD_OUTPUT and CLASSICAL SOLUTION, revise input list.

!  open file, call standard input write, close file

      IF ( DO_WRITE_INPUT ) THEN
        IUNIT = VLIDORT_INUNIT
        OPEN(IUNIT,FILE=INPUT_WRITE_FILENAME,STATUS='UNKNOWN')
        CALL VLIDORT_WRITEINPUT ( &
          IUNIT, DO_USER_VZANGLES, DO_NO_AZIMUTH,                       &
          DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, DO_RAYLEIGH_ONLY,  &
          DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_LAMBERTIAN_SURFACE,    &
          DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_DIRECT_BEAM,     &
          NSTREAMS, NLAYERS, N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS,    &
          N_USER_LEVELS, NGREEK_MOMENTS_INPUT, VLIDORT_ACCURACY, FLUX_FACTOR )
        CLOSE(IUNIT)
      ENDIF

!  Get derived inputs
!  ==================

!  Miscellaneous and layer input.
!  Version 2.8 revision 7/8/16, I/O list. USER_VZANGLES_ADJUST no longer present
!  Version 2.8 revision 3/1/17, I/O list for FOCORR flags
!mick fix 9/19/2017 - added DO_FOCORR_NADIR & DO_FOCORR_OUTGOING to input

      CALL VLIDORT_DERIVE_INPUT ( &
        DO_FULLRAD_MODE, DO_FOCORR, DO_FOCORR_ALONE,                           & ! Input Boolean
        DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_SSCORR_USEFMAT,                & ! Input Boolean
        DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_RAYLEIGH_ONLY,            & ! Input Boolean
        DO_UPWELLING, DO_DNWELLING, DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY, & ! Input Boolean
        DO_THERMAL_TRANSONLY, DO_SPECIALIST_OPTION_2, DO_SOLUTION_SAVING,      & ! Input Boolean
        DO_BVP_TELESCOPING, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_NO_AZIMUTH, & ! Input Boolean
        DO_DOUBLE_CONVTEST, DO_ALL_FOURIER, DO_DIRECT_BEAM, DO_DBCORRECTION,   & ! Input Boolean
        NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, N_SZANGLES,  NLAYERS_NOMS,  & ! Input Integer
        NGREEK_MOMENTS_INPUT, N_USER_RELAZMS, N_USER_VZANGLES,  N_OUT_STREAMS, & ! Input Integer
        SZANGLES, USER_VZANGLES, USER_LEVELS,                                  & ! Input Floating point
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,            & ! Input Floating point
        COS_SZANGLES, SIN_SZANGLES, QUAD_STREAMS,                              & ! OUTPUT
        QUAD_WEIGHTS,QUAD_STRMWTS, QUAD_HALFWTS, QUAD_SINES, QUAD_ANGLES,      & ! OUTPUT
        DO_MSMODE_VLIDORT, NMOMENTS, NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, & ! OUTPUT
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, NBEAMS, NPARTICSOLS, NSTOKES_SQ,         & ! OUTPUT
        FLUXVEC, DMAT, MUELLER_INDEX, GREEKMAT_INDEX, DO_REAL_EIGENSOLVER,     & ! OUTPUT
        BVP_REGULAR_FLAG, LAYER_MAXMOMENTS, DO_LAYER_SCATTERING,               & ! OUTPUT
        N_CONVTESTS, N_CONV_STREAMS, N_DIRECTIONS, WHICH_DIRECTIONS,           & ! OUTPUT
        USER_STREAMS, USER_SINES, USER_SECANTS, PARTLAYERS_OUTFLAG,            & ! OUTPUT
        PARTLAYERS_OUTINDEX, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,           & ! OUTPUT
        DO_PARTLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,   & ! OUTPUT
        N_LAYERSOURCE_UP, N_LAYERSOURCE_DN, N_ALLLAYERS_UP, N_ALLLAYERS_DN,    & ! OUTPUT
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DFLUX, TAUGRID_INPUT,          & ! OUTPUT
        STATUS_SUB, MESSAGE )

!  Set thermal MS flag
!mick fix 9/19/2017 - changed def of DO_MSMODE_THERMAL flag based on implementation
!                     of new FO code.  When DO_FOCORR set, both solar AND THERMAL
!                     direct now come from the FO code.

      !DO_MSMODE_THERMAL = (.not.DO_FULLRAD_MODE) .and. ( DO_SURFACE_EMISSION .and.DO_THERMAL_EMISSION )
      DO_MSMODE_THERMAL = DO_MSMODE_VLIDORT .and. ( DO_SURFACE_EMISSION .and. DO_THERMAL_EMISSION )

!  Exception handling

!  Rob Fix 2011. Output from DERIVE_INPUTS are checks, not execution failures
!  Old code---------
!      IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
!        TRACE_1 = ''
!        TRACE_2 = 'Derive_Input Call in VLIDORT_L_MASTER'
!        TRACE_3 = ' ** VLIDORT_L_MASTER'
!        STATUS_INPUTCHECK = VLIDORT_WARNING
!        VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
!        VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
!        VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!        VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
!        VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!      ENDIF
!  New code-----------
      IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = VLIDORT_WARNING
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        NCHECKMESSAGES = NCHECKMESSAGES + 1
        CHECKMESSAGES(NCHECKMESSAGES)  = TRIM(ADJUSTL(MESSAGE))
        ACTIONS(NCHECKMESSAGES) = ' Action taken in VLIDORT_DERIVE_INPUT to set internal default'
      ENDIF

!  Rob Fix, 3/17/15. Serious Error check on use of FO calculation, must exit
!   --> Plan for Version 2.8 is to relax this condition
!      IF (DO_FO_CALC .AND. DO_PARTLAYERS) THEN
!        STATUS_INPUTCHECK = VLIDORT_SERIOUS
!        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
!        NCHECKMESSAGES = NCHECKMESSAGES + 1
!        CHECKMESSAGES(NCHECKMESSAGES) = ' Internal SS calculation using FO code: ONLY LAYER-BOUNDARY output (no partials)'
!        ACTIONS(NCHECKMESSAGES)       = ' Set USER-LEVEL output only for layer boundaries'
!        VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
!        VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
!        VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
!        RETURN
!      ENDIF

!  Chapman function calculation of TAUTHICK_INPUT
!  ----------------------------------------------

!  Calling statement revised I/O, Version 2.8 7/7/16.

!mick fix - comment out 1st IF
!mick fix 9/19/2017 - added input N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES
!                   - added output PARTIAL_CHAPFACS
!                   - moved the call to VLIDORT_CHAPMAN from before VLIDORT_DERIVE_INPUT
!                     to here to resolve an I/O contradiction

        !IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_CHAPMAN_FUNCTION ) THEN

          CALL VLIDORT_CHAPMAN ( &
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,            & !  Input
            NLAYERS, N_SZANGLES, FINEGRID, SZANGLES,              & !  Input
            N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, & !  Input
            EARTH_RADIUS, RFINDEX_PARAMETER,                      & !  Input
            HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID,         & !  Input
            CHAPMAN_FACTORS, PARTIAL_CHAPFACS, SZA_LOCAL_INPUT,   & !  Output
            SUN_SZA_COSINES, FAIL, MESSAGE, TRACE_1 )               !  Output

          IF (FAIL) THEN
            TRACE_2 = 'Direct call in VLIDORT_MASTER'
            TRACE_3 = ' ** VLIDORT_MASTER '
            STATUS_CALCULATION = VLIDORT_SERIOUS
            VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
            VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
          ENDIF

        ENDIF
      !ENDIF

!  Single scatter correction: flux multiplier. Change to DO_SSCORR.
!      IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN

      IF ( DO_FOCORR ) THEN
        SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4
      ELSE
        SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4
      ENDIF

!  if there's no azimuth dependence, just do one value in azimuth loop

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_N_USERAZM = 1
      ELSE
        LOCAL_N_USERAZM = N_USER_RELAZMS
      ENDIF

!  Number of sources

      IF ( DO_SOLAR_SOURCES ) THEN
        NSOURCES = N_SZANGLES
      ELSE
        NSOURCES = 1
        NBEAMS   = 1
      ENDIF

!  save some offsets for indexing geometries
!   This section revised for the Observational Geometry option

      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        N_VIEWING    = N_USER_VZANGLES * LOCAL_N_USERAZM
        N_GEOMETRIES = NSOURCES * N_VIEWING
        DO IBEAM = 1, NBEAMS
          SZA_OFFSETS(IBEAM) = N_VIEWING * ( IBEAM - 1 )
          DO UM = 1, N_USER_VZANGLES
            VZA_OFFSETS(IBEAM,UM) = SZA_OFFSETS(IBEAM) + LOCAL_N_USERAZM * (UM - 1)
          END DO
        END DO
      ELSE
        N_VIEWING    = N_USER_OBSGEOMS
        N_GEOMETRIES = N_USER_OBSGEOMS
        SZA_OFFSETS = 0 ; VZA_OFFSETS = 0
      ENDIF

!  mick add 3/25/2015 - geometry order check (leave for debug)
!      WRITE(*,*)
!      WRITE(*,'(1X,A,I2)') 'LOCAL_N_USERAZM = ',LOCAL_N_USERAZM
!      WRITE(*,'(1X,A,I2)') 'N_USER_VZANGLES = ',N_USER_VZANGLES
!      WRITE(*,'(1X,A,I2)') 'NBEAMS          = ',NBEAMS
!      WRITE(*,*)
!      DO IBEAM = 1, NBEAMS
!        DO UM = 1, N_USER_VZANGLES
!          DO UA = 1, LOCAL_N_USERAZM
!            V = VZA_OFFSETS(IBEAM,UM) + UA
!            WRITE(*,'(1X,A,I2,3(A,F6.2))') 'FOR V = ',V,&
!              '  RAA(UA) = ',USER_RELAZMS(UA), &
!              '  UZA(UM) = ',USER_VZANGLES(UM),&
!              '  SZA(IB) = ',SZANGLES(IBEAM)
!          ENDDO
!        ENDDO
!      ENDDO
!      STOP

!  #################
!  Set up operations
!  #################

!  Setups for Fourier = 0
!  ======================

!  Each call is foolowed by a packing routine

!  1.  VLIDORT_MISCSETUPS  : Delta-M, average-secant formulation, transmittances
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

      CALL VLIDORT_MISCSETUPS ( &
        DO_SOLAR_SOURCES, DO_DELTAM_SCALING, DO_PLANE_PARALLEL,        & ! Input flags
        DO_REFRACTIVE_GEOMETRY,                                        & ! Input flags
        DO_PARTLAYERS, DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY,      & ! Input flags
        DO_SOLUTION_SAVING, DO_SPECIALIST_OPTION_3, DO_TOA_CONTRIBS,   & ! Input flags
        NSTOKES, NLAYERS, NSTREAMS, N_USER_VZANGLES, N_USER_LEVELS,    & ! Input Numbers
        NBEAMS, NMOMENTS, N_PARTLAYERS, NLAYERS_CUTOFF,                & ! Input Numbers
        MUELLER_INDEX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,         & ! Input level/Stokes indices
        PARTLAYERS_OUTFLAG,  PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,   & ! Input Partlayers
        QUAD_STREAMS, USER_SECANTS, COS_SZANGLES, SUN_SZA_COSINES,     & ! Input streams, SZA cosines
        OMEGA_TOTAL_INPUT, DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT,    & ! Input Optical
        TAUGRID_INPUT, CHAPMAN_FACTORS, PARTIAL_CHAPFACS,              & ! Input Chapman
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
        CUMTRANS, ITRANS_USERM )                                         ! Output PREPTRANS (auxiliary)       
      
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
                DO_UPWELLING, DO_DNWELLING,                                   & ! Input flags
                NLAYERS, NBEAMS, N_USER_VZANGLES, TAYLOR_ORDER, N_PARTLAYERS, & ! Input numbers
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,  & ! Input Level control
                USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                       & ! Input optical, streams   
                BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,      & ! Input solar beam
                T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, ITRANS_USERM,       & ! Input User-stream Trans.
                EMULT_HOPRULE, SIGMA_M, SIGMA_P,                              & ! Output Multipliers
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )                  ! Output Multipliers
          ELSE
            CALL EMULT_MASTER_OBSGEO ( &
                DO_UPWELLING, DO_DNWELLING,                                  & ! Input flags
                NLAYERS, NBEAMS, TAYLOR_ORDER, N_PARTLAYERS,                 & ! Input numbers
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
!           ( NSTOKES, NSTREAMS, NLAYERS, NBEAMS, NMOMENTS, NSTREAMS_2,           & ! Input
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
!         do ibeam = 1, nbeams
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
!         TRANS_ATMOS(1:nbeams) = FLUX_DIFFUSE_TF(1:nbeams,1) / flux_factor / local_csza(nlayers,1:nbeams)
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
!          NSTOKES, NSTREAMS, NLAYERS, NBEAMS, NMOMENTS, NSTREAMS_2,                 & ! Input Numbers
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
!         do ibeam = 1, nbeams
!            SLTERM_ISOTROPIC(1,IBEAM) = SLTERM_ISOTROPIC(1,IBEAM) * TRANS_ATMOS_FINAL(IBEAM)
!         enddo
!  Scale the Non-isotropic terms, 24 December 2015, if flagged
!         if ( .not. DO_SL_ISOTROPIC ) THEN
!           do ibeam = 1, nbeams
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

      !IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_USER_VZANGLES ) THEN
          IF ( DO_FOCORR ) THEN

!  Call to interface. Updated, 17 September 2016.
!   -- 4/9/19. Added output FO Surface-leaving assignation + cumulative transmittances.
!              Added input, water-leaving control

            CALL VFO_MASTER_INTERFACE ( &
                DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                         & ! Input Sources flags
                DO_PLANE_PARALLEL, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_DELTAM_SCALING,          & ! Input SS control flags
                DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY, DO_PARTLAYERS,                 & ! input Model control flags
                DO_SSCORR_USEFMAT, DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_WATER_LEAVING,     & ! Input Optical/Surface flags
                DO_SL_ISOTROPIC, NSTOKES, NLAYERS, NFINELAYERS, NGREEK_MOMENTS_INPUT,               & ! Input numbers
                N_SZANGLES, SZANGLES, N_USER_VZANGLES, USER_VZANGLES, N_USER_RELAZMS, USER_RELAZMS, & ! Input geometry
                N_USER_LEVELS, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, N_PARTLAYERS,                & ! Input levels  control
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,    & ! Input partial control
                EARTH_RADIUS, HEIGHT_GRID, SS_FLUX_MULTIPLIER, FLUXVEC,                             & ! Input Flux/Heights
                DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                         & ! Inputs (Optical - Regular)
                DELTAU_VERT, FMATRIX_UP, FMATRIX_DN, TRUNC_FACTOR, THERMAL_BB_INPUT,                & ! Inputs (Optical - Regular)
                ALBEDO, EXACTDB_BRDFUNC, SURFBB, USER_EMISSIVITY,                                   & ! Inputs (Optical - Surface)
                SLTERM_ISOTROPIC, SLTERM_USERANGLES,                                                & ! Inputs (Optical - Surface)
                FO_STOKES_SS, FO_STOKES_DB, FO_STOKES_DTA, FO_STOKES_DTS,                           & ! Output Stokes vectors
                FO_STOKES_ATMOS, FO_STOKES_SURF, FO_STOKES, FO_CUMTRANS, FO_SLTERM,                 & ! Output Stokes vectors
                FAIL, MESSAGE, TRACE_1, TRACE_2 )                                                     ! Exception handling

!  Exception handling

            IF ( FAIL ) THEN
               TRACE_3 = 'VFO_MASTER_INTERFACE failed, VLIDORT_MASTER'
               STATUS_CALCULATION = VLIDORT_SERIOUS
               VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
               VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
               VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
               VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
               VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
               RETURN
            ENDIF

!  Copy FO results to VLIDORT arrays
!mick note 9/19/2017 - Important! STOKES_SS & STOKES_DB contain BOTH solar AND thermal
!                      direct terms when computing in the crossover region!

            STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)   = &
              FO_STOKES_ATMOS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)
            STOKES_DB(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)     = &
              FO_STOKES_SURF(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)

          ENDIF

!  End user-angle and solar-source if blocks

        ENDIF
      !ENDIF

!  ####################
!   MAIN FOURIER LOOP
!  ####################

!  Initialise Fourier loop
!  =======================

!  Set Number of Fourier terms (NMOMENTS = Maximum).
!    ( Starting from 0 = Fundamental )

      SAVE_DO_NO_AZIMUTH  = DO_NO_AZIMUTH
      LOCAL_DO_NO_AZIMUTH = DO_NO_AZIMUTH

!  No azimuth dependency for following cases
!    Thermal-only treatment is isotropic, no azimuth

!      IF ( DO_TRANSMITTANCE_ONLY ) THEN
!        LOCAL_DO_NO_AZIMUTH = .TRUE.
!      ENDIF

!mick mod 3/30/2015 - consolidated if conditions
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

      LOCAL_ITERATION   = .TRUE.
      FOURIER = -1
      TESTCONV          = 0

!  set up solar beam flags. Required in all cases.
!   ---Even required for the thermal.....

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
            NSTOKES, NSTREAMS, NLAYERS, NBEAMS, N_USER_VZANGLES, N_USER_LEVELS, N_THERMAL_COEFFS, NMOMENTS, NSTREAMS_2, & !Input
            NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2, N_PARTLAYERS, N_DIRECTIONS, TAYLOR_ORDER, & !Input
            N_ALLLAYERS_UP, N_ALLLAYERS_DN, FLUX_FACTOR, FLUXVEC, COS_SZANGLES, SZA_LOCAL_INPUT, SUN_SZA_COSINES, & !Input
            USER_STREAMS, USER_SECANTS, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS,                   & !Input
            MUELLER_INDEX, DMAT, BVP_REGULAR_FLAG, LOCAL_UM_START, WHICH_DIRECTIONS, SOLARBEAM_BOATRANS,          & !Input
            UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & !Input
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DFLUX, SURFBB, EMISSIVITY, USER_EMISSIVITY,                   & !Input
            Misc, Therm, Mult,                                                                                    & !Input
            ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,  & !Input
            PIMM_11, PIMM_KM, BVTEL_FOURIER, DO_INCLUDE_THERMEMISS, DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS,     & !Output
            STOKES_F, MEANST_DIFFUSE, FLUX_DIFFUSE, DNMEANST_DIRECT, DNFLUX_DIRECT, MS_CONTRIBS_F,                & !Output MAIN
            TRANS_ATMOS_FINAL, TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,                 & !Output 4/26/19
            STATUS_SUB, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )                                                        !Output Status

!  Useful debug
!        if ( FOURIER.eq.0) write(*,*)FOURIER,TRANS_ATMOS_FINAL(1),STOKES_F(1:2,1,1,1,1),STOKES_F(2,1,1,1,2)
         
!  error handling

         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            TRACE_3 = ' Called by VLIDORT_MASTER '
            STATUS_CALCULATION = VLIDORT_SERIOUS
            VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
            VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
         ENDIF

!  4/9/19. FO-DB adjustment for Fourier 0, using self-consistent water-leaving term
!          Intensity only, Surface-leaving not polarized.         

         if ( FOURIER .eq. 0 .and. DO_FOCORR .and. DO_WATER_LEAVING ) then
            O1 = 1
            if ( do_OBSERVATION_GEOMETRY ) THEN
               DO IB = 1, NBEAMS
                  SLTERM_LOCAL = FO_SLTERM(O1,ib) * TRANS_ATMOS_FINAL(ib)
                  do uta = 1, n_user_levels
                     CUMSOURCE_DB = FO_CUMTRANS(UTA,IB) * SLTERM_LOCAL
                     STOKES_DB(UTA,IB,O1) = STOKES_DB(UTA,IB,O1) + CUMSOURCE_DB
                  enddo
               enddo
            else
               DO IB = 1, NBEAMS ; DO UM = 1, N_USER_VZANGLES
                  g1 = VZA_OFFSETS(IB,UM) + 1 ; g2 = VZA_OFFSETS(IB,UM) + n_user_relazms
                  SLTERM_LOCAL = FO_SLTERM(O1,g1) * TRANS_ATMOS_FINAL(ib)
                  do uta = 1, n_user_levels
                     do g = g1, g2
                        CUMSOURCE_DB = FO_CUMTRANS(UTA,g) * SLTERM_LOCAL
                        STOKES_DB(UTA,g,O1) = STOKES_DB(UTA,g,O1) + CUMSOURCE_DB
                     enddo
                  enddo
               enddo; enddo
            Endif
         endif

!  Useful debug
!           if ( FOURIER.eq.0) write(*,*)'FO_Rev',STOKES_SS(1,1,1,1),STOKES_DB(1:2,1:2,1)

!  Output for WLADJUSTED Water-Leaving . Introduced 4/22/19 for Version 3.8.1
      ! Sleave Results need to be modified from their original inputs.
      ! Debug results - USE WITH CAUTION. Note the preliminary zeroing to avoid unassigned arrays.

         if ( FOURIER .eq. 0 .and. DO_WATER_LEAVING .and. DO_WLADJUSTED_OUTPUT ) then
            VLIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC = zero
            VLIDORT_Out%WLOut%TS_WLADJUSTED_DIRECT    = zero
            VLIDORT_Out%WLOut%TS_WLADJUSTED_F_Ords_0  = zero
            VLIDORT_Out%WLOut%TS_WLADJUSTED_F_User_0  = zero
            DO IB = 1, NBEAMS
               TFACTOR = TRANS_ATMOS_FINAL(ib) ; O1 = 1
               VLIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC(O1,IB) = SLTERM_ISOTROPIC(O1,IB) * TFACTOR
               VLIDORT_Out%WLOut%TS_WLADJUSTED_F_Ords_0(0,O1,1:NSTREAMS,IB) = TFACTOR * SLTERM_F_0(0,O1,1:NSTREAMS,IB)
               IF ( DO_USER_VZANGLES ) THEN
                  VLIDORT_Out%WLOut%TS_WLADJUSTED_DIRECT(O1,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IB) = &
                                          TFACTOR * SLTERM_USERANGLES(O1,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IB)
                  VLIDORT_Out%WLOut%TS_WLADJUSTED_F_User_0(0,O1,1:N_USER_VZANGLES,IB) = &
                                          TFACTOR * USER_SLTERM_F_0(0,O1,1:N_USER_VZANGLES,IB)
               ENDIF
            ENDDO
         ENDIF

! #################### INNOVATIONS 5/5/20 ##############
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
! #################### INNOVATIONS 5/5/20 ##############

!  Fourier summation and Convergence examination
!  ---------------------------------------------

!mick fix 3/30/2015 - added if condition and moved azimuth block from
!                     before call to VLIDORT_FOURIER to here

!  Begin convergence if block

         IF ( .NOT.DO_MVOUT_ONLY ) THEN

!  azimuth cosine/sine factors
!    - Use of adjusted geometries retired for Version 2.8
!mick eff 3/22/2017 - both IF and ELSE sections

          IF ( FOURIER .GT. 0 ) THEN
            DFC = DBLE(FOURIER)
            IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
              DO UA = 1, LOCAL_N_USERAZM
                AZM_ARGUMENT = USER_RELAZMS(UA) * DFC * DEG_TO_RAD
                AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,1)  = COS(AZM_ARGUMENT)
                AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,3)  = SIN(AZM_ARGUMENT)
                AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,2)  = AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,1)
                AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,4)  = AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,3)
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

              IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
                CALL VLIDORT_CONVERGE ( &
                  DO_FOCORR, DO_FOCORR_EXTERNAL, DO_FOCORR_ALONE, DO_DBCORRECTION,           & ! Input flags
                  DO_DOUBLE_CONVTEST, DO_RAYLEIGH_ONLY, DO_UPWELLING,                        & ! Input flags
                  DO_TOA_CONTRIBS, DO_ALL_FOURIER,  DO_NO_AZIMUTH,                           & ! Input flags
                  FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_USER_RELAZMS, N_USER_LEVELS, & ! Input control numbers
                  N_CONVTESTS, LOCAL_UM_START, LOCAL_N_USERAZM, N_DIRECTIONS, N_OUT_STREAMS, & ! Input derived numbers
                  WHICH_DIRECTIONS, VZA_OFFSETS, VLIDORT_ACCURACY, AZMFAC,                   & ! Input bookkeeping
                  STOKES_SS, STOKES_DB, SS_CONTRIBS, STOKES_F, MS_CONTRIBS_F,                & ! Input Fourier and SS fields
                  STOKES, FOURIER_SAVED, CONTRIBS, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) ) ! Output Stokes and diagnostics
              ELSE
                CALL VLIDORT_CONVERGE_OBSGEO ( &
                  DO_FOCORR, DO_FOCORR_EXTERNAL, DO_FOCORR_ALONE, DO_DBCORRECTION,           & ! Input flags
                  DO_DOUBLE_CONVTEST, DO_RAYLEIGH_ONLY, DO_UPWELLING,                        & ! Input flags
                  DO_TOA_CONTRIBS, DO_ALL_FOURIER,  DO_NO_AZIMUTH,                           & ! Input flags
                  FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS,  N_USER_LEVELS,                & ! Input control numbers
                  N_CONVTESTS, N_DIRECTIONS,  WHICH_DIRECTIONS, VLIDORT_ACCURACY, AZMFAC,    & ! Input bookkeeping
                  STOKES_SS, STOKES_DB, SS_CONTRIBS, STOKES_F, MS_CONTRIBS_F,                & ! Input Fourier and SS fields
                  STOKES, FOURIER_SAVED, CONTRIBS, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) ) ! Output Stokes and diagnostics
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
          RUNIT, DO_FULLRAD_MODE, DO_FOCORR_NADIR,   &
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
          STOKES, MEANST_DIFFUSE,             &
          FLUX_DIFFUSE, FOURIER_SAVED )

        CLOSE(RUNIT)
      ENDIF

!  Geophysical input (scenario) write
!  ----------------------------------

      IF ( DO_WRITE_SCENARIO ) THEN
        SUNIT = VLIDORT_SCENUNIT
        OPEN(SUNIT,FILE=SCENARIO_WRITE_FILENAME,STATUS='UNKNOWN')

        CALL VLIDORT_WRITESCEN ( &
          SUNIT, DO_DELTAM_SCALING, NSTREAMS, &
          NLAYERS, NGREEK_MOMENTS_INPUT,      &
          N_SZANGLES, SZANGLES,               &
          N_USER_RELAZMS, USER_RELAZMS,       &
          N_USER_VZANGLES, USER_VZANGLES,     &
          N_USER_LEVELS, USER_LEVELS,         &
          OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,  &
          DO_LAMBERTIAN_SURFACE, ALBEDO, &
          DO_THERMAL_EMISSION, N_THERMAL_COEFFS,    &
          DO_SURFACE_EMISSION, SURFBB, &
          QUAD_STREAMS, QUAD_WEIGHTS,  &
          QUAD_ANGLES, DO_NO_AZIMUTH,  &
          NMOMENTS, GREEKMAT_INDEX,    &
          TAUGRID_INPUT, OMEGA_TOTAL,  &
          GREEKMAT_TOTAL, TAUGRID )

        CLOSE(SUNIT)
      ENDIF

!  ========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (IN/OUT variables)
!  ========================================================

!  FOCORR Booleans reorganized for Version 2.8
!mick mod 9/19/2017 - reordered FO variables to conform to newly modified input type structure
!                   - DO_FOCORR_ALONE now defined internally

      VLIDORT_ModIn%MBool%TS_DO_FOCORR              = DO_FOCORR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL     = DO_FOCORR_EXTERNAL
      !VLIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE        = DO_FOCORR_ALONE
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR        = DO_FOCORR_NADIR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING     = DO_FOCORR_OUTGOING

      VLIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION   = DO_SSCORR_TRUNCATION
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT      = DO_SSCORR_USEFMAT

!  Additional Control for Externalized input (SLEAVE). Introduced 3/18/19 for Version 2.8.1

      VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE     = DO_EXTERNAL_WLEAVE
      
!  Solar control

      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = DO_SOLAR_SOURCES
      VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = DO_REFRACTIVE_GEOMETRY
      VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = DO_CHAPMAN_FUNCTION

      !VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY      = DO_ISOTROPIC_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH          = DO_NO_AZIMUTH
      !VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER         = DO_ALL_FOURIER

!  Performance control

      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY       = DO_RAYLEIGH_ONLY
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING      = DO_DELTAM_SCALING
      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST     = DO_DOUBLE_CONVTEST
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING     = DO_SOLUTION_SAVING
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING     = DO_BVP_TELESCOPING

!  RT Model control

      VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES       = DO_USER_VZANGLES
      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = DO_ADDITIONAL_MVOUT
      VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY          = DO_MVOUT_ONLY
      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY   = DO_THERMAL_TRANSONLY
      VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY= DO_OBSERVATION_GEOMETRY

!  Modified control inputs

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT   = NGREEK_MOMENTS_INPUT

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

      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS)             = USER_LEVELS(1:N_USER_LEVELS)

      VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS     = N_USER_OBSGEOMS
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,:) = &
        USER_OBSGEOMS(1:N_USER_OBSGEOMS,:)

!  Modified Chapman function inputs

      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS        = EARTH_RADIUS

!  Modified optical inputs

      VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS) = OMEGA_TOTAL_INPUT(1:NLAYERS)

!  ==========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (pure OUT variables)
!  ==========================================================

!  Radiances and fluxes
!    Output direct Mean/Flux added 17 May 2012

      VLIDORT_Out%Main%TS_STOKES(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)    = &
        STOKES(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)

      VLIDORT_Out%Main%TS_MEANST_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        MEANST_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)
      VLIDORT_Out%Main%TS_FLUX_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        FLUX_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)
      VLIDORT_Out%Main%TS_DNMEANST_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = &
        DNMEANST_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)
      VLIDORT_Out%Main%TS_DNFLUX_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = &
        DNFLUX_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)

!  Verify TF calculation. Debug 2/3/16. Worked out.
!      do ibeam = 1, n_szangles
!         write(*,*)'MK2 Result = ',FLUX_DIFFUSE(2,1,1,2)
!         write(778,'(a,i2,1p3e15.6,2x)')'Original', ibeam,MEANST_DIFFUSE(2,ibeam,1,2),FLUX_DIFFUSE(2,ibeam,1,2),&
!         DNFLUX_DIRECT(2,ibeam,1)
!      enddo

!  4/26/19. Media properties output.

      VLIDORT_Out%Main%TS_ALBMED_USER(1:NSTOKES,1:N_USER_VZANGLES) = ALBMED_USER(1:NSTOKES,1:N_USER_VZANGLES)
      VLIDORT_Out%Main%TS_TRNMED_USER(1:NSTOKES,1:N_USER_VZANGLES) = TRNMED_USER(1:NSTOKES,1:N_USER_VZANGLES)
      VLIDORT_Out%Main%TS_ALBMED_FLUXES(1:NSTOKES,1:2) = ALBMED_FLUXES(1:NSTOKES,1:2)
      VLIDORT_Out%Main%TS_TRNMED_FLUXES(1:NSTOKES,1:2) = TRNMED_FLUXES(1:NSTOKES,1:2)

!  4/28/19. Planetary problem output
!  ---------------------------------

!  5/5/20. Version 2.8.1 Upgrades 
!    ==> Revised code to use only the first index of TRANSBEAM ==> makes the Q-problem valid

!  Here is the Old code (pre 5/5/20 Upgrade) 
!       IF ( DO_PLANETARY_PROBLEM ) THEN
!         VLIDORT_Out%Main%TS_PLANETARY_SBTERM = TRNMED_FLUXES(1,2)
!         IF ( DO_OBSERVATION_GEOMETRY ) THEN
!            DO IB =1, N_GEOMETRIES ; DO O1 = 1, NSTOKES
!               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,IB) = TRANSBEAM(O1,IB) * TRNMED_USER(O1,IB) / PIE
!            ENDDO ; ENDDO
!         ELSE
!            DO IB =1, NBEAMS ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
!               OFF = VZA_OFFSETS(IB,UM) ; TRANS = TRANSBEAM(O1,IB) * TRNMED_USER(O1,UM) / PIE
!               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,OFF+1:OFF+N_USER_RELAZMS) = TRANS
!            ENDDO ; ENDDO ; ENDDO
!         ENDIF   
!      ENDIF

!  Here is the New code (5/5/20 Upgrade)  

      IF ( DO_PLANETARY_PROBLEM ) THEN
         VLIDORT_Out%Main%TS_PLANETARY_SBTERM = TRNMED_FLUXES(1,2)
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            DO IB =1, N_GEOMETRIES ; DO O1 = 1, NSTOKES
               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,IB) = TRANSBEAM(1,IB) * TRNMED_USER(O1,IB) / PIE
            ENDDO ; ENDDO
         ELSE
            DO IB =1, NBEAMS ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
               OFF = VZA_OFFSETS(IB,UM) ; TRANS = TRANSBEAM(1,IB) * TRNMED_USER(O1,UM) / PIE
               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,OFF+1:OFF+N_USER_RELAZMS) = TRANS
            ENDDO ; ENDDO ; ENDDO
         ENDIF   
      ENDIF

!  new 12 March 2012
!   IF SS results already available, no need to copy them !

      IF ( .NOT. DO_FOCORR_EXTERNAL ) THEN
         VLIDORT_Sup%SS%TS_STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:) = &
           STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)
         VLIDORT_Sup%SS%TS_STOKES_DB(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES) = &
           STOKES_DB(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
      ENDIF

!  Additional Control for Externalized input (SLEAVE). Introduced 3/18/19 for Version 2.8.1
      ! Sleave Results May have been modified from their external inputs.
      ! Allows you to come out with modified SLEAVE.

      IF ( DO_EXTERNAL_WLEAVE ) THEN
        VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(1:NSTOKES,1:N_SZANGLES) = &
                SLTERM_ISOTROPIC(1:NSTOKES,1:N_SZANGLES)
        VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0(0:NMOMS,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES) = &
                SLTERM_F_0(0:NMOMS,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES)
        IF ( DO_USER_VZANGLES ) THEN
          VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
                SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
          VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(0:NMOMS,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES) = &
            USER_SLTERM_F_0(0:NMOMS,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES)
        ENDIF
      ENDIF
      
!  Bookkeeping

      VLIDORT_Out%Main%TS_FOURIER_SAVED(1:N_SZANGLES) = FOURIER_SAVED(1:N_SZANGLES)
      VLIDORT_Out%Main%TS_N_GEOMETRIES = N_GEOMETRIES
      VLIDORT_Out%Main%TS_SZA_OFFSETS(1:N_SZANGLES)   = SZA_OFFSETS(1:N_SZANGLES)
      VLIDORT_Out%Main%TS_VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES) = VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES)

!  Solar Beam Transmittance to BOA
!  rob fix 11/17/2014, for diagnostic use only

      VLIDORT_Out%Main%TS_SOLARBEAM_BOATRANS (1:NBEAMS) = SOLARBEAM_BOATRANS (1:NBEAMS)

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
        
      CALL VLIDORT_WRITE_STD_INPUT ( DO_DEBUG_WRITE, &
        DO_SOLAR_SOURCES,  DO_TOAFLUX, DO_BOAFLUX, DO_ALBTRN_MEDIA,                  & ! Sources
        DO_THERMAL_EMISSION, DO_THERMAL_TRANSONLY, DO_SURFACE_EMISSION,                                & ! Sources
        DO_FULLRAD_MODE,  DO_FOCORR,           DO_FOCORR_EXTERNAL,   DO_PLANETARY_PROBLEM,             & ! FOCORR/Planetary
        DO_FOCORR_NADIR,  DO_FOCORR_OUTGOING,  DO_SSCORR_TRUNCATION, DO_SSCORR_USEFMAT,                & ! FOCORR
        DO_RAYLEIGH_ONLY, DO_DELTAM_SCALING,   DO_DOUBLE_CONVTEST, DO_SOLUTION_SAVING,     DO_BVP_TELESCOPING, & ! Performance
        DO_UPWELLING,     DO_DNWELLING,        DO_PLANE_PARALLEL,  DO_REFRACTIVE_GEOMETRY, DO_CHAPMAN_FUNCTION,& ! RT model
        DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY, DO_ADDITIONAL_MVOUT,    DO_MVOUT_ONLY,              & ! RT model
        DO_TOA_CONTRIBS,  DO_SPECIALIST_OPTION_1,  DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3,     & ! Specialist
        DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_EXTERNAL_WLEAVE,                & ! Surface
        DO_WATER_LEAVING, DO_FLUORESCENCE,     DO_TF_ITERATION, DO_WLADJUSTED_OUTPUT,                  & ! Surface
        TAYLOR_ORDER, NSTOKES, NSTREAMS, NLAYERS, NFINELAYERS, N_THERMAL_COEFFS, NGREEK_MOMENTS_INPUT, & ! Main numbers
        NLAYERS_NOMS, NLAYERS_CUTOFF, VLIDORT_ACCURACY, FLUX_FACTOR, EARTH_RADIUS,                     & ! Flux/Acc/Radius
        N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES, N_USER_OBSGEOMS, N_USER_LEVELS,                   & ! Geo & lev cont
        SZANGLES,   USER_RELAZMS,   USER_VZANGLES,   USER_OBSGEOMS,   USER_LEVELS,                     & ! Geo & lev cont
        HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, FINEGRID, RFINDEX_PARAMETER, GEOMETRY_SPECHEIGHT,& ! RF
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, FMATRIX_UP, FMATRIX_DN,            & ! Optical
        THERMAL_BB_INPUT, ALBEDO, SURFBB, TOAFLUX, BOAFLUX,               & ! Optical
        ATMOS_WAVELENGTH, TF_MAXITER, TF_CRITERION,                                                    & ! Misc
        DO_WRITE_INPUT,    INPUT_WRITE_FILENAME,    DO_WRITE_FOURIER,       DO_WRITE_RESULTS,          & ! Debug control
        DO_WRITE_SCENARIO, SCENARIO_WRITE_FILENAME, FOURIER_WRITE_FILENAME, RESULTS_WRITE_FILENAME )     ! Debug control

      IF (.NOT. DO_LAMBERTIAN_SURFACE) THEN
!mick mod 9/19/2017 - added NMOMENTS to input
        CALL VLIDORT_WRITE_SUP_BRDF_INPUT ( &
          DO_USER_VZANGLES, DO_SURFACE_EMISSION, &
          NSTOKES,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,NSTREAMS,NMOMS,&
          EXACTDB_BRDFUNC,BRDF_F_0,BRDF_F,USER_BRDF_F_0,USER_BRDF_F,&
          EMISSIVITY,USER_EMISSIVITY)
      END IF

      IF (DO_FOCORR_EXTERNAL) THEN
        CALL VLIDORT_WRITE_SUP_SS_INPUT ( &
          NSTOKES,N_USER_LEVELS,STOKES_SS,STOKES_DB)
      END IF

      IF (DO_SURFACE_LEAVING) THEN
        CALL VLIDORT_WRITE_SUP_SLEAVE_INPUT ( &
          DO_USER_VZANGLES, &
          NSTOKES,N_SZANGLES,N_USER_VZANGLES,N_USER_RELAZMS,NSTREAMS,NMOMS,&
          SLTERM_ISOTROPIC,SLTERM_USERANGLES,SLTERM_F_0,USER_SLTERM_F_0)
      END IF

      END SUBROUTINE VLIDORT_DEBUG_INPUT_MASTER

      END SUBROUTINE VLIDORT_MASTER

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
            NSTOKES, NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_LEVELS, N_THERMAL_COEFFS, NMOMENTS, NSTREAMS_2, & !Input
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

      INTEGER, INTENT (IN) ::          NBEAMS
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

      DOUBLE PRECISION, INTENT (IN) :: BRDF_F        ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F_0      ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F   ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F_0 ( 0:MAXMOMENTS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )
      DOUBLE PRECISION, INTENT (IN) :: EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  New surface-leaving stuff 17 May 2012

      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_F_0 ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   USER_SLTERM_F_0 ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  Disabled

      !DOUBLE PRECISION, INTENT (IN) :: EXACTDB_BRDFUNC    ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: SLTERM_USERANGLES  ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

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
        NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS, N_USER_STREAMS, NBEAMS, NMOMENTS, & ! Input Numbers
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
            NLAYERS, N_PARTLAYERS, NBEAMS, N_USER_STREAMS,  & ! Input numbers
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
      !      DO IBEAM = 1, NBEAMS
      !        DO_REFLECTED_DIRECTBEAM(IBEAM) = .TRUE.
      !      ENDDO
      !    ELSE
      !      DO IBEAM = 1, NBEAMS
      !        DO_REFLECTED_DIRECTBEAM(IBEAM) = .FALSE.
      !      ENDDO
      !    ENDIF
      !  ENDIF
      !ENDIF

      IF ( DO_SOLAR_SOURCES .and. DO_INCLUDE_SURFACE ) THEN
        DO_LOCALBEAM(1:NBEAMS) = DO_REFLECTED_DIRECTBEAM(1:NBEAMS)
      ELSE
        DO_LOCALBEAM(1:NBEAMS) = .FALSE.
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
                FOURIER, NSTOKES, NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS,               & ! Input numbers
                MUELLER_INDEX, FLUX_FACTOR, FLUXVEC, DELTA_FACTOR,                         & ! Input flux-factors, indices
                SZA_LOCAL_INPUT, COS_SZANGLES, TRANS_SOLAR_BEAM, DO_LOCALBEAM,             & ! Input solar beam
                ALBEDO, BRDF_F_0, USER_BRDF_F_0,                                           & ! Input surface reflectance
                SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                             & ! Input surface leaving
                ATMOS_ATTN, RF_DIRECT_BEAM, RF_USER_DIRECT_BEAM, SL_QUADTERM, SL_USERTERM )  ! Output reflected direct beams.
        ENDIF
     ENDIF

!  Useful debug
!      do IBEAM = 1, nbeams
!         if ( FOURIER.eq.0 ) write(*,*)IBEAM,ATMOS_ATTN(IBEAM),RF_DIRECT_BEAM(1,IBEAM,1), RF_USER_DIRECT_BEAM(1,IBEAM,1)
!         if ( FOURIER.eq.0 ) write(*,*)IBEAM,SL_QUADTERM(1,IBEAM,1), SL_USERTERM(1,IBEAM,1)
!      enddo

!  ###################
!  Spherical functions
!  ###################

!  Get Pi Matrices for this Fourier componentIncluding all beam angles !!
!   - Version 2.7. Use "OMP" routine for thread-safety. Mick fix.
!   - Version 2.8. Clean up arguments 7/8/16

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
        CALL VLIDORT_PIMATRIX_SETUP_OMP ( &
          DO_REFRACTIVE_GEOMETRY, DO_USER_STREAMS, FOURIER,             & ! Input flags and Fourier
          NSTOKES, NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS,           & ! Input numbers
          QUAD_STREAMS, USER_STREAMS, COS_SZANGLES, SUN_SZA_COSINES,    & ! Input streams and angles
          NMOMENTS, MUELLER_INDEX, DMAT,                                & ! auxiliary inputs
          PIMM_11, PIMM_KM, PI_XQP, PI_XQM, PI_XUP, PI_XUM, PI_X0P,     & ! Output Pi-Matrix stuff
          PI_XQM_POST, PI_XQM_PRE, PI_XQP_PRE, PI_XUM_POST, PI_XUP_PRE )  ! Output Pi-Matrix stuff
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

      DO IBEAM = 1, NBEAMS

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
                        SL_USERTERM(1:N_USER_STREAMS,IBEAM,O1) = USER_SLTERM_F_0(FOURIER,O1,1:N_USER_STREAMS,IBEAM) * TFACTOR
                     ELSE
                        SL_USERTERM(LUI,IBEAM,O1) = USER_SLTERM_F_0(FOURIER,O1,LUI,IBEAM) * TFACTOR
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

      END MODULE vlidort_masters_m

