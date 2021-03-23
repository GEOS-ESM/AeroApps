
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

      module VLIDORT_Setups_def_m

!  Version 2.8, Inclusion of phase matrices, 03 March 2017

!  Version 2.8.1, September 2019
!      --- New Water-leaving inputs DO_WLADJUSTED_OUTPUT, DO_EXTERNAL_WLEAVE
!      --- New TOA/BOA Isotropic illumination: control and fluxes
!      --- New flags for Planetary problem and mediaproperties

!  Version 2.8.2. 15 April 2020. Summary of Changes
!     1. Removed several Modified Structures in lidort_inputs_def.f90
!     2. This module of Setups is new, uses Derive_Input variables (bookkeeping), and geometries

!  This Module contains the following VLIDORT Structures, with Intents INOUT :

!             VLIDORT_Bookkeep      ==> Bookkeeping stuff
!             VLIDORT_Geometry_FO   ==> FO geometry variables
!             VLIDORT_Geometry_MS   ==> MS geometry variables

      use VLIDORT_PARS_m, only : MAXSTREAMS, MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAXFINELAYERS,   &
                                 MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_STREAMS, MAX_GEOMETRIES, &
                                 MAXBEAMS, MAX_DIRECTIONS, MAX_USER_LEVELS, MAXMOMENTS_INPUT

      implicit none

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Bookkeeping

!  Mode of operation

      LOGICAL   :: DO_MSMODE_VLIDORT, DO_MSMODE_THERMAL, DO_FOCORR_ALONE

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER   :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS, NSTKS_NSTRMS_2 = NSTOKES * NSTREAMS_2
!  total number of layers and streams NTOTAL = NLAYERS * NSTKS_NSTRMS_2
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER   :: NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2
      INTEGER   :: NTOTAL, N_SUBDIAG, N_SUPDIAG

!  Number of directions (1 or 2) and directional array

      INTEGER   :: N_DIRECTIONS
      INTEGER   :: WHICH_DIRECTIONS (MAX_DIRECTIONS)

!  output optical depth masks and indices

      LOGICAL   :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER   :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER   :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER   :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)

!  off-grid optical depths (values, masks, indices)

      LOGICAL   :: DO_PARTLAYERS
      INTEGER   :: N_PARTLAYERS
      INTEGER   :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)
      DOUBLE PRECISION :: PARTLAYERS_VALUES (MAX_PARTLAYERS)

!  Partial heights

      DOUBLE PRECISION :: PARTLAYERS_HEIGHTS (MAX_PARTLAYERS)

!  Number of convergences

      INTEGER   :: N_OUT_STREAMS
      INTEGER   :: N_CONVTESTS

!  Layer masks and limits for doing integrated source terms

      LOGICAL   :: STERM_LAYERMASK_UP (MAXLAYERS)
      LOGICAL   :: STERM_LAYERMASK_DN (MAXLAYERS)

      INTEGER   :: N_ALLLAYERS_UP, N_ALLLAYERS_DN

!  Offsets for geometry indexing
!   4/15/20. Version 2.8.2. Add the Doublet-geometry offset array

      INTEGER   :: LOCAL_N_USERAZM, N_GEOMETRIES
      INTEGER   :: SZA_OFFSETS (MAX_SZANGLES)
      INTEGER   :: VZA_OFFSETS (MAX_SZANGLES, MAX_USER_VZANGLES)
      INTEGER   :: SZD_OFFSETS (MAX_SZANGLES)

!  Post-processing masks. Introduced 4/9/19.

      INTEGER   :: N_PPSTREAMS, PPSTREAM_MASK (MAX_USER_VZANGLES,MAX_SZANGLES)

!  Mueller index

      INTEGER   :: MUELLER_INDEX (MAXSTOKES,MAXSTOKES)

!  Greekmat indices, DMAT, FLUXVEC and related quantities (VLIDORT only)

      INTEGER          :: GREEKMAT_INDEX (6)
      DOUBLE PRECISION :: FLUXVEC (MAXSTOKES), DFLUX (MAXSTOKES)
      DOUBLE PRECISION :: DMAT (MAXSTOKES,MAXSTOKES)

!  Number of particular solutions (not enabled yet)

      INTEGER          :: NPARTICSOLS
      DOUBLE PRECISION :: DMAT_PSOLS (MAXSTOKES,MAXSTOKES,2)
      DOUBLE PRECISION :: DMAT_PSOLS_FLUXVEC (MAXSTOKES,2)

!  Misc control flags

      LOGICAL   :: DO_ALL_FOURIER
      LOGICAL   :: DO_DBCORRECTION

      End TYPE VLIDORT_Bookkeeping

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Geometry_MS

!  Local solar zenith angles Cosines (regular case)

      DOUBLE PRECISION  :: COS_SZANGLES (MAX_SZANGLES)
      DOUBLE PRECISION  :: SIN_SZANGLES (MAX_SZANGLES)

!  Quadrature weights and abscissae, and product

      DOUBLE PRECISION  :: QUAD_STREAMS (MAXSTREAMS)
      DOUBLE PRECISION  :: QUAD_HALFWTS (MAXSTREAMS)
      DOUBLE PRECISION  :: QUAD_WEIGHTS (MAXSTREAMS)
      DOUBLE PRECISION  :: QUAD_STRMWTS (MAXSTREAMS)
      DOUBLE PRECISION  :: QUAD_SINES   (MAXSTREAMS)
      DOUBLE PRECISION  :: QUAD_ANGLES  (MAXSTREAMS)

!  User stream cosines, secants

      DOUBLE PRECISION  :: USER_STREAMS ( MAX_USER_STREAMS )
      DOUBLE PRECISION  :: USER_SECANTS ( MAX_USER_STREAMS )
      DOUBLE PRECISION  :: USER_SINES   ( MAX_USER_STREAMS )
      DOUBLE PRECISION  :: USER_ANGLES  ( MAX_USER_STREAMS )

!  Chapman Factors (Level and (1/9/18) partials, Local solar zenith angles at nadir, average cosines
!    These will be calculated either from FO values, or stand-alone if FO not called.

      DOUBLE PRECISION  :: CHAPMAN_FACTORS  (MAXLAYERS,     MAXLAYERS,MAX_SZANGLES)
      DOUBLE PRECISION  :: PARTIAL_CHAPFACS (MAX_PARTLAYERS,MAXLAYERS,MAX_SZANGLES)

      DOUBLE PRECISION  :: SZA_LEVEL_OUTPUT (0:MAXLAYERS,MAX_SZANGLES)
      DOUBLE PRECISION  :: SUNLAYER_COSINES (MAXLAYERS  ,MAX_SZANGLES)  ! Average cosines for refractive geometry case

!  Azimuth cosines

      END TYPE VLIDORT_Geometry_MS

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Geometry_FO

!  SOLAR Variables
!  ===============

!  Flag for the Nadir case

      LOGICAL           :: doNadir (MAX_GEOMETRIES)

!  Ray constants

      DOUBLE PRECISION  :: Raycon_up (MAX_GEOMETRIES)
      DOUBLE PRECISION  :: Raycon_dn (MAX_GEOMETRIES)

!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      DOUBLE PRECISION  :: Mu0_up (MAX_GEOMETRIES), Mu1_up (MAX_GEOMETRIES)
      DOUBLE PRECISION  :: Mu0_dn (MAX_GEOMETRIES), Mu1_dn (MAX_GEOMETRIES)

!  Cosine scattering angles
   
      DOUBLE PRECISION  :: COSSCAT_UP (MAX_GEOMETRIES)
      DOUBLE PRECISION  :: COSSCAT_DN (MAX_GEOMETRIES)

!  Generalized spherical functions.
!    Rotations(1-4)    = C1, S1, C2, S2

      DOUBLE PRECISION  :: GENSPHER_UP (0:MAXMOMENTS_INPUT,4,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: GENSPHER_DN (0:MAXMOMENTS_INPUT,4,MAX_GEOMETRIES)

      DOUBLE PRECISION  :: ROTATIONS_UP (4,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: ROTATIONS_DN (4,MAX_GEOMETRIES)

!  LOS Quadratures for Enhanced PS. Partials added 8/17/16.

      INTEGER           :: nfinedivs (MAXLAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: xfine     (MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: wfine     (MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)

      INTEGER           :: nfinedivs_p_up (MAX_PARTLAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: xfine_p_up     (MAX_PARTLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: wfine_p_up     (MAX_PARTLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)

      INTEGER           :: nfinedivs_p_dn (MAX_PARTLAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: xfine_p_dn     (MAX_PARTLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: wfine_p_dn     (MAX_PARTLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)

!  LOS path lengths

      DOUBLE PRECISION  :: LosW_paths (MAXLAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: LosP_paths (MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Chapman factor outputs

      DOUBLE PRECISION  :: chapfacs    (MAXLAYERS,     MAXLAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: chapfacs_p  (MAX_PARTLAYERS,MAXLAYERS,MAX_GEOMETRIES)

!  Level-boundary solar paths

      INTEGER           :: ntraverse_up     (0:MAXLAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: sunpaths_up      (0:MAXLAYERS,MAXLAYERS,MAX_GEOMETRIES)
      INTEGER           :: ntraversefine_up (MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: sunpathsfine_up  (MAXLAYERS,MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)

      INTEGER           :: ntraverse_dn     (0:MAXLAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: sunpaths_dn      (0:MAXLAYERS,MAXLAYERS,MAX_GEOMETRIES)
      INTEGER           :: ntraversefine_dn (MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: sunpathsfine_dn  (MAXLAYERS,MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)

!  Partial-level outputs (Sunpaths, fine-sunpaths)

      INTEGER           :: ntraverse_p_up     (MAX_PARTLAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: sunpaths_p_up      (MAX_PARTLAYERS,MAXLAYERS,MAX_GEOMETRIES)
      INTEGER           :: ntraversefine_p_up (MAX_PARTLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: sunpathsfine_p_up  (MAX_PARTLAYERS,MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)

      INTEGER           :: ntraverse_p_dn     (MAX_PARTLAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: sunpaths_p_dn      (MAX_PARTLAYERS,MAXLAYERS,MAX_GEOMETRIES)
      INTEGER           :: ntraversefine_p_dn (MAX_PARTLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)
      DOUBLE PRECISION  :: sunpathsfine_p_dn  (MAX_PARTLAYERS,MAXLAYERS,MAXFINELAYERS,MAX_GEOMETRIES)

!  THERMAL GEOMETRY variables.
!  ===========================

!  Geometry. Los paths added, 8/25/16. Partials added 8/22/16

      DOUBLE PRECISION  :: Mu1_LOS (MAX_USER_STREAMS)
      DOUBLE PRECISION  :: LosW_paths_LOS (MAXLAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION  :: LosP_paths_LOS (MAX_PARTLAYERS,MAX_USER_STREAMS)

!  LOS Quadratures for Enhanced PS. Partials added 8/25/16.
!  5/22/20. Version 2.8.2 Upgrades.
!     ==> Must use vertical distances in Thermal source terms, introduced here (hfine/hfine_p)

      INTEGER           :: nfinedivs_LOS (MAXLAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION  :: xfine_LOS     (MAXLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION  :: wfine_LOS     (MAXLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION  :: hfine_LOS_up     (MAXLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION  :: hfine_LOS_dn     (MAXLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)

      INTEGER           :: nfinedivs_p_LOS_up (MAX_PARTLAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION  :: xfine_p_LOS_up     (MAX_PARTLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION  :: wfine_p_LOS_up     (MAX_PARTLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION  :: hfine_p_LOS_up     (MAX_PARTLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)

      INTEGER           :: nfinedivs_p_LOS_dn (MAX_PARTLAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION  :: xfine_p_LOS_dn     (MAX_PARTLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION  :: wfine_p_LOS_dn     (MAX_PARTLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)
      DOUBLE PRECISION  :: hfine_p_LOS_dn     (MAX_PARTLAYERS,MAXFINELAYERS,MAX_USER_STREAMS)

      END TYPE VLIDORT_Geometry_FO

! #####################################################################
! #####################################################################

!  Three nested sub-structures
      
      TYPE VLIDORT_Setups
         TYPE(VLIDORT_Bookkeeping) ::  Bookkeep
         TYPE(VLIDORT_Geometry_MS) ::  MSGeom
         TYPE(VLIDORT_Geometry_FO) ::  FOGeom
      END TYPE VLIDORT_Setups

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VLIDORT_Bookkeeping, &
                VLIDORT_Geometry_MS, &
                VLIDORT_Geometry_FO, &
                VLIDORT_Setups

end module VLIDORT_Setups_def_m



