
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
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
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
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
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2021.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

      program BRDFplus_Tester

!  This is the Version 2.8 driver. Created 9/18/16 from 2.7 Driver
!     R. Spurr. RT Solutions Inc.

!  Module files for VBRDF and VLIDORT. Strict usage, Version 2.8

      USE VBRDF_SUP_AUX_m, Only : VBRDF_READ_ERROR
      USE VBRDF_SUP_MOD_m
      USE VBRDF_LINSUP_MOD_m

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_LIN_IO_DEFS_m

      USE VLIDORT_VBRDF_SUP_ACCESSORIES_m
      USE VLIDORT_VBRDF_LINSUP_ACCESSORIES_m

      USE VLIDORT_AUX_m,    Only : VLIDORT_READ_ERROR, VLIDORT_WRITE_STATUS

      USE VLIDORT_INPUTS_m, Only : VLIDORT_INPUT_MASTER, VLIDORT_SLEAVE_Sup_Init, VLIDORT_SS_Sup_Init
      USE VLIDORT_MASTERS_m

      USE VLIDORT_L_INPUTS_m, Only : VLIDORT_L_INPUT_MASTER, VLIDORT_SLEAVE_LinSup_Init, VLIDORT_SS_LinSup_Init
      USE VLIDORT_LCS_MASTERS_m

      IMPLICIT NONE

!  NOTES 11 January 2011
!  THIS is a test of the Main Code with VBRDF supplement first.

!  File inputs status structure

      TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus

!  VLIDORT debug input control

      LOGICAL :: DO_DEBUG_INPUT

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn

!  VLIDORT supplements i/o structure

      TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

!  VLIDORT linearized input structures

      TYPE(VLIDORT_Fixed_LinInputs)          :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs)       :: VLIDORT_LinModIn

!  VLIDORT linearized supplements i/o structure

      TYPE(VLIDORT_LinSup_InOut)             :: VLIDORT_LinSup

!  VLIDORT linearized output structure

      TYPE(VLIDORT_LinOutputs)               :: VLIDORT_LinOut

!  VBRDF supplement file inputs status structure

      TYPE(VBRDF_Input_Exception_Handling)   :: VBRDF_Sup_InputStatus

!  VBRDF supplement input structure

      TYPE(VBRDF_Sup_Inputs)                 :: VBRDF_Sup_In

!  VBRDF supplement output structure

      TYPE(VBRDF_Sup_Outputs)                :: VBRDF_Sup_Out
      TYPE(VBRDF_Output_Exception_Handling)  :: VBRDF_Sup_OutputStatus

!  VBRDF supplement linearized input structure

      TYPE(VBRDF_LinSup_Inputs)              :: VBRDF_LinSup_In

!  VBRDF supplement linearized output structure

      TYPE(VBRDF_LinSup_Outputs)             :: VBRDF_LinSup_Out

!  VBRDF supplement / VLIDORT VBRDF-related inputs consistency check status

      TYPE(VLIDORT_Exception_Handling)       :: VLIDORT_VBRDFCheck_Status

!  Local Variables
!  ===============

!  Flag for opening error output file

      LOGICAL ::          OPENFILEFLAG, GCM

!  Help variables
!  mick chg 4/3/2015 - renamed thread --> task

      INTEGER ::          NGREEKMAT_ENTRIES, L, LDUM, LCHECK
      INTEGER ::          K, CK, GK, RK, CMASK(8), GMASK(8), SMASK(8), RMASK(8)

      INTEGER ::          N, N6, UTA, NDUM, V, T, NTASKS, TASK, O1
      INTEGER ::          W, Q, I, J, IB, UM, IA, M, NSTOKESSQ, STATUS_MBCHECK
      DOUBLE PRECISION :: KD, GAER, WAER, TAER, PARCEL, RAYWT, AERWT, EPS
      DOUBLE PRECISION :: AERSCA, AEREXT, MOLSCA, TOTSCA, TOTEXT, EPSFAC
      DOUBLE PRECISION :: MOLOMG(MAXLAYERS), MOLEXT(MAXLAYERS), PAR(6)
      DOUBLE PRECISION :: AMOMS(6), AERMOMS(6,0:MAXMOMENTS_INPUT)
      DOUBLE PRECISION :: RAYMOMS(0:2,MAXLAYERS), DEPOL, BETA2
      DOUBLE PRECISION :: PROBLEM_RAY(6,0:2), SK, CUTOFF, MOM

      INTEGER, PARAMETER :: LUM=1, LIA=1
      CHARACTER(LEN=7)   :: C1
      CHARACTER(LEN=4)   :: C2
      CHARACTER (LEN=120) :: TRACE

      INTEGER ::          OutUnit

!  VBRDF supplement variables

      LOGICAL ::          BS_DO_BRDF_SURFACE
      LOGICAL ::          BS_DO_USER_STREAMS
      LOGICAL ::          BS_DO_SURFACE_EMISSION

      INTEGER ::          BS_NSTOKES

      INTEGER ::          BS_NSTREAMS
      INTEGER ::          BS_NMOMENTS_INPUT

      INTEGER ::          BS_NBEAMS
      DOUBLE PRECISION :: BS_BEAM_SZAS ( MAXBEAMS )

      INTEGER ::          BS_N_USER_RELAZMS
      DOUBLE PRECISION :: BS_USER_RELAZMS (MAX_USER_RELAZMS)

      INTEGER ::          BS_N_USER_STREAMS
      DOUBLE PRECISION :: BS_USER_ANGLES_INPUT (MAX_USER_STREAMS)

      LOGICAL ::          DO_DEBUG_RESTORATION

!  VLIDORT standard input preparation
!   Rob 3/22/15. Introduced FO_CALC local variable, renamed FOCORR 3/1/17

      LOGICAL ::          DO_LAMBERTIAN_SURFACE
      LOGICAL ::          DO_THERMAL_EMISSION, DO_SURFACE_EMISSION
      LOGICAL ::          DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY
      LOGICAL ::          DO_FOCORR, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING
      LOGICAL ::          DO_DELTAM_SCALING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING
      INTEGER ::          NLAYERS, NSTOKES, NFINELAYERS, N_THERMAL_COEFFS

!  Height proxy, # moments

      INTEGER ::          NGREEK_MOMENTS_INPUT
      DOUBLE PRECISION :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Optical proxies. Fmatrix proxies are new for Version 2.8

      DOUBLE PRECISION :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      DOUBLE PRECISION :: DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION :: FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 ) 
      DOUBLE PRECISION :: FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 ) 

      DOUBLE PRECISION :: SURFBB, THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION :: FLUX_FACTOR

!  VLIDORT standard output preparation

      INTEGER ::            N_GEOMETRIES
      INTEGER ::            N_SZANGLES
      INTEGER ::            N_USER_LEVELS
      DOUBLE PRECISION ::   USER_LEVELS ( MAX_USER_LEVELS )
      LOGICAL ::            DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY

!  VLIDORT linearized input preparation

      LOGICAL ::            DO_SIMULATION_ONLY
      LOGICAL ::            DO_COLUMN_LINEARIZATION
      LOGICAL ::            DO_PROFILE_LINEARIZATION
      LOGICAL ::            DO_ATMOS_LINEARIZATION
      LOGICAL ::            DO_SURFACE_LINEARIZATION
      LOGICAL ::            DO_LINEARIZATION
      INTEGER ::            N_TOTALCOLUMN_WFS
      INTEGER ::            N_TOTALPROFILE_WFS
      INTEGER ::            LAYER_VARY_NUMBER ( MAXLAYERS )
      LOGICAL ::            LAYER_VARY_FLAG ( MAXLAYERS )
      DOUBLE PRECISION ::   L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION ::   L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION ::   L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      INTEGER ::            N_SURFACE_WFS

!  Saved baseline results

      DOUBLE PRECISION :: STOKES_BAS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: MEANST_BAS &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FLUX_BAS &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: SURFWF_BAS &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: MEANST_SURFWF_BAS &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FLUX_SURFWF_BAS &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Saved Perturbation results

      INTEGER, PARAMETER :: MAX_TASKS = 10

      DOUBLE PRECISION :: STOKES_PT &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION :: MEANST_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION :: FLUX_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )

!  Define some arrays
!  @@@ Rob Fix 6/13/13. CMASK is for the AEROSOLS, RMASK for RAYLEIGH.

      GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)
      RMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)    !  This   Rayleigh
!      CMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)
      CMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)    !  OUTPUT OF RTS MIE
      SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /)

!  Initialize error file output flag

      OPENFILEFLAG = .false.

!  BASELINE CALCULATION
!  ********************

!  VBRDF SUPPLEMENT CALCULATION
!  ============================

!  Get the VBRDF inputs

      CALL VBRDF_LIN_INPUTMASTER ( &
        'vlidort_v_test/VBRDF_ReadInput.cfg', & ! Input
        VBRDF_Sup_In,         & ! Outputs
        VBRDF_LinSup_In,      & ! Outputs
        VBRDF_Sup_InputStatus ) ! Outputs

      IF ( VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VBRDF_READ_ERROR ( 'V2p8p3_VBRDF_ReadInput.log', VBRDF_Sup_InputStatus )

!  Baseline parameters

      par(1) = VBRDF_Sup_In%BS_BRDF_FACTORS(1)
      par(2) = VBRDF_Sup_In%BS_BRDF_FACTORS(2)
      par(3) = VBRDF_Sup_In%BS_BRDF_PARAMETERS(2,1)
      par(4) = VBRDF_Sup_In%BS_BRDF_PARAMETERS(2,2)
      par(5) = VBRDF_Sup_In%BS_BRDF_FACTORS(3)
      par(6) = VBRDF_Sup_In%BS_BRDF_PARAMETERS(3,1)

!  Do not want debug restoration

      DO_DEBUG_RESTORATION = .false.

!  Set output file designation for later use

      do_observation_geometry = VBRDF_Sup_In%BS_DO_USER_OBSGEOMS
      do_doublet_geometry     = VBRDF_Sup_In%BS_DO_DOUBLET_GEOMETRY

      if (do_observation_geometry) then
        c1 = 'obsgeo'
      else if (do_doublet_geometry) then
        c1 = 'doublet'
      else
        c1 = 'lattice'
      endif

!  A normal calculation will require

      BS_NMOMENTS_INPUT = 2 * VBRDF_Sup_In%BS_NSTREAMS - 1

!  Linearized VBRDF call
!    The output will now be used for the Main VLIDORT calculation

      CALL VBRDF_LIN_MAINMASTER ( &
        DO_DEBUG_RESTORATION,     & ! Inputs
        BS_NMOMENTS_INPUT,        & ! Inputs
        VBRDF_Sup_In,             & ! Inputs
        VBRDF_LinSup_In,          & ! Inputs
        VBRDF_Sup_Out,            & ! Outputs
        VBRDF_LinSup_Out,         & ! Outputs
        VBRDF_Sup_OutputStatus )    ! Outputs

      if ( VBRDF_Sup_OutputStatus%BS_STATUS_OUTPUT .ne. vlidort_success ) THEN
         TRACE = 'VBRDF_LIN_MAINMASTER failed' ; go to 678
      endif

!  MAIN Calculation
!  ================

!  Read Main VLIDORT inputs
!  ------------------------

!  (EITHER USE THE SAVED VBRDF INPUTS and COPY
!    OR READ NEW INPUTS and CHECK SAME)

!  VLIDORT control Read input, abort if failed

      CALL VLIDORT_L_INPUT_MASTER ( &
        'vlidort_v_test/V2p8p3_VLIDORT_ReadInput.cfg', & ! Input
        VLIDORT_FixIn,      & ! Outputs
        VLIDORT_ModIn,      & ! Outputs
        VLIDORT_LinFixIn,   & ! Outputs
        VLIDORT_LinModIn,   & ! Outputs
        VLIDORT_InputStatus ) ! Outputs

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( 'V2p8p3_VLIDORT_ReadInput.log', VLIDORT_InputStatus )

!  1/31/21. Version 2.8.3, New variables must be set by hand

      VLIDORT_FixIn%Cont%TS_ASYMTX_TOLERANCE = 1.0d-20

!  1/31/21. Version 2.8.3, Set the DO_MSSTS flag to generate output for the MS sphericity corrections

      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .false.

!  There are very strict conditions for this specialist option, as follows :==>
!     1. Either DO_UPWELLING or DO_DNWELLING must be set, Not Both !!!!
!     2. DO_FULLRAD_MODE must be set
!     3. DO_OBSERVATION_GEOMETRY must be set, with N_USER_OBSGEOMS = 3
!     4. DO_FOCORR and DO_FOCORR_OUTGOING must both be set
!     5a. Upwelling  : N_USER_LEVELS = 1, and USER_LEVELS(1) = 0.0            [ TOA output only ]
!     5b. Downwelling: N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(nlayers)  [ BOA output only ]
!  These checks have been implemented inside VLIDORT.

!  1/31/21. Version 2.8.3, Set the NSTOKES2 FOURIER0 flag by hand

      VLIDORT_FixIn%Bool%TS_DO_FOURIER0_NSTOKES2 = .true.

!  mick add 3/31/2015 - add wavelength diagnostic (just a dummy here)
!   Version 2.8, now a configuration-file read. Overrides the value in there!!!
!   Needs to be the same as the VBRDF wavelength

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = 0.333d0

!  Define some input variables not handled by VLIDORT_L_INPUT_MASTER

      CALL VLIDORT_SLEAVE_Sup_Init ( VLIDORT_Sup )
      CALL VLIDORT_SS_Sup_Init ( VLIDORT_Sup )
      CALL VLIDORT_SLEAVE_LinSup_Init ( VLIDORT_LinSup )
      CALL VLIDORT_SS_LinSup_Init ( VLIDORT_LinSup )

!  Set some local variables

      nstokes = VLIDORT_FixIn%Cont%TS_nstokes

!   IMPORTANT - Check compatibility of VBRDF and Main VLIDORT inputs
!   ----------------------------------------------------------------

!  VBRDF surface

      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = .false.

!  Set the number of weighting functions

      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS = VBRDF_LinSup_In%BS_N_SURFACE_WFS

!  This is the checking routine

      CALL VLIDORT_VBRDF_INPUT_CHECK ( &
        VBRDF_Sup_In,             & ! Inputs
        VLIDORT_FixIn,            & ! Inputs
        VLIDORT_ModIn,            & ! Inputs
        VLIDORT_VBRDFCheck_Status ) ! Outputs

      IF ( VLIDORT_VBRDFCheck_Status%TS_STATUS_INPUTCHECK .ne. VLIDORT_SUCCESS ) THEN
        write(*,'(/1X,A)') 'VLIDORT/VBRDF baseline check:'
        CALL VLIDORT_VBRDF_INPUT_CHECK_ERROR ( &
          'V2p8p3_VLIDORT_VBRDFcheck.log', VLIDORT_VBRDFCheck_Status )
      ENDIF

!  Alternative: Force inputs to agree
!  ----------------------------------

!      NSTOKES  = BS_NSTOKES
!      NSTREAMS = BS_NSTREAMS

!      DO_MVOUT_ONLY         = .not. BS_DO_USER_STREAMS
!      DO_LAMBERTIAN_SURFACE = .not. BS_DO_BRDF_SURFACE
!      DO_SURFACE_EMISSION   = BS_DO_SURFACE_EMISSION
!      DO_USER_VZANGLES  = BS_DO_USER_STREAMS

!      N_SZANGLES = BS_NBEAMS
!      DO I = 1, N_SZANGLES
!        SZANGLES(I) = BS_BEAM_SZAS(I)
!      ENDDO

!      N_USER_VZANGLES = BS_N_USER_STREAMS
!      DO I = 1, N_USER_VZANGLES
!        USER_VZANGLES(I) = BS_USER_ANGLES_INPUT(I)
!      ENDDO

!      N_USER_RELAZMS = BS_N_USER_RELAZMS
!      DO I = 1,N_USER_RELAZMS
!        USER_RELAZMS(I) = BS_USER_RELAZMS(I)
!      ENDDO

!  Copy VBRDF Sup outputs to VLIDORT's VBRDF Sup inputs (std & lin)

      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .TRUE.
      CALL SET_VLIDORT_L_VBRDF_INPUTS ( &
        VBRDF_Sup_Out, VBRDF_LinSup_Out,      & !Inputs
        VLIDORT_FixIn, VLIDORT_ModIn,         & !Inputs
        VLIDORT_LinFixIn, VLIDORT_LinModIn,   & !Inputs
        VLIDORT_Sup, VLIDORT_LinSup)            !Outputs

!  mick note 9/19/2017 - adjusted indices on some of the BRDF output arrays to promote
!                        consistency in their use
!                      - added this table for clarity
!
!                        Index     Task for BRDF output
!                        ----- ----------------------------
!                         ib   incident solar zenith angle
!                         i    incident quadrature angle
!                         j    reflected quadrature angle
!                         um   reflected user zenith angle
!                         ia   reflected user azimuth angle

!  Write normal VBRDF supplement results (as a check)

      GCM = (VBRDF_Sup_In%BS_which_brdf(1) .eq. 10 .or. &
             VBRDF_Sup_In%BS_which_brdf(2) .eq. 10 .or. &
             VBRDF_Sup_In%BS_which_brdf(3) .eq. 10)

      if ( GCM ) then
         NSTOKESSQ = VBRDF_Sup_In%BS_NSTOKES * VBRDF_Sup_In%BS_NSTOKES
      else
         NSTOKESSQ = 1
      endif

      open(36, file='vlidort_v_test/results_brdf_supcheck_' // trim(c1) // '.res', &
           status = 'unknown')
      write(36,'(/a/)') ' Number of Beams/Quad-streams/User-streams/' // &
                        'User-relazms/ObsGeoms/nstokessq/GCM flag'
      write(36,'(6i5,2x,L2)') &
        VBRDF_Sup_In%BS_nbeams,         VBRDF_Sup_In%BS_nstreams, &
        VBRDF_Sup_In%BS_n_user_streams, VBRDF_Sup_In%BS_n_user_relazms, &
        VBRDF_Sup_In%BS_n_user_obsgeoms, &
        NSTOKESSQ, GCM

      write(36,'(/T20,a/)') 'Exact BRDF, all geometries'
      if (do_observation_geometry) then
        do ib = 1, VBRDF_Sup_In%BS_nbeams
          write(36,'(3i4,1p16e20.11)') ib, lum, lia, &
            (VLIDORT_Sup%BRDF%TS_exactdb_brdfunc(q,lum,lia,ib), &
             q = 1, nstokessq)
        enddo
      elseif (do_doublet_geometry) then
        do ib = 1, VBRDF_Sup_In%BS_nbeams
          do um = 1, VBRDF_Sup_In%BS_n_user_doublets
             write(36,'(3i4,1p16e20.11)') ib, um, lia, &
               (VLIDORT_Sup%BRDF%TS_exactdb_brdfunc(q,um,lia,ib), &
                q = 1, nstokessq)
          enddo
        enddo
      else
        do ib = 1, VBRDF_Sup_In%BS_nbeams
          do um = 1, VBRDF_Sup_In%BS_n_user_streams
            do ia = 1, VBRDF_Sup_In%BS_n_user_relazms
             write(36,'(3i4,1p16e20.11)') ib, um, ia, &
               (VLIDORT_Sup%BRDF%TS_exactdb_brdfunc(q,um,ia,ib), &
                q = 1, nstokessq)
            enddo
          enddo
        enddo
      endif

      DO M = 0, BS_NMOMENTS_INPUT
!mick mod 9/19/2017 - swapped function of i and j indices in both write statements
        write(36,'(/T20,a,i3/)') 'Fourier component # ', M
        do j = 1, VBRDF_Sup_In%BS_nstreams
          write(36,99) J, ((VLIDORT_Sup%BRDF%TS_BRDF_F(M,Q,J,I), &
            I = 1, VBRDF_Sup_In%BS_nstreams), q = 1, nstokessq)
          write(36,99) J, ((VLIDORT_Sup%BRDF%TS_BRDF_F_0(M,Q,J,IB), &
            IB = 1, VBRDF_Sup_In%BS_nbeams), q = 1, nstokessq)
        enddo

!mick mod 9/19/2017 - swapped index j for i in both USER_BRDF_F write statements
        write(36,*)' '
        if (do_observation_geometry) then
          do um = 1, VBRDF_Sup_In%BS_n_user_streams
            write(36,99) UM, ((VLIDORT_Sup%BRDF%TS_USER_BRDF_F(M,Q,um,I), &
              I = 1, VBRDF_Sup_In%BS_nstreams), q = 1, nstokessq)
          enddo
          write(36,99) LUM, ((VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0(M,Q,lum,IB), &
            IB = 1, VBRDF_Sup_In%BS_nbeams), q = 1, nstokessq)
        elseif (do_doublet_geometry) then
          do um = 1, VBRDF_Sup_In%BS_n_user_streams
            write(36,99) UM, ((VLIDORT_Sup%BRDF%TS_USER_BRDF_F(M,Q,um,I), &
              I = 1, VBRDF_Sup_In%BS_nstreams), q = 1, nstokessq)
            write(36,99) UM, ((VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0(M,Q,um,IB), &
              IB = 1, VBRDF_Sup_In%BS_nbeams), q = 1, nstokessq)
          enddo
        else
          do um = 1, VBRDF_Sup_In%BS_n_user_streams
            write(36,99) UM, ((VLIDORT_Sup%BRDF%TS_USER_BRDF_F(M,Q,um,I), &
              I = 1, VBRDF_Sup_In%BS_nstreams), q = 1, nstokessq)
            write(36,99) UM, ((VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0(M,Q,um,IB), &
              IB = 1, VBRDF_Sup_In%BS_nbeams), q = 1, nstokessq)
          enddo
        endif

      ENDDO
      close(36)

98    format(3I4,1p16e15.7)
99    format(I4,1p500e15.7)

!  Write linearized VBRDF supplement results

      if ( .not. do_debug_restoration ) then
      open(36, file = 'vlidort_v_test/results_brdf_supcheck_' // trim(c1) // '.wfs', &
           status = 'unknown')
      write(36,'(/a/)') ' Number of SurfaceWFs/Beams/Quad-streams/' // &
                        'User-streams/User-relazms/ObsGeoms/nstokessq/GCM flag'
      write(36,'(7i5,2x,L2)') &
        VBRDF_LinSup_In%BS_N_SURFACE_WFS,  VBRDF_Sup_In%BS_nbeams, &
        VBRDF_Sup_In%BS_nstreams,          VBRDF_Sup_In%BS_n_user_streams, &
        VBRDF_Sup_In%BS_n_user_relazms,    VBRDF_Sup_In%BS_n_user_obsgeoms, &
        NSTOKESSQ, GCM

      do w = 1, VBRDF_LinSup_In%BS_N_SURFACE_WFS
        write(36,'(/T20,a,i2/)') 'Weighting function # ', W
        write(36,'(/T20,a/)') 'Exact BRDF, all geometries'
        if (do_observation_geometry) then
          do ib = 1, VBRDF_Sup_In%BS_nbeams
            write(36,98) ib, lum, lia, &
              (VLIDORT_LinSup%BRDF%TS_LS_exactdb_brdfunc(w,q,lum,lia,ib), &
               q = 1, nstokessq)
          enddo
        elseif (do_doublet_geometry) then
          do ib = 1, VBRDF_Sup_In%BS_nbeams
           do um = 1, VBRDF_Sup_In%BS_n_user_streams
             write(36,98) ib, um, lia, &
               (VLIDORT_LinSup%BRDF%TS_LS_exactdb_brdfunc(w,q,um,lia,ib), &
                q = 1, nstokessq)
           enddo
          enddo
        else
          do ib = 1, VBRDF_Sup_In%BS_nbeams
           do um = 1, VBRDF_Sup_In%BS_n_user_streams
            do ia = 1, VBRDF_Sup_In%BS_n_user_relazms
             write(36,98) ib, um, ia, &
               (VLIDORT_LinSup%BRDF%TS_LS_exactdb_brdfunc(w,q,um,ia,ib), &
                q = 1, nstokessq)
            enddo
           enddo
          enddo
        endif

        do M = 0, 2*VBRDF_Sup_In%BS_nstreams - 1
!mick mod 9/19/2017 - swapped function of i and j indices in both write statements
         write(36,'(/T20,a,i3/)') 'Fourier component # ', M
         do j = 1, VBRDF_Sup_In%BS_nstreams
          write(36,99) J, ((VLIDORT_LinSup%BRDF%TS_LS_BRDF_F(w,M,Q,J,I), &
            I = 1, VBRDF_Sup_In%BS_nstreams), q = 1, nstokessq)
          write(36,99) J, ((VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0(w,M,Q,J,IB), &
            IB = 1, VBRDF_Sup_In%BS_nbeams),  q = 1, nstokessq)
         enddo

!mick mod 9/19/2017 - swapped index j for i in both LS_USER_BRDF_F write statements
         write(36,*) ' '
         if (do_observation_geometry) then
           do um = 1, VBRDF_Sup_In%BS_n_user_streams
            write(36,99) UM, &
              ((VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F(w,M,Q,um,I), &
                I = 1, VBRDF_Sup_In%BS_nstreams), q = 1, nstokessq)
           enddo
           write(36,99) LUM, &
             ((VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0(w,M,Q,lum,IB), &
               IB = 1, VBRDF_Sup_In%BS_nbeams), q = 1, nstokessq)
         elseif (do_doublet_geometry) then
           do um = 1, VBRDF_Sup_In%BS_n_user_streams
            write(36,99) UM, &
              ((VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F(w,M,Q,um,I), &
                I = 1, VBRDF_Sup_In%BS_nstreams), q = 1, nstokessq)
            write(36,99) UM, &
              ((VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0(w,M,Q,um,IB), &
                IB = 1, VBRDF_Sup_In%BS_nbeams), q = 1, nstokessq)
           enddo
         else
           do um = 1, VBRDF_Sup_In%BS_n_user_streams
            write(36,99) UM, &
              ((VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F(w,M,Q,um,I), &
                I = 1, VBRDF_Sup_In%BS_nstreams), q = 1, nstokessq)
            write(36,99) UM, &
              ((VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0(w,M,Q,um,IB), &
                IB = 1, VBRDF_Sup_In%BS_nbeams), q = 1, nstokessq)
           enddo
         endif
        enddo
      enddo
      close(36)
      endif

!  Get the pre-prepared atmosphere

      height_grid = 0.0d0
      nlayers = VLIDORT_FixIn%Cont%TS_nlayers
      open(45,file='vlidort_v_test/input_atmos.dat',status='old' )
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(1,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0) = 60.0d0
      do n = 1, nlayers
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n), molext(n),molomg(n),kd,kd,kd,kd
      enddo
      close(45)

!  Set masking limits

      if ( nstokes .eq. 1 ) ngreekmat_entries = 1
      if ( nstokes .eq. 3 ) ngreekmat_entries = 5
      if ( nstokes .eq. 4 ) ngreekmat_entries = 8

!  Add Aerosols bottom 6 layers, spread evenly
!     Gamma-distributed, 0.782 Microns, n = 1.43, reff = 1.05, veff = 0.07

      open(1,file='vlidort_v_test/ProblemIII.Moms',status = 'unknown')
      do L = 1, 5
        read(1,*)
      enddo
      LCHECK = 0
      L      = -1
      cutoff = 1.0d-06
      DO while (L.lt.106 .and. LCHECK.lt.2)
         read(1,*)LDUM,(amoms(k),k=1,6)
         if ( amoms(1).gt.cutoff ) then
           L = L + 1
           do k = 1, 6
             aermoms(k,L) = amoms(k)
           enddo
           if ( Lcheck.eq.1 ) Lcheck = Lcheck + 1
         else
           Lcheck = Lcheck + 1
           if ( Lcheck .eq. 1 ) then
             L = L + 1
             do k = 1, 6
               aermoms(k,L) = amoms(k)
             enddo
           endif
         endif
      enddo
      close(1)
      ngreek_moments_input = L

!  Rayleigh scattering law

!  @@@ Rob Fix 6/13/13. Formula for depol was wrong, Index was wrong.
!                       However, this only affects the circular-polarization
!                       Also the (3,1) instead of (4,1) will not affect result

!  Old code
!      depol   = ( 1.0d0 - 2.0d0*raymoms(2,NLAYERS) ) / &
!                ( 1.0d0 + 2.0d0*raymoms(2,NLAYERS) )
      depol   = ( 1.0d0 - 2.0d0*raymoms(2,NLAYERS) ) / &
                ( 1.0d0 + raymoms(2,NLAYERS) )

      beta2 = raymoms(2,NLAYERS)
      PROBLEM_RAY = 0.0d0
      PROBLEM_RAY(1,0) =  1.D0

!  Old code, index was wrong
!      PROBLEM_RAY(3,1) =  3.D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
      PROBLEM_RAY(4,1) =  3.D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
      PROBLEM_RAY(1,2) =  beta2
      PROBLEM_RAY(5,2) =  - DSQRT(6.D0) * beta2
      PROBLEM_RAY(2,2) =  6.D0 * beta2

!  Initialise optical proxies

      deltau_vert_input    = 0.0d0
      omega_total_input    = 0.0d0
      greekmat_total_input = 0.0d0
      Fmatrix_up           = zero ! 2p8
      Fmatrix_dn           = zero ! 2p8

!  Fix the non-aerosol layers

      n6 = nlayers - 6
      do n = 1, n6
        deltau_vert_input(n) = molext(n)
        omega_total_input(n) = molomg(n)
        DO K = 1, ngreekmat_entries
          GK = gmask(k) ; rk = rmask(k) ! ; ck = cmask(k)
          DO L = 0, 2
! Old       GREEKMAT_TOTAL_INPUT(L,n,gk) = PROBLEM_RAY(ck,L)
            GREEKMAT_TOTAL_INPUT(L,n,gk) = PROBLEM_RAY(rk,L)
          ENDDO
          GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0
       ENDDO
      enddo

!  Fix layers with aerosol

      waer = 0.99999d0 ; taer = 0.5d0
      parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
      do n = n6 + 1, nlayers
        aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
        aersca = aerext * waer
        molsca = molomg(n) * molext(n)
        totext = molext(n) + aerext
        totsca = molsca    + aersca
        raywt  = molsca / totsca
        aerwt  = aersca / totsca
        deltau_vert_input(n) = totext
        omega_total_input(n) = totsca / totext
        DO K = 1, ngreekmat_entries
          GK = gmask(k); sk = dble(smask(k)) ; ck = cmask(k) ; rk = rmask(k)
          DO L = 0, 2
! Old       MOM = RAYWT*PROBLEM_RAY(ck,L) + AERWT*SK*AERMOMS(CK,L)
            MOM = RAYWT*PROBLEM_RAY(rk,L) + AERWT*SK*AERMOMS(CK,L)
            GREEKMAT_TOTAL_INPUT(L,N,GK) = MOM
          ENDDO
          DO L = 3, NGREEK_MOMENTS_INPUT
            MOM = AERWT * SK * AERMOMS(CK,L)
            GREEKMAT_TOTAL_INPUT(L,n,gk) = MOM
          ENDDO
        ENDDO
        GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0
      enddo

!  Linearization input

      L_deltau_vert_input = 0.0d0
      L_omega_total_input = 0.0d0
      L_greekmat_total_input = 0.0d0

      DO_SIMULATION_ONLY       = .FALSE.
      DO_COLUMN_LINEARIZATION  = .FALSE.
      DO_PROFILE_LINEARIZATION = .FALSE.
      DO_ATMOS_LINEARIZATION   = .FALSE.
      DO_SURFACE_LINEARIZATION = .TRUE.
      DO_LINEARIZATION         = .TRUE.

!  2p6 variables superceded
!      DO_LTE_LINEARIZATION     = .FALSE.
!      DO_SURFBB_LINEARIZATION  = .FALSE.

      n_totalcolumn_wfs        = 0
      n_totalprofile_wfs       = 0

      layer_vary_number = 0
      layer_vary_flag   = .false.

!  Baseline calculation,
!      In/Outgoing  Single-scatter correction, With delta-M scaling. Upgrade 2.8, 3/1/17

      DO_FOCORR          = .true.
      DO_FOCORR_NADIR    = .true.
      DO_FOCORR_OUTGOING = .false.
      !DO_FOCORR_NADIR    = .false.
      !DO_FOCORR_OUTGOING = .true.

      DO_DELTAM_SCALING  = .true.
      DO_SOLUTION_SAVING = .false.
      DO_BVP_TELESCOPING = .false.
      !DO_SOLUTION_SAVING = .true.
      !DO_BVP_TELESCOPING = .true.
      NFINELAYERS        = 4

!  Copy local inputs to type-structure inputs:

!  Misc

      VLIDORT_FixIn%Cont%TS_nlayers               = nlayers
      VLIDORT_ModIn%MCont%TS_ngreek_moments_input = ngreek_moments_input
      VLIDORT_FixIn%Chapman%TS_height_grid        = height_grid

!  Linearization control

      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = DO_SIMULATION_ONLY
      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = DO_ATMOS_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = DO_SURFACE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = DO_LINEARIZATION

      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS         = N_TOTALCOLUMN_WFS
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = N_TOTALPROFILE_WFS

      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER         = LAYER_VARY_NUMBER
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG           = LAYER_VARY_FLAG

      VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES  = ' '
      VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES = ' '

!  Normal and linearized optical properties

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input
      VLIDORT_FixIn%Optical%TS_FMATRIX_UP           = fmatrix_up ! zero here
      VLIDORT_FixIn%Optical%TS_FMATRIX_DN           = fmatrix_dn ! zero here

      VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input    = L_deltau_vert_input
      VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input = L_greekmat_total_input
      VLIDORT_LinFixIn%Optical%TS_L_omega_total_input    = L_omega_total_input

!  Additional Control. FO upgrade for V2.8, 3/1/17

      VLIDORT_ModIn%MBool%TS_DO_FOCORR          = DO_FOCORR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR    = DO_FOCORR_NADIR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING = DO_FOCORR_OUTGOING

      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING  = DO_DELTAM_SCALING
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = DO_SOLUTION_SAVING
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = DO_BVP_TELESCOPING

      VLIDORT_FixIn%Cont%TS_NFINELAYERS         = NFINELAYERS

!  Progress

      write(*,*)'Doing Baseline calculation'

!  VLIDORT call

      do_debug_input = .false.
      CALL VLIDORT_LCS_MASTER ( do_debug_input, &
        VLIDORT_FixIn, &
        VLIDORT_ModIn, &
        VLIDORT_Sup, &
        VLIDORT_Out, &
        VLIDORT_LinFixIn, &
        VLIDORT_LinModIn, &
        VLIDORT_LinSup, &
        VLIDORT_LinOut )

!  Exception handling, write-up (optional)
!    Will generate file only if errors or warnings are encountered

      CALL VLIDORT_WRITE_STATUS ( &
        'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Set some variables

      n_szangles    = VLIDORT_ModIn%MSunrays%TS_n_szangles
      n_user_levels = VLIDORT_FixIn%UserVal%TS_n_user_levels
      user_levels   = VLIDORT_ModIn%MUserVal%TS_user_levels
      n_geometries  = VLIDORT_Out%Main%TS_n_geometries
      n_surface_wfs = VLIDORT_LinFixIn%Cont%TS_n_surface_wfs

!  Save Baseline results

      do uta = 1, n_user_levels
        do o1 = 1, nstokes

          do v = 1, n_geometries
            stokes_BAS(uta,v,o1,upidx) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,upidx)
            stokes_BAS(uta,v,o1,dnidx) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,dnidx)
            do q = 1, n_surface_wfs
              SURFWF_BAS(q,uta,v,o1,upidx)= VLIDORT_LinOut%Surf%TS_surfacewf(q,uta,v,o1,upidx)
              SURFWF_BAS(q,uta,v,o1,dnidx)= VLIDORT_LinOut%Surf%TS_surfacewf(q,uta,v,o1,dnidx)
            enddo
          enddo

          do v = 1, n_szangles
            MEANST_BAS(uta,v,o1,upidx) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,upidx)
            MEANST_BAS(uta,v,o1,dnidx) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,dnidx) &
                                       + VLIDORT_Out%Main%TS_DNMEANST_DIRECT(uta,v,o1)
            FLUX_BAS(uta,v,o1,upidx)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,upidx)
            FLUX_BAS(uta,v,o1,dnidx)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,dnidx) &
                                       + VLIDORT_Out%Main%TS_DNFLUX_DIRECT(uta,v,o1)

            do q = 1, n_surface_wfs
              MEANST_SURFWF_BAS(q,uta,v,o1,upidx) = VLIDORT_LinOut%Surf%TS_MEANST_DIFFUSE_SURFWF(q,uta,v,o1,upidx)
              MEANST_SURFWF_BAS(q,uta,v,o1,dnidx) = VLIDORT_LinOut%Surf%TS_MEANST_DIFFUSE_SURFWF(q,uta,v,o1,dnidx)
              FLUX_SURFWF_BAS(q,uta,v,o1,upidx)   = VLIDORT_LinOut%Surf%TS_FLUX_DIFFUSE_SURFWF(q,uta,v,o1,upidx)
              FLUX_SURFWF_BAS(q,uta,v,o1,dnidx)   = VLIDORT_LinOut%Surf%TS_FLUX_DIFFUSE_SURFWF(q,uta,v,o1,dnidx)
            enddo
          enddo

        enddo
      enddo

!  FINITE DIFFERENCE SECTION
!  *************************

!  Reset VLIDORT's linearized surface inputs

      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .FALSE.
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 0

!  Define perturbation for FD computations

      eps = 1.0d-04
      epsfac = 1.0d0 + eps

!  Number of tasks

      ntasks = 6

!  For each task

      do t = 1, ntasks

!  progress

        write(*,*)'Doing task # ',t

!  Read the VBRDF information again

        CALL VBRDF_INPUTMASTER ( &
          'vlidort_v_test/VBRDF_ReadInput.cfg', &
          VBRDF_Sup_In, &
          VBRDF_Sup_InputStatus )

      IF ( VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VBRDF_READ_ERROR ( 'V2p8p3_VBRDF_ReadInput.log', VBRDF_Sup_InputStatus )

!  Perturbation

        if ( t .eq. 1 ) VBRDF_Sup_In%BS_brdf_factors(1) = &
                        VBRDF_Sup_In%BS_brdf_factors(1) * epsfac
        if ( t .eq. 2 ) VBRDF_Sup_In%BS_brdf_factors(2) = &
                        VBRDF_Sup_In%BS_brdf_factors(2) * epsfac
        if ( t .eq. 3 ) VBRDF_Sup_In%BS_BRDF_PARAMETERS(2,1) = &
                        VBRDF_Sup_In%BS_BRDF_PARAMETERS(2,1) * epsfac
        if ( t .eq. 4 ) VBRDF_Sup_In%BS_BRDF_PARAMETERS(2,2) = &
                        VBRDF_Sup_In%BS_BRDF_PARAMETERS(2,2) * epsfac
        if ( t .eq. 5 ) VBRDF_Sup_In%BS_brdf_factors(3) = &
                        VBRDF_Sup_In%BS_brdf_factors(3) * epsfac
        if ( t .eq. 6 ) VBRDF_Sup_In%BS_BRDF_PARAMETERS(3,1) = &
                        VBRDF_Sup_In%BS_BRDF_PARAMETERS(3,1) * epsfac

!  A normal calculation will require

        BS_NMOMENTS_INPUT = 2 * VBRDF_Sup_In%BS_NSTREAMS - 1

!  Standard VBRDF call

        CALL VBRDF_MAINMASTER ( &
          DO_DEBUG_RESTORATION,    & ! Inputs
          BS_NMOMENTS_INPUT,       & ! Inputs
          VBRDF_Sup_In,            & ! Inputs
          VBRDF_Sup_Out,           & ! Outputs
          VBRDF_Sup_OutputStatus )   ! Output exception handling

        if ( VBRDF_Sup_OutputStatus%BS_STATUS_OUTPUT .ne. vlidort_success ) THEN
           write(c1,'(i1)')t ; TRACE = 'VBRDF_MAINMASTER failed, perturbation # '//trim(c1) ; go to 678
        endif

!  GCM flag

        GCM = (VBRDF_Sup_In%BS_which_brdf(1) .eq. 10 .or. &
               VBRDF_Sup_In%BS_which_brdf(2) .eq. 10 .or. &
               VBRDF_Sup_In%BS_which_brdf(3) .eq. 10)

        if ( GCM ) then
           NSTOKESSQ = VBRDF_Sup_In%BS_NSTOKES * VBRDF_Sup_In%BS_NSTOKES
        else
           NSTOKESSQ = 1
        endif

!   IMPORTANT - Check compatibility of VBRDF and MAIN Inputs
!   --------------------------------------------------------

        CALL VLIDORT_VBRDF_INPUT_CHECK ( &
          VBRDF_Sup_In,             & ! Inputs
          VLIDORT_FixIn,            & ! Inputs
          VLIDORT_ModIn,            & ! Inputs
          VLIDORT_VBRDFCheck_Status ) ! Outputs

        IF ( VLIDORT_VBRDFCheck_Status%TS_STATUS_INPUTCHECK .ne. vlidort_success ) THEN
          write(*,'(/1X,A)') 'VLIDORT/VBRDF FD check:'
          CALL VLIDORT_VBRDF_INPUT_CHECK_ERROR ( &
             'V2p8p3_VLIDORT_VBRDFcheck.log', VLIDORT_VBRDFCheck_Status )
        ENDIF

!  Copy VBRDF Sup outputs to VLIDORT's VBRDF Sup inputs (std only)

        CALL SET_VLIDORT_VBRDF_INPUTS ( &
          VBRDF_Sup_Out, VLIDORT_FixIn, VLIDORT_ModIn, & !Inputs
          VLIDORT_Sup )                                  !Outputs

!  VLIDORT call

        do_debug_input = .false.
        CALL VLIDORT_MASTER ( do_debug_input, &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

!  Exception handling, write-up (optional)
!    Will generate file only if errors or warnings are encountered

        CALL VLIDORT_WRITE_STATUS ( &
          'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Copy some variables

        n_szangles    = VLIDORT_ModIn%MSunrays%TS_n_szangles
        n_user_levels = VLIDORT_FixIn%UserVal%TS_n_user_levels
        user_levels   = VLIDORT_ModIn%MUserVal%TS_user_levels
        n_geometries  = VLIDORT_Out%Main%TS_n_geometries

!  Save results

        do uta = 1, n_user_levels
          do o1 = 1, nstokes
            do v = 1, n_geometries
              stokes_PT(uta,v,o1,upidx,t) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,upidx)
              stokes_PT(uta,v,o1,dnidx,t) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,dnidx)
            enddo
            do v = 1, n_szangles
              MEANST_PT(uta,v,o1,upidx,t) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,upidx)
              MEANST_PT(uta,v,o1,dnidx,t) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,dnidx) &
                                          + VLIDORT_Out%Main%TS_DNMEANST_DIRECT(uta,v,o1)
              FLUX_PT(uta,v,o1,upidx,t)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,upidx)
              FLUX_PT(uta,v,o1,dnidx,t)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,dnidx) &
                                          + VLIDORT_Out%Main%TS_DNFLUX_DIRECT(uta,v,o1)
            enddo
          enddo
        enddo

!  End task loop

      enddo

!  Write file of Surface WF test results
!  =====================================

!  Title of output file

      do_observation_geometry = VBRDF_Sup_In%BS_DO_USER_OBSGEOMS
      do_doublet_geometry     = VBRDF_Sup_In%BS_DO_DOUBLET_GEOMETRY

      if (do_observation_geometry) then
        c1 = 'obsgeo'
      else if (do_doublet_geometry) then
        c1 = 'doublet'
      else
        c1 = 'lattice'
      endif

      if ( nstokes .eq. 1 ) then
        c2 = 'I000'
      else if ( nstokes .eq. 3 ) then
        c2 = 'IQU0'
      else if ( nstokes .eq. 4 ) then
        c2 = 'IQUV'
      endif

      OutUnit = 36
      OPEN(OutUnit,file='vlidort_v_test/results_brdfplus_tester_' &
           // trim(c1) // '_' // c2 // '.all',status='unknown')

!  Write Stokes Intensity and Mean values

      call write_stokes_brdf ( &
        OutUnit, 1, 'I', max_tasks, n_geometries, n_szangles, n_user_levels, &
        user_levels, par, eps, &
        STOKES_BAS,    MEANST_BAS,    FLUX_BAS, &
        STOKES_PT,     MEANST_PT,     FLUX_PT, &
        SURFWF_BAS, MEANST_SURFWF_BAS, FLUX_SURFWF_BAS )

!  Write Stokes Q and U and Mean values

      if ( NSTOKES .ge. 3 ) then
        call write_stokes_brdf ( &
          OutUnit, 2, 'Q', max_tasks, n_geometries, n_szangles, n_user_levels, &
          user_levels, par, eps, &
          STOKES_BAS,    MEANST_BAS,    FLUX_BAS, &
          STOKES_PT,     MEANST_PT,     FLUX_PT, &
          SURFWF_BAS, MEANST_SURFWF_BAS, FLUX_SURFWF_BAS )

        call write_stokes_brdf ( &
          OutUnit, 3, 'U', max_tasks, n_geometries, n_szangles, n_user_levels, &
          user_levels, par, eps, &
          STOKES_BAS,    MEANST_BAS,    FLUX_BAS, &
          STOKES_PT,     MEANST_PT,     FLUX_PT, &
          SURFWF_BAS, MEANST_SURFWF_BAS, FLUX_SURFWF_BAS )
      endif

!  Write Stokes V and Mean values

      if ( NSTOKES .eq. 4 ) then
        call write_stokes_brdf ( &
          OutUnit, 4, 'V', max_tasks, n_geometries, n_szangles, n_user_levels, &
          user_levels, par, eps, &
          STOKES_BAS,    MEANST_BAS,    FLUX_BAS, &
          STOKES_PT,     MEANST_PT,     FLUX_PT, &
          SURFWF_BAS, MEANST_SURFWF_BAS, FLUX_SURFWF_BAS )
      endif

!  Close file

      close(OutUnit)

!  Close VLIDORT error file if it has been opened

      write(*,*)
      IF ( OPENFILEFLAG ) then
        close( VLIDORT_ERRUNIT )
        write(*,*)'Main program executed with internal defaults,'// &
                  ' warnings in "V2p8p3_VLIDORT_Execution.log"'
      ELSE
        write(*,*)'Main program finished successfully'
      ENDIF

      stop

!  BRDF error finish

678   continue

      write(*,*)
      write(*,'(1x,a)') trim(TRACE)
      write(*,*)'Number of error messages from VBRDF calculation = ',VBRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES
      do i = 1, VBRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES
         write(*,'(a,i2,a,a)')' Message # ', i, ': ',adjustl(trim(VBRDF_Sup_OutputStatus%BS_OUTPUTMESSAGES(i)))
      enddo

      stop
      end program BRDFplus_Tester

!**********************************************************************
      subroutine write_stokes_brdf ( &
        un, o1, a1, max_tasks, n_geometries, nbeams, n_user_levels, &
        user_levels, par, eps, &
        STOKES_BAS,    MEANST_BAS,    FLUX_BAS, &
        STOKES_PT,     MEANST_PT,     FLUX_PT, &
        SURFWF_BAS, MEANST_SURFWF_BAS, FLUX_SURFWF_BAS )

      use vlidort_pars_m, Only : max_surfacewfs, max_user_levels, max_geometries, &
                                 max_szangles, maxstokes, max_directions, &
                                 upidx, dnidx

      implicit none

!  input

      integer, intent(in) ::           max_tasks, n_geometries, nbeams, n_user_levels
      integer, intent(in) ::           un, o1
      character (len=1), intent(in) :: a1
      double precision, intent(in) ::  user_levels ( max_user_levels )
      double precision, intent(in) ::  par(6)
      double precision, intent(in) ::  eps

      DOUBLE PRECISION, INTENT(IN) :: STOKES_BAS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: MEANST_BAS &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: FLUX_BAS &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT(IN) :: SURFWF_BAS &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: MEANST_SURFWF_BAS &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: FLUX_SURFWF_BAS &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT(IN) :: STOKES_PT &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION, INTENT(IN) :: MEANST_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION, INTENT(IN) :: FLUX_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )

!  local variables

      integer :: n, v

!  Write stokes component

      write(un,'(/T32,a/T32,a/)') &
       'STOKES-'//a1//', 6 SURFACE JACOBIANS', &
       '============================='
      write(un,'(a,T32,a/)')'Geometry    Level/Output', &
        'Stokes-'//a1//'       Surface AJ1   Surface FD1     '// &
                              'Surface AJ2   Surface FD2     '// &
                              'Surface AJ3   Surface FD3     '// &
                              'Surface AJ4   Surface FD4     '// &
                              'Surface AJ5   Surface FD5     '// &
                              'Surface AJ6   Surface FD6'

      do v = 1, n_geometries
         do n = 1, n_user_levels
           write(un,367)v,'Upwelling @',user_levels(n), &
           STOKES_BAS(n,v,o1,upidx), &
           par(1) * SURFWF_BAS(1,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,1)-STOKES_BAS(n,v,o1,upidx))/eps, &
           par(2) * SURFWF_BAS(2,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,2)-STOKES_BAS(n,v,o1,upidx))/eps, &
           par(3) * SURFWF_BAS(3,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,3)-STOKES_BAS(n,v,o1,upidx))/eps, &
           par(4) * SURFWF_BAS(4,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,4)-STOKES_BAS(n,v,o1,upidx))/eps, &
           par(5) * SURFWF_BAS(5,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,5)-STOKES_BAS(n,v,o1,upidx))/eps, &
           par(6) * SURFWF_BAS(6,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,6)-STOKES_BAS(n,v,o1,upidx))/eps
         enddo

         do n = 1, n_user_levels
           write(un,367)v,'Dnwelling @',user_levels(n), &
           STOKES_BAS(n,v,o1,dnidx), &
           par(1) * SURFWF_BAS(1,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,1)-STOKES_BAS(n,v,o1,dnidx))/eps, &
           par(2) * SURFWF_BAS(2,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,2)-STOKES_BAS(n,v,o1,dnidx))/eps, &
           par(3) * SURFWF_BAS(3,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,3)-STOKES_BAS(n,v,o1,dnidx))/eps, &
           par(4) * SURFWF_BAS(4,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,4)-STOKES_BAS(n,v,o1,dnidx))/eps, &
           par(5) * SURFWF_BAS(5,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,5)-STOKES_BAS(n,v,o1,dnidx))/eps, &
           par(6) * SURFWF_BAS(6,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,6)-STOKES_BAS(n,v,o1,dnidx))/eps
         enddo
         write(un,*)' '
      enddo

!  Write actinic fluxes

      write(un,'(/T32,a/T32,a/)') &
       'ACT.FLUX-'//a1//', 6 SURFACE JACOBIANS', &
       '==============================='
      write(un,'(a,T32,a/)')'  SunZA #   Level/Output', &
        'Act.Flux-'//a1//'     Surface AJ1   Surface FD1     '// &
                              'Surface AJ2   Surface FD2     '// &
                              'Surface AJ3   Surface FD3     '// &
                              'Surface AJ4   Surface FD4     '// &
                              'Surface AJ5   Surface FD5     '// &
                              'Surface AJ6   Surface FD6'

      do v = 1, nbeams
         do n = 1, n_user_levels
           write(un,367) v,'Upwelling @',user_levels(n), &
           MEANST_BAS(n,v,o1,upidx), &
           par(1) * MEANST_SURFWF_BAS(1,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,1)-MEANST_BAS(n,v,o1,upidx))/eps, &
           par(2) * MEANST_SURFWF_BAS(2,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,2)-MEANST_BAS(n,v,o1,upidx))/eps, &
           par(3) * MEANST_SURFWF_BAS(3,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,3)-MEANST_BAS(n,v,o1,upidx))/eps, &
           par(4) * MEANST_SURFWF_BAS(4,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,4)-MEANST_BAS(n,v,o1,upidx))/eps, &
           par(5) * MEANST_SURFWF_BAS(5,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,5)-MEANST_BAS(n,v,o1,upidx))/eps, &
           par(6) * MEANST_SURFWF_BAS(6,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,6)-MEANST_BAS(n,v,o1,upidx))/eps
         enddo

         do n = 1, n_user_levels
           write(un,367)v,'Dnwelling @',user_levels(n), &
           MEANST_BAS(n,v,o1,upidx), &
           par(1) * MEANST_SURFWF_BAS(1,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,1)-MEANST_BAS(n,v,o1,dnidx))/eps, &
           par(2) * MEANST_SURFWF_BAS(2,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,2)-MEANST_BAS(n,v,o1,dnidx))/eps, &
           par(3) * MEANST_SURFWF_BAS(3,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,3)-MEANST_BAS(n,v,o1,dnidx))/eps, &
           par(4) * MEANST_SURFWF_BAS(4,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,4)-MEANST_BAS(n,v,o1,dnidx))/eps, &
           par(5) * MEANST_SURFWF_BAS(5,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,5)-MEANST_BAS(n,v,o1,dnidx))/eps, &
           par(6) * MEANST_SURFWF_BAS(6,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,6)-MEANST_BAS(n,v,o1,dnidx))/eps
         enddo
         write(un,*)' '
      enddo

!  Write regular fluxes

      write(un,'(/T32,a/T32,a/)') &
       'REG.FLUX-'//a1//' , 6 SURFACE JACOBIANS', &
       '================================'
      write(un,'(a,T32,a/)')'  SunZA #   Level/Output', &
        'Reg.Flux-'//a1//'     Surface AJ1   Surface FD1     '// &
                              'Surface AJ2   Surface FD2     '// &
                              'Surface AJ3   Surface FD3     '// &
                              'Surface AJ4   Surface FD4     '// &
                              'Surface AJ5   Surface FD5     '// &
                              'Surface AJ6   Surface FD6'

      do v = 1, nbeams
         do n = 1, n_user_levels
           write(un,367) v,'Upwelling @',user_levels(n), &
           FLUX_BAS(n,v,o1,upidx), &
           par(1) * flux_SURFWF_BAS(1,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,1)-FLUX_BAS(n,v,o1,upidx))/eps, &
           par(2) * flux_SURFWF_BAS(2,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,2)-FLUX_BAS(n,v,o1,upidx))/eps, &
           par(3) * flux_SURFWF_BAS(3,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,3)-FLUX_BAS(n,v,o1,upidx))/eps, &
           par(4) * flux_SURFWF_BAS(4,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,4)-FLUX_BAS(n,v,o1,upidx))/eps, &
           par(5) * flux_SURFWF_BAS(5,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,5)-FLUX_BAS(n,v,o1,upidx))/eps, &
           par(6) * flux_SURFWF_BAS(6,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,6)-FLUX_BAS(n,v,o1,upidx))/eps
         enddo

         do n = 1, n_user_levels
           write(un,367)v,'Dnwelling @',user_levels(n), &
           FLUX_BAS(n,v,o1,upidx), &
           par(1) * flux_SURFWF_BAS(1,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,1)-FLUX_BAS(n,v,o1,dnidx))/eps, &
           par(2) * flux_SURFWF_BAS(2,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,2)-FLUX_BAS(n,v,o1,dnidx))/eps, &
           par(3) * flux_SURFWF_BAS(3,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,3)-FLUX_BAS(n,v,o1,dnidx))/eps, &
           par(4) * flux_SURFWF_BAS(4,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,4)-FLUX_BAS(n,v,o1,dnidx))/eps, &
           par(5) * flux_SURFWF_BAS(5,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,5)-FLUX_BAS(n,v,o1,dnidx))/eps, &
           par(6) * flux_SURFWF_BAS(6,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,6)-FLUX_BAS(n,v,o1,dnidx))/eps
         enddo
         write(un,*)' '
      enddo

 367  format(i5,T11,a,f6.2,2x,1pe13.5,12(2x,1p2e14.5))

      end subroutine write_stokes_brdf
