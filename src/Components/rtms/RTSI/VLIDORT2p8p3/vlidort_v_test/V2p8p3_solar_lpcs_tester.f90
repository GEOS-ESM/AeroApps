
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

      program Solar_lpcs_Tester

!  This is the Version 2.8 driver. Created 9/18/16 from 2.7 Driver
!     R. Spurr. RT Solutions Inc.

!  Upgrade for Version 2.8.1, August 2019
!  ---------------------------------------

!  Module files for VLIDORT. Strict usage, Version 2.8 upwards

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_LIN_IO_DEFS_m

      USE VLIDORT_AUX_m,      Only : VLIDORT_READ_ERROR, VLIDORT_WRITE_STATUS
      USE VLIDORT_INPUTS_m,   Only : VLIDORT_Sup_Init
      USE VLIDORT_L_INPUTS_m, Only : VLIDORT_L_INPUT_MASTER, VLIDORT_LinSup_Init
      USE VLIDORT_LPS_MASTERS_m
      USE VLIDORT_LCS_MASTERS_m

      IMPLICIT NONE

!  VLIDORT file inputs status structure

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

!  Local Variables
!  ===============

!  Flag for opening error output file

      LOGICAL ::          OPENFILEFLAG, PART1, PART2

!  Help variables
!  mick chg 4/3/2015 - renamed thread --> task

      INTEGER ::          NGREEKMAT_ENTRIES, L, LDUM, LCHECK
      INTEGER ::          K, CK, GK, RK, CMASK(8), GMASK(8), SMASK(8), RMASK(8)

      INTEGER ::          N,N6,NDUM,V,T,T1,T2,T3,T4,T5,T6,T7
      INTEGER ::          O1,Q,TASK, UTA,K3,K4,K5, KN
      DOUBLE PRECISION :: KD, GAER, WAER, TAER, PARCEL, RAYWT, AERWT
      DOUBLE PRECISION :: AERSCA, AEREXT, MOLSCA, TOTSCA, TOTEXT, OMEGA
      DOUBLE PRECISION :: MOLOMG(MAXLAYERS),MOLEXT(MAXLAYERS),C0
      DOUBLE PRECISION :: EPS, EPSFAC, RATIO1, RATIO2, MOLABS
      DOUBLE PRECISION :: AMOMS(6), AERMOMS(6,0:MAXMOMENTS_INPUT)
      DOUBLE PRECISION :: RAYMOMS(0:2,MAXLAYERS), DEPOL, BETA2
      DOUBLE PRECISION :: PROBLEM_RAY(6,0:2), SK, CUTOFF, MOM

      CHARACTER(LEN=7) :: C1
      CHARACTER(LEN=4) :: C2

      INTEGER ::          OutUnit

!  Saved optical property inputs

      DOUBLE PRECISION :: OMEGA_TOTAL_INPUT_SAVE ( MAXLAYERS )
      DOUBLE PRECISION :: DELTAU_VERT_INPUT_SAVE ( MAXLAYERS )
      DOUBLE PRECISION :: GREEKMAT_TOTAL_INPUT_SAVE ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION :: LAMBERTIAN_ALBEDO_SAVE

!  VLIDORT standard input preparation
!   Rob 3/22/15. Introduced FO_CALC local variable, renamed FOCORR 3/1/17

      LOGICAL ::          DO_LAMBERTIAN_SURFACE
      LOGICAL ::          DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY
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

      DOUBLE PRECISION :: FLUX_FACTOR

!  VLIDORT standard output preparation

      INTEGER ::          N_SZANGLES
      INTEGER ::          N_USER_LEVELS
      DOUBLE PRECISION :: USER_LEVELS ( MAX_USER_LEVELS )
      INTEGER ::          N_GEOMETRIES

!  VLIDORT linearized input preparation

      LOGICAL ::          DO_SIMULATION_ONLY
      LOGICAL ::          DO_PROFILE_LINEARIZATION
      LOGICAL ::          DO_COLUMN_LINEARIZATION
      LOGICAL ::          DO_ATMOS_LINEARIZATION
      INTEGER ::          N_TOTALPROFILE_WFS
      INTEGER ::          N_TOTALCOLUMN_WFS
      LOGICAL ::          LAYER_VARY_FLAG ( MAXLAYERS )
      INTEGER ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      DOUBLE PRECISION :: L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      LOGICAL ::          DO_SURFACE_LINEARIZATION
      INTEGER ::          N_SURFACE_WFS

!  Saved baseline results

      DOUBLE PRECISION :: STOKES_BAS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: MEANST_BAS &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FLUX_BAS &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )

      DOUBLE PRECISION :: PROFWF_SAVE &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: MEANST_PROFWF_SAVE &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FLUX_PROFWF_SAVE &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: COLWF_SAVE &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: MEANST_COLWF_SAVE &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FLUX_COLWF_SAVE &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: SURFWF_SAVE &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: MEANST_SURFWF_SAVE &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FLUX_SURFWF_SAVE &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Saved Perturbation results

      INTEGER, PARAMETER :: MAX_TASKS = 10

      DOUBLE PRECISION :: STOKES_PT &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION :: MEANST_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION :: FLUX_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )

!  Start of code
!  =============

!  Define some arrays
!  @@@ Rob Fix 6/13/13. CMASK is for the AEROSOLS, RMASK for RAYLEIGH.

      GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)
      RMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)    !  This   Rayleigh
!      CMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)
      CMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)    !  OUTPUT OF RTS MIE
      SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /)

!  Define some control variables

!  Profile WF calculation
!      PART1 = .FALSE.
      PART1 = .TRUE.

!  Column WF calculation
!      PART2 = .FALSE.
      PART2 = .TRUE.

!  Only interested in 1 test now. Version 2.8. 3/1/17

!mick add 3/25/2015 - added SS/FO test loop
!  Start single-scatter test loop
!    loop 1 - performs tests using original SSCOR code
!    loop 2 - performs tests using first-order (FO) code
!      DO test_ss=1,2
!      IF (test_ss == 1) THEN
!        write(*,*) 'test using SSCorr code' ; write(*,*) '---------------------'
!      ELSE IF (test_ss == 2) THEN
!        write(*,*) 'test using FO code' ; write(*,*) '------------------'
!      ENDIF
!      write(*,*)

!  Initialize error file output flag

      OPENFILEFLAG = .false.

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

!  Rob 3/22/15. Add wavelength Diagnostic (just a dummy here)
!   Versino 2.8, actually a file-read quantity, overwritten here

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = 0.333d0

!  Initialize some input variables not handled by VLIDORT_L_INPUT_MASTER

      CALL VLIDORT_Sup_Init ( VLIDORT_Sup )
      CALL VLIDORT_LinSup_Init ( VLIDORT_LinSup )

!  Set Output-level Proxies (saved values)

      nstokes       = VLIDORT_FixIn%Cont%TS_nstokes
      n_user_levels = VLIDORT_FixIn%UserVal%TS_n_user_levels
      user_levels   = VLIDORT_ModIn%MUserVal%TS_user_levels
      n_szangles    = VLIDORT_ModIn%MSunrays%TS_n_szangles

!  Set output file designation for later use

      do_observation_geometry = VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY
      do_doublet_geometry     = VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY

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

!  Set masking limits

      if ( nstokes .eq. 1 ) ngreekmat_entries = 1
      if ( nstokes .eq. 3 ) ngreekmat_entries = 5
      if ( nstokes .eq. 4 ) ngreekmat_entries = 8

!  Standard inputs - cross-over w thermal and surface emission
!   Updated, 2.8, 3/1/17.
!mick note 9/19/2107 - unlike the scalar version of the solar lpcs test, inputs in
!                      this section ARE active (override config file settings)

      DO_SOLAR_SOURCES      = .true.

      DO_FOCORR             = .true.
      DO_FOCORR_NADIR       = .true.
      DO_FOCORR_OUTGOING    = .false.

      DO_THERMAL_EMISSION   = .false.
      DO_THERMAL_TRANSONLY  = .false.
      DO_DELTAM_SCALING     = .true.
      DO_SOLUTION_SAVING    = .false.
      DO_BVP_TELESCOPING    = .false.
      DO_LAMBERTIAN_SURFACE = .true.
      DO_SURFACE_EMISSION   = .false.

      NFINELAYERS           = 4 ! Non zero here
      N_THERMAL_COEFFS      = 2

      FLUX_FACTOR           = 1.0d-03

!mick add 3/25/2015 - added FO control (by do loop)
!  This section removed, 3/1/17   for Version 2.8
!  Additional control for testing FO code
!      if (test_ss == 1) then
!         VLIDORT_ModIn%MBool%TS_DO_FO_CALC = .false.
!      else if (test_ss == 2) then
!        VLIDORT_ModIn%MBool%TS_DO_FO_CALC = .true.
!        n_user_levels = 3
!        user_levels(1:2) = VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:2)
!        user_levels(3)   = VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(5)
!        VLIDORT_ModIn%MUserVal%TS_USER_LEVELS  = zero
!        VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS = n_user_levels
!        VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:3) = user_levels(1:3)
!      endif

!  Skip part 1 of the tests if desired

      IF ( .NOT. PART1 ) GO TO 567

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!                     P A R T    I

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Baseline Profile WF calculation
!  ================================

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

!  Initialise

      deltau_vert_input    = 0.0d0
      omega_total_input    = 0.0d0
      greekmat_total_input = 0.0d0
      L_deltau_vert_input  = 0.0d0
      L_omega_total_input  = 0.0d0
      L_greekmat_total_input = 0.0d0

      write(*,*)'Doing PROFILEWF Baseline: '// &
      'Baseline, Intensities + 3 profile Jacobians + 1 Surface Jacobian'

!  Linearized inputs

      DO_COLUMN_LINEARIZATION  = .FALSE.
      DO_PROFILE_LINEARIZATION = .TRUE.
      DO_ATMOS_LINEARIZATION   = .TRUE.
      DO_SURFACE_LINEARIZATION = .TRUE.
      DO_SIMULATION_ONLY       = .FALSE.

      n_totalprofile_wfs       = 2
      do n = 1, nlayers
        layer_vary_number(n) = 2
        layer_vary_flag(n)   = .true.
      enddo
      n_surface_wfs            = 1

!  Fix the non-aerosol layers

      n6 = nlayers - 6
      do n = 1, n6
        deltau_vert_input(n) = molext(n)
        omega_total_input(n) = molomg(n)
        ratio1 = 1.0d0 - molomg(n)
        L_deltau_vert_input(1,n) =   ratio1
        L_omega_total_input(1,n) = - ratio1
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
        molabs = molext(n) * ( 1.0d0 - molomg(n) )
        molsca = molomg(n) * molext(n)
        totext = molext(n) + aerext
        totsca = molsca    + aersca
        raywt  = molsca / totsca
        aerwt  = aersca / totsca
        omega  = totsca / totext
        deltau_vert_input(n) = totext
        omega_total_input(n) = totsca / totext
        ratio1 = molabs / totext
        L_deltau_vert_input(1,n) =   ratio1
        L_omega_total_input(1,n) = - ratio1
        ratio2 = aerext / totext 
        L_deltau_vert_input(2,n) = ratio2
        L_omega_total_input(2,n) = ratio2 * ((waer/omega) - 1.0d0)
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
          do L = 1, ngreek_moments_input
            if ( greekmat_total_input(L,n,GK) .ne. 0.0d0 ) then
              ratio2 = aersca / totsca
              L_greekmat_total_input(2,L,n,GK) = ratio2 * &
            ( (SK*AERMOMS(CK,L)/greekmat_total_input(L,n,GK)) - 1.0d0 )
            endif
          enddo
        ENDDO
        GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0
      enddo

!  Surface

      lambertian_albedo = 0.05d0

!  Save inputs
!  -----------

      do n = 1, nlayers
       deltau_vert_input_save(n) = deltau_vert_input(n) 
       omega_total_input_save(n) = omega_total_input(n)
       do L = 0, ngreek_moments_input
        do k = 1, 16
         greekmat_total_input_save(L,n,k) = greekmat_total_input(L,n,k)
        enddo
       enddo
      enddo
      lambertian_albedo_save = lambertian_albedo

!  Copy local variables. Upgraded for Version 2.8. 3/1/17

      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION   = DO_THERMAL_EMISSION
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION   = DO_SURFACE_EMISSION
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = DO_LAMBERTIAN_SURFACE

      VLIDORT_ModIn%MBool%TS_DO_FOCORR            = DO_FOCORR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR      = DO_FOCORR_NADIR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING   = DO_FOCORR_OUTGOING

      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES     = DO_SOLAR_SOURCES
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING    = DO_DELTAM_SCALING
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING   = DO_SOLUTION_SAVING
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING   = DO_BVP_TELESCOPING
      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY = DO_THERMAL_TRANSONLY

      VLIDORT_FixIn%Cont%TS_NFINELAYERS           = NFINELAYERS
      VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS      = N_THERMAL_COEFFS
      VLIDORT_FixIn%Cont%TS_NLAYERS               = NLAYERS

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = NGREEK_MOMENTS_INPUT

      VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR        = FLUX_FACTOR

      VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID        = HEIGHT_GRID

!  Copy local linearization control

      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = DO_ATMOS_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = DO_SURFACE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = DO_SIMULATION_ONLY

      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)   = LAYER_VARY_FLAG(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS) = LAYER_VARY_NUMBER(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS           = N_TOTALPROFILE_WFS
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS                = N_SURFACE_WFS

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_lambertian_albedo    = lambertian_albedo
      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input    = L_deltau_vert_input
      VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input = L_greekmat_total_input
      VLIDORT_LinFixIn%Optical%TS_L_omega_total_input    = L_omega_total_input

!  VLIDORT call

      do_debug_input = .false.
      CALL VLIDORT_LPS_MASTER ( do_debug_input, &
        VLIDORT_FixIn, &
        VLIDORT_ModIn, &
        VLIDORT_Sup, &
        VLIDORT_Out, &
        VLIDORT_LinFixIn, &
        VLIDORT_LinModIn, &
        VLIDORT_LinSup, &
        VLIDORT_LinOut )

!  Exception handling, write-up (optional)

      CALL VLIDORT_WRITE_STATUS ( &
        'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Set some variables

      n_geometries  = VLIDORT_Out%Main%TS_n_geometries

!  Save baseline results

      do uta = 1, n_user_levels
        do o1 = 1, nstokes

          do v = 1, n_geometries
            stokes_BAS(uta,v,o1,upidx) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,upidx)
            stokes_BAS(uta,v,o1,dnidx) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,dnidx)
            do n = 1, nlayers
              do q = 1, n_totalprofile_wfs
                PROFWF_SAVE(q,n,uta,v,o1,upidx) = VLIDORT_LinOut%Prof%TS_profilewf(q,n,uta,v,o1,upidx)
                PROFWF_SAVE(q,n,uta,v,o1,dnidx) = VLIDORT_LinOut%Prof%TS_profilewf(q,n,uta,v,o1,dnidx)
              enddo
            enddo
            do q = 1, n_surface_wfs
              SURFWF_SAVE(q,uta,v,o1,upidx) = VLIDORT_LinOut%Surf%TS_surfacewf(q,uta,v,o1,upidx)
              SURFWF_SAVE(q,uta,v,o1,dnidx) = VLIDORT_LinOut%Surf%TS_surfacewf(q,uta,v,o1,dnidx)
            enddo
          enddo

          do v = 1, n_szangles
            MEANST_BAS(uta,v,o1,upidx) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,upidx)
            MEANST_BAS(uta,v,o1,dnidx) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,dnidx) &
                                       + VLIDORT_Out%Main%TS_DNMEANST_DIRECT(uta,v,o1)
            FLUX_BAS(uta,v,o1,upidx)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,upidx)
            FLUX_BAS(uta,v,o1,dnidx)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,dnidx) &
                                       + VLIDORT_Out%Main%TS_DNFLUX_DIRECT(uta,v,o1)

            do n = 1, nlayers
              do q = 1, n_totalprofile_wfs
                MEANST_PROFWF_save(q,n,uta,v,o1,upidx) = VLIDORT_LinOut%Prof%TS_MEANST_DIFFUSE_PROFWF(q,n,uta,v,o1,upidx)
                MEANST_PROFWF_save(q,n,uta,v,o1,dnidx) = VLIDORT_LinOut%Prof%TS_MEANST_DIFFUSE_PROFWF(q,n,uta,v,o1,dnidx) &
                                                       + VLIDORT_LinOut%Prof%TS_DNMEANST_DIRECT_PROFWF(q,n,uta,v,o1)
                FLUX_PROFWF_SAVE(q,n,uta,v,o1,upidx)   = VLIDORT_LinOut%Prof%TS_FLUX_DIFFUSE_PROFWF(q,n,uta,v,o1,upidx)
                FLUX_PROFWF_SAVE(q,n,uta,v,o1,dnidx)   = VLIDORT_LinOut%Prof%TS_FLUX_DIFFUSE_PROFWF(q,n,uta,v,o1,dnidx) &
                                                       + VLIDORT_LinOut%Prof%TS_DNFLUX_DIRECT_PROFWF(q,n,uta,v,o1)
              enddo
            enddo

            do q = 1, n_surface_wfs
              MEANST_SURFWF_SAVE(q,uta,v,o1,upidx) = VLIDORT_LinOut%Surf%TS_MEANST_DIFFUSE_SURFWF(q,uta,v,o1,upidx)
              MEANST_SURFWF_SAVE(q,uta,v,o1,dnidx) = VLIDORT_LinOut%Surf%TS_MEANST_DIFFUSE_SURFWF(q,uta,v,o1,dnidx)
              FLUX_SURFWF_SAVE(q,uta,v,o1,upidx)   = VLIDORT_LinOut%Surf%TS_FLUX_DIFFUSE_SURFWF(q,uta,v,o1,upidx)
              FLUX_SURFWF_SAVE(q,uta,v,o1,dnidx)   = VLIDORT_LinOut%Surf%TS_FLUX_DIFFUSE_SURFWF(q,uta,v,o1,dnidx)
            enddo
          enddo

        enddo
      enddo

      !stop 'baseline profile WFS only'

!  Initialize FD section
!  =====================

!  Fd perturbation

      eps = 1.0d-03
      epsfac = 1.0d0 + eps

!  Initialize linearized inputs

      L_deltau_vert_input = 0.0d0
      L_omega_total_input = 0.0d0
      L_greekmat_total_input = 0.0d0

      DO_COLUMN_LINEARIZATION  = .FALSE.
      DO_PROFILE_LINEARIZATION = .FALSE.
      DO_ATMOS_LINEARIZATION   = .FALSE.
      DO_SURFACE_LINEARIZATION = .FALSE.
      DO_SIMULATION_ONLY       = .TRUE.

      n_totalprofile_wfs      = 0
      n_surface_wfs           = 0
      do n = 1, nlayers
        layer_vary_number(n) = 0
        layer_vary_flag(n)   = .false.
      enddo

!  Copy local linearization control

      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = DO_ATMOS_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = DO_SURFACE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = DO_SIMULATION_ONLY

      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)   = LAYER_VARY_FLAG(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS) = LAYER_VARY_NUMBER(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS           = N_TOTALPROFILE_WFS
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS                = N_SURFACE_WFS

!  TASKS for FD testing
!  ====================

!  Task 1, FD Perturbation on Surface Albedo
!  Task 2, FD Perturbation on trace gas optical depth in layer 6
!  Task 3, FD Perturbation on trace gas optical depth in layer 21
!  Task 4, FD Perturbation on aerosol optical depth   in layer 23

      do task = 1, 4
!      do task = 1, 1

        if ( task.eq.2 ) kn = 6
        if ( task.eq.3 ) kn = 21
        if ( task.eq.4 ) kn = 23

!  task counter

        t = task
        deltau_vert_input    = 0.0d0
        omega_total_input    = 0.0d0
        greekmat_total_input = 0.0d0

!  Task 1

        if ( task .eq. 1 ) then
          write(*,*)'Doing task 1: '// 'Finite Difference, perturb Lambertian albedo'

          lambertian_albedo = lambertian_albedo_save * epsfac
          do n = 1, nlayers
            deltau_vert_input(n) = deltau_vert_input_save(n)
            omega_total_input(n) = omega_total_input_save(n)
            do L = 0, ngreek_moments_input
              do k = 1, 16
                greekmat_total_input(L,n,k) = greekmat_total_input_save(L,n,k)
              enddo
            enddo
          enddo
        endif

!  Task 2 or 3

        if ( task.eq.2 .or. task.eq.3  ) then
          write(*,'(a,i2,a,i2)')' Doing task',task, ': Finite Difference, perturb molecular absorption Layer ',kn

          !  Surface
          lambertian_albedo = lambertian_albedo_save

          !  rayleigh layers
          do n = 1, n6
            molabs = molext(n) * ( 1.0d0 - molomg(n) )
            if ( n.eq.kn ) molabs = molabs * epsfac
            molsca = molomg(n) * molext(n)
            deltau_vert_input(n) = molabs + molsca
            omega_total_input(n) = molsca / ( molabs + molsca )
            do L = 0, ngreek_moments_input
              do k = 1, 16
                greekmat_total_input(L,n,k) = greekmat_total_input_save(L,n,k)
              enddo
            enddo
          enddo

          !  aerosol layers
          do n = n6 + 1, nlayers
            aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
            aersca = aerext * waer
            molabs = molext(n) * ( 1.0d0 - molomg(n) )
            if ( n.eq.kn ) molabs = molabs * epsfac
            molsca = molomg(n) * molext(n)
            totext = molabs + molsca + aerext
            totsca = molsca + aersca
            raywt  = molsca / totsca
            aerwt  = aersca / totsca
            omega  = totsca / totext
            deltau_vert_input(n) = totext
            omega_total_input(n) = omega
            do L = 0, ngreek_moments_input
              do k = 1, 16
                greekmat_total_input(L,n,k) = greekmat_total_input_save(L,n,k)
              enddo
            enddo
          enddo
        endif

!  Task 4

        if ( task.eq.4  ) then
          write(*,'(a,i2,a,i2)')' Doing task',task, &
            ': Finite Difference, perturb aerosol optical depth Layer ',kn

          !  Surface and Rayleigh layers
          lambertian_albedo = lambertian_albedo_save
          do n = 1, nlayers
            deltau_vert_input(n) = deltau_vert_input_save(n)
            omega_total_input(n) = omega_total_input_save(n)
            do L = 0, ngreek_moments_input
              do k = 1, 16
                greekmat_total_input(L,n,k) = greekmat_total_input_save(L,n,k)
              enddo
            enddo
          enddo

          !  aerosol layers
          parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
          do n = n6 + 1, nlayers
            aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
            if ( kn .eq. n ) aerext = aerext * epsfac
            aersca = aerext * waer
            molabs = molext(n) * ( 1.0d0 - molomg(n) )
            molsca = molomg(n) * molext(n)
            totext = molabs + molsca + aerext
            totsca = molsca + aersca
            raywt  = molsca / totsca
            aerwt  = aersca / totsca
            omega  = totsca / totext
            deltau_vert_input(n) = totext
            omega_total_input(n) = omega
            DO K = 1, ngreekmat_entries
              GK = gmask(k); sk = dble(smask(k)) ; ck = cmask(k) ; rk = rmask(k)
              DO L = 0, 2
! Old           MOM = RAYWT*PROBLEM_RAY(ck,L) + AERWT*SK*AERMOMS(CK,L)
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
        endif

!  Copy to optical property type-structure inputs

        VLIDORT_FixIn%Optical%TS_lambertian_albedo    = lambertian_albedo
        VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
        VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
        VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

!  VLIDORT call

        do_debug_input = .false.
        CALL VLIDORT_LPS_MASTER ( do_debug_input, &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup, &
          VLIDORT_Out, &
          VLIDORT_LinFixIn, &
          VLIDORT_LinModIn, &
          VLIDORT_LinSup, &
          VLIDORT_LinOut )

!  Exception handling, write-up (optional)

        CALL VLIDORT_WRITE_STATUS ( &
          'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Set some variables

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

!  Write files of Profile WF test results
!  ======================================

!  Shorthand indices

      k3 = 6 ; k4 = 21 ; k5 = 23
      t1 = 1 ; t2 = 2  ; t3 = 3 ; t4 = 4

!  Superceded code from Version 2.7
!      c3 = 'FO'
!      if ( .not. VLIDORT_ModIn%MBool%TS_DO_FO_CALC ) then
!        OPEN(OutUnit,file='vlidort_v_test/results_solar_lps_tester_' &
!             // c1 // '_' // c2 // '.all',status='unknown')
!      else
!        OPEN(OutUnit,file='vlidort_v_test/results_solar_lps_tester_' &
!             // c1 // '_' // c2 // '_' // c3 //'.all',status='unknown')
!      end if

!  Version 2.8 only the FO calculation

      OutUnit = 36
      OPEN(OutUnit,file='vlidort_v_test/results_solar_lps_tester_' &
           // trim(c1) // '_' // c2 // '.all',status='unknown')

!  Write Stokes Intensity and Mean values

      call write_stokes_lps ( &
        OutUnit, 1, 'I', k3, k4, k5, t1, t2, t3, t4, &
        max_tasks, n_geometries, n_szangles,  n_user_levels, &
        user_levels, lambertian_albedo_save, eps, &
        STOKES_BAS, MEANST_BAS, FLUX_BAS, &
        STOKES_PT,  MEANST_PT,  FLUX_PT, &
        PROFWF_SAVE, SURFWF_SAVE, &
        MEANST_PROFWF_SAVE, FLUX_PROFWF_SAVE, &
        MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )

!  Write Stokes Q and U and Mean values

      if ( NSTOKES .ge. 3 ) then
        call write_stokes_lps ( &
          OutUnit, 2, 'Q', k3, k4, k5, t1, t2, t3, t4, &
          max_tasks, n_geometries, n_szangles,  n_user_levels, &
          user_levels, lambertian_albedo_save, eps, &
          STOKES_BAS, MEANST_BAS, FLUX_BAS, &
          STOKES_PT,  MEANST_PT,  FLUX_PT, &
          PROFWF_SAVE, SURFWF_SAVE, &
          MEANST_PROFWF_SAVE, FLUX_PROFWF_SAVE, &
          MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )

        call write_stokes_lps ( &
          OutUnit, 3, 'U', k3, k4, k5, t1, t2, t3, t4, &
          max_tasks, n_geometries, n_szangles,  n_user_levels, &
          user_levels, lambertian_albedo_save, eps, &
          STOKES_BAS, MEANST_BAS, FLUX_BAS, &
          STOKES_PT,  MEANST_PT,  FLUX_PT, &
          PROFWF_SAVE, SURFWF_SAVE, &
          MEANST_PROFWF_SAVE, FLUX_PROFWF_SAVE, &
          MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )
      endif

!  Write Stokes V and Mean values

      if ( NSTOKES .eq. 4 ) then
        call write_stokes_lps ( &
          OutUnit, 4, 'V', k3, k4, k5, t1, t2, t3, t4, &
          max_tasks, n_geometries, n_szangles,  n_user_levels, &
          user_levels, lambertian_albedo_save, eps, &
          STOKES_BAS, MEANST_BAS, FLUX_BAS, &
          STOKES_PT,  MEANST_PT,  FLUX_PT, &
          PROFWF_SAVE,      SURFWF_SAVE, &
          MEANST_PROFWF_SAVE, FLUX_PROFWF_SAVE, &
          MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )
      endif

      close(OutUnit)

!  Continuation point

 567  continue

      IF ( .NOT. PART2 ) STOP 'Program ended because PART 2 tests are not flagged'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!                     P A R T    I I

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Baseline Column WF calculation
!  ==============================

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

!  Initialize

      deltau_vert_input      = 0.0d0
      omega_total_input      = 0.0d0
      greekmat_total_input   = 0.0d0
      L_deltau_vert_input    = 0.0d0
      L_omega_total_input    = 0.0d0
      L_greekmat_total_input = 0.0d0

      IF ( PART1 ) write(*,*)
      write(*,*)'Doing COLUMNWF Baseline: '// &
       'Baseline, Intensities + 2 Column Jacobians + 1 Surface Jacobian'

!  Initialize linearized inputs

      DO_COLUMN_LINEARIZATION  = .TRUE.
      DO_PROFILE_LINEARIZATION = .FALSE.
      DO_ATMOS_LINEARIZATION   = .TRUE.
      DO_SURFACE_LINEARIZATION = .TRUE.
      DO_SIMULATION_ONLY       = .FALSE.

      n_totalcolumn_wfs        = 2
      do n = 1, nlayers
        layer_vary_number(n) = 2
        layer_vary_flag(n)   = .true.
      enddo
      n_surface_wfs            = 1

!  Total column

      c0 = 0.0d0
      do n = 1, nlayers
        molabs = molext(n) * ( 1.0d0 - molomg(n) )
        c0 = c0 + molabs
      enddo

!  Fix the non-aerosol layers

      n6 = nlayers - 6
      do n = 1, n6
        deltau_vert_input(n) = molext(n)
        omega_total_input(n) = molomg(n)
        ratio1 = 1.0d0 - molomg(n)
        L_deltau_vert_input(1,n) =   ratio1
        L_omega_total_input(1,n) = - ratio1
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
        molabs = molext(n) * ( 1.0d0 - molomg(n) )
        molsca = molomg(n) * molext(n)
        totext = molext(n) + aerext
        totsca = molsca    + aersca
        raywt  = molsca / totsca
        aerwt  = aersca / totsca
        omega  = totsca / totext
        deltau_vert_input(n) = totext
        omega_total_input(n) = omega
        ratio1 = molabs / totext
        L_deltau_vert_input(1,n) =   ratio1
        L_omega_total_input(1,n) = - ratio1
        ratio2 = aerext / totext 
        L_deltau_vert_input(2,n) = ratio2
        L_omega_total_input(2,n) = ratio2 * ((waer/omega) - 1.0d0)
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
          do L = 1, ngreek_moments_input
            if ( greekmat_total_input(L,n,GK) .ne. 0.0d0 ) then
              ratio2 = aersca / totsca
              L_greekmat_total_input(2,L,n,GK) = ratio2 * &
            ( (SK*AERMOMS(CK,L)/greekmat_total_input(L,n,GK)) - 1.0d0 )
            endif
          enddo
        ENDDO
        GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0
      enddo

!  Surface

      lambertian_albedo = 0.05d0

!  Save inputs
!  -----------

      do n = 1, nlayers
       deltau_vert_input_save(n) = deltau_vert_input(n) 
       omega_total_input_save(n) = omega_total_input(n)
       do L = 0, ngreek_moments_input
        do k = 1, 16
         greekmat_total_input_save(L,n,k) = greekmat_total_input(L,n,k)
        enddo
       enddo
      enddo
      lambertian_albedo_save = lambertian_albedo

!  Copy local variables. List of FOCORR updated for Version 2.8, 3/1/17

      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION   = DO_THERMAL_EMISSION
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION   = DO_SURFACE_EMISSION
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = DO_LAMBERTIAN_SURFACE

      VLIDORT_ModIn%MBool%TS_DO_FOCORR            = DO_FOCORR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR      = DO_FOCORR_NADIR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING   = DO_FOCORR_OUTGOING

      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES     = DO_SOLAR_SOURCES

      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING    = DO_DELTAM_SCALING
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING   = DO_SOLUTION_SAVING
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING   = DO_BVP_TELESCOPING
      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY = DO_THERMAL_TRANSONLY

      VLIDORT_FixIn%Cont%TS_NFINELAYERS           = NFINELAYERS
      VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS      = N_THERMAL_COEFFS
      VLIDORT_FixIn%Cont%TS_NLAYERS               = NLAYERS

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = NGREEK_MOMENTS_INPUT

      VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR        = FLUX_FACTOR

      VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID        = HEIGHT_GRID

!  Copy local linearization control

      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = DO_ATMOS_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = DO_SURFACE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = DO_SIMULATION_ONLY

      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)   = LAYER_VARY_FLAG(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS) = LAYER_VARY_NUMBER(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS            = N_TOTALCOLUMN_WFS
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS                = N_SURFACE_WFS

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_lambertian_albedo    = lambertian_albedo
      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input    = L_deltau_vert_input
      VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input = L_greekmat_total_input
      VLIDORT_LinFixIn%Optical%TS_L_omega_total_input    = L_omega_total_input

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

      CALL VLIDORT_WRITE_STATUS ( &
        'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Set some variables

      n_geometries  = VLIDORT_Out%Main%TS_n_geometries

!  Save baseline results

      do uta = 1, n_user_levels
        do o1 = 1, nstokes

          do v = 1, n_geometries
            stokes_BAS(uta,v,o1,upidx) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,upidx)
            stokes_BAS(uta,v,o1,dnidx) = VLIDORT_Out%Main%TS_stokes(uta,v,o1,dnidx)
            do q = 1, n_totalcolumn_wfs
              COLWF_SAVE(q,uta,v,o1,upidx) = VLIDORT_LinOut%Col%TS_columnwf(q,uta,v,o1,upidx)
              COLWF_SAVE(q,uta,v,o1,dnidx) = VLIDORT_LinOut%Col%TS_columnwf(q,uta,v,o1,dnidx)
            enddo
            do q = 1, n_surface_wfs
              SURFWF_SAVE(q,uta,v,o1,upidx) = VLIDORT_LinOut%Surf%TS_surfacewf(q,uta,v,o1,upidx)
              SURFWF_SAVE(q,uta,v,o1,dnidx) = VLIDORT_LinOut%Surf%TS_surfacewf(q,uta,v,o1,dnidx)
            enddo
          enddo

          do v = 1, n_szangles
            MEANST_BAS(uta,v,o1,upidx) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,upidx)
            MEANST_BAS(uta,v,o1,dnidx) = VLIDORT_Out%Main%TS_MEANST_DIFFUSE(uta,v,o1,dnidx) &
                                       + VLIDORT_Out%Main%TS_DNMEANST_DIRECT(uta,v,o1)
            FLUX_BAS(uta,v,o1,upidx)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,upidx)
            FLUX_BAS(uta,v,o1,dnidx)   = VLIDORT_Out%Main%TS_FLUX_DIFFUSE(uta,v,o1,dnidx) &
                                       + VLIDORT_Out%Main%TS_DNFLUX_DIRECT(uta,v,o1)

            do q = 1, n_totalcolumn_wfs
              MEANST_COLWF_save(q,uta,v,o1,upidx) = VLIDORT_LinOut%Col%TS_MEANST_DIFFUSE_COLWF(q,uta,v,o1,upidx)
              MEANST_COLWF_save(q,uta,v,o1,dnidx) = VLIDORT_LinOut%Col%TS_MEANST_DIFFUSE_COLWF(q,uta,v,o1,dnidx) &
                                                  + VLIDORT_LinOut%Col%TS_DNMEANST_DIRECT_COLWF(q,uta,v,o1)
              FLUX_COLWF_SAVE(q,uta,v,o1,upidx)   = VLIDORT_LinOut%Col%TS_FLUX_DIFFUSE_COLWF(q,uta,v,o1,upidx)
              FLUX_COLWF_SAVE(q,uta,v,o1,dnidx)   = VLIDORT_LinOut%Col%TS_FLUX_DIFFUSE_COLWF(q,uta,v,o1,dnidx) &
                                                  + VLIDORT_LinOut%Col%TS_DNFLUX_DIRECT_COLWF(q,uta,v,o1)
            enddo

            do q = 1, n_surface_wfs
              MEANST_SURFWF_SAVE(q,uta,v,o1,upidx) = VLIDORT_LinOut%Surf%TS_MEANST_DIFFUSE_SURFWF(q,uta,v,o1,upidx)
              MEANST_SURFWF_SAVE(q,uta,v,o1,dnidx) = VLIDORT_LinOut%Surf%TS_MEANST_DIFFUSE_SURFWF(q,uta,v,o1,dnidx)
              FLUX_SURFWF_SAVE(q,uta,v,o1,upidx)   = VLIDORT_LinOut%Surf%TS_FLUX_DIFFUSE_SURFWF(q,uta,v,o1,upidx)
              FLUX_SURFWF_SAVE(q,uta,v,o1,dnidx)   = VLIDORT_LinOut%Surf%TS_FLUX_DIFFUSE_SURFWF(q,uta,v,o1,dnidx)
            enddo
          enddo

        enddo
      enddo

!  Initialize FD section
!  =====================

!  Fd perturbation

      eps    = 1.0d-02
      epsfac = 1.0d0 + eps

!  Initialize linearized inputs

      L_deltau_vert_input = 0.0d0
      L_omega_total_input = 0.0d0
      L_greekmat_total_input = 0.0d0

      DO_COLUMN_LINEARIZATION  = .FALSE.
      DO_PROFILE_LINEARIZATION = .FALSE.
      DO_ATMOS_LINEARIZATION   = .FALSE.
      DO_SURFACE_LINEARIZATION = .FALSE.
      DO_SIMULATION_ONLY       = .TRUE.

      n_totalcolumn_wfs = 0
      n_surface_wfs = 0
      do n = 1, nlayers
        layer_vary_number(n) = 0
        layer_vary_flag(n)   = .false.
      enddo

!  Copy local linearization control

      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = DO_ATMOS_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = DO_SURFACE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = DO_SIMULATION_ONLY

      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)   = LAYER_VARY_FLAG(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS) = LAYER_VARY_NUMBER(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS            = N_TOTALCOLUMN_WFS
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS                = N_SURFACE_WFS

!  TASKS for FD testing
!  ====================

!  Task 5, FD Perturbation on Surface Albedo
!  Task 6, FD Perturbation on Total trace gas optical depth C0
!  Task 7, FD Perturbation on Total aerosol optical depth TAU

      do task = 5, 7
!      do task = 6, 6
!      do task = 7, 7

!  task counter

        t = task
        deltau_vert_input    = 0.0d0
        omega_total_input    = 0.0d0
        greekmat_total_input = 0.0d0

!  Task 5

        if ( task .eq. 5 ) then
          write(*,*)'Doing task 5: '// 'Finite Difference, perturb Lambertian albedo'

          lambertian_albedo = lambertian_albedo_save * epsfac
          do n = 1, nlayers
            deltau_vert_input(n) = deltau_vert_input_save(n) 
            omega_total_input(n) = omega_total_input_save(n) 
            do L = 0, ngreek_moments_input
              do k = 1, 16
                greekmat_total_input(L,n,k) = greekmat_total_input_save(L,n,k)
              enddo
            enddo
          enddo
        endif

!  Task 6

        if ( task .eq. 6 ) then
          write(*,*)'Doing task 6: '// 'Finite Difference, perturb total column absorber'

          !  Surface and rayleigh layers
          lambertian_albedo = lambertian_albedo_save
          do n = 1, n6
            molabs = molext(n) * ( 1.0d0 - molomg(n) ) * epsfac
            molsca = molomg(n) * molext(n)
            deltau_vert_input(n) = molabs + molsca
            omega_total_input(n) = molsca / ( molabs + molsca )
            do L = 0, ngreek_moments_input
              do k = 1, 16
                greekmat_total_input(L,n,k) = greekmat_total_input_save(L,n,k)
              enddo
            enddo
          enddo

          !  aerosol layers
          do n = n6 + 1, nlayers
            aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
            aersca = aerext * waer
            molabs = molext(n) * ( 1.0d0 - molomg(n) ) * epsfac
            molsca = molomg(n) * molext(n)
            totext = molabs + molsca + aerext
            totsca = molsca + aersca
            raywt  = molsca / totsca
            aerwt  = aersca / totsca
            omega  = totsca / totext
            deltau_vert_input(n) = totext
            omega_total_input(n) = omega
            do L = 0, ngreek_moments_input
              do k = 1, 16
                greekmat_total_input(L,n,k) = greekmat_total_input_save(L,n,k)
              enddo
            enddo
          enddo
        endif

!  Task 7

        if ( task .eq. 7 ) then
          write(*,*)'Doing task 7: '// 'Finite Difference, perturb total aerosol tau'

          !  Surface and rayleigh layers
          lambertian_albedo = lambertian_albedo_save
          do n = 1, n6
            deltau_vert_input(n) = deltau_vert_input_save(n) 
            omega_total_input(n) = omega_total_input_save(n) 
            do L = 0, ngreek_moments_input
              do k = 1, 16
                greekmat_total_input(L,n,k) = greekmat_total_input_save(L,n,k)
              enddo
            enddo
          enddo

          !  aerosol layers
          taer = taer * epsfac
          parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
          do n = n6 + 1, nlayers
            aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
            aersca = aerext * waer
            molabs = molext(n) * ( 1.0d0 - molomg(n) )
            molsca = molomg(n) * molext(n)
            totext = molabs + molsca + aerext
            totsca = molsca + aersca
            raywt  = molsca / totsca
            aerwt  = aersca / totsca
            omega  = totsca / totext
            deltau_vert_input(n) = totext
            omega_total_input(n) = omega
            DO K = 1, ngreekmat_entries
              GK = gmask(k); sk = dble(smask(k)) ; ck = cmask(k) ; rk = rmask(k)
              DO L = 0, 2
! Old           MOM = RAYWT*PROBLEM_RAY(ck,L) + AERWT*SK*AERMOMS(CK,L)
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
        endif

!  Copy to optical property type-structure inputs

        VLIDORT_FixIn%Optical%TS_lambertian_albedo    = lambertian_albedo
        VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
        VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
        VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

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

        CALL VLIDORT_WRITE_STATUS ( &
          'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Set some variables

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

!  Write files of Column WF test results
!  =====================================

!  Shorthand indices

      t5 = 5 ; t6 = 6  ; t7 = 7

!   Version 2.7 code superceded
!      c3 = 'FO'
!      if ( .not. VLIDORT_ModIn%MBool%TS_DO_FO_CALC ) then
!        OPEN(OutUnit,file='vlidort_v_test/results_solar_lcs_tester_' &
!             // c1 // '_' // c2 // '.all',status='unknown')
!      else
!        OPEN(OutUnit,file='vlidort_v_test/results_solar_lcs_tester_' &
!             // c1 // '_' // c2 // '_' // c3 //'.all',status='unknown')
!      end if

!  Version 2.8 file

      OutUnit = 36
      OPEN(OutUnit,file='vlidort_v_test/results_solar_lcs_tester_' &
             // trim(c1) // '_' // c2 // '.all',status='unknown')

!  Write Stokes Intensity and Mean values

      call write_stokes_lcs ( &
        OutUnit, 1, 'I', t5, t6, t7, &
        max_tasks, n_geometries, n_szangles,  n_user_levels, &
        user_levels, lambertian_albedo_save, eps, &
        STOKES_BAS, MEANST_BAS, FLUX_BAS, &
        STOKES_PT,  MEANST_PT,  FLUX_PT, &
        COLWF_SAVE, SURFWF_SAVE, &
        MEANST_COLWF_SAVE,  FLUX_COLWF_SAVE, &
        MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )

!  Write Stokes Q and U and Mean values

      if ( NSTOKES .ge. 3 ) then
        call write_stokes_lcs ( &
          OutUnit, 2, 'Q', t5, t6, t7, &
          max_tasks, n_geometries, n_szangles,  n_user_levels, &
          user_levels, lambertian_albedo_save, eps, &
          STOKES_BAS, MEANST_BAS, FLUX_BAS, &
          STOKES_PT,  MEANST_PT,  FLUX_PT, &
          COLWF_SAVE, SURFWF_SAVE, &
          MEANST_COLWF_SAVE,  FLUX_COLWF_SAVE, &
          MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )

        call write_stokes_lcs ( &
          OutUnit, 3, 'U', t5, t6, t7, &
          max_tasks, n_geometries, n_szangles,  n_user_levels, &
          user_levels, lambertian_albedo_save, eps, &
          STOKES_BAS, MEANST_BAS, FLUX_BAS, &
          STOKES_PT,  MEANST_PT,  FLUX_PT, &
          COLWF_SAVE, SURFWF_SAVE, &
          MEANST_COLWF_SAVE,  FLUX_COLWF_SAVE, &
          MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )
      endif

!  Write Stokes V and Mean values

      if ( NSTOKES .eq. 4 ) then
        call write_stokes_lcs ( &
          OutUnit, 4, 'V', t5, t6, t7, &
          max_tasks, n_geometries, n_szangles,  n_user_levels, &
          user_levels, lambertian_albedo_save, eps, &
          STOKES_BAS, MEANST_BAS, FLUX_BAS, &
          STOKES_PT,  MEANST_PT,  FLUX_PT, &
          COLWF_SAVE, SURFWF_SAVE, &
          MEANST_COLWF_SAVE,  FLUX_COLWF_SAVE, &
          MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )
      endif

      close(OutUnit)

!  This is no longer requried for Version 2.8
!  End single-scatter test loop
!      write(*,*)
!      enddo

!  Finish

      write(*,*)
      write(*,*)'Main program finished successfully'

      stop

      end program Solar_lpcs_Tester

!  ****************************************************************************

      subroutine write_stokes_lps ( &
        un, o1, a1, k3, k4, k5, t1, t2, t3, t4, &
        max_tasks, n_geometries, nbeams,  n_user_levels, &
        user_levels, lambertian_albedo_save, eps, &
        STOKES_BAS, MEANST_BAS, FLUX_BAS, &
        STOKES_PT,  MEANST_PT,  FLUX_PT, &
        PROFWF_SAVE, SURFWF_SAVE, &
        MEANST_PROFWF_SAVE, FLUX_PROFWF_SAVE, &
        MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )

      use vlidort_pars_m, Only : max_user_levels, max_geometries, max_szangles, maxstokes,  &
                                 max_directions, max_atmoswfs, max_surfacewfs, maxlayers, &
                                 upidx, dnidx

      implicit none

!  input

      integer, intent(in) ::           max_tasks, n_geometries, nbeams, n_user_levels
      integer, intent(in) ::           k3, k4, k5, t1, t2, t3, t4
      integer, intent(in) ::           un, o1
      character (len=1), intent(in) :: a1
      double precision, intent(in) ::  user_levels ( max_user_levels )
      double precision, intent(in) ::  lambertian_albedo_save, eps

      DOUBLE PRECISION, INTENT(IN) :: STOKES_BAS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: MEANST_BAS &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: FLUX_BAS &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT(IN) :: PROFWF_SAVE &
         ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: SURFWF_SAVE &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT(IN) :: MEANST_PROFWF_SAVE &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: MEANST_SURFWF_SAVE &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT(IN) :: FLUX_PROFWF_SAVE &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: FLUX_SURFWF_SAVE &
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
        'STOKES-'//a1//', 3 PROFILE JACOBIANS, 1 SURFACE JACOBIAN', &
        '================================================='
      write(un,'(a,T32,a,a/)')'Geometry    Level/Output', &
        'Stokes-'//a1// &
        '        Profile AJ1   Profile FD1     Profile AJ2', &
        '   Profile FD2     Profile AJ3   Profile FD3     '// &
        'Surface AJ1   Surface FD1'

      do v = 1, n_geometries
         do n = 1, n_user_levels
           write(un,366)v,'Upwelling @',user_levels(n), &
           STOKES_BAS(n,v,o1,upidx), &
           PROFWF_SAVE(1,k3,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,t2)-STOKES_BAS(n,v,o1,upidx))/eps, &
           PROFWF_SAVE(1,k4,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,t3)-STOKES_BAS(n,v,o1,upidx))/eps, &
           PROFWF_SAVE(2,k5,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,t4)-STOKES_BAS(n,v,o1,upidx))/eps, &
           lambertian_albedo_save*SURFWF_SAVE(1,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,t1)-STOKES_BAS(n,v,o1,upidx))/eps
         enddo

         do n = 1, n_user_levels
           write(un,366)v,'Dnwelling @',user_levels(n), &
           STOKES_BAS(n,v,o1,dnidx), &
           PROFWF_SAVE(1,k3,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,t2)-STOKES_BAS(n,v,o1,dnidx))/eps, &
           PROFWF_SAVE(1,k4,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,t3)-STOKES_BAS(n,v,o1,dnidx))/eps, &
           PROFWF_SAVE(2,k5,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,t4)-STOKES_BAS(n,v,o1,dnidx))/eps, &
           lambertian_albedo_save*SURFWF_SAVE(1,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,t1)-STOKES_BAS(n,v,o1,dnidx))/eps
         enddo
         write(un,*)' '
      enddo

!  Write actinic fluxes

      write(un,'(/T32,a/T32,a/)') &
        'ACT.FLUX-'//a1//', 3 PROFILE JACOBIANS, 1 SURFACE JACOBIAN', &
        '==================================================='
      write(un,'(a,T32,a,a/)')' Sun SZA    Level/Output', &
        'ActFlux-'//a1// &
        '      Profile AJ1   Profile FD1     Profile AJ2', &
        '   Profile FD2     Profile AJ3   Profile FD3    '// &
        'Surface AJ1   Surface FD1'

      do v = 1, nbeams
         do n = 1, n_user_levels
           write(un,366)v,'Upwelling @',user_levels(n), &
           MEANST_BAS(n,v,o1,upidx), &
           MEANST_PROFWF_save(1,k3,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,t2)-MEANST_BAS(n,v,o1,upidx))/eps, &
           MEANST_PROFWF_save(1,k4,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,t3)-MEANST_BAS(n,v,o1,upidx))/eps, &
           MEANST_PROFWF_save(2,k5,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,t4)-MEANST_BAS(n,v,o1,upidx))/eps, &
           lambertian_albedo_save*MEANST_SURFWF_save(1,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,t1)-MEANST_BAS(n,v,o1,upidx))/eps
         enddo

         do n = 1, n_user_levels
           write(un,366)v,'Dnwelling @',user_levels(n), &
           MEANST_BAS(n,v,o1,dnidx), &
           MEANST_PROFWF_save(1,k3,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,t2)-MEANST_BAS(n,v,o1,dnidx))/eps, &
           MEANST_PROFWF_save(1,k4,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,t3)-MEANST_BAS(n,v,o1,dnidx))/eps, &
           MEANST_PROFWF_save(2,k5,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,t4)-MEANST_BAS(n,v,o1,dnidx))/eps, &
           lambertian_albedo_save*MEANST_SURFWF_save(1,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,t1)-MEANST_BAS(n,v,o1,dnidx))/eps
         enddo
         write(un,*)' '
      enddo

!  Write regular fluxes

      write(un,'(/T32,a/T32,a/)') &
        'REG.FLUX-'//a1//', 3 PROFILE JACOBIANS, 1 SURFACE JACOBIAN', &
        '==================================================='
      write(un,'(a,T32,a,a/)')' Sun SZA    Level/Output', &
        'RegFlux-'//a1// &
        '      Profile AJ1   Profile FD1     Profile AJ2', &
        '   Profile FD2     Profile AJ3   Profile FD3    '// &
        'Surface AJ1   Surface FD1'

      do v = 1, nbeams
         do n = 1, n_user_levels
           write(un,366)v,'Upwelling @',user_levels(n), &
           FLUX_BAS(n,v,o1,upidx), &
           flux_PROFWF_SAVE(1,k3,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,t2)-FLUX_BAS(n,v,o1,upidx))/eps, &
           flux_PROFWF_SAVE(1,k4,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,t3)-FLUX_BAS(n,v,o1,upidx))/eps, &
           flux_PROFWF_SAVE(2,k5,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,t4)-FLUX_BAS(n,v,o1,upidx))/eps, &
           lambertian_albedo_save*flux_SURFWF_SAVE(1,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,t1)-FLUX_BAS(n,v,o1,upidx))/eps
         enddo

         do n = 1, n_user_levels
           write(un,366)v,'Dnwelling @',user_levels(n), &
           FLUX_BAS(n,v,o1,dnidx), &
           flux_PROFWF_SAVE(1,k3,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,t2)-FLUX_BAS(n,v,o1,dnidx))/eps, &
           flux_PROFWF_SAVE(1,k4,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,t3)-FLUX_BAS(n,v,o1,dnidx))/eps, &
           flux_PROFWF_SAVE(2,k5,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,t4)-FLUX_BAS(n,v,o1,dnidx))/eps, &
           lambertian_albedo_save*flux_SURFWF_SAVE(1,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,t1)-FLUX_BAS(n,v,o1,dnidx))/eps
         enddo
         write(un,*)' '
      enddo

366   format(i5,T11,a,f6.2,2x,1pe13.6,4(2x,1p2e14.6))

      end subroutine write_stokes_lps

!

      subroutine write_stokes_lcs ( &
        un, o1, a1, t5, t6, t7, &
        max_tasks, n_geometries, nbeams,  n_user_levels, &
        user_levels, lambertian_albedo_save, eps, &
        STOKES_BAS, MEANST_BAS, FLUX_BAS, &
        STOKES_PT,  MEANST_PT,  FLUX_PT, &
        COLWF_SAVE, SURFWF_SAVE, &
        MEANST_COLWF_SAVE,  FLUX_COLWF_SAVE, &
        MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )

      use vlidort_pars_m, Only : max_user_levels, max_geometries, max_szangles, maxstokes, &
                                 max_directions, max_atmoswfs, max_surfacewfs, &
                                 upidx, dnidx

      implicit none

!  input

      integer, intent(in) ::           max_tasks, n_geometries, nbeams, n_user_levels
      integer, intent(in) ::           t5, t6, t7
      integer, intent(in) ::           un, o1
      character (len=1), intent(in) :: a1
      double precision, intent(in) ::  user_levels ( max_user_levels )
      double precision, intent(in) ::  lambertian_albedo_save, eps

      DOUBLE PRECISION, INTENT(IN) :: STOKES_BAS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: MEANST_BAS &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: FLUX_BAS &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT(IN) :: COLWF_SAVE &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: SURFWF_SAVE &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT(IN) :: MEANST_COLWF_SAVE &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: MEANST_SURFWF_SAVE &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT(IN) :: FLUX_COLWF_SAVE &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(IN) :: FLUX_SURFWF_SAVE &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

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
        'STOKES-'//a1//' , 2 COLUMN JACOBIANS, 1 SURFACE JACOBIAN', &
        '================================================='
      write(un,'(a,T32,a/)')'Geometry    Level/Output', &
        'Stokes-'//a1//'        Column AJ1    Column FD1      '// &
        'Column AJ2    Column FD2      Surface AJ1   Surface FD1'

      do v = 1, n_geometries
         do n = 1, n_user_levels
           write(un,367)v,'Upwelling @',user_levels(n), &
           STOKES_BAS(n,v,o1,upidx), &
           COLWF_SAVE(1,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,t6)-STOKES_BAS(n,v,o1,upidx))/eps, &
           COLWF_SAVE(2,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,t7)-STOKES_BAS(n,v,o1,upidx))/eps, &
           lambertian_albedo_save*SURFWF_SAVE(1,n,v,o1,upidx), &
           (STOKES_PT(n,v,o1,upidx,t5)-STOKES_BAS(n,v,o1,upidx))/eps
         enddo

         do n = 1, n_user_levels
           write(un,367)v,'Dnwelling @',user_levels(n), &
           STOKES_BAS(n,v,o1,dnidx), &
           COLWF_SAVE(1,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,t6)-STOKES_BAS(n,v,o1,dnidx))/eps, &
           COLWF_SAVE(2,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,t7)-STOKES_BAS(n,v,o1,dnidx))/eps, &
           lambertian_albedo_save*SURFWF_SAVE(1,n,v,o1,dnidx), &
           (STOKES_PT(n,v,o1,dnidx,t5)-STOKES_BAS(n,v,o1,dnidx))/eps
         enddo
         write(un,*)' '
      enddo

!  Write actinic fluxes

      write(un,'(/T32,a/T32,a/)') &
        'ACT.FLUX-'//a1//', 2 COLUMN JACOBIANS, 1 SURFACE JACOBIAN', &
        '=================================================='
      write(un,'(a,T32,a/)')'  SunZA #   Level/Output', &
        'Act.Flux-'//a1//'       Column AJ1    Column FD1      '// &
        'Column AJ2    Column FD2      Surface AJ1   Surface FD1'

      do v = 1, nbeams
         do n = 1, n_user_levels
           write(un,367)v,'Upwelling @',user_levels(n), &
           MEANST_BAS(n,v,o1,upidx), &
           MEANST_COLWF_save(1,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,t6)-MEANST_BAS(n,v,o1,upidx))/eps, &
           MEANST_COLWF_save(2,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,t7)-MEANST_BAS(n,v,o1,upidx))/eps, &
           lambertian_albedo_save*MEANST_SURFWF_save(1,n,v,o1,upidx), &
           (MEANST_PT(n,v,o1,upidx,t5)-MEANST_BAS(n,v,o1,upidx))/eps
         enddo

         do n = 1, n_user_levels
           write(un,367)v,'Dnwelling @',user_levels(n), &
           MEANST_BAS(n,v,o1,dnidx), &
           MEANST_COLWF_save(1,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,t6)-MEANST_BAS(n,v,o1,dnidx))/eps, &
           MEANST_COLWF_save(2,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,t7)-MEANST_BAS(n,v,o1,dnidx))/eps, &
           lambertian_albedo_save*MEANST_SURFWF_save(1,n,v,o1,dnidx), &
           (MEANST_PT(n,v,o1,dnidx,t5)-MEANST_BAS(n,v,o1,dnidx))/eps
         enddo
         write(un,*)' '
      enddo

!  Write regular fluxes

      write(un,'(/T32,a/T32,a/)') &
        'REG.FLUX-'//a1//', 2 COLUMN JACOBIANS, 1 SURFACE JACOBIAN', &
        '=================================================='
      write(un,'(a,T32,a/)')'  SunZA #   Level/Output', &
        'Reg.Flux-'//a1//'      Column AJ1    Column FD1      '// &
        'Column AJ2    Column FD2      Surface AJ1   Surface FD1'

      do v = 1, nbeams
         do n = 1, n_user_levels
           write(un,367)v,'Upwelling @',user_levels(n), &
           FLUX_BAS(n,v,o1,upidx), &
           flux_COLWF_SAVE(1,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,t6)-FLUX_BAS(n,v,o1,upidx))/eps, &
           flux_COLWF_SAVE(2,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,t7)-FLUX_BAS(n,v,o1,upidx))/eps, &
           lambertian_albedo_save*flux_SURFWF_SAVE(1,n,v,o1,upidx), &
           (FLUX_PT(n,v,o1,upidx,t5)-FLUX_BAS(n,v,o1,upidx))/eps
         enddo

         do n = 1, n_user_levels
           write(un,367)v,'Dnwelling @',user_levels(n), &
           FLUX_BAS(n,v,o1,dnidx), &
           flux_COLWF_SAVE(1,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,t6)-FLUX_BAS(n,v,o1,dnidx))/eps, &
           flux_COLWF_SAVE(2,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,t7)-FLUX_BAS(n,v,o1,dnidx))/eps, &
           lambertian_albedo_save*flux_SURFWF_SAVE(1,n,v,o1,dnidx), &
           (FLUX_PT(n,v,o1,dnidx,t5)-FLUX_BAS(n,v,o1,dnidx))/eps
         enddo
         write(un,*)' '
      enddo

 367  format(i5,T11,a,f6.2,2x,1pe13.6,3(2x,1p2e14.6))

      end subroutine write_stokes_lcs
