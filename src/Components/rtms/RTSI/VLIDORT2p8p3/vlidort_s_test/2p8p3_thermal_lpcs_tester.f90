
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

      program Thermal_lpcs_Tester

!   Upgrade for Version 2.8 driver. Created 9/18/16 from 2.7 Driver
!     R. Spurr. RT Solutions Inc.
        
!  Upgrade for Version 2.8.1, August 2019
!  ---------------------------------------

!  Module files for VLIDORT. Strict usage, Version 2.8 upwards

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_LIN_IO_DEFS_m
      USE VLIDORT_GETPLANCK_m

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
!   Rob 3/22/15. renamed thread --> task

      INTEGER ::          N,N6,L,NDUM,LDUM,V,T,T1,T2,T3,T4,T5,T6,T7
      INTEGER ::          O1,Q,TASK, UTA,K3,K4,K5
      DOUBLE PRECISION :: KD, GAER, WAER, TAER, PARCEL, RAYWT, AERWT
      DOUBLE PRECISION :: AERSCA, AEREXT, MOLSCA, TOTSCA, TOTEXT, OMEGA
      DOUBLE PRECISION :: MOLOMG(MAXLAYERS),MOLEXT(MAXLAYERS),C0
      DOUBLE PRECISION :: AERMOMS(0:MAXMOMENTS_INPUT), EPS, EPSFAC
      DOUBLE PRECISION :: RAYMOMS(0:2,MAXLAYERS), RATIO1, RATIO2, MOLABS

      INTEGER ::          OutUnit

!  Saved optical property inputs

      DOUBLE PRECISION :: OMEGA_TOTAL_INPUT_SAVE ( MAXLAYERS )
      DOUBLE PRECISION :: DELTAU_VERT_INPUT_SAVE ( MAXLAYERS )
      DOUBLE PRECISION :: GREEKMAT_TOTAL_INPUT_SAVE ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      DOUBLE PRECISION :: LAMBERTIAN_ALBEDO_SAVE

!  Thermal emission setups

      LOGICAL ::            THERMFAIL
      CHARACTER (LEN=70) :: THERMALMESSAGE
      DOUBLE PRECISION ::   SURFTEMP, AIRTEMPS(0:MAXLAYERS)
      DOUBLE PRECISION ::   WNUMLO, WNUMHI
      INTEGER ::            NTEMPS, SMALLV

!  VLIDORT standard input preparation
!   Rob 3/22/15. Introduced FO_CALC local variable, renamed FOCORR 3/1/17

      LOGICAL ::          DO_LAMBERTIAN_SURFACE
      LOGICAL ::          DO_THERMAL_EMISSION, DO_SURFACE_EMISSION
      LOGICAL ::          DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY
      LOGICAL ::          DO_FOCORR, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING
      LOGICAL ::          DO_DELTAM_SCALING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING
      LOGICAL ::          DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY
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

      INTEGER ::          NSTREAMS
      INTEGER ::          N_SZANGLES
      DOUBLE PRECISION :: SZANGLES ( MAX_SZANGLES )
      LOGICAL ::          DO_USER_VZANGLES
      INTEGER ::          N_USER_VZANGLES
      DOUBLE PRECISION :: USER_VZANGLES ( MAX_USER_VZANGLES )
      INTEGER ::          N_USER_RELAZMS
      DOUBLE PRECISION :: USER_RELAZMS ( MAX_USER_RELAZMS )
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
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS  )
      DOUBLE PRECISION :: MEANST_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS  )
      DOUBLE PRECISION :: FLUX_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS, MAX_TASKS  )

!  Start of code
!  =============

!  Initialize error file output flag

      OPENFILEFLAG = .false.

!  Define some control variables

!      PART1 = .FALSE.
      PART1 = .TRUE.

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

!  VLIDORT control Read input, abort if failed

      CALL VLIDORT_L_INPUT_MASTER ( &
        'vlidort_s_test/2p8p3_VLIDORT_ReadInput.cfg', & ! Input
        VLIDORT_FixIn,      & ! Outputs
        VLIDORT_ModIn,      & ! Outputs
        VLIDORT_LinFixIn,   & ! Outputs
        VLIDORT_LinModIn,   & ! Outputs
        VLIDORT_InputStatus ) ! Outputs

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( '2p8p3_VLIDORT_ReadInput.log', VLIDORT_InputStatus )

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

!  Initialize some input variables not handled by VLIDORT_L_INPUT_MASTER

      CALL VLIDORT_Sup_Init ( VLIDORT_Sup )
      CALL VLIDORT_LinSup_Init ( VLIDORT_LinSup )

!  Define some driver variables

      NSTOKES              = VLIDORT_FixIn%Cont%TS_NSTOKES
      !NSTREAMS             = VLIDORT_FixIn%Cont%TS_NSTREAMS
      N_USER_LEVELS        = VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS
      N_SZANGLES           = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
      !SZANGLES             = VLIDORT_ModIn%MSunrays%TS_SZANGLES
      !DO_USER_VZANGLES     = VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES
      !N_USER_VZANGLES      = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
      !USER_VZANGLES        = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT
      !N_USER_RELAZMS       = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      !USER_RELAZMS         = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS
      USER_LEVELS          = VLIDORT_ModIn%MUserVal%TS_USER_LEVELS

!  Standard inputs - cross-over w thermal and surface emission
!   Updated Version 2.8

      DO_SOLAR_SOURCES      = .true.

      DO_FOCORR             = .true.
      DO_FOCORR_NADIR       = .true.
      DO_FOCORR_OUTGOING    = .false.
      !DO_FOCORR_NADIR       = .false.
      !DO_FOCORR_OUTGOING    = .true.

      DO_THERMAL_EMISSION   = .true.
      DO_THERMAL_TRANSONLY  = .false.
      DO_DELTAM_SCALING     = .true.
      DO_SOLUTION_SAVING    = .false.
      DO_BVP_TELESCOPING    = .false.
      DO_LAMBERTIAN_SURFACE = .true.
      DO_SURFACE_EMISSION   = .true.

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
      open(45,file='vlidort_s_test/input_atmos.dat',status='old' )
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(1,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0) = 60.0d0
      do n = 1, nlayers
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n), molext(n),molomg(n),kd,kd,kd,kd
      enddo
      close(45)

!  Surface, with saved albedo

      lambertian_albedo_save = 0.05d0

!  Add Aerosols bottom 6 layers, spread evenly

      !do_rayleigh_only = .false.
      n6 = nlayers - 6; ngreek_moments_input = 80
      gaer = 0.8d0 ; waer = 0.95d0 ; taer = 0.5d0 ; aermoms(0) = 1.0d0
      do l = 1, ngreek_moments_input
        aermoms(l) = dble(2*L+1) * gaer ** dble(L)
      enddo

!  Thermal stuff
!  -------------

!   Array of temperatures should be 24 levels

      open(1,file= 'vlidort_s_test/input_temp23.dat', status='old')
      read(1,*)SURFTEMP
      read(1,*)NTEMPS
      if ( NTEMPS .ne. NLAYERS ) then
         write(*,*)'Number of layers in temperature data set not equal to input NLAYERS'
         stop 'Test program 2p8p3_lpcs_thermal_tester stopped'
      endif
      do n = 0, NTEMPS
         read(1,*)ndum, AIRTEMPS(NTEMPS-N)
      enddo
      close(1)

!  Do the BB inputs

      wnumhi = 2500.0d0 + 1.0d0
      wnumlo = 2500.0d0 - 1.0d0
      call get_planckfunction &
         ( wnumlo,wnumhi, SURFTEMP, &
           SURFBB, SMALLV, THERMFAIL, THERMALMESSAGE )
      do n = 0, nlayers
         call get_planckfunction &
         ( wnumlo,wnumhi, AIRTEMPS(N), &
           THERMAL_BB_INPUT(N), SMALLV, THERMFAIL, THERMALMESSAGE )
      enddo
      if ( THERMFAIL ) THEN
         write(*,'(A)')Trim(THERMALMESSAGE)
         stop 'Failure from call to get_planckfunction in program 2p8p3_lpcs_thermal_tester'
      endif

!  Set some VLIDORT variables

      VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT = SURFBB
      VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT(0:NLAYERS) = THERMAL_BB_INPUT(0:NLAYERS)

!     pause'thermal check'

!  Initialise

      deltau_vert_input = 0.0d0
      omega_total_input = 0.0d0
      greekmat_total_input = 0.0d0
      L_deltau_vert_input = 0.0d0
      L_omega_total_input = 0.0d0
      L_greekmat_total_input = 0.0d0

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_lambertian_albedo    = lambertian_albedo_save
      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input    = L_deltau_vert_input
      VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input = L_greekmat_total_input
      VLIDORT_LinFixIn%Optical%TS_L_omega_total_input    = L_omega_total_input

      write(*,*)'Doing PROFILEWF Baseline: '// &
      'Baseline, Intensities + 3 profile Jacobians + 1 Surface Jacobian'

!  Initialize linearized inputs

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

!  Surface

      n_surface_wfs   = 1
      lambertian_albedo = lambertian_albedo_save

!  rayleigh layers

      do n = 1, n6
        deltau_vert_input(n) = molext(n)
        omega_total_input(n) = molomg(n)
        greekmat_total_input(0,n,1) = raymoms(0,n)
        greekmat_total_input(1,n,1) = raymoms(1,n)
        greekmat_total_input(2,n,1) = raymoms(2,n)
        ratio1 = 1.0d0 - molomg(n)
        L_deltau_vert_input(1,n) =   ratio1
        L_omega_total_input(1,n) = - ratio1
      enddo

!  aerosol layers

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
        greekmat_total_input(0,n,1) = 1.0d0
        do L = 1, 2
           greekmat_total_input(L,n,1) = raywt * raymoms(L,n) + aerwt * aermoms(L)
        enddo
        do L = 3, ngreek_moments_input
           greekmat_total_input(L,n,1) = aerwt * aermoms(L)
        enddo
        do L = 1, ngreek_moments_input
           ratio2 = aersca / totsca
           L_greekmat_total_input(2,L,n,1) = ratio2 * ( (aermoms(L)/greekmat_total_input(L,n,1)) - 1.0d0 )
        enddo
      enddo

!  Save inputs
!  -----------

      do n = 1, nlayers
        deltau_vert_input_save(n) = deltau_vert_input(n)
        omega_total_input_save(n) = omega_total_input(n)
        do L = 0, ngreek_moments_input
          greekmat_total_input_save(L,n,1) = greekmat_total_input(L,n,1)
        enddo
      enddo
      lambertian_albedo_save = lambertian_albedo

!  Copy local variables. Upgraded for Version 2.8. 3/1/17

      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION    = DO_THERMAL_EMISSION
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION    = DO_SURFACE_EMISSION
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE  = DO_LAMBERTIAN_SURFACE

      VLIDORT_ModIn%MBool%TS_DO_FOCORR             = DO_FOCORR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR       = DO_FOCORR_NADIR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING    = DO_FOCORR_OUTGOING

      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES      = DO_SOLAR_SOURCES
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING     = DO_DELTAM_SCALING
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING    = DO_SOLUTION_SAVING
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING    = DO_BVP_TELESCOPING
      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY  = DO_THERMAL_TRANSONLY

      VLIDORT_FixIn%Cont%TS_NFINELAYERS            = NFINELAYERS
      VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS       = N_THERMAL_COEFFS
      VLIDORT_FixIn%Cont%TS_NLAYERS                = NLAYERS

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT  = NGREEK_MOMENTS_INPUT

      VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR         = FLUX_FACTOR

      VLIDORT_FixIn%Chapman%TS_height_grid         = HEIGHT_GRID

!  Copy local linearization variables

      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)   = LAYER_VARY_FLAG(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS) = LAYER_VARY_NUMBER(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS           = N_TOTALPROFILE_WFS
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS                = N_SURFACE_WFS

      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = DO_ATMOS_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = DO_SURFACE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = DO_SIMULATION_ONLY

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_FixIn%Optical%TS_lambertian_albedo = lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_emissivity             = 1.0d0 - lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_user_emissivity        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_emissivity(1,:)        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_user_emissivity(1,:)   = 1.0d0 - lambertian_albedo

      VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input    = L_deltau_vert_input
      VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input = L_greekmat_total_input
      VLIDORT_LinFixIn%Optical%TS_L_omega_total_input    = L_omega_total_input

      VLIDORT_LinSup%BRDF%TS_ls_emissivity(1,1,:)      = -1.0d0
      VLIDORT_LinSup%BRDF%TS_ls_user_emissivity(1,1,:) = -1.0d0

!  For special MV-only check

      if ( VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY ) then
         VLIDORT_ModIn%MBool%TS_DO_FOCORR             = .FALSE.
         VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR       = .FALSE.
         VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING    = .FALSE.
         !VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING     = .FALSE.
      endif

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
        '2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

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

!      stop 'baseline profile WFS only'

!  Initialize FD section
!  =====================

!  Fd perturbation

      eps = 1.0d-03
      epsfac = 1.0d0 + eps

!  Initialize some inputs

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

!  Task 1, FD Perturbation on Surface Albedo
!  =========================================

      deltau_vert_input    = 0.0d0
      omega_total_input    = 0.0d0
      greekmat_total_input = 0.0d0

!  Task

      task = 1; t = task; t1 = t
      write(*,*)'Doing task 1: '// 'Finite Difference, perturb Lambertian albedo'

!  Surface - Perturb albedo

      lambertian_albedo = lambertian_albedo_save * epsfac

!  Atmosphere - Copy saved values

      do n = 1, nlayers
        deltau_vert_input(n) = deltau_vert_input_save(n)
        omega_total_input(n) = omega_total_input_save(n)
        do L = 0, ngreek_moments_input
          greekmat_total_input(L,n,1) = greekmat_total_input_save(L,n,1)
        enddo
      enddo

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_FixIn%Optical%TS_lambertian_albedo = lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_emissivity             = 1.0d0 - lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_user_emissivity        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_emissivity             = zero
      VLIDORT_Sup%BRDF%TS_user_emissivity        = zero
      VLIDORT_Sup%BRDF%TS_emissivity(1,:)        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_user_emissivity(1,:)   = 1.0d0 - lambertian_albedo

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
        '2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

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

!      pause' only interested in surf WF albedo'

!  Task 2, FD Perturbation on trace gas optical depth in layer 1
!  =============================================================

      deltau_vert_input    = 0.0d0
      omega_total_input    = 0.0d0
      greekmat_total_input = 0.0d0

!  Task

      k3 = 1
      task = 2; t = task; t2 = t
      write(*,*)'Doing task 2: '// 'Finite Difference, perturb molecular absorption in Layer 1'

!  Surface

      lambertian_albedo = lambertian_albedo_save

!  rayleigh layers

      do n = 1, n6
        molabs = molext(n) * ( 1.0d0 - molomg(n) )
        if ( n.eq.k3 ) molabs = molabs * epsfac
        molsca = molomg(n) * molext(n)
        deltau_vert_input(n) = molabs + molsca
        omega_total_input(n) = molsca / ( molabs + molsca )
        greekmat_total_input(0:2,n,1) = raymoms(0:2,n)
      enddo

!  aerosol layers

      do n = n6 + 1, nlayers
        aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
        aersca = aerext * waer
        molabs = molext(n) * ( 1.0d0 - molomg(n) )
        if ( n.eq.k3 ) molabs = molabs * epsfac
        molsca = molomg(n) * molext(n)
        totext = molabs + molsca + aerext
        totsca = molsca + aersca
        raywt  = molsca / totsca
        aerwt  = aersca / totsca
        omega  = totsca / totext
        deltau_vert_input(n) = totext
        omega_total_input(n) = omega
        greekmat_total_input(0,n,1) = 1.0d0
        do L = 1, 2
           greekmat_total_input(L,n,1) = raywt * raymoms(L,n) + aerwt * aermoms(L)
        enddo
        do L = 3, ngreek_moments_input
           greekmat_total_input(L,n,1) = aerwt * aermoms(L)
        enddo
      enddo

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_FixIn%Optical%TS_lambertian_albedo = lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_emissivity             = 1.0d0 - lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_user_emissivity        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_emissivity             = zero
      VLIDORT_Sup%BRDF%TS_user_emissivity        = zero
      VLIDORT_Sup%BRDF%TS_emissivity(1,:)        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_user_emissivity(1,:)   = 1.0d0 - lambertian_albedo

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
        '2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

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

!  Task 3, FD Perturbation on trace gas optical depth in layer 21
!  ==============================================================

      deltau_vert_input    = 0.0d0
      omega_total_input    = 0.0d0
      greekmat_total_input = 0.0d0

!  Task

      k4 = 21
      task = 3; t = task; t3 = t
      write(*,*)'Doing task 3: '// 'Finite Difference, perturb molecular absorption in Layer 21'

!  Surface

      lambertian_albedo = lambertian_albedo_save

!  rayleigh layers

      do n = 1, n6
        molabs = molext(n) * ( 1.0d0 - molomg(n) )
        if ( n.eq.k4 ) molabs = molabs * epsfac
        molsca = molomg(n) * molext(n)
        deltau_vert_input(n) = molabs + molsca
        omega_total_input(n) = molsca / ( molabs + molsca )
        greekmat_total_input(0:2,n,1) = raymoms(0:2,n)
      enddo

!  aerosol layers

      do n = n6 + 1, nlayers
        aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
        aersca = aerext * waer
        molabs = molext(n) * ( 1.0d0 - molomg(n) )
        if ( n.eq.k4 ) molabs = molabs * epsfac
        molsca = molomg(n) * molext(n)
        totext = molabs + molsca + aerext
        totsca = molsca + aersca
        raywt  = molsca / totsca
        aerwt  = aersca / totsca
        omega  = totsca / totext
        deltau_vert_input(n) = totext
        omega_total_input(n) = omega
        greekmat_total_input(0,n,1) = 1.0d0
        do L = 1, 2
           greekmat_total_input(L,n,1) = raywt * raymoms(L,n) + aerwt * aermoms(L)
        enddo
        do L = 3, ngreek_moments_input
           greekmat_total_input(L,n,1) = aerwt * aermoms(L)
        enddo
      enddo

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_FixIn%Optical%TS_lambertian_albedo = lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_emissivity             = 1.0d0 - lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_user_emissivity        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_emissivity             = zero
      VLIDORT_Sup%BRDF%TS_user_emissivity        = zero
      VLIDORT_Sup%BRDF%TS_emissivity(1,:)        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_user_emissivity(1,:)   = 1.0d0 - lambertian_albedo

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
        '2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

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

!  Task 4, FD Perturbation on aerosol optical depth in layer 23
!  ============================================================

      deltau_vert_input    = 0.0d0
      omega_total_input    = 0.0d0
      greekmat_total_input = 0.0d0

!  Task

      k5 = 23
      task = 4; t = task; t4 = t
      write(*,*)'Doing task 4: '// &
        'Finite Difference, perturb aerosol optical depth in layer 23'

!  Surface

      lambertian_albedo = lambertian_albedo_save

!  rayleigh layers

      do n = 1, n6
        molabs = molext(n) * ( 1.0d0 - molomg(n) )
        molsca = molomg(n) * molext(n)
        deltau_vert_input(n) = molabs + molsca
        omega_total_input(n) = molsca / ( molabs + molsca )
        greekmat_total_input(0:2,n,1) = raymoms(0:2,n)
      enddo

!  aerosol layers

      parcel = taer / ( height_grid(n6) - height_grid(nlayers) )
      do n = n6 + 1, nlayers
        aerext = Parcel * ( height_grid(n-1) - height_grid(n) )
        if ( k5 .eq. n ) aerext = aerext * epsfac
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
        greekmat_total_input(0,n,1) = 1.0d0
        do L = 1, 2
           greekmat_total_input(L,n,1) = raywt * raymoms(L,n) + aerwt * aermoms(L)
        enddo
        do L = 3, ngreek_moments_input
           greekmat_total_input(L,n,1) = aerwt * aermoms(L)
        enddo
      enddo

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_FixIn%Optical%TS_lambertian_albedo = lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_emissivity             = 1.0d0 - lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_user_emissivity        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_emissivity             = zero
      VLIDORT_Sup%BRDF%TS_user_emissivity        = zero
      VLIDORT_Sup%BRDF%TS_emissivity(1,:)        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_user_emissivity(1,:)   = 1.0d0 - lambertian_albedo

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
        '2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

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

!  Write files of Profile WF test results
!  ======================================

!  Shorthand indices

      k3 = 1 ; k4 = 21 ; k5 = 23
      t1 = 1 ; t2 = 2  ; t3 = 3 ; t4 = 4

!  Set some variables

      do_additional_mvout = VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT
      do_mvout_only       = VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

!  Title of output file

!  Superceded code from Version 2.7
!      c3 = 'FO'
!      if ( .not. VLIDORT_ModIn%MBool%TS_DO_FO_CALC ) then
!        OPEN(36,file = 'vlidort_s_test/results_thermal_lps_tester.all', status = 'unknown')
!      else
!        OPEN(36,file = 'vlidort_s_test/results_thermal_lps_tester' // '_' // c3 // '.all', status = 'unknown')
!      end if

 !  Version 2.8 only the FO calculation

      OutUnit = 36
      OPEN(OutUnit,file = 'vlidort_s_test/results_thermal_lps_tester.all', status = 'unknown')

!  Write Stokes Intensity and Mean values

      call write_stokes_lps_thermal ( &
        OutUnit, 1, 'I', k3, k4, k5, t1, t2, t3, t4, &
        max_tasks, n_geometries, n_szangles,  n_user_levels, &
        user_levels, lambertian_albedo_save, eps, &
        do_additional_mvout, do_mvout_only, &
        STOKES_BAS, MEANST_BAS, FLUX_BAS, &
        STOKES_PT,  MEANST_PT,  FLUX_PT, &
        PROFWF_SAVE, SURFWF_SAVE, &
        MEANST_PROFWF_SAVE, FLUX_PROFWF_SAVE, &
        MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )

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

!  Get the preprepared atmosphere

      height_grid = 0.0d0
      nlayers = VLIDORT_FixIn%Cont%TS_nlayers
      open(45,file='vlidort_s_test/input_atmos.dat',status='old' )
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(1,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0) = 60.0d0
      do n = 1, nlayers
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n), molext(n),molomg(n),kd,kd,kd,kd
      enddo
      close(45)

!  Surface, with saved albedo

      lambertian_albedo_save = 0.05d0

!  Add Aerosols bottom 6 layers, spread evenly

      !do_rayleigh_only = .false.
      n6 = nlayers - 6; ngreek_moments_input = 80
      gaer = 0.8d0 ; waer = 0.95d0 ; taer = 0.5d0 ; aermoms(0) = 1.0d0
      do l = 1, ngreek_moments_input
        aermoms(l) = dble(2*L+1) * gaer ** dble(L)
      enddo

!  Total column

      c0 = 0.0d0
      do n = 1, nlayers
        molabs = molext(n) * ( 1.0d0 - molomg(n) )
        c0 = c0 + molabs
      enddo

!  Initialize

      deltau_vert_input      = 0.0d0
      omega_total_input      = 0.0d0
      greekmat_total_input   = 0.0d0
      L_deltau_vert_input    = 0.0d0
      L_omega_total_input    = 0.0d0
      L_greekmat_total_input = 0.0d0

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_lambertian_albedo    = lambertian_albedo_save
      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input    = L_deltau_vert_input
      VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input = L_greekmat_total_input
      VLIDORT_LinFixIn%Optical%TS_L_omega_total_input    = L_omega_total_input

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

!  Surface

      n_surface_wfs   = 1
      lambertian_albedo = lambertian_albedo_save

!  rayleigh layers

      do n = 1, n6
        deltau_vert_input(n) = molext(n)
        omega_total_input(n) = molomg(n)
        greekmat_total_input(0,n,1) = raymoms(0,n)
        greekmat_total_input(1,n,1) = raymoms(1,n)
        greekmat_total_input(2,n,1) = raymoms(2,n)
        ratio1 = 1.0d0 - molomg(n)
        L_deltau_vert_input(1,n) =   ratio1
        L_omega_total_input(1,n) = - ratio1
      enddo

!  aerosol layers

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
        greekmat_total_input(0,n,1) = 1.0d0
        do L = 1, 2
           greekmat_total_input(L,n,1) = raywt * raymoms(L,n) + aerwt * aermoms(L)
        enddo
        do L = 3, ngreek_moments_input
           greekmat_total_input(L,n,1) = aerwt * aermoms(L)
        enddo
        do L = 1, ngreek_moments_input
           ratio2 = aersca / totsca
           L_greekmat_total_input(2,L,n,1) = ratio2 * &
             ( (aermoms(L)/greekmat_total_input(L,n,1)) - 1.0d0 )
        enddo
      enddo

!  Save inputs
!  -----------

      do n = 1, nlayers
        deltau_vert_input_save(n) = deltau_vert_input(n)
        omega_total_input_save(n) = omega_total_input(n)
        do L = 0, ngreek_moments_input
          greekmat_total_input_save(L,n,1) = greekmat_total_input(L,n,1)
        enddo
      enddo
      lambertian_albedo_save = lambertian_albedo

!  Copy local control integers

      VLIDORT_FixIn%Cont%TS_NLAYERS               = NLAYERS
      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = NGREEK_MOMENTS_INPUT
      VLIDORT_FixIn%Chapman%TS_height_grid        = HEIGHT_GRID

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

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_FixIn%Optical%TS_lambertian_albedo = lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_emissivity                = 1.0d0 - lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_user_emissivity           = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_emissivity                = zero
      VLIDORT_Sup%BRDF%TS_user_emissivity           = zero
      VLIDORT_Sup%BRDF%TS_emissivity(1,:)           = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_user_emissivity(1,:)      = 1.0d0 - lambertian_albedo

      VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input    = L_deltau_vert_input
      VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input = L_greekmat_total_input
      VLIDORT_LinFixIn%Optical%TS_L_omega_total_input    = L_omega_total_input

      VLIDORT_LinSup%BRDF%TS_ls_emissivity(1,1,:)      = -1.0d0
      VLIDORT_LinSup%BRDF%TS_ls_user_emissivity(1,1,:) = -1.0d0

!  For special MV-only check

      if ( VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY ) then
         VLIDORT_ModIn%MBool%TS_DO_FOCORR             = .FALSE.
         VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR       = .FALSE.
         VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING    = .FALSE.
         !VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING     = .FALSE.
      endif

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
        '2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

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

      eps    = 1.0d-03
      epsfac = 1.0d0 + eps

!  Initialize some inputs

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

!  Task 5, FD Perturbation on Surface Albedo
!  =========================================

      deltau_vert_input = 0.0d0
      omega_total_input = 0.0d0
      greekmat_total_input = 0.0d0

!  Task

      task = 5; t = task; t5 = t
      write(*,*)'Doing task 5: '// 'Finite Difference, perturb Lambertian albedo'

!  Surface - Perturb albedo

      lambertian_albedo = lambertian_albedo_save * epsfac

!  Atmosphere - Copy saved values

      do n = 1, nlayers
        deltau_vert_input(n) = deltau_vert_input_save(n)
        omega_total_input(n) = omega_total_input_save(n)
        do L = 0, ngreek_moments_input
          greekmat_total_input(L,n,1) = greekmat_total_input_save(L,n,1)
        enddo
      enddo

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_FixIn%Optical%TS_lambertian_albedo = lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_emissivity             = 1.0d0 - lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_user_emissivity        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_emissivity             = zero
      VLIDORT_Sup%BRDF%TS_user_emissivity        = zero
      VLIDORT_Sup%BRDF%TS_emissivity(1,:)        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_user_emissivity(1,:)   = 1.0d0 - lambertian_albedo

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
        '2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

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

!  Task 6, FD Perturbation on Total trace gas optical depth C0
!  =============================================================

      deltau_vert_input = 0.0d0
      omega_total_input = 0.0d0
      greekmat_total_input = 0.0d0

!  Task

      task = 6; t = task; t6 = t
      write(*,*)'Doing task 6: '// 'Finite Difference, perturb total molecular absoprtion'

!  Surface

      lambertian_albedo = lambertian_albedo_save

!  rayleigh layers

      do n = 1, n6
        molabs = molext(n) * ( 1.0d0 - molomg(n) ) * epsfac
        molsca = molomg(n) * molext(n)
        deltau_vert_input(n) = molabs + molsca
        omega_total_input(n) = molsca / ( molabs + molsca )
        greekmat_total_input(0:2,n,1) = raymoms(0:2,n)
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
        greekmat_total_input(0,n,1) = 1.0d0
        do L = 1, 2
           greekmat_total_input(L,n,1) = raywt * raymoms(L,n) + aerwt * aermoms(L)
        enddo
        do L = 3, ngreek_moments_input
           greekmat_total_input(L,n,1) = aerwt * aermoms(L)
        enddo
      enddo

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_FixIn%Optical%TS_lambertian_albedo = lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_emissivity             = 1.0d0 - lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_user_emissivity        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_emissivity             = zero
      VLIDORT_Sup%BRDF%TS_user_emissivity        = zero
      VLIDORT_Sup%BRDF%TS_emissivity(1,:)        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_user_emissivity(1,:)   = 1.0d0 - lambertian_albedo

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
        '2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

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

!  Task 7, FD Perturbation on Total aerosol optical depth TAU
!  =============================================================

      deltau_vert_input = 0.0d0
      omega_total_input = 0.0d0
      greekmat_total_input = 0.0d0

!  Task

      task = 7; t = task; t7 = t
      write(*,*)'Doing task 7: '// 'Finite Difference, perturb total aerosol'

!  Surface

      lambertian_albedo = lambertian_albedo_save

!  rayleigh layers

      do n = 1, n6
        molabs = molext(n) * ( 1.0d0 - molomg(n) )
        molsca = molomg(n) * molext(n)
        deltau_vert_input(n) = molabs + molsca
        omega_total_input(n) = molsca / ( molabs + molsca )
        greekmat_total_input(0:2,n,1) = raymoms(0:2,n)
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
        greekmat_total_input(0,n,1) = 1.0d0
        do L = 1, 2
           greekmat_total_input(L,n,1) = raywt * raymoms(L,n) + aerwt * aermoms(L)
        enddo
        do L = 3, ngreek_moments_input
           greekmat_total_input(L,n,1) = aerwt * aermoms(L)
        enddo
      enddo

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_FixIn%Optical%TS_lambertian_albedo = lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_emissivity             = 1.0d0 - lambertian_albedo
      !VLIDORT_Sup%BRDF%TS_user_emissivity        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_emissivity             = zero
      VLIDORT_Sup%BRDF%TS_user_emissivity        = zero
      VLIDORT_Sup%BRDF%TS_emissivity(1,:)        = 1.0d0 - lambertian_albedo
      VLIDORT_Sup%BRDF%TS_user_emissivity(1,:)   = 1.0d0 - lambertian_albedo

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
        '2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

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

!  Write files of Column WF test results
!  =====================================

!  Shorthand indices

      t5 = 5 ; t6 = 6  ; t7 = 7

!  Set some variables

      do_additional_mvout = VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT
      do_mvout_only       = VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

!  Title of output file

!   Version 2.7 code superceded
!      c3 = 'FO'
!      if ( .not. VLIDORT_ModIn%MBool%TS_DO_FO_CALC ) then
!        OPEN(36,file = 'vlidort_s_test/results_thermal_lcs_tester.all', status = 'unknown')
!      else
!        OPEN(36,file = 'vlidort_s_test/results_thermal_lcs_tester' // '_' // c3 // '.all', status = 'unknown')
!      end if

!  Version 2.8 file

      OutUnit = 36
      OPEN(36,file = 'vlidort_s_test/results_thermal_lcs_tester.all', status = 'unknown')

!  Write Stokes Intensity and Mean values

      call write_stokes_lcs_thermal ( &
        OutUnit, 1, 'I', t5, t6, t7, &
        max_tasks, n_geometries, n_szangles,  n_user_levels, &
        user_levels, lambertian_albedo_save, eps, &
        do_additional_mvout, do_mvout_only, &
        STOKES_BAS, MEANST_BAS, FLUX_BAS, &
        STOKES_PT,  MEANST_PT,  FLUX_PT, &
        COLWF_SAVE, SURFWF_SAVE, &
        MEANST_COLWF_SAVE,  FLUX_COLWF_SAVE, &
        MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )

      close(OutUnit)

!  This is no longer requried for Version 2.8
!  End single-scatter test loop
!      write(*,*)
!      enddo

!  Finish

      write(*,*)
      write(*,*)'Main program finished successfully'

      stop

      end program Thermal_lpcs_Tester

!  ****************************************************************************

      subroutine write_stokes_lps_thermal ( &
        un, o1, a1, k3, k4, k5, t1, t2, t3, t4, &
        max_tasks, n_geometries, nbeams,  n_user_levels, &
        user_levels, lambertian_albedo_save, eps, &
        do_additional_mvout, do_mvout_only, &
        STOKES_BAS, MEANST_BAS, FLUX_BAS, &
        STOKES_PT,  MEANST_PT,  FLUX_PT, &
        PROFWF_SAVE,      SURFWF_SAVE, &
        MEANST_PROFWF_SAVE, FLUX_PROFWF_SAVE, &
        MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )

      use vlidort_pars_m, Only : max_user_levels, max_geometries, max_szangles, maxstokes, &
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
      logical, intent(in) ::           do_additional_mvout, do_mvout_only

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

      if ( .not. do_mvout_only ) then
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
      endif

      if ( do_additional_mvout .or. do_mvout_only ) then

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
             MEANST_BAS(n,v,o1,1), &
             MEANST_PROFWF_SAVE(1,k3,n,v,o1,1), &
             (MEANST_PT(n,v,o1,1,t2)-MEANST_BAS(n,v,o1,1))/eps, &
             MEANST_PROFWF_SAVE(1,k4,n,v,o1,1), &
             (MEANST_PT(n,v,o1,1,t3)-MEANST_BAS(n,v,o1,1))/eps, &
             MEANST_PROFWF_SAVE(2,k5,n,v,o1,1), &
             (MEANST_PT(n,v,o1,1,t4)-MEANST_BAS(n,v,o1,1))/eps, &
             lambertian_albedo_save*MEANST_SURFWF_SAVE(1,n,v,o1,1), &
             (MEANST_PT(n,v,o1,1,t1)-MEANST_BAS(n,v,o1,1))/eps
           enddo

           do n = 1, n_user_levels
             write(un,366)v,'Dnwelling @',user_levels(n), &
             MEANST_BAS(n,v,o1,2), &
             MEANST_PROFWF_SAVE(1,k3,n,v,o1,2), &
             (MEANST_PT(n,v,o1,2,t2)-MEANST_BAS(n,v,o1,2))/eps, &
             MEANST_PROFWF_SAVE(1,k4,n,v,o1,2), &
             (MEANST_PT(n,v,o1,2,t3)-MEANST_BAS(n,v,o1,2))/eps, &
             MEANST_PROFWF_SAVE(2,k5,n,v,o1,2), &
             (MEANST_PT(n,v,o1,2,t4)-MEANST_BAS(n,v,o1,2))/eps, &
             lambertian_albedo_save*MEANST_SURFWF_SAVE(1,n,v,o1,2), &
             (MEANST_PT(n,v,o1,2,t1)-MEANST_BAS(n,v,o1,2))/eps
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
             FLUX_BAS(n,v,o1,1), &
             flux_PROFWF_SAVE(1,k3,n,v,o1,1), &
             (FLUX_PT(n,v,o1,1,t2)-FLUX_BAS(n,v,o1,1))/eps, &
             flux_PROFWF_SAVE(1,k4,n,v,o1,1), &
             (FLUX_PT(n,v,o1,1,t3)-FLUX_BAS(n,v,o1,1))/eps, &
             flux_PROFWF_SAVE(2,k5,n,v,o1,1), &
             (FLUX_PT(n,v,o1,1,t4)-FLUX_BAS(n,v,o1,1))/eps, &
             lambertian_albedo_save*flux_SURFWF_SAVE(1,n,v,o1,1), &
             (FLUX_PT(n,v,o1,1,t1)-FLUX_BAS(n,v,o1,1))/eps
           enddo

           do n = 1, n_user_levels
             write(un,366)v,'Dnwelling @',user_levels(n), &
             FLUX_BAS(n,v,o1,2), &
             flux_PROFWF_SAVE(1,k3,n,v,o1,2), &
             (FLUX_PT(n,v,o1,2,t2)-FLUX_BAS(n,v,o1,2))/eps, &
             flux_PROFWF_SAVE(1,k4,n,v,o1,2), &
             (FLUX_PT(n,v,o1,2,t3)-FLUX_BAS(n,v,o1,2))/eps, &
             flux_PROFWF_SAVE(2,k5,n,v,o1,2), &
             (FLUX_PT(n,v,o1,2,t4)-FLUX_BAS(n,v,o1,2))/eps, &
             lambertian_albedo_save*flux_SURFWF_SAVE(1,n,v,o1,2), &
             (FLUX_PT(n,v,o1,2,t1)-FLUX_BAS(n,v,o1,2))/eps
           enddo
           write(un,*)' '
        enddo
      endif

366   format(i5,T11,a,f6.2,2x,1pe13.6,4(2x,1p2e14.6))

      end subroutine write_stokes_lps_thermal

!

      subroutine write_stokes_lcs_thermal ( &
        un, o1, a1, t5, t6, t7, &
        max_tasks, n_geometries, nbeams,  n_user_levels, &
        user_levels, lambertian_albedo_save, eps, &
        do_additional_mvout, do_mvout_only, &
        STOKES_BAS, MEANST_BAS, FLUX_BAS, &
        STOKES_PT,  MEANST_PT,  FLUX_PT, &
        COLWF_SAVE,       SURFWF_SAVE, &
        MEANST_COLWF_SAVE,  FLUX_COLWF_SAVE, &
        MEANST_SURFWF_SAVE, FLUX_SURFWF_SAVE )

      use vlidort_pars_m, Only : max_user_levels, max_geometries, max_szangles, maxstokes, &
                                 max_directions, max_atmoswfs, max_surfacewfs, upidx, dnidx

      implicit none

!  input

      integer, intent(in) ::           max_tasks, n_geometries, nbeams, n_user_levels
      integer, intent(in) ::           t5, t6, t7
      integer, intent(in) ::           un, o1
      character (len=1), intent(in) :: a1
      double precision, intent(in) ::  user_levels ( max_user_levels )
      double precision, intent(in) ::  lambertian_albedo_save, eps
      logical, intent(in) ::           do_additional_mvout, do_mvout_only

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
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
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

      if ( .not. do_mvout_only ) then
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
      endif

      if ( do_additional_mvout .or. do_mvout_only ) then

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
             MEANST_BAS(n,v,o1,1), &
             MEANST_COLWF_SAVE(1,n,v,o1,1), &
             (MEANST_PT(n,v,o1,1,t6)-MEANST_BAS(n,v,o1,1))/eps, &
             MEANST_COLWF_SAVE(2,n,v,o1,1), &
             (MEANST_PT(n,v,o1,1,t7)-MEANST_BAS(n,v,o1,1))/eps, &
             lambertian_albedo_save*MEANST_SURFWF_SAVE(1,n,v,o1,1), &
             (MEANST_PT(n,v,o1,1,t5)-MEANST_BAS(n,v,o1,1))/eps
           enddo

           do n = 1, n_user_levels
             write(un,367)v,'Dnwelling @',user_levels(n), &
             MEANST_BAS(n,v,o1,2), &
             MEANST_COLWF_SAVE(1,n,v,o1,2), &
             (MEANST_PT(n,v,o1,2,t6)-MEANST_BAS(n,v,o1,2))/eps, &
             MEANST_COLWF_SAVE(2,n,v,o1,2), &
             (MEANST_PT(n,v,o1,2,t7)-MEANST_BAS(n,v,o1,2))/eps, &
             lambertian_albedo_save*MEANST_SURFWF_SAVE(1,n,v,o1,2), &
             (MEANST_PT(n,v,o1,2,t5)-MEANST_BAS(n,v,o1,2))/eps
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
             FLUX_BAS(n,v,o1,1), &
             flux_COLWF_SAVE(1,n,v,o1,1), &
             (FLUX_PT(n,v,o1,1,t6)-FLUX_BAS(n,v,o1,1))/eps, &
             flux_COLWF_SAVE(2,n,v,o1,1), &
             (FLUX_PT(n,v,o1,1,t7)-FLUX_BAS(n,v,o1,1))/eps, &
             lambertian_albedo_save*flux_SURFWF_SAVE(1,n,v,o1,1), &
             (FLUX_PT(n,v,o1,1,t5)-FLUX_BAS(n,v,o1,1))/eps
           enddo

           do n = 1, n_user_levels
             write(un,367)v,'Dnwelling @',user_levels(n), &
             FLUX_BAS(n,v,o1,2), &
             flux_COLWF_SAVE(1,n,v,o1,2), &
             (FLUX_PT(n,v,o1,2,t6)-FLUX_BAS(n,v,o1,2))/eps, &
             flux_COLWF_SAVE(2,n,v,o1,2), &
             (FLUX_PT(n,v,o1,2,t7)-FLUX_BAS(n,v,o1,2))/eps, &
             lambertian_albedo_save*flux_SURFWF_SAVE(1,n,v,o1,2), &
             (FLUX_PT(n,v,o1,2,t5)-FLUX_BAS(n,v,o1,2))/eps
           enddo
           write(un,*)' '
        enddo
      endif

 367  format(i5,T11,a,f6.2,2x,1pe13.6,3(2x,1p2e14.6))

      end subroutine write_stokes_lcs_thermal
