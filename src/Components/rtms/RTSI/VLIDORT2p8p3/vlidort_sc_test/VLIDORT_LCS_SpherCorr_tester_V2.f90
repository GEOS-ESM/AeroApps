
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

      PROGRAM VLIDORT_LCS_SPHERCORR_TESTER_V1

!  HERE IS THE HISTORY
!  ===================

!  2-point correction history --
!    V1: 11/20/19. Originally Coded 20-29 November 2019 for VLIDORT Version 2.8.1
!    V2: 12/18/19. Extension to include BOA downwelling situation as well as TOA upwelling.
!    V3: 03/01/20. Renamed, along with addition of 3pt Correction

!  3-point correction history --
!    V1: 03/01/20. New 3pt Correction

!  Multi-point correction history --
!    V1: 03/01/20. New  Correction

!  Reprise November 2020.
!    Bring together tests for all corrections in one driver

!  1/31/21. Developed for Official VLIDORT Version 2.8.3 Example
!     R. Spurr. RT Solutions Inc.
        
!  THIS IS THE VLIDORT TESTER
!  ==========================

!  Module files for VLIDORT.

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_LIN_IO_DEFS_m

      USE VLIDORT_AUX_m,      Only : VLIDORT_READ_ERROR, VLIDORT_WRITE_STATUS
      USE VLIDORT_INPUTS_m,   Only : VLIDORT_Sup_Init
      USE VLIDORT_L_INPUTS_m, Only : VLIDORT_L_INPUT_MASTER, VLIDORT_LinSup_Init

!  MS 2-, 3- and Multi-point sphericity correction modules

      USE VLIDORT_LCS_SPHERCORR_ROUTINES_m

!  Implicit none

      IMPLICIT NONE

!  VLIDORT file inputs status structure

      TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus

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

!  Flags for opening error output file, debug input

      LOGICAL ::  OPENFILEFLAG, DO_DEBUG_INPUT, Do_Aerosols, Do_Albedo

!  Help variables

      INTEGER   :: NGREEKMAT_ENTRIES, L, LDUM, LCHECK, N_AERMOMS
      INTEGER   :: K, CK, GK, RK, CMASK(8), GMASK(8), SMASK(8), RMASK(8)
      INTEGER   :: N,N6,NDUM
      REAL(FPK) :: KD, WAER, TAER, PARCEL, RAYWT, AERWT
      REAL(FPK) :: AERSCA, AEREXT, MOLSCA, MOLABS, TOTSCA, TOTEXT
      REAL(FPK) :: MOLOMG ( MAXLAYERS ), MOLEXT ( MAXLAYERS )
      REAL(FPK) :: AMOMS(6), AERMOMS(6,0:MAXMOMENTS_INPUT)
      REAL(FPK) :: RAYMOMS(0:2,MAXLAYERS), PROBLEM_RAY(6,0:2), BETA2, DEPOL
      REAL(FPK) :: CUTOFF, MOM, RATIO1, RATIO2, OMEGA, SK

!  Control linearization

      LOGICAL         ::  DO_COLUMN_LINEARIZATION
      LOGICAL         ::  DO_SURFACE_LINEARIZATION
      LOGICAL         ::  LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER         ::  LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER         ::  N_TOTALCOLUMN_WFS
      INTEGER         ::  N_SURFACE_WFS

!  Height proxy, Optical proxies

      INTEGER   :: NLAYERS, NSTOKES, NS, NFINE, NGREEK_MOMENTS_INPUT
      REAL(FPK) :: HEIGHT_GRID ( 0:MAXLAYERS )
      REAL(FPK) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      REAL(FPK) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      REAL(FPK) :: greekmat_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, 16 )

!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (F-matrix  variation) input

      REAL(fpk) :: L_OMEGA_TOTAL_INPUT(MAX_ATMOSWFS,MAXLAYERS)
      REAL(fpk) :: L_DELTAU_VERT_INPUT(MAX_ATMOSWFS,MAXLAYERS)
      REAL(fpk) :: L_greekmat_TOTAL_INPUT (MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,MAXLAYERS,16)

!  Main MS sphericity correction module output with explanations
!    Applies equally to the Satellite view (upwelling) or the ground view (downwelling)

      INTEGER, PARAMETER :: MAXGEO = 150

      REAL(FPK) :: FO_STOKES_BOAGEOM   ( MAXSTOKES, MAXGEO ) ! 1. FO (SS/DB) Stokes vector for BOA geometry
      REAL(FPK) :: MS_STOKES_BOAGEOM   ( MAXSTOKES, MAXGEO ) ! 2. MS Stokes vector for BOA geometry
      REAL(FPK) :: MS_STOKES_SPHERCORR ( MAXSTOKES, MAXGEO ) ! 3. MS Stokes vector with MS2pt corr
      REAL(FPK) :: STOKES_BOAGEOM      ( MAXSTOKES, MAXGEO ) ! 4. Full Stokes vector for BOA geometry      4 = 1 + 2
      REAL(FPK) :: STOKES_SPHERCORR    ( MAXSTOKES, MAXGEO ) ! 5. Full Stokes vector with Sphericity corr  5 = 1 + 3

!  Column Jacobian output with explanations

      REAL(fpk) :: FO_COLUMNWF_BOAGEOM   ( MAX_ATMOSWFS, MAXSTOKES, MAXGEO ) ! 1. FO (SS/DB) LC Jac for BOA geometry
      REAL(fpk) :: MS_COLUMNWF_BOAGEOM   ( MAX_ATMOSWFS, MAXSTOKES, MAXGEO ) ! 2. MS LC Jac for BOA geometry
      REAL(fpk) :: MS_COLUMNWF_SPHERCORR ( MAX_ATMOSWFS, MAXSTOKES, MAXGEO ) ! 3. MS LC Jac with 2-point corr
      REAL(fpk) :: COLUMNWF_BOAGEOM      ( MAX_ATMOSWFS, MAXSTOKES, MAXGEO ) ! 4. Full LC Jac for BOA geometry      4 = 1 + 2
      REAL(fpk) :: COLUMNWF_SPHERCORR    ( MAX_ATMOSWFS, MAXSTOKES, MAXGEO ) ! 5. Full LC Jac with Sphericity corr  5 = 1 + 3

!  Surface Jacobian output with explanations

      REAL(fpk) :: FO_SURFACEWF_BOAGEOM   ( MAX_SURFACEWFS, MAXSTOKES, MAXGEO ) ! 1. FO (SS/DB) LC Jac for BOA geometry
      REAL(fpk) :: MS_SURFACEWF_BOAGEOM   ( MAX_SURFACEWFS, MAXSTOKES, MAXGEO ) ! 2. MS LC Jac BOA geometry
      REAL(fpk) :: MS_SURFACEWF_SPHERCORR ( MAX_SURFACEWFS, MAXSTOKES, MAXGEO ) ! 3. MS LC Jac with 2-point corr
      REAL(fpk) :: SURFACEWF_BOAGEOM      ( MAX_SURFACEWFS, MAXSTOKES, MAXGEO ) ! 4. Full LC Jac for BOA geometry      4 = 1 + 2
      REAL(fpk) :: SURFACEWF_SPHERCORR    ( MAX_SURFACEWFS, MAXSTOKES, MAXGEO ) ! 5. Full LC Jac with Sphericity corr  5 = 1 + 3

!  Twilight flag and exception handling
!    ==> Twilight flag active when TOA-geometry SZA > 90 degrees (VLIDORT can't handle this)
      
      LOGICAL       :: Twilight_Flag, FAIL
      CHARACTER*256 :: Local_Message, Local_Action

!  Other help variables

      INTEGER          :: I, IA, IS, IV, ICASE, IGEO, NGEO, CorrType, &
                          NSZA, NVZA, NRAA
      REAL(FPK)        :: SZA, SZA_SET(4), VZA, VZA_SET(12), RAA, RAA_SET(3), &
                          ALBEDO
      CHARACTER*1      :: MSType              
      CHARACTER*2      :: CALB, CFINE
      CHARACTER*6      :: Title
      CHARACTER*8      :: C8
      CHARACTER*100    :: Input_Filename, Output_FileName

!  Start of code
!  =============

!  Define some arrays
!  @@@ Rob Fix 6/13/13. CMASK is for the AEROSOLS, RMASK for RAYLEIGH.

      GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)
      RMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)    ! This Rayleigh
!      CMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)
      CMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)    ! OUTPUT OF RTS MIE
      SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /)

      SZA_SET = (/  0.0_fpk, 40.0_fpk, 70.0_fpk, 80.0_fpk /)
      VZA_SET = (/ 20.0_fpk, 30.0_fpk, 40.0_fpk, 50.0_fpk, 60.0_fpk, 65.0_fpk, &
                   70.0_fpk, 72.0_fpk, 74.0_fpk, 76.0_fpk, 78.0_fpk, 80.0_fpk /)
      RAA_SET = (/  0.0_fpk, 90.0_fpk, 180.0_fpk /)

!  Initialize error file output flag

      OPENFILEFLAG = .false.

!  Some driver control options

      Do_Aerosols  = .true.
      Do_Albedo    = .true.

      Albedo       = 0.2d0

!  Begin type of multiple-scatter correction loop

      DO CorrType = 1,3

        IF (CorrType .eq. 1) THEN
          MSType = '2'
        ELSEIF (CorrType .eq. 2) THEN
          MSType = '3'
        ELSEIF (CorrType .eq. 3) THEN
          MSType = 'M'
        ENDIF

!  Begin case loop [ Case 1 (Satellite) or Case 2 (Ground-based) ]

      DO ICASE = 1,2

      if ( ICASE .eq. 1 ) Input_Filename = 'vlidort_sc_test/V2p8p3_MSST_ReadInput_ToaUp.cfg'
      if ( ICASE .eq. 2 ) Input_Filename = 'vlidort_sc_test/V2p8p3_MSST_ReadInput_BoaDn.cfg'

!  VLIDORT control Read input, abort if failed

      CALL VLIDORT_L_INPUT_MASTER ( Trim(Input_Filename), & ! Input
          VLIDORT_FixIn,        & ! Outputs
          VLIDORT_ModIn,        & ! Outputs
          VLIDORT_LinFixIn,     & ! Outputs
          VLIDORT_LinModIn,     & ! Outputs
          VLIDORT_InputStatus )   ! Outputs

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( 'VLIDORT_LCS_SpherCorr_ReadInput.log', VLIDORT_InputStatus )

!  Initialize some input variables not handled by VLIDORT_L_INPUT_MASTER
!     These are the same for all Local tasks

      CALL VLIDORT_Sup_Init    ( VLIDORT_Sup )
      CALL VLIDORT_LinSup_Init ( VLIDORT_LinSup )

!  Set Proxies (saved values)

      nlayers = VLIDORT_FixIn%Cont%TS_nlayers
      nstokes = VLIDORT_FixIn%Cont%TS_nstokes
      nfine   = VLIDORT_FixIn%Cont%TS_nfinelayers

!  1/31/21. Version 2.8.3, New variables must be set by hand

      VLIDORT_FixIn%Cont%TS_ASYMTX_TOLERANCE = 1.0d-20

!  1/31/21. Version 2.8.3, Set the DO_MSSTS flag to generate output for the MS sphericity corrections

      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .true.

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

!  Identify extension for Rayleigh vs. Aerosols

      if ( Do_Aerosols ) then
        C8 = 'Aero6lyr'
      else
        C8 = 'RayOnly'
      endif

!  Identify NFINE file extension

      write(CFINE,'(I2.2)') NFINE

!  Identify Albedo file extension

      if ( Do_Albedo ) then
        write(CALB,'(I2.2)') INT(100.0d0*Albedo)
      else
        write(CALB,'(I2.2)') INT(100.0d0*VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO)
      endif

!  Write to screen

      if ( MSType.eq.'2' .or. MSType.eq.'3' ) then
        write(*,'(/a/)')'This is a ' // MSType // '-pt sphericity correction test for a ' // &
                        ' Scenario, Albedo '//CALB//'%'
      else
        write(*,'(/A/)')'This is a MULTI-pt sphericity correction test for the ' // &
                        ' Scenario, Albedo '//CALB//'%'
      endif

!  Open output files for normal type of results
!mick fix 1/5/2021 - added IF conditions for Q & U (as in "VLIDORT_SpherCorr_tester_V2.f90")

      if ( ICASE .eq. 1 ) then
        Output_FileName = 'vlidort_sc_test/results_LCS_SpherCorr_test_' // &
                           MSType // '-pt_' // trim(C8) // '_NF_' // CFINE // '_Alb_' // CALB
        if (nstokes.eq.1) then
          Open(701,file=Trim(Output_FileName)//'_I000.all',status='replace' )
        elseif (nstokes.gt.1) then
          Open(701,file=Trim(Output_FileName)//'_IQU0.all',status='replace' )
        else
          write(*,*) 'Error: test set up for nstokes 1 or 3'
          stop
        endif
      endif

!  Define files headers for the current output file
!    -- 1. Standard  VLIDORT: Single   scatter result with BOA geometry input
!    -- 2. Standard  VLIDORT: Multiple scatter result with BOA geometry input
!    -- 3. Standard  VLIDORT: SS + MS complete result with BOA geometry input (3 = 1 + 2)
!    -- 4. SpherCorr VLIDORT: Multiple scatter result
!    -- 5. SpherCorr VLIDORT: SS + MS complete result (5 = 1 + 4)
!    -- 6. relative difference (%) between SpherCorr and Standard VLIDORT results (MS FIELD ONLY)
!    -- 7. relative difference (%) between SpherCorr and Standard VLIDORT results (TOTAL FIELD)

      if ( ICASE .eq. 1 ) Title = 'ToaUp:'
      if ( ICASE .eq. 2 ) Title = 'BoaDn:'

      write(701,'(/a)') Title
      write(701,'(a,31x,7a)') '    ','                Std VLIDORT                     '      ,&
                              '         SC VLIDORT           ' ,'          Rel diff (%)'

      write(701,'(12a)')      ' Geo','   SZA  ','   VZA  ','   RAA  ',' STK ',&
                                     '       SS       ','       MS       ','     SS+MS      ',&
                                                        '       MS       ','     SS+MS      ',&
                                                        '       MS       ','     SS+MS'

!  Linearization proxies

      DO_COLUMN_LINEARIZATION  = .TRUE. ; layer_vary_number = 0
      DO_SURFACE_LINEARIZATION = .TRUE. ; layer_vary_flag   = .false.
      layer_vary_flag(1:nlayers) = .true.
      if ( do_aerosols ) then
         layer_vary_number(1:nlayers) = 2
         n_totalcolumn_wfs = 2
      else
         layer_vary_number(1:nlayers) = 1
         n_totalcolumn_wfs = 1
      endif
      n_surface_wfs     = 1

!  Set masking limits

      if ( nstokes .eq. 1 ) ngreekmat_entries = 1
      if ( nstokes .eq. 3 ) ngreekmat_entries = 5
      if ( nstokes .eq. 4 ) ngreekmat_entries = 8

!  Get the pre-prepared atmosphere
!  -------------------------------

      height_grid = 0.0d0
      open(45,file='vlidort_sc_test/input_atmos.dat',status='old' )
      read(45,'(i5,e18.9)') ldum, raymoms(0,1)
      read(45,'(i5,e18.9)') ldum, raymoms(1,1)
      read(45,'(i5,e18.9)') ldum, raymoms(2,1)
      do n = 1, nlayers
        read(45,'(i4,f12.5,5e16.7)') ndum,height_grid(n),molext(n),molomg(n),kd,kd,kd
        raymoms(0,n)=raymoms(0,1)
        raymoms(1,n)=raymoms(1,1)
        raymoms(2,n)=raymoms(2,1)
        if(molomg(n).gt.0.99999999d0) molomg(n)=0.99999999d0
      enddo
      height_grid(0) = 60.0d0
      close(45)

!  Aerosols
!  --------

!  Add Aerosols bottom 6 layers, spread evenly
!     Gamma-distributed, 0.782 Microns, n = 1.43, reff = 1.05, veff = 0.07

      if ( do_aerosols ) then
        open(1,file='vlidort_sc_test/ProblemIII.Moms',status = 'old')
        do L = 1, 5
          read(1,*)
        enddo
        LCHECK = 0 ; L = -1 ; cutoff = 1.0d-06
        DO while (L.lt.106 .and. LCHECK.lt.2)
          read(1,*)LDUM,amoms(1:6)
          if ( amoms(1).gt.cutoff ) then
            L = L + 1 ; aermoms(1:6,L) = amoms(1:6)
            if ( Lcheck.eq.1 ) Lcheck = Lcheck + 1
          else
            Lcheck = Lcheck + 1
            if ( Lcheck .eq. 1 ) then
              L = L + 1 ; aermoms(1:6,L) = amoms(1:6)
            endif
          endif
        enddo
        close(1)
        N_Aermoms = L
      endif

!  Rayleigh scattering
!  -------------------

!  Rayleigh scattering law

      depol   = ( 1.0d0 - 2.0d0*raymoms(2,NLAYERS) ) / &
                ( 1.0d0 + raymoms(2,NLAYERS) )
      beta2 = raymoms(2,NLAYERS)
      PROBLEM_RAY = 0.0d0
      PROBLEM_RAY(1,0) =  1.D0
      PROBLEM_RAY(4,1) =  3.D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
      PROBLEM_RAY(1,2) =  beta2
      PROBLEM_RAY(5,2) =  - SQRT(6.D0) * beta2
      PROBLEM_RAY(2,2) =  6.D0 * beta2

!  Initialise optical properties

      deltau_vert_input    = zero
      omega_total_input    = zero
      greekmat_total_input = zero

!  Initialize linearized inputs (local)

      L_deltau_vert_input    = zero
      L_omega_total_input    = zero
      L_greekmat_total_input = zero

!  Wavelength. Override the value in the configuration file (dummy, not needed)

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = 0.333d0

!  Set albedo

      if ( Do_Albedo ) VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = ALBEDO
 
!  Aerosols occupy the lowest 6 layers

      if ( do_aerosols ) then
        n6 = nlayers - 6
      else
        n6 = nlayers              ! RAYLEIGH ONLY
      endif

!  Fix optical properties for the non-aerosol layers

      do n = 1, n6
         deltau_vert_input(n) = molext(n)
         omega_total_input(n) = molomg(n)
         DO K = 1, ngreekmat_entries
           GK = gmask(k) ; rk = rmask(k) ! ; ck = cmask(k)
           DO L = 0, 2
             GREEKMAT_TOTAL_INPUT(L,n,gk) = PROBLEM_RAY(rk,L)
           ENDDO
           GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0
         ENDDO
         greekmat_total_input(0,n,1) = 1.0d0
         ratio1 = 1.0d0 - molomg(n)
         L_deltau_vert_input(1,n) =   ratio1
         L_omega_total_input(1,n) = - ratio1
      enddo      

!  Fix layers with aerosol

      if ( Do_Aerosols ) then
        waer = 0.99d0 ; taer = 0.5d0
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
            GK = gmask(k) ; rk = rmask(k) ; ck = cmask(k) ; sk = real(smask(k),fpk)
            DO L = 0, 2
              MOM = RAYWT*PROBLEM_RAY(RK,L) + AERWT * SK * AERMOMS(CK,L)
              greekmat_total_input(L,N,GK) = MOM
            ENDDO
            DO L = 3, N_AERMOMS
              MOM = AERWT * SK * AERMOMS(CK,L)
              greekmat_total_input(L,n,GK) = MOM
            ENDDO
          ENDDO
          greekmat_total_input(0,n,1) = 1.0d0
          DO K = 1, ngreekmat_entries
            GK = gmask(k) ; rk = rmask(k) ; ck = cmask(k) ; sk = real(smask(k),fpk)
            do L = 1, N_AERMOMS
              ratio2 = aersca / totsca
              MOM = greekmat_total_input(L,N,GK)
              IF ( MOM .ne.zero ) L_greekmat_total_input(2,L,n,GK) = ratio2 * ( (SK * aermoms(CK,L)/mom) - 1.0d0 )
            enddo
          enddo
        enddo
      Endif

!  1/31/20. Version 2.8.3.
!    -- set the Rayleigh-only, delta-m Scaling and Double-convergence flags

      if ( Do_Aerosols ) then
         ngreek_moments_input = N_Aermoms
         VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY    = .false.
         VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST  = .true.
         VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING   = .true.
      else
         ngreek_moments_input = 2
         VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY    = .true.
         VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST  = .false.
         VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING   = .false.
      endif

!  Copy to type-structure inputs

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = ngreek_moments_input
      VLIDORT_FixIn%Chapman%TS_height_grid        = height_grid

!  Copy to optical property type-structure inputs

      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input

      VLIDORT_FixIn%Optical%TS_FMATRIX_UP = zero ! zero here
      VLIDORT_FixIn%Optical%TS_FMATRIX_DN = zero ! zero here

!  Copy Linearization inputs

      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION  = .false.
      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION   = DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION  = DO_SURFACE_LINEARIZATION
      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS          = n_totalcolumn_wfs
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS              = n_surface_wfs
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG            = layer_vary_flag
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER          = layer_vary_number

      VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input     = L_deltau_vert_input
      VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input  = L_greekmat_total_input
      VLIDORT_LinFixIn%Optical%TS_L_omega_total_input     = L_omega_total_input

!  Need DO_FOCORR, DO_FOCORR_OUTGOING flags, always.
!    This is checked inside VLIDORT when the MSST flag is set.

      VLIDORT_ModIn%MBool%TS_DO_FOCORR          = .true.     
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR    = .false.
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING = .true.

!  debug input flag, for Dumping VLIDORT inputs

      do_debug_input = .false.

!  Start geometry loops
!  ====================

      igeo = 0

!  Trim local # of geometries here!

      NRAA = 3
      NSZA = 4
      NVZA = 12

!  Start RAA loop

      DO ia = 1,NRAA

!  Set RAA according to loop value

      raa = raa_set(ia)
      if( ia.eq.1 ) raa = 0.001d0 
      if( ia.eq.NRAA ) raa = 179.999d0
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,3) = RAA
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1)          = RAA

!  Start SZA loop

      DO is = 1,NSZA

!  Set SZA according to loop value

      sza = sza_set(is)
      if( is.eq.1 ) sza = 0.001d0
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1) = SZA
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(1)              = SZA

!  Start VZA loop

      DO iv= 1,NVZA

!  Set VZA according to loop value

      vza = vza_set(iv) !; if( iv.eq.1 ) vza = 0.001d0 
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,2) = VZA
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1)   = VZA

!  Progress

      igeo = igeo + 1
      if (iv .eq. 1) write(*,*)
      write(*,'(3(a,f6.2),a,i3)') ' Doing  SZA = ',sza,'  VZA = ',vza,'  RAA = ',raa,'  ==> IGEO = ',igeo

! Find critical VZA for a given SZA
!   sza = 90.0d0 - 85.0d0
!   write(*,*)'Critang = ',atan(sin(sza*deg_to_rad)/(cos(sza*deg_to_rad)-(6371.0d0/6471.0d0)))/deg_to_rad
!   stop

!  Call to dedicated sphericity correction routine (either 2-point, 3-point or multi-point)

      if ( MSType .eq. '2' ) then
        CALL MS2pt_LCS_SPHERCORR_V1 ( do_debug_input, &
            VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                        & ! VLIDORT Standard   Inputs
            VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup,               & ! VLIDORT Linearized Inputs
            FO_STOKES_BOAGEOM     (:,igeo),   MS_STOKES_BOAGEOM(:,igeo),      & ! Standard Component outputs (FO/MS)
            MS_STOKES_SPHERCORR   (:,igeo),                                   & ! Standard Component outputs (FO/MS)
            FO_COLUMNWF_BOAGEOM   (:,:,igeo), MS_COLUMNWF_BOAGEOM(:,:,igeo),  & ! LC Jacobian Component outputs (FO/MS)
            MS_COLUMNWF_SPHERCORR (:,:,igeo),                                 & ! LC Jacobian Component outputs (FO/MS)
            FO_SURFACEWF_BOAGEOM  (:,:,igeo), MS_SURFACEWF_BOAGEOM(:,:,igeo), & ! LS Jacobian Component outputs (FO/MS)
            MS_SURFACEWF_SPHERCORR(:,:,igeo),                                 & ! LS Jacobian Component outputs (FO/MS)
            VLIDORT_Out, VLIDORT_LinOut,                                      & ! VLIDORT and Standard SpherCorr output
            STOKES_BOAGEOM        (:,igeo),   STOKES_SPHERCORR    (:,igeo),   & ! VLIDORT and Standard SpherCorr output
            COLUMNWF_BOAGEOM      (:,:,igeo), COLUMNWF_SPHERCORR  (:,:,igeo), & ! Linearized SpherCorr output
            SURFACEWF_BOAGEOM     (:,:,igeo), SURFACEWF_SPHERCORR (:,:,igeo), & ! Linearized SpherCorr output
            Twilight_Flag, Fail, Local_Message, Local_Action )                  ! Exceptions and errors

      else if ( MSType .eq. '3' ) then
        CALL MS3pt_LCS_SPHERCORR_V1 ( do_debug_input, &
            VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                        & ! VLIDORT Standard   Inputs
            VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup,               & ! VLIDORT Linearized Inputs
            FO_STOKES_BOAGEOM     (:,igeo),   MS_STOKES_BOAGEOM(:,igeo),      & ! Standard Component outputs (FO/MS)
            MS_STOKES_SPHERCORR   (:,igeo),                                   & ! Standard Component outputs (FO/MS)
            FO_COLUMNWF_BOAGEOM   (:,:,igeo), MS_COLUMNWF_BOAGEOM(:,:,igeo),  & ! LC Jacobian Component outputs (FO/MS)
            MS_COLUMNWF_SPHERCORR (:,:,igeo),                                 & ! LC Jacobian Component outputs (FO/MS)
            FO_SURFACEWF_BOAGEOM  (:,:,igeo), MS_SURFACEWF_BOAGEOM(:,:,igeo), & ! LS Jacobian Component outputs (FO/MS)
            MS_SURFACEWF_SPHERCORR(:,:,igeo),                                 & ! LS Jacobian Component outputs (FO/MS)
            VLIDORT_Out, VLIDORT_LinOut,                                      & ! VLIDORT and Standard SpherCorr output
            STOKES_BOAGEOM        (:,igeo),   STOKES_SPHERCORR    (:,igeo),   & ! VLIDORT and Standard SpherCorr output
            COLUMNWF_BOAGEOM      (:,:,igeo), COLUMNWF_SPHERCORR  (:,:,igeo), & ! Linearized SpherCorr output
            SURFACEWF_BOAGEOM     (:,:,igeo), SURFACEWF_SPHERCORR (:,:,igeo), & ! Linearized SpherCorr output
            Twilight_Flag, Fail, Local_Message, Local_Action )                  ! Exceptions and errors

      else if ( MSType .eq. 'M' ) then
        CALL MSMpt_LCS_SPHERCORR_V1 ( do_debug_input, & 
            VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                        & ! VLIDORT Standard   Inputs
            VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup,               & ! VLIDORT Linearized Inputs
            FO_STOKES_BOAGEOM     (:,igeo),   MS_STOKES_BOAGEOM(:,igeo),      & ! Standard Component outputs (FO/MS)
            MS_STOKES_SPHERCORR   (:,igeo),                                   & ! Standard Component outputs (FO/MS)
            FO_COLUMNWF_BOAGEOM   (:,:,igeo), MS_COLUMNWF_BOAGEOM(:,:,igeo),  & ! LC Jacobian Component outputs (FO/MS)
            MS_COLUMNWF_SPHERCORR (:,:,igeo),                                 & ! LC Jacobian Component outputs (FO/MS)
            FO_SURFACEWF_BOAGEOM  (:,:,igeo), MS_SURFACEWF_BOAGEOM(:,:,igeo), & ! LS Jacobian Component outputs (FO/MS)
            MS_SURFACEWF_SPHERCORR(:,:,igeo),                                 & ! LS Jacobian Component outputs (FO/MS)
            VLIDORT_Out, VLIDORT_LinOut,                                      & ! VLIDORT and Standard SpherCorr output
            STOKES_BOAGEOM        (:,igeo),   STOKES_SPHERCORR    (:,igeo),   & ! VLIDORT and Standard SpherCorr output
            COLUMNWF_BOAGEOM      (:,:,igeo), COLUMNWF_SPHERCORR  (:,:,igeo), & ! Linearized SpherCorr output
            SURFACEWF_BOAGEOM     (:,:,igeo), SURFACEWF_SPHERCORR (:,:,igeo), & ! Linearized SpherCorr output
            Twilight_Flag, Fail, Local_Message, Local_Action )                  ! Exceptions and errors
      endif

!  Exception handling, write-up (optional)

      IF ( Fail ) THEN
        write(*,'(A)')'Failure from Call to MS'//MStype//'pt_LCS_SPHERCORR_V1, here are 2 local messages - '
        write(*,'(A)')'  --'//Trim(Local_Message)
        write(*,'(A)')'  --'//Trim(Local_Action)
        write(*,'(A)')'Now look at "VLIDORT_LCS_SpherCorr_tester_V1.log" file for error details'
        CALL VLIDORT_WRITE_STATUS ('VLIDORT_LCS_SpherCorr_tester_V1.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )
        STOP 'Program stopped'
      ENDIF

!  End geometry loops

      ENDDO !VZA
      ENDDO !SZA
      ENDDO !RAA

      ngeo = igeo
      write(*,'(/a,i3.3/)') 'ngeo = ',ngeo

!  Write results to file
!  ---------------------

!  One line for each SZA/VZA/RAA BOA-geometry calculation

!    -- 1. Standard  VLIDORT: Single   scatter result with BOA geometry input
!    -- 2. Standard  VLIDORT: Multiple scatter result with BOA geometry input
!    -- 3. Standard  VLIDORT: SS + MS complete result with BOA geometry input (3 = 1 + 2)
!    -- 4. SpherCorr VLIDORT: Multiple scatter result
!    -- 5. SpherCorr VLIDORT: SS + MS complete result (5 = 1 + 4)
!    -- 6. relative difference (%) between SpherCorr and Standard VLIDORT results (MS FIELD ONLY)
!    -- 7. relative difference (%) between SpherCorr and Standard VLIDORT results (TOTAL FIELD)
!    -- No result for Twilight cases (SZA > 90 at TOA or BOA destination). Set default value 99
!mick fix 1/5/2021 - added IF conditions for Q & U (as in "VLIDORT_SpherCorr_tester_V3.f90")
!                  - added IF condition for case with zero U-component

      if ( .not. Twilight_Flag ) then
        write(701,*)
        igeo = 0
        DO ia = 1,NRAA
        DO is = 1,NSZA
        DO iv = 1,NVZA
          if (iv .eq. 1) write(701,*)
          igeo = igeo + 1
          write(701,155) igeo,sza_set(is),vza_set(iv),raa_set(ia),'  I  ',&
              FO_COLUMNWF_BOAGEOM(1,1,igeo), MS_COLUMNWF_BOAGEOM(1,1,igeo), COLUMNWF_BOAGEOM(1,1,igeo), &
              MS_COLUMNWF_SPHERCORR(1,1,igeo), COLUMNWF_SPHERCORR(1,1,igeo),&
              100.0d0*(1.0-(MS_COLUMNWF_SPHERCORR(1,1,igeo)/MS_COLUMNWF_BOAGEOM(1,1,igeo))),&
              100.0d0*(1.0-(COLUMNWF_SPHERCORR(1,1,igeo)/COLUMNWF_BOAGEOM(1,1,igeo)))
        ENDDO ; ENDDO ; ENDDO

        if (nstokes.gt.1) then
          write(701,*)
          igeo = 0
          DO ia = 1,NRAA
          DO is = 1,NSZA
          DO iv = 1,NVZA
            if (iv .eq. 1) write(701,*)
            igeo = igeo + 1
            write(701,155) igeo,sza_set(is),vza_set(iv),raa_set(ia),'  Q  ',&
                FO_COLUMNWF_BOAGEOM(1,2,igeo), MS_COLUMNWF_BOAGEOM(1,2,igeo), COLUMNWF_BOAGEOM(1,2,igeo), &
                MS_COLUMNWF_SPHERCORR(1,2,igeo), COLUMNWF_SPHERCORR(1,2,igeo),&
                100.0d0*(1.0-(MS_COLUMNWF_SPHERCORR(1,2,igeo)/MS_COLUMNWF_BOAGEOM(1,2,igeo))),&
                100.0d0*(1.0-(COLUMNWF_SPHERCORR(1,2,igeo)/COLUMNWF_BOAGEOM(1,2,igeo)))
          ENDDO ; ENDDO ; ENDDO
        endif

        if (nstokes.eq.3) then
          write(701,*)
          igeo = 0
          DO ia = 1,NRAA
          DO is = 1,NSZA
          DO iv = 1,NVZA
            if (iv .eq. 1) write(701,*)
            igeo = igeo + 1
            if (abs(COLUMNWF_BOAGEOM(1,3,igeo)) .gt. ZERO) then
              write(701,155) igeo,sza_set(is),vza_set(iv),raa_set(ia),'  U  ',&
                  FO_COLUMNWF_BOAGEOM(1,3,igeo), MS_COLUMNWF_BOAGEOM(1,3,igeo), COLUMNWF_BOAGEOM(1,3,igeo), &
                  MS_COLUMNWF_SPHERCORR(1,3,igeo), COLUMNWF_SPHERCORR(1,3,igeo),&
                  100.0d0*(1.0-(MS_COLUMNWF_SPHERCORR(1,3,igeo)/MS_COLUMNWF_BOAGEOM(1,3,igeo))),&
                  100.0d0*(1.0-(COLUMNWF_SPHERCORR(1,3,igeo)/COLUMNWF_BOAGEOM(1,3,igeo)))
            else
              write(701,155) igeo,sza_set(is),vza_set(iv),raa_set(ia),'  U  ', &
                  FO_COLUMNWF_BOAGEOM(1,3,igeo), MS_COLUMNWF_BOAGEOM(1,3,igeo), COLUMNWF_BOAGEOM(1,3,igeo), &
                  MS_COLUMNWF_SPHERCORR(1,3,igeo), COLUMNWF_SPHERCORR(1,3,igeo),&
                  ZERO,ZERO
            endif
          ENDDO ; ENDDO ; ENDDO
        endif
        write(701,*)

      else
        write(*,'(/a/)') 'Warning: twilight flag set - geometry adjustment required!'
      endif
155   format(1x,i3,3f8.3,a5,7es16.8)

!  End case loop (TOA Up & BOA Dn)

      ENDDO

!  Close file

      close(701)

!  End type of MS correction loop (2-pt, 3-pt, & M-pt)

      ENDDO

!  Finish

      write(*,'(/a/)')'Main program VLIDORT_LCS_SPHERCORR_TESTER_V1 finished successfully'

      END PROGRAM VLIDORT_LCS_SPHERCORR_TESTER_V1
