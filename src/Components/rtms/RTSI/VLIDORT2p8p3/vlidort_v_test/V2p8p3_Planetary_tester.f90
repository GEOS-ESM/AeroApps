
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

      program Planetary_Tester

!  THIS IS THE VECTOR TEST, 6/18/19.
        
!  Introduced for Version 2.8.1. April 2019
!  ----------------------------------------
        
!  Module files for VLIDORT. Strict usage, Version 2.8

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_LIN_IO_DEFS_m

      USE VLIDORT_AUX_m,      Only : VLIDORT_READ_ERROR,     VLIDORT_WRITE_STATUS
      USE VLIDORT_INPUTS_m,   Only : VLIDORT_INPUT_MASTER,   VLIDORT_Sup_Init
      USE VLIDORT_L_INPUTS_m, Only : VLIDORT_L_INPUT_MASTER, VLIDORT_LinSup_Init

      USE VLIDORT_MASTERS_m
      USE VLIDORT_LCS_MASTERS_m
      USE VLIDORT_LPS_MASTERS_m
 
!  Implicit none

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

!  Control shorthand
!  =================

!  VLIDORT standard input preparation

      LOGICAL ::   DO_FOCORR, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING
      LOGICAL ::   DO_DELTAM_SCALING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING
      LOGICAL ::   DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY

!  Proxies

      INTEGER ::   NLAYERS, NSTOKES, NFINELAYERS, NGREEK_MOMENTS_INPUT
      INTEGER   :: N_GEOMETRIES, N_USER_LEVELS, NBEAMS
      DOUBLE PRECISION :: USER_LEVELS ( MAX_USER_LEVELS )

!  Control linearization

      LOGICAL         ::  DO_SIMULATION_ONLY
      LOGICAL         ::  DO_COLUMN_LINEARIZATION
      LOGICAL         ::  DO_PROFILE_LINEARIZATION
      LOGICAL         ::  DO_ATMOS_LINEARIZATION
      LOGICAL         ::  DO_SURFACE_LINEARIZATION
      LOGICAL         ::  LAYER_VARY_FLAG  (MAXLAYERS)
      INTEGER         ::  LAYER_VARY_NUMBER(MAXLAYERS)
      INTEGER         ::  N_SURFACE_WFS

!  Local Optical input Variables
!  =============================

!  Number of tasks

      INTEGER, PARAMETER :: MAXTASKS = 5

!  multilayer Height inputs

      DOUBLE PRECISION :: HEIGHT_GRID( 0:MAXLAYERS )

!  multilayer optical property (bulk) inputs

      DOUBLE PRECISION :: OMEGA_TOTAL_INPUT  ( MAXLAYERS, MAXTASKS )
      DOUBLE PRECISION :: DELTAU_VERT_INPUT  ( MAXLAYERS, MAXTASKS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations
!      Phasfunc proxies are new for Version 2.8. Not needed here

      DOUBLE PRECISION :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, 16, MAXTASKS )

!  Lambertian Surface control

      DOUBLE PRECISION :: LAMBERTIAN_ALBEDO (MAXTASKS)

!  Fmatrix stuff non needed here

      DOUBLE PRECISION :: FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 ) 
      DOUBLE PRECISION :: FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 ) 

!  Optical property linearizations

      DOUBLE PRECISION :: L_OMEGA_TOTAL_INPUT(MAX_ATMOSWFS,MAXLAYERS,MAXTASKS)
      DOUBLE PRECISION :: L_DELTAU_VERT_INPUT(MAX_ATMOSWFS,MAXLAYERS,MAXTASKS)
      DOUBLE PRECISION :: L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,MAXLAYERS,16,MAXTASKS)
 
!  Radiances

      DOUBLE PRECISION, dimension ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTASKS ) :: INTENSITY
      DOUBLE PRECISION, dimension ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTASKS ) :: PINTENSITY

!  Column weighting functions at all angles and optical depths

      DOUBLE PRECISION, dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES,  MAX_DIRECTIONS, MAXTASKS ) :: COLUMNWF

!  Profile weighting functions at all angles and optical depths

      DOUBLE PRECISION, dimension (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS,MAXTASKS) :: PROFILEWF

!  Other local variables
!  ---------------------

!  Task number and timing

      INTEGER          :: TASK, IRUN, NRUNS, JOB, NTASKS
      REAL             :: E1, E2, TIMING(MAXTASKS)
      
!  Flag for opening error output file

      LOGICAL ::          OPENFILEFLAG

!  Help variables

      integer                 :: n,n6,l,ndum,ldum, v, t, ndirs, nmi, ncwfs, npwfs, q, qfd, o1
      INTEGER                 :: K, CK, GK, RK, CMASK(8), GMASK(8), SMASK(8), RMASK(8)
      INTEGER                 :: NGREEKMAT_ENTRIES, LCHECK
      
      DOUBLE PRECISION        :: kd, gaer, waer, taer, parcel, raywt, aerwt
      double precision        :: aersca, aerext, molabs, molsca, totsca, totext, omega, delta, rayod
      double precision        :: molomg(maxlayers),molext(maxlayers)
      double precision        :: AMOMS(6), aermoms(6,0:maxmoments_input), lamb, eps, epsfac
      double precision        :: raymoms(0:2,maxlayers), ratio1, ratio2, beta2, depol

      
      DOUBLE PRECISION        :: PROBLEM_RAY(6,0:2), SK, CUTOFF, MOM
      
      DOUBLE PRECISION        :: R1, R2, R3, Z0, Z1, Z2, Z3, Z4, Z5, ZDIFF, W1, W2, Y0, Y1, Y2
      DOUBLE PRECISION        :: BZ0, BZ1, BZ2, BZDIFF, BZ3(max_geometries), BZ4(max_geometries), BZ5(max_geometries)
      DOUBLE PRECISION        :: PZ0, PZ1, PZ2, PZDIFF, PZ3, PZ4, PZ5

      DOUBLE PRECISION        :: TRANS, SPHER, SPHER1, TRANS1, PTRANS, PSPHER, PSPHER1, PTRANS1
      DOUBLE PRECISION        :: BTRANS(max_geometries), BSPHER(max_geometries), BSPHER1(max_geometries), BTRANS1(max_geometries)
      DOUBLE PRECISION        :: LC_TRANS (max_Atmoswfs,max_geometries), LC_SPHER (max_Atmoswfs,max_geometries)
      DOUBLE PRECISION        :: LC_SPHER1(max_Atmoswfs,max_geometries), LC_TRANS1(max_Atmoswfs,max_geometries)
      DOUBLE PRECISION        :: LC_Z3(max_Atmoswfs,max_geometries)
      DOUBLE PRECISION        :: LC_Z4(max_Atmoswfs,max_geometries), LC_Z5(max_Atmoswfs,max_geometries)
      DOUBLE PRECISION        :: LP_TRANS(maxlayers,max_Atmoswfs,max_geometries),  LP_SPHER (maxlayers,max_Atmoswfs,max_geometries)
      DOUBLE PRECISION        :: LP_SPHER1(maxlayers,max_Atmoswfs,max_geometries), LP_TRANS1(maxlayers,max_Atmoswfs,max_geometries)
      DOUBLE PRECISION        :: LP_Z3(maxlayers,max_Atmoswfs,max_geometries)
      DOUBLE PRECISION        :: LP_Z4(maxlayers,max_Atmoswfs,max_geometries), LP_Z5(maxlayers,max_Atmoswfs,max_geometries)

!  Only interested in Intensity

      O1 = 1 ; ntasks = maxtasks ; nruns = 1
      
!  Define some arrays
!  @@@ Rob Fix 6/13/13. CMASK is for the AEROSOLS, RMASK for RAYLEIGH.

      GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)
      RMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)    !  This   Rayleigh
!      CMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)
      CMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)    !  OUTPUT OF RTS MIE
      SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /)

!  ***********************************************
!  ***********************************************
!  PART 1. BASELINE CALCULATIONS WITH NO JACOBIANS
!  ***********************************************
!  ***********************************************
      
!  Initialize error file output flag

      OPENFILEFLAG = .false.

!  VLIDORT control Read input, abort if failed

      CALL VLIDORT_INPUT_MASTER ( &
        'vlidort_v_test/V2p8p3_VLIDORT_Planetary.cfg', & ! Input
        VLIDORT_FixIn,      & ! Outputs
        VLIDORT_ModIn,      & ! Outputs
        VLIDORT_InputStatus ) ! Outputs

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( 'V2p8p3_VLIDORT_Planetary.log', VLIDORT_InputStatus )

!  Initialize some input variables not handled by VLIDORT_INPUT_MASTER
!     These are the same for all Local tasks

      CALL VLIDORT_Sup_Init ( VLIDORT_Sup )

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
!   Version 2.8, now a configuration-file read. Override the value in there!!!

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = 0.333d0

!  Shorthand

      nstokes       = VLIDORT_FixIn%Cont%TS_nstokes
      nbeams        = VLIDORT_ModIn%MSunrays%TS_n_szangles
      n_user_levels = VLIDORT_FixIn%UserVal%TS_n_user_levels
      user_levels   = VLIDORT_ModIn%MUserVal%TS_user_levels

      if ( nstokes.eq.1) ngreekmat_entries = 1
      if ( nstokes.eq.3) ngreekmat_entries = 5

!  Get the pre-prepared atmosphere. NO AEROSOLS

      height_grid = zero
      nlayers = VLIDORT_FixIn%Cont%TS_nlayers
      open(45,file='vlidort_v_test/input_atmos.dat',status='old' )
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(1,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)') ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0) = 60.0d0
      do n = 1, nlayers
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n),molext(n),molomg(n),kd,kd,kd,kd
      enddo
      close(45)

!  Set the aerosol. if N6 = nlayers, it will be Rayleigh only.

      n6 = nlayers ; ngreek_moments_input = 2
!      n6 = nlayers - 6; ngreek_moments_input = 80
      
!  Add Aerosols bottom 6 layers, spread evenly
!     Gamma-distributed, 0.782 Microns, n = 1.43, reff = 1.05, veff = 0.07

      if ( n6 .lt. nlayers ) then
        open(1,file='vlidort_v_test/ProblemIII.Moms',status = 'unknown')
        read(1,*) ; read(1,*) ; read(1,*) ; read(1,*) ; read(1,*)
        LCHECK = 0 ; L = -1 ; cutoff = 1.0d-06
        DO while (L.lt.106 .and. LCHECK.lt.2)
          read(1,*)LDUM,(amoms(k),k=1,6)
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
        ngreek_moments_input = L
     endif
     
!  Rayleigh scattering law

      depol   = ( 1.0d0 - 2.0d0*raymoms(2,NLAYERS) ) / ( 1.0d0 + raymoms(2,NLAYERS) )
      beta2 = raymoms(2,NLAYERS)
      PROBLEM_RAY      = zero
      PROBLEM_RAY(1,0) = one

!  Old code, index was wrong
!      PROBLEM_RAY(3,1) =  3.D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
      PROBLEM_RAY(4,1) =  3.0D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
      PROBLEM_RAY(1,2) =  beta2
      PROBLEM_RAY(5,2) =  - sqrt(6.D0) * beta2
      PROBLEM_RAY(2,2) =  6.0D0 * beta2

!  Initialise optical proxies

      deltau_vert_input    = zero
      omega_total_input    = zero
      greekmat_total_input = zero
      Fmatrix_up           = zero ! 2p8
      Fmatrix_dn           = zero ! 2p8

!  Fill up optical properties, All tasks
!  -------------------------------------

!  Fix the non-aerosol layers

      do n = 1, n6
        deltau_vert_input(n,1:ntasks) = molext(n)
        omega_total_input(n,1:ntasks) = molomg(n)
        DO K = 1, ngreekmat_entries
          GK = gmask(k) ; rk = rmask(k)
          DO L = 0, 2
            GREEKMAT_TOTAL_INPUT(L,n,gk,1:ntasks) = PROBLEM_RAY(rk,L)
          ENDDO
          GREEKMAT_TOTAL_INPUT(0,n,1,1:ntasks) = 1.0d0
       ENDDO
      enddo

!  Fix layers with aerosol

      if ( n6 .lt. nlayers ) then
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
          deltau_vert_input(n,1:ntasks) = totext
          omega_total_input(n,1:ntasks) = totsca / totext
          DO K = 1, ngreekmat_entries
            GK = gmask(k); sk = dble(smask(k)); ck = cmask(k) ; rk = rmask(k)
            DO L = 0, 2
! Old         MOM = RAYWT*PROBLEM_RAY(ck,L) + AERWT*SK*AERMOMS(CK,L)
              MOM = RAYWT*PROBLEM_RAY(rk,L) + AERWT*SK*AERMOMS(CK,L)
              GREEKMAT_TOTAL_INPUT(L,N,GK,1:ntasks) = MOM
            ENDDO
            DO L = 3, NGREEK_MOMENTS_INPUT
              MOM = AERWT * SK * AERMOMS(CK,L)
              GREEKMAT_TOTAL_INPUT(L,n,gk,1:ntasks) = MOM
            ENDDO
          ENDDO
          GREEKMAT_TOTAL_INPUT(0,n,1,1:ntasks) = 1.0d0
        enddo
      endif
     
!  Surface

      lambertian_albedo(1) = 0.0d0
      lambertian_albedo(2) = 0.1d0
      lambertian_albedo(3) = 0.2d0
      lambertian_albedo(4) = 0.3d0
      lambertian_albedo(5) = 0.0d0

!  Copy to type-structure inputs, zero type structure opticals

      VLIDORT_FixIn%Cont%TS_nlayers                 = nlayers
      VLIDORT_ModIn%MCont%TS_ngreek_moments_input   = ngreek_moments_input
      VLIDORT_FixIn%Chapman%TS_height_grid          = height_grid
      VLIDORT_ModIn%MUserVal%TS_user_levels(1)      = zero
      VLIDORT_ModIn%MUserVal%TS_user_levels(2)      = dble(nlayers)

      VLIDORT_FixIn%Optical%TS_FMATRIX_UP           = zero ! zero here
      VLIDORT_FixIn%Optical%TS_FMATRIX_DN           = zero ! zero here
      VLIDORT_FixIn%Optical%TS_deltau_vert_input    = zero
      VLIDORT_ModIn%MOptical%TS_omega_total_input   = zero
      VLIDORT_FixIn%Optical%TS_greekmat_total_input = zero
      
!  Define additional proxies

      ndirs = max_directions
      nmi   = ngreek_moments_input
      
!  Start task loop
!  ===============

      do task = 1, ntasks
!      do task = 1, 1
!      do task = 5, 5

!  Copy to optical property type-structure inputs
!    This must now be done for each task

        VLIDORT_FixIn%Optical%TS_deltau_vert_input(1:nlayers)            = deltau_vert_input(1:nlayers,task)
        VLIDORT_ModIn%MOptical%TS_omega_total_input(1:nlayers)           = omega_total_input(1:nlayers,task)
        VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:nmi,1:nlayers,1:16) = greekmat_total_input(0:nmi,1:nlayers,1:16,task)
        VLIDORT_FixIn%Optical%TS_lambertian_albedo                       = lambertian_albedo(task)

!  Copy local variables to LIDORT (override config file)

        VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING  = .false.
        VLIDORT_ModIn%MBool%TS_DO_FOCORR          = .false.
        VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR    = .false.
        VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING = .false.

!  This is not set here, but may be invoked internally
        
        VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA = .false.

!  Task 1-4: 3 calls with different albedos, No FO correction, no delta-M scaling
!  Task 5:   1 call with no albedo but with Planetary problem flag set

        if ( task .lt. 5 ) then
           VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .false.
        else
           VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .true.
        endif

!  VLidort Call

        call CPU_TIME(e1)
        do_debug_input = .false.
        do irun = 1, nruns
           CALL VLIDORT_MASTER ( do_debug_input, &
             VLIDORT_FixIn, &
             VLIDORT_ModIn, &
             VLIDORT_Sup,   &
             VLIDORT_Out )
        enddo
        call CPU_TIME(e2)
        timing(task) = e2 - e1
        
!  Exception handling, write-up (optional)
!    Will generate file only if errors or warnings are encountered

        CALL VLIDORT_WRITE_STATUS ( &
          'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Write VLIDORT intensity to Local output (all tasks)

        n_geometries = VLIDORT_Out%Main%TS_n_geometries
        Intensity(1,1:n_geometries,1,task) = VLIDORT_Out%Main%TS_Stokes(1,1:n_geometries,1,1)
!        write(*,*)'Doing task # ', task,Intensity(1,1:n_geometries,1,task)

!  End task loop
        
      ENDDO

!  FIRST RESULTS Planetary experiment # 1 (NO JACOBIANS)
      
      open(1,file='vlidort_v_test/results_planetary_tester.all',status='unknown')

      write(1,55)' PART I : BASELINE RESULTS and timings for 10 calls to VLIDORT',&
                 ' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      write(1,56)'  Old Planetary Method    New Planetary Method     Direct Calculation                S and T values ',&
                 '  3 calls, 2 Equations     1 call, direct S/T      1 call with albedo',&
         'Timing(s), Scaling ',Sum(Timing(1:3)),Sum(Timing(1:3))/Timing(4),Timing(5), Timing(5)/Timing(4),Timing(4),1.0,&
                 ' # Geometry'

      do v = 1, n_geometries
        Z1 = Intensity(1,v,1,2) - Intensity(1,v,1,1)
        Z2 = Intensity(1,v,1,3) - Intensity(1,v,1,1)
        R1 = lambertian_albedo(2)
        R2 = lambertian_albedo(3)
        SPHER = ( (Z1/R1) - (Z2/R2) ) / ( Z1 - Z2 )
        TRANS = Z1 * ( one - R1 * SPHER ) / R1
        R3 = lambertian_albedo(4)
        Z3 = Intensity(1,v,1,1) + (R3 * TRANS / (ONe - R3 * SPHER))      
        TRANS1 = VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(o1,v)
        SPHER1 = VLIDORT_Out%Main%TS_PLANETARY_SBTERM
        Z4 = Intensity(1,v,1,1) + R3 * TRANS1 / (ONe - R3 * SPHER1)
        Z5 = Intensity(1,v,1,4)
        write(1,57)v, Z3, Z4, Z5, SPHER1, TRANS1
      enddo
55    format(/T30,A/T30,A)
56    format(/T20,A/T20,A//A,T20,3(3x,f8.5,2x,f8.5,3x)//A/)
57    format(2x,i3,3x,T20,3(4x,1pe16.7,4x),2(2x,1pe16.7,2x))

!  temporary stop/skip

!Go TO 6000

!      close(1)
!      STOP 'AFter first baseline test (No Jacobians)'
      
!  ***************************************
!  ***************************************
!  PART 2. CALCULATIONS WITH LCS JACOBIANS
!  ***************************************
!  ***************************************

!  Initialize

      OPENFILEFLAG = .false.

!  VLIDORT control Read input, abort if failed

      CALL VLIDORT_L_INPUT_MASTER ( &
        'vlidort_v_test/V2p8p3_VLIDORT_Planetary.cfg', & ! Input
        VLIDORT_FixIn,        & ! Outputs
        VLIDORT_ModIn,        & ! Outputs
        VLIDORT_LinFixIn,     & ! Outputs
        VLIDORT_LinModIn,     & ! Outputs
        VLIDORT_InputStatus )   ! Outputs

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( 'V2p8p3_VLIDORT_ReadInput.log', VLIDORT_InputStatus )

!  Initialize some input variables not handled by VLIDORT_L_INPUT_MASTER

      CALL VLIDORT_Sup_Init ( VLIDORT_Sup )
      CALL VLIDORT_LinSup_Init ( VLIDORT_LinSup )

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
!   Version 2.8, now a configuration-file read. Override the value in there!!!

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = 0.333d0

!  Shorthand

      nbeams        = VLIDORT_ModIn%MSunrays%TS_n_szangles
      n_user_levels = VLIDORT_FixIn%UserVal%TS_n_user_levels
      user_levels   = VLIDORT_ModIn%MUserVal%TS_user_levels

!  Set the aerosol. if N6 = nlayers, it will be Rayleigh only.

      n6 = nlayers ; ngreek_moments_input = 2
!      n6 = nlayers - 6; ngreek_moments_input = 80
      
!  Get the pre-prepared atmosphere

      nlayers = VLIDORT_FixIn%Cont%TS_nlayers
      open(45,file='vlidort_v_test/input_atmos.dat',status='old' )
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(1,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0) = 60.0d0
      do n = 1, nlayers
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n),molext(n),molomg(n),kd,kd,kd,kd
      enddo
      close(45)
      
!  Set the aerosol. if N6 = nlayers, it will be Rayleigh only.

      n6 = nlayers ; ngreek_moments_input = 2
!      n6 = nlayers - 6; ngreek_moments_input = 80
      
!  Add Aerosols bottom 6 layers, spread evenly
!     Gamma-distributed, 0.782 Microns, n = 1.43, reff = 1.05, veff = 0.07

      if ( n6 .lt. nlayers ) then
        open(1,file='vlidort_v_test/ProblemIII.Moms',status = 'unknown')
        read(1,*) ; read(1,*) ; read(1,*) ; read(1,*) ; read(1,*)
        LCHECK = 0 ; L = -1 ; cutoff = 1.0d-06
        DO while (L.lt.106 .and. LCHECK.lt.2)
          read(1,*)LDUM,(amoms(k),k=1,6)
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
        ngreek_moments_input = L
     endif

!  Define some additional proxies

      nmi   = ngreek_moments_input
      ncwfs = 1
      n_surface_wfs  = 0
      DO_SURFACE_LINEARIZATION = .FALSE.

!  Surface

      lambertian_albedo(1) = 0.0d0
      lambertian_albedo(2) = 0.1d0
      lambertian_albedo(3) = 0.2d0
      lambertian_albedo(4) = 0.3d0
      lambertian_albedo(5) = 0.0d0

!  Initialize linearized inputs (local)

      deltau_vert_input    = zero
      omega_total_input    = zero
      greekmat_total_input = zero
      L_deltau_vert_input    = zero
      L_omega_total_input    = zero
      L_greekmat_total_input = zero
      
!  Copy local control integers, height grid.
!      SAME for all tasks

      VLIDORT_FixIn%Cont%TS_nlayers                   = nlayers
      VLIDORT_ModIn%MCont%TS_ngreek_moments_input     = ngreek_moments_input
      VLIDORT_FixIn%Chapman%TS_height_grid(0:nlayers) = height_grid(0:nlayers)
      VLIDORT_ModIn%MUserVal%TS_user_levels(1) = zero
      VLIDORT_ModIn%MUserVal%TS_user_levels(2) = dble(nlayers)

!  JOB 1, baseline calculation
!  ===========================

      Job = 1 ; ntasks = 5
      write(*,*) 'Doing Job 1: ' // 'Baseline, Intensity + 2 column Jacobians'

!  Initialize linearized inputs

      DO_COLUMN_LINEARIZATION  = .TRUE.
      DO_PROFILE_LINEARIZATION = .FALSE.
      DO_ATMOS_LINEARIZATION   = .TRUE.
      DO_SIMULATION_ONLY       = .FALSE.
      do n = 1, nlayers
        layer_vary_number(n) = ncwfs
        layer_vary_flag(n)   = .true.
      enddo

!  Rayleigh layers

      do n = 1, n6
        deltau_vert_input(n,1:ntasks) = molext(n)
        omega_total_input(n,1:ntasks) = molomg(n)
        DO K = 1, ngreekmat_entries
          GK = gmask(k) ; rk = rmask(k) ! ; ck = cmask(k)
          DO L = 0, 2
            GREEKMAT_TOTAL_INPUT(L,n,gk,1:ntasks) = PROBLEM_RAY(rk,L)
          ENDDO
          GREEKMAT_TOTAL_INPUT(0,n,1,1:ntasks) = 1.0d0
        ENDDO
        ratio1 = 1.0d0 - molomg(n)
        L_deltau_vert_input(1,n,1:ntasks) =   ratio1
        L_omega_total_input(1,n,1:ntasks) = - ratio1
      enddo

!  Aerosol layers
      
      if ( n6 .lt. nlayers ) then
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
          deltau_vert_input(n,1:ntasks) = totext
          omega_total_input(n,1:ntasks) = omega
          ratio1 = molabs / totext
          L_deltau_vert_input(1,n,1:ntasks) =   ratio1
          L_omega_total_input(1,n,1:ntasks) = - ratio1
          ratio2 = aerext / totext 
          L_deltau_vert_input(2,n,1:ntasks) = ratio2
          L_omega_total_input(2,n,1:ntasks) = ratio2 * ((waer/omega) - 1.0d0)
          DO K = 1, ngreekmat_entries
            GK = gmask(k); sk = dble(smask(k)) ; ck = cmask(k) ; rk = rmask(k)
            DO L = 0, 2
              MOM = RAYWT*PROBLEM_RAY(rk,L) + AERWT*SK*AERMOMS(CK,L)
              GREEKMAT_TOTAL_INPUT(L,N,GK,1:ntasks) = MOM
            ENDDO
            DO L = 3, NGREEK_MOMENTS_INPUT
              MOM = AERWT * SK * AERMOMS(CK,L)
              GREEKMAT_TOTAL_INPUT(L,n,gk,1:ntasks) = MOM
            ENDDO
            do L = 1, ngreek_moments_input
              if ( greekmat_total_input(L,n,GK,1) .ne. 0.0d0 ) then
                ratio2 = aersca / totsca
                L_greekmat_total_input(2,L,n,GK,1:ntasks) = ratio2 * &
                  ( (SK*AERMOMS(CK,L)/greekmat_total_input(L,n,GK,1)) - 1.0d0 )
              endif
            enddo
          ENDDO
          GREEKMAT_TOTAL_INPUT(0,n,1,1:ntasks) = 1.0d0
        enddo
      endif
     
!  Copy local linearization control (override config file)

      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY          = .false.
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION    = .false.
      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION     = DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION      = DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION    = .false. 
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION            =  DO_COLUMN_LINEARIZATION
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)   = LAYER_VARY_FLAG(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS) = LAYER_VARY_NUMBER(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS                = 0
      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS            = ncwfs

!  Start task loop
!  ===============

      do task = 1, ntasks
!      do task = 2, 2
!      do task = 5, 5

!  Copy to optical property type-structure inputs
!    This must now be done for each task

        VLIDORT_FixIn%Optical%TS_deltau_vert_input(1:nlayers)          = deltau_vert_input(1:nlayers,task)
        VLIDORT_ModIn%MOptical%TS_omega_total_input(1:nlayers)         = omega_total_input(1:nlayers,task)
        VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:nmi,1:nlayers,1:16) = greekmat_total_input(0:nmi,1:nlayers,1:16,task)
        VLIDORT_FixIn%Optical%TS_lambertian_albedo                     = lambertian_albedo(task)
        VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input(1:ncwfs,1:nlayers) = &
                                   L_deltau_vert_input(1:ncwfs,1:nlayers,task)
        VLIDORT_LinFixIn%Optical%TS_L_omega_total_input(1:ncwfs,1:nlayers) = &               
                                   L_omega_total_input(1:ncwfs,1:nlayers,task)
        VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input(1:ncwfs,0:nmi,1:nlayers,1:16) = &
                                   L_greekmat_total_input(1:ncwfs,0:nmi,1:nlayers,1:16,task)

!  Copy local variables to VLIDORT (override config file)

        VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING  = .false.
        VLIDORT_ModIn%MBool%TS_DO_FOCORR          = .false.
        VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR    = .false.
        VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING = .false.

!  This is not set here, but may be invoked internally
        
        VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA = .false.

!  Task 1-4: 3 calls with different albedos, No FO correction, no delta-M scaling
!  Task 5:   1 call with no albedo but with Planetary problem flag set

        if ( task .lt. 5 ) then
           VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .false.
        else
           VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .true.
        endif

!  VLidort Call

        call CPU_TIME(e1)
        do_debug_input = .false.
        do irun = 1, nruns
          CALL VLIDORT_LCS_MASTER ( do_debug_input, &
            VLIDORT_FixIn,    & ! INPUTS
            VLIDORT_ModIn,    & ! INPUTS (possibly modified)
            VLIDORT_Sup,      & ! INPUTS/OUTPUTS
            VLIDORT_Out,      & ! OUTPUTS
            VLIDORT_LinFixIn, & ! INPUTS
            VLIDORT_LinModIn, & ! INPUTS (possibly modified)
            VLIDORT_LinSup,   & ! INPUTS/OUTPUTS
            VLIDORT_LinOut )    ! OUTPUTS
        enddo
        call CPU_TIME(e2)
        timing(task) = e2 - e1

!  Exception handling, write-up (optional)
!    Will generate file only if errors or warnings are encountered

        CALL VLIDORT_WRITE_STATUS ( &
          'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Copy to Local Baseline output
!  Write Baseline intensity to Local output (all tasks)

        n_geometries = VLIDORT_Out%Main%TS_n_geometries
        Intensity(1,1:n_geometries,1,task) = VLIDORT_Out%Main%TS_Stokes(1,1:n_geometries,1,1)
        Columnwf (1:ncwfs,1,1:n_geometries,1,task) =  &
             VLIDORT_LinOut%Col%TS_Columnwf(1:ncwfs,1,1:n_geometries,1,1)
        write(*,*)' ---- Done Baseline LCS task # ', task
        
!write(*,*)task,Columnwf (1,1,1,1,task),Intensity(1,1,1,task)

!  End task loop
        
      ENDDO
      
!  Geometry loop for LC linearized planetary experiment test
      
      do v = 1, n_geometries

!  Base intensity and S/T values

         BZ0 = Intensity(1,v,1,1)
         BZ1 = Intensity(1,v,1,2) - Intensity(1,v,1,1)
         BZ2 = Intensity(1,v,1,3) - Intensity(1,v,1,1) ; BZDIFF = BZ1 - BZ2
         R1 = lambertian_albedo(2) ; R2 = lambertian_albedo(3) ; R3 = lambertian_albedo(4)
         BSPHER(v) = ( (BZ1/R1) - (BZ2/R2) ) / BZDIFF
         BTRANS(v) = BZ1 * ( one - R1 * BSPHER(v) ) / R1
         BZ3(v) = Intensity(1,v,1,1) + R3 * BTRANS(v)  / (ONe - R3 * BSPHER(v))      
         BTRANS1(v) = VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(o1,v)
         BSPHER1(v) = VLIDORT_Out%Main%TS_PLANETARY_SBTERM
         BZ4(v) = Intensity(1,v,1,1) + R3 * BTRANS1(v) / (ONe - R3 * BSPHER1(v))
         BZ5(v) = Intensity(1,v,1,4)

!  Linearized Planetary experiment
!  Baseline LC Jacobian values of S/T and Intensity
 
         do q = 1, ncwfs
            Y0 =  Columnwf(q,1,v,1,1)
            Y1 =  Columnwf(q,1,v,1,2) - Y0
            Y2 =  Columnwf(q,1,v,1,3) - Y0
            W1 = BTRANS(v) * Y1 / Bz1 ; W2 = BTRANS(v) * Y2 / Bz2
            LC_SPHER(q,v) = ( W1 - W2 ) / BZDIFF 
            LC_TRANS(q,v) = W1 - LC_SPHER(q,v) * BZ1 
            LC_Z3(q,v) = Y0 + (BZ3(v)-BZ0) * ( LC_TRANS(q,v) + LC_SPHER(q,v) * (BZ3(v)-BZ0) ) / BTRANS(v)
            LC_TRANS1(q,v) = VLIDORT_LinOut%Col%TS_PLANETARY_TRANSTERM_COLWF(o1,v,q)
!write(*,*)v,LC_TRANS1(q,v)
            LC_SPHER1(q,v) = VLIDORT_LinOut%Col%TS_PLANETARY_SBTERM_COLWF(q)
            LC_Z4(q,v) = Y0 + (BZ4(v)-BZ0) * ( LC_TRANS1(q,v) + LC_SPHER1(q,v) * (BZ4(v)-BZ0) ) / BTRANS1(v)
            LC_Z5(q,v) = Columnwf(q,1,v,1,4)
         enddo

      enddo      

!STOP 'BASELINE'

!  FD CALCULATION, 1 COLUMN JACOBIAN ONLY
!  --------------------------------------

!  Write header for this part
      
      write(1,55)' PART II : LCS JACOBIAN RESULTS: Analytic WFS vs. FD JACOBIANS',&
                 ' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      write(1,66)'    Old Planetary Method          New Planetary Method           Direct Calculation',&
                 '                            S and T values ',&
                 '    3 calls, 2 Equations           1 call, direct S/T            1 call with albedo',&
                 '   AnalyticWF   FDJacobian       AnalyticWF   FDJacobian       AnalyticWF   FDJacobian', &
                 '         AnalyticWF  FDJacobian      AnalyticWF  FDJacobian',&
                 ' # Geometry'
      
!  Copy local linearization control (override config file)

      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY          = .true.
      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION     = .false.
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION      = .false.
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION            = .false.
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)   = .false.
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS) = 0
      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS            = 0

!  Initialize linearized inputs (local)

      L_deltau_vert_input    = zero
      L_omega_total_input    = zero
      L_greekmat_total_input = zero

      eps = 0.001_fpk ; epsfac = one + eps

      do n = 1, nlayers
         rayod = molext(n) * molomg(n)
         delta = molext(n) * epsfac - rayod * eps
         omega= rayod / delta
         deltau_vert_input(n,1:ntasks) = delta
         omega_total_input(n,1:ntasks) = omega
      enddo

!  Start task loop
!  ---------------

      do task = 1, 5

!  Copy to optical property type-structure inputs
!    This must now be done for each task

        VLIDORT_FixIn%Optical%TS_deltau_vert_input(1:nlayers)            = deltau_vert_input(1:nlayers,task)
        VLIDORT_ModIn%MOptical%TS_omega_total_input(1:nlayers)           = omega_total_input(1:nlayers,task)
        VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:nmi,1:nlayers,1:16) = greekmat_total_input(0:nmi,1:nlayers,1:16,task)
        VLIDORT_FixIn%Optical%TS_lambertian_albedo                       = lambertian_albedo(task)
      
!  Task 1-4: 3 calls with different albedos, No FO correction, no delta-M scaling
!  Task 5:   1 call with no albedo but with Planetary problem flag set

        VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA = .false.
        if ( task .lt. 5 ) then
           VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .false.
        else
           VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .true.
        endif

!  VLidort Call

        do_debug_input = .false.
        CALL VLIDORT_LCS_MASTER ( do_debug_input, &
            VLIDORT_FixIn,    & ! INPUTS
            VLIDORT_ModIn,    & ! INPUTS (possibly modified)
            VLIDORT_Sup,      & ! INPUTS/OUTPUTS
            VLIDORT_Out,      & ! OUTPUTS
            VLIDORT_LinFixIn, & ! INPUTS
            VLIDORT_LinModIn, & ! INPUTS (possibly modified)
            VLIDORT_LinSup,   & ! INPUTS/OUTPUTS
            VLIDORT_LinOut )    ! OUTPUTS

!  Exception handling, write-up (optional)
!    Will generate file only if errors or warnings are encountered

        CALL VLIDORT_WRITE_STATUS ( &
          'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Copy to Local Baseline output
!  Write Baseline intensity to Local output (all tasks)

        n_geometries = VLIDORT_Out%Main%TS_n_geometries
        PIntensity(1,1:n_geometries,1,task) = VLIDORT_Out%Main%TS_Stokes(1,1:n_geometries,1,1)
        
!  End task loop for FD perturbations
        
!write(*,*)task,PIntensity(1,1,1,task)

      ENDDO

!  Compute FD Jacobians and write-up results
      
      q = 1
      do v = 1, n_geometries
         PZ0 = PIntensity(1,v,1,1)
         PZ1 = PIntensity(1,v,1,2) - PIntensity(1,v,1,1)
         PZ2 = PIntensity(1,v,1,3) - PIntensity(1,v,1,1) ; PZDIFF = PZ1 - PZ2
         R1 = lambertian_albedo(2) ; R2 = lambertian_albedo(3) ; R3 = lambertian_albedo(4)
         PSPHER = ( (PZ1/R1) - (PZ2/R2) ) / PZDIFF
         PTRANS = PZ1 * ( one - R1 * PSPHER ) / R1
         PZ3 = PIntensity(1,v,1,1) + R3 * PTRANS  / (ONe - R3 * PSPHER)      
         PTRANS1 = VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(o1,v)
         PSPHER1 = VLIDORT_Out%Main%TS_PLANETARY_SBTERM
         PZ4 = PIntensity(1,v,1,1) + R3 * PTRANS1 / (ONe - R3 * PSPHER1)
         PZ5 = PIntensity(1,v,1,4)
         write(1,67)v,  LC_Z3(q,v), (PZ3-BZ3(v))/eps,  LC_Z4(q,v), (PZ4-BZ4(v))/eps, LC_Z5(q,v), (PZ5-BZ5(v))/eps,&
              LC_SPHER1(q,v), (PSPHER1-BSPHER1(v))/eps,  LC_TRANS1(q,v), (PTRANS1-BTRANS1(v))/eps
      enddo
66    format(/T20,A,A/T20,A/T20,A,A//A/)
67    format(2x,i3,3x,T19,3(1x,1p2e14.6,1x),2x,2(1x,1p2e13.5,1x))

!
!     close(1)
!     Stop 'Tempo LCS----'

!  ***************************************
!  ***************************************
!  PART 3. CALCULATIONS WITH LPS JACOBIANS
!  ***************************************
!  ***************************************

!  Skip point
!6000 continue

!  Initialize

      OPENFILEFLAG = .false.

!  VLIDORT control Read input, abort if failed

      CALL VLIDORT_L_INPUT_MASTER ( &
        'vlidort_v_test/V2p8p3_VLIDORT_Planetary.cfg', & ! Input
        VLIDORT_FixIn,        & ! Outputs
        VLIDORT_ModIn,        & ! Outputs
        VLIDORT_LinFixIn,     & ! Outputs
        VLIDORT_LinModIn,     & ! Outputs
        VLIDORT_InputStatus )   ! Outputs

      IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) &
        CALL VLIDORT_READ_ERROR ( 'V2p8p3_VLIDORT_ReadInput.log', VLIDORT_InputStatus )

!  Initialize some input variables not handled by VLIDORT_L_INPUT_MASTER

      CALL VLIDORT_Sup_Init ( VLIDORT_Sup )
      CALL VLIDORT_LinSup_Init ( VLIDORT_LinSup )

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
!   Version 2.8, now a configuration-file read. Override the value in there!!!

      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = 0.333d0

!  Shorthand

      nbeams        = VLIDORT_ModIn%MSunrays%TS_n_szangles
      n_user_levels = VLIDORT_FixIn%UserVal%TS_n_user_levels
      user_levels   = VLIDORT_ModIn%MUserVal%TS_user_levels

!  Get the pre-prepared atmosphere

      nlayers = VLIDORT_FixIn%Cont%TS_nlayers
      open(45,file='vlidort_v_test/input_atmos.dat',status='old' )
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(1,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0) = 60.0d0
      do n = 1, nlayers
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n),molext(n),molomg(n),kd,kd,kd,kd
      enddo
      close(45)
      
!  Set the aerosol. if N6 = nlayers, it will be Rayleigh only.

      n6 = nlayers ; ngreek_moments_input = 2
!      n6 = nlayers - 6; ngreek_moments_input = 80
      
!  Add Aerosols bottom 6 layers, spread evenly
!     Gamma-distributed, 0.782 Microns, n = 1.43, reff = 1.05, veff = 0.07

      if ( n6 .lt. nlayers ) then
        open(1,file='vlidort_v_test/ProblemIII.Moms',status = 'unknown')
        read(1,*) ; read(1,*) ; read(1,*) ; read(1,*) ; read(1,*)
        LCHECK = 0 ; L = -1 ; cutoff = 1.0d-06
        DO while (L.lt.106 .and. LCHECK.lt.2)
          read(1,*)LDUM,(amoms(k),k=1,6)
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
        ngreek_moments_input = L
     endif
     
!  Rayleigh scattering law

      depol   = ( 1.0d0 - 2.0d0*raymoms(2,NLAYERS) ) / ( 1.0d0 + raymoms(2,NLAYERS) )
      beta2 = raymoms(2,NLAYERS)
      PROBLEM_RAY      = zero
      PROBLEM_RAY(1,0) = one

!  Old code, index was wrong
!      PROBLEM_RAY(3,1) =  3.D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
      PROBLEM_RAY(4,1) =  3.0D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
      PROBLEM_RAY(1,2) =  beta2
      PROBLEM_RAY(5,2) =  - sqrt(6.D0) * beta2
      PROBLEM_RAY(2,2) =  6.0D0 * beta2

!  Define some additional proxies

      nmi   = ngreek_moments_input
      npwfs = 1
      n_surface_wfs  = 0
      DO_SURFACE_LINEARIZATION = .FALSE.
      
!  Surface

      lambertian_albedo(1) = 0.0d0
      lambertian_albedo(2) = 0.1d0
      lambertian_albedo(3) = 0.2d0
      lambertian_albedo(4) = 0.3d0
      lambertian_albedo(5) = 0.0d0

!  Initialize linearized inputs (local)

      L_deltau_vert_input    = zero
      L_omega_total_input    = zero
      L_greekmat_total_input = zero
      
!  Copy local control integers, height grid.
!      SAME for all tasks

      VLIDORT_FixIn%Cont%TS_nlayers                   = nlayers
      VLIDORT_ModIn%MCont%TS_ngreek_moments_input     = ngreek_moments_input
      VLIDORT_FixIn%Chapman%TS_height_grid(0:nlayers) = height_grid(0:nlayers)
      VLIDORT_ModIn%MUserVal%TS_user_levels(1) = zero
      VLIDORT_ModIn%MUserVal%TS_user_levels(2) = dble(nlayers)

!  JOB 1, baseline calculation
!  ===========================

      Job = 1 ; ntasks = 5
      write(*,*) 'Doing Job 1: ' // 'Baseline, Intensity + 1 profile Jacobians'

!  Initialize linearized inputs

      DO_PROFILE_LINEARIZATION = .TRUE.
      do n = 1, nlayers
        layer_vary_number(n) = npwfs
        layer_vary_flag(n)   = .true.
      enddo

!  Rayleigh layers

      do n = 1, n6
        deltau_vert_input(n,1:ntasks) = molext(n)
        omega_total_input(n,1:ntasks) = molomg(n)
        DO K = 1, ngreekmat_entries
          GK = gmask(k) ; rk = rmask(k) ! ; ck = cmask(k)
          DO L = 0, 2
            GREEKMAT_TOTAL_INPUT(L,n,gk,1:ntasks) = PROBLEM_RAY(rk,L)
          ENDDO
          GREEKMAT_TOTAL_INPUT(0,n,1,1:ntasks) = 1.0d0
        ENDDO
        ratio1 = 1.0d0 - molomg(n)
        L_deltau_vert_input(1,n,1:ntasks) =   ratio1
        L_omega_total_input(1,n,1:ntasks) = - ratio1
      enddo

!  Aerosol layers
      
      if ( n6 .lt. nlayers ) then
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
          deltau_vert_input(n,1:ntasks) = totext
          omega_total_input(n,1:ntasks) = omega
          ratio1 = molabs / totext
          L_deltau_vert_input(1,n,1:ntasks) =   ratio1
          L_omega_total_input(1,n,1:ntasks) = - ratio1
          ratio2 = aerext / totext 
          L_deltau_vert_input(2,n,1:ntasks) = ratio2
          L_omega_total_input(2,n,1:ntasks) = ratio2 * ((waer/omega) - 1.0d0)
          DO K = 1, ngreekmat_entries
            GK = gmask(k); sk = dble(smask(k)) ; ck = cmask(k) ; rk = rmask(k)
            DO L = 0, 2
              MOM = RAYWT*PROBLEM_RAY(rk,L) + AERWT*SK*AERMOMS(CK,L)
              GREEKMAT_TOTAL_INPUT(L,N,GK,1:ntasks) = MOM
            ENDDO
            DO L = 3, NGREEK_MOMENTS_INPUT
              MOM = AERWT * SK * AERMOMS(CK,L)
              GREEKMAT_TOTAL_INPUT(L,n,gk,1:ntasks) = MOM
            ENDDO
            do L = 1, ngreek_moments_input
              if ( greekmat_total_input(L,n,GK,1) .ne. 0.0d0 ) then
                ratio2 = aersca / totsca
                L_greekmat_total_input(2,L,n,GK,1:ntasks) = ratio2 * &
                  ( (SK*AERMOMS(CK,L)/greekmat_total_input(L,n,GK,1)) - 1.0d0 )
              endif
            enddo
          ENDDO
          GREEKMAT_TOTAL_INPUT(0,n,1,1:ntasks) = 1.0d0
        enddo
      endif
     
!  Copy local linearization control (override config file)

      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY          = .false.
      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION     = .false.
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION    = DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION      = DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION    = .false. 
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION            = DO_PROFILE_LINEARIZATION
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)   = LAYER_VARY_FLAG(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS) = LAYER_VARY_NUMBER(1:NLAYERS)
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS                = 0
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS           = npwfs
      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS            = 0

!  Start task loop
!  ===============

      do task = 1, ntasks
!      do task = 2, 2
!      do task = 5, 5

!  Copy to optical property type-structure inputs
!    This must now be done for each task

        VLIDORT_FixIn%Optical%TS_deltau_vert_input(1:nlayers)          = deltau_vert_input(1:nlayers,task)
        VLIDORT_ModIn%MOptical%TS_omega_total_input(1:nlayers)         = omega_total_input(1:nlayers,task)
        VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:nmi,1:nlayers,1:16) = greekmat_total_input(0:nmi,1:nlayers,1:16,task)
        VLIDORT_FixIn%Optical%TS_lambertian_albedo                     = lambertian_albedo(task)
        VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input(1:npwfs,1:nlayers) = &
                                   L_deltau_vert_input(1:npwfs,1:nlayers,task)
        VLIDORT_LinFixIn%Optical%TS_L_omega_total_input(1:npwfs,1:nlayers) = &               
                                   L_omega_total_input(1:npwfs,1:nlayers,task)
        VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input(1:npwfs,0:nmi,1:nlayers,1:16) = &
                                   L_greekmat_total_input(1:npwfs,0:nmi,1:nlayers,1:16,task)

!  Copy local variables to VLIDORT (override config file)

        VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING  = .false.
        VLIDORT_ModIn%MBool%TS_DO_FOCORR          = .false.
        VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR    = .false.
        VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING = .false.

!  This is not set here, but may be invoked internally
        
        VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA = .false.

!  Task 1-4: 3 calls with different albedos, No FO correction, no delta-M scaling
!  Task 5:   1 call with no albedo but with Planetary problem flag set

        if ( task .lt. 5 ) then
           VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .false.
        else
           VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .true.
        endif

!  VLidort Call

        call CPU_TIME(e1)
        do_debug_input = .false.
        do irun = 1, nruns
          CALL VLIDORT_LPS_MASTER ( do_debug_input, &
            VLIDORT_FixIn,    & ! INPUTS
            VLIDORT_ModIn,    & ! INPUTS (possibly modified)
            VLIDORT_Sup,      & ! INPUTS/OUTPUTS
            VLIDORT_Out,      & ! OUTPUTS
            VLIDORT_LinFixIn, & ! INPUTS
            VLIDORT_LinModIn, & ! INPUTS (possibly modified)
            VLIDORT_LinSup,   & ! INPUTS/OUTPUTS
            VLIDORT_LinOut )    ! OUTPUTS
        enddo
        call CPU_TIME(e2)
        timing(task) = e2 - e1

!  Exception handling, write-up (optional)
!    Will generate file only if errors or warnings are encountered

        CALL VLIDORT_WRITE_STATUS ( &
          '2p8p1_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Copy to Local Baseline output
!  Write Baseline intensity to Local output (all tasks)

        n_geometries = VLIDORT_Out%Main%TS_n_geometries
        Intensity(1,1:n_geometries,1,task) = VLIDORT_Out%Main%TS_Stokes(1,1:n_geometries,1,1)
        Profilewf (1:npwfs,1:nlayers,1,1:n_geometries,1,task) =  &
             VLIDORT_LinOut%Prof%TS_Profilewf(1:npwfs,1:nlayers,1,1:n_geometries,1,1)
!        write(*,*)'Doing task # ', task,Intensity(1,1,1,task),Profilewf (1,1:nlayers,1,1:n_geometries,1,task)
        write(*,*)' ---- Done Baseline LPS task # ', task
        
!  End task loop
        
      ENDDO
      
!  Geometry loop for LP linearized planetary experiment test
      
      do v = 1, n_geometries

!  Base intensity and S/T values

         BZ0 = Intensity(1,v,1,1)
         BZ1 = Intensity(1,v,1,2) - Intensity(1,v,1,1)
         BZ2 = Intensity(1,v,1,3) - Intensity(1,v,1,1) ; BZDIFF = BZ1 - BZ2
         R1 = lambertian_albedo(2) ; R2 = lambertian_albedo(3) ; R3 = lambertian_albedo(4)
         BSPHER(v) = ( (BZ1/R1) - (BZ2/R2) ) / BZDIFF
         BTRANS(v) = BZ1 * ( one - R1 * BSPHER(v) ) / R1
         BZ3(v) = Intensity(1,v,1,1) + R3 * BTRANS(v)  / (ONe - R3 * BSPHER(v))      
         BTRANS1(v) = VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(o1,v)
         BSPHER1(v) = VLIDORT_Out%Main%TS_PLANETARY_SBTERM
         BZ4(v) = Intensity(1,v,1,1) + R3 * BTRANS1(v) / (ONe - R3 * BSPHER1(v))
         BZ5(v) = Intensity(1,v,1,4)
         
!  Baseline LP Jacobian values of S/T and Intensity
               
         do n = 1, nlayers
           do q = 1, layer_vary_number(n)
             Y0 =  Profilewf(q,n,1,v,1,1)
             Y1 =  Profilewf(q,n,1,v,1,2) - Y0
             Y2 =  Profilewf(q,n,1,v,1,3) - Y0
             W1 = BTRANS(v) * Y1 / Bz1 ; W2 = BTRANS(v) * Y2 / Bz2
             LP_SPHER(n,q,v) = ( W1 - W2 ) / BZDIFF 
             LP_TRANS(n,q,v) = W1 - LP_SPHER(n,q,v) * BZ1 
             LP_Z3(n,q,v) = Y0 + (BZ3(v)-BZ0) * ( LP_TRANS(n,q,v) + LP_SPHER(n,q,v) * (BZ3(v)-BZ0) ) / BTRANS (v) 
             LP_TRANS1(n,q,v) = VLIDORT_LinOut%Prof%TS_PLANETARY_TRANSTERM_PROFWF(o1,v,n,q)
             LP_SPHER1(n,q,v) = VLIDORT_LinOut%Prof%TS_PLANETARY_SBTERM_PROFWF(n,q)
             LP_Z4(n,q,v) = Y0 + (BZ4(v)-BZ0) * ( LP_TRANS1(n,q,v) + LP_SPHER1(n,q,v) * (BZ4(v)-BZ0) ) / BTRANS1(v)
             LP_Z5(n,q,v) = Profilewf(q,n,1,v,1,4)
           enddo
         enddo

!  end geometry
         
      enddo

      !STOP 'BASELINE'
      
!  FD CALCULATION, ALL LAYERS
!  --------------------------
      
!  Copy local linearization control (override config file)

      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY          = .true.
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION    = .false.
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION      = .false.
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION            = .false.
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)   = .false.
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS) = 0
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS           = 0

!  Write header for this part
      
      write(1,55)' PART III : LPS JACOBIAN RESULTS: Analytic WFS vs. FD JACOBIANS',&
                 ' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      write(1,66)'    Old Planetary Method          New Planetary Method           Direct Calculation',&
                 '                            S and T values ',&
                 '    3 calls, 2 Equations           1 call, direct S/T            1 call with albedo',&
                 '   AnalyticWF   FDJacobian       AnalyticWF   FDJacobian       AnalyticWF   FDJacobian', &
                 '         AnalyticWF  FDJacobian      AnalyticWF  FDJacobian',&
                 ' # Geometry/Layer '

!  Initialize linearized inputs (local)

      L_deltau_vert_input    = zero
      L_omega_total_input    = zero
      L_greekmat_total_input = zero

      eps = 0.001_fpk ; epsfac = one + eps
      q = 1

!  Start profile Jacobian loop

      do qfd = 1, nlayers

        do n = 1, nlayers
          if ( qfd.eq.n ) then
            rayod = molext(n) * molomg(n)
            delta = molext(n) * epsfac - rayod * eps
            omega= rayod / delta
          else
            delta = molext(n)
            omega = molomg(n)
          endif
          deltau_vert_input(n,1:ntasks) = delta
          omega_total_input(n,1:ntasks) = omega
        enddo

!  Start task loop
!  ===============

        do task = 1, 5

!  Copy to optical property type-structure inputs
!    This must now be done for each task

          VLIDORT_FixIn%Optical%TS_deltau_vert_input(1:nlayers)            = deltau_vert_input(1:nlayers,task)
          VLIDORT_ModIn%MOptical%TS_omega_total_input(1:nlayers)           = omega_total_input(1:nlayers,task)
          VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:nmi,1:nlayers,1:16) = greekmat_total_input(0:nmi,1:nlayers,1:16,task)
          VLIDORT_FixIn%Optical%TS_lambertian_albedo                       = lambertian_albedo(task)
      
!  Task 1-4: 3 calls with different albedos, No FO correction, no delta-M scaling
!  Task 5:   1 call with no albedo but with Planetary problem flag set        

          VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA = .false.
          if ( task .lt. 5 ) then
            VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .false.
          else
            VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .true.
          endif

!  VLidort Call

          do_debug_input = .false.
          CALL VLIDORT_LPS_MASTER ( do_debug_input, &
            VLIDORT_FixIn,    & ! INPUTS
            VLIDORT_ModIn,    & ! INPUTS (possibly modified)
            VLIDORT_Sup,      & ! INPUTS/OUTPUTS
            VLIDORT_Out,      & ! OUTPUTS
            VLIDORT_LinFixIn, & ! INPUTS
            VLIDORT_LinModIn, & ! INPUTS (possibly modified)
            VLIDORT_LinSup,   & ! INPUTS/OUTPUTS
            VLIDORT_LinOut )    ! OUTPUTS

!  Exception handling, write-up (optional)
!    Will generate file only if errors or warnings are encountered

          CALL VLIDORT_WRITE_STATUS ( &
            'V2p8p3_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )

!  Copy to Local Baseline output
!  Write Baseline intensity to Local output (all tasks)

          n_geometries = VLIDORT_Out%Main%TS_n_geometries
          PIntensity(1,1:n_geometries,1,task) = VLIDORT_Out%Main%TS_Stokes(1,1:n_geometries,1,1)
!          write(*,*)'Doing WFTest task # ', task,qfd
                
!  End task loop
        
        ENDDO

!  Compute FD Jacobians and write-up results, layer by layer
        
        q = 1      
        do v = 1, n_geometries
          PZ0 = PIntensity(1,v,1,1)
          PZ1 = PIntensity(1,v,1,2) - PIntensity(1,v,1,1)
          PZ2 = PIntensity(1,v,1,3) - PIntensity(1,v,1,1) ; PZDIFF = PZ1 - PZ2
          R1 = lambertian_albedo(2) ; R2 = lambertian_albedo(3) ; R3 = lambertian_albedo(4)
          PSPHER = ( (PZ1/R1) - (PZ2/R2) ) / PZDIFF
          PTRANS = PZ1 * ( one - R1 * PSPHER ) / R1
          PZ3 = PIntensity(1,v,1,1) + R3 * PTRANS  / (ONe - R3 * PSPHER)      
          PTRANS1 = VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(o1,v)
          PSPHER1 = VLIDORT_Out%Main%TS_PLANETARY_SBTERM
          PZ4 = PIntensity(1,v,1,1) + R3 * PTRANS1 / (ONe - R3 * PSPHER1)
          PZ5 = PIntensity(1,v,1,4)
          write(1,68)v,qfd, LP_Z3(qfd,q,v),(PZ3-BZ3(v))/eps,LP_Z4(qfd,q,v),(PZ4-BZ4(v))/eps,LP_Z5(qfd,q,v),(PZ5-BZ5(v))/eps,&
              LP_SPHER1(qfd,q,v), (PSPHER1-BSPHER1(v))/eps,  LP_TRANS1(qfd,q,v), (PTRANS1-BTRANS1(v))/eps
        enddo
        write(1,*)' '
68      format(1x,i2,2x,i2,T19,3(1x,1p2e14.6,1x),2x,2(1x,1p2e13.5,1x))

!  End layer profile Jacobian loop

      enddo

!  Close file

      close(1)
      
!  Close error file if it has been opened

      write(*,*)
      IF ( OPENFILEFLAG ) then
        close( VLIDORT_ERRUNIT )
        write(*,*)'Main program executed with internal defaults,'// &
                  ' warnings in "V2p8p3_VLIDORT_Execution.log"'
      ELSE
        write(*,*)'Main program finished successfully'
      ENDIF

!  Finish

      stop
      end program Planetary_Tester

