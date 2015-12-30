      program V2OS_Full_LPCS_OMPTester

!  This is a Full tester for the combined linearized VLIDORT/2OS code
!   [ Uses the same setups and VLIDORT Inputs as the TwoOS_Self testers !!! ]

!  Module files for VLIDORT

      USE VLIDORT_PARS
      USE VLIDORT_IO_DEFS
      USE VLIDORT_LIN_IO_DEFS

      USE VLIDORT_AUX
      USE VLIDORT_L_INPUTS
      USE VLIDORT_LCS_MASTERS
      USE VLIDORT_LPS_MASTERS

      USE vlidort_2OScorr_master_m
      USE vlidort_2OScorr_lps_master_m
      USE vlidort_2OScorr_lcs_master_m
!mick fix 5/1/2015 - added new 2OS module
      USE vlidort_2OS_corr_aux

      IMPLICIT NONE

!  Dimensioning integer parameters

      INTEGER, PARAMETER :: MAXWVN = 4!50!100!1000

!  VLIDORT file inputs status structure

      TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn, VLIDORT_FixIn_SAVE
      TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn, VLIDORT_ModIn_SAVE

!  VLIDORT supplements i/o structure

      TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup, VLIDORT_Sup_SAVE

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

!  VLIDORT linearized input structures

      TYPE(VLIDORT_Fixed_LinInputs)          :: VLIDORT_LinFixIn, VLIDORT_LinFixIn_SAVE
      TYPE(VLIDORT_Modified_LinInputs)       :: VLIDORT_LinModIn, VLIDORT_LinModIn_SAVE

!  VLIDORT linearized supplements i/o structure

      TYPE(VLIDORT_LinSup_InOut)             :: VLIDORT_LinSup, VLIDORT_LinSup_SAVE

!  VLIDORT linearized output structure

      TYPE(VLIDORT_LinOutputs)               :: VLIDORT_LinOut

!  Local Variables
!  ===============

!  Dump variables

      integer          :: nlay,nmug,nmoms,nstokes,nfoumax,ncoefs(maxlayers)
      logical          :: regular_ps,enhanced_ps
      double precision :: epsilon,phi,emu0,emu,spars,xmu(maxstreams),w(maxstreams),&
                          Sun_Chapman2(maxlayers,maxlayers),omega
      double precision :: ssa(maxlayers),opd(maxlayers),coefs(0:maxmoments_input,maxlayers,6),&
                          eradius,heights(0:maxlayers)

!  Flag for opening error output file

      LOGICAL ::          OPENFILEFLAG, PART1, PART2

!  local

      integer :: nlayers, n_totalcolumn_wfs, n_totalprofile_wfs, n_surface_wfs      
      integer :: n, n1, k, l, nd, kd, ld, t, task, ks, nstokes_saved, o1, v, q
      double precision :: eps, epsfac, gaspert, scat, taupert, coalb
      double precision :: ssa_save(maxlayers),opd_save(maxlayers),spars_save

!  Saved Baseline results

      DOUBLE PRECISION :: STOKES_V_BAS ( MAX_GEOMETRIES, MAXSTOKES, MAXWVN )
      DOUBLE PRECISION :: STOKES_S_BAS ( MAX_GEOMETRIES, MAXSTOKES, MAXWVN )
      DOUBLE PRECISION :: STOKES_C_BAS ( MAX_GEOMETRIES, MAXSTOKES, MAXWVN )
!mick mod 5/1/2015 - add Icorr
      DOUBLE PRECISION :: Icorr ( MAX_GEOMETRIES, MAXWVN )

      DOUBLE PRECISION :: PROFILEWF_V_BAS ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, MAXSTOKES, MAXWVN )
      DOUBLE PRECISION :: PROFILEWF_S_BAS ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, MAXSTOKES, MAXWVN )
      DOUBLE PRECISION :: PROFILEWF_C_BAS ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, MAXSTOKES, MAXWVN )

      DOUBLE PRECISION :: SURFACEWF_V_BAS ( MAX_SURFACEWFS, MAX_GEOMETRIES, MAXSTOKES, MAXWVN )
      DOUBLE PRECISION :: SURFACEWF_S_BAS ( MAX_SURFACEWFS, MAX_GEOMETRIES, MAXSTOKES, MAXWVN )
      DOUBLE PRECISION :: SURFACEWF_C_BAS ( MAX_SURFACEWFS, MAX_GEOMETRIES, MAXSTOKES, MAXWVN )

      DOUBLE PRECISION :: COLUMNWF_V_BAS ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES, MAXWVN )
      DOUBLE PRECISION :: COLUMNWF_S_BAS ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES, MAXWVN )
      DOUBLE PRECISION :: COLUMNWF_C_BAS ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES, MAXWVN )

!  Saved Perturbation results. Note the nubmer of tasks !

      DOUBLE PRECISION :: STOKES_V_PT ( MAX_GEOMETRIES, MAXSTOKES, 0:70, MAXWVN )
      DOUBLE PRECISION :: STOKES_S_PT ( MAX_GEOMETRIES, MAXSTOKES, 0:70, MAXWVN )
      DOUBLE PRECISION :: STOKES_C_PT ( MAX_GEOMETRIES, MAXSTOKES, 0:70, MAXWVN )

!  Local shared array

!mick fix 5/1/2015 - added VLIDORT_STATUS_INPUTCHECK
      INTEGER :: VLIDORT_STATUS_INPUTCHECK
      INTEGER, DIMENSION ( MAXWVN ) :: VLIDORT_TID_SAVE

!  OpenMP tests (general)

      INTEGER       :: TEST, N_OMP_TESTS, OMP_MAXTHREADS
      INTEGER       :: TID, OMP_NTHREADS
      INTEGER       :: wvn, numwvn, ompt

      CHARACTER (LEN=1)  :: charTID
      CHARACTER (LEN=1)  :: CH1 
      CHARACTER (LEN=5)  :: CH2,CH3
      CHARACTER (LEN=30) :: TAIL

!  OpenMP tests (timing)

      INTEGER :: n_core, time_divider
      REAL    :: omp_e1, omp_e2

      LOGICAL :: do_vlid_ns3_timecheck_only
      LOGICAL :: do_vlid_ns1_v2os_timecheck_only

!  OpenMP functions

      INTEGER :: OMP_GET_NUM_THREADS, &
                 OMP_GET_THREAD_NUM

!  Debug

      LOGICAL, parameter :: debug = .true.
      !LOGICAL, parameter :: debug = .false.

!  Flag for standard output file

      !LOGICAL, parameter :: do_std_output_file = .true.
      LOGICAL, parameter :: do_std_output_file = .false.

!  Start of code
!  =============

!  Define some control variables:

!  (1) test control
!      * LPS - part 1
!      * LCS - part 2

!      PART1 = .FALSE.
      PART1 = .TRUE.

!      PART2 = .FALSE.
      PART2 = .TRUE.

!  (2) individual timing check control
!      * NS3     - VLIDORT nstokes=3 time check only (no output files created)
!      * NS1 2OS - VLIDORT nstokes=1 with V2OS code time check only (no output files created)

!      Mick note: The default settings for both of these flags is FALSE for a normal lpcs test.
!                 Only turn on ONE at a time to do either a VLID-NS3 or a VLID-NS1/V2OS time check

      !do_vlid_ns3_timecheck_only = .true.
      do_vlid_ns3_timecheck_only = .false.

      !do_vlid_ns1_v2os_timecheck_only = .true.
      do_vlid_ns1_v2os_timecheck_only = .false.

!  Set test parameters

      N_CORE = 2!4
      N_OMP_TESTS = 2!3

      numwvn = maxwvn
!      ntasks = maxtasks

!  Echo test setup

      write(*,*)
      if ( do_vlid_ns3_timecheck_only ) then
        write(*,*) 'Doing VLID-NS3 time check (no output files created)'
      else if ( do_vlid_ns1_v2os_timecheck_only ) then
        write(*,*) 'Doing VLID-NS1/V2OS time check (no output files created)'
      else if ( .not.do_vlid_ns3_timecheck_only .and. .not.do_vlid_ns1_v2os_timecheck_only ) then
        write(*,*) 'Doing normal lpcs test'
      else if ( do_vlid_ns3_timecheck_only .and. do_vlid_ns1_v2os_timecheck_only ) then
        write(*,*) 'Error: setup to perform ONE time check at a time only'
        stop
      end if

!  Begin test loop

      DO TEST = 1, N_OMP_TESTS
      !DO TEST = 2, N_OMP_TESTS
      !DO TEST = 1, 1
      !DO TEST = 2, 2

      write(*,*)
      write(*,*) '******************************************'
      write(*,'(1x,a,i1)') 'Doing test ',TEST

!  Set total number of tasks (optional)
!      Note: if not set, OpenMP will set the number of threads
!            to the number of cores on the system
!      Note: if not set, OpenMP will default to dividing up the
!            wvn loop below equally among the available threads

      if ( test == 1 ) then
        OMP_MAXTHREADS = 1
      elseif ( test == 2 ) then
        OMP_MAXTHREADS = 2
      elseif ( test == 3 ) then
        OMP_MAXTHREADS = 4
      endif

      !write(*,*)
      !write(*,*) 'OMP_MAXTHREADS = ',OMP_MAXTHREADS

!  Set total number of threads

      !write(*,*)
      !write(*,*) 'Setting total number of threads'
      CALL OMP_SET_NUM_THREADS(OMP_MAXTHREADS)

!  Fixed Boolean inputs
!  --------------------

      VLIDORT_FixIn_Save%Bool%TS_DO_FULLRAD_MODE      = .true.    ! V2OS Full - Must be TRUE for the Full test
      VLIDORT_FixIn_Save%Bool%TS_DO_SSCORR_TRUNCATION = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_SS_EXTERNAL       = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_SSFULL            = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_THERMAL_EMISSION  = .false.!.true.
      VLIDORT_FixIn_Save%Bool%TS_DO_SURFACE_EMISSION  = .false.!.true.

      VLIDORT_FixIn_Save%Bool%TS_DO_PLANE_PARALLEL      = .false.  ! Set later
      VLIDORT_FixIn_Save%Bool%TS_DO_LAMBERTIAN_SURFACE  = .true.

      VLIDORT_FixIn_Save%Bool%TS_DO_UPWELLING           = .true.
      VLIDORT_FixIn_Save%Bool%TS_DO_DNWELLING           = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_SURFACE_LEAVING     = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_SL_ISOTROPIC        = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_QUAD_OUTPUT         = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_TOA_CONTRIBS        = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_SPECIALIST_OPTION_1 = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_SPECIALIST_OPTION_2 = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_SPECIALIST_OPTION_3 = .false.

!  Fixed control
!  -------------

      VLIDORT_FixIn_Save%Cont%TS_NSTOKES          = 0 ! Set later
      VLIDORT_FixIn_Save%Cont%TS_NSTREAMS         = 0 ! set later
      VLIDORT_FixIn_Save%Cont%TS_NLAYERS          = 0 ! set later
      VLIDORT_FixIn_Save%Cont%TS_NFINELAYERS      = 2 ! V2OS Full - Should be set for the Full test
      VLIDORT_FixIn_Save%Cont%TS_N_THERMAL_COEFFS = 0!2 ! V2OS Full - Should be set for the Full test

      VLIDORT_FixIn_Save%Cont%TS_TAYLOR_ORDER     = 3 ! V2OS Full - Should be set for the Full test

      VLIDORT_FixIn_Save%Cont%TS_VLIDORT_ACCURACY  = zero ! Set later

      VLIDORT_FixIn_Save%Cont%TS_NLAYERS_NOMS   = 0
      VLIDORT_FixIn_Save%Cont%TS_NLAYERS_CUTOFF = 0

!  Fixed Beam Userval inputs
!  --------------------------

      VLIDORT_FixIn_Save%Sunrays%TS_FLUX_FACTOR     = ONE
      VLIDORT_FixIn_Save%UserVal%TS_N_USER_LEVELS   = 1

!  Chapman inputs
!  --------------

      VLIDORT_FixIn_Save%Chapman%TS_HEIGHT_GRID      = ZERO ! V2OS Full - must be assigned from input (see later)

      VLIDORT_FixIn_Save%Chapman%TS_PRESSURE_GRID    = ZERO
      VLIDORT_FixIn_Save%Chapman%TS_TEMPERATURE_GRID = ZERO
      VLIDORT_FixIn_Save%Chapman%TS_FINEGRID         = 0
      VLIDORT_FixIn_Save%Chapman%TS_RFINDEX_PARAMETER = ZERO

!  Fixed Optical inputs
!  --------------------

      VLIDORT_FixIn_Save%Optical%TS_THERMAL_BB_INPUT = ZERO
      VLIDORT_FixIn_Save%Optical%TS_SURFACE_BB_INPUT = ZERO

      VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT = zero ! Set later
      VLIDORT_FixIn_Save%Optical%TS_DELTAU_VERT_INPUT    = zero ! Set later
      VLIDORT_FixIn_Save%Optical%TS_LAMBERTIAN_ALBEDO    = zero ! Set later

!  Fixed write
!  -----------

      VLIDORT_FixIn_Save%Write%TS_DO_DEBUG_WRITE    = .false.
      VLIDORT_FixIn_Save%Write%TS_DO_WRITE_INPUT    = .false.
      VLIDORT_FixIn_Save%Write%TS_DO_WRITE_SCENARIO = .false.
      VLIDORT_FixIn_Save%Write%TS_DO_WRITE_FOURIER  = .false.
      VLIDORT_FixIn_Save%Write%TS_DO_WRITE_RESULTS  = .false.

      VLIDORT_FixIn_Save%Write%TS_INPUT_WRITE_FILENAME     = ' '
      VLIDORT_FixIn_Save%Write%TS_SCENARIO_WRITE_FILENAME  = ' '
      VLIDORT_FixIn_Save%Write%TS_FOURIER_WRITE_FILENAME   = ' '
      VLIDORT_FixIn_Save%Write%TS_RESULTS_WRITE_FILENAME   = ' '

!  Modified Booleans
!  -----------------

      VLIDORT_ModIn_Save%MBool%TS_DO_SOLAR_SOURCES   = .true.

!  Fixed defaults
!    -- Full calculation with outgoing sphericity single-scatter

      VLIDORT_ModIn_Save%MBool%TS_DO_SSCORR_NADIR        = .false. !  Set later
      VLIDORT_ModIn_Save%MBool%TS_DO_SSCORR_OUTGOING     = .false. !  Set later
!mick fix 5/1/2015 - added DO_FO_CALC initialization
      VLIDORT_ModIn_Save%MBool%TS_DO_FO_CALC             = .false.

      VLIDORT_ModIn_Save%MBool%TS_DO_REFRACTIVE_GEOMETRY = .false.
      VLIDORT_ModIn_Save%MBool%TS_DO_CHAPMAN_FUNCTION    = .true.  !  V2OS Full : Must be set now

      VLIDORT_ModIn_Save%MBool%TS_DO_USER_VZANGLES       = .true.
      VLIDORT_ModIn_Save%MBool%TS_DO_ADDITIONAL_MVOUT    = .false.
      VLIDORT_ModIn_Save%MBool%TS_DO_MVOUT_ONLY          = .false.
      VLIDORT_ModIn_Save%MBool%TS_DO_THERMAL_TRANSONLY   = .false.

      VLIDORT_ModIn_Save%MBool%TS_DO_SOLUTION_SAVING     = .false.
      VLIDORT_ModIn_Save%MBool%TS_DO_BVP_TELESCOPING     = .false.

!  From the input (Observation Geometry is new)

      VLIDORT_ModIn_Save%MBool%TS_DO_OBSERVATION_GEOMETRY = .true.

      VLIDORT_ModIn_Save%MBool%TS_DO_RAYLEIGH_ONLY        = .false.
      VLIDORT_ModIn_Save%MBool%TS_DO_DOUBLE_CONVTEST      = .true.
      VLIDORT_ModIn_Save%MBool%TS_DO_DELTAM_SCALING       = .false. ! Status uncertain

      VLIDORT_ModIn_Save%MBool%TS_DO_2OS_CORRECTION       = .false. ! Will be set later in the Full tests

!  Modified Control
!  ----------------

      VLIDORT_ModIn_Save%MCont%TS_NGREEK_MOMENTS_INPUT = 0 ! Set later

!  Modified Sunrays
!  ----------------

      VLIDORT_ModIn_Save%MSunrays%TS_N_SZANGLES       = 1
      VLIDORT_ModIn_Save%MSunrays%TS_SZANGLES         = zero ! set later

!  Modified Uservals
!  -----------------

      VLIDORT_ModIn_Save%MUserVal%TS_N_USER_VZANGLES      = 1
      VLIDORT_ModIn_Save%MUserVal%TS_USER_VZANGLES_INPUT  = zero ! set later
      VLIDORT_ModIn_Save%MUserVal%TS_N_USER_RELAZMS       = 1
      VLIDORT_ModIn_Save%MUserVal%TS_USER_RELAZMS         = zero  ! Set later
!mick fix 4/21/2015 - initialize all elements of USER_LEVELS
      !VLIDORT_ModIn_Save%MUserVal%TS_USER_LEVELS(1)       = zero
      VLIDORT_ModIn_Save%MUserVal%TS_USER_LEVELS          = zero
      VLIDORT_ModIn_Save%MUserVal%TS_GEOMETRY_SPECHEIGHT  = zero
      VLIDORT_ModIn_Save%MUserVal%TS_N_USER_OBSGEOMS      = 1
      VLIDORT_ModIn_Save%MUserVal%TS_USER_OBSGEOMS_INPUT  = zero !  set later

!  Modified Chapman inputs
!  -----------------------

      VLIDORT_ModIn_Save%MChapman%TS_EARTH_RADIUS         = zero ! V2OS_Full : Now assiged from input file. See later
      VLIDORT_ModIn_Save%MChapman%TS_CHAPMAN_FACTORS      = zero ! V2OS_Full : Now Calculated Internally

!  Modified Optical
!  ----------------

      VLIDORT_ModIn_Save%MOptical%TS_OMEGA_TOTAL_INPUT    = zero ! Set later

!  Linearized
!  ----------

      VLIDORT_LinFixIn_Save%Cont%TS_LAYER_VARY_FLAG   = .false.  ! Set later
      VLIDORT_LinFixIn_Save%Cont%TS_LAYER_VARY_NUMBER = 0        ! Set later
!mick fix 5/1/2015 - added these initializations
      VLIDORT_LinFixIn_Save%Cont%TS_N_TOTALCOLUMN_WFS  = 0       ! Set later
      VLIDORT_LinFixIn_Save%Cont%TS_N_TOTALPROFILE_WFS = 0       ! Set later
      VLIDORT_LinFixIn_Save%Cont%TS_N_SURFACE_WFS      = 0       ! Set later
      VLIDORT_LinFixIn_Save%Cont%TS_N_SLEAVE_WFS       = 0       ! Set later
      VLIDORT_LinFixIn_Save%Cont%TS_COLUMNWF_NAMES     = ' '
      VLIDORT_LinFixIn_Save%Cont%TS_PROFILEWF_NAMES    = ' '

      VLIDORT_LinFixIn_Save%Optical%TS_L_deltau_vert_input    = zero ! Set later
      VLIDORT_LinFixIn_Save%Optical%TS_L_greekmat_total_input = zero
      VLIDORT_LinFixIn_Save%Optical%TS_L_omega_total_input    = zero ! Set later

!mick fix 5/1/2015 - added these initializations
      VLIDORT_LinModIn_Save%MCont%TS_DO_COLUMN_LINEARIZATION  = .false. ! Set later
      VLIDORT_LinModIn_Save%MCont%TS_DO_PROFILE_LINEARIZATION = .false. ! Set later
      VLIDORT_LinModIn_Save%MCont%TS_DO_ATMOS_LINEARIZATION   = .false. ! Set later
      VLIDORT_LinModIn_Save%MCont%TS_DO_SURFACE_LINEARIZATION = .false. ! Set later
      VLIDORT_LinModIn_Save%MCont%TS_DO_LINEARIZATION         = .false.
      VLIDORT_LinModIn_Save%MCont%TS_DO_SIMULATION_ONLY       = .false.
      VLIDORT_LinModIn_Save%MCont%TS_DO_ATMOS_LBBF            = .false.
      VLIDORT_LinModIn_Save%MCont%TS_DO_SURFACE_LBBF          = .false.
      VLIDORT_LinModIn_Save%MCont%TS_DO_SLEAVE_WFS            = .false.

! Supplemental Initialization
! ---------------------------

!  BRDF contributions, may be set later, in CODEC_VLIDORT_Optical

      VLIDORT_Sup_Save%BRDF%TS_EXACTDB_BRDFUNC  = ZERO
      VLIDORT_Sup_Save%BRDF%TS_BRDF_F_0         = ZERO
      VLIDORT_Sup_Save%BRDF%TS_BRDF_F           = ZERO
      VLIDORT_Sup_Save%BRDF%TS_USER_BRDF_F_0    = ZERO
      VLIDORT_Sup_Save%BRDF%TS_USER_BRDF_F      = ZERO

!  Emissitivity may be set later, in CODEC_VLIDORT_Optical

      VLIDORT_Sup_Save%BRDF%TS_EMISSIVITY       = ZERO
      VLIDORT_Sup_Save%BRDF%TS_USER_EMISSIVITY  = ZERO

!  Surface leaving stuff may be set later, in CODEC_VLIDORT_Optical

      VLIDORT_Sup_Save%SLEAVE%TS_SLTERM_ISOTROPIC  = ZERO
      VLIDORT_Sup_Save%SLEAVE%TS_SLTERM_USERANGLES = ZERO
      VLIDORT_Sup_Save%SLEAVE%TS_SLTERM_F_0        = ZERO
      VLIDORT_Sup_Save%SLEAVE%TS_USER_SLTERM_F_0   = ZERO

!mick fix 5/1/2015 - added SUP initializations

      VLIDORT_Sup_Save%SS%TS_STOKES_SS  = ZERO
      VLIDORT_Sup_Save%SS%TS_STOKES_DB  = ZERO
      VLIDORT_Sup_Save%SS%TS_FO_ZMATRIX = ZERO
      VLIDORT_Sup_Save%SS%TS_FO_R1SAVED = ZERO

!  Linearized (not required)

      VLIDORT_LinSup_Save%BRDF%TS_LS_EXACTDB_BRDFUNC = ZERO
      VLIDORT_LinSup_Save%BRDF%TS_LS_BRDF_F_0        = ZERO
      VLIDORT_LinSup_Save%BRDF%TS_LS_BRDF_F          = ZERO
      VLIDORT_LinSup_Save%BRDF%TS_LS_USER_BRDF_F_0   = ZERO
      VLIDORT_LinSup_Save%BRDF%TS_LS_USER_BRDF_F     = ZERO
      VLIDORT_LinSup_Save%BRDF%TS_LS_EMISSIVITY      = ZERO
      VLIDORT_LinSup_Save%BRDF%TS_LS_USER_EMISSIVITY = ZERO

      VLIDORT_LinModIn_Save%MCont%TS_DO_SLEAVE_WFS = .FALSE.
      VLIDORT_LinFixIn_Save%Cont%TS_N_SLEAVE_WFS   = 0
   
      VLIDORT_LinSup_Save%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC  = ZERO
      VLIDORT_LinSup_Save%SLEAVE%TS_LSSL_SLTERM_USERANGLES = ZERO
      VLIDORT_LinSup_Save%SLEAVE%TS_LSSL_SLTERM_F_0        = ZERO
      VLIDORT_LinSup_Save%SLEAVE%TS_LSSL_USER_SLTERM_F_0   = ZERO

!mick fix 5/1/2015 - added linearized SUP initializations
      VLIDORT_LinSup_Save%SS%Prof%TS_PROFILEWF_SS = ZERO
      VLIDORT_LinSup_Save%SS%Prof%TS_PROFILEWF_DB = ZERO
      VLIDORT_LinSup_Save%SS%Col%TS_COLUMNWF_SS   = ZERO
      VLIDORT_LinSup_Save%SS%Col%TS_COLUMNWF_DB   = ZERO
      VLIDORT_LinSup_Save%SS%Surf%TS_SURFACEWF_DB = ZERO

!  Now open the files
!  ------------------

      open(87,file='fort.87_Rg_h',status='old')
      read(87,*)nlay,nmug,nstokes,nmoms,nfoumax,regular_ps,enhanced_ps
      read(87,*)epsilon,phi,emu0,emu,spars,eradius,heights(0) ; spars_save = spars
      do k = 1, nmug
         read(87,*)kd,xmu(k),w(k) ! Not used
      enddo
      do n = 1, nlay
         read(87,*)heights(n),ssa(n),opd(n),ncoefs(n)
         opd_save(n) = opd(n) ; ssa_save(n) = ssa(n)
         do l = 0, 2*nmug-1
            read(87,*)nd,ld,(coefs(l,n,k),k=1,6)
         enddo
      enddo
      close(87)

!      Sun_Chapman2 = 0.0d0
!      open(88,file='fort.88',status='old')
!         do n = 1, nlay
!            do k = 1, n
!               read(88,*)nd,kd,Sun_Chapman2(n,k)
!            enddo
!         enddo
!      close(88)

!  Now set the Remaining VLIDORT inputs
!  ------------------------------------

!  Control

      VLIDORT_ModIn_Save%MBool%TS_DO_SSCORR_NADIR        = regular_ps
      VLIDORT_ModIn_Save%MBool%TS_DO_SSCORR_OUTGOING     = enhanced_ps
      VLIDORT_FixIn_Save%Bool%TS_DO_PLANE_PARALLEL       = .not. regular_ps .and. .not. enhanced_ps

      VLIDORT_FixIn_Save%Cont%TS_NSTOKES          = nstokes
      VLIDORT_FixIn_Save%Cont%TS_NSTREAMS         = nmug
      VLIDORT_FixIn_Save%Cont%TS_NLAYERS          = nlay
      VLIDORT_FixIn_Save%Cont%TS_VLIDORT_ACCURACY = epsilon

      VLIDORT_ModIn_Save%MCont%TS_NGREEK_MOMENTS_INPUT = nmoms

!  Geometry

      VLIDORT_ModIn_Save%MSunrays%TS_SZANGLES             = acos(emu0)/deg_to_rad
      VLIDORT_ModIn_Save%MUserVal%TS_USER_VZANGLES_INPUT  = acos(emu) /deg_to_rad
      VLIDORT_ModIn_Save%MUserVal%TS_USER_RELAZMS         = phi

      VLIDORT_ModIn_Save%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1) = acos(emu0)/deg_to_rad
      VLIDORT_ModIn_Save%MUserVal%TS_USER_OBSGEOMS_INPUT(1,2) = acos(emu) /deg_to_rad
      VLIDORT_ModIn_Save%MUserVal%TS_USER_OBSGEOMS_INPUT(1,3) = phi

!  Chapman (heights and eradius)
!    Factors are no longer copied from the file-read, Now Internal calculation!
!        VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS(:,:,1)   = Sun_Chapman2(:,:)

      VLIDORT_ModIn_Save%MChapman%TS_EARTH_RADIUS     = eradius
      do n = 0, nlay
        VLIDORT_FixIn_Save%Chapman%TS_HEIGHT_GRID(n)  = heights(n)
      enddo

!  Optical

      do n = 1, nlay
        n1 = nlay + 1 - n
        VLIDORT_FixIn_Save%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
        VLIDORT_ModIn_Save%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,1)  =   coefs(0:nmoms,n,1)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,2)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,5)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,6)  =   coefs(0:nmoms,n,2)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,11) =   coefs(0:nmoms,n,3)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,12) =   coefs(0:nmoms,n,6)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,15) = - coefs(0:nmoms,n,6)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,16) =   coefs(0:nmoms,n,4)
      enddo
      VLIDORT_FixIn_Save%Optical%TS_LAMBERTIAN_ALBEDO =   spars_save

!  Skip part 1 of the tests if desired

      IF ( .NOT. PART1 ) GO TO 567

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!                     P A R T    I

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      write(*,*)
      write(*,*)'Doing PROFILEWF part: '// &
        'Intensities + 60 profile Jacobians + 1 Surface Jacobian'

!  Start timing

      call cpu_time(omp_e1)

!     Only master thread does this

      !OMP_NTHREADS = OMP_GET_NUM_THREADS()
      !write(*,*)
      !write(*,'(1x,a,i1)') 'Total number of threads (to start) = ', OMP_NTHREADS

!     Initialize output

      STOKES_V_BAS = ZERO
      STOKES_S_BAS = ZERO
      STOKES_C_BAS = ZERO
      Icorr = ZERO

      PROFILEWF_V_BAS = ZERO
      PROFILEWF_S_BAS = ZERO
      PROFILEWF_C_BAS = ZERO

      SURFACEWF_V_BAS = ZERO
      SURFACEWF_S_BAS = ZERO
      SURFACEWF_C_BAS = ZERO

      STOKES_V_PT  = ZERO
      STOKES_S_PT  = ZERO
      STOKES_C_PT  = ZERO

!  Initialize error file output flag

      OPENFILEFLAG = .false.

!  Initialize shared input check flag

      VLIDORT_STATUS_INPUTCHECK = 0

!     Begin parallel region

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (OMP_NTHREADS, &
!$OMP     DO_VLID_NS3_TIMECHECK_ONLY, DO_VLID_NS1_V2OS_TIMECHECK_ONLY, &
!$OMP     NUMWVN, NLAY, &
!$OMP     opd_save, ssa_save, spars_save, eps, &
!$OMP     VLIDORT_FixIn_SAVE,    VLIDORT_ModIn_SAVE,    VLIDORT_Sup_SAVE, &
!$OMP     VLIDORT_LinFixIn_SAVE, VLIDORT_LinModIn_SAVE, VLIDORT_LinSup_SAVE, &
!$OMP     STOKES_V_BAS, STOKES_S_BAS, STOKES_C_BAS, Icorr, &
!$OMP     PROFILEWF_V_BAS, PROFILEWF_S_BAS, PROFILEWF_C_BAS, &
!$OMP     SURFACEWF_V_BAS, SURFACEWF_S_BAS, SURFACEWF_C_BAS, &
!$OMP     STOKES_V_PT, STOKES_S_PT, STOKES_C_PT, &
!$OMP     VLIDORT_STATUS_INPUTCHECK, VLIDORT_TID_SAVE ) &
!$OMP   FIRSTPRIVATE (OPENFILEFLAG)

!     Obtain thread number

      TID = OMP_GET_THREAD_NUM()
      !TID = 0

!     Obtain and display total number of threads and local thread

      IF (TID == 0) THEN
        OMP_NTHREADS = OMP_GET_NUM_THREADS()
        !OMP_NTHREADS = 1
        !write(*,*)
        !write(*,'(1x,a,i1)') 'Total number of threads (inside parallel region) = ', OMP_NTHREADS
        write(*,*)
        write(*,'(1x,a,i1,a)') 'Running V2OS driver with ',OMP_NTHREADS,' thread(s)'
      END IF

      write(*,'(1x,a,i1,a)') 'Thread ',TID,' beginning its wavenumber loop'

!$OMP DO

!  Begin "wavenumber" loop

      do wvn = 1, numwvn
      !do wvn = 1, 1

        !write(*,'(1x,a,i1,a,i5)') 'Thread  = ',TID,' wvn = ',wvn

!  Initialize some inputs

        VLIDORT_FixIn = VLIDORT_FixIn_SAVE
        VLIDORT_ModIn = VLIDORT_ModIn_SAVE
        VLIDORT_Sup   = VLIDORT_Sup_SAVE

        VLIDORT_LinFixIn = VLIDORT_LinFixIn_SAVE
        VLIDORT_LinModIn = VLIDORT_LinModIn_SAVE
        VLIDORT_LinSup   = VLIDORT_LinSup_SAVE

!  linearization control

        VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = .false.
        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAY)   = .true.
        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAY) = 1
        VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = 1
        VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS         = 0
        VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 1

        VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = .false.
        VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = .true.
        VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .true.
        VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .true.
!mick fix 5/1/2015 - add control variable
        VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = .true.

!  Additional proxies

        nlayers            = nlay
        n_surface_wfs      = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS
        n_totalprofile_wfs = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
        n_totalcolumn_wfs  = VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS

!  Linearization inputs. 1 Profile WF, WRT Gas absorption

        do n = 1, nlay
          n1 = nlay + 1 - n ; coalb = one - ssa_save(n)
          VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(1,n1)  =  coalb
          VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(1,n1)  = -coalb
        enddo

!  Precautions for VLIDORT

        do n = 1, nlay
           VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0,n,1) = one
           omega = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)
           if ( omega.gt.0.99999999d0 ) VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n) = 0.99999999d0
        enddo

!  First Call. VECTOR VLIDORT
!  ==========================

!  Set conditions

        nstokes_saved = VLIDORT_FixIn%Cont%TS_nstokes 
        VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .false.

!  Begin VLID-NS3 control if block - turn off VLID-NS3 only during VLID-NS1/V20S time check
!  mick note: if block started here as the above two conditions need set for VLID-NS1/V2OS also

        if (.not.do_vlid_ns1_v2os_timecheck_only) then

!  VLIDORT call

          CALL VLIDORT_LPS_MASTER ( &
            VLIDORT_FixIn, &
            VLIDORT_ModIn, &
            VLIDORT_Sup, &
            VLIDORT_Out, &
            VLIDORT_LinFixIn, &
            VLIDORT_LinModIn, &
            VLIDORT_LinSup, &
            VLIDORT_LinOut )

          if ((debug).and.(TID == 0)) then
            !write(*,*)
            write(*,*)'Done Vector VLIDORT, Baseline LPS calculation'
          end if

!  Exception handling, write-up (optional)

          write(charTID,'(i1)') TID
          CALL VLIDORT_WRITE_STATUS ( &
              'V2OS_full_LPS_OMPExecution_TID_' // charTID // '.log', &
              35+TID, OPENFILEFLAG, &
              VLIDORT_Out%Status )

          IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
               VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save baseline results

          do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
             do v = 1, VLIDORT_Out%Main%TS_n_geometries
                stokes_V_BAS(v,o1,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
                do q = 1, n_totalprofile_wfs
                   profilewf_V_BAS(q,1:nlayers,v,o1,wvn) = &
                     VLIDORT_LinOut%Prof%TS_profilewf(q,1:nlayers,1,v,o1,1)
                enddo
                do q = 1, n_surface_wfs
                   surfacewf_V_BAS(q,v,o1,wvn) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,1)
                enddo
             enddo
          enddo

!  End VLID-NS3 control if block

        endif


!  Begin VLID-NS1/V20S control if block -
!       turn off VLID-NS1/V20S only during VLID-NS3 time check

        if (.not.do_vlid_ns3_timecheck_only) then

!  Second Call. SCALAR VLIDORT + 2OS CORRECTION
!  ============================================

!  Set conditions

          VLIDORT_FixIn%Cont%TS_nstokes            = 1
          VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .true.

!  VLIDORT call

          CALL VLIDORT_LPS_MASTER ( &
            VLIDORT_FixIn, &
            VLIDORT_ModIn, &
            VLIDORT_Sup, &
            VLIDORT_Out, &
            VLIDORT_LinFixIn, &
            VLIDORT_LinModIn, &
            VLIDORT_LinSup, &
            VLIDORT_LinOut )

          if ((debug).and.(TID == 0)) then
            !write(*,*)
            write(*,*)'Done Scalar VLIDORT, Baseline LPS calculation'
          end if

!  Exception handling, write-up (optional)

          CALL VLIDORT_WRITE_STATUS ( &
              'V2OS_full_LPS_OMPExecution_TID_' // charTID // '.log', &
              35+TID, OPENFILEFLAG, &
              VLIDORT_Out%Status )

          IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
               VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

          do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
             do v = 1, VLIDORT_Out%Main%TS_n_geometries
                stokes_S_BAS(v,o1,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
                do q = 1, n_totalprofile_wfs
                   profilewf_S_BAS(q,1:nlayers,v,o1,wvn) = &
                     VLIDORT_LinOut%Prof%TS_profilewf(q,1:nlayers,1,v,o1,1)
                enddo
                do q = 1, n_surface_wfs
                  surfacewf_S_BAS(q,v,o1,wvn) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,1)
                enddo
             enddo
          enddo

!  Perform 2OS correction
!  ----------------------

!  Set conditions

          VLIDORT_FixIn%Cont%TS_nstokes = nstokes_saved

!  Call 2OS correction routine

          CALL VLIDORT_2OSCORR_LPS_MASTER ( &
            VLIDORT_FixIn, &
            VLIDORT_ModIn, &
            VLIDORT_Sup, &
            VLIDORT_Out, &
            VLIDORT_LinFixIn, &
            VLIDORT_LinModIn, &
            VLIDORT_LinSup, &
            VLIDORT_LinOut )

          if ((debug).and.(TID == 0)) then
            !write(*,*)
            write(*,*)'Done 2OS correction, Baseline LPS calculation'
          end if

!  Exception handling, write-up (optional)

!mick fix 5/1/2015 - created new status subroutine for 2OS code 
          CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
            'V2OS_full_LPS_OMPExecution_TID_' // charTID // '.log', &
            35+TID, OPENFILEFLAG, &
            VLIDORT_Out%Status )

          IF ( VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
               VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

          do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
             do v = 1, VLIDORT_Out%Main%TS_n_geometries
                stokes_C_BAS(v,o1,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
                if (o1 == 1) Icorr(v,wvn) = VLIDORT_Out%Main%TS_2OSCORR_ICORR(v)
                do q = 1, n_totalprofile_wfs
                   profilewf_C_BAS(q,1:nlayers,v,o1,wvn) = &
                    VLIDORT_LinOut%Prof%TS_profilewf(q,1:nlayers,1,v,o1,upidx)
                enddo
                do q = 1, n_surface_wfs
                  surfacewf_C_BAS(q,v,o1,wvn) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,upidx)
                enddo
             enddo
          enddo

!  Write results

!        write(*,*)'LPS,R2(1) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,1),VLIDORT_Out%Main%TS_stokes(1,1,1,upidx)
!        write(*,*)'LPS,R2(2) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,2),VLIDORT_Out%Main%TS_stokes(1,1,2,upidx)
!        write(*,*)'LPS,R2(3) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,3),VLIDORT_Out%Main%TS_stokes(1,1,3,upidx)
!        write(*,*)'LPS,Icorr = ',VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

!        write(*,*)'LPS,LS_R2(1) = ',VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,1)
!        write(*,*)'LPS,LS_R2(2) = ',VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,2)
!        write(*,*)'LPS,LS_R2(3) = ',VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,3)
!        write(*,*)'LPS,LS_Icorr = ',VLIDORT_LinOut%Surf%TS_2OSCORR_LS_Icorr(1,1)

!        ks = 36
!        write(*,*)'LPS 36,LP_R2(1) = ',VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(1,ks,1,1)
!        write(*,*)'LPS 36,LP_R2(2) = ',VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(1,ks,1,2)
!        write(*,*)'LPS 36,LP_R2(3) = ',VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(1,ks,1,3)
!        write(*,*)'LPS 36,LP_Icorr = ',VLIDORT_LinOut%Prof%TS_2OSCORR_LP_Icorr(1,ks,1)

!        stop'baseline profile WFS only'


!  End VLID-NS1/V20S control if block

        endif


!  Begin FD control if block - turn off FD during either VLID-NS3 or VLID-NS1/V20S time check

        if (.not.do_vlid_ns3_timecheck_only .and. .not.do_vlid_ns1_v2os_timecheck_only) then

!  Initialize FD section
!  =====================

!  Fd perturbation

          eps = 1.0d-02
          epsfac = 1.0d0 + eps

!  Initialize linearized inputs

          VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = .true.
          VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAY)   = .false.
          VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAY) = 0
          VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = 0
          VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 0

          VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = .false.
          VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .false.
          VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .false.
!mick fix 5/1/2015 - add control variable
          VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = .false.

          VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT = zero
          VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT = zero

!  Additional proxies

          n_surface_wfs      = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS
          n_totalprofile_wfs = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS

!  TASKS for FD testing
!  ======================

!  Task 0, FD Perturbation on Surface Albedo
!  Task t, FD Perturbation on trace gas optical depth in layer t  (2OS = lowest layer)

          do task = 0, nlay
!        do task = 0,0

!  task counter

           t = task

!  Task 0 - surface jacobian

           if ( task .eq. 0 ) then
              if ((debug).and.(TID == 0)) &
                write(*,*)'Doing task 0: Finite Difference, perturb Lambertian albedo'
              VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =   spars_save * epsfac
              do n = 1, nlay
                 n1 = nlay + 1 - n
                 VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
                 VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
              enddo
           endif

!  Tasks 1 to NLAYERS - profile jacobians

           if ( task.gt.0  ) then
              if ((debug).and.(TID == 0)) &
                write(*,'(a,i2,a,i2)')' Doing task ',task, &
                  ': Finite Difference, perturb molecular absorption Layer ',t
              !if (task.eq.1) write(*,'(a,i2,a,i2)')' Doing tasks 1 -',nlay, &
              !   ': Finite Difference, perturb molecular absorption layers 1 -',nlay
              VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =   spars_save
              do n = 1, nlay
                 n1 = nlay + 1 - n
                 VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
                 VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
              enddo
              do n = 1, nlay
                 n1 = nlay + 1 - n
                 if ( n1.eq.t ) then
                    scat = opd_save(n) * ssa_save(n)
                    gaspert = ( opd_save(n) - scat ) * epsfac
                    taupert = scat + gaspert
                    VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = taupert
                    VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = scat / taupert
                 endif
              enddo
           endif

!  Precautions for VLIDORT

           do n = 1, nlay
              omega = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)
              if ( omega.gt.0.99999999d0 ) VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n) = 0.99999999d0
           enddo

!  First Call. VECTOR VLIDORT
!  ==========================

!  Set conditions

           nstokes_saved = VLIDORT_FixIn%Cont%TS_nstokes 
           VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .false.

!  VLIDORT call

           CALL VLIDORT_LPS_MASTER ( &
             VLIDORT_FixIn, &
             VLIDORT_ModIn, &
             VLIDORT_Sup, &
             VLIDORT_Out, &
             VLIDORT_LinFixIn, &
             VLIDORT_LinModIn, &
             VLIDORT_LinSup, &
             VLIDORT_LinOut )

           if ((debug).and.(TID == 0)) then
             !write(*,*)
             write(*,*)'Done Vector VLIDORT for FD task #', t
           end if

!  Exception handling, write-up (optional)

           CALL VLIDORT_WRITE_STATUS ( &
             'V2OS_full_LPS_OMPExecution_TID_' // charTID // '.log', &
             35+TID, OPENFILEFLAG, &
             VLIDORT_Out%Status )

           IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
                VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

           do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
             do v = 1, VLIDORT_Out%Main%TS_n_geometries
               stokes_V_PT(v,o1,t,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
             enddo
           enddo

!  Second Call. SCALAR VLIDORT + 2OS CORRECTION
!  ============================================

!  Set conditions

           VLIDORT_FixIn%Cont%TS_nstokes            = 1
           VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .true.

!  VLIDORT call

           CALL VLIDORT_LPS_MASTER ( &
             VLIDORT_FixIn, &
             VLIDORT_ModIn, &
             VLIDORT_Sup, &
             VLIDORT_Out, &
             VLIDORT_LinFixIn, &
             VLIDORT_LinModIn, &
             VLIDORT_LinSup, &
             VLIDORT_LinOut )

           if ((debug).and.(TID == 0)) then
             !write(*,*)
             write(*,*)'Done Scalar VLIDORT for FD task #', t
           end if

!  Exception handling, write-up (optional)

           CALL VLIDORT_WRITE_STATUS ( &
             'V2OS_full_LPS_OMPExecution_TID_' // charTID // '.log', &
             35+TID, OPENFILEFLAG, &
             VLIDORT_Out%Status )

           IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
                VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

           do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
             do v = 1, VLIDORT_Out%Main%TS_n_geometries
               stokes_S_PT(v,o1,t,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
             enddo
           enddo

!  Perform 2OS correction
!  ----------------------

!  Set conditions

           VLIDORT_FixIn%Cont%TS_nstokes = nstokes_saved

!  Call 2OS correction routine

           CALL VLIDORT_2OSCORR_LPS_MASTER ( &
               VLIDORT_FixIn, &
               VLIDORT_ModIn, &
               VLIDORT_Sup, &
               VLIDORT_Out, &
               VLIDORT_LinFixIn, &
               VLIDORT_LinModIn, &
               VLIDORT_LinSup, &
               VLIDORT_LinOut )

           if ((debug).and.(TID == 0)) then
             !write(*,*)
             write(*,*)'Done 2OS correction for FD task #', t
           end if

!  Exception handling, write-up (optional)

!mick fix 5/1/2015 - created new status subroutine for 2OS code 
           CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
              'V2OS_full_LPS_OMPExecution_TID_' // charTID // '.log', &
              35+TID, OPENFILEFLAG, &
              VLIDORT_Out%Status )

           IF ( VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
                VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

           do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
             do v = 1, VLIDORT_Out%Main%TS_n_geometries
               stokes_C_PT(v,o1,t,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
             enddo
           enddo

!  End task loop

          enddo

!  End FD control if block

        endif

!  Save OMP task information

        VLIDORT_TID_SAVE(wvn) = TID + 1

!  End "wavenumber" loop

      enddo

!$OMP END DO

      !IF (TID == 0) call cpu_time(omp_e2)

!  End parallel region

!$OMP END PARALLEL

!  Timing tests (OpenMP)
!  =====================

      call cpu_time(omp_e2)

      if (omp_nthreads <= n_core) then
        time_divider = omp_nthreads
      else
        time_divider = n_core
      endif

      write(*,*)
      write(*,'(15x,a)')                    'Timing report'
      write(*,'(4(1x,a))') 'Numwvn','# OMP Threads','Time (sec)'
      write(*,'(4(1x,a))') '------','-------------','----------'
      write(*,'(1x,i5,8x,i1,7x,f6.2)') &
        numwvn, omp_nthreads, (omp_e2 - omp_e1)/real(time_divider)

!  Output file control if block - turn off output file writing during either VLID-NS3 or VLID-NS1/V20S time check

      if (.not.do_vlid_ns3_timecheck_only .and. .not.do_vlid_ns1_v2os_timecheck_only) then

!  Prepare output file

        if ( do_std_output_file ) then
          !place all output in one file
          write(ch1,'(i1)')   omp_nthreads
          write(ch2,'(i5.5)') numwvn 
          tail = '_nt' // ch1 // '_nwn' // ch2

          OPEN(36,file='V2OS_full_LPS_results.all' &
                  // '_OMP' // trim(tail),status='replace')
        endif

!  Write output. (Just a sampling of Results)

        !do wvn = 1, numwvn
        do wvn = 1, numwvn, numwvn/omp_nthreads

          if ( .not. do_std_output_file ) then
            !place one set of outputs per file
            write(ch1,'(i1)')   omp_nthreads
            write(ch2,'(i5.5)') numwvn 
            write(ch3,'(i5.5)') wvn
            tail = '_nt' // ch1 // '_nwn' // ch2 // '_wvn' // ch3

            OPEN(36,file='V2OS_full_LPS_results.all' &
                 // '_OMP' // trim(tail),status='replace')
          endif

          v = 1
          write(36,77)'Baseline Stokes I (V,S,C) = ',stokes_V_BAS(v,1,wvn),stokes_S_BAS(v,1,wvn),stokes_C_BAS(v,1,wvn)
          write(36,77)'Baseline Stokes Q (V,S,C) = ',stokes_V_BAS(v,2,wvn),stokes_S_BAS(v,2,wvn),stokes_C_BAS(v,2,wvn)
          write(36,77)'Baseline Stokes U (V,S,C) = ',stokes_V_BAS(v,3,wvn),stokes_S_BAS(v,3,wvn),stokes_C_BAS(v,3,wvn)
          write(36,78)'Baseline Icorrection      = ',Icorr(v,wvn)

!mick fix 5/1/2015 - define q!
          q = 1
          do t = 0, nlay
             if ((t.eq.0) .or. (t.eq.1)) write(36,*)
             if ( t.eq.0 ) then
                 do o1 = 1, 3
                   write(36,56)'Surface WF (V,S,C), Stokes # = ',O1, &
                        spars_save*surfacewf_V_BAS(q,v,o1,wvn),(STOKES_V_PT(v,o1,t,wvn)-stokes_V_BAS(v,o1,wvn))/eps, &
                        spars_save*surfacewf_S_BAS(q,v,o1,wvn),(STOKES_S_PT(v,o1,t,wvn)-stokes_S_BAS(v,o1,wvn))/eps, &
                        spars_save*surfacewf_C_BAS(q,v,o1,wvn),(STOKES_C_PT(v,o1,t,wvn)-stokes_C_BAS(v,o1,wvn))/eps
                 enddo
             else
                 write(36,55)'Layer = ',t
                 do o1 = 1, 3
                   write(36,56)'Profile WF (V,S,C), Stokes # = ',O1, &
                        profilewf_V_BAS(q,t,v,o1,wvn),(STOKES_V_PT(v,o1,t,wvn)-stokes_V_BAS(v,o1,wvn))/eps,   &
                        profilewf_S_BAS(q,t,v,o1,wvn),(STOKES_S_PT(v,o1,t,wvn)-stokes_S_BAS(v,o1,wvn))/eps, &
                        profilewf_C_BAS(q,t,v,o1,wvn),(STOKES_C_PT(v,o1,t,wvn)-stokes_C_BAS(v,o1,wvn))/eps
                 enddo
             endif
          enddo

          if ( .not. do_std_output_file ) close(36)

!  End file write loop

        enddo

        if ( do_std_output_file ) close(36)

!  End output file writing control if block

      end if

!  Continuation point

 567  continue

      IF ( .NOT. PART2 ) STOP

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!                     P A R T    I I

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Baseline Column WF calculation
!  ==============================

      write(*,*)
      write(*,*)'Doing COLUMNWF part: '// &
        'Intensities + 1 Column Jacobian + 1 Surface Jacobian'

!  Start timing

      call cpu_time(omp_e1)

!     Only master thread does this

      !OMP_NTHREADS = OMP_GET_NUM_THREADS()
      !write(*,*)
      !write(*,'(1x,a,i1)') 'Total number of threads (to start) = ', OMP_NTHREADS

!     Initialize output

      STOKES_V_BAS = ZERO
      STOKES_S_BAS = ZERO
      STOKES_C_BAS = ZERO
      Icorr = ZERO

      COLUMNWF_V_BAS  = ZERO
      COLUMNWF_S_BAS  = ZERO
      COLUMNWF_C_BAS  = ZERO

      SURFACEWF_V_BAS = ZERO
      SURFACEWF_S_BAS = ZERO
      SURFACEWF_C_BAS = ZERO

      STOKES_V_PT  = ZERO
      STOKES_S_PT  = ZERO
      STOKES_C_PT  = ZERO

!  Initialize error file output flag

      OPENFILEFLAG = .false.

!  Initialize shared input check flag

      VLIDORT_STATUS_INPUTCHECK = 0

!     Begin parallel region

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (OMP_NTHREADS, &
!$OMP     DO_VLID_NS3_TIMECHECK_ONLY, DO_VLID_NS1_V2OS_TIMECHECK_ONLY, &
!$OMP     NUMWVN, NLAY, &
!$OMP     opd_save, ssa_save, spars_save, eps, &
!$OMP     VLIDORT_FixIn_SAVE,    VLIDORT_ModIn_SAVE,    VLIDORT_Sup_SAVE, &
!$OMP     VLIDORT_LinFixIn_SAVE, VLIDORT_LinModIn_SAVE, VLIDORT_LinSup_SAVE, &
!$OMP     STOKES_V_BAS, STOKES_S_BAS, STOKES_C_BAS, Icorr, &
!$OMP     COLUMNWF_V_BAS,  COLUMNWF_S_BAS,  COLUMNWF_C_BAS,  &
!$OMP     SURFACEWF_V_BAS, SURFACEWF_S_BAS, SURFACEWF_C_BAS, &
!$OMP     STOKES_V_PT, STOKES_S_PT, STOKES_C_PT, &
!$OMP     VLIDORT_STATUS_INPUTCHECK, VLIDORT_TID_SAVE ) &
!$OMP   FIRSTPRIVATE (OPENFILEFLAG)

!     Obtain thread number

      TID = OMP_GET_THREAD_NUM()
      !TID = 0

!     Obtain and display total number of threads and local thread

      IF (TID == 0) THEN
        OMP_NTHREADS = OMP_GET_NUM_THREADS()
        !OMP_NTHREADS = 1
        !write(*,*)
        !write(*,'(1x,a,i1)') 'Total number of threads (inside parallel region) = ', OMP_NTHREADS
        write(*,*)
        write(*,'(1x,a,i1,a)') 'Running V2OS driver with ',OMP_NTHREADS,' thread(s)'
      END IF

      write(*,'(1x,a,i1,a)') 'Thread ',TID,' beginning its wavenumber loop'

!$OMP DO

!  Begin "wavenumber" loop

      do wvn = 1, numwvn
      !do wvn = 1, 1

        write(*,'(1x,a,i1,a,i5)') 'Thread  = ',TID,' wvn = ',wvn

!  Initialize some inputs

        VLIDORT_FixIn = VLIDORT_FixIn_SAVE
        VLIDORT_ModIn = VLIDORT_ModIn_SAVE
        VLIDORT_Sup   = VLIDORT_Sup_SAVE

        VLIDORT_LinFixIn = VLIDORT_LinFixIn_SAVE
        VLIDORT_LinModIn = VLIDORT_LinModIn_SAVE
        VLIDORT_LinSup   = VLIDORT_LinSup_SAVE

!  Re-set the two basic inputs for PART II (VERY IMPORTANT !!!!!!!!!!!!!!!!!!! )
!  Mick note: the two variables below now reset in the triplets of type structure
!             initialization statements above

        !VLIDORT_FixIn%Cont%TS_NSTOKES            = nstokes
        !VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .false.

!  linearization control

        VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = .false.
        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAY)   = .true.
        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAY) = 1
        VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = 0
        VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS         = 1
        VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 1

        VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = .true.
        VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = .false.
        VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .true.
        VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .true.
!mick fix 5/1/2015 - add control variable
        VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = .true.

!  Additional proxies

        nlayers            = nlay
        n_surface_wfs      = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS
        n_totalprofile_wfs = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
        n_totalcolumn_wfs  = VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS

!  Linearization inputs. 1 Column WF, WRT Total Gas absorption

        do n = 1, nlay
          n1 = nlay + 1 - n ; coalb = one - ssa_save(n)
          VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(1,n1)  =  coalb
          VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(1,n1)  = -coalb
        enddo

!  Precautions for VLIDORT

        do n = 1, nlay
           VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0,n,1) = one
           omega = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)
           if ( omega.gt.0.99999999d0 ) VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n) = 0.99999999d0
        enddo

!  First Call. VECTOR VLIDORT
!  ==========================

!  Set conditions

        nstokes_saved = VLIDORT_FixIn%Cont%TS_nstokes 
        VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .false.

!  Begin VLID-NS3 control if block - turn off VLID-NS3 only during VLID-NS1/V20S time check
!  mick note: if block started here as the above two conditions need set for VLID-NS1/V2OS also

        if (.not. do_vlid_ns1_v2os_timecheck_only) then

!  VLIDORT call

          CALL VLIDORT_LCS_MASTER ( &
            VLIDORT_FixIn, &
            VLIDORT_ModIn, &
            VLIDORT_Sup, &
            VLIDORT_Out, &
            VLIDORT_LinFixIn, &
            VLIDORT_LinModIn, &
            VLIDORT_LinSup, &
            VLIDORT_LinOut )

          if ((debug).and.(TID == 0)) then
            !write(*,*)
            write(*,*)'Done Vector VLIDORT, Baseline LCS calculation'
          end if

!  Exception handling, write-up (optional)

          write(charTID,'(i1)') TID
          CALL VLIDORT_WRITE_STATUS ( &
              'V2OS_full_LCS_OMPExecution_TID_' // charTID // '.log', &
              35+TID, OPENFILEFLAG, &
              VLIDORT_Out%Status )

          IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
               VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save baseline results

          do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
             do v = 1, VLIDORT_Out%Main%TS_n_geometries
                stokes_V_BAS(v,o1,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
                do q = 1, n_totalcolumn_wfs
                   columnwf_V_BAS(q,v,o1,wvn) = VLIDORT_LinOut%Col%TS_columnwf(q,1,v,o1,upidx)
                enddo
                do q = 1, n_surface_wfs
                   surfacewf_V_BAS(q,v,o1,wvn) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,upidx)
                enddo
             enddo
          enddo

!  End VLID-NS3 control if block

        endif


!  Begin VLID-NS1/V20S control if block - turn off VLID-NS1/V20S only during VLID-NS3 time check

        if (.not. do_vlid_ns3_timecheck_only) then

!  Second Call. SCALAR VLIDORT + 2OS CORRECTION
!  ============================================

!  Set conditions

          VLIDORT_FixIn%Cont%TS_nstokes            = 1
          VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .true.

!  VLIDORT call

          CALL VLIDORT_LCS_MASTER ( &
            VLIDORT_FixIn, &
            VLIDORT_ModIn, &
            VLIDORT_Sup, &
            VLIDORT_Out, &
            VLIDORT_LinFixIn, &
            VLIDORT_LinModIn, &
            VLIDORT_LinSup, &
            VLIDORT_LinOut )

          if ((debug).and.(TID == 0)) then
            !write(*,*)
            write(*,*)'Done Scalar VLIDORT, Baseline LCS calculation'
          end if

!  Exception handling, write-up (optional)

          CALL VLIDORT_WRITE_STATUS ( &
            'V2OS_full_LCS_OMPExecution_TID_' // charTID // '.log', &
            35+TID, OPENFILEFLAG, &
            VLIDORT_Out%Status )

          IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
               VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

          do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
             do v = 1, VLIDORT_Out%Main%TS_n_geometries
                stokes_S_BAS(v,o1,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
                do q = 1, n_totalcolumn_wfs
                   columnwf_S_BAS(q,v,o1,wvn) = VLIDORT_LinOut%Col%TS_columnwf(q,1,v,o1,upidx)
                enddo
                do q = 1, n_surface_wfs
                  surfacewf_S_BAS(q,v,o1,wvn) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,upidx)
                enddo
             enddo
          enddo

!  Perform 2OS correction
!  ----------------------

!  Set conditions

          VLIDORT_FixIn%Cont%TS_nstokes = nstokes_saved

!  Call 2OS LCS correction routine

          CALL VLIDORT_2OSCORR_LCS_MASTER ( &
            VLIDORT_FixIn, &
            VLIDORT_ModIn, &
            VLIDORT_Sup, &
            VLIDORT_Out, &
            VLIDORT_LinFixIn, &
            VLIDORT_LinModIn, &
            VLIDORT_LinSup, &
            VLIDORT_LinOut )

          if ((debug).and.(TID == 0)) then
            !write(*,*)
            write(*,*)'Done 2OS correction, Baseline LCS calculation'
          end if

!  Exception handling, write-up (optional)

!mick fix 5/1/2015 - created new status subroutine for 2OS code 
          CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
            'V2OS_full_LCS_OMPExecution_TID_' // charTID // '.log', &
            35+TID, OPENFILEFLAG, &
            VLIDORT_Out%Status )

          IF ( VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
               VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

          do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
             do v = 1, VLIDORT_Out%Main%TS_n_geometries
                stokes_C_BAS(v,o1,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
                do q = 1, n_totalcolumn_wfs
                   columnwf_C_BAS(q,v,o1,wvn) = VLIDORT_LinOut%Col%TS_columnwf(q,1,v,o1,upidx)
                enddo
                do q = 1, n_surface_wfs
                  surfacewf_C_BAS(q,v,o1,wvn) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,upidx)
                enddo
             enddo
          enddo

!          stop'baseline Column WFS only'

!  End VLID-NS1/V20S control if block

        endif


!  Begin FD control if block - turn off FD during either VLID-NS3 or VLID-NS1/V20S time check

        if (.not.do_vlid_ns3_timecheck_only .and. .not.do_vlid_ns1_v2os_timecheck_only) then

!  Initialize FD section
!  =====================

!  Fd perturbation

          eps = 1.0d-03
          epsfac = 1.0d0 + eps

!  Initialize linearized inputs

          VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = .true.
          VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAY)   = .false.
          VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAY) = 0
          VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS         = 0
          VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 0

          VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = .false.
          VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .false.
          VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .false.
!mick fix 5/1/2015 - add control variable
          VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = .false.

          VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT = zero
          VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT = zero

!  Additional proxies

          n_surface_wfs      = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS
          n_totalcolumn_wfs  = VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS

!  TASKS for FD testing
!  ======================

!  Task 0, FD Perturbation on Surface Albedo
!  Task 1, FD Perturbation on trace gas total optical depth

          do task = 0, 1
!          do task = 1, 1

!  Task counter

             t = task

!  Task 0 - surface jacobian

             if ( task .eq. 0 ) then
                if ((debug).and.(TID == 0)) &
                  write(*,*)'Doing task 0: Finite Difference, perturb Lambertian albedo'
                VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =   spars_save * epsfac
                do n = 1, nlay
                   n1 = nlay + 1 - n
                   VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
                   VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
                enddo
             endif

!  Task 1 - column jacobian

             if ( task.gt.0  ) then
                if ((debug).and.(TID == 0)) &
                  write(*,'(a,i2,a)')' Doing task',task, &
                    ': Finite Difference, perturb total molecular absorption'
                VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =   spars_save
                do n = 1, nlay
                   n1 = nlay + 1 - n
                   scat = opd_save(n) * ssa_save(n)
                   gaspert = ( opd_save(n) - scat ) * epsfac
                   taupert = scat + gaspert
                   VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = taupert
                   VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = scat / taupert
                enddo
             endif

!  Precautions for VLIDORT

             do n = 1, nlay
                omega = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)
                if ( omega.gt.0.99999999d0 ) VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n) = 0.99999999d0
             enddo

!  First Call. VECTOR VLIDORT
!  ==========================

!  Set conditions

             nstokes_saved = VLIDORT_FixIn%Cont%TS_nstokes 
             VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .false.

!  VLIDORT call

             CALL VLIDORT_LCS_MASTER ( &
               VLIDORT_FixIn, &
               VLIDORT_ModIn, &
               VLIDORT_Sup, &
               VLIDORT_Out, &
               VLIDORT_LinFixIn, &
               VLIDORT_LinModIn, &
               VLIDORT_LinSup, &
               VLIDORT_LinOut )

             if ((debug).and.(TID == 0)) then
               !write(*,*)
               write(*,*)'Done Vector VLIDORT for FD task #', t
             end if

!  Exception handling, write-up (optional)

             CALL VLIDORT_WRITE_STATUS ( &
               'V2OS_full_LCS_OMPExecution_TID_' // charTID // '.log', &
               35+TID, OPENFILEFLAG, &
               VLIDORT_Out%Status )

             IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
                  VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

             do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
               do v = 1, VLIDORT_Out%Main%TS_n_geometries
                 stokes_V_PT(v,o1,t,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
               enddo
             enddo

!  Second Call. SCALAR VLIDORT + 2OS CORRECTION
!  ============================================

!  Set conditions

             VLIDORT_FixIn%Cont%TS_nstokes            = 1
             VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .true.

!  VLIDORT call

             CALL VLIDORT_LCS_MASTER ( &
               VLIDORT_FixIn, &
               VLIDORT_ModIn, &
               VLIDORT_Sup, &
               VLIDORT_Out, &
               VLIDORT_LinFixIn, &
               VLIDORT_LinModIn, &
               VLIDORT_LinSup, &
               VLIDORT_LinOut )

             if ((debug).and.(TID == 0)) then
               !write(*,*)
               write(*,*)'Done Scalar VLIDORT for FD task #', t
             end if

!  Exception handling, write-up (optional)

             CALL VLIDORT_WRITE_STATUS ( &
               'V2OS_full_LCS_OMPExecution_TID_' // charTID // '.log', &
               35+TID, OPENFILEFLAG, &
               VLIDORT_Out%Status )

             IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
                  VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

             do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
               do v = 1, VLIDORT_Out%Main%TS_n_geometries
                 stokes_S_PT(v,o1,t,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
               enddo
             enddo

!  Perform 2OS correction
!  ----------------------

!  Set conditions

             VLIDORT_FixIn%Cont%TS_nstokes = nstokes_saved

!  Call 2OS correction routine

             CALL VLIDORT_2OSCORR_LCS_MASTER ( &
                 VLIDORT_FixIn, &
                 VLIDORT_ModIn, &
                 VLIDORT_Sup, &
                 VLIDORT_Out, &
                 VLIDORT_LinFixIn, &
                 VLIDORT_LinModIn, &
                 VLIDORT_LinSup, &
                 VLIDORT_LinOut )

             if ((debug).and.(TID == 0)) then
               !write(*,*)
               write(*,*)'Done 2OS correction for FD task #', t
             end if

!  Exception handling, write-up (optional)

!mick fix 5/1/2015 - created new status subroutine for 2OS code 
             CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
               'V2OS_full_LCS_OMPExecution_TID_' // charTID // '.log', &
               35+TID, OPENFILEFLAG, &
               VLIDORT_Out%Status )

             IF ( VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
                  VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

             do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
               do v = 1, VLIDORT_Out%Main%TS_n_geometries
                 stokes_C_PT(v,o1,t,wvn) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
               enddo
             enddo

!  End task loop

          enddo

!  End FD control if block

        endif

!  Save OMP thread information

        VLIDORT_TID_SAVE(wvn) = TID + 1

!  End "wavenumber" loop

      enddo

!$OMP END DO

      !IF (TID == 0) call cpu_time(omp_e2)

!  End parallel region

!$OMP END PARALLEL

!  Timing tests (OpenMP)
!  =====================

      call cpu_time(omp_e2)

      if (omp_nthreads <= n_core) then
        time_divider = omp_nthreads
      else
        time_divider = n_core
      endif

      write(*,*)
      write(*,'(15x,a)')                    'Timing report'
      write(*,'(4(1x,a))') 'Numwvn','# OMP Threads','Time (sec)'
      write(*,'(4(1x,a))') '------','-------------','----------'
      write(*,'(1x,i5,8x,i1,7x,f6.2)') &
        numwvn, omp_nthreads, (omp_e2 - omp_e1)/real(time_divider)

!  Prepare output file

      if ( do_std_output_file ) then
        !place all output in one file
        write(ch1,'(i1)')   omp_nthreads
        write(ch2,'(i5.5)') numwvn 
        tail = '_nt' // ch1 // '_nwn' // ch2

        OPEN(36,file='V2OS_full_LCS_results.all' &
               // '_OMP' // trim(tail),status='replace')
      endif

!  Output file control if block - turn off output file writing during either VLID-NS3 or VLID-NS1/V20S time check

      if (.not.do_vlid_ns3_timecheck_only .and. .not.do_vlid_ns1_v2os_timecheck_only) then

!  Write output. (Just a sampling of Results)

        !do wvn = 1, numwvn
        do wvn = 1, numwvn, numwvn/omp_nthreads

          if ( .not. do_std_output_file ) then
            !place one set of outputs per file
            write(ch1,'(i1)')   omp_nthreads
            write(ch2,'(i5.5)') numwvn 
            write(ch3,'(i5.5)') wvn
            tail = '_nt' // ch1 // '_nwn' // ch2 // '_wvn' // ch3

            OPEN(36,file='V2OS_full_LCS_results.all' &
                 // '_OMP' // trim(tail),status='replace')
          endif

          v = 1
          write(36,77)'Baseline Stokes I (V,S,C) = ',stokes_V_BAS(v,1,wvn),stokes_S_BAS(v,1,wvn),stokes_C_BAS(v,1,wvn)
          write(36,77)'Baseline Stokes Q (V,S,C) = ',stokes_V_BAS(v,2,wvn),stokes_S_BAS(v,2,wvn),stokes_C_BAS(v,2,wvn)
          write(36,77)'Baseline Stokes U (V,S,C) = ',stokes_V_BAS(v,3,wvn),stokes_S_BAS(v,3,wvn),stokes_C_BAS(v,3,wvn)
          write(36,78)'Baseline Icorrection      = ',Icorr(v,wvn)

!mick fix 5/1/2015 - changed hardcoded jac value of 1 in arrays below to use q now; defined q here
          q = 1
          do t = 0, 1
             write(36,*)
             if ( t.eq.0 ) then
                 do o1 = 1, 3
                   write(36,56)'Surface WF (V,S,C), Stokes # = ',O1, &
                        spars_save*surfacewf_V_BAS(q,v,o1,wvn),(STOKES_V_PT(v,o1,t,wvn)-stokes_V_BAS(v,o1,wvn))/eps, &
                        spars_save*surfacewf_S_BAS(q,v,o1,wvn),(STOKES_S_PT(v,o1,t,wvn)-stokes_S_BAS(v,o1,wvn))/eps, &
                        spars_save*surfacewf_C_BAS(q,v,o1,wvn),(STOKES_C_PT(v,o1,t,wvn)-stokes_C_BAS(v,o1,wvn))/eps
                 enddo
             else
                 do o1 = 1, 3
                   write(36,56)'Column  WF (V,S,C), Stokes # = ',O1, &
                        columnwf_V_BAS(q,v,o1,wvn),(STOKES_V_PT(v,o1,t,wvn)-stokes_V_BAS(v,o1,wvn))/eps, &
                        columnwf_S_BAS(q,v,o1,wvn),(STOKES_S_PT(v,o1,t,wvn)-stokes_S_BAS(v,o1,wvn))/eps, &
                        columnwf_C_BAS(q,v,o1,wvn),(STOKES_C_PT(v,o1,t,wvn)-stokes_C_BAS(v,o1,wvn))/eps
                 enddo
             endif
          enddo

          if ( .not. do_std_output_file ) close(36)

!  End file write loop

        enddo

        if ( do_std_output_file ) close(36)

!  End output file writing control if block

      endif

!  End tests loop

      enddo

!  Format statements

55    format(A,i2)
56    format(A,i3, 4x, 3(1p2e17.8,2x))
77    format(A,1p3e17.8)
78    format(A,1pe17.8)

!  Finish

      write(*,*)
      write(*,*) '******************************************'
      write(*,*)

      if ( VLIDORT_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) THEN
         stop 'Successful run, with Warning: look in ' // &
              'V2OS_full_LPS_OMPExecution_TID_*.log or V2OS_full_LCS_OMPExecution_TID_*.log'
      else
         stop 'Successful run'
      endif

      end program V2OS_full_LPCS_OMPTester
