      program TWOOS_Self_LPCS_OMPTester

!  This is a Self tester for the linearized 2OS code, validation against Vijay's latest.
!  3 March 2015 (OMP testing)

!  Module files for VLIDORT

      USE VLIDORT_PARS
      USE VLIDORT_IO_DEFS
      USE VLIDORT_LIN_IO_DEFS

      USE VLIDORT_AUX
      USE vlidort_2OScorr_master_m
      USE vlidort_2OScorr_lps_master_m
      USE vlidort_2OScorr_lcs_master_m
!mick fix 4/21/2015 - added new 2OS module
      USE vlidort_2OS_corr_aux

      IMPLICIT NONE

!  Dimensioning integer parameters

      INTEGER, PARAMETER :: MAXWVN = 10!4!50!100!1000

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
                          Sun_Chapman2(maxlayers,maxlayers)
      double precision :: ssa(maxlayers),opd(maxlayers),coefs(0:maxmoments_input,maxlayers,6)


!  Saved results

      !INTEGER, DIMENSION(:), ALLOCATABLE :: OMP_N_GEOMETRIES
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: R2_PTSave
      DOUBLE PRECISION, DIMENSION(:,:)  , ALLOCATABLE :: Icorr_PTSave

      DOUBLE PRECISION, DIMENSION(:,:)  , ALLOCATABLE :: R2_Baseline
      DOUBLE PRECISION, DIMENSION(:)    , ALLOCATABLE :: Icorr_Baseline
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: R2_LPBaseline
      DOUBLE PRECISION, DIMENSION(:,:)  , ALLOCATABLE :: Icorr_LPBaseline
      DOUBLE PRECISION, DIMENSION(:,:)  , ALLOCATABLE :: R2_LCBaseline
      DOUBLE PRECISION, DIMENSION(:)    , ALLOCATABLE :: Icorr_LCBaseline
      DOUBLE PRECISION, DIMENSION(:,:)  , ALLOCATABLE :: R2_LSBaseline
      DOUBLE PRECISION, DIMENSION(:)    , ALLOCATABLE :: Icorr_LSBaseline

!  Local shared array

!mick fix 4/21/2015 - added VLIDORT_STATUS_INPUTCHECK
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

      INTEGER       :: n_core, time_divider
      REAL          :: omp_e1, omp_e2

!  OpenMP functions

      INTEGER :: OMP_GET_NUM_THREADS, &
                 OMP_GET_THREAD_NUM

!  Flag for opening error output file

      LOGICAL :: OPENFILEFLAG, PART1, PART2

!  flag for standard output file

      !LOGICAL, parameter :: do_std_output_file = .true.
      LOGICAL, parameter :: do_std_output_file = .false.

!  local

      integer :: n, n1, k, l, nd, kd, ld, t, task, o1
      double precision :: eps, epsfac, gaspert, scat, taupert, coalb
      double precision :: ssa_save(maxlayers),opd_save(maxlayers),spars_save

!  Start of code
!  =============

!  Define some control variables

!      PART1 = .FALSE.
      PART1 = .TRUE.

!      PART2 = .FALSE.
      PART2 = .TRUE.

!  Set test parameters

      N_CORE = 2!4
      N_OMP_TESTS = 2!3

      numwvn = maxwvn
!      ntasks = maxtasks

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

!  Start timing

      call cpu_time(omp_e1)

!  Fixed Boolean inputs
!  --------------------

      VLIDORT_FixIn_Save%Bool%TS_DO_FULLRAD_MODE      = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_SSCORR_TRUNCATION = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_SS_EXTERNAL       = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_SSFULL            = .false.
      VLIDORT_FixIn_Save%Bool%TS_DO_THERMAL_EMISSION  = .true.
      VLIDORT_FixIn_Save%Bool%TS_DO_SURFACE_EMISSION  = .true.

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
      VLIDORT_FixIn_Save%Cont%TS_NFINELAYERS      = 0
      VLIDORT_FixIn_Save%Cont%TS_N_THERMAL_COEFFS = 0

      VLIDORT_FixIn_Save%Cont%TS_TAYLOR_ORDER     = 0

      VLIDORT_FixIn_Save%Cont%TS_VLIDORT_ACCURACY  = zero ! Set later

      VLIDORT_FixIn_Save%Cont%TS_NLAYERS_NOMS   = 0
      VLIDORT_FixIn_Save%Cont%TS_NLAYERS_CUTOFF = 0

!  Fixed Beam Userval inputs
!  --------------------------

      VLIDORT_FixIn_Save%Sunrays%TS_FLUX_FACTOR     = ONE
      VLIDORT_FixIn_Save%UserVal%TS_N_USER_LEVELS   = 1

!  Chapman inputs
!  --------------

!mick fix 4/21/2015 - added HEIGHT_GRID initialization
      VLIDORT_FixIn_Save%Chapman%TS_HEIGHT_GRID       = ZERO
      VLIDORT_FixIn_Save%Chapman%TS_PRESSURE_GRID     = ZERO
      VLIDORT_FixIn_Save%Chapman%TS_TEMPERATURE_GRID  = ZERO
      VLIDORT_FixIn_Save%Chapman%TS_FINEGRID          = 0
      VLIDORT_FixIn_Save%Chapman%TS_RFINDEX_PARAMETER = Zero

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
!mick fix 4/21/2015 - added DO_FO_CALC initialization
      VLIDORT_ModIn_Save%MBool%TS_DO_FO_CALC             = .false.

      VLIDORT_ModIn_Save%MBool%TS_DO_REFRACTIVE_GEOMETRY = .false.
      VLIDORT_ModIn_Save%MBool%TS_DO_CHAPMAN_FUNCTION    = .false.   ! Not required

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

      VLIDORT_ModIn_Save%MBool%TS_DO_2OS_CORRECTION       = .true.  ! Must be set

!  Modified Control
!  ----------------

      VLIDORT_ModIn_Save%MCont%TS_NGREEK_MOMENTS_INPUT = 0 ! Set later

!  Modified Sunrays
!  ----------------

      VLIDORT_ModIn_Save%MSunrays%TS_N_SZANGLES        = 1
      VLIDORT_ModIn_Save%MSunrays%TS_SZANGLES          = zero ! set later

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

      VLIDORT_ModIn_Save%MChapman%TS_EARTH_RADIUS        = zero
      VLIDORT_ModIn_Save%MChapman%TS_CHAPMAN_FACTORS     = zero ! Set later

!  Modified Optical
!  ----------------

      VLIDORT_ModIn_Save%MOptical%TS_OMEGA_TOTAL_INPUT = zero ! Set later

!  Linearized
!  ----------

      VLIDORT_LinFixIn_Save%Cont%TS_LAYER_VARY_FLAG    = .false. ! Set later
      VLIDORT_LinFixIn_Save%Cont%TS_LAYER_VARY_NUMBER  = 0       ! Set later
!mick fix 4/21/2015 - added these initializations
      VLIDORT_LinFixIn_Save%Cont%TS_N_TOTALCOLUMN_WFS  = 0       ! Set later
      VLIDORT_LinFixIn_Save%Cont%TS_N_TOTALPROFILE_WFS = 0       ! Set later
      VLIDORT_LinFixIn_Save%Cont%TS_N_SURFACE_WFS      = 0       ! Set later
      VLIDORT_LinFixIn_Save%Cont%TS_N_SLEAVE_WFS       = 0       ! Set later
      VLIDORT_LinFixIn_Save%Cont%TS_COLUMNWF_NAMES     = ' '
      VLIDORT_LinFixIn_Save%Cont%TS_PROFILEWF_NAMES    = ' '

      VLIDORT_LinFixIn_Save%Optical%TS_L_deltau_vert_input    = zero ! Set later
      VLIDORT_LinFixIn_Save%Optical%TS_L_greekmat_total_input = zero
      VLIDORT_LinFixIn_Save%Optical%TS_L_omega_total_input    = zero ! Set later

!mick fix 4/21/2015 - added these initializations
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

!mick fix 4/21/2015 - added SUP initializations
!  SS and DB stuff may be set later

      VLIDORT_Sup_SAVE%SS%TS_STOKES_SS  = ZERO
      VLIDORT_Sup_SAVE%SS%TS_STOKES_DB  = ZERO
      VLIDORT_Sup_SAVE%SS%TS_FO_ZMATRIX = ZERO
      VLIDORT_Sup_SAVE%SS%TS_FO_R1SAVED = ZERO

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

!mick fix 4/21/2015 - added linearized SS & DB initializations
      VLIDORT_LinSup_Save%SS%Prof%TS_PROFILEWF_SS = ZERO
      VLIDORT_LinSup_Save%SS%Prof%TS_PROFILEWF_DB = ZERO
      VLIDORT_LinSup_Save%SS%Col%TS_COLUMNWF_SS   = ZERO
      VLIDORT_LinSup_Save%SS%Col%TS_COLUMNWF_DB   = ZERO
      VLIDORT_LinSup_Save%SS%Surf%TS_SURFACEWF_DB = ZERO

!  Now open the files

      open(87,file='fort.87_Rg',status='old')
      read(87,*)nlay,nmug,nstokes,nmoms,nfoumax,regular_ps,enhanced_ps
      read(87,*)epsilon,phi,emu0,emu,spars ; spars_save = spars
      do k = 1, nmug
         read(87,*)kd,xmu(k),w(k)
      enddo
      do n = 1, nlay
         read(87,*)ssa(n),opd(n),ncoefs(n)
         opd_save(n) = opd(n) ; ssa_save(n) = ssa(n)
         do l = 0, 2*nmug-1
            read(87,*)nd,ld,(coefs(l,n,k),k=1,6)
         enddo
      enddo
      close(87)

!mick fix 4/21/2015 - added Sun_Chapman2 initialization
      Sun_Chapman2 = 0.0d0
      open(88,file='fort.88',status='old')
      do n = 1, nlay
         do k = 1, n
         read(88,*)nd,kd,Sun_Chapman2(n,k)
         enddo
      enddo
      close(88)

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

      VLIDORT_ModIn_Save%MSunrays%TS_SZANGLES             = acos(emu0)/deg_to_rad
      VLIDORT_ModIn_Save%MUserVal%TS_USER_VZANGLES_INPUT  = acos(emu) /deg_to_rad
      VLIDORT_ModIn_Save%MUserVal%TS_USER_RELAZMS         = phi

      VLIDORT_ModIn_Save%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1) = acos(emu0)/deg_to_rad
      VLIDORT_ModIn_Save%MUserVal%TS_USER_OBSGEOMS_INPUT(1,2) = acos(emu) /deg_to_rad
      VLIDORT_ModIn_Save%MUserVal%TS_USER_OBSGEOMS_INPUT(1,3) = phi

!mick fix 4/20/2015 - trimmed argument passing
      !VLIDORT_ModIn_Save%MChapman%TS_CHAPMAN_FACTORS(:,:,1) = Sun_Chapman2(:,:)
      VLIDORT_ModIn_Save%MChapman%TS_CHAPMAN_FACTORS(1:nlay,1:nlay,1) = &
        Sun_Chapman2(1:nlay,1:nlay)

!  Optical

      do n = 1, nlay
        n1 = nlay + 1 - n
        VLIDORT_FixIn_Save%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd(n)
        VLIDORT_ModIn_Save%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa(n)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,1)  =   coefs(0:nmoms,n,1)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,2)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,5)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,6)  =   coefs(0:nmoms,n,2)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,11) =   coefs(0:nmoms,n,3)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,12) =   coefs(0:nmoms,n,6)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,15) = - coefs(0:nmoms,n,6)
        VLIDORT_FixIn_Save%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,16) =   coefs(0:nmoms,n,4)
      enddo
      VLIDORT_FixIn_Save%Optical%TS_LAMBERTIAN_ALBEDO =   spars

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

!     Only master thread does this

      !OMP_NTHREADS = OMP_GET_NUM_THREADS()
      !write(*,*)
      !write(*,'(1x,a,i1)') 'Total number of threads (to start) = ', OMP_NTHREADS

!     Allocate output memory

      ALLOCATE ( &
      !ALLOCATE ( OMP_N_GEOMETRIES ( MAXWVN ), &
        R2_PTSave    ( 0:MAXLAYERS, MAXSTOKES, MAXWVN ), Icorr_PTSave( 0:MAXLAYERS, MAXWVN )   , &
        R2_Baseline  ( MAXSTOKES, MAXWVN )             , Icorr_Baseline   ( MAXWVN )           , &
        R2_LSBaseline( MAXSTOKES, MAXWVN )             , Icorr_LSBaseline ( MAXWVN )           , &
        R2_LPBaseline( MAXLAYERS, MAXSTOKES, MAXWVN )  , Icorr_LPBaseline ( MAXLAYERS, MAXWVN )  )

!     Initialize output

      !OMP_N_GEOMETRIES = 0

      R2_PTSave        = ZERO
      Icorr_PTSave     = ZERO

      R2_Baseline      = ZERO
      Icorr_Baseline   = ZERO

      R2_LSBaseline    = ZERO
      Icorr_LSBaseline = ZERO

      R2_LPBaseline    = ZERO
      Icorr_LPBaseline = ZERO

!  Initialize error file output flag

      OPENFILEFLAG = .false.

!  Initialize shared input check variable

      VLIDORT_STATUS_INPUTCHECK = 0

!     Begin parallel region

!mick fix 4/21/2015 -
!Added:   ssa, opd_save, ssa_save, spars_save, eps,
!         VLIDORT_STATUS_INPUTCHECK
!Removed: OMP_N_GEOMETRIES

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (OMP_NTHREADS, &
!$OMP     NUMWVN, NLAY, &
!$OMP     ssa, opd_save, ssa_save, spars_save, eps, &
!$OMP     VLIDORT_FixIn_SAVE,    VLIDORT_ModIn_SAVE,    VLIDORT_Sup_SAVE, &
!$OMP     VLIDORT_LinFixIn_SAVE, VLIDORT_LinModIn_SAVE, VLIDORT_LinSup_SAVE, &
!$OMP     R2_PTSave,     Icorr_PTSave, &
!$OMP     R2_Baseline,   Icorr_Baseline, &
!$OMP     R2_LPBaseline, Icorr_LPBaseline, &
!$OMP     R2_LSBaseline, Icorr_LSBaseline, &
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
        write(*,'(1x,a,i1,a)') 'Running TwoOS driver with ',OMP_NTHREADS,' thread(s)'
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

!mick fix 4/21/2015 - added initialization
!  Must zero the Main output, since now intent(INOUT)

        VLIDORT_Out%Main%TS_STOKES       = zero
        VLIDORT_LinOut%Prof%TS_PROFILEWF = zero
        VLIDORT_LinOut%Surf%TS_SURFACEWF = zero

!  linearization control

        VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = .false.
        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAY)   = .true.
        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAY) = 1
        VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = 1
        VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS         = 0
        VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 1

        VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .true.
        VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = .false.
        VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = .true.
        VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .true.

!  Linearization inputs. 1 Profile WF, WRT Gas absorption

        do n = 1, nlay
          n1 = nlay + 1 - n ; coalb = one - ssa(n)
          VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(1,n1)  =  coalb
          VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(1,n1)  = -coalb
        enddo

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

        !write(*,*)
        !write(*,*)'Done 2OS LPS Baseline correction'

!  Exception handling, write-up (optional)

        write(charTID,'(i1)') TID
!mick fix 4/21/2015 - created new status subroutine for 2OS code 
        CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
          'TwoOS_self_LPS_OMPExecution_TID_' // charTID // '.log', &
          35+TID, OPENFILEFLAG, &
          VLIDORT_Out%Status )

        IF ( VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
          VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

        !!OMP_N_GEOMETRIES(wvn) = VLIDORT_Out%Main%TS_n_geometries
        !R2_Baseline(1:3,wvn)  = VLIDORT_Out%Main%TS_2OSCORR_R2(1,1:3)
        !Icorr_Baseline (wvn)  = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)
        !do n = 1, nlay
        !   R2_LPBaseline(n,1:3,wvn) = VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(1,n,1,1:3)
        !   Icorr_LPBaseline (n,wvn) = VLIDORT_LinOut%Prof%TS_2OSCORR_LP_ICORR(1,n,1)
        !enddo
        !R2_LSBaseline(1:3,wvn) = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,1:3)
        !Icorr_LSBaseline (wvn) = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_ICORR(1,1)

        do o1=1,3
          R2_Baseline(o1,wvn) = VLIDORT_Out%Main%TS_2OSCORR_R2(1,o1)
          do n = 1, nlay
            R2_LPBaseline(n,o1,wvn) = VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(1,n,1,o1)
          enddo
          R2_LSBaseline(o1,wvn) = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,o1)
        enddo

        Icorr_Baseline(wvn) = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)
        do n = 1, nlay
          Icorr_LPBaseline(n,wvn) = VLIDORT_LinOut%Prof%TS_2OSCORR_LP_ICORR(1,n,1)
        enddo
        Icorr_LSBaseline(wvn) = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_ICORR(1,1)

!  Write results

!      write(*,*)'LPS,R2(1) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,1),VLIDORT_Out%Main%TS_stokes(1,1,1,upidx)
!      write(*,*)'LPS,R2(2) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,2),VLIDORT_Out%Main%TS_stokes(1,1,2,upidx)
!      write(*,*)'LPS,R2(3) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,3),VLIDORT_Out%Main%TS_stokes(1,1,3,upidx)
!      write(*,*)'LPS,Icorr = ',VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

!      write(*,*)'LPS,LS_R2(1) = ',VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,1)
!      write(*,*)'LPS,LS_R2(2) = ',VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,2)
!      write(*,*)'LPS,LS_R2(3) = ',VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,3)
!      write(*,*)'LPS,LS_Icorr = ',VLIDORT_LinOut%Surf%TS_2OSCORR_LS_Icorr(1,1)

!      ks = 36
!      write(*,*)'LPS 36,LP_R2(1) = ',VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(1,ks,1,1)
!      write(*,*)'LPS 36,LP_R2(2) = ',VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(1,ks,1,2)
!      write(*,*)'LPS 36,LP_R2(3) = ',VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(1,ks,1,3)
!      write(*,*)'LPS 36,LP_Icorr = ',VLIDORT_LinOut%Prof%TS_2OSCORR_LP_Icorr(1,ks,1)

!       stop'baseline profile WFS only'

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

        VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .false.
        VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = .false.
        VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .false.

        VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT    = zero
        VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT    = zero

!  TASKS for FD testing
!  ======================

!  Task 0, FD Perturbation on Surface Albedo
!  Task t, FD Perturbation on trace gas optical depth in layer t  (2OS = lowest layer)

        do task = 0, nlay
!        do task = 1, 1

!  Task counter

         t = task

!  Task 0

         if ( task .eq. 0 ) then
            !write(*,*)'Doing task 0: '//'Finite Difference, perturb Lambertian albedo'
            VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = spars_save * epsfac
            do n = 1, nlay
               n1 = nlay + 1 - n
               VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
               VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
            enddo
         endif

!  Tasks 1 to NLAYERS

         if ( task.gt.0  ) then
            !write(*,'(a,i2,a,i2)')' Doing task ',task, &
            ! ': Finite Difference, perturb molecular absorption Layer ',t
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

!  Must zero the Main output, since now intent(INOUT)

         VLIDORT_Out%Main%TS_STOKES = zero

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

!  Exception handling, write-up (optional)

!mick fix 4/21/2015 - created new status subroutine for 2OS code 
         CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
           'TwoOS_self_LPS_OMPExecution_TID_' // charTID // '.log', &
           35+TID, OPENFILEFLAG, &
           VLIDORT_Out%Status )

         IF ( VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
           VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

         !!OMP_N_GEOMETRIES(wvn) = VLIDORT_Out%Main%TS_n_geometries
         !R2_PTSave(t,1:3,wvn)  = VLIDORT_Out%Main%TS_2OSCORR_R2(1,1:3)
         !Icorr_PTSave(t,wvn)   = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

         do o1=1,3
           R2_PTSave(t,o1,wvn) = VLIDORT_Out%Main%TS_2OSCORR_R2(1,o1)
         enddo
         Icorr_PTSave(t,wvn) = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

!  End task loop

      enddo

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

        OPEN(36,file='TwoOS_self_LPS_results.all' &
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

          OPEN(36,file='TwoOS_self_LPS_results.all' &
               // '_OMP' // trim(tail),status='replace')
        endif

        write(36,55)'Baseline R2/Icorr',R2_Baseline(1:3,wvn), Icorr_Baseline(wvn)
        do t = 0, nlay
          if ( t.eq.0 ) then
            write(36,56)'Surface WF R2/Icorr          ', &
                     spars_save*R2_LSBaseline(1:3,wvn), spars_save*Icorr_LSBaseline(wvn), &
                    (R2_PTSave(t,1,wvn)-R2_Baseline(1,wvn))/eps, &
                    (R2_PTSave(t,2,wvn)-R2_Baseline(2,wvn))/eps, &
                    (R2_PTSave(t,3,wvn)-R2_Baseline(3,wvn))/eps, &
                    (Icorr_PTSave(t,wvn)-Icorr_Baseline(wvn))/eps
          else
            write(36,57)'Profile WF R2/Icorr, Layer ',t,& 
                     R2_LPBaseline(t,1:3,wvn), Icorr_LPBaseline(t,wvn), &
                    (R2_PTSave(t,1,wvn)-R2_Baseline(1,wvn))/eps, &
                    (R2_PTSave(t,2,wvn)-R2_Baseline(2,wvn))/eps, &
                    (R2_PTSave(t,3,wvn)-R2_Baseline(3,wvn))/eps, &
                    (Icorr_PTSave(t,wvn)-Icorr_Baseline(wvn))/eps
          endif
        enddo

        if ( .not. do_std_output_file ) close(36)

!  End file write loop

      enddo

      if ( do_std_output_file ) close(36)

!     Deallocate output memory

      !DEALLOCATE ( OMP_N_GEOMETRIES, R2_PTSave, Icorr_PTSave, R2_Baseline, Icorr_Baseline, &
      !             R2_LPBaseline, Icorr_LPBaseline, R2_LSBaseline, Icorr_LSBaseline )
      DEALLOCATE ( R2_PTSave, Icorr_PTSave, R2_Baseline, Icorr_Baseline, &
                   R2_LPBaseline, Icorr_LPBaseline, R2_LSBaseline, Icorr_LSBaseline )

!  Format

55    format(A,    1p3e15.6,2x,1pe15.6)
56    format(A,    1p3e15.6,2x,1pe15.6, 4x, 1p3e15.6,2x,1pe15.6 )
57    format(A,I2, 1p3e15.6,2x,1pe15.6, 4x, 1p3e15.6,2x,1pe15.6 )

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

!     Only master thread does this

      !OMP_NTHREADS = OMP_GET_NUM_THREADS()
      !write(*,*)
      !write(*,'(1x,a,i1)') 'Total number of threads (to start) = ', OMP_NTHREADS

!     Allocate output memory

      ALLOCATE ( & 
        !OMP_N_GEOMETRIES ( MAXWVN ), &
        R2_PTSave    ( 0:1, MAXSTOKES, MAXWVN ), Icorr_PTSave ( 0:1, MAXWVN )      , &
        R2_Baseline  ( MAXSTOKES, MAXWVN )       , Icorr_Baseline   ( MAXWVN )     , &
        R2_LSBaseline( MAXSTOKES, MAXWVN )       , Icorr_LSBaseline ( MAXWVN )     , &
        R2_LCBaseline( MAXSTOKES, MAXWVN )       , Icorr_LCBaseline ( MAXWVN )  )

!     Initialize output

      !OMP_N_GEOMETRIES = 0

      R2_PTSave            = ZERO
      Icorr_PTSave         = ZERO

      R2_Baseline          = ZERO
      Icorr_Baseline       = ZERO

      R2_LSBaseline        = ZERO
      Icorr_LSBaseline     = ZERO

      R2_LCBaseline        = ZERO
      Icorr_LCBaseline     = ZERO

!  Initialize error file output flag

      OPENFILEFLAG = .false.

!  Initialize shared input check variable

      VLIDORT_STATUS_INPUTCHECK = 0

!     Begin parallel region

!mick fix 4/21/2015 -
!Added:   nlay, ssa, opd_save, ssa_save, spars_save, eps,
!         VLIDORT_STATUS_INPUTCHECK
!Removed: OMP_N_GEOMETRIES

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (OMP_NTHREADS, &
!$OMP     NUMWVN, NLAY, &
!$OMP     ssa, opd_save, ssa_save, spars_save, eps, &
!$OMP     VLIDORT_FixIn_SAVE,    VLIDORT_ModIn_SAVE,    VLIDORT_Sup_SAVE, &
!$OMP     VLIDORT_LinFixIn_SAVE, VLIDORT_LinModIn_SAVE, VLIDORT_LinSup_SAVE, &
!$OMP     R2_PTSave,     Icorr_PTSave, &
!$OMP     R2_Baseline,   Icorr_Baseline, &
!$OMP     R2_LCBaseline, Icorr_LCBaseline, &
!$OMP     R2_LSBaseline, Icorr_LSBaseline, &
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
        write(*,'(1x,a,i1,a)') 'Running TwoOS driver with ',OMP_NTHREADS,' thread(s)'
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
!mick fix 4/21/2015 - added initialization
        VLIDORT_Out%Main%TS_STOKES = ZERO

        VLIDORT_LinFixIn = VLIDORT_LinFixIn_SAVE
        VLIDORT_LinModIn = VLIDORT_LinModIn_SAVE
        VLIDORT_LinSup   = VLIDORT_LinSup_SAVE
!mick fix 4/21/2015 - added initialization
        VLIDORT_LinOut%Col%TS_COLUMNWF   = ZERO
        VLIDORT_LinOut%Surf%TS_SURFACEWF = ZERO

!  linearization control

        VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = .false.
        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAY)   = .true.
        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAY) = 1
        VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = 0
        VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS         = 1
        VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 1

        VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .true.
        VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = .true.
        VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = .false.
        VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .true.

!  Linearization inputs. 1 Column WF, WRT Total Gas absorption

        do n = 1, nlay
          n1 = nlay + 1 - n ; coalb = one - ssa_save(n)
          VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(1,n1)  =  coalb
          VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(1,n1)  = -coalb
        enddo

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

        !write(*,*)
        !write(*,*)'Done 2OS LCS Baseline correction'

        write(charTID,'(i1)') TID
!mick fix 4/21/2015 - created new status subroutine for 2OS code 
        CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
          'TwoOS_self_LCS_OMPExecution_TID_' // charTID // '.log', &
          35+TID, OPENFILEFLAG, &
          VLIDORT_Out%Status )

        IF ( VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
          VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

        !!OMP_N_GEOMETRIES(wvn)  = VLIDORT_Out%Main%TS_n_geometries
        !R2_Baseline(1:3,wvn)   = VLIDORT_Out%Main%TS_2OSCORR_R2(1,1:3)
        !Icorr_Baseline (wvn)   = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)
        !R2_LCBaseline(1:3,wvn) = VLIDORT_LinOut%Col%TS_2OSCORR_LC_R2(1,1,1:3)
        !Icorr_LCBaseline (wvn) = VLIDORT_LinOut%Col%TS_2OSCORR_LC_ICORR(1,1)
        !R2_LSBaseline(1:3,wvn) = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,1:3)
        !Icorr_LSBaseline (wvn) = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_ICORR(1,1)

        do o1=1,3
          R2_Baseline(o1,wvn)   = VLIDORT_Out%Main%TS_2OSCORR_R2(1,o1)
          R2_LCBaseline(o1,wvn) = VLIDORT_LinOut%Col%TS_2OSCORR_LC_R2(1,1,o1)
          R2_LSBaseline(o1,wvn) = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,o1)
        enddo
        Icorr_Baseline(wvn)   = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)
        Icorr_LCBaseline(wvn) = VLIDORT_LinOut%Col%TS_2OSCORR_LC_ICORR(1,1)
        Icorr_LSBaseline(wvn) = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_ICORR(1,1)

!      stop'baseline Column WFS only'

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

      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .false.
      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = .false.
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .false.

      VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT = zero
      VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT = zero

!  TASKS for FD testing
!  ======================

!  Task 0, FD Perturbation on Surface Albedo
!  Task 1, FD Perturbation on trace gas total optical depth

      do task = 0, 1
!      do task = 1, 1

!  Task counter

         t = task

!  Task 0

         if ( task .eq. 0 ) then
            !write(*,*)'Doing task  0: '//'Finite Difference, perturb Lambertian albedo'
            VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =  spars_save * epsfac
            do n = 1, nlay
               n1 = nlay + 1 - n
               VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
               VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
            enddo
         endif

!  Task 1

         if ( task.gt.0  ) then
            !write(*,'(a,i2,a)')' Doing task ',task, &
            ! ': Finite Difference, perturb total molecular absorption '
            VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = spars_save
            do n = 1, nlay
               n1 = nlay + 1 - n
               scat = opd_save(n) * ssa_save(n)
               gaspert = ( opd_save(n) - scat ) * epsfac
               taupert = scat + gaspert
               VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = taupert
               VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = scat / taupert
            enddo
         endif

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

!  Exception handling, write-up (optional)

!mick fix 4/21/2015 - created new status subroutine for 2OS code 
         CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
           'TwoOS_self_LCS_OMPExecution_TID_' // charTID // '.log', &
           35+TID, OPENFILEFLAG, &
           VLIDORT_Out%Status )

         IF ( VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
           VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Save results

         !R2_PTSave(t,1:3,wvn) = VLIDORT_Out%Main%TS_2OSCORR_R2(1,1:3)
         !Icorr_PTSave (t,wvn) = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

         do o1=1,3
           R2_PTSave(t,o1,wvn) = VLIDORT_Out%Main%TS_2OSCORR_R2(1,o1)
         enddo
         Icorr_PTSave (t,wvn) = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

!  End task loop

      enddo

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

        OPEN(36,file='TwoOS_self_LCS_results.all' &
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

          OPEN(36,file='TwoOS_self_LCS_results.all' &
               // '_OMP' // trim(tail),status='replace')
        endif

        write(36,55)'Baseline R2/Icorr',R2_Baseline(1:3,wvn), Icorr_Baseline(wvn)
        do t = 0, 1
          if ( t.eq.0 ) then
            write(36,56)'Surface WF R2/Icorr ', &
                     spars_save*R2_LSBaseline(1:3,wvn), spars_save*Icorr_LSBaseline(wvn), &
                    (R2_PTSave(t,1,wvn)-R2_Baseline(1,wvn))/eps, &
                    (R2_PTSave(t,2,wvn)-R2_Baseline(2,wvn))/eps, &
                    (R2_PTSave(t,3,wvn)-R2_Baseline(3,wvn))/eps, &
                    (Icorr_PTSave(t,wvn)-Icorr_Baseline(wvn))/eps
          else
            write(36,56)'Column WF R2/Icorr  ', & 
                     R2_LCBaseline(1:3,wvn), Icorr_LCBaseline(wvn), &
                    (R2_PTSave(t,1,wvn)-R2_Baseline(1,wvn))/eps, &
                    (R2_PTSave(t,2,wvn)-R2_Baseline(2,wvn))/eps, &
                    (R2_PTSave(t,3,wvn)-R2_Baseline(3,wvn))/eps, &
                    (Icorr_PTSave(t,wvn)-Icorr_Baseline(wvn))/eps
          endif
        enddo

        if ( .not. do_std_output_file ) close(36)

!  End file write loop

      enddo

      if ( do_std_output_file ) close(36)

!     Deallocate output memory

      !DEALLOCATE ( OMP_N_GEOMETRIES, R2_PTSave, Icorr_PTSave, R2_Baseline, Icorr_Baseline, &
      !             R2_LCBaseline, Icorr_LCBaseline, R2_LSBaseline, Icorr_LSBaseline )
      DEALLOCATE ( R2_PTSave, Icorr_PTSave, R2_Baseline, Icorr_Baseline, &
                   R2_LCBaseline, Icorr_LCBaseline, R2_LSBaseline, Icorr_LSBaseline )

!  End tests loop

      enddo

!  Finish

      write(*,*)
      write(*,*) '******************************************'
      write(*,*)

      if ( VLIDORT_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) THEN
         stop 'Successful run, with Warning: look in ' // &
           'TwoOS_self_LPS_OMPExecution_TID_*.log or TwoOS_self_LCS_OMPExecution_TID_*.log'
      else
         stop 'Successful run'
      endif

      end program TWOOS_Self_LPCS_OMPTester
