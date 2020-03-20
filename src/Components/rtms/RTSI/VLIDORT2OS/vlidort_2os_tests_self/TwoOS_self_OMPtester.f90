      program TWOOS_OMP_Self_Tester

!  This is a Self tester for the 2OS code, validation against Vijay's latest.
!   25 february 2015, 3 March 2015 (OMP testing)

!  Module files for VLIDORT

      USE VLIDORT_PARS
      USE VLIDORT_IO_DEFS

      USE VLIDORT_AUX

      USE vlidort_2OScorr_master_m
!mick fix 4/21/2015 - added new 2OS module
      USE vlidort_2OS_corr_aux

      IMPLICIT NONE

!  Dimensioning integer parameters

      INTEGER, PARAMETER :: MAXWVN = 100!10!4!50!100!1000

!  VLIDORT file inputs status structure

      TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn, VLIDORT_FixIn_SAVE
      TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn, VLIDORT_ModIn_SAVE

!  VLIDORT supplements i/o structure

      TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup, VLIDORT_Sup_SAVE

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

!  Local Variables
!  ===============

!  Dump variables

      integer :: nlay,nmug,nmoms,nstokes,nfoumax,ncoefs(maxlayers)
      logical :: regular_ps,enhanced_ps
      double precision :: epsilon,phi,emu0,emu,spars,xmu(maxstreams),w(maxstreams), &
                          Sun_Chapman2(maxlayers,maxlayers)
      double precision :: ssa(maxlayers),opd(maxlayers),coefs(0:maxmoments_input,maxlayers,6)

!  Saved results

      !INTEGER, DIMENSION(:), ALLOCATABLE :: OMP_N_GEOMETRIES
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: R2_PT
      DOUBLE PRECISION, DIMENSION(:,:)  , ALLOCATABLE :: Icorr_PT

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

      LOGICAL :: OPENFILEFLAG

!  flag for standard output file

      !LOGICAL, parameter :: do_std_output_file = .true.
      LOGICAL, parameter :: do_std_output_file = .false.

!  local

      integer :: n, n1, k, l, nd, kd, ld, o1

!  Start of code
!  =============

!  Set test parameters

      N_CORE = 2!4
      N_OMP_TESTS = 2!3

      numwvn = maxwvn

!  Begin test loop

      DO TEST = 1, N_OMP_TESTS
      !DO TEST = 3, N_OMP_TESTS
      !DO TEST = 1, 1
      !DO TEST = 2, 2

      write(*,*)
      write(*,*) '******************************************'
      write(*,'(1x,a,i1)') 'Doing test ',TEST

!  Set total number of threads (optional)
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

!  Hardwired input

!  Fixed Boolean inputs
!  --------------------

      VLIDORT_FixIn_SAVE%Bool%TS_DO_FULLRAD_MODE      = .false.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_SSCORR_TRUNCATION = .false.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_SS_EXTERNAL       = .false.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_SSFULL            = .false.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_THERMAL_EMISSION  = .true.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_SURFACE_EMISSION  = .true.

      VLIDORT_FixIn_SAVE%Bool%TS_DO_PLANE_PARALLEL      = .false.  ! Set later
      VLIDORT_FixIn_SAVE%Bool%TS_DO_LAMBERTIAN_SURFACE  = .true.

      VLIDORT_FixIn_SAVE%Bool%TS_DO_UPWELLING           = .true.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_DNWELLING           = .false.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_SURFACE_LEAVING     = .false.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_SL_ISOTROPIC        = .false.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_QUAD_OUTPUT         = .false.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_TOA_CONTRIBS        = .false.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_SPECIALIST_OPTION_1 = .false.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_SPECIALIST_OPTION_2 = .false.
      VLIDORT_FixIn_SAVE%Bool%TS_DO_SPECIALIST_OPTION_3 = .false.

!  Fixed control
!  -------------

      VLIDORT_FixIn_SAVE%Cont%TS_NSTOKES          = 0 ! Set later
      VLIDORT_FixIn_SAVE%Cont%TS_NSTREAMS         = 0 ! set later
      VLIDORT_FixIn_SAVE%Cont%TS_NLAYERS          = 0 ! set later
      VLIDORT_FixIn_SAVE%Cont%TS_NFINELAYERS      = 0
      VLIDORT_FixIn_SAVE%Cont%TS_N_THERMAL_COEFFS = 0

      VLIDORT_FixIn_SAVE%Cont%TS_TAYLOR_ORDER     = 0

      VLIDORT_FixIn_SAVE%Cont%TS_VLIDORT_ACCURACY = zero ! Set later

      VLIDORT_FixIn_SAVE%Cont%TS_NLAYERS_NOMS   = 0
      VLIDORT_FixIn_SAVE%Cont%TS_NLAYERS_CUTOFF = 0

!  Fixed Beam Userval inputs
!  --------------------------

      VLIDORT_FixIn_SAVE%Sunrays%TS_FLUX_FACTOR     = ONE
      VLIDORT_FixIn_SAVE%UserVal%TS_N_USER_LEVELS   = 1

!  Chapman inputs
!  --------------

!mick fix 4/21/2015 - added HEIGHT_GRID initialization
      VLIDORT_FixIn_SAVE%Chapman%TS_HEIGHT_GRID       = ZERO
      VLIDORT_FixIn_SAVE%Chapman%TS_PRESSURE_GRID     = ZERO
      VLIDORT_FixIn_SAVE%Chapman%TS_TEMPERATURE_GRID  = ZERO
      VLIDORT_FixIn_SAVE%Chapman%TS_FINEGRID          = 0
      VLIDORT_FixIn_SAVE%Chapman%TS_RFINDEX_PARAMETER = Zero

!  Fixed Optical inputs
!  --------------------

      VLIDORT_FixIn_SAVE%Optical%TS_THERMAL_BB_INPUT = ZERO
      VLIDORT_FixIn_SAVE%Optical%TS_SURFACE_BB_INPUT = ZERO

      VLIDORT_FixIn_SAVE%Optical%TS_GREEKMAT_TOTAL_INPUT = zero ! Set later
      VLIDORT_FixIn_SAVE%Optical%TS_DELTAU_VERT_INPUT    = zero ! Set later
      VLIDORT_FixIn_SAVE%Optical%TS_LAMBERTIAN_ALBEDO    = zero ! Set later

!  Fixed write
!  -----------

      VLIDORT_FixIn_SAVE%Write%TS_DO_DEBUG_WRITE    = .false.
      VLIDORT_FixIn_SAVE%Write%TS_DO_WRITE_INPUT    = .false.
      VLIDORT_FixIn_SAVE%Write%TS_DO_WRITE_SCENARIO = .false.
      VLIDORT_FixIn_SAVE%Write%TS_DO_WRITE_FOURIER  = .false.
      VLIDORT_FixIn_SAVE%Write%TS_DO_WRITE_RESULTS  = .false.

      VLIDORT_FixIn_SAVE%Write%TS_INPUT_WRITE_FILENAME     = ' '
      VLIDORT_FixIn_SAVE%Write%TS_SCENARIO_WRITE_FILENAME  = ' '
      VLIDORT_FixIn_SAVE%Write%TS_FOURIER_WRITE_FILENAME   = ' '
      VLIDORT_FixIn_SAVE%Write%TS_RESULTS_WRITE_FILENAME   = ' '

!  Modified Booleans
!  -----------------

      VLIDORT_ModIn_SAVE%MBool%TS_DO_SOLAR_SOURCES   = .true.

!  Fixed defaults
!    -- Full calculation with outgoing sphericity single-scatter

      VLIDORT_ModIn_SAVE%MBool%TS_DO_SSCORR_NADIR        = .false. !  Set later
      VLIDORT_ModIn_SAVE%MBool%TS_DO_SSCORR_OUTGOING     = .false. !  Set later
!mick fix 4/21/2015 - added DO_FO_CALC initialization
      VLIDORT_ModIn_SAVE%MBool%TS_DO_FO_CALC             = .false.

      VLIDORT_ModIn_SAVE%MBool%TS_DO_REFRACTIVE_GEOMETRY = .false.
      VLIDORT_ModIn_SAVE%MBool%TS_DO_CHAPMAN_FUNCTION    = .false.   ! Not required

      VLIDORT_ModIn_SAVE%MBool%TS_DO_USER_VZANGLES       = .true.
      VLIDORT_ModIn_SAVE%MBool%TS_DO_ADDITIONAL_MVOUT    = .false.
      VLIDORT_ModIn_SAVE%MBool%TS_DO_MVOUT_ONLY          = .false.
      VLIDORT_ModIn_SAVE%MBool%TS_DO_THERMAL_TRANSONLY   = .false.

      VLIDORT_ModIn_SAVE%MBool%TS_DO_SOLUTION_SAVING     = .false.
      VLIDORT_ModIn_SAVE%MBool%TS_DO_BVP_TELESCOPING     = .false.

!  From the input (Observation Geometry is new)

      VLIDORT_ModIn_SAVE%MBool%TS_DO_OBSERVATION_GEOMETRY = .true.

      VLIDORT_ModIn_SAVE%MBool%TS_DO_RAYLEIGH_ONLY        = .false.
      VLIDORT_ModIn_SAVE%MBool%TS_DO_DOUBLE_CONVTEST      = .true.
      VLIDORT_ModIn_SAVE%MBool%TS_DO_DELTAM_SCALING       = .false. ! Status uncertain

      VLIDORT_ModIn_SAVE%MBool%TS_DO_2OS_CORRECTION       = .true.  ! Must be set

!  Modified Control
!  ----------------

      VLIDORT_ModIn_SAVE%MCont%TS_NGREEK_MOMENTS_INPUT = 0 ! Set later

!  Modified Sunrays
!  ----------------

      VLIDORT_ModIn_SAVE%MSunrays%TS_N_SZANGLES        = 1
      VLIDORT_ModIn_SAVE%MSunrays%TS_SZANGLES          = zero ! set later

!  Modified Uservals
!  -----------------

      VLIDORT_ModIn_SAVE%MUserVal%TS_N_USER_VZANGLES      = 1
      VLIDORT_ModIn_SAVE%MUserVal%TS_USER_VZANGLES_INPUT  = zero ! set later
      VLIDORT_ModIn_SAVE%MUserVal%TS_N_USER_RELAZMS       = 1
      VLIDORT_ModIn_SAVE%MUserVal%TS_USER_RELAZMS         = zero  ! Set later
!mick fix 4/21/2015 - initialize all elements of USER_LEVELS
      !VLIDORT_ModIn_SAVE%MUserVal%TS_USER_LEVELS(1)       = zero
      VLIDORT_ModIn_SAVE%MUserVal%TS_USER_LEVELS          = zero
      VLIDORT_ModIn_SAVE%MUserVal%TS_GEOMETRY_SPECHEIGHT  = zero
      VLIDORT_ModIn_SAVE%MUserVal%TS_N_USER_OBSGEOMS      = 1
      VLIDORT_ModIn_SAVE%MUserVal%TS_USER_OBSGEOMS_INPUT  = zero !  set later

!  Modified Chapman inputs
!  -----------------------

      VLIDORT_ModIn_SAVE%MChapman%TS_EARTH_RADIUS        = zero
      VLIDORT_ModIn_SAVE%MChapman%TS_CHAPMAN_FACTORS     = zero ! Set later

!  Modified Optical
!  ----------------

      VLIDORT_ModIn_SAVE%MOptical%TS_OMEGA_TOTAL_INPUT   = zero ! Set later

! Supplemental Initialization
! ---------------------------

!  BRDF contributions, may be set later, in CODEC_VLIDORT_Optical

      VLIDORT_Sup_SAVE%BRDF%TS_EXACTDB_BRDFUNC  = ZERO
      VLIDORT_Sup_SAVE%BRDF%TS_BRDF_F_0         = ZERO
      VLIDORT_Sup_SAVE%BRDF%TS_BRDF_F           = ZERO
      VLIDORT_Sup_SAVE%BRDF%TS_USER_BRDF_F_0    = ZERO
      VLIDORT_Sup_SAVE%BRDF%TS_USER_BRDF_F      = ZERO

!  Emissitivity may be set later, in CODEC_VLIDORT_Optical

      VLIDORT_Sup_SAVE%BRDF%TS_EMISSIVITY       = ZERO
      VLIDORT_Sup_SAVE%BRDF%TS_USER_EMISSIVITY  = ZERO

!  Surface leaving stuff may be set later, in CODEC_VLIDORT_Optical

      VLIDORT_Sup_SAVE%SLEAVE%TS_SLTERM_ISOTROPIC  = ZERO
      VLIDORT_Sup_SAVE%SLEAVE%TS_SLTERM_USERANGLES = ZERO
      VLIDORT_Sup_SAVE%SLEAVE%TS_SLTERM_F_0        = ZERO
      VLIDORT_Sup_SAVE%SLEAVE%TS_USER_SLTERM_F_0   = ZERO

!mick fix 4/21/2015 - added SUP initializations
!  SS and DB stuff may be set later

      VLIDORT_Sup_SAVE%SS%TS_STOKES_SS  = ZERO
      VLIDORT_Sup_SAVE%SS%TS_STOKES_DB  = ZERO
      VLIDORT_Sup_SAVE%SS%TS_FO_ZMATRIX = ZERO
      VLIDORT_Sup_SAVE%SS%TS_FO_R1SAVED = ZERO

!  Now open the files

      open(87,file='fort.87_Rg',status='old')
      read(87,*)nlay,nmug,nstokes,nmoms,nfoumax,regular_ps,enhanced_ps
      read(87,*)epsilon,phi,emu0,emu,spars
      do k = 1, nmug
         read(87,*)kd,xmu(k),w(k)
      enddo
      do n = 1, nlay
         read(87,*)ssa(n),opd(n),ncoefs(n)
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

      VLIDORT_ModIn_SAVE%MBool%TS_DO_SSCORR_NADIR    = regular_ps
      VLIDORT_ModIn_SAVE%MBool%TS_DO_SSCORR_OUTGOING = enhanced_ps
      VLIDORT_FixIn_SAVE%Bool%TS_DO_PLANE_PARALLEL   = .not. regular_ps .and. .not. enhanced_ps

      VLIDORT_FixIn_SAVE%Cont%TS_NSTOKES          = nstokes
      VLIDORT_FixIn_SAVE%Cont%TS_NSTREAMS         = nmug
      VLIDORT_FixIn_SAVE%Cont%TS_NLAYERS          = nlay
      VLIDORT_FixIn_SAVE%Cont%TS_VLIDORT_ACCURACY = epsilon

      VLIDORT_ModIn_SAVE%MCont%TS_NGREEK_MOMENTS_INPUT = nmoms

      VLIDORT_ModIn_SAVE%MSunrays%TS_SZANGLES             = acos(emu0)/deg_to_rad
      VLIDORT_ModIn_SAVE%MUserVal%TS_USER_VZANGLES_INPUT  = acos(emu) /deg_to_rad
      VLIDORT_ModIn_SAVE%MUserVal%TS_USER_RELAZMS         = phi

      VLIDORT_ModIn_SAVE%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1) = acos(emu0)/deg_to_rad
      VLIDORT_ModIn_SAVE%MUserVal%TS_USER_OBSGEOMS_INPUT(1,2) = acos(emu) /deg_to_rad
      VLIDORT_ModIn_SAVE%MUserVal%TS_USER_OBSGEOMS_INPUT(1,3) = phi

!mick fix 4/20/2015 - trimmed argument passing
      !VLIDORT_ModIn_SAVE%MChapman%TS_CHAPMAN_FACTORS(:,:,1)   = Sun_Chapman2(:,:)
      VLIDORT_ModIn_SAVE%MChapman%TS_CHAPMAN_FACTORS(1:nlay,1:nlay,1) = &
        Sun_Chapman2(1:nlay,1:nlay)

      do n = 1, nlay
        n1 = nlay + 1 - n
        VLIDORT_FixIn_SAVE%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd(n)
        VLIDORT_ModIn_SAVE%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa(n)
        VLIDORT_FixIn_SAVE%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,1)  =   coefs(0:nmoms,n,1)
        VLIDORT_FixIn_SAVE%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,2)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn_SAVE%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,5)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn_SAVE%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,6)  =   coefs(0:nmoms,n,2)
        VLIDORT_FixIn_SAVE%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,11) =   coefs(0:nmoms,n,3)
        VLIDORT_FixIn_SAVE%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,12) =   coefs(0:nmoms,n,6)
        VLIDORT_FixIn_SAVE%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,15) = - coefs(0:nmoms,n,6)
        VLIDORT_FixIn_SAVE%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,16) =   coefs(0:nmoms,n,4)
      enddo
      VLIDORT_FixIn_SAVE%Optical%TS_LAMBERTIAN_ALBEDO =   spars

!     Only master thread does this

      !OMP_NTHREADS = OMP_GET_NUM_THREADS()
      !write(*,*)
      !write(*,'(1x,a,i1)') 'Total number of threads (to start) = ', OMP_NTHREADS

!     Allocate output memory

      ALLOCATE ( & 
        !OMP_N_GEOMETRIES ( MAXWVN ), &
        R2_PT            ( MAX_GEOMETRIES, MAXSTOKES, MAXWVN ), &
        Icorr_PT         ( MAX_GEOMETRIES, MAXWVN ) )

!     Initialize output

        !OMP_N_GEOMETRIES = 0
        R2_PT            = ZERO
        Icorr_PT         = ZERO

!  Initialize error file output flag

        OPENFILEFLAG = .false.

!  Initialize shared input check variable

        VLIDORT_STATUS_INPUTCHECK = 0

!     Begin parallel region

!mick fix 4/21/2015 -
!Added:   VLIDORT_STATUS_INPUTCHECK
!Removed: OMP_N_GEOMETRIES

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (OMP_NTHREADS, &
!$OMP     NUMWVN, &
!$OMP     VLIDORT_FixIn_SAVE, VLIDORT_ModIn_SAVE, VLIDORT_Sup_SAVE, &
!$OMP     R2_PT, Icorr_PT, &
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

      !write(*,'(1x,a,i1,a,i5)') 'Thread  = ',TID,' wvn = ',wvn

!  Initialize some inputs

      VLIDORT_FixIn = VLIDORT_FixIn_SAVE
      VLIDORT_ModIn = VLIDORT_ModIn_SAVE
      VLIDORT_Sup   = VLIDORT_Sup_SAVE

!mick fix 4/21/2015 - added initialization
!  Must zero the Main output, since now intent(INOUT)

      VLIDORT_Out%Main%TS_STOKES = ZERO

!  Call 2OS correction routine

      CALL VLIDORT_2OSCORR_MASTER ( &
           VLIDORT_FixIn, &
           VLIDORT_ModIn, &
           VLIDORT_Sup,   &
           VLIDORT_Out )

      !write(*,*)
      !write(*,*)'Done 2OS correction'

!  Exception handling, write-up (optional)

        write(charTID,'(i1)') TID
!mick fix 4/21/2015 - created new status subroutine for 2OS code 
        CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
          'TwoOS_self_OMPExecution_TID_' // charTID // '.log', &
          35+TID, OPENFILEFLAG, &
          VLIDORT_Out%Status )

        IF ( VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) &
          VLIDORT_STATUS_INPUTCHECK = VLIDORT_WARNING

!  Saved results

        !OMP_N_GEOMETRIES(wvn) = VLIDORT_Out%Main%TS_n_geometries
        !R2_PT(1,1:3,wvn)      = VLIDORT_Out%Main%TS_2OSCORR_R2(1,1:3)
        !Icorr_PT(1,wvn)       = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

        do o1=1,3
          R2_PT(1,o1,wvn)      = VLIDORT_Out%Main%TS_2OSCORR_R2(1,o1)
        end do
        Icorr_PT(1,wvn)       = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

!  Write results

!      write(*,*)'R2(1) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,1),VLIDORT_Out%Main%TS_stokes(1,1,1,upidx)
!      write(*,*)'R2(2) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,2),VLIDORT_Out%Main%TS_stokes(1,1,2,upidx)
!      write(*,*)'R2(3) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,3),VLIDORT_Out%Main%TS_stokes(1,1,3,upidx)
!      write(*,*)'Icorr = ',VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

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

        OPEN(36,file='TwoOS_self_results.all' &
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

          OPEN(36,file='TwoOS_self_results.all' &
               // '_OMP' // trim(tail),status='replace')
        endif

!        write(36,*)wvn,R2_PT(1,1:3,wvn),Icorr_PT(1,wvn)
        write(36,77)'Baseline R2(1) = ',R2_PT(1,1,wvn)
        write(36,77)'Baseline R2(2) = ',R2_PT(1,2,wvn)
        write(36,77)'Baseline R2(3) = ',R2_PT(1,3,wvn)
        write(36,78)'Baseline Icorr = ',Icorr_PT(1,wvn)
        close(36)
77      format(A,1p2e17.8)
78      format(A,1pe17.8)

        if ( .not. do_std_output_file ) close(36)

!  End file write loop

      enddo

      if ( do_std_output_file ) close(36)

!     Deallocate output memory

      !DEALLOCATE ( OMP_N_GEOMETRIES, R2_PT, Icorr_PT )
      DEALLOCATE ( R2_PT, Icorr_PT )

!  End test loop

      enddo

!  Finish

      write(*,*)
      write(*,*) '******************************************'
      write(*,*)

      if ( VLIDORT_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) THEN
         stop 'Successful run, with Warning: look in TwoOS_self_OMPExecution_TID_*.log'
      else
         stop 'Successful run'
      endif

      end program TWOOS_OMP_Self_Tester
