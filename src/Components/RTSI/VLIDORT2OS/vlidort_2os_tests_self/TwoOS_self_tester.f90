      program TWOOS_Self_Tester

!  This is a Self tester for the 2OS code, validation against Vijay's latest.
!   25 february 2015

!  Module files for VLIDORT

      USE VLIDORT_PARS
      USE VLIDORT_IO_DEFS

      USE VLIDORT_AUX

      USE vlidort_2OScorr_master_m
!mick fix 4/21/2015 - added new 2OS module
      USE vlidort_2OS_corr_aux

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

!  Local Variables
!  ===============

!  Dump variables

      integer :: nlay,nmug,nmoms,nstokes,nfoumax,ncoefs(maxlayers)
      logical :: regular_ps,enhanced_ps
      double precision :: epsilon,phi,emu0,emu,spars,xmu(maxstreams),w(maxstreams),&
                          Sun_Chapman2(maxlayers,maxlayers)
      double precision :: ssa(maxlayers),opd(maxlayers),&
                          coefs(0:maxmoments_input,maxlayers,6)

!  Flag for opening error output file

      LOGICAL ::          OPENFILEFLAG

!  local

      integer :: n, n1, k, l, nd, kd, ld

!  Start of code
!  =============

!  Hardwired input

!  Fixed Boolean inputs
!  --------------------

      VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE      = .false.
      VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION = .false.
      VLIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL       = .false.
      VLIDORT_FixIn%Bool%TS_DO_SSFULL            = .false.
      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION  = .true.
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION  = .true.

      VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL      = .false.  ! Set later
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE  = .true.

      VLIDORT_FixIn%Bool%TS_DO_UPWELLING           = .true.
      VLIDORT_FixIn%Bool%TS_DO_DNWELLING           = .false.
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING     = .false.
      VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC        = .false.
      VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT         = .false.
      VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS        = .false.
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_1 = .false.
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_2 = .false.
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_3 = .false.

!  Fixed control
!  -------------

      VLIDORT_FixIn%Cont%TS_NSTOKES          = 0 ! Set later
      VLIDORT_FixIn%Cont%TS_NSTREAMS         = 0 ! set later
      VLIDORT_FixIn%Cont%TS_NLAYERS          = 0 ! set later
      VLIDORT_FixIn%Cont%TS_NFINELAYERS      = 0
      VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = 0

      VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER     = 0

      VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY  = zero ! Set later

      VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS   = 0
      VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF = 0

!  Fixed Beam Userval inputs
!  --------------------------

      VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR     = ONE
      VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS   = 1

!  Chapman inputs
!  --------------

!mick fix 4/21/2015 - added HEIGHT_GRID initialization
      VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID       = ZERO
      VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID     = ZERO
      VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID  = ZERO
      VLIDORT_FixIn%Chapman%TS_FINEGRID          = 0
      VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = Zero

!  Fixed Optical inputs
!  --------------------

      VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT = ZERO
      VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT = ZERO

      VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT = zero ! Set later
      VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT    = zero ! Set later
      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO    = zero ! Set later

!  Fixed write
!  -----------

      VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE    = .false.
      VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT    = .false.
      VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO = .false.
      VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER  = .false.
      VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS  = .false.

      VLIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME     = ' '
      VLIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME  = ' '
      VLIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME   = ' '
      VLIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME   = ' '

!  Modified Booleans
!  -----------------

      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES   = .true.

!  Fixed defaults
!    -- Full calculation with outgoing sphericity single-scatter

      VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR        = .false. !  Set later
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING     = .false. !  Set later
!mick fix 4/21/2015 - added DO_FO_CALC initialization
      VLIDORT_ModIn%MBool%TS_DO_FO_CALC             = .false.

      VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = .false.
      VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = .false.   ! Not required

      VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES       = .true.
      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = .false.
      VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY          = .false.
      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY   = .false.

      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING     = .false.
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING     = .false.

!  From the input (Observation Geometry is new)

      VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = .true.

      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY        = .false.
      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST      = .true.
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING       = .false. ! Status uncertain

      VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION       = .true.  ! Must be set

!  Modified Control
!  ----------------

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 0 ! Set later

!  Modified Sunrays
!  ----------------

      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES        = 1
      VLIDORT_ModIn%MSunrays%TS_SZANGLES          = zero ! set later

!  Modified Uservals
!  -----------------

      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES      = 1
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT  = zero ! set later
      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS       = 1
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS         = zero  ! Set later
!mick fix 4/21/2015 - initialize all elements of USER_LEVELS
      !VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1)       = zero
      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS          = zero
      VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT  = zero
      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS      = 1
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT  = zero !  set later

!  Modified Chapman inputs
!  -----------------------

      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS        = zero
      VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS     = zero ! Set later

!  Modified Optical
!  ----------------

      VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT   = zero ! Set later

! Supplemental Initialization
! ---------------------------

!  BRDF contributions, may be set later, in CODEC_VLIDORT_Optical

      VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC  = ZERO
      VLIDORT_Sup%BRDF%TS_BRDF_F_0         = ZERO
      VLIDORT_Sup%BRDF%TS_BRDF_F           = ZERO
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0    = ZERO
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F      = ZERO

!  Emissitivity may be set later, in CODEC_VLIDORT_Optical

      VLIDORT_Sup%BRDF%TS_EMISSIVITY       = ZERO
      VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY  = ZERO

!  Surface leaving stuff may be set later, in CODEC_VLIDORT_Optical

      VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC  = ZERO
      VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES = ZERO
      VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0        = ZERO
      VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0   = ZERO

!mick fix 4/21/2015 - added SS & DB initializations
!  SS and DB stuff may be set later

      VLIDORT_Sup%SS%TS_STOKES_SS = ZERO
      VLIDORT_Sup%SS%TS_STOKES_DB = ZERO

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

      VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR    = regular_ps
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING = enhanced_ps
      VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL   = .not. regular_ps .and. .not. enhanced_ps

      VLIDORT_FixIn%Cont%TS_NSTOKES          = nstokes
      VLIDORT_FixIn%Cont%TS_NSTREAMS         = nmug
      VLIDORT_FixIn%Cont%TS_NLAYERS          = nlay
      VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY = epsilon

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = nmoms

      VLIDORT_ModIn%MSunrays%TS_SZANGLES             = acos(emu0)/deg_to_rad
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT  = acos(emu) /deg_to_rad
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS         = phi

      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1) = acos(emu0)/deg_to_rad
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,2) = acos(emu) /deg_to_rad
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,3) = phi

!mick fix 4/20/2015 - trimmed argument passing
      !VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS(:,:,1) = Sun_Chapman2(:,:)
      VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS(1:nlay,1:nlay,1) = &
        Sun_Chapman2(1:nlay,1:nlay)

      do n = 1, nlay
        n1 = nlay + 1 - n
        VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd(n)
        VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa(n)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,1)  =   coefs(0:nmoms,n,1)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,2)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,5)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,6)  =   coefs(0:nmoms,n,2)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,11) =   coefs(0:nmoms,n,3)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,12) =   coefs(0:nmoms,n,6)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,15) = - coefs(0:nmoms,n,6)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,16) =   coefs(0:nmoms,n,4)
      enddo
      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =   spars

!  Must zero the Main output, since now intent(INOUT)

      VLIDORT_Out%Main%TS_STOKES = ZERO

!  Set status flag

      OPENFILEFLAG = .false.

!  Call 2OS correction routine

      CALL VLIDORT_2OSCORR_MASTER ( &
           VLIDORT_FixIn, &
           VLIDORT_ModIn, &
           VLIDORT_Sup,   &
           VLIDORT_Out )

      write(*,*)
      write(*,*)'Done 2OS correction'

!  Exception handling, write-up (optional)

!mick fix 4/21/2015 - created new status subroutine for 2OS code 
      CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
          'TwoOS_self_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Write results

      open(36,file='TwoOS_self_results.all',status='replace')
      write(36,77)'Baseline R2(1) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,1),VLIDORT_Out%Main%TS_stokes(1,1,1,upidx)
      write(36,77)'Baseline R2(2) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,2),VLIDORT_Out%Main%TS_stokes(1,1,2,upidx)
      write(36,77)'Baseline R2(3) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,3),VLIDORT_Out%Main%TS_stokes(1,1,3,upidx)
      write(36,78)'Baseline Icorr = ',VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)
      close(36)
77    format(A,1p2e17.8)
78    format(A,1pe17.8)

!  Finish

      write(*,*)
      if (VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) THEN
         stop 'Successful run, with Warning: look in TwoOS_self_Execution.log'
      else
         stop 'Successful run'
      endif

      end program TWOOS_Self_Tester
