      program TWOOS_Self_LPCS_Tester

!  This is a Self tester for the linearized 2OS code, validation against Vijay's latest.
!   25 february 2015, 3 March 2015

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

!  Dump variables

      integer          :: nlay,nmug,nmoms,nstokes,nfoumax,ncoefs(maxlayers)
      logical          :: regular_ps,enhanced_ps
      double precision :: epsilon,phi,emu0,emu,spars,xmu(maxstreams),w(maxstreams),&
                          Sun_Chapman2(maxlayers,maxlayers)
      double precision :: ssa(maxlayers),opd(maxlayers),coefs(0:maxmoments_input,maxlayers,6)

!  Saved results output

      double precision :: R2_Baseline(3),Icorr_Baseline
      double precision :: R2_LPBaseline(maxlayers,3),Icorr_LPBaseline(maxlayers)
      double precision :: R2_LCBaseline(3),Icorr_LCBaseline
      double precision :: R2_LSBaseline(3),Icorr_LSBaseline
      double precision :: R2_PTSave(0:62,3),Icorr_PTSave(0:62)

!  Flag for opening error output file

      LOGICAL ::          OPENFILEFLAG, PART1, PART2

!  local

      integer :: n, n1, k, l, nd, kd, ld, t, task, ks, o1
      double precision :: eps, epsfac, gaspert, scat, taupert, coalb
      double precision :: ssa_save(maxlayers),opd_save(maxlayers),spars_save

!  Start of code
!  =============

!  Define some control variables

!      PART1 = .FALSE.
      PART1 = .TRUE.

!      PART2 = .FALSE.
      PART2 = .TRUE.

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

      VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT = zero ! Set later

!  Linearized
!  ----------

      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG    = .false. ! Set later
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER  = 0       ! Set later
!mick fix 4/21/2015 - added these initializations
      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS  = 0       ! Set later
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = 0       ! Set later
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS      = 0       ! Set later
      VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS       = 0       ! Set later
      VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES     = ' '
      VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES    = ' '

      VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input    = zero ! Set later
      VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input = zero
      VLIDORT_LinFixIn%Optical%TS_L_omega_total_input    = zero ! Set later

!mick fix 4/21/2015 - added these initializations
      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = .false. ! Set later
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = .false. ! Set later
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .false. ! Set later
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .false. ! Set later
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = .false.
      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = .false.
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF            = .false.
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF          = .false.
      VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS            = .false.

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

!  Linearized (not required)

      VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC = ZERO
      VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0        = ZERO
      VLIDORT_LinSup%BRDF%TS_LS_BRDF_F          = ZERO
      VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0   = ZERO
      VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F     = ZERO
      VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY      = ZERO
      VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY = ZERO

      VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS = .FALSE.
      VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS   = 0

      VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC  = ZERO
      VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES = ZERO
      VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0        = ZERO
      VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0   = ZERO

!mick fix 4/21/2015 - added linearized SS & DB initializations
      VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_SS = ZERO
      VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_DB = ZERO
      VLIDORT_LinSup%SS%Col%TS_COLUMNWF_SS   = ZERO
      VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB   = ZERO
      VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB = ZERO
   
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

      VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR        = regular_ps
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING     = enhanced_ps
      VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL       = .not. regular_ps .and. .not. enhanced_ps

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

!  Optical

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

!  Must zero the Main output, since now intent(INOUT)

      VLIDORT_Out%Main%TS_STOKES       = zero
      VLIDORT_LinOut%Prof%TS_PROFILEWF = zero
      VLIDORT_LinOut%Surf%TS_SURFACEWF = zero

!  Set status flag

      OPENFILEFLAG = .false.

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

      write(*,*)'Done 2OS LPS Baseline correction'

!  Exception handling, write-up (optional)

!mick fix 4/21/2015 - created new status subroutine for 2OS code 
      CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
          'TwoOS_self_LPS_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save results

      !R2_Baseline   (1:3) = VLIDORT_Out%Main%TS_2OSCORR_R2(1,1:3)
      !Icorr_Baseline      = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)
      !do n = 1, nlay
      !   R2_LPBaseline   (n,1:3) = VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(1,n,1,1:3)
      !   Icorr_LPBaseline(n)     = VLIDORT_LinOut%Prof%TS_2OSCORR_LP_ICORR(1,n,1)
      !enddo
      !R2_LSBaseline   (1:3) = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,1:3)
      !Icorr_LSBaseline      = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_ICORR(1,1)

      do o1=1,3
        R2_Baseline(o1) = VLIDORT_Out%Main%TS_2OSCORR_R2(1,o1)
        do n = 1, nlay
          R2_LPBaseline(n,o1) = VLIDORT_LinOut%Prof%TS_2OSCORR_LP_R2(1,n,1,o1)
        enddo
        R2_LSBaseline(o1) = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,o1)
      enddo

      Icorr_Baseline = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)
      do n = 1, nlay
        Icorr_LPBaseline(n) = VLIDORT_LinOut%Prof%TS_2OSCORR_LP_ICORR(1,n,1)
      enddo
      Icorr_LSBaseline = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_ICORR(1,1)

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

!      stop'baseline profile WFS only'

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

      VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT = zero
      VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT = zero

!  Tasks for FD testing
!  ====================

!  Task 0, FD Perturbation on Surface Albedo
!  Task t, FD Perturbation on trace gas optical depth in layer t  (2OS = lowest layer)

      do task = 0, nlay
!      do task = 1, 1

!  task counter

         t = task

!  Task 0

         if ( task .eq. 0 ) then
            write(*,*)'Doing task  0: '//'Finite Difference, perturb Lambertian albedo'
            VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =   spars_save * epsfac
            do n = 1, nlay
               n1 = nlay + 1 - n
               VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
               VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
            enddo
         endif

!  Tasks 1 to NLAYERS

         if ( task.gt.0  ) then
            VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =   spars_save
            write(*,'(a,i2,a,i2)')' Doing task ',task, &
             ': Finite Difference, perturb molecular absorption Layer ',t
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
             'TwoOS_self_LPS_Execution.log', &
             35, OPENFILEFLAG, &
             VLIDORT_Out%Status )

!  Save results

         !R2_PTSave   (t,1:3) = VLIDORT_Out%Main%TS_2OSCORR_R2(1,1:3)
         !Icorr_PTSave(t)     = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

         do o1=1,3
           R2_PTSave(t,o1) = VLIDORT_Out%Main%TS_2OSCORR_R2(1,o1)
         enddo
         Icorr_PTSave(t) = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

!  End task loop

      enddo

!  write results

      open(36,file='TwoOS_self_LPS_results.all',status='replace')
      write(36,55)'Baseline R2/Icorr',R2_Baseline(1:3), Icorr_Baseline
      do t = 0, nlay
         if ( t.eq.0 ) then
            write(36,56)'Surface WF R2/Icorr          ', &
                     spars_save*R2_LSBaseline(1:3), spars_save*Icorr_LSBaseline, &
                    (R2_PTSave(t,1)-R2_Baseline(1))/eps, &
                    (R2_PTSave(t,2)-R2_Baseline(2))/eps, &
                    (R2_PTSave(t,3)-R2_Baseline(3))/eps, &
                    (Icorr_PTSave(t)-Icorr_Baseline)/eps
         else
            write(36,57)'Profile WF R2/Icorr, Layer ',t,& 
                     R2_LPBaseline(t,1:3), Icorr_LPBaseline(t), &
                    (R2_PTSave(t,1)-R2_Baseline(1))/eps, &
                    (R2_PTSave(t,2)-R2_Baseline(2))/eps, &
                    (R2_PTSave(t,3)-R2_Baseline(3))/eps, &
                    (Icorr_PTSave(t)-Icorr_Baseline)/eps
         endif
      enddo
      close(36)
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

!  Must zero the Main output, since now intent(INOUT)

      VLIDORT_Out%Main%TS_STOKES       = zero
      VLIDORT_LinOut%Col%TS_COLUMNWF   = zero
      VLIDORT_LinOut%Surf%TS_SURFACEWF = zero

!  Set status flag

      OPENFILEFLAG = .false.

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

      write(*,*)'Done 2OS LCS Baseline correction'

!  Exception handling, write-up (optional)

!mick fix 4/21/2015 - created new status subroutine for 2OS code 
      CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
          'TwoOS_self_LCS_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save results

      R2_LCBaseline  (1:3) = VLIDORT_LinOut%Col%TS_2OSCORR_LC_R2(1,1,1:3)
      Icorr_LCBaseline     = VLIDORT_LinOut%Col%TS_2OSCORR_LC_ICORR(1,1)
      R2_LSBaseline  (1:3) = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_R2(1,1,1:3)
      Icorr_LSBaseline     = VLIDORT_LinOut%Surf%TS_2OSCORR_LS_ICORR(1,1)
      R2_Baseline    (1:3) = VLIDORT_Out%Main%TS_2OSCORR_R2(1,1:3)
      Icorr_Baseline       = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

!  Write results (debug)
!      write(*,*)'LCS,R2(1) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,1),VLIDORT_Out%Main%TS_stokes(1,1,1,upidx)
!      write(*,*)'LCS,R2(2) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,2),VLIDORT_Out%Main%TS_stokes(1,1,2,upidx)
!      write(*,*)'LCS,R2(3) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,3),VLIDORT_Out%Main%TS_stokes(1,1,3,upidx)
!      write(*,*)'LCS,Icorr = ',VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

!      stop'baseline Column WFS only'

!  Initialize FD section
!  =====================

!  Fd perturbation

      eps = 1.0d-02
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

!  Tasks for FD testing
!  ====================

!  Task 0, FD Perturbation on Surface Albedo
!  Task 1, FD Perturbation on trace gas total optical depth

      do task = 0, 1
!      do task = 1, 1

!  task counter

         t = task

!  Task 1

         if ( task .eq. 0 ) then
            write(*,*)'Doing task  0: '//'Finite Difference, perturb Lambertian albedo'
            VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =   spars_save * epsfac
            do n = 1, nlay
               n1 = nlay + 1 - n
               VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
               VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
            enddo
         endif

!  Task 1

         if ( task.gt.0  ) then
            VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =   spars_save
            write(*,'(a,i2,a)')' Doing task ',task, &
             ': Finite Difference, perturb total molecular absorption '
            do n = 1, nlay
               n1 = nlay + 1 - n
               scat = opd_save(n) * ssa_save(n)
               gaspert = ( opd_save(n) - scat ) * epsfac
               taupert = scat + gaspert
               VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = taupert
               VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = scat / taupert
            enddo
         endif

!  Must zero the Main output, since now intent(INOUT)

         VLIDORT_Out%Main%TS_STOKES       = zero

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
             'TwoOS_self_LCS_Execution.log', &
             35, OPENFILEFLAG, &
             VLIDORT_Out%Status )

!  Save results

         R2_PTSave   (t,1:3) = VLIDORT_Out%Main%TS_2OSCORR_R2(1,1:3)
         Icorr_PTSave(t)     = VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

!  write results debug
!      write(*,*)'T-pert,R2(1) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,1),VLIDORT_Out%Main%TS_stokes(1,1,1,upidx)
!      write(*,*)'T-pert,R2(2) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,2),VLIDORT_Out%Main%TS_stokes(1,1,2,upidx)
!      write(*,*)'T-pert,R2(3) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,3),VLIDORT_Out%Main%TS_stokes(1,1,3,upidx)
!      write(*,*)'T-oert,Icorr = ',VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

!  End task loop

      enddo

!  write results

      open(36,file='TwoOS_self_LCS_results.all',status='replace')
      write(36,55)'Baseline R2/Icorr',R2_Baseline(1:3), Icorr_Baseline
      do t = 0, 1
         if ( t.eq.0 ) then
            write(36,56)'Surface WF R2/Icorr ',&
                     spars_save*R2_LSBaseline(1:3), spars_save*Icorr_LSBaseline, &
                    (R2_PTSave(t,1)-R2_Baseline(1))/eps, &
                    (R2_PTSave(t,2)-R2_Baseline(2))/eps, &
                    (R2_PTSave(t,3)-R2_Baseline(3))/eps, &
                    (Icorr_PTSave(t)-Icorr_Baseline)/eps
         else
            write(36,56)'Column WF R2/Icorr  ',& 
                     R2_LCBaseline(1:3), Icorr_LCBaseline, &
                    (R2_PTSave(t,1)-R2_Baseline(1))/eps, &
                    (R2_PTSave(t,2)-R2_Baseline(2))/eps, &
                    (R2_PTSave(t,3)-R2_Baseline(3))/eps, &
                    (Icorr_PTSave(t)-Icorr_Baseline)/eps
         endif
      enddo
      close(36)

!  Finish

      write(*,*)
      if (VLIDORT_Out%Status%TS_2OSCORR_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) THEN
         stop 'Successful run, with Warning: look in TwoOS_self_LPS_Execution.log or TwoOS_self_LCS_Execution.log'
      else
         stop 'Successful run'
      endif

      end program TWOOS_Self_LPCS_Tester
