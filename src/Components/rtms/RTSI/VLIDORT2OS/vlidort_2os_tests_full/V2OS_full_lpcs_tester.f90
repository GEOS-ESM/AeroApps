      program V2OS_Full_LPCS_Tester

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

      DOUBLE PRECISION :: STOKES_V_BAS ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: STOKES_S_BAS ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: STOKES_C_BAS ( MAX_GEOMETRIES, MAXSTOKES )
!mick mod 5/1/2015 - add Icorr
      DOUBLE PRECISION :: Icorr ( MAX_GEOMETRIES )

      DOUBLE PRECISION :: PROFILEWF_V_BAS ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: PROFILEWF_S_BAS ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: PROFILEWF_C_BAS ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION :: SURFACEWF_V_BAS ( MAX_SURFACEWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: SURFACEWF_S_BAS ( MAX_SURFACEWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: SURFACEWF_C_BAS ( MAX_SURFACEWFS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION :: COLUMNWF_V_BAS ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: COLUMNWF_S_BAS ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: COLUMNWF_C_BAS ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXSTOKES )

!  Saved Perturbation results. Note the number of tasks!

      DOUBLE PRECISION :: STOKES_V_PT ( MAX_GEOMETRIES, MAXSTOKES, 0:70 )
      DOUBLE PRECISION :: STOKES_S_PT ( MAX_GEOMETRIES, MAXSTOKES, 0:70 )
      DOUBLE PRECISION :: STOKES_C_PT ( MAX_GEOMETRIES, MAXSTOKES, 0:70 )

!  Debug

      LOGICAL, parameter :: debug = .true.
      !LOGICAL, parameter :: debug = .false.

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

      VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE      = .true.    ! V2OS Full - Must be TRUE for the Full test
      VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION = .false.
      VLIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL       = .false.
      VLIDORT_FixIn%Bool%TS_DO_SSFULL            = .false.
      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION  = .false.!.true.
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION  = .false.!.true.

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
      VLIDORT_FixIn%Cont%TS_NFINELAYERS      = 2 ! V2OS Full - Should be set for the Full test
      VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = 0!2 ! V2OS Full - Should be set for the Full test

      VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER     = 3 ! V2OS Full - Should be set for the Full test

      VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY  = zero ! Set later

      VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS   = 0
      VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF = 0

!  Fixed Beam Userval inputs
!  --------------------------

      VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR     = ONE
      VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS   = 1

!  Chapman inputs
!  --------------

      VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID      = ZERO ! V2OS Full - must be assigned from input (see later)

      VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID    = ZERO
      VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID = ZERO
      VLIDORT_FixIn%Chapman%TS_FINEGRID         = 0
      VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = ZERO

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
!mick fix 5/1/2015 - added DO_FO_CALC initialization
      VLIDORT_ModIn%MBool%TS_DO_FO_CALC             = .false.

      VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = .false.
      VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = .true.  !  V2OS Full : Must be set now

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

      VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION       = .false.  ! Will be set later in the Full tests

!  Modified Control
!  ----------------

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 0 ! Set later

!  Modified Sunrays
!  ----------------

      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES       = 1
      VLIDORT_ModIn%MSunrays%TS_SZANGLES         = zero ! set later

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

      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS         = zero ! V2OS_Full : Now assiged from input file. See later
      VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS      = zero ! V2OS_Full : Now Calculated Internally

!  Modified Optical
!  ----------------

      VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT    = zero ! Set later

!  Linearized
!  ----------

      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG   = .false.  ! Set later
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER = 0        ! Set later
!mick fix 5/1/2015 - added these initializations
      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS  = 0       ! Set later
      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = 0       ! Set later
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS      = 0       ! Set later
      VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS       = 0       ! Set later
      VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES     = ' '
      VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES    = ' '

      VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input    = zero ! Set later
      VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input = zero
      VLIDORT_LinFixIn%Optical%TS_L_omega_total_input    = zero ! Set later

!mick fix 5/1/2015 - added these initializations
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

!mick fix 5/1/2015 - added SUP initializations

      VLIDORT_Sup%SS%TS_STOKES_SS  = ZERO
      VLIDORT_Sup%SS%TS_STOKES_DB  = ZERO
      VLIDORT_Sup%SS%TS_FO_ZMATRIX = ZERO
      VLIDORT_Sup%SS%TS_FO_R1SAVED = ZERO

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

!mick fix 5/1/2015 - added linearized SUP initializations
      VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_SS = ZERO
      VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_DB = ZERO
      VLIDORT_LinSup%SS%Col%TS_COLUMNWF_SS   = ZERO
      VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB   = ZERO
      VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB = ZERO

!  Now open the files
!  ------------------

!  Fort 88 (Chapman Factors) not required now.
!  Fort 87, now contains the height grid and eradius

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
!      do n = 1, nlay
!         do k = 1, n
!            read(88,*)nd,kd,Sun_Chapman2(n,k)
!         enddo
!      enddo
!      close(88)

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

!  Geometry

      VLIDORT_ModIn%MSunrays%TS_SZANGLES             = acos(emu0)/deg_to_rad
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT  = acos(emu) /deg_to_rad
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS         = phi

      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1) = acos(emu0)/deg_to_rad
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,2) = acos(emu) /deg_to_rad
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,3) = phi

!  Chapman (heights and eradius)
!    Factors are no longer copied from the file-read, Now Internal calculation!
!        VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS(:,:,1)   = Sun_Chapman2(:,:)

      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS     = eradius
      do n = 0, nlay
        VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(n)  = heights(n)
      enddo

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

!  Optical

      do n = 1, nlay
        n1 = nlay + 1 - n
        VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
        VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,1)  =   coefs(0:nmoms,n,1)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,2)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,5)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,6)  =   coefs(0:nmoms,n,2)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,11) =   coefs(0:nmoms,n,3)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,12) =   coefs(0:nmoms,n,6)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,15) = - coefs(0:nmoms,n,6)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,16) =   coefs(0:nmoms,n,4)
      enddo
      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = spars_save

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

!  Set status flag

      OPENFILEFLAG = .false.

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

      !write(*,*)
      write(*,*)'Done Vector VLIDORT, Baseline LPS calculation'

!  Exception handling, write-up (optional)

      CALL VLIDORT_WRITE_STATUS ( &
          'V2OS_full_LPS_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save baseline results

      do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
         do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_V_BAS(v,o1) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
            do q = 1, n_totalprofile_wfs
               profilewf_V_BAS(q,1:nlayers,v,o1) = &
                 VLIDORT_LinOut%Prof%TS_profilewf(q,1:nlayers,1,v,o1,upidx)
            enddo
            do q = 1, n_surface_wfs
               surfacewf_V_BAS(q,v,o1) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,upidx)
            enddo
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

      !write(*,*)
      write(*,*)'Done Scalar VLIDORT, Baseline LPS calculation'

!  Exception handling, write-up (optional)

      CALL VLIDORT_WRITE_STATUS ( &
          'V2OS_full_LPS_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save results

      do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
         do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_S_BAS(v,o1) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
            do q = 1, n_totalprofile_wfs
               profilewf_S_BAS(q,1:nlayers,v,o1) = &
                 VLIDORT_LinOut%Prof%TS_profilewf(q,1:nlayers,1,v,o1,upidx)
            enddo
            do q = 1, n_surface_wfs
              surfacewf_S_BAS(q,v,o1) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,upidx)
            enddo
         enddo
      enddo

!  Perform 2OS correction
!  ----------------------

!  Set conditions

      VLIDORT_FixIn%Cont%TS_nstokes = nstokes_saved

!  Call correction routine

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
      write(*,*)'Done 2OS correction, Baseline LPS calculation'

!  Exception handling, write-up (optional)

!mick fix 5/1/2015 - created new status subroutine for 2OS code 
      CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
          'V2OS_full_LPS_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save results

      do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
         do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_C_BAS(v,o1) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
            if (o1 == 1) Icorr(v) = VLIDORT_Out%Main%TS_2OSCORR_ICORR(v)
            do q = 1, n_totalprofile_wfs
               profilewf_C_BAS(q,1:nlayers,v,o1) = &
                 VLIDORT_LinOut%Prof%TS_profilewf(q,1:nlayers,1,v,o1,upidx)
            enddo
            do q = 1, n_surface_wfs
              surfacewf_C_BAS(q,v,o1) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,upidx)
            enddo
         enddo
      enddo

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
!      do task = 0,0

!  task counter

         t = task

!  Task 0 - surface jacobian

         if ( task .eq. 0 ) then
            if (debug) write(*,*)'Doing task 0: Finite Difference, perturb Lambertian albedo'
            VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =   spars_save * epsfac
            do n = 1, nlay
               n1 = nlay + 1 - n
               VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
               VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
            enddo
         endif

!  Tasks 1 to NLAYERS - profile jacobians

         if ( task.gt.0  ) then
            if (debug) write(*,'(a,i2,a,i2)')' Doing task',task, &
              ': Finite Difference, perturb molecular absorption Layer ',t
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
            VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0,n,1) = one
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

        !write(*,*)
        write(*,*)'Done Vector VLIDORT for FD task #', t

!  Exception handling, write-up (optional)

        CALL VLIDORT_WRITE_STATUS ( &
          'V2OS_full_LPS_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save results

        do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
          do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_V_PT(v,o1,t) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
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

        !write(*,*)
        write(*,*)'Done Scalar VLIDORT for FD task #', t

!  Exception handling, write-up (optional)

        CALL VLIDORT_WRITE_STATUS ( &
          'V2OS_full_LPS_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save results

        do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
          do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_S_PT(v,o1,t) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
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

        !write(*,*)
        write(*,*)'Done 2OS correction for FD task #', t

!  Exception handling, write-up (optional)

!mick fix 5/1/2015 - created new status subroutine for 2OS code 
         CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
             'V2OS_full_LPS_Execution.log', &
             35, OPENFILEFLAG, &
             VLIDORT_Out%Status )

!  Save results

        do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
          do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_C_PT(v,o1,t) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
          enddo
        enddo

!  End task loop

      enddo

!  write results

      open(36,file='V2OS_full_LPS_results.all',status='replace')
      v = 1
      write(36,77)'Baseline Stokes I (V,S,C) = ',stokes_V_BAS(v,1),stokes_S_BAS(v,1),stokes_C_BAS(v,1)
      write(36,77)'Baseline Stokes Q (V,S,C) = ',stokes_V_BAS(v,2),stokes_S_BAS(v,2),stokes_C_BAS(v,2)
      write(36,77)'Baseline Stokes U (V,S,C) = ',stokes_V_BAS(v,3),stokes_S_BAS(v,3),stokes_C_BAS(v,3)
      write(36,78)'Baseline Icorrection      = ',Icorr(v)

!mick fix 5/1/2015 - define q!
      q = 1
      do t = 0, nlay
         if ((t.eq.0) .or. (t.eq.1)) write(36,*)
         if ( t.eq.0 ) then
             do o1 = 1, 3
               write(36,56)'Surface WF (V,S,C), Stokes # = ',O1, &
                    spars_save*surfacewf_V_BAS(q,v,o1),(STOKES_V_PT(v,o1,t)-stokes_V_BAS(v,o1))/eps, &
                    spars_save*surfacewf_S_BAS(q,v,o1),(STOKES_S_PT(v,o1,t)-stokes_S_BAS(v,o1))/eps, &
                    spars_save*surfacewf_C_BAS(q,v,o1),(STOKES_C_PT(v,o1,t)-stokes_C_BAS(v,o1))/eps
             enddo
         else
             write(36,55)'Layer = ',t
             do o1 = 1, 3
               write(36,56)'Profile WF (V,S,C), Stokes # = ',O1, &
                    profilewf_V_BAS(q,t,v,o1),(STOKES_V_PT(v,o1,t)-stokes_V_BAS(v,o1))/eps, &
                    profilewf_S_BAS(q,t,v,o1),(STOKES_S_PT(v,o1,t)-stokes_S_BAS(v,o1))/eps, &
                    profilewf_C_BAS(q,t,v,o1),(STOKES_C_PT(v,o1,t)-stokes_C_BAS(v,o1))/eps
             enddo
         endif
      enddo
      close(36)

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

!  Re-set the two basic inputs for PART II (VERY IMPORTANT !!!!!!!!!!!!!!!!!!! )

      VLIDORT_FixIn%Cont%TS_NSTOKES            = nstokes
      VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .false.

!  Optical

      do n = 1, nlay
        n1 = nlay + 1 - n
        VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
        VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,1)  =   coefs(0:nmoms,n,1)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,2)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,5)  = - coefs(0:nmoms,n,5)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,6)  =   coefs(0:nmoms,n,2)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,11) =   coefs(0:nmoms,n,3)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,12) =   coefs(0:nmoms,n,6)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,15) = - coefs(0:nmoms,n,6)
        VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:nmoms,n1,16) =   coefs(0:nmoms,n,4)
      enddo
      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = spars_save

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

!  Set status flag

      OPENFILEFLAG = .false.

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

      !write(*,*)
      write(*,*)'Done Vector VLIDORT, Baseline LCS calculation'

!  Exception handling, write-up (optional)

      CALL VLIDORT_WRITE_STATUS ( &
          'V2OS_full_LCS_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save baseline results

      do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
         do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_V_BAS(v,o1) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
            do q = 1, n_totalcolumn_wfs
               columnwf_V_BAS(q,v,o1) = VLIDORT_LinOut%Col%TS_columnwf(q,1,v,o1,upidx)
            enddo
            do q = 1, n_surface_wfs
               surfacewf_V_BAS(q,v,o1) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,upidx)
            enddo
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

      !write(*,*)
      write(*,*)'Done Scalar VLIDORT, Baseline LCS calculation'

!  Exception handling, write-up (optional)

      CALL VLIDORT_WRITE_STATUS ( &
          'V2OS_full_LCS_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save results

      do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
         do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_S_BAS(v,o1) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
            do q = 1, n_totalcolumn_wfs
               columnwf_S_BAS(q,v,o1) = VLIDORT_LinOut%Col%TS_columnwf(q,1,v,o1,upidx)
            enddo
            do q = 1, n_surface_wfs
              surfacewf_S_BAS(q,v,o1) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,upidx)
            enddo
         enddo
      enddo

!  Perform 2OS correction
!  ----------------------

!  Set conditions

      VLIDORT_FixIn%Cont%TS_nstokes = nstokes_saved

!  Call correction routine

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
      write(*,*)'Done 2OS correction, Baseline LCS calculation'

!  Exception handling, write-up (optional)

!mick fix 5/1/2015 - created new status subroutine for 2OS code 
      CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
          'V2OS_full_LCS_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save results

      do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
         do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_C_BAS(v,o1) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
            do q = 1, n_totalcolumn_wfs
               columnwf_C_BAS(q,v,o1) = VLIDORT_LinOut%Col%TS_columnwf(q,1,v,o1,upidx)
            enddo
            do q = 1, n_surface_wfs
              surfacewf_C_BAS(q,v,o1) = VLIDORT_LinOut%Surf%TS_surfacewf(q,1,v,o1,upidx)
            enddo
         enddo
      enddo

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
!      do task = 1, 1

!  task counter

         t = task

!  Task 0 - surface jacobian

         if ( task .eq. 0 ) then
            if (debug) write(*,*)'Doing task 0: Finite Difference, perturb Lambertian albedo'
            VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO =   spars_save * epsfac
            do n = 1, nlay
               n1 = nlay + 1 - n
               VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n1)  = opd_save(n)
               VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n1) = ssa_save(n)
            enddo
         endif

!  Task 1 - column jacobian

         if ( task.gt.0  ) then
            if (debug) write(*,'(a,i2,a)')' Doing task',task, &
              ': Finite Difference, perturb total molecular absorption '
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
            VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0,n,1) = one
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

         !write(*,*)
         write(*,*)'Done Vector VLIDORT for FD task #', t

!  Exception handling, write-up (optional)

         CALL VLIDORT_WRITE_STATUS ( &
           'V2OS_full_LCS_Execution.log', &
           35, OPENFILEFLAG, &
           VLIDORT_Out%Status )

!  Save results

         do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
           do v = 1, VLIDORT_Out%Main%TS_n_geometries
             stokes_V_PT(v,o1,t) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
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

         !write(*,*)
         write(*,*)'Done Scalar VLIDORT for FD task #', t

!  Exception handling, write-up (optional)

         CALL VLIDORT_WRITE_STATUS ( &
           'V2OS_full_LCS_Execution.log', &
           35, OPENFILEFLAG, &
           VLIDORT_Out%Status )

!  Save results

         do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
           do v = 1, VLIDORT_Out%Main%TS_n_geometries
             stokes_S_PT(v,o1,t) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
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

        !write(*,*)
        write(*,*)'Done 2OS correction for FD task #', t

!  Exception handling, write-up (optional)

!mick fix 5/1/2015 - created new status subroutine for 2OS code 
         CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
             'V2OS_full_LCS_Execution.log', &
             35, OPENFILEFLAG, &
             VLIDORT_Out%Status )

!  Save results

         do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
           do v = 1, VLIDORT_Out%Main%TS_n_geometries
             stokes_C_PT(v,o1,t) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
           enddo
         enddo

!  End task loop

      enddo

!  Write results

      open(36,file='V2OS_full_LCS_results.all',status='replace')
      v = 1
      write(36,77)'Baseline Stokes I (V,S,C) = ',stokes_V_BAS(v,1),stokes_S_BAS(v,1),stokes_C_BAS(v,1)
      write(36,77)'Baseline Stokes Q (V,S,C) = ',stokes_V_BAS(v,2),stokes_S_BAS(v,2),stokes_C_BAS(v,2)
      write(36,77)'Baseline Stokes U (V,S,C) = ',stokes_V_BAS(v,3),stokes_S_BAS(v,3),stokes_C_BAS(v,3)
      write(36,78)'Baseline Icorrection      = ',Icorr(v)

!mick fix 5/1/2015 - changed hardcoded jac value of 1 in arrays below to use q now; defined q here
      q = 1
      do t = 0, 1
         write(36,*)
         if ( t.eq.0 ) then
             do o1 = 1, 3
               write(36,56)'Surface WF (V,S,C), Stokes # = ',O1, &
                    spars_save*surfacewf_V_BAS(q,v,o1),(STOKES_V_PT(v,o1,t)-stokes_V_BAS(v,o1))/eps, &
                    spars_save*surfacewf_S_BAS(q,v,o1),(STOKES_S_PT(v,o1,t)-stokes_S_BAS(v,o1))/eps, &
                    spars_save*surfacewf_C_BAS(q,v,o1),(STOKES_C_PT(v,o1,t)-stokes_C_BAS(v,o1))/eps
             enddo
         else
             do o1 = 1, 3
               write(36,56)'Column  WF (V,S,C), Stokes # = ',O1, &
                    columnwf_V_BAS(q,v,o1),(STOKES_V_PT(v,o1,t)-stokes_V_BAS(v,o1))/eps, &
                    columnwf_S_BAS(q,v,o1),(STOKES_S_PT(v,o1,t)-stokes_S_BAS(v,o1))/eps, &
                    columnwf_C_BAS(q,v,o1),(STOKES_C_PT(v,o1,t)-stokes_C_BAS(v,o1))/eps
             enddo
         endif
      enddo
      close(36)

!  Format statements

55    format(A,i2)
56    format(A,i3, 4x, 3(1p2e17.8,2x))
77    format(A,1p3e17.8)
78    format(A,1pe17.8)

!  Finish

      write(*,*)
      if (VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) THEN
         stop 'Successful run, with Warning: look in V2p7_2OS_Execution.log'
      else
         stop 'Successful run'
      endif

      end program V2OS_Full_LPCS_Tester
