      program V2OS_Full_Tester

!  This is a Full tester for the combined VLIDORT/2OS code
!   [ Uses the same setups and VLIDORT Inputs as the TwoOS_Self testers !!! ]

!  Module files for VLIDORT

      USE VLIDORT_PARS
      USE VLIDORT_IO_DEFS

      USE VLIDORT_AUX
      USE VLIDORT_INPUTS
      USE VLIDORT_MASTERS

      USE vlidort_2OScorr_master_m
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

!  Local Variables
!  ===============

!  Dump variables

      integer :: nlay,nmug,nmoms,nstokes,nfoumax,ncoefs(maxlayers)
      logical :: regular_ps,enhanced_ps
      double precision :: epsilon,phi,emu0,emu,spars,xmu(maxstreams),w(maxstreams),&
                          Sun_Chapman2(maxlayers,maxlayers),omega
      double precision :: ssa(maxlayers),opd(maxlayers),&
                          coefs(0:maxmoments_input,maxlayers,6),&
                          eradius,heights(0:maxlayers)

!  Flag for opening error output file

      LOGICAL :: OPENFILEFLAG

!  Local

      integer :: n, n1, k, l, nd, kd, ld, nstokes_saved, o1, v

!  Saved results

      DOUBLE PRECISION :: STOKES_V ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: STOKES_S ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: STOKES_C ( MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: Icorr ( MAX_GEOMETRIES )

!  Start of code
!  =============

!  Hardwired input

!  Fixed Boolean inputs
!  --------------------

      VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE      = .true.    ! V2OS Full - Must be TRUE for the Full test
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
      VLIDORT_FixIn%Cont%TS_NFINELAYERS      = 2 ! V2OS Full - Should be set for the Full test
      VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = 2 ! V2OS Full - Should be set for the Full test

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

!      VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION       = .true.
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
!mick fix 5/1/2015 - initialize all elements of USER_LEVELS
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

!  Now open the files
!  ------------------

!  Fort 88 (Chapman Factors) not required now.
!  Fort 87, now contains the height grid and eradius

      open(87,file='fort.87_Rg_h',status='old')
      read(87,*)nlay,nmug,nstokes,nmoms,nfoumax,regular_ps,enhanced_ps
      read(87,*)epsilon,phi,emu0,emu,spars,eradius,heights(0)
      do k = 1, nmug
         read(87,*)kd,xmu(k),w(k)  ! Not used
      enddo
      do n = 1, nlay
         read(87,*)heights(n),ssa(n),opd(n),ncoefs(n)
         do l = 0, 2*nmug-1
            read(87,*)nd,ld,(coefs(l,n,k),k=1,6)
         enddo
      enddo
      close(87)

!      Sun_Chapman2 = 0.0d0
!      open(88,file='fort.88',status='old')
!         do n = 1, nlay
!            do k = 1, n
!              read(88,*)nd,kd,Sun_Chapman2(n,k)
!            enddo
!       enddo
!       close(88)

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

!  Precautions for VLIDORT

      do n = 1, nlay
         VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0,n,1) = one
         omega = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)
         if ( omega.gt.0.99999999d0 ) VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n) = 0.99999999d0
      enddo

!  Initialize output

      STOKES_V = ZERO
      STOKES_S = ZERO
      STOKES_C = ZERO
      Icorr = ZERO

!  Set status flag

      OPENFILEFLAG = .false.

!  First Call. VECTOR VLIDORT
!  ==========================

!  Set conditions. Full Vector, No 2OS

      nstokes_saved = VLIDORT_FixIn%Cont%TS_nstokes 
      VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .false.

!  VLIDORT call

      CALL VLIDORT_MASTER ( &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

      write(*,*)
      write(*,*)'Done Vector VLIDORT'

!  Exception handling, write-up (optional)

      CALL VLIDORT_WRITE_STATUS ( &
          'V2OS_full_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save results

      do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
         do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_V(v,o1) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
         enddo
         write(*,*)'Vector RT SS,SS+MS',&
                    VLIDORT_Sup%SS%TS_stokes_ss(1,1,o1,upidx),&
                    VLIDORT_Out%Main%TS_stokes(1,1,o1,upidx)
      enddo

!  Second Call. SCALAR VLIDORT + 2OS CORRECTION
!  ============================================

!  Set conditions. Scalar Call, 2OS correction on.

      VLIDORT_FixIn%Cont%TS_nstokes            = 1
      VLIDORT_ModIn%MBool%TS_DO_2OS_CORRECTION = .true.

!  VLIDORT call

      CALL VLIDORT_MASTER ( &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

      write(*,*)
      write(*,*)'Done Scalar VLIDORT'

!  Exception handling, write-up (optional)

      CALL VLIDORT_WRITE_STATUS ( &
          'V2OS_full_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save results

      do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
         do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_S(v,o1) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
         enddo
         write(*,*)'Scalar RT SS,SS+MS',&
           VLIDORT_Sup%SS%TS_stokes_ss(1,1,o1,upidx),&
           VLIDORT_Out%Main%TS_stokes(1,1,o1,upidx)
      enddo

!  Perform 2OS correction
!  ----------------------

!  Set conditions

      VLIDORT_FixIn%Cont%TS_nstokes = nstokes_saved

!  Call 2OS correction routine

      CALL VLIDORT_2OSCORR_MASTER ( &
           VLIDORT_FixIn, &
           VLIDORT_ModIn, &
           VLIDORT_Sup,   &
           VLIDORT_Out )

      write(*,*)
      write(*,*)'Done 2OS correction'

!  Exception handling, write-up (optional)

!mick fix 5/1/2015 - created new status subroutine for 2OS code 
      CALL VLIDORT_2OSCORR_WRITE_STATUS ( &
          'V2OS_full_Execution.log', &
          35, OPENFILEFLAG, &
          VLIDORT_Out%Status )

!  Save results

      do o1 = 1, VLIDORT_FixIn%Cont%TS_nstokes
         do v = 1, VLIDORT_Out%Main%TS_n_geometries
            stokes_C(v,o1) = VLIDORT_Out%Main%TS_stokes(1,v,o1,upidx)
            if (o1 == 1) Icorr(v) = VLIDORT_Out%Main%TS_2OSCORR_ICORR(v)
         enddo
         write(*,*)'Corrected RT SS,SS+MS',&
                    VLIDORT_Sup%SS%TS_stokes_ss(1,1,o1,upidx),&
                    VLIDORT_Out%Main%TS_stokes(1,1,o1,upidx)
      enddo

!  Debug output

      write(*,*)
      write(*,*)'R2(1) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,1)
      write(*,*)'R2(2) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,2)
      write(*,*)'R2(3) = ',VLIDORT_Out%Main%TS_2OSCORR_R2(1,3)
      write(*,*)'Icorr = ',VLIDORT_Out%Main%TS_2OSCORR_ICORR(1)

!  Write results

      open(36,file='V2OS_full_results.all',status='replace')
      v = 1
      write(36,77)'Baseline Stokes I (V,S,C) = ',stokes_V(v,1),stokes_S(v,1),stokes_C(v,1)
      write(36,77)'Baseline Stokes Q (V,S,C) = ',stokes_V(v,2),stokes_S(v,2),stokes_C(v,2)
      write(36,77)'Baseline Stokes U (V,S,C) = ',stokes_V(v,3),stokes_S(v,3),stokes_C(v,3)
      write(36,78)'Baseline Icorrection      = ',Icorr(v)
      close(36)
77    format(A,1p3e17.8)
78    format(A,1pe17.8)

!  Finish

      write(*,*)
      if (VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .eq. VLIDORT_WARNING ) THEN
         stop 'Successful run, with Warning: look in V2OS_full_Execution.log'
      else
         stop 'Successful run'
      endif

      end program V2OS_Full_Tester
