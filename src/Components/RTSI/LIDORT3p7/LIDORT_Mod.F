      module LIDORT_Mod

      USE LIDORT_PARS
      
      USE LIDORT_IO_DEFS
      USE BRDF_SUP_MOD
      
      implicit NONE

      PUBLIC  LIDORT_Init ! Initialize once for all pixels
      PUBLIC  LIDORT_Fin  ! Finalize once
       
      !  LIDORT input/output structures

      TYPE LIDORT_IO
          !  LIDORT input structures
          TYPE(LIDORT_Fixed_Inputs)             :: LIDORT_FixIn
          TYPE(LIDORT_Modified_Inputs)          :: LIDORT_ModIn



          !  LIDORT output structures
!!          TYPE(Outputs_Main_def)             :: LIDORT_Outputs
          TYPE(LIDORT_Outputs)                  :: LIDORT_Out

!!          TYPE(Exception_Handling_def)       :: LIDORT_Status
          TYPE(BRDF_Output_Exception_Handling)   :: BRDF_Sup_OutputStatus


          !  BRDF input and output structure
!!          TYPE(BRDF_Surface_def)             :: LIDORT_BRDF_inputs
!!          TYPE(BRDF_Sup_Inputs)                 :: BRDF_Sup_In
          TYPE(BRDF_Sup_Outputs)                :: BRDF_Sup_Out

          !  LIDORT supplements i/o structure

          TYPE(LIDORT_Sup_InOut)                :: LIDORT_Sup

      END TYPE LIDORT_IO
      

      TYPE LIDORT
        logical     :: initialized = .false.  
        integer     :: NSTREAMS = 6        ! Number of half-space streams
        integer     :: NBEAMS = 1          ! Number of solar zenith angles
        integer     :: N_USER_STREAMS = 1  ! Number of Viewing zenith angles
        integer     :: N_USER_RELAZMS = 1  ! Number of relative azimuth angles
        integer     :: N_USER_LEVELS  = 1  ! Number of user-defined vertical output levels
        integer     :: N_USER_OBSGEOMS = 1 ! Number of azimuth angles calculated by surface supplement

        type(LIDORT_IO) :: VIO

      END TYPE LIDORT

      contains

!.............................................................................

      subroutine LIDORT_Init (self, km, rc)
                                   
      USE LIDORT_PARS
      USE LIDORT_Inputs_def

      USE LIDORT_INPUTS
     
      type(LIDORT),              intent(inout) :: self

      integer,                    intent(in)    :: km     ! number of atmospheric layers
      integer,                    intent(out)   :: rc     ! error code
      
!                     ---

      logical :: vector


!  Local variables

      logical                                :: DO_FULLRAD_MODE
      logical                                :: DO_SSCORR_NADIR
      logical                                :: DO_SSCORR_OUTGOING
      logical                                :: DO_SSCORR_TRUNCATION
      logical                                :: DO_SSFULL
      logical                                :: DO_SS_EXTERNAL
      logical                                :: DO_DOUBLE_CONVTEST
      logical                                :: DO_SOLAR_SOURCES
      logical                                :: DO_PLANE_PARALLEL
      logical                                :: DO_FO_CALC
      logical                                :: DO_OBSERVATION_GEOMETRY
      logical                                :: DO_REFRACTIVE_GEOMETRY
      logical                                :: DO_CHAPMAN_FUNCTION
      logical                                :: DO_RAYLEIGH_ONLY
      logical                                :: DO_DELTAM_SCALING
      logical                                :: DO_SOLUTION_SAVING
      logical                                :: DO_BVP_TELESCOPING
      logical                                :: DO_UPWELLING
      logical                                :: DO_DNWELLING
      logical                                :: DO_QUAD_OUTPUT
      logical                                :: DO_USER_VZANGLES
      logical                                :: DO_ADDITIONAL_MVOUT
      logical                                :: DO_MVOUT_ONLY
      logical                                :: DO_DEBUG_WRITE
      logical                                :: DO_WRITE_INPUT
      logical                                :: DO_WRITE_SCENARIO
      logical                                :: DO_WRITE_FOURIER
      logical                                :: DO_WRITE_RESULTS
      character (LEN=60)                     :: INPUT_WRITE_FILENAME
      character (LEN=60)                     :: SCENARIO_WRITE_FILENAME
      character (LEN=60)                     :: FOURIER_WRITE_FILENAME
      character (LEN=60)                     :: RESULTS_WRITE_FILENAME
      integer                                :: TAYLOR_ORDER
      integer                                :: NSTOKES
      integer                                :: NSTREAMS
      integer                                :: NLAYERS
      integer                                :: NFINELAYERS
      integer                                :: NGREEK_MOMENTS_INPUT
      real*8                                 :: LIDORT_ACCURACY
      real*8                                 :: FLUX_FACTOR
      integer                                :: N_SZANGLES
      real*8, dimension( MAX_SZANGLES )      :: SZANGLES 
      real*8                                 :: EARTH_RADIUS
      real*8                                 :: RFINDEX_PARAMETER
      real*8                                 :: GEOMETRY_SPECHEIGHT
      integer                                :: N_USER_OBSGEOMS
      real*8, dimension( MAX_USER_OBSGEOMS,3):: USER_OBSGEOMS
      integer                                :: N_USER_RELAZMS
      real*8, dimension( MAX_USER_RELAZMS )  :: USER_RELAZMS  
      integer                                :: N_USER_VZANGLES
      real*8, dimension( MAX_USER_VZANGLES ) :: USER_VZANGLES 
      integer                                :: N_USER_LEVELS
      real*8, dimension( MAX_USER_LEVELS )   :: USER_LEVELS 
      logical                                :: DO_LAMBERTIAN_SURFACE
      real*8                                 :: LAMBERTIAN_ALBEDO
      logical                                :: DO_THERMAL_EMISSION
      integer                                :: N_THERMAL_COEFFS
      real*8, dimension( 0:MAXLAYERS )       :: THERMAL_BB_INPUT 
      logical                                :: DO_SURFACE_EMISSION
      logical                                :: DO_SURFACE_LEAVING
      logical                                :: DO_SL_ISOTROPIC
      real*8                                 :: SURFBB
      logical                                :: DO_THERMAL_TRANSONLY
      logical                                :: DO_SPECIALIST_OPTION_1
      logical                                :: DO_SPECIALIST_OPTION_2
      logical                                :: DO_SPECIALIST_OPTION_3
      logical                                :: DO_TOA_CONTRIBS

      rc = 0

     
!    initialize variables to false or zero
!    -------------------------------------    

      call LIDORT_INIT_CONTROL_VARS (DO_FULLRAD_MODE, DO_SSCORR_NADIR,
     &     DO_SSCORR_OUTGOING, DO_FO_CALC, DO_SSCORR_TRUNCATION,
     &     DO_SS_EXTERNAL, DO_SSFULL,
     &     DO_DOUBLE_CONVTEST, DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, 
     &     DO_CHAPMAN_FUNCTION, 
     &     DO_RAYLEIGH_ONLY, DO_DELTAM_SCALING, DO_SOLUTION_SAVING, 
     &     DO_BVP_TELESCOPING, DO_UPWELLING, DO_DNWELLING, DO_QUAD_OUTPUT, 
     &     DO_USER_VZANGLES, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, 
     &     DO_DEBUG_WRITE, DO_WRITE_INPUT, DO_WRITE_SCENARIO, 
     &     DO_WRITE_FOURIER, DO_WRITE_RESULTS, DO_LAMBERTIAN_SURFACE, 
     &     DO_OBSERVATION_GEOMETRY,DO_SPECIALIST_OPTION_1, 
     &     DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3, DO_TOA_CONTRIBS,
     &     DO_SURFACE_LEAVING, DO_SL_ISOTROPIC)


      call LIDORT_INIT_THERMAL_VARS(DO_THERMAL_EMISSION, N_THERMAL_COEFFS, 
     &     THERMAL_BB_INPUT, DO_SURFACE_EMISSION, 
     &     SURFBB, DO_THERMAL_TRANSONLY)

      call LIDORT_INIT_MODEL_VARS(TAYLOR_ORDER,NSTOKES, 
     &     NSTREAMS, NLAYERS, NFINELAYERS, 
     &     NGREEK_MOMENTS_INPUT, LIDORT_ACCURACY,  
     &     FLUX_FACTOR, N_SZANGLES, SZANGLES,       
     &     EARTH_RADIUS, RFINDEX_PARAMETER,       
     &     GEOMETRY_SPECHEIGHT, N_USER_RELAZMS,     
     &     USER_RELAZMS, N_USER_VZANGLES,           
     &     USER_VZANGLES, N_USER_LEVELS, USER_LEVELS,
     &     N_USER_OBSGEOMS, USER_OBSGEOMS,
     &     LAMBERTIAN_ALBEDO)


      !    Customize those parameters that do not vary from pixel to pixel
!     ---------------------------------------------------------------


!                         Modes of Operation
!                         ------------------

      DO_FULLRAD_MODE    = .true.  ! Do full Stokes vector calculation?
      DO_SSCORR_NADIR    = .false. ! Do nadir single scatter correction?
      DO_SSCORR_OUTGOING = .true.  ! Do outgoing single scatter correction?
      DO_FO_CALC         = .false. ! Do outgoung single Nakajima scatter correction?
      DO_SSFULL          = .false. ! Do Full-up single scatter calculation?
!      DO_DBCORRECTION    = .true.  ! Do direct beam correction?
      DO_DOUBLE_CONVTEST = .true.  ! Perform double convergence test?
       
!                            Solar Sources
!                            -------------

      DO_SOLAR_SOURCES    = .true.     ! Include solar sources?
      DO_PLANE_PARALLEL   = .false.    ! Plane-parallel treatment of direct beam?
      DO_CHAPMAN_FUNCTION = .true.     ! Perform internal Chapman function calculation?
      DO_REFRACTIVE_GEOMETRY = .false. ! Beam path with refractive atmosphere?
     
!                         Performance Control
!                         -------------------

      DO_RAYLEIGH_ONLY     = .false. ! Rayleigh atmosphere only?
      DO_DELTAM_SCALING    = .true.  ! Include Delta-M scaling?
      DO_SSCORR_TRUNCATION = .true. ! Additional Delta-M scaling for SS correction?
      DO_SOLUTION_SAVING   = .false. ! Solution saving mode?
      DO_BVP_TELESCOPING   = .false. ! Boundary value problem telescoping mode?
      
!                      User-defined output control
!                      ---------------------------

      DO_UPWELLING = .true.     ! Upwelling output?
      DO_DNWELLING = .false.    ! Downwelling output?
      DO_USER_VZANGLES = .true. ! User-defined Viewing zenith angles?
      DO_OBSERVATION_GEOMETRY = .false. ! Do Observation Geometry?
      

      DO_ADDITIONAL_MVOUT = .true.  ! Generate mean value output additionally?
      DO_MVOUT_ONLY       = .false. ! Generate only mean value output?
      
!                           Write Control
!                           -------------

      DO_DEBUG_WRITE    = .false. ! Debug write?
      DO_WRITE_INPUT    = .false. ! Input control write?
      DO_WRITE_SCENARIO = .false. ! Input scenario write?
      DO_WRITE_FOURIER  = .false. ! Fourier component output write?
      DO_WRITE_RESULTS  = .false. ! Results write?
      
      INPUT_WRITE_FILENAME    = '/dev/null' ! filename for input write
      SCENARIO_WRITE_FILENAME = '/dev/null' ! filename for scenario write
      FOURIER_WRITE_FILENAME  = '/dev/null' ! Fourier output filename
      RESULTS_WRITE_FILENAME  = '/dev/null' ! filename for main output

      
      
! WAS NOT THERE in the preVIOus version (f77), do we need them??????

!                Specialist options. Should always be initialized here
!                -----------------------------------------------------
      DO_SPECIALIST_OPTION_1 = .false.    
      DO_SPECIALIST_OPTION_2 = .false.
      DO_SPECIALIST_OPTION_3 = .false.
      
!                  TOA contributions flag
!             --------------------------
      DO_TOA_CONTRIBS = .false.
      
!                     Quadrature output is a debug flag only.
!                    ----------------------------------------
      DO_QUAD_OUTPUT  = .false.
      

     
      NSTREAMS = self%NSTREAMS         ! Number of half-space streams
      NLAYERS = km                    ! Number of atmospheric layers
      NFINELAYERS = 3                 ! Number of fine layers (outgoing sphericity correction)
      NGREEK_MOMENTS_INPUT = 300     ! Number of scattering matrix expansion coefficients
      TAYLOR_ORDER = 3                ! Number of small-number terms in Taylor series expansions
      N_USER_OBSGEOMS = 1             ! Number of observation Geometry inputs
      

!                            accuracy input
!                            --------------

      LIDORT_ACCURACY = 0.0001  ! Fourier series convergence
      FLUX_FACTOR = 1.00d0      ! Solar flux constant; =1 if no solar sources.

!                        Pseudo-spherical inputs
!                        -----------------------

      EARTH_RADIUS = 6371.0 ! Earth radius (km)
      RFINDEX_PARAMETER = 0.000288 ! Refractive index parameter
      GEOMETRY_SPECHEIGHT = 0.0 ! Input geometry specification height [km]


!                          Thermal controls
!                          ---------------- 

      DO_THERMAL_EMISSION  = .false.  ! Do thermal emission?


      DO_THERMAL_TRANSONLY = .false.  ! Do thermal emission, transmittance only?
      N_THERMAL_COEFFS     = 2        ! Number of thermal coefficients
      DO_SURFACE_EMISSION  =  .false. ! Do Surface emission?



      !  normal execution: Copy all variables and return
!  -----------------------------------------------

!  Copy data to structure variables:

!  Fixed Boolean inputs

      self%VIO%LIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE        = DO_FULLRAD_MODE
      self%VIO%LIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION   = DO_SSCORR_TRUNCATION
      self%VIO%LIDORT_FixIn%Bool%TS_DO_SSFULL              = DO_SSFULL
      self%VIO%LIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION    = DO_THERMAL_EMISSION
      self%VIO%LIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION    = DO_SURFACE_EMISSION
      self%VIO%LIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL      = DO_PLANE_PARALLEL
      !self%VIO%LIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE       = DO_BRDF_SURFACE
      self%VIO%LIDORT_FixIn%Bool%TS_DO_UPWELLING           = DO_UPWELLING
      self%VIO%LIDORT_FixIn%Bool%TS_DO_DNWELLING           = DO_DNWELLING
      self%VIO%LIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT         = DO_QUAD_OUTPUT
      self%VIO%LIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS        = DO_TOA_CONTRIBS
      !self%VIO%LIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = DO_LAMBERTIAN_SURFACE
      self%VIO%LIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_1 = DO_SPECIALIST_OPTION_1
      self%VIO%LIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_2 = DO_SPECIALIST_OPTION_2
      self%VIO%LIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_3 = DO_SPECIALIST_OPTION_3


      self%VIO%LIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING     = DO_SURFACE_LEAVING
      self%VIO%LIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC        = DO_SL_ISOTROPIC


!  Modified Boolean inputs

      self%VIO%LIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR        = DO_SSCORR_NADIR
      self%VIO%LIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING     = DO_SSCORR_OUTGOING
      self%VIO%LIDORT_ModIn%MBool%TS_DO_FO_CALC             = DO_FO_CALC
      self%VIO%LIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST     = DO_DOUBLE_CONVTEST
      self%VIO%LIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = DO_SOLAR_SOURCES
      self%VIO%LIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = DO_REFRACTIVE_GEOMETRY
      self%VIO%LIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = DO_CHAPMAN_FUNCTION
      self%VIO%LIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY       = DO_RAYLEIGH_ONLY
      !LIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY              = DO_ISOTROPIC_ONLY
      !LIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH                  = DO_NO_AZIMUTH
      !LIDORT_ModIn%MBool%TS_DO_ALL_FOURIER                 = DO_ALL_FOURIER
      self%VIO%LIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING      = DO_DELTAM_SCALING
      self%VIO%LIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING     = DO_SOLUTION_SAVING
      self%VIO%LIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING     = DO_BVP_TELESCOPING
      !LIDORT_ModIn%MBool%TS_DO_USER_STREAMS                = DO_USER_STREAMS
      self%VIO%LIDORT_ModIn%MBool%TS_DO_USER_VZANGLES       = DO_USER_VZANGLES
      self%VIO%LIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = DO_ADDITIONAL_MVOUT
      self%VIO%LIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY          = DO_MVOUT_ONLY
      self%VIO%LIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY   = DO_THERMAL_TRANSONLY

!  Fixed control inputs
!      self%VIO%LIDORT_FixIn%Cont%TS_NSTOKES               = NSTOKES
      self%VIO%LIDORT_FixIn%Cont%TS_TAYLOR_ORDER           = TAYLOR_ORDER
      self%VIO%LIDORT_FixIn%Cont%TS_NSTREAMS               = self%NSTREAMS
      self%VIO%LIDORT_FixIn%Cont%TS_NLAYERS                = NLAYERS
      self%VIO%LIDORT_FixIn%Cont%TS_NFINELAYERS            = NFINELAYERS
      self%VIO%LIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS       = N_THERMAL_COEFFS
      self%VIO%LIDORT_FixIn%Cont%TS_LIDORT_ACCURACY        = LIDORT_ACCURACY

!  Modified control inputs

      self%VIO%LIDORT_ModIn%Mcont%TS_NGREEK_MOMENTS_INPUT   = NGREEK_MOMENTS_INPUT

!  Beam inputs

      self%VIO%LIDORT_FixIn%SunRays%TS_FLUX_FACTOR               = FLUX_FACTOR      
      self%VIO%LIDORT_ModIn%MSunRays%TS_N_SZANGLES               = self%NBEAMS
      self%VIO%LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS           = self%N_USER_RELAZMS     
      self%VIO%LIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES          = self%N_USER_STREAMS
      self%VIO%LIDORT_FixIn%UserVal%TS_N_USER_LEVELS             = self%N_USER_LEVELS
      
      self%VIO%LIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT      = GEOMETRY_SPECHEIGHT

      self%VIO%LIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS          = N_USER_OBSGEOMS     
      self%VIO%LIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT      = USER_OBSGEOMS
      

!  Fixed Chapman function inputs

      !LIDORT_Chapman_inputs%TS_HEIGHT_GRID                     = HEIGHT_GRID
      !LIDORT_Chapman_inputs%TS_PRESSURE_GRID                   = PRESSURE_GRID
      !LIDORT_Chapman_inputs%TS_TEMPERATURE_GRID                = TEMPERATURE_GRID
      !LIDORT_Chapman_inputs%TS_FINEGRID                        = FINEGRID
      self%VIO%LIDORT_ModIn%MChapman%TS_EARTH_RADIUS            = EARTH_RADIUS
      self%VIO%LIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER       = RFINDEX_PARAMETER

!  Fixed optical inputs

      !LIDORT_Optical_inputs%TS_DELTAU_VERT_INPUT               = DELTAU_VERT_INPUT
      !LIDORT_Optical_inputs%TS_OMEGA_TOTAL_INPUT               = OMEGA_TOTAL_INPUT
      !LIDORT_Optical_inputs%TS_GREEKMAT_TOTAL_INPUT            = GREEKMAT_TOTAL_INPUT
      self%VIO%LIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT        = THERMAL_BB_INPUT
      self%VIO%LIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO       = LAMBERTIAN_ALBEDO
      self%VIO%LIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT        = SURFBB

!  Fixed write inputs

      self%VIO%LIDORT_FixIn%Write%TS_DO_DEBUG_WRITE            = DO_DEBUG_WRITE
  
      self%VIO%LIDORT_FixIn%Write%TS_DO_WRITE_INPUT            = DO_WRITE_INPUT
      self%VIO%LIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME      = INPUT_WRITE_FILENAME

      self%VIO%LIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO         = DO_WRITE_SCENARIO
      self%VIO%LIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME   = SCENARIO_WRITE_FILENAME

      self%VIO%LIDORT_FixIn%Write%TS_DO_WRITE_FOURIER          = DO_WRITE_FOURIER
      self%VIO%LIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME    = FOURIER_WRITE_FILENAME

      self%VIO%LIDORT_FixIn%Write%TS_DO_WRITE_RESULTS          = DO_WRITE_RESULTS
      self%VIO%LIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME    = RESULTS_WRITE_FILENAME

!                            Error Checking
!                            --------------

!      if ( NSTOKES  .GT. MAXSTOKES  )                   rc = 1
      if ( NSTREAMS .GT. MAXSTREAMS )                   rc = 2
      if ( NLAYERS  .GT. MAXLAYERS  )                   rc = 3
      if ( NFINELAYERS .GT. MAXFINELAYERS )             rc = 4
      if ( NGREEK_MOMENTS_INPUT .GT. MAXMOMENTS_INPUT)  rc = 5
      if ( N_SZANGLES .GT. MAX_SZANGLES )               rc = 6
      if ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS )       rc = 7
      if ( N_USER_VZANGLES .GT. MAX_USER_VZANGLES )     rc = 8
      if ( N_USER_LEVELS .GT. MAX_USER_LEVELS )         rc = 9 
      if ( N_THERMAL_COEFFS .GT. MAX_THERMAL_COEFFS )   rc =  10
     
      
!      print*, 'test...', NSTREAMS, MAXSTREAMS, rc

!     All done
!     --------
      self%initialized = .true.

      end subroutine LIDORT_Init

!.............................................................................

      subroutine LIDORT_Fin (self, rc)
      type(LIDORT),  intent(inout) :: self
      integer,        intent(out) :: rc
      self%initialized = .false.
      
      rc = 0
      
      end subroutine LIDORT_Fin

!.............................................................................

      end module LIDORT_Mod
