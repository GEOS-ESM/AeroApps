      module VLIDORT90_LambMod
 
      USE VLIDORT_PARS
      
      USE VLIDORT_BRDFSup_def
      USE VLIDORT_Inputs_def
      USE VLIDORT_Outputs_def

      USE VLIDORT_AUX
      USE VLIDORT_INPUTS
      USE VLIDORT_MASTERS


      implicit NONE

!     PRIVATE

      PUBLIC  VLIDORT_LambInit ! Initialize once
      PUBLIC  VLIDORT_LambRun  ! Run for each profile
      PUBLIC  VLIDORT_LambFin  ! Finalize once
      PUBLIC  VLIDORT_LER      ! Lambertian Equivalent Reflectivity
      PUBLIC  VLIDORT_AI       ! Aerosol Index
      
      
!  VLIDORT input structures

      TYPE VLIDORT_IO
 
          TYPE(Fixed_Inputs_Boolean_def)     :: VLIDORT_FBoolean_inputs
          TYPE(Fixed_Inputs_Control_def)     :: VLIDORT_FControl_inputs
          TYPE(Modified_Inputs_Boolean_def)  :: VLIDORT_MBoolean_inputs
          TYPE(Modified_Inputs_Control_def)  :: VLIDORT_MControl_inputs

          TYPE(Fixed_Inputs_Sunrays_def)     :: VLIDORT_Beam_inputs
          TYPE(Fixed_Inputs_UserVAL_def)     :: VLIDORT_User_inputs
          TYPE(Fixed_Inputs_Chapman_def)     :: VLIDORT_Chapman_inputs

          TYPE(Fixed_Inputs_Optical_def)     :: VLIDORT_Optical_inputs

          TYPE(Fixed_Inputs_Write_def)       :: VLIDORT_Write_inputs

!  VLIDORT output structures

          TYPE(Input_Exception_Handling_def) :: VLIDORT_InputStatus
          TYPE(Outputs_Main_def)             :: VLIDORT_Outputs
          TYPE(Exception_Handling_def)       :: VLIDORT_Status
!  BRDF input structure

          TYPE(BRDF_Surface_def)             :: VLIDORT_BRDF_inputs

      END TYPE VLIDORT_IO

      type VLIDORT_Lamb

         logical         :: initialized = .false.
         real*8          :: wavelength          ! in [nm]
         real*8          :: albedo              ! surface albedo
         real*8          :: solar_zenith        ! in degrees
         real*8          :: relative_azymuth    ! in degrees
         real*8          :: sensor_zenith       ! in degrees
         integer         :: nMom                ! number of momemts read (phase function) 
         real*8          :: MISSING             ! MISSING VALUE
         real*8, pointer :: tau(:)              ! aerosol tau
         real*8, pointer :: ssa(:)              ! aerosol ssa
         real*8, pointer ::   g(:)              ! aerosol asymmetry factor
         real*8, pointer ::  pe(:)              ! pressure    at layer edges [Pa]
         real*8, pointer ::  ze(:)              ! height      at layer edges [m]
         real*8, pointer ::  te(:)              ! temperature at layer edges [K]
         real*8, pointer:: pmom(:,:,:)          ! components of the scattering phase matrix

         type(VLIDORT_IO) :: VIO

       end type VLIDORT_Lamb

      Contains

!.............................................................................

      subroutine VLIDORT_LambInit (self, km, rc, scalar)
                                   
      USE VLIDORT_PARS
      USE VLIDORT_Inputs_def
      USE VLIDORT_Outputs_def

      USE VLIDORT_AUX
      USE VLIDORT_INPUTS
     
      type(VLIDORT_Lamb),         intent(inout) :: self

      integer,                    intent(in)    :: km     ! number of atmospheric layers
      integer,                    intent(out)   :: rc     ! error code
      

      logical, optional,    intent(in)  :: scalar ! if true, perform scalar calculation
                                                  !  default is VECTOR calculation

!                     ---

      logical :: vector


!  Local variables

      LOGICAL ::             DO_FULLRAD_MODE
      LOGICAL ::             DO_SSCORR_NADIR
      LOGICAL ::             DO_SSCORR_OUTGOING
      LOGICAL ::             DO_SSCORR_TRUNCATION
      LOGICAL ::             DO_SSFULL
      LOGICAL ::             DO_DOUBLE_CONVTEST
      LOGICAL ::             DO_SOLAR_SOURCES
      LOGICAL ::             DO_PLANE_PARALLEL
      LOGICAL ::             DO_REFRACTIVE_GEOMETRY
      LOGICAL ::             DO_CHAPMAN_FUNCTION
      LOGICAL ::             DO_RAYLEIGH_ONLY
      LOGICAL ::             DO_DELTAM_SCALING
      LOGICAL ::             DO_SOLUTION_SAVING
      LOGICAL ::             DO_BVP_TELESCOPING
      LOGICAL ::             DO_UPWELLING
      LOGICAL ::             DO_DNWELLING
      LOGICAL ::             DO_QUAD_OUTPUT
      LOGICAL ::             DO_USER_VZANGLES
      LOGICAL ::             DO_ADDITIONAL_MVOUT
      LOGICAL ::             DO_MVOUT_ONLY
      LOGICAL ::             DO_DEBUG_WRITE
      LOGICAL ::             DO_WRITE_INPUT
      LOGICAL ::             DO_WRITE_SCENARIO
      LOGICAL ::             DO_WRITE_FOURIER
      LOGICAL ::             DO_WRITE_RESULTS
      CHARACTER (LEN=60) ::  INPUT_WRITE_FILENAME
      CHARACTER (LEN=60) ::  SCENARIO_WRITE_FILENAME
      CHARACTER (LEN=60) ::  FOURIER_WRITE_FILENAME
      CHARACTER (LEN=60) ::  RESULTS_WRITE_FILENAME
      INTEGER ::             NSTOKES
      INTEGER ::             NSTREAMS
      INTEGER ::             NLAYERS
      INTEGER ::             NFINELAYERS
      INTEGER ::             NGREEK_MOMENTS_INPUT
      DOUBLE PRECISION ::    VLIDORT_ACCURACY
      DOUBLE PRECISION ::    FLUX_FACTOR
      INTEGER ::             N_SZANGLES
      DOUBLE PRECISION ::    SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION ::    EARTH_RADIUS
      DOUBLE PRECISION ::    RFINDEX_PARAMETER
      DOUBLE PRECISION ::    GEOMETRY_SPECHEIGHT
      INTEGER ::             N_USER_RELAZMS
      DOUBLE PRECISION ::    USER_RELAZMS  ( MAX_USER_RELAZMS )
      INTEGER ::             N_USER_VZANGLES
      DOUBLE PRECISION ::    USER_VZANGLES ( MAX_USER_VZANGLES )
      INTEGER ::             N_USER_LEVELS
      DOUBLE PRECISION ::    USER_LEVELS ( MAX_USER_LEVELS )
      LOGICAL ::             DO_LAMBERTIAN_SURFACE
      DOUBLE PRECISION ::    LAMBERTIAN_ALBEDO
      LOGICAL ::             DO_THERMAL_EMISSION
      INTEGER ::             N_THERMAL_COEFFS
      DOUBLE PRECISION ::    THERMAL_BB_INPUT ( 0:MAXLAYERS )
      LOGICAL ::             DO_SURFACE_EMISSION
      DOUBLE PRECISION ::    SURFBB
      LOGICAL ::             DO_THERMAL_TRANSONLY
      LOGICAL ::             DO_SPECIALIST_OPTION_1
      LOGICAL ::             DO_SPECIALIST_OPTION_2
      LOGICAL ::             DO_SPECIALIST_OPTION_3
      LOGICAL ::             DO_TOA_CONTRIBS


      if ( present(scalar) ) then
         vector = .not. scalar
      else
         vector = .true. 
      end if
      rc = 0



!    initialize variables to false or zero
!    -------------------------------------
      call VLIDORT_INIT_CONTROL_VARS (DO_FULLRAD_MODE, DO_SSCORR_NADIR,
     &     DO_SSCORR_OUTGOING, DO_SSCORR_TRUNCATION, DO_SSFULL,
     &     DO_DOUBLE_CONVTEST, DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, 
     &     DO_CHAPMAN_FUNCTION, 
     &     DO_RAYLEIGH_ONLY, DO_DELTAM_SCALING, DO_SOLUTION_SAVING, 
     &     DO_BVP_TELESCOPING, DO_UPWELLING, DO_DNWELLING, DO_QUAD_OUTPUT, 
     &     DO_USER_VZANGLES, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, 
     &     DO_DEBUG_WRITE, DO_WRITE_INPUT, DO_WRITE_SCENARIO, 
     &     DO_WRITE_FOURIER, DO_WRITE_RESULTS, DO_LAMBERTIAN_SURFACE, 
     &     DO_SPECIALIST_OPTION_1, 
     &     DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3, DO_TOA_CONTRIBS)


      call VLIDORT_INIT_THERMAL_VARS(DO_THERMAL_EMISSION, N_THERMAL_COEFFS, 
     &     THERMAL_BB_INPUT, DO_SURFACE_EMISSION, 
     &     SURFBB, DO_THERMAL_TRANSONLY)

      call VLIDORT_INIT_MODEL_VARS(NSTOKES, NSTREAMS, NLAYERS, NFINELAYERS, 
     &     NGREEK_MOMENTS_INPUT, VLIDORT_ACCURACY,  
     &     FLUX_FACTOR, N_SZANGLES, SZANGLES,       
     &     EARTH_RADIUS, RFINDEX_PARAMETER,       
     &     GEOMETRY_SPECHEIGHT, N_USER_RELAZMS,     
     &     USER_RELAZMS, N_USER_VZANGLES,           
     &     USER_VZANGLES, N_USER_LEVELS, USER_LEVELS,
     &     LAMBERTIAN_ALBEDO)



!     Customize those parameters that do not vary from pixel to pixel
!     ---------------------------------------------------------------


!                         Modes of Operation
!                         ------------------

      DO_FULLRAD_MODE    = .true.  ! Do full Stokes vector calculation?
      DO_SSCORR_NADIR    = .false. ! Do nadir single scatter correction?
      DO_SSCORR_OUTGOING = .true.  ! Do outgoing single scatter correction?
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
      DO_DELTAM_SCALING    = .false.  ! Include Delta-M scaling?
      DO_SSCORR_TRUNCATION = .false. ! Additional Delta-M scaling for SS correction?
      DO_SOLUTION_SAVING   = .false. ! Solution saving mode?
      DO_BVP_TELESCOPING   = .false. ! Boundary value problem telescoping mode?

!                      User-defined output control
!                      ---------------------------

      DO_UPWELLING = .true.     ! Upwelling output?
      DO_DNWELLING = .false.    ! Downwelling output?
      DO_USER_VZANGLES = .true. ! User-defined viewing zenith angles?

      DO_ADDITIONAL_MVOUT = .false.  ! Generate mean value output additionally?
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

! WAS NOT THERE in the previous version (f77), do we need them??????

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


!                     Stokes/streams/layers/moments
!                     -----------------------------

      if ( vector ) then
         NSTOKES = 3        ! Number of Stokes vector components
      else
         NSTOKES = 1
      end if
      NSTREAMS = 6         ! Number of half-space streams
      NLAYERS = km          ! Number of atmospheric layers
      NFINELAYERS = 2       ! Number of fine layers (outgoing sphericity correction)
      NGREEK_MOMENTS_INPUT = 500 ! Number of scattering matrix expansion coefficients
 
!      print *, 'NSTOKES', NSTOKES

!                            accuracy input
!                            --------------

      VLIDORT_ACCURACY = 0.0001  ! Fourier series convergence
      FLUX_FACTOR = 1.00d0      ! Solar flux constant; =1 if no solar sources.

!                        Pseudo-spherical inputs
!                        -----------------------

      EARTH_RADIUS = 6371.0 ! Earth radius (km)
      RFINDEX_PARAMETER = 0.000288 ! Refractive index parameter
      GEOMETRY_SPECHEIGHT = 0.0 ! Input geometry specification height [km]

!                       Solar/Viewing Geometry 
!                       ----------------------

      N_SZANGLES = 1       ! Number of BOA solar zenith angles 
      N_USER_RELAZMS = 1   ! Number of user defined relative azimuth angles    
      N_USER_VZANGLES = 1  ! Number of user-defined stream angles     
      N_USER_LEVELS = 1    ! Number of user-defined vertical output levels

!                          Lambertian Only
!                          ---------------

      DO_LAMBERTIAN_SURFACE = .true. ! Lambertian surface?
      LAMBERTIAN_ALBEDO = 0.0  ! Lambertian (isotropic) input albedo, will be reset
                                     !  during run
      

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

      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_FULLRAD_MODE        = DO_FULLRAD_MODE
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_SSCORR_TRUNCATION   = DO_SSCORR_TRUNCATION
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_SSFULL              = DO_SSFULL
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_THERMAL_EMISSION    = DO_THERMAL_EMISSION
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_SURFACE_EMISSION    = DO_SURFACE_EMISSION
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_PLANE_PARALLEL      = DO_PLANE_PARALLEL
      !self%VIO%VLIDORT_FBoolean_inputs%TS_DO_BRDF_SURFACE        = DO_BRDF_SURFACE
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_UPWELLING           = DO_UPWELLING
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_DNWELLING           = DO_DNWELLING
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_QUAD_OUTPUT         = DO_QUAD_OUTPUT
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_TOA_CONTRIBS        = DO_TOA_CONTRIBS
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_LAMBERTIAN_SURFACE  = DO_LAMBERTIAN_SURFACE
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_SPECIALIST_OPTION_1 = DO_SPECIALIST_OPTION_1
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_SPECIALIST_OPTION_2 = DO_SPECIALIST_OPTION_2
      self%VIO%VLIDORT_FBoolean_inputs%TS_DO_SPECIALIST_OPTION_3 = DO_SPECIALIST_OPTION_3

!  Modified Boolean inputs

      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_SSCORR_NADIR        = DO_SSCORR_NADIR
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_SSCORR_OUTGOING     = DO_SSCORR_OUTGOING
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_DOUBLE_CONVTEST     = DO_DOUBLE_CONVTEST
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_SOLAR_SOURCES       = DO_SOLAR_SOURCES
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_REFRACTIVE_GEOMETRY = DO_REFRACTIVE_GEOMETRY
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_CHAPMAN_FUNCTION    = DO_CHAPMAN_FUNCTION
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_RAYLEIGH_ONLY       = DO_RAYLEIGH_ONLY
      !VLIDORT_MBoolean_inputs%TS_DO_ISOTROPIC_ONLY      = DO_ISOTROPIC_ONLY
      !VLIDORT_MBoolean_inputs%TS_DO_NO_AZIMUTH          = DO_NO_AZIMUTH
      !VLIDORT_MBoolean_inputs%TS_DO_ALL_FOURIER         = DO_ALL_FOURIER
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_DELTAM_SCALING      = DO_DELTAM_SCALING
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_SOLUTION_SAVING     = DO_SOLUTION_SAVING
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_BVP_TELESCOPING     = DO_BVP_TELESCOPING
      !VLIDORT_MBoolean_inputs%TS_DO_USER_STREAMS        = DO_USER_STREAMS
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_USER_VZANGLES       = DO_USER_VZANGLES
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_ADDITIONAL_MVOUT    = DO_ADDITIONAL_MVOUT
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_MVOUT_ONLY          = DO_MVOUT_ONLY
      self%VIO%VLIDORT_MBoolean_inputs%TS_DO_THERMAL_TRANSONLY   = DO_THERMAL_TRANSONLY

!  Fixed control inputs

      self%VIO%VLIDORT_FControl_inputs%TS_NSTOKES          = NSTOKES
      self%VIO%VLIDORT_FControl_inputs%TS_NSTREAMS         = NSTREAMS
      self%VIO%VLIDORT_FControl_inputs%TS_NLAYERS          = NLAYERS
      self%VIO%VLIDORT_FControl_inputs%TS_NFINELAYERS      = NFINELAYERS
      self%VIO%VLIDORT_FControl_inputs%TS_N_THERMAL_COEFFS = N_THERMAL_COEFFS
      self%VIO%VLIDORT_FControl_inputs%TS_VLIDORT_ACCURACY = VLIDORT_ACCURACY

!  Modified control inputs

      self%VIO%VLIDORT_MControl_inputs%TS_NGREEK_MOMENTS_INPUT = NGREEK_MOMENTS_INPUT

!  Beam inputs

      self%VIO%VLIDORT_Beam_inputs%TS_FLUX_FACTOR = FLUX_FACTOR
      !VLIDORT_Beam_inputs%TS_NBEAMS      = NBEAMS
      self%VIO%VLIDORT_Beam_inputs%TS_N_SZANGLES  = N_SZANGLES
      !VLIDORT_Beam_inputs%TS_BEAM_SZAS   = BEAM_SZAS
      self%VIO%VLIDORT_Beam_inputs%TS_SZANGLES    = SZANGLES

!  User value inputs

      self%VIO%VLIDORT_User_inputs%TS_N_USER_RELAZMS      = N_USER_RELAZMS
      self%VIO%VLIDORT_User_inputs%TS_USER_RELAZMS        = USER_RELAZMS

      !VLIDORT_User_inputs%TS_N_USER_STREAMS      = N_USER_STREAMS
      self%VIO%VLIDORT_User_inputs%TS_N_USER_VZANGLES     = N_USER_VZANGLES
      !VLIDORT_User_inputs%TS_USER_ANGLES_INPUT   = USER_ANGLES
      self%VIO%VLIDORT_User_inputs%TS_USER_VZANGLES_INPUT = USER_VZANGLES

      self%VIO%VLIDORT_User_inputs%TS_N_USER_LEVELS       = N_USER_LEVELS
      self%VIO%VLIDORT_User_inputs%TS_USER_LEVELS         = USER_LEVELS

      self%VIO%VLIDORT_User_inputs%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT

!  Fixed Chapman function inputs

      !VLIDORT_Chapman_inputs%TS_HEIGHT_GRID       = HEIGHT_GRID
      !VLIDORT_Chapman_inputs%TS_PRESSURE_GRID     = PRESSURE_GRID
      !VLIDORT_Chapman_inputs%TS_TEMPERATURE_GRID  = TEMPERATURE_GRID
      !VLIDORT_Chapman_inputs%TS_FINEGRID          = FINEGRID
      self%VIO%VLIDORT_Chapman_inputs%TS_EARTH_RADIUS      = EARTH_RADIUS
      self%VIO%VLIDORT_Chapman_inputs%TS_RFINDEX_PARAMETER = RFINDEX_PARAMETER

!  Fixed optical inputs

      !VLIDORT_Optical_inputs%TS_DELTAU_VERT_INPUT    = DELTAU_VERT_INPUT
      !VLIDORT_Optical_inputs%TS_OMEGA_TOTAL_INPUT    = OMEGA_TOTAL_INPUT
      !VLIDORT_Optical_inputs%TS_GREEKMAT_TOTAL_INPUT = GREEKMAT_TOTAL_INPUT
      self%VIO%VLIDORT_Optical_inputs%TS_THERMAL_BB_INPUT     = THERMAL_BB_INPUT
      self%VIO%VLIDORT_Optical_inputs%TS_LAMBERTIAN_ALBEDO    = LAMBERTIAN_ALBEDO
      self%VIO%VLIDORT_Optical_inputs%TS_SURFACE_BB_INPUT     = SURFBB

!  Fixed write inputs

      self%VIO%VLIDORT_Write_inputs%TS_DO_DEBUG_WRITE          = DO_DEBUG_WRITE

      self%VIO%VLIDORT_Write_inputs%TS_DO_WRITE_INPUT          = DO_WRITE_INPUT
      self%VIO%VLIDORT_Write_inputs%TS_INPUT_WRITE_FILENAME    = INPUT_WRITE_FILENAME

      self%VIO%VLIDORT_Write_inputs%TS_DO_WRITE_SCENARIO       = DO_WRITE_SCENARIO
      self%VIO%VLIDORT_Write_inputs%TS_SCENARIO_WRITE_FILENAME = SCENARIO_WRITE_FILENAME

      self%VIO%VLIDORT_Write_inputs%TS_DO_WRITE_FOURIER        = DO_WRITE_FOURIER
      self%VIO%VLIDORT_Write_inputs%TS_FOURIER_WRITE_FILENAME  = FOURIER_WRITE_FILENAME

      self%VIO%VLIDORT_Write_inputs%TS_DO_WRITE_RESULTS        = DO_WRITE_RESULTS
      self%VIO%VLIDORT_Write_inputs%TS_RESULTS_WRITE_FILENAME  = RESULTS_WRITE_FILENAME



!      RETURN
!      print *, 'maxSTOKES', MAXSTOKES

!                            Error Checking
!                            --------------

      if ( NSTOKES  .GT. MAXSTOKES  )                   rc = 1
      if ( NSTREAMS .GT. MAXSTREAMS )                   rc = 2
      if ( NLAYERS  .GT. MAXLAYERS  )                   rc = 3
      if ( NFINELAYERS .GT. MAXFINELAYERS )             rc = 4
      if ( NGREEK_MOMENTS_INPUT .GT. MAXMOMENTS_INPUT)  rc = 5
      if ( N_SZANGLES .GT. MAX_SZANGLES )               rc = 6
      if ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS )       rc = 7
      if ( N_USER_VZANGLES .GT. MAX_USER_VZANGLES )     rc = 8
      if ( N_USER_LEVELS .GT. MAX_USER_LEVELS )         rc = 9 
      if ( N_THERMAL_COEFFS .GT. MAX_THERMAL_COEFFS )   rc =  10
     

!     All done
!     --------
      self%initialized = .true.

      end subroutine VLIDORT_LambInit

!.............................................................................
 
      subroutine VLIDORT_LambRun (self,radiance, reflectance, AOT, rc, 
     &                            scalar, aerosol)
!
!     Computes radiances for a single wavelength, pixel. Optical properties
!     and met fields in self are assumed to have been updated with the
!     apropriate values.
!
      USE VLIDORT_PARS
      
      USE VLIDORT_BRDFSup_def
      USE VLIDORT_Inputs_def
      USE VLIDORT_Outputs_def

      USE VLIDORT_AUX
      USE VLIDORT_INPUTS
      USE VLIDORT_MASTERS
     

      type(VLIDORT_Lamb),         intent(inout)   :: self        ! Contains most input
      real*8,             intent(out)     :: radiance    ! TOA radiance
      real*8,             intent(out)     :: reflectance ! TOA reflectance
      real*8,             intent(out)     :: AOT         ! Aerosol Optical Thickness
      integer,            intent(out)     :: rc

      logical, optional,  intent(in)      :: scalar  ! If True, do scalar calculation
      logical,            intent(in)      :: aerosol ! if False, Rayleigh only

!                           ----

      integer :: STATUS_INPUTCHECK, STATUS_CALCULATION 
      logical :: vector

!                           ----

!     local variables
!     ---------------
      integer                                          :: i, j, k, l, m, n
      integer                                          :: ncoeffs, k1, IDR, ierror
      integer                                          :: NLAYERS
      real*8                                           :: ray_l      
      real*8                                           :: tau_l 
      real*8                                           :: ssa_l 
      real*8                                           :: g_l 
      real*8                                           :: tau_ext
      real*8                                           :: tau_scat
      real*8                                           :: ssa_tot
      real*8                                           :: raysmom2
      real*8                                           :: gammamom2
      real*8                                           :: alphamom2 
      real*8                                           :: deltamom1 
      real*8                                           :: aerswt 
      real*8                                           :: rayswt
      real*8                                           :: factor      
      real*8                                           :: COEF_DEPOL
      real*8                                           :: wmicron

   
      real*8                                           :: x
      real*8, dimension(MAXLAYERS)                       :: Ray
      real*8, dimension(0:MAXMOMENTS_INPUT)            :: aersmom 
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: aervmoms 
      real*8, dimension(0:2, 16)                       :: rayvmoms
      real*8                                           :: difz
      real*8                                           :: somray
      logical                                          :: DO_UPWELLING
      integer                                          :: NGREEK_MOMENTS_INPUT
      real*8, dimension(0:MAXMOMENTS_INPUT,MAXLAYERS,16) :: greekmat_total_input
      real*8, dimension(0:MAXLAYERS)                     :: height_grid                     
      real*8, dimension(0:MAXLAYERS)                     :: pressure_grid 
      real*8, dimension(0:MAXLAYERS)                     :: temperature_grid

      real*8, dimension(MAXLAYERS)                       :: deltau_vert_input
      real*8, dimension(MAXLAYERS)                       :: omega_total_input

      real*8, dimension(MAX_USER_LEVELS, MAX_GEOMETRIES
     &, MAXSTOKES, MAX_DIRECTIONS)                     :: STOKES
      real*8                                           :: FLUX_FACTOR

!     Volume Rayleigh scattering coefficient (depends closely on the 
!     thermodynamic conditions of the atm and varies with height, 
!     product of the molecular number density of air by the total ray 
!     scattering cross section).
      real*8, dimension(MAXLAYERS+1)                     :: Vol 
      real*8, dimension(MAXLAYERS+1)                     :: sect 
      real*8, parameter                                :: pi = 3.14159
      real*8, parameter                                :: DEPOL_RATIO = 0.030
      
      
       rc = 0
 
       if ( .not. self%initialized ) then
         rc = 1
         return
       end if
      
       if ( present(scalar) ) then
         vector = .not. scalar
      else
         vector = .true.
      end if

      if ( self%VIO%VLIDORT_FBoolean_inputs%TS_DO_UPWELLING ) then
         IDR = 1
      else
         IDR = 2
      end if
     
!      print *, 'IDR',IDR

      self%VIO%VLIDORT_Beam_inputs%TS_SZANGLES(1) = self%solar_zenith
      self%VIO%VLIDORT_User_inputs%TS_USER_RELAZMS(1) = self%relative_azymuth
      self%VIO%VLIDORT_User_inputs%TS_USER_VZANGLES_INPUT(1) = self%sensor_zenith
      self%VIO%VLIDORT_User_inputs%TS_USER_LEVELS(1) = 0.0
      self%VIO%VLIDORT_Optical_inputs%TS_LAMBERTIAN_ALBEDO = self%albedo  
      
      NGREEK_MOMENTS_INPUT = self%VIO%VLIDORT_MControl_inputs%TS_NGREEK_MOMENTS_INPUT
      NLAYERS = self%VIO%VLIDORT_FControl_inputs%TS_NLAYERS     
      print*,'ALBEDO', self%wavelength,self%VIO%VLIDORT_Optical_inputs%TS_LAMBERTIAN_ALBEDO
      print*,'SZA', self%wavelength,self%VIO%VLIDORT_Beam_inputs%TS_SZANGLES(1)
!      print *, 'NLAYERS',NLAYERS, self%VIO%VLIDORT_Optical_inputs%TS_LAMBERTIAN_ALBEDO

!               Calculation of the Rayleigh-Scattering Optical Depth
!               ----------------------------------------------------
   
      height_grid(0) = self%ze(1) * 1.E-3  ! en km
      pressure_grid(0) = self%pe(1) * 1.E-2 ! en hPa
      temperature_grid(0) = self%te(1)
!      print *, 'height, pres, t,',height_grid(0),pressure_grid(0),temperature_grid(0)

      do i = 1, NLAYERS 
         height_grid(i) = self%ze(i+1) * 1.E-3 ! en km
         pressure_grid(i) = self%pe(i+1) * 1.E-2 ! en hPa
         temperature_grid(i) = self%te(i+1)
      end do
            
      self%VIO%VLIDORT_Chapman_inputs%TS_height_grid = height_grid
      self%VIO%VLIDORT_Chapman_inputs%TS_pressure_grid = pressure_grid
      self%VIO%VLIDORT_Chapman_inputs%TS_temperature_grid = temperature_grid


!     Rayleigh extinction profile from Bodhaine et al., (1999) 
!     and Tomasi et al., (2005)
!     (wavelength in micrometer, pressure in hpa)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      wmicron = self%wavelength * 1.E-3  ! micrometer
      
      do j = 0, NLAYERS 
        
         Vol(j) = 3.69296E-18 * A(wmicron) * pressure_grid(j) 
     &          / temperature_grid(j)
         Vol(j) = Vol(j) * 1.E5 ! en km-1

        sect(j) = Vol(j)/2.546899E19 * 1013.25/pressure_grid(j) 
     &          * temperature_grid(j)/288.15

!         print*, "sect", j,height_grid(j), sect(j) ! test de compa with Bodhaine
           
      end do
       
!     logarithmique interpolation procedure of the Rayleigh profile 
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      somray = 0.
     
      do j = 0, NLAYERS-1
         difz = height_grid(j) - height_grid(j+1)
!         ray(j+1) = difz * Vol(j) * (Vol(j+1)/Vol(j))**0.5 
         ray(j+1) = difz * (Vol(j) + Vol(j+1))/2.
         somray = somray + ray(j+1) ! Total Tau Ray
      end do

!      print*, 'somray', somray
!     greek moments for Rayleigh only( only if vector )
!     ----------------------------- 
      if ( vector ) then

         gammamom2 = -SQRT(6.) * (1 - DEPOL_RATIO) / (2 + DEPOL_RATIO)
         alphamom2 = 6 * (1 - DEPOL_RATIO) / (2 + DEPOL_RATIO)
         deltamom1 = 3 * (1 - 2 * DEPOL_RATIO) / (2 + DEPOL_RATIO)
         raysmom2 = (1.0 - DEPOL_RATIO)/(2.0 + DEPOL_RATIO) 

!DEPOL_RATIO = (6 * COEF_DEPOL - 6) / (3 + 7 * COEF_DEPOL) 
      
         rayvmoms(0,1) = 1.0
         rayvmoms(1,1) = 0.0
         rayvmoms(2,1) = raysmom2
         rayvmoms(0,2) = 0.0
         rayvmoms(1,2) = 0.0
         rayvmoms(2,2) = gammamom2
         do k = 3, 4
            do l = 0, 2
               rayvmoms(l,k) = 0.0
            end do
         end do
         rayvmoms(0,5) = 0.0
         rayvmoms(1,5) = 0.0
         rayvmoms(2,5) = gammamom2 
         rayvmoms(0,6) = 0.0
         rayvmoms(1,6) = 0.0
         rayvmoms(2,6) = alphamom2 
         do k = 7, 15
            do l = 0, 2
               rayvmoms(l,k) = 0.0
            end do
         end do
         rayvmoms(0,16) = 0.0
         rayvmoms(1,16) = deltamom1
         rayvmoms(2,16) = 0.0 
     
      end if

!     Loop over the layers:
!     ---------------------
      AOT = 0.0
      PRINT*, 'AOD', self%tau(1)
      do i = 1, NLAYERS  
         ray_l = ray(i)         ! indice l for  each layer

         if ( .not. aerosol ) then
            tau_l = 0.0 
            ssa_l = 0.0 
            g_l = 0.0 
         else
            tau_l = self%tau(i)
            ssa_l = self%ssa(i) 
            g_l = self%g(i) 
            AOT = AOT + tau_l
         end if

!        total optical depths for extinction and scattering 
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         tau_ext = ray_l + tau_l
         tau_scat = ray_l +  ssa_l * tau_l

!        single scattering albedo total
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ssa_tot = tau_scat / tau_ext
         if ( ssa_tot > 0.99999 ) then
            ssa_tot = 0.99999
         end if
     
         deltau_vert_input(i) = tau_ext
         omega_total_input(i) = ssa_tot 

!            Compute Henyey-Greenstein phase function, including Rayleigh
!            ------------------------------------------------------------

!     SCALAR testing Rayleigh second moment
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      wmicron = self%wavelength * 1.E-3  ! micrometer
     
      raysmom2 = (1.0 - DEPOL_RATIO)/(2.0 + DEPOL_RATIO) 
      if (scalar) then    
    
!     Phase function moments (Rayleigh only)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (tau_l == 0.0) then
     
            greekmat_total_input(0,i,1) = 1.0
            greekmat_total_input(1,i,1) = 0.0
            greekmat_total_input(2,i,1) = raysmom2
            do l = 3, NGREEK_MOMENTS_INPUT        
               greekmat_total_input(l,i,1) = 0.0
            end do 
         end if
     

!     SCALAR phase function moments (aerosol + Rayleigh)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if (tau_l /= 0.0) then
             aerswt = ssa_l * tau_l / tau_scat  
             rayswt = ray_l / tau_scat    
             aersmom(0) = 1.0
             aersmom(1) = 3.0 * g_l
             aersmom(2) = 5.0 * g_l * aersmom(1) / 3.0

             greekmat_total_input(0,i,1) = 1.0
             greekmat_total_input(1,i,1) = aersmom(1) * aerswt
             greekmat_total_input(2,i,1) = raysmom2 * rayswt + aersmom(2) * aerswt
         
             do l = 3, NGREEK_MOMENTS_INPUT         
                factor = REAL(2*l+1) / REAL(2*l-1) 
                aersmom(l) = factor * g_l * aersmom(l-1)
                greekmat_total_input(l,i,1) = aersmom(l) * aerswt
             end do  
          end if ! end if tau_l /= 0.0
      end if  ! end scalar
      
      
!     VECTOR phase function moments 
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if ( vector ) then  
      
!     Phase function moments (Rayleigh only)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (tau_l == 0.0) then
   
           do k=1, 16
             do l=0, 2
              greekmat_total_input(l,i,k) = rayvmoms(l,k)
             end do
           end do

            do k = 1, 16 
              do l = 3, NGREEK_MOMENTS_INPUT
             greekmat_total_input(l,i,k) = 0.0
             end do
           end do

      end if
      
      
!     VECTOR phase function moments (aerosol + Rayleigh)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (tau_l /= 0.0) then 
             
         aerswt = ssa_l * tau_l / tau_scat  
         rayswt = ray_l / tau_scat 

!         if (g_l > 0.9.and.g_l<=1.0) then
!            self%nmom = 500
!         else if (g_l > 0.8.and.g_l<=0.9) then
!           self%nmom = 200
!         else 
!            self%nmom = 100
!         end if
!        print*, 'ok1', self%nmom
             
         if ( self%nmom -1 < NGREEK_MOMENTS_INPUT ) then 
            do l= 0, self%nmom-1
               
               aervmoms(l,i,1)  = self%pmom(i,l+1,1) ! P11             
               aervmoms(l,i,6)  = self%pmom(i,l+1,1) ! P11               
               aervmoms(l,i,2)  = self%pmom(i,l+1,2) ! P12                
               aervmoms(l,i,5)  = self%pmom(i,l+1,2) ! P12             
               aervmoms(l,i,11) = self%pmom(i,l+1,3) ! P33           
               aervmoms(l,i,16) = self%pmom(i,l+1,3) ! P33          
               aervmoms(l,i,12) = self%pmom(i,l+1,4) ! P34            
               aervmoms(l,i,15) = -self%pmom(i,l+1,4) ! - P34 = P44
            end do
         
            do l = self%nmom, NGREEK_MOMENTS_INPUT
               aervmoms(l,i,1)  = 0.0 ! P11
               aervmoms(l,i,6)  = 0.0 ! P11 = P22
               aervmoms(l,i,2)  = 0.0 ! P12
               aervmoms(l,i,5)  = 0.0 ! P12
               aervmoms(l,i,11) = 0.0 ! P33
               aervmoms(l,i,16) = 0.0 ! P33
               aervmoms(l,i,12) = 0.0 ! P34
               aervmoms(l,i,15) = 0.0 ! P34 = P44
        
            end do
        
         else 
         
            do l= 0, NGREEK_MOMENTS_INPUT
               aervmoms(l,i,1)  = self%pmom(i,l+1,1) ! P11
               aervmoms(l,i,6)  = self%pmom(i,l+1,1) ! P11
               aervmoms(l,i,2)  = self%pmom(i,l+1,2) ! P12
               aervmoms(l,i,5)  = self%pmom(i,l+1,2) ! P12
               aervmoms(l,i,11) = self%pmom(i,l+1,3) ! P33
               aervmoms(l,i,16) = self%pmom(i,l+1,3) ! P33
               aervmoms(l,i,12) = self%pmom(i,l+1,4) ! P34
               aervmoms(l,i,15) = -self%pmom(i,l+1,4) ! - P34 = P44
            end do
       
         end if

         do k = 1, 16
            do l = 0,2
               greekmat_total_input(l,i,k) = rayvmoms(l,k) * rayswt 
     &                                     + aervmoms(l,i,k) * aerswt
            end do
           
            do l = 3, NGREEK_MOMENTS_INPUT
               greekmat_total_input(l,i,k) = aervmoms(l,i,k) * aerswt
               
            end do
                  
         end do
         greekmat_total_input(0,i,1) = 1.0
      
     
      end if                    ! end if tau_l /= 0.0

      end if                    ! end if vector
      
!     end layer loop
!     ---------------
      end do
 
      self%VIO%VLIDORT_Optical_inputs%TS_DELTAU_VERT_INPUT = deltau_vert_input
      self%VIO%VLIDORT_Optical_inputs%TS_OMEGA_TOTAL_INPUT = omega_total_input
      self%VIO%VLIDORT_Optical_inputs%TS_GREEKMAT_TOTAL_INPUT = greekmat_total_input
       
!      print*, 'DELTAU',deltau_vert_input 
!      print*, 'OMEGA',omega_total_input

      do i=1, NGREEK_MOMENTS_INPUT
      print*,  'GREEK',i,greekmat_total_input(i,72,1)
      enddo

!      g_tot = g_tot/NLAYERS
     

!     Call the Lambertian driver for doing the actual calculation
!     -----------------------------------------------------------
      call VLIDORT_MASTER (self%VIO%VLIDORT_FBoolean_inputs,
     &     self%VIO%VLIDORT_MBoolean_inputs, 
     &     self%VIO%VLIDORT_FControl_inputs, 
     &     self%VIO%VLIDORT_MControl_inputs, 
     &     self%VIO%VLIDORT_Beam_inputs,     
     &     self%VIO%VLIDORT_User_inputs,     
     &     self%VIO%VLIDORT_Chapman_inputs,  
     &     self%VIO%VLIDORT_Optical_inputs,  
     &     self%VIO%VLIDORT_Write_inputs,    
     &     self%VIO%VLIDORT_BRDF_inputs,     
     &     self%VIO%VLIDORT_Outputs,         
     &     self%VIO%VLIDORT_Status )

      if ( self%VIO%VLIDORT_Status%TS_STATUS_INPUTCHECK /= 0 ) rc = 1
      if ( self%VIO%VLIDORT_Status%TS_STATUS_CALCULATION /= 0 ) rc = 2
      if ( rc /= 0 ) return
       
      STOKES = self%VIO%VLIDORT_Outputs%TS_STOKES ! output of VLIDORT_MASTER subroutine
      FLUX_FACTOR = self%VIO%VLIDORT_Beam_inputs%TS_FLUX_FACTOR
      print *, 'FF', FLUX_FACTOR
!     Return TOA radiance
!     -------------------
      RADIANCE = 0.0
      REFLECTANCE = 0.0
      if ( scalar ) then
       RADIANCE = STOKES(1, 1, 1, IDR)
              
       REFLECTANCE = (pi * RADIANCE) / ( cos(self%VIO%VLIDORT_Beam_inputs%TS_SZANGLES(1)*pi/180.0) * FLUX_FACTOR ) 
      end if
  
      if ( vector ) then
       RADIANCE = STOKES(1, 1, 1, IDR)
!       print *, 'RAD', RADIANCE
       REFLECTANCE = (pi * RADIANCE) / ( cos(self%VIO%VLIDORT_Beam_inputs%TS_SZANGLES(1)*pi/180.0) * FLUX_FACTOR )   
      end if
      print *, 'RADIANCE',self%wavelength, RADIANCE 

      STOP 
      Contains  

!.................................................................................
! determination of the function A(wmicron)
!-----------------------------------      
    
         function A(X)
!  pour le Rayleigh, X = lambda en micron
       implicit none
       real*8 :: A
       real*8, intent(in) :: X
       real*8 :: depol_ratio
       real*8 :: XX, XX2, C
       real*8 :: depol1, depol2, depol3, depol4
       real*8 :: coef_depol
       real*8 :: RI

      XX=1./X

      XX2=XX*XX
     
!  ns-1  sans approximation
      RI=8060.77+2481070/(132.274-XX2)+17456.3/(39.32957-XX2)


! autre maniere de determiner coef_depol
        depol1 = 1.034 + 3.17E-4 * XX2
        depol2 = 1.096 + 1.385E-3 * XX2 + 1.448E-4 * XX2 * XX2
        depol3 = 1.
        depol4 = 1.15
!depol5 = 1.001 if we 
        C = 0.030 ! for standard air , concentration of CO2 = 300 ppmv

        coef_depol= (78.084 * depol1 + 20.946 * depol2 + 0.934 * depol3 +  C * depol4) / (78.084 + 20.946 + 0.934 + C)

         
      A=coef_depol * XX2 * XX2 * RI * RI
      
      return

      end function A      
    
      end subroutine VLIDORT_LambRun

!.............................................................................

      subroutine VLIDORT_LambFin (self, rc)
      type(VLIDORT_Lamb), intent(inout) :: self
      integer,              intent(out) :: rc
      self%initialized = .false.
      
      rc = 0
      
      end subroutine VLIDORT_LambFin

!..............................................................................

      subroutine VLIDORT_LER (self, rad, refl, rc)
 
      ! First compute the spherical albedo and the transmission of the 
      ! atmosphere purely rayleigh and
      ! the Lambertian Equivalent Reflectivity R*
      !-------------------------------------------------------------------------------------

       type(VLIDORT_Lamb), intent(inout)           :: self 
       real*8,  intent(in)               :: rad       !Radiance mesure (atm Ray + aer)
       real*8,  intent(out)              :: refl      !surface reflectivity R*           
       integer, intent(out)              :: rc
       real*8                            :: spher_alb !spherical albedo
       real*8                            :: trans     !transmission
       real*8                            :: reflectance 
       real*8                            :: beta
       real*8                            :: rad0      !Ray atmo and no surface reflection
       real*8                            :: rad1      !Ray atmo and surface albedo = 0.1
       real*8                            :: rad2      !Ray atmo and surface albedo = 0.2
       real*8                            :: rad_aer
       real*8                            :: rad_ray
       real*8                            :: delta_rad
       real*8                            :: AOT
      
       
      
       ! atmosphere Rayleigh and surface albedo = 0.0
       ! --------------------------------------------
       self%albedo = 0.0
       call VLIDORT_LambRun (self, rad0, reflectance, AOT, rc, .false., .false.)

       ! atmosphere Rayleigh and surface albedo = 0.1
       ! --------------------------------------------
       self%albedo = 0.1
       call VLIDORT_LambRun (self, rad1, reflectance, AOT, rc, .false., .false.)

       ! atmosphere Rayleigh and surface albedo = 0.2
       ! --------------------------------------------
       self%albedo = 0.2
       call VLIDORT_LambRun (self, rad2, reflectance, AOT, rc, .false., .false.)

       beta =  ( rad1 - rad0 ) / ( rad2 - rad0 )
       spher_alb = (( beta * 0.2 ) - 0.1 ) / ( 0.1 * 0.2 * ( beta - 1 )) ! Spherical albedo
       trans = ( 1 - ( spher_alb * 0.1 )) * ( rad1 - rad0 ) / 0.1        ! transmission
       
!       print*, 'spher_alb', spher_alb, trans
 
       ! Compute the Reflectivity R*
       ! --------------------------
       ! Compute the radiance for an atmosphere purely Ray and surface albedo = 0.0
       !----------------------------------------------------------------------------
       self%albedo = 0.0
       call VLIDORT_LambRun (self, rad_ray, reflectance, AOT, rc, .false., .false.)
!       print*, 'rad_ray', rad_ray      

       delta_rad =  rad - rad_ray    
        
       refl = delta_rad / ( trans + ( spher_alb * delta_rad )) ! Reflectivity R*

!       print*, 'refl', refl

       end subroutine VLIDORT_LER
       
!.........................................................................

      subroutine VLIDORT_AI (self,rad, refl, AI, rc)   

      ! Compute the Aerosol Index from the reflectivity calculate 
      ! from the previous subroutine
      ! ---------------------------------------------------------
       type(VLIDORT_Lamb),         intent(inout)   :: self
       real*8,             intent(in)     :: rad
       real*8,             intent(in)     :: refl      !reflectivity R*
       integer,            intent(out)    :: rc
       real*8,             intent(out)    :: AI        !Aerosol Index
       real*8                             :: rad_nouv  !calcul with R* as the surface albedo param
       real*8                             :: reflectance 
       real*8                             :: alpha  
       
       real*8                             :: AOT 
       real*8                             :: g_tot
      ! Compute the radiance for an atmosphere Ray with the R* as the surface alb parameter
      !-----------------------------------------------------------------------------------

      self%albedo = refl  
      
      call VLIDORT_LambRun (self, rad_nouv, reflectance, AOT, rc, .false., .false.)

      if (rad > 0.0) then
      alpha = rad / rad_nouv
      AI = -100 * log10(alpha) ! Calcul of the Aerosol Index
      else
      AI = 999.999
      end if
     
      
      end subroutine VLIDORT_AI


!................................................................................
! function coef_depolarisation ( KING FACTOR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      function COEF_DEPOL(X)
!      implicit none
!      real*8              :: COEF_DEPOL ! King Factor
!      real*8              :: DEPOL_RATIO ! depolarisation ratio
!      real*8, intent(in)  :: X
!      real*8              :: XX, XX2, C
!      real*8              :: DEPOL1, DEPOL2, DEPOL3, DEPOL4
      
!      XX=1./X

!      XX2=XX*XX
     
!     Calculation of the King Factor
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      DEPOL1 = 1.034 + 3.17E-4 * XX2
!      DEPOL2 = 1.096 + 1.385E-3 * XX2 + 1.448E-4 * XX2 * XX2
!      DEPOL3 = 1.
!      DEPOL4 = 1.15
!      C = 0.030 ! for standard air , concentration of CO2 = 300 ppmv

!      COEF_DEPOL= (78.084 * DEPOL1 + 20.946 * DEPOL2 + 0.934 * DEPOL3 +  C * DEPOL4) / (78.084 + 20.946 + 0.934 + C)
    
!      return

!      end function COEF_DEPOL

      end module VLIDORT90_LambMod
