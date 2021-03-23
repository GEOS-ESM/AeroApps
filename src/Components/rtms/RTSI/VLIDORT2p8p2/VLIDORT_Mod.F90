      module VLIDORT_Mod

      USE VLIDORT_PARS_m
      
      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_SETUPS_DEF_m
      USE VBRDF_SUP_MOD_m
      USE VSLEAVE_SUP_MOD_m
      
      implicit NONE

	      PUBLIC  VLIDORT_Init ! Initialize once for all pixels
	      PUBLIC  VLIDORT_Fin  ! Finalize once
	       
	      !  VLIDORT input/output structures

	      TYPE VLIDORT_IO
		  !  VLIDORT input structures
		  TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn
		  TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn



		  !  VLIDORT output structures
	!!          TYPE(Outputs_Main_def)             :: VLIDORT_Outputs
		  TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

	!!          TYPE(Exception_Handling_def)       :: VLIDORT_Status
		  TYPE(VBRDF_Output_Exception_Handling)   :: VBRDF_Sup_OutputStatus


		  !  BRDF input and output structure
	!!          TYPE(BRDF_Surface_def)             :: VLIDORT_BRDF_inputs
	!!          TYPE(VBRDF_Sup_Inputs)                 :: VBRDF_Sup_In
		  TYPE(VBRDF_Sup_Outputs)                :: VBRDF_Sup_Out

		  ! SLEAVE output structure
		  TYPE(VSLEAVE_Sup_Outputs)                :: VSLEAVE_Sup_Out                      

		  !  VLIDORT supplements i/o structure

		  TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup

          !  VLIDORT Bookkeeping and FO/MS Geometry outputs
          TYPE(VLIDORT_Bookkeeping) :: VLIDORT_Bookkeep
          TYPE(VLIDORT_Geometry_FO) :: VLIDORT_FOGeom
          TYPE(VLIDORT_Geometry_MS) :: VLIDORT_MSGeom
	      END TYPE VLIDORT_IO
	      

	      TYPE VLIDORT
		logical     :: initialized = .false.  
		integer     :: NSTREAMS = 6        ! Number of half-space streams
		integer     :: NBEAMS = 1          ! Number of solar zenith angles
		integer     :: N_USER_STREAMS = 1  ! Number of Viewing zenith angles
		integer     :: N_USER_RELAZMS = 1  ! Number of relative azimuth angles
		integer     :: N_USER_LEVELS  = 1  ! Number of user-defined vertical output levels
		integer     :: N_USER_OBSGEOMS = 1 ! Number of azimuth angles calculated by surface supplement
		logical     :: DO_PLANE_PARALLEL = .false.
		type(VLIDORT_IO) :: VIO

	      END TYPE VLIDORT

	      contains

!.............................................................................

      subroutine VLIDORT_Init (self, km, rc, DO_BOA)
                                   
      USE VLIDORT_PARS_m
      USE VLIDORT_Inputs_def_m

      USE VLIDORT_INPUTS_m
     
      type(VLIDORT),              intent(inout) :: self

      integer,                    intent(in)    :: km     ! number of atmospheric layers
      integer,                    intent(out)   :: rc     ! error code

    ! !OPTIONAL PARAMETERS
      logical, optional, intent(in) :: DO_BOA      
!                     ---

      logical :: vector


!  Local variables

      logical                                :: DO_FULLRAD_MODE
      logical                                :: DO_FOCORR_NADIR
      logical                                :: DO_FOCORR_OUTGOING
      logical                                :: DO_FOCORR
      logical                                :: DO_FOCORR_EXTERNAL
      logical                                :: DO_FOCORR_DOCRIT
      real*8                                 :: DO_FOCORR_ACRIT
      logical                                :: DO_DOUBLE_CONVTEST
      logical                                :: DO_SOLAR_SOURCES
      logical                                :: DO_PLANE_PARALLEL
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
      real*8                                 :: VLIDORT_ACCURACY
      real*8                                 :: FLUX_FACTOR
      integer                                :: TF_MAXITER
      real*8                                 :: TF_CRITERION
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
      real*8                                 :: LAMBERTIAN_ALBEDO
      logical                                :: DO_THERMAL_EMISSION
      integer                                :: N_THERMAL_COEFFS
      real*8, dimension( 0:MAXLAYERS )       :: THERMAL_BB_INPUT 
      logical                                :: DO_SURFACE_EMISSION
      real*8                                 :: SURFBB
      logical                                :: DO_THERMAL_TRANSONLY

      rc = 0

     
!     initialize variables to false or zero
!     -------------------------------------    

      call VLIDORT_INIT_INPUTS(self%VIO%VLIDORT_FixIn,self%VIO%VLIDORT_ModIn)

      !    Customize those parameters that do not vary from pixel to pixel
!     ---------------------------------------------------------------

!     Increase N_USER_LEVELS if BOA Calculation is requested
      if (present(DO_BOA)) then
        if (DO_BOA) then
          self%N_USER_LEVELS = 2
        end if
      end if

!                         Modes of Operation
!                         ------------------

      DO_FULLRAD_MODE    = .true.  ! Do full Stokes vector calculation?  If false, and DO_FOCORR is true, returns single scatter
      DO_FOCORR_NADIR    = .false. ! Do nadir single scatter correction?
      DO_FOCORR_OUTGOING = .true.  ! Do outgoing single scatter correction?
      DO_FOCORR          = .true.  ! Do First-Order correction?  Must be set tu use exact single scatter instead of the truncated phase function
      DO_FOCORR_EXTERNAL = .false. ! Use First-Order results computed externally
      DO_FOCORR_DOCRIT   = .false. ! Outgoing Sphericity criticality flag (not implemented yet) 
      DO_FOCORR_ACRIT    = 0.0     ! Outgoing Sphericity attenuation criterion (not implemented yet)
!      DO_DBCORRECTION    = .true.  ! Do direct beam correction?
      DO_DOUBLE_CONVTEST = .true.  ! Perform double convergence test?
       
!                            Solar Sources
!                            -------------

      DO_SOLAR_SOURCES    = .true.     ! Include solar sources?
      DO_PLANE_PARALLEL   = self%DO_PLANE_PARALLEL    ! Plane-parallel treatment of direct beam?
      if (DO_PLANE_PARALLEL) then
          DO_FOCORR_NADIR = .false.
          DO_FOCORR_OUTGOING = .false.
      end if
      DO_CHAPMAN_FUNCTION = .true.     ! Perform internal Chapman function calculation?
      DO_REFRACTIVE_GEOMETRY = .false. ! Beam path with refractive atmosphere?
     
!                         Performance Control
!                         -------------------

      DO_RAYLEIGH_ONLY     = .false. ! Rayleigh atmosphere only?
      DO_DELTAM_SCALING    = .true.  ! Include Delta-M scaling?
      DO_SOLUTION_SAVING   = .false. ! Solution saving mode?
      DO_BVP_TELESCOPING   = .false. ! Boundary value problem telescoping mode?
      
!                      User-defined output control
!                      ---------------------------

      DO_UPWELLING = .true.     ! Upwelling output?
      DO_DNWELLING = .false.    ! Downwelling output?
      DO_USER_VZANGLES = .true. ! User-defined Viewing zenith angles?
      DO_OBSERVATION_GEOMETRY = .true. ! Do Observation Geometry?
      

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

      NSTREAMS = self%NSTREAMS         ! Number of half-space streams
      NLAYERS = km                    ! Number of atmospheric layers
      NFINELAYERS = 3                 ! Number of fine layers (outgoing sphericity correction)
      NGREEK_MOMENTS_INPUT = 300     ! Number of scattering matrix expansion coefficients
      TAYLOR_ORDER = 3                ! Number of small-number terms in Taylor series expansions
      N_USER_OBSGEOMS = self%N_USER_OBSGEOMS             ! Number of observation Geometry inputs 

!                            accuracy input
!                            --------------

      VLIDORT_ACCURACY = 0.0001  ! Fourier series convergence
      FLUX_FACTOR = 1.00d0      ! Solar flux constant; =1 if no solar sources.
      TF_MAXITER = 10 ! Maximum number of iterations in calculation of Water-leaving transmittance
      TF_CRITERION = 0.000001 ! Convergence criterion for iterative calculation of Water-leaving transmittance 
!                        Pseudo-spherical inputs
!                        -----------------------

      EARTH_RADIUS = 6371.0 ! Earth radius (km)
      RFINDEX_PARAMETER = 0.000288 ! Refractive index parameter
      GEOMETRY_SPECHEIGHT = 0.0 ! Input geometry specification height [km]


!                          Thermal controls
!                          ---------------- 

      DO_THERMAL_EMISSION  = .false.  ! Do thermal emission?
      SURFBB               = 0.0      ! surface black body emissionnn

      DO_THERMAL_TRANSONLY = .false.  ! Do thermal emission, transmittance only?
      N_THERMAL_COEFFS     = 2        ! Number of thermal coefficients
      DO_SURFACE_EMISSION  =  .false. ! Do Surface emission?



      !  normal execution: Copy all variables and return
!  -----------------------------------------------

!  Copy data to structure variables:

!  Fixed Boolean inputs

      self%VIO%VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE        = DO_FULLRAD_MODE
      self%VIO%VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION    = DO_THERMAL_EMISSION
      self%VIO%VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION    = DO_SURFACE_EMISSION
      self%VIO%VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL      = DO_PLANE_PARALLEL
      !self%VIO%VLIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE       = DO_BRDF_SURFACE
      self%VIO%VLIDORT_FixIn%Bool%TS_DO_UPWELLING           = DO_UPWELLING
      self%VIO%VLIDORT_FixIn%Bool%TS_DO_DNWELLING           = DO_DNWELLING
      !self%VIO%VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = DO_LAMBERTIAN_SURFACE


!  Modified Boolean inputs

      self%VIO%VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR        = DO_FOCORR_NADIR
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING     = DO_FOCORR_OUTGOING
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_FOCORR              = DO_FOCORR
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL     = DO_FOCORR_EXTERNAL
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST     = DO_DOUBLE_CONVTEST
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = DO_SOLAR_SOURCES
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = DO_REFRACTIVE_GEOMETRY
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = DO_CHAPMAN_FUNCTION
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY       = DO_RAYLEIGH_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY              = DO_ISOTROPIC_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH                  = DO_NO_AZIMUTH
      !VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER                 = DO_ALL_FOURIER
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING      = DO_DELTAM_SCALING
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING     = DO_SOLUTION_SAVING
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING     = DO_BVP_TELESCOPING
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY= DO_OBSERVATION_GEOMETRY
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES       = DO_USER_VZANGLES
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = DO_ADDITIONAL_MVOUT
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY          = DO_MVOUT_ONLY
      self%VIO%VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY   = DO_THERMAL_TRANSONLY

!  Fixed control inputs
!      self%VIO%VLIDORT_FixIn%Cont%TS_NSTOKES               = NSTOKES
      self%VIO%VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER           = TAYLOR_ORDER
      self%VIO%VLIDORT_FixIn%Cont%TS_NSTREAMS               = self%NSTREAMS
      self%VIO%VLIDORT_FixIn%Cont%TS_NLAYERS                = NLAYERS
      self%VIO%VLIDORT_FixIn%Cont%TS_NFINELAYERS            = NFINELAYERS
      self%VIO%VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS       = N_THERMAL_COEFFS
      self%VIO%VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY       = VLIDORT_ACCURACY
      self%VIO%VLIDORT_FixIn%Cont%TS_TF_MAXITER             = TF_MAXITER
      self%VIO%VLIDORT_FixIn%Cont%TS_TF_CRITERION           = TF_CRITERION
!  Modified control inputs

      self%VIO%VLIDORT_ModIn%Mcont%TS_NGREEK_MOMENTS_INPUT   = NGREEK_MOMENTS_INPUT
      self%VIO%VLIDORT_ModIn%Mcont%TS_GEOMETRY_SPECHEIGHT    = GEOMETRY_SPECHEIGHT
      self%VIO%VLIDORT_ModIn%MCont%TS_EARTH_RADIUS           = EARTH_RADIUS
      self%VIO%VLIDORT_ModIn%MCont%TS_DO_FOCORR_ACRIT        = DO_FOCORR_ACRIT 

!  Beam inputs

      self%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR               = FLUX_FACTOR      
      self%VIO%VLIDORT_ModIn%MSunRays%TS_N_SZANGLES               = self%NBEAMS
      self%VIO%VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS           = self%N_USER_RELAZMS     
      self%VIO%VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES          = self%N_USER_STREAMS
      self%VIO%VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS             = self%N_USER_LEVELS
      
      self%VIO%VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS          = N_USER_OBSGEOMS     
      self%VIO%VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT      = USER_OBSGEOMS
      

!  Fixed Chapman function inputs

      !VLIDORT_Chapman_inputs%TS_HEIGHT_GRID                     = HEIGHT_GRID
      !VLIDORT_Chapman_inputs%TS_PRESSURE_GRID                   = PRESSURE_GRID
      !VLIDORT_Chapman_inputs%TS_TEMPERATURE_GRID                = TEMPERATURE_GRID
      !VLIDORT_Chapman_inputs%TS_FINEGRID                        = FINEGRID
      self%VIO%VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER       = RFINDEX_PARAMETER

!  Fixed optical inputs

      !VLIDORT_Optical_inputs%TS_DELTAU_VERT_INPUT               = DELTAU_VERT_INPUT
      !VLIDORT_Optical_inputs%TS_OMEGA_TOTAL_INPUT               = OMEGA_TOTAL_INPUT
      !VLIDORT_Optical_inputs%TS_GREEKMAT_TOTAL_INPUT            = GREEKMAT_TOTAL_INPUT
      self%VIO%VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT        = THERMAL_BB_INPUT
      self%VIO%VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO       = LAMBERTIAN_ALBEDO
      self%VIO%VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT        = SURFBB

!  Fixed write inputs

      self%VIO%VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE            = DO_DEBUG_WRITE
  
      self%VIO%VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT            = DO_WRITE_INPUT
      self%VIO%VLIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME      = INPUT_WRITE_FILENAME

      self%VIO%VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO         = DO_WRITE_SCENARIO
      self%VIO%VLIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME   = SCENARIO_WRITE_FILENAME

      self%VIO%VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER          = DO_WRITE_FOURIER
      self%VIO%VLIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME    = FOURIER_WRITE_FILENAME

      self%VIO%VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS          = DO_WRITE_RESULTS
      self%VIO%VLIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME    = RESULTS_WRITE_FILENAME

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

      end subroutine VLIDORT_Init

!.............................................................................

      subroutine VLIDORT_Fin (self, rc)
      type(VLIDORT),  intent(inout) :: self
      integer,        intent(out) :: rc
      self%initialized = .false.
      
      rc = 0
      
      end subroutine VLIDORT_Fin

!.............................................................................

      end module VLIDORT_Mod
