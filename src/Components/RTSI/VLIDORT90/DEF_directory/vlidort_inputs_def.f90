
! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6                                  #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

      MODULE VLIDORT_Inputs_def

!  This module contains the following VLIDORT input structures,
!  with intents :

!        VLIDORT_Fixed_Boolean    nested in VLIDORT_Fixed_Inputs
!        VLIDORT_Fixed_Control    nested in VLIDORT_Fixed_Inputs
!        VLIDORT_Fixed_Sunrays    nested in VLIDORT_Fixed_Inputs
!     VLIDORT_Fixed_UserValues    nested in VLIDORT_Fixed_Inputs
!        VLIDORT_Fixed_Chapman    nested in VLIDORT_Fixed_Inputs
!        VLIDORT_Fixed_Optical    nested in VLIDORT_Fixed_Inputs
!          VLIDORT_Fixed_Write    nested in VLIDORT_Fixed_Inputs
!         VLIDORT_Fixed_Inputs    Intent(In)

!     VLIDORT_Modified_Boolean    nested in VLIDORT_Modified_Inputs
!     VLIDORT_Modified_Control    nested in VLIDORT_Modified_Inputs
!     VLIDORT_Modified_Sunrays    nested in VLIDORT_Modified_Inputs
!  VLIDORT_Modified_UserValues    nested in VLIDORT_Modified_Inputs
!     VLIDORT_Modified_Chapman    nested in VLIDORT_Modified_Inputs
!     VLIDORT_Modified_Optical    nested in VLIDORT_Modified_Inputs
!      VLIDORT_Modified_Inputs    Intent(InOut)

      USE VLIDORT_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Boolean


!  Full Radiance calculation

      LOGICAL     :: TS_DO_FULLRAD_MODE

!  Flag for SSCORR truncation

      LOGICAL     :: TS_DO_SSCORR_TRUNCATION

!  Flag for Use of Externally-derived single scatter results

      LOGICAL     :: TS_DO_SS_EXTERNAL

!  Flag for Full-up single scatter calculation

      LOGICAL     :: TS_DO_SSFULL

!  Surface and thermal emission flags

      LOGICAL     :: TS_DO_THERMAL_EMISSION
      LOGICAL     :: TS_DO_SURFACE_EMISSION

!  Beam particular solution pseudo-spherical options

      LOGICAL     :: TS_DO_PLANE_PARALLEL

!  Surface control (New, 23 March 2010)

!      LOGICAL     :: TS_DO_BRDF_SURFACE

!  directional control

      LOGICAL     :: TS_DO_UPWELLING
      LOGICAL     :: TS_DO_DNWELLING

!  stream angle flags

      LOGICAL     :: TS_DO_QUAD_OUTPUT

!  Contributions (RT Solutions Use Only)

      LOGICAL     :: TS_DO_TOA_CONTRIBS

!  Surface

      LOGICAL     :: TS_DO_LAMBERTIAN_SURFACE

!  Special Options (RT Solutions Use Only)

      LOGICAL     :: TS_DO_SPECIALIST_OPTION_1
      LOGICAL     :: TS_DO_SPECIALIST_OPTION_2
      LOGICAL     :: TS_DO_SPECIALIST_OPTION_3

!  surface leaving Control. New 17 May 2012

      LOGICAL     :: TS_DO_SURFACE_LEAVING
      LOGICAL     :: TS_DO_SL_ISOTROPIC

      END TYPE VLIDORT_Fixed_Boolean

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Control


!  Number of Stokes parameters

      INTEGER   :: TS_NSTOKES

!  Number of discrete ordinate streams

      INTEGER   :: TS_NSTREAMS

!  Number of computational layers

      INTEGER   :: TS_NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing spherical correction algorithm)

      INTEGER   :: TS_NFINELAYERS

!  Number of thermal coefficients (2 should be the default)

      INTEGER   :: TS_N_THERMAL_COEFFS

!  Accuracy for convergence of Fourier series

      REAL(fpk) :: TS_VLIDORT_ACCURACY

!  Special Options (RT Solutions Use Only)

      INTEGER   :: TS_NLAYERS_NOMS
      INTEGER   :: TS_NLAYERS_CUTOFF


      END TYPE VLIDORT_Fixed_Control

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Sunrays


!  Bottom-of-atmosphere solar zenith angles, DEGREES

!      REAL(fpk), dimension (MAXBEAMS) :: TS_BEAM_SZAS
      REAL(fpk), dimension (MAX_SZANGLES) :: TS_SZANGLES


      END TYPE VLIDORT_Fixed_Sunrays

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_UserValues


!  User-defined zenith angle input

!      INTEGER                                  :: TS_N_USER_STREAMS
      INTEGER                                  :: TS_N_USER_VZANGLES

!  User-defined vertical level output. From Top-of-atmosphere. Example for 4 outputs:
!     USER_LEVELS(1) = 0.0           --> Top-of-atmosphere
!     USER_LEVELS(2) = 1.1           --> One tenth of the way down into Layer 2
!     USER_LEVELS(3) = 3.5           --> One half  of the way down into Layer 4
!     USER_LEVELS(4) = dble(NLAYERS) --> Bottom of atmosphere

      INTEGER                                 :: TS_N_USER_LEVELS
      REAL(fpk), dimension (MAX_USER_LEVELS)  :: TS_USER_LEVELS


      END TYPE VLIDORT_Fixed_UserValues

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Chapman


!  multilayer Height inputsm in [km]
!   Required for the Chapman function calculations

      REAL(fpk), dimension (0:MAXLAYERS) :: TS_HEIGHT_GRID

!  multilayer atmospheric inputs. Pressures in [mb], temperatures in [K]
!   Required for the Chapman function calculations, refractive geometry

      REAL(fpk), dimension (0:MAXLAYERS) :: TS_PRESSURE_GRID
      REAL(fpk), dimension (0:MAXLAYERS) :: TS_TEMPERATURE_GRID

!  Number of fine layer gradations
!   Required for the Chapman function calculations, refractive geometry

      INTEGER,   dimension (MAXLAYERS)   :: TS_FINEGRID

!  Refractive index parameter
!    (only for Chapman function calculations with refractive index)

      REAL(fpk) :: TS_RFINDEX_PARAMETER


      END TYPE VLIDORT_Fixed_Chapman

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Optical


!  multilayer optical property (bulk)

      REAL(fpk), dimension ( MAXLAYERS ) :: TS_DELTAU_VERT_INPUT

!  Phase function Legendre-polynomial expansion coefficients
!    Include all that you require for exact single scatter calculations

      REAL(fpk), dimension ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ ) :: &
        TS_GREEKMAT_TOTAL_INPUT

!  Thermal Black Body functions

      REAL(fpk), dimension ( 0:MAXLAYERS ) :: TS_THERMAL_BB_INPUT

!  Lambertian Surface control

      REAL(fpk) :: TS_LAMBERTIAN_ALBEDO

!  Surface Black body inputs

      REAL(fpk) :: TS_SURFACE_BB_INPUT

!  Special LTE variables (RT Solutions Use Only)

      REAL(fpk), dimension ( 2, MAXLAYERS ) :: TS_LTE_DELTAU_VERT_INPUT
      REAL(fpk), dimension ( 0:MAXLAYERS )  :: TS_LTE_THERMAL_BB_INPUT


      END TYPE VLIDORT_Fixed_Optical

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Write


!  Debug control

      LOGICAL ::            TS_DO_DEBUG_WRITE

!  Input file

      LOGICAL ::            TS_DO_WRITE_INPUT
      CHARACTER (LEN=60) :: TS_INPUT_WRITE_FILENAME

!  Scene file

      LOGICAL ::            TS_DO_WRITE_SCENARIO
      CHARACTER (LEN=60) :: TS_SCENARIO_WRITE_FILENAME

!  Fourier component results file

      LOGICAL ::            TS_DO_WRITE_FOURIER
      CHARACTER (LEN=60) :: TS_FOURIER_WRITE_FILENAME

!  Results file

      LOGICAL ::            TS_DO_WRITE_RESULTS
      CHARACTER (LEN=60) :: TS_RESULTS_WRITE_FILENAME


      END TYPE VLIDORT_Fixed_Write

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Fixed_Inputs


      TYPE(VLIDORT_Fixed_Boolean)    :: Bool
      TYPE(VLIDORT_Fixed_Control)    :: Cont
      TYPE(VLIDORT_Fixed_Sunrays)    :: Sunrays
      TYPE(VLIDORT_Fixed_UserValues) :: UserVal
      TYPE(VLIDORT_Fixed_Chapman)    :: Chapman
      TYPE(VLIDORT_Fixed_Optical)    :: Optical
      TYPE(VLIDORT_Fixed_Write)      :: Write


      END TYPE VLIDORT_Fixed_Inputs

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Boolean


!  single scatter and direct beam corrections

      LOGICAL    :: TS_DO_SSCORR_NADIR
      LOGICAL    :: TS_DO_SSCORR_OUTGOING

!  double convergence test flag

      LOGICAL    :: TS_DO_DOUBLE_CONVTEST

!  Basic top-level solar beam control

      LOGICAL    :: TS_DO_SOLAR_SOURCES

!  Pseudo-spherical input control

      LOGICAL    :: TS_DO_REFRACTIVE_GEOMETRY
      LOGICAL    :: TS_DO_CHAPMAN_FUNCTION

!  scatterers and phase function control

      LOGICAL    :: TS_DO_RAYLEIGH_ONLY
!      LOGICAL    :: TS_DO_ISOTROPIC_ONLY
!      LOGICAL    :: TS_DO_NO_AZIMUTH
!      LOGICAL    :: TS_DO_ALL_FOURIER

!  Deltam scaling flag

      LOGICAL    :: TS_DO_DELTAM_SCALING

!  2 new flags in Version 3.0

      LOGICAL    :: TS_DO_SOLUTION_SAVING
      LOGICAL    :: TS_DO_BVP_TELESCOPING

!  stream angle flags

!      LOGICAL    :: TS_DO_USER_STREAMS
      LOGICAL    :: TS_DO_USER_VZANGLES

!  mean value control

      LOGICAL    :: TS_DO_ADDITIONAL_MVOUT
      LOGICAL    :: TS_DO_MVOUT_ONLY

!  Transmittance only for thermal mode.

      LOGICAL    :: TS_DO_THERMAL_TRANSONLY


      END TYPE VLIDORT_Modified_Boolean

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Control


!  number of Legendre phase function expansion moments
!       May be re-set after Checking

      INTEGER   :: TS_NGREEK_MOMENTS_INPUT


      END TYPE VLIDORT_Modified_Control

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Sunrays


!  Flux factor ( should be 1 or pi ). Same for all solar beams.

      REAL(fpk)  :: TS_FLUX_FACTOR

!  number of solar beams to be processed

!      INTEGER    :: TS_NBEAMS
      INTEGER    :: TS_N_SZANGLES


      END TYPE VLIDORT_Modified_Sunrays

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_UserValues


!  User-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER                                 :: TS_N_USER_RELAZMS
      REAL(fpk), dimension (MAX_USER_RELAZMS) :: TS_USER_RELAZMS

!  User-defined zenith angle input

!      REAL(fpk), dimension (MAX_USER_STREAMS)  :: TS_USER_ANGLES_INPUT
      REAL(fpk), dimension (MAX_USER_VZANGLES) :: TS_USER_VZANGLES_INPUT

!  Geometry specification height

      REAL(fpk)                               :: TS_GEOMETRY_SPECHEIGHT


      END TYPE VLIDORT_Modified_UserValues

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Chapman


!  Output from Chapman function calculations - can also be input (not used)

      !REAL(fpk), dimension (MAXLAYERS,MAXLAYERS,MAXBEAMS) :: &
      !  TS_CHAPMAN_FACTORS

!  Earth radius in [km] for Chapman function calculation of TAUTHICK_INPUT

      REAL(fpk) :: TS_EARTH_RADIUS


      END TYPE VLIDORT_Modified_Chapman

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Optical


!  multilayer optical property (bulk)

      REAL(fpk), dimension ( MAXLAYERS ) :: TS_OMEGA_TOTAL_INPUT


      END TYPE VLIDORT_Modified_Optical

! #####################################################################
! #####################################################################

      TYPE VLIDORT_Modified_Inputs


      TYPE(VLIDORT_Modified_Boolean)    :: MBool
      TYPE(VLIDORT_Modified_Control)    :: MCont
      TYPE(VLIDORT_Modified_Sunrays)    :: MSunrays
      TYPE(VLIDORT_Modified_UserValues) :: MUserVal
      TYPE(VLIDORT_Modified_Chapman)    :: MChapman
      TYPE(VLIDORT_Modified_Optical)    :: MOptical


      END TYPE VLIDORT_Modified_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VLIDORT_Fixed_Boolean, &
                VLIDORT_Fixed_Control, &
                VLIDORT_Fixed_Sunrays, &
                VLIDORT_Fixed_UserValues, &
                VLIDORT_Fixed_Chapman, &
                VLIDORT_Fixed_Optical, &
                VLIDORT_Fixed_Write, &
                VLIDORT_Fixed_Inputs, &
                VLIDORT_Modified_Boolean, &
                VLIDORT_Modified_Control, &
                VLIDORT_Modified_Sunrays, &
                VLIDORT_Modified_UserValues, &
                VLIDORT_Modified_Chapman, &
                VLIDORT_Modified_Optical, &
                VLIDORT_Modified_Inputs

      END MODULE VLIDORT_Inputs_def
