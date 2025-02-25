!
! CRTM_RTSolution_Define
!
! Module defining the CRTM RTSolution structure and containing routines
! to manipulate it.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 13-May-2004
!                       paul.vandelst@noaa.gov
!

MODULE CRTM_RTSolution_Define


  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use statements
  USE Type_Kinds           , ONLY: fp
  USE Message_Handler      , ONLY: SUCCESS, FAILURE, WARNING, INFORMATION, Display_Message
  USE Compare_Float_Numbers, ONLY: DEFAULT_N_SIGFIG, &
                                   OPERATOR(.EqualTo.), &
                                   Compares_Within_Tolerance
  USE File_Utility         , ONLY: File_Open, File_Exists
  USE Binary_File_Utility  , ONLY: Open_Binary_File      , &
                                   WriteGAtts_Binary_File, &
                                   ReadGAtts_Binary_File
  USE SensorInfo_Parameters, ONLY: INVALID_SENSOR, &
                                   INVALID_WMO_SATELLITE_ID, &
                                   INVALID_WMO_SENSOR_ID
  USE CRTM_Parameters      , ONLY: STRLEN
  ! Disable all implicit typing
  IMPLICIT NONE


  ! --------------------
  ! Default visibilities
  ! --------------------
  ! Everything private by default
  PRIVATE
  ! Datatypes
  PUBLIC :: CRTM_RTSolution_type
  ! Operators
  PUBLIC :: OPERATOR(==)
  PUBLIC :: OPERATOR(-)
  ! Public procedures
  PUBLIC :: CRTM_RTSolution_Associated
  PUBLIC :: CRTM_RTSolution_Destroy
  PUBLIC :: CRTM_RTSolution_Create
  PUBLIC :: CRTM_RTSolution_Zero
  PUBLIC :: CRTM_RTSolution_Inspect
  PUBLIC :: CRTM_RTSolution_DefineVersion
  PUBLIC :: CRTM_RTSolution_Compare
  PUBLIC :: CRTM_RTSolution_InquireFile
  PUBLIC :: CRTM_RTSolution_ReadFile
  PUBLIC :: CRTM_RTSolution_WriteFile


  ! ---------------------
  ! Procedure overloading
  ! ---------------------
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE CRTM_RTSolution_Equal
  END INTERFACE OPERATOR(==)

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE CRTM_RTSolution_Subtract
  END INTERFACE OPERATOR(-)

  INTERFACE CRTM_RTSolution_Inspect
    MODULE PROCEDURE Scalar_Inspect
    MODULE PROCEDURE Rank2_Inspect
  END INTERFACE CRTM_RTSolution_Inspect


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id$'
  ! Literal constants
  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  ! Message string length
  INTEGER, PARAMETER :: ML = 256
  ! File status on close after write error
  CHARACTER(*), PARAMETER :: WRITE_ERROR_STATUS = 'DELETE'


  ! -------------------------------
  ! RTSolution data type definition
  ! -------------------------------
  !:tdoc+:
  TYPE :: CRTM_RTSolution_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Dimensions
    INTEGER :: n_Layers = 0  ! K
    ! Sensor information
    CHARACTER(STRLEN) :: Sensor_ID        = ''
    INTEGER           :: WMO_Satellite_ID = INVALID_WMO_SATELLITE_ID
    INTEGER           :: WMO_Sensor_ID    = INVALID_WMO_SENSOR_ID
    INTEGER           :: Sensor_Channel   = 0
    ! RT algorithm information
    CHARACTER(STRLEN) :: RT_Algorithm_Name = ''
    ! Internal variables. Users do not need to worry about these.
    LOGICAL :: Scattering_Flag = .TRUE.
    INTEGER :: n_Full_Streams  = 0
    INTEGER :: n_Stokes        = 0
    ! Forward radiative transfer intermediate results for a single channel
    !    These components are not defined when they are used as TL, AD
    !    and K variables
    REAL(fp) :: SOD                     = ZERO  ! Scattering Optical Depth
    REAL(fp) :: Surface_Emissivity      = ZERO
    REAL(fp) :: Up_Radiance             = ZERO
    REAL(fp) :: Down_Radiance           = ZERO
    REAL(fp) :: Down_Solar_Radiance     = ZERO
    REAL(fp) :: Surface_Planck_Radiance = ZERO
    REAL(fp), ALLOCATABLE :: Upwelling_Radiance(:)   ! K
    REAL(fp), ALLOCATABLE :: Layer_Optical_Depth(:)  ! K
    ! Radiative transfer results for a single channel/node
    REAL(fp) :: Radiance               = ZERO
    REAL(fp) :: Brightness_Temperature = ZERO
  END TYPE CRTM_RTSolution_type
  !:tdoc-:


CONTAINS


!##################################################################################
!##################################################################################
!##                                                                              ##
!##                           ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_RTSolution_Associated
!
! PURPOSE:
!       Elemental function to test the status of the allocatable components
!       of a CRTM RTSolution object.
!
! CALLING SEQUENCE:
!       Status = CRTM_RTSolution_Associated( RTSolution )
!
! OBJECTS:
!       RTSolution:   RTSolution structure which is to have its member's
!                     status tested.
!                     UNITS:      N/A
!                     TYPE:       CRTM_RTSolution_type
!                     DIMENSION:  Scalar or any rank
!                     ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       Status:       The return value is a logical value indicating the
!                     status of the RTSolution members.
!                       .TRUE.  - if the array components are allocated.
!                       .FALSE. - if the array components are not allocated.
!                     UNITS:      N/A
!                     TYPE:       LOGICAL
!                     DIMENSION:  Same as input RTSolution argument
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL FUNCTION CRTM_RTSolution_Associated( RTSolution ) RESULT( Status )
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: RTSolution
    LOGICAL :: Status
    Status = RTSolution%Is_Allocated
  END FUNCTION CRTM_RTSolution_Associated


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_RTSolution_Destroy
!
! PURPOSE:
!       Elemental subroutine to re-initialize CRTM RTSolution objects.
!
! CALLING SEQUENCE:
!       CALL CRTM_RTSolution_Destroy( RTSolution )
!
! OBJECTS:
!       RTSolution:   Re-initialized RTSolution structure.
!                     UNITS:      N/A
!                     TYPE:       CRTM_RTSolution_type
!                     DIMENSION:  Scalar OR any rank
!                     ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE CRTM_RTSolution_Destroy( RTSolution )
    TYPE(CRTM_RTSolution_type), INTENT(OUT) :: RTSolution
    RTSolution%Is_Allocated = .FALSE.
    RTSolution%n_Layers = 0
  END SUBROUTINE CRTM_RTSolution_Destroy


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_RTSolution_Create
!
! PURPOSE:
!       Elemental subroutine to create an instance of the CRTM RTSolution object.
!
! CALLING SEQUENCE:
!       CALL CRTM_RTSolution_Create( RTSolution, n_Layers )
!
! OBJECTS:
!       RTSolution:   RTSolution structure.
!                     UNITS:      N/A
!                     TYPE:       CRTM_RTSolution_type
!                     DIMENSION:  Scalar or any rank
!                     ATTRIBUTES: INTENT(OUT)
!
! INPUTS:
!       n_Layers:     Number of layers for which there is RTSolution data.
!                     Must be > 0.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Same as RTSolution object
!                     ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE CRTM_RTSolution_Create( RTSolution, n_Layers )
    ! Arguments
    TYPE(CRTM_RTSolution_type), INTENT(OUT) :: RTSolution
    INTEGER,                    INTENT(IN)  :: n_Layers
    ! Local variables
    INTEGER :: alloc_stat

    ! Check input
    IF ( n_Layers < 1 ) RETURN

    ! Perform the allocation
    ALLOCATE( RTSolution%Upwelling_Radiance(n_Layers), &
              RTSolution%Layer_Optical_Depth(n_Layers), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN

    ! Initialise
    ! ...Dimensions
    RTSolution%n_Layers = n_Layers
    ! ...Arrays
    RTSolution%Upwelling_Radiance  = ZERO
    RTSolution%Layer_Optical_Depth = ZERO

    ! Set allocation indicator
    RTSolution%Is_Allocated = .TRUE.

  END SUBROUTINE CRTM_RTSolution_Create


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_RTSolution_Zero
!
! PURPOSE:
!       Elemental subroutine to zero out the data components
!       in a CRTM RTSolution object.
!
! CALLING SEQUENCE:
!       CALL CRTM_RTSolution_Zero( rts )
!
! OUTPUTS:
!       rts:          CRTM RTSolution structure in which the data components
!                     are to be zeroed out.
!                     UNITS:      N/A
!                     TYPE:       CRTM_RTSolution_type
!                     DIMENSION:  Scalar or any rank
!                     ATTRIBUTES: INTENT(IN OUT)
!
! COMMENTS:
!       - The dimension components of the structure are *NOT* set to zero.
!       - The sensor infomration and RT algorithm components are
!         *NOT* reset in this routine.
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE CRTM_RTSolution_Zero( RTSolution )
    TYPE(CRTM_RTSolution_type), INTENT(IN OUT) :: RTSolution

    ! Zero out the scalar data components
    RTSolution%SOD                     = ZERO
    RTSolution%Surface_Emissivity      = ZERO
    RTSolution%Up_Radiance             = ZERO
    RTSolution%Down_Radiance           = ZERO
    RTSolution%Down_Solar_Radiance     = ZERO
    RTSolution%Surface_Planck_Radiance = ZERO
    RTSolution%Radiance                = ZERO
    RTSolution%Brightness_Temperature  = ZERO

    ! Zero out the array data components
    IF ( CRTM_RTSolution_Associated(RTSolution) ) THEN
      RTSolution%Upwelling_Radiance  = ZERO
      RTSolution%Layer_Optical_Depth = ZERO
    END IF

  END SUBROUTINE CRTM_RTSolution_Zero


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_RTSolution_Inspect
!
! PURPOSE:
!       Subroutine to print the contents of a CRTM RTSolution object to stdout.
!
! CALLING SEQUENCE:
!       CALL CRTM_RTSolution_Inspect( RTSolution )
!
! INPUTS:
!       RTSolution:    CRTM RTSolution object to display.
!                      UNITS:      N/A
!                      TYPE:       CRTM_RTSolution_type
!                      DIMENSION:  Scalar or Rank-2 (n_channels x n_profiles)
!                      ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE Scalar_Inspect( RTSolution )
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: RTSolution
    WRITE(*,'(1x,"RTSolution OBJECT")')
    ! Display components
    WRITE(*,'(3x,"Sensor Id                : ",a )') TRIM(RTSolution%Sensor_ID)
    WRITE(*,'(3x,"WMO Satellite Id         : ",i0)') RTSolution%WMO_Satellite_ID
    WRITE(*,'(3x,"WMO Sensor Id            : ",i0)') RTSolution%WMO_Sensor_ID
    WRITE(*,'(3x,"Channel                  : ",i0)') RTSolution%Sensor_Channel
    WRITE(*,'(3x,"RT Algorithm Name        : ",a )') RTSolution%RT_Algorithm_Name
    WRITE(*,'(3x,"Scattering Optical Depth : ",es13.6)') RTSolution%SOD
    WRITE(*,'(3x,"Surface Emissivity       : ",es13.6)') RTSolution%Surface_Emissivity
    WRITE(*,'(3x,"Up Radiance              : ",es13.6)') RTSolution%Up_Radiance
    WRITE(*,'(3x,"Down Radiance            : ",es13.6)') RTSolution%Down_Radiance
    WRITE(*,'(3x,"Down Solar Radiance      : ",es13.6)') RTSolution%Down_Solar_Radiance
    WRITE(*,'(3x,"Surface Planck Radiance  : ",es13.6)') RTSolution%Surface_Planck_Radiance
    IF ( CRTM_RTSolution_Associated(RTSolution) ) THEN
      WRITE(*,'(3x,"n_Layers : ",i0)') RTSolution%n_Layers
      WRITE(*,'(3x,"Upwelling Radiance       :")')
      WRITE(*,'(5(1x,es13.6,:))') RTSolution%Upwelling_Radiance
      WRITE(*,'(3x,"Layer Optical Depth      :")')
      WRITE(*,'(5(1x,es13.6,:))') RTSolution%Layer_Optical_Depth
    END IF
    WRITE(*,'(3x,"Radiance                 : ",es13.6)') RTSolution%Radiance
    WRITE(*,'(3x,"Brightness Temperature   : ",es13.6)') RTSolution%Brightness_Temperature
  END SUBROUTINE Scalar_Inspect


  SUBROUTINE Rank2_Inspect( RTSolution )
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: RTSolution(:,:)
    INTEGER :: i, n_channels
    INTEGER :: j, n_profiles

    n_channels = SIZE(RTSolution,1)
    n_profiles = SIZE(RTSolution,2)
    DO j = 1, n_profiles
      DO i = 1, n_channels
        WRITE(*, FMT='(1x,"PROFILE INDEX:",i0,", CHANNEL INDEX:",i0," - ")', ADVANCE='NO') j,i
        CALL Scalar_Inspect(RTSolution(i,j))
      END DO
    END DO
  END SUBROUTINE Rank2_Inspect



!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_RTSolution_DefineVersion
!
! PURPOSE:
!       Subroutine to return the module version information.
!
! CALLING SEQUENCE:
!       CALL CRTM_RTSolution_DefineVersion( Id )
!
! OUTPUTS:
!       Id:            Character string containing the version Id information
!                      for the module.
!                      UNITS:      N/A
!                      TYPE:       CHARACTER(*)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_RTSolution_DefineVersion( Id )
    CHARACTER(*), INTENT(OUT) :: Id
    Id = MODULE_VERSION_ID
  END SUBROUTINE CRTM_RTSolution_DefineVersion

!------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       CRTM_RTSolution_Compare
!
! PURPOSE:
!       Elemental function to compare two CRTM_RTSolution objects to within
!       a user specified number of significant figures.
!
! CALLING SEQUENCE:
!       is_comparable = CRTM_RTSolution_Compare( x, y, n_SigFig=n_SigFig )
!
! OBJECTS:
!       x, y:          Two CRTM RTSolution objects to be compared.
!                      UNITS:      N/A
!                      TYPE:       CRTM_RTSolution_type
!                      DIMENSION:  Scalar or any rank
!                      ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUTS:
!       n_SigFig:      Number of significant figure to compare floating point
!                      components.
!                      UNITS:      N/A
!                      TYPE:       INTEGER
!                      DIMENSION:  Conformable with inputs
!                      ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       is_comparable: Logical value indicating whether the inputs are
!                      comparable.
!                      UNITS:      N/A
!                      TYPE:       LOGICAL
!                      DIMENSION:  Same as inputs.
!:sdoc-:
!------------------------------------------------------------------------------

  ELEMENTAL FUNCTION CRTM_RTSolution_Compare( &
    x, &
    y, &
    n_SigFig ) &
  RESULT( is_comparable )
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: x, y
    INTEGER,          OPTIONAL, INTENT(IN) :: n_SigFig
    LOGICAL :: is_comparable
    ! Variables
    INTEGER :: n

    ! Set up
    is_comparable = .FALSE.
    IF ( PRESENT(n_SigFig) ) THEN
      n = ABS(n_SigFig)
    ELSE
      n = DEFAULT_N_SIGFIG
    END IF

    ! Check the structure association status
    IF ( CRTM_RTSolution_Associated(x) .NEQV. CRTM_RTSolution_Associated(y) ) RETURN

    ! Check the sensor information
    IF ( (x%Sensor_ID        /= y%Sensor_ID       ) .OR. &
         (x%WMO_Satellite_ID /= y%WMO_Satellite_ID) .OR. &
         (x%WMO_Sensor_ID    /= y%WMO_Sensor_ID   ) .OR. &
         (x%Sensor_Channel   /= y%Sensor_Channel  ) ) RETURN

    ! Check the RT algorithm name
    IF ( x%RT_Algorithm_Name /= y%RT_Algorithm_Name ) RETURN

    ! Check the scalar components
    IF ( .NOT. Compares_Within_Tolerance(x%SOD                    , y%SOD                    , n) .OR. &
         .NOT. Compares_Within_Tolerance(x%Surface_Emissivity     , y%Surface_Emissivity     , n) .OR. &
         .NOT. Compares_Within_Tolerance(x%Up_Radiance            , y%Up_Radiance            , n) .OR. &
         .NOT. Compares_Within_Tolerance(x%Down_Radiance          , y%Down_Radiance          , n) .OR. &
         .NOT. Compares_Within_Tolerance(x%Down_Solar_Radiance    , y%Down_Solar_Radiance    , n) .OR. &
         .NOT. Compares_Within_Tolerance(x%Surface_Planck_Radiance, y%Surface_Planck_Radiance, n) .OR. &
         .NOT. Compares_Within_Tolerance(x%Radiance               , y%Radiance               , n) .OR. &
         .NOT. Compares_Within_Tolerance(x%Brightness_Temperature , y%Brightness_Temperature , n) ) RETURN

    ! Check the array components
    IF ( CRTM_RTSolution_Associated(x) .AND. CRTM_RTSolution_Associated(y) ) THEN
      IF ( (.NOT. ALL(Compares_Within_Tolerance(x%Upwelling_Radiance ,y%Upwelling_Radiance ,n))) .OR. &
           (.NOT. ALL(Compares_Within_Tolerance(x%Layer_Optical_Depth,y%Layer_Optical_Depth,n))) ) RETURN
    END IF

    ! If we get here, the structures are comparable
    is_comparable = .TRUE.

  END FUNCTION CRTM_RTSolution_Compare


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_RTSolution_InquireFile
!
! PURPOSE:
!       Function to inquire CRTM RTSolution object files.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_RTSolution_InquireFile( Filename               , &
!                                                   n_Channels = n_Channels, &
!                                                   n_Profiles = n_Profiles  )
!
! INPUTS:
!       Filename:       Character string specifying the name of a
!                       CRTM RTSolution data file to read.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(*)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OPTIONAL OUTPUTS:
!       n_Channels:     The number of spectral channels for which there is
!                       data in the file.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: OPTIONAL, INTENT(OUT)
!
!       n_Profiles:     The number of profiles in the data file.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: OPTIONAL, INTENT(OUT)
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS, the file inquire was successful
!                          == FAILURE, an unrecoverable error occurred.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_RTSolution_InquireFile( &
    Filename   , &  ! Input
    n_Channels , &  ! Optional output
    n_Profiles ) &  ! Optional output
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN)  :: Filename
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Channels
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Profiles
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_RTSolution_InquireFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat
    INTEGER :: fid
    INTEGER :: l, m

    ! Set up
    err_stat = SUCCESS
    ! Check that the file exists
    IF ( .NOT. File_Exists( TRIM(Filename) ) ) THEN
      msg = 'File '//TRIM(Filename)//' not found.'
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Open the file
    err_stat = Open_Binary_File( Filename, fid )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error opening '//TRIM(Filename)
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Read the number of channels,profiles
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) l, m
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading dimensions from '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Close the file
    CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
    IF ( io_stat /= 0 ) THEN
      msg = 'Error closing '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Set the return arguments
    IF ( PRESENT(n_Channels) ) n_Channels = l
    IF ( PRESENT(n_Profiles) ) n_Profiles = m

  CONTAINS

    SUBROUTINE Inquire_CleanUp()
      IF ( File_Open( Filename ) ) THEN
        CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
        IF ( io_stat /= SUCCESS ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup - '//TRIM(io_msg)
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Inquire_CleanUp

  END FUNCTION CRTM_RTSolution_InquireFile


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_RTSolution_ReadFile
!
! PURPOSE:
!       Function to read CRTM RTSolution object files.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_RTSolution_ReadFile( Filename                , &
!                                                RTSolution              , &
!                                                Quiet      = Quiet      , &
!                                                n_Channels = n_Channels , &
!                                                n_Profiles = n_Profiles , &
!
! INPUTS:
!       Filename:     Character string specifying the name of an
!                     RTSolution format data file to read.
!                     UNITS:      N/A
!                     TYPE:       CHARACTER(*)
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       RTSolution:   CRTM RTSolution object array containing the RTSolution
!                     data.
!                     UNITS:      N/A
!                     TYPE:       CRTM_RTSolution_type
!                     DIMENSION:  Rank-2 (n_Channels x n_Profiles)
!                     ATTRIBUTES: INTENT(OUT)
!
! OPTIONAL INPUTS:
!       Quiet:        Set this logical argument to suppress INFORMATION
!                     messages being printed to stdout
!                     If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                        == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                     If not specified, default is .FALSE.
!                     UNITS:      N/A
!                     TYPE:       LOGICAL
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(IN), OPTIONAL
!
! OPTIONAL OUTPUTS:
!       n_Channels:   The number of channels for which data was read.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: OPTIONAL, INTENT(OUT)
!
!       n_Profiles:   The number of profiles for which data was read.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: OPTIONAL, INTENT(OUT)
!
!
! FUNCTION RESULT:
!       Error_Status: The return value is an integer defining the error status.
!                     The error codes are defined in the Message_Handler module.
!                     If == SUCCESS, the file read was successful
!                        == FAILURE, an unrecoverable error occurred.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_RTSolution_ReadFile( &
    Filename   , &  ! Input
    RTSolution , &  ! Output
    Quiet      , &  ! Optional input
    n_Channels , &  ! Optional output
    n_Profiles , &  ! Optional output
    Debug      ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),               INTENT(IN)  :: Filename
    TYPE(CRTM_RTSolution_type), INTENT(OUT) :: RTSolution(:,:)
    LOGICAL,          OPTIONAL, INTENT(IN)  :: Quiet
    INTEGER,          OPTIONAL, INTENT(OUT) :: n_Channels
    INTEGER,          OPTIONAL, INTENT(OUT) :: n_Profiles
    LOGICAL,          OPTIONAL, INTENT(IN)  :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_RTSolution_ReadFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat
    LOGICAL :: noisy
    INTEGER :: fid
    INTEGER :: l, n_file_channels, n_input_channels
    INTEGER :: m, n_file_profiles, n_input_profiles


    ! Set up
    err_stat = SUCCESS
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Override Quiet settings if debug set.
    IF ( PRESENT(Debug) ) noisy = Debug
    ! ...Check that the file exists
    IF ( .NOT. File_Exists( TRIM(Filename) ) ) THEN
      msg = 'File '//TRIM(Filename)//' not found.'
      CALL Read_Cleanup(); RETURN
    END IF


    ! Open the file
    err_stat = Open_Binary_File( Filename, fid )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error opening '//TRIM(Filename)
      CALL Read_Cleanup(); RETURN
    END IF


    ! Read the dimensions
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) n_file_channels, n_file_profiles
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading dimensions from '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Check if n_Channels in file is > size of output array
    n_input_channels = SIZE(RTSolution,DIM=1)
    IF ( n_file_channels > n_input_channels ) THEN
      WRITE( msg,'("Number of channels, ",i0," > size of the output RTSolution", &
                  &" array dimension, ",i0,". Only the first ",i0, &
                  &" channels will be read.")' ) &
                  n_file_channels, n_input_channels, n_input_channels
      CALL Display_Message( ROUTINE_NAME, msg, WARNING )
    END IF
    n_input_channels = MIN(n_input_channels, n_file_channels)
    ! ...Check if n_Profiles in file is > size of output array
    n_input_profiles = SIZE(RTSolution,DIM=2)
    IF ( n_file_profiles > n_input_profiles ) THEN
      WRITE( msg,'( "Number of profiles, ",i0," > size of the output RTSolution", &
                   &" array dimension, ",i0,". Only the first ",i0, &
                   &" profiles will be read.")' ) &
                   n_file_profiles, n_input_profiles, n_input_profiles
      CALL Display_Message( ROUTINE_NAME, msg, WARNING )
    END IF
    n_input_profiles = MIN(n_input_profiles, n_file_profiles)


    ! Loop over all the profiles and channels
    Profile_Loop: DO m = 1, n_input_profiles
      Channel_Loop: DO l = 1, n_input_channels
        err_stat = Read_Record( fid, RTSolution(l,m) )
        IF ( err_stat /= SUCCESS ) THEN
          WRITE( msg,'("Error reading RTSolution element (",i0,",",i0,") from ",a)' ) &
                 l, m, TRIM(Filename)
          CALL Read_Cleanup(); RETURN
        END IF
      END DO Channel_Loop
    END DO Profile_Loop


    ! Close the file
    CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
    IF ( io_stat /= 0 ) THEN
      msg = 'Error closing '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF


    ! Set the return values
    IF ( PRESENT(n_Channels) ) n_Channels = n_input_channels
    IF ( PRESENT(n_Profiles) ) n_Profiles = n_input_profiles


    ! Output an info message
    IF ( noisy ) THEN
      WRITE( msg,'("Number of channels and profiles read from ",a,": ",i0,1x,i0)' ) &
             TRIM(Filename), n_input_channels, n_input_profiles
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Read_CleanUp()
      IF ( File_Open( Filename ) ) THEN
        CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
        IF ( io_stat /= 0 ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup - '//TRIM(io_msg)
      END IF
      CALL CRTM_RTSolution_Destroy( RTSolution )
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Read_CleanUp

  END FUNCTION CRTM_RTSolution_ReadFile


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_RTSolution_WriteFile
!
! PURPOSE:
!       Function to write CRTM RTSolution object files.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_RTSolution_WriteFile( Filename     , &
!                                                 RTSolution   , &
!                                                 Quiet = Quiet  )
!
! INPUTS:
!       Filename:     Character string specifying the name of the
!                     RTSolution format data file to write.
!                     UNITS:      N/A
!                     TYPE:       CHARACTER(*)
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(IN)
!
!       RTSolution:   CRTM RTSolution object array containing the RTSolution
!                     data.
!                     UNITS:      N/A
!                     TYPE:       CRTM_RTSolution_type
!                     DIMENSION:  Rank-2 (n_Channels x n_Profiles)
!                     ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUTS:
!       Quiet:        Set this logical argument to suppress INFORMATION
!                     messages being printed to stdout
!                     If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                        == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                     If not specified, default is .FALSE.
!                     UNITS:      N/A
!                     TYPE:       LOGICAL
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status: The return value is an integer defining the error status.
!                     The error codes are defined in the Message_Handler module.
!                     If == SUCCESS, the file write was successful
!                        == FAILURE, an unrecoverable error occurred.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       - If the output file already exists, it is overwritten.
!       - If an error occurs during *writing*, the output file is deleted before
!         returning to the calling routine.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_RTSolution_WriteFile( &
    Filename   , &  ! Input
    RTSolution , &  ! Input
    Quiet      , &  ! Optional input
    Debug      ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),               INTENT(IN) :: Filename
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: RTSolution(:,:)
    LOGICAL,          OPTIONAL, INTENT(IN) :: Quiet
    LOGICAL,          OPTIONAL, INTENT(IN) :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_RTSolution_WriteFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat
    LOGICAL :: noisy
    INTEGER :: fid
    INTEGER :: l, n_output_channels
    INTEGER :: m, n_output_profiles

    ! Set up
    err_stat = SUCCESS
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Override Quiet settings if debug set.
    IF ( PRESENT(Debug) ) noisy = Debug
    n_output_channels = SIZE(RTSolution,DIM=1)
    n_output_profiles = SIZE(RTSolution,DIM=2)


    ! Open the file
    err_stat = Open_Binary_File( Filename, fid, For_Output = .TRUE. )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error opening '//TRIM(Filename)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Write the dimensions
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) n_output_channels, n_output_profiles
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing dimensions to '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Write the data
    Profile_Loop: DO m = 1, n_output_profiles
      Channel_Loop: DO l = 1, n_output_channels
        err_stat = Write_Record( fid, RTSolution(l,m) )
        IF ( err_stat /= SUCCESS ) THEN
          WRITE( msg,'("Error writing RTSolution element (",i0,",",i0,") to ",a)' ) &
                 l, m, TRIM(Filename)
          CALL Write_Cleanup(); RETURN
        END IF
      END DO Channel_Loop
    END DO Profile_Loop


    ! Close the file (if error, no delete)
    CLOSE( fid,STATUS='KEEP',IOSTAT=io_stat,IOMSG=io_msg )
    IF ( io_stat /= 0 ) THEN
      msg = 'Error closing '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Output an info message
    IF ( noisy ) THEN
      WRITE( msg,'("Number of channels and profiles written to ",a,": ",i0,1x,i0 )' ) &
             TRIM(Filename), n_output_channels, n_output_profiles
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Write_CleanUp()
      IF ( File_Open( Filename ) ) THEN
        CLOSE( fid,STATUS=WRITE_ERROR_STATUS,IOSTAT=io_stat,IOMSG=io_msg )
        IF ( io_stat /= 0 ) &
          msg = TRIM(msg)//'; Error deleting output file during error cleanup - '//TRIM(io_msg)
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Write_CleanUp

  END FUNCTION CRTM_RTSolution_WriteFile



!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################

!------------------------------------------------------------------------------
!
! NAME:
!       CRTM_RTSolution_Equal
!
! PURPOSE:
!       Elemental function to test the equality of two CRTM_RTSolution objects.
!       Used in OPERATOR(==) interface block.
!
! CALLING SEQUENCE:
!       is_equal = CRTM_RTSolution_Equal( x, y )
!
!         or
!
!       IF ( x == y ) THEN
!         ...
!       END IF
!
! OBJECTS:
!       x, y:          Two CRTM RTSolution objects to be compared.
!                      UNITS:      N/A
!                      TYPE:       CRTM_RTSolution_type
!                      DIMENSION:  Scalar or any rank
!                      ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       is_equal:      Logical value indicating whether the inputs are equal.
!                      UNITS:      N/A
!                      TYPE:       LOGICAL
!                      DIMENSION:  Same as inputs.
!
!------------------------------------------------------------------------------

  ELEMENTAL FUNCTION CRTM_RTSolution_Equal( x, y ) RESULT( is_equal )
    TYPE(CRTM_RTSolution_type) , INTENT(IN)  :: x, y
    LOGICAL :: is_equal

    ! Setup
    is_equal = .FALSE.

    ! Check the structure association status
    IF ( CRTM_RTSolution_Associated(x) .NEQV. CRTM_RTSolution_Associated(y) ) RETURN

    ! Check scalars
    IF ( (x%n_Layers == y%n_Layers) .AND. &
         (x%Sensor_ID         == y%Sensor_ID        ) .AND. &
         (x%WMO_Satellite_ID  == y%WMO_Satellite_ID ) .AND. &
         (x%WMO_Sensor_ID     == y%WMO_Sensor_ID    ) .AND. &
         (x%Sensor_Channel    == y%Sensor_Channel   ) .AND. &
         (x%RT_Algorithm_Name == y%RT_Algorithm_Name) .AND. &
         (x%SOD                     .EqualTo. y%SOD                    ) .AND. &
         (x%Surface_Emissivity      .EqualTo. y%Surface_Emissivity     ) .AND. &
         (x%Up_Radiance             .EqualTo. y%Up_Radiance            ) .AND. &
         (x%Down_Radiance           .EqualTo. y%Down_Radiance          ) .AND. &
         (x%Down_Solar_Radiance     .EqualTo. y%Down_Solar_Radiance    ) .AND. &
         (x%Surface_Planck_Radiance .EqualTo. y%Surface_Planck_Radiance) .AND. &
         (x%Radiance                .EqualTo. y%Radiance               ) .AND. &
         (x%Brightness_Temperature  .EqualTo. y%Brightness_Temperature ) ) &
      is_equal = .TRUE.

    ! Check arrays (which may or may not be allocated)
    IF ( CRTM_RTSolution_Associated(x) .AND. CRTM_RTSolution_Associated(y) ) THEN
      is_equal = is_equal .AND. &
                 ALL(x%Upwelling_Radiance  .EqualTo. y%Upwelling_Radiance ) .AND. &
                 ALL(x%Layer_Optical_Depth .EqualTo. y%Layer_Optical_Depth)
    END IF

  END FUNCTION CRTM_RTSolution_Equal


!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_RTSolution_Subtract
!
! PURPOSE:
!       Pure function to subtract two CRTM RTSolution objects.
!       Used in OPERATOR(-) interface block.
!
! CALLING SEQUENCE:
!       rtsdiff = CRTM_RTSolution_Subtract( rts1, rts2 )
!
!         or
!
!       rtsdiff = rts1 - rts2
!
!
! INPUTS:
!       rts1, rts2: The RTSolution objects to difference.
!                   UNITS:      N/A
!                   TYPE:       CRTM_RTSolution_type
!                   DIMENSION:  Scalar
!                   ATTRIBUTES: INTENT(IN OUT)
!
! RESULT:
!       rtsdiff:    RTSolution object containing the differenced components.
!                   UNITS:      N/A
!                   TYPE:       CRTM_RTSolution_type
!                   DIMENSION:  Scalar
!
!--------------------------------------------------------------------------------

  ELEMENTAL FUNCTION CRTM_RTSolution_Subtract( rts1, rts2 ) RESULT( rtsdiff )
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: rts1, rts2
    TYPE(CRTM_RTSolution_type) :: rtsdiff
    INTEGER :: k

    ! Check input
    ! ...If input structure association status different, do nothing
    IF ( CRTM_RTSolution_Associated(rts1) .NEQV. CRTM_RTSolution_Associated(rts2) ) RETURN
    ! ...If input structure for different sensors, do nothing
    IF ( (rts1%Sensor_ID         /= rts2%Sensor_ID        ) .AND. &
         (rts1%WMO_Satellite_ID  /= rts2%WMO_Satellite_ID ) .AND. &
         (rts1%WMO_Sensor_ID     /= rts2%WMO_Sensor_ID    ) .AND. &
         (rts1%Sensor_Channel    /= rts2%Sensor_Channel   ) ) RETURN

    ! Copy the first structure
    rtsdiff = rts1

    ! And subtract the second one's components from it
    ! ...Handle RT_Algorithm_Name
    rtsdiff%RT_Algorithm_Name = TRIM(rtsdiff%RT_Algorithm_Name)//' - '//TRIM(rts2%RT_Algorithm_Name)
    ! ...The scalar values
    rtsdiff%SOD                     = rtsdiff%SOD                     - rts2%SOD
    rtsdiff%Surface_Emissivity      = rtsdiff%Surface_Emissivity      - rts2%Surface_Emissivity
    rtsdiff%Up_Radiance             = rtsdiff%Up_Radiance             - rts2%Up_Radiance
    rtsdiff%Down_Radiance           = rtsdiff%Down_Radiance           - rts2%Down_Radiance
    rtsdiff%Down_Solar_Radiance     = rtsdiff%Down_Solar_Radiance     - rts2%Down_Solar_Radiance
    rtsdiff%Surface_Planck_Radiance = rtsdiff%Surface_Planck_Radiance - rts2%Surface_Planck_Radiance
    rtsdiff%Radiance                = rtsdiff%Radiance                - rts2%Radiance
    rtsdiff%Brightness_Temperature  = rtsdiff%Brightness_Temperature  - rts2%Brightness_Temperature
    ! ...The arrays (which may or may not be allocated)
    IF ( CRTM_RTSolution_Associated(rts1) .AND. CRTM_RTSolution_Associated(rts2) ) THEN
      k = rts1%n_Layers
      rtsdiff%Upwelling_Radiance(1:k)  = rtsdiff%Upwelling_Radiance(1:k)  - rts2%Upwelling_Radiance(1:k)
      rtsdiff%Layer_Optical_Depth(1:k) = rtsdiff%Layer_Optical_Depth(1:k) - rts2%Layer_Optical_Depth(1:k)
    END IF

  END FUNCTION CRTM_RTSolution_Subtract


!
! NAME:
!       Read_Record
!
! PURPOSE:
!       Utility function to read a single RTSolution data record
!

  FUNCTION Read_Record( &
    fid, &  ! Input
    rts) &  ! Output
  RESULT( err_stat )
    ! Arguments
    INTEGER,                    INTENT(IN)  :: fid
    TYPE(CRTM_RTSolution_type), INTENT(OUT) :: rts
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_RTSolution_ReadFile(Record)'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat
    INTEGER :: n_layers

    ! Set up
    err_stat = SUCCESS


    ! Read the dimensions
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) n_layers
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading dimensions - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF


    ! Allocate the RTSolution structure if necessary
    IF ( n_layers > 0 ) THEN
      CALL CRTM_RTSolution_Create( rts, n_layers )
      IF ( .NOT. CRTM_RTSolution_Associated( rts ) ) THEN
        msg = 'Error creating output object.'
        CALL Read_Record_Cleanup(); RETURN
      END IF
    END IF


    ! Read the sensor info
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      rts%Sensor_Id       , &
      rts%WMO_Satellite_Id, &
      rts%WMO_Sensor_Id   , &
      rts%Sensor_Channel
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading sensor information - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF


    ! Read the RT algorithm name
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      rts%RT_Algorithm_Name
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading RT Algorithm Name'//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF


    ! Read the forward radiative transfer intermediate results
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      rts%SOD                    , &
      rts%Surface_Emissivity     , &
      rts%Up_Radiance            , &
      rts%Down_Radiance          , &
      rts%Down_Solar_Radiance    , &
      rts%Surface_Planck_Radiance
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading scalar intermediate results - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF
    IF ( n_Layers > 0 ) THEN
      READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
        rts%Upwelling_Radiance , &
        rts%Layer_Optical_Depth
      IF ( io_stat /= 0 ) THEN
        msg = 'Error reading array intermediate results - '//TRIM(io_msg)
        CALL Read_Record_Cleanup(); RETURN
      END IF
    END IF


    ! Read the radiative transfer results
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      rts%Radiance              , &
      rts%Brightness_Temperature
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading result data - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF

  CONTAINS

    SUBROUTINE Read_Record_Cleanup()
      CALL CRTM_RTSolution_Destroy( rts )
      CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
      IF ( io_stat /= SUCCESS ) &
        msg = TRIM(msg)//'; Error closing file during error cleanup - '//TRIM(io_msg)
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Read_Record_Cleanup

  END FUNCTION Read_Record


!
! NAME:
!       Write_Record
!
! PURPOSE:
!       Function to write a single RTSolution data record
!

  FUNCTION Write_Record( &
    fid, &  ! Input
    rts) &  ! Input
  RESULT( err_stat )
    ! Arguments
    INTEGER,                    INTENT(IN) :: fid
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: rts
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_RTSolution_WriteFile(Record)'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat

    ! Set up
    err_stat = SUCCESS


    ! Write the data dimensions
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) rts%n_Layers
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing dimensions - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF


    ! Write the sensor info
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      rts%Sensor_Id       , &
      rts%WMO_Satellite_Id, &
      rts%WMO_Sensor_Id   , &
      rts%Sensor_Channel
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing sensor information - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF


    ! Write the sensor info
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      rts%RT_Algorithm_Name
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing RT Algorithm Name'//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF


    ! Write the forward radiative transfer intermediate results
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      rts%SOD                    , &
      rts%Surface_Emissivity     , &
      rts%Up_Radiance            , &
      rts%Down_Radiance          , &
      rts%Down_Solar_Radiance    , &
      rts%Surface_Planck_Radiance
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing scalar intermediate results - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF
    IF ( rts%n_Layers > 0 ) THEN
      WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
        rts%Upwelling_Radiance , &
        rts%Layer_Optical_Depth
      IF ( io_stat /= 0 ) THEN
        msg = 'Error writing array intermediate results - '//TRIM(io_msg)
        CALL Write_Record_Cleanup(); RETURN
      END IF
    END IF


    ! Write the radiative transfer results
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      rts%Radiance              , &
      rts%Brightness_Temperature
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing result data - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF

  CONTAINS

    SUBROUTINE Write_Record_Cleanup()
      CLOSE( fid,STATUS=WRITE_ERROR_STATUS,IOSTAT=io_stat,IOMSG=io_msg )
      IF ( io_stat /= 0 ) &
        msg = TRIM(msg)//'; Error closing file during error cleanup - '//TRIM(io_msg)
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Write_Record_Cleanup

  END FUNCTION Write_Record

END MODULE CRTM_RTSolution_Define
