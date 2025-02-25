!
! CRTM_Options_Define
!
! Module defining the CRTM Options optional argument data structure
! and containing routines to manipulate it.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 25-Sep-2004
!                       paul.vandelst@noaa.gov
!

MODULE CRTM_Options_Define

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use statements
  USE Type_Kinds           , ONLY: fp, Long, Double
  USE Message_Handler      , ONLY: SUCCESS, FAILURE, WARNING, INFORMATION, Display_Message
  USE Compare_Float_Numbers, ONLY: OPERATOR(.EqualTo.)
  USE File_Utility         , ONLY: File_Open, File_Exists
  USE Binary_File_Utility  , ONLY: Open_Binary_File        , &
                                   WriteGAtts_Binary_File  , &
                                   ReadGAtts_Binary_File   , &
                                   WriteLogical_Binary_File, &
                                   ReadLogical_Binary_File
  USE CRTM_Parameters      , ONLY: RT_ADA
  USE SSU_Input_Define     , ONLY: SSU_Input_type, &
                                   OPERATOR(==), &
                                   SSU_Input_IsValid, &
                                   SSU_Input_Inspect, &
                                   SSU_Input_GetValue, &
                                   SSU_Input_SetValue, &
                                   SSU_Input_ReadFile, &
                                   SSU_Input_WriteFile
  USE Zeeman_Input_Define  , ONLY: Zeeman_Input_type, &
                                   OPERATOR(==), &
                                   Zeeman_Input_IsValid, &
                                   Zeeman_Input_Inspect, &
                                   Zeeman_Input_GetValue, &
                                   Zeeman_Input_SetValue, &
                                   Zeeman_Input_ReadFile, &
                                   Zeeman_Input_WriteFile
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Datatypes
  PUBLIC :: CRTM_Options_type
  ! Operators
  PUBLIC :: OPERATOR(==)
  ! Public procedures
  PUBLIC :: CRTM_Options_Associated
  PUBLIC :: CRTM_Options_Destroy
  PUBLIC :: CRTM_Options_Create
  PUBLIC :: CRTM_Options_IsValid
  PUBLIC :: CRTM_Options_Inspect
  PUBLIC :: CRTM_Options_DefineVersion
  PUBLIC :: CRTM_Options_InquireFile
  PUBLIC :: CRTM_Options_ReadFile
  PUBLIC :: CRTM_Options_WriteFile
  ! ...Inherited procedures
  PUBLIC :: SSU_Input_GetValue
  PUBLIC :: SSU_Input_SetValue
  PUBLIC :: Zeeman_Input_GetValue
  PUBLIC :: Zeeman_Input_SetValue



  ! -------------------
  ! Procedure overloads
  ! -------------------
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE CRTM_Options_Equal
  END INTERFACE OPERATOR(==)


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PRIVATE, PARAMETER :: MODULE_VERSION_ID = &
  '$Id$'
  ! Literal constants
  REAL(Double), PARAMETER :: ZERO = 0.0_Double
  REAL(Double), PARAMETER :: ONE  = 1.0_Double
  ! Integer "logicals" for I/O
  INTEGER(Long), PARAMETER :: FALSE = 0_Long
  INTEGER(Long), PARAMETER :: TRUE  = 1_Long
  ! Message string length
  INTEGER, PARAMETER :: ML = 256
  ! File status on close after write error
  CHARACTER(*), PARAMETER :: WRITE_ERROR_STATUS = 'DELETE'


  ! ----------------------------
  ! Options data type definition
  ! ----------------------------
  !:tdoc+:
  TYPE :: CRTM_Options_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.

    ! Input checking on by default
    LOGICAL :: Check_Input = .TRUE.

    ! User defined MW water emissivity algorithm
    LOGICAL :: Use_Old_MWSSEM = .FALSE.

    ! Antenna correction application
    LOGICAL :: Use_Antenna_Correction = .FALSE.

    ! NLTE radiance correction is ON by default
    LOGICAL :: Apply_NLTE_Correction = .TRUE.

    ! RT Algorithm is set to ADA by default
    INTEGER(Long) :: RT_Algorithm_Id = RT_ADA

    ! Aircraft flight level pressure
    ! Value > 0 turns "on" the aircraft option
    REAL(Double) :: Aircraft_Pressure = -ONE

    ! User defined number of RT solver streams (streams up + streams down)
    LOGICAL       :: Use_n_Streams = .FALSE.
    INTEGER(Long) :: n_Streams = 0

    ! Scattering switch. Default is for
    ! Cloud/Aerosol scattering to be included.
    LOGICAL :: Include_Scattering = .TRUE.

    ! User defined emissivity/reflectivity
    ! ...Dimensions
    INTEGER(Long) :: n_Channels = 0  ! L dimension
    ! ...Index into channel-specific components
    INTEGER(Long) :: Channel = 0
    ! ...Emissivity optional arguments
    LOGICAL :: Use_Emissivity = .FALSE.
    REAL(Double), ALLOCATABLE :: Emissivity(:)  ! L
    ! ...Direct reflectivity optional arguments
    LOGICAL :: Use_Direct_Reflectivity = .FALSE.
    REAL(Double), ALLOCATABLE :: Direct_Reflectivity(:) ! L

    ! SSU instrument input
    TYPE(SSU_Input_type) :: SSU

    ! Zeeman-splitting input
    TYPE(Zeeman_Input_type) :: Zeeman
  END TYPE CRTM_Options_type
  !:tdoc-:


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Options_Associated
!
! PURPOSE:
!       Elemental function to test the status of the allocatable components
!       of a CRTM Options object.
!
! CALLING SEQUENCE:
!       Status = CRTM_Options_Associated( Options )
!
! OBJECTS:
!       Options:      Options structure which is to have its member's
!                     status tested.
!                     UNITS:      N/A
!                     TYPE:       CRTM_Options_type
!                     DIMENSION:  Scalar or any rank
!                     ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       Status:       The return value is a logical value indicating the
!                     status of the Options members.
!                       .TRUE.  - if the array components are allocated.
!                       .FALSE. - if the array components are not allocated.
!                     UNITS:      N/A
!                     TYPE:       LOGICAL
!                     DIMENSION:  Same as input Options argument
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL FUNCTION CRTM_Options_Associated( self ) RESULT( Status )
    TYPE(CRTM_Options_type), INTENT(IN) :: self
    LOGICAL :: Status
    Status = self%Is_Allocated
  END FUNCTION CRTM_Options_Associated


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Options_Destroy
!
! PURPOSE:
!       Elemental subroutine to re-initialize CRTM Options objects.
!
! CALLING SEQUENCE:
!       CALL CRTM_Options_Destroy( Options )
!
! OBJECTS:
!       Options:      Re-initialized Options structure.
!                     UNITS:      N/A
!                     TYPE:       CRTM_Options_type
!                     DIMENSION:  Scalar OR any rank
!                     ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE CRTM_Options_Destroy( self )
    TYPE(CRTM_Options_type), INTENT(OUT) :: self
    self%Is_Allocated = .FALSE.
  END SUBROUTINE CRTM_Options_Destroy


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Options_Create
!
! PURPOSE:
!       Elemental subroutine to create an instance of the CRTM Options object.
!
! CALLING SEQUENCE:
!       CALL CRTM_Options_Create( Options, n_Channels )
!
! OBJECTS:
!       Options:      Options structure.
!                     UNITS:      N/A
!                     TYPE:       CRTM_Options_type
!                     DIMENSION:  Scalar or any rank
!                     ATTRIBUTES: INTENT(OUT)
!
! INPUTS:
!       n_Channels:   Number of channels for which there is Options data.
!                     Must be > 0.
!                     This dimension only applies to the emissivity-related
!                     components.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Same as Options object
!                     ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE CRTM_Options_Create( self, n_Channels )
    ! Arguments
    TYPE(CRTM_Options_type), INTENT(OUT) :: self
    INTEGER,                 INTENT(IN)  :: n_Channels
    ! Local variables
    INTEGER :: alloc_stat

    ! Check input
    IF ( n_Channels < 1 ) RETURN

    ! Perform the allocation
    ALLOCATE( self%Emissivity(n_Channels), &
              self%Direct_Reflectivity(n_Channels), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN

    ! Initialise
    ! ...Dimensions
    self%n_Channels = n_Channels
    ! ...Arrays
    self%Emissivity          = ZERO
    self%Direct_Reflectivity = ZERO

    ! Set allocation indicator
    self%Is_Allocated = .TRUE.

  END SUBROUTINE CRTM_Options_Create


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Options_IsValid
!
! PURPOSE:
!       Non-pure function to perform some simple validity checks on a
!       CRTM Options object.
!
!       If invalid data is found, a message is printed to stdout.
!
! CALLING SEQUENCE:
!       result = CRTM_Options_IsValid( opt )
!
!         or
!
!       IF ( CRTM_Options_IsValid( opt ) ) THEN....
!
! OBJECTS:
!       opt:       CRTM Options object which is to have its
!                  contents checked.
!                  UNITS:      N/A
!                  TYPE:       CRTM_Options_type
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       result:    Logical variable indicating whether or not the input
!                  passed the check.
!                  If == .FALSE., Options object is unused or contains
!                                 invalid data.
!                     == .TRUE.,  Options object can be used in CRTM.
!                  UNITS:      N/A
!                  TYPE:       LOGICAL
!                  DIMENSION:  Scalar
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION CRTM_Options_IsValid( self ) RESULT( IsValid )
    TYPE(CRTM_Options_type), INTENT(IN) :: self
    LOGICAL :: IsValid
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Options_IsValid'
    CHARACTER(ML) :: msg

    ! Setup
    IsValid = .TRUE.

    ! Check emissivity options
    IF ( self%Use_Emissivity .OR. self%Use_Direct_Reflectivity ) THEN
      IsValid = CRTM_Options_Associated(self)
      IF ( .NOT. IsValid ) THEN
        msg = 'Options structure not allocated'
        CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
        RETURN
      ENDIF
      IF ( self%Use_Emissivity ) THEN
        IF ( ANY(self%Emissivity < ZERO) .OR. ANY(self%Emissivity > ONE) ) THEN
          msg = 'Invalid emissivity'
          CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
          IsValid = .FALSE.
        END IF
      END IF
      IF ( self%Use_Direct_Reflectivity ) THEN
        IF ( ANY(self%Direct_Reflectivity < ZERO) .OR. ANY(self%Direct_Reflectivity > ONE) ) THEN
          msg = 'Invalid direct reflectivity'
          CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
          IsValid = .FALSE.
        END IF
      END IF
    END IF

    ! Check SSU input options
    IsValid = SSU_Input_IsValid( self%SSU ) .AND. IsValid

    ! Check Zeeman input options
    IsValid = Zeeman_Input_IsValid( self%Zeeman ) .AND. IsValid

  END FUNCTION CRTM_Options_IsValid


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Options_Inspect
!
! PURPOSE:
!       Subroutine to print the contents of a CRTM Options object to stdout.
!
! CALLING SEQUENCE:
!       CALL CRTM_Options_Inspect( Options )
!
! INPUTS:
!       Options:       CRTM Options object to display.
!                      UNITS:      N/A
!                      TYPE:       CRTM_Options_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_Options_Inspect( self )
    TYPE(CRTM_Options_type), INTENT(IN) :: self
    WRITE(*,'(1x,"Options OBJECT")')
    ! Display components
    WRITE(*,'(3x,"Check input flag            :",1x,l1)') self%Check_Input
    WRITE(*,'(3x,"Use old MWSSEM flag         :",1x,l1)') self%Use_Old_MWSSEM
    WRITE(*,'(3x,"Use antenna correction flag :",1x,l1)') self%Use_Antenna_Correction
    WRITE(*,'(3x,"Apply NLTE correction flag  :",1x,l1)') self%Apply_NLTE_Correction
    WRITE(*,'(3x,"Aircraft pressure altitude  :",1x,es13.6)') self%Aircraft_Pressure
    WRITE(*,'(3x,"RT algorithm Id             :",1x,i0)') self%RT_Algorithm_Id
    WRITE(*,'(3x,"Include scattering flag     :",1x,l1)') self%Include_Scattering
    WRITE(*,'(3x,"Use n_Streams flag          :",1x,l1)') self%Use_N_Streams
    WRITE(*,'(3x,"n_Streams                   :",1x,i0)') self%n_Streams
    ! ...Emissivity component
    IF ( CRTM_Options_Associated(self) ) THEN
      WRITE(*,'(3x,"Emissivity component")')
      WRITE(*,'(5x,"n_Channels                   :",1x,i0)') self%n_Channels
      WRITE(*,'(5x,"Channel index                :",1x,i0)') self%Channel
      WRITE(*,'(5x,"Use emissivity flag          :",1x,l1)') self%Use_Emissivity
      WRITE(*,'(5x,"Use direct reflectivity flag :",1x,l1)') self%Use_Direct_Reflectivity
      WRITE(*,'(5x,"Emissivity :")')
      WRITE(*,'(5(1x,es13.6,:))') self%Emissivity
      WRITE(*,'(5x,"Direct reflectivity :")')
      WRITE(*,'(5(1x,es13.6,:))') self%Direct_Reflectivity
    END IF
    ! ...SSU input
    CALL SSU_Input_Inspect( self%SSU )
    ! ...Zeeman input
    CALL Zeeman_Input_Inspect( self%Zeeman )

  END SUBROUTINE CRTM_Options_Inspect


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Options_DefineVersion
!
! PURPOSE:
!       Subroutine to return the module version information.
!
! CALLING SEQUENCE:
!       CALL CRTM_Options_DefineVersion( Id )
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

  SUBROUTINE CRTM_Options_DefineVersion( Id )
    CHARACTER(*), INTENT(OUT) :: Id
    Id = MODULE_VERSION_ID
  END SUBROUTINE CRTM_Options_DefineVersion


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Options_InquireFile
!
! PURPOSE:
!       Function to inquire CRTM Options object files.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Options_InquireFile( &
!                        Filename               , &
!                        n_Profiles = n_Profiles  )
!
! INPUTS:
!       Filename:       Character string specifying the name of a
!                       CRTM Options data file to read.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(*)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OPTIONAL OUTPUTS:
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

  FUNCTION CRTM_Options_InquireFile( &
    Filename   , &  ! Input
    n_Profiles ) &  ! Optional output
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN)  :: Filename
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Profiles
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Options_InquireFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat
    INTEGER :: fid
    INTEGER :: m

    ! Set up
    err_stat = SUCCESS
    ! ...Check that the file exists
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


    ! Read the number of profiles dimension
    READ( fid, IOSTAT=io_stat,IOMSG=io_msg ) m
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


    ! Set the optional return arguments
    IF ( PRESENT(n_Profiles) ) n_Profiles = m

  CONTAINS

    SUBROUTINE Inquire_CleanUp()
      IF ( File_Open(fid) ) THEN
        CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
        IF ( io_stat /= SUCCESS ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup - '//TRIM(io_msg)
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Inquire_CleanUp

  END FUNCTION CRTM_Options_InquireFile


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Options_ReadFile
!
! PURPOSE:
!       Function to read CRTM Options object files.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Options_ReadFile( &
!                        Filename               , &
!                        Options                , &
!                        Quiet      = Quiet     , &
!                        n_Profiles = n_Profiles  )
!
! INPUTS:
!       Filename:     Character string specifying the name of an
!                     Options format data file to read.
!                     UNITS:      N/A
!                     TYPE:       CHARACTER(*)
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Options:      CRTM Options object array containing the Options
!                     data.
!                     UNITS:      N/A
!                     TYPE:       CRTM_Options_type
!                     DIMENSION:  Rank-1 (n_Profiles)
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

  FUNCTION CRTM_Options_ReadFile( &
    Filename  , &  ! Input
    Options   , &  ! Output
    Quiet     , &  ! Optional input
    n_Profiles, &  ! Optional output
    Debug     ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),            INTENT(IN)  :: Filename
    TYPE(CRTM_Options_type), INTENT(OUT) :: Options(:)  ! n_Profiles
    LOGICAL,       OPTIONAL, INTENT(IN)  :: Quiet
    INTEGER,       OPTIONAL, INTENT(OUT) :: n_Profiles
    LOGICAL,       OPTIONAL, INTENT(IN)  :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Options_ReadFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat
    LOGICAL :: noisy
    INTEGER :: fid
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
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) n_file_profiles
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading dimensions from '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Check if n_Profiles in file is > size of output array
    n_input_profiles = SIZE(Options)
    IF ( n_file_profiles > n_input_profiles ) THEN
      WRITE( msg,'("Number of profiles, ",i0,", > size of the output Options", &
                  &" array, ",i0,". Only the first ",i0, &
                  &" profiles will be read.")' ) &
                  n_file_profiles, n_input_profiles, n_input_profiles
      CALL Display_Message( ROUTINE_NAME, msg, WARNING )
    END IF
    n_input_profiles = MIN(n_input_profiles, n_file_profiles)


    ! Loop over all the profiles
    Profile_Loop: DO m = 1, n_input_profiles
      err_stat = Read_Record( fid, Options(m), &
                              Quiet = Quiet, &
                              Debug = Debug  )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error reading Options element (",i0,") from ",a)' ) m, TRIM(Filename)
        CALL Read_Cleanup(); RETURN
      END IF
    END DO Profile_Loop


    ! Close the file
    CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
    IF ( io_stat /= 0 ) THEN
      msg = 'Error closing '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF


    ! Set the optional return values
    IF ( PRESENT(n_Profiles) ) n_Profiles = n_input_profiles


    ! Output an info message
    IF ( noisy ) THEN
      WRITE( msg,'("Number of profiles read from ",a,": ",i0)' ) TRIM(Filename), n_input_profiles
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Read_CleanUp()
      IF ( File_Open( Filename ) ) THEN
        CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
        IF ( io_stat /= 0 ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup - '//TRIM(io_msg)
      END IF
      CALL CRTM_Options_Destroy( Options )
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Read_CleanUp

  END FUNCTION CRTM_Options_ReadFile


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Options_WriteFile
!
! PURPOSE:
!       Function to write CRTM Options object files.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Options_WriteFile( Filename     , &
!                                              Options      , &
!                                              Quiet = Quiet  )
!
! INPUTS:
!       Filename:     Character string specifying the name of the
!                     Options format data file to write.
!                     UNITS:      N/A
!                     TYPE:       CHARACTER(*)
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(IN)
!
!       Options:      CRTM Options object array containing the Options
!                     data.
!                     UNITS:      N/A
!                     TYPE:       CRTM_Options_type
!                     DIMENSION:  Rank-1 (n_Profiles)
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

  FUNCTION CRTM_Options_WriteFile( &
    Filename, &  ! Input
    Options , &  ! Input
    Quiet   , &  ! Optional input
    Debug   ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),            INTENT(IN) :: Filename
    TYPE(CRTM_Options_type), INTENT(IN) :: Options(:)  ! n_Profiles
    LOGICAL,       OPTIONAL, INTENT(IN) :: Quiet
    LOGICAL,       OPTIONAL, INTENT(IN) :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Options_WriteFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat
    LOGICAL :: noisy
    INTEGER :: fid
    INTEGER :: m, n_output_profiles

    ! Setup
    err_stat = SUCCESS
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Override Quiet settings if debug set.
    IF ( PRESENT(Debug) ) noisy = Debug


    ! Any valid profiles?
    n_output_profiles = SIZE(Options)
    IF ( n_output_profiles == 0 ) THEN
      msg = 'Zero dimension profiles in input!'
      CALL Write_Cleanup(); RETURN
    END IF


    ! Open the file
    err_stat = Open_Binary_File( Filename, fid, For_Output = .TRUE. )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error opening '//TRIM(Filename)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Write the dimensions
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) n_output_profiles
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing dimensions to '//TRIM(Filename)//'- '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Write the data
    Profile_Loop: DO m = 1, n_output_profiles
      err_stat = Write_Record( fid, Options(m), &
                               Quiet = Quiet, &
                               Debug = Debug )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error writing Options element (",i0,") to ",a)' ) m, TRIM(Filename)
        CALL Write_Cleanup(); RETURN
      END IF
    END DO Profile_Loop


    ! Close the file (if error, no delete)
    CLOSE( fid,STATUS='KEEP',IOSTAT=io_stat,IOMSG=io_msg )
    IF ( io_stat /= 0 ) THEN
      msg = 'Error closing '//TRIM(Filename)//'- '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Output an info message
    IF ( noisy ) THEN
      WRITE( msg,'("Number of profiles written to ",a,": ",i0)' ) TRIM(Filename), n_output_profiles
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

  END FUNCTION CRTM_Options_WriteFile



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
!       CRTM_Options_Equal
!
! PURPOSE:
!       Elemental function to test the equality of two CRTM_Options objects.
!       Used in OPERATOR(==) interface block.
!
!       Note: Only the dimensionality and radiance/brightness temperatures
!             are checked for equality.
!
! CALLING SEQUENCE:
!       is_equal = CRTM_Options_Equal( x, y )
!
!         or
!
!       IF ( x == y ) THEN
!         ...
!       END IF
!
! OBJECTS:
!       x, y:          Two CRTM Options objects to be compared.
!                      UNITS:      N/A
!                      TYPE:       CRTM_Options_type
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

  ELEMENTAL FUNCTION CRTM_Options_Equal( x, y ) RESULT( is_equal )
    TYPE(CRTM_Options_type) , INTENT(IN) :: x, y
    LOGICAL :: is_equal

    is_equal = (x%Check_Input              .EQV.   y%Check_Input           ) .AND. &
               (x%Use_Old_MWSSEM           .EQV.   y%Use_Old_MWSSEM        ) .AND. &
               (x%Use_Antenna_Correction   .EQV.   y%Use_Antenna_Correction) .AND. &
               (x%Apply_NLTE_Correction    .EQV.   y%Apply_NLTE_Correction ) .AND. &
               (x%RT_Algorithm_Id           ==     y%RT_Algorithm_Id       ) .AND. &
               (x%Aircraft_Pressure      .EqualTo. y%Aircraft_Pressure     ) .AND. &
               (x%Use_n_Streams            .EQV.   y%Use_n_Streams         ) .AND. &
               (x%n_Streams                 ==     y%n_Streams             ) .AND. &
               (x%Include_Scattering       .EQV.   y%Include_Scattering    )

    ! Emissivity component
    is_equal = is_equal .AND. &
               ( (x%n_Channels == y%n_Channels) .AND. &
                 (x%Channel    == y%Channel   ) .AND. &
                 (x%Use_Emissivity           .EQV. y%Use_Emissivity          ) .AND. &
                 (x%Use_Direct_Reflectivity  .EQV. y%Use_Direct_Reflectivity ) .AND. &
                 (CRTM_Options_Associated(x) .EQV. CRTM_Options_Associated(y)) )
    IF ( CRTM_Options_Associated(x) .AND. CRTM_Options_Associated(y) ) &
      is_equal = is_equal .AND. &
                 ALL(x%Emissivity          .EqualTo. y%Emissivity         ) .AND. &
                 ALL(x%Direct_Reflectivity .EqualTo. y%Direct_Reflectivity)

    ! SSU input
    is_equal = is_equal .AND. &
               (x%SSU == y%SSU)

    ! Zeeman input
    is_equal = is_equal .AND. &
               (x%Zeeman == y%Zeeman)

  END FUNCTION CRTM_Options_Equal


!
! NAME:
!       Read_Record
!
! PURPOSE:
!       Utility function to read a single options data record
!

  FUNCTION Read_Record( &
    fid        , &  ! Input
    opt        , &  ! Output
    Quiet      , &  ! Optional input
    Debug      ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    INTEGER,                 INTENT(IN)  :: fid
    TYPE(CRTM_Options_type), INTENT(OUT) :: opt
    LOGICAL,       OPTIONAL, INTENT(IN)  :: Quiet
    LOGICAL,       OPTIONAL, INTENT(IN)  :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Options_ReadFile(Record)'
    ! Function variables
    CHARACTER(ML) :: fname
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat
    INTEGER :: n_channels
    LOGICAL :: emissivity_data_present

    ! Set up
    err_stat = SUCCESS


    ! Read the dimensions
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) n_channels
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading dimensions - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ... No emissivity data if n_channels == 0
    emissivity_data_present = (n_channels > 0)



    ! Allocate the Options structure if necessary
    IF ( emissivity_data_present ) THEN
      CALL CRTM_Options_Create( opt, n_channels )
      IF ( .NOT. CRTM_Options_Associated( opt ) ) THEN
        msg = 'Error creating output object.'
        CALL Read_Record_Cleanup(); RETURN
      END IF
    END IF


    ! Read the optional values
    ! ...Input checking logical
    err_stat = ReadLogical_Binary_File( fid, opt%Check_Input )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading input checking option'
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ...Old MWSSEM logical
    err_stat = ReadLogical_Binary_File( fid, opt%Use_Old_MWSSEM )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading old MW water emissivity algorithm switch option'
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ...Antenna correction logical
    err_stat = ReadLogical_Binary_File( fid, opt%Use_Antenna_Correction )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading antenna correction option'
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ...NLTE correction logical
    err_stat = ReadLogical_Binary_File( fid, opt%Apply_NLTE_Correction )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading NLTE correction option'
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ...RT algorithm ID
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) opt%RT_Algorithm_Id
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading RT algorithm id option - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ...Aircraft flight level pressure
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) opt%Aircraft_Pressure
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading aircraft flight level pressure option - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ...Number of RT streams options
    err_stat = ReadLogical_Binary_File( fid, opt%Use_n_Streams )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading n_Streams option'
      CALL Read_Record_Cleanup(); RETURN
    END IF
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) opt%n_Streams
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading n_Streams optional value - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ...Scattering options
    err_stat = ReadLogical_Binary_File( fid, opt%Include_Scattering )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading include scattering option'
      CALL Read_Record_Cleanup(); RETURN
    END IF


    ! Read the emissivity/reflectivity data
    IF ( emissivity_data_present ) THEN
      ! Read the emissivity option
      ! ...The switch...
      err_stat = ReadLogical_Binary_File( fid, opt%Use_Emissivity )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error reading emissivity option'
        CALL Read_Record_Cleanup(); RETURN
      END IF
      ! ...and the data
      READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) opt%Emissivity
      IF ( io_stat /= 0 ) THEN
        msg = 'Error reading emissivity data - '//TRIM(io_msg)
        CALL Read_Record_Cleanup(); RETURN
      END IF

      ! Read the direct reflectivity option
      ! ...The switch...
      err_stat = ReadLogical_Binary_File( fid, opt%Use_Direct_Reflectivity )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error reading direct reflectivity option'
        CALL Read_Record_Cleanup(); RETURN
      END IF
      ! ...and the data
      READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) opt%Direct_Reflectivity
      IF ( io_stat /= 0 ) THEN
        msg = 'Error reading direct reflectivity data - '//TRIM(io_msg)
        CALL Read_Record_Cleanup(); RETURN
      END IF
    END IF


    ! Read the contained object data
    INQUIRE( UNIT=fid,NAME=fname )
    ! ...The SSU input data
    err_stat = SSU_Input_ReadFile( &
                 opt%SSU, &
                 fname, &
                 Quiet    = Quiet, &
                 No_Close = .TRUE., &
                 Debug    = Debug )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading SSU input data'
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ...The Zeeman input data
    err_stat = Zeeman_Input_ReadFile( &
                 opt%Zeeman, &
                 fname, &
                 Quiet    = Quiet, &
                 No_Close = .TRUE., &
                 Debug    = Debug )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading Zeeman input data'
      CALL Read_Record_Cleanup(); RETURN
    END IF

  CONTAINS

    SUBROUTINE Read_Record_Cleanup()
      CALL CRTM_Options_Destroy( opt )
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
!       Utility function to write a single options data record
!

  FUNCTION Write_Record( &
    fid        , &  ! Input
    opt        , &  ! Input
    Quiet      , &  ! Optional input
    Debug      ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    INTEGER,                 INTENT(IN) :: fid
    TYPE(CRTM_Options_type), INTENT(IN) :: opt
    LOGICAL,       OPTIONAL, INTENT(IN) :: Quiet
    LOGICAL,       OPTIONAL, INTENT(IN) :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Options_WriteFile(Record)'
    ! Function variables
    CHARACTER(ML) :: fname
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat

    ! Set up
    err_stat = SUCCESS


    ! Write the dimensions
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) opt%n_channels
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing dimensions - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF


    ! Write the optional values
    ! ...Input checking logical
    err_stat = WriteLogical_Binary_File( fid, opt%Check_Input )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error writing input checking option'
      CALL Write_Record_Cleanup(); RETURN
    END IF
    ! ...Old MWSSEM logical
    err_stat = WriteLogical_Binary_File( fid, opt%Use_Old_MWSSEM )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error writing old MW water emissivity algorithm switch option'
      CALL Write_Record_Cleanup(); RETURN
    END IF
    ! ...Antenna correction logical
    err_stat = WriteLogical_Binary_File( fid, opt%Use_Antenna_Correction )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error writing antenna correction option'
      CALL Write_Record_Cleanup(); RETURN
    END IF
    ! ...NLTE correction logical
    err_stat = WriteLogical_Binary_File( fid, opt%Apply_NLTE_Correction )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error writing NLTE correction option'
      CALL Write_Record_Cleanup(); RETURN
    END IF
    ! ...RT algorithm ID
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) opt%RT_Algorithm_Id
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing RT algorithm id option - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF
    ! ...Aircraft flight level pressure
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) opt%Aircraft_Pressure
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing aircraft flight level pressure option - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF
    ! ...Number of RT streams options
    err_stat = WriteLogical_Binary_File( fid, opt%Use_n_Streams )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error writing n_Streams option'
      CALL Write_Record_Cleanup(); RETURN
    END IF
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) opt%n_Streams
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing n_Streams optional value - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF
    ! ...Scattering options
    err_stat = WriteLogical_Binary_File( fid, opt%Include_Scattering )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error writing include scattering option'
      CALL Write_Record_Cleanup(); RETURN
    END IF


    ! Write the emissivity/reflectivity data
    IF ( CRTM_Options_Associated(opt) ) THEN
      ! Write the emissivity option
      ! ...The switch...
      err_stat = WriteLogical_Binary_File( fid, opt%Use_Emissivity )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error writing emissivity option'
        CALL Write_Record_Cleanup(); RETURN
      END IF
      ! ...and the data
      WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) opt%Emissivity
      IF ( io_stat /= 0 ) THEN
        msg = 'Error writing emissivity data - '//TRIM(io_msg)
        CALL Write_Record_Cleanup(); RETURN
      END IF

      ! Write the direct reflectivity option
      ! ...The switch...
      err_stat = WriteLogical_Binary_File( fid, opt%Use_Direct_Reflectivity )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error writing direct reflectivity option'
        CALL Write_Record_Cleanup(); RETURN
      END IF
      ! ...and the data
      WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) opt%Direct_Reflectivity
      IF ( io_stat /= 0 ) THEN
        msg = 'Error writing direct reflectivity data - '//TRIM(io_msg)
        CALL Write_Record_Cleanup(); RETURN
      END IF
    END IF


    ! Write the contained object data
    INQUIRE( UNIT=fid,NAME=fname )
    ! ...The SSU input data
    err_stat = SSU_Input_WriteFile( &
                 opt%SSU, &
                 fname, &
                 Quiet    = Quiet, &
                 No_Close = .TRUE., &
                 Debug    = Debug )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error writing SSU input data'
      CALL Write_Record_Cleanup(); RETURN
    END IF
    ! ...The Zeeman input data
    err_stat = Zeeman_Input_WriteFile( &
                 opt%Zeeman, &
                 fname, &
                 Quiet    = Quiet, &
                 No_Close = .TRUE., &
                 Debug    = Debug )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error writing Zeeman input data'
      CALL Write_Record_Cleanup(); RETURN
    END IF

  CONTAINS

    SUBROUTINE Write_Record_Cleanup()
      CLOSE( fid,STATUS=WRITE_ERROR_STATUS,IOSTAT=io_stat,IOMSG=io_msg )
      IF ( io_stat /= SUCCESS ) &
        msg = TRIM(msg)//'; Error closing file during error cleanup - '//TRIM(io_msg)
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Write_Record_Cleanup

  END FUNCTION Write_Record

END MODULE CRTM_Options_Define
