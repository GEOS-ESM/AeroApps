

*...................................................................


      subroutine ODS_CheckI ( ncid, varid, nval, values, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_CheckI
! 
! !DESCRIPTION: 
!     This routine checks native integers to determine whether all
!     values are within the required range as specified by the
!     NetCDF file attributes.  If a missing value is also defined,
!     then any value equal to the missing value is not considered
!     to be out of range.  If a value is determined to be out of
!     range then ierr is given a value of NCEInVal as defined in
!     the header file, netcdf.inc
!   
! !INTERFACE:
!     ODS_CheckI ( ncid, varid, nval, values, ierr )
!
! !INPUT PARAMETERS:
      integer    ncid              ! NetCDF file id
      integer    varid             ! NetCDF variable id
      integer    nval              ! Number of values to scaled
      integer    values  ( nval )  ! Values to be checked
!
! !OUTPUT PARAMETER:
      integer    ierr              ! Returned error (status) code
!
! !SEE ALSO:
!     ODS_CheckR     ( The parallel routine for real numbers )  
!
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!     ods_worksp.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures for work space
!
! !LIBRARIES ACCESSED:
!     NetCDF
!
! !REVISION HISTORY: 
!     16May96   C. Redder   Origional version
!     16Feb2000 R. Todling  Rename stdio.h to ods_stdio.h
!     02Mar2005 D. Dee      More informative error message
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_stdio.h'
      include 'ods_worksp.h'

*     variables defining the allowed range of each integer
*     ----------------------------------------------------
      integer   valid_max             ! maximim value
      integer   valid_min             ! minimum value
      integer   valid_range ( 2 )     ! valid range of values
                                      !   The first element contains
                                      !   the minimum value and the
                                      !   second contains the maximum
                                      !   value.
      integer   IMin                  ! temp for minimum value
      integer   IMax                  ! temp for maximum value

*     variables defining the missing value
*     ------------------------------------
      integer   missing_val           ! code for missing value
      integer   miss_val              ! temporary storage

*     variables containing information about the NetCDF variable
*     and its attributes
*     ----------------------------------------------------------
      character VarNam * ( MaxNCNam ) ! name of variable
      integer   NC_VarType            ! type of NetCDF variable
      integer   NVDims                ! number of NetCDF dimensions
      integer   VDims    ( MaxVDims ) ! NetCDF variable dimensions
      integer   NVAtts                ! number of variable attribute
      integer   AttLen                ! length of an attribute

*     Other variables
*     ---------------
      integer   ival                  ! index variable
      integer   val                   ! temporary storage for values
      integer   ncopts                ! NetCDF error handling options

*     Get information about the variable
*     ----------------------------------
      call ncvinq ( ncid,   varid,
     .              VarNam, NC_VarType,
     .              NVDims, VDims,
     .              NVAtts, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Default error code
*     ------------------
      ierr = NCNoErr

*     Save current options of error handling and turn off
*     error messages and set errors to be non-fatal.
*     ---------------------------------------------------
      call ncgopt ( ncopts )
      call ncpopt ( 0 )

*     Get valid range
*     ---------------
      call ODS_NCAGTI ( ncid, varid, 'valid_range',
     .                        AttLen, valid_range, ierr  )
      if ( ierr  .ne. NCNoErr .and.
     .     ierr  .ne. NCENoAtt ) then
         write ( stderr, 901 ) 'valid_range'
         return

      end if

*     If the attribute, valid_range, exists then ...
*     ----------------------------------------------
      if ( ierr .eq. NCNoErr ) then

*        Check to determine if the maximum and minimum
*        values stored in valid_range are in reverse order
*        -------------------------------------------------
         if ( valid_range ( 2 ) .lt. valid_range ( 1 ) ) then
            valid_min         = valid_range ( 2 )
            valid_max         = valid_range ( 1 )
            valid_range ( 1 ) = valid_min
            valid_range ( 2 ) = valid_max
         end if

      else
*     ----

*        Set defaults for the range of valid values to the
*        minimum and maximum values for native integers
*        -------------------------------------------------
         valid_range ( 1 ) = I_Min
         valid_range ( 2 ) = I_Max

*        Get the attribute, valid_min
*        ----------------------------
         call ODS_NCAGTI ( ncid, varid, 'valid_min',
     .                           AttLen, IMin,      ierr  )
         if ( ierr  .ne. NCNoErr .and.
     .        ierr  .ne. NCENoAtt ) then
            write ( stderr, 901 ) 'valid_min'
            return
         end if
         if ( ierr  .eq. NCNoErr  ) valid_range ( 1 ) = IMin

*        Get the attribute, valid_max
*        ----------------------------
         call ODS_NCAGTI ( ncid, varid, 'valid_max',
     .                           AttLen, IMax,      ierr  )
         if ( ierr  .ne. NCNoErr .and.
     .        ierr  .ne. NCENoAtt ) then
            write ( stderr, 901 ) 'valid_max'
            return
         end if
         if ( ierr  .eq. NCNoErr  ) valid_range ( 2 ) = IMax

      end if

*     Set the valid max and min
*     -------------------------
      valid_min   = valid_range ( 1 )
      valid_max   = valid_range ( 2 )

*     Set defaults for the missing value so that if the attribute
*     specifying the missing value is absent, then the value
*     will have no affect on the results of the check.
*     -----------------------------------------------------------
      missing_val = valid_max

*     Get the attribute, missing_value, if it exits
*     ---------------------------------------------
      call ODS_NCAGTI ( ncid, varid, 'missing_value',
     .                        AttLen, miss_val,     ierr  )
      if ( ierr  .ne. NCNoErr .and.
     .     ierr  .ne. NCENoAtt ) then
         write ( stderr, 901 ) 'missing_val'
         return
      end if
      if ( ierr  .eq. NCNoErr  ) missing_val = miss_val

*     Return error handling option to their previous values
*     -----------------------------------------------------
      call ncpopt ( ncopts )

*     Set default for returned error code
*     -----------------------------------
      ierr = NCNoErr

*     Check the the attribute values to determine
*     if the values are within the user-specified range
*     -----------------------------------------------------
      do 10, ival = 1, nval
         val = values ( ival )
         if ( val  .lt. valid_min .or.
     .        val  .gt. valid_max ) then

*           A value can fail the check only if the
*           value is not missing
*           --------------------------------------
            if ( val .ne. missing_val ) ierr = NCEInVal

         end if
 10   continue

*     Return with an error message if any value is out of range
*     ---------------------------------------------------------
      if ( ierr .ne. NCNoErr ) then
         write ( stderr, 902 ) VarNam, val
         return
      end if

      return
*     ------

 901  format ( /, ' ODS_CheckI: Error in extracting an attribute. ',
     .         /, '             Attribute name is ', a )
 902  format ( /, ' ODS_CheckI: Value to be written is out of ',
     .         /, '             the range as specified by the',
     .         /, '             NetCDF file attributes. ',
     .         /, '             Variable name is ', a,
     .         /, '             Value is ', i12 )

      end
