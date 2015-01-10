
*..............................................................


      character * ( * ) function ODS_VerTag ( version_index, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:           ODS_VerTag
! 
! !DESCRIPTION: 
!     Returns the appropriate version tag from a list given the
!     array list index number.  The list is sorted so that the
!     index of 1 corresponds to the earliest version.  If
!     version index number = 0, then the function returns
!     the index tag corresponding to the latest version.  
!
! !INTERFACE: version_tag = ODS_VerTag ( version_index, ierr )
!
! !INPUT PARAMETER:
      integer           version_index ! version index

! !OUTPUT PARAMETERS:
!     character * ( * ) version_tag   ! version tag
      integer           ierr          ! The return error code
!
! !SEE ALSO:
!     ODS_VerIndex - Returns the ODS/HDF version array list
!                    index number
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!
! !REVISION HISTORY:
!     06Apr98   C. Redder   Origional version
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Default error code is valid
*     ---------------------------
      ierr           = NCNoErr

*     If desired, extract the tag for
*     the latest version and return
*     -------------------------------
      if ( version_index .eq. 0 ) then
         ODS_VerTag = version_list ( NVersions )
         return

      end if

*     Perform an array bound check on
*     the input version list index number
*     -----------------------------------
      if ( version_index .lt. 1 .or.
     .     version_index .gt. NVersions ) then
         ierr = NCSysErr
         write ( stderr, 901 ) version_index, NVersions
         return

      end if

*     Extract the desired version tag
*     -------------------------------
      ODS_VerTag = version_list ( version_index )

 901  format ( /, ' ODS_VerTag: The version list index number, ',
     .         /, '             ', i6, ' is out of bounds. ',
     .         /, '             The bounds are 0 to ', i4, '.' )

      return
      end
