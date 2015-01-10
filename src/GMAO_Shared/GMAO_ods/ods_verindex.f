
*..............................................................

      integer function ODS_VerIndex ( version_tag, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:           ODS_VerIndex
! 
! !DESCRIPTION: 
!     Performs a search in a list of valid ODS version tags and
!     returns the returns the appropriate version index number.
!     The list is sorted so that the index of 1 corresponds
!     to the earliest version.  If version_tag = 'latest version',
!     then the function returns the index number corresponding to
!     the latest version.  If the tag in not found in the list
!     then the function returns a zero and a VALID error code.
!
! !INTERFACE: version_index = ODS_VerIndex ( version_tag, ierr )
!
! !INPUT PARAMETER:
      character * (*)   version_tag   ! version tag

! !OUTPUT PARAMETERS:
      integer           version_index ! version index
      integer           ierr          ! The return error code
!
! !SEE ALSO:
!     ODS_VerTag - Returns the ODS/HDF version tag
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
!     06Apr1998  C. Redder   Origional version
!     01Nov1999  C. Redder   Revised code to prevent subscript errors
!                            in character strings
!     24Feb2000  C. Redder   Minor correction in documentation.
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Fuctions referenced
*     -------------------
      integer  ODS_StrSize   ! Determines the string size
                             !   excluding trailing blanks
      integer  ODS_StrSearch ! Return list entry number from
                             !   a string search

*     Other variables
*     ---------------
      character * ( max_strlen )
     .         tag_temp      ! temporary storage for version tag 
      integer  tag_len       ! string length for version tag 
      integer  tag_index     ! temporary storage for the
                             !   version index number

*     Default output parameters to their defaults
*     -------------------------------------------
      ODS_VerIndex = 0
      ierr         = NCNoErr

*     Check to determine if the input string length is larger enough
*     --------------------------------------------------------------
      tag_len = ODS_StrSize ( version_tag )
      if ( tag_len .gt. max_strlen ) then
         ierr = NCSysErr
         write ( stderr, 901 ) tag_len, max_strlen
         return

      end if


*     If the index number for the latest version is desired
*     then set the index value to the maximum and return
*     -----------------------------------------------------
      tag_temp = version_tag
      tag_len  = ODS_StrSize ( tag_temp )
      if ( tag_temp ( : tag_len ) .eq. 'latest version' ) then
         ODS_VerIndex = NVersions
         return

      end if

*     Search for the tag among the list of valid versions
*     ---------------------------------------------------
      tag_index = ODS_StrSearch ( tag_temp,
     .                            NVersions,
     .                            version_list,
     .                            ierr )

*     If an error condition occurred then
*     set the index number and return
*     -----------------------------------
      if ( ierr .ne. NCNoErr ) then
         ODS_VerIndex = 0
         return

      end if

*     If the version index number is not found
*     then return with an error message
*     ----------------------------------------
      if    ( tag_index .eq. 0 ) then
         ODS_VerIndex = 0
         if ( tag_len   .gt. 0 ) then
            write ( stderr, 902 ) tag_temp ( : tag_len )
         else
            write ( stderr, 902 ) '(blank name)'
         end if
         return

      end if

      ODS_VerIndex = tag_index

 901  format ( /, ' ODS_VerIndex: The length of the string ',
     .         /, '               supplied by the calling routine ',
     .         /, '               (= ', i4, ') is greater than the ',
     .         /, '               allowed maximum (= ', i4, ')' )

 902  format ( /, ' ODS_VerIndex: The version tag, ', a,
     .                         ', not valid. ' )

      return
      end
