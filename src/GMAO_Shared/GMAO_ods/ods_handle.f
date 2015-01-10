
*..............................................................


      integer function ODS_Handle ( filename, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:           ODS_Handle
! 
! !DESCRIPTION: 
!     Returns a ODS file handle number that is available for opening
!     or creating another file.  If no handles are found or if
!     if the file is already open, then the function returns a
!     zero and an invalid error code.
!
! !INTERFACE: id = ODS_Handle ( filename, ierr )
!
! !INPUT PARAMETER:
      character * (*)   filename   ! file name

! !OUTPUT PARAMETERS:
      integer           id         ! available file handle
      integer           ierr       ! The return error code
!
! !SEE ALSO:
!     ODS_Create, ODS_NCCreate,  ODS_Open, ODS_NCOpen
!     ( opens/creates the ODS file )
!     ODS_Close.  Closes a file.
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
!     28Mar1995   C. Redder   Origional version
!     28Jun1998   C. Redder   Fixed bug in checking names of files
!                             already opened.
!     01Nov1999   C. Redder   Revised code to prevent subscript errors
!                             in character strings
!     16Feb2000   R. Todling  Rename stdio.h to ods_stdio.h
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Other variables
*     ---------------
      integer  id_not_used  ! available file handle
      integer  lf           ! number of characters in filename

      ierr        = NCNoErr
      id_not_used = Not_Found

*     Find handle not being used
*     --------------------------
      do 10, id = id_max, 1, -1
         if ( IOMode ( id ) .eq. CLOSED ) then
            id_not_used      = id
            filenames ( id ) = BLANK   ! Clear entries in list of
                                       !   files that are not
                                       !   opened
         end if
 10   continue

      ODS_Handle = id_not_used

*     Return with error message if all handles are in use
*     ---------------------------------------------------
      if ( id_not_used .eq. Not_Found ) then
         lf = index ( filename, ' ' ) - 1
         write ( stderr, * )
         if ( lf .gt. 0 ) then
            write ( stderr, * )
     .       ' ODS_Handle: No additional file handles are ',
     .                    'available for file, ', filename ( :lf ) 
         else
            write ( stderr, * )
     .       ' ODS_Handle: No additional file handles are ',
     .                    'available for file, ', '(blank name)'
         end if

         ierr = NCENFile
         return
      end if

*     Check to determine if file with same filename is already open
*     -------------------------------------------------------------
      do 20, id = 1, id_max
         if ( filenames ( id ) .eq. filenames ( id_not_used ) .and.
     .           IOMode ( id ) .ne. CLOSED ) then

            lf = index ( filenames ( id_not_used ), ' ' ) - 1
            write ( stderr, * )
            if ( lf .gt. 0 ) then
               write ( stderr, * )
     .            ' ODS_Handle: ODS/HDF file, ',
     .              filenames ( id_not_used ) ( 1:lf ),
     .           ', is already open. '
            else
               write ( stderr, * )
     .            ' ODS_Handle: ODS/HDF file, ', '(blank name)',
     .           ', is already open. '

            end if

            ODS_Handle                = Not_Found
            filenames ( id_not_used ) = BLANK
            ierr                      = NCEName
            return
         end if

 20   continue

      filenames ( id_not_used ) = BLANK

      return
      end
