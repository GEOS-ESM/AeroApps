

*..............................................................


      subroutine ODS_NCOpen ( id, filename, mode, ierr)

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:          ODS_NCOpen
! 
! !DESCRIPTION:      Opens a ODS file, returns the ODS file handle id
!                    and reads pointer data from file.
!
! !INTERFACE: call ODS_NCOpen ( id, filename, mode, ierr )
!
! !INPUT PARAMETERS:
      character * (*)  filename         ! name of ODS file
      character * (*)  mode             ! mode = 'w' open for writing
                                        ! mode = 'r' open for reading
!
! !OUTPUT PARAMETERS:
      integer          id               ! ODS file handle
      integer          ierr             ! error code
!
! !SEE ALSO:
!     ODS_NCCreate ( creates the ODS file )
!     ODS_Close    ( closes  the ODS file )
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_hdf.h, a header file, for defining hardwired constants and
!            defining global variables and setting up data
!            structures
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!
! !LIBRARIES ACCESSED:
!     NetCDF
!
! !REVISION HISTORY:
!     17May1996   C. Redder   Origional version
!     04Jun1998   C. Redder   Patch to a bug found in the NetCDF/HDF
!                             library routine.
!     01Nov1999   C. Redder   Revised code to prevent subscript errors
!                             in character strings
!     16Feb2000   R. Todling  Rename stdio.h to ods_stdio.h
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     functions referenced
*     --------------------
      integer     ODS_Handle  ! Obtains an ODS file handle id
      character * ( 1 )
     .            ODS_Case    ! Sets the case for a string
                              !   of characters
      integer     ODS_StrSize ! Returns the string size or length
                              !   excluding trailing blanks

*     Other variables
*     ---------------
      integer     nc_id       ! temporary storage for NetCDF
                              !   file id
      integer     dimid       ! dimension id
      character * ( max_strlen )
     .            DimName     ! name of dimension
      integer     DimSz       ! dimension size
      integer     varid       ! temporary storage for NetCDF
                              !   variable id
      integer     NameSz      ! number of characters in filename
      integer     AttLen      ! string length for attribute or
                              !   attribute length
      integer     itemp       ! temporary storage for attribute
                              !   value
      character * ( 1 )       ! temporary storage for the input
     .            mode_       !   parameter, mode
      integer     IO_Mode     ! temporary storage for the
                              !   array, IOMode ( defined in
                              !   the header file, ods_hdf.h )
      integer     ierr_temp   ! temporary storage for the
                              !   returned error code


*     Set ierr code to valid 
*     ----------------------
      ierr    = NCNoErr

*     Get the number of an used file handle
*     -------------------------------------
      id      = ODS_Handle ( filename, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Set error message flag to ON    ( i.e. prints error mesages )
*     and the fatal error flag to OFF ( i.e. errors are not fatal )
*     -------------------------------------------------------------
      call NCPOPT ( NCVERBOS )          ! verbose mode

*     ensure that any character element in the string
*     parameter, mode, are in lower case
*     -----------------------------------------------
      mode_   = ODS_Case ( mode, 'lower' )
      NameSz  = ODS_StrSize ( filename )

*     Open the file
*     -------------
      if  ( mode_ ( 1:1 ) .eq. 'w' ) then
         nc_id         = NCOPN ( filename, NCWRITE,  ierr )
         if ( ierr .ne. NCNoErr ) return

         IO_Mode       = NCWRITE        ! set IO mode to write 

      else
         nc_id         = NCOPN ( filename, NCNOWRIT, ierr )
         if ( ierr .ne.  NCNoErr ) return

         IO_Mode       = NCNOWRIT       ! set IO mode to read only

      end if

*     A patch to a bug found in the netcdf routine, NCOPN.  The
*     newest version of the routine would print no error message
*     and return a invalid error code for a nonexistant file
*     ----------------------------------------------------------
      if ( nc_id .eq. -1 ) then
         if ( NameSz .gt. 0 ) then
            write ( stderr, 901 ) filename ( : NameSz ), nc_id
         else
            write ( stderr, 901 ) '(blank name)',        nc_id
         end if
         ierr = NCSysErr
         return

      end if

*     Store important file information in common
*     ------------------------------------------
      ncid      ( id )  = nc_id
      IOMode    ( id )  = IO_Mode
      filenames ( id )  = filename

*     Get pointers and related data and store them in common 
*     ------------------------------------------------------
      call ODS_ReSetP ( id, ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      call ODS_ReadP  ( id, ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     Leave file open
*     ---------------

      return
*     ------

*     Clean up by closing the file and
*     setting io mode to CLOSED
*     --------------------------------
 801  continue
      call ncclos ( nc_id, ierr_temp )
      IOMode ( id ) = CLOSED

      return
*     ------

 901  format ( /, ' ODS_NCOpen: Invalid file handle id from the ',
     .         /, '             NetCDF library routine, NCOPN. ',
     .         /, '             file                  = ', a,
     .         /, '             NetCDF file handle id = ', i4 )

      end
