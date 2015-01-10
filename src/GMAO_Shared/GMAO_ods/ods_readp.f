
*..............................................................


      subroutine ODS_ReadP ( id, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:    ODS_ReadP
! 
! !DESCRIPTION: 
!     Reads all tables of pointers to a file identified by the ODS
!     file handle, id.  These tables contain all pointers that are
!     necessary to locate data within the file organized according
!     to the julian day and synoptic hour.  The tables are stored
!     in common as defined in the header file, ods_hdf.h.
!
! !INTERFACE: call ODS_ReadP ( id, ierr )
!
! !INPUT PARAMETER:
      integer          id          ! File handle
!
! !OUTPUT PARAMETER: 
      integer          ierr        ! The return error code
!
!     note: The routine checks to ensure that ndays and
!           nsyn do not exceed mdays and msyn, respectively.
!           If they do, then the routine returns with an error
!           message.
!
! !SEE ALSO:
!     ODS_UpdateP ( Updates the pointer data )
!     ODS_ReSetP  ( Reset pointers and other relevant data )
!     ODS_WriteP  ( Write all pointer data to file. )
!     ODS_Get P   ( Gets pointer data for reading/writing
!                   a block of data. )
!     ODS_JulHr   ( Returns the number of hours from the 
!                   Julian hour of the last block of data
!                   written to file )
!     ODS_GetAttP ( Gets the attribute of the pointers )
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
! !LIBRARIES ACCESSED:
!     NetCDF
!
! !REVISION HISTORY: 
!     15Apr96   C. Redder   Origional version
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     corners and edge lengths
*     ------------------------
      integer     corner ( MAXNCDIM )
      integer     edges  ( MAXNCDIM )

*     Other variables
*     ---------------
      integer    nc_id       ! temporary storage for NetCDF file id
      integer    varid       ! temporary storage for NetCDF variable id
      integer    iptr        ! location of pointer data in storage
      integer    isyn, iday  ! index variables for do loops 

*     Default value for error code
*     ----------------------------
      ierr       = NCNoErr

*     If the nsyn is too large, then print error message and return
*     -------------------------------------------------------------
      if ( nsyn ( id ) .gt. msyn ) then
         write ( stderr, * )
         write ( stderr, * )
     .    ' ODS_ReadP: The allocated work space is too small '
         write ( stderr, * )
     .    '            Increase msyn to at least ', nsyn ( id )
         ierr = NCSysErr
         return
      end if

*     If the ndays is too large, then print error message and return
*     --------------------------------------------------------------
      if ( ndays ( id ) .gt. mdays ) then
         write ( stderr, * )
         write ( stderr, * )
     .    ' ODS_ReadP: The allocated work space is too small. '
         write ( stderr, * )
     .    '            Increase mdays to at least ', ndays ( id )
         ierr = NCSysErr
         return
      end if

*     Get NetCDF file handle
*     ----------------------
      nc_id      = ncid      ( id )

*     Set corner and edges for days
*     -----------------------------
      corner (1) = 1
      corner (2) = 1
      edges  (1) = ndays     ( id )
      edges  (2) = 1

*     Read data from file and store into scatch space
*     -----------------------------------------------
      varid      = ncvid ( nc_id, 'days', ierr )
      if ( ierr .ne. NCNoErr ) return
      call    ODS_NCVGTI ( nc_id, varid, corner,  edges,
     .                                   ptrtemp, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Copy data in days to scatch space
*     ---------------------------------
      do 10, iday = 1, ndays ( id )
         days ( iday, id ) = ptrtemp ( iday )
 10   continue

*     Set edges for syn_beg and syn_len
*     ---------------------------------
      edges  (1) = nsyn      ( id )
      edges  (2) = ndays     ( id )

*     Read data from file and store into scatch space
*     -----------------------------------------------
      varid      = ncvid ( nc_id, 'syn_beg', ierr )
      call ODS_NCVGTI    ( nc_id, varid, corner,  edges,
     .                                   ptrtemp, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Copy data in syn_beg to scatch space
*     ------------------------------------
      iptr = CLEAR
      do 30, iday = 1, ndays ( id )
      do 20, isyn = 1, nsyn  ( id )
         iptr                       = iptr + 1
         syn_beg ( isyn, iday, id ) = ptrtemp ( iptr )
 20   continue
 30   continue

*     Read data from file and store into scatch space
*     -----------------------------------------------
      varid      = ncvid ( nc_id, 'syn_len', ierr )
      call    ODS_NCVGTI ( nc_id, varid, corner,  edges,
     .                                   ptrtemp, ierr)
      if ( ierr .ne. NCNoErr ) return

*     Copy data in scatch space to syn_len
*     ------------------------------------
      iptr = CLEAR
      do 50, iday = 1, ndays ( id )
      do 40, isyn = 1, nsyn  ( id )
         iptr                       = iptr + 1
         syn_len ( isyn, iday, id ) = ptrtemp ( iptr )
 40   continue
 50   continue

*     Re-initialize the arrays containing the pointers
*     for append mode
*     ------------------------------------------------
      append_mode ( id ) = OFF
      call ODS_GetP ( id, 'read',
     .                     latest_day  ( id ),
     .                     latest_hour ( id ),
     .                     append_beg  ( id ),
     .                     append_len  ( id ), ierr )

      return
      end
