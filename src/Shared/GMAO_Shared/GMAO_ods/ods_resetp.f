
*..............................................................

      subroutine ODS_ReSetP ( id, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:          ODS_ReSetP
! 
! !DESCRIPTION:
!     Initializes or resets all tables of pointers for a NetCDF
!     file (identified by the ODS file handle, id ).  These tables
!     contain all pointers that are necessary to locate a block
!     of data within a file ( specified by the file handle, id )
!     partitioned according to the julian day and synoptic hour.
!     The tables are stored in common defined in the header file,
!     ods_hdf.h.  In addition, the routine extracts relavent
!     NetCDF attributes from the ODS file if the file is opened.
!
! !INTERFACE: call ODS_ReSetP ( id, ierr )
!
! !INPUT PARAMETER:
      integer      id        ! ODS file handle
!
! !OUTPUT PARAMETERS: 
      integer      ierr      ! The return error code
!
! !SEE ALSO:
!     ODS_UpdateP ( Updates the pointer data )
!     ODS_ReadP   ( Read all pointer data from file. )
!     ODS_WriteP  ( Write all pointer data to file. )
!     ODS_GetP    ( Gets pointer data for reading/writing
!                   a block of data. )
!     ODS_GetAttP ( Gets the attribute of the pointers )
!     ODS_JulHr   ( Returns the number of hours from the 
!                   Julian hour of the last block of data
!                   written to file )
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
!     24Apr96   Redder   Origional version
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Variable storing information about the ...
*     ------------------------------------------

*     NetCDF file id
*     --------------
      integer     nc_id

*     NetCDF dimension
*     ----------------
      integer     dimid       ! dimension id
      character * ( max_strlen )
     .            DimName     ! name of dimension
      integer     DimSz       ! dimension size

*     NetCDF variable
*     ---------------
      integer     varid       ! variable id
      character * ( max_strlen )
     .            VarName     ! name of variable
      integer     VarType     ! variable type
      integer     NDims       ! number of associated dimensions
      integer     DimIDs ( 2 )! dimension ids
      integer     NAtts       ! number of associated attributes

*     NetCDF attribute
*     ----------------
      integer     AttLen      ! string length for attribute or
                              !   attribute length

*     Other varibles
*     --------------
      integer     itemp       ! temporary storage for attribute
                              !   value
      integer     iday, isyn  ! index variable for do loops

*     Reset all tables of pointers
*     ----------------------------
      do 10, iday = 1, mdays
         days ( iday, id ) = CLEAR
 10   continue

      do 30, isyn = 1, msyn
      do 20, iday = 1, mdays
         syn_beg ( isyn, iday, id ) = CLEAR + 1
         syn_len ( isyn, iday, id ) = CLEAR
 20   continue
 30   continue

*     Reset the arrays containing the pointers for append mode
*     --------------------------------------------------------
      append_mode ( id ) = OFF
      append_beg  ( id ) = CLEAR + 1
      append_len  ( id ) = 0

*     If file is not opened, then return
*     ----------------------------------
      if ( IOMode ( id ) .eq. CLOSED ) return

*     else ...
*     --------

*     Get NetCDF file id
*     ------------------
      nc_id = ncid ( id )

*     Get information abount the NetCDF variable
*     ------------------------------------------
      varid = NCVID   ( nc_id,  'syn_beg', ierr )
      if ( ierr .ne. NCNoErr ) return
      call ncvinq     ( nc_id,   varid,
     .                  VarName, VarType,
     .                  NDims,   DimIDs,
     .                  NAtts,   ierr )
      if ( ierr .ne. NCNoErr ) return

*     Get NetCDF Dimension
*     --------------------
      dimid        = DimIDs ( 1 )
      call NCDINQ ( nc_id, dimid, DimName, DimSz, ierr )
      if ( ierr .ne. NCNoErr ) return
      nsyn  ( id ) = DimSz

      dimid        = DimIDs ( 2 )
      call NCDINQ ( nc_id, dimid, DimName, DimSz, ierr )
      if ( ierr .ne. NCNoErr ) return
      ndays ( id ) = DimSz

*     Get attribute values and information
*     related to the use of pointers
*     ------------------------------------

*     Julian day offset
*     -----------------
      call ODS_NCAGTI ( nc_id,  varid,
     .                         'first_julian_day',
     .                          AttLen,
     .                          itemp,
     .                          ierr )
      julian_offset ( id ) = itemp - 1
      if ( ierr .ne. NCNoErr ) return

*     latest Julian day for which there is data
*     -----------------------------------------
      call ODS_NCAGTI ( nc_id,  varid,
     .                         'latest_julian_day',
     .                          AttLen,
     .                          itemp,
     .                          ierr )
      if ( ierr .ne. NCNoErr ) return
      latest_day    ( id ) = itemp

*     latest synoptic hour for which there is data
*     --------------------------------------------
      call ODS_NCAGTI ( nc_id,  varid,
     .                         'latest_synoptic_hour',
     .                          AttLen,
     .                          itemp,
     .                          ierr )
      if ( ierr .ne. NCNoErr ) return
      latest_hour   ( id ) = itemp

*     no more attributes
*     ------------------

      return
      end
