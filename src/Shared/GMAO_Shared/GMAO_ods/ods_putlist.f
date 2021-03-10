
*..............................................................


      subroutine ODS_PutList ( id, ListName, ListSz, List, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:    ODS_PutList
! 
! !DESCRIPTION: 
!     Writes the user-defined list or array of strings ( ident-
!     ified by the character string, ListName) to the NetCDF data
!     file ( identified by the ODS file handle id ).  Any allocated 
!     space in the NetCDF file not overwritten are filled with
!     blanks.
!
! !INTERFACE:
!     call ODS_PutList ( id, ListName, ListSz, List, ierr )
!
! !INPUT PARAMETERS:
      integer          id          ! ODS file handle
      character  * (*) ListName    ! The name of the list to be
                                   !   saved.  The case of each
                                   !   letter is significant.
      integer          ListSz      ! The number of entries in the
                                   !   list (including blanks )
      character  * (*) List  ( * ) ! The entries of the list.
!
! !OUTPUT PARAMETER: 
      integer          ierr        ! Error code. If non-zero, an error 
                                   !  has occurred. For a list of
                                   !  possible values see Table 8
                                   !  of da Silva and Redder (1995).
!
!     note: The routine checks to ensure that the size of the
!           list does not exceed the space allocated in the 
!           NetCDF file.
!
! !SEE ALSO: 
!     ODS_GetList ( Read the user-defined list string or array
!                   of strings to the NetCDF data file. )
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
!     08Mar1996   C. Redder   Origional version
!     09Apr1998   C. Redder   Added checks to input arguments
!     20Apr1998   C. Redder   Correct bug in bounds check
!     26May1999   C. Redder   Fixed bug in handling the case when ListSz
!                             is zero or less.
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     corners and edge lengths
*     ------------------------
      integer   corner ( MAXNCDIM )
      integer   edges  ( MAXNCDIM )

*     Other variables
*     ---------------
      character VarNam * ( MaxNCNam )
                            ! name of variable (i.e. list )
                            !   returned from NetCDF routine
                            !   ncvinq
      integer   nc_id       ! temporary storage for NetCDF
                            !   file id
      integer   varid       ! temporary storage for NetCDF
                            !   variable id
      integer   NC_VarType  ! type of NetCDF variable
      integer   NVDims      ! number of NetCDF dimensions
      integer   VDims    ( MaxVDims )
                            ! NetCDF variable dimension ids
      integer   DimSz       ! Dimension size
      integer   DimID       ! Dimension id
      integer   NVAtts      ! number of variable attribute
      integer   NC_NStr     ! the allocated number of list
                            !   entries allocated in the file
      integer   NC_StrSz    ! the allocated size of each string
                            !   or list entry allocated 
      integer   iList       ! list entry index number
      integer   iChar0      ! location of substring in temporary
                            !   storage space
      integer   NBlocks     ! number of blocks or chunks necessary
                            !   to write the data to the file
                            !   by using in scratch space in
                            !   the variable, strtemp, ( defined
                            !   in the header file, ods_hdf.h )
      integer   iBlock      ! block index number for do loop
      integer   NChar       ! total number of characters in a
                            !   block of list entries
      integer   NSectors    ! number of sectors in the scatch
                            !   space.  Each sector is NC_StrSz
                            !   in length
      integer   iSector     ! sector index number for do loop

*     Default value for error code
*     ----------------------------
      ierr         = NCNoErr

*     Check to determine if the file handle id is valid
*     -------------------------------------------------
      if ( id            .lt. 1      .or.
     .     id            .gt. id_max .or.
     .     IOMode ( id ) .eq. CLOSED ) then
         write ( stderr, 901 )
         ierr = NCEBadID
         return
      end if

*     Nothing to do if the list size (ListSz) is zero, ...
*     ----------------------------------------------------
      if ( ListSz .eq. 0 ) then
         return

*     ... but if the list size is negative then return with error status
*     ------------------------------------------------------------------
      else if ( ListSz .lt. 0 ) then
         ierr = ODS_DimErr
         write ( stderr, 902 ) ListName
         return

      end if

*     Set NetCDF file id
*     ------------------
      nc_id        = ncid  ( id )

*     Set NetCDF variable id
*     ----------------------
      varid        = ncvid ( nc_id, ListName, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Get information about the NetCDF variable
*     -----------------------------------------
      call ncvinq ( nc_id,  varid,
     .              VarNam, NC_VarType,
     .              NVDims, VDims,
     .              NVAtts, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Determine the hyperslab dimension in the NetCDF file
*     ----------------------------------------------------
      DimID        = VDims ( 1 )
      call NCDINQ ( nc_id, DimID, VarNam, DimSz, ierr )
      if ( ierr .ne. NCNoErr ) return
      NC_StrSz     = DimSz

      DimID        = VDims ( 2 )
      call NCDINQ ( nc_id, DimID, VarNam, DimSz, ierr )
      if ( ierr .ne. NCNoErr ) return
      NC_NStr      = DimSz

*     If the declared length of each string in the array, String,
*     is larger than allocted scratch space, then print error
*     message and return
*     -----------------------------------------------------------
      if ( len ( List ( 1 ) ) .gt. MChar ) then
         write ( stderr, 903 ) len ( List ( 1 ) ), ListName
         ierr      = NCSysErr
         return
      end if

*     If the list size is larger than the number of string
*     allocated in the NetCDF file, then print error message
*     and return
*     ------------------------------------------------------
      if ( ListSz .gt. NC_NStr ) then
         write ( stderr, 904 ) ListName
         ierr      = NCECOORD
         return
      end if

*     Define a block of data and allocate
*     the work space, strtemp
*     -----------------------------------
      NChar        = NC_StrSz * NC_NStr
      NSectors     = min ( NC_StrSz * NC_NStr, MChar ) / NC_StrSz
      NBlocks      = NC_NStr / NSectors
      if ( NBlocks * NSectors .ne. NC_NStr )
     .   NBlocks   = NBlocks + 1

*     Initialize do loop index and pointers
*     -------------------------------------
      iList        = 0
      corner ( 1 ) = 1
      corner ( 2 ) = 1
      edges  ( 1 ) = NC_StrSz
      edges  ( 2 ) = 1

*     for each block ...
*     ------------------
      do 20, iBlock = 1, NBlocks

*        Copy the sublist of NetCDF entries to scatrch space
*        ---------------------------------------------------
         iChar0 = CLEAR
         do 10, iSector = 1, NSectors
            iList = iList + 1

*           If the NetCDF list index number less
*           than the size of the input list ...
*           ------------------------------------
            if ( iList .le. ListSz ) then

*              then write the input list entry to file
*              ---------------------------------------
               strtemp ( iChar0 + 1 : iChar0 + NC_StrSz )
     .            = List ( iList )

            else
*              else write an entry of blanks
*              -----------------------------
               strtemp ( iChar0 + 1 : iChar0 + NC_StrSz )
     .            = ' '

            end if

            iChar0 = iChar0 + NC_StrSz

 10      continue

*        Define the hyperslab of the block of list entries
*        -------------------------------------------------
         edges ( 2 ) = NSectors
         if ( iBlock .eq. NBlocks )
     .      edges ( 2 ) = mod ( NC_NStr - 1, NSectors ) + 1
         NChar       = edges ( 1 ) * edges ( 2 )

*        Write the block of entries to file
*        ----------------------------------
         call ncvptc ( nc_id, varid, corner, edges,
     .                 strtemp, NChar, ierr )
         if ( ierr .ne. NCNoErr ) return

*        Define the corner of the hyperslab
*        for the next block of entries
*        ----------------------------------
         corner ( 2 ) = corner ( 2 ) + edges ( 2 )

 20   continue

      return
*     ------

 901  format ( /, ' ODS_PutList: File handle id number does not ',
     .         /, '              correspond to an opened ODS file' )
 902  format ( /, ' ODS_PutList: The size of the input list ',
     .         /, '              must be a non-negative integer. ',
     .         /, '              Name of list is: ', a )
 903  format ( /, ' ODS_PutList: The character string length of ',
     .         /, '              each list entry exceeds the ',
     .         /, '              size of the scratch space. ',
     .         /, '              Increase MChar to at least ', i9,
     .         /, '              Name of list is: ', a )
 904  format ( /, ' ODS_PutList: The number of entries in the ',
     .         /, '              input list exceeds the space '
     .         /, '              allocated in the NetCDF file. ',
     .         /, '              Name of list is: ', a )

      end
