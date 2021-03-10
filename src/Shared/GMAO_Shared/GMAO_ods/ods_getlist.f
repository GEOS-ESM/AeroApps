
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:    ODS_GetList --- Reads a list from an ODS file
! 
! !DESCRIPTION: 
!    \label{ODS:GetList}
!     Reads the user-defined list or array of strings
!     (identified by the character string, ListName) from a
!     NetCDF data file (identified by the ODS file handle id).
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
! !INTERFACE:
!
      subroutine ODS_GetList ( id, ListName, ListSz, List, ierr )
!
! !INPUT PARAMETERS:
      implicit         NONE
      integer          id          ! ODS file handle
      character  * (*) ListName    ! The name of the list to be
                                   !   obtained.  The case of each
                                   !   letter is significant.
!
! !INPUT/OUTPUT PARAMETERS:
      integer         ListSz       ! on input: maximum number of
                                   !   values (usually determined
                                   !   by the space allocated for
                                   !   storage)
                                   ! on output: the number of
                                   !   values obtained from the
                                   !   file 
                                   ! note: If the input value is
                                   !   smaller than the number
                                   !   of values to be read, then
                                   !   the routine exits with a
                                   !   system error message and
                                   !   code of ODS_DimErr ( = -3 ).
                                   !   Also, nval is set to the
                                   !   minimum value required in
                                   !   order for the routine to
                                   !   execute successfully.
!
! !OUTPUT PARAMETERS: 
      character  * (*) List  ( * ) ! The entries of the list.
      integer    ierr              ! Error code. If non-zero, an 
                                   !  error has occurred. For a list
                                   !  of possible values, see the
                                   !  description section of this
                                   !  prologue.
!
!     NOTE:  Each file in the present version of ODS ( version 2.00 )
!            contains the lists with the following names
!
!            kt_names    Name of GEOS/DAS data types
!            kt_units    Units for each GEOS/DAS data type
!            kx_names    Name of GEOS/DAS data sources
!            kx_meta     kx specific metadata information.  This 
!                        information is used to specify the meaning
!                        of the ODS variable "xm" or to specify the
!                        name of the external OMS file name.  This
!                        list is available only in versions 2.00
!                        or later.
!            qcx_names   information about the meaning of the ODS 
!                        variable, "qcexcl".  This list is available
!                        only in versions 2.00 or later.
!
! !SEE ALSO: 
!     ODS_PutList ( Put the user-defined list string or array
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
!     08Mar1996   C. Redder   Original version
!     13Apr1998   C. Redder   Made routine to be an interface routine
!                             in ODS version 2.00.  Added checks made
!                             to input arguments.
!     20Apr1998   C. Redder   Fixed bug in bounds check
!     02Nov1999   C. Redder   Fixed bug in storing the list entries
!                             in the last block of data.
!     19Nov1999   C. Redder   Added a latex label in and moved the
!                             subroutine statement into the prologue.
!                             Modified the comments for the return
!                             status code.
!     06Dec1999   C. Redder   Corrections to the documentation in the
!                             prologue.
!     16Feb2000   R. Todling  Rename stdio.h to ods_stdio.h
!
!EOP
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

*     Return with an error message if the number of
*     observations to be read exceeds the user-defined maximum
*     --------------------------------------------------------
      if ( NC_NStr .gt. ListSz ) then
         write ( stderr, * )
         write ( stderr, 902 ) ListSz, NC_NStr
         ierr = ODS_DimErr
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

*        Define the hyperslab of the block of list entries
*        -------------------------------------------------
         edges ( 2 ) = NSectors
         if ( iBlock .eq. NBlocks ) then   ! for the last block
                                           ! ------------------
            edges ( 2 ) = mod ( NC_NStr - 1, NSectors ) + 1
            NSectors    = edges ( 2 )

         end if
         NChar       = edges ( 1 ) * edges ( 2 )

*        Write the block of entries to file
*        ----------------------------------
         call ncvgtc ( nc_id,   varid, corner, edges,
     .                 strtemp, NChar, ierr )
         if ( ierr .ne. NCNoErr ) return

*        Define the corner of the hyperslab
*        for the next block of entries
*        ----------------------------------
         corner ( 2 ) = corner ( 2 ) + edges ( 2 )

*        Copy the block of data in scatch space to the output array
*        ----------------------------------------------------------
         iChar0 = CLEAR
         do 10, iSector = 1, NSectors
            iList  = iList + 1
            List ( iList )
     .          = strtemp ( iChar0 + 1 : iChar0 + NC_StrSz )

            iChar0 = iChar0 + NC_StrSz

 10      continue

 20   continue

*     return the list size
*     --------------------
      ListSz = NC_NStr

      return
*     ------

 901  format ( /, ' ODS_GetList: File handle id number does not ',
     .         /, '              correspond to an opened ODS file' )
 902  format ( /, ' ODS_GetList: The number of list entries ',
     .         /, '              to be read exceeds the user ',
     .         /, '              defined maximum of ', i10, '.',
     .         /, '              Number of list entries to be ',
     .         /, '              read = ', i10, '.' )

      end
