
*...................................................................


      subroutine ODS_NCVGTC ( ncid,   varid,
     .                        start,  count, 
     .                        values, ierr ) 

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE: ODS_NCVGTC
! 
! !DESCRIPTION:
!     This routine reads character data from a file using NetCDF interface
!     routines and then stores the data into an array of character strings.
!     If the string length of each array element is larger than the number
!     of characters read for each element, then each element is padded with
!     blanks.  If the string length is smaller than the number of 
!     characters read, then the number of characters read is reduced.
!
! !INTERFACE: 
!     call ODS_NCVGTC ( ncid, varid, start, count, values, ierr )
!
      integer   ncid         ! NetCDF file id
      integer   varid        ! NetCDF variable id
      integer   start  ( * ) ! NetCDF file indicies specifying the
                             !   location of the corner of the
                             !   hyperslab where the first of the
                             !   data values will be written.  A
                             !   hyperslab is a multidimensional
                             !   block of data within a NetCDF
                             !   file.
      integer   count  ( * ) ! The edge lengths of the NetCDF
                             !   file hyperslab.
!
! !OUTPUT PARAMETERS:
      character      * ( * )
     .          values ( * ) ! Strings to be read.
      integer   ierr         ! Return error code
!
!     NOTE: The first element in the input array, count, is assumed to be
!           the number of characters to be read for each element in the
!           output array, values. 
!
! !LIBRARIES ACCESSED:
!     NetCDF
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_worksp.h, a header file, for defining hardwired
!            constants and defining global variables and
!            setting up data structures for work space
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!
! !REVISION HISTORY: 
!     16Sep2002   C. Redder   Origional version
!     10Jan2003   C. Redder   Added helpful comments and removed 
!                             unnecessary check
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_worksp.h'
      include 'ods_stdio.h'

*     NetCDF Variable information
*     ---------------------------
      character VarNam * ( MaxNCNam ) ! name of variable
      integer   NC_VarType            ! type of NetCDF variable
      integer   NVDims                ! number of NetCDF dimensions
      integer   VDims    ( MaxVDims ) ! NetCDF dimensions
      integer   NVAtts                ! number of attributes

*     Parameters determining the usage of work space
*     ----------------------------------------------
      integer   nblocks               ! number of blocks of values
                                      !   that are required to
                                      !   read all data via the
                                      !   work space.  This
                                      !   parameter is also the
                                      !   number of NetCDF I/O
                                      !   routines calls
      integer   iblock                ! block index number for
                                      !   do loop
      integer   nchar                 ! number of characters in the block
      integer   ichar                 ! character index number for do loop
      integer   ival                  ! value index number for array
      integer   iworksp               ! index number for work space
      integer   LVStr                 ! number of characters in each
                                      !   string that is to be extracted

*     Default value for error code
*     ----------------------------
      ierr    = NCNoErr

*     Get information about the NetCDF variable
*     -----------------------------------------
      call ncvinq ( ncid,   varid,
     .              VarNam, NC_VarType,
     .              NVDims, VDims,
     .              NVAtts, ierr )
      if ( ierr .ne. NCNoErr ) return

*     The data is being partitioned into subslabs ( or blocks or
*     subsets).  Each subslab is processed using the work space
*     as defined in the header file, ods_worksp.h.  Set up work
*     space and determine NBlocks.
*     ----------------------------------------------------------
      call ODS_DefWSp ( NVDims, start, count, Max_CWorkSp, nblocks )

*     For each block of values...
*     ---------------------------
      do 20, iblock = 1, nblocks

*        Get indecies defining block in 1-D array and NetCDF hyperslab
*        -------------------------------------------------------------
         call ODS_ValWSp ( iblock,         ichar_beg,   ichar_end   )
         call ODS_NCVWSp ( iblock, NVDims, start_block, count_block )

*        Determine the number of values in the block
*        -------------------------------------------
         nchar = ichar_end - ichar_beg + 1

*        Use NetCDF routine to read the strings from the file
*        and save the values in work space
*        ----------------------------------------------------
         if ( NC_VarType .eq. NCChar ) then
            call NCVGTC ( ncid,  varid, 
     .                    start_block,
     .                    count_block,
     .                    C_Val, nchar, ierr )
            if ( ierr .ne. NCNoErr ) return

         else

            write ( stderr, 901 )
            ierr = NCEBadTy

         end if

*        If there is an error, print out variable name
*        --------------------------------------------- 
         if ( ierr .ne. NCNoErr ) then
            write ( stderr, 902 ) VarNam
            return
         end if

*        Transfer data in block to the output array
*        ------------------------------------------
         LVStr = count ( 1 )   ! Size of count and LVStr must be >= 1
                               !   or output from the routine ODS_DefWSp
                               !   and ODS_NVal would prevent the
                               !   current loop and these statements
                               !   from being executed.
         ival  = ( ichar_beg - 1 ) / LVStr
         do 10, ichar = ichar_beg, ichar_end, LVStr
            iworksp         = ichar - ichar_beg + 1
            ival            = ival  + 1
            values ( ival ) = C_Val ( iworksp : iworksp + LVStr - 1 )
 10      continue

 20   continue

      return
*     ------

 901  format ( /, ' ODS_NCVGTC : not a recognized NetCDF ',
     .                          'variable type for a character ',
     .                          'string ' )
 902  format ( /, ' ODS_NCVGTC : variable name is ', a )

      end
