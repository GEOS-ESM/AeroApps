
*...................................................................


      subroutine ODS_NCVGTR ( ncid,   varid,
     .                        start,  count, 
     .                        values, ierr ) 

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE: ODS_NCVGTR
! 
! !DESCRIPTION:
!     This routine reads from a file using NetCDF interface
!     routines and then converts the data to native real.
!     format.  The routine invokes other modules to read the
!     values from the file and then convert the values from a
!     format determined by the NetCDF routines to native format.
!     If the attributes, scale_factor and add_offset, are defined
!     in the NetCDF file, then the values are also scaled using
!     these attributes.
!
! !INTERFACE: 
!     call ODS_NCVGTR ( ncid, varid, start, count, values, ierr )
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
      real      values ( * ) ! Values to be read
      integer   ierr         ! Return error code
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
!     05Jun1996   C. Redder   Origional version
!     13Sep2002   C. Redder   Made changes in calling ods_defwsp in order
!                             process and store one-byte integers in a
!                             multi-character string rather than in an 
!                             array of one-character strings
!
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
                                      !   write all data via the
                                      !   work space.  This
                                      !   parameter is also the
                                      !   number of NetCDF I/O
                                      !   routines calls
      integer   iblock                ! block index number for do
                                      !   loop
      integer   nval                  ! number of values in the block
      integer   ival                  ! value index number for do loop
      integer   iworksp               ! index number for work space

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
*     space and determine NBlocks for ...
*     ----------------------------------------------------------
      if ( NC_VarType .eq. NCByte ) then   ! ... one byte integers
         call ODS_DefWSp ( NVDims, start, count, Max_CWorkSp, nblocks )

      else                                 ! ... for all other types
         call ODS_DefWSp ( NVDims, start, count, Max_WorkSp,  nblocks )

      end if

*     For each block of values...
*     ---------------------------
      do 20, iblock = 1, nblocks

*        Get indecies defining block in 1-D array and NetCDF hyperslab
*        -------------------------------------------------------------
         call ODS_ValWSp ( iblock,         ival_beg,    ival_end    )
         call ODS_NCVWSp ( iblock, NVDims, start_block, count_block )

*        Determine the number of values in the block
*        -------------------------------------------
         nval = ival_end - ival_beg + 1

*        Convert and read array of values from the
*        file if the NetCDF variable type is ...
*        -----------------------------------------

*        one byte integer
*        -----------------
         if      ( NC_VarType .eq. NCByte  ) then

            call ODS_NCG_I1toR ( ncid,        varid,
     .                           start_block, count_block,
     .                           nval,        R_Val,
     .                           ierr )

*        two byte integers
*        -----------------
         else if ( NC_VarType .eq. NCShort ) then

            call ODS_NCG_I2toR ( ncid,        varid,
     .                           start_block, count_block,
     .                           nval,        R_Val,
     .                           ierr )

*        four byte integers
*        ------------------
         else if ( NC_VarType .eq. NCLong  ) then

            call ODS_NCG_I4toR ( ncid,        varid,
     .                           start_block, count_block,
     .                           nval,        R_Val,
     .                           ierr )

*        four byte real numbers
*        ----------------------
         else if ( NC_VarType .eq. NCFloat ) then

            call ODS_NCG_R4toR ( ncid,        varid,
     .                           start_block, count_block,
     .                           nval,        R_Val,
     .                           ierr )

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

*        Scale the values
*        ----------------
         call ODS_ScaleR ( ncid, varid, nval, R_Val, ierr )

*        Transfer data in block to the output array
*        ------------------------------------------
         do 10, ival = ival_beg, ival_end
            iworksp            = ival - ival_beg + 1
            values ( ival )    = R_Val ( iworksp )
 10      continue

 20   continue

      return
*     ------

 901  format ( /, ' ODS_NCVGTR : not a recognized NetCDF ',
     .                          'variable type ' )
 902  format ( /, ' ODS_NCVGTR : variable name is ', a )

      end
