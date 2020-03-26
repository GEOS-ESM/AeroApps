

*...................................................................


      subroutine ODS_NCVPTI ( ncid,   varid,
     .                        start,  count, 
     .                        values, ierr ) 

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE: ODS_NCVPTI
! 
! !DESCRIPTION:
!     This routine writes an array of native integers to a
!     file using NetCDF interface routines. The routine invokes
!     other modules to convert the values from native format to
!     a format required by the NetCDF routines, and then writes
!     the data to the file.  If the attributes, scale_factor and
!     add_offset, are defined in the NetCDF file, then the values
!     are also scaled using these attributes.
!
!     NOTE: This routine does perform checks to verify whether
!     the range of values on input is consistent with the
!     internal variable type in the file. For example, an
!     integer variable here could be stored as a two byte integer
!     variable on file.  Therefore, all elements in the array
!     to be stored as two byte integers should have values
!     between -32,767 and 32,767.  Be careful.  If the software
!     does detect a number inconsistent with the variable type,
!     then the routine returns with an error message without
!     saving at least some of the values.  This check prevents
!     overflows from occurring which would produce unexpected
!     results and probable abnormal program termination.  The 
!     range of values for the internal variable types are as
!     follows:
!
!       type                      minimum value    maximum value
!       ----                      -------------    -------------
!       byte                               -128              255
!       short (2 byte) integer          -32,767          -32,767
!       long  (4 byte) integer   -2,147,483,643    2,147,483,643
!       real  (4 bytes )               -3.40e38          3.40e38
!
!     A value of 256 is added to any negative integer greater than
!     -128 that is to be stored as a byte.

! !INTERFACE: 
!     call ODS_NCVPTI ( ncid, varid, start, count, values, ierr )
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
      integer   values ( * ) ! Values to be written
!
! !OUTPUT PARAMETER:
      integer   ierr         ! Return error code
!
!     NOTE: This routine does perform checks to verify whether
!     the range of values on input is consistent with the
!     internal variable type in the file. For example, an
!     integer variable here could be stored as a two byte integer
!     variable on file.  Therefore, all elements in the array
!     to be stored as two byte integers should have values
!     between -32,767 and 32,767.  Be careful.  If the software
!     does detect a number inconsistent with the variable type,
!     then the routine returns with an error message without
!     saving at least some of the values.  This check prevents
!     overflows from occurring which would produce unexpected
!     results and probable abnormal program termination.  The 
!     range of values for the internal variable types are as
!     follows:
!
!       type                      minimum value    maximum value
!       ----                      -------------    -------------
!       byte                               -128              255
!       short (2 byte) integer          -32,767          -32,767
!       long  (4 byte) integer   -2,147,483,643    2,147,483,643
!       real  (4 bytes )               -3.40e38          3.40e38
!
!     This routine treats each byte as unsigned integers
!     with values within the range from 0 to 255, inclusively.
!     An offset value of 256 is added to any negative value.
!
! !LIBRARIES ACCESSED:
!     NetCDF
!
! !Files USED:
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
      integer   iblock                ! block index number for
                                      !   do loop
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

*        Transfer all values in the block to the work space
*        --------------------------------------------------
         do 10, ival = ival_beg, ival_end
            iworksp           = ival - ival_beg + 1
            I_Val ( iworksp ) = values ( ival )
 10      continue

*        Check to determine if the values are within the range
*        as specified by the NetCDF file attributes
*        -----------------------------------------------------
         call ODS_CheckI    ( ncid, varid, nval, I_Val, ierr )
         if ( ierr .ne. NCNoErr ) return

*        Scale the array of values
*        -------------------------
         call ODS_ScaleIRev ( ncid, varid, nval, I_Val, ierr )
         if ( ierr .ne. NCNoErr ) return

*        Convert and write array of values to
*        file if the NetCDF variable type is ...
*        ---------------------------------------

*        one byte integer
*        -----------------
         if      ( NC_VarType .eq. NCByte  ) then

            call ODS_NCP_ItoI1 ( ncid,        varid,
     .                           start_block, count_block,
     .                           nval,        I_Val,
     .                           ierr )

*        two byte integers
*        -----------------
         else if ( NC_VarType .eq. NCShort ) then

            call ODS_NCP_ItoI2 ( ncid,        varid,
     .                           start_block, count_block,
     .                           nval,        I_Val,
     .                           ierr )

*        four byte integers
*        ------------------
         else if ( NC_VarType .eq. NCLong  ) then

            call ODS_NCP_ItoI4 ( ncid,        varid,
     .                           start_block, count_block,
     .                           nval,        I_Val,
     .                           ierr )

*        four byte real numbers
*        ----------------------
         else if ( NC_VarType .eq. NCFloat ) then

            call ODS_NCP_ItoR4 ( ncid,        varid,
     .                           start_block, count_block,
     .                           nval,        I_Val,
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

 20   continue

*     Insure that disk copy of file is updated
*     ----------------------------------------
*      call NCSNC ( ncid, ierr )
*      if ( ierr .ne. NCNoErr ) return

      return
*     ------

 901  format ( /, ' ODS_NCVPTI : not a recognized NetCDF ',
     .                          'variable type ' )
 902  format ( /, ' ODS_NCVPTI : variable name is ', a )

      end
