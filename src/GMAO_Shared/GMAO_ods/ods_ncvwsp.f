
*...................................................................


      subroutine ODS_NCVWSp ( iBlock, NDim, start_blk, count_blk )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_NCVWSp
! 
! !DESCRIPTION: 
!     This routine returns the NetCDF indices and edge lengths
!     that define the given data block as a subslab.  The
!     subslab is a subset within the hyperslab as defined
!     by the arrays, start_HSlab and count_HSlab, that are 
!     stored in common and defined in the header file,
!     ods_worksp.h
!
! !INTERFACE:  ODS_NCVWSp ( iBlock, NDim, start_blk, count_blk )
!
! !INPUT PARAMETERS
      integer  iBlock          ! Block index number ( = ith routine
                               !   call )
      integer  NDim            ! number of NetCDF dimensions
!
! !OUTPUT PARAMETERS
      integer  start_blk ( * ) ! NetCDF file indicies specifying the
                               !   location of the corner of the
                               !   subslab where the first of the
                               !   data values will be written in the
                               !   NetCDF file
      integer  count_blk ( * ) ! The NetCDF edge lengths for the data
                               !   block.
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_worksp.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures for work space
!
! !SEE ALSO:
!     ODS_DefWSp ( Defines the workspace used to process data )
!
! !REVISION HISTORY:
!     31May96   C. Redder   Origional version
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_worksp.h'

*     Other variables
*     ---------------
      integer  iDim          ! dimension index variable for do loop
      integer  iBlock_Cycle  ! subslab (i.e. block) index number
                             !    for the cycle
      integer  iBlock_Dim    ! subslab (i.e. block) index number
                             !    for dimension, iDim

*     Set the parameters for the lowest NetCDF dimensions
*     of the subslab.  These parameters are to be identical
*     to those for the hyperslab
*     -----------------------------------------------------
      do 10, iDim = 1, NDim_SubSlab - 1
         start_blk ( iDim ) = start_HSlab ( iDim )
         count_blk ( iDim ) = count_HSlab ( iDim )
 10   continue

*     Define the starting position parameter for the dimension,
*     NDim_SubSlab, of the subslab.
*     ---------------------------------------------------------
      iBlock_Cycle  =   mod ( iBlock - 1, NBlocks_Cycle ) + 1
      start_blk ( NDim_SubSlab )
     .              =   start_HSlab ( iDim )
     .              + ( iBlock_Cycle - 1 ) * Max_SlabLen

*     Define the edge length for the same dimension
*     ---------------------------------------------
      if ( iBlock_Cycle .lt. NBlocks_Cycle ) then
         count_blk ( NDim_SubSlab )
     .              =   Max_SlabLen
      else  ! The edge length for the last subslab in each
*     ----  ! cycle will be equal or less than the lengths
            ! for each of the preceeding subslabs
            ! --------------------------------------------
         count_blk ( NDim_SubSlab )
     .              =   Count_HSlab   ( NDim_SubSlab )
     .              -   Max_SlabLen * ( NBlocks_Cycle - 1 )
      end if

*     Define the starting position of the subslab for the 
*     highest dimensions.  The edge lengths are set to 1.
*     ---------------------------------------------------
      iBlock_Dim    = ( iBlock - 1 ) / NBlocks_Cycle + 1
      do 20, iDim = NDim_SubSlab + 1, NDim
         start_blk ( iDim ) =   start_HSlab ( iDim ) 
     .                      +   iBlock_Dim - 1
         count_blk ( iDim ) =   1
         iBlock_Dim         = ( iBlock_Dim - 1 )
     .                      /   Count_HSlab ( iDim ) + 1
 20   continue

      return
      end
