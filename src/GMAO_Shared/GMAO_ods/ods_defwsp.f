

*...................................................................


      subroutine ODS_DefWSp ( NDim, start, count, Size_WorkSp, NBlocks )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_DefWSp
! 
! !DESCRIPTION:
!     This routine determines how the calling routine will
!     partition the data and allocate the work space that will
!     be used to process a NetCDF hyperslab ( as defined by the
!     arrays, start and count ).  This routine defines each 
!     subset of the partioned data by setting the indices of
!     each one-dimensional array block and NetCDF sub-hyperslab.
!     With the exception of the variable, NBlocks, all results
!     are stored in common ( defined in the header file,
!     ods_worksp.h ).
!
! !INTERFACE:  ODS_DefWSp ( NDim, start, count, Size_WorkSp, NBlocks )
!
! !INPUT PARAMETERS
      integer  NDim         ! number of NetCDF dimensions
      integer  start  ( * ) ! NetCDF file indicies specifying the
                            !   location of the corner of the
                            !   hyperslab where the first of the
                            !   data values will be written.  A
                            !   hyperslab is multidimensional 
                            !   block of data within a NetCDF file
                            !   For this routine, the indices 
                            !   correspond to the entire 
                            !   hyperslab that is to be written.
      integer  count  ( * ) ! The edge lengths of the NetCDF
                            !   file hyperslab.
      integer  Size_WorkSp  ! Size of work space

! !OUTPUT PARAMETER
      integer  NBlocks      ! Number of blocks of values that are
                            !   required to process all data via
                            !   the work space.
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_worksp.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures for the work space
!
! !REVISION HISTORY: 
!     31May1996  C. Redder   Origional version
!     16Sep2002  C. Redder   Added input argument, Max_BlockSz.  Replaced
!                            the variable Max_WorkSp with Max_BlockSz.
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_worksp.h'

*     functions refereced
*     -------------------
      integer ODS_NVal    ! calculates the total number of values
                          ! to be processed

*     Other variables
*     ---------------
      integer iDim        ! index variable for do loop
      integer NVal        ! total number of values to be processed
      integer Min_NBlocks ! minimum value for NBlocks
      integer Min_BlockSz ! minimum size of each block

*     Determine the total number of values to be processed
*     ----------------------------------------------------
      NVal        = ODS_NVal ( NDim, count )

*     Initialize variables
*     --------------------
      Min_NBlocks = min ( NVal, Size_WorkSp )
      Min_BlockSz = 1

*     Return and set NBlocks to zero if
*     Min_NBlocks ( i.e. NVal ) is zero
*     -----------------------------------
      if ( Min_NBlocks .eq. 0 ) then
         NBlocks = 0
         return
      end if

*     for each NetCDF dimension
*     -------------------------
      do 10, iDim = 1, NDim

*        The parameters necessary for partitioning the data are
*        determined when the minimim size of each block (i.e.
*        Min_BlockSz) exceeds the available work space.  Until
*        that occurs ...
*        ------------------------------------------------------
         if ( Min_BlockSz .le. Size_WorkSp ) then

*           Save iDim as NDim_SubSlabs, the number
*           of dimensions in a NetCDF sub-hyperslab.  
*           ----------------------------------------
            NDim_SubSlab      = iDim

*           Determine NBlocks_Cycle for the dimension, iDim 
*           -----------------------------------------------
            NBlocks_Cycle     = count ( iDim ) / Min_NBlocks
            if ( NBlocks_Cycle * Min_NBLocks .ne. count ( iDim ) )
     .         NBlocks_Cycle  = NBlocks_Cycle + 1

*           Initialize NBlocks and determine the
*           the associated value for Max_BlockSz
*           and Max_SlabLen
*           ------------------------------------
            NBlocks           = NBlocks_Cycle
            Max_BlockSz       = Min_NBlocks * Min_BlockSz
            Max_SlabLen       = Min_NBlocks

*        After the sector size is determined ...
*        ---------------------------------------
         else

*           Update NBlocks
*           --------------
            NBlocks           = NBlocks  * count ( iDim )

         end if

*        Store the values of start and count in common.
*        ----------------------------------------------
         start_HSlab ( iDim ) = start ( iDim )
         count_HSlab ( iDim ) = count ( iDim ) 

*        Update the number of sectors and sector size for the
*        next dimension.
*        ----------------------------------------------------
         Min_NBlocks          = Min_NBlocks / count ( iDim )
         Min_BlockSz          = Min_BlockSz * count ( iDim )

 10   continue

*     Determine the number of values associated with, NBlocks_Cycle
*     -------------------------------------------------------------
      NVal_Cycle = ODS_NVal ( NDim_SubSlab, count )

      return
      end
