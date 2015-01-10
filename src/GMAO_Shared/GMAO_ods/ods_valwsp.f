
*...................................................................


      subroutine ODS_ValWsp ( iBlock, iValBeg, iValEnd )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_ValWSp
! 
! !DESCRIPTION: 
!     This routine returns the indicies of the first and last 
!     elements of the data block ( in a one-dimensional array )
!     as defined by the subeoutine, ODS_DefWSp
!
! !INTERFACE:  ODS_ValWSp ( iBlock, iValBeg, iValEnd )
!
! !INPUT PARAMETER
      integer iBlock   ! Block index number ( = ith routine
                       !   call )
!
! !OUTPUT PARAMETERS
      integer iValBeg  ! Index number defining the beginning
                       !   of the block of values in a one
                       !   dimensional array
      integer iValEnd  ! Index number defining the end
                       !   of the block of values in a one
                       !   dimensional array
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
!     04Jun96   C. Redder   Origional version
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_worksp.h'

*     Other variables 
*     ---------------
      integer  iCycle       ! cycle index number
      integer  iBlock_Cycle ! block index number within 
                            !   a cycle

*     note: For a discussion of a cycle, see the file, ods_worksp.h

*     Determine iCycle and iBlock_Cycle
*     ---------------------------------
      iCycle       =     ( iBlock - 1 ) / NBlocks_Cycle   + 1
      iBlock_Cycle = mod ( iBlock - 1,    NBlocks_Cycle ) + 1

*     Determine the index number of the first value in the block
*     ----------------------------------------------------------
      iValBeg      = ( iCycle       - 1 ) * NVal_Cycle
     .             + ( iBlock_Cycle - 1 ) * Max_BlockSz + 1

*     Determine the index number of the last value
*     --------------------------------------------
      if ( iBlock_Cycle .lt. NBlocks_Cycle ) then
         iValEnd   =   iValBeg + Max_BlockSz - 1

      else  ! The size for the last block in each cycle
*     ----  ! will be equal or less than the size for
            ! each of the preceeding block
            ! -----------------------------------------
         iValEnd   =   iValBeg - 1
     .             +   NVal_Cycle 
     .             - ( NBlocks_Cycle - 1 ) * Max_BlockSz

      end if

      return
      end
