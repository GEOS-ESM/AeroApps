
*...................................................................


      integer function ODS_NVal ( ndim, count )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_NVal
! 
! !DESCRIPTION: 
!     This routine determines the total number of values written to
!     or read from a NetCDF file.  The number of values is a 
!     product of the number for each NetCDF dimension for a variable.
!
! !INTERFACE: nval = ODS_NVal ( ndim, count )
!
! !INPUT PARAMETERS:
      integer ndim         ! Number of NetCDF dimensions
      integer count ( * )  ! Number of values for each dimension
!
! !REVISION HISTORY: 
!     13Oct95   C. Redder   Origional version
!
!-------------------------------------------------------------------------

*     Other variables
*     ---------------
      integer nval   ! temporary storage for the number of values
      integer idim   ! index variable for do loop

*     Initialize nval if ndim is reasonable
*     -------------------------------------
      if ( ndim .ge. 1 ) then
         nval = count ( 1 )

*     Otherwise, set numbers of values to zero and return
*     ---------------------------------------------------
      else
         ODS_NVal = 0
         return

      end if

*     Determine nval
*     --------------
      do 10, idim = 2, ndim
         nval = nval * count ( idim )
 10   continue

      ODS_NVal = nval

      return
      end
