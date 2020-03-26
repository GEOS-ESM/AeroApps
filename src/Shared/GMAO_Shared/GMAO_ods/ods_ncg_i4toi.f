
*...................................................................


      subroutine ODS_NCG_I4toI ( ncid,  varid,
     .                           start, count,
     .                           nval,  values, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_NCG_I4toI
! 
! !DESCRIPTION:
!     This routine calls the NetCDF routine, NCVGT, to read four
!     byte integer data from an ODS/HDF file and converts the 
!     data to native integers.  This routine insures that the
!     data are correctly returned from the NetCDF routine,
!     irrespective of the machine on which the software is running. 
!
! !INTERFACE:  ODS_NCG_I4toI ( ncid,  varid,
!                                     start, count,
!                                     nval,  values, ierr )
!
! !INPUT PARAMETERS: 
      integer            ncid         ! NetCDF file id
      integer            varid        ! NetCDF variable id
      integer            start  ( * ) ! NetCDF file indicies specifying the
                                      !   location of the corner of the
                                      !   hyperslab where the first of the
                                      !   data values will be read.  A
                                      !   hyperslab is a space within
                                      !   ODS/HDF file where the data are
                                      !   to be read.
      integer            count  ( * ) ! The edge lengths of the NetCDF
                                      !   file hyperslab.
      integer            nval         ! Number of values to be read
!
! !OUTPUT PARAMETERS:
      integer            values ( * ) ! Values to be read
      integer            ierr         ! Return error code
!
! !LIBRARIES ACCESSED:
!     NetCDF
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_worksp.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures for work space
!
! !SEE ALSO:
!     ODS_NCP_ItoI4  ( The inverse companion of this routine )
!
! !REVISION HISTORY: !
!     29Mar96    C. Redder   Origional code
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_worksp.h'

*     Other variables
*     ---------------
      integer     ival     ! Index variable for do loop

*     Default value for error code
*     ----------------------------
      ierr    = NCNoErr

*     Use NetCDF routine to read the values from the file
*     and save the values in the work space
*     ---------------------------------------------------
      call NCVGT  ( ncid,   varid, 
     .              start,  count,
     .              I4_Val, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Convert to native format
*     ------------------------
      do 10, ival = 1, nval
         values ( ival ) = I4_Val ( ival )
 10   continue

      return
      end
