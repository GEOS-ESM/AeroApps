
*...................................................................


      subroutine ODS_NCG_I1toI ( ncid,   varid,
     .                           start,  count,
     .                           nval,   values, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_NCG_I1toI
! 
! !DESCRIPTION:
!     This routine calls the NetCDF routine, NCVGTC, to read one
!     byte integer data from an ODS/HDF file and converts the 
!     data to native integers.  The ODS package treats each byte
!     as unsigned integers with values within the range from 0 to
!     to 255, inclusively.  An offset value of 256 is added to any
!     negative value.  This routine insures that the data are
!     correctly returned from the NetCDF routine, irrespective of
!     the machine on which the software is running. 
!
! !INTERFACE:  ODS_NCG_I1toI ( ncid,   varid,
!                                     start,  count,
!                                     nval,   values, ierr )
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
!     ODS_NCP_ItoI1  ( The inverse companion of this routine )
!
! !REVISION HISTORY: !
!     29Mar1996   C. Redder   Original code
!     13Sep2002   C. Redder   Changes made to handle one-byte integers
!                             stored in a long character string 
!     27Feb2009   D. Nadeau   Fixed variable type mismatch (char/int1)
!-------------------------------------------------------------------------
      include 'netcdf.inc'
      include 'ods_worksp.h'

*     Other variables
*     ---------------
      integer     ival     ! Index variable for do loop
      integer     itemp    ! Temporary storage in native integers
      integer     i1max    ! One more than the maximum value for
                           !   an integer in unsigned byte format
      parameter ( i1max    = 256 )

*     Default value for error code
*     ----------------------------
      ierr    = NCNoErr

*     Use NetCDF routine to read then values from the file
*     and save the values in work space
*     ----------------------------------------------------
      call NCVGT ( ncid,   varid, 
     .              start,  count,
     .              I1_Val, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Convert from unsigned integer byte to native format
*     ---------------------------------------------------
      do 10, ival = 1, nval
         itemp              = I1_Val ( ival ) 
         if ( itemp .lt. 0 ) itemp = itemp + i1max
         values ( ival )    = itemp
 10   continue

      return
      end
