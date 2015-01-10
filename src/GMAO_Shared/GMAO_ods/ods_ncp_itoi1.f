
*...................................................................


      subroutine ODS_NCP_ItoI1 ( ncid,   varid,
     .                           start,  count,
     .                           nval,   values, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_NCP_ItoI1
! 
! !DESCRIPTION: 
!     This routine converts native integers to one byte integers
!     and calls the NetCDF routine, NCVPTC, to write the one byte
!     integers to a ODS/HDF file.  The ODS package treats each
!     byte as unsigned integers with values within the range 
!     from 0 and 255, inclusively.  This routine insures that the
!     data are correctly passed into the NetCDF routine,
!     irrespective of the machine on which the software is running. 
!
! !INTERFACE:  ODS_NCP_ItoI1 ( ncid,   varid,
!                                     start,  count,
!                                     nval,   values, ierr )
!
! !INPUT PARAMETERS: 
      integer            ncid         ! NetCDF file id
      integer            varid        ! NetCDF variable id
      integer            start  ( * ) ! NetCDF file indicies specifying the
                                      !   location of the corner of the
                                      !   hyperslab where the first of the
                                      !   data values will be written.  A
                                      !   hyperslab is a space within
                                      !   ODS/HDF file where the data are
                                      !   to be written.
      integer            count  ( * ) ! The edge lengths of the NetCDF
                                      !   file hyperslab.
      integer            nval         ! Number of values to be written
      integer            values ( * ) ! Values to be written
!
!     note: An offset value of 256 will be added to any negative
!           integer in the array, values.
!
!
! !OUTPUT PARAMETER:
      integer            ierr         ! Return error code
!
! !LIBRARIES ACCESSED:
!     NetCDF
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!     ods_worksp.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures for work space
!
! !SEE ALSO:
!     ODS_NCG_I1toI  ( The inverse companion of this routine )
!
! !REVISION HISTORY: !
!     02Apr1996  C. Redder   Origional code
!     16Feb2000  R. Todling  Rename stdio.h to ods_stdio.h
!     13Sep2002  C. Redder   Changes made to handle one-byte integers
!                            stored in a long character string 
!     27Feb2009  D. Nadeau   Fixed variable type mismatch (char/int1)
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_stdio.h'
      include 'ods_worksp.h'

*     Other variables
*     ---------------
      integer      val     ! Temporary storage
      integer     ival     ! Index variable for do loop

*     Default value for error code
*     ----------------------------
      ierr    = NCNoErr

*     Check to determine if the values are within the range
*     of NetCDF data type
*     -----------------------------------------------------
      do 10, ival = 1, nval
         val = values ( ival )
         if ( val  .lt. I1_Min .or.
     .        val  .gt. I1_Max ) ierr = NCSysErr
 10   continue

*     If a value is out of range then return with an error message
*     ------------------------------------------------------------
      if ( ierr .ne. NCNoErr  ) then
         write ( stderr, 901) 'NCByte'
         return
      end if

*     Convert to 1-byte integers and store into work space
*     ----------------------------------------------------
      do 20, ival = 1, nval
         I1_Val ( ival ) = values ( ival ) 
 20   continue

*     Use NetCDF routine to write the values to file
*     ----------------------------------------------
      call NCVPT ( ncid,   varid, 
     .              start,  count,
     .              I1_Val, ierr )
      if ( ierr .ne. NCNoErr ) return

      return
*     ------

 901  format ( /, ' ODS_NCP_ItoI1: Value to be written as the ',
     .         /, '                NetCDF variable type, ', a,
     .         /, '                is out of range. ' )

      end
