
*...................................................................


      subroutine ODS_NCG_R4toI ( ncid,  varid,
     .                           start, count,
     .                           nval,  values, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_NCG_R4toI
! 
! !DESCRIPTION:
!     This routine calls the NetCDF routine, NCVGT, to read four
!     byte integer data from an ODS/HDF file and converts the 
!     data to native integers.  This routine insures that the
!     data are correctly returned from the NetCDF routine,
!     irrespective of the machine on which the software is running. 
!
! !INTERFACE:  ODS_NCG_R4toI ( ncid,  varid,
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
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!     ods_worksp.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures for work space
!
! !SEE ALSO:
!     ODS_NCG_ItoR4  ( The inverse companion of this routine )
!
! !REVISION HISTORY: !
!     14May96    C. Redder   Origional code
!     16Feb2000  R. Todling  Rename stdio.h to ods_stdio.h
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_stdio.h'
      include 'ods_worksp.h'

*     Other variables
*     ---------------
      double precision
     .             val     ! Temporary storage
      integer     ival     ! Index variable for do loop

*     Default value for error code
*     ----------------------------
      ierr    = NCNoErr

*     Use NetCDF routine to read the values from the file
*     and save the values in the work space
*     ---------------------------------------------------
      call NCVGT  ( ncid,   varid, 
     .              start,  count,
     .              R4_Val, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Default value for error code
*     ----------------------------
      ierr    = NCNoErr

*     Check to determine if the values are within the range
*     of the native format
*     -----------------------------------------------------
      do 10, ival = 1, nval
         val      = dble ( R4_Val ( ival ) )
         if ( val  .lt. I_MinD .or.
     .        val  .gt. I_MaxD ) ierr = NCSysErr
 10   continue

*     If a value is out of range then return with an error message
*     ------------------------------------------------------------
      if ( ierr .ne. NCNoErr  ) then
         write ( stderr, 901) 'NCFloat'
         return
      end if

*     Convert to native format
*     ------------------------
      do 20, ival = 1, nval
         values ( ival ) = nint ( R4_Val ( ival ) )
 20   continue

      return
*     ------

 901  format ( /, ' ODS_NCG_R4toI: A value read as the ',
     .         /, '                NetCDF variable type, ', a,
     .         /, '                and to be converted to a ',
     .         /, '                native integer is out of ',
     .         /, '                range. ' )

      end
