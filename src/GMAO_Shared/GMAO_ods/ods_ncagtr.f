

*...................................................................


      subroutine ODS_NCAGTR ( ncid,   varid,
     .                        AttNam, AttLen,
     .                        values, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_NCAGTR
! 
! !DESCRIPTION: 
!     This routine calls the NetCDF routines, NCAINQ and 
!     NCAGT to extract the NetCDF attribute values for
!     a user specified variable id and attribute name.  The
!     attribute values are return as native real numbers.
!     This routine was design to extract a limited number
!     of values as specified by the variable, Max_WorkSp,
!     as defined in the header file, ods_worksp.h
!
! !INTERFACE:  ODS_NCAGTR ( ncid,   varid,
!                                  AttNam, AttLen,
!                                  values, ierr )
!
! !INPUT PARAMETERS:
      integer            ncid          ! NetCDF file id
      integer            varid         ! NetCDF variable id
      character * ( * )  AttNam        ! Name of attribute values
!
! !INPUT/OUTPUT PARAMETERS:
      integer            AttLen        ! Number of attribute values
      real               values  ( * ) ! Attribute values
      integer            ierr          ! Returned error code
!
!     note: This routine treats each byte as unsigned integers
!           with values within the range from 0 to 255, inclusively.
!           An offset value of 256 is added to any negative value.
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
! !LIBRARIES ACCESSED:
!     NetCDF
!
! !REVISION HISTORY: 
!     05Jun1996   C. Redder   Origional version
!     16Feb2000   R. Todling  Rename stdio.h to ods_stdio.h
!     13Sep2002   C. Redder   Changes made to handle one-byte integers
!                             stored in a long character string 
!     27Feb2009   D. Nadeau   Fixed variable type mismatch (char/int1)
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_stdio.h'
      include 'ods_worksp.h'

*     one more than the maximum value for
*     any integer in unsigned byte format
*     -----------------------------------
      integer            I1Max
      parameter        ( I1Max      = 256 )

*     scratch space for ...
*     ----------------------
      integer            itemp   ! native integers

*     Other variables
*     ---------------
      integer            AtType  ! attribute type
      integer            iAtt    ! index variable for do loop

*     Default value for error code
*     ----------------------------
      ierr    = NCNoErr

*     Get information about attribute
*     -------------------------------
      call ncainq ( ncid, varid, AttNam, AtType, AttLen, ierr )

*     If the error code is not valid, return.
*     ---------------------------------------
      if ( ierr .ne. NCNoErr ) return

*     If the number of attribute values is too large ( i.e. insuf-
*     ficient scatch space), then print message and return
*     ------------------------------------------------------------
      if ( AttLen .gt. Max_WorkSp ) then

         write ( stderr, 901 ) AttLen
         ierr = NCSysErr
         return

      end if

*     Extract the attribute values and convert to real numbers
*     if the type of values as stored in the NetCDF file are 
*     --------------------------------------------------------

*     one byte integers
*     -----------------
      if      ( AtType .eq. NCByte  ) then

         call ncagt  ( ncid, varid, AttNam, I1_Val, ierr )

         do 10, iAtt = 1, AttLen
            itemp           = I1_Val ( iAtt )
            if ( itemp .lt. 0 ) itemp = itemp + i1max
            values ( iAtt ) = real ( itemp )
 10      continue

*     two byte integers
*     -----------------
      else if ( AtType .eq. NCShort ) then

         call ncagt  ( ncid, varid, AttNam, I2_Val, ierr )

         do 20, iAtt = 1, AttLen
            values ( iAtt ) = real ( I2_Val ( iAtt ) )
 20      continue

*     four byte integers
*     ------------------
      else if ( AtType .eq. NCLong  ) then

         call ncagt  ( ncid, varid, AttNam, I4_Val, ierr )

         do 30, iAtt = 1, AttLen
            values ( iAtt ) = real ( I4_Val ( iAtt ) )
 30      continue

*     four byte real numbers
*     ----------------------
      else if ( AtType .eq. NCFloat ) then

         call ncagt  ( ncid, varid, AttNam, R4_Val, ierr )

         do 40, iAtt = 1, AttLen
            values ( iAtt ) =        R4_Val ( iAtt )
 40      continue

*     Return with an error message if the attribute type
*     is not recognized as valid.
*     --------------------------------------------------
      else 

         write ( stderr, 902 )
         ierr = NCEBadTy
         return

      end if

      return
*     ------

 901  format ( /, ' ODS_NCAGTR: Attribute length exceeds the ',
     .         /, '             maximum allowed.  Increase ',
     .         /, '             Max_WorkSp to at least ', i9 )
 902  format ( /, ' ODS_NCAGTR: Not an appropriate attribute ',
     .         /, '             type for this routine. ' )

      end
