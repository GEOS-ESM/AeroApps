
*...................................................................


      subroutine ODS_NCAGTC ( ncid,   varid,
     .                        AttNam, AttLen,
     .                        string, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_NCAGTC
! 
! !DESCRIPTION: 
!     This routine calls the NetCDF routines, NCAINQ and 
!     NCAGTG to extract the NetCDF attribute character
!     string for a user specified variable id and attribute
!     name.
!
! !INTERFACE:  ODS_NCAGTC ( ncid,  varid,
!                                  AttNam, AttLen,
!                                  string, ierr )
!
! !INPUT PARAMETERS:
      integer            ncid    ! NetCDF file id
      integer            varid   ! NetCDF variable id
      character * ( * )  AttNam  ! Name of attribute string
!
! !INPUT/OUTPUT PARAMETERS:
      integer            AttLen  ! Length of attribute character string
      character * ( * )  string  ! Attribute character string
      integer            ierr    ! Returned error code
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!
! !LIBRARIES ACCESSED:
!     NetCDF
!
! !REVISION HISTORY: 
!     05Jun96   C. Redder   Origional version
!     15Mar99   C. Redder   Fixed bug in error handling.  The routine
!                           now checks the status code returned from
!                           ncainq before checking the scratch space.
!     16Feb2000 R. Todling  Rename stdio.h to ods_stdio.h
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_stdio.h'

*     Other variables
*     ---------------
      integer            AtType  ! attribute type
      integer            strlen  ! length of String

*     Default value for error code
*     ----------------------------
      ierr    = NCNoErr

*     Get information about attribute
*     -------------------------------
      call ncainq ( ncid, varid, AttNam, AtType, AttLen, ierr )

*     If the error code is not valid, return.
*     ---------------------------------------
      if ( ierr .ne. NCNoErr ) return

*     Return and print error message if space is insufficient.
*     --------------------------------------------------------
      strlen = len ( String )
      if ( AttLen .gt. strlen ) then
         write ( stderr, * )
         write ( stderr, * )
     .     ' ODS_NCAGTC: the space allocated for the string',  
     .                 ' attribute is '
         write ( stderr, * )
     .     '             insufficient.  Increase max_strlen',
     .                 ' to at least ',
     .                   AttLen
                          
         ierr = NCSysErr
         return

         end if

*     If the attribute is not a character string return
*     -------------------------------------------------
      if ( AtType .ne. NCChar ) then

         write ( stderr, * )
         write ( stderr, * )
     .      ' ODS_NCAGTC: Attribute is not a character string. '

         ierr = NCEBadTy
         return

      end if

*     Get Attribute
*     -------------
      call ncagtc ( ncid, varid, AttNam, String, AttLen, ierr )
      String = String ( : AttLen )   ! Insure that any unused portion of 
     .                               !   of String is filled with blanks

      return
      end
