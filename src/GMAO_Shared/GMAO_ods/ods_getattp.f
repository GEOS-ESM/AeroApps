
*..............................................................


      subroutine ODS_GetAttP ( id, AttName, AttVal, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:      ODS_GetAttP
! 
! !DESCRIPTION: 
!     Returns an attribute of a pointer in a ODS/HDF file.  The
!     attribute is identified by the character string, AttName.
!     The current valid names are latest_julian_day and
!     latest_synoptic_hour.  The calling routine attempts to
!     extract an invalid parameter, then then the routine returns
!     an error code of NCEInval (as defined in the header file,
!     netcdf.inc) but does not print an error message.  Otherwise
!     the error code is set to NCNoErr (also defined in the header
!     file).  This routine is especially designed for internal use
!     in routines such as ODS_IGet.
!
! !INTERFACE: call ODS_GetAttP  ( id, AttName, AttVal, ierr )
!
! !INPUT PARAMETERS:
      integer         id      ! ODS file handle
      character * (*) AttName ! The name of the attribute to be 
                              !   obtained.  The case of each
                              !   letter is significant.
!
! !OUTPUT PARAMETERS: 
      integer         AttVal  ! The integer value of the attribute
      integer         ierr    ! The return error code
!
! !SEE ALSO: 
!     ODS_UpdateP ( Updates the pointer data )
!     ODS_ReSetP  ( Reset pointers and other relevant data )
!     ODS_ReadP   ( Read all pointer data from file. )
!     ODS_WriteP  ( Write all pointer data to file. )
!     ODS_Get P   ( Gets pointer data for reading/writing
!                   a block of data. )
!     ODS_JulHr   ( Returns the number of hours from the 
!                   Julian hour of the last block of data
!                   written to file )
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!
! !REVISION HISTORY: 
!     15Apr1996   Redder   Origional version
!     01Nov1999   Redder   Revised code to prevent subscript errors
!                          in character strings
!     16Feb2000   Todling  Rename stdio.h to ods_stdio.h
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Function referenced
*     -------------------
      integer  ODS_StrSize  ! String length excluding trailing blanks
                            ! 

*     Local variable
*     --------------
      integer  LAttName     ! String length of attribute name

      ierr     = NCNoErr
      LAttName = ODS_StrSize ( AttName )

      if ( LAttName .eq. 25 ) then
      if ( AttName ( 1:LAttName ) .eq.
     .         'syn_beg:latest_julian_day'    ) then
         AttVal = latest_day   ( id )
         return
      end if
      end if

      if ( LAttName .eq. 28 ) then
      if ( AttName ( 1:LAttName ) .eq.
     .         'syn_beg:latest_synoptic_hour' ) then
         AttVal = latest_hour  ( id )
         return

      end if
      end if

      ierr   = NCEInval

      return
      end
