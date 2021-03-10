

*..............................................................


      integer function ODS_IMonth ( Month_Name, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:      ODS_IMonth
! 
! !DESCRIPTION: 
!     Returns the integer representation of the month as specified
!     by the first three characters of the name.
!
! !INTERFACE: IMonth = ODS_IMonth  ( Month_Name, ierr )
!
! !INPUT PARAMETER:
      character * (*) Month_Name ! The name of the month as determined
                                 !   by the first three characters.
                                 !   These characters can be either
                                 !   upper or lower case.  The remaining
                                 !   characters are ignored by the
                                 !   routine

! !OUTPUT PARAMETER:
      integer         ierr       ! Returned error code
!
! !SEE ALSO:
!     ODS_NameMonth ( Returns the name of the month as specified
!                     by its integer representation )
! 
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers

! !REVISION HISTORY: 
!     14Feb1996   Redder   Origional version
!     01Nov1999   Redder   Revised code to prevent subscript errors
!                          in character strings
!     16Feb2000   Todling  Rename stdio.h to ods_stdio.h
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_stdio.h'

      character * ( 3 ) Month
      character * ( 3 ) ODS_Case

      ierr  = NCNoErr

*     specify the month using only the first
*     three characters in upper case
*     ---------------------------------------
      Month = Month_Name
      Month = ODS_Case ( Month, 'upper' )

*     search for the appropriate name and
*     specify its integer representation
*     -----------------------------------
      if      ( Month .eq. 'JAN' ) then
         ODS_IMonth = 1

      else if ( Month .eq. 'FEB' ) then
         ODS_IMonth = 2

      else if ( Month .eq. 'MAR' ) then
         ODS_IMonth = 3

      else if ( Month .eq. 'APR' ) then
         ODS_IMonth = 4

      else if ( Month .eq. 'MAY' ) then
         ODS_IMonth = 5

      else if ( Month .eq. 'JUN' ) then
         ODS_IMonth = 6

      else if ( Month .eq. 'JUL' ) then
         ODS_IMonth = 7

      else if ( Month .eq. 'AUG' ) then
         ODS_IMonth = 8

      else if ( Month .eq. 'SEP' ) then
         ODS_IMonth = 9

      else if ( Month .eq. 'OCT' ) then
         ODS_IMonth = 10

      else if ( Month .eq. 'NOV' ) then
         ODS_IMonth = 11

      else if ( Month .eq. 'DEC' ) then
         ODS_IMonth = 12

      else ! print error message and return with a invalid error code 
           ! --------------------------------------------------------
         ODS_IMonth = 0
         ierr       = NCEInVal
         write ( stderr, 901 ) Month_Name
      end if

      return

 901  format ( /, ' ODS_IMonth: Invalid input name for month ',
     .         /, '             Input name is: ', a ) 

      end
