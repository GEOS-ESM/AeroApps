
*..............................................................


      character * ( * ) function ODS_NameMonth ( IMonth, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:      ODS_NameMonth
! 
! !DESCRIPTION: 
!     Returns the name of the month as specified by its integer
!     representation.  The first letter of the name is in upper
!     case, and the remaining characters are in lower case.
!
! !INTERFACE: NameMonth = ODS_NameMonth ( IMonth, ierr )
!
! !INPUT PARAMETER:
      integer   IMonth  ! The integer representation of the month

! !OUTPUT PARAMETER:
      integer   ierr    ! Returned error code
!
! !SEE ALSO:
!     ODS_IMonth ( Returns the integer representation of the month
!                  as specified by the first three characters of 
!                  the name. )
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers

! !REVISION HISTORY: 
!     14Feb96   Redder   Origional version
!     16Feb2000 Todling  Rename stdio.h to ods_stdio.h
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_stdio.h'

*     Table for the names of all months
*     ---------------------------------
      character * ( 9 ) Months ( 12 )

      data  Months ( 1 )  / 'January' /,
     .      Months ( 2 )  / 'February' /,
     .      Months ( 3 )  / 'March' /,
     .      Months ( 4 )  / 'April' /,
     .      Months ( 5 )  / 'May' /,
     .      Months ( 6 )  / 'June' /,
     .      Months ( 7 )  / 'July' /,
     .      Months ( 8 )  / 'August' /,
     .      Months ( 9 )  / 'September' /,
     .      Months ( 10 ) / 'October' /,
     .      Months ( 11 ) / 'November' /,
     .      Months ( 12 ) / 'December' /

      ierr  = NCNoErr

*     if integer value is valid ...
*     -----------------------------
      if ( IMonth .ge. 1 .and.
     .     IMonth .le. 12 ) then

*        specify the name of the month
*        -----------------------------
         ODS_NameMonth = Months ( IMonth )

      else ! print error message and return with a invalid error code 
           ! --------------------------------------------------------

         ODS_NameMonth = ' '
         ierr          = NCEInVal
         write ( stderr, 901 ) IMonth
      end if

      return

 901  format ( /, ' ODS_NameMonth: Invalid input integer for month ',
     .         /, '                Input integer is: ', i7 ) 

      end
