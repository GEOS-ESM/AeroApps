

*..............................................................


      integer function ODS_JulHr ( id, julian_day, syn_hour, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE: ODS_JulHr
! 
! !DESCRIPTION: 
!     Calculates the Julian hour using, as a reference, the synoptic
!     hour and julian day that correspond to the last block of data
!     written to file.  For example, if the Julian day is 5001 and
!     the synoptic hour is 12 and the reference day and hour is 5000
!     and 6, then the routine returns a value of 30.
!
! !INTERFACE:
!     JulHr = ODS_JulHr ( id, julian_day, syn_hour, ierr )
!
! !INPUT PARAMETERS:
      integer           id         ! ODS file handle
      integer           julian_day ! Julian day
      integer           syn_hour   ! hours since 0:00 GMT    
!
! !OUTPUT PARAMETER:
      integer           ierr     ! return error code
!
! !SEE ALSO:
!     ODS_UpdateP ( Updates the pointer data )
!     ODS_ReSetP  ( Reset pointers and other relevant data )
!     ODS_ReadP   ( Read all pointer data from file. )
!     ODS_WriteP  ( Write all pointer data to file. )
!     ODS_Get P   ( Gets pointer data for reading/writing
!                   a block of data. )
!     ODS_GetAttP ( Gets the attribute of the pointers )
!
! !FILES USED:  
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!
! !REVISION HISTORY:
!
!  11Apr96   C. Redder   Origional Code
!
!-------------------------------------------------------------------------
      include 'netcdf.inc'
      include 'ods_hdf.h'

      ODS_JulHr =   julian_day        * 24 + syn_hour
     .          - ( latest_day ( id ) * 24 + latest_hour ( id ) )

      return
      end
