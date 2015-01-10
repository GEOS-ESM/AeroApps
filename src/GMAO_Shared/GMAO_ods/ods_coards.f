
*..............................................................

      subroutine ODS_COARDS ( CalDate, CalTime, COARDS, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE: ODS_COARDS() --- Returns the time and date in COARDS
!                            format.
! 
! !DESCRIPTION:
!     The routine returns the time and date in COARDS format.
!     (e.g. 1992-10-08 15:15:42.5 ).  The routine checks the
!     length of the string to ensure that sufficient space
!     exists.  Trailing blanks are inserted if the string
!     length exceeds the required minimum.  The date begins
!     in the first character in the string.
!
! !INTERFACE: 
!     call ODS_COARDS ( CalDate, CalTime, COARDS, ierr )
!
! !INPUT PARAMETERS:
      integer           CalDate       ! Date in the format YYYYMMDD
                                      !  where YYYY is the year, MM is
                                      !  the month and DD is the day.
      integer           CalTime       ! Time in the format HHMMSS 
                                      !  where HH is the hour, MM is
                                      !  the minute and SS is the second.
! !OUTPUT PARAMETERS:
      character * ( * ) COARDS        ! Date and time in COARDS format
      integer           ierr          ! Error code. If non-zero, an error 
                                      !  has occurred. For a list of
                                      !  possible values see Table 8
                                      !  of da Silva and Redder (1995).

!
! !SEE ALSO:
!
! !REVISION HISTORY: 
!      3Apr1998   C. Redder   Original code.
!     16Feb2000   Todling     Rename stdio.h to ods_stdio.h
!
!-------------------------------------------------------------------------

      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Other variables
*     ---------------
      integer       Year
      integer       Month
      integer       Day
      integer       Hour
      integer       Minute
      integer       Second
      integer       LenStr
      integer       MinLStr
      parameter   ( MinLStr = 21 )
      character * ( MinLStr )
     .              Temp

*     Check the length of the string
*     ------------------------------
      LenStr = len ( COARDS )
      if ( LenStr .lt. MinLStr ) then
         ierr = ODS_DimErr
         write ( stderr, 901 )  LenStr, MinLStr

      end if

      Year   =       CalDate / 10000
      Month  = mod ( CalDate,  10000 ) / 100
      Day    = mod ( CalDate,    100 ) 
      Hour   =       CalTime / 10000
      Minute = mod ( CalTime,  10000 ) / 100
      Second = mod ( CalTime,    100 )

      write ( Temp, 801 ) Year, Month,  Day,
     .                    Hour, Minute, Second

      COARDS = Temp

      return

 801  format ( i4,   '-',
     .         i2.2, '-',
     .         i2.2, ' ',
     .         i2.2, ':',
     .         i2.2, ':',
     .         i2.2, '.0' )
 901  format ( /, ' ODS_COARDS: The length of the string ',
     .         /, '             supplied by the calling routine ',
     .         /, '             (= ', i4, ') is less than the ',
     .         /, '             required minimum (= ', i4, ')' )

      end
