
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_Min2Cal() --- Convert minutes since a given reference to "calendar" date and time
! 
! !DESCRIPTION:
!    \label{ODS:Min2Cal}
!     The routine converts the minutes since a given reference
!     date and time to calendar format.  The reference date and
!     time is also in calendar format.  The format for the
!     calendar date and time is year 2000 compliant.
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
! !INTERFACE: 
!
      subroutine ODS_Min2Cal ( NVal,    RefDate, RefTime,
     .                         Minutes, CalDate, CalTime )
!
! !INPUT PARAMETERS:
      implicit        NONE
      integer         NVal              ! number of values to be 
                                        !  converted
      integer         RefDate           ! Reference date in the format
                                        !  YYYYMMDD where YYYY is the
                                        !  year, MM is the month and
                                        !  DD is the day.  A negative
                                        !  number implies that the
                                        !  year is B.C.
      integer         RefTime           ! Reference time in the format
                                        !  HHMMSS where HH is the hour,
                                        !  MM is the minute and SS is
                                        !  the second.
      integer         Minutes  ( NVal ) ! Minutes since given reference
                                        !  date and time
!
! !OUTPUT PARAMETERS:
      integer         CalDate  ( NVal ) ! Calendar date in the format
                                        !  YYYYMMDD where YYYY is the
                                        !  year, MM is the month and
                                        !  DD is the day.  A negative
                                        !  number implies that the
                                        !  year is B.C.
      integer         CalTime  ( NVal ) ! Time in the format HHMMSS 
                                        !  where HH is the hour, MM
                                        !  is the minute and SS is
                                        !  the second.
!
! !SEE ALSO:
!     ODS_Cal2Min()   - Convert "calendar" date and time to minutes
!                       since a given reference
!     ODS_Time2Cal()  - converts ODS "time" attribute to "calendar"
!                       date and time
!     ODS_Cal2Time()  - converts "calendar" date and time to ODS
!                       time attribute
!     ODS_Julian()    - calculates the Julian day from date and
!                       time in "calendar" format
!     ODS_CalDat()    - calculates the "calendar" date and time
!                       from the Julian day
!
! !REVISION HISTORY: 
!     10May2000 C. Redder  Original code.
!
!EOP
!-------------------------------------------------------------------------

*     Function referenced
*     -------------------
      integer     ODS_CalDat
      integer     ODS_Julian

*     Other variables
*     ---------------
      integer     Hour,    RHour
      integer     Minute,  RMinute
      integer     Mins,    RMinutes
      integer     Second,  RSecond
      integer     RJDay,   days_offset
      integer     iVal
      integer     MinDay
      parameter ( MinDay = 60 * 24 )

*     Find the reference Julian minutes since reference Julian day
*     ------------------------------------------------------------ 
      RJDay    = ODS_Julian ( RefDate )
      RHour    =       RefTime / 10000
      RMinute  = mod ( RefTime,  10000 ) / 100
      RSecond  = mod ( RefTime,    100 )
      RMinutes = 60 *  RHour + RMinute
      if ( RSecond .ge. 30 ) RMinutes = RMinutes + 1

*     For each value ...
*     ------------------
      do iVal = 1, NVal

*        Determine the number of days from first_jday
*        --------------------------------------------
         Mins    = Minutes ( iVal ) + RMinutes
         days_offset
     .           = Mins / MinDay
         Mins    = mod ( Mins, MinDay )

*        Special handling for time prior to reference time
*        -------------------------------------------------
         if ( Mins .lt. 0 ) then
            days_offset = days_offset - 1
            Mins        = Mins        + MinDay

         end if

*        Determine the hour, minute and second of the day
*        ------------------------------------------------
         Hour    = mod ( Mins, MinDay ) / 60 
         Minute  = mod ( Mins, 60 )
         Second  = 0

*        Determine the calendar time and date
*        ------------------------------------
         CalTime ( iVal )
     .           = Hour   * 10000
     .           + Minute * 100
     .           + Second
         CalDate ( iVal )
     .           = ODS_CalDat ( RJDay + days_offset )

      end do

      return
      end
