
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_Cal2Min() --- Convert "calendar" date and time to minutes since a given reference
!
! !DESCRIPTION:
!    \label{ODS:Cal2Min}
!     The routine converts date and time in calendar format to
!     minutes since a given reference date and time (also in
!     calendar format).  The format for the calendar date and
!     time is year 2000 compliant.
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
! !INTERFACE: 
!
      subroutine ODS_Cal2Min ( NVal, 
     .                         CalDate, CalTime,
     .                         RefDate, RefTime, Minutes )
!
!
! !INPUT PARAMETERS:
      implicit        NONE
      integer         NVal              ! number of values to be 
                                        !  converted
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
! !OUTPUT PARAMETERS:
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
! !SEE ALSO:
!     ODS_Min2Cal()   - Convert minutes since a given reference to
!                       "calendar" date and time
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
      integer     ODS_Julian

*     Other variables
*     ---------------
      integer     Time
      integer     Hour,    RHour
      integer     Minute,  RMinute
      integer              RMinutes
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

*        Determine the hour, minute and second of the day
*        ------------------------------------------------
         Time    = CalTime ( iVal )
         Hour    =       Time / 10000
         Minute  = mod ( Time,  10000 ) / 100
         Second  = mod ( Time,    100 )
         if ( Second .ge. 30 ) Minute = Minute + 1 

*        Determine the number of days from first_jday
*        --------------------------------------------
         days_offset
     .           = ODS_Julian ( CalDate ( iVal ) ) - RJDay

*        Determine the time in ODS foramt
*        --------------------------------
         Minutes ( iVal )
     .           = days_offset * MinDay 
     .           + Hour * 60
     .           + Minute
     .           - RMinutes

      end do

      return
      end
