
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_CalDat() --- Convert Julian day to calendar date
! 
! !DESCRIPTION:
!    \label{ODS:CalDat}
!     The routine returns the calendar day based on the Julian
!     day number.  The algorithm was adopted from Press et al.
!     The Julian day number on May 23, 1968 is 2,440,000.  The
!     return date is in the format YYYYMMDD where YYYY is the
!     year, MM is the month and DD is the day.  A negative
!     number implies that the year is B.C.
!
!    \bigskip {\bf Reference:}
!    \begin{description}
!    \item Press, William H., Saul A. Teukolsky, William T.
!          Vetterling and Brian P. Flannery, 1992: {\em Numerical
!          Recipes in Fortran}, 2nd Ed. Cambridge University
!          Press, New York, NY, 963pp.
!    \end{description}
!
! !INTERFACE: 
!
      integer function ODS_CalDat ( JulDay )
!
! !INPUT PARAMETERS:
      implicit NONE
      integer  JulDay   ! Julian day number
!
! !OUTPUT PARAMETERS:
!     integer  CalDat    Date in "calendar" format
!
! !SEE ALSO:
!     ODS_Min2Cal()  - Convert minutes since a given reference to
!                      "calendar" date and time
!     ODS_Cal2Min()  - Convert "calendar" date and time to minutes
!                      since a given reference
!     ODS_Time2Cal() - converts ODS "time" attribute to "calendar"
!                      date and time
!     ODS_Cal2Time() - converts "calendar" date and time to ODS
!                      time attribute
!     ODS_Julian()   - calculates the Julian day from date and
!                      time in "calendar" format
!
! !REVISION HISTORY: 
!     13Apr1998  C. Redder  Original code.  Routine was developed to
!                           create ODS version 2.00
!     19Nov1999  C. Redder  Added a latex label in and moved the
!                           subroutine statement into the prologue.
!     10May2000  C. Redder  Updated prologue to include the routines
!                           ODS_Min2Cal and ODS_Cal2Min in the 
!                           "See Also" list
!
!EOP
!-------------------------------------------------------------------------

*     Other variables
*     ---------------
      integer       Year
      integer       Month
      integer       Day
      integer      iGreg  ! The Julian day number of the Gregorian
      parameter  ( iGreg = 2299161 )  ! Calendar adopted Oct 12, 1582
      integer       Alpha
      integer       ja, jb, jc, jd, je

*     Cross-over to Gregorian Calendar produces this correction
*     ---------------------------------------------------------
      if ( JulDay .ge. iGreg ) then
         Alpha = int ((( JulDay - 1867216 ) - 0.25 ) / 36524.25 )
         ja    = JulDay + 1 + Alpha - int ( 0.25 * Alpha )

      else ! no correction
*     --------------------
         ja    = JulDay

      endif

      jb    = ja + 1524
      jc    = int ( 6680. + (( jb - 2439870 ) - 122.1 )
     .      / 365.25)
      jd    = 365 * jc + int ( 0.25 * jc )
      je    = int (( jb - jd ) / 30.6001 )

      Day   = jb - jd - int ( 30.6001 * je )

      Month = je - 1
      if ( Month .gt. 12 ) Month = Month - 12

      Year  = jc - 4715
      if ( Month .gt. 2 ) Year = Year - 1
      if ( Year  .le. 0 ) Year = Year - 1

      ODS_CalDat = Year * 10000 + Month * 100 + Day 

      return
      end
