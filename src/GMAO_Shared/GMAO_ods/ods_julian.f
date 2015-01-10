*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_Julian() --- Returns the Julian day
! 
! !DESCRIPTION:
!    \label{ODS:Julian}
!     The routine returns the Julian day number based on the
!     calendar date.  The algorithm was adopted from Press et al.
!     The Julian day number on May 23, 1968 is 2,440,000.  The
!     year zero is treated as the year 1 since the year 0 does
!     not exist.
!
!    \bigskip {\bf Reference:}
!    \begin{description}
!    \item Press, William H., Saul A. Teukolsky, William T.
!        Vetterling and Brian P. Flannery, 1992: {\em Numerical
!        Recipes in Fortran}, 2nd Ed. Cambridge University
!        Press, New York, NY, 963pp.
!    \end{description}
!
! !INTERFACE: 
!
      integer function ODS_Julian ( CalDate )
!
! !INPUT PARAMETERS:
      implicit NONE
      integer  CalDate  ! Calendar date in the format YYYYMMDD
                        !   where YYYY is the year, MM is the
                        !   month and DD is the day.  A negative
                        !   number implies that the year is B.C.
!
! !OUTPUT PARAMETERS:
!     integer  JulDay     Julian day number
!
! !SEE ALSO:
!     ODS_Min2Cal()   - Convert minutes since a given reference to
!                       "calendar" date and time
!     ODS_Cal2Min()   - Convert "calendar" date and time to minutes
!                       since a given reference
!     ODS_Time2Cal()  - converts ODS "time" attribute to "calendar"
!                       date and time
!     ODS_Cal2Time()  - converts "calendar" date and time to ODS
!                       time attribute
!     ODS_CalDat()    - calculates the "calendar" date and time
!                       from the Julian day
!
! !REVISION HISTORY: 
!     13Apr1998  C. Redder   Original code.  Routine was developed to
!                            create ODS version 2.00
!     19Nov1999  C. Redder   Added a latex label in and moved the
!                            subroutine statement into the prologue.
!     06Dec1999  C. Redder   Corrections to the documentation in the
!                            prologue.
!     10May2000  C. Redder   Updated prologue to include the routines
!                            ODS_Min2Cal and ODS_Cal2Min in the 
!                            "See Also" list
!
!EOP
!-------------------------------------------------------------------------

*     Other variables
*     ---------------
      integer     Year
      integer     Month
      integer     Day
      integer     iGreg  ! Gregorian Calendar adopted Oct 12, 1582
      parameter ( iGreg = 15 + 31 * ( 10 + 12 * 1582 ) )
      integer     JulDay
      integer     jy, jm, ja

      Year   =       CalDate / 10000
      Month  = mod ( CalDate,  10000 ) / 100
      Day    = mod ( CalDate,    100 )
 
*     Change year 0 to year 1
*     -----------------------
      if ( Year  .eq. 0 ) Year = 1

*     Account for the nonexisting year 0
*     ----------------------------------
      if ( Year  .lt. 0 ) Year = Year + 1

      if ( Month .gt. 2 ) then
         jy = Year
         jm = Month + 1

      else
         jy = Year  - 1
         jm = Month + 13

      endif

      JulDay = int ( 365.25  * jy )
     .       + int ( 30.6001 * jm )
     .       + Day + 1720995

*     Test whether to change to Gregorian Celendar
*     --------------------------------------------
      if ( Day + 31 * ( Month + 12 * Year ) .ge. iGreg) then
        ja     = int ( 0.01 * jy )
        Julday = JulDay + 2 - ja + int ( 0.25 * ja )

      endif

      ODS_Julian = JulDay

      return
      end
