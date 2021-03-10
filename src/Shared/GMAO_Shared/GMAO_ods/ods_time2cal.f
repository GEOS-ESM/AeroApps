
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_Time2Cal() --- Convert ODS to calendar date and time
! 
! !DESCRIPTION:
!    \label{ODS:Time2Cal}
!     The routine converts the ODS "time" to calendar format.
!     The ODS time is the number of minutes since 0:00 GMT of
!     the first Julian day (stored in common) on file.  The
!     Julian day in this routine is assumed to begin at 0:00
!     GMT.  The format for the calendar date and time is year
!     2000 compliant.
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
! !INTERFACE: 
!
      subroutine ODS_Time2Cal ( id, nval, ODSTime,
     .                          CalDate, CalTime, ierr )
!
! !INPUT PARAMETERS:
      implicit        NONE
      integer         id                ! ODS file handle as returned
                                        !  from ODS_Create() or
                                        !  ODS_Open().
      integer         nval              ! number of values to be 
                                        !  converted
      integer         ODSTime  ( nval ) ! ODS time
!
! !OUTPUT PARAMETERS:
      integer         CalDate  ( nval ) ! Calendar date in the format
                                        !  YYYYMMDD where YYYY is the
                                        !  year, MM is the month and
                                        !  DD is the day.  A negative
                                        !  number implies that the
                                        !  year is B.C.
      integer         CalTime  ( nval ) ! Time in the format HHMMSS 
                                        !  where HH is the hour, MM
                                        !  is the minute and SS is
                                        !  the second.
      integer         ierr              ! Error code. If non-zero, an 
                                        !  error has occurred. For a list
                                        !  of possible values, see the
                                        !  description section of this
                                        !  prologue.
!
! !SEE ALSO:
!     ODS_Min2Cal()   - Convert minutes since a given reference to
!                       "calendar" date and time
!     ODS_Cal2Min()   - Convert "calendar" date and time to minutes
!                       since a given reference
!     ODS_Cal2Time()  - converts "calendar" date and time to ODS
!                       time attribute
!     ODS_Julian()    - calculates the Julian day from date and
!                       time in "calendar" format
!     ODS_CalDat()    - calculates the "calendar" date and time
!                       from the Julian day
!
! !REVISION HISTORY: 
!     13Apr1998 C. Redder  Original code.  Routine was developed to
!                          create ODS version 2.00
!     19Nov1999 C. Redder  Added a latex label in and moved the
!                          subroutine statement into the prologue.
!                          Modified the comments for the return status
!                          code.
!     10May2000 C. Redder  Updated prologue to include the routines
!                          ODS_Min2Cal and ODS_Cal2Min in the 
!                          "See Also" list.  Fixed bug for handling
!                          negative values for the ODS time
!
!EOP
!-------------------------------------------------------------------------

      include    'netcdf.inc'
      include    'ods_hdf.h'
      include    'ods_stdio.h'

*     Function referenced
*     -------------------
      integer     ODS_CalDat

*     Other variables
*     ---------------
      integer     Hour
      integer     Minute, Minutes
      integer     Second
      integer     days_offset
      integer     ival
      integer     first_jday
      integer     MinDay
      parameter ( MinDay = 60 * 24 )

*     Set ierr code to valid 
*     ----------------------
      ierr    = NCNoErr

*     Check to determine if the file handle id is valid
*     -------------------------------------------------
      if ( id            .lt. 1      .or.
     .     id            .gt. id_max .or.
     .     IOMode ( id ) .eq. CLOSED ) then
         write ( stderr, 901 )
         ierr = NCEBadID
         return

      end if

*     Extract the first Julian day on file from common
*     (as declared in the header file ods_hdf.h)
*     ------------------------------------------------
      first_jday = julian_offset ( id ) + 1


*     For each value ...
*     ------------------
      do ival = 1, nval

*        Determine the number of days from first_jday
*        --------------------------------------------
         Minutes = ODSTime ( ival )
         days_offset    
     .           = Minutes / MinDay
         Minutes = mod ( Minutes, MinDay )

*        Special handling for negative values of ODS time
*        ------------------------------------------------
         if ( Minutes .lt. 0 ) then
            days_offset = days_offset - 1
            Minutes     = Minutes     + MinDay

         end if

*        Determine the hour, minute and second of the day
*        ------------------------------------------------
         Hour    = mod ( Minutes, MinDay ) / 60 
         Minute  = mod ( Minutes, 60 )
         Second  = 0

*        Determine the calendar time and date
*        ------------------------------------
         CalTime ( ival )
     .           = Hour   * 10000
     .           + Minute * 100
     .           + Second
         CalDate ( ival )
     .           = ODS_CalDat ( first_jday + days_offset )

      end do


 901  format ( /, ' ODS_Time2Cal: File handle id number does ',
     .         /, '               not correspond to an opened ',
     .         /, '               ODS file' )

      return
      end
