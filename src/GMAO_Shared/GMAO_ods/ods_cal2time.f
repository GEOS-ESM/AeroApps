
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_Cal2Time() --- Convert calendar to ODS "time" format
! 
! !DESCRIPTION:
!    \label{ODS:Cal2Time}
!     The routine converts data and time in "calendar" format
!     to the ODS "time" format.  The ODS time is the number of
!     minutes since 0:00 GMT of the first Julian day (stored in
!     common) on file.  The Julian day in this routine is assumed
!     to begin at 0:00 GMT. The format for the calendar date and
!     time must be year 2000 compliant.
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
! !INTERFACE: 
!
      subroutine ODS_Cal2Time ( id, nval, CalDate, CalTime,
     .                          ODSTime, ierr )
!
!
! !INPUT PARAMETERS:
      implicit        NONE
      integer         id                ! ODS file handle as returned
                                        !  from ODS_Create() or
                                        !  ODS_Open().
      integer         nval              ! number of values to be 
                                        !  converted
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
!
! !OUTPUT PARAMETERS:
      integer         ODSTime  ( nval ) ! ODS time
      integer         ierr              ! Error code. If non-zero, an
                                        !  error has occurred. For a
                                        !  list of possible values
                                        !  see Table 8 of da Silva
                                        !  and Redder (1995).
!
! !SEE ALSO:
!     ODS_Min2Cal()   - Convert minutes since a given reference to
!                       "calendar" date and time
!     ODS_Cal2Min()   - Convert "calendar" date and time to minutes
!                       since a given reference
!     ODS_Time2Cal()  - converts ODS "time" attribute to "calendar"
!                       date and time
!     ODS_Julian()    - calculates the Julian day from date and
!                       time in "calendar" format
!     ODS_CalDat()    - calculates the "calendar" date and time
!                       from the Julian day
!
! !REVISION HISTORY: 
!     13Apr1998  C. Redder  Original code.  Routine was developed to
!                           create ODS version 2.00
!     19Nov1999  C. Redder  Added a latex label in and moved the
!                           subroutine statement into the prologue.
!                           Modified the comments for the return status
!                           code.
!     06Dec1999  C. Redder  Minor change in case of variable name.
!     16Feb2000  R. Todling Rename stdio.h to ods_stdio.h
!     10May2000  C. Redder  Updated prologue to include the routines
!                           ODS_Min2Cal and ODS_Cal2Min in the 
!                           "See Also" list
!
!EOP
!-------------------------------------------------------------------------

      include    'netcdf.inc'
      include    'ods_hdf.h'
      include    'ods_stdio.h'

*     Function referenced
*     -------------------
      integer     ODS_Julian

*     Other variables
*     ---------------
      integer     Time
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

*        Determine the hour, minute and second of the day
*        ------------------------------------------------
         Time    = CalTime ( ival )
         Hour    = Time / 10000
         Minute  = mod ( Time,  10000 ) / 100
         Second  = mod ( Time,    100 )
         if ( Second .ge. 30 ) Minute = Minute + 1 

*        Determine the number of days from first_jday
*        --------------------------------------------
         days_offset
     .           = ODS_Julian ( CalDate ( ival ) ) - first_jday

*        Determine the time in ODS foramt
*        --------------------------------
         ODSTime ( ival )
     .           = days_offset * MinDay + Hour * 60 + Minute

      end do


 901  format ( /, ' ODS_Cal2Time: File handle id number does ',
     .         /, '               not correspond to an opened ',
     .         /, '               ODS file' )

      return
      end
