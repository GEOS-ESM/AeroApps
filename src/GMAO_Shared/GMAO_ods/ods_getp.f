

*..............................................................


      subroutine ODS_GetP   ( id,         rwmode, 
     .                        julian_day, syn_hour,
     .                        start,      count,   ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE: ODS_GetP
!
! !DESCRIPTION:
!     Gets pointer data for NetCDF read and write routines based
!     on day and synoptic day.  These tables contain all pointers
!     that are necessary to locate a block of data within a file
!     partitioned according to the julian day and synoptic hour.
!     The pointers are extracted from tables stored in common.
!
! !INTERFACE: ODS_GetP ( id,    rwmode, 
!                               julian_day, time,
!                               start,      count, ierr )
!
! !INPUT PARAMETERS:
      integer           id         ! ODS file handle
      character * ( * ) rwmode     ! = 'write' if the calling
                                   !    routine is writing data
                                   !    to file
                                   ! = 'read' if the calling
                                   !    routine is read data from
                                   !    file
                                   !  (case insensitive)
      integer           julian_day ! Julian day
      integer           syn_hour   ! minutes since 0:00 GMT    
!
! ! OUTPUT PARAMETERS:
      integer           start      ! NetCDF index for first observation
      integer           count      ! number of observation for a day,
                                   !    and synptic period
      integer           ierr       ! return error code
!
! !SEE ALSO:
!     ODS_UpdateP ( Updates the pointer data )
!     ODS_ReSetP  ( Reset pointers and other relevant data )
!     ODS_ReadP   ( Read all pointer data from file. )
!     ODS_WriteP  ( Write all pointer data to file. )
!     ODS_JulHr   ( Returns the number of hours from the 
!                   Julian hour of the last block of data
!                   written to file )
!     ODS_GetAttP ( Gets the attribute of the pointers )
!
! !FILES USED:  
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!
! !REVISION HISTORY: 
!
!  15Apr1996   C. Redder   Origional Code
!  01Nov1999   C. Redder   Revised code to produce more stringent 
!                          comparison character strings.
!  28Nov1999   C. Redder   Fixed bug in determining the length of
!                          the character string, rw_mode.
!  16Feb2000   R. Todling  Rename stdio.h to ods_stdio.h
!
!-------------------------------------------------------------------------

      include   'netcdf.inc'
      include   'ods_hdf.h'
      include   'ods_stdio.h'

*     Functions referenced
*     --------------------
      character   ODS_Case   * ( 5 ) ! set string to upper or
                                     !   lower case
      integer     ODS_JulHr          ! determine the Julian hour
                                     !   from to the Julian hour
                                     !   of the last data block
                                     !   written to file
      integer     ODS_StrSize        ! String length excluding trailing
                                     !   blanks

*     Other variables
*     ---------------
      integer     day                ! number of days from julian
                                     !   offset
      integer     syn                ! syn period of day
      character   rw_mode    * ( 6 ) ! temporary storage for input
                                     !   argument rwmode
      integer     LMode              ! String length of rw_mode
      integer     JulHr              ! temporary storage for
                                     !    returned from the
                                     !    function, ODS_JulHr
      integer     JulHr_Latest       ! JulHr for the last data
                                     !    block written to file
      parameter ( JulHr_Latest = 0 ) ! Set to 0 for clarity and
                                     !    convenience

*     Default value for error code
*     ----------------------------
      ierr     = NCNoErr

*     Determine day and hour to be used as table indicies
*     --------------------------------------------------
      day      = julian_day - julian_offset  ( id )
      syn      = syn_hour * nsyn ( id ) / hour_max + 1

*     If the table indicies are not within range specified
*     by the NetCDF file, return with an error message
*     ---------------------------------------------------
      if ( day .lt. 1 .or.
     .     day .gt. ndays ( id ) ) then

         write ( stderr, 901 ) julian_day
         ierr = NCEInVal
         return

      end if

      if ( syn .lt. 1 .or.
     .     syn .gt. nsyn  ( id ) ) then

         write ( stderr, 902 ) syn_hour
         ierr = NCEInVal
         return

      end if

*     Get the pointer data from the tables
*     ------------------------------------
      start = syn_beg ( syn, day, id )
      count = syn_len ( syn, day, id )

*     If  append mode is ON then set start and count using
*     the pointer values for the last data block written to file.
*     -----------------------------------------------------------
      if ( append_mode ( id ) .eq. ON ) then
         rw_mode = ODS_Case    ( rwmode, 'lower' )
         LMode   = ODS_StrSize ( rw_mode )
         JulHr   = ODS_JulHr   ( id, julian_day, syn_hour, ierr ) 

*        This section of code is applicable only if the 
*        calling routine is writing data ( i.e rwmode is set
*        to 'write' ) and the data is being added to the last
*        block written to file ( i.e. JulHr = JulHr_Latest )
*        ---------------------------------------------------
         if ( JulHr           .eq.  JulHr_Latest .and.
     .        LMode           .eq.  5            .and.
     .        rw_mode ( : 5 ) .eq. 'write' )  then
            start = append_beg ( id )
            count = append_len ( id )
         end if
      end if

      return
*     ------

 901  format ( /, ' ODS_GetP : Invalid Julian day, ', i10 )
 902  format ( /, ' ODS_GetP : Invalid synoptic hour, ', i10 )


      end
