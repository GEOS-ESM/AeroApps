
*..............................................................


      subroutine ODS_UpdateP ( id,         update_option,
     .                         julian_day, syn_hour,
     .                         nval,       ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE: ODS_UpdateP
! 
! !DESCRIPTION: 
!     This routine performs the required modifications on a table of
!     pointers and related data when data for a specified julian day
!     and synoptic hour is being added (or deleted).  These tables
!     contain all pointers that are necessary to locate a block of
!     data within a file ( specified by the ODS file handle, id )
!     partitioned according to the julian day and synoptic hour.
!     The tables are stored in common as defined in the header file,
!     ods_hdf.h.  In addition, this routine updates NetCDF attribute
!     values also stored in common.
!
! !INTERFACE:
!      call ODS_UpdateP ( id,         update_option,
!                         julian_day, syn_hour,
!                         nval,       ierr )
!
! !INPUT PARAMETERS:
      integer           id             ! ODS file handle
      character * ( * ) update_option  ! = 'no_option' if none is
                                       !    desired.
                                       ! = 'append' to set up the
                                       !    pointers for appending
                                       !    data to a block
      integer           julian_day     ! Julian day
      integer           syn_hour       ! hours since 0:00 GMT    
      integer           nval           ! number of values corresponding
                                       !   to the day and synoptic
                                       !   period
!
! !OUTPUT PARAMETER:
      integer           ierr           ! return error code
!
!     note : The arguments julian_day and syn_hour are not accessed
!            if the update_option is set to 'append'
!
! !SEE ALSO:
!     ODS_ReSetP  ( Reset pointers and other relevant data )
!     ODS_ReadP   ( Read all pointer data from file. )
!     ODS_WriteP  ( Write all pointer data to file. )
!     ODS_Get P   ( Gets pointer data for reading/writing
!                   a block of data. )
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
!  25Apr1996   C. Redder   Origional Code
!  01Nov1999   C. Redder   Revised code to produce more stringent 
!                          comparison character strings.
!
!-------------------------------------------------------------------------

      include   'netcdf.inc'
      include   'ods_hdf.h'
      include   'ods_stdio.h'

*     Functions referenced
*     --------------------
      character * (20)
     .            ODS_Case        ! Set all character in string
                                  !   to upper or lower case
      integer     ODS_JulHr       ! Julian hour using, as the
                                  !   reference, the synoptic hour
                                  !   and julian day that corresponds
                                  !   to the last block of data
                                  !   written to file
      integer     ODS_StrSize     ! String length excluding trailing
                                  !   blanks

*     Other variables
*     ---------------
      integer     julday          ! temporary storage for the input
                                  !   argument, julian_day
      integer     synhour         ! temporary storage for the input
                                  !   argument, syn_hour
      integer     day             ! number of days from julian offset
      integer     syn             ! syn period of day
      integer     iday, isyn      ! index variables for do loop
      integer     delta_jhour     ! change of julian hour since the
                                  !   last block of data was written
                                  !   to file
      character * ( 20 )
     .            Option          ! temporary storage for the input
                                  !   argument, update_option
      integer     LOption         ! length of string

*     Default value for error code
*     ----------------------------
      ierr    = NCNoErr

*     Ensure that all characters in option are in lower case
*     ------------------------------------------------------
      Option  = ODS_Case    ( update_option, 'lower' )
      LOption = ODS_StrSize ( Option )

*     If option is 'append' ... 
*     -------------------------
      if ( LOption .eq. 6 .and.
     .     Option ( : 6 ) .eq. 'append' ) then

*        use the julian day and synoptic hour for
*        the last block of data written to file
*        ----------------------------------------
         julday  = latest_day  ( id )
         synhour = latest_hour ( id )
      else

*        use the input arguments
*        -----------------------
         julday  = julian_day
         synhour = syn_hour

      end if

*     Determine day and syn to be used as array indicies
*     --------------------------------------------------
      day      = julday  - julian_offset  ( id )
      syn      = synhour * nsyn ( id ) / hour_max + 1

*     If the table indicies are not within range specified
*     by the NetCDF file, return with an error message
*     ----------------------------------------------------
      if ( day .lt. 1 .or.
     .     day .gt. ndays ( id ) ) then

         write ( stderr, 901 ) julday
         ierr = NCEInVal
         return

      end if

      if ( syn .lt. 1 .or.
     .     syn .gt. nsyn  ( id ) ) then

         write ( stderr, 902 ) syn_hour
         ierr = NCEInVal
         return

      end if

*     Update the pointer table syn_beg for all
*     successive synoptic periods on the current day
*     ----------------------------------------------
      do 10, isyn = syn + 1, nsyn  ( id )
         syn_beg ( isyn, day,  id )
     .      = syn_beg ( isyn, day,   id ) + nval
 10   continue

*     Update the pointer table for all successive
*     days and synoptic periods
*     -------------------------------------------
      do 30, iday = day + 1, ndays ( id )
      do 20, isyn = 1,       nsyn  ( id )
         syn_beg ( isyn, iday, id )
     .      = syn_beg ( isyn, iday,  id ) + nval
 20   continue
 30   continue

*     Update the table, syn_len, only for current day
*     and synoptic period
*     -----------------------------------------------
      syn_len ( syn, day,id )
     .   = syn_len ( syn, day, id ) + nval

*     Determine the diffence of julian hour since the last block
*     of data was written to file
*     ----------------------------------------------------------
      delta_jhour = ODS_JulHr ( id, julday, synhour, ierr )
      if ( ierr .ne. NCNoErr ) return

*     If the hour has changed, then ...
*     ---------------------------------
      if ( delta_jhour .eq. 0 ) then

*        Set append mode to ON and set the
*        pointers to their appropriate values
*        ------------------------------------
         append_mode ( id ) = ON
         append_beg  ( id ) = syn_beg ( syn, day, id )
     .                      + syn_len ( syn, day, id ) - nval  
         append_len  ( id ) = nval

         else

*        Set the append mode to OFF and
*        ReSet the corresponding pointers
*        --------------------------------
         append_mode ( id ) = OFF
         append_beg  ( id ) = syn_beg ( syn, day, id )  
         append_len  ( id ) = syn_len ( syn, day, id )


      end if

*     Update the arrays containing the latest julian
*     day and hour
*     ----------------------------------------------
      latest_day     ( id ) = julday
      latest_hour    ( id ) = synhour

      return
*     ------

 901  format ( /, ' ODS_UpdateP : Invalid Julian day, ', i10 )
 902  format ( /, ' ODS_UpdateP : Invalid synoptic hour, ', i10 )

      end
