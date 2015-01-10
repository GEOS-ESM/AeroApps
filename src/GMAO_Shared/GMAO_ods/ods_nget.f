
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_NGet --- Returns the number of observation reports
! 
! !DESCRIPTION: 
!    \label{ODS:NGet}
!     Returns the number of observation reports for a given
!     julian day and synoptic hour
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
! !INTERFACE:
!
      integer function ODS_NGet ( id, julian_day, syn_hour, ierr )
!
! !INPUT PARAMETERS:
      implicit   NONE
      integer    id                ! ODS file handle
      integer    julian_day        ! Julian day.  Use the function
                                   !   ODS_Julian to obtain this
                                   !   number.
      integer    syn_hour          ! hour of synoptic time since
                                   !   0:00 GMT (e.g., 0, 6, 12,
                                   !   18)
!
! !OUTPUT PARAMETERS:
!     integer    NObs              ! Number of observation report
      integer    ierr              ! Error code. If non-zero, an 
                                   !  error has occurred. For a list
                                   !  of possible values, see the
                                   !  description section of this
                                   !  prologue.
!
! !SEE ALSO: 
!     ODS_IGet()    gets the integer of a user-selected
!                   NetCDF file parameter
!     ODS_RGet()    gets the floating point value of a user-selected
!                   NetCDF file parameter
!     ODS_CGet()    gets the character string of a user-selected
!                   NetCDF file parameter
!     ODS_Julian()  converts the "calendar" date to the Julian day
!
! !REVISION HISTORY: 
!     13Apr1996   Redder   Original version.  Routine developed to
!                          create ODS version 2.00
!     20Apr1998   Redder   Correct bug in error message
!     19Nov1999   Redder   Added a latex label in and moved the
!                          subroutine statement into the prologue.
!                          Modified the comments for the return status
!                          code.
!     06Dec1999   Redder   Corrections to the documentation in the
!                          prologue.
!
!EOP
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Indices defining the set of observation data 
*     corresponding to a synoptic time in the ODS file
*     ------------------------------------------------
      integer  synbeg      ! starting index for the data block
      integer  synlen      ! number of observation reports in
                           !   the data block

*     Set output parameters to default
*     --------------------------------
      ODS_NGet = 0
      ierr     = NCNoErr

*     Check to determine if the file handle id is valid
*     -------------------------------------------------
      if ( id            .lt. 1      .or.
     .     id            .gt. id_max .or.
     .     IOMode ( id ) .eq. CLOSED ) then
         write ( stderr, 901 )
         ierr = NCEBadID
         return
      end if

*     Get pointer data from tables stored in common
*     ---------------------------------------------
      call ODS_GetP    ( id,        'write',
     .                   julian_day, syn_hour,
     .                   synbeg,     synlen,  ierr )
      if ( ierr .ne. NCNoErr ) return

*     The number of observation reports
*     be can extracted from pointer data
*     ----------------------------------
      ODS_NGet = synlen

      return
*     ------

 901  format ( /, ' ODS_NGet: File handle id number does not ',
     .         /, '           correspond to an opened ODS file' )

      end
