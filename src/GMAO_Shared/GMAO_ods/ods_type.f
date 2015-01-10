
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:      ODS_Type
!
! !DESCRIPTION: 
!
!     Returns the ODS type.  The valid types are 'pre-analysis' 
!     and 'post-analysis'
!
      character * (*) function ODS_Type ( id, ierr )
!
! !INPUT PARAMETERS:
      implicit   NONE
      integer    id      ! ODS file handle
      integer    ierr    ! return status code
!
! !REVISION HISTORY: 
!     24Sep98   Redder   Origional version
!     15Mar99   Redder   Incorporated routine into ODS library
!
!-------------------------------------------------------------------------

      include    'netcdf.inc'
      include    'ods_hdf.h'
      include    'ods_stdio.h'

*     NetCDF file id
*     --------------
      integer     nc_id

*     NetCDF variable
*     ---------------
      integer     varid       ! variable id

*     Temporary storage
*     -----------------
      integer     AttLen      ! returned attribute length
      integer     ierr_omf    ! return status code for inquiring
                              !   about the variable, 'omf'
      integer     ierr_oma    ! return status code fpr inquiring
                              !   about the variable, 'oma'
      integer     ncopts      ! NetCDF error handling options

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

*     Extract the NetCDF file id
*     --------------------------
      nc_id = ncid ( id )

*     Save current options of error handling and turn off
*     error messages and set errors to be non-fatal.
*     ---------------------------------------------------
      call NCGOPT     ( ncopts )
      call NCPOPT     ( 0 )

*     Retrieve ODS type
*     -----------------      
      call ODS_NCAGTC ( nc_id,    NCGLOBAL,
     .                 'type',    AttLen,
     .                  ODS_Type, ierr )

*     Return error handling option to their previous values
*     -----------------------------------------------------
      call NCPOPT     ( ncopts )

*     Check status codes and ...
*     --------------------------
      if      ( ierr .eq. NCNoErr  ) then
         return                ! ... return if successful
                               ! ------------------------

      else if ( ierr .eq. NCENoAtt ) then
         ierr = NCNoErr        ! ... continue and determine the type by
                               ! checking the ODS variables, oma and omf,
                               ! the the attribute does not exist
                               ! --------------------------------------

      else if ( ierr .eq. NCSysErr ) then
         ODS_Type = ' '        ! ... return with no error message since 
         return                ! a message should already have been 
                               ! printed by the routine, ODS_NCAGTC
                               ! --------------------------------------

      else if ( ierr .eq. NCEBadTy ) then
         ODS_Type = ' '        ! ... return with no error message since
         return                ! a message should already have been
                               ! printed by the routine, ODS_NCAGTC
                               ! --------------------------------------

      else
         ODS_Type = ' '        ! ... return with error message
         write ( stderr, 902 ) ! -----------------------------
         return

      end if

*     Save current options of error handling and turn off
*     error messages and set errors to be non-fatal.
*     ---------------------------------------------------
      call NCGOPT     ( ncopts )
      call NCPOPT     ( 0 )

*     Inquire about the variables, omf and oma.
*     -----------------------------------------
      varid = NCVID   ( nc_id,  'omf', ierr_omf )
      varid = NCVID   ( nc_id,  'oma', ierr_oma )

*     Return error handling option to their previous values
*     -----------------------------------------------------
      call NCPOPT     ( ncopts )

*     Check return status codes for ...
*     ---------------------------------
      if      ( ierr_omf .ne. NCNoErr   .and.
     .          ierr_omf .ne. NCENotVr ) then
         ODS_Type  = ' '             ! ... inquiring about 'omf'
         ierr      =  ierr_omf       ! -------------------------
         write ( stderr, 903 ) 'omf'
         return

      else if ( ierr_oma .ne. NCNoErr   .and.
     .          ierr_oma .ne. NCENotVr ) then
         ODS_Type  = ' '             ! ... inquiring about 'oma'
         ierr      =  ierr_oma       ! -------------------------
         write ( stderr, 903 ) 'oma'
         return

      else if ( ierr_omf .ne. ierr_oma ) then
         ODS_Type  = ' '             ! ... inconsistent status codes
         ierr      =  NCSysErr       ! -----------------------------
         write ( stderr, 904 ) ierr_omf, ierr_oma
         return

      end if

*     If the both variables exists ...
*     --------------------------------
      if      ( ierr_omf .eq. NCNoErr  ) then
          ODS_Type = 'post_analysis' ! ... then the file is post-analysis
                                     ! ----------------------------------

      else if ( ierr_omf .eq. NCENotVr ) then
          ODS_Type = 'pre_analysis'  ! ... else the file is pre-analysis
                                     ! if the error code is NCENotVr
                                     ! ---------------------------------

      end if

 901  format ( /, ' ODS_Type: File handle id number does not ',
     .         /, '           correspond to an opened ODS file' )
 902  format ( /, ' ODS_Type: Error in retrieving the ODS  ',
     .                       'global attribute, tupe.' )
 903  format ( /, ' ODS_Type: Error in inquiring about the ',
     .                       'variable, ', a, '.' )
 904  format ( /, ' ODS_Type: Inconsistent status codes from the ',
     .         /, '           the NetCDF routine for the ',
     .         /, '           variables, omf and oma. ',
     .         /, '           Status codes for the variable ... ',
     .         /, '           ... omf = ', i5, 
     .         /, '           ... oma = ', i5 )


      return
      end 
