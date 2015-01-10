

*....................................................................


      subroutine ODS_Info ( id,        filename,   mode, version,
     .                      FirstJDay, LatestJDay, LatestHour,
     .                      n_kt,      kt_names,   kt_units,
     .                      n_kx,      kx_names,   kx_meta,
     .                      n_qcx,     qcx_names,  ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:  ODS_Info() --- Returns information about an opened
!                           ODS file
! 
! !DESCRIPTION:
!     Returns the NetCDF and file status parameters
!
! !INTERFACE: 
!     call ODS_Info  ( id,        filename,   mode, version,
!                      FirstJDay, LatestJDay, LatestHour
!                      n_kt,      kt_names,   kt_units,
!                      n_kx,      kx_names,   kx_meta,
!                      n_qcx,     qcx_names,
!                      ierr )
!
!
! !INPUT PARAMETERS:
      integer          id             ! ODS file handle (file id)
                                      !   use this to refer to this
                                      !   file later on.

! !INPUT/OUTPUT PARAMETERS:
      integer          n_kt           ! on input: the maximum value
                                      !   allowed for this integer
                                      !   (usually determined by
                                      !   the space allocated for
                                      !   storage)
                                      ! on output: for maximum
                                      !   number of GEOS/DAS data
     .                                !   types
      integer          n_kx           ! on input: the maximum value
                                      !   allowed for this integer
                                      !   (usually determined by
                                      !   the space allocated for
                                      !   storage)
                                      ! on output: for maximum
                                      !   number of GEOS/DAS data
     .                                !   sources
      integer          n_qcx          ! on input: the maximum value
                                      !   allowed for this integer
                                      !   (usually determined by
                                      !   the space allocated for
                                      !   storage)
                                      ! on output: for maximum
                                      !   number of possible 
                                      !   values for the quality
                                      !   control exclusion mark
                                      !   (qcexcl)
                                      ! note: If the input values
                                      !   for either nkt, kx or
                                      !   nqcx is smaller than
                                      !   the number of values to
                                      !   be read, then the routine
                                      !   exits with a system error
                                      !   message and a code of
                                      !   ODS_DimErr ( = -3 ).
                                      !   Also, nkt, nkx and/or
                                      !   nqcx are set to the
                                      !   minimum value required
                                      !   in order for the routine
                                      !   to execute successfully.
!
! !OUTPUT PARAMETERS:
      character * (*)  filename       ! name of ODS file
      character * (*)  mode           ! mode = 'w' open for writing
                                      ! mode = 'r' open for reading 
      character * (*)  version        ! version tag
      integer          FirstJDay      ! first julian day on file
      integer          LatestJDay     ! latest julian day for which
                                      !   data exists
      integer          LatestHour     ! latest julian hour for which
                                      !   data exists
      character * (*)  kt_names  ( n_kt )
                                      ! name for each data type
      character * (*)  kt_units  ( n_kt )
                                      ! units for each data type
      character * (*)  kx_names  ( n_kx )
                                      ! name for each data source
      character * (*)  kx_meta   ( n_kx )
                                      ! kx specific metadata information
                                      !   Use this to specify the
                                      !   meaning of the parameter "xm"
                                      !   or to specify the name of the
                                      !   external OMS file name.
                                      !   e.g. "oms_file:myinstr.oms"
      character * (*)  qcx_names ( n_qcx ) 
                                      ! information about the meaning
                                      !   of the variable, "qcexcl".
      integer          ierr           ! Error code. If non-zero, an 
                                      !   error has occurred. For a list
                                      !   of possible values see Table 8
                                      !   of da Silva and Redder (1995).
                                      !   If an error has occurred, then
                                      !   file was closed.
!
!     NOTE: No more than id_max ( as defined in header file ods_hdf.h )
!           files can be opened at any one time.  Use the function
!           ODS_Julian to obtain the input parameters, FirstJDay, and
!           LatestJDay.
!
! !SEE ALSO:
!     ODS_Create()   creates the ODS file, sets the dimensions
!                    and saves the text data for kt_names, kt_units
!                    and kx_names. 
!
!
! !REVISION HISTORY:
!     09Apr1998   Redder    Original code.
!     16Feb2000   Todling   Rename stdio.h to ods_stdio.h
!
!-------------------------------------------------------------------------

      include  'ods_hdf.h'
      include  'netcdf.inc'
      include  'ods_stdio.h'

*     Function referenced
*     -------------------
      integer   ODS_VerIndex  ! determines the ODS version
                              !  index number

*     Other variables
*     ---------------
      integer   str_len       ! string length for each name
     .                        !  and unit
      integer   ierr_temp     ! temporary storage for the
                              !  returned error code
      integer   version_index ! version index number
      integer   version_2     ! index number for version 2.00

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

*     Get ODS file name
*     -----------------
      call ODS_CGet    ( id, ':filename',   filename,    ierr )
      if ( ierr .ne. NCNoErr ) return

*     Get mode (read or write)
*     ------------------------
      call ODS_CGet    ( id, ':mode',       mode,        ierr )
      if ( ierr .ne. NCNoErr ) return

*     Get ODS version number
*     ----------------------
      call ODS_CGet    ( id, ':version',    version,     ierr )
      if ( ierr .ne. NCNoErr ) return
      version_index = ODS_VerIndex ( version, ierr )

*     Get some NetCDF file dimensions and parameters from ODS
*     -------------------------------------------------------
      version_2     = ODS_VerIndex ( '2.00', ierr )
      if ( ierr .ne. NCNoErr ) return

      if ( version_index .ge. version_2 ) then
         call ODS_IGet    ( id, 'nkt',      n_kt,        ierr )
         if ( ierr .ne. NCNoErr ) return
         call ODS_IGet    ( id, 'nkx',      n_kx,        ierr )
         if ( ierr .ne. NCNoErr ) return
         call ODS_IGet    ( id, 'nqcx',     n_qcx,       ierr )
         if ( ierr .ne. NCNoErr ) return

      else
         call ODS_IGet    ( id, 'ktmax',    n_kt,        ierr )
         if ( ierr .ne. NCNoErr ) return
         call ODS_IGet    ( id, 'kxmax',    n_kx,        ierr )
         if ( ierr .ne. NCNoErr ) return

      end if

      call ODS_IGet    ( id, 'strlen',      str_len,     ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_IGet    ( id, 'syn_beg:first_julian_day',  
     .                                      FirstJDay,   ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_IGet    ( id, 'syn_beg:latest_julian_day',
     .                                      LatestJDay,  ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_IGet    ( id, 'syn_beg:latest_synoptic_hour',
     .                                      LatestHour,  ierr )
      if ( ierr .ne. NCNoErr ) return

*     Read the text labels for data types and sources currently
*     available and quality control codes currently implemented
*     ---------------------------------------------------------
      call ODS_GetList ( id, 'kt_names', n_kt, kt_names, ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_GetList ( id, 'kt_units', n_kt, kt_units, ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_GetList ( id, 'kx_names', n_kx, kx_names, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Read these tables only if the file corre-
*     sponds to one of the later versions
*     -----------------------------------------
      if ( version_index .ge. version_2 ) then
         call ODS_GetList
     .               ( id, 'kx_meta',   n_kx,  kx_meta,   ierr )
         if ( ierr .ne. NCNoErr ) return
         call ODS_GetList
     .               ( id, 'qcx_names', n_qcx, qcx_names, ierr )
         if ( ierr .ne. NCNoErr ) return

      end if

      return
*     ------

 901  format ( /, ' ODS_Info: File handle id number does not ',
     .         /, '           correspond to an opened ODS file' )

      end
