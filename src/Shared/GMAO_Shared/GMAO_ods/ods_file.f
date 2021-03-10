* ods_hdf.f - last update: Nov 03, 1999
*
*                  OBSERVATION DATA STREAM VERSION 2.03
*
* This file contains routines implementing DAO's Observation Data
* Stream (ODS), Version 2.00. A description of this package is given in
*
*  da Silva, A. and C. Redder, 1995: Documentation of the GEOS/DAS 
*        Observation Data Stream (ODS) Version 2.00,  DAO Office Note 95-01,
*        NASA Goddard Space Flight Center, Greenbelt, MD, 30pp.
*
*  On-line versions of this document are available from
*
*  ftp://dao.gsfc.nasa.gov/pub/office_notes/on9501.ps.Z  (postscript)
*
*  ftp://niteroi.gsfc.nasa.gov/www/on9501/ods.html       (HTML)
*
*....................................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  ODS_Create() --- Creates an ODS file
! 
! !DESCRIPTION:
!    \label{ODS:Create}
!     Creates an ODS file (including setting up the data structure
!     and contents), sets some internal file parameters and saves
!     text information associated with GEOS/DAS standard data types,
!     data sources and quality control exclusion marks.
!
!     Two types of ODS can be created:
!
! \begin{description}
!
! \item[Pre-analysis: ] this file contains the observed value, along
!   with space-time coordinates and a meta-data indicator. At DAO
!   this file is usually produced by the REPACK pre-processing system.
!
! \item[Post-analysis: ] In addition to the same information present
!  in the pre-analysis files, additional information produced during
!  the assimilation process is included. Examples are the difference
!  between 6 hour forecast and the observed values, and quality
!  control flags assigned to the observations.
!
! \end{description}
!
!  Note: For a list of error codes, see Table~\ref{tab:errors}.
!  
!
! !INTERFACE: 
      subroutine ODS_Create ( id,    filename, ods_type,
     .                        FirstJDay,
     .                        n_kt,  kt_names, kt_units,
     .                        n_kx,  kx_names, kx_meta,
     .                        n_qcx, qcx_names,
     .                        ierr )

!
! !INPUT PARAMETERS:
      implicit         NONE
      character * (*)  filename             ! name of ODS file
      character * (*)  ods_type             ! = 'pre_anal', if a pre-analysis 
                                            !   ODS file is to be created.
                                            ! = 'post_anal', if a post-analysis
                                            !   ODS file is to be created.
      integer          FirstJDay            ! first julian day to be written
                                            !  to the file.  Use the
                                            !  function ODS_Julian to
                                            !  obtain this number.
      integer          n_kt                 ! maximum number of GEOS/DAS 
                                            !   data types
      character * (*)  kt_names  ( n_kt )   ! names for each data type
      character * (*)  kt_units  ( n_kt )   ! units for each data type
      integer          n_kx                 ! maximum number of GEOS/DAS
                                            !   data sources
      character * (*)  kx_names  ( n_kx )   ! names fore each data source
      character * (*)  kx_meta   ( n_kx )   ! kx specific meta-data
                                            !   information Use this to specify
                                            !   the meaning of the parameter
                                            !   "xm" or to specify the name of
                                            !    the external OMS file name.
                                            !   e.g. "oms_file:myinstr.oms"
      integer          n_qcx                ! maximum number of possible 
                                            !   values for the quality
                                            !   control exclusion mark
                                            !   (qcexcl)
      character * (*)  qcx_names ( n_qcx )  ! information about the meaning
                                            !   of the variable, "qcexcl".
!
! !OUTPUT PARAMETERS:
      integer          id                   ! ODS file handle (or file id);
                                            !  use this to refer to this
                                            !  file later on.
      integer          ierr                 ! Error code. If non-zero, an 
                                            !  error has occurred. For a list
                                            !  of possible values, see the
                                            !  description section of this
                                            !  prologue.  If an error has
                                            !  occurred, then file was either
                                            !  closed or not created
!
!     NOTE: No more than id_max ( as defined in header file ods_hdf.h )
!           files can be opened at any one time.
!
! !SEE ALSO:
!     ODS_Open()    opens a pre-existing ODS file and returns the
!                   dimensions and text data for kt_names, kt_units
!                   and kx_names.
!     ODS_Close()   closes the ODS file
!     ODS_Julian()  convert the "calendar" date to the Julian day
!
! !REVISION HISTORY:
!     02Oct1995   da Silva  Specification.
!     17May1996   Redder    Original code.
!     13Apr1998   Redder    Changed the variable names from kt_max and
!                           kx_max to nkt and nkx.  Added code to
!                           process according to version number.
!                           Added code to retrieve the tables kx_meta
!                           and qcx_names and netcdf dimension variable
!                           nqcx.  These changes were made to create ODS
!                           version 2.00
!     19Nov1999   Redder    Added a latex label in and moved the
!                           subroutine statement into the prologue.
!                           Modified the comments for the return status
!                           code.
!     06Dec1999   Redder    Corrections to the documentation in the
!                           prologue.
!
!EOP
!-------------------------------------------------------------------------

*     Other variables
*     ---------------
      integer   str_len       ! string length of each text
                              !   label for data types and 
                              !   source
      integer   ierr_temp     ! temporary storage for the
                              !   returned error code

      include 'netcdf.inc'

*     Set ierr code to valid 
*     ----------------------
      ierr    = NCNoErr

*     Determine the NetCDF dimension, str_len
*     ---------------------------------------
      str_len =       len ( kt_names  (1) )
      str_len = max ( len ( kt_units  (1) ), str_len )
      str_len = max ( len ( kx_names  (1) ), str_len )
      str_len = max ( len ( kx_meta   (1) ), str_len )
      str_len = max ( len ( qcx_names (1) ), str_len )

*     Set some NetCDF file dimensions for ODS
*     ---------------------------------------
      call ODS_SetParmI ( 'nkt',     n_kt,     ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_SetParmI ( 'nkx',     n_kx,     ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_SetParmI ( 'nqcx',    n_qcx,    ierr )
      if ( ierr .ne. NCNoErr ) return
      call ODS_SetParmI ( 'strlen',  str_len,  ierr )
      if ( ierr .ne. NCNoErr ) return

*     Create the NetCDF file for ODS
*     ------------------------------
      call ODS_NCCreate ( id, filename, ods_type, FirstJDay, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Write the text labels for data types and
*     sources currently available
*     ----------------------------------------
      call ODS_PutList  ( id, 'kt_names', n_kt, kt_names,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      call ODS_PutList  ( id, 'kt_units', n_kt, kt_units,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      call ODS_PutList  ( id, 'kx_names', n_kx, kx_names,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      call ODS_PutList  ( id, 'kx_meta',  n_kx, kx_meta,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      call ODS_PutList  ( id, 'qcx_names',
     .                                   n_qcx, qcx_names,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

      return
*     ------

*     Clean up by closing the file and
*     setting io mode to CLOSED
*     --------------------------------
 801  continue
      call ODS_Close    ( id, ' ', ierr_temp )

      return
*     ------

      end

*....................................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  ODS_Open() --- Opens an ODS file
! 
! !DESCRIPTION:
!    \label{ODS:Open}
!     Opens an existing ODS file for reading or writing. 
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
! !INTERFACE: 
      subroutine ODS_Open ( id, filename, mode, ierr )
!
! !INPUT PARAMETERS:
      implicit         NONE
      character * (*)  filename       ! name of ODS file
      character * (*)  mode           ! mode = 'w' open for writing
                                      ! mode = 'r' open for reading 
!
! !OUTPUT PARAMETERS:
      integer          id             ! ODS file handle (file id)
                                      !   use this to refer to this
                                      !   file later on.
      integer          ierr           ! Error code. If non-zero, an 
                                      !  error has occurred. For a list
                                      !  of possible values, see the
                                      !  description section of this
                                      !  prologue.  If an error has
                                      !  occurred, then file was closed.
!
!     NOTE: No more than id_max ( as defined in header file ods_hdf.h )
!           files can be opened at any one time.
!
! !SEE ALSO:
!     ODS_Create()   creates the ODS file, sets the dimensions
!                    and saves the text data for kt_names, kt_units,
!                    kx_names, kx_meta, and qcx_names. 
!
!
! !REVISION HISTORY:
!     02Oct1995   da Silva  Specification.
!     17May1996   Redder    Original code.
!     13Apr1998   Redder    Changed the variable names from ktmax and
!                           kxmax to nkt and nkx.  Added code to process
!                           according to version number.  Added code to
!                           retrieve the tables kx_meta and qcx_names
!                           and netcdf dimension variable nqcx.  Removed
!                           input and output arguments to create a
!                           minimalistic version of this routine.  These
!                           changes were made to create ODS version 2.00
!     19Nov1999   Redder    Added a latex label in and moved the
!                           subroutine statement into the prologue.
!                           Modified the comments for the return status
!                           code.
!
!EOP
!-------------------------------------------------------------------------

      include  'ods_hdf.h'
      include  'netcdf.inc'

*     Set ierr code to valid 
*     ----------------------
      ierr    = NCNoErr

*     Open the NetCDF file for ODS
*     ----------------------------
      call ODS_NCOpen  ( id, filename, mode, ierr )
      if ( ierr .ne. NCNoErr ) return

      return
*     ------
      end

*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_Append() --- Sets up software to append data to file
! 
! !DESCRIPTION:
!    \label{ODS:Append}
!     Sets up the software package to add additional observation
!     reports to the data segment specified by the latest julian
!     day and synoptic time for which data was written ( stored
!     in common ).  The set up does not affect the subsequent
!     execution of the routines, ODS\_GetI and ODS\_GetR.
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
! !INTERFACE: 
      subroutine ODS_Append ( id, nval, ierr )
!
! !INPUT PARAMETERS:
      implicit        NONE
      integer         id              ! ODS file handle as returned
                                      !  from ODS_Create() or
                                      !  ODS_Open().
                                      !   
      integer         nval            ! the number of values to be
                                      !   added.
! !OUTPUT PARAMETERS:
      integer          ierr           ! Error code. If non-zero, an 
                                      !  error has occurred. For a list
                                      !  of possible values, see the
                                      !  description section of this
                                      !  prologue.

!
! !SEE ALSO:
!     ODS_PutC()    write a list of character strings to   an ODS file
!     ODS_GetC()    get   a list of character strings from an ODS file
!     ODS_PutI()    write a list of integer values    to   an ODS file
!     ODS_GetI()    get   a list of integer values    from an ODS file
!     ODS_PutR()    write a list of real    values    to   an ODS file
!     ODS_GetR()    get   a list of real    values    to   an ODS file
!     ODS_Julian()  convert the "calendar" date to the Julian day
!     ods_hdf.h     include file defining internal constants
!                    and global variables, and setting up data
!                    structures
!
! !REVISION HISTORY: 
!     02Oct1995   da Silva    Specification.
!     17May1996   C. Redder   Original code.
!     13Apr1998   C. Redder   Updated documentation to reflection
!                             changes made in other routines to
!                             create ODS version 2.00
!     19Nov1999   C. Redder   Added a latex label in and moved the
!                             subroutine statement into the prologue.
!                             Modified the comments for the return status
!                             code.
!     06Dec1999   C. Redder   Corrections to the documentation in the
!                             prologue.
!     16Feb2000   Todling     Rename stdio.h to ods_stdio.h
!     16Sep2002   C. Redder   Updated prologue
!
!EOP
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Other variables
*     ---------------
      integer  julian_day ! julian_day
      integer  syn_hour   ! hour of synoptic time

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

*     If mode set IO mode to read only,
*     then exit with an error message
*     --------------------------------
      if ( IOMode ( id ) .ne. NCWRITE ) then
         write ( stderr, 902 ) filenames ( id )
         ierr = NCEPerm
         return

      end if

*     Set variables as filler arguments for the routine,
*     ODS_UpdateP ( These variables will not be accessed ) 
*     ----------------------------------------------------
      julian_day  = 0
      syn_hour    = 0

*     Update pointer tables and variables ( tables stored
*     in common ) to be used to append the data.
*     ---------------------------------------------------
      call ODS_UpdateP ( id,        'append',
     .                   julian_day, syn_hour,
     .                   nval,       ierr )
      if ( ierr .ne. NCNoErr ) return

      return

 901  format ( /, ' ODS_Append: File handle id number does not ',
     .         /, '             correspond to an opened ODS file' )
 902  format ( /, ' ODS_Append: Cannot append to file in read ',
     .         /, '             only mode.  Filename is: ', a )

      end

*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:           ODS_Close() --- Closes an ODS file
! 
! !DESCRIPTION: 
!    \label{ODS:Close}
!     Closes an ODS file and saves the pointer data necessary to
!     observation data stored in segments according Julian day
!     and synoptic time. If the file has been opened for writing
!     the event tag and the updated pointer data are saved
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
! !INTERFACE:  
!
      subroutine ODS_Close ( id, event, ierr )
!
! !INPUT PARAMETERS:
      implicit          NONE
      integer           id         ! ODS file handle as returned from
                                   !  ODS_Create() or ODS_Open().
      character * ( * ) event      ! program event tag - usually the
                                   !  name of the program which modified
                                   !  the file and a date.  The
                                   !  length of the string should be
                                   !  as small as possible with no
                                   !  trailing blanks
!
! !OUTPUT PARAMETERS: 
      integer           ierr       ! Error code. If non-zero, an 
                                   !  error has occurred. For a list
                                   !  of possible values, see the
                                   !  description section of this
                                   !  prologue.  If an error has
                                   !  occurred, then the pointer data
                                   !  and thus the observation data will
                                   !  likely not be updated.
!      
! !SEE ALSO:
!     ODS_Create()   creates an ODS file 
!     ODS_Open()     opens an existing ODS file
!     ods_hdf.h      an include file defining internal constants
!                     and global variables, and setting up data
!                     structures
!     ods_stdio.h        an include file defining standard input/output
!                     unit numbers
!
! !SYSTEM ROUTINES:
!     The NetCDF API is used to access the HDF/ODS file.
!
! !REVISION HISTORY:
!     02Oct1995   da Silva    Specification.
!     17May1996   Redder      Original code.
!     07Jul1998   Redder      Fixed bug that truncates the last
!                             character in the history tag when an
!                             event tag is appended.
!     19Nov1999   Redder      Added a latex label in and moved the
!                             subroutine statement into the prologue.
!                             Modified the comments for the return
!                             status code.
!     06Dec1999   Redder      Corrections to the documentation in the
!                             prologue.
!     16Feb2000   Todling     Rename stdio.h to ods_stdio.h
!
!EOP
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Functions referenced
*     --------------------
      integer     ODS_StrSize    ! index variable for last non-blank
                                 !   character

*     Other variables
*     ---------------
      integer     nc_id          ! temporary storage for NetCDF file id
      integer     event_len      ! length of event tag
      integer     ohist_len      ! length of old version of history tag
      integer     hist_len       ! length of updated history tag
      integer     att_len        ! length of attribute
      integer     varid          ! variable id
      integer     NameSz         ! number of characters in the filename
      integer     ierr_temp      ! temporary storage for returned error
                                 !   code

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

*     Get NetCDF id
*     -------------
      nc_id   = ncid ( id )

*     if io mode set for writing ...
*     ------------------------------
      if ( IOMode ( id ) .eq. NCWRITE ) then

*        Get file history tag
*        --------------------
         call ODS_NCAGTC ( nc_id, NCGLOBAL,
     .                           'history', att_len,
     .                            history,  ierr )
         if ( ierr .ne. NCNoErr ) go to 10  ! Clean up


*        Determine the lengths of all tags
*        ---------------------------------
         ohist_len = ODS_StrSize ( history )
         event_len = ODS_StrSize ( event )
         hist_len  = ohist_len + event_len
     .             + 2                  ! account the event
                                        !   separator, ";"
                                        !   and blank
  
*        If necessary, adjust the lengths, hist_len and event_len.
*        Any adjustments imply that the event tag will be truncated
*        to fit into the space allocated for the history tag.
*        ----------------------------------------------------------
         if ( hist_len .gt. max_strlen .or. 
     .        hist_len .gt. att_len ) then

            hist_len  = min ( hist_len, max_strlen )
            hist_len  = min ( hist_len, att_len )

*           print warning message
*           ---------------------
            NameSz    = index ( filenames ( id ), ' ' ) - 1
            if ( NameSz .le. 0 )
     .         NameSz = len   ( filenames ( id ) )

            write ( stderr, 902 ) filenames ( id ) ( : NameSz )

         end if

*        If event tag is not empty ...
*        -----------------------------
         if    ( event_len .ge. 1 ) then

*           If history tag is already empty, initialize the tag
*           ---------------------------------------------------
            if ( ohist_len .lt. 1 ) then
               history =  event   ( : event_len )

            else ! update the history tag
*           -----------------------------
               history =  history ( : ohist_len ) // '; '
     .                 // event   ( : event_len )
            end if

*           Save history tag
*           ----------------
            call ncaptc ( nc_id, NCGLOBAL, 'history', NCCHAR,   
     .                           att_len,
     .                           history ( : hist_len ),
     .                           ierr )
            if ( ierr .ne. NCNoErr ) go to 10  ! Clean up

         end if

*        Save pointers stored in tables and in common 
*        --------------------------------------------
         call ODS_WriteP ( id, ierr )
         if ( ierr .ne. NCNoErr ) go to 10  ! Clean up

      end if

 10   continue

*     Close NetCDF data file
*     ----------------------
      call ncclos ( nc_id, ierr_temp )
      if ( ierr .eq. NCNoErr ) ierr  = ierr_temp

*     IOMode is set to CLOSED
*     -----------------------
      IOMode ( id ) = CLOSED

      return
*     ------

 901  format ( /, ' ODS_Close: File handle id number does not ',
     .         /, '            correspond to an opened ODS file' )
 902  format ( /, ' ODS_Close: The event tag was truncated to ',
     .         /, '            fit into the available allocated ',
     .         /, '            space for the history tag',
     .         /, '            filename is: ',a )

      end

*..............................................................

      block data ODS_VerList

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:     ODS_VerList ( data block )
! 
! !DESCRIPTION:
!     Initializes version list and stores the list in 
!     common; defined in header file ods_hdf.h)
!
! !INTERFACE: non executable
!
! !FILES USED:
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!
! !REVISION HISTORY:
!     06Apr1098   Redder   Origional version
!     26May1998   Redder   Added version 2.02 to list
!     02Nov1999   Redder   Added version 2.03 to list
!     19Nov1999   Redder   Added version 2.10 to list
!     12Dec2000   Redder   Added version 2.11-2.15 to list
!     16Feb2001   Todling  Added version 2.16
!     16Feb2001   Redder   Added version 2.17
!     24Sep2002   Redder   Added version 2.18
!
!-------------------------------------------------------------------------

      include  'ods_hdf.h'

      data NVersions           /  NVer  /
      data version_list (  1 ) / '1.01' /,
     .     version_list (  2 ) / '1.02' /,
     .     version_list (  3 ) / '2.00' /,
     .     version_list (  4 ) / '2.01' /,
     .     version_list (  5 ) / '2.02' /,
     .     version_list (  6 ) / '2.03' /,
     .     version_list (  7 ) / '2.04' /,
     .     version_list (  8 ) / '2.10' /,
     .     version_list (  9 ) / '2.11' /,
     .     version_list ( 10 ) / '2.12' /,
     .     version_list ( 11 ) / '2.13' /,
     .     version_list ( 12 ) / '2.14' /,
     .     version_list ( 13 ) / '2.15' /,
     .     version_list ( 14 ) / '2.16' /,
     .     version_list ( 15 ) / '2.17' /,
     .     version_list ( 16 ) / '2.18' /

      end 

*..............................................................


      block data ODS_IOStat 

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:     ODS_IOStat ( data block )
! 
! !DESCRIPTION:
!     Initializes integers specifying input / output mode (stored in
!     common; defined in header file ods_hdf.h)
!
! !INTERFACE: non executable
!
! !FILES USED:
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!
! !REVISION HISTORY:
!     13Oct95   Redder   Origional version
!
!-------------------------------------------------------------------------

      include 'ods_hdf.h'

*     Other variables
*     ---------------
      integer   id   ! index variable for data statement

      data ( IOMode ( id ), id = 1, id_max ) / id_max * CLOSED /

      end


*..............................................................


      block data ODS__Parm

      implicit NONE


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:     ODS__Parm ( data block )
!
! !DESCRIPTION:
!     Initializes variables used to set the NetCDF file parameters.
!
! !INTERFACE: non executable
!
! !FILES USED:
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!
! !REVISION HISTORY:
!     10Apr96   Redder   Origional version
!
!-------------------------------------------------------------------------

      include 'ods_hdf.h'

      data   ParmI_nlist / CLEAR /
      data   ParmC_nlist / CLEAR /

      end

*..............................................................
