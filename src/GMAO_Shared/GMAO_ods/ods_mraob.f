!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 610.1, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:          MRAOB_Create - Create OMS file
! 
! !DESCRIPTION:      
!     Creates an OMS file and initializes the pointers
!     specifying the location and size of data blocks that
!     are partitioned according to julian day and synoptic
!     hour.  
!
! !INTERFACE:
      
      subroutine MRAOB_Create ( id, filename, first_jday, ierr )

      implicit NONE

! !INPUT PARAMETERS:
      character * (*)  filename             ! name of ODS file
      integer          first_jday           ! first julian day
!
! !OUTPUT PARAMETERS: 
      integer          id                   ! ODS file handle
      integer          ierr                 ! error code
!
! !REVISION HISTORY:
!    07Apr2002  da Silva    Based on ODS_NCCreate.
!    23Apr2002  Lucchesi    Added support for integer station ID
!    24Sep2002  C. Redder   Corrected string length in email-address,
!                           and the number of dimensions for the NetCDF
!                           variable, istatn.  Added dimension for
!                           string length of station ID and added the
!                           variables for the station ID (Id), station
!                           elevation (elev) and launch time (ltime).
!                           Added Julian day at reference date attribute
!                           for the variable, syn_beg.  Fixed bug in
!                           defining the attribute, 
!                           syn_beg:value_at_reference_date
!    10Oct2002  C. Redder   Added code to define the lists, it_names
!                           (names of WMO instrument types), and rc_descr
!                           (desription of WMO radiation correction codes).
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

!     dimension ids for ...
!     ---------------------
      integer     nbatchesdim ! number of batches of
                              !   observations stored in file
      integer     batchlendim ! length of each batch of
                              !   observations stored in file
      integer     idlendim    ! string length of each WNO station id
      integer     ntypesdim   ! number of WMO rawinsondes types
      integer     nradcordim  ! number of WMO radiation correction codes
      integer     ndaysdim    ! numbers of days
      integer     nsyndim     ! number of synoptic times
      integer     strlendim   ! string length of kt and kx
                              !   names and units

!     scatch space for variable shapes (i.e. dimension )
!     --------------------------------------------------
      integer  dims   ( MAXNCDIM )

!     scratch space for attributes
!     ----------------------------
      integer   intval    (2)
      real      floatval  (2)
      character charstr * ( max_strlen )

!     functions referenced
!     --------------------
      integer     ODS_CalDat  ! Determine calender date in YYYYMMDD
                              !   format
      character * ( max_strlen )
     .            ODS_Case    ! Sets the case for a string
                              !   of characters
      integer     ODS_Handle  ! Obtains the ODS file handle id
      character * ( max_strlen )
     .            ODS_ParmC   ! Sets string  NetCDF file parameters
      integer     ODS_ParmI   ! Sets integer NetCDF file parameters
      real        ODS_ParmR   ! Sets floating point NetCDF file
                              !   parameters
      integer     ODS_StrSize ! Returns the character string length
                              !   excluding trailing blanks
      character * ( max_strlen )
     .            ODS_VerTag  ! Returns the character string length
                              !   excluding trailing blanks

!     Other variables
!     ---------------
      integer  NameSz         ! number of characters in filename
      integer  DimSz          ! temporary storage for NetCDF
                              !   dimension size
      integer  nc_id          ! temporary storage for NetCDF
                              !   file id
      integer  varid          ! temporary storage for NetCDF
                              !   variable id
      integer  charstrsz      ! string length for character
                              !   string attribute
      integer  type_strlen    ! string length for ODS type
                              !   attribute
      integer  version_strlen ! string length for version
                              !   attribute
      character * ( 20 )      ! temporary storage for the input
     .         ods__type      !   parameter, ods_type
      integer  ntypes         ! temporary storage for the
                              !   dimension, ntypes
      integer  nradcor        ! temporary storage for the
                              !   dimension, nradcor
      integer  ndayz          ! temporary storage for the
                              !   dimension, ndays
      integer  n_syn          ! temporary storage for the
                              !   dimension, nsyn
      integer  ierr_temp      ! temporary storage for the
                              !   returned error code
      integer  lstr           ! string length

!.............................................................

!     Set ierr code to valid 
!     ----------------------
      ierr    = NCNoErr

!     Get the number of an used file handle
!     -------------------------------------
      id      = ODS_Handle ( filename, ierr )
      if ( ierr .ne. NCNoErr ) return

!     Set error message flag to ON    ( i.e. prints error mesages )
!     and the fatal error flag to OFF ( i.e. errors are not fatal )
!     -------------------------------------------------------------
      call ncpopt ( NCVERBOS )          ! verbose mode

!     Enter define mode and use NetCDF routines to set the structure
!     and contents of the ODS/HDF file in terms of the CDL
!     --------------------------------------------------------------
      nc_id             = nccre ( filename, NCCLOB, ierr )
      if ( ierr .ne. NCNoErr ) return

!     Store important file information in common
!     ------------------------------------------
      ncid      ( id )  = nc_id
      IOMode    ( id )  = NCWRITE
      filenames ( id )  = filename

!     From this point on, if an error is detected then perform
!     clean up operations before exiting this routine
!     --------------------------------------------------------

!     define NetCDF dimensions
!     ------------------------

!     nbatches     - number of batches
!     --------------------------------
      nbatchesdim  = ncddef  ( nc_id, 'nbatches', NCUNLIM,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     batchlen     - size of each batch
!     ---------------------------------
      DimSz        = ODS_ParmI      ( 'batchlen', Dbatchlen, ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      batchlendim  = ncddef  ( nc_id, 'batchlen', DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     idlen        - string length of each station id
!     -----------------------------------------------
      DimSz        = ODS_ParmI      ( 'Idlen',    Didlen,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      idlendim     = ncddef  ( nc_id, 'Idlen',    DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     ntypes       - number of WMO instrument types
!     ---------------------------------------------
      DimSz        = ODS_ParmI      ( 'ntypes',   Dntypes,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      ntypesdim    = ncddef  ( nc_id, 'ntypes',   DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      ntypes       = DimSz

!     nradcor      - number of WMO radiation correction codes
!     -------------------------------------------------------
      DimSz        = ODS_ParmI      ( 'nradcor',  Dnradcor,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      nradcordim   = ncddef  ( nc_id, 'nradcor',  DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      nradcor      = DimSz

!     ndays        - maximum number of days stored in file
!     ----------------------------------------------------
      DimSz        = ODS_ParmI      ( 'ndays',    Dndays,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      ndaysdim     = ncddef  ( nc_id, 'ndays',    DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      ndayz        = DimSz ! store this value for later access

!     nsyn - number of synoptic times per day
!     ---------------------------------------
      DimSz        = ODS_ParmI      ( 'nsyn',     Dnsyn,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      nsyndim      = ncddef  ( nc_id, 'nsyn',     DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      n_syn        = DimSz ! store this value for later access

!     strlen - the size of the string in each list entry
!     --------------------------------------------------
      DimSz        = ODS_ParmI     ( 'strlen',   Dstrlen,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      strlendim    = ncddef ( nc_id, 'strlen',   DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     Define NetDCF variables and their attributes for text
!     describing data types and souces currently available.
!     -----------------------------------------------------
!
!     note: Comments that begin with "*** attribute" refer
!           to NetCDF attributes
!     ----------------------------------------------------

!     it_names - name for each rawinsonde index number
!     ------------------------------------------------
      dims  ( 1 ) = strlendim
      dims  ( 2 ) = ntypesdim
      varid       = ncvdef (  nc_id, 'it_names',  NCCHAR,   2,
     .                        dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,    varid, 'name',      NCCHAR,  28,
     .             'Name of WMO instrument types',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     rc_descr - descr for each rawinsonde index number
!     ------------------------------------------------
      dims  ( 1 ) = strlendim
      dims  ( 2 ) = nradcordim
      varid       = ncvdef (  nc_id, 'rc_descr',  NCCHAR,   2,
     .                        dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,    varid, 'name',      NCCHAR,  29,
     .             'WMO solar/IR correction flags',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     end of text describing OMS attributes
!     -------------------------------------

!     For every NetCDF variable defined below with the
!     NetCDF dimension nbatches, the data will be read
!     and/or written in batches. That is each batch
!     contains a subset of the values to be processed.
!     The maximum size of each subset is dbatchlen as
!     defined in the header file, ods_hdf.h or possibly a
!     user specified number.  (through the use of the
!     routines, ODS_SetParmI and ODS_ParmI ) The values
!     are being processed in this manner in order to
!     prevent the execution of the NetCDF library routines
!     from being prohibitively slow, especially when
!     read/writing in HDF format.  However, the number of
!     values that can be stored in the file remains
!     unlimited since the NetCDF dimension nbatches will
!     be set to unlimited.
!     ---------------------------------------------------

!     Define NetDCF variables and their attributes for
!     space-time coordinates including ...
!     ------------------------------------------------
!
!     note: Comments that begin with "*** attribute" refer
!           to NetCDF attributes
!     ----------------------------------------------------


!     -------------------
!     xm - metadata index
!     -------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'xm',        NCFLOAT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,   14,
     .             'Metadata Index',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - missing value set to DAO standard
!     -------------------------------------------------
      floatval      ( 1 ) =  missing_ob
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'missing_value',NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801



!     -------------------
!     Id - WMO station ID
!     -------------------
      dims        ( 1 ) = idlendim
      dims        ( 2 ) = batchlendim
      dims        ( 3 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'Id',        NCChar,   3,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCChar,  14,
     .             'WMO station ID',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801



!     --------
!     latitude
!     --------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'lat',       NCSHORT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,   8,
     .             'Latitude',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - units are degrees north
!     ---------------------------------------
      call ncaptc ( nc_id,     varid,  'units',        NCCHAR,  13,
     .             'degrees north',
     .              ierr)
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - valid range is from -90 to 90
!     ---------------------------------------------
      floatval    ( 1 ) = -90
      floatval    ( 2 ) =  90
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'valid_range',  NCFLOAT,  2,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - scale factor selected to maximize resolution
!                     when stored as a two byte integer
!     ------------------------------------------------------------
c      floatval    ( 1 ) =  90.0 / ( 2.0 ** 15 - 1.0 )
      floatval    ( 1 ) =  1.0 / 100.0
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'scale_factor', NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801



!     ---------
!     longitude
!     ---------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'lon',       NCSHORT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,   9,
     .             'Longitude',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - units are degree east
!     -------------------------------------
      call ncaptc ( nc_id,     varid,  'units',        NCCHAR,  12,
     .             'degrees east',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - value at 90 degrees east is 90.0
!     ------------------------------------------------
      floatval    ( 1 ) =   90
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'value_at_90E', NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - value at 90 degrees west is -90.0
!     -------------------------------------------------
      floatval    ( 1 ) =  -90
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'value_at_90W', NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - valid range is from -180 to 180 degrees
!     -------------------------------------------------------
      floatval    ( 1 ) = -180
      floatval    ( 2 ) =  180
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'valid_range',  NCFLOAT,  2,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - scale factor selected to maximize resolution
!                     when stored as a two byte integer
!     ------------------------------------------------------------
      floatval    ( 1 ) =  1.0 / 100.0d0
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'scale_factor', NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801



!     ------------------------
!     elev - Station elevation
!     ------------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'elev',      NCSHORT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  17,
     .             'Station elevation',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - units
!     ---------------------
      call ncaptc ( nc_id,     varid,  'units',        NCCHAR,  27,
     .             'meters above mean sea level',
     .              ierr)
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - valid minimum
!     -----------------------------
      intval      ( 1 ) =  -2000
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'valid_min',    NCLONG,   1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - missing value
!     -----------------------------
      intval      ( 1 ) =  -10000
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'missing_value',NCLONG,   1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801


!     -----------------------
!     itype - Instrument type
!     -----------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'itype',     NCBYTE,   2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  25,
     .             'WMO instrument type index',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

      intval      ( 1 ) =  0
      intval      ( 2 ) =  ntypes - 1
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'valid_range',  NCLONG,   2,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801


!     ---------------------------------------
!     radcor - solar/infrared correction flag
!     ---------------------------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'radcor',    NCBYTE,   2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  28,
     .             'WMO solar/IR correction flag',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

      intval      ( 1 ) =  0
      intval      ( 2 ) =  nradcor - 1
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'valid_range',  NCLONG,   2,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801


!     -------------------
!     ltime - launch time
!     -------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'ltime',     NCSHORT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  11,
     .             'Launch time',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - units
!     ---------------------
      call ncaptc ( nc_id,     varid,  'units',        NCCHAR,  27,
     .             'minutes since synoptic time',
     .              ierr)
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - missing value
!     -----------------------------
      intval      ( 1 ) =  -10000
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'missing_value',NCLONG,   1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801


!     Define NetDCF variables and their attributes for
!     Tables of pointers including ...
!     ------------------------------------------------
!
!     days
!     ----
      dims         ( 1 ) = ndaysdim
      varid              = ncvdef ( nc_id, 'days',     NCLONG,   1,
     .                              dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     syn_beg
!     -------
      dims         ( 1 ) = nsyndim
      dims         ( 2 ) = ndaysdim
      varid              = ncvdef ( nc_id, 'syn_beg',  NCLONG,   2,
     .                              dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  39,
     .             'Beginning of synoptic hour for each day',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - reference date which must be identical 
!                     to the same attribute for the variable,
!                     julian.  Thus, the date can be identified
!                     by using either one of two parameter names.
!     -----------------------------------------------------------
      charstr            =  ODS_ParmC ( 'julian:reference_date',
     .                                   Dref_date, ierr )
*     *** attribute - the julian day at the reference date 
*                     which must be identical to the same
*                     attribute for the variable, julian.
*                     Thus, this attribute can be identified
*                     by using either one of two parameter
*                     names.
*     ---------------------------------------------------
      charstr            =  ODS_ParmC ( 'syn_beg:reference_date',
     .                                   charstr,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      charstrsz          =  index ( charstr, ' ' ) - 1
      if ( charstrsz .le. 0 ) charstrsz = len ( charstr )
      call ncaptc ( nc_id,     varid,           'reference_date', 
     .              NCCHAR,
     .              charstrsz,
     .              charstr ( 1 : charstrsz ),
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801


!     *** attribute - Julian day at reference date
!     --------------------------------------------
      intval       ( 1 ) =  ODS_ParmI
     .                       ( 'syn_beg:value_at_reference_date',
     .                          Dref_jday,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      call ODS_NCAPTI 
     .            ( nc_id,     varid,  'value_at_reference_date', 
     .              NCLONG,    1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - julian day corresponding to the first
!                     block of data in the file
!     -----------------------------------------------------
      intval       ( 1 ) = first_jday
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'first_julian_day',
     .              NCLONG,    1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - julian day corresponding to the last
!                     block of data written to file
!     ----------------------------------------------------
      intval      ( 1 )  =  first_jday - 1
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'latest_julian_day',
     .              NCLONG,    1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - synoptic hour corresponding to the
!                     last block of data written to file
!     --------------------------------------------------
      intval      ( 1 )  =  0
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'latest_synoptic_hour',  
     .               NCLONG,   1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     syn_len
!     -------
      dims       (  1 )  = nsyndim
      dims       (  2 )  = ndaysdim
      varid              = ncvdef ( nc_id, 'syn_len',   NCLONG,  2,
     .                              dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     *** attribute - name
!     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  36,
     .             'Number of observations for syn. time',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     end of pointers
!     ---------------

!     end of definitions of NetCDF variables and their attributes
!     -----------------------------------------------------------

!     define NetCDF global attributes for ...
!     ---------------------------------------

!     source of this ODS data file
!     ----------------------------
      call ncaptc ( nc_id,  NCGLOBAL,  'source',       NCCHAR,  62,
     .             'Global Modeling and Assimilation Office, Code 610.1, NASA/GSFC',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     title of the data file
!     ----------------------
      call ncaptc ( nc_id,  NCGLOBAL,  'title',        NCCHAR,  45,
     .             'RAOB Observational Metadata Stream (OMS) File',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801


!     ODS version number
!     ------------------
      charstr        = ODS_VerTag ( 0, ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      version_strlen = ODS_StrSize ( charstr )
      call ncaptc ( nc_id,  NCGLOBAL,  'version',      NCCHAR,
     .              version_strlen,
     .              charstr ( 1:version_strlen ),
     .              ierr ) 
      if ( ierr .ne. NCNoErr ) go to 801

!     email address
!     -------------
      call ncaptc ( nc_id,  NCGLOBAL,  'data_info',    NCCHAR,  31,
     .             'Contact data@gmao.gsfc.nasa.gov',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     history tag
!     -----------
      history               = ' '
      call ncaptc ( nc_id,  NCGLOBAL,  'history',      NCCHAR,
     .              max_strlen,
     .              history,
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

!     end of global attributes
!     ------------------------

!     leave define mode
!     -----------------
      call ncendf ( nc_id, ierr)
      if ( ierr .ne. NCNoErr ) go to 801

!     Intialize pointers ( stored in common )
!     ---------------------------------------
      call ODS_ReSetP ( id, ierr )
      if ( ierr .ne. NCNoErr ) go to 802

!     Leave file open for writing
!     ---------------------------

      return
!     ------

!     Clean up by deleting the file and
!     setting io mode to CLOSED
!     ---------------------------------
 801  continue
      call ncabor ( nc_id, ierr_temp )
      IOMode ( id ) = CLOSED

      return
!     ------

!     Clean up by closing the file and
!     setting io mode to CLOSED
!     --------------------------------
 802  continue
      call ncclos ( nc_id, ierr_temp )
      IOMode ( id ) = CLOSED

      return
!     ------

      end
