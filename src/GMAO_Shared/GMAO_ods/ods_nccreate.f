
*....................................................................


      subroutine ODS_NCCreate ( id, filename, ods_type, 
     .                                        first_jday, ierr )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 610.1, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:          ODS\_NCCreate
! 
! !DESCRIPTION:      
!     Creates an ODS file and initializes the pointers
!     specifying the location and size of data blocks that
!     are partitioned according to julian day and synoptic
!     hour.  Most of this routine is dedicated to calling
!     NetCDF interface routines to determine the structure
!     of the ODS file.  This structure is described in
!     terms of the network common data form language (CDL)
!     and in the header as produced with the NetCDF interface
!     to HDF utility ncdump (see Appendix C of da Silva and
!     Redder, 1996 ).  For further information regarding
!     NetCDF routines and CDL see Rew et al. (1993).
!
! !REFERENCE: 
!
!     Rew, Russ, Glenn Davis and Steve Emmerson, 1993:
!        NetCDF User's Guide, Unidata Program Center,
!        University Corporation for Atmospheric Research,
!        National Science Foundation, Boulder, CO.
!
!
! !INTERFACE:
!     call ODS_NCCreate  ( id, filename, ods_type, 
!                                        first_jday, ierr )
!
! !INPUT PARAMETERS:
      character * (*)  filename             ! name of ODS file
      character * (*)  ods_type             ! = 'pre_anal', if only
                                            !    repacked obs are to be
                                            !    stored.
                                            ! = 'post_anal', if
                                            !    data from assimilation
                                            !    is to be stored
      integer          first_jday           ! first julian day
!
! !OUTPUT PARAMETERS: 
      integer          id                   ! ODS file handle
      integer          ierr                 ! error code
!
! !SEE ALSO:
!     ODS_NCOpen  ( opens  the ODS file )
!     ODS_Close   ( closes the ODS file )!
!     ODS_Julian()  convert the "calendar" date to the Julian day

! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!
! !LIBRARIES ACCESSED:
!     NetCDF
!
! !REVISION HISTORY:
!     3 Sept 1996  Redder   Origional version
!    19 Aug  1997  Redder   Removed the ODS variable for the
!                           modification flag
!    17 Sept 1997  Redder   Modified code in order store the
!                           QC flag as a 4-byte integer
!    13 Apr  1998  Redder   Changed variable name "level" to
!                          "lev"; removed the variables, qc_flag,
!                           julian and km, added the variables,
!                           qcexcl, qchist, xm; modified the
!                           variable, time; hanged variable names,
!                           kt_max and kx_max, to nkt and nkx;
!                           modified code to handle the version
!                           tag with greater flexibility.  These
!                           changes were made to create ODS
!                           version 2.00
!     22 Apr 1998  Redder   Changed the offset and valid range of
!                           the NetCDF variable, time, to allow for  
!                           for negative values.  Time values can
!                           be negative when, for example, the
!                           sample time is 21Z of the synptic time
!                           0Z of the first Julian day on file.
!     16 Jun 1998  Redder   Modified scale factor for the variables
!                           lat and lon
!     07 Jul 1998  Redder   Made minor corrections to the strings
!                           for the NetCDF global attributes,
!                           data_info and title.
!     11 Mar 1999  Redder   Made changes to store ODS observation
!                           attribute, ks (data sounding index), as
!                           a 4-byte instead of 2-byte integer. 
!                           The maximum acceptable value was
!                           increased from ~64,000 to ~1 billion.
!     15 Mar 1999  Redder   Added global attribute, type.
!     16 Feb 2000  Todling  Rename stdio.h to ods_stdio.h
!     16 Feb 2000  RT/CR    Added xvec observation attribute
!     04 Mar 2002  Redder   Changed the lower limit of the time attribute
!                           attribute by -60 min (to allow the raob launch
!                           time to be stored)
!     16 Sep 2002  Redder   Minor correction in name of attribute syn_beg
!     20 Jul 2005  RT/CR    Changed TimeMin to accommodate 10-day time window
!
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     dimension ids for ...
*     ---------------------
      integer     nbatchesdim ! number of batches of
                              !   observations stored in file
      integer     batchlendim ! length of each batch of
                              !   observations stored in file
      integer     nktdim      ! list of kt names and units
      integer     nkxdim      ! list of kx names
      integer     nqcxdim     ! meaning of each possible value of the
                              !   quality control exclusion mark
      integer     ndaysdim    ! numbers of days
      integer     nsyndim     ! number of synoptic times
      integer     strlendim   ! string length of kt and kx
                              !   names and units

*     scatch space for variable shapes (i.e. dimension )
*     --------------------------------------------------
      integer  dims   ( MAXNCDIM )

*     scratch space for attributes
*     ----------------------------
      integer   intval    (2)
      real      floatval  (2)
      character charstr * ( max_strlen )

*     functions referenced
*     --------------------
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

*     Other variables
*     ---------------
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
      integer  nkt            ! temporary storage for the
                              !   dimension, nkt
      integer  nkx            ! temporary storage for the
                              !   dimension, nkx
      integer  nqcx           ! temporary storage for the
                              !   dimension, nqcx
      integer  ndayz          ! temporary storage for the
                              !   dimension, ndays
      integer  n_syn          ! temporary storage for the
                              !   dimension, nsyn
      integer  ierr_temp      ! temporary storage for the
                              !   returned error code
      integer  lstr           ! string length
      integer  TimeMin        ! minimum sample time

*.............................................................

*     Set ierr code to valid 
*     ----------------------
      ierr    = NCNoErr

*     Get the number of an used file handle
*     -------------------------------------
      id      = ODS_Handle ( filename, ierr )
      if ( ierr .ne. NCNoErr ) return

*     ensure that any character element in the string
*     parameter, ods_type, are in lower case
*     -----------------------------------------------
      ods__type = ODS_Case ( ods_type, 'lower' )

*     Give full name to ods file type.
*     If ods type is not valid, return
*     --------------------------------
      if      ( ods__type ( :8 ) .eq. 'pre_anal'  ) then
         ods__type = 'pre_analysis'

      else if ( ods__type ( :9 ) .eq. 'post_anal' ) then
         ods__type = 'post_analysis'

      else
         write ( stderr, * )
         write ( stderr, * )
     .    ' ODS_NCCreate: not a valid ods type: ',
     .                    ods__type
         ierr = NCEInVal
         return

      end if 

*     Set error message flag to ON    ( i.e. prints error mesages )
*     and the fatal error flag to OFF ( i.e. errors are not fatal )
*     -------------------------------------------------------------
      call ncpopt ( NCVERBOS )          ! verbose mode

*     Enter define mode and use NetCDF routines to set the
*     structure and contents of the ODS/HDF file in terms of
*     the CDL
*     ------------------------------------------------------
      nc_id             = nccre ( filename, NCCLOB, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Store important file information in common
*     ------------------------------------------
      ncid      ( id )  = nc_id
      IOMode    ( id )  = NCWRITE
      filenames ( id )  = filename

*     From this point on, if an error is detected then perform
*     clean up operations before exiting this routine
*     --------------------------------------------------------

*     define NetCDF dimensions
*     ------------------------

*     nbatches - number of batches
*     ----------------------------
      nbatchesdim = ncddef  ( nc_id, 'nbatches', NCUNLIM,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     batchlen - size of each batch
*     -----------------------------
      DimSz       = ODS_ParmI      ( 'batchlen', Dbatchlen, ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      batchlendim = ncddef  ( nc_id, 'batchlen', DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     nkt - maximum value for the data type index
*     -------------------------------------------
      DimSz        = ODS_ParmI     ( 'nkt',      Dnkt,      ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      nktdim       = ncddef ( nc_id, 'nkt',      DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      nkt          = DimSz ! store this value for later access

*     nkx - maximum value for the data source index
*     ---------------------------------------------
      DimSz        = ODS_ParmI     ( 'nkx',      Dnkx,      ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      nkxdim       = ncddef ( nc_id, 'nkx',      DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      nkx          = DimSz ! store this value for later access

*     nqcx - maximum number of possible values for the
*            quality control exclusion mark (qcexcl)
*     ------------------------------------------------
      DimSz        = ODS_ParmI     ( 'nqcx',     Dnqcx,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      nqcxdim      = ncddef ( nc_id, 'nqcx',     DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      nqcx         = DimSz ! store this value for later access

*     ndays - maximum number of days stored in file
*     ---------------------------------------------
      DimSz        = ODS_ParmI     ( 'ndays',    Dndays,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      ndaysdim     = ncddef ( nc_id, 'ndays',    DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      ndayz        = DimSz ! store this value for later access

*     nsyn - number of synoptic times per day
*     ---------------------------------------
      DimSz        = ODS_ParmI     ( 'nsyn',     Dnsyn,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      nsyndim      = ncddef ( nc_id, 'nsyn',     DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      n_syn        = DimSz ! store this value for later access

*     strlen - the size of the string in each list entry
*     --------------------------------------------------
      DimSz        = ODS_ParmI     ( 'strlen',   Dstrlen,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      strlendim    = ncddef ( nc_id, 'strlen',   DimSz,     ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     Define NetDCF variables and their attributes for text
*     describing data types and souces currently available.
*     -----------------------------------------------------
*
*     note: Comments that begin with "*** attribute" refer
*           to NetCDF attributes
*     ----------------------------------------------------


*     kt_names for each GEOS/DAS data type index number
*     -------------------------------------------------
      dims  ( 1 ) = strlendim
      dims  ( 2 ) = nktdim
      varid       = ncvdef (  nc_id, 'kt_names',  NCCHAR,  2,
     .                        dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid, 'name',     NCCHAR,  27,
     .             'Name of GEOS/DAS data types',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     kt_units - unit for each GEOS/DAS data type index number
*     --------------------------------------------------------
      dims  ( 1 ) = strlendim
      dims  ( 2 ) = nktdim
      varid       = ncvdef (  nc_id, 'kt_units',  NCCHAR,  2,
     .                        dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,    varid, 'name',      NCCHAR,  33,
     .             'Units for each GEOS/DAS data type',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     kx_names - name for each GEOS/DAS data source index number
*     ----------------------------------------------------------
      dims  ( 1 ) = strlendim
      dims  ( 2 ) = nkxdim
      varid       = ncvdef (  nc_id, 'kx_names',  NCCHAR,  2,
     .                        dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,    varid, 'name',      NCCHAR,  29,
     .             'Name of GEOS/DAS data sources',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     kx_meta - kx specific metadata information
*     ------------------------------------------
      dims  ( 1 ) = strlendim
      dims  ( 2 ) = nkxdim
      varid       = ncvdef (  nc_id, 'kx_meta',   NCCHAR,  2,
     .                        dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,    varid, 'name',      NCCHAR,  32,
     .             'kx specific metadata information',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     qcx_names - meaning of each possible value of the
*                 quality control exclusion mark
*     -------------------------------------------------
      dims  ( 1 ) = strlendim
      dims  ( 2 ) = nqcxdim
      varid       = ncvdef (  nc_id, 'qcx_names', NCCHAR,  2,
     .                        dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      charstr     =  'Meaning of each possible value of '
     .            // 'the quality control exclusion mark'
      charstrsz   = ODS_StrSize ( charstr )
      call ncaptc ( nc_id,    varid, 'name',      NCCHAR,
     .              charstrsz,
     .              charstr ( : charstrsz ),
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     end of text describing data types and souces
*     --------------------------------------------

*     For every NetCDF variable defined below with the
*     NetCDF dimension nbatches, the data will be read
*     and/or written in batches. That is each batch
*     contains a subset of the values to be processed.
*     The maximum size of each subset is dbatchlen as
*     defined in the header file, ods_hdf.h or possibly a
*     user specified number.  (through the use of the
*     routines, ODS_SetParmI and ODS_ParmI ) The values
*     are being processed in this manner in order to
*     prevent the execution of the NetCDF library routines
*     from being prohibitively slow, especially when
*     read/writing in HDF format.  However, the number of
*     values that can be stored in the file remains
*     unlimited since the NetCDF dimension nbatches will
*     be set to unlimited.
*     ---------------------------------------------------

*     Define NetDCF variables and their attributes for
*     space-time coordinates including ...
*     ------------------------------------------------
*
*     note: Comments that begin with "*** attribute" refer
*           to NetCDF attributes
*     ----------------------------------------------------

*     latitude
*     --------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'lat',       NCSHORT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,   8,
     .             'Latitude',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - units are degrees north
*     ---------------------------------------
      call ncaptc ( nc_id,     varid,  'units',        NCCHAR,  13,
     .             'degrees north',
     .              ierr)
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - valid range is from -90 to 90
*     ---------------------------------------------
      floatval    ( 1 ) = -90
      floatval    ( 2 ) =  90
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'valid_range',  NCFLOAT,  2,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - scale factor selected to maximize resolution
*                     when stored as a two byte integer
*     ------------------------------------------------------------
c      floatval    ( 1 ) =  90.0 / ( 2.0 ** 15 - 1.0 )
      floatval    ( 1 ) =  1.0 / 100.0
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'scale_factor', NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801


*     longitude
*     ---------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'lon',       NCSHORT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,   9,
     .             'Longitude',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - units are degree east
*     -------------------------------------
      call ncaptc ( nc_id,     varid,  'units',        NCCHAR,  12,
     .             'degrees east',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801


*     *** attribute - value at 90 degrees east is 90.0
*     ------------------------------------------------
      floatval    ( 1 ) =   90
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'value_at_90E', NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - value at 90 degrees west is -90.0
*     -------------------------------------------------
      floatval    ( 1 ) =  -90
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'value_at_90W', NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - valid range is from -180 to 180 degrees
*     -------------------------------------------------------
      floatval    ( 1 ) = -180
      floatval    ( 2 ) =  180
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'valid_range',  NCFLOAT,  2,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - scale factor selected to maximize resolution
*                     when stored as a two byte integer
*     ------------------------------------------------------------
c      floatval    ( 1 ) =  180.0 / ( 2.0 ** 15 - 1.0 )
      floatval    ( 1 ) =  1.0 / 100.0d0
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'scale_factor', NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     level
*     -----
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'lev',       NCFLOAT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  25,
     .             'Pressure level or channel',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - units are hPa or none
*                    (for satellite channel numbers ) 
*     -----------------------------------------------
      call ncaptc ( nc_id,     varid,  'units',        NCCHAR,  11,
     .             'hPa or none',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801


*     time - minutes since 0:00 GMT of first julian day on file
*     ---------------------------------------------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'time',      NCSHORT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  50,
     .        'minutes since 0:00 GMT of first julian day on file',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - units in minutes since first julian day on
*                     file ( in COARDS; e.g. 1992-10-8 00:00:00.0 )
*     -------------------------------------------------------------
      intval ( 1 )      =  ODS_CalDat ( first_jday )
      intval ( 2 )      =  0
      call ODS_COARDS ( intval ( 1 ), intval ( 2 ), charstr, ierr )
      charstr           = 'minutes since ' // charstr
      charstrsz         =  ODS_StrSize ( charstr )
      call ncaptc ( nc_id,     varid,  'units',        NCCHAR,
     .              charstrsz,
     .              charstr ( : charstrsz ),
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - add offset which is selected to double
*                     the range without exceeding the limits
*                     of a two byte integer
*     --------------------------------------------------------
c      TimeMin           =  -1440 / ( 2 * n_syn ) - 60
!     TimeMin           =  -1440 / n_syn - 60
      TimeMin           =  -14400
      intval      ( 1 ) =  2 ** 15 - 1 + TimeMin
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'add_offset',   NCLONG,   1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - valid range as determined by the limits
*                     of a two byte integer
*     -------------------------------------------------------
      intval      ( 1 ) =  TimeMin
      intval      ( 2 ) =  2 ** 16 - 2 + TimeMin
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'valid_range',  NCLONG,   2,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     end of space-time coordinates
*     -----------------------------

*     Define NetDCF variables and their attributes for
*     data types and sources including ...
*     ------------------------------------------------
*
*     note: Comments that begin with "*** attribute" refer
*           to NetCDF attributes
*     ----------------------------------------------------

*     kt - GEOS/DAS data type
*     -----------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'kt',        NCBYTE,   2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  15,
     .             'Data type index',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - add offset set to 1
*     -----------------------------------
      intval      ( 1 ) =  1
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'add_offset',   NCLONG,   1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

c      intval      ( 1 ) =  2
c      call ODS_NCAPTI
c     .            ( nc_id,     varid,  'scale_factor', NCLONG,  1,
c     .              intval,    ierr )
c      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - valid range that depends on the dimension nkt
*     -------------------------------------------------------------
      intval      ( 1 ) =  1
      intval      ( 2 ) =  nkt
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'valid_range',  NCLONG,   2,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     kx - GEOS/DAS data source
*     -------------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'kx',        NCSHORT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,   17,
     .             'Data source index',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - add offset which is selected to double
*                     the range without exceeding the limits
*                     of a two byte integer
*     --------------------------------------------------------
      intval      ( 1 ) =  2 ** 15
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'add_offset',   NCLONG,   1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - valid range that depends on the dimension nkx
*     -------------------------------------------------------------
      intval      ( 1 ) =  1
      intval      ( 2 ) =  nkx
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'valid_range',  NCLONG,   2,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     ks - GEOS/DAS sounding index
*     ----------------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'ks',        NCLONG,   2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  14,
     .             'Sounding index',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - add offset which is selected to double
*                     the range without exceeding the limits
*                     of a two byte integer
*     --------------------------------------------------------
      intval      ( 1 ) =  0
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'add_offset',   NCLONG,   1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - valid range as determined by the limits
*                     of a two byte integer
*     -------------------------------------------------------
      intval      ( 1 ) =  1
      intval      ( 2 ) =  2 ** 30
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'valid_range',  NCLONG,   2,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     xm - GEOS/DAS metadata index
*     ----------------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'xm',        NCFLOAT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,   8,
     .             'Metadata',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - missing value set to DAO standard
*     -------------------------------------------------
      floatval    ( 1 ) =  missing_ob
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'missing_value',NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801


*     end of data types and sources
*     -----------------------------

*     Define NetDCF variables and their attributes for
*     observations values
*     ------------------------------------------------
*
*     note: Comments that begin with "*** attribute" refer
*           to NetCDF attributes
*     ----------------------------------------------------

*     obs - observation
*     -----------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'obs',       NCFLOAT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  17,
     .             'Observation value',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - missing value set to DAO standard
*                   ( see header file, ods_hdf.h )
*     -------------------------------------------------
      floatval    ( 1 ) =  missing_ob
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'missing_value',NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     Define the variables, omf and oma, only if the 
*     the ods type is set to 'post_nal'
*     ----------------------------------------------
      if ( ods__type ( 1:9 ) .eq. 'post_anal' ) then

*     omf - observation minus forecast
*     --------------------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'omf',       NCFLOAT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  29,
     .             'Observation minus 6h forecast',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - missing value set to DAO standard
*                   ( see header file, ods_hdf.h )
*     -------------------------------------------------
      floatval    ( 1 ) =  missing_ob
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'missing_value',NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     oma - observation minus analysis
*     --------------------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'oma',       NCFLOAT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  26,
     .             'Observation minus analysis',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - missing value set to DAO standard
*                   ( see header file, ods_hdf.h )
*     -------------------------------------------------
      floatval    ( 1 ) =  missing_ob
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'missing_value',NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801


*     xvec - PSAS CG solution vector
*     ------------------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'xvec',      NCFLOAT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  23,
     .             'PSAS CG solution vector',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - missing value set to DAO standard
*                   ( see header file, ods_hdf.h )
*     -------------------------------------------------
      floatval    ( 1 ) =  missing_ob
      call ODS_NCAPTR
     .            ( nc_id,     varid,  'missing_value',NCFLOAT,  1,
     .              floatval,  ierr )
      if ( ierr .ne. NCNoErr ) go to 801


      end if
*     ------

*     Define NetDCF variables and their attributes for
*     quality control flags including ...
*     ------------------------------------------------
*
*     note: Comments that begin with "*** attribute" refer
*           to NetCDF attributes
*     ----------------------------------------------------

*     qcexcl - quality control exclusion mark
*     ---------------------------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'qcexcl',    NCBYTE,   2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  30,
     .             'Quality control exclusion mark',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - valid range as determined by the limits
*                     of a one byte integer
*     -------------------------------------------------------
      intval      ( 1 ) =  0
      intval      ( 2 ) =  2 ** 8 - 1
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'valid_range',  NCLONG,   2,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801


*     qchist - quality control history mark
*     -------------------------------------
      dims        ( 1 ) = batchlendim
      dims        ( 2 ) = nbatchesdim
      varid             = ncvdef ( nc_id, 'qchist',    NCSHORT,  2,
     .                             dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  28,
     .             'Quality control history mark',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - add offset which is selected to double
*                     the range without exceeding the limits
*                     of a one byte integer
*     ------------------------------------------------------
      intval      ( 1 ) =  2 ** 15 - 1
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'add_offset',    NCLONG,  1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - valid range as determined by the limits
*                     of a one byte integer
*     -------------------------------------------------------
      intval      ( 1 ) =  0
      intval      ( 2 ) =  2 ** 16 - 2
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'valid_range',   NCLONG,  2,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     modification flag
*     -----------------
c      dims        ( 1 ) = batchlendim
c      dims        ( 2 ) = nbatchesdim
c      varid             = ncvdef ( nc_id, 'mod_flag',  NCBYTE,  2,
c     .                             dims,   ierr )
c      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
c      call ncaptc ( nc_id,     varid,  'name',         NCCHAR, 17,
c     .             'Modification flag',
c     .              ierr )
c      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - add offset which is selected to double
*                     the range without exceeding the limits
*                     of a two byte integer
*    -------------------------------------------------------
c      intval      ( 1 ) =  0
c      call ODS_NCAPTI
c     .            ( nc_id,     varid,  'add_offset',   NCLONG,  1,
c     .              intval,    ierr )
c      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - valid range as determined by the limits
*                     of a one byte integer
*     -------------------------------------------------------
c      intval      ( 1 ) =  0
c      intval      ( 2 ) =  255
c      call ODS_NCAPTI
c     .            ( nc_id,     varid,  'valid_range',  NCLONG,  2,
c     .              intval,    ierr )
c      if ( ierr .ne. NCNoErr ) go to 801


*     end of quality control flags
*     ----------------------------

*     Define NetDCF variables and their attributes for
*     Tables of pointers including ...
*     ------------------------------------------------
*
*     note: Comments that begin with "*** attribute" refer
*           to NetCDF attributes
*     ----------------------------------------------------

*     days
*     ----
      dims         ( 1 ) = ndaysdim
      varid              = ncvdef ( nc_id, 'days',     NCLONG,   1,
     .                              dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     syn_beg
*     -------
      dims         ( 1 ) = nsyndim
      dims         ( 2 ) = ndaysdim
      varid              = ncvdef ( nc_id, 'syn_beg',  NCLONG,   2,
     .                              dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  39,
     .             'Beginning of synoptic hour for each day',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - reference date which must be identical 
*                     to the same attribute for the variable,
*                     julian.  Thus, the date can be identified
*                     by using either one of two parameter names.
*     -----------------------------------------------------------
      charstr            =  ODS_ParmC ( 'julian:reference_date',
     .                                   Dref_date, ierr )
      if ( ierr .ne. NCNoErr ) go to 801
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

*     *** attribute - the julian day at the reference date 
*                     which must be identical to the same
*                     attribute for the variable, julian.
*                     Thus, this attribute can be identified
*                     by using either one of two parameter
*                     names.
*     ---------------------------------------------------
      intval       ( 1 ) =  ODS_ParmI
     .                       (  'julian:value_at_reference_date',
     .                          Dref_jday, ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      intval       ( 1 ) =  ODS_ParmI
     .                       ( 'syn_beg:value_at_reference_date',
     .                          intval,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      call ODS_NCAPTI 
     .            ( nc_id,     varid,  'value_at_reference_date', 
     .              NCLONG,    1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - julian day corresponding to the first
*                     block of data in the file
*     -----------------------------------------------------
      intval       ( 1 ) = first_jday
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'first_julian_day',
     .              NCLONG,    1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - julian day corresponding to the last
*                     block of data written to file
*     ----------------------------------------------------
      intval      ( 1 )  =  first_jday - 1
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'latest_julian_day',
     .              NCLONG,    1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - synoptic hour corresponding to the
*                     last block of data written to file
*     --------------------------------------------------
      intval      ( 1 )  =  0
      call ODS_NCAPTI
     .            ( nc_id,     varid,  'latest_synoptic_hour',  
     .               NCLONG,   1,
     .              intval,    ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     syn_len
*     -------
      dims       (  1 )  = nsyndim
      dims       (  2 )  = ndaysdim
      varid              = ncvdef ( nc_id, 'syn_len',   NCLONG,  2,
     .                              dims,   ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     *** attribute - name
*     --------------------
      call ncaptc ( nc_id,     varid,  'name',         NCCHAR,  36,
     .             'Number of observations for syn. time',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     end of pointers
*     ---------------

*     end of definitions of NetCDF variables and their attributes
*     -----------------------------------------------------------

*     define NetCDF global attributes for ...
*     ---------------------------------------

*     source of this ODS data file
*     ----------------------------
      call ncaptc ( nc_id,  NCGLOBAL,  'source',       NCCHAR,  62,
     .             'Global Modeling and Assimilation Office, Code 610.1, NASA/GSFC',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     title of the data file
*     ----------------------
      call ncaptc ( nc_id,  NCGLOBAL,  'title',        NCCHAR,  45,
     .             'GEOS DAS Observational Data Stream (ODS) File',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     ODS type
*     --------
      type_strlen = ODS_StrSize ( ods__type )
      call ncaptc ( nc_id,  NCGLOBAL,  'type',         NCCHAR,
     .              type_strlen,
     .              ods__type ( 1:type_strlen ),
     .              ierr ) 
      if ( ierr .ne. NCNoErr ) go to 801

*     ODS version number
*     ------------------
      charstr        = ODS_VerTag ( 0, ierr )
      if ( ierr .ne. NCNoErr ) go to 801
      version_strlen = ODS_StrSize ( charstr )
      call ncaptc ( nc_id,  NCGLOBAL,  'version',      NCCHAR,
     .              version_strlen,
     .              charstr ( 1:version_strlen ),
     .              ierr ) 
      if ( ierr .ne. NCNoErr ) go to 801

*     email address
*     -------------
      call ncaptc ( nc_id,  NCGLOBAL,  'data_info',    NCCHAR,  31,
     .             'Contact data@gmao.gsfc.nasa.gov',
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     history tag
*     -----------
      history               = ' '
      call ncaptc ( nc_id,  NCGLOBAL,  'history',      NCCHAR,
     .              max_strlen,
     .              history,
     .              ierr )
      if ( ierr .ne. NCNoErr ) go to 801

*     end of global attributes
*     ------------------------

*     leave define mode
*     -----------------
      call ncendf ( nc_id, ierr)
      if ( ierr .ne. NCNoErr ) go to 801

*     Intialize pointers ( stored in common )
*     ---------------------------------------
      call ODS_ReSetP ( id, ierr )
      if ( ierr .ne. NCNoErr ) go to 802

*     Leave file open for writing
*     ---------------------------

      return
*     ------

*     Clean up by deleting the file and
*     setting io mode to CLOSED
*     ---------------------------------
 801  continue
      call ncabor ( nc_id, ierr_temp )
      IOMode ( id ) = CLOSED

      return
*     ------

*     Clean up by closing the file and
*     setting io mode to CLOSED
*     --------------------------------
 802  continue
      call ncclos ( nc_id, ierr_temp )
      IOMode ( id ) = CLOSED

      return
*     ------

      end
