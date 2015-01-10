* file: ods_hdf.h -
*
*     Include file for ODS 
*
*................................................................

*     array and character string dimensions
*     -------------------------------------
      integer     id_max           ! maximum number of available 
                                   !   ODS/HDF file handles
      parameter ( id_max      = 10   )
      integer     mdays            ! maximum number of days
      parameter ( mdays       = 255  )
      integer     msyn             ! maximum synoptic times per day
      parameter ( msyn        = 24   )
      integer     MPtr             ! maximum number of pointers
                                   !   allowed in temporary storage
      parameter ( MPtr        = mdays * msyn )
      integer     mstrings         ! maximum number of strings
      parameter ( mstrings    = 200  )
      integer     max_strlen       ! maximum string length
      parameter ( max_strlen  = 180  )
      integer     MChar            ! maximum number of characters 
                                   !   allowed in temporary storage
      parameter ( MChar       = max_strlen * mstrings  )

*     symbolic constants for ...
*     --------------------------
      character   Blank * ( 1 )     ! ' ', i.e. a blank
      parameter ( Blank       = ' ' )
      integer     hour_max          ! number of hours in a day
      parameter ( hour_max    = 24  ) 
      integer     ON                ! switch is on
      parameter ( ON          = 1   )
      integer     OFF               ! switch is off
      parameter ( OFF         = 0   )
      integer     Not_Found         ! seach status of not found
      parameter ( Not_Found   = 0   )

*     parameters for ODS/HDF files
*     ----------------------------
      integer     ncid      ( id_max ) ! NetCDF file id
      character * ( max_strlen ) 
     .            filenames ( id_max ) ! list of HDF/ODS file names
      integer     IOMode    ( id_max ) ! integers specifying input
                                       !   / output status status
      integer     CLOSED               ! file is closed ( initial
      parameter ( CLOSED      =  -1 )  !   or default value )

*     error codes supplementing those in the header file, netcdf.inc
*     --------------------------------------------------------------
      integer     ODS_DimErr        ! inappropriate dimension size
      parameter ( ODS_DimErr  =  -3 )
      integer     ODS_BadInt        ! inappropriate FORTRAN interface
      parameter ( ODS_BadInt =   -4 )
      integer     ODS_BadLen        ! inappropriate length of each input
      parameter ( ODS_BadLen =   -5 )  ! string

*     constants used to initialize or reset variables
*     -----------------------------------------------
      integer     CLEAR             ! initialize to zero
      parameter ( CLEAR       =   0 ) 

*     Storage for setting NetCDF file parameters
*     ------------------------------------------
      integer     Parm_mlist          ! maximum size of the
                                      !   parameter list
      parameter ( Parm_mlist = 200 )
      integer     ParmI_nlist         ! size of the list of integer
                                      !   parameters
      integer     ParmR_nlist         ! size of the list of floating
                                      !   point parameters
      integer     ParmC_nlist         ! size of the list of characters
                                      !   string parameters
      character   ParmI_list   ( Parm_mlist ) * ( max_strlen )
                                      ! names of each integer parameter
      character   ParmR_list   ( Parm_mlist ) * ( max_strlen )
                                      ! names of each floating point
                                      !   parameter
      character   ParmC_list   ( Parm_mlist ) * ( max_strlen )
                                      ! names of each character string
                                      !   parameter
      integer     ParmI_val    ( Parm_mlist )
                                      ! specified values for each
                                      !   integer parameter
      real        ParmR_val    ( Parm_mlist )
                                      ! specified values for each
                                      !   floating point parameter
      character   ParmC_string ( Parm_mlist ) * ( max_strlen )
                                      ! specified values for each
                                      !   character string parameter

*     defaults for NetCDF file dimension sizes
*     ----------------------------------------
      integer     Dbatchlen           ! size of each batch
      parameter ( Dbatchlen = 1000 )
      integer     Didlen              ! string length of WMO station id
      parameter ( Didlen    = 8  )
      integer     Dntypes             ! number of WMO instrument types
      parameter ( Dntypes   = 256 )
      integer     Dnradcor            ! number of WMO radiation corr codes
      parameter ( Dnradcor  = 16 )
      integer     Dnkt                ! maximum data type index
      parameter ( Dnkt      = 7  )
      integer     Dnkx                ! maximum data source index
      parameter ( Dnkx      = 87 )  
      integer     Dnqcx               ! maximum number of possible
                                      !   values for the quality
                                      !   control exclusion mark
      parameter ( Dnqcx     = 100 )  
      integer     Dndays              ! maximum number of days
      parameter ( Dndays    = 255 )
      integer     Dnsyn               ! maximum number of 
      parameter ( Dnsyn     = 4  )    !   synoptic times
      integer     Dstrlen             ! size of each string in
      parameter ( Dstrlen   = max_strlen ) ! a list entry

*     defaults for NetCDF attributes
*     ------------------------------
      character   Dref_date *   ( max_strlen )   ! reference date
      parameter ( Dref_date =    '1968-05-23' )
      integer     Dref_jday                      ! reference julian day
      parameter ( Dref_jday =     2440000 )
      real        missing_ob                     ! missing observation
      parameter ( missing_ob  =   1.0e15 )       !   DAO standard

*     history and version attributes ( The list of
*     versions is set in block data ODS_VerList )
*     --------------------------------------------
      integer     NVer
      parameter ( NVer = 16 )
      integer     NVersions
      character   history                 * ( max_strlen )
      character   version_list   ( NVer ) * ( max_strlen )

*     Tables of pointers for blocks of data in a file partitioned
*     according to julian day and synoptic hour.
*     -----------------------------------------------------------
      integer     days          ( mdays, id_max ) ! 
      integer     syn_beg ( msyn, mdays, id_max ) ! beginning of block
      integer     syn_len ( msyn, mdays, id_max ) ! block size ( i.e. the 
                                                  !   number of values )

*     Julian day offset ( used in determining table indices )
*     -------------------------------------------------------
      integer     julian_offset ( id_max )  !

*     Integers defining the dimension of the tables
*     ---------------------------------------------
      integer     nsyn          ( id_max )  ! number of synoptic
                                            !   times per day
      integer     ndays         ( id_max )  ! maximum number of
                                            !   days on file

*     data corresponding to the last block of data written fo file
*     ------------------------------------------------------------
      integer     latest_day    ( id_max )  ! julian day
      integer     latest_hour   ( id_max )  ! synoptic hour

*     pointer data for writing in append mode
*     ---------------------------------------
      integer     append_mode   ( id_max )  ! = ON if in append mode
                                            ! = OFF otherwise
      integer     append_beg    ( id_max )  ! beginning of block
      integer     append_len    ( id_max )  ! block size

*     scratch (or work) space for ...
*     -------------------------------
      integer     ptrtemp       ( MPtr  )   ! pointers
      character   strtemp *     ( MChar )   ! array of character strings

*     common storage for ...
*     ----------------------

*     ODS/HDF file parameters
*     -----------------------
      common / ods_io   / ncid,        IOMode,      filenames

*     ODS/HDF version and history attributes
*     --------------------------------------
      common / ods_verh / NVersions, 
     .                    version_list,
     .                    history

*     NetCDF file parameters
*     ----------------------
      common / ods_parm / ParmI_nlist, ParmR_nlist, ParmC_nlist,
     .                    ParmI_val,   ParmR_val,   ParmC_string,
     .                    ParmI_list,  ParmR_list,  ParmC_list

*     pointers and current values
*     ---------------------------
      common / ods_ptr  / days,
     .                    syn_beg,
     .                    syn_len,
     .                    ndays,
     .                    nsyn,
     .                    julian_offset,
     .                    latest_day,
     .                    latest_hour,
     .                    append_mode,
     .                    append_beg,
     .                    append_len

*     scratch space
*     -------------
      common / ods_scr  / ptrtemp,    strtemp


