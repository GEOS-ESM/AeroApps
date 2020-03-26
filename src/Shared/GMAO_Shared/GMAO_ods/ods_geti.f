
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_GetI() --- Reads data from file to native integers
! 
! !DESCRIPTION:
!    \label{ODS:GetI}
!     Reads the contents of an integer variable for a full synoptic
!     time from an opened ODS file.
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
! !INTERFACE: 
      subroutine ODS_GetI ( id, VarName, julian_day, syn_hour,
     .                          nval,    values,     ierr )

!
! !INPUT PARAMETERS:
      implicit        NONE
      integer         id              ! ODS file handle as returned
                                      !  from ODS_Create() or
                                      !  ODS_Open().
      character * (*) VarName         ! name of variable
      integer         julian_day      ! Julian Day.  Use the function
                                      !  ODS_Julian to obtain this
                                      !  number.
      integer         syn_hour        ! hour of synoptic time since
                                      !  0:00 GMT (e.g., 0, 6, 12, 18)

! !INPUT/OUTPUT PARAMETERS:
      integer         nval            ! on input: maximum number of
                                      !  values (usually determined
                                      !  by the space allocated for
                                      !  storage)
                                      ! on output: the number of
                                      !  values obtained from the
                                      !  file 
                                      ! note: If the input value is
                                      !  smaller than the number
                                      !  of values to be read, then
                                      !  the routine exits with a
                                      !  system error message and
                                      !  code of ODS_DimErr ( = -3 ).
                                      !  Also, nval is set to the
                                      !  minimum value required in
                                      !  order for the routine to
                                      !  execute successfully.

! !OUTPUT PARAMETERS:
      integer         values ( nval ) ! Values read from file
      integer         ierr            ! Error code. If non-zero, an 
                                      !  error has occurred. For a list
                                      !  of possible values, see the
                                      !  description section of this
                                      !  prologue.
!
! !SEE ALSO:
!     ODS_PutC()    write a list of character  strings to   an ODS file
!     ODS_PutI()    write a list of integer    values  to   an ODS file
!     ODS_PutR()    write a list of real       values  to   an ODS file
!     ODS_GetC()    get   a list of characters strings from an ODS file
!     ODS_GetR()    get   a list of real       values  from an ODS file
!     ODS_NGet()    get the number of values for a given julian day
!                    and synoptic hour
!     ODS_Julian()  convert the "calendar" date to the Julian day
!     ods_hdf.h     include file defining internal constants
!                    and global variables, and setting up data
!                    structures
!
! !REVISION HISTORY: 
!     02Oct1995   da Silva    Specification.
!     16May1996   C. Redder   Original code
!     13Apr1998   C. Redder   Updated documentation to reflect changes
!                             made in other routines to create ODS
!                             version 2.00
!     19Nov1999   C. Redder   Added a latex label in and moved the
!                             subroutine statement into the prologue.
!                             Modified the comments for the return status
!                             code.
!     16Feb2000   Todling     Rename stdio.h to ods_stdio.h
!     30Nov2001   C. Redder   Fixed bug that prevents last batch from 
!                             being processed twice when nval = 0
!     16Sep2002   C. Redder   Updated prologue and added error handling
!                             for call to NetCDF routine, NCDINQ
!EOP
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Functions referenced
*     --------------------
      character * ( 10 )
     .         ODS_VarType ! recommended FORTRAN variable type for data

*     Arguments for NetCDF library routines
*     -------------------------------------
      integer  nc_id       ! file id
      integer  varid       ! variable id
      integer  start  ( 2 )! starting file index for hyperslab
      integer  count  ( 2 )! length of edges in hyperslab
      character * ( MAXNCNAM )
     .         Var_Name    ! variable name ( returned from NetCDF
                           !   routine NCVINQ )
      integer  Var_Type    ! variable type
      integer  NVDims      ! number of dimensions
      integer  DimIDs ( 2 )! dimension ids
      integer  DimID       ! dimension id
      integer  NVAtts      ! number of attributes
      character * ( MAXNCNAM )
     .         DimName     ! variable name ( returned from NetCDF
                           !   routine NCDINQ )

*     Indices defining the set of observation data 
*     corresponding to a synoptic time in the ODS file
*     ------------------------------------------------
      integer  synbeg      ! starting index for the data block
      integer  synlen      ! number of observation reports in
                           !   the data block

*     Parameters for defining the observation
*     data batches to be written to the ODS file 
*     ------------------------------------------
      integer  batch_len   ! number of values in each batch
      integer  first_batch ! index for the first batch
      integer  last_batch  ! index for the last batch
      integer  first_ival  ! data value index of the first
                           !   value in the first batch
      integer  first_nval  ! number of values in the first batch
      integer  last_nval   ! number of values in the last batch
      integer  ibatch      ! batch number index for do loop

*     Other variables
*     ---------------
      integer  ival        ! index value for the input array values
      integer  nval_max    ! maximum number of values
      character * ( 10 ) 
     .         VarType     ! recommended FORTRAN variable type for
                           !   data

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

*     Check to determine if variable is appropriate for this routine
*     --------------------------------------------------------------
      VarType    = ODS_VarType ( VarName )
      if ( VarType ( :7  ) .ne. 'integer'     .and.
     .     VarType ( :10 ) .ne. 'indefinite' ) then
         write ( stderr, * )
         write ( stderr, 902 ) VarName
         ierr = ODS_BadInt
         return

      end if

*     Get pointer data from tables stored in common
*     ---------------------------------------------
      call ODS_GetP    ( id,        'read',
     .                   julian_day, syn_hour,
     .                   synbeg,     synlen,  ierr )
      if ( ierr .ne. NCNoErr ) return

*     Get NetCDF file and variable ids
*     --------------------------------
      nc_id = ncid ( id )
      varid = NCVID ( nc_id, VarName, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Set nval and nval_max
*     ---------------------
      nval_max = nval
      nval     = synlen

*     Return with an error message if the number of
*     observations to be read exceeds the user-defined maximum
*     --------------------------------------------------------
      if ( nval .gt. nval_max ) then
         write ( stderr, * )
         write ( stderr, 903 ) nval_max, nval
         ierr = ODS_DimErr
         return

      end if  

*     Get the length of each batch of observation data
*     in an ODS file
*     ------------------------------------------------
      call NCVINQ ( nc_id,    varid,
     .              Var_Name, Var_Type,
     .              NVDims,   DimIDs,
     .              NVatts,   ierr ) 
      if ( ierr .ne. NCNoErr ) return
      DimID       = DimIDs ( 1 )
      call NCDINQ ( nc_id,    DimID,
     .              DimName,  batch_len, ierr )
      if ( ierr .ne. NCNoErr ) return

*     The remaining code uses the ODS IO routines to read
*     the data in batches.  Each batch contains a subset
*     of the values to be read.  The maximum size of
*     each subset is batchlen as defined in the header
*     file, ods_hdf.h.  The values are being processed in
*     this manner in order to prevent the execution of
*     the NetCDF library routines from being
*     prohibitively slow.
*     ---------------------------------------------------

*     Determine the indices of the first and last batches
*     of data to be read from the ODS file
*     ---------------------------------------------------
      first_batch = ( synbeg - 1 ) / batch_len + 1
      last_batch  = ( synbeg + synlen - 2 ) / batch_len + 1

*     Read the first batch of values
*     ------------------------------
      first_ival  =   mod ( synbeg - 1, batch_len ) + 1
      first_nval  =   min ( synlen, batch_len - first_ival + 1 )
      start ( 1 ) =   first_ival
      start ( 2 ) =   first_batch
      count ( 1 ) =   first_nval
      count ( 2 ) =   1
      call ODS_NCVGTI ( nc_id, varid, start, count, values, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Read data not in the first and last batch of values
*     ---------------------------------------------------
      ival = first_nval + 1
      do 10, ibatch = first_batch + 1, last_batch - 1
         start ( 1 ) = 1
         start ( 2 ) = ibatch
         count ( 1 ) = batch_len
         count ( 2 ) = 1
         call ODS_NCVGTI ( nc_id, varid, start, count,
     .                     values ( ival ), ierr )
         if ( ierr .ne. NCNoErr ) return
         ival        = ival + batch_len
 10   continue

*     Read the last batch of values, if necessary
*     -------------------------------------------
      start ( 1 ) =   1
      start ( 2 ) =   last_batch
      count ( 1 ) =   mod ( synbeg + synlen - 2, batch_len ) + 1
      count ( 2 ) =   1
      if ( last_batch .gt. first_batch )
     .   call ODS_NCVGTI ( nc_id, varid, start, count,
     .                     values ( ival ), ierr )
      if ( ierr .ne. NCNoErr ) return

      return
*     ------

 901  format ( /, ' ODS_GetI: File handle id number does not ',
     .         /, '           correspond to an opened ODS file' )
 902  format ( /, ' ODS_GetI: Inappropriate variable for this ',
     .         /, '           routine.  Variable name is: ', a ) 
 903  format ( /, ' ODS_GetI: The number of values to be read ',
     .         /, '           exceeds the user defined maximum',
     .         /, '           of ', i10, '. Number of values ',
     .         /, '           to be read = ', i10,'.' )

      end
