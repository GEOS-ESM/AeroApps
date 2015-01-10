
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_PutR() --- Writes data in native real format to file
! 
! !DESCRIPTION:
!    \label{ODS:PutR}
!     Stores the contents of a {\tt real} variable for a full
!     synoptic hour on an opened ODS file.\\
!
!  \noindent {\bf Notes: }\\ 
!  \begin{enumerate}
!  \item For a list of error codes, see Table~\ref{tab:errors}.
!  \item This routine does perform checks to verify whether
!     the range of values on input is consistent with the internal
!     variable type in the ODS file.  For example, a real variable
!     here could be stored as a two byte variable on file.  
!     Therefore, all elements in the array to be stored as two byte
!     integers should have values between -32,767 and 32,767.
!     If the software does detect a number inconsistent with the
!     variable type, then the routine returns with an error
!     message without saving at least some of the values.  This
!     check prevents overflows from occurring which would produce
!     unexpected results and probably abnormal program termination.  
!
!  \item The input values should be within the range specified in 
!     Table~\ref{tab:ObAtt} or by the ODS file annotated header of
!     the ODS file produced with the MFHDF utility, ncdump.  If no
!     range is specified, then the maximum and minimum values for
!     the internal variable types are as follows:
!
!  \begin{verbatim}
!       type                      minimum value    maximum value
!       ----                      -------------    -------------
!       byte                                  0              255
!       short (2 byte) integer          -32,767           32,767
!       long  (4 byte) integer   -2,147,483,643    2,147,483,643
!       real  (4 bytes )               -3.40e38          3.40e38
!  \end{verbatim}
!
!  \item Once data from the specified date and hour have been
!     written, the routine will allow overwriting of data
!     from a previous date and hour.  The only restriction
!     is that the number to be written does not exceed the
!     space already allocated in the ODS file.  Only the first
!     {\tt nval} values will be overwritten.
!
!  \end{enumerate}
!
! !INTERFACE: 
!
      subroutine ODS_PutR ( id, VarName, julian_day, syn_hour,
     .                          nval,    values,     ierr )
!
! !INPUT PARAMETERS:
      implicit        NONE
      integer         id               ! ODS file handle as returned
                                       !  from ODS_Create() or
                                       !  ODS_Open().
      character * (*) VarName          ! name of variable
      integer         julian_day       ! Julian Day.  Use the function
                                       !  ODS_Julian to obtain this
                                       !  number.
      integer         syn_hour         ! hour of synoptic time since
                                       !  0:00 GMT (e.g., 0, 6, 12, 18)
      integer         nval             ! number of values for this
                                       !  synoptic time
      real            values ( nval )  ! values to be written
!
! !OUTPUT PARAMETERS:
      integer         ierr             ! Error code. If non-zero, an 
                                       !  error has occurred. For a list
                                       !  of possible values, see the
                                       !  description section of this
                                       !  prologue.
!
! !SEE ALSO:
!     ODS_PutC()    write a list of character  strings to   an ODS file
!     ODS_PutI()    write a list of integer    values  to   an ODS file
!     ODS_GetC()    get   a list of characters strings from an ODS file
!     ODS_GetI()    get   a list of integer    values  from an ODS file
!     ODS_GetR()    get   a list of real       values  from an ODS file
!     ODS_Julian()  convert the "calendar" date to the Julian day
!     ods_hdf.h     include file defining internal constants
!                    and global variables, and setting up data
!                    structures
!
! !REVISION HISTORY: 
!     02Oct1995   da Silva    Specification.
!     17May1996   Redder      Original code
!     13Apr1998   C. Redder   Updated documentation to reflect changes
!                             made in other routines to create ODS
!                             version 2.00
!     26May1999   C. Redder   Fixed bug in handling the case when nval
!                             is zero or less.
!     19Nov1999   C. Redder   Added a latex label in and moved the
!                             subroutine statement into the prologue.
!                             Modified the comments for the return status
!                             code.  Removed references to the old version
!                             of the ODS documentation.
!     06Dec1999   C. Redder   Corrections to the documentation in the
!                             prologue.
!     30Nov2001   C. Redder   Fixed bug that prevents last batch from 
!                             being processed twice when nval = 0
!     16Sep2002   C. Redder   Updated prologue and added error handling
!                             for call to NetCDF routine, NCDINQ
!     17Oct2002   C. Redder   Fixed bug that can causes an array index
!                             error when only part of the data for the 
!                             latest date and time on file is being
!                             overwritten.
!
!EOP
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Functions referenced
*     --------------------
      character * ( 10 )
     .         ODS_VarType ! recommended FORTRAN variable
                           !   type for data
      integer  ODS_JulHr   ! Julian hour using, as the reference,
                           !   the synoptic hour and julian day
                           !   that corresponds to the last block
                           !   of data written to file

*     Arguments for NetCDF library routines
*     -------------------------------------
      integer  nc_id       ! file id
      integer  varid       ! variable id
      integer  start  ( 3 )! starting file index for hyperslab
      integer  count  ( 3 )! length of edges in hyperslab
      character * ( MAXNCNAM )
     .         Var_Name    ! variable name ( returned from NetCDF
                           !   routine NCVINQ )
      integer  Var_Type    ! variable type
      integer  NVDims      ! number of dimensions
      integer  DimIDs ( 3 )! dimension ids
      integer  DimID       ! dimension id
      integer  NVAtts      ! number of attributes
      character * ( MAXNCNAM )
     .         DimName     ! dimension name ( returned from NetCDF
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
      character * ( 10 ) 
     .         VarType     ! recommended FORTRAN variable type for
                           !   data
      integer JulHr        ! Julian hour corresponding to julian_day
                           !   and syn_hour
      integer JulHr_Latest ! Julian hour corresponding to the synoptic
                           !   date and hour of the last block of data
                           !   that was written to file 
      parameter 
     .      ( JulHr_Latest = 0 )
                           ! Set to zero for clarity and convenience


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

*     If IO mode is not set to write,
*     then exit with an error message
*     -------------------------------
      if ( IOMode ( id ) .ne. NCWRITE ) then
         write ( stderr, 902 ) filenames ( id )
         ierr = NCEPerm
         return

      end if

*     Check to determine if variable is appropriate for this routine
*     --------------------------------------------------------------
      VarType    = ODS_VarType ( VarName )
      if ( VarType ( :4  ) .ne. 'real'     .and.
     .     VarType ( :10 ) .ne. 'indefinite' ) then
         write ( stderr, 903 ) VarName
         ierr = ODS_BadInt
         return

      end if

*     Return with an error message if the number of values
*     to be written is a negative integer
*     ----------------------------------------------------
      if ( nval .lt. 0 ) then
         ierr = ODS_DimErr
         write ( stderr, 904 ) VarName
         return

      end if

*     Calculate the Julian hour ( using, as the reference,
*     the date and hour of the last block of data that
*     was written to file )
*     ---------------------------------------------------
      JulHr = ODS_JulHr ( id, julian_day, syn_hour, ierr )

*     For synoptic times earlier than the synoptic time
*     for the last data block written to file ...
*     ------------------------------------------------------
      if      ( JulHr .lt. JulHr_Latest ) then

*        Get pointer data from tables stored in common
*        ---------------------------------------------
         call ODS_GetP    ( id,        'write',
     .                      julian_day, syn_hour,
     .                      synbeg,     synlen, ierr )
         if ( ierr .ne. NCNoErr ) return

*        Return with an error message if the number of values
*        exceeds the space already allocated in the ODS file
*        ----------------------------------------------------
         if ( nval .gt. synlen ) then
            write ( stderr, 905 )
            ierr  = NCEInVal
            return

         end if

         synlen  = nval ! Set syn_len to the user specified value

*     For synoptic times at the same time for the last data
*     block written to file ...
*     -----------------------------------------------------
      else if ( JulHr .eq. JulHr_Latest ) then

*        Get pointer data from tables stored in common
*        ---------------------------------------------
         call ODS_GetP    ( id,        'write',
     .                      julian_day, syn_hour,
     .                      synbeg,     synlen, ierr )
         if ( ierr .ne. NCNoErr ) return

*        Return with an error message if the number of values
*        exceeds the space already allocated in the ODS file
*        ----------------------------------------------------
         if ( nval .gt. synlen ) then
            write ( stderr, 906 )
            ierr  = NCEInVal
            return

         end if

         synlen  = nval ! Set syn_len to the user specified value

*     For synoptic times later than the time for the
*     last data block written to file ...
*     -------------------------------------------------
      else if ( JulHr .gt. JulHr_Latest ) then

*        Update pointer tables ( tables stored in common )
*        -------------------------------------------------
         call ODS_UpdateP ( id,        'no_option',
     .                      julian_day, syn_hour,
     .                      nval,       ierr )
         if ( ierr .ne. NCNoErr ) return

*        Get pointer data from tables stored in common
*        ---------------------------------------------
         call ODS_GetP    ( id,        'write',
     .                      julian_day, syn_hour,
     .                      synbeg,     synlen, ierr )
         if ( ierr .ne. NCNoErr ) return

      end if

*     Nothing more to do if no values are to be written
*     -------------------------------------------------
      if ( nval .le. 0 ) return

*     Get NetCDF file and variable ids
*     --------------------------------
      nc_id = ncid ( id )
      varid = NCVID ( nc_id, VarName, ierr )
      if ( ierr .ne. NCNoErr ) return

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

*     The remaining code uses the ODS IO routines to write
*     the data in batches.  Each batch contains a subset
*     of the values to be written.  The maximum size of
*     each subset is batchlen as defined in the header
*     file, ods_hdf.h.  The values are being processed in
*     this manner in order to prevent the execution of
*     the NetCDF library routines from being
*     prohibitively slow.
*     ---------------------------------------------------

*     Determine the indices of the first and last batches
*     of data to be written to the ODS file
*     ---------------------------------------------------
      first_batch = ( synbeg - 1 ) / batch_len + 1
      last_batch  = ( synbeg + synlen - 2 ) / batch_len + 1

*     Write the first batch of values
*     -------------------------------
      first_ival  =   mod ( synbeg - 1, batch_len ) + 1
      first_nval  =   min ( synlen, batch_len - first_ival + 1 )
      start ( 1 ) =   first_ival
      start ( 2 ) =   first_batch
      count ( 1 ) =   first_nval
      count ( 2 ) =   1
      call ODS_NCVPTR ( nc_id, varid, start, count, values, ierr )
      if ( ierr .ne. NCNoErr ) return

*     Write data not in the first and last batch of values,
*     if necessary
*     -----------------------------------------------------
      ival = first_nval + 1
      do 10, ibatch = first_batch + 1, last_batch - 1
         start ( 1 ) = 1
         start ( 2 ) = ibatch
         count ( 1 ) = batch_len
         count ( 2 ) = 1
         call ODS_NCVPTR ( nc_id, varid, start, count,
     .                     values ( ival ), ierr )
         if ( ierr .ne. NCNoErr ) return
         ival        = ival + batch_len
 10   continue

*     Write the last batch of values, if necessary
*     --------------------------------------------
      start ( 1 ) =   1
      start ( 2 ) =   last_batch
      count ( 1 ) =   mod ( synbeg + synlen - 2, batch_len ) + 1
      count ( 2 ) =   1
      if ( last_batch .gt. first_batch )
     .   call ODS_NCVPTR ( nc_id, varid, start, count,
     .                     values ( ival ), ierr )
      if ( ierr .ne. NCNoErr ) return

      return
*     ------

 901  format ( /, ' ODS_PutR: File handle id number does not ',
     .         /, '           correspond to an opened ODS file' )
 902  format ( /, ' ODS_PutR: Cannot write to file unless mode ',
     .         /, '           is set to write.  Filename is: ', a )
 903  format ( /, ' ODS_PutR: Inappropriate variable for this ',
     .         /, '           routine.  Variable name is: ', a ) 
 904  format ( /, ' ODS_PutR: The number of values to be written ',
     .         /, '           must be a non-negative integer. ',
     .         /, '           Variable name is: ', a ) 
 905  format ( /, ' ODS_PutR: Cannot increase the size of a ',
     .         /, '           block of data that corresponds ',
     .         /, '           to a date and time earlier than ',
     .         /, '           the latest date and time for ',
     .         /, '           which there is data. ' )
 906  format ( /, ' ODS_PutR: Must use ODS_Append to increase the '
     .         /, '           size of a block of data that corre-'
     .         /, '           sponds to the latest date and time '
     .         /, '           for which there is data. ' )

      end
