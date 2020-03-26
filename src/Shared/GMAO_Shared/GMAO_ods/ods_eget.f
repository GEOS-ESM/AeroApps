
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_EGet --- Tests for the existence of an ODS file parameter
! 
! !DESCRIPTION: 
!    \label{ODS:EGet}
!     Returns true if a the user-selected parameter of a NetCDF file
!     ( identified by the character string, parm\_name, and the file
!     handle, id ) exist and false otherwise.  The file parameters are
!     the NetCDF dimensions and the NetCDF attributes as defined in
!     the routine, ODS\_Create.  These parameters are identified by
!     the parameter name, parm\_name, and the ODS file handle, id.
!     The parameter name must be consistent with the notation of
!     network Common Data form Language (CDL) as described in the
!     Rew et al. (1993).  For this routine, if the parameter is a
!     NetCDF attribute, then parm\_name must consist of the NetCDF
!     variable name followed by the delimiter character, ':', and
!     then the attribute name with no blanks.  If the attribute is
!     global, then the first character in parm\_name must be the
!     delimiter.  If the parameter is a NetCDF variable, then the
!     variable must be followed by the delimiter character.  If the
!     parameter is a NetCDF dimension, then no delimiter character
!     must be present.  Examples for parm\_name:
!\begin{verbatim}
!
!           kt:valid_max  is the attribute, 'valid_max', for
!                         variable, 'kt'
!           kt:           is the varialbe kt
!           :type         is the global attribute, type
!           nsyn          if the dimension nsyn.
!
!\end{verbatim}
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
!  \bigskip {\bf Reference:}
!  \begin{description}
!  \item Rew, Russ, Glenn Davis and Steve Emmerson, 1993:
!        {\em NetCDF User's Guide}, Unidata Program Center,
!        University Corporation for Atmospheric Research,
!        National Science Foundation, Boulder, CO.
!  \end{description}
!
! !INTERFACE:
!
      logical function ODS_EGet ( id, parm_name, ierr )
!
! !INPUT PARAMETERS:
      implicit   NONE
      integer    id                ! ODS file handle
      character  parm_name   * (*) ! The name of the parameter
                                   !   to be obtained.  The case of
                                   !   each letter is significant.
!
! !OUTPUT PARAMETERS: 
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
!
! !REVISION HISTORY: 
!     20Sep2002   Redder   Original version.
!EOP
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Function referenced
*     -------------------
      integer     ODS_StrSize ! Returns the character string length
                              !   excluding trailing blanks

*     Additional variables for ...
*     ----------------------------

*     NetCDF dimension information
*     ----------------------------
      integer  dimid          ! dimension id
      character * ( MAXNCNAM )
     .         DimName        ! dimension name

*     NetCDF variable information
*     ---------------------------
      integer  nc_id          ! file id
      integer  varid          ! variable id
      character * ( MAXNCNAM )
     .         VarName        ! variable name

*     NetCDF attribute information
*     ----------------------------
      character * ( MAXNCNAM )
     .         AttName        ! attribute name
      integer  AttLen         ! attribute length
      integer  AtType         ! attribute type

*     Other variables
*     -----------------
      integer  delimiter_loc
     .                        ! location of delimiter in  
      integer  ParmNameSz     ! number of character in the
                              !   parameter name
      integer     ncopts      ! NetCDF error handling options

*     Set error code to default
*     -------------------------
      ierr  = NCNoErr

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

*     Determine the length of the NetCDF file parameter name
*     ------------------------------------------------------
      ParmNameSz = ODS_StrSize ( parm_name )
      if ( ParmNameSz .le. 0 ) then
         write ( stderr, 902 ) '(blank name)'
         ierr = NCENOTVR
         return

      end if

*     Locate the delimiter character, ':', in the NetCDF file
*     parameter name.
*     -------------------------------------------------------
      delimiter_loc = index ( parm_name, ':' )
      if ( ParmNameSz .eq. delimiter_loc .and.
     .     ParmNameSz .eq. 1 ) then
         write ( stderr, 902 ) parm_name ( : ParmNameSz )
         ierr = NCENOTVR
         return

      end if

*     Save current options of error handling and turn off
*     error messages and set errors to be non-fatal.
*     ---------------------------------------------------
      call NCGOPT     ( ncopts )
      call NCPOPT     ( 0 )

*     If the requested parameter is a NetCDF dimension ( i.e.
*     no delimiter character is found in the parameter name
*     string ) then ...
*     -------------------------------------------------------
      if      ( delimiter_loc .eq. Not_Found ) then

*        ... determine if the dimension exists
*        -------------------------------------
         DimName     = parm_name
         dimid       = ncdid ( nc_id,  DimName,     ierr )
         ODS_EGet    = .true.
         call  NCPOPT ( ncopts )           ! Return to previous error handling
         if    ( ierr .ne. NCNoErr  ) then ! If error status is ...
            ODS_EGet = .false.
            if ( ierr .eq. NCEBADD  ) then ! ... variable not found then
               ierr  = NCNoErr             ! ... reset error status to valid
            else                           ! ... otherwise recall the NetCDF 
               dimid = ncdid ( nc_id,  DimName, ierr )
            end if                         !     routine to perform the 
         end if                            !     error handling

*     If the requested parameter is a NetCDF global attribute
*     ( i.e. the delimiter is the first character in the
*     parameter name string ) then ...
*     -------------------------------------------------------
      else if ( delimiter_loc .eq. 1         ) then

*        determine if the global attribute exists
*        ----------------------------------------
         varid    = NCGLOBAL
         AttName  = parm_name ( 2:ParmNameSz )
         call NCAINQ ( nc_id,  varid,  AttName,
     .                  AtType, AttLen, ierr )
         ODS_EGet = .true.
         call  NCPOPT ( ncopts )           ! Return to previous error handling
         if    ( ierr .ne. NCNoErr  ) then ! If error status is ...
            ODS_EGet = .false.
            if ( ierr .eq. NCENoAtt ) then ! ... attribute not found then
               ierr  = NCNoErr             ! ... reset error status to valid
            else                           ! ... otherwise recall the NetCDF 
               call NCAINQ ( nc_id,  varid,  AttName,
     .                       AtType, AttLen, ierr )
            end if                         !     routine to perform the 
         end if                            !     error handling

*     If the requested parameter is a NetCDF variable
*     ( i.e. the delimiter is the last character in the
*     parameter name string ) then ...
*     -------------------------------------------------
      else if ( delimiter_loc .eq. ParmNameSz ) then

*        determine if the variable exists ...
*        ------------------------------------
         VarName  = parm_name ( 1:ParmNameSz - 1 )
         varid    = NCVID ( nc_id, VarName, ierr )
         ODS_EGet = .true.
         call  NCPOPT ( ncopts )           ! Return to previous error handling
         if    ( ierr .ne. NCNoErr  ) then ! If error status is ...
            ODS_EGet = .false.
            if ( ierr .eq. NCENotVr ) then ! ... variable not found then
               ierr  = NCNoErr             ! ... reset error status to valid
            else                           ! ... otherwise recall the NetCDF 
               call NCAINQ ( nc_id,  varid,  AttName,
     .                       AtType, AttLen, ierr )
            end if                         !     routine to perform the 
         end if                            !     error handling

*     If the requested parameter is a NetCDF nonglobal attribute
*     ( i.e. the delimiter is located in the middle of the 
*     parameter name string ) then ...
*     ----------------------------------------------------------
      else if ( delimiter_loc .gt. 1         ) then

*        get the NetCDF variable id
*        --------------------------
         VarName = parm_name (  :delimiter_loc - 1 )
         varid   = ncvid ( nc_id,  VarName,     ierr )

         if ( ierr .ne. NCNoErr ) then     ! If error status is ...
            ods_eget = .false.
            call  NCPOPT ( ncopts )        ! Recall the NetCDF routine to
            varid    = NCVID ( nc_id, VarName, ierr )
            return
         end if                            !    perform the error handling

*        ... and determine if the NetCDF attribute exists 
*        ------------------------------------------------
         AttName  = parm_name ( 2:ParmNameSz )
         AttName  = parm_name (   delimiter_loc + 1 : ParmNameSz )
         call NCAINQ ( nc_id,  varid,  AttName,
     .                 AtType, AttLen, ierr )

         ODS_EGet = .true.
         call  NCPOPT ( ncopts )           ! Return to previous error handling
         if    ( ierr .ne. NCNoErr  ) then ! If error status is ...
            ODS_EGet = .false.
            if ( ierr .eq. NCENoAtt ) then ! ... attribute not found then
               ierr  = NCNoErr             ! ... reset error status to valid
            else                           ! ... otherwise recall the NetCDF 
               call NCAINQ ( nc_id,  varid,  AttName,
     .                       AtType, AttLen, ierr )
            end if                         !     routine to perform the 
         end if                            !     error handling
      end if

 901  format ( /, ' ODS_EGet: File handle id number does not ',
     .         /, '           correspond to an opened ODS file' )
 902  format ( /, ' ODS_EGet: Parameter name is inappropriate. ',
     .         /, '           Name = ', a )

      end
