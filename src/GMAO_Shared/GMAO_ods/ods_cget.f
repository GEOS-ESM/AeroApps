

*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_CGet --- Returns a character ODS file parameter
! 
! !DESCRIPTION: 
!    \label{ODS:CGet}
!     Returns the character string of the user-selected parameter
!     of a NetCDF file ( identified by the character string,
!     parm\_name, and the file handle, id ).  If the input string is
!     larger than is necessary, then trailing blanks is inserted
!     into the unused portion of the string.  The file parameters are
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
!     delimiter.  Examples for parm\_name:
!\begin{verbatim}
!
!           kt:valid_max  is the attribute, 'valid_max', for
!                         variable, 'kt'
!           :type         is the global attribute, type
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
      subroutine ODS_CGet ( id, parm_name, parm_string, ierr )
!
! !INPUT PARAMETERS:
      implicit   NONE
      integer    id                ! ODS file handle
      character  parm_name   * (*) ! The name of the parameter
                                   !   to be obtained.  The case of
                                   !   each letter is significant.
!
! !OUTPUT PARAMETERS: 
      character  parm_string * (*) ! The parameter values
      integer    ierr              ! Error code. If non-zero, an 
                                   !  error has occurred. For a list
                                   !  of possible values, see the
                                   !  description section of this
                                   !  prologue.
!
!     NOTES: If insufficient space is allocated for parm_string,
!            the parameter string is truncated. This routine is
!            not appropriate for extracting the dimension sizes
!            To obtain the dimension sizes, use the routine,
!            ODS_IGet.  Some of the parameters including
!           ':filename' and ':mode' (read or write mode) are
!            stored in common rather than the netcdf file.
!            The names of some useful parameters:
!
!              :filename      name of ODS file
!              :mode          mode = 'w' open for write and read
!                             mode = 'r' open for read only
!              :version       version tag
!              :history       history tag
!              :type          file type ("pre-analysis" or 
!                                        "post-analysis")
!
!
! !SEE ALSO:
!     ODS_IGet() Gets the integer of a user-selected NetCDF file
!                parameter
!     ODS_RGet() Gets the floating point value of a user-selected
!                NetCDF file parameter
!     ODS_NGet() Gets the number of observation reports for a given
!                julian day and synoptic hour
!
! !REVISION HISTORY: 
!     05Apr1996   Redder   Original version
!     13Apr1998   Redder   Made routine to be an interface routine
!                          in ODS version 2.00.  Added checks made
!                          to input arguments.
!     15Mar1999   Redder   Added code to retrieve the ODS type
!     01Nov1999   Redder   Revised code to prevent subscript errors
!                          in character strings
!     19Nov1999   Redder   Added a latex label in and moved the
!                          subroutine statement into the prologue.
!                          Modified the comments for the return status
!                          code.
!     06Dec1999   Redder   Corrections to the documentation in the
!                          prologue.
!     16Feb2000   Todling  Rename stdio.h to ods_stdio.h
!
!EOP
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

*     Function referenced
*     -------------------
      integer     ODS_StrSize ! Returns the character string length
                              !   excluding trailing blanks
      character * (max_strlen)
     .            ODS_Type    ! Returns the ODS type

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
      integer  AttLen         ! attribute

*     Other variables
*     -----------------
      integer  delimiter_loc
     .                        ! location of delimiter in  
      integer  ParmNameSz     ! number of character in the
                              !   parameter name

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
      if ( delimiter_loc .ge. ParmNameSz ) then
         write ( stderr, 902 ) parm_name ( : ParmNameSz )
         ierr = NCENOTVR
         return

      end if

*     For variables found in common ...
*     ---------------------------------

*     extract the file name
*     ---------------------
      if ( ParmNameSz        .eq. 9        ) then
      if ( parm_name ( : 9 ) .eq. ':filename' ) then
         parm_string = filenames ( id )
         return

      end if
      end if

*     extract the mode 
*     ----------------
      if ( ParmNameSz        .eq. 5        ) then
      if ( parm_name ( : 5 ) .eq. ':mode'  ) then
         if ( IOMode ( id )  .eq. NCWRITE  )
     .      parm_string = 'w'
         if ( IOMode ( id )  .eq. NCNOWRIT )
     .      parm_string = 'r'
         return

      end if
      end if

*     ... with special handling required for retrieval
*     ------------------------------------------------

*     extract the ODS file type
*     -------------------------
      if ( ParmNameSz        .eq. 5        ) then
      if ( parm_name ( : 5 ) .eq. ':type'  ) then
         parm_string = ODS_Type ( id, ierr )
         return

      end if
      end if

*     If the delimiter string is not found, then the program
*     will assume that the calling program is requesting a
*     NetCDF dimension.  If no delimiter is found, then ...
*     ------------------------------------------------------
      if      ( delimiter_loc .eq. Not_Found ) then

*        print error message and return
*        ------------------------------
         ierr = NCEInVal
         write ( stderr, 903 ) parm_name
         return

*     If the requested parameter is a NetCDF global attribute
*     ( i.e. the delimiter is the first character in the
*     parameter name string ) then ...
*     -------------------------------------------------------
      else if ( delimiter_loc .eq. 1         ) then

*        ... or extract the NetCDF global attribute
*        ------------------------------------------
         varid       = NCGLOBAL
         AttName     = parm_name ( 2:ParmNameSz )
         call ODS_NCAGTC ( nc_id,  varid,       AttName,
     .                     AttLen, parm_string, ierr )
         parm_string = parm_string ( : AttLen )
         if ( ierr .ne. NCNoErr  )  return

*     If the requested parameter is a NetCDF nonglobal attribute
*     ( i.e. the delimiter is located after the first character
*     in the parameter name string ) then ...
*     ----------------------------------------------------------
      else if ( delimiter_loc .gt. 1         ) then

*        extract the NetCDF attribute for the NetCDF 
*        variable, VarName
*        -------------------------------------------
         VarName     = parm_name (  :delimiter_loc - 1 )
         AttName     = parm_name (   delimiter_loc + 1 : ParmNameSz )
         varid       = ncvid ( nc_id,  VarName,     ierr )
         if ( ierr .ne. NCNoErr )  return
         call ODS_NCAGTC ( nc_id,  varid,       AttName,
     .                     AttLen, parm_string, ierr )
         parm_string = parm_string ( : AttLen )
         if ( ierr .ne. NCNoErr )  return
      end if

      return
*     ------

 901  format ( /, ' ODS_CGet: File handle id number does not ',
     .         /, '           correspond to an opened ODS file' )
 902  format ( /, ' ODS_CGet: Parameter name is inappropriate. ',
     .         /, '           Name = ', a )
 903  format ( /, ' ODS_CGet: No delimiter string was found ',
     .         /, '           in the name of the parameter, ',
     .         /, '           and, therefore, this routine ',
     .         /, '           assumes that the dimension ',
     .         /, '           size is requested.  This ',
     .         /, '           routine is inappropriate for ',
     .         /, '           obtaining this parameter. ',
     .         /, '           Parameter name is ', a )

      end
