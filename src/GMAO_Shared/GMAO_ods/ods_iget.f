

*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:      ODS_IGet --- Returns an integer ODS file parameter
! 
! !DESCRIPTION: 
!    \label{ODS:IGet}
!     Returns the integer value of the user-selected parameter of
!     a NetCDF file ( identified by the character string, parm\_name,
!     and the file handle, id ).  The file parameters are the
!     NetCDF dimensions and the NetCDF attributes as defined in the
!     routine, ODS\_Create.  These parameters are identified by the
!     parameter name, parm\_name, and the ODS file handle, id.  The
!     parameter name must be consistent with the notation of network
!     Common Data form Language (CDL) as described in the Rew et al.
!     (1993).  For this routine, if the parameter is an NetCDF
!     attribute, then parm\_name must consist of the NetCDF variable
!     name followed by the delimiter character, ':', and then the
!     attribute name with no blanks.  If the attribute is global,
!     then the first character in parm\_name must be the delimiter.
!     If the parameter is a NetCDF dimension, then parm\_name must
!     not contain any delimiters.  Examples for parm\_name:
!\begin{verbatim}
!
!           nkt           is the dimension, 'nkt'
!           kt:valid_max  is the attribute, 'valid_max', for
!                         the variable, 'kt'
!           :type         is the global attribute, type
!
!\end{verbatim}
!
!     Note: For a list of error codes, see Table~\ref{tab:errors}.
!
!    \bigskip {\bf Reference:}
!    \begin{description}
!    \item Rew, Russ, Glenn Davis and Steve Emmerson, 1993:
!        {\em NetCDF User's Guide}, Unidata Program Center,
!        University Corporation for Atmospheric Research,
!        National Science Foundation, Boulder, CO.
!    \end{description}
!
! !INTERFACE:
!
      subroutine ODS_IGet ( id, parm_name, parm_values, ierr )
!
! !INPUT PARAMETERS:
      implicit   NONE
      integer    id              ! ODS file handle
      character  parm_name * (*) ! The name of the NetCDF parameter
                                 !   to be obtained.  The case of
                                 !   each letter is significant.
!
! !OUTPUT PARAMETERS: 
      integer    parm_values (*) ! The parameter values
      integer    ierr            ! Error code. If non-zero, an 
                                 !  error has occurred. For a list
                                 !  of possible values, see the
                                 !  description section of this
                                 !  prologue.
!
!     NOTES: This routine is designed to extract NetCDF 
!            attributes of limited length and assumes that
!            the calling routine has already predetermined
!            the number the values that this routine will
!            extract.  The names of some useful parameters:
!
!            nkt                       maximum number of GEOS/DAS 
!                                       data types  This dimension
!                                       is named ktmax for versions
!                                       earlier than 2.00
!            nkx                       maximum number of GEOS/DAS
!                                       data sources.  This dimension
!                                       is named kxmax for versions
!                                       earlier than 2.00
!            nqcx                      maximum number of possible 
!                                       values for the quality control
!                                       exclusion mark (qcexcl).  This
!                                       variable is available only for
!                                      versions 2.00 or later
!            syn_beg:first_julian_day  first julian day on file
!                                       to the file.
!            syn_beg:latest_julian_day latest julian day for which
!                                       data exists
!            syn_beg:latest_synoptic_hour
!                                      latest julian hour for which
!                                       data exists
!
! !SEE ALSO: 
!     ODS_RGet() gets the floating point value of a user-selected
!                NetCDF file parameter
!     ODS_CGet() gets the character string of a user-selected
!                NetCDF file parameter
!     ODS_NGet() get the number of observation reports for a given
!                julian day and synoptic hour
!
! !REVISION HISTORY: 
!     15Apr1996   Redder   Original version
!     13Apr1998   Redder   Made routine to be an interface routine
!                          in ODS version 2.00.  Added checks made
!                          to input arguments.
!     01Nov1999   Redder   Revised code to prevent subscript errors
!                          in character strings
!     19Nov1999   Redder   Added a latex label in and moved the
!                          subroutine statement into the prologue.
!                          Modified the comments for the return status
!                          code.
!     06Dec1999   Redder   Corrections to the documentation in the
!                          prologue.
!     16Feb2000   Todling  Rename stdio.h to ods_stdio.h
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
      integer  dimid       ! dimension id
      character * ( MAXNCNAM )
     .         DimName     ! dimension name

*     NetCDF variable information
*     ---------------------------
      integer  nc_id       ! file id
      integer  varid       ! variable id
      character * ( MAXNCNAM )
     .         VarName     ! variable name

*     NetCDF attribute information
*     ----------------------------
      character * ( MAXNCNAM )
     .         AttName     ! attribute name
      integer  AttLen      ! attribute

*     Other variables
*     -----------------
      integer  delimiter_loc
     .                     ! location of delimiter in  
      integer  ParmNameSz  ! number of character in the parameter name
      integer  parm_value  ! temporary storage for a 
                           !   single value in parm_values

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

*     Some NetCDF file parameters may change during program
*     execution and therefore are stored in common.  If the
*     calling routine requests any of these parameters, then
*     extract them directly from common storage rather by
*     calling internal routines rather than use NetCDF routines.
*     ----------------------------------------------------------

      call ODS_GetAttP ( id, parm_name, parm_value, ierr )

*     If the NetCDF file parameter is stored in common storage,
*     then extract the value and return
*     ---------------------------------------------------------
      if ( ierr .eq. NCNoErr ) then
         parm_values ( 1 ) = parm_value
         return

      else ! Else reset the return error code and continue
*     ----------------------------------------------------
         ierr = NCNoErr
         continue

      end if

*     If the requested parameter is a NetCDF dimension ( i.e.
*     no delimiter character is found in the parameter name
*     string ) then ...
*     -------------------------------------------------------
      if      ( delimiter_loc .eq. Not_Found ) then

*        extract the NetCDF dimension size
*        ---------------------------------
         DimName = parm_name
         dimid   = NCDID ( nc_id,  DimName,     ierr )
         if ( ierr .ne. NCNoErr )  return
         call ncdinq     ( nc_id,  dimid,       DimName,
     .                             parm_values, ierr )
         if ( ierr .ne. NCNoErr )  return

*     If the requested parameter is a NetCDF global attribute
*     ( i.e. the delimiter is the first character in the
*     parameter name string ) then ...
*     -------------------------------------------------------
      else if ( delimiter_loc .eq. 1         ) then

*        extract the NetCDF global attribute
*        -----------------------------------
         varid   = NCGLOBAL
         AttName = parm_name ( 2:ParmNameSz )
         call ODS_NCAGTI ( nc_id,  varid,       AttName,
     .                     AttLen, parm_values, ierr )
         if ( ierr .ne. NCNoErr  )  return

*     If the requested parameter is a NetCDF nonglobal attribute
*     ( i.e. the delimiter is located after the first character
*     in the parameter name string ) then ...
*     ----------------------------------------------------------
      else if ( delimiter_loc .gt. 1         ) then

*        extract the attribute for the NetCDF variable, VarName
*        ------------------------------------------------------
         VarName = parm_name (  :delimiter_loc - 1 )
         AttName = parm_name (   delimiter_loc + 1 : ParmNameSz )
         varid   = ncvid ( nc_id,  VarName,     ierr )
         if ( ierr .ne. NCNoErr )  return
         call ODS_NCAGTI ( nc_id,  varid,       AttName,
     .                     Attlen, parm_values, ierr )
         if ( ierr .ne. NCNoErr )  return
      end if

 901  format ( /, ' ODS_IGet: File handle id number does not ',
     .         /, '           correspond to an opened ODS file' )
 902  format ( /, ' ODS_IGet: Parameter name is inappropriate. ',
     .         /, '           Name = ', a )

      return
      end
