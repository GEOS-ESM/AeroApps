
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:      ODS_ParmC --- Gets a NetCDF string parameter for creating files 
! 
! !DESCRIPTION: 
!    \label{ODS:ParmC}
!     Returns an ODS parameter for creating files by subsequent
!     calls to the interface routine, \verb|ODS_Create|.  The
!     returned value is determined by a default value or by a
!     call to the companion routine \verb|ODS_SetParmC| (or the
!     routine \verb|ODS_Create| for the parameters, 'nkt', 'nkx'
!     and 'nqcx').  The parameter name must be consistent with
!     the notation of network Common Data form Language (CDL) as
!     described in the Rew et al. (1993).  For this routine, if
!     the parameter is an NetCDF attribute, then the parameter
!     name must consist of the NetCDF variable name followed by
!     the delimiter character, ':', and then the attribute name
!     with no blanks.  If the attribute is global, then the
!     first character in the name must be the delimiter.  If the
!     parameter is a NetCDF dimension, then the name must not
!     contain any delimters.
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
      character * (*) function ODS_ParmC ( ParmName,
     .                                     DFault_String,
     .                                     ierr )

      implicit NONE
!
! !INPUT PARAMETERS:
      character * (*) ParmName      ! The name of the NetCDF file
                                    !   parameter to be set.  The
                                    !   case of each letter is
                                    !   significant.
      character * (*) DFault_String ! The default string in case the
                                    !   parameter name is not on the
                                    !   list of specified integers
!
! !OUTPUT PARAMETER:
      integer         ierr          ! The returned error code
!
! !SEE ALSO: 
!     ODS_ParmI    ( Sets the integer value of a NetCDF file parameter )
!     ODS_ParmR    ( Sets the floating point value of a NetCDF file
!                    parameter )
!     ODS_SetParmC ( Saves the character string of a NetCDF file parameter )
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!
! !REVISION HISTORY: 
!     15Apr1996   Redder   Origional version
!     25Oct2000   Redder   Made routine callable.
!
!EOP
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'

*     Functions referened
*     -------------------
      integer       ODS_StrSearch ! Search for the entry among a list of 
                                  !   character strings containing the
                                  !   string being searched

*     Other variables
*     ---------------
      integer       list_entry    ! Index number of the list entry
                                  !   returned from the search string
                                  !   function, ODS_StrSearch.  If no
                                  !   string is found then the value
                                  !   is equal to the parameter,
                                  !   Not_Found, as defined in the
                                  !   header file, ods_hdf.h

*     Search for the NetCDF parameter name, ParmName, among the
*     list of file parameters that are being set
*     ---------------------------------------------------------
      list_entry = ODS_StrSearch ( ParmName,
     .                             ParmC_nlist,
     .                             ParmC_list,
     .                             ierr )
      if ( ierr .ne. NCNoErr ) return

*     If the NetCDF file parameter name is not found ...
*     --------------------------------------------------
      if ( list_entry .eq . Not_Found ) then

         ODS_ParmC = DFault_String ! use default string

      else

         ODS_ParmC = ParmC_string ( list_entry ) 
                                   ! use the string corresponding to
                                   !   the index number, list_entry
      end if

      return
      end

*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:      ODS_ParmI --- Gets a NetCDF integer parameter for creating files 
! 
! !DESCRIPTION: 
!    \label{ODS:ParmI}
!     Returns an ODS parameter for creating files by subsequent
!     calls to the interface routine, \verb|ODS_Create|.  The
!     returned value is determined by a default value or by a
!     call to the companion routine \verb|ODS_SetParmI| (or the
!     routine \verb|ODS_Create| for the parameters, 'nkt', 'nkx'
!     and 'nqcx').  The parameter name must be consistent with
!     the notation of network Common Data form Language (CDL) as
!     described in the Rew et al. (1993).  For this routine, if
!     the parameter is an NetCDF attribute, then the parameter
!     name must consist of the NetCDF variable name followed by
!     the delimiter character, ':', and then the attribute name
!     with no blanks.  If the attribute is global, then the
!     first character in the name must be the delimiter.  If the
!     parameter is a NetCDF dimension, then the name must not
!     contain any delimters.  Some useful parameters are as
!     follows:
! \bv
!        ndays - maximum number of days stored in file
!        nsyn  - number of synoptic times per day
! \ev
!
!    \bigskip {\bf Reference:}
!    \begin{description}
!    \item Rew, Russ, Glenn Davis and Steve Emmerson, 1993:
!        {\em NetCDF User's Guide}, Unidata Program Center,
!        University Corporation for Atmospheric Research,
!        National Science Foundation, Boulder, CO.
!    \end{description}
!
!
! !INTERFACE:
      integer function ODS_ParmI ( ParmName, DFault_Val, ierr )
!
! !INPUT PARAMETERS:
      character * (*) ParmName   ! The name of the NetCDF file
                                 !   parameter to be set.  The
                                 !   case of each letter is
                                 !   significant.
      integer         DFault_Val ! The default value in case the
                                 !   parameter name is not on the
                                 !   list of specified integers
!
! !OUTPUT PARAMETER:
      integer         ierr       ! The returned error code
!
! !SEE ALSO: 
!     ODS_ParmR    ( Sets the floating point value of a NetCDF file
!                    parameter )
!     ODS_ParmC    ( Sets the character string of a NetCDF file parameter )
!     ODS_SetParmC ( Saves the character string of a NetCDF file parameter )
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!
! !REVISION HISTORY: 
!     04Sep1996   Redder   Origional version
!     25Oct2000   Redder   Made routine callable.
!
!EOP
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'

*     Functions referened
*     -------------------
      integer       ODS_StrSearch ! Search for the entry among a list of 
                                  !   character strings containing the
                                  !   string being searched

*     Other variables
*     ---------------
      integer       list_entry    ! Index number of the list entry
                                  !   returned from the search string
                                  !   function, ODS_StrSearch.  If no
                                  !   string is found then the value
                                  !   is equal to the parameter,
                                  !   Not_Found, as defined in the
                                  !   header file, ods_hdf.h

*     Search for the NetCDF file parameter name, ParmName, among
*     the list of file parameters that are being set
*     ----------------------------------------------------------
      list_entry = ODS_StrSearch ( ParmName,
     .                             ParmI_nlist,
     .                             ParmI_list,
     .                             ierr )

*     If an error is detected ...
*     ---------------------------
      if ( ierr .ne. NCNoErr ) then

         ODS_ParmI = DFault_Val   ! use default value
         return

      end if

*     If the NetCDF parameter name is not found ...
*     ---------------------------------------------
      if ( list_entry .eq . Not_Found ) then

         ODS_ParmI = DFault_Val   ! use default value

      else

         ODS_ParmI = ParmI_val ( list_entry ) 
                                  ! use the value corresponding to
                                  !   the index number, list_entry
      end if

      return
      end

*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:      ODS_ParmR --- Gets a NetCDF real parameter for creating files 
! 
! !DESCRIPTION: 
!    \label{ODS:ParmR}
!     Returns an ODS parameter for creating files by subsequent
!     calls to the interface routine, \verb|ODS_Create|.  The
!     returned value is determined by a default value or by a
!     call to the companion routine \verb|ODS_SetParmR| (or the
!     routine \verb|ODS_Create| for the parameters, 'nkt', 'nkx'
!     and 'nqcx').  The parameter name must be consistent with
!     the notation of network Common Data form Language (CDL) as
!     described in the Rew et al. (1993).  For this routine, if
!     the parameter is an NetCDF attribute, then the parameter
!     name must consist of the NetCDF variable name followed by
!     the delimiter character, ':', and then the attribute name
!     with no blanks.  If the attribute is global, then the
!     first character in the name must be the delimiter.  If the
!     parameter is a NetCDF dimension, then the name must not
!     contain any delimters.  Some useful parameters are as
!     follows:
! \bv
!        ndays - maximum number of days stored in file
!        nsyn  - number of synoptic times per day
! \ev
!
!    \bigskip {\bf Reference:}
!    \begin{description}
!    \item Rew, Russ, Glenn Davis and Steve Emmerson, 1993:
!        {\em NetCDF User's Guide}, Unidata Program Center,
!        University Corporation for Atmospheric Research,
!        National Science Foundation, Boulder, CO.
!    \end{description}
!
!
! !INTERFACE:
      real function ODS_ParmR ( ParmName, DFault_Val, ierr )
!
! !INPUT PARAMETERS:
      character * (*) ParmName   ! The name of the NetCDF file
                                 !   parameter to be set.  The case
                                 !   of each letter is significant.
      real            DFault_Val ! The default value in case the
                                 !   parameter name is not on the
                                 !   list of specified integers
!
! !OUTPUT PARAMETER:
      integer         ierr       ! The returned error code
!
! !SEE ALSO: 
!     ODS_ParmI    ( Sets the integer value of a NetCDF file parameter )
!     ODS_ParmC    ( Sets the character string of a NetCDF file parameter )
!     ODS_SetParmC ( Saves the character string of a NetCDF file parameter )
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!
! !REVISION HISTORY: 
!     04Sep1996   Redder   Origional version
!     25Oct2000   Redder   Made routine callable.
!
!EOP
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'

*     Functions referened
*     -------------------
      integer       ODS_StrSearch ! Search for the entry among a list of 
                                  !   character strings containing the
                                  !   string being searched

*     Other variables
*     ---------------
      integer       list_entry    ! Index number of the list entry
                                  !   returned from the search string
                                  !   function, ODS_StrSearch.  If no
                                  !   string is found then the value
                                  !   is equal to the parameter,
                                  !   Not_Found, as defined in the
                                  !   header file, ods_hdf.h

*     Search for the NetCDF file parameter name, ParmName, among
*     the list of file parameters that are being set
*     ----------------------------------------------------------
      list_entry = ODS_StrSearch ( ParmName,
     .                             ParmR_nlist,
     .                             ParmR_list,
     .                             ierr )

*     If an error is detected ...
*     ---------------------------
      if ( ierr .ne. NCNoErr ) then

         ODS_ParmR = DFault_Val   ! use default value
         return

      end if

*     If the NetCDF parameter name is not found ...
*     ---------------------------------------------
      if ( list_entry .eq . Not_Found ) then

         ODS_ParmR = DFault_Val    ! use default value

      else

         ODS_ParmR = ParmR_val ( list_entry ) 
                                   ! use the value corresponding to
                                   !   the index number, list_entry
      end if

      return
      end
*..............................................................
