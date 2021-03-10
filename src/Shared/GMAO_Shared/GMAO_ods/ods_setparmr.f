
*..............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:      ODS_SetParmR --- Sets a NetCDF real parameter for creating files
! 
! !DESCRIPTION: 
!    \label{ODS:SetParmR}
!     Sets an ODS parameter for creating files by subsequent
!     calls to the interface routine, \verb|ODS_Create|.  The
!     new value will be in affect until this parameter is
!     changed by this routine (or the routine \verb|ODS_Create|
!     for the parameters, 'nkt', 'nkx' and 'nqcx').  The
!     parameter name must be consistent with the notation of
!     network Common Data form Language (CDL) as described in
!     the Rew et al. (1993).  For this routine, if the parameter
!     is an NetCDF attribute, then the parameter name must
!     consist of the NetCDF variable name followed by the
!     delimiter character, ':', and then the attribute name with
!     no blanks.  If the attribute is global, then the first
!     character in the name must be the delimiter.  If the
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
      subroutine ODS_SetParmR ( ParmName, ParmVal, ierr )
!
! !INPUT PARAMETERS:
      implicit        NONE
      character * (*) ParmName   ! The name of the NetCDF file
                                 !   parameter to be set.  The
                                 !   case of each letter is
                                 !   significant.
      real            ParmVal    ! The value corresponding to the
                                 !   parameter name
!
! !OUTPUT PARAMETER:
      integer         ierr       ! The returned error code
!
!     NOTE: If the parameter name is not valid, then this routine
!           will neither print an error message nor return an
!           error status.  In addition, the routine ODS_Create
!           will ignore the value associated with the invalid
!           name.
!
! !SEE ALSO: 
!     ODS_SetParmI ( Saves the integer of a NetCDF file parameter )
!     ODS_SetParmC ( Saves the character string of a NetCDF file
!                    parameter )
!
! !FILES USED:
!     netcdf.inc, a header file, for defining NetCDF library 
!            parameters
!     ods_hdf.h, a header file, for defining hardwired constants
!            and defining global variables and setting up data
!            structures
!     ods_stdio.h, a header file, for defining standard input/output
!            unit numbers
!
! !REVISION HISTORY:
!     15Apr1996   Redder   Origional version
!     06Oct2000   Redder   Made routine callable.
!
!EOP
!-------------------------------------------------------------------------

      include 'netcdf.inc'
      include 'ods_hdf.h'
      include 'ods_stdio.h'

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

      ierr        = NCNoErr       ! Default error code

*     Search for the NetCDF file parameter name, ParmName, among
*     the list of file parameters that are being set
*     ----------------------------------------------------------
      list_entry = ODS_StrSearch ( ParmName,
     .                             ParmR_nlist,
     .                             ParmR_list,
     .                             ierr )
      if ( ierr .ne. NCNoErr ) return

*     If the NetCDF file parameter name is not found ...
*     --------------------------------------------------
      if ( list_entry .eq. Not_Found ) then

*        Increment the size of the list
*        ------------------------------
         ParmR_nlist = ParmR_nlist + 1

*        Return with an error message if the list size
*        exceeds the allocated space
*        ---------------------------------------------
         if ( ParmR_nlist .gt. Parm_mlist ) then
            ierr = NCSysErr
            write ( stderr, 901 )
            return
         end if

*        Store the NetCDF file parameter name and value
*        ----------------------------------------------
         ParmR_list ( ParmR_nlist ) = ParmName
         ParmR_val  ( ParmR_nlist ) = ParmVal

      else

*        Store the NetCDF file parameter value
*        -------------------------------------
         ParmR_val  ( list_entry  ) = ParmVal

      end if

      return
*     ------

 901  format ( /, ' ODS_SetParmR : The number of parameters to ',
     .         /, '                be set has exceeded the ',
     .         /, '                allocated space. ' )

      end
