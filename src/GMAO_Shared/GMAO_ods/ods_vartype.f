
*..............................................................


      character * ( * ) function ODS_VarType ( VarName )

      implicit NONE

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
! !ROUTINE:      ODS_VarType
!
! !DESCRIPTION: 
!     Returns a character string indicating the FORTRAN type
!     that should be used when using the ODS library to
!     read or write data from "del" files
!
! !INTERFACE: VarType = ODS_VarType  ( VarName )
!
! !INPUT PARAMETER:
      character * (*) VarName ! The variable name 
!
! !REVISION HISTORY: 
!     20Feb96   Redder   Origional version
!     29May97   Redder   Modifications implemented to prevent
!                        character subscript errors
!     17Sep97   Redder   Added code to prevent using the variable,
!                        Var_Name, before it is defined.
!      2Apr98   Redder   Changed variable name "level" to "lev",
!                        Added the variables, qcexcl, qchist, xm
!                        Modified the variable, time
!     16Feb01   Redder   Added error checks for omf, oma and xvec
!
!-------------------------------------------------------------------------

      include 'ods_hdf.h'

      character * ( max_strlen ) Var_Name

      Var_Name = VarName

      if      ( Var_Name ( :3 ) .eq. 'lat'      ) then
         ODS_VarType = 'real'

      else if ( Var_Name ( :3 ) .eq. 'lon'      ) then
         ODS_VarType = 'real'

      else if ( Var_Name ( :3 ) .eq. 'lev'      ) then
         ODS_VarType = 'real'

      else if ( Var_Name ( :5 ) .eq. 'level'    ) then
         ODS_VarType = 'real'

      else if ( Var_Name ( :6 ) .eq. 'julian'   ) then
         ODS_VarType = 'integer'

      else if ( Var_Name ( :4 ) .eq. 'time'     ) then
         ODS_VarType = 'integer'

      else if ( Var_Name ( :2 ) .eq. 'kt'       ) then
         ODS_VarType = 'integer'

      else if ( Var_Name ( :2 ) .eq. 'kx'       ) then
         ODS_VarType = 'integer'

      else if ( Var_Name ( :2 ) .eq. 'ks'       ) then
         ODS_VarType = 'integer'

      else if ( Var_Name ( :2 ) .eq. 'km'       ) then
         ODS_VarType = 'integer'

      else if ( Var_Name ( :2 ) .eq. 'xm'       ) then
         ODS_VarType = 'real'

      else if ( Var_Name ( :3 ) .eq. 'obs'      ) then
         ODS_VarType = 'real'

      else if ( Var_Name ( :3 ) .eq. 'omf'      ) then
         ODS_VarType = 'real'

      else if ( Var_Name ( :3 ) .eq. 'oma'      ) then
         ODS_VarType = 'real'

      else if ( Var_Name ( :4 ) .eq. 'xvec'     ) then
         ODS_VarType = 'real'

      else if ( Var_Name ( :7 ) .eq. 'qc_flag'  ) then
         ODS_VarType = 'integer'

      else if ( Var_Name ( :8 ) .eq. 'mod_flag' ) then
         ODS_VarType = 'integer'

      else if ( Var_Name ( :6 ) .eq. 'qcexcl'   ) then
         ODS_VarType = 'integer'

      else if ( Var_Name ( :6 ) .eq. 'qchist'   ) then
         ODS_VarType = 'integer'

      else
         ODS_VarType = 'indefinite'

      end if

      return
      end
