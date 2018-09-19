!
!  Simple fortran wrapper for parallel netcdf read
!  
!
!  Patricia Castellanos
!  May 2017
!.............................................................................

subroutine py_mp_readvar1D(varname, filename, var)

    use mp_netcdf_Mod, only: mp_readvar1d  
    implicit None

  ! !INPUT PARAMETERS:

    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:), intent(inout)      ::  var


    call mp_readvar1d (varname, filename, var )  


end subroutine py_mp_readvar1D

subroutine py_mp_readvar2D(varname, filename, var)

    use mp_netcdf_Mod, only: mp_readvar2d  
    implicit None

  ! !INPUT PARAMETERS:

    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:,:), intent(inout)    ::  var


    call mp_readvar2d (varname, filename, var )  


end subroutine py_mp_readvar2D

subroutine py_mp_readvar3D(varname, filename, var)

    use mp_netcdf_Mod, only: mp_readvar3d  
    implicit None

  ! !INPUT PARAMETERS:

    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real, dimension(:,:,:), intent(inout)  ::  var


    call mp_readvar3d (varname, filename, var )  


end subroutine py_mp_readvar3D

subroutine py_mp_readvar4D(varname, filename, var)

    use mp_netcdf_Mod, only: mp_readvar4d  
    implicit None

  ! !INPUT PARAMETERS:

    character(len=*), intent(in)             ::  varname
    character(len=*), intent(in)             ::  filename
    real, dimension(:,:,:,:), intent(inout)  ::  var


    call mp_readvar4d (varname, filename, var )  


end subroutine py_mp_readvar4D

