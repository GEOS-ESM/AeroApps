module netcdf_Mod

  use netcdf
  !use mpi
  implicit none 

  private

  public :: readvar4D
  public :: readvar3D
  public :: readvar2D
  public :: readvar1D
  public :: readGattr
  public :: readGattr_char  
  public :: readVattr
  public :: readDim
  public :: check
  contains
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readDim
! PURPOSE
!     reads a dimension from a netcdf file
! INPUT
!     dimname   : string of dimension name
!     filename  : file to be read
!     dimvar    : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine readDim(dimname, filename, dimvar)
    character(len=*), intent(in)           ::  dimname
    character(len=*), intent(in)           ::  filename
    integer,  intent(inout)                ::  dimvar

    integer                                :: ncid, dimid

    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_dimid(ncid, dimname, dimid), "getting dimid for dimension " // dimname)
    call check( nf90_inquire_dimension(ncid, dimid, len = dimvar), "getting value for length of dimension " // dimname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readDim


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar3D
! PURPOSE
!     reads a 3D variable from a netcdf file all at once
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine readvar3D(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real, dimension(:,:,:), intent(inout)  ::  var

    integer                                :: ncid, varid


    call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar3D


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar4D
! PURPOSE
!     reads a 4D variable from a netcdf file all at once
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine readvar4D(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real, dimension(:,:,:,:), intent(inout)  ::  var

    integer                                :: ncid, varid


    call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar4D


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar2D
! PURPOSE
!     reads a 2D variable from a netcdf file all at once
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine readvar2D(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:,:), intent(inout)    ::  var

    integer                                :: ncid, varid


    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)

  end subroutine readvar2D

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar1D
! PURPOSE
!     reads a 1D variable from a netcdf file all at once
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine readvar1D(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:), intent(inout)      ::  var

    integer                                :: ncid, varid

    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar1D



!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readGattr
! PURPOSE
!     reads a global attribute from a netcdf file
! INPUT
!     attrname  : string of attribute name
!     filename  : file to be read
!     attrvar   : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine readGattr(attrname, filename, attrvar)
    character(len=*), intent(in)           ::  attrname
    character(len=*), intent(in)           ::  filename
    real,  intent(inout)                   ::  attrvar

    integer                                :: ncid

    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_get_att(ncid, NF90_GLOBAL, attrname, attrvar), "getting value for global attribute" // attrname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readGattr

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readGattr_char
! PURPOSE
!     reads a global attribute from a netcdf file
! INPUT
!     attrname  : string of attribute name
!     filename  : file to be read
!     attrvar   : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine readGattr_char(attrname, filename, attrvar)
    character(len=*), intent(in)           ::  attrname
    character(len=*), intent(in)           ::  filename
    character(len=*), intent(inout)        ::  attrvar

    integer                                :: ncid

    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_get_att(ncid, NF90_GLOBAL, attrname, attrvar), "getting value for global attribute" // attrname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readGattr_char  

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readVattr
! PURPOSE
!     reads a variable attribute from a netcdf file
! INPUT
!     attrname  : string of attribute name
!     filename  : file to be read
!     varname   : variable to be read
!     attrvar   : the output variable to be read to
! OUTPUT
!     None
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine readVattr(attrname, filename, varname, attrvar)
    character(len=*), intent(in)           ::  attrname
    character(len=*), intent(in)           ::  filename
    character(len=*), intent(in)           ::  varname
    real,  intent(inout)                   ::  attrvar

    integer                                :: ncid, varid

    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    call check( nf90_get_att(ncid, varid, attrname, attrvar), "getting value for variable attribute" // attrname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readVattr  

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     check
! PURPOSE
!     tests the return value of an NF90 call
!     prints a message (loc) if the return value indicates an error
! INPUT
!     status : NF90 return value to be checked
!     loc    : use character string indicating where in the code the NF90 call is
! OUTPUT
!     Writes to the standard output the loc and the NF90 error
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine check(status, loc)

    integer, intent(in) :: status
    character(len=*), intent(in) :: loc 

    if(status /= NF90_NOERR) then
      write (*,*) "Error at ", loc
      write (*,*) NF90_STRERROR(status)
      stop 2
    end if

  end subroutine check

end module netcdf_Mod  