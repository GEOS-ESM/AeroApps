module netcdf_Mod

  use netcdf
  !use mpi
  implicit none 

  private

  public :: readvar5D
  public :: readvar4D
  public :: readvar3D
  public :: readvar2D
  public :: readvar2Dgrp  
  public :: readvar1Dgrp    
  public :: readvar1D
  public :: readGattr
  public :: readGattr_char  
  public :: readVattr
  public :: readVattrgrp
  public :: readDim
  public :: check
  public :: readvar3Dslice

  interface readvar2Dgrp
       module procedure readvar2DgrpR4
       module procedure readvar2DgrpR8
  end interface

  interface readvar2D
       module procedure readvar2DR4
       module procedure readvar2DR8
  end interface  

  interface readvar3D
       module procedure readvar3DR4
       module procedure readvar3DR8
  end interface  

  interface readvar4D
       module procedure readvar4DR4
       module procedure readvar4DR8
  end interface  

  interface readvar5D
       module procedure readvar5DR4
       module procedure readvar5DR8
  end interface


  interface readvar1Dgrp
       module procedure readvar1DgrpR4
       module procedure readvar1DgrpR8
  end interface


  interface readvar1D
       module procedure readvar1DR4
       module procedure readvar1DR8
  end interface
  
  contains

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar3Dslice
! PURPOSE
!     reads a slice of a 3D variable from a netcdf file all at once
!     can only be called by one processor
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     2018 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine readvar3Dslice(varname, filename, countsize, n, startl, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    integer, dimension(:), intent(in)      ::  countsize
    integer, intent(in)                    ::  n
    integer, intent(in)                    ::  startl
    real, dimension(:,:,:), intent(inout)  ::  var


    integer                                :: ncid, varid


    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    if (n == 1) then
      call check( nf90_get_var(ncid, varid,var, start = (/ startl, 1, 1 /), count=countsize), "reading " // varname)
    else if (n == 2) then
      call check( nf90_get_var(ncid, varid,var, start = (/ 1, startl, 1 /), count=countsize), "reading " // varname)
    else 
      call check( nf90_get_var(ncid, varid,var, start = (/ 1, 1, startl /), count=countsize), "reading " // varname)
    end if                 

    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar3Dslice


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
!    readvar3DR4
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
  subroutine readvar3DR4(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real, dimension(:,:,:), intent(inout)  ::  var

    integer                                :: ncid, varid


    call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar3DR4

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar3DR8
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
  subroutine readvar3DR8(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:,:,:), intent(inout)  ::  var

    integer                                :: ncid, varid


    call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar3DR8

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar4DR4
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
  subroutine readvar4DR4(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real, dimension(:,:,:,:), intent(inout)  ::  var

    integer                                :: ncid, varid


    call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar4DR4

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar4DR8
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
  subroutine readvar4DR8(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:,:,:,:), intent(inout)  ::  var

    integer                                :: ncid, varid


    call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar4DR8  


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar5DR4
! PURPOSE
!     reads a 5D variable from a netcdf file all at once
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine readvar5DR4(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real, dimension(:,:,:,:,:), intent(inout)  ::  var

    integer                                :: ncid, varid


    call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar5DR4

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar5DR8
! PURPOSE
!     reads a 5D variable from a netcdf file all at once
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine readvar5DR8(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:,:,:,:,:), intent(inout)  ::  var

    integer                                :: ncid, varid


    call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar5DR8  

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar2DR4
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
  subroutine readvar2DR4(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real, dimension(:,:), intent(inout)    ::  var

    integer                                :: ncid, varid


    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)

  end subroutine readvar2DR4

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar2DR8
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
  subroutine readvar2DR8(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:,:), intent(inout)    ::  var

    integer                                :: ncid, varid


    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)

  end subroutine readvar2DR8

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar2DgrpR4
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
  subroutine readvar2DgrpR4(varname, groupname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  groupname
    character(len=*), intent(in)           ::  filename
    real, dimension(:,:), intent(inout)    ::  var

    integer                                :: ncid, grp_ncid, varid


    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_ncid(ncid,groupname,grp_ncid), "getting ncid for " // groupname)    
    call check( nf90_inq_varid(grp_ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(grp_ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)

  end subroutine readvar2DgrpR4

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar2DgrpR8
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
  subroutine readvar2DgrpR8(varname, groupname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  groupname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:,:), intent(inout)    ::  var

    integer                                :: ncid, grp_ncid, varid


    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_ncid(ncid,groupname,grp_ncid), "getting ncid for " // groupname)    
    call check( nf90_inq_varid(grp_ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(grp_ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)

  end subroutine readvar2DgrpR8

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar1DgrpR4
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
  subroutine readvar1DgrpR4(varname, groupname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  groupname
    character(len=*), intent(in)           ::  filename
    real, dimension(:), intent(inout)      ::  var

    integer                                :: ncid, grp_ncid, varid


    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_ncid(ncid,groupname,grp_ncid), "getting ncid for " // groupname)    
    call check( nf90_inq_varid(grp_ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(grp_ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)

  end subroutine readvar1DgrpR4

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar1DgrpR8
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
  subroutine readvar1DgrpR8(varname, groupname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  groupname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:), intent(inout)    ::  var

    integer                                :: ncid, grp_ncid, varid


    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_ncid(ncid,groupname,grp_ncid), "getting ncid for " // groupname)    
    call check( nf90_inq_varid(grp_ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(grp_ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)

  end subroutine readvar1DgrpR8

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar1DR8
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
  subroutine readvar1DR8(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:), intent(inout)      ::  var

    integer                                :: ncid, varid

    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar1DR8

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    readvar1DR4
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
  subroutine readvar1DR4(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real, dimension(:), intent(inout)      ::  var

    integer                                :: ncid, varid

    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar1DR4


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
!    readVattrgrp
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
  subroutine readVattrgrp(attrname, groupname, filename, varname, attrvar)
    character(len=*), intent(in)           ::  attrname
    character(len=*), intent(in)           ::  groupname    
    character(len=*), intent(in)           ::  filename
    character(len=*), intent(in)           ::  varname
    real,  intent(inout)                   ::  attrvar

    integer                                :: ncid, grp_ncid, varid

    call check( nf90_open(filename,nf90_nowrite,ncid), "opening file " // filename)
    call check( nf90_inq_ncid(ncid,groupname,grp_ncid), "getting ncid for " // groupname)   
    call check( nf90_inq_varid(grp_ncid, varname, varid), "getting varid for " // varname)
    call check( nf90_get_att(grp_ncid, varid, attrname, attrvar), "getting value for variable attribute" // attrname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readVattrgrp  


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
