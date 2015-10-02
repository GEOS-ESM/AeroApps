module mp_netcdf_Mod

  use netcdf
  !use mpi
  implicit none 

  private

  include 'mpif.h'

  public :: mp_readDim
  public :: mp_readvar4D
  public :: mp_readvar3D
  public :: mp_readvar2D
  public :: mp_readvar1D
  public :: mp_readGattr
  public :: mp_readGattr_char  
  public :: mp_readVattr
  public :: mp_colreadvar
  public :: mp_layreadvar
  public :: mp_readvar4Dchunk
  public :: mp_readvar3Dchunk
  public :: mp_readvar2Dchunk
  public :: mp_readvar1Dchunk
  public :: mp_readvarReduced3Dprofile
  public :: mp_readvarReduced2D
  public :: mp_readvarReduced1D
  public :: mp_check
  contains
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readDim
! PURPOSE
!     reads a dimension from a netcdf file
! INPUT
!     dimname   : string of dimension name
!     filename  : file to be read
!     dimvar    : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readDim(dimname, filename, dimvar)
    character(len=*), intent(in)           ::  dimname
    character(len=*), intent(in)           ::  filename
    integer,  intent(inout)                ::  dimvar

    integer                                :: ncid, dimid

    call mp_check( nf90_open(filename,IOR(nf90_nowrite, nf90_mpiio),ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_dimid(ncid, dimname, dimid), "getting dimid for dimension " // dimname)
    call mp_check( nf90_inquire_dimension(ncid, dimid, len = dimvar), "getting value for length of dimension " // dimname)
    call mp_check( nf90_close(ncid), "closing " // filename)
  end subroutine mp_readDim


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvar3D
! PURPOSE
!     reads a 3D variable from a netcdf file all at once
!     can only be called by one processor
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readvar3D(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real, dimension(:,:,:), intent(inout)  ::  var

    integer                                :: ncid, varid


    call mp_check( nf90_open(filename,IOR(nf90_nowrite, nf90_mpiio),ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call mp_check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call mp_check( nf90_close(ncid), "closing " // filename)
  end subroutine mp_readvar3D


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvar4D
! PURPOSE
!     reads a 4D variable from a netcdf file all at once
!     can only be called by one processor
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readvar4D(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real, dimension(:,:,:,:), intent(inout)  ::  var

    integer                                :: ncid, varid


    call mp_check( nf90_open(filename,IOR(nf90_nowrite, nf90_mpiio),ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call mp_check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call mp_check( nf90_close(ncid), "closing " // filename)
  end subroutine mp_readvar4D


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvar2D
! PURPOSE
!     reads a 2D variable from a netcdf file all at once
!     can only be called by one processor
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readvar2D(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:,:), intent(inout)    ::  var

    integer                                :: ncid, varid


    call mp_check( nf90_open(filename,IOR(nf90_nowrite, nf90_mpiio),ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call mp_check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call mp_check( nf90_close(ncid), "closing " // filename)

  end subroutine mp_readvar2D

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvar1D
! PURPOSE
!     reads a 1D variable from a netcdf file all at once
!     can only be called by one processor
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readvar1D(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real*8, dimension(:), intent(inout)      ::  var

    integer                                :: ncid, varid

    call mp_check( nf90_open(filename,IOR(nf90_nowrite, nf90_mpiio),ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call mp_check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call mp_check( nf90_close(ncid), "closing " // filename)
  end subroutine mp_readvar1D

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readGattr
! PURPOSE
!     reads a global attribute from a netcdf file
! INPUT
!     attrname  : string of attribute name
!     filename  : file to be read
!     attrvar   : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     4 May P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readGattr(attrname, filename, attrvar)
    character(len=*), intent(in)           ::  attrname
    character(len=*), intent(in)           ::  filename
    real,  intent(inout)                   ::  attrvar

    integer                                :: ncid

    call mp_check( nf90_open(filename,IOR(nf90_nowrite, nf90_mpiio),ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_get_att(ncid, NF90_GLOBAL, attrname, attrvar), "getting value for global attribute" // attrname)
    call mp_check( nf90_close(ncid), "closing " // filename)
  end subroutine mp_readGattr

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readGattr_char
! PURPOSE
!     reads a global attribute from a netcdf file
! INPUT
!     attrname  : string of attribute name
!     filename  : file to be read
!     attrvar   : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     4 May P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readGattr_char(attrname, filename, attrvar)
    character(len=*), intent(in)           ::  attrname
    character(len=*), intent(in)           ::  filename
    character(len=*), intent(inout)        ::  attrvar

    integer                                :: ncid

    call mp_check( nf90_open(filename,IOR(nf90_nowrite, nf90_mpiio),ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_get_att(ncid, NF90_GLOBAL, attrname, attrvar), "getting value for global attribute" // attrname)
    call mp_check( nf90_close(ncid), "closing " // filename)
  end subroutine mp_readGattr_char  

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readVattr
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
!     13 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readVattr(attrname, filename, varname, attrvar)
    character(len=*), intent(in)           ::  attrname
    character(len=*), intent(in)           ::  filename
    character(len=*), intent(in)           ::  varname
    real,  intent(inout)                   ::  attrvar

    integer                                :: ncid, varid

    call mp_check( nf90_open(filename,IOR(nf90_nowrite, nf90_mpiio),ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    call mp_check( nf90_get_att(ncid, varid, attrname, attrvar), "getting value for variable attribute" // attrname)
    call mp_check( nf90_close(ncid), "closing " // filename)
  end subroutine mp_readVattr  

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_colreadvar
! PURPOSE
!     uses npet processors to read a variable from a netcdf file in chunks of columns
!     is called by multiple processors
!     variable to be read must have dimensions (ew,ns,time=1)
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     7 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_colreadvar(varname, filename, npet, myid, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, intent(in)                       ::  npet, myid
    real, dimension(:,:), intent(inout)       ::  var

    integer                       :: p, startc, countc, endc
    integer                       :: ncid, varid
    integer, dimension(npet)      :: ncol                ! how many columns each processor reads
    integer                       :: im, jm
    integer, dimension(2)         :: dims


    ! Size of domain
    !---------------------
    dims = shape(var)
    im   = dims(1)
    jm   = dims(2)

    ! Everyone Figure out how many columns each PE has to read
    ! -----------------------------
    ncol = 0
    if (npet >= jm) then
      ncol(1:npet) = 1
    else if (npet < jm) then
      ncol(1:npet) = jm/npet
      ncol(npet)   = ncol(npet) + mod(jm,npet)
    end if 


    call mp_check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    do p = 0, npet-1
      if (myid == p) then
        if (p == 0) then
          startc = 1
        else
          startc = sum(ncol(1:p))+1
        end if
        countc = ncol(p+1)
        endc   = startc + countc - 1

        call mp_check( nf90_get_var(ncid, varid,var(:,startc:endc), start = (/ 1, startc, 1 /), count=(/im,countc,1/)), "reading " // varname)

      end if
    end do
    call mp_check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_colreadvar

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_layreadvar
! PURPOSE
!     uses npet processors to read a variable from a netcdf file in chunks of layers
!     is called by multiple processors
!     variable to be read must have dimensions (ew,ns,lev,time=1)
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_layreadvar(varname, filename, npet, myid, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, intent(in)                       ::  npet, myid
    real, dimension(:,:,:), intent(inout)     ::  var

    integer                       :: p, startl, countl, endl
    integer                       :: ncid, varid
    integer, dimension(npet)      :: nlayer                ! how many layers each processor reads
    integer                       :: im, jm, km
    integer, dimension(3)         :: dims

    ! Size of domain
    !---------------------
    dims = shape(var)
    im   = dims(1)
    jm   = dims(2)
    km   = dims(3)

    ! Everyone Figure out how many layers each PE has to read
    ! -----------------------------
    nlayer = 0
    if (npet >= km) then
      nlayer(1:npet) = 1
    else if (npet < km) then
      nlayer(1:npet) = km/npet
      nlayer(npet)   = nlayer(npet) + mod(km,npet)
    end if 


    call mp_check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    do p = 0, npet-1
      if (myid == p) then
        if (p == 0) then
          startl = 1
        else
          startl = sum(nlayer(1:p))+1
        end if
        countl = nlayer(p+1)
        endl   = startl + countl - 1

        call mp_check( nf90_get_var(ncid, varid,var(:,:,startl:endl), start = (/ 1, 1, startl, 1 /), count=(/im,jm,countl,1/)), "reading " // varname)

      end if
    end do
    call mp_check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_layreadvar


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvar4Dchunk
! PURPOSE
!     General code to uses npet processors to read a variable from a netcdf file in chunks across dimenion n
!     is called by multiple processors
!     variable to be read must have 4 dimensions 
! INPUT
!     varname         : string of variable name
!     filename        : file to be read
!     dim = e.g. [im, jm, km, tm]  : size of array to be read
!     n               : dimension to split up
!     var             : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     11 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readvar4Dchunk(varname, filename, dim, n, npet, myid, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, dimension(:), intent(in)         ::  dim
    integer, intent(in)                       ::  n
    integer, intent(in)                       ::  npet, myid
    real, dimension(:,:,:,:), intent(inout)   ::  var

    integer                       :: p, startl, countl, endl
    integer                       :: ncid, varid
    integer, dimension(npet)      :: nlayer                ! how many layers each processor reads
    integer                       :: km
    integer, allocatable          :: countsize(:)

    ! allocate count array
    !---------------------------------------
    allocate(countsize(size(dim)))


    ! Everyone Figure out how many indeces each PE has to read
    ! -----------------------------
    km     = dim(n)
    nlayer = 0
    if (npet >= km) then
      nlayer(1:npet) = 1
    else if (npet < km) then
      nlayer(1:npet) = km/npet
      nlayer(npet)   = nlayer(npet) + mod(km,npet)
    end if 

    call mp_check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    do p = 0, npet-1
      if (myid == p) then
        if (p == 0) then
          startl = 1
        else
          startl = sum(nlayer(1:p))+1
        end if
        countl = nlayer(p+1)
        endl   = startl + countl - 1

        countsize = dim
        countsize(n) = countl

        if (n == 1) then
          call mp_check( nf90_get_var(ncid, varid,var(startl:endl,:,:,:), start = (/ startl, 1, 1, 1 /), count=countsize), "reading " // varname)
        else if (n == 2) then
          call mp_check( nf90_get_var(ncid, varid,var(:,startl:endl,:,:), start = (/ 1, startl, 1, 1 /), count=countsize), "reading " // varname)
        else if (n == 3) then
          call mp_check( nf90_get_var(ncid, varid,var(:,:,startl:endl,:), start = (/ 1, 1, startl, 1 /), count=countsize), "reading " // varname)
        else
          call mp_check( nf90_get_var(ncid, varid,var(:,:,:,startl:endl), start = (/ 1, 1, 1, startl /), count=countsize), "reading " // varname)
        end if                 

      end if
    end do
    call mp_check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_readvar4Dchunk

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvar3Dchunk
! PURPOSE
!     General code to uses npet processors to read a variable from a netcdf file in chunks across dimenion n
!     is called by multiple processors
!     variable to be read must have 3 dimensions 
! INPUT
!     varname         : string of variable name
!     filename        : file to be read
!     dim = e.g. [im, jm, km]  : size of array to be read
!     n               : dimension to split up
!     var             : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     11 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readvar3Dchunk(varname, filename, dim, n, npet, myid, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, dimension(:), intent(in)         ::  dim
    integer, intent(in)                       ::  n  !dimensions to read over in chunks
    integer, intent(in)                       ::  npet, myid
    real, dimension(:,:,:), intent(inout)     ::  var

    integer                       :: p, startl, countl, endl
    integer                       :: ncid, varid
    integer, dimension(npet)      :: nlayer                ! how many layers each processor reads
    integer                       :: km
    integer, allocatable          :: countsize(:)


    ! allocate count array
    !---------------------------------------
    allocate(countsize(size(dim)))


    ! Everyone Figure out how many indeces each PE has to read
    ! -----------------------------
    km     = dim(n)
    nlayer = 0
    if (npet >= km) then
      nlayer(1:npet) = 1
    else if (npet < km) then
      nlayer(1:npet) = km/npet
      nlayer(npet)   = nlayer(npet) + mod(km,npet)
    end if 

    call mp_check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    do p = 0, npet-1
      if (myid == p) then
        if (p == 0) then
          startl = 1
        else
          startl = sum(nlayer(1:p))+1
        end if
        countl = nlayer(p+1)
        endl   = startl + countl - 1

        countsize = dim
        countsize(n) = countl

        if (n == 1) then
          call mp_check( nf90_get_var(ncid, varid,var(startl:endl,:,:), start = (/ startl, 1, 1 /), count=countsize), "reading " // varname)
        else if (n == 2) then
          call mp_check( nf90_get_var(ncid, varid,var(:,startl:endl,:), start = (/ 1, startl, 1 /), count=countsize), "reading " // varname)
        else 
          call mp_check( nf90_get_var(ncid, varid,var(:,:,startl:endl), start = (/ 1, 1, startl /), count=countsize), "reading " // varname)
        end if                 

      end if
    end do
    call mp_check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_readvar3Dchunk  

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvar2Dchunk
! PURPOSE
!     General code to uses npet processors to read a variable from a netcdf file in chunks across dimenion n
!     is called by multiple processors
!     variable to be read must have 3 dimensions 
! INPUT
!     varname         : string of variable name
!     filename        : file to be read
!     dim = e.g. [im, jm]  : size of array to be read
!     n               : dimension to split up
!     var             : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     11 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readvar2Dchunk(varname, filename, dim, n, npet, myid, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, dimension(:), intent(in)         ::  dim
    integer, intent(in)                       ::  n  !dimensions to read over in chunks
    integer, intent(in)                       ::  npet, myid    
    real*8, dimension(:,:), intent(inout)     ::  var

    integer                       :: p, startl, countl, endl
    integer                       :: ncid, varid
    integer, dimension(npet)      :: nlayer                ! how many layers each processor reads
    integer                       :: km
    integer, allocatable          :: countsize(:)


    ! allocate count array
    !---------------------------------------
    allocate(countsize(size(dim)))


    ! Everyone Figure out how many indeces each PE has to read
    ! -----------------------------
    km     = dim(n)
    nlayer = 0
    if (npet >= km) then
      nlayer(1:npet) = 1
    else if (npet < km) then
      nlayer(1:npet) = km/npet
      nlayer(npet)   = nlayer(npet) + mod(km,npet)
    end if 

    call mp_check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    do p = 0, npet-1
      if (myid == p) then
        if (p == 0) then
          startl = 1
        else
          startl = sum(nlayer(1:p))+1
        end if
        countl = nlayer(p+1)
        endl   = startl + countl - 1

        countsize = dim
        countsize(n) = countl

        if (n == 1) then
          call mp_check( nf90_get_var(ncid, varid,var(startl:endl,:), start = (/ startl, 1 /), count=countsize), "reading " // varname)
        else 
          call mp_check( nf90_get_var(ncid, varid,var(:,startl:endl), start = (/ 1, startl /), count=countsize), "reading " // varname)
        end if                 

      end if
    end do
    call mp_check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_readvar2Dchunk  


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvar1Dchunk
! PURPOSE
!     General code to uses npet processors to read a variable from a netcdf file in chunks across dimenion n
!     is called by multiple processors
!     variable to be read must have 3 dimensions 
! INPUT
!     varname         : string of variable name
!     filename        : file to be read
!     dim = e.g. [im]  : size of array to be read
!     n               : dimension to split up
!     var             : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     11 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readvar1Dchunk(varname, filename, dim, n, npet, myid, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, dimension(:), intent(in)         ::  dim
    integer, intent(in)                       ::  n  !dimensions to read over in chunks
    integer, intent(in)                       ::  npet, myid    
    real*8, dimension(:), intent(inout)         ::  var

    integer                       :: p, startl, countl, endl
    integer                       :: ncid, varid
    integer, dimension(npet)      :: nlayer                ! how many layers each processor reads
    integer                       :: km
    integer, allocatable          :: countsize(:)


    ! allocate count array
    !---------------------------------------
    allocate(countsize(size(dim)))


    ! Everyone Figure out how many indeces each PE has to read
    ! -----------------------------
    km     = dim(n)
    nlayer = 0
    if (npet >= km) then
      nlayer(1:npet) = 1
    else if (npet < km) then
      nlayer(1:npet) = km/npet
      nlayer(npet)   = nlayer(npet) + mod(km,npet)
    end if 

    call mp_check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    do p = 0, npet-1
      if (myid == p) then
        if (p == 0) then
          startl = 1
        else
          startl = sum(nlayer(1:p))+1
        end if
        countl = nlayer(p+1)
        endl   = startl + countl - 1

        countsize = dim
        countsize(n) = countl

        call mp_check( nf90_get_var(ncid, varid,var(startl:endl), start = (/ startl /), count=countsize), "reading " // varname)

      end if
    end do
    call mp_check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_readvar1Dchunk    

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvarReduced3Dprofile
! PURPOSE
!     General code to use npet processors to read profiles from a netcdf file  
!     user provides array of indices to read 
! INPUT
!     varname         : string of variable name
!     filename        : file to be read
!     i_work, j_work  : indices to be read
!     km              : size of profile
!     npet            : number of processors
!     myid            : processor id
!     var             : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     28 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readvarReduced3Dprofile(varname, filename, i_work, j_work, km, npet, myid, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, dimension(:), intent(in)         ::  i_work, j_work
    integer, intent(in)                       ::  km
    integer, intent(in)                       ::  npet, myid    
    real, dimension(:,:), intent(inout)       ::  var

    integer                       :: p, c, startl, countl, endl
    integer                       :: ncid, varid
    integer, dimension(npet)      :: nlayer                ! how many layers each processor reads
    integer                       :: clrm
    integer                       :: countsize(3) 

    countsize = (/1,1,km/)

    ! Everyone Figure out how many indeces each PE has to read
    ! -----------------------------
    clrm     = size(i_work)

    nlayer = 0
    if (npet >= clrm) then
      nlayer(1:npet) = 1
    else if (npet < clrm) then
      nlayer(1:npet) = clrm/npet
      nlayer(npet)   = nlayer(npet) + mod(clrm,npet)
    end if 

    call mp_check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    do p = 0, npet-1
      if (myid == p) then
        if (p == 0) then
          startl = 1
        else
          startl = sum(nlayer(1:p))+1
        end if
        countl = nlayer(p+1)
        endl   = startl + countl - 1

        do c = startl, endl
          call mp_check( nf90_get_var(ncid, varid, var(c,:), start = (/i_work(c), j_work(c),1/), count=countsize), "reading " // varname)
        end do
      end if
    end do
    call mp_check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_readvarReduced3Dprofile      


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvarReduced2D
! PURPOSE
!     General code to use npet processors to read a variable from a netcdf file  
!     user provides array of indices to read 
! INPUT
!     varname         : string of variable name
!     filename        : file to be read
!     i_work, j_work  : indices to be read
!     npet            : number of processors
!     myid            : processor id
!     var             : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     28 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readvarReduced2D(varname, filename, i_work, j_work, npet, myid, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, dimension(:), intent(in)         ::  i_work, j_work
    integer, intent(in)                       ::  npet, myid    
    real, dimension(:), intent(inout)         ::  var

    integer                       :: p, c, startl, countl, endl
    integer                       :: ncid, varid
    integer, dimension(npet)      :: nlayer                ! how many layers each processor reads
    integer                       :: km
    integer,parameter             :: countsize(2) = (/1,1/)
    real                          :: temp(1,1) 

    ! Everyone Figure out how many indeces each PE has to read
    ! -----------------------------
    km     = size(i_work)

    nlayer = 0
    if (npet >= km) then
      nlayer(1:npet) = 1
    else if (npet < km) then
      nlayer(1:npet) = km/npet
      nlayer(npet)   = nlayer(npet) + mod(km,npet)
    end if 

    call mp_check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    do p = 0, npet-1
      if (myid == p) then
        if (p == 0) then
          startl = 1
        else
          startl = sum(nlayer(1:p))+1
        end if
        countl = nlayer(p+1)
        endl   = startl + countl - 1

        do c = startl, endl
          call mp_check( nf90_get_var(ncid, varid, temp, start = (/i_work(c), j_work(c)/), count=countsize), "reading " // varname)
          var(c) = temp(1,1)
        end do
      end if
    end do
    call mp_check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_readvarReduced2D      

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvarReduced1D
! PURPOSE
!     General code to use npet processors to read a variable from a netcdf file  
!     user provides array of indices to read 
! INPUT
!     varname         : string of variable name
!     filename        : file to be read
!     i_work          : indices to be read
!     npet            : number of processors
!     myid            : processor id
!     var             : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     28 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_readvarReduced1D(varname, filename, i_work, npet, myid, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, dimension(:), intent(in)         ::  i_work
    integer, intent(in)                       ::  npet, myid    
    real, dimension(:), intent(inout)         ::  var

    integer                       :: p, c, startl, countl, endl
    integer                       :: ncid, varid
    integer, dimension(npet)      :: nlayer                ! how many layers each processor reads
    integer                       :: km
    integer,parameter             :: countsize(1) = (/1/)
    real                          :: temp(1) 

    ! Everyone Figure out how many indeces each PE has to read
    ! -----------------------------
    km     = size(i_work)

    nlayer = 0
    if (npet >= km) then
      nlayer(1:npet) = 1
    else if (npet < km) then
      nlayer(1:npet) = km/npet
      nlayer(npet)   = nlayer(npet) + mod(km,npet)
    end if 

    call mp_check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call mp_check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    do p = 0, npet-1
      if (myid == p) then
        if (p == 0) then
          startl = 1
        else
          startl = sum(nlayer(1:p))+1
        end if
        countl = nlayer(p+1)
        endl   = startl + countl - 1

        do c = startl, endl
          call mp_check( nf90_get_var(ncid, varid, temp, start = (/i_work(c)/), count=countsize), "reading " // varname)
          var(c) = temp(1)
        end do
      end if
    end do
    call mp_check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_readvarReduced1D      



!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     mp_check
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

  subroutine mp_check(status, loc)

    integer, intent(in) :: status
    character(len=*), intent(in) :: loc 

    if(status /= NF90_NOERR) then
      write (*,*) "Error at ", loc
      write (*,*) NF90_STRERROR(status)
    end if

  end subroutine mp_check

end module mp_netcdf_Mod  