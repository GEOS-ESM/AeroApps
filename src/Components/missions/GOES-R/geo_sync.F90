!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    geo_sync
! PURPOSE
!     synchronize one geostationary swath with another
!     Heuristic for GOES-R and TEMPO
! INPUT
!     inst_a  : instrument name for which you have scan times
!     inst_b  : instrument name that needs to be synchronized
!     indir : main directory for input data
!     outdir: directory for output data
!     nccs  : flag to indicate nccs directory convention
! OUTPUT
!     None
!  HISTORY
!     Aug 2015 P. Castellanos 
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#  include "MAPL_Generic.h"
#  include "MAPL_ErrLogMain.h"
program geo_sync

  use ESMF                         ! ESMF modules
  use MAPL_Mod
  use netcdf                       ! for reading the NR files
  use mp_netcdf_Mod
  use netcdf_Mod

  implicit none
  include "mpif.h"

! ESMF Objects
!----------------
  type(ESMF_Config)       :: cf
  type(ESMF_VM)           :: vm 

! Information to be retrieved from resource file 
! ------------------------------------------------
  character(len=256)                    :: arg                    ! command line rc file argument
  character(len=256)                    :: inst_a, inst_b, indir, outdir
  logical                               :: nccs

! File names
! ----------
  character(len=256)                    :: INVa_file, INVb_file, TIME_file, OUT_file

! Main data arrays
! ---------------------------------------------
  real*8, allocatable                   :: CLONa(:,:), CLATa(:,:)
  real*8, allocatable                   :: CLONb(:,:), CLATb(:,:)
  real*8, allocatable                   :: TIME(:) 
  real*8, allocatable                   :: SYNCTIME(:,:) 
  real*8, allocatable                   :: dsq(:,:)


! Satellite domain variables
!------------------------------
  integer                               :: im, jm                                    ! size of TEMPO domain
  integer                               :: xm, ym                                    ! size of GOES-R domain
  integer                               :: i, j, x, y                                ! domain working variables
  integer                               :: counti                                    ! array indices and counts for each processor
  integer, allocatable                  :: nclr(:)                                   ! how many pixels each processor works on
  integer                               :: c                                         ! n-pixel working variable
  integer                               :: left

! netcdf variables
!----------------------  
  integer                               :: ncid                                      ! netcdf file id
  integer                               :: TimeID 

! Miscellaneous
! -------------
  integer                               :: ierr, rc, status                            ! MPI error message
  integer                               :: status_mpi(MPI_STATUS_SIZE)                 ! MPI status
  integer                               :: myid, npet, CoresPerNode                    ! MPI dimensions and processor id
  integer                               :: p                                           ! i-processor
  character(len=100)                    :: msg, msgr                                   ! message to be printed
  integer                               :: progress                                    ! 
  real*8                                :: MISSING 
  real*8                                :: lon, lat
  integer, dimension(2)                 :: loc

! System tracking variables
! -----------------------------
  integer*8                             :: t1, t2, clock_max
  real*8                                :: clock_rate
  character(len=*), parameter           :: Iam = 'geo_vlidort'

!                               END OF VARIABLE DECLARATIONS
!----------------------------------------------------------------------------------------------------------
  
! Start Timing
! ------------
  call system_clock ( t1, clock_rate, clock_max )
  progress = -1

! Initialize MPI with ESMF
! ------------------------
  call ESMF_Initialize (logkindflag=ESMF_LOGKIND_NONE, vm=vm, __RC__)

  call ESMF_VMGet(vm, localPET=myid, PETcount=npet) 
  if ( MAPL_am_I_root() ) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'


! Parse Resource file provided at command line for input info 
! -----------------------------------------------------------
  call getarg(1, arg)
  call get_config(arg)

! Write out settings to use
! --------------------------
  if (myid == 0) then
    write(*,*) 'TEMPO lg1 file ', trim(INVa_file)
    write(*,*) 'GOES-R lg1 file ', trim(INVb_file)
    write(*,*) 'Scantimes ', trim(TIME_file)
    write(*,*) 'OUTFILE ', trim(OUT_file)
  end if 

! Query for domain dimensions and missing value
!----------------------------------------------
  call mp_readDim("ew", INVa_file, im)
  call mp_readDim("ns", INVa_file, jm)
  call mp_readDim("ew", INVb_file, xm)
  call mp_readDim("ns", INVb_file, ym)  

! Split up domain among processors
!----------------------------------------------
  allocate (nclr(npet)) 
  nclr = 0
  if (npet >= xm) then
    nclr(1:xm) = 1
  else 
    nclr = xm/npet
    left = mod(xm,npet)
    if (npet >= left) then
      nclr(1:left) = nclr(1:left) + 1
    else 
      do while (left > 0) 
        if (npet >= left) then
          nclr(1:left) = nclr(1:left) + 1
          left = 0
        else
          nclr = nclr + 1
          left = left - npet
        end if
      end do
    end if
  end if   

! Allocate arrays that will be copied on each processor - unshared
! -----------------------------------------------------------------
  call allocate_unshared()

! Create OUTFILE
! --------------
  call create_outfile(TimeID)

! Read in data
! ---------------
  call mp_colreader("clon", INVb_file, nclr, CLONb)
  call mp_colreader("clat", INVb_file, nclr, CLATb)
  call mp_readvar2D("clon", INVa_file, CLONa)
  call mp_readvar2D("clat", INVa_file, CLATa)
  call mp_readvar1D("scanTime", TIME_file, TIME)
! Initialize outputs to be safe
! -------------------------------
  MISSING          = 1e15
  SYNCTIME         = MISSING
  counti = nclr(myid+1)

  do x = 1, counti
    do y=1, ym
      lon = CLONb(x,y)
      lat = CLATb(x,y)
      if ((.not. IS_MISSING(lon,MISSING)) .and. (.not. IS_MISSING(lat,MISSING))) then
        dsq = ((lon - CLONa)**2 + (lat - CLATa)**2)
        loc = minloc(dsq)
        SYNCTIME(x,y) = TIME(loc(1))
      end if
    end do
  end do   

!                             Write to main OUT_File
!                             ----------------------
  call check( nf90_open(OUT_file, nf90_write, ncid), "opening file " // OUT_file )

  call check(nf90_put_var(ncid, TimeID, SYNCTIME), "writing out Time")
  call check( nf90_close(ncid), "close outfile" )


  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)
! -----------------------------------------------------------------------------------

  contains


  logical function IS_MISSING(x,MISSING)
    real*8, intent(in)     :: x
    real*8, intent(in)     :: MISSING

    IS_MISSING = abs(x/MISSING-1)<0.001
    return
  end function IS_MISSING

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_colreader
! PURPOSE
!     uses npet processors to read a variable from a netcdf file in chunks of columns
!     is called by multiple processors
!     variable to be read must have dimensions (ew,ns)
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     7 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine mp_colreader(varname, filename, nclr, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, dimension(:), intent(in)         ::  nclr
    real*8, dimension(:,:), intent(inout)       ::  var

    integer                       :: starti, counti, endi
    integer                       :: ncid, varid
    integer                       :: im, jm
    integer, dimension(2)         :: dims, dimIDS
    integer                       :: num_1, num_2


    ! Size of domain
    !---------------------
    dims = shape(var)
    im   = dims(1)
    jm   = dims(2)


    call check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)

    call check( nf90_inquire_variable(ncid, varid, dimids = dimIDs), 'get dimids')
    call check( nf90_inquire_dimension(ncid, dimIDs(1), len = num_1), 'get num_1')
    call check( nf90_inquire_dimension(ncid, dimIDs(2), len = num_2), 'get num_2')

    if (myid == 0) then
      starti = 1
    else
      starti = sum(nclr(1:myid))+1
    end if
    counti = nclr(myid+1)
    endi   = starti + counti - 1
    call check( nf90_get_var(ncid, varid,var, start = (/ starti, 1 /), count=(/counti,jm/)), "reading " // varname)
    call check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_colreader

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     create_outfile
! PURPOSE
!     create netcdf outfile
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     Aug 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine create_outfile(TimeID)
    character(len=256)                  :: filename
    character(len=2000)                 :: comment
    integer                             :: starti, endi, counti
    integer                             :: ncid
    integer                             :: ewDimID, nsDimID
    integer,intent(out)                 :: TimeID                          

    if (myid == 0) then
      starti = 1
    else
      starti = sum(nclr(1:myid))+1
    end if
    counti = nclr(myid+1)
    endi   = starti + counti - 1
    ! Open File
    write(filename,'(A,I4.4,A,I4.4,A)') trim(OUT_file),starti,'-',endi,'.nc4'
    OUT_file = trim(filename)
    call check(nf90_create(OUT_file, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // filename)

    ! Create dimensions
    call check(nf90_def_dim(ncid, "ew", counti, ewDimID), "creating ew dimension")
    call check(nf90_def_dim(ncid, "ns", ym, nsDimID), "creating ns dimension")

    ! Global Attributes
    write(comment,'(A)') 'GOES-R intermediate synchronized scan time file'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'title',trim(comment)),"title attr")

    ! Define Variables
!                                     Data
!                                     -----    
    call check(nf90_def_var(ncid,'scanTime',nf90_float,(/ewDimID,nsDimID/),TimeID),"create scanTime var")    

    ! Variable Attributes
!                                     Time
!                                     ----      
    write(comment,'(A)') 'Time of Scan'
    call check(nf90_put_att(ncid,TimeID,'long_name',trim(adjustl(comment))),"long_name attr")
    call check(nf90_put_att(ncid,TimeID,'units','seconds since 2000-01-01 00:00:00'),"units attr")

    !Leave define mode
    call check(nf90_enddef(ncid),"leaving define mode")

    call check( nf90_close(ncid), "close outfile" )
  end subroutine create_outfile


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     filenames
! PURPOSE
!     populate file name varaibles
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     26 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
subroutine filenames()

  if (nccs) then
    ! INFILES
    write(TIME_file,'(6A)') trim(indir),'/',lower_to_upper(trim(inst_a)),'/DATA/LevelB/Y2005/M12/D31/',trim(inst_a),'-g5nr.lb2.met_Nv.20051231_01z.nc4'
    write(INVa_file,'(6A)') trim(indir),'/',lower_to_upper(trim(inst_a)),'/DATA/LevelG/invariant/',trim(inst_a),'.lg1.invariant.nc4'
    write(INVb_file,'(6A)') trim(indir),'/',lower_to_upper(trim(inst_b)),'/DATA/LevelG/invariant/',trim(inst_b),'.lg1.invariant.nc4'
  else
    write(TIME_file,'(6A)') trim(indir),'/',lower_to_upper(trim(inst_a)),'/LevelB/Y2005/M12/D31/',trim(inst_a),'-g5nr.lb2.met_Nv.20051231_01z.nc4'
    write(INVa_file,'(6A)') trim(indir),'/',lower_to_upper(trim(inst_a)),'/LevelG/invariant/',trim(inst_a),'.lg1.invariant.nc4'
    write(INVb_file,'(6A)') trim(indir),'/',lower_to_upper(trim(inst_b)),'/LevelG/invariant/',trim(inst_b),'.lg1.invariant.nc4'  
  end if  
! OUTFILES
  write(OUT_file,'(4A)') trim(outdir),'/',trim(inst_b),'.lg1.scantime.'

end subroutine filenames


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     lower_to_upper
! PURPOSE
!     converts lower case characters to uppercase
! INPUT
!     word
! OUTPUT
!     lower_to_upper
!  HISTORY
!     18 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   

  function lower_to_upper ( word )
    implicit none

    character(len=*), intent(in)   :: word
    character(len=:),allocatable   :: lower_to_upper
    integer                        :: i,n

    character(*), parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
    character(*), parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    allocate(character(len=len(word)) :: lower_to_upper)

    lower_to_upper = word

    ! Loop over string elements
    do i = 1, LEN(word)
    ! Find location of letter in lower case constant string
      n = INDEX(LOWER_CASE, word( i:i ))
    ! If current substring is a lower case letter, make it upper case
      if ( n /= 0 ) lower_to_upper( i:i ) = UPPER_CASE( n:n )
    end do

  end function lower_to_upper


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     reducedProfile
! PURPOSE
!     reduces a 3D array of profiles into a 2D array of profiles according to a mask
! INPUT
!     var: variables to reduce
!     mask: 2-D mask
! OUTPUT
!     none
!  HISTORY
!     29 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine reduceProfile(var,mask,reducedProfile)
    real,intent(in),dimension(:,:,:)               :: var
    logical,intent(in),dimension(:,:)              :: mask
    real,intent(inout),dimension(:,:)              :: reducedProfile

    integer                                        :: im, jm, km, k

    im = size(var,1)
    jm = size(var,2)
    km = size(var,3)

    do k = 1, km
      reducedProfile(:,k) = pack(reshape(var(:,:,k),(/im,jm/)),mask)
    end do

  end subroutine reduceProfile




!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     allocate_unshared
! PURPOSE
!     allocates all the arrays that I need, but each processor has its own copy
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine allocate_unshared()
  ! Needed for vlidort
  ! -----------------------
    allocate (CLONa(im,jm))
    allocate (CLATa(im,jm))
    allocate (CLONb(nclr(myid+1),ym))
    allocate (CLATb(nclr(myid+1),ym))

    allocate (TIME(im))
    allocate (SYNCTIME(nclr(myid+1),ym))

    allocate (dsq(im,jm))
  end subroutine allocate_unshared

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    get_config()
! PURPOSE
!     Read and parse resource file
!     Set filenames
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     29 May P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine get_config(rcfile)
    character(len=*),intent(in)      :: rcfile

    cf = ESMF_ConfigCreate()
    call ESMF_ConfigLoadFile(cf, fileName=trim(rcfile), __RC__)

    ! Read in variables
    ! -----------------------
    call ESMF_ConfigGetAttribute(cf, inst_a, label = 'INSTNAMEA:',__RC__)
    call ESMF_ConfigGetAttribute(cf, inst_b, label = 'INSTNAMEB:',__RC__)
    call ESMF_ConfigGetAttribute(cf, indir, label = 'INDIR:',__RC__)
    call ESMF_ConfigGetAttribute(cf, outdir, label = 'OUTDIR:',default=indir)
    call ESMF_ConfigGetAttribute(cf, nccs, label = 'NCCS:',default=.true.)

    ! INFILES & OUTFILE set names
    ! -------------------------------
    call filenames()

    call ESMF_ConfigDestroy(cf)

  end subroutine get_config




end program geo_sync
