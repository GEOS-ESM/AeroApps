!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    cloud_calculator
! PURPOSE
!     Writes out cloud properties
! INPUT
!     rcFile: contains all the input parameters
! OUTPUT
!     None
!  HISTORY
!     Feb 2019 P. Castellanos adapted from leo_vlidort_cloud.F90
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#  include "MAPL_Generic.h"
#  include "MAPL_ErrLogMain.h"
program cloud_calculator

  use ESMF                         ! ESMF modules
  use ESMF_CFIOMod, only: ESMF_CFIODownBit
  use MAPL_Mod
  use MAPL_ShmemMod                ! The SHMEM infrastructure
  use netcdf                       ! for reading the NR files
  use mp_netcdf_Mod
  use netcdf_Mod
  use cloud_MieMod

  implicit none
  include "mpif.h"

! ESMF Objects
!----------------
  type(ESMF_Config)       :: cf
  type(ESMF_VM)           :: vm 

  character(len=256)                    :: arg                    ! command line rc file argument
  integer,parameter                     :: nMom = 300, nPol =6, nband = 3, nobs = 1
  integer                               :: nch 
  character(len=4)                      :: bandname
  integer                               :: clrm
  integer, allocatable                  :: iIndex(:),jIndex(:)                       ! indices of pixels to run  
  integer, allocatable                  :: nclr(:)                                   ! how many clear pixels each processor works on  
! Information to be retrieved from resource file 
! ------------------------------------------------
  character(len=256)                    :: IcldTable, LcldTable
  character(len=256)                    :: ANG_file, CLD_file, OUT_file

! Global, 3D inputs to be allocated using SHMEM
! ---------------------------------------------
  real, pointer                         :: REI(:,:,:) => null()  
  real, pointer                         :: REL(:,:,:) => null()  
  real, pointer                         :: TAUI(:,:,:) => null()  
  real, pointer                         :: TAUL(:,:,:) => null() 
  integer, pointer                      :: indices(:) => null()


!                                  Final Shared Arrays
!                                  -------------------
  real, pointer                         :: LTAU(:,:,:) => null()                 ! liquid cloud optical depth
  real, pointer                         :: LSSA(:,:,:) => null()                 ! liquid cloud single scattering albedo
  real, pointer                         :: LG(:,:,:) => null()                   ! liquid cloud asymmetry factor
  real, pointer                         :: ITAU(:,:,:) => null()                 ! ice cloud optical depth
  real, pointer                         :: ISSA(:,:,:) => null()                 ! ice cloud single scattering albedo
  real, pointer                         :: IG(:,:,:) => null()                   ! ice cloud asymmetry factor

! working variables
!------------------------------
  real, allocatable                     :: channels(:)                              ! channels to simulate
  integer                               :: ch                                       ! i-channel  
  integer                               :: iband                                    ! i-band
  real, allocatable                     :: VtauLcl(:,:,:)          ! cloud aerosol optical depth
  real, allocatable                     :: VssaIcl(:,:,:)          ! cloud aerosol optical depth
  real, allocatable                     :: VssaLcl(:,:,:)          ! cloud single scattering albedo
  real, allocatable                     :: VtauIcl(:,:,:)          ! cloud single scattering albedo
  real, allocatable                     :: VgLcl(:,:,:)            ! cloud asymmetry factor
  real, allocatable                     :: VgIcl(:,:,:)            ! cloud asymmetry factor
  real,allocatable                      :: VpmomIcl(:,:,:,:,:)                      ! elements of scattering phase matrix for vector calculations
  real,allocatable                      :: VpmomLcl(:,:,:,:,:)                      ! elements of scattering phase matrix for vector calculations

! Satellite domain variables
!------------------------------
  integer                               :: im, jm, km, tm                            ! size of satellite domain
  integer                               :: i, j, k, n                                ! satellite domain working variable
  integer                               :: starti, counti, endi                      ! array indices and counts for each processor
  integer                               :: c, cc

! Miscellaneous
! -------------
  integer                               :: ierr, rc, status                            ! MPI error message
  integer                               :: status_mpi(MPI_STATUS_SIZE)                 ! MPI status
  integer                               :: myid, npet, CoresPerNode                    ! MPI dimensions and processor id
  integer                               :: p                                           ! i-processor
  character(len=100)                    :: msg                                         ! message to be printed
  real                                  :: progress                                    ! 
  real                                  :: g5nr_missing  

! System tracking variables
! -----------------------------
  character(len=*), parameter           :: Iam = 'cloud_calculator'

!                               END OF VARIABLE DECLARATIONS
!----------------------------------------------------------------------------------------------------------
  
! Initialize MPI with ESMF
! ------------------------
  call ESMF_Initialize (logkindflag=ESMF_LOGKIND_NONE, vm=vm, __RC__)

  call ESMF_VMGet(vm, localPET=myid, PETcount=npet) 
  if ( MAPL_am_I_root() ) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'

! Initialize SHMEM
! ----------------
  CoresPerNode = MAPL_CoresPerNodeGet(MPI_COMM_WORLD,rc=ierr) ! a must
  call MAPL_InitializeShmem(rc=ierr)

! Parse Resource file provided at command line for input info 
! -----------------------------------------------------------
  call getarg(1, arg)
  call get_config(arg)

  if (MAPL_am_I_root()) then
    write(*,*) '<> Reading Data from'
    write(*,*) 'ANG_file: ',trim(ANG_file)
    write(*,*) 'CLD_file: ',trim(CLD_file)
  end if

! Query for domain dimensions and missing value
!----------------------------------------------
  call mp_readDim("ccd_pixels", CLD_file, im)
  call mp_readDim("number_of_scans", CLD_file, jm)
  call mp_readDim("lev", CLD_file, km)
  tm = 1
  
  call mp_readVattr("missing_value", CLD_FILE, "REI", g5nr_missing)
  
! Allocate arrays that will be shared
! -----------------------------------------------------------------
  call allocate_shared()
  call MAPL_SyncSharedMemory(rc=ierr)

! Read in the global arrays
! ------------------------------
  call read_cld_Tau() 
  call MAPL_SyncSharedMemory(rc=ierr)

! Create out file 
! ----------------------
  if ( MAPL_am_I_root() ) call create_cldfile()
  call MAPL_SyncSharedMemory(rc=ierr)

! Collapse indeces
!------------------------------------------
  clrm = im*jm
  allocate(iIndex(clrm))
  allocate(jIndex(clrm))
  k = 1
  do i=1,im
    do j=1,jm
        iIndex(k) = i
        jIndex(k) = j
        k = k + 1      
    end do
  end do

! Split up filtered domain among processors
!----------------------------------------------
  allocate (nclr(npet))
  nclr = 0
  if (npet >= clrm) then
    nclr(1:clrm) = 1
  else if (npet < clrm) then
    nclr(1:npet) = clrm/npet
    nclr(npet)   = nclr(npet) + mod(clrm,npet)
  end if 

  if (myid == 0) then
    starti  = 1
  else
    starti = sum(nclr(1:myid)) + 1
  end if    
  counti = nclr(myid+1)
  endi   = starti + counti - 1 

! Shuffle indices 1-clrm to minimize load imbalance on the node
  call MAPL_AllocNodeArray(indices,(/clrm/),rc=ierr)
  if (MAPL_am_I_root())  indices = knuthshuffle(clrm)
  call MAPL_SyncSharedMemory(rc=ierr)

! Loop Through Bands
! ------------------------------  
do iband = 1, nband

  ! Read in channels for band
  ! --------------------------
    if (iband == 1) then 
      bandname = 'blue'
      nch = 60
    else if (iband == 2) then
      bandname = 'red'
      nch = 60
    else if (iband == 3) then
      bandname = 'SWIR'
      nch = 9
    end if
    if (MAPL_am_I_root()) then
      write(*,*) '<> Working on ',bandname,' band'
    end if
    allocate(channels(nch))
    call readvar1Dgrp(trim(bandname)//"_wavelength", "sensor_band_parameters", ANG_file, channels)

  ! Allocate arrays that will be shared
  ! -----------------------------------------------------------------
    call allocate_shared_band()
    call MAPL_SyncSharedMemory(rc=ierr)

  ! Allocate arrays that will be copied on each processor - unshared
  ! -----------------------------------------------------------------
    call allocate_unshared()


  ! Initialize outputs to be safe
  ! -------------------------------
    LTAU           = g5nr_missing
    LSSA           = g5nr_missing
    LG             = g5nr_missing
    ITAU           = g5nr_missing
    ISSA           = g5nr_missing
    IG             = g5nr_missing
    call MAPL_SyncSharedMemory(rc=ierr)
  ! Main do loop over the part of the shuffled domain assinged to each processor
    write(*,*)'<> Work on all channels ',myid
    progress = 0
    do cc = starti, endi  
      c = indices(cc)
      i  = iIndex(c)
      j  = jIndex(c)

  !   Cloud Optical Properties
  !   ------------------------
      ! nch = 1
      call getCOPvector(IcldTable, km, nobs, nch, nMom, nPol, channels, REI(i,j,:), TAUI(i,j,:), VssaIcl, VgIcl, VpmomIcl, VtauIcl)
      call getCOPvector(LcldTable, km, nobs, nch, nMom, nPol, channels, REL(i,j,:), TAUL(i,j,:), VssaLcl, VgLcl, VpmomLcl, VtauLcl)

  !   Save some variables on the 2D Grid for Writing Later
  !   ------------------------
                             !km,nch,nobs
      LTAU(i,j,:) = SUM(VtauLcl(:,:,nobs),DIM=1)
      LSSA(i,j,:) = SUM(VssaLcl(:,:,nobs)*VtauLcl(:,:,nobs),DIM=1)
      LG(i,j,:)   = SUM(VgLcl(:,:,nobs)*VtauLcl(:,:,nobs),DIM=1)      

      ITAU(i,j,:) = SUM(VtauIcl(:,:,nobs),DIM=1)
      ISSA(i,j,:) = SUM(VssaIcl(:,:,nobs)*VtauIcl(:,:,nobs),DIM=1)
      IG(i,j,:)   = SUM(VgIcl(:,:,nobs)*VtauIcl(:,:,nobs),DIM=1)    
      do ch = 1,nch
        if (LTAU(i,j,ch) > 0) then
          LSSA(i,j,ch) = LSSA(i,j,ch)/LTAU(i,j,ch)
          LG(i,j,ch)   = LG(i,j,ch)/LTAU(i,j,ch)
        else
          LSSA(i,j,ch) = g5nr_missing
          LG(i,j,ch)   = g5nr_missing
        end if
        
        if (ITAU(i,j,ch) > 0) then        
          ISSA(i,j,ch) = ISSA(i,j,ch)/ITAU(i,j,ch)
          IG(i,j,ch)   = IG(i,j,ch)/ITAU(i,j,ch)
        else
          ISSA(i,j,ch) = g5nr_missing
          IG(i,j,ch)   = g5nr_missing
        end if
      end do
     
  !   Keep track of progress of each processor
  !   -----------------------------------------        
      if (nint(100.*real(cc-starti)/real(counti)) > progress) then         
        write(*,'(A,I,A,I,A,I2,A,I3,A)') 'Pixel: ',cc,'  End Pixel: ',endi,'  ID:',myid,'  Progress:', nint(progress),'%'
        progress = progress + 10       
      end if
                  
    end do ! do clear pixels
    ! Sync processors before shutting down
    ! ------------------------------------
    call MAPL_SyncSharedMemory(rc=ierr)

    !Write cloud data
    !----------------------
    call write_cldfile()

    ! Sync processors before shutting down
    ! ------------------------------------
    call MAPL_SyncSharedMemory(rc=ierr)

    if (MAPL_am_I_root()) then
      write(*,*) '<> Wrote Out File for band ',bandname
      write(*,*) ' '
    end if

    ! All done band
    ! -----------------
    deallocate(channels)
    call deallocate_shared_band()
    call deallocate_unshared()

    call MAPL_SyncSharedMemory(rc=ierr)
  end do !nband
500  call MAPL_SyncSharedMemory(rc=ierr)
  call MAPL_FinalizeShmem (rc=ierr)
  call ESMF_Finalize(__RC__)

! -----------------------------------------------------------------------------------

  contains

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     write_cldfile()
! PURPOSE
!     writes data to cld output file
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     Oct 2018 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

  subroutine write_cldfile()
    integer                            :: ncid, varid, found, length 
    character(len=256)                 :: filename
    integer, dimension(npet)           :: nlayer  ! how many layers each processor writes
    integer                            :: startl,countl,endl
! Write to CLDO_file
! -------------------
    nlayer = 0
    if (npet >= jm) then
      nlayer(1:npet) = 1
    else if (npet < jm) then
      nlayer(1:npet) = jm/npet
      nlayer(npet)   = nlayer(npet) + mod(jm,npet)
    end if 

    if (myid == 0) then
      startl = 1
    else
      startl = sum(nlayer(1:myid))+1
    end if
    countl = nlayer(myid+1)
    endl   = startl + countl - 1  


    found = index(OUT_file,'_{}')
    length = len(trim(OUT_file))
    ! filename = OUT_file(1:found-1)//trim(bandname)//OUT_file(found+2:length)
    filename = OUT_file(1:found-1)//OUT_file(found+2:length)
    call check( nf90_open(filename, IOR(nf90_write, nf90_mpiio), ncid,comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename )
      
    !                             Write to Cloud Outputs File
    !                             --------------------------------
    ! ! Liquid Cloud Stuff
    call check(nf90_inq_varid(ncid, 'lcot_' // trim(bandname), varid), "get lcot vaird")    
    call check(nf90_put_var(ncid, varid, LTAU(:,startl:endl,:),start = (/1,startl,1/), count = (/im,countl,nch/)), "writing out ltau")

    ! ! call check(nf90_inq_varid(ncid, 'lc_g_' // trim(bandname), varid), "get lc_g vaird")
    ! call check(nf90_put_var(ncid, varid, LG), "writing out lg")

    call check(nf90_inq_varid(ncid, 'lc_ssa_' // trim(bandname), varid), "get lc_ssa vaird")
    call check(nf90_put_var(ncid, varid, LSSA(:,startl:endl,:), start = (/1,startl,1/), count = (/im,countl,nch/)), "writing out lssa")

    ! Liquid Cloud Stuff
    call check(nf90_inq_varid(ncid, 'icot_' // trim(bandname), varid), "get icot vaird")
    call check(nf90_put_var(ncid, varid, ITAU(:,startl:endl,:), start = (/1,startl,1/), count = (/im,countl,nch/)), "writing out itau")

    ! ! call check(nf90_inq_varid(ncid, 'ic_g_' // trim(bandname), varid), "get ic_g vaird")
    ! ! call check(nf90_put_var(ncid, varid, IG), "writing out ig")

    call check(nf90_inq_varid(ncid, 'ic_ssa_' // trim(bandname), varid), "get ic_ssa vaird")
    call check(nf90_put_var(ncid, varid, ISSA(:,startl:endl,:), start = (/1,startl,1/), count = (/im,countl,nch/)), "writing out issa")

    if (myid == 0) then
      call check(nf90_inq_varid(ncid, trim(bandname) // '_wavelength', varid), "get wavelength vaird")
      call check(nf90_put_var(ncid,varid,channels), "writing out channels")  
    end if
        
    call check( nf90_close(ncid), "close cldofile" )

  end subroutine write_cldfile


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     when
! PURPOSE
!     return a string with current date nad time
! INPUT
!     NONE
! OUTPUT
!     when: date string with format Date 04/01/1995; time 14:07:40
! HISTORY
!     Oct 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  function when()
        integer,dimension(3)      :: now
        integer,dimension(4)      :: today
        character(len=256)        :: when

        call idate(today(1),today(2),today(3))   ! today(1)=month, (2)=day, (3)=year
        call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
        write ( when, 1000 )  today(1), today(2), today(3)    !, now
        ! 1000 format ( 'Date ', i2.2, '/', i2.2, '/', i2.2, '; time ',&
        !                i2.2, ':', i2.2, ':', i2.2 )

        1000 format ( 'Date ', i2.2, '/', i2.2, '/', i2.2)
        
  end function when



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
!     read_cld_Tau
! PURPOSE
!     read in cloud optical depth
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     August 2017 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine read_cld_Tau()

    call mp_readvar3Dchunk("REI"    , CLD_file, (/im,jm,km/), 1, npet, myid, REI) 
    call mp_readvar3Dchunk("REL"    , CLD_file, (/im,jm,km/), 1, npet, myid, REL) 
    call mp_readvar3Dchunk("TAUI"   , CLD_file, (/im,jm,km/), 1, npet, myid, TAUI) 
    call mp_readvar3Dchunk("TAUL"   , CLD_file, (/im,jm,km/), 1, npet, myid, TAUL) 

    call MAPL_SyncSharedMemory(rc=ierr)    
    if (MAPL_am_I_root()) then
      write(*,*) '<> Read cloud data to shared memory'
    end if      

  end subroutine read_cld_Tau


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
!     allocate_shared
! PURPOSE
!     allocates all the shared memory arrays that I need
! INPUT
!     None
! OUTPUT
!     none
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine allocate_shared()
    if (MAPL_am_I_root()) then
      write(*,*) '<> Allocating shared memory variables'
    end if
    call MAPL_AllocNodeArray(REI,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(REL,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(TAUI,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(TAUL,(/im,jm,km/),rc=ierr)
  end subroutine allocate_shared

  subroutine allocate_shared_band()
    if (MAPL_am_I_root()) then
      write(*,*) '<> Allocating shared memory band variables'
    end if
    call MAPL_AllocNodeArray(LTAU,(/im,jm,nch/),rc=ierr)
    call MAPL_AllocNodeArray(LSSA,(/im,jm,nch/),rc=ierr)
    call MAPL_AllocNodeArray(LG,(/im,jm,nch/),rc=ierr)
    call MAPL_AllocNodeArray(ITAU,(/im,jm,nch/),rc=ierr)
    call MAPL_AllocNodeArray(ISSA,(/im,jm,nch/),rc=ierr)
    call MAPL_AllocNodeArray(IG,(/im,jm,nch/),rc=ierr)

  end subroutine allocate_shared_band

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     deallocate_shared_band
! PURPOSE
!     deallocates all the shared memory arrays that I need
! INPUT
!     None
! OUTPUT
!     none
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine deallocate_shared_band()
    if (MAPL_am_I_root()) then
      write(*,*) '<> Dellocating shared band memory variables'
    end if
    call MAPL_DeallocNodeArray(LTAU,rc=ierr)
    call MAPL_DeallocNodeArray(LSSA,rc=ierr)
    call MAPL_DeallocNodeArray(LG,rc=ierr)
    call MAPL_DeallocNodeArray(ITAU,rc=ierr)
    call MAPL_DeallocNodeArray(ISSA,rc=ierr)
    call MAPL_DeallocNodeArray(IG,rc=ierr)

  end subroutine deallocate_shared_band


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
   !nch = 1
   allocate (VtauIcl(km,nch,nobs))
   allocate (VtauLcl(km,nch,nobs))
   allocate (VssaIcl(km,nch,nobs))
   allocate (VssaLcl(km,nch,nobs))    
   allocate (VgIcl(km,nch,nobs))    
   allocate (VgLcl(km,nch,nobs))    
   allocate (VpmomLcl(km,nch,nobs,nMom,nPol))
   allocate (VpmomIcl(km,nch,nobs,nMom,nPol))

  end subroutine allocate_unshared

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     deallocate_unshared
! PURPOSE
!     allocates all the arrays that I need, but each processor has its own copy
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine deallocate_unshared()
   deallocate (VtauIcl)
   deallocate (VtauLcl)
   deallocate (VssaIcl)
   deallocate (VssaLcl)    
   deallocate (VgIcl)  
   deallocate (VgLcl)    
   deallocate (VpmomLcl)
   deallocate (VpmomIcl)

  end subroutine deallocate_unshared

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     CREATE_CLDFILE
! PURPOSE
!     creates outfile for wiriting to
! INPUT
!     filename
! OUTPUT
!     none
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

  subroutine create_cldfile()

    integer                            :: LssaVarID, LtauVarID, LgVarID
    integer                            :: IssaVarID, ItauVarID, IgVarID
    
    integer                            :: ncid
    integer                            :: ewDimID, nsDimID, levDimID, chDimID 
    integer                            :: scantimeVarID, clonVarID, clatVarID
    integer                            :: chVarID, levVarID, ewVarID, nsVarID
    integer                            :: leveVarID, e

    real*8,allocatable,dimension(:,:)  :: clon, clat
    real*8,allocatable,dimension(:)    :: scantime, ew, ns, tyme, lev

    character(len=2000)                :: comment
    character(len=256)                 :: filename
    integer                            :: found, length

    found = index(OUT_file,'_{}')
    length = len(trim(OUT_file))
    ! filename = OUT_file(1:found-1)//trim(bandname)//OUT_file(found+2:length)
    filename = OUT_file(1:found-1)//OUT_file(found+2:length)
    write(*,*) 'Creating file: ',trim(filename)

!                           CLOUD OUTPUTS
!                           ------------------
      ! Open File
      call check(nf90_create(filename, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // filename)

      ! Create dimensions
      ! call check(nf90_def_dim(ncid, "lev", km, levDimID), "creating ns dimension") !km
      call check(nf90_def_dim(ncid, "ccd_pixels", im, ewDimID), "creating ew dimension") !im
      call check(nf90_def_dim(ncid, "number_of_scans", jm, nsDimID), "creating ns dimension") !jm

      ! Global Attributes
      write(comment,'(A)') 'Atmospheric inputs for VLIDORT Simulation of GEOS-5 PACE Sampler'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'title',trim(comment)),"title attr")

      write(comment,'(A)') 'NASA/Goddard Space Flight Center'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'institution',trim(comment)),"institution attr")

      write(comment,'(A)') 'Global Model and Assimilation Office'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'source',trim(comment)),"source attr")

      write(comment,'(A)') 'n/a'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'references',trim(comment)),"references attr") 

      write(comment,'(A)') 'This file contains intermediate input data for the VLIDORT simulation'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'comment',trim(comment)),"comment attr")   


      call check(nf90_put_att(ncid,NF90_GLOBAL,"contact","Patricia Castellanos <patricia.castellanos@nasa.gov>"),"contact attr")
      call check(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","cf"),"conventions attr")

      ! Define Variables
  !                                     Dimensions
  !                                     ----------    
      ! call check(nf90_def_var(ncid,'lev',nf90_float,(/levDimID/),levVarID),"create lev var")
      call check(nf90_def_var(ncid,'ccd_pixels',nf90_float,(/ewDimID/),ewVarID),"create ew var")
      call check(nf90_def_var(ncid,'number_of_scans',nf90_float,(/nsDimID/),nsVarID),"create ns var")

      call check(nf90_def_var(ncid,'ev_mid_time',nf90_float,(/ewDimID/),scantimeVarID),"create scanTime var")
      call check(nf90_def_var(ncid,'longitude',nf90_float,(/ewDimID,nsDimID/),clonVarID),"create clon var")
      call check(nf90_def_var(ncid,'latitude',nf90_float,(/ewDimID,nsDimID/),clatVarID),"create clat var")

    do iband = 1, nband

      ! Read in channels for band
      ! --------------------------
        if (iband == 1) then 
          bandname = 'blue'
          nch = 60
        else if (iband == 2) then
          bandname = 'red'
          nch = 60
        else if (iband == 3) then
          bandname = 'SWIR'
          nch = 9
        end if

        
        call check(nf90_def_dim(ncid, trim(bandname)//"_wavelength", nch, chDimID), "creating ch dimension") !nch
        call check(nf90_def_var(ncid,trim(bandname)//"_wavelength",nf90_float,(/chDimID/),chVarID),"create ch var")

  !                                     Data
  !                                     ----
        call check(nf90_def_var(ncid, 'lcot_' // trim(bandname) ,nf90_float,(/ewDimID,nsDimID,chDimID/),LtauVarID),"create lcot var")      
        call check(nf90_def_var(ncid, 'lc_ssa_' // trim(bandname) ,nf90_float,(/ewDimID,nsDimID,chDimID/),LssaVarID),"create lc_ssa var")
  !      call check(nf90_def_var(ncid, 'lc_g_' // trim(bandname),nf90_float,(/ewDimID,nsDimID,levDimID,chDimID/),LgVarID),"create lc_g var")

        call check(nf90_def_var(ncid, 'icot_' // trim(bandname) ,nf90_float,(/ewDimID,nsDimID,chDimID/),ItauVarID),"create icot var")      
        call check(nf90_def_var(ncid, 'ic_ssa_' // trim(bandname) ,nf90_float,(/ewDimID,nsDimID,chDimID/),IssaVarID),"create ic_ssa var")
  !      call check(nf90_def_var(ncid, 'ic_g_' // trim(bandname) ,nf90_float,(/ewDimID,nsDimID,levDimID,chDimID/),IgVarID),"create ic_g var")


    ! Variable Attributes
    !                                          Cloud Data
    !                                          -----------------  
        ! Liquid Cloud Stuff
        write(comment,'(A,A)') bandname, ' wavelengths Liquid COT'
        call check(nf90_put_att(ncid,LtauVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(A,A)') bandname, ' wavelengths layer Liquid Cloud Optical Thickness'
        call check(nf90_put_att(ncid,LtauVarID,'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,LtauVarID,'missing_value',real(g5nr_missing)),"missing_value attr")
        call check(nf90_put_att(ncid,LtauVarID,'units','none'),"units attr")
        call check(nf90_put_att(ncid,LtauVarID,"_FillValue",real(g5nr_missing)),"_Fillvalue attr")

        write(comment,'(A,A)') bandname, ' wavelengths Liquid Cloud SSA'
        call check(nf90_put_att(ncid,LssaVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(A,A)') bandname, ' wavelengths layer Liquid Cloud Single Scattering Albedo'
        call check(nf90_put_att(ncid,LssaVarID,'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,LssaVarID,'missing_value',real(g5nr_missing)),"missing_value attr")
        call check(nf90_put_att(ncid,LssaVarID,'units','none'),"units attr")
        call check(nf90_put_att(ncid,LssaVarID,"_FillValue",real(g5nr_missing)),"_Fillvalue attr")   

    !       write(comment,'(A,A)') bandname, ' wavelengths Liquid Cloud g'
    !       call check(nf90_put_att(ncid,LgVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
    !       write(comment,'(A,A)') bandname, ' wavelengths layer Liquid Cloud Asymmetry Parameter'
    !       call check(nf90_put_att(ncid,LgVarID,'long_name',trim(adjustl(comment))),"long_name attr")
    !       call check(nf90_put_att(ncid,LgVarID,'missing_value',real(g5nr_missing)),"missing_value attr")
    !       call check(nf90_put_att(ncid,LgVarID,'units','none'),"units attr")
    !       call check(nf90_put_att(ncid,LgVarID,"_FillValue",real(g5nr_missing)),"_Fillvalue attr")  

        ! Ice Cloud Stuff
        write(comment,'(A,A)') bandname, ' wavelengths Ice COT'
        call check(nf90_put_att(ncid,ItauVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(A,A)') bandname, ' wavelengths layer Ice Cloud Optical Thickness'
        call check(nf90_put_att(ncid,ItauVarID,'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,ItauVarID,'missing_value',real(g5nr_missing)),"missing_value attr")
        call check(nf90_put_att(ncid,ItauVarID,'units','none'),"units attr")
        call check(nf90_put_att(ncid,ItauVarID,"_FillValue",real(g5nr_missing)),"_Fillvalue attr")

        write(comment,'(A,A)') bandname, ' wavelengths Ice Cloud SSA'
        call check(nf90_put_att(ncid,IssaVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(A,A)') bandname, ' wavelengths layer Ice Cloud Single Scattering Albedo'
        call check(nf90_put_att(ncid,IssaVarID,'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,IssaVarID,'missing_value',real(g5nr_missing)),"missing_value attr")
        call check(nf90_put_att(ncid,IssaVarID,'units','none'),"units attr")
        call check(nf90_put_att(ncid,IssaVarID,"_FillValue",real(g5nr_missing)),"_Fillvalue attr")   

    !       write(comment,'(A,A)') bandname, ' wavelengths Ice Cloud g'
    !       call check(nf90_put_att(ncid,IgVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
    !       write(comment,'(A,A)') bandname, ' wavelengths layer Ice Cloud Asymmetry Parameter'
    !       call check(nf90_put_att(ncid,IgVarID,'long_name',trim(adjustl(comment))),"long_name attr")
    !       call check(nf90_put_att(ncid,IgVarID,'missing_value',real(g5nr_missing)),"missing_value attr")
    !       call check(nf90_put_att(ncid,IgVarID,'units','none'),"units attr")
    !       call check(nf90_put_att(ncid,IgVarID,"_FillValue",real(g5nr_missing)),"_Fillvalue attr")  
    end do !bands

  !                                          scanTime
  !                                          -------  
      call check(nf90_put_att(ncid,scantimeVarID,'long_name','Earth view mid time (seconds of day)'),"long_name attr")
      call check(nf90_put_att(ncid,scantimeVarID,'units','seconds'),"units attr")

  !                                          EW, NS, LEV, TIME
  !                                          -----------------------  
      call check(nf90_put_att(ncid,ewVarID,'long_name','pseudo longitude'),"long_name attr")
      call check(nf90_put_att(ncid,ewVarID,'units','degrees_east'),"units attr")
      call check(nf90_put_att(ncid,nsVarID,'long_name','pseudo latitude'),"long_name attr")
      call check(nf90_put_att(ncid,nsVarID,'units','degrees_north'),"units attr")   
      ! call check(nf90_put_att(ncid,levVarID,'long_name','Vertical Level'),"long_name attr")
      ! call check(nf90_put_att(ncid,levVarID,'units','layer'),"units attr")
      ! call check(nf90_put_att(ncid,levVarID,'positive','down'),"positive attr")
      ! call check(nf90_put_att(ncid,levVarID,'axis','z'),"axis attr")

  !                                          clon & clat
  !                                          -------  
      call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
      call check(nf90_put_att(ncid,clonVarID,'missing_value',real(g5nr_missing)),"missing_value attr")
      call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
      call check(nf90_put_att(ncid,clatVarID,'missing_value',real(g5nr_missing)),"missing_value attr")

      call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
      call check(nf90_put_att(ncid,clonVarID,'missing_value',real(g5nr_missing)),"missing_value attr")
      call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
      call check(nf90_put_att(ncid,clatVarID,'missing_value',real(g5nr_missing)),"missing_value attr")  

      call check(nf90_put_att(ncid,chVarID,'long_name','Band center wavelengths'),"long_name attr")
      call check(nf90_put_att(ncid,chVarID,'units','nm'),"units attr")
        
      !Leave define mode
      call check(nf90_enddef(ncid),"leaving define mode")

      ! write out ew, ns, lev, time, clon, clat, & scantime
      allocate (scantime(im))
      allocate (clon(im, jm))
      allocate (clat(im, jm))
      allocate (ew(im))
      allocate (ns(jm))    
      ! allocate (lev(km))

      call readvar1D("ev_mid_time", CLD_file, scantime)
      call check(nf90_put_var(ncid,scantimeVarID,scantime), "writing out scantime")

      call readvar2D("longitude", CLD_file, clon)
      call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

      call readvar2D("latitude", CLD_file, clat)
      call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

      ! call readvar1D("lev", CLD_file, lev)
      ! call check(nf90_put_var(ncid,levVarID,lev), "writing out lev")

      call readvar1D("ccd_pixels", CLD_file, ew)
      call check(nf90_put_var(ncid,ewVarID,ew), "writing out ew")

      call readvar1D("number_of_scans", CLD_file, ns)
      call check(nf90_put_var(ncid,nsVarID,ns), "writing out ns")   

      deallocate (clon)
      deallocate (clat)  
      deallocate (scantime)
      deallocate (ns)
      deallocate (ew)
      ! deallocate (lev)

      call check( nf90_close(ncid), "close cldofile" )
  end subroutine create_cldfile  




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
    ! Clouds
    call ESMF_ConfigGetAttribute(cf, IcldTable, label = 'ICLDTABLE:',__RC__)       
    call ESMF_ConfigGetAttribute(cf, LcldTable, label = 'LCLDTABLE:',__RC__)

    ! Input Files
    call ESMF_ConfigGetAttribute(cf, OUT_file, label = 'OUT_file:',__RC__)
    call ESMF_ConfigGetAttribute(cf, CLD_file, label = 'CLD_file:',__RC__)
    call ESMF_ConfigGetAttribute(cf, ANG_file, label = 'ANG_file:',__RC__)

    call ESMF_ConfigDestroy(cf)

  end subroutine get_config

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    knuthshuffle
! PURPOSE
!     Implements Knuth Shuffle for array indices 1-N
! INPUT
!     N : max index
! OUTPUT
!     knuthshuffle: array of shuffled indices
!  HISTORY
!     November 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  function knuthshuffle(N)
    integer, intent(in)               :: N
    integer,dimension(N)              :: knuthshuffle
    integer                           :: randpos, i, temp
    real                              :: r

    knuthshuffle = (/(i,i=1,N)/)

    do i = N, 2, -1
      call random_number(r)
      randpos = int(r*i) + 1
      temp    = knuthshuffle(randpos)
      knuthshuffle(randpos) = knuthshuffle(i)
      knuthshuffle(i)       = temp
    end do

  end function knuthshuffle




end program cloud_calculator
