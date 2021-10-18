!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    leo_vlidort_cloud
! PURPOSE
!     Reads in parallel the model data (in a netcdf file) interpolated to satellite grid 
!     The variables are needed as input to vlidort
!     A shared memory array created with a call to MAPL_ShmemMod is used for the variables
!     Do some filtering and run the vlidort code
! INPUT
!     rcFile: contains all the input parameters
! OUTPUT
!     None
!  HISTORY
!     October 2018 P. Castellanos adapted from geo_vlidort_cloud.F90
!     27 April 2015 P. Castellanos adapted from A. da Silva shmem_reader.F90
!     Aug 2017 P. Castellanos add clouds
! NOTE
!     VLIDORT should be configured to run with the same number of streams found 
!     in the ice cloud optical properties file.  This is because a truncation
!     factor is used to scale ice optical thickness.  This truncation factor corrects
!     for the loss of the forward peak in the ice cloud scatterimg matrix introduced by
!     using a finite number of streams.
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#  define I_AM_MAIN
#  include "MAPL_Generic.h"
#  include "MAPL_ErrLogMain.h"
program pace_vlidort

  use ESMF                         ! ESMF modules
  use MAPL_Mod
  use MAPL_ShmemMod                ! The SHMEM infrastructure
  use netcdf                       ! for reading the NR files
  use vlidort_brdf_modis           ! Module to run VLIDORT with MODIS BRDF surface supplement
  use vlidort_brdf_modis_bpdf      ! Module to run VLIDORT with MODIS BRDF+BPDF surface supplement  
  use vlidort_lambert              ! Module to run VLIDORT with lambertian surface  
  use vlidort_brdf_cx              ! Module to run VLIDORT with cox munk surface
  use vlidort_rot                  ! Module to calculate Rayleigh
  use Chem_MieMod
  use mp_netcdf_Mod
  use netcdf_Mod

  implicit none
  include "accp_polar_vlidort_pars.F90"
  include "mpif.h"

! ESMF Objects
!----------------
  type(ESMF_Config)       :: cf
  type(ESMF_VM)           :: vm 

! Information to be retrieved from resource file 
! ------------------------------------------------
  character(len=256)                    :: arg                    ! command line rc file argument
  character(len=256)                    :: nodenumarg             ! command line node number argument
  character(len=8)                      :: date         
  character(len=6)                      :: time 
  character(len=256)                    :: instname
  character(len=256)                    :: landname, landmodel, watername, watermodel
  integer                               :: landbandmBRDF, landbandmLER              ! number of wavelength bands or channels in surface reflectance data file
  real, allocatable                     :: landband_cBRDF(:), landband_cLER(:)      ! modis band center wavelength
  logical                               :: scalar
  integer                               :: nstreams               ! number of half space streams, default = 12
  logical                               :: plane_parallel    
  logical                               :: aerosol_free
  real, allocatable                     :: channels(:)            ! channels to simulate
  real*8, allocatable                   :: mr(:)                  ! water real refractive index    
  integer                               :: nch                    ! number of channels  
  real                                  :: szamax                 ! Geomtry filtering
  character(len=256)                    :: rcfile                 ! resource file

! Test flag
! -----------
  logical                               :: test_shmem = .False.

! File names
! ----------
  character(len=256)                    :: AER_file, ANG_file, INV_file, BRDF_file, LER_file, OUT_file
  character(len=256)                    :: MET_file, NDVI_file, BPDF_file 

! Global, 3D inputs to be allocated using SHMEM
! ---------------------------------------------
  real, pointer                         :: AIRDENS(:,:) => null()
  real, pointer                         :: RH(:,:) => null()
  real, pointer                         :: DELP(:,:) => null()
  real, pointer                         :: PS(:) => null()  
  real, pointer                         :: DU001(:,:) => null()
  real, pointer                         :: DU002(:,:) => null()
  real, pointer                         :: DU003(:,:) => null()
  real, pointer                         :: DU004(:,:) => null()
  real, pointer                         :: DU005(:,:) => null()
  real, pointer                         :: SS001(:,:) => null()
  real, pointer                         :: SS002(:,:) => null()
  real, pointer                         :: SS003(:,:) => null()
  real, pointer                         :: SS004(:,:) => null()
  real, pointer                         :: SS005(:,:) => null()
  real, pointer                         :: BCPHOBIC(:,:) => null()
  real, pointer                         :: BCPHILIC(:,:) => null()
  real, pointer                         :: OCPHOBIC(:,:) => null()
  real, pointer                         :: OCPHILIC(:,:) => null()
  real, pointer                         :: SO4(:,:) => null()  
  real, pointer                         :: KISO(:,:) => null()
  real, pointer                         :: KVOL(:,:) => null()
  real, pointer                         :: KGEO(:,:) => null() 
  real, pointer                         :: LER(:,:) => null()   
  real, pointer                         :: FRLAND(:) => null()
  real, pointer                         :: SZA(:,:) => null()
  real, pointer                         :: VZA(:,:) => null()
  real, pointer                         :: SAA(:,:) => null()
  real, pointer                         :: VAA(:,:) => null()
  real, pointer                         :: RAA(:,:) => null()
  real, pointer                         :: V10M(:) => null()  
  real, pointer                         :: U10M(:) => null()  
  real, pointer                         :: NDVI(:) => null()  
  real, pointer                         :: BPDFcoef(:) => null()  
  integer, pointer                      :: indices(:) => null()
  real, pointer                         :: READER2D(:,:) => null()
  real, pointer                         :: READER1D(:) => null()


! VLIDORT input arrays
! ---------------------------
  real, allocatable                     :: Vpe(:,:)                ! edge pressure [Pa]
  real, allocatable                     :: Vze(:,:)                ! edge height above sfc [m]
  real, allocatable                     :: Vte(:,:)                ! edge Temperature [K]
  real, allocatable                     :: Vqm(:,:,:)              ! (mixing ratio) * delp/g
  real, allocatable                     :: Vtau(:,:,:)             ! aerosol optical depth
  real, allocatable                     :: Vssa(:,:,:)             ! single scattering albedo
  real, allocatable                     :: Vg(:,:,:)               ! asymmetry factor
  real*8, allocatable                   :: Valbedo(:,:,:)            ! surface albedo

! VLIDORT output arrays
!-------------------------------
!                                  Intermediate Unshared Arrays
!                                  -----------------------------
  real*8, allocatable                   :: radiance_VL_int(:,:,:)                   ! TOA normalized radiance from VLIDORT
  real*8, allocatable                   :: reflectance_VL_int(:,:,:)                ! TOA reflectance from VLIDORT  
  real*8, allocatable                   :: Q_int(:,:,:)                             ! Q Stokes component
  real*8, allocatable                   :: U_int(:,:,:)                             ! U Stokes component
  real*8, allocatable                   :: ROT_int(:,:,:)                         ! rayleigh optical thickness
  real*8, allocatable                   :: depol(:)                               ! rayleigh depolarization ratio
  real*8, allocatable                   :: alpha(:,:,:)                           ! trace gas absorption
  real*8, allocatable                   :: flux_factor(:,:)                       ! solar irradiance (F0)
  real*8, allocatable                   :: BR_Q_int(:,:,:)                          ! surface albedo Q
  real*8, allocatable                   :: BR_U_int(:,:,:)                          ! surface albedo U


!                                  Final Shared Arrays
!                                  -------------------
  real*8, pointer                       :: radiance_VL(:,:) => null()             ! TOA normalized radiance from VLIDORT
  real*8, pointer                       :: reflectance_VL(:,:) => null()          ! TOA reflectance from VLIDORT
  real*8, pointer                       :: Q(:,:) => null()                      ! Q Stokes component
  real*8, pointer                       :: U(:,:) => null()                      ! U Stokes component
  real*8, pointer                       :: ROT(:,:) => null()                  ! rayleigh optical depth
  real*8, pointer                       :: ALBEDO(:,:) => null()                 ! bi-directional surface reflectance
  real*8, pointer                       :: BR_Q(:,:) => null()                   ! bi-directional surface reflectance Q
  real*8, pointer                       :: BR_U(:,:) => null()                   ! bi-directional surface reflectance U


  real, pointer                         :: TAU(:) => null()                  ! aerosol optical depth
  real, pointer                         :: SSA(:) => null()                  ! single scattering albedo
  real, pointer                         :: G(:) => null()                    ! asymmetry factor

  real, pointer                         :: PE(:,:) => null()
  real, pointer                         :: ZE(:,:) => null()  
  real, pointer                         :: TE(:,:) => null()

! VLIDORT working variables
!------------------------------
  integer                               :: ch                                       ! i-channel  
  integer                               :: iband                                    ! i-surfaceband
  real,allocatable                      :: Vpmom(:,:,:,:,:)                         ! elements of scattering phase matrix for vector calculations

! MODIS Kernel variables
!--------------------------
  real*8, allocatable                   :: kernel_wt(:,:,:)                         ! kernel weights (/fiso,fgeo,fvol/)
  real*8, allocatable                   :: param(:,:,:)                             ! Li-Sparse parameters 
                                                                                    ! param1 = crown relative height (h/b)
                                                                                    ! param2 = shape parameter (b/r)
  real                                  :: land_missing                                                                 
  real*8, allocatable                   :: BPDFparam(:,:,:)                         ! BPDF model parameters (/1.5,NDVI,BPDFcoef/)

! Mie Table Stucture
!---------------------
  type(Chem_Mie)                        :: mieTables

! Satellite domain variables
!------------------------------
  integer                               :: nalong, km, tm                            ! size of satellite domain
  integer                               :: i, j, k, n, ialong                        ! satellite domain working variable
  integer                               :: starti, counti, endi                      ! array indices and counts for each processor
  integer, allocatable                  :: nclr(:)                                   ! how many clear pixels each processor works on
  integer                               :: clrm                                      ! number of clear pixels for this part of decomposed domain
  integer                               :: c, cc                                     ! clear pixel working variable
  real, allocatable                     :: SOLAR_ZENITH(:,:)                         ! solar zenith angles used for data filtering
  logical, allocatable                  :: clmask(:)                                 ! cloud-land mask

! netcdf variables
!----------------------  
  integer                               :: ncid                                           ! netcdf file id
  integer                               :: varid

! Miscellaneous
! -------------
  integer                               :: ierr, rc, status                            ! MPI error message
  integer                               :: status_mpi(MPI_STATUS_SIZE)                 ! MPI status
  integer                               :: myid, npet, CoresPerNode                    ! MPI dimensions and processor id
  logical                               :: amOnFirstNode
  integer                               :: pp                                          ! i-processor
  logical,allocatable,dimension(:)      :: FirstNode
  character(len=100)                    :: msg                                         ! message to be printed
  character(len=256)                    :: line                                        ! dummy variables to read line of file 
  real                                  :: progress                                    ! 
  real                                  :: g5nr_missing  

! System tracking variables
! -----------------------------
  integer*8                             :: t1, t2, clock_max
  real*8                                :: clock_rate
  character(len=*), parameter           :: Iam = 'accp_vlidort'

!                               END OF VARIABLE DECLARATIONS
!----------------------------------------------------------------------------------------------------------
  
  progress = 10

! Initialize MPI with ESMF
! ------------------------
  call ESMF_Initialize (logkindflag=ESMF_LOGKIND_NONE, vm=vm, __RC__)

  call ESMF_VMGet(vm, localPET=myid, PETcount=npet) 
  if ( MAPL_am_I_root() ) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'

! Initialize SHMEM
! ----------------
  CoresPerNode = MAPL_CoresPerNodeGet(MPI_COMM_WORLD,rc=ierr) ! a must
  call MAPL_InitializeShmem(rc=ierr)
  amOnFirstNode = MAPL_ShmemAmOnFirstNode(comm=MPI_COMM_WORLD,rc=ierr)

  if (MAPL_am_I_root()) then
    allocate (FirstNode(npet-1))
    do pp = 1,npet-1
      call mpi_recv(FirstNode(pp), 1, MPI_LOGICAL, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)
    end do
  else
    call mpi_send(amOnFirstNode, 1, MPI_LOGICAL, 0, 2001, MPI_COMM_WORLD, ierr)
  end if  

! Parse Resource file provided at command line for input info 
! -----------------------------------------------------------
  call getarg(1, arg)
  call get_config(arg)

! Write out settings to use
! --------------------------
  if (MAPL_am_I_root()) then
    write(*,*) 'Simulating ', lower_to_upper(trim(instname)),' domain on ',date,' ', time, 'Z'
    write(*,*) 'Channels [nm]: ',channels
    write(*,*) 'SZA < ', szamax
    if (scalar) write(*,*) 'Scalar calculations'
    if (.not. scalar) write(*,*) 'Vector calculations'
    write(*,*) ' '
  end if 

! Query for domain dimensions and missing value
!----------------------------------------------  
  call mp_readDim("time", AER_file,tm)
  call mp_readDim("lev", AER_file, km)
  call mp_readDim("nalong", ANG_file, nalong)


  call mp_readVattr("missing_value", AER_FILE, "DELP", g5nr_missing)
  if ((index(lower_to_upper(landmodel),'RTLS') > 0)) then
    call mp_readVattr("missing_value", BRDF_file, "Riso470", land_missing) 
  else if (lower_to_upper(landmodel) == 'LAMBERTIAN') then
    call mp_readVattr("missing_value", LER_file, "SRFLER354", land_missing) 
  end if

! Allocate arrays that will be copied on each processor - unshared
! -----------------------------------------------------------------
  call allocate_unshared()

! Create OUTFILE
! --------------
  if ( MAPL_am_I_root() )  call create_outfile(date, time)


  call MAPL_SyncSharedMemory(rc=ierr)

! Read the land and angle data 
! -------------------------------------
  call read_sza()

! Create mask
! ------------------------
  ! first check that you can do anything!!!!!  
  if (.not. ANY(SOLAR_ZENITH < szamax)) then   
    if (MAPL_am_I_root()) then
      write(*,*) 'The sun has set, nothing to do'
      write(*,*) 'Exiting.....'
    end if
    GOTO 500
  end if

! Figure out how many indices to work on
!------------------------------------------
  clrm = 0
  clmask = .False.
  do i=1,tm
    if (ALL(SOLAR_ZENITH(:,i) < szamax)) then
      clrm = clrm + 1
      clmask(i) = .True.
    end if
  end do

  if (clrm == 0) then
    if (MAPL_am_I_root()) then
      write(*,*) 'No good pixels, nothing to do'
      write(*,*) 'Exiting.....'
    end if
    GOTO 500
  end if


! Unshared angles no longer needed
!---------------------------------------
  deallocate (SOLAR_ZENITH)

  if (MAPL_am_I_root()) then
    write(*,*) '<> Created Mask'
    write(*,'(A,I3,A)') '       ',nint(100.*clrm/(tm)),'% of the domain will be simulated'
    write(*,'(A,I,A)')  '       Simulating ',clrm,' pixels'
    write(*,*) ' '
  end if   

! Allocate the Global arrays using SHMEM
! It will be available on all processors
! ---------------------------------------------------------
  if (amOnFirstNode) then
    call allocate_shared(clrm)
  end if
  call MAPL_SyncSharedMemory(rc=ierr)
  
! Read in the global arrays
! ------------------------------
  call read_land()
  call read_aer_Nv()
  call read_surf_land()
  call read_angles()
  call read_wind()

! Wait for everyone to finish reading 
! ------------------------------------------------------------------  
  call MAPL_SyncSharedMemory(rc=ierr)
   

! Split up filtered domain among processors
! Root processor does not work. Reserved for message passing.
! Split up among npet-1 processors
!----------------------------------------------
  nclr = 0
  if (npet-1 >= clrm) then
    nclr(1:clrm) = 1
  else if (npet-1 < clrm) then
    nclr(1:npet-1) = clrm/(npet-1)
    do i=1,mod(clrm,npet-1)
      nclr(i)   = nclr(i) + 1
    end do
  end if 

! Create unshared arrays on other nodes
! --------------------------------------
  if (.not. amOnFirstNode) then
    call allocate_multinode(nclr(myid))
  end if

! Initialize outputs to be safe
! -------------------------------
  TAU            = dble(MISSING)
  SSA            = dble(MISSING)
  G              = dble(MISSING)  
  ROT            = dble(MISSING)
  PE             = dble(MISSING)
  ZE             = dble(MISSING)
  TE             = dble(MISSING) 
  ALBEDO         = dble(MISSING)   
  radiance_VL    = dble(MISSING)
  reflectance_VL = dble(MISSING)
  if (.not. scalar) then
    Q = dble(MISSING)
    U = dble(MISSING)
    BR_Q = dble(MISSING)
    BR_U = dble(MISSING)    
  end if
  call MAPL_SyncSharedMemory(rc=ierr)

! Prepare inputs and run VLIDORT
! -----------------------------------
  call strarr_2_chararr(vnames_string,nq,16,vnames)

! Create the Mie Tables
! ---------------------
  mieTables = Chem_MieCreate(trim(rcfile),rc)
  if ( rc /= 0 ) then
    print *, 'Cannot create Mie tables from '//trim(rcfile)
    call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
  end if

  if ( nMom > mieTables%nMom ) then ! mieTables%nMom is writen in Aod_EOS.rc file
    print *, 'mieTables do not have enough moments', nMom, mieTables%nMom
    call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
  end if

! Get start and end index for each processor
! start from process 1 as root (0) does not work
! -------------------------------------------
  if (myid == 1) then
    starti  = 1
  else
    starti = sum(nclr(1:myid-1)) + 1
  end if    
  counti = nclr(myid)
  endi   = starti + counti - 1 


! Processors not on the root node need to ask for data
! -----------------------------------------------------
  if (MAPL_am_I_root()) then
    do pp = 1,npet-1
      if (.not. FirstNode(pp)) then
        starti = sum(nclr(1:pp-1))+1
        counti = nclr(pp)
        endi   = starti + counti - 1  
        if (counti > 0) then
          call mpi_send(AIRDENS(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(RH(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(DELP(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(DU001(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(DU002(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(DU003(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(DU004(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(DU005(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(SS001(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(SS002(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(SS003(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(SS004(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(SS005(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)    
          call mpi_send(BCPHOBIC(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr) 
          call mpi_send(BCPHILIC(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr) 
          call mpi_send(OCPHOBIC(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr) 
          call mpi_send(OCPHILIC(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)
          call mpi_send(SO4(starti:endi,:), counti*km, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)   
          if ((index(lower_to_upper(landmodel),'RTLS') > 0)) then
            call mpi_send(KISO(starti:endi,:), counti*landbandmBRDF, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr) 
            call mpi_send(KVOL(starti:endi,:), counti*landbandmBRDF, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr) 
            call mpi_send(KGEO(starti:endi,:), counti*landbandmBRDF, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr) 
          end if
          if ((lower_to_upper(landmodel) == 'LAMBERTIAN') .or. (index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0)) then
            call mpi_send(LER(starti:endi,:), counti*landbandmLER, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)   
          end if
          if ( (index(lower_to_upper(landmodel),'BPDF') > 0) ) then
            if ( (index(lower_to_upper(landname),'MAIGNAN') > 0) ) then
              call mpi_send(NDVI(starti:endi), counti, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)  
              call mpi_send(BPDFcoef(starti:endi), counti, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)  
            end if
          end if 
          call mpi_send(SZA(starti:endi,:), counti*nalong, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)  
          call mpi_send(VZA(starti:endi,:), counti*nalong, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)  
          call mpi_send(RAA(starti:endi,:), counti*nalong, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)  

          call mpi_send(FRLAND(starti:endi), counti, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)    
          call mpi_send(U10M(starti:endi), counti, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr)   
          call mpi_send(V10M(starti:endi), counti, MPI_REAL, pp, 2001, MPI_COMM_WORLD, ierr) 

        end if
      end if 
    end do
  end if 

  if (.not. amOnFirstNode) then
    if (counti > 0) then
      call mpi_recv(AIRDENS, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(RH, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(DELP, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(DU001, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(DU002, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(DU003, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(DU004, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(DU005, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(SS001, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(SS002, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(SS003, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(SS004, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
      call mpi_recv(SS005, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)    
      call mpi_recv(BCPHOBIC, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr) 
      call mpi_recv(BCPHILIC, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)  
      call mpi_recv(OCPHOBIC, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr) 
      call mpi_recv(OCPHILIC, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr) 
      call mpi_recv(SO4, counti*km, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)      
      if ((index(lower_to_upper(landmodel),'RTLS') > 0)) then
        call mpi_recv(KISO, counti*landbandmBRDF, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)  
        call mpi_recv(KVOL, counti*landbandmBRDF, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)
        call mpi_recv(KGEO, counti*landbandmBRDF, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr) 
      end if
      if ((lower_to_upper(landmodel) == 'LAMBERTIAN') .or. (index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0)) then
        call mpi_recv(LER, counti*landbandmLER, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr) 
      end if
      if ( (index(lower_to_upper(landmodel),'BPDF') > 0) ) then
        if ( (index(lower_to_upper(landname),'MAIGNAN') > 0) ) then
          call mpi_recv(NDVI, counti, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr) 
          call mpi_recv(BPDFcoef, counti, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr) 
        end if
      end if

      call mpi_recv(SZA, counti*nalong, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)   
      call mpi_recv(VZA, counti*nalong, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)   
      call mpi_recv(RAA, counti*nalong, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)      

      call mpi_recv(FRLAND, counti, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)    
      call mpi_recv(U10M, counti, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)   
      call mpi_recv(V10M, counti, MPI_REAL, 0, 2001, MPI_COMM_WORLD, status_mpi, ierr)   

    end if
  end if  


! Main do loop over the part of the domain assinged to each processor
  if ((.not. MAPL_am_I_root()) .and. (counti > 0)) then
    if (.not. amOnFirstNode) then
      starti = 1
      endi = counti
    end if
    do c = starti, endi  

      call getEdgeVars ( km, nobs, reshape(AIRDENS(c,:),(/km,nobs/)), &
                       reshape(DELP(c,:),(/km,nobs/)), ptop, &
                       Vpe, Vze, Vte )   

      PE(c,:) = Vpe(:,nobs)
      ZE(c,:) = Vze(:,nobs)
      TE(c,:) = Vte(:,nobs)

      write(msg,'(A,I)') 'getEdgeVars ', myid
      call write_verbose(msg)
      
      call calc_qm()

      write(msg,'(A,I)') 'calc_qm ', myid
      call write_verbose(msg)  

  !   Rayleigh Optical Thickness
  !   --------------------------
      call VLIDORT_ROT_CALC (km, nch, nobs, dble(channels), dble(Vpe), dble(Vze), dble(Vte), &
                                     dble(MISSING),verbose, &
                                     ROT_int, depol, ierr ) 


  !   Aerosol Optical Properties
  !   --------------------------
      call VLIDORT_getAOPvector ( mieTables, km, nobs, nch, nq, channels, vnames, verbose, &
                        Vqm, reshape(RH(c,:),(/km,nobs/)),&
                        nMom,nPol, Vtau, Vssa, Vg, Vpmom, ierr )
      if ( ierr /= 0 ) then
        print *, 'cannot get aerosol optical properties'
        call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)      
      end if
      ! Multiply by -1 to go from Mischenko convention to VLIDORT
      Vpmom(:,:,:,:,2) = -1.*Vpmom(:,:,:,:,2)
      Vpmom(:,:,:,:,4) = -1.*Vpmom(:,:,:,:,4)

  !   Save some variables for Writing Later
  !   ----------------------------------------
      ! Rayleigh
      ROT(c,:) = ROT_int(:,nobs,nch)                                      
      
      ! Aerosols
      TAU(c) = SUM(Vtau(:,nch,nobs))
      SSA(c) = SUM(Vssa(:,nch,nobs)*Vtau(:,nch,nobs))
      G(c)   = SUM(Vg(:,nch,nobs)*Vtau(:,nch,nobs))

      if (TAU(c) > 0) then
        SSA(c) = SSA(c)/TAU(c)
        G(c)   = G(c)/TAU(c)
      else
        SSA(c) = dble(MISSING)
        G(c)   = dble(MISSING)
      end if

      write(msg,*) 'getAOP ', myid
      call write_verbose(msg)


  !   Surface Parameters
  !   -------------------
      if ( FRLAND(c) >= 0.99 ) then
        call get_surf_params()
      end if

  !   Call VlIDORT
  !   ------------
      if ( FRLAND(c) >= 0.99 ) then
        call DO_LAND()
      else
        call DO_OCEAN()
      end if

      
    ! Store Outputs in Shared Arrays
    ! ----------------------------------------------------    
      radiance_VL(c,:)    = radiance_VL_int(nobs,nch,:)
      reflectance_VL(c,:) = reflectance_VL_int(nobs,nch,:)
      ALBEDO(c,:) = Valbedo(nobs,nch,:)

      if (.not. scalar) then
        Q(c,:)      = Q_int(nobs,nch,:)
        U(c,:)      = U_int(nobs,nch,:)
        BR_Q(c,:)   = BR_Q_int(nobs,nch,:)
        BR_U(c,:)   = BR_U_int(nobs,nch,:)
      end if

      write(msg,*) 'VLIDORT Calculations DONE', myid, ierr
      call write_verbose(msg)

  !   Keep track of progress of each processor
  !   -----------------------------------------        
      if (nint(100.*real(c-starti)/real(counti)) > progress) then
        progress = progress + 10
        write(*,'(A,I,A,I,A,I4,A,I3,A)') 'Pixel: ',c,'  End Pixel: ',endi,'  ID:',myid,'  Progress:', nint(progress),'%'           
      end if
                  
    end do ! do clear pixels
  end if ! not Root processor
! Wait for everyone to finish calculations
! ----------------------------------------
  call MAPL_SyncSharedMemory(rc=ierr)

! Processors not on root node need to send data back to be put into shared memory
! ----------------------------------------
  if (.not. amOnFirstNode) then
    if (counti > 0) then
      call mpi_send(PE, counti*km+1, MPI_REAL, 0, 2001, MPI_COMM_WORLD, ierr)
      call mpi_send(ZE, counti*km+1, MPI_REAL, 0, 2001, MPI_COMM_WORLD, ierr)
      call mpi_send(TE, counti*km+1, MPI_REAL, 0, 2001, MPI_COMM_WORLD, ierr)
      
      call mpi_send(ROT, counti*km, MPI_REAL8, 0, 2001, MPI_COMM_WORLD, ierr)
      call mpi_send(TAU, counti, MPI_REAL, 0, 2001, MPI_COMM_WORLD, ierr)
      call mpi_send(SSA, counti, MPI_REAL, 0, 2001, MPI_COMM_WORLD, ierr)
      call mpi_send(G, counti, MPI_REAL, 0, 2001, MPI_COMM_WORLD, ierr)

      call mpi_send(ALBEDO, counti*nalong, MPI_REAL8, 0, 2001, MPI_COMM_WORLD, ierr)
      if (.not. scalar) then
        call mpi_send(Q, counti*nalong, MPI_REAL8, 0, 2001, MPI_COMM_WORLD, ierr)
        call mpi_send(U, counti*nalong, MPI_REAL8, 0, 2001, MPI_COMM_WORLD, ierr)
        call mpi_send(BR_Q, counti*nalong, MPI_REAL8, 0, 2001, MPI_COMM_WORLD, ierr)
        call mpi_send(BR_U, counti*nalong, MPI_REAL8, 0, 2001, MPI_COMM_WORLD, ierr)
      end if

      call mpi_send(radiance_VL, counti*nalong, MPI_REAL8, 0, 2001, MPI_COMM_WORLD, ierr)
      call mpi_send(reflectance_VL, counti*nalong, MPI_REAL8, 0, 2001, MPI_COMM_WORLD, ierr)
    end if
  end if

  if (MAPL_am_I_root()) then
    do pp = 1,npet-1
      if (.not. FirstNode(pp)) then
        starti = sum(nclr(1:pp-1))+1
        counti = nclr(pp)
        endi   = starti + counti - 1  
        if (counti > 0) then
          call mpi_recv(PE(starti:endi,:), counti*km+1, MPI_REAL, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)
          call mpi_recv(ZE(starti:endi,:), counti*km+1, MPI_REAL, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)
          call mpi_recv(TE(starti:endi,:), counti*km+1, MPI_REAL, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)

          call mpi_recv(ROT(starti:endi,:), counti*km, MPI_REAL8, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)
          call mpi_recv(TAU(starti:endi), counti, MPI_REAL, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)
          call mpi_recv(SSA(starti:endi), counti, MPI_REAL, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)
          call mpi_recv(G(starti:endi), counti, MPI_REAL, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)

          call mpi_recv(ALBEDO(starti:endi,:), counti*nalong, MPI_REAL8, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr) 
          if (.not. scalar) then
            call mpi_recv(Q(starti:endi,:), counti*nalong, MPI_REAL8, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)
            call mpi_recv(U(starti:endi,:), counti*nalong, MPI_REAL8, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr) 
            call mpi_recv(BR_Q(starti:endi,:), counti*nalong, MPI_REAL8, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)
            call mpi_recv(BR_U(starti:endi,:), counti*nalong, MPI_REAL8, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr) 
          end if

          call mpi_recv(radiance_VL(starti:endi,:), counti*nalong, MPI_REAL8, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)
          call mpi_recv(reflectance_VL(starti:endi,:), counti*nalong, MPI_REAL8, pp, 2001, MPI_COMM_WORLD, status_mpi, ierr)
        end if
      end if
    end do
  end if
  call MAPL_SyncSharedMemory(rc=ierr)
! Write to outfile
! Expand radiance to im x jm using mask
! Write output to correct position in file 
! -----------------
  if (myid .eq. 1) then
    call write_outfile()
  end if

! Sync processors before shutting down
! ------------------------------------
  call MAPL_SyncSharedMemory(rc=ierr)
  if (MAPL_am_I_root()) then
    write(*,*) '<> Finished VLIDORT Simulation of '//trim(lower_to_upper(instname))//' domain'
    write(*,*) ' '
  end if

! All done
! --------
500  call MAPL_SyncSharedMemory(rc=ierr)
  call MAPL_FinalizeShmem (rc=ierr)
  call ESMF_Finalize(__RC__)

! -----------------------------------------------------------------------------------

  contains

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     write_outfile()
! PURPOSE
!     writes data to outfile
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     Oct 2018 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

  subroutine write_outfile()
    real*8, dimension(tm)            :: field 
    real, dimension(tm)              :: field4  

    field = dble(MISSING)
    field4 = MISSING
!                             Write to main OUT_File
!                             ----------------------
    call check( nf90_open(OUT_file, nf90_write, ncid), "opening file " // OUT_file )
    write(msg,'(F10.2)') channels(nch)
    
    call check(nf90_inq_varid(ncid, 'ROT', varid), "get ref vaird")
    do k=1,km
      call check(nf90_put_var(ncid, varid, real(unpack(reshape(ROT(:,k),(/clrm/)),clmask,field)), &
                  start = (/k,1/), count = (/1,tm/)), "writing out ROT")
    enddo

    call check(nf90_inq_varid(ncid, 'ZE', varid), "get ref vaird")
    do k=1,km+1
      call check(nf90_put_var(ncid, varid, real(unpack(reshape(ZE(:,k),(/clrm/)),clmask,field)), &
                  start = (/k,1/), count = (/1,tm/)), "writing out ZE")
    enddo


    do ialong=1,nalong
      call check(nf90_inq_varid(ncid, 'toa_reflectance', varid), "get ref vaird")
      call check(nf90_put_var(ncid, varid, real(unpack(reshape(reflectance_VL(:,ialong),(/clrm/)),clmask,field)), &
                  start = (/ialong,1/), count = (/1,tm/)), "writing out reflectance")

      call check(nf90_inq_varid(ncid, 'I' , varid), "get rad vaird")
      call check(nf90_put_var(ncid, varid, real(unpack(reshape(radiance_VL(:,ialong),(/clrm/)),clmask,field)), &
                    start = (/ialong,1/), count = (/1,tm/)), "writing out radiance")

      call check(nf90_inq_varid(ncid, 'surf_reflectance', varid), "get ref vaird")
      call check(nf90_put_var(ncid, varid, real(unpack(reshape(ALBEDO(:,ialong),(/clrm/)),clmask,field)), &
                    start = (/ialong,1/), count = (/1,tm/)), "writing out albedo")

      if (.not. scalar) then
        call check(nf90_inq_varid(ncid, 'Q', varid), "get q vaird")
        call check(nf90_put_var(ncid, varid, real(unpack(reshape(Q(:,ialong),(/clrm/)),clmask,field)), &
                    start = (/ialong,1/), count = (/1,tm/)), "writing out Q")

        call check(nf90_inq_varid(ncid, 'U', varid), "get u vaird")
        call check(nf90_put_var(ncid, varid, real(unpack(reshape(U(:,ialong),(/clrm/)),clmask,field)), &
                    start = (/ialong,1/), count = (/1,tm/)), "writing out U")

        call check(nf90_inq_varid(ncid, 'surf_reflectance_Q', varid), "get ref vaird")
        call check(nf90_put_var(ncid, varid, real(unpack(reshape(BR_Q(:,ialong),(/clrm/)),clmask,field)), &
                      start = (/ialong,1/), count = (/1,tm/)), "writing out albedo Q")

        call check(nf90_inq_varid(ncid, 'surf_reflectance_U', varid), "get ref vaird")
        call check(nf90_put_var(ncid, varid, real(unpack(reshape(BR_U(:,ialong),(/clrm/)),clmask,field)), &
                      start = (/ialong,1/), count = (/1,tm/)), "writing out albedo Q")

      endif        

      call check(nf90_inq_varid(ncid, 'solar_zenith', varid), "get sza vaird")
      call check(nf90_put_var(ncid, varid, unpack(reshape(SZA(:,ialong),(/clrm/)),clmask,field4), &
                    start = (/ialong,1/), count = (/1,tm/)), "writing out sza")

      call check(nf90_inq_varid(ncid, 'sensor_zenith', varid), "get vza vaird")
      call check(nf90_put_var(ncid, varid, unpack(reshape(VZA(:,ialong),(/clrm/)),clmask,field4), &
                    start = (/ialong,1/), count = (/1,tm/)), "writing out vza")

      call check(nf90_inq_varid(ncid, 'solar_azimuth', varid), "get saa vaird")
      call check(nf90_put_var(ncid, varid, unpack(reshape(SAA(:,ialong),(/clrm/)),clmask,field4), &
                    start = (/ialong,1/), count = (/1,tm/)), "writing out saa")

      call check(nf90_inq_varid(ncid, 'sensor_azimuth', varid), "get vaa vaird")
      call check(nf90_put_var(ncid, varid, unpack(reshape(VAA(:,ialong),(/clrm/)),clmask,field4), &
                    start = (/ialong,1/), count = (/1,tm/)), "writing out vaa")


    end do
    call check(nf90_inq_varid(ncid, 'rayleigh_depol_ratio', varid), "get depol vaird")
    call check(nf90_put_var(ncid, varid, real(depol)), "writing out depol")

    call check(nf90_inq_varid(ncid, 'ocean_refractive_index', varid), "get ori vaird")
    call check(nf90_put_var(ncid, varid, real(mr)), "writing out ori")


    call check( nf90_close(ncid), "close outfile" )
    if (MAPL_am_I_root()) then
      write(*,*) '<> Wrote Outfile '
    end if


  end subroutine write_outfile


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     DO_OCEAN()
! PURPOSE
!     call the correct VLIDORT Ocean Subroutine
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     Oct 2018 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine DO_OCEAN()
!   WATER PIXELS
!   ------------
    if ( (lower_to_upper(watermodel) .eq. 'CX') ) then

    ! Cox Munk Model
    ! -------------------------------
      ! Call to vlidort vector code
      call VLIDORT_Vector_CX (km, nch, nobs ,nalong, dble(channels), nstreams, plane_parallel, nMom,   &
             nPol, ROT_int, depol, dble(Vtau), dble(Vssa), dble(Vpmom), &
             dble(Vpe), dble(Vze), dble(Vte), &
             (/dble(U10M(c))/), &
             (/dble(V10M(c))/), &
             mr, &
             reshape(dble(SZA(c,:)),(/nobs,nalong/)), &
             reshape(dble(RAA(c,:)),(/nobs,nalong/)), &
             reshape(dble(VZA(c,:)),(/nobs,nalong/)), &
             dble(MISSING),verbose, &
             radiance_VL_int,reflectance_VL_int, Q_int, U_int, Valbedo, BR_Q_int, BR_U_int, ierr)
      !   Check VLIDORT Status
      !   ----------------------------------------------------    
      call mp_check_vlidort(radiance_VL_int,reflectance_VL_int)  

    else
!       Save code for pixels that were not gap filled
!       ---------------------------------------------    
      radiance_VL_int = MISSING
      reflectance_VL_int = MISSING
      Valbedo = MISSING
      if (.not. scalar) then
        Q_int = MISSING
        U_int = MISSING
        BR_Q_int = MISSING
        BR_U_int = MISSING
      end if
      ierr = 0

    end if

  end subroutine DO_OCEAN


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     get_surf_params()
! PURPOSE
!     get surface parameters
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     Oct 2018 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine get_surf_params()
  real,dimension(2)    :: temp_iso, temp_geo, temp_vol
  integer,dimension(1) :: brdf_i, ler_i

  !   Surface Reflectance Parameters
  !    ----------------------------------
      if ( (index(lower_to_upper(landmodel),'RTLS') > 0) ) then
      !  Get BRDF Kernel Weights
      !  ----------------------------
        do ch = 1, nch
          if (landbandmBRDF == 1) then
            kernel_wt(:,ch,nobs) = (/dble(KISO(c,1)),&
                                dble(KGEO(c,1)),&
                                dble(KVOL(c,1))/)
          else
            ! > 2130 uses highest MODIS wavelength band     
            if (channels(ch) >= maxval(landband_cBRDF)) then  
              iband = minloc(abs(landband_cBRDF - channels(ch)), dim = 1)
              kernel_wt(:,ch,nobs) = (/dble(KISO(c,iband)),&
                                dble(KGEO(c,iband)),&
                                dble(KVOL(c,iband))/)
            end if
            
            if (channels(ch) < maxval(landband_cBRDF)) then
              ! nearest neighbor interpolation of kernel weights to wavelength
              if ( (index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0)  .and. (channels(ch) < minval(landband_cBRDF)) ) then
                !  relax kernel weights below 470 to lamberitan 
                !  ----------------------------
                brdf_i = minloc(landband_cBRDF)
                ler_i  = maxloc(landband_cLER)

                temp_iso(1) = LER(c,ler_i(1))
                temp_iso(2) = KISO(c,brdf_i(1))

                temp_geo(1) = 0
                temp_geo(2) = KGEO(c,brdf_i(1))

                temp_vol(1) = 0
                temp_vol(2) = KVOL(c,brdf_i(1))

                kernel_wt(1,ch,nobs) = dble(nn_interp((/landband_cLER(ler_i(1)),landband_cBRDF(brdf_i(1))/),temp_iso,channels(ch),land_missing,MISSING))
                kernel_wt(2,ch,nobs) = dble(nn_interp((/landband_cLER(ler_i(1)),landband_cBRDF(brdf_i(1))/),temp_geo,channels(ch),land_missing,MISSING))
                kernel_wt(3,ch,nobs) = dble(nn_interp((/landband_cLER(ler_i(1)),landband_cBRDF(brdf_i(1))/),temp_vol,channels(ch),land_missing,MISSING))

              else
                kernel_wt(1,ch,nobs) = dble(nn_interp(landband_cBRDF,reshape(KISO(c,:),(/landbandmBRDF/)),channels(ch),land_missing,MISSING))
                kernel_wt(2,ch,nobs) = dble(nn_interp(landband_cBRDF,reshape(KGEO(c,:),(/landbandmBRDF/)),channels(ch),land_missing,MISSING))
                kernel_wt(3,ch,nobs) = dble(nn_interp(landband_cBRDF,reshape(KVOL(c,:),(/landbandmBRDF/)),channels(ch),land_missing,MISSING))    
              end if      
            end if
          end if
          param(:,ch,nobs)     = (/dble(2),dble(1)/)
        end do

        if ( (index(lower_to_upper(landmodel),'BPDF') > 0) ) then
        ! Get BPDF parameters
        ! ----------------------
          if ( (index(lower_to_upper(landname),'MAIGNAN') > 0) ) then
            BPDFparam(1,:,:) = 1.5   ! water refractive index
            BPDFparam(2,:,:) = dble(NDVI(c))
            BPDFparam(3,:,:) = dble(BPDFcoef(c))
          end if
        end if

      else if ( (lower_to_upper(landmodel) == 'LAMBERTIAN') ) then
      !  Get Lambertian Albedo
      !  ------------------------------
        do ch = 1, nch
          if (landbandmLER == 1) then
            Valbedo(nobs,ch,:) = dble(LER(c,1))
          else
            Valbedo(nobs,ch,:) = dble(nn_interp(landband_cLER,reshape(LER(c,:),(/landbandmLER/)),channels(ch),land_missing,MISSING))
          end if
        end do
      end if 

  end subroutine get_surf_params

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     DO_LAND()
! PURPOSE
!     call the correct VLIDORT Land Subroutine
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     Oct 2018 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

  subroutine DO_LAND()

  !     LAND PIXELS
  !     -----------

    if ( (lower_to_upper(landmodel) .eq. 'LAMBERTIAN') ) then
    ! Simple lambertian surface model
    ! -------------------------------
      ! Call to vlidort vector code
      call VLIDORT_Vector_Lambert (km, nch, nobs , nalong, dble(channels), nstreams, plane_parallel, nMom,   &
             nPol, ROT_int, depol, alpha, dble(Vtau), dble(Vssa), dble(Vpmom), &
             dble(Vpe), dble(Vze), dble(Vte), Valbedo,&
             reshape(dble(SZA(c,:)),(/nobs,nalong/)), &
             reshape(dble(RAA(c,:)),(/nobs,nalong/)), &
             reshape(dble(VZA(c,:)),(/nobs,nalong/)), &
             flux_factor, &
             dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, Q_int, U_int, ierr)
      BR_Q_int = 0
      BR_U_int = 0
      !   Check VLIDORT Status
      !   ----------------------------------------------------    
      call mp_check_vlidort(radiance_VL_int,reflectance_VL_int)  

    else if ( (index(lower_to_upper(landmodel),'RTLS') > 0) .and. (index(lower_to_upper(landmodel),'BPDF') .eq. 0) ) then
      if ( ANY(kernel_wt == MISSING) ) then
!       Save code for pixels that were not gap filled
!       ---------------------------------------------    
        radiance_VL_int = MISSING
        reflectance_VL_int = MISSING
        Valbedo = MISSING
        if (.not. scalar) then
          Q_int = MISSING
          U_int = MISSING
          BR_Q_int = MISSING
          BR_U_int = MISSING
        end if
        ierr = 0

      else   
!       MODIS BRDF Surface Model
!       ------------------------------
        ! Call to vlidort vector code
        call VLIDORT_Vector_LandMODIS (km, nch, nobs, nalong, dble(channels), nstreams, plane_parallel, nMom, &
                nPol, ROT_int, depol, dble(Vtau), dble(Vssa), dble(Vpmom), &
                dble(Vpe), dble(Vze), dble(Vte), &
                kernel_wt, param, &
                reshape(dble(SZA(c,:)),(/nobs,nalong/)), &
                reshape(dble(RAA(c,:)),(/nobs,nalong/)), &
                reshape(dble(VZA(c,:)),(/nobs,nalong/)), &
                dble(MISSING),verbose, &
                radiance_VL_int,reflectance_VL_int, Valbedo, Q_int, U_int, BR_Q_int, BR_U_int, ierr )  
    !   Check VLIDORT Status
    !   ----------------------------------------------------    
        call mp_check_vlidort(radiance_VL_int,reflectance_VL_int)  

      end if

    else if ( (index(lower_to_upper(landmodel),'RTLS') > 0) .and. (index(lower_to_upper(landmodel),'BPDF') > 0) ) then  
      if ( ANY(kernel_wt == MISSING) .or. ANY(BPDFparam == land_missing)) then
!       Save code for pixels that were not gap filled
!       ---------------------------------------------    
        radiance_VL_int = MISSING
        reflectance_VL_int = MISSING
        Valbedo = MISSING
        if (.not. scalar) then
          Q_int = MISSING
          U_int = MISSING
        end if
        ierr = 0

!      else   
!!       MODIS BPDF Surface Model
!!       ------------------------------

!        ! Call to vlidort vector code
!        call VLIDORT_Vector_LandMODIS_BPDF (km, nch, nobs, dble(channels), nstreams, plane_parallel, nMom, &
!                nPol, ROT_int, depol, dble(Vtau), dble(Vssa), dble(Vpmom), &
!                dble(Vpe), dble(Vze), dble(Vte), &
!                kernel_wt, param, BPDFparam, &
!                (/dble(SZA(c,:))/), &
!                (/dble(RAA(c,:))/), &
!                (/dble(VZA(c,:))/), &
!                dble(MISSING),verbose, &
!                radiance_VL_int,reflectance_VL_int, Valbedo, Q_int, U_int, BR_Q_int, BR_U_int, ierr )  

!    !   Check VLIDORT Status
!    !   ----------------------------------------------------    
!        call mp_check_vlidort(radiance_VL_int,reflectance_VL_int)  

      end if

    end if   

  end subroutine DO_LAND


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     nn_interp
! PURPOSE
!     nearest neighbor interpolation
! INPUT
!     x     : x-values
!     y     : y-values
!     xint  : x-value to interpolate to
! OUTPUT
!     nn_interp  : interpolated y-value
!  HISTORY
!     10 June 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

  function nn_interp(x,y,xint,ymissing,missing)
    real,intent(in),dimension(:)    :: x, y
    real,intent(in)                 :: xint
    real                            :: nn_interp
    integer                         :: below, above
    real                            :: top, bottom
    real                            :: ymissing,missing

    if (xint .GE. maxval(x)) then
        below = maxloc(x, dim = 1)
        nn_interp = y(below)
    else if (xint .LT. minval(x)) then
        above = minloc(x, dim = 1)
        nn_interp = y(above)
    else
      above = minloc(abs(xint - x), dim = 1, mask = (xint - x) .LT. 0)
      below = minloc(abs(xint - x), dim = 1, mask = (xint - x) .GE. 0)

      if (.not. ANY((/y(above),y(below)/) == ymissing)) then
        top = y(above) - y(below)
        bottom = x(above) - x(below)
        nn_interp = y(below) + (xint-x(below)) * top / bottom
      else
        nn_interp  = missing
      end if
    end if
  end function nn_interp

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
!     read_land
! PURPOSE
!     allocates land variables and reads in data
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     28 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine read_land()
    if (MAPL_am_I_root()) then
      call readvar1D("FRLAND", INV_file, READER1D)
      FRLAND = pack(READER1D,clmask)
      write(*,*) '<> Read fraction land data to shared memory'
    end if       

    call MAPL_SyncSharedMemory(rc=ierr)
  end subroutine read_land

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     read_sza
! PURPOSE
!     allocates sza variable and reads in data
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     28 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine read_sza()

    allocate (SOLAR_ZENITH(nalong,tm))

    call readvar2D("sza_ss", ANG_file, SOLAR_ZENITH)
 
  end subroutine read_sza


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
!     mp_check_vlidort
! PURPOSE
!     exits if VLIDORT does not give good output
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine mp_check_vlidort(radiance,reflectance)
    real*8, dimension(:,:,:)    :: radiance, reflectance


    if (ANY(radiance == MISSING) .or. ANY(reflectance == MISSING)) then
      write(*,*) 'VLIDORT returned a missing value'
      write(*,*) 'Exiting......'
      write(*,*) 'Pixel is ',c
      call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
    else if (ierr > 0) then
      write(*,*) 'VLIDORT returned rc code ', ierr
      write(*,*) 'Exiting......'
      write(*,*) 'Pixel is ',c
      call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
    end if 

  end subroutine mp_check_vlidort

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     read_aer_Nv
! PURPOSE
!     read in all the aerosol variables
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine read_aer_Nv()
    if (MAPL_am_I_root()) then
      call readvar2D("DU001", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,DU001)

      call readvar2D("DU002", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,DU002)

      call readvar2D("DU003", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,DU003)

      call readvar2D("DU004", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,DU004)

      call readvar2D("DU005", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,DU005)

      call readvar2D("SS001", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,SS001)

      call readvar2D("SS002", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,SS002)

      call readvar2D("SS003", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,SS003)

      call readvar2D("SS004", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,SS004)

      call readvar2D("SS005", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,SS005)

      call readvar2D("BCPHOBIC", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,BCPHOBIC)

      call readvar2D("BCPHILIC", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,BCPHILIC)

      call readvar2D("OCPHOBIC", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,OCPHOBIC)

      call readvar2D("OCPHILIC", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,OCPHILIC)

      call readvar2D("SO4", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,SO4)

      call readvar2D("RH", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,RH)

      call readvar2D("DELP", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,DELP)

      call readvar2D("AIRDENS", AER_file, READER2D) 
      call reduceProfile(READER2D,clmask,AIRDENS)

      call readvar1D("PS", AER_file, READER1D) 
      PS = pack(READER1D,clmask)

      write(*,*) '<> Read aerosol data to shared memory'
    end if

    call MAPL_SyncSharedMemory(rc=ierr)    

  end subroutine read_aer_Nv


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
    real,intent(in),dimension(:,:)               :: var
    logical,intent(in),dimension(:)              :: mask
    real,intent(inout),dimension(:,:)            :: reducedProfile

    integer                                        :: tm, km, k

    km = size(var,1)
    tm = size(var,2)

    do k = 1, km
      reducedProfile(:,k) = pack(reshape(var(k,:),(/tm/)),mask)
    end do

  end subroutine reduceProfile


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     read_surf_land
! PURPOSE
!     read in all the surface reflectance variables
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!     Jul 2015 P. Castellanos - added OMI LER option
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine read_surf_land()
    if (MAPL_am_I_root()) then
      ! RTLS Kernel
      if ((index(lower_to_upper(landmodel),'RTLS') > 0)) then
        call read_RTLS()
        write(*,*) '<> Read BRDF data to shared memory'

        if ((index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0)) then
          call read_LER()   
          write(*,*) '<> Read LER data to shared memory'
        end if
      end if

      ! BPDF
      if ( (index(lower_to_upper(landmodel),'BPDF') > 0) ) then
        call read_BPDF()
        write(*,*) '<> Read BPDF data to shared memory'
      end if

      ! Lambertian
      if (lower_to_upper(landmodel) == 'LAMBERTIAN') then
        call read_LER()     
        write(*,*) '<> Read LER data to shared memory'
      end if 
    end if
    call MAPL_SyncSharedMemory(rc=ierr) 
  end subroutine read_surf_land


  subroutine read_RTLS()
    character(len=100)                 :: sds

    do ch = 1,landbandmBRDF
      if ( landband_cBRDF(ch) < 1000 ) then
        write(sds,'(A4,I3)') "Riso",int(landband_cBRDF(ch))
      else
        write(sds,'(A4,I4)') "Riso",int(landband_cBRDF(ch))
      end if

      call readvar1D(trim(sds), BRDF_file, READER1D) 
      KISO(:,ch) = pack(READER1D,clmask)

      if ( landband_cBRDF(ch) < 1000 ) then
        write(sds,'(A4,I3)') "Rgeo",int(landband_cBRDF(ch))
      else
        write(sds,'(A4,I4)') "Rgeo",int(landband_cBRDF(ch))
      end if
      
      call readvar1D(trim(sds), BRDF_file, READER1D) 
      KGEO(:,ch) = pack(READER1D,clmask)

      if ( landband_cBRDF(ch) < 1000 ) then
        write(sds,'(A4,I3)') "Rvol",int(landband_cBRDF(ch))
      else
        write(sds,'(A4,I4)') "Rvol",int(landband_cBRDF(ch))
      end if

      call readvar1D(trim(sds), BRDF_file, READER1D) 
      KVOL(:,ch) = pack(READER1D,clmask)      
    end do
    
  end subroutine read_RTLS

  subroutine read_BPDF()

    call readvar1D('NDVI', NDVI_file, READER1D) 
    NDVI = pack(READER1D,clmask)
    call readvar1D('BPDFcoef', BPDF_file, READER1D) 
    BPDFcoef= pack(READER1D,clmask)
  end subroutine read_BPDF

  subroutine read_LER()
    real, pointer                      :: temp(:,:)
    character(len=100)                 :: sds

    allocate(temp(1,tm))
    do ch = 1,landbandmLER
      if ( landband_cLER(ch) < 1000 ) then
        write(sds,'(A6,I3)') "SRFLER",int(landband_cLER(ch))
      end if
      
      call readvar2D(trim(sds), LER_file, temp) 
      LER(:,ch) = pack(reshape(temp(1,:),(/tm/)),clmask)
    end do
    deallocate(temp)
  end subroutine read_LER


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     read_angles
! PURPOSE
!     read in all the sensor and solar position angles
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine read_angles()
    real, allocatable          :: saa_(:,:),angle(:,:)
    integer                    :: i,n

    if (MAPL_am_I_root()) then  
      allocate (saa_(clrm,nalong))
      allocate (angle(nalong,tm))

      call readvar2D("sza_ss", ANG_file, angle)
      call reduceProfile(angle,clmask,SZA)

      call readvar2D("vza_ss", ANG_file, angle)
      call reduceProfile(angle,clmask,VZA)

      call readvar2D("saa_ss", ANG_file, angle)
      call reduceProfile(angle,clmask,SAA)

      call readvar2D("vaa_ss", ANG_file, angle)
      call reduceProfile(angle,clmask,VAA)

     
      ! define according to photon travel direction
      saa_ = SAA + 180.0
      do i = 1, clrm
        do n = 1,nalong
          if (saa_(i,n) >= 360.0) then
            saa_(i,n) = saa_(i,n) - 360.0
          end if
        end do
      end do

      RAA = VAA - saa_
      do i = 1, clrm
        do n = 1,nalong
          if (RAA(i,n) < 0) then
            RAA(i,n) = RAA(i,n) + 360.0
          end if
        end do
      end do
      deallocate (saa_)
      write(*,*) '<> Read angle data to shared memory' 

    end if
    call MAPL_SyncSharedMemory(rc=ierr)    


  end subroutine read_angles  

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     read_wind
! PURPOSE
!     read in wind speed
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine read_wind()

    if (MAPL_am_I_root()) then
      call readvar1D("U10M",   MET_file, READER1D) 
      U10M = pack(READER1D,clmask)

      call readvar1D("V10M",  MET_file, READER1D) 
      V10M = pack(READER1D,clmask)
      write(*,*) '<> Read wind data to shared memory' 
    end if

    call MAPL_SyncSharedMemory(rc=ierr)    
    
  end subroutine read_wind


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     allocate_shared
! PURPOSE
!     allocates all the shared memory arrays that I need
! INPUT
!     clrm   :: number of pixels after cloud filtering
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine allocate_shared(clrm)
    integer, intent(in)      :: clrm

    if (MAPL_am_I_root()) then
      write(*,*) '<> Allocating shared memory variables',tm,km, clrm
    end if
    call MAPL_AllocNodeArray(READER2D,(/km,tm/),rc=ierr)
    call MAPL_AllocNodeArray(READER1D,(/tm/),rc=ierr)
    call MAPL_AllocNodeArray(AIRDENS,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(RH,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DELP,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(PS,(/clrm/),rc=ierr)
    call MAPL_AllocNodeArray(DU001,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU002,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU003,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU004,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU005,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS001,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS002,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS003,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS004,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS005,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(BCPHOBIC,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(BCPHILIC,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(OCPHOBIC,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(OCPHILIC,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SO4,(/clrm,km/),rc=ierr)
    if ((index(lower_to_upper(landmodel),'RTLS') > 0)) then
      call MAPL_AllocNodeArray(KISO,(/clrm,landbandmBRDF/),rc=ierr)
      call MAPL_AllocNodeArray(KVOL,(/clrm,landbandmBRDF/),rc=ierr)
      call MAPL_AllocNodeArray(KGEO,(/clrm,landbandmBRDF/),rc=ierr)
    end if
    if ((lower_to_upper(landmodel) == 'LAMBERTIAN') .or. (index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0)) then
      call MAPL_AllocNodeArray(LER,(/clrm,landbandmLER/),rc=ierr)
    end if 
    if ( (index(lower_to_upper(landmodel),'BPDF') > 0) ) then
      if ( (index(lower_to_upper(landname),'MAIGNAN') > 0) ) then
        call MAPL_AllocNodeArray(NDVI,(/clrm/),rc=ierr)
        call MAPL_AllocNodeArray(BPDFcoef,(/clrm/),rc=ierr)
      end if
    end if

    call MAPL_AllocNodeArray(SZA,(/clrm,nalong/),rc=ierr)
    call MAPL_AllocNodeArray(VZA,(/clrm,nalong/),rc=ierr)
    call MAPL_AllocNodeArray(RAA,(/clrm,nalong/),rc=ierr)
    call MAPL_AllocNodeArray(SAA,(/clrm,nalong/),rc=ierr)
    call MAPL_AllocNodeArray(VAA,(/clrm,nalong/),rc=ierr) 

    call MAPL_AllocNodeArray(FRLAND,(/clrm/),rc=ierr) 

    call MAPL_AllocNodeArray(U10M,(/clrm/),rc=ierr) 
    call MAPL_AllocNodeArray(V10M,(/clrm/),rc=ierr) 

    call MAPL_AllocNodeArray(TAU,(/clrm/),rc=ierr)
    call MAPL_AllocNodeArray(SSA,(/clrm/),rc=ierr)
    call MAPL_AllocNodeArray(G,(/clrm/),rc=ierr)

    call MAPL_AllocNodeArray(ROT,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(PE,(/clrm,km+1/),rc=ierr)
    call MAPL_AllocNodeArray(ZE,(/clrm,km+1/),rc=ierr)
    call MAPL_AllocNodeArray(TE,(/clrm,km+1/),rc=ierr)
    
    if (.not. scalar) then
      call MAPL_AllocNodeArray(Q,(/clrm,nalong/),rc=ierr)
      call MAPL_AllocNodeArray(U,(/clrm,nalong/),rc=ierr)
      call MAPL_AllocNodeArray(BR_Q,(/clrm,nalong/),rc=ierr)
      call MAPL_AllocNodeArray(BR_U,(/clrm,nalong/),rc=ierr)

    end if

    call MAPL_AllocNodeArray(ALBEDO,(/clrm,nalong/),rc=ierr)
    call MAPL_AllocNodeArray(radiance_VL,(/clrm,nalong/),rc=ierr)
    call MAPL_AllocNodeArray(reflectance_VL,(/clrm,nalong/),rc=ierr)

  end subroutine allocate_shared

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     allocate_multinode
! PURPOSE
!     allocates all the unshared arrays that I need on the nodes that are not root
! INPUT
!     clrm   :: number of pixels this processor needs to work on
! OUTPUT
!     none
!  HISTORY
!     
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine allocate_multinode(clrm)
    integer, intent(in)    :: clrm

    if (clrm > 0) then
      allocate (AIRDENS(clrm,km))
      allocate (RH(clrm,km))
      allocate (DELP(clrm,km))
      allocate (PS(clrm))
      allocate (DU001(clrm,km))
      allocate (DU002(clrm,km))
      allocate (DU003(clrm,km))
      allocate (DU004(clrm,km))
      allocate (DU005(clrm,km))
      allocate (SS001(clrm,km))
      allocate (SS002(clrm,km))
      allocate (SS003(clrm,km))
      allocate (SS004(clrm,km))
      allocate (SS005(clrm,km))
      allocate (BCPHOBIC(clrm,km))
      allocate (BCPHILIC(clrm,km))
      allocate (OCPHOBIC(clrm,km))
      allocate (OCPHILIC(clrm,km))
      allocate (SO4(clrm,km))
      if ((index(lower_to_upper(landmodel),'RTLS') > 0)) then
        allocate (KISO(clrm,landbandmBRDF))
        allocate (KVOL(clrm,landbandmBRDF))
        allocate (KGEO(clrm,landbandmBRDF))
      end if
      if ((lower_to_upper(landmodel) == 'LAMBERTIAN') .or. (index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0)) then
        allocate (LER(clrm,landbandmLER))
      end if 
      if ( (index(lower_to_upper(landmodel),'BPDF') > 0) ) then
        if ( (index(lower_to_upper(landname),'MAIGNAN') > 0) ) then  
          allocate (NDVI(clrm))    
          allocate (BPDFcoef(clrm))
        end if
      end if

      allocate (SZA(clrm,nalong))
      allocate (VZA(clrm,nalong))
      allocate (RAA(clrm,nalong))
      allocate (SAA(clrm,nalong))
      allocate (VAA(clrm,nalong))

      allocate (FRLAND(clrm))

      allocate (U10M(clrm))
      allocate (V10M(clrm))

      allocate (TAU(clrm))
      allocate (SSA(clrm))
      allocate (G(clrm))

      allocate (ROT(clrm,km))      
      allocate (PE(clrm,km+1))
      allocate (ZE(clrm,km+1))
      allocate (TE(clrm,km+1))
      
      if (.not. scalar) then
        allocate (Q(clrm,nalong))
        allocate (U(clrm,nalong))
        allocate (BR_Q(clrm,nalong))
        allocate (BR_U(clrm,nalong))
      end if

      allocate (ALBEDO(clrm,nalong))
      allocate (radiance_VL(clrm,nalong))
      allocate (reflectance_VL(clrm,nalong))
    end if
  end subroutine allocate_multinode


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
    allocate (Vpe(km+1,nobs))
    allocate (Vze(km+1,nobs))
    allocate (Vte(km+1,nobs))
    allocate (Vqm(km,nq,nobs))
    allocate (Vtau(km,nch,nobs))
    allocate (Vssa(km,nch,nobs))
    allocate (Vg(km,nch,nobs))
    allocate (Valbedo(nobs,nch,nalong))

    allocate (radiance_VL_int(nobs,nch,nalong))
    allocate (reflectance_VL_int(nobs, nch,nalong))    

    if ((index(lower_to_upper(landmodel),'RTLS') > 0) ) then
      allocate (kernel_wt(nkernel,nch,nobs))
      allocate (param(nparam,nch,nobs))
    end if

    if ( (index(lower_to_upper(landmodel),'BPDF') > 0) ) then  
      if ( (index(lower_to_upper(landname),'MAIGNAN') > 0) ) then    
        allocate (BPDFparam(3,nch,nobs)) 
      end if
    end if
    

    allocate (ROT_int(km,nobs,nch))
    allocate (depol(nch))
    allocate (alpha(km,nobs,nch))
    alpha = 0.0
    allocate (Vpmom(km,nch,nobs,nMom,nPol))
    allocate (flux_factor(nch,nobs))
    flux_factor = 1.0

    if (.not. scalar) then      
      allocate (Q_int(nobs, nch, nalong))
      allocate (U_int(nobs, nch, nalong))
      allocate (BR_Q_int(nobs, nch, nalong))
      allocate (BR_U_int(nobs, nch, nalong))
    end if

  ! Needed for reading
  ! ----------------------
    allocate (nclr(npet-1)) 
    nclr = 0

    allocate(clmask(tm))
    clmask = .False.

  end subroutine allocate_unshared


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     write_verbose
! PURPOSE
!     manages verbose messaging
! INPUT
!     msg  : text to output
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

  subroutine write_verbose(msg)
    character(len=*), intent(in)      :: msg

    if (verbose==1 .and. verbose_mp==1) then
      call mp_write_msg(msg) 
    else if (verbose==1) then
      write(*,*) trim(msg)
    end if
  end subroutine write_verbose

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     mp_write_msg
! PURPOSE
!     writes stdio messages from parallel processors in order
! INPUT
!     msg  : text to output
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine mp_write_msg(msg)
    character(len=100), intent(in)    :: msg
    character(len=100)                :: msgr
    integer                           :: pp

    if ( MAPL_am_I_root() ) then
      do pp = 1,npet-1
        call mpi_recv(msgr, 100, MPI_CHARACTER, pp, 1, MPI_COMM_WORLD, status_mpi, ierr)
        write(*,*) trim(msgr)
      end do
      write(*,*) trim(msg)
    else
      call mpi_send(msg, 100, MPI_CHARACTER, 0, 1, MPI_COMM_WORLD, ierr)
    end if 

  end subroutine mp_write_msg

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     create_outfile
! PURPOSE
!     creates the output netcdf4 file 
! INPUT
!     date, time  : from RC file
! OUTPUT
!     None
!  HISTORY
!     6 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine create_outfile(date, time)
    character(len=*) ,intent(in)       :: date, time

    integer                            :: radVarID, refVarID
    integer                            :: qVarID, uVarID, albVarID, brUVarID, brQVarID      
    
    integer                            :: ncid
    integer                            :: timeDimID, xDimID, yDimID, levDimID, lsDimID, leveDimID 
    integer                            :: nvzaDimID, chDimID     
    integer                            :: lonVarID, latVarID
    integer                            :: timeVarID, levVarID, xVarID, yVarID
    integer                            :: vzaVarID, vaaVarID, szaVarID, saaVarID
    integer                            :: rotVarID, depolVarID, oriVarID, zeVarID
    integer                            :: ch

    real*8,allocatable,dimension(:)    :: lon, lat, sza, vza, raa
    real*8,allocatable,dimension(:)    :: scantime, x, y, lev

    character(len=2000)                :: comment

!                                MAIN LevelC2 OUT_FILE 
!                                ----------------------
    ! Open File
    call check(nf90_create(OUT_file, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // OUT_file)

    ! Create dimensions
    call check(nf90_def_dim(ncid, "channel", 1, chDimID), "creating channel dimension") 
    call check(nf90_def_dim(ncid, "view_angles", nalong, nvzaDimID), "creating view_angles dimension") 
    call check(nf90_def_dim(ncid, "time", tm, timeDimID), "creating time dimension")
    call check(nf90_def_dim(ncid, "lev", km, levDimID), "creating lev dimension") !km
    call check(nf90_def_dim(ncid, "leve", km+1, leveDimID), "creating leve dimension") !km+1
    call check(nf90_def_dim(ncid, "x", 1, xDimID), "creating x dimension") 
    call check(nf90_def_dim(ncid, "y", 1, yDimID), "creating y dimension") 

    ! Global Attributes
    write(comment,'(A)') "VLIDORT Simulation of GEOS-5 multiangle polarized reflectance"
    call check(nf90_put_att(ncid,NF90_GLOBAL,'title',trim(comment)),"title attr")

    write(comment,'(A)') 'NASA/Goddard Space Flight Center'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'institution',trim(comment)),"institution attr")

    write(comment,'(A)') 'Global Model and Assimilation Office'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'source',trim(comment)),"source attr")

    write(comment,'(A)') "VLIDORT simulation run on sampled GEOS-5"
    call check(nf90_put_att(ncid,NF90_GLOBAL,'history',trim(comment)),"history attr")

    write(comment,'(A)') "n/a"
    call check(nf90_put_att(ncid,NF90_GLOBAL,'references',trim(comment)),"reference attr")  

    write(comment,'(A)') "Patricia Castellanos <patricia.castellanos@nasa.gov>"
    call check(nf90_put_att(ncid,NF90_GLOBAL,'contact',trim(comment)),"contact attr")      

    write(comment,'(A)') "CF"
    call check(nf90_put_att(ncid,NF90_GLOBAL,'conventions',trim(comment)),"conventions attr")  

    ! Define Variables
!                                     Dimensions
!                                     ----------    
    call check(nf90_def_var(ncid,'trjLon',nf90_float,(/timeDimID/),lonVarID),"create lon var")
    call check(nf90_def_var(ncid,'trjLat',nf90_float,(/timeDimID/),latVarID),"create lat var")
    call check(nf90_def_var(ncid,'time',nf90_int,(/timeDimID/),timeVarID),"create time var")
    call check(nf90_def_var(ncid,'lev',nf90_float,(/levDimID/),levVarID),"create lev var")
    call check(nf90_def_var(ncid,'x',nf90_float,(/xDimID/),xVarID),"create x var")
    call check(nf90_def_var(ncid,'y',nf90_float,(/yDimID/),yVarID),"create y var")

!                                     Angles
!                                     ------
    call check(nf90_def_var(ncid,'sensor_zenith',nf90_float,(/nvzaDimID,timeDimID/),vzaVarID),"create vza var")
    call check(nf90_def_var(ncid,'sensor_azimuth',nf90_float,(/nvzaDimID,timeDimID/),vaaVarID),"create vaa var")
    call check(nf90_def_var(ncid,'solar_zenith',nf90_float,(/nvzaDimID,timeDimID/),szaVarID),"create sza var")
    call check(nf90_def_var(ncid,'solar_azimuth',nf90_float,(/nvzaDimID,timeDimID/),saaVarID),"create saa var")

!                                     Data
!                                     ----
    call check(nf90_def_var(ncid, 'toa_reflectance',nf90_float,(/nvzaDimID,timeDimID/),refVarID),"create reflectance var")
    call check(nf90_def_var(ncid, 'I',nf90_float,(/nvzaDimID,timeDimID/),radVarID),"create radiance var")              
    if (.not. scalar) then
      call check(nf90_def_var(ncid, 'Q',nf90_float,(/nvzaDimID,timeDimID/),qVarID),"create Q var")
      call check(nf90_def_var(ncid, 'U',nf90_float,(/nvzaDimID,timeDimID/),uVarID),"create U var")
    end if        

    call check(nf90_def_var(ncid, 'surf_reflectance',nf90_float,(/nvzaDimID,timeDimID/),albVarID),"create albedo var")
    if (.not. scalar) then
      call check(nf90_def_var(ncid, 'surf_reflectance_Q',nf90_float,(/nvzaDimID,timeDimID/),brQVarID),"create Q var")
      call check(nf90_def_var(ncid, 'surf_reflectance_U',nf90_float,(/nvzaDimID,timeDimID/),brUVarID),"create U var")
    end if        

    call check(nf90_def_var(ncid,'ROT',nf90_float,(/levDimID,timeDimID/),rotVarID),"create rot var")
    call check(nf90_def_var(ncid,'rayleigh_depol_ratio',nf90_float,(/chDimID/),depolVarID),"create depol var")
    call check(nf90_def_var(ncid,'ocean_refractive_index',nf90_float,(/chDimID/),oriVarID),"create ori var")

    call check(nf90_def_var(ncid,'ZE',nf90_float,(/leveDimID,timeDimID/),zeVarID),"create ze var")

    ! Variable Attributes
!                                          Reflectance TOA & Surface
!                                          -------------------------  
    ch = 1
    write(comment,'(F10.2,A)') channels(ch), ' nm TOA Reflectance'
    call check(nf90_put_att(ncid,refVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
    write(comment,'(F10.2,A)') channels(ch), ' nm reflectance at the top of the atmosphere'
    call check(nf90_put_att(ncid,refVarID,'long_name',trim(adjustl(comment))),"long_name attr")
    call check(nf90_put_att(ncid,refVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,refVarID,'units','None'),"units attr")
    call check(nf90_put_att(ncid,refVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

    write(comment,'(F10.2,A)') channels(ch), ' nm TOA I'
    call check(nf90_put_att(ncid,radVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
    write(comment,'(F10.2,A)') channels(ch), ' nm intensity at the top of the atmosphere'
    call check(nf90_put_att(ncid,radVarID,'long_name',trim(adjustl(comment))),"long_name attr")
    call check(nf90_put_att(ncid,radVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,radVarID,'units','W m-2 sr-1 nm-1'),"units attr")
    call check(nf90_put_att(ncid,radVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

    if (.not. scalar) then
      write(comment,'(F10.2,A)') channels(ch), ' nm TOA Q'
      call check(nf90_put_att(ncid,qVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm TOA Q Component of the Stokes Vector'
      call check(nf90_put_att(ncid,qVarID,'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,qVarID,'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,qVarID,'units','W m-2 sr-1 nm-1'),"units attr")
      call check(nf90_put_att(ncid,qVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

      write(comment,'(F10.2,A)') channels(ch), ' nm TOA U'
      call check(nf90_put_att(ncid,uVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm TOA U Component of the Stokes Vector'
      call check(nf90_put_att(ncid,uVarID,'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,uVarID,'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,uVarID,'units','W m-2 sr-1 nm-1'),"units attr")
      call check(nf90_put_att(ncid,uVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")
    end if

    write(comment,'(F10.2,A)') channels(ch), ' nm Surface Reflectance'
    call check(nf90_put_att(ncid,albVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
    write(comment,'(F10.2,A)') channels(ch), ' nm Bi-Directional Surface Reflectance'
    call check(nf90_put_att(ncid,albVarID,'long_name',trim(adjustl(comment))),"long_name attr")
    call check(nf90_put_att(ncid,albVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,albVarID,'units','None'),"units attr")
    call check(nf90_put_att(ncid,albVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

    if (.not. scalar) then
      write(comment,'(F10.2,A)') channels(ch), ' nm Surface Reflectance Q'
      call check(nf90_put_att(ncid,brQVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm Bi-Directional Surface Reflectance Q'
      call check(nf90_put_att(ncid,brQVarID,'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,brQVarID,'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,brQVarID,'units','None'),"units attr")
      call check(nf90_put_att(ncid,brQVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

      write(comment,'(F10.2,A)') channels(ch), ' nm Surface Reflectance U'
      call check(nf90_put_att(ncid,brUVarID,'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm Bi-Directional Surface Reflectance U'
      call check(nf90_put_att(ncid,brUVarID,'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,brUVarID,'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,brUVarID,'units','None'),"units attr")
      call check(nf90_put_att(ncid,brUVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

    end if


!                                          scanTime
!                                          -------  
    call check(nf90_put_att(ncid,timeVarID,'long_name','Time'),"long_name attr")
    write(comment,*) 'seconds since ',date(1:4),'-',date(5:6),'-',date(7:8),' ',time(1:2),':',time(3:4),':',time(5:6)
    call check(nf90_put_att(ncid,timeVarID,'units',trim(adjustl(comment))),"units attr")



!                                          EW, NS, LEV, TIME
!                                          -----------------------  
    call check(nf90_put_att(ncid,xVarID,'long_name','Fake Longitude for GrADS Compatibility'),"long_name attr")
    call check(nf90_put_att(ncid,xVarID,'units','degrees_east'),"units attr")
    call check(nf90_put_att(ncid,yVarID,'long_name','Fake Latitude for GrADS Compatibility'),"long_name attr")
    call check(nf90_put_att(ncid,yVarID,'units','degrees_north'),"units attr")   
    call check(nf90_put_att(ncid,levVarID,'long_name','Vertical Level'),"long_name attr")
    call check(nf90_put_att(ncid,levVarID,'units','layer'),"units attr")
    call check(nf90_put_att(ncid,levVarID,'positive','down'),"positive attr")
    call check(nf90_put_att(ncid,levVarID,'axis','z'),"axis attr")

!                                          clon & clat
!                                          -------  
    call check(nf90_put_att(ncid,lonVarID,'long_name','Trajectory Longitude'),"long_name attr")
    call check(nf90_put_att(ncid,lonVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,lonVarID,'units','degrees_east'),"units attr")
    call check(nf90_put_att(ncid,latVarID,'long_name','Trajectory Latitude'),"long_name attr")
    call check(nf90_put_att(ncid,latVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,latVarID,'units','degrees_north'),"units attr")     

!                                          geometry
!                                          -------  
    call check(nf90_put_att(ncid,vzaVarID,'long_name','sensor viewing zenith angle (VZA)'),"long_name attr")
    call check(nf90_put_att(ncid,vzaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,vzaVarID,'units','degrees'),"units attr")

    call check(nf90_put_att(ncid,vaaVarID,'long_name','sensor viewing azimuth angle (VAA)'),"long_name attr")
    call check(nf90_put_att(ncid,vaaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,vaaVarID,'units','degrees clockwise from North (0-360)'),"units attr")

    call check(nf90_put_att(ncid,szaVarID,'long_name','solar viewing zenith angle (SZA)'),"long_name attr")
    call check(nf90_put_att(ncid,szaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,szaVarID,'units','degrees'),"units attr")

    call check(nf90_put_att(ncid,saaVarID,'long_name','solar viewing azimuth angle (SAA)'),"long_name attr")
    call check(nf90_put_att(ncid,saaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,saaVarID,'units','degrees clockwise from North (0-360)'),"units attr")

!                                          ROT
!                                          -------  
    write(comment,'(F10.2,A)') channels(ch), ' nm Rayleigh Optical Thickness'
    call check(nf90_put_att(ncid,rotVarID,'long_name',trim(adjustl(comment))),"standard_name attr")
    call check(nf90_put_att(ncid,rotVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,rotVarID,'units','None'),"units attr")
    call check(nf90_put_att(ncid,rotVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

    write(comment,'(F10.2,A)') channels(ch), ' nm Rayleigh Depolrization Ratio'
    call check(nf90_put_att(ncid,depolVarID,'long_name',trim(adjustl(comment))),"standard_name attr")
    call check(nf90_put_att(ncid,depolVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,depolVarID,'units','None'),"units attr")
    call check(nf90_put_att(ncid,depolVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

    write(comment,'(F10.2,A)') channels(ch), ' nm ocean refreactive index'
    call check(nf90_put_att(ncid,oriVarID,'long_name',trim(adjustl(comment))),"standard_name attr")
    call check(nf90_put_att(ncid,oriVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,oriVarID,'units','None'),"units attr")
    call check(nf90_put_att(ncid,oriVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

!                                          ZE
!                                          -------
    write(comment,'(A)') 'Layer Edge Altitude Above Surface'
    call check(nf90_put_att(ncid,zeVarID,'long_name',trim(adjustl(comment))),"standard_name attr")
    call check(nf90_put_att(ncid,zeVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,zeVarID,'units','m'),"units attr")
    call check(nf90_put_att(ncid,zeVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")



    !Leave define mode
    call check(nf90_enddef(ncid),"leaving define mode")

    ! write out ew, ns, lev, time, clon, clat, & scantime
    allocate (scantime(tm))
    allocate (lon(tm))
    allocate (lat(tm))
    allocate (x(1))
    allocate (y(1))
    allocate (lev(km))


    call readvar1D("time", INV_file, scantime)
    call check(nf90_put_var(ncid,timeVarID,int(scantime)), "writing out time")

    call readvar1D("trjLon", INV_file, lon)
    call check(nf90_put_var(ncid,lonVarID,lon), "writing out lon")

    call readvar1D("trjLat", INV_file, lat)
    call check(nf90_put_var(ncid,latVarID,lat), "writing out lat")

    call readvar1D("lev", AER_file, lev)
    call check(nf90_put_var(ncid,levVarID,lev), "writing out lev")

    call readvar1D("x", INV_file, x)
    call check(nf90_put_var(ncid,xVarID,x), "writing out x")

    call readvar1D("y", INV_file, y)
    call check(nf90_put_var(ncid,yVarID,y), "writing out y")   

    call check( nf90_close(ncid), "close outfile" )

  end subroutine create_outfile


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    cal_qm
! PURPOSE
!     calculates qm needed by VLIDORT
! INPUT
!     var             : aerosol mixin ratio
!     q               : species index
!     n               : obs index
!     i, j            : ew, ns index
! OUTPUT
!     None
!  HISTORY
!     6 May P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine calc_qm()
    integer                                :: k

    do k = 1, km
      Vqm(k,1,nobs) = DU001(c,k)*DELP(c,k)/grav
      Vqm(k,2,nobs) = DU002(c,k)*DELP(c,k)/grav
      Vqm(k,3,nobs) = DU003(c,k)*DELP(c,k)/grav
      Vqm(k,4,nobs) = DU004(c,k)*DELP(c,k)/grav
      Vqm(k,5,nobs) = DU005(c,k)*DELP(c,k)/grav
      Vqm(k,6,nobs) = SS001(c,k)*DELP(c,k)/grav
      Vqm(k,7,nobs) = SS002(c,k)*DELP(c,k)/grav
      Vqm(k,8,nobs) = SS003(c,k)*DELP(c,k)/grav
      Vqm(k,9,nobs) = SS004(c,k)*DELP(c,k)/grav
      Vqm(k,10,nobs) = SS005(c,k)*DELP(c,k)/grav
      Vqm(k,11,nobs) = BCPHOBIC(c,k)*DELP(c,k)/grav
      Vqm(k,12,nobs) = BCPHILIC(c,k)*DELP(c,k)/grav
      Vqm(k,13,nobs) = OCPHOBIC(c,k)*DELP(c,k)/grav
      Vqm(k,14,nobs) = OCPHILIC(c,k)*DELP(c,k)/grav
      Vqm(k,15,nobs) = SO4(c,k)*DELP(c,k)/grav
    end do
   end subroutine calc_qm

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    strarr_2_chararr
! PURPOSE
!     converts a 1-D array of strings to a 2-D array of characters
! INPUT
!     strings  : 1D string array
!     ns       : size of string array
!     maxlen   : max length of strings
! OUTPUT
!     chars
!  HISTORY
!     6 MayP. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine strarr_2_chararr(strings, ns, maxlen, chars)
    character(len=maxlen), intent(in)     :: strings(ns)
    integer, intent(in)                   :: ns
    integer, intent(in)                   :: maxlen
    character, intent(inout)              :: chars(ns,maxlen)

    integer                               :: s, c

    do s = 1,ns         
      do c = 1,maxlen
        chars(s,c) = strings(s)(c:c)
      end do
    end do
  end subroutine strarr_2_chararr

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
  subroutine get_config(configfile)
    character(len=*),intent(in)      :: configfile
    character(len=256)               :: file

    cf = ESMF_ConfigCreate()
    call ESMF_ConfigLoadFile(cf, fileName=trim(configfile), __RC__)

    ! Read in variables
    ! -----------------------
    ! General
    call ESMF_ConfigGetAttribute(cf, date, label = 'DATE:',__RC__)
    call ESMF_ConfigGetAttribute(cf, time, label = 'TIME:',__RC__)
    call ESMF_ConfigGetAttribute(cf, instname, label = 'INSTNAME:',__RC__)
    call ESMF_ConfigGetAttribute(cf, scalar, label = 'SCALAR:',default=.TRUE.)
    call ESMF_ConfigGetAttribute(cf, plane_parallel, label = 'PLANE_PARALLEL:',default=.TRUE.)
    call ESMF_ConfigGetAttribute(cf, aerosol_free, label = 'AEROSOL_FREE:',default=.TRUE.)
    call ESMF_ConfigGetAttribute(cf, nstreams, label = 'NSTREAMS:',default=12)
    call ESMF_ConfigGetAttribute(cf, szamax, label = 'SZAMAX:',default=80.0)
    call ESMF_ConfigGetAttribute(cf, rcfile, label = 'RCFILE:',default='Aod_EOS.rc')
    
    ! Land Surface
    call ESMF_ConfigGetAttribute(cf, landname, label = 'LANDNAME:',default='MAIACRTLS')
    call ESMF_ConfigGetAttribute(cf, landmodel, label = 'LANDMODEL:',default='RTLS')
    if ( (index(lower_to_upper(landmodel),'RTLS') > 0) ) then
      landbandmBRDF =  ESMF_ConfigGetLen(cf, label = 'LANDBAND_C_BRDF:',__RC__)
      allocate (landband_cBRDF(landbandmBRDF))    
      call ESMF_ConfigGetAttribute(cf, landband_cBRDF, label = 'LANDBAND_C_BRDF:', __RC__)
    end if
    if ( (index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0) .or. (lower_to_upper(landmodel) == 'LAMBERTIAN') ) then
      landbandmLER =  ESMF_ConfigGetLen(cf, label = 'LANDBAND_C_LER:',__RC__)
      allocate (landband_cLER(landbandmLER))    
      call ESMF_ConfigGetAttribute(cf, landband_cLER, label = 'LANDBAND_C_LER:', __RC__)
    end if

    ! Water Surface
    call ESMF_ConfigGetAttribute(cf, watername, label = 'WATERNAME:',default='CX')
    call ESMF_ConfigGetAttribute(cf, watermodel, label = 'WATERMODEL:',default='CX')

    ! Input Files
    call ESMF_ConfigGetAttribute(cf, AER_file, label = 'AER_file:',__RC__)
    call ESMF_ConfigGetAttribute(cf, MET_file, label = 'MET_file:',__RC__)
    call ESMF_ConfigGetAttribute(cf, ANG_file, label = 'ANG_file:',__RC__)
    call ESMF_ConfigGetAttribute(cf, INV_file, label = 'INV_file:',__RC__)
    if ( (index(lower_to_upper(landmodel),'RTLS') > 0) ) then
      call ESMF_ConfigGetAttribute(cf, BRDF_file, label = 'BRDF_file:',__RC__)
    end if
    if ( (index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0) .or. (lower_to_upper(landmodel) == 'LAMBERTIAN') ) then
      call ESMF_ConfigGetAttribute(cf, LER_file, label = 'LER_file:',__RC__)
    end if

    call ESMF_ConfigGetAttribute(cf, OUT_file, label = 'OUT_file:',__RC__)
    
    if ( (index(lower_to_upper(landmodel),'BPDF') > 0) ) then
      if ( (index(lower_to_upper(landname),'MAIGNAN') > 0) ) then    
        call ESMF_ConfigGetAttribute(cf, NDVI_file, label = 'NDVI_file:',__RC__)
        call ESMF_ConfigGetAttribute(cf, BPDF_file, label = 'BPDF_file:',__RC__)    
      end if
    end if


    ! Figure out number of channels and read into vector
    !------------------------------------------------------
    nch =  ESMF_ConfigGetLen(cf, label = 'CHANNELS:',__RC__)
    allocate (channels(nch))
    call ESMF_ConfigGetAttribute(cf, channels, label = 'CHANNELS:', default=550.)

    allocate (mr(nch))
    call ESMF_ConfigGetAttribute(cf, mr, label = 'WATERMR:',__RC__) 
    
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

end program pace_vlidort
