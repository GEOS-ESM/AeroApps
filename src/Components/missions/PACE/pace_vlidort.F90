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
  use cloud_MieMod

  implicit none
  include "pace_vlidort_pars.F90"
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
  integer                               :: nstreams               ! number of half space streams, default = 6
  logical                               :: plane_parallel    
  real, allocatable                     :: channels(:)            ! channels to simulate
  real, allocatable                     :: mr(:)                  ! water real refractive index    
  integer                               :: nch                    ! number of channels  
  real                                  :: szamax, vzamax         ! Geomtry filtering
  logical                               :: additional_output      ! does user want additional output
  logical                               :: aerosol_output         ! does user want aerosol output  
  logical                               :: cloud_output           ! does user want cloud output
  integer                               :: nodemax                ! number of nodes requested
  integer                               :: nodenum                ! which node is this?
  character(len=256)                    :: version
  character(len=256)                    :: layout 
  character(len=256)                    :: IcldTable, LcldTable
  integer                               :: idxCld

! Test flag
! -----------
  logical                               :: test_shmem = .False.

! File names
! ----------
  character(len=256)                    :: AER_file, ANG_file, INV_file, BRDF_file, LER_file, OUT_file
  character(len=256)                    :: ADD_file, CLD_file, MET_file, WAT_file, NDVI_file, BPDF_file 
  character(len=256)                    :: AERO_file, CLDO_file

! Global, 3D inputs to be allocated using SHMEM
! ---------------------------------------------
  real, pointer                         :: AIRDENS(:,:,:) => null()
  real, pointer                         :: RH(:,:,:) => null()
  real, pointer                         :: DELP(:,:,:) => null()
  real, pointer                         :: DU001(:,:,:) => null()
  real, pointer                         :: DU002(:,:,:) => null()
  real, pointer                         :: DU003(:,:,:) => null()
  real, pointer                         :: DU004(:,:,:) => null()
  real, pointer                         :: DU005(:,:,:) => null()
  real, pointer                         :: SS001(:,:,:) => null()
  real, pointer                         :: SS002(:,:,:) => null()
  real, pointer                         :: SS003(:,:,:) => null()
  real, pointer                         :: SS004(:,:,:) => null()
  real, pointer                         :: SS005(:,:,:) => null()
  real, pointer                         :: BCPHOBIC(:,:,:) => null()
  real, pointer                         :: BCPHILIC(:,:,:) => null()
  real, pointer                         :: OCPHOBIC(:,:,:) => null()
  real, pointer                         :: OCPHILIC(:,:,:) => null()
  real, pointer                         :: SO4(:,:,:) => null()  
  real, pointer                         :: KISO(:,:,:) => null()
  real, pointer                         :: KVOL(:,:,:) => null()
  real, pointer                         :: KGEO(:,:,:) => null() 
  real, pointer                         :: LER(:,:,:) => null()   
  real, pointer                         :: FRLAND(:,:) => null()
  real, pointer                         :: SZA(:,:) => null()
  real, pointer                         :: VZA(:,:) => null()
  real, pointer                         :: SAA(:,:) => null()
  real, pointer                         :: VAA(:,:) => null()
  real, pointer                         :: RAA(:,:) => null()
  real, pointer                         :: REI(:,:,:) => null()  
  real, pointer                         :: REL(:,:,:) => null()  
  real, pointer                         :: TAUI(:,:,:) => null()  
  real, pointer                         :: TAUL(:,:,:) => null() 
  real, pointer                         :: V10M(:,:) => null()  
  real, pointer                         :: U10M(:,:) => null()  
  real, pointer                         :: SLEAVE(:,:,:,:) => null()
  real, pointer                         :: WATER_CH(:) => null()
  real, pointer                         :: NDVI(:,:) => null()  
  real, pointer                         :: BPDFcoef(:,:) => null()  
  integer, pointer                      :: indices(:) => null()


! VLIDORT input arrays
! ---------------------------
  real, allocatable                     :: Vpe(:,:)                ! edge pressure [Pa]
  real, allocatable                     :: Vze(:,:)                ! edge height above sfc [m]
  real, allocatable                     :: Vte(:,:)                ! edge Temperature [K]
  real, allocatable                     :: Vqm(:,:,:)              ! (mixing ratio) * delp/g
  real, allocatable                     :: Vtau(:,:,:)             ! aerosol optical depth
  real, allocatable                     :: Vssa(:,:,:)             ! single scattering albedo
  real, allocatable                     :: Vg(:,:,:)               ! asymmetry factor
  real*8, allocatable                   :: Valbedo(:,:)            ! surface albedo
  real*8, allocatable                   :: Vsleave(:,:)            ! normalized water leaving radiance [Pa]

  real, allocatable                     :: VtauLcl(:,:,:)          ! cloud aerosol optical depth
  real, allocatable                     :: VssaIcl(:,:,:)          ! cloud aerosol optical depth
  real, allocatable                     :: VssaLcl(:,:,:)          ! cloud single scattering albedo
  real, allocatable                     :: VtauIcl(:,:,:)          ! cloud single scattering albedo
  real, allocatable                     :: VgLcl(:,:,:)            ! cloud asymmetry factor  
  real, allocatable                     :: VgIcl(:,:,:)            ! cloud asymmetry factor    

! VLIDORT output arrays
!-------------------------------
!                                  Intermediate Unshared Arrays
!                                  -----------------------------
  real*8, allocatable                   :: radiance_VL_int(:,:)                   ! TOA normalized radiance from VLIDORT
  real*8, allocatable                   :: reflectance_VL_int(:,:)                ! TOA reflectance from VLIDORT  
  real*8, allocatable                   :: Q_int(:,:)                             ! Q Stokes component
  real*8, allocatable                   :: U_int(:,:)                             ! U Stokes component
  real*8, allocatable                   :: ROT(:,:,:)                             ! rayleigh optical thickness
  real*8, allocatable                   :: depol(:)                               ! rayleigh depolarization ratio
  real*8, allocatable                   :: BR_Q_int(:,:)                          ! surface albedo Q
  real*8, allocatable                   :: BR_U_int(:,:)                          ! surface albedo U


!                                  Final Shared Arrays
!                                  -------------------
  real*8, pointer                       :: radiance_VL(:,:) => null()             ! TOA normalized radiance from VLIDORT
  real*8, pointer                       :: reflectance_VL(:,:) => null()          ! TOA reflectance from VLIDORT
  real*8, pointer                       :: Q(:,:) => null()                      ! Q Stokes component
  real*8, pointer                       :: U(:,:) => null()                      ! U Stokes component
  real*8, pointer                       :: ROD(:,:) => null()                    ! rayleigh optical depth
  real*8, pointer                       :: ALBEDO(:,:) => null()                 ! bi-directional surface reflectance
  real*8, pointer                       :: BR_Q(:,:) => null()                   ! bi-directional surface reflectance Q
  real*8, pointer                       :: BR_U(:,:) => null()                   ! bi-directional surface reflectance U


  real, pointer                         :: TAU(:,:) => null()                  ! aerosol optical depth
  real, pointer                         :: SSA(:,:) => null()                  ! single scattering albedo
  real, pointer                         :: G(:,:) => null()                    ! asymmetry factor
  real, pointer                         :: LTAU(:,:) => null()                 ! liquid cloud optical depth
  real, pointer                         :: LSSA(:,:) => null()                 ! liquid cloud single scattering albedo
  real, pointer                         :: LG(:,:) => null()                   ! liquid cloud asymmetry factor
  real, pointer                         :: ITAU(:,:) => null()                 ! ice cloud optical depth
  real, pointer                         :: ISSA(:,:) => null()                 ! ice cloud single scattering albedo
  real, pointer                         :: IG(:,:) => null()                   ! ice cloud asymmetry factor

  real, pointer                         :: PE(:,:,:) => null()
  real, pointer                         :: ZE(:,:,:) => null()  
  real, pointer                         :: TE(:,:,:) => null()

  real*8,allocatable                    :: AOD(:,:)                                 ! Temporary variable to add up AOD
! VLIDORT working variables
!------------------------------
  integer                               :: ch                                       ! i-channel  
  integer                               :: iband                                    ! i-surfaceband
  real,allocatable                      :: Vpmom(:,:,:,:,:)                         ! elements of scattering phase matrix for vector calculations
  real,allocatable                      :: VpmomIcl(:,:,:,:,:)                      ! elements of scattering phase matrix for vector calculations
  real,allocatable                      :: VpmomLcl(:,:,:,:,:)                      ! elements of scattering phase matrix for vector calculations

! MODIS Kernel variables
!--------------------------
  real*8, allocatable                   :: kernel_wt(:,:,:)                         ! kernel weights (/fiso,fgeo,fvol/)
  real*8, allocatable                   :: param(:,:,:)                             ! Li-Sparse parameters 
                                                                                    ! param1 = crown relative height (h/b)
                                                                                    ! param2 = shape parameter (b/r)
  real                                  :: land_missing                                                                 
  real                                  :: sleave_missing                                                                 
  real*8, allocatable                   :: BPDFparam(:,:,:)                         ! BPDF model parameters (/1.5,NDVI,BPDFcoef/)

! Mie Table Stucture
!---------------------
  type(Chem_Mie)                        :: mieTables

! Satellite domain variables
!------------------------------
  integer                               :: im, jm, km, tm                            ! size of satellite domain
  integer                               :: wnch                                      ! number of wavelength in water leaving file
  integer                               :: i, j, k, n                                ! satellite domain working variable
  integer                               :: starti, counti, endi                      ! array indices and counts for each processor
  integer, allocatable                  :: nclr(:)                                   ! how many clear pixels each processor works on
  integer                               :: clrm                                      ! number of clear pixels for this part of decomposed domain
  integer                               :: clrm_total                                ! number of clear pixels 
  integer                               :: c, cc                                     ! clear pixel working variable
  real, allocatable                     :: SOLAR_ZENITH(:,:)                         ! solar zenith angles used for data filtering
  real,allocatable                      :: SENSOR_ZENITH(:,:)                        ! SENSOR zenith angles used for data filtering  
  logical, allocatable                  :: clmask(:,:)                               ! cloud-land mask
  integer, allocatable                  :: iIndex(:),jIndex(:)                       ! indices of pixels to run

! netcdf variables
!----------------------  
  integer                               :: ncid                                           ! netcdf file id
  integer                               :: varid

! Miscellaneous
! -------------
  integer                               :: ierr, rc, status                            ! MPI error message
  integer                               :: status_mpi(MPI_STATUS_SIZE)                 ! MPI status
  integer                               :: myid, npet, CoresPerNode                    ! MPI dimensions and processor id
  integer                               :: p                                           ! i-processor
  character(len=100)                    :: msg                                         ! message to be printed
  real                                  :: progress                                    ! 
  real                                  :: g5nr_missing  
  logical                               :: do_cxonly, do_cx_sleave                              !

! System tracking variables
! -----------------------------
  integer*8                             :: t1, t2, clock_max
  real*8                                :: clock_rate
  character(len=*), parameter           :: Iam = 'pace_vlidort'

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

! Parse Resource file provided at command line for input info 
! -----------------------------------------------------------
  call getarg(1, arg)
  call getarg(2, nodenumarg)
  call get_config(arg, nodenumarg,ierr)
  if (ierr /= 0) then 
    if (MAPL_am_I_root()) then
        write(*,*) 'Problem reading nodenum'
        if (ierr == 1) then
          write(*,*) 'NODEMAX is > 1.  Must provide nodenum at command line.  Exiting'       
        else if (ierr == 2) then
          write(*,*) 'nodenum > NODEMAX'
        end if                
    end if
    call MAPL_SyncSharedMemory(rc=ierr)
    call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
  end if

! Write out settings to use
! --------------------------
  if (MAPL_am_I_root()) then
    write(*,*) 'Simulating ', lower_to_upper(trim(instname)),' domain on ',date,' ', time, 'Z'
    write(*,*) 'Channels [nm]: ',channels
    write(*,*) 'SZA < ', szamax
    write(*,*) 'VZA < ', vzamax
    if (scalar) write(*,*) 'Scalar calculations'
    if (.not. scalar) write(*,*) 'Vector calculations'
    write(*,*) 'Additional Output: ',additional_output
    write(*,*) 'Aerosol Output: ',aerosol_output
    write(*,*) 'Cloud Output: ',cloud_output
    if (trim(layout) /= '111') write(*,*) 'layout: ',trim(layout)
    write(*,*) ' '
  end if 

! Query for domain dimensions and missing value
!----------------------------------------------
  call mp_readDim("ccd_pixels", CLD_file, im)
  call mp_readDim("number_of_scans", CLD_file, jm)
  call mp_readDim("lev", CLD_file, km)
  call mp_readDim("time", AER_file,tm)

  call mp_readVattr("missing_value", AER_FILE, "DELP", g5nr_missing)
  if ((index(lower_to_upper(landmodel),'RTLS') > 0)) then
    call mp_readVattr("missing_value", BRDF_file, "Riso470", land_missing) 
  else if (lower_to_upper(landmodel) == 'LAMBERTIAN') then
    call mp_readVattr("missing_value", LER_file, "SRFLER354", land_missing) 
  end if

  if ( (index(lower_to_upper(watername),'NOBM') > 0) ) then
    call mp_readDim("wavelength", WAT_file, wnch)
    call mp_readVattr("missing_value", WAT_file, "lwn",sleave_missing)
  end if

! Allocate arrays that will be copied on each processor - unshared
! -----------------------------------------------------------------
  call allocate_unshared()

! Create OUTFILE
! --------------
  if ( MAPL_am_I_root() )  call create_outfile(date, time)

! Create additional data
! ----------------------
  if (additional_output) then
    if ( MAPL_am_I_root() ) call create_addfile(date, time)
  end if

! Createaerosol data
! ----------------------
  if (aerosol_output) then
    if ( MAPL_am_I_root() ) call create_aerfile(date, time)
  end if

! Create cloud data
! ----------------------
  if (cloud_output) then
    if ( MAPL_am_I_root() ) call create_cldfile(date, time)
  end if


  call MAPL_SyncSharedMemory(rc=ierr)
! Read the land and angle data 
! -------------------------------------
  call read_sza()
  call read_vza()

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

  if (.not. ANY(SENSOR_ZENITH < vzamax)) then   
    if (MAPL_am_I_root()) then
      write(*,*) 'View angle too slant, nothing to do'
      write(*,*) 'Exiting.....'
    end if
    GOTO 500
  end if  


! Figure out how many indices to work on
!------------------------------------------
  clrm = 0
  clmask = .False.
  do i=1,im
    do j=1,jm
      if (SOLAR_ZENITH(i,j) < szamax) then
        if (SENSOR_ZENITH(i,j) < vzamax) then
          clrm = clrm + 1
          clmask(i,j) = .True.
        end if
      end if
    end do
  end do

  if (clrm == 0) then
    if (MAPL_am_I_root()) then
      write(*,*) 'No good pixels, nothing to do'
      write(*,*) 'Exiting.....'
    end if
    GOTO 500
  end if

  allocate(iIndex(clrm))
  allocate(jIndex(clrm))
  k = 1
  do i=1,im
    do j=1,jm
      if (clmask(i,j)) then
        iIndex(k) = i
        jIndex(k) = j
        k = k + 1
      end if
    end do
  end do

! Unshared angles no longer needed
!---------------------------------------
  deallocate (SOLAR_ZENITH)
  deallocate (SENSOR_ZENITH)

  if (MAPL_am_I_root()) then
    write(*,*) '<> Created Mask'
    write(*,'(A,I3,A)') '       ',nint(100.*clrm/(im*jm)),'% of the domain will be simulated'
    write(*,'(A,I,A)')  '       Simulating ',clrm,' pixels'
    write(*,*) ' '
  end if   

! Allocate the Global arrays using SHMEM
! It will be available on all processors
! ---------------------------------------------------------
 call allocate_shared()
 call MAPL_SyncSharedMemory(rc=ierr)
  
! ! Read in the global arrays
! ! ------------------------------
  call read_land()
  call read_aer_Nv()
  call read_surf_land()
  call read_water()
  call read_angles()
  call read_cld_Tau()
  call read_wind()

! Wait for everyone to finish reading 
! ------------------------------------------------------------------  
  call MAPL_SyncSharedMemory(rc=ierr)
   

! Split up filtered domain among nodes
! User must verify that there are not more nodes than pixels!
!--------------------------------------
  clrm_total = clrm
  if (nodemax > 1) then
    if (nodenum == nodemax) then
      clrm = clrm/nodemax + mod(clrm,nodemax)
    else
      clrm = clrm/nodemax
    end if
  end if

! Split up filtered domain among processors
!----------------------------------------------
  nclr = 0
  if (npet >= clrm) then
    nclr(1:clrm) = 1
  else if (npet < clrm) then
    nclr(1:npet) = clrm/npet
    nclr(npet)   = nclr(npet) + mod(clrm,npet)
  end if 

! Initialize outputs to be safe
! -------------------------------
  TAU            = dble(MISSING)
  SSA            = dble(MISSING)
  G              = dble(MISSING)
  LTAU           = dble(MISSING)
  LSSA           = dble(MISSING)
  LG             = dble(MISSING)
  ITAU           = dble(MISSING)
  ISSA           = dble(MISSING)
  IG             = dble(MISSING)
  ALBEDO         = dble(MISSING)
  ROT            = dble(MISSING)
  PE             = dble(MISSING)
  ZE             = dble(MISSING)
  TE             = dble(MISSING)    
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
  mieTables = Chem_MieCreate(rcfile,rc)
  if ( rc /= 0 ) then
    print *, 'Cannot create Mie tables from '//trim(rcfile)
    call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
  end if

  if ( nMom > mieTables%nMom ) then ! mieTables%nMom is writen in Aod_EOS.rc file
    print *, 'mieTables do not have enough moments', nMom, mieTables%nMom
    call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
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

! Main do loop over the part of the shuffled domain assinged to each processor
  do cc = starti, endi  
    c = indices(cc)
    c = c + (clrm_total/nodemax)*(nodenum-1)

    i  = iIndex(c)
    j  = jIndex(c)

    call getEdgeVars ( km, nobs, reshape(AIRDENS(i,j,:),(/km,nobs/)), &
                       reshape(DELP(i,j,:),(/km,nobs/)), ptop, &
                       Vpe, Vze, Vte )   
    PE(i,j,:) = Vpe(:,nobs)
    ZE(i,j,:) = Vze(:,nobs)
    TE(i,j,:) = Vte(:,nobs)

    write(msg,'(A,I)') 'getEdgeVars ', myid
    call write_verbose(msg)
    
    call calc_qm()

    write(msg,'(A,I)') 'calc_qm ', myid
    call write_verbose(msg)  

!   Rayleigh Optical Thickness
!   --------------------------
    call VLIDORT_ROT_CALC (km, nch, nobs, dble(channels), dble(Vpe), dble(Vze), dble(Vte), &
                                   dble(MISSING),verbose, &
                                   ROT, depol, ierr )  


!   Aerosol Optical Properties
!   --------------------------
    call VLIDORT_getAOPvector ( mieTables, km, nobs, nch, nq, channels, vnames, verbose, &
                        Vqm, reshape(RH(i,j,:),(/km,nobs/)),&
                        nMom,nPol, Vtau, Vssa, Vg, Vpmom, ierr )
    if ( ierr /= 0 ) then
      print *, 'cannot get aerosol optical properties'
      call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)      
    end if

!   Cloud Optical Properties
!   ------------------------
    call getCOPvector(IcldTable, km, nobs, nch, nMom, nPol, channels, REI(i,j,:), TAUI(i,j,:), VssaIcl, VgIcl, VpmomIcl, VtauIcl)
    call getCOPvector(LcldTable, km, nobs, nch, nMom, nPol, channels, REL(i,j,:), TAUL(i,j,:), VssaLcl, VgLcl, VpmomLcl, VtauLcl)

!   Save some variables on the 2D Grid for Writing Later
!   ------------------------
    ! Rayleigh
    ROD(i,j) = SUM(ROT(:,nobs,nch))

    ! Aerosols
    TAU(i,j) = SUM(Vtau(:,nch,nobs))
    SSA(i,j) = SUM(Vssa(:,nch,nobs)*Vtau(:,nch,nobs))
    G(i,j)   = SUM(Vg(:,nch,nobs)*Vtau(:,nch,nobs))

    ! Clouds
    LTAU(i,j) = SUM(VtauLcl(:,nch,nobs))
    LSSA(i,j) = SUM(VssaLcl(:,nch,nobs)*VtauLcl(:,nch,nobs))
    LG(i,j)   = SUM(VgLcl(:,nch,nobs)*VtauLcl(:,nch,nobs))

    ITAU(i,j) = SUM(VtauIcl(:,nch,nobs))
    ISSA(i,j) = SUM(VssaIcl(:,nch,nobs)*VtauIcl(:,nch,nobs))
    IG(i,j)   = SUM(VgIcl(:,nch,nobs)*VtauIcl(:,nch,nobs))

    if (TAU(i,j) > 0) then
      SSA(i,j) = SSA(i,j)/TAU(i,j)
      G(i,j)   = G(i,j)/TAU(i,j)
    else
      SSA(i,j) = dble(MISSING)
      G(i,j)   = dble(MISSING)
    end if

    if (LTAU(i,j) > 0) then
      LSSA(i,j) = LSSA(i,j)/LTAU(i,j)
      LG(i,j)   = LG(i,j)/LTAU(i,j)
    else
      LSSA(i,j) = dble(MISSING)
      LG(i,j)   = dble(MISSING)
    end if

    if (ITAU(i,j) > 0) then
      ISSA(i,j) = ISSA(i,j)/ITAU(i,j)
      IG(i,j)   = IG(i,j)/ITAU(i,j)
    else
      ISSA(i,j) = dble(MISSING)
      IG(i,j)   = dble(MISSING)
    end if


    write(msg,*) 'getAOP ', myid
    call write_verbose(msg)

!   Call VlIDORT
!   ------------
    if ( FRLAND(i,j) >= 0.99 ) then
      call get_surf_params()
      call DO_LAND()
    else
      call get_ocean_params()
      call DO_OCEAN()
    end if

    
!   Check VLIDORT Status, Store Outputs in Shared Arrays
!   ----------------------------------------------------    
    call mp_check_vlidort(radiance_VL_int,reflectance_VL_int)  
    radiance_VL(i,j)    = radiance_VL_int(nobs,nch)
    reflectance_VL(i,j) = reflectance_VL_int(nobs,nch)
    ALBEDO(i,j) = Valbedo(nobs,nch)
    
    if (.not. scalar) then
      Q(i,j)      = Q_int(nobs,nch)
      U(i,j)      = U_int(nobs,nch)
      BR_Q(i,j)   = BR_Q_int(nobs,nch)
      BR_U(i,j)   = BR_U_int(nobs,nch)
    end if

    write(msg,*) 'VLIDORT Calculations DONE', myid, ierr
    call write_verbose(msg)

!   Keep track of progress of each processor
!   -----------------------------------------        
    if (nint(100.*real(cc-starti)/real(counti)) > progress) then
      progress = progress + 10
      write(*,'(A,I,A,I,A,I2,A,I3,A)') 'Pixel: ',cc,'  End Pixel: ',endi,'  ID:',myid,'  Progress:', nint(progress),'%'           
    end if
                
  end do ! do clear pixels

! Wait for everyone to finish calculations
! ----------------------------------------
  call MAPL_SyncSharedMemory(rc=ierr)

! Write to outfile
! -----------------
  call write_outfile()

! Write additional data
! ----------------------
  if (additional_output) then
    call write_addfile()
  end if

! Write aerosol data
! ----------------------
  if (aerosol_output) then
    call write_aerfile()
  end if

! Write cloud data
! ----------------------
  if (cloud_output) then
    call write_cldfile()
  end if

! Sync processors before shutting down
! ------------------------------------
  call MAPL_SyncSharedMemory(rc=ierr)
  if (MAPL_am_I_root()) then
    write(*,*) '<> Finished VLIDORT Simulation of '//trim(lower_to_upper(instname))//' domain'
    write(*,*) ' '
  end if

! ! All done
! ! --------
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
    integer, dimension(npet)           :: nlayer  ! how many layers each processor writes
    integer                            :: startl,countl,endl  
! Write to OUT_file 
! ------------------     
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

!                             Write to main OUT_File
!                             ----------------------
    call check( nf90_open(OUT_file, IOR(nf90_write, nf90_mpiio), ncid,comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // OUT_file )
    write(msg,'(F10.2)') channels(nch)
    
    call check(nf90_inq_varid(ncid, 'ref_' // trim(adjustl(msg)), varid), "get ref vaird")
    call check(nf90_put_var(ncid, varid, reflectance_VL(:,startl:endl), &
                start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out reflectance")

    call check(nf90_inq_varid(ncid, 'I_' // trim(adjustl(msg)), varid), "get rad vaird")
    call check(nf90_put_var(ncid, varid, radiance_VL(:,startl:endl), &
                  start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out radiance")

    call check(nf90_inq_varid(ncid, 'surf_ref_I_' // trim(adjustl(msg)), varid), "get ref vaird")
    call check(nf90_put_var(ncid, varid, ALBEDO(:,startl:endl), &
                  start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out albedo")

    if (.not. scalar) then
      call check(nf90_inq_varid(ncid, 'Q_' // trim(adjustl(msg)), varid), "get q vaird")
      call check(nf90_put_var(ncid, varid, Q(:,startl:endl), &
                  start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out Q")

      call check(nf90_inq_varid(ncid, 'U_' // trim(adjustl(msg)), varid), "get u vaird")
      call check(nf90_put_var(ncid, varid, U(:,startl:endl), &
                  start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out U")

      call check(nf90_inq_varid(ncid, 'surf_ref_Q_' // trim(adjustl(msg)), varid), "get ref vaird")
      call check(nf90_put_var(ncid, varid, BR_Q(:,startl:endl), &
                    start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out albedo Q")

      call check(nf90_inq_varid(ncid, 'surf_ref_U_' // trim(adjustl(msg)), varid), "get ref vaird")
      call check(nf90_put_var(ncid, varid, BR_U(:,startl:endl), &
                    start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out albedo U")

    endif        

    call check( nf90_close(ncid), "close outfile" )
    if (MAPL_am_I_root()) then
      write(*,*) '<> Wrote Outfile '
    end if


  end subroutine write_outfile

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     write_addfile()
! PURPOSE
!     writes data to addfile
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     Oct 2018 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

  subroutine write_addfile
    integer, dimension(npet)           :: nlayer  ! how many layers each processor writes
    integer                            :: startl,countl,endl    
! Write to ADD_file
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
     
  !                             Write to Additional Outputs File
  !                             --------------------------------

    call check( nf90_open(ADD_file, IOR(nf90_write, nf90_mpiio), ncid,comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // ADD_file )
    write(msg,'(F10.2)') channels(nch)
    call check(nf90_inq_varid(ncid, 'rod_' // trim(adjustl(msg)), varid), "get rod vaird")
    call check(nf90_put_var(ncid, varid, ROD(:,startl:endl), &
              start = (/1,startl,nobs/), count = (/im,countl,nobs/)), "writing out rod")

    call check( nf90_close(ncid), "close addfile" )
    if (MAPL_am_I_root()) then
      write(*,*) '<> Wrote Additional Data file '
    end if
      

  end subroutine write_addfile

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     write_aerfile()
! PURPOSE
!     writes data to aer outfile
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     Oct 2018 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

  subroutine write_aerfile
    integer, dimension(npet)           :: nlayer  ! how many layers each processor writes
    integer                            :: startl,countl,endl     
! Write to AERO_file
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
     
  !                             Write to Aerosol Outputs File
  !                             --------------------------------
    call check( nf90_open(AERO_file, IOR(nf90_write, nf90_mpiio), ncid,comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // AERO_file )
    write(msg,'(F10.2)') channels(nch)

    ! Aerosol Stuff
    call check(nf90_inq_varid(ncid, 'aod_' // trim(adjustl(msg)), varid), "get aot vaird")
    call check(nf90_put_var(ncid, varid, TAU(:,startl:endl), &
              start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out tau")

    call check(nf90_inq_varid(ncid, 'g_' // trim(adjustl(msg)), varid), "get g vaird")
    call check(nf90_put_var(ncid, varid, G(:,startl:endl), &
              start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out g")

    call check(nf90_inq_varid(ncid, 'ssa_' // trim(adjustl(msg)), varid), "get ssa vaird")
    call check(nf90_put_var(ncid, varid, SSA(:,startl:endl), &
              start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out ssa")


    call check( nf90_close(ncid), "close aerofile" )
    if (MAPL_am_I_root()) then
      write(*,*) '<> Wrote Aerosol file '
    end if
      
  end subroutine write_aerfile
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

  subroutine write_cldfile
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
     
  !                             Write to Cloud Outputs File
  !                             --------------------------------
    call check( nf90_open(CLDO_file, IOR(nf90_write, nf90_mpiio), ncid,comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // CLDO_file )
    write(msg,'(F10.2)') channels(nch)

    ! Liquid Cloud Stuff
    call check(nf90_inq_varid(ncid, 'lcod_' // trim(adjustl(msg)), varid), "get lcot vaird")
    call check(nf90_put_var(ncid, varid, LTAU(:,startl:endl), &
              start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out ltau")

    ! call check(nf90_inq_varid(ncid, 'lc_g_' // trim(adjustl(msg)), varid), "get lc_g vaird")
    ! call check(nf90_put_var(ncid, varid, LG(:,startl:endl,k), &
    !           start = (/1,startl,k,nobs/), count = (/im,countl,1,nobs/)), "writing out lg")

    call check(nf90_inq_varid(ncid, 'lc_ssa_' // trim(adjustl(msg)), varid), "get lc_ssa vaird")
    call check(nf90_put_var(ncid, varid, LSSA(:,startl:endl), &
              start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out lssa")

    ! Ice Cloud Stuff
    call check(nf90_inq_varid(ncid, 'icod_' // trim(adjustl(msg)), varid), "get icot vaird")
    call check(nf90_put_var(ncid, varid, ITAU(:,startl:endl), &
              start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out itau")

    ! call check(nf90_inq_varid(ncid, 'ic_g_' // trim(adjustl(msg)), varid), "get ic_g vaird")
    ! call check(nf90_put_var(ncid, varid, IG(:,startl:endl,k), &
    !           start = (/1,startl,k,nobs/), count = (/im,countl,1,nobs/)), "writing out ig")

    call check(nf90_inq_varid(ncid, 'ic_ssa_' // trim(adjustl(msg)), varid), "get ic_ssa vaird")
    call check(nf90_put_var(ncid, varid, ISSA(:,startl:endl), &
              start = (/1,startl,1,nobs/), count = (/im,countl,1,nobs/)), "writing out issa")


    call check( nf90_close(ncid), "close cldofile" )
    
    if (MAPL_am_I_root()) then
      write(*,*) '<> Wrote Cloud file '
    end if

  end subroutine write_cldfile

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
  subroutine get_ocean_params()
    integer              :: below


    if ( (index(lower_to_upper(watername),'NOBM') > 0) ) then
      do ch = 1, nch
        below = minloc(abs(channels(ch) - WATER_CH), dim = 1, mask = (channels(ch) - WATER_CH) .GE. 0)
        if (channels(ch) .eq. maxval(WATER_CH)) then
          Vsleave(ch,nobs) = dble(SLEAVE(i,j,ch,1))
          if (SLEAVE(i,j,ch,1) .eq. sleave_missing) then
            Vsleave(ch,nobs) = dble(MISSING)
          end if
        else
          Vsleave(ch,nobs) = dble(nn_interp(WATER_CH(below:below+1),reshape(SLEAVE(i,j,ch,:),(/2/)),channels(ch),sleave_missing,MISSING))
        end if
      end do

    end if
  end subroutine get_ocean_params



  subroutine DO_OCEAN()
!   WATER PIXELS
!   ------------
    if ( (lower_to_upper(watermodel) .eq. 'CX') ) then

      if  ( (lower_to_upper(watername) .eq. 'CX') ) then
        do_cxonly = .true.
        do_cx_sleave = .false.
      else if ( (index(lower_to_upper(watername),'NOBM') > 0) ) then
        do_cxonly = .false.
        do_cx_sleave = .true.
        if (ANY(Vsleave .eq. dble(MISSING))) then
          do_cx_sleave = .false.
        end if
      end if

      if  ( do_cxonly ) then
      ! GISS Cox Munk Model
      ! -------------------------------
        if (scalar) then
          ! Call to vlidort scalar code       
          call VLIDORT_Scalar_GissCX_Cloud (km, nch, nobs ,dble(channels), nstreams, plane_parallel, nMom,      &
                  nPol, ROT, depol, dble(Vtau), dble(Vssa), dble(Vg), dble(Vpmom), &
                  dble(VtauIcl), dble(VssaIcl), dble(VgIcl), dble(VpmomIcl),&
                  dble(VtauLcl), dble(VssaLcl), dble(VgLcl), dble(VpmomLcl),&                
                  dble(Vpe), dble(Vze), dble(Vte), &
                  (/dble(U10M(i,j))/), &
                  (/dble(V10M(i,j))/), &
                  dble(mr), &
                  (/dble(SZA(i,j))/), &
                  (/dble(abs(RAA(i,j)))/), &
                  (/dble(VZA(i,j))/), &
                  dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, Valbedo, ierr)
        else
          ! Call to vlidort vector code
          call VLIDORT_Vector_GissCX_Cloud (km, nch, nobs ,dble(channels), nstreams, plane_parallel, nMom,   &
                 nPol, ROT, depol, dble(Vtau), dble(Vssa), dble(Vpmom), &
                 dble(VtauIcl), dble(VssaIcl), dble(VpmomIcl),&
                 dble(VtauLcl), dble(VssaLcl), dble(VpmomLcl),&               
                 dble(Vpe), dble(Vze), dble(Vte), &
                 (/dble(U10M(i,j))/), &
                 (/dble(V10M(i,j))/), &
                 dble(mr), &
                 (/dble(SZA(i,j))/), &
                 (/dble(abs(RAA(i,j)))/), &
                 (/dble(VZA(i,j))/), &
                 dble(MISSING),verbose, &
                 radiance_VL_int,reflectance_VL_int, Q_int, U_int, Valbedo, BR_Q_int, BR_U_int, ierr)
        end if
      else if  ( do_cx_sleave ) then
      ! GISS Cox Munk Model with NOBM Water Leaving Reflectance
      ! ----------------------------------------------------------
        if (scalar) then
          ! Call to vlidort scalar code       
          call VLIDORT_Scalar_GissCX_NOBM_Cloud (km, nch, nobs ,dble(channels), nstreams, plane_parallel, nMom,      &
                  nPol, ROT, depol, dble(Vtau), dble(Vssa), dble(Vg), dble(Vpmom), &
                  dble(VtauIcl), dble(VssaIcl), dble(VgIcl), dble(VpmomIcl),&
                  dble(VtauLcl), dble(VssaLcl), dble(VgLcl), dble(VpmomLcl),&                
                  dble(Vpe), dble(Vze), dble(Vte), &
                  (/dble(U10M(i,j))/), &
                  (/dble(V10M(i,j))/), &
                  dble(mr), &
                  Vsleave, &
                  (/dble(SZA(i,j))/), &
                  (/dble(abs(RAA(i,j)))/), &
                  (/dble(VZA(i,j))/), &
                  dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, Valbedo, ierr)
        else
          ! Call to vlidort vector code
          call VLIDORT_Vector_GissCX_NOBM_Cloud (km, nch, nobs ,dble(channels), nstreams, plane_parallel, nMom,   &
                 nPol, ROT, depol, dble(Vtau), dble(Vssa), dble(Vpmom), &
                 dble(VtauIcl), dble(VssaIcl), dble(VpmomIcl),&
                 dble(VtauLcl), dble(VssaLcl), dble(VpmomLcl),&               
                 dble(Vpe), dble(Vze), dble(Vte), &
                 (/dble(U10M(i,j))/), &
                 (/dble(V10M(i,j))/), &
                 dble(mr), &
                 Vsleave, &
                 (/dble(SZA(i,j))/), &
                 (/dble(abs(RAA(i,j)))/), &
                 (/dble(VZA(i,j))/), &
                 dble(MISSING),verbose, &
                 radiance_VL_int,reflectance_VL_int, Q_int, U_int, Valbedo, BR_Q_int, BR_U_int, ierr)
        end if
      else
!       Save code for pixels that were not gap filled
!       ---------------------------------------------    
        radiance_VL_int(nobs,:) = -500
        reflectance_VL_int(nobs,:) = -500
        Valbedo = -500
        if (.not. scalar) then
          Q_int = -500
          U_int = -500
        end if
        ierr = 0

      end if
    end if          

  end subroutine DO_OCEAN


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
            kernel_wt(:,ch,nobs) = (/dble(KISO(i,j,1)),&
                                dble(KGEO(i,j,1)),&
                                dble(KVOL(i,j,1))/)
          else
            ! > 2130 uses highest MODIS wavelength band     
            if (channels(ch) >= maxval(landband_cBRDF)) then  
              iband = minloc(abs(landband_cBRDF - channels(ch)), dim = 1)
              kernel_wt(:,ch,nobs) = (/dble(KISO(i,j,iband)),&
                                dble(KGEO(i,j,iband)),&
                                dble(KVOL(i,j,iband))/)
            end if
            
            if (channels(ch) < maxval(landband_cBRDF)) then
              ! nearest neighbor interpolation of kernel weights to wavelength
              if ( (index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0)  .and. (channels(ch) < minval(landband_cBRDF)) ) then
                !  relax kernel weights below 470 to lamberitan 
                !  ----------------------------
                brdf_i = minloc(landband_cBRDF)
                ler_i  = maxloc(landband_cLER)

                temp_iso(1) = LER(i,j,ler_i(1))
                temp_iso(2) = KISO(i,j,brdf_i(1))

                temp_geo(1) = 0
                temp_geo(2) = KGEO(i,j,brdf_i(1))

                temp_vol(1) = 0
                temp_vol(2) = KVOL(i,j,brdf_i(1))

                kernel_wt(1,ch,nobs) = dble(nn_interp((/landband_cLER(ler_i(1)),landband_cBRDF(brdf_i(1))/),temp_iso,channels(ch),land_missing,MISSING))
                kernel_wt(2,ch,nobs) = dble(nn_interp((/landband_cLER(ler_i(1)),landband_cBRDF(brdf_i(1))/),temp_geo,channels(ch),land_missing,MISSING))
                kernel_wt(3,ch,nobs) = dble(nn_interp((/landband_cLER(ler_i(1)),landband_cBRDF(brdf_i(1))/),temp_vol,channels(ch),land_missing,MISSING))

              else
                kernel_wt(1,ch,nobs) = dble(nn_interp(landband_cBRDF,reshape(KISO(i,j,:),(/landbandmBRDF/)),channels(ch),land_missing,MISSING))
                kernel_wt(2,ch,nobs) = dble(nn_interp(landband_cBRDF,reshape(KGEO(i,j,:),(/landbandmBRDF/)),channels(ch),land_missing,MISSING))
                kernel_wt(3,ch,nobs) = dble(nn_interp(landband_cBRDF,reshape(KVOL(i,j,:),(/landbandmBRDF/)),channels(ch),land_missing,MISSING))    
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
            BPDFparam(2,:,:) = dble(NDVI(i,j))
            BPDFparam(3,:,:) = dble(BPDFcoef(i,j))
          end if
        end if

      else if ( (lower_to_upper(landmodel) == 'LAMBERTIAN') ) then
      !  Get Lambertian Albedo
      !  ------------------------------
        do ch = 1, nch
          if (landbandmLER == 1) then
            Valbedo(nobs,ch) = dble(LER(i,j,1))
          else
            Valbedo(nobs,ch) = dble(nn_interp(landband_cLER,reshape(LER(i,j,:),(/landbandmLER/)),channels(ch),land_missing,MISSING))
          end if
        end do
      end if 

  end subroutine get_surf_params

  subroutine DO_LAND()

  !     LAND PIXELS
  !     -----------

    if ( (lower_to_upper(landmodel) .eq. 'LAMBERTIAN') ) then
    ! Simple lambertian surface model
    ! -------------------------------
      if (scalar) then
          ! Call to vlidort scalar code       
          call VLIDORT_Scalar_Lambert_Cloud (km, nch, nobs ,dble(channels), nstreams, plane_parallel, nMom,      &
                  nPol, ROT, depol, dble(Vtau), dble(Vssa), dble(Vg), dble(Vpmom),&
                  dble(VtauIcl), dble(VssaIcl), dble(VgIcl), dble(VpmomIcl),&
                  dble(VtauLcl), dble(VssaLcl), dble(VgLcl), dble(VpmomLcl),&
                  dble(Vpe), dble(Vze), dble(Vte), Valbedo,&
                  (/dble(SZA(i,j))/), &
                  (/dble(abs(RAA(i,j)))/), &
                  (/dble(VZA(i,j))/), &
                  dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ierr)
      else
        ! Call to vlidort vector code
        call VLIDORT_Vector_Lambert_Cloud (km, nch, nobs ,dble(channels), nstreams, plane_parallel, nMom,   &
               nPol, ROT, depol, dble(Vtau), dble(Vssa), dble(Vpmom), &
               dble(VtauIcl), dble(VssaIcl), dble(VpmomIcl), &
               dble(VtauLcl), dble(VssaLcl), dble(VpmomLcl), &
               dble(Vpe), dble(Vze), dble(Vte), Valbedo,&
               (/dble(SZA(i,j))/), &
               (/dble(abs(RAA(i,j)))/), &
               (/dble(VZA(i,j))/), &
               dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, Q_int, U_int, ierr)
        BR_Q_int = 0
        BR_U_int = 0
      end if

    else if ( (index(lower_to_upper(landmodel),'RTLS') > 0) .and. (index(lower_to_upper(landmodel),'BPDF') .eq. 0) ) then
      if ( ANY(kernel_wt == MISSING) ) then
!       Save code for pixels that were not gap filled
!       ---------------------------------------------    
        radiance_VL_int(nobs,:) = -500
        reflectance_VL_int(nobs,:) = -500
        Valbedo = -500
        if (.not. scalar) then
          Q_int = -500
          U_int = -500
        end if
        ierr = 0

      else   
!       MODIS BRDF Surface Model
!       ------------------------------

        if (scalar) then 
            ! Call to vlidort scalar code            
            call VLIDORT_Scalar_LandMODIS_Cloud (km, nch, nobs, dble(channels), nstreams, plane_parallel, nMom,  &
                    nPol, ROT, depol, dble(Vtau), dble(Vssa), dble(Vg), dble(Vpmom), &
                    dble(VtauIcl), dble(VssaIcl), dble(VgIcl), dble(VpmomIcl), &
                    dble(VtauLcl), dble(VssaLcl), dble(VgLcl), dble(VpmomLcl), &
                    dble(Vpe), dble(Vze), dble(Vte), &
                    kernel_wt, param, &
                    (/dble(SZA(i,j))/), &
                    (/dble(abs(RAA(i,j)))/), &
                    (/dble(VZA(i,j))/), &
                    dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, Valbedo, ierr )  
        else
          ! Call to vlidort vector code
          call VLIDORT_Vector_LandMODIS_cloud (km, nch, nobs, dble(channels), nstreams, plane_parallel, nMom, &
                  nPol, ROT, depol, dble(Vtau), dble(Vssa), dble(Vpmom), &
                  dble(VtauIcl), dble(VssaIcl), dble(VpmomIcl), &
                  dble(VtauLcl), dble(VssaLcl), dble(VpmomLcl), &                
                  dble(Vpe), dble(Vze), dble(Vte), &
                  kernel_wt, param, &
                  (/dble(SZA(i,j))/), &
                  (/dble(abs(RAA(i,j)))/), &
                  (/dble(VZA(i,j))/), &
                  dble(MISSING),verbose, &
                  radiance_VL_int,reflectance_VL_int, Valbedo, Q_int, U_int, BR_Q_int, BR_U_int, ierr )  
        end if    
      end if

    else if ( (index(lower_to_upper(landmodel),'RTLS') > 0) .and. (index(lower_to_upper(landmodel),'BPDF') > 0) ) then  
      if ( ANY(kernel_wt == MISSING) .or. ANY(BPDFparam == land_missing)) then
!       Save code for pixels that were not gap filled
!       ---------------------------------------------    
        radiance_VL_int(nobs,:) = -500
        reflectance_VL_int(nobs,:) = -500
        Valbedo = -500
        if (.not. scalar) then
          Q_int = -500
          U_int = -500
        end if
        ierr = 0

      else   
!       MODIS BPDF Surface Model
!       ------------------------------

        ! Call to vlidort vector code
        call VLIDORT_Vector_LandMODIS_BPDF_cloud (km, nch, nobs, dble(channels), nstreams, plane_parallel, nMom, &
                nPol, ROT, depol, dble(Vtau), dble(Vssa), dble(Vpmom), &
                dble(VtauIcl), dble(VssaIcl), dble(VpmomIcl), &
                dble(VtauLcl), dble(VssaLcl), dble(VpmomLcl), &                
                dble(Vpe), dble(Vze), dble(Vte), &
                kernel_wt, param, BPDFparam, &
                (/dble(SZA(i,j))/), &
                (/dble(abs(RAA(i,j)))/), &
                (/dble(VZA(i,j))/), &
                dble(MISSING),verbose, &
                radiance_VL_int,reflectance_VL_int, Valbedo, Q_int, U_int, BR_Q_int, BR_U_int, ierr )  
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

    call mp_readvar2Dchunk("FRLAND", INV_file, (/im,jm/), 1, npet, myid, FRLAND)

    call MAPL_SyncSharedMemory(rc=ierr)    
    if (MAPL_am_I_root()) then
      write(*,*) '<> Read fraction land data to shared memory'
    end if      


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
    real :: add_offset, scale_factor

    allocate (SOLAR_ZENITH(im,jm))

    call readvar2Dgrp("solar_zenith", "geolocation_data", ANG_file, SOLAR_ZENITH)
    call readVattrgrp("add_offset", "geolocation_data", ANG_file, "solar_zenith", add_offset)
    call readVattrgrp("scale_factor", "geolocation_data", ANG_file, "solar_zenith", scale_factor)

    SOLAR_ZENITH = SOLAR_ZENITH*scale_factor + add_offset

  end subroutine read_sza

  subroutine read_vza()
    real :: add_offset, scale_factor

    allocate (SENSOR_ZENITH(im,jm))

    call readvar2Dgrp("sensor_zenith", "geolocation_data", ANG_file, SENSOR_ZENITH)
    call readVattrgrp("add_offset", "geolocation_data", ANG_file, "sensor_zenith", add_offset)
    call readVattrgrp("scale_factor", "geolocation_data", ANG_file, "sensor_zenith", scale_factor)

    SENSOR_ZENITH = SENSOR_ZENITH*scale_factor + add_offset

  end subroutine read_vza


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
    real*8, dimension(:,:)    :: radiance, reflectance


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

    call mp_readvar3Dchunk("DU001", AER_file, (/im,jm,km/), 1, npet, myid, DU001) 
    call mp_readvar3Dchunk("DU002", AER_file, (/im,jm,km/), 1, npet, myid, DU002) 
    call mp_readvar3Dchunk("DU003", AER_file, (/im,jm,km/), 1, npet, myid, DU003) 
    call mp_readvar3Dchunk("DU004", AER_file, (/im,jm,km/), 1, npet, myid, DU004) 
    call mp_readvar3Dchunk("DU005", AER_file, (/im,jm,km/), 1, npet, myid, DU005) 
    call mp_readvar3Dchunk("SS001", AER_file, (/im,jm,km/), 1, npet, myid, SS001) 
    call mp_readvar3Dchunk("SS002", AER_file, (/im,jm,km/), 1, npet, myid, SS002) 
    call mp_readvar3Dchunk("SS003", AER_file, (/im,jm,km/), 1, npet, myid, SS003) 
    call mp_readvar3Dchunk("SS004", AER_file, (/im,jm,km/), 1, npet, myid, SS004) 
    call mp_readvar3Dchunk("SS005", AER_file, (/im,jm,km/), 1, npet, myid, SS005) 
    call mp_readvar3Dchunk("BCPHOBIC", AER_file, (/im,jm,km/), 1, npet, myid, BCPHOBIC) 
    call mp_readvar3Dchunk("BCPHILIC", AER_file, (/im,jm,km/), 1, npet, myid, BCPHILIC) 
    call mp_readvar3Dchunk("OCPHOBIC", AER_file, (/im,jm,km/), 1, npet, myid, OCPHOBIC) 
    call mp_readvar3Dchunk("OCPHILIC", AER_file, (/im,jm,km/), 1, npet, myid, OCPHILIC) 
    call mp_readvar3Dchunk("SO4", AER_file, (/im,jm,km/), 1, npet, myid, SO4) 
    call mp_readvar3Dchunk("AIRDENS", AER_file, (/im,jm,km/), 1, npet, myid, AIRDENS)  
    call mp_readvar3Dchunk("RH"     , AER_file, (/im,jm,km/), 1, npet, myid, RH) 
    call mp_readvar3Dchunk("DELP"   , AER_file, (/im,jm,km/), 1, npet, myid, DELP) 

    call MAPL_SyncSharedMemory(rc=ierr)    
    if (MAPL_am_I_root()) then
      write(*,*) '<> Read aeorosl data to shared memory'
    end if      

  end subroutine read_aer_Nv

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
!     read_water
! PURPOSE
!     read in water leaving reflectance
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine read_water()

    if ( (index(lower_to_upper(watername),'NOBM') > 0) ) then

      call read_SLEAVE()
      if (MAPL_am_I_root()) then
        write(*,*) '<> Read SLEAVE data to shared memory'
      end if
      call MAPL_SyncSharedMemory(rc=ierr) 

    end if 
  end subroutine read_water

  subroutine read_SLEAVE()
    integer                            :: below
    integer                            :: ncid, varid
    real, dimension(im,jm,2)           :: temp
    real, dimension(im,jm,1)           :: tempmax

    call mp_readvar1Dchunk('wavelength', WAT_file, (/wnch/), 1, npet, myid, WATER_CH)
    if (MAPL_am_I_root()) then
      
      do ch = 1, nch
        ! get channel below
        below = minloc(abs(channels(ch) - WATER_CH), dim = 1, mask = (channels(ch) - WATER_CH) .GE. 0)
        if (channels(ch) .eq. maxval(WATER_CH)) then
          call readvar3Dslice('lwn', WAT_file, (/im,jm,1/), 3, below, tempmax)
          temp = 0
          temp(:,:,1) = tempmax(:,:,1)
        else
          call readvar3Dslice('lwn', WAT_file, (/im,jm,2/), 3, below, temp)          
        end if
        SLEAVE(:,:,ch,:) = temp
        
      end do
    end if
    call MAPL_SyncSharedMemory(rc=ierr) 

  end subroutine read_SLEAVE


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
    ! RTLS Kernel
    if ((index(lower_to_upper(landmodel),'RTLS') > 0)) then
      call read_RTLS()
      if (MAPL_am_I_root()) then
        write(*,*) '<> Read BRDF data to shared memory'
      end if
      call MAPL_SyncSharedMemory(rc=ierr) 

      if ((index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0)) then
        call read_LER()   
        if (MAPL_am_I_root()) then               
          write(*,*) '<> Read LER data to shared memory'
        end if
        call MAPL_SyncSharedMemory(rc=ierr) 
      end if
    end if

    ! BPDF
    if ( (index(lower_to_upper(landmodel),'BPDF') > 0) ) then
      call read_BPDF()
      if (MAPL_am_I_root()) then
        write(*,*) '<> Read BPDF data to shared memory'
      end if
    call MAPL_SyncSharedMemory(rc=ierr) 

    end if

    ! Lambertian
    if (lower_to_upper(landmodel) == 'LAMBERTIAN') then
      call read_LER()     
      if (MAPL_am_I_root()) then           
        write(*,*) '<> Read LER data to shared memory'
      end if
      call MAPL_SyncSharedMemory(rc=ierr) 
    end if 
  end subroutine read_surf_land


  subroutine read_RTLS()
    integer                            :: Nx, Ny, ntile
    integer                            :: xstart, xend, ystart, yend
    integer                            :: x, y
    character(len=100)                 :: sds
    real, pointer                      :: temp(:,:)

    call MAPL_AllocNodeArray(temp,(/im,jm/),rc=ierr)
    call MAPL_SyncSharedMemory(rc=ierr) 
    do ch = 1,landbandmBRDF
      if ( landband_cBRDF(ch) < 1000 ) then
        write(sds,'(A4,I3)') "Riso",int(landband_cBRDF(ch))
      else
        write(sds,'(A4,I4)') "Riso",int(landband_cBRDF(ch))
      end if

      call mp_readvar2Dchunk(trim(sds), BRDF_file, (/im,jm/), 1, npet, myid, temp) 
      call MAPL_SyncSharedMemory(rc=ierr) 
      if (MAPL_am_I_root()) then
        KISO(:,:,ch) = temp
      end if
      call MAPL_SyncSharedMemory(rc=ierr) 

      if ( landband_cBRDF(ch) < 1000 ) then
        write(sds,'(A4,I3)') "Rgeo",int(landband_cBRDF(ch))
      else
        write(sds,'(A4,I4)') "Rgeo",int(landband_cBRDF(ch))
      end if
      
      call mp_readvar2Dchunk(trim(sds), BRDF_file, (/im,jm/), 1, npet, myid, temp) 
      call MAPL_SyncSharedMemory(rc=ierr)
      if (MAPL_am_I_root()) then
        KGEO(:,:,ch) = temp
      end if
      call MAPL_SyncSharedMemory(rc=ierr)         

      if ( landband_cBRDF(ch) < 1000 ) then
        write(sds,'(A4,I3)') "Rvol",int(landband_cBRDF(ch))
      else
        write(sds,'(A4,I4)') "Rvol",int(landband_cBRDF(ch))
      end if

      call mp_readvar2Dchunk(trim(sds), BRDF_file, (/im,jm/), 1, npet, myid, temp) 
      call MAPL_SyncSharedMemory(rc=ierr)
      if (MAPL_am_I_root()) then
        KVOL(:,:,ch) = temp
      end if
      call MAPL_SyncSharedMemory(rc=ierr)         
      
    end do
    call MAPL_DeallocNodeArray(temp,rc=ierr)
  end subroutine read_RTLS

  subroutine read_BPDF()

    call mp_readvar2Dchunk('NDVI', NDVI_file, (/im,jm/), 1, npet, myid, NDVI) 
    call mp_readvar2Dchunk('BPDFcoef', BPDF_file, (/im,jm/), 1, npet, myid, BPDFcoef) 

  end subroutine read_BPDF

  subroutine read_LER()
    !real, dimension(im,jm,1,1)         :: temp
    real, pointer                      :: temp(:,:,:,:)
    character(len=100)                 :: sds

    call MAPL_AllocNodeArray(temp,(/im,jm,1,1/),rc=ierr)
    call MAPL_SyncSharedMemory(rc=ierr)   
    do ch = 1,landbandmLER
      if ( landband_cLER(ch) < 1000 ) then
        write(sds,'(A6,I3)') "SRFLER",int(landband_cLER(ch))
      end if
      !call readvar4d(trim(sds), LER_file, temp)
      call mp_readvar4Dchunk(trim(sds), LER_file, (/im,jm,1,1/), 1, npet, myid, temp) 
      call MAPL_SyncSharedMemory(rc=ierr)  
      if (MAPL_am_I_root()) then 
        LER(:,:,ch) = temp(:,:,1,1)
      end if
      call MAPL_SyncSharedMemory(rc=ierr)     
    end do
    call MAPL_DeallocNodeArray(temp,rc=ierr)
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
    real, dimension(im,jm)     :: temp
    real, allocatable          :: saa_(:,:),angle(:,:)
    integer                    :: i,j
    real                       :: add_offset, scale_factor

    if (MAPL_am_I_root()) then  
      allocate (saa_(im,jm))
      allocate (angle(im,jm))

      call readvar2Dgrp("solar_zenith", "geolocation_data", ANG_file, angle)
      call readVattrgrp("add_offset", "geolocation_data", ANG_file, "solar_zenith", add_offset)
      call readVattrgrp("scale_factor", "geolocation_data", ANG_file, "solar_zenith", scale_factor)

      SZA = angle*dble(scale_factor) + dble(add_offset)


      call readvar2Dgrp("sensor_zenith", "geolocation_data", ANG_file, angle)
      call readVattrgrp("add_offset", "geolocation_data", ANG_file, "sensor_zenith", add_offset)
      call readVattrgrp("scale_factor", "geolocation_data", ANG_file, "sensor_zenith", scale_factor)

      VZA = angle*dble(scale_factor) + dble(add_offset)

      call readvar2Dgrp("solar_azimuth", "geolocation_data", ANG_file, angle)
      call readVattrgrp("add_offset", "geolocation_data", ANG_file, "solar_azimuth", add_offset)
      call readVattrgrp("scale_factor", "geolocation_data", ANG_file, "solar_azimuth", scale_factor)

      SAA = angle*dble(scale_factor) + dble(add_offset)

      call readvar2Dgrp("sensor_azimuth", "geolocation_data", ANG_file, angle)
      call readVattrgrp("add_offset", "geolocation_data", ANG_file, "sensor_azimuth", add_offset)
      call readVattrgrp("scale_factor", "geolocation_data", ANG_file, "sensor_azimuth", scale_factor)

      VAA = angle*dble(scale_factor) + dble(add_offset)

      ! make azimiuths clockwise from north
      do i = 1, im
        do j = 1, jm
          if (SAA(i,j) < 0) then
            SAA(i,j) = 360.0 + SAA(i,j)
          end if
          if (VAA(i,j) < 0) then
            VAA(i,j) = 360.0 + VAA(i,j)
          end if
        end do
      end do
     
      ! define according to photon travel direction
      saa_ = SAA + 180.0
      do i = 1, im
        do j = 1, jm
          if (saa_(i,j) >= 360.0) then
            saa_(i,j) = saa_(i,j) - 360.0
          end if
        end do
      end do

      RAA = VAA - saa_

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


    call mp_readvar2Dchunk("U10M",   MET_file, (/im,jm/), 1, npet, myid, U10M) 
    call mp_readvar2Dchunk("V10M",  MET_file, (/im,jm/), 1, npet, myid, V10M) 
    call MAPL_SyncSharedMemory(rc=ierr)    

    call MAPL_SyncSharedMemory(rc=ierr)    
    if (MAPL_am_I_root()) then
        
      write(*,*) '<> Read wind data to shared memory' 

    end if
    
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
  subroutine allocate_shared()
    if (MAPL_am_I_root()) then
      write(*,*) '<> Allocating shared memory variables'
    end if
    call MAPL_AllocNodeArray(AIRDENS,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(RH,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DELP,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(REI,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(REL,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(TAUI,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(TAUL,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU001,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU002,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU003,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU004,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU005,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS001,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS002,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS003,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS004,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS005,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(BCPHOBIC,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(BCPHILIC,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(OCPHOBIC,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(OCPHILIC,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SO4,(/im,jm,km/),rc=ierr)
    if ((index(lower_to_upper(landmodel),'RTLS') > 0)) then
      call MAPL_AllocNodeArray(KISO,(/im,jm,landbandmBRDF/),rc=ierr)
      call MAPL_AllocNodeArray(KVOL,(/im,jm,landbandmBRDF/),rc=ierr)
      call MAPL_AllocNodeArray(KGEO,(/im,jm,landbandmBRDF/),rc=ierr)
    end if
    if ((lower_to_upper(landmodel) == 'LAMBERTIAN') .or. (index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0)) then
      call MAPL_AllocNodeArray(LER,(/im,jm,landbandmLER/),rc=ierr)
    end if 
    if ( (index(lower_to_upper(landmodel),'BPDF') > 0) ) then
      if ( (index(lower_to_upper(landname),'MAIGNAN') > 0) ) then
        call MAPL_AllocNodeArray(NDVI,(/im,jm/),rc=ierr)
        call MAPL_AllocNodeArray(BPDFcoef,(/im,jm/),rc=ierr)
      end if
    end if

    call MAPL_AllocNodeArray(SZA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(VZA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(RAA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(SAA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(VAA,(/im,jm/),rc=ierr) 

    call MAPL_AllocNodeArray(FRLAND,(/im,jm/),rc=ierr) 

    call MAPL_AllocNodeArray(U10M,(/im,jm/),rc=ierr) 
    call MAPL_AllocNodeArray(V10M,(/im,jm/),rc=ierr) 

    call MAPL_AllocNodeArray(TAU,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(SSA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(G,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(LTAU,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(LSSA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(LG,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(ITAU,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(ISSA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(IG,(/im,jm/),rc=ierr)

    call MAPL_AllocNodeArray(ROD,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(ALBEDO,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(PE,(/im,jm,km+1/),rc=ierr)
    call MAPL_AllocNodeArray(ZE,(/im,jm,km+1/),rc=ierr)
    call MAPL_AllocNodeArray(TE,(/im,jm,km+1/),rc=ierr)
    
    if (.not. scalar) then
      call MAPL_AllocNodeArray(Q,(/im,jm/),rc=ierr)
      call MAPL_AllocNodeArray(U,(/im,jm/),rc=ierr)
      call MAPL_AllocNodeArray(BR_Q,(/im,jm/),rc=ierr)
      call MAPL_AllocNodeArray(BR_U,(/im,jm/),rc=ierr)

    end if

    if ( (index(lower_to_upper(watername),'NOBM') > 0) ) then
      call MAPL_AllocNodeArray(SLEAVE,(/im,jm,nch,2/),rc=ierr)
      call MAPL_AllocNodeArray(WATER_CH,(/wnch/),rc=ierr)
    end if

    call MAPL_AllocNodeArray(radiance_VL,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(reflectance_VL,(/im,jm/),rc=ierr)

  end subroutine allocate_shared

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
    allocate (VtauIcl(km,nch,nobs))
    allocate (VtauLcl(km,nch,nobs))
    allocate (VssaIcl(km,nch,nobs))
    allocate (VssaLcl(km,nch,nobs))    
    allocate (VgIcl(km,nch,nobs))    
    allocate (VgLcl(km,nch,nobs))    
    allocate (Valbedo(nobs,nch))

    allocate (radiance_VL_int(nobs,nch))
    allocate (reflectance_VL_int(nobs, nch))    

    if ((index(lower_to_upper(landmodel),'RTLS') > 0) ) then
      allocate (kernel_wt(nkernel,nch,nobs))
      allocate (param(nparam,nch,nobs))
    end if

    if ( (index(lower_to_upper(landmodel),'BPDF') > 0) ) then  
      if ( (index(lower_to_upper(landname),'MAIGNAN') > 0) ) then    
        allocate (BPDFparam(3,nch,nobs)) 
      end if
    end if
    
    if ( (index(lower_to_upper(watername),'NOBM') > 0) ) then
      allocate (Vsleave(nch,nobs))
    end if

    allocate (ROT(km,nobs,nch))
    allocate (depol(nch))
    allocate (Vpmom(km,nch,nobs,nMom,nPol))
    allocate (VpmomLcl(km,nch,nobs,nMom,nPol))
    allocate (VpmomIcl(km,nch,nobs,nMom,nPol))

    if (.not. scalar) then      
      allocate (Q_int(nobs, nch))
      allocate (U_int(nobs, nch))
      allocate (BR_Q_int(nobs, nch))
      allocate (BR_U_int(nobs, nch))
    end if

  ! Needed for reading
  ! ----------------------
    allocate (nclr(npet)) 
    nclr = 0

    allocate(clmask(im,jm))
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

    integer,dimension(nch)             :: radVarID, refVarID, aotVarID   
    integer,dimension(nch)             :: qVarID, uVarID, albVarID, brUVarID, brQVarID      
    integer,dimension(nch)             :: ssaVarID, tauVarID, gVarID
    integer,dimension(nch)             :: LssaVarID, LtauVarID, LgVarID
    integer,dimension(nch)             :: IssaVarID, ItauVarID, IgVarID
    integer                            :: peVarID, zeVarID, teVarID     
    
    integer                            :: ncid
    integer                            :: timeDimID, ewDimID, nsDimID, levDimID, chaDimID 
    integer                            :: leveDimID      
    integer                            :: scantimeVarID, clonVarID, clatVarID
    integer                            :: timeVarID, levVarID, ewVarID, nsVarID
    integer                            :: leveVarID, e
    integer                            :: ch

    real*8,allocatable,dimension(:,:)  :: clon, clat, sza, vza, raa
    real*8,allocatable,dimension(:)    :: scantime, ew, ns, tyme

    character(len=2000)                :: comment

!                                MAIN LevelC2 OUT_FILE 
!                                ----------------------
    ! Open File
    call check(nf90_create(OUT_file, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // OUT_file)

    ! Create dimensions
    call check(nf90_def_dim(ncid, "time", tm, timeDimID), "creating time dimension")
    call check(nf90_def_dim(ncid, "lev", 1, levDimID), "creating ns dimension") !km
    call check(nf90_def_dim(ncid, "ccd_pixels", im, ewDimID), "creating ew dimension") !im
    call check(nf90_def_dim(ncid, "number_of_scans", jm, nsDimID), "creating ns dimension") !jm

    ! Global Attributes
    write(comment,'(A)') 'VLIDORT Simulation of GEOS-5 '//lower_to_upper(trim(instname))//' Sampler'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'title',trim(comment)),"title attr")

    write(comment,'(A)') 'NASA/Goddard Space Flight Center'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'institution',trim(comment)),"institution attr")

    write(comment,'(A)') 'Global Model and Assimilation Office'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'source',trim(comment)),"source attr")

    write(comment,'(A)') 'VLIDORT simulation run from pace_vlidort.x'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'history',trim(comment)),"history attr")

    call check(nf90_put_att(ncid,NF90_GLOBAL,'grid_inputs',trim(INV_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'angle_inputs',trim(ANG_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'aerosol_inputs',trim(AER_file)),"input files attr")
    if ( (index(lower_to_upper(landmodel),'RTLS') > 0) ) then
      call check(nf90_put_att(ncid,NF90_GLOBAL,'surface_inputs',trim(BRDF_file)),"input files attr")
    else if ( (index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0) ) then
      call check(nf90_put_att(ncid,NF90_GLOBAL,'surface_inputs',trim(BRDF_file)),"input files attr")
      call check(nf90_put_att(ncid,NF90_GLOBAL,'surface_inputs',trim(LER_file)),"input files attr")
    else
      call check(nf90_put_att(ncid,NF90_GLOBAL,'surface_inputs',trim(LER_file)),"input files attr")
    end if
    call check(nf90_put_att(ncid,NF90_GLOBAL,'cloud_inputs',trim(CLD_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'wind_inputs',trim(MET_file)),"input files attr")

    call check(nf90_put_att(ncid,NF90_GLOBAL,'inputs',trim(comment)),"input files attr")
    write(comment,'(A)') 'n/a'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'references',trim(comment)),"references attr") 

    write(comment,'(A)') 'This file contains VLIDORT simulated surface reflectance (albedo) ' // &
                         'and top of the atmosphere ' // &
                         'radiance and reflectance from GEOS-5 parameters sampled ' // &
                         ' on the '//lower_to_upper(trim(instname))//' swath '
    call check(nf90_put_att(ncid,NF90_GLOBAL,'comment',trim(comment)),"comment attr")   

    if ( scalar ) then
      write(comment,'(A)') 'Scalar calculations'
    else
      write(comment,'(A)') 'Vector calculations'
    end if
    call check(nf90_put_att(ncid,NF90_GLOBAL,'vlidort_comment',trim(comment)),"vlidort_comment")

    if (.not. scalar) then
      write(comment,'(I1,A)') nPol, ' components of the scattering matrix '
      call check(nf90_put_att(ncid,NF90_GLOBAL,'scat_comment',trim(comment)),"scat_comment") 
    end if

    write(comment,'(I3,A)') nMom,' phase function moments'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'mie_comment',trim(comment)),"mie_comment")

    write(comment,*) channels
    call check(nf90_put_att(ncid,NF90_GLOBAL,'channels',trim(adjustl(comment))),"channels_comment") 

    call check(nf90_put_att(ncid,NF90_GLOBAL,"contact","Patricia Castellanos <patricia.castellanos@nasa.gov>"),"contact attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","cf"),"conventions attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,"Version",trim(version)),"version attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,"Created",when()),"created attr")    


    write(comment,*) date(1:4),'-',date(5:6),'-',date(7:8),'T',time(1:2),':',time(3:4),':',time(5:6)
    call check(nf90_put_att(ncid,NF90_GLOBAL,'time_coverage_start',trim(comment)),"created time_coverage_start")

    ! Define Variables
!                                     Dimensions
!                                     ----------    
    call check(nf90_def_var(ncid,'lev',nf90_float,(/levDimID/),levVarID),"create lev var")
    call check(nf90_def_var(ncid,'ccd_pixels',nf90_float,(/ewDimID/),ewVarID),"create ew var")
    call check(nf90_def_var(ncid,'number_of_scans',nf90_float,(/nsDimID/),nsVarID),"create ns var")

    call check(nf90_def_var(ncid,'ev_mid_time',nf90_float,(/ewDimID/),scantimeVarID),"create scanTime var")
    call check(nf90_def_var(ncid,'longitude',nf90_float,(/ewDimID,nsDimID/),clonVarID),"create clon var")
    call check(nf90_def_var(ncid,'latitude',nf90_float,(/ewDimID,nsDimID/),clatVarID),"create clat var")

!                                     Data
!                                     ----
    do ch=1,nch
      write(comment,'(F10.2)') channels(ch)
      call check(nf90_def_var(ncid, 'ref_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),refVarID(ch)),"create reflectance var")
      call check(nf90_def_var(ncid, 'surf_ref_I_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),albVarID(ch)),"create albedo var")
      call check(nf90_def_var(ncid, 'I_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),radVarID(ch)),"create radiance var")              

      if (.not. scalar) then
        call check(nf90_def_var(ncid, 'Q_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),qVarID(ch)),"create Q var")
        call check(nf90_def_var(ncid, 'U_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),uVarID(ch)),"create U var")
        call check(nf90_def_var(ncid, 'surf_ref_Q_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),brQVarID(ch)),"create Q var")
        call check(nf90_def_var(ncid, 'surf_ref_U_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),brUVarID(ch)),"create U var")
      end if        

    end do

    ! Variable Attributes
!                                          Reflectance TOA & Surface
!                                          -------------------------  
    do ch=1,size(channels)
      write(comment,'(F10.2,A)') channels(ch), ' nm TOA Reflectance'
      call check(nf90_put_att(ncid,refVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm Top of Atmosphere Reflectance'
      call check(nf90_put_att(ncid,refVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,refVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,refVarID(ch),'units','None'),"units attr")
      call check(nf90_put_att(ncid,refVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

      write(comment,'(F10.2,A)') channels(ch), ' nm Surface Reflectance I'
      call check(nf90_put_att(ncid,albVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm Bi-Directional Surface Reflectance'
      call check(nf90_put_att(ncid,albVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,albVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,albVarID(ch),'units','None'),"units attr")
      call check(nf90_put_att(ncid,albVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

      write(comment,'(F10.2,A)') channels(ch), ' nm TOA I'
      call check(nf90_put_att(ncid,radVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm TAO I Component of the Stokes Vector'
      call check(nf90_put_att(ncid,radVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,radVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,radVarID(ch),'units','W m-2 sr-1 nm-1'),"units attr")
      call check(nf90_put_att(ncid,radVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

      if (.not. scalar) then
        write(comment,'(F10.2,A)') channels(ch), ' nm TOA Q'
        call check(nf90_put_att(ncid,qVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm TOA Q Component of the Stokes Vector'
        call check(nf90_put_att(ncid,qVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,qVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,qVarID(ch),'units','W m-2 sr-1 nm-1'),"units attr")
        call check(nf90_put_att(ncid,qVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

        write(comment,'(F10.2,A)') channels(ch), ' nm TOA U'
        call check(nf90_put_att(ncid,uVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm TOA U Component of the Stokes Vector'
        call check(nf90_put_att(ncid,uVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,uVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,uVarID(ch),'units','W m-2 sr-1 nm-1'),"units attr")
        call check(nf90_put_att(ncid,uVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

        write(comment,'(F10.2,A)') channels(ch), ' nm Surface Reflectance Q'
        call check(nf90_put_att(ncid,brQVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm Bi-Directional Surface Reflectance'
        call check(nf90_put_att(ncid,brQVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,brQVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,brQVarID(ch),'units','None'),"units attr")
        call check(nf90_put_att(ncid,brQVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

        write(comment,'(F10.2,A)') channels(ch), ' nm Surface Reflectance U'
        call check(nf90_put_att(ncid,brUVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm Bi-Directional Surface Reflectance'
        call check(nf90_put_att(ncid,brUVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,brUVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,brUVarID(ch),'units','None'),"units attr")
        call check(nf90_put_att(ncid,brUVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

      end if

    end do

!                                          scanTime
!                                          -------  
    call check(nf90_put_att(ncid,scantimeVarID,'long_name','Earth view midtime (seconds of day)'),"long_name attr")
    call check(nf90_put_att(ncid,scantimeVarID,'units','seconds'),"units attr")

!                                          EW, NS, LEV, TIME
!                                          -----------------------  
    call check(nf90_put_att(ncid,ewVarID,'long_name','pseudo longitude'),"long_name attr")
    call check(nf90_put_att(ncid,ewVarID,'units','degrees_east'),"units attr")
    call check(nf90_put_att(ncid,nsVarID,'long_name','pseudo latitude'),"long_name attr")
    call check(nf90_put_att(ncid,nsVarID,'units','degrees_north'),"units attr")   
    call check(nf90_put_att(ncid,levVarID,'long_name','Vertical Level'),"long_name attr")
    call check(nf90_put_att(ncid,levVarID,'units','layer'),"units attr")
    call check(nf90_put_att(ncid,levVarID,'positive','down'),"positive attr")
    call check(nf90_put_att(ncid,levVarID,'axis','z'),"axis attr")

!                                          clon & clat
!                                          -------  
    call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
    call check(nf90_put_att(ncid,clonVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
    call check(nf90_put_att(ncid,clatVarID,'missing_value',real(MISSING)),"missing_value attr")

    call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
    call check(nf90_put_att(ncid,clonVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
    call check(nf90_put_att(ncid,clatVarID,'missing_value',real(MISSING)),"missing_value attr")  
      

    !Leave define mode
    call check(nf90_enddef(ncid),"leaving define mode")

    ! write out ew, ns, lev, time, clon, clat, & scantime
    allocate (scantime(im))
    allocate (clon(im, jm))
    allocate (clat(im, jm))
    allocate (ew(im))
    allocate (ns(jm))    

    call readvar1D("ev_mid_time", INV_file, scantime)
    call check(nf90_put_var(ncid,scantimeVarID,scantime), "writing out scantime")

    call readvar2D("longitude", INV_file, clon)
    call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

    call readvar2D("latitude", INV_file, clat)
    call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

    call check(nf90_put_var(ncid,levVarID,(/1/)), "writing out lev")

    call readvar1D("ccd_pixels", INV_file, ew)
    call check(nf90_put_var(ncid,ewVarID,ew), "writing out ew")

    call readvar1D("number_of_scans", INV_file, ns)
    call check(nf90_put_var(ncid,nsVarID,ns), "writing out ns")   

    deallocate (clon)
    deallocate (clat)  
    deallocate (scantime)
    deallocate (ns)
    deallocate (ew)

    call check( nf90_close(ncid), "close outfile" )

  end subroutine create_outfile

  subroutine create_aerfile(date, time)
    character(len=*) ,intent(in)       :: date, time

    integer,dimension(nch)             :: ssaVarID, tauVarID, gVarID
    
    integer                            :: ncid
    integer                            :: timeDimID, ewDimID, nsDimID, levDimID, chaDimID 
    integer                            :: scantimeVarID, clonVarID, clatVarID
    integer                            :: timeVarID, levVarID, ewVarID, nsVarID
    integer                            :: leveVarID, e
    integer                            :: ch

    real*8,allocatable,dimension(:,:)  :: clon, clat, sza, vza, raa
    real*8,allocatable,dimension(:)    :: scantime, ew, ns, tyme, lev

    character(len=2000)                :: comment

!                           AEROSOL OUTPUTS
!                           ------------------
      ! Open File
      call check(nf90_create(AERO_file, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // AERO_file)

      ! Create dimensions
      call check(nf90_def_dim(ncid, "time", tm, timeDimID), "creating time dimension")
      call check(nf90_def_dim(ncid, "lev", 1, levDimID), "creating ns dimension") !km
      call check(nf90_def_dim(ncid, "ccd_pixels", im, ewDimID), "creating ew dimension") !im
      call check(nf90_def_dim(ncid, "number_of_scans", jm, nsDimID), "creating ns dimension") !jm

      ! Global Attributes
      write(comment,'(A)') 'Atmospheric inputs for VLIDORT Simulation of GEOS-5 '//lower_to_upper(trim(instname))//' Sampler'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'title',trim(comment)),"title attr")

      write(comment,'(A)') 'NASA/Goddard Space Flight Center'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'institution',trim(comment)),"institution attr")

      write(comment,'(A)') 'Global Model and Assimilation Office'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'source',trim(comment)),"source attr")

      write(comment,'(A)') 'VLIDORT simulation run from pace_vlidort.x'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'history',trim(comment)),"history attr")

      write(comment,'(A)') 'n/a'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'references',trim(comment)),"references attr") 

      write(comment,'(A)') 'This file contains intermediate input data for the VLIDORT simulation'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'comment',trim(comment)),"comment attr")   

      if ( scalar ) then
        write(comment,'(A)') 'Scalar calculations'
      else
        write(comment,'(A)') 'Vector calculations'
      end if
      call check(nf90_put_att(ncid,NF90_GLOBAL,'vlidort_comment',trim(comment)),"vlidort_comment")

      if (.not. scalar) then
        write(comment,'(I1,A)') nPol, ' components of the scattering matrix '
        call check(nf90_put_att(ncid,NF90_GLOBAL,'scat_comment',trim(comment)),"scat_comment") 
      end if

      write(comment,'(I3,A)') nMom,' phase function moments'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'mie_comment',trim(comment)),"mie_comment")

      write(comment,*) channels
      call check(nf90_put_att(ncid,NF90_GLOBAL,'channels',trim(adjustl(comment))),"channels_comment") 

      call check(nf90_put_att(ncid,NF90_GLOBAL,"contact","Patricia Castellanos <patricia.castellanos@nasa.gov>"),"contact attr")
      call check(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","cf"),"conventions attr")

      write(comment,*) date(1:4),'-',date(5:6),'-',date(7:8),'T',time(1:2),':',time(3:4),':',time(5:6)
      call check(nf90_put_att(ncid,NF90_GLOBAL,'time_coverage_start',trim(comment)),"created time_coverage_start")

      ! Define Variables
  !                                     Dimensions
  !                                     ----------    
      call check(nf90_def_var(ncid,'lev',nf90_float,(/levDimID/),levVarID),"create lev var")
      call check(nf90_def_var(ncid,'ccd_pixels',nf90_float,(/ewDimID/),ewVarID),"create ew var")
      call check(nf90_def_var(ncid,'number_of_scans',nf90_float,(/nsDimID/),nsVarID),"create ns var")

      call check(nf90_def_var(ncid,'ev_mid_time',nf90_float,(/ewDimID/),scantimeVarID),"create scanTime var")
      call check(nf90_def_var(ncid,'longitude',nf90_float,(/ewDimID,nsDimID/),clonVarID),"create clon var")
      call check(nf90_def_var(ncid,'latitude',nf90_float,(/ewDimID,nsDimID/),clatVarID),"create clat var")

  !                                     Data
  !                                     ----
      do ch=1,nch
        write(comment,'(F10.2)') channels(ch)

        call check(nf90_def_var(ncid, 'aod_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),tauVarID(ch)),"create aot var")      
        call check(nf90_def_var(ncid, 'ssa_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),ssaVarID(ch)),"create ssa var")
        call check(nf90_def_var(ncid, 'g_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),gVarID(ch)),"create g var")

      end do
      ! Variable Attributes
  !                                          Aerosol Data
  !                                          -----------------  
      do ch=1,size(channels)

        ! Aerosol Stuff
        write(comment,'(F10.2,A)') channels(ch), ' nm AOD'
        call check(nf90_put_att(ncid,tauVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm Aerosol Optical Depth'
        call check(nf90_put_att(ncid,tauVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,tauVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,tauVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,tauVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

        write(comment,'(F10.2,A)') channels(ch), ' nm SSA'
        call check(nf90_put_att(ncid,ssaVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm Aerosol Single Scattering Albedo'
        call check(nf90_put_att(ncid,ssaVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,ssaVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,ssaVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,ssaVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")   

        write(comment,'(F10.2,A)') channels(ch), ' nm g'
        call check(nf90_put_att(ncid,gVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm Aerosol Asymmetry Parameter'
        call check(nf90_put_att(ncid,gVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,gVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,gVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,gVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")  

      end do


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
      call check(nf90_put_att(ncid,levVarID,'long_name','Vertical Level'),"long_name attr")
      call check(nf90_put_att(ncid,levVarID,'units','layer'),"units attr")
      call check(nf90_put_att(ncid,levVarID,'positive','down'),"positive attr")
      call check(nf90_put_att(ncid,levVarID,'axis','z'),"axis attr")

  !                                          clon & clat
  !                                          -------  
      call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
      call check(nf90_put_att(ncid,clonVarID,'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
      call check(nf90_put_att(ncid,clatVarID,'missing_value',real(MISSING)),"missing_value attr")

      call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
      call check(nf90_put_att(ncid,clonVarID,'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
      call check(nf90_put_att(ncid,clatVarID,'missing_value',real(MISSING)),"missing_value attr")  
        

      !Leave define mode
      call check(nf90_enddef(ncid),"leaving define mode")

      ! write out ew, ns, lev, time, clon, clat, & scantime
      allocate (scantime(im))
      allocate (clon(im, jm))
      allocate (clat(im, jm))
      allocate (ew(im))
      allocate (ns(jm))    
      allocate (lev(km))

      call readvar1D("ev_mid_time", INV_file, scantime)
      call check(nf90_put_var(ncid,scantimeVarID,scantime), "writing out scantime")

      call readvar2D("longitude", INV_file, clon)
      call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

      call readvar2D("latitude", INV_file, clat)
      call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

      call readvar1D("lev", AER_file, lev)
      call check(nf90_put_var(ncid,levVarID,(/1/)), "writing out lev")

      call readvar1D("ccd_pixels", INV_file, ew)
      call check(nf90_put_var(ncid,ewVarID,ew), "writing out ew")

      call readvar1D("number_of_scans", INV_file, ns)
      call check(nf90_put_var(ncid,nsVarID,ns), "writing out ns")   

      deallocate (clon)
      deallocate (clat)  
      deallocate (scantime)
      deallocate (ns)
      deallocate (ew)
      deallocate (lev)

      call check( nf90_close(ncid), "close aerfile" )
  end subroutine create_aerfile  

  subroutine create_cldfile(date, time)
    character(len=*) ,intent(in)       :: date, time

    integer,dimension(nch)             :: LssaVarID, LtauVarID, LgVarID
    integer,dimension(nch)             :: IssaVarID, ItauVarID, IgVarID
    
    integer                            :: ncid
    integer                            :: timeDimID, ewDimID, nsDimID, levDimID, chaDimID 
    integer                            :: scantimeVarID, clonVarID, clatVarID
    integer                            :: timeVarID, levVarID, ewVarID, nsVarID
    integer                            :: leveVarID, e
    integer                            :: ch

    real*8,allocatable,dimension(:,:)  :: clon, clat, sza, vza, raa
    real*8,allocatable,dimension(:)    :: scantime, ew, ns, tyme, lev

    character(len=2000)                :: comment

!                           CLOUD OUTPUTS
!                           ------------------
      ! Open File
      call check(nf90_create(CLDO_file, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // CLDO_file)

      ! Create dimensions
      call check(nf90_def_dim(ncid, "time", tm, timeDimID), "creating time dimension")
      call check(nf90_def_dim(ncid, "lev", 1, levDimID), "creating ns dimension") !km
      call check(nf90_def_dim(ncid, "ccd_pixels", im, ewDimID), "creating ew dimension") !im
      call check(nf90_def_dim(ncid, "number_of_scans", jm, nsDimID), "creating ns dimension") !jm

      ! Global Attributes
      write(comment,'(A)') 'Atmospheric inputs for VLIDORT Simulation of GEOS-5 '//lower_to_upper(trim(instname))//' Sampler'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'title',trim(comment)),"title attr")

      write(comment,'(A)') 'NASA/Goddard Space Flight Center'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'institution',trim(comment)),"institution attr")

      write(comment,'(A)') 'Global Model and Assimilation Office'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'source',trim(comment)),"source attr")

      write(comment,'(A)') 'VLIDORT simulation run from pace_vlidort.x'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'history',trim(comment)),"history attr")

      write(comment,'(A)') 'n/a'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'references',trim(comment)),"references attr") 

      write(comment,'(A)') 'This file contains intermediate input data for the VLIDORT simulation'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'comment',trim(comment)),"comment attr")   

      if ( scalar ) then
        write(comment,'(A)') 'Scalar calculations'
      else
        write(comment,'(A)') 'Vector calculations'
      end if
      call check(nf90_put_att(ncid,NF90_GLOBAL,'vlidort_comment',trim(comment)),"vlidort_comment")

      if (.not. scalar) then
        write(comment,'(I1,A)') nPol, ' components of the scattering matrix '
        call check(nf90_put_att(ncid,NF90_GLOBAL,'scat_comment',trim(comment)),"scat_comment") 
      end if

      write(comment,'(I3,A)') nMom,' phase function moments'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'mie_comment',trim(comment)),"mie_comment")

      write(comment,*) channels
      call check(nf90_put_att(ncid,NF90_GLOBAL,'channels',trim(adjustl(comment))),"channels_comment") 

      call check(nf90_put_att(ncid,NF90_GLOBAL,"contact","Patricia Castellanos <patricia.castellanos@nasa.gov>"),"contact attr")
      call check(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","cf"),"conventions attr")

      write(comment,*) date(1:4),'-',date(5:6),'-',date(7:8),'T',time(1:2),':',time(3:4),':',time(5:6)
      call check(nf90_put_att(ncid,NF90_GLOBAL,'time_coverage_start',trim(comment)),"created time_coverage_start")

      ! Define Variables
  !                                     Dimensions
  !                                     ----------    
      call check(nf90_def_var(ncid,'lev',nf90_float,(/levDimID/),levVarID),"create lev var")
      call check(nf90_def_var(ncid,'ccd_pixels',nf90_float,(/ewDimID/),ewVarID),"create ew var")
      call check(nf90_def_var(ncid,'number_of_scans',nf90_float,(/nsDimID/),nsVarID),"create ns var")

      call check(nf90_def_var(ncid,'ev_mid_time',nf90_float,(/ewDimID/),scantimeVarID),"create scanTime var")
      call check(nf90_def_var(ncid,'longitude',nf90_float,(/ewDimID,nsDimID/),clonVarID),"create clon var")
      call check(nf90_def_var(ncid,'latitude',nf90_float,(/ewDimID,nsDimID/),clatVarID),"create clat var")

  !                                     Data
  !                                     ----
      do ch=1,nch
        write(comment,'(F10.2)') channels(ch)

        call check(nf90_def_var(ncid, 'lcod_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),LtauVarID(ch)),"create lcot var")      
        call check(nf90_def_var(ncid, 'lc_ssa_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),LssaVarID(ch)),"create lc_ssa var")
        call check(nf90_def_var(ncid, 'lc_g_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),LgVarID(ch)),"create lc_g var")

        call check(nf90_def_var(ncid, 'icod_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),ItauVarID(ch)),"create icot var")      
        call check(nf90_def_var(ncid, 'ic_ssa_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),IssaVarID(ch)),"create ic_ssa var")
        call check(nf90_def_var(ncid, 'ic_g_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),IgVarID(ch)),"create ic_g var")

      end do
      ! Variable Attributes
  !                                          Cloud Data
  !                                          -----------------  
      do ch=1,size(channels)
        ! Liquid Cloud Stuff
        write(comment,'(F10.2,A)') channels(ch), ' nm Liquid COD'
        call check(nf90_put_att(ncid,LtauVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm Liquid Cloud Optical Depth'
        call check(nf90_put_att(ncid,LtauVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,LtauVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,LtauVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,LtauVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

        write(comment,'(F10.2,A)') channels(ch), ' nm Liquid Cloud SSA'
        call check(nf90_put_att(ncid,LssaVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm Liquid Cloud Single Scattering Albedo'
        call check(nf90_put_att(ncid,LssaVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,LssaVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,LssaVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,LssaVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")   

        ! write(comment,'(F10.2,A)') channels(ch), ' nm Liquid Cloud g'
        ! call check(nf90_put_att(ncid,LgVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        ! write(comment,'(F10.2,A)') channels(ch), ' nm layer Liquid Cloud Asymmetry Parameter'
        ! call check(nf90_put_att(ncid,LgVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        ! call check(nf90_put_att(ncid,LgVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        ! call check(nf90_put_att(ncid,LgVarID(ch),'units','none'),"units attr")
        ! call check(nf90_put_att(ncid,LgVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")  

        ! Ice Cloud Stuff
        write(comment,'(F10.2,A)') channels(ch), ' nm Ice COD'
        call check(nf90_put_att(ncid,ItauVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm Ice Cloud Optical Depth'
        call check(nf90_put_att(ncid,ItauVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,ItauVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,ItauVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,ItauVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

        write(comment,'(F10.2,A)') channels(ch), ' nm Ice Cloud SSA'
        call check(nf90_put_att(ncid,IssaVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm Ice Cloud Single Scattering Albedo'
        call check(nf90_put_att(ncid,IssaVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,IssaVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,IssaVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,IssaVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")   

        ! write(comment,'(F10.2,A)') channels(ch), ' nm Ice Cloud g'
        ! call check(nf90_put_att(ncid,IgVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        ! write(comment,'(F10.2,A)') channels(ch), ' nm layer Ice Cloud Asymmetry Parameter'
        ! call check(nf90_put_att(ncid,IgVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        ! call check(nf90_put_att(ncid,IgVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        ! call check(nf90_put_att(ncid,IgVarID(ch),'units','none'),"units attr")
        ! call check(nf90_put_att(ncid,IgVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")  

      end do


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
      call check(nf90_put_att(ncid,levVarID,'long_name','Vertical Level'),"long_name attr")
      call check(nf90_put_att(ncid,levVarID,'units','layer'),"units attr")
      call check(nf90_put_att(ncid,levVarID,'positive','down'),"positive attr")
      call check(nf90_put_att(ncid,levVarID,'axis','z'),"axis attr")

  !                                          clon & clat
  !                                          -------  
      call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
      call check(nf90_put_att(ncid,clonVarID,'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
      call check(nf90_put_att(ncid,clatVarID,'missing_value',real(MISSING)),"missing_value attr")

      call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
      call check(nf90_put_att(ncid,clonVarID,'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
      call check(nf90_put_att(ncid,clatVarID,'missing_value',real(MISSING)),"missing_value attr")  
        

      !Leave define mode
      call check(nf90_enddef(ncid),"leaving define mode")

      ! write out ew, ns, lev, time, clon, clat, & scantime
      allocate (scantime(im))
      allocate (clon(im, jm))
      allocate (clat(im, jm))
      allocate (ew(im))
      allocate (ns(jm))    
      allocate (lev(km))

      call readvar1D("ev_mid_time", INV_file, scantime)
      call check(nf90_put_var(ncid,scantimeVarID,scantime), "writing out scantime")

      call readvar2D("longitude", INV_file, clon)
      call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

      call readvar2D("latitude", INV_file, clat)
      call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

      call readvar1D("lev", AER_file, lev)
      call check(nf90_put_var(ncid,levVarID,(/1/)), "writing out lev")

      call readvar1D("ccd_pixels", INV_file, ew)
      call check(nf90_put_var(ncid,ewVarID,ew), "writing out ew")

      call readvar1D("number_of_scans", INV_file, ns)
      call check(nf90_put_var(ncid,nsVarID,ns), "writing out ns")   

      deallocate (clon)
      deallocate (clat)  
      deallocate (scantime)
      deallocate (ns)
      deallocate (ew)
      deallocate (lev)

      call check( nf90_close(ncid), "close cldofile" )
  end subroutine create_cldfile  


  subroutine create_addfile(date, time)
    character(len=*) ,intent(in)       :: date, time

    integer,dimension(nch)             :: rodVarID
    integer                            :: peVarID, zeVarID, teVarID     
    
    integer                            :: ncid
    integer                            :: timeDimID, ewDimID, nsDimID, levDimID, chaDimID 
    integer                            :: scantimeVarID, clonVarID, clatVarID
    integer                            :: timeVarID, levVarID, ewVarID, nsVarID
    integer                            :: leveVarID, e
    integer                            :: ch

    real*8,allocatable,dimension(:,:)  :: clon, clat, sza, vza, raa
    real*8,allocatable,dimension(:)    :: scantime, ew, ns, tyme, lev

    character(len=2000)                :: comment

!                           ADDITIONAL OUTPUTS
!                           ------------------
      ! Open File
      call check(nf90_create(ADD_file, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // ADD_file)

      ! Create dimensions
      call check(nf90_def_dim(ncid, "time", tm, timeDimID), "creating time dimension")
      call check(nf90_def_dim(ncid, "lev", km, levDimID), "creating ns dimension") !km
      call check(nf90_def_dim(ncid, "ccd_pixels", im, ewDimID), "creating ew dimension") !im
      call check(nf90_def_dim(ncid, "number_of_scans", jm, nsDimID), "creating ns dimension") !jm

      ! Global Attributes
      write(comment,'(A)') 'Atmospheric inputs for VLIDORT Simulation of GEOS-5 '//lower_to_upper(trim(instname))//' Sampler'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'title',trim(comment)),"title attr")

      write(comment,'(A)') 'NASA/Goddard Space Flight Center'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'institution',trim(comment)),"institution attr")

      write(comment,'(A)') 'Global Model and Assimilation Office'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'source',trim(comment)),"source attr")

      write(comment,'(A)') 'VLIDORT simulation run from pace_vlidort.x'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'history',trim(comment)),"history attr")

      write(comment,'(A)') 'n/a'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'references',trim(comment)),"references attr") 

      write(comment,'(A)') 'This file contains intermediate input data for the VLIDORT simulation'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'comment',trim(comment)),"comment attr")   

      if ( scalar ) then
        write(comment,'(A)') 'Scalar calculations'
      else
        write(comment,'(A)') 'Vector calculations'
      end if
      call check(nf90_put_att(ncid,NF90_GLOBAL,'vlidort_comment',trim(comment)),"vlidort_comment")

      if (.not. scalar) then
        write(comment,'(I1,A)') nPol, ' components of the scattering matrix '
        call check(nf90_put_att(ncid,NF90_GLOBAL,'scat_comment',trim(comment)),"scat_comment") 
      end if

      write(comment,'(I3,A)') nMom,' phase function moments'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'mie_comment',trim(comment)),"mie_comment")

      write(comment,*) channels
      call check(nf90_put_att(ncid,NF90_GLOBAL,'channels',trim(adjustl(comment))),"channels_comment") 

      call check(nf90_put_att(ncid,NF90_GLOBAL,"contact","Patricia Castellanos <patricia.castellanos@nasa.gov>"),"contact attr")
      call check(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","cf"),"conventions attr")

      write(comment,*) date(1:4),'-',date(5:6),'-',date(7:8),'T',time(1:2),':',time(3:4),':',time(5:6)
      call check(nf90_put_att(ncid,NF90_GLOBAL,'time_coverage_start',trim(comment)),"created time_coverage_start")

      ! Define Variables
  !                                     Dimensions
  !                                     ----------    
      call check(nf90_def_var(ncid,'lev',nf90_float,(/levDimID/),levVarID),"create lev var")
      call check(nf90_def_var(ncid,'ccd_pixels',nf90_float,(/ewDimID/),ewVarID),"create ew var")
      call check(nf90_def_var(ncid,'number_of_scans',nf90_float,(/nsDimID/),nsVarID),"create ns var")

      call check(nf90_def_var(ncid,'ev_mid_time',nf90_float,(/ewDimID/),scantimeVarID),"create scanTime var")
      call check(nf90_def_var(ncid,'longitude',nf90_float,(/ewDimID,nsDimID/),clonVarID),"create clon var")
      call check(nf90_def_var(ncid,'latitude',nf90_float,(/ewDimID,nsDimID/),clatVarID),"create clat var")

  !                                     Data
  !                                     ----
      do ch=1,nch
        write(comment,'(F10.2)') channels(ch)
        call check(nf90_def_var(ncid, 'rod_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,timeDimID/),rodVarID(ch)),"create rod var")
      end do
      ! Variable Attributes
  !                                          Additional Data
  !                                          -----------------  
      do ch=1,size(channels)
        write(comment,'(F10.2,A)') channels(ch), ' nm ROD'
        call check(nf90_put_att(ncid,rodVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm column Rayliegh Optical Depth'
        call check(nf90_put_att(ncid,rodVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,rodVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,rodVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,rodVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

      end do


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
      call check(nf90_put_att(ncid,levVarID,'long_name','Vertical Level'),"long_name attr")
      call check(nf90_put_att(ncid,levVarID,'units','layer'),"units attr")
      call check(nf90_put_att(ncid,levVarID,'positive','down'),"positive attr")
      call check(nf90_put_att(ncid,levVarID,'axis','z'),"axis attr")

  !                                          clon & clat
  !                                          -------  
      call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
      call check(nf90_put_att(ncid,clonVarID,'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
      call check(nf90_put_att(ncid,clatVarID,'missing_value',real(MISSING)),"missing_value attr")

      call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
      call check(nf90_put_att(ncid,clonVarID,'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
      call check(nf90_put_att(ncid,clatVarID,'missing_value',real(MISSING)),"missing_value attr")  
        

      !Leave define mode
      call check(nf90_enddef(ncid),"leaving define mode")

      ! write out ew, ns, lev, time, clon, clat, & scantime
      allocate (scantime(im))
      allocate (clon(im, jm))
      allocate (clat(im, jm))
      allocate (ew(im))
      allocate (ns(jm))    
      allocate (lev(km))

      call readvar1D("ev_mid_time", INV_file, scantime)
      call check(nf90_put_var(ncid,scantimeVarID,scantime), "writing out scantime")

      call readvar2D("longitude", INV_file, clon)
      call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

      call readvar2D("latitude", INV_file, clat)
      call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

      call readvar1D("lev", AER_file, lev)
      call check(nf90_put_var(ncid,levVarID,lev), "writing out lev")

      call readvar1D("ccd_pixels", INV_file, ew)
      call check(nf90_put_var(ncid,ewVarID,ew), "writing out ew")

      call readvar1D("number_of_scans", INV_file, ns)
      call check(nf90_put_var(ncid,nsVarID,ns), "writing out ns")   

      deallocate (clon)
      deallocate (clat)  
      deallocate (scantime)
      deallocate (ns)
      deallocate (ew)
      deallocate (lev)

      call check( nf90_close(ncid), "close addfile" )
  end subroutine create_addfile  

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
      Vqm(k,1,nobs) = DU001(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,2,nobs) = DU002(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,3,nobs) = DU003(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,4,nobs) = DU004(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,5,nobs) = DU005(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,6,nobs) = SS001(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,7,nobs) = SS002(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,8,nobs) = SS003(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,9,nobs) = SS004(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,10,nobs) = SS005(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,11,nobs) = BCPHOBIC(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,12,nobs) = BCPHILIC(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,13,nobs) = OCPHOBIC(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,14,nobs) = OCPHILIC(i,j,k)*DELP(i,j,k)/grav
      Vqm(k,15,nobs) = SO4(i,j,k)*DELP(i,j,k)/grav
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
  subroutine get_config(rcfile, nodenumarg, ierr)
    character(len=*),intent(in)      :: rcfile
    character(len=*),intent(in)      :: nodenumarg
    integer,intent(out)              :: ierr
    character(len=256)               :: file

    cf = ESMF_ConfigCreate()
    call ESMF_ConfigLoadFile(cf, fileName=trim(rcfile), __RC__)

    ! Read in variables
    ! -----------------------
    ! General
    call ESMF_ConfigGetAttribute(cf, date, label = 'DATE:',__RC__)
    call ESMF_ConfigGetAttribute(cf, time, label = 'TIME:',__RC__)
    call ESMF_ConfigGetAttribute(cf, instname, label = 'INSTNAME:',__RC__)
    call ESMF_ConfigGetAttribute(cf, scalar, label = 'SCALAR:',default=.TRUE.)
    call ESMF_ConfigGetAttribute(cf, plane_parallel, label = 'PLANE_PARALLEL:',default=.FALSE.)
    call ESMF_ConfigGetAttribute(cf, nstreams, label = 'NSTREAMS:',default=6)
    call ESMF_ConfigGetAttribute(cf, szamax, label = 'SZAMAX:',default=80.0)
    call ESMF_ConfigGetAttribute(cf, vzamax, label = 'VZAMAX:',default=80.0)
    call ESMF_ConfigGetAttribute(cf, additional_output, label = 'ADDITIONAL_OUTPUT:',default=.false.)
    call ESMF_ConfigGetAttribute(cf, aerosol_output, label = 'AEROSOL_OUTPUT:',default=.false.)
    call ESMF_ConfigGetAttribute(cf, cloud_output, label = 'CLOUD_OUTPUT:',default=.false.)
    call ESMF_ConfigGetAttribute(cf, nodemax, label = 'NODEMAX:',default=1) 
    call ESMF_ConfigGetAttribute(cf, version, label = 'VERSION:',default='1.0') 
    call ESMF_ConfigGetAttribute(cf, layout, label = 'LAYOUT:',default='111')    

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
    if ( (index(lower_to_upper(watername),'NOBM') > 0) ) then
      call ESMF_ConfigGetAttribute(cf, WAT_file, label = 'WAT_file:',__RC__)
    end if

    ! Clouds
    call ESMF_ConfigGetAttribute(cf, IcldTable, label = 'ICLDTABLE:',__RC__)       
    call ESMF_ConfigGetAttribute(cf, LcldTable, label = 'LCLDTABLE:',__RC__)

    ! Input Files
    call ESMF_ConfigGetAttribute(cf, AER_file, label = 'AER_file:',__RC__)
    call ESMF_ConfigGetAttribute(cf, ANG_file, label = 'ANG_file:',__RC__)
    call ESMF_ConfigGetAttribute(cf, INV_file, label = 'INV_file:',__RC__)
    if ( (index(lower_to_upper(landmodel),'RTLS') > 0) ) then
      call ESMF_ConfigGetAttribute(cf, BRDF_file, label = 'BRDF_file:',__RC__)
    end if
    if ( (index(lower_to_upper(landmodel),'RTLS-HYBRID') > 0) .or. (lower_to_upper(landmodel) == 'LAMBERTIAN') ) then
      call ESMF_ConfigGetAttribute(cf, LER_file, label = 'LER_file:',__RC__)
    end if

    call ESMF_ConfigGetAttribute(cf, OUT_file, label = 'OUT_file:',__RC__)
    if (additional_output) then
      call ESMF_ConfigGetAttribute(cf, ADD_file, label = 'ADD_file:',__RC__)
    end if 
    if (aerosol_output) then
      call ESMF_ConfigGetAttribute(cf, AERO_file, label = 'AERO_file:',__RC__)
    end if 
    if (cloud_output) then
      call ESMF_ConfigGetAttribute(cf, CLDO_file, label = 'CLDO_file:',__RC__)
    end if 

    call ESMF_ConfigGetAttribute(cf, CLD_file, label = 'CLD_file:',__RC__)
    call ESMF_ConfigGetAttribute(cf, MET_file, label = 'MET_file:',__RC__)
    
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


    ! Get the node number from the rcfile name
    ! ----------------------------------------
    ierr = 0
    if (nodemax == 1) then
      nodenum = 1
    else if (trim(nodenumarg) == '') then
      ierr = 1
      return
    else
      read(nodenumarg,*) nodenum
    end if

    if (nodenum > nodemax) then
      ierr = 2
      return
    end if

    ! Change out_file and add_file names if nodemax > 1
    if (nodemax > 1) then
      if (nodenum < 10) then
        write(file,'(A,A,I1)') trim(OUT_file),'_',nodenum
        OUT_file = file
        if (additional_output) then
          write(file,'(A,A,I1)') trim(ADD_file),'_',nodenum
          ADD_file = file 
        end if   
        if (aerosol_output) then
          write(file,'(A,A,I1)') trim(AERO_file),'_',nodenum
          AERO_file = file 
        end if  
        if (cloud_output) then
          write(file,'(A,A,I1)') trim(CLDO_file),'_',nodenum
          CLDO_file = file 
        end if                        
      else if (nodenum >= 10 .and. nodenum < 100) then
        write(file,'(A,A,I2)') trim(OUT_file),'_',nodenum
        OUT_file = file
        if (additional_output) then
          write(file,'(A,A,I2)') trim(ADD_file),'_',nodenum
          ADD_file = file  
        end if      
        if (aerosol_output) then
          write(file,'(A,A,I2)') trim(AERO_file),'_',nodenum
          AERO_file = file  
        end if        
        if (cloud_output) then
          write(file,'(A,A,I2)') trim(CLDO_file),'_',nodenum
          CLDO_file = file  
        end if        

      end if
    end if    

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
