!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    geo_vlidort
! PURPOSE
!     Reads in parallel the model data (in a netcdf file) interpolated to satellite grid 
!     The variables are needed as input to vlidort
!     A shared memory array created with a call to MAPL_ShmemMod is used for the variables
!     Do some filtering and run the vlidort code
! INPUT
!     date  : string of variable name
!     time  : file to be read
!     inst  : instrument name
!     indir : main directory for input data
!     outdir: directory for output data
! OUTPUT
!     None
!  HISTORY
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
program geo_vlidort_cloud

  use ESMF                         ! ESMF modules
  use MAPL_Mod
  use MAPL_ShmemMod                ! The SHMEM infrastructure
  use netcdf                       ! for reading the NR files
  use vlidort_brdf_modis           ! Module to run VLIDORT with MODIS BRDF surface supplement
  use vlidort_lambert              ! Module to run VLIDORT with lambertian surface  
  use lidort_brdf_modis
  use Chem_MieMod
!  use netcdf_helper                ! Module with netcdf routines
  use mp_netcdf_Mod
  use netcdf_Mod
  use GeoAngles                    ! Module with geostationary satellite algorithms for scene geometry
  use cloud_MieMod

  implicit none
  include "geo_vlidort_pars.F90"
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
  character(len=7)                      :: surfdate
  character(len=2)                      :: time 
  character(len=256)                    :: instname, angname, indir, outdir
  character(len=256)                    :: surfname, surfmodel
  character(len=256)                    :: surfband               ! flag to use nearest-neighbor interpolation or an exact value given
  integer                               :: surfbandm              ! number of wavelength bands or channels in surface reflectance data file
  integer, allocatable                  :: surfband_i(:)          ! surface band indeces that overlap with vlidort channels
  real, allocatable                     :: surfband_c(:)          ! modis band center wavelength
  logical                               :: scalar
  real, allocatable                     :: channels(:)            ! channels to simulate
  integer                               :: nch                    ! number of channels  
  real                                  :: szamax, vzamax         ! Geomtry filtering
  logical                               :: additional_output      ! does user want additional output
  integer                               :: nodemax                ! number of nodes requested
  integer                               :: nodenum                ! which node is this?
  character(len=256)                    :: version
  character(len=256)                    :: surf_version
  character(len=256)                    :: layout 
  logical                               :: vlidort
  character(len=256)                    :: IcldTable, LcldTable
  integer                               :: idxCld

! Test flag
! -----------
  logical                               :: test_shmem = .False.

! File names
! ----------
  character(len=256)                    :: AER_file, ANG_file, INV_file, SURF_file, OUT_file, ATMOS_file, CLD_file 

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
  real, pointer                         :: BAND1(:,:,:) => null()
  real, pointer                         :: BAND2(:,:,:) => null()
  real, pointer                         :: BAND3(:,:,:) => null()  
  real, pointer                         :: BAND4(:,:,:) => null()
  real, pointer                         :: BAND5(:,:,:) => null()
  real, pointer                         :: BAND6(:,:,:) => null()  
  real, pointer                         :: BAND7(:,:,:) => null()
  real, pointer                         :: BAND8(:,:,:) => null()
  real, pointer                         :: KISO(:,:,:) => null()
  real, pointer                         :: KVOL(:,:,:) => null()
  real, pointer                         :: KGEO(:,:,:) => null() 
  real, pointer                         :: LER(:,:,:) => null()     
  real*8, pointer                       :: SZA(:,:) => null()
  real*8, pointer                       :: VZA(:,:) => null()
  real*8, pointer                       :: SAA(:,:) => null()
  real*8, pointer                       :: VAA(:,:) => null()
  real*8, pointer                       :: RAA(:,:) => null()
  real, pointer                         :: REI(:,:,:) => null()  
  real, pointer                         :: REL(:,:,:) => null()  
  real, pointer                         :: TAUI(:,:,:) => null()  
  real, pointer                         :: TAUL(:,:,:) => null()  
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
  real*8, allocatable                   :: Q_int(:,:)                                 ! Q Stokes component
  real*8, allocatable                   :: U_int(:,:)                                 ! U Stokes component
  real*8, allocatable                   :: ROT_int(:,:,:)                             ! rayleigh optical thickness

!                                  Final Shared Arrays
!                                  -------------------
  real*8, pointer                       :: radiance_VL(:,:) => null()             ! TOA normalized radiance from VLIDORT
  real*8, pointer                       :: reflectance_VL(:,:) => null()          ! TOA reflectance from VLIDORT
  real*8, pointer                       :: Q(:,:) => null()                      ! Q Stokes component
  real*8, pointer                       :: U(:,:) => null()                      ! U Stokes component
  real*8, pointer                       :: ROT(:,:,:) => null()                  ! rayleigh optical thickness
  real*8, pointer                       :: ALBEDO(:,:) => null()                 ! bi-directional surface reflectance

  real, pointer                         :: TAU(:,:,:) => null()                  ! aerosol optical depth
  real, pointer                         :: SSA(:,:,:) => null()                  ! single scattering albedo
  real, pointer                         :: G(:,:,:) => null()                    ! asymmetry factor
  real, pointer                         :: PE(:,:,:) => null()

  real*8,allocatable                    :: AOD(:,:)                                 ! Temporary variable to add up AOD
! VLIDORT working variables
!------------------------------
  integer                               :: ch                                       ! i-channel  
  integer                               :: iband                                    ! i-surfaceband
  real,allocatable                      :: Vpmom(:,:,:,:,:)                         ! elements of scattering phase matrix for vector calculations
  real,allocatable                      :: VpmomIcl(:,:,:,:,:)                      ! elements of scattering phase matrix for vector calculations
  real,allocatable                      :: VpmomLcl(:,:,:,:,:)                      ! elements of scattering phase matrix for vector calculations
  real,allocatable                      :: betaIcl(:,:,:), betaLcl(:,:,:)           ! cloud optical thickness scaling factors 
  real,allocatable                      :: truncIcl(:,:,:), truncLcl(:,:,:)          ! cloud optical thickness scaling factors
! MODIS Kernel variables
!--------------------------
  real*8, allocatable                   :: kernel_wt(:,:,:)                         ! kernel weights (/fiso,fgeo,fvol/)
  real*8, allocatable                   :: param(:,:,:)                             ! Li-Sparse parameters 
                                                                                    ! param1 = crown relative height (h/b)
                                                                                    ! param2 = shape parameter (b/r)
  real                                  :: surf_missing                                                                 

! Mie Table Stucture
!---------------------
  type(Chem_Mie)                        :: mieTables

! Satellite domain variables
!------------------------------
  integer                               :: im, jm, km, tm                            ! size of satellite domain
  integer                               :: imC, jmC                                  ! size of coarse satellite domain  
  integer                               :: i, j, k, n                                ! satellite domain working variable
  integer                               :: iC, jC                                    ! coarse satellite domain working variable
  integer                               :: starti, counti, endi                      ! array indices and counts for each processor
  integer, allocatable                  :: nclr(:)                                   ! how many clear pixels each processor works on
  integer                               :: clrm                                      ! number of clear pixels for this part of decomposed domain
  integer                               :: clrm_total                                ! number of clear pixels 
  integer                               :: c, cc                                     ! clear pixel working variable
  real*8, allocatable                   :: CLDTOT(:,:)                               ! GEOS-5 cloud fraction
  real*8, allocatable                   :: FRLAND(:,:)                               ! GEOS-5 land fraction
  real*8, allocatable                   :: SOLAR_ZENITH(:,:)                         ! solar zenith angles used for data filtering
  real*8,allocatable                    :: SENSOR_ZENITH(:,:)                        ! SENSOR zenith angles used for data filtering  
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
  real                                  :: g5nr_missing                                !

! System tracking variables
! -----------------------------
  integer*8                             :: t1, t2, clock_max
  real*8                                :: clock_rate
  character(len=*), parameter           :: Iam = 'geo_vlidort_cloud'

!                               END OF VARIABLE DECLARATIONS
!----------------------------------------------------------------------------------------------------------
  
  progress = -1

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
    write(*,*) 'Input directory: ',trim(indir)
    write(*,*) 'Output directory: ',trim(outdir)
    write(*,*) 'Instrument geometry: ',trim(angname)
    write(*,*) 'BRDF dataset: ',trim(surfname),' ',trim(surfdate)
    write(*,*) 'Channels [nm]: ',channels
    if (lower_to_upper(surfband) == 'EXACT') write(*,*) 'Using exact surface relflectance parameters on bands : ',surfband_i
    if (lower_to_upper(surfband) == 'INTERPOLATE') write(*,*) 'Using interpolated surface reflectance parameters' 
    write(*,*) 'SZA < ', szamax
    write(*,*) 'VZA < ', vzamax
    if (scalar) write(*,*) 'Scalar calculations'
    if (.not. scalar) write(*,*) 'Vector calculations'
    write(*,*) 'Additional Output: ',additional_output
    if (trim(layout) /= '111') write(*,*) 'layout: ',trim(layout)
    write(*,*) ' '
  end if 

! Query for domain dimensions and missing value
!----------------------------------------------
  call mp_readDim("ew", CLD_file, im)
  call mp_readDim("ns", CLD_file, jm)
  call mp_readDim("lev", CLD_file, km)
  call mp_readDim("time", AER_file,tm)
  call mp_readDim("ew", AER_file, imC)
  call mp_readDim("ns", AER_file, jmC)

  if (lower_to_upper(surfmodel) == 'RTLS') then
    call mp_readVattr("missing_value", SURF_file, "Band1", surf_missing) 
    surf_missing = -99999.0
  else
    call mp_readVattr("missing_value", SURF_file, "SRFLER354", surf_missing) 
  end if
  call mp_readVattr("missing_value", AER_FILE, "DELP", g5nr_missing)

! Allocate arrays that will be copied on each processor - unshared
! -----------------------------------------------------------------
  call allocate_unshared()

! Create OUTFILE
! --------------
  if ( MAPL_am_I_root() )  call create_outfile(date, time)

  call MAPL_SyncSharedMemory(rc=ierr)
! Read the land and angle data 
! -------------------------------------
  call read_land()
  call read_sza()
  call read_vza()

! Create cloud-land mask
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
      if ((FRLAND(i,j) .ne. g5nr_missing) .and. (FRLAND(i,j) >= 0.99))  then
        if (SOLAR_ZENITH(i,j) < szamax) then
          if (SENSOR_ZENITH(i,j) < vzamax) then
            clrm = clrm + 1
            clmask(i,j) = .True.
          end if
        end if
      end if 
    end do
  end do

  if (clrm == 0) then
    if (MAPL_am_I_root()) then
      write(*,*) 'No good pixels over land, nothing to do'
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

! Land data and unshared angles no longer needed
!---------------------------------------
  deallocate (FRLAND)
  deallocate (SOLAR_ZENITH)
  deallocate (SENSOR_ZENITH)

  if (MAPL_am_I_root()) then
    write(*,*) '<> Created Land Mask'
    write(*,'(A,I3,A)') '       ',nint(100.*clrm/(im*jm)),'% of the domain is land and sunlit'
    write(*,'(A,I,A)')  '       Simulating ',clrm,' pixels'
    write(*,*) ' '
  end if   

! Allocate the Global arrays using SHMEM
! It will be available on all processors
! ---------------------------------------------------------
 call allocate_shared()
 call MAPL_SyncSharedMemory(rc=ierr)
  
! Read in the global arrays
! ------------------------------
 call read_aer_Nv()
 call read_surf()
 call read_angles()
 call read_cld_Tau()

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
  TAU           = dble(MISSING)
  SSA           = dble(MISSING)
  G             = dble(MISSING)
  ALBEDO        = dble(MISSING)
  ROT           = dble(MISSING)
  PE            = dble(MISSING)
  radiance_VL    = dble(MISSING)
  reflectance_VL = dble(MISSING)
  if (.not. scalar) then
    Q = dble(MISSING)
    U = dble(MISSING)
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
    iC = iCoarse(i,layout,im,imC)
    jC = jCoarse(j,layout,jm,jmC)

    call getEdgeVars ( km, nobs, reshape(AIRDENS(i,j,:),(/km,nobs/)), &
                       reshape(DELP(i,j,:),(/km,nobs/)), ptop, &
                       Vpe, Vze, Vte )   
    PE(i,j,:) = Vpe(:,nobs)

    write(msg,'(A,I)') 'getEdgeVars ', myid
    call write_verbose(msg)
    
    call calc_qm()

    write(msg,'(A,I)') 'calc_qm ', myid
    call write_verbose(msg)  

!   Surface Reflectance Parameters
!   ----------------------------------
    if ( (lower_to_upper(surfmodel) == 'RTLS') ) then
!     Get BRDF Kernel Weights
!     ----------------------------
      do ch = 1, nch
        if (lower_to_upper(surfband) == 'EXACT') then
          kernel_wt(:,ch,nobs) = (/dble(KISO(i,j,surfband_i(ch))),&
                              dble(KGEO(i,j,surfband_i(ch))),&
                              dble(KVOL(i,j,surfband_i(ch)))/)
        else
          ! > 2130 uses highest MODIS wavelength band     
          if (channels(ch) >= 2130) then  
            iband = minloc(abs(surfband_c - channels(ch)), dim = 1)
            kernel_wt(:,ch,nobs) = (/dble(KISO(i,j,iband)),&
                              dble(KGEO(i,j,iband)),&
                              dble(KVOL(i,j,iband))/)
          end if
          
          if (channels(ch) < 2130) then
            ! nearest neighbor interpolation of kernel weights to wavelength
            ! ******channel has to fall above available range, user is responsible to verify
            kernel_wt(1,ch,nobs) = dble(nn_interp(surfband_c,reshape(KISO(i,j,:),(/surfbandm/)),channels(ch)))
            kernel_wt(2,ch,nobs) = dble(nn_interp(surfband_c,reshape(KGEO(i,j,:),(/surfbandm/)),channels(ch)))
            kernel_wt(3,ch,nobs) = dble(nn_interp(surfband_c,reshape(KVOL(i,j,:),(/surfbandm/)),channels(ch)))          
          end if
        end if
        param(:,ch,nobs)     = (/dble(2),dble(1)/)
      end do
    else
      ! Use a provided albedo file
      do ch = 1, nch
        if (lower_to_upper(surfband) == 'EXACT') then
          Valbedo(nobs,ch) = dble(LER(iC,jC,surfband_i(ch)))
        else
          Valbedo(nobs,ch) = dble(nn_interp(surfband_c,reshape(LER(iC,jC,:),(/surfbandm/)),channels(ch)))
        end if
      end do
    end if 

!   Aerosol Optical Properties
!   --------------------------
    call VLIDORT_getAOPvector ( mieTables, km, nobs, nch, nq, channels, vnames, verbose, &
                        Vqm, reshape(RH(i,j,:),(/km,nobs/)),&
                        nMom,nPol, Vtau, Vssa, Vg, Vpmom, ierr )

!   Cloud Optical Properties
!   ------------------------
    call getCOPvector(IcldTable, km, nobs, nch, nMom, nPol, idxCld, REI(i,j,:), VssaIcl, VgIcl, VpmomIcl, betaIcl, truncIcl)
    call getCOPvector(LcldTable, km, nobs, nch, nMom, nPol, idxCld, REL(i,j,:), VssaLcl, VgLcl, VpmomLcl, betaLcl, truncLcl)

!   Scale Cloud Optical Thickness from reference wavelength of 0.65um
!   ------------------------
    VtauIcl(:,nch,nobs) = TAUI(i,j,:)
    VtauIcl             = VtauIcl*betaIcl
    VtauLcl(:,nch,nobs) = TAUL(i,j,:)
    VtauLcl             = VtauLcl*betaLcl

!   Scale Cloud optical thickness for phase function truncation            
!   ------------------------
    VtauIcl = VtauIcl*(1. - truncIcl*VssaIcl)
    VtauLcl = VtauLcl*(1. - truncLcl*VssaLcl)
    VssaIcl = VssaIcl*(1. - truncIcl)/(1. - VssaIcl*truncIcl)
    VssaLcl = VssaLcl*(1. - truncLcl)/(1. - VssaLcl*truncLcl)

!   Save some variables on the 2D Grid for Writing Later
!   ------------------------
    TAU(i,j,:) = Vtau(:,nch,nobs)
    SSA(i,j,:) = Vssa(:,nch,nobs)
    G(i,j,:)   = Vg(:,nch,nobs)

    write(msg,*) 'getAOP ', myid
    call write_verbose(msg)


!   Call VlIDORT
!   ------------
    if ( (lower_to_upper(surfmodel) /= 'RTLS') ) then
!     Simple lambertian surface model
!     -------------------------------
      if (scalar) then
        if (vlidort) then
          ! Call to vlidort scalar code       
          call VLIDORT_Scalar_Lambert_Cloud (km, nch, nobs ,dble(channels), nMom,      &
                  nPol, dble(Vtau), dble(Vssa), dble(Vg), dble(Vpmom),&
                  dble(VtauIcl), dble(VssaIcl), dble(VgIcl), dble(VpmomIcl),&
                  dble(VtauLcl), dble(VssaLcl), dble(VgLcl), dble(VpmomLcl),&
                  dble(Vpe), dble(Vze), dble(Vte), Valbedo,&
                  (/dble(SZA(i,j))/), &
                  (/dble(abs(RAA(i,j)))/), &
                  (/dble(VZA(i,j))/), &
                  dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT_int, ierr)
        else 
          call LIDORT_Scalar_Lambert_Cloud (km, nch, nobs ,dble(channels), nMom,      &
                  nPol, dble(Vtau), dble(Vssa), dble(Vg), dble(Vpmom),&
                  dble(VtauIcl), dble(VssaIcl), dble(VgIcl), dble(VpmomIcl),&
                  dble(VtauLcl), dble(VssaLcl), dble(VgLcl), dble(VpmomLcl),&
                  dble(Vpe), dble(Vze), dble(Vte), Valbedo,&
                  (/dble(SZA(i,j))/), &
                  (/dble(abs(RAA(i,j)))/), &
                  (/dble(VZA(i,j))/), &
                  dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT_int, ierr)
        end if 
      else
        ! Call to vlidort vector code
        call VLIDORT_Vector_Lambert_Cloud (km, nch, nobs ,dble(channels), nMom,   &
               nPol, dble(Vtau), dble(Vssa), dble(Vpmom), &
               dble(VtauIcl), dble(VssaIcl), dble(VpmomIcl), &
               dble(VtauLcl), dble(VssaLcl), dble(VpmomLcl), &
               dble(Vpe), dble(Vze), dble(Vte), Valbedo,&
               (/dble(SZA(i,j))/), &
               (/dble(abs(RAA(i,j)))/), &
               (/dble(VZA(i,j))/), &
               dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT_int, Q_int, U_int, ierr)
      end if

    else if ( ANY(kernel_wt == surf_missing) ) then
!     Save code for pixels that were not gap filled
!     ---------------------------------------------    
      radiance_VL_int(nobs,:) = -500
      reflectance_VL_int(nobs,:) = -500
      Valbedo = -500
      ROT_int = -500
      if (.not. scalar) then
        Q_int = -500
        U_int = -500
      end if
      ierr = 0
    else   
!     MODIS BRDF Surface Model
!     ------------------------------
      if (scalar) then 
        if (vlidort) then
          ! Call to vlidort scalar code            
          call VLIDORT_Scalar_LandMODIS_Cloud (km, nch, nobs, dble(channels), nMom,  &
                  nPol, dble(Vtau), dble(Vssa), dble(Vg), dble(Vpmom), &
                  dble(VtauIcl), dble(VssaIcl), dble(VgIcl), dble(VpmomIcl), &
                  dble(VtauLcl), dble(VssaLcl), dble(VgLcl), dble(VpmomLcl), &
                  dble(Vpe), dble(Vze), dble(Vte), &
                  kernel_wt, param, &
                  (/dble(SZA(i,j))/), &
                  (/dble(abs(RAA(i,j)))/), &
                  (/dble(VZA(i,j))/), &
                  dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT_int, Valbedo, ierr )  
        else
        ! Call to vlidort scalar code            
          call LIDORT_Scalar_LandMODIS_Cloud (km, nch, nobs, dble(channels), nMom,  &
                  nPol, dble(Vtau), dble(Vssa), dble(Vg), dble(Vpmom), &
                  dble(VtauIcl), dble(VssaIcl), dble(VgIcl), dble(VpmomIcl), &
                  dble(VtauLcl), dble(VssaLcl), dble(VgLcl), dble(VpmomLcl), &                  
                  dble(Vpe), dble(Vze), dble(Vte), &
                  kernel_wt, param, &
                  (/dble(SZA(i,j))/), &
                  (/dble(abs(RAA(i,j)))/), &
                  (/dble(VZA(i,j))/), &
                  dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT_int, Valbedo, ierr )  
        end if
      else
        ! Call to vlidort vector code
        call VLIDORT_Vector_LandMODIS_cloud (km, nch, nobs, dble(channels), nMom, &
                nPol, dble(Vtau), dble(Vssa), dble(Vpmom), &
                dble(VtauIcl), dble(VssaIcl), dble(VpmomIcl), &
                dble(VtauLcl), dble(VssaLcl), dble(VpmomLcl), &                
                dble(Vpe), dble(Vze), dble(Vte), &
                kernel_wt, param, &
                (/dble(SZA(i,j))/), &
                (/dble(abs(RAA(i,j)))/), &
                (/dble(VZA(i,j))/), &
                dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT_int, Valbedo, Q_int, U_int, ierr )  
      end if      
    end if          
    
!   Check VLIDORT Status, Store Outputs in Shared Arrays
!   ----------------------------------------------------    
    call mp_check_vlidort(radiance_VL_int,reflectance_VL_int)  
    radiance_VL(i,j)    = radiance_VL_int(nobs,nch)
    reflectance_VL(i,j) = reflectance_VL_int(nobs,nch)
    ALBEDO(i,j) = Valbedo(nobs,nch)

    ROT(i,j,:) = ROT(:,nobs,nch)
    
    if (.not. scalar) then
      Q(i,j)      = Q_int(nobs,nch)
      U(i,j)      = U_int(nobs,nch)
    end if
    
    write(msg,*) 'VLIDORT Calculations DONE', myid, ierr
    call write_verbose(msg)

!   Keep track of progress of each processor
!   -----------------------------------------        
    if (nint(100.*real(cc-starti)/real(counti)) > progress) then
      progress = nint(100.*real(cc-starti)/real(counti))
      write(*,'(A,I,A,I,A,I2,A,I3,A)') 'Pixel: ',cc,'  End Pixel: ',endi,'  ID:',myid,'  Progress:', nint(progress),'%'           
    end if
                
  end do ! do clear pixels

! Wait for everyone to finish calculations
! ----------------------------------------
  call MAPL_SyncSharedMemory(rc=ierr)


! Write to OUT_file and ATMOS_file
! -----------------------------------------
  if (MAPL_am_I_root()) then
   
    allocate (AOD(im,jm))

!                             Write to main OUT_File
!                             ----------------------
    call check( nf90_open(OUT_file, nf90_write, ncid), "opening file " // OUT_file )
    do ch = 1, nch
      write(msg,'(F10.2)') channels(ch)
      
      call check(nf90_inq_varid(ncid, 'ref_' // trim(adjustl(msg)), varid), "get ref vaird")
      call check(nf90_put_var(ncid, varid, reflectance_VL, &
                  start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out reflectance")

      call check(nf90_inq_varid(ncid, 'surf_ref_' // trim(adjustl(msg)), varid), "get ref vaird")
      call check(nf90_put_var(ncid, varid, ALBEDO, &
                    start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out albedo")

      AOD = 0
      do k=1,km 
        AOD = AOD + TAU(:,:,k)
      end do

      where(ALBEDO .eq. MISSING)  AOD = MISSING
      

      call check(nf90_inq_varid(ncid, 'aod_' // trim(adjustl(msg)), varid), "get aod vaird")
      call check(nf90_put_var(ncid, varid, AOD, &
                    start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out aod")
    end do
    call check( nf90_close(ncid), "close outfile" )

    if (additional_output) then
!                             Write to Additional Outputs File
!                             --------------------------------
      call check( nf90_open(ATMOS_file, nf90_write, ncid), "opening file " // ATMOS_file )
      do ch = 1, nch
        write(msg,'(F10.2)') channels(ch)
        call check(nf90_inq_varid(ncid, 'rad_' // trim(adjustl(msg)), varid), "get rad vaird")
        call check(nf90_put_var(ncid, varid, radiance_VL, &
                      start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out radiance")

        if (.not. scalar) then
          call check(nf90_inq_varid(ncid, 'q_' // trim(adjustl(msg)), varid), "get q vaird")
          call check(nf90_put_var(ncid, varid, Q, &
                      start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out Q")

          call check(nf90_inq_varid(ncid, 'u_' // trim(adjustl(msg)), varid), "get u vaird")
          call check(nf90_put_var(ncid, varid, U, &
                      start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out U")
        endif        

        do k=1,km 

          call check(nf90_inq_varid(ncid, 'aot_' // trim(adjustl(msg)), varid), "get aot vaird")
          call check(nf90_put_var(ncid, varid, TAU(:,:,k), &
                    start = (/1,1,k,nobs/), count = (/im,jm,1,nobs/)), "writing out tau")

          call check(nf90_inq_varid(ncid, 'g_' // trim(adjustl(msg)), varid), "get g vaird")
          call check(nf90_put_var(ncid, varid, G(:,:,k), &
                    start = (/1,1,k,nobs/), count = (/im,jm,1,nobs/)), "writing out g")

          call check(nf90_inq_varid(ncid, 'ssa_' // trim(adjustl(msg)), varid), "get ssa vaird")
          call check(nf90_put_var(ncid, varid, SSA(:,:,k), &
                    start = (/1,1,k,nobs/), count = (/im,jm,1,nobs/)), "writing out ssa")

          call check(nf90_inq_varid(ncid, 'rot_' // trim(adjustl(msg)), varid), "get rot vaird")
          call check(nf90_put_var(ncid, varid, ROT(:,:,k), &
                    start = (/1,1,k,nobs/), count = (/im,jm,1,nobs/)), "writing out rot")
        end do
      end do

      call check(nf90_inq_varid(ncid, 'pe', varid), "get pe vaird")
      do k=1,km+1
        call check(nf90_put_var(ncid, varid, PE(:,:,k), &
                   start = (/1,1,k,nobs/), count = (/im,jm,1,nobs/)), "writing out pe")
      end do
      call check( nf90_close(ncid), "close atmosfile" )
    end if  !additional_output
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
!  call deallocate_shared()
500  call MAPL_SyncSharedMemory(rc=ierr)
  call MAPL_FinalizeShmem (rc=ierr)
  call ESMF_Finalize(__RC__)

! -----------------------------------------------------------------------------------

  contains

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

function nn_interp(x,y,xint)
  real,intent(in),dimension(:)    :: x, y
  real,intent(in)                 :: xint
  real                            :: nn_interp
  integer                         :: below, above
  real                            :: top, bottom

  above = minloc(abs(xint - x), dim = 1, mask = (xint - x) .LT. 0)
  below = minloc(abs(xint - x), dim = 1, mask = (xint - x) .GE. 0)

  if (.not. ANY((/y(above),y(below)/) == surf_missing)) then
    top = y(above) - y(below)
    bottom = x(above) - x(below)
    nn_interp = y(below) + (xint-x(below)) * top / bottom
  else
    nn_interp  = surf_missing
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
!  call MAPL_AllocNodeArray(FROCEAN,(/im,jm/),rc=ierr)
!  call MAPL_AllocNodeArray(FRLAKE,(/im,jm/),rc=ierr)
  allocate(FRLAND(im,jm))
!  call MAPL_AllocNodeArray(FRLANDICE,(/im,jm/),rc=ierr)

!  call mp_colreadvar("FROCEAN", LAND_file, npet, myid, FROCEAN)
!  call mp_colreadvar("FRLAKE", LAND_file,  npet, myid, FRLAKE)
  call readvar2D("FRLAND", CLD_file, FRLAND)
!  call mp_colreadvar("FRLANDICE", LAND_file,  npet, myid, FRLANDICE)  

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
  allocate (SOLAR_ZENITH(im,jm))

  call readvar2D("solar_zenith", ANG_file, SOLAR_ZENITH)
end subroutine read_sza

subroutine read_vza()
  allocate (SENSOR_ZENITH(im,jm))

  call readvar2D("sensor_zenith", ANG_file, SENSOR_ZENITH)
end subroutine read_vza

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
  
  ! INFILES
  write(AER_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/', &
                          trim(instname),'-g5nr.lb2.aer_Nv.',date,'_',time,'z.nc4'

  if (trim(layout) == '111') then
    write(CLD_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/', &
                            trim(instname),'-g5nr-icacl-TOTWPDF-GCOP-SKEWT.',date,'_',time,'z.nc4'                                     
    write(ANG_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/', &
                          trim(angname),'.cloud.lb2.angles.',date,'_',time,'z.nc4'
  else
    write(CLD_file,'(16A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/', &
                            trim(instname),'-g5nr-icacl-TOTWPDF-GCOP-SKEWT.',date,'_',time,'z.',trim(layout),'.nc4'  
    write(ANG_file,'(16A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/', &
                          trim(angname),'.cloud.lb2.angles.',date,'_',time,'z.',trim(layout),'.nc4'                              
  end if  

  if ( lower_to_upper(surfmodel) == 'RTLS' ) then
    if ( lower_to_upper(surfname) == 'MAIACRTLS' ) then
      write(SURF_file,'(8A)') trim(indir),'/BRDF/v',trim(surf_version),'/',trim(surfname),'.',surfdate,'.hdf'
    else
      write(SURF_file,'(10A)') trim(indir),'/BRDF/v',trim(surf_version),'/',trim(surfname),'.',surfdate,'.',trim(layout),'.nc4'
    end if
  else
    write(SURF_file,'(6A)') trim(indir),'/SurfLER/',trim(instname),'-omi.SurfLER.',date(5:6),'.nc4'
  end if
  if (trim(layout) == '111') then
    write(INV_file,'(4A)')  trim(indir),'/LevelG/invariant/',trim(instname),'.lg1.invariant.nc4'
  else
    write(INV_file,'(6A)')  trim(indir),'/LevelG/invariant/',trim(instname),'.lg1.invariant.',trim(layout),'.nc4'
  end if

! OUTFILES
  if (vlidort) then
    write(OUT_file,'(4A)') trim(outdir),'/',trim(angname),'-g5nr.cloud.lc2.vlidort.'
  else
    write(OUT_file,'(4A)') trim(outdir),'/',trim(angname),'-g5nr.cloud.lc2.lidort.'
  end if 
  call outfile_extname(OUT_file)

  if (additional_output) then
    if (vlidort) then
      write(ATMOS_file,'(4A)') trim(outdir),'/',trim(angname),'-g5nr.cloud.add.vlidort.'
    else
      write(ATMOS_file,'(4A)') trim(outdir),'/',trim(angname),'-g5nr.cloud.add.lidort.'
    end if
    call outfile_extname(ATMOS_file)
  end if

end subroutine filenames

subroutine outfile_extname(file)
  character(len=256),intent(inout)     :: file
  character(len=256)                   :: chmax, chmin
  integer                              :: i


  if (scalar) then 
    write(file,'(2A)') trim(file),'scalar.'
  else
    write(file,'(2A)') trim(file),'vector.'
  end if 
  
  if (lower_to_upper(surfmodel) /= 'RTLS') then
    if (lower_to_upper(surfband) == 'EXACT') then
      write(file,'(2A)') trim(file),'lambertian.'
    else
      write(file,'(2A)') trim(file),'ilambertian.'
    end if
  else
    if (lower_to_upper(surfband) == 'EXACT') then
      write(file,'(3A)') trim(file),trim(surfname),'.'
    else
      write(file,'(4A)') trim(file),'i',trim(surfname),'.'
    end if 
  end if

  write(chmax,'(F10.2)') maxval(channels)
  write(chmin,'(F10.2)') minval(channels)

  i = index(chmax,'.')
  chmax(i:i) = 'd'
  i = index(chmin,'.')
  chmin(i:i) = 'd'

  if (nch == 1) then    
    write(file,'(7A)') trim(file),date,'_',time,'z_',trim(adjustl(chmin)),'nm.'
  else 
    write(file,'(9A)') trim(file),date,'_',time,'z_',trim(adjustl(chmin)),'-',trim(adjustl(chmax)),'nm.'
  end if

  if (trim(layout) /= '111' ) then
    write(file,'(3A)') trim(file),trim(layout),'.'
  end if
 
  if (nodemax > 1) then
    if (nodenum < 10) then
      write(file,'(A,I1,A)') trim(file),nodenum,'.nc4'
    else if (nodenum >= 10 .and. nodenum < 100) then
      write(file,'(A,I2,A)') trim(file),nodenum,'.nc4'
    end if
  else
    write(file,'(2A)') trim(file),'nc4'
  end if
end subroutine outfile_extname

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

    call mp_readvar3Dchunk("DU001", AER_file, (/imC,jmC,km/), 1, npet, myid, DU001) 
    call mp_readvar3Dchunk("DU002", AER_file, (/imC,jmC,km/), 1, npet, myid, DU002) 
    call mp_readvar3Dchunk("DU003", AER_file, (/imC,jmC,km/), 1, npet, myid, DU003) 
    call mp_readvar3Dchunk("DU004", AER_file, (/imC,jmC,km/), 1, npet, myid, DU004) 
    call mp_readvar3Dchunk("DU005", AER_file, (/imC,jmC,km/), 1, npet, myid, DU005) 
    call mp_readvar3Dchunk("SS001", AER_file, (/imC,jmC,km/), 1, npet, myid, SS001) 
    call mp_readvar3Dchunk("SS002", AER_file, (/imC,jmC,km/), 1, npet, myid, SS002) 
    call mp_readvar3Dchunk("SS003", AER_file, (/imC,jmC,km/), 1, npet, myid, SS003) 
    call mp_readvar3Dchunk("SS004", AER_file, (/imC,jmC,km/), 1, npet, myid, SS004) 
    call mp_readvar3Dchunk("SS005", AER_file, (/imC,jmC,km/), 1, npet, myid, SS005) 
    call mp_readvar3Dchunk("BCPHOBIC", AER_file, (/imC,jmC,km/), 1, npet, myid, BCPHOBIC) 
    call mp_readvar3Dchunk("BCPHILIC", AER_file, (/imC,jmC,km/), 1, npet, myid, BCPHILIC) 
    call mp_readvar3Dchunk("OCPHOBIC", AER_file, (/imC,jmC,km/), 1, npet, myid, OCPHOBIC) 
    call mp_readvar3Dchunk("OCPHILIC", AER_file, (/imC,jmC,km/), 1, npet, myid, OCPHILIC) 
    call mp_readvar3Dchunk("SO4", AER_file, (/imC,jmC,km/), 1, npet, myid, SO4) 

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

    call mp_readvar3Dchunk("AIRDENS", CLD_file, (/im,jm,km/), 1, npet, myid, AIRDENS)  
    call mp_readvar3Dchunk("RH"     , CLD_file, (/im,jm,km/), 1, npet, myid, RH) 
    call mp_readvar3Dchunk("DELP"   , CLD_file, (/im,jm,km/), 1, npet, myid, DELP) 
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
!     read_surf
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
  subroutine read_surf()
    real, dimension(im,jm,surfbandm)   :: temp
    integer                            :: Nx, Ny, ntile
    integer                            :: xstart, xend, ystart, yend
    integer                            :: x, y

    if (lower_to_upper(surfmodel) == 'RTLS') then
      
      call mp_readvar3Dchunk("Band1", SURF_file, (/im,jm,nkernel/), 1, npet, myid, Band1) 
      call mp_readvar3Dchunk("Band2", SURF_file, (/im,jm,nkernel/), 1, npet, myid, Band2) 
      call mp_readvar3Dchunk("Band3", SURF_file, (/im,jm,nkernel/), 1, npet, myid, Band3) 
      call mp_readvar3Dchunk("Band4", SURF_file, (/im,jm,nkernel/), 1, npet, myid, Band4) 
      call mp_readvar3Dchunk("Band5", SURF_file, (/im,jm,nkernel/), 1, npet, myid, Band5) 
      call mp_readvar3Dchunk("Band6", SURF_file, (/im,jm,nkernel/), 1, npet, myid, Band6) 
      call mp_readvar3Dchunk("Band7", SURF_file, (/im,jm,nkernel/), 1, npet, myid, Band7) 
      if (lower_to_upper(surfname) == 'MAIACRTLS') then
        call mp_readvar3Dchunk("Band8", SURF_file, (/im,jm,nkernel/), 1, npet, myid, Band8) 
      end if

      call MAPL_SyncSharedMemory(rc=ierr)    

      if (MAPL_am_I_root()) then
        KISO(:,:,1) = Band1(:,:,1)
        KISO(:,:,2) = Band2(:,:,1)
        KISO(:,:,3) = Band3(:,:,1)
        KISO(:,:,4) = Band4(:,:,1)
        KISO(:,:,5) = Band5(:,:,1)
        KISO(:,:,6) = Band6(:,:,1)
        KISO(:,:,7) = Band7(:,:,1)
        if (lower_to_upper(surfname) == 'MAIACRTLS') KISO(:,:,8) = Band8(:,:,1)

        KVOL(:,:,1) = Band1(:,:,2)
        KVOL(:,:,2) = Band2(:,:,2)
        KVOL(:,:,3) = Band3(:,:,2)
        KVOL(:,:,4) = Band4(:,:,2)
        KVOL(:,:,5) = Band5(:,:,2)
        KVOL(:,:,6) = Band6(:,:,2)
        KVOL(:,:,7) = Band7(:,:,2)
        if (lower_to_upper(surfname) == 'MAIACRTLS') KVOL(:,:,8) = Band8(:,:,2)

        KGEO(:,:,1) = Band1(:,:,3)
        KGEO(:,:,2) = Band2(:,:,3)
        KGEO(:,:,3) = Band3(:,:,3)
        KGEO(:,:,4) = Band4(:,:,3)
        KGEO(:,:,5) = Band5(:,:,3)
        KGEO(:,:,6) = Band6(:,:,3)
        KGEO(:,:,7) = Band7(:,:,3)
        if (lower_to_upper(surfname) == 'MAIACRTLS') KGEO(:,:,8) = Band8(:,:,3)

        write(*,*) '<> Read BRDF data to shared memory' 
      end if 

    else
      if (MAPL_am_I_root()) then
        call read_LER(LER)        
        write(*,*) '<> Read LER data to shared memory'
      end if
    end if

  end subroutine read_surf

  subroutine read_LER(indata)
    real, intent(inout),dimension(imC,jmC,surfbandm)         :: indata
    real, dimension(imC,jmC,1,1)                             :: temp

    call readvar4d("SRFLER354", SURF_file, temp)
    indata(:,:,1) = temp(:,:,1,1)

    call readvar4d("SRFLER388", SURF_file, temp)
    indata(:,:,2) = temp(:,:,1,1)
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
    real, allocatable          :: saa_(:,:)
    integer                    :: i,j


    call mp_readvar2Dchunk("solar_zenith",   ANG_file, (/im,jm/), 1, npet, myid, SZA) 
    call mp_readvar2Dchunk("sensor_zenith",  ANG_file, (/im,jm/), 1, npet, myid, VZA) 
    call mp_readvar2Dchunk("solar_azimuth",  ANG_file, (/im,jm/), 1, npet, myid, SAA)
    call mp_readvar2Dchunk("sensor_azimuth", ANG_file, (/im,jm/), 1, npet, myid, VAA)
    call MAPL_SyncSharedMemory(rc=ierr)    
    if (MAPL_am_I_root()) then        
      allocate (saa_(im,jm))
      
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
    ! call MAPL_DeallocNodeArray(SZA,rc=ierr) 
    ! call MAPL_DeallocNodeArray(VZA,rc=ierr) 
    ! call MAPL_DeallocNodeArray(SAA,rc=ierr) 
    ! call MAPL_DeallocNodeArray(VAA,rc=ierr) 

  end subroutine read_angles  


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
    call MAPL_AllocNodeArray(DU001,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU002,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU003,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU004,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU005,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS001,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS002,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS003,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS004,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS005,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(BCPHOBIC,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(BCPHILIC,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(OCPHOBIC,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(OCPHILIC,(/imC,jmC,km/),rc=ierr)
    call MAPL_AllocNodeArray(SO4,(/imC,jmC,km/),rc=ierr)
    if (lower_to_upper(surfmodel) == 'RTLS') then
      call MAPL_AllocNodeArray(KISO,(/im,jm,surfbandm/),rc=ierr)
      call MAPL_AllocNodeArray(KVOL,(/im,jm,surfbandm/),rc=ierr)
      call MAPL_AllocNodeArray(KGEO,(/im,jm,surfbandm/),rc=ierr)
    else
      call MAPL_AllocNodeArray(LER,(/imC,jmC,surfbandm/),rc=ierr)
    end if 

    call MAPL_AllocNodeArray(SZA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(VZA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(RAA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(SAA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(VAA,(/im,jm/),rc=ierr) 

    call MAPL_AllocNodeArray(TAU,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SSA,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(G,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(ROT,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(ALBEDO,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(PE,(/im,jm,km+1/),rc=ierr)
    
    if (.not. scalar) then
      call MAPL_AllocNodeArray(Q,(/im,jm,nch/),rc=ierr)
      call MAPL_AllocNodeArray(U,(/im,jm,nch/),rc=ierr)
    end if

    call MAPL_AllocNodeArray(radiance_VL,(/im,jm,nch/),rc=ierr)
    call MAPL_AllocNodeArray(reflectance_VL,(/im,jm,nch/),rc=ierr)

    if (lower_to_upper(surfmodel) == 'RTLS') then
      call MAPL_AllocNodeArray(BAND1,(/im,jm,nkernel/),rc=ierr)
      call MAPL_AllocNodeArray(BAND2,(/im,jm,nkernel/),rc=ierr)
      call MAPL_AllocNodeArray(BAND3,(/im,jm,nkernel/),rc=ierr)
      call MAPL_AllocNodeArray(BAND4,(/im,jm,nkernel/),rc=ierr)
      call MAPL_AllocNodeArray(BAND5,(/im,jm,nkernel/),rc=ierr)
      call MAPL_AllocNodeArray(BAND6,(/im,jm,nkernel/),rc=ierr)
      call MAPL_AllocNodeArray(BAND7,(/im,jm,nkernel/),rc=ierr)
      if (lower_to_upper(surfname) == 'MAIACRTLS') then
        call MAPL_AllocNodeArray(BAND8,(/im,jm,nkernel/),rc=ierr)
      end if
    end if

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

    if (lower_to_upper(surfmodel) == 'RTLS') then
      allocate (kernel_wt(nkernel,nch,nobs))
    end if
    allocate (param(nparam,nch,nobs))

    allocate (ROT_int(km,nobs,nch))
    allocate (Vpmom(km,nch,nobs,nMom,nPol))
    allocate (VpmomLcl(km,nch,nobs,nMom,nPol))
    allocate (VpmomIcl(km,nch,nobs,nMom,nPol))
    allocate (betaIcl(km,nch,nobs))
    allocate (betaLcl(km,nch,nobs))
    allocate (truncIcl(km,nch,nobs))
    allocate (truncLcl(km,nch,nobs))

    if (.not. scalar) then      
      allocate (Q_int(nobs, nch))
      allocate (U_int(nobs, nch))
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
    integer,dimension(nch)             :: qVarID, uVarID, albVarID      
    integer,dimension(nch)             :: ssaVarID, tauVarID, gVarID, rotVarID
    integer                            :: peVarID     
    
    integer                            :: ncid
    integer                            :: timeDimID, ewDimID, nsDimID, levDimID, chaDimID 
    integer                            :: leveDimID      
    integer                            :: scantimeVarID, clonVarID, clatVarID
    integer                            :: timeVarID, levVarID, ewVarID, nsVarID
    integer                            :: leveVarID, e
    integer                            :: ch

    real*8,allocatable,dimension(:,:)  :: clon, clat, sza, vza, raa
    real*8,allocatable,dimension(:)    :: scantime, ew, ns, tyme, lev

    character(len=2000)                :: comment

!                                MAIN LevelC2 OUT_FILE 
!                                ----------------------
    ! Open File
    call check(nf90_create(OUT_file, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // OUT_file)

    ! Create dimensions
    call check(nf90_def_dim(ncid, "time", tm, timeDimID), "creating time dimension")
    call check(nf90_def_dim(ncid, "lev", 1, levDimID), "creating ns dimension") !km
    call check(nf90_def_dim(ncid, "ew", im, ewDimID), "creating ew dimension") !im
    call check(nf90_def_dim(ncid, "ns", jm, nsDimID), "creating ns dimension") !jm

    ! Global Attributes
    write(comment,'(A)') 'VLIDORT Simulation of GEOS-5 '//lower_to_upper(trim(instname))//' Sampler'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'title',trim(comment)),"title attr")

    write(comment,'(A)') 'NASA/Goddard Space Flight Center'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'institution',trim(comment)),"institution attr")

    write(comment,'(A)') 'Global Model and Assimilation Office'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'source',trim(comment)),"source attr")

    write(comment,'(A)') 'VLIDORT simulation run from geo_vlidort_cloud.x'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'history',trim(comment)),"history attr")

    call check(nf90_put_att(ncid,NF90_GLOBAL,'grid_inputs',trim(INV_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'angle_inputs',trim(ANG_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'aerosol_inputs',trim(AER_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'surface_inputs',trim(SURF_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'cloud_inputs',trim(CLD_file)),"input files attr")

    call check(nf90_put_att(ncid,NF90_GLOBAL,'inputs',trim(comment)),"input files attr")
    write(comment,'(A)') 'n/a'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'references',trim(comment)),"references attr") 

    write(comment,'(A)') 'This file contains VLIDORT simulated surface reflectance (albedo) ' // &
                         'and top of the atmosphere ' // &
                         'radiance and reflectance from GEOS-5 parameters sampled ' // &
                         ' on the '//lower_to_upper(trim(instname))//' geostationary grid '
    call check(nf90_put_att(ncid,NF90_GLOBAL,'comment',trim(comment)),"comment attr")   

    if (lower_to_upper(surfmodel) /= 'RTLS') then
      if (lower_to_upper(surfband) == 'INTERPOLATE') then
        write(comment,'(A)') 'Lambertian surface reflectance interpolated to channel'
      else
        write(comment,'(A)') 'Lambertian surface reflectance without interpolation to channel'
      end if
    else
      if (lower_to_upper(surfband) == 'INTERPOLATE') then
        write(comment,'(2A)') lower_to_upper(trim(surfname)),' surface BRDF kernel weights interpolated to channel'
      else
        write(comment,'(2A)') lower_to_upper(trim(surfname)),' surface BRDF kernel weights without inerpolation to channel'
      end if
    end if

    call check(nf90_put_att(ncid,NF90_GLOBAL,'surface_comment',trim(comment)),"surface_comment")

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

    ! Define Variables
!                                     Dimensions
!                                     ----------    
    call check(nf90_def_var(ncid,'time',nf90_int,(/timeDimID/),timeVarID),"create time var")
    call check(nf90_def_var(ncid,'lev',nf90_float,(/levDimID/),levVarID),"create lev var")
    call check(nf90_def_var(ncid,'ew',nf90_float,(/ewDimID/),ewVarID),"create ew var")
    call check(nf90_def_var(ncid,'ns',nf90_float,(/nsDimID/),nsVarID),"create ns var")

    call check(nf90_def_var(ncid,'scanTime',nf90_float,(/ewDimID/),scantimeVarID),"create scanTime var")
    call check(nf90_def_var(ncid,'clon',nf90_float,(/ewDimID,nsDimID/),clonVarID),"create clon var")
    call check(nf90_def_var(ncid,'clat',nf90_float,(/ewDimID,nsDimID/),clatVarID),"create clat var")

!                                     Data
!                                     ----
    do ch=1,nch
      write(comment,'(F10.2)') channels(ch)
      call check(nf90_def_var(ncid, 'ref_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),refVarID(ch)),"create reflectance var")
      call check(nf90_def_var(ncid, 'surf_ref_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),albVarID(ch)),"create albedo var")
      call check(nf90_def_var(ncid, 'aod_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),aotVarID(ch)),"create aot var")

    end do

    ! Variable Attributes
!                                          Reflectance and AOD
!                                          ------------------------  
    do ch=1,size(channels)
      write(comment,'(F10.2,A)') channels(ch), ' nm TOA Reflectance'
      call check(nf90_put_att(ncid,refVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm Top of Atmosphere Reflectance'
      call check(nf90_put_att(ncid,refVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,refVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,refVarID(ch),'units','None'),"units attr")
      call check(nf90_put_att(ncid,refVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

      write(comment,'(F10.2,A)') channels(ch), ' nm Surface Reflectance'
      call check(nf90_put_att(ncid,albVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm Bi-Directional Surface Reflectance'
      call check(nf90_put_att(ncid,albVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,albVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,albVarID(ch),'units','None'),"units attr")
      call check(nf90_put_att(ncid,albVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

      write(comment,'(F10.2,A)') channels(ch), ' nm AOD'
      call check(nf90_put_att(ncid,aotVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm Aerosol Optical Depth'
      call check(nf90_put_att(ncid,aotVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,aotVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,aotVarID(ch),'units','None'),"units attr")
      call check(nf90_put_att(ncid,aotVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")      
    end do

!                                          scanTime
!                                          -------  
    call check(nf90_put_att(ncid,scantimeVarID,'long_name','Initial Time of Scan'),"long_name attr")
    call check(nf90_put_att(ncid,scantimeVarID,'units','seconds since '//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '// &
                                                  time//':00:00'),"units attr")

!                                          EW, NS, LEV, TIME
!                                          -----------------------  
    call check(nf90_put_att(ncid,ewVarID,'long_name','pseudo longitude'),"long_name attr")
    call check(nf90_put_att(ncid,ewVarID,'units','degrees_east'),"units attr")
    call check(nf90_put_att(ncid,nsVarID,'long_name','pseudo latitude'),"long_name attr")
    call check(nf90_put_att(ncid,nsVarID,'units','degrees_north'),"units attr")   
    call check(nf90_put_att(ncid,timeVarID,'long_name','Initial Time of Scan'),"long_name attr")
    call check(nf90_put_att(ncid,timeVarID,'units','seconds since '//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '// &
                                                  time//':00:00'),"units attr")   
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
    allocate (lev(1))
    allocate (tyme(tm))

    call readvar1D("scanTime", INV_file, scantime)
    call check(nf90_put_var(ncid,scantimeVarID,scantime), "writing out scantime")

    call readvar2D("clon", INV_file, clon)
    call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

    call readvar2D("clat", INV_file, clat)
    call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

    call readvar1D("time", AER_file, tyme)
    call check(nf90_put_var(ncid,timeVarID,tyme), "writing out time")

    call readvar1D("lev", AER_file, lev)
    call check(nf90_put_var(ncid,levVarID,lev), "writing out lev")

    call readvar1D("ew", INV_file, ew)
    call check(nf90_put_var(ncid,ewVarID,ew), "writing out ew")

    call readvar1D("ns", INV_file, ns)
    call check(nf90_put_var(ncid,nsVarID,ns), "writing out ns")   

    deallocate (clon)
    deallocate (clat)  
    deallocate (scantime)
    deallocate (ns)
    deallocate (ew)
    deallocate (tyme)
    deallocate (lev)

    call check( nf90_close(ncid), "close outfile" )

    if (additional_output) then
!                           ADDITIONAL OUTPUTS
!                           ------------------
      ! Open File
      call check(nf90_create(ATMOS_file, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // ATMOS_file)

      ! Create dimensions
      call check(nf90_def_dim(ncid, "time", tm, timeDimID), "creating time dimension")
      call check(nf90_def_dim(ncid, "lev", km, levDimID), "creating ns dimension") !km
      call check(nf90_def_dim(ncid, "leve", km+1, leveDimID), "creating edge level dimension") !km+1
      call check(nf90_def_dim(ncid, "ew", im, ewDimID), "creating ew dimension") !im
      call check(nf90_def_dim(ncid, "ns", jm, nsDimID), "creating ns dimension") !jm

      ! Global Attributes
      write(comment,'(A)') 'Atmospheric inputs for VLIDORT Simulation of GEOS-5 '//lower_to_upper(trim(instname))//' Sampler'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'title',trim(comment)),"title attr")

      write(comment,'(A)') 'NASA/Goddard Space Flight Center'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'institution',trim(comment)),"institution attr")

      write(comment,'(A)') 'Global Model and Assimilation Office'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'source',trim(comment)),"source attr")

      write(comment,'(A)') 'VLIDORT simulation run from geo_vlidort_cloud.x'
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

      ! Define Variables
  !                                     Dimensions
  !                                     ----------    
      call check(nf90_def_var(ncid,'time',nf90_int,(/timeDimID/),timeVarID),"create time var")
      call check(nf90_def_var(ncid,'lev',nf90_float,(/levDimID/),levVarID),"create lev var")
      call check(nf90_def_var(ncid,'leve',nf90_float,(/leveDimID/),leveVarID),"create leve var")
      call check(nf90_def_var(ncid,'ew',nf90_float,(/ewDimID/),ewVarID),"create ew var")
      call check(nf90_def_var(ncid,'ns',nf90_float,(/nsDimID/),nsVarID),"create ns var")

      call check(nf90_def_var(ncid,'scanTime',nf90_float,(/ewDimID/),scantimeVarID),"create scanTime var")
      call check(nf90_def_var(ncid,'clon',nf90_float,(/ewDimID,nsDimID/),clonVarID),"create clon var")
      call check(nf90_def_var(ncid,'clat',nf90_float,(/ewDimID,nsDimID/),clatVarID),"create clat var")

  !                                     Data
  !                                     ----
      do ch=1,nch
        write(comment,'(F10.2)') channels(ch)
        call check(nf90_def_var(ncid, 'rad_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),radVarID(ch)),"create radiance var")              
        call check(nf90_def_var(ncid, 'aot_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),tauVarID(ch)),"create aot var")      
        call check(nf90_def_var(ncid, 'rot_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),rotVarID(ch)),"create rot var")
        call check(nf90_def_var(ncid, 'ssa_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),ssaVarID(ch)),"create ssa var")
        call check(nf90_def_var(ncid, 'g_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),gVarID(ch)),"create g var")
        if (.not. scalar) then
          call check(nf90_def_var(ncid, 'q_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),qVarID(ch)),"create Q var")
          call check(nf90_def_var(ncid, 'u_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),uVarID(ch)),"create U var")
        end if        
      end do
      call check(nf90_def_var(ncid,'pe',nf90_float,(/ewDimID,nsDimID,leveDimID,timeDimID/),peVarID),"create pe var")
      ! Variable Attributes
  !                                          Additional Data
  !                                          -----------------  
      do ch=1,size(channels)
        write(comment,'(F10.2,A)') channels(ch), ' nm TOA Radiance'
        call check(nf90_put_att(ncid,radVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm Top of Atmosphere Radiance'
        call check(nf90_put_att(ncid,radVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        write(comment,'(A)') 'only l=1 level contains data'
        call check(nf90_put_att(ncid,radVarID(ch),'comment',trim(adjustl(comment))),"long_name attr")        
        call check(nf90_put_att(ncid,radVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,radVarID(ch),'units','W m-2 sr-1 nm-1'),"units attr")
        call check(nf90_put_att(ncid,radVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

        write(comment,'(F10.2,A)') channels(ch), ' nm AOT'
        call check(nf90_put_att(ncid,tauVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm layer Aerosol Optical Thickness'
        call check(nf90_put_att(ncid,tauVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,tauVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,tauVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,tauVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

        write(comment,'(F10.2,A)') channels(ch), ' nm ROT'
        call check(nf90_put_att(ncid,rotVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm layer Rayliegh Optical Thickness'
        call check(nf90_put_att(ncid,rotVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,rotVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,rotVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,rotVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

        write(comment,'(F10.2,A)') channels(ch), ' nm SSA'
        call check(nf90_put_att(ncid,ssaVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm layer Aerosol Single Scattering Albedo'
        call check(nf90_put_att(ncid,ssaVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,ssaVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,ssaVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,ssaVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")   

        write(comment,'(F10.2,A)') channels(ch), ' nm g'
        call check(nf90_put_att(ncid,gVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
        write(comment,'(F10.2,A)') channels(ch), ' nm layer Aerosol Asymmetry Parameter'
        call check(nf90_put_att(ncid,gVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
        call check(nf90_put_att(ncid,gVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
        call check(nf90_put_att(ncid,gVarID(ch),'units','none'),"units attr")
        call check(nf90_put_att(ncid,gVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")  

        if (.not. scalar) then
          write(comment,'(F10.2,A)') channels(ch), ' nm TOA Q'
          call check(nf90_put_att(ncid,qVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
          write(comment,'(F10.2,A)') channels(ch), ' nm TOA Q Component of the Stokes Vector'
          call check(nf90_put_att(ncid,qVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
          write(comment,'(A)') 'only l=1 level contains data'
          call check(nf90_put_att(ncid,qVarID(ch),'comment',trim(adjustl(comment))),"long_name attr")                
          call check(nf90_put_att(ncid,qVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
          call check(nf90_put_att(ncid,qVarID(ch),'units','W m-2 sr-1 nm-1'),"units attr")
          call check(nf90_put_att(ncid,qVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

          write(comment,'(F10.2,A)') channels(ch), ' nm TOA U'
          call check(nf90_put_att(ncid,uVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
          write(comment,'(F10.2,A)') channels(ch), ' nm TOA U Component of the Stokes Vector'
          call check(nf90_put_att(ncid,uVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
          write(comment,'(A)') 'only l=1 level contains data'
          call check(nf90_put_att(ncid,uVarID(ch),'comment',trim(adjustl(comment))),"long_name attr")                          
          call check(nf90_put_att(ncid,uVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
          call check(nf90_put_att(ncid,uVarID(ch),'units','W m-2 sr-1 nm-1'),"units attr")
          call check(nf90_put_att(ncid,uVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")
        end if

      end do

      call check(nf90_put_att(ncid,peVarID,'standard_name','Edge Pressure'),"standard_name attr")
      call check(nf90_put_att(ncid,peVarID,'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,peVarID,'units','Pa'),"units attr")
      call check(nf90_put_att(ncid,peVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")  

  !                                          scanTime
  !                                          -------  
      call check(nf90_put_att(ncid,scantimeVarID,'long_name','Initial Time of Scan'),"long_name attr")
      call check(nf90_put_att(ncid,scantimeVarID,'units','seconds since '//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '// &
                                                    time//':00:00'),"units attr")

  !                                          EW, NS, LEV, TIME
  !                                          -----------------------  
      call check(nf90_put_att(ncid,ewVarID,'long_name','pseudo longitude'),"long_name attr")
      call check(nf90_put_att(ncid,ewVarID,'units','degrees_east'),"units attr")
      call check(nf90_put_att(ncid,nsVarID,'long_name','pseudo latitude'),"long_name attr")
      call check(nf90_put_att(ncid,nsVarID,'units','degrees_north'),"units attr")   
      call check(nf90_put_att(ncid,timeVarID,'long_name','Initial Time of Scan'),"long_name attr")
      call check(nf90_put_att(ncid,timeVarID,'units','seconds since '//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '// &
                                                    time//':00:00'),"units attr")   
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
      allocate (tyme(tm))

      call readvar1D("scanTime", INV_file, scantime)
      call check(nf90_put_var(ncid,scantimeVarID,scantime), "writing out scantime")

      call readvar2D("clon", INV_file, clon)
      call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

      call readvar2D("clat", INV_file, clat)
      call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

      call readvar1D("time", AER_file, tyme)
      call check(nf90_put_var(ncid,timeVarID,tyme), "writing out time")

      call readvar1D("lev", AER_file, lev)
      call check(nf90_put_var(ncid,levVarID,lev), "writing out lev")

      call check(nf90_put_var(ncid,leveVarID,(/(real(e), e = 1, km+1)/)), "writing out leve")      

      call readvar1D("ew", INV_file, ew)
      call check(nf90_put_var(ncid,ewVarID,ew), "writing out ew")

      call readvar1D("ns", INV_file, ns)
      call check(nf90_put_var(ncid,nsVarID,ns), "writing out ns")   

      deallocate (clon)
      deallocate (clat)  
      deallocate (scantime)
      deallocate (ns)
      deallocate (ew)
      deallocate (tyme)
      deallocate (lev)

      call check( nf90_close(ncid), "close atmosfile" )
    end if ! additional_output
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
      Vqm(k,1,nobs) = DU001(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,2,nobs) = DU002(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,3,nobs) = DU003(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,4,nobs) = DU004(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,5,nobs) = DU005(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,6,nobs) = SS001(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,7,nobs) = SS002(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,8,nobs) = SS003(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,9,nobs) = SS004(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,10,nobs) = SS005(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,11,nobs) = BCPHOBIC(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,12,nobs) = BCPHILIC(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,13,nobs) = OCPHOBIC(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,14,nobs) = OCPHILIC(iC,jC,k)*DELP(i,j,k)/grav
      Vqm(k,15,nobs) = SO4(iC,jC,k)*DELP(i,j,k)/grav
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
!    shmem_test3D
! PURPOSE
!     test that shared memory arrays have same values across all processors
! INPUT
!     varname: string of variable name
!     var    : the variable to be checked
! OUTPUT
!     Writes to the file shmem_test.txt the min and max value of the variable as 
!     reported by each processor
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine do_testing()

    if ( MAPL_am_I_root() ) then
      open (unit = 2, file="shmem_test.txt")
      write(2,'(A)') '--- Array Statistics ---'
    end if

    call shmem_test3D('RH',RH)
    call shmem_test3D('AIRDENS',AIRDENS)
    call shmem_test3D('DELP',DELP)
    call shmem_test3D('DU001',DU001)
    call shmem_test3D('DU002',DU002)
    call shmem_test3D('DU003',DU003)
    call shmem_test3D('DU004',DU004)
    call shmem_test3D('DU005',DU005)
    call shmem_test3D('SS001',SS001)
    call shmem_test3D('SS002',SS002)
    call shmem_test3D('SS003',SS003)
    call shmem_test3D('SS004',SS004)
    call shmem_test3D('SS005',SS005)
    call shmem_test3D('BCPHOBIC',BCPHOBIC)
    call shmem_test3D('BCPHILIC',BCPHILIC)
    call shmem_test3D('OCPHOBIC',OCPHOBIC)
    call shmem_test3D('OCPHILIC',OCPHILIC)
    call shmem_test3D('SO4',SO4)

    !   Wait for everyone to finish and print max memory used
    !   -----------------------------------------------------------  
    call MAPL_SyncSharedMemory(rc=ierr)
    if (MAPL_am_I_root()) then  
      write(*,*) 'Tested shared memory' 
      call sys_tracker()   
    end if   

  end subroutine do_testing

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    shmem_test3D
! PURPOSE
!     test that shared memory arrays have same values across all processors
! INPUT
!     varname: string of variable name
!     var    : the variable to be checked
! OUTPUT
!     Writes to the file shmem_test.txt the min and max value of the variable as 
!     reported by each processor
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine shmem_test3D(varname,var)
    character(len=*), intent(in)  :: varname
    real,dimension(:,:,:)         :: var
    character(len=61)             :: msg

    if ( MAPL_am_I_root() ) then
      open (unit = 2, file="shmem_test.txt",position="append")
      do p = 1,npet-1
        call mpi_recv(msg, 61, MPI_CHARACTER, p, 1, MPI_COMM_WORLD, status_mpi, ierr)
        write(2,*) msg
      end do
      write(msg,'(A9,I4,E24.17,E24.17)') varname,myid,maxval(var),minval(var)
      write(2,*) msg
      write(2,*) 'These should all have the same min/max values!'
      close(2)
    else
      write(msg,'(A9,I4,E24.17,E24.17)') varname,myid,maxval(var),minval(var)
      call mpi_send(msg, 61, MPI_CHARACTER, 0, 1, MPI_COMM_WORLD, ierr)
    end if

  end subroutine shmem_test3D

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    shmem_test2D
! PURPOSE
!     test that shared memory arrays have same values across all processors
! INPUT
!     varname: string of variable name
!     var    : the variable to be checked
! OUTPUT
!     Writes to the file shmem_test.txt the min and max value of the variable as 
!     reported by each processor
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine shmem_test2D(varname,var)
    character(len=*), intent(in)  :: varname
    real,dimension(:,:)      :: var

    if ( MAPL_am_I_root() ) then
      open (unit = 2, file="shmem_test.txt",position="append")
      do p = 1,npet-1
        call mpi_recv(msg, 61, MPI_CHARACTER, p, 1, MPI_COMM_WORLD, status_mpi, ierr)
        write(2,*) msg
      end do
        write(msg,'(A9,I4,E24.17,E24.17)') varname,myid,maxval(var),minval(var)
        write(2,*) msg
        write(2,*) 'These should all have the same min/max values!'
        close(2)
    else
        write(msg,'(A9,I4,E24.17,E24.17)') varname,myid,maxval(var),minval(var)
        call mpi_send(msg, 61, MPI_CHARACTER, 0, 1, MPI_COMM_WORLD, ierr)
    end if

  end subroutine shmem_test2D

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

    cf = ESMF_ConfigCreate()
    call ESMF_ConfigLoadFile(cf, fileName=trim(rcfile), __RC__)

    ! Read in variables
    ! -----------------------
    call ESMF_ConfigGetAttribute(cf, date, label = 'DATE:',__RC__)
    call ESMF_ConfigGetAttribute(cf, time, label = 'TIME:',__RC__)
    call ESMF_ConfigGetAttribute(cf, instname, label = 'INSTNAME:',__RC__)
    call ESMF_ConfigGetAttribute(cf, angname, label = 'ANGNAME:',default='NONE')
    call ESMF_ConfigGetAttribute(cf, indir, label = 'INDIR:',__RC__)
    call ESMF_ConfigGetAttribute(cf, outdir, label = 'OUTDIR:',default=indir)
    call ESMF_ConfigGetAttribute(cf, surfname, label = 'SURFNAME:',default='MAIACRTLS')
    call ESMF_ConfigGetAttribute(cf, surfmodel, label = 'SURFMODEL:',default='RTLS')
    call ESMF_ConfigGetAttribute(cf, surfdate, label = 'SURFDATE:',__RC__)
    call ESMF_ConfigGetAttribute(cf, scalar, label = 'SCALAR:',default=.TRUE.)
    call ESMF_ConfigGetAttribute(cf, szamax, label = 'SZAMAX:',default=80.0)
    call ESMF_ConfigGetAttribute(cf, vzamax, label = 'VZAMAX:',default=80.0)
    call ESMF_ConfigGetAttribute(cf, surfband, label = 'SURFBAND:', default='INTERPOLATE')
    call ESMF_ConfigGetAttribute(cf, surfbandm, label = 'SURFBANDM:',__RC__)
    call ESMF_ConfigGetAttribute(cf, additional_output, label = 'ADDITIONAL_OUTPUT:',default=.false.)
    call ESMF_ConfigGetAttribute(cf, nodemax, label = 'NODEMAX:',default=1) 
    call ESMF_ConfigGetAttribute(cf, version, label = 'VERSION:',default='1.0') 
    call ESMF_ConfigGetAttribute(cf, surf_version, label = 'SURF_VERSION:',default='1.0')    
    call ESMF_ConfigGetAttribute(cf, layout, label = 'LAYOUT:',default='111')    
    call ESMF_ConfigGetAttribute(cf, vlidort, label = 'VLIDORT:',default=.true.) 
    call ESMF_ConfigGetAttribute(cf, IcldTable, label = 'ICLDTABLE:',__RC__)       
    call ESMF_ConfigGetAttribute(cf, LcldTable, label = 'LCLDTABLE:',__RC__)
    call ESMF_ConfigGetAttribute(cf, idxCld, label = 'IDXCLD:',__RC__)

  ! angname is instname if not given in the rcfile
  !------------------------------------------
    if (trim(lower_to_upper(angname)) == 'NONE') then
      angname = instname
    end if

    ! Check that LER configuration is correct
    !----------------------------------------
    if (lower_to_upper(surfmodel) /= 'RTLS' ) then
      if (surfbandm /= 2) then
        surfbandm = 2
        if (MAPL_am_I_root()) then
          write(*,*) 'Wrong number of surface bands.  If not RTLS surface, SURFBANDM must equal 2'
          write(*,*) 'Forcing surfbandm = 2'          
        end if
      end if
    end if

    ! Figure out number of channels and read into vector
    !------------------------------------------------------
    nch =  ESMF_ConfigGetLen(cf, label = 'CHANNELS:',__RC__)
    allocate (channels(nch))
    call ESMF_ConfigGetAttribute(cf, channels, label = 'CHANNELS:', default=550.)

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


    ! INFILES & OUTFILE set names
    ! -------------------------------
    call filenames()

    if (lower_to_upper(surfband) == 'EXACT' ) then
      allocate (surfband_i(nch))
      call ESMF_ConfigGetAttribute(cf, surfband_i, label = 'SURFBAND_I:', __RC__)
    else
      allocate (surfband_c(surfbandm))
      call ESMF_ConfigGetAttribute(cf, surfband_c, label = 'SURFBAND_C:', __RC__)
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

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    iCoarse
! PURPOSE
!     Finds i-index of high-res tiled cloud data on coarse domain
! INPUT
!     layout : layout of tiles
!     i      : index on high-res tile
! OUTPUT
!     iC     : coarse domain indices
!  HISTORY
!     August 2017 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  function iCoarse(i,layout,im,imC)
    integer, intent(in)               :: i, im, imC
    character(len=*)                  :: layout 
    integer                           :: iCoarse
    integer                           :: Nx, ntile
    integer                           :: i_, ifactor

    read(layout(1:1),*) Nx
    read(layout(3:) ,*) ntile

    ifactor = Nx*im/imC

    i_ = mod(ntile,nX)

    iCoarse = int((i + im*i_ - 1)/ifactor) + 1

  end function iCoarse

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    jCoarse
! PURPOSE
!     Finds j-index of high-res tiled cloud data on coarse domain
! INPUT
!     layout : layout of tiles
!     j      : index on high-res tile
! OUTPUT
!     jC     : coarse domain indices
!  HISTORY
!     August 2017 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  function jCoarse(j,layout,jm,jmC)
    integer, intent(in)               :: j, jm, jmC
    character(len=*)                  :: layout 
    integer                           :: jCoarse
    integer                           :: Nx, Ny, ntile
    integer                           :: j_, jfactor

    read(layout(1:1),*) Nx
    read(layout(2:2),*) Ny
    read(layout(3:) ,*) ntile

    jfactor = Ny*jm/jmC

    j_ = int(ntile/nX) 

    jCoarse = int((j + jm*j_ - 1)/jfactor) + 1

  end function jCoarse

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    deallocate_shared
! PURPOSE
!     Clean up allocated arrays 
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine deallocate_shared()         

    ! shmem must deallocate shared memory arrays
    call MAPL_DeallocNodeArray(AIRDENS,rc=ierr)
    call MAPL_DeallocNodeArray(RH,rc=ierr)
    call MAPL_DeallocNodeArray(DELP,rc=ierr)
    call MAPL_DeallocNodeArray(DU001,rc=ierr)
    call MAPL_DeallocNodeArray(DU002,rc=ierr)
    call MAPL_DeallocNodeArray(DU003,rc=ierr)
    call MAPL_DeallocNodeArray(DU004,rc=ierr)                           
    call MAPL_DeallocNodeArray(DU005,rc=ierr)
    call MAPL_DeallocNodeArray(SS001,rc=ierr) 
    call MAPL_DeallocNodeArray(SS002,rc=ierr) 
    call MAPL_DeallocNodeArray(SS003,rc=ierr) 
    call MAPL_DeallocNodeArray(SS004,rc=ierr) 
    call MAPL_DeallocNodeArray(SS005,rc=ierr) 
    call MAPL_DeallocNodeArray(BCPHOBIC,rc=ierr) 
    call MAPL_DeallocNodeArray(BCPHILIC,rc=ierr) 
    call MAPL_DeallocNodeArray(OCPHOBIC,rc=ierr) 
    call MAPL_DeallocNodeArray(OCPHILIC,rc=ierr) 
    call MAPL_DeallocNodeArray(SO4,rc=ierr) 
    if (lower_to_upper(surfmodel) == 'RTLS') then
      call MAPL_DeallocNodeArray(KISO,rc=ierr) 
      call MAPL_DeallocNodeArray(KVOL,rc=ierr) 
      call MAPL_DeallocNodeArray(KGEO,rc=ierr) 
    else
      call MAPL_DeallocNodeArray(LER,rc=ierr) 
    end if
    call MAPL_DeallocNodeArray(SZA,rc=ierr) 
    call MAPL_DeallocNodeArray(VZA,rc=ierr) 
    call MAPL_DeallocNodeArray(RAA,rc=ierr) 

    call MAPL_DeallocNodeArray(TAU,rc=ierr)
    call MAPL_DeallocNodeArray(SSA,rc=ierr)
    call MAPL_DeallocNodeArray(G,rc=ierr)
    call MAPL_DeallocNodeArray(ROT,rc=ierr)
    call MAPL_DeallocNodeArray(ALBEDO,rc=ierr)
    call MAPL_DeallocNodeArray(PE,rc=ierr)
    if (.not. scalar) then
      call MAPL_DeallocNodeArray(Q,rc=ierr)
      call MAPL_DeallocNodeArray(U,rc=ierr)
    end if

    call MAPL_DeallocNodeArray(radiance_VL,rc=ierr) 
    call MAPL_DeallocNodeArray(reflectance_VL,rc=ierr) 

  end subroutine deallocate_shared


end program geo_vlidort_cloud
