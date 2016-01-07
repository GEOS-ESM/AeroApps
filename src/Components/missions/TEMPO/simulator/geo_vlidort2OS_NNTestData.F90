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
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#  include "MAPL_Generic.h"
#  include "MAPL_ErrLogMain.h"
program geo_vlidort2OS_NNTestData

  use ESMF                         ! ESMF modules
  use MAPL_Mod
  use MAPL_ShmemMod                ! The SHMEM infrastructure
  use netcdf                       ! for reading the NR files
  use vlidort_brdf_modis           ! Module to run VLIDORT with MODIS BRDF surface supplement
  use lidort_brdf_modis
  use Chem_MieMod
!  use netcdf_helper                ! Module with netcdf routines
  use mp_netcdf_Mod
  use netcdf_Mod
  use GeoAngles                    ! Module with geostationary satellite algorithms for scene geometry

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
  character(len=256)                    :: instname, indir, outdir
  character(len=256)                    :: surfname, surfmodel
  character(len=256)                    :: surfband               ! flag to use nearest-neighbor interpolation or an exact value given
  integer                               :: surfbandm              ! number of wavelength bands or channels in surface reflectance data file
  integer, allocatable                  :: surfband_i(:)          ! surface band indeces that overlap with vlidort channels
  real, allocatable                     :: surfband_c(:)          ! modis band center wavelength
  logical                               :: scalar
  real, allocatable                     :: channels(:)            ! channels to simulate
  integer                               :: nch                    ! number of channels  
  real                                  :: cldmax                 ! Cloud Filtering  
  real                                  :: szamax, vzamax         ! Geomtry filtering
  logical                               :: additional_output      ! does user want additional output
  integer                               :: nodemax                ! number of nodes requested
  integer                               :: nodenum                ! which node is this?
  character(len=256)                    :: version
  character(len=256)                    :: surf_version
  character(len=256)                    :: layout 
  logical                               :: vlidort
  logical                               :: DO_2OS_CORRECTION

! Test flag
! -----------
  logical                               :: test_shmem = .False.

! File names
! ----------
  character(len=256)                    :: MET_file, AER_file, ANG_file, INV_file, SURF_file, LAND_file, OUT_file, ATMOS_file 

! Global, 3D inputs to be allocated using SHMEM
! ---------------------------------------------
  real, pointer                         :: AIRDENS(:,:) => null()
  real, pointer                         :: RH(:,:) => null()
  real, pointer                         :: DELP(:,:) => null()
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
  real*8, pointer                       :: SZA(:) => null()
  real*8, pointer                       :: VZA(:) => null()
  real*8, pointer                       :: RAA(:) => null()  

  real, pointer                         :: AIRDENS_(:,:,:) => null()
  real, pointer                         :: RH_(:,:,:) => null()
  real, pointer                         :: DELP_(:,:,:) => null()
  real, pointer                         :: DU001_(:,:,:) => null()
  real, pointer                         :: DU002_(:,:,:) => null()
  real, pointer                         :: DU003_(:,:,:) => null()
  real, pointer                         :: DU004_(:,:,:) => null()
  real, pointer                         :: DU005_(:,:,:) => null()
  real, pointer                         :: SS001_(:,:,:) => null()
  real, pointer                         :: SS002_(:,:,:) => null()
  real, pointer                         :: SS003_(:,:,:) => null()
  real, pointer                         :: SS004_(:,:,:) => null()
  real, pointer                         :: SS005_(:,:,:) => null()
  real, pointer                         :: BCPHOBIC_(:,:,:) => null()
  real, pointer                         :: BCPHILIC_(:,:,:) => null()
  real, pointer                         :: OCPHOBIC_(:,:,:) => null()
  real, pointer                         :: OCPHILIC_(:,:,:) => null()
  real, pointer                         :: SO4_(:,:,:) => null()  
  real, pointer                         :: BAND1(:,:,:) => null()
  real, pointer                         :: BAND2(:,:,:) => null()
  real, pointer                         :: BAND3(:,:,:) => null()  
  real, pointer                         :: BAND4(:,:,:) => null()
  real, pointer                         :: BAND5(:,:,:) => null()
  real, pointer                         :: BAND6(:,:,:) => null()  
  real, pointer                         :: BAND7(:,:,:) => null()
  real, pointer                         :: BAND8(:,:,:) => null()
  real, pointer                         :: KISO_(:,:,:) => null()
  real, pointer                         :: KVOL_(:,:,:) => null()
  real, pointer                         :: KGEO_(:,:,:) => null()    
  real*8, pointer                       :: SZA_(:,:) => null()
  real*8, pointer                       :: VZA_(:,:) => null()
  real*8, pointer                       :: SAA_(:,:) => null()
  real*8, pointer                       :: VAA_(:,:) => null()
  integer, pointer                      :: indices(:) => null()


! VLIDORT input arrays
! ---------------------------
  real, allocatable                     :: pe(:,:)                ! edge pressure [Pa]
  real, allocatable                     :: ze(:,:)                ! edge height above sfc [m]
  real, allocatable                     :: te(:,:)                ! edge Temperature [K]
  real, allocatable                     :: qm(:,:,:)              ! (mixing ratio) * delp/g
  real, allocatable                     :: tau(:,:,:)             ! aerosol optical depth
  real, allocatable                     :: ssa(:,:,:)             ! single scattering albedo
  real, allocatable                     :: g(:,:,:)               ! asymmetry factor
  real*8, allocatable                   :: albedo(:,:)            ! surface albedo

! VLIDORT output arrays
!-------------------------------
!                                  Intermediate Unshared Arrays
!                                  -----------------------------
  real*8, allocatable                   :: radiance_VL_int(:,:)                   ! TOA normalized radiance from VLIDORT
  real*8, allocatable                   :: reflectance_VL_int(:,:)                ! TOA reflectance from VLIDORT  
  real*8, allocatable                   :: Q(:,:)                                 ! Q Stokes component
  real*8, allocatable                   :: U(:,:)                                 ! U Stokes component
  real*8, allocatable                   :: ROT(:,:,:)                             ! rayleigh optical thickness

!                                  Final Shared Arrays
!                                  -------------------
  real*8, pointer                       :: radiance_VL(:,:) => null()             ! TOA normalized radiance from VLIDORT
  real*8, pointer                       :: reflectance_VL(:,:) => null()          ! TOA reflectance from VLIDORT
  real*8, pointer                       :: Q_(:,:) => null()                      ! Q Stokes component
  real*8, pointer                       :: U_(:,:) => null()                      ! U Stokes component
  real*8, pointer                       :: ROT_(:,:,:) => null()                  ! rayleigh optical thickness
  real*8, pointer                       :: ALBEDO_(:,:) => null()                 ! bi-directional surface reflectance

  real, pointer                         :: TAU_(:,:,:) => null()                  ! aerosol optical depth
  real, pointer                         :: SSA_(:,:,:) => null()                  ! single scattering albedo
  real, pointer                         :: G_(:,:,:) => null()                    ! asymmetry factor

  real*8,allocatable                    :: field(:,:)                             ! Template for unpacking shared arrays
  real*8,allocatable                    :: AOD(:)                                 ! Temporary variable to add up AOD
! VLIDORT working variables
!------------------------------
  integer                               :: ch                                       ! i-channel  
  integer                               :: iband                                    ! i-surfaceband
  real,allocatable                      :: pmom(:,:,:,:,:)                          ! elements of scattering phase matrix for vector calculations

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
  integer                               :: i, j, k, n                                ! satellite domain working variable
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
    write(*,*) 'BRDF dataset: ',trim(surfname),' ',trim(surfdate)
    write(*,*) 'Channels [nm]: ',channels
    if (lower_to_upper(surfband) == 'EXACT') write(*,*) 'Using exact surface relflectance parameters on bands : ',surfband_i
    if (lower_to_upper(surfband) == 'INTERPOLATE') write(*,*) 'Using interpolated surface reflectance parameters' 
    write(*,*) 'Cloud Fraction <= ', cldmax
    write(*,*) 'SZA < ', szamax
    write(*,*) 'VZA < ', vzamax
    if (DO_2OS_CORRECTION) then
      write(*,*) 'Scalar + 2OS Correction Calculations'
    else
      if (scalar) write(*,*) 'Scalar calculations'
      if (.not. scalar) write(*,*) 'Vector calculations'
    end if
    write(*,*) 'Additional Output: ',additional_output
    if (trim(layout) /= '111') write(*,*) 'layout: ',trim(layout)
    write(*,*) ' '
  end if 

! Query for domain dimensions and missing value
!----------------------------------------------
  call mp_readDim("ew", MET_file, im)
  call mp_readDim("ns", MET_file, jm)
  call mp_readDim("lev", MET_file, km)
  call mp_readDim("time", MET_file,tm)
  if (lower_to_upper(surfmodel) == 'RTLS') then
    call mp_readVattr("missing_value", SURF_file, "Band1", surf_missing) 
    surf_missing = -99999.0
  else
    call mp_readVattr("missing_value", SURF_file, "SRFLER354", surf_missing) 
  end if
  call mp_readVattr("missing_value", MET_FILE, "CLDTOT", g5nr_missing)

! Allocate arrays that will be copied on each processor - unshared
! -----------------------------------------------------------------
  call allocate_unshared()

! Create OUTFILE
! --------------
  if ( MAPL_am_I_root() )  call create_outfile(date, time)

  call MAPL_SyncSharedMemory(rc=ierr)
! Read the cloud, land, and angle data 
! -------------------------------------
  call read_cloud()
  call read_land()
  call read_sza()
  call read_vza()

! Create cloud-land mask
! ------------------------
  ! first check that you can do anything!!!!!  
  if (.not. ANY(CLDTOT <= cldmax)) then   
    if (MAPL_am_I_root()) then
      write(*,*) 'The domain is too cloudy, nothing to do'
      write(*,*) 'Exiting.....'
    end if
    GOTO 500
  end if

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
  do i=1,im
    do j=1,jm
      if ((FRLAND(i,j) .ne. g5nr_missing) .and. (FRLAND(i,j) >= 0.99))  then
        if (CLDTOT(i,j) <= cldmax) then
          if (SOLAR_ZENITH(i,j) < szamax) then
            if (SENSOR_ZENITH(i,j) < vzamax) then
              clrm = clrm + 1
              clmask(i,j) = .True.
            end if
          end if
        end if
      end if 
    end do
  end do

  if (clrm == 0) then
    if (MAPL_am_I_root()) then
      write(*,*) 'No clear pixels over land, nothing to do'
      write(*,*) 'Exiting.....'
    end if
    GOTO 500
  end if

! Cloud and Land data no longer needed
!---------------------------------------
  deallocate (CLDTOT)
  deallocate (FRLAND)
  deallocate (SOLAR_ZENITH)
  deallocate (SENSOR_ZENITH)

  if (MAPL_am_I_root()) then
    write(*,*) '<> Created Cloud Mask'
    write(*,'(A,I3,A)') '       ',nint(100.*clrm/(im*jm)),'% of the domain is clear and sunlit'
    write(*,'(A,I,A)')  '       Simulating ',clrm,' pixels'
    write(*,*) ' '
  end if   

! Allocate the Global arrays using SHMEM
! It will be available on all processors
! ---------------------------------------------------------
 call allocate_shared(clrm)
 call MAPL_SyncSharedMemory(rc=ierr)
  
! Read in the global arrays
! ------------------------------
 call read_aer_Nv()
 call read_surf()
 call read_angles()

! Wait for everyone to finish reading and print max memory used
! ******This needs to loop through all processors....
! ------------------------------------------------------------------  
  call MAPL_SyncSharedMemory(rc=ierr)
   
  if (MAPL_am_I_root()) then 
    write(*,*) 'Read all variables' 
    call sys_tracker()   
    write(*,*) ' '
    call system_clock ( t2, clock_rate, clock_max )
    write ( *, * ) 'Elapsed real time = ', real ( t2 - t1 ) / real ( clock_rate )
    write(*,*) ' '
  end if   

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
  if (npet >= clrm) then
    nclr(1:clrm) = 1
  else if (npet < clrm) then
    nclr(1:npet) = clrm/npet
    nclr(npet)   = nclr(npet) + mod(clrm,npet)
  end if 

! Although read on individual PEs, all shared variables should have the same
! data in all PEs. Let's verify that.
! ----------------------------------------------------------- 
  if (test_shmem) call do_testing()

! Initialize outputs to be safe
! -------------------------------
  TAU_           = dble(MISSING)
  SSA_           = dble(MISSING)
  G_             = dble(MISSING)
  ALBEDO_        = dble(MISSING)
  ROT_           = dble(MISSING)
  radiance_VL    = dble(MISSING)
  reflectance_VL = dble(MISSING)
  if (.not. scalar) then
    Q_ = dble(MISSING)
    U_ = dble(MISSING)
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
    !starti = (clrm_total/nodemax)*(nodenum-1) +  1
    starti  = 1
  else
    !starti = sum(nclr(1:myid)) + (clrm_total/nodemax)*(nodenum-1) + 1
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

    call getEdgeVars ( km, nobs, reshape(AIRDENS(c,:),(/km,nobs/)), &
                       reshape(DELP(c,:),(/km,nobs/)), ptop, &
                       pe, ze, te )   

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
          kernel_wt(:,ch,nobs) = (/dble(KISO(c,surfband_i(ch))),&
                              dble(KGEO(c,surfband_i(ch))),&
                              dble(KVOL(c,surfband_i(ch)))/)
        else
          ! > 2130 uses highest MODIS wavelength band     
          if (channels(ch) >= 2130) then  
            iband = minloc(abs(surfband_c - channels(ch)), dim = 1)
            kernel_wt(:,ch,nobs) = (/dble(KISO(c,iband)),&
                              dble(KGEO(c,iband)),&
                              dble(KVOL(c,iband))/)
          end if
          
          if (channels(ch) < 2130) then
            ! nearest neighbor interpolation of kernel weights to wavelength
            ! ******channel has to fall above available range, user is responsible to verify
            kernel_wt(1,ch,nobs) = dble(nn_interp(surfband_c,reshape(KISO(c,:),(/surfbandm/)),channels(ch)))
            kernel_wt(2,ch,nobs) = dble(nn_interp(surfband_c,reshape(KGEO(c,:),(/surfbandm/)),channels(ch)))
            kernel_wt(3,ch,nobs) = dble(nn_interp(surfband_c,reshape(KVOL(c,:),(/surfbandm/)),channels(ch)))          
          end if
        end if
        param(:,ch,nobs)     = (/dble(2),dble(1)/)
      end do
    else
      ! Use a provided albedo file
      do ch = 1, nch
        if (lower_to_upper(surfband) == 'EXACT') then
          albedo(nobs,ch) = dble(LER(c,surfband_i(ch)))
        else
          albedo(nobs,ch) = dble(nn_interp(surfband_c,reshape(LER(c,:),(/surfbandm/)),channels(ch)))
        end if
      end do
    end if 

!   Aerosol Optical Properties
!   --------------------------
    ! if (scalar) then
    !   call VLIDORT_getAOPscalar ( mieTables, km, nobs, nch, nq, channels, vnames, verbose, &
    !                       qm, reshape(RH(c,:),(/km,nobs/)), &
    !                       tau, ssa, g, ierr )
    ! else
      call VLIDORT_getAOPvector ( mieTables, km, nobs, nch, nq, channels, vnames, verbose, &
                          qm, reshape(RH(c,:),(/km,nobs/)),&
                          nMom,nPol, tau, ssa, g, pmom, ierr )
    ! end if

    TAU_(c,:,:) = tau(:,:,nobs)
    SSA_(c,:,:) = ssa(:,:,nobs)
    G_(c,:,:)   = g(:,:,nobs)

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
          call VLIDORT_Scalar_Lambert (km, nch, nobs ,dble(channels), nMom,      &
                  nPol, dble(tau), dble(ssa), dble(g), dble(pmom), dble(pe), dble(ze), dble(te), albedo,&
                  (/dble(SZA(c))/), &
                  (/dble(abs(RAA(c)))/), &
                  (/dble(VZA(c))/), &
                  dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT, ierr)
        else 
          call LIDORT_Scalar_Lambert (km, nch, nobs ,dble(channels), nMom,      &
                  nPol, dble(tau), dble(ssa), dble(g), dble(pmom), dble(pe), dble(ze), dble(te), albedo,&
                  (/dble(SZA(c))/), &
                  (/dble(abs(RAA(c)))/), &
                  (/dble(VZA(c))/), &
                  dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT, ierr)
        end if 
      else
        ! Call to vlidort vector code
        call VLIDORT_Vector_Lambert (km, nch, nobs ,dble(channels), nMom,   &
               nPol, dble(tau), dble(ssa), dble(pmom), dble(pe), dble(ze), dble(te), albedo,&
               (/dble(SZA(c))/), &
               (/dble(abs(RAA(c)))/), &
               (/dble(VZA(c))/), &
               dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT, Q, U, ierr, DO_2OS_CORRECTION)
      end if

    else if ( ANY(kernel_wt == surf_missing) ) then
!     Save code for pixels that were not gap filled
!     ---------------------------------------------    
      radiance_VL_int(nobs,:) = -500
      reflectance_VL_int(nobs,:) = -500
      albedo = -500
      ROT = -500
      if (.not. scalar) then
        Q = -500
        U = -500
      end if
      ierr = 0
    else   
!     MODIS BRDF Surface Model
!     ------------------------------
      if (scalar) then 
        if (vlidort) then
          ! Call to vlidort scalar code            
          call VLIDORT_Scalar_LandMODIS (km, nch, nobs, dble(channels), nMom,  &
                  nPol, dble(tau), dble(ssa), dble(g), dble(pmom), dble(pe), dble(ze), dble(te), &
                  kernel_wt, param, &
                  (/dble(SZA(c))/), &
                  (/dble(abs(RAA(c)))/), &
                  (/dble(VZA(c))/), &
                  dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT, albedo, ierr )  
        else
        ! Call to vlidort scalar code            
          call LIDORT_Scalar_LandMODIS (km, nch, nobs, dble(channels), nMom,  &
                  nPol, dble(tau), dble(ssa), dble(g), dble(pmom), dble(pe), dble(ze), dble(te), &
                  kernel_wt, param, &
                  (/dble(SZA(c))/), &
                  (/dble(abs(RAA(c)))/), &
                  (/dble(VZA(c))/), &
                  dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT, albedo, ierr )  
        end if
      else
   
        ! Call to vlidort vector code
        call VLIDORT_Vector_LandMODIS (km, nch, nobs, dble(channels), nMom, &
                nPol, dble(tau), dble(ssa), dble(pmom), dble(pe), dble(ze), dble(te), &
                kernel_wt, param, &
                (/dble(SZA(c))/), &
                (/dble(abs(RAA(c)))/), &
                (/dble(VZA(c))/), &
                dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT, albedo, Q, U, ierr, DO_2OS_CORRECTION )  
      end if      
    end if          
    
!   Check VLIDORT Status, Store Outputs in Shared Arrays
!   ----------------------------------------------------    
    call mp_check_vlidort(radiance_VL_int,reflectance_VL_int)  
    radiance_VL(c,:)    = radiance_VL_int(nobs,:)
    reflectance_VL(c,:) = reflectance_VL_int(nobs,:)
    ALBEDO_(c,:) = albedo(nobs,:)

    do ch=1,nch      
        ROT_(c,:,ch) = ROT(:,nobs,ch)
    end do

    if (.not. scalar) then
      Q_(c,:)      = Q(nobs,:)
      U_(c,:)      = U(nobs,:)
    end if
    
    write(msg,*) 'VLIDORT Calculations DONE', myid, ierr
    call write_verbose(msg)

!   Keep track of progress of each processor
!   -----------------------------------------        
    if (nint(100.*real(cc-starti)/real(counti)) > progress) then
      progress = nint(100.*real(cc-starti)/real(counti))
      write(*,'(A,I,A,I,A,I2,A,I3,A)') 'Pixel: ',c,'  End Pixel: ',endi,'  ID:',myid,'  Progress:', nint(progress),'%'           
    end if
                
  end do ! do clear pixels

! Wait for everyone to finish calculations
! ----------------------------------------
  call MAPL_SyncSharedMemory(rc=ierr)


! Write to OUT_file and ATMOS_file
! Expand radiance to im x jm using mask
! Write output to correct position in file 
! -----------------------------------------
  if (MAPL_am_I_root()) then
   
    allocate (field(im,jm))
    field = g5nr_missing
    allocate (AOD(clrm_total))

!                             Write to main OUT_File
!                             ----------------------
    call check( nf90_open(OUT_file, nf90_write, ncid), "opening file " // OUT_file )
    do ch = 1, nch
      write(msg,'(F10.2)') channels(ch)
      
      call check(nf90_inq_varid(ncid, 'ref_' // trim(adjustl(msg)), varid), "get ref vaird")
      call check(nf90_put_var(ncid, varid, unpack(reshape(reflectance_VL(:,ch),(/clrm_total/)),clmask,field), &
                  start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out reflectance")

      call check(nf90_inq_varid(ncid, 'surf_ref_' // trim(adjustl(msg)), varid), "get ref vaird")
      call check(nf90_put_var(ncid, varid, unpack(reshape(ALBEDO_(:,ch),(/clrm_total/)),clmask,field), &
                    start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out albedo")

      AOD = 0
      do k=1,km 
        AOD = AOD + reshape(TAU_(:,k,ch),(/clrm_total/))
      end do

      where(reshape(ALBEDO_(:,ch),(/clrm_total/)) .eq. MISSING)  AOD = MISSING
      

      call check(nf90_inq_varid(ncid, 'aod_' // trim(adjustl(msg)), varid), "get aod vaird")
      call check(nf90_put_var(ncid, varid, unpack(AOD,clmask,field), &
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
        call check(nf90_put_var(ncid, varid, unpack(reshape(radiance_VL(:,ch),(/clrm_total/)),clmask,field), &
                      start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out radiance")

        if (.not. scalar) then
          call check(nf90_inq_varid(ncid, 'q_' // trim(adjustl(msg)), varid), "get q vaird")
          call check(nf90_put_var(ncid, varid, unpack(reshape(Q_(:,ch),(/clrm_total/)),clmask,field), &
                      start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out Q")

          call check(nf90_inq_varid(ncid, 'u_' // trim(adjustl(msg)), varid), "get u vaird")
          call check(nf90_put_var(ncid, varid, unpack(reshape(U_(:,ch),(/clrm_total/)),clmask,field), &
                      start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out U")
        endif        

        do k=1,km 

          call check(nf90_inq_varid(ncid, 'aot_' // trim(adjustl(msg)), varid), "get aot vaird")
          call check(nf90_put_var(ncid, varid, unpack(reshape(TAU_(:,k,ch),(/clrm_total/)),clmask,field), &
                    start = (/1,1,k,nobs/), count = (/im,jm,1,nobs/)), "writing out tau")

          call check(nf90_inq_varid(ncid, 'g_' // trim(adjustl(msg)), varid), "get g vaird")
          call check(nf90_put_var(ncid, varid, unpack(reshape(G_(:,k,ch),(/clrm_total/)),clmask,field), &
                    start = (/1,1,k,nobs/), count = (/im,jm,1,nobs/)), "writing out g")

          call check(nf90_inq_varid(ncid, 'ssa_' // trim(adjustl(msg)), varid), "get ssa vaird")
          call check(nf90_put_var(ncid, varid, unpack(reshape(SSA_(:,k,ch),(/clrm_total/)),clmask,field), &
                    start = (/1,1,k,nobs/), count = (/im,jm,1,nobs/)), "writing out ssa")

          call check(nf90_inq_varid(ncid, 'rot_' // trim(adjustl(msg)), varid), "get rot vaird")
          call check(nf90_put_var(ncid, varid, unpack(reshape(ROT_(:,k,ch),(/clrm_total/)),clmask,field), &
                    start = (/1,1,k,nobs/), count = (/im,jm,1,nobs/)), "writing out rot")
        end do
      end do
      call check( nf90_close(ncid), "close atmosfile" )
    end if  !additional_output
    deallocate(field)
  end if

! Sync processors before shutting down
! ------------------------------------
  call MAPL_SyncSharedMemory(rc=ierr)
  if (MAPL_am_I_root()) then
    write(*,*) '<> Finished VLIDORT Simulation of '//trim(lower_to_upper(instname))//' domain'
    call sys_tracker()
    write(*,*) ' '
    call system_clock ( t2, clock_rate, clock_max )
    write ( *, * ) 'Elapsed real time = ', real ( t2 - t1 ) / real ( clock_rate )
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
!     read_cloud
! PURPOSE
!     allocates cloud variables and reads in data
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     28 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
subroutine read_cloud()
  allocate (CLDTOT(im,jm))

  call readvar2D("CLDTOT", MET_file, CLDTOT)
end subroutine read_cloud

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
  call readvar2D("FRLAND", LAND_file, FRLAND)
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
  if (trim(layout) == '111') then
    write(MET_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/', &
                            trim(instname),'-g5nr.lb2.met_Nv.',date,'_',time,'z.nc4'
    write(AER_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/', &
                            trim(instname),'-g5nr.lb2.aer_Nv.',date,'_',time,'z.nc4'
    write(ANG_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/', &
                            trim(instname),'.lb2.angles.',date,'_',time,'z.nc4'

    write(LAND_file,'(4A)') trim(indir),'/LevelB/invariant/',trim(instname),'-g5nr.lb2.asm_Nx.nc4'                                 
  else
    write(MET_file,'(16A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/', &
                            trim(instname),'-g5nr.lb2.met_Nv.',date,'_',time,'z_',trim(layout),'.nc4'
    write(AER_file,'(16A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/', &
                            trim(instname),'-g5nr.lb2.aer_Nv.',date,'_',time,'z_',trim(layout),'.nc4'
    write(ANG_file,'(16A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/', &
                            trim(instname),'.lb2.angles.',date,'_',time,'z_',trim(layout),'.nc4'
    write(LAND_file,'(6A)') trim(indir),'/LevelB/invariant/',trim(instname),'-g5nr.lb2.asm_Nx_',trim(layout),'.nc4'
  end if  

  if ( lower_to_upper(surfmodel) == 'RTLS' ) then
    if ( lower_to_upper(surfname) == 'MAIACRTLS' ) then
      write(SURF_file,'(8A)') trim(indir),'/BRDF/v',trim(surf_version),'/',trim(surfname),'.',surfdate,'.hdf'
    else
      write(SURF_file,'(8A)') trim(indir),'/BRDF/v',trim(surf_version),'/',trim(surfname),'.',surfdate,'.nc4'
    end if
  else
    write(SURF_file,'(6A)') trim(indir),'/SurfLER/',trim(instname),'-omi.SurfLER.',date(5:6),'.nc4'
  end if

  write(INV_file,'(4A)')  trim(indir),'/LevelG/invariant/',trim(instname),'.lg1.invariant.nc4'

! OUTFILES
  if (vlidort) then
    write(OUT_file,'(4A)') trim(outdir),'/',trim(instname),'-g5nr.lc2.vlidort.'
  else
    write(OUT_file,'(4A)') trim(outdir),'/',trim(instname),'-g5nr.lc2.lidort.'
  end if 
  call outfile_extname(OUT_file)

  if (additional_output) then
    if (vlidort) then
      write(ATMOS_file,'(4A)') trim(outdir),'/',trim(instname),'-g5nr.lc.vlidort.'
    else
      write(ATMOS_file,'(4A)') trim(outdir),'/',trim(instname),'-g5nr.lc.lidort.'
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
    else if (nodenum > 10 .and. nodenum < 100) then
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

    call mp_readvar3Dchunk("AIRDENS", AER_file, (/im,jm,km/), 1, npet, myid, AIRDENS_)  
    call mp_readvar3Dchunk("RH", AER_file, (/im,jm,km/), 1, npet, myid, RH_) 
    call mp_readvar3Dchunk("DELP", AER_file, (/im,jm,km/), 1, npet, myid, DELP_) 
    call mp_readvar3Dchunk("DU001", AER_file, (/im,jm,km/), 1, npet, myid, DU001_) 
    call mp_readvar3Dchunk("DU002", AER_file, (/im,jm,km/), 1, npet, myid, DU002_) 
    call mp_readvar3Dchunk("DU003", AER_file, (/im,jm,km/), 1, npet, myid, DU003_) 
    call mp_readvar3Dchunk("DU004", AER_file, (/im,jm,km/), 1, npet, myid, DU004_) 
    call mp_readvar3Dchunk("DU005", AER_file, (/im,jm,km/), 1, npet, myid, DU005_) 
    call mp_readvar3Dchunk("SS001", AER_file, (/im,jm,km/), 1, npet, myid, SS001_) 
    call mp_readvar3Dchunk("SS002", AER_file, (/im,jm,km/), 1, npet, myid, SS002_) 
    call mp_readvar3Dchunk("SS003", AER_file, (/im,jm,km/), 1, npet, myid, SS003_) 
    call mp_readvar3Dchunk("SS004", AER_file, (/im,jm,km/), 1, npet, myid, SS004_) 
    call mp_readvar3Dchunk("SS005", AER_file, (/im,jm,km/), 1, npet, myid, SS005_) 
    call mp_readvar3Dchunk("BCPHOBIC", AER_file, (/im,jm,km/), 1, npet, myid, BCPHOBIC_) 
    call mp_readvar3Dchunk("BCPHILIC", AER_file, (/im,jm,km/), 1, npet, myid, BCPHILIC_) 
    call mp_readvar3Dchunk("OCPHOBIC", AER_file, (/im,jm,km/), 1, npet, myid, OCPHOBIC_) 
    call mp_readvar3Dchunk("OCPHILIC", AER_file, (/im,jm,km/), 1, npet, myid, OCPHILIC_) 
    call mp_readvar3Dchunk("SO4", AER_file, (/im,jm,km/), 1, npet, myid, SO4_) 

    call MAPL_SyncSharedMemory(rc=ierr)    
    if (MAPL_am_I_root()) then
      call reduceProfile(AIRDENS_,clmask,AIRDENS) 
      call reduceProfile(RH_,clmask,RH) 
      call reduceProfile(DELP_,clmask,DELP) 
      call reduceProfile(DU001_,clmask,DU001) 
      call reduceProfile(DU002_,clmask,DU002)
      call reduceProfile(DU003_,clmask,DU003) 
      call reduceProfile(DU004_,clmask,DU004) 
      call reduceProfile(DU005_,clmask,DU005)   
      call reduceProfile(SS001_,clmask,SS001)  
      call reduceProfile(SS002_,clmask,SS002) 
      call reduceProfile(SS003_,clmask,SS003) 
      call reduceProfile(SS004_,clmask,SS004) 
      call reduceProfile(SS005_,clmask,SS005) 
      call reduceProfile(BCPHOBIC_,clmask,BCPHOBIC) 
      call reduceProfile(BCPHILIC_,clmask,BCPHILIC) 
      call reduceProfile(OCPHOBIC_,clmask,OCPHOBIC) 
      call reduceProfile(OCPHILIC_,clmask,OCPHILIC)  
      call reduceProfile(SO4_,clmask,SO4)  
      write(*,*) '<> Read aeorosl data to shared memory'
    end if      

    ! call MAPL_DeallocNodeArray(AIRDENS_,rc=ierr)
    ! call MAPL_DeallocNodeArray(RH_,rc=ierr)
    ! call MAPL_DeallocNodeArray(DELP_,rc=ierr)
    ! call MAPL_DeallocNodeArray(DU001_,rc=ierr)
    ! call MAPL_DeallocNodeArray(DU002_,rc=ierr)
    ! call MAPL_DeallocNodeArray(DU003_,rc=ierr)
    ! call MAPL_DeallocNodeArray(DU004_,rc=ierr)                           
    ! call MAPL_DeallocNodeArray(DU005_,rc=ierr)
    ! call MAPL_DeallocNodeArray(SS001_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(SS002_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(SS003_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(SS004_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(SS005_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(BCPHOBIC_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(BCPHILIC_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(OCPHOBIC_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(OCPHILIC_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(SO4_,rc=ierr)       
    ! write(*,*) 'finished read_aer_Nv'  
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

      if (trim(layout) == '111') then
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
      else
        !read only the chunk that belongs to the tile
        read(layout(1:1),*) Nx
        read(layout(2:2),*) Ny
        read(layout(3:) ,*) ntile

        ! figure out x,y position of tile
        x = mod(ntile, Nx)
        y = int(ntile/Nx)

        xstart = x*(im) + 1
        xend   = xstart + im - 1
        ystart = y*(jm) + 1
        yend   = ystart + jm - 1

        call mp_readTilevar3Dchunk("Band1", SURF_file, (/im,jm,nkernel/), &
                                  (/xstart,ystart,1/), 1, npet, myid, Band1) 
        call mp_readTilevar3Dchunk("Band2", SURF_file, (/im,jm,nkernel/), &
                                  (/xstart,ystart,1/), 1, npet, myid, Band2) 
        call mp_readTilevar3Dchunk("Band3", SURF_file, (/im,jm,nkernel/), &
                                  (/xstart,ystart,1/), 1, npet, myid, Band3) 
        call mp_readTilevar3Dchunk("Band4", SURF_file, (/im,jm,nkernel/), &
                                  (/xstart,ystart,1/), 1, npet, myid, Band4) 
        call mp_readTilevar3Dchunk("Band5", SURF_file, (/im,jm,nkernel/), &
                                  (/xstart,ystart,1/), 1, npet, myid, Band5) 
        call mp_readTilevar3Dchunk("Band6", SURF_file, (/im,jm,nkernel/), &
                                  (/xstart,ystart,1/), 1, npet, myid, Band6) 
        call mp_readTilevar3Dchunk("Band7", SURF_file, (/im,jm,nkernel/), &
                                  (/xstart,ystart,1/), 1, npet, myid, Band7) 
        if (lower_to_upper(surfname) == 'MAIACRTLS') then
          call mp_readTilevar3Dchunk("Band8", SURF_file, (/im,jm,nkernel/), &
                                  (/xstart,ystart,1/), 1, npet, myid, Band8)
        end if
      end if

      call MAPL_SyncSharedMemory(rc=ierr)    

      if (MAPL_am_I_root()) then
        KISO_(:,:,1) = Band1(:,:,1)
        KISO_(:,:,2) = Band2(:,:,1)
        KISO_(:,:,3) = Band3(:,:,1)
        KISO_(:,:,4) = Band4(:,:,1)
        KISO_(:,:,5) = Band5(:,:,1)
        KISO_(:,:,6) = Band6(:,:,1)
        KISO_(:,:,7) = Band7(:,:,1)
        if (lower_to_upper(surfname) == 'MAIACRTLS') KISO_(:,:,8) = Band8(:,:,1)

        KVOL_(:,:,1) = Band1(:,:,2)
        KVOL_(:,:,2) = Band2(:,:,2)
        KVOL_(:,:,3) = Band3(:,:,2)
        KVOL_(:,:,4) = Band4(:,:,2)
        KVOL_(:,:,5) = Band5(:,:,2)
        KVOL_(:,:,6) = Band6(:,:,2)
        KVOL_(:,:,7) = Band7(:,:,2)
        if (lower_to_upper(surfname) == 'MAIACRTLS') KVOL_(:,:,8) = Band8(:,:,2)

        KGEO_(:,:,1) = Band1(:,:,3)
        KGEO_(:,:,2) = Band2(:,:,3)
        KGEO_(:,:,3) = Band3(:,:,3)
        KGEO_(:,:,4) = Band4(:,:,3)
        KGEO_(:,:,5) = Band5(:,:,3)
        KGEO_(:,:,6) = Band6(:,:,3)
        KGEO_(:,:,7) = Band7(:,:,3)
        if (lower_to_upper(surfname) == 'MAIACRTLS') KGEO_(:,:,8) = Band8(:,:,3)

        call reduceProfile(KISO_,clmask,KISO) 
        call reduceProfile(KVOL_,clmask,KVOL)
        call reduceProfile(KGEO_,clmask,KGEO)  
        write(*,*) '<> Read BRDF data to shared memory' 
      end if 

    else
      if (MAPL_am_I_root()) then
        call read_LER(temp)        
        call reduceProfile(temp,clmask,LER)
        write(*,*) '<> Read LER data to shared memory'
      end if
    end if

  end subroutine read_surf

  subroutine read_LER(indata)
    real, intent(inout),dimension(im,jm,surfbandm)         :: indata
    real, dimension(im,jm,1,1)                             :: temp

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
    real, allocatable          :: saa(:), vaa(:)
    integer                    :: i


    call mp_readvar2Dchunk("solar_zenith",   ANG_file, (/im,jm/), 1, npet, myid, SZA_) 
    call mp_readvar2Dchunk("sensor_zenith",  ANG_file, (/im,jm/), 1, npet, myid, VZA_) 
    call mp_readvar2Dchunk("solar_azimuth",  ANG_file, (/im,jm/), 1, npet, myid, SAA_)
    call mp_readvar2Dchunk("sensor_azimuth", ANG_file, (/im,jm/), 1, npet, myid, VAA_)
    call MAPL_SyncSharedMemory(rc=ierr)    
    if (MAPL_am_I_root()) then
      SZA = pack(SZA_,clmask)
      VZA = pack(VZA_,clmask)
        
      allocate (saa(clrm))
      allocate (vaa(clrm))
      
      saa = pack(SAA_,clmask)
      ! define according to photon travel direction
      saa = saa + 180.0
      do i = 1, clrm
        if (saa(i) >= 360.0) then
          saa(i) = saa(i) - 360.0
        end if
      end do

      vaa = pack(VAA_,clmask)

      RAA = vaa - saa

      deallocate (saa)
      deallocate (vaa)
      write(*,*) '<> Read angle data to shared memory' 

    end if
    ! call MAPL_DeallocNodeArray(SZA_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(VZA_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(SAA_,rc=ierr) 
    ! call MAPL_DeallocNodeArray(VAA_,rc=ierr) 

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
  subroutine allocate_shared(clrm)
    integer, intent(in)    :: clrm

    call MAPL_AllocNodeArray(AIRDENS,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(RH,(/clrm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DELP,(/clrm,km/),rc=ierr)
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
    if (lower_to_upper(surfmodel) == 'RTLS') then
      call MAPL_AllocNodeArray(KISO,(/clrm,surfbandm/),rc=ierr)
      call MAPL_AllocNodeArray(KVOL,(/clrm,surfbandm/),rc=ierr)
      call MAPL_AllocNodeArray(KGEO,(/clrm,surfbandm/),rc=ierr)
    else
      call MAPL_AllocNodeArray(LER,(/clrm,surfbandm/),rc=ierr)
    end if 

    call MAPL_AllocNodeArray(SZA,(/clrm/),rc=ierr)
    call MAPL_AllocNodeArray(VZA,(/clrm/),rc=ierr)
    call MAPL_AllocNodeArray(RAA,(/clrm/),rc=ierr)

    call MAPL_AllocNodeArray(TAU_,(/clrm,km,nch/),rc=ierr)
    call MAPL_AllocNodeArray(SSA_,(/clrm,km,nch/),rc=ierr)
    call MAPL_AllocNodeArray(G_,(/clrm,km,nch/),rc=ierr)
    call MAPL_AllocNodeArray(ROT_,(/clrm,km,nch/),rc=ierr)
    call MAPL_AllocNodeArray(ALBEDO_,(/clrm,nch/),rc=ierr)
    
    if (.not. scalar) then
      call MAPL_AllocNodeArray(Q_,(/clrm,nch/),rc=ierr)
      call MAPL_AllocNodeArray(U_,(/clrm,nch/),rc=ierr)
    end if

    call MAPL_AllocNodeArray(radiance_VL,(/clrm,nch/),rc=ierr)
    call MAPL_AllocNodeArray(reflectance_VL,(/clrm,nch/),rc=ierr)

    call MAPL_AllocNodeArray(AIRDENS_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(RH_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DELP_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU001_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU002_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU003_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU004_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(DU005_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS001_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS002_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS003_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS004_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SS005_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(BCPHOBIC_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(BCPHILIC_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(OCPHOBIC_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(OCPHILIC_,(/im,jm,km/),rc=ierr)
    call MAPL_AllocNodeArray(SO4_,(/im,jm,km/),rc=ierr)
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
      call MAPL_AllocNodeArray(KISO_,(/im,jm,surfbandm/),rc=ierr)
      call MAPL_AllocNodeArray(KVOL_,(/im,jm,surfbandm/),rc=ierr)
      call MAPL_AllocNodeArray(KGEO_,(/im,jm,surfbandm/),rc=ierr)
    end if

    call MAPL_AllocNodeArray(SZA_,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(VZA_,(/im,jm/),rc=ierr)    
    call MAPL_AllocNodeArray(SAA_,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(VAA_,(/im,jm/),rc=ierr) 
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
    allocate (pe(km+1,nobs))
    allocate (ze(km+1,nobs))
    allocate (te(km+1,nobs))
    allocate (qm(km,nq,nobs))
    allocate (tau(km,nch,nobs))
    allocate (ssa(km,nch,nobs))
    allocate (g(km,nch,nobs))
    allocate (albedo(nobs,nch))

    allocate (radiance_VL_int(nobs,nch))
    allocate (reflectance_VL_int(nobs, nch))    

    if (lower_to_upper(surfmodel) == 'RTLS') then
      allocate (kernel_wt(nkernel,nch,nobs))
    end if
    allocate (param(nparam,nch,nobs))

    allocate (ROT(km,nobs,nch))
    allocate (pmom(km,nch,nobs,nMom,nPol))

    if (.not. scalar) then      
      allocate (Q(nobs, nch))
      allocate (U(nobs, nch))
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
    
    integer                            :: ncid
    integer                            :: timeDimID, ewDimID, nsDimID, levDimID, chaDimID       
    integer                            :: scantimeVarID, clonVarID, clatVarID
    integer                            :: timeVarID, levVarID, ewVarID, nsVarID
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

    write(comment,'(A)') 'VLIDORT simulation run from geo_vlidort.x'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'history',trim(comment)),"history attr")

    call check(nf90_put_att(ncid,NF90_GLOBAL,'grid_inputs',trim(INV_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'angle_inputs',trim(ANG_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'land_inputs',trim(LAND_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'met_inputs',trim(MET_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'aerosol_inputs',trim(AER_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'surface_inputs',trim(SURF_file)),"input files attr")

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

    call readvar1D("scanTime", MET_file, scantime)
    call check(nf90_put_var(ncid,scantimeVarID,scantime), "writing out scantime")

    call readvar2D("clon", MET_file, clon)
    call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

    call readvar2D("clat", MET_file, clat)
    call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

    call readvar1D("time", MET_file, tyme)
    call check(nf90_put_var(ncid,timeVarID,tyme), "writing out time")

    call readvar1D("lev", MET_file, lev)
    call check(nf90_put_var(ncid,levVarID,lev), "writing out lev")

    call readvar1D("ew", MET_file, ew)
    call check(nf90_put_var(ncid,ewVarID,ew), "writing out ew")

    call readvar1D("ns", MET_file, ns)
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
      call check(nf90_def_dim(ncid, "ew", im, ewDimID), "creating ew dimension") !im
      call check(nf90_def_dim(ncid, "ns", jm, nsDimID), "creating ns dimension") !jm

      ! Global Attributes
      write(comment,'(A)') 'Atmospheric inputs for VLIDORT Simulation of GEOS-5 '//lower_to_upper(trim(instname))//' Sampler'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'title',trim(comment)),"title attr")

      write(comment,'(A)') 'NASA/Goddard Space Flight Center'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'institution',trim(comment)),"institution attr")

      write(comment,'(A)') 'Global Model and Assimilation Office'
      call check(nf90_put_att(ncid,NF90_GLOBAL,'source',trim(comment)),"source attr")

      write(comment,'(A)') 'VLIDORT simulation run from geo_vlidort.x'
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

      call readvar1D("scanTime", MET_file, scantime)
      call check(nf90_put_var(ncid,scantimeVarID,scantime), "writing out scantime")

      call readvar2D("clon", MET_file, clon)
      call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

      call readvar2D("clat", MET_file, clat)
      call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

      call readvar1D("time", MET_file, tyme)
      call check(nf90_put_var(ncid,timeVarID,tyme), "writing out time")

      call readvar1D("lev", MET_file, lev)
      call check(nf90_put_var(ncid,levVarID,lev), "writing out lev")

      call readvar1D("ew", MET_file, ew)
      call check(nf90_put_var(ncid,ewVarID,ew), "writing out ew")

      call readvar1D("ns", MET_file, ns)
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
      qm(k,1,nobs) = DU001(c,k)*DELP(c,k)/grav
      qm(k,2,nobs) = DU002(c,k)*DELP(c,k)/grav
      qm(k,3,nobs) = DU003(c,k)*DELP(c,k)/grav
      qm(k,4,nobs) = DU004(c,k)*DELP(c,k)/grav
      qm(k,5,nobs) = DU005(c,k)*DELP(c,k)/grav
      qm(k,6,nobs) = SS001(c,k)*DELP(c,k)/grav
      qm(k,7,nobs) = SS002(c,k)*DELP(c,k)/grav
      qm(k,8,nobs) = SS003(c,k)*DELP(c,k)/grav
      qm(k,9,nobs) = SS004(c,k)*DELP(c,k)/grav
      qm(k,10,nobs) = SS005(c,k)*DELP(c,k)/grav
      qm(k,11,nobs) = BCPHOBIC(c,k)*DELP(c,k)/grav
      qm(k,12,nobs) = BCPHILIC(c,k)*DELP(c,k)/grav
      qm(k,13,nobs) = OCPHOBIC(c,k)*DELP(c,k)/grav
      qm(k,14,nobs) = OCPHILIC(c,k)*DELP(c,k)/grav
      qm(k,15,nobs) = SO4(c,k)*DELP(c,k)/grav
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

    call shmem_test2D('RH',RH)
    call shmem_test2D('AIRDENS',AIRDENS)
    call shmem_test2D('DELP',DELP)
    call shmem_test2D('DU001',DU001)
    call shmem_test2D('DU002',DU002)
    call shmem_test2D('DU003',DU003)
    call shmem_test2D('DU004',DU004)
    call shmem_test2D('DU005',DU005)
    call shmem_test2D('SS001',SS001)
    call shmem_test2D('SS002',SS002)
    call shmem_test2D('SS003',SS003)
    call shmem_test2D('SS004',SS004)
    call shmem_test2D('SS005',SS005)
    call shmem_test2D('BCPHOBIC',BCPHOBIC)
    call shmem_test2D('BCPHILIC',BCPHILIC)
    call shmem_test2D('OCPHOBIC',OCPHOBIC)
    call shmem_test2D('OCPHILIC',OCPHILIC)
    call shmem_test2D('SO4',OCPHILIC)

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
    call ESMF_ConfigGetAttribute(cf, indir, label = 'INDIR:',__RC__)
    call ESMF_ConfigGetAttribute(cf, outdir, label = 'OUTDIR:',default=indir)
    call ESMF_ConfigGetAttribute(cf, surfname, label = 'SURFNAME:',default='MAIACRTLS')
    call ESMF_ConfigGetAttribute(cf, surfmodel, label = 'SURFMODEL:',default='RTLS')
    call ESMF_ConfigGetAttribute(cf, surfdate, label = 'SURFDATE:',__RC__)
    call ESMF_ConfigGetAttribute(cf, scalar, label = 'SCALAR:',default=.TRUE.)
    call ESMF_ConfigGetAttribute(cf, szamax, label = 'SZAMAX:',default=80.0)
    call ESMF_ConfigGetAttribute(cf, vzamax, label = 'VZAMAX:',default=80.0)
    call ESMF_ConfigGetAttribute(cf, cldmax, label = 'CLDMAX:',default=0.01)
    call ESMF_ConfigGetAttribute(cf, surfband, label = 'SURFBAND:', default='INTERPOLATE')
    call ESMF_ConfigGetAttribute(cf, surfbandm, label = 'SURFBANDM:',__RC__)
    call ESMF_ConfigGetAttribute(cf, additional_output, label = 'ADDITIONAL_OUTPUT:',default=.false.)
    call ESMF_ConfigGetAttribute(cf, nodemax, label = 'NODEMAX:',default=1) 
    call ESMF_ConfigGetAttribute(cf, version, label = 'VERSION:',default='1.0') 
    call ESMF_ConfigGetAttribute(cf, surf_version, label = 'SURF_VERSION:',default='1.0')    
    call ESMF_ConfigGetAttribute(cf, layout, label = 'LAYOUT:',default='111')    
    call ESMF_ConfigGetAttribute(cf, vlidort, label = 'VLIDORT:',default=.true.)  
    call ESMF_ConfigGetAttribute(cf, DO_2OS_CORRECTION, label = 'DO_2OS_CORRECTION:',default=.false.)        


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

    call MAPL_DeallocNodeArray(TAU_,rc=ierr)
    call MAPL_DeallocNodeArray(SSA_,rc=ierr)
    call MAPL_DeallocNodeArray(G_,rc=ierr)
    call MAPL_DeallocNodeArray(ROT_,rc=ierr)
    call MAPL_DeallocNodeArray(ALBEDO_,rc=ierr)
    if (.not. scalar) then
      call MAPL_DeallocNodeArray(Q_,rc=ierr)
      call MAPL_DeallocNodeArray(U_,rc=ierr)
    end if

    call MAPL_DeallocNodeArray(radiance_VL,rc=ierr) 
    call MAPL_DeallocNodeArray(reflectance_VL,rc=ierr) 

  end subroutine deallocate_shared


end program geo_vlidort2OS_NNTestData
