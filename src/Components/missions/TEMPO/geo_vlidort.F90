!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    geo_vlidort
! PURPOSE
!     Reads in parallel the model data (in a netcdf file) interpolated to TEMPO grid 
!     The variables are needed as input to vlidort
!     A shared memory array created with a call to MAPL_ShmemMod is used for the variables
! INPUT
!     varname  : string of variable name
!     filename : file to be read
!     var      : the variable to be read to
! OUTPUT
!     None
!  HISTORY
!     27 April P. Castellanos adapted from A. da Silva shmem_reader.F90
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

program shmem_reader

  use MAPL_ShmemMod    ! The SHMEM infrastructure
  use netcdf           ! for reading the NR files

  implicit none
  include "mpif.h"

! Date 
! --------
  character(len=8)      :: date 
  character(len=2)      :: time 

! Test flag
! -----------
  logical               :: test_shmem = .False.

! File names
! ----------
  character(len=256)    :: MET_file, AER_file, INV_file, BRDF_file, OUT_file 

! Global, 3D arrays to be allocated using SHMEM
! ---------------------------------------------
  real, pointer         :: CLON(:,:) => null()
  real, pointer         :: CLAT(:,:) => null()
  real, pointer         :: CLDTOT(:,:) => null()
  real, pointer         :: AIRDENS(:,:,:) => null()
  real, pointer         :: RH(:,:,:) => null()
  real, pointer         :: DELP(:,:,:) => null()
  real, pointer         :: DU001(:,:,:) => null()
  real, pointer         :: DU002(:,:,:) => null()
  real, pointer         :: DU003(:,:,:) => null()
  real, pointer         :: DU004(:,:,:) => null()
  real, pointer         :: DU005(:,:,:) => null()
  real, pointer         :: SS001(:,:,:) => null()
  real, pointer         :: SS002(:,:,:) => null()
  real, pointer         :: SS003(:,:,:) => null()
  real, pointer         :: SS004(:,:,:) => null()
  real, pointer         :: SS005(:,:,:) => null()
  real, pointer         :: BCPHOBIC(:,:,:) => null()
  real, pointer         :: BCPHILIC(:,:,:) => null()
  real, pointer         :: OCPHOBIC(:,:,:) => null()
  real, pointer         :: OCPHILIC(:,:,:) => null()
  real*8, pointer         :: KISO(:,:,:) => null()
  real*8, pointer         :: KVOL(:,:,:) => null()
  real*8, pointer         :: KGEO(:,:,:) => null()


! VLIDORT input arrays
! ---------------------------
  real, pointer                         :: pe(:,:) => null()      ! edge pressure [Pa]
  real, pointer                         :: ze(:,:) => null()      ! edge height above sfc [m]
  real, pointer                         :: te(:,:) => null()      ! edge Temperature [K]
  integer, parameter                    :: nobs = 1                   ! number of profiles VLIDORT will work on
  real, parameter                       :: ptop = 1.0             ! top (edge) pressure [Pa]
  integer, parameter                    :: nch = 1                    ! number of channels
  integer, parameter                    :: nq = 14                ! number of tracers
  integer, parameter                    :: verbose = 0
  real, pointer                         :: qm(:,:,:) => null()    ! (mixing ratio) * delp/g
  real, pointer                         :: tau(:,:,:) => null()   ! aerosol optical depth
  real, pointer                         :: ssa(:,:,:) => null()   ! single scattering albedo
  real, pointer                         :: g(:,:,:) => null()     ! asymmetry factor
  real, pointer                         :: albedo(:,:) => null()  ! surface albedo
  real, pointer                         :: solar_zenith(:) => null() 
  real, pointer                         :: relat_azimuth(:) => null()
  real, pointer                         :: sensor_zenith(:) => null()
  integer, parameter                    :: nkernel = 3
  integer, parameter                    :: nparam  = 2
  real*8, pointer                       :: kernel_wt(:,:,:) => null()  ! kernel weights (/fiso,fgeo,fvol/)
  real*8, pointer                       :: param(:,:,:) => null()      ! Li-Sparse parameters 
                                                                 ! param1 = crown relative height (h/b)
                                                                 ! param2 = shape parameter (b/r)

  real                                  :: sensor_azimuth
  real                                  :: solar_azimuth
  real, dimension(nch),parameter        :: channels = (/405.0/)    ! channels to simulate
  real, parameter                       :: MISSING = -999
  real*8,pointer                        :: radiance_VL(:,:) => null()      ! TOA normalized radiance from VLIDORT
  real*8,pointer                        :: reflectance_VL(:,:) => null()  ! TOA reflectance from VLIDORT
  real*8,pointer                        :: radiance_VL_LandMODIS(:,:) => null()
  real*8,pointer                        :: reflectance_VL_LandMODIS(:,:) => null()

  character(len=*), parameter           :: rcfile = 'Aod_EOS.rc'  ! resource file
  character(len=16), parameter          :: vnames_string(nq) = (/'du001', 'du002', 'du003', 'du004', 'du005', &
                                                             'ss001', 'ss002', 'ss003', 'ss004', 'ss005', &
                                                             'BCphobic', 'BCphilic',                      &
                                                             'OCphobic', 'OCphilic'/) ! array of variable name strings
  character                             :: vnames(nq,16)          ! character array of variable names


! Miscellaneous
! -------------
  integer               :: ierr                      ! MPI error message
  integer               :: status(MPI_STATUS_SIZE)   ! MPI status
  integer               :: myid, npet, CoresPerNode  ! MPI dimensions and processor id
  integer               :: im, jm, km                ! size of TEMPO domain
  integer               :: i, j, k, n                ! size of TEMPO domain
  integer, parameter    :: nbands = 7                ! number of modis bands
  integer, dimension(nch), parameter    :: bands = (/3/)             ! modis band that overlaps with vlidort channel
  integer               :: ch                        ! i-channel
  integer               :: clrm, c                   ! number of clear 
  integer               :: ncid, radVarID, refVarID  ! netcdf ids
  integer               :: p                         ! i-processor
  integer               :: doy                       ! day of year
  integer               :: starti, counti, endi      ! array indices and counts for e-w columns to be read
  integer, allocatable  :: nclr(:)                   ! how many clear pixels each processor reads

  character(len=100)    :: msg                       ! message printed to stdout

  real                  :: sat_lon                   ! satellite longitude
  real                  :: sat_lat                   ! satellite latitude
  real                  :: sat_alt                   ! satellite altitude
  real                  :: Earth_rad                 ! earth radius
  real, parameter       :: grav = 9.81               ! gravity
  real                  :: progress 

  integer               :: t1, t2, clock_max
  real*8                :: clock_rate

  integer,allocatable,dimension(:)    :: i_work, j_work     ! working indeces
  real, parameter                     :: cldmax = 0.01

  call system_clock ( t1, clock_rate, clock_max )

  progress = 0.1

! Initialize MPI
! --------------
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,npet,ierr)
  if (myid == 0) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'

! Initialize SHMEM
! ----------------
  CoresPerNode = MAPL_CoresPerNodeGet(MPI_COMM_WORLD,rc=ierr) ! a must
  call MAPL_InitializeShmem(rc=ierr)

! For now hard code file names and dimensions, could query file for dimensions
! ----------------------------------------------------------------------------
  date = '20051231'
  time = '18'
  im = 1250
  jm = 2000
  km = 72

! INFILES
  write(MET_file,'(A,A,A,A,A,A,A,A,A)') '/nobackup/TEMPO/met_Nv/Y',date(1:4),'/M',date(5:6),'/tempo-g5nr.lb2.met_Nv.',date,'_',time,'z.nc4'
  write(AER_file,'(A,A,A,A,A,A,A,A,A)') '/nobackup/TEMPO/aer_Nv/Y',date(1:4),'/M',date(5:6),'/tempo-g5nr.lb2.aer_Nv.',date,'_',time,'z.nc4'
  INV_file = "/nobackup/TEMPO/LevelG/invariant/tempo.lg1.invariant.nc4"
  BRDF_file = "/nobackup/TEMPO/BRDF/raw/MAIACRTLS.2006008.hdf"

! OUTFILES
  write(OUT_file,'(A,A,A,A,A)') 'tempo-g5nr.lb2.vlidort.',date,'_',time,'z.nc4'
  call mp_create_outfile(ncid,radVarID,refVarID)

  call date2doy()
 
! Allocate the Global arraya using SHMEM
! It will be available on all processors
! ---------------------------------------------------------
  call MAPL_AllocNodeArray(CLON,(/im,jm/),rc=ierr)
  call MAPL_AllocNodeArray(CLAT,(/im,jm/),rc=ierr)
  call MAPL_AllocNodeArray(CLDTOT,(/im,jm/),rc=ierr)
  call MAPL_AllocNodeArray(CLDTOT,(/im,jm/),rc=ierr)
  call MAPL_AllocNodeArray(AIRDENS,(/im,jm,km/),rc=ierr)
  call MAPL_AllocNodeArray(RH,(/im,jm,km/),rc=ierr)
  call MAPL_AllocNodeArray(DELP,(/im,jm,km/),rc=ierr)
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
  call MAPL_AllocNodeArray(KISO,(/im,jm,nbands/),rc=ierr)
  call MAPL_AllocNodeArray(KVOL,(/im,jm,nbands/),rc=ierr)
  call MAPL_AllocNodeArray(KGEO,(/im,jm,nbands/),rc=ierr)

! Allocate arrays that will be copied on each processor - unshared
! ---------------------------------------------------------

! Needed for vlidort
! -----------------------
  allocate (pe(km+1,nobs))
  allocate (ze(km+1,nobs))
  allocate (te(km+1,nobs))
  allocate (qm(km,nq,nobs))
  allocate (tau(km,nch,nobs))
  allocate (ssa(km,nch,nobs))
  allocate (g(km,nch,nobs))
  allocate (albedo(nch,nobs))

  allocate (solar_zenith(nobs))
  allocate (relat_azimuth(nobs))
  allocate (sensor_zenith(nobs))

  allocate (radiance_VL(nobs,nch))
  allocate (reflectance_VL(nobs, nch))

  allocate (radiance_VL_LandMODIS(nobs,nch))
  allocate (reflectance_VL_LandMODIS(nobs, nch))

  allocate (kernel_wt(nkernel,nch,nobs))
  allocate (param(nparam,nch,nobs))

! Needed for reading
! ----------------------
  allocate (nclr(npet)) 
  nclr = 0

! Read the cloud and geographic data
! -------------------------------------
  call mp_colreadvar("CLDTOT", MET_file, CLDTOT)
  call mp_colreadvar("clat", INV_FILE, CLAT)
  call mp_colreadvar("clon", INV_FILE, CLON)
  call mp_readGattr("sat_lat", INV_FILE, sat_lat)
  call mp_readGattr("sat_lon", INV_FILE, sat_lon)
  call mp_readGattr("sat_alt", INV_FILE, sat_alt)
  call mp_readGattr("Earth_radius", INV_FILE, earth_rad) 
  
  if (myid == 0) then
    write(*,*) 'Allocated all shared memory'
    write(*,*) 'Read cloud and geographic information'
    call sys_tracker()
  end if

! Wait for everyone to have access to what's been read into shared memory
! -----------------------------------------------------------   
  call MAPL_SyncSharedMemory(rc=ierr)

! Create cloud mask
! ------------------------
  ! first check that you can!!!!!
  if (.not. ANY(CLDTOT <= cldmax)) then
    write(*,*) 'domain is too cloudy, nothing to do'
    write(*,*) 'Exiting.....'
    call MPI_ABORT(MPI_COMM_WORLD,myid)
  endif

  ! figure out how many indices to work on
  clrm = 0
  do i=1,im
    do j=1,jm
      if (CLDTOT(i,j) <= cldmax) then
        clrm = clrm + 1
      end if 
    end do
  end do

  ! store them
  allocate(i_work(clrm))
  allocate(j_work(clrm))

  clrm = 0
  do i=1,im
    do j=1,jm
      if (CLDTOT(i,j) <= cldmax) then
        clrm = clrm + 1
        i_work(clrm) = i
        j_work(clrm) = j
      end if 
    end do
  end do

  if (myid == 0) then
    write(*,*) 100.*clrm/(im*jm),'% of the domain is clear'
    write(*,*) 'simulating ',clrm,' pixels'
  end if 

! Split up domain among processors
  if (npet >= clrm) then
    nclr(1:npet) = 1
  else if (npet < clrm) then
    nclr(1:npet) = clrm/npet
    nclr(npet)   = nclr(npet) + mod(clrm,npet)
  end if 

! Read the netcdf variables layer-by-layer in parallel
! ------------------------------
  call mp_layreadvar("AIRDENS", AER_file, AIRDENS)
  call mp_layreadvar("RH", AER_file, RH)
  call mp_layreadvar("DELP", AER_file, DELP)
  call mp_layreadvar("DU001", AER_file, DU001)
  call mp_layreadvar("DU002", AER_file, DU002)
  call mp_layreadvar("DU003", AER_file, DU003)
  call mp_layreadvar("DU004", AER_file, DU004)
  call mp_layreadvar("DU005", AER_file, DU005)
  call mp_layreadvar("SS001", AER_file, SS001)
  call mp_layreadvar("SS002", AER_file, SS002)
  call mp_layreadvar("SS003", AER_file, SS003)
  call mp_layreadvar("SS004", AER_file, SS004)
  call mp_layreadvar("SS005", AER_file, SS005)
  call mp_layreadvar("BCPHOBIC", AER_file, BCPHOBIC)
  call mp_layreadvar("BCPHILIC", AER_file, BCPHILIC)
  call mp_layreadvar("OCPHOBIC", AER_file, OCPHOBIC)
  call mp_layreadvar("OCPHILIC", AER_file, OCPHILIC)
  call mp_readvar3D("Kiso", BRDF_file, (/im,jm,nbands/), 1, KISO)
  call mp_readvar3D("Kvol", BRDF_file, (/im,jm,nbands/), 1, KVOL)
  call mp_readvar3D("Kgeo", BRDF_file, (/im,jm,nbands/), 1, KGEO)

! Wait for everyone to finish and print max memory used
! -----------------------------------------------------------  
  call MAPL_SyncSharedMemory(rc=ierr)
  if (myid == 0) then 
    write(*,*) 'Read all variables' 
    call sys_tracker()   
  end if     


  if (test_shmem) then
! Although read on individual PEs, all shared variables should have the same
! data in all PEs. Let's verify that.
! -----------------------------------------------------------     
    call do_testing()
  end if 


! Prepare inputs and run VLIDORT
! ------------------------------
  call strarr_2_chararr(vnames_string,nq,16,vnames)
  
  do p = 7,7 !npet-1
    if (myid == p) then
      if (p == 0) then
        starti = 1
      else
        starti = sum(nclr(1:p))+1
      end if
      counti = nclr(p+1)
      endi   = starti + counti - 1

      do c = 571506, 571506 !starti, endi
        call getEdgeVars ( km, nobs, reshape(AIRDENS(i_work(c),j_work(c),:),(/km,nobs/)), reshape(DELP(i_work(c),j_work(c),:),(/km,nobs/)), ptop, &
                           pe, ze, te )   

        if (verbose==1) then
          write(*,'(A,I)') 'getEdgeVars ', myid
          ! write(msg,'(A,I)') 'getEdgeVars ', myid
          ! call mp_write_msg(msg)                          
        endif

        do n = 1, nobs
          call calc_qm(DU001,1,n,i_work(c),j_work(c))
          call calc_qm(DU002,2,n,i_work(c),j_work(c))
          call calc_qm(DU003,3,n,i_work(c),j_work(c))
          call calc_qm(DU004,4,n,i_work(c),j_work(c))
          call calc_qm(DU005,5,n,i_work(c),j_work(c))
          call calc_qm(SS001,6,n,i_work(c),j_work(c))
          call calc_qm(SS002,7,n,i_work(c),j_work(c))
          call calc_qm(SS003,8,n,i_work(c),j_work(c))
          call calc_qm(SS004,9,n,i_work(c),j_work(c))
          call calc_qm(SS005,10,n,i_work(c),j_work(c))
          call calc_qm(BCPHOBIC,11,n,i_work(c),j_work(c))
          call calc_qm(BCPHILIC,12,n,i_work(c),j_work(c))
          call calc_qm(OCPHOBIC,13,n,i_work(c),j_work(c))
          call calc_qm(OCPHILIC,14,n,i_work(c),j_work(c))

          if (verbose==1) then
            write(*,*) 'calc_qm ', myid
            ! write(msg,'(A,I)') 'calc_qm ', myid
            ! call mp_write_msg(msg) 
          endif

          call sensor_geometry(CLAT(i_work(c),j_work(c)),CLON(i_work(c),j_work(c)),n)


          if (verbose==1) then
            write(*,*) 'sensor_geometry ', myid, sensor_zenith, sensor_azimuth
            ! write(msg,'(A,I)') 'sensor_geometry ', myid
            ! call mp_write_msg(msg) 
          endif

          call solar_geometry(CLAT(i_work(c),j_work(c)),CLON(i_work(c),j_work(c)),n)

          if (verbose==1) then
            write(*,*) 'solar_geometry ', myid, solar_zenith, solar_azimuth
            ! write(msg,*) 'solar_geometry ', myid, solar_zenith, solar_azimuth
            ! call mp_write_msg(msg) 
          end if 

          if (solar_zenith(n) < 90.0) then
            relat_azimuth(n) = abs(sensor_azimuth - solar_azimuth)        
          end if

          do ch = 1, nch
            kernel_wt(:,ch,n) = (/KISO(i,j,bands(ch)),KGEO(i,j,bands(ch)),KVOL(i,j,bands(ch))/)
            param(:,ch,n)     = (/2,1/)
          end do
        end do           

        if (solar_zenith(nobs) < 90.0) then
 
          call getAOPscalar ( km, nobs, nch, nq, rcfile, channels, vnames, verbose, &
                         qm, reshape(RH(i_work(c),j_work(c),:),(/km,nobs/)), &
                         tau, ssa, g, ierr )

          if (verbose==1) then
            write(*,*) 'getAOPscalar ', myid
            ! write(msg,'(A,I)') 'getAOPscalar ', myid
            ! call mp_write_msg(msg) 
          endif

          albedo(:,:)      = 0.05

          ! Call to vlidort scalar code
          ! Simple lambertian surface
          ! ------------------------------
          call Scalar (km, nch, nobs ,dble(channels),        &
                  dble(tau), dble(ssa), dble(g), dble(pe), dble(ze), dble(te), dble(albedo),&
                  dble(solar_zenith), dble(relat_azimuth), dble(sensor_zenith), &
                  dble(MISSING),verbose,radiance_VL,reflectance_VL, ierr)

          call Scalar_LandMODIS (km, nch, nobs, dble(channels),        &
                   dble(tau), dble(ssa), dble(g), dble(pe), dble(ze), dble(te), &
                   kernel_wt, param, &
                   dble(solar_zenith), dble(relat_azimuth), dble(sensor_zenith), &
                   dble(MISSING),verbose,radiance_VL_LandMODIS,reflectance_VL_LandMODIS,ierr )

          if (verbose==1) then
            write(*,*) 'VLIDORT scalar ', myid, ierr
            ! write(msg,'(A,I,I)') 'VLIDORT scalar ', myid, ierr
            ! call mp_write_msg(msg) 
          endif

          if (ANY(radiance_VL == MISSING) .or. ANY(reflectance_VL == MISSING) &
              .or. ANY(radiance_VL_LandMODIS == MISSING) .or. ANY(reflectance_VL_LandMODIS == MISSING)) then
            write(*,*) 'VLIDORT returned a missing value'
            write(*,*) 'Exiting......'
            call MPI_ABORT(MPI_COMM_WORLD,myid)
          else if (ierr > 0) then
            write(*,*) 'VLIDORT returned rc code ', ierr
            write(*,*) 'Exiting......'
            call MPI_ABORT(MPI_COMM_WORLD,myid)
          end if 

          ! Write output to correct position in array
          ! --------------------------------------------
          call mp_write_outfile(ncid,radVarID,refVarID,(/i_work(c),j_work(c),1,1/),(/1,1,nobs,nch/),radiance_VL,reflectance_VL)

          ! Keep track of progress of each processor
          ! -----------------------------------------
          ! progress = 100.*real(c-starti)/real(counti)
          ! write(msg,'(A,I,A,I,F,A)') 'geo_vlidort ',myid,' ',c, progress,'%'
          ! call mp_write_msg(msg)

        end if 
        !if (100.*real(c-starti)/real(counti) >= progress) then
          progress = 100.*real(c-starti)/real(counti)
          write(*,*) 'geo_vlidort ',myid,' ',c,endi, progress,'%'
          !progress = progress + 0.1
        !end if
      
      end do
    end if
  end do 

! Everyone close OUT_file and sync up
! -----------------------------------------

  call check( nf90_close(ncid), "close outfile" )
  call MAPL_SyncSharedMemory(rc=ierr)


  if (myid == 0) then
    write(*,*) 'Called VLIDORT Scalar' 
    call sys_tracker()
    call system_clock ( t2, clock_rate, clock_max )
    write ( *, * ) 'Elapsed real time = ', real ( t2 - t1 ) / real ( clock_rate )
  end if

! All done
! --------
  call MAPL_SyncSharedMemory(rc=ierr)
  call shutdown()

! -----------------------------------------------------------------------------------

  contains

  subroutine mp_write_msg(msg)
    character(len=100), intent(in)    :: msg
    character(len=100)                :: msgr
    integer                           :: pp

    if ( myid  == 0 ) then
      do pp = 1,npet-1
        call mpi_recv(msgr, 100, MPI_CHARACTER, pp, 1, MPI_COMM_WORLD, status, ierr)
        write(*,*) trim(msgr)
      end do
      write(*,*) trim(msg)
    else
      call mpi_send(msg, 100, MPI_CHARACTER, 0, 1, MPI_COMM_WORLD, ierr)
    end if 

  end subroutine mp_write_msg
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     mp_create_outfile
! PURPOSE
!     creates a netcdf4 file that can be written to my multiple processors independently
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     6 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine mp_create_outfile(ncid,radVarID,refVarID)
    integer,intent(out)                :: ncid
    integer,intent(out)                :: radVarID, refVarID    
    
    integer, dimension(4)              :: chunk_size
    integer, dimension(4)              :: dimids
    integer                            :: timeDimID, ewDimID, nsDimID, chaDimID
    integer                            :: chaVarID



    call check(nf90_create(OUT_file, IOR(IOR(nf90_netcdf4, nf90_clobber),nf90_mpiio), ncid, &
                               comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "creating file " // OUT_file)
    call check(nf90_def_dim(ncid, "time", 1, timeDimID), "creating time dimension")
    call check(nf90_def_dim(ncid, "ew", im, ewDimID), "creating ew dimension") !im
    call check(nf90_def_dim(ncid, "ns", jm, nsDimID), "creating ns dimension") !jm
    call check(nf90_def_dim(ncid, "ch", nch, chaDimID), "creating nch dimension")

    call check(nf90_put_att(ncid,NF90_GLOBAL,'title','VLIDORT Simulation of GEOS-5 TEMPO Sampler'),"title attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'institution','NASA/Goddard Space Flight Center'),"institution attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'source','Global Model and Assimilation Office'),"source attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'history','VLIDORT simulation run from geo_vlidort.x'),"history attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'references','n/a'),"references attr")    
    call check(nf90_put_att(ncid,NF90_GLOBAL,'comment','This file contains VLIDORT simulated top of the atmosphere ' // &
                                                       'radiance and reflectance from GEOS-5 parameters sampled ' // &
                                                       ' on the TEMPO geostationary grid '),"comment attr")   
    call check(nf90_put_att(ncid,NF90_GLOBAL,"contact","Arlindo da Silva <arlindo.dasilva@nasa.gov>"),"contact attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","cf"),"conventions attr")

    dimids = (/ewDimID,nsDimID,timeDimID,chaDimID/)
    chunk_size = (/1,1,1,nch/)

    call check(nf90_def_var(ncid,'channel',nf90_float,(/chaDimID/),chaVarID),"create channel var")
    call check(nf90_def_var(ncid,'radiance',nf90_float,dimids,radVarID,chunksizes=chunk_size),"create radiance var")
    call check(nf90_def_var(ncid,'reflectance',nf90_float,dimids,refVarID,chunksizes=chunk_size),"create reflectance var")

    call check(nf90_put_att(ncid,chaVarID,'standard_name','Channel'),"standard_name attr")
    call check(nf90_put_att(ncid,chaVarID,'long_name','Channel Wavelength'),"long_name attr")
    call check(nf90_put_att(ncid,chaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,chaVarID,'units','nm'),"units attr")
    call check(nf90_put_att(ncid,chaVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

    call check(nf90_put_att(ncid,radVarID,'standard_name','Top of Atmosphere Radiance'),"standard_name attr")
    call check(nf90_put_att(ncid,radVarID,'long_name','Top of Atmosphere Radiance'),"long_name attr")
    call check(nf90_put_att(ncid,radVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,radVarID,'units','W m-2 sr-1 nm-1'),"units attr")
    call check(nf90_put_att(ncid,radVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

    call check(nf90_put_att(ncid,refVarID,'standard_name','Top of Atmosphere Reflectance'),"standard_name attr")
    call check(nf90_put_att(ncid,refVarID,'long_name','Top of Atmosphere Reflectance'),"long_name attr")
    call check(nf90_put_att(ncid,refVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,refVarID,'units','None'),"units attr")
    call check(nf90_put_att(ncid,refVarID,"_FillValue",real(MISSING)),"_Fillvalue attr")

    !Leave define mode
    call check(nf90_enddef(ncid),"leaving define mode")

    ! one processor writes out channels
    if (myid == 0) then
      call check(nf90_put_var(ncid, chaVarID, channels), "writing out channels")
    endif

    ! force independent write 
    call check(nf90_var_par_access(ncid, refVarID, nf90_independent),"reflectance writes independently")
    call check(nf90_var_par_access(ncid, radVarID, nf90_independent),"radiance writes independently")

!    call check( nf90_close(ncid), "close outfile" )

  end subroutine mp_create_outfile  


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     mp_write_outfile
! PURPOSE
!     write reflectance and radiance data to outfile
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     6 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine mp_write_outfile(ncid,radVarID,refVarID,start,count,radData,refData)
    integer,intent(in)                  :: ncid
    integer,intent(in)                  :: radVarID, refVarID    
    real*8,dimension(:,:),intent(in)    :: refData, radData
    integer,dimension(4),intent(in)     :: start, count

    call check(nf90_put_var(ncid, radVarID, radData, start = start, count = count), "writing out radiance")
    call check(nf90_put_var(ncid, refVarID, refData, start = start, count = count), "writing out reflectance")

  end subroutine mp_write_outfile


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     date2doy
! PURPOSE
!     calculate the day of the year from the date (Jan 1 = 1)
! INPUT
!     none
! OUTPUT
!     doy  : day of year
!  HISTORY
!     5 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine date2doy()
    integer                                  :: yyyy
    integer                                  :: mm 
    integer                                  :: dd
    integer, dimension(12)                   :: daysInMonth 
    integer                                  :: i               

    read(date(1:4),'(i)') yyyy
    read(date(5:6),'(i)') mm
    read(date(7:8),'(i)') dd

    daysInMonth = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/) 
    if (leap(yyyy)) then
      daysInMonth(1) = 29
    end if 

    doy = 0
    do i = 1, mm-1
      doy = doy + daysInMonth(i)
    end do
    doy = doy + dd

  end subroutine date2doy

!***********************************************************************************************************************************
!  LEAP
!
!  Input:
!     year  -  integer
!  Output:
!     Function return value = .true. if year is a leap year, and .false. otherwise.
!***********************************************************************************************************************************

  FUNCTION LEAP (YEAR) RESULT (LEAPFLAG)

    IMPLICIT NONE

    INTEGER :: YEAR
    LOGICAL :: LEAPFLAG

    LEAPFLAG = .FALSE.
    IF (MOD(YEAR,4) .EQ. 0)   LEAPFLAG = .TRUE.
    IF (MOD(YEAR,100) .EQ. 0) LEAPFLAG = .FALSE.
    IF (MOD(YEAR,400) .EQ. 0) LEAPFLAG = .TRUE.
    RETURN
  END FUNCTION

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     sensor_geometry
! PURPOSE
!     calculate the sensor viewing zenith and azimuth angles
! INPUT
!     lon : longitude of pixel center
!     n   : observation number
! OUTPUT
!     solar_zenith 
!     solar_azimuth
!  HISTORY
!     5 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine sensor_geometry(lat,lon,n)
    real, intent(in)                  :: lat, lon
    integer, intent(in)               :: n

    real                              :: delta_L
    real                              :: lambda
    real                              :: sin_rho
    real                              :: sensor_elev
    real                              :: cos_az
    real, parameter                   :: pi  = 2.0*acos(0.0)
    real, parameter                   :: d2r = pi/180.0
    real, parameter                   :: r2d = 180.0/pi

 
    delta_L = d2r*abs(sat_lon - lon)
    lambda = acos(sin(d2r*sat_lat)*sin(d2r*lat) + cos(d2r*sat_lat)*cos(d2r*lat)*cos(delta_L))
 
    cos_az = (sin(d2r*lat) - cos(lambda)*sin(d2r*sat_lat))/(sin(lambda)*cos(d2r*sat_lat))
    
    ! hack to deal with nadir pixels where numerical errors give slightly more than 1
    if (cos_az > 1.0)  cos_az = 1.0
    if (cos_az < -1.0) cos_az = -1.0

    sensor_azimuth = r2d*acos(cos_az)

    sin_rho = earth_rad/(earth_rad+sat_alt)

    sensor_zenith(n) = r2d*atan(sin_rho*sin(lambda)/(1-sin_rho*cos(lambda)))

  end subroutine sensor_geometry

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    solar_geometry
! PURPOSE
!     calculated solar geomtery for the pixel
! INPUT
!     lat, lon  : pixel location
! OUTPUT
!     solar_azimuth
!     solar_zenith
!  HISTORY
!     4 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine solar_geometry(lat,lon,n)
    real, intent(in)                  :: lat, lon
    integer, intent(in)               :: n            

    real, parameter                   :: pi  = 2.0*acos(0.0)
    real, parameter                   :: d2r = pi/180.0
    real, parameter                   :: r2d = 180.0/pi
    real                              :: lst
    real                              :: hr_ang
    real                              :: solar_declin
    real                              :: solar_elev
    real                              :: tt
    real                              :: D
    real                              :: ET
    real                              :: cos_az


    read(time,*) tt
    D      = 360*(doy-81)/365 
    ET     = 9.87*sin(2*D*d2r) - 7.53*cos(D*d2r) - 1.5*sin(D*d2r)  ! Equation of time in minutes needed for true solar time
    lst    = lon/15.0 + tt + ET/60.0

    if (lst < 0.0) then
      lst = 24.0 + lst
    else if (lst > 24.0) then
      lst = lst - 24.0
    end if 

    hr_ang = d2r*360*(lst - 12.0)/24.0

    solar_declin = 23.45*d2r*sin(D*d2r)

    solar_zenith(n) = r2d*acos(sin(abs(d2r*lat))*sin(solar_declin) + cos(d2r*abs(lat))*cos(solar_declin)*cos(hr_ang))

    solar_elev   = 90.0 - solar_zenith(n)

    cos_az = ((cos(d2r*lat)*sin(solar_declin)) - cos(hr_ang)*cos(solar_declin)*sin(d2r*lat))/(cos(d2r*solar_elev))

    ! hack to deal with nadir pixels where numerical errors give slightly more than 1
    if (cos_az > 1.0)  cos_az = 1.0
    if (cos_az < -1.0) cos_az = -1.0

    if (solar_elev > 0 ) then
      solar_azimuth = r2d*acos(cos_az)

      if (hr_ang > 0) then
        solar_azimuth = solar_azimuth + 180.0
      end if
    end if

  end subroutine solar_geometry

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

  subroutine calc_qm(var,q,n,i,j)
    real, intent(in), dimension(:,:,:)     :: var
    integer, intent(in)                    :: q, n, i, j

    integer                                :: k

    do k = 1, km
      qm(k,q,n) = var(i,j,k)*DELP(i,j,k)/grav
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
!    readvar2D
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
  subroutine readvar2D(varname, filename, var)
    character(len=*), intent(in)           ::  varname
    character(len=*), intent(in)           ::  filename
    real, dimension(:,:), intent(inout)    ::  var

    integer                                :: ncid, varid


!         write(*,'(A,A,A,I4)')'Reading ',trim(varname), ' on PE ', myid
    call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
    call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
    call check( nf90_get_var(ncid,varid,var), "reading " // varname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine readvar2D

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

    integer                                :: ncid, varid

    call check( nf90_open(filename,IOR(nf90_nowrite, nf90_mpiio),ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call check( nf90_get_att(ncid, NF90_GLOBAL, attrname, attrvar), "getting value for global attribute" // attrname)
    call check( nf90_close(ncid), "closing " // filename)
  end subroutine mp_readGattr


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
  subroutine mp_colreadvar(varname, filename, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    real, dimension(im,jm), intent(inout)  ::  var

    integer                       :: p, startc, countc, endc
    integer                       :: ncid, varid
    integer, dimension(npet)      :: ncol                ! how many columns each processor reads

    ! Everyone Figure out how many columns each PE has to read
    ! -----------------------------
    ncol = 0
    if (npet >= jm) then
      ncol(1:npet) = 1
    else if (npet < jm) then
      ncol(1:npet) = jm/npet
      ncol(npet)   = ncol(npet) + mod(jm,npet)
    end if 


    call check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    do p = 0, npet-1
      if (myid == p) then
        if (p == 0) then
          startc = 1
        else
          startc = sum(ncol(1:p))+1
        end if
        countc = ncol(p+1)
        endc   = startc + countc - 1

        call check( nf90_get_var(ncid, varid,var(:,startc:endc), start = (/ 1, startc, 1 /), count=(/im,countc,1/)), "reading " // varname)

      end if
    end do
    call check( nf90_close(ncid), "closing "// filename)
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
  subroutine mp_layreadvar(varname, filename, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    real, dimension(im,jm,km), intent(inout)  ::  var

    integer                       :: p, startl, countl, endl
    integer                       :: ncid, varid
    integer, dimension(npet)      :: nlayer                ! how many layers each processor reads

    ! Everyone Figure out how many layers each PE has to read
    ! -----------------------------
    nlayer = 0
    if (npet >= km) then
      nlayer(1:npet) = 1
    else if (npet < km) then
      nlayer(1:npet) = km/npet
      nlayer(npet)   = nlayer(npet) + mod(km,npet)
    end if 


    call check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
    do p = 0, npet-1
      if (myid == p) then
        if (p == 0) then
          startl = 1
        else
          startl = sum(nlayer(1:p))+1
        end if
        countl = nlayer(p+1)
        endl   = startl + countl - 1

        call check( nf90_get_var(ncid, varid,var(:,:,startl:endl), start = (/ 1, 1, startl, 1 /), count=(/im,jm,countl,1/)), "reading " // varname)

      end if
    end do
    call check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_layreadvar


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvar4D
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
  subroutine mp_readvar4D(varname, filename, dim, n, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, dimension(:), intent(in)         ::  dim
    integer, intent(in)                       ::  n
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

    call check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
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
          call check( nf90_get_var(ncid, varid,var(startl:endl,:,:,:), start = (/ startl, 1, 1, 1 /), count=countsize), "reading " // varname)
        else if (n == 2) then
          call check( nf90_get_var(ncid, varid,var(:,startl:endl,:,:), start = (/ 1, startl, 1, 1 /), count=countsize), "reading " // varname)
        else if (n == 3) then
          call check( nf90_get_var(ncid, varid,var(:,:,startl:endl,:), start = (/ 1, 1, startl, 1 /), count=countsize), "reading " // varname)
        else
          call check( nf90_get_var(ncid, varid,var(:,:,:,startl:endl), start = (/ 1, 1, 1, startl /), count=countsize), "reading " // varname)
        end if                 

      end if
    end do
    call check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_readvar4D

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    mp_readvar3D
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
  subroutine mp_readvar3D(varname, filename, dim, n, var)
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  filename
    integer, dimension(:), intent(in)         ::  dim
    integer, intent(in)                       ::  n  !dimensions to read over in chunks
    real*8, dimension(:,:,:), intent(inout)     ::  var

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

    call check( nf90_open(filename, IOR(nf90_nowrite, nf90_mpiio), ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL), "opening file " // filename)
    call check( nf90_inq_varid(ncid, varname, varid), "getting varid for " // varname)
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
          call check( nf90_get_var(ncid, varid,var(startl:endl,:,:), start = (/ startl, 1, 1 /), count=countsize), "reading " // varname)
        else if (n == 2) then
          call check( nf90_get_var(ncid, varid,var(:,startl:endl,:), start = (/ 1, startl, 1 /), count=countsize), "reading " // varname)
        else 
          call check( nf90_get_var(ncid, varid,var(:,:,startl:endl), start = (/ 1, 1, startl /), count=countsize), "reading " // varname)
        end if                 

      end if
    end do
    call check( nf90_close(ncid), "closing "// filename)
    !
  end subroutine mp_readvar3D  

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

    if ( myid  == 0 ) then
      open (unit = 2, file="shmem_test.txt")
      write(2,'(A)') '--- Array Statistics ---'
    end if

    call shmem_test2D('CLDTOT',CLDTOT)
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

    !   Wait for everyone to finish and print max memory used
    !   -----------------------------------------------------------  
    call MAPL_SyncSharedMemory(rc=ierr)
    if (myid == 0) then  
      write(*,*) 'Tested shared memory' 
      call sys_tracker()   
    end if   

  end subroutine do_testing



  subroutine shmem_test3D(varname,var)
    character(len=*), intent(in)  :: varname
    real,dimension(:,:,:)         :: var
    character(len=61)             :: msg

    if ( myid  == 0 ) then
      open (unit = 2, file="shmem_test.txt",position="append")
      do p = 1,npet-1
        call mpi_recv(msg, 61, MPI_CHARACTER, p, 1, MPI_COMM_WORLD, status, ierr)
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

    if ( myid  == 0 ) then
      open (unit = 2, file="shmem_test.txt",position="append")
      do p = 1,npet-1
        call mpi_recv(msg, 61, MPI_CHARACTER, p, 1, MPI_COMM_WORLD, status, ierr)
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
    end if

  end subroutine check

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    shutdown
! PURPOSE
!     Clean up allocated arrays and finalize MPI
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine shutdown()         

    deallocate (pe)
    deallocate (ze)
    deallocate (te)
    deallocate (qm)
    deallocate (tau)
    deallocate (ssa)
    deallocate (g)

    deallocate (nclr)

    ! shmem must deallocate shared memory arrays
    ! call MAPL_DeallocNodeArray(CLDTOT,rc=ierr)
    ! call MAPL_DeallocNodeArray(AIRDENS,rc=ierr)
    ! call MAPL_DeallocNodeArray(RH,rc=ierr)
    ! call MAPL_DeallocNodeArray(DELP,rc=ierr)
    ! call MAPL_DeallocNodeArray(DU001,rc=ierr)
    ! call MAPL_DeallocNodeArray(DU002,rc=ierr)
    ! call MAPL_DeallocNodeArray(DU003,rc=ierr)
    ! call MAPL_DeallocNodeArray(DU004,rc=ierr)                           
    ! call MAPL_DeallocNodeArray(DU005,rc=ierr)
    ! call MAPL_DeallocNodeArray(SS001,rc=ierr) 
    ! call MAPL_DeallocNodeArray(SS002,rc=ierr) 
    ! call MAPL_DeallocNodeArray(SS003,rc=ierr) 
    ! call MAPL_DeallocNodeArray(SS004,rc=ierr) 
    ! call MAPL_DeallocNodeArray(SS005,rc=ierr) 
    ! call MAPL_DeallocNodeArray(BCPHOBIC,rc=ierr) 
    ! call MAPL_DeallocNodeArray(BCPHILIC,rc=ierr) 
    ! call MAPL_DeallocNodeArray(OCPHOBIC,rc=ierr) 
    ! call MAPL_DeallocNodeArray(OCPHILIC,rc=ierr) 

    !call MAPL_FinalizeShmem (rc=ierr)

    call MPI_Finalize(ierr)

  end subroutine shutdown


end program shmem_reader
