!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    geo_vlidort
! PURPOSE
!     Reads in parallel the model data (in a netcdf file) interpolated to TEMPO grid 
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

program geo_vlidort

  use MAPL_ShmemMod                ! The SHMEM infrastructure
  use netcdf                       ! for reading the NR files
  use mpi
  use vlidort_brdf_modis_py        ! Module to run VLIDORT with MODIS BRDF surface supplement
  use geo_vlidort_netcdf

  implicit none
  include "geo_vlidort_pars.F90"
  !include "mpif.h"

! Date 
! --------
  character(len=8)      :: date 
  character(len=2)      :: time 

! Test flag
! -----------
  logical               :: test_shmem = .False.

! File names
! ----------
  character(len=256)    :: MET_file, AER_file, INV_file, BRDF_file, LAND_file, OUT_file 

! Global, 3D arrays to be allocated using SHMEM
! ---------------------------------------------
  real, pointer                         :: CLON(:,:) => null()
  real, pointer                         :: CLAT(:,:) => null()
  real, pointer                         :: CLDTOT(:,:) => null()
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
  real, pointer                         :: KISO(:,:,:) => null()
  real, pointer                         :: KVOL(:,:,:) => null()
  real, pointer                         :: KGEO(:,:,:) => null()
  real, pointer                         :: FROCEAN(:,:) => null()
  real, pointer                         :: FRLAKE(:,:) => null()
  real, pointer                         :: FRLAND(:,:) => null()
  real, pointer                         :: FRLANDICE(:,:) => null()  

! VLIDORT input arrays
! ---------------------------
  real, pointer                         :: pe(:,:) => null()      ! edge pressure [Pa]
  real, pointer                         :: ze(:,:) => null()      ! edge height above sfc [m]
  real, pointer                         :: te(:,:) => null()      ! edge Temperature [K]

  real, pointer                         :: qm(:,:,:) => null()    ! (mixing ratio) * delp/g
  real, pointer                         :: tau(:,:,:) => null()   ! aerosol optical depth
  real, pointer                         :: ssa(:,:,:) => null()   ! single scattering albedo
  real, pointer                         :: g(:,:,:) => null()     ! asymmetry factor
  real, pointer                         :: albedo(:,:) => null()  ! surface albedo
  real, pointer                         :: solar_zenith(:) => null() 
  real, pointer                         :: relat_azimuth(:) => null()
  real, pointer                         :: sensor_zenith(:) => null()

! VLIDORT output arrays
!-------------------------------
  real*8,pointer                        :: radiance_VL(:,:) => null()      ! TOA normalized radiance from VLIDORT
  real*8,pointer                        :: reflectance_VL(:,:) => null()  ! TOA reflectance from VLIDORT
  real*8,pointer                        :: radiance_VL_Surface(:,:) => null()      ! TOA normalized radiance from VLIDORT
  real*8,pointer                        :: reflectance_VL_Surface(:,:) => null()  ! TOA reflectance from VLIDORT  

! VLIDORT working variables
!------------------------------
  real                                  :: sensor_azimuth
  real                                  :: solar_azimuth
  real                                  :: sat_lon                   ! satellite longitude
  real                                  :: sat_lat                   ! satellite latitude
  real                                  :: sat_alt                   ! satellite altitude
  real                                  :: Earth_rad                 ! earth radius
  integer                               :: doy                       ! day of year
  integer                               :: ch                        ! i-channel  

! MODIS Kernel  arrays
!--------------------------
  real*8, pointer                       :: kernel_wt(:,:,:) => null()  ! kernel weights (/fiso,fgeo,fvol/)
  real*8, pointer                       :: param(:,:,:) => null()      ! Li-Sparse parameters 
                                                                 ! param1 = crown relative height (h/b)
                                                                 ! param2 = shape parameter (b/r)
  real                                  :: modis_missing                                                                 


! Satellite domain variables
!------------------------------
  integer                               :: im, jm, km                ! size of TEMPO domain
  integer                               :: i, j, k, n                ! TEMPO domain working variable
  integer                               :: starti, counti, endi      ! array indices and counts for e-w columns to be read
  integer, allocatable                  :: nclr(:)                   ! how many clear pixels each processor reads
  integer                               :: clrm                      ! number of clear pixels
  integer                               :: c                         ! clear pixel working variables
  integer,allocatable,dimension(:)      :: i_work, j_work            ! clear pixel working indeces

! netcdf variables
!----------------------  
  integer                               :: ncid, radVarID, refVarID  ! netcdf ids

! Miscellaneous
! -------------
  integer                               :: ierr                      ! MPI error message
  integer                               :: status(MPI_STATUS_SIZE)   ! MPI status
  integer                               :: myid, npet, CoresPerNode  ! MPI dimensions and processor id
  integer                               :: p                         ! i-processor
  character(len=100)                    :: msg                       ! message to be printed
  real                                  :: progress 
  real                                  :: g5nr_missing

! System tracking variables
! -----------------------------
  integer                               :: t1, t2, clock_max
  real*8                                :: clock_rate


! Start Timing
! -----------------------------
  call system_clock ( t1, clock_rate, clock_max )


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

! For now hard code file names 
! ----------------------------------------------------------------------------
  date = '20051231'
  time = '18'

! INFILES
  write(MET_file,'(A,A,A,A,A,A,A,A,A)') '/nobackup/TEMPO/met_Nv/Y',date(1:4),'/M',date(5:6),'/tempo-g5nr.lb2.met_Nv.',date,'_',time,'z.nc4'
  write(AER_file,'(A,A,A,A,A,A,A,A,A)') '/nobackup/TEMPO/aer_Nv/Y',date(1:4),'/M',date(5:6),'/tempo-g5nr.lb2.aer_Nv.',date,'_',time,'z.nc4'
  INV_file = "/nobackup/TEMPO/LevelG/invariant/tempo.lg1.invariant.nc4"
  BRDF_file = "/nobackup/TEMPO/BRDF/raw/MAIACRTLS.2006008.hdf"
  LAND_file = "/nobackup/TEMPO/const_2d_asm_Nx/Y2005/M12/tempo-g5nr.lb2.asm_Nx.20051231_00z.nc4"

! Query for domain dimensions, satellite constants, and missing value
!--------------------------------------------------------------
  call mp_readDim("ew", MET_file, im)
  call mp_readDim("ns", MET_file, jm)
  call mp_readDim("lev", MET_file, km)
  call mp_readGattr("sat_lat", INV_FILE, sat_lat)
  call mp_readGattr("sat_lon", INV_FILE, sat_lon)
  call mp_readGattr("sat_alt", INV_FILE, sat_alt)
  call mp_readGattr("Earth_radius", INV_FILE, earth_rad) 
  call mp_readVattr("missing_value", BRDF_FILE, "Kiso", modis_missing) 
  call mp_readVattr("missing_value", MET_FILE, "CLDTOT", g5nr_missing)

! OUTFILES
  write(OUT_file,'(A,A,A,A,A)') 'tempo-g5nr.lb2.vlidort.',date,'_',time,'z.nc4'
  call mp_create_outfile(ncid,radVarID,refVarID)

  call date2doy()

 
! Allocate the Global arraya using SHMEM
! It will be available on all processors
! ---------------------------------------------------------
  call allocate_shared()


! Allocate arrays that will be copied on each processor - unshared
! ---------------------------------------------------------
  call allocate_unshared()


! Read the cloud, land and BRDF data
! -------------------------------------
  call mp_colreadvar("CLDTOT", MET_file, npet, myid, CLDTOT)
  call mp_colreadvar("clat", INV_FILE, npet, myid, CLAT)
  call mp_colreadvar("clon", INV_FILE, npet, myid, CLON)
  call mp_colreadvar("FROCEAN", LAND_file, npet, myid, FROCEAN)
  call mp_colreadvar("FRLAKE", LAND_file,  npet, myid, FRLAKE)
  call mp_colreadvar("FRLAND", LAND_file,  npet, myid, FRLAND)
  call mp_colreadvar("FRLANDICE", LAND_file,  npet, myid, FRLANDICE)  
  call mp_readvar3D("Kiso", BRDF_file, (/im,jm,nbands/), 1, npet, myid, KISO)
  call mp_readvar3D("Kvol", BRDF_file, (/im,jm,nbands/), 1, npet, myid, KVOL)
  call mp_readvar3D("Kgeo", BRDF_file, (/im,jm,nbands/), 1, npet, myid, KGEO)
  
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
  ! first check that you can do anything!!!!!
  if (.not. ANY(CLDTOT <= cldmax)) then
    write(*,*) 'domain is too cloudy, nothing to do'
    write(*,*) 'Exiting.....'
    call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
  endif

  ! figure out how many indices to work on
  ! right now skipping over erroneous missing values
  clrm = 0
  do i=1,im
    do j=1,jm
      if (FROCEAN(i,j) .ne. g5nr_missing) then
        if (CLDTOT(i,j) <= cldmax) then
          clrm = clrm + 1
        end if
      end if 
    end do
  end do

  if (clrm == 0) then
    write(*,*) 'no clear pixels over land, nothing to do'
    write(*,*) 'Exiting.....'
    call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
  endif

  ! store them
  allocate(i_work(clrm))
  allocate(j_work(clrm))

  clrm = 0
  do i=1,im
    do j=1,jm
      if (FROCEAN(i,j) .ne. g5nr_missing) then
        if (CLDTOT(i,j) <= cldmax) then
          clrm = clrm + 1
          i_work(clrm) = i
          j_work(clrm) = j
        end if
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

! Read the aerosol variables layer-by-layer in parallel
! ------------------------------
  call read_aer_Nv()

! Wait for everyone to finish reading and print max memory used
! -----------------------------------------------------------  
  call MAPL_SyncSharedMemory(rc=ierr)
  if (myid == 0) then 
    write(*,*) 'Read all variables' 
    call sys_tracker()   
  end if     

! Although read on individual PEs, all shared variables should have the same
! data in all PEs. Let's verify that.
! ----------------------------------------------------------- 
  if (test_shmem) call do_testing()

! Prepare inputs and run VLIDORT
! ------------------------------
  call strarr_2_chararr(vnames_string,nq,16,vnames)
  
  do p = 1, npet-1
    if (myid == p) then
      if (p == 0) then
        starti = 1
      else
        starti = sum(nclr(1:p))+1
      end if
      counti = nclr(p+1)
      endi   = starti + counti - 1

      do c = starti, endi
        call getEdgeVars ( km, nobs, reshape(AIRDENS(i_work(c),j_work(c),:),(/km,nobs/)), &
                           reshape(DELP(i_work(c),j_work(c),:),(/km,nobs/)), ptop, &
                           pe, ze, te )   

        write(msg,'(A,I)') 'getEdgeVars ', myid
        call write_verbose(msg)
        
        call calc_qm(DU001,1,nobs,i_work(c),j_work(c))
        call calc_qm(DU002,2,nobs,i_work(c),j_work(c))
        call calc_qm(DU003,3,nobs,i_work(c),j_work(c))
        call calc_qm(DU004,4,nobs,i_work(c),j_work(c))
        call calc_qm(DU005,5,nobs,i_work(c),j_work(c))
        call calc_qm(SS001,6,nobs,i_work(c),j_work(c))
        call calc_qm(SS002,7,nobs,i_work(c),j_work(c))
        call calc_qm(SS003,8,nobs,i_work(c),j_work(c))
        call calc_qm(SS004,9,nobs,i_work(c),j_work(c))
        call calc_qm(SS005,10,nobs,i_work(c),j_work(c))
        call calc_qm(BCPHOBIC,11,nobs,i_work(c),j_work(c))
        call calc_qm(BCPHILIC,12,nobs,i_work(c),j_work(c))
        call calc_qm(OCPHOBIC,13,nobs,i_work(c),j_work(c))
        call calc_qm(OCPHILIC,14,nobs,i_work(c),j_work(c))

        write(msg,'(A,I)') 'calc_qm ', myid
        call write_verbose(msg)

        call sensor_geometry(CLAT(i_work(c),j_work(c)),CLON(i_work(c),j_work(c)),nobs)

        write(msg,*) 'sensor_geometry ', myid, sensor_zenith, sensor_azimuth
        call write_verbose(msg)

        call solar_geometry(CLAT(i_work(c),j_work(c)),CLON(i_work(c),j_work(c)),nobs)

        write(msg,*) 'solar_geometry ', myid, solar_zenith, solar_azimuth
        call write_verbose(msg)

        if (solar_zenith(nobs) < 90.0) then
          relat_azimuth(nobs) = abs(sensor_azimuth - solar_azimuth)        

          do ch = 1, nch
            kernel_wt(:,ch,nobs) = (/dble(KISO(i_work(c),j_work(c),bands(ch))),&
                                  dble(KGEO(i_work(c),j_work(c),bands(ch))),&
                                  dble(KVOL(i_work(c),j_work(c),bands(ch)))/)
            param(:,ch,nobs)     = (/dble(2),dble(1)/)
          end do
 
          call getAOPscalar ( km, nobs, nch, nq, rcfile, channels, vnames, verbose, &
                         qm, reshape(RH(i_work(c),j_work(c),:),(/km,nobs/)), &
                         tau, ssa, g, ierr )

          write(msg,*) 'getAOPscalar ', myid
          call write_verbose(msg)

          ! Figure out what kinds of surfaces to simulate
          ! Weight the radiances according to the surface area fractions
          !------------------------------------------------------------------

          ! WATER
          !---------------------
          !if (FROCEAN(i_work(c),j_work(c)) + FRLAKE(i_work(c),j_work(c)) .ne. 0 ) then
            ! here i calculate water leaving surface calculations
            ! also need to do something about FRSEAICE
            !if ((FRLAKE(i_work(c),j_work(c)) > 0) .or. ((FROCEAN(i_work(c),j_work(c))-FRSEAICE(i_work(c),j_work(c))) > 0) then
              ! do some water calcs!
            !end if
            !if (FRSEAICE(i_work(c),j_work(c)) > 0)
              !do some ice calcs!
            !end if 
            ! radiance_VL = radiance_VL + radiance_VL_Surface*(FROCEAN(i_work(c),j_work(c))-FRSEAICE(i_work(c),j_work(c))) 
            ! radiance_VL = radiance_VL + radiance_VL_Surface*FRLAKE(i_work(c),j_work(c))

            ! reflectance_VL = reflectance_VL + reflectance_VL_Surface*(FROCEAN(i_work(c),j_work(c))-FRSEAICE(i_work(c),j_work(c)))
            ! reflectance_VL = reflectance_VL + reflectance_VL_Surface*FRLAKE(i_work(c),j_work(c))

          !end if 

          ! LAND
          !---------------------          
          if (FRLAND(i_work(c),j_work(c)) .ne. 0 ) then
            ! do land leaving surface calculations
            ! need to implement snow cover (FRSNO)
            if (ANY(kernel_wt == modis_missing)) then
              albedo(:,:)      = 0.05  !this needs to be a climatology(?)

              ! Call to vlidort scalar code
              ! Simple lambertian surface model
              !------------------------------
              ! call Scalar (km, nch, nobs ,dble(channels),        &
              !         dble(tau), dble(ssa), dble(g), dble(pe), dble(ze), dble(te), dble(albedo),&
              !         dble(solar_zenith), dble(relat_azimuth), dble(sensor_zenith), &
              !         dble(MISSING),verbose,radiance_VL_Surface,reflectance_VL_Surface, ierr)
            else 
              ! Call to vlidort scalar code
              ! MODIS BRDF Surface Model
              !------------------------------
              call Scalar_LandMODIS (km, nch, nobs, dble(channels),        &
                      dble(tau), dble(ssa), dble(g), dble(pe), dble(ze), dble(te), &
                      kernel_wt, param, &
                      dble(solar_zenith), dble(relat_azimuth), dble(sensor_zenith), &
                      dble(MISSING),verbose,radiance_VL_Surface,reflectance_VL_Surface,ierr )  

            end if

            ! Done Calling VLIDORT - check that everything exits OK
            !--------------------------------
            call mp_check_vlidort(radiance_VL_Surface,reflectance_VL_Surface)

            ! radiance_VL    = radiance_VL + radiance_VL_Surface*FRLAND(i_work(c),j_work(c))
            ! reflectance_VL = reflectance_VL + reflectance_VL_Surface*FRLAND(i_work(c),j_work(c))

            radiance_VL    = radiance_VL_Surface
            reflectance_VL = reflectance_VL_Surface
          end if  


          ! LAND PERMANENTLY COVERED BY ICE
          !---------------------          
          !if (FRLANDICE(i_work(c),j_work(c)) .ne. 0) then
          !  write(*,*) 'We found landice!!!!???'
          !end if
          
          ! Done with VLIDORT Simulations
          !--------------------------------
          write(msg,*) 'VLIDORT scalar DONE', myid, ierr
          call write_verbose(msg)

          ! Write output to correct position in file
          ! --------------------------------------------
          call mp_write_outfile(ncid,radVarID,refVarID,(/i_work(c),j_work(c),1,1/),(/1,1,nobs,nch/),radiance_VL,reflectance_VL)

          ! Keep track of progress of each processor
          ! -----------------------------------------        
          progress = 100.*real(c-starti)/real(counti)
          write(*,*) 'geo_vlidort ',myid,' ',c,endi, nint(progress),'%'
        
        end if      
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
      call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
    else if (ierr > 0) then
      write(*,*) 'VLIDORT returned rc code ', ierr
      write(*,*) 'Exiting......'
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
    call mp_layreadvar("AIRDENS", AER_file, npet, myid, AIRDENS)
    call mp_layreadvar("RH", AER_file, npet, myid, RH)
    call mp_layreadvar("DELP", AER_file, npet, myid, DELP)
    call mp_layreadvar("DU001", AER_file, npet, myid, DU001)
    call mp_layreadvar("DU002", AER_file, npet, myid, DU002)
    call mp_layreadvar("DU003", AER_file, npet, myid, DU003)
    call mp_layreadvar("DU004", AER_file, npet, myid, DU004)
    call mp_layreadvar("DU005", AER_file, npet, myid, DU005)
    call mp_layreadvar("SS001", AER_file, npet, myid, SS001)
    call mp_layreadvar("SS002", AER_file, npet, myid, SS002)
    call mp_layreadvar("SS003", AER_file, npet, myid, SS003)
    call mp_layreadvar("SS004", AER_file, npet, myid, SS004)
    call mp_layreadvar("SS005", AER_file, npet, myid, SS005)
    call mp_layreadvar("BCPHOBIC", AER_file, npet, myid, BCPHOBIC)
    call mp_layreadvar("BCPHILIC", AER_file, npet, myid, BCPHILIC)
    call mp_layreadvar("OCPHOBIC", AER_file, npet, myid, OCPHOBIC)
    call mp_layreadvar("OCPHILIC", AER_file, npet, myid, OCPHILIC)
  end subroutine read_aer_Nv

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     allocate_shared
! PURPOSE
!     allocates all the shared memory arrays that I need
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine allocate_shared()

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
    call MAPL_AllocNodeArray(FROCEAN,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(FRLAKE,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(FRLAND,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(FRLANDICE,(/im,jm/),rc=ierr)

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
    allocate (albedo(nch,nobs))

    allocate (solar_zenith(nobs))
    allocate (relat_azimuth(nobs))
    allocate (sensor_zenith(nobs))

    allocate (radiance_VL(nobs,nch))
    allocate (reflectance_VL(nobs, nch))

    allocate (radiance_VL_Surface(nobs,nch))
    allocate (reflectance_VL_Surface(nobs, nch))    

    radiance_VL = 0
    reflectance_VL = 0

    allocate (kernel_wt(nkernel,nch,nobs))
    allocate (param(nparam,nch,nobs))

  ! Needed for reading
  ! ----------------------
    allocate (nclr(npet)) 
    nclr = 0
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


end program geo_vlidort
