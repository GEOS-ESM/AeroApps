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
#  include "MAPL_Generic.h"
#  include "MAPL_ErrLogMain.h"
program geo_vlidort

  use ESMF                         ! ESMF modules
  use MAPL_Mod
  use MAPL_ShmemMod                ! The SHMEM infrastructure
  use netcdf                       ! for reading the NR files
  use vlidort_brdf_modis           ! Module to run VLIDORT with MODIS BRDF surface supplement
  use geo_vlidort_netcdf           ! Module with netcdf routines
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
  character(len=8)                      :: date         
  character(len=7)                      :: brdfdate
  character(len=2)                      :: time 
  character(len=256)                    :: instname, indir, outdir, brdfname
  character(len=256)                    :: brdfband
  integer, dimension(:),allocatable     :: brdfband_i                  ! modis band indeces that overlap with vlidort channels
  real, dimension(:),allocatable        :: brdfband_c                  ! modis band center wavelength
  logical                               :: scalar
  real, dimension(:), allocatable       :: channels               ! channels to simulate
  integer                               :: nch                    ! number of channels  
  real                                  :: cldmax                 ! Cloud Filtering  
  real                                  :: szamax                 ! Geomtry filtering

! Test flag
! -----------
  logical                               :: test_shmem = .False.

! File names
! ----------
  character(len=256)                    :: MET_file, AER_file, ANG_file, INV_file, BRDF_file, LAND_file, OUT_file 

! Global, 3D arrays to be allocated using SHMEM
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
  real, pointer                         :: KISO(:,:) => null()
  real, pointer                         :: KVOL(:,:) => null()
  real, pointer                         :: KGEO(:,:) => null()
  real, pointer                         :: SZA(:) => null()
  real, pointer                         :: VZA(:) => null()
  real, pointer                         :: RAA(:) => null()  

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

! VLIDORT output arrays
!-------------------------------
  real*8,pointer                        :: radiance_VL(:,:) => null()             ! TOA normalized radiance from VLIDORT
  real*8,pointer                        :: reflectance_VL(:,:) => null()          ! TOA reflectance from VLIDORT
  real*8,pointer                        :: radiance_VL_Surface(:,:) => null()     ! TOA normalized radiance from VLIDORT
  real*8,pointer                        :: reflectance_VL_Surface(:,:) => null()  ! TOA reflectance from VLIDORT  
  real*8,pointer                        :: Q(:,:) => null()                       ! Q Stokes component
  real*8,pointer                        :: U(:,:) => null()                       ! U Stokes component
  real*8,pointer                        :: field(:,:) => null()

! VLIDORT working variables
!------------------------------
  integer                               :: ch                        ! i-channel  
  integer                               :: iband                     ! i-brdfband
  real,allocatable                      :: pmom(:,:,:,:,:)           ! elements of scattering phase matrix for vector calculations

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
  real, pointer                         :: CLDTOT(:,:) => null()     ! GEOS-5 cloud fraction
  real, pointer                         :: FRLAND(:,:) => null()     ! GEOS-5 land fraction
  real, pointer                         :: SOLAR_ZENITH(:,:) => null() 
  integer,allocatable,dimension(:)      :: i_work, j_work            ! clear pixel working indeces
  logical,allocatable,dimension(:,:)    :: clmask                    ! cloud-land mask

! netcdf variables
!----------------------  
  integer                               :: ncid, radVarID, refVarID  ! netcdf ids

! Miscellaneous
! -------------
  integer                               :: ierr, rc, status          ! MPI error message
  integer                               :: status_mpi(MPI_STATUS_SIZE)   ! MPI status
  integer                               :: myid, npet, CoresPerNode  ! MPI dimensions and processor id
  integer                               :: p                         ! i-processor
  character(len=100)                    :: msg                       ! message to be printed
  real                                  :: progress 
  real                                  :: g5nr_missing

! System tracking variables
! -----------------------------
  integer*8                             :: t1, t2, clock_max
  real*8                                :: clock_rate
  character(len=*), parameter           :: Iam = 'geo_vlidort'

  
! Start Timing
! -----------------
  call system_clock ( t1, clock_rate, clock_max )
  progress = -1

! Initialize MPI with ESMF
! --------------
  call ESMF_Initialize (logkindflag=ESMF_LOGKIND_NONE, vm=vm, __RC__)

  call ESMF_VMGet(vm, localPET=myid, PETcount=npet) 
  if ( MAPL_am_I_root() ) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'

! Initialize SHMEM
! ----------------
  CoresPerNode = MAPL_CoresPerNodeGet(MPI_COMM_WORLD,rc=ierr) ! a must
  call MAPL_InitializeShmem(rc=ierr)

! Parse Resource file provided at command line for input info 
! ----------------------------------------------------------------------------
  call getarg(1, arg)
  call get_config(arg)

! Write out settings to use
! -------------------------------
  if (myid == 0) then
    write(*,*) 'Simulating ', lower_to_upper(trim(instname)),' domain on ',date,' ', time, 'Z'
    write(*,*) 'Input directory: ',trim(indir)
    write(*,*) 'Output directory: ',trim(outdir)
    write(*,*) 'BRDF dataset: ',trim(brdfname),' ',trim(brdfdate)
    write(*,*) 'Channels [nm]: ',channels
    if (lower_to_upper(brdfband) == 'EXACT') write(*,*) 'Using exact BRDF on bands : ',brdfband_i
    if (lower_to_upper(brdfband) == 'INTERPOLATE') write(*,*) 'Using interpolated BRDF' 
    write(*,*) 'Cloud Fraction <= ', cldmax
    write(*,*) 'SZA < ', szamax
    if (scalar) write(*,*) 'Scalar calculations'
    if (.not. scalar) write(*,*) 'Vector calculations'
    write(*,*) ' '
  end if 


! Query for domain dimensions and missing value
!--------------------------------------------------------------
  call mp_readDim("ew", MET_file, im)
  call mp_readDim("ns", MET_file, jm)
  call mp_readDim("lev", MET_file, km)
  call mp_readVattr("missing_value", BRDF_FILE, "Kiso", modis_missing) 
  call mp_readVattr("missing_value", MET_FILE, "CLDTOT", g5nr_missing)

! Create OUTFILE
! --------------------
  if ( MAPL_am_I_root() )  call create_outfile(date,time,ncid,radVarID,refVarID)

! Allocate arrays that will be copied on each processor - unshared
! ---------------------------------------------------------
  call allocate_unshared()

! Read the cloud, land, and angle data 
! -------------------------------------
  call read_cloud()
  call read_land()
  call read_sza()

! Create cloud-land mask
! ------------------------
  ! first check that you can do anything!!!!!  
  if (.not. ANY(CLDTOT <= cldmax)) then   
    if (myid == 0) then
      write(*,*) 'The domain is too cloudy, nothing to do'
      write(*,*) 'Exiting.....'
    end if
    GOTO 500
  end if

  if (.not. ANY(SOLAR_ZENITH < szamax)) then   
    if (myid == 0) then
      write(*,*) 'The sun has set, nothing to do'
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
            clrm = clrm + 1
            clmask(i,j) = .True.
          end if
        end if
      end if 
    end do
  end do

  if (clrm == 0) then
    if (myid == 0) then
      write(*,*) 'no clear pixels over land, nothing to do'
      write(*,*) 'Exiting.....'
    end if
    GOTO 500
  end if

! Cloud and Land data no longer needed
!---------------------------------------
  deallocate (CLDTOT)
  deallocate (FRLAND)
  deallocate (SOLAR_ZENITH)

! Store indices
!-------------------
  allocate(i_work(clrm))
  allocate(j_work(clrm))

  clrm = 0
  do i=1,im
    do j=1,jm
      if (clmask(i,j)) then
        clrm = clrm + 1
        i_work(clrm) = i
        j_work(clrm) = j
      end if 
    end do
  end do       

  if (myid == 0) then
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
 call read_BRDF()
 call read_angles()

! Wait for everyone to finish reading and print max memory used
! ------------------------------------------------------------------  
  call MAPL_SyncSharedMemory(rc=ierr)
   
  if (myid == 0) then 
    write(*,*) 'Read all variables' 
    call sys_tracker()   
    write(*,*) ' '
    call system_clock ( t2, clock_rate, clock_max )
    write ( *, * ) 'Elapsed real time = ', real ( t2 - t1 ) / real ( clock_rate )
    write(*,*) ' '
  end if   

! Split up filtered domain among processors
!----------------------------------------------
  if (npet >= clrm) then
    nclr(1:npet) = 1
  else if (npet < clrm) then
    nclr(1:npet) = clrm/npet
    nclr(npet)   = nclr(npet) + mod(clrm,npet)
  end if 

! Although read on individual PEs, all shared variables should have the same
! data in all PEs. Let's verify that.
! ----------------------------------------------------------- 
  if (test_shmem) call do_testing()

! Prepare inputs and run VLIDORT
! -----------------------------------
  call strarr_2_chararr(vnames_string,nq,16,vnames)
  
  if (myid == 0) then
    starti = 1
  else
    starti = sum(nclr(1:myid))+1
  end if
  counti = nclr(myid+1)
  endi   = starti + counti - 1

  do c = starti, endi
    call getEdgeVars ( km, nobs, reshape(AIRDENS(c,:),(/km,nobs/)), &
                       reshape(DELP(c,:),(/km,nobs/)), ptop, &
                       pe, ze, te )   

    write(msg,'(A,I)') 'getEdgeVars ', myid
    call write_verbose(msg)
    
    call calc_qm()

    write(msg,'(A,I)') 'calc_qm ', myid
    call write_verbose(msg)  

    ! Get BRDF Kernel Weights
    !----------------------------
    do ch = 1, nch
      if (lower_to_upper(brdfband) == 'EXACT') then
        kernel_wt(:,ch,nobs) = (/dble(KISO(c,brdfband_i(ch))),&
                            dble(KGEO(c,brdfband_i(ch))),&
                            dble(KVOL(c,brdfband_i(ch)))/)
      else
        ! > 2130 uses highest MODIS wavelength band     
        if (channels(ch) >= 2130) then  
          iband = minloc(abs(brdfband_c - channels(ch)), dim = 1)
          kernel_wt(:,ch,nobs) = (/dble(KISO(c,iband)),&
                            dble(KGEO(c,iband)),&
                            dble(KVOL(c,iband))/)
        end if
        
        if (channels(ch) < 2130) then
          ! nearest neighbor interpolation of kernel weights to wavelength
          ! ******has to fall above available range
          kernel_wt(1,ch,nobs) = dble(nn_interp(brdfband_c,reshape(KISO(c,:),(/nbands/)),channels(ch)))
          kernel_wt(2,ch,nobs) = dble(nn_interp(brdfband_c,reshape(KGEO(c,:),(/nbands/)),channels(ch)))
          kernel_wt(3,ch,nobs) = dble(nn_interp(brdfband_c,reshape(KVOL(c,:),(/nbands/)),channels(ch)))          
        end if
      end if
      param(:,ch,nobs)     = (/dble(2),dble(1)/)
    end do

    if (scalar) then
      call getAOPscalar ( km, nobs, nch, nq, rcfile, channels, vnames, verbose, &
                          qm, reshape(RH(c,:),(/km,nobs/)), &
                          tau, ssa, g, ierr )
    else
      call getAOPvector ( km, nobs, nch, nq, rcfile, channels, vnames, verbose, &
                          qm, reshape(RH(c,:),(/km,nobs/)),&
                          nMom,nPol, tau, ssa, g, pmom, ierr )
    end if

    write(msg,*) 'getAOP ', myid
    call write_verbose(msg)

    if ((lower_to_upper(BRDFNAME) /= 'MAIACRTLS') .or. (ANY(kernel_wt == modis_missing))) then
      ! Simple lambertian surface model
      !------------------------------
      albedo(:,:)      = 0.05  !this needs to be a climatology(?)
      if (scalar) then
        ! Call to vlidort scalar code
        
        call Scalar_Lambert (km, nch, nobs ,dble(channels),        &
                dble(tau), dble(ssa), dble(g), dble(pe), dble(ze), dble(te), dble(albedo),&
                (/dble(SZA(c))/), &
                (/dble(abs(RAA(c)))/), &
                (/dble(VZA(c))/), &
                dble(MISSING),verbose,radiance_VL_Surface,reflectance_VL_Surface, ierr)
      else
        ! Call to vlidort vector code
        call Vector_Lambert (km, nch, nobs ,dble(channels), nMom,   &
               nPol, dble(tau), dble(ssa), dble(g), dble(pmom), dble(pe), dble(ze), dble(te), dble(albedo),&
               (/dble(SZA(c))/), &
               (/dble(abs(RAA(c)))/), &
               (/dble(VZA(c))/), &
               dble(MISSING),verbose,radiance_VL_Surface,reflectance_VL_Surface, Q, U, ierr)
      end if
    else             
      ! MODIS BRDF Surface Model
      !------------------------------
      if (scalar) then 
        ! Call to vlidort scalar code            
        call Scalar_LandMODIS (km, nch, nobs, dble(channels),        &
                dble(tau), dble(ssa), dble(g), dble(pe), dble(ze), dble(te), &
                kernel_wt, param, &
                (/dble(SZA(c))/), &
                (/dble(abs(RAA(c)))/), &
                (/dble(VZA(c))/), &
                dble(MISSING),verbose,radiance_VL_Surface,reflectance_VL_Surface,ierr )  
      else
        ! Call to vlidort vector code
        write(msg,*) 'getting ready to do vector calculations', myid, ierr
        call write_verbose(msg)

        call Vector_LandMODIS (km, nch, nobs, dble(channels), nMom, &
                nPol, dble(tau), dble(ssa), dble(g), dble(pmom), dble(pe), dble(ze), dble(te), &
                kernel_wt, param, &
                (/dble(SZA(c))/), &
                (/dble(abs(RAA(c)))/), &
                (/dble(VZA(c))/), &
                dble(MISSING),verbose,radiance_VL_Surface,reflectance_VL_Surface, Q, U, ierr )  

      end if

      call mp_check_vlidort(radiance_VL_Surface,reflectance_VL_Surface)

    end if          
      
    radiance_VL(c,:)    = radiance_VL_Surface(nobs,:)
    reflectance_VL(c,:) = reflectance_VL_Surface(nobs,:)

    write(msg,*) 'VLIDORT Calculations DONE', myid, ierr
    call write_verbose(msg)

    ! Keep track of progress of each processor
    ! -----------------------------------------        
    if (nint(100.*real(c-starti)/real(counti)) > progress) then
      progress = nint(100.*real(c-starti)/real(counti))
      write(*,'(A,I,A,I,A,I2,A,I3,A)') 'Pixel: ',c,'  End Pixel: ',endi,'  ID:',myid,'  Progress:', nint(progress),'%'           
    end if
                
  end do

  ! Wait for everyone to finish calculations
  ! ------------------------------------------
  call MAPL_SyncSharedMemory(rc=ierr)

  if (MAPL_am_I_root()) then
    ! Expand radiance to im x jm using mask
    !-------------------------------
    allocate (field(im,jm))
    field = g5nr_missing

    call check( nf90_open(OUT_file, nf90_write, ncid), "opening file " // OUT_file )
    ! Write output to correct position in file
    !  --------------------------------------------
    do ch = 1, nch
      call check(nf90_put_var(ncid, radVarID, unpack(reshape(radiance_VL(:,ch),(/clrm/)),clmask,field), &
                  start = (/1,1,1,ch/), count = (/im,jm,nobs,1/)), "writing out radiance")
      call check(nf90_put_var(ncid, refVarID, unpack(reshape(reflectance_VL(:,ch),(/clrm/)),clmask,field), &
                  start = (/1,1,1,ch/), count = (/im,jm,nobs,1/)), "writing out reflectance")
    end do

    call check( nf90_close(ncid), "close outfile" )
    deallocate(field)
  end if

500 call MAPL_SyncSharedMemory(rc=ierr)
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
  call MAPL_SyncSharedMemory(rc=ierr)
  call shutdown()

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

  below = minloc(abs(xint - x), dim = 1, mask = (xint - x) .LE. 0)
  above = minloc(abs(xint - x), dim = 1, mask = (xint - x) .GT. 0)

  
  if (.not. ANY((/y(above),y(below)/) == modis_missing)) then
    top = y(above) - y(below)
    bottom = x(above) - x(below)
    nn_interp = y(below) + (xint-x(below)) * top / bottom
  else
    nn_interp  = modis_missing
  end if
end function nn_interp

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
  character(len=256)          :: chmax, chmin

  ! INFILES
  write(MET_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/',trim(instname),'-g5nr.lb2.met_Nv.',date,'_',time,'z.nc4'
  write(AER_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/',trim(instname),'-g5nr.lb2.aer_Nv.',date,'_',time,'z.nc4'
  write(ANG_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/',trim(instname),'.lb2.angles.',date,'_',time,'z.nc4'
  write(INV_file,'(4A)')  trim(indir),'/LevelG/invariant/',trim(instname),'.lg1.invariant.nc4'
  write(BRDF_file,'(6A)') trim(indir),'/BRDF/raw/',trim(brdfname),'.',brdfdate,'.hdf'
  write(LAND_file,'(4A)') trim(indir),'/LevelB/invariant/',trim(instname),'-g5nr.lb2.asm_Nx.nc4' 

  write(chmax,'(I4)') int(maxval(channels))
  write(chmin,'(I4)') int(minval(channels))

! OUTFILES
  write(OUT_file,'(4A)') trim(outdir),'/',trim(instname),'-g5nr.lb2.vlidort.'

  if (scalar) then 
    write(OUT_file,'(2A)') trim(OUT_file),'scalar.'
  else
    write(OUT_file,'(2A)') trim(OUT_file),'vector.'
  end if 
  
  if (lower_to_upper(brdfname) == 'MAIACRTLS' .and. lower_to_upper(brdfband) == 'INTERPOLATE') then
    write(OUT_file,'(2A)') trim(OUT_file),'iMAIACRTLS.'
  else if (lower_to_upper(brdfname) == 'MAIACRTLS' .and. lower_to_upper(brdfband) == 'EXACT') then
    write(OUT_file,'(2A)') trim(OUT_file),'MAIACRTLS.'
  else if (lower_to_upper(brdfname) /= 'MAIACRTLS') then
    write(OUT_file,'(2A)') trim(OUT_file),'lambertian.'
  end if

  if (nch == 1) then    
    write(OUT_file,'(7A)') trim(OUT_file),date,'_',time,'z_',trim(adjustl(chmin)),'nm.nc4'
  else 
    write(OUT_file,'(9A)') trim(OUT_file),date,'_',time,'z_',trim(adjustl(chmin)),'-',trim(adjustl(chmax)),'nm.nc4'
  end if

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
    real, dimension(im,jm,km)   :: temp

    if (myid == 0) then
      call readvar3D("AIRDENS", AER_file, temp)
      call reduceProfile(temp,clmask,AIRDENS)      

      call readvar3D("RH", AER_file, temp)
      call reduceProfile(temp,clmask,RH)

      call readvar3D("DELP", AER_file, temp)
      call reduceProfile(temp,clmask,DELP)

      call readvar3D("DU001", AER_file, temp)
      call reduceProfile(temp,clmask,DU001)

      call readvar3D("DU002", AER_file, temp)
      call reduceProfile(temp,clmask,DU002)

      call readvar3D("DU003", AER_file, temp)
      call reduceProfile(temp,clmask,DU003)

      call readvar3D("DU004", AER_file, temp)
      call reduceProfile(temp,clmask,DU004)

      call readvar3D("DU005", AER_file, temp)
      call reduceProfile(temp,clmask,DU005)

      call readvar3D("SS001", AER_file, temp)
      call reduceProfile(temp,clmask,SS001)

      call readvar3D("SS002", AER_file, temp)
      call reduceProfile(temp,clmask,SS002)

      call readvar3D("SS003", AER_file, temp)
      call reduceProfile(temp,clmask,SS003)

      call readvar3D("SS004", AER_file, temp)
      call reduceProfile(temp,clmask,SS004)

      call readvar3D("SS005", AER_file, temp)
      call reduceProfile(temp,clmask,SS005)

      call readvar3D("BCPHOBIC", AER_file, temp)
      call reduceProfile(temp,clmask,BCPHOBIC) 

      call readvar3D("BCPHILIC", AER_file, temp)
      call reduceProfile(temp,clmask,BCPHILIC)

      call readvar3D("OCPHOBIC", AER_file, temp)
      call reduceProfile(temp,clmask,OCPHOBIC) 

      call readvar3D("OCPHILIC", AER_file, temp)
      call reduceProfile(temp,clmask,OCPHILIC)    

      write(*,*) '<> Read aeorosl data to shared memory'  
    end if
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
!     read_BRDF
! PURPOSE
!     read in all the BRDF variables
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     15 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine read_BRDF()
    real, dimension(im,jm,nbands)   :: temp

    if (MAPL_am_I_root()) then
      call readvar3D("Kiso", BRDF_file, temp)
      call reduceProfile(temp,clmask,KISO)

      call readvar3D("Kvol", BRDF_file, temp)
      call reduceProfile(temp,clmask,KVOL)

      call readvar3D("Kgeo", BRDF_file, temp)
      call reduceProfile(temp,clmask,KGEO)

      write(*,*) '<> Read BRDF data to shared memory' 
    end if 
  end subroutine read_BRDF  

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

    if (myid == 0) then
      call readvar2D("solar_zenith", ANG_file, temp)
      SZA = pack(temp,clmask)

      call readvar2D("sensor_zenith", ANG_file, temp)
      VZA = pack(temp,clmask)

      call readvar2D("relat_azimuth", ANG_file, temp)
      RAA = pack(temp,clmask)
      write(*,*) '<> Read angle data to shared memory' 
    end if

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
    call MAPL_AllocNodeArray(KISO,(/clrm,nbands/),rc=ierr)
    call MAPL_AllocNodeArray(KVOL,(/clrm,nbands/),rc=ierr)
    call MAPL_AllocNodeArray(KGEO,(/clrm,nbands/),rc=ierr)
    call MAPL_AllocNodeArray(SZA,(/clrm/),rc=ierr)
    call MAPL_AllocNodeArray(VZA,(/clrm/),rc=ierr)
    call MAPL_AllocNodeArray(RAA,(/clrm/),rc=ierr)

    call MAPL_AllocNodeArray(radiance_VL,(/clrm,nch/),rc=ierr)
    call MAPL_AllocNodeArray(reflectance_VL,(/clrm,nch/),rc=ierr)

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

    allocate (radiance_VL_Surface(nobs,nch))
    allocate (reflectance_VL_Surface(nobs, nch))    

    allocate (kernel_wt(nkernel,nch,nobs))
    allocate (param(nparam,nch,nobs))

    if (.not. scalar) then
      allocate (pmom(km,nch,nobs,nMom,nPol))
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

    if ( myid  == 0 ) then
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
!     ncid        : netcdf file id
!     radVarID, refVarID : radiance and reflectance variable IDs
!  HISTORY
!     6 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine create_outfile(date, time, ncid, radVarID, refVarID)
    character(len=*)                   :: date, time
    integer,intent(out)                :: ncid
    integer,intent(out)                :: radVarID, refVarID    
    
    integer, dimension(4)              :: chunk_size
    integer, dimension(4)              :: dimids
    integer                            :: timeDimID, ewDimID, nsDimID, chaDimID
    integer                            :: szaVarID, vzaVarID, raaVarID    
    integer                            :: timeVarID, clonVarID, clatVarID, chaVarID
    integer                            :: nPolVarID, nMomVarID, ewVarID, nsVarID
    real,allocatable,dimension(:,:)    :: clon, clat, sza, vza, raa
    real,allocatable,dimension(:)      :: scantime, ew, ns


    call check(nf90_create(OUT_file, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // OUT_file)
    call check(nf90_def_dim(ncid, "time", 1, timeDimID), "creating time dimension")
    call check(nf90_def_dim(ncid, "ew", im, ewDimID), "creating ew dimension") !im
    call check(nf90_def_dim(ncid, "ns", jm, nsDimID), "creating ns dimension") !jm
    call check(nf90_def_dim(ncid, "ch", nch, chaDimID), "creating nch dimension")

    call check(nf90_put_att(ncid,NF90_GLOBAL,'title','VLIDORT Simulation of GEOS-5 '//lower_to_upper(trim(instname))//' Sampler'),"title attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'institution','NASA/Goddard Space Flight Center'),"institution attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'source','Global Model and Assimilation Office'),"source attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'history','VLIDORT simulation run from geo_vlidort.x'),"history attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'references','n/a'),"references attr")    
    call check(nf90_put_att(ncid,NF90_GLOBAL,'comment','This file contains VLIDORT simulated top of the atmosphere ' // &
                                                       'radiance and reflectance from GEOS-5 parameters sampled ' // &
                                                       ' on the '//lower_to_upper(trim(instname))//' geostationary grid '),"comment attr")   
    if (lower_to_upper(brdfname) == 'MAIACRTLS' .and. lower_to_upper(brdfband) == 'INTERPOLATE') then
      call check(nf90_put_att(ncid,NF90_GLOBAL,'surface_comment','MAIAC RTLS surface BRDF kernel weights interpolated to channels'),"surface_comment")
    else if (lower_to_upper(brdfname) == 'MAIACRTLS' .and. lower_to_upper(brdfband) == 'EXACT') then
      call check(nf90_put_att(ncid,NF90_GLOBAL,'surface_comment','MAIAC RTLS surface BRDF kernel weights without inerpolation to channel'),"surface_comment")
    else if (lower_to_upper(brdfname) /= 'MAIACRTLS') then
      call check(nf90_put_att(ncid,NF90_GLOBAL,'surface_comment','Lambertian surface reflectance'),"surface_comment")
    end if

    if ( scalar ) then
      call check(nf90_put_att(ncid,NF90_GLOBAL,'vlidort_comment','Scalar calculations'),"vlidort_comment")
    else
      call check(nf90_put_att(ncid,NF90_GLOBAL,'vlidort_comment','Vector calculations'),"surface_comment")
      call check(nf90_def_var(ncid,'nPol',nf90_float,(/chaDimID/),nPolVarID),"create scanTime var")
      call check(nf90_def_var(ncid,'nMom',nf90_float,(/chaDimID/),nMomVarID),"create scanTime var")

      call check(nf90_put_att(ncid,nPolVarID,'long_name','number of components of the scattering matrix'),"long_name attr")
      call check(nf90_put_att(ncid,nMomVarID,'long_name','number of phase function moments'),"long_name attr")
    end if

    call check(nf90_put_att(ncid,NF90_GLOBAL,"contact","Patricia Castellanos <patricia.castellanos@nasa.gov>"),"contact attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","cf"),"conventions attr")

    dimids = (/ewDimID,nsDimID,timeDimID,chaDimID/)
    chunk_size = (/1,1,1,nch/)

    call check(nf90_def_var(ncid,'scanTime',nf90_float,(/ewDimID/),timeVarID),"create scanTime var")
    call check(nf90_def_var(ncid,'ew',nf90_float,(/ewDimID/),ewVarID),"create ew var")
    call check(nf90_def_var(ncid,'ns',nf90_float,(/nsDimID/),nsVarID),"create ns var")
    call check(nf90_def_var(ncid,'clon',nf90_float,(/ewDimID,nsDimID/),clonVarID),"create clon var")
    call check(nf90_def_var(ncid,'clat',nf90_float,(/ewDimID,nsDimID/),clatVarID),"create clat var")
    call check(nf90_def_var(ncid,'channel',nf90_float,(/chaDimID/),chaVarID),"create channel var")
    call check(nf90_def_var(ncid,'radiance',nf90_float,dimids,radVarID,chunksizes=chunk_size),"create radiance var")
    call check(nf90_def_var(ncid,'reflectance',nf90_float,dimids,refVarID,chunksizes=chunk_size),"create reflectance var")
    call check(nf90_def_var(ncid,'solar_zenith',nf90_float,(/ewDimID,nsDimID/),szaVarID),"create solar_zenith var")
    call check(nf90_def_var(ncid,'sensor_zenith',nf90_float,(/ewDimID,nsDimID/),vzaVarID),"create sensor_zenith var")
    call check(nf90_def_var(ncid,'relat_azimuth',nf90_float,(/ewDimID,nsDimID/),raaVarID),"create relat_azimuth var")

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


    call check(nf90_put_att(ncid,timeVarID,'long_name','Initial Time of Scan'),"long_name attr")
    call check(nf90_put_att(ncid,timeVarID,'units','seconds since '//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '// &
                                                  time//':00:00'),"units attr")

    call check(nf90_put_att(ncid,ewVarID,'long_name','pseudo longitude'),"long_name attr")
    call check(nf90_put_att(ncid,ewVarID,'units','degrees_east'),"units attr")
    call check(nf90_put_att(ncid,nsVarID,'long_name','pseudo latitude'),"long_name attr")
    call check(nf90_put_att(ncid,nsVarID,'units','degrees_north'),"units attr")   

    call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
    call check(nf90_put_att(ncid,clonVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
    call check(nf90_put_att(ncid,clatVarID,'missing_value',real(MISSING)),"missing_value attr")  

    call check(nf90_put_att(ncid,szaVarID,'long_name','solar zenith angle'),"long_name attr")
    call check(nf90_put_att(ncid,szaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,szaVarID,'units','degrees'),"units attr")

    call check(nf90_put_att(ncid,vzaVarID,'long_name','sensor zenith angle'),"long_name attr")
    call check(nf90_put_att(ncid,vzaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,vzaVarID,'units','degrees'),"units attr") 

    call check(nf90_put_att(ncid,raaVarID,'long_name','relative azimuth angle'),"long_name attr")
    call check(nf90_put_att(ncid,raaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,raaVarID,'units','degrees'),"units attr")       

    !Leave define mode
    call check(nf90_enddef(ncid),"leaving define mode")

    ! write out channels and clon, clat, scantime
    call check(nf90_put_var(ncid, chaVarID, channels), "writing out channels")

    allocate (scantime(im))
    allocate (clon(im, jm))
    allocate (clat(im, jm))
    allocate (ew(im))
    allocate (ns(jm))    

    call readvar1D("scanTime", MET_file, scantime)
    call check(nf90_put_var(ncid,timeVarID,scantime), "writing out scantime")

    call readvar2D("clon", INV_file, clon)
    call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

    call readvar2D("clat", INV_file, clat)
    call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

    call readvar1D("ew", MET_file, ew)
    call check(nf90_put_var(ncid,ewVarID,ew), "writing out ew")

    call readvar1D("ns", MET_file, ns)
    call check(nf90_put_var(ncid,nsVarID,ns), "writing out ns")   

    deallocate (clon)
    deallocate (clat)  
    deallocate (scantime)
    deallocate (ns)
    deallocate (ew)

    allocate (sza(im,jm))
    allocate (vza(im,jm))
    allocate (raa(im,jm))

    call readvar2D("solar_zenith", ANG_file, sza)
    call check(nf90_put_var(ncid,szaVarID,sza), "writing out sza")

    call readvar2D("sensor_zenith", ANG_file, vza)
    call check(nf90_put_var(ncid,vzaVarID,vza), "writing out vza")

    call readvar2D("relat_azimuth", ANG_file, raa)
    call check(nf90_put_var(ncid,raaVarID,raa), "writing out raa")

    deallocate (sza)
    deallocate (vza)
    deallocate (raa)

    if (.not. scalar) then
      call check(nf90_put_var(ncid,nPolVarID,nPol), "writing out nPol")
      call check(nf90_put_var(ncid,nMomVarID,nMom), "writing out nMom")
    end if

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

    !   Wait for everyone to finish and print max memory used
    !   -----------------------------------------------------------  
    call MAPL_SyncSharedMemory(rc=ierr)
    if (myid == 0) then  
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

    if ( myid  == 0 ) then
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

    if ( myid  == 0 ) then
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
  subroutine get_config(rcfile)
    character(len=*),intent(in)      :: rcfile

    cf = ESMF_ConfigCreate()
    call ESMF_ConfigLoadFile(cf, fileName=trim(rcfile), __RC__)

    ! Read in variables
    ! -----------------------
    call ESMF_ConfigGetAttribute(cf, date, label = 'DATE:',__RC__)
    call ESMF_ConfigGetAttribute(cf, time, label = 'TIME:',__RC__)
    call ESMF_ConfigGetAttribute(cf, instname, label = 'INSTNAME:',__RC__)
    call ESMF_ConfigGetAttribute(cf, indir, label = 'INDIR:',__RC__)
    call ESMF_ConfigGetAttribute(cf, outdir, label = 'OUTDIR:',default=indir)
    call ESMF_ConfigGetAttribute(cf, brdfname, label = 'BRDFNAME:',default='MAIACRTLS')
    call ESMF_ConfigGetAttribute(cf, brdfdate, label = 'BRDFDATE:',__RC__)
    call ESMF_ConfigGetAttribute(cf, scalar, label = 'SCALAR:',default=.TRUE.)
    call ESMF_ConfigGetAttribute(cf, szamax, label = 'SZAMAX:',default=90.0)
    call ESMF_ConfigGetAttribute(cf, cldmax, label = 'CLDMAX:',default=0.01)
    call ESMF_ConfigGetAttribute(cf, brdfband, label = 'BRDFBAND:', default='INTERPOLATE')

    ! Figure out number of channels and read into vector
    !------------------------------------------------------
    nch =  ESMF_ConfigGetLen(cf, label = 'CHANNELS:',__RC__)
    allocate (channels(nch))
    call ESMF_ConfigGetAttribute(cf, channels, label = 'CHANNELS:', default=550.)

    ! INFILES & OUTFILE set names
    ! -------------------------------
    call filenames()

    if (lower_to_upper(brdfband) == 'EXACT' ) then
      allocate (brdfband_i(nch))
      call ESMF_ConfigGetAttribute(cf, brdfband_i, label = 'BRDFBAND_I:', __RC__)
    else
      allocate (brdfband_c(nbands))
      call ESMF_ConfigGetAttribute(cf, brdfband_c, label = 'BRDFBAND_C:', __RC__)
    end if

    call ESMF_ConfigDestroy(cf)

  end subroutine get_config

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
     call MAPL_DeallocNodeArray(KISO,rc=ierr) 
     call MAPL_DeallocNodeArray(KVOL,rc=ierr) 
     call MAPL_DeallocNodeArray(KGEO,rc=ierr) 
     call MAPL_DeallocNodeArray(SZA,rc=ierr) 
     call MAPL_DeallocNodeArray(VZA,rc=ierr) 
     call MAPL_DeallocNodeArray(RAA,rc=ierr) 
     call MAPL_DeallocNodeArray(radiance_VL,rc=ierr) 
     call MAPL_DeallocNodeArray(reflectance_VL,rc=ierr) 

    call MAPL_FinalizeShmem (rc=ierr)

    call ESMF_Finalize(__RC__)

  end subroutine shutdown


end program geo_vlidort
