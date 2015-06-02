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
  use vlidort_brdf_modis           ! Module to run VLIDORT with MODIS BRDF surface supplement
  use geo_vlidort_netcdf           ! Module with netcdf routines
  use GeoAngles                    ! Module with geostationary satellite algorithms for scene geometry

  implicit none
  include "geo_vlidort_pars.F90"
  !include "mpif.h"

! rcfile variables
!----------------
  character(len=256)                    :: arg, rcline, rcvar, rcvalue  
  integer                               :: io

! RC Inputs 
! --------
  character(len=8)                      :: date 
  character(len=7)                      :: brdfdate
  character(len=2)                      :: time 
  character(len=256)                    :: satname, indir, outdir, brdfname
  logical                               :: scalar

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

! VLIDORT working variables
!------------------------------
  integer                               :: ch                        ! i-channel  
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
  integer                               :: ierr                      ! MPI error message
  integer                               :: status(MPI_STATUS_SIZE)   ! MPI status
  integer                               :: myid, npet, CoresPerNode  ! MPI dimensions and processor id
  integer                               :: p                         ! i-processor
  character(len=100)                    :: msg                       ! message to be printed
  real                                  :: progress 
  real                                  :: g5nr_missing

! System tracking variables
! -----------------------------
  integer*8                             :: t1, t2, clock_max
  real*8                                :: clock_rate


! Start Timing
! -----------------
  call system_clock ( t1, clock_rate, clock_max )
  progress = -1


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

! Parse Resource file for input info 
! ----------------------------------------------------------------------------
  call getarg(1, arg)

  call get_config()


! Write out settings to use
! -------------------------------
  if (myid == 0) then
    write(*,*) 'Simulating ', lower_to_upper(trim(satname)),' domain on ',date,' ', time, 'Z'
    write(*,*) 'Input directory: ',trim(indir)
    write(*,*) 'Output directory: ',trim(outdir)
    write(*,*) 'BRDF dataset: ',trim(brdfname),' ',trim(brdfdate)
    if (scalar) write(*,*) 'Scalar calculations'
    if (.not. scalar) write(*,*) 'Vector calculations'
    write(*,*) ' '
  end if 

! INFILES & OUTFILE set names
! -------------------------------
  call filenames()

! Query for domain dimensions and missing value
!--------------------------------------------------------------
  call mp_readDim("ew", MET_file, im)
  call mp_readDim("ns", MET_file, jm)
  call mp_readDim("lev", MET_file, km)
  call mp_readVattr("missing_value", BRDF_FILE, "Kiso", modis_missing) 
  call mp_readVattr("missing_value", MET_FILE, "CLDTOT", g5nr_missing)


! Create OUTFILE
! --------------------
  call mp_create_outfile(date,time,ncid,radVarID,refVarID)

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

  if (.not. ANY(SOLAR_ZENITH <= szamax)) then   
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
          if (SOLAR_ZENITH(i,j) <= szamax) then
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
  
  do p = 0, npet-1
    if (myid == p) then
      if (p == 0) then
        starti = 1
      else
        starti = sum(nclr(1:p))+1
      end if
      counti = nclr(p+1)
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

        do ch = 1, nch
          kernel_wt(:,ch,nobs) = (/dble(KISO(c,bands(ch))),&
                                dble(KGEO(c,bands(ch))),&
                                dble(KVOL(c,bands(ch)))/)
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

        if ((lower_to_upper(BRDFNAME) == 'MAIACRTLS') .or. (ANY(kernel_wt == modis_missing))) then
          ! Simple lambertian surface model
          !------------------------------
          albedo(:,:)      = 0.05  !this needs to be a climatology(?)
          if (scalar) then
            ! Call to vlidort scalar code
            
            call Scalar_Lambert (km, nch, nobs ,dble(channels),        &
                    dble(tau), dble(ssa), dble(g), dble(pe), dble(ze), dble(te), dble(albedo),&
                    (/dble(SZA(c))/), &
                    (/dble(RAA(c))/), &
                    (/dble(VZA(c))/), &
                    dble(MISSING),verbose,radiance_VL_Surface,reflectance_VL_Surface, ierr)
          else
            ! Call to vlidort vector code
            call Vector_Lambert (km, nch, nobs ,dble(channels), nMom,   &
                   nPol, dble(tau), dble(ssa), dble(g), dble(pmom), dble(pe), dble(ze), dble(te), dble(albedo),&
                   (/dble(SZA(c))/), &
                   (/dble(RAA(c))/), &
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
                    (/dble(RAA(c))/), &
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
                    (/dble(RAA(c))/), &
                    (/dble(VZA(c))/), &
                    dble(MISSING),verbose,radiance_VL_Surface,reflectance_VL_Surface, Q, U, ierr )  

          end if

          call mp_check_vlidort(radiance_VL_Surface,reflectance_VL_Surface)

          ! radiance_VL    = radiance_VL + radiance_VL_Surface*FRLAND(i_work(c),j_work(c))
          ! reflectance_VL = reflectance_VL + reflectance_VL_Surface*FRLAND(i_work(c),j_work(c))
        end if          
          
        radiance_VL    = radiance_VL_Surface
        reflectance_VL = reflectance_VL_Surface

        write(msg,*) 'VLIDORT Calculations DONE', myid, ierr
        call write_verbose(msg)

        ! Write output to correct position in file
        ! --------------------------------------------
        call mp_write_vlidort(ncid,radVarID,refVarID,(/i_work(c),j_work(c),1,1/),(/1,1,nobs,nch/),radiance_VL,reflectance_VL)

        ! Keep track of progress of each processor
        ! -----------------------------------------        
        if (nint(100.*real(c-starti)/real(counti)) > progress) then
          progress = nint(100.*real(c-starti)/real(counti))
          write(*,'(A,I,A,I,A,I2,A,I3,A)') 'Pixel: ',c,'  End Pixel: ',endi,'  ID:',myid,'  Progress:', nint(progress),'%'           
        end if
                    
      end do
    end if
  end do 

! Everyone close OUT_file and sync up
! -----------------------------------------
500  call check( nf90_close(ncid), "close outfile" )

  if (myid == 0) then
    write(*,*) '<> Finished VLIDORT Simulation of '//trim(lower_to_upper(satname))//' domain'
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

  ! INFILES
  write(MET_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/',trim(satname),'-g5nr.lb2.met_Nv.',date,'_',time,'z.nc4'
  write(AER_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/',trim(satname),'-g5nr.lb2.aer_Nv.',date,'_',time,'z.nc4'
  write(ANG_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/',trim(satname),'.lb2.angles.',date,'_',time,'z.nc4'
  write(INV_file,'(4A)')  trim(indir),'/LevelG/invariant/',trim(satname),'.lg1.invariant.nc4'
  write(BRDF_file,'(6A)') trim(indir),'/BRDF/raw/',trim(brdfname),'.',brdfdate,'.hdf'
  write(LAND_file,'(4A)') trim(indir),'/LevelB/invariant/',trim(satname),'-g5nr.lb2.asm_Nx.nc4' 

! OUTFILES
  write(OUT_file,'(8A)') trim(outdir),'/',trim(satname),'-g5nr.lb2.vlidort.',date,'_',time,'z.nc4'
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
!     upper_to_lower
! PURPOSE
!     converts uppercase characters to lowercase
! INPUT
!     word
! OUTPUT
!     upper_to_lower
!  HISTORY
!     18 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   

  function upper_to_lower ( word )
    implicit none

    character(len=*), intent(in)   :: word
    character(len=:),allocatable   :: upper_to_lower
    integer                        :: i,n

    character(*), parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
    character(*), parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    allocate(character(len=len(word)) :: upper_to_lower)

    upper_to_lower = word

    ! Loop over string elements
    do i = 1, LEN(word)
    ! Find location of letter in lower case constant string
      n = INDEX(UPPER_CASE, word( i:i ))
    ! If current substring is a lower case letter, make it upper case
      if ( n /= 0 ) upper_to_lower( i:i ) = LOWER_CASE( n:n )
    end do

  end function upper_to_lower  

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

    if (myid == 0) then
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

    allocate (radiance_VL(nobs,nch))
    allocate (reflectance_VL(nobs, nch))

    allocate (radiance_VL_Surface(nobs,nch))
    allocate (reflectance_VL_Surface(nobs, nch))    

    radiance_VL = 0
    reflectance_VL = 0

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

  subroutine mp_create_outfile(date, time, ncid, radVarID, refVarID)
    character(len=*)                   :: date, time
    integer,intent(out)                :: ncid
    integer,intent(out)                :: radVarID, refVarID    
    
    integer, dimension(4)              :: chunk_size
    integer, dimension(4)              :: dimids
    integer                            :: timeDimID, ewDimID, nsDimID, chaDimID
        integer                        :: szaVarID, vzaVarID, raaVarID    
    integer                            :: timeVarID, clonVarID, clatVarID, chaVarID
    real,allocatable,dimension(:,:)    :: clon, clat, sza, vza, raa
    real,allocatable,dimension(:)      :: scantime


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

    call check(nf90_def_var(ncid,'scanTime',nf90_float,(/ewDimID/),timeVarID),"create scanTime var")
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

    ! one processor writes out channels and clon, clat, scantime
    if (myid == 0) then
      call check(nf90_put_var(ncid, chaVarID, channels), "writing out channels")

      allocate (scantime(im))
      allocate (clon(im, jm))
      allocate (clat(im, jm))

      call readvar1D("scanTime", MET_file, scantime)
      call check(nf90_put_var(ncid,timeVarID,scantime), "writing out scantime")

      call readvar2D("clon", INV_file, clon)
      call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

      call readvar2D("clon", INV_file, clat)
      call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

      deallocate (clon)
      deallocate (clat)  
      deallocate (scantime)

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

    endif

    ! force independent write 
    call check(nf90_var_par_access(ncid, refVarID, nf90_independent),"reflectance writes independently")
    call check(nf90_var_par_access(ncid, radVarID, nf90_independent),"radiance writes independently")

    !call check( nf90_close(ncid), "close outfile" )

  end subroutine mp_create_outfile  


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     mp_write_vlidort
! PURPOSE
!     write reflectance and radiance data to outfile
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     6 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine mp_write_vlidort(ncid,radVarID,refVarID,start,count,radData,refData)
    integer,intent(in)                  :: ncid
    integer,intent(in)                  :: radVarID, refVarID    
    real*8,dimension(:,:),intent(in)    :: refData, radData
    integer,dimension(4),intent(in)     :: start, count

    call check(nf90_put_var(ncid, radVarID, radData, start = start, count = count), "writing out radiance")
    call check(nf90_put_var(ncid, refVarID, refData, start = start, count = count), "writing out reflectance")

  end subroutine mp_write_vlidort

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     write_outfile3D
! PURPOSE
!     write general data to outfile by one processor
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     18 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine write_outfile3D(ncid,VarID,myData)
    integer,intent(in)                  :: ncid
    integer,intent(in)                  :: VarID    
    real,dimension(:,:,:),intent(in)    :: myData

    call check(nf90_put_var(ncid, VarID, myData ), "writing out myData")

  end subroutine write_outfile3D

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     write_outfile1D
! PURPOSE
!     write general data to outfile by one processor
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     18 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine write_outfile1D(ncid,VarID,myData)
    integer,intent(in)                  :: ncid
    integer,intent(in)                  :: VarID    
    real,dimension(:),intent(in)        :: myData

    call check(nf90_put_var(ncid, VarID, myData ), "writing out myData")

  end subroutine write_outfile1D

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
!    get_config()
! PURPOSE
!     Read and parse resource file
!     fixed format
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     29 May P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine get_config()
    open(unit=5, file = arg, IOSTAT=io, mode = 'READ', status = 'OLD')
    if (io > 0) then
      if (myid == 0) then
        write(*,*) 'Problem opening geo_vlidort rcfile'
        write(*,*) 'Exiting.....'
        call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
      end if
    end if
    do 
      read(5,'(A)',IOSTAT=io) rcline
      if (io > 0) then
        if (myid == 0) then
          write(*,*) 'Problem reading geo_vlidort rcfile'
          write(*,*) 'Exiting.....'
          call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)
        end if
      else if (io < 0) then
        ! end of file
        exit
      else 
        !parse line
        i = index(rcline, ':')  
        if (i > 0) then    
          rcvar = lower_to_upper(rcline(1:i-1))
          rcvalue = adjustl(rcline(i+1:))

          if (trim(rcvar) .eq. 'DATE') date = trim(rcvalue)
          if (trim(rcvar) .eq. 'TIME') time = trim(rcvalue)
          if (trim(rcvar) .eq. 'SATNAME') satname = upper_to_lower(trim(rcvalue))
          if (trim(rcvar) .eq. 'INDIR') indir = trim(rcvalue) 
          if (trim(rcvar) .eq. 'OUTDIR') outdir = trim(rcvalue)
          if (trim(rcvar) .eq. 'BRDFNAME') brdfname = trim(rcvalue)
          if (trim(rcvar) .eq. 'BRDFDATE') brdfdate = trim(rcvalue)
          if (trim(rcvar) .eq. 'SCALAR') then
            if (lower_to_upper(trim(rcvalue)) .eq. 'TRUE') scalar = .TRUE.
            if (lower_to_upper(trim(rcvalue)) .eq. 'FALSE') scalar = .FALSE.
          end if
        end if 
      end if
    end do 
    close(5)

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

    ! deallocate (pe)
    ! deallocate (ze)
    ! deallocate (te)
    ! deallocate (qm)
    ! deallocate (tau)
    ! deallocate (ssa)
    ! deallocate (g)

    ! deallocate (nclr)

    ! shmem must deallocate shared memory arrays
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
