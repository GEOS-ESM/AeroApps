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
program geo_vlidort

  use ESMF                         ! ESMF modules
  use MAPL_Mod
  use MAPL_ShmemMod                ! The SHMEM infrastructure
  use netcdf                       ! for reading the NR files
  use vlidort_lambert              ! Module to run VLIDORT with lambertian surface  
  use vlidort_rot                  ! Module to calculate Rayleigh  
  use Chem_MieMod
!  use netcdf_helper                ! Module with netcdf routines
  use mp_netcdf_Mod
  use netcdf_Mod

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
  character(len=2)                      :: time 
  character(len=256)                    :: instname, indir, outdir
  logical                               :: scalar
  logical                               :: plane_parallel
  integer                               :: nstreams               ! number of half space streams, default = 6
  real, allocatable                     :: channels(:)            ! channels to simulate
  integer                               :: nch                    ! number of channels  
  real                                  :: cldmax                 ! Cloud Filtering  
  real                                  :: szamax, vzamax         ! Geomtry filtering
  integer                               :: nodemax                ! number of nodes requested
  integer                               :: nodenum                ! which node is this?
  character(len=256)                    :: version
  character(len=256)                    :: layout 


! Test flag
! -----------
  logical                               :: test_shmem = .False.

! File names
! ----------
  character(len=256)                    :: MET_file, AER_file, ANG_file, INV_file, LAND_file, OUT_file 

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
  real*8, allocatable                   :: depol(:)                               ! rayleigh depolarization ratio
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
  real, pointer                         :: PE_(:,:) => null()

  real*8,allocatable                    :: field(:,:)                             ! Template for unpacking shared arrays
  real*8,allocatable                    :: AOD(:)                                 ! Temporary variable to add up AOD
! VLIDORT working variables
!------------------------------
  integer                               :: ch                                       ! i-channel  
  integer                               :: iband                                    ! i-surfaceband
  real,allocatable                      :: pmom(:,:,:,:,:)                          ! elements of scattering phase matrix for vector calculations

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
    write(*,*) 'Channels [nm]: ',channels
    write(*,*) 'Cloud Fraction <= ', cldmax
    write(*,*) 'SZA < ', szamax
    write(*,*) 'VZA < ', vzamax
    if (scalar) write(*,*) 'Scalar calculations'
    if (.not. scalar) write(*,*) 'Vector calculations'
    if (trim(layout) /= '111') write(*,*) 'layout: ',trim(layout)
    write(*,*) ' '
  end if 

! Query for domain dimensions and missing value
!----------------------------------------------
  call mp_readDim("ew", MET_file, im)
  call mp_readDim("ns", MET_file, jm)
  call mp_readDim("lev", MET_file, km)
  call mp_readDim("time", MET_file,tm)
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
  ROT_           = dble(MISSING)
  PE_            = dble(MISSING)
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

    call getEdgeVars ( km, nobs, reshape(AIRDENS(c,:),(/km,nobs/)), &
                       reshape(DELP(c,:),(/km,nobs/)), ptop, &
                       pe, ze, te )   

    PE_(c,:) = pe(:,nobs)

    write(msg,'(A,I)') 'getEdgeVars ', myid
    call write_verbose(msg)
    
    call calc_qm()

    write(msg,'(A,I)') 'calc_qm ', myid
    call write_verbose(msg)  

!   Surface Reflectance Parameters
!   ----------------------------------
    albedo = 0

!   Rayleigh Optical Thickness
!   --------------------------
    call VLIDORT_ROT_CALC (km, nch, nobs, dble(channels), dble(pe), dble(ze), dble(te), &
                                   dble(MISSING),verbose, &
                                   ROT, depol, ierr )  


!   Aerosol Optical Properties
!   --------------------------
    call VLIDORT_getAOPvector ( mieTables, km, nobs, nch, nq, channels, vnames, verbose, &
                          qm, reshape(RH(c,:),(/km,nobs/)),&
                          nMom,nPol, tau, ssa, g, pmom, ierr )

    if ( ierr /= 0 ) then
      print *, 'cannot get aerosol optical properties'
      call MPI_ABORT(MPI_COMM_WORLD,myid,ierr)      
    end if

    TAU_(c,:,:) = tau(:,:,nobs)
    SSA_(c,:,:) = ssa(:,:,nobs)
    G_(c,:,:)   = g(:,:,nobs)
    do ch=1,nch      
        ROT_(c,:,ch) = ROT(:,nobs,ch)
    end do

    write(msg,*) 'getAOP ', myid
    call write_verbose(msg)

!   Call VlIDORT
!   ------------

!   Simple lambertian surface model
!   -------------------------------
    ! Call to vlidort vector code
    call VLIDORT_Vector_Lambert (km, nch, nobs ,dble(channels), nstreams, plane_parallel, nMom,   &
               nPol, ROT, depol, dble(tau), dble(ssa), dble(pmom), dble(pe), dble(ze), dble(te), albedo,&
               (/dble(SZA(c))/), &
               (/dble(RAA(c))/), &
               (/dble(VZA(c))/), &
               dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, Q, U, ierr)
   

!   Check VLIDORT Status, Store Outputs in Shared Arrays
!   ----------------------------------------------------    
    call mp_check_vlidort(radiance_VL_int,reflectance_VL_int)  
    radiance_VL(c,:)    = radiance_VL_int(nobs,:)
    reflectance_VL(c,:) = reflectance_VL_int(nobs,:)

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


      AOD = 0
      do k=1,km 
        AOD = AOD + reshape(TAU_(:,k,ch),(/clrm_total/))
      end do

      where(reshape(reflectance_VL(:,ch),(/clrm_total/)) .eq. MISSING)  AOD = MISSING
      

      call check(nf90_inq_varid(ncid, 'aod_' // trim(adjustl(msg)), varid), "get aod vaird")
      call check(nf90_put_var(ncid, varid, unpack(AOD,clmask,field), &
                    start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out aod")
    end do
    call check( nf90_close(ncid), "close outfile" )

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


  write(INV_file,'(4A)')  trim(indir),'/LevelG/invariant/',trim(instname),'.lg1.invariant.nc4'

! OUTFILES
  write(OUT_file,'(4A)') trim(outdir),'/',trim(instname),'-g5nr.lc2.vlidort_blacksurface.'

  call outfile_extname(OUT_file)

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
      do i = 1, clrm
        if (RAA(i) < 0) then
          RAA(i) = RAA(i) + 360.0
        end if
      end do
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

    call MAPL_AllocNodeArray(SZA,(/clrm/),rc=ierr)
    call MAPL_AllocNodeArray(VZA,(/clrm/),rc=ierr)
    call MAPL_AllocNodeArray(RAA,(/clrm/),rc=ierr)

    call MAPL_AllocNodeArray(TAU_,(/clrm,km,nch/),rc=ierr)
    call MAPL_AllocNodeArray(SSA_,(/clrm,km,nch/),rc=ierr)
    call MAPL_AllocNodeArray(G_,(/clrm,km,nch/),rc=ierr)
    call MAPL_AllocNodeArray(ROT_,(/clrm,km,nch/),rc=ierr)
    call MAPL_AllocNodeArray(PE_,(/clrm,km+1/),rc=ierr)
    
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

    allocate (ROT(km,nobs,nch))
    allocate (depol(nch))
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
    integer,dimension(nch)             :: qVarID, uVarID     
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

    write(comment,'(A)') 'VLIDORT simulation run from geo_vlidort.x'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'history',trim(comment)),"history attr")

    call check(nf90_put_att(ncid,NF90_GLOBAL,'grid_inputs',trim(INV_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'angle_inputs',trim(ANG_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'land_inputs',trim(LAND_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'met_inputs',trim(MET_file)),"input files attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'aerosol_inputs',trim(AER_file)),"input files attr")


    call check(nf90_put_att(ncid,NF90_GLOBAL,'inputs',trim(comment)),"input files attr")
    write(comment,'(A)') 'n/a'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'references',trim(comment)),"references attr") 

    write(comment,'(A)') 'This file contains VLIDORT ' // &
                         'top of the atmosphere ' // &
                         'reflectance over a black surface from GEOS-5 parameters sampled ' // &
                         ' on the '//lower_to_upper(trim(instname))//' geostationary grid '
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
    call ESMF_ConfigGetAttribute(cf, scalar, label = 'SCALAR:',default=.TRUE.)
    call ESMF_ConfigGetAttribute(cf, plane_parallel, label = 'PLANE_PARALLEL:',default=.FALSE.)
    call ESMF_ConfigGetAttribute(cf, nstreams, label = 'NSTREAMS:',default=6)
    call ESMF_ConfigGetAttribute(cf, szamax, label = 'SZAMAX:',default=80.0)
    call ESMF_ConfigGetAttribute(cf, vzamax, label = 'VZAMAX:',default=80.0)
    call ESMF_ConfigGetAttribute(cf, cldmax, label = 'CLDMAX:',default=0.01)
    call ESMF_ConfigGetAttribute(cf, nodemax, label = 'NODEMAX:',default=1) 
    call ESMF_ConfigGetAttribute(cf, version, label = 'VERSION:',default='1.0') 
    call ESMF_ConfigGetAttribute(cf, layout, label = 'LAYOUT:',default='111')    


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
    call MAPL_DeallocNodeArray(SZA,rc=ierr) 
    call MAPL_DeallocNodeArray(VZA,rc=ierr) 
    call MAPL_DeallocNodeArray(RAA,rc=ierr) 

    call MAPL_DeallocNodeArray(TAU_,rc=ierr)
    call MAPL_DeallocNodeArray(SSA_,rc=ierr)
    call MAPL_DeallocNodeArray(G_,rc=ierr)
    call MAPL_DeallocNodeArray(ROT_,rc=ierr)
    call MAPL_DeallocNodeArray(PE_,rc=ierr)
    if (.not. scalar) then
      call MAPL_DeallocNodeArray(Q_,rc=ierr)
      call MAPL_DeallocNodeArray(U_,rc=ierr)
    end if

    call MAPL_DeallocNodeArray(radiance_VL,rc=ierr) 
    call MAPL_DeallocNodeArray(reflectance_VL,rc=ierr) 

  end subroutine deallocate_shared


end program geo_vlidort
