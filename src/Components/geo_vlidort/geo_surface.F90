!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    geo_surface
! PURPOSE
!     Reads in parallel the model data (in a netcdf file) interpolated to TEMPO grid 
!     The variables are needed as input to vlidort
!     A shared memory array created with a call to MAPL_ShmemMod is used for the variables
!     Do some filtering and run the vlidort code 
!     **--> Only the vlidort surface supplement is run and outputs are stored.
! INPUT
!      Resource File
! OUTPUT
!     None
!  HISTORY
!     27 April 2015 P. Castellanos adapted from A. da Silva shmem_reader.F90
!     July 2015 P. Castellanos adapted rrom geo_vlidort.F90
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#  include "MAPL_Generic.h"
#  include "MAPL_ErrLogMain.h"
program geo_vlidort_surface

  use ESMF                         ! ESMF modules
  use MAPL_Mod
  use MAPL_ShmemMod                ! The SHMEM infrastructure
  use SURFACE                      ! modules for calculating surface reflectance 
  use netcdf                       ! for reading the NR files
  use netcdf_helper                ! Module with netcdf routines
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
  character(len=7)                      :: surfdate
  character(len=2)                      :: time 
  character(len=256)                    :: instname, indir, outdir, surfname
  character(len=256)                    :: surfband               ! flag to use nearest-neighbor interpolation or an exact value given
  integer                               :: surfbandm              ! number of wavelength bands or channels in surface reflectance data file
  integer, allocatable                  :: surfband_i(:)          ! surface band indeces that overlap with vlidort channels
  real, allocatable                     :: surfband_c(:)          ! modis band center wavelength
  real, allocatable                     :: channels(:)            ! channels to simulate
  integer                               :: nch                    ! number of channels  
  real                                  :: cldmax                 ! Cloud Filtering  
  real                                  :: szamax                 ! Geomtry filtering

! Test flag
! -----------
  logical                               :: test_shmem = .False.

! File names
! ----------
  character(len=256)                    :: MET_file, ANG_file, INV_file, SURF_file, LAND_file, OUT_file 

! Global, 3D inputs to be allocated using SHMEM
! ---------------------------------------------
  real, pointer                         :: KISO(:,:) => null()
  real, pointer                         :: KVOL(:,:) => null()
  real, pointer                         :: KGEO(:,:) => null()
  real, pointer                         :: LER(:,:) => null()
  real, pointer                         :: SZA(:) => null()
  real, pointer                         :: VZA(:) => null()
  real, pointer                         :: RAA(:) => null()  


! VLIDORT output arrays
!-------------------------------
!                                  Intermediate Unshared Arrays
!                                  -----------------------------
  real*8, allocatable                   :: albedo(:,:)                            ! bi-directional surface reflectance
  real*8, allocatable                   :: k_vol(:,:)                             ! Ross-thick kernel
  real*8, allocatable                   :: k_geo(:,:)                             ! Li-sparse kernel  

!                                  Final Shared Arrays
!                                  -------------------
  real*8, pointer                       :: Kvol_(:,:) => null()                   ! Ross-thick kernel
  real*8, pointer                       :: Kgeo_(:,:) => null()                   ! Li-sparse kernel
  real*8, pointer                       :: ALBEDO_(:,:) => null()                 ! bi-directional surface reflectance

  real*8,allocatable                    :: field(:,:)                             ! Template for unpacking shared arrays

! VLIDORT working variables
!------------------------------
  integer                               :: ch                                       ! i-channel  
  integer                               :: iband                                    ! i-surfaceband

! MODIS Kernel variables
!--------------------------
  real*8, allocatable                   :: kernel_wt(:,:,:)                         ! kernel weights (/fiso,fgeo,fvol/)
  real*8, allocatable                   :: param(:,:,:)                             ! Li-Sparse parameters 
                                                                                    ! param1 = crown relative height (h/b)
                                                                                    ! param2 = shape parameter (b/r)
  real                                  :: surf_missing                                                                 


! Satellite domain variables
!------------------------------
  integer                               :: im, jm, km, tm                            ! size of TEMPO domain
  integer                               :: i, j, k, n                                ! TEMPO domain working variable
  integer                               :: starti, counti, endi                      ! array indices and counts for each processor
  integer, allocatable                  :: nclr(:)                                   ! how many clear pixels each processor works on
  integer                               :: clrm                                      ! number of clear pixels
  integer                               :: c                                         ! clear pixel working variable
  real, allocatable                     :: CLDTOT(:,:)                               ! GEOS-5 cloud fraction
  real, allocatable                     :: FRLAND(:,:)                               ! GEOS-5 land fraction
  real, allocatable                     :: SOLAR_ZENITH(:,:)                         ! solar zenith angles used for data filtering
  logical, allocatable                  :: clmask(:,:)                               ! cloud-land mask

! netcdf variables
!----------------------  
  integer                               :: ncid                                           ! netcdf file id
  integer,allocatable,dimension(:)      :: volVarID, geoVarID, albVarID                   ! netcdf OUT_file variable IDs

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
  call get_config(arg)

! Write out settings to use
! --------------------------
  if (myid == 0) then
    write(*,*) 'Simulating ', lower_to_upper(trim(instname)),' domain on ',date,' ', time, 'Z'
    write(*,*) 'Input directory: ',trim(indir)
    write(*,*) 'Output directory: ',trim(outdir)
    write(*,*) 'BRDF dataset: ',trim(surfname),' ',trim(surfdate)
    write(*,*) 'Channels [nm]: ',channels
    if (lower_to_upper(surfband) == 'EXACT') write(*,*) 'Using exact surface relflectance parameters on bands : ',surfband_i
    if (lower_to_upper(surfband) == 'INTERPOLATE') write(*,*) 'Using interpolated surface reflectance parameters' 
    write(*,*) 'Cloud Fraction <= ', cldmax
    write(*,*) 'SZA < ', szamax
    write(*,*) ' '
  end if 


! Query for domain dimensions and missing value
!----------------------------------------------
  call mp_readDim("ew", MET_file, im)
  call mp_readDim("ns", MET_file, jm)
  call mp_readDim("lev", MET_file, km)
  call mp_readDim("time", MET_file,tm)
  if (lower_to_upper(surfname) == 'MAIACRTLS') then
    call mp_readVattr("missing_value", SURF_file, "Kiso", surf_missing) 
  else
    call mp_readVattr("missing_value", SURF_file, "SRFLER354", surf_missing) 
  end if
  call mp_readVattr("missing_value", MET_FILE, "CLDTOT", g5nr_missing)

! Allocate arrays that will be copied on each processor - unshared
! -----------------------------------------------------------------
  call allocate_unshared()

! Create OUTFILE
! --------------
  if ( MAPL_am_I_root() )  call create_outfile(date, time, volVarID, geoVarID, albVarID)
  call MAPL_SyncSharedMemory(rc=ierr)
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
 call read_surf()
 call read_angles()

! Wait for everyone to finish reading and print max memory used
! ******This needs to loop through all processors....
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

! Initialize outputs to be safe
! -------------------------------
  Kvol_           = dble(MISSING)
  Kgeo_           = dble(MISSING)
  ALBEDO_         = dble(MISSING)
  call MAPL_SyncSharedMemory(rc=ierr)
! Prepare inputs and run VLIDORT
! -----------------------------------
  
  if (myid == 0) then
    starti = 1
  else
    starti = sum(nclr(1:myid))+1
  end if
  counti = nclr(myid+1)
  endi   = starti + counti - 1

  do c = starti, endi

!   Surface Reflectance Parameters
!   ----------------------------------
    if ( (lower_to_upper(surfname) == 'MAIACRTLS') ) then
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


!   Call VlIDORT
!   ------------
    if ( (lower_to_upper(surfname) /= 'MAIACRTLS') ) then
      ! call VLIDORT_Scalar_Lambert (km, nch, nobs ,dble(channels),        &
      !         dble(tau), dble(ssa), dble(g), dble(pe), dble(ze), dble(te), albedo,&
      !         (/dble(SZA(c))/), &
      !         (/dble(RAA(c))/), &
      !         (/dble(VZA(c))/), &
      !         dble(MISSING),verbose,radiance_VL_int,reflectance_VL_int, ROT, ierr)
      write(*,*) 'other surfaces not implemented yet'
    else if ( ANY(kernel_wt == surf_missing) ) then
!     Save code for pixels that were not gap filled
!     ---------------------------------------------    
      k_vol    = -500
      k_geo    = -500
      albedo   = -500
      ierr     = 0
    else             
!     MODIS BRDF Surface Model
!     ------------------------------
      call SURFACE_LandMODIS (km, nch, nobs, dble(channels),        &
              kernel_wt, param, &
              (/dble(SZA(c))/), &
              (/dble(RAA(c))/), &
              (/dble(VZA(c))/), &
              dble(MISSING),verbose, albedo, k_vol, k_geo, ierr )  

    end if          
    
!   Check VLIDORT Status, Store Outputs in Shared Arrays
!   ----------------------------------------------------    
    call mp_check_vlidort(k_vol, k_geo, albedo)  
    Kvol_(c,:)    = k_vol(nobs,:)
    Kgeo_(c,:)    = k_geo(nobs,:)
    ALBEDO_(c,:)  = albedo(nobs,:)
    
    write(msg,*) 'VLIDORT Calculations DONE', myid, ierr
    call write_verbose(msg)

!   Keep track of progress of each processor
!   -----------------------------------------        
    if (nint(100.*real(c-starti)/real(counti)) > progress) then
      progress = nint(100.*real(c-starti)/real(counti))
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

!                             Write to main OUT_File
!                             ----------------------
    call check( nf90_open(OUT_file, nf90_write, ncid), "opening file " // OUT_file )

    do ch = 1, nch
      call check(nf90_put_var(ncid, volVarID(ch), unpack(reshape(Kvol_(:,ch),(/clrm/)),clmask,field), &
                  start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out kvol")
      call check(nf90_put_var(ncid, geoVarID(ch), unpack(reshape(Kgeo_(:,ch),(/clrm/)),clmask,field), &
                  start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out kgeo")
      call check(nf90_put_var(ncid, albVarID(ch), unpack(reshape(ALBEDO_(:,ch),(/clrm/)),clmask,field), &
                    start = (/1,1,1,nobs/), count = (/im,jm,1,nobs/)), "writing out albedo")
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
!  call MAPL_SyncSharedMemory(rc=ierr)
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

  below = minloc(abs(xint - x), dim = 1, mask = (xint - x) .LE. 0)
  above = minloc(abs(xint - x), dim = 1, mask = (xint - x) .GT. 0)

  
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
  write(MET_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/',trim(instname),'-g5nr.lb2.met_Nv.',date,'_',time,'z.nc4'
  write(ANG_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/',trim(instname),'.lb2.angles.',date,'_',time,'z.nc4'
  write(INV_file,'(4A)')  trim(indir),'/LevelG/invariant/',trim(instname),'.lg1.invariant.nc4'
  if ( lower_to_upper(surfname) == 'MAIACRTLS' ) then
    write(SURF_file,'(6A)') trim(indir),'/BRDF/',trim(surfname),'.',surfdate,'.hdf'
  else
    write(SURF_file,'(4A)') trim(indir),'/SurfLER/tempo-omi.SurfLER.',date(5:6),'.nc4'
  end if

  write(LAND_file,'(4A)') trim(indir),'/LevelB/invariant/',trim(instname),'-g5nr.lb2.asm_Nx.nc4' 

! OUTFILES
  write(OUT_file,'(4A)') trim(outdir),'/',trim(instname),'-g5nr.lb2.surface.'

  call outfile_extname(OUT_file)

end subroutine filenames

subroutine outfile_extname(file)
  character(len=256),intent(inout)     :: file
  character(len=256)                   :: chmax, chmin
  integer                              :: i

  
  if (lower_to_upper(surfname) == 'MAIACRTLS' .and. lower_to_upper(surfband) == 'INTERPOLATE') then
    write(file,'(2A)') trim(file),'iMAIACRTLS.'
  else if (lower_to_upper(surfname) == 'MAIACRTLS' .and. lower_to_upper(surfband) == 'EXACT') then
    write(file,'(2A)') trim(file),'MAIACRTLS.'
  else if (lower_to_upper(surfname) /= 'MAIACRTLS' .and. lower_to_upper(surfband) == 'INTERPOLATE') then
    write(file,'(2A)') trim(file),'ilambertian.'
  else if (lower_to_upper(surfname) /= 'MAIACRTLS' .and. lower_to_upper(surfband) == 'EXACT') then
    write(file,'(2A)') trim(file),'lambertian.'
  end if

  write(chmax,'(F10.2)') maxval(channels)
  write(chmin,'(F10.2)') minval(channels)

  i = index(chmax,'.')
  chmax(i:i) = 'd'
  i = index(chmin,'.')
  chmin(i:i) = 'd'

  if (nch == 1) then    
    write(file,'(7A)') trim(file),date,'_',time,'z_',trim(adjustl(chmin)),'nm.nc4'
  else 
    write(file,'(9A)') trim(file),date,'_',time,'z_',trim(adjustl(chmin)),'-',trim(adjustl(chmax)),'nm.nc4'
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
  subroutine mp_check_vlidort(k_vol, k_geo, albedo)
    real*8, dimension(:,:)    :: k_vol, k_geo, albedo


    if (ANY(k_vol == MISSING) .or. ANY(k_geo == MISSING) .or. ANY(albedo == MISSING)) then
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

    if (MAPL_am_I_root()) then
      if (lower_to_upper(surfname) == 'MAIACRTLS') then
 
        call readvar3D("Kiso", SURF_file, temp)
        call reduceProfile(temp,clmask,KISO)

        call readvar3D("Kvol", SURF_file, temp)
        call reduceProfile(temp,clmask,KVOL)

        call readvar3D("Kgeo", SURF_file, temp)
        call reduceProfile(temp,clmask,KGEO)

        write(*,*) '<> Read BRDF data to shared memory' 
      else
        call read_LER(temp)        
        call reduceProfile(temp,clmask,LER)
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
    real, allocatable          :: SAA(:), VAA(:)
    integer                    :: clrm

    if (myid == 0) then
      allocate(SAA(clrm))
      allocate(VAA(clrm))

      call readvar2D("solar_zenith", ANG_file, temp)
      SZA = pack(temp,clmask)

      call readvar2D("sensor_zenith", ANG_file, temp)
      VZA = pack(temp,clmask)

      call readvar2D("solar_azimuth", ANG_file, temp)
      SAA = pack(temp,clmask)

      call readvar2D("sensor_azimuth", ANG_file, temp)
      VAA = pack(temp,clmask)

      ! define according to photon travel direction
      SAA = SAA + 180.0
      do i = 1, clrm
        if (SAA(i) >= 360.0) then
            SAA(i) = SAA(i) - 360.0
        end if
      end do

      RAA = VAA - SAA
      do i = 1, clrm
        if (RAA(i) < 0) then
          RAA(i) = RAA(i) + 360.0
        end if
      end do
      
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

    if (lower_to_upper(surfname) == 'MAIACRTLS') then
      call MAPL_AllocNodeArray(KISO,(/clrm,surfbandm/),rc=ierr)
      call MAPL_AllocNodeArray(KVOL,(/clrm,surfbandm/),rc=ierr)
      call MAPL_AllocNodeArray(KGEO,(/clrm,surfbandm/),rc=ierr)
    else
      call MAPL_AllocNodeArray(LER,(/clrm,surfbandm/),rc=ierr)
    end if 
    call MAPL_AllocNodeArray(SZA,(/clrm/),rc=ierr)
    call MAPL_AllocNodeArray(VZA,(/clrm/),rc=ierr)
    call MAPL_AllocNodeArray(RAA,(/clrm/),rc=ierr)

    call MAPL_AllocNodeArray(Kvol_,(/clrm,nch/),rc=ierr)
    call MAPL_AllocNodeArray(Kgeo_,(/clrm,nch/),rc=ierr)
    call MAPL_AllocNodeArray(ALBEDO_,(/clrm,nch/),rc=ierr)
    
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
    allocate (albedo(nobs,nch))
    allocate (k_vol(nobs,nch))
    allocate (k_geo(nobs, nch))    

    if (lower_to_upper(surfname) == 'MAIACRTLS') then
      allocate (kernel_wt(nkernel,nch,nobs))
    end if
    allocate (param(nparam,nch,nobs))

  ! Needed for reading
  ! ----------------------
    allocate (nclr(npet)) 
    nclr = 0

    allocate(clmask(im,jm))
    clmask = .False.

  ! Netcdf IDs
  ! ----------
    allocate (volVarID(nch))
    allocate (geoVarID(nch))
    allocate (albVarID(nch))

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

  subroutine create_outfile(date, time, volVarID, geoVarID, albVarID)
    character(len=*) ,intent(in)       :: date, time
    integer,dimension(:),intent(out)   :: volVarID, geoVarID, albVarID  !OUT_file variables    
    
    integer                            :: ncid
    integer                            :: timeDimID, ewDimID, nsDimID, levDimID, chaDimID       
    integer                            :: scantimeVarID, clonVarID, clatVarID
    integer                            :: timeVarID, levVarID, ewVarID, nsVarID
    integer                            :: ch

    real,allocatable,dimension(:,:)    :: clon, clat, sza, vza, raa
    real,allocatable,dimension(:)      :: scantime, ew, ns, tyme, lev

    character(len=2000)                :: comment

!                                 OUT_FILE 
!                                ----------
    ! Open File
    call check(nf90_create(OUT_file, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // OUT_file)

    ! Create dimensions
    call check(nf90_def_dim(ncid, "time", tm, timeDimID), "creating time dimension")
    call check(nf90_def_dim(ncid, "lev", 1, levDimID), "creating ns dimension") !km
    call check(nf90_def_dim(ncid, "ew", im, ewDimID), "creating ew dimension") !im
    call check(nf90_def_dim(ncid, "ns", jm, nsDimID), "creating ns dimension") !jm

    ! Global Attributes
    write(comment,'(A)') 'Surface Reflectance Simulation of GEOS-5 '//lower_to_upper(trim(instname))//' Sampler'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'title',trim(comment)),"title attr")

    write(comment,'(A)') 'NASA/Goddard Space Flight Center'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'institution',trim(comment)),"institution attr")

    write(comment,'(A)') 'Global Model and Assimilation Office'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'source',trim(comment)),"source attr")

    write(comment,'(A)') 'Simulation run from geo_surface.x'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'history',trim(comment)),"history attr")

    write(comment,'(A)') trim(INV_file)  // CHAR(13) // &
                         trim(ANG_file)  // CHAR(13) //  &
                         trim(LAND_file) // CHAR(13) // &
                         trim(MET_file)  // CHAR(13) //  &
                         trim(SURF_file) 
    call check(nf90_put_att(ncid,NF90_GLOBAL,'inputs',trim(comment)),"input files attr")
    write(comment,'(A)') 'n/a'
    call check(nf90_put_att(ncid,NF90_GLOBAL,'references',trim(comment)),"references attr") 

    write(comment,'(A)') 'This file contains simulated surface reflectance (albedo) ' // &
                         ' on the '//lower_to_upper(trim(instname))//' geostationary grid '
    call check(nf90_put_att(ncid,NF90_GLOBAL,'comment',trim(comment)),"comment attr")   

    if (lower_to_upper(surfname) == 'MAIACRTLS' .and. lower_to_upper(surfband) == 'INTERPOLATE') then
      write(comment,'(A)') 'MAIAC RTLS surface BRDF kernel weights interpolated to channel'
    else if (lower_to_upper(surfname) == 'MAIACRTLS' .and. lower_to_upper(surfband) == 'EXACT') then
      write(comment,'(A)') 'MAIAC RTLS surface BRDF kernel weights without inerpolation to channel'
    else if (lower_to_upper(surfname) /= 'MAIACRTLS' .and. lower_to_upper(surfband) == 'INTERPOLATE') then
      write(comment,'(A)') 'Lambertian surface reflectance interpolated to channel'
    else if (lower_to_upper(surfname) /= 'MAIACRTLS' .and. lower_to_upper(surfband) == 'EXACT') then
      write(comment,'(A)') 'Lambertian surface reflectance without interpolation to channel'
    end if
    call check(nf90_put_att(ncid,NF90_GLOBAL,'surface_comment',trim(comment)),"surface_comment")

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
      call check(nf90_def_var(ncid, 'kvol_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),volVarID(ch)),"create kvol var")      
      call check(nf90_def_var(ncid, 'kgeo_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),geoVarID(ch)),"create kgeo var")
      call check(nf90_def_var(ncid, 'albedo_' // trim(adjustl(comment)) ,nf90_float,(/ewDimID,nsDimID,levDimID,timeDimID/),albVarID(ch)),"create albedo var")
    end do

    ! Variable Attributes
!                                            Surface Reflectance
!                                          ------------------------  
    do ch=1,size(channels)
      write(comment,'(F10.2,A)') channels(ch), ' nm Kvol'
      call check(nf90_put_att(ncid,volVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm Ross-Thick Kernel'
      call check(nf90_put_att(ncid,volVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,volVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,volVarID(ch),'units','None'),"units attr")
      call check(nf90_put_att(ncid,volVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

      write(comment,'(F10.2,A)') channels(ch), ' nm Kgeo'
      call check(nf90_put_att(ncid,geoVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm Li-sparse Kernel'
      call check(nf90_put_att(ncid,geoVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,geoVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,geoVarID(ch),'units','None'),"units attr")
      call check(nf90_put_att(ncid,geoVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")

      write(comment,'(F10.2,A)') channels(ch), ' nm Albedo'
      call check(nf90_put_att(ncid,albVarID(ch),'standard_name',trim(adjustl(comment))),"standard_name attr")
      write(comment,'(F10.2,A)') channels(ch), ' nm Bi-Directional Surface Reflectance'
      call check(nf90_put_att(ncid,albVarID(ch),'long_name',trim(adjustl(comment))),"long_name attr")
      call check(nf90_put_att(ncid,albVarID(ch),'missing_value',real(MISSING)),"missing_value attr")
      call check(nf90_put_att(ncid,albVarID(ch),'units','None'),"units attr")
      call check(nf90_put_att(ncid,albVarID(ch),"_FillValue",real(MISSING)),"_Fillvalue attr")
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

    call readvar2D("clon", INV_file, clon)
    call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

    call readvar2D("clat", INV_file, clat)
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

    call shmem_test2D('KIO',KISO)
    call shmem_test2D('KVOL',KVOL)
    call shmem_test2D('KGEO',KGEO)

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
    call ESMF_ConfigGetAttribute(cf, surfname, label = 'SURFNAME:',default='MAIACRTLS')
    call ESMF_ConfigGetAttribute(cf, surfdate, label = 'SURFDATE:',__RC__)
    call ESMF_ConfigGetAttribute(cf, szamax, label = 'SZAMAX:',default=90.0)
    call ESMF_ConfigGetAttribute(cf, cldmax, label = 'CLDMAX:',default=0.01)
    call ESMF_ConfigGetAttribute(cf, surfband, label = 'SURFBAND:', default='INTERPOLATE')
    call ESMF_ConfigGetAttribute(cf, surfbandm, label = 'SURFBANDM:',__RC__)

    ! Check that LER configuration is correct
    !----------------------------------------
    if (lower_to_upper(surfname) /= 'MAIACRTLS' ) then
      if (surfbandm /= 2) then
        surfbandm = 2
        if (MAPL_am_I_root()) then
          write(*,*) 'Wrong number of surface bands.  If not MAIACRTLS, SURFBANDM must equal 2'
          write(*,*) 'Forcing surfbandm = 2'          
        end if
      end if
    end if

    ! Figure out number of channels and read into vector
    !------------------------------------------------------
    nch =  ESMF_ConfigGetLen(cf, label = 'CHANNELS:',__RC__)
    allocate (channels(nch))
    call ESMF_ConfigGetAttribute(cf, channels, label = 'CHANNELS:', default=550.)

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
    if (lower_to_upper(surfname) == 'MAIACRTLS') then
      call MAPL_DeallocNodeArray(KISO,rc=ierr) 
      call MAPL_DeallocNodeArray(KVOL,rc=ierr) 
      call MAPL_DeallocNodeArray(KGEO,rc=ierr) 
    else
      call MAPL_DeallocNodeArray(LER,rc=ierr) 
    end if
    call MAPL_DeallocNodeArray(SZA,rc=ierr) 
    call MAPL_DeallocNodeArray(VZA,rc=ierr) 
    call MAPL_DeallocNodeArray(RAA,rc=ierr) 

    call MAPL_DeallocNodeArray(Kvol_,rc=ierr)
    call MAPL_DeallocNodeArray(Kgeo_,rc=ierr)
    call MAPL_DeallocNodeArray(ALBEDO_,rc=ierr)

  end subroutine deallocate_shared


end program geo_vlidort_surface
