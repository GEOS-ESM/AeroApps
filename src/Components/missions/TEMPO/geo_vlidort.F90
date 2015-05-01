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

!  Test flat
!  -----------
   logical               :: test_shmem = .False.

!  File names
!  ----------
   character(len=256)    :: MET_file, AER_file 

!  Global, 3D arrays to be allocated using SHMEM
!  ---------------------------------------------
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

!  VLIDORT input arrays
!  ---------------------------
   real, pointer         :: pe(:,:) => null()      ! edge pressure [Pa]
   real, pointer         :: ze(:,:) => null()      ! edge height above sfc [m]
   real, pointer         :: te(:,:) => null()      ! edge Temperature [K]
   integer               :: nobs                   ! number of profiles VLIDORT will work on
   real, parameter       :: ptop = 1.0             ! top (edge) pressure [Pa]
   integer               :: nch                    ! number of channels
   integer, parameter    :: nq = 14                ! number of tracers
   integer               :: verbose 
   real, pointer         :: qm(:,:,:) => null()    ! (mixing ratio) * delp/g
!   real, pointer         :: rh(:,:) => null()      ! relative humidity
   real, pointer         :: tau(:,:,:) => null()   ! aerosol optical depth
   real, pointer         :: ssa(:,:,:) => null()   ! single scattering albedo
   real, pointer         :: g(:,:,:) => null()     ! asymmetry factor
   real, pointer         :: albedo(:,:) => null()  ! surface albedo
   real, pointer         :: solar_zenith(:) => null() 
   real, pointer         :: relat_azymuth(:) => null()
   real, pointer         :: sensor_zenith(:) => null()
   real, parameter       :: channels(1) = 550.0    ! channels to simulate
   real, parameter       :: MISSING = -999
   real*8,pointer        :: radiance_VL(:,:) => null()      ! TOA normalized radiance from VLIDORT
   real*8,pointer        :: reflectance_VL(:,:) => null()  ! TOA reflectance from VLIDORT

   character(len=*), parameter      :: rcfile = 'Aod_EOS.rc'  ! resource file
   character(len=16), parameter     :: vnames_string(nq) = (/'du001', 'du002', 'du003', 'du004', 'du005', &
                                                             'ss001', 'ss002', 'ss003', 'ss004', 'ss005', &
                                                             'BCphobic', 'BCphilic',                      &
                                                             'OCphobic', 'OCphilic'/) ! array of variable name strings
   character             :: vnames(nq,16)          ! character array of variable names

!   real, allocatable     :: ps(:,:)                  ! saturation vapor pressure temporary variable
!   real, allocatable     :: rh(:,:)                  ! relative humidity temporary variable


!  Miscellaneous
!  -------------
   integer               :: ierr                      ! MPI error message
   integer               :: status(MPI_STATUS_SIZE)   ! MPI status
   integer               :: myid, npet, CoresPerNode  ! MPI dimensions and processor id
   integer               :: im, jm, km                ! size of TEMPO domain
   integer               :: i, j, k, n                ! size of TEMPO domain
   integer               :: ncid, varid, rcid         ! netcdf ids
   integer               :: p                         ! i-processor
   integer               :: startl, countl, endl      ! array indices and counts for layers to be read
   integer, allocatable  :: nlayer(:)                 ! how many layers each processor reads
   character(len=61)     :: msg                       ! message printed by shmen_test
   real, parameter       :: grav = 9.81                  ! gravity


!  Initialize MPI
!  --------------
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,npet,ierr)
   if (myid == 0) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'

!  Initialize SHMEM
!  ----------------
   CoresPerNode = MAPL_CoresPerNodeGet(MPI_COMM_WORLD,rc=ierr) ! a must
   call MAPL_InitializeShmem(rc=ierr)

!  For now hard code file names and dimensions, could query file for dimensions
!  ----------------------------------------------------------------------------
   im = 1250
   jm = 2000
   km = 72
   MET_file = "/nobackup/TEMPO/met_Nv/Y2005/M12/tempo-g5nr.lb2.met_Nv.20051231_00z.nc4"
   AER_file = "/nobackup/TEMPO/aer_Nv/Y2005/M12/tempo-g5nr.lb2.aer_Nv.20051231_00z.nc4"

 
!  Allocate the Global arraya using SHMEM
!  It will be available on all processors
!  ---------------------------------------------------------
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

!  Allocate arrays that will be copied on each processor - unshared
!  ---------------------------------------------------------

!  Needed for vlidort
!  -----------------------
   verbose = 1
   nobs    = 1
   nch     = 1
   allocate (pe(km+1,nobs))
   allocate (ze(km+1,nobs))
   allocate (te(km+1,nobs))
   allocate (qm(km,nq,nobs))
!   allocate (rh(km,nobs))
   allocate (tau(km,nch,nobs))
   allocate (ssa(km,nch,nobs))
   allocate (g(km,nch,nobs))
   allocate (albedo(nch,nobs))

   allocate (solar_zenith(nobs))
   allocate (relat_azymuth(nobs))
   allocate (sensor_zenith(nobs))

   allocate (radiance_VL(nobs,nch))
   allocate (reflectance_VL(nobs, nch))

!   allocate (rh(km,nobs))
!   allocate (ps(km,nobs))

!  Needed for reading
!  ----------------------
   allocate (nlayer(npet))

!  Read the cloud data
!  ------------------------------
   if (myid == 0) then  
      write(*,*) 'Allocated all shared memory'
      call sys_tracker()    
       call readvar2D("CLDTOT", MET_file, CLDTOT)
      write(*,*) 'Read cloud information'
      call sys_tracker()
    end if

!  Wait for everyone to have access to CLDTOT and shared memory
!  -----------------------------------------------------------   
   call MAPL_SyncSharedMemory(rc=ierr)

!  Everyone Figure out how many layers each PE has to read
!  -----------------------------
   if (npet >= km) then
      nlayer(1:npet) = 1
   else if (npet < km) then
      nlayer(1:npet) = km/npet
      nlayer(npet)   = nlayer(npet) + mod(km,npet)
   end if 

!  Read the netcdf variables layer-by-layer in parallel
!  ------------------------------
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

!  Wait for everyone to finish and print max memory used
!  -----------------------------------------------------------  
   call MAPL_SyncSharedMemory(rc=ierr)
   if (myid == 0) then 
      write(*,*) 'Read all variables' 
      call sys_tracker()   
   end if     


   if (test_shmem) then
!  Although read on individual PEs, all shared variables should have the same
!  data in all PEs. Let's verify that.
!  -----------------------------------------------------------   
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

!  Wait for everyone to finish and print max memory used
!  -----------------------------------------------------------  
      call MAPL_SyncSharedMemory(rc=ierr)
      if (myid == 0) then  
         write(*,*) 'Tested shared memory' 
         call sys_tracker()   
      end if   
    end if 

!  Prepare inputs for VLIDORT
!  split up domain among processors
!  ------------------------------
   if (myid == 0) then  
      call getEdgeVars ( km, nobs, reshape(AIRDENS(1,1,:),(/km,nobs/)), reshape(DELP(1,1,:),(/km,nobs/)), ptop, &
                          pe, ze, te )

!  this is a temporary hack just to be able to test vlidort
      ! do n = 1, nobs
      !    do k = 1, km
      !       ps(k,n) = 6.11 * 10 ** ((7.5*(te(k,n)-273.15))/(237.3+(te(k,n)-273.15)))   ! hPa
      !       ps(k,n) = 100.*ps(k,n)  !Pa
      !       rh(k,n) = 100.*QV(1,1,k)*pe(k,n)/ps(k,n)  ! %
      !    end do
      ! end do

      call strarr_2_chararr(vnames_string,nq,16,vnames)

      do n = 1, nobs
         call calc_qm(DU001,1,n)
         call calc_qm(DU002,2,n)
         call calc_qm(DU003,3,n)
         call calc_qm(DU004,4,n)
         call calc_qm(DU005,5,n)
         call calc_qm(SS001,6,n)
         call calc_qm(SS002,7,n)
         call calc_qm(SS003,8,n)
         call calc_qm(SS004,9,n)
         call calc_qm(SS005,10,n)
         call calc_qm(BCPHOBIC,11,n)
         call calc_qm(BCPHILIC,12,n)
         call calc_qm(OCPHOBIC,13,n)
         call calc_qm(OCPHILIC,14,n)
      end do
 
      call getAOPscalar ( km, nobs, nch, nq, rcfile, channels, vnames, verbose, &
                           qm, reshape(RH(1,1,:),(/km,nobs/)), &
                           tau, ssa, g, ierr )

      albedo(:,:) = 0.05
      sensor_zenith(:) = 2.0
      solar_zenith(:)  = 60.0
      relat_azymuth    = 20.0


       call Scalar (km, nch, nobs ,dble(channels),        &
                    dble(tau), dble(ssa), dble(g), dble(pe), dble(ze), dble(te), dble(albedo),&
                    dble(solar_zenith), dble(relat_azymuth), dble(sensor_zenith), &
                    dble(MISSING),verbose,radiance_VL,reflectance_VL, ierr)

       write(*,*) reflectance_VL
       write(*,*) 'Called VLIDORT Scalar' 
       call sys_tracker()
   end if 

!  All done
!  --------
   call MAPL_SyncSharedMemory(rc=ierr)
   call shutdown()

!-----------------------------------------------------------------------------------
   contains

   subroutine calc_qm(var,n,iobs)
      real, intent(in), dimension(:,:,:)     :: var
      integer, intent(in)                    :: n, iobs

      integer                                :: k

      do k = 1, km
            qm(k,n,iobs) = var(1,1,k)*DELP(1,1,k)/grav
      end do
   end subroutine calc_qm


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

         integer                       :: ncid, varid


!         write(*,'(A,A,A,I4)')'Reading ',trim(varname), ' on PE ', myid
         call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
         call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
         call check( nf90_get_var(ncid,varid,var), "reading " // varname)
         call check( nf90_close(ncid), "closing " // filename)
      end subroutine readvar2D

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


         do p = 0, npet-1
            if (myid == p) then
               startl = myid*nlayer(p+1)+1
               countl = nlayer(p+1)
               endl   = startl + countl
!               write(*,'(A,A,A,I2,A,I2,A,I4)')'Reading ',trim(varname),' ', countl, ' layers starting on ', startl, ' on PE ', myid
               call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
               call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
               call check( nf90_get_var(ncid,varid,var(:,:,startl:endl), start = (/ 1, 1, startl, 1 /), count=(/im,jm,countl,1/)), "reading " // varname)
               call check( nf90_close(ncid), "closing MET file")
            end if
         end do
      end subroutine mp_layreadvar

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
         real,dimension(:,:,:)      :: var

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
!         deallocate (rh)
         deallocate (tau)
         deallocate (ssa)
         deallocate (g)

         ! deallocate (rh)
         ! deallocate(ps)

         deallocate(nlayer)


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
