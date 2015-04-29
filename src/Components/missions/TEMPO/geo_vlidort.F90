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
   character, pointer    :: vnames(:,:) => null() ! variable name
   real, pointer         :: qm(:,:,:) => null()    ! (mixing ratio) * delp/g
   real, pointer         :: rh(:,:) => null()      ! relative humidity
   real, pointer         :: tau(:,:,:) => null()   ! aerosol optical depth
   real, pointer         :: ssa(:,:,:) => null()   ! single scattering albedo
   real, pointer         :: g(:,:,:) => null()     ! asymmetry factor
   real, pointer         :: channels(:) => null()  ! channels to simulate
   character(len=*), parameter      :: rcfile = 'Aod_EOS.rc'  ! resource file


!  Miscellaneous
!  -------------
   integer               :: ierr                      ! MPI error message
   integer               :: status(MPI_STATUS_SIZE)   ! MPI status
   integer               :: myid, npet, CoresPerNode  ! MPI dimensions and processor id
   integer               :: im, jm, km                ! size of TEMPO domain
   integer               :: ncid, varid, rcid         ! netcdf ids
   integer               :: p                         ! i-processor
   integer               :: startl, countl, endl      ! array indices and counts for layers to be read
   integer, allocatable  :: nlayer(:)                 ! how many layers each processor reads
   character(len=61)     :: msg                       ! message printed by shmen_test


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
   AER_file = "/nobackup/TEMPO/aer_nv/Y2005/M12/tempo-g5nr.lb2.aer_Nv.20051231_00z.nc4"

!  Allocate the Global arraya using SHMEM
!  It will be available on all processors
!  ---------------------------------------------------------
   call MAPL_AllocNodeArray(CLDTOT,(/im,jm/),rc=ierr)
   call MAPL_AllocNodeArray(AIRDENS,(/im,jm,km/),rc=ierr)
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
   allocate (vnames(nq,16))
   allocate (pe(km+1,nobs))
   allocate (ze(km+1,nobs))
   allocate (te(km+1,nobs))
   allocate (qm(km,nq,nobs))
   allocate (rh(km,nobs))
   allocate (tau(km,nch,nobs))
   allocate (ssa(km,nch,nobs))
   allocate (g(km,nch,nobs))
   allocate (channels(nch))

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
  call par_layreadvar("AIRDENS", AER_file, AIRDENS)
  call par_layreadvar("DELP", AER_file, DELP)
  call par_layreadvar("DU001", AER_file, DU001)
  call par_layreadvar("DU002", AER_file, DU002)
  call par_layreadvar("DU003", AER_file, DU003)
  call par_layreadvar("DU004", AER_file, DU004)
  call par_layreadvar("DU005", AER_file, DU005)
  call par_layreadvar("SS001", AER_file, SS001)
  call par_layreadvar("SS002", AER_file, SS002)
  call par_layreadvar("SS003", AER_file, SS003)
  call par_layreadvar("SS004", AER_file, SS004)
  call par_layreadvar("SS005", AER_file, SS005)
  call par_layreadvar("BCPHOBIC", AER_file, BCPHOBIC)
  call par_layreadvar("BCPHILIC", AER_file, BCPHILIC)
  call par_layreadvar("OCPHOBIC", AER_file, OCPHOBIC)
  call par_layreadvar("OCPHILIC", AER_file, OCPHILIC)

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
!  ------------------------------
   if (myid == 0) then  
      call getEdgeVars ( km, nobs, reshape(airdens(1,1,:),(/km,nobs/)), reshape(delp(1,1,:),(/km,nobs/)), ptop, &
                           pe, ze, te )

      ! call getAOPscalar ( km, nobs, nch, nq, rcfile, channels, vnames, verbose, &
      !                      qm, rh, &
      !                      tau, ssa, g, ierr )
   end if 

!  All done
!  --------
   call shutdown()

   contains
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
         character(len=*), intent(in)  ::  varname
         character(len=*), intent(in)  ::  filename
         real, dimension(:,:)        ::  var

         integer                       :: ncid, varid


!         write(*,'(A,A,A,I4)')'Reading ',trim(varname), ' on PE ', myid
         call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
         call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
         call check( nf90_get_var(ncid,varid,var), "reading " // varname)
         call check( nf90_close(ncid), "closing " // filename)
      end subroutine readvar2D

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    par_layreadvar
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
      subroutine par_layreadvar(varname, filename, var)
         character(len=*), intent(in)  ::  varname
         character(len=*), intent(in)  ::  filename
         real, dimension(im,jm,km)        ::  var

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
      end subroutine par_layreadvar

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
         deallocate (vnames)
         deallocate (qm)
         deallocate (rh)
         deallocate (tau)
         deallocate (ssa)
         deallocate (g)

         ! shmem must deallocate shared memory arrays
         call MAPL_DeallocNodeArray(CLDTOT,rc=ierr)
         call MAPL_DeallocNodeArray(AIRDENS,rc=ierr)
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

         call MAPL_FinalizeShmem (rc=ierr)

         call MPI_Finalize(ierr)

      end subroutine shutdown


end program shmem_reader
