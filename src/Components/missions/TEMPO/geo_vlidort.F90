! program to read several variables from a netcdf file into a
! a shared memory array using MAPL_ShmemMod

program shmem_reader

   use MAPL_ShmemMod    ! The SHMEM infrastructure
   use netcdf           ! for reading the NR files

   implicit none
   include "mpif.h"

!  File names
!  ----------
   character(len=256) :: MET_file, AER_file 

!  Global, 3D arrays to be allocated using SHMEM
!  ---------------------------------------------
   real, pointer :: CLDTOT(:,:) => null()
   real, pointer :: AIRDENS(:,:,:) => null()
   real, pointer :: DELP(:,:,:) => null()

!  Miscellaneous
!  -------------
   integer               :: ierr
   integer               :: myid, npet, CoresPerNode
   integer               :: im, jm, lm
   integer               :: ncid, varid, rcid
   integer               :: p
   integer               :: startl, countl, endl
   integer, pointer      :: nlayer(:) => null()
   real                  :: memusage
   character(len=61)    :: msg
   integer               :: status(MPI_STATUS_SIZE)

!  Initialize MPI
!  --------------
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,npet,ierr)
   allocate (nlayer(npet))
   if (myid == 0) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'

!  Initialize SHMEM
!  ----------------
   CoresPerNode = MAPL_CoresPerNodeGet(MPI_COMM_WORLD,rc=ierr) ! a must
   call MAPL_InitializeShmem(rc=ierr)

!  For now hard code file names and dimensions, could query file for dimensions
!  ----------------------------------------------------------------------------
   im = 1250
   jm = 2000
   lm = 72
   MET_file = "/nobackup/TEMPO/met_Nv/Y2005/M12/tempo-g5nr.lb2.met_Nv.20051231_00z.nc4"
   AER_file = "/nobackup/TEMPO/aer_nv/Y2005/M12/tempo-g5nr.lb2.aer_Nv.20051231_00z.nc4"

!  Allocate the Global arraya using SHMEM
!  It will be available on all processors
!  ---------------------------------------------------------
   call MAPL_AllocNodeArray(CLDTOT,(/im,jm/),rc=ierr)
   call MAPL_AllocNodeArray(AIRDENS,(/im,jm,lm/),rc=ierr)
   call MAPL_AllocNodeArray(DELP,(/im,jm,lm/),rc=ierr)

!  Read the cloud data
!  ------------------------------
   if (myid == 0) then  
      call sys_tracker()    
      call readvar2D("CLDTOT", MET_file, CLDTOT)
      call sys_tracker()
    end if

!  Wait for everyone to have access to CLDTOT and shared memory
   call MAPL_SyncSharedMemory(rc=ierr)

!  Everyone Figure out how many layers each PE has to read
!  -----------------------------
   if (npet >= lm) then
      nlayer(1:npet) = 1
   else if (npet < lm) then
      nlayer(1:npet) = lm/npet
      nlayer(npet)   = nlayer(npet) + mod(lm,npet)
   end if 

!  Read the AIRDENS variable layer-by-layer in parallel
!  ------------------------------
  call par_layreadvar("AIRDENS", AER_file, AIRDENS)
  call par_layreadvar("DELP", AER_file, DELP)


!  Although read on individual PEs, all shared variables should have the same
!  data in all PEs. Let's verify that.
!  -----------------------------------------------------------   
   !call MPI_FILE_OPEN(MPI_COMM_WORLD, 'test_write', MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
   !                    MPI_INFO_NULL, rcid, ierr)


   if ( myid  == 0 ) then
         open (unit = 2, file="test_write")
         write(2,'(A)') '--- Array Statistics ---'
   end if

   call shmem_test3D('AIRDENS',AIRDENS)

!  Wait for everyone to finish
   call MAPL_SyncSharedMemory(rc=ierr)
   if (myid == 0) then  
      call sys_tracker()   
   endif   

!  All done
!  --------
   call shutdown()

   contains

      subroutine readvar2D(varname, filename, var)
         character(len=*), intent(in)  ::  varname
         character(len=*), intent(in)  ::  filename
         real, dimension(im,jm)        ::  var

         integer                       :: ncid, varid


         write(*,'(A,A,A,I4)')'Reading ',trim(varname), ' on PE ', myid
         call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
         call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
         call check( nf90_get_var(ncid,varid,var), "reading " // varname)
         call check( nf90_close(ncid), "closing " // filename)
      end subroutine readvar2D

      subroutine par_layreadvar(varname, filename, var)
         character(len=*), intent(in)  ::  varname
         character(len=*), intent(in)  ::  filename
         real, dimension(im,jm,lm)        ::  var

         integer                       :: p, startl, countl, endl
         integer                       :: ncid, varid


         do p = 0, npet-1
            if (myid == p) then
               startl = myid*nlayer(p+1)+1
               countl = nlayer(p+1)
               endl   = startl + countl
               write(*,'(A,A,A,I2,A,I2,A,I4)')'Reading ',trim(varname),' ', countl, ' layers starting on ', startl, ' on PE ', myid
               call check( nf90_open(filename,NF90_NOWRITE,ncid), "opening file " // filename)
               call check( nf90_inq_varid(ncid,varname,varid), "getting varid for " // varname)
               call check( nf90_get_var(ncid,varid,var(:,:,startl:endl), start = (/ 1, 1, startl, 1 /), count=(/im,jm,countl,1/)), "reading " // varname)
               call check( nf90_close(ncid), "closing MET file")
            end if
         end do
      end subroutine par_layreadvar

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    shmem_test
! PURPOSE
!     test that shared memory arrays have same values across all processors
! INPUT
!     varname: string of variable name
!     var    : the variable to be checked
! OUTPUT
!     Writes to the file test_write the min and max value of the variable as 
!     reported by each processor
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      subroutine shmem_test3D(varname,var)
         character(len=*), intent(in)  :: varname
         real,dimension(im,jm,lm)      :: var

         if ( myid  == 0 ) then
            open (unit = 2, file="test_write",position="append")
            do p = 1,npet-1
               call mpi_recv(msg, 61, MPI_CHARACTER, p, 1, MPI_COMM_WORLD, status, ierr)
               write(2,*) msg
            end do
         else
            write(msg,'(A9,I4,E24.17,E24.17)') varname,myid,maxval(var),minval(var)
            call mpi_send(msg, 61, MPI_CHARACTER, 0, 1, MPI_COMM_WORLD, ierr)
         end if
         if (myid == 0) then
            write(2,*) 'These should all have the same min/max values!'
            close(2)
         end if

      end subroutine shmem_test3D


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

         deallocate(nlayer)

         ! shmem must deallocate shared memory arrays
         call MAPL_DeallocNodeArray(CLDTOT,rc=ierr)
         call MAPL_DeallocNodeArray(AIRDENS,rc=ierr)
         call MAPL_DeallocNodeArray(DELP,rc=ierr)

         call MAPL_FinalizeShmem (rc=ierr)

         call MPI_Finalize(ierr)

      end subroutine shutdown

end program shmem_reader
