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
   integer :: ierr
   integer :: myid, npet, CoresPerNode
   integer :: im, jm, lm
   integer :: ncid, varid
   integer :: p
   integer :: startl, countl, endl
   integer, pointer :: nlayer(:) => null()
   real    :: memusage

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

!  Read the data
!  ------------------------------
   if (myid == 0) then  
      call sys_tracker()    
      write(*,*) 'Reading the cloud fraction on PE', myid
      call check( nf90_open(MET_file,NF90_NOWRITE,ncid), "opening MET file")
      call check( nf90_inq_varid(ncid,"CLDTOT",varid), "getting CLDTOT varid")
      call check( nf90_get_var(ncid,varid,CLDTOT, start = (/ 1, 1, 1 /), count=(/im,jm,1/)), "reading CLDTOT")
      call check( nf90_close(ncid), "closing MET file")
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
  call par_layreader("AIRDENS", AER_file, AIRDENS)
  call par_layreader("DELP", AER_file, DELP)

   ! do p = 0, npet-1
   !    if (myid == p) then
   !       startl = myid*nlayer(p+1)+1
   !       countl = nlayer(p+1)
   !       endl   = startl + countl
   !       write(*,*)'Reading ', countl, ' layers starting on ', startl, ' on PE ', myid
   !       call check( nf90_open(AER_file,NF90_NOWRITE,ncid), "opening AER file")
   !       call check( nf90_inq_varid(ncid,"AIRDENS",varid), "getting AIRDENS varid")
   !       call check( nf90_get_var(ncid,varid,AIRDENS(:,:,startl:endl), start = (/ 1, 1, startl, 1 /), count=(/im,jm,countl,1/)), "reading CLDTOT")
   !       call check( nf90_close(ncid), "closing MET file")
   !    endif
   ! end do

!  Wait for everyone to finish reading
!   call MPI_Barrier(MPI_COMM_WORLD, ierr)
!   if (myid == 0) then  
!      call sys_tracker()   
!   endif



!  Although read on individual PEs, U,V,T should have the same
!  data in all PEs. Let's verify that.
!  -----------------------------------------------------------
   ! call MPI_Barrier(MPI_COMM_WORLD,ierr)    
   ! if ( myid == 0 ) then
   !    print *, '--- Array Statistics ---'
   ! end if
   ! write(*,*)'CLDTOT: ',myid,maxval(CLDTOT),minval(CLDTOT)
   ! call MPI_Barrier(MPI_COMM_WORLD,ierr)
   ! write(*,*)'AIRDENS: ',myid,maxval(AIRDENS),minval(AIRDENS)
   ! call MPI_Barrier(MPI_COMM_WORLD,ierr)
   ! write(*,*)'DELP: ',myid,maxval(DELP),minval(DELP)
   ! call MPI_Barrier(MPI_COMM_WORLD,ierr)
   ! if ( myid == 0 ) then
   !    open (unit = 2, file="test_write")
   !    write(2,*) "these values should all have the same min max values"
   !    close(2)
   ! end if

!  Wait for everyone to finish
   call MAPL_SyncSharedMemory(rc=ierr)

!  All done
!  --------
   call shutdown()

   contains

      subroutine par_layreader(varname, filename, var)
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
      end subroutine par_layreader

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
