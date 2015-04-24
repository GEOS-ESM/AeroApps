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
   real, pointer :: AIRDENS(:,:,:,:) => null()
!   real, pointer :: T(:,:,:) => null()

!  Miscellaneous
!  -------------
   integer :: ierr
   integer :: myid, npet, CoresPerNode
   integer :: im, jm, lm
   integer :: ncid, varid
   real    :: memusage

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
   lm = 72
   MET_file = "/nobackup/TEMPO/met_Nv/Y2005/M12/tempo-g5nr.lb2.met_Nv.20051231_00z.nc4"
   AER_file = "/nobackup/TEMPO/aer_nv/Y2005/M12/tempo-g5nr.lb2.aer_Nv.20051231_00z.nc4"
   !T_file = "/home/adasilva/opendap/c1440_NR/DATA/0.0625_deg/inst/inst30mn_3d_T_Nv/Y2005/M07/D12/c1440_NR.inst30mn_3d_T_Nv.20050712_1200z.nc4"

!  Allocate the Global CLDTOT array using SHMEM
!  It will be available on all processors
!  ---------------------------------------------------------
   call MAPL_AllocNodeArray(CLDTOT,(/im,jm/),rc=ierr)
   call MAPL_AllocNodeArray(AIRDENS,(/im,jm,lm/),rc=ierr)
!   call MAPL_AllocNodeArray(T,(/im,jm,lm/),rc=ierr)

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
!  Read the data files in parallel
!  ------------------------------
   ! if (npet >= 3) then
   !    if (myid == 0) then
   !       write(*,*)'Reading U on PE ', myid
   !       call check( nf90_open(U_file,NF90_NOWRITE,ncid), "opening U file")
   !       call check( nf90_inq_varid(ncid,"U",varid), "getting U varid")
   !       call check( nf90_get_var(ncid,varid,U), "reading U")
   !       call check( nf90_close(ncid), "closing U file")
   !    else if (myid == 1) then
   !       write(*,*)'Reading V on PE ', myid
   !       call check( nf90_open(V_file,NF90_NOWRITE,ncid), "opening V file")
   !       call check( nf90_inq_varid(ncid,"V",varid), "getting V varid")
   !       call check( nf90_get_var(ncid,varid,V), "reading V")
   !       call check( nf90_close(ncid), "closing V file")
   !    else if (myid == 2) then
   !       write(*,*)'Reading T on PE ', myid
   !       call check( nf90_open(T_file,NF90_NOWRITE,ncid), "opening T file")
   !       call check( nf90_inq_varid(ncid,"T",varid), "getting T varid")
   !       call check( nf90_get_var(ncid,varid,T), "reading T")
   !       call check( nf90_close(ncid), "closing T file")
   !    end if
   ! else
   !    call shutdown() ! should never happen
   ! end if

!  Although read on individual PEs, U,V,T should have the same
!  data in all PEs. Let's verify that.
!  -----------------------------------------------------------
   ! call MPI_Barrier(MPI_COMM_WORLD,ierr)    
   ! if ( myid == 0 ) then
   !      print *
   !      print *, '--- Array Statistics ---'
   ! end if
   ! write(*,*)'U: ',myid,maxval(U),minval(U)
   ! call MPI_Barrier(MPI_COMM_WORLD,ierr)
   ! write(*,*)'V: ',myid,maxval(V),minval(V)
   ! call MPI_Barrier(MPI_COMM_WORLD,ierr)
   ! write(*,*)'T: ',myid,maxval(T),minval(T)
   ! call MPI_Barrier(MPI_COMM_WORLD,ierr)
   ! if (myid == 0) write(*,*)'all processors should have same min/max values!'

!  All done
!  --------
   call shutdown()

   contains

      subroutine check(status, loc)

         integer, intent(in) :: status
         character(len=*), intent(in) :: loc

         if(status /= NF90_NOERR) then
            write (*,*) "Error at ", loc
            write (*,*) NF90_STRERROR(status)
         end if

      end subroutine check

      subroutine shutdown()

         ! shmem must deallocate shared memory arrays
         call MAPL_DeallocNodeArray(CLDTOT,rc=ierr)
         !call MAPL_DeallocNodeArray(V,rc=ierr)
         !call MAPL_DeallocNodeArray(T,rc=ierr)

         call MAPL_FinalizeShmem (rc=ierr)

         call MPI_Finalize(ierr)

      end subroutine shutdown

end program shmem_reader
