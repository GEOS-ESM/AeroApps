!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SparseComm
!
! !DESCRIPTION:
!
! !INTERFACE:

Module m_SparseComm

#ifdef USE_SHMEM_PTR
  Use m_SparseComm_SHMEM_PTR, only : SparseComm
  Use m_SparseComm_SHMEM_PTR, only : Analyze
  Use m_SparseComm_SHMEM_PTR, only : PartialReduceScatter
  Use m_SparseComm_SHMEM_PTR, only : PartialGather
  Use m_SparseComm_SHMEM_PTR, only : DescribeComm
  Use m_SparseComm_SHMEM_PTR, only : Clean
  Use m_SparseComm_SHMEM_PTR, only : ReleaseSymmetricMemory
  Use m_SparseComm_SHMEM_PTR, only : GATHER
  Use m_SparseComm_SHMEM_PTR, only : REDUCE_SCATTER
#endif

#ifdef USE_SHMEM_CRAY
  Use m_SparseComm_SHMEM_CRAY, only : SparseComm
  Use m_SparseComm_SHMEM_CRAY, only : Analyze
  Use m_SparseComm_SHMEM_CRAY, only : PartialReduceScatter
  Use m_SparseComm_SHMEM_CRAY, only : PartialGather
  Use m_SparseComm_SHMEM_CRAY, only : DescribeComm
  Use m_SparseComm_SHMEM_CRAY, only : Clean
  Use m_SparseComm_SHMEM_CRAY, only : ReleaseSymmetricMemory
  Use m_SparseComm_SHMEM_CRAY, only : GATHER
  Use m_SparseComm_SHMEM_CRAY, only : REDUCE_SCATTER
#endif

#ifdef USE_MPI_1
  Use m_SparseComm_MPI_1, only : SparseComm
  Use m_SparseComm_MPI_1, only : Analyze
  Use m_SparseComm_MPI_1, only : PartialReduceScatter
  Use m_SparseComm_MPI_1, only : PartialGather
  Use m_SparseComm_MPI_1, only : DescribeComm
  Use m_SparseComm_MPI_1, only : Clean
  Use m_SparseComm_MPI_1, only : ReleaseSymmetricMemory
  Use m_SparseComm_MPI_1, only : GATHER
  Use m_SparseComm_MPI_1, only : REDUCE_SCATTER
#endif

#ifdef USE_MPI_2
  Use m_SparseComm_MPI_2, only : SparseComm
  Use m_SparseComm_MPI_2, only : Analyze
  Use m_SparseComm_MPI_2, only : PartialReduceScatter
  Use m_SparseComm_MPI_2, only : PartialGather
  Use m_SparseComm_MPI_2, only : DescribeComm
  Use m_SparseComm_MPI_2, only : Clean
  Use m_SparseComm_MPI_2, only : ReleaseSymmetricMemory
  Use m_SparseComm_MPI_2, only : GATHER
  Use m_SparseComm_MPI_2, only : REDUCE_SCATTER
#endif

  Implicit None
  Private

  Public :: SparseComm              ! public datatype
  Public :: Analyze                 ! analyze communication pattern
  Public :: PartialReduceScatter    ! implements collective operation
  Public :: PartialGather           ! implements collective operation
  Public :: DescribeComm            ! writes a description of communication pattern
  Public :: Clean                   ! deallocates SparseComm objects
  Public :: ReleaseSymmetricMemory  ! release symmetric heap memory for complete clean-up
  Public :: GATHER
  Public :: REDUCE_SCATTER

! !REVISION HISTORY:
!       11Feb02 - Tom Clune <clune@sgi.com>
!               . Split module into separate module for each paradigm
!       29Oct01 - Tom Clune
!               . Finally added this intro section

!EOP ___________________________________________________________________

End Module m_SparseComm
