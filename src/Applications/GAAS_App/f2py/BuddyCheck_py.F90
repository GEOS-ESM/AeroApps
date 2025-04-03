!***********
! Fortran wrapper for the Python interface
! Use of openMP
!***********

subroutine Icosahedron_regions(nobs, xobs, yobs, zobs, kr)

use m_Spherical_Partition, only : NumberOfRegions, Initialize, Clean
use m_Spherical_Partition, only : Spherical_Partition, GetRegion
use m_Spherical_Partition, only : xyz2reg

   implicit none

! Arguments:
! ----------

integer, intent(in) :: nobs
real*8, intent(in) :: xobs(nobs)
real*8, intent(in) :: yobs(nobs)
real*8, intent(in) :: zobs(nobs)
integer, intent(out) :: kr(nobs)


! Region pointers
! ---------------
Type (Spherical_Partition) :: partition
integer, allocatable :: iregbeg(:)     ! Pointer to first ob in region
integer, allocatable :: ireglen(:)     ! No. of obs in region

integer:: n_levels, nregions

n_levels = 0

!Initialize Partitioner
!     ----------------------
call Initialize ( n_levels=0, partition=partition, compress=.true. )
nregions = NumberOfRegions(partition)  

! Compute index region from (x,y,z) after bin obs
! -------------------------
call XYZ2REG ( nobs, xobs, yobs, zobs, kr, partition )  ! bin obs into regions

end subroutine Icosahedron_regions
