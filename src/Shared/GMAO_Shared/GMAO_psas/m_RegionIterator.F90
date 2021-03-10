!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_RegionIterator
!
! !DESCRIPTION:
!
!   This module implements an iterator class for use with the Spherical_Partitions class,
!   defined below.
!  
! !INTERFACE:
!

#include "assert.H"

! This module implements at iterator class for the regions in
! a partition.  In particular, there is a need to loop over regions
! at a given refinement level.  In C++, these two classes would be
! friends, but with F90 I am simply making more things public than 
! otherwise ought to be.

Module m_RegionIterator
  Use m_die,only : assert_
  Implicit None
  Private

  Public :: RegionIterator
  Public :: New
  Public :: Clean
  Public :: Reset
  Public :: Next
  Public :: base
  Public :: IsLeaf
  Public :: GetNumberOfSubregions
  Public :: GetLevel
  Public :: Clone
  Public :: IsNested

  Public :: rotate
  Public :: reflect

  Public :: INVALID, DONE
  Public :: base_refinement, max_refinement
  Public :: n_refine
  Public :: coef
  Public :: RefinementToRegionIndex

  Integer, Parameter :: n_refine = 4 ! subregions per region
  Integer, Parameter :: INVALID = -9999
  Integer, Parameter :: DONE    = -1
  Integer, Parameter :: base_refinement = -1
  Integer, Parameter :: max_refinement = 7

#ifndef F95
  Type RegionIterator
     Integer, Pointer :: refinement_path(:)
     Integer :: index
  End Type RegionIterator
#else
  Type RegionIterator
     Integer, Pointer :: refinement_path(:) => Null()
     Integer :: index = INVALID
  End Type RegionIterator
#endif
  Integer, Parameter :: coef(base_refinement:max_refinement) = &
       & (/ 1, 5, 21, 85, 341, 1365, 5461, 21845, 87381 /)

  Integer :: err ! used for error conditions

! !REVISION HISTORY:
!        1Dec00 - Tom CLune and Peter Lyster <lys@dao.gsfc.nasa.gov>
!                 Install Clune's new algorithms in support of
!                 advanced geometry and Icosahedral refinement scheme etc.
!                 The initial (birthday) files are in root molotov1.gsfc.nasa.gov:/home/jguo/Gcvs/ in 
!                 repository mppi under tag new_refinement
!EOP ___________________________________________________________________

Contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 
  ! Given an integer list specifying the refinement path from the
  ! icosahedron, return a unique index for that path.
  ! Refinement is specified with integers 0..n_refine, where 0 means
  ! no refinement and 1..n_refine correspond to the usual triangular 
  ! refinement.
  !
  ! The ordering chosen for index is such that all descendants are 
  ! contiguous.  This simplifies resorting binned data if a given region is
  ! subsequently refined.
  !
  ! Note that no assumption about the base partition is made here.
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Function RefinementToRegionIndex(iter) Result (idx)
    Implicit None
    Type (RegionIterator), Intent(In) :: iter
    Integer :: idx

    Integer :: i, j, n, idx_

    ASSERT_NOMSG(Associated(iter%refinement_path))

    n = Ubound(iter%refinement_path,1)

    idx_ = 0

    Do i = base_refinement, n
       j = iter%refinement_path(i)
       If (j == 0) Exit
       idx_ = idx_ + 1 + coef(n-i-1)*(j-1)
    End Do

    ASSERT(idx_ >= 0)

    idx = idx_

  End Function RefinementToRegionIndex

  Subroutine new(iter, level)
    Implicit None
    Type (RegionIterator), Intent(Out) :: iter
    Integer, Intent(In) :: level

!_JG    ASSERT_NOMSG(.not.Associated(iter%refinement_path))
!	At this point, new_ter%refinement_path is undefined.
!	Checking if it is associated is violating Fortran 90
!	standard.

    Allocate(iter%refinement_path(base_refinement:level), STAT = err)
    ALWAYS_ASSERT(err == 0)
    
    iter%refinement_path(base_refinement) = 0
    iter%index = DONE

  End Subroutine new

  Subroutine clean(iter)
    Type (RegionIterator), Intent(InOut) :: iter

    ALWAYS_ASSERT_NOMSG(Associated(iter%refinement_path))

    Deallocate(iter%refinement_path, STAT = err)
    ASSERT(err == 0)
    Nullify(iter%refinement_path)

  End Subroutine clean

  Integer Function reset(iter, level, index)
    Implicit None
    Type (RegionIterator), Intent(InOut) :: iter
    Integer, Intent(In), Optional :: level
    Integer, Intent(In), Optional :: index
    Integer :: level_
    Integer :: idx, i, n

    ASSERT_NOMSG(Associated(iter%refinement_path))
    ASSERT_NOMSG(.not. (Present(index) .and. Present(level)))

    level_ = base_refinement
    If (Present(level)) level_ = level

    n = Ubound(iter%refinement_path,1)
    ASSERT(level_ >= -1)
    ASSERT(level_ <= n)

    iter%refinement_path = 0

    If (Present(index)) Then
       idx = index
       Do i = base_refinement, Ubound(iter%refinement_path,1)
	  idx = idx - 1
	  iter%refinement_path(i) = 1 + idx / coef(n - 1 - i)
	  idx = idx - (iter%refinement_path(i)-1) * coef(n - 1 - i)
	  If (idx == 0) Exit
       End Do
    Else
       iter%refinement_path(base_refinement) = 1
       iter%refinement_path(base_refinement+1:level_) = 1
       iter%refinement_path(level_+1:) = 0
    End If

    iter%index = RefinementToRegionIndex(iter)
#ifndef NDEBUG
    If (Present(index)) Then
       If (iter%index /= index) Then
	  print*,index, iter%index
	  print*,iter%refinement_path
       End If
       ASSERT(iter%index == index)
    End If
#endif

    reset = iter%index

  End Function reset
    
  Function next(iter, level) Result (idx)
    Use m_Icosahedron, only : n_faces
    Implicit None
    Type (RegionIterator), Intent(InOut) :: iter
    Integer, Optional, Intent(In) :: level
    Integer :: idx

    Logical :: done_

    Integer :: level_, n, i

    ASSERT_NOMSG(Associated(iter%refinement_path))
    n = Ubound(iter%refinement_path,1)

    If (Present(level)) Then
       level_ = level
    Else
       level_ = n
       Do i =  base_refinement, n
	  If (iter%refinement_path(i) == 0) Then
	     level_ = i
	     Exit
	  End If
       End Do
    End If

    ASSERT_NOMSG(level_ >= base_refinement)
    ASSERT(level_ <= n)

    done_ = .false.

    iter%refinement_path(level_+1:) = 0
    iter%refinement_path(:level_-1) = max(1,iter%refinement_path(:level_-1))

    Do 
       iter%refinement_path(level_) = iter%refinement_path(level_) + 1

       Select Case (level_)
       Case (base_refinement+1:)
	  If (iter%refinement_path(level_) <= n_refine) Exit
       Case (base_refinement)
	  If (iter%refinement_path(level_) == n_faces + 1) iter%index = DONE
	  Exit
       End Select

       If (Present(level)) Then
	  iter%refinement_path(level_) = 1 
       Else
	  iter%refinement_path(level_) = 0
       End If
       level_ = level_ - 1

    End Do

    If (iter%index /= DONE) Then
       iter%index = RefinementToRegionIndex(iter)
    End If

    ASSERT_NOMSG(All(iter%refinement_path(base_refinement+1:) <= n_refine))
    idx = iter%index

  End Function next
  
  Integer Function base(iter)
    Implicit None
    Type (RegionIterator), Intent(In) :: iter

    ASSERT_NOMSG(associated(iter%refinement_path))
    base = iter%refinement_path(base_refinement)

  End Function base

  Logical Function IsLeaf(iter)
    Implicit None
    Type (RegionIterator), Intent(In) :: iter

    If (iter%refinement_path(Ubound(iter%refinement_path,1)) > 0)  Then
       IsLeaf = .true.
    Else
       IsLeaf = .false.
    End If

  End Function IsLeaf

  Function GetNumberOfSubregions(iter) Result(n)
    Implicit None
    Type (RegionIterator), Intent(In) :: iter
    Integer :: n
    Integer :: i, ub

    ASSERT_NOMSG(Associated(iter%refinement_path))
    ub = ubound(iter%refinement_path,1)
    n = 0
    Do i = base_refinement, ub
       If (iter%refinement_path(i) == 0) Then
	  n = coef(ub-i) - 1
	  exit
       End If
    End Do
    
  End Function GetNumberOfSubregions

  Function GetLevel(iter) Result(n)
    Implicit None
    Type (RegionIterator), Intent(In) :: iter
    Integer :: n
    Integer :: i, ub

    ub = ubound(iter%refinement_path,1)
    n = ub
    Do i = base_refinement + 1, ub
       If (iter%refinement_path(i) == 0) Then
	  n = i - 1
	  Exit
       End If
    End Do
    
  End Function GetLevel

  Subroutine Clone(old_iter, new_iter)
    Implicit None
    Type (RegionIterator), Intent(In) :: old_iter
    Type (RegionIterator), Intent(Out) :: new_iter

!_JG    ASSERT_NOMSG(.not.associated(new_iter%refinement_path))
!	At this point, new_ter%refinement_path is undefined.
!	Checking if it is associated is violating Fortran 90
!	standard.
    Call new(new_iter, level = ubound(old_iter%refinement_path,1))
    new_iter%refinement_path = old_iter%refinement_path
    new_iter%index = old_iter%index

  End Subroutine Clone

  Logical Function IsNested(iter_1, iter_2)
    Implicit None
    Type (RegionIterator), Intent(In) :: iter_1, iter_2
    Integer :: i
    Logical :: flag

    flag = .true.
    Do i = base_refinement, ubound(iter_1%refinement_path,1)
       If ((iter_1%refinement_path(i) == 0) .or. (iter_2%refinement_path(i) == 0)) Exit
       If (iter_1%refinement_path(i) /= iter_2%refinement_path(i)) Then
	  flag = .false.
	  Exit
       End If
    End Do
    
    IsNested = flag
	  
  End Function IsNested

  Function Rotate(iter) Result (idx)
    Type (RegionIterator), Intent(InOut) :: iter
    Integer :: idx, i, j

    ASSERT_NOMSG(associated(iter%refinement_path))

    Do i = base_refinement + 1, ubound(iter%refinement_path,1)
       j = iter%refinement_path(i)
       Select Case(j)
       Case(1:3)
	  iter%refinement_path(i) = 1 + mod(j,3)
	  exit
       Case(0)
	  exit
       Case (4)
	  ! continue
       End Select

    End Do
    iter%index = RefinementToRegionIndex(iter)

    idx = iter%index

  End Function Rotate
    
  Function Reflect(iter) Result (idx)
    Type (RegionIterator), Intent(InOut) :: iter
    Integer :: idx, i, j

    Do i = base_refinement + 1, ubound(iter%refinement_path,1)
       j = iter%refinement_path(i)
       Select Case(j)
       Case(2:3)
	  iter%refinement_path(i) = 5-j ! swap2 and 3
       Case(0)
	  exit
       End Select
    End Do
    
    iter%index = RefinementToRegionIndex(iter)

    idx = iter%index

  End Function Reflect
    
End Module m_RegionIterator
