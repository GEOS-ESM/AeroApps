#include "assert.H"
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Spherical_Partition
!
! !DESCRIPTION:
!
!   This module implements a class representing triangular partitions of the
!   unit sphere.  A strong assumption is made that the partition is based
!   upon refinements of the icosahedron accessed via m_Icosahedron, though
!   much of the structure is more general.  
!  
!   A few of the members are obsolete, but provided for backward compatability.
!   
! !INTERFACE:
!

Module m_Spherical_Partition
  Use m_die,only : assert_
  Use m_Icosahedron, Only : n_faces
  Use m_RegionIterator, Only : base_refinement, max_refinement, n_refine, INVALID
  Use m_Spherical_Triangle, only : Spherical_Triangle
  Implicit None
  Private

  Public :: Spherical_Partition
  Public :: Initialize
  Public :: Clean
  Public :: MAXREG_TRADITIONAL
  Public :: BuildMask
  Public :: xyz2reg
  Public :: LatLon_to_Region

  Public :: Get
  Public :: GetRegion
  Public :: NumberOfLeaves
  Public :: NumberOfRegions
  Public :: NumberOfNodes
  Public :: NumberOfEdges
  Public :: ListCenters
  Public :: ListAreas
!!$  Public :: ListNodes
  Public :: BaseRegion
  Public :: Base20Region
  Public :: Base80Region

  Public :: krname ! names of original 80 regions

  Interface Get
     Module Procedure Get_by_partition
     Module Procedure Get_by_index
  End Interface
!________________________

    public :: new	! create an object as a pointer
    public :: delete	! delete an object as a pointer

    interface new   ; module procedure new_   ; end interface
    interface delete; module procedure delete_; end interface

#ifndef F95
  ! f95 allows default components in derived types and poiner initialization.
  Type Region
     Type (Spherical_Triangle)     :: triangle
     Real, Dimension(3,n_refine-1) :: test
     Integer :: parent_index      
     Integer :: child_index(n_refine)
  End Type Region

  Type Spherical_Partition
     Private
     Integer                :: refinement_level
     Type (Region), Pointer :: region(:)       
     Integer                :: base_indices(n_faces)
     Logical                :: Compress
     Integer, Pointer       :: leaf(:)
  End Type Spherical_Partition
#else
  Type Region
     Type (Spherical_Triangle)     :: triangle
     Real, Dimension(3,n_refine-1) :: test
     Integer :: parent_index          = INVALID
     Integer :: child_index(n_refine) = INVALID
  End Type Region

  Type Spherical_Partition
     Private
     Integer                :: refinement_level = INVALID
     Type (Region), Pointer :: region(:)        => Null()
     Integer                :: base_indices(n_faces) = INVALID
     Logical                :: Compress = .false.
     Integer, Pointer       :: leaf(:) => Null()
  End Type Spherical_Partition
#endif

  ! Probably only one instance will ever be used, and keeping a
  ! private object here reduces the burden on the class user
  ! in most instances.  Methods will all have an optional dummy
  ! argument to override the default use of this private object.
  Type (Spherical_Partition), Target :: private_partition ! default

  Character(Len=*), Parameter :: myname='m_Spherical_Partition'

  Integer, Parameter :: MAXREG_TRADITIONAL = 80

  Integer :: err ! used for internal error handling

! !REVISION HISTORY:
!       27Mar02 - Tom Clune
!                 . Added headers to routines.
!                 . Added new interface 'BaseRegion'
!        1Dec00 - Tom CLune and Peter Lyster <lys@dao.gsfc.nasa.gov>
!                 Install Clune's new algorithms in support of
!                 advanced geometry and Icosahedral refinement scheme etc.
!                 The initial (birthday) files are in root molotov1.gsfc.nasa.gov:/home/jguo/Gcvs/ in 
!                 repository mppi under tag new_refinement
!EOP ___________________________________________________________________

!.......................................................................
!.... Icosahedral region names.

	Character(len=8), Parameter :: krname(80) =  (/ &
     &	  'npolkr01', 'eurokr02', 'natlkr03', 'eurokr04', 'npolkr05', &
     &	  'asiakr06', 'asiakr07', 'asiakr08', 'npolkr09', 'npackr10', &
     &	  'asiakr11', 'asiakr12', 'npolkr13', 'namrkr14', 'npackr15', &
     &	  'npackr16', 'npolkr17', 'natlkr18', 'namrkr19', 'namrkr20', &
     &	  'eatlkr21', 'eqafkr22', 'eqafkr23', 'eqafkr24', 'iocnkr25', &
     &	  'eqafkr26', 'seaskr27', 'iocnkr28', 'austkr29', 'seaskr30', &
     &	  'wpackr31', 'wpackr32', 'epackr33', 'wpackr34', 'epackr35', &
     &	  'epackr36', 'samrkr37', 'namrkr38', 'eatlkr39', 'samrkr40', &
     &	  'eqafkr41', 'iocnkr42', 'eqafkr43', 'eqafkr44', 'seaskr45', &
     &	  'austkr46', 'iocnkr47', 'iocnkr48', 'wpackr49', 'wpackr50', &
     &	  'austkr51', 'wpackr52', 'epackr53', 'epackr54', 'epackr55', &
     &	  'epackr56', 'eatlkr57', 'eatlkr58', 'samrkr59', 'samrkr60', &
     &	  'spolkr61', 'satlkr62', 'siockr63', 'siockr64', 'spolkr65', &
     &	  'siockr66', 'austkr67', 'siockr68', 'spolkr69', 'austkr70', &
     &	  'spackr71', 'spackr72', 'spolkr73', 'spackr74', 'spackr75', &
     &	  'spackr76', 'spolkr77', 'samrkr78', 'satlkr79', 'satlkr80' /)

Contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Initialize - Initialize a partition of the unit sphere
!
! !DESCRIPTION:
!
! This routine initializes a partition.  For convenience, 
! if the optional dummy argument is missing, the "private_partition" 
! is used instead.
!
! !INTERFACE:

  Subroutine Initialize(n_levels, partition, compress)
    Use m_Icosahedron, only : n_faces, init_Icosahedron => Initialize
    Use m_Icosahedron, only : GetFaceTriangles
    Use m_RegionIterator, only : NewIter => new, reset, next, RegionIterator, base, DONE
    Use m_Spherical_Triangle, only : Refine
    Use m_mall, only : mall_ci, mall_co, mall_ison

    Integer, Intent(In) :: n_levels ! refinement level
    Type (Spherical_Partition), Optional, Target :: partition
    Logical, Intent(In), Optional :: compress

! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::Initialize'
    Type (Spherical_Partition), Pointer  :: p_

    Integer :: n_coarse_regions, n_regions_all, n_fine_regions
    Integer :: idx, idxp, idx2, idx_mod
    Integer :: level, i, leaf_idx

    Type (RegionIterator) :: iter, iter2

    ASSERT_NOMSG(n_levels >= base_refinement)

    p_ => DefaultPartition(partition)
! should check, but doesn't work here because in fact
! p_%refinement_level = 0 is set by the F90 compiler at runtime (should
! really be undefined).
!   ALWAYS_ASSERT(.not. InitializedQ(p_))

    p_%refinement_level = n_levels
    P_%compress = .true.
    If (Present(compress)) p_%compress = compress

    ! We must ensure that m_icosahedron has been initialized.
    Call init_icosahedron()

    ! The number of regions on the finest partition
    n_fine_regions = n_refine**(n_levels+1) * n_faces

    ! The number of regions on all the coarser partitions combined
    n_coarse_regions  = (n_fine_regions - n_faces)/3

    ! The number of regions including all partitions - use
    ! simple geometric serios formula
    n_regions_all  = n_fine_regions + n_coarse_regions

        Allocate(p_%region(n_regions_all), p_%leaf(n_fine_regions), STAT = err)
    ALWAYS_ASSERT(err == 0)
	If (mall_ison()) Call mall_ci(44 * size(p_%region), myname)

    ! Use icosahedron to fill coarsest partition
    ! Recall that the numbering scheme is a bit odd.
    
    Call NewIter(iter, n_levels)
    idx = reset(iter, level = base_refinement)
    
    Do i = 1, n_faces
       p_%base_indices(i) = idx
       idx = next(iter, level = base_refinement)
       If (idx == DONE) exit
       ASSERT_NOMSG((idx > 0) .and. (idx <= n_regions_all))
    End Do
    
    p_%region(p_%base_indices)%triangle = GetFaceTriangles()

    If (n_levels == base_refinement) Then
       ! leaves are the regions
       p_%leaf = p_%base_indices
    End If

    leaf_idx = 0
    Do level = base_refinement, n_levels - 1 ! loop over refinement levels
       ! index on the coarser (and finer) regions
       idx  = reset(iter, level = level)

       Do ! loop over coarse regions
	  Do i = 1, n_refine
	     idxp = next(iter, level = level+1) ! next fine region
	     p_%region(idx)%child_index(i) = idxp
	     ASSERT_NOMSG((idxp > 0) .and. (idxp <= size(p_%region)))
	  End Do

	  p_%region(p_%region(idx)%child_index)%triangle = Refine(p_%region(idx)%triangle, &
	       & tests = p_%region(idx)%test)
	  p_%region(p_%region(idx)%child_index)%parent_index = idx

	  If (level == n_levels-1) Then
	     p_%leaf(leaf_idx+1:leaf_idx+n_refine) = p_%region(idx)%child_index(:)
	     leaf_idx = leaf_idx + n_refine
	  End If

	  idx = next(iter, level) ! next coarse region
	  If (idx == DONE) Exit ! Done

	  ASSERT(idx > 0)

       End Do

    End Do

  End Subroutine Initialize

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean - clean a partition
!
! !DESCRIPTION:
!    This subroutine cleans up what initialize creates
!
! !INTERFACE:

  Subroutine clean(partition)
    Use m_mall, only: mall_co, mall_ison
    Type (Spherical_Partition), Intent(InOut), Target, Optional :: partition
! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________
    Character(Len=*), Parameter :: myname_=myname//'::clean'
    Type (Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)

	If (mall_ison()) Call mall_co(44 * size(p_%region),myname)
    Deallocate(p_%region, STAT = err)
        ALWAYS_ASSERT(err == 0)
    p_%refinement_level = INVALID

  End Subroutine clean

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: BuildMask
!
! !DESCRIPTION:
!  This subroutine constructs a logical mask based upon the separation
!  angle between two regions.  The mask is true if at least some portion
!  of the regions are within the angle given by "threshold".  threshold
!  is specified in degrees _not_ radians.
!
! !INTERFACE:

  Subroutine BuildMask(threshold, mask, partition)
    Use m_RegionIterator, only : GetNumberOfSubregions, RegionIterator, clone
    Use m_RegionIterator, only : reset, clean, NewIter => new, next, DONE
    Use m_RegionIterator, only : Rotate, Reflect
    Use m_Icosahedron, only : n_faces, archetype, rotatetostandard
    Use m_mall, only : mall_mci, mall_ci, mall_mco, mall_co, mall_ison

    Real, Intent(In) :: threshold ! in degrees
    Logical(Kind=1), Target :: mask(:,:)
    Type (Spherical_Partition), Intent(In), Optional, Target :: partition

! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________
    Character(Len=*), Parameter :: myname_=myname//'::BuildMask'

    Type (Spherical_Partition), Pointer :: p_

    Real :: cos_threshold
    Integer :: idx_1, idx_2, n_dcnts
    Type (RegionIterator) :: iter_1, iter_2
    Integer :: i, j
    Integer :: f1, f2, f1p, f2p, rot_1, rot_2
    Logical :: sym
    Integer :: idx_1p, idx_2p
    Logical :: scratch(n_faces, n_faces)
    Integer, Allocatable :: permutations(:,:,:)
    Logical(Kind=1), Pointer :: ref(:,:)

    p_ => DefaultPartition(partition)

    ASSERT_NOMSG(p_%refinement_level >= base_refinement)
    ASSERT_NOMSG(size(mask,1) == NumberOfRegions(partition))
    ASSERT_NOMSG(size(mask,2) == NumberOfRegions(partition))

    ! cosine of threshold given in degrees
    cos_threshold = Cos(4*atan(1.)*threshold/180)

    ! Build the permutations for the 6 symmetries of a triangle.
    ! 3 rotations x 2 reflections - Identity
    ! n_dcnts is the number of subregions per icosahedral face
    Call NewIter(iter_1, p_%refinement_level)
    idx_1 = reset(iter_1)
    n_dcnts = GetNumberOfSubregions(iter_1)

    Allocate(permutations(n_dcnts+1,0:2,0:1), STAT = err)
        ALWAYS_ASSERT(err == 0)
	If (mall_ison()) Call mall_mci(permutations,myname)

    If (p_%compress) Then
       Allocate(ref(size(p_%region),size(p_%region)), STAT=err)
          ALWAYS_ASSERT(err ==  0)
          If (mall_ison()) Call mall_ci(size(ref)/4, myname)
    Else
       ref => mask
    End If

    Do
       Call clone(iter_1, iter_2)
       idx_2 = idx_1
       Do rot_1 = 0, 2

	  permutations(idx_1,rot_1,0) = idx_2 
	  idx_2 = Reflect(iter_2)
	  permutations(idx_1,rot_1,1) = idx_2 

	  idx_2 = Reflect(iter_2) ! undo the reflection
	  idx_2 = Rotate(iter_2)  ! and rotate
       End Do
       Call Clean(iter_2)
       idx_1 = next(iter_1)
       If (iter_1%refinement_path(base_refinement) == 2) Exit

    End Do
    Call Clean(iter_1)

    ! First build the archetype cases for face 1
    f1 = 1
    Call NewIter(iter_1,p_%refinement_level)
    Call NewIter(iter_2,p_%refinement_level)
    scratch = .false.
#ifndef NDEBUG
       ref = .false.
#endif
    Do j = 1, 6
       idx_1 = reset(iter_1, index = 1)
       f2 = archetype(j)
       idx_2 = 1 + (f2-1)* (n_dcnts+1)
       idx_2 = reset(iter_2, index = idx_2)
       Call BuildSubMask(iter_1, iter_2)
       scratch(1,f2) = .true.
    End Do

    Do f1 = 1, n_faces
       idx_1 = 1 + (f1-1)*(n_dcnts+1)

       Do f2 = 1, n_faces
	  idx_2 = 1 + (f2-1)*(n_dcnts + 1)

	  If (scratch(f1,f2)) Cycle
	     
	  Call RotateToStandard(f1, f2, f2p, rot_1, rot_2, sym)

	  ASSERT_NOMSG(Count(f2p==archetype) == 1)
	  Call Fill()
	  
	  scratch(f1,f2) = .true.

       End Do
    End Do
    
    ASSERT(all(scratch))

       If (mall_ison()) call mall_mco(permutations,myname)
    Deallocate(permutations, STAT = err)
       ALWAYS_ASSERT(err == 0)

    Call Clean(iter_1)
    Call Clean(iter_2)

    If (p_%compress) Then
       mask = ref(p_%leaf,p_%leaf)
          If (mall_ison()) call mall_co(size(ref)/4,myname)
       Deallocate(ref, STAT = err)
          ALWAYS_ASSERT(err == 0)
    End If

    ASSERT_NOMSG(All(mask .eqv. Transpose(mask)))

       

  Contains

    Subroutine Fill()
      Integer :: idx_1p, idx_2p, isym

      idx_1p = 1
      idx_2p = 1 + (f2p-1)*(n_dcnts+1)

      isym = 0
      IF (sym) isym = 1

      ref(idx_2:idx_2+n_dcnts,idx_1:idx_1+n_dcnts) = &
	   & ref(idx_2p-1+permutations(:,rot_2,isym),idx_1p-1+permutations(:,rot_1,isym))

    End Subroutine Fill

    Recursive Subroutine BuildSubMask(iter_1, iter_2, terminate)
      Use m_Spherical_Triangle, only : CosSeparation, SEPANG_MAX, SEPANG_MIN, GetVertices, Separation
      Use m_RegionIterator, only : IsLeaf, GetLevel
      Type (RegionIterator), Intent(InOut) :: iter_1, iter_2
      Logical, Optional, Intent(In) :: terminate
      Type (RegionIterator)                :: iter_1p, iter_2p
      Logical :: test_maxsep, test_minsep
      Integer :: n_dcnts_1, n_dcnts_2
      Integer :: level_1, level_2, idx_1, idx_2
      Integer :: i, j
      Logical, Parameter :: FALSE = .false.
      Logical, Parameter :: TRUE = .true.
      Logical :: terminate_

      idx_1 = iter_1%index
      idx_2 = iter_2%index

      level_1 = GetLevel(iter_1)
      level_2 = GetLevel(iter_2)

      terminate_ = .false.
      If (Present(terminate)) terminate_ = terminate

      test_maxsep = CosSeparation(p_%region(idx_1)%triangle, p_%region(idx_2)%triangle, SEPANG_MAX,cos_threshold)

      If (.not. test_maxsep) Then
	 ! all descendents are within the bounds

	 n_dcnts_1 = GetNumberOfSubregions(iter_1)
	 n_dcnts_2 = GetNumberOfSubregions(iter_2)
	 If (terminate_) n_dcnts_2 = 0

#ifndef NDEBUG
!!$	 If (f2 == 8) Then
!!$	    Print*,'Build: A',idx_1,idx_2,n_dcnts_1,n_dcnts_2
!!$	    Print*,'ref 1:',iter_1%refinement_path
!!$	    Print*,'ref 2:',iter_2%refinement_path
!!$	    Print*,separation(GetRegion(12),GetRegion(149),method = SEPANG_MIN)
!!$	    Print*,separation(GetRegion(149),GetRegion(12),method = SEPANG_MIN)
!!$	    Print*,separation(GetRegion(1),GetRegion(148),method = SEPANG_MAX)
!!$	 End If
#endif
	 ref(idx_2:idx_2+n_dcnts_2, idx_1:idx_1+n_dcnts_1) = TRUE

      Else
	 test_minsep = CosSeparation(p_%region(idx_1)%triangle, p_%region(idx_2)%triangle, SEPANG_MIN, cos_threshold)

	 If (.not. test_minsep) Then
	    ! no descendants are within the bounds
	    n_dcnts_1 = GetNumberOfSubregions(iter_1)
	    n_dcnts_2 = GetNumberOfSubregions(iter_2)
	    If (terminate_) n_dcnts_2 = 0
#ifndef NDEBUG
!!$	 If (f2 == 8) Then
!!$	    Print*,'Build: B',idx_1,idx_2,n_dcnts_1,n_dcnts_2
!!$	    Print*,'ref 1:',iter_1%refinement_path
!!$	    Print*,'ref 2:',iter_2%refinement_path
!!$	 End If
#endif
	    ref(idx_2:idx_2+n_dcnts_2, idx_1:idx_1+n_dcnts_1) = FALSE
	 Else
	    ! possibly some are in each category
	    ! construct mutually exclusive list of subpartitions 

	    ASSERT(.not. ref(idx_2,idx_1))
#ifndef NDEBUG
!!$	 If (f2 == 8) Then
!!$	    Print*,'Build: C',idx_1,idx_2
!!$	    Print*,'ref 1:',iter_1%refinement_path
!!$	    Print*,'ref 2:',iter_2%refinement_path
!!$	 End If
#endif

	    ref(idx_2,idx_1) = TRUE
	    
	    Call clone(iter_1, iter_1p)
	    If (.not. IsLeaf(iter_1)) Then
	       Do i = 1, n_refine
		  idx_1 = next(iter_1p, level = level_1 + 1)
		  Call BuildSubMask(iter_1p, iter_2, terminate = .true.)
	       End Do
	    End If
	    Call clean(iter_1p)
	    
	    If ((.not. terminate_) .and. (.not. IsLeaf(iter_2))) Then
	       Call clone(iter_2, iter_2p)
	       Do i = 1, n_refine
		  idx_2 = next(iter_2p, level = level_2 + 1)
		  Call BuildSubMask(iter_1, iter_2p)
	       End Do
	       Call clean(iter_2p)
	    End If

	 End If
      End If

    End Subroutine BuildSubMask

  End Subroutine BuildMask

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: NumberOfRegions
!
! !DESCRIPTION:
!
! !INTERFACE:

  Function NumberOfRegions(partition,atlevel) Result(n)
    Type (Spherical_Partition), Intent(In), Optional, Target :: partition
    integer,optional,intent(in) :: atlevel
    Integer :: n

! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________
    Character(Len=*), Parameter :: myname_=myname//'::NumberOfRegions'

    Integer :: level_

    Type (Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)

    if(present(atlevel)) then
      n = n_faces * (n_refine ** (atlevel - base_refinement))
      return
    endif

    If (p_%compress) Then
       n = n_faces * (n_refine ** (p_%refinement_level - base_refinement))
    Else
       n = size(p_%region)
    End If

  End Function NumberOfRegions
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 
  ! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Function NumberOfNodes(partition) Result(n)
    Use m_Icosahedron, only : n_edges, n_vertices
    Type (Spherical_Partition), Intent(In), Optional, Target :: partition
    Integer :: n

! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::NumberOfNodes'
    Type (Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)

    ! Geometric series
    n = n_vertices + &
	 & n_edges * (4**(p_%refinement_level-base_refinement) - 1)/(4-1)

  End Function NumberOfNodes

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 
  ! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Function NumberOfLeaves(partition) Result(n)
    Type (Spherical_Partition), Intent(In), Optional, Target :: partition
    Integer :: n

! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::NumberOfLeaves'
    Type (Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)

    ! Geometric series
    n = n_faces * n_refine ** (p_%refinement_level-base_refinement)

  End Function NumberOfLeaves

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GetLeaves
!
! !DESCRIPTION:
!    This routine bins xyz data into regien indices.
!
! !INTERFACE:

  Subroutine GetLeaves(indices, partition)
    Use m_RegionIterator, only : RegionIterator, IsLeaf, next, reset, clean, NewIter => new
    Integer, Intent(Out) :: indices(:)
    Type (Spherical_Partition), Intent(In), Optional, Target :: partition
! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::GetLeaves'

    Type (Spherical_Partition), Pointer :: p_
    Integer :: n, idx
    Type (RegionIterator) :: iter

    p_ => DefaultPartition(partition)

    Call NewIter(iter,p_%refinement_level)
    idx = reset(iter, level = p_%refinement_level)
    n = 0

    Do
       If (IsLeaf(iter)) Then
	  n = n + 1
	  ASSERT(n < size(indices))
	  indices(n) = idx
       End If

       idx = next(iter)
	  
       If (idx == INVALID) Exit
       ASSERT_NOMSG(idx > 0 .and. idx <= size(p_%region))
       
    End Do

    Call clean(iter)
    ASSERT_NOMSG(n == NumberOfLeaves(partition))

  End Subroutine GetLeaves

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ListCenters
!
! !DESCRIPTION:
!    This routine bins xyz data into regien indices.
!
! !INTERFACE:

  Subroutine ListCenters(centers, partition)
    Use m_Spherical_Triangle, only : Circumcenter
    Real, Intent(Out) :: centers(:,:)
    Type (Spherical_Partition), Intent(In), Optional, Target :: partition
! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::ListCenters'

    Integer :: i
    Type (Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)

    If (p_%compress) Then

       Do i = 1, NumberOfRegions(partition=p_)
	  centers(:,i) = Circumcenter(p_%region(p_%leaf(i))%triangle)
       End Do

    Else

       Do i = 1, NumberOfRegions(partition=p_)
	  centers(:,i) = Circumcenter(p_%region(i)%triangle)
       End Do

    End If

  End Subroutine ListCenters

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ListAreas
!
! !DESCRIPTION:
!    This routine bins xyz data into regien indices.
!
! !INTERFACE:

  Subroutine ListAreas(areas, partition)
    Use m_Spherical_Triangle, only : Area
    Real, Intent(Out) :: areas(:)
    Type (Spherical_Partition), Intent(In), Optional, Target :: partition
! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::ListAreas'

    Integer :: i
    Integer, Pointer :: indices_(:)
    Type (Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)

    If (p_%compress) Then
       Do i = 1, NumberOfRegions(partition)
	  areas(i) = Area(p_%region(p_%leaf(i))%triangle)
       End Do
    Else
       Do i = 1, NumberOfRegions(partition)
	  areas(i) = Area(p_%region(i)%triangle)
       End Do
    End If

  End Subroutine ListAreas


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ListNodes
!
! !DESCRIPTION:
!    This routine bins xyz data into regien indices.
!
! !INTERFACE:


!!$  Function ListNodes(partition) Result(nodes)
!!$    Type (Spherical_Partition), Intent(In), Optional, Target :: partition
!!$    Real :: nodes(3,NumberOfNodes(partition))
!!$
! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________
!!$
!!$    Character(Len=*), Parameter :: myname_=myname//'::ListNodes'
!!$    Type (Spherical_Partition), Pointer :: p_
!!$    Real, Allocatable :: vc(:,:)
!!$    Integer, Allocatable :: v(:,:)
!!$    Integer, Allocatable :: e(:,:)
!!$    Integer, Allocatable :: f(:,:)
!!$    
!!$    If (Present(partition)) Then
!!$       p_ => partition
!!$    Else
!!$       p_ => private_partition
!!$    End If
!!$
!!$    Allocate(vc(3,NumberOfNodes(partition)), STAT = err) ! vertices
!!$    Allocate(v(2,NumberOfNodes(partition)), STAT = err) ! vertices
!!$    Allocate(e(2,NumberOfEdges(partition)), STAT = err) ! edges
!!$    Allocate(f(3,NumberOfEdges(partition)), STAT = err) ! nodes
!!$
!!$
!!$    vc(:,:n_vertices) = GetVertices()
!!$    v(:,:)    = Reshape((/ (i,i,i = 1, n_vertices) /), SHAPE = (/2, n_vertices/))
!!$    e(:,:n_edges)    = GetEdges()
!!$    f(:,:n_faces)    = GetFaces()
!!$
!!$    nv = n_vertices
!!$    ne = n_edges
!!$    nf = n_faces
!!$
!!$    Do level = base_refinement+1, p_%refinement_level
!!$       nv_old = nv
!!$       ! Add midpoints to the vert. list and new "interior" edges
!!$       Do i = 1, ne
!!$	  v(:,nv+1) = (/ e(1,i), e(2,i) /)
!!$	  vc(:,nv+1) = UnitVector(vc(:,e(1,i)) + vc(:,e(2,i)))
!!$	  nv = nv + 1
!!$	  ASSERT_NOMSG(nv < NumberOfNodes(partition))
!!$       End Do
!!$
!!$       nf_old = nf
!!$       Do i = 1, nf_old
!!$	  Do j = 1, 3 ! loop over edges
!!$	     jp = 1 + mod(j+3-2,3) ! next edge
!!$
!!$	     ! Find the new vertex that is the midpoint
!!$	     k = nv_old + 1
!!$	     Do 
!!$		If ((v(1,k) == Min(f(j,i),f(jp,i))) .and. (v(2,k) == Max(f(j,i),f(jp,i)))) Then
!!$		   midpt(j) = k
!!$		   Exit
!!$		End If
!!$		k = k + 1
!!$		ASSERT(k <= nv)
!!$	     End Do
!!$	  End Do
!!$
!!$	  ! Add "interior" edges
!!$	  Do j = 1, 3
!!$	     jp = 1 + mod(j+3-2,3) ! next edge
!!$	     e(1,ne+1) = Min(midpt(j),midpt(jp))
!!$	     e(2,ne+1) = Max(midpt(j),midpt(jp))
!!$	     ne = ne + 1
!!$	     ASSERT_NOMSG(ne <= NumberOfEdges(partition))
!!$	  End Do
!!$	  
!!$	  ! Split original edges
!!$	  Do j = 1, 3
!!$	     jp = 1 + mod(j+3-2,3) ! next edge
!!$	     If (f(j,i) < f(jp, i)) Then
!!$		! Find original edge
!!$		Do
!!$		   If ((e(1,k) == f(j,i)) .and (e(2,k) == f(jp,i))) Exit
!!$		   ASSERT(k <= ne_old)
!!$		End Do
!!$		
!!$		e(1,ne+1) = Min(f(jp,i),mipdt(j))
!!$		e(2,ne+1) = Max(f(jp,i),mipdt(j))
!!$		ne = ne + 1
!!$
!!$		e(1,k) = Min(e(1,k), mipdt(j)
!!$		e(2,k) = Max(e(1,k), midpt(j))
!!$		! not a new one, so do not increment ne
!!$	     End If
!!$	  End Do
!!$
!!$	  ! new faces
!!$	  Do j = 1, 3
!!$	     jp = 1 + mod(j+3-2,3) ! next edge
!!$	     f(:,nf + 1) = (/ f(j,i), midpt(j), midpt(jp) /)
!!$	     nf = nf + 1
!!$	     ASSERT_NOMSG(nf <= NumberOfRegions(partition))
!!$	  End Do
!!$	  f(:,i) = midpt
!!$
!!$       End Do
!!$
!!$    End Do
!!$
!!$    nodes = vc
!!$
!!$    Deallocate(vc, v, e, f, STAT = err)
!!$    ASSERT(err == 0)
!!$
!!$  End Function ListNodes

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: xyz2reg - bin points into regions
!
! !DESCRIPTION:
!    This routine bins xyz data into regien indices.
!
! !INTERFACE:

  Subroutine xyz2reg(n, x, y, z, idx, partition, atlevel)
    Use m_Icosahedron, only : icos_xyz2reg => xyz2reg
    Use m_RegionIterator, only : coef
    Use m_die, only : die
    Integer, Intent(In) :: n
    Real, Intent(In) :: x(n), y(n), z(n)
    Integer, Intent(Out) :: idx(n)
    Type (Spherical_Partition), Target, Intent(In), Optional :: partition
    Integer, Optional, Intent(In) :: atlevel

! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::xyz2reg'

    Type (Spherical_Partition), Pointer :: part_
    Integer :: i, level, ireg, f, idx_old, iregp
    Real :: p(3), x_old, y_old, z_old, test
    Integer :: atlevel_

    part_ => DefaultPartition(partition)
    ASSERT(InitializedQ(part_))

    ! No loop over refinement levels to smaller and smaller bins
    atlevel_ = part_%refinement_level
    If (Present(atlevel)) Then
       If (atlevel > atlevel_) Call die(myname_,'insufficient refinement for level',atlevel)
       atlevel_ = atlevel
    End If

    x_old = -2. ! impossible value to initialize
    y_old = -2.
    z_old = -2.
    idx_old = INVALID

    Do i = 1, n
       ! First check to see if this obs is in the same position
       ! as the old one.
       If ((x(i)-x_old)**2 + (y(i)-y_old)**2 + (z(i)-z_old)**2 == 0) Then
	  idx(i) = idx_old
	  Cycle
       End If

       ! Construct a vector from the coordinates for convenience.
       p = (/ x(i), y(i), z(i) /)

       ! Now determine which icosahedral region the point belongs to:
       f = icos_xyz2reg(p)

       ! Unfortunately, we must now convert to the awkward numbering scheme.
       ireg = 1 + (f-1) * coef(part_%refinement_level)
       iregp = f - 1

       Do level = base_refinement, atlevel_ - 1

	  ! could be made into a loop, but unrolled for efficiency.
	  test = dot_product(p,part_%region(ireg)%test(:,1))
	  If (test >= 0) Then
	     ireg = part_%region(ireg)%child_index(1)
	     iregp = n_refine*iregp + 0
	  Else
	     test = dot_product(p,part_%region(ireg)%test(:,2))
	     If (test >= 0) Then
		ireg = part_%region(ireg)%child_index(2)
		iregp = n_refine*iregp + 1
	     Else
		test = dot_product(p,part_%region(ireg)%test(:,3))
		If (test >= 0) Then
		   ireg = part_%region(ireg)%child_index(3)
		   iregp = n_refine*iregp + 2
		Else
		   ireg = part_%region(ireg)%child_index(4)
		   iregp = n_refine*iregp + 3
		End If
	     End If
	  End If
       End Do

       If (part_%compress) Then
	  idx(i) = iregp + 1
	  idx_old = iregp + 1
       Else
	  idx(i) = ireg
	  idx_old = ireg
       End If

       x_old = x(i)
       y_old = y(i)
       z_old = z(i)

    End Do

  End Subroutine xyz2reg
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: LatLon_to_Region - bin points into regions
!
! !DESCRIPTION:
!    This routine bins LatLon data into regien indices.
!
! !INTERFACE:

  Subroutine LatLon_to_Region(n, lat, lon, idx, ier, partition, atlevel)
    Use m_Icosahedron, only : icos_xyz2reg => xyz2reg
    Use m_RegionIterator, only : coef
    Use m_stdio, only : stderr
    Use m_die, only : die,perr

    Integer, Intent(In) :: n
    Real, Intent(In) :: lat(n), lon(n)
    Integer, Intent(Out) :: idx(n)
    Integer, Intent(Out) :: ier
    Type (Spherical_Partition), Target, Intent(In), Optional :: partition
    Integer, Optional,Intent(In) :: atlevel

! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::LatLon_to_Region'
    Type (Spherical_Partition), Pointer :: part_
    Integer :: i, level, ireg, f, idx_old, iregp
    Real :: p(3), lat_old,lon_old,test
    Integer :: atlevel_
    Real :: rlon, rlat,clon,slon,clat,slat,deg2rad

!.......................................................................
!.... Statement function.

      logical::      llbad
      llbad(n) = abs(lat(n)).gt.90.0

    part_ => DefaultPartition(partition)
    ASSERT(InitializedQ(part_))

    
    
    lat_old = -1000. ! impossible value to initialize
    lon_old = -1000. ! impossible value to initialize
    idx_old = INVALID
    deg2rad=4.*atan(1.)/180.

    ! No loop over refinement levels to smaller and smaller bins
    atlevel_ = part_%refinement_level
    If (Present(atlevel)) Then
       If (atlevel > atlevel_) Call die(myname_,	&
       		'insufficient refinement for level',atlevel)
       atlevel_ = atlevel
    End If

    Do i = 1, n

       if( llbad(i) ) then
          Call perr(myname_,'latitude or longitude out of range',i)
          write(stderr,*) myname_,':  n = ', i,	&
	  	'  rlats(n) = ', lat(i),'  rlons(n) = ', lon(i)
            ier = 1
            Nullify(part_)
          return
       endif

       ! First check to see if this obs is in the same position
       ! as the old one.

       If (i>1 .and. lon_old == lon(i) .and. lat_old == lat(i)) Then
          idx(i) = idx_old
          Cycle
       End If

       ! Convert to x, y, z
       rlon=lon(i)*deg2rad
       rlat=lat(i)*deg2rad
       clon=cos(rlon)
       slon=sin(rlon)
       clat=cos(rlat)
       slat=sin(rlat)

       ! Construct a vector from the coordinates for convenience.
       p = (/ clat*clon, clat*slon, slat /)


       ! Now determine which icosahedral region the point belongs to:
       f = icos_xyz2reg(p)

       ! Unfortunately, we must now convert to the awkward numbering scheme.
       ireg = 1 + (f-1) * coef(part_%refinement_level)
       iregp = f - 1

       Do level = base_refinement, atlevel_ - 1

	  ! could be made into a loop, but unrolled for efficiency.
	  test = dot_product(p,part_%region(ireg)%test(:,1))
	  If (test >= 0) Then
	     ireg = part_%region(ireg)%child_index(1)
	     iregp = n_refine*iregp + 0
	  Else
	     test = dot_product(p,part_%region(ireg)%test(:,2))
	     If (test >= 0) Then
		ireg = part_%region(ireg)%child_index(2)
		iregp = n_refine*iregp + 1
	     Else
		test = dot_product(p,part_%region(ireg)%test(:,3))
		If (test >= 0) Then
		   ireg = part_%region(ireg)%child_index(3)
		   iregp = n_refine*iregp + 2
		Else
		   ireg = part_%region(ireg)%child_index(4)
		   iregp = n_refine*iregp + 3
		End If
	     End If
	  End If
       End Do

       If (part_%compress) Then
	  idx(i) = iregp + 1
	  idx_old = iregp + 1
       Else
	  idx(i) = ireg
	  idx_old = ireg
       End If

       lat_old = lat(i)
       lon_old = lon(i)

    End Do

    Nullify(part_)
    ier = 0

  End Subroutine LatLon_to_Region
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GetRegion
!
! !DESCRIPTION:
!    Return the triangle of the specified region index
!
! !INTERFACE:

  Function GetRegion(idx, partition) Result(tri)
    Integer, Intent(In) :: idx
    Type (Spherical_Partition), Target, Intent(In), Optional :: partition
    Type (Spherical_Triangle) :: tri

! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::GetRegion'
    Type (Spherical_Partition), Pointer :: p_
    Integer :: idxp

    p_ => DefaultPartition(partition)
    ASSERT(InitializedQ(p_))
    ASSERT_NOMSG((idx > 0) .and. idx <= NumberOfRegions(partition))

    If (p_%compress) Then
       idxp = p_%leaf(idx)
       tri = p_%region(idxp)%triangle
    Else
       tri = p_%region(idx)%triangle
    End If

  End Function GetRegion

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: DefaultPartition
!
! !DESCRIPTION:
! A   very useful utility to allow optional "partition" arguments.
!
! !INTERFACE:

  Function DefaultPartition(partition) Result (p)
    Type (Spherical_Partition), Target, Intent(In), Optional :: partition
    Type (Spherical_Partition), Pointer :: p

! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::DefaultPartition'
    If (Present(partition)) Then
       p => partition
    Else
       p => private_partition
    End If

  End Function DefaultPartition


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: BaseRegion
!
! !DESCRIPTION:
!    Return the base region specified by level
!
! !INTERFACE:

  Function BaseRegion(idx, level, partition) Result(base_idx)
    Use m_RegionIterator, only : coef,base_refinement
    Integer, Intent(In) :: idx
    Integer, Intent(In) :: level
    Type (Spherical_Partition), Intent(In), Optional, Target :: partition
    Integer :: base_idx

! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(len=*), Parameter :: myname_ = myname//'::BaseRegion'
    Type (Spherical_Partition), Pointer :: p_
    Integer :: nlevel,ilevel,idigit,ipower,iremin
!------------------------------------------------------------------------------
    p_ => DefaultPartition(partition)

    base_idx=0
    if(level<base_refinement .or. level>p_%refinement_level) return

	! The returned base_idx value is 0, if either the level or the
	! region index is invalid.  For example, the specified input
	! level may be higher than base_refinement level or lower than
	! %refinement_level defining the DefaultPartition(), which is
	! handled there.  Or, an input region index idx may have a value
	! for a region on a level higher (smaller level value) then the
	! specified level, which is handled below.

    If (p_%compress) Then
      iremin=idx-1	! Shift the lower-bound of idx from 1 to 0.
      nlevel=p_%refinement_level+1
      do ilevel=nlevel,nlevel-level+base_refinement,-1
        ipower=n_refine**ilevel
        idigit=iremin/ipower
        iremin=mod(iremin,ipower)
        base_idx=base_idx*n_refine+idigit
      end do

    Else
      iremin=idx
      nlevel=p_%refinement_level
      do ilevel=nlevel,nlevel-level+base_refinement,-1

	if(iremin==0) then
		! This is a special case, where the index is on a level
		! higher than then specified level.
	  base_idx=0
	  return
	endif
		! Remove the place reserved for the parent
		! level.  For the first iteration of ilevel,
		! it shifts the lower-bound of idx from 1 to 0.
	iremin=iremin-1

        idigit=iremin/coef(ilevel)
        iremin=mod(iremin,coef(ilevel))

		! Add back as an index for a compressed partition.
        base_idx=base_idx*n_refine+idigit
      end do

    End If

    base_idx=base_idx+1	! The range of idx is [1..N], not [0..N-1]

  End Function BaseRegion


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Base20Region
!
! !DESCRIPTION:
!   Return the face of the underlying icosahedron that this region lies in
!
! !INTERFACE:

  Function Base20Region(idx, partition) Result(base_idx)
    Use m_RegionIterator, only : coef
    Integer, Intent(In) :: idx
    Type (Spherical_Partition), Intent(In), Optional, Target :: partition
    Integer :: base_idx

! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::Base20Region'
    Type (Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)

    If (p_%compress) Then
       base_idx = 1 + (idx-1)/n_refine**(p_%refinement_level-base_refinement)
    Else
       base_idx = 1 + (idx-1) / coef(p_%refinement_level)
    End If

  End Function Base20Region


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Base80Region
!
! !DESCRIPTION:
!    Return the region number corresponding to the original (Jim Fentner numbering).
!
! !INTERFACE:

  Function Base80Region(idx, partition) Result(base_idx)
    Use m_RegionIterator, only : coef
    Integer, Intent(In) :: idx
    Type (Spherical_Partition), Intent(In), Optional, Target :: partition
    Integer :: base_idx
! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::Base80Region'
    Integer :: face, subreg

    Type (Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)
    If (p_%compress) Then
       face = 1 + (idx-1)/(n_refine ** (p_%refinement_level - base_refinement))
       subreg = 1 + mod((idx-1)/(n_refine ** (p_%refinement_level - 1 - base_refinement)),n_refine)
    Else
       face = 1 + (idx-1)/coef(p_%refinement_level)
       subreg = 1 + (mod(idx-1,coef(p_%refinement_level))-1)/coef(p_%refinement_level-1)
    End If
       base_idx = 4*(face-1) + subreg

  End Function Base80Region

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: NumberOfEdges
!
! !DESCRIPTION:
!
! !INTERFACE:

  Function NumberOfEdges(partition) Result(n)
    Type (Spherical_Partition), Intent(In), Target, Optional :: partition
    Integer :: n
! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::NumberOfEdges'

    Type (Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)

    n = 30 * (n_refine ** (p_%refinement_level-base_refinement))

  End Function NumberOfEdges

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: InitializedQ
!
! !DESCRIPTION:
!
! !INTERFACE:

  Logical Function InitializedQ(partition)
    Type (Spherical_Partition), Intent(In), Target, Optional :: partition
! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::InitializedQ'

    Type (Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)

    InitializedQ = (p_%refinement_level >= base_refinement)

  End Function InitializedQ

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Get_by_partition
!
! !DESCRIPTION:
!
! !INTERFACE:

  Subroutine Get_by_partition(partition, refinement_level)
    Type (Spherical_Partition), Target, Intent(In), Optional :: partition
    Integer, Intent(Out) :: refinement_level
! !REVISION HISTORY:
! 	
!EOP ___________________________________________________________________

    Character(Len=*), Parameter :: myname_=myname//'::Get_by_partition'
    
    Type (Spherical_Partition), Pointer :: p_
    p_ => DefaultPartition(partition)
    ASSERT(InitializedQ(p_))

    refinement_level = p_%refinement_level

  End Subroutine Get_by_partition

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Get_by_index
!
! !DESCRIPTION:
!
! !INTERFACE:

  Recursive Subroutine Get_by_index(idx, parent, siblings, depth, part)
  ! Return the index of the parent region for the
  ! region specified by idx

    Integer, Intent(In) :: idx 
    Integer, Intent(Out), Optional :: parent
    Integer, Intent(Out), Optional :: siblings(:)
    Integer, Intent(Out), Optional :: depth
    Type (Spherical_Partition), Optional :: part

    Type (Spherical_Partition), Pointer :: p_
    Integer :: parent_, depth_

    p_ => DefaultPartition(part)
    ASSERT(InitializedQ(p_))
    ASSERT(idx > 0)
    ASSERT(idx <= NumberOfRegions())

    parent_ = p_%region(idx)%parent_index
    If (Present(parent)) parent = parent_

    If (parent_ /= INVALID) Then

       If (Present(siblings)) siblings = p_%region(parent_)%child_index
       If (Present(depth)) Then
          Call  Get_by_index(parent_, depth = depth_, part = part)
          depth = depth_ + 1 ! one deeper than parent
       End If
       
       
    Else

       If (Present(siblings)) siblings = -1
       If (Present(depth)) depth = -1 ! base is -1 level for historical reasons

    End If

  End Subroutine Get_by_index

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: new_ - create an object as a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    function new_(mold,stat)
      use m_die, only : die,perr
      use m_mall,only : mall_ison,mall_ci
      implicit none
      type(Spherical_Partition),pointer :: mold
      integer,optional,intent(out) :: stat
      type(Spherical_Partition),pointer :: new_

! !REVISION HISTORY:
! 	11Feb02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		. initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::new_'
  type(Spherical_Partition),pointer :: obj
  integer :: ier

  if(present(stat)) stat=0
  allocate(obj,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'allocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	if(mall_ison()) call mall_ci(1,myname)

  new_ => obj
  nullify(obj)		! to prevent the compiler touching the memory.
end function new_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: delete_ - delete an object as a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:
    subroutine delete_(obj,stat)
      use m_die, only : die,perr
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(Spherical_Partition),pointer :: obj
      integer,optional,intent(out) :: stat


! !REVISION HISTORY:
! 	11Feb02	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		. initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::delete_'
  integer :: ier

  if(present(stat)) stat=0
	if(mall_ison()) call mall_co(1,myname)

  deallocate(obj,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine delete_
End Module m_Spherical_Partition
