Module m_Spherical_Partition_Simplified
  Implicit None
  Private

  ! PUBLIC INTERFACES
  Public :: Spherical_Partition
  Public :: Initialize
  Public :: Clean
  Public :: NumberOfRegions
  Public :: xyz2reg
  Public :: InitializedQ
  Public :: GetRegion
  Public :: Separation  ! Added Separation to public interface
 
  Public :: SEPANG_MIN 
  ! CONSTANTS
  ! From m_Icosahedron
  Integer, Parameter :: n_faces     = 20
  Integer, Parameter :: n_edges     = 30
  Integer, Parameter :: n_vertices  = 12
  Integer, Parameter :: PrimaryFace = 7 ! arbitrary (but unfortunate choice)
  
  ! From m_RegionIterator
  Integer, Parameter :: n_refine = 4 ! subregions per region
  Integer, Parameter :: INVALID = -9999
  Integer, Parameter :: DONE    = -1
  Integer, Parameter :: base_refinement = -1
  Integer, Parameter :: max_refinement = 7
  
  ! From m_Spherical_Triangle - separation methods
  Integer, Parameter :: SEPANG_MIN            = 1
  Integer, Parameter :: SEPANG_MAX            = 2
  Integer, Parameter :: SEPANG_CIRCUMCENTER   = 3
  Integer, Parameter :: SEPANG_CENTER_OF_MASS = 4
  Integer, Parameter :: sep_methods(4) = &
       &(/ SEPANG_MIN, SEPANG_MAX, SEPANG_CIRCUMCENTER, SEPANG_CENTER_OF_MASS /)

  ! DERIVED TYPES
  ! From m_Spherical_Triangle
  Type Spherical_Triangle
     Private
     SEQUENCE
     Real :: A(3) = 0 ! vertex A
     Real :: B(3) = 0 ! vertex B
     Real :: C(3) = 0 ! vertex C
  End Type Spherical_Triangle
  
  ! From m_RegionIterator
  Type RegionIterator
     Integer, Pointer :: refinement_path(:) => Null()
     Integer :: index = INVALID
  End Type RegionIterator
  
  ! From m_Spherical_Partition
  Type Region
     Type (Spherical_Triangle)     :: triangle
     Real, Dimension(3,n_refine-1) :: test
     Integer :: parent_index = INVALID
     Integer :: child_index(n_refine) = INVALID
  End Type Region

  Type Spherical_Partition
     Private
     Integer                :: refinement_level = INVALID
     Type (Region), Pointer :: region(:) => Null()
     Integer                :: base_indices(n_faces) = INVALID
     Logical                :: Compress = .false.
     Integer, Pointer       :: leaf(:) => Null()
  End Type Spherical_Partition

  ! Default partition for when no partition is specified
  Type (Spherical_Partition), Target :: private_partition

  ! Helper arrays (from m_RegionIterator)
  Integer, Parameter :: coef(base_refinement:max_refinement) = &
       & (/ 1, 5, 21, 85, 341, 1365, 5461, 21845, 87381 /)
  
  ! Internal variables
  Integer :: err ! used for internal error handling
  Logical :: icosahedron_initialized = .false.
  Real :: vertices(3,n_vertices) ! Store icosahedron vertices
  
  ! Variables for xyz2reg
  Real :: octahedral_axis_1(3)
  Real :: octahedral_axis_2(3)
  Real :: octahedral_axis_3(3)
  Real :: refine_a(3)
  Real :: refine_b(3)
  Real :: refine_c(3)

  ! Internal array for face vertices (from m_Icosahedron)
  Integer, Parameter :: face_vertices(3,n_faces) =  &
       & RESHAPE(SOURCE = &
       & (/  1,  2,  6,     1,  3,  2,     1,  4,  3,     1,  5,  4, &
       &     1,  6,  5,     7,  6,  2,     8,  2,  3,     9,  3,  4, &
       &    10,  4,  5,    11,  5,  6,     2,  8,  7,     3,  9,  8, &
       &     4, 10,  9,     5, 11, 10,     6,  7, 11,    12,  7,  8, &
       &    12,  8,  9,    12,  9, 10,    12, 10, 11,    12, 11,  7 /), &
       & SHAPE = (/ 3, n_faces /))

Contains

  !=====================================================================
  ! FUNCTIONS FROM m_Spherical_Triangle
  !=====================================================================

  ! Create a new triangle from three points on the unit sphere
  Function Construct_Triangle(v1, v2, v3) Result(tri)
    Real, Intent(In) :: v1(3), v2(3), v3(3)
    Type(Spherical_Triangle) :: tri
    
    tri%A = v1
    tri%B = v2
    tri%C = v3
  End Function Construct_Triangle

  ! Get the vertices of a triangle
  Function GetVertices_Triangle(tri) Result(v)
    Type(Spherical_Triangle), Intent(In) :: tri
    Real :: v(3,3)
    
    v(:,1) = tri%A
    v(:,2) = tri%B
    v(:,3) = tri%C
  End Function GetVertices_Triangle
  
  ! Normalize a vector to unit length
  Function UnitVector(v) Result(u)
    Real, Intent(In) :: v(:)
    Real :: u(size(v))
    Real :: norm
    
    norm = sqrt(sum(v**2))
    If (norm > 0.0) Then
      u = v / norm
    Else
      u = 0.0
    End If
  End Function UnitVector
  
  ! Cross product of two 3D vectors
  Function Cross(a, b) Result(c)
    Real, Intent(In) :: a(3), b(3)
    Real :: c(3)
    
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  End Function Cross
  
  ! Safe version of acos to handle small numerical errors
  Function Safe_Acos(x) Result(angle)
    Real, Intent(In) :: x
    Real :: angle
    
    If (x >= 1.0) Then
      angle = 0.0
    Else If (x <= -1.0) Then
      angle = 4.0 * atan(1.0)  ! pi
    Else
      angle = acos(x)
    End If
  End Function Safe_Acos

  ! Get the center of mass of a triangle
  Function CenterOfMass(tri) Result(center)
    Type(Spherical_Triangle), Intent(In) :: tri
    Real :: center(3)
    
    center = tri%A + tri%B + tri%C
    center = UnitVector(center)
  End Function CenterOfMass

  ! Get the circumcenter of a triangle
  Function Circumcenter(tri, cosrad) Result(center)
    Type(Spherical_Triangle), Intent(In) :: tri
    Real, Optional :: cosrad   ! cosine of the radius of the circumcircle
    Real :: center(3)
    
    center = Cross(tri%B, tri%A) + Cross(tri%C, tri%B) + Cross(tri%A, Tri%C)
    center = UnitVector(center)
    
    If (Present(cosrad)) cosrad = dot_product(tri%A, center)
  End Function Circumcenter

  ! Refine a triangle into four subtriangles
  Function Refine(tri, tests) Result(daughters)
    Type(Spherical_Triangle), Intent(In) :: tri
    Real, Optional, Intent(Out) :: tests(3,3)
    Type(Spherical_Triangle) :: daughters(4)
    
    Real :: AB(3), BC(3), CA(3)
    Real :: A(3), B(3), C(3)
    
    A = tri%A
    B = tri%B
    C = tri%C
    
    ! Find midpoints on great circles
    AB = UnitVector(A+B)
    BC = UnitVector(B+C)
    CA = UnitVector(C+A)
    
    ! Create four subtriangles
    daughters(4) = Spherical_Triangle(BC, CA, AB)
    daughters(1) = Spherical_Triangle(A, AB, CA)
    daughters(2) = Spherical_Triangle(B, BC, AB)
    daughters(3) = Spherical_Triangle(C, CA, BC)
    
    If (Present(tests)) Then
       tests(:,1) = Cross(CA, AB)
       tests(:,2) = Cross(AB, BC)
       tests(:,3) = Cross(BC, CA)
    End If
  End Function Refine
  
  ! Trivial function to select a vertex from a triangle
  Function vertex(tri,i) Result(v)
    Type(Spherical_Triangle), Intent(In) :: tri
    Integer, Intent(In) :: i
    Real :: v(3)
    
    Select Case(i)
    Case(1)
       v = tri%A
    Case(2)
       v = tri%B
    Case(3)
       v = tri%C
    Case Default
       v = 0.0  ! Invalid index
    End Select
  End Function vertex
  
  ! Function to select an arc (specified by opposing vertex) from a triangle
  Function arc(tri,i) Result(vv)
    Type(Spherical_Triangle), Intent(In) :: tri
    Integer, Intent(In) :: i
    Real :: vv(3,2)
    
    Select Case(i)
    Case(1)
       vv(:,1) = tri%B
       vv(:,2) = tri%C
    Case(2)
       vv(:,1) = tri%C
       vv(:,2) = tri%A
    Case(3)
       vv(:,1) = tri%A
       vv(:,2) = tri%B
    Case Default
       vv = 0.0  ! Invalid index
    End Select
  End Function arc
  
  ! Compute minimum separation between a point and an arc
  Function arc_sep_min(p, arc_ab, inside) Result(cosang)
    Real, Intent(In) :: p(3)       ! The point
    Real, Intent(In) :: arc_ab(3,2) ! The arc
    Real, Intent(InOut) :: inside
    Real :: cosang
    
    Real :: n(3)    ! The normal vector to the great circle containing the arc
    
    n = Cross(arc_ab(:,1), arc_ab(:,2))
    inside = inside + dot_product(p, n)
    
    If (dot_product(arc_ab(:,1), Cross(p, n)) > 0) Then
       ! A is on the "right" side 
       If (dot_product(arc_ab(:,2), Cross(p, n)) < 0) Then
          ! B is on the "left" side and the
          ! perpendicular distance to the arc is the shortest
          cosang = sqrt(max(0., 1 - dot_product(p, UnitVector(n))**2)) 
       Else ! B must be the closer point
          cosang = dot_product(p, arc_ab(:,2)) 
       End If
    Else ! A is on the "left" side
       If (dot_product(arc_ab(:,2), Cross(p, n)) < 0) Then
          ! B is on the "right" side
          cosang = dot_product(p, arc_ab(:,1)) ! A must be the closer point
       Else ! arc is on the far side of the sphere 
          ! must check both endpoints
          cosang = max(dot_product(p, arc_ab(:,1)), dot_product(p, arc_ab(:,2)))
       End If
    End If
    
    cosang = min(cosang, +1.0) ! prevent roundoff errors
  End Function arc_sep_min
  
  ! Compute maximum separation between a point and an arc
  Function arc_sep_max(p, arc_ab, inside) Result(cosang)
    Real, Intent(In) :: p(3)       ! The point
    Real, Intent(In) :: arc_ab(3,2) ! The arc
    Real, Intent(InOut) :: inside
    Real :: cosang
    
    Real :: n(3)    ! The normal vector to the great circle containing the arc
    
    n = Cross(arc_ab(:,1), arc_ab(:,2))
    inside = inside + dot_product(p, n)
    
    If (dot_product(arc_ab(:,1), Cross(p, n)) > 0) Then
       ! A is on the "right" side 
       If (dot_product(arc_ab(:,2), Cross(p, n)) < 0) Then
          ! B is on the "left" side and the
          ! the farthest distance is one of the endpoints
          cosang = min(dot_product(p, arc_ab(:,1)), dot_product(p, arc_ab(:,2)))
       Else ! A must be the further point
          cosang = dot_product(p, arc_ab(:,1)) 
       End If
    Else ! A is on the "left" side
       If (dot_product(arc_ab(:,2), Cross(p, n)) < 0) Then
          ! B is on the "right" side
          cosang = dot_product(p, arc_ab(:,2)) ! B must be the furthest point
       Else ! arc is on the far side of the sphere 
          ! perpendicular distance to the arc is the furthest
          cosang = -sqrt(max(0., 1 - dot_product(p, UnitVector(n))**2))
       End If
    End If
    
    cosang = max(cosang, -1.0) ! prevent roundoff errors
  End Function arc_sep_max
  
  ! Calculate area of a spherical triangle
  Function Area_(tri) Result(A)
    Type(Spherical_Triangle), Intent(In) :: tri
    Real :: A
    
    Real, Dimension(3) :: U_ab, U_ac, U_ba, U_bc, U_ca, U_cb
    Real :: pi
    
    ! U_ij is the unit vector at vertex i that is tangent to the 
    ! equator passing through arc ij.
    U_ab = UnitVector(Cross(tri%A, Cross(tri%B, tri%A)))
    U_ac = UnitVector(Cross(tri%A, Cross(tri%C, tri%A)))
    
    U_ba = UnitVector(Cross(tri%B, Cross(tri%A, tri%B)))
    U_bc = UnitVector(Cross(tri%B, Cross(tri%C, tri%B)))
    
    U_ca = UnitVector(Cross(tri%C, Cross(tri%A, tri%C)))
    U_cb = UnitVector(Cross(tri%C, Cross(tri%B, tri%C)))
    
    ! The area of a spherical triangle is 2π minus the sum of the angles
    pi = 4.0 * atan(1.0)
    
    A = Safe_Acos(dot_product(U_ab, U_ac)) + &
        Safe_Acos(dot_product(U_ba, U_bc)) + &
        Safe_Acos(dot_product(U_ca, U_cb)) - pi
  End Function Area_
  
  ! Calculate separation angle between two triangles
  Function Separation(tri_1, tri_2, method) Result(angle)
    Type(Spherical_Triangle), Intent(In) :: tri_1
    Type(Spherical_Triangle), Intent(In) :: tri_2
    Integer, Intent(In) :: method
    Real :: angle
    
    Real :: c1(3), c2(3)
    Real :: cosang, cosang_max, cosang_min
    Real :: inside(3,2)
    Real :: pi
    Integer :: i, j
    
    ! Check if method is valid
    If (.not. any(method == sep_methods)) Then
      Print *, "Error: Invalid separation method"
      angle = -999.0
      Return
    End If
    
    Select Case(method)
    Case(SEPANG_MIN)
       cosang_max = -1.0 ! worst case
       inside = 0.0
       Do j = 1, 3
          Do i = 1, 3
             ! Compute distance between i-th vertex of tri_1 and
             ! side (arc) j on tri_2
             cosang = arc_sep_min(vertex(tri_1, i), arc(tri_2, j), inside(j,1))
             cosang_max = max(cosang_max, cosang)
             
             ! Same with triangles swapped
             cosang = arc_sep_min(vertex(tri_2, i), arc(tri_1, j), inside(j,2))
             cosang_max = max(cosang_max, cosang)
          End Do
       End Do
       
       ! Check if triangles are nested
       If (Any(All(inside < 0.0, 1))) cosang_max = +1.0
       
       cosang = cosang_max
       
    Case(SEPANG_MAX)
       cosang_min = +1.0 ! worst case
       inside = 0.0
       Do j = 1, 3
          Do i = 1, 3
             ! Compute distance between i-th vertex of tri_1 and
             ! side (arc) j on tri_2
             cosang = arc_sep_max(vertex(tri_1, i), arc(tri_2, j), inside(j,1))
             cosang_min = min(cosang_min, cosang)
             
             ! Same with triangles swapped
             cosang = arc_sep_max(vertex(tri_2, i), arc(tri_1, j), inside(j,2))
             cosang_min = min(cosang_min, cosang)
          End Do
       End Do
       
       ! Check if triangles are nested
       If (Any(All(inside < 0.0, 1))) cosang_min = +1.0
       
       cosang = cosang_min
       
    Case(SEPANG_CIRCUMCENTER)
       c1 = Circumcenter(tri_1)
       c2 = Circumcenter(tri_2)
       cosang = dot_product(c1, c2)
       
    Case(SEPANG_CENTER_OF_MASS)
       c1 = CenterOfMass(tri_1)
       c2 = CenterOfMass(tri_2)
       cosang = dot_product(c1, c2)
    End Select
    
    pi = 4.0 * atan(1.0)
    angle = 180.0 * Safe_Acos(cosang) / pi
  End Function Separation

  !=====================================================================
  ! FUNCTIONS FROM m_RegionIterator
  !=====================================================================

  ! Create a new iterator
  Subroutine NewIter(iter, level)
    Type(RegionIterator), Intent(Out) :: iter
    Integer, Intent(In) :: level
    
    Allocate(iter%refinement_path(base_refinement:level), STAT=err)
    If (err /= 0) Then
      Print *, "Error: Failed to allocate iterator memory"
      Return
    End If
    
    iter%refinement_path(base_refinement) = 0
    iter%index = DONE
  End Subroutine NewIter
  
  ! Reset an iterator
  Function reset_iter(iter, level, index) Result(idx)
    Type(RegionIterator), Intent(InOut) :: iter
    Integer, Intent(In), Optional :: level
    Integer, Intent(In), Optional :: index
    Integer :: idx
    Integer :: level_
    Integer :: i, n
    
    level_ = base_refinement
    If (Present(level)) level_ = level
    
    n = Ubound(iter%refinement_path, 1)
    
    iter%refinement_path = 0
    
    If (Present(index)) Then
       idx = index
       Do i = base_refinement, Ubound(iter%refinement_path, 1)
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
    idx = iter%index
  End Function reset_iter
  
  ! Get next region in iterator
  Function next(iter, level) Result(idx)
    Type(RegionIterator), Intent(InOut) :: iter
    Integer, Optional, Intent(In) :: level
    Integer :: idx
    
    Integer :: level_, n, i
    
    n = Ubound(iter%refinement_path, 1)
    
    If (Present(level)) Then
       level_ = level
    Else
       level_ = n
       Do i = base_refinement, n
          If (iter%refinement_path(i) == 0) Then
             level_ = i
             Exit
          End If
       End Do
    End If
    
    iter%refinement_path(level_+1:) = 0
    iter%refinement_path(:level_-1) = max(1, iter%refinement_path(:level_-1))
    
    Do 
       iter%refinement_path(level_) = iter%refinement_path(level_) + 1
       
       Select Case(level_)
       Case(base_refinement+1:)
          If (iter%refinement_path(level_) <= n_refine) Exit
       Case(base_refinement)
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
    
    idx = iter%index
  End Function next
  
  ! Convert refinement path to region index
  Function RefinementToRegionIndex(iter) Result(idx)
    Type(RegionIterator), Intent(In) :: iter
    Integer :: idx
    
    Integer :: i, j, n, idx_
    
    n = Ubound(iter%refinement_path, 1)
    
    idx_ = 0
    
    Do i = base_refinement, n
       j = iter%refinement_path(i)
       If (j == 0) Exit
       idx_ = idx_ + 1 + coef(n-i-1)*(j-1)
    End Do
    
    idx = idx_
  End Function RefinementToRegionIndex

  !=====================================================================
  ! FUNCTIONS FROM m_Icosahedron
  !=====================================================================

  ! Initialize the icosahedron
  Subroutine init_icosahedron()
    Real :: dphi, phi, lat
    Integer :: i
    
    If (icosahedron_initialized) Return
    
    ! Set north and south poles
    vertices(:,1) = (/ 0.0, 0.0, +1.0 /)  ! North pole
    vertices(:,12) = (/ 0.0, 0.0, -1.0 /) ! South pole
    
    ! Set remaining vertices
    dphi = (4.0*atan(1.0))/5.0  ! pi/5
    phi = 0.0
    lat = atan(0.5) 
    
    Do i = 1, 5
       vertices(:,6+i) = (/ cos(lat)*cos(phi), cos(lat)*sin(phi), -sin(lat) /)
       phi = phi + dphi
       vertices(:,1+i) = (/ cos(lat)*cos(phi), cos(lat)*sin(phi), +sin(lat) /)
       phi = phi + dphi
    End Do
    
    ! Set up octahedral axes and refinement vectors needed for xyz2reg
    octahedral_axis_1 = UnitVector(Cross(vertices(:,9), vertices(:,3)))
    octahedral_axis_2 = UnitVector(Cross(vertices(:,1), vertices(:,2)))
    octahedral_axis_3 = UnitVector(Cross(vertices(:,7), vertices(:,8)))
    
    refine_a = Cross(vertices(:,2), vertices(:,3))
    refine_b = Cross(vertices(:,3), vertices(:,8))
    refine_c = Cross(vertices(:,8), vertices(:,2))
    
    icosahedron_initialized = .true.
  End Subroutine init_icosahedron
  
  ! Get triangular faces of the icosahedron
  Function GetFaceTriangles() Result(facelist)
    Type(Spherical_Triangle) :: facelist(n_faces)
    Integer :: face
    
    Do face = 1, n_faces
       facelist(face) = Construct_Triangle(&
         vertices(:,face_vertices(1,face)),&
         vertices(:,face_vertices(2,face)),&
         vertices(:,face_vertices(3,face)))
    End Do
  End Function GetFaceTriangles

  ! Determine which icosahedral face contains a point
  Function icos_xyz2reg(p) Result(region)
    Real, Intent(In) :: p(3)
    Integer :: region
    
    Real :: c(3), pp(3)
    Real :: test_a, test_b, test_c
    Integer :: reg
    
    ! This table specifies which icosahedron faces are associated
    ! with which octahedral face
    Integer, Parameter :: lookup(4,8) = Reshape(Source = (/ &
       &         2, 12, 11,  7,            2,  8,  4,  3, &
       &         1, 10,  4,  5,            1, 15, 11,  6, &
       &        19, 15, 16, 20,           19, 10,  9, 14, &
       &        18,  8,  9, 13,           18, 12, 16, 17 /), &        
       & Shape = (/ 4, 8 /))
    
    ! Determine which octant the point lies in
    reg = 1
    c(1) = dot_product(p, octahedral_axis_1)
    c(2) = dot_product(p, octahedral_axis_2)
    c(3) = dot_product(p, octahedral_axis_3)
    
    If (c(3) < 0) Then ! "Southern" hemisphere
       reg = reg + 4   
       ! rotate 180 degrees about "x" axis
       c(2) = -c(2)
       c(3) = -c(3)
    End If
    
    If (c(2) < 0) Then ! "Western" hemisphere
       reg = reg + 2
       ! rotate 180 degrees about "north pole"
       c(1) = -c(1)
       c(2) = -c(2)
    End If
    
    If (c(1) < 0) Then ! Reflect about x axis
       reg = reg + 1
    End If
    
    ! Project the point p back onto the primary octant
    pp = c(2) * octahedral_axis_2 + c(3) * octahedral_axis_3 + abs(c(1)) * octahedral_axis_1
    
    ! Refine octant onto icosahedral regions
    test_a = dot_product(pp, refine_a) ! >0 in corner A
    test_b = dot_product(pp, refine_b) ! >0 in corner B
    
    If ((test_a >= 0) .or. (test_b >= 0)) Then
       If (test_a >= 0) Then 
          region = lookup(1, reg)
       Else ! (test_b >= 0)
          region = lookup(2, reg)
       End If
    Else
       test_c = dot_product(pp, refine_c) ! >0 in corner C
       If (test_c > 0) Then  
          region = lookup(3, reg)
       Else ! center region
          region = lookup(4, reg)
       End If
    End If
  End Function icos_xyz2reg

  !=====================================================================
  ! MAIN FUNCTIONS FROM m_Spherical_Partition
  !=====================================================================

  ! Initialize a partition of the unit sphere
  Subroutine Initialize(n_levels, partition, compress)
    Integer, Intent(In) :: n_levels ! refinement level
    Type(Spherical_Partition), Optional, Target :: partition
    Logical, Intent(In), Optional :: compress

    Type(Spherical_Partition), Pointer :: p_
    Integer :: n_coarse_regions, n_regions_all, n_fine_regions
    Integer :: idx, idxp, i, leaf_idx, level
    Type(RegionIterator) :: iter

    ! Check validity of refinement level
    If (n_levels < base_refinement) Then
      Print *, "Error: n_levels must be >= base_refinement"
      Return
    End If

    p_ => DefaultPartition(partition)
    p_%refinement_level = n_levels
    p_%compress = .true.
    If (Present(compress)) p_%compress = compress

    ! Initialize icosahedron
    Call init_icosahedron()

    ! Calculate number of regions at different levels
    n_fine_regions = n_refine**(n_levels+1) * n_faces
    n_coarse_regions = (n_fine_regions - n_faces)/3
    n_regions_all = n_fine_regions + n_coarse_regions

    ! Allocate memory for regions and leaves
    Allocate(p_%region(n_regions_all), p_%leaf(n_fine_regions), STAT=err)
    If (err /= 0) Then
      Print *, "Error: Memory allocation failed"
      Return
    End If

    ! Use icosahedron to fill coarsest partition
    Call NewIter(iter, n_levels)
    idx = reset_iter(iter, level=base_refinement)
    
    ! Set up base indices
    Do i = 1, n_faces
       p_%base_indices(i) = idx
       idx = next(iter, level=base_refinement)
       If (idx == DONE) Exit
    End Do
    
    ! Set base triangles
    p_%region(p_%base_indices)%triangle = GetFaceTriangles()

    ! For single-level case
    If (n_levels == base_refinement) Then
       ! leaves are the regions
       p_%leaf = p_%base_indices
    End If

    ! Loop through levels and refine triangles
    leaf_idx = 0
    Do level = base_refinement, n_levels - 1 
       idx = reset_iter(iter, level=level)

       Do ! loop over coarse regions
          ! For each coarse region, get its children
          Do i = 1, n_refine
             idxp = next(iter, level=level+1) ! next fine region
             p_%region(idx)%child_index(i) = idxp
          End Do

          ! Refine triangles
          p_%region(p_%region(idx)%child_index)%triangle = Refine(p_%region(idx)%triangle, &
               tests=p_%region(idx)%test)
          p_%region(p_%region(idx)%child_index)%parent_index = idx

          ! If at the finest level, set leaf indices
          If (level == n_levels-1) Then
             p_%leaf(leaf_idx+1:leaf_idx+n_refine) = p_%region(idx)%child_index(:)
             leaf_idx = leaf_idx + n_refine
          End If

          idx = next(iter, level) ! next coarse region
          If (idx == DONE) Exit ! Done
       End Do
    End Do
  End Subroutine Initialize
  
   ! Clean up a partition and free associated memory
  Subroutine Clean(partition)
    Type(Spherical_Partition), Intent(InOut), Target, Optional :: partition
    
    Type(Spherical_Partition), Pointer :: p_
    
    p_ => DefaultPartition(partition)
    
    ! Deallocate memory
    If (Associated(p_%region)) Then
      Deallocate(p_%region, STAT=err)
      If (err /= 0) Then
        Print *, "Warning: Error deallocating region array"
      End If
      Nullify(p_%region)
    End If
    
    If (Associated(p_%leaf)) Then
      Deallocate(p_%leaf, STAT=err)
      If (err /= 0) Then
        Print *, "Warning: Error deallocating leaf array"
      End If
      Nullify(p_%leaf)
    End If
    
    ! Reset the partition to uninitialized state
    p_%refinement_level = INVALID
    p_%base_indices = INVALID
    p_%compress = .false.
  End Subroutine Clean
  
  ! Return total number of regions or regions at a specific level
  Function NumberOfRegions(partition, atlevel) Result(n)
    Type(Spherical_Partition), Intent(In), Optional, Target :: partition
    Integer, Optional, Intent(In) :: atlevel
    Integer :: n

    Type(Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)

    If (Present(atlevel)) Then
      n = n_faces * (n_refine ** (atlevel - base_refinement))
      Return
    End If

    If (p_%compress) Then
       n = n_faces * (n_refine ** (p_%refinement_level - base_refinement))
    Else
       n = size(p_%region)
    End If
  End Function NumberOfRegions

  ! Helper function to get the partition pointer
  Function DefaultPartition(partition) Result(p)
    Type(Spherical_Partition), Target, Intent(In), Optional :: partition
    Type(Spherical_Partition), Pointer :: p

    If (Present(partition)) Then
       p => partition
    Else
       p => private_partition
    End If
  End Function DefaultPartition
  
  ! Function to bin points into regions
  Subroutine xyz2reg(n, x, y, z, idx, partition, atlevel)
    Integer, Intent(In) :: n
    Real*8, Intent(In) :: x(n), y(n), z(n)
    Integer, Intent(Out) :: idx(n)
    Type(Spherical_Partition), Target, Intent(In), Optional :: partition
    Integer, Optional, Intent(In) :: atlevel

    Type(Spherical_Partition), Pointer :: part_
    Integer :: i, level, ireg, f, idx_old, iregp
    Real :: p(3), x_old, y_old, z_old, test
    Integer :: atlevel_

    part_ => DefaultPartition(partition)
    
    ! Check if partition is properly initialized
    If (part_%refinement_level < base_refinement) Then
      Print *, "Error: Partition not initialized"
      Return
    End If

    ! Set refinement level
    atlevel_ = part_%refinement_level
    If (Present(atlevel)) Then
       If (atlevel > atlevel_) Then
          Print *, "Error: Insufficient refinement for requested level"
          Return
       End If
       atlevel_ = atlevel
    End If

    x_old = -2. ! impossible value to initialize
    y_old = -2.
    z_old = -2.
    idx_old = INVALID

    Do i = 1, n
       ! Check if this point is at the same location as previous
       If ((x(i)-x_old)**2 + (y(i)-y_old)**2 + (z(i)-z_old)**2 == 0) Then
          idx(i) = idx_old
          Cycle
       End If

       ! Construct a vector from the coordinates
       p = (/ x(i), y(i), z(i) /)

       ! Determine which icosahedral region the point belongs to
       f = icos_xyz2reg(p)

       ! Convert to the numbering scheme used by the partition
       ireg = 1 + (f-1) * coef(part_%refinement_level)
       iregp = f - 1

       ! Refine to the appropriate level
       Do level = base_refinement, atlevel_ - 1
          ! Test which subregion contains the point
          test = dot_product(p, part_%region(ireg)%test(:,1))
          If (test >= 0) Then
             ireg = part_%region(ireg)%child_index(1)
             iregp = n_refine*iregp + 0
          Else
             test = dot_product(p, part_%region(ireg)%test(:,2))
             If (test >= 0) Then
                ireg = part_%region(ireg)%child_index(2)
                iregp = n_refine*iregp + 1
             Else
                test = dot_product(p, part_%region(ireg)%test(:,3))
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

       ! Return the appropriate index
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

  Logical Function InitializedQ(partition)
    Type (Spherical_Partition), Intent(In), Target, Optional :: partition
    Type (Spherical_Partition), Pointer :: p_

    p_ => DefaultPartition(partition)
    InitializedQ = (p_%refinement_level >= base_refinement)
  End Function InitializedQ

  ! Return the triangle of the specified region index
  Function GetRegion(idx, partition) Result(tri)
    Integer, Intent(In) :: idx
    Type(Spherical_Partition), Target, Intent(In), Optional :: partition
    Type(Spherical_Triangle) :: tri
    
    Type(Spherical_Partition), Pointer :: p_
    Integer :: idxp
    Logical :: error_detected = .false.

    ! Initialize triangle to prevent returning uninitialized values
    tri%A = (/0.0, 0.0, 0.0/)
    tri%B = (/0.0, 0.0, 0.0/)
    tri%C = (/0.0, 0.0, 0.0/)

    p_ => DefaultPartition(partition)
     
    ! Check if partition is properly initialized (replacing ASSERT)
    If (.not. InitializedQ(p_)) Then
       Print *, "Error in GetRegion: Partition not properly initialized"
       error_detected = .true.
    End If
  
    ! Check if index is valid (replacing ASSERT_NOMSG)
    If (.not. error_detected) Then
       If (.not. ((idx > 0) .and. idx <= NumberOfRegions(partition))) Then
          Print *, "Error in GetRegion: Invalid region index:", idx
          error_detected = .true.
       End If
    End If
  
  ! Only proceed if no errors were detected
  If (.not. error_detected) Then
    If (p_%compress) Then
      
      idxp = p_%leaf(idx)
      
      tri = p_%region(idxp)%triangle
    else
      tri = p_%region(idx)%triangle
    end If
  end If

  End Function GetRegion

End Module m_Spherical_Partition_Simplified
