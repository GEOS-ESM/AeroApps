!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Spherical_Triangle
!
! !DESCRIPTION:
!
!   This module implements a class of triangles on the unit sphere.
!   All of the methods are geometric (and in my opinion fairly transparent).
!   This module is primarily used by m_Spherical_Partition.
!   
! !INTERFACE:
!

#include "assert.H"
Module m_Spherical_Triangle
  Use m_die,only : assert_
  Implicit None
  Private

  Public :: Spherical_Triangle
  Public :: Construct
  Public :: GetVertices
  Public :: Circumcenter
  Public :: CenterOfMass
  Public :: Separation
  Public :: CosSeparation
  Public :: Refine
  Public :: Area

  Public :: SEPANG_MIN
  Public :: SEPANG_MAX
  Public :: SEPANG_CIRCUMCENTER
  Public :: SEPANG_CENTER_OF_MASS

  Public :: n_refine

#ifndef F95
  Type Spherical_Triangle
     Private
     SEQUENCE
     Real :: A(3) ! vertex A
     Real :: B(3) ! vertex B
     Real :: C(3) ! vertex C
  End Type Spherical_Triangle
#else
  Type Spherical_Triangle
     Private
     SEQUENCE
     Real :: A(3) = 0 ! vertex A
     Real :: B(3) = 0 ! vertex B
     Real :: C(3) = 0 ! vertex C
  End Type Spherical_Triangle
#endif
  Character(len=*), Parameter :: myname = 'm_Spherical_Triangle'

  ! Would be nice to use a "enum" type here ...
  Integer, Parameter :: SEPANG_MIN            = 1
  Integer, Parameter :: SEPANG_MAX            = 2
  Integer, Parameter :: SEPANG_CIRCUMCENTER   = 3
  Integer, Parameter :: SEPANG_CENTER_OF_MASS = 4
  Integer, Parameter :: sep_methods(4) = &
       &(/ SEPANG_MIN, SEPANG_MAX, SEPANG_CIRCUMCENTER, SEPANG_CENTER_OF_MASS /)

  ! The number of regions created by a refinement.
  Integer, Parameter :: n_refine = 4

  Interface Area
     Module Procedure Area_
  End Interface

! !REVISION HISTORY:
!        1Dec00 - Tom CLune and Peter Lyster <lys@dao.gsfc.nasa.gov>
!                 Install Clune's new algorithms in support of
!                 advanced geometry and Icosahedral refinement scheme etc.
!                 The initial (birthday) files are in root molotov1.gsfc.nasa.gov:/home/jguo/Gcvs/ in 
!                 repository mppi under tag new_refinement
!EOP ___________________________________________________________________
     
Contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Area of a spherical triangle (shudder)
  ! This is expensive, since it uses trigonometric functions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Function Area_(tri) Result (A)
    Use m_geometry, only : Cross, UnitVector, Safe_Acos
    Implicit None
    Type (Spherical_Triangle), Intent(In) :: tri
    Real :: A

    Real, Dimension(3)   :: U_ab, U_ac, U_ba, U_bc, U_ca, U_cb
    Real, Dimension(3,3) :: r
    Real :: pi
    Character(len=*), Parameter :: myname_ = myname//'::Area_scalar'

    ! U_ij is the unit vector at vertex i that is tangent to the 
    ! equator passing through arc ij.  Note that U_ij is perpendicular
    ! to r_i and to r_i X r_j 
    ! This uses twice as many "sqrt" calls as strictly necessary,
    ! but is clearer, and the "acos" dominates the work anyway.

    U_ab = UnitVector( cross(tri%A, cross(tri%B,tri%A)) )
    U_ac = UnitVector( cross(tri%A, cross(tri%C,tri%A)) )

    U_ba = UnitVector( cross(tri%B, cross(tri%A,tri%B)) )
    U_bc = UnitVector( cross(tri%B, cross(tri%C,tri%B)) )

    U_ca = UnitVector( cross(tri%C, cross(tri%A,tri%C)) )
    U_cb = UnitVector( cross(tri%C, cross(tri%B,tri%C)) )

    ! The area of a spherical triangle is 2 pi minus the sum of the 
    ! angles.

    pi = 4 * atan(1.0)

    A = safe_acos(dot_product(U_ab,U_ac)) + safe_acos(dot_product(U_ba,U_bc)) + &
	 & safe_acos(dot_product(U_ca,U_cb)) - pi

  End Function Area_
        

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The center of mass is simply the (normalized) average
  ! of the vertices.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Function CenterOfMass(tri) Result (center)
    Use m_geometry, only : UnitVector
    Implicit None
    Type (Spherical_Triangle), Intent(In) :: tri
    Real :: center(3)
    Character(len=*), Parameter :: myname_ = myname//'::CenterOfMass'


    center = tri%A + tri%B + tri%C
    center = UnitVector(center)

  End Function CenterOfMass

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function computes the circumcenter of a 
  ! spherical triangle. 
  Function Circumcenter(tri, cosrad) Result(center)
    Use m_geometry, only : UnitVector, cross
    Implicit None
    Type(Spherical_Triangle), Intent(In) :: tri
    Real, Optional :: cosrad   ! cosine of the radius of the circumcircle
    Real :: center(3)
    Character(len=*), Parameter :: myname_ = myname//'::Circumcenter'

    center = cross(tri%B,tri%A) + cross(tri%C,tri%B) + cross(tri%A,Tri%C)
    center = UnitVector(center)

    If (Present(cosrad)) cosrad = dot_product(tri%A,center)

  End Function Circumcenter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Used to create new triangles by specifying the vertices
  ! on the unit sphere.   (Perhaps the lentgths should be verified?)
  ! This would be  a trivial operation were it not for the private nature
  ! of the derived datatype.

  Function Construct(v1,v2,v3) Result (tri)
    Implicit None
    Real, Intent(In) :: v1(3), v2(3), v3(3)
    Type (Spherical_Triangle) :: tri
    Character(len=*), Parameter :: myname_ = myname//'::Construct'

    tri%A = v1
    tri%B = v2
    tri%C = v3

  End Function Construct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Returns the 4 daughter triangles based upon subdividing
  ! a triangle at its midpoints.

  Function Refine(tri, tests) Result (daughters)
    Use m_geometry, only : UnitVector, Cross
    Implicit None
    Type (Spherical_Triangle), Intent(In) :: tri
    Real, Optional, Intent(Out) :: tests(3,3)
    Type (Spherical_Triangle) :: daughters(4)

    Real :: AB(3),BC(3),CA(3)
    Real :: a(3),B(3),C(3)
    Character(len=*), Parameter :: myname_ = myname//'::Refine'


    ! center triangle
    A = tri%A
    B = tri%B
    C = tri%C
    AB = UnitVector(A+B)
    BC = UnitVector(B+C)
    CA = UnitVector(C+A)

    daughters(4) = Spherical_Triangle(BC,CA,AB)
    daughters(1) = Spherical_Triangle(A,AB,CA)
    daughters(2) = Spherical_Triangle(B,BC,AB)
    daughters(3) = Spherical_Triangle(C,CA,BC)

    If (Present(tests)) Then
       tests(:,1) = Cross(CA,AB)
       tests(:,2) = Cross(AB,BC)
       tests(:,3) = Cross(BC,CA)
    End If

  End Function Refine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Separation returns separation angle between two
  ! triangles by any of several measures.
  ! The maximum distance may be useful in higher level refinement procedures
  ! by indicating whethe further refinement will take a measure out of
  ! some separation window.
  ! Triangles are assumed to be non-intersecting, but may be nested.
  !
  ! The implemented metrics are:
  !
  ! SEPANG_MIN             minimum angle between any pair of points
  ! SEPANG_MAX             maximum angle between any pair of points
  ! SEPANG_CIRCUMCENTER    angle between circumcenters
  ! SEPANG_CENTER_OF_MASS  angle between COM
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Function Separation(tri_1, tri_2, method) Result(angle)
    Use m_geometry, only : safe_acos
    Implicit None
    Type (Spherical_Triangle), Intent(In) :: tri_1
    Type (Spherical_Triangle), Intent(In) :: tri_2
    Integer, Intent(In) :: method
    Real :: angle

    Real :: c1(3), c2(3)
    Real :: cosang, cosang_max, cosang_min
    Integer :: i, j
    Character(len=*), Parameter :: myname_ = myname//'::Separation'
    Real :: inside(3,2), pi

    ASSERT_MSG(Any(method==sep_methods), myname_)
    
    Select Case (method)
    Case (SEPANG_MIN)
       cosang_max = -1.0 ! worst case
       inside = 0
       Do j = 1, 3
	  Do i = 1, 3

	     ! Compute (cos) distance between i_th vertex of tri_1 and
	     ! side (arc) j on tri_2
	     cosang = arc_sep_min(vertex(tri_1,i),arc(tri_2,j),inside(j,1))
	     cosang_max = max(cosang_max, cosang)

	     ! Same with triangles swapped
	     cosang = arc_sep_min(vertex(tri_2,i),arc(tri_1,j),inside(j,2))
	     cosang_max = max(cosang_max, cosang)

	  End Do
       End Do
       ! Check to see if triangles are nested, in which case the separation
       ! is zero. (cosang = 1)
       If (Any(All(inside<0,1))) cosang_max = +1

       cosang = cosang_max

    Case (SEPANG_MAX)

       cosang_min = +1.0 ! worst case
       inside = 0
       Do j = 1, 3
	  Do i = 1, 3

	     ! Compute (cos) distance between i_th vertex of tri_1 and
	     ! side (arc) j on tri_2
	     cosang = arc_sep_max(vertex(tri_1,i),arc(tri_2,j),inside(j,1))
	     cosang_min = min(cosang_min, cosang)

	     ! Same with triangles swapped
	     cosang = arc_sep_max(vertex(tri_2,i),arc(tri_1,j),inside(j,2))
	     cosang_min = min(cosang_min, cosang)

	  End Do
       End Do

       ! Check to see if triangles are nested, in which case the separation
       ! is zero. (cosang = 1)
       If (Any(All(inside<0,1))) cosang_min = +1
       cosang = cosang_min

    Case (SEPANG_CIRCUMCENTER)

       c1 = Circumcenter(tri_1)
       c2 = Circumcenter(tri_2)
       cosang = dot_product(c1,c2)

    Case (SEPANG_CENTER_OF_MASS)

       c1 = CenterOfMass(tri_1)
       c2 = CenterOfMass(tri_2)
       cosang = dot_product(c1,c2)

    End Select

    pi = 4*atan(1.)
    angle = 180*safe_acos(cosang)/pi

  End Function Separation

!---------------------------------------------------------------------------
! If the regions A and B are allowed to be nested, but not intersecting,
! we have one of four possibilities:
!
!   I. Every vertex of A is outside of region B, and every vertex of B is
!      outside of region A.
!  II. Every vertex of A is inside of region B, and every vertex of B is
!      outside of region A.
! III. Every vertex of A is outside of region B, and every vertex of B is
!      inside of region A.
!  IV. Every vertex of A is inside of region B, and every vertex of B is
!      inside of region A.
! 
! (We refer to points on the boundary of a region as being both inside and
! outside.)
! 
! For case IV, region A and B are identical, and the arc based formulae
! are correct.
! 
! For case I the arc based formulae are correct.
! 
! Fore cases II and III, the arc-based formulae are clearly incorrect,
! and the sepang min and max should be zero.  One must merely detect
! these cases.  A simple test can determine if a point is on the
! "correct" side of the great-circle of the other region.  This test
! produces a real quantity (the variable "inside") whose sign specifies
! inside (+) or outside (-).  We can sum over all the vertices of a
! region, because there sign is guaranteed to be constant for
! non-intersecting triangles.  (I.e. all three vertices of a region are
! either inside or all three vertices are outside.)  If the sum is
! positive for all arcs of the partner region, then the triangle is
! nested. We thus form 3 sums for Case II and 3 sums for Case III,
! and set cos_sepang = 1 if either set is all poistive.
! 
! The logic contained here fails for intersecting non-nested triangles.
!-----------------------------------------------------------------------

  Function CosSeparation(tri_1, tri_2, method, threshold) Result(test)
    Implicit None
    Type (Spherical_Triangle), Intent(In) :: tri_1
    Type (Spherical_Triangle), Intent(In) :: tri_2
    Integer, Intent(In) :: method
    Real, Intent(In) :: threshold
    Logical :: test

    Real :: angle
    Real :: inside(3,2)
    Real :: c1(3), c2(3)
    Real :: cosang, cosang_max, cosang_min
    Integer :: i, j
    Character(len=*), Parameter :: myname_ = myname//'::CosSeparation'

    ASSERT_MSG(Any(method==(/SEPANG_MIN,SEPANG_MAX/)), myname_)
    test = .false.

    Select Case (method)
    Case (SEPANG_MIN)
       cosang_max = -1.0 ! worst case
       inside = 0
       Do j = 1, 3
	  Do i = 1, 3

	     ! Compute (cos) distance between i_th vertex of tri_1 and
	     ! side (arc) j on tri_2
	     cosang = arc_sep_min(vertex(tri_1,i),arc(tri_2,j),inside(j,1))
	     cosang_max = max(cosang_max, cosang)
	     If (cosang_max > threshold) Then
		test = .true.
		Exit
	     End If

	     ! Same with triangles swapped
	     cosang = arc_sep_min(vertex(tri_2,i),arc(tri_1,j),inside(j,2))
	     cosang_max = max(cosang_max, cosang)
	     If (cosang_max > threshold) Then
		test = .true.
		Exit
	     End If

	  End Do
	  If (test) exit
       End Do

       If (.not. test) Then
	  ! Check to see if triangles are nested, in which case the separation
	  ! is zero. (cosang = 1) and the test is false.
	  If (Any(All(inside<0,1))) test = .true.
       End If

    Case (SEPANG_MAX)

       cosang_min = +1.0 ! worst case

       inside = 0
       Do j = 1, 3
	  Do i = 1, 3

	     ! Compute (cos) distance between i_th vertex of tri_1 and
	     ! side (arc) j on tri_2
	     cosang = arc_sep_max(vertex(tri_1,i),arc(tri_2,j),inside(j,1))
	     cosang_min = min(cosang_min, cosang)

	     If (cosang_min < threshold) Then
		test = .true.
		Exit
	     End If

	     ! Same with triangles swapped
	     cosang = arc_sep_max(vertex(tri_2,i),arc(tri_1,j),inside(j,2))
	     cosang_min = min(cosang_min, cosang)
	     If (cosang_min < threshold) Then
		test = .true.
		Exit
	     End If

	  End Do
	  If (test) exit
       End Do

       If (.not. test) Then
	  ! Check to see if triangles are nested, in which case the separation
	  ! is zero. (cosang = 1) and the test is true.
	  If (Any(All(inside<0,1))) test = .false.
       End If

    Case (SEPANG_CIRCUMCENTER)


       c1 = Circumcenter(tri_1)
       c2 = Circumcenter(tri_2)
       cosang = dot_product(c1,c2)
       test = (cosang > threshold)

    Case (SEPANG_CENTER_OF_MASS)

       c1 = CenterOfMass(tri_1)
       c2 = CenterOfMass(tri_2)
       cosang = dot_product(c1,c2)
       test = (cosang > threshold)

    End Select


  End Function CosSeparation

  ! Compute the cosine of the minimimum separation between a point and
  ! an arc.  The minimum distance may be to either endpoint, or to
  ! the perpendicular path from the point to the arc.
  Function arc_sep_min(p,arc_ab,inside) Result(cosang)
    Use m_geometry, only : Cross, UnitVector
    Implicit None
    Real, Intent(In) :: p(3)    ! The point
    Real, Intent(In) :: arc_ab(3,2) ! The arc
    Real, Intent(InOut) :: inside
    Real :: cosang

    Real :: n(3)    ! The north pole for the arc AB
    Character(len=*), Parameter :: myname_ = myname//'::arc_sep_min'
    
    n = cross(arc_ab(:,1),arc_ab(:,2))
    inside = inside + dot_product(p,n)
    
    If (dot_product(arc_ab(:,1),cross(p,n)) > 0) Then
       ! A is on the "right" side 
       If (dot_product(arc_ab(:,2),cross(p,n)) < 0) Then
	  ! B is on the "left" side and the
	  ! perpendicular distance to the arc is the shortest
	  cosang = Sqrt(max(0.,1 - dot_product(p,UnitVector(n))**2)) 
       Else ! B must be the closer point
	  cosang = dot_product(p,arc_ab(:,2)) 
       End If
    Else ! A is on the "left" side
       If (dot_product(arc_ab(:,2),cross(p,n)) < 0) Then
	  ! B is on the "right" side
	  cosang = dot_product(p,arc_ab(:,1)) ! A must be the closer point
       Else ! arc is on the far side of the sphere 
	  ! must check both endpoints
	  cosang = max(dot_product(p,arc_ab(:,1)), dot_product(p,arc_ab(:,2)))
       End If
    End If
    
    cosang = Min(cosang,+1.) ! to prevent roundoff causing illegal values
    
  End Function arc_sep_min
  
  ! Compute the cosine of the maximum separation between a point and
  ! an arc.  The maximum distance may be to either endpoint, or to
  ! the perpendicular path from the point to the arc.
  
  Function arc_sep_max(p,arc_ab,inside) Result(cosang)
    Use m_geometry, only : Cross, UnitVector
    Implicit None
    Real, Intent(In) :: p(3)    ! The point
    Real, Intent(In) :: arc_ab(3,2) ! The arc
    Real, Intent(InOut) :: inside
    Real :: cosang
    
    Real :: n(3)    ! The north pole for the arc AB
    Character(len=*), Parameter :: myname_ = myname//'::arc_sep_max'
    
    n = cross(arc_ab(:,1),arc_ab(:,2))
    inside = inside + dot_product(p,n)

    If (dot_product(arc_ab(:,1),cross(p,n)) > 0) Then
       ! A is on the "right" side 
       If (dot_product(arc_ab(:,2),cross(p,n)) < 0) Then
	  ! B is on the "left" side and the
	  ! the farthest distance is one of the endpoints
	  cosang = min(dot_product(p,arc_ab(:,1)), dot_product(p,arc_ab(:,2)))
       Else ! A must be the further point
	  cosang = dot_product(p,arc_ab(:,1)) 
       End If
    Else ! A is on the "left" side
       If (dot_product(arc_ab(:,2),cross(p,n)) < 0) Then
	  ! B is on the "right" side
	  cosang = dot_product(p,arc_ab(:,1)) ! B must be the furthest point arc
       Else ! arc is on the far side of the sphere 
	  ! perpendicular distance to the arc is the furthest
	  ! must check both endpoints
	  cosang = -Sqrt(max(0.,1 - dot_product(p,UnitVector(n))**2))
       End If
    End If
    
    cosang = Max(cosang,-1.) ! to prevent roundoff causing illegal values
    
  End Function arc_sep_max
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This trivial function selects an arc (specified by the opposing vertex)
  ! of a spherical triangle.
  
  Function arc(tri,i) Result (vv)
    Implicit None
    Type (Spherical_Triangle), Intent(In) :: tri
    Integer, Intent(In) :: i
    Real :: vv(3,2)
    Character(len=*), Parameter :: myname_ = myname//'::arc'
    
    ASSERT_MSG((Any(i==(/1,2,3/))),myname_)
    Select Case (i)
    Case (1)
       vv(:,1) = tri%B
       vv(:,2) = tri%C
    Case (2)
       vv(:,1) = tri%C
       vv(:,2) = tri%A
    Case (3)
       vv(:,1) = tri%A
       vv(:,2) = tri%B
    End Select
    
  End Function arc
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Trivial function to select a lone vertex from a triangle.
  ! This enables looping over vertices, but really points to
  ! the weakness in the implementation of the derived type.
  ! A singel 3x3 matrix would have been better, but more opaque
  ! in most of the methods.

  Function vertex(tri,i) Result (v)
    Implicit None
    Type (Spherical_Triangle), Intent(In) :: tri
    Integer, Intent(In) :: i
    Real :: v(3)
    Character(len=*), Parameter :: myname_ = myname//'::vertex'
    
    ASSERT_MSG((Any(i==(/1,2,3/))),myname_)
    Select Case (i)
    Case (1)
       v = tri%A
    Case (2)
       v = tri%B
    Case (3)
       v = tri%C
    End Select

  End Function vertex

  Function GetVertices(tri) Result(v)
    Implicit None
    Type (Spherical_Triangle), Intent(In) :: tri
    Real :: v(3,3)
    
    v(:,1) = tri%A
    v(:,2) = tri%B
    v(:,3) = tri%C

  End Function GetVertices

End Module m_Spherical_Triangle
