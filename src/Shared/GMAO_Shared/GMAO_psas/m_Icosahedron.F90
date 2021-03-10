!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Icosahedron
!
! !DESCRIPTION:
!
!   This module implements a singleton class of an icosahedron inscribed in the
!   unit sphere (otherwise known as the earth).   A few of the methods are necessary
!   for non-intuitive reasons related to the m_Spherical_Partition module which
!   sits on top of this module.
! 
!   None of the methods may be used without first using "Initialize"
!   but, for performance reasons, this is not verified by the methods.
!   (Big argument in favor of using ASSERT, but disabling it for production.)
!   
! !INTERFACE:
!

#include "assert.H"
Module m_Icosahedron
  Use m_die,only : assert_
  Implicit None
  Private

  Public :: Initialize

  Public :: GetFaces
  Public :: GetFaceTriangles
  Public :: GetVertices
  Public :: GetEdges
  Public :: xyz2reg
  Public :: AreaOfFace
  Public :: LengthOfEdge
  Public :: RotateToStandard
  Public :: Archetype

  Public :: n_faces
  Public :: n_edges
  Public :: n_vertices

  Integer, Parameter :: n_faces     = 20
  Integer, Parameter :: n_edges     = 30
  Integer, Parameter :: n_vertices  = 12
  Integer, Parameter :: PrimaryFace = 7 ! arbitrary (but unfortunate choice)

  ! This table labels the vertices of each face of the icosahedron.
  ! I.e. face #1 is defined by the vertices numbered (1,2,6) and so on.
  Integer, Parameter :: face_vertices(3,n_faces) =  &
       & RESHAPE(SOURCE = &
       & (/  1,  2,  6,     1,  3,  2,     1,  4,  3,     1,  5,  4, &
       &     1,  6,  5,     7,  6,  2,     8,  2,  3,     9,  3,  4, &
       &    10,  4,  5,    11,  5,  6,     2,  8,  7,     3,  9,  8, &
       &     4, 10,  9,     5, 11, 10,     6,  7, 11,    12,  7,  8, &
       &    12,  8,  9,    12,  9, 10,    12, 10, 11,    12, 11,  7 /), &
       & SHAPE = (/ 3, n_faces /))

  Real :: vertices(3,n_vertices)
  Integer, Parameter :: archetype(6) = (/ 1,  2,  3,  8,  13,  18 /)

  ! A flag to test if the module has been initialized
  Logical :: init = .false.

  ! The primary axes for the octahedron on which the icosahedron
  ! is based.  These are initialized in the "initialize routine"
  Real :: octahedral_axis_1(3)
  Real :: octahedral_axis_2(3)
  Real :: octahedral_axis_3(3)

  ! These vectors and matrices are used to refine the octahedron
  ! onto the icosahedron.  
  Real :: refine_a(3)
  Real :: refine_b(3)
  Real :: refine_c(3)
  Real, Dimension(3,3) :: rotate_a, rotate_ap
  Real, Dimension(3,3) :: rotate_b, rotate_bp
  Real, Dimension(3,3) :: rotate_c, rotate_cp

#ifndef F95
  Type Element
     Integer :: shift            ! 1-20
     Integer :: rotation         ! 0,1,2 clockwise
  End Type Element
#else
  Type Element
     Integer :: shift = -1       ! 1-20
     Integer :: rotation         ! 0,1,2 clockwise
  End Type Element
#endif

  Type (Element) :: translationtable(n_faces, n_faces)

  Character(Len=*), Parameter :: myname='m_Icosahedron'

  ! Use * for group operations
  Interface Operator(*)
     Module Procedure operate
  End Interface
  
  Real :: LengthOfEdge_ = -1
  Real :: AreaOfFace_ = -1
  Real :: pi

! !REVISION HISTORY:
!        1Dec00 - Tom CLune and Peter Lyster <lys@dao.gsfc.nasa.gov>
!                 Install Clune's new algorithms in support of
!                 advanced geometry and Icosahedral refinement scheme etc.
!                 The initial (birthday) files are in root molotov1.gsfc.nasa.gov:/home/jguo/Gcvs/ in 
!                 repository mppi under tag new_refinement
!EOP ___________________________________________________________________

Contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This module really implements a "singleton" class of icosahedra.
  ! I.e. since all icosahedra are identical, there is only one instance.
  ! This routine fills in some essential tables that would be parameters,
  ! if trig functions were allowed in Fortran initialization expressions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine Initialize()
    Use m_geometry, only : MakeRotationMatrix, Cross, UnitVector
    Implicit None
    Real :: v_main(3,3)
    Real :: pi, angle, anglep
    Character(Len=*), Parameter :: myname_=myname//'::Initialize'

    If (init) Return
    init = .true.

    ! get the coordinates of the vertices
    Call ConstructVertices()

    ! get the primary axes of the base octahedron

    octahedral_axis_1 = UnitVector(Cross(vertices(:,9), vertices(:,3)))
    octahedral_axis_2 = UnitVector(Cross(vertices(:,1), vertices(:,2)))
    octahedral_axis_3 = UnitVector(Cross(vertices(:,7), vertices(:,8)))

    pi = 4*atan(1.0)
    angle = 2*pi/5
    anglep = 4*pi/5

    v_main = vertices(:,face_vertices(:,PrimaryFace))

    refine_a = cross(vertices(:,2),vertices(:,3))
    refine_b = cross(vertices(:,3),vertices(:,8))
    refine_c = cross(vertices(:,8),vertices(:,2))

    Call FillTranslationTable()

    LengthOfEdge_ = Sqrt(2 - 2/Sqrt(5.))
    AreaOfFace_ = 4*pi/n_faces

  End Subroutine Initialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The translation table is used to enable symmetry considerations to reduce
  ! the size of separation tables constructed in m_Spherical_Partition.
  ! Rotations about face #1 and rotations about the north pole are the
  ! 2 generators used to construct this translator.
  Subroutine FillTranslationTable()
    Implicit None

    Integer :: face
    Type (Element), Parameter :: Identity(n_faces) =  &
	 (/ (Element(face, 0), face = 1, n_faces) /)

    ! Clockwise rotation about the north pole
    Type (Element), Parameter :: C5(n_faces) = (/     & 
	 &  (Element( 1+5*((face-1)/5) + Mod(face+3,5),0),face= 1,20)/)


    ! Clockwise rotation about face 1
    Type (Element), Parameter :: C3(n_faces) = (/    &
	 & Element( 1,1), Element( 6,2), Element(11,0), Element( 7,1), Element( 2,2), &
	 & Element( 5,2), Element(15,2), Element(16,0), Element(12,1), Element( 3,1), &
	 & Element(10,2), Element(20,2), Element(17,1), Element( 8,1), Element( 4,0), &
	 & Element(14,2), Element(19,1), Element(18,2), Element(13,1), Element( 9,0)  &
	 & /)
	 
    Character(Len=*), Parameter :: myname_=myname//'::FillTranslationTable'

    TranslationTable(:,1) = Identity

    Do face = 2, 5
       TranslationTable(:,face) = C5 * TranslationTable(:,face-1)
    End Do

    TranslationTable(:,6) = C3 * (TranslationTable(:,5) * C3)
    Do face = 7, 10
       TranslationTable(:,face) = TranslationTable(:,face-1) * C5
    End Do

    TranslationTable(:,11) = C5 * (C5 * (C3 * C3))
    Do face = 12, 15
       TranslationTable(:,face) = TranslationTable(:,face-1) * C5
    End Do

    TranslationTable(:,16) = C3 * (TranslationTable(:,11) * (TranslationTable(:,4) * C3))
    Do face = 17, 20
       TranslationTable(:,face) = TranslationTable(:,face-1) * C5
    End Do

#ifndef NDEBUG
    Do face = 1, n_faces
       ASSERT_NOMSG(TranslationTable(face,face)%shift == 1)
       ASSERT_NOMSG(TranslationTable(face,face)%rotation == 0)
    End Do
#endif

  End Subroutine FillTranslationTable

  Function operate(gen1, gen2) Result(gen)
    Implicit None
    Type (Element) :: gen(n_faces)
    Type (Element), Intent(In) :: gen1(n_faces), gen2(n_faces)
    
    Integer :: i
    Type (Element) :: e1, e2

    Do i = 1, n_faces
       e1 = gen2(i)
       e2 = gen1(e1%shift)
       gen(i) = Element(e2%shift, mod(e2%rotation + e1%rotation, 3))
     !	      &   (e1%parity_flip .xor. e2%parity_flip))
      End Do

    End Function operate


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! A trivial function, since by symmetry each face has the same area
  ! and hence is equal to 1/20th of the entire sphere (4 pi).
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Real Function AreaOfFace()
    Implicit None
    Real :: four_pi ! area of unit sphere
    Character(Len=*), Parameter :: myname_=myname//'::AreaOfFace'

    ASSERT(init)
    AreaOfFace = areaofface_

  End Function AreaOfFace

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Not a public method, but used for the initialization.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine ConstructVertices()
    Implicit None
    Character(Len=*), Parameter :: myname_=myname//'::ConstructVertices'

    Integer :: i
    Real :: dphi ! pi/10 is a common angle in this geometry
    Real :: phi, lat  ! used for useful intermediate angles

    vertices(:,1)  = (/ 0, 0, +1 /) ! The north pole
    vertices(:,12) = (/ 0, 0, -1 /) ! The south pole

    dphi = (4*atan(1.0))/5  ! pi/5
    phi  = 0
    lat = Atan(0.5)

    Do i = 1, 5
       vertices(:,6+i) = (/ cos(lat) * cos(phi), cos(lat) * sin(phi), - sin(lat) /)
       phi = phi + dphi
       vertices(:,1+i) = (/ cos(lat) * cos(phi), cos(lat) * sin(phi), + sin(lat) /)
       phi = phi + dphi
    End Do
    
  End Subroutine ConstructVertices

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! The public method which simply copies values from an internal table.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Function GetVertices() Result(v)
    Implicit None
    Real :: v(3,n_vertices)
    Character(Len=*), Parameter :: myname_=myname//'::GetVertices'

    ASSERT(init)
    v = vertices

  End Function GetVertices

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! GetFaces returns the indices of the vertices of the faces of the
  ! icosahedra.  Note that each vertex is listed 5 times, since it occurs
  ! on 5 different faces.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Function GetFaces() Result(facelist)
    Implicit None
    Integer facelist(3,n_faces)
    Character(Len=*), Parameter :: myname_=myname//'::GetFaces'
    
    ASSERT(init)
    facelist = face_vertices

  End Function GetFaces

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! GetFaceTriangles returns the coordinates of the vertices of the faces of the
  ! icosahedra.  Note that each vertex is listed 5 times, since it occurs
  ! on 5 different faces.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Function GetFaceTriangles() Result(facelist)
    Use m_Spherical_Triangle, only : Spherical_Triangle, Construct
    Implicit None
    Type (Spherical_Triangle) :: facelist(n_faces)

    Real :: v(3,n_vertices) 
    Integer :: face
    Character(Len=*), Parameter :: myname_=myname//'::GetFaceTriangles'
    
    ASSERT(init)
    v = GetVertices()
    Do face = 1, n_faces
       facelist(face) = Construct(v(:,face_vertices(1,face)),v(:,face_vertices(2,face)), &
	    & v(:,face_vertices(3,face)))
    End Do


  End Function GetFaceTriangles

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! GetEdges returns the coordinates of the vertices of the edges of the
  ! icosahedra.  Note that each vertex is listed 5 times, since it occurs
  ! on 5 different edges.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Function GetEdges() Result(edgelist)
    Implicit None
    Integer :: edgelist(2,n_edges)
    Integer :: face, n, i, ip
    Character(Len=*), Parameter :: myname_=myname//'::GetEdges'

    ! Rely on the ordering is such that each edge
    ! is listed once in ascending order and once is descending
    ! order, where order is based on vertex numbering scheme.
    ! This actually only requires that all faces follow the 
    ! right-hand rule in the ordering of vertices.

    ASSERT(init)
    n = 0
    Do face = 1, n_faces ! loop over faces
       Do i = 1, 3       ! loop over edges
	  ip = 1 + mod(3+(i-2),3)
	  If (face_vertices(i,face) > face_vertices(ip,face)) Then
	     n = n + 1
	     edgelist(1,n) = face_vertices(ip,face)
	     edgelist(2,n) = face_vertices(i,face)
	  End If
       End Do
    End Do

  End Function GetEdges

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Length of edge could actually be a parameter, if only trig functions
  ! coud be used for initialization expressions.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Real Function LengthOfEdge()
    Implicit None

    ASSERT(init)
    LengthOfEdge = lengthofedge_ ! basic geometry

  End Function LengthOfEdge

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Given a point on the unit sphere, return the face which contains that 
  ! point.  This is a critical method for binning data.  The
  ! point itself is also translated (by symmetry) to the equivalent point
  ! on the default "PrimaryFace", thereby minimizing any external tables
  ! used for further triangular refinements.
  ! 
  ! This algorithm is substantially different from the original.
  ! It uses octahedral symmetry for the first cut and then triangular
  ! refinement to bin onto the icosahedron.  This _slightly_
  ! reduces the number of logical tests that must be performed, and
  ! as of this version is probably a bit more transparent.  The previous
  ! version made heavy use of the m4 preprocessor, and was quite opaque.
  ! The method of translating the point p into a narrowing set of regions,
  ! enables clarity without sacrificing performance.  Of course that
  ! version went further and refined to 80 regions rather than 20.
  ! That functionality is now split, with the remainder lying in
  ! m_Spherical_Partition.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Function xyz2reg(p) Result(region)
    Use m_geometry, only : Cross
    Implicit None
    Integer :: region
    Real, Intent(In) :: p(3)

    Real :: c(3), pp(3), temp
    Real :: test_a, test_b, test_c
    Integer :: reg

    ! This table specifies which icosahedron faces are associated
    ! with which octahedral face.  It is inseparable from the
    ! choices of axes for the octahedra themselves.

    Integer, Parameter :: lookup(4,8) = Reshape(Source = (/ &
	 &         2, 12, 11,  7,            2,  8,  4,  3, &
	 &         1, 10,  4,  5,            1, 15, 11,  6, &
	 &        19, 15, 16, 20,           19, 10,  9, 14, &
	 &        18,  8,  9, 13,           18, 12, 16, 17 /), &        
	 & Shape = (/ 4, 8 /))

    Logical :: reflect
    Character(Len=*), Parameter :: myname_=myname//'::xyz2reg'

    ASSERT(init)

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
       reflect = .true.
    Else
       reflect = .false.
    End If

    ! Using the new values for "c" project the point p back
    ! onto the primary octant.  This allows a the same refinement tests
    ! to be applied independent of octant.

    pp = c(2) * octahedral_axis_2 + c(3) * octahedral_axis_3 + abs(c(1)) * octahedral_axis_1

    ! Now refine octant onto icosahedral regions.
    ! Use vector tests to see if in corner subregion:

    test_a = dot_product(pp, refine_a) ! >0 in corner A
    test_b = dot_product(pp, refine_b) ! >0 in corner B
    
    If ((test_a >= 0).or. (test_b >= 0)) Then
       If (test_a >= 0) Then 
	  region = lookup(1,reg)
	  ! rotate to the "center" region
       Else ! (test_b >= 0)
	  region = lookup(2,reg)
	  ! rotate to the "center" region
       End If
    Else
       test_c = dot_product(pp, refine_c) ! >0 in corner C
       If (test_c > 0 ) Then  
	  region = lookup(3,reg)
	  ! rotate to the "center" region
       Else ! center region
	  region = lookup(4,reg)
	  ! rotate to the "center" region, which is a null op in this case.
       End If
    End If

  End Function xyz2reg

  !
  ! This function rotates faces f1 and f2 such that f1
  ! is in the position (and orientation) of PrimaryFace
  ! and returns the new position and orientation of f2.

  Subroutine RotateToStandard(f1, f2, f2p, rotation_1, rotation_2, reflect)
    Implicit None
    Integer, Intent(In) :: f1, f2
    Integer, Intent(Out) :: f2p, rotation_1, rotation_2
    Logical, Intent(Out) :: reflect

    Type sym
       Integer :: face
       Integer :: rotation_1
       Integer :: rotation_2
       Logical :: reflect
    End Type sym
    Logical, Parameter :: T = .true. , F = .false.

    Type (sym), Parameter :: ToArchetype(n_faces) = (/ &
	 & sym( 1,0,0,F), sym( 2,0,0,F), sym( 3,0,0,F), sym( 3,0,0,T),  &
	 & sym( 2,0,0,T), sym( 2,2,1,F), sym( 3,2,2,T), sym( 8,0,0,F),  &
         & sym( 8,0,0,T), sym( 3,1,1,F), sym( 3,2,0,F), sym( 8,2,2,T),  &
	 & sym(13,0,0,F), sym( 8,1,1,F), sym( 3,1,0,T), sym( 8,2,0,F),  &
	 & sym(13,2,2,F), sym(18,0,0,F), sym(13,1,1,F), sym( 8,1,0,T) /)

    ASSERT(init)

    ASSERT((f1 > 0) .and. (f1 <= n_faces))
    ASSERT((f2 > 0) .and. (f2 <= n_faces))

    f2p = TranslationTable(f2,f1)%shift
    rotation_2 = TranslationTable(f2,f1)%rotation

    rotation_1 = ToArchetype(f2p)%rotation_1
    rotation_2 = Mod(ToArchetype(f2p)%rotation_2 + rotation_2,3)
    reflect = ToArchetype(f2p)%reflect
    f2p   = ToArchetype(f2p)%face

    ASSERT_NOMSG( (f2p > 0) .and. (f2p <= n_faces) )
    ASSERT_NOMSG( (rotation_1 >= 0) .and. (rotation_1 < 3) )
    ASSERT_NOMSG( (rotation_2 >= 0) .and. (rotation_2 < 3) )

  End Subroutine RotateToStandard
  
End Module m_Icosahedron
