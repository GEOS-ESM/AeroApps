! F90 module encapsulating:
!     ll2xyz.F    02/03/2000
!     qtrig.f     02/03/2000
!     qtrig0.f    02/03/2000

Module m_spherical_geometry
  Implicit None
  Private 

  Public :: ll2xyz
  Public :: circumcenter
  Public :: normalize
  Public :: cross

  Interface circumcenter
     Module Procedure circumcenter_
  End Interface

  Interface normalize
     Module Procedure normalize_vec
     Module Procedure normalize_vec_list
  End Interface

  Interface cross
     Module Procedure cross_
     Module Procedure cyclic_cross_
  End Interface

!(trigtab.h)
  Integer, Parameter :: NN = 360
  Integer, Parameter :: lentab = 1 + 5*NN/4
  Real   , save :: trigtab(lentab)
  Logical, save :: initialized_=.false.
!.

  character(len=*),parameter :: myname='m_spherical_geometry'

Contains

!  03Feb2000 - Tom C. - converted to f90 free-format
!  02Feb95  - Jing G.	- Changed CRAY to _UNICOS for consistency and
!			  to follow the guide lines.
!  modification for dynamic storage on CRAY                     05/28/93
!  03oct94 - Implemented CRAY specifics with IFDEFs (A. da S.)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!.... Compute (x,y,z) coordinates from (rlons,rlats).

  subroutine ll2xyz( rlons, rlats, length, x, y, z, ierr )
    
    use m_stdio
    use m_mall,only : mall_ison,mall_mci,mall_mco
!.......................................................................
!.... Argument declarations.
    
    integer ::   length
    real    ::   rlons(length)
    real    ::   rlats(length)
    real    ::   x(length)
    real    ::   y(length)
    real    ::   z(length)
    integer ::   ierr


!.......................................................................
!.... Local storage.
    
    Real, allocatable :: coslat(:)
    
    Character(len=6), Parameter :: myname = 'll2xyz'


!.......................................................................
!.... Statement function.
    
    logical :: llbad
    Integer :: n
    llbad(n) = (abs(rlons(n)) > 180.0).or.(abs(rlats(n)) > 90.0)
    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   ..Consistency check
    
    allocate(coslat(length),stat=ierr)
    if(ierr.ne.0) then
       write(stderr,'(2a,i5)') myname, &
       &		': allocate() error, stat =',ierr
       return
    endif
    if(mall_ison()) call mall_mci(coslat,myname)
    
    ierr = 0

!.......................................................................
!.... Loop over the points.


    do n = 1, length
       
!.......................................................................
!....... Check that longitude and latitude are in range.
       
       if( llbad(n) ) then
	  print *, '  ll2xyz:  latitude or longitude out of range'
	  print *, '  ll2xyz:  n = ', n, '  rlats(n) = ', rlats(n), &
	       &   '  rlons(n) = ', rlons(n)
	  ierr = 1
	  if(mall_ison()) call mall_mco(coslat,myname)
	  deallocate(coslat)
	  return
       endif
       
    end do

!.......................................................................
!.... Compute sines and cosines of longitude and latitude.
       
    call qtrig( length, rlons, x, y )
    call qtrig( length, rlats, coslat, z )
       
!.......................................................................
!.... Compute x and y coordinates, ( z = sin(lat) ).
       
    do n = 1, length
       x(n) = coslat(n)*x(n)
       y(n) = coslat(n)*y(n)
    end do
	  
	  if(mall_ison()) call mall_mco(coslat,myname)
    deallocate(coslat)
    return
	  
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  end subroutine ll2xyz



!!$  subroutine qtrig( n, dangs, qcos, qsin )
!!$  Quick approximation to cosine and sine functions.

  subroutine qtrig( n, dangs, qcos, qsin )

!!$  trigtab    < input, real vector dimensioned trigtab(lentab) >
!!$          Tabulated array of sines for every one degree.
!!$          trigtab(i) must be initialized with Sin((i-181)*pi/180).
!!$          Note that trigtab can be initialized using routine Qtrig0.

!!$  n       < input, integer >
!!$          Number of cosines and sines to be approximated. Used
!!$          to dimension the vectors dangs, qcos, and qsin.

!!$  dangs   < input, real vector dimensioned n >
!!$          Angles for which cosines and sines are to be computed.
!!$          Units are degrees, and (-180.0.le.dangs(i).lt.180.0)
!!$          must be satisfied.

!!$  qcos    < output, real vector dimensioned n >
!!$          On return qcos(i) approximates cos(dangs(i)*pi/180).

!!$  qsin    < output, real vector dimensioned n >
!!$          On return qsin(i) approximates sin(dangs(i)*pi/180).

!!$  Jim Pfaendtner, July 26, 1992.

    Real, Parameter :: pi     = 3.141592654
    Real, Parameter :: pid180 = pi / 180.0 
    Real, Parameter :: r181p5 = NN/2 + 1.5 
    Real, Parameter :: conj   = NN / 360.0 
    Real, Parameter :: qconj  = 1.0 / conj 
    Real, Parameter :: r180   = 180.0      
    Integer, Parameter :: n181   = 181        
    Real, Parameter :: one    = 1.0        
    Real, Parameter :: half   = 0.5        
    Integer, Parameter :: n90    = NN/4       

    Integer :: n
    Real :: dangs(n), qcos(n), qsin(n)
    Integer :: i, j
    Real :: del, sdel, cdel

    if(.not.initialized_) call qtrig0()

    if( NN.eq.360 ) then
       do i = 1, n
	  j         = ifix( r181p5 + dangs(i) )
	  del       = dangs(i) - float( j - n181 )
	  sdel      = pid180 * del
	  cdel      = one - half * sdel**2
	  qcos(i)   = trigtab(j+n90) * cdel - trigtab(    j) * sdel
	  qsin(i)   = trigtab(j    ) * cdel + trigtab(j+n90) * sdel
       end do
    else
       do i = 1, n
	  j         = ifix( r181p5 + conj * dangs(i) )
	  del       = dangs(i) - qconj* float( j - 1 ) + r180
	  sdel      = pid180 * del
	  cdel      = one - half * sdel**2
	  qcos(i)   = trigtab(j+n90) * cdel - trigtab(    j) * sdel
	  qsin(i)   = trigtab(j    ) * cdel + trigtab(j+n90) * sdel
       end do
    endif

    return
  end subroutine qtrig

  subroutine qtrig0()
    use m_die,only : die

!!$  Fills sine table for qtrig.

!!$  trigtab    < output, real vector dimensioned trigtab(lentab) >
!!$          Tabulated array of sines for use by qtrig.
!!$          trigtab(i) is set to Sin((2*(i-1)-NN)*(pi/NN)) for
!!$          i = 1,2,3,...,lentab.

!!$     Jim Pfaendtner, July 26, 1992.

    Real, Parameter :: c00=1.000000000,     s00=0.000000000 
    Real, Parameter :: c01=0.999847695,     s01=0.017452406 
    Real, Parameter :: c02=0.999390827,     s02=0.034899497 
    Real, Parameter :: c03=c02*c01-s02*s01, s03=s02*c01+c02*s01 
    Real, Parameter :: c04=0.997564050,     s04=0.069756474 
    Real, Parameter :: c05=c04*c01-s04*s01, s05=s04*c01+c04*s01 
    Real, Parameter :: c06=c04*c02-s04*s02, s06=s04*c02+c04*s02 
    Real, Parameter :: c07=c04*c03-s04*s03, s07=s04*c03+c04*s03 
    Real, Parameter :: c08=0.990268069,     s08=0.139173101 
    Real, Parameter :: c09=c08*c01-s08*s01, s09=s08*c01+c08*s01 
    Real, Parameter :: c10=c08*c02-s08*s02, s10=s08*c02+c08*s02 
    Real, Parameter :: c11=c08*c03-s08*s03, s11=s08*c03+c08*s03 
    Real, Parameter :: c12=c08*c04-s08*s04, s12=s08*c04+c08*s04 
    Real, Parameter :: c13=c08*c05-s08*s05, s13=s08*c05+c08*s05 
    Real, Parameter :: c14=c08*c06-s08*s06, s14=s08*c06+c08*s06 
    Real, Parameter :: c15=c08*c07-s08*s07, s15=s08*c07+c08*s07 
    Real, Parameter :: c16=0.961261696,     s16=0.275637356 
    Real, Parameter :: c17=c16*c01-s16*s01, s17=s16*c01+c16*s01 
    Real, Parameter :: c18=c16*c02-s16*s02, s18=s16*c02+c16*s02 
    Real, Parameter :: c19=c16*c03-s16*s03, s19=s16*c03+c16*s03 
    Real, Parameter :: c20=c16*c04-s16*s04, s20=s16*c04+c16*s04 
    Real, Parameter :: c21=c16*c05-s16*s05, s21=s16*c05+c16*s05 
    Real, Parameter :: c22=c16*c06-s16*s06, s22=s16*c06+c16*s06 
    Real, Parameter :: c23=c16*c07-s16*s07, s23=s16*c07+c16*s07 
    Real, Parameter :: c24=c16*c08-s16*s08, s24=s16*c08+c16*s08 
    Real, Parameter :: c25=c16*c09-s16*s09, s25=s16*c09+c16*s09 
    Real, Parameter :: c26=c16*c10-s16*s10, s26=s16*c10+c16*s10 
    Real, Parameter :: c27=c16*c11-s16*s11, s27=s16*c11+c16*s11 
    Real, Parameter :: c28=c16*c12-s16*s12, s28=s16*c12+c16*s12 
    Real, Parameter :: c29=c16*c13-s16*s13, s29=s16*c13+c16*s13 
    Real, Parameter :: c30=c16*c14-s16*s14, s30=s16*c14+c16*s14 
    Real, Parameter :: c31=c16*c15-s16*s15, s31=s16*c15+c16*s15 
    Real, Parameter :: c32=0.848048096,     s32=0.529919264 
    Real, Parameter :: c33=c32*c01-s32*s01, s33=s32*c01+c32*s01 
    Real, Parameter :: c34=c32*c02-s32*s02, s34=s32*c02+c32*s02 
    Real, Parameter :: c35=c32*c03-s32*s03, s35=s32*c03+c32*s03 
    Real, Parameter :: c36=c32*c04-s32*s04, s36=s32*c04+c32*s04 
    Real, Parameter :: c37=c32*c05-s32*s05, s37=s32*c05+c32*s05 
    Real, Parameter :: c38=c32*c06-s32*s06, s38=s32*c06+c32*s06 
    Real, Parameter :: c39=c32*c07-s32*s07, s39=s32*c07+c32*s07 
    Real, Parameter :: c40=c32*c08-s32*s08, s40=s32*c08+c32*s08 
    Real, Parameter :: c41=c32*c09-s32*s09, s41=s32*c09+c32*s09 
    Real, Parameter :: c42=c32*c10-s32*s10, s42=s32*c10+c32*s10 
    Real, Parameter :: c43=c32*c11-s32*s11, s43=s32*c11+c32*s11 
    Real, Parameter :: c44=c32*c12-s32*s12, s44=s32*c12+c32*s12 
    Real, Parameter :: c45=c32*c13-s32*s13                      
    
    Integer, Parameter :: n90  = 90           
    Integer, Parameter :: i0   = 181          
    Integer, Parameter :: i9   = i0 + n90      
    Integer, Parameter :: m9   = i0 - n90     
    Integer, Parameter :: i8   = i9 + n90     
    Real, Parameter :: pi   = 3.141592654  
    Real, Parameter :: cang = (pi+pi) / NN 
    Integer :: i

    character(len=*),parameter :: myname_=myname//'::qtrig0'

    if(initialized_) call die(myname_,'multiple object definition')
    initialized_=.true.

    if ( NN == 360 ) then

       trigtab(i9- 0)  = c00
       trigtab(i0+ 0)  = s00
       trigtab(i9- 1)  = c01
       trigtab(i0+ 1)  = s01
       trigtab(i9- 2)  = c02
       trigtab(i0+ 2)  = s02
       trigtab(i9- 3)  = c03
       trigtab(i0+ 3)  = s03
       trigtab(i9- 4)  = c04
       trigtab(i0+ 4)  = s04
       trigtab(i9- 5)  = c05
       trigtab(i0+ 5)  = s05
       trigtab(i9- 6)  = c06
       trigtab(i0+ 6)  = s06
       trigtab(i9- 7)  = c07
       trigtab(i0+ 7)  = s07
       trigtab(i9- 8)  = c08
       trigtab(i0+ 8)  = s08
       trigtab(i9- 9)  = c09
       trigtab(i0+ 9)  = s09
       trigtab(i9-10)  = c10
       trigtab(i0+10)  = s10
       trigtab(i9-11)  = c11
       trigtab(i0+11)  = s11
       trigtab(i9-12)  = c12
       trigtab(i0+12)  = s12
       trigtab(i9-13)  = c13
       trigtab(i0+13)  = s13
       trigtab(i9-14)  = c14
       trigtab(i0+14)  = s14
       trigtab(i9-15)  = c15
       trigtab(i0+15)  = s15
       trigtab(i9-16)  = c16
       trigtab(i0+16)  = s16
       trigtab(i9-17)  = c17
       trigtab(i0+17)  = s17
       trigtab(i9-18)  = c18
       trigtab(i0+18)  = s18
       trigtab(i9-19)  = c19
       trigtab(i0+19)  = s19
       trigtab(i9-20)  = c20
       trigtab(i0+20)  = s20
       trigtab(i9-21)  = c21
       trigtab(i0+21)  = s21
       trigtab(i9-22)  = c22
       trigtab(i0+22)  = s22
       trigtab(i9-23)  = c23
       trigtab(i0+23)  = s23
       trigtab(i9-24)  = c24
       trigtab(i0+24)  = s24
       trigtab(i9-25)  = c25
       trigtab(i0+25)  = s25
       trigtab(i9-26)  = c26
       trigtab(i0+26)  = s26
       trigtab(i9-27)  = c27
       trigtab(i0+27)  = s27
       trigtab(i9-28)  = c28
       trigtab(i0+28)  = s28
       trigtab(i9-29)  = c29
       trigtab(i0+29)  = s29
       trigtab(i9-30)  = c30
       trigtab(i0+30)  = s30
       trigtab(i9-31)  = c31
       trigtab(i0+31)  = s31
       trigtab(i9-32)  = c32
       trigtab(i0+32)  = s32
       trigtab(i9-33)  = c33
       trigtab(i0+33)  = s33
       trigtab(i9-34)  = c34
       trigtab(i0+34)  = s34
       trigtab(i9-35)  = c35
       trigtab(i0+35)  = s35
       trigtab(i9-36)  = c36
       trigtab(i0+36)  = s36
       trigtab(i9-37)  = c37
       trigtab(i0+37)  = s37
       trigtab(i9-38)  = c38
       trigtab(i0+38)  = s38
       trigtab(i9-39)  = c39
       trigtab(i0+39)  = s39
       trigtab(i9-40)  = c40
       trigtab(i0+40)  = s40
       trigtab(i9-41)  = c41
       trigtab(i0+41)  = s41
       trigtab(i9-42)  = c42
       trigtab(i0+42)  = s42
       trigtab(i9-43)  = c43
       trigtab(i0+43)  = s43
       trigtab(i9-44)  = c44
       trigtab(i0+44)  = s44
       trigtab(i9-45)  = c45

       do i = 0, n90
	  trigtab( i8 - i ) = + trigtab( i0 + i )
	  trigtab( i8 + i ) = - trigtab( i0 + i ) 
	  trigtab( m9 - i ) = - trigtab( i9 - i )
	  trigtab( m9 + i ) = - trigtab( i9 - i )
       end do

    else

       do  i = 1, lentab
	  trigtab(i) = sin(cang*(i-1)-pi)
       end do

    endif

    return
  end subroutine qtrig0


  ! This function computes the circumcenter of a 
  ! spherical triangle. 
  Function circumcenter_(v, cosrad) Result(center)
    ! Inputs
    Real, Intent(In) :: v(3,3) ! vertices of triangle in clockwise order
    ! Output
    Real :: center(3)
    Real, Optional :: cosrad   ! cosine of the radius of the circumcircle

    center = cross(v(:,2),v(:,1)) + cross(v(:,3),v(:,2)) + cross(v(:,1),v(:,3))
    center = normalize(center)

    If (Present(cosrad)) cosrad = dot_product(v(:,1),center)

  End Function circumcenter_

  ! returns unit vector in the direction of v
  function normalize_vec(v) result(nv)
    Real, Intent(In) :: v(3)
    Real :: nv(3)
    nv = v / sqrt(sum(v**2))
  end function normalize_vec

  ! normalize a list of vectors
  Function normalize_vec_list(v) result(nv)
    Real, Intent(In) :: v(:,:)
    Real :: nv(3,size(v,2))
    Real :: temp
    Integer :: i
    
    Do i = 1, size(v,2)
       nv(:,i) = normalize(v(:,i))
    End Do
  End Function normalize_vec_list

  ! Vector cross product
  Function cross_(a, b) result(axb)
    Real, Intent(In) :: a(3), b(3)
    Real :: axb(3)
    
    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)
    
  End Function cross_


  ! Vector dot product
  Real Function dot(p_a, p_b)
    Real, Intent(In) :: p_a(3), p_b(3)
    
    dot = p_a(1)*p_b(1) + p_a(2)*p_b(2) + p_a(3)*p_b(3)
    
  End Function dot
  
  Function cyclic_cross_(v123) Result(tri)
    Real, Intent(in) :: v123(3,3)
    Real :: tri(3,3)

    tri(:,1) = cross(v123(:,2),v123(:,3))
    tri(:,2) = cross(v123(:,3),v123(:,1))
    tri(:,3) = cross(v123(:,1),v123(:,2))

  End Function cyclic_cross_

End Module m_spherical_geometry

