module hcorfuns
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: hcorfuns - (to do)
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	17May96 - J. Guo	- (to do)
!       16dec99 - da Silva  Added WIN_POWER-LAW-3000km
!_______________________________________________________________________

#ifndef _REAL_
#  define	_REAL_	real*8
#endif

private	! except

public	:: i_DAMP_COSINE
public	:: i_POWER_LAW
public	:: i_GAUSSIAN
public	:: i_EXPONENTIAL
public	:: i_GASPARI_COHN
public	:: i_WIN_POWERLAW
public	:: i_WIN_POWERLAW_6000km
public	:: i_WIN_POWERLAW_3000km

!public	:: damp_cosine
!public	:: power_law
!public	:: gaussian
!public	:: exponential
!public	:: gaspari_cohn

public	:: corfun
public	:: richxtr

public	:: i_funcname

integer,          parameter	:: i_DAMP_COSINE  = 1
character(len=*), parameter	:: c_DAMP_COSINE  = 'DAMP-COSINE'

integer,          parameter	:: i_POWER_LAW    = 2
character(len=*), parameter	:: c_POWER_LAW    = 'POWER-LAW'

integer,          parameter	:: i_GAUSSIAN     = 3
character(len=*), parameter	:: c_GAUSSIAN     = 'GAUSSIAN'

integer,          parameter	:: i_EXPONENTIAL  = 4
character(len=*), parameter	:: c_EXPONENTIAL  = 'EXPONENTIAL'

integer,          parameter	:: i_GASPARI_COHN = 5
character(len=*), parameter	:: c_GASPARI_COHN = 'GASPARI-COHN'

integer,          parameter	:: i_WIN_POWERLAW = 6
character(len=*), parameter	:: c_WIN_POWERLAW = 'WIN-POWERLAW'

integer,          parameter	:: i_WIN_POWERLAW_6000km = 6
character(len=*), parameter	:: c_WIN_POWERLAW_6000km = 'WIN-POWERLAW-6000km'

integer,          parameter	:: i_WIN_POWERLAW_3000km = 7
character(len=*), parameter	:: c_WIN_POWERLAW_3000km = 'WIN-POWERLAW-300'

! Note: on rc file, you can also specify: WIN-POWERLAW-3000km

character(len=*), parameter :: myname='hcorfuns'

contains
!=======================================================================
integer function i_funcname(c_funcname)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: i_funcname - return an integer code for a give function name
!
! !INTERFACE:
!	<@interface
!
! !DESCRIPTION:
!
!	For a user specified function name (c_funcname), i_funcname()
!	returns an integer as the function value from a set of function
!	names known to the module.  If the name is not known, -1 is
!	returned.
!
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	24May96 - J. Guo	- (to do)
!_______________________________________________________________________
implicit none
character(len=*),intent(in)	:: c_funcname

  character(len=len(c_funcname)) :: UCType
  integer i,l,ifun

  l=len_trim(c_funcname)
  UCType=c_funcname
  do i=1,l
    if(UCType(i:i).ge.'a'.and.UCType(i:i).le.'z')	&
	& UCType(i:i) = char(ichar(UCType(i:i))+ichar('A')-ichar('a'))
  end do

  ifun=-1

  select case(UCType(1:l))
  case (c_DAMP_COSINE)
    ifun=i_DAMP_COSINE

  case (c_POWER_LAW)
    ifun=i_POWER_LAW

  case (c_GAUSSIAN)
    ifun=i_GAUSSIAN

  case (c_EXPONENTIAL)
    ifun=i_EXPONENTIAL

  case (c_GASPARI_COHN)
    ifun=i_GASPARI_COHN

  case (c_WIN_POWERLAW)
    ifun=i_WIN_POWERLAW

  case (c_WIN_POWERLAW_6000km)
    ifun=i_WIN_POWERLAW        ! yes, this is the same as c_WIN_POWERLAW

  case (c_WIN_POWERLAW_3000km)
    ifun=i_WIN_POWERLAW_3000km

  end select

  i_funcname=ifun

!@interface
end function i_funcname
!@end/interface
!=======================================================================
!@interface
subroutine damp_cosine(nctau,ctaus,dctau,cors,npar,pars,ierr)
!@end/interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: damp_cosine - (to do)
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	17May96 - J. Guo	- (to do)
!_______________________________________________________________________
!@usemodule
use const, only : R_EARTH
!@end/usemodule

!@interface
implicit none

  integer, intent(in)	:: nctau	! size of ctaus[]
  _REAL_,  intent(in)	:: ctaus(nctau)	! function argument table
  _REAL_,  intent(in)	:: dctau	! offset from ctaus(:)
  _REAL_,  intent(out)	:: cors(nctau)	! function values
  integer, intent(in)	:: npar		! number of parameters
  real,    intent(in)	:: pars(npar)	! function parameter table
  integer, intent(out)	:: ierr		! status code

!@end/interface
!-----------------------------------------------------------------------
	! Local parameters
  integer, parameter	:: mpar = 6	! minimum parameters
  _REAL_,  parameter	:: rade = .001*R_EARTH
  character(len=*),parameter :: myname='damp_cosine'
	!----------------------------------------
	! Local vars.
  _REAL_ dm,c1,c2,c3,c4,c5
  _REAL_ s
  integer i
!-----------------------------------------------------------------------
	! Checking

  ierr=0
  if(npar.lt.mpar) then
    ierr=-1
    return
  endif
	!----------------------------------------
	! define parameters

  dm=pars(1)	! not used
  c1=pars(2)
  c2=pars(3)*rade
  c3=pars(4)
  c4=pars(5)*rade
  c5=pars(6)

	!----------------------------------------
  do i=1,nctau
    s=acos(1.-ctaus(i)-dctau)
    cors(i) = (c1*cos(c2*s)+c3)/(1.0+(c4*s)*(c4*s))**c5
  end do

end subroutine damp_cosine
!=======================================================================
!@interface
subroutine power_law(nctau,ctaus,dctau,cors,npar,pars,ierr)
!@end/interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: power_law - (to do)
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	17May96 - J. Guo	- (to do)
!_______________________________________________________________________
!@usemodule
use const, only : R_EARTH
!@end/usemodule

!@interface
implicit none

  integer, intent(in)	:: nctau	! size of ctaus[]
  _REAL_,  intent(in)	:: ctaus(nctau)	! function argument table
  _REAL_,  intent(in)	:: dctau	! offset from ctaus(:)
  _REAL_,  intent(out)	:: cors(nctau)	! function values
  integer, intent(in)	:: npar		! number of parameters
  real,    intent(in)	:: pars(npar)	! function parameter table
  integer, intent(out)	:: ierr		! status code

!@end/interface
!-----------------------------------------------------------------------
	! Local parameters
  integer, parameter	:: mpar = 2	! minimum parameters
  _REAL_,  parameter	:: rade = .001*R_EARTH
  character(len=*),parameter :: myname='power_law'
	!----------------------------------------
	! Local vars.
  _REAL_ dm,L,b
  _REAL_ r
  integer i
!-----------------------------------------------------------------------
	! Checking

  ierr=0
  if(npar.lt.mpar) then
    ierr=-1
    return
  endif
	!----------------------------------------
	! define parameters

  dm=pars(1)	! not used
  L =pars(2)/rade
  b =1.
  if(npar.ge.3) b =pars(3)

	!----------------------------------------
  if(npar.eq.2) then
    do i=1,nctau
      ! as suggested by Gaspari, to be tested
      ! r=sqrt(2.*(ctaus(i)+dctau))/L

      r=2.*sin(.5*acos(1.-ctaus(i)-dctau))/L
      cors(i) = 1./(1.+r*r)
    end do
  else
    do i=1,nctau
      ! as suggested by Gaspari, to be tested
      ! r=sqrt(2.*(ctaus(i)+dctau))/L

      r=2.*sin(.5*acos(1.-ctaus(i)-dctau))/L
      cors(i) = (1.+r*r)**(-b)
    end do
  endif

end subroutine power_law
!=======================================================================
!@interface
subroutine gaussian(nctau,ctaus,dctau,cors,npar,pars,ierr)
!@end/interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: gaussian - (to do)
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	17May96 - J. Guo	- (to do)
!_______________________________________________________________________
!@usemodule
use const, only : R_EARTH
!@end/usemodule

!@interface
implicit none

  integer, intent(in)	:: nctau	! size of ctaus[]
  _REAL_,  intent(in)	:: ctaus(nctau)	! function argument table
  _REAL_,  intent(in)	:: dctau	! offset from ctaus(:)
  _REAL_,  intent(out)	:: cors(nctau)	! function values
  integer, intent(in)	:: npar		! number of parameters
  real,    intent(in)	:: pars(npar)	! function parameter table
  integer, intent(out)	:: ierr		! status code

!@end/interface
!-----------------------------------------------------------------------
	! Local parameters
  integer, parameter	:: mpar = 2	! minimum parameters
  _REAL_,  parameter	:: rade = .001*R_EARTH
  character(len=*),parameter :: myname='gaussian'
	!----------------------------------------
	! Local vars.
  _REAL_ dm,L
  _REAL_ s
  integer i
!-----------------------------------------------------------------------
	! Checking

  ierr=0
  if(npar.lt.mpar) then
    ierr=-1
    return
  endif
	!----------------------------------------
	! define parameters

  dm=pars(1)	! not used
  L =pars(2)/rade

	!----------------------------------------
  do i=1,nctau
    s=acos(1.-ctaus(i)-dctau)/L
    cors(i) = exp(-s*s)
  end do

end subroutine gaussian
!=======================================================================
!@interface
subroutine exponential(nctau,ctaus,dctau,cors,npar,pars,ierr)
!@end/interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: exponential - (to do)
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	17May96 - J. Guo	- (to do)
!_______________________________________________________________________
!@usemodule
use const, only : R_EARTH
!@end/usemodule

!@interface
implicit none

  integer, intent(in)	:: nctau	! size of ctaus[]
  _REAL_,  intent(in)	:: ctaus(nctau)	! function argument table
  _REAL_,  intent(in)	:: dctau	! offset from ctaus(:)
  _REAL_,  intent(out)	:: cors(nctau)	! function values
  integer, intent(in)	:: npar		! number of parameters
  real,    intent(in)	:: pars(npar)	! function parameter table
  integer, intent(out)	:: ierr		! status code

!@end/interface
!-----------------------------------------------------------------------
	! Local parameters
  integer, parameter	:: mpar = 2	! minimum parameters
  _REAL_,  parameter	:: rade = .001*R_EARTH
  character(len=*),parameter :: myname='exponential'
	!----------------------------------------
	! Local vars.
  _REAL_ dm,L,b
  _REAL_ s
  integer i
!-----------------------------------------------------------------------
	! Checking

  ierr=0
  if(npar.lt.mpar) then
    ierr=-1
    return
  endif
	!----------------------------------------
	! define parameters

  dm=pars(1)	! not used
  L =pars(2)/rade

	!----------------------------------------
  do i=1,nctau
    s=acos(1.-ctaus(i)-dctau)/L
    cors(i) = exp(-s)
  end do

end subroutine exponential
!=======================================================================
!@interface
subroutine gaspari_cohn(nctau,ctaus,dctau,cors,npar,pars,ierr)
!@end/interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: gaspari_cohn - (to do)
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	17May96 - J. Guo	- (to do)
!_______________________________________________________________________
!@usemodule
use const, only : R_EARTH
!@end/usemodule

!@interface
implicit none

  integer, intent(in)	:: nctau	! size of ctaus[]
  _REAL_,  intent(in)	:: ctaus(nctau)	! function argument table
  _REAL_,  intent(in)	:: dctau	! offset from ctaus(:)
  _REAL_,  intent(out)	:: cors(nctau)	! function values
  integer, intent(in)	:: npar		! number of parameters
  real,    intent(in)	:: pars(npar)	! function parameter table
  integer, intent(out)	:: ierr		! status code

!@end/interface
!-----------------------------------------------------------------------
	! Local parameters
  integer, parameter	:: mpar = 2	! minimum parameters
  _REAL_,  parameter	:: rade = .001*R_EARTH
  character(len=*),parameter :: myname='gaspari_cohn'

	!----------------------------------------

  _REAL_, parameter :: a0= 1.
  _REAL_, parameter :: a2=-5./3.
  _REAL_, parameter :: a3= 5./8.
  _REAL_, parameter :: a4= 1./2.
  _REAL_, parameter :: a5=-1./4.

  _REAL_, parameter :: bm=-2./3.
  _REAL_, parameter :: b0= 4.
  _REAL_, parameter :: b1=-5.
  _REAL_, parameter :: b2= 5./3.
  _REAL_, parameter :: b3= 5./8.
  _REAL_, parameter :: b4=-1./2.
  _REAL_, parameter :: b5= 1./12.

	!----------------------------------------
	! Local vars.
  _REAL_ dm,L
  _REAL_ r,c,z,w
  integer i
!-----------------------------------------------------------------------
	! Checking

  ierr=0
  if(npar.lt.mpar) then
    ierr=-1
    return
  endif

#ifdef	_DEBUG
  _DEBUG  write(*,*) myname,': #1 - ',nctau,npar,ierr,ctaus(1)
#endif
	!----------------------------------------
	! define parameters

  dm=pars(1)	! not used
  L =pars(2)/rade
  c =sqrt(10.d0/3.d0)*L

#ifdef	_DEBUG
  _DEBUG  write(*,*) myname,': #2 - ',nctau,npar,ierr,ctaus(1)
#endif
	!----------------------------------------
  do i=1,nctau

!  The follow #ifndef structure is put in for reproducing "etest14"
!  result.

#ifndef _ETEST14_
    w=sqrt(2.*ctaus(i))/c		! reference argument
    z=sqrt(2.*(ctaus(i)+dctau))/c	! actual argument
#else
    w=2.*sin(.5*acos(1.-ctaus(i)))/c		! reference argument
    z=2.*sin(.5*acos(1.-ctaus(i)-dctau))/c	! actual argument
#endif

    if(w.le.1.) then
      cors(i) = z*z*(z*(z*(z*a5+a4)+a3)+a2) + a0
    elseif(w.lt.2.) then
      cors(i) = z*(z*(z*(z*(z*b5+b4)+b3)+b2)+b1) + b0 + bm/z
    else
      cors(i) = 0.
    endif

  end do

#ifdef	_DEBUG
  _DEBUG  write(*,*) myname,': #3 - ',nctau,npar,ierr,ctaus(1)
#endif
end subroutine gaspari_cohn

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: win_powerlaw_ - interpreting a windowed-power-law function
!
! !DESCRIPTION:
!
!   This windowed-powerlaw correlation function uses a compact-
!   supported correlation function (see subroutine gaspary_cohn() above)
!   for the window function with a "cutoff" distance (Cw=3000./Re by default).
!   Notice the definition of the power-law function used here is also
!   different from the one used in subroutine power_law().
!
! !INTERFACE:

    subroutine win_powerlaw_(nctau,ctaus,dctau,cors,npar,pars,ierr, &
                             cutoff ) ! optional
      use const, only : R_EARTH
      implicit none

      integer,             intent(in)  :: nctau	! size of ctaus(:)
      _REAL_, dimension(:),intent(in)  :: ctaus	! argument table
      _REAL_,              intent(in)  :: dctau	! offset from ctaus(:)
      _REAL_, OPTIONAL,    intent(in)  :: cutoff ! wind cuttof [km]
      _REAL_, dimension(:),intent(out) :: cors	! function values

      integer,             intent(in)  :: npar	! number of parameters
      real,   dimension(:),intent(in)  :: pars	! parameter table

      integer,             intent(out) :: ierr	! status code

! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	13Mar97 - Jing Guo <guo@eramus> - defined template
!       15dec99 - da Silva - added cutoff optional parameter
!_______________________________________________________________________

	! Local parameters
  integer, parameter	:: mpar = 2	! minimum parameters
  _REAL_,  parameter	:: rade = .001*R_EARTH

  _REAL_                :: Cw   ! cutoff parameter


		! 

  character(len=*),parameter :: myname_=myname//'::win_powerlaw_'

	!----------------------------------------

  _REAL_, parameter :: a0= 1.
  _REAL_, parameter :: a2=-5./3.
  _REAL_, parameter :: a3= 5./8.
  _REAL_, parameter :: a4= 1./2.
  _REAL_, parameter :: a5=-1./4.

  _REAL_, parameter :: bm=-2./3.
  _REAL_, parameter :: b0= 4.
  _REAL_, parameter :: b1=-5.
  _REAL_, parameter :: b2= 5./3.
  _REAL_, parameter :: b3= 5./8.
  _REAL_, parameter :: b4=-1./2.
  _REAL_, parameter :: b5= 1./12.

	!----------------------------------------
	! Local vars.
  _REAL_ dm,L,Lw
  _REAL_ r,z,w
  integer i
!-----------------------------------------------------------------------
	! Checking

! Set cuttoff parameter Cw. Notice the "Cw" is actually the half distance
! where the correlation vanishes. 
! -----------------------------------------------------------------------
  if ( present(cutoff) ) then
       Cw = 0.5 * cutoff / rade
   else
       Cw = 0.5 * 6000. / rade
   end if

       
  ierr=0
  if(npar.lt.mpar) then
    ierr=-1
    return
  endif

#ifdef	_DEBUG
  _DEBUG  write(*,*) myname,': #1 - ',nctau,npar,ierr,ctaus(1)
#endif
	!----------------------------------------
	! define parameters

  dm=pars(1)	! not used
  L =pars(2)/rade

  Lw=Cw/sqrt(10.d0/3.d0)	! the length scale of the window

#ifdef	_DEBUG
  _DEBUG  write(*,*) myname,': #2 - ',nctau,npar,ierr,ctaus(1)
#endif
	!----------------------------------------
  do i=1,nctau
		! Computer dimesionless cordal distances

    w=sqrt(2.*ctaus(i))/Cw		! reference argument
    z=sqrt(2.*(ctaus(i)+dctau))/Cw	! actual argument
    r=sqrt(2.*(ctaus(i)+dctau))/L

		! Evaluate the function values
    if(w.le.1.) then
      cors(i)=(z* z*(z*(z*(z*a5+a4)+a3)+a2)    +a0     ) /(1.+.5*r*r)
    elseif(w.lt.2.) then
      cors(i)=(z*(z*(z*(z*(z*b5+b4)+b3)+b2)+b1)+b0+bm/z) /(1.+.5*r*r)
    else
      cors(i)=0.
    endif

  end do

#ifdef	_DEBUG
  _DEBUG  write(*,*) myname,': #3 - ',nctau,npar,ierr,ctaus(1)
#endif
end subroutine win_powerlaw_
!=======================================================================
subroutine richxtr(nctau,ctaus,dctau,hcorf,dcorf,qcorf,		&
  &	i_CorType,npar,pars,niter,ierr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! !ROUTINE: richxtr - evaluate derivatives of an analytic function
!
! !DESCRIPTION:
!	richxtr() uses repeated Richardson extrapolation algorithm to
!	evaluate derivatives of an analytic function.
!_______________________________________________________________________
use m_stdio, only : stderr

!@interface

implicit none
integer, intent(in)	:: nctau	! size of the tables
_REAL_,  intent(in)	:: ctaus(nctau)	! function argument table
_REAL_,  intent(in)	:: dctau	! starting step size
_REAL_,  intent(in)	:: hcorf(nctau)	! function value table (done)
_REAL_,  intent(out)	:: dcorf(nctau)	! first derivatives
_REAL_,  intent(out)	:: qcorf(nctau)	! second derivatives
	!----------------------------------------
integer, intent(in)	:: i_CorType	! Function name
integer, intent(in)	:: npar		! size of the parameter table
real,    intent(in)	:: pars(npar)	! a parameter table
integer, intent(in)	:: niter	! number of iterations
	!----------------------------------------
integer, intent(out)	:: ierr		! status code
!@end/interface
!-----------------------------------------------------------------------
	! Local parameters

  integer,             parameter	:: np2 = 2	! power of 2

	! The value of np2 is determined by the reduction of the step h
	! as the base (2), as well as the power series step of the
	! evaluated function.  In this case, the power series step is
	! taken as 1 for generality.

  character(len=*), parameter	:: myname = 'richxtr'

!-----------------------------------------------------------------------
! Additional information from my older version:
!
!	..Numerical approximation using "Richardson's deferred approch
!	to the limit" algorithm with any given function form.  For
!	reference, see,
!
!	+ Dahlquist, G. and A. Bjorck, 1974, Numerical Methods,
!	Prentic-Hall.
!
!	..The function name given at the calling of the function is
!	used to evaluate a series of the given function, which is
!	determined at a given x=t and a series of value of h (=dt at
!	the begining).  The value of h is halfed at every iteration,
!	until a condition about the error is meet.
!-----------------------------------------------------------------------
	! Local workspace and variables

  integer		:: mb,i,k,nm
  _REAL_		:: h
  _REAL_		:: fq(1)
  _REAL_,allocatable	:: fm(:)
  _REAL_,allocatable	:: fp(:)

  _REAL_,allocatable	:: at_d(:),am_d(:,:)
  _REAL_,allocatable	:: at_q(:),am_q(:,:)

	!----------------------------------------
	! The initial step should not be too small!

	h=dctau

	!----------------------------------------
	! Allocate workspace

	allocate( fm(nctau),  fp(nctau),			&
	  &	  at_d(nctau),am_d(nctau,0:niter),		&
	  &	  at_q(nctau),am_q(nctau,0:niter),		&
	  &	  stat=ierr					)

	if(ierr.ne.0) then
	  write(stderr,'(2a,i5)') myname,			&
	    &	': allocate() error, stat =',ierr
	  return
	endif

	!----------------------------------------
!	..First differentiation,

	call corfun(      1,ctaus(1), 2.*h,fq(1),	&
	  & i_CorType,npar,pars, ierr			)
	if(ierr.ne.0) then
	  write(stderr,'(2a,i5)') myname,			&
	    &	': fq=corfun(tq) error, stat =',ierr
	  deallocate(fm,fp,at_d,at_q,am_d,am_q)
	  return
	endif

	call corfun(nctau-1,ctaus(2),-1.*h,fm(2),	&
	  & i_CorType,npar,pars, ierr			)
	if(ierr.ne.0) then
	  write(stderr,'(2a,i5)') myname,			&
	    &	': fm=corfun(tm) error, stat =',ierr
	  deallocate(fm,fp,at_d,at_q,am_d,am_q)
	  return
	endif

	call corfun(nctau  ,ctaus(1), 1.*h,fp(1),	&
	  & i_CorType,npar,pars, ierr			)
	if(ierr.ne.0) then
	  write(stderr,'(2a,i5)') myname,			&
	    &	': fp=corfun(tp) error, stat =',ierr
	  deallocate(fm,fp,at_d,at_q,am_d,am_q)
	  return
	endif

!	..Store data for extrapolations.  Note that dm or am_x(0) is the
!	directly computed value using the formula, not the returned
!	derivative, which is obtained from using the Richardson's
!	extrapolation.

		! Boundary values are handled specially

	am_d(1,0)=(3.*fp(1)-2.*hcorf(1)-fq(1))/h
	am_q(1,0)=(fq(1)-2.*fp(1)+hcorf(1))/(h*h)

		! Other are handled using regular difference expressions
	do i=2,nctau
	  am_d(i,0)=.5*(fp(i)-fm(i))/h
	  am_q(i,0)=(fp(i)-2.*hcorf(i)+fm(i))/(h*h)
	end do


!	..Print a table to verify the result
!	print'(i2,1p,e10.3,0p,10(f10.6,i5,1x),f10.6)',nm,h,
!     &	  (am(k-1)*scale,nint((am(k)-am(k-1))*scale*1.e+6),k=1,nm),
!     &	  am(nm)*scale
	
	nm=0
	do while(nm.lt.niter)
	  nm=nm+1

!	  ..New values of h=h/2 and derivative.
	  h=h/2

	  call corfun(      1,ctaus(1), 2.*h,fq(1),	&
	    &	i_CorType,npar,pars, ierr		)
	  if(ierr.ne.0) then
	    write(stderr,'(2a,i2,a,i5)') myname,		&
	      &	': fq=corfun(tq,',nm,') error, stat =',ierr
	    deallocate(fm,fp,at_d,at_q,am_d,am_q)
	    return
	  endif

	  call corfun(nctau-1,ctaus(2),-1.*h,fm(2),	&
	    &	i_CorType,npar,pars, ierr		)
	  if(ierr.ne.0) then
	    write(stderr,'(2a,i2,a,i5)') myname,		&
	      &	': fm=corfun(tm,',nm,') error, stat =',ierr
	    deallocate(fm,fp,at_d,at_q,am_d,am_q)
	    return
	  endif
	  
	  call corfun(nctau  ,ctaus(1), 1.*h,fp(1),	&
	    &	i_CorType,npar,pars, ierr		)
	  if(ierr.ne.0) then
	    write(stderr,'(2a,i2,a,i5)') myname,		&
	      &	': fp=corfun(tp,',nm,') error, stat =',ierr
	    deallocate(fm,fp,at_d,at_q,am_d,am_q)
	    return
	  endif

!-----------------------------------------------------------------------
!	  ..Extrapolate with Richardson extrapolation.

		!------------------------------------------------
		! am_d(:,0:nm) = ..., using at_q(:) as a buffer

	  at_d(1)=(3.*fp(1)-2.*hcorf(1)-fq(1))/h
	  do i=2,nctau
	    at_d(i)=.5*(fp(i)-fm(i))/h
	  end do

	  mb=1
	  do k=1,nm
	    mb=mb*np2

	    at_q(:)=at_d(:)+(at_d(:)-am_d(:,k-1))/(mb-1)
	    am_d(:,k-1)=at_d(:)
	    at_d(:)=at_q(:)
	  end do

	  am_d(:,nm)=at_d(:)

		!------------------------------------------------
		! am_q(:,0:nm) = ..., using at_d(:) as a buffer

	  at_q(1)=(fq(1)-2.*fp(1)+hcorf(1))/(h*h)
	  do i=2,nctau
	    at_q(i)=(fp(i)-2.*hcorf(i)+fm(i))/(h*h)
	  end do

	  mb=1
	  do k=1,nm
	    mb=mb*np2

	    at_d(:)=at_q(:)+(at_q(:)-am_q(:,k-1))/(mb-1)
	    am_q(:,k-1)=at_q(:)
	    at_q(:)=at_d(:)

	  end do
	  am_q(:,nm)=at_q(:)

!	  ..If the correction term is not greater than the minimum
!	  fraction permited by numerical computation |era*a|, or it is
!	  less than the least significant magnitude given by the user
!	  before the computation |erf|, there is no point to continue
!	  the computation.  Don't let erf = 0!  Also, |h| ~= 0 is not
!	  checked.

!	  ..Print a table to verify the result
!	  print'(i2,1p,e10.3,0p,10(f10.6,i5,1x),f10.6)',nm,h,
!     &	    (am(k-1)*scale,nint((am(k)-am(k-1))*scale*1.e+6),k=1,nm),
!     &	    am(nm)*scale

	end do

!	print'(1x,a,5x,a,4x,10(5x,i2,9x),7x,i2)','m','h',
!     &	  (k-1,k=1,nm),nm

!	da=am(nm)-am(nm-1)
!	print'(a,2(a,1p,e10.3,a,l1),a,i2,x,l1,a,1p,e10.3,a,l1)',
!     &	  'Condition: ','<=era(',abs(era*a),') ',abs(da).le.abs(era*a),
!     &	  ' <erf(',erf,') ',abs(da).lt.erf,' nm=20=',nm,nm.eq.mmax

	dcorf(:)=am_d(:,nm)
	qcorf(:)=am_q(:,nm)

	deallocate(fm,fp,at_d,at_q,am_d,am_q)
!-----------------------------------------------------------------------
end subroutine richxtr
!=======================================================================
subroutine corfun(nctau,ctaus,dctau,hcors,i_CorType,npar,pars,ierr)
use m_stdio, only : stderr
  implicit none
  integer, intent(in)	:: nctau	! size of ctaus(:)
  _REAL_,  intent(in)	:: ctaus(nctau)	! argument table
  _REAL_,  intent(in)	:: dctau	! offset from ctaus(:)
  _REAL_,  intent(out)	:: hcors(nctau)	! function value table
  integer, intent(in)	:: i_CorType	! function type code
  integer, intent(in)	:: npar		! size of pars[]
  real,    intent(in)	:: pars(npar)	! function parameter table
  integer, intent(out)	:: ierr		! status code
!-----------------------------------------------------------------------
  character(len=*), parameter	:: myname='corfun'
!-----------------------------------------------------------------------
!#define _DEBUG

  _REAL_ cutoff

  ierr=0	! should be initialized
  if(nctau.le.0) return

#ifdef	_DEBUG
  write(*,*) myname,': #1 - ',nctau,npar,ierr,ctaus(1:2),i_CorType
#endif

  select case(i_CorType)
    case ( i_DAMP_COSINE)
      call damp_cosine(	nctau,ctaus,dctau,hcors, npar,pars, ierr )

    case (    i_GAUSSIAN)
      call gaussian(	nctau,ctaus,dctau,hcors, npar,pars, ierr )

    case ( i_EXPONENTIAL)
      call exponential(	nctau,ctaus,dctau,hcors, npar,pars, ierr )

    case (   i_POWER_LAW)
      call power_law(	nctau,ctaus,dctau,hcors, npar,pars, ierr )

    case (i_GASPARI_COHN)
      call gaspari_cohn(nctau,ctaus,dctau,hcors, npar,pars, ierr )

    case (i_WIN_POWERLAW)
      call win_powerlaw_(nctau,ctaus,dctau,hcors, npar,pars, ierr )

    case (i_WIN_POWERLAW_3000km)
      cutoff = 3000.
      call win_powerlaw_(nctau,ctaus,dctau,hcors, npar,pars, ierr, cutoff )

    case default
      write(stderr,'(2a,i2)') myname,		&
	& ': unknown function code, ',i_CorType
      ierr=-1
      return
    end select

    if(ierr /= 0) then
      write(stderr,'(2a,i2,a,i5)') myname,		&
	& ': error with function code ',i_CorType, ', ierr =',ierr
      ierr=-2
      return
    endif
#ifdef	_DEBUG
  write(*,*) myname,': #2 - ',nctau,npar,ierr,ctaus(1:2),i_CorType
#endif
  end subroutine corfun
!=======================================================================
end module hcorfuns
!.
