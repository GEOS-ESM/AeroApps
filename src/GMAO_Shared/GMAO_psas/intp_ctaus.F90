subroutine intp_ctaus(nctau1,dctau1,nctau2,dctau2,nctau,ctaus)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: intp_ctaus - (to do)
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	25Jan96 - J. Guo	- (to do)
!_______________________________________________________________________
#ifndef _REAL_
# define	_REAL_	real*8
#endif

use m_stdio,only : stderr
use m_die,  only : die
implicit none
integer, intent(in)	:: nctau1
real,    intent(in)	:: dctau1
integer, intent(in)	:: nctau2
real,    intent(in)	:: dctau2
integer, intent(in)	:: nctau
_REAL_,    intent(out)	:: ctaus(0:nctau-1)

!  local vars.
  integer i

  character(len=*), parameter	:: myname='intp_ctaus'

	! verify the argument

  if(nctau.ne.nctau1+nctau2+1) then
    write(stderr,'(2a,i5)') myname,	&
      ': invalid argument, $5=$1+$3+1 expected but',nctau
    call die(myname)
  endif
  
  ctaus(0)=0.		! the first value is given
  
  	! region 1, i= 0,nctau1 = int(.5 + (ctau)/dctau1)
  do i=1,nctau1
    ctaus(i)=dctau1*i
  end do
  
  	! region 2, i= nctau1+1,nctau2 = int(.5 + (ctau-ctau1)/dctau2)
  do i=1,nctau2
    ctaus(i+nctau1)=dctau2*i+ctaus(nctau1)
  end do
end
!.
