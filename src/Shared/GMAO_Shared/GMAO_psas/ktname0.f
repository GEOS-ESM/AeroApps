
c                                                   02/05/93 - ktname0.f
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Initialize array ktname with the known data type names.

      subroutine ktname0
      use m_die,  only : die
      use m_stdio,only : stdout

c.......................................................................
c.... Table ktname to hold data type names.

	use m_output, only : output_ObsErrCov
	use m_stdio
	implicit none

      include     "ktmax.h"
      include     "kttabl.h"

c.......................................................................
c.... Data type names.

c     DATA ktname /
c    1 'u_SeaLev' , 'v_SeaLev' , 'p_SeaLev' , 'u_UprAir' , 'v_UprAir',
c    2 'H_UprAir' , 'q_UprAir' /


	character*7 myname
	parameter  (myname='ktname0')

	integer istat

	integer ikt,jkt

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	call rdkttbl(RC_kt,ktmax,ktname,ktunit,ktdesc,ktmvar,istat)
	if(istat.ne.0) call die(myname,'rdkttbl()',istat)

c	..Echo input in its original form
c	=================================
	if(output_ObsErrCov) then
	  write(stdout,'(/a)') RC_kt
	  do ikt=1,ktmax
	    write(stdout,'(i2,4(1x,a),1x,7l2)') ikt,
     &	      ktname(ikt),ktunit(ikt),ktdesc(ikt),'$',
     &	      (ktmvar(ikt,jkt),jkt=1,ikt)
	  end do
	  write(stdout,'(a)') '::'
	endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	end
c.
