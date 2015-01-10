	subroutine cgnorm(CGname,criter,mxpass,npass,nvecs,sizerr,ndat)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: cgnorm - print a table of CG solution residuals
!
! !SYNOPSIS: (to do)
! !INPUT PARAMETERS: (to do)
! !OUTPUT PARAMETERS: (to do)
!
! !DESCRIPTION:
!	print a table of saved CG solution residuals.
!
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
!
! !REVISION HISTORY:
! 	25Jan96 - J. Guo	- prolog added to a previously written
!				  subroutine
!_______________________________________________________________________
	use m_mpout,only : mpout_ison
	use m_mpout,only : mpout
	implicit none

!!interface

	character*(*) CGname	! name of the routin
	real criter		! convergence criterion
	integer mxpass		! the first mension
	integer npass		! number of values may have stored
	integer nvecs		! number of vectors
	real sizerr(0:mxpass,nvecs)	! see conjgr() for the meaning
	integer ndat		! local size of a vector

!!end_interface
	include "realvals.h"

c   ..Local variables
c  ===================
	real scale		! the order of magnitude of sizerr
	real pwr		! log10(scale)
	integer ipwr		! the power of scale
	integer ivec		! vector index
	integer k,kn		! pass index

c   ..Local parameters
c   ================-=

	integer mxprt
	parameter(mxprt=16)	! total number of values to list

	integer lnblnk,ln	! string size
	external lnblnk

c   ===========================================
	if(.not.mpout_ison()) return

c   ..Check the magnitudes of all sizerr values
c   ===========================================
	scale=rsmlval
	do ivec=1,nvecs
	  do k=0,min(mxpass,npass)
	    if(scale.lt.sizerr(k,ivec)) scale=sizerr(k,ivec)
	  end do
	end do

c   ..Extract the power of the magnitude and the scale
c   ==================================================
	pwr=log(scale)/log(10.)
	ipwr=int(pwr)
	if(ipwr.gt.pwr) ipwr=ipwr-1
	pwr=ipwr
	scale=10.**(pwr)

c   ..Write the header (label)
c   ==========================
	ln=max(lnblnk(CGname),1)
	write(mpout,'(2a)') CGname(1:ln),'::'

c   ..Write related values
c   ======================
	write(mpout,'(t5,a,t17,sp,i10)')	'Scale',	ipwr
	write(mpout,'(t5,a,t17,1p,e10.3)')	'Criterion',	criter
	write(mpout,'(t5,a,t17,i10)')		'localSize',	ndat
	write(mpout,'(t5,a,t17,i10)')		'Iterations',	npass
	if(npass.gt.mxpass)
     &	  write(mpout,'(a,i3,a)')		'# the last',	mxpass,
     &	    ' iterations only'
	write(mpout,'(t5,a,t17,i3)')		'Vectors',	nvecs
	if(nvecs.gt.mxprt)
     &	  write(mpout,'(a,i3,a)')		'# the first',	mxprt,
     &	    ' vectors only'

c   ..Write the table
c   =================
	write(mpout,'(i3,16(f8.4))') 0,
     &	  (sizerr(0,ivec)/scale,ivec=1,min(nvecs,mxprt))

c   ..List only the last mxpass entris, if npass is greater than mxpass
c   ===================================================================
	do k=max(1,npass-mxpass+1),npass
	  kn=mod(k-1,mxpass)+1
	  write(mpout,'(i3,16(f8.4))') k,
     &	    (sizerr(kn,ivec)/scale,ivec=1,min(nvecs,mxprt))
	end do

c   ..List the values determined by criter
c   ======================================
	write(mpout,'(t2,a,16(f8.4))') '**',
     &	  (criter*sizerr(0,ivec)/scale,ivec=1,min(nvecs,mxprt))

c   ..Write the end of table mark
c   =============================
	write(mpout,'(2a)') '::      # ',CGname(1:ln)

c=======================================================================
	end
c.
