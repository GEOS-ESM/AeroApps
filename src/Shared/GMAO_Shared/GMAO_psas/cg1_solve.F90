
	subroutine cg1_solve( verbose,kind_cov,			&
     &		kr, kt,						&
     &		nnobs,  kx,  ks,ktab,jtab,sigU,sigO,sigF,	&
     &		      qrPhi,qmPhi,qlPhi,klPhi,			&
     &		      qrOth,qmOth,qlOth,klOth,			&
     &		nvecs,						&
     &		     ldy,  y,					&
     &		     ldx,  x,       ierr  )

!.......................................................................
!.... Argument declarations.
	use m_stdio,only : stderr,stdout
	use m_diagcorM,only : diagcorM

	implicit none

	logical,intent(in) :: verbose
	integer,intent(in) :: kind_cov

	integer,intent(in) :: kr
	integer,intent(in) :: kt

	integer,intent(in) :: nnobs
	integer,intent(in) ::   kx(nnobs)
	integer,intent(in) ::   ks(nnobs)
	integer,intent(in) :: ktab(nnobs)
	integer,intent(in) :: jtab(nnobs)

	real   ,intent(in) :: sigU(nnobs)
	real   ,intent(in) :: sigO(nnobs)
	real   ,intent(in) :: sigF(nnobs)

	real   ,intent(in) :: qrPhi(3,nnobs)
	real   ,intent(in) :: qmPhi(3,nnobs)
	real   ,intent(in) :: qlPhi(3,nnobs)
	integer,intent(in) :: klPhi(  nnobs)

	real   ,intent(in) :: qrOth(3,nnobs)
	real   ,intent(in) :: qmOth(3,nnobs)
	real   ,intent(in) :: qlOth(3,nnobs)
	integer,intent(in) :: klOth(  nnobs)


	integer,intent(in) :: nvecs
	integer,intent(in) :: ldy
	   real,intent(in) :: y(ldy,nvecs)
	integer,intent(in) :: ldx
	   real,intent(in) :: x(ldx,nvecs)

	integer,intent(out) :: ierr

!.......................................................................
!  Local workspace

	real, allocatable :: corrI(:),corrM(:)
	real, allocatable :: xx(:,:),rr(:,:),zz(:,:),pp(:,:),Cpp(:,:)
	real, allocatable :: sizerr(:,:),rrzznew(:)


!.......................................................................
!.... Local storage.

	logical      time_to_quit
	integer      begin_blk, next_blk, begin_sav
	real         endqrx
	logical	next
	integer ivec

	integer ij

	real	rrzzold,ppCpp
	real	alphap,alpham,beta
	integer	m,i,j,kn,km,iter
	logical cnvrged,cnvrging,passmin

!.......................................................................
!.... Convergence control parameters.

#include "lapack.H"

	include	  "mxpass.h"
	include     "bands.h"
	include	  "realvals.h"

	character*7  myname
	parameter   (myname='cg1_solve')

	real _SDOT,_SNRM2
	external _SDOT,_SNRM2	! BLAS functions

	integer lnblnk,luavail
	external lnblnk,luavail

!.......................................................................

	allocate(	xx(nnobs,nvecs), rr(nnobs,nvecs),	&
     &			zz(nnobs,nvecs), pp(nnobs,nvecs),	&
     &			Cpp(nnobs,nvecs),			&
     &			sizerr(0:mxpass,nvecs),			&
     &			rrzznew(nvecs),				&
     &			corrM(nnobs*(nnobs+1)/2),		&
     &			corrI(nnobs*(nnobs+1)/2), stat=ierr )

	if(ierr.ne.0) then
	  write(stderr,'(2a,i3)') myname,	&
     &	    ': allocate() error, stat =',ierr
	  return
	endif

!.......................................................................
!.... Initialization.

	do ivec=1,nvecs
	  do i=1,nnobs
	    xx(i,ivec)=0.
	  end do
	end do

	do ivec=1,nvecs
	  call _SCOPY(nnobs,y(1,ivec),1,rr(1,ivec),1)	! rr=y
	  call _SCOPY(nnobs,xx(1,ivec),1,x(1,ivec),1)	! x=xx
	end do

	call diagcorM(kind_cov,				&
     &	  kr,kt,nnobs,kx,ks,ktab,jtab,sigU,sigO,sigF,	&
     &	  qrPhi,qmPhi,qlPhi,klPhi,			&
     &	  qrOth,qmOth,qlOth,klOth,			&
     &	  corrM, ierr)

	if(ierr.ne.0) then
	  write(stderr,'(a,3(a,i3))') myname,		&
     &	    ': unexpected variable type for diagcorM(), (kr,kt) = (',&
     &	    kr,',',kt,'), ierr =',ierr

	  deallocate(xx,rr,zz,pp,Cpp,sizerr,rrzznew,corrM,corrI)
	  return
	endif

!.......................................................................
!.... Start iterations.

	iter = 0
	kn=0

	CONJGR1_LOOP: do

!.......................................................................
!.... Test for convergence.

!	..iter is the actual iteration index for the norm.  The initial
!	norm is stored in sizerr(0,ivec).  The later norms are stored in
!	sizerr(1:mxpass,ivec) in a cyclic list (kn=1..mxpass,1..mxpass).
!	Kn is where the iter-th entry is in sizerr() list.  Km is the
!	previous kn entry.

		! Creating a circulating pointer kn

	  km=kn
	  if(iter.eq.0) then
	    kn=0
	  else
	    kn=mod(iter-1,mxpass)+1
	  endif

		! compute current norm
	  do ivec=1,nvecs
	    sizerr(kn,ivec) = _SNRM2(nnobs,rr(1,ivec),1)
	  end do

		! are the solutions all converging?
	  cnvrging=.true.
	  ivec=0
	  do while(iter.gt.0 .and. cnvrging .and. ivec.lt.nvecs)
	    ivec=ivec+1
	    cnvrging = sizerr(kn,ivec).le.sizerr(km,ivec)
	  end do

		! have the solutions converged?
	  cnvrged=.true.
	  ivec=0
	  do while(cnvrged .and. ivec.lt.nvecs)
	    ivec=ivec+1
	    cnvrged = sizerr(kn,ivec).le.sizerr(0,ivec)*criter(1,1)
	  end do

		! have the solutions pass the threashold value

	  passmin=iter.gt.0
	  ivec=0
	  do while(passmin .and. ivec.lt.nvecs)
	    ivec=ivec+1
	    passmin = sizerr(km,ivec).le.sizerr(0,ivec)*criter(1,2)
	  end do

		! may be it is the time to quit?

	  time_to_quit = cnvrged
	  if(.not.time_to_quit) time_to_quit = iter.gt.maxpass(1)
	  if(.not.time_to_quit) time_to_quit =	&
     &	    (.not. cnvrging) .and. (iter.gt.minpass(1) .or. passmin)


		! return with the saved solution

	  if(time_to_quit.and..not.cnvrging) exit CONJGR1_LOOP

		! Save this solution before the next iteration

	  do ivec=1,nvecs
	    call _SCOPY(nnobs,xx(1,ivec),1,x(1,ivec),1)	! x=xx
	  end do

	  if(time_to_quit) exit CONJGR1_LOOP

!.......... Preconditioner for level 1 (one region, one kt) is direct 
!.......... solver on diagonal sub-blocks.
!.......... Try to keep soundings together (group by qr_x)

	  do ivec=1,nvecs
	    call _SCOPY(nnobs,rr(1,ivec),1,zz(1,ivec),1)	! zz=rr
	  end do

	  begin_blk = 1
	  begin_sav = 1

	  do while ( begin_blk .le. nnobs )           

!	    ..It (next_blk) is actually the end-of-this-block here.

	    next_blk = min(begin_blk+msmall-1,nnobs)
	    endqrx = qrPhi(1,next_blk)

!.......................................................................
!.......... search for end of this sounding (at end of msmall sized 
!.......... block) and set block break where soundings change

!	    ..Tests are made in sequence to avoid qr_x(nnobs+1) ever
!	    being referenced.

	    next=.true.
	    do while ( next )
	      next_blk = next_blk + 1

	      next=next_blk.le.nnobs
	      if(next) next=qrPhi(1,next_blk).eq.endqrx
	    enddo

	    m = next_blk - begin_blk                               

	    if(iter.eq.0) then

	      call smex(corrM,nnobs,begin_blk,m,corrI(begin_sav))

	      call _SPPTRF('U',m,corrI(begin_sav),ierr)
	      if(ierr.ne.0) then
		write(stderr,'(2a,i4)') myname,		&
     &		  ': SPPTRF() error, stat =',ierr

	        deallocate(xx,rr,zz,pp,Cpp,sizerr,rrzznew,corrM,corrI)
		return
	      endif

	    endif	! iter.eq.0

	    do ivec=1,nvecs
	      call _SPPTRS('U',m,1,corrI(begin_sav),	&
     &		zz(begin_blk,ivec),nnobs,ierr)
	    end do

	    if(ierr.ne.0) then	! if it ever happens.
	      write(stderr,'(2a,i3,a,i5,a,i3))') myname,	&
     &		': SPPTRS(',m,',',begin_blk,') error, stat =',ierr

	      deallocate(xx,rr,zz,pp,Cpp,sizerr,rrzznew,corrM,corrI)
	      return
	    endif

	    begin_blk = next_blk
	    begin_sav = begin_sav+m*(m+1)/2

	  end do	! next block (starting from begin_blk)?

!.......................................................................
!.......... Set search direction, pp.

	  if( iter.eq.0 ) then
	    do ivec=1,nvecs
	      rrzznew(ivec) = _SDOT(nnobs,rr(1,ivec),1,zz(1,ivec),1)
	      call _SCOPY(nnobs,zz(1,ivec),1,pp(1,ivec),1)
	    end do

	  else

	    do ivec=1,nvecs
	      rrzzold = rrzznew(ivec)
	      rrzznew(ivec) = _SDOT(nnobs,rr(1,ivec),1,zz(1,ivec),1)
	      beta = rrzznew(ivec) / rrzzold
	      call _SAXPY(nnobs,beta,pp(1,ivec),1,zz(1,ivec),1)
	      call _SCOPY(nnobs,zz(1,ivec),1,pp(1,ivec),1)
	    end do

	  endif

!.......................................................................
!.......... Update solution estimate, xx, and residual, rr.

	  do ivec=1,nvecs
	    call _SSPMV('U',nnobs, 1.,corrM,pp(1,ivec),1,	&
     &		                  0.,     Cpp(1,ivec),1)

	    ppCpp=_SDOT(nnobs,pp(1,ivec),1,Cpp(1,ivec),1)

	    alphap = rrzznew(ivec) / ppCpp
	    alpham = -alphap
	    call _SAXPY(nnobs,alphap, pp(1,ivec),1,xx(1,ivec),1)
	    call _SAXPY(nnobs,alpham,Cpp(1,ivec),1,rr(1,ivec),1)
	  end do

	  iter=iter+1

	end do CONJGR1_LOOP

!.......................................................................
!....... time to quit, give some messages

	
	if(verbose) then

	  if(.not.cnvrging) then
	    write(stdout,'(2a,3(a,i2))') myname,': premature exit, ', &
     &		'iter = ',iter-1,', kr = ',kr,', kt = ',kt
	  elseif(.not.cnvrged) then
	    write(stdout,'(2a,3(a,i2))') myname,': max-pass exceeded, ', &
     &		'iter = ',iter,', kr = ',kr,', kt = ',kt
	  else
	    write(stdout,'(2a,3(a,i2))') myname,	&
     &		': convergence achieved, ',		&
     &		'iter = ',iter,', kr = ',kr,', kt = ',kt
	  endif

	  call cgnorm(myname,criter(1,1),mxpass,iter,nvecs,sizerr, &
     &		nnobs)
	endif

	deallocate(xx,rr,zz,pp,Cpp,sizerr,rrzznew,corrM,corrI)
!***********************************************************************
	end subroutine cg1_solve
