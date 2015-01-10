!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_diagcorM - forming a diagonal block matrix of given (kr,kt)
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_diagcorM
      implicit none
      private	! except

      public :: diagcorM_solve		! The class data structure
      public :: cg1_solve
      interface diagcorM_solve; module procedure	&
	cg1_solve; end interface

! !REVISION HISTORY:
! 	06Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_diagcorM'

#include "lapack.H"
#include "assert.H"

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: cg1_solve - the good-old conjgr1() solve
!
! !DESCRIPTION:
!
! !INTERFACE:

	subroutine cg1_solve( verbose,kind_cov,		&
     &		handle, kr, kt,					&
     &		nnobs,sigOu,sigOc,kx,ks,qr,kl,		&
     &		ktab,jtab,sigFi,			&
     &		      qrPhi,qmPhi,qlPhi,klPhi,		&
     &		      qrPsi,qmPsi,qlPsi,klPsi,		&
     &		nvecs,y,x, ierr)

!.......................................................................
!.... Argument declarations.
	use m_stdio,only : stderr,stdout
	use m_die  ,only : perr, assert_
	use m_block_storage,only : get_handle, release_handle
	use m_mall ,only : mall_ison,mall_mci,mall_mco

	implicit none

	logical,intent(in) :: verbose
	integer,intent(in) :: kind_cov

	integer,intent(in) :: handle
	integer,intent(in) :: kr
	integer,intent(in) :: kt

	integer,intent(in) :: nnobs

	real   ,intent(in) :: sigOu(nnobs)
	real   ,intent(in) :: sigOc(nnobs)
	integer,intent(in) ::   kx(nnobs)
	integer,intent(in) ::   ks(nnobs)
	real   ,intent(in) :: qr(3,nnobs)
	integer,intent(in) ::   kl(nnobs)

	integer,intent(in) :: ktab(nnobs)
	integer,intent(in) :: jtab(nnobs)
	real   ,intent(in) :: sigFi(nnobs)

	real   ,intent(in) :: qrPhi(3,nnobs)
	real   ,intent(in) :: qmPhi(3,nnobs)
	real   ,intent(in) :: qlPhi(3,nnobs)
	integer,intent(in) :: klPhi(  nnobs)

	real   ,intent(in) :: qrPsi(3,nnobs)
	real   ,intent(in) :: qmPsi(3,nnobs)
	real   ,intent(in) :: qlPsi(3,nnobs)
	integer,intent(in) :: klPsi(  nnobs)

	integer,intent(in) :: nvecs
	   real,intent(in) :: y(nvecs,nnobs)
	   real,intent(out):: x(nvecs,nnobs)

	integer,intent(out) :: ierr

! !REVISION HISTORY:
! 	06Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::cg1_solve'

!.......................................................................
!  Local workspace

	real, pointer :: corrI(:)
	real, pointer :: corrM(:)
	real, Allocatable :: tmp_I(:)
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

	include	  "mxpass.h"
	include     "bands.h"
	include	  "realvals.h"

	real _SDOT,_SNRM2
	external _SDOT,_SNRM2	! BLAS functions

	logical :: computeM, releaseM
	logical :: computeI, releaseI
	character(len=1), pointer :: mtyp
	Integer :: nn_elem, corrI_size

!!$	Integer, Parameter :: FUSE_SIZE = 10000
!!$	Integer, Parameter :: FUSE_SIZE = 24
	Integer, Parameter :: FUSE_SIZE = 1

!.......................................................................

	ASSERT(nnobs >=0)
	ASSERT(nvecs >=0)

	if(nnobs==0) return
	if(nvecs==0) return

	allocate(	xx(nnobs,nvecs), rr(nnobs,nvecs),	&
     &			zz(nnobs,nvecs), pp(nnobs,nvecs),	&
     &			Cpp(nnobs,nvecs),			&
     &			sizerr(0:mxpass,nvecs),			&
     &			rrzznew(nvecs),	stat=ierr )

	if(ierr.ne.0) then
	  call perr(myname_,'allocate()',ierr)
	  return
	endif

	if(mall_ison()) then
	  call mall_mci(xx,myname)
	  call mall_mci(rr,myname)
	  call mall_mci(zz,myname)
	  call mall_mci(pp,myname)
	  call mall_mci(Cpp,myname)
	  call mall_mci(sizerr,myname)
	  call mall_mci(rrzznew,myname)
	endif

!.......................................................................
!.... Initialization.

	xx(:,:) =0.
	rr(:,:) = transpose( y(:,1:nnobs))
	 x(:,:) = transpose(xx(1:nnobs,:))

	 nn_elem = nnobs*(nnobs+1)/2
	 corrM => get_handle(2*handle-1,nn_elem,computeM,releaseM, ierr)
         if(ierr.ne.0) return


	 If (computeM) Then

	    begin_blk = 1
	    begin_sav = 1
	    corrI_size = 0

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
		  if (next)  &
		       next= (qrPhi(1,next_blk).eq.endqrx) .or. (next_blk - begin_blk < FUSE_SIZE)
	       enddo

	       m = next_blk - begin_blk                               
	       corrI_size = corrI_size + m*(m+1)/2
	       begin_blk = next_blk
	       begin_sav = begin_sav+m*(m+1)/2
	    end do	! next block (starting from begin_blk)?

	 End If

	 corrI => get_handle(2*handle, corrI_size, computeI, releaseI, ierr)
         if(ierr.ne.0) return

	 if (computeM) Then
	call diagcorM(kind_cov,			&
     &	  kr,kt,nnobs,sigOu,sigOc,kx,ks,qr,kl,	&
     &		ktab,jtab,sigFi,		&
     &		qrPhi,qmPhi,qlPhi,klPhi,	&
     &		qrPsi,qmPsi,qlPsi,klPsi,	&
     &	  corrM, ierr)

	if(ierr.ne.0) then
	  write(stderr,'(a,3(a,i3))') myname_,		&
     &	    ': unexpected variable type for diagcorM(), (kr,kt) = (',&
     &	    kr,',',kt,'), ierr =',ierr

	  If (releaseM .or. releaseI) Then
             Call release_handle(2*handle-1)
             Call release_handle(2*handle)
          End if

		if(mall_ison()) then
		  call mall_mco(xx,myname)
		  call mall_mco(rr,myname)
		  call mall_mco(zz,myname)
		  call mall_mco(pp,myname)
		  call mall_mco(Cpp,myname)
		  call mall_mco(sizerr,myname)
		  call mall_mco(rrzznew,myname)
		endif
	  deallocate(xx,rr,zz,pp,Cpp,sizerr,rrzznew)

	  return
	endif
     end if


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

	  x(:,:) = transpose(xx(1:nnobs,:))

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
	      if (next)  &
		 next= (qrPhi(1,next_blk).eq.endqrx) .or. (next_blk - begin_blk < FUSE_SIZE)
	    enddo

	    m = next_blk - begin_blk                               

	    if(computeI .and. iter.eq.0) then

	      call smex(corrM,nnobs,begin_blk,m,corrI(begin_sav))

	      call _SPPTRF('U',m,corrI(begin_sav),ierr)
	      if(ierr.ne.0) then
		write(stderr,'(2a,i4)') myname_,	&
     &		  ': SPPTRF() error, stat =',ierr

         	If (releaseM .or. releaseI) Then
                   Call release_handle(2*handle-1)
	           Call release_handle(2*handle)
                End if

			if(mall_ison()) then
			  call mall_mco(xx,myname)
			  call mall_mco(rr,myname)
			  call mall_mco(zz,myname)
			  call mall_mco(pp,myname)
			  call mall_mco(Cpp,myname)
			  call mall_mco(sizerr,myname)
			  call mall_mco(rrzznew,myname)
			endif
	        deallocate(xx,rr,zz,pp,Cpp,sizerr,rrzznew)
		return
	      endif

	    endif	! iter.eq.0

	      call _SPPTRS('U',m,nvecs,corrI(begin_sav),	&
     &		zz(begin_blk,1),nnobs,ierr)

	    if(ierr.ne.0) then	! if it ever happens.
	      write(stderr,'(2a,i3,a,i5,a,i3)') myname_,	&
     &		': SPPTRS(',m,',',begin_blk,') error, stat =',ierr


         	If (releaseM .or. releaseI) Then
                   Call release_handle(2*handle-1)
	           Call release_handle(2*handle)
                End If

			if(mall_ison()) then
			  call mall_mco(xx,myname)
			  call mall_mco(rr,myname)
			  call mall_mco(zz,myname)
			  call mall_mco(pp,myname)
			  call mall_mco(Cpp,myname)
			  call mall_mco(sizerr,myname)
			  call mall_mco(rrzznew,myname)
			endif
	      deallocate(xx,rr,zz,pp,Cpp,sizerr,rrzznew)

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
	    write(stdout,'(2a,3(a,i2))') myname_,	&
     &		': premature exit, ', &
     &		'iter = ',iter-1,', kr = ',kr,', kt = ',kt
	  elseif(.not.cnvrged) then
	    write(stdout,'(2a,3(a,i2))') myname_,	&
     &		': max-pass exceeded, ', &
     &		'iter = ',iter,', kr = ',kr,', kt = ',kt
	  else
	    write(stdout,'(2a,3(a,i2))') myname_,	&
     &		': convergence achieved, ',		&
     &		'iter = ',iter,', kr = ',kr,', kt = ',kt
	  endif

	  call cgnorm(myname_,criter(1,1),mxpass,iter,nvecs,sizerr, &
     &		nnobs)
	endif

        If (releaseM .or. releaseI) Then
           Call release_handle(2*handle-1)
           Call release_handle(2*handle)
        End If

	Nullify(corrM)
	Nullify(corrI)
	
		if(mall_ison()) then
		  call mall_mco(xx,myname)
		  call mall_mco(rr,myname)
		  call mall_mco(zz,myname)
		  call mall_mco(pp,myname)
		  call mall_mco(Cpp,myname)
		  call mall_mco(sizerr,myname)
		  call mall_mco(rrzznew,myname)
		endif
	deallocate(xx,rr,zz,pp,Cpp,sizerr,rrzznew)
	

!***********************************************************************
end subroutine cg1_solve

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !IROUTINE: diagcorM - form a diagonal block matrix of given (kr,kt)
!
! !DESCRIPTION: (to do)
!
! !INTERFACE:

    subroutine diagcorM(kind_cov,		&
	kr,kt,n_x,sigOu,sigOc,kx,ks,qr,kl,	&
		  ktab,jtab,sigFi,		&
		  qrPhi,qmPhi,qlPhi,klPhi,	&
		  qrPsi,qmPsi,qlPsi,klPsi,	&
		corM, istat)

      use config,only : kind_covU,kind_covC
      use config,only : kind_covF,kind_covS,kind_covV
      use config,only : ktuu,ktus,ktvv,ktvs
      use FEsigW_imat ,only : FEsigS_imat,FEsigV_imat
      use m_block_corU,only : diagcor
      use m_block_corO,only : diagcor
      use m_block_corF,only : diagcor,offdcor
      use m_block_corD,only : diagcorD
      use m_stdio,only : stderr
      use m_die,  only : die,perr
      use m_mall, only : mall_ison,mall_mci,mall_mco
      implicit none

  integer, intent(in)	:: kind_cov	! which correlation matrices
  integer, intent(in)	:: kr		! regional index
  integer, intent(in)	:: kt		! data type index
  integer, intent(in)	:: n_x		! size

	! R = R_u + R_c

  real,dimension(n_x),intent(in) :: sigOu ! uncorrelated obs. err. stdv.
  real,dimension(n_x),intent(in) :: sigOc ! correlated obs. err. stdv.

  integer,dimension(  n_x),intent(in) :: kx
  integer,dimension(  n_x),intent(in) :: ks
  real   ,dimension(3,n_x),intent(in) :: qr
  integer,dimension(  n_x),intent(in) :: kl

	! Pf = ...

  integer,dimension(  n_x),intent(in) :: ktab
  integer,dimension(  n_x),intent(in) :: jtab
  real,dimension(n_x),intent(in) :: sigFi ! fcst. err. stdv.

  real   ,dimension(3,n_x),intent(in) :: qrPhi
  real   ,dimension(3,n_x),intent(in) :: qmPhi
  real   ,dimension(3,n_x),intent(in) :: qlPhi
  integer,dimension(  n_x),intent(in) :: klPhi

  real   ,dimension(3,n_x),intent(in) :: qrPsi
  real   ,dimension(3,n_x),intent(in) :: qmPsi
  real   ,dimension(3,n_x),intent(in) :: qlPsi
  integer,dimension(  n_x),intent(in) :: klPsi

  real,dimension(n_x*(n_x+1)/2),intent(out) :: corM ! packed matrix

  integer,intent(out) :: istat	! status

! !REVISION HISTORY:
! 	21Feb96 - J. Guo	- (to do)
! 	 6Jul96	- J. Larson - Replaced addition of variables kind_covF, 
!                 et cetera with intrinsic ior() calls.
!_______________________________________________________________________

	! Local workspace

  real, allocatable	:: corr(:)
  real, allocatable	:: sigW(:)

	! local vars.

  integer		:: i
  character(len=1)	:: Mtyp

  character(len=*),parameter :: myname_=myname//'::diagcorM'
  integer, parameter	:: kind_covX=	&
	ior(ior(kind_covU,kind_covC),	&
	    ior(kind_covF,ior(kind_covS,kind_covV)))
!=======================================================================
	! option checking, switches to be implemented

  if(kind_cov .ne. kind_covX)	&
	call die(myname_,'unexpected kind_cov',kind_cov)
!-----------------------------------------------------------------------

	allocate(corr(n_x*(n_x+1)/2),sigW(n_x), stat=istat)
	if(istat.ne.0) then
	  call perr(myname_,'allocate()',istat)
	  return
	endif
	if(mall_ison()) then
	  call mall_mci(corr,myname)
	  call mall_mci(sigW,myname)
	endif

	!--------------------------------
	! Initialize the output

  istat=0		! a scalar
  corM(:)=0.		! a packed matrix
!-----------------------------------------------------------------------

  Mtyp='Z'
  call diagcor(kind_covC,kt,n_x,kx,qr,kl,	&
		Mtyp,corr,istat)

  	call check_istat_('diagcorO')
  	if(istat/=0) return

  call add_corM_(sigOc)
	if(istat/=0) return

!-----------------------------------------------------------------------
  Mtyp='Z'
  call diagcor(kind_covU,kt,n_x,kx,ks,kl,Mtyp,corr,istat)

  	call check_istat_('diagcorU')
  	if(istat/=0) return

  call add_corM_(sigOu)
	if(istat/=0) return

!-----------------------------------------------------------------------
  select case(kt)
  case (ktuu,ktus,ktvv,ktvs)

    Mtyp='Z'
    call diagcorD(kind_covF,kt,n_x,		&
	ktab,jtab, qrPhi,qmPhi,qlPhi,klPhi,	&
	Mtyp,corr,istat)

  case default

    Mtyp='Z'
    call diagcor(kind_covF,kt,n_x,qrPhi,qmPhi,klPhi,	&
	Mtyp,corr,istat)

  end select

	call check_istat_('diagcorF')
  	if(istat/=0) then
	  return
	endif

		! Calculation in this routing is not always
		! correct, it applies to only if sigFi for both
		! u and v is the same.

  call add_corM_(sigFi)
	if(istat/=0) return

!-----------------------------------------------------------------------
  if(kt/=ktuu .and. kt/=ktvv .and. kt/=ktus .and. kt/=ktvs) then
		if(mall_ison()) then
		  call mall_mco(corr,myname)
		  call mall_mco(sigW,myname)
		endif
	deallocate(corr, sigW)
    return
  endif

	!----------------------------------------
  Mtyp='Z'
  call diagcorD(kind_covS,kt,n_x,		&
	qrPsi,qmPsi,qlPsi,klPsi,		&
	Mtyp,corr,istat)

	call check_istat_('diagcorS')
	if(istat/=0) return

  do i=1,n_x
    sigW(i)=FEsigS_imat(ktab(i),jtab(i))
  end do
  call add_corM_(sigW)
	if(istat/=0) return

	!----------------------------------------
  Mtyp='Z'
  call diagcorD(kind_covV,kt,n_x,		&
	qrPsi,qmPsi,qlPsi,klPsi,		&
	Mtyp,corr,istat)

	call check_istat_('diagcorV')
	if(istat/=0) return

  do i=1,n_x
    sigW(i)=FEsigV_imat(ktab(i),jtab(i))
  end do
  call add_corM_(sigW)
	if(istat/=0) return

!-----------------------------------------------------------------------
		if(mall_ison()) then
		  call mall_mco(corr,myname)
		  call mall_mco(sigW,myname)
		endif
	deallocate(corr, sigW)
!=======================================================================
contains
!-----------------------------------------------------------------------
subroutine check_istat_(name)	! check istat values
use m_stdio,only : stderr
implicit none
	!--------------------------------
  character(len=*), intent(in)	:: name
	!--------------------------------
  if(istat/=0) then
    write(stderr,'(3a,3(a,i3))') myname_,		&
      ': unexpected variable type for ',name,		&
      '(), (kr,kt) = (',kr,',',kt,'), istat =',istat

		if(mall_ison()) then
		  call mall_mco(corr,myname)
		  call mall_mco(sigW,myname)
		endif
    deallocate(corr, sigW)
    return
  endif

#ifdef	_DEBUG
	_DEBUG	write(*,'(3a,i2,a)') myname_,	&
	_DEBUG	  ': diagcorM(',Mtyp,kt,')'
#endif

end subroutine check_istat_
!-----------------------------------------------------------------------
subroutine add_corM_(sigm)
  use m_stdio,only : stderr

	!--------------------------------
	! add a covariance matrix to corM()

  implicit none
  real, intent(in)	:: sigm(n_x)

	! Locals
  integer i,j,ij

#ifdef	_DEBUG
_DEBUG	write(*,*) n_x
_DEBUG	do i=1,min(3,n_x)
_DEBUG	  write(*,'(i2,x,f6.2,x,i3,x,i7,x,i3,x,i3,$)') i,	&
_DEBUG		sigm(i),kx(i),ks(i),ktab(i),jtab(i)
_DEBUG	  do j=1,min(3,n_x)
_DEBUG	    ij=min(i,j)+(max(i,j)-1)*max(i,j)/2
_DEBUG	    write(*,'(1x,f5.3,$)') corr(ij)
_DEBUG	  end do
_DEBUG	  write(*,'(x,$)')
_DEBUG	  do j=1,min(3,n_x)
_DEBUG	    ij=min(i,j)+(max(i,j)-1)*max(i,j)/2
_DEBUG	    write(*,'(1x,1p,e9.3,$)') corM(ij)
_DEBUG	  end do
_DEBUG	  write(*,*)
_DEBUG	end do
#endif
	!--------------------------------
  select case(Mtyp)
  case ('U','u')		! a regular packed matrix
    call add_smat_(sigm)

  case ('I','i')		! an identity matrix
    call add_iden_(sigm)

  case ('D','d')		! a diagonal matrix
    call add_diag_(sigm)

  case ('Z','z')		! a zero matrix
	! do nothing with a zero matrix

  case default			! likely an error
    istat=-1
    write(stderr,'(4a)') myname_,	&
      ': unexpected matrix, Mtyp=(',Mtyp,')'

		if(mall_ison()) then
		  call mall_mco(corr,myname)
		  call mall_mco(sigW,myname)
		endif
    deallocate(corr, sigW)
    return
  end select
end subroutine add_corM_
!-----------------------------------------------------------------------
subroutine add_smat_(sigm)

	! add a packed symmetric covariance matrix to corM()

  implicit none
  real, intent(in)	:: sigm(n_x)

	! Locals
  integer	:: i,j,ij

  do j=1,n_x		! for all column
    ij=j*(j-1)/2	! where the last column end
    do i=1,j
      corM(ij+i)=sigm(i)*corr(ij+i)*sigm(j) + corM(ij+i)
    end do
  end do
end subroutine add_smat_
!-----------------------------------------------------------------------
subroutine add_diag_(sigm)

	! add a diagonal covariance matrix to corM()

  implicit none
  real, intent(in)	:: sigm(n_x)

	! Locals
  integer	:: j,ij

  do j=1,n_x		! for all column
    ij=j*(j+1)/2	! the diagonal element in the matrix only

    corM(ij)=sigm(j)*corr(ij)*sigm(j) + corM(ij)
  end do
end subroutine add_diag_
!-----------------------------------------------------------------------
subroutine add_iden_(sigm)

	! add a covariance matrix of identity correlation to corM()

  implicit none
  real, intent(in)	:: sigm(n_x)

	! Locals
  integer	:: j,ij

  do j=1,n_x		! for all column
    ij=j*(j+1)/2	! the diagonal element in the matrix only

    corM(ij)=sigm(j)*sigm(j) + corM(ij)
  end do
end subroutine add_iden_
!-----------------------------------------------------------------------
#ifdef	_SAVE
_SAVE	subroutine savevectI_(name,kr,kt,nx,ivec)
_SAVE	  character(len=*), intent(in)	:: name
_SAVE	  integer,	    intent(in)	:: kr,kt,nx
_SAVE	  integer,	    intent(in)	:: ivec(nx)
_SAVE	
_SAVE	  if(kr.eq.4.and.kt.eq.ktHH.and.nx.gt.0) then
_SAVE	    call savevect(name,kr,kt,nx,ivec)
_SAVE	  endif
_SAVE	end subroutine savevectI_
!-----------------------------------------------------------------------
_SAVE	subroutine savevectR_(name,kr,kt,nx,xvec)
_SAVE	  character(len=*), intent(in)	:: name
_SAVE	  integer,	    intent(in)	:: kr,kt,nx
_SAVE	  real,		    intent(in)	:: xvec(nx)
_SAVE	
_SAVE	  if(kr.eq.4.and.kt.eq.ktuu.and.nx.gt.0) then
_SAVE	    call savevect(name,kr,kt,nx,xvec)
_SAVE	  endif
_SAVE	end subroutine savevectR_
!-----------------------------------------------------------------------
_SAVE	subroutine savemat1_(name,kr,kt,nx,corM)
_SAVE	  character(len=*), intent(in)	:: name
_SAVE	  integer,	    intent(in)	:: kr,kt,nx
_SAVE	  real,		    intent(in)	:: corM(nx*(nx+1)/2)
_SAVE	
_SAVE	  if(kr.eq.4.and.kt.eq.ktuu.and.nx.gt.0) then
_SAVE	    call savemat1(name,kr,kt,nx,corM)
_SAVE	  endif
_SAVE	end subroutine savemat1_
!-----------------------------------------------------------------------
_SAVE	subroutine savecorrd_(kr,kt,nx,qr_x,qr_y,qr_z)
_SAVE	  implicit none
_SAVE	  integer, intent(in) :: kr,kt,nx
_SAVE	  real,    intent(in) :: qr_x(nx),qr_y(nx),qr_z(nx)
_SAVE
_SAVE	  integer ierr,i
_SAVE	  real coslat
_SAVE	  real,allocatable :: rlon(:),rlat(:)
_SAVE
_SAVE	  allocate(rlat(nx),rlon(nx),stat=ierr)
_SAVE	  if(ierr.ne.0) then
_SAVE	    write(stderr,'(2a,$)') myname_,	&
_SAVE		': savecorrd_()/allocate() error, stat = '
_SAVE	    write(stderr,*) ierr
_SAVE	    call die(myname_)
_SAVE	  endif
_SAVE	  if(kr.eq.4.and.kt.eq.ktuu.and.nx.gt.0) then
_SAVE	    do i=1,nx
_SAVE	      coslat=sqrt(qr_x(i)*qr_x(i)+qr_y(i)*qr_y(i))
_SAVE
_SAVE	      rlon(i)=modulo(atan2(qr_y(i)/coslat,qr_x(i)/coslat),360.)
_SAVE	      rlat(i)=atan2(qr_z(i),coslat)
_SAVE	    end do
_SAVE	  endif
_SAVE
_SAVE	  call savevect('rlat',kr,kt,nx,rlat)
_SAVE	  call savevect('rlon',kr,kt,nx,rlon)
_SAVE	  deallocate(rlat,rlon)
_SAVE	end subroutine savecorrd_
#endif
!-----------------------------------------------------------------------
end subroutine diagcorM
end module m_diagcorM
!.
