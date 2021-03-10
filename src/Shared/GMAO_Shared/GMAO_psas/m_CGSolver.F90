!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_CGSolver - Conjugate Gradient solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_CGSolver
      use m_InnovCovMatx,only : InnovCovMatx
      implicit none
      private	! except

      public :: CGSolver
      public :: CGSolver_init
      public :: CGSolver_solve
      public :: CGSolver_clean
      public :: nband
      public :: nncov

      interface CGSolver_init; module procedure init_; end interface
      interface CGSolver_solve; module procedure	&
	recurSolve_; end interface
      interface CGSolver_clean; module procedure clean_; end interface
      interface nband; module procedure nband_; end interface
      interface nncov; module procedure nncov_; end interface

      type CGSolver
	private
	integer :: iband
        type(InnovCovMatx) :: incov
	type(CGSolver),pointer :: next
      end type CGSolver

! !REVISION HISTORY:
! 	28Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_CGSolver'

#include "assert.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(pcg,corObs,corPhi,corPsi,root,comm,rsrc)
      use m_die,    only : die
      use m_psasrc, only : psasrc_open,psasrc_close
      use m_InnovCovMatx,only : InnovCovMatx_init
      use m_CorAttrF,only : CorAttrF
      use m_CorAttrX,only : CorAttrX
      use m_block_storage, only : block_storage_clean => clean, reserve_storage

      implicit none
      type(CGSolver),intent(out) :: pcg
      type(CorAttrX),intent(in) :: corObs	! C^o_c & C^o_u
      type(CorAttrF),intent(in) :: corPhi	! C^phi
      type(CorAttrF),intent(in) :: corPsi	! C^psi, C^chi
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	28Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  include "kind_mats.h"

  logical :: success
  integer :: ier

  call psasrc_open(root,comm,rsrc=rsrc,stat=ier)
	if(ier/=0) call die(myname_,'psasrc_open()',ier)

  call bands0()

  call psasrc_close(ier)
	if(ier/=0) call die(myname_,'psasrc_close()',ier)

  call recurPush_(pcg,(/kind_5mat,kind_Rmat/),	&
	corObs,corPhi,corPsi,root,comm)

  success = reserve_storage(kind_Umat,10**7)

end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nband_ - top band number
!
! !DESCRIPTION:
!
! !INTERFACE:

    function nband_(pcg)
      implicit none
      type(CGSolver),intent(in) :: pcg
      integer :: nband_

! !REVISION HISTORY:
! 	30May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nband_'
  nband_ = pcg%iband

end function nband_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nncov - return contents of incov
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine nncov_(pcg,nncov)
      implicit none
      type(CGSolver),intent(in) :: pcg
      type(InnovCovMatx),intent(out) :: nncov

! !REVISION HISTORY:
! 	31Aug01	- Todling, initial code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nncov_'
  nncov = pcg%incov

end subroutine nncov_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recurPush_ - recursively specify precoditioners
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine recurPush_(pcg,ibands,		&
	corObs,corPhi,corPsi,root,comm)
      use m_mall,only : mall_ison,mall_ci
      use m_die ,only : die
      use m_mpout,only : mpout
      use m_CorAttrF,only : CorAttrF
      use m_CorAttrX,only : CorAttrX
      use m_InnovCovMatx,only : InnovCovMatx_init
      use m_InnovCovMatx,only : InnovCovMatx_showCosts
      implicit none
      type(CGSolver),intent(out) :: pcg
      integer,dimension(:),intent(in) :: ibands
      type(CorAttrX),intent(in) :: corObs
      type(CorAttrF),intent(in) :: corPhi
      type(CorAttrF),intent(in) :: corPsi
      integer,intent(in) :: root
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recurPush_'
  integer :: ier

  if(size(ibands)==0) call die(myname_,'invalid argument')

  pcg%iband=ibands(1)

  call InnovCovMatx_init(pcg%incov,		&
	corObs,corPhi,corPsi,comm,pcg%iband)

	call InnovCovMatx_showCosts(pcg%incov,mpout,root,comm)

  if(size(ibands)==1) then
    nullify(pcg%next)

  else

    allocate(pcg%next,stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)
	if(mall_ison()) call mall_ci(1,myname)

    call recurPush_(pcg%next,ibands(2:),	&
	corObs,corPhi,corPsi,root,comm)
  endif

end subroutine recurPush_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(pcg)
      use m_InnovCovMatx,only : InnovCovMatx_clean
      use m_die,only : die
      use m_block_storage, only : block_storage_clean => clean

      implicit none
      type(CGSolver),intent(inout) :: pcg

! !REVISION HISTORY:
! 	18May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  call recurPop_(pcg)
  call block_storage_clean()

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recurPop_ - recursively remove precondioners
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine recurPop_(pcg)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_co
      use m_InnovCovMatx,only : InnovCovMatx_clean
      implicit none
      type(CGSolver),intent(inout) :: pcg

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recurPop_'
  integer :: ier

  pcg%iband=-1
  call InnovCovMatx_clean(pcg%incov)

  if(associated(pcg%next)) then
    call recurPop_(pcg%next)
	if(mall_ison()) call mall_co(1,myname)
    deallocate(pcg%next,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)
  endif

end subroutine recurPop_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recurSolve_ - a recursive solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine recurSolve_(pcg,comm,	&
	obNav,krNav,ktNav,			&
	nvi,sigU,sigC,corObs,			&
	ktab,jtab,sigPhi,corPhi,corPsi,		&
	y,x,ierr)

      use m_parDOT,only : parDOT	! SDOT()
      use m_parDOT,only : parNRM2	! SNRM2()
      use m_die   ,only : perr,die,assert_
      use m_mpout ,only : mpout,mpout_ison
      use m_mall  ,only : mall_ison,mall_mci,mall_mco

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize

      use m_InnovCovMatx,only : InnovCovMatx_Cx

      use m_sigmaPhi,only : sigmaPhi
      use m_CorAttrF,only : CorAttrF
      use m_CorAttrX,only : CorAttrX
      use m_zeit, only : zeit_ci, zeit_co

      implicit none

      type(CGSolver),intent(inout) :: pcg
      integer ,intent(in) :: comm

      type(Navigator),intent(in) :: obNav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer ,intent(in) :: nvi

      real   ,dimension(:),intent(in) :: sigU
      real   ,dimension(:),intent(in) :: sigC
      type(CorAttrX),intent(in) :: corObs

      integer,dimension(:),intent(in) :: ktab
      integer,dimension(:),intent(in) :: jtab

      type(sigmaPhi),intent(in) :: sigPhi
      type(CorAttrF),intent(in) :: corPhi
      type(CorAttrF),intent(in) :: corPsi

	! Not (!) defined as assumed-shape

      real   ,dimension(:,:),intent(in)  :: y	! (nvi,nvecs)
      real   ,dimension(:,:),intent(out) :: x	! (nvi,nvecs)
      integer,intent(out) :: ierr

! !REVISION HISTORY:
! 	18May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recurSolve_'
  integer :: ibandcg

#include "lapack.H"

!.... Convergence control parameters.

  include "mxpass.h"
  include "bands.h"

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!.......................................................................

  real,allocatable,dimension(:,:) :: xx,rr,zz,pp,Cpp
  real,allocatable,dimension(:,:) :: sizerr
  real,allocatable,dimension(:)   :: rrzznew,rrzzold

!  Local vars.

	logical      time_to_quit
	real alphap,alpham,beta
	integer iter,km,kn,ivec
	integer kr,ln,kt
	integer i
	logical cnvrged,cnvrging,passmin
	integer :: nvecs,nvseg
!________________________________________

	ASSERT( nvi <= size(ktab) )
	ASSERT( nvi <= size(jtab) )
	ASSERT( nvi <= size(sigU) )
	ASSERT( nvi <= size(sigC) )

	ASSERT( nvi <= size(y,1) )
	ASSERT( nvi == size(x,1) )	! x is the output
	ASSERT( nvi >= 0 )

  nvecs=size(y,2)
	ASSERT( nvecs == size(x,2) )

  nvseg=lsize(obNav)
	ASSERT( nvseg <= size(krNav,1) )
	ASSERT( nvseg <= size(ktNav,1) )
!________________________________________

	ierr=0
	if(nvecs<=0) return
!________________________________________

  ibandcg=pcg%iband

	allocate( xx(nvi,nvecs),rr(nvi,nvecs),	&
		  zz(nvi,nvecs),pp(nvi,nvecs),	&
		  Cpp(nvi,nvecs),			&
		  sizerr(0:mxpass,nvecs),		&
		  rrzznew(nvecs),rrzzold(nvecs),	&
		  stat=ierr		)

		if(ierr/=0) then
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
		  call mall_mci(rrzzold,myname)
		endif

!.......................................................................
!.... Initialization.

	xx(:,:) =0.

	rr(:,:) = y(1:nvi,1:nvecs)
	 x(:,1:nvecs) =xx(:,:)

!.......................................................................
!.... Start next iteration.

	iter = 0
	kn=0

	CONJGR_LOOP: do

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

		! compute current norm globally

	  sizerr(kn,1:nvecs)= parNRM2(rr(1:nvi,1:nvecs),comm)

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
	    cnvrged=sizerr(kn,ivec).le.sizerr(0,ivec)*criter(ibandcg,1)
	  end do

		! have the solutions pass the threashold value

	  passmin=iter.gt.0
	  ivec=0
	  do while(passmin .and. ivec.lt.nvecs)
	    ivec=ivec+1
	    passmin=sizerr(km,ivec).le.sizerr(0,ivec)*criter(ibandcg,2)
	  end do

		! may be it is the time to quit?

	  time_to_quit = cnvrged
	  if(.not.time_to_quit) time_to_quit = iter.gt.maxpass(ibandcg)
	  if(.not.time_to_quit) time_to_quit =	&
	    (.not. cnvrging) .and. (iter.gt.minpass(ibandcg).or.passmin)

		! return with the saved solution
	  if(time_to_quit.and..not.cnvrging) exit CONJGR_LOOP

		! Save this solution before the next iteration

	  x(1:nvi,1:nvecs) = xx(:,:)

	  if(time_to_quit) exit CONJGR_LOOP

!.......................................................................
!....... Test for maximum iterations.

!         Preconditioner for level ibandcg is conjugate gradient on 
!         level 2 done at regional diagonal blocks.
!         ---------------------------------------------------

	  if(.not.associated(pcg%next)) then

		! This is a local solver.  Note that the without
		! finally remove the direct solver as the low level
		! preconditioner, corPhi can not be really implemented.
		! (What I meant is that forming block matrices is a
		! headache.
	     call zeit_ci('localSolve_')
	    call localSolve_(obNav,krNav,ktNav,		&
		nvi,sigU,sigC,corObs,			&
		ktab,jtab,sigPhi,corPhi,corPsi,		&
		rr,zz,ierr)

		if(ierr/=0) then
		  call perr(myname_,'localSolve_()',ierr)
			if(mall_ison()) then
			  call mall_mco(xx,myname)
			  call mall_mco(rr,myname)
			  call mall_mco(zz,myname)
			  call mall_mco(pp,myname)
			  call mall_mco(Cpp,myname)
			  call mall_mco(sizerr,myname)
			  call mall_mco(rrzznew,myname)
			  call mall_mco(rrzzold,myname)
			endif
		  deallocate(xx,rr,zz,pp,Cpp,sizerr,rrzznew,rrzzold)
		  return
		endif
	     call zeit_co('localSolve_')

	  else

		! This is a global solver

	    call recurSolve_(pcg%next,comm,		&
		obNav,krNav,ktNav,			&
		nvi,sigU,sigC,corObs,			&
		ktab,jtab,sigPhi,corPhi,corPsi,		&
		rr,zz,ierr)

		if(ierr/=0) then
		  call perr(myname_,'recurSolve_()',ierr)
			if(mall_ison()) then
			  call mall_mco(xx,myname)
			  call mall_mco(rr,myname)
			  call mall_mco(zz,myname)
			  call mall_mco(pp,myname)
			  call mall_mco(Cpp,myname)
			  call mall_mco(sizerr,myname)
			  call mall_mco(rrzznew,myname)
			  call mall_mco(rrzzold,myname)
			endif
		  deallocate(xx,rr,zz,pp,Cpp,sizerr,rrzznew,rrzzold)
		  return
		endif

	  endif

!.......................................................................
!.......... Set search direction, pp.

	  if(iter==0) then

	    rrzznew(1:nvecs) = parDOT(rr(1:nvi,1:nvecs),	&
				      zz(1:nvi,1:nvecs), comm	)
	    pp(:,:) = zz(:,:)

	  else

	    rrzzold(1:nvecs) = rrzznew(1:nvecs)
	    rrzznew(1:nvecs) = parDOT(rr(1:nvi,1:nvecs),	&
				      zz(1:nvi,1:nvecs), comm	)

	    if(nvi>0) then
	      do ivec=1,nvecs
	        beta = rrzznew(ivec)/rrzzold(ivec)

	        call _SAXPY(nvi,beta,pp(1,ivec),1,zz(1,ivec),1)
	        call _SCOPY(nvi,zz(1,ivec),1,pp(1,ivec),1)
	      end do
	    endif
	  endif

!.......................................................................
!.......... Update solution estimate, xx, and residual, rr.

  call InnovCovMatx_Cx(pcg%incov,comm,	&
	obNav,krNav,ktNav,		&
	nvi,sigU,sigC,corObs,		&
	ktab,jtab,sigPhi,corPhi,corPsi,	&
	pp,Cpp)

	  rrzzold(1:nvecs) = parDOT (pp(1:nvi,1:nvecs),	&
				    Cpp(1:nvi,1:nvecs), comm	)
	  if (nvi>0) then
	    do ivec=1,nvecs
	      alphap = rrzznew(ivec) / rrzzold(ivec)
	      alpham = -alphap

	      call _SAXPY(nvi,alphap, pp(1,ivec),1,xx(1,ivec),1)
	      call _SAXPY(nvi,alpham,Cpp(1,ivec),1,rr(1,ivec),1)
	    end do
	  endif

	  iter=iter+1

	end do CONJGR_LOOP

!.......................................................................
!....... time to quit, give some messages

	if(cgverb(ibandcg)) then

	  if(mpout_ison()) then
	  if(.not.cnvrging) then
	    write(mpout,'(2a,i2,a,i3)') myname_,	&
		': band ',ibandcg,', premature exit, iter = ',iter-1

	  elseif(.not.cnvrged) then
	    write(mpout,'(2a,i2,a,i3)') myname_,	&
		': band ',ibandcg,', max-pass exceeded, iter = ',iter

	  else
	    write(mpout,'(2a,i2,a,i3)') myname_,	&
		': band ',ibandcg,', convergence achieved, iter = ',iter
	  endif
	  endif

	  call cgnorm(cgname(ibandcg),criter(ibandcg,1),	&
		mxpass,iter,nvecs,sizerr,nvi)
	endif

			if(mall_ison()) then
			  call mall_mco(xx,myname)
			  call mall_mco(rr,myname)
			  call mall_mco(zz,myname)
			  call mall_mco(pp,myname)
			  call mall_mco(Cpp,myname)
			  call mall_mco(sizerr,myname)
			  call mall_mco(rrzznew,myname)
			  call mall_mco(rrzzold,myname)
			endif
	deallocate(xx,rr,zz,pp,Cpp,sizerr,rrzznew,rrzzold)
!***********************************************************************

end subroutine recurSolve_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: mt_Solve_ - a multitask solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine mt_Solve_(nSegi,			&
	krNav,ktNav,lcNav,lnNav,		&
	nVeci,sigU,sigC,kx,ks,qr,kl,		&
	ktab,jtab,sigFi,			&
		qrPhi,qmPhi,qlPhi,klPhi,	&
		qrPsi,qmPsi,qlPsi,klPsi,	&
	nvecs,y,x,ierrs)

      use m_die ,only : perr,die
      use m_Navigator,only : Navigator
      use m_Navigator,only : get
      use m_InnovCovMatx,only : kind_INCov
      use m_diagcorM,only : diagcorM_solve => cg1_solve
!!$      use m_diagcorM,only : diagcorM_solve => cg1_solve_1vec
      implicit none

      integer,intent(in) :: nSegi
      integer,dimension(nSegi),intent(in) :: krNav
      integer,dimension(nSegi),intent(in) :: ktNav
      integer,dimension(nSegi),intent(in) :: lcNav
      integer,dimension(nSegi),intent(in) :: lnNav

      integer ,intent(in) :: nVeci
      real   ,dimension(nVeci),intent(in) :: sigU
      real   ,dimension(nVeci),intent(in) :: sigC

      integer,dimension(nVeci),intent(in) :: kx
      integer,dimension(nVeci),intent(in) :: ks
      real ,dimension(3,nVeci),intent(in) :: qr
      integer,dimension(nVeci),intent(in) :: kl

      integer,dimension(nVeci),intent(in) :: ktab
      integer,dimension(nVeci),intent(in) :: jtab
      real   ,dimension(nVeci),intent(in) :: sigFi

      real ,dimension(3,nVeci),intent(in) :: qrPhi
      real ,dimension(3,nVeci),intent(in) :: qmPhi
      real ,dimension(3,nVeci),intent(in) :: qlPhi
      integer,dimension(nVeci),intent(in) :: klPhi

      real ,dimension(3,nVeci),intent(in) :: qrPsi
      real ,dimension(3,nVeci),intent(in) :: qmPsi
      real ,dimension(3,nVeci),intent(in) :: qlPsi
      integer,dimension(nVeci),intent(in) :: klPsi

	! Not (!) defined as assumed-shape

      integer,intent(in) :: nvecs
      real   ,dimension(nvecs,nVeci),intent(in ) :: y
      real   ,dimension(nvecs,nVeci),intent(out) :: x

      integer,dimension(nSegi),intent(out) :: ierrs

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::mt_Solve_'

  include "kind_mats.h"
  include "bands.h"

  integer :: iSeg
  integer :: kr,kt,lc,le,ln,ier

  do iSeg=1,nSegi
    kr=krNav(iSeg)
    kt=ktNav(iSeg)
    lc=lcNav(iSeg)
    ln=lnNav(iSeg)
    le=lc+ln-1
    ierrs(iSeg)=0

    call diagcorM_solve(cgverb(kind_Umat),kind_INCov,iSeg,kr,kt,	&
	ln, sigU(lc:le),sigC(lc:le),			&
	kx(lc:le),ks(lc:le),qr(:,lc:le),kl(lc:le),	&
	ktab(lc:le),jtab(lc:le), sigFi(lc:le),		&
	qrPhi(:,lc:le),qmPhi(:,lc:le),			&
	qlPhi(:,lc:le),klPhi(  lc:le),			&
	qrPsi(:,lc:le),qmPsi(:,lc:le),			&
	qlPsi(:,lc:le),klPsi(  lc:le),			&
     nvecs,y(:,lc:le),x(:,lc:le),ier)

    ierrs(iSeg)=ier
  end do

end subroutine mt_Solve_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: localSolve_ - a local solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine localSolve_(obNav,krNav,ktNav,		&
	nVeci, sigU,sigC,corObs,			&
	ktab,jtab, sigPhi,corPhi,corPsi, y,x,ierr)

      use m_die ,only : perr,die,assert_
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_InnovCovMatx,only : kind_INCov

      use m_sigmaPhi,only : sigmaPhi
      use m_sigmaPhi,only : ptr_sigma

      use m_CorAttrF,only : CorAttrF
      use m_CorAttrF,only : ptr_qr
      use m_CorAttrF,only : ptr_qm
      use m_CorAttrF,only : ptr_ql
      use m_CorAttrF,only : ptr_kl

      use m_CorAttrX,only : CorAttrX
      use m_CorAttrX,only : ptr_kx
      use m_CorAttrX,only : ptr_ks
      use m_CorAttrX,only : ptr_qr
      use m_CorAttrX,only : ptr_kl

      implicit none

      type(Navigator),intent(in) :: obNav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer ,intent(in) :: nVeci

      real   ,dimension(:),intent(in) :: sigU
      real   ,dimension(:),intent(in) :: sigC
      type(CorAttrX),intent(in) :: corObs

      integer,dimension(:),intent(in) :: ktab
      integer,dimension(:),intent(in) :: jtab

      type(sigmaPhi),intent(in) :: sigPhi
      type(CorAttrF),intent(in) :: corPhi
      type(CorAttrF),intent(in) :: corPsi

	! Not (!) defined as assumed-shape

      real   ,dimension(:,:),intent(in)  :: y	! (nveci,nvecs)
      real   ,dimension(:,:),intent(out) :: x	! (nveci,nvecs)

      integer,intent(out) :: ierr

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::localSolve_'

  integer :: lc,ln
  integer :: ier
  integer :: nSegi,iSeg
  character(len=4) :: cSeg
  integer,allocatable,dimension(:) :: lcNav,lnNav,iErrs
  real   ,allocatable,dimension(:,:) :: tx
  integer :: nvecs

	! This constraint may be removed only if the low level
	! preconditioner is redefined as an iterative solver.

	! allocate workspaces

  nvecs= size(x,2)
	ASSERT( nvecs <= size(y,2) )

  nSegi=lsize(obNav)

	allocate(lcNav(nSegi),lnNav(nSegi),iErrs(nSegi),	&
		tx(nvecs,nVeci),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

		if(mall_ison()) then
		  call mall_mci(lcNav,myname)
		  call mall_mci(lnNav,myname)
		  call mall_mci(iErrs,myname)
		  call mall_mci(tx,myname)
		endif

  do iSeg=1,nSegi
    call get(obNav,iSeg,lc=lc,ln=ln)
    lcNav(iSeg)=lc
    lnNav(iSeg)=ln
    iErrs(iSeg)=0
  end do

		! This is a local solver.  Note that the without
		! finally remove the direct solver as the low level
		! preconditioner, corPhi can not be really implemented.
		! (What I meant is that forming block matrices is a
		! headache.

    call mt_Solve_(nSegi,krNav,ktNav,lcNav,lnNav,	&
	nVeci,sigU,sigC,				&
		ptr_kx(corObs),ptr_ks(corObs),		&
		ptr_qr(corObs),ptr_kl(corObs),		&
	ktab,jtab,ptr_sigma(sigPhi),			&
	ptr_qr(corPhi),ptr_qm(corPhi),			&
	ptr_ql(corPhi),ptr_kl(corPhi),			&
	ptr_qr(corPsi),ptr_qm(corPsi),			&
	ptr_ql(corPsi),ptr_kl(corPsi),			&
	nvecs,transpose(y),tx,iErrs)

	! Verify all blocks

  ierr=0
  do iSeg=1,nSegi
    ier=iErrs(iSeg)
    if(ier/=0) then
      write(cSeg,'(i4.4)') iSeg
      call perr(myname_,'mt_Solve('//cSeg//')',ier)

      if(ierr==0) ierr=iSeg
    endif
  end do

  x(1:nVeci,1:nvecs)=transpose(tx(:,:))

		if(mall_ison()) then
		  call mall_mco(lcNav,myname)
		  call mall_mco(lnNav,myname)
		  call mall_mco(iErrs,myname)
		  call mall_mco(tx,myname)
		endif
	deallocate(lcNav,lnNav,iErrs,tx,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine localSolve_

end module m_CGSolver
