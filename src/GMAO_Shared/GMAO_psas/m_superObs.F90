!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_superObs - Super-observation processor
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_superObs
      implicit none
      private	! except

      public :: superObs_prox

      interface superObs_prox; module procedure	prox_; end interface

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_superObs'

  character(len=*),parameter :: PROXRC="superobs.rc"

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize this superObs processor
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_()
      use m_inpak90,only : i90_loadf,i90_release
      use m_mpout,  only : mpout_log
      use m_die,    only : die
      use m_boxes,  only : boxes_init
      implicit none

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

  integer :: ier

  call i90_loadf(PROXRC,ier)
	if(ier.ne.0) call die(myname_,'i90_loadf("'//PROXRC//'")',ier)

	call mpout_log(myname_,		&
	  ': using "'//PROXRC//'" for the runtime configuration')

  call kxname0()	! load the resource
  call proxel0()
  call boxes_init()

  call i90_release(ier)
	if(ier/=0) call die(myname_,'i90_release("'//PROXRC//'")',ier)

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean resources
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      use m_boxes,only : boxes_clean
      implicit none

! !REVISION HISTORY:
! 	20Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  call boxes_clean()

end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: qualify_ - qualify observations
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine qualify_(nobs,qlst,rlat,rlon,rlev,tims,	&
	kxs,kts,ktLst)
      use m_die, only : die
      use m_boxes,only : boxes_restrict
      implicit none

      integer,intent(in) :: nobs
      logical,dimension(:),intent(out):: qlst
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon
      real   ,dimension(:),intent(in) :: rlev
      real   ,dimension(:),intent(in) :: tims
      integer,dimension(:),intent(in) :: kxs
      integer,dimension(:),intent(in) :: kts

      integer,dimension(:),optional,intent(in) :: ktLst

! !REVISION HISTORY:
! 	07Feb00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::qualify_'

  integer,allocatable,dimension(:) :: nkt
  real   ,allocatable,dimension(:) :: sigF,sigU,sigO
  integer :: ier
  integer :: i,k,kt

!________________________________________

	! Pre-partitioning processing, making sure the data
	! to be partitioned are all qualified 
	!________________________________________

	allocate(  nkt(nobs),sigF(nobs),sigU(nobs),sigO(nobs),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

	! Work on a copy of the input data.  This copy will become
	! the output.  It is assumed the sizes of the output(=workspace)
	! arrays are at least as large as the input arrays.

	! All data are assumed to be qualified at the start.

  do i=1,nobs
    qlst(i)=.true.
  end do

	! Mark observations with not listed kt values

  if(present(ktLst)) then	! only data with listed kt values

    do i=1,nobs		! none has been verified
      nkt(i)=-1
    end do

    do k=1,size(ktLst)
      if(ktLst(k) <= 0) call die(myname_,	&
		'invalide ktLst value',ktLst(k)	)
    end do

    do k=1,size(ktLst)
      kt=ktLst(k)

	! Mark the flag with its true kt value, if this kt value is
	! listed.

      do i=1,nobs
	if(kts(i)==kt) nkt(i)=kt
      end do
    end do

	! Any observation with a non-listed kt value would not be
	! wanted.

    do i=1,nobs
      if(nkt(i)==-1) qlst(i)=.false.
    end do

  endif

  call boxes_restrict(nobs,qlst,rlat,rlon,rlev,tims,kxs,kts)
	!________________________________________

  ! call psas_sigs(nobs,rlat,rlon,kxs,kts,sigF,sigO,sigU)
  ! call mark_nsig(nobs,sigU,sigO,qlst)
	!________________________________________

	deallocate(nkt,sigF,sigU,sigO, stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine qualify_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: proxel_ 
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine proxel_(nxkr,nxkt,lc_kr,ln_kr,lc_kt,ln_kt,	&
	nobs,rlat,rlon,rlev,tims,kxs,kts,xmUs,dels	)

      use m_mpout,only : mpout
      use m_die,  only : die
      implicit none

      integer,intent(in) :: nxkr
      integer,intent(in) :: nxkt
      integer,dimension(:)  ,intent(inout) :: lc_kr
      integer,dimension(:)  ,intent(inout) :: ln_kr
      integer,dimension(:,:),intent(inout) :: lc_kt
      integer,dimension(:,:),intent(inout) :: ln_kt

      integer,intent(inout) :: nobs
      real   ,dimension(:),intent(inout) :: rlat
      real   ,dimension(:),intent(inout) :: rlon
      real   ,dimension(:),intent(inout) :: rlev
      real   ,dimension(:),intent(inout) :: tims
      integer,dimension(:),intent(inout) :: kxs
      integer,dimension(:),intent(inout) :: kts
      real   ,dimension(:),intent(inout) :: xmUs
      real   ,dimension(:),intent(inout) :: dels

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::proxel_'

  logical,parameter :: verbose=.false.

  integer :: ier
  integer :: n,nprox
  integer :: kr,kt,lc
  logical,allocatable,dimension(:) :: kls
  real   ,allocatable,dimension(:) :: sigO,sigF

!________________________________________

	allocate(kls(nobs),sigO(nobs),sigF(nobs), stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

	  sigO(1:nobs)=0.
	  sigF(1:nobs)=0.

!	..Remove "duplicates" in the observations and adjust lc_kr,
!	ln_kr, lc_kt, and ln_kt accordingly.

  call dupelim(verbose, mpout,				&
	nobs, kxs, kts, kls,rlat,rlon,rlev,dels,	&
	     xmUs,sigO,sigF,tims,			&
	nxkr,lc_kr,ln_kr,nxkt,ln_kt			)

!	.."Superob" observations that are within a given range.  Quit
!	searching loop if nothing to "superob", or have looped 5 times.
!	lc_kr, ln_kr, and ln_kt arrays are adjusted accordingly.

  nprox=0
  n=1
  do while(n.eq.1 .or. nprox.ne.0.and.n.le.5)
    call proxel ( verbose, mpout,			&
	nobs, kxs, kts, kls,rlat,rlon,rlev,dels,	&
	     sigO,xmUs,sigF,tims,			&
	nxkr,lc_kr,ln_kr,nxkt,ln_kt, nprox		)
    n=n+1
  end do

  do kr=1,nxkr
    lc=0
    do kt=1,nxkt
      lc_kt(kt,kr)=lc
      lc=lc+ln_kt(kt,kr)
    end do
  end do

	deallocate(kls,sigO,sigF,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine proxel_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: prox_ - creating super-Obs.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine prox_(	&
	nout,olat,olon,olev,otim,okx,oks,okt,odel,oxmU,	&
	ninp,rlat,rlon,rlev,tims,kxs,kts,dels,ktLst	)

      use m_ktList, only : KTMAX
      use m_Sndx,only : setSndx
      use m_Regioner,only : Regioner_init
      use m_Regioner,only : Regioner_clean
      use m_Regioner,only : Regioner_part
      use m_die,  only : die
      use m_ioutil,only : luflush
      implicit none

      integer,intent(out) :: nout
      real   ,dimension(:),intent(out) :: olat
      real   ,dimension(:),intent(out) :: olon
      real   ,dimension(:),intent(out) :: olev
      real   ,dimension(:),intent(out) :: otim
      integer,dimension(:),intent(out) :: okx
      integer,dimension(:),intent(out) :: oks
      integer,dimension(:),intent(out) :: okt
      real   ,dimension(:),intent(out) :: odel
      real   ,dimension(:),intent(out) :: oxmU

      integer,intent(in) :: ninp
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon
      real   ,dimension(:),intent(in) :: rlev
      real   ,dimension(:),intent(in) :: tims
      integer,dimension(:),intent(in) :: kxs
      integer,dimension(:),intent(in) :: kts
      real   ,dimension(:),intent(in) :: dels

      integer,dimension(:),optional,intent(in) :: ktLst

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::prox_'

  integer :: i,k,ier
  integer :: nxkr,nxkt
  integer,allocatable,dimension(:)   :: ln_kr,lc_kr
  integer,allocatable,dimension(:,:) :: ln_kt,lc_kt

  logical,allocatable,dimension(:)   :: qlst
!________________________________________

  call init_()
	!________________________________________

	allocate(qlst(ninp),stat=ier)
		if(ier/=0) call die(myname_,'allocate(qlst)',ier)
	!________________________________________

  call qualify_(ninp,qlst,rlat,rlon,rlev,tims,kxs,kts)

	! Collect all elements with a "qualified" data.  See qualify_()
	! for the meaning of "qualified" data.

    k=0
    do i=1,ninp
      if( qlst(i) ) then
	k=k+1

	olat(k)=rlat(i)
	olon(k)=rlon(i)
	olev(k)=rlev(i)
	otim(k)=tims(i)
	odel(k)=dels(i)
	 okx(k)= kxs(i)
	 okt(k)= kts(i)

      endif
    end do
    nout=k
	!________________________________________

	deallocate(qlst,stat=ier)
		if(ier/=0) call die(myname_,'deallocate(qlst)',ier)
	!________________________________________

  if(nout==0) then
    call clean_()
    return
  endif
!________________________________________

  nxkt=KTMAX
  call Regioner_init(nxkr)

	allocate( lc_kr(     nxkr),ln_kr(     nxkr),	&
		  lc_kt(nxkt,nxkr),ln_kt(nxkt,nxkr), stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call Regioner_part(nxkr,nxkt,lc_kr,ln_kr,lc_kt,ln_kt,	&
	nout,olat,olon,olev,otim,okx,okt,odel	)
	!________________________________________

	! Post-partitioning processing:
	!________________________________________

	oxmU(1:nout)=1.	! Assumed a local constant sigU for the
			! purpose of superObs only.

  call proxel_(nxkr,nxkt,lc_kr,ln_kr,lc_kt,ln_kt,	&
	nout,olat,olon,olev,otim,okx,okt,oxmU,odel	)

		! In this procedure, lc_kt is not used at all.
		! Neverthless, it is still defined by Regioner_part()
		! and redefined by proxel_().

	deallocate(lc_kr,ln_kr,lc_kt,ln_kt,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

  call setSndx(oks(1:nout),okx,olat,olon)

  call Regioner_clean()

  call clean_()
!________________________________________

end subroutine prox_
end module m_superObs
!.
