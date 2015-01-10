!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_subBlocks - Sub-blocking a PSAS vector strcuture
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_subBlocks
      use m_AttrVect, only : AttrVect	! The base data structure

      implicit none
      private	! except

      public :: subBlocks	! the data structure
      public :: init  		! construct a (subBlocks)
      public :: blockAttr	! index all the block attributes
      public :: clean		! clean a (subBlocks)
      public :: lsize		! determine the size
      public :: listing		! A debugging utility
      public :: nIAttr		! Number of attributes in the object

      type subBlocks
	type(AttrVect) :: av
	integer :: ikr,ikt,ikx,iks,ilc,iln	! Key-indices
      end type subBlocks

      interface init
	module procedure init_asis_,	&! as given kr-kt blocks
			 init_any_,	&! with a nominal block size
			 init_snd_,	&! do not break a sounding
			 init_nav_,	&! convert from a navigator
			 init_kt_	 ! block the same kt
      end interface

      interface blockAttr
	module procedure blockAttr_	! return specified attributes
      end interface

      interface clean; module procedure clean_; end interface
      interface lsize; module procedure lsize_; end interface
      interface listing;  module procedure list_;  end interface
      interface nIAttr;module procedure nIAttr_;end interface

! !REVISION HISTORY:
! 	30Oct98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	16Nov98 - Jing Guo <guo@thunder> - added init_snd_()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_subBlock'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_kt_ - block a sorted integer vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_kt_(subB,kt)
      use m_AttrVect, only : init
      use m_AttrVect, only : indexIA
      use m_AttrVect, only : ptr_iAttr
      implicit none
      type(subBlocks),intent(out) :: subB
      integer,intent(in) :: kt(:)

! !REVISION HISTORY:
! 	04Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_kt_'
  integer :: ln
  integer :: nblk
  integer :: iSav,i
  integer :: ikt,ilc,iln
  integer,pointer,dimension(:,:) :: iAttr

  ln=size(kt)
  nblk=0
  if(ln > 0) then
    nblk=1
    iSav=kt(1)

    do i=2,ln
      if(iSav /= kt(i)) then
	nblk=nblk+1
	iSav=kt(i)
      endif
    end do
  endif

  call init(subB%av,'kt:lc:ln','',nblk)

  if(nblk == 0) return

  ikt=indexIA(subB%av,'kt')
  ilc=indexIA(subB%av,'lc')
  iln=indexIA(subB%av,'ln')

  subB%ikr=-1
  subB%ikt=ikt
  subB%ikx=-1
  subB%iks=-1
  subB%ilc=ilc
  subB%iln=iln

  nblk=1
  iSav=kt(1)

  iAttr => ptr_iAttr(subB%av)

  iAttr(ikt,nblk)=iSav
  iAttr(ilc,nblk)=1

  do i=2,ln
    if(iSav /= kt(i)) then
      iAttr(iln,nblk)=i-iAttr(ilc,nblk)

      nblk=nblk+1
      iSav=kt(i)

      iAttr(ikt,nblk)=iSav
      iAttr(ilc,nblk)=i
    endif
  end do
  iAttr(iln,nblk)=ln-iAttr(ilc,nblk)+1

  nullify(iAttr)
end subroutine init_kt_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_asis_ - define a (subBlocks) as given kr-kt blocks
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_asis_(subB, nkr,kr_loc,kr_len,nkt,kt_loc,kt_len)
      use m_AttrVect, only : init
      use m_AttrVect, only : indexIA
      use m_AttrVect, only : ptr_iAttr
      implicit none

      type(subBlocks), intent(out) :: subB ! sub-block attributes

      integer,intent(in) :: nkr		! number of regions
      integer,intent(in) :: kr_loc(:)	! region locations
      integer,intent(in) :: kr_len(:)	! region lenths
      integer,intent(in) :: nkt		! number of (variable) types
      integer,intent(in) :: kt_loc(:,:)	! type loc-offsets (to a rgion)
      integer,intent(in) :: kt_len(:,:)	! type lenths

! !REVISION HISTORY:
! 	30Oct98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_asis_'

  integer :: nSub,iSub
  integer :: ikr,ikt,ilc,iln
  integer ::  kr, kt, lc, ln
  integer,pointer,dimension(:,:) :: iAttr

	! Count the number of non-zero size segments

  nSub=0
  do kr=1,nkr
    do kt=1,nkt
      if(kt_len(kt,kr) > 0) nSub=nSub+1
    end do
  end do

	! Define (AttrVect) variable

  call init(subB%av,'kr:kt:lc:ln','',nSub)

	! Index the attributes

  ikr=indexIA(subB%av,'kr')
  ikt=indexIA(subB%av,'kt')
  ilc=indexIA(subB%av,'lc')
  iln=indexIA(subB%av,'ln')

  subB%ikr=ikr
  subB%ikt=ikt
  subB%ikx=-1
  subB%iks=-1
  subB%ilc=ilc
  subB%iln=iln

  iAttr => ptr_iAttr(subB%av)

	! Set the values of block-attributes

  iSub=0
  do kr=1,nkr
    do kt=1,nkt
      if(kt_len(kt,kr) > 0) then
	iSub=iSub+1

	lc = kt_loc(kt,kr) + kr_loc(kr)
	ln = kt_len(kt,kr)

	iAttr(ikr,iSub) = kr
	iAttr(ikt,iSub) = kt
	iAttr(ilc,iSub) = lc
	iAttr(iln,iSub) = ln
      endif
    end do
  end do

  nullify(iAttr)

end subroutine init_asis_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_nav_ - define a (subBlocks) as given kr-kt blocks
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_nav_(subB, nav,krNav,ktNav)
      use m_AttrVect, only : init
      use m_AttrVect, only : indexIA
      use m_AttrVect, only : ptr_iAttr
      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      implicit none

      type(subBlocks), intent(out) :: subB ! sub-block attributes
      type(Navigator), intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

! !REVISION HISTORY:
!	31May00	- Jing Guo
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_nav_'

  integer :: inav
  integer :: ikr,ikt,ilc,iln
  integer ::  kr, kt, lc, ln
  integer,pointer,dimension(:,:) :: iAttr

	! Define (AttrVect) variable

  call init(subB%av,'kr:kt:lc:ln','',lsize(nav))

	! Index the attributes

  ikr=indexIA(subB%av,'kr')
  ikt=indexIA(subB%av,'kt')
  ilc=indexIA(subB%av,'lc')
  iln=indexIA(subB%av,'ln')

  subB%ikr=ikr
  subB%ikt=ikt
  subB%ikx=-1
  subB%iks=-1
  subB%ilc=ilc
  subB%iln=iln

  iAttr => ptr_iAttr(subB%av)

	! Set the values of block-attributes

  do inav=1,lsize(nav)
    call get(nav,inav,lc=lc,ln=ln)
    kr=krNav(inav)
    kt=ktNav(inav)

    iAttr(ikr,inav) = kr
    iAttr(ikt,inav) = kt
    iAttr(ilc,inav) = lc
    iAttr(iln,inav) = ln
  end do

  nullify(iAttr)

end subroutine init_nav_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_any_ - initialize sub-blocks with a given size
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_any_(subB, nkr,kr_loc,kr_len,nkt,kt_loc,kt_len, &
	BlockSize )
      use m_AttrVect, only : init
      use m_AttrVect, only : indexIA
      use m_AttrVect, only : ptr_iAttr
      use m_stdio,    only : stderr
      use m_die,      only : die
      implicit none

      type(subBlocks),intent(out) :: subB ! sub-block attributes

      integer,intent(in) :: nkr	! number of regions
      integer,intent(in) :: kr_loc(:)	! region locations
      integer,intent(in) :: kr_len(:)	! region lenths
      integer,intent(in) :: nkt	! number of (variable) types
      integer,intent(in) :: kt_loc(:,:)	! type locations (in the rgion)
      integer,intent(in) :: kt_len(:,:)	! type lenths

      integer,intent(in) :: BlockSize	! a nominal block size

! !REVISION HISTORY:
! 	30Oct98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_any_'

  integer :: lnB	! The nominal size of sub-blocks used in code
  integer :: nBlocks,nSubs,nB,iB
  integer :: ikr,ikt,ilc,iln
  integer ::  kr, kt, lc_i, ln, ln_i
  integer :: i
  integer,pointer,dimension(:,:) :: iAttr

	! Count the number of non-zero size sub-blocks

  lnB=BlockSize
  if(lnB<=0) then
    write(stderr,'(2a,i4)') myname_,': Invalid, BlockSize = ',lnB
    call die(myname_)
  endif

  nBlocks=0
  do kr=1,nkr
    do kt=1,nkt
      ln=kt_len(kt,kr)
      if(ln > 0) then
        nSubs=1
        if(ln > lnB) nSubs=(ln+lnB-1)/lnB
        nBlocks=nBlocks+nSubs
      endif
    enddo
  enddo

	! Define (AttrVect) variable

  call init(subB%av,'kr:kt:lc:ln','',nBlocks)

	! Index the attributes

  ikr=indexIA(subB%av,'kr')
  ikt=indexIA(subB%av,'kt')
  ilc=indexIA(subB%av,'lc')
  iln=indexIA(subB%av,'ln')

  subB%ikr=ikr
  subB%ikt=ikt
  subB%ikx=-1
  subB%iks=-1
  subB%ilc=ilc
  subB%iln=iln

  iAttr => ptr_iAttr(subB%av)

	! Set the values of block-attributes

  iB=0
  do kr=1,nkr
    do kt=1,nkt
      ln = kt_len(kt,kr)

      if(ln > 0) then

		! Location of the current block

	lc_i = kt_loc(kt,kr)+kr_loc(kr)

        nSubs=1
        if(ln > lnB) nSubs=(ln+lnB-1)/lnB

	nB=nSubs

		! Sizes of each sub-blocks are computed in a sequential
		! algorithm that adjusts the size every step to make
		! the sub-blocks more evenly divided.  e.g. instead of
		! using (3,3,3,1), this algorithm will make them as
		! (3,3,2,2).

	do i=1,nSubs
	  iB=iB+1
	  ln_i=(ln+nB-1)/nB

	  iAttr(ikr,iB) = kr
	  iAttr(ikt,iB) = kt
	  iAttr(ilc,iB) = lc_i
	  iAttr(iln,iB) = ln_i

	  lc_i=lc_i+ln_i
	  ln=ln-ln_i
	  nB=nB-1
	end do

      endif
    end do
  end do

	! Verify the check-sum, in case of anyone has changed the
	! program with inconsitant 

  if (iB/=nBlocks) then
    write(stderr,'(a,2(a,i5))') myname_,	&
	': failed invariant check, iB =', iB,	&
	' and nBlocks = ',nBlocks
    call die(myname_)
  endif

  nullify(iAttr)
end subroutine init_any_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_snd_ - initialize sub-blocks according to soundings
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_snd_(subB, nkr,kr_loc,kr_len,nkt,kt_loc,kt_len, &
	nx,kx,ks,BlockSize )
      use m_AttrVect, only : init
      use m_AttrVect, only : indexIA
      use m_AttrVect, only : ptr_iAttr
      use m_stdio,    only : stderr
      use m_die,      only : die
      implicit none

      type(subBlocks),intent(out) :: subB ! sub-block attributes

      integer,intent(in) :: nkr	! number of regions
      integer,intent(in) :: kr_loc(:)	! region locations
      integer,intent(in) :: kr_len(:)	! region lenths
      integer,intent(in) :: nkt	! number of (variable) types
      integer,intent(in) :: kt_loc(:,:)	! type locations (in the rgion)
      integer,intent(in) :: kt_len(:,:)	! type lenths

      integer,intent(in) :: nx ! number of total data
      integer,dimension(:),intent(in) :: kx	! instrument indices
      integer,dimension(:),intent(in) :: ks	! sounding indices

      integer,intent(in) :: BlockSize	! a nominal block size

! !REVISION HISTORY:
! 	16Nov98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_snd_'

  integer :: lnB	! The nominal size of sub-blocks used in code
  integer :: nBlocks,nSubs,nB,iB
  integer :: ikr,ikt,ikx,iks,ilc,iln
  integer ::  kr, kt, lc_i, ln_i
  integer :: lc,le,ln
  integer :: i
  integer,pointer,dimension(:,:) :: iAttr

	! Count the number of non-zero size sub-blocks

  lnB=BlockSize
  if(lnB<=0) then
    write(stderr,'(2a,i4)') myname_,': Invalid, BlockSize = ',lnB
    call die(myname_)
  endif

  nBlocks=0
  do kr=1,nkr
    do kt=1,nkt
      ln=kt_len(kt,kr)

      if(ln > 0) then

        lc=kt_loc(kt,kr)+kr_loc(kr)
	le=lc+ln-1

	ln_i=1

	do i=lc+1,le

	  if( ln_i >= lnB      .and.	&	! at least this long,
	      ( kx(i)/=kx(i-1) .or.	&	! also broken only where
	        ks(i)/=ks(i-1) )	) then	! kx-ks switch

	    nBlocks=nBlocks+1	! count as a segment
	    ln_i=0		! start the next segment
	  endif

	  ln_i=ln_i+1
	end do

	if(ln_i > 0) nBlocks=nBlocks+1	! count the rest as a segment

      endif
    enddo
  enddo

	! Define (AttrVect) variable

  call init(subB%av,'kr:kt:kx:ks:lc:ln','',nBlocks)

	! Index the attributes

  ikr=indexIA(subB%av,'kr')
  ikt=indexIA(subB%av,'kt')
  ikx=indexIA(subB%av,'kx')
  iks=indexIA(subB%av,'ks')
  ilc=indexIA(subB%av,'lc')
  iln=indexIA(subB%av,'ln')

  subB%ikr=ikr
  subB%ikt=ikt
  subB%ikx=ikx
  subB%iks=iks
  subB%ilc=ilc
  subB%iln=iln

  iAttr => ptr_iAttr(subB%av)

	! Set the values of block-attributes

  iB=0
  do kr=1,nkr
    do kt=1,nkt
      ln = kt_len(kt,kr)

      if(ln > 0) then

		! Location of the current block

	lc = kt_loc(kt,kr)+kr_loc(kr)
	le = lc+ln-1

	ln_i=1
	lc_i=lc
	do i=lc+1,le

	  if( ln_i >= lnB      .and.	&	! at least this long,
	      ( kx(i)/=kx(i-1) .or.	&	! also broken only where
	        ks(i)/=ks(i-1) )	) then	! kx-ks switch

	    iB=iB+1	! count as a segment

	    iAttr(ikr,iB) = kr
	    iAttr(ikt,iB) = kt
	    iAttr(ikx,iB) = kx(lc_i)
	    iAttr(iks,iB) = ks(lc_i)
	    iAttr(ilc,iB) = lc_i
	    iAttr(iln,iB) = ln_i

	    lc_i=lc_i+ln_i
	    ln_i=0		! start the next segment
	  endif

	  ln_i=ln_i+1
	end do

	if(ln_i > 0) then
	  iB=iB+1

	  iAttr(ikr,iB) = kr
	  iAttr(ikt,iB) = kt
	  iAttr(ikx,iB) = kx(lc_i)
	  iAttr(iks,iB) = ks(lc_i)
	  iAttr(ilc,iB) = lc_i
	  iAttr(iln,iB) = ln_i

	endif

      endif
    end do
  end do

	! Verify the check-sum, in case of anyone has changed the
	! program with inconsitant 

  if (iB/=nBlocks) then
    write(stderr,'(a,2(a,i5))') myname_,	&
	': failed invariant check, iB =', iB,	&
	' and nBlocks = ',nBlocks
    call die(myname_)
  endif

  nullify(iAttr)
end subroutine init_snd_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: blockAttr_ - return specified block attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine blockAttr_(subB,iBlk,kr,kt,kx,ks,lc,ln)
      use m_stdio, only : stderr
      use m_die,      only : die
      use m_AttrVect,only : ptr_iAttr
      use m_AttrVect,only : lsize
      implicit none
      type(subBlocks), intent(in)  :: subB
      integer,         intent(in)  :: iBlk
      integer,optional,intent(out) :: kr
      integer,optional,intent(out) :: kt
      integer,optional,intent(out) :: kx
      integer,optional,intent(out) :: ks
      integer,optional,intent(out) :: lc
      integer,optional,intent(out) :: ln

! !REVISION HISTORY:
! 	09Nov98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::blockAttr_'
  integer,pointer,dimension(:,:) :: iAttr

  if(iBlk <= 0 .or. iBlk > lsize(subB%av)) then
    write(stderr,'(2a,i6)') myname_,	&
	': the second dimension overflow, iBlk =',iBlk
    call die(myname_)
  endif

  if( (present(kx).and.subB%ikx<=0)	.or.	&
      (present(ks).and.subB%iks<=0)	)	then

    if( present(kx).and.subB%ikx<=0 )	&
	write(stderr,'(2a)') myname_, ': undefined %ikx.'
    if( present(ks).and.subB%iks<=0 )	&
	write(stderr,'(2a)') myname_, ': undefined %iks.'

    call die(myname_)
  endif

  iAttr => ptr_iAttr(subB%av)
  if(present(kr).and.subB%ikr>0) kr=iAttr(subB%ikr,iBlk)
  if(present(kt).and.subB%ikt>0) kt=iAttr(subB%ikt,iBlk)
  if(present(kx).and.subB%ikx>0) kx=iAttr(subB%ikx,iBlk)
  if(present(ks).and.subB%iks>0) ks=iAttr(subB%iks,iBlk)
  if(present(lc).and.subB%ilc>0) lc=iAttr(subB%ilc,iBlk)
  if(present(ln).and.subB%iln>0) ln=iAttr(subB%iln,iBlk)
  nullify(iAttr)

end subroutine blockAttr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a (subBlocks) attribute vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(subB)
      use m_AttrVect, only : clean
      implicit none
      type(subBlocks), intent(inout) :: subB

! !REVISION HISTORY:
! 	09Nov98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  call clean(subB%av)
  subB%ikr=-1
  subB%ikt=-1
  subB%ikx=-1
  subB%iks=-1
  subB%ilc=-1
  subB%iln=-1
end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - the size of the given (subBlocks)
!
! !DESCRIPTION:
!
! !INTERFACE:

    function lsize_(subB)
      use m_AttrVect, only : lsize
      implicit none
      type(subBlocks),intent(in) :: subB
      integer :: lsize_

! !REVISION HISTORY:
! 	09Nov98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_=lsize(subB%av)
end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nIAttr_ - the number of attributes of in the object
!
! !DESCRIPTION:
!
! !INTERFACE:

    function nIAttr_(subB)
      use m_AttrVect, only : nIAttr
      implicit none
      type(subBlocks),intent(in) :: subB
      integer :: nIAttr_

! !REVISION HISTORY:
! 	09Nov98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nIAttr_'

  nIAttr_=nIAttr(subB%av)
end function nIAttr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: list_ - List the contents of a (subBlocks)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine list_(subB)
      use m_stdio, only : stdout
      use m_AttrVect,only : ptr_iAttr
      implicit none
      type(subBlocks),intent(in) :: subB

! !REVISION HISTORY:
! 	30Nov98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::list_'

  integer :: nSubs
  integer :: ikr,ikt,ilc,iln
  integer ::  kr, kt, lc, ln
  integer :: i,n
  integer,pointer,dimension(:,:) :: iAttr

  nSubs=lsize_(subB)
  if(nSubs <= 0) return

  ikr=subB%ikr
  ikt=subB%ikt
  ilc=subB%ilc
  iln=subB%iln

  iAttr => ptr_iAttr(subB%av)

  n=1
  kr=iAttr(ikr,1)
  kt=iAttr(ikt,1)
  lc=iAttr(ilc,1)
  ln=iAttr(iln,1)

  do i=2,nSubs
    if( iAttr(ikr,i) /= kr      .or.    &
        iAttr(ikt,i) /= kt      )       then

      write(stdout,'(3i4,2i8)') n,kr,kt,lc,ln

      n=0
      kr=iAttr(ikr,i)
      kt=iAttr(ikt,i)
      lc=iAttr(ilc,i)
      ln=0
    endif
    n=n+1
    ln=ln+iAttr(iln,i)
  end do
  write(stdout,'(3i4,2i8)') n,kr,kt,lc,ln

  nullify(iAttr)

end subroutine list_
end module m_subBlocks
