!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_gammaPhi - wind-mass balance operator 
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_gammaPhi
      implicit none
      private	! except

      public :: gammaPhi_matvecpy
      public :: gammaPhi_adjvec

      interface gammaPhi_matvecpy; module procedure	&
	matvecpy_; end interface
      interface gammaPhi_adjvec; module procedure	&
	adjvec_; end interface

! !REVISION HISTORY:
! 	24May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_gammaPhi'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
!
! !IROUTINE: matvecpy_ -  convert (H,Hm,Hl) vectors to (H,U,V) vectors
!
! !DESCRIPTION:
!
!	matvecpy_() converts a (H,Hm,Hl,..) vector to a (H,U,V,..)
!	vector, by applying a wind-mass-balance operator to the input
!	Hvec ([H,U,V] = A[H,Hm,Hl]).  Both Hvec and Xvec share the same
!	navigator, except that the space for wind components (H,U,V) in
!	Xvec corresponds to the space for gradient components (H,Hm,Hl)
!	in Hvec.  The location attributes (ktab,jtab) are defined
!	through the attributes of Xvec.
!
! !INTERFACE:

    subroutine matvecpy_(nav,krNav,ktNav,	&
	nsize,ktab,jtab,nvecs,Hvec,Xvec		)

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_die,  only : die,perr
      use m_mall, only : mall_ison,mall_mci,mall_mco
      use FEalpha_imat

      implicit none
      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer,intent(in) :: nsize
      integer,dimension(:),intent(in) :: ktab
      integer,dimension(:),intent(in) :: jtab

      integer,intent(in) :: nvecs
      real   ,dimension(:,:),intent(in   ) :: Hvec
      real   ,dimension(:,:),intent(inout) :: Xvec

! !REVISION HISTORY:
! 	24May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- rewrote with new arguments and as a module
!	10Jul98 - J. Guo
!		- changed the intent of Xvec(:,:) from
!		  (out) to (inout) to fix an apperent bug.  There
!		  are no apperent effect on the results by this
!		  fix.  The final solution is yet to be studied.
! 	07Feb96 - J. Guo
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::matvecpy_'
  include "ktmax.h"

  real,allocatable :: a_u(:) ! (nsize) loosely dimensioned (a_um a_ul)
  real,allocatable :: a_v(:) ! (nsize) loosely dimensioned (a_vm a_vl)

  integer :: kt,i,ivec,ktx
  integer :: lc_m,le_m,ln_m
  integer :: lc_l,le_l,ln_l
  integer :: lc  ,le  ,ln
  integer :: ier
  integer :: inav

  logical :: invalid	! a condition flag

!-----------------------------------------------------------------------
  
	! Dimension consistency checking

  invalid=.false.
  if(nsize>size(Hvec,1)) then
    call perr(myname_,'nsize',nsize,'size(Hvec,1)',size(Hvec,1))
    invalid=.true.
  endif

  if(nsize>size(Xvec,1)) then
    call perr(myname_,'nsize',nsize,'size(Xvec,1)',size(Xvec,1))
    invalid=.true.
  endif

  if(nvecs>size(Hvec,2)) then
    call perr(myname_,'nvecs',nvecs,'size(Hvec,2)',size(Hvec,2))
    invalid=.true.
  endif

  if(nvecs>size(Xvec,2)) then
    call perr(myname_,'nvecs',nvecs,'size(Xvec,2)',size(Xvec,2))
    invalid=.true.
  endif

  if(invalid) call die(myname_)

!-----------------------------------------------------------------------
	! allocate workspace

	allocate(a_u(nsize),a_v(nsize),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(a_u,myname)
		  call mall_mci(a_v,myname)
		endif

!-----------------------------------------------------------------------
	! loop over all (X) partitions to convert (Hm,Hl) to (u,v)

  do inav=1,lsize(nav)

    kt=ktNav(inav)

    select case(kt)
    case (ktUU,ktUs)
!________________________________________

      if(inav+1>lsize(nav))    &
	call die(myname_,'inav+1',inav+1,'lsize(nav)',lsize(nav))

      ktx=ktNav(inav+1)

      select case(kt)
      case (ktUU)
        if(ktx/=ktVV) call die(myname_,'kt_m',kt,'kt_l',ktx)
      case (ktUs)
        if(ktx/=ktVs) call die(myname_,'kt_m',kt,'kt_l',ktx)
      end select

      call get(nav,inav  ,lc=lc  ,le=le  ,ln=ln  )
      call get(nav,inav+1,lc=lc_l,le=le_l,ln=ln_l)

      lc_m=lc
      le_m=le
      ln_m=ln

      invalid=ln_m/=ln_l
      if(.not.invalid) invalid =	&
	any(ktab(lc_m:le_m)/=ktab(lc_l:le_l)) .or.	&
	any(jtab(lc_m:le_m)/=jtab(lc_l:le_l))

      if(invalid) then
	call perr(myname_,'ln_m',ln_m,'ln_l',ln_l)
	call die(myname_,'unpaired u-v',kt)
      endif
!________________________________________

      a_u(lc_m:le_m) =	&
	(/ (Aum_imat(ktab(i),jtab(i)), i=lc,le) /)

      a_u(lc_l:le_l) =	&
	(/ (Aul_imat(ktab(i),jtab(i)), i=lc,le) /)

      do ivec=1,nvecs
	Xvec(lc:le, ivec) = Xvec(lc:le, ivec)    +	&
	  a_u(lc_m:le_m) * Hvec(lc_m:le_m, ivec) +	&
	  a_u(lc_l:le_l) * Hvec(lc_l:le_l, ivec)
      end do

    case (ktVV,ktVs)

      if(inav-1<1) call die(myname_,'inav-1',inav-1)

      ktx=ktNav(inav-1)

      select case(kt)
      case (ktVV)
        if(ktx/=ktUU) call die(myname_,'kt_m',ktx,'kt_l',kt)
      case (ktVs)
        if(ktx/=ktUs) call die(myname_,'kt_m',ktx,'kt_l',kt)
      end select

      call get(nav,inav  ,lc=lc  ,le=le  ,ln=ln  )
      call get(nav,inav-1,lc=lc_m,le=le_m,ln=ln_m)

      lc_l=lc
      le_l=le
      ln_l=ln

      invalid=ln_m/=ln_l
      if(.not.invalid) invalid =	&
	any(ktab(lc_m:le_m)/=ktab(lc_l:le_l)) .or.	&
	any(jtab(lc_m:le_m)/=jtab(lc_l:le_l))

      if(invalid) then
	call perr(myname_,'ln_m',ln_m,'ln_l',ln_l)
	call die(myname_,'unpaired u-v',kt)
      endif

      a_v(lc_m:le_m) =	&
	(/ (Avm_imat(ktab(i),jtab(i)), i=lc,le) /)

      a_v(lc_l:le_l) =	&
	(/ (Avl_imat(ktab(i),jtab(i)), i=lc,le) /)

      do ivec=1,nvecs
	Xvec(lc:le, ivec) = Xvec(lc:le, ivec)    +	&
	  a_v(lc_m:le_m) * Hvec(lc_m:le_m, ivec) +	&
	  a_v(lc_l:le_l) * Hvec(lc_l:le_l, ivec)
      end do

    case default
      call get(nav,inav,lc=lc,le=le)
      do ivec=1,nvecs
        Xvec(lc:le,ivec)=Xvec(lc:le,ivec)+Hvec(lc:le,ivec)
      end do
    end select

  end do

		if(mall_ison()) then
		  call mall_mco(a_u,myname)
		  call mall_mco(a_v,myname)
		endif
	deallocate(a_u,a_v,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine matvecpy_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
!
! !IROUTINE: adjvec_ - convert (H,U,V) vectors to (H,Hm,Hl) vectors
!
! !DESCRIPTION:
!
!	adjvec_() converts a (H,U,V,..) vector to a (H,Hm,Hl,..)
!	vector, by applying the adjoint of a wind-mass-balance operator
!	to the input Xvec ([H,Hm,Hl] = A'[H,U,V]).  Both Hvec and Xvec
!	share the same navigator, except that the space for wind-mass
!	components (H,U,V) in Xvec corresponds to the space for mass-
!	gradient components (H,Hm, Hl) in Hvec.  The location
!	attributes (ktab, jtab) are defined through the attributes of
!	Xvec.
!
! !INTERFACE:

    subroutine adjvec_(nav,krNav,ktNav,	&
	nsize,ktab,jtab,nvecs,Xvec,Hvec		)

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_die,  only : die,perr
      use m_mall, only : mall_ison,mall_mci,mall_mco
      use FEalpha_imat

      implicit none
      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer,intent(in) :: nsize
      integer,dimension(:),intent(in) :: ktab
      integer,dimension(:),intent(in) :: jtab

      integer,intent(in) :: nvecs
      real   ,dimension(:,:),intent(in ) :: Xvec
      real   ,dimension(:,:),intent(out) :: Hvec

! !REVISION HISTORY:
! 	24May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- rewrote with new arguments and as a module
!		- redefined the intent of Hvec from inout to out.
!	10Jul98 - J. Guo
!		- changed the intent of Hvec(:,:) from
!		  (out) to (inout) to fix an apperent bug.  There
!		  are no apperent effect on the results by this
!		  fix.  The final solution is yet to be studied.
! 	07Feb96 - J. Guo
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::adjvec_'
  include "ktmax.h"

  real,allocatable :: a_m(:) ! (nsize) loosely dimensioned (a_um a_ul)
  real,allocatable :: a_l(:) ! (nsize) loosely dimensioned (a_vm a_vl)

  integer :: kt,i,ivec,ktx
  integer :: lc_u,le_u,ln_u
  integer :: lc_v,le_v,ln_v
  integer :: lc  ,le  ,ln
  integer :: ier
  integer :: inav

  logical :: invalid	! a condition flag

!-----------------------------------------------------------------------
  
	! Dimension consistency checking

  invalid=.false.
  if(nsize>size(Hvec,1)) then
    call perr(myname_,'nsize',nsize,'size(Hvec,1)',size(Hvec,1))
    invalid=.true.
  endif

  if(nsize>size(Xvec,1)) then
    call perr(myname_,'nsize',nsize,'size(Xvec,1)',size(Xvec,1))
    invalid=.true.
  endif

  if(nvecs>size(Hvec,2)) then
    call perr(myname_,'nvecs',nvecs,'size(Hvec,2)',size(Hvec,2))
    invalid=.true.
  endif

  if(nvecs>size(Xvec,2)) then
    call perr(myname_,'nvecs',nvecs,'size(Xvec,2)',size(Xvec,2))
    invalid=.true.
  endif

  if(invalid) call die(myname_)

!-----------------------------------------------------------------------
	! allocate workspace

	allocate(a_m(nsize),a_l(nsize),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(a_m,myname)
		  call mall_mci(a_l,myname)
		endif

!-----------------------------------------------------------------------
	! loop over all (H) partitions to convert (Hm,Hl) to (u,v)

  do inav=1,lsize(nav)

    kt=ktNav(inav)

    select case(kt)
    case (ktUU,ktUs)
!________________________________________

      if(inav+1>lsize(nav))    &
	call die(myname_,'inav+1',inav+1,'lsize(nav)',lsize(nav))

      ktx=ktNav(inav+1)

      select case(kt)
      case (ktUU)
        if(ktx/=ktVV) call die(myname_,'kt_u',kt,'kt_v',ktx)
      case (ktUs)
        if(ktx/=ktVs) call die(myname_,'kt_u',kt,'kt_v',ktx)
      end select

      call get(nav,inav  ,lc=lc  ,le=le  ,ln=ln  )
      call get(nav,inav+1,lc=lc_v,le=le_v,ln=ln_v)

      lc_u=lc
      le_u=le
      ln_u=ln

      invalid=ln_u/=ln_v
      if(.not.invalid) invalid =	&
	any(ktab(lc_u:le_u)/=ktab(lc_v:le_v)) .or.	&
	any(jtab(lc_u:le_u)/=jtab(lc_v:le_v))

      if(invalid) then
	call perr(myname_,'ln_u',ln_u,'ln_v',ln_v)
	call die(myname_,'unpaired u/v',kt)
      endif
!________________________________________

      a_m(lc_u:le_u) =	&
	(/ (Aum_imat(ktab(i),jtab(i)), i=lc_u,le_u) /)

      a_m(lc_v:le_v) =	&
	(/ (Avm_imat(ktab(i),jtab(i)), i=lc_v,le_v) /)

      do ivec=1,nvecs
	Hvec(lc:le, ivec) =				&
	  a_m(lc_u:le_u) * Xvec(lc_u:le_u, ivec) +	&
	  a_m(lc_v:le_v) * Xvec(lc_v:le_v, ivec)
      end do

    case (ktVV,ktVs)
!________________________________________

      if(inav-1<1) call die(myname_,'inav-1',inav-1)

      ktx=ktNav(inav-1)

      select case(kt)
      case (ktVV)
        if(ktx/=ktUU) call die(myname_,'kt_u',ktx,'kt_v',kt)
      case (ktVs)
        if(ktx/=ktUs) call die(myname_,'kt_u',ktx,'kt_v',kt)
      end select

      call get(nav,inav  ,lc=lc  ,le=le  ,ln=ln  )
      call get(nav,inav-1,lc=lc_u,le=le_u,ln=ln_u)

      lc_v=lc
      le_v=le
      ln_v=ln

      invalid=ln_u/=ln_v
      if(.not.invalid) invalid =	&
	any(ktab(lc_u:le_u)/=ktab(lc_v:le_v)) .or.	&
	any(jtab(lc_u:le_u)/=jtab(lc_v:le_v))

      if(invalid) then
	call perr(myname_,'ln_u',ln_u,'ln_v',ln_v)
	call die(myname_,'unpaired u/v',kt)
      endif
!________________________________________

      a_l(lc_u:le_u) =	&
	(/ (Aul_imat(ktab(i),jtab(i)), i=lc_u,le_u) /)

      a_l(lc_v:le_v) =	&
	(/ (Avl_imat(ktab(i),jtab(i)), i=lc_v,le_v) /)

      do ivec=1,nvecs
	Hvec(lc:le, ivec) =				&
	  a_l(lc_u:le_u) * Xvec(lc_u:le_u, ivec) +	&
	  a_l(lc_v:le_v) * Xvec(lc_v:le_v, ivec)
      end do

    case default
      call get(nav,inav,lc=lc,le=le)
      do ivec=1,nvecs
        Hvec(lc:le,ivec)=Xvec(lc:le,ivec)
      end do
    end select

  end do

		if(mall_ison()) then
		  call mall_mco(a_m,myname)
		  call mall_mco(a_l,myname)
		endif
	deallocate(a_m,a_l,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine adjvec_

end module m_gammaPhi
