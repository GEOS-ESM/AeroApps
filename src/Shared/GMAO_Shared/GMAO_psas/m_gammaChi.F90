!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_gammaChi - wind balancer: Velosity potential
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_gammaChi
      implicit none
      private	! except

      public :: gammaChi_matvecpy
      public :: gammaChi_adjvec

      interface gammaChi_matvecpy; module procedure	&
	matvecpy_; end interface
      interface gammaChi_adjvec; module procedure	&
	adjvec_; end interface

! !REVISION HISTORY:
!
! 	24May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_gammaChi'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
!
! !IROUTINE: matvecpy_ - convert (chi-m,chi-l) vectors to (U,V) vectors
!
! !DESCRIPTION:
!
!	matvecpy_() converts a (chi-m,chi-l) vector to a (U,V) vector,
!	by applying a wind-balance (through velosity potential, chi)
!	operator to the input argument Hvec ([U,V] = G[chi-m,chi-l]).
!	Both Xvec and Hvec share the same navigator, except that the
!	space for wind components (U,V) in Xvec corresponds to the
!	space for gradiant components (chi-m, chi-l) in Hvec.  The
!	location attributes (ktab, jtab) are defined through the
!	attributes of Xvec.
!
! !INTERFACE:

    subroutine matvecpy_(nav,krNav,ktNav,	&
	nsize,ktab,jtab,nvecs,Hvec,Xvec		)

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_die,  only : die,perr

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
!
! 	24May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- rewrote with new arguments and as a module.
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

  integer :: kt,i,ivec
  integer :: lc  ,le
  integer :: lc_m,le_m
  integer :: lc_l,le_l
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
	! loop over all (X) partitions to convert (Hm,Hl) to (u,v)

  do inav=1,lsize(nav)

    kt=ktNav(inav)

    select case(kt)
    case (ktUU,ktUs)

      call get(nav,inav  ,lc=lc  ,le=le  )
      call get(nav,inav+1,lc=lc_l,le=le_l)
      do ivec=1,nvecs
	Xvec(lc:le,ivec) = -Hvec(lc_l:le_l,ivec) + Xvec(lc:le,ivec)
      end do

    case (ktVV,ktVs)

      call get(nav,inav  ,lc=lc  ,le=le  )
      call get(nav,inav-1,lc=lc_m,le=le_m)
      do ivec=1,nvecs
	Xvec(lc:le,ivec) = -Hvec(lc_m:le_m,ivec) + Xvec(lc:le,ivec)
      end do

    end select

  end do

end subroutine matvecpy_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
!
! !IROUTINE: adjvec_ - convert (U,V) vectors to (chi-m,chi-l) vectors
!
! !DESCRIPTION:
!
!	adjvec_() converts a (U,V) vector to a (chi-m,chi-l) vector,
!	by applying the _adjoint_ of a wind-balance (through velosity
!	potential, chi) operator to the input argument Xvec
!	([U,V] = G'[chi-m,chi-l]).  Both Xvec and Hvec share the same
!	navigator, except that the space for gradient components
!	(chi-m,chi-l) in Hvec corresponds to the space for wind
!	components (U,V) in Xvec.  The location attributes (ktab, jtab)
!	are defined through the attributes of Xvec.
!
! !INTERFACE:

    subroutine adjvec_(nav,krNav,ktNav,	&
	nsize,ktab,jtab,nvecs,Xvec,Hvec		)

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_die,  only : die,perr

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
!
! 	24May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- redefined the intent of Hvec from inout to out;
!		- rewrote with new arguments and as a module.
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

  integer :: kt,i,ivec
  integer :: lc  ,le
  integer :: lc_u,le_u
  integer :: lc_v,le_v
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
	! loop over all (H) partitions to convert (Hm,Hl) to (u,v)

  do inav=1,lsize(nav)

    kt=ktNav(inav)

    select case(kt)
    case (ktUU,ktUs)

      call get(nav,inav  ,lc=lc  ,le=le  )
      call get(nav,inav+1,lc=lc_v,le=le_v)

      do ivec=1,nvecs
	Hvec(lc:le,ivec) = -Xvec(lc_v:le_v,ivec)
      end do

    case (ktVV,ktVs)

      call get(nav,inav  ,lc=lc  ,le=le  )
      call get(nav,inav-1,lc=lc_u,le=le_u)

      do ivec=1,nvecs
	Hvec(lc:le,ivec) = -Xvec(lc_u:le_u,ivec)
      end do

    case default
      call get(nav,inav,lc=lc,le=le)
      do ivec=1,nvecs
        Hvec(lc:le,ivec)=0.
      end do
    end select

  end do

end subroutine adjvec_

end module m_gammaChi
