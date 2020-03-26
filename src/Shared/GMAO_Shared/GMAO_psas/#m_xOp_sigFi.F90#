!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_xOp_sigFi - Constructing sigFi operaters from sigFi_lookups
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_xOp_sigFi
      implicit none
      private	! except

      public :: intp_sigFi	! construct a sigFi operater

      interface intp_sigFi; module procedure	&
	intp_,		&
	intp_subB_,	&
	intp_kkNav_,	&
	intp_var1_
      end interface

! !REVISION HISTORY:
! 	04Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_xOp_sigFi'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intp_ - construct a sigFi operator
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intp_(ndat,kts,rlat,rlon,rlev,sigFi,stat)

      use m_die,   only : die
      use m_stdio, only : stderr
      use m_SortingTools,only : indexSet, indexSort
      use m_subBlocks, only : subBlocks,init,clean

      implicit none

		! Attributes of the output sigF

      integer, intent(in) :: ndat
      integer,dimension(:), intent(in) :: kts
      real,   dimension(:), intent(in) :: rlat
      real,   dimension(:), intent(in) :: rlon
      real,   dimension(:), intent(in) :: rlev
      real,   dimension(:), intent(out) :: sigFi

      integer, optional, intent(out) :: stat

! !REVISION HISTORY:
! 	18Mar97 - Jing Guo <guo@eramus> - initial prototyping and coding
! 	04Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intp_'

  integer,dimension(:),allocatable :: iD
  real,   dimension(:),allocatable :: wk

  type(subBlocks) :: ktList
  integer :: ier


  if(present(stat)) stat=0
  if(ndat == 0) return		! nothing else happens

  allocate(iD(ndat),wk(ndat), stat=ier)
  if(ier /= 0) then
    write(stderr,'(2a,i4)') myname_,': allocate() error, stat =',ier
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  call indexSet(ndat,iD)
  call indexSort(ndat,iD,kts,descend=.false.)

	! Build a table out of a _sorted_ kts(:).  Note that kts(:)
	! as the input to init_block(), is sorted at the entry via
	! kts(iD(1:ndat)), where iD(1:ndat) is a sorting index table.

  call init(ktList,kts(iD(1:ndat)))

	! Enter the subroutine for each kt with sorted input arrays

  call intp_subB_(ktList,					&
	rlat(iD(1:ndat)),rlon(iD(1:ndat)),rlev(iD(1:ndat)),	&
	sigFi, stat=ier						)

	! _Unsort_ the result sigFi[]

  if(ier == 0) sigFi(iD(1:ndat))=sigFi(1:ndat)

  if(ier /= 0) then
    write(stderr,'(2a,i4)') myname_,	&
	': intp_subB_() error, stat =',ier
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  call clean(ktList)
  deallocate(iD,wk,stat=ier)

end subroutine intp_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intp_kkNav_ - interpolate with a segmented vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intp_kkNav_(nav,krNav,ktNav,rlat,rlon,rlev,sigFi,stat)

      use m_die,   only : die
      use m_stdio, only : stderr
      use m_Navigator,only : Navigator,lsize,get

      implicit none

      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav
      real,intent(in) :: rlat(:) ! It has been sorted w.r.t kt
      real,intent(in) :: rlon(:) ! It has been sorted w.r.t kt
      real,intent(in) :: rlev(:) ! It has been sorted w.r.t kt
      real,intent(out):: sigFi(:)	! return sig_Phi
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	04Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intp_kkNav_'

  integer :: kt,lc,ln,le
  integer :: l
  integer :: ier

  if(present(stat)) stat=0

	! This is a loop over ktList entries.  Inside the loop, the
	! elements of the working array dsigF(lc:ln) is accessed
	! without overlapping, therefore, maybe accessed in parallel.

  do l=1,lsize(nav)
    call get(nav,l,lc=lc,le=le,ln=ln)

    if(ln <= 0) cycle

    kt=ktNav(l)

	! Segment (lc:le) is a collection of the data of the same
	! variable type

    call intp_var1_(kt,rlat(lc:le),rlon(lc:le),rlev(lc:le),	&
	sigFi(lc:le), stat=ier					)

	if(ier /= 0) then
	  write(stderr,'(2a,i3,a,i3)') myname_,	&
	    ': intp_var1_(kt=',kt,') error, stat =',ier
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  end do	! l=1,lsize(nav)

end subroutine intp_kkNav_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intp_subB_ - interpolate with a segmented vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intp_subB_(ktList,rlat,rlon,rlev,sigFi,stat)

      use m_die,   only : die
      use m_stdio, only : stderr
      use m_subBlocks,only : subBlocks,lsize,blockAttr

      implicit none

      type(subBlocks),intent(in) :: ktList
      real,intent(in) :: rlat(:) ! It has been sorted w.r.t kt
      real,intent(in) :: rlon(:) ! It has been sorted w.r.t kt
      real,intent(in) :: rlev(:) ! It has been sorted w.r.t kt
      real,intent(out):: sigFi(:)	! return sig_Phi
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	04Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intp_subB_'

  integer :: kt,lc,ln,le
  integer :: l
  integer :: ier

  if(present(stat)) stat=0

	! This is a loop over ktList entries.  Inside the loop, the
	! elements of the working array dsigF(lc:ln) is accessed
	! without overlapping, therefore, maybe accessed in parallel.

  do l=1,lsize(ktList)
    call blockAttr(ktList,l,kt=kt,lc=lc,ln=ln)

    if(ln <= 0) cycle

    le=lc+ln-1

	! Segment (lc:le) is a collection of the data of the same
	! variable type

    call intp_var1_(kt,rlat(lc:le),rlon(lc:le),rlev(lc:le),	&
	sigFi(lc:le), stat=ier					)

	if(ier /= 0) then
	  write(stderr,'(2a,i3,a,i3)') myname_,	&
	    ': intp_var1_(kt=',kt,') error, stat =',ier
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  end do	! l=1,lsize_block(ktList)

end subroutine intp_subB_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intp_var1_ - interpolate for a given variable type
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intp_var1_(kt,rlat,rlon,rlev,sigFi,stat)

      use const,  only : G_EARTH  ! acc. due to the Earth gravity, m/ss
      use const,  only : R_EARTH  ! radius of Earth in meters
      use const,  only : OMEGA
      use const,  only : RHOBAR	  ! mean sea-level density

      use m_ktList, only : ktHH,ktslp,ktQQ
      use m_ktList, only : ktHm => ktUU
      use m_ktList, only : ktHl => ktVV
      use m_ktList, only : ktPm => ktUs
      use m_ktList, only : ktPl => ktVs

      use m_sigFi_lookups
      use m_normCor, only : intp_normCor
      use m_intpAP,  only : intpAP

      use m_die,   only : die
      use m_stdio, only : stderr

      implicit none

      integer,intent(in) :: kt
      real,dimension(:),intent(in) :: rlat,rlon,rlev
      real,dimension(:),intent(out):: sigFi
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	14Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intp_var1_'

  real, parameter :: scalar_p = .01*G_Earth*RHOBAR	! Pa to mb
  real, parameter :: scalar_d = G_EARTH/(2.*OMEGA*R_EARTH)

  integer :: ier,ln
  real,dimension(:),allocatable :: anorm

  if(present(stat)) stat=0

    ier = -1
    ln = size(sigFi)

    if(	ln /= size(rlat) .or.	&
	ln /= size(rlon) .or.	&
	ln /= size(rlev)	) then

      write(stderr,'(2a,i4)') myname_,': size(sigFi) = ',ln
      write(stderr,'(2a,i4)') myname_,': size(rlat ) = ',size(rlat)
      write(stderr,'(2a,i4)') myname_,': size(rlon ) = ',size(rlon)
      write(stderr,'(2a,i4)') myname_,': size(rlev ) = ',size(rlev)

      if(.not.present(stat)) call die(myname_)
      stat=ier
      return
    endif

      select case(kt)
      case(ktslp,ktPm,ktPl,ktHH,ktHm,ktHl)

        call intpAP(ln,rlat,rlon,rlev,sigFi,		&
	  hghte_dim%nlon, hghte_dim%nlat,		&
	  hghte_dim%nlev, hghte_dim%vlev,		&
	  hghte_fld,hghte_fill,	stat=ier		)

	if(ier == 0) then
	  select case(kt)
	  case(ktslp)
            sigFi(:) = scalar_p*sigFi(:)

	  case(ktHm,ktHl,ktPm,ktPl)

		! Define the normalization scalar

	    allocate(anorm(ln),stat=ier)
	    if(ier /= 0) then
	      write(stderr,'(2a,i4)') myname_,	&
		': allocate() error, stat =',ier
	      if(.not.present(stat)) call die(myname_)
	      stat=ier
	      return
	    endif

	    call intp_normCor(ln,rlev,anorm)
            sigFi(:) = scalar_d*sigFi(:)/anorm(:)

	    deallocate(anorm)
	  end select
	endif

      case(ktQQ)

        call intpAP(ln,rlat,rlon,rlev,sigFi,		&
	  mixre_dim%nlon, mixre_dim%nlat,		&
	  mixre_dim%nlev, mixre_dim%vlev,		&
	  mixre_fld,mixre_fill,	stat=ier		)

      end select

  if(ier /= 0) then
    write(stderr,'(2a,i3,a,i4)') myname_,	&
	': intpAP(',kt,') error, stat =',ier
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

end subroutine intp_var1_

!=======================================================================
end module m_xOp_sigFi
!.
