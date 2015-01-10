!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AttrVect - a distributed Innovation vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_AttrVect
      use m_List, only : List
      implicit none
      private	! except

      public :: AttrVect		! The class data structure

      public :: init,AttrVect_init	! create a local vector
      public :: clean		! clean the local vector
      public :: lsize		! the actual size of the local vector
      public :: msize		! the maximum size of the local vector
      public :: nIAttr		! number of integer attributes on local
      public :: nRAttr		! number of real attributes on local
      public :: resize		! adjust the actual size

      public :: indexIA		! index the integer attributes
      public :: indexRA		! index the real attributes

      public :: getIList	! get %iList
      public :: getRList	! get %rList
      public :: ptr_iAttr	! addressing %iAttr
      public :: ptr_rAttr	! addressing %rAttr

    type AttrVect
      private
      type(List) :: iList
      type(List) :: rList
      integer    :: lsize		! the actual size
      integer,dimension(:,:),pointer :: iAttr
      real   ,dimension(:,:),pointer :: rAttr
    end type AttrVect

    interface AttrVect_init; module procedure	&
	init__,		&
	initas__,	&
	initStr__
    end interface
    interface init   ; module procedure		&
	init__,		&
	initas__,	&
	initStr__
    end interface

    interface clean  ; module procedure clean__  ; end interface
    interface lsize  ; module procedure lsize__  ; end interface
    interface msize  ; module procedure msize__  ; end interface
    interface resize ; module procedure resize__ ; end interface
    interface nIAttr ; module procedure nIAttr__ ; end interface
    interface nRAttr ; module procedure nRAttr__ ; end interface

    interface indexIA; module procedure indexIA__; end interface
    interface indexRA; module procedure indexRA__; end interface

    interface getIList; module procedure getIList__; end interface
    interface getRList; module procedure getRList__; end interface

    interface ptr_iAttr; module procedure ptr_iAttr__; end interface
    interface ptr_rAttr; module procedure ptr_rAttr__; end interface

! !REVISION HISTORY:
! 	10Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AttrVect'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init__ - initialize with given iList, rList, and the size
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init__(aV,iList,rList,lsize)
      use m_List, only : init,nitem
      use m_die , only : die
      use m_mall, only : mall_mci,mall_ison

      implicit none
      type(AttrVect),intent(out) :: aV
      character(len=*),optional,intent(in) :: iList
      character(len=*),optional,intent(in) :: rList
      integer,         optional,intent(in) :: lsize

! !REVISION HISTORY:
! 	09Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname__=myname//'::init_'
  integer :: nIA,nRA,n,ier

  if(present(iList)) then
    call init(aV%iList,iList)	! init.List()
  else
    call init(aV%iList,'')	! init.List()
  endif

  if(present(rList)) then
    call init(aV%rList,rList)	! init.List()
  else
    call init(aV%rList,' ')	! init.List()
  endif

  nIA=nitem(aV%iList)		! nitem.List()
  nRA=nitem(aV%rList)		! nitem.List()

  n=0
  if(present(lsize)) n=lsize

  allocate( aV%iAttr(nIA,n),aV%rAttr(nRA,n),	stat=ier)
  if(ier /= 0) call die(myname__,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(aV%rAttr,myname)
	  call mall_mci(aV%iAttr,myname)
	endif

  aV%lsize=n
end subroutine init__
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initStr__ - initialize with given iList, rList, and the size
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initStr__(aV,iList,rList,lsize)
      use m_List  ,only : init
      use m_String,only : String
      use m_die   ,only : die
      use m_mall  ,only : mall_mci,mall_ison

      implicit none
      type(AttrVect),intent(out) :: aV
      type(String),dimension(:),intent(in) :: iList
      type(String),dimension(:),intent(in) :: rList
      integer,     intent(in) :: lsize

! !REVISION HISTORY:
! 	09Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname__=myname//'::initStr_'
  integer :: nIA,nRA,n,ier

  call init(aV%iList,iList)	! init.List()
  call init(aV%rList,rList)	! init.List()

  nIA=size(iList)
  nRA=size(rList)

  n=lsize

  allocate( aV%iAttr(nIA,n),aV%rAttr(nRA,n),	stat=ier)
  if(ier /= 0) call die(myname__,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(aV%rAttr,myname)
	  call mall_mci(aV%iAttr,myname)
	endif

  aV%lsize=n
end subroutine initStr__
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initas__ - initialize as the given vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initas__(aV,bV,lsize)
      use m_String, only : String,toChar,clean
      use m_List,   only : get
      implicit none
      type(AttrVect),intent(out) :: aV
      type(AttrVect),intent(in)  :: bV
      integer,optional,intent(in)  :: lsize

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::initas_'
  type(String) :: iLStr,rLStr
  integer :: isize_

	! Convert the two Lists to two Strings

  call get(iLStr,bv%iList)
  call get(rLStr,bv%rList)

  isize_=lsize__(bv)
  if(present(lsize)) isize_=lsize

  call init__(aV,iList=toChar(iLStr),rList=toChar(rLStr),lsize=isize_)

  call clean(iLStr)
  call clean(rLStr)

end subroutine initas__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean__ - clean a vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean__(aV)
      use m_die,only : die
      use m_mall,only : mall_mco,mall_ison
      use m_List, only : clean
      implicit none
      type(AttrVect),intent(inout) :: aV

! !REVISION HISTORY:
! 	09Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::clean_'
  integer :: ier

	if(mall_ison()) then
	  call mall_mco(aV%iAttr,myname)
	  call mall_mco(aV%rAttr,myname)
	endif

  deallocate(aV%iAttr,aV%rAttr,stat=ier)
  if(ier /= 0) call die(myname__,'deallocte()',ier)

  call clean(aV%iList)
  call clean(aV%rList)

  aV%lsize=0
end subroutine clean__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize__ - the actual local size of the vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    function lsize__(aV)
      implicit none
      type(AttrVect), intent(in) :: aV
      integer :: lsize__

! !REVISION HISTORY:
! 	09Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::lsize_'

  lsize__=aV%lsize

end function lsize__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: msize__ - the maximum local size of the vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    function msize__(aV)
      implicit none
      type(AttrVect), intent(in) :: aV
      integer :: msize__

! !REVISION HISTORY:
! 	09Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::msize_'

  msize__=size(aV%iAttr,2)

end function msize__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: resize__ - adjust the actual size
!
! !DESCRIPTION:
!
!	It is user's responsibility to ensure lsize is no greater than
!   the maximum size of the vector.
!
!	If lsize is not specified, the size of the vector is adjusted
!   to its original size.
!
! !INTERFACE:

    subroutine resize__(aV,lsize)
      implicit none
      type(AttrVect)  ,intent(inout) :: aV
      integer,optional,intent(in   ) :: lsize	! a new size

! !REVISION HISTORY:
! 	03May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::resize_'
  integer :: m

  m=msize__(aV)

  aV%lsize=m
  if(present(lsize)) aV%lsize=lsize

end subroutine resize__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nIAttr__ - number of INTEGER type attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    function nIAttr__(aV)
      use m_List, only : nitem
      implicit none
      type(AttrVect),intent(in) :: aV
      integer :: nIAttr__

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::nIAttr_'
  nIAttr__=nitem(aV%iList)

end function nIAttr__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nRAttr__ - number of REAL type attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    function nRAttr__(aV)
      use m_List, only : nitem
      implicit none
      type(AttrVect),intent(in) :: aV
      integer :: nRAttr__

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::nRAttr_'
  nRAttr__=nitem(aV%rList)

end function nRAttr__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getIList__ - get an item from iList
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getIList__(item,ith,aVect)
      use m_String, only : String
      use m_List,   only : get
      implicit none
      type(String),intent(out) :: item
      integer,     intent(in)  :: ith
      type(AttrVect),intent(in) :: aVect

! !REVISION HISTORY:
! 	24Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::getIList_'

  call get(item,ith,aVect%iList)

end subroutine getIList__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getRList__ - get an item from rList
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getRList__(item,ith,aVect)
      use m_String, only : String
      use m_List,   only : get
      implicit none
      type(String),intent(out) :: item
      integer,     intent(in)  :: ith
      type(AttrVect),intent(in) :: aVect

! !REVISION HISTORY:
! 	24Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::getRList_'

  call get(item,ith,aVect%rList)

end subroutine getRList__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexIA__ - index the integer attribute List
!
! !DESCRIPTION:
!
! !INTERFACE:

    function indexIA__(aV,item)
      use m_List, only : index
      implicit none
      type(AttrVect), intent(in) :: aV
      character(len=*),intent(in) :: item
      integer :: indexIA__

! !REVISION HISTORY:
! 	27Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::indexIA_'

  indexIA__=index(aV%iList,item)

end function indexIA__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexRA__ - index the integer attribute List
!
! !DESCRIPTION:
!
! !INTERFACE:

    function indexRA__(aV,item)
      use m_List, only : index
      implicit none
      type(AttrVect), intent(in) :: aV
      character(len=*),intent(in) :: item
      integer :: indexRA__

! !REVISION HISTORY:
! 	27Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::indexRA_'

  indexRA__=index(aV%rList,item)

end function indexRA__

! Sorting:
!
!	aV%iVect(:,:) =		&
!		aV%iVect((/(indx(i),i=1,lsize(aV))/),:)
!
!	aV%iVect((/(indx(i),i=1,lsize(aV))/),:) =	&
!		aV%iVect(:,:)
!
!	aV%iVect(:,ikx),aV%iVect(:,iks)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_iAttr__ - addressing of %iAttr
!
! !DESCRIPTION:
!
!	It is the user's responsibility to ensure lbnd:ubnd are in
!   the range of iAttr(:,:)
!
! !INTERFACE:

    function ptr_iAttr__(aVect,lbnd,ubnd)
      implicit none
      type(AttrVect),intent(in) :: aVect
      integer,pointer,dimension(:,:) :: ptr_iAttr__
      integer,optional,intent(in) :: lbnd,ubnd

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::ptr_iAttr_'
  integer :: lc,le,ln

  if(present(lbnd).or.present(ubnd)) then
    lc=lbound(aVect%iAttr,2)
    if(present(lbnd)) lc=lbnd
    le=ubound(aVect%iAttr,2)
    if(present(ubnd)) le=ubnd
    ptr_iAttr__ => aVect%iAttr(:,lc:le)	! lbound()==1
  else
    ln=aVect%lsize
    ptr_iAttr__ => aVect%iAttr(:,1:ln)	! lbound()==lbound(aVect%iAttr)
  endif

end function ptr_iAttr__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_rAttr__ - addressing of %rAttr
!
! !DESCRIPTION:
!
!	It is the user's responsibility to ensure lbnd:ubnd are in
!   the range of iAttr(:,:)
!
! !INTERFACE:

    function ptr_rAttr__(aVect,lbnd,ubnd)
      implicit none
      type(AttrVect),intent(in) :: aVect
      real,pointer,dimension(:,:) :: ptr_rAttr__
      integer,optional,intent(in) :: lbnd,ubnd

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname__=myname//'::ptr_rAttr_'
  integer :: lc,le,ln

  if(present(lbnd).or.present(ubnd)) then
    lc=lbound(aVect%rAttr,2)
    if(present(lbnd)) lc=lbnd
    le=ubound(aVect%rAttr,2)
    if(present(ubnd)) le=ubnd
    ptr_rAttr__ => aVect%rAttr(:,lc:le)	! lbound()==1
  else
    ln=aVect%lsize
    ptr_rAttr__ => aVect%rAttr(:,1:ln)	! lbound()==lbound(aVect%rAttr)
  endif

end function ptr_rAttr__

end module m_AttrVect
!.
