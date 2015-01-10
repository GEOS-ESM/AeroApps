!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_LVGrid - vertical grid combined with variable types
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_LVGrid
      implicit none
      private	! except

      public :: LVGrid		! The class data structure
      public :: LVGrid_init	! initialize a (LVGrid)
      public :: LVGrid_clean	! clean a (LVGrid)
      public :: LVGrid_index	! locate a kt or kt-kp from a (LVGrid)
      public :: LVGrid_range	! find the range of all or a given kt
      public :: LVGrid_size	! find the full size of LVGrid

    type LVGrid
      private
      integer :: nLV
      integer :: nkt
      integer :: mlev

      integer,pointer,dimension(:) :: kt	! size() is nLV
      integer,pointer,dimension(:) :: kp	! size() is nLV

      integer,pointer,dimension(:) :: nlev	! numbers of levels
      integer,pointer,dimension(:) :: ilev	! level locations

      real,   pointer,dimension(:) :: levs	! size() is mlev
    end type LVGrid

    interface LVGrid_init;  module procedure	&
      init_,	&
      init2_
    end interface
    interface LVGrid_clean; module procedure clean_; end interface
    interface LVGrid_index; module procedure	&
      index2_,	&
      index3_
    end interface
    interface LVGrid_range; module procedure range_; end interface
    interface LVGrid_size ; module procedure size_ ; end interface
!________________________

    public :: ptr_kt
    public :: ptr_kp
    public :: ptr_levs

    interface ptr_kt  ; module procedure ptr_kt_  ; end interface
    interface ptr_kp  ; module procedure ptr_kp_  ; end interface
    interface ptr_levs; module procedure ptr_levs_; end interface
!________________________

    public :: new	! create an object as a pointer
    public :: delete	! delete an object as a pointer

    interface new   ; module procedure new_   ; end interface
    interface delete; module procedure delete_; end interface

! !REVISION HISTORY:
! 	04Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_LVGrid'
contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - make a level times variable grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(lvG,mlev,plevs,nkt,ktList,nlev)
      use m_ktList,only : ktList_dim
      use m_mall,only : mall_ison,mall_mci
      use m_die, only : die
      implicit none
      type(LVGrid),intent(out) :: lvG
      integer,intent(in) :: mlev
      real,dimension(:),intent(in) :: plevs
      integer,intent(in) :: nkt
      integer,dimension(:),intent(in) :: ktList		! (nkt)
      integer,dimension(:),optional,intent(in) :: nlev	! (nkt)

! !REVISION HISTORY:
! 	04Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ndim,nlv,ikt,kt,k,ier,llev

	! Get a count of "level" grid points

  nlv=0
  do ikt=1,nkt
    ndim=ktList_dim(ktList(ikt))

    select case(ndim)
    case(2)
      llev=0
      if(present(nlev)) llev=nlev(ikt)
      if(llev/=0) call die(myname_,	&
		'invalid nlev for a 2-d variable',llev)

      nlv=nlv+1

    case(3)
      llev=mlev
      if(present(nlev)) llev=nlev(ikt)

      nlv=nlv+llev

    case default
      call die(myname_,'invalide kt',kt)
    end select
  end do

  lvG%nlv=nlv
  lvG%nkt=nkt
  lvG%mlev=mlev

  allocate(lvG%kt(nlv),lvG%kp(nlv),			&
	lvG%nlev(nkt),lvG%ilev(nkt),lvG%levs(mlev),	&
	stat=ier					)

	if(ier/=0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(lvG%kt,myname)
	  call mall_mci(lvG%kp,myname)
	  call mall_mci(lvG%nlev,myname)
	  call mall_mci(lvG%ilev,myname)
	  call mall_mci(lvG%levs,myname)
	endif

	! Save or "make" the tables of the grid values

  lvG%levs(1:mlev)=plevs(1:mlev)

  nlv=0
  do ikt=1,nkt
    kt=ktList(ikt)
    ndim=ktList_dim(kt)

    select case(ndim)
    case(2)

      lvG%nlev(ikt)=0
      lvG%ilev(ikt)=nlv+1

      nlv=nlv+1
      lvG%kp(nlv)=-1
      lvG%kt(nlv)=kt

    case(3)
      llev=mlev
      if(present(nlev)) llev=nlev(ikt)

      lvG%nlev(ikt)=llev
      lvG%ilev(ikt)=nlv+1

      do k=1,llev
	nlv=nlv+1
	lvG%kt(nlv)=kt
	lvG%kp(nlv)=k
      end do

    case default
      call die(myname_,'invalide kt',kt)
    end select
  end do

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init2_ - make a level times variable grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init2_(lvG,nkt,ktList)
      use m_ktList,only : ktList_dim
      use m_mall,only : mall_ison,mall_mci
      use m_die, only : die
      implicit none
      type(LVGrid),intent(out) :: lvG
      integer,intent(in) :: nkt
      integer,dimension(:),intent(in) :: ktList

! !REVISION HISTORY:
! 	04Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init2_'
  integer :: ndim,nlv,ikt,kt,ier

	! Get a count of "level" grid points

  nlv=0
  do ikt=1,nkt
    ndim=ktList_dim(ktList(ikt))

    select case(ndim)
    case(2)
      nlv=nlv+1
    case default
      call die(myname_,'invalide kt',kt)
    end select
  end do

  lvG%nlv=nlv
  lvG%mlev=0

  allocate(lvG%kt(nlv),lvG%kp(nlv),	&
	lvG%nlev(nkt),lvG%ilev(nkt),lvG%levs(0),stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(lvG%kt,myname)
	  call mall_mci(lvG%kp,myname)
	  call mall_mci(lvG%nlev,myname)
	  call mall_mci(lvG%ilev,myname)
	  call mall_mci(lvG%levs,myname)
	endif

	! Save or "make" the tables of the grid values

  nlv=0
  do ikt=1,nkt
    kt=ktList(ikt)
    ndim=ktList_dim(kt)

    select case(ndim)
    case(2)
      lvG%nlev(ikt)=0
      lvG%ilev(ikt)=nlv+1

      nlv=nlv+1
      lvG%kp(nlv)=-1
      lvG%kt(nlv)=kt

    case default
      call die(myname_,'invalide kt',kt)
    end select
  end do

end subroutine init2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a levelXvariable grid table
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(lvG)
      use m_mall,only : mall_ison,mall_mco
      use m_die, only : die
      implicit none
      type(LVGrid),intent(inout) :: lvG

! !REVISION HISTORY:
! 	04Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

	if(mall_ison()) then
	  call mall_mco(lvG%kt,myname)
	  call mall_mco(lvG%kp,myname)
	  call mall_mco(lvG%nlev,myname)
	  call mall_mco(lvG%ilev,myname)
	  call mall_mco(lvG%levs,myname)
	endif

  deallocate(lvG%kt,lvG%kp,lvG%nlev,lvG%ilev,lvG%levs,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

  lvG%nlv=0
  lvG%mlev=0

end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: size_ - return the total size
!
! !DESCRIPTION:
!
! !INTERFACE:

    function size_(lvG,mlev)
      implicit none
      type(LVGrid),intent(in) :: lvG
      integer,optional,intent(in) :: mlev
      integer :: size_

! !REVISION HISTORY:
! 	25Feb00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::size_'

  size_=lvG%nLV
  if(present(mlev)) size_=lvG%mlev

end function size_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: range_ - find the storage range of a given kt
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine range_(lvG,kt,lower,upper)
      use m_ktList,only : ktList_dim
      implicit none
      type(LVGrid),    intent(in)  :: lvG
      integer,         intent(in)  :: kt
      integer,         intent(out) :: lower,upper

! !REVISION HISTORY:
! 	04Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::range_'
  integer :: lv

  lower=1	! a non-existent variable corresponds to a zero-sized
  upper=0	! array.

  do lv=1,lvG%nkt
    if(lvG%kt(lv) == kt) then
      lower=lvG%ilev(lv)
      upper=lower+lvG%nlev(lv)-1
      exit
    endif
  end do

end subroutine range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: index2_ - index a given value of kt-lev
!
! !DESCRIPTION:
!
! !INTERFACE:

    function index2_(lvG,kt)
      use m_ktList,only : ktList_dim
      implicit none
      type(LVGrid),intent(in) :: lvG
      integer,intent(in) :: kt
      integer :: index2_

! !REVISION HISTORY:
! 	04Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::index2_'
  integer :: lv,ndim

  index2_=0

  ndim=ktList_dim(kt)
  select case(ndim)
  case(2)
    do lv=1,lvG%nlv
      if(lvG%kt(lv) == kt .and. lvG%kp(lv)<=0) then
	index2_=lv
	return
      endif
    end do
  end select

end function index2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: index3_ index a give value of kt-lev
!
! !DESCRIPTION:
!
! !INTERFACE:

    function index3_(lvG,kt,kp)
      use m_ktList,only : ktList_dim
      implicit none
      type(LVGrid),intent(in) :: lvG
      integer,intent(in) :: kt
      integer,intent(in) :: kp
      integer :: index3_

! !REVISION HISTORY:
! 	04Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::index3_'
  integer :: lv,ndim

  index3_=0

  if(kp<=0) return

  ndim=ktList_dim(kt)
  select case(ndim)
  case(3)
    do lv=1,lvG%nlv
      if(lvG%kt(lv)==kt .and. lvG%kp(lv)==kp) then
	index3_=lv
	return
      endif
    end do
  end select

end function index3_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: new_ - create an object as a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    function new_(mold,stat)
      use m_die, only : die,perr
      use m_mall,only : mall_ison,mall_ci
      implicit none
      type(LVGrid),pointer :: mold
      integer,optional,intent(out) :: stat
      type(LVGrid),pointer :: new_

! !REVISION HISTORY:
! 	18Sep01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		. initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::new_'
  type(LVGrid),pointer :: obj
  integer :: ier

  if(present(stat)) stat=0
  allocate(obj,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'allocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	if(mall_ison()) call mall_ci(1,myname)

  new_ => obj
  nullify(obj)		! to prevent the compiler touching the memory.
end function new_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: delete_ - delete an object as a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine delete_(obj,stat)
      use m_die, only : die,perr
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(LVGrid),pointer :: obj
      integer,optional,intent(out) :: stat


! !REVISION HISTORY:
! 	18Sep01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		. initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::delete_'
  integer :: ier

  if(present(stat)) stat=0

	if(mall_ison()) call mall_co(1,myname)

  deallocate(obj,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine delete_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_kt_ - Component %kt of all LV grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_kt_(obj)
      implicit none
      type(LVGrid),intent(in) :: obj
      integer,pointer,dimension(:) :: ptr_kt_

! !REVISION HISTORY:
! 	18Sep01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_kt_'

  ptr_kt_ => obj%kt(:)
end function ptr_kt_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_kp_ - Component %kp of all LV grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_kp_(obj)
      implicit none
      type(LVGrid),intent(in) :: obj
      integer,pointer,dimension(:) :: ptr_kp_

! !REVISION HISTORY:
! 	18Sep01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_kp_'

  ptr_kp_ => obj%kp(:)
end function ptr_kp_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_levs_ - Component %lev indexed by values of %kp
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_levs_(obj)
      implicit none
      type(LVGrid),intent(in) :: obj
      real,pointer,dimension(:) :: ptr_levs_

! !REVISION HISTORY:
! 	18Sep01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_levs_'

  ptr_levs_ => obj%levs(:)
end function ptr_levs_

end module m_LVGrid
