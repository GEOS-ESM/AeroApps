!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_CorAttrF - Attributes of <Phi,Phi> correlation
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_CorAttrF
      use m_Navigator,only : Navigator
      use m_Collector,only : Collector

      implicit none
      private	! except

      public :: CorAttrF		! The class data structure
      public :: CorAttrF_init,init	! initialize an object
      public :: clean			! clean an object

      public :: ptr_Navigator 
      public :: ptr_krNav
      public :: ptr_ktNav
      public :: ptr_qr
      public :: ptr_qd
      public :: ptr_qm
      public :: ptr_ql
      public :: ptr_kl

      public :: ptr_collNav
      public :: ptr_collVec

    type CorAttrF
      private

		! The Navigator for the "scattered" attributes

      type(Navigator),pointer :: nav
      integer,pointer,dimension(:) :: krNav
      integer,pointer,dimension(:) :: ktNav

		! The "scattered" attributes

      real   ,pointer,dimension(:,:) :: qr
      real   ,pointer,dimension(:,:) :: qd
      real   ,pointer,dimension(:,:) :: qm
      real   ,pointer,dimension(:,:) :: ql
      integer,pointer,dimension(  :) :: kl

		! The Navigator for the "all-gathered" attributes

      type(Navigator),pointer :: nav_all
      integer,pointer,dimension(:) :: krNav_all
      integer,pointer,dimension(:) :: ktNav_all

		! The "all-gathered" attributes

      real   ,pointer,dimension(:,:) :: qr_all
      real   ,pointer,dimension(:,:) :: qd_all
      integer,pointer,dimension(  :) :: kl_all

		! Collectors

      type(Collector),pointer :: collNav ! the collector of navigators
      type(Collector),pointer :: collVec ! the collector of attributes

    end type CorAttrF

      interface CorAttrF_init; module procedure	&
	init_; end interface
      interface init; module procedure	&
	init_; end interface
      interface clean; module procedure clean_; end interface

      interface ptr_Navigator; module procedure	&
	ptr_Navigator_; end interface
      interface ptr_krNav; module procedure	&
	ptr_krNav_; end interface
      interface ptr_ktNav; module procedure	&
	ptr_ktNav_; end interface
      interface ptr_qr; module procedure	&
	ptr_qr_; end interface
      interface ptr_qd; module procedure	&
	ptr_qd_; end interface
      interface ptr_qm; module procedure	&
	ptr_qm_; end interface
      interface ptr_ql; module procedure	&
	ptr_ql_; end interface
      interface ptr_kl; module procedure	&
	ptr_kl_; end interface

      interface ptr_collNav; module procedure	&
	ptr_collNav_; end interface
      interface ptr_collVec; module procedure	&
	ptr_collVec_; end interface

! !REVISION HISTORY:
! 	06Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_CorAttrF'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initqvec_ - initialize spherical unit vectors
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initqvec_(qr,qm,ql,qd,nav,ktNav,rlon,rlat)
      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_die,only : die
      implicit none

      real,dimension(:,:),intent(out) :: qr
      real,dimension(:,:),intent(out) :: qm
      real,dimension(:,:),intent(out) :: ql
      real,dimension(:,:),intent(out) :: qd

      type(Navigator)     ,intent(in)  :: nav
      integer,dimension(:),intent(in)  :: ktNav
      real   ,dimension(:),intent(in)  :: rlon
      real   ,dimension(:),intent(in)  :: rlat

! !REVISION HISTORY:
! 	06Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initqvec_'
  include "ktmax.h"
  integer,parameter :: ktHm=ktUU
  integer,parameter :: ktHl=ktVV
  integer,parameter :: ktPm=ktUs
  integer,parameter :: ktPl=ktVs

  integer :: inav
  integer :: lc,le,kt

  do inav=1,lsize(nav)
    call get(nav,inav,lc=lc,le=le)

	! Horizontal distance indexing through qr

    call setqr_(rlon(lc:le),rlat(lc:le),qr(:,lc:le))
    call setqm_(rlon(lc:le),rlat(lc:le),qm(:,lc:le))
    call setql_(rlon(lc:le),rlat(lc:le),ql(:,lc:le))

    kt=ktNav(inav)

    select case(kt)
    case(ktHH,ktslp,ktQQ)
      qd(:,lc:le)=0.

    case(ktHm,ktPm)
      qd(:,lc:le)=qm(:,lc:le)

    case(ktHl,ktPl)
      qd(:,lc:le)=ql(:,lc:le)

    case default
      call die(myname_,'unknown kt',kt)
    end select
  end do

end subroutine initqvec_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize the object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(qattr,nav,krNav,ktNav,rlon,rlat,rlev,comm)

      use m_mall,only : mall_ison,mall_mci
      use m_die ,only : die

      use m_Navigator,only : Navigator
      use m_Navigator,only : new
      use m_Navigator,only : lsize
      use m_Navigator,only : allgather

      use m_Collector,only : Collector_init
      use m_Collector,only : new
      use m_Collector,only : globalSize
      use m_Collector,only : get

      use m_CollectorComm,only : allgatherv

      use m_xTab_levs,only : index_levs

      use m_zeit, only : zeit_ci, zeit_co

      implicit none

      type(CorAttrF)       ,intent(out) :: qattr

      type(Navigator),target,intent(in)  :: nav
      integer,target,dimension(:),intent(in) :: krNav
      integer,target,dimension(:),intent(in) :: ktNav
      real   ,dimension(:),intent(in)  :: rlon
      real   ,dimension(:),intent(in)  :: rlat
      real   ,dimension(:),intent(in)  :: rlev
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	06Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

  integer :: nsct,ngat
  integer :: lsct,lgat,ld,lc,le
  integer :: ier

  call zeit_ci('attr_1')
	! Note that nxxx are in general counters for Navigators, while
	! lxxx are in general counters for vectors or attributes.

  qattr%nav   => nav	! make shallow copies
  qattr%krNav => krNav
  qattr%ktNav => ktNav

		! define the navigator Collector, collnav

  nsct=lsize(nav)	! the size of the local navigator

  	qattr%collnav => new(qattr%collnav)
  call Collector_init(qattr%collnav,nsct,comm)

  call zeit_co('attr_1')
  call zeit_ci('attr_2')

		! define the vector Collector, collvec

  lsct=size(rlon)	! the size of the local vector

  		! the vector size is specified by the input argument.
		! Note that all sizes must be consistent.  Specifically,
		! the data in argument _nav_ must be consistant with
		! the value of _lsct_ to avoid posible memory violation.


  	qattr%collvec => new(qattr%collvec)
  call Collector_init(qattr%collvec,lsct,comm)
  call zeit_co('attr_2')
  call zeit_ci('attr_3')

		! get the global location references of the local
		! segement.  displ=ld is used later as the difference
		! between the local location references and global
		! location references.

  call get(qattr%collvec,displ=ld,lbound=lc,ubound=le)

		! Collect local navigators to define an "all-gathered"
		! navigator

  	qattr%nav_all => new(qattr%nav_all)

  call allgather(nav,ld,qattr%nav_all,qattr%collnav,comm)
  call zeit_co('attr_3')
  call zeit_ci('attr_4')

	! Navigator attributes allocated

  ngat=globalSize(qattr%collnav)
  allocate(qattr%krNav_all(ngat),qattr%ktNav_all(ngat),stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)
	if(mall_ison()) then
	  call mall_mci(qattr%krNav_all,myname)
	  call mall_mci(qattr%ktNav_all,myname)
	endif

  call allgatherv(krNav,qattr%krNav_all,qattr%collnav,comm)
  call allgatherv(ktNav,qattr%ktNav_all,qattr%collnav,comm)

		! define the buffer for "all-gathered" attrbutes

  lgat=globalSize(qattr%collvec)

  allocate( qattr%qr_all(3,lgat),	&
	    qattr%qd_all(3,lgat),	&
	    qattr%kl_all(  lgat),	&
	    qattr%qm(3,lsct),		&
	    qattr%ql(3,lsct),	stat=ier)

	if(ier/=0) call die(myname_,'allocate()',ier)
	if(mall_ison()) then
	  call mall_mci(qattr%qr_all,myname)
	  call mall_mci(qattr%qd_all,myname)
	  call mall_mci(qattr%kl_all,myname)
	  call mall_mci(qattr%qm    ,myname)
	  call mall_mci(qattr%ql    ,myname)
	endif

	! Localize the "local" segment

  qattr%qr => qattr%qr_all(:,lc:le)
  qattr%qd => qattr%qd_all(:,lc:le)
  qattr%kl => qattr%kl_all(  lc:le)
  call zeit_co('attr_4')

	! Define spherical coordinate unit vectors to the local
	! segment.  The pointers are used as buffers

  call zeit_ci('attr_5')
  call initqvec_(qattr%qr,qattr%qm,qattr%ql,qattr%qd,	&
	nav,ktNav,rlon,rlat)
  call zeit_co('attr_5')

	! This approach for level indexing must be changed in the
	! future obs.op. version.  Indexing should be done against
	! individual level tables for different covariance components
	! and even the tables for different correlation components.

  call zeit_ci('attr_6')
  call index_levs(rlev(1:lsct),qattr%kl)
  call zeit_co('attr_6')

		! define the "all-gathered" attributes.  Note that
		! local segment pointers are put in "()" to ensure
		! their "intent(in)" status.

  call zeit_ci('attr_7')
  call allgatherv((qattr%qr),qattr%qr_all,qattr%collvec,comm)
  call allgatherv((qattr%qd),qattr%qd_all,qattr%collvec,comm)
  call allgatherv((qattr%kl),qattr%kl_all,qattr%collvec,comm)
  call zeit_co('attr_7')

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(qattr)
      use m_mall,only : mall_ison,mall_mco
      use m_die ,only : die
      use m_Collector,only : clean
      use m_Collector,only : delete
      use m_Navigator,only : clean
      use m_Navigator,only : delete
      implicit none
      type(CorAttrF),intent(inout) :: qattr

! !REVISION HISTORY:
! 	06Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  call  clean(qattr%collnav)
  call delete(qattr%collnav)

  call  clean(qattr%collvec)
  call delete(qattr%collvec)

  call  clean(qattr%nav_all)
  call delete(qattr%nav_all)

	if(mall_ison()) then
	  call mall_mco(qattr%qr_all   ,myname)
	  call mall_mco(qattr%qd_all   ,myname)
	  call mall_mco(qattr%kl_all   ,myname)
	  call mall_mco(qattr%qm       ,myname)
	  call mall_mco(qattr%ql       ,myname)
	  call mall_mco(qattr%krNav_all,myname)
	  call mall_mco(qattr%ktNav_all,myname)
	endif

  deallocate( qattr%qr_all,	&
	      qattr%qd_all,	&
	      qattr%kl_all,	&
	      qattr%qm    ,	&
	      qattr%ql    ,	&
	      qattr%krNav_all,	&
	      qattr%ktNav_all,	stat=ier)

	if(ier/=0) call die(myname_,'deallocate()',ier)

  nullify(qattr%nav)
  nullify(qattr%krNav)
  nullify(qattr%ktNav)

  nullify(qattr%qr)
  nullify(qattr%qd)
  nullify(qattr%kl)

end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setqr_ - define qr(x,y,z)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setqr_(rlon,rlat,qr)
      implicit none
      real,dimension(:),intent(in ) :: rlon
      real,dimension(:),intent(in ) :: rlat
      real,dimension(:,:),intent(out) :: qr

! !REVISION HISTORY:
! 	06Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setqr_'
  real :: deg
  real :: xlon,xlat
  real :: coslon,coslat
  real :: sinlon,sinlat
  integer :: i

  deg=4.*atan(1.)/180.

  do i=1,size(qr,2)

   xlon=rlon(i)*deg
   xlat=rlat(i)*deg

   coslon=cos(xlon)
   coslat=cos(xlat)
   sinlon=sin(xlon)
   sinlat=sin(xlat)

   qr(1,i)=coslat*coslon	! qr_x
   qr(2,i)=coslat*sinlon	! qr_y
   qr(3,i)=sinlat		! qr_z

  end do

end subroutine setqr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setqm_ - define qm(x,y,z)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setqm_(rlon,rlat,qm)
      implicit none
      real,dimension(:),intent(in ) :: rlon
      real,dimension(:),intent(in ) :: rlat
      real,dimension(:,:),intent(out) :: qm

! !REVISION HISTORY:
! 	06Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setqm_'

  real :: deg
  real :: xlon,xlat
  real :: coslon,coslat
  real :: sinlon,sinlat
  integer :: i

  deg=4.*atan(1.)/180.

  do i=1,size(qm,2)

   xlon=rlon(i)*deg
   xlat=rlat(i)*deg

   coslon=cos(xlon)
   coslat=cos(xlat)
   sinlon=sin(xlon)
   sinlat=sin(xlat)

   qm(1,i)=-sinlat*coslon	! qm_x
   qm(2,i)=-sinlat*sinlon	! qm_y
   qm(3,i)= coslat		! qm_z

  end do

end subroutine setqm_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setql_ - define ql(x,y,z)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setql_(rlon,rlat,ql)
      implicit none
      real,dimension(:),intent(in ) :: rlon
      real,dimension(:),intent(in ) :: rlat
      real,dimension(:,:),intent(out) :: ql

! !REVISION HISTORY:
! 	06Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setql_'

  real :: deg
  real :: xlon,xlat
  real :: coslon,coslat
  real :: sinlon,sinlat
  integer :: i

  deg=4.*atan(1.)/180.

  do i=1,size(ql,2)

   xlon=rlon(i)*deg
   xlat=rlat(i)*deg

   coslon=cos(xlon)
   coslat=cos(xlat)
   sinlon=sin(xlon)
   sinlat=sin(xlat)

   ql(1,i)=-sinlon		! ql_x
   ql(2,i)= coslon		! ql_y
   ql(3,i)= 0.			! ql_z

  end do

end subroutine setql_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_Navigator_ - referencing %navigator
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_Navigator_(qattr,whole)
      use m_Navigator,only : Navigator
      implicit none
      type(CorAttrF)  ,intent(in) :: qattr
      logical,optional,intent(in) :: whole
      type(Navigator),pointer :: ptr_Navigator_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_Navigator_'

  ptr_Navigator_ => qattr%nav

  if(present(whole)) then
    if(whole) ptr_Navigator_ => qattr%nav_all
  endif

end function ptr_Navigator_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_krNav_ - referencing %krNav
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_krNav_(qattr,whole)
      implicit none
      type(CorAttrF)  ,intent(in) :: qattr
      logical,optional,intent(in) :: whole
      integer,pointer,dimension(:) :: ptr_krNav_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_krNav_'

  ptr_krNav_ => qattr%krNav
  if(present(whole)) then
    if(whole) ptr_krNav_ => qattr%krNav_all
  endif

end function ptr_krNav_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_ktNav_ - referencing %ktNav
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_ktNav_(qattr,whole)
      implicit none
      type(CorAttrF)  ,intent(in) :: qattr
      logical,optional,intent(in) :: whole
      integer,pointer,dimension(:) :: ptr_ktNav_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_ktNav_'

  ptr_ktNav_ => qattr%ktNav
  if(present(whole)) then
    if(whole) ptr_ktNav_ => qattr%ktNav_all
  endif

end function ptr_ktNav_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_qr_ - referencing %qr(3,:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_qr_(qattr,whole)
      implicit none
      type(CorAttrF)  ,intent(in) :: qattr
      logical,optional,intent(in) :: whole
      real,pointer,dimension(:,:) :: ptr_qr_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_qr_'

  ptr_qr_ => qattr%qr
  if(present(whole)) then
    if(whole) ptr_qr_ => qattr%qr_all
  endif

end function ptr_qr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_qd_ - referencing %qd(3,:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_qd_(qattr,whole)
      implicit none
      type(CorAttrF)  ,intent(in) :: qattr
      logical,optional,intent(in) :: whole
      real,pointer,dimension(:,:) :: ptr_qd_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_qd_'

  ptr_qd_ => qattr%qd
  if(present(whole)) then
    if(whole) ptr_qd_ => qattr%qd_all
  endif

end function ptr_qd_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_qm_ - referencing %qm(3,:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_qm_(qattr)
      implicit none
      type(CorAttrF)  ,intent(in) :: qattr
      real,pointer,dimension(:,:) :: ptr_qm_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_qm_'

  ptr_qm_ => qattr%qm

end function ptr_qm_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_ql_ - referencing %ql(3,:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_ql_(qattr)
      implicit none
      type(CorAttrF)  ,intent(in) :: qattr
      real,pointer,dimension(:,:) :: ptr_ql_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_ql_'

  ptr_ql_ => qattr%ql

end function ptr_ql_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_kl_ - referencing %kl(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_kl_(qattr,whole)
      implicit none
      type(CorAttrF)  ,intent(in) :: qattr
      logical,optional,intent(in) :: whole
      integer,pointer,dimension(:) :: ptr_kl_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_kl_'

  ptr_kl_ => qattr%kl
  if(present(whole)) then
    if(whole) ptr_kl_ => qattr%kl_all
  endif

end function ptr_kl_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_collNav_ - referencing %collNav(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_collNav_(qattr)
      implicit none
      type(CorAttrF),intent(in) :: qattr
      type(Collector),pointer :: ptr_collNav_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_collNav_'

  ptr_collNav_ => qattr%collNav

end function ptr_collNav_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_collVec_ - referencing %collVec(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_collVec_(qattr)
      implicit none
      type(CorAttrF),intent(in) :: qattr
      type(Collector),pointer :: ptr_collVec_

! !REVISION HISTORY:
! 	08Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_collVec_'

  ptr_collVec_ => qattr%collVec

end function ptr_collVec_

end module m_CorAttrF
