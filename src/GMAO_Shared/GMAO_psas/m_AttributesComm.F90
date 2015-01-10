!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AttributesComm - communication methods of Attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_AttributesComm
      implicit none
      private	! except

      public :: gatherv
      public :: allgatherv
      public :: indexedTrans
      public :: indexedUntrans
      public :: trans
      public :: untrans

      interface gatherv; module procedure gatherv__; end interface
      interface allgatherv; module procedure allgatherv__; end interface
      interface trans; module procedure	&
	transposev__; end interface
      interface untrans; module procedure	&
	untransposev__; end interface
      interface indexedTrans; module procedure	&
	indexedtransposev__; end interface
      interface indexedUntrans; module procedure	&
	indexeduntransposev__; end interface

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AttributesComm'
contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gatherv__ - gather a subset Attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gatherv__(sattr,rattr,coll,root,comm,stat)
      use m_Attributes,only : Attributes
      use m_Attributes,only : Attributes_init
      use m_Attributes,only : key
      use m_Attributes,only : ptr_kr
      use m_Attributes,only : ptr_ks
      use m_Attributes,only : ptr_kx
      use m_Attributes,only : ptr_kt
      use m_Attributes,only : ptr_lat
      use m_Attributes,only : ptr_lon
      use m_Attributes,only : ptr_lev

      use m_Attributes,only : KR
      use m_Attributes,only : KS
      use m_Attributes,only : KX
      use m_Attributes,only : KT
      use m_Attributes,only : LAT
      use m_Attributes,only : LON
      use m_Attributes,only : LEV

      use m_Collector    ,only : Collector
      use m_Collector    ,only : globalSize
      use m_CollectorComm,only : gatherv

      use m_die,only : die,perr
      implicit none

      type(Attributes),intent(in ) :: sattr	! send buffer
      type(Attributes),intent(out) :: rattr	! receive buffer
      type(Collector) ,intent(in)  :: coll	! a Collector
      integer,intent(in) :: root		! root PE
      integer,intent(in) :: comm		! communicator
      integer,optional,intent(out) :: stat	! return code

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gatherv__'
  integer :: kset
  integer :: gsize_
  integer :: ier
  integer,pointer,dimension(:) :: iattr
  real   ,pointer,dimension(:) :: xattr

	! This is the code for the subset Attributes

  if(present(stat)) stat=0

  kset=key(sattr)

	! Define a skeleton of Attributes

  gsize_=globalSize(coll)
  call Attributes_init(rattr,kset,gsize_)

	! Start message-passing

  ier=-1

  if(iand(kset,KR)/=0) then
	iattr => ptr_kr(rattr)
	call gatherv(ptr_kr(sattr),iattr,coll,root,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'gatherv(kr)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif
  if(iand(kset,KT)/=0) then
	iattr => ptr_kt(rattr)
	call gatherv(ptr_kt(sattr),iattr,coll,root,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'gatherv(kt)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif
  if(iand(kset,KX)/=0) then
	iattr => ptr_kx(rattr)
	call gatherv(ptr_kx(sattr),iattr,coll,root,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'gatherv(kx)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif
  if(iand(kset,KS)/=0) then
	iattr => ptr_ks(rattr)
	call gatherv(ptr_ks(sattr),iattr,coll,root,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'gatherv(ks)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif
  if(iand(kset,LAT)/=0) then
	xattr => ptr_lat(rattr)
	call gatherv(ptr_lat(sattr),xattr,coll,root,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'gatherv(lat)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif
  if(iand(kset,LON)/=0) then
	xattr => ptr_lon(rattr)
	call gatherv(ptr_lon(sattr),xattr,coll,root,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'gatherv(lon)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif
  if(iand(kset,LEV)/=0) then
	xattr => ptr_lev(rattr)
	call gatherv(ptr_lev(sattr),xattr,coll,root,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'gatherv(lev)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif

  	if(ier/=0) then
	  call perr(myname_,'invalid subset',kset)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine gatherv__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgatherv__ - allgather a subset Attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgatherv__(sattr,rattr,coll,comm,stat)
      use m_Attributes,only : Attributes
      use m_Attributes,only : Attributes_init
      use m_Attributes,only : key
      use m_Attributes,only : ptr_kr
      use m_Attributes,only : ptr_ks
      use m_Attributes,only : ptr_kx
      use m_Attributes,only : ptr_kt
      use m_Attributes,only : ptr_lat
      use m_Attributes,only : ptr_lon
      use m_Attributes,only : ptr_lev

      use m_Attributes,only : KR
      use m_Attributes,only : KS
      use m_Attributes,only : KX
      use m_Attributes,only : KT
      use m_Attributes,only : LAT
      use m_Attributes,only : LON
      use m_Attributes,only : LEV

      use m_Collector    ,only : Collector
      use m_Collector    ,only : globalSize
      use m_CollectorComm,only : allgatherv

      use m_die,only : die,perr
      implicit none

      type(Attributes),intent(in ) :: sattr	! send buffer
      type(Attributes),intent(out) :: rattr	! receive buffer
      type(Collector) ,intent(in)  :: coll	! a Collector
      integer,intent(in) :: comm		! communicator
      integer,optional,intent(out) :: stat	! return code

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgatherv__'
  integer :: kset
  integer :: gsize_
  integer :: ier
  integer,pointer,dimension(:) :: iattr
  real   ,pointer,dimension(:) :: xattr

	! This is the code for the subset Attributes

  if(present(stat)) stat=0

  kset=key(sattr)

	! Define a skeleton of Attributes

  gsize_=globalSize(coll)
  call Attributes_init(rattr,kset,gsize_)

	! Start message-passing

  ier=-1

  if(iand(kset,KR)/=0) then
	iattr => ptr_kr(rattr)
	call allgatherv(ptr_kr(sattr),iattr,coll,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allgatherv(kr)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif
  if(iand(kset,KT)/=0) then
	iattr => ptr_kt(rattr)
	call allgatherv(ptr_kt(sattr),iattr,coll,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allgatherv(kt)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif
  if(iand(kset,KX)/=0) then
	iattr => ptr_kx(rattr)
	call allgatherv(ptr_kx(sattr),iattr,coll,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allgatherv(kx)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif
  if(iand(kset,KS)/=0) then
	iattr => ptr_ks(rattr)
	call allgatherv(ptr_ks(sattr),iattr,coll,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allgatherv(ks)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif
  if(iand(kset,LAT)/=0) then
	xattr => ptr_lat(rattr)
	call allgatherv(ptr_lat(sattr),xattr,coll,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allgatherv(lat)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif
  if(iand(kset,LON)/=0) then
	xattr => ptr_lon(rattr)
	call allgatherv(ptr_lon(sattr),xattr,coll,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allgatherv(lon)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif
  if(iand(kset,LEV)/=0) then
	xattr => ptr_lev(rattr)
	call allgatherv(ptr_lev(sattr),xattr,coll,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allgatherv(lev)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif

  	if(ier/=0) then
	  call perr(myname_,'invalid subset',kset)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine allgatherv__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedtransposev__ - all-to-all send indexed data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexedtransposev__(indx,sattr,rattr,tran,comm,stat)
      use m_Attributes,only : Attributes
      use m_Attributes,only : Attributes_init
      use m_Attributes,only : key
      use m_Attributes,only : ptr_kr
      use m_Attributes,only : ptr_ks
      use m_Attributes,only : ptr_kx
      use m_Attributes,only : ptr_lat
      use m_Attributes,only : ptr_lon

      use m_Attributes,only : KR_SUBSET
      use m_Attributes,only : KS_SUBSET
      use m_Attributes,only : LL_SUBSET
      use m_Attributes,only : KX_SUBSET

      use m_Transposer    ,only : Transposer
      use m_Transposer    ,only : recvSize
      use m_TransposerComm,only : indexedTrans

      use m_die,only : die,perr
      implicit none

      integer,dimension(:),intent(in) :: indx	! sendSize(tran)
      type(Attributes),intent(in)  :: sattr	! sendSize(tran)
      type(Attributes),intent(out) :: rattr	! recvSize(tran)
      type(Transposer),intent(in)  :: tran	! see m_Transposer
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexedtransposev__'
  integer :: ier
  integer :: kset
  integer :: rsize_
  integer,pointer,dimension(:) :: iattr
  real   ,pointer,dimension(:) :: xattr

  if(present(stat)) stat=0

  kset=key(sattr)
  rsize_=recvSize(tran)
  call Attributes_init(rattr,kset,rsize_)

  ier=0

  select case(kset)
  case(KR_SUBSET,LL_SUBSET,KS_SUBSET,KX_SUBSET)

	iattr => ptr_kr(rattr)
	call indexedTrans(indx,ptr_kr(sattr),iattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'indexedTrans(kr)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)

    select case(kset)
    case(KS_SUBSET,KX_SUBSET)

	iattr => ptr_ks(rattr)
	call indexedTrans(indx,ptr_ks(sattr),iattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'indexedTrans(ks)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)

      select case(kset)
      case(KX_SUBSET)

	iattr => ptr_kx(rattr)
	call indexedTrans(indx,ptr_kx(sattr),iattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'indexedTrans(kx)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)

      end select

    case(LL_SUBSET)

	xattr => ptr_lat(rattr)
	call indexedTrans(indx,ptr_lat(sattr),xattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'indexedTrans(lat)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)

	xattr => ptr_lon(rattr)
	call indexedTrans(indx,ptr_lon(sattr),xattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'indexedTrans(lon)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)

    end select

  case default
    ier=-1
  end select

	if(ier/=0) then
	  call perr(myname_,'unknown subset',kset)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine indexedtransposev__
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexeduntransposev__ - all-to-all send indexed data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexeduntransposev__(sattr,rattr,tran,indx,comm,stat)
      use m_Attributes,only : Attributes
      use m_Attributes,only : Attributes_init
      use m_Attributes,only : key
      use m_Attributes,only : ptr_kr
      use m_Attributes,only : ptr_ks
      use m_Attributes,only : ptr_kx
      use m_Attributes,only : ptr_lat
      use m_Attributes,only : ptr_lon

      use m_Attributes,only : KR_SUBSET
      use m_Attributes,only : KS_SUBSET
      use m_Attributes,only : LL_SUBSET
      use m_Attributes,only : KX_SUBSET

      use m_Transposer    ,only : Transposer
      use m_Transposer    ,only : sendSize
      use m_TransposerComm,only : indexedUntrans

      use m_die,only : die,perr
      implicit none

      type(Attributes),intent(in)  :: sattr	! recvSize(tran)
      type(Attributes),intent(out) :: rattr	! sendSize(tran)
      type(Transposer),intent(in)  :: tran	! see m_Tranposer
      integer,dimension(:),intent(in) :: indx	! sendSize(tran)
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexeduntransposev__'
  integer :: ier
  integer :: kset
  integer :: ssize_
  integer,pointer,dimension(:) :: iattr
  real   ,pointer,dimension(:) :: xattr

  if(present(stat)) stat=0

			! Note that the size of send buffer is defined
			! by recvSize() and the size of recv buffer is
			! defined by sendSize().  Also, indx(:) now is
			! defined for the recv buffer.
  kset=key(sattr)
  ssize_=sendSize(tran)
  call Attributes_init(rattr,kset,ssize_)

  ier=0

  select case(kset)
  case(KR_SUBSET,LL_SUBSET,KS_SUBSET,KX_SUBSET)

	iattr => ptr_kr(rattr)
	call indexedUntrans(ptr_kr(sattr),iattr,tran,indx,	&
	  comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'indexedUntrans(kr)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)

    select case(kset)
    case(KS_SUBSET,KX_SUBSET)

	iattr => ptr_ks(rattr)
	call indexedUntrans(ptr_ks(sattr),iattr,tran,indx,	&
	  comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'indexedUntrans(ks)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)

      select case(kset)
      case(KX_SUBSET)

	iattr => ptr_kx(rattr)
	call indexedUntrans(ptr_kx(sattr),iattr,tran,indx,	&
	  comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'indexedUntrans(kx)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)

      end select

    case(LL_SUBSET)

	xattr => ptr_lat(rattr)
	call indexedUntrans(ptr_lat(sattr),xattr,tran,indx,	&
	  comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'indexedUntrans(lat)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)

	xattr => ptr_lon(rattr)
	call indexedUntrans(ptr_lon(sattr),xattr,tran,indx,	&
	  comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'indexedUntrans(lon)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)

    end select

  case default
    ier=-1
  end select

	if(ier/=0) then
	  call perr(myname_,'unknown subset',kset)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine indexeduntransposev__
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: transposev__ - all-to-all send data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine transposev__(sattr,rattr,tran,comm,stat)
      use m_Attributes,only : Attributes
      use m_Attributes,only : Attributes_init
      use m_Attributes,only : key
      use m_Attributes,only : ptr_kr
      use m_Attributes,only : ptr_ks
      use m_Attributes,only : ptr_kx
      use m_Attributes,only : ptr_kt
      use m_Attributes,only : ptr_lat
      use m_Attributes,only : ptr_lon
      use m_Attributes,only : ptr_lev

      use m_Attributes,only : KR
      use m_Attributes,only : KS
      use m_Attributes,only : KX
      use m_Attributes,only : KT
      use m_Attributes,only : LAT
      use m_Attributes,only : LON
      use m_Attributes,only : LEV

      use m_Transposer    ,only : Transposer
      use m_Transposer    ,only : recvSize
      use m_TransposerComm,only : trans

      use m_die,only : die,perr
      implicit none

      type(Attributes),intent(in)  :: sattr	! sendSize(tran)
      type(Attributes),intent(out) :: rattr	! recvSize(tran)
      type(Transposer),intent(in)  :: tran	! see m_Transposer
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::transposev__'
  integer :: ier
  integer :: kset
  integer :: rsize_
  integer,pointer,dimension(:) :: iattr
  real   ,pointer,dimension(:) :: xattr

  if(present(stat)) stat=0

  kset=key(sattr)
  rsize_=recvSize(tran)
  call Attributes_init(rattr,kset,rsize_)

  ier=-1

  if(iand(kset,KR)/=0) then
	iattr => ptr_kr(rattr)
	call trans(ptr_kr(sattr),iattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'trans(kr)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif

  if(iand(kset,KS)/=0) then
	iattr => ptr_ks(rattr)
	call trans(ptr_ks(sattr),iattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'trans(ks)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif

  if(iand(kset,KX)/=0) then
	iattr => ptr_kx(rattr)
	call trans(ptr_kx(sattr),iattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'trans(kx)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif

  if(iand(kset,KT)/=0) then
	iattr => ptr_kt(rattr)
	call trans(ptr_kt(sattr),iattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'trans(kt)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif

  if(iand(kset,LAT)/=0) then
	xattr => ptr_lat(rattr)
	call trans(ptr_lat(sattr),xattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'trans(lat)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif

  if(iand(kset,LON)/=0) then
	xattr => ptr_lon(rattr)
	call trans(ptr_lon(sattr),xattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'trans(lon)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif

  if(iand(kset,LEV)/=0) then
	xattr => ptr_lev(rattr)
	call trans(ptr_lev(sattr),xattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'trans(lev)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif

	if(ier/=0) then
	  call perr(myname_,'invalid subset',kset)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine transposev__
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: untransposev__ - all-to-all send  data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine untransposev__(sattr,rattr,tran,comm,stat)
      use m_Attributes,only : Attributes
      use m_Attributes,only : Attributes_init
      use m_Attributes,only : key
      use m_Attributes,only : ptr_kr
      use m_Attributes,only : ptr_ks
      use m_Attributes,only : ptr_kx
      use m_Attributes,only : ptr_kt
      use m_Attributes,only : ptr_lat
      use m_Attributes,only : ptr_lon
      use m_Attributes,only : ptr_lev

      use m_Attributes,only : KR
      use m_Attributes,only : KS
      use m_Attributes,only : KX
      use m_Attributes,only : KT
      use m_Attributes,only : LAT
      use m_Attributes,only : LON
      use m_Attributes,only : LEV

      use m_Transposer    ,only : Transposer
      use m_Transposer    ,only : sendSize
      use m_TransposerComm,only : untrans

      use m_die,only : die,perr
      implicit none

      type(Attributes),intent(in)  :: sattr	! recvSize(tran)
      type(Attributes),intent(out) :: rattr	! sendSize(tran)
      type(Transposer),intent(in)  :: tran	! see m_Tranposer
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::untransposev__'
  integer :: ier
  integer :: kset
  integer :: ssize_
  integer,pointer,dimension(:) :: iattr
  real   ,pointer,dimension(:) :: xattr

  if(present(stat)) stat=0

			! Note that the size of send buffer is defined
			! by recvSize() and the size of recv buffer is
			! defined by sendSize().
  kset=key(sattr)
  ssize_=sendSize(tran)
  call Attributes_init(rattr,kset,ssize_)

  ier=-1

  if(iand(kset,KR)/=0) then
	iattr => ptr_kr(rattr)
	call untrans(ptr_kr(sattr),iattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'untrans(kr)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif

  if(iand(kset,KS)/=0) then
	iattr => ptr_ks(rattr)
	call untrans(ptr_ks(sattr),iattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'untrans(ks)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif

  if(iand(kset,KX)/=0) then
	iattr => ptr_kx(rattr)
	call untrans(ptr_kx(sattr),iattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'untrans(kx)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif

  if(iand(kset,KT)/=0) then
	iattr => ptr_kt(rattr)
	call untrans(ptr_kt(sattr),iattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'untrans(kt)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(iattr)
  endif

  if(iand(kset,LAT)/=0) then
	xattr => ptr_lat(rattr)
	call untrans(ptr_lat(sattr),xattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'untrans(lat)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif

  if(iand(kset,LON)/=0) then
	xattr => ptr_lon(rattr)
	call untrans(ptr_lon(sattr),xattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'untrans(lon)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif

  if(iand(kset,LEV)/=0) then
	xattr => ptr_lev(rattr)
	call untrans(ptr_lev(sattr),xattr,tran,comm,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'untrans(lev)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
	nullify(xattr)
  endif

	if(ier/=0) then
	  call perr(myname_,'invalid subset',kset)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine untransposev__
end module m_AttributesComm
