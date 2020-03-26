!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_TransposerComm - Communication through a Transposer
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_TransposerComm
      implicit none
      private	! except

      public :: trans
      public :: indexedTrans
      public :: untrans
      public :: indexedUntrans

      interface trans; module procedure	&
	transposev_,	&
	transposevr2_,	&
	transposevi2_,	&
	transposevr_,	&
	transposevi_
      end interface

      interface untrans; module procedure	&
	untransposev_,		&
	untransposevr2_,	&
	untransposevi2_,	&
	untransposevr_,		&
	untransposevi_
      end interface

      interface indexedTrans; module procedure	&
	indexedtransposev_,	&
	indexedtransposevr2_,	&
	indexedtransposevi2_,	&
	indexedtransposevr_,	&
	indexedtransposevi_
      end interface

      interface indexedUntrans; module procedure	&
	indexeduntransposev_,	&
	indexeduntransposevr2_,	&
	indexeduntransposevi2_,	&
	indexeduntransposevr_,	&
	indexeduntransposevi_
      end interface

! !REVISION HISTORY:
! 	30Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_TransposerComm'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: transposev_ - transpose a distrbuted AttrVect
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine transposev_(attr_send,attr_recv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : AttrVect_init
      use m_AttrVect,only : ptr_iAttr,ptr_rAttr
      use m_AttrVect,only : nIAttr
      use m_AttrVect,only : nRAttr
      use m_Transposer,only : Transposer
      use m_Transposer,only : recvSize
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      type(AttrVect)  ,intent(in)  :: attr_send
      type(AttrVect)  ,intent(out) :: attr_recv
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::transposev_'

  integer :: nIA,nRA
  integer :: rsize
  integer,pointer,dimension(:,:) :: iAttrSend
  real   ,pointer,dimension(:,:) :: rAttrSend
  integer,pointer,dimension(:,:) :: iAttrRecv
  real   ,pointer,dimension(:,:) :: rAttrRecv

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls

  integer :: iAttrType,rAttrType
  integer :: ier

  if(present(stat)) stat=0

  rsize=recvSize(trans)

  scounts => ptr_sendCounts(trans)
  sdispls => ptr_sendDispls(trans)
  rcounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)

  call AttrVect_init(attr_recv,attr_send,rsize)

  nIA=nIAttr(attr_send)
  nRA=nRAttr(attr_send)

  iAttrSend =>ptr_iAttr(attr_send)
  rAttrSend =>ptr_rAttr(attr_send)
  iAttrRecv =>ptr_iAttr(attr_recv)
  rAttrRecv =>ptr_rAttr(attr_recv)

  iAttrType=MP_type(iAttrSend)
  rAttrType=MP_type(rAttrSend)

	! Send INTEGER data first

  call mpi_alltoallv(						&
	iAttrSend,nIA*scounts(:),nIA*sdispls(:),iAttrType,	&
	iAttrRecv,nIA*rcounts(:),nIA*rdispls(:),iAttrType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv(iAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	! Send REAL data second

  call mpi_alltoallv(						&
	rAttrSend,nRA*scounts(:),nRA*sdispls(:),rAttrType,	&
	rAttrRecv,nRA*rcounts(:),nRA*rdispls(:),rAttrType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv(rAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  
  nullify(iAttrSend)
  nullify(rAttrSend)
  nullify(iAttrRecv)
  nullify(rAttrRecv)

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine transposev_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: untransposev_ - untranspose a distrbuted AttrVect
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine untransposev_(attr_send,attr_recv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : AttrVect_init
      use m_AttrVect,only : ptr_iAttr,ptr_rAttr
      use m_AttrVect,only : nIAttr
      use m_AttrVect,only : nRAttr
      use m_Transposer,only : Transposer
      use m_Transposer,only : sendSize
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      type(AttrVect)  ,intent(in)  :: attr_send
      type(AttrVect)  ,intent(out) :: attr_recv
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::untransposev_'

  integer :: nIA,nRA
  integer :: ssize
  integer,pointer,dimension(:,:) :: iAttrSend
  real   ,pointer,dimension(:,:) :: rAttrSend
  integer,pointer,dimension(:,:) :: iAttrRecv
  real   ,pointer,dimension(:,:) :: rAttrRecv

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls

  integer :: iAttrType,rAttrType
  integer :: ier

  if(present(stat)) stat=0

  ssize=sendSize(trans)

  scounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)
  rcounts => ptr_sendCounts(trans)
  rdispls => ptr_recvDispls(trans)

  call AttrVect_init(attr_recv,attr_send,ssize)

  nIA=nIAttr(attr_send)
  nRA=nRAttr(attr_send)

  iAttrSend =>ptr_iAttr(attr_send)
  rAttrSend =>ptr_rAttr(attr_send)
  iAttrRecv =>ptr_iAttr(attr_recv)
  rAttrRecv =>ptr_rAttr(attr_recv)

  iAttrType=MP_type(iAttrSend)
  rAttrType=MP_type(rAttrSend)

	! Send INTEGER data first

  call mpi_alltoallv(						&
	iAttrSend,nIA*scounts(:),nIA*sdispls(:),iAttrType,	&
	iAttrRecv,nIA*rcounts(:),nIA*rdispls(:),iAttrType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv(iAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	! Send REAL data second

  call mpi_alltoallv(						&
	rAttrSend,nRA*scounts(:),nRA*sdispls(:),rAttrType,	&
	rAttrRecv,nRA*rcounts(:),nRA*rdispls(:),rAttrType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv(rAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  
  nullify(iAttrSend)
  nullify(rAttrSend)
  nullify(iAttrRecv)
  nullify(rAttrRecv)

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine untransposev_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: transposevr2_ - transpose a distrbuted REAL vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine transposevr2_(rSend,rRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      real  ,dimension(:,:),intent(in)  :: rSend
      real  ,dimension(:,:),intent(out) :: rRecv
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::transposevr2_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ldm
  integer :: ier

  if(present(stat)) stat=0

  scounts => ptr_sendCounts(trans)
  sdispls => ptr_sendDispls(trans)
  rcounts => ptr_recvCounts(trans)
  rdispls => ptr_recvDispls(trans)

  ldm=size(rSend,1)
  mesgType=MP_type(rSend)

  call mpi_alltoallv(				&
	rSend,ldm*scounts(:),ldm*sdispls(:),mesgType,	&
	rRecv,ldm*rcounts(:),ldm*rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine transposevr2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: transposevi2_ - transpose a distrbuted INTEGER vecter
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine transposevi2_(iSend,iRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      integer,dimension(:,:),intent(in)  :: iSend
      integer,dimension(:,:),intent(out) :: iRecv
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::transposevi2_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ldm
  integer :: ier

  if(present(stat)) stat=0

  scounts => ptr_sendCounts(trans)
  sdispls => ptr_sendDispls(trans)
  rcounts => ptr_recvCounts(trans)
  rdispls => ptr_recvDispls(trans)

  ldm=size(iSend,1)
  mesgType=MP_type(iSend)

  call mpi_alltoallv(				&
	iSend,ldm*scounts(:),ldm*sdispls(:),mesgType,	&
	iRecv,ldm*rcounts(:),ldm*rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine transposevi2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: transposevr_ - transpose a distrbuted REAL vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine transposevr_(rSend,rRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      real  ,dimension(:),intent(in)  :: rSend
      real  ,dimension(:),intent(out) :: rRecv
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::transposevr_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  scounts => ptr_sendCounts(trans)
  sdispls => ptr_sendDispls(trans)
  rcounts => ptr_recvCounts(trans)
  rdispls => ptr_recvDispls(trans)

  mesgType=MP_type(rSend)

  call mpi_alltoallv(				&
	rSend,scounts(:),sdispls(:),mesgType,	&
	rRecv,rcounts(:),rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine transposevr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: transposevi_ - transpose a distrbuted INTEGER vecter
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine transposevi_(iSend,iRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      integer,dimension(:),intent(in)  :: iSend
      integer,dimension(:),intent(out) :: iRecv
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::transposevi_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  scounts => ptr_sendCounts(trans)
  sdispls => ptr_sendDispls(trans)
  rcounts => ptr_recvCounts(trans)
  rdispls => ptr_recvDispls(trans)

  mesgType=MP_type(iSend)

  call mpi_alltoallv(				&
	iSend,scounts(:),sdispls(:),mesgType,	&
	iRecv,rcounts(:),rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine transposevi_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: untransposevr2_ - untranspose a distrbuted REAL vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine untransposevr2_(rSend,rRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      real  ,dimension(:,:),intent(in)  :: rSend
      real  ,dimension(:,:),intent(out) :: rRecv
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::untransposevr2_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ldm
  integer :: ier

  if(present(stat)) stat=0

  scounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)
  rcounts => ptr_sendCounts(trans)
  rdispls => ptr_sendDispls(trans)

  ldm=size(rSend,1)
  mesgType=MP_type(rSend)

  call mpi_alltoallv(					&
	rSend,ldm*scounts(:),ldm*sdispls(:),mesgType,	&
	rRecv,ldm*rcounts(:),ldm*rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine untransposevr2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: untransposevi2_ - untranspose a INTEGER vecter
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine untransposevi2_(iSend,iRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      integer,dimension(:,:),intent(in)  :: iSend
      integer,dimension(:,:),intent(out) :: iRecv
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::untransposevi2_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ldm
  integer :: ier

  if(present(stat)) stat=0

  scounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)
  rcounts => ptr_sendCounts(trans)
  rdispls => ptr_sendDispls(trans)

  ldm=size(iSend,1)
  mesgType=MP_type(iSend)

  call mpi_alltoallv(					&
	iSend,ldm*scounts(:),ldm*sdispls(:),mesgType,	&
	iRecv,ldm*rcounts(:),ldm*rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine untransposevi2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: untransposevr_ - untranspose a distrbuted REAL vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine untransposevr_(rSend,rRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      real  ,dimension(:),intent(in)  :: rSend
      real  ,dimension(:),intent(out) :: rRecv
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::untransposevr_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  scounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)
  rcounts => ptr_sendCounts(trans)
  rdispls => ptr_sendDispls(trans)

  mesgType=MP_type(rSend)

  call mpi_alltoallv(				&
	rSend,scounts(:),sdispls(:),mesgType,	&
	rRecv,rcounts(:),rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine untransposevr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: untransposevi_ - untranspose a INTEGER vecter
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine untransposevi_(iSend,iRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      integer,dimension(:),intent(in)  :: iSend
      integer,dimension(:),intent(out) :: iRecv
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::untransposevi_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  scounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)
  rcounts => ptr_sendCounts(trans)
  rdispls => ptr_sendDispls(trans)

  mesgType=MP_type(iSend)

  call mpi_alltoallv(				&
	iSend,scounts(:),sdispls(:),mesgType,	&
	iRecv,rcounts(:),rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine untransposevi_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedtransposev_ - transpose a distrbuted AttrVect
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexedtransposev_(indx,attr_send,attr_recv,	&
	trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : AttrVect_init
      use m_AttrVect,only : ptr_iAttr,ptr_rAttr
      use m_AttrVect,only : nIAttr
      use m_AttrVect,only : nRAttr
      use m_Transposer,only : Transposer
      use m_Transposer,only : sendSize
      use m_Transposer,only : recvSize
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      integer,dimension(:),intent(in) :: indx	! sendSize(trans)
      type(AttrVect)  ,intent(in)  :: attr_send	! :,sendSize(trans)
      type(AttrVect)  ,intent(out) :: attr_recv	! :,recvSize(trans)
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexedtransposev_'

  integer :: nIA,nRA
  integer :: ssize,rsize
  integer,pointer,dimension(:,:) :: iAttrSend
  real   ,pointer,dimension(:,:) :: rAttrSend
  integer,pointer,dimension(:,:) :: iAttrRecv
  real   ,pointer,dimension(:,:) :: rAttrRecv

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls

  integer :: iAttrType,rAttrType
  integer :: ier

  if(present(stat)) stat=0

  ssize=sendSize(trans)
  rsize=recvSize(trans)

  scounts => ptr_sendCounts(trans)
  sdispls => ptr_sendDispls(trans)
  rcounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)

  call AttrVect_init(attr_recv,attr_send,rsize)

  nIA=nIAttr(attr_send)
  nRA=nRAttr(attr_send)

  iAttrSend =>ptr_iAttr(attr_send)
  rAttrSend =>ptr_rAttr(attr_send)
  iAttrRecv =>ptr_iAttr(attr_recv)
  rAttrRecv =>ptr_rAttr(attr_recv)

  iAttrType=MP_type(iAttrSend)
  rAttrType=MP_type(rAttrSend)

	! Send INTEGER data first

  call mpi_alltoallv(						&
	iAttrSend(:,indx(1:ssize)),				&
		  nIA*scounts(:),nIA*sdispls(:),iAttrType,	&
	iAttrRecv,nIA*rcounts(:),nIA*rdispls(:),iAttrType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv(iAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	! Send REAL data second

  call mpi_alltoallv(						&
	rAttrSend(:,indx(1:ssize)),				&
		  nRA*scounts(:),nRA*sdispls(:),rAttrType,	&
	rAttrRecv,nRA*rcounts(:),nRA*rdispls(:),rAttrType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv(rAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  
  nullify(iAttrSend)
  nullify(rAttrSend)
  nullify(iAttrRecv)
  nullify(rAttrRecv)

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine indexedtransposev_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexeduntransposev_ - untranspose a distrbuted AttrVect
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexeduntransposev_(attr_send,attr_recv,	&
	trans,indx,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : AttrVect_init
      use m_AttrVect,only : ptr_iAttr,ptr_rAttr
      use m_AttrVect,only : nIAttr
      use m_AttrVect,only : nRAttr
      use m_Transposer,only : Transposer
      use m_Transposer,only : sendSize
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      type(AttrVect)  ,intent(in)  :: attr_send	! :,recvSize(trans)
      type(AttrVect)  ,intent(out) :: attr_recv	! :,sendSize(trans)
      type(Transposer),intent(in)  :: trans
      integer,dimension(:),intent(in) :: indx	! sendSize(trans)
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter::myname_=myname//'::indexeduntransposev_'

  integer :: nIA,nRA
  integer :: rsize
  integer,pointer,dimension(:,:) :: iAttrSend
  real   ,pointer,dimension(:,:) :: rAttrSend
  integer,pointer,dimension(:,:) :: iAttrRecv
  real   ,pointer,dimension(:,:) :: rAttrRecv

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls

  integer :: iAttrType,rAttrType
  integer :: ier

  if(present(stat)) stat=0

  rsize=sendSize(trans)

  scounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)
  rcounts => ptr_sendCounts(trans)
  rdispls => ptr_recvDispls(trans)

  call AttrVect_init(attr_recv,attr_send,rsize)

  nIA=nIAttr(attr_send)
  nRA=nRAttr(attr_send)

  iAttrSend =>ptr_iAttr(attr_send)
  rAttrSend =>ptr_rAttr(attr_send)
  iAttrRecv =>ptr_iAttr(attr_recv)
  rAttrRecv =>ptr_rAttr(attr_recv)

  iAttrType=MP_type(iAttrSend)
  rAttrType=MP_type(rAttrSend)

	! Send INTEGER data first

  call mpi_alltoallv(						&
	iAttrSend,nIA*scounts(:),nIA*sdispls(:),iAttrType,	&
	iAttrRecv,nIA*rcounts(:),nIA*rdispls(:),iAttrType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv(iAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  iAttrRecv(:,indx(1:rsize))=iAttrRecv(:,1:rsize)

	! Send REAL data second

  call mpi_alltoallv(						&
	rAttrSend,nRA*scounts(:),nRA*sdispls(:),rAttrType,	&
	rAttrRecv,nRA*rcounts(:),nRA*rdispls(:),rAttrType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv(rAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  rAttrRecv(:,indx(1:rsize))=rAttrRecv(:,1:rsize)
  
  nullify(iAttrSend)
  nullify(rAttrSend)
  nullify(iAttrRecv)
  nullify(rAttrRecv)

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine indexeduntransposev_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedtransposevr2_ - transpose a distrbuted REAL vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexedtransposevr2_(indx,rSend,rRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : sendSize
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      integer,dimension(:) ,intent(in)  :: indx	  ! sendSize(trans)
      real  ,dimension(:,:),intent(in)  :: rSend  ! :,sendSize(trans)
      real  ,dimension(:,:),intent(out) :: rRecv  ! :,recvSize(trans)
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter::myname_=myname//'::indexedtransposevr2_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ldm
  integer :: ier
  integer :: ssize

  if(present(stat)) stat=0

  scounts => ptr_sendCounts(trans)
  sdispls => ptr_sendDispls(trans)
  rcounts => ptr_recvCounts(trans)
  rdispls => ptr_recvDispls(trans)

  ldm=size(rSend,1)
  ssize=sendSize(trans)
  mesgType=MP_type(rSend)

  call mpi_alltoallv(					&
	rSend(:,indx(1:ssize)),				&
	      ldm*scounts(:),ldm*sdispls(:),mesgType,	&
	rRecv,ldm*rcounts(:),ldm*rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine indexedtransposevr2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedtransposevi2_ - transpose a distrbuted INTEGER vecter
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexedtransposevi2_(indx,iSend,iRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      use m_Transposer,only : sendSize
      implicit none

      integer,dimension(:)  ,intent(in)  :: indx    ! SendSize(trans)
      integer,dimension(:,:),intent(in)  :: iSend   ! :,sendSize(trans)
      integer,dimension(:,:),intent(out) :: iRecv   ! :,recvSize(trans)
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter::myname_=myname//'::indexedtransposevi2_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ldm
  integer :: ier
  integer :: ssize

  if(present(stat)) stat=0

  scounts => ptr_sendCounts(trans)
  sdispls => ptr_sendDispls(trans)
  rcounts => ptr_recvCounts(trans)
  rdispls => ptr_recvDispls(trans)

  ldm=size(iSend,1)
  ssize=sendSize(trans)
  mesgType=MP_type(iSend)

  call mpi_alltoallv(					&
	iSend(:,indx(1:ssize)),				&
	      ldm*scounts(:),ldm*sdispls(:),mesgType,	&
	iRecv,ldm*rcounts(:),ldm*rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine indexedtransposevi2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedtransposevr_ - transpose a distrbuted REAL vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexedtransposevr_(indx,rSend,rRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : sendSize
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      integer,dimension(:),intent(in) :: indx	! sendSize(trans)
      real  ,dimension(:),intent(in)  :: rSend	! sendSize(trans)
      real  ,dimension(:),intent(out) :: rRecv	! recvSize(trans)
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter::myname_=myname//'::indexedtransposevr_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ier
  integer :: ssize

  if(present(stat)) stat=0

  scounts => ptr_sendCounts(trans)
  sdispls => ptr_sendDispls(trans)
  rcounts => ptr_recvCounts(trans)
  rdispls => ptr_recvDispls(trans)

  ssize=sendSize(trans)
  mesgType=MP_type(rSend)

  call mpi_alltoallv(				&
	rSend(indx(1:ssize)),			&
	      scounts(:),sdispls(:),mesgType,	&
	rRecv,rcounts(:),rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine indexedtransposevr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedtransposevi_ - transpose a distrbuted INTEGER vecter
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexedtransposevi_(indx,iSend,iRecv,trans,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : sendSize
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      integer,dimension(:),intent(in)  :: indx	! sendSize(trans)
      integer,dimension(:),intent(in)  :: iSend	! sendSize(trans)
      integer,dimension(:),intent(out) :: iRecv	! recvSize(trans)
      type(Transposer),intent(in)  :: trans
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter::myname_=myname//'::indexedtransposevi_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ier
  integer :: ssize

  if(present(stat)) stat=0

  scounts => ptr_sendCounts(trans)
  sdispls => ptr_sendDispls(trans)
  rcounts => ptr_recvCounts(trans)
  rdispls => ptr_recvDispls(trans)

  ssize=sendSize(trans)
  mesgType=MP_type(iSend)

  call mpi_alltoallv(				&
	iSend(indx(1:ssize)),			&
	      scounts(:),sdispls(:),mesgType,	&
	iRecv,rcounts(:),rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine indexedtransposevi_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexeduntransposevr2_ - untranspose a distrbuted REAL vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexeduntransposevr2_(rSend,rRecv,trans,indx,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : sendSize
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      real  ,dimension(:,:),intent(in)  :: rSend ! :,recvSize(trans)
      real  ,dimension(:,:),intent(out) :: rRecv ! :,sendSize(trans)
      type(Transposer),intent(in)  :: trans
      integer,dimension(:),intent(in) :: indx	 ! sendSize(trans)
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter::myname_=myname//'::indexeduntransposevr2_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ldm
  integer :: ier
  integer :: rsize

  if(present(stat)) stat=0

  scounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)
  rcounts => ptr_sendCounts(trans)
  rdispls => ptr_sendDispls(trans)

  rsize=sendSize(trans)
  ldm=size(rSend,1)
  mesgType=MP_type(rSend)

  call mpi_alltoallv(					&
	rSend,ldm*scounts(:),ldm*sdispls(:),mesgType,	&
	rRecv,ldm*rcounts(:),ldm*rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  rRecv(:,indx(1:rsize))=rRecv(:,1:rsize)

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine indexeduntransposevr2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexeduntransposevi2_ - untranspose a INTEGER vecter
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexeduntransposevi2_(iSend,iRecv,trans,indx,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : sendSize
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      integer,dimension(:,:),intent(in)  :: iSend ! :,recvSize(trans)
      integer,dimension(:,:),intent(out) :: iRecv ! :,sendSize(trans)
      type(Transposer),intent(in)  :: trans
      integer,dimension(:),intent(in) :: indx	  ! sendSize(trans)
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter::myname_=myname//'::indexeduntransposevi2_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ldm
  integer :: ier
  integer :: rsize

  if(present(stat)) stat=0

  scounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)
  rcounts => ptr_sendCounts(trans)
  rdispls => ptr_sendDispls(trans)

  rsize=sendSize(trans)
  ldm=size(iSend,1)
  mesgType=MP_type(iSend)

  call mpi_alltoallv(					&
	iSend,ldm*scounts(:),ldm*sdispls(:),mesgType,	&
	iRecv,ldm*rcounts(:),ldm*rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  iRecv(:,indx(1:rsize))=iRecv(:,1:rsize)

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine indexeduntransposevi2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexeduntransposevr_ - untranspose a distrbuted REAL vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexeduntransposevr_(rSend,rRecv,trans,indx,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : sendSize
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      real  ,dimension(:),intent(in)  :: rSend	! recvSize(trans)
      real  ,dimension(:),intent(out) :: rRecv	! sendSize(trans)
      type(Transposer),intent(in)  :: trans
      integer,dimension(:),intent(in) :: indx	! sendSize(trans)
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter::myname_=myname//'::indexeduntransposevr_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ier
  integer :: rsize

  if(present(stat)) stat=0

  scounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)
  rcounts => ptr_sendCounts(trans)
  rdispls => ptr_sendDispls(trans)

  rsize=sendSize(trans)
  mesgType=MP_type(rSend)

  call mpi_alltoallv(				&
	rSend,scounts(:),sdispls(:),mesgType,	&
	rRecv,rcounts(:),rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  rRecv(indx(1:rsize))=rRecv(1:rsize)

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine indexeduntransposevr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexeduntransposevi_ - untranspose a INTEGER vecter
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexeduntransposevi_(iSend,iRecv,trans,indx,comm,stat)
      use m_mpif90,only : MP_type,MP_perr
      use m_die,   only : die
      use m_Transposer,only : Transposer
      use m_Transposer,only : sendSize
      use m_Transposer,only : ptr_sendCounts,ptr_sendDispls
      use m_Transposer,only : ptr_recvCounts,ptr_recvDispls
      implicit none

      integer,dimension(:),intent(in)  :: iSend	! recvSize(trans)
      integer,dimension(:),intent(out) :: iRecv	! sendSize(trans)
      type(Transposer),intent(in)  :: trans
      integer,dimension(:),intent(in) :: indx	! sendSize(trans)
      integer         ,intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	22Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter::myname_=myname//'::indexeduntransposevi_'

  integer,pointer,dimension(:) :: scounts,sdispls
  integer,pointer,dimension(:) :: rcounts,rdispls
  integer :: mesgType
  integer :: ier
  integer :: rsize

  if(present(stat)) stat=0

  scounts => ptr_recvCounts(trans)
  sdispls => ptr_recvDispls(trans)
  rcounts => ptr_sendCounts(trans)
  rdispls => ptr_sendDispls(trans)

  rsize=sendSize(trans)
  mesgType=MP_type(iSend)

  call mpi_alltoallv(				&
	iSend,scounts(:),sdispls(:),mesgType,	&
	iRecv,rcounts(:),rdispls(:),mesgType,	&
	comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'mpi_alltoallv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  iRecv(indx(1:rsize))=iRecv(1:rsize)

  nullify(scounts)
  nullify(sdispls)
  nullify(rcounts)
  nullify(rdispls)

end subroutine indexeduntransposevi_

end module m_TransposerComm
