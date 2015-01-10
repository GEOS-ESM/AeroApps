!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_CollectorComm - Communications of a Collector
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_CollectorComm
      implicit none
      private	! except

      public :: gatherv
      public :: scatterv
      public :: allgatherv
      public :: indexedAllGatherv

      interface gatherv; module procedure	&
	gathervi2_,	&
	gathervr2_,	&
	gathervi_,	&
	gathervr_
      end interface

      interface scatterv; module procedure	&
	scattervi2_,	&
	scattervr2_,	&
	scattervi_,	&
	scattervr_
      end interface

      interface allgatherv; module procedure	&
	allgatherv_,	&
	allgathervi2_,	&
	allgathervr2_,	&
	allgathervi_,	&
	allgathervr_
      end interface

      interface indexedAllGatherv; module procedure	&
	indexedallgatherv_,	&
	indexedallgathervi2_,	&
	indexedallgathervr2_,	&
	indexedallgathervi_,	&
	indexedallgathervr_
      end interface

! !REVISION HISTORY:
! 	30Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_CollectorComm'

#include "assert.H"

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgatherv_ - all-gatherv of attrvect
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgatherv_(attr_send,attr_recv,coll,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : AttrVect_init => init
      use m_AttrVect,only : nIAttr
      use m_AttrVect,only : nRAttr
      use m_AttrVect,only : ptr_iAttr,ptr_rAttr
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : globalSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      type(AttrVect),  intent(in)  :: attr_send
      type(AttrVect),  intent(out) :: attr_recv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgatherv_'

  integer,pointer,dimension(:,:) :: iAttrSend,iAttrRecv
  real   ,pointer,dimension(:,:) :: rAttrSend,rAttrRecv
  integer,pointer,dimension(:)   :: counts,displs

  integer :: nIA,nRA
  integer :: lsize,gsize
  integer :: iAttrType,rAttrType
  integer :: ier

  if(present(stat)) stat=0

  nIA=nIAttr(attr_send)
  nRA=nRAttr(attr_send)
  lsize=localSize(coll)
  gsize=globalSize(coll)

  call AttrVect_init(attr_recv,attr_send,gsize)

  iAttrSend=>ptr_iAttr(attr_send)
  rAttrSend=>ptr_rAttr(attr_send)
  iAttrRecv=>ptr_iAttr(attr_recv)
  rAttrRecv=>ptr_rAttr(attr_recv)

  iAttrType=MP_type(iAttrSend)
  rAttrType=MP_type(rAttrSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_allgatherv(					&
	iAttrSend,nIA*lsize,		iAttrType,	&
	iAttrRecv,nIA*counts,nIA*displs,iAttrType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv(iAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  call MPI_allgatherv(					&
	rAttrSend,nRA*lsize,		rAttrType,	&
	rAttrRecv,nRA*counts,nRA*displs,rAttrType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv(rAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(iAttrSend)
  nullify(rAttrSend)
  nullify(iAttrRecv)
  nullify(rAttrRecv)

  nullify(counts)
  nullify(displs)

contains
subroutine write_()
end subroutine write_
end subroutine allgatherv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgathervi_ - all_gatherv() of integer(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgathervi_(iSend,iRecv,coll,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      integer,dimension(:),intent(in)  :: iSend
      integer,dimension(:),intent(out) :: iRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgathervi_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  lsize=localSize(coll)
  mesgType=MP_type(iSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_allgatherv(			&
	iSend,lsize,	    mesgType,	&
	iRecv,counts,displs,mesgType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine allgathervi_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgathervr_ - all_gatherv() of real(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgathervr_(rSend,rRecv,coll,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      real   ,dimension(:),intent(in)  :: rSend
      real   ,dimension(:),intent(out) :: rRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgathervr_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  lsize=localSize(coll)
  mesgType=MP_type(rSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_allgatherv(			&
	rSend,lsize,	    mesgType,	&
	rRecv,counts,displs,mesgType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine allgathervr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgathervi2_ - all_gatherv() of integer(:,:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgathervi2_(iSend,iRecv,coll,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      integer,dimension(:,:),intent(in)  :: iSend
      integer,dimension(:,:),intent(out) :: iRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgathervi2_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize,ldm
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  ldm  =size(iSend,1)
  lsize=localSize(coll)
  mesgType=MP_type(iSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_allgatherv(				&
	iSend,ldm*lsize,	    mesgType,	&
	iRecv,ldm*counts,ldm*displs,mesgType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine allgathervi2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgathervr2_ - all_gatherv() of real(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgathervr2_(rSend,rRecv,coll,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      real   ,dimension(:,:),intent(in)  :: rSend
      real   ,dimension(:,:),intent(out) :: rRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgathervr2_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize,ldm
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  ldm  =size(rSend,1)
  lsize=localSize(coll)
  mesgType=MP_type(rSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_allgatherv(				&
	rSend,ldm*lsize,	    mesgType,	&
	rRecv,ldm*counts,ldm*displs,mesgType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine allgathervr2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedallgatherv_ - all-gatherv of attrvect
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexedallgatherv_(indx,attr_send,	&
	attr_recv,coll,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : AttrVect_init => init
      use m_AttrVect,only : nIAttr
      use m_AttrVect,only : nRAttr
      use m_AttrVect,only : ptr_iAttr,ptr_rAttr
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : globalSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      integer,dimension(:),intent(in) :: indx	! pre-comm sorting
      type(AttrVect),  intent(in)  :: attr_send
      type(AttrVect),  intent(out) :: attr_recv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexedallgatherv_'

  integer,pointer,dimension(:,:) :: iAttrSend,iAttrRecv
  real   ,pointer,dimension(:,:) :: rAttrSend,rAttrRecv
  integer,pointer,dimension(:)   :: counts,displs

  integer :: nIA,nRA
  integer :: lsize,gsize
  integer :: iAttrType,rAttrType
  integer :: ier

  if(present(stat)) stat=0

  nIA=nIAttr(attr_send)
  nRA=nRAttr(attr_send)
  lsize=localSize(coll)
  gsize=globalSize(coll)

  call AttrVect_init(attr_recv,attr_send,gsize)

  iAttrSend=>ptr_iAttr(attr_send)
  rAttrSend=>ptr_rAttr(attr_send)
  iAttrRecv=>ptr_iAttr(attr_recv)
  rAttrRecv=>ptr_rAttr(attr_recv)

  iAttrType=MP_type(iAttrSend)
  rAttrType=MP_type(rAttrSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_allgatherv(					&
	iAttrSend(:,indx(1:lsize)),			&
	          nIA*lsize,		iAttrType,	&
	iAttrRecv,nIA*counts,nIA*displs,iAttrType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv(iAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  call MPI_allgatherv(					&
	rAttrSend(:,indx(1:lsize)),			&
	          nRA*lsize,		rAttrType,	&
	rAttrRecv,nRA*counts,nRA*displs,rAttrType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv(rAttr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(iAttrSend)
  nullify(rAttrSend)
  nullify(iAttrRecv)
  nullify(rAttrRecv)

  nullify(counts)
  nullify(displs)

end subroutine indexedallgatherv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedallgathervi_ - all_gatherv() of integer(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexedallgathervi_(indx,iSend,iRecv,coll,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      integer,dimension(:),intent(in)  :: indx	! pre-comm sorting
      integer,dimension(:),intent(in)  :: iSend
      integer,dimension(:),intent(out) :: iRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexedallgathervi_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  lsize=localSize(coll)
  mesgType=MP_type(iSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_allgatherv(			&
	iSend(indx(1:lsize)),		&
	      lsize,	    mesgType,	&
	iRecv,counts,displs,mesgType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine indexedallgathervi_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedallgathervr_ - all_gatherv() of real(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexedallgathervr_(indx,rSend,rRecv,coll,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      integer,dimension(:),intent(in)  :: indx	! pre-comm sorting
      real   ,dimension(:),intent(in)  :: rSend
      real   ,dimension(:),intent(out) :: rRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexedallgathervr_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  lsize=localSize(coll)
  mesgType=MP_type(rSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_allgatherv(			&
	rSend(indx(1:lsize)),		&
	      lsize,	    mesgType,	&
	rRecv,counts,displs,mesgType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine indexedallgathervr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedallgathervi2_ - all_gatherv() of integer(:,:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexedallgathervi2_(indx,iSend,iRecv,coll,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      integer,dimension(:)  ,intent(in)  :: indx  ! pre-comm sorting
      integer,dimension(:,:),intent(in)  :: iSend
      integer,dimension(:,:),intent(out) :: iRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexedallgathervi2_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize,ldm
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  ldm  =size(iSend,1)
  lsize=localSize(coll)
  mesgType=MP_type(iSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_allgatherv(				&
	iSend(:,indx(1:lsize)),			&
	      ldm*lsize,	    mesgType,	&
	iRecv,ldm*counts,ldm*displs,mesgType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine indexedallgathervi2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexedallgathervr2_ - all_gatherv() of real(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine indexedallgathervr2_(indx,rSend,rRecv,coll,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      integer,dimension(:)  ,intent(in)  :: indx  ! pre-comm sorting
      real   ,dimension(:,:),intent(in)  :: rSend
      real   ,dimension(:,:),intent(out) :: rRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexedallgathervr2_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize,ldm
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  ldm  =size(rSend,1)
  lsize=localSize(coll)
  mesgType=MP_type(rSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_allgatherv(				&
	rSend(:,indx(1:lsize)),			&
	      ldm*lsize,	    mesgType,	&
	rRecv,ldm*counts,ldm*displs,mesgType,	&
	comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_allgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine indexedallgathervr2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gathervi_ - gatherv() of integer(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gathervi_(iSend,iRecv,coll,root,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      integer,dimension(:),intent(in)  :: iSend
      integer,dimension(:),intent(out) :: iRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: root
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gathervi_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  lsize=localSize(coll)

  mesgType=MP_type(iSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_gatherv(			&
	iSend,lsize,	    mesgType,	&
	iRecv,counts,displs,mesgType,	&
	root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_gatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine gathervi_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gathervr_ - gatherv() of real(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gathervr_(rSend,rRecv,coll,root,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      real   ,dimension(:),intent(in)  :: rSend
      real   ,dimension(:),intent(out) :: rRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: root
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gathervr_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  lsize=localSize(coll)
  mesgType=MP_type(rSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_gatherv(			&
	rSend,lsize,	    mesgType,	&
	rRecv,counts,displs,mesgType,	&
	root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_gatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine gathervr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gathervi2_ - gatherv() of integer(:,:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gathervi2_(iSend,iRecv,coll,root,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die,assert_
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      integer,dimension(:,:),intent(in)  :: iSend
      integer,dimension(:,:),intent(out) :: iRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: root
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gathervi2_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize,ldm
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  ldm  =size(iSend,1)

	ASSERT(size(iSend,1)==size(iRecv,1))

  lsize=localSize(coll)
  mesgType=MP_type(iSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_gatherv(				&
	iSend,ldm*lsize,	    mesgType,	&
	iRecv,ldm*counts,ldm*displs,mesgType,	&
	root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_gatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine gathervi2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gathervr2_ - gatherv() of real(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gathervr2_(rSend,rRecv,coll,root,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die,assert_
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      real   ,dimension(:,:),intent(in)  :: rSend
      real   ,dimension(:,:),intent(out) :: rRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: root
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gathervr2_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize,ldm
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  ldm  =size(rSend,1)

	ASSERT(size(rSend,1)==size(rRecv,1))

  lsize=localSize(coll)
  mesgType=MP_type(rSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_gatherv(				&
	rSend,ldm*lsize,	    mesgType,	&
	rRecv,ldm*counts,ldm*displs,mesgType,	&
	root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_gatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine gathervr2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scattervi_ - scatterv() of integer(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine scattervi_(iSend,iRecv,coll,root,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      integer,dimension(:),intent(in)  :: iSend
      integer,dimension(:),intent(out) :: iRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: root
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scattervi_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  lsize=localSize(coll)

  mesgType=MP_type(iSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_scatterv(			&
	iSend,counts,displs,mesgType,	&
	iRecv,lsize,	    mesgType,	&
	root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_scatterv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine scattervi_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scattervr_ - scatterv() of real(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine scattervr_(rSend,rRecv,coll,root,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      real   ,dimension(:),intent(in)  :: rSend
      real   ,dimension(:),intent(out) :: rRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: root
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scattervr_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: lsize
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  lsize=localSize(coll)
  mesgType=MP_type(rSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_scatterv(			&
	rSend,counts,displs,mesgType,	&
	rRecv,lsize,	    mesgType,	&
	root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_scatterv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine scattervr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scattervi2_ - scatterv() of integer(:,:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine scattervi2_(iSend,iRecv,coll,root,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die, assert_
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      integer,dimension(:,:),intent(in)  :: iSend
      integer,dimension(:,:),intent(out) :: iRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: root
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scattervi2_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: ldm
  integer :: lsize
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  ldm  =size(iSend,1)

	ASSERT(size(iSend,1)==size(iRecv,1))

  lsize=localSize(coll)
  mesgType=MP_type(iSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_scatterv(			&
	iSend,ldm*counts,ldm*displs,mesgType,	&
	iRecv,ldm*lsize,	    mesgType,	&
	root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_scatterv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine scattervi2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scattervr2_ - scatterv() of real(:,:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine scattervr2_(rSend,rRecv,coll,root,comm,stat)  
      use m_mpif90,only : MP_type,MP_perr
      use m_die,only : die,assert_
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : ptr_counts,ptr_displs
      implicit none

      real   ,dimension(:,:),intent(in)  :: rSend
      real   ,dimension(:,:),intent(out) :: rRecv
      type(Collector), intent(in)  :: coll
      integer,         intent(in)  :: root
      integer,         intent(in)  :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scattervr2_'

  integer,pointer,dimension(:)   :: counts,displs
  integer :: ldm
  integer :: lsize
  integer :: mesgType
  integer :: ier

  if(present(stat)) stat=0

  ldm  =size(rSend,1)

	ASSERT(size(rSend,1)==size(rRecv,1))

  lsize=localSize(coll)
  mesgType=MP_type(rSend)

  counts =>ptr_counts(coll)
  displs =>ptr_displs(coll)

  call MPI_scatterv(			&
	rSend,ldm*counts,ldm*displs,mesgType,	&
	rRecv,ldm*lsize,	    mesgType,	&
	root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_scatterv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nullify(counts)
  nullify(displs)

end subroutine scattervr2_
end module m_CollectorComm
