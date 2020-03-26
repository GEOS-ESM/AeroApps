!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_PartitionComm - Partition-triple (indx,tran,coll)
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_PartitionComm
      implicit none
      private	! except

      public :: allgather
      public :: a2aForward
      public :: a2aBackward
      
      interface allgather; module procedure	&
	allgatherav_,	&
	allgatheri2_,	&
	allgatherr2_,	&
	allgatheriv_,	&
	allgatherrv_
      end interface

      interface a2aForward; module procedure	&
	a2aFwdav_,	&
	a2aFwdi2_,	&
	a2aFwdr2_,	&
	a2aFwdiv_,	&
	a2aFwdrv_
      end interface

      interface a2aBackward; module procedure	&
	a2aBwdav_,	&
	a2aBwdi2_,	&
	a2aBwdr2_,	&
	a2aBwdiv_,	&
	a2aBwdrv_
      end interface

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_PartitionComm'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgatherav_ - Allgather local vectors for a share vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgatherav_(local,shared, part,comm,stat)
      use m_die,only : perr,die
      use m_AttrVect,only : AttrVect
      use m_CollectorComm,only : indexedallgatherv
      use m_Partition,only : Partition
      use m_Partition,only : ptr_coll
      use m_Partition,only : ptr_indx

      implicit none
      type(AttrVect)    ,intent(in)  :: local  ! to be gathered
      type(AttrVect)    ,intent(out) :: shared ! the gathered
      type(Partition)   ,intent(in)  :: part
      integer           ,intent(in)  :: comm   ! Communicator
      integer,optional  ,intent(out) :: stat

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgatherav_'
  integer :: ier
  
  if(present(stat)) stat=0

  call indexedallgatherv((ptr_indx(part)),local,shared,	&
	(ptr_coll(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedallgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine allgatherav_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgatheri2_ - Allgather local vectors for a share vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgatheri2_(local,shared, part,comm,stat)
      use m_die,only : perr,die
      use m_CollectorComm,only : indexedallgatherv
      use m_Partition,only : Partition
      use m_Partition,only : ptr_coll
      use m_Partition,only : ptr_indx

      implicit none
      integer,dimension(:,:),intent(in)  :: local  ! :,localSize(part)
      integer,dimension(:,:),intent(out) :: shared ! :,globalSize(part)
      type(Partition)   ,intent(in)  :: part
      integer           ,intent(in)  :: comm   ! Communicator
      integer,optional  ,intent(out) :: stat

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgatheri2_'
  integer :: ier
  
  if(present(stat)) stat=0

  call indexedallgatherv((ptr_indx(part)),local,shared,	&
	(ptr_coll(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedallgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine allgatheri2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgatherr2_ - Allgather local vectors for a share vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgatherr2_(local,shared, part,comm,stat)
      use m_die,only : perr,die
      use m_CollectorComm,only : indexedallgatherv
      use m_Partition,only : Partition
      use m_Partition,only : ptr_coll
      use m_Partition,only : ptr_indx

      implicit none
      real,dimension(:,:),intent(in)  :: local  ! :,localSize(part)
      real,dimension(:,:),intent(out) :: shared ! :,globalSize(part)
      type(Partition)   ,intent(in)  :: part
      integer           ,intent(in)  :: comm   ! Communicator
      integer,optional  ,intent(out) :: stat

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgatherr2_'
  integer :: ier
  
  if(present(stat)) stat=0

  call indexedallgatherv((ptr_indx(part)),local,shared,	&
	(ptr_coll(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedallgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine allgatherr2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgatheriv_ - Allgather local vectors for a share vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgatheriv_(local,shared, part,comm,stat)
      use m_die,only : perr,die
      use m_CollectorComm,only : indexedallgatherv
      use m_Partition,only : Partition
      use m_Partition,only : ptr_coll
      use m_Partition,only : ptr_indx

      implicit none
      integer,dimension(:),intent(in)  :: local  ! localSize(part)
      integer,dimension(:),intent(out) :: shared ! globalSize(part)
      type(Partition)   ,intent(in)  :: part
      integer           ,intent(in)  :: comm   ! Communicator
      integer,optional  ,intent(out) :: stat

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgatheriv_'
  integer :: ier
  
  if(present(stat)) stat=0

  call indexedallgatherv((ptr_indx(part)),local,shared,	&
	(ptr_coll(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedallgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine allgatheriv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgatherrv_ - Allgather local vectors for a share vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgatherrv_(local,shared, part,comm,stat)
      use m_die,only : perr,die
      use m_CollectorComm,only : indexedallgatherv
      use m_Partition,only : Partition
      use m_Partition,only : ptr_coll
      use m_Partition,only : ptr_indx

      implicit none
      real,dimension(:),intent(in)  :: local  ! localSize(part)
      real,dimension(:),intent(out) :: shared ! globalSize(part)
      type(Partition)   ,intent(in)  :: part
      integer           ,intent(in)  :: comm   ! Communicator
      integer,optional  ,intent(out) :: stat

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgatherrv_'
  integer :: ier
  
  if(present(stat)) stat=0

  call indexedallgatherv((ptr_indx(part)),local,shared,	&
	(ptr_coll(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedallgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine allgatherrv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2afwdav_ - transpose forward
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2afwdav_(local,distr,part,comm,stat)
      use m_die,only : perr,die
      use m_AttrVect,only : AttrVect
      use m_TransposerComm,only : indexedTrans
      use m_Partition,only : Partition
      use m_Partition,only : ptr_tran
      use m_Partition,only : ptr_indx
      implicit none
      type(AttrVect) ,intent(in) :: local ! :,localSize(part)
      type(AttrVect) ,intent(out):: distr ! :,distrSize(part)
      type(Partition),intent(in) :: part
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2afwdav_'
  integer :: ier

  if(present(stat)) stat=0

  call indexedTrans((ptr_indx(part)),local,distr,	&
	(ptr_tran(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedTrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine a2afwdav_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2afwdi2_ - transpose forward
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2afwdi2_(local,distr,part,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : indexedTrans
      use m_Partition,only : Partition
      use m_Partition,only : ptr_tran
      use m_Partition,only : ptr_indx
      implicit none
      integer,dimension(:,:),intent(in)  :: local  ! :,localSize(part)
      integer,dimension(:,:),intent(out) :: distr  ! :,distrSize(part)
      type(Partition),intent(in) :: part
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2afwdi2_'
  integer :: ier

  if(present(stat)) stat=0

  call indexedTrans((ptr_indx(part)),local,distr,	&
	(ptr_tran(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedTrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine a2afwdi2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2afwdr2_ - transpose forward
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2afwdr2_(local,distr,part,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : indexedTrans
      use m_Partition,only : Partition
      use m_Partition,only : ptr_tran
      use m_Partition,only : ptr_indx
      implicit none
      real,dimension(:,:),intent(in)  :: local	! :,localSize(part)
      real,dimension(:,:),intent(out) :: distr	! :,distrSize(part)
      type(Partition),intent(in) :: part
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2afwdr2_'
  integer :: ier

  if(present(stat)) stat=0

  call indexedTrans((ptr_indx(part)),local,distr,	&
	(ptr_tran(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedTrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine a2afwdr2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2afwdiv_ - transpose forward
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2afwdiv_(local,distr,part,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : indexedTrans
      use m_Partition,only : Partition
      use m_Partition,only : ptr_tran
      use m_Partition,only : ptr_indx
      implicit none
      integer,dimension(:),intent(in)  :: local  ! localSize(part)
      integer,dimension(:),intent(out) :: distr  ! distrSize(part)
      type(Partition),intent(in) :: part
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2afwdiv_'
  integer :: ier

  if(present(stat)) stat=0

  call indexedTrans((ptr_indx(part)),local,distr,	&
	(ptr_tran(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedTrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine a2afwdiv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2afwdrv_ - transpose forward
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2afwdrv_(local,distr,part,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : indexedTrans
      use m_Partition,only : Partition
      use m_Partition,only : ptr_tran
      use m_Partition,only : ptr_indx
      implicit none
      real,dimension(:),intent(in)  :: local	! localSize(part)
      real,dimension(:),intent(out) :: distr	! distrSize(part)
      type(Partition),intent(in) :: part
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2afwdrv_'
  integer :: ier

  if(present(stat)) stat=0

  call indexedTrans((ptr_indx(part)),local,distr,	&
	(ptr_tran(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedTrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine a2afwdrv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2abwdav_ - transpose backward
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2abwdav_(distr,local,part,comm,stat)
      use m_die,only : perr,die
      use m_AttrVect,only : AttrVect
      use m_TransposerComm,only : indexedUntrans
      use m_Partition,only : Partition
      use m_Partition,only : ptr_tran
      use m_Partition,only : ptr_indx
      implicit none
      type(AttrVect) ,intent(in) :: distr ! :,distrSize(part)
      type(AttrVect) ,intent(out):: local ! :,localSize(part)
      type(Partition),intent(in) :: part
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2abwdav_'
  integer :: ier

  if(present(stat)) stat=0

  call indexedUntrans(distr,local,(ptr_tran(part)),	&
	(ptr_indx(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedUntrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine a2abwdav_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2abwdi2_ - transpose backward
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2abwdi2_(distr,local,part,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : indexedUntrans
      use m_Partition,only : Partition
      use m_Partition,only : ptr_tran
      use m_Partition,only : ptr_indx
      implicit none
      integer,dimension(:,:),intent(in) :: distr ! :,distrSize(part)
      integer,dimension(:,:),intent(out):: local ! :,localSize(part)
      type(Partition),intent(in) :: part
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2abwdi2_'
  integer :: ier

  if(present(stat)) stat=0

  call indexedUntrans(distr,local,(ptr_tran(part)),	&
	(ptr_indx(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedUntrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine a2abwdi2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2abwdr2_ - transpose backward
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2abwdr2_(distr,local,part,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : indexedUntrans
      use m_Partition,only : Partition
      use m_Partition,only : ptr_tran
      use m_Partition,only : ptr_indx
      implicit none
      real,dimension(:,:),intent(in) :: distr ! :,distrSize(part)
      real,dimension(:,:),intent(out):: local ! :,localSize(part)
      type(Partition),intent(in) :: part
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2abwdr2_'
  integer :: ier

  if(present(stat)) stat=0

  call indexedUntrans(distr,local,(ptr_tran(part)),	&
	(ptr_indx(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedUntrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine a2abwdr2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2abwdiv_ - transpose backward
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2abwdiv_(distr,local,part,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : indexedUntrans
      use m_Partition,only : Partition
      use m_Partition,only : ptr_tran
      use m_Partition,only : ptr_indx
      implicit none
      integer,dimension(:),intent(in) :: distr ! distrSize(part)
      integer,dimension(:),intent(out):: local ! localSize(part)
      type(Partition),intent(in) :: part
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2abwdiv_'
  integer :: ier

  if(present(stat)) stat=0

  call indexedUntrans(distr,local,(ptr_tran(part)),	&
	(ptr_indx(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedUntrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine a2abwdiv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2abwdrv_ - transpose backward
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2abwdrv_(distr,local,part,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : indexedUntrans
      use m_Partition,only : Partition
      use m_Partition,only : ptr_tran
      use m_Partition,only : ptr_indx
      implicit none
      real,dimension(:),intent(in) :: distr ! distrSize(part)
      real,dimension(:),intent(out):: local ! localSize(part)
      type(Partition),intent(in) :: part
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2abwdrv_'
  integer :: ier

  if(present(stat)) stat=0

  call indexedUntrans(distr,local,(ptr_tran(part)),	&
	(ptr_indx(part)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'indexedUntrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine a2abwdrv_

end module m_PartitionComm
