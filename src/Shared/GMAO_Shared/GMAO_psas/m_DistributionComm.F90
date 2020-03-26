!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_DistributionComm - Distribution-triple (indx,tran,coll)
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_DistributionComm
      implicit none
      private	! except

      public :: distribute
      public :: undistribute
      
      interface distribute; module procedure	&
	a2aFwdav_,	&
	a2aFwdi2_,	&
	a2aFwdr2_,	&
	a2aFwdiv_,	&
	a2aFwdrv_
      end interface

      interface undistribute; module procedure	&
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

  character(len=*),parameter :: myname='m_DistributionComm'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2afwdav_ - transpose forward
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2afwdav_(local,distr,dstr,comm,stat)
      use m_die,only : perr,die
      use m_Attributes,only : Attributes
      use m_Attributes,only : permute
      use m_Attributes,only : unpermute
      use m_Attributes,only : clean

      use m_AttributesComm,only : trans

      use m_Distribution,only : Distribution
      use m_Distribution,only : ptr_tran
      use m_Distribution,only : ptr_aindx
      use m_Distribution,only : ptr_prank
      implicit none
      type(Attributes)  ,intent(in) :: local
      type(Attributes)  ,intent(out):: distr
      type(Distribution),intent(in) :: dstr
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2afwdav_'
  integer :: ier
  type(Attributes) :: sattr	! send buffer
  type(Attributes) :: rattr	! recv buffer

  if(present(stat)) stat=0

  call permute(local,sattr,ptr_aindx(dstr))

  call trans(sattr,rattr,(ptr_tran(dstr)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'trans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  call unpermute(rattr,distr,ptr_prank(dstr))

	call clean(rattr)
	call clean(sattr)

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

    subroutine a2afwdi2_(local,distr,dstr,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : trans
      use m_Distribution,only : Distribution
      use m_Distribution,only : distrSize
      use m_Distribution,only : ptr_tran
      use m_Distribution,only : ptr_aindx
      use m_Distribution,only : ptr_prank
      implicit none
      integer,dimension(:,:),intent(in)  :: local  ! :,localSize(dstr)
      integer,dimension(:,:),intent(out) :: distr  ! :,distrSize(dstr)
      type(Distribution),intent(in) :: dstr
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2afwdi2_'
  integer :: ier
  integer,pointer,dimension(:) :: indx
  integer,pointer,dimension(:) :: rank
  integer :: dsize

  if(present(stat)) stat=0

  dsize=distrSize(dstr)

	indx => ptr_aindx(dstr)
	rank => ptr_prank(dstr)

  call trans(local(:,indx(:)),distr,ptr_tran(dstr),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'trans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  distr(:,rank(:)) = distr(:,1:dsize)

	nullify(rank)
	nullify(indx)

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

    subroutine a2afwdr2_(local,distr,dstr,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : trans
      use m_Distribution,only : Distribution
      use m_Distribution,only : distrSize
      use m_Distribution,only : ptr_tran
      use m_Distribution,only : ptr_aindx
      use m_Distribution,only : ptr_prank
      implicit none
      real,dimension(:,:),intent(in)  :: local	! :,localSize(dstr)
      real,dimension(:,:),intent(out) :: distr	! :,distrSize(dstr)
      type(Distribution),intent(in) :: dstr
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2afwdr2_'
  integer :: ier
  integer,pointer,dimension(:) :: indx
  integer,pointer,dimension(:) :: rank
  integer :: dsize

  if(present(stat)) stat=0

  dsize=distrSize(dstr)

	indx => ptr_aindx(dstr)
	rank => ptr_prank(dstr)

  call trans(local(:,indx(:)),distr,ptr_tran(dstr),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'trans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  distr(:,rank(:)) = distr(:,1:dsize)

	nullify(rank)
	nullify(indx)

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

    subroutine a2afwdiv_(local,distr,dstr,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : trans
      use m_Distribution,only : Distribution
      use m_Distribution,only : distrSize
      use m_Distribution,only : ptr_tran
      use m_Distribution,only : ptr_aindx
      use m_Distribution,only : ptr_prank
      implicit none
      integer,dimension(:),intent(in)  :: local  ! localSize(dstr)
      integer,dimension(:),intent(out) :: distr  ! distrSize(dstr)
      type(Distribution),intent(in) :: dstr
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2afwdiv_'
  integer :: ier
  integer,pointer,dimension(:) :: indx
  integer,pointer,dimension(:) :: rank
  integer :: dsize

  if(present(stat)) stat=0

  dsize=distrSize(dstr)

	indx => ptr_aindx(dstr)
	rank => ptr_prank(dstr)

  call trans(local(indx(:)),distr,ptr_tran(dstr),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'trans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  distr(rank(:)) = distr(1:dsize)

	nullify(rank)
	nullify(indx)

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

    subroutine a2afwdrv_(local,distr,dstr,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : trans
      use m_Distribution,only : Distribution
      use m_Distribution,only : distrSize
      use m_Distribution,only : ptr_tran
      use m_Distribution,only : ptr_aindx
      use m_Distribution,only : ptr_prank
      implicit none
      real,dimension(:),intent(in)  :: local	! localSize(dstr)
      real,dimension(:),intent(out) :: distr	! distrSize(dstr)
      type(Distribution),intent(in) :: dstr
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2afwdrv_'
  integer :: ier
  integer,pointer,dimension(:) :: indx
  integer,pointer,dimension(:) :: rank
  integer :: dsize

  if(present(stat)) stat=0

  dsize=distrSize(dstr)

	indx => ptr_aindx(dstr)
	rank => ptr_prank(dstr)

  call trans(local(indx(:)),distr,ptr_tran(dstr),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'trans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  distr(rank(:)) = distr(1:dsize)

	nullify(rank)
	nullify(indx)

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

    subroutine a2abwdav_(distr,local,dstr,comm,stat)
      use m_die,only : perr,die
      use m_Attributes,only : Attributes
      use m_Attributes,only : permute
      use m_Attributes,only : unpermute
      use m_Attributes,only : clean
      use m_AttributesComm,only : untrans

      use m_Distribution,only : Distribution
      use m_Distribution,only : ptr_tran
      use m_Distribution,only : ptr_aindx
      use m_Distribution,only : ptr_prank
      implicit none
      type(Attributes) ,intent(in) :: distr ! :,distrSize(dstr)
      type(Attributes) ,intent(out):: local ! :,localSize(dstr)
      type(Distribution),intent(in) :: dstr
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2abwdav_'
  integer :: ier
  type(Attributes) :: sattr	! send buffer
  type(Attributes) :: rattr	! recv buffer

  if(present(stat)) stat=0

  call permute(distr,sattr,ptr_prank(dstr))

  call untrans(sattr,rattr,(ptr_tran(dstr)),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'untrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  call unpermute(rattr,local,ptr_aindx(dstr))

	call clean(sattr)
	call clean(rattr)

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

    subroutine a2abwdi2_(distr,local,dstr,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : untrans
      use m_Distribution,only : Distribution
      use m_Distribution,only : localSize
      use m_Distribution,only : ptr_tran
      use m_Distribution,only : ptr_aindx
      use m_Distribution,only : ptr_prank
      implicit none
      integer,dimension(:,:),intent(in) :: distr ! :,distrSize(dstr)
      integer,dimension(:,:),intent(out):: local ! :,localSize(dstr)
      type(Distribution),intent(in) :: dstr
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2abwdi2_'
  integer :: ier
  integer,pointer,dimension(:) :: rank
  integer,pointer,dimension(:) :: indx
  integer :: lsize

  if(present(stat)) stat=0

  lsize=localSize(dstr)

	rank => ptr_prank(dstr)
	indx => ptr_aindx(dstr)

  call untrans(distr(:,rank(:)),local,ptr_tran(dstr),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'untrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  local(:,indx(:)) = local(:,1:lsize)

	nullify(indx)
	nullify(rank)

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

    subroutine a2abwdr2_(distr,local,dstr,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : untrans
      use m_Distribution,only : Distribution
      use m_Distribution,only : localSize
      use m_Distribution,only : ptr_tran
      use m_Distribution,only : ptr_aindx
      use m_Distribution,only : ptr_prank
      implicit none
      real,dimension(:,:),intent(in) :: distr ! :,distrSize(dstr)
      real,dimension(:,:),intent(out):: local ! :,localSize(dstr)
      type(Distribution),intent(in) :: dstr
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2abwdr2_'
  integer :: ier
  integer,pointer,dimension(:) :: rank
  integer,pointer,dimension(:) :: indx
  integer :: lsize

  if(present(stat)) stat=0

  lsize=localSize(dstr)

	rank => ptr_prank(dstr)
	indx => ptr_aindx(dstr)

  call untrans(distr(:,rank(:)),local,ptr_tran(dstr),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'untrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  local(:,indx(:))=local(:,1:lsize)

	nullify(indx)
	nullify(rank)

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

    subroutine a2abwdiv_(distr,local,dstr,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : untrans
      use m_Distribution,only : Distribution
      use m_Distribution,only : localSize
      use m_Distribution,only : ptr_tran
      use m_Distribution,only : ptr_aindx
      use m_Distribution,only : ptr_prank
      implicit none
      integer,dimension(:),intent(in) :: distr ! distrSize(dstr)
      integer,dimension(:),intent(out):: local ! localSize(dstr)
      type(Distribution),intent(in) :: dstr
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2abwdiv_'
  integer :: ier
  integer,pointer,dimension(:) :: rank
  integer,pointer,dimension(:) :: indx
  integer :: lsize

  if(present(stat)) stat=0

  lsize=localSize(dstr)

	rank => ptr_prank(dstr)
	indx => ptr_aindx(dstr)

  call untrans(distr(rank(:)),local,ptr_tran(dstr),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'untrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  local(indx(:))=local(1:lsize)

	nullify(indx)
	nullify(rank)

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

    subroutine a2abwdrv_(distr,local,dstr,comm,stat)
      use m_die,only : perr,die
      use m_TransposerComm,only : untrans
      use m_Distribution,only : Distribution
      use m_Distribution,only : localSize
      use m_Distribution,only : ptr_tran
      use m_Distribution,only : ptr_aindx
      use m_Distribution,only : ptr_prank
      implicit none
      real,dimension(:),intent(in) :: distr ! distrSize(dstr)
      real,dimension(:),intent(out):: local ! localSize(dstr)
      type(Distribution),intent(in) :: dstr
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2abwdrv_'
  integer :: ier
  integer,pointer,dimension(:) :: rank
  integer,pointer,dimension(:) :: indx
  integer :: lsize

  if(present(stat)) stat=0

  lsize=localSize(dstr)

	rank => ptr_prank(dstr)
	indx => ptr_aindx(dstr)

  call untrans(distr(rank(:)),local,ptr_tran(dstr),comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'untrans()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  local(indx(:))=local(1:lsize)

	nullify(indx)
	nullify(rank)

end subroutine a2abwdrv_

end module m_DistributionComm
