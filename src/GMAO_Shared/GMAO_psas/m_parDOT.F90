!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_parDOT - parallel dot-product
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_parDOT
      implicit none
      private	! except

      public :: parDOT		! The class data structure
      public :: parNRM2

      public :: parL1_norm	! l1-norm
      public :: parL2_norm	! l2-norm
      public :: parL3_norm	! l3-norm
      public :: parLi_norm	! l-infinity-norm

      interface parDOT ; module procedure	&
	vdotm_,		&
	vdot_ ; end interface
      interface parNRM2; module procedure	&
	vnrm2m_,	&
	vnrm2_; end interface

      interface parL1_norm; module procedure	&
	vnrm1m_,	&
	vnrm1_; end interface
      interface parL2_norm; module procedure	&
	vnrm2m_,	&
	vnrm2_; end interface
      interface parL3_norm; module procedure	&
	vnrm3m_,	&
	vnrm3_; end interface
      interface parLi_norm; module procedure	&
	vnrmim_,	&
	vnrmi_; end interface

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_parDOT'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vdot_ - parallel dot-product
!
! !DESCRIPTION:
!
! !INTERFACE:

    function vdot_(x,y,comm)
      use m_mpif90,only : MP_type
      use m_mpif90,only : MP_SUM
      use m_die   ,only : MP_die
      implicit none
      real,dimension(:),intent(in) :: x
      real,dimension(:),intent(in) :: y
      integer,intent(in) :: comm
      real :: vdot_

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vdot_'
  integer :: ier
  real :: sdot,rdot

  sdot=dot_product(x,y)

  call MPI_allreduce(sdot,rdot,1,MP_TYPE(sdot),MP_SUM,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_allreduce()',ier)

  vdot_=rdot

end function vdot_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vdotm_ - parallel dot-product for multiple vectors
!
! !DESCRIPTION:
!
! !INTERFACE:

    function vdotm_(x,y,comm)
      use m_mpif90,only : MP_type
      use m_mpif90,only : MP_SUM
      use m_die   ,only : MP_die,die
      use m_mall  ,only : mall_ison,mall_mci,mall_mco
      implicit none
      real,dimension(:,:),intent(in) :: x	! (n,nvecs)
      real,dimension(:,:),intent(in) :: y	! (n,nvecs)
      integer,intent(in) :: comm
      real,dimension(size(x,2)) :: vdotm_

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vdotm_'
  integer :: i,ier
  integer :: nvecs
  real,allocatable,dimension(:) :: sdot,rdot

  nvecs=size(x,2)
  if(nvecs==0) return

	allocate(sdot(nvecs),rdot(nvecs),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(sdot,myname)
		  call mall_mci(rdot,myname)
		endif

  do i=1,nvecs
    sdot(i)=dot_product(x(:,i),y(:,i))
  end do

  call MPI_allreduce(sdot,rdot,nvecs,MP_TYPE(sdot),MP_SUM,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_allreduce()',ier)

  vdotm_(:)=rdot(:)

		if(mall_ison()) then
		  call mall_mco(sdot,myname)
		  call mall_mco(rdot,myname)
		endif
	deallocate(sdot,rdot,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end function vdotm_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vnrm1_ - parallel l1-norm with a single vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    function vnrm1_(x,comm)
      use m_mpif90,only : MP_type
      use m_mpif90,only : MP_SUM
      use m_die   ,only : MP_die
      implicit none
      real,dimension(:),intent(in) :: x
      integer,intent(in) :: comm
      real :: vnrm1_

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vnrm1_'
  real :: sumx
  integer :: ier

  sumx=0.
  if(size(x)>0) sumx=sum(abs(x))

  call MPI_allreduce((sumx),sumx,1,MP_type(sumx),MP_SUM,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_allreduce(sumx)',ier)

  vnrm1_=sumx
end function vnrm1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vnrm1m_ - parallel l1-norm with multiple vectors
!
! !DESCRIPTION:
!
! !INTERFACE:

    function vnrm1m_(x,comm)
      use m_mpif90,only : MP_type
      use m_mpif90,only : MP_SUM
      use m_die   ,only : MP_die,die
      use m_mall  ,only : mall_ison,mall_mci,mall_mco
      implicit none
      real,dimension(:,:),intent(in) :: x	! (n,nvecs)
      integer,intent(in) :: comm
      real,dimension(size(x,2)) :: vnrm1m_

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vnrm1m_'
  integer :: nvecs
  integer :: i,ier
  real,allocatable,dimension(:) :: sumx

  nvecs=size(x,2)
  if(nvecs==0) return

	allocate(sumx(nvecs),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(sumx,myname)

  sumx(:)=0.
  if(size(x)>0) then
    do i=1,nvecs
      sumx(i)=sum(abs(x(:,i)))
    end do
  endif

  call MPI_allreduce((sumx),sumx,nvecs,MP_type(sumx),MP_SUM,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_allreduce(sumx)',ier)

  vnrm1m_(:)=sumx(:)

		if(mall_ison()) call mall_mco(sumx,myname)
	deallocate(sumx,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end function vnrm1m_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vnrm2_ - parallel l2-norm with a single vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    function vnrm2_(x,comm)
      implicit none
      real,dimension(:),intent(in) :: x
      integer,intent(in) :: comm
      real :: vnrm2_

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vnrm2_'

  vnrm2_=sqrt(vdot_(x,x,comm))

end function vnrm2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vnrm2m_ - parallel l2-norm with multiple vectors
!
! !DESCRIPTION:
!
! !INTERFACE:

    function vnrm2m_(x,comm)
      implicit none
      real,dimension(:,:),intent(in) :: x	! (n,nvecs)
      integer,intent(in) :: comm
      real,dimension(size(x,2)) :: vnrm2m_

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vnrm2m_'

  vnrm2m_(:)=sqrt(vdotm_(x,x,comm))

end function vnrm2m_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vnrm3_ - parallel l1-norm with a single vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    function vnrm3_(x,comm)
      use m_mpif90,only : MP_type
      use m_mpif90,only : MP_SUM
      use m_die   ,only : MP_die
      implicit none
      real,dimension(:),intent(in) :: x
      integer,intent(in) :: comm
      real :: vnrm3_

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vnrm3_'
  real :: sumx
  integer :: ier

  sumx=0.
  if(size(x)>0) sumx=sum(abs(x))

  call MPI_allreduce((sumx),sumx,1,MP_type(sumx),MP_SUM,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_allreduce(sumx)',ier)

  vnrm3_=sumx
end function vnrm3_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vnrm3m_ - parallel l1-norm with multiple vectors
!
! !DESCRIPTION:
!
! !INTERFACE:

    function vnrm3m_(x,comm)
      use m_mpif90,only : MP_type
      use m_mpif90,only : MP_SUM
      use m_die   ,only : MP_die,die
      use m_mall  ,only : mall_ison,mall_mci,mall_mco
      implicit none
      real,dimension(:,:),intent(in) :: x	! (n,nvecs)
      integer,intent(in) :: comm
      real,dimension(size(x,2)) :: vnrm3m_

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vnrm3m_'
  integer :: nvecs
  integer :: i,ier
  real,allocatable,dimension(:) :: sumx

  nvecs=size(x,2)
  if(nvecs==0) return

	allocate(sumx(nvecs),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(sumx,myname)

  sumx(:)=0.
  if(size(x)>0) then
    do i=1,nvecs
      sumx(i)=sum(abs( x(:,i)*x(:,i)*x(:,i) ))
    end do
  endif

  call MPI_allreduce((sumx),sumx,nvecs,MP_type(sumx),MP_SUM,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_allreduce(sumx)',ier)

  vnrm3m_(:)=sumx(:)**(1./3.)

		if(mall_ison()) call mall_mco(sumx,myname)
	deallocate(sumx,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end function vnrm3m_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vnrmi_ - parallel l<inf>-norm with a single vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    function vnrmi_(x,comm,iPE,loc)
      use m_mpif90,only : MP_2type,MP_comm_rank,MP_MAXLOC
      use m_die   ,only : MP_die
      implicit none
      real,dimension(:),intent(in) :: x
      integer,intent(in) :: comm
      integer,optional,intent(out) :: iPE
      integer,optional,intent(out) :: loc
      real :: vnrmi_

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vnrmi_'
  real,dimension(2) :: sumx
  integer :: ier,myID

	call MP_comm_rank(comm,myID,ier)
		if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  sumx(1)=0.
  if(size(x)>0) sumx(1)=maxval(abs(x))
  sumx(2)=myID

  call MPI_allreduce((sumx),sumx,1,MP_2type(sumx(1)),MP_MAXLOC,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_allreduce(sumx)',ier)

  if(present(loc)) loc=maxloc(abs(x),dim=1)
  if(present(iPE)) iPE=nint(sumx(2))
  vnrmi_=sumx(1)
end function vnrmi_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: vnrmim_ - parallel l<inf>-norm with multiple vectors
!
! !DESCRIPTION:
!
! !INTERFACE:

    function vnrmim_(x,comm,iPE,loc)
      use m_mpif90,only : MP_2type,MP_comm_rank,MP_MAXLOC
      use m_die   ,only : MP_die,die
      use m_mall  ,only : mall_ison,mall_mci,mall_mco
      implicit none
      real,dimension(:,:),intent(in) :: x	! (n,nvecs)
      integer,intent(in) :: comm
      integer,dimension(size(x,2)),optional,intent(out) :: iPE
      integer,dimension(size(x,2)),optional,intent(out) :: loc
      real,dimension(size(x,2)) :: vnrmim_

! !REVISION HISTORY:
! 	23May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::vnrmim_'
  integer :: nvecs
  integer :: i,ier
  real,allocatable,dimension(:,:) :: sumx

  nvecs=size(x,2)
  if(nvecs==0) return

	allocate(sumx(2,nvecs),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(sumx,myname)

  sumx(1,:)=0.
  if(size(x)>0) then
    do i=1,nvecs
      sumx(1,i)=maxval(abs(x(:,i)))
    end do
  endif

  call MPI_allreduce((sumx),sumx,nvecs,MP_2type(sumx(1,1)),	&
    MP_MAXLOC,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_allreduce(sumx)',ier)

  if(present(loc)) loc(:)=maxloc(abs(x(:,:)),dim=1)
  if(present(iPE)) iPE(:)=nint(sumx(2,:))
  vnrmim_(:)=sumx(1,:)

		if(mall_ison()) call mall_mco(sumx,myname)
	deallocate(sumx,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end function vnrmim_

end module m_parDOT
