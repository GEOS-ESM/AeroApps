!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_sigDCWindErr - Fcst. err. cov. models of decoupled winds
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_sigDCWindErr
      implicit none
      private	! except

      public :: sigDCWindErr_store
      public :: sigDCWindErr_remove
      public :: WPOTE

      interface sigDCWindErr_store; module procedure	&
	mp_store_,	&
	store_
      end interface

      interface sigDCWindErr_remove; module procedure	&
	remove_
      end interface

! !REVISION HISTORY:
! 	08Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_sigDCWindErr'

  character(len=*),parameter :: WPOTE='WPOTE'
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: mp_store_ - store_() on multiple processors.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine mp_store_(comm,root,mtyp,nlev,levs,npar,tabl,desc,stat)

      use FEsigW_tabl,only : FEsigW_type
      use FEsigW_tabl,only : FEsigW_desc
      use FEsigW_tabl,only : FEsigW_nlev
      use FEsigW_tabl,only : FEsigW_plev
      use FEsigW_tabl,only : FEsigW_npar
      use FEsigW_tabl,only : FEsigW_pars
      use FEsigW_tabl,only : FEsigW_Mpar
      use FEsigW_tabl,only : FEsigW_Mlev
      use FEsigW_tabl,only : FEsigW_defined
      use FEsigW_tabl,only : FEsigW_setbyRC

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_type,MP_perr
      use m_stdio,only : stderr
      use m_die,  only : die
      implicit none

		! Commnication arguments required on all PEs

      integer,intent(in) :: comm
      integer,intent(in) :: root

		! Data arguments are only significant on the root PE

      character(len=*), intent(in) :: mtyp
      integer,          intent(in) :: nlev
      real,dimension(:),intent(in) :: levs
      integer,          intent(in) :: npar
      real,dimension(:,:),intent(in) :: tabl
      character(len=*),optional,intent(in) :: desc

      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	08Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::mp_store_'

  integer,dimension(2) :: nsize
  integer :: myID
  integer :: ier

  if(present(stat)) stat=0

	! Exchange the size information

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif

  if(myID==root) then
    nsize(1)=nlev
    nsize(2)=npar
  endif

  call MPI_Bcast(nsize,size(nsize),MP_type(nsize),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(nsize)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif

  FEsigW_nlev=nsize(1)
  FEsigW_npar=nsize(2)

	! Verify the arguments, only on the root

if(myID==root) then
  if(	nlev>size(tabl,2) .or.	&
	npar>size(tabl,1)	) then

    if(	nlev>size(tabl,2) ) write(stderr,'(2a,2(a,i5))')	&
	myname_,': invalid arguments',		&
	', nlev ='        ,nlev,		&
	', size(tabl,2) =',size(tabl,2)

    if(	npar>size(tabl,1) ) write(stderr,'(2a,2(a,i5))')	&
	myname_,': invalid arguments',		&
	', npar ='        ,npar,		&
	', size(tabl,1) =',size(tabl,1)

    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif
endif

	! Verify the storage, on all PEs.

  if(	FEsigW_nlev>FEsigW_Mlev .or.	&
	FEsigW_npar>FEsigW_Mpar	) then

    if(	FEsigW_nlev>FEsigW_Mlev ) write(stderr,'(2a,2(a,i5))')	&
	myname_,': invalid arguments',				&
	', nlev ='       ,FEsigW_nlev,				&
	', FEsigW_Mlev =',FEsigW_Mlev

    if(	FEsigW_npar>FEsigW_Mpar ) write(stderr,'(2a,2(a,i5))')	&
	myname_,': invalid arguments',				&
	', npar ='       ,FEsigW_npar,				&
	', FEsigW_Mpar =',FEsigW_Mpar

    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

if(myID==root) then

	! Save the information on the root PE.

  FEsigW_type=mtyp
  FEsigW_plev(1:nlev)=levs(1:nlev)
  FEsigW_pars(1:npar,1:nlev)=tabl(1:npar,1:nlev)

  FEsigW_desc=""
  if(present(desc)) FEsigW_desc=desc
endif

	! Broadcast the information to all PEs.

  call MPI_bcast(FEsigW_type,len(FEsigW_type),		&
    MP_type(FEsigW_type),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(type)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif

  call MPI_bcast(FEsigW_plev(              1:FEsigW_nlev),	&
    FEsigW_nlev            ,MP_type(FEsigW_plev),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(plev)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif

  call MPI_bcast(FEsigW_pars(1:FEsigW_npar,1:FEsigW_nlev),	&
    FEsigW_nlev*FEsigW_npar,MP_type(FEsigW_pars),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(pars)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif

  call MPI_bcast(FEsigW_desc,len(FEsigW_desc),		&
     MP_type(FEsigW_desc),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(desc)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif

  FEsigW_defined=.true.
  FEsigW_setbyRC=.false.

end subroutine mp_store_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: store_ - store the mass-decoupled wind er. parameter table
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine store_(mtyp,nlev,levs,npar,tabl,desc,stat)

      use FEsigW_tabl,only : FEsigW_type
      use FEsigW_tabl,only : FEsigW_desc
      use FEsigW_tabl,only : FEsigW_nlev
      use FEsigW_tabl,only : FEsigW_plev
      use FEsigW_tabl,only : FEsigW_npar
      use FEsigW_tabl,only : FEsigW_pars
      use FEsigW_tabl,only : FEsigW_Mpar
      use FEsigW_tabl,only : FEsigW_Mlev
      use FEsigW_tabl,only : FEsigW_defined
      use FEsigW_tabl,only : FEsigW_setbyRC

      use m_stdio,only : stderr
      use m_die,  only : die
      implicit none

      character(len=*), intent(in) :: mtyp
      integer,          intent(in) :: nlev
      real,dimension(:),intent(in) :: levs
      integer,          intent(in) :: npar
      real,dimension(:,:),intent(in) :: tabl
      character(len=*),optional,intent(in) :: desc

      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	08Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::store_'

  if(present(stat)) stat=0

	! Verify the arguments

  if(	nlev>size(tabl,2) .or.	&
	npar>size(tabl,1)	) then

    if(	nlev>size(tabl,2) ) write(stderr,'(2a,2(a,i5))')	&
	myname_,': invalid arguments',		&
	', nlev ='        ,nlev,		&
	', size(tabl,2) =',size(tabl,2)

    if(	npar>size(tabl,1) ) write(stderr,'(2a,2(a,i5))')	&
	myname_,': invalid arguments',		&
	', npar ='        ,npar,		&
	', size(tabl,1) =',size(tabl,1)

    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

	! Verify the storage

  if(	nlev>FEsigW_Mlev .or.	&
	npar>FEsigW_Mpar	) then

    if(	nlev>FEsigW_Mlev ) write(stderr,'(2a,2(a,i5))')	&
	myname_,': invalid arguments',		&
	', nlev ='       ,nlev,			&
	', FEsigW_Mlev =',FEsigW_Mlev

    if(	npar>FEsigW_Mpar ) write(stderr,'(2a,2(a,i5))')	&
	myname_,': invalid arguments',		&
	', npar ='       ,npar,			&
	', FEsigW_Mpar =',FEsigW_Mpar

    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

	! Save the information

  FEsigW_type=mtyp
  FEsigW_nlev=nlev
  FEsigW_plev(1:nlev)=levs(1:nlev)

  FEsigW_npar=npar
  FEsigW_pars(1:npar,1:nlev)=tabl(1:npar,1:nlev)

  FEsigW_desc=""
  if(present(desc)) FEsigW_desc=desc

  FEsigW_defined=.true.
  FEsigW_setbyRC=.false.

end subroutine store_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: remove_ - remove the component
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine remove_(stat)
      use FEsigW_tabl,only : FEsigW_type
      use FEsigW_tabl,only : FEsigW_desc
      use FEsigW_tabl,only : FEsigW_nlev
      use FEsigW_tabl,only : FEsigW_npar
      use FEsigW_tabl,only : FEsigW_defined
      use FEsigW_tabl,only : FEsigW_setbyRC
      implicit none
      integer,optional,intent(out):: stat

! !REVISION HISTORY:
! 	08Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::remove_'

  if(present(stat)) stat=0

  FEsigW_defined=.false.
  FEsigW_setbyRC=.false.

  FEsigW_nlev=0
  FEsigW_npar=0
  FEsigW_type=""
  FEsigW_desc=""

end subroutine remove_

end module m_sigDCWindErr
