!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_xTab_levs - Tabulate a shared pressure level table
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_xTab_levs
      use rlev_imat
      implicit none
      private	! except

      public :: tab_levs		! The class data structure
      public :: rmtab_levs
      public :: index_levs
      public :: levs_xtab

      public :: MXveclev
      public :: nveclev
      public :: pveclev

      interface tab_levs; module procedure	&
	tab0_,	&
	tab1_
      end interface
      interface levs_xtab;  module procedure	&
	parxtab_,	&
	xtab_;   end interface
      interface rmtab_levs; module procedure rmtab_; end interface
      interface index_levs; module procedure	&
	index0_,	&
	index1_
      end interface

! !REVISION HISTORY:
! 	07Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_xTab_levs'

  integer,parameter :: NMAND = 30
  real,   dimension(NMAND),parameter :: PMAND=	&
	(/ 1000.,925.,850.,700.,600.,500.,400.,300.,250.,200.,	&
	    150.,100., 70., 50., 40., 30., 20., 10.,  7.,  5.,  &
	      4.,  3.,  2.,  1.,  .7,  .5,  .4,  .3,  .2,  .1	/)

  logical,save :: xTab_lev_defined=.false.

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: tab0_ - tabulate with a default table
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine tab0_(nlev)
      use config, only : pres4slp
      use m_stdio,only : stderr
      use m_die,  only : die
      implicit none
      integer,optional,intent(in) :: nlev ! by pass a compiler bug

! !REVISION HISTORY:
! 	07Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::tab0_'

  call merg_plevs(NMAND,PMAND, 1,pres4slp,	&
	MXveclev,nveclev,pveclev		)

end subroutine tab0_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: tab1_ - tabulate level entries
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine tab1_(ndat,rlev)
      use config, only : pres4slp
      implicit none

      integer,          intent(in) :: ndat
      real,dimension(:),intent(in) :: rlev

! !REVISION HISTORY:
! 	07Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::tab1_'

  call merg_plevs(NMAND,PMAND, 1,pres4slp,	&
	MXveclev,nveclev,pveclev		)
  call merg_plevs(ndat,rlev,nveclev,pveclev,	&
	MXveclev,nveclev,pveclev		)

end subroutine tab1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: xtab_ - extent the table
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine xtab_(ndat,rlev)
      use m_die,only : die
      implicit none
      integer,intent(in) :: ndat
      real,dimension(:),intent(in) :: rlev

! !REVISION HISTORY:
! 	28Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::xtab_'

  if(size(rlev)<ndat) call die(myname_,		&
	'size(rlev)',size(rlev),'ndat',ndat	)

  call merg_plevs(ndat,rlev,nveclev,pveclev,	&
	MXveclev,nveclev,pveclev		)

end subroutine xtab_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: parxtab_ - extent the table from all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine parxtab_(rlev,comm)
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_Collector,only : Collector
      use m_Collector,only : Collector_init
      use m_Collector,only : clean
      use m_Collector,only : globalSize
      use m_CollectorComm,only : allgatherv

      use m_mpout, only : mpout_log

      implicit none
      real,dimension(:),intent(in) :: rlev
      integer	       ,intent(in) :: comm

! !REVISION HISTORY:
! 	28Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::parxtab_'
  type(Collector) :: coll
  real,allocatable,dimension(:) :: alev, tlev
  integer :: ndat,nall
  integer :: ier

  allocate(tlev(size(rlev)),stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)

  Call unique(rlev, size(rlev), tlev, ndat)
  Call mpout_log(myname_,'ndat original: ',size(rlev))
  Call mpout_log(myname_,'ndat final:    ',ndat)

  call Collector_init(coll,ndat,comm)

  nall=globalSize(coll)
	allocate(alev(nall),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(alev,myname)

  call allgatherv(tlev,alev,coll,comm,stat=ier)
	if(ier/=0) call die(myname_,'allgatherv()',ier)
  
  call merg_plevs(nall,alev,nveclev,pveclev,	&
	MXveclev,nveclev,pveclev		)

		if(mall_ison()) call mall_mco(alev,myname)
		if(mall_ison()) call mall_mco(tlev,myname)
	deallocate(alev,tlev,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

  call clean(coll)

contains

  ! This routine selects the unique values of the array "x".
  ! First a sort is performed to make detection simple.
  Subroutine unique(x, nx, y, ny)
    Integer, Intent(In)  :: nx
    Real, Intent(In)     :: x(nx)
    Real, Intent(Out)    :: y(nx)
    Integer, Intent(Out) :: ny

    Integer :: index(nx), i
    Real :: y_next
    Real, Parameter :: sfix = 0.01 * 0.01

    call indexxr(nx,x,index)
    call permutr(x,index,nx,y)

    ny = 0
    y_next = -1
    If (nx > 0) ny = 1
    Do i = 1, nx - 1
       If (y(i+1) >= y_next) Then ! a new value
          ny = ny + 1
          y(ny) = y(i+1)
	  y_next = sfix * ((y(ny)/sfix) + 1)
       End If
    End Do

  End Subroutine unique

end subroutine parxtab_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rmtab_ - remove the defined table (a place holder)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine rmtab_()
      implicit none

! !REVISION HISTORY:
! 	07Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rmtab_'

  nveclev=0

end subroutine rmtab_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: index1_ - index w.r.t. pveclev
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine index1_(pres,ilev,wlev)
      use m_stdio,only : stderr
      use m_die,  only : die
      use m_mall, only : mall_ison,mall_mci,mall_mco
      implicit none
      real,   dimension(:),intent(in) :: pres
      integer,dimension(:),intent(out):: ilev
      real,   dimension(:),intent(out),optional :: wlev

! !REVISION HISTORY:
! 	15Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::index1_'
  logical,parameter :: NEAREST=.true.

  real,dimension(:),allocatable :: wtmp
  integer :: ier

  if(present(wlev)) then
    call slogtab(.not.NEAREST,nveclev,pveclev,	&
	size(pres),pres,ilev,wlev)
  else
    allocate(wtmp(size(pres)),stat=ier)
    if(ier /= 0) then
      write(stderr,'(2a,i4)') myname_,	&
	': allocate() error, stat =',ier
     call die(myname_)
    endif
    if(mall_ison()) call mall_mci(wtmp,myname)

    call slogtab(     NEAREST,nveclev,pveclev,	&
	size(pres),pres,ilev,wtmp)

    if(mall_ison()) call mall_mco(wtmp,myname)
    deallocate(wtmp)
  endif

end subroutine index1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: index0_ - index w.r.t. pveclev
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine index0_(pres,ilev,wlev)
      implicit none
      real,   intent(in) :: pres
      integer,intent(out):: ilev
      real,   intent(out),optional :: wlev

! !REVISION HISTORY:
! 	15Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::index0_'
  logical,parameter :: NEAREST=.true.

  real :: wtmp
  integer :: ier

  if(present(wlev)) then
    call slogtab(.not.NEAREST,nveclev,pveclev,	&
	1,pres,ilev,wlev)
  else
    call slogtab(     NEAREST,nveclev,pveclev,	&
	1,pres,ilev,wtmp)
  endif

end subroutine index0_
end module m_xTab_levs
