!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_xTab_lats - tabulate latitude entries
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_xTab_lats
      use rlat_imat
      implicit none
      private	! except

      public :: tab_lats
      public :: rmtab_lats
      public :: index_lats
      public :: lats_xtab

      public :: MXveclat
      public :: nveclat
      public :: veclats

      interface tab_lats; module procedure	&
	tab0_,	&
	tab1_
      end interface
      interface lats_xtab;  module procedure	&
	parxtab_,	&
	xtab_;  end interface
      interface rmtab_lats; module procedure rmtab_; end interface
      interface index_lats; module procedure	&
	index0_,	&
	index1_
      end interface

! !REVISION HISTORY:
! 	07Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_xTab_lats'

  integer,parameter :: JNP=91
  real,   parameter :: NLAT=+90.
  real,   parameter :: SLAT=-90.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: tab0_ - tabulate a default table
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine tab0_(nlev)
      implicit none
      integer,optional,intent(in) :: nlev ! by pass a compiler bug

! !REVISION HISTORY:
! 	07Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::tab0_'

  real :: dlat
  integer :: j

  nveclat=max(2,JNP)
  dlat=(NLAT-SLAT)/(nveclat-1)

  do j=1,nveclat
    veclats(j)=(j-1)*dlat+SLAT
  end do

end subroutine tab0_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: tab1_ - tabulate latitude entries
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine tab1_(ndat,rlat)
      implicit none
      integer,intent(in) :: ndat
      real,dimension(:),intent(in) :: rlat

! !REVISION HISTORY:
! 	07Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::tab1_'

  nveclat=0
  call merg_lats(JNP,ndat,rlat,MXveclat,nveclat,veclats)

end subroutine tab1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: xtab_ - extent the latitude table
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine xtab_(ndat,rlat)
      use m_die, only : die
      implicit none
      integer,intent(in) :: ndat
      real,dimension(:),intent(in) :: rlat

! !REVISION HISTORY:
! 	28Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::xtab_'

  if(size(rlat)<ndat) call die(myname_,		&
	'size(rlat)',size(rlat),'ndat',ndat	)

  call merg_lats(JNP,ndat,rlat,MXveclat,nveclat,veclats)

end subroutine xtab_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: parxtab_ - extent the latitude table from all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine parxtab_(rlat,comm)
      use m_die, only : die
      use m_Collector,only : Collector
      use m_Collector,only : Collector_init
      use m_Collector,only : clean
      use m_Collector,only : globalSize
      use m_CollectorComm,only : allgatherv
      implicit none
      real,dimension(:),intent(in) :: rlat
      integer          ,intent(in) :: comm

! !REVISION HISTORY:
! 	28Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::parxtab_'
  type(Collector) :: coll
  real,allocatable,dimension(:) :: alat
  integer :: ndat,nall
  integer :: ier

  ndat=size(rlat)
  call Collector_init(coll,ndat,comm)

  nall=globalSize(coll)
	allocate(alat(nall),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call allgatherv(rlat,alat,coll,comm,stat=ier)
		if(ier/=0) call die(myname_,'allgatherv()',ier)

  call merg_lats(JNP,nall,alat,MXveclat,nveclat,veclats)

	deallocate(alat,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

  call clean(coll)

end subroutine parxtab_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rmtab_ - clean the table (a place holder)
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

  nveclat=0

end subroutine rmtab_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: index1_ - index w.r.t. veclats
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine index1_(rlat,ilat,wlat)
      use m_stdio,only : stderr
      use m_die,  only : die
      implicit none
      real,   dimension(:),intent(in) :: rlat
      integer,dimension(:),intent(out):: ilat
      real,   dimension(:),intent(out),optional :: wlat

! !REVISION HISTORY:
! 	15Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::index1_'
  logical,parameter :: NEAREST=.true.

  real,dimension(:),allocatable :: wtmp
  integer :: ier

  if(present(wlat)) then
    call slintab(.not.NEAREST,nveclat,veclats,	&
	size(rlat),rlat,ilat,wlat)
  else
    allocate(wtmp(size(rlat)),stat=ier)
    if(ier /= 0) then
      write(stderr,'(2a,i4)') myname_,	&
	': allocate() error, stat =',ier
     call die(myname_)
    endif

    call slintab(     NEAREST,nveclat,veclats,	&
	size(rlat),rlat,ilat,wtmp)

    deallocate(wtmp)
  endif

end subroutine index1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: index0_ - index w.r.t. veclats
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine index0_(rlat,ilat,wlat)
      implicit none
      real,   intent(in) :: rlat
      integer,intent(out):: ilat
      real,   intent(out),optional :: wlat

! !REVISION HISTORY:
! 	15Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::index0_'
  logical,parameter :: NEAREST=.true.
  real    :: wtmp

  if(present(wlat)) then
    call slintab(.not.NEAREST,nveclat,veclats,	&
	1,rlat,ilat,wlat)
  else
    call slintab(     NEAREST,nveclat,veclats,	&
	1,rlat,ilat,wtmp)
  endif

end subroutine index0_
end module m_xTab_lats
