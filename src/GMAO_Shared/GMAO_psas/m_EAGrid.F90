!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_EAGrid - Q.E.A. 2 dimensional grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_EAGrid
      implicit none
      private	! except

      public :: EAGrid		! The class data structure
      public :: EAGrid_init
      public :: EAGrid_clean
      public :: EAGrid_ea2ll

    type EAGrid
      integer :: nEA
      real,pointer,dimension(:) :: lat
      real,pointer,dimension(:) :: lon
    end type EAGrid

    interface EAGrid_init ; module procedure init_ ; end interface
    interface EAGrid_clean; module procedure clean_; end interface
    interface EAGrid_ea2ll; module procedure ea2ll_; end interface

! !REVISION HISTORY:
! 	08Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_EAGrid'

	! For the time being, only a singleton is supported, because
	! of the invoking of the routines in gridxx.f, where a common
	! block (qea.h and now m_qea) is used.

  logical,save :: singleton_defined=.false.
  character(len=*),parameter :: EAGrid_rsrc="eagrid.rc"

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialzed an (EAGrid)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(eaG,nlon,nlat,rsrc)
      use m_inpak90,only : i90_loadf,i90_release
      use m_mall,only : mall_ison,mall_mci
      use m_die, only : die
      implicit none
      type(EAGrid),intent(out) :: eaG
      integer,intent(in) :: nlon
      integer,intent(in) :: nlat
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	08Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier
  integer :: mgrid

  if(singleton_defined) call die(myname_,'multiply defined')

  if(present(rsrc)) then
    call i90_LoadF(rsrc,ier)
    if(ier/=0) call die(myname_,'i90_LoadF("'//trim(rsrc)//'")',ier)
  else
    call i90_LoadF(EAGrid_rsrc,ier)
    if(ier/=0) call die(myname_,'i90_LoadF("'//EAGrid_rsrc//'")',ier)
  endif

  call gridxx0()	! which does qea_init().  See qea_clean() below.

  call i90_release(ier)
	if(ier/=0) then
	  if(present(rsrc)) then
	    call die(myname_,'i90_release("'//trim(rsrc)//'")',ier)
	  else
	    call die(myname_,'i90_release("'//EAGrid_rsrc//'")',ier)
	  endif
	endif

  mgrid=nlon*nlat
  allocate(eaG%lat(mgrid),eaG%lon(mgrid),stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(eaG%lat,myname)
	  call mall_mci(eaG%lon,myname)
	endif

  call eagrid_set(ier,eaG%lon,eaG%lat,mgrid,eaG%nEA,nlon,nlat)
	if(ier/=0) call die(myname_,'eagrid_set()',ier)

  singleton_defined=.true.

end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean an EAGrid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(eaG)
      use m_mall,only : mall_ison,mall_mco
      use m_die, only : die
      use m_qea, only : qea_clean
      implicit none
      type(EAGrid),intent(inout) :: eaG

! !REVISION HISTORY:
! 	08Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  if(.not.singleton_defined) call die(myname_,'undefined')

	if(mall_ison()) then
	  call mall_mco(eaG%lat,myname)
	  call mall_mco(eaG%lon,myname)
	endif

  deallocate(eaG%lat,eaG%lon,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

  eaG%nEA=0

  call qea_clean()

  singleton_defined=.false.

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ea2ll_ - interface ea2ll (gridxx.f)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine ea2ll_(zLL,nlon,nlat,zEA,nEA)
      use m_die,only : die
      implicit none
      real,dimension(:,:),intent(out) :: zLL
      integer,intent(in) :: nlon
      integer,intent(in) :: nlat
      real,dimension(:),intent(in) :: zEA
      integer,intent(in) :: nEA

! !REVISION HISTORY:
! 	08Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ea2ll_'

  if(.not.singleton_defined) call die(myname_,'undefined')

  call ea2ll(zLL,nlon,nlat,zEA,nEA)

end subroutine ea2ll_

end module m_EAGrid
