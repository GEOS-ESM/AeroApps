!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ob_Navigator - Navigator of the obs. vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ob_Navigator
      use m_subBlocks, only : subBlocks
      implicit none
      private	! except

      public :: ob_Nav
      public :: ob_Indx

      public :: ob_krlc
      public :: ob_krln
      public :: ob_ktlc
      public :: ob_ktln

      public :: ob_Navigator_init
      public :: ob_Navigator_clean

      interface ob_Navigator_init; module procedure	&
	init_
      end interface
      interface ob_Navigator_clean; module procedure	&
	clean_
      end interface


! !REVISION HISTORY:
! 	03Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ob_Navigator'

	! The navigator data

  type(subBlocks) :: ob_Nav
  integer,allocatable,dimension(:)   :: ob_Indx	! (nobs)

	! Original navigator data (nxkr) or (nxkt,nxkr)

  integer,allocatable,dimension(:)   :: ob_krlc,ob_krln
  integer,allocatable,dimension(:,:) :: ob_ktlc,ob_ktln

  logical,save :: ob_Navigator_defined=.false.

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize the navigator
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(nxkt,nxkr,nobs)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      integer,intent(in) :: nxkt
      integer,intent(in) :: nxkr
      integer,intent(in) :: nobs

! !REVISION HISTORY:
! 	06Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier

  if(ob_Navigator_defined) call die(myname_,'multiple definitions')

  allocate( ob_krlc(nxkr),ob_ktlc(nxkt,nxkr),	&
	    ob_krln(nxkr),ob_ktln(nxkt,nxkr),	&
	    ob_Indx(nobs),	 stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(ob_krlc,myname)
	  call mall_mci(ob_ktlc,myname)
	  call mall_mci(ob_krln,myname)
	  call mall_mci(ob_ktln,myname)
	  call mall_mci(ob_Indx,myname)
	endif

  ob_Navigator_defined=.true.

end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the data object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mco
      implicit none

! !REVISION HISTORY:
! 	06Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  if(.not.ob_Navigator_defined) call die(myname_,'object undefined')

	if(mall_ison()) then
	  call mall_mco(ob_Indx,myname)
	  call mall_mco(ob_ktln,myname)
	  call mall_mco(ob_krln,myname)
	  call mall_mco(ob_ktlc,myname)
	  call mall_mco(ob_krlc,myname)
	endif

  deallocate(ob_krlc,ob_ktlc,ob_krln,ob_ktln,ob_Indx,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

  ob_Navigator_defined=.false.

end subroutine clean_

end module m_ob_Navigator
