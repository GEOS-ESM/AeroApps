!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ai_Navigator - Navigator of the analysis increment vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ai_Navigator
      use m_subBlocks, only : subBlocks
      implicit none
      private	! except

      public :: ai_Nav
      public :: ai_Indx

      public :: ai_krlc
      public :: ai_krln
      public :: ai_ktlc
      public :: ai_ktln

      public :: ai_Navigator_init
      public :: ai_Navigator_clean

      interface ai_Navigator_init; module procedure	&
	init_
      end interface
      interface ai_Navigator_clean; module procedure	&
	clean_
      end interface

! !REVISION HISTORY:
! 	03Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ai_Navigator'

	! The navigator data

  type(subBlocks) :: ai_Nav
  integer,allocatable,dimension(:)   :: ai_Indx	! (ninc)

	! Original navigator data (nxkr) or (nxkt,nxkr)

  integer,allocatable,dimension(:)   :: ai_krlc,ai_krln
  integer,allocatable,dimension(:,:) :: ai_ktlc,ai_ktln

  logical,save :: ai_Navigator_defined=.false.

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize the data object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(nxkt,nxkr,ninc)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      integer,intent(in) :: nxkt
      integer,intent(in) :: nxkr
      integer,intent(in) :: ninc

! !REVISION HISTORY:
! 	06Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier

  if(ai_Navigator_defined) call die(myname_,'multiple definitions')

  allocate( ai_krlc(nxkr),ai_ktlc(nxkt,nxkr),	&
	    ai_krln(nxkr),ai_ktln(nxkt,nxkr),	&
	    ai_Indx(ninc),	stat=ier	)
	if(ier/=0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(ai_krlc,myname)
	  call mall_mci(ai_ktlc,myname)
	  call mall_mci(ai_krln,myname)
	  call mall_mci(ai_ktln,myname)
	  call mall_mci(ai_Indx,myname)
	endif

  ai_Navigator_defined=.true.

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

  if(.not.ai_Navigator_defined) call die(myname_,'undefined object')

	if(mall_ison()) then
	  call mall_mco(ai_krlc,myname)
	  call mall_mco(ai_ktlc,myname)
	  call mall_mco(ai_krln,myname)
	  call mall_mco(ai_ktln,myname)
	  call mall_mco(ai_Indx,myname)
	endif

  deallocate( ai_krlc,ai_ktlc,ai_krln,ai_ktln,ai_Indx,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

  ai_Navigator_defined=.false.

end subroutine clean_

end module m_ai_Navigator
