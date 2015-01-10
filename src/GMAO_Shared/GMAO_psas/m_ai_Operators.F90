!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ai_Operators - Coveriance operators of the a.inc. vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ai_Operators
      implicit none
      private	! except

      public :: ai_Operators_init
      public :: ai_Operators_clean

      public :: ai_qrx,ai_qry,ai_qrz
      public :: ai_qmx,ai_qmy,ai_qmz
      public :: ai_qlx,ai_qly
      public :: ai_sigFi
      public :: ai_xlev
      public :: ai_ktab,ai_jtab

      interface ai_Operators_init; module procedure	&
	init_
      end interface

      interface ai_Operators_clean; module procedure	&
	clean_
      end interface

! !REVISION HISTORY:
! 	03Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ai_Operators'

	! covariance attributes (ninc)

  real   ,save,allocatable,dimension(:) :: ai_qrx,ai_qry,ai_qrz
  real   ,save,allocatable,dimension(:) :: ai_qmx,ai_qmy,ai_qmz
  real   ,save,allocatable,dimension(:) :: ai_qlx,ai_qly
  real   ,save,allocatable,dimension(:) :: ai_sigFi
  real   ,save,allocatable,dimension(:) :: ai_xlev
  integer,save,allocatable,dimension(:) :: ai_ktab,ai_jtab

  logical,save :: ai_Operators_defined=.false.

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize the data objects
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(ninc)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      integer,intent(in) :: ninc

! !REVISION HISTORY:
! 	06Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier

  if(ai_Operators_defined) call die(myname_,'multiple definitions')

  allocate( ai_qrx(ninc),ai_qry(ninc),ai_qrz(ninc),	&
	    ai_qmx(ninc),ai_qmy(ninc),ai_qmz(ninc),	&
	    ai_qlx(ninc),ai_qly(ninc),			&
	    ai_sigFi(ninc),				&
	    ai_xlev(ninc),				&
	    ai_ktab(ninc),ai_jtab(ninc), stat=ier	)

	if(ier/=0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(ai_qrx,myname)
	  call mall_mci(ai_qry,myname)
	  call mall_mci(ai_qrz,myname)
	  call mall_mci(ai_qmx,myname)
	  call mall_mci(ai_qmy,myname)
	  call mall_mci(ai_qmz,myname)
	  call mall_mci(ai_qlx,myname)
	  call mall_mci(ai_qly,myname)
	  call mall_mci(ai_sigFi,myname)
	  call mall_mci(ai_xlev,myname)
	  call mall_mci(ai_ktab,myname)
	  call mall_mci(ai_jtab,myname)
	endif

  ai_Operators_defined=.true.

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the data objects
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

  if(.not.ai_Operators_defined) call die(myname_,'object undefined')

	if(mall_ison()) then
	  call mall_mco(ai_qrx,myname)
	  call mall_mco(ai_qry,myname)
	  call mall_mco(ai_qrz,myname)
	  call mall_mco(ai_qmx,myname)
	  call mall_mco(ai_qmy,myname)
	  call mall_mco(ai_qmz,myname)
	  call mall_mco(ai_qlx,myname)
	  call mall_mco(ai_qly,myname)
	  call mall_mco(ai_sigFi,myname)
	  call mall_mco(ai_xlev,myname)
	  call mall_mco(ai_ktab,myname)
	  call mall_mco(ai_jtab,myname)
	endif

  deallocate( ai_qrx,ai_qry,ai_qrz,	&
	      ai_qmx,ai_qmy,ai_qmz,	&
	      ai_qlx,ai_qly,		&
	      ai_sigFi,			&
	      ai_xlev,			&
	      ai_ktab,ai_jtab,		&
	      stat=ier			)

	if(ier/=0) call die(myname_,'deallocate()',ier)

  ai_Operators_defined=.false.

end subroutine clean_
end module m_ai_Operators
