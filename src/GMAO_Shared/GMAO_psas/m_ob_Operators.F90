!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ob_Operators - Coveriance operators of the obs. vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ob_Operators
      implicit none
      private	! except

      public :: ob_Operators_init
      public :: ob_Operators_clean

      public :: ob_qrx,ob_qry,ob_qrz
      public :: ob_qmx,ob_qmy,ob_qmz
      public :: ob_qlx,ob_qly
      public :: ob_sigU,ob_sigO,ob_sigFi
      public :: ob_xlev
      public :: ob_ktab,ob_jtab

      interface ob_Operators_init;   module procedure	&
	init_
      end interface
      interface ob_Operators_clean; module procedure	&
	clean_
      end interface

! !REVISION HISTORY:
! 	03Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ob_Operators'

	! covariance attributes (nobs)

  real,allocatable,dimension(:) :: ob_qrx,ob_qry,ob_qrz
  real,allocatable,dimension(:) :: ob_qmx,ob_qmy,ob_qmz
  real,allocatable,dimension(:) :: ob_qlx,ob_qly

  real,allocatable,dimension(:) :: ob_sigU,ob_sigO,ob_sigFi

  real,target,allocatable,dimension(:) :: ob_xlev	! extended levels
  integer,allocatable,dimension(:) :: ob_ktab,ob_jtab

  logical,save :: ob_Operators_defined=.false.

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize operators
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(nobs)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      integer,intent(in) :: nobs

! !REVISION HISTORY:
! 	06Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier

		! Sorting observations (-forecasts) is taking place now.

  if(ob_Operators_defined) call die(myname_,'multiple definitions')

  allocate( ob_qrx(nobs),ob_qry(nobs),ob_qrz(nobs),	&
	    ob_qmx(nobs),ob_qmy(nobs),ob_qmz(nobs),	&
	    ob_qlx(nobs),ob_qly(nobs),			&
	    ob_sigU(nobs),ob_sigO(nobs),		&
	    ob_sigFi(nobs),				&
	    ob_xlev(nobs),				&
	    ob_ktab(nobs),ob_jtab(nobs), stat=ier	)

	if(ier/=0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(ob_qrx, myname)
	  call mall_mci(ob_qry, myname)
	  call mall_mci(ob_qrz, myname)
	  call mall_mci(ob_qmx, myname)
	  call mall_mci(ob_qmy, myname)
	  call mall_mci(ob_qmz, myname)
	  call mall_mci(ob_qlx, myname)
	  call mall_mci(ob_qly, myname)
	  call mall_mci(ob_sigU,myname)
	  call mall_mci(ob_sigO,myname)
	  call mall_mci(ob_sigFi,myname)
	  call mall_mci(ob_xlev,myname)
	  call mall_mci(ob_ktab,myname)
	  call mall_mci(ob_jtab,myname)
	endif

  ob_Operators_defined=.true.
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

  if(.not.ob_Operators_defined) call die(myname_,'object undefined')

	if(mall_ison()) then
	  call mall_mco(ob_qrx, myname)
	  call mall_mco(ob_qry, myname)
	  call mall_mco(ob_qrz, myname)
	  call mall_mco(ob_qmx, myname)
	  call mall_mco(ob_qmy, myname)
	  call mall_mco(ob_qmz, myname)
	  call mall_mco(ob_qlx, myname)
	  call mall_mco(ob_qly, myname)
	  call mall_mco(ob_sigU,myname)
	  call mall_mco(ob_sigO,myname)
	  call mall_mco(ob_sigFi,myname)
	  call mall_mco(ob_xlev,myname)
	  call mall_mco(ob_ktab,myname)
	  call mall_mco(ob_jtab,myname)
	endif

  deallocate( ob_qrx,ob_qry,ob_qrz,	&
	      ob_qmx,ob_qmy,ob_qmz,	&
	      ob_qlx,ob_qly,		&
	      ob_sigU,ob_sigO,		&
	      ob_sigFi,			&
	      ob_ktab,ob_jtab,		&
	      ob_xlev,			&
	      stat=ier			)

	if(ier/=0) call die(myname_,'deallocate()',ier)

  ob_Operators_defined=.false.

end subroutine clean_
end module m_ob_Operators
