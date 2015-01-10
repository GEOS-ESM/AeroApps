!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_LLGridDistribution - Distribution of a lat-long grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_LLGridDistribution
      implicit none
      private	! except

      public :: LLGridDistribution_init

      interface LLGridDistribution_init; module procedure	&
	init_; end interface

! !REVISION HISTORY:
!	31Oct01	- Jing Guo
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_LLGridDistribution'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - Initialize an object and other by-products
!
! !DESCRIPTION:
!
! !INTERFACE:
#include "assert.H"
    subroutine init_(dstr,attr_recv,man, attr_send,	&
	comm,wpower,nozero)

      use m_Distribution,only : Distribution
      use m_Attributes  ,only : Attributes, ptr_kr
      use m_MultiAccessNavigator,only : MultiAccessNavigator

      use m_Distribution,only : Distribution_init

      use m_Attributes  ,only : clean

      use m_Dictionary  ,only : Dictionary
      use m_Dictionary80,only : Dictionary80_init
      use m_Dictionary  ,only : Dictionary_clean

      use m_mpif90,      only : MP_TYPE, MP_SUM, MP_MAX
      use m_mall,        only : mall_mci, mall_mco, mall_ison
      use m_die,         only : assert_, MP_die

      implicit none

			! targetted output and by-products

      type(Distribution),intent(out) :: dstr	    ! targetted object
      type(Attributes)  ,intent(out) :: attr_recv   ! distributed data
      type(MultiAccessNavigator),intent(out) :: man ! block Navigators

			! Input to be used to determine the Distribution
			! and to be distributed.

      type(Attributes),intent(in) :: attr_send

      integer         ,intent(in) :: comm	   ! communicator

			! The partitioner may be defined by balancing
			! the (wpower) power of the counts.

      real   ,optional,intent(in) :: wpower

			! The partitioner may be required to make sure
			! there is no PE left with no data.

      logical,optional,intent(in) :: nozero

! !REVISION HISTORY:
!	31Oct01	- Jing Guo
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

  type(Dictionary),target :: dict80
  type(Attributes),target :: attr80

  Integer, Allocatable :: loc_pop(:), glo_pop(:)
  Integer :: i
  Integer :: ier
  Integer, Pointer :: kr(:)
  Integer :: mxkr ! maximum number of regions


     kr => ptr_kr(attr_send)

		! Allocate an array to contain the region populations and fill
  mxkr = maxval(kr)
  call MPI_allreduce((mxkr),mxkr,1,MP_TYPE(mxkr),MP_MAX,comm,ier)
     if(ier/=0) call MP_die(myname_,'MPI_allreduce()',ier)

     
  Allocate(loc_pop(mxkr),glo_pop(mxkr), STAT = ier)
     ALWAYS_ASSERT(ier == 0)
     If (mall_ison()) Then
        Call mall_mci(loc_pop, myname)
        Call mall_mci(glo_pop, myname)
     End If
     loc_pop = 0

     Do i = 1, Size(kr)
        loc_pop(kr(i)) = loc_pop(kr(i)) + 1
     End Do

     Call MPI_ALLreduce(loc_pop, glo_pop, Size(loc_pop), MP_TYPE(loc_pop),	&
          &	     MP_SUM, comm, ier)
          ALWAYS_ASSERT(ier == 0)

		! Construct a Dictionary with population information
	call Dictionary80_init(dict80,attr80, glo_pop)

           If (mall_ison()) Then
              Call mall_mco(loc_pop, myname)
              Call mall_mco(glo_pop, myname)
           End If
        Deallocate(loc_pop, glo_pop, STAT = ier)
           ALWAYS_ASSERT(ier == 0)

        
		! Determine the Distribution with this pre-determined
		! Dictionary.

  call Distribution_init(dstr,attr_recv,man, attr_send,	&
	dict80,attr80, comm,wpower=wpower,nozero=nozero)

		! Remove the Dictionary.

	call Dictionary_clean(dict80)
	call clean(attr80)

end subroutine init_

end module m_LLGridDistribution
