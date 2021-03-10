!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ai_Attributes - Attributes of the analysis-increment vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ai_Attributes
      implicit none
      private	! except

      public :: ai_Attributes_init
      public :: ai_Attributes_clean

      public :: ai_lat
      public :: ai_lon
      public :: ai_lev
      public :: ai_kt
      public :: ai_vinc

      interface ai_Attributes_init; module procedure	&
	init_
      end interface

      interface ai_Attributes_clean; module procedure	&
	clean_
      end interface

! !REVISION HISTORY:
! 	03Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ai_Attributes'

	! Observation attributes (ninc)

  real   ,save,allocatable,dimension(:) :: ai_lat,ai_lon,ai_lev
  integer,save,allocatable,dimension(:) :: ai_kt
  real   ,save,allocatable,dimension(:,:) :: ai_vinc

  logical,save :: ai_Attributes_defined=.false.

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

    subroutine init_(ninc,indx,vlat,vlon,vlev,vkts)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      integer,intent(in) :: ninc
      integer,dimension(:),intent(in) :: indx
      real   ,dimension(:),intent(in) :: vlat
      real   ,dimension(:),intent(in) :: vlon
      real   ,dimension(:),intent(in) :: vlev
      integer,dimension(:),intent(in) :: vkts

! !REVISION HISTORY:
! 	06Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier
  integer :: i,l

  if(ai_Attributes_defined) call die(myname_,'multiple definitions')

  allocate(ai_lat(ninc),ai_lon(ninc),ai_lev(ninc),ai_kt(ninc),	&
	ai_vinc(ninc,1),stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(ai_lat,myname)
	  call mall_mci(ai_lon,myname)
	  call mall_mci(ai_lev,myname)
	  call mall_mci(ai_kt ,myname)
	  call mall_mci(ai_vinc,myname)
	endif

  do i=1,ninc
    l=indx(i)
    ai_lat(i)=vlat(l)
    ai_lon(i)=vlon(l)
    ai_lev(i)=vlev(l)
    ai_kt (i)=vkts(l)
  end do

  ai_Attributes_defined=.true.

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

  if(.not.ai_Attributes_defined) call die(myname_,'object undefined')

	if(mall_ison()) then
	  call mall_mco(ai_lat,myname)
	  call mall_mco(ai_lon,myname)
	  call mall_mco(ai_lev,myname)
	  call mall_mco(ai_kt ,myname)
	  call mall_mco(ai_vinc,myname)
	endif

  deallocate(ai_lat,ai_lon,ai_lev,ai_kt,ai_vinc, stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

  ai_Attributes_defined=.false.

end subroutine clean_

end module m_ai_Attributes
