!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Dictionary80 - Hardwired 80 region Dictionary
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Dictionary80
      implicit none
      private	! except

      public :: Dictionary80_init	! create an object

      interface Dictionary80_init; module procedure init__; end interface

! !REVISION HISTORY:
! 	31Oct01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Dictionary80'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init__ - create a 80 region Dictionary
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init__(dict80,attr80,pseudo_population)
      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : Dictionary_init, Dictionary_compile
      use m_Dictionary,only : ptr_counts
      use m_Attributes,only : Attributes
      use m_Attributes,only : KR_SUBSET
      use m_Attributes,only : Attributes_init
      use m_Attributes,only : ptr_kr

      implicit none

      type(Dictionary),intent(out) :: dict80
      type(Attributes),intent(out) :: attr80
      Integer, Intent(In) :: pseudo_population(:)

! !REVISION HISTORY:
! 	31Oct01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init__'

  integer,pointer,dimension(:) :: p,q
  integer :: i,k

  call Attributes_init(attr80,KR_SUBSET,Size(pseudo_population))
  call Dictionary_init(dict80,KR_SUBSET,Size(pseudo_population))

	q => ptr_counts(dict80)
	p => ptr_kr(attr80)

  do i=1,size(pseudo_population)
    q(i)=pseudo_population(i)	! "volumes" of the regions
    p(i)=i			! pseudo region numbers
  end do

	nullify(p)
	nullify(q)

	! Building a link list out of an ordered (ascending) table.

  call Dictionary_compile(dict80)
  
end subroutine init__
end module m_Dictionary80
