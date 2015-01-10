!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_DictionaryTable - Dictionary generated from a table
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_DictionaryTable
      implicit none
      private	! except

      public :: DictionaryTable_init	! create an object

      interface DictionaryTable_init; module procedure	&
      	init__; end interface

! !REVISION HISTORY:
! 	31Oct01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_DictionaryTable'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init__ - create a Dictionary from a simple table
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init__(dict,attr,pseudo_population)
      use m_Dictionary,only : Dictionary
      use m_Dictionary,only : Dictionary_init, Dictionary_compile
      use m_Dictionary,only : ptr_counts
      use m_Attributes,only : Attributes
      use m_Attributes,only : KR_SUBSET
      use m_Attributes,only : Attributes_init
      use m_Attributes,only : ptr_kr

      implicit none

      type(Dictionary),intent(out) :: dict
      type(Attributes),intent(out) :: attr
      Integer, Intent(In) :: pseudo_population(:)

! !REVISION HISTORY:
! 	31Oct01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init__'

  integer,pointer,dimension(:) :: p,q
  integer :: i,k

  call Attributes_init(attr,KR_SUBSET,Size(pseudo_population))
  call Dictionary_init(dict,KR_SUBSET,Size(pseudo_population))

	q => ptr_counts(dict)
	p => ptr_kr(attr)

  do i=1,size(pseudo_population)
    q(i)=pseudo_population(i)	! "volumes" of the regions
    p(i)=i			! pseudo region numbers
  end do

	nullify(p)
	nullify(q)

	! Building a link list out of an ordered (ascending) table.

  call Dictionary_compile(dict)
  
end subroutine init__
end module m_DictionaryTable
