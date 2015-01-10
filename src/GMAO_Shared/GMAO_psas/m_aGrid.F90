!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_aGrid - A simple grid data structure
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_aGrid
      use m_String, only : String
      implicit none
      private	! except

      public :: aGrid		! The class data structure
      public :: clean

      public :: LINEAR
      public :: LEVELS
      public :: tag_LLAGrid

      public :: aGrid_name

    type aGrid
      type(String) :: tag

      integer		:: nlon		! number of grid points
      integer		:: mapx		! the grid interval mapping
      logical		:: perx		! periodic?
      real,dimension(:),pointer :: vlon	! grid values

      integer		:: nlat		! numberof grid points
      integer		:: mapy		! the grid interval mapping
      logical		:: pery		! periodic?
      real,dimension(:),pointer :: vlat	! grid values

      integer		:: nlev		! number of grid points
      integer		:: mapz		! the grid interval mapping
      real,dimension(:),pointer	:: vlev	! grid values
    end type aGrid

    interface clean; module procedure clean_; end interface

! !REVISION HISTORY:
! 	08Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_aGrid'

  character(len=*),parameter :: aGrid_name=myname
  integer,parameter :: LINEAR=0
  integer,parameter :: LEVELS=1
  integer,parameter :: TWO_D =-2

  character(len=*),parameter :: tag_LLAGrid='LLAGrid'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean an aGrid object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(aG,stat)
      use m_stdio, only : stderr
      use m_die,   only : die
      use m_mall,  only : mall_mco,mall_ison
      use m_String,only : String_clean => clean
      implicit none
      type(aGrid),intent(inout) :: aG
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	08Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  if(present(stat)) stat=0

	if(mall_ison()) then
	  call mall_mco(aG%vlon,myname)
	  call mall_mco(aG%vlat,myname)
	  if(aG%nlev/=TWO_D) call mall_mco(aG%vlev,myname)
	endif

  deallocate(aG%vlon,aG%vlat,stat=ier)
  if(aG%nlev/=TWO_D) deallocate(aG%vlev,stat=ier)

  if(ier /= 0) then
    write(stderr,'(2a,i4)') myname_,': deallocate() error, ier =',ier
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  call String_clean(aG%tag)

  aG%nlon=-1
  aG%nlat=-1
  aG%nlev=-1

end subroutine clean_
end module m_aGrid
