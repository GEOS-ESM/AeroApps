!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_PartMatx - Parallel (Matx) partitioner
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_PartMatx
      implicit none
      private	! except

      public :: PartMatx		! The class data structure
      public :: ChooseMethod            ! Determine partition method from env

      public :: PartMatx_SCATTER	! spread expensive blocks (optimal load-balance)
      public :: PartMatx_COMPACT	! crude locality generator
      public :: PartMatx_RECTANG	! optimal for asymmetric matvec (good cache-reuse)
      public :: PartMatx_HILBERT	! Hilbert (space-filling) schedule for compact and cache-reuse

      interface PartMatx; module procedure part_; end interface
      interface ChooseMethod; module procedure choose_; end interface

! !REVISION HISTORY:
! 	18Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_PartMatx'

  integer,parameter :: PartMatx_SCATTER=0	! default, sorted costs
  integer,parameter :: PartMatx_COMPACT=1	! compact but tangled
  integer,parameter :: PartMatx_RECTANG=2	! clear cuts in row/col
  integer,parameter :: PartMatx_HILBERT=3	! clear cuts in Hilbert curve
  integer,parameter :: METHODS(4) = (/ PartMatx_SCATTER, PartMatx_COMPACT, &
	&                              PartMatx_RECTANG, PartMatx_HILBERT /)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: part_ - a parallel (Matx) partitioner driver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine part_(Matx,iMatx,nMatx,comm,TotalCost,myCost,method)

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : ptr_rAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_stdio,    only : stderr
      use m_die,      only : die

      use m_scatPartMatx,only : scatPartMatx => partMatx
      use m_compPartMatx,only : compPartmatx => partMatx
      use m_rectPartMatx,only : rectPartMatx => partMatx
      use m_HilbertPartMatx,only : HilbertPartMatx => partMatx

      implicit none

      type(AttrVect),intent(in) :: Matx	! a list of matrix blocks

      integer,dimension(:),intent(inout) :: iMatx ! indices to Matx
      integer,             intent(inout) :: nMatx ! <=size(iMatx)

      integer,intent(in) :: comm

      real,   optional,intent(out) :: TotalCost
      real,   optional,intent(out) :: myCost
      integer,optional,intent(in)  :: method	! 0, 1, 3, 4

! !REVISION HISTORY:
! 	17Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::part_'

  real,pointer,dimension(:,:) :: rAttr

  integer :: i,l
  integer :: icost
  real    :: TotalCost_,myCost_
  integer :: method_

	! Set the default partitioning algorithm
  method_=PartMatx_SCATTER
        ! If optional method is specified, use it ...
  if(present(method)) method_=method
       

	if( .not. Any(method_ == (/ METHODS /))) Then
    
	  write(stderr,'(2a,i4)') myname_,	&
		': unknown value, method =',method_
	  call die(myname_)
	endif

	! Check indices

  icost=AttrVect_indexRA(Matx,'cost')

	if(icost<=0) then
	  write(stderr,'(2a,i4)') myname_,': invalid, icost =',icost
	  call die(myname_)
	endif

	! Estimate the _total_ cost to be divided

		rAttr => ptr_rAttr(Matx)

  TotalCost_=0.
  do i=1,nMatx
    l=iMatx(i)
    TotalCost_=TotalCost_+rAttr(icost,l)
  end do
  if(present(TotalCost)) TotalCost=TotalCost_

		nullify(rAttr)

  myCost_=TotalCost_

  select case(method_)
  case(PartMatx_SCATTER)

	! Use an algorithm that the _costlier_ tasks will be 
	! distributed first, without compactness consideration.

    call scatPartMatx(Matx,iMatx,nMatx, myCost_, comm)

  case(PartMatx_COMPACT)

	! Use an algorithm that produces work partitions with
	! their _surfaces_ of loosely minimized.

    call compPartMatx(Matx,iMatx,nMatx, myCost_, comm)

  case(PartMatx_RECTANG)

	! Use an algorithm that preserves order

    call rectPartMatx(Matx,iMatx,nMatx, myCost_, comm)

  case(PartMatx_HILBERT)

	! Use an algorithm that preserves order in Hilbert sense

    call HilbertPartMatx(Matx,iMatx,nMatx, myCost_, comm)

  end select

  if(present(myCost)) myCost=myCost_

end subroutine part_

#include "assert.H"
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: choose_ - select method for partitioning
!
! !DESCRIPTION:
!
! !INTERFACE:

Integer Function choose_(env_var, default)
  Use m_chars, only : uppercase
  use m_die, only : assert_, die
  Character(Len=*), Intent(In) :: env_var
  Integer, Intent(In) :: default

! !REVISION HISTORY:
! 	05Nov01 - Tom Clune <clune@sgi.com> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::choose_'
  Character(Len=255) :: env_val
  Integer :: method

     ! Default method
  ASSERT(Any(default == METHODS))
  method = default

  Call getenv(env_var, env_val)
  if ( len_trim(env_val) .gt. 0 ) then
     env_val = uppercase(env_val)
     Select Case(env_val)
     Case ('SCATTER')
        method = PartMatx_SCATTER
     Case ('COMPACT')
        method = PartMatx_COMPACT
     Case ('RECTANG')
        method = PartMatx_RECTANG
     Case ('HILBERT')
        method = PartMatx_HILBERT
     Case Default
        call die(myname_,'no such method: '//env_val)
     End Select
 End if
 choose_ = method

End Function Choose_

  
end module m_PartMatx
!.
