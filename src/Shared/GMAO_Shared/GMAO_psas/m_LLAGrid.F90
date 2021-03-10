!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_LLAGrid - A lat-long A-grid data structure
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_LLAGrid
      use m_aGrid, only : LLAGrid => aGrid
      use m_aGrid, only : tag_LLAGrid
      use m_aGrid, only : LINEAR,LEVELS
      implicit none
      private	! except

      public :: LLAGrid		! The class data structure
      public :: LLAGrid_init
      public :: LLAGrid_clean

      interface LLAGrid_init; module procedure	&
	init3d_
      end interface
      interface LLAGrid_clean; module procedure	&
	clean_
      end interface

! !REVISION HISTORY:
! 	08Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_LLAGrid'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init3d_ - initialize an aGrid structure
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init3d_(aG, mlon,mlat,mlev,levs,stat)
      use m_die,  only : die
      use m_stdio,only : stderr
      use m_String, only : init
      use m_mall, only : mall_mci,mall_ison
      use m_aGrid,only : aGrid_name
      implicit none
      type(LLAGrid),intent(out) :: aG

      integer,intent(in) :: mlon
      integer,intent(in) :: mlat
      integer,intent(in) :: mlev
      real,   intent(in) :: levs(:)
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	08Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init3d_'
  integer :: ier

  if(present(stat)) stat=0

  allocate(aG%vlon(3),aG%vlat(3),aG%vlev(mlev),stat=ier)
  if(ier /= 0) then
    write(stderr,'(2a,i4)') myname_,': allocate() error, stat =',ier
    if(.not.present(stat)) call die(myname_)
    stat=-1
  endif

	if(mall_ison()) then
	  call mall_mci(aG%vlon,aGrid_name)
	  call mall_mci(aG%vlat,aGrid_name)
	  call mall_mci(aG%vlev,aGrid_name)
	endif

  call init(aG%tag,tag_LLAGrid)

  aG%nlon=mlon
  aG%perx=.true.
  aG%mapx=LINEAR
  aG%vlon(1)=-180.
  aG%vlon(2)= 360.
  aG%vlon(3)=aG%vlon(2)/max(1,mlon)
  if(aG%perx) aG%vlon(3)= aG%vlon(2)/max(1,mlon-1)

  aG%nlat=mlat
  aG%pery=.false.
  aG%mapy=LINEAR
  aG%vlat(1)=-90.
  aG%vlat(2)=180.
  aG%vlat(3)=aG%vlat(2)/max(1,mlat)
  if(aG%pery) aG%vlat(3)= aG%vlat(2)/max(1,mlat-1)

  aG%nlev=mlev
  aG%mapz=LEVELS
  aG%vlev(:)=levs(:)

end subroutine init3d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init2d_ - initialize an aGrid structure
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init2d_(aG, mlon,mlat,stat)
      use m_die,  only : die
      use m_stdio,only : stderr
      use m_String, only : init
      use m_aGrid,only : aGrid_name
      use m_mall, only : mall_mci,mall_ison
      implicit none
      type(LLAGrid),intent(out) :: aG

      integer,intent(in) :: mlon
      integer,intent(in) :: mlat
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	08Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init2d_'
  integer :: ier

  if(present(stat)) stat=0

  allocate(aG%vlon(3),aG%vlat(3),stat=ier)
  if(ier /= 0) then
    write(stderr,'(2a,i4)') myname_,': allocate() error, stat =',ier
    if(.not.present(stat)) call die(myname_)
    stat=-1
  endif

	if(mall_ison()) then
	  call mall_mci(aG%vlon,aGrid_name)
	  call mall_mci(aG%vlat,aGrid_name)
	endif

  call init(aG%tag,tag_LLAGrid)

  aG%nlon=mlon
  aG%perx=.true.
  aG%mapx=LINEAR
  aG%vlon(1)=-180.
  aG%vlon(2)= 360.
  aG%vlon(3)=aG%vlon(2)/max(1,mlon)
  if(aG%perx) aG%vlon(3)= aG%vlon(2)/max(1,mlon-1)

  aG%nlat=mlat
  aG%pery=.false.
  aG%mapy=LINEAR
  aG%vlat(1)=-90.
  aG%vlat(2)=180.
  aG%vlat(3)=aG%vlat(2)/max(1,mlat)
  if(aG%pery) aG%vlat(3)= aG%vlat(2)/max(1,mlat-1)

  aG%nlev=-1	! should be TWO_D if this routine is to be used
  aG%mapz=-1

end subroutine init2d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a LLAGrid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(aG,stat)
      use m_die,   only : perr,die
      use m_aGrid, only : clean
      implicit none
      type(LLAGrid),   intent(inout) :: aG
      integer,optional,intent(out)   :: stat

! !REVISION HISTORY:
! 	03Jan00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  if(present(stat)) stat=0

	! Note that mall_mco() calls are made in aGrid_clean()

  call clean(aG,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'clean()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine clean_
end module m_LLAGrid
