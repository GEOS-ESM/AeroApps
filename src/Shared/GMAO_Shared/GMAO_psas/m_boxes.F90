!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_boxes - Data selection
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_boxes
      implicit none
      private	! except

      public :: boxes_init
      public :: boxes_restrict
      public :: boxes_clean

    interface boxes_init;  module procedure init_;  end interface
    interface boxes_clean; module procedure clean_; end interface
    interface boxes_restrict; module procedure	&
	restrict_
    end interface

! !REVISION HISTORY:
! 	20Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_boxes'

  integer,parameter :: NBMAX=100
  integer,save      :: nboxes=0
  real,dimension(2,6,NBMAX),save :: boxes

  logical,save :: boxes_defined = .false.

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize boxes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_()
      use m_die,only : die
      implicit none

! !REVISION HISTORY:
! 	20Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

  if(boxes_defined) call die(myname_,'multiple object definition')

  call setbox(NBMAX,nboxes,boxes)

  boxes_defined=.true.

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean boxes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      use m_die,only : die
      implicit none

! !REVISION HISTORY:
! 	20Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  if(.not.boxes_defined) call die(myname_,'undefined object')

  boxes_defined=.false.

end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: restrict_ - limit data through boxes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine restrict_(nobs,qlst,rlat,rlon,rlev,tims,kxs,kts)
      use m_mpout,only : mpout
      use m_die,  only : die
      implicit none
      integer,intent(in) :: nobs
      logical,dimension(:),intent(inout) :: qlst

      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon
      real   ,dimension(:),intent(in) :: rlev
      real   ,dimension(:),intent(in) :: tims
      integer,dimension(:),intent(in) :: kxs
      integer,dimension(:),intent(in) :: kts

! !REVISION HISTORY:
! 	20Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::restrict_'

  logical,parameter :: verbose=.false.

  integer :: ier

  if(.not.boxes_defined) call die(myname_,'undefined object')

	! ..Restrict observations only to those "within" at least
	! one of "hyper-boxes", defined by lat/lon/pres/kx/kt/time.
	! Remove data outside the "hyper-boxes" by push them to the
	! end of the list and reset `nobs' to the size of the front
	! part of the list.

!     Restrict the choice to data within a union of lat/lon/lev boxes.
!     We want to set  kl(n) = kl(n).and.(box1.or.box2.or....or.boxk).
!     ---------------------------------------------------------------

  call LLBOXES ( verbose,mpout,nboxes,boxes,	&
	nobs,rlon,rlat,rlev,kxs,kts,tims,qlst	)

end subroutine restrict_

end module m_boxes
