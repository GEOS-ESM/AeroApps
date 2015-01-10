!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_simplePart - simple partitioner
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_simplePart
      implicit none
      private	! except

      public :: simplePart	! Partition parameters
      public :: simplePart_comm	! Partition parameters

! !REVISION HISTORY:
! 	17Sep01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_simplePart'

#include "assert.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: simplePart - a simple chunked partition
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine simplePart(total,nParts,myPart,count,displ)
      use m_die,only : assert_
      implicit none
      integer,intent(in) :: total	! the total size
      integer,intent(in) :: nParts	! number of chunks
      integer,intent(in) :: myPart	! specify my part (0:nParts-1)
      integer,intent(out) :: count	! size of this part
      integer,intent(out) :: displ	! displacement of this part

! !REVISION HISTORY:
! 	17Sep01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::simplePart'
  integer :: resid

	ASSERT(nParts> 0)
	ASSERT(myPart>=0 .and. myPart< nParts)

  resid=mod(total,nParts)

  count=total/nParts
  if(myPart< resid) count=count+1

  displ=count*myPart
  if(myPart>=resid) displ=displ+resid

end subroutine simplePart
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: simplePart_comm - a simple chunked partition
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine simplePart_comm(total,comm,count,displ)
      use m_die,only : MP_die
      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      implicit none
      integer,intent(in) :: total	! the total size
      integer,intent(in) :: comm	! communicator
      integer,intent(out) :: count	! size of this part
      integer,intent(out) :: displ	! displacement of this part

! !REVISION HISTORY:
! 	17Sep01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::simplePart_comm'
  integer :: nPEs,myID
  integer :: resid
  integer :: ier

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)
  call MP_comm_size(comm,nPEs,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)

  resid=mod(total,nPEs)

  count=total/nPEs
  if(myID< resid) count=count+1

  displ=count*myID
  if(myID>=resid) displ=displ+resid

end subroutine simplePart_comm

end module m_simplePart
