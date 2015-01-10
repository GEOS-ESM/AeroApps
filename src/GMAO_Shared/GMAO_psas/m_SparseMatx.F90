!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_CmprSprsRowCoord - Compressed Sparse Row Coordinate scheme
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_CmprSprsRowCoord
      implicit none
      private	! except

      public :: CmprSprsRowCoord	! The class data structure

    type CmprSprsRowCoord
      private
      integer :: ni
      integer,pointer,dimension(:) :: icoord
      type(Navigator) :: iNav_j
      integer :: nij
      integer,pointer,dimension(:) :: jcoord
    end type CmprSprsRowCoord

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_CmprSprsRowCoord'

subroutine toCSCC(cscc,csrc)

  pjcoord => ptr_jcoord(csrc)
  countj(:)=0
  do i=1,nrow(csrc)
    call get(ptr_navigator(csrc),i,row=ic,lc=jlc,le=jle)

    do j=jlc,jle
      jc=pjcoord(j)			! ordered by icoord()
      counti(jc)=counti(jc)+1
    end do
  end do

	! Sorting

  nelem=lsize(csrc)
  call indexSet(indx(1:nelem))
  call indexSort(indx(1:nelem),ptr_jcoord(csrc))
  icoord(1:nelem)=icoord(indx(1:nelem))
  pjcoord => ptr_jcoord(csrc)
  jcoord(1:nelem)=pjcoord(indx(1:nelem))

  ncol=?
  do jcol=1,ncol
    if(counti(jcol)
  end do

end module m_CmprSprsRowCoord
