!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_xtenlev - extended levels (meta-data)
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_xtenlev
      implicit none
      private	! except

      public :: xtenlev_set

      interface xtenlev_set; module procedure set_; end interface

! !REVISION HISTORY:
! 	28Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_xtenlev'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: set_ - set extended level meta-data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine set_(n,xlev,rlev,kts)
      use m_ktList,only : ktUs,ktVs,ktslp
      use config,  only : pres4slp
      use m_die,only : die,perr
      implicit none
      integer,intent(in) :: n
      real,   dimension(:),intent(out) :: xlev
      real,   dimension(:),intent(in ) :: rlev
      integer,dimension(:),intent(in ) :: kts

! !REVISION HISTORY:
! 	28Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::set_'
  integer :: i

  if( size(xlev) < n	.or.	&
      size(rlev) < n	.or.	&
      size(kts)  < n	) then

    if(size(xlev)<n) call perr(myname_,'size(xlev)',size(xlev),'n',n)
    if(size(rlev)<n) call perr(myname_,'size(rlev)',size(rlev),'n',n)
    if(size(kts )<n) call perr(myname_,'size(kts )',size(kts ),'n',n)

    call die(myname_)
  endif

  do i=1,n
    select case(kts(i))
    case (ktUs,ktVs,ktslp)
      xlev(i)=pres4slp
    case default
      xlev(i)=rlev(i)
    end select
  end do

end subroutine set_

end module m_xtenlev
