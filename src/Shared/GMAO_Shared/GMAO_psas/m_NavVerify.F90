!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_NavVerify
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_NavVerify
      implicit none
      private	! except

      public :: NavVerify

      interface NavVerify; module procedure verify_; end interface

! !REVISION HISTORY:
! 	24May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_NavVerify'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: verify_ - verify the attributes for this operator
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine verify_(nav,krNav,ktNav)
      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_stdio,only : stderr
      use m_die  ,only : die
      implicit none
      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

! !REVISION HISTORY:
! 	24May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::verify_'
  include "ktmax.h"

  integer :: inav,kt
  integer ::  im, ip
  integer :: krm,krp
  integer :: ktm,ktp
  integer :: lnm,lnp
  logical :: invalid,failed

	! Consistency checking: sizes of u and v in each partition
	! must be the same

  failed=.false.
  do inav=1,lsize(nav)
    kt=ktNav(inav)

    select case(kt)
    case (ktUU)
      ktm=kt
      ktp=ktVV
      im=inav
      ip=inav+1
    case (ktUs)
      ktm=kt
      ktp=ktVs
      im=inav
      ip=inav+1
    case (ktVV)
      ktm=ktUU
      ktp=kt
      im=inav-1
      ip=inav
    case (ktVs)
      ktm=ktUs
      ktp=kt
      im=inav-1
      ip=inav
    case default
      ktm=kt
      ktp=kt
    end select

    if(ktm/=ktp) then
      invalid = ktNav(im)/=ktm .or. ktNav(ip)/=ktp
      if(.not.invalid) invalid = krNav(im)/=krNav(ip)
      if(.not.invalid) then
        call get(nav,im,ln=lnm)
	call get(nav,ip,ln=lnp)
	invalid = lnm/=lnp
      endif

      if(invalid) then

	write(stderr,'(2a,4(a,i5,a,i5,a))')		&
	  myname_,': unmatched blocks',			&
	       ' (',      im ,',',      ip ,')',	&
	  ' kr = (',krNav(im),',',krNav(ip),')',	&
	  ' kt = (',ktNav(im),',',ktNav(ip),')',	&
	  ' ln = (',     lnm ,',',     lnp ,')'

	failed=.true.
      endif
    endif
  end do

  if(failed) call die(myname_)

end subroutine verify_
end module m_NavVerify
