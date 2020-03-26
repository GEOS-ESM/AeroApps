!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Checks - check attributes of an operator
!
! !DESCRIPTION:
!
! !INTERFACE:
!#include "regime.H"

    module m_Checks
      implicit none
      private	! except

!      public :: CheckAttr	! Check block attributes
      public :: CheckNorm	! Check vector l2-norm

      interface CheckNorm; module procedure	&
	l2norm1_loc,	&
	l2norm2_loc,	&
	l2norm1_,	&
	l2norm2_; end interface

! !REVISION HISTORY:
! 	23Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Checks'

contains
#ifdef _TODO_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: checkAttr - check type-blocks
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine checkAttr(name,nav,krnav,ktnav,rlon,rlat,rlev)
      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_Navigator2,only : getExtent
      use m_mpout,only : mpout,mpout_log,mpout_flush

      implicit none

      character(len=*)    ,intent(in) :: name
      type(Navigator)     ,intent(in) :: nav
      integer,dimension(:),intent(in) :: krnav
      integer,dimension(:),intent(in) :: ktnav
      real   ,dimension(:),intent(in) :: rlon
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlev

! !REVISION HISTORY:
! 	23Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::checkAttr'
  integer :: nAtom,nType
  integer :: inav,lc,le,ln
  integer :: inext,idspl,irest,istep
  integer :: lnext,ldspl

  nType=lsize(nav)
  call getExtent(nav,count=nAtom)

  call mpout_log(myname,'lsize('//trim(name)//')=',nAtom,	&
	showrank=.true.)
  call mpout_log(myname,'nType('//trim(name)//')=',nType,	&
	showrank=.true.)

  call mpout_log(myname,'('//trim(name)//')',showrank=.true.)

	! Show a header

!23.|....1....|....2....|....3....|....4....|....5....|....6....|....7..
!.inav......lc......le...ln...kr...kt.......l....rlon....rlat....rlev

  write(mpout,'(1x,a)')	&
' iblk      lc      le   ln   kr   kt       l    rlon    rlat    rlev'

	! List blocks as well vector contents with a specifiable
	! interval.

  istep=1
  idspl=0
  inext=idspl+istep
  irest=nAtom-1

  do inav=1,nType
    call get(nav,inav,lc=lc,le=le,ln=ln,displ=ldspl)

    if(lc+ln/=le+1) then

      write(mpout,'(i6,2i8,i5,a)') inav,lc,le,ln,'**'

    else
      if(ln<=0) then
      
	write(mpout,'(i6,2i8,i5)') inav,lc,le,ln

      else

	lnext=ldspl+inext-idspl
	if(lnext > le) then
		! If no listing for this segment, list the block.

	  write(mpout,'(i6,2i8,3i5)') inav,		&
		lc,le,ln,krNav(inav),ktNav(inav)

	else
		! otherwise, list data.

	  do while(lnext<=le)
	    write(mpout,'(i6,2i8,3i5,i8,3f8.2)') inav,	&
		lc,le,ln,krNav(inav),ktNav(inav),	&
		lnext,rlon(lnext),rlat(lnext),rlev(lnext)

		! Adjust istep.  If there are already many, increase
		! istep.  If there are less left, decrease istep.
		! However, don't become an infinit loop.

	    if(inext==istep*10) istep=istep*10
	    do while(irest<istep)
	      istep=istep/10
	    end do
	    if(istep<1) istep=1

	    inext=inext+istep
	    irest=irest-istep
	    lnext=lnext+istep
	  end do
	endif

	idspl=idspl+ln

      endif
    endif

  end do

  call mpout_log(myname_,'('//trim(name)//')',showrank=.true.)
end subroutine checkAttr
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: l2norm1_ - check type-blocks
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine l2norm1_(name,vect,comm)
      use m_parDOT,only : parNRM2
      use m_mpout ,only : mpout_log
      implicit none

      character(len=*) ,intent(in) :: name
      real,dimension(:),intent(in) :: vect
      integer          ,intent(in) :: comm

! !REVISION HISTORY:
! 	23Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::l2norm1_'
  real :: nrm2

  nrm2=parNRM2(vect,comm)
  call mpout_log(myname,'l2-norm('//trim(name)//')=',	&
	nrm2,showrank=.false.)

end subroutine l2norm1_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: l2norm2_ - check type-blocks
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine l2norm2_(name,vect,comm)
      use m_parDOT,only : parNRM2
      use m_mpout ,only : mpout_log
      implicit none

      character(len=*)   ,intent(in) :: name
      real,dimension(:,:),intent(in) :: vect
      integer            ,intent(in) :: comm

! !REVISION HISTORY:
! 	23Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::l2norm2_'
  integer :: i
  real,dimension(size(vect,2)) :: nrm2

  nrm2(:)=parNRM2(vect,comm)

  do i=1,size(vect,2)
    call mpout_log(myname,'l2-norm('//trim(name)//')=',	&
	nrm2(i),showrank=.false.)
  end do

end subroutine l2norm2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: l2norm1_loc - check type-blocks
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine l2norm1_loc(name,vect)
      use m_mpout ,only : mpout_log
      implicit none

      character(len=*) ,intent(in) :: name
      real,dimension(:),intent(in) :: vect

! !REVISION HISTORY:
! 	23Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::l2norm1_loc'

  call mpout_log(myname,'l2-norm('//trim(name)//')=',	&
	sqrt(dot_product(vect(:),vect(:))), showrank=.true.)

end subroutine l2norm1_loc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: l2norm2_loc - check type-blocks
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine l2norm2_loc(name,vect)
      use m_parDOT,only : parNRM2
      use m_mpout ,only : mpout_log
      implicit none

      character(len=*)   ,intent(in) :: name
      real,dimension(:,:),intent(in) :: vect

! !REVISION HISTORY:
! 	23Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::l2norm2_loc'
  integer :: i

  do i=1,size(vect,2)
    call mpout_log(myname,'l2-norm('//trim(name)//')=',	&
	sqrt(dot_product(vect(:,i),vect(:,i))), showrank=.true.)
  end do

end subroutine l2norm2_loc

end module m_Checks
