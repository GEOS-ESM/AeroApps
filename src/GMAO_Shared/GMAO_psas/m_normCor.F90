!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_normCor - Normalization factors of a multi-var correlation
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_normCor
      use config, only	: kmat_HGHT
      implicit none
      private	! except

      public :: kmat_HGHT
      public :: intp_normCor

      interface intp_normCor; module procedure intp_; end interface

! !REVISION HISTORY:
! 	04Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_normCor'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intp_ - interpolation for a vector of level values
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intp_(nlev,rlev,norm)

      use m_stdio,only : stderr
      use m_die,  only : die
      use m_mall, only : mall_ison,mall_mci,mall_mco
      use rlev_imat, only : MXveclev,nveclev
      use m_xTab_levs, only : index_levs
      use m_mvcorF_blox,only : norm_DD

      implicit none
      integer,intent(in)  :: nlev	! the size of rlev[] and norm[]
      real,   intent(in)  :: rlev(:)	! levels
      real,   intent(out) :: norm(:)	! normalization factors

		! Accessing the common block

! !REVISION HISTORY:
! 	23Feb96 - J. Guo	- as dervsigF_upD()
! 	04Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intp_'

	! Local workspace

  integer,allocatable	:: ilev(:)
  real,   allocatable	:: wlev(:)

	! Loval vars.

  real		:: wt
  integer	:: i,lv,ierr

  logical, parameter	:: nearest =.true.
!-----------------------------------------------------------------------
	! Allocate workspace

	allocate(ilev(nlev),wlev(nlev), stat=ierr )
	if(ierr.ne.0) then
	  write(stderr,'(2a,i3)') myname_,		&
	    ': allocate() error, stat =',ierr
	  call die(myname_)
	endif
	if(mall_ison()) then
	  call mall_mci(ilev,myname)
	  call mall_mci(wlev,myname)
	endif

!-----------------------------------------------------------------------
	! Index rlev(1:nlev) on norm_DD(:) vs.  pveclev(:) table
	! with log-linear weights

  call index_levs(rlev,ilev,wlev)

  do i=1,nlev
    lv=ilev(i)
    wt=wlev(i)

		! conditional linear interpolation

    norm(i)=norm_DD(lv,kmat_HGHT)
    if(lv.lt.nveclev)	&
	norm(i)=norm(i) + wt*(norm_DD(lv+1,kmat_HGHT)-norm(i))

  end do

	if(mall_ison()) then
	  call mall_mco(ilev,myname)
	  call mall_mco(wlev,myname)
	endif
  deallocate(ilev,wlev)

end subroutine intp_
end module m_normCor
!.
