!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: PSAS_Sigs - define error stdv. of observations
!
! !INTERFACE:

  subroutine PSAS_Sigs(	nobs,rlat,rlon,rlev,kx,kt,		&
			sigF,sigOc,sigOu			)

  use m_xRSRC_sigFi, only : fetch_sigFi
	use m_FEsigFi_tabl, only : clean_sigFi => FEsigFi_clean

  use m_xTab_sigFi,  only : tab_sigFi,rmtab_sigFi

  use m_xTab_levs,   only : tab_levs, rmtab_levs
  use m_xTab_lats,   only : tab_lats, rmtab_lats
  use m_die,	     only : die
  use m_zeit,	     only : zeit_ci,zeit_co

  use config,        only : pres4slp
  use config,        only : ktslp,ktUs,ktVs,ktQQ
  implicit none

        ! Observations

  integer, intent(in) :: nobs           ! no. of observations

  real,    intent(in) :: rlat(nobs)     ! lat. of observations
  real,    intent(in) :: rlon(nobs)     ! lon. of observations
  real,    intent(in) :: rlev(nobs)     ! lev. (pres) of observations
  integer, intent(in) ::   kx(nobs)     ! instrument code
  integer, intent(in) ::   kt(nobs)     ! variable type code

  real,    intent(out):: sigF(nobs)     ! Forecast error stdv. of obs.
  real,    intent(out):: sigOc(nobs)    ! Correlated Obs. error stdv.
  real,    intent(out):: sigOu(nobs)    ! Uncorrelated Obs. error stdv.

! !DESCRIPTION:
!
!   PSAS\_Sigs() returns error standard deviations for a given
!   observation data set.  The returned values includes values for the
!   forecast errors and two components of the observation errors.
!
!   In case that the observation error covariance information for a
!   given observations error class (C for _correlated_ or U for
!   _uncorrelated_) is not available, a negative value will be returned.
!   It is the user's reponsibility to check if such values exist.
!
! !SYSTEM ROUTINES:
!
!   getenv(3f) is required to check if there is an user specified
!   resource file to use, through an environment variable PSASRC.
!   If the variable has not been set, a default filename (psas.rc) is
!   selected.
!
! DESIGN ISSUES:
!
!   - Should error variances be returned instead of error stdv.?
!     I am hoping in GEOS-DAS 2.1, PSAS will provide calls return
!     variances not stdv., such that, we can name the interfaces as
!     diag(R,..) for the diagonal of an observation error covariance
!     matrix R, etc..  We will let users to take the square-root.
!
!     In the future interface design:
!
!       diag(Ru,obs), diag(Rc,obs), diag(R,obs), diag(Pf,obs), and
!       diag(Pf,grid) through interface overloading.
!
!     FORTRAN 77 binding with different interface names will be
!     provided accordingly.
!
! !REVISION HISTORY:
! 	12Dec96 - J. Guo	- initial interface and implementation
! 	01Oct98 - Jing Guo <guo@thunder> - revised the prolog for the
!					  latest protex() rules
!EOP ___________________________________________________________________

character(len=*),parameter :: myname='PSAS_Sigs'

	! Local workspace

real,dimension(:),allocatable :: tlev

	! Nothing needs to be done.  Not even considered an error

  if(nobs == 0) return

	! Initializes tables from the resource file

	call zeit_ci(myname)
  call init_()

  call derive_OEstdv_()	! interpolating from the sigO tables
  call derive_FEstdv_()	! interpolating from the sigF tables

  call clean_()
	call zeit_co(myname)

!=======================================================================
contains
!----------------------------------------
subroutine init_()
  use m_die,only : die
  use m_psasrc,only : psasrc_open,psasrc_close
  use m_mall, only : mall_ison,mall_mci
  implicit none
integer :: i,ierr
character(len=*),parameter :: myname_=myname//'::init_'

		! Loading the resource handle

  call psasrc_open(stat=ierr)
	if(ierr/=0) call die(myname_,'psasrc_open()',ierr)

	! Define those common blocks (A part of OE information?)

  call ktname0()		! initialize kttabl.h

!-----------------------------------------------------------------------
	! Define observation error information tables

  call kxname0()		! initialize kxtabl.h
  call set_OEclas()	! OEclass_tbl

!-----------------------------------------------------------------------
	! Define forecast error information tables

  call fetch_sigFi()
  call set_FEhCor()
  call set_FEvCor()
  call tabl_FEsigW()
  call tabl_FEalpha()

!-----------------------------------------------------------------------

	! Release PSASRC
  call psasrc_close(stat=ierr)
	if(ierr /= 0) call die(myname_,'psasrc_close()',ierr)

	! Making a temporary pressure level array

  allocate(tlev(1:nobs), stat=ierr)
	if(ierr /= 0) call die(myname_,'allocate()',ierr)

	if(mall_ison()) call mall_mci(tlev,myname)

  do i=1,nobs
    tlev(i)=rlev(i)
    if(kt(i)==ktslp.or.kt(i)==ktUs.or.kt(i)==ktVs) tlev(i)=pres4slp
  end do

end subroutine init_
!----------------------------------------
subroutine derive_OEstdv_()
  implicit none

  character(len=*), parameter :: myname_=myname//'::derive_OEstdv_'

  call intp_sigO(nobs,kx,kt,tlev,sigOc,sigOu)

end subroutine derive_OEstdv_
!----------------------------------------
subroutine derive_FEstdv_()
  use m_stdv_FE, only : stdv_FE
  use m_mvcorF_bldr,only : mvcorF_init
  use m_mvcorF_bldr,only : mvcorF_clean
  implicit none

  character(len=*), parameter :: myname_=myname//'::derive_FEstdv_'

	! Tabulizing FE cov. objects

  call tab_levs(nobs,rlev)
  call tab_lats(nobs,rlat)

  call mvcorF_init()	! was "call set_fecHH()"
  call imat_alpha()
  call imat_sigW()
  call tab_sigFi()

	! Deriving stdv, from FE cov. lookup-table objects

  call stdv_FE(nobs,kt,rlat,rlon,tlev, sigF)

	! Remove all the lookup tables (at least trying to)

  call rmtab_sigFi()

	! Other rmtab_() calls yet to be defined

  call rmtab_lats()
  call rmtab_levs()
  call mvcorF_clean()

end subroutine derive_FEstdv_

subroutine clean_()
  use m_mall,only : mall_ison,mall_mco
  implicit none

  character(len=*),parameter :: myname_=myname//'::clean_'

	if(mall_ison()) call mall_mco(tlev,myname)

  deallocate(tlev)

  call clean_sigFi()

  ! There is nothing you can do to those common blocks
end subroutine clean_

!=======================================================================
end subroutine PSAS_Sigs
!.
