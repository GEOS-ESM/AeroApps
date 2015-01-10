!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: stdv_OEqc - update an observation err. stdv. table for QC
!
! !INTERFACE:

  subroutine stdv_OEqc( n_kx,n_lt,mxlev, mlev,plevs, sigO_tabl)

    use config, only : ktmax,ktus,ktvs,ktslp,ktHH,ktUU,ktVV,ktQQ
    use config, only : pres4slp
    use m_die,  only : die
    implicit none

  integer, intent(in) :: n_kx		! Number of insturments (kxmax)
  integer, intent(in) :: n_lt		! Number of variable columes
  integer, intent(in) :: mxlev		! table size parameter
  integer, intent(in) :: mlev		! number of levels
  real,    intent(in) :: plevs(mlev)	! presure levels

  real,    intent(inout) :: sigO_tabl(n_kx,n_lt,mxlev+1) ! OE. stdv.

	! Define the output dependent table references for QC.  See
	! OI "obserr.data" for detail.

    integer, parameter :: ltps =2
    integer, parameter :: ltus =3
    integer, parameter :: ltvs =4

    integer, parameter :: ltQQ =1
    integer, parameter :: ltHH =2
    integer, parameter :: ltUU =3
    integer, parameter :: ltVV =4

! !DESCRIPTION:
!
!   stdv_OEqc() returns observation error standard deviations for
!   GEOS/DAS-QC.
!
! !SYSTEM ROUTINES:
!
!   getenv(3f) is required to check if there is an user specified
!   resource file to use, through an environment variable PSASRC.
!   If the variable has not been set, a default filename (psas.rc) is
!   selected.
!
! !REVISION HISTORY:
!
! 	15Nov96 - J. Guo	- prototyped.  And confirmed the 
!		interface design with Genia Brin and David Lamich.
!	17Nov96 - J. Guo	- finished code and tested.
!	18Nov96 - J. Guo	- redefined the interface to fit
!		setobse() of OI.
!
!_______________________________________________________________________
!#define _TRACE	!

character(len=*),parameter :: myname="stdv_OEqc"

	! Local workspace

integer             :: n_lst
integer,allocatable :: kx_lst(:),kt_lst(:),lv_lst(:),lt_lst(:)
real,   allocatable :: pr_lst(:), sigRc(:),sigRu(:)

	! ... code goes here

  call check_args_()	! arguments consistancy checking
  call init_tabl_()	! reading tables for the resource file
  call set_list_()	! creating a list of "observation"
  call derive_stdv_()	! interpolating the tables
!=======================================================================
contains
!----------------------------------------
subroutine check_args_()
  use m_die,only : die
  use m_stdio,only : stderr
  implicit none

character(len=*),parameter :: myname_=myname//"::check_args_"

  if(n_lt.lt.max(ltps,ltus,ltvs,ltHH,ltUU,ltVV,ltQQ)) then
    write(stderr,*) myname_,": expecting ",	&
      max(ltps,ltus,ltvs,ltHH,ltUU,ltVV,ltQQ),	&
      ", but n_lt = ",n_lt
    call die(myname_)
  endif

  if(mxlev.lt.mlev) then
    write(stderr,*) myname_,": expecting ",mlev,", but mxlev = ",mxlev
    call die(myname_)
  endif

end subroutine check_args_
!----------------------------------------
subroutine init_tabl_()
  use m_die,  only : die
  use m_psasrc,only : psasrc_open,psasrc_close
  implicit none
integer ier
character(len=*),parameter :: myname_=myname//"::init_tabl_"

		! Loading the resource handle

  call psasrc_open(stat=ier)
	if(ier/=0) call die(myname,'psasrc_open()',ier)

  call ktname0()		! initialize kttabl.h
  call kxname0()		! initialize kxtabl.h
  call set_OEclas()	! OEclass_tbl

  call psasrc_close(stat=ier)
	if(ier/=0) call die(myname,'psasrc_close()',ier)

end subroutine init_tabl_
!----------------------------------------
subroutine set_list_()
  use m_die,  only : die
  implicit none

  integer istat,n,kx,kt,lv
  character(len=*), parameter :: myname_=myname//"::set_list_"

  n_lst=n_kx*ktmax*mlev		! it will create more memory than needed,
				! but is blind to the actual values of
				! the variable type codes (ktslp, etc.).

  allocate(kx_lst(n_lst),kt_lst(n_lst),lv_lst(n_lst),lt_lst(ktmax), &
	   pr_lst(n_lst), sigRc(n_lst), sigRu(n_lst), stat=istat    )
	if(istat/=0) call die(myname_,'allocate()',istat)

		! Configuration dependent indexing (for QC only)

  do kt=1,ktmax
    select case(kt)
      case(ktslp)
	lt_lst(kt)=ltps
      case(ktus)
	lt_lst(kt)=ltus
      case(ktvs)
	lt_lst(kt)=ltvs
      case(ktHH)
	lt_lst(kt)=ltHH
      case(ktuu)
	lt_lst(kt)=ltUU
      case(ktvv)
	lt_lst(kt)=ltVV
      case(ktqq)
	lt_lst(kt)=ltQQ
      case default
	call die(myname_,'invalid kt',kt)
    end select
  end do

  n=0
  do kx=1,n_kx
    do kt=1,ktmax
      if(kt.eq.ktslp.or.kt.eq.ktus.or.kt.eq.ktvs) then
	n=n+1
	kx_lst(n)=kx
	kt_lst(n)=kt
	lv_lst(n)=mxlev+1	! Output table reference
	pr_lst(n)=pres4slp
      else
	do lv=1,mlev
	  n=n+1
	  kx_lst(n)=kx
	  kt_lst(n)=kt
	  lv_lst(n)=lv		! Output table reference
	  pr_lst(n)=plevs(lv)
	end do
      endif
    end do
  end do
  
  n_lst=n
end subroutine set_list_
!----------------------------------------
subroutine derive_stdv_()
  implicit none
  integer i,kx,kt,lv,lt

  call intp_sigO(n_lst,kx_lst,kt_lst,pr_lst,sigRc,sigRu)

		! Selective filling the table
  do i=1,n_lst
    if(sigRc(i).ge.0. .or. sigRu(i).ge.0.) then

      kx=kx_lst(i)
      kt=kt_lst(i)
      lv=lv_lst(i)
      lt=lt_lst(kt)

      sigO_tabl(kx,lt,lv)=sqrt( sigRc(i)*sigRc(i) + sigRu(i)*sigRu(i) )

    endif
  end do

  deallocate(kx_lst,kt_lst,lv_lst,lt_lst,pr_lst,sigRc,sigRu)
end subroutine derive_stdv_
!----------------------------------------
end subroutine stdv_OEqc
!.
