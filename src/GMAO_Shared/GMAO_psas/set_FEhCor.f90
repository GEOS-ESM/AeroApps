subroutine set_FEhCor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: set_FEhCor - set hfecH_tbl and hfecQ_tbl
!
! !INTERFACE:
!
! !DESCRIPTION:
!
! !EXAMPLES:
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
! 	10Jan96 - J. Guo	- programmed and added the prolog
!_______________________________________________________________________

use config, only : lvmax_hc, MXpar_hc, MX_fecH
use m_mpout,only : mpout,mpout_ison
use m_die,  only : die
use hfecH_tbl
use hfecQ_tbl

implicit none

character(len=*), parameter	:: myname='set_FEhCor'

integer i_CorType
integer i,l,km
integer istat

!-----------------------------------------------------------------------
! processing the FE correlation function from the resource file
	!---------------------------------------------------------------
	! hfecH_tbl

do km=1,MX_fecH
  call rdpars(	name_hfecH(km),type_hfecH(km),desc_hfecH(km),	&
		lvmax_hc,nlev_hfecH(km),plev_hfecH(1,km),	&
		MXpar_hc,npar_hfecH(km),pars_hfecH(1,1,km),	&
		istat						)

  if(istat.ne.0)	&
    call die(myname,'rdpars("'//trim(name_hfecH(km))//'")',istat)
end do

	!---------------------------------------------------------------
	! hfecQ_tbl

  call rdpars(	name_hfecQ,type_hfecQ,desc_hfecQ,		&
		lvmax_hc,nlev_hfecQ,plev_hfecQ,			&
		MXpar_hc,npar_hfecQ,pars_hfecQ,istat		)

  if(istat.ne.0)	&
    call die(myname,'rdpars("'//trim(name_hfecQ)//'")',istat)

!-----------------------------------------------------------------------
! echo the input

if(.not.mpout_ison()) return

	!---------------------------------------------------------------
	! hfecH_tbl

do km=1,MX_fecH
  write(mpout,'(a)') name_hfecH(km)

	! a description

  l=max(len_trim(desc_hfecH(km)),1)
  write(mpout,'(2x,a,x,a)') type_hfecH(km),desc_hfecH(km)(1:l)

	! table content

  do l=1,nlev_hfecH(km)
    if(abs(plev_hfecH(l,km)) >= 1.) then
      write(mpout,'(2x,i5,  2x,$)') nint(plev_hfecH(l,km))
    else
      write(mpout,'(2x,f5.3,2x,$)')      plev_hfecH(l,km)
    endif

    do i=1,npar_hfecH(km)
      if(abs(pars_hfecH(i,l,km)) > 10.) then
	write(mpout,'(x,f7.1,$)') pars_hfecH(i,l,km)
      else
	write(mpout,'(x,f7.4,$)') pars_hfecH(i,l,km)
      endif
    end do

    write(mpout,*)
  end do

	! table ends
  write(mpout,'(a)') '::'
end do

	!---------------------------------------------------------------
	! hfecH_tbl

  l=max(len_trim(name_hfecQ),1)
  write(mpout,'(a)') name_hfecQ(1:l)

	! a description

  l=max(len_trim(desc_hfecQ),1)
  write(mpout,'(2x,a,x,a)') type_hfecQ,desc_hfecQ(1:l)

	! table content

  do l=1,nlev_hfecQ

    if(abs(plev_hfecQ(l)) >= 1.) then
      write(mpout,'(2x,i5,  2x,$)') nint(plev_hfecQ(l))
    else
      write(mpout,'(2x,f5.2,2x,$)')      plev_hfecQ(l)
    endif

    do i=1,npar_hfecQ
      if(abs(pars_hfecQ(i,l)) > 10.) then
	write(mpout,'(x,f7.1,$)') pars_hfecQ(i,l)
      else
	write(mpout,'(x,f7.4,$)') pars_hfecQ(i,l)
      endif
    end do

    write(mpout,*)
  end do

	! table ends
  write(mpout,'(a)') '::'
!-----------------------------------------------------------------------
end subroutine set_FEhCor
!.
