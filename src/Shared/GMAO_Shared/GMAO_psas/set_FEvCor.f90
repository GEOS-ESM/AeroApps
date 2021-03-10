subroutine set_FEvCor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: set_FEvCor - set vfecH_tbl and vfecQ_tbl
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

use config, only : lvmax_vc,MX_fecH
use m_mpout,only : mpout,mpout_ison
use m_die,  only : die
use vfecH_tbl
use vfecQ_tbl

implicit none

character(len=*), parameter	:: myname='set_FEvCor'

integer i,l,istat,km

!-----------------------------------------------------------------------
! processing all function types from the resource file
	!---------------------------------------------------------------
	! vfecH_tbl

do km=1,MX_fecH
  call rdvctbl(	name_vfecH(km),type_vfecH(km),desc_vfecH(km),	&
		lvmax_vc,nlev_vfecH(km),plev_vfecH(1,km),	&
		corr_vfecH(1,1,km),istat	)

  if(istat.ne.0)	&
    call die(myname,'rdvctbl("'//trim(name_vfecH(km))//'")',istat)
end do

	!---------------------------------------------------------------
	! vfecQ_tbl

  call rdvctbl(	name_vfecQ,type_vfecQ,desc_vfecQ,		&
		lvmax_vc,nlev_vfecQ,plev_vfecQ,corr_vfecQ,istat	)

  if(istat.ne.0)	&
    call die(myname,'rdvctbl("'//trim(name_vfecQ)//'")',istat)

!-----------------------------------------------------------------------
! echo the input before the final "mirroring".

if(mpout_ison()) then

	!---------------------------------------------------------------
	! vfecH_tbl label

do km=1,MX_fecH
  write(mpout,'(a)') name_vfecH(km)

	! a description

  l=max(len_trim(desc_vfecH(km)),1)
  write(mpout,'(2x,a,x,a)') type_vfecH(km),desc_vfecH(km)(1:l)

	! table content

  do l=1,nlev_vfecH(km)
    if(abs(plev_vfecH(l,km)) >= 1.) then
      write(mpout,'(2x,i5,  2x,$)') nint(plev_vfecH(l,km))
    else
      write(mpout,'(2x,f5.3,2x,$)')      plev_vfecH(l,km)
    endif

    do i=1,l
      write(mpout,'(x,f4.2,$)') corr_vfecH(i,l,km)
    end do

    write(mpout,*)
  end do

	! table ends
  write(mpout,'(a)') '::'
end do

	!---------------------------------------------------------------
	! vfecQ_tbl label

  write(mpout,'(a)') name_vfecQ

	! a description

  l=max(len_trim(desc_vfecQ),1)
  write(mpout,'(2x,a,x,a)') type_vfecQ,desc_vfecQ(1:l)

	! table content

  do l=1,nlev_vfecQ
    if(abs(plev_vfecQ(l)) >= 1.) then
      write(mpout,'(2x,i5,  2x,$)') nint(plev_vfecQ(l))
    else
      write(mpout,'(2x,f5.3,2x,$)')      plev_vfecQ(l)
    endif

    do i=1,l
      write(mpout,'(x,f4.2,$)') corr_vfecQ(i,l)
    end do

    write(mpout,*)
  end do

	! table ends
  write(mpout,'(a)') '::'
endif

!-----------------------------------------------------------------------
! Make the correlation matrices easier to interpolate by making them
! linear around their main diagonals.
	!---------------------------------------------------------------
	! vfecH_tbl

do km=1,MX_fecH
  do l=1,nlev_vfecH(km)
    do i=1,l-1
      corr_vfecH(l,i,km)=2.-corr_vfecH(i,l,km)
    end do
  end do
end do

	!---------------------------------------------------------------
	! vfecQ_tbl

  do l=1,nlev_vfecQ
    do i=1,l-1
      corr_vfecQ(l,i)=2.-corr_vfecQ(i,l)
    end do
  end do
!-----------------------------------------------------------------------
end subroutine set_FEvCor
!.
