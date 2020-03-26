subroutine set_OEvCor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: set_OEvCor - set voecH_tbl
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

use config, only : kxmax,lvmax,MX_voecH
use m_mpout,only : mpout,mpout_ison
use m_die,  only : die
use OEclass_tbl
use voecH_tbl
use m_output, only : output_ObsErrCov

implicit none

integer i,l,k,nvc
integer kx,kc
integer istat

character(len=len(voecH_c)) kxvoecH(max(kxmax,MX_voecH))
character(len=len('ObsErr*vCor_HH-::')+len(name_voecH)) rc_tmp

character(len=*), parameter	:: myname='set_OEvCor'

include "kxtabl.h"

!-----------------------------------------------------------------------
	! Remove undefined first

  nvc=0
  do i=1,nOEclas
    if(voecH_c(i).ne.' '.and.voecH_c(i).ne.'-') then
      nvc=nvc+1
      kxvoecH(nvc)=voecH_c(i)
    endif

    if(voecH_u(i).ne.' '.and.voecH_u(i).ne.'-') then
      nvc=nvc+1
      kxvoecH(nvc)=voecH_u(i)
    endif
  end do

	! Tabulate all types of voecH_c and voecH_u.

  n_voecH=0
  call tabSlist(nvc,kxvoecH,MX_voecH,n_voecH,name_voecH)

!-----------------------------------------------------------------------
	! Index the function types to the table

  do kx=1,kxmax
    kc=i_kxclas(kx)		! its OE class index
    kxvoecH(kx)='-'
    if(kc.gt.0) kxvoecH(kx)=voecH_c(kc)	! its voecH_c name
  end do
  call inxSlist(n_voecH,name_voecH,kxmax,kxvoecH,i_voecHc)

  !-do kx=1,kxmax
  !-  kc=i_voecHc(kx)
  !-  if(kc.gt.0) write(mpout,*) kx,kc,':',kxvoecH(kx),name_voecH(kc)
  !-end do

  do kx=1,kxmax
    kc=i_kxclas(kx)		! its OE class index
    kxvoecH(kx)='-'
    if(kc.gt.0) kxvoecH(kx)=voecH_u(kc)	! its voecH_u name
  end do
  call inxSlist(n_voecH,name_voecH,kxmax,kxvoecH,i_voecHu)

  !-do kx=1,kxmax
  !-  kc=i_voecHu(kx)
  !-  if(kc.gt.0) write(mpout,*) kx,kc,':',kxvoecH(kx),name_voecH(kc)
  !-end do

!-----------------------------------------------------------------------
	! processing all function types from the resource file

  do i=1,n_voecH
    l=len_trim(name_voecH(i))
    write(rc_tmp,'(3a)') 'ObsErr*vCor_HH-',name_voecH(i)(1:l),'::'

    call rdvctbl(rc_tmp,type_voecH(i),desc_voecH(i),	&
	lvmax_vc,nlev_voecH(i),plev_voecH(1,i),		&
	corr_voecH(1,1,i),istat				)

    if(istat.ne.0)	&
      call die(myname,'rdvctbl("'//trim(rc_tmp)//'")',istat)
  end do

!-----------------------------------------------------------------------
	! Echo the input

if(output_ObsErrCov.and.mpout_ison()) then
  do i=1,n_voecH
		! table label
    l=len_trim(name_voecH(i))
    write(rc_tmp,'(3a)') 'ObsErr*vCor_HH-',name_voecH(i)(1:l),'::'

    l=max(len_trim(rc_tmp),1)
    write(mpout,'(a)') rc_tmp(1:l)

		! a description

    l=max(len_trim(desc_voecH(i)),1)
    write(mpout,'(2x,a,x,a)') type_voecH(i),desc_voecH(i)(1:l)

		! table content

    do l=1,nlev_voecH(i)
      if(abs(plev_voecH(l,i)) >= 1.) then
        write(mpout,'(2x,i5,  2x,$)') nint(plev_voecH(l,i))
      else
        write(mpout,'(2x,f5.3,2x,$)')      plev_voecH(l,i)
      endif

      do k=1,l
        write(mpout,'(x,f4.2,$)') corr_voecH(k,l,i)
      end do

      write(mpout,*)
    end do

		! table ends
    write(mpout,'(a)') '::'
  end do
endif

!-----------------------------------------------------------------------
	! Let the correlation matrices easier to interpolate by making
	! them linear around their main diagonals.

  do i=1,n_voecH
    do k=1,nlev_voecH(i)
      do l=1,k-1
        corr_voecH(k,l,i)=2.-corr_voecH(l,k,i)
      end do
    end do
  end do
!-----------------------------------------------------------------------
end
!.
