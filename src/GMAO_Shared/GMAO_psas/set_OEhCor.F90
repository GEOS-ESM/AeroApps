!#define	_TRACE	!
subroutine set_OEhCor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: set_OEhCor - set hoecH_tbl
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

use config, only : kxmax, lvmax, MX_hoecH
use m_mpout,only : mpout,mpout_ison
use m_die,  only : die
use OEclass_tbl
use hoecH_tbl
use m_output, only : output_ObsErrCov

implicit none

integer nhc,i,l,ip
integer kx,kc
integer i_CorType
integer istat

character(len=len(hoecH_c)) kxhoecH(kxmax)
character(len=len('ObsErr*hCor_HH-::')+len(name_hoecH)) rc_tmp

character(len=*), parameter	:: myname='set_OEhCor'

include "kxtabl.h"

#ifdef	_TRACE
	_TRACE	write(mpout,'(2a)') myname,': enterend'
#endif

	! Remove undefined first from hoecH_c
  nhc=0
  do i=1,nOEclas
    if(hoecH_c(i).ne.' ' .and. hoecH_c(i).ne.'-') then
      nhc=nhc+1
      kxhoecH(nhc)=hoecH_c(i)
    endif
  end do

	! Tabulate all types of hoecH_c.

#ifdef	_TRACE
	_TRACE	write(mpout,'(2a)') myname,': calling tabSlist()'
#endif
  n_hoecH=0
  call tabSlist(nhc,kxhoecH,MX_hoecH,n_hoecH,name_hoecH)

	! Index the function types to the table

  do kx=1,kxmax
    kc=i_kxclas(kx)		! its OE class index
    kxhoecH(kx)='-'
    if(kc.gt.0) kxhoecH(kx)=hoecH_c(kc)	! its hoecH_c name
  end do

#ifdef	_TRACE
	_TRACE	write(mpout,'(2a)') myname,': calling inxSlist()'
#endif
  call inxSlist(n_hoecH,name_hoecH,kxmax,kxhoecH,i_hoecHc)

  !-do kx=1,kxmax
  !-  kc=i_hoecHc(kx)
  !-  if(kc.gt.0) write(mpout,*) kx,kc,':',kxhoecH(kx),name_hoecH(kc)
  !-end do

	! processing all function types from the resource file

  do i=1,n_hoecH
    l=len_trim(name_hoecH(i))
    write(rc_tmp,'(3a)') 'ObsErr*hCor_HH-',name_hoecH(i)(1:l),'::'

#ifdef	_TRACE
	_TRACE	write(mpout,'(4a)') myname,			&
	_TRACE	  ': calling rdpars("',name_hoecH(i)(1:l),'")'
#endif

    call rdpars(rc_tmp,type_hoecH(i),desc_hoecH(i),	&
	lvmax,nlev_hoecH(i),plev_hoecH(1,i),		&
	MXpar_hc,npar_hoecH(i),pars_hoecH(1,1,i),istat	)

    if(istat.ne.0)	&
      call die(myname,'rdpars("'//trim(rc_tmp)//'")',istat)
  end do

  if(.not.output_ObsErrCov) return
  if(.not.mpout_ison()) return

	! echo the input

  do i=1,n_hoecH
		! table label
    l=len_trim(name_hoecH(i))
    write(rc_tmp,'(3a)') 'ObsErr*hCor_HH-',name_hoecH(i)(1:l),'::'

    l=max(len_trim(rc_tmp),1)
    write(mpout,'(a)') rc_tmp(1:l)

		! a description

    l=max(len_trim(desc_hoecH(i)),1)
    write(mpout,'(2x,a,x,a)') type_hoecH(i),desc_hoecH(i)(1:l)

		! table content

    write(mpout,*)
    do l=1,nlev_hoecH(i)
			! the level
      if(abs(plev_hoecH(l,i)) >= 1.) then
        write(mpout,'(2x,i5,  2x,$)') nint(plev_hoecH(l,i))
      else
        write(mpout,'(2x,f5.3,2x,$)')      plev_hoecH(l,i)
      endif

      do ip=1,npar_hoecH(i)
        if(abs(pars_hoecH(ip,l,i)) > 10.) then
	  write(mpout,'(x,f7.1,$)') pars_hoecH(ip,l,i)
        else
	  write(mpout,'(x,f7.4,$)') pars_hoecH(ip,l,i)
        endif
      end do

      write(mpout,*)
    end do

		! table ends
    write(mpout,'(a)') '::'
  end do
end
