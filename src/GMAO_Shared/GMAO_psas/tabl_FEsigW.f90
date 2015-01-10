subroutine tabl_FEsigW()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: tabl_FEsigW - initialize FEsigW_tabl from the resource file
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
!	08Oct99	- Jing Guo
!		. Added the controls of _defined and _setbyRC, to
!		  support module m_sigDCWindErr.
! 	23Apr96 - J. Guo	- (to do)
!_______________________________________________________________________

use config, only : lvmax
use m_mpout,only : mpout,mpout_ison,mpout_log
use m_die,  only : die
use m_output,only : output_FcstErrCov

use FEsigW_tabl

implicit none

character(len=*), parameter	:: myname='tabl_FEsigW'

integer i,l
integer istat

!-----------------------------------------------------------------------

  if(FEsigW_defined .and. .not.FEsigW_setbyRC) then

    call mpout_log(myname,	&
	'using a user specified "'//trim(FEsigW_rsrc)//'"')

  else
	! Search the resource file for the table, FEsigW_tabl

    call rdpars( FEsigW_rsrc,FEsigW_type,FEsigW_desc,	&
	lvmax, FEsigW_nlev,FEsigW_plev,			&
	FEsigW_Mpar, FEsigW_npar,FEsigW_pars, istat		)

    if(istat.ne.0)	&
      call die(myname,'rdpars("'//trim(FEsigW_rsrc)//'")',istat)

    FEsigW_defined=.true.
    FEsigW_setbyRC=.true.

    call mpout_log(myname,	&
	'using a resource based "'//trim(FEsigW_rsrc)//'"')
  endif

!-----------------------------------------------------------------------
	! Echo the input

  if(.not.output_FcstErrCov) return
  if(.not.mpout_ison()) return

  write(mpout,'(a)') trim(FEsigW_rsrc)

	! a description

  l=len_trim(FEsigW_type)
  if(l.gt.0) then
    l=max(len_trim(FEsigW_type),1)
    write(mpout,'(2x,a,x,a)') FEsigW_type,FEsigW_desc(1:l)
  endif

	! Table content.  The format is based on the fact that their
	! values are in the unit of wind errors.

  do l=1,FEsigW_nlev
    write(mpout,'(2x,f7.2,a,$)') FEsigW_plev(l),'  '
    do i=1,FEsigW_npar
      write(mpout,'(x,f7.2,$)') FEsigW_pars(i,l)
    end do
    write(mpout,*)
  end do
  write(mpout,'(a)') '::'		! table ends
!-----------------------------------------------------------------------
end subroutine tabl_FEsigW
!.
