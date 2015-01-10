subroutine tabl_FEalpha()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: tabl_FEalpha - (to do)
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	05Mar96 - J. Guo	- (to do)
!_______________________________________________________________________

use config,	only	: lvmax
use m_die,	only	: die
use m_stdio,	only	: stdout
use m_output,only : output_FcstErrCov
use FEalpha_tabl
implicit none
!-----------------------------------------------------------------------
	! define local variables
  integer			:: istat	! status code
  character(len=*), parameter	:: myname='tabl_FEalpha'
  integer			:: i,l
!=======================================================================
  call rdPars(FEalpha_rsrc,FEalpha_type,FEalpha_desc,		&
		lvmax,FEalpha_nlev,FEalpha_plev,		&
		FEalpha_Mpar,FEalpha_npar, FEalpha_pars, istat	)

	! No default model.  A model must be specified

  if(istat /= 0) call die(myname,'rdPars()',istat)
!-----------------------------------------------------------------------
	! Echo the input

  if(.not.output_FcstErrCov) return

  l=max(len_trim(FEalpha_rsrc),1)
  write(stdout,'(a)') FEalpha_rsrc(1:l)

	! a description

  l=len_trim(FEalpha_type)
  if(l.gt.0) then
    l=max(len_trim(FEalpha_type),1)
    write(stdout,'(2x,a,x,a)') FEalpha_type,FEalpha_desc(1:l)
  endif

	! Table content.  The format is based on the fact that their
	! values are in the unit of wind errors.

  do l=1,FEalpha_nlev
    write(stdout,'(2x,f7.2,a,$)') FEalpha_plev(l),'  '
    do i=1,FEalpha_npar
      write(stdout,'(x,f7.2,$)') FEalpha_pars(i,l)
    end do
    write(stdout,*)
  end do
  write(stdout,'(a)') '::'		! table ends
!_______________________________________________________________________
end subroutine tabl_FEalpha
