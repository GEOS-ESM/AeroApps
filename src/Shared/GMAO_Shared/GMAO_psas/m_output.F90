!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_output - A control module of the PSAS stdout (log)
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_output
      implicit none
      private	! except

      public :: output_ObsErrCov
      public :: output_FcstErrCov
      public :: output_CGSolver
      public :: output_zeit

      public :: output_ON
      public :: output_OFF

      public :: output_fetch
      public :: output_reset

      interface output_fetch; module procedure fetch_;end interface
      interface output_reset; module procedure reset_;end interface

! !REVISION HISTORY:
! 	19Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_output'

  logical,parameter :: output_ON = .true.
  logical,parameter :: output_OFF= .false.

  logical,parameter :: dft_ObsErrCov =output_OFF
  logical,parameter :: dft_FcstErrCov=output_OFF
  logical,parameter :: dft_CGSolver  =output_ON
  logical,parameter :: dft_zeit	     =output_OFF

  logical,save :: output_ObsErrCov =dft_ObsErrCov
  logical,save :: output_FcstErrCov=dft_FcstErrCov
  logical,save :: output_CGSolver  =dft_CGSolver
  logical,save :: output_zeit      =dft_zeit

  character(len=*),parameter :: outputRC_ObsErrCov ='output*ObsErrCov:'
  character(len=*),parameter :: outputRC_FcstErrCov='output*FcstErrCov:'
  character(len=*),parameter :: outputRC_CGSolver  ='output*CGSolver:'
  character(len=*),parameter :: outputRC_zeit      ='output*zeit:'
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: fetch_ - initialize the control module data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine fetch_()
      implicit none

! !REVISION HISTORY:
! 	19Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::fetch_'

  call fetch_setting_(outputRC_ObsErrCov ,output_ObsErrCov )
  call fetch_setting_(outputRC_FcstErrCov,output_FcstErrCov)
  call fetch_setting_(outputRC_CGSolver  ,output_CGSolver  )
  call fetch_setting_(outputRC_zeit      ,output_zeit      )

contains
subroutine fetch_setting_(rsrc,lval)
  use m_inpak90,only : i90_label, i90_gtoken
  use m_chars,  only : uppercase
  use m_stdio,  only : stderr
  use m_die,    only : die
  character(len=*),intent(in) :: rsrc
  logical,intent(inout) :: lval

  integer :: ier
  character(len=32) :: cval

  call i90_label(rsrc,ier)
  if(ier == 0) then
    call i90_gtoken(cval,ier)
    if(ier /= 0) then
      write(stderr,'(4a)') myname_,	&
	': missing value, "',rsrc,'"'
      call die(myname_)
    endif

    select case(uppercase(cval))
    case('YES')
      lval=.true.
    case('NO')
      lval=.false.
    case default
      write(stderr,'(6a)') myname_,	&
	': invalid "',rsrc,'" value, "',trim(cval),'"'
      call die(myname_)
    end select

  endif
end subroutine fetch_setting_
end subroutine fetch_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: reset_ - reset to their default values
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine reset_()
      implicit none

! !REVISION HISTORY:
! 	19Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::reset_'

  output_ObsErrCov =dft_ObsErrCov
  output_FcstErrCov=dft_FcstErrCov
  output_CGSolver  =dft_CGSolver
  output_zeit      =dft_zeit

end subroutine reset_
end module m_output
!.
