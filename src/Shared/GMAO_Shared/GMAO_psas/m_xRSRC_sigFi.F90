!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_xRSRC_sigFi - a place holder for an extented sigFi resource
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_xRSRC_sigFi
      use m_FEsigFi_tabl, only : FEsigFi_Name
      use m_FEsigFi_tabl, only : clean_sigFi => FEsigFi_clean
      implicit none
      private	! except

      public :: fetch_sigFi	! Fetch a sigFi resource
      public :: clean_sigFi

      interface fetch_sigFi; module procedure fetch_; end interface

! !REVISION HISTORY:
! 	31Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_xRSRC_sigFi'

  character(len=*),parameter :: sigFi_rc='fcst_err_grads_descr_file:'

  ! A pending resource handle change:
  ! character(len=*),parameter :: sigFi_rc='FcstErrStdv_MassVars:'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: fetch_ - Fetch a sigFi resource from the resource database
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine fetch_(stat)
      use m_inpak90, only : i90_label, i90_Gtoken
      use m_stdio, only : stderr
      use m_die,   only : die
      implicit none

      integer, optional,intent(out) :: stat

! !REVISION HISTORY:
! 	17Mar97 - Jing Guo <guo@eramus> - initial prototyping and coding
! 	31Dec98 - Jing Guo <guo@thunder> - restructured after GFname()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::fetch_'

  integer :: ierr

  if(present(stat)) stat=0

	! Locate the resource

  call i90_label(sigFi_rc,ierr)
  if(ierr /= 0) then
    write(stderr,'(2a,i4)') myname_,': i90_label() error, irt =',ierr
    if(.not.present(stat)) call die(myname_)
    stat=ierr
    return
  endif

	! Check the token

  call i90_Gtoken(FEsigFi_Name, ierr)
  if(ierr /= 0) then
    write(stderr,'(2a,i3)') myname_,': i90_Gtoken() error, irt =',ierr
    if(.not.present(stat)) call die(myname_)
    stat=ierr
    return
  endif

end subroutine fetch_
end module m_xRSRC_sigFi
!.
