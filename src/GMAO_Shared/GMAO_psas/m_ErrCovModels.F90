!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ErrCovModels - PSAS error covariance models
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ErrCovModels
      implicit none
      private	! except

      public :: ErrCovModels_update	! update models

      interface ErrCovModels_update; module procedure	&
	update_
      end interface

! !REVISION HISTORY:
! 	22Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ErrCovModels'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: update_ - Update error covariance models
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine update_(root,comm,rsrc)
      use m_die,only : die
      use m_xRSRC_sigFi,only : fetch_sigFi
      use m_xTab_sigFi, only : mptab_sigFi
      use m_psasrc, only : psasrc_open,psasrc_close
      use m_psasrc, only : psasrc_name
      implicit none
      integer,intent(in) :: root	! root PE
      integer,intent(in) :: comm	! MPI communicator
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	22Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::update_'
  integer :: ier

  call psasrc_open(root,comm,rsrc=rsrc,stat=ier)
	if(ier/=0) call die(myname_,'psasrc_open()',ier)

	! Observation error covariance models

  call ktname0()
  call kxname0()
  call set_OEclas()
  call set_OEhCor()
  call set_OEvCor()

	! Forecast error covariance models

  call set_FEhCor()
  call set_FEvCor()
  call tabl_FEsigW()
  call tabl_FEalpha()
  call fetch_sigFi()

  call psasrc_close(stat=ier)
	if(ier/=0) call die(myname_,'psasrc_close()',ier)

	! tab_sigFi() needs to be parallelized.

  call mptab_sigFi(comm=comm,root=root)

end subroutine update_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize all error covariance models
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_()
      implicit none

! !REVISION HISTORY:
! 	22Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

  ! Actions yet to be defined

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean all error covariance models
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      implicit none

! !REVISION HISTORY:
! 	22Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  ! Actions yet to be defined

end subroutine clean_
end module m_ErrCovModels
