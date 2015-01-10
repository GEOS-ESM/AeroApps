!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_FcstErrCovModels - A front desk for fcst. err. cov. models.
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_FcstErrCovModels
      use m_sigFi_lookups, only : HGHTE
      use m_sigFi_lookups, only : MIXRE
      use m_sigDCWindErr,  only : WPOTE
      implicit none
      private	! except

      public :: FcstErrCovModels_store
      public :: FcstErrCovModels_remove

      public :: FcstErrCovModels_HGHTE
      public :: FcstErrCovModels_MIXRE
      public :: FcstErrCovModels_WPOTE

      interface FcstErrCovModels_store; module procedure	&
	store3d_,	&
	storePT_
      end interface

      interface FcstErrCovModels_remove; module procedure	&
	remove_
      end interface

! !REVISION HISTORY:
! 	08Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_FcstErrCovModels'

  character(len=*),parameter :: FcstErrCovModels_HGHTE=HGHTE
  character(len=*),parameter :: FcstErrCovModels_MIXRE=MIXRE
  character(len=*),parameter :: FcstErrCovModels_WPOTE=WPOTE

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: store3d_ - store a fcst. err. cov. model component
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine store3d_(varn,nlon,nlat,nlev,levs,sigf,fill,	&
	comm,root,stat)
      use m_sigFi_lookups, only : HGHTE,MIXRE
      use m_sigFi_lookups, only : sigFi_store
      use m_die,  only : perr,die
      use m_stdio,only : stderr
      implicit none

      character(len=*),intent(in) :: varn	! component name

      integer,intent(in) :: nlon	! no. of longitude. grid points.
      integer,intent(in) :: nlat	! no. of latitude. grid points.
      integer,intent(in) :: nlev	! no. of level grid points.
      real,dimension(:),    intent(in) :: levs	! level grid values
      real,dimension(:,:,:),intent(in) :: sigf	! 3-d field values.
      real,optional,intent(in) :: fill	! a "missing" value in the field

      integer,optional,intent(in)  :: comm	! communicator
      integer,optional,intent(in)  :: root	! which is the root PE
      integer,optional,intent(out) :: stat	! status of the call

! !REVISION HISTORY:
! 	08Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::store3d_'
  logical :: comm_
  integer :: ier

  if(present(stat)) stat=0

	! Verify input arguments

  comm_=.false.
  if(present(comm).or.present(root)) then
    if(.not.present(comm) .or. .not.present(root)) then
      if(.not.present(comm)) call perr(myname_,'missing argmument comm')
      if(.not.present(root)) call perr(myname_,'missing argmument root')
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    comm_=.true.
  endif

  select case(varn)
  case (FcstErrCovModels_HGHTE)
    if(comm_) then
      call sigFi_store(comm,root,HGHTE,		&
	nlon,nlat,nlev,levs,sigf,fill=fill,stat=ier)
    else
      call sigFi_store(		 HGHTE,		&
	nlon,nlat,nlev,levs,sigf,fill=fill,stat=ier)
    endif
	if(ier/=0) then
	  call perr(myname_,'sigFi_store("'//HGHTE//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  case (FcstErrCovModels_MIXRE)
    if(comm_) then
      call sigFi_store(comm,root,MIXRE,		&
	nlon,nlat,nlev,levs,sigf,fill=fill,stat=ier)
    else
      call sigFi_store(		 MIXRE,		&
	nlon,nlat,nlev,levs,sigf,fill=fill,stat=ier)
    endif
	if(ier/=0) then
	  call perr(myname_,'sigFi_store("'//MIXRE//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  case default
    write(stderr,'(4a)') myname_,': unknown varn, "',varn,'"'
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return

  end select

end subroutine store3d_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: storePT_ - store a fcst. err. cov. model component
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine storePT_(varn,mtyp,nlev,levs,npar,tabl,desc,	&
	comm,root,stat)
      use m_sigDCWindErr, only : sigDCWindErr_store
      use m_die,  only : perr,die
      use m_stdio,only : stderr
      implicit none
      character(len=*), intent(in) :: varn	! component name

      character(len=*), intent(in) :: mtyp	! model type
      integer,          intent(in) :: nlev	! no. of level grid pt.
      real,dimension(:),intent(in) :: levs	! level grid values
      integer,          intent(in) :: npar	! no. of parameters
      real,dimension(:,:),intent(in) :: tabl	! a parameter table
      character(len=*),optional,intent(in) :: desc	! a short desc.

      integer, optional,intent(in)  :: comm	! communcator
      integer, optional,intent(in)  :: root	! which is the root PE
      integer, optional,intent(out) :: stat	! status of the call

! !REVISION HISTORY:
! 	08Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::storePT_'
  logical :: comm_
  integer :: ier

  if(present(stat)) stat=0

	! Verify input arguments

  comm_=.false.
  if(present(comm).or.present(root)) then
    if(.not.present(comm) .or. .not.present(root)) then
      if(.not.present(comm)) call perr(myname_,'missing argmument comm')
      if(.not.present(root)) call perr(myname_,'missing argmument root')
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    comm_=.true.
  endif

  select case(varn)
  case (FcstErrCovModels_WPOTE)
    if(comm_) then
      call sigDCWindErr_store(comm,root,	&
	mtyp,nlev,levs,npar,tabl,desc=desc,stat=ier)
    else
      call sigDCWindErr_store(			&
	mtyp,nlev,levs,npar,tabl,desc=desc,stat=ier)
    endif
	if(ier/=0) then
	  call perr(myname_,'sigDCWindErr_store()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  case default
    write(stderr,'(4a)') myname_,': unknown varn, "',varn,'"'
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return

  end select

end subroutine storePT_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: remove_ - remove a fcst. err. cov. model component
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine remove_(varn,stat)

      use m_sigFi_lookups, only : HGHTE,MIXRE
      use m_sigFi_lookups, only : sigFi_remove

      use m_sigDCWindErr,  only : sigDCWindErr_remove

      use m_die,  only : perr,die
      use m_stdio,only : stderr

      implicit none
      character(len=*),intent(in) :: varn	! component name
      integer,optional,intent(out):: stat	! status of the call

! !REVISION HISTORY:
! 	08Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::remove_'
  integer :: ier

  if(present(stat)) stat=0

  select case(varn)
  case (FcstErrCovModels_HGHTE)
    call sigFi_remove(HGHTE,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'sigFi_remove("'//HGHTE//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  case (FcstErrCovModels_MIXRE)
    call sigFi_remove(MIXRE,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'sigFi_remove("'//MIXRE//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  case (FcstErrCovModels_WPOTE)
    call sigDCWindErr_remove(stat=ier)
	if(ier/=0) then
	  call perr(myname_,'sigDCWindErr_remove()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  case default
    write(stderr,'(4a)') myname_,': unknown varn, "',varn,'"'
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return

  end select

end subroutine remove_

end module m_FcstErrCovModels
