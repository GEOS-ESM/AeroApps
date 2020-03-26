!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_xTab_sigFi - Tabulate sigFi
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_xTab_sigFi
      implicit none
      private	! except

      public :: tab_sigFi	! tabulate sigFi from a resource
      public :: mptab_sigFi	! tabulate distributed sigFi
      public :: rmtab_sigFi	! remove the tabulated sigFi

      interface mptab_sigFi;   module procedure	&
	mptab_sigFi_
      end interface
      interface tab_sigFi;   module procedure	&
	tab_sigFi_
      end interface
      interface rmtab_sigFi; module procedure	&
	rmtab_sigFi_
      end interface

! !REVISION HISTORY:
! 	31Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_xTab_sigFi'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: tab_sigFi_ - setting up the (m_sigFi_lookups) module
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine tab_sigFi_(tick,stat)

      use m_GrADSfiles, only : GrADS_open
      use m_GrADSfiles, only : GrADS_close
      use m_GrADSfiles, only : GrADSfiles

      use m_mpout,      only : mpout_log
      use m_die,        only : die,perr

      use m_FEsigFi_tabl,  only : FEsigFi_Name
      use m_sigFi_lookups, only : HGHTE,hghte_defined,hghte_setbyRC
      use m_sigFi_lookups, only : MIXRE,mixre_defined,mixre_setbyRC

      implicit none
      integer,optional,intent(in)  :: tick
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	09Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!	06Oct99	- Jing Guo
!		. Initialized stat, if the argument presents
!		. Allowed initializing HGHTE or MIXRE at a different
!		  place.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::tab_sigFi_'
  type(GrADSfiles) :: gs
  integer :: ier
  integer :: llev

  if(present(stat)) stat=0

  if( hghte_defined.and.hghte_setbyRC	.or.	&
      mixre_defined.and.mixre_setbyRC	)	then

    if( hghte_defined.and.hghte_setbyRC )	&
	call perr(myname_,'multiple "'//HGHTE//'" definition')
    if( mixre_defined.and.mixre_setbyRC )	&
	call perr(myname_,'multiple "'//MIXRE//'" definition')

    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

	! If both HGHTE and MIXRE are defined (not _setbyRC), there
	! is not need to contidue

  if(hghte_defined .and. mixre_defined) return

  llev=1
  if(present(tick)) llev=tick

  call mpout_log(myname_,'open "'//trim(FEsigFi_Name)//'" for input')

  call GrADS_open(gs,FEsigFi_Name,stat=ier)
	if(ier/=0) then
	  call perr(myname_,	&
		'GrADS_open("'//trim(FEsigFi_Name)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  if(.not.hghte_defined) then
	! If HGHTE has not been defined (not _setbyRC), read from
	! the GrADS file.

    call get_sigFi_var_(gs,HGHTE,tick=llev,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'get_sigFi_var_("'//HGHTE//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

    hghte_defined=.true.
    hghte_setbyRC=.true.
  endif

  if(.not.mixre_defined) then
	! If MIXRE has not been defined (not _setbyRC), read from
	! the GrADS file.

    call get_sigFi_var_(gs,MIXRE,tick=llev,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'get_sigFi_var_("'//MIXRE//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

    mixre_defined=.true.
    mixre_setbyRC=.true.
  endif

  call GrADS_close(gs,stat=ier)
	if(ier/=0) then
	  call perr(myname_,	&
		'GrADS_close("'//trim(FEsigFi_Name)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine tab_sigFi_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: mptab_sigFi_ - setting up the (m_sigFi_lookups) module
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine mptab_sigFi_(comm,root,tick,stat)

      use m_die,        only : die,perr
      use m_mpif90,     only : MP_comm_rank
      use m_mpif90,     only : MP_type,MP_perr

      use m_sigFi_lookups, only : HGHTE,hghte_defined,hghte_setbyRC
      use m_sigFi_lookups, only : MIXRE,mixre_defined,mixre_setbyRC

      implicit none
      integer,intent(in)  :: comm
      integer,intent(in)  :: root
      integer,optional,intent(in)  :: tick
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	09Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!	06Oct99	- Jing Guo
!		. Initialized stat, if the argument presents
!		. Allowed initializing HGHTE or MIXRE at a different
!		  place.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::mptab_sigFi_'
  integer :: ier
  integer :: myID
  logical,dimension(2) :: lbufr

  if(present(stat)) stat=0

	! Before proceed, check the status on all PEs.

  if( hghte_defined.and.hghte_setbyRC	.or.	&
      mixre_defined.and.mixre_setbyRC	)	then

    if( hghte_defined.and.hghte_setbyRC )	&
	call perr(myname_,'multiple "'//HGHTE//'" definition')
    if( mixre_defined.and.mixre_setbyRC )	&
	call perr(myname_,'multiple "'//MIXRE//'" definition')

    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  if(myID==root) then
    call tab_sigFi_(tick=tick,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'tab_sigFi_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  endif

	! Broadcast the information before take another step.

  lbufr(1)=hghte_setbyRC
  lbufr(2)=mixre_setbyRC
  call MPI_bcast(lbufr,size(lbufr),MP_type(lbufr),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(lbufr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  hghte_setbyRC=lbufr(1)
  mixre_setbyRC=lbufr(2)

  if(hghte_setbyRC) call bcast_(HGHTE,myID,root,comm)
  if(mixre_setbyRC) call bcast_(MIXRE,myID,root,comm)

end subroutine mptab_sigFi_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcast_ - broadcast a sigFi field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine bcast_(varn,myID,root,comm)
      use m_sigFi_lookups, only : HGHTE,hghte_defined,hghte_setbyRC
      use m_sigFi_lookups, only : hghte_dim,hghte_fill,hghte_fld
      use m_sigFi_lookups, only : MIXRE,mixre_defined,mixre_setbyRC
      use m_sigFi_lookups, only : mixre_dim,mixre_fill,mixre_fld
      use m_sigFi_lookups, only : sigFi_lookups_init

      use m_die   ,only : MP_die,die
      use m_mpif90,only : MP_type
      implicit none
      character(len=*),intent(in) :: varn
      integer,intent(in) :: myID,root,comm

! !REVISION HISTORY:
! 	25Jan01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bcast_'
  integer,dimension(3) :: ibufr
  real,allocatable,dimension(:) :: rbufr
  integer :: nlon,nlat,nlev
  integer :: ier

  select case(varn)
  case(HGHTE)
    ibufr(1)=hghte_dim%nlon
    ibufr(2)=hghte_dim%nlat
    ibufr(3)=hghte_dim%nlev
  case(MIXRE)
    ibufr(1)=mixre_dim%nlon
    ibufr(2)=mixre_dim%nlat
    ibufr(3)=mixre_dim%nlev
  case default
    call die(myname_,'unknown varn, "'//varn//'"')
  end select

  call MPI_bcast(ibufr,size(ibufr),MP_type(ibufr),root,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_bcast(ibufr)',ier)

  nlon=ibufr(1)
  nlat=ibufr(2)
  nlev=ibufr(3)

  allocate(rbufr(0:nlev),stat=ier)
	if(ier/=0) call die(myname_,'allocate(rbufr)',ier)

  if(myID==root) then
    select case(varn)
    case(HGHTE)
      rbufr(1:nlev)=hghte_dim%vlev(1:nlev)
      rbufr(0)=hghte_fill
    case(MIXRE)
      rbufr(1:nlev)=mixre_dim%vlev(1:nlev)
      rbufr(0)=mixre_fill
    end select
  endif

  call MPI_bcast(rbufr,size(rbufr),MP_type(rbufr),root,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_bcast(rbufr)',ier)

  if(myID/=root) then
    call sigFi_lookups_init(varn,	&
	nlon,nlat,nlev,rbufr(1:nlev),fill=rbufr(0),stat=ier)

	if(ier/=0) call die(myname_,	&
		'sigFi_lookups_init("'//varn//'")',ier)
  endif

  deallocate(rbufr,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(rbufr)',ier)

  select case(varn)
  case(HGHTE)
    hghte_defined=.true.
    call MPI_bcast(hghte_fld,size(hghte_fld),	&
	MP_type(hghte_fld),root,comm,ier)
  case(MIXRE)
    mixre_defined=.true.
    call MPI_bcast(mixre_fld,size(mixre_fld),	&
	MP_type(mixre_fld),root,comm,ier)
  end select

	if(ier/=0) call MP_die(myname_,'MPI_bcast("'//varn//'")',ier)

end subroutine bcast_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_sigFi_var_ - setting up a sigFi lookup table
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_sigFi_var_(gs,varn,tick,stat)
      use m_GrADSfiles,    only : GrADS_getdims
      use m_GrADSfiles,    only : GrADS_zdef
      use m_GrADSfiles,    only : GrADS_input
      use m_GrADSfiles,    only : GrADSfiles
      use m_sigFi_lookups, only : sigFi_lookups_init
      use m_sigFi_lookups, only : HGHTE,MIXRE
      use m_sigFi_lookups, only : hghte_fld,mixre_fld
      use m_die,           only : die,perr

      implicit none

      type(GrADSfiles),intent(inout) :: gs
      character(len=*),intent(in)  :: varn
      integer,optional,intent(in)  :: tick
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	09Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!	06Oct99	- Jing Guo
!		. Initialized stat, if the argument presents
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_sigFi_var_'
  integer :: nlon,nlat,nlev,ier
  integer :: llev
  real    :: fill

  if(present(stat)) stat=0

  llev=1
  if(present(tick)) llev=tick

  call GrADS_getdims(gs,varn,nlon,nlat,nlev,udef=fill,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'GrADS_getdims("'//trim(varn)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  call sigFi_lookups_init(varn,nlon,nlat,nlev,	&
    GrADS_zdef(gs,nlev),fill=fill,stat=ier)
	if(ier/=0) then
	  call perr(myname_,	&
		'sigFi_lookups_init("'//trim(varn)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  select case(varn)
  case(HGHTE)
    call GrADS_input(gs,varn,llev,hghte_fld,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'GrADS_input("'//trim(varn)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  case(MIXRE)
    call GrADS_input(gs,varn,llev,mixre_fld,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'GrADS_input("'//trim(varn)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  case default
	call perr(myname_,'unknown tag, "'//trim(varn)//'"')
	if(.not.present(stat)) call die(myname_)
	stat=-1
	return
  end select

end subroutine get_sigFi_var_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rmtab_sigFi_ - remove the table defined by tab_
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine rmtab_sigFi_(stat)
      use m_sigFi_lookups,only : sigFi_lookups_clean
      use m_sigFi_lookups,only : HGHTE,hghte_defined,hghte_setbyRC
      use m_sigFi_lookups,only : MIXRE,mixre_defined,mixre_setbyRC
      use m_die,  only : die,perr
      implicit none
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	06Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rmtab_sigFi_'
  integer :: ier

  if(present(stat)) stat=0

	! If HGHTE or MIXRE has not been defined, this call
	! has not been paired by tab_sigFi()_, which should
	! ensure both HGHTE and MIXRE to be defined.

  if(.not.(hghte_defined.and.mixre_defined)) then
    call perr(myname_,'undefined object')
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

  if(hghte_defined.and.hghte_setbyRC) then

	! If HGHTE was defined by tab_sigFi_(), ...
  
    call sigFi_lookups_clean(HGHTE,stat=ier)
	if(ier /= 0) then
	  call perr(myname_,'sigFi_lookups_clean("'//HGHTE//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

    hghte_defined=.false.
    hghte_setbyRC=.false.

  endif

  if(mixre_defined.and.mixre_setbyRC) then
  
	! If MIXRE was defined by tab_sigFi_(), ...

    call sigFi_lookups_clean(MIXRE,stat=ier)
	if(ier /= 0) then
	  call perr(myname_,'sigFi_lookups_clean("'//MIXRE//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

    mixre_defined=.false.
    mixre_setbyRC=.false.

  endif

end subroutine rmtab_sigFi_
end module m_xTab_sigFi
!.
