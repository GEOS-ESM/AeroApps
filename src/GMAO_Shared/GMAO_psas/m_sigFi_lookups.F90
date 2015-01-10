!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_sigFi_lookups - State 2(lookup tables) of mass-coupled sigF
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_sigFi_lookups

      use m_LLAGrid, only : LLAGrid
      use config,    only : dft_FILL => FILL

      implicit none
      private	! except

      public :: sigFi_lookups_init
      public :: sigFi_lookups_clean

      public :: sigFi_store
      public :: sigFi_remove

      public :: HGHTE,hghte_defined,hghte_setbyRC
      public :: MIXRE,mixre_defined,mixre_setbyRC
      public :: hghte_dim,hghte_fld,hghte_fill
      public :: mixre_dim,mixre_fld,mixre_fill

      interface sigFi_lookups_init
	module procedure init_
      end interface

      interface sigFi_lookups_clean; module procedure	&
	clean_all_,	&
	clean__
      end interface

      interface sigFi_store; module procedure	&
	mp_store_,		&
	store_
      end interface

      interface sigFi_remove; module procedure	&
	clean__
      end interface

! !REVISION HISTORY:
! 	09Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_sigFi_lookups'

  character(len=*),parameter :: HGHTE='HGHTE'
  logical,save          :: hghte_defined=.false.
  logical,save		:: hghte_setbyRC=.false.
  type(LLAGrid),save    :: hghte_dim
  real,save,allocatable :: hghte_fld(:,:,:)
  real,save		:: hghte_fill=dft_FILL

  character(len=*),parameter :: MIXRE='MIXRE'
  logical,save          :: mixre_defined=.false.
  logical,save		:: mixre_setbyRC=.false.
  type(LLAGrid),save    :: mixre_dim
  real,save,allocatable :: mixre_fld(:,:,:)
  real,save		:: mixre_fill=dft_FILL

  integer,parameter :: LEN_VARN=max(len(HGHTE),len(MIXRE))

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: mp_store_ - initialized on all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine mp_store_(comm,root,varn,nlon,nlat,nlev,levs,	&
	sigf,fill,stat)
      use m_die,  only : die,perr
      use m_mpif90,only : MP_perr
      use m_mpif90,only : MP_type
      use m_mpif90,only : MP_comm_rank
      implicit none
      integer,intent(in) :: comm
      integer,intent(in) :: root

      character(len=*),intent(in) :: varn
      integer,intent(in) :: nlon,nlat,nlev
      real,intent(in) :: levs(:)
      real,dimension(:,:,:),intent(in) :: sigf
      real,   optional,intent(in)  :: fill
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	06Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::mp_store_'
  character(len=LEN_VARN) :: varn_
  real :: fill_
  integer :: nlon_,nlat_,nlev_
  integer,dimension(3) :: ibufr
  real,allocatable,dimension(:) :: levs_
  integer :: myID
  integer :: ier

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif

  if(myID==root) then
			! Arguments of sigF are only significant on
			! the root PE.
    call check_varn_(myname_,varn,ier,stat=stat)
	if(ier/=0) return

    call check_3d_(myname_,nlon,nlat,nlev,sigf,ier,stat=stat)
	if(ier/=0) return

    call check_1d_(myname_,nlev,levs,ier,stat=stat)
	if(ier/=0) return

  endif
!________________________________________

  fill_=dft_FILL
  if(present(fill)) fill_=fill

  varn_=varn	! if varn is qualified on the root PE, the size of
		! varn_ is also qualified.

  nlon_=nlon
  nlat_=nlat
  nlev_=nlev
!________________________________________

  call MPI_bcast(varn_,LEN_VARN,MP_type(varn_),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(varn)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif
	!________________________________

	ibufr(1)=nlon_
	ibufr(2)=nlat_
	ibufr(3)=nlev_

  call MPI_bcast(ibufr,size(ibufr),MP_type(ibufr),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(ibufr)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif

	nlon_=ibufr(1)
	nlat_=ibufr(2)
	nlev_=ibufr(3)
	!________________________________

  call MPI_bcast(fill_,1,MP_type(fill_),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(fill)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif
	!________________________________

	allocate(levs_(nlev_),stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allocate(levs_)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=-1
		  return
		endif

	if(myID==root) levs_(1:nlev_)=levs(1:nlev)

  call MPI_bcast(levs_,nlev_,MP_type(levs_),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast(levs)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif
	!________________________________

		! Allocate/define local storages

  call init_(varn_,nlon_,nlat_,nlev_,levs_,fill=fill_,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'init_("'//trim(varn_)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	deallocate(levs_,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'deallocate(levs_)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=-1
		  return
		endif
!________________________________________

  select case(varn_)
  case (HGHTE)
  	if(myID==root) hghte_fld(:,:,:)=sigf(:,:,:)

    call MPI_bcast(hghte_fld,size(hghte_fld),	&
	MP_type(hghte_fld),root,comm,ier)
		if(ier/=0) then
		  call MP_perr(myname_,'MPI_bcast(hghte_fld)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=-1
		  return
		endif

    hghte_defined=.true.
    hghte_setbyRC=.false.

  case (MIXRE)
  	if(myID==root) mixre_fld(:,:,:)=sigf(:,:,:)

    call MPI_bcast(mixre_fld,size(mixre_fld),	&
	MP_type(mixre_fld),root,comm,ier)
		if(ier/=0) then
		  call MP_perr(myname_,'MPI_bcast(mixre_fld)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=-1
		  return
		endif

    mixre_defined=.true.
    mixre_setbyRC=.false.

  case default
    call perr(myname_,'unknown varn, "'//trim(varn_)//'"')
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  end select

end subroutine mp_store_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: check_varn_ - check varn
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine check_varn_(where,varn,ier,stat)
      use m_die,only : perr,die
      implicit none
      character(len=*),intent(in) :: where
      character(len=*),intent(in) :: varn
      integer,intent(out) :: ier
      integer,optional,intent(inout) :: stat

! !REVISION HISTORY:
! 	27Dec00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::check_varn_'

  ier=0

  select case(varn)
  case (HGHTE)
    if(hghte_defined) then
      call perr(where,'multiple "'//trim(varn)//'" definition')
      if(.not.present(stat)) call die(where)
      stat=-1
      ier=-1
      return
    endif

  case (MIXRE)
    if(mixre_defined) then
      call perr(where,'multiple "'//trim(varn)//'" definition')
      if(.not.present(stat)) call die(where)
      stat=-1
      ier=-1
      return
    endif

  case default
    call perr(where,'unknown varn, "'//trim(varn)//'"')
    if(.not.present(stat)) call die(where)
    stat=-1
    ier=-1
    return
  end select

end subroutine check_varn_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: check_3d_ - check sigf(nlon,nlat,nlev)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine check_3d_(where,nlon,nlat,nlev,sigf,ier,stat)
      use m_die,only : perr,die
      implicit none
      character(len=*),intent(in) :: where
      integer,intent(in) :: nlon
      integer,intent(in) :: nlat
      integer,intent(in) :: nlev
      real,dimension(:,:,:),intent(in) :: sigf
      integer,intent(out) :: ier
      integer,optional,intent(inout) :: stat

! !REVISION HISTORY:
! 	27Dec00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::check_3d_'

  ier=0

  if(	nlon/=size(sigf,1)	.or.	&
	nlat/=size(sigf,2)	.or.	&
	nlev/=size(sigf,3)	) then

    call perr(where,'inconsistent size(sigf) error.')

    if(	nlon/=size(sigf,1) ) call perr(		&
	where,': nlon=',nlon,', size(sigf,1) =',size(sigf,1))
    if(	nlat/=size(sigf,2) ) call perr(		&
	where,': nlat=',nlat,', size(sigf,2) =',size(sigf,2))
    if(	nlev/=size(sigf,3) ) call perr(		&
	where,': nlev=',nlev,', size(sigf,3) =',size(sigf,3))

    if(.not.present(stat)) call die(where)
    stat=-1
    ier=-1
    return
  endif
end subroutine check_3d_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: check_1d_ - check a 1-d array dimension
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine check_1d_(where,nlev,levs,ier,stat)
      use m_die,only : perr,die
      implicit none
      character(len=*),intent(in) :: where
      integer,intent(in) :: nlev
      real,dimension(:),intent(in) :: levs
      integer,intent(out) :: ier
      integer,optional,intent(inout) :: stat

! !REVISION HISTORY:
! 	27Dec00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::check_1d_'

  ier=0
  if(nlev/=size(levs,1) ) then
    call perr(where,'inconsistent size(levs) error.')
    call perr(where,': nlev=',nlev,', size(levs) =',size(levs))
    if(.not.present(stat)) call die(where)
    stat=-1
    ier=-1
    return
  endif
end subroutine check_1d_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: store_ - initialized both the datastructure and the field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine store_(varn,nlon,nlat,nlev,levs,sigf,fill,stat)
      use m_die,  only : die,perr
      implicit none
      character(len=*),intent(in) :: varn
      integer,intent(in) :: nlon,nlat,nlev
      real,intent(in) :: levs(:)
      real,dimension(:,:,:),intent(in) :: sigf
      real,   optional,intent(in)  :: fill
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	06Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::store_'
  integer :: ier

  if(present(stat)) stat=0

  call check_varn_(myname_,varn,ier,stat=stat)
	if(ier/=0) return
  call check_3d_(myname_,nlon,nlat,nlev,sigf,ier,stat=stat)
	if(ier/=0) return
  call check_1d_(myname_,nlev,levs,ier,stat=stat)
	if(ier/=0) return

  call init_(varn,nlon,nlat,nlev,levs,fill=fill,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'init_("'//trim(varn)//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  select case(varn)
  case (HGHTE)
    hghte_fld(:,:,:)=sigf(:,:,:)
    hghte_defined=.true.
    hghte_setbyRC=.false.

  case (MIXRE)
    mixre_fld(:,:,:)=sigf(:,:,:)
    mixre_defined=.true.
    mixre_setbyRC=.false.

  case default
    call perr(myname_,'unknown varn, "'//trim(varn)//'"')
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  end select

end subroutine store_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(varn,nlon,nlat,nlev,levs,fill,stat)
      use m_LLAGrid, only : LLAGrid_init
      use m_die,     only : die,perr
      use m_mall,    only : mall_mci,mall_ison
      implicit none
      character(len=*),intent(in) :: varn
      integer,intent(in) :: nlon,nlat,nlev
      real,intent(in) :: levs(:)
      real,   optional,intent(in)  :: fill
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	09Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier

  if(present(stat)) stat = 0
  select case(varn)
  case (HGHTE)
    if(hghte_defined) then
      call perr(myname_,'multiple "'//trim(varn)//'" definition',ier)
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    call LLAGrid_init(hghte_dim,nlon,nlat,nlev,levs,stat=ier)
    if(ier/=0) then
      call perr(myname_,'LLAGrid_init("'//trim(varn)//'")',ier)
      if(.not.present(stat)) call die(myname_)
      stat=ier
      return
    endif
    allocate(hghte_fld(nlon,nlat,nlev),stat=ier)
    if(ier/=0) then
      call perr(myname_,'allocate("'//trim(varn)//'")',ier)
      if(.not.present(stat)) call die(myname_)
      stat=ier
      return
    endif

	if(mall_ison()) call mall_mci(hghte_fld,myname)

    hghte_defined=.true.
    hghte_fill   =dft_FILL
    if(present(fill)) hghte_fill=fill

  case (MIXRE)
    if(mixre_defined) then
      call perr(myname_,'multiple "'//trim(varn)//'" definition',ier)
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    call LLAGrid_init(mixre_dim,nlon,nlat,nlev,levs,stat=ier)
    if(ier/=0) then
      call perr(myname_,'LLAGrid_init("'//trim(varn)//'")',ier)
      if(.not.present(stat)) call die(myname_)
      stat=ier
      return
    endif
    allocate(mixre_fld(nlon,nlat,nlev),stat=ier)
    if(ier/=0) then
      call perr(myname_,'allocate("'//trim(varn)//'")',ier)
      if(.not.present(stat)) call die(myname_)
      stat=ier
      return
    endif

	if(mall_ison()) call mall_mci(mixre_fld,myname)

    mixre_defined=.true.
    mixre_fill   =dft_FILL
    if(present(fill)) mixre_fill=fill

  case default
    call perr(myname_,'unknown varn, "'//trim(varn)//'"')
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  end select

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: clean_all_ - clean all sigFi_lookups
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_all_(stat)
      use m_die,     only : die,perr
      implicit none
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	10Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_all_'
  integer :: ier

  if(present(stat)) stat=0

  call clean__(HGHTE,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'clean__("'//HGHTE//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  call clean__(MIXRE,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'clean__("'//MIXRE//'")',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine clean_all_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean__ - clean a sigFi_lookups
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean__(varn,stat)
      use m_LLAGrid, only : LLAGrid_clean
      use m_die,     only : die,perr
      use m_mall,    only : mall_mco,mall_ison
      implicit none
      character(len=*),intent(in) :: varn
      integer,optional,intent(out):: stat

! !REVISION HISTORY:
! 	09Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean__'
  integer :: ier

  if(present(stat)) stat=0

  select case(varn)
  case(HGHTE)
    if(.not.hghte_defined) then
      call perr(myname_,'undefined object, "'//trim(varn)//'"')
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    call LLAGrid_clean(hghte_dim,stat=ier)
    if(ier /= 0) then
      call perr(myname_,'LLAGrid_clean("'//trim(varn)//'"_dim)',ier)
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

	if(mall_ison()) call mall_mco(hghte_fld,myname)

    deallocate(hghte_fld,stat=ier)
    if(ier/=0) then
      call perr(myname_,'deallocate("'//trim(varn)//'"_fld)',ier)
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    hghte_defined=.false.
    hghte_setbyRC=.false.
    hghte_fill	 =dft_FILL

  case(MIXRE)
    if(.not.mixre_defined) then
      call perr(myname_,'undefined object, "'//trim(varn)//'"')
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    call LLAGrid_clean(mixre_dim,stat=ier)
    if(ier /= 0) then
      call perr(myname_,'LLAGrid_clean("'//trim(varn)//'"_dim)',ier)
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

	if(mall_ison()) call mall_mco(mixre_fld,myname)

    deallocate(mixre_fld,stat=ier)
    if(ier/=0) then
      call perr(myname_,'deallocate("'//trim(varn)//'"_fld)',ier)
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    mixre_defined=.false.
    mixre_setbyRC=.false.
    mixre_fill	 =dft_FILL

  case default
    call perr(myname_,'unknown varn, "'//trim(varn)//'"')
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  end select

end subroutine clean__
end module m_sigFi_lookups
