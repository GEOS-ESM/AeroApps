!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AISolver - MPI PSAS solver interfaces
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_AISolver
      use m_ktList,only : ktslp,ktus,ktvs,ktHH,ktUU,ktVV,ktQQ

      implicit none
      private	! except

      public :: solve
      public :: multiSolve

      interface solve; module procedure	&
	solve_mix_,	&
	solve_puv_,	&
	solve_zuv_,	&
	solve_all_,	&
	solve_vec_
      end interface

      interface multiSolve; module procedure	&
	solve_multi_
      end interface

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AISolver'

  integer,dimension(1),parameter :: KTmix=(/ktQQ/)
  integer,dimension(3),parameter :: KTzuv=(/ktUU,ktVV,ktHH/)
  integer,dimension(3),parameter :: KTpuv=(/ktUs,ktVs,ktslp/)
  integer,dimension(7),parameter :: KTall=	&
	(/ktUs,ktVs,ktslp,ktUU,ktVV,ktHH,ktQQ/)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_mix_ - A univariate upper-air ana. incr. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_mix_( im,jnp,mlev,plevs,qAinc,		&
	nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,dels,		&
	comm,root,psas_rsrc,eagrid_rsrc)

      use m_rootedAIGrid,only : rootedAIGrid
      use m_rootedAIGrid,only : rootedAIGrid_init
      use m_rootedAIGrid,only : distrSize
      use m_rootedAIGrid,only : ptr_lat
      use m_rootedAIGrid,only : ptr_lon
      use m_rootedAIGrid,only : ptr_lev
      use m_rootedAIGrid,only : ptr_kt
      use m_rootedAIGrid,only : gatherv
      use m_rootedAIGrid,only : clean
      use m_rootedAIGrid,only : rootedAIGrid_intp

      use m_die,only : die,perr,MP_die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : gdstat
      use m_AE,    only : AE_solve

      implicit none

      integer,intent(in) :: im		! no. of longitude grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:,:),intent(out) :: qAinc ! gridded analysis

      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(out):: xvec	! solution vector "x"
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU
      real   ,dimension(:),intent(in) :: dels	! an O-F vector

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_mix_'

	! The workspace for a local copy of A-F vector and its
	! attributes.

  type(rootedAIGrid) :: aig

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: inc
  integer :: ier
!________________________________________

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  call verifyObsVec_(nobs,xvec,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,dels,myname_)
		!________________________________
  invalid=.false.
  if(myID==root) then
    if(.not.invalid) invalid = size(plevs)  < mlev
    if(.not.invalid) invalid = size(qAinc,1) < im
    if(.not.invalid) invalid = size(qAinc,2) < jnp
    if(.not.invalid) invalid = size(qAinc,3) < mlev

    if(invalid) then
      if(size(plevs)  < mlev) call perr(myname_,	&
	'mlev',mlev,'size(plevs)',size(plevs)	)

      if(size(qAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(qAinc,1)',size(qAinc,1)	)
      if(size(qAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(qAinc,2)',size(qAinc,2)	)
      if(size(qAinc,3) < mlev) call perr(myname_,		&
	'mlev',mlev,'size(qAinc,3)',size(qAinc,3)	)

      call die(myname_,'invalid argument(s)')
    endif

  endif
!________________________________________

	! Create a distributed AIGrid vector based on its specification
	! on the root PE.

  if(present(eagrid_rsrc)) then
    call rootedAIGrid_init(aig,im,jnp,mlev,plevs,KTmix,	&
	root,comm,rsrc=eagrid_rsrc)
  else
    call rootedAIGrid_init(aig,im,jnp,mlev,plevs,KTmix,	&
	root,comm,rsrc=psas_rsrc)
  endif

	ninc=distrSize(aig)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate(vinc)',ier)
		if(mall_ison()) call mall_mci(vinc,myname_)

  call AE_solve(ninc,vinc,		&
	ptr_lat(aig),ptr_lon(aig),ptr_lev(aig),ptr_kt(aig),	&
	nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,dels,		&
	comm,root,rsrc=psas_rsrc)

		! Cache distributed vinc on the root PE

  call gatherv(vinc,aig)

		! Note that qAinc is only defined on the root PE.

  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktQQ,qAinc)

		if(mall_ison()) call mall_mco(vinc,myname_)
	deallocate(vinc, stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vinc)',ier)

  call clean(aig)

  if(myID==root.and.mpout_ison()) then
    inc=1
    if(mlev > 0) then
      if(plevs(1)>plevs(mlev)) inc=-1
    endif

    call gdstat(mpout,im,jnp,mlev,qAinc,plevs,'MIXR','PRES',1.e+15,  &
	'Upper-air water vapor mixing-ratio analysis increments',inc)
  endif

end subroutine solve_mix_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_zuv_ - A multivariate upper-air analysis eqn. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_zuv_(im,jnp,mlev,plevs,zAinc,uAinc,vAinc,	&
	nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,dels,		&
	comm,root,psas_rsrc,eagrid_rsrc)

      use m_rootedAIGrid,only : rootedAIGrid
      use m_rootedAIGrid,only : rootedAIGrid_init
      use m_rootedAIGrid,only : distrSize
      use m_rootedAIGrid,only : ptr_lat
      use m_rootedAIGrid,only : ptr_lon
      use m_rootedAIGrid,only : ptr_lev
      use m_rootedAIGrid,only : ptr_kt
      use m_rootedAIGrid,only : gatherv
      use m_rootedAIGrid,only : clean
      use m_rootedAIGrid,only : rootedAIGrid_intp

      use m_die,only : die,perr,MP_die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : gdstat
      use m_AE,    only : AE_solve

      use m_zeit

      implicit none

      integer,intent(in) :: im		! no. of longitude grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:,:),intent(out) :: zAinc ! Gridded h anaysis
      real,dimension(:,:,:),intent(out) :: uAinc ! gridded u analysis
      real,dimension(:,:,:),intent(out) :: vAinc ! gridded v analysis

      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(out):: xvec	! solution vector "x"
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU
      real   ,dimension(:),intent(in) :: dels	! an O-F vector

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_zuv_'

	! The workspace for a local copy of A-F vector and its
	! attributes.

  type(rootedAIGrid) :: aig

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: inc
  integer :: ier

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  call verifyObsVec_(nobs,xvec,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,dels,myname_)

  if(myID==root) then
    invalid=.false.
    if(.not.invalid) invalid = size(plevs)   < mlev

    if(.not.invalid) invalid = size(zAinc,1) < im
    if(.not.invalid) invalid = size(zAinc,2) < jnp
    if(.not.invalid) invalid = size(zAinc,3) < mlev

    if(.not.invalid) invalid = size(uAinc,1) < im
    if(.not.invalid) invalid = size(uAinc,2) < jnp
    if(.not.invalid) invalid = size(uAinc,3) < mlev

    if(.not.invalid) invalid = size(vAinc,1) < im
    if(.not.invalid) invalid = size(vAinc,2) < jnp
    if(.not.invalid) invalid = size(vAinc,3) < mlev

    if(invalid) then
      if(size(plevs)  < mlev) call perr(myname_,	&
		'mlev',mlev,'size(plevs)',size(plevs)	)

      if(size(zAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(zAinc,1)',size(zAinc,1)	)
      if(size(zAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(zAinc,2)',size(zAinc,2)	)
      if(size(zAinc,3) < mlev) call perr(myname_,	&
	'mlev',mlev,'size(zAinc,3)',size(zAinc,3)	)

      if(size(uAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(uAinc,1)',size(uAinc,1)	)
      if(size(uAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(uAinc,2)',size(uAinc,2)	)
      if(size(uAinc,3) < mlev) call perr(myname_,	&
	'mlev',mlev,'size(uAinc,3)',size(uAinc,3)	)

      if(size(vAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(vAinc,1)',size(vAinc,1)	)
      if(size(vAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(vAinc,2)',size(vAinc,2)	)
      if(size(vAinc,3) < mlev) call perr(myname_,	&
	'mlev',mlev,'size(vAinc,3)',size(vAinc,3)	)

      call die(myname_,'invalid argument(s)')
    endif

  endif
!________________________________________

	! Create a distributed AIGrid vector based on its specification
	! on the root PE.

	call zeit_ci('rootedAIGrid')

  if(present(eagrid_rsrc)) then
    call rootedAIGrid_init(aig,im,jnp,mlev,plevs,KTzuv,	&
	root,comm,rsrc=eagrid_rsrc)
  else
    call rootedAIGrid_init(aig,im,jnp,mlev,plevs,KTzuv,	&
	root,comm,rsrc=psas_rsrc)
  endif
	call zeit_co('rootedAIGrid')

	ninc=distrSize(aig)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate(vinc)',ier)
		if(mall_ison()) call mall_mci(vinc,myname_)

	call zeit_ci('AE_solve')

  call AE_solve(ninc,vinc,		&
	ptr_lat(aig),ptr_lon(aig),ptr_lev(aig),ptr_kt(aig),	&
	nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,dels,		&
	comm,root,rsrc=psas_rsrc)
	call zeit_co('AE_solve')

		! Cache distributed vinc on the root PE.

	call zeit_ci('rootedAIGrid_intp')
  call gatherv(vinc,aig)

		! Note that h-u-v ainc. is only defined on root PE.

  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktHH,     zAinc)
  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktUU,ktVV,uAinc,vAinc)
	call zeit_co('rootedAIGrid_intp')

		if(mall_ison()) call mall_mco(vinc,myname_)
	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vinc)',ier)

  call clean(aig)

  if(myID==root.and.mpout_ison()) then
    inc=1
    if(mlev > 0) then
      if(plevs(1)>plevs(mlev)) inc=-1
    endif

    call gdstat(mpout,im,jnp,mlev,uAinc,plevs,'WIND','PRES',	&
	1.e+15,'Upper-air u-wind analysis increments',inc	)
    call gdstat(mpout,im,jnp,mlev,vAinc,plevs,'WIND','PRES',	&
	1.e+15,'Upper-air v-wind analysis increments',inc	)
    call gdstat(mpout,im,jnp,mlev,zAinc,plevs,'HGHT','PRES',	&
	1.e+15,'Upper-air height analysis increments',inc	)
  endif

end subroutine solve_zuv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_puv_ - A multivariate sea-level analysis eqn. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_puv_( im,jnp,pAinc,uAinc,vAinc,		&
	nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,dels,		&
	comm,root,psas_rsrc,eagrid_rsrc)

      use m_rootedAIGrid,only : rootedAIGrid
      use m_rootedAIGrid,only : rootedAIGrid_init
      use m_rootedAIGrid,only : distrSize
      use m_rootedAIGrid,only : ptr_lat
      use m_rootedAIGrid,only : ptr_lon
      use m_rootedAIGrid,only : ptr_lev
      use m_rootedAIGrid,only : ptr_kt
      use m_rootedAIGrid,only : gatherv
      use m_rootedAIGrid,only : clean
      use m_rootedAIGrid,only : rootedAIGrid_intp

      use m_die,only : die,perr,MP_die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : lvstat
      use m_AE,    only : AE_solve

      implicit none

      integer,intent(in) :: im		! no. of longitude grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      real,dimension(:,:),intent(out) :: pAinc	! gridded p_sl analysis
      real,dimension(:,:),intent(out) :: uAinc	! gridded u_sl analysis
      real,dimension(:,:),intent(out) :: vAinc	! gridded v_sl analysis

      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(out):: xvec	! solution vector "x"
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU
      real   ,dimension(:),intent(in) :: dels	! an O-F vector

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_puv_'

	! The workspace for a local copy of the O-F vector and its
	! attributes.

  type(rootedAIGrid) :: aig

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: ier
!________________________________________

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  call verifyObsVec_(nobs,xvec,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,dels,myname_)

  if(myID==root) then
    invalid=.false.
    if(.not.invalid) invalid = size(pAinc,1) < im
    if(.not.invalid) invalid = size(pAinc,2) < jnp
    if(.not.invalid) invalid = size(uAinc,1) < im
    if(.not.invalid) invalid = size(uAinc,2) < jnp
    if(.not.invalid) invalid = size(vAinc,1) < im
    if(.not.invalid) invalid = size(vAinc,2) < jnp

    if(invalid) then
      if(size(pAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(pAinc,1)',size(pAinc,1)	)
      if(size(pAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(pAinc,2)',size(pAinc,2)	)

      if(size(uAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(uAinc,1)',size(uAinc,1)	)
      if(size(uAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(uAinc,2)',size(uAinc,2)	)

      if(size(vAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(vAinc,1)',size(vAinc,1)	)
      if(size(vAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(vAinc,2)',size(vAinc,2)	)

      call die(myname_,'invalid argument(s)')
    endif

  endif
!________________________________________

	! Create a distributed AIGrid vector based on its specification
	! on the root PE.

  if(present(eagrid_rsrc)) then
    call rootedAIGrid_init(aig,im,jnp,KTpuv,root,comm,rsrc=eagrid_rsrc)
  else
    call rootedAIGrid_init(aig,im,jnp,KTpuv,root,comm,rsrc=psas_rsrc)
  endif

	ninc=distrSize(aig)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate(vinc)',ier)
		if(mall_ison()) call mall_mci(vinc,myname_)

  call AE_solve(ninc,vinc,		&
	ptr_lat(aig),ptr_lon(aig),ptr_lev(aig),ptr_kt(aig),	&
	nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,dels,		&
	comm,root,rsrc=psas_rsrc)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,aig)

		! Note that p-u-v ainc. is only defined on root PE.

  call rootedAIGrid_intp(aig,im,jnp,ktslp,    pAinc)
  call rootedAIGrid_intp(aig,im,jnp,ktUs,ktVs,uAinc,vAinc)

		if(mall_ison()) call mall_mco(vinc,myname_)
	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vinc)',ier)

  call clean(aig)

  if(myID==root.and.mpout_ison()) then
    call lvstat(mpout,im,jnp,uAinc,0.,'WIND','HGHT',1.e+15,'USL')
    call lvstat(mpout,im,jnp,vAinc,0.,'WIND','HGHT',1.e+15,'VSL')
    call lvstat(mpout,im,jnp,pAinc,0.,'PRES','HGHT',1.e+15,'PSL')
  endif

end subroutine solve_puv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_all_ - A multivariate analysis eqn. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_all_(im,jnp,mlev,plevs,			&
	pslAinc,uslAinc,vslAinc,zAinc,uAinc,vAinc,qAinc,	&
	nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,dels,		&
	comm,root,psas_rsrc,eagrid_rsrc)

      use m_rootedAIGrid,only : rootedAIGrid
      use m_rootedAIGrid,only : rootedAIGrid_init
      use m_rootedAIGrid,only : distrSize
      use m_rootedAIGrid,only : ptr_lat
      use m_rootedAIGrid,only : ptr_lon
      use m_rootedAIGrid,only : ptr_lev
      use m_rootedAIGrid,only : ptr_kt
      use m_rootedAIGrid,only : gatherv
      use m_rootedAIGrid,only : clean
      use m_rootedAIGrid,only : rootedAIGrid_intp

      use m_die,only : die,perr,MP_die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat
      use m_AE,    only : AE_solve

      implicit none

      integer,intent(in) :: im		! no. of longutide grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:),intent(out) :: pslAinc ! gridded p_sl analysis
      real,dimension(:,:),intent(out) :: uslAinc ! gridded u_sl analysis
      real,dimension(:,:),intent(out) :: vslAinc ! gridded v_sl analysis
      real,dimension(:,:,:),intent(out) :: zAinc ! gridded z analysis
      real,dimension(:,:,:),intent(out) :: uAinc ! gridded u analysis
      real,dimension(:,:,:),intent(out) :: vAinc ! gridded v analysis
      real,dimension(:,:,:),intent(out) :: qAinc ! gridded q analysis

      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(out):: xvec	! solution vector "x"
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU
      real   ,dimension(:),intent(in) :: dels	! an O-F vector

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_all_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedAIGrid) :: aig

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: inc
  integer :: ier

!________________________________________

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  call verifyObsVec_(nobs,xvec,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,dels,myname_)

  if(myID==root) then
    invalid=.false.
    if(.not.invalid) invalid = size(plevs)  < mlev

    if(.not.invalid) invalid = size(pslAinc,1) < im
    if(.not.invalid) invalid = size(pslAinc,2) < jnp

    if(.not.invalid) invalid = size(uslAinc,1) < im
    if(.not.invalid) invalid = size(uslAinc,2) < jnp

    if(.not.invalid) invalid = size(vslAinc,1) < im
    if(.not.invalid) invalid = size(vslAinc,2) < jnp

    if(.not.invalid) invalid = size(zAinc  ,1) < im
    if(.not.invalid) invalid = size(zAinc  ,2) < jnp
    if(.not.invalid) invalid = size(zAinc  ,3) < mlev

    if(.not.invalid) invalid = size(uAinc  ,1) < im
    if(.not.invalid) invalid = size(uAinc  ,2) < jnp
    if(.not.invalid) invalid = size(uAinc  ,3) < mlev

    if(.not.invalid) invalid = size(vAinc  ,1) < im
    if(.not.invalid) invalid = size(vAinc  ,2) < jnp
    if(.not.invalid) invalid = size(vAinc  ,3) < mlev

    if(.not.invalid) invalid = size(qAinc  ,1) < im
    if(.not.invalid) invalid = size(qAinc  ,2) < jnp
    if(.not.invalid) invalid = size(qAinc,  3) < mlev

    if(invalid) then
      if(size(plevs)  < mlev) call perr(myname_,	&
	'mlev',mlev,'size(plevs)',size(plevs)	)

      if(size(pslAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(pslAinc,1)',size(pslAinc,1)	)
      if(size(pslAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(pslAinc,2)',size(pslAinc,2)	)

      if(size(uslAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(uslAinc,1)',size(uslAinc,1)	)
      if(size(uslAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(uslAinc,2)',size(uslAinc,2)	)

      if(size(vslAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(vslAinc,1)',size(vslAinc,1)	)
      if(size(vslAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(vslAinc,2)',size(vslAinc,2)	)

      if(size(zAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(zAinc,1)',size(zAinc,1)	)
      if(size(zAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(zAinc,2)',size(zAinc,2)	)
      if(size(zAinc,3) < mlev) call perr(myname_,	&
	'mlev',mlev,'size(zAinc,3)',size(zAinc,3)	)

      if(size(uAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(uAinc,1)',size(uAinc,1)	)
      if(size(uAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(uAinc,2)',size(uAinc,2)	)
      if(size(uAinc,3) < mlev) call perr(myname_,	&
	'mlev',mlev,'size(uAinc,3)',size(uAinc,3)	)

      if(size(vAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(vAinc,1)',size(vAinc,1)	)
      if(size(vAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(vAinc,2)',size(vAinc,2)	)
      if(size(vAinc,3) < mlev) call perr(myname_,	&
	'mlev',mlev,'size(vAinc,3)',size(vAinc,3)	)

      if(size(qAinc,1) < im  ) call perr(myname_,	&
	'im'  ,im  ,'size(qAinc,1)',size(qAinc,1)	)
      if(size(qAinc,2) < jnp ) call perr(myname_,	&
	'jnp' ,jnp ,'size(qAinc,2)',size(qAinc,2)	)
      if(size(qAinc,3) < mlev) call perr(myname_,	&
	'mlev',mlev,'size(qAinc,3)',size(qAinc,3)	)

      call die(myname_,'invalid argument(s)')
    endif

  endif
!________________________________________

	! Create a distributed AIGrid vector based on its specification
	! on the root PE.

  if(present(eagrid_rsrc)) then
    call rootedAIGrid_init(aig,im,jnp,mlev,plevs,KTall,	&
	root,comm,rsrc=eagrid_rsrc)
  else
    call rootedAIGrid_init(aig,im,jnp,mlev,plevs,KTall,	&
	root,comm,rsrc=psas_rsrc)
  endif

	ninc=distrSize(aig)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate(vinc)',ier)
		if(mall_ison()) call mall_mci(vinc,myname_)

  call AE_solve(ninc,vinc,		&
	ptr_lat(aig),ptr_lon(aig),ptr_lev(aig),ptr_kt(aig),	&
	nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,dels,		&
	comm,root,rsrc=psas_rsrc)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,aig)

		! Note that the ainc. is only defined on root PE.

  call rootedAIGrid_intp(aig,im,jnp,ktUs,ktVs,uslAinc,vslAinc)
  call rootedAIGrid_intp(aig,im,jnp,ktslp,    pslAinc)

  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktUU,ktVV,uAinc,vAinc)
  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktHH,     zAinc)
  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktQQ,     qAinc)

		if(mall_ison()) call mall_mco(vinc,myname_)
	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vinc)',ier)

  call clean(aig)

  if(myID==root.and.mpout_ison()) then

    call lvstat(mpout,im,jnp,uslAinc,0.,'WIND','HGHT',1.e+15,'U_SL')
    call lvstat(mpout,im,jnp,vslAinc,0.,'WIND','HGHT',1.e+15,'V_SL')
    call lvstat(mpout,im,jnp,pslAinc,0.,'PRES','HGHT',1.e+15,'P_SL')

    inc=1
    if(mlev > 0) then
      if(plevs(1)>plevs(mlev)) inc=-1
    endif

    call gdstat(mpout,im,jnp,mlev,uAinc,plevs,'WIND','PRES',	&
	1.e+15,'Upperair u-wind analysis increments',inc	)
    call gdstat(mpout,im,jnp,mlev,vAinc,plevs,'WIND','PRES',	&
	1.e+15,'Upperair v-wind analysis increments',inc	)
    call gdstat(mpout,im,jnp,mlev,zAinc,plevs,'HGHT','PRES',	&
	1.e+15,'Upperair height analysis increments',inc	)

    call gdstat(mpout,im,jnp,mlev,qAinc,plevs,'MIXR','PRES',1.e+15,  &
	'Upperair water vapor mixing-ratio analysis increments',inc)
  endif

end subroutine solve_all_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_vec_ - A multivariate analysis eqn. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_vec_(ninc,vinc,vlat,vlon,vlev,vkts,	&
	nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,dels,		&
	comm,root,rsrc)

      use m_die,only : die,perr
      use m_AE, only : AE_solve
      implicit none

      integer,intent(in) :: ninc	! size of the increment vector
      real   ,dimension(:),intent(out):: vinc	! an increment vector
      real   ,dimension(:),intent(in) :: vlat	! latitudes
      real   ,dimension(:),intent(in) :: vlon	! longitudes
      real   ,dimension(:),intent(in) :: vlev	! levels
      integer,dimension(:),intent(in) :: vkts	! variable types

      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(out):: xvec	! solution vector "x"
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU
      real   ,dimension(:),intent(in) :: dels	! an O-F vector

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_vec_'

  logical :: invalid
  integer :: ier

  call verifyObsVec_(nobs,xvec,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,dels,myname_)
		!________________________________

  invalid=.false.
  if(.not.invalid) invalid = size(vinc) < ninc
  if(.not.invalid) invalid = size(vlat) < ninc
  if(.not.invalid) invalid = size(vlon) < ninc
  if(.not.invalid) invalid = size(vlev) < ninc
  if(.not.invalid) invalid = size(vkts) < ninc

  if(invalid) then
    if(size(vinc) < ninc) call perr(myname_,	&
	'ninc',ninc,'size(vinc)',size(vinc)	)
    if(size(vlat) < ninc) call perr(myname_,	&
	'ninc',ninc,'size(vlat)',size(vlat)	)
    if(size(vlon) < ninc) call perr(myname_,	&
	'ninc',ninc,'size(vlon)',size(vlon)	)
    if(size(vlev) < ninc) call perr(myname_,	&
	'ninc',ninc,'size(vlev)',size(vlev)	)
    if(size(vkts ) < ninc) call perr(myname_,	&
	'ninc',ninc,'size(vkts)',size(vkts)	)

    call die(myname_,'invalid argument(s)')
  endif

!________________________________________

  call AE_solve( ninc,vinc,vlat,vlon,vlev,vkts,			&
	nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,dels,		&
	comm,root,rsrc=rsrc)

end subroutine solve_vec_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_multi_ - A solver with multiiple variale type classes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_multi_(ninc,vinc,vlat,vlon,vlev,vkts,	&
	nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,dels,		&
	comm,root,ktClasses,rsrcFiles)

      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_die ,only : die,perr
      use m_AE  ,only : AE_solve
      implicit none

      integer,intent(in) :: ninc	! size of the increment vector
      real   ,dimension(:),intent(out):: vinc	! an increment vector
      real   ,dimension(:),intent(in) :: vlat	! latitudes
      real   ,dimension(:),intent(in) :: vlon	! longitudes
      real   ,dimension(:),intent(in) :: vlev	! levels
      integer,dimension(:),intent(in) :: vkts	! variable types

      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(out):: xvec	! solution vector "x"
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU
      real   ,dimension(:),intent(in) :: dels	! an O-F vector

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! Resource files for multiple types

      integer,         dimension(:),intent(in) :: ktClasses
      character(len=*),dimension(:),intent(in) :: rsrcFiles

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_multi_'
!________________________________________
!
! This interface is similar to solve_vec_(), except that the optional
! resource filename argumentin solve_vec_() becomes a pair of required
! arguments, to specify the resource filenames for each user specified
! analysis variable classes.  For example,
!
!	integer,dimension(4) :: ktClasses	=	&
!	  (/ 0,	& ! the base atmospheric variables class [1-7], and
!	     8,	& ! other univariate variable classes [7,:].
!	    12,	&
!	    16	/)
!
!	character(len=80),dimension(4) :: rsrcFiles	=	&
!	  (/"psas.rc",		&
!	    "psas-008.rc",	&
!	    "psas-012.rc",	&
!	    "psas-016.rc"	/)
!
!	call multiSolve(..., ktClasses=ktClasses,rsrcFiles=rsrcFiles)
!________________________________________

  integer,allocatable,dimension(:) :: kinc	! selected vinc(:)
  integer,allocatable,dimension(:) :: kobs	! selected xvec(:)
  real   ,allocatable,dimension(:) :: vtmp	! selected vinc
  real   ,allocatable,dimension(:) :: xtmp	! selected xvec
  integer :: minc,mobs
  logical :: invalid
  integer :: ikt,k
  integer :: ier
		!________________________________

  call verifyObsVec_(nobs,xvec,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,dels,myname_)
		!________________________________

  invalid=.false.
  if(.not.invalid) invalid = size(vinc) < ninc
  if(.not.invalid) invalid = size(vlat) < ninc
  if(.not.invalid) invalid = size(vlon) < ninc
  if(.not.invalid) invalid = size(vlev) < ninc
  if(.not.invalid) invalid = size(vkts) < ninc

  if(invalid) then
    if(size(vinc) < ninc) call perr(myname_,	&
	'ninc',ninc,'size(vinc)',size(vinc)	)
    if(size(vlat) < ninc) call perr(myname_,	&
	'ninc',ninc,'size(vlat)',size(vlat)	)
    if(size(vlon) < ninc) call perr(myname_,	&
	'ninc',ninc,'size(vlon)',size(vlon)	)
    if(size(vlev) < ninc) call perr(myname_,	&
	'ninc',ninc,'size(vlev)',size(vlev)	)
    if(size(vkts ) < ninc) call perr(myname_,	&
	'ninc',ninc,'size(vkts)',size(vkts)	)

    call die(myname_,'invalid argument(s)')
  endif
		!________________________________

  do ikt=1,size(ktClasses)
    select case(ktClasses(ikt))
    case(1,2,3,4,5,6,7)

		! Any memeber in Class 0 should not be specified
		! seperately.

      call die(myname_,'invalid kt class ID',ktClasses(ikt))
    end select
  end do
!________________________________________

	allocate(kinc(ninc),kobs(nobs),stat=ier)
		if(ier/=0) call die(myname_,'allocate(kinc,kobs)',ier)
		if(mall_ison()) then
		  call mall_mci(kinc,myname_)
		  call mall_mci(kobs,myname_)
		endif

  do ikt=1,size(ktClasses)
		! Mark kt values in a specified class

    call sel_types_(ktClasses(ikt), vkts(1:ninc), minc, kinc)
    call sel_types_(ktClasses(ikt), kts (1:nobs), mobs, kobs)

	allocate(vtmp(minc),xtmp(mobs),stat=ier)
		if(ier/=0) call die(myname_,'allocate(vtmp,xtmp)',ier)
		if(mall_ison()) then
		  call mall_mci(vtmp,myname_)
		  call mall_mci(xtmp,myname_)
		endif

		! Perform analysis with marked kt values only

    select case(ktClasses(ikt))
    case(0)
		! This class is for kt = [1,7], as the standard
		! configuration for the atmosphere.

      call AE_solve( minc,vtmp(:),		&
	vlat(kinc(1:minc)),vlon(kinc(1:minc)),	&
	vlev(kinc(1:minc)),vkts(kinc(1:minc)),	&
		     mobs,xtmp(:),		&
	rlat(kobs(1:mobs)),rlon(kobs(1:mobs)),	&
	rlev(kobs(1:mobs)),			&
	kxs (kobs(1:mobs)),kss (kobs(1:mobs)),	&
	kts (kobs(1:mobs)),			&
	xmUs(kobs(1:mobs)),dels(kobs(1:mobs)),	&
	comm,root,rsrc=rsrcFiles(ikt))

    case(1,2,3,4,5,6,7)

		! Any memeber in Class 0 should not be specified
		! seperately.

      call die(myname_,'invalid kt class ID',ktClasses(ikt))

    case default

		! Any other member will be handled as a univariate case
		! with its resource file specified in a way similar to
		! the mixing ratio variable (kt=7).

      call AE_solve( minc,vtmp(:),			&
	vlat(kinc(1:minc) ),vlon(kinc(1:minc) ),	&
	vlev(kinc(1:minc) ),   (/(7,k=1,minc)/),	&
		     mobs,xtmp(:),			&
	rlat(kobs(1:mobs) ),rlon(kobs(1:mobs) ),	&
	rlev(kobs(1:mobs) ),				&
	kxs (kobs(1:mobs) ),kss (kobs(1:mobs) ),	&
	   (/(7,k=1,mobs)/),				&
	xmUs(kobs(1:mobs) ),dels(kobs(1:mobs)),		&
	comm,root,rsrc=rsrcFiles(ikt))

    end select

		! Store the result for marked kt values only

    vinc(kinc(1:minc)) = vtmp(1:minc)
    xvec(kinc(1:mobs)) = xtmp(1:mobs)

		if(mall_ison()) then
		  call mall_mco(vtmp,myname_)
		  call mall_mco(xtmp,myname_)
		endif
	deallocate(vtmp,xtmp,stat=ier)
		if(ier/=0) call die(myname_,'deallocate(vtmp,xtmp)',ier)

  end do

		if(mall_ison()) then
		  call mall_mco(kinc,myname_)
		  call mall_mco(kobs,myname_)
		endif
	deallocate(kinc,kobs,stat=ier)
		if(ier/=0) call die(myname_,'deallocate(kinc,kobs)',ier)

end subroutine solve_multi_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sel_types_ - select kts in a specified class
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine sel_types_(ktClass, vkts, minc, kinc)
      use m_die,only : die
      implicit none
      integer,intent(in) :: ktClass
      integer,dimension(:),intent(in)  :: vkts	! all kt values
      integer,intent(out) :: minc		! size of selected
      integer,dimension(:),intent(out) :: kinc	! selected elements

! !REVISION HISTORY:
! 	14Dec01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sel_types_'
  integer :: nout
  integer :: i

  nout=size(kinc)
  minc=0

  do i=1,size(vkts)
    select case(ktClass)
    case(0)
      if( vkts(i)>=1 .and. vkts(i)<=7 ) then
	minc=minc+1
	if(minc<=nout) kinc(minc)=i
      endif
    case(1,2,3,4,5,6,7)
    case default
      if(vkts(i)==ktClass) then
	minc=minc+1
	if(minc<=nout) kinc(minc)=i
      endif
    end select
  end do

  if(minc>nout) call die(myname_,'size(kinc)',nout,'minc',minc)

end subroutine sel_types_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: verifyObsVec_ - verify obs. vec. dimensions
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine verifyObsVec_(nobs,xvec,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,dels, where)

      use m_die,only : die,perr
      implicit none

      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(in) :: xvec	! solution vector "x"
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU
      real   ,dimension(:),intent(in) :: dels	! an O-F vector

      character(len=*),intent(in) :: where

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::verifyObsVec_'

  logical :: invalid
  integer :: ier

  invalid=.false.

  if(.not.invalid) invalid = size(xvec) < nobs
  if(.not.invalid) invalid = size(rlat) < nobs
  if(.not.invalid) invalid = size(rlon) < nobs
  if(.not.invalid) invalid = size(rlev) < nobs
  if(.not.invalid) invalid = size(kxs ) < nobs
  if(.not.invalid) invalid = size(kss ) < nobs
  if(.not.invalid) invalid = size(kts ) < nobs
  if(.not.invalid) invalid = size(xmUs) < nobs
  if(.not.invalid) invalid = size(dels) < nobs

  if(invalid) then
    if(size(xvec) < nobs) call perr(where,	&
	'nobs',nobs,'size(xvec)',size(xvec)	)
    if(size(rlat) < nobs) call perr(where,	&
	'nobs',nobs,'size(rlat)',size(rlat)	)
    if(size(rlon) < nobs) call perr(where,	&
	'nobs',nobs,'size(rlon)',size(rlon)	)
    if(size(rlev) < nobs) call perr(where,	&
	'nobs',nobs,'size(rlev)',size(rlev)	)
    if(size(kxs ) < nobs) call perr(where,	&
	'nobs',nobs,'size(kxs)',size(kxs)	)
    if(size(kss ) < nobs) call perr(where,	&
	'nobs',nobs,'size(kss)',size(kss)	)
    if(size(kts ) < nobs) call perr(where,	&
	'nobs',nobs,'size(kts)',size(kts)	)
    if(size(xmUs) < nobs) call perr(where,	&
	'nobs',nobs,'size(xmUs)',size(xmUs)	)
    if(size(dels) < nobs) call perr(where,	&
	'nobs',nobs,'size(dels)',size(dels)	)

    call die(where,'invalid argument(s)')
  endif
end subroutine verifyObsVec_

end module m_AISolver
