!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_PHxSolver - PSAS solver interfaces
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_PHxSolver
      use m_ktList,only : ktslp,ktus,ktvs,ktHH,ktUU,ktVV,ktQQ
      implicit none
      private	! except

      public :: PHx

      interface PHx; module procedure	&
	PHx_mix_,	&
	PHx_puv_,	&
	PHx_zuv_,	&
	PHx_all_
      end interface

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_PHxSolver'

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
! !IROUTINE: PHx_mix_ - A univariate upper-air FcstErrCov eq. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine PHx_mix_(nobs,rlat,rlon,rlev,kts,xvec,	&
	im,jnp,mlev,plevs,qAinc,comm,root,rsrc)

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
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
      use m_FcstErrCov,only : FcstErrCov_Cx
      use m_mpout, only : mpout,mpout_ison
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat

      implicit none
      integer,intent(in) :: nobs		! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels (pressure)
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xvec	! an O-F vector

		! These output avariables are only significant on the
		! root PE.

      integer,intent(in) :: im		! no. of longitude grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:,:),intent(out) :: qAinc ! gridded analysis

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! who the root PE is.
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::PHx_mix_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedAIGrid) :: aig

  integer :: ninc
  real,allocatable,dimension(:) :: vinc
  real   ,pointer ,dimension(:) :: vlat
  real   ,pointer ,dimension(:) :: vlon
  real   ,pointer ,dimension(:) :: vlev
  integer,pointer ,dimension(:) :: vkts

  integer :: myID
  logical :: invalid
  integer :: inc
  integer :: ier

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  invalid=.false.
  if(.not.invalid) invalid = size(rlat) < nobs
  if(.not.invalid) invalid = size(rlon) < nobs
  if(.not.invalid) invalid = size(rlev) < nobs
  if(.not.invalid) invalid = size(xvec) < nobs
  if(.not.invalid) invalid = size(kts ) < nobs

  if(myID==root) then
    if(.not.invalid) invalid = size(plevs)  < mlev

    if(.not.invalid) invalid = size(qAinc,1) < im
    if(.not.invalid) invalid = size(qAinc,2) < jnp
    if(.not.invalid) invalid = size(qAinc,3) < mlev
  endif

  if(invalid) then
    if(size(rlat) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlat)',size(rlat)	)
    if(size(rlon) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlon)',size(rlon)	)
    if(size(rlev) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlev)',size(rlev)	)
    if(size(xvec) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(xvec)',size(xvec)	)
    if(size(kts ) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(kts)',size(kts)	)

    if(myID==root) then
      if(size(plevs)  < mlev) call perr(myname_,	&
	'mlev',mlev,'size(plevs)',size(plevs)	)

      if(size(qAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(qAinc,1)',size(qAinc,1)	)
      if(size(qAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(qAinc,2)',size(qAinc,2)	)
      if(size(qAinc,3) < mlev) call perr(myname_,		&
	'mlev',mlev,'size(qAinc,3)',size(qAinc,3)	)
    endif

    call die(myname_,'invalid argument(s)')
  endif

!________________________________________

	! Created a distributed AIGrid vector based on its specification
	! on the root PE.

  call rootedAIGrid_init(aig,im,jnp,mlev,plevs,KTmix,root,comm)

		! Refer to local attribute data

	vlat => ptr_lat(aig)
	vlon => ptr_lon(aig)
	vlev => ptr_lev(aig)
	vkts => ptr_kt(aig)

	ninc=distrSize(aig)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call FcstErrCov_Cx(	&
	ninc, vinc, vlat, vlon,	vlev, vkts,	&
	nobs, xvec, rlat, rlon, rlev,  kts, comm,root,rsrc=rsrc)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,aig)

		! Note that qAinc is only defined on root PE.

  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktQQ,qAinc)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

	nullify(vlat)
	nullify(vlon)
	nullify(vlev)
	nullify(vkts)

  call clean(aig)

  if(myID==root.and.mpout_ison()) then
    inc=1
    if(mlev > 0) then
      if(plevs(1)>plevs(mlev)) inc=-1
    endif

    call gdstat(mpout,im,jnp,mlev,qAinc,plevs,'MIXR','PRES',1.e+15,  &
	'Upperair water vapor mixing-ratio analysis increments',inc)
  endif

end subroutine PHx_mix_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PHx_zuv_ - A multivariate upper-air FcstErrCov eq. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine PHx_zuv_(nobs,rlat,rlon,rlev,kts,xvec,	&
	im,jnp,mlev,plevs,zAinc,uAinc,vAinc,comm,root,rsrc)

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
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
      use m_FcstErrCov,only : FcstErrCov_Cx
      use m_mpout, only : mpout,mpout_ison
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat

      implicit none
      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels (pressure)
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xvec	! an O-F vector

		! These output avariables are only significant on the
		! root PE.

      integer,intent(in) :: im		! no. of longitude grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:,:),intent(out) :: zAinc ! Gridded h anaysis
      real,dimension(:,:,:),intent(out) :: uAinc ! gridded u analysis
      real,dimension(:,:,:),intent(out) :: vAinc ! gridded v analysis

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! who the root PE is.
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::PHx_zuv_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedAIGrid) :: aig

  integer :: ninc
  real,allocatable,dimension(:) :: vinc
  real   ,pointer ,dimension(:) :: vlat
  real   ,pointer ,dimension(:) :: vlon
  real   ,pointer ,dimension(:) :: vlev
  integer,pointer ,dimension(:) :: vkts

  integer :: myID
  logical :: invalid
  integer :: inc
  integer :: ier

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  invalid=.false.
  if(.not.invalid) invalid = size(rlat) < nobs
  if(.not.invalid) invalid = size(rlon) < nobs
  if(.not.invalid) invalid = size(rlev) < nobs
  if(.not.invalid) invalid = size(xvec) < nobs
  if(.not.invalid) invalid = size(kts ) < nobs

  if(myID==root) then
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
  endif

  if(invalid) then
    if(size(rlat) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(rlat)',size(rlat)	)
    if(size(rlon) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(rlon)',size(rlon)	)
    if(size(rlev) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(rlev)',size(rlev)	)
    if(size(xvec) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(xvec)',size(xvec)	)
    if(size(kts ) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(kts)',size(kts)	)

    if(myID==root) then
      if(size(plevs)  < mlev) call perr(myname_,	&
		'mlev',mlev,'size(plevs)',size(plevs)	)

      if(size(zAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(zAinc,1)',size(zAinc,1)	)
      if(size(zAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(zAinc,2)',size(zAinc,2)	)
      if(size(zAinc,3) < mlev) call perr(myname_,		&
	'mlev',mlev,'size(zAinc,3)',size(zAinc,3)	)

      if(size(uAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(uAinc,1)',size(uAinc,1)	)
      if(size(uAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(uAinc,2)',size(uAinc,2)	)
      if(size(uAinc,3) < mlev) call perr(myname_,		&
	'mlev',mlev,'size(uAinc,3)',size(uAinc,3)	)

      if(size(vAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(vAinc,1)',size(vAinc,1)	)
      if(size(vAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(vAinc,2)',size(vAinc,2)	)
      if(size(vAinc,3) < mlev) call perr(myname_,		&
	'mlev',mlev,'size(vAinc,3)',size(vAinc,3)	)
    endif

    call die(myname_,'invalid argument(s)')
  endif

!________________________________________

	! Created a distributed AIGrid vector based on its specification
	! on the root PE.

  call rootedAIGrid_init(aig,im,jnp,mlev,plevs,KTzuv,root,comm)

		! Refer to local attribute data

	vlat => ptr_lat(aig)
	vlon => ptr_lon(aig)
	vlev => ptr_lev(aig)
	vkts => ptr_kt(aig)

	ninc=distrSize(aig)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call FcstErrCov_Cx(	&
	ninc, vinc, vlat, vlon,	vlev, vkts,	&
	nobs, xvec, rlat, rlon, rlev,  kts, comm,root,rsrc=rsrc)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,aig)

		! Note that h-u-v ainc. is only defined on root PE.

  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktHH,     zAinc)
  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktUU,ktVV,uAinc,vAinc)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

	nullify(vlat)
	nullify(vlon)
	nullify(vlev)
	nullify(vkts)

  call clean(aig)

  if(myID==root.and.mpout_ison()) then
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
  endif

end subroutine PHx_zuv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PHx_puv_ - A multivariate sea-level FcstErrCov eq. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine PHx_puv_(nobs,rlat,rlon,rlev,kts,xvec,	&
	im,jnp,pAinc,uAinc,vAinc,comm,root,rsrc)

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
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
      use m_FcstErrCov,only : FcstErrCov_Cx
      use m_mpout, only : mpout,mpout_ison
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat

      implicit none
      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels (pressures)
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xvec	! an O-F vector

		! These output avariables are only significant on the
		! root PE.

      integer,intent(in) :: im		! no. of longitude grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      real,dimension(:,:),intent(out) :: pAinc	! gridded p_sl analysis
      real,dimension(:,:),intent(out) :: uAinc	! gridded u_sl analysis
      real,dimension(:,:),intent(out) :: vAinc	! gridded v_sl analysis

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! who the root PE is.
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::PHx_puv_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedAIGrid) :: aig

  integer :: ninc
  real,allocatable,dimension(:) :: vinc
  real   ,pointer ,dimension(:) :: vlat
  real   ,pointer ,dimension(:) :: vlon
  real   ,pointer ,dimension(:) :: vlev
  integer,pointer ,dimension(:) :: vkts

  integer :: myID
  logical :: invalid
  integer :: ier

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  invalid=.false.
  if(.not.invalid) invalid = size(rlat) < nobs
  if(.not.invalid) invalid = size(rlon) < nobs
  if(.not.invalid) invalid = size(rlev) < nobs
  if(.not.invalid) invalid = size(xvec) < nobs
  if(.not.invalid) invalid = size(kts ) < nobs

  if(myID==root) then
    if(.not.invalid) invalid = size(pAinc,1) < im
    if(.not.invalid) invalid = size(pAinc,2) < jnp
    if(.not.invalid) invalid = size(uAinc,1) < im
    if(.not.invalid) invalid = size(uAinc,2) < jnp
    if(.not.invalid) invalid = size(vAinc,1) < im
    if(.not.invalid) invalid = size(vAinc,2) < jnp
  endif

  if(invalid) then
    if(size(rlat) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlat)',size(rlat)	)
    if(size(rlon) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlon)',size(rlon)	)
    if(size(rlev) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlev)',size(rlev)	)
    if(size(xvec) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(xvec)',size(xvec)	)
    if(size(kts ) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(kts)',size(kts)	)

    if(myID==root) then
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
    endif

    call die(myname_,'invalid argument(s)')
  endif

!________________________________________

	! Created a distributed AIGrid vector based on its specification
	! on the root PE.

  call rootedAIGrid_init(aig,im,jnp,KTpuv,root,comm)

		! Refer to local attribute data

	vlat => ptr_lat(aig)
	vlon => ptr_lon(aig)
	vlev => ptr_lev(aig)
	vkts => ptr_kt(aig)

	ninc=distrSize(aig)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call FcstErrCov_Cx(	&
	ninc, vinc, vlat, vlon,	vlev, vkts,	&
	nobs, xvec, rlat, rlon, rlev,  kts, comm,root,rsrc=rsrc)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,aig)

		! Note that p-u-v ainc. is only defined on root PE.

  call rootedAIGrid_intp(aig,im,jnp,ktslp,    pAinc)
  call rootedAIGrid_intp(aig,im,jnp,ktUs,ktVs,uAinc,vAinc)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

	nullify(vlat)
	nullify(vlon)
	nullify(vlev)
	nullify(vkts)

  call clean(aig)

  if(myID==root.and.mpout_ison()) then
    call lvstat(mpout,im,jnp,uAinc,0.,'WIND','HGHT',1.e+15,'USL')
    call lvstat(mpout,im,jnp,vAinc,0.,'WIND','HGHT',1.e+15,'VSL')
    call lvstat(mpout,im,jnp,pAinc,0.,'PRES','HGHT',1.e+15,'PSL')
  endif

end subroutine PHx_puv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PHx_all_ - A multivariate FcstErrCov eq. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine PHx_all_(nobs,rlat,rlon,rlev,kts,xvec,	&
	im,jnp,mlev,plevs,pslAinc,uslAinc,vslAinc,	&
			  zAinc,uAinc,vAinc,qAinc,comm,root,rsrc)

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
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
      use m_FcstErrCov,only : FcstErrCov_Cx
      use m_mpout, only : mpout,mpout_ison
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat

      implicit none
      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels (pressures)
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xvec	! an O-F vector

		! These output avariables are only significant on the
		! root PE.

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

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! who the root PE is.
      ! character(len=*),optional,intent(in) :: rsrc
      character(len=*),intent(in) :: rsrc

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::PHx_all_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedAIGrid) :: aig

  integer :: ninc
  real,allocatable,dimension(:) :: vinc
  real   ,pointer ,dimension(:) :: vlat
  real   ,pointer ,dimension(:) :: vlon
  real   ,pointer ,dimension(:) :: vlev
  integer,pointer ,dimension(:) :: vkts

  integer :: myID
  logical :: invalid
  integer :: inc
  integer :: ier

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  invalid=.false.
  if(.not.invalid) invalid = size(rlat) < nobs
  if(.not.invalid) invalid = size(rlon) < nobs
  if(.not.invalid) invalid = size(rlev) < nobs
  if(.not.invalid) invalid = size(xvec) < nobs
  if(.not.invalid) invalid = size(kts ) < nobs

  if(myID==root) then
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
  endif

  if(invalid) then
    if(size(rlat) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlat)',size(rlat)	)
    if(size(rlon) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlon)',size(rlon)	)
    if(size(rlev) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlev)',size(rlev)	)
    if(size(xvec) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(xvec)',size(xvec)	)
    if(size(kts ) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(kts)',size(kts)	)

    if(myID==root) then
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

      if(size(zAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(zAinc,1)',size(zAinc,1)	)
      if(size(zAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(zAinc,2)',size(zAinc,2)	)
      if(size(zAinc,3) < mlev) call perr(myname_,		&
	'mlev',mlev,'size(zAinc,3)',size(zAinc,3)	)

      if(size(uAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(uAinc,1)',size(uAinc,1)	)
      if(size(uAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(uAinc,2)',size(uAinc,2)	)
      if(size(uAinc,3) < mlev) call perr(myname_,		&
	'mlev',mlev,'size(uAinc,3)',size(uAinc,3)	)

      if(size(vAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(vAinc,1)',size(vAinc,1)	)
      if(size(vAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(vAinc,2)',size(vAinc,2)	)
      if(size(vAinc,3) < mlev) call perr(myname_,		&
	'mlev',mlev,'size(vAinc,3)',size(vAinc,3)	)

      if(size(qAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(qAinc,1)',size(qAinc,1)	)
      if(size(qAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(qAinc,2)',size(qAinc,2)	)
      if(size(qAinc,3) < mlev) call perr(myname_,		&
	'mlev',mlev,'size(qAinc,3)',size(qAinc,3)	)
    endif

    call die(myname_,'invalid argument(s)')
  endif

!________________________________________

	! Created a distributed AIGrid vector based on its specification
	! on the root PE.

  call rootedAIGrid_init(aig,im,jnp,mlev,plevs,KTall,root,comm)

		! Refer to local attribute data

	vlat => ptr_lat(aig)
	vlon => ptr_lon(aig)
	vlev => ptr_lev(aig)
	vkts => ptr_kt(aig)

	ninc=distrSize(aig)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call FcstErrCov_Cx(	&
	ninc, vinc, vlat, vlon,	vlev, vkts,	&
	nobs, xvec, rlat, rlon, rlev,  kts, comm,root,rsrc=rsrc)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,aig)

		! Note that the ainc. is only defined on root PE.

  call rootedAIGrid_intp(aig,im,jnp,ktUs,ktVs,uslAinc,vslAinc)
  call rootedAIGrid_intp(aig,im,jnp,ktslp,    pslAinc)

  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktUU,ktVV,uAinc,vAinc)
  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktHH,     zAinc)
  call rootedAIGrid_intp(aig,im,jnp,mlev,plevs,ktQQ,     qAinc)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

	nullify(vlat)
	nullify(vlon)
	nullify(vlev)
	nullify(vkts)

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

end subroutine PHx_all_

end module m_PHxSolver
