!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AISolverNoSuperOb - PSAS solver interfaces (w/o superob)
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_AISolverNoSuperOb
      use m_ktList,only : ktslp,ktus,ktvs,ktHH,ktUU,ktVV,ktQQ
      implicit none
      private	! except

      public :: getAInoSuperOb

      interface getAInoSuperOb; module procedure	&
	getAImix_,	&
	getAIpuv_,	&
	getAIzuv_,	&
	getAIall_
      end interface

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	01Jul00 - Guo/RT, bug fix, ks now made a dummy argument
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AISolverNoSuperOb'

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
! !IROUTINE: getAImix_ - A univariate upper-air analysis equation solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getAImix_(nobs,rlat,rlon,rlev,tims,kxs,kts,kss,dels,	&
	im,jnp,mlev,plevs,qAinc,xvec				)

      use m_die,only : die,perr
      use m_AIGrid,only : AIGrid_init
      use m_AIGrid,only : AIGrid_clean
      use m_AIGrid,only : AIGrid_bind
      use m_AIGrid,only : AIGrid_intp
      use m_AE,    only : AE_solve
      use m_ObsSmry,only : ObsSmry
      use m_mpout, only : mpout,mpout_ison

      implicit none
      integer,intent(in) :: nobs		! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels (pressure)
      real   ,dimension(:),intent(in) :: tims	! time offsets
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      integer,dimension(:),intent(in) :: kss    ! sounding index
      real   ,dimension(:),intent(in) :: dels	! an O-F vector

      integer,intent(in) :: im		! no. of longitude grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:,:),intent(out) :: qAinc ! gridded analysis
      real ,dimension(:),intent(out) :: xvec	 ! solution of conjgr

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	09May00 - Todling - Removed superobbing
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getAImix_'

	! The workspace for a local copy of the O-F vector and its
	! artributes.

  real,allocatable,dimension(:) :: xmUs_

	! The workspace for a local copy of A-F vector and its
	! artributes.

  integer :: ninc
  real,allocatable,dimension(:) :: vinc
  real,allocatable,dimension(:) :: vlat
  real,allocatable,dimension(:) :: vlon
  real,allocatable,dimension(:) :: vlev
  integer,allocatable,dimension(:) :: vkts

  logical :: invalid
  integer :: inc
  integer :: ier
  integer :: i

  invalid=.false.
  if(.not.invalid) invalid = size(rlat) < nobs
  if(.not.invalid) invalid = size(rlon) < nobs
  if(.not.invalid) invalid = size(rlev) < nobs
  if(.not.invalid) invalid = size(tims) < nobs
  if(.not.invalid) invalid = size(dels) < nobs
  if(.not.invalid) invalid = size(xvec) < nobs
  if(.not.invalid) invalid = size(kxs ) < nobs
  if(.not.invalid) invalid = size(kts ) < nobs
  if(.not.invalid) invalid = size(kss ) < nobs

  if(.not.invalid) invalid = size(plevs)  < mlev

  if(.not.invalid) invalid = size(qAinc,1) < im
  if(.not.invalid) invalid = size(qAinc,2) < jnp
  if(.not.invalid) invalid = size(qAinc,3) < mlev

  if(invalid) then
    if(size(rlat) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlat)',size(rlat)	)
    if(size(rlon) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlon)',size(rlon)	)
    if(size(rlev) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlev)',size(rlev)	)
    if(size(tims) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(tims)',size(tims)	)
    if(size(dels) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(dels)',size(dels)	)
    if(size(kxs ) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(kxs)',size(kxs)	)
    if(size(kts ) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(kts)',size(kts)	)
    if(size(xvec) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(xvec)',size(xvec)	)	

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

!________________________________________

	allocate( xmUs_(nobs), stat=ier	)
		if(ier/=0) call die(myname_,'allocate(obs)',ier)

	if(mpout_ison()) then
	  write(mpout,'(2a,i2)') myname_,': on entry()'
	  call ObsSmry(mpout,nobs,kxs,kts)
	endif

!  Fill up arrays that normally would come from superObs
!  -----------------------------------------------------
  xmUs_(1:nobs) = 1.0  			! Scale sigO factor 

  call AIGrid_init(ninc,im,jnp,mlev,plevs,size(KTmix),KTmix )

	allocate( vinc(ninc),vlat(ninc),vlon(ninc),	&
		  vlev(ninc),vkts(ninc), stat=ier	)
		if(ier/=0) call die(myname_,'allocate(vec)',ier)

  call AIGrid_bind(ninc, vlat,vlon,vlev,vkts )

  call AE_solve(	&
	ninc, vinc, vlat, vlon,	vlev, vkts,	&
	nobs, xvec, rlat, rlon, rlev, tims,	&
	       kxs,  kss,  kts,xmUs_,dels	)

  call AIGrid_intp(ktQQ,im,jnp,mlev,plevs,qAinc, ninc,vinc )

	deallocate(vinc,vlat,vlon,vlev,vkts, stat=ier )

		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

  call AIGrid_clean()

	deallocate(xmUs_,	stat=ier	)
		if(ier/=0) call die(myname_,'deallocate(obs)',ier)

  inc=1
  if(mlev > 0) then
    if(plevs(1)>plevs(mlev)) inc=-1
  endif

  call gdstat(mpout,im,jnp,mlev,qAinc,plevs,'MIXR','PRES',1.e+15,  &
	'Upperair water vapor mixing-ratio analysis increments',inc)

end subroutine getAImix_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getAIzuv_ - A multivariate upper-air analysis eqn. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getAIzuv_(nobs,rlat,rlon,rlev,tims,kxs,kts,kss,dels,	&
	im,jnp,mlev,plevs,zAinc,uAinc,vAinc,xvec		)

      use m_die,only : die,perr
      use m_AIGrid,only : AIGrid_init
      use m_AIGrid,only : AIGrid_clean
      use m_AIGrid,only : AIGrid_bind
      use m_AIGrid,only : AIGrid_intp
      use m_AE,    only : AE_solve
      use m_ObsSmry,only : ObsSmry
      use m_mpout, only : mpout,mpout_ison

      implicit none
      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels (pressure)
      real   ,dimension(:),intent(in) :: tims	! time off-sets
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      integer,dimension(:),intent(in) :: kss    ! sounding index
      real   ,dimension(:),intent(in) :: dels	! an O-F vector

      integer,intent(in) :: im		! no. of longitude grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:,:),intent(out) :: zAinc ! Gridded h anaysis
      real,dimension(:,:,:),intent(out) :: uAinc ! gridded u analysis
      real,dimension(:,:,:),intent(out) :: vAinc ! gridded v analysis
      real,dimension(:),intent(out) :: xvec	 ! solution of conjgr

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
! 	09May00 - Todling - Removed superobbing
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getAIzuv_'

	! The workspace for a local copy of the O-F vector and its
	! artributes.

  real,allocatable,dimension(:) :: xmUs_

	! The workspace for a local copy of A-F vector and its
	! artributes.

  integer :: ninc
  real,allocatable,dimension(:) :: vinc
  real,allocatable,dimension(:) :: vlat
  real,allocatable,dimension(:) :: vlon
  real,allocatable,dimension(:) :: vlev
  integer,allocatable,dimension(:) :: vkts

  logical :: invalid
  integer :: inc
  integer :: ier
  integer :: i

  invalid=.false.
  if(.not.invalid) invalid = size(rlat) < nobs
  if(.not.invalid) invalid = size(rlon) < nobs
  if(.not.invalid) invalid = size(rlev) < nobs
  if(.not.invalid) invalid = size(tims) < nobs
  if(.not.invalid) invalid = size(dels) < nobs
  if(.not.invalid) invalid = size(xvec) < nobs
  if(.not.invalid) invalid = size(kxs ) < nobs
  if(.not.invalid) invalid = size(kts ) < nobs
  if(.not.invalid) invalid = size(kss ) < nobs

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
    if(size(rlat) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(rlat)',size(rlat)	)
    if(size(rlon) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(rlon)',size(rlon)	)
    if(size(rlev) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(rlev)',size(rlev)	)
    if(size(tims) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(tims)',size(tims)	)
    if(size(dels) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(dels)',size(dels)	)
    if(size(kxs ) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(kxs)',size(kxs)	)
    if(size(kts ) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(kts)',size(kts)	)
    if(size(xvec) < nobs) call perr(myname_,	&
		'nobs',nobs,'size(xvec)',size(xvec)	)		

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

    call die(myname_,'invalid argument(s)')
  endif

!________________________________________

	allocate(  xmUs_(nobs), stat=ier	)
		if(ier/=0) call die(myname_,'allocate(obs)',ier)

	if(mpout_ison()) then
	  write(mpout,'(2a,i2)') myname_,': on entry'
	  call ObsSmry(mpout,nobs,kxs,kts)
	endif

!  Fill up arrays that normally would come from superObs
!  -----------------------------------------------------
  xmUs_(1:nobs) = 1.0  			! Scale sigO factor
	
  call AIGrid_init(ninc,im,jnp,mlev,plevs,size(KTzuv),KTzuv )

	allocate( vinc(ninc),vlat(ninc),vlon(ninc),	&
		  vlev(ninc),vkts(ninc), stat=ier	)
		if(ier/=0) call die(myname_,'allocate(vec)',ier)

  call AIGrid_bind(ninc, vlat,vlon,vlev,vkts )

  call AE_solve(	&
	ninc, vinc, vlat, vlon,	vlev, vkts,	&
	nobs, xvec, rlat, rlon, rlev, tims,	&
	       kxs,  kss,  kts,xmUs_,dels 	)

  call AIGrid_intp(ktHH,     im,jnp,mlev,plevs,zAinc,       ninc,vinc)
  call AIGrid_intp(ktUU,ktVV,im,jnp,mlev,plevs,uAinc,vAinc, ninc,vinc)

	deallocate(vinc,vlat,vlon,vlev,vkts, stat=ier )

		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

  call AIGrid_clean()

	deallocate(xmUs_,	stat=ier	)

		if(ier/=0) call die(myname_,'deallocate(obs)',ier)

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

end subroutine getAIzuv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getAIpuv_ - A multivariate sea-level analysis eqn. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getAIpuv_(nobs,rlat,rlon,rlev,tims,kxs,kts,kss,dels,	&
	im,jnp,pAinc,uAinc,vAinc,xvec				)

      use m_die,only : die,perr
      use m_AIGrid,only : AIGrid_init
      use m_AIGrid,only : AIGrid_clean
      use m_AIGrid,only : AIGrid_bind
      use m_AIGrid,only : AIGrid_intp
      use m_AE,    only : AE_solve
      use m_ObsSmry,only : ObsSmry
      use m_mpout, only : mpout,mpout_ison

      implicit none
      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels (pressures)
      real   ,dimension(:),intent(in) :: tims	! time off-sets
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      integer,dimension(:),intent(in) :: kss    ! sounding index
      real   ,dimension(:),intent(in) :: dels	! an O-F vector

      integer,intent(in) :: im		! no. of longitude grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      real,dimension(:,:),intent(out) :: pAinc	! gridded p_sl analysis
      real,dimension(:,:),intent(out) :: uAinc	! gridded u_sl analysis
      real,dimension(:,:),intent(out) :: vAinc	! gridded v_sl analysis
      real,dimension(:),intent(out) :: xvec	! solution of conjgr

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
! 	09May00 - Todling - Removed superobbing
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getAIpuv_'

	! The workspace for a local copy of the O-F vector and its
	! artributes.

  real,allocatable,dimension(:) :: xmUs_

	! The workspace for a local copy of A-F vector and its
	! artributes.

  integer :: ninc
  real,allocatable,dimension(:) :: vinc
  real,allocatable,dimension(:) :: vlat
  real,allocatable,dimension(:) :: vlon
  real,allocatable,dimension(:) :: vlev
  integer,allocatable,dimension(:) :: vkts

  logical :: invalid
  integer :: ier
  integer :: i

  invalid=.false.
  if(.not.invalid) invalid = size(rlat) < nobs
  if(.not.invalid) invalid = size(rlon) < nobs
  if(.not.invalid) invalid = size(rlev) < nobs
  if(.not.invalid) invalid = size(tims) < nobs
  if(.not.invalid) invalid = size(dels) < nobs
  if(.not.invalid) invalid = size(xvec) < nobs
  if(.not.invalid) invalid = size(kxs ) < nobs
  if(.not.invalid) invalid = size(kts ) < nobs
  if(.not.invalid) invalid = size(kss ) < nobs
  if(.not.invalid) invalid = size(pAinc,1) < im
  if(.not.invalid) invalid = size(pAinc,2) < jnp
  if(.not.invalid) invalid = size(uAinc,1) < im
  if(.not.invalid) invalid = size(uAinc,2) < jnp
  if(.not.invalid) invalid = size(vAinc,1) < im
  if(.not.invalid) invalid = size(vAinc,2) < jnp

  if(invalid) then
    if(size(rlat) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlat)',size(rlat)	)
    if(size(rlon) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlon)',size(rlon)	)
    if(size(rlev) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlev)',size(rlev)	)
    if(size(tims) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(tims)',size(tims)	)
    if(size(dels) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(dels)',size(dels)	)
    if(size(kxs ) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(kxs)',size(kxs)	)
    if(size(kts ) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(kts)',size(kts)	)
    if(size(xvec) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(xvec)',size(xvec)	)	

    if(size(pAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(pAinc,1)',size(pAinc,1)	)
    if(size(pAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(pAinc,2)',size(pAinc,2)	)

    if(size(uAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(uAinc,1)',size(uAinc,1)	)
    if(size(uAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(uAinc,2)',size(uAinc,2)	)

    if(size(vAinc,1) < im  ) call perr(myname_,		&
	'im'  ,im  ,'size(vAinc,1)',size(vAinc,1)	)
    if(size(vAinc,2) < jnp ) call perr(myname_,		&
	'jnp' ,jnp ,'size(vAinc,2)',size(vAinc,2)	)

    call die(myname_,'invalid argument(s)')
  endif

!________________________________________

	allocate( xmUs_(nobs), stat=ier	)
		if(ier/=0) call die(myname_,'allocate(obs)',ier)

	if(mpout_ison()) then
	  write(mpout,'(2a,i2)') myname_,': on entry'
	  call ObsSmry(mpout,nobs,kxs,kts)
	endif

!  Fill up arrays that normally would come from superObs
!  -----------------------------------------------------
  xmUs_(1:nobs) = 1.0  			! Scale sigO factor

  call AIGrid_init(ninc,im,jnp, size(KTpuv),KTpuv  )

	allocate( vinc(ninc),vlat(ninc),vlon(ninc),	&
		  vlev(ninc),vkts(ninc), stat=ier	)
		if(ier/=0) call die(myname_,'allocate(vec)',ier)

  call AIGrid_bind(ninc, vlat,vlon,vlev,vkts )

  call AE_solve(	&
	ninc, vinc, vlat, vlon,	vlev, vkts,	&
	nobs, xvec, rlat, rlon, rlev, tims,	&
	       kxs,  kss,  kts,xmUs_,dels	)

  call AIGrid_intp(ktslp,    im,jnp,pAinc,       ninc,vinc)
  call AIGrid_intp(ktUs,ktVs,im,jnp,uAinc,vAinc, ninc,vinc)

	deallocate(vinc,vlat,vlon,vlev,vkts, stat=ier )

		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

  call AIGrid_clean()

	deallocate(xmUs_,	stat=ier	)

		if(ier/=0) call die(myname_,'deallocate(obs)',ier)

  call lvstat(mpout,im,jnp,uAinc,0.,'WIND','HGHT',1.e+15,'USL')
  call lvstat(mpout,im,jnp,vAinc,0.,'WIND','HGHT',1.e+15,'VSL')
  call lvstat(mpout,im,jnp,pAinc,0.,'PRES','HGHT',1.e+15,'PSL')

end subroutine getAIpuv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getAIall_ - A multivariate analysis eqn. solver
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getAIall_(nobs,rlat,rlon,rlev,tims,kxs,kts,kss,dels,	&
	im,jnp,mlev,plevs,pslAinc,uslAinc,vslAinc,		&
			  zAinc,uAinc,vAinc,qAinc,xvec		)

      use m_die,only : die,perr
      use m_AIGrid,only : AIGrid_init
      use m_AIGrid,only : AIGrid_clean
      use m_AIGrid,only : AIGrid_bind
      use m_AIGrid,only : AIGrid_intp
      use m_AE,    only : AE_solve
      use m_ObsSmry,only : ObsSmry
      use m_mpout, only : mpout,mpout_ison

      implicit none
      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels (pressures)
      real   ,dimension(:),intent(in) :: tims	! time off-sets (min)
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      integer,dimension(:),intent(in) :: kss    ! sounding index
      real   ,dimension(:),intent(in) :: dels	! an O-F vector

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
      real,dimension(:),intent(out) :: xvec	 ! solution of conjgr

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
! 	09May00 - Todling - Removed superobbing
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getAIall_'

	! The workspace for a local copy of the O-F vector and its
	! artributes.

  real,allocatable,dimension(:) :: xmUs_

	! The workspace for a local copy of A-F vector and its
	! artributes.

  integer :: ninc
  real,allocatable,dimension(:) :: vinc
  real,allocatable,dimension(:) :: vlat
  real,allocatable,dimension(:) :: vlon
  real,allocatable,dimension(:) :: vlev
  integer,allocatable,dimension(:) :: vkts

  logical :: invalid
  integer :: inc
  integer :: ier
  integer :: i

  invalid=.false.
  if(.not.invalid) invalid = size(rlat) < nobs
  if(.not.invalid) invalid = size(rlon) < nobs
  if(.not.invalid) invalid = size(rlev) < nobs
  if(.not.invalid) invalid = size(tims) < nobs
  if(.not.invalid) invalid = size(dels) < nobs
  if(.not.invalid) invalid = size(xvec) < nobs
  if(.not.invalid) invalid = size(kxs ) < nobs
  if(.not.invalid) invalid = size(kts ) < nobs
  if(.not.invalid) invalid = size(kss ) < nobs

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
    if(size(rlat) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlat)',size(rlat)	)
    if(size(rlon) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlon)',size(rlon)	)
    if(size(rlev) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(rlev)',size(rlev)	)
    if(size(tims) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(tims)',size(tims)	)
    if(size(dels) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(dels)',size(dels)	)
    if(size(kxs ) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(kxs)',size(kxs)	)
    if(size(kts ) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(kts)',size(kts)	)
    if(size(xvec) < nobs) call perr(myname_,	&
	'nobs',nobs,'size(xvec)',size(xvec)	)	

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

    call die(myname_,'invalid argument(s)')
  endif

!________________________________________

	allocate( xmUs_(nobs), stat=ier	)
		if(ier/=0) call die(myname_,'allocate(obs)',ier)

	if(mpout_ison()) then
	  write(mpout,'(2a,i2)') myname_,': on entry'
	  call ObsSmry(mpout,nobs,kxs,kts)
	endif

!  Fill up arrays that normally would come from superObs
!  -----------------------------------------------------
  xmUs_(1:nobs) = 1.0  			! Scale sigO factor

  call AIGrid_init(ninc,im,jnp,mlev,plevs,size(KTall),KTall )

	allocate( vinc(ninc),vlat(ninc),vlon(ninc),	&
		  vlev(ninc),vkts(ninc), stat=ier	)
		if(ier/=0) call die(myname_,'allocate(vec)',ier)

  call AIGrid_bind(ninc, vlat,vlon,vlev,vkts )

  call AE_solve(	&
	ninc, vinc, vlat, vlon,	vlev, vkts,	&
	nobs, xvec, rlat, rlon, rlev, tims,	&
	       kxs,  kss,  kts,xmUs_,dels	)

  call AIGrid_intp(ktUs,ktVs,im,jnp,uslAinc,vslAinc, ninc,vinc)
  call AIGrid_intp(ktslp,    im,jnp,pslAinc,         ninc,vinc)

  call AIGrid_intp(ktUU,ktVV,im,jnp,mlev,plevs, uAinc,vAinc, ninc,vinc)
  call AIGrid_intp(ktHH,     im,jnp,mlev,plevs, zAinc,       ninc,vinc)

  call AIGrid_intp(ktQQ,     im,jnp,mlev,plevs, qAinc,       ninc,vinc)

	deallocate(vinc,vlat,vlon,vlev,vkts, stat=ier )

		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

  call AIGrid_clean()

	deallocate(xmUs_,	stat=ier	)

		if(ier/=0) call die(myname_,'deallocate(obs)',ier)

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

end subroutine getAIall_

end module m_AISolverNoSuperOb
