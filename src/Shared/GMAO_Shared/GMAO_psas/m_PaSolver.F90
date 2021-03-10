!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_PaSolver - MPI Pa times vector module
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_PaSolver
      use m_ktList,only : ktslp,ktus,ktvs,ktHH,ktUU,ktVV,ktQQ

      implicit none
      private	! except

      public :: PaAll
      public :: PfAll   ! should be in a separate file

      interface PaAll; module procedure	&
	solve_mix_, &
	solve_puv_, &
	solve_zuv_, &
	solve_all_
      end interface
      interface PfAll; module procedure	&
	solve_mix2_, &
	solve_puv2_, &
	solve_zuv2_, &
	solve_all2_
      end interface

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!		- provided framework for this implementation
!	15Mar02	- Todling  - implemented this MPI Pa operator
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_PaSolver'

  integer,dimension(1),parameter :: KTmix=(/ktQQ/)
  integer,dimension(3),parameter :: KTpuv=(/ktUs,ktVs,ktslp/)
  integer,dimension(3),parameter :: KTzuv=(/ktUU,ktVV,ktHH/)
  integer,dimension(7),parameter :: KTall=	&
	(/ktUs,ktVs,ktslp,ktUU,ktVV,ktHH,ktQQ/)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_mix_ - A univariate Pa operation on a vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_mix_(im,jnp,mlev,plevs,			&
	qAinc,							&
	nobs,rlat,rlon,rlev,kxs,kss,kts,xmUs,			&
	comm,root,psas_rsrc,eagrid_rsrc,what)

      use m_rootedLLGrid,only : rootedLLGrid
      use m_rootedLLGrid,only : rootedLLGrid_init
      use m_rootedLLGrid,only : ptr_lat
      use m_rootedLLGrid,only : ptr_lon
      use m_rootedLLGrid,only : ptr_lev
      use m_rootedLLGrid,only : ptr_kt
      use m_rootedLLGrid,only : distrSize
      use m_rootedLLGrid,only : clean
      use m_rootedLLGrid,only : scatterv
      use m_rootedLLGrid,only : gatherv

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat
      use m_PaE,   only : Pa_Cx
      use m_ObsSmry,only : ObsSmry

      implicit none

      integer,intent(in) :: im		! no. of longutide grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:,:),intent(out) :: qAinc ! gridded q analysis

      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc
      integer         ,optional,intent(in) :: what ! operator defined according to
                                                   !  what=-1 ==> op(Pf)
                                                   !  what= 0 ==> op(inc(Pa))
                                                   !  what= 1 ==> op(Pa) default

! !REVISION HISTORY:
!	15Mar02	- RT/JG - implemented this MPI Pa operator
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_mix_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedLLGrid) :: llg

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: inc
  integer :: ier

!________________________________________

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  call verifyObsVec_(nobs,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,myname_)

  if(myID==root) then
    invalid=.false.
    if(.not.invalid) invalid = size(plevs)  < mlev

    if(.not.invalid) invalid = size(qAinc  ,1) < im
    if(.not.invalid) invalid = size(qAinc  ,2) < jnp
    if(.not.invalid) invalid = size(qAinc,  3) < mlev

    if(invalid) then
      if(size(plevs)  < mlev) call perr(myname_,	&
	'mlev',mlev,'size(plevs)',size(plevs)	)

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

  call rootedLLGrid_init(llg,im,jnp,plevs,KTmix,root,comm)

	ninc=distrSize(llg)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call scatterv(ktQQ ,  qAinc,vinc,llg,root,comm)

  call Pa_Cx(ninc,vinc,						&
	ptr_lat(llg),ptr_lon(llg),ptr_lev(llg),ptr_kt(llg),	&
	nobs,rlat,rlon,rlev,kxs,kss,kts,xmUs,			&
	comm,root,rsrc=psas_rsrc,what=what)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,ktQQ,   qAinc,llg,root,comm)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

  call clean(llg)

  if(myID==root.and.mpout_ison()) then

    inc=1
    if(mlev > 0) then
      if(plevs(1)>plevs(mlev)) inc=-1
    endif

    call gdstat(mpout,im,jnp,mlev,qAinc,plevs,'MIXR','PRES',1.e+15,  &
	'Upperair water vapor mixing-ratio analysis increments',inc)
  endif

end subroutine solve_mix_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_puv_ - A multivariate 2-d Pa operation on a vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_puv_(im,jnp,mlev,plevs,			&
	pslAinc,uslAinc,vslAinc,				&
	nobs,rlat,rlon,rlev,kxs,kss,kts,xmUs,			&
	comm,root,psas_rsrc,eagrid_rsrc,what)

      use m_rootedLLGrid,only : rootedLLGrid
      use m_rootedLLGrid,only : rootedLLGrid_init
      use m_rootedLLGrid,only : ptr_lat
      use m_rootedLLGrid,only : ptr_lon
      use m_rootedLLGrid,only : ptr_lev
      use m_rootedLLGrid,only : ptr_kt
      use m_rootedLLGrid,only : distrSize
      use m_rootedLLGrid,only : clean
      use m_rootedLLGrid,only : scatterv
      use m_rootedLLGrid,only : gatherv

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat
      use m_PaE,   only : Pa_Cx
      use m_ObsSmry,only : ObsSmry

      implicit none

      integer,intent(in) :: im		! no. of longutide grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:),intent(out) :: pslAinc ! gridded p_sl analysis
      real,dimension(:,:),intent(out) :: uslAinc ! gridded u_sl analysis
      real,dimension(:,:),intent(out) :: vslAinc ! gridded v_sl analysis

      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc
      integer         ,optional,intent(in) :: what ! operator defined according to
                                                   !  what=-1 ==> op(Pf)
                                                   !  what= 0 ==> op(inc(Pa))
                                                   !  what= 1 ==> op(Pa) default

! !REVISION HISTORY:
!	15Mar02	- RT/JG - implemented this MPI Pa operator
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_puv_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedLLGrid) :: llg

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: inc
  integer :: ier

!________________________________________

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  call verifyObsVec_(nobs,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,myname_)

  if(myID==root) then
    invalid=.false.
    if(.not.invalid) invalid = size(plevs)  < mlev

    if(.not.invalid) invalid = size(pslAinc,1) < im
    if(.not.invalid) invalid = size(pslAinc,2) < jnp

    if(.not.invalid) invalid = size(uslAinc,1) < im
    if(.not.invalid) invalid = size(uslAinc,2) < jnp

    if(.not.invalid) invalid = size(vslAinc,1) < im
    if(.not.invalid) invalid = size(vslAinc,2) < jnp

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

      call die(myname_,'invalid argument(s)')
    endif

!	if(mpout_ison()) then
!	  call mpout_log(myname_,'A summary of the input observations:')
!	  call ObsSmry(mpout,nobs,kxs,kts)
!	endif
  endif
!________________________________________

	! Create a distributed AIGrid vector based on its specification
	! on the root PE.

  call rootedLLGrid_init(llg,im,jnp,plevs,KTpuv,root,comm)

	ninc=distrSize(llg)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call scatterv(ktSLP,pslAinc,vinc,llg,root,comm)
  call scatterv(ktUs ,uslAinc,vinc,llg,root,comm)
  call scatterv(ktVs ,vslAinc,vinc,llg,root,comm)

  call Pa_Cx(ninc,vinc,						&
	ptr_lat(llg),ptr_lon(llg),ptr_lev(llg),ptr_kt(llg),	&
	nobs,rlat,rlon,rlev,kxs,kss,kts,xmUs,			&
	comm,root,rsrc=psas_rsrc,what=what)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,ktSLP,pslAinc,llg,root,comm)
  call gatherv(vinc,ktUs ,uslAinc,llg,root,comm)
  call gatherv(vinc,ktVs ,vslAinc,llg,root,comm)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

  call clean(llg)

  if(myID==root.and.mpout_ison()) then

    call lvstat(mpout,im,jnp,uslAinc,0.,'WIND','HGHT',1.e+15,'U_SL')
    call lvstat(mpout,im,jnp,vslAinc,0.,'WIND','HGHT',1.e+15,'V_SL')
    call lvstat(mpout,im,jnp,pslAinc,0.,'PRES','HGHT',1.e+15,'P_SL')

  endif

end subroutine solve_puv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_zuv_ - A multivariate Pa operation on a vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_zuv_(im,jnp,mlev,plevs,			&
	zAinc,uAinc,vAinc,					&
	nobs,rlat,rlon,rlev,kxs,kss,kts,xmUs,			&
	comm,root,psas_rsrc,eagrid_rsrc,what)

      use m_rootedLLGrid,only : rootedLLGrid
      use m_rootedLLGrid,only : rootedLLGrid_init
      use m_rootedLLGrid,only : ptr_lat
      use m_rootedLLGrid,only : ptr_lon
      use m_rootedLLGrid,only : ptr_lev
      use m_rootedLLGrid,only : ptr_kt
      use m_rootedLLGrid,only : distrSize
      use m_rootedLLGrid,only : clean
      use m_rootedLLGrid,only : scatterv
      use m_rootedLLGrid,only : gatherv

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat
      use m_PaE,   only : Pa_Cx
      use m_ObsSmry,only : ObsSmry

      implicit none

      integer,intent(in) :: im		! no. of longutide grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:,:),intent(out) :: zAinc ! gridded z analysis
      real,dimension(:,:,:),intent(out) :: uAinc ! gridded u analysis
      real,dimension(:,:,:),intent(out) :: vAinc ! gridded v analysis

      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc
      integer         ,optional,intent(in) :: what ! operator defined according to
                                                   !  what=-1 ==> op(Pf)
                                                   !  what= 0 ==> op(inc(Pa))
                                                   !  what= 1 ==> op(Pa) default

! !REVISION HISTORY:
!	15Mar02	- RT/JG - implemented this MPI Pa operator
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_zuv_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedLLGrid) :: llg

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: inc
  integer :: ier

!________________________________________

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  call verifyObsVec_(nobs,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,myname_)

  if(myID==root) then
    invalid=.false.
    if(.not.invalid) invalid = size(plevs)  < mlev

    if(.not.invalid) invalid = size(zAinc  ,1) < im
    if(.not.invalid) invalid = size(zAinc  ,2) < jnp
    if(.not.invalid) invalid = size(zAinc  ,3) < mlev

    if(.not.invalid) invalid = size(uAinc  ,1) < im
    if(.not.invalid) invalid = size(uAinc  ,2) < jnp
    if(.not.invalid) invalid = size(uAinc  ,3) < mlev

    if(.not.invalid) invalid = size(vAinc  ,1) < im
    if(.not.invalid) invalid = size(vAinc  ,2) < jnp
    if(.not.invalid) invalid = size(vAinc  ,3) < mlev

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

!	if(mpout_ison()) then
!	  call mpout_log(myname_,'A summary of the input observations:')
!	  call ObsSmry(mpout,nobs,kxs,kts)
!	endif
  endif
!________________________________________

	! Create a distributed AIGrid vector based on its specification
	! on the root PE.

  call rootedLLGrid_init(llg,im,jnp,plevs,KTzuv,root,comm)

	ninc=distrSize(llg)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call scatterv(ktHH ,  zAinc,vinc,llg,root,comm)
  call scatterv(ktUU ,  uAinc,vinc,llg,root,comm)
  call scatterv(ktVV ,  vAinc,vinc,llg,root,comm)

  call Pa_Cx(ninc,vinc,						&
	ptr_lat(llg),ptr_lon(llg),ptr_lev(llg),ptr_kt(llg),	&
	nobs,rlat,rlon,rlev,kxs,kss,kts,xmUs,			&
	comm,root,rsrc=psas_rsrc,what=what)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,ktHH,   zAinc,llg,root,comm)
  call gatherv(vinc,ktUU,   uAinc,llg,root,comm)
  call gatherv(vinc,ktVV,   vAinc,llg,root,comm)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

  call clean(llg)

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

end subroutine solve_zuv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_all_ - A multivariate Pa operation on a vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_all_(im,jnp,mlev,plevs,			&
	pslAinc,uslAinc,vslAinc,zAinc,uAinc,vAinc,qAinc,	&
	nobs,rlat,rlon,rlev,kxs,kss,kts,xmUs,			&
	comm,root,psas_rsrc,eagrid_rsrc,what)

      use m_rootedLLGrid,only : rootedLLGrid
      use m_rootedLLGrid,only : rootedLLGrid_init
      use m_rootedLLGrid,only : ptr_lat
      use m_rootedLLGrid,only : ptr_lon
      use m_rootedLLGrid,only : ptr_lev
      use m_rootedLLGrid,only : ptr_kt
      use m_rootedLLGrid,only : distrSize
      use m_rootedLLGrid,only : clean
      use m_rootedLLGrid,only : scatterv
      use m_rootedLLGrid,only : gatherv

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat
      use m_PaE,   only : Pa_Cx
      use m_ObsSmry,only : ObsSmry

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
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc
      integer         ,optional,intent(in) :: what ! operator defined according to
                                                   !  what=-1 ==> op(Pf)
                                                   !  what= 0 ==> op(inc(Pa))
                                                   !  what= 1 ==> op(Pa) default

! !REVISION HISTORY:
!	15Mar02	- RT/JG - implemented this MPI Pa operator
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_all_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedLLGrid) :: llg

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: inc
  integer :: ier

!________________________________________

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  call verifyObsVec_(nobs,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,myname_)

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

!	if(mpout_ison()) then
!	  call mpout_log(myname_,'A summary of the input observations:')
!	  call ObsSmry(mpout,nobs,kxs,kts)
!	endif
  endif
!________________________________________

	! Create a distributed AIGrid vector based on its specification
	! on the root PE.

  call rootedLLGrid_init(llg,im,jnp,plevs,KTall,root,comm)

	ninc=distrSize(llg)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  call scatterv(ktHH ,  zAinc,vinc,llg,root,comm)
  call scatterv(ktUU ,  uAinc,vinc,llg,root,comm)
  call scatterv(ktVV ,  vAinc,vinc,llg,root,comm)
  call scatterv(ktQQ ,  qAinc,vinc,llg,root,comm)
  call scatterv(ktSLP,pslAinc,vinc,llg,root,comm)
  call scatterv(ktUs ,uslAinc,vinc,llg,root,comm)
  call scatterv(ktVs ,vslAinc,vinc,llg,root,comm)

  call Pa_Cx(ninc,vinc,						&
	ptr_lat(llg),ptr_lon(llg),ptr_lev(llg),ptr_kt(llg),	&
	nobs,rlat,rlon,rlev,kxs,kss,kts,xmUs,			&
	comm,root,rsrc=psas_rsrc,what=what)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,ktHH,   zAinc,llg,root,comm)
  call gatherv(vinc,ktUU,   uAinc,llg,root,comm)
  call gatherv(vinc,ktVV,   vAinc,llg,root,comm)
  call gatherv(vinc,ktQQ,   qAinc,llg,root,comm)
  call gatherv(vinc,ktSLP,pslAinc,llg,root,comm)
  call gatherv(vinc,ktUs ,uslAinc,llg,root,comm)
  call gatherv(vinc,ktVs ,vslAinc,llg,root,comm)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vec)',ier)

  call clean(llg)

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
! !IROUTINE: solve_mix2_ - A univariate Pf operation on a vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_mix2_(im,jnp,mlev,plevs,			&
	qAinc,							&
	comm,root,psas_rsrc,eagrid_rsrc)

      use m_rootedLLGrid,only : rootedLLGrid
      use m_rootedLLGrid,only : rootedLLGrid_init
      use m_rootedLLGrid,only : distrSize
      use m_rootedLLGrid,only : ptr_lat
      use m_rootedLLGrid,only : ptr_lon
      use m_rootedLLGrid,only : ptr_lev
      use m_rootedLLGrid,only : ptr_kt
      use m_rootedLLGrid,only : scatterv
      use m_rootedLLGrid,only : gatherv
      use m_rootedLLGrid,only : clean

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat
      use m_PaE,   only : Pf_Cx
      use m_ObsSmry,only : ObsSmry

      implicit none

      integer,intent(in) :: im		! no. of longutide grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:,:),intent(out) :: qAinc ! gridded q analysis

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc

! !REVISION HISTORY:
!	15Mar02	- RT/JG - implemented this MPI Pa operator
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_mix2_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedLLGrid) :: llg

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: inc
  integer :: ier

!________________________________________

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  if(myID==root) then
    invalid=.false.
    if(.not.invalid) invalid = size(plevs)  < mlev

    if(.not.invalid) invalid = size(qAinc  ,1) < im
    if(.not.invalid) invalid = size(qAinc  ,2) < jnp
    if(.not.invalid) invalid = size(qAinc,  3) < mlev

    if(invalid) then
      if(size(plevs)  < mlev) call perr(myname_,	&
	'mlev',mlev,'size(plevs)',size(plevs)	)

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

  call rootedLLGrid_init(llg,im,jnp,plevs,KTmix,root,comm)

	ninc=distrSize(llg)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

		! Cache distributed vinc from the root PE.

  call scatterv(ktQQ ,  qAinc,vinc,llg,root,comm)

  call Pf_Cx(ninc,vinc,						&
	ptr_lat(llg),ptr_lon(llg),ptr_lev(llg),ptr_kt(llg),	&
	comm,root,rsrc=psas_rsrc)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,ktQQ,   qAinc,llg,root,comm)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vinc)',ier)

  call clean(llg)

  if(myID==root.and.mpout_ison()) then

    inc=1
    if(mlev > 0) then
      if(plevs(1)>plevs(mlev)) inc=-1
    endif

    call gdstat(mpout,im,jnp,mlev,qAinc,plevs,'MIXR','PRES',1.e+15,  &
	'Upperair water vapor mixing-ratio analysis increments',inc)
  endif

end subroutine solve_mix2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_puv2_ - A multivariate 2-d Pf operation on a vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_puv2_(im,jnp,mlev,plevs,			&
	pslAinc,uslAinc,vslAinc,				&
	comm,root,psas_rsrc,eagrid_rsrc)

      use m_rootedLLGrid,only : rootedLLGrid
      use m_rootedLLGrid,only : rootedLLGrid_init
      use m_rootedLLGrid,only : distrSize
      use m_rootedLLGrid,only : ptr_lat
      use m_rootedLLGrid,only : ptr_lon
      use m_rootedLLGrid,only : ptr_lev
      use m_rootedLLGrid,only : ptr_kt
      use m_rootedLLGrid,only : scatterv
      use m_rootedLLGrid,only : gatherv
      use m_rootedLLGrid,only : clean

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat
      use m_PaE,   only : Pf_Cx
      use m_ObsSmry,only : ObsSmry

      implicit none

      integer,intent(in) :: im		! no. of longutide grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:),intent(out) :: pslAinc ! gridded p_sl analysis
      real,dimension(:,:),intent(out) :: uslAinc ! gridded u_sl analysis
      real,dimension(:,:),intent(out) :: vslAinc ! gridded v_sl analysis

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc

! !REVISION HISTORY:
!	15Mar02	- RT/JG - implemented this MPI Pa operator
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_all2_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedLLGrid) :: llg

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: inc
  integer :: ier

!________________________________________

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  if(myID==root) then
    invalid=.false.
    if(.not.invalid) invalid = size(plevs)  < mlev

    if(.not.invalid) invalid = size(pslAinc,1) < im
    if(.not.invalid) invalid = size(pslAinc,2) < jnp

    if(.not.invalid) invalid = size(uslAinc,1) < im
    if(.not.invalid) invalid = size(uslAinc,2) < jnp

    if(.not.invalid) invalid = size(vslAinc,1) < im
    if(.not.invalid) invalid = size(vslAinc,2) < jnp

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

      call die(myname_,'invalid argument(s)')
    endif

  endif
!________________________________________

	! Create a distributed AIGrid vector based on its specification
	! on the root PE.

  call rootedLLGrid_init(llg,im,jnp,plevs,KTpuv,root,comm)

	ninc=distrSize(llg)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

		! Cache distributed vinc from the root PE.

  call scatterv(ktSLP,pslAinc,vinc,llg,root,comm)
  call scatterv(ktUs ,uslAinc,vinc,llg,root,comm)
  call scatterv(ktVs ,vslAinc,vinc,llg,root,comm)

  call Pf_Cx(ninc,vinc,						&
	ptr_lat(llg),ptr_lon(llg),ptr_lev(llg),ptr_kt(llg),	&
	comm,root,rsrc=psas_rsrc)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,ktSLP,pslAinc,llg,root,comm)
  call gatherv(vinc,ktUs ,uslAinc,llg,root,comm)
  call gatherv(vinc,ktVs ,vslAinc,llg,root,comm)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vinc)',ier)

  call clean(llg)

  if(myID==root.and.mpout_ison()) then

    call lvstat(mpout,im,jnp,uslAinc,0.,'WIND','HGHT',1.e+15,'U_SL')
    call lvstat(mpout,im,jnp,vslAinc,0.,'WIND','HGHT',1.e+15,'V_SL')
    call lvstat(mpout,im,jnp,pslAinc,0.,'PRES','HGHT',1.e+15,'P_SL')

  endif

end subroutine solve_puv2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_zuv2_ - A multivariate (upper-air) Pf operation on a vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_zuv2_(im,jnp,mlev,plevs,			&
	zAinc,uAinc,vAinc,					&
	comm,root,psas_rsrc,eagrid_rsrc)

      use m_rootedLLGrid,only : rootedLLGrid
      use m_rootedLLGrid,only : rootedLLGrid_init
      use m_rootedLLGrid,only : distrSize
      use m_rootedLLGrid,only : ptr_lat
      use m_rootedLLGrid,only : ptr_lon
      use m_rootedLLGrid,only : ptr_lev
      use m_rootedLLGrid,only : ptr_kt
      use m_rootedLLGrid,only : scatterv
      use m_rootedLLGrid,only : gatherv
      use m_rootedLLGrid,only : clean

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat
      use m_PaE,   only : Pf_Cx
      use m_ObsSmry,only : ObsSmry

      implicit none

      integer,intent(in) :: im		! no. of longutide grid points
      integer,intent(in) :: jnp		! no. of latitude grid points
      integer,intent(in) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	 ! level grid points
      real,dimension(:,:,:),intent(out) :: zAinc ! gridded z analysis
      real,dimension(:,:,:),intent(out) :: uAinc ! gridded u analysis
      real,dimension(:,:,:),intent(out) :: vAinc ! gridded v analysis

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc

! !REVISION HISTORY:
!	15Mar02	- RT/JG - implemented this MPI Pa operator
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_zuv2_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedLLGrid) :: llg

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: inc
  integer :: ier

!________________________________________

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  if(myID==root) then
    invalid=.false.
    if(.not.invalid) invalid = size(plevs)  < mlev

    if(.not.invalid) invalid = size(zAinc  ,1) < im
    if(.not.invalid) invalid = size(zAinc  ,2) < jnp
    if(.not.invalid) invalid = size(zAinc  ,3) < mlev

    if(.not.invalid) invalid = size(uAinc  ,1) < im
    if(.not.invalid) invalid = size(uAinc  ,2) < jnp
    if(.not.invalid) invalid = size(uAinc  ,3) < mlev

    if(.not.invalid) invalid = size(vAinc  ,1) < im
    if(.not.invalid) invalid = size(vAinc  ,2) < jnp
    if(.not.invalid) invalid = size(vAinc  ,3) < mlev

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

  call rootedLLGrid_init(llg,im,jnp,plevs,KTzuv,root,comm)

	ninc=distrSize(llg)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

		! Cache distributed vinc from the root PE.

  call scatterv(ktHH ,  zAinc,vinc,llg,root,comm)
  call scatterv(ktUU ,  uAinc,vinc,llg,root,comm)
  call scatterv(ktVV ,  vAinc,vinc,llg,root,comm)

  call Pf_Cx(ninc,vinc,						&
	ptr_lat(llg),ptr_lon(llg),ptr_lev(llg),ptr_kt(llg),	&
	comm,root,rsrc=psas_rsrc)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,ktHH,   zAinc,llg,root,comm)
  call gatherv(vinc,ktUU,   uAinc,llg,root,comm)
  call gatherv(vinc,ktVV,   vAinc,llg,root,comm)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vinc)',ier)

  call clean(llg)

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

end subroutine solve_zuv2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_all2_ - A multivariate Pf operation on a vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_all2_(im,jnp,mlev,plevs,			&
	pslAinc,uslAinc,vslAinc,zAinc,uAinc,vAinc,qAinc,	&
	comm,root,psas_rsrc,eagrid_rsrc)

      use m_rootedLLGrid,only : rootedLLGrid
      use m_rootedLLGrid,only : rootedLLGrid_init
      use m_rootedLLGrid,only : distrSize
      use m_rootedLLGrid,only : ptr_lat
      use m_rootedLLGrid,only : ptr_lon
      use m_rootedLLGrid,only : ptr_lev
      use m_rootedLLGrid,only : ptr_kt
      use m_rootedLLGrid,only : scatterv
      use m_rootedLLGrid,only : gatherv
      use m_rootedLLGrid,only : clean

      use m_die,only : die,perr,MP_die
      use m_mpif90,only : MP_comm_rank
      use m_mpout, only : mpout,mpout_ison,mpout_log
      use m_gdstat,only : gdstat
      use m_gdstat,only : lvstat
      use m_PaE,   only : Pf_Cx
      use m_ObsSmry,only : ObsSmry

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

      integer,intent(in) :: comm	! communicator
      integer,intent(in) :: root	! where is the root PE

				! A resource file for m_AE
      character(len=*),optional,intent(in) :: psas_rsrc
				! A resource file for m_EAGrid
      character(len=*),optional,intent(in) :: eagrid_rsrc

! !REVISION HISTORY:
!	15Mar02	- RT/JG - implemented this MPI Pa operator
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_all2_'

	! The workspace for a local copy of A-F vector and its
	! artributes.

  type(rootedLLGrid) :: llg

  integer :: ninc
  real,allocatable,dimension(:) :: vinc

  logical :: invalid
  integer :: myID
  integer :: inc
  integer :: ier

!________________________________________

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

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

  call rootedLLGrid_init(llg,im,jnp,plevs,KTall,root,comm)

	ninc=distrSize(llg)
	allocate(vinc(ninc),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

		! Cache distributed vinc from the root PE.

  call scatterv(ktHH ,  zAinc,vinc,llg,root,comm)
  call scatterv(ktUU ,  uAinc,vinc,llg,root,comm)
  call scatterv(ktVV ,  vAinc,vinc,llg,root,comm)
  call scatterv(ktQQ ,  qAinc,vinc,llg,root,comm)
  call scatterv(ktSLP,pslAinc,vinc,llg,root,comm)
  call scatterv(ktUs ,uslAinc,vinc,llg,root,comm)
  call scatterv(ktVs ,vslAinc,vinc,llg,root,comm)

  call Pf_Cx(ninc,vinc,						&
	ptr_lat(llg),ptr_lon(llg),ptr_lev(llg),ptr_kt(llg),	&
	comm,root,rsrc=psas_rsrc)

		! Cache distributed vinc on the root PE.

  call gatherv(vinc,ktHH,   zAinc,llg,root,comm)
  call gatherv(vinc,ktUU,   uAinc,llg,root,comm)
  call gatherv(vinc,ktVV,   vAinc,llg,root,comm)
  call gatherv(vinc,ktQQ,   qAinc,llg,root,comm)
  call gatherv(vinc,ktSLP,pslAinc,llg,root,comm)
  call gatherv(vinc,ktUs ,uslAinc,llg,root,comm)
  call gatherv(vinc,ktVs ,vslAinc,llg,root,comm)

	deallocate(vinc,stat=ier )
		if(ier/=0) call die(myname_,'deallocate(vinc)',ier)

  call clean(llg)

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

end subroutine solve_all2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: verifyObsVec_ - verify obs. vec. dimensions
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine verifyObsVec_(nobs,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs, where)

      use m_die,only : die,perr
      implicit none

      integer,intent(in) :: nobs	! size of the O-F vector
      real   ,dimension(:),intent(in) :: rlat	! latitudes
      real   ,dimension(:),intent(in) :: rlon	! longitudes
      real   ,dimension(:),intent(in) :: rlev	! levels
      integer,dimension(:),intent(in) :: kxs	! instrument IDs
      integer,dimension(:),intent(in) :: kss	! sounding IDs
      integer,dimension(:),intent(in) :: kts	! variable types
      real   ,dimension(:),intent(in) :: xmUs	! "meta" scalar of sigU

      character(len=*),intent(in) :: where

! !REVISION HISTORY:
! 	26Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::verifyObsVec_'

  logical :: invalid
  integer :: ier

  invalid=.false.

  if(.not.invalid) invalid = size(rlat) < nobs
  if(.not.invalid) invalid = size(rlon) < nobs
  if(.not.invalid) invalid = size(rlev) < nobs
  if(.not.invalid) invalid = size(kxs ) < nobs
  if(.not.invalid) invalid = size(kss ) < nobs
  if(.not.invalid) invalid = size(kts ) < nobs
  if(.not.invalid) invalid = size(xmUs) < nobs

  if(invalid) then
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

    call die(where,'invalid argument(s)')
  endif
end subroutine verifyObsVec_

end module m_PaSolver
