!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AIGrid - Analysis grid vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_AIGrid
      use m_LVGrid,only : LVGrid
      use m_EAGrid,only : EAGrid
      implicit none
      private	! except

      public :: AIGrid_init
      public :: AIGrid_bind
      public :: AIGrid_intp
      public :: AIGrid_clean

      interface AIGrid_init; module procedure	&
	init__,	&
	init2__
      end interface
      interface AIGrid_clean; module procedure clean__; end interface

      interface AIGrid_bind; module procedure bind_; end interface
      interface AIGrid_intp; module procedure	&
	intps_,		&
	intpv_,		&
	intps2_,	&
	intpv2_
      end interface

! !REVISION HISTORY:
! 	03Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AIGrid'

  type(LVGrid),save :: AIGrid_lvG
  type(EAGrid),save :: AIGrid_eaG

  logical,save :: AIGrid_defined=.false.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init__ - initialized a grid description (internal)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init__(ninc,nlon,nlat,nlev,plevs,nkt,ktList,rsrc)
      use m_LVGrid,only : LVGrid_init
      use m_LVGrid,only : LVGrid_size
      use m_EAGrid,only : EAGrid_init
      use m_die,   only : die
      implicit none
      integer,intent(out) :: ninc	! the size of the AIGrid vector
      integer,intent(in ) :: nlon	! no. of longitude grid points
      integer,intent(in ) :: nlat	! no. of latitude grid points
      integer,intent(in ) :: nlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	! level grid points
      integer,intent(in ) :: nkt	! no. of variable types
      integer,dimension(:),intent(in) :: ktList	! variable types
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	03Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init__'

  if(AIGrid_defined) call die(myname_,'multiple object definition')

  call LVGrid_init(AIGrid_lvG,nlev,plevs,nkt,ktList)
  call EAGrid_init(AIGrid_eaG,nlon,nlat,rsrc=rsrc)

  ninc=LVGrid_size(AIGrid_lvG)*AIGrid_eaG%nEA
	
  AIGrid_defined=.true.
end subroutine init__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init2__ - initialized a grid description (internal)
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init2__(ninc,nlon,nlat,nkt,ktList,rsrc)
      use m_LVGrid,only : LVGrid_init
      use m_LVGrid,only : LVGrid_size
      use m_EAGrid,only : EAGrid_init
      use m_die,   only : die
      implicit none
      integer,intent(out) :: ninc	! the size of the AIGrid vector
      integer,intent(in ) :: nlon	! no. of longitude grid points
      integer,intent(in ) :: nlat	! no. of latitude grid points
      integer,intent(in ) :: nkt	! no. of variable types
      integer,dimension(:),intent(in) :: ktList	! variable types
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	03Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init2__'

  if(AIGrid_defined) call die(myname_,'multiple object definition')

  call LVGrid_init(AIGrid_lvG,nkt,ktList)
  call EAGrid_init(AIGrid_eaG,nlon,nlat,rsrc=rsrc)

  ninc=LVGrid_size(AIGrid_lvG)*AIGrid_eaG%nEA
	
  AIGrid_defined=.true.
end subroutine init2__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bind_ - return attributes of the AIGrid vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine bind_(ninc, vlat,vlon,vlev,vkts )
      use m_die, only : die
      use m_LVGrid,only : LVGrid_size
      use m_LVGrid,only : ptr_kt,ptr_kp,ptr_levs
      implicit none
      integer,intent(in) :: ninc	! size of the AIGrid vector
      real,   dimension(:),intent(out) :: vlat	! latitudes
      real,   dimension(:),intent(out) :: vlon	! longitudes
      real,   dimension(:),intent(out) :: vlev	! levels
      integer,dimension(:),intent(out) :: vkts	! variable types

! !REVISION HISTORY:
! 	03Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bind_'

  if(.not.AIGrid_defined) call die(myname_,'object undefined')

  if(ninc /= LVGrid_size(AIGrid_lvG)*AIGrid_eaG%nEA )	&
	call die(myname_,'ninc',ninc,'nLV*nEA',	&
		LVGrid_size(AIGrid_lvG)*AIGrid_eaG%nEA)

	! should have more error checking here

		! Would this shield improve optimization?

  call bindEALV_(vlat,vlon,vlev,vkts,	&
	AIGrid_eaG%nEA,AIGrid_eaG%lat,AIGrid_eaG%lon,	&
	LVGrid_size(AIGrid_lvG),ptr_kt(AIGrid_lvG),	&
	ptr_kp(AIGrid_lvG),ptr_levs(AIGrid_lvG))

contains

subroutine bindEALV_(vlat,vlon,vlev,vkts,	&
	nEA,glat,glon,nLV,kt,kp,levs)
  implicit none
  real,   dimension(:),intent(out) :: vlat	! latitudes
  real,   dimension(:),intent(out) :: vlon	! longitudes
  real,   dimension(:),intent(out) :: vlev	! levels
  integer,dimension(:),intent(out) :: vkts	! variable types
  integer,intent(in) :: nEA
  real,   dimension(:),intent(in)  :: glat,glon
  integer,intent(in) :: nLV
  integer,dimension(:),intent(in)  :: kt,kp
  real,   dimension(:),intent(in)  :: levs

	! Local variables
  integer :: i,iv,ih
  integer :: kti,kpi
  real    :: rlev

	! The vectors are ordered as (1:nEA,1:nLV)

  do iv=1,nLV
    kti=kt(iv)
    kpi=kp(iv)
    rlev=-1.			! This could be a level-less grid
    if(kpi>0) rlev=levs(kpi)

    i=(iv-1)*nEA
    do ih=1,nEA
      vlat(i+ih)=glat(ih)
      vlon(i+ih)=glon(ih)

      vlev(i+ih)=rlev
      vkts(i+ih)=kti
    end do
  end do

end subroutine bindEALV_
end subroutine bind_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intps_ - interpolate a scalar field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intps_(kt,im,jnp,mlev,plevs,zinc, ninc,vinc )
      use m_LevCache,only : LevCache_init
      use m_LevCache,only : LevCache_locate
      use m_LevCache,only : LevCache_clean
      use m_LVGrid,  only : LVGrid_index
      use m_LVGrid,  only : ptr_levs
      use m_die,     only : die
      use m_mall,    only : mall_ison,mall_mci,mall_mco
      implicit none
      integer,intent(in ) :: kt	! the variable type to be interpolated
      integer,intent(in ) :: im		! no. of longitude grid points
      integer,intent(in ) :: jnp	! no. of latitude grid points
      integer,intent(in ) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	! level grid points
      real,dimension(:,:,:),intent(out) :: zinc ! interpolated field
      integer,intent(in) :: ninc	! size of the AIGrid vector
      real,dimension(:),intent(in) :: vinc ! a scattered field

! !REVISION HISTORY:
! 	03Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intps_'

  logical,parameter :: nearest=.true.

  real,pointer,dimension(:,:) :: zm,zp
  real,pointer,dimension(:) :: p_levs
  integer :: m_levs

  real   ,allocatable,dimension(:) :: wlev
  integer,allocatable,dimension(:) :: klev

  integer :: i,j,k,kvm,kvp
  integer :: ier
  logical :: found
  integer :: llev
  real    :: wght

  if(.not.AIGrid_defined) call die(myname_,'object undefined')

  call LevCache_init(im,jnp,nVar=1)

	allocate(klev(mlev),wlev(mlev),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

		if(mall_ison()) then
		  call mall_mci(klev,myname)
		  call mall_mci(wlev,myname)
		endif

	p_levs => ptr_levs(AIGrid_lvG)
	m_levs = size(p_levs)

  call slogtab(.not.nearest, m_levs,p_levs,mlev,plevs,klev,wlev)

	nullify(p_levs)

  do k=1,mlev
    kvm=LVGrid_index(AIGrid_lvG,kt,klev(k))

	! zm will be pointing to the level with kvm
    call LevCache_locate(kvm,found,zm)
    if(.not.found) call intpslev_(zm,AIGrid_eaG,kvm,vinc)

    do j=1,jnp
      do i=1,im
	zinc(i,j,k)=zm(i,j)
      end do
    end do

    if( klev(k) < m_levs ) then
      llev=klev(k)+1
      wght=wlev(k)

      kvp=LVGrid_index(AIGrid_lvG,kt,llev)

      call LevCache_locate(kvp,found,zp)
      if(.not.found) call intpslev_(zp,AIGrid_eaG,kvp,vinc)

      do j=1,jnp
      do i=1,im
	zinc(i,j,k)=zinc(i,j,k)+wght*(zp(i,j)-zinc(i,j,k))
      end do
      end do
    endif
  end do

		if(mall_ison()) then
		  call mall_mco(klev,myname)
		  call mall_mco(wlev,myname)
		endif
	deallocate(klev,wlev,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

  call LevCache_clean()
end subroutine intps_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intpv_ - interpolate from a vector field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intpv_(ktu,ktv,im,jnp,mlev,plevs,uAinc,vAinc, ninc,vinc)
      use m_LevCache,only : LevCache_init
      use m_LevCache,only : LevCache_locate
      use m_LevCache,only : LevCache_clean
      use m_LVGrid,  only : LVGrid_index
      use m_LVGrid,  only : ptr_levs
      use m_die,     only : die
      use m_mall,    only : mall_ison,mall_mci,mall_mco
      implicit none
      integer,intent(in ) :: ktu,ktv ! variables to be interpolated
      integer,intent(in ) :: im		! no. of longitude grid points
      integer,intent(in ) :: jnp	! no. of latitude grid points
      integer,intent(in ) :: mlev	! no. of level grid points
      real,dimension(:),intent(in) :: plevs	! level grid points
      real,dimension(:,:,:),intent(out) :: uAinc ! interpolated field
      real,dimension(:,:,:),intent(out) :: vAinc ! interpolated field
      integer,intent(in) :: ninc	! size of the AIGrid vector
      real,dimension(:),intent(in) :: vinc ! a scattered field

! !REVISION HISTORY:
! 	03Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intpv_'

  logical,parameter :: nearest=.true.

  real,pointer,dimension(:,:) :: um,up
  real,pointer,dimension(:,:) :: vm,vp
  real,pointer,dimension(:) :: p_levs
  integer :: m_levs

  real   ,allocatable,dimension(:) :: wlev
  integer,allocatable,dimension(:) :: klev

  integer :: i,j,k,kum,kup,kvm,kvp
  integer :: ier
  logical :: found
  integer :: llev
  real    :: wght

  if(.not.AIGrid_defined) call die(myname_,'object undefined')

  call LevCache_init(im,jnp,nVar=4)

	allocate(klev(mlev),wlev(mlev),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

		if(mall_ison()) then
		  call mall_mci(klev,myname)
		  call mall_mci(wlev,myname)
		endif

	p_levs => ptr_levs(AIGrid_lvG)
	m_levs = size(p_levs)

  call slogtab(.not.nearest, m_levs,p_levs,mlev,plevs,klev,wlev)

	nullify(p_levs)

  do k=1,mlev
    kum=LVGrid_index(AIGrid_lvG,ktu,klev(k))
    kvm=LVGrid_index(AIGrid_lvG,ktv,klev(k))

    call LevCache_locate(kum,kvm,found,um,vm)
    if(.not.found) call intpvlev_(um,vm,AIGrid_eaG,kum,kvm,vinc)

    do j=1,jnp
    do i=1,im
      uAinc(i,j,k)=um(i,j)
      vAinc(i,j,k)=vm(i,j)
    end do
    end do

    if( klev(k) < m_levs ) then
      llev=klev(k)+1
      wght=wlev(k)

      kup=LVGrid_index(AIGrid_lvG,ktu,llev)
      kvp=LVGrid_index(AIGrid_lvG,ktv,llev)

      call LevCache_locate(kup,kvp,found,up,vp)
      if(.not.found) call intpvlev_(up,vp,AIGrid_eaG,kup,kvp,vinc)

      do j=1,jnp
      do i=1,im
	uAinc(i,j,k)=uAinc(i,j,k)+wght*(up(i,j)-uAinc(i,j,k))
	vAinc(i,j,k)=vAinc(i,j,k)+wght*(vp(i,j)-vAinc(i,j,k))
      end do
      end do
    endif
  end do

		if(mall_ison()) then
		  call mall_mco(klev,myname)
		  call mall_mco(wlev,myname)
		endif
	deallocate(klev,wlev,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

  call LevCache_clean()
end subroutine intpv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intps2_ - interpolate a scalar field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intps2_(kt,im,jnp,zinc, ninc,vinc )
      use m_LVGrid,  only : LVGrid_index
      use m_die, only : perr,die
      implicit none
      integer,intent(in ) :: kt	! the variable type to be interpolated
      integer,intent(in ) :: im		! no. of longitude grid points
      integer,intent(in ) :: jnp	! no. of latitude grid points
      real,dimension(:,:),intent(out) :: zinc ! interpolated field
      integer,intent(in) :: ninc	! size of the AIGrid vector
      real,dimension(:),intent(in) :: vinc ! a scattered field

! !REVISION HISTORY:
! 	03Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intps2_'

  integer :: kv

  if(.not.AIGrid_defined) call die(myname_,'object undefined')

    kv=LVGrid_index(AIGrid_lvG,kt)
    if(kv<=0) call die(myname_,'LVGrid_index(kt)',kt)

    if(im>size(zinc,1) .or. jnp>size(zinc,2)) then

      if(  im>size(zinc,1))	&
	call perr(myname_, 'im', im,'size(zinc,1)',size(zinc,1))
      if( jnp>size(zinc,2))	&
	call perr(myname_,'jnp',jnp,'size(zinc,2)',size(zinc,2))

      call die(myname_)
    endif

    call intpslev_(zinc(1:im,1:jnp),AIGrid_eaG,kv,vinc)

end subroutine intps2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intpv2_ - interpolate from a vector field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intpv2_(ktu,ktv,im,jnp,uAinc,vAinc, ninc,vinc)
      use m_LVGrid,  only : LVGrid_index
      use m_die, only : perr,die
      implicit none
      integer,intent(in ) :: ktu,ktv ! variables to be interpolated
      integer,intent(in ) :: im		! no. of longitude grid points
      integer,intent(in ) :: jnp	! no. of latitude grid points
      real,dimension(:,:),intent(out) :: uAinc ! interpolated field
      real,dimension(:,:),intent(out) :: vAinc ! interpolated field
      integer,intent(in) :: ninc	! size of the AIGrid vector
      real,dimension(:),intent(in) :: vinc ! a scattered field

! !REVISION HISTORY:
! 	03Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intpv2_'

  integer :: ku,kv

  if(.not.AIGrid_defined) call die(myname_,'object undefined')

  ku=LVGrid_index(AIGrid_lvG,ktu)
  kv=LVGrid_index(AIGrid_lvG,ktv)

  if(ku<=0 .or. kv<=0) then
    if(ku<=0) call perr(myname_,'LVGrid_index(ktu)',ktu)
    if(kv<=0) call perr(myname_,'LVGrid_index(ktv)',ktv)
    call die(myname_)
  endif

  if(	 im>size(uAinc,1) .or.	 im>size(vAinc,1) .or.	&
	jnp>size(uAinc,2) .or.	jnp>size(vAinc,2)	) then

    if(	 im>size(uAinc,1))	&
	call perr(myname_, 'im', im,'size(uAinc,1)',size(uAinc,1))
    if(	 im>size(vAinc,1))	&
	call perr(myname_, 'im', im,'size(vAinc,1)',size(vAinc,1))
    if( jnp>size(uAinc,2))	&
	call perr(myname_,'jnp',jnp,'size(uAinc,2)',size(uAinc,2))
    if( jnp>size(vAinc,2))	&
	call perr(myname_,'jnp',jnp,'size(vAinc,2)',size(vAinc,2))
    call die(myname_)
  endif

  call intpvlev_( uAinc(1:im,1:jnp),	&
		  vAinc(1:im,1:jnp),	&
		  AIGrid_eaG,ku,kv,vinc	)

end subroutine intpv2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intpslev_ - interpolate a scalar field to a 2-d grid
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intpslev_(zll,eaG,kv,vinc)
      use m_EAGrid,only : EAGrid,EAGrid_ea2ll
      implicit none
      real,dimension(:,:),intent(out) :: zll	! a lon-lat grid
      type(EAGrid),intent(in) :: eaG	! an EAGrid
      integer,intent(in) :: kv		! index of the level
      real,dimension(:),intent(in) :: vinc	! an all-var. vector

! !REVISION HISTORY:
! 	05Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intpslev_'
  integer :: lc,le
  integer :: nx,ny

  nx=size(zll,1)
  ny=size(zll,2)

  le=eaG%nEA*kv
  lc=le-eaG%nEA+1

  call EAGrid_ea2ll(zll,nx,ny,vinc(lc:le),eaG%nEA)

end subroutine intpslev_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intpvlev_ interplate a spherical vector field
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intpvlev_(ull,vll,eaG,ku,kv,vinc)
      use m_EAGrid,only : EAGrid,EAGrid_ea2ll
      use m_die,only : perr,die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none
      real,dimension(:,:),intent(out) :: ull	! u lon-lat grid
      real,dimension(:,:),intent(out) :: vll	! v lon-lat grid
      type(EAGrid),intent(in) :: eaG	! an EAGrid
      integer,intent(in) :: ku		! index of u level
      integer,intent(in) :: kv		! index of v level
      real,dimension(:),intent(in) :: vinc	! an all-var. vector

! !REVISION HISTORY:
! 	05Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intpvlev_'
  integer :: nx,ny,nEA
  integer :: lcu,leu,lcv,lev
  integer :: ier

  real,allocatable,dimension(:)   :: wx,wy,wz
  real,allocatable,dimension(:,:) :: vx,vy,vz

  nx=size(ull,1)
  ny=size(ull,2)
  if(nx/=size(vll,1) .or. ny/=size(vll,2)) then
    if(nx/=size(vll,1)) call perr(myname_,		&
	'size(ull,1)',nx,'size(vll,1)',size(vll,1)	)
    if(ny/=size(vll,2)) call perr(myname_,		&
	'size(ull,2)',ny,'size(vll,2)',size(vll,2)	)
    call die(myname_)
  endif

  nEA=eaG%nEA

  leu=nEA*ku
  lcu=leu-nEA+1
  lev=nEA*kv
  lcv=lev-nEA+1

	! This is a vector interpolation algorithm. It should
	! conserve the vector magnitude along the spherical flow.
	!
	! Compare to the simple scalar interpolation algorithm:
	!
	!   call EAGrid_ea2ll(ull,nx,ny,vinc(lcu:leu),nEA)
	!   call EAGrid_ea2ll(vll,nx,ny,vinc(lcv:lev),nEA)
	!
	! or the 3-d vector interpolation algorithm that treats
	! a spherical flow (u,v) field as an ordinary 3-d vector.

	allocate(wx(nEA),  wy(nEA),  wz(nEA),	&
		 vx(nx,ny),vy(nx,ny),vz(nx,ny),	stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

		if(mall_ison()) then
		  call mall_mci(wx,myname)
		  call mall_mci(wy,myname)
		  call mall_mci(wz,myname)
		  call mall_mci(vx,myname)
		  call mall_mci(vy,myname)
		  call mall_mci(vz,myname)
		endif

	! Conver linear vectors to unit angular momentums

  call l2avec_(nEA,eaG%lat,eaG%lon,wx,wy,wz,	&
	vinc(lcu:leu),vinc(lcv:lev)		)

  call EAGrid_ea2ll(vx,nx,ny,wx,nEA)
  call EAGrid_ea2ll(vy,nx,ny,wy,nEA)
  call EAGrid_ea2ll(vz,nx,ny,wz,nEA)

	! Conver unit angular momentums to linear vectors

  call a2lvec_(nx,ny,ull,vll,vx,vy,vz)

		if(mall_ison()) then
		  call mall_mco(wx,myname)
		  call mall_mco(wy,myname)
		  call mall_mco(wz,myname)
		  call mall_mco(vx,myname)
		  call mall_mco(vy,myname)
		  call mall_mco(vz,myname)
		endif
	deallocate(wx,wy,wz,vx,vy,vz,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)
end subroutine intpvlev_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: l2avec_ - convert linear vectors (u,v) to r|.cross.vectors.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine l2avec_(n,vlat,vlon,wx,wy,wz,uVec,vVec)
      implicit none
      integer,intent(in) :: n
      real,dimension(:),intent(in)  :: vlat,vlon
      real,dimension(:),intent(out) :: wx,wy,wz
      real,dimension(:),intent(in)  :: uVec,vVec

! !REVISION HISTORY:
! 	09Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::l2avec_'

  real :: deg		! radiance per degree
  real :: rlat,clat,slat
  real :: rlon,clon,slon
  real :: emx,emy,emz
  real :: elx,ely
  integer :: i

  deg=4.*atan(1.)/180.		! pi/180.

  do i=1,n
    rlat=vlat(i)*deg	! convert to radiance
    rlon=vlon(i)*deg

    clat=cos(rlat)
    slat=sin(rlat)
    clon=cos(rlon)
    slon=sin(rlon)

    emx=-slat*clon
    emy=-slat*slon
    emz= clat
    elx=-slon
    ely= clon

	! Let
	!
	!   w| = er| .cross. (u*el| + v*em|)
	!
	! one has
	!
	!   w|= u*em| - v*el|

    wx(i)=uVec(i)*emx - vVec(i)*elx
    wy(i)=uVec(i)*emy - vVec(i)*ely
    wz(i)=uVec(i)*emz
  end do

end subroutine l2avec_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: a2lvec_ - convert r|.cross.vectors to linear vectors (u,v).
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine a2lvec_(im,jnp,uVec,vVec,wx,wy,wz)
      implicit none
      integer,intent(in) :: im,jnp
      real,dimension(:,:),intent(out) :: uVec,vVec
      real,dimension(:,:),intent(in)  :: wx,wy,wz

! !REVISION HISTORY:
! 	09Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::a2lvec_'

  real,parameter :: SPlat= -90.,NPlat= 90.
  real,parameter :: WSlon=-180.,ESlon=180.

  real :: deg		! radiance per degree
  real :: dlat,dlon
  real :: rlat,clat,slat
  real :: rlon,clon,slon
  real :: emx,emy,emz
  real :: elx,ely
  integer :: i,j

  deg=4.*atan(1.)/180.		! pi/180.

  dlat=(NPlat-SPlat)/(jnp-1)
  dlon=(ESlon-WSlon)/im

  do j=1,jnp
    rlat=((j-1)*dlat + SPlat)*deg
    clat=cos(rlat)
    slat=sin(rlat)

    do i=1,im
      rlon=((i-1)*dlon + WSlon)*deg
      clon=cos(rlon)
      slon=sin(rlon)

      emx=-slat*clon
      emy=-slat*slon
      emz= clat
      elx=-slon
      ely= clon

	! From the definition of w| in l2avec_(), one has
	!
	!   u = w| .dot. em|
	!
	! and
	!
	!   v = w| .dot. el|

      uVec(i,j)=  wx(i,j)*emx + wy(i,j)*emy + wz(i,j)*emz
      vVec(i,j)=-(wx(i,j)*elx + wy(i,j)*ely)
    end do
  end do

end subroutine a2lvec_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean__ - clean all stored Q.E.A grid information
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean__()
      use m_LVGrid,only : LVGrid_clean
      use m_EAGrid,only : EAGrid_clean
      use m_die,   only : die
      implicit none

! !REVISION HISTORY:
! 	03Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean__'

  if(.not.AIGrid_defined) call die(myname_,'object undefined')

  call LVGrid_clean(AIGrid_lvG)
  call EAGrid_clean(AIGrid_eaG)

  AIGrid_defined=.false.
end subroutine clean__
end module m_AIGrid
