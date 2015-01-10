!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_rootedLLGrid - an implementation of a distributed LLGrid
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_rootedLLGrid
      use m_LVGrid      ,only : LVGrid
      use m_Collector   ,only : Collector
      implicit none
      private	! except

      public :: rootedLLGrid		! The class data structure
      public :: rootedLLGrid_init,init	! initialize an object
      public :: clean			! clean a distributed grid

      public :: distrSize	! size of the distrbuted vector
      public :: scatterv	! create a scattered vector
      public :: gatherv		! create a gathered vector

      public :: ptr_lat		! Component referencing, latitudes
      public :: ptr_lon		! Component referencing, longitudes
      public :: ptr_lev		! Component referencing, levels
      public :: ptr_kt		! Component referencing, var-types

      type rootedLLGrid
	private

		! In this implementation, a distributed vector is
		! stored in the following order (in case one needs to
		! know):
		!
		!   1) full lev-kt (v-dim),
		!   2) distributed lon-lat (h-dim).
		!
		! "v-dim" is the faster running index than "h-dim"

		! -%lvG- speicifies "v-dim", in the ordered of
		!
		!   a) pressure levels (lev), if kt is a 3-d variable.
		!   b) variable types (kt),
		!
		! "lev" is the faster index than "kt"
		!
		! See m_LVGrid for additional details

	type(LVGrid),pointer :: lvG

		! -%im- and -%jnp- specifies the information of "h-dim'
		! as a global grid.  -%icount- and -%idispl- specifies
		! the distribution information.  A h-dim is ordered
		! globally as
		!
		!   a) longitudes (lon).
		!   b) latitudes (lat),
		!
		! "lon" is the faster index than "lat"

	integer :: im
	integer :: jnp
	type(Collector),pointer :: coll

		! In a notation similar to the one for Fortran arrays,
		! the order may be denoted in the following form:
		!
		!   vect[[lev,kt],[lon,lat]]

		! Attributes of all local vector elements

	real   ,pointer,dimension(:) :: lat
	real   ,pointer,dimension(:) :: lon

	integer,pointer,dimension(:) :: kt
	real   ,pointer,dimension(:) :: lev

      end type rootedLLGrid

      interface rootedLLGrid_init; module procedure	&
	init_,		&
	init2_; end interface

      interface init; module procedure	&
	init_,		&
	init2_; end interface

      interface clean   ; module procedure clean_   ; end interface
      interface distrSize; module procedure dsize_  ; end interface
      interface gatherv ; module procedure	&
	gatherv2_,	&
	gatherv3_ ; end interface
      interface scatterv; module procedure	&
	scatterv2_,	&
	scatterv3_; end interface

      interface ptr_lat; module procedure ptr_lat_; end interface
      interface ptr_lon; module procedure ptr_lon_; end interface
      interface ptr_lev; module procedure ptr_lev_; end interface
      interface ptr_kt ; module procedure ptr_kt_ ; end interface

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_rootedLLGrid'
#include "assert.H"

  real,parameter :: WLON=-180.
  real,parameter :: SLAT=-90.

! Usecase:
!
!	use m_rootedLLGrid,only : rootedLLGrid
!	use m_rootedLLGrid,only : rootedLLGrid_init
!	use m_rootedLLGrid,only : distrSize
!	use m_rootedLLGrid,only : clean
!	use m_rootedLLGrid,only : scatterv
!	use m_rootedLLGrid,only : gatherv
!	use m_ktList,only : ktHH,ktUU,ktVV,ktQQ
!	use m_ktList,only : ktUs,ktVs,ktslp
!	use m_xyz,only : solve
!	[...]
!
!	implicit none
!
!		! Although must defined on all PEs, grid dimension data
!		! are significant only on the root PE.,
!
!	integer,parameter :: im=144,jnp=91,mlev=6
!	real,parameter,dimension(mlev) ::	&
!		plevs = (/ 850.,700.,500.,300.,200.,100./)
!	integer,parameter,dimension(7) ::	&
!		ktList= (/ ktHH,ktUU,ktVV,ktQQ,ktslp,ktUs,ktVs /)
!
!	integer,parameter :: root=0
!	integer :: comm
!
!	type(rootedLLGrid) :: llGrid 
!
!		! Grid variables are expected to be dimensioned
!		! as (im,jnp,mlev) on the root PE, defined but any
!		! dimension on other PEs.
!
!	real,allocatable,dimension(:,:,:) :: hght,uwnd,vwnd,mixr
!	real,allocatable,dimension(:,:)   :: slp ,slu ,slv
!
!		! Vector variables are explicitly dimensioned in the
!		! code below.
!
!	real,allocatable,dimension(:) :: vectin
!	real,allocatable,dimension(:) :: vectout
!
!	integer :: lvec
!
!		! Grid variables must be properly allocated and defined
!		! on all PEs before passed as the input to scatterv().
!		!
!		! Grid variables must be properly allocated on all PEs
!		! before used to store the output from gatherv().
!
!	[...]
!
!	comm=MPI_comm_world
!
!   call rootedLLGrid_init(llGrid,im,jnp,plevs,ktList,root,comm)
!
!	lvec=distrSize(llGrid)
!	allocate(vectin(lvec))
!
!   call scatterv(ktHH,hght,vectin,llGrid,root,comm)
!   call scatterv(ktUU,uwnd,vectin,llGrid,root,comm)
!   call scatterv(ktVV,vwnd,vectin,llGrid,root,comm)
!   call scatterv(ktQQ,mixr,vectin,llGrid,root,comm)
!   call scatterv(ktSLP,slp,vectin,llGrid,root,comm)
!   call scatterv(ktUs ,slu,vectin,llGrid,root,comm)
!   call scatterv(ktVs ,slv,vectin,llGrid,root,comm)
!
!	allocate(vectout(lvec))
!
!	call solve(ptr_lat(llGrid),ptr_lon(llGrid),	&
!		   plt_lev(llGrid),ptr_kts(llGrid),	&
!		   vectin,vectout,comm)
!
!	deallocate(vectin(lvec))
!
!   call gatherv(vectout,ktHH,hght,llGrid,root,comm)
!   call gatherv(vectout,ktUU,uwnd,llGrid,root,comm)
!   call gatherv(vectout,ktVV,vwnd,llGrid,root,comm)
!   call gatherv(vectout,ktQQ,mixr,llGrid,root,comm)
!   call gatherv(vectout,ktSLP,slp,llGrid,root,comm)
!   call gatherv(vectout,ktUs ,slu,llGrid,root,comm)
!   call gatherv(vectout,ktVs ,slv,llGrid,root,comm)
!
!	deallocate(vectout(lvec))
!
!   call clean(llGrid)
!

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a LLGrid with scattered attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(obj,im,jnp,plevs,KTList,root,comm)
      use m_LVGrid    ,only : new,LVGrid_init
      use m_LVGrid    ,only : LVGrid_size
      use m_LVGrid    ,only : ptr_kt
      use m_LVGrid    ,only : ptr_kp
      use m_LVGrid    ,only : ptr_levs
      use m_simplePart,only : simplePart
      use m_Collector ,only : new,Collector_init
      use m_Collector ,only : get

      use m_mpif90,only : MP_comm_rank
      use m_mpif90,only : MP_comm_size
      use m_mpif90,only : MP_type

      use m_die ,only : MP_die,die
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none

      type(rootedLLGrid),intent(out) :: obj	! output arguments

		! input arguments used only on root

      integer,intent(in) :: im
      integer,intent(in) :: jnp
      real   ,dimension(:),intent(in) :: plevs
      integer,dimension(:),intent(in) :: KTList

      integer,intent(in) :: root
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

  integer :: ier
  integer :: myID,nPEs
  integer :: kt,kp
  integer :: lc,le,icount,idispl
  integer :: im_,jnp_,mlev,nkt
  integer :: l,i,j,lv,ij,lij,nLV,lsize
  real    :: dlon,dlat,pr,rlonij,rlatij

  integer,dimension(4) :: ibuf
  integer,allocatable,dimension(:) :: kbuf
  real   ,allocatable,dimension(:) :: rbuf

  real   ,pointer,dimension(:) :: p_pr
  integer,pointer,dimension(:) :: p_kt,p_kp
!________________________________________

  call MP_comm_size(comm,nPEs,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_size()',ier)

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)
!________________________________________

  if(myID==root) then
    ibuf(1)=im
    ibuf(2)=jnp
    ibuf(3)=size(plevs)
    ibuf(4)=size(ktList)
  endif

  call MPI_bcast(ibuf,size(ibuf),MP_type(ibuf),root,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_bcast(ibuf)',ier)

	! Localize the variable names

  im_ =ibuf(1)
  jnp_=ibuf(2)
  mlev=ibuf(3)
  nkt =ibuf(4)
!________________________________________

	! Define "v-dim"

	allocate(kbuf(nkt),rbuf(mlev),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(kbuf,myname_)
		  call mall_mci(rbuf,myname_)
		endif
		 
  if(myID==root) kbuf(:)=ktList(:)
  call MPI_bcast(kbuf,size(kbuf),MP_type(kbuf),root,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_bcast(kbuf)',ier)

  if(myID==root) rbuf(:)=plevs(:)
  call MPI_bcast(rbuf,size(rbuf),MP_type(rbuf),root,comm,ier)
	if(ier/=0) call MP_die(myname_,'MPI_bcast(rbuf)',ier)

	obj%lvG => new(obj%lvG)
  call LVGrid_init(obj%lvG,mlev,rbuf,nkt,kbuf)

		if(mall_ison()) then
		  call mall_mco(kbuf,myname_)
		  call mall_mco(rbuf,myname_)
		endif
	deallocate(kbuf,rbuf,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)
!________________________________________

	! Define "h-dim"

  obj%im  =im_
  obj%jnp =jnp_

	! Determine the size of a local section in a given kt-lev grid
	! level.  Note that the size of a local section (icount) is the
	! same for all kt-lev grid level on a given PE.  However, this
	! size is not the same on different PEs.

  call simplePart(obj%im*obj%jnp,nPEs,myID,count=icount,displ=idispl)

	obj%coll => new(obj%coll)
  call Collector_init(obj%coll,icount,comm)
!________________________________________

	! Define local vector attributes

  dlon=360./max(1,obj%im)
  dlat=180./max(1,obj%jnp-1)

	! local vector storage size for all sections

  nLV=LVGrid_size(obj%lvG)
  lsize=icount*nLV

	allocate(obj%lat(lsize),obj%lon(lsize),	&
		 obj%kt (lsize),obj%lev(lsize), stat=ier)

		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(obj%lat,myname)
		  call mall_mci(obj%lon,myname)
		  call mall_mci(obj%kt ,myname)
		  call mall_mci(obj%lev,myname)
		endif

	p_kt => ptr_kt(obj%lvG)		! (1:nLV)
	p_kp => ptr_kp(obj%lvG)		! (1:nLV)
	p_pr => ptr_levs(obj%lvG)	! (1:mlev)

  call get(obj%coll,lbound=lc,ubound=le)

  do ij=lc,le

		! Find out (i,j) for a local section.  Note lat-lon
		! grid values are the same for all local vector
		! sections.

    i=1+mod(ij-1 ,obj%im)
    j=1+   (ij-1)/obj%im

    rlonij=(i-1)*dlon + WLON
    rlatij=(j-1)*dlat + SLAT

    lij=(ij-lc)*nLV
    do lv=1,nLV

		! Find out the vector storage range.  Note that the
		! size of a (lv) section is the same as in both %coll
		! and %navi

      kt=p_kt(lv)
      kp=p_kp(lv)
      pr=-1
      if(kp>0)pr=p_pr(kp)

      obj%lon(lij+lv)=rlonij
      obj%lat(lij+lv)=rlatij
      obj%lev(lij+lv)=pr
      obj%kt (lij+lv)=kt
    end do
  end do

	nullify(p_kt)
	nullify(p_kp)
	nullify(p_pr)

end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init2_ - initialize a 2-d LLGrid with scattered attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init2_(obj,im,jnp,KTList,root,comm)
      use m_mpif90,only : MP_comm_rank
      use m_die ,only : MP_die,die
      implicit none

      type(rootedLLGrid),intent(out) :: obj	! output arguments

		! input arguments used only on root

      integer,intent(in) :: im
      integer,intent(in) :: jnp
      integer,dimension(:),intent(in) :: KTList

      integer,intent(in) :: root
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init2_'
  real,dimension(0) :: plevs

  call init_(obj,im,jnp,plevs,KTList,root,comm)

end subroutine init2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gatherv3_ - gather vinc from all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gatherv3_(vect,kt,grid,obj,root,comm)
      use m_Collector    ,only : localSize
      use m_CollectorComm,only : gatherv

      use m_LVGrid,only : LVGrid_range
      use m_LVGrid,only : LVGrid_size

      use m_mpif90,only : MP_comm_rank

      use m_die ,only : MP_die,die,assert_
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none

      real,dimension(:),intent(in) :: vect	! distributed grid
      integer          ,intent(in) :: kt	! var-type of the grid
      real,dimension(:,:,:),intent(out) :: grid	! output grid
      type(rootedLLGrid),intent(in) :: obj

      integer,intent(in) :: root	! use an root
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gatherv3_'

  real,allocatable,dimension(:,:) :: bvect
  real,allocatable,dimension(:,:) :: bgrid
  integer :: nLV,lbnd,ubnd,mlev
  integer :: nIJ,mIJ
  integer :: myID
  integer :: l,lv,i,j
  integer :: ier
!________________________________________

! Algorithm:
!
!   Given a distributed vector, sending a block of the vector for a
! specified variable type (kt) to -root- for a collected (gathered)
! grid requires following steps:
!
!   1) Buffering local sections of specified kt sections of the
! distributed vector (1-d to 2-d copying),
!
!   2) Collecting the distributed buffers to the root PE for a copy of
! the full grid (2-d to 2-d gatherv()),
!
!   3) Un-buffering the collected buffer to the output storage (2-d to
! 3-d copying and transpose).
!
!   Note that for each gatherv() call, only partial data of the whole
! input vector -vect- is used to produce the output vector -grid-.
!________________________________________

  nIJ=localSize(obj%coll)		! local size of "h-dim"

  nLV=LVGrid_size(obj%lvG)		! total "v-dim" size
  call LVGrid_range(obj%lvG,kt,lbnd,ubnd) ! "v-dim" grid range for -kt-
  mlev=ubnd-lbnd+1			! no. of "v-dim" grid points

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

#ifndef NDEBUG
			! Verify the input argument

	ASSERT(nIJ*nLV==size(vect))

			! Verify the output argument
  if(myID==root) then
	ASSERT(obj%im ==size(grid,1))
	ASSERT(obj%jnp==size(grid,2))
	ASSERT(mlev   ==size(grid,3))
  endif
#endif

  if(mlev==0) return
!_______________________________________________________________________

		! Determine h-dim local size of the vector rooted at
		! -root-.  i.e. mIJ is zero if this is not the root PE,
		! is im*jnp if this is the root PE.
	mIJ=0
	if(myID==root) mIJ=obj%im*obj%jnp

		! -bvect- and -bgrid- are message send/receive buffers.
		! They are 2-d arrays, such that a common -%coll- can
		! be defined and used for variable types (kt) with
		! different number of levels (mlev).

	allocate(bvect(mlev,nIJ),bgrid(mlev,mIJ),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(bvect,myname_)
		  call mall_mci(bgrid,myname_)
		endif
!________________________________________

	! Create a send-buffer one colume at a time.

  do l=1,nIJ		! for all h-dim points
			! for vect(:), starting from
    lv=(l-1)*nLV
			! Copy a colume in v-dim

    bvect(:,l)=vect(lv+lbnd:lv+ubnd)
  end do
!________________________________________

	! Collect distributed buffers to root PE.

  call gatherv(bvect,bgrid,obj%coll,root,comm)
!________________________________________

	! Un-buffer to fill a tranposed output gridded data.  Note
	! that "mIJ>0" implies "myID==root".

  do l=1,mIJ
    i=1+mod(l-1 ,obj%im)
    j=1+   (l-1)/obj%im

		! Copy a column in v-dim

    grid(i,j,:)=bgrid(:,l)
  end do
!________________________________________

		if(mall_ison()) then
		  call mall_mco(bvect,myname_)
		  call mall_mco(bgrid,myname_)
		endif
	deallocate(bvect,bgrid,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine gatherv3_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gatherv2_ - gather vinc from all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gatherv2_(vect,kt,grid,obj,root,comm)
      use m_Collector    ,only : localSize
      use m_CollectorComm,only : gatherv

      use m_LVGrid,only : LVGrid_range
      use m_LVGrid,only : LVGrid_size

      use m_mpif90,only : MP_comm_rank

      use m_die ,only : MP_die,die,assert_
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none

      real,dimension(:),intent(in) :: vect	! distributed grid
      integer          ,intent(in) :: kt	! var-type of the grid
      real,dimension(:,:),intent(out) :: grid	! output grid
      type(rootedLLGrid),intent(in) :: obj

      integer,intent(in) :: root	! use an root
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gatherv2_'

  real,allocatable,dimension(:) :: bvect
  real,allocatable,dimension(:) :: bgrid
  integer :: nLV,lbnd,ubnd,mlev
  integer :: nIJ,mIJ
  integer :: myID
  integer :: l,lv,i,j
  integer :: ier
!________________________________________

! Algorithm:
!
!   Given a distributed vector, sending a block of the vector for a
! specified variable type (kt) to -root- for a collected (gathered)
! grid requires following steps:
!
!   1) Buffering local sections of specified kt sections of the
! distributed vector (1-d to 1-d copying),
!
!   2) Collecting the distributed buffers to the root PE for a copy of
! the full grid (1-d to 1-d gatherv()),
!
!   3) Un-buffering the collected buffer to the output storage (1-d to
! 2-d copying).
!
!   Note that for each gatherv() call, only partial data of the whole
! input vector -vect- is used to produce the output vector -grid-.
!________________________________________

  nIJ=localSize(obj%coll)		! local size of "h-dim"

  nLV=LVGrid_size(obj%lvG)		! total "v-dim" size
  call LVGrid_range(obj%lvG,kt,lbnd,ubnd) ! "v-dim" grid range for -kt-
  mlev=ubnd-lbnd+1			! no. of "v-dim" grid points

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

#ifndef NDEBUG
			! Verify the input argument

	ASSERT(nIJ*nLV==size(vect))
	ASSERT(mlev   ==1)

			! Verify the output argument
  if(myID==root) then
	ASSERT(obj%im ==size(grid,1))
	ASSERT(obj%jnp==size(grid,2))
  endif
#endif

!_______________________________________________________________________

		! Determine h-dim local size of the vector rooted at
		! -root-.  i.e. mIJ is zero if this is not the root PE,
		! is im*jnp if this is the root PE.
	mIJ=0
	if(myID==root) mIJ=obj%im*obj%jnp

		! -bvect- and -bgrid- are message send/receive buffers.
		! They are 2-d arrays, such that a common -%coll- can
		! be defined and used for variable types (kt) with
		! different number of levels (mlev).

	allocate(bvect(nIJ),bgrid(mIJ),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(bvect,myname_)
		  call mall_mci(bgrid,myname_)
		endif
!________________________________________

	! Create a send-buffer for all h-dim points

  ubnd=(nIJ-1)*nLV+lbnd
  bvect(1:nIJ)=vect(lbnd:ubnd:nLV)
!________________________________________

	! Collect distributed buffers to root PE.

  call gatherv(bvect,bgrid,obj%coll,root,comm)
!________________________________________

	! Un-buffer to fill a tranposed output gridded data.  Note
	! that "mIJ>0" implies "myID==root".

  do l=0,mIJ-1,obj%im
    j=l/obj%im+1
    grid(:,j)=bgrid(l+1:l+obj%im)
  end do
!________________________________________

		if(mall_ison()) then
		  call mall_mco(bvect,myname_)
		  call mall_mco(bgrid,myname_)
		endif
	deallocate(bvect,bgrid,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine gatherv2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scatterv3_ - scatter vinc from all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine scatterv3_(kt,grid,vect,obj,root,comm)
      use m_Collector    ,only : localSize
      use m_CollectorComm,only : scatterv

      use m_LVGrid,only : LVGrid_range
      use m_LVGrid,only : LVGrid_size

      use m_mpif90,only : MP_comm_rank

      use m_die ,only : MP_die,die,assert_
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none

      integer              ,intent(in) :: kt	! var-type of the grid
      real,dimension(:,:,:),intent(in) :: grid	! output grid
      real,dimension(:),intent(inout)  :: vect	! distributed grid
      type(rootedLLGrid),intent(in) :: obj

      integer,intent(in) :: root	! use an root
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scatterv3_'

  real,allocatable,dimension(:,:) :: bvect
  real,allocatable,dimension(:,:) :: bgrid
  integer :: nLV,lbnd,ubnd,mlev
  integer :: nIJ,mIJ
  integer :: myID
  integer :: l,lv,i,j
  integer :: ier
!________________________________________

! Algorithm:
!
!   This is an inverse operation of gatherv() above.  Note that for
! each scatterv() call, only partial data of the whole output vector
! -vect- is assigned by the input vector -grid-.
!________________________________________

  nIJ=localSize(obj%coll)		! local size of "h-dim"

  nLV=LVGrid_size(obj%lvG)		! total "v-dim" size
  call LVGrid_range(obj%lvG,kt,lbnd,ubnd) ! "v-dim" grid range for -kt-
  mlev=ubnd-lbnd+1			! no. of "v-dim" grid points

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

#ifndef NDEBUG
			! Verify the output argument

	ASSERT(nIJ*nLV==size(vect))

			! Verify the input argument
  if(myID==root) then
	ASSERT(obj%im ==size(grid,1))
	ASSERT(obj%jnp==size(grid,2))
	ASSERT(mlev   ==size(grid,3))
  endif
#endif

  if(mlev==0) return
!_______________________________________________________________________

		! Determine h-dim local size of the vector rooted at
		! -root-.  i.e. mIJ is zero if this is not the root PE,
		! is im*jnp if this is the root PE.
	mIJ=0
	if(myID==root) mIJ=obj%im*obj%jnp

		! -bvect- and -bgrid- are message send/receive buffers.
		! They are 2-d arrays, such that a common -%coll- can
		! be defined and used for variable types (kt) with
		! different number of levels (mlev).

	allocate(bvect(mlev,nIJ),bgrid(mlev,mIJ),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(bvect,myname_)
		  call mall_mci(bgrid,myname_)
		endif
!________________________________________

	! Create a send-buffer one colume at a time.  Note that
	! "mIJ>0" implies "myID==root"

  do l=1,mIJ
    i=1+mod(l-1 ,obj%im)
    j=1+   (l-1)/obj%im

		! Copy a column in v-dim

    bgrid(:,l)=grid(i,j,:)
  end do
!________________________________________

	! Scatter buffers on root PE, to all PEs.

  call scatterv(bgrid,bvect,obj%coll,root,comm)
!________________________________________

	! Un-buffer to selectively fill an output vector.

  do l=1,nIJ		! for all h-dim points
			! for vect(:), starting from
    lv=(l-1)*nLV
			! Copy a colume in v-dim

    vect(lv+lbnd:lv+ubnd)=bvect(:,l)
  end do
!________________________________________

		if(mall_ison()) then
		  call mall_mco(bvect,myname_)
		  call mall_mco(bgrid,myname_)
		endif
	deallocate(bvect,bgrid,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine scatterv3_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scatterv2_ - scatter vinc from all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine scatterv2_(kt,grid,vect,obj,root,comm)
      use m_Collector    ,only : localSize
      use m_CollectorComm,only : scatterv

      use m_LVGrid,only : LVGrid_range
      use m_LVGrid,only : LVGrid_size

      use m_mpif90,only : MP_comm_rank

      use m_die ,only : MP_die,die,assert_
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none

      integer              ,intent(in) :: kt	! var-type of the grid
      real,dimension(:,:)  ,intent(in) :: grid	! output grid
      real,dimension(:) ,intent(inout) :: vect	! distributed grid
      type(rootedLLGrid)   ,intent(in) :: obj

      integer,intent(in) :: root	! use an root
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scatterv2_'

  real,allocatable,dimension(:) :: bvect
  real,allocatable,dimension(:) :: bgrid
  integer :: nLV,lbnd,ubnd,mlev
  integer :: nIJ,mIJ
  integer :: myID
  integer :: l,lv,i,j
  integer :: ier
!________________________________________

! Algorithm:
!
!   This is an inverse operation of gatherv() above.  Note that for
! each scatterv() call, only partial data of the whole output vector
! -vect- is assigned by the input vector -grid-.
!________________________________________

  nIJ=localSize(obj%coll)		! local size of "h-dim"

  nLV=LVGrid_size(obj%lvG)		! total "v-dim" size
  call LVGrid_range(obj%lvG,kt,lbnd,ubnd) ! "v-dim" grid range for -kt-
  mlev=ubnd-lbnd+1			! no. of "v-dim" grid points

  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

#ifndef NDEBUG
			! Verify the output argument

	ASSERT(nIJ*nLV==size(vect))
	ASSERT(mlev   ==1)

			! Verify the input argument
  if(myID==root) then
	ASSERT(obj%im ==size(grid,1))
	ASSERT(obj%jnp==size(grid,2))
  endif
#endif

!_______________________________________________________________________

		! Determine h-dim local size of the vector rooted at
		! -root-.  i.e. mIJ is zero if this is not the root PE,
		! is im*jnp if this is the root PE.
	mIJ=0
	if(myID==root) mIJ=obj%im*obj%jnp

		! -bvect- and -bgrid- are message send/receive buffers.
		! They are 2-d arrays, such that a common -%coll- can
		! be defined and used for variable types (kt) with
		! different number of levels (mlev).

	allocate(bvect(nIJ),bgrid(mIJ),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(bvect,myname_)
		  call mall_mci(bgrid,myname_)
		endif
!________________________________________

	! Create a send-buffer one row at a time.  Note that "mIJ>0"
	! implies "myID==root"

  do l=0,mIJ-1,obj%im
    j=l/obj%im+1
    bgrid(l+1:l+obj%im)=grid(:,j)
  end do
!________________________________________

	! Scatter buffers on root PE, to all PEs.

  call scatterv(bgrid,bvect,obj%coll,root,comm)
!________________________________________

	! Un-buffer to selectively fill the output vector.

  ubnd=(nIJ-1)*nLV+lbnd
  vect(lbnd:ubnd:nLV)=bvect(1:nIJ)
!________________________________________

		if(mall_ison()) then
		  call mall_mco(bvect,myname_)
		  call mall_mco(bgrid,myname_)
		endif
	deallocate(bvect,bgrid,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine scatterv2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(obj)
      use m_LVGrid   ,only : clean => LVGrid_clean
      use m_LVGrid   ,only : delete
      use m_Collector,only : clean,delete
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mco
      implicit none
      type(rootedLLGrid),intent(inout) :: obj

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  call clean(obj%lvG)
  call delete(obj%lvG)

  call clean(obj%coll)
  call delete(obj%coll)

  obj%im=0
  obj%jnp=0

	if(mall_ison()) then
	  call mall_mco(obj%lat,myname)
	  call mall_mco(obj%lon,myname)
	  call mall_mco(obj%kt ,myname)
	  call mall_mco(obj%lat,myname)
	endif

  deallocate(obj%lat,obj%lon,obj%kt,obj%lev,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_lat_ - referencing latitudes
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_lat_(obj)
      implicit none
      type(rootedLLGrid),intent(in) :: obj
      real,pointer,dimension(:) :: ptr_lat_

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_lat_'

  ptr_lat_ => obj%lat

end function ptr_lat_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_lon_ - referencing longitudes
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_lon_(obj)
      implicit none
      type(rootedLLGrid),intent(in) :: obj
      real,pointer,dimension(:) :: ptr_lon_

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_lon_'

  ptr_lon_ => obj%lon

end function ptr_lon_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_lev_ - referencing levels
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_lev_(obj)
      implicit none
      type(rootedLLGrid),intent(in) :: obj
      real,pointer,dimension(:) :: ptr_lev_

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_lev_'

  ptr_lev_ => obj%lev

end function ptr_lev_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_kt_ - referencing var-types
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_kt_(obj)
      implicit none
      type(rootedLLGrid),intent(in) :: obj
      integer,pointer,dimension(:) :: ptr_kt_

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_kt_'

  ptr_kt_ => obj%kt

end function ptr_kt_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dsize_ - distributed vector size
!
! !DESCRIPTION:
!
! !INTERFACE:

    function dsize_(obj)
      use m_LVGrid   ,only : LVGrid_size
      use m_Collector,only : get
      implicit none
      type(rootedLLGrid),intent(in) :: obj
      integer :: dsize_

! !REVISION HISTORY:
! 	04Dec00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::dsize_'
  integer :: nIJ,nLV

	nLV=LVGrid_size(obj%lvG)	! v-dim size
	call get(obj%coll,count=nIJ)	! local h-dim size

  dsize_=nIJ*nLV
end function dsize_
end module m_rootedLLGrid

