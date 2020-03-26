!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_rootedAIGrid - an implementation of a distributed AIGrid
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_rootedAIGrid
      use m_GlobalPartition,only : GlobalPartition
      use m_Distribution,only : Distribution
      use m_Attributes  ,only : Attributes
      implicit none
      private	! except

      public :: rootedAIGrid		! The class data structure
      public :: rootedAIGrid_init,init	! initialize an object
      public :: clean			! clean a distributed grid

      public :: distrSize		! size of the distrbuted vector
      public :: gatherv			! create a gathered vector
      public :: cleangv			! clean the gathered vector
      public :: rootedAIGrid_intp,intp	! interpolate to a LL-grid

      public :: ptr_lat		! Component referencing, latitudes
      public :: ptr_lon		! Component referencing, longitudes
      public :: ptr_lev		! Component referencing, levels
      public :: ptr_kt		! Component referencing, var-types

      public :: ptr_attr	! referencing Attributes vector
      public :: ptr_dstr	! referencing a Distribution

      type rootedAIGrid
	private
	integer :: comm
	logical :: root
	type(GlobalPartition),pointer :: gpart
	type(Distribution),pointer :: dstr
	type(Attributes)  ,pointer :: attr
	real,pointer,dimension(:)  :: vinc0	! a rooted grid vector
      end type rootedAIGrid

      interface rootedAIGrid_init; module procedure	&
	init_,		&
	init2_; end interface

      interface init; module procedure	&
	init_,		&
	init2_; end interface

      interface clean  ; module procedure clean__  ; end interface
      interface distrSize; module procedure dsize_  ; end interface
      interface gatherv; module procedure gatherv_; end interface
      interface cleangv; module procedure cleanv_ ; end interface

      interface rootedAIGrid_intp; module procedure	&
	intps_,		&
	intpv_,		&
	intps2_,	&
	intpv2_; end interface

      interface intp; module procedure	&
	intps_,		&
	intpv_,		&
	intps2_,	&
	intpv2_; end interface

      interface ptr_lat; module procedure ptr_lat__; end interface
      interface ptr_lon; module procedure ptr_lon__; end interface
      interface ptr_lev; module procedure ptr_lev__; end interface
      interface ptr_kt ; module procedure ptr_kt__ ; end interface
      interface ptr_dstr; module procedure ptr_dstr_; end interface
      interface ptr_attr; module procedure ptr_attr_; end interface

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_rootedAIGrid'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a AIGrid with scattered attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(aig,im,jnp,mlev,plevs,KTList,root,comm,rsrc)
      use m_mpif90,only : MP_comm_rank
      use m_die ,only : MP_die,die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_AIGrid,only : AIGrid_init
      use m_AIGrid,only : AIGrid_bind
      use m_GlobalPartition,only : new,GlobalPartition_init
      use m_Distribution,only : distribution_init
      use m_Distribution,only : new
      use m_Distribution,only : distrSize
      use m_Attributes,only : Attributes
      use m_Attributes,only : wrap
      use m_Attributes,only : clean
      use m_Attributes,only : new
      use m_Attributes,only : KR_SUBSET
      use m_MultiAccessNavigator,only : MultiAccessNavigator
      use m_MultiAccessNavigator,only : clean
      use m_zeit ,only : zeit_ci,zeit_co
      use m_mpout,only : mpout
      use m_showDistrib,only : showDistrib
      implicit none

      type(rootedAIGrid),intent(out) :: aig	! output arguments

		! input arguments used only on root

      integer,intent(in) :: im
      integer,intent(in) :: jnp
      integer,intent(in) :: mlev
      real,dimension(:),intent(in) :: plevs
      integer,dimension(:),intent(in) :: KTList

      integer,intent(in) :: root
      integer,intent(in) :: comm

      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: myID
  integer :: ier
  integer :: ninc0
  real   ,allocatable,dimension(:) :: vlat0,vlon0,vlev0
  integer,allocatable,dimension(:) :: vkts0
  type(Attributes) :: ai_attr0
  type(MultiAccessNavigator) :: ai_man

  aig%comm=comm
  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  if(myID==root) then
    call AIGrid_init(ninc0,im,jnp,mlev,plevs,	&
	size(KTList),KTList,rsrc=rsrc )
    aig%root=.true.
  else
    ninc0=0
    aig%root=.false.
  endif

	call showDistrib(mpout,'pre-distrib-ninc',ninc0,	&
		root,comm,listall=.true.)

		! Allocate the memories

	allocate( vlat0(ninc0),vlon0(ninc0),	&
		  vlev0(ninc0),vkts0(ninc0), stat=ier	)
		if(ier/=0) call die(myname_,'allocate(vinc0)',ier)
		if(mall_ison()) then
		  call mall_mci(vlat0,myname)
		  call mall_mci(vlon0,myname)
		  call mall_mci(vlev0,myname)
		  call mall_mci(vkts0,myname)
		endif

  if(aig%root) call AIGrid_bind(ninc0, vlat0,vlon0,vlev0,vkts0 )

	aig%gpart => new(aig%gpart)
  call GlobalPartition_init(aig%gpart,comm=comm,root=root,rsrc=rsrc)

	! Attributes of the AIGrid vector are defined only on the
	! root PE.  On other PEs, the attributes are defined as
	! size zero.  ai_attr0 is a handle to those attributes.

	call wrap(ai_attr0, ninc0,vlat0,vlon0,vlev0,vkts0,	&
		aig%gpart,comm)

	! Create a Distribution for analysis increments
	! (aig%dstr), a local MultiAccssNavigator (ai_man), and
	! a copy of distributed Attributes (aig%attr).

  aig%dstr => new(aig%dstr)
  aig%attr => new(aig%attr)

	call zeit_ci('AIGrid_distr_3d')
  call distribution_init(aig%dstr,aig%attr,ai_man,	&
	(/KR_SUBSET/),ai_attr0,comm)
	call zeit_co('AIGrid_distr_3d')

  nullify(aig%vinc0)

		! Remove the handle, then the actual memories.

	call clean(ai_attr0)
		if(mall_ison()) then
		  call mall_mco(vlat0,myname)
		  call mall_mco(vlon0,myname)
		  call mall_mco(vlev0,myname)
		  call mall_mco(vkts0,myname)
		endif
	deallocate(vlat0,vlon0,vlev0,vkts0,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

	call clean(ai_man)

	call showDistrib(mpout,'post-distrib-ninc',	&
		distrSize(aig%dstr),root,comm,listall=.true.)

end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init2_ - initialize a 2-d AIGrid with scattered attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init2_(aig,im,jnp,KTList,root,comm,rsrc)
      use m_mpif90,only : MP_comm_rank
      use m_die ,only : MP_die,die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_AIGrid,only : AIGrid_init
      use m_AIGrid,only : AIGrid_bind
      use m_GlobalPartition,only : new,GlobalPartition_init
      use m_Distribution,only : distribution_init
      use m_Distribution,only : new
      use m_Distribution,only : distrSize
      use m_Attributes,only : Attributes
      use m_Attributes,only : wrap
      use m_Attributes,only : clean
      use m_Attributes,only : new
      use m_Attributes,only : KR_SUBSET
      use m_MultiAccessNavigator,only : MultiAccessNavigator
      use m_MultiAccessNavigator,only : clean
      use m_mpout,only : mpout
      use m_showDistrib,only : showDistrib
      implicit none

      type(rootedAIGrid),intent(out) :: aig	! output arguments

		! input arguments used only on root

      integer,intent(in) :: im
      integer,intent(in) :: jnp
      integer,dimension(:),intent(in) :: KTList

      integer,intent(in) :: root
      integer,intent(in) :: comm

      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init2_'
  integer :: myID
  integer :: ier
  integer :: ninc0
  real   ,allocatable,dimension(:) :: vlat0,vlon0,vlev0
  integer,allocatable,dimension(:) :: vkts0
  type(Attributes) :: ai_attr0
  type(MultiAccessNavigator) :: ai_man

  aig%comm=comm
  call MP_comm_rank(comm,myID,ier)
	if(ier/=0) call MP_die(myname_,'MP_comm_rank()',ier)

  if(myID==root) then
    call AIGrid_init(ninc0,im,jnp,size(KTList),KTList,rsrc=rsrc)
    aig%root=.true.
  else
    ninc0=0
    aig%root=.false.
  endif

	call showDistrib(mpout,'pre-distrib-ninc',ninc0,	&
		root,comm,listall=.true.)

		! Allocate the memories

	allocate( vlat0(ninc0),vlon0(ninc0),	&
		  vlev0(ninc0),vkts0(ninc0), stat=ier	)
		if(ier/=0) call die(myname_,'allocate(vinc0)',ier)
		if(mall_ison()) then
		  call mall_mci(vlat0,myname)
		  call mall_mci(vlon0,myname)
		  call mall_mci(vlev0,myname)
		  call mall_mci(vkts0,myname)
		endif

  if(aig%root) call AIGrid_bind(ninc0, vlat0,vlon0,vlev0,vkts0 )

	aig%gpart => new(aig%gpart)
  call GlobalPartition_init(aig%gpart,comm=comm,root=root,rsrc=rsrc)

	! Attributes of the AIGrid vector are defined only on the
	! root PE.  On other PEs, the attributes are defined as
	! size zero.  ai_attr0 is a handle to those attributes.

	call wrap(ai_attr0, ninc0,vlat0,vlon0,vlev0,vkts0,	&
		aig%gpart,comm)

	! Create a Distribution for analysis increments
	! (aig%dstr), a local MultiAccssNavigator (ai_man), and
	! a copy of distributed Attributes (aig%attr).

  aig%dstr => new(aig%dstr)
  aig%attr => new(aig%attr)

  call distribution_init(aig%dstr,aig%attr,ai_man,	&
	(/KR_SUBSET/),ai_attr0,comm)

  nullify(aig%vinc0)

		! Remove the handle, then the actual memories.

	call clean(ai_attr0)
		if(mall_ison()) then
		  call mall_mco(vlat0,myname)
		  call mall_mco(vlon0,myname)
		  call mall_mco(vlev0,myname)
		  call mall_mco(vkts0,myname)
		endif
	deallocate(vlat0,vlon0,vlev0,vkts0,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

	call clean(ai_man)

	call showDistrib(mpout,'post-distrib-ninc',	&
		distrSize(aig%dstr),root,comm,listall=.true.)

end subroutine init2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gatherv_ - gather vinc from all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine gatherv_(vinc,aig)
      use m_Distribution    ,only : localSize
      use m_DistributionComm,only : undistribute

      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci

      implicit none

      real,dimension(:) ,intent(in   ) :: vinc	! distributed vinc
      type(rootedAIGrid),intent(inout) :: aig	! a rooted grid info.

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gatherv_'
  integer :: ninc
  integer :: ier

  if(associated(aig%vinc0))	&
	call die(myname_,'associated(aig%vinc0)')

  ninc=localSize(aig%dstr)	! the source is (was) either size 0 or
				! the full size

  allocate(aig%vinc0(ninc),stat=ier)

	if(ier/=0) call die(myname_,'allocate()',ier)
	if(mall_ison()) call mall_mci(aig%vinc0,myname)

  call undistribute(vinc,aig%vinc0,aig%dstr,aig%comm)

end subroutine gatherv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: cleanv_ - clean %vinc only
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine cleanv_(aig)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_mco
      implicit none
      type(rootedAIGrid),intent(inout) :: aig

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::cleanv_'
  integer :: ier

  if(.not.associated(aig%vinc0))	&
	call die(myname_,'.not.associated(aig%vinc0)')

	if(mall_ison()) call mall_mco(aig%vinc0,myname)

  deallocate(aig%vinc0,stat=ier)

	if(ier/=0) call die(myname_,'deallocate()',ier)

  nullify(aig%vinc0)
end subroutine cleanv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean__ - clean an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean__(aig)
      use m_GlobalPartition,only : clean,delete
      use m_Distribution,only : clean,delete
      use m_Attributes  ,only : clean,delete
      use m_mpif90,only : MP_comm_null
      use m_AIGrid,only : AIGrid_clean
      implicit none
      type(rootedAIGrid),intent(inout) :: aig

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean__'

  if(associated(aig%vinc0)) call cleanv_(aig)

  call clean(aig%dstr)
  call delete(aig%dstr)

  call clean(aig%attr)
  call delete(aig%attr)

  call clean(aig%gpart)
  call delete(aig%gpart)

  aig%comm=MP_comm_null
  if(aig%root) then
    call AIGrid_clean()
    aig%root=.false.
  endif


end subroutine clean__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intps_ - interpolate vinc to a scalar grid on the root PE
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intps_(aig,im,jnp,mlev,plevs,ktz,z)
      use m_die ,only : die
      use m_AIGrid,only : AIGrid_intp
      implicit none

      type(rootedAIGrid),intent(in) :: aig

		! these arguments are significant only on the root PE,
		! which was defined when the (rootedAIGrid)aig is
		! initialized.

      integer,intent(in) :: im
      integer,intent(in) :: jnp
      integer,intent(in) :: mlev
      real,dimension(:),intent(in) :: plevs
      integer,intent(in) :: ktz
      real,dimension(:,:,:),intent(out) :: z

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intps_'
  integer :: ninc0

  if(.not.associated(aig%vinc0))		&
	call die(myname_,'.not.associated(aig%vinc0)')

  if(aig%root) then
    ninc0=size(aig%vinc0)
    call AIGrid_intp(ktz,    im,jnp,mlev,plevs,z,  ninc0,aig%vinc0)
  endif

end subroutine intps_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intpv_ - interpolate vinc to a scalar grid on the root PE
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intpv_(aig,im,jnp,mlev,plevs,ktu,ktv,u,v)
      use m_die ,only : die
      use m_AIGrid,only : AIGrid_intp
      implicit none

      type(rootedAIGrid),intent(in) :: aig

		! these arguments are significant only on the root PE,
		! which was defined when the (rootedAIGrid)aig is
		! initialized.

      integer,intent(in) :: im
      integer,intent(in) :: jnp
      integer,intent(in) :: mlev
      real,dimension(:),intent(in) :: plevs

      integer,intent(in) :: ktu,ktv
      real,dimension(:,:,:),intent(out) :: u
      real,dimension(:,:,:),intent(out) :: v

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intpv_'
  integer :: ninc0

  if(.not.associated(aig%vinc0))		&
	call die(myname_,'.not.associated(aig%vinc0)')

  if(aig%root) then
    ninc0=size(aig%vinc0)
    call AIGrid_intp(ktu,ktv,im,jnp,mlev,plevs,u,v,ninc0,aig%vinc0)
  endif

end subroutine intpv_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intps2_ - interpolate vinc to a scalar grid on the root PE
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intps2_(aig,im,jnp,ktz,z)
      use m_die ,only : die
      use m_AIGrid,only : AIGrid_intp
      implicit none

      type(rootedAIGrid),intent(in) :: aig

		! these arguments are significant only on the root PE,
		! which was defined when the (rootedAIGrid)aig is
		! initialized.

      integer,intent(in) :: im
      integer,intent(in) :: jnp
      integer,intent(in) :: ktz
      real,dimension(:,:),intent(out) :: z

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intps2_'
  integer :: ninc0

  if(.not.associated(aig%vinc0))		&
	call die(myname_,'.not.associated(aig%vinc0)')

  if(aig%root) then
    ninc0=size(aig%vinc0)
    call AIGrid_intp(ktz,    im,jnp,z,  ninc0,aig%vinc0)
  endif

end subroutine intps2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: intpv2_ - interpolate vinc to a scalar grid on the root PE
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine intpv2_(aig,im,jnp,ktu,ktv,u,v)
      use m_die ,only : die
      use m_AIGrid,only : AIGrid_intp
      implicit none

      type(rootedAIGrid),intent(in) :: aig

		! these arguments are significant only on the root PE,
		! which was defined when the (rootedAIGrid)aig is
		! initialized.

      integer,intent(in) :: im
      integer,intent(in) :: jnp
      integer,intent(in) :: ktu,ktv
      real,dimension(:,:),intent(out) :: u
      real,dimension(:,:),intent(out) :: v

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::intpv2_'
  integer :: ninc0

  if(.not.associated(aig%vinc0))		&
	call die(myname_,'.not.associated(aig%vinc0)')

  if(aig%root) then
    ninc0=size(aig%vinc0)
    call AIGrid_intp(ktu,ktv,im,jnp,u,v,ninc0,aig%vinc0)
  endif

end subroutine intpv2_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_lat__ - referencing latitudes
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_lat__(aig)
      use m_Attributes,only : ptr_lat
      implicit none
      type(rootedAIGrid),intent(in) :: aig
      real,pointer,dimension(:) :: ptr_lat__

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_lat__'

  ptr_lat__ => ptr_lat(aig%attr)

end function ptr_lat__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_lon__ - referencing longitudes
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_lon__(aig)
      use m_Attributes,only : ptr_lon
      implicit none
      type(rootedAIGrid),intent(in) :: aig
      real,pointer,dimension(:) :: ptr_lon__

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_lon__'

  ptr_lon__ => ptr_lon(aig%attr)

end function ptr_lon__
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_lev__ - referencing levels
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_lev__(aig)
      use m_Attributes,only : ptr_lev
      implicit none
      type(rootedAIGrid),intent(in) :: aig
      real,pointer,dimension(:) :: ptr_lev__

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_lev__'

  ptr_lev__ => ptr_lev(aig%attr)

end function ptr_lev__
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_kt__ - referencing var-types
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_kt__(aig)
      use m_Attributes,only : ptr_kt
      implicit none
      type(rootedAIGrid),intent(in) :: aig
      integer,pointer,dimension(:) :: ptr_kt__

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_kt__'

  ptr_kt__ => ptr_kt(aig%attr)

end function ptr_kt__
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_dstr_ - referencing %dstr
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_dstr_(aig)
      use m_Distribution,only : Distribution
      implicit none
      type(rootedAIGrid),intent(in) :: aig
      type(Distribution),pointer    :: ptr_dstr_

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_dstr_'

  ptr_dstr_ => aig%dstr

end function ptr_dstr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_attr_ - referencing %attr
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_attr_(aig)
      use m_Attributes,only : Attributes
      implicit none
      type(rootedAIGrid),intent(in) :: aig
      type(Attributes)  ,pointer    :: ptr_attr_

! !REVISION HISTORY:
! 	07Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_attr_'

  ptr_attr_ => aig%attr

end function ptr_attr_
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
      use m_Distribution,only : distrSize
      implicit none
      type(rootedAIGrid),intent(in) :: obj
      integer :: dsize_

! !REVISION HISTORY:
! 	04Dec00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::dsize_'

  dsize_=distrSize(obj%dstr)

end function dsize_
end module m_rootedAIGrid

