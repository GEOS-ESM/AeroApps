!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Attributes - Data attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Attributes
      implicit none
      private	! except

      public :: Attributes		! The class data structure

      public :: Attributes_init,init	! create a skeleton object
      public :: wrap			! wrap arrays to an object or
      public :: subset			! create a subset

      public :: permute			! permute and unpermute
      public :: unpermute

      public :: clean			! clean an object

      public :: lsize		! attribute vector sizes
      public :: key		! attribute subset key
      public :: order		! the order of two attributes entires
      public :: aSuperset	! arg1 contains arg2

      public :: KR,KS,KX,KT,LAT,LON,LEV
      public :: KR_SUBSET
      public :: KS_SUBSET
      public :: LL_SUBSET
      public :: KX_SUBSET
      public :: OBS_SUBSET
      public :: FLD_SUBSET

      public :: LESS
      public :: GREATER

      public :: new
      public :: delete

      public :: ptr_lat		! object component of latitudes
      public :: ptr_lon		! object component of longitudes
      public :: ptr_lev		! object component of levels
      public :: ptr_kr		! object component of regions
      public :: ptr_kt		! object component of types
      public :: ptr_ks		! object component of sounding indices
      public :: ptr_kx		! object component of data sources

    type Attributes
      private

      integer :: n			! vector size
      integer :: key			! code for subset

      real   ,pointer,dimension(:) :: lon	! longitudes
      real   ,pointer,dimension(:) :: lat	! latitidues
      real   ,pointer,dimension(:) :: lev	! levels
      integer,pointer,dimension(:) :: kt	! variable types
      integer,pointer,dimension(:) :: kr	! regions
      integer,pointer,dimension(:) :: ks	! sounding indices
      integer,pointer,dimension(:) :: kx	! data sources

    end type Attributes

    interface Attributes_init; module procedure init_; end interface
    interface init  ; module procedure init_  ; end interface
    interface subset; module procedure subset_; end interface
    interface permute; module procedure permute_; end interface
    interface unpermute; module procedure unpermute_; end interface

    interface wrap; module procedure	&
      wrap_ob_,	&
      wrap_
    end interface

    interface clean; module procedure clean_; end interface

    interface key  ; module procedure key_  ; end interface
    interface lsize; module procedure lsize_; end interface
    interface order; module procedure order_; end interface
    interface aSuperset; module procedure	&
	aSuper_,	&
	aSuperAttr_; end interface

    interface ptr_lat; module procedure ptr_lat_; end interface
    interface ptr_lon; module procedure ptr_lon_; end interface
    interface ptr_lev; module procedure ptr_lev_; end interface
    interface ptr_kt ; module procedure ptr_kt_ ; end interface
    interface ptr_kr ; module procedure ptr_kr_ ; end interface
    interface ptr_ks ; module procedure ptr_ks_ ; end interface
    interface ptr_kx ; module procedure ptr_kx_ ; end interface

    interface new   ; module procedure new_   ; end interface
    interface delete; module procedure delete_; end interface

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Attributes'

  integer,parameter :: KR = 1	! 00000001
  integer,parameter :: KX = 2	! 00000010
  integer,parameter :: KS = 4	! 00000100
  integer,parameter :: KT = 8	! 00001000
  integer,parameter :: LAT=16	! 00010000
  integer,parameter :: LON=32	! 00100000
  integer,parameter :: LEV=64	! 01000000

  integer,parameter :: KR_SUBSET = KR		! code of kr only
  integer,parameter :: KS_SUBSET = KR+KS	! code for kr-ks
  integer,parameter :: KX_SUBSET = KR+KX+KS	! code for kr-kx-ks
  integer,parameter :: LL_SUBSET = KR+LON+LAT	! code for kr-lon-lat
  integer,parameter :: OBS_SUBSET= LEV+LON+LAT+KT+KS+KX+KR
  integer,parameter :: FLD_SUBSET= LEV+LON+LAT+KT+KR
  integer,parameter :: OBS_ALIAS = -OBS_SUBSET
  integer,parameter :: FLD_ALIAS = -FLD_SUBSET

  integer,parameter :: LESS    = -1
  integer,parameter :: GREATER = +1

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: wrap_ob_ - initialize an object of observation attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine wrap_ob_(obj, ndat,rlat,rlon,rlev,kx,ks,kt,gpart,comm)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mci
      use m_GlobalPartition,only : GlobalPartition
      use m_GlobalPartition,only : setkr
      implicit none

      type(Attributes),intent(out) :: obj

      integer,intent(in)  :: ndat
      real   ,target,dimension(:),intent(in) :: rlat
      real   ,target,dimension(:),intent(in) :: rlon
      real   ,target,dimension(:),intent(in) :: rlev
      integer,target,dimension(:),intent(in) :: kx
      integer,target,dimension(:),intent(in) :: ks
      integer,target,dimension(:),intent(in) :: kt
      type(GlobalPartition),intent(in) :: gpart
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::wrap_ob_'
  integer :: ier

  obj%n = ndat
  obj%key = OBS_ALIAS

  obj%lat => rlat(1:ndat)
  obj%lon => rlon(1:ndat)
  obj%lev => rlev(1:ndat)
  obj%kt  =>   kt(1:ndat)

  allocate(obj%kr(ndat),stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)
	if(mall_ison()) call mall_mci(obj%kr,myname)

  call setkr(obj%kr,gpart,kx,ks,rlat,rlon,comm)

  obj%kx => kx(1:ndat)
  obj%ks => ks(1:ndat)

end subroutine wrap_ob_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: wrap_ - initialize a object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine wrap_(obj, ndat,rlat,rlon,rlev,kt,gpart,comm)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mci
      use m_GlobalPartition,only : GlobalPartition
      use m_GlobalPartition,only : setkr
      implicit none

      type(Attributes),intent(out) :: obj

      integer,intent(in)  :: ndat
      real   ,target,dimension(:),intent(in) :: rlat
      real   ,target,dimension(:),intent(in) :: rlon
      real   ,target,dimension(:),intent(in) :: rlev
      integer,target,dimension(:),intent(in) :: kt
      type(GlobalPartition),intent(in) :: gpart
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::wrap_'
  integer :: ier

  obj%n = ndat
  obj%key = FLD_ALIAS

  obj%lat => rlat(1:ndat)
  obj%lon => rlon(1:ndat)
  obj%lev => rlev(1:ndat)
  obj%kt  =>   kt(1:ndat)

  allocate(obj%kr(ndat),stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)
	if(mall_ison()) call mall_mci(obj%kr,myname)

  call setkr(obj%kr,gpart,rlat,rlon,comm)

  nullify(obj%kx)
  nullify(obj%ks)

end subroutine wrap_

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
      use m_die, only : die
      implicit none
      type(Attributes),intent(inout) :: obj

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  ier=0
  select case(obj%key)		! can not use key_()
  case(FLD_ALIAS,OBS_ALIAS)
    call clean_alias_(obj)
  case(FLD_SUBSET,OBS_SUBSET,KR_SUBSET,LL_SUBSET,KS_SUBSET,KX_SUBSET)
    call clean_subset_(obj)
  case default
    ier=-1
  end select

  if(ier/=0) call die(myname_,'unknown subset',obj%key)
end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_alias_ - clean an object as formed by wrap_()
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_alias_(obj)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mco
      implicit none
      type(Attributes),intent(inout) :: obj

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_alias_'
  integer :: ier
  integer :: key

  key=-obj%key		! Note this procedure should only be called
			! from clean_() with %key checked for FLD_ALIAS
			! or OBS_ALIAS.

	if(mall_ison()) call mall_mco(obj%kr,myname)
  deallocate(obj%kr,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

  nullify(obj%lon)
  nullify(obj%lat)
  nullify(obj%lev)
  nullify(obj%kt)
  nullify(obj%kr)
  nullify(obj%ks)
  nullify(obj%kx)

  obj%n=-1
  obj%key=0
end subroutine clean_alias_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_subset_ - clean an object as a SUBSET
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_subset_(obj)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mco
      implicit none
      type(Attributes),intent(inout) :: obj

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_subset_'
  integer :: ier
  integer :: key

  obj%n=-1
		! Note this routine is call by clean_() only with
		! the obj%key value already verified by clean_().
  key=obj%key

!	if(mall_ison()) then
!	  if(iand(key,KR )/=0) call mall_mco(obj%kr ,myname)
!	  if(iand(key,KS )/=0) call mall_mco(obj%ks ,myname)
!	  if(iand(key,KX )/=0) call mall_mco(obj%kx ,myname)
!	  if(iand(key,KT )/=0) call mall_mco(obj%kt ,myname)
!	  if(iand(key,LAT)/=0) call mall_mco(obj%lat,myname)
!	  if(iand(key,LON)/=0) call mall_mco(obj%lon,myname)
!	  if(iand(key,LEV)/=0) call mall_mco(obj%lev,myname)
!	endif

  ier=-1
  select case(key)
  case(KR_SUBSET)
		if(mall_ison()) call mall_mco(obj%kr ,myname)
	deallocate(obj%kr,stat=ier)
  case(KS_SUBSET)
		if(mall_ison()) then
		  call mall_mco(obj%kr ,myname)
		  call mall_mco(obj%ks ,myname)
		endif
	deallocate(obj%kr,obj%ks,stat=ier)
  case(LL_SUBSET)
		if(mall_ison()) then
		  call mall_mco(obj%kr ,myname)
		  call mall_mco(obj%lon,myname)
		  call mall_mco(obj%lat,myname)
		endif
	deallocate(obj%kr,obj%lon,obj%lat,stat=ier)
  case(KX_SUBSET)
		if(mall_ison()) then
		  call mall_mco(obj%kr ,myname)
		  call mall_mco(obj%ks ,myname)
		  call mall_mco(obj%kx ,myname)
		endif
	deallocate(obj%kr,obj%ks,obj%kx,stat=ier)
  case(OBS_SUBSET)
		if(mall_ison()) then
		  call mall_mco(obj%kr ,myname)
		  call mall_mco(obj%ks ,myname)
		  call mall_mco(obj%kx ,myname)
		  call mall_mco(obj%kt ,myname)
		  call mall_mco(obj%lat,myname)
		  call mall_mco(obj%lon,myname)
		  call mall_mco(obj%lev,myname)
		endif
	deallocate( obj%kr ,obj%ks ,obj%kx ,obj%kt,	&
		    obj%lat,obj%lon,obj%lev,stat=ier)
  case(FLD_SUBSET)
		if(mall_ison()) then
		  call mall_mco(obj%kr ,myname)
		  call mall_mco(obj%kt ,myname)
		  call mall_mco(obj%lat,myname)
		  call mall_mco(obj%lon,myname)
		  call mall_mco(obj%lev,myname)
		endif
	deallocate( obj%kr ,                obj%kt,	&
		    obj%lat,obj%lon,obj%lev,stat=ier)
  end select

	if(ier/=0) call die(myname_,'deallocate()',ier)

  obj%key=0

  nullify(obj%lon)
  nullify(obj%lat)
  nullify(obj%lev)
  nullify(obj%kt )
  nullify(obj%kr )
  nullify(obj%ks )
  nullify(obj%kx )

end subroutine clean_subset_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a skeleton object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(obj,key,lsize)
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      type(Attributes),intent(out) :: obj
      integer,intent(in) :: key
      integer,intent(in) :: lsize

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier

  ier=0

  nullify(obj%lon)
  nullify(obj%lat)
  nullify(obj%lev)
  nullify(obj%kt)
  nullify(obj%kr)
  nullify(obj%ks)
  nullify(obj%kx)

  obj%n=0
  obj%key=0

  select case(abs(key))
  case(KR_SUBSET)

	allocate(obj%kr(lsize),stat=ier)
		if(ier/=0) call die(myname_,'allocate(KR_SUBSET)',ier)
		if(mall_ison()) call mall_mci(obj%kr,myname)

    obj%n=lsize
    obj%key=KR_SUBSET

  case(KS_SUBSET)

	allocate(obj%kr(lsize),obj%ks(lsize),stat=ier)
		if(ier/=0) call die(myname_,'allocate(KS_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(obj%kr,myname)
		  call mall_mci(obj%ks,myname)
		endif

    obj%n=lsize
    obj%key=KS_SUBSET

  case(LL_SUBSET)

	allocate(obj%kr(lsize),obj%lon(lsize),obj%lat(lsize),stat=ier)
		if(ier/=0) call die(myname_,'allocate(LL_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(obj%kr,myname)
		  call mall_mci(obj%lon,myname)
		  call mall_mci(obj%lat,myname)
		endif

    obj%n=lsize
    obj%key=LL_SUBSET

  case(KX_SUBSET)

	allocate( obj%kr(lsize),obj%ks(lsize),obj%kx(lsize),stat=ier)
		if(ier/=0) call die(myname_,'allocate(KX_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(obj%kr,myname)
		  call mall_mci(obj%ks,myname)
		  call mall_mci(obj%kx,myname)
		endif

    obj%n=lsize
    obj%key=KX_SUBSET

  case(OBS_SUBSET)

	allocate( obj%lat(lsize),obj%lon(lsize),obj%lev(lsize),	&
		  obj%kr (lsize),obj%ks (lsize),obj%kx (lsize),	&
		  obj%kt (lsize),	stat=ier)
		if(ier/=0) call die(myname_,'allocate(OBS_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(obj%kr ,myname)
		  call mall_mci(obj%ks ,myname)
		  call mall_mci(obj%kx ,myname)
		  call mall_mci(obj%kt ,myname)
		  call mall_mci(obj%lat,myname)
		  call mall_mci(obj%lon,myname)
		  call mall_mci(obj%lev,myname)
		endif

    obj%n=lsize
    obj%key=OBS_SUBSET

  case(FLD_SUBSET)

	allocate( obj%lat(lsize),obj%lon(lsize),obj%lev(lsize),	&
		  obj%kr (lsize),obj%kt (lsize), stat=ier)
		if(ier/=0) call die(myname_,'allocate(FLD_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(obj%kr ,myname)
		  call mall_mci(obj%kt ,myname)
		  call mall_mci(obj%lat,myname)
		  call mall_mci(obj%lon,myname)
		  call mall_mci(obj%lev,myname)
		endif

    obj%n=lsize
    obj%key=FLD_SUBSET

  case default
    ier=-1
  end select

	if(ier/=0) call die(myname_,'unknown subset',key)

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: subset_ - create a subset of attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine subset_(subset,key,indices,data)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      type(Attributes),intent(out) :: subset
      integer,intent(in) :: key
      integer,dimension(:),intent(in) :: indices
      type(Attributes),intent(in) :: data

! !REVISION HISTORY:
! 	23Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::subset_'
  integer :: ier
  integer :: lsize_
  integer :: key_data

	! If input Attributes(data ) does not have all components
	! specified for subset by key, the subset will not be
	! generated.

  key_data=key_(data)	! to make sure the sign is not seen.

  if(.not.aSuper_(key_data,key))	&
	call die(myname_,'key(data)',key_data,'key',key)

  lsize_=size(indices)

  ier=0

  nullify(subset%lon)
  nullify(subset%lat)
  nullify(subset%lev)
  nullify(subset%kt)
  nullify(subset%kr)
  nullify(subset%ks)
  nullify(subset%kx)

  select case(key)
  case(KR_SUBSET)

	allocate(subset%kr(lsize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate(KR_SUBSET)',ier)
		if(mall_ison()) call mall_mci(subset%kr,myname)

    subset%kr(:)=data%kr(indices(:))

  case(KS_SUBSET)

	allocate(subset%kr(lsize_),subset%ks(lsize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate(KS_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(subset%kr,myname)
		  call mall_mci(subset%ks,myname)
		endif

    subset%kr(:)=data%kr(indices(:))
    subset%ks(:)=data%ks(indices(:))

  case(LL_SUBSET)

	allocate(subset%kr (lsize_),subset%lon(lsize_),	&
		 subset%lat(lsize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate(LL_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(subset%kr,myname)
		  call mall_mci(subset%lon,myname)
		  call mall_mci(subset%lat,myname)
		endif

    subset%kr(:)=data%kr(indices(:))
    subset%lon(:)=data%lon(indices(:))
    subset%lat(:)=data%lat(indices(:))

  case(KX_SUBSET)

	allocate( subset%kr(lsize_),subset%ks(lsize_),	&
		  subset%kx(lsize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate(KX_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(subset%kr,myname)
		  call mall_mci(subset%ks,myname)
		  call mall_mci(subset%kx,myname)
		endif

    subset%kr(:)=data%kr(indices(:))
    subset%ks(:)=data%ks(indices(:))
    subset%kx(:)=data%kx(indices(:))

  case(OBS_SUBSET)

	allocate( subset%kr (lsize_),subset%ks (lsize_),	&
		  subset%kx (lsize_),subset%kt (lsize_),	&
		  subset%lat(lsize_),subset%lon(lsize_),	&
		  subset%lev(lsize_), stat=ier)
		if(ier/=0) call die(myname_,'allocate(OBS_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(subset%kr ,myname)
		  call mall_mci(subset%ks ,myname)
		  call mall_mci(subset%kx ,myname)
		  call mall_mci(subset%kt ,myname)
		  call mall_mci(subset%lat,myname)
		  call mall_mci(subset%lon,myname)
		  call mall_mci(subset%lev,myname)
		endif

    subset%kr (:)=data%kr (indices(:))
    subset%ks (:)=data%ks (indices(:))
    subset%kx (:)=data%kx (indices(:))
    subset%kt (:)=data%kt (indices(:))
    subset%lat(:)=data%lat(indices(:))
    subset%lon(:)=data%lon(indices(:))
    subset%lev(:)=data%lev(indices(:))

  case(FLD_SUBSET)

	allocate( subset%kr (lsize_),subset%kt (lsize_),	&
		  subset%lat(lsize_),subset%lon(lsize_),	&
		  subset%lev(lsize_), stat=ier)
		if(ier/=0) call die(myname_,'allocate(FLD_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(subset%kr ,myname)
		  call mall_mci(subset%kt ,myname)
		  call mall_mci(subset%lat,myname)
		  call mall_mci(subset%lon,myname)
		  call mall_mci(subset%lev,myname)
		endif

    subset%kr (:)=data%kr (indices(:))
    subset%kt (:)=data%kt (indices(:))
    subset%lat(:)=data%lat(indices(:))
    subset%lon(:)=data%lon(indices(:))
    subset%lev(:)=data%lev(indices(:))

  case default
    ier=-1
  end select

	if(ier/=0) call die(myname_,'unknown subset',key)

  subset%key=key
  subset%n=lsize_

end subroutine subset_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unsubset_ - create a subset of attributes in order
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unsubset_(indices,subset,key,data)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      integer,dimension(:),intent(in) :: indices
      type(Attributes),intent(out) :: subset
      integer,intent(in) :: key
      type(Attributes),intent(in) :: data

! !REVISION HISTORY:
! 	23Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unsubset_'
  integer :: ier
  integer :: lsize_
  integer :: key_data

	! If input Attributes(data ) does not have all components
	! specified for subset by key, the subset will not be
	! generated.

  key_data=key_(data)	! to make sure the sign is not seen.

  if(.not.aSuper_(key_data,key))	&
	call die(myname_,'key(data)',key_data,'key',key)

  lsize_=size(indices)

  ier=0

  nullify(subset%lon)
  nullify(subset%lat)
  nullify(subset%lev)
  nullify(subset%kt)
  nullify(subset%kr)
  nullify(subset%ks)
  nullify(subset%kx)

  select case(key)
  case(KR_SUBSET)

	allocate(subset%kr(lsize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate(KR_SUBSET)',ier)
		if(mall_ison()) call mall_mci(subset%kr,myname)

    subset%kr(indices(:))=data%kr(:)

  case(KS_SUBSET)

	allocate(subset%kr(lsize_),subset%ks(lsize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate(KS_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(subset%kr,myname)
		  call mall_mci(subset%ks,myname)
		endif

    subset%kr(indices(:))=data%kr(:)
    subset%ks(indices(:))=data%ks(:)

  case(LL_SUBSET)

	allocate(subset%kr (lsize_),subset%lon(lsize_),	&
		 subset%lat(lsize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate(LL_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(subset%kr,myname)
		  call mall_mci(subset%lon,myname)
		  call mall_mci(subset%lat,myname)
		endif

    subset%kr(indices(:))=data%kr(:)
    subset%lon(indices(:))=data%lon(:)
    subset%lat(indices(:))=data%lat(:)

  case(KX_SUBSET)

	allocate( subset%kr(lsize_),subset%ks(lsize_),	&
		  subset%kx(lsize_),stat=ier)
		if(ier/=0) call die(myname_,'allocate(KX_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(subset%kr,myname)
		  call mall_mci(subset%ks,myname)
		  call mall_mci(subset%kx,myname)
		endif

    subset%kr(indices(:))=data%kr(:)
    subset%ks(indices(:))=data%ks(:)
    subset%kx(indices(:))=data%kx(:)

  case(OBS_SUBSET)

	allocate( subset%kr (lsize_),subset%ks (lsize_),	&
		  subset%kx (lsize_),subset%kt (lsize_),	&
		  subset%lat(lsize_),subset%lon(lsize_),	&
		  subset%lev(lsize_), stat=ier)
		if(ier/=0) call die(myname_,'allocate(OBS_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(subset%kr ,myname)
		  call mall_mci(subset%ks ,myname)
		  call mall_mci(subset%kx ,myname)
		  call mall_mci(subset%kt ,myname)
		  call mall_mci(subset%lat,myname)
		  call mall_mci(subset%lon,myname)
		  call mall_mci(subset%lev,myname)
		endif

    subset%kr(indices(:))=data%kr(:)
    subset%ks(indices(:))=data%ks(:)
    subset%kx(indices(:))=data%kx(:)
    subset%kt(indices(:))=data%kt(:)
    subset%lon(indices(:))=data%lon(:)
    subset%lat(indices(:))=data%lat(:)
    subset%lev(indices(:))=data%lev(:)

  case(FLD_SUBSET)

	allocate( subset%kr (lsize_),subset%kt (lsize_),	&
		  subset%lat(lsize_),subset%lon(lsize_),	&
		  subset%lev(lsize_), stat=ier)
		if(ier/=0) call die(myname_,'allocate(FLD_SUBSET)',ier)
		if(mall_ison()) then
		  call mall_mci(subset%kr ,myname)
		  call mall_mci(subset%kt ,myname)
		  call mall_mci(subset%lat,myname)
		  call mall_mci(subset%lon,myname)
		  call mall_mci(subset%lev,myname)
		endif

    subset%kr(indices(:))=data%kr(:)
    subset%kt(indices(:))=data%kt(:)
    subset%lon(indices(:))=data%lon(:)
    subset%lat(indices(:))=data%lat(:)
    subset%lev(indices(:))=data%lev(:)

  case default
    ier=-1
  end select

	if(ier/=0) call die(myname_,'unknown subset',key)

  subset%key=key
  subset%n=lsize_

end subroutine unsubset_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: permute_ - create a new subset ordered by given indices
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine permute_(data,sorted,indices)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      type(Attributes)    ,intent(in ) :: data		! input
      type(Attributes)    ,intent(out) :: sorted	! output subset
      integer,dimension(:),intent(in ) :: indices	! which input

! !REVISION HISTORY:
! 	23Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::permute_'
  integer :: key_data

  key_data=key_(data)	! to make sure the sign is not seen.
  call subset_(sorted,key_data,indices,data)

end subroutine permute_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: unpermute_ - create a new subset ordered by given indices
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine unpermute_(data,unsorted,indices)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      type(Attributes)    ,intent(in ) :: data		! input
      type(Attributes)    ,intent(out) :: unsorted	! output subset
      integer,dimension(:),intent(in ) :: indices	! which output

! !REVISION HISTORY:
! 	23Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::unpermute_'
  integer :: key_data

  key_data=key_(data)	! to make sure the sign is not seen.
  call unsubset_(indices,unsorted,key_data,data)

end subroutine unpermute_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - local size
!
! !DESCRIPTION:
!
! !INTERFACE:

    function lsize_(obj)
      implicit none
      type(Attributes),intent(in) :: obj
      integer :: lsize_

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_=obj%n

end function lsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: key_ - get the subset key of a given object
!
! !DESCRIPTION:
!
! !INTERFACE:

    function key_(obj)
      implicit none
      type(Attributes),intent(in) :: obj
      integer :: key_

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::key_'

  key_=abs(obj%key)	! note that the sign of the key is only known
			! internally.

end function key_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_lat_ - referencing component %lat
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_lat_(obj)
      use m_die,only : die
      implicit none
      type(Attributes),intent(in) :: obj
      real,pointer,dimension(:) :: ptr_lat_

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_lat_'
  integer :: key

  key=key_(obj)
  if(iand(key,LAT)==0) call die(myname_,'component %lat undefined')

  ptr_lat_ => obj%lat

end function ptr_lat_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_lon_ - referencing component %lon
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_lon_(obj)
      use m_die,only : die
      implicit none
      type(Attributes),intent(in) :: obj
      real,pointer,dimension(:) :: ptr_lon_

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_lon_'
  integer :: key

  key=key_(obj)
  if(iand(key,LON)==0) call die(myname_,'component %lon undefined')

  ptr_lon_ => obj%lon

end function ptr_lon_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_lev_ - referencing component %lev
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_lev_(obj)
      use m_die,only : die
      implicit none
      type(Attributes),intent(in) :: obj
      real,pointer,dimension(:) :: ptr_lev_

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_lev_'
  integer :: key

  key=key_(obj)
  if(iand(key,LEV)==0) call die(myname_,'component %lev undefined')

  ptr_lev_ => obj%lev

end function ptr_lev_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_kt_ - referencing component %kt
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_kt_(obj)
      use m_die,only : die
      implicit none
      type(Attributes),intent(in) :: obj
      integer,pointer,dimension(:) :: ptr_kt_

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_kt_'
  integer :: key

  key=key_(obj)
  if(iand(key,KT)==0) call die(myname_,'component %kt undefined')

  ptr_kt_ => obj%kt

end function ptr_kt_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_kr_ - referencing component %kr
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_kr_(obj)
      use m_die,only : die
      implicit none
      type(Attributes),intent(in) :: obj
      integer,pointer,dimension(:) :: ptr_kr_

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_kr_'
  integer :: key

  key=key_(obj)
  if(iand(key,KR)==0) call die(myname_,'component %kr undefined')

  ptr_kr_ => obj%kr

end function ptr_kr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_kx_ - referencing component %kx
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_kx_(obj)
      use m_die,only : die
      implicit none
      type(Attributes),intent(in) :: obj
      integer,pointer,dimension(:) :: ptr_kx_

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_kx_'
  integer :: key

  key=key_(obj)
  if(iand(key,KX) ==0) call die(myname_,'component %kx undefined')

  ptr_kx_ => obj%kx

end function ptr_kx_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_ks_ - referencing component %ks
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_ks_(obj)
      use m_die,only : die
      implicit none
      type(Attributes),intent(in) :: obj
      integer,pointer,dimension(:) :: ptr_ks_

! !REVISION HISTORY:
! 	20Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_ks_'
  integer :: key

  key=key_(obj)
  if(iand(key,KS) ==0) call die(myname_,'component %ks undefined')

  ptr_ks_ => obj%ks

end function ptr_ks_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: order_ - order two entries
!
! !DESCRIPTION:
!
! !INTERFACE:

    function order_(i,attri,j,attrj,key)
      use m_die,only : die
      implicit none
      integer,intent(in) :: i
      type(Attributes),intent(in) :: attri
      integer,intent(in) :: j
      type(Attributes),intent(in) :: attrj
      integer,intent(in) :: key
      integer :: order_

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::order_'
  integer :: keyi,keyj

  keyi=key_(attri)
  keyj=key_(attrj)

  order_=0
  select case(key)
  case(KR_SUBSET,KS_SUBSET,LL_SUBSET,KX_SUBSET)

    if(attri%kr(i)<attrj%kr(j)) order_=LESS
    if(attri%kr(i)>attrj%kr(j)) order_=GREATER
    if(order_/=0) return

    select case(key)
    case(KS_SUBSET)

      if(attri%ks(i)<attrj%ks(j)) order_=LESS
      if(attri%ks(i)>attrj%ks(j)) order_=GREATER

    case(LL_SUBSET)

      if(attri%lat(i)<attrj%lat(j)) order_=LESS
      if(attri%lat(i)>attrj%lat(j)) order_=GREATER
      if(order_/=0) return
      if(attri%lon(i)<attrj%lon(j)) order_=LESS
      if(attri%lon(i)>attrj%lon(j)) order_=GREATER

    case(KX_SUBSET)

      if(attri%kx(i)<attrj%kx(j)) order_=LESS
      if(attri%kx(i)>attrj%kx(j)) order_=GREATER
      if(order_/=0) return
      if(attri%ks(i)<attrj%ks(j)) order_=LESS
      if(attri%ks(i)>attrj%ks(j)) order_=GREATER

    end select

  case default
    call die(myname_,'unknown subset',key)
  end select

end function order_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: aSuper_ - key 1 contains (is a superset of) key 2
!
! !DESCRIPTION:
!
! !INTERFACE:

    function aSuper_(key1,key2)
      implicit none
      integer,intent(in) :: key1
      integer,intent(in) :: key2
      logical :: aSuper_

! !REVISION HISTORY:
! 	26Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::aSuper_'

  aSuper_=iand(key1,key2)==key2

end function aSuper_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: aSuperAttr_ - Attribute 1 contains (is a superset of) 2
!
! !DESCRIPTION:
!
! !INTERFACE:

    function aSuperAttr_(attr1,attr2)
      implicit none
      type(Attributes),intent(in) :: attr1
      type(Attributes),intent(in) :: attr2
      logical :: aSuperAttr_

! !REVISION HISTORY:
! 	26Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::aSuperAttr_'
  integer :: key1,key2

  key1=key_(attr1)
  key2=key_(attr2)
  aSuperAttr_=iand(key1,key2)==key2

end function aSuperAttr_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: new_ - create an object for a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    function new_(mold,stat)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_ci
      implicit none
      type(Attributes),pointer :: mold
      integer,optional,intent(out) :: stat
      type(Attributes),pointer :: new_

! !REVISION HISTORY:
! 	28Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::new_'
  type(Attributes),pointer :: obj
  integer :: ier

  if(present(stat)) stat=0
  allocate(obj,stat=ier)
	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_,'allocate()',ier)
	  stat=ier
	  return
	endif

	if(mall_ison()) call mall_ci(1,myname)

  new_ => obj
  nullify(obj)		! to prevent the compiler touching the memory.
end function new_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: delete_ - delete an object for a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine delete_(obj,stat)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(Attributes),pointer :: obj
      integer,optional,intent(out) :: stat


! !REVISION HISTORY:
! 	28Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::delete_'
  integer :: ier

  if(present(stat)) stat=0

	if(mall_ison()) call mall_co(1,myname)

  deallocate(obj,stat=ier)
	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_,'deallocate()',ier)
	  stat=ier
	  return
	endif

end subroutine delete_

end module m_Attributes
