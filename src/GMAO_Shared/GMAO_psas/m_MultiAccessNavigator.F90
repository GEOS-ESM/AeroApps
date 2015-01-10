!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_MultiAccessNavigator - access data blocks in multiple ways
!
! !DESCRIPTION:
!
!	A MultiAccessNavigator allows an application of accessing a
!   given data set in multiple ways, as type blocks (kr-kt), as profile
!   blocks, or as sounding blocks.  During the initialization of a
!   MultiAccessNavigator, the module procedure will also return an
!   index array to let an application to sort the given data set to
!   make it accessible by the returned MultiAccessNavigator.
!
! !INTERFACE:

    module m_MultiAccessNavigator
      use m_Navigator,only : Navigator
      implicit none
      private	! except

      public :: MultiAccessNavigator	! The class data structure

      public :: MultiAccessNavigator_init
      public :: init
      public :: clean

      public :: new
      public :: delete

      public :: build_fromProf
      public :: build_fromType
      public :: build_fromSorted
      public :: build_fromCounts

      public :: getType,nType
      public :: getSndx,nSndx
      public :: getProf,nProf

		! a type block accessing navigator

      public :: ptr_typeNav

		! attributes of type blocks

	public :: ptr_krType
	public :: ptr_ktType
	public :: ptr_kxType

		! a sndx+prof block accessing naviagtor

      public :: ptr_sndxNav
      public :: ptr_profNav

		! attributes of sndx+prof blocks

	public :: ptr_krProf
	public :: ptr_kxProf
	public :: ptr_ksProf
	public :: ptr_latProf
	public :: ptr_lonProf
	public :: ptr_ktProf

    type MultiAccessNavigator
      private

			! Orderly access to blocks with data of
			! the same type

      type(Navigator),pointer :: typeNav

	integer,pointer,dimension(:) :: krType
	integer,pointer,dimension(:) :: ktType
	integer,pointer,dimension(:) :: kxType

			! Randomly access to blocks with data from
			! the same profile

      type(Navigator),pointer :: profNav

	integer,pointer,dimension(:) :: krProf
	integer,pointer,dimension(:) :: ktProf
	integer,pointer,dimension(:) :: kxProf
	integer,pointer,dimension(:) :: ksProf
	real   ,pointer,dimension(:) :: latProf
	real   ,pointer,dimension(:) :: lonProf

			! Orderly access to profNav with profiles
			! from the same sounding

      type(Navigator),pointer :: sndxNav

			! Information about the sounding may be
			! obtained from (krProf,ktProf,kxProf,ksProf)

    end type MultiAccessNavigator

    interface build_fromProf; module procedure		&
      build_fromProf_; end interface
    interface build_fromType; module procedure		&
      build_fromType_; end interface
    interface build_fromSorted; module procedure	&
      build_fromSorted_; end interface
    interface build_fromCounts; module procedure	&
      build_fromCounts_; end interface

    interface MultiAccessNavigator_init; module procedure	&
      init_; end interface

    interface init; module procedure	&
      init_; end interface

    interface clean; module procedure clean__; end interface

    interface    new; module procedure    new_; end interface
    interface delete; module procedure delete__; end interface

    interface getType; module procedure getType_; end interface
    interface getSndx; module procedure getSndx_; end interface
    interface getProf; module procedure getProf_; end interface

    interface nType; module procedure nType_; end interface
    interface nSndx; module procedure nSndx_; end interface
    interface nProf; module procedure nProf_; end interface

    interface ptr_typeNav; module procedure ptr_typeNav_; end interface
    interface ptr_krType; module procedure ptr_krType_; end interface
    interface ptr_ktType; module procedure ptr_ktType_; end interface
    interface ptr_kxType; module procedure ptr_kxType_; end interface

    interface ptr_sndxNav; module procedure ptr_sndxNav_; end interface

    interface ptr_profNav; module procedure ptr_profNav_; end interface
    interface ptr_krProf ; module procedure ptr_krProf_ ; end interface
    interface ptr_kxProf ; module procedure ptr_kxProf_ ; end interface
    interface ptr_ksProf ; module procedure ptr_ksProf_ ; end interface
    interface ptr_ktProf ; module procedure ptr_ktProf_ ; end interface
    interface ptr_latProf; module procedure ptr_latProf_; end interface
    interface ptr_lonProf; module procedure ptr_lonProf_; end interface


! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_MultiAccessNavigator'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - instantiate components of an object
!
! !DESCRIPTION:
!
!	Note that interfaces like this one is needed if any additional
!   module is expected to be developed to actually initialize an object
!   with components of pointer attribute.  This procedure instantiate
!   all dynamic components with fixed sizes.
!
! !INTERFACE:

    subroutine init_(man,nType,nProf,obs,kxblocking,stat)
      use m_die ,only : die,perr
      use m_mall,only : mall_ison,mall_mci
      use m_Navigator,only : new
      implicit none
      type(MultiAccessNavigator),intent(inout) :: man
      integer,optional,intent(in) :: nType	! initialize typeNav
      integer,optional,intent(in) :: nProf	! initialize profNav
      logical,optional,intent(in) :: obs	! one of two datatypes
      logical,optional,intent(in) :: kxblocking	! one of two datatypes
      integer,optional,intent(out) :: stat	! return status code

! !REVISION HISTORY:
! 	13Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier

! !REMARKS:
! Information about the sounding may be obtained from 
! (krProf,ktProf,kxProf,ksProf), therefore, we do not
! have nSndx as an input argument to init_.
! init_ is called in the build_from* subroutines.
! Check different ways of usage therein.

  if(present(stat)) stat=0
!________________________________________

  if( .not.present(nType) .and.	&
      .not.present(nProf) ) then

    if(present(obs).or.present(kxblocking)) then
      if(present(obs))	call perr(myname_,'unexpected obs')
      if(present(kxblocking))	&
			call perr(myname_,'unexpected kxblocking')
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

	! Instantiate all navigators

    man%typeNav => new(man%typeNav,stat=ier)
    if(ier==0) man%profNav => new(man%profNav,stat=ier)
    if(ier==0) man%sndxNav => new(man%sndxNav,stat=ier)
    if(ier/=0) then
      call perr(myname_,'new()',ier)
      if(.not.present(stat)) call die(myname_)
      stat=ier
      return
    endif
  endif
!________________________________________

  if(present(nType)) then

	! Instantiate type-block attributes

    if(present(obs)) then
      call perr(myname_,'unexpected obs for type-blocking')
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    if(.not.present(kxblocking)) then
      call perr(myname_,'unspecified kxblocking for type-blocking')
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    if(kxblocking) then
      allocate(man%krType(nType),man%ktType(nType),	&
	man%kxType(nType),stat=ier)
    else
      allocate(man%krType(nType),man%ktType(nType),stat=ier)
      nullify(man%kxType)	!related to 'if(associated(...))' in
				!clean__()
    endif

	if(ier/=0) then
	  call perr(myname_,'allocate(type)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
	if(mall_ison()) then
	  call mall_mci(man%krType,myname)
	  call mall_mci(man%ktType,myname)
	  if(kxblocking) call mall_mci(man%kxType,myname)
	endif
  endif
!________________________________________

  if(present(nProf)) then

	! Instantiate prof-block attributes

    if(present(kxblocking)) then
      call perr(myname_,'unexpected kxblocking for prof-blocking')
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    if(.not.present(obs)) then
      call perr(myname_,'unspecified obs for prof-blocking')
      if(.not.present(stat)) call die(myname_)
      stat=-1
      return
    endif

    if(obs) then
      nullify(man%latProf)
      nullify(man%lonProf)
      allocate(man%krProf(nProf),man%ktProf(nProf),	&
	man%ksProf(nProf),man%kxProf(nProf), stat=ier)
    else
      nullify(man%kxProf)
      nullify(man%ksProf)
      allocate(man%krProf(nProf),man%ktProf(nProf),	&
	man%latProf(nProf),man%lonProf(nProf), stat=ier)
    endif

	if(ier/=0) then
	  call perr(myname_,'allocate(prof)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
	if(mall_ison()) then
	  call mall_mci(man%krProf,myname)
	  call mall_mci(man%ktProf,myname)
	  if(obs) then
	    call mall_mci(man%kxProf,myname)
	    call mall_mci(man%ksProf,myname)
	  else
	    call mall_mci(man%latProf,myname)
	    call mall_mci(man%lonProf,myname)
	  endif
	endif
  endif

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean__ - clean an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean__(man,stat)
      use m_Navigator,only : clean
      use m_Navigator,only : delete
      use m_mall,only : mall_ison,mall_mco
      use m_die ,only : perr,die
      implicit none
      type(MultiAccessNavigator),intent(inout) :: man
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean__'
  integer :: ier

  if(present(stat)) stat=0	! to be extended

	! clean their contents

  call clean(man%typeNav)
  call clean(man%sndxNav)
  call clean(man%profNav)

	! actually deallocate them

  call delete(man%typeNav)
  call delete(man%sndxNav)
  call delete(man%profNav)

	if(mall_ison()) then

		! clean typeNav stuff

	  call mall_mco(man%krType,myname)
	  call mall_mco(man%ktType,myname)
	  if(associated(man%kxType))	&
	  	call mall_mco(man%kxType,myname)

		! clean profNav stuff

	  call mall_mco(man%krProf,myname)
	  if(associated(man%kxProf)) then
	    call mall_mco(man%kxProf,myname)
	    call mall_mco(man%ksProf,myname)
	  else
	    call mall_mco(man%latProf,myname)
	    call mall_mco(man%lonProf,myname)
	  endif
	  call mall_mco(man%ktProf,myname)
	endif

  if(associated(man%kxType)) then
    deallocate(man%krType,man%ktType,man%kxType,stat=ier)
  else
    deallocate(man%krType,man%ktType,stat=ier)
  endif
	if(ier/=0) then
	  call perr(myname_,'deallocate(type)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif

  if(associated(man%kxProf)) then
    deallocate(man%krProf,man%kxProf,man%ksProf,man%ktProf,stat=ier)
  else
    deallocate(man%krProf,man%ktProf,man%latProf,man%lonProf,stat=ier)
  endif
	if(ier/=0) then
	  call perr(myname_,'deallocate(prof)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=-1
	  return
	endif

end subroutine clean__

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: build_fromSorted_ - build an object by sorted kr-kt(-kx)-ks
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine build_fromSorted_(man,kr,kt,kx,ks,lat,lon,kxblocking)

      use m_Navigator,only : Navigator
      use m_Navigator,only : Navigator_init
      use m_Navigator,only : get
      use m_Navigator,only : lsize
      use m_Navigator,only : ptr_displs
      use m_Navigator,only : ptr_counts

      use m_die ,only : die
      use m_mall,only : mall_ison
      use m_mall,only : mall_mci,mall_mco

      use m_SortingTools,only : IndexSet
      use m_SortingTools,only : IndexSort

      implicit none

      type(MultiAccessNavigator),intent(out) :: man

		! data attributes.  They should have the same size.

      integer,optional,dimension(:),intent(in) :: kr
      integer,optional,dimension(:),intent(in) :: kt
      integer,optional,dimension(:),intent(in) :: kx
      integer,optional,dimension(:),intent(in) :: ks
      real   ,optional,dimension(:),intent(in) :: lat
      real   ,optional,dimension(:),intent(in) :: lon

      logical,optional,intent(in) :: kxblocking

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::build_fromSorted_'
  integer :: ier
  integer :: iType,nType
  integer :: iProf,nProf
  integer :: i,n,lc
  integer,allocatable,dimension(:) :: ixProf
  integer,pointer,dimension(:) :: displs
  integer,pointer,dimension(:) :: counts

  logical :: setobs
  logical :: setfld
  logical :: kxblox

	! Arguments kr, kt, etc., are not really optional.
        ! They must each present when calling this routine.

  setobs=present( kr).and.present( kt) .and.	&
	 present( kx).and.present( ks)

  setfld=present( kr).and.present( kt) .and.	&
	 present(lat).and.present(lon)

  if(.not.(setobs.or.setfld))	&	! at least one
		call die(myname_,'invalid arguments')
  if(setobs.and.setfld)		&	! at most one
		call die(myname_,'invalid arguments')

	! This test checks if the type-blocking should take into
	! account of separating kx blocks.

  kxblox=.false.
  if(setobs.and.present(kxblocking)) kxblox=kxblocking

	! The smallest array is used to define the logical size of the
	! input.

  n=min(size(kr),size(kt))
  if(setobs) then
    n=min(n,size(kx))
    n=min(n,size(ks))
  else
    n=min(n,size(lat))
    n=min(n,size(lon))
  endif
!________________________________________

		! instantiate navigators

	call init_(man,stat=ier)
		if(ier/=0) call die(myname_,'init_()',ier)
!________________________________________

	! Construct the type-Navigator, by two/three sorted keys
        ! in the order of kr-kt(-kx).

  if(kxblox) then
    call Navigator_init(man%typeNav,kr(1:n),kt(1:n),kx(1:n),stat=ier)
  else
    call Navigator_init(man%typeNav,kr(1:n),kt(1:n),        stat=ier)
  endif
	if(ier/=0) call die(myname_,'Navigator_init(typeNav)',ier)

  nType=lsize(man%typeNav)
		! instantiate attributes of type-blocks

	call init_(man,nType=nType,kxblocking=kxblox,stat=ier)
		if(ier/=0) call die(myname_,'init_(nType)',ier)

		! Collect type-block attributes

  do iType=1,nType
    call get(man%typeNav,iType,lc=lc)

    man%krType(iType)=kr(lc)
    man%ktType(iType)=kt(lc)
    if(kxblox) man%kxType(iType)=kx(lc)
  end do
!________________________________________

	! Construct the prof-Navigator, by three/four sorted keys
        ! in the order of kr-kt(-ks)-kx.

  if(setobs) then
    call Navigator_init(man%profNav,kr(1:n),kt(1:n),	&
				    kx(1:n),ks(1:n),stat=ier)
  else
    call Navigator_init(man%profNav, kr(1:n), kt(1:n),	&
	nint(1000.*lat(1:n)),nint(1000.*lon(1:n)),stat=ier)
  endif
	if(ier/=0) call die(myname_,'Navigator_init(profNav)',ier)

  nProf=lsize(man%profNav)
		! instantiate attributes of prof-blocks

	call init_(man,nProf=nProf,obs=setobs,stat=ier)
		if(ier/=0) call die(myname_,'init_(nProf)',ier)

		! Collect prof-block attributes

  do iProf=1,nProf
    call get(man%profNav,iProf,lc=lc)

    man%krProf(iProf)=kr(lc)
    man%ktProf(iProf)=kt(lc)
    if(setobs) then
      man%kxProf(iProf)= kx(lc)
      man%ksProf(iProf)= ks(lc)
    else
      man%latProf(iProf)=lat(lc)
      man%lonProf(iProf)=lon(lc)
    endif
  end do

		! Sort the profile navigator data by their sounding
		! indices, i.e. as kr(-kx)-ks-kt.  Notice that this sorting
		! order is different from the data order, which will
		! result in a reordered prof-blocks with a random access
		! profNav.

	allocate(ixProf(nProf),stat=ier)
		if(ier/=0) call die(myname_,'allocate(ix)',ier)
		if(mall_ison()) call mall_mci(ixProf,myname)

  call IndexSet(ixProf)
                ! Pay attention to the calling sequence:
  call IndexSort(ixProf,man%ktProf)
  if(setobs) then
    call IndexSort(ixProf,man%ksProf )
    call IndexSort(ixProf,man%kxProf )
  else
    call IndexSort(ixProf,man%lonProf)
    call IndexSort(ixProf,man%latProf)
  endif
  call IndexSort(ixProf,man%krProf)

		! These prof-block attributes are reordered.

  man%ktProf(:)=man%ktProf(ixProf(:))
  if(setobs) then
    man%ksProf(:)=man%ksProf(ixProf(:))
    man%kxProf(:)=man%kxProf(ixProf(:))
  else
    man%lonProf(:)=man%lonProf(ixProf(:))
    man%latProf(:)=man%latProf(ixProf(:))
  endif
  man%krProf(:)=man%krProf(ixProf(:))

		! The navigator data are also reordered, including
		! %displs(:).  This will result in a random access
		! navigator.

	displs => ptr_displs(man%profNav)
	counts => ptr_counts(man%profNav)

  displs(:)=displs(ixProf(:))
  counts(:)=counts(ixProf(:))

	nullify(displs)		! make sure no dangling pointers
	nullify(counts)

		if(mall_ison()) call mall_mco(ixProf,myname)
	deallocate(ixProf,stat=ier)
		if(ier/=0) call die(myname_,'deallocate(ix)',ier)
!________________________________________

	! Construct a sndx-block navigator of the sorted profNav, by
	! two/three sorted keys in the order of kr(-kx)-ks.

  if(setobs) then
    call Navigator_init(man%sndxNav,man%krProf,man%kxProf,	&
				    man%ksProf,stat=ier)
  else
    call Navigator_init(man%sndxNav,man%krProf,			&
	nint(1000.*man%latProf),nint(1000.*man%lonProf),stat=ier)
  endif
	if(ier/=0) call die(myname_,'Navigator_init(sndxNav)',ier)

end subroutine build_fromSorted_

subroutine showSndx_(sndxNav,profNav,mSndx)
  use m_Navigator,only : Navigator
  use m_Navigator,only : lsize
  use m_Navigator,only : get
  use m_stdio,only : stdout
  implicit none
  type(Navigator),intent(in) :: sndxNav
  type(Navigator),intent(in) :: profNav
  integer,intent(in) :: mSndx

  integer :: iSndx
  integer :: iProf
  integer :: nSndx
  integer :: lc,lcProf
  integer :: le,leProf

  !write out some sndx and prof values:

  write(stdout,*)' iSndx iProf   lc    le'
  write(stdout,*)' ----- -----   --    --'

  nSndx=lsize(sndxNav)
  do iSndx=1,min(mSndx,nSndx)
    call get(sndxNav,iSndx,lc=lcProf,le=leProf)
    do iProf=lcProf,leProf
      call get(profNav,iProf,lc=lc,le=le)
      write(stdout,23) iSndx,iProf,lc,le
23    format(1x,8(i5,1x))
    enddo
  enddo
end subroutine showSndx_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: build_fromProf_ - build from prof-blocks
!
! !DESCRIPTION:
!
!	It is assumed that the input data structure (profNav,krProf,
!   kxProf,ksProf,ktProf) is already sorted as prof-blocking (kr-kx-
!   ks-kt), and that no two blocks have the same key values (kr-
!   kx-ks-kt).  Although it may not be a problem to this routine, 
!   without the latter constraint, there could be a problem for 
!   some applications of this module.
!
!	The prof-block navigator of the MAN created by this subroutine
!   (man%profNav) has the same accessing order as the input profNav.
!   i.e. the attributes of the prof-blocks, krProf, ksProf, kxProf,
!   ktProf, including %counts(:) of man%profNav, are the same as that of
!   the input profNav.  However, %displs(:) of man%profNav are expected
!   to be different from of the input profNav, because the storage
!   order of the data points is expected to be reordered.
!
! !INTERFACE:

    subroutine build_fromProf_(man,indx,profNav,	&
	krProf,ksProf,ktProf,kxProf,latProf,lonProf,kxblocking)

      use m_SortingTools,only : indexSet,indexSort

      use m_Navigator,only : Navigator
      use m_Navigator,only : Navigator_init
      use m_Navigator,only : clean
      use m_Navigator,only : define
      use m_Navigator,only : weighted_define
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_Navigator,only : ptr_counts
      use m_Navigator,only : ptr_displs

      use m_die ,only : die
      use m_mall,only : mall_ison
      use m_mall,only : mall_mci,mall_mco

      implicit none

      type(MultiAccessNavigator),intent(out) :: man
      integer,dimension(:),intent(out) :: indx ! a permutation operator

      type(Navigator),intent(in) :: profNav
      integer,optional,dimension(:),intent(in) :: krProf ! region
      integer,optional,dimension(:),intent(in) :: ksProf ! sounding
      integer,optional,dimension(:),intent(in) :: ktProf ! variable
      integer,optional,dimension(:),intent(in) :: kxProf ! instrument
      real   ,optional,dimension(:),intent(in) :: latProf ! latitudes
      real   ,optional,dimension(:),intent(in) :: lonProf ! longitudes

      logical,optional,intent(in) :: kxblocking

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::build_fromProf_'

  integer :: lco,leo,lno
  integer :: lcx,lex,lnx
  integer :: i,ier
  type(Navigator) :: tmpNav

  integer,pointer,dimension(:) :: counts
  integer,pointer,dimension(:) :: displs

  integer,allocatable,dimension(:) :: ixProf

  integer :: nProf,iProf
  integer :: nType,iType
  integer :: msize
  logical :: setobs
  logical :: setfld
  logical :: kxblox

	! Arguments kr, kt, ks are not really optional.

  setobs=present( krProf).and.present( ktProf) .and.	&
	 present( kxProf).and.present( ksProf)

  setfld=present( krProf).and.present( ktProf) .and.	&
	 present(latProf).and.present(lonProf)

  if(.not.(setobs.or.setfld))	&	! at least one
		call die(myname_,'invalid arguments')
  if(setobs.and.setfld)		&	! at most one
		call die(myname_,'invalid arguments')

	! This test checks if the type-blocking should take into
	! account of separating kx blocks.

  kxblox=.false.
  if(setobs.and.present(kxblocking)) kxblox=kxblocking
!________________________________________

	! The size of (profNav) is a size well defined.

  nProf=lsize(profNav)
	if( nProf > size(krProf) ) call die(myname_,	&
		'nProf',nProf, 'size(krProf)',size(krProf))

	if(setobs) then
	  if( nProf > size(kxProf) ) call die(myname_,	&
		'nProf',nProf, 'size(kxProf)',size(kxProf))
	  if( nProf > size(ksProf) ) call die(myname_,	&
		'nProf',nProf, 'size(ksProf)',size(ksProf))
	else
	  if( nProf > size(latProf) ) call die(myname_,	&
		'nProf',nProf, 'size(latProf)',size(latProf))
	  if( nProf > size(lonProf) ) call die(myname_,	&
		'nProf',nProf, 'size(lonProf)',size(lonProf))
	endif

	if( nProf > size(ktProf) ) call die(myname_,	&
		'nProf',nProf, 'size(ktProf)',size(ktProf))

	allocate( ixProf(nProf), stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(ixProf,myname)
!________________________________________

	! Construct a permutation operator to sort block attribute data
	! from the input order (kr-kx-ks-kt) to the MAN order (kr-kt-kx-
	! ks), by prof-blocks.

  call indexSet(ixProf)
  if(setobs) then
    call indexSort(ixProf,ksProf)
    call indexSort(ixProf,kxProf)
  else
    call indexSort(ixProf,lonProf)
    call indexSort(ixProf,latProf)
  endif
  call indexSort(ixProf,ktProf)
  call indexSort(ixProf,krProf)
!________________________________________

	! Get the data structure ready.  This init_() instantiate
	! %profNav, %sndxNav, and %typeNav.

	call init_(man,stat=ier)
		if(ier/=0) call die(myname_,'init_()',ier)
!________________________________________

	! Construct the prof-block and the sndx-block navigators, with
	!
	!   nProf, krProf(), kxProf(), ksProf(), ktProf(),
	!   profNav%counts(ixProf)

  call Navigator_init(man%profNav,nProf)

		! Sort the navigator

	call Navigator_init(tmpNav,nProf)

		! Use the input %counts(:) to define the intermediate
		! prof-block navigator, tmpNav.  %displs(:) of tmpNav
		! is now defined according to the new order.

		counts => ptr_counts(profNav)
  call define(tmpNav,counts=counts(ixProf(:)))
		nullify(counts)

		! Use both %counts(:) of input profNav and %displs(:)
		! of the intermediate prof-block navigator tmpNav, to
		! define a new prof-block navigator with only its
		! %displs(:) are different from the input prof-block
		! navigator, profNav.

		displs => ptr_displs(man%profNav)
		counts => ptr_counts(man%profNav)

		! Get the counts(:) from the navigator of the input data
		! order. (kr-ks-kt)

  counts(:) = ptr_counts(profNav)	  ! Note this is not "=>"

		! Get the displs(:) from the navigator of the targetted
		! data order. (kr-kt-ks)

  displs(ixProf(:)) = ptr_displs(tmpNav)  ! Note this is not "=>"

		nullify(counts)
		nullify(displs)

	call clean(tmpNav)
	!________________________________________

	! Store the prof-block attributes

	call init_(man,nProf=nProf,obs=setobs,stat=ier)
		if(ier/=0) call die(myname_,'init_(nProf)',ier)

  man%krProf(:)=krProf(1:nProf)
  if(setobs) then
    man%kxProf(:)=kxProf(1:nProf)
    man%ksProf(:)=ksProf(1:nProf)
  else
    man%latProf(:)=latProf(1:nProf)
    man%lonProf(:)=lonProf(1:nProf)
  endif
  man%ktProf(:)=ktProf(1:nProf)
!________________________________________

	! Construct the sndx-block navigator from the prof-block
	! attributes

  if(setobs) then
    call Navigator_init(man%sndxNav,krProf,kxProf,ksProf)
  else
    call Navigator_init(man%sndxNav,krProf,	&
	nint(1000.*latProf),nint(1000.*lonProf))
  endif
!________________________________________

	! Construct the type-block navigator, with
	!
	!   krProf(ixProf),ktProf(ixProf),kxProf(ixProf),
	!   profNav%counts(ixProf)
	!
	! Note that nType is unknown until the new navigator is formed.

  if(kxblox) then
    call Navigator_init(man%typeNav,krProf(ixProf(:)),	&
		  ktProf(ixProf(:)),kxProf(ixProf(:)))
  else
    call Navigator_init(man%typeNav,krProf(ixProf(:)),	&
		  ktProf(ixProf(:)))
  endif

  nType=lsize(man%typeNav)

		!________________________________________
		!
		! Storing attributes kr-kt-kx needs to be done before
		! the unit of the counts/displs are redefined to data
		! points.  Right now they are profiles.

	call init_(man,nType=nType,kxblocking=kxblox,stat=ier)
		if(ier/=0) call die(myname_,'init_(nType)',ier)

  do iType=1,nType
		! for each type-block, get the first prof-block no. (l).

    call get(man%typeNav,iType,lc=iProf)

		! store the corresponding (pre-sorting) attributes.

    man%krType(iType)=krProf(ixProf(iProf))
    man%ktType(iType)=ktProf(ixProf(iProf))
    if(kxblox)		&
	man%kxType(iType)=kxProf(ixProf(iProf))
  end do

		! Redefine the type-block navigator with weights.  This
		! procedure converts the unit of the counts/displs from
		! the current per-profile to the required per-data
		! points.

	counts => ptr_counts(man%profNav)
  call weighted_define(man%typeNav,weights=counts(ixProf(:)))
	nullify(counts)

!________________________________________

	! Verify the size of indx(:)

  msize=0
  if(nType>0) call get(man%typeNav,nType,le=msize)

	if(msize>size(indx)) call die(myname_,	&
		'msize',msize,'size(indx)',size(indx))
!________________________________________

	! Construct the point-to-point permutation operator, with
	!
	!   profNav and man%profNav

  do iProf=1,nProf
		! for each prof-block, get the original segment.

    call get(    profNav,iProf,lc=lco,le=leo,ln=lno)

		! get the targetted locations of the same segment.

    call get(man%profNav,iProf,lc=lcx,le=lex,ln=lnx)

	if(lno/=lnx) call die(myname_,'lno',lno,'lnx',lnx)

	! lco:leo is where the data are as the input
	! lcx:lex is where the sorted data should be

    indx(lcx:lex)= (/(i,i=lco,leo)/)
  end do

		if(mall_ison()) call mall_mco(ixProf,myname)
	deallocate(ixProf,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine build_fromProf_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: build_fromCounts_ - build from counted profile-blocks
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine build_fromCounts_(man,indx,	&
	lnProf,krProf,latProf,lonProf,ktProf,kxProf,kxblocking)

      use m_Navigator,only : Navigator
      use m_Navigator,only : Navigator_init
      use m_Navigator,only : clean
      use m_Navigator,only : define

      implicit none

      type(MultiAccessNavigator),intent(out) :: man
      integer,dimension(:),intent(out) :: indx ! a permutation operator

      integer,optional,dimension(:),intent(in) :: lnProf ! counts
      integer,optional,dimension(:),intent(in) :: krProf ! region
      real,optional,dimension(:),intent(in) :: latProf ! sounding
      real,optional,dimension(:),intent(in) :: lonProf ! sounding
      integer,optional,dimension(:),intent(in) :: ktProf ! variable
      integer,optional,dimension(:),intent(in) :: kxProf ! instrument

      logical,optional,intent(in) :: kxblocking

! !REVISION HISTORY:
!       10Aug01 - L.P. Chang
!               - Systematic replacement of 'ks' with 'lat and lon' such as
!                 <    lnProf,krProf,ksProf,ktProf,kxProf,kxblocking)
!                 >    lnProf,krProf,latProf,lonProf,ktProf,kxProf,kxblocking)
!
!                 <    integer,optional,dimension(:),intent(in) :: ksProf ! sounding
!                 >    real,optional,dimension(:),intent(in) :: latProf ! sounding
!                 >    real,optional,dimension(:),intent(in) :: lonProf ! sounding
!
!                 <    ksProf=ksProf,ktProf=ktProf,kxProf=kxProf,      &
!                 >    latProf=latProf,lonProf=lonProf,ktProf=ktProf,kxProf=kxProf,&
! 	20Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::build_fromCounts_'

  type(Navigator) :: tmpNav
  integer :: nprof
  
                ! Construct a temporary prof-block navigator for
                ! column vectors.

  nprof=size(lnProf)

  call Navigator_init(tmpNav,nprof)
  call define(tmpNav,counts=lnProf)

  call build_fromProf_(man,indx,tmpNav,krProf=krProf,	&
	latProf=latProf,lonProf=lonProf,ktProf=ktProf,kxProf=kxProf,	&
	kxblocking=kxblocking)

  call clean(tmpNav)

end subroutine build_fromCounts_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: build_fromType_ - build from type-blocks
!
! !DESCRIPTION:
!
!	It is assumed that the input data structure (typeNav,krType,
!   ktType, kxType) is already sorted as type-blocking (kr-kt-kx), and
!   that no two blocks have the same key values (kr-kt-kx).
!   Without the latter constraint, there could be a problem for the
!   applications of this module.
!
! !INTERFACE:

    subroutine build_fromType_(man,indx,typeNav,	&
	krType,ktType,kxType,kx,ks,lat,lon)
 
      use m_Navigator,only : Navigator
      use m_Navigator,only : Navigator_init
      use m_Navigator,only : get
      use m_Navigator,only : lsize
      use m_Navigator,only : define
      use m_Navigator,only : clean
      use m_Navigator,only : ptr_counts
      use m_Navigator,only : ptr_displs
      use m_SortingTools,only : IndexSet
      use m_SortingTools,only : IndexSort
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none

      type(MultiAccessNavigator),intent(out) :: man

      integer,dimension(:),intent(out) :: indx	! has a packed size

      type(Navigator),intent(in) :: typeNav
      integer,optional,dimension(:),intent(in) :: krType
      integer,optional,dimension(:),intent(in) :: ktType
      integer,optional,dimension(:),intent(in) :: kxType

      integer,optional,dimension(:),intent(in) :: kx	! instruments
      integer,optional,dimension(:),intent(in) :: ks	! Sndx
      real   ,optional,dimension(:),intent(in) :: lat	! latitudes
      real   ,optional,dimension(:),intent(in) :: lon	! longitudes

! !REVISION HISTORY:
! 	13Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::build_fromType_'
  integer :: iType,nType
  integer :: iProf,nProf
  integer :: msize
  integer :: lc,le
  integer :: lco,leo,lno
  integer :: lcx,lex,lnx
  integer :: i,lnp
  integer :: displ
  integer :: ier
  type(Navigator) :: tmpNav
  integer,pointer,dimension(:) :: counts
  integer,pointer,dimension(:) :: displs
  integer,allocatable,dimension(:) :: ixProf
  integer,allocatable,dimension(:) :: icount
  integer,allocatable,dimension(:) :: itypes
  logical :: setobs
  logical :: setfld
  logical :: kxblox

	! Arguments krType, ktType, etc., are not really optional

  if( .not.present(krType) .or.	&
      .not.present(ktType) ) call die(myname_,'invalid arguments')

  if(present(kxType).and.present(kx))	&	! too many arguments
		call die(myname_,'presence of both kxType and kx')

  setobs=(present(kxType).or.present(kx)).and.present(ks)
  setfld=present(lat).and.present(lon)

  if(.not.(setobs.or.setfld))	&	! at least one
		call die(myname_,'not enough arguments')
  if(setobs.and.setfld)		&	! at most one
		call die(myname_,'too many arguments')

  kxblox=setobs.and.present(kxType)
!________________________________________

	call init_(man,stat=ier)
		if(ier/=0) call die(myname_,'init_()',ier)

	! Construct a typeNav with the same order, but packed and
	! ordered.

  nType=lsize(typeNav)
  call Navigator_init(man%typeNav,nType)

  call define(man%typeNav,counts=ptr_counts(typeNav))

	call init_(man,nType=nType,kxblocking=kxblox,stat=ier)
		if(ier/=0) call die(myname_,'init_(nType)',ier)

  man%krType(:)=krType(1:nType)
  man%ktType(:)=ktType(1:nType)
  if(kxblox) man%kxType(:)=kxType(1:nType)

!________________________________________

	! Construct the initial indx(:)

		! this is to find the memsize for the output

  msize=0
  if(nType>0) call get(man%typeNav,nType,le=msize)

	if(msize>size(indx)) call die(myname_,	&
		'msize',msize,'size(indx)',size(indx))

  do iType=1,nType
    call get(    typeNav,iType,lc=lco,le=leo)	! original locations
    call get(man%typeNav,iType,lc=lcx,le=lex)	! targetted locations

    indx(lcx:lex)=(/(i,i=lco,leo)/)
  end do
!________________________________________

	! Construct a profNav, from type-blocks with "local" sortings.

	allocate(icount(msize),itypes(msize),stat=ier)
		if(ier/=0) call die(myname_,'allocate(ii)',ier)
		if(mall_ison()) then
		  call mall_mci(icount,myname)
		  call mall_mci(itypes,myname)
		endif

  nProf=0
  do iType=1,nType

    call get(man%typeNav,iType,lc=lc,le=le)

	! Sectional index-sorting.  The full size ks(:) must be given,
	! since the indices were defined for the full array.

    if(setobs) then
      if(kxblox) then
        call IndexSort(indx(lc:le),ks)
        call Navigator_init(tmpNav,ks(indx(lc:le)))
      else
        call IndexSort(indx(lc:le),ks)
	call IndexSort(indx(lc:le),kx)
        call Navigator_init(tmpNav,kx(indx(lc:le)),ks(indx(lc:le)))
      endif
    else
      call IndexSort(indx(lc:le),lon)
      call IndexSort(indx(lc:le),lat)
      call Navigator_init(tmpNav,nint(1000.*lat(indx(lc:le))),	&
				 nint(1000.*lon(indx(lc:le)))	)
    endif

	lnp=lsize(tmpNav)	! number of prof-blocks for iType
	iProf=nProf+1
	nProf=nProf+lnp

    itypes(iProf:nProf)=iType
    icount(iProf:nProf)=ptr_counts(tmpNav)

	call clean(tmpNav)
  end do

  call Navigator_init(man%profNav,nProf)
  call define(man%profNav,counts=icount(1:nProf))

	call init_(man,nProf=nProf,obs=setobs,stat=ier)
		if(ier/=0) call die(myname_,'init_(nProf)',ier)

  do iProf=1,nProf

    iType=itypes(iProf)

    call get(man%profNav,iProf,lc=lc)
    man%krProf(iProf)=krType(iType)
    man%ktProf(iProf)=ktType(iType)
    if(setobs) then
      if(kxblox) then
	man%kxProf(iProf)=kxType(iType)
      else
	man%kxProf(iProf)=kx(indx(lc))
      endif
      man%ksProf(iProf)=ks(indx(lc))
    else
      man%latProf(iProf)=lat(indx(lc))
      man%lonProf(iProf)=lon(indx(lc))
    endif
  end do

		if(mall_ison()) then
		  call mall_mco(icount,myname)
		  call mall_mco(itypes,myname)
		endif
	deallocate(icount,itypes,stat=ier)
		if(ier/=0) call die(myname_,'deallocate(ii)',ier)

	! Index the order of prof-blocks as kr-(kx-)ks-kt

	allocate(ixProf(nProf),stat=ier)
		if(ier/=0) call die(myname_,'allocate(ix)',ier)
		if(mall_ison()) call mall_mci(ixProf,myname)

  call IndexSet(ixProf)
  call IndexSort(ixProf,man%ktProf)
  if(setobs) then
    call IndexSort(ixProf,man%ksProf)
    call IndexSort(ixProf,man%kxProf)
  else
    call IndexSort(ixProf,man%lonProf)
    call IndexSort(ixProf,man%latProf)
  endif
  call IndexSort(ixProf,man%krProf)

	! Sort the navigator and all attributes to the order of
	! kr-kx-ks-kt

	counts => ptr_counts(man%profNav)
	displs => ptr_displs(man%profNav)

  counts(:)=counts(ixProf(:))	! in place sorting
  displs(:)=displs(ixProf(:))	! in place sorting

	nullify(counts)
	nullify(displs)

  man%krProf(:)=man%krProf(ixProf(:))
  if(setobs) then
    man%kxProf(:)=man%kxProf(ixProf(:))
    man%ksProf(:)=man%ksProf(ixProf(:))
  else
    man%latProf(:)=man%latProf(ixProf(:))
    man%lonProf(:)=man%lonProf(ixProf(:))
  endif
  man%ktProf(:)=man%ktProf(ixProf(:))

		if(mall_ison()) call mall_mco(ixProf,myname)
	deallocate(ixProf,stat=ier)
		if(ier/=0) call die(myname_,'deallocate(ix)',ier)

	! Make a sndx-block navigator of %profNav

  if(setobs) then
    call Navigator_init(man%sndxNav,man%krProf(:),	&
		      man%kxProf(:),man%ksProf(:)	)
  else
    call Navigator_init(man%sndxNav,man%krProf(:),	&
		nint(1000.*man%latProf(:)),		&
		nint(1000.*man%lonProf(:))		)
  endif

end subroutine build_fromType_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getType_ - get information of a type-block
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getType_(man,i,displ,count,lc,ln,le,kr,kt,kx)
      use m_Navigator,only : get
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer,intent(in) :: i
      integer,optional,intent(out) :: displ
      integer,optional,intent(out) :: count
      integer,optional,intent(out) :: lc
      integer,optional,intent(out) :: ln
      integer,optional,intent(out) :: le
      integer,optional,intent(out) :: kr
      integer,optional,intent(out) :: kt
      integer,optional,intent(out) :: kx

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getType_'

	if(i<=0) call die(myname_,'invalid i',i)

  call get(man%typeNav,i,displ=displ,count=count,lc=lc,ln=ln,le=le)
  if(present(kr)) kr=man%krType(i)
  if(present(kt)) kt=man%ktType(i)
  if(present(kx)) then
    if(.not.associated(man%kxType))	&
	call die(myname_,'%kxType undefined')
    kx=man%kxType(i)
  endif

end subroutine getType_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getProf_ - get information of a prof-block
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getProf_(man,i,displ,count,lc,ln,le,kr,kt,kx,ks,lat,lon)
      use m_Navigator,only : get
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer,intent(in) :: i
      integer,optional,intent(out) :: displ
      integer,optional,intent(out) :: count
      integer,optional,intent(out) :: lc
      integer,optional,intent(out) :: ln
      integer,optional,intent(out) :: le
      integer,optional,intent(out) :: kr
      integer,optional,intent(out) :: kt
      integer,optional,intent(out) :: kx
      integer,optional,intent(out) :: ks
      real   ,optional,intent(out) :: lat
      real   ,optional,intent(out) :: lon

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getProf_'

	if(i<=0) call die(myname_,'invalid i',i)

  call get(man%profNav,i,displ=displ,count=count,lc=lc,ln=ln,le=le)
  if(present(kr)) kr=man%krProf(i)
  if(present(kt)) kt=man%ktProf(i)
  if(present(kx).or.present(ks)) then
    if(.not.associated(man%kxProf))	&
	call die(myname_,'%kxProf undefined')
    if(present(kx)) kx=man%kxProf(i)
    if(present(ks)) ks=man%ksProf(i)
  endif
  if(present(lat).or.present(lon)) then
    if(.not.associated(man%latProf))	&
	call die(myname_,'%latProf undefined')
    if(present(lat)) lat=man%latProf(i)
    if(present(lon)) lon=man%lonProf(i)
  endif

end subroutine getProf_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getSndx_ - get information of a sndx-block
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getSndx_(man,i,displ,count,lc,ln,le,kr,kx,ks,lat,lon)
      use m_Navigator,only : get
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer,intent(in) :: i
      integer,optional,intent(out) :: displ
      integer,optional,intent(out) :: count
      integer,optional,intent(out) :: lc
      integer,optional,intent(out) :: ln
      integer,optional,intent(out) :: le
      integer,optional,intent(out) :: kr
      integer,optional,intent(out) :: kx
      integer,optional,intent(out) :: ks
      real   ,optional,intent(out) :: lat
      real   ,optional,intent(out) :: lon

! !REVISION HISTORY:
!       10Aug01 - LPC
!               - Fix bug: two 'if(present(lat)) lat=man%latProf(lc_)'
!                 statements become one 'if(present(lat)) lat=man%latProf(lc_)'
!                 and one 'if(present(lon)) lon=man%lonProf(lc_)'.
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getSndx_'
  integer :: lc_

	if(i<=0) call die(myname_,'invalid i',i)

  call get(man%sndxNav,i,displ=displ,count=count,lc=lc_,ln=ln,le=le)
  if(present(lc)) lc=lc_

  if(present(kr)) kr=man%krProf(lc_)

  if(present(kx).or.present(ks)) then
    if(.not.associated(man%kxProf))	&
	call die(myname_,'%kxProf undefined')
    if(present(kx)) kx=man%kxProf(lc_)
    if(present(ks)) ks=man%ksProf(lc_)
  endif
  if(present(lat).or.present(lon)) then
    if(.not.associated(man%latProf))	&
	call die(myname_,'%latProf undefined')
    if(present(lat)) lat=man%latProf(lc_)
    if(present(lon)) lon=man%lonProf(lc_)
  endif

end subroutine getSndx_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nType_ - number of type-blocks
!
! !DESCRIPTION:
!
! !INTERFACE:

    function nType_(man)
      use m_Navigator,only : lsize
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer :: nType_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nType_'

  nType_=lsize(man%typeNav)

end function nType_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nProf_ - number of prof-blocks
!
! !DESCRIPTION:
!
! !INTERFACE:

    function nProf_(man)
      use m_Navigator,only : lsize
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer :: nProf_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nProf_'

  nProf_=lsize(man%profNav)

end function nProf_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nSndx_ - number of sndx-blocks
!
! !DESCRIPTION:
!
! !INTERFACE:

    function nSndx_(man)
      use m_Navigator,only : lsize
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer :: nSndx_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nSndx_'

  nSndx_=lsize(man%sndxNav)

end function nSndx_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_typeNav_ - referencing the type-block navigator
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_typeNav_(man)
      use m_Navigator,only : Navigator
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      type(Navigator),pointer :: ptr_typeNav_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_typeNav_'

  ptr_typeNav_ => man%typeNav

end function ptr_typeNav_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_krType_ - referencing a type-block attribute array
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_krType_(man)
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer,pointer,dimension(:) :: ptr_krType_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_krType_'

  if(.not.associated(man%krType))	&
	call die(myname_,'%krType undefined')

  ptr_krType_ => man%krType

end function ptr_krType_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_ktType_ - referencing a type-block attribute array
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_ktType_(man)
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer,pointer,dimension(:) :: ptr_ktType_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_ktType_'

  if(.not.associated(man%ktType))	&
	call die(myname_,'%ktType undefined')

  ptr_ktType_ => man%ktType

end function ptr_ktType_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_kxType_ - referencing a type-block attribute array
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_kxType_(man)
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer,pointer,dimension(:) :: ptr_kxType_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_kxType_'

  if(.not.associated(man%kxType))	&
	call die(myname_,'%kxType undefined')

  ptr_kxType_ => man%kxType

end function ptr_kxType_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_profNav_ - referencing the prof-block navigator
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_profNav_(man)
      use m_Navigator,only : Navigator
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      type(Navigator),pointer :: ptr_profNav_


! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_profNav_'

  ptr_profNav_ => man%profNav

end function ptr_profNav_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_krProf_ - referencing a prof-block attribute array
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_krProf_(man)
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer,pointer,dimension(:) :: ptr_krProf_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_krProf_'

  if(.not.associated(man%krProf))	&
	call die(myname_,'%krProf undefined')

  ptr_krProf_ => man%krProf

end function ptr_krProf_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_ktProf_ - referencing a prof-block attribute array
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_ktProf_(man)
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer,pointer,dimension(:) :: ptr_ktProf_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_ktProf_'

  if(.not.associated(man%ktProf))	&
	call die(myname_,'%ktProf undefined')

  ptr_ktProf_ => man%ktProf

end function ptr_ktProf_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_kxProf_ - referencing a prof-block attribute array
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_kxProf_(man)
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer,pointer,dimension(:) :: ptr_kxProf_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_kxProf_'

  if(.not.associated(man%kxProf))	&
	call die(myname_,'%kxProf undefined')

  ptr_kxProf_ => man%kxProf

end function ptr_kxProf_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_ksProf_ - referencing a prof-block attribute array
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_ksProf_(man)
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      integer,pointer,dimension(:) :: ptr_ksProf_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_ksProf_'

  if(.not.associated(man%ksProf))	&
	call die(myname_,'%ksProf undefined')

  ptr_ksProf_ => man%ksProf

end function ptr_ksProf_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_latProf_ - referencing a prof-block attribute array
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_latProf_(man)
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      real,pointer,dimension(:) :: ptr_latProf_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_latProf_'

  if(.not.associated(man%latProf))	&
	call die(myname_,'%latProf undefined')

  ptr_latProf_ => man%latProf

end function ptr_latProf_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_lonProf_ - referencing a prof-block attribute array
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_lonProf_(man)
      use m_die,only : die
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      real,pointer,dimension(:) :: ptr_lonProf_

! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_lonProf_'

  if(.not.associated(man%lonProf))	&
	call die(myname_,'%lonProf undefined')

  ptr_lonProf_ => man%lonProf

end function ptr_lonProf_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_sndxNav_ - referencing the sndx-block navigator
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_sndxNav_(man)
      use m_Navigator,only : Navigator
      implicit none
      type(MultiAccessNavigator),intent(in) :: man
      type(Navigator),pointer :: ptr_sndxNav_


! !REVISION HISTORY:
! 	08Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_sndxNav_'

  ptr_sndxNav_ => man%sndxNav

end function ptr_sndxNav_

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
      type(MultiAccessNavigator),pointer :: mold
      integer,optional,intent(out) :: stat
      type(MultiAccessNavigator),pointer :: new_

! !REVISION HISTORY:
! 	28Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::new_'
  type(MultiAccessNavigator),pointer :: obj
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
! !IROUTINE: delete__ - delete an object for a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine delete__(obj,stat)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(MultiAccessNavigator),pointer :: obj
      integer,optional,intent(out) :: stat


! !REVISION HISTORY:
! 	28Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::delete__'
  integer :: ier

  if(present(stat)) stat=0

	if(mall_ison()) call mall_co(1,myname)

  deallocate(obj,stat=ier)
	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_,'deallocate()',ier)
	  stat=ier
	  return
	endif

end subroutine delete__

end module m_MultiAccessNavigator

