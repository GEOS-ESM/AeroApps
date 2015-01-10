!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AttributesMAN - MAN builder from given Attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_AttributesMAN
      implicit none
      private	! except

      public :: build_fromAttr		! build a MAN from Attributes

      interface build_fromAttr; module procedure build_; end interface

! !REVISION HISTORY:
! 	06Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AttributesMAN'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: build_ - build a MAN from given Attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine build_(man,indx,attr)
      use m_die,only : die
      use m_MultiAccessNavigator,only : MultiAccessNavigator
      use m_MultiAccessNavigator,only : build_fromSorted

      use m_Attributes,only : Attributes
      use m_Attributes,only : key
      use m_Attributes,only : OBS_SUBSET
      use m_Attributes,only : FLD_SUBSET

      use m_Attributes,only : ptr_kr
      use m_Attributes,only : ptr_kt
      use m_Attributes,only : ptr_kx
      use m_Attributes,only : ptr_ks
      use m_Attributes,only : ptr_lat
      use m_Attributes,only : ptr_lon
      use m_Attributes,only : ptr_lev

      use m_SortingTools,only : indexSet
      use m_SortingTools,only : indexSort

      implicit none
      type(MultiAccessNavigator),intent(out) :: man
      integer,dimension(:)      ,intent(out) :: indx
      type(Attributes)          ,intent(in ) :: attr

! !REVISION HISTORY:
! 	06Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::build_'
  integer,pointer,dimension(:) :: kr
  integer,pointer,dimension(:) :: kt
  integer,pointer,dimension(:) :: kx
  integer,pointer,dimension(:) :: ks
  real   ,pointer,dimension(:) :: lat
  real   ,pointer,dimension(:) :: lon
  integer :: key_

	! Initialize sorting indices

  call indexSet(indx)

	! Sort the index by attr%lev in an ascending order

  call indexSort(indx,ptr_lev(attr),descend=.true.)

  kr => ptr_kr(attr)
  kt => ptr_kt(attr)

  key_=key(attr)

  select case(key_)
  case(OBS_SUBSET)

    kx => ptr_kx(attr)
    ks => ptr_ks(attr)

		! Sort in kr-kt-kx-ks order

    call indexSort(indx,ks,descend=.false.)
    call indexSort(indx,kx,descend=.false.)
    call indexSort(indx,kt,descend=.false.)
    call indexSort(indx,kr,descend=.false.)

    call build_fromSorted(man,			&
	 kr= kr(indx(:)), kt= kt(indx(:)),	&
	 kx= kx(indx(:)), ks= ks(indx(:))	)

    nullify(ks)
    nullify(kx)

  case(FLD_SUBSET)

    lat => ptr_lat(attr)
    lon => ptr_lon(attr)

		! Sort in kr-kt-lat-lon order

    call indexSort(indx,lon,descend=.false.)
    call indexSort(indx,lat,descend=.false.)
    call indexSort(indx, kt,descend=.false.)
    call indexSort(indx, kr,descend=.false.)

    call build_fromSorted(man,			&
	 kr= kr(indx(:)), kt= kt(indx(:)),	&
	lat=lat(indx(:)),lon=lon(indx(:))	)

    nullify(lon)
    nullify(lat)

  case default
    call die(myname_,'invalid attributes key',key_)
  end select

  nullify(kt)
  nullify(kr)

end subroutine build_

end module m_AttributesMAN
