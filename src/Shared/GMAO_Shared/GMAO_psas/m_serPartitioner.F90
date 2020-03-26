!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_serPartitioner - A partitioner
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_serPartitioner
      use m_Navigator,only : Navigator
      implicit none
      private	! except

      public :: serPartitioner		! The class data structure

      public :: serPartitioner_init
      public :: serPartitioner_clean

      public :: ptr_Navigator
      public :: ptr_krNav
      public :: ptr_ktNav
      public :: ptr_indx

    type serPartitioner
      private
      type(Navigator),pointer :: nav
      integer,pointer,dimension(:) :: krNav
      integer,pointer,dimension(:) :: ktNav

      integer,pointer,dimension(:) :: indx
    end type serPartitioner

      interface serPartitioner_init ; module procedure	&
	initob_,	&
	initai_; end interface
      interface serPartitioner_clean; module procedure	&
	clean_; end interface

      interface ptr_Navigator; module procedure	&
	ptr_Navigator_; end interface
      interface ptr_krNav; module procedure ptr_krNav_; end interface
      interface ptr_ktNav; module procedure ptr_ktNav_; end interface
      interface ptr_indx ; module procedure ptr_indx_ ; end interface

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_serPartitioner'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initob_ - Initialize an object for observations
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initob_(part, ndat,rlat,rlon,rlev,rtim,kx,ks,kt)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none

      type(serPartitioner),intent(out) :: part

      integer,intent(in)  :: ndat
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon
      real   ,dimension(:),intent(in) :: rlev
      real   ,dimension(:),intent(in) :: rtim
      integer,dimension(:),intent(in) :: kx
      integer,dimension(:),intent(in) :: ks
      integer,dimension(:),intent(in) :: kt

! !REVISION HISTORY:
!	26May00	- Jing Guo
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initob_'
  integer :: nxkr
  integer :: ier
  integer,allocatable,dimension(:) :: kr

	allocate(kr(ndat),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(kr,myname)

  allocate(part%indx(ndat),stat=ier)
	if(ier/=0) call die(myname_,'allocate(%indx)',ier)
	if(mall_ison()) call mall_mci(part%indx,myname)

  call sortob_(ndat,rlat,rlon,rlev,rtim,kx,ks,kt, kr,part%indx)

  call locateix_(part%nav,part%krNav,part%ktNav,	&
	kr,kt,part%indx)

		if(mall_ison()) call mall_mco(kr,myname)
	deallocate(kr,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine initob_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initai_ - Initialize an object for analysis increments
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initai_(part, ndat,rlat,rlon,rlev,kt)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none

      type(serPartitioner),intent(out) :: part

      integer,intent(in)  :: ndat
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon
      real   ,dimension(:),intent(in) :: rlev
      integer,dimension(:),intent(in) :: kt

! !REVISION HISTORY:
!	26May00	- Jing Guo
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initai_'
  integer :: nxkr
  integer :: ier
  integer,allocatable,dimension(:) :: kr

	allocate(kr(ndat),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(kr,myname)

  allocate(part%indx(ndat),stat=ier)
	if(ier/=0) call die(myname_,'allocate(%indx)',ier)
	if(mall_ison()) call mall_mci(part%indx,myname)

  call sortai_(ndat,rlat,rlon,rlev,kt, kr,part%indx)

  call locateix_(part%nav,part%krNav,part%ktNav,	&
	kr,kt,part%indx)

		if(mall_ison()) call mall_mco(kr,myname)
	deallocate(kr,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine initai_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - In case the region definition becoming complicate
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(part)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mco,mall_co
      use m_Navigator,only : clean
      implicit none
      type(serPartitioner),intent(inout) :: part

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  call clean(part%nav)

	if(mall_ison()) then
	  call mall_mco(part%krNav,myname)
	  call mall_mco(part%ktNav,myname)
	  call mall_mco(part%indx,myname)
	  call mall_co(1,myname)	! for part%nav
	endif

  deallocate(part%krNav,part%ktNav,part%indx,part%nav,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sortob_ - An "as-is" sorting scheme with a new data struct.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine sortob_(ndat,rlat,rlon,rlev,rtim,kx,ks,kt, kr,indx)

      use m_die,     only : die
      use m_SortingTools, only : indexSet
      use m_SortingTools, only : indexSort
      use m_Regioner,only : Regioner_set
      implicit none

      integer,       intent(in)  :: ndat
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon
      real   ,dimension(:),intent(in) :: rlev
      real   ,dimension(:),intent(in) :: rtim
      integer,dimension(:),intent(in) :: kx
      integer,dimension(:),intent(in) :: ks
      integer,dimension(:),intent(in) :: kt

      integer,dimension(:),intent(out) :: kr
      integer,dimension(:),intent(out) :: indx

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sortob_'

  call Regioner_set(ndat,kr,kx,ks,rlat,rlon)

  call indexSet(ndat,indx)

	!..Sort data in the order of:
	!
	!	region(lat,lon)-kt-kx-ks-pres

  call indexSort(ndat,indx,rlev,descend=.true. )
  call indexSort(ndat,indx,ks,  descend=.false.)
  call indexSort(ndat,indx,kx,  descend=.false.)
  call indexSort(ndat,indx,kt,  descend=.false.)
  call indexSort(ndat,indx,kr,  descend=.false.)

  

end subroutine sortob_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sortai_ - An "as-is" sorting scheme with a new data struct.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine sortai_(ndat,rlat,rlon,rlev,kt, kr,indx)

      use m_die,     only : die
      use m_SortingTools, only : indexSet
      use m_SortingTools, only : indexSort
      use m_Regioner,only : Regioner_set
      implicit none

      integer,       intent(in)  :: ndat
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon
      real   ,dimension(:),intent(in) :: rlev
      integer,dimension(:),intent(in) :: kt

      integer,dimension(:),intent(out) :: kr
      integer,dimension(:),intent(out) :: indx

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sortai_'

  call Regioner_set(ndat,kr,rlat,rlon)

  call indexSet(ndat,indx)

	!..Sort data in the order of:
	!
	!	region(lat,lon)-kt-lat-lon-pres

  call indexSort(ndat,indx,rlev,descend=.true. )
  call indexSort(ndat,indx,rlon,descend=.false. )
  call indexSort(ndat,indx,rlat,descend=.false. )
  call indexSort(ndat,indx,kt,  descend=.false.)
  call indexSort(ndat,indx,kr,  descend=.false.)

end subroutine sortai_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_Navigator_ - referering to the Navigator
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_Navigator_(part)
      use m_Navigator,only : Navigator
      implicit none
      type(serPartitioner),intent(in) :: part
      type(Navigator),pointer :: ptr_Navigator_

! !REVISION HISTORY:
! 	26May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_Navigator_'

  ptr_Navigator_ => part%nav

end function ptr_Navigator_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_krNav_ - referering to the kr-table
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_krNav_(part)
      implicit none
      type(serPartitioner),intent(in) :: part
      integer,pointer,dimension(:) :: ptr_krNav_

! !REVISION HISTORY:
! 	26May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_krNav_'

  ptr_krNav_ => part%krNav

end function ptr_krNav_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_ktNav_ - referering to the kt-table
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_ktNav_(part)
      implicit none
      type(serPartitioner),intent(in) :: part
      integer,pointer,dimension(:) :: ptr_ktNav_

! !REVISION HISTORY:
! 	26May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_ktNav_'

  ptr_ktNav_ => part%ktNav

end function ptr_ktNav_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_indx_ - refereing to the permutation index table
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_indx_(part)
      implicit none
      type(serPartitioner),intent(in) :: part
      integer,pointer,dimension(:) :: ptr_indx_

! !REVISION HISTORY:
! 	26May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_indx_'

  ptr_indx_ => part%indx

end function ptr_indx_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: locateix_ - locate indexed kr-kt segments
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine locateix_(nav,kxNav,kyNav,kx,ky,indx)
      use m_Navigator,only : Navigator
      use m_Navigator,only : Navigator_init
      use m_Navigator,only : ptr_displs
      use m_Navigator,only : ptr_counts
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci,mall_ci,mall_mco
      implicit none
      type(Navigator),pointer :: nav
      integer,pointer,dimension(:) :: kxNav
      integer,pointer,dimension(:) :: kyNav

      integer,dimension(:),intent(in) :: kx
      integer,dimension(:),intent(in) :: ky
      integer,dimension(:),intent(in) :: indx

! !REVISION HISTORY:
! 	23Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::locateix_'

  integer :: ln,i,l,n,kxi,kyi,nseg,ier
  integer,allocatable,dimension(:) :: lns
  integer,pointer,dimension(:) :: counts
  integer,pointer,dimension(:) :: displs

  n=size(indx)
  nseg=0

	allocate(lns(n),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(lns,myname)

  if(n>0) then

    l=indx(1)
    kxi=kx(l)
    kyi=ky(l)

    ln=1
    nseg=nseg+1

    do i=2,n
      l=indx(i)

      if(kxi/=kx(l) .or. kyi/=ky(l)) then	! this is a new segment
        lns(nseg)=ln	! Save the size of the last segment
        nseg=nseg+1	! Count in the new segment
        ln=0

        kxi=kx(l)
        kyi=ky(l)
      endif

      ln=ln+1
    end do

    lns(nseg)=ln
  endif

  allocate(nav,stat=ier)
	if(ier/=0) call die(myname_,'allocate(nav)',ier)
	if(mall_ison()) call mall_ci(1,myname)

  call Navigator_init(nav,nseg,stat=ier)
	if(ier/=0) call die(myname_,'Navigator_init()',ier)

  allocate(kxNav(nseg),kyNav(nseg),stat=ier)
	if(ier/=0) call die(myname_,'allocate(kxNav,kyNav)',ier)
	if(mall_ison()) then
	  call mall_mci(kxNav,myname)
	  call mall_mci(kyNav,myname)
	endif

  counts => ptr_counts(nav)
  displs => ptr_displs(nav)

  counts(:) = lns(1:nseg)
  displs(1)=0
  l=indx(1)
  kxNav(1)=kx(l)
  kyNav(1)=ky(l)
  do i=2,nseg
    displs(i)=displs(i-1)+counts(i-1)
    l=indx(displs(i)+1)
    kxNav(i)=kx(l)
    kyNav(i)=ky(l)
  end do

  nullify(counts)
  nullify(displs)

		if(mall_ison()) call mall_mco(lns,myname)
	deallocate(lns,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine locateix_
end module m_serPartitioner
