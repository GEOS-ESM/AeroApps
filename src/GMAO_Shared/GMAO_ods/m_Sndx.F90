!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Sndx - Sounding index generator
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Sndx
      use m_odsmeta, only : kxmax
      use m_odsmeta, only : kxmod
      use m_die,     only : perr
      implicit none
      private	! except

      public :: setSndx		! Set sounding index values
      public :: Sndx2pix	! another sounding index (pix)
      public :: dumpks

      interface setSndx; module procedure	&
	setstn_,	&
	setkx_,		&
	setll_,		&
	set_; end interface


! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Sndx'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setstn_ - set sounding index by stn. IDs.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setstn_(ks,kx,stn)
      use m_SortingTools,only : indexSet,indexSort
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none
      integer         ,dimension(:),intent(out) :: ks	! sounding index
      integer         ,dimension(:),intent(in ) :: kx	! instrument IDs
      character(len=*),dimension(:),intent(in ) :: stn	! station IDs

! REVISION HISTORY:
!      22Apr04	- Jing Guo <guo@dao.gsfc.nasa.gov>
!      	- initial prototype/prolog/code
!EP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setstn_'
  integer :: ikx,iks,i,j
  character(len=len(stn)) :: xID
  integer,allocatable,dimension(:) :: indx
  integer :: lobs
  integer :: ier

  lobs=size(kx)
	if(lobs/=size(ks) .or. lobs/=size(stn)) then
	  call perr(myname_,'size( ks)',size( ks))
	  call perr(myname_,'size( kx)',size( kx))
	  call perr(myname_,'size(stn)',size(stn))
	  call die(myname_,'mismatched argument sizes',ier)
	endif

  if(lobs==0) return

	allocate(indx(lobs), stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(indx,myname)

  	! local sorting by kx and by station IDs.

  call indexSet(indx)
  call indexSort(indx,stn)
  call indexSort(indx,kx)

  iks=1

  i=1
  j=indx(i)

  ikx=kx(j)
  xid=stn(j)
  ks(j)=iks

  	! Set all sounding index values.  All soundings of the same
	! instrument are assumed to be on the same PE.

  do i=2,lobs
    j=indx(i)	! refer to the actual location

    if(kx(j)/=ikx) then
    	! if this is a new intrument type (kx), reset sounding index
	! counter (iks) and instrument marker (ikx) to the new value.

      iks=1
      ikx=kx(j)
      xID=stn(j)

    else
    	! if this is the same intrument type as the previous entry,

      if(stn(j)/=xID) then
      	! if this is a new stn ID, reset sounding index (ks).
        xID=stn(j)
        iks=iks+1
      endif

    endif

    ks(j)=iks
  end do

		if(mall_ison()) call mall_mco(indx,myname)
	deallocate(indx, stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine setstn_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setkx_ - assign sounding indices
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setkx_(ks,kx,rlat,rlon)
      use m_SortingTools,only : indexSet,indexSort
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none
      integer,dimension(:),intent(out) :: ks
      integer,dimension(:),intent(in ) :: kx
      real   ,dimension(:),intent(in ) :: rlat
      real   ,dimension(:),intent(in ) :: rlon

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setkx_'
  integer :: ndat
  integer :: ier
  integer :: i,l
  integer :: kx_i
  integer :: kx_l
  real    :: rlon_l
  real    :: rlat_l
  integer :: ks_l
  integer,allocatable,dimension(:) :: indx

  ndat=size(ks)
  if(ndat<=0) return

	allocate(indx(ndat),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(indx,myname)

  call indexSet(indx)
  call indexSort(indx,rlon)
  call indexSort(indx,rlat)
  call indexSort(indx,kx)

!________________________________________
	! Segmentation using indices

  l=indx(1)	! Process the initial data point.
  kx_i=kx(l)

		! Count as another sounding and save the current data.
  ks_l=1

  kx_l=kx_i
  rlon_l=rlon(l)
  rlat_l=rlat(l)

  ks(l)=ks_l	! Save the current sounding index.

		! Process rest data points.
  do i=2,ndat
    l=indx(i)
    kx_i=kx(l)

		! If the saved data are not the same as the current,
		! count as another sounding and save the current data.

    if(kx_l/=kx_i .or. rlon_l/=rlon(l) .or. rlat_l/=rlat(l)) then

		! The following statement relies on the sorting order
		! of the data to be primerily by kx.  The purpose is
		! to restart a counter for every kx.

      if(kx_l/=kx_i) ks_l=0
      ks_l=ks_l+1

      kx_l=kx_i
      rlon_l=rlon(l)
      rlat_l=rlat(l)
    endif

		! Save the current sounding index.

    ks(l)=ks_l
  end do

		if(mall_ison()) call mall_mco(indx,myname)
	deallocate(indx,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine setkx_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: set_ - assign sounding indices
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine set_(ks,rlat,rlon)
      use m_SortingTools,only : indexSet,indexSort
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none
      integer,dimension(:),intent(out) :: ks
      real   ,dimension(:),intent(in ) :: rlat
      real   ,dimension(:),intent(in ) :: rlon

! !REVISION HISTORY:
! 	05Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::set_'
  integer :: ndat
  integer :: ier
  integer :: i,l
  real    :: rlon_l
  real    :: rlat_l
  integer :: ks_l
  integer,allocatable,dimension(:) :: indx

  ndat=size(ks)
  if(ndat<=0) return

	allocate(indx(ndat),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(indx,myname)

  call indexSet(indx)
  call indexSort(indx,rlon)
  call indexSort(indx,rlat)

!________________________________________
	! Segmentation using indices

  l=indx(1)	! Process the initial data point.

		! Count as another sounding and save the current data.
  ks_l=1

  rlon_l=rlon(l)
  rlat_l=rlat(l)

  ks(l)=ks_l	! Save the current sounding index.

		! Process rest data points.
  do i=2,ndat
    l=indx(i)

		! If the saved data are not the same as the current,
		! count as another sounding and save the current data.

    if(rlon_l/=rlon(l) .or. rlat_l/=rlat(l)) then

      ks_l=ks_l+1

      rlon_l=rlon(l)
      rlat_l=rlat(l)
    endif

		! Save the current sounding index.

    ks(l)=ks_l
  end do

		if(mall_ison()) call mall_mco(indx,myname)
	deallocate(indx,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine set_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setll_ - set sounding indices locally from lat-lon values
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setll_(ks,rlat,dlat,rlon,dlon)
      use m_die,only : die
      implicit none
      integer,dimension(:),intent(out) :: ks

      real   ,dimension(:),intent(in) :: rlat
      real		  ,intent(in) :: dlat

      real   ,dimension(:),intent(in) :: rlon
      real		  ,intent(in) :: dlon

! !REVISION HISTORY:
! 	03Nov00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setll_'
  integer :: i
  integer :: nlon,ilon,ilat
  real :: alat,alon

  alat=abs(dlat)
  alon=abs(dlon)

  if(alat*alon < .0001) call die(myname_,'unreliable configuration')

  nlon=360./alon
  do i=1,size(ks)
    ilat=nint(rlat(i)/alat)
    ilon=modulo(nint(rlon(i)/alon),nlon)
    ks(i)=ilat*nlon+ilon
  end do

end subroutine setll_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Sndx2pix -
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine Sndx2pix(ks,kx,rlat,rlon)
      implicit none
      integer,dimension(:),intent(inout) :: ks
      integer,dimension(:),intent(in) :: kx
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon

! !REVISION HISTORY:
! 	22Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!	04Jun04 - Todling  removed kxmax.h (using m_odsmeta)
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Sndx2pix'

  integer :: i,nobs
  integer :: ks_i,kx_i
  real    :: rlat_i,rlon_i
  integer :: knt(kxmax)
  integer :: pix_now
  logical :: begin_prof
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  nobs=size(ks)
  if (nobs.le.0) return

  knt(1:kxmax)=-1		! none yet

	! this part needs to be serially processed, since it is
	! it is upto the previous entry to determine the pix
	! value of the current entry.

  kx_i=kx(1)
  ks_i=ks(1)
  rlat_i=rlat(1)
  rlon_i=rlon(1)

  knt(kx_i)=knt(kx_i)+1
  pix_now=kx_i+knt(kx_i)*kxmod

  ks(1)=pix_now

  do i=2,nobs

    begin_prof=	ks_i.ne.ks(i)		.or.		&
		kx_i.ne.kx(i)

    kx_i=kx(i)
    ks_i=ks(i)
    rlat_i=rlat(i)
    rlon_i=rlon(i)

    if(begin_prof) then		! If a new profile, count it
      knt(kx_i)=knt(kx_i)+1
      pix_now=kx_i+knt(kx_i)*kxmod
    endif

    ks(i)=pix_now
  end do
!_______________________________________________________________________
end subroutine Sndx2pix
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dumpks -
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine dumpks(fname,kx,ks,rlat,rlon)
      use m_ioutil,only : luavail,opntext,clstext
      implicit none
      character(len=*),intent(in) :: fname
      integer,dimension(:),intent(in) :: kx
      integer,dimension(:),intent(in) :: ks
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon

! !REVISION HISTORY:
! 	22Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::dumpks'

  integer :: lu,ier,i

  lu=luavail()
  call opntext(lu,fname,'unknown',ier)

  do i=1,size(kx)
    write(lu,'(3i8,2f8.2)') i,kx(i),ks(i),rlat(i),rlon(i)
  end do

  call clstext(lu,ier)

end subroutine dumpks
end module m_Sndx
