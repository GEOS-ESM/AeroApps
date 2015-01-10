!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Regioner - Everything about the "regions"
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Regioner
      implicit none
      private	! except

      public :: Regioner_init
      public :: Regioner_clean
      public :: Regioner_part
      public :: Regioner_set

      interface Regioner_init ; module procedure init_ ; end interface
      interface Regioner_clean; module procedure clean_; end interface
      interface Regioner_part
         module procedure part_ 
         module procedure part_all_
      end interface
      interface Regioner_set ; module procedure	&
	setks_,	&
	setll_
      end interface

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!       31Oct01 - Tom Clune <clune@sgi.com>
!               . Modified to use m_Spherical_Partition
!       05mar02 - Tom Clune
!               . brought m_Regioner::part_all_ from m_Regioner.f90
!                 (note the lower case "f")
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Regioner'

  logical,save :: Regioner_defined=.false.
  logical,save :: Partition_init=.false.
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - Initialize the definition of the "regions"
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(nxkr)
      use m_Spherical_Partition, only : Initialize_Partition => Initialize
      use m_Spherical_Partition, only : NumberOfRegions
      use m_die, only : die
      use m_mpout, only : mpout_log, mpout
      implicit none
      integer,intent(out) :: nxkr

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

  integer :: refinement_level, ier
  Character(len=255) :: env_val ! buffer for environment variables

  if(Regioner_defined) call die(myname_,'multiple object definition')

  If(.not. Partition_init) Then ! only do this once

!  Read environment to determine the refinement level for regions
!  -------------------------------------------------------------------
     Call getenv('REFINEMENT_LEVEL',env_val)

     if ( len_trim(env_val) .gt. 0 ) then
        read(env_val,*,IOSTAT=ier) refinement_level
        if(ier /= 0) call die(myname_,'multiple object definition')
     else
        refinement_level = 0 ! 80 regions by default
     end if

     Call mpout_log(myname_,'Using refinement level',refinement_level)
     Call Initialize_Partition(n_levels = refinement_level, compress=.false.)
     Partition_init = .true.
     call seticos()


  End If

  nxkr=NumberOfRegions()


	! That is easy ..

  Regioner_defined=.true.

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - In case the region definition becoming complicate
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      use m_die, only : die
      implicit none

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  if(.not.Regioner_defined) call die(myname_,'object undefined')

  Regioner_defined=.false.
end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: part_ - An "as-is" partition scheme
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine part_(nxkr,nxkt,lc_kr,ln_kr,lc_kt,ln_kt,	&
	nobs,rlat,rlon,rlev,tims,kxs,kts,dels		)

      use m_mpout,   only : mpout, mpout_log
      use m_die,     only : die
      use m_mall,    only : mall_ison,mall_mci,mall_mco
      implicit none

      integer,intent(in) :: nxkr
      integer,intent(in) :: nxkt
      integer,dimension(:)  ,intent(out) :: lc_kr
      integer,dimension(:)  ,intent(out) :: ln_kr
      integer,dimension(:,:),intent(out) :: lc_kt
      integer,dimension(:,:),intent(out) :: ln_kt

      integer,intent(in) :: nobs
      real   ,dimension(:),intent(inout) :: rlat
      real   ,dimension(:),intent(inout) :: rlon
      real   ,dimension(:),intent(inout) :: rlev
      real   ,dimension(:),intent(inout) :: tims
      integer,dimension(:),intent(inout) :: kxs
      integer,dimension(:),intent(inout) :: kts
      real   ,dimension(:),intent(inout) :: dels

! !REVISION HISTORY:
! 	27Oct99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::part_'

  logical,parameter :: verbose=.false.
  integer :: ier
  integer :: kr,kt,lc

  real,allocatable,dimension(:) :: sigU,sigO,sigF

  if(.not.Regioner_defined) call die(myname_,'object undefined')

	allocate(sigU(nobs),sigO(nobs),sigF(nobs), stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(sigU,myname)
		  call mall_mci(sigO,myname)
		  call mall_mci(sigF,myname)
		endif

	!..Sort observations in the order of:
	!
	!	region(lat,lon)-kt-kx-lat-lon-pres
	!
	!  Also, define pointer/size information of each region and type
	!  by set arrays iregbeg, ireglen, and ityplen.

  call sort (	myname_,verbose,mpout,			&
	nobs,	rlat,	rlon,	rlev,	kxs,	kts,	&
	dels,	sigU,	sigO,	sigF,	tims,		&
	nxkr,	nxkt,	lc_kr,	ln_kr,	ln_kt		)

	! define lc_kt(:,:) for the offset of each type block (kt)
	! with respect to the initial location of the regions (kr).

  do kr=1,nxkr
    lc=0
    do kt=1,nxkt
      lc_kt(kt,kr)=lc
      lc=lc+ln_kt(kt,kr)
    end do
  end do

		if(mall_ison()) then
		  call mall_mco(sigU,myname)
		  call mall_mco(sigO,myname)
		  call mall_mco(sigF,myname)
		endif
	deallocate(sigU,sigO,sigF,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine part_

      subroutine part_all_(nxkr,nxkt,lc_kr,ln_kr,lc_kt,ln_kt,       &
      nobs,rlat,rlon,rlev,kxs,kts,kid                               )
!******************************************************************************
!*                                                                            *
!*                                                                            *
!*                                                                            *
!*                                                                            *
!******************************************************************************
      use m_mpout,   only : mpout
      use m_die,     only : die
      

      integer,intent(in) :: nxkr
      integer,intent(in) :: nxkt
      integer,dimension(:)  ,intent(out) :: lc_kr
      integer,dimension(:)  ,intent(out) :: ln_kr
      integer,dimension(:,:),intent(out) :: lc_kt
      integer,dimension(:,:),intent(out) :: ln_kt

      integer,intent(in) :: nobs
      real   ,dimension(:),intent(inout) :: rlat
      real   ,dimension(:),intent(inout) :: rlon
      real   ,dimension(:),intent(inout) :: rlev
      integer,dimension(:),intent(inout) :: kxs
      integer,dimension(:),intent(inout) :: kts
      integer,dimension(:),intent(inout) :: kid

      character(len=*),parameter :: myname_=myname//'::part_'

      logical,parameter :: verbose=.false.
      integer :: ier
      integer :: kr,kt,lc

      if(.not.Regioner_defined) call die(myname_,'object undefined')

!..Sort observations 
      call sort_all(myname_,verbose,mpout,nobs,rlat,rlon,rlev,kxs,kts, &
                    kid,nxkr,nxkt,lc_kr,ln_kr,ln_kt)

! define lc_kt(:,:) for the offset of each type block (kt)
! with respect to the initial location of the regions (kr).

      do kr=1,nxkr
      lc=0
      do kt=1,nxkt
      lc_kt(kt,kr)=lc
      lc=lc+ln_kt(kt,kr)
      end do
      end do

      end subroutine part_all_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setks_ - set region IDs for soundings
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setks_(n,kr,kx,ks,rlat,rlon)
      use m_geometry, only : ll2xyz
      use m_Spherical_Partition, only : LatLon_to_region
      use m_die,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none

      integer,intent(in) :: n
      integer,dimension(:),intent(out) :: kr
      integer,dimension(:),intent(in) :: kx
      integer,dimension(:),intent(in) :: ks
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon

! !REVISION HISTORY:
! 	23Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setks_'

  integer :: nsnd,i,ier,ksi,kxi
  integer,allocatable,dimension(:) :: ksnd,kreg
  real,   allocatable,dimension(:) :: slat,slon
  integer :: nxkr
  logical :: local_Regioner

  if(n<=0) return

  local_Regioner = .not. Regioner_defined
  if(local_Regioner) call init_(nxkr)

	allocate( ksnd(n),kreg(n),slat(n),slon(n), stat=ier	)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(ksnd,myname)
		  call mall_mci(kreg,myname)
		  call mall_mci(slat,myname)
		  call mall_mci(slon,myname)
		endif

  nsnd=0
  i=1
  nsnd=nsnd+1
  slat(nsnd)=rlat(i)
  slon(nsnd)=rlon(i)
  kxi=kx(i)
  ksi=ks(i)
  ksnd(i)=nsnd

  do i=2,n
    if(kx(i)/=kxi .or. ks(i)/=ksi) then
      nsnd=nsnd+1
      slat(nsnd)=rlat(i)
      slon(nsnd)=rlon(i)
      kxi=kx(i)
      ksi=ks(i)
    endif
    ksnd(i)=nsnd
  end do

  call LatLon_to_Region(nsnd,slat,slon,kreg,ier,atlevel=0)
	if(ier/=0) call die(myname_,'LatLon_to_Region()',ier)

  do i=1,n
    kr(i)=kreg(ksnd(i))
  end do

		if(mall_ison()) then
		  call mall_mco(ksnd,myname)
		  call mall_mco(kreg,myname)
		  call mall_mco(slat,myname)
		  call mall_mco(slon,myname)
		endif
	deallocate(ksnd,kreg,slat,slon,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

  if(local_Regioner) call clean_()

end subroutine setks_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setll_ - set region IDs for lat-long datasets
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setll_(n,kr,rlat,rlon)
      use m_Spherical_Partition, only : xyz2reg
      use m_die,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none
      integer,intent(in) :: n
      integer,dimension(:),intent(out) :: kr
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon

! !REVISION HISTORY:
! 	23Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setll_'

  integer :: ier
  real,   allocatable,dimension(:) :: x,y,z
  integer :: nxkr
  logical :: local_Regioner

  if(n<=0) return

  local_Regioner = .not. Regioner_defined
  if(local_Regioner) call init_(nxkr)

	allocate( x(n),y(n),z(n),	stat=ier	)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(x,myname)
		  call mall_mci(y,myname)
		  call mall_mci(z,myname)
		endif

  call ll2xyz(rlon,rlat, n, x,y,z, ier)
	if(ier/=0) call die(myname_,'ll2xyz()',ier)

!   Only region onto the base 80 regions for now
!   Further refinement comes later.
!   --------------------------------------------
  call xyz2reg(n,x,y,z,kr,atlevel=0)

		if(mall_ison()) then
		  call mall_mco(x,myname)
		  call mall_mco(y,myname)
		  call mall_mco(z,myname)
		endif
	deallocate(x,y,z,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

  if(local_Regioner) call clean_()

end subroutine setll_
end module m_Regioner
