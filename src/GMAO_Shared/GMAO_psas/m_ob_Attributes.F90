!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ob_Attributes - Attributes of the obs. vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ob_Attributes
      implicit none
      private	! except

      public :: ob_Attributes_init
      public :: ob_Attributes_clean

      public :: ob_kx, ob_ks, ob_kt
      public :: ob_lat,ob_lon,ob_lev
      public :: ob_del,ob_xvec

      interface ob_Attributes_init; module procedure	&
	init_
      end interface

      interface ob_Attributes_clean; module procedure	&
	clean_
      end interface

! !REVISION HISTORY:
! 	03Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ob_Attributes'

	! Observation attributes (nobs)

  integer,target,allocatable,dimension(:) :: ob_kx,ob_ks,ob_kt
  real   ,target,allocatable,dimension(:) :: ob_lat,ob_lon,ob_lev
  real   ,target,allocatable,dimension(:,:) :: ob_del
  real   ,target,allocatable,dimension(:,:) :: ob_xvec

  logical,save :: ob_Attributes_defined=.false.

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize the data object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(nobs,indx,kxs,kss,kts,rlat,rlon,rlev,dels)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      integer,intent(in) :: nobs
      integer,dimension(:),intent(in) :: indx
      integer,dimension(:),intent(in) :: kxs
      integer,dimension(:),intent(in) :: kss
      integer,dimension(:),intent(in) :: kts
      real   ,dimension(:),intent(in) :: rlat
      real   ,dimension(:),intent(in) :: rlon
      real   ,dimension(:),intent(in) :: rlev
      real   ,dimension(:),intent(in) :: dels

! !REVISION HISTORY:
! 	06Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier
  integer :: i,l

		! Sorting observations (-forecasts) is taking place now.

  if(ob_Attributes_defined) call die(myname_,'multiple definitions')

  allocate( ob_kx(nobs),  ob_ks(nobs),  ob_kt(nobs),	&
	    ob_lat(nobs), ob_lon(nobs), ob_lev(nobs),	&
	    ob_del(nobs,1),ob_xvec(nobs,1), stat=ier	)

	if(ier/=0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(ob_kx, myname)
	  call mall_mci(ob_ks, myname)
	  call mall_mci(ob_kt, myname)
	  call mall_mci(ob_lat,myname)
	  call mall_mci(ob_lon,myname)
	  call mall_mci(ob_lev,myname)
	  call mall_mci(ob_del,myname)
	  call mall_mci(ob_xvec,myname)
	endif

  do i=1,nobs
    l=indx(i)
    ob_kx(i)=kxs(l)
    ob_ks(i)=kss(l)
    ob_kt(i)=kts(l)

    ob_lat(i)=rlat(l)
    ob_lon(i)=rlon(l)
    ob_lev(i)=rlev(l)

    ob_del(i,1)=dels(l)
  end do

  ob_Attributes_defined=.true.

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the data object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      use m_die, only : die
      use m_mall,only : mall_ison,mall_mco
      implicit none

! !REVISION HISTORY:
! 	06Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  if(.not.ob_Attributes_defined) call die(myname_,'object undefined')

	if(mall_ison()) then
	  call mall_mco(ob_kx, myname)
	  call mall_mco(ob_ks, myname)
	  call mall_mco(ob_kt, myname)
	  call mall_mco(ob_lat,myname)
	  call mall_mco(ob_lon,myname)
	  call mall_mco(ob_lev,myname)
	  call mall_mco(ob_del,myname)
	  call mall_mco(ob_xvec,myname)
	endif

  deallocate(ob_kx,ob_ks,ob_kt,ob_lat,ob_lon,ob_lev,	&
	ob_del,ob_xvec,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

  ob_Attributes_defined=.false.

end subroutine clean_
end module m_ob_Attributes
