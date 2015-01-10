!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_mvcorF_nsep - Data module of non-sperable mvcorF
!
! !DESCRIPTION:
!
!   m_mvcorF_nsep defines the storage components required for the windowed
!   correlation function and its first two derivatives.  It also defines
!   storage components for the squared length scales rLav2 and the
!   multiplier CLfact required for the cross-correlation.
!
!   m_mvcorF_nsep also defines the storage components for the vertical
!   correlations.  Notice vfecHD and vfecDD are not normalized by
!   themself, but factorized with horizontal correlation components to
!   normalize the total correlation coefficients.  However, the
!   assumption is that vfecHH and hfecHH are determined from normalized
!   functions.
!
! !INTERFACE:

    module m_mvcorF_nsep 
      implicit none
      private	! except

      public :: mvcorF_nsep_setid
      public :: mvcorF_nsep_getid

      public :: mvcorF_nsep_alloc
      public :: mvcorF_nsep_dealloc

      public :: Hcoslim
      public :: qxHtb
      public ::   winHH
      public ::  DwinHH
      public :: DDwinHH
      public ::   winH
      public ::  DwinH
      public :: CLfact
      public :: rLav2

      public :: vfecHH
      public :: vfecHD
      public :: vfecDD
      public :: normHH
      public :: normDD

    interface mvcorF_nsep_setid; module procedure setid_; end interface
    interface mvcorF_nsep_getid; module procedure getid_; end interface
    interface mvcorF_nsep_alloc; module procedure alloc_; end interface
    interface mvcorF_nsep_dealloc; module procedure	&
	dealloc_; end interface

! !REVISION HISTORY:
!
! 	16Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- adapted from hfecHH.h and vfecHH.h for non-separable
!		  correlation implementation.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_mvcorF_nsep'
!_______________________________________________________________________
!
! HISTORY from "hfecHH.h":
!
!       21Feb2001- G. Gaspari  Redesigned from hfecHH.h
!       11Dec98 - J. Larson  Loop Blocking for Cache Optimization
!       27Nov98 - J. Larson  Incorporated Tom Clune's Cache 
!                            Optimizations:
!                              1) Fused Imat Tables hfecRR and 
!                                 hfecTT into hfecRRTT.
!                              2) Switched order of IMAT tau 
!                                 and level indices.
! 	21Nov95 - J. Guo
!		- added the prolog
!		- made adjustment of resolutions easier for direct
!		  table reference without interpolation.
!
!     Feb 93  - Jim Pfaendtner 	- original version 
!    2Feb 94  - Meta Sienkiewicz- declare coslim as parameter  
!   20Mar 94  - Meta Sienkiewicz- change coslim to variable (simplifies
!                                  expansion of tables if needed)
!   28Apr 94  - Meta Sienkiewicz- new tables based on cos(dist/rade) 
!                                  (cosine of separation angle)
!________________________________________
!
! HISTORY from "vfecHH.h":
!
!  18Sep95  - Jing G.	- Created multiple copies of statistical tables
!  22Feb96  - Jing G.	- Added vfecHD/vfecDD to normalize hfecHD/hfecDD
!
!  18Sep95  - Jing G.	- Created multiple copies of statistical tables
!			- Added the prolog
!  24Mar95  - Jing G.	- created for text based data table
!  02/20/93 -		- an old vfecH.h
!_______________________________________________________________________

			! Number of hfecH tables (see "MX_hfecH.h")

  include "MX_hfecH.h"

	! Book keeping data

  integer,save :: n_fecH=0
  integer,save,dimension(MX_fecH) :: kmat_tbl = -1

	! Data storage

  logical,save :: initialized_=.false.

  integer,save :: nHHtab  = 0		! number of ctau mesh points
  integer,save :: nveclev = 0		! number of level mesh points

					! upper mesh limit in tau
  real,save,allocatable,dimension(:    ) :: Hcoslim

					! mesh resolution in 1/dtau
  real,save,allocatable,dimension(:    ) :: qxHtb

			! G.G. 7/6/00, len scales & window fcn 

  real,save,allocatable,dimension(:    ) :: winHH
  real,save,allocatable,dimension(:    ) :: DwinHH
  real,save,allocatable,dimension(:    ) :: DDwinHH

  real,save,allocatable,dimension(:,:,:) :: CLfact
  real,save,allocatable,dimension(:,:,:) :: rLav2

		! winH and DwinH are used to set TT[itau=1] := 0

  real,save,allocatable,dimension(:) :: winH
  real,save,allocatable,dimension(:) :: DwinH

			!________________________________________
			! Vertical Forecast-Error correlation tables
			! for H related variables.

  real,save,allocatable,dimension(:,:,:) :: vfecHH	! H-H ver.cor.
  real,save,allocatable,dimension(:,:,:) :: vfecHD	! H-D ver.cor.
  real,save,allocatable,dimension(:,:,:) :: vfecDD	! D-D ver.cor.

				! norm_HH=1/sqrt(diag(fecHH))
				! norm_DD=1/sqrt(diag(fecDD))

  real,save,allocatable,dimension(  :,:) :: normHH
  real,save,allocatable,dimension(  :,:) :: normDD

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: setid_ - signing in an instance
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine setid_(kindcov,reset)
      use m_die, only : die
      implicit none
      integer,intent(in ) :: kindcov
      logical,optional,intent(in) :: reset

! !REVISION HISTORY:
! 	16Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::setid_'
  include "kind_covs.h"

  logical :: reset_

  reset_=.false.
  if(present(reset)) reset_=reset

  if(reset_) then
    n_fecH=0
    kmat_tbl(:)=-1
  endif

  n_fecH=n_fecH+1

  if(n_fecH>MX_fecH) call die(myname_,	&
	'n_fecH',n_fecH,'MX_fecH',MX_fecH)

  select case(kindcov)
  case(kind_covF)
    kmat_tbl(kmat_HGHT)=n_fecH

  case(kind_covS)
    kmat_tbl(kmat_STRM)=n_fecH

  case(kind_covV)
    kmat_tbl(kmat_VELP)=n_fecH

  case default
    call die(myname_,'unknown kind_cov',kindcov)
  end select

end subroutine setid_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getid_ - inquire previously _set_ ID.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getid_(kindcov,kmat)
      use m_die, only : die
      implicit none
      integer,intent(in ) :: kindcov
      integer,intent(out) :: kmat

! !REVISION HISTORY:
! 	16Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getid_'
  include "kind_covs.h"

  select case(kindcov)
  case(kind_covF)
    kmat=kmat_tbl(kmat_HGHT)

  case(kind_covS)
    kmat=kmat_tbl(kmat_STRM)

  case(kind_covV)
    kmat=kmat_tbl(kmat_VELP)

  case default
    call die(myname_,'unknown kind_cov',kindcov)
  end select

  if(kmat<=0 .or. kmat>n_fecH) call die(myname_,	&
	'kind_cov not registered',kindcov)
    
end subroutine getid_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: alloc_ - allocate the whole storage
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine alloc_(ntau,nlev)
      use m_mall,only : mall_ison,mall_mci
      use m_die ,only : die
      implicit none
      integer,intent(in) :: ntau	! number of tau mesh points
      integer,intent(in) :: nlev	! number of level mesh points

! !REVISION HISTORY:
! 	16Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::alloc_'
  integer :: ier

  if(initialized_) call die(myname_,'multiple module initialization')

  initialized_=.true.
  nHHtab =ntau
  nveclev=nlev

  allocate( Hcoslim(n_fecH),qxHtb(n_fecH),			&
	    winHH(nHHtab),DwinHH(nHHtab),DDwinHH(nHHtab),	&
	    winH (nHHtab),DwinH (nHHtab),			&
	    CLfact(nveclev,nveclev,n_fecH),			&
	    rLav2 (nveclev,nveclev,n_fecH), stat=ier)

	if(ier/=0) call die(myname_,'allocate(hfecHH)',ier)

	if(mall_ison()) then
	  call mall_mci(Hcoslim,myname)
	  call mall_mci(qxHtb  ,myname)
	  call mall_mci(  winHH,myname)
	  call mall_mci( DwinHH,myname)
	  call mall_mci(DDwinHH,myname)
	  call mall_mci(  winH ,myname)
	  call mall_mci( DwinH ,myname)
	  call mall_mci(CLfact ,myname)
	  call mall_mci(rLav2  ,myname)
	endif

  allocate( vfecHH(nveclev,nveclev,n_fecH),		&
	    vfecHD(nveclev,nveclev,n_fecH),		&
	    vfecDD(nveclev,nveclev,n_fecH),		&
	    normHH(        nveclev,n_fecH),		&
	    normDD(        nveclev,n_fecH), stat=ier)

	if(ier/=0) call die(myname_,'allocate(vfecHH)',ier)

	if(mall_ison()) then
	  call mall_mci(vfecHH,myname)
	  call mall_mci(vfecHD,myname)
	  call mall_mci(vfecDD,myname)
	  call mall_mci(normHH,myname)
	  call mall_mci(normDD,myname)
	endif

end subroutine alloc_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dealloc_ - deallocate the storage
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine dealloc_()
      use m_mall,only : mall_ison,mall_mco
      use m_die ,only : die
      implicit none

! !REVISION HISTORY:
! 	16Nov01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::dealloc_'
  integer :: ier

  if(.not.initialized_) call die(myname_,'uninitialized module')

  initialized_=.false.
  nHHtab =0
  nveclev=0

	if(mall_ison()) then
	  call mall_mco(Hcoslim,myname)
	  call mall_mco(qxHtb  ,myname)
	  call mall_mco(  winHH,myname)
	  call mall_mco( DwinHH,myname)
	  call mall_mco(DDwinHH,myname)
	  call mall_mco(  winH ,myname)
	  call mall_mco( DwinH ,myname)
	  call mall_mco(CLfact ,myname)
	  call mall_mco(rLav2  ,myname)
	endif

  deallocate(Hcoslim,qxHtb,winHH,DwinHH,DDwinHH,winH,DwinH,	&
	CLfact,rLav2,stat=ier)
	
	if(ier/=0) call die(myname_,'deallocate(hfecHH)',ier)

	if(mall_ison()) then
	  call mall_mco(vfecHH,myname)
	  call mall_mco(vfecHD,myname)
	  call mall_mco(vfecDD,myname)
	  call mall_mco(normHH,myname)
	  call mall_mco(normDD,myname)
	endif

  deallocate(vfecHH,vfecHD,vfecDD,normHH,normDD, stat=ier)

	if(ier/=0) call die(myname_,'deallocate(vfecHH)',ier)

end subroutine dealloc_

end module m_mvcorF_nsep
