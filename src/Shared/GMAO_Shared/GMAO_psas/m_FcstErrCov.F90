!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_FcstErrCov - A high level fcst.err.cov. operator
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_FcstErrCov
      implicit none
      private	! except

      public :: FcstErrCov_Cx
!      public :: FcstErrCov_Cxpy

      interface FcstErrCov_Cx  ; module procedure recCx_; end interface
!	symCx_,		&
!	recCx_; end interface

!      interface FcstErrCov_Cxpy; module procedure	&
!	symCxpy_,	&
!	recCxpy_; end interface

! !REVISION HISTORY:
! 	22Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_FcstErrCov'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recCx_ - a rectangular fcst.err.cov. matrix operator
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine recCx_(ninc,vinc,vlat,vlon,vlev,vkts,	&
		      nobs,xvec,rlat,rlon,rlev, kts,	&
	comm,root,rsrc	)

      use m_Navigator,only : Navigator

      use m_MultiAccessNavigator,only : MultiAccessNavigator
      use m_MultiAccessNavigator,only : clean
      use m_MultiAccessNavigator,only : ptr_typeNav
      use m_MultiAccessNavigator,only : ptr_krType
      use m_MultiAccessNavigator,only : ptr_ktType

      use m_GlobalPartition,only : GlobalPartition
      use m_GlobalPartition,only : GlobalPartition_init
      use m_GlobalPartition,only : clean

      use m_redwin,only : redwin_init
      use m_redwin,only : redwin_clean

      use m_sparse,only : sparse_init
      use m_sparse,only : sparse_clean

      use m_ob_Operators
      use m_ai_Operators

      use m_Distribution,only : Distribution
      use m_Distribution,only : Distribution_init
      use m_Distribution,only : clean
      use m_Distribution,only : distrSize

      use m_DistributionComm,only : distribute
      use m_DistributionComm,only : undistribute

      use m_Attributes,only : Attributes
      use m_Attributes,only : KR_SUBSET
      use m_Attributes,only : wrap
      use m_Attributes,only : clean
      use m_Attributes,only : ptr_kx
      use m_Attributes,only : ptr_ks
      use m_Attributes,only : ptr_kt
      use m_Attributes,only : ptr_lat
      use m_Attributes,only : ptr_lon
      use m_Attributes,only : ptr_lev

      use m_FcstErrCovMatx,only : FcstErrCovMatx
      use m_FcstErrCovMatx,only : FcstErrCovMatx_init
      use m_FcstErrCovMatx,only : FcstErrCovMatx_Cx
      use m_FcstErrCovMatx,only : FcstErrCovMatx_clean

      use m_sigmaPhi,only : sigmaPhi
      use m_sigmaPhi,only : sigmaPhi_init
      use m_sigmaPhi,only : col_Navigator
      use m_sigmaPhi,only : col_krNav
      use m_sigmaPhi,only : col_ktNav
      use m_sigmaPhi,only : col_rlat
      use m_sigmaPhi,only : col_rlon
      use m_sigmaPhi,only : col_rlev
      use m_sigmaPhi,only : clean

      use m_CorAttrF,only : CorAttrF
      use m_CorAttrF,only : CorAttrF_init
      use m_CorAttrF,only : clean

      use m_StatLevels,only : nStatLevel,StatLevels !,listStat

      use m_ErrCovModels,only : ErrCovModels_update
      use m_mvcorF_bldr,only : mvcorF_init
      use m_mvcorF_bldr,only : mvcorF_clean

      use m_xtenlev,  only : xtenlev_set

      use m_xTab_levs,only : tab_levs,levs_xtab,index_levs,rmtab_levs
      use m_xTab_lats,only : tab_lats,lats_xtab,index_lats,rmtab_lats
      use m_xTab_sigFi  ,only : rmtab_sigFi
      use m_FEsigFi_tabl,only : clean_sigFi => FEsigFi_clean
      use m_psasrc,   only : PSAS_DEFRC

      use m_costs ,only : costs_init, costs_clean, print_cost_statistics
      use m_parDOT,only : parDOT

      use m_die,  only : die
      use m_mall, only : mall_ison,mall_mci,mall_mco
      use m_zeit

      implicit none

      integer,intent(in) :: ninc
      real,dimension(:),intent(out) :: vinc	! vinc = (Cov)*xvec
      real,dimension(:),intent(in ) :: vlat	! latitudes
      real,dimension(:),intent(in ) :: vlon	! longitudes
      real,dimension(:),intent(in ) :: vlev	! levels
      integer,dimension(:),intent(in ) :: vkts	! types

      integer,intent(in) :: nobs
      real,dimension(:),intent(in ) :: xvec	! x
      real,target,dimension(:),intent(in ) :: rlat
      real,target,dimension(:),intent(in ) :: rlon
      real,target,dimension(:),intent(in ) :: rlev
      integer,target,dimension(:),intent(in ) :: kts
  
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
! 	22Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Cx_'

  include "kind_mats.h"

  integer :: nxkr
  integer :: ier
  integer :: i,l

  type(Navigator),pointer :: ob_Nav
  integer,pointer,dimension(:) :: ob_krNav
  integer,pointer,dimension(:) :: ob_ktNav

  type(Navigator),pointer :: ai_Nav
  integer,pointer,dimension(:) :: ai_krNav
  integer,pointer,dimension(:) :: ai_ktNav

  type(FcstErrCovMatx) :: fecov	! a fcst.err.cov. matrix operator

  type(sigmaPhi) :: ob_sigPhi
  type(sigmaPhi) :: ai_sigPhi

  type(CorAttrF) :: ob_corPhi
  type(CorAttrF) :: ai_corPhi

  type(Attributes) :: ob_attr0	! a temporary vector wrapper
  type(Attributes) :: ob_attr	! distributed vectors
  type(Distribution) :: ob_dstr
  type(MultiAccessNavigator) :: ob_man
  integer :: lobs
  integer,pointer,dimension(:) :: ob_kt
  real   ,pointer,dimension(:) :: ob_lat,ob_lon,ob_lev
  real,allocatable,dimension(:,:) :: ob_xvec

  type(Attributes) :: ai_attr0	! a temporary vector wrapper
  type(Attributes) :: ai_attr	! distributed vectors
  type(Distribution) :: ai_dstr
  type(MultiAccessNavigator) :: ai_man
  integer :: linc
  integer,pointer,dimension(:) :: ai_kt
  real   ,pointer,dimension(:) :: ai_lat,ai_lon,ai_lev
  real,allocatable,dimension(:,:) :: ai_vinc

  type(GlobalPartition) :: gp
!________________________________________

  if(present(rsrc)) then
    call GlobalPartition_init(gp,comm=comm,root=root,rsrc=rsrc)
  else
    call GlobalPartition_init(gp,comm=comm,root=root,rsrc=PSAS_DEFRC)
  endif
  call redwin_init(comm=comm,root=root,rsrc=rsrc)
  call sparse_init(gpart=gp,comm=comm,root=root)
!________________________________________

	call wrap(ob_attr0, nobs,rlat,rlon,rlev,kts,gp,comm)

		! Create a Distribution for observations (ob_dstr),
		! a local MultiAccssNavigator (ob_man), and
		! a copy of distributed Attributes (ob_attr).

  call distribution_init(ob_dstr,ob_attr,ob_man,	&
	(/KR_SUBSET/),ob_attr0,comm)

	call clean(ob_attr0)

		! This is the local size

  lobs = distrSize(ob_dstr)

	allocate( ob_xvec(lobs,1), stat=ier)
		if(ier/=0) call die(myname_,'allocate(ob_xvec)',ier)

		if(mall_ison()) call mall_mci(ob_xvec,myname)

  call distribute(xvec,ob_xvec(:,1),ob_dstr,comm)

  ob_Nav   => ptr_typeNav(ob_man)
  ob_krNav => ptr_krType(ob_man)
  ob_ktNav => ptr_ktType(ob_man)

  ob_kt  => ptr_kt(ob_attr)
  ob_lat => ptr_lat(ob_attr)
  ob_lon => ptr_lon(ob_attr)
  ob_lev => ptr_lev(ob_attr)

  call ob_Operators_init(lobs)
!________________________________________

	call wrap(ai_attr0, ninc,vlat,vlon,vlev,vkts,gp,comm)

		! Create a Distribution for analysis increments
		! (ai_dstr), a local MultiAccssNavigator (ai_man), and
		! a copy of distributed Attributes (ai_attr).

  call distribution_init(ai_dstr,ai_attr,ai_man,	&
	(/KR_SUBSET/),ai_attr0,comm)

	call clean(ai_attr0)

		! This is the local size

  linc = distrSize(ai_dstr)

	allocate(ai_vinc(linc,1),stat=ier)
		if(ier/=0) call die(myname_,'allocate(ai_vinc)',ier)

		if(mall_ison()) call mall_mci(ai_vinc,myname)

  ai_Nav   => ptr_typeNav(ai_man)
  ai_krNav => ptr_krType(ai_man)
  ai_ktNav => ptr_ktType(ai_man)

  ai_kt  => ptr_kt(ai_attr)
  ai_lat => ptr_lat(ai_attr)
  ai_lon => ptr_lon(ai_attr)
  ai_lev => ptr_lev(ai_attr)

  call ai_Operators_init(linc)
!________________________________________
  call ErrCovModels_update(root,comm,rsrc=rsrc)
!________________________________________

!  call ErrCovOperators_init() {

		! Merge in observation levels.  Notice that the level
		! for sea-level/surface observation or analysis is
		! assumed to be 1000 mb, and should be included in list
		! pres_list(:) already if the level would ever been
		! used in later computations.

	call xtenlev_set(linc,ai_xlev,ai_lev,ai_kt)
	call xtenlev_set(lobs,ob_xlev,ob_lev,ob_kt)

	call tab_levs()
	call levs_xtab(ai_xlev,comm)
	call levs_xtab(ob_xlev,comm)

	call tab_lats()
	call lats_xtab(ai_lat ,comm)
	call lats_xtab(ob_lat ,comm)

		! Set imat structures of correlation functions.  The
		! resultant imat structures are dependent on levels,
		! as well as distances (in cosine-radians).
!________________________________________
!  call FcstErrCov_init(..) {

	call mvcorF_init()	! was "call set_fecHH()"
	call set_fecQQ()

		! Set imat structures of alpha-operators (wind-mass
		! balancing scheme).  The resultant imat structures are
		! dependent on levels and latitudes.

	call imat_alpha()
	call imat_sigW()

!  }
!________________________________________

	call ll2qvec(lobs,ob_lat,ob_lon,	&
		ob_qrx,ob_qry,ob_qrz,		&
		ob_qmx,ob_qmy,ob_qmz,		&
		ob_qlx,ob_qly			)

		! Index rlev(:) to the common pressure level lookup
		! table.  Note that _nearest_ option is chosen.  Future
		! changes should make this index table operator
		! dependent.

        call index_levs(ob_xlev,ob_ktab)

		! Index rlat(:) to the common latitude lookup table.
		! Note that _nearest_ option is chosen.  Future changes
		! should make this index table operator dependent.

        call index_lats(ob_lat,ob_jtab)

	call sigmaPhi_init(ob_sigPhi,ob_Nav,ob_krNav,ob_ktNav,	&
		ob_lon,ob_lat,ob_xlev)
!________________________________________

	call CorAttrF_init(ob_corPhi,		&
		col_Navigator(ob_sigPhi),	&
		col_krNav(ob_sigPhi),		&
		col_ktNav(ob_sigPhi),		&
		col_rlon(ob_sigPhi),		&
		col_rlat(ob_sigPhi),		&
		col_rlev(ob_sigPhi),comm)
!________________________________________

	call ll2qvec(linc,ai_lat,ai_lon,	&
		ai_qrx,ai_qry,ai_qrz,		&
		ai_qmx,ai_qmy,ai_qmz,		&
		ai_qlx,ai_qly			)

		! Index rlev(:) to the common pressure level lookup
		! table.  Note that _nearest_ option is chosen.  Future
		! changes should make this index table operator
		! dependent.

        call index_levs(ai_xlev,ai_ktab)

		! Index rlat(:) to the common latitude lookup table.
		! Note that _nearest_ option is chosen.  Future changes
		! should make this index table operator dependent.

        call index_lats(ai_lat,ai_jtab)

	call sigmaPhi_init(ai_sigPhi,ai_Nav,ai_krNav,ai_ktNav,	&
		ai_lon,ai_lat,ai_xlev)
!________________________________________

	call CorAttrF_init(ai_corPhi,		&
		col_Navigator(ai_sigPhi),	&
		col_krNav(ai_sigPhi),		&
		col_ktNav(ai_sigPhi),		&
		col_rlon(ai_sigPhi),		&
		col_rlat(ai_sigPhi),		&
		col_rlev(ai_sigPhi),comm	)
!________________________________________

    call costs_init(comm=comm,root=root)
!_______________________________________________________________________

!	..Project the solution to grid points through [P^f][H]`[x] to
!	get Analysis Increments (AInc).
!	=============================================================

!	..Compute "normalized" AInc.
!	=============================

  call FcstErrCovMatx_init(fecov,ai_corPhi,ai_corPhi,	&
	ob_corPhi,ob_corPhi,comm,kind_5mat)

  call zeit_ci('FcstErrCovMatx_Cx')
  call FcstErrCovMatx_Cx(fecov,comm,			&
	ai_Nav,ai_krNav,ai_ktNav,			&
	  linc,ai_ktab,ai_jtab,ai_sigPhi,ai_corPhi,ai_corPhi,	&
	ob_Nav,ob_krNav,ob_ktNav,			&
	  lobs,ob_ktab,ob_jtab,ob_sigPhi,ob_corPhi,ob_corPhi,	&
	ob_xvec,ai_vinc)
  call zeit_co('FcstErrCovMatx_Cx')

  call FcstErrCovMatx_clean(fecov)
  
! Call print_cost_statistics()
  call costs_clean()
!________________________________________
		! "Undistribute" the outputs

  call undistribute(ai_vinc(1:linc,1),vinc,ai_dstr,comm)

!________________________________________

!  call ErrCovOperators_clean() {

	call clean(ob_sigPhi)
	call clean(ob_corPhi)
	call clean(ai_sigPhi)
	call clean(ai_corPhi)

	call ai_Operators_clean()
	call ob_Operators_clean()

	call rmtab_levs
	call rmtab_lats
	call mvcorF_clean()

!  }

!  call ErrCovModels_clean() {

	call rmtab_sigFi	! remove cached sigH/sigQ grids (file)
	call clean_sigFi	! remove sigFi filename

!  }
!______________________

	nullify(ob_Nav)
	nullify(ob_krNav)
	nullify(ob_ktNav)

	nullify(ob_kt)
	nullify(ob_lat)
	nullify(ob_lon)
	nullify(ob_lev)

	call clean(ob_dstr)
	call clean(ob_attr)
	call clean(ob_man)

		if(mall_ison()) call mall_mco(ob_xvec,myname)
	deallocate(ob_xvec,stat=ier)
		if(ier/=0) call die(myname_,'deallocate(ob_xvec)',ier)
!______________________

	nullify(ai_Nav)
	nullify(ai_krNav)
	nullify(ai_ktNav)

	nullify(ai_kt)
	nullify(ai_lat)
	nullify(ai_lon)
	nullify(ai_lev)

	call clean(ai_dstr)
	call clean(ai_attr)
	call clean(ai_man)

		if(mall_ison()) call mall_mco(ai_vinc,myname)
	deallocate(ai_vinc,stat=ier)
		if(ier/=0) call die(myname_,'deallocate(ai_vinc)',ier)

!______________________

	call sparse_clean()
	call redwin_clean()
	call clean(gp)
end subroutine recCx_
end module m_FcstErrCov
