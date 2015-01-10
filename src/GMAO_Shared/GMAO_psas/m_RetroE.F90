!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_RetroE - solves eq. for single retro-analysis 
!
! !DESCRIPTION:
!          Calculates:   xvec =>  M * x = H * Pf * y
!
! !INTERFACE:

    module m_RetroE
      implicit none
      private	! except

      public :: RetroE_solve	! Solve a PSAS analysis equation

      interface RetroE_solve; module procedure solve_; end interface

! !REVISION HISTORY:
!	31Aug01 - R. Todling - Adapted m_AE to present case.
!	12Jan01	- Jing Guo
!		. Updated the prolog and minor things for a prelease
!		  clean-up.
! 	22Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_RetroE'

! Note:
!
!  Name _solve() is for complicate equations, while _Cx() or _Cxpy()
!  are for solvers with simple meanings.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: solve_ - solves eq. for single retro-analysis
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine solve_( nobs,xvec,rlat,rlon,rlev,kxs,kss,kts,xmUs,&
	ninc,vinc,vlat,vlon,vlev,vkts,	&
	comm,root,rsrc	)

      use m_Navigator,only : Navigator

      use m_MultiAccessNavigator,only : MultiAccessNavigator
      use m_MultiAccessNavigator,only : clean
      use m_MultiAccessNavigator,only : ptr_typeNav
      use m_MultiAccessNavigator,only : ptr_krType
      use m_MultiAccessNavigator,only : ptr_ktType
      use m_MultiAccessNavigator,only : nType,getType

      use m_Regioner,only : Regioner_init
      use m_Regioner,only : Regioner_clean

      use m_ob_Operators
      use m_ai_Operators

      use m_Distribution,only : Distribution
      use m_Distribution,only : Distribution_init
      use m_Distribution,only : clean
      use m_Distribution,only : distrSize
      use m_Distribution,only : globalSize

      use m_ClusterDistribution,only : ClusterDistribution_init

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

      use m_CGSolver,only : CGSolver
      use m_CGSolver,only : CGSolver_init
      use m_CGSolver,only : CGSolver_solve
      use m_CGSolver,only : nband
      use m_CGSolver,only : CGSolver_clean

      use m_FcstErrCovMatx,only : FcstErrCovMatx
      use m_FcstErrCovMatx,only : FcstErrCovMatx_init
      use m_FcstErrCovMatx,only : FcstErrCovMatx_Cx
      use m_FcstErrCovMatx,only : FcstErrCovMatx_clean
      use m_FcstErrCovMatx,only : FcstErrCovMatx_showCosts

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

      use m_CorAttrX,only : CorAttrX
      use m_CorAttrX,only : CorAttrX_init
      use m_CorAttrX,only : clean

      use m_StatLevels,only : nStatLevel,StatLevels !,listStat

      use m_ErrCovModels,only : ErrCovModels_update
      use m_mvcorF_bldr,only : mvcorF_init
      use m_mvcorF_bldr,only : mvcorF_clean

      use m_xtenlev,  only : xtenlev_set

      use m_xTab_levs,only : tab_levs,levs_xtab,index_levs,rmtab_levs
      use m_xTab_lats,only : tab_lats,lats_xtab,index_lats,rmtab_lats
      use m_xTab_sigFi,only : rmtab_sigFi
      use m_FEsigFi_tabl,only : clean_sigFi => FEsigFi_clean

      use m_costs,only : costs_init,costs_clean,print_cost_statistics

      use m_parDOT,only : parDOT

      use m_showDistrib,only : showDistrib
      use m_die,  only : die
      use m_mall, only : mall_ison,mall_mci,mall_mco
      use m_mpout,only : mpout,mpout_ison,mpout_log

      use m_zeit

      implicit none

      integer,intent(in) :: ninc
      real,dimension(:),intent(in ) :: vinc	! perturbation vector
      real,dimension(:),intent(in ) :: vlat	! latitudes
      real,dimension(:),intent(in ) :: vlon	! longitudes
      real,dimension(:),intent(in ) :: vlev	! levels
      integer,dimension(:),intent(in ) :: vkts	! types

      integer,intent(in) :: nobs
      real,dimension(:),intent(out) :: xvec	! x (retro CG solution)
      real,target,dimension(:),intent(in ) :: rlat
      real,target,dimension(:),intent(in ) :: rlon
      real,target,dimension(:),intent(in ) :: rlev
      integer,target,dimension(:),intent(in ) :: kxs
      integer,target,dimension(:),intent(in ) :: kss
      integer,target,dimension(:),intent(in ) :: kts
      real,dimension(:),intent(in ) :: xmUs
  
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
!	31Aug01 - R. Todling - Adapted m_AE to present case.
!	12Jan01	- Jing Guo
!		. Updated the prolog and minor things for a prelease
!		  clean-up.
! 	22Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::solve_'

  include "bands.h"

  integer :: nxkr
  integer :: ier
  integer :: i,l

  type(Navigator),pointer :: ob_Nav
  integer,pointer,dimension(:) :: ob_krNav
  integer,pointer,dimension(:) :: ob_ktNav

  type(Navigator),pointer :: ai_Nav
  integer,pointer,dimension(:) :: ai_krNav
  integer,pointer,dimension(:) :: ai_ktNav

  type(CGSolver) :: pcg		! a nested preconditioned CG solver
  type(FcstErrCovMatx) :: fecov	! a fcst.err.cov. matrix operator

  type(sigmaPhi) :: ob_sigPhi
  type(sigmaPhi) :: ai_sigPhi

  type(CorAttrF) :: ob_corPhi
  type(CorAttrF) :: ai_corPhi
  type(CorAttrX) :: ob_corObs

  integer :: nbandMatx

  type(Attributes) :: ob_attr0		! a temporary wrapper
  type(Attributes) :: ob_attr		! distributed vectors
  type(Distribution) :: ob_dstr
  type(MultiAccessNavigator) :: ob_man
  integer :: lobs,mobs
  integer,pointer,dimension(:) :: ob_kx ,ob_kt ,ob_ks
  real   ,pointer,dimension(:) :: ob_lat,ob_lon,ob_lev
  real,allocatable,dimension(:) :: ob_xmUs
  real,allocatable,dimension(:,:) :: ob_yvec,ob_xvec

  type(Attributes) :: ai_attr0		! a temporary wrapper
  type(Attributes) :: ai_attr		! distributed vectors
  type(Distribution) :: ai_dstr
  type(MultiAccessNavigator) :: ai_man
  integer :: linc,minc
  integer,pointer,dimension(:) :: ai_kt
  real   ,pointer,dimension(:) :: ai_lat,ai_lon,ai_lev
  real,allocatable,dimension(:,:) :: ai_vinc

  integer :: kr,kt,kx,ln,itype
!________________________________________

  call Regioner_init(nxkr)
!________________________________________
	call zeit_ci('distribute_ob')

	call showDistrib(mpout,'pre-distrib-nobs',nobs,	&
		root,comm,listall=.true.)

	call wrap(ob_attr0, nobs,rlat,rlon,rlev,kxs,kss,kts)

		! Create a Distribution for observations (ob_dstr),
		! a local MultiAccssNavigator (ob_man), and
		! a copy of distributed Attributes (ob_attr).

  call ClusterDistribution_init(ob_dstr,ob_attr,ob_man,ob_attr0,root,comm,symmetric=.true.)

	call clean(ob_attr0)

		! Show the local and global sizes of the distributed
		! vector.

  lobs = distrSize(ob_dstr)
  mobs = globalSize(ob_dstr)

	call showDistrib(mpout,'post-distrib-nobs',lobs,	&
		root,comm,listall=.true.)

	call mpout_log(myname_,'globalSize(ob_dstr) =',mobs)
	call mpout_log(myname_,' localSize(ob_dstr) =',lobs)

	allocate( ob_yvec(lobs,1),ob_xvec(lobs,1),	&
		  ob_xmUs(lobs), stat=ier)
		if(ier/=0) call die(myname_,'allocate(ob_yvec)',ier)

		if(mall_ison()) then
		  call mall_mci(ob_yvec ,myname)
		  call mall_mci(ob_xvec,myname)
		  call mall_mci(ob_xmUs,myname)
		endif

  call distribute(xmUs,ob_xmUs    ,ob_dstr,comm)

  ob_Nav   => ptr_typeNav(ob_man)
  ob_krNav => ptr_krType(ob_man)
  ob_ktNav => ptr_ktType(ob_man)

  ob_kx  => ptr_kx(ob_attr)
  ob_ks  => ptr_ks(ob_attr)
  ob_kt  => ptr_kt(ob_attr)
  ob_lat => ptr_lat(ob_attr)
  ob_lon => ptr_lon(ob_attr)
  ob_lev => ptr_lev(ob_attr)

  call ob_Operators_init(lobs)

	call zeit_co('distribute_ob')
!________________________________________
	call zeit_ci('distribute_ai')

	call showDistrib(mpout,'pre-distrib-ninc',ninc,	&
		root,comm,listall=.true.)

	call wrap(ai_attr0, ninc,vlat,vlon,vlev,vkts)

		! Create a Distribution for analysis increments
		! (ai_dstr), a local MultiAccssNavigator (ai_man), and
		! a copy of distributed Attributes (ai_attr).

  call ClusterDistribution_init(ai_dstr,ai_attr,ai_man,ai_attr0,root,comm,symmetric=.false.)

	call clean(ai_attr0)

		! Show the local and global sizes of the distributed
		! vector.

  linc = distrSize(ai_dstr)
  minc = globalSize(ai_dstr)

	call showDistrib(mpout,'post-distrib-ninc',linc,	&
		root,comm,listall=.true.)

	call mpout_log(myname_,'globalSize(ai_dstr) =',minc)
	call mpout_log(myname_,' localSize(ai_dstr) =',linc)

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

  call distribute(vinc,ai_vinc(:,1),ai_dstr,comm)  ! initialize with input vector

  call ai_Operators_init(linc)

	call zeit_co('distribute_ai')
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

	! It looks awful.  However, this part code builds xtab_levs and
	! xtab_lats based on all data.  This communication may be
	! removed if the information can be determined locally by
	! 1) remove the dependency of sigW and alpha on IMATs; and
	! 2) determined IMATs for correlation functions with different
	!   algorithms.

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
!  call ObsErrCov_init(..) {

	call set_oecHH()
!-----------------------------------------------------------------------
	! Recreate sigO_list[] with sigU_list[] based the tables.  The
	! observation error deviations are created first, such that they
	! can be processed by restrict(), where a negative value is
	! considered for an invalid data point.

	call intp_sigO(lobs,ob_kx,ob_kt,ob_xlev,		&
		ob_sigO,ob_sigU				)

	if(mpout_ison()) then
	  call obstat(mpout, lobs,ob_kx,ob_kt,ob_xlev,ob_sigO,	&
	    nStatLevel,StatLevels(1:nStatLevel),myname//'*ObsErr*sigO')

	  call obstat(mpout, lobs,ob_kx,ob_kt,ob_xlev,ob_sigU,	&
	    nStatLevel,StatLevels(1:nStatLevel),myname//'*ObsErr*sigU')
	endif

	do i=1,lobs
	  if(ob_sigU(i)<0. .or. ob_sigO(i)<0.)	&
		call die(myname_,'incomplete ObsErrCov')
	end do

	do i=1,lobs
	  ob_sigU(i)=ob_sigU(i)*ob_xmUs(i)
	end do

	if(mpout_ison()) then
	  call obstat(mpout, lobs,ob_kx,ob_kt,ob_xlev,ob_sigU,	&
	    nStatLevel,StatLevels(1:nStatLevel),		&
	    myname//'*ObsErr*sigU(*xmUs)')
	endif

!  }
!________________________________________
!  call FcstErrCov_init(..) {

	call mvcorF_init()
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

	call CorAttrX_init(ob_corObs,		&
		ob_Nav,ob_krNav,ob_ktNav,	&
		ob_kx,ob_ks,			&
		ob_lon,ob_lat,ob_xlev,comm)
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

!  }

	call costs_init(comm=comm,root=root)

!_______________________________________________________________________

!	Initilize conjugate gradient to get nbandmatx
!	=============================================


  call CGSolver_init(pcg,ob_corObs,ob_corPhi,ob_corPhi,	&
	root,comm,rsrc=rsrc)

	! get the top level %iband of pcg for use in
	! (FcstErrCovMatx)fecov

  nbandmatx=nband(pcg)

!________________________________________
!	..Project the solution to grid points through [P^f][H]`[x] to
!	get Analysis Increments (AInc).
!	=============================================================

!	..Compute "normalized" AInc.
!	=============================

  call FcstErrCovMatx_init(fecov,ob_corPhi,ob_corPhi,	&
	ai_corPhi,ai_corPhi,comm,nbandmatx)

  call FcstErrCovMatx_showCosts(fecov,mpout,root,comm)

  call zeit_ci('FcstErrCovMatx_Cx')
  call FcstErrCovMatx_Cx(fecov,comm,			&
	ob_Nav,ob_krNav,ob_ktNav,			&
	  lobs,ob_ktab,ob_jtab,ob_sigPhi,ob_corPhi,ob_corPhi,	&
	ai_Nav,ai_krNav,ai_ktNav,			&
	  linc,ai_ktab,ai_jtab,ai_sigPhi,ai_corPhi,ai_corPhi,	&
	ai_vinc,ob_yvec)
  call zeit_co('FcstErrCovMatx_Cx')

  call FcstErrCovMatx_clean(fecov)


	if(mpout_ison()) then
	  call obstat(mpout, lobs,ob_kx,ob_kt,ob_xlev,ob_yvec(:,1),&
	    nStatLevel,StatLevels(1:nStatLevel),		&
	    myname//'*HP-vector')
	endif

!_______________________________________________________________________

!	Solve ([H][P^f][H]^t+[R])[x] = [r] for [x]
!	===============================================

	call zeit_ci('CGSolver_solve')

  call CGSolver_solve(pcg,comm,				&
	ob_Nav,ob_krNav,ob_ktNav,			&
	lobs,ob_sigU,ob_sigO,ob_corObs,			&
	ob_ktab,ob_jtab,ob_sigPhi,ob_corPhi,ob_corPhi,	&
	ob_yvec,ob_xvec,ier				)

	if(ier/=0) call die(myname_,'CGSolver_solve()',ier)

  call CGSolver_clean(pcg)

	call zeit_co('CGSolver_solve')

	if(mpout_ison()) then
	  call obstat(mpout, lobs,ob_kx,ob_kt,ob_xlev,ob_xvec(:,1), &
	    nStatLevel,StatLevels(1:nStatLevel),	&
	    myname//'*Solutions')
	endif

! call print_cost_statistics()
	call costs_clean()

!________________________________________
		! "Undistribute" the outputs

	call undistribute(ob_xvec(1:lobs,1),xvec,ob_dstr,comm)

!________________________________________

!  call ErrCovOperators_clean() {

	call clean(ob_corObs)

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

	nullify(ob_kx)
	nullify(ob_ks)
	nullify(ob_kt)
	nullify(ob_lat)
	nullify(ob_lon)
	nullify(ob_lev)

	call clean(ob_dstr)
	call clean(ob_attr)
	call clean(ob_man)

		if(mall_ison()) then
		  call mall_mco(ob_yvec,myname)
		  call mall_mco(ob_xvec,myname)
		  call mall_mco(ob_xmUs,myname)
		endif
	deallocate(ob_yvec,ob_xvec,ob_xmUs,stat=ier)
		if(ier/=0) call die(myname_,'deallocate(ob_yvec)',ier)
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

	call Regioner_clean()
end subroutine solve_
end module m_RetroE
