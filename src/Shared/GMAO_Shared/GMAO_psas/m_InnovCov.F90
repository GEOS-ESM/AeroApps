!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_InnovCov - an innovation covariance operator
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_InnovCov
      implicit none
      private	! except

      public :: InnovCov_Cx	! Solve an innov. cov. equation

      interface InnovCov_cx; module procedure symCx_; end interface

! !REVISION HISTORY:
!	12Jan01	- Jing Guo
!		- Converted from m_AE.F90
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_InnovCov'

! Note:
!
!  Name _solve() is for complicate equations, while _Cx() or _Cxpy()
!  for solvers with simple meanings.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: symCx_ - solve a PSAS analysis equation
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symCx_(nobs,xvec,rlat,rlon,rlev,	&
	kxs,kss,kts,xmUs,yvec,comm,root,rsrc	)

      use m_Navigator,only : Navigator

      use m_MultiAccessNavigator,only : MultiAccessNavigator
      use m_MultiAccessNavigator,only : clean
      use m_MultiAccessNavigator,only : ptr_typeNav
      use m_MultiAccessNavigator,only : ptr_krType
      use m_MultiAccessNavigator,only : ptr_ktType

      use m_GlobalPartition,only : GlobalPartition
      use m_GlobalPartition,only : GlobalPartition_init
      use m_GlobalPartition,only : clean

      use m_amat    ,only : amat_init
      use m_amat    ,only : amat_clean
      use m_amat    ,only : MAXBAND

      use m_sparse  ,only : sparse_init
      use m_sparse  ,only : sparse_clean

      use m_ob_Operators

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

      use m_InnovCovMatx,only : InnovCovMatx
      use m_InnovCovMatx,only : InnovCovMatx_init
      use m_InnovCovMatx,only : InnovCovMatx_Cx
      use m_InnovCovMatx,only : InnovCovMatx_clean

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
      use m_mvcorF_bldr,only : mvcorF_init
      use m_mvcorF_bldr,only : mvcorF_clean

      use m_ErrCovModels,only : ErrCovModels_update

      use m_xtenlev,  only : xtenlev_set

      use m_xTab_levs,only : tab_levs,levs_xtab,index_levs,rmtab_levs
      use m_xTab_lats,only : tab_lats,lats_xtab,index_lats,rmtab_lats
      use m_xTab_sigFi,only : rmtab_sigFi
      use m_FEsigFi_tabl,only : clean_sigFi => FEsigFi_clean

      use m_die   , only : die
      use m_mall  , only : mall_ison,mall_mci,mall_mco
      use m_mpout , only : mpout,mpout_ison
      use m_zeit  , only : zeit_ci,zeit_co
      use m_psasrc, only : psasrc_open,psasrc_close

      implicit none

      integer,intent(in) :: nobs
      real,dimension(:),intent(in ) :: xvec	! x
      real,target,dimension(:),intent(in ) :: rlat
      real,target,dimension(:),intent(in ) :: rlon
      real,target,dimension(:),intent(in ) :: rlev
      integer,target,dimension(:),intent(in ) :: kxs
      integer,target,dimension(:),intent(in ) :: kss
      integer,target,dimension(:),intent(in ) :: kts
      real,dimension(:),intent(in ) :: xmUs
      real,dimension(:),intent(out) :: yvec	! y=(IM)x
  
      integer,intent(in) :: root
      integer,intent(in) :: comm
      character(len=*),optional,intent(in) :: rsrc

! !REVISION HISTORY:
!	12Jan01	- Jing Guo
!		. Updated the prolog and minor things for a prelease
!		  clean-up.
! 	22Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::symCx_'

  integer :: ier
  integer :: i,l

  type(Navigator),pointer :: ob_Nav
  integer,pointer,dimension(:) :: ob_krNav
  integer,pointer,dimension(:) :: ob_ktNav

  type(InnovCovMatx) :: incov	! an innovation cov. matrix operator

  type(sigmaPhi) :: ob_sigPhi

  type(CorAttrF) :: ob_corPhi
  type(CorAttrX) :: ob_corObs

  type(Attributes) :: ob_attr0		! a temporary wrapper
  type(Attributes) :: ob_attr		! distributed vectors
  type(Distribution) :: ob_dstr
  type(MultiAccessNavigator) :: ob_man
  integer :: lobs
  integer,pointer,dimension(:) :: ob_kx ,ob_kt ,ob_ks
  real   ,pointer,dimension(:) :: ob_lat,ob_lon,ob_lev
  real,allocatable,dimension(:) :: ob_xmUs
  real,allocatable,dimension(:,:) :: ob_yvec,ob_xvec

  type(GlobalPartition) :: gp
!________________________________________

  call GlobalPartition_init(gp,comm=comm,root=root,rsrc=rsrc)
  call amat_init(gpart=gp,comm=comm,root=root,rsrc=rsrc)
  call sparse_init(gpart=gp,comm=comm,root=root)
!________________________________________

	call wrap(ob_attr0, nobs,rlat,rlon,rlev,kxs,kss,kts,gp,comm)

		! Create a Distribution for observations (ob_dstr),
		! a local MultiAccssNavigator (ob_man), and
		! a copy of distributed Attributes (ob_attr).

  call distribution_init(ob_dstr,ob_attr,ob_man,	&
	(/KR_SUBSET/),ob_attr0,comm)

	call clean(ob_attr0)

		! This is the local size

  lobs = distrSize(ob_dstr)

	allocate( ob_yvec(lobs,1),ob_xvec(lobs,1),	&
		  ob_xmUs(lobs), stat=ier)
		if(ier/=0) call die(myname_,'allocate(ob_yvec)',ier)

		if(mall_ison()) then
		  call mall_mci(ob_yvec,myname)
		  call mall_mci(ob_xvec,myname)
		  call mall_mci(ob_xmUs,myname)
		endif

  call distribute(xmUs,ob_xmUs     ,ob_dstr,comm)
  call distribute(xvec,ob_xvec(:,1),ob_dstr,comm)

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
!________________________________________

  call ErrCovModels_update(root,comm,rsrc=rsrc)
!________________________________________

!  call ErrCovOperators_init() {

		! Merge in observation levels.  Notice that the level
		! for sea-level/surface observation or analysis is
		! assumed to be 1000 mb, and should be included in list
		! pres_list(:) already if the level would ever been
		! used in later computations.

	call xtenlev_set(lobs,ob_xlev,ob_lev,ob_kt)

	! It looks awful.  However, this part code builds xtab_levs and
	! xtab_lats based on all data.  This communication may be
	! removed if the information can be determined locally by
	! 1) remove the dependency of sigW and alpha on IMATs; and
	! 2) determined IMATs for correlation functions with different
	!   algorithms.

	call tab_levs()
	call levs_xtab(ob_xlev,comm)

	call tab_lats()
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
	  nStatLevel,StatLevels(1:nStatLevel), myname//'*ObsErr*sigO')

	call obstat(mpout, lobs,ob_kx,ob_kt,ob_xlev,ob_sigU,	&
	  nStatLevel,StatLevels(1:nStatLevel), myname//'*ObsErr*sigU')
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
	  nStatLevel,StatLevels(1:nStatLevel),			&
	  myname//'*ObsErr*sigU(*xmUs)')

	call obstat(mpout, lobs,ob_kx,ob_kt,ob_xlev,ob_xvec(:,1),&
	  nStatLevel,StatLevels(1:nStatLevel), myname//'*xvec')
      endif

!  }
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

	call CorAttrX_init(ob_corObs,		&
		ob_Nav,ob_krNav,ob_ktNav,	&
		ob_kx,ob_ks,			&
		ob_lon,ob_lat,ob_xlev,comm)
!________________________________________
!  }

!_______________________________________________________________________

!	Solve y=([H][P^f][H]^t+[R]) x
!	=============================

  call InnovCovMatx_init(incov,ob_corObs,ob_corPhi,ob_corPhi,	&
	comm,MAXBAND)

	call zeit_ci('InnovCovMatx_Cx')

  call InnovCovMatx_Cx(incov,comm,ob_Nav,ob_krNav,ob_ktNav,	&
	lobs,ob_sigU,ob_sigO,ob_corObs,			&
	ob_ktab,ob_jtab,ob_sigPhi,ob_corPhi,ob_corPhi,	&
	ob_xvec,ob_yvec)
	
	call zeit_co('InnovCovMatx_Cx')

  call InnovCovMatx_clean(incov)

      if(mpout_ison()) then
	call obstat(mpout, lobs,ob_kx,ob_kt,ob_xlev,ob_yvec(:,1), &
	  nStatLevel,StatLevels(1:nStatLevel), myname//'*yvec')
      endif
!________________________________________

		! "Undistribute" the outputs

	call undistribute(ob_yvec(1:lobs,1),yvec,ob_dstr,comm)
!________________________________________

!  call ErrCovOperators_clean() {

	call clean(ob_corObs)

	call clean(ob_sigPhi)
	call clean(ob_corPhi)

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

	call sparse_clean()
	call amat_clean()
	call clean(gp)
end subroutine symCx_
end module m_InnovCov
