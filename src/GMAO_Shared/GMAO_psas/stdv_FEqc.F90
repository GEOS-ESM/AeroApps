!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: stdv_FEqc - define forecast err. stdv. at grid points
!
! !INTERFACE:

  subroutine stdv_FEqc( im,jnp,mlev,plev,	&
	slp_sigF,slu_sigF,slv_sigF,		&
	hght_sigF,uwnd_sigF,vwnd_sigF,mixr_sigF	)

    use config, only : ktmax
    use m_zeit, only : zeit_ci,zeit_co
    implicit none

	! Grid point information

  integer, intent(in) :: im	! no. of lon. grid points (circular)
  integer, intent(in) :: jnp	! no. of lat. grid points
  integer, intent(in) :: mlev	! no. of lev. grid points
  real,    intent(in) :: plev(mlev)	! pres. lev. values

	! surface variables at (2-d) grid points

  real,    intent(inout) :: slp_sigF(im,jnp) ! Fcst. err. stdv. of slp
  real,    intent(out)   :: slu_sigF(im,jnp) ! Fcst. err. stdv. of slu
  real,    intent(out)   :: slv_sigF(im,jnp) ! Fcst. err. stdv. of slv

	! upper-air variables at (3-d) grid points

  real,    intent(inout) :: hght_sigF(im,jnp,mlev) ! FE. stdv. of hght
  real,    intent(out)   :: uwnd_sigF(im,jnp,mlev) ! FE. stdv. of uwnd
  real,    intent(out)   :: vwnd_sigF(im,jnp,mlev) ! FE. stdv. of vwnd
  real,    intent(inout) :: mixr_sigF(im,jnp,mlev) ! FE. stdv. of mixr

! !DESCRIPTION:
!
!   stdv\_FEqc() returns forecast error standard deviations for a given
!   analysis grid point setting.
!
! !SYSTEM ROUTINES:
!
!   getenv(3f) is required to check if there is an user specified
!   resource file to use, through an environment variable PSASRC.
!   If the variable has not been set, a default filename (psas.rc) is
!   selected.
!
! !REVISION HISTORY:
!
! 	15Nov96 - J. Guo	- prototyped.  And confirmed the 
!		interface design with Genia Brin and David Lamich.
!	17Nov96 - J. Guo	- finished code and tested.
! 	01Oct98 - Jing Guo <guo@thunder> - revised the prolog for the
!				latest protex() rules.
!EOP ___________________________________________________________________

	! Local variables and parameters
character(len=*), parameter :: myname="stdv_FEqc"

	! ... code goes here

	call zeit_ci(myname)
  call init_()		! read the resources
  call derive_FEstdv_()	! derive gridded forecast error stdv.
  call clean_()		! clean the cached resources whenver possible
	call zeit_co(myname)

!=======================================================================
contains
!-------------------------------
subroutine init_()
  use m_die,     only : die
  use m_xRSRC_sigFi,only : fetch_sigFi

  use m_psasrc,  only : psasrc_open,psasrc_close
  implicit none

	! local variables

character(len=*),parameter :: myname_=myname//"::init_"
integer :: ier


		! Loading the resource handle

  call psasrc_open(stat=ier)
	if(ier/=0) call die(myname_,'psasrc_open()',ier)

  call fetch_sigFi()
  call set_FEhCor()	! FcstErr*hCor_xx/(hfecH_tbl,hfecQ_tbl)
  call set_FEvCor()	! FcstErr*vCor_xx/(vfecH_tbl,vfecQ_tbl)
  call tabl_FEsigW()	! FcstErr*Sigma_Wind/FEsigW_tabl
  call tabl_FEalpha()

		! release PSASRC

  call psasrc_close(stat=ier)
	if(ier/=0) call die(myname_,'psasrc_close()',ier)

end subroutine init_

subroutine clean_()
  ! use m_xRSRC_sigFi,only : clean_sigFi
	use m_FEsigFi_tabl, only : clean_sigFi => FEsigFi_clean
  implicit none

  character(len=*),parameter :: myname_=myname//"::clean_"
  call clean_sigFi()

end subroutine clean_

subroutine derive_FEstdv_()
  use m_xTab_sigFi,only : tab_sigFi,rmtab_sigFi
  use m_xTab_levs, only : tab_levs, rmtab_levs
  use m_xTab_lats, only : tab_lats, rmtab_lats
  use m_mvcorF_bldr,only : mvcorF_init
  use m_mvcorF_bldr,only : mvcorF_clean
  implicit none

character(len=*), parameter :: myname_=myname//"::derive_FEstdv_"


	! Set level/latitude table entries with built-in resolutions

  call tab_levs()
  call tab_lats()

	! Set imat structures of correlation functions.  The
	! resultant imat structures are dependent on levels,
	! as well as distances (in cosine-radians).

  call mvcorF_init()	! was "call set_fecHH()"

	! Set imat structures of alpha-operators (wind-mass
	! balancing scheme).  The resultant imat structures are
	! dependent on levels and latitudes.

  call imat_alpha()

	! Set imat structures of sig_Psi/sig_Chi operators (height-
	! decoupled winds).  The resultant imat structures are
	! dependent on levels and latitudes.

  call imat_sigW()

	! Load the sigFi file

  call tab_sigFi()

	! Derive gridded FcstErrStdv.

  call stdv_( im,jnp,mlev,plev,		&
	slp_sigF,slu_sigF,slv_sigF,		&
	hght_sigF,uwnd_sigF,vwnd_sigF,mixr_sigF	)

	! Clean up the lookup tables

  call rmtab_sigFi()
  call rmtab_lats()
  call rmtab_levs()
  call mvcorF_clean()

end subroutine derive_FEstdv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: stdv_ - interpolate to gridded variables
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine stdv_(mlon,mlat,mlev,plev,		&
	slp_sigF,slu_sigF,slv_sigF,		&
	hght_sigF,uwnd_sigF,vwnd_sigF,mixr_sigF	)

      use config, only : ktslp,ktHH,ktUU,ktVV,ktUs,ktVs,ktQQ
      use config, only : pres4slp

      use m_xTab_lats, only : index_lats
      use m_xTab_levs, only : index_levs

      use m_stdv_FE,   only : stdv_FE

      use m_die,  only : die
      use m_mall, only : mall_ison,mall_mci,mall_mco

      implicit none

      integer,intent(in) :: mlon
      integer,intent(in) :: mlat
      integer,intent(in) :: mlev
      real,dimension(:),intent(in) :: plev
      real,dimension(:,:),  target,intent(out) ::  slp_sigF
      real,dimension(:,:),  target,intent(out) ::  slu_sigF
      real,dimension(:,:),  target,intent(out) ::  slv_sigF
      real,dimension(:,:,:),target,intent(out) :: hght_sigF
      real,dimension(:,:,:),target,intent(out) :: uwnd_sigF
      real,dimension(:,:,:),target,intent(out) :: vwnd_sigF
      real,dimension(:,:,:),target,intent(out) :: mixr_sigF

! !REVISION HISTORY:
! 	13Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::stdv_'

!  logical,parameter :: NEAREST=.true.

  integer :: islp
  real    :: wslp
  integer,dimension(:),allocatable :: ilev
  real,   dimension(:),allocatable :: wlev
  integer,dimension(:),allocatable :: ilat
  real,   dimension(:),allocatable :: wlat
  real,   dimension(:),allocatable :: rlat

  integer :: j,l,ier
  real    :: dlat

  if(mlon==0.or.mlat==0.or.mlev==0) return

  allocate( rlat(mlat),	ilat(mlat),wlat(mlat),	&
			ilev(mlev),wlev(mlev),	&
			stat=ier		)

	if(ier /= 0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(rlat,myname)
	  call mall_mci(ilat,myname)
	  call mall_mci(wlat,myname)
	  call mall_mci(ilev,myname)
	  call mall_mci(wlev,myname)
	endif

		! Indices to the loopup tables of norm_DD[],
		! sigWs[], and sigWv[].

  dlat=180./max(1,mlat-1)
  do j=1,mlat
    rlat(j)=-90.+(j-1)*dlat
  end do

  call index_levs(plev,     ilev,wlev)
  call index_levs(pres4slp, islp,wslp)
  call index_lats(rlat,ilat,wlat)

		! Interpolate for the surface variables

  call stdv_FE(ktslp,mlon,mlat,pres4slp,	&
	ilat,wlat,islp,wslp,slp_sigF(:,:)	)

  call stdv_FE(ktUs, mlon,mlat,pres4slp,	&
	ilat,wlat,islp,wslp,slu_sigF(:,:)	)

  call stdv_FE(ktVs, mlon,mlat,pres4slp,	&
	ilat,wlat,islp,wslp,slv_sigF(:,:)	)

		! Interpolate for the upper-air variables
  do l=1,mlev

    call stdv_FE(ktHH,mlon,mlat,plev(l),	&
	ilat,wlat,ilev(l),wlev(l),		&
	hght_sigF(:,:,l)			)

    call stdv_FE(ktUU,mlon,mlat,plev(l),	&
	ilat,wlat,ilev(l),wlev(l),		&
	uwnd_sigF(:,:,l)			)

    call stdv_FE(ktVV,mlon,mlat,plev(l),	&
	ilat,wlat,ilev(l),wlev(l),		&
	vwnd_sigF(:,:,l)			)

    call stdv_FE(ktQQ,mlon,mlat,plev(l),	&
	ilat,wlat,ilev(l),wlev(l),		&
	mixr_sigF(:,:,l)			)

  end do

	if(mall_ison()) then
	  call mall_mco(rlat,myname)
	  call mall_mco(ilat,myname)
	  call mall_mco(wlat,myname)
	  call mall_mco(ilev,myname)
	  call mall_mco(wlev,myname)
	endif
  deallocate(rlat,ilat,wlat,ilev,wlev)
end subroutine stdv_
!----------------------------------------
end subroutine stdv_FEqc
!.
