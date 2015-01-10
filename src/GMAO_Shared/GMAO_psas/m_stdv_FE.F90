!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_stdv_FE - stdv. of fcst.err, given FE operators
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_stdv_FE
      implicit none
      private	! except

      public :: stdv_FE

      interface stdv_FE; module procedure	&
	stdv_,	&
	stdv_subB_, &
	stdvlev_
      end interface

! !REVISION HISTORY:
! 	04Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_stdv_FE'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: stdv_ - compute stdv. of fcst.err.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine stdv_(ndat,kts,rlat,rlon,rlev,sigF)
      use m_die,       only : die
      use m_stdio,     only : stderr
      use m_SortingTools,only : indexSet,indexSort
      use m_subBlocks, only : subBlocks,init,clean
      implicit none

		! Attributes of the output sigF

      integer, intent(in) :: ndat
      integer,dimension(ndat), intent(in) :: kts
      real,   dimension(ndat), intent(in) :: rlat
      real,   dimension(ndat), intent(in) :: rlon
      real,   dimension(ndat), intent(in) :: rlev

      real,   dimension(ndat), intent(out):: sigF

! !REVISION HISTORY:
! 	04Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::stdv_'
  integer :: ier
  integer :: i
  integer,allocatable,dimension(:) :: indx
  type(subBlocks) :: ktList

	! Initializing the _sorter_

  allocate(indx(ndat),stat=ier)
  if(ier /= 0) then
    write(stderr,'(2a,i4)') myname_,	&
	': allocate() error, stat =',ier
    call die(myname_)
  endif

  call indexSet(ndat,indx)
  call indexSort(ndat,indx,kts,descend=.false.)

	! Blocking.  Sorting is implied through the argument in a
	! array expression.

  call init(ktList,kts(indx(1:ndat)))

	! Computing.  Sorting is implied through the arguments in
	! array expressions.

  call stdv_subB_(ktList,	&
	rlat(indx(1:ndat)),	&
	rlon(indx(1:ndat)),	&
	rlev(indx(1:ndat)),	&
	sigF			)

	! Un-sorting

  sigF(indx(1:ndat))=sigF(1:ndat)

	! Cleaning up

  call clean(ktList)
  deallocate(indx)

end subroutine stdv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: stdv_subB_ - fcst.err.stdv. with a blocked vector.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine stdv_subB_(vAttr,rlat,rlon,rlev,sigF)
      use m_die,      only : die
      use m_stdio,    only : stderr
      use m_subBlocks,only : subBlocks,blockAttr,lsize
      use m_AttrVect, only : AttrVect
      implicit none

      type(subBlocks),  intent(in) :: vAttr
      real,dimension(:),intent(in) :: rlat,rlon,rlev
      real,dimension(:),intent(out):: sigF

! !REVISION HISTORY:
! 	04Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::stdv_subB_'

  integer :: kt,lc,ln,le
  integer :: l,ier

	! Compute the fcst.err.stdv. of wind variables

  do l=1,lsize(vAttr)

    call blockAttr(vAttr,l,kt=kt,lc=lc,ln=ln)
    
    if(ln <= 0) cycle

    le=lc+ln-1

    call stdv_var1_(kt,rlat(lc:le),rlon(lc:le),rlev(lc:le),sigF(lc:le))

  end do	! l=1,lsize(vAttr)

end subroutine stdv_subB_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: stdv_var1_ - fcst.err.stdv. with a blocked vector.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine stdv_var1_(kt,rlat,rlon,rlev,sigF)
      use m_die,    only : die
      use m_stdio,  only : stderr
      use config,   only : ktUU,ktUS,ktVV,ktVS

      use rlat_imat,only : nveclat,veclats
      use rlev_imat,only : nveclev,pveclev

      use m_xTab_lats,only : index_lats
      use m_xTab_levs,only : index_levs

      use FEalpha_imat,only : Aum_imat,Aul_imat,Avm_imat,Avl_imat
      use FEsigW_imat, only : FEsigS_imat,FEsigV_imat

      use m_xOp_sigFi, only : intp_sigHD => intp_sigFi

      implicit none
      integer,          intent(in) :: kt
      real,dimension(:),intent(in) :: rlat,rlon,rlev
      real,dimension(:),intent(out):: sigF

! !REVISION HISTORY:
! 	04Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::stdv_var1_'

  logical,parameter :: NEAREST=.true.

  integer :: ln
  real,   dimension(:),allocatable :: wlat,wlev
  integer,dimension(:),allocatable :: ilat,ilev

  integer :: i,ier
  real    :: am,al,bm,bl
  real    :: sigWs,sigWv,sigXs,sigXv
  integer :: lj,lv
  integer :: wj,wv

  ln=size(rlat)
  allocate( ilat(ln),wlat(ln), ilev(ln),wlev(ln), stat=ier)
	if(ier /= 0) then
	  write(stderr,'(2a,i4)') myname_,	&
		': allocate() error, stat =',ier
	  call die(myname_)
	endif

	! Compute the stdv. of missing-ratio and height-gradient
	! coupled fcst.err.

  call intp_sigHD(kt,rlat,rlon,rlev,sigF,stat=ier)

	if(ier /= 0) then
	  write(stderr,'(2a,i4)') myname_,	&
		': intp_sigHD() error, stat =',ier
	  call die(myname_)
	endif

	! Compute the fcst.err.stdv. of wind variables

    select case(kt)
    case(ktUU,ktUs,ktVV,ktVs)

      call index_levs(rlev,ilev,wlev)
      call index_lats(rlat,ilat,wlat)

      do i=1,ln
	lj=ilat(i)
	wj=wlat(i)
	lv=ilev(i)
	wv=wlev(i)

	select case(kt)
	case(ktUU,ktUs)

	  am=Aum_imat(lv,lj)
	  al=Aul_imat(lv,lj)
	  if(lj<nveclat) then
	    am=am+wj*(Aum_imat(lv,lj+1)-am)
	    al=al+wj*(Aul_imat(lv,lj+1)-al)
	  endif

	  if(lv<nveclev) then
	    bm=Aum_imat(lv+1,lj)
	    bl=Aul_imat(lv+1,lj)
	    if(lj<nveclat) then
	      bm=bm+wj*(Aum_imat(lv+1,lj+1)-bm)
	      bl=bl+wj*(Aul_imat(lv+1,lj+1)-bl)
	    endif

	    am=am+wv*(bm-am)
	    al=al+wv*(bl-al)
	  endif

	case(ktVV,ktVs)

	  am=Avm_imat(lv,lj)
	  al=Avl_imat(lv,lj)
	  if(lj<nveclat) then
	    am=am+wj*(Avm_imat(lv,lj+1)-am)
	    al=al+wj*(Avl_imat(lv,lj+1)-al)
	  endif

	  if(lv<nveclev) then
	    bm=Avm_imat(lv+1,lj)
	    bl=Avl_imat(lv+1,lj)
	    if(lj<nveclat) then
	      bm=bm+wj*(Avm_imat(lv+1,lj+1)-bm)
	      bl=bl+wj*(Avl_imat(lv+1,lj+1)-bl)
	    endif

	    am=am+wv*(bm-am)
	    al=al+wv*(bl-al)
	  endif

	end select

	sigWs=FEsigS_imat(lv,lj)	! already in wind units
	sigWv=FEsigV_imat(lv,lj)	! already in wind units
	if(lj.lt.nveclat) then
	  sigWs=sigWs + wj*(FEsigS_imat(lv,lj+1)-sigWs)
	  sigWv=sigWv + wj*(FEsigV_imat(lv,lj+1)-sigWv)
	endif
	if(lv.lt.nveclev) then
	  sigXs=FEsigS_imat(lv+1,lj)	! already in wind units
	  sigXv=FEsigV_imat(lv+1,lj)	! already in wind units
	  if(lj.lt.nveclat) then
	    sigXs=sigXs + wj*(FEsigS_imat(lv+1,lj+1)-sigXs)
            sigXv=sigXv + wj*(FEsigV_imat(lv+1,lj+1)-sigXv)
	  endif

	  sigWs=sigWs + wv*(sigXs-sigWs)
	  sigWv=sigWv + wv*(sigXv-sigWv)
	endif

	sigF(i)= sqrt( sigF(i)*sigF(i)*(am*am+al*al) + 	&
			  sigWs*sigWs + sigWv*sigWv )
      end do

    end select

  deallocate(ilat,wlat,ilev,wlev)
end subroutine stdv_var1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: stdvlev_ - interpolate with each level
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine stdvlev_(kt, mlon,mlat,pres,ilat,wlat,ilev,wlev,	&
    	sigF,kind_cov)

      use const,  only : G_EARTH  ! acc. due to the Earth gravity, m/ss
      use const,  only : R_EARTH  ! radius of Earth in meters
      use const,  only : OMEGA
      use const,  only : RHOBAR	  ! mean sea-level density
      use config, only : pres4slp

      use config,only : ktslp,ktHH,ktUU,ktVV,ktUs,ktVs,ktQQ

      use rlat_imat,only : nveclat
      use rlev_imat,only : nveclev

      use FEalpha_imat,only : Aum_imat,Aul_imat,Avm_imat,Avl_imat
      use FEsigW_imat, only : FEsigS_imat,FEsigV_imat

      use m_normCor,only : intp_normCor

      use m_sigFi_lookups
      use m_aGrid,only : LINEAR,LEVELS

      use m_intpAP, only : intpAP

      use m_stdio,only : stderr
      use m_die,  only : die

      implicit none
      integer,            intent(in)  :: kt
      integer,            intent(in)  :: mlon
      integer,            intent(in)  :: mlat
      real,               intent(in)  :: pres
      integer,dimension(:),intent(in) :: ilat
      real,   dimension(:),intent(in) :: wlat
      integer,            intent(in)  :: ilev
      real,               intent(in)  :: wlev
      real,dimension(:,:),intent(out) :: sigF
      integer,optional   ,intent(in)  :: kind_cov

      include "kind_covs.h"

! !REVISION HISTORY:
! 	13Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::stdvlev_'

  real, parameter :: scalar_p = .01*G_Earth*RHOBAR	! Pa to mb
  real, parameter :: scalar_d = G_EARTH/(2.*OMEGA*R_EARTH)

  integer :: i,j,ier
  integer :: lv,lj
  real    :: wv,wj
  real    :: anorm(1)
  real :: am,al,bm,bl
  real :: sigWs,sigWv,sigXs,sigXv
  real :: rlev
  integer :: kcov

  kcov=ior(kind_covF,ior(kind_covS,kind_covV))
  if(present(kind_cov)) kcov=kind_cov

  ier = -1

  if( .not.mixre_dim%perx .or. mixre_dim%mapx /=LINEAR	.or.	&
           mixre_dim%pery .or. mixre_dim%mapy /=LINEAR	.or.	&
                               mixre_dim%mapz /=LEVELS	.or.	&
      .not.hghte_dim%perx .or. hghte_dim%mapx /=LINEAR	.or.	&
           hghte_dim%pery .or. hghte_dim%mapy /=LINEAR	.or.	&
                               hghte_dim%mapz /=LEVELS	) then

    write(stderr,'(2a)') myname_,	&
	': can not handel mixre_dim%lon/lat.'

	write(stderr,*) mixre_dim%perx,mixre_dim%mapx,LINEAR
	write(stderr,*) mixre_dim%pery,mixre_dim%mapy,LINEAR
	write(stderr,*) .false.,       mixre_dim%mapz,LEVELS

	write(stderr,*) hghte_dim%perx,hghte_dim%mapx,LINEAR
	write(stderr,*) hghte_dim%pery,hghte_dim%mapy,LINEAR
	write(stderr,*) .false.,       hghte_dim%mapz,LEVELS

    call die(myname_)
  endif

  select case(kt)
  case(ktslp,ktHH,ktUU,ktVV,ktUs,ktVs)

    rlev=pres
    select case(kt)
    case(ktslp,ktUs,ktVs)
      rlev=pres4slp
    end select

    call intpAP(mlon,mlat,rlev,sigF,		&
	hghte_dim%nlon, hghte_dim%nlat,		&
	hghte_dim%nlev,	hghte_dim%vlev,		&
	hghte_fld,hghte_fill, stat=ier		)

    if(ier == 0) then
      select case(kt)
      case(ktslp)

	sigF(:,:) = scalar_p*sigF(:,:)

      case(ktUU,ktVV,ktUs,ktVs)

	call intp_normCor(1,(/pres/),anorm)
	sigF(:,:) = scalar_d*sigF(:,:)/anorm(1)

	lv=ilev
	wv=wlev

	do j=1,mlat
	  lj=ilat(j)
	  wj=wlat(j)

	  select case(kt)
	  case(ktUU,ktUs)

	    am=Aum_imat(lv,lj)
	    al=Aul_imat(lv,lj)
	    if(lj<nveclat) then
	      am=am+wj*(Aum_imat(lv,lj+1)-am)
	      al=al+wj*(Aul_imat(lv,lj+1)-al)
	    endif

	    if(lv<nveclev) then
	      bm=Aum_imat(lv+1,lj)
	      bl=Aul_imat(lv+1,lj)
	      if(lj<nveclat) then
	        bm=bm+wj*(Aum_imat(lv+1,lj+1)-bm)
	        bl=bl+wj*(Aul_imat(lv+1,lj+1)-bl)
	      endif

	      am=am+wv*(bm-am)
	      al=al+wv*(bl-al)
	    endif

	  case(ktVV,ktVs)

	    am=Avm_imat(lv,lj)
	    al=Avl_imat(lv,lj)
	    if(lj<nveclat) then
	      am=am+wj*(Avm_imat(lv,lj+1)-am)
	      al=al+wj*(Avl_imat(lv,lj+1)-al)
	    endif

	    if(lv<nveclev) then
	      bm=Avm_imat(lv+1,lj)
	      bl=Avl_imat(lv+1,lj)
	      if(lj<nveclat) then
	        bm=bm+wj*(Avm_imat(lv+1,lj+1)-bm)
	        bl=bl+wj*(Avl_imat(lv+1,lj+1)-bl)
	      endif

	      am=am+wv*(bm-am)
	      al=al+wv*(bl-al)
	    endif

	  end select

	  sigWs=FEsigS_imat(lv,lj)	! already in wind units
	  sigWv=FEsigV_imat(lv,lj)	! already in wind units
	  if(lj.lt.nveclat) then
	    sigWs=sigWs + wj*(FEsigS_imat(lv,lj+1)-sigWs)
	    sigWv=sigWv + wj*(FEsigV_imat(lv,lj+1)-sigWv)
	  endif
	  if(lv.lt.nveclev) then
	    sigXs=FEsigS_imat(lv+1,lj)	! already in wind units
	    sigXv=FEsigV_imat(lv+1,lj)	! already in wind units
	    if(lj.lt.nveclat) then
	      sigXs=sigXs + wj*(FEsigS_imat(lv+1,lj+1)-sigXs)
              sigXv=sigXv + wj*(FEsigV_imat(lv+1,lj+1)-sigXv)
	    endif

	    sigWs=sigWs + wv*(sigXs-sigWs)
	    sigWv=sigWv + wv*(sigXv-sigWv)
	  endif

	  if(iand(kcov,kind_covF)==0) sigF(:,j)=0.
	  if(iand(kcov,kind_covS)==0) sigWs=0.
	  if(iand(kcov,kind_covV)==0) sigWv=0.

	  do i=1,mlon
	    sigF(i,j)=sqrt( sigF(i,j)*sigF(i,j)*(am*am+al*al) + &
			    sigWs*sigWs + sigWv*sigWv )
	  end do
	end do	! j=1,mlat

      end select
    endif

  case(ktQQ)

    if(iand(kcov,kind_covF)==0) then
      sigF(:,:)=0.
    else
      call intpAP(mlon,mlat,pres,sigF,		&
	mixre_dim%nlon, mixre_dim%nlat,		&
	mixre_dim%nlev,	mixre_dim%vlev,		&
	mixre_fld,mixre_fill, stat=ier		)
    endif

  end select

  if(ier /= 0) then
    write(stderr,'(2a,i2.2,a,i4)') myname_,	&
	': intp_lev(',kt,') error, stat =',ier
    call die(myname_)
  endif

end subroutine stdvlev_

end module m_stdv_FE
!.
