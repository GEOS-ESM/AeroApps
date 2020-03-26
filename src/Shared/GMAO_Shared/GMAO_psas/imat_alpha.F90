!#define	_DEBUG	!
subroutine imat_alpha()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: imat_alpha - (to do)
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	01Mar96 - J. Guo	- (to do)
!_______________________________________________________________________
use rlat_imat
use rlev_imat
use FEalpha_tabl
use FEalpha_imat
use m_chars, only: uppercase
use m_stdio,only : stderr
use m_die,  only : die
implicit none

	! parameters

  character(len=*), parameter	:: myname='imat_alpha'
  logical,	    parameter	:: nearest=.true.	! a flag

	! Local vars.
  integer :: i,j,ln,istat
!-----------------------------------------------------------------------
 
  select case (uppercase(FEalpha_type))
  case ('PSAS:X1')

    if(FEalpha_nlev.le.0.or.FEalpha_npar.lt.5) then
      ln=max(len_trim(FEalpha_type),1)

      if(FEalpha_nlev.le.0) write(stderr,'(4a,i3)') myname,	&
	  ': invalid table for "',FEalpha_type(1:ln),		&
	  '", FEalpha_nlev =',FEalpha_nlev

      if(FEalpha_npar.lt.5) write(stderr,'(4a,i3)') myname,	&
	  ': invalid table for "',FEalpha_type(1:ln),		&
	  '", expecting 5 parameters, but only ',FEalpha_npar

      call die(myname)
    endif

    call EXP_ALPHAS_model_x1

  case ('PSAS:X2')

    if(FEalpha_nlev.le.0.or.FEalpha_npar.lt.5) then
      ln=max(len_trim(FEalpha_type),1)

      if(FEalpha_nlev.le.0) write(stderr,'(4a,i3)') myname,	&
	  ': invalid table for "',FEalpha_type(1:ln),		&
	  '", FEalpha_nlev =',FEalpha_nlev

      if(FEalpha_npar.lt.5) write(stderr,'(4a,i3)') myname,	&
	  ': invalid table for "',FEalpha_type(1:ln),		&
	  '", expecting 5 parameters, but only ',FEalpha_npar

      call die(myname)
    endif

    call EXP_ALPHAS_model_x2

  case ('EXP-ALPHAS','PSAS:X0')

	! A expotentially adjusted ageostrophically balanced model

    if(FEalpha_nlev.le.0.or.FEalpha_npar.lt.3) then
      ln=max(len_trim(FEalpha_type),1)

      if(FEalpha_nlev.le.0) write(stderr,'(4a,i3)') myname,	&
	  ': invalid table for "',FEalpha_type(1:ln),		&
	  '", FEalpha_nlev =',FEalpha_nlev

      if(FEalpha_npar.lt.3) write(stderr,'(4a,i3)') myname,	&
	  ': invalid table for "',FEalpha_type(1:ln),		&
	  '", expecting 3 parameters, but only ',FEalpha_npar

      call die(myname)
    endif

    call EXP_ALPHAS_model_

  case ('PSAS:X3')

	! A expotentially adjusted ageostrophically balanced model
        ! A slight modification of X0, allowing a latitudinal
        ! independent frictionalo component (for use near the sfc)

    if(FEalpha_nlev.le.0.or.FEalpha_npar.lt.4) then
      ln=max(len_trim(FEalpha_type),1)

      if(FEalpha_nlev.le.0) write(stderr,'(4a,i3)') myname,	&
	  ': invalid table for "',FEalpha_type(1:ln),		&
	  '", FEalpha_nlev =',FEalpha_nlev

      if(FEalpha_npar.lt.4) write(stderr,'(4a,i3)') myname,	&
	  ': invalid table for "',FEalpha_type(1:ln),		&
	  '", expecting 3 parameters, but only ',FEalpha_npar

      call die(myname)

    endif

    call EXP_ALPHAS_model_x3

  case ('EPSILONS')

	! A four-parameter model

    if(FEalpha_nlev.le.0.or.FEalpha_npar.lt.4) then
      ln=max(len_trim(FEalpha_type),1)

      if(FEalpha_nlev.le.0) write(stderr,'(4a,i3)') myname,	&
	  ': invalid table for "',FEalpha_type(1:ln),		&
	  '", FEalpha_nlev =',FEalpha_nlev

      if(FEalpha_npar.lt.4) write(stderr,'(4a,i3)') myname,	&
	  ': invalid table for "',FEalpha_type(1:ln),		&
	  '", expecting 4 parameters, but only ',FEalpha_npar

      call die(myname)
    endif

    call EPSILONS_model_

  case ('GEOS/DAS-OI')

	! A model similar to the GEOS/DAS-OI mass-wind balance model
	! which expecting no parameter

    call OI_model_

  case default

	! Others are considered an error

    ln=max(len_trim(FEalpha_type),1)
    write(stderr,'(4a)') myname,': unknown model, "',	&
	FEalpha_type(1:ln),'"'
    call die(myname)
  end select
!=======================================================================

#ifdef	_DEBUG
	      	write(*,*) myname,'::'
		write(*,*) nveclev,nveclat
		write(*,*) 'Aum, Aul, Avm, Avl'
		do i=1,nveclev,(nveclev-1)/2
		  write(*,*) 'level = ',i,pveclev(i)
		  do j=1,nveclat
		    write(*,'(f7.2,4f8.3)') veclats(j),		&
			Aum_imat(i,j),Aul_imat(i,j),		&
			Avm_imat(i,j),Avl_imat(i,j)
		  end do
		end do
#endif

contains
!-----------------------------------------------------------------------
  subroutine OI_model_

	! OI_model tries to make a matrix similar to the effect by the
	! OI mass-wind balanced covariance model.  It would not be the
	! same, but may be close.

    use const, only	: OMEGA, ONEPI
    implicit none

    integer j,lv
    real Aum_upa,Aul_upa,Avm_upa,Avl_upa
    real Aum_sfc,Aul_sfc,Avm_sfc,Avl_sfc
    real wlev

    real Alf,rlat_r,wt,Hfac

    real, parameter	:: RAD_PER_DEG	= ONEPI/180.
    real, parameter	:: WDRG		= 1.e-5/(2.*OMEGA)

    real, parameter	:: reflat = 25.
    real, parameter	:: dellat =  5.

!-----------------------------------------------------------------------
	! OI has a special treatment for the location within 10 degrees
	! around the eauator.  This special handling is not implemented
	! because with the PSAS implementation approach, it would result
	! three-dimentional (Aum, Aul, Avm, Avl) IMATs.

    do j=1,nveclat

      rlat_r=veclats(j)*RAD_PER_DEG

      if(abs(veclats(j)).lt.reflat-dellat) then

	Alf=sign(sin((reflat-dellat)*RAD_PER_DEG),rlat_r)

      elseif(abs(veclats(j)).lt.reflat+dellat) then

	wt=.5*(1.+sin(ONEPI*(abs(veclats(j))-reflat)/(2.*dellat)))
	Alf=(1.-wt)*sign(sin((reflat-dellat)*RAD_PER_DEG),rlat_r) + &
		wt*sin(rlat_r)

      else

	Alf=sin(rlat_r)

      endif

      Hfac=1.-exp(-abs(veclats(j))/reflat)

      Aum_upa=-Hfac/Alf
      Aul_upa= 0.
      Avm_upa= 0.
      Avl_upa= Hfac/Alf


      Aum_sfc =- Alf*Hfac/(Alf*Alf+WDRG*WDRG)
      Aul_sfc =-WDRG*Hfac/(Alf*Alf+WDRG*WDRG)
      Avm_sfc =-WDRG*Hfac/(Alf*Alf+WDRG*WDRG)
      Avl_sfc =  Alf*Hfac/(Alf*Alf+WDRG*WDRG)

      do lv=1,nveclev
	if(pveclev(lv).le.850.) then
          Aum_imat(lv,j)=Aum_upa
	  Aul_imat(lv,j)=Aul_upa
	  Avm_imat(lv,j)=Avm_upa
	  Avl_imat(lv,j)=Avl_upa
	elseif(pveclev(lv).ge.1000.) then
	  Aum_imat(lv,j)=Aum_sfc
	  Aul_imat(lv,j)=Aul_sfc
	  Avm_imat(lv,j)=Avm_sfc
	  Avl_imat(lv,j)=Avl_sfc
	else
	  wlev=log(pveclev(lv)/1000.)/log(850./1000.)
	  Aum_imat(lv,j)=Aum_sfc + wlev*(Aum_upa-Aum_sfc)
	  Aul_imat(lv,j)=Aul_sfc + wlev*(Aul_upa-Aul_sfc)
	  Avm_imat(lv,j)=Avm_sfc + wlev*(Avm_upa-Avm_sfc)
	  Avl_imat(lv,j)=Avl_sfc + wlev*(Avl_upa-Avl_sfc)
        endif
      end do

    end do

  end subroutine OI_model_
!-----------------------------------------------------------------------
  subroutine EPSILONS_model_

	! This model assume a latitude independent epsilon value that
	! varies with respect to pressure levels only.

    use const, only	: OMEGA, ONEPI
    use m_die, only	: die
    use m_mall,only	: mall_ison,mall_mci,mall_mco
    implicit none

    integer j,k,il,istat
    real wt,slat
    real Exuk,Exvk,Eyuk,Eyvk,Edet
    real,    allocatable :: Exu(:),Exv(:),Eyu(:),Eyv(:)
    integer, allocatable :: ilev(:)
    real,    allocatable :: wlev(:)

    real, parameter	:: RAD_PER_DEG	= ONEPI/180.

	! index and assign log-linear weights to the table

    allocate(ilev(1:nveclev),wlev(1:nveclev),		&
	      Exu(1:nveclev), Exv(1:nveclev),		&
	      Eyu(1:nveclev), Eyv(1:nveclev), stat=istat)
    if(istat.ne.0)	&
      call die(myname,'EPSILONS_model_ allocate()',istat)

	if(mall_ison()) then
	  call mall_mci(ilev,myname)
	  call mall_mci(wlev,myname)
	  call mall_mci( Exu,myname)
	  call mall_mci( Exv,myname)
	  call mall_mci( Eyu,myname)
	  call mall_mci( Eyv,myname)
	endif

    call slogtab(.not.nearest, FEalpha_nlev,FEalpha_plev,	&
		nveclev,pveclev, ilev,wlev		)

    do k=1,nveclev
      il=ilev(k)		! conditioned log-linear interpolation
      wt=wlev(k)

      Exuk=FEalpha_pars(1,il)
      Eyuk=FEalpha_pars(2,il)
      Exvk=FEalpha_pars(3,il)
      Eyvk=FEalpha_pars(4,il)
      if(il.lt.FEalpha_nlev) then
	Exuk=Exuk+wt*(FEalpha_pars(1,il+1)-Exuk)
	Eyuk=Eyuk+wt*(FEalpha_pars(2,il+1)-Eyuk)
	Exvk=Exvk+wt*(FEalpha_pars(3,il+1)-Exvk)
	Eyvk=Eyvk+wt*(FEalpha_pars(4,il+1)-Eyvk)
      endif

      Exu(k)=Exuk
      Eyu(k)=Eyuk
      Exv(k)=Exvk
      Eyv(k)=Eyvk
    end do

    do j=1,nveclat
      slat=sin(veclats(j)*RAD_PER_DEG)

      do k=1,nveclev

	Edet=Exu(k)*Eyv(k)+Exv(k)*Eyu(k)*slat*slat

	Aum_imat(k,j) = -slat*Exv(k)/Edet
	Aul_imat(k,j) =      -Eyv(k)/Edet
	Avm_imat(k,j) =      -Exu(k)/Edet
	Avl_imat(k,j) =  slat*Eyu(k)/Edet

      end do
    end do

	if(mall_ison()) then
	  call mall_mco(ilev,myname)
	  call mall_mco(wlev,myname)
	  call mall_mco( Exu,myname)
	  call mall_mco( Exv,myname)
	  call mall_mco( Eyu,myname)
	  call mall_mco( Eyv,myname)
	endif

    deallocate(ilev,wlev,Exu,Exv,Eyu,Eyv)
  end subroutine EPSILONs_model_
!-----------------------------------------------------------------------
  subroutine EXP_ALPHAS_model_

	! This model assume a latitude dependent epsilon value that
	! varies with respect to pressure levels too.

    use const, only	: OMEGA, ONEPI
    use m_die, only 	: die
    use m_mall,only	: mall_ison,mall_mci,mall_mco
    implicit none

    integer j,k,il,istat
    real wt,slat,rlat
    real Ampk,scLk,Betk,Edet,eps
    real,    allocatable :: Amp(:),scL(:),Bet(:)
    integer, allocatable :: ilev(:)
    real,    allocatable :: wlev(:)

    real, parameter	:: RAD_PER_DEG	= ONEPI/180.

	! index and assign log-linear weights to the table

    allocate(ilev(1:nveclev),wlev(1:nveclev),		&
	      Amp(1:nveclev), scL(1:nveclev),		&
	      Bet(1:nveclev), stat=istat		)
    if(istat.ne.0)	&
      call die(myname,'EPSILONS_model_ allocate()',istat)

	if(mall_ison()) then
	  call mall_mci(ilev,myname)
	  call mall_mci(wlev,myname)
	  call mall_mci( Amp,myname)
	  call mall_mci( scL,myname)
	  call mall_mci( Bet,myname)
	endif

    call slogtab(.not.nearest, FEalpha_nlev,FEalpha_plev,	&
		nveclev,pveclev, ilev,wlev		)

    do k=1,nveclev
      il=ilev(k)		! conditioned log-linear interpolation
      wt=wlev(k)

      Ampk=FEalpha_pars(1,il)
      scLk=FEalpha_pars(2,il)
      Betk=FEalpha_pars(3,il)
      if(il.lt.FEalpha_nlev) then
	Ampk=Ampk+wt*(FEalpha_pars(1,il+1)-Ampk)
	scLk=scLk+wt*(FEalpha_pars(2,il+1)-scLk)
	Betk=Betk+wt*(FEalpha_pars(3,il+1)-Betk)
      endif

      Amp(k)=Ampk
      scL(k)=scLk
      Bet(k)=Betk
    end do

    do j=1,nveclat
      rlat=veclats(j)*RAD_PER_DEG
      slat=sin(rlat)

      do k=1,nveclev

	eps=Amp(k)*exp(-(rlat/scL(k))*(rlat/scL(k)))
	Edet=eps*eps+slat*slat

	Aum_imat(k,j) =       -slat/Edet
	Aul_imat(k,j) = -Bet(k)*eps/Edet
	Avm_imat(k,j) = -Bet(k)*eps/Edet
	Avl_imat(k,j) =        slat/Edet

      end do
    end do

	if(mall_ison()) then
	  call mall_mco(ilev,myname)
	  call mall_mco(wlev,myname)
	  call mall_mco( Amp,myname)
	  call mall_mco( scL,myname)
	  call mall_mco( Bet,myname)
	endif

    deallocate(ilev,wlev,Amp,scL,Bet)
  end subroutine EXP_ALPHAS_model_
!_______________________________________________________________________
  subroutine EXP_ALPHAS_model_x1

	! This model assume a latitude dependent epsilon value that
	! varies with respect to pressure levels too.

    use const, only	: OMEGA, ONEPI
    use m_die, only 	: die
    use m_mall,only	: mall_ison,mall_mci,mall_mco
    implicit none

    integer j,k,il,istat
    real wt,slat,rlat
    real Ampk,scLk,Betk,Awtk,Bwtk,Edet,eps
    real,    allocatable :: Amp(:),scL(:),Bet(:),Awt(:),Bwt(:)
    integer, allocatable :: ilev(:)
    real,    allocatable :: wlev(:)

    real, parameter	:: RAD_PER_DEG	= ONEPI/180.

	! index and assign log-linear weights to the table

    allocate(ilev(1:nveclev),wlev(1:nveclev),		&
	      Amp(1:nveclev), scL(1:nveclev),		&
	      Bet(1:nveclev), Awt(1:nveclev),		&
	      Bwt(1:nveclev), stat=istat		)
    if(istat.ne.0)	&
      call die(myname,'EPSILONS_model_x1 allocate()',istat)

	if(mall_ison()) then
	  call mall_mci(ilev,myname)
	  call mall_mci(wlev,myname)
	  call mall_mci( Amp,myname)
	  call mall_mci( scL,myname)
	  call mall_mci( Bet,myname)
	  call mall_mci( Awt,myname)
	  call mall_mci( Bwt,myname)
	endif

    call slogtab(.not.nearest, FEalpha_nlev,FEalpha_plev,	&
		nveclev,pveclev, ilev,wlev		)

    do k=1,nveclev
      il=ilev(k)		! conditioned log-linear interpolation
      wt=wlev(k)

      Ampk=FEalpha_pars(1,il)
      scLk=FEalpha_pars(2,il)
      Betk=FEalpha_pars(3,il)
      Awtk=FEalpha_pars(4,il)
      Bwtk=FEalpha_pars(5,il)

      if(il.lt.FEalpha_nlev) then
	Ampk=Ampk+wt*(FEalpha_pars(1,il+1)-Ampk)
	scLk=scLk+wt*(FEalpha_pars(2,il+1)-scLk)
	Betk=Betk+wt*(FEalpha_pars(3,il+1)-Betk)
	Awtk=Awtk+wt*(FEalpha_pars(4,il+1)-Awtk)
	Bwtk=Bwtk+wt*(FEalpha_pars(5,il+1)-Bwtk)
      endif

      Amp(k)=Ampk
      scL(k)=scLk
      Bet(k)=Betk
      Awt(k)=Awtk
      Bwt(k)=Bwtk
    end do

    do j=1,nveclat
      rlat=veclats(j)*RAD_PER_DEG
      slat=sin(rlat)

      do k=1,nveclev

	eps=Amp(k)*exp(-(rlat/scL(k))*(rlat/scL(k)))
	Edet=eps*eps+slat*slat

			! For test purpose only

	Aum_imat(k,j) = -slat*(Awt(k) + Bwt(k)/Edet)
	Aul_imat(k,j) = -Bet(k)*eps/Edet
	Avm_imat(k,j) = -Bet(k)*eps/Edet
	Avl_imat(k,j) =  slat*(Awt(k) + Bwt(k)/Edet)

      end do
    end do

	if(mall_ison()) then
	  call mall_mco(ilev,myname)
	  call mall_mco(wlev,myname)
	  call mall_mco( Amp,myname)
	  call mall_mco( scL,myname)
	  call mall_mco( Bet,myname)
	  call mall_mco( Awt,myname)
	  call mall_mco( Bwt,myname)
	endif
    deallocate(ilev,wlev,Amp,scL,Bet,Awt,Bwt)
  end subroutine EXP_ALPHAS_model_x1
!_______________________________________________________________________
  subroutine EXP_ALPHAS_model_x2

	! This model assumes a latitude dependent epsilon value that
	! varies with respect to pressure levels too.
        ! The formula are:
        ! 
        !  G(x)=b*(1-exp(-x^2/K^2))/sin(x), |x| > 0, 
        !  G(x)=0.0,                         x  = 0,
        !  e(x)=A + B*exp(-x^2/L^2),
        !  D(x)=e(x).
        !  
        !  u = -G(x)(dh/dy) - D(x)(dh/dx),
        !  v = -D(x)(dh/dy) + G(x)(dh/dx).
        !  
        !  Guang Ping Lou, Oct.1, 1998  
        !  Guang Ping Lou, Nov.5, 1998  
        !  
!_______________________________________________________________________
    use const, only	: OMEGA, ONEPI
    use m_die, only 	: die
    use m_mall,only	: mall_ison,mall_mci,mall_mco
    implicit none

    integer j,k,il,istat
    real wt,slat,rlat
    real Ampk,scLk,Betk,scKk,Bwtk,G,eps
    real,    allocatable :: Amp(:),scL(:),Bet(:),scK(:),Bwt(:)
    integer, allocatable :: ilev(:)
    real,    allocatable :: wlev(:)

    real, parameter	:: RAD_PER_DEG	= ONEPI/180.

	! index and assign log-linear weights to the table

    allocate(ilev(1:nveclev),wlev(1:nveclev),		&
	      Amp(1:nveclev), scL(1:nveclev),		&
	      Bet(1:nveclev), scK(1:nveclev),		&
	      Bwt(1:nveclev), stat=istat		)
    if(istat.ne.0)	&
      call die(myname,'EPSILONS_model_x1 allocate()',istat)

	if(mall_ison()) then
	  call mall_mci(ilev,myname)
	  call mall_mci(wlev,myname)
	  call mall_mci( Amp,myname)
	  call mall_mci( scL,myname)
	  call mall_mci( Bet,myname)
	  call mall_mci( scK,myname)
	  call mall_mci( Bwt,myname)
	endif

    call slogtab(.not.nearest, FEalpha_nlev,FEalpha_plev,	&
		nveclev,pveclev, ilev,wlev		)

    do k=1,nveclev
      il=ilev(k)		! conditioned log-linear interpolation
      wt=wlev(k)

      Ampk=FEalpha_pars(1,il)
      scLk=FEalpha_pars(2,il)
      Betk=FEalpha_pars(3,il)
      scKk=FEalpha_pars(4,il)
      Bwtk=FEalpha_pars(5,il)

      if(il.lt.FEalpha_nlev) then
	Ampk=Ampk+wt*(FEalpha_pars(1,il+1)-Ampk)
	scLk=scLk+wt*(FEalpha_pars(2,il+1)-scLk)
	Betk=Betk+wt*(FEalpha_pars(3,il+1)-Betk)
	scKk=scKk+wt*(FEalpha_pars(4,il+1)-scKk)
	Bwtk=Bwtk+wt*(FEalpha_pars(5,il+1)-Bwtk)
      endif

      Amp(k)=Ampk
      scL(k)=scLk
      Bet(k)=Betk
      scK(k)=scKk
      Bwt(k)=Bwtk
    end do

    do j=1,nveclat
      rlat=veclats(j)*RAD_PER_DEG
      slat=sin(rlat)

        !  G(x)=b*(1-exp(-x^2/K^2))/sin(x), |x| > 0, 
        !  G(x)=0.0,                         x  = 0,
        !  e(x)=A + B*exp(-x^2/L^2),
        !  D(x)=e(x).
        !  
        !  u = -Aum(x)(dh/dy) - Aul(x)(dh/dx),
        !  v = -Avm(x)(dh/dy) + Avl(x)(dh/dx).
        !  The "-" signs are moved to inside of the parameters  
        !  in the following calculations.
        !  In the calculations of u and v, they will be
        !  multiplied by a factor of g/(2*omega).
        !  
      do k=1,nveclev

	eps=Amp(k) + Bet(k)*exp(-(rlat/scL(k))*(rlat/scL(k)))
	G=Bwt(k)*(1.0-exp(-(rlat/scK(k))*(rlat/scK(k))))

       ! Compute the geostrophic(ageostropic) coefficients

       if ( abs(slat) .ge. 0.005 ) then
       
          Aum_imat(k,j) = -G/slat
          Avl_imat(k,j) =  G/slat
         else
       
          Aum_imat(k,j) = 0.0
          Avl_imat(k,j) = 0.0
        endif
	Aul_imat(k,j) = -eps
	Avm_imat(k,j) = -eps

      end do
    end do

	if(mall_ison()) then
	  call mall_mco(ilev,myname)
	  call mall_mco(wlev,myname)
	  call mall_mco( Amp,myname)
	  call mall_mco( scL,myname)
	  call mall_mco( Bet,myname)
	  call mall_mco( scK,myname)
	  call mall_mco( Bwt,myname)
	endif
    deallocate(ilev,wlev,Amp,scL,Bet,scK,Bwt)
  end subroutine EXP_ALPHAS_model_x2
!-----------------------------------------------------------------------
  subroutine EXP_ALPHAS_model_x3

	! A model similar to PSAS:X0, but allowing the frictional component
	! to have a latitudinal independent portion.
        ! Note: See ATBD 1998 for documentation of PSAS:X0

        !  u = -G(x)(dh/dy) - D(x)(dh/dx),
        !  v = -D(x)(dh/dy) + G(x)(dh/dx).

        !  e(x) = A * exp ( -x^2/L^2 ) 
        !  G(x) = sin(y) / ( sin(x)^2 + e(x)^2 )
        !  D(x) = C + B * e(x) / ( sin(x)^2 + e(x)^2 )
        !


    use const, only	: OMEGA, ONEPI
    use m_die, only 	: die
    implicit none

    integer j,k,il,istat
    real wt,slat,rlat
    real Ampk,scLk,Betk,Cetk,Edet,eps
    real,    allocatable :: Amp(:),scL(:),Bet(:), Cet(:)
    integer, allocatable :: ilev(:)
    real,    allocatable :: wlev(:)

    real, parameter	:: RAD_PER_DEG	= ONEPI/180.

	! index and assign log-linear weights to the table

    allocate(ilev(1:nveclev),wlev(1:nveclev),		&
	      Amp(1:nveclev), scL(1:nveclev),		&
	      Bet(1:nveclev), Cet(1:nveclev),  stat=istat		)
    if(istat.ne.0)	&
      call die(myname,'EPSILONS_model_ allocate()',istat)

    call slogtab(.not.nearest, FEalpha_nlev,FEalpha_plev,	&
		nveclev,pveclev, ilev,wlev		)

    do k=1,nveclev
      il=ilev(k)		! conditioned log-linear interpolation
      wt=wlev(k)

      Ampk=FEalpha_pars(1,il)
      scLk=FEalpha_pars(2,il)
      Betk=FEalpha_pars(3,il)
      Cetk=FEalpha_pars(4,il)
      if(il.lt.FEalpha_nlev) then
	Ampk=Ampk+wt*(FEalpha_pars(1,il+1)-Ampk)
	scLk=scLk+wt*(FEalpha_pars(2,il+1)-scLk)
	Betk=Betk+wt*(FEalpha_pars(3,il+1)-Betk)
	Cetk=Cetk+wt*(FEalpha_pars(4,il+1)-Cetk)
      endif

      Amp(k)=Ampk
      scL(k)=scLk
      Bet(k)=Betk
      Cet(k)=Cetk

    end do

    do j=1,nveclat
      rlat=veclats(j)*RAD_PER_DEG
      slat=sin(rlat)

      do k=1,nveclev

	eps=Amp(k)*exp(-(rlat/scL(k))*(rlat/scL(k)))
	Edet=eps*eps+slat*slat

	Aum_imat(k,j) =       -slat/Edet
	Aul_imat(k,j) = - ( Cet(k) + Bet(k)*eps/Edet )
	Avm_imat(k,j) = - ( Cet(k) + Bet(k)*eps/Edet )
	Avl_imat(k,j) =        slat/Edet

      end do
    end do

    deallocate(ilev,wlev,Amp,scL,Bet)
  end subroutine EXP_ALPHAS_model_x3
!_______________________________________________________________________
end subroutine imat_alpha
!.
