!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: imat_sigW - determine mass-decoupled wind forecast errors
!
! !INTERFACE:

    subroutine imat_sigW()
      use config, only	: kmat_STRM, kmat_VELP
      use m_mvcorF_blox,only : norm_DD
      use m_chars, only : uppercase
      use rlat_imat
      use rlev_imat
      use FEsigW_tabl
      use FEsigW_imat
      use m_die,only : die
      use m_stdio, only : stderr
      implicit none

! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	24Mar96 - J. Guo	- (to do)
!_______________________________________________________________________
!#define	_DEBUG	!

	! parameters

  character(len=*), parameter	:: myname='imat_sigW'
  logical,	    parameter	:: nearest=.true.	! a flag

	! Local vars.
  integer :: j,k,ln,istat,mpar
!-----------------------------------------------------------------------
	! Consistency checking

  select case (uppercase(FEsigW_type))
  case (' ','PSAS1.1','PSAS1.1A')
    mpar=2
  case ('GEOS/DAS-OI')
    mpar=0
  case ('PSAS/S&L')
    mpar=3
  case ('NULL')
    mpar=0
  case default
    write(stderr,'(5a)') myname,': unknown type of "',	&
      trim(FEsigW_rsrc),'", "',trim(FEsigW_type),'"'
    call die(myname)
  end select

  if(FEsigW_npar < mpar) then
    write(stderr,'(2a,i1,3a,i4)') myname,': expecting ',mpar,	&
      ' parameters for "',FEsigW_type,'", found only =',FEsigW_npar
    call die(myname)
  endif
!-----------------------------------------------------------------------
	! Interpolation

  select case (uppercase(FEsigW_type))
  case ('PSAS/S&L')
		! sigS and sigV are lat & level dependent

    call PSAS_sigW_SL()

  case ('PSAS1.1A')
		! sigS and sigV are lat & level dependent

    call PSAS11_sigW_a()

  case (' ','PSAS1.1')
		! sigS and sigV are only level dependent

    call PSAS11_sigW_()

  case ('GEOS/DAS-OI')

    call OI_sigW_()

  case ('NULL')

    FEsigS_imat(1:nveclev,1:nveclat)=0.
    FEsigV_imat(1:nveclev,1:nveclat)=0.

  end select
!=======================================================================

#ifdef	_DEBUG
_DEBUG	write(*,'(2a)') myname,'*DEBUG::'
_DEBUG	ln=max(len_trim(FEsigW_desc),1)
_DEBUG	write(*,'(a,2x,a)') FEsigW_type,FEsigW_desc(1:ln)
_DEBUG	write(*,'(a)') ' lv  plev  FEsigS  FEsigV'
_DEBUG	do k=1,nveclev,(nveclev-1)/2
_DEBUG	  write(*,*) 'level = ',k,pveclev(k),'  FEsigW  FEsigV'
_DEBUG	  do j=1,nveclat
_DEBUG	    write(*,'(f7.2,2f8.2)') veclats(j),		&
_DEBUG		FEsigS_imat(k,j),FEsigV_imat(k,j)
_DEBUG	  end do
_DEBUG	end do
_DEBUG	write(*,'(a)')  '::'
#endif

contains
!-----------------------------------------------------------------------
  subroutine OI_sigW_()

	! This model is made similar to what used in OI

    use const, only	: R_EARTH, G_EARTH, OMEGA, ONEPI
    use m_die, only	: die
    use m_mall,only	: mall_ison,mall_mci,mall_mco
    implicit none

    integer j,lv,istat
    real coriolis,wt
    real rlat_r
    real bR_sfc,bD_sfc,bC_sfc,bF_sfc,bC_Eqs,bF_Eqs
    real bR_upa,bD_upa
    real Alf,Hfac

    integer, parameter	:: Trop_nlev = 11
    real, parameter	:: Trop_plev(Trop_nlev) =	&
	(/.01, 40., 60.,175.,225.,275.,350.,450.,600.,925.,1000. /)
    real, parameter	:: Trop_sigW(Trop_nlev) =	&
	(/7.1, 8.0, 7.3, 7.4, 7.0, 6.2, 5.2, 4.5, 4.4, 4.3, 5.0  /)

    real, parameter	:: reflat	= 25.
    real, parameter	:: dellat	=  5.

    real, parameter	:: RAD_PER_DEG	= ONEPI/180.
    real, parameter	:: WDRG		= 1.e-5/(2.*OMEGA)

    integer, allocatable :: ilev(:)
    real,    allocatable :: wlev(:)
!-----------------------------------------------------------------------

	! index and assign log-linear weights to the table

    allocate(ilev(1:nveclev),wlev(1:nveclev), stat=istat)
    if(istat.ne.0)	&
      call die(myname,'OI_sigW_() allocate()',istat)

	if(mall_ison()) then
	  call mall_mci(ilev,myname)
	  call mall_mci(wlev,myname)
	endif

    call slogtab(.not.nearest, Trop_nlev,Trop_plev,	&
		nveclev,pveclev, ilev,wlev		)

    do j=1,nveclat

      rlat_r=abs(veclats(j))*RAD_PER_DEG

      if(abs(veclats(j)).lt.reflat-dellat) then

	Alf=sin((reflat-dellat)*RAD_PER_DEG)

      elseif(abs(veclats(j)).lt.reflat+dellat) then

	wt=.5*(1.+sin(ONEPI*(abs(veclats(j))-reflat)/(2.*dellat)))
	Alf=(1.-wt)*sin((reflat-dellat)*RAD_PER_DEG) +	&
		wt *sin(rlat_r)

      else

	Alf=sin(rlat_r)

      endif

      Hfac=1.-exp(-abs(veclats(j))/reflat)

      bR_upa=sin((reflat-dellat)*RAD_PER_DEG)*sqrt(1.-Hfac*Hfac)/Alf
      bD_upa = 0.

      bC_sfc = Alf/(Alf*Alf+WDRG*WDRG)
      bF_sfc =WDRG/(Alf*Alf+WDRG*WDRG)

      Alf    =sin((reflat-dellat)*RAD_PER_DEG)
      bC_Eqs = Alf/(Alf*Alf+WDRG*WDRG)
      bF_Eqs =WDRG/(Alf*Alf+WDRG*WDRG)

      bR_sfc=bC_sfc/sqrt(bC_Eqs*bC_Eqs+bF_Eqs*bF_Eqs)
      bR_sfc=sqrt(1.-Hfac*Hfac)*bR_sfc

      bD_sfc=bF_sfc/sqrt(bC_Eqs*bC_Eqs+bF_Eqs*bF_Eqs)
      bD_sfc=sqrt(1.-Hfac*Hfac)*bD_sfc

      do lv=1,nveclev
	if(pveclev(lv).le.850.) then
	  FEsigS_imat(lv,j)=bR_upa*Trop_sigW(ilev(lv))
	  FEsigV_imat(lv,j)=bD_upa*Trop_sigW(ilev(lv))
	elseif(pveclev(lv).lt.1000.) then
	  wt=log(pveclev(lv)/1000.)/log(850./1000.)
	  FEsigS_imat(lv,j)=(bR_sfc+wt*(bR_upa-bR_sfc))*	&
	    Trop_sigW(ilev(lv))
	  FEsigV_imat(lv,j)=(bD_sfc+wt*(bD_upa-bD_sfc))*	&
	    Trop_sigW(ilev(lv))
	else
	  FEsigS_imat(lv,j)=bR_sfc*Trop_sigW(ilev(lv))
	  FEsigV_imat(lv,j)=bD_sfc*Trop_sigW(ilev(lv))
	endif
      end do

!     do lv=1,nveclev
!	FEsigS_imat(lv,j)=bR_upa*Trop_sigW(ilev(lv))
!	FEsigV_imat(lv,j)=0.
!     end do

    end do

	if(mall_ison()) then
	  call mall_mco(ilev,myname)
	  call mall_mco(wlev,myname)
	endif
    deallocate(ilev,wlev)
    
  end subroutine OI_sigW_
!-----------------------------------------------------------------------
  subroutine PSAS11_sigW_

	! This model assume a latitude independent sigS and sigV

    use const, only	: R_EARTH, G_EARTH, OMEGA
    use m_die, only	: die
    use m_mall,only	: mall_ison,mall_mci,mall_mco
    implicit none

    real, parameter	:: scalar = G_EARTH/(2.*OMEGA*R_EARTH)

    integer k,lv,istat
    real wt,sigS,sigV
    integer, allocatable :: ilev(:)
    real,    allocatable :: wlev(:)

	! index and assign log-linear weights to the table

    allocate(ilev(1:nveclev),wlev(1:nveclev), stat=istat)
    if(istat.ne.0)	&
      call die(myname,'PSAS11_sigW_() allocate()',istat)
	if(mall_ison()) then
	  call mall_mci(ilev,myname)
	  call mall_mci(wlev,myname)
	endif

    call slogtab(.not.nearest, FEsigW_nlev,FEsigW_plev,	&
		nveclev,pveclev, ilev,wlev		)

    do k=1,nveclev
      lv=ilev(k)		! conditioned log-linear interpolation
      wt=wlev(k)

	! For given sigS and sigV, interpolate to table levels.

      sigS=FEsigW_pars(1,lv)
      sigV=FEsigW_pars(2,lv)

      if(lv.lt.FEsigW_nlev) then
	sigS=sigS+wt*(FEsigW_pars(1,lv+1)-sigS)
	sigV=sigV+wt*(FEsigW_pars(2,lv+1)-sigV)
      endif

	! Convert to wind units.  norm_DD is already on table levels.

      FEsigS_imat(k,1:nveclat)=sigS*scalar/norm_DD(k,kmat_STRM)
      FEsigV_imat(k,1:nveclat)=sigV*scalar/norm_DD(k,kmat_VELP)
    end do

	if(mall_ison()) then
	  call mall_mco(ilev,myname)
	  call mall_mco(wlev,myname)
	endif
    deallocate(ilev,wlev)
  end subroutine PSAS11_sigW_
!_______________________________________________________________________
  subroutine PSAS11_sigW_a

	! This model assume a latitude independent sigS and sigV

    use const, only	: R_EARTH, G_EARTH, OMEGA
    use m_die, only	: die
    use m_mall,only	: mall_ison,mall_mci,mall_mco
    implicit none

    real, parameter	:: scalar = G_EARTH/(2.*OMEGA*R_EARTH)

    integer k,lv,istat
    real wt,sigS,sigV
    integer, allocatable :: ilev(:)
    real,    allocatable :: wlev(:)

    real x,x2,x4,scalar_a

    real,parameter :: pbltop=700.
    real,parameter :: pblsfc=1000.

	! index and assign log-linear weights to the table

    allocate(ilev(1:nveclev),wlev(1:nveclev), stat=istat)
    if(istat.ne.0)	&
      call die(myname,'PSAS11_sigW_a() allocate()',istat)

	if(mall_ison()) then
	  call mall_mci(ilev,myname)
	  call mall_mci(wlev,myname)
	endif

    call slogtab(.not.nearest, FEsigW_nlev,FEsigW_plev,	&
		nveclev,pveclev, ilev,wlev		)

    do j=1,nveclat
      do k=1,nveclev
        lv=ilev(k)		! conditioned log-linear interpolation
        wt=wlev(k)

	! For given sigS and sigV, interpolate to table levels.

        sigS=FEsigW_pars(1,lv)
        sigV=FEsigW_pars(2,lv)

        if(lv.lt.FEsigW_nlev) then
	  sigS=sigS+wt*(FEsigW_pars(1,lv+1)-sigS)
	  sigV=sigV+wt*(FEsigW_pars(2,lv+1)-sigV)
        endif

	! Convert to wind units.  norm_DD is already on table levels.

	x=veclats(j)/25.
	x2=x*x
	x4=x2*x2
	wt=max(0.,min(log(pblsfc/pveclev(k))/log(pblsfc/pbltop),1.))
        scalar_a=scalar*( 1.-wt*(1.-exp(-x4)) )

        FEsigS_imat(k,j)=sigS*scalar_a/norm_DD(k,kmat_STRM)
        FEsigV_imat(k,j)=sigV*scalar_a/norm_DD(k,kmat_VELP)

      end do
    end do

	if(mall_ison()) then
	  call mall_mco(ilev,myname)
	  call mall_mco(wlev,myname)
	endif
    deallocate(ilev,wlev)
  end subroutine PSAS11_sigW_a
!_______________________________________________________________________
  subroutine PSAS_sigW_SL

	! This model assume a latitude independent sigS and sigV
! The height decoupled winds are computed via:
!
! ud = + dChi/dx - dPsi/dy
! vd = + dPsi/dx + dChi/dx
! Cov(Chi) = sigChi.sigChi.corChi
! sigChi = sigWv.(scaling).exp(-x^2/J^2)
! sigPsi = sigWs.(scaling).exp(-x^2/J^2)
! where sigWs, sigWv, and J are read in parameters.


    use const, only	: R_EARTH, G_EARTH, OMEGA, ONEPI
    use m_die, only	: die
    use m_mall,only	: mall_ison,mall_mci,mall_mco
    implicit none

    real, parameter	:: scalar = G_EARTH/(2.*OMEGA*R_EARTH)
    real, parameter	:: RAD_PER_DEG	= ONEPI/180.

    integer k,lv,istat
    real wt,sigS,sigV,sigJ
    integer, allocatable :: ilev(:)
    real,    allocatable :: wlev(:)

    real x,scalar_a

	! index and assign log-linear weights to the table

    allocate(ilev(1:nveclev),wlev(1:nveclev), stat=istat)
    if(istat.ne.0)	&
      call die(myname,'PSAS11_sigW_a() allocate()',istat)

	if(mall_ison()) then
	  call mall_mci(ilev,myname)
	  call mall_mci(wlev,myname)
	endif

    call slogtab(.not.nearest, FEsigW_nlev,FEsigW_plev,	&
		nveclev,pveclev, ilev,wlev		)

    do j=1,nveclat
      do k=1,nveclev
        lv=ilev(k)		! conditioned log-linear interpolation
        wt=wlev(k)

	! For given sigS and sigV, interpolate to table levels.

        sigS=FEsigW_pars(1,lv)
        sigV=FEsigW_pars(2,lv)
        sigJ=FEsigW_pars(3,lv)

        if(lv.lt.FEsigW_nlev) then
	  sigS=sigS+wt*(FEsigW_pars(1,lv+1)-sigS)
	  sigV=sigV+wt*(FEsigW_pars(2,lv+1)-sigV)
	  sigJ=sigJ+wt*(FEsigW_pars(3,lv+1)-sigJ)
        endif

	! Convert to wind units.  norm_DD is already on table levels.

        x=veclats(j)*RAD_PER_DEG
        scalar_a=scalar*exp(-x*x/(sigJ*sigJ))

        FEsigS_imat(k,j)=sigS*scalar_a/norm_DD(k,kmat_STRM)
        FEsigV_imat(k,j)=sigV*scalar_a/norm_DD(k,kmat_VELP)

      end do
    end do

	if(mall_ison()) then
	  call mall_mco(ilev,myname)
	  call mall_mco(wlev,myname)
	endif

    deallocate(ilev,wlev)
  end subroutine PSAS_sigW_SL
!_______________________________________________________________________
end subroutine imat_sigW
!.
