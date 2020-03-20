! file: mod_simTauRe.f90
! simulate an emsemble of tauRe subcolumns for a GSM

module mod_simTauRe

  use mod_utils
  use mod_StaticGridcolumnStatisticalModel
  implicit none

  private
  public simTauRe
  public simTauRe_create, simTauRe_destroy
  public simTauRe_generate
  public getftauVis

  ! ~~~~~~~~~~~
  type simTauRe
  ! ~~~~~~~~~~~
    !private

    ! pointer to the underlying GSM
    type (GridcolumnStatisticalModel), pointer :: gsm

    ! pre-calcs for all layers (1:gsm%nk)
    real*4, allocatable, dimension(:) :: &
      rLeff, rIeff, delpog, cgTf1VisLiq, cgTf1VisIce

    ! pre-calcs for top (non-variable) layers (1:gsm%klm1)
    real*4, allocatable, dimension(:) :: &
      top_dcotLiq, top_dcotIce  ! lyr cld opt thicknesses

  end type simTauRe

contains

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine simTauRe_create (s, gsm, pref, u, qcils, qcian)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! set temperature dependent pre-calcs, etc.

    use mod_consts
    use mod_reff, only: get_tempor, reff
    implicit none

    type (simTauRe), intent(out) :: s
    type (GridcolumnStatisticalModel), intent(in), target :: gsm

    ! these four (pref, u, qcils, qcian) are used in eff. radius calc. only
    real*4, intent(in) :: pref (gsm%nk+1)  ! reference edge pressures [Pa]
    real*4, intent(in), dimension(gsm%nk) :: &
      u,      & ! layer zonal (eastward) velocity [m/s]
      qcils,  & ! large scale ice water content component [kg/kg]
      qcian     ! anvil       ice water content component [kg/kg]

    integer :: k
    real*4 :: tempor
    real*4, dimension(gsm%nk) :: &
      ftauLqdVis, ftauIceVis
!     ftauLqdIR,  ftauIceIR

    ! load basics
    call simTauRe_destroy (s)
    s%gsm => gsm

    ! effective radii
    call get_tempor (gsm%nk, pref, u, tempor)
    allocate(s%rLeff(gsm%nk),s%rIeff(gsm%nk))
!   call getreff_simple (gsm%pMid, gsm%T, gsm%nk, s%rLeff, s%rIeff)
    do k = 1, gsm%nk
      call reff (gsm%T(k), gsm%pMid(k)/100., qcils(k), qcian(k), tempor, &
        s%rLeff(k), s%rIeff(k))
    end do

    ! Visible COT per condensed water path
    call getftauVis (s%rLeff, s%rIeff, gsm%nk, ftauLqdVis, ftauIceVis)

!   ! IR COT per condensed water path
!   call getftauIR  (s%rLeff, s%rIeff, gsm%nk, ftauLqdIR,  ftauIceIR)

    ! cloud optical thickness:
    !    dCOT ~ ftau * dCWP,
    !    dCWP = rhoc * dz = scon * rho * dz = scon * delp / grav,
    ! => dCOT = scon * cgTf1,  cgTf1 = ftau * delp / grav
    ! Note: Because the phase split is only currently a function
    ! of T (fice v. fliq in gsm), we may apply the phase split in
    ! cgTf1 instead of scon.
    allocate(s%delpog(gsm%nk))
    s%delpog = gsm%delp / grav
    allocate(s%cgTf1VisLiq(gsm%nk),s%cgTf1VisIce(gsm%nk))
    s%cgTf1VisLiq = ftaulqdVis * gsm%fliq * s%delpog
    s%cgTf1VisIce = ftauiceVis * gsm%fice * s%delpog

    ! top dcot
    if (gsm%klm1 .ge. 1) then
      allocate(s%top_dcotLiq(gsm%klm1),s%top_dcotIce(gsm%klm1))
      s%top_dcotLiq = gsm%top_scon * s%cgTf1VisLiq(1:gsm%klm1)
      s%top_dcotIce = gsm%top_scon * s%cgTf1VisIce(1:gsm%klm1)
    end if

  end subroutine simTauRe_create


  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine simTauRe_destroy (s)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    type (simTauRe), intent(out) :: s
    s%gsm => NULL()
    if (allocated (s%delpog))       deallocate (s%delpog)
    if (allocated (s%rLeff))        deallocate (s%rLeff)
    if (allocated (s%rIeff))        deallocate (s%rIeff)
    if (allocated (s%cgTf1VisLiq))  deallocate (s%cgTf1VisLiq)
    if (allocated (s%cgTf1VisIce))  deallocate (s%cgTf1VisIce)
    if (allocated (s%top_dcotLiq))  deallocate (s%top_dcotLiq)
    if (allocated (s%top_dcotIce))  deallocate (s%top_dcotIce)
  end subroutine simTauRe_destroy


  ! tauRe simulator
  ! produces ngen subcolumns of tau and re
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine simTauRe_generate(s, ngen, taul, taui, rel, rei)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    implicit none

    ! inputs ...
    type (simTauRe), intent(in) :: s
    integer, intent(in) :: ngen

    ! outputs ...
    real*4, dimension(ngen,s%gsm%nk), intent(out) :: &
      taul, taui, rel, rei

    ! locals ...
    integer :: k
    type (MoistureSubcolumnEnsemble) :: mse

    ! copy static effective radii
    rel = spread(s%rLeff,1,ngen)
    rei = spread(s%rIeff,1,ngen)

    ! generate subcolumns
    call GridcolumnStatisticalModel_generate (s%gsm, ngen, mse)

    ! Calculate subcolumn COTs
    if (s%gsm%klm1 .ge. 1) then
      taul(:,1:s%gsm%klm1) = spread(s%top_dcotLiq,1,ngen)
      taui(:,1:s%gsm%klm1) = spread(s%top_dcotIce,1,ngen)
    end if
    do k = 1, s%gsm%nkv
      taul(:,s%gsm%klm1+k) =  mse%scon(:,k) * s%cgTf1VisLiq(s%gsm%klm1+k)
      taui(:,s%gsm%klm1+k) =  mse%scon(:,k) * s%cgTf1VisIce(s%gsm%klm1+k)
    end do

  end subroutine simTauRe_generate

  ! estimate liquid and ice effective radii
  ! (simplified from GEOS-5 version, and p now in Pa)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine getreff_simple(p, T, nk, rLeff, rIeff)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    implicit none
    integer, intent(in) :: nk
    ! pressure p in Pa, temperature T in K
    real*4, intent(in),  dimension(nk) :: p, T
    ! effective radii in m
    real*4, intent(out), dimension(nk) :: rLeff, rIeff

    ! limits
    real*4, parameter :: rImin   = 20.0e-6
    real*4, parameter :: rImax   = 40.0e-6
    real*4, parameter :: rLmin   = 10.0e-6
    real*4, parameter :: rLmax   = 21.0e-6

    ! Ice
    where (p <= 150.e2)
      rIeff = rImax
    elsewhere
      rIeff = rImax * 150.e2 / p
    end where
    where (rIeff < rImin) rIeff = rImin

    ! Liquid
    where (p <= 300.e2)
      rLeff = rLmax
    elsewhere
      rLeff = rLmax * 300.e2 / p
    end where
    where (rLeff < rLmin) rLeff = rLmin

  end subroutine getreff_simple


  ! visible COT per layer per phase per condensate water path.
  ! Multiply externally by cwp [in kg/m2] in each phase to get COTs.
  ! NB: users responsibility to make sure the reff > 0.
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine getftauVis(re_lqd, re_ice, nk, ftau_lqd, ftau_ice)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    implicit none
    integer, intent(in) :: nk
    real*4, intent(in),  dimension(nk) :: &
       re_lqd,  re_ice    ! effective radii in [m]
    real*4, intent(out), dimension(nk) :: &
      ftau_lqd, ftau_ice  ! (tau / cwp) in [m2/kg]
  
    ! exctinction coefficient constants
    real*4, parameter :: awb1 = -6.59  ! m2/kg
    real*4, parameter :: awb2 =  1.65e-3 ! m3/kg
    real*4, parameter :: aib  =  1.64e-3 ! m3/kg

    ! optical depths per cwp
    ftau_lqd = awb1 + awb2 / re_lqd
    ftau_ice = aib / re_ice
  
    ! make sure non-negative
    where(ftau_ice < 0.) ftau_ice = 0.
    where(ftau_lqd < 0.) ftau_lqd = 0.

  end subroutine getftauVis


  ! 10.5 um longwave COT per layer per phase per condensate water path.
  ! Multiply externally by cwp [in kg/m2] in each phase to get COTs.
  ! PS: ice and lqd only, no rain.
  ! References are to Chou et al, 2001, NASA/TM-2001-104606, Vol. 19.
  ! NB: users responsibility to make sure the reff > 0.
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine getftauIR(re_lqd, re_ice, nk, ftau_lqd, ftau_ice)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    implicit none
    integer, intent(in) :: nk
    real*4, intent(in),  dimension(nk) :: &
      re_lqd, re_ice      ! effective radii in [m]
    real*4, intent(out), dimension(nk) :: &
      ftau_lqd, ftau_ice  ! (tau / cwp) in [m2/kg]

    ! coefficients for computing the extinction coefficient
    ! for cloud ice particles (Table 11a, Eq. 6.4a).
    real*4, parameter :: aib(3) = &
      (/ -0.01896, 0.78955, 0.69493 /)

    ! coefficients for computing the extinction coefficient
    ! for cloud liquid drops. (Table 11b, Eq. 6.4b)
    real*4, parameter :: awb(4) = &
      (/ 0.15587, 0.00371, -7.7705e-4, 2.0547e-5 /)

    ! locals
    real*4, dimension(nk) :: rl, ri

    ! convert effective radii to microns
    ! (as equations below require. See Chou et al.)
    ri = re_ice * 1.e6; rl = re_lqd * 1.e6

    ! compute cloud optical thickness -- Eqs. (6.4a,b) and (6.7)
    ! original formula gave m2/g, so multiply by 1000. to get m2/kg
    ftau_ice = 1.e3 * (aib(1)+aib(2)/ri**aib(3))
    ftau_lqd = 1.e3 * (awb(1)+(awb(2)+(awb(3)+awb(4)*rl)*rl)*rl)

  end subroutine getftauIR

end module mod_simTauRe
