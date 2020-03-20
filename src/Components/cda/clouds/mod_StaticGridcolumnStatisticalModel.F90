! file: mod_StaticGridcolumnStatisticalModel.f90
! A StaticGridcolumn Statistical Model (GSM) class.
! A GSM is a factory to produce MoistureSubcolumnEnsemble (MSE).
! Pressure, temperature and their children are fixed.
! An upper atmosphere section is also completely fixed.
! Moisture in the lower atmosphere is statistically variable.
! Use create() to make a GSM with profiles of moisture PDFs.
! Use generate() to produce MSEs.

module mod_StaticGridcolumnStatisticalModel

  use mod_mkl
  use mod_utils
  use mod_consts
  use mod_GaussianCopula
  use MAPL_SatVaporMod
  implicit none

  private
  public GridcolumnStatisticalModel
  public MoistureSubcolumnEnsemble
  public GridcolumnStatisticalModel_create
  public GridcolumnStatisticalModel_destroy
  public GridcolumnStatisticalModel_generate
  public OVERLAP_GCOP, MARGIN_TRIANGLE

  ! key flags:
  integer, parameter :: OVERLAP_GCOP = 1
  integer, parameter :: MARGIN_TRIANGLE = 1

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type GridcolumnStatisticalModel
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Notes:
    ! (1) parameter variability is forced to zero for p <= pLim
    ! (2) the q??? are moisture contents ala Norris et al., QJRMS (2008),
    ! i.e., densities normalized by the virtual density ala Norris.

    !private

    integer :: nk	   ! number of layers (1=top, nk=bot)

    real*4 :: pTop         ! top pressure of model [Pa]
    real*4 :: pLim         ! largest press at which parameter
                           !   variability is zero [hPa]

    ! maximum allowed total saturation ratios in different phases.
    ! For mixed phase regions, the Smax will scale between the two
    ! based on icefrac(T).
    real*4 :: Smax_liq, Smax_ice

    ! layer profiles, all dim(nk)
    real*4, allocatable, dimension (:) :: &
      ! ~~~~~
      ! these first set are fixed since Temperature and
      ! pressure do not participate in the variability
      ! ~~~~~
      delp,        &! press thickness [Pa]
      pTope,       &! top edge press [Pa]
      pBote,       &! bot edge press [Pa]
      pMid,        &! midpoint press [Pa]
      T,           &! temperature [K]
      fice, fliq,  &! frac of ice & liq condensate 
                    !   (currently uses icefrac(T))
      qs,          &! saturation vapor content
                    !   (function of T and p only)
      rqs,         &! 1/qs
      cgTf2,       &! rqs-reps
      Smax,        &! maximum allowed S=qt/qs
                    !   anywhere in PDF
      ! ~~~~~
      ! Reference PDF parameter set for whole column (see "x" params below).
      ! ~~~~~
      Smean,       &! mean totl satn ratio, qtmean/qs
      PL,          &! triangle skewness parameter
      alpha         ! normalized triangle width

    ! klim is highest of moisture variable layers (klim:nk)
    integer :: klim
    ! number of top (non-variable) layers = klim-1 = "klm1"
    ! also k = kk + klm1 converts variable layer
    !   index kk=1...nkv to full layer index.
    integer :: klm1
    ! num of variable lyrs (nk-klm1)
    integer :: nkv
   
    ! 'variable layer' profiles, use (1:g%nkv)
    real*4, allocatable, dimension (:) :: &
      xSmean,  &! Smean
      xPL,     &! PL
      xalpha    ! alpha

    ! interlayer correlations
    real*4, allocatable :: CorM (:,:)

    ! copula for layer overlap
    type(GaussianCopula) :: copula

    ! top (non-variable) MSE-type output
    real*4, dimension(:), allocatable :: &
      top_svap, top_scon

  end type GridcolumnStatisticalModel

  type MoistureSubcolumnEnsemble
    !private
    integer :: nSubcolumns            ! "ncols"
    ! specific (mass / mass of parcel) quantities
    real*4, allocatable :: svap(:,:)  ! vapor      (ncols,nkv)
    real*4, allocatable :: scon(:,:)  ! condensate (ncols,nkv)
  end type MoistureSubcolumnEnsemble

  ! a very small S used to convert the open range (0, Smax),
  ! for acceptable Smean, to a closed range [Seps, Smax-Seps].
  ! Should have
  !   0 < Seps << 1,
  ! and also
  !   Seps << min(Smax_liq,Smax_ice) - 1.
  real*4, parameter :: Seps = 1.0e-3

  ! type of margin
  integer, parameter :: margin = MARGIN_TRIANGLE

  ! type of overlap
  integer, parameter :: overlap = OVERLAP_GCOP

  ! specify the acceptable parameter ranges used for the triangles.
  real*4, parameter :: PLmin = 0.1, PLmax = 0.9
  real*4, parameter :: almin = 0., almax = 1.
  ! Smean uses [Seps,Smax-Seps], where Smax
  !   is calculated based on Smax_liq, Smax_ice.

  ! PDF initialization from background state ...
  ! Estimate a PDF initial condition based on bkg qvmean and f
  logical, parameter :: ICpdfEst = .true.
  ! Maximum permitted modifications to f for PDF initialization
  integer, parameter :: maxfmods = 100
  ! f tolerance relative to fhi-flo, e.g., 0.01
  real*4, parameter :: ftol = 0.01
  ! nominal half width of pdf(S)
  ! used in special cases, e.g., 0.2
  real*4, parameter :: hDe = 0.2

contains

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GridcolumnStatisticalModel_create( &
    g, pTop, pLim, Smax_liq, Smax_ice, &
    nk, delp, T, qtmean_in, qvmean_in, f_in, &
    delpCorr, delSmeanRiish)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Notes:
  ! (1) pLim is the largest pressure at which parameter
  !   variability is forced to zero.
  ! (2) delpCorr is a pressure 'decorrelation length'
  !   for building inter-layer correlations.
  ! (3) delSmeanRiish: if > 0, use Riishojgaard flow-dependent
  !   modification to inter-layer S correlation.

    use mod_icefrac
    use mod_eswi
    use mod_triangle
    implicit none

    type (GridcolumnStatisticalModel), intent(out) :: g
    integer, intent(in) :: nk
    real*4,  intent(in) :: &
      pTop, pLim, Smax_liq, Smax_ice, delpCorr, delSmeanRiish
    real*4,  intent(in), dimension (nk) :: &
      delp, T, qtmean_in, qvmean_in, f_in

    real*4, parameter :: &
      tpi = 8.*atan(1.), nhlogtpi = -log(tpi)/2.

    ! locals
    integer :: i, j, k, klim, nmods, rc
    logical :: found, fInRange, low, validSolution
    real*4 :: PLx, De, qc_qs, qv_qs, ftolc, SL, SH, PHx, &
      SM, cR, qv0h_qs, qch_qs, cRlo, cRhi, flo, fhi, cRx, &
      fgoal, lof, hif, hlogdet
    real*4, dimension (nk) :: &
      facq, qsl, qsi, qtmean, qvmean, qcmean, f, Svmean
    real*4, allocatable :: pMidv(:)
    real*4, dimension(:,:), allocatable :: &
      HcorM, Sva, HcovS

    ! basic checks
    if (PLmin <= 0.)    call myError('require PLmin >  0.')
    if (PLmax >= 1.)    call myError('require PLmax <  1.')
    if (PLmin >= 0.5)   call myError('require PLmin <  0.5')
    if (PLmax <= 0.5)   call myError('require PLmax >  0.5')
    if (almin < 0.)     call myError('require almin >=  0.')
    if (almax > 1.)     call myError('require almax <=  1.')
    if (almin >= almax) call myError('require almin < almax')
    if (ftol <= 0. .or. ftol >= 0.5) call myError('need ftol in (0,0.5)')
    if (hDe  <= 0. .or. hDe  >= 0.5) call myError('need hDe in (0,0.5)')

    ! destroy any existing g
    call GridcolumnStatisticalModel_destroy(g)

    ! ~~~~~~~~~~~~~~~~~~~~~~~
    ! Establish mean profiles
    ! ~~~~~~~~~~~~~~~~~~~~~~~

    ! allocate and populate input profiles
    g%nk = nk
    allocate (g%delp(nk)); g%delp = delp
    allocate (g%T   (nk)); g%T    = T
    qtmean = qtmean_in; qvmean = qvmean_in; f = f_in

    ! evaluate pBote, pTope, pMid
    g%pTop = pTop
    allocate (g%pTope(nk), g%pBote(nk), g%pMid(nk))
    g%pTope(1) = pTop
    g%pBote(1) = pTop + delp(1)
    do k = 2, nk
      g%pTope(k) = g%pBote(k-1)
      g%pBote(k) = g%pTope(k) + delp(k)
    end do
    g%pMid = (g%pBote + g%pTope) / 2.

    ! basic rectifications (qvmean later)
    where (qtmean < 0.) qtmean = 0.
    where (f < 0.)
      f = 0.
    elsewhere (f > 1.)
      f = 1.
    endwhere

    ! Per Norris et al., QJRMS 2008, qx = rhox / rhostar,
    ! where rhostar = p / (rgas * T), the virtual density,
    ! and rgas = 287.04 (dry air). In this system, the
    ! saturation vapor content is qs = eps * esat(T) / p.

    ! notation / meanings:
    ! qtmean = mean total water
    ! qs = constant vapor saturation
    ! Smean = qtmean/qs = mean total saturation ratio
    ! Smax = max value of S
    ! qvmean = mean water vapor
    ! qv0h = mean vapor in clear portion
    ! qcmean = mean condensate
    ! qch = mean in-cloud condensate

    ! store Smax_liq, Smax_ice in g
    g%Smax_liq = Smax_liq
    g%Smax_ice = Smax_ice

    ! qs and Smax calculation
    ! (qs using rhostar normalization)
    facq = eps / g%pMid
!   qsl = facq * esw_gcet(T)
!   qsi = facq * esi_gcet(T)
    qsl = facq * MAPL_EQsat(T,OverIce=.false.)
    qsi = facq * MAPL_EQsat(T,OverIce=.true. )
    allocate (g%fice(nk), g%fliq(nk))
    call icefrac (nk, T, g%fice, g%fliq)
    allocate (g%qs(nk), g%Smax(nk))
    g%qs    = g%fliq * qsl      + g%fice * qsi
    g%Smax  = g%fliq * Smax_liq + g%fice * Smax_ice

    ! Note: At some point, we should add a gain factor to qsi
    ! to account for ice supersaturation. Then qsi would be
    ! like an "effective ice vapor saturation". This would
    ! then allow Smax_ice to come down to similar to a value
    ! closer to Smax_liq.

    ! specific humidity
    ! svap = qv / (1 - (reps-1)*qv + qc)
    !      = qv / (1 + qt - reps*qv)
    ! By definition and by bulk assumption
    !   qt = S * qs;  qv = min(S,1) * qs
    ! svap =  min(S,1) / (rqs + S - reps*min(S,1))
    !   = 1 / (rqs/min(S,1) + S/min(S,1) - reps)
    !   = 1 / (rqs/S + omreps) for S <= 1
    !   = 1 / (cgTf2 + S)      for S >  1
    ! where rqs = 1/qs, omreps = 1-reps,
    !   and cgTf2 = rqs - reps
    allocate (g%rqs(nk), g%cgTf2(nk))
    g%rqs   = 1. / g%qs
    g%cgTf2 = g%rqs - reps

    ! Seps checks
    if (Seps <= 0. .or. Seps > 0.01) &
      call myError('need Seps in (0,0.01]')
    if (Seps > 0.01*(min(Smax_liq,Smax_ice)-1.)) &
      call myError('need Seps <= 0.01*(min(Smax_liq,Smax_ice)-1)')
      ! (since otherwise clouds will not properly form)

    ! mean total water saturation ratio clipped in (0,Smax)
    !   (open because must leave room for tails)
    allocate (g%Smean(nk))
    g%Smean = qtmean / g%qs
    where (g%Smean < Seps)
      g%Smean = Seps
    elsewhere (g%Smean > g%Smax - Seps)
      g%Smean = g%Smax - Seps
    endwhere
    qtmean = g%Smean * g%qs

    ! enforce physicality on mean qv:
    !   qcmean = qtmean - qvmean >= 0 so qvmean <= qtmean
    where (qvmean < 0)
      qvmean = 0.
    elsewhere (qvmean > qtmean)
      qvmean = qtmean
    endwhere

    ! enforce bulk assumption on mean vapor
    !   qvmean = (1-f) * qv0h + f * qs <= qs,
    ! since qv0h <= qs.
    where (qvmean > g%qs) qvmean = g%qs

    !  mean condensate
    qcmean = qtmean - qvmean

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Establish variable parameter region
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! parameter variation for index (klim:)
    g%pLim = pLim
    found = .false.
    do klim = 1, nk
      ! must not include equality:
      ! there must be some variability at index klim
      found = (g%pMid(klim) > g%pLim)
      if (found) exit
    end do
    if (.not.found) call myError('cant find pLim')
    g%klim = klim
    g%klm1 = klim - 1
    g%nkv = nk - g%klm1
    allocate (pMidv(g%nkv))
    pMidv = g%pMid(klim:nk)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Diagnosis of initial PDFs
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~

    ! set up initial margin parameters
    if (margin == MARGIN_TRIANGLE) then

      !!!!!!!!!!!!!!!!!!!!!!!!!
      !  SKEWED TRIANGLE PDF  !
      !!!!!!!!!!!!!!!!!!!!!!!!!

      ! Triangle with layers of PL, alpha and Smean
      allocate (g%PL(nk), g%alpha(nk))

      ! =======================
      ! non-variable atmosphere
      ! =======================

      ! Currently, for (1:klm1), we just assume a delta function at Smean.
      ! This is mainly for speed. May later diagnose a triangle there
      ! too, in which case this "non-variable" section will disappear.
      if (g%klm1 .ge. 1) then
        g%PL(1:g%klm1) = 0.5
        g%alpha(1:g%klm1) = 0.
      end if

      ! ===================
      ! variable atmosphere
      ! ===================

      if (.not.ICpdfEst) then

        g%PL(klim:nk) = 0.5  ! all symmetric
        g%alpha(klim:nk) = 0.5  ! half max width

      else

        ! Here we would like to diagnose a triangular PDF from Smean, qvmean and f
        ! We will trust these input variables in that order (i.e., Smean most)

        ! diagnose a layer at a time
        do k = klim, nk

          ! ............................................
          ! (A) No condensate in background (qcmean = 0)
          ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          ! This is taken as sufficient evidence of a clear case,
          ! regardless of f (which is considered less reliable).
          ! Note: qcmean = 0 implies qtmean = qvmean <= qs, so Smean <= 1

          ! Clear case
          if (qcmean(k) == 0.) then

            ! NB: f will be forced clear by update_fSv

            ! Want to use a symmetric triangle (PL=0.5, SM=Smean)
            ! of half width hDe, if it gives endpoints in [0,1]
            call triangleSymPlus(0.,1.,hDe,g%Smean(k),PLmin,PLmax,PLx,De)

            ! set parameters
            g%PL(k) = PLx
            g%alpha(k) = min(De/DeMaxScalar(g%Smean(k),g%Smax(k),PLx,1.-PLx),1.)

          ! .........................................
          ! (B) Condensate in background (qcmean > 0)
          ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          ! This is taken as sufficient evidence of a cloudy case,
          ! regardless of f (which is considered less reliable).
          ! In cloudy parcels, 
          !   qc = qt - qv = qt - qs > 0., so S = qt / qs > 1.

          ! Overcast case
          else if (qvmean(k) == g%qs(k)) then

            ! NB: f will be forced overcast by update_fSv

            ! Want to use a symmetric triangle (PL=0.5, SM=Smean)
            ! of half width hDe, if it gives endpoints in [1,Smax]
            call triangleSymPlus(1.,g%Smax(k),hDe,g%Smean(k),PLmin,PLmax,PLx,De)

            ! set parameters
            g%PL(k) = PLx
            g%alpha(k) = min(De/DeMaxScalar(g%Smean(k),g%Smax(k),PLx,1.-PLx),1.)

          ! Partially cloudy case
          else

            ! setup
            qc_qs = qcmean(k) / g%qs(k)
            qv_qs = qvmean(k) / g%qs(k)
            ftolc = 1. - ftol
            ! f mod routine wont work if f == 0 or 1
            f(k) = min(max(f(k),0.001),0.999) ! for dble prec

            ! iterate until in range
            fInRange = .false. ! by default
            do nmods = 1, maxfmods+1

              ! (1) cR consistent with qv_qs, qc_qs, f,
              !     i.e., along an Rm contour.
              qch_qs = qc_qs / f(k)
              qv0h_qs = (qv_qs - f(k)) / (1. - f(k))
              cR = qch_qs / (1. - qv0h_qs)
              cR = min(max(cR, 1.e-6), 1.e+6) ! for dble prec

              ! (2) Plim type flo, fhi
              ! [for reference, open limits for PLmin=0, PLmax=1
              ! flo = Finv(1./cR); fhi = 1. - Finv(cR)]
              cRlo = (1. - PLmax) / PLmax
              cRhi = (1. - PLmin) / PLmin
              if (cR <= cRlo) then
                ! both lower and upper in SM <= 1
                flo = efcontScalar(1-PLmin,1./cR)
                fhi = efcontScalar(1-PLmax,1./cR)
              else if (cR >= cRhi) then
                ! both lower and upper in SM >= 1
                flo = 1. - efcontScalar(PLmin,cR)
                fhi = 1. - efcontScalar(PLmax,cR)
              else
                ! lower in SM <= 1, upper in SM >= 1
                flo = efcontScalar(1-PLmin,1./cR)
                fhi = 1. - efcontScalar(PLmax,cR)
              end if

              ! (3) first time through set goal
              if (nmods == 1) then
                if (f(k) < flo) then
                  ! initial point is below range, so goal is to slide up
                  ! (higher f) and to left (lower cR) on the Rm contour
                  ! until we get to the correct f = flo. But during this
                  ! slide, the local flo decreases. Hence, the initial
                  ! f is the lower bound, and the initial flo is the
                  ! upper bound.
                  low = .true.; lof = f(k); hif = flo
                else if (f(k) > fhi) then
                  ! initial point is above range, so goal is to slide down
                  ! (lower f) and to right (higher cR) on the Rm contour
                  ! until we get to the correct f = fhi. But during this
                  ! slide, the local fhi increases. Hence, the initial
                  ! f is the upper bound, and the initial fhi is the
                  ! lower bound.
                  low = .false.; hif = f(k); lof = fhi
                else
                  ! in range first time through
                  fInRange = .true.
                  exit
                end if

              ! (4) subsequent iterations
              else
                ! (a) Not first time through, so goal already set.
                ! Put goal slightly (fraction ftol) inside the
                ! fhi - flo range to allow for convergence error
                if (low) then
                  fgoal = ftolc*flo+ftol*fhi
                else
                  fgoal = ftol*flo+ftolc*fhi
                end if
                ! (b) convergence check
                if (abs(f(k)-fgoal) < ftol*(fhi-flo)) then
                  fInRange = .true.
                  exit
                ! (c) New bounds on goal.
                ! again use monotonic properties of flo|fhi and Rm
                ! New bounds must not be outside existing bounds,
                ! but this is automatic for f, since it bisects
                ! the old bounds.
                else
                  if (f(k) < fgoal) then
                    lof = f(k); hif = min(hif,fgoal)
                  else  ! f(k) > fgoal
                    hif = f(k); lof = max(lof,fgoal)
                  end if
                end if
              end if

              ! (5) try a bisection of current range
              f(k) = 0.5 * (lof + hif)

            end do ! nmods

            ! get a solution
            validSolution = .false.
            if (fInRange) then
              call diagTriangle(f(k), qv0h_qs, qch_qs, SL, SM, SH, De, PLx, PHx, cRx, rc)
              validSolution = (rc == 0)
            end if
            if (validSolution) &
              validSolution = (SL >= 0. .and. SH <= g%Smax(k))

            ! Process a valid solution
            if (validSolution) then
              g%PL(k) = PLx
              g%alpha(k) = min(De/DeMaxScalar(g%Smean(k),g%Smax(k),PLx,PHx),1.)

            ! No solution or a bad solution are all covered
            ! by this "fail-all" branch, which attempts a
            ! [0,Smax] triangle, with modifications as
            ! needed to honor Smean(k).
            else
              call triangleSLSH(0.,g%Smax(k),g%Smean(k),PLmin,PLmax,PLx,De)
              g%PL(k) = PLx
              g%alpha(k) = min(De/DeMaxScalar(g%Smean(k),g%Smax(k),PLx,1.-PLx),1.)
            end if

          end if  ! partially cloudy

        end do  ! layer(k)

      end if  ! ICpdfEst

    else
      call myError('unrecognized margin')
    end if

!   ! update qvmean and f from S, Smax, PL, and al
!   ! really only need in get_
!   call update_fSv(g%nk, g%Smean, g%Smax, g%PL, g%alpha, f, Svmean)
!   qvmean = Svmean * g%qs

    ! make a copy of current variable params
    allocate (g%xSmean(g%nkv)); g%xSmean = g%Smean(klim:nk)
    allocate (g%xPL   (g%nkv)); g%xPL    = g%PL   (klim:nk) 
    allocate (g%xalpha(g%nkv)); g%xalpha = g%alpha(klim:nk)
     
    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    ! Inter-layer Correlations
    ! ~~~~~~~~~~~~~~~~~~~~~~~~

    ! establish correlation matrix based on corrDelp
    ! a. scaled layer separation
    allocate (g%CorM(g%nkv,g%nkv))
    g%CorM = spread(pMidv,2,g%nkv)
    g%CorM = abs(g%CorM-transpose(g%CorM))/delpCorr
    ! b. correlation model
   !g%CorM = exp(-g%CorM)                          ! FOAR (Exponential)
    g%CorM = (1+g%CorM)*exp(-g%CorM)               ! SOAR
   !g%CorM = (1+g%CorM*(1+g%CorM/3.))*exp(-g%CorM) ! TOAR
    ! Riishojgaard flow-dependent correlation
    if (delSmeanRiish > 0.) then
      allocate (Sva(g%nkv,g%nkv))
      Sva = spread(g%Smean(klim:nk),2,g%nkv)
      Sva = (Sva - transpose(Sva)) / delSmeanRiish
      g%CorM  = g%CorM * exp(-0.5 * Sva**2)
      deallocate (Sva)
    end if
    ! c. Upper Cholesky of correlation
    allocate (HcorM(g%nkv,g%nkv))
    HcorM = g%CorM  ! spotrf overwrites input
    call spotrf('U', g%nkv, HcorM, g%nkv, rc)
    if (rc /= 0) call myError('CorM: not positive definite')
    do j=1,g%nkv-1
      do i=j+1,g%nkv
        HcorM(i,j)=0. ! set lower to zero
      end do
    end do
    call GaussianCopula_create (g%copula, g%nkv, HcorM)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Pre-calculate top (non-variable) MSE-type output
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Specific quantity (mass / mass of parcel),
    !   i.e., svap = rhov/rho, scon = rhoc/rho.
    ! For svap, see notes on rqs and cgTf2 earlier.
    ! For scon in cloudy air (where qv = qs),
    !   scon = (rhot/rhov-1)*rhov/rho
    !     = (qt/qv-1)*svap = (qt/qs-1)*svap
    !     = (S-1)*svap
    if (g%klm1 .ge. 1) then
      allocate (g%top_svap(g%klm1))
      allocate (g%top_scon(g%klm1))
      where (g%Smean(1:g%klm1) > 1.)
        ! cloudy
        g%top_svap = 1. / (g%cgTf2(1:g%klm1) + g%Smean(1:g%klm1))
        g%top_scon = (g%Smean(1:g%klm1) - 1.) * g%top_svap
      elsewhere
        ! clear
        where (g%Smean(1:g%klm1) > 0.)
          g%top_svap = 1. / (g%rqs(1:g%klm1) / g%Smean(1:g%klm1) + omreps)
        elsewhere
          g%top_svap = 0.
        end where
        g%top_scon = 0.
      end where
    end if

  end subroutine GridcolumnStatisticalModel_create

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine triangleSLSH (SL, SH, Smean, PLmin, PLmax, PL, De)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Attempt to use a [SL,SH] triangle than honors Smean.
  ! If necessary, make endpoint adjustments to honor Smean.

    real*4, intent(in)  :: SL, SH, Smean, PLmin, PLmax
    real*4, intent(out) :: PL, De

    real*4 :: SM

    ! argument validity
    if (.not.(SL < Smean .and. Smean < SH)) &
      call myError('require SL < Smean < SH')

    ! default width
    De = SH - SL

    ! A [SL,SH] triangle cannot give an Smean this small,
    ! so adjust SH downwards, keeping SL fixed, and PL at
    ! the lower bound PLmin, to honor Smean.
    if (Smean < SL + De/3.*(1.+PLmin)) then
      PL = PLmin; De = 3.*(Smean-SL)/(1.+PLmin)

    ! A [SL,SH] triangle cannot give an Smean this big,
    ! so adjust SL upwards, keeping SH fixed, and PL at
    ! the upper bound PLmax, to honor Smean.
    else if (Smean > SL + De/3.*(1.+PLmax)) then
      PL = PLmax; De = 3.*(SH-Smean)/(2.-PLmax)

    ! Here, endpoints can be fixed at SL and SH.
    ! and SM diagnosed honoring the Smean.
    else
      ! use default De = SH - SL
      SM = 3.*Smean-(SL+SH)
      PL = (SM-SL)/De
    end if

  end subroutine triangleSLSH

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine triangleSymPlus(S0, S1, hDe, Smean, PLmin, PLmax, PL, De)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Try a symmetric triangle (PL=0.5, SM=Smean) of half
  ! width hDe, if it gives endpoints in [S0,S1]. If
  ! necessary, make endpoint adjustments to honor Smean.

    real*4, intent(in)  :: S0, S1, hDe, Smean, PLmin, PLmax
    real*4, intent(out) :: PL, De

    ! locals
    real*4 :: Smean0, Smean1, TSmean0, TSmean1, fac

    ! argument validity
    if (.not.(S0 < Smean .and. Smean < S1)) &
      call myError('require S0 < Smean < S1')
    if (hDe <= 0) call myError('require hDe > 0')

    ! useful quantities
    Smean0 = Smean - S0
    Smean1 = S1 - Smean

    ! (1) both ends violated, try a [S0,S1] triangle
    if (Smean0 < hDe .and. Smean1 < hDe) then
      call triangleSLSH(S0, S1, Smean, PLmin, PLmax, PL, De)

    ! (2) only S0 violated, try a [S0,SM+hDe] triangle
    else if (Smean0 < hDe) then

      TSmean0 = 3.*Smean0; fac = 1.+PLmin
      if (TSmean0 >= hDe*fac/(1.-PLmin)) then
        ! An [S0,SM+hDe] triangle can honor Smean
        PL = (TSmean0-hDe)/(TSmean0+hDe)
        De = hDe/(1.-PL)
      else
        ! Need to shift dn SH, while keeping SL=S0
        ! and PL=PLmin, in order to honor Smean.
        PL = PLmin; De = TSmean0/fac
      end if

    ! (3) only S1 violated, try a [SM-hDe,S1] triangle
    else if (Smean1 < hDe) then

      ! note: PHmin = 1 - PLmax
      TSmean1 = 3.*Smean1; fac = 2.-PLmax
      if (TSmean1 >= hDe*fac/PLmax) then
        ! An [SM-hDe,S1] triangle can honor Smean
        PL = 1. - (TSmean1-hDe)/(TSmean1+hDe)
        De = hDe/PL
      else
        ! Need to shift up SL, while keeping SH=S1
        ! and PL=PLmax, in order to honor Smean.
        PL = PLmax; De = TSmean1/fac
      end if

    ! (4) neither violated, use a sym(hDe) triangle
    else

      PL = 0.5; De = hDe * 2.

    end if

  end subroutine triangleSymPlus


  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GridcolumnStatisticalModel_destroy (g)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    type (GridcolumnStatisticalModel), intent(out) :: g
    g%nk = -1
    if (allocated (g%delp  )) deallocate (g%delp  )
    if (allocated (g%pTope )) deallocate (g%pTope )
    if (allocated (g%pBote )) deallocate (g%pBote )
    if (allocated (g%pMid  )) deallocate (g%pMid  )
    if (allocated (g%T     )) deallocate (g%T     )
    if (allocated (g%fice  )) deallocate (g%fice  )
    if (allocated (g%fliq  )) deallocate (g%fliq  )
    if (allocated (g%qs    )) deallocate (g%qs    )
    if (allocated (g%Smax  )) deallocate (g%Smax  )
    if (allocated (g%rqs   )) deallocate (g%rqs   )
    if (allocated (g%cgTf2 )) deallocate (g%cgTf2 )
    if (allocated (g%Smean )) deallocate (g%Smean )
    if (allocated (g%PL    )) deallocate (g%PL    )
    if (allocated (g%alpha )) deallocate (g%alpha )
    if (allocated (g%xSmean)) deallocate (g%xSmean)
    if (allocated (g%xPL   )) deallocate (g%xPL   )
    if (allocated (g%xalpha)) deallocate (g%xalpha)
    if (allocated (g%CorM  )) deallocate (g%CorM  )
    call GaussianCopula_destroy(g%copula)
    if (allocated (g%top_svap)) deallocate (g%top_svap)
    if (allocated (g%top_scon)) deallocate (g%top_scon)
  end subroutine GridcolumnStatisticalModel_destroy


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! generate MoistureSubcolumnEnsemble of nsamples from the GSM
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GridcolumnStatisticalModel_generate ( &
    g, nSubcols, mse)

    use MKL_VSL_TYPE
    use MKL_VSL
    use mod_GaussianCopula
    use mod_triangle
    implicit none

    character(*), parameter :: loc = 'GSM_generate' // ': '

    type (GridcolumnStatisticalModel), intent(in) :: g 
    type (MoistureSubcolumnEnsemble), intent(out) :: mse
    integer, intent(in) :: nSubcols
     
    integer :: k, rc
    type (VSL_STREAM_STATE) :: brngStream
    integer, parameter :: seed = 7777777
    integer, parameter :: brng = VSL_BRNG_MCG31
    real*4, dimension(nSubcols, g%nkv) :: rank, S
    real*4 :: U(nsubcols)
    real*4, dimension (g%nkv) :: &
      xPH, xDe, xSmode, xDeL, xDeH

    ! make space for MSE
    call MoistureSubcolumnEnsemble_destroy (mse)
    mse%nSubcolumns = nSubcols
    allocate (mse%svap (nSubcols, g%nkv))
    allocate (mse%scon (nSubcols, g%nkv))

    ! generate the ranks
    if (overlap == OVERLAP_GCOP) then
      rc = vslnewstream (brngStream, brng, seed)
      if (rc /= VSL_STATUS_OK) &
        call myError (loc//'vslnewstream problem')
      call GaussianCopula_generate ( &
        g%copula, nSubcols, brngStream, rank)
      rc = vsldeletestream (brngStream)
      if (rc /= VSL_STATUS_OK) &
        call myError (loc//'vsldeletestream problem')
    else
      call myError(loc//'unknown overlap')
    end if

    ! generation of S
    if (margin == MARGIN_TRIANGLE) then
      ! augmentation of xparams for subcol generation
      xPH = 1. - g%xPL
      xDe = g%xalpha * DeMax(g%nkv, g%xSmean, g%Smax(g%klim:g%nk), g%xPL, xPH)
      xDeL = xDe * g%xPL; xDeH = xDe * xPH
      xSmode = g%xSmean + (xDeL - xDeH) / 3.
      ! sample from the distribution
      do k = 1, g%nkv
        U = rank(:,k)
        where (U <= g%xPL(k))
          U = xSmode(k) - xDeL(k) * (1. - sqrt(  U    / g%xPL(k)))
        elsewhere
          U = xSmode(k) + xDeH(k) * (1. - sqrt((1.-U) /   xPH(k)))
        end where
        S(:,k) = U
      end do
    else
      call myError (loc//'unknown margin')
    end if

    ! mse specific quantity (mass / mass of parcel)
    do k = 1, g%nkv
      where (S(:,k) > 1.)                               !SPEED
        ! cloudy
        mse%svap(:,k) = 1. / (g%cgTf2(g%klm1+k) + S(:,k))
        mse%scon(:,k) = (S(:,k) - 1.) * mse%svap(:,k)
      elsewhere
        ! clear
        mse%svap(:,k) = 1. / (g%rqs(g%klm1+k) / S(:,k) + omreps)
        mse%scon(:,k) = 0.
      end where
    end do

  end subroutine GridcolumnStatisticalModel_generate


  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MoistureSubcolumnEnsemble_destroy (m)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    type (MoistureSubcolumnEnsemble), intent(out) :: m
    m%nSubcolumns = -1
    if (allocated (m%svap)) deallocate (m%svap)
    if (allocated (m%scon)) deallocate (m%scon)
  end subroutine MoistureSubcolumnEnsemble_destroy


end module mod_StaticGridcolumnStatisticalModel
