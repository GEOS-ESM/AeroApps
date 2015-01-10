! file: mod_triangle.f90
! triangle PDF functions

module mod_triangle

  use mod_utils
  implicit none

  private
  public pdf, cdf, diagTriangle, DeMax, DeMaxScalar, efcont, efcontScalar, update_fSv

contains

  subroutine InvalidSoln (context, variable, value, problem)
    character(*), intent(in) :: context, variable, problem
    real*4, intent(in) :: value
    character(255) :: err
    write(err,'(a,": ",a," (",g,") ",a)') &
      context, variable, value, problem
    ! call myError(err)
    print *, 'InvalidSolution: ', trim(err)
  end subroutine InvalidSoln

  ! PDF 
  ! a worker function:
  !   all argument error and consistency
  !   checking must be done outside
  ! expect:
  !   SL < SM < SH
  !   and DeL, DeH, De > 0
  !   and SL = SM - DeL,
  !       SH = SM + DeH,
  !       De = DeL + DeH.
  function pdf (n, s, SL, SM, SH, DeL, DeH, De) result(p)
    integer, intent(in) :: n
    real*4, intent(in), dimension(n) :: &
      s, SL, SM, SH, DeL, DeH, De
    real*4 :: p(n)
    p = 0.
    where (SL < s .and. s <= SM)
      p = 2. / De * (s - SL) / DeL
    elsewhere (SM < s .and. s < SH)
      p = 2. / De * (SH - s) / DeH
    endwhere
  end function pdf

  ! CDF
  ! a worker function:
  !   all argument error and consistency
  !   checking must be done outside
  ! expect:
  !   SL < SM < SH
  !   and DeL, DeH > 0
  !   and SL = SM - DeL,
  !       SH = SM + DeH,
  !   and DeL = PL * De,
  !       DeH = PH * De,
  !   and PL in (0,1)
  !   and PH = 1 - PL.
  function cdf (n, s, SL, SM, SH, DeL, DeH, PL, PH) result (r)
    integer, intent(in) :: n
    real*4, intent(in), dimension(n) :: &
      s, SL, SM, SH, DeL, DeH, PL, PH
    real*4 :: r(n)
    r = 0.
    where (SL < s .and. s <= SM)
      r = PL * ((s - SL) / DeL) ** 2
    elsewhere (SM < s .and. s < SH)
      r = 1. - PH * ((SH - s) / DeH) ** 2
    elsewhere (SH <= s) 
      r = 1. 
    endwhere
  end function cdf

  ! integral of S over clear part of gridbox
  ! a worker function: 
  !   all argument error and consistency
  !   checking must be done outside
  ! expect:
  !   SL < SM < SH
  !   and DeL, DeH, De > 0
  !   and SL = SM - DeL,
  !       SH = SM + DeH,
  !       De = DeL + DeH.
  !   and fd = PS(1,SL,SM,SH,DeL,DeH,PL,PH)
  !     as above.
  function funIS0 (n, SL, SM, SH, DeL, DeH, De, fd) result (x)
    integer, intent(in) :: n
    real*4, intent(in), dimension(n) :: &
      SL, SM, SH, DeL, DeH, De, fd
    real*4 :: x(n), SM1(n)
    SM1 = SM - 1.; x = 0.
    where (SL < 1. .and. 1. <= SM)
      x = (SM1**2 * (1. - 2./3. * SM1 / DeL) - DeL**2 / 3.) / De
    elsewhere (SM < 1. .and. 1. < SH)
      x = (SM1**2 * (1. + 2./3. * SM1 / DeH) - DeL**2 / 3.) / De
    elsewhere (1. >= SH)
      x = (DeH - DeL) / 3.
    endwhere
    x = x + SM * fd
  end function funIS0

  ! solve for sigma
  ! returns an sg solution or sg < 0 if none
  ! valid solutions in the range [0,1)
  ! << scalar version >>
  function sgSolve (ef, eR) result(sg)
    real*4, intent(in) :: ef, eR
    real*4 :: sg
    real*4 :: tau, D
    ! triangle.tex: ef = \f, eR = \R
    ! Note: not necessary to check for efmin explic-
    !   itly as D > 0 handles existence condition
    tau = ef*(3.-ef+(1.-ef)*eR)/2.
    D = tau**2 - ef
    if (D > 0. .and. ef <= 1./(1.+eR)) then
      sg = 1. - tau - sqrt(D)
    else
      sg = -1.
    end if
  end function sgSolve

  ! diagnose triangle for partial cloudy case
  ! only knows the three input parameters: f,qv0h_qs,qch_qs
  ! returns the PDF parameters: SL, SM, SH, De, PL, PH, plus cR and return code, rc
  ! Does NOT check for endpoint [0, Smax] violations
  ! << scalar version >>
  subroutine diagTriangle (f, qv0h_qs, qch_qs, SL, SM, SH, De, PL, PH, cR, rc)
    real*4, intent(in)  :: f, qv0h_qs, qch_qs
    real*4, intent(out) ::  SL, SM, SH, De, PL, PH, cR
    integer, intent(out) :: rc
    real*4 :: fd, Q, sgL, sgH, DeL, DeH
    if (.not. (0 < f .and. f < 1)) call myError('require f in (0,1)')
    if (qch_qs <= 0.)              call myError('require qch_qs > 0')
    if (qv0h_qs >= 1.)             call myError('require qv0h_qs < 1')
    rc = -1  ! default failure
    fd = 1. - f; cR = qch_qs / (1. - qv0h_qs)
    Q = fd * (qv0h_qs - 1.) + f * qch_qs
    if (f >= cR/(1.+cR)) then
      ! (a) SM >= 1 branch
      sgL = sgSolve(fd,cR)
      if (.not.(0. <= sgL .and. sgL < 1.)) then
        call InvalidSoln('SM >= 1','sgL',sgL,'not in [0,1)')
        return
      end if
      PL = fd / (1. - sgL) ** 2
      if (.not.(0. < PL .and. PL < 1.)) then
        call InvalidSoln('SM >= 1','PL',PL,'not in (0,1)')
        return
      end if
      PH = 1. - PL
      De = Q / (PL * sgL + (PH - PL) / 3.)
      if (.not.(De > 0.)) then
        call InvalidSoln('SM >= 1','De',De,'not > 0')
        return
      end if
      DeL = PL * De; SM = 1. + sgL * DeL
      SL = SM - DeL; SH = SL + De
    else 
      ! (b) SM <= 1 branch
      sgH = sgSolve(f,1./cR)
      if (.not.(0. <= sgH .and. sgH < 1.)) then
        call InvalidSoln('SM <= 1','sgH',sgH,'not in [0,1)')
        return
      end if
      PH = f / (1. - sgH) ** 2
      if (.not.(0. < PH .and. PH < 1.)) then
        call InvalidSoln('SM <= 1','PH',PH,'not in (0,1)')
        return
      end if
      PL = 1. - PH
      De = Q / ((PH - PL) / 3. - PH * sgH)
      if (.not.(De > 0.)) then
        call InvalidSoln('SM <= 1','De',De,'not > 0')
        return
      end if
      DeH = PH * De; SM = 1. - sgH * DeH
      SH = SM + DeH; SL = SH - De
    end if
    rc = 0  ! success
  end subroutine diagTriangle

  ! the maximum width permitted without either
  !   (a) the lower tail crossing 0
  !   (b) the upper tail crossing Smax
  ! a worker function:
  !   all argument error and consistency
  !   checking must be done outside
  ! expect:
  !   PL in (0,1), PH = 1 - PL
  !   Smean in [0,Smax], Smax >= 0
  function DeMax (n, Smean, Smax, PL, PH) result (De)
    integer, intent(in) :: n
    real*4, intent(in), dimension(n) :: Smean, Smax, PL, PH
    real*4 :: De(n), Dee(n,2)
    Dee(:,1) = Smean / (1. + PL) 
    Dee(:,2) = (Smax - Smean) / (1. + PH)
    De = 3. * minval(transpose(Dee), dim=1)
  end function DeMax

  function DeMaxScalar (Smean, Smax, PL, PH) result (De)
    real*4, intent(in) :: Smean, Smax, PL, PH
    real*4 :: De, De1(1)
    De1 = DeMax (1, [Smean], [Smax], [PL], [PH])
    De = De1(1)
  end function DeMaxScalar

  ! F function
  function F(y) result (x)
    real*4, intent(in) :: y(:)
    real*4 :: x(size(y))
    if (any(y <= 0. .or. y >= 1.)) &
      call myError('require y in (0,1)')
    x = sqrt(y)
    x = 2. / (x * (1. + x)) - 1.
  end function F

  ! inverse F function
  function Finv(x) result (y)
    real*4, intent(in) :: x(:)
    real*4 :: y(size(x))
    if (any(x <= 0)) call myError('require x > 0')
    y = 0.25 * (sqrt(1.+8./(x+1.))-1.)**2
  end function Finv

  ! ef contour function
  function efcont(n, p, eR) result(ef)
    integer, intent(in) :: n
    real*4, intent(in) :: p(n), eR(n)
    real*4, dimension(n) :: ef, &
      eRmin, kpp, fac, psi, X, c1, c2
    if (any(p <= 0. .or. p >= 1.)) then
      print *, 'p: ', p
      call myError('require p in (0,1)')
    end if
    eRmin = (1.-p)/p
    if (any(eR < eRmin)) then
      print *, 'eR, eRmin: ', eR, eRmin
      call myError('require eR >= eRmin')
    end if
    kpp = sqrt(p); kpp = kpp + 1./kpp
    fac = 1.+eR/3.; psi = fac/(1.+eR)
    X = -0.5*kpp/fac/sqrt(psi)
    c1 = cos(acos(-X)/3.)
    c2 = cos(acos( X)/3.)
    ef = 4.*psi*(c1-c2)**2
  end function efcont

  function efcontScalar(p,eR) result(ef)
    real*4, intent(in) :: p, eR
    real*4 :: ef, ef1(1)
    ef1 = efcont(1, [p],[eR])
    ef = ef1(1)
  end function efcontScalar

  ! update f and Sv given:
  !   S, Smax, PL, al
  ! a worker function:
  !   all argument error and consistency
  !   checking must be done outside
  ! expect:
  !   PL in (0,1), al in [0,1]
  !   S, Smax as needed for DeMax
  subroutine update_fSv (n, S, Smax, PL, al, f, Sv)
    integer, intent(in) :: n
    real*4, intent(in), dimension(n) :: S, Smax, PL, al
    real*4, intent(out), dimension(n) :: f, Sv
    logical :: ix(n)
    real*4, dimension(count(al > 0)) :: &
      PLx, PH, Sx, De, DeL, DeH, SM, SL, SH, fd, fx
    integer :: m 
    ! regular function cases
    ix = (al > 0.)
    if (any(ix)) then
      m = count(ix)
      PLx = pack(PL,mask=ix)
      PH = 1. - PLx
      Sx = pack(S,mask=ix)
      De = pack(al,mask=ix) * DeMax(m,Sx,pack(Smax,mask=ix),PLx,PH)
      DeL = PLx * De
      DeH = PH  * De
      SM = Sx - (DeH - DeL) / 3.
      SL = SM - DeL
      SH = SM + DeH
      Sx = 1.; fd = cdf(m,Sx,SL,SM,SH,DeL,DeH,PLx,PH)
      fx = 1. - fd; f = unpack(fx,ix,f)
      Sv = unpack(fx+funIS0(m,SL,SM,SH,DeL,DeH,De,fd),ix,Sv)
    end if
    ! delta function cases
    where (.not. ix)
      where (S <= 1.)
        f = 0.; Sv = S  ! clear
      elsewhere
        f = 1.; Sv = 1. ! overcast
      endwhere
    endwhere
  end subroutine update_fSv

end module mod_triangle
