
  subroutine getTauSfc ( tau, ch, freq, nch, T, q, o3, pm, km, &
                         ps, lat, zen, rc)
    use m_bt_drv
    implicit NONE   
    integer, intent(in)  :: km 
    integer, intent(in)  :: nch
    integer, intent(in)  :: ch(nch)
    real,    intent(in)  :: freq(nch)
    real,    intent(in)  :: T(km)
    real,    intent(in)  :: q(km)
    real,    intent(in)  :: o3(km)
    real,    intent(in)  :: pm(km)
    real,    intent(in)  :: ps
    real,    intent(in)  :: lat
    real,    intent(in)  :: zen  ! zenith angle of sounding
    
    real,    intent(out) :: tau(nch)
    integer, intent(out) :: rc

!
!   Uses bt_drv to compute surface transmittances.
!

    real*8  :: r8_freq(nch)
    real*8  :: r8_T(km)
    real*8  :: r8_q(km)
    real*8  :: r8_o3(km)
    real*8  :: r8_pm(km)
    real*8  :: r8_ps
    real*8  :: r8_lat
    real*8  :: r8_zen
    real*8  :: r8_tau(nch)

    integer :: mype=0, satnumber=1, iprt=0
    integer :: day=0, land=1, tovsrt(nch) 
    real*8  :: mw_zen=0.0, su_zen=0.0, emissmw=1.0, rho=0.0, tground=300.
    real*8  :: ssu_eff_cell(nch), secsun=1.0
    logical :: comp_trans = .true.
 
    real*8    :: btnew(nch)

    r8_freq = freq
    r8_T    = T
    r8_q    = q
    r8_o3   = o3
    r8_pm   = pm
    r8_ps   = ps
    r8_lat  = lat
    r8_zen  = zen
    tovsrt  = 1
    ssu_eff_cell = 1.0

    call bt_drv(mype, r8_T, r8_pm, r8_q, r8_pm, r8_o3, r8_pm,  r8_ps, r8_lat,   &
         ch, mw_zen, su_zen, r8_zen, day, secsun, emissmw, &
         land, rho, tground, btnew, rc, comp_trans, tovsrt,    &
         satnumber, iprt, r8_freq, ssu_eff_cell, tausfc=r8_tau)

    tau = r8_tau

  end subroutine getTauSfc
