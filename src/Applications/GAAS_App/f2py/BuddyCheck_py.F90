! Direct f2py interface for buddy check functionality

subroutine py_icosahedron_regions(nobs, xobs, yobs, zobs, kr, maxreg)
    use BuddyCheck_Mod
    implicit none
    
    ! Arguments:
    ! ----------
    integer, intent(in) :: nobs
    real*8, intent(in) :: xobs(nobs)
    real*8, intent(in) :: yobs(nobs)
    real*8, intent(in) :: zobs(nobs)
    integer, intent(out) :: kr(nobs)
    integer, intent(out) :: maxreg
    
    ! Call the actual implementation
    call icosahedron_regions(nobs, xobs, yobs, zobs, kr, maxreg)
    
end subroutine py_icosahedron_regions

subroutine py_find_buddies(nobs, n_susp, ki_susp, kr_susp, &
                          xobs, yobs, zobs, lats, lons,lev, omf, varF, varO, qcx,&
                          ls_h, ls_v, search_rad, tau_buddy, single_level, nbuddy_max, &
                          iregbeg, ireglen, maxreg, seplim, &
                          reaccept)
    use BuddyCheck_Mod
    implicit none
    
    ! Input parameters
    integer, intent(in) :: nobs
    integer, intent(in) :: n_susp
    integer, intent(in) :: ki_susp(n_susp)
    integer, intent(in) :: kr_susp(n_susp)
    real*8, intent(in) :: xobs(nobs)
    real*8, intent(in) :: yobs(nobs)
    real*8, intent(in) :: zobs(nobs)
    real*8, intent(in) :: lats(nobs)
    real*8, intent(in) :: lons(nobs)
    real*8, intent(in) :: lev(nobs)
    real*8, intent(in) :: omf(nobs)
    real*8, intent(in) :: varF(nobs)
    real*8, intent(in) :: varO(nobs)
    integer, intent(in) :: qcx(nobs)
    real*8, intent(in) :: ls_h
    real*8, intent(in) :: ls_v
    real*8, intent(in) :: search_rad
    real*8, intent(in) :: tau_buddy
    real*8, intent(in) :: seplim
    logical, intent(in) :: single_level
    integer, intent(in) :: nbuddy_max
    integer, intent(in) :: iregbeg(maxreg)
    integer, intent(in) :: ireglen(maxreg)
    integer, intent(in) :: maxreg
    ! Output parameters
    logical, intent(out) :: reaccept(nobs)
    
    ! Call the actual implementation
    call find_buddies(nobs, n_susp, ki_susp, kr_susp, &
                     xobs, yobs, zobs, lats, lons, lev, omf, varF, varO,qcx, & 
                     ls_h, ls_v, search_rad, tau_buddy,single_level, nbuddy_max, &
                     iregbeg, ireglen, maxreg, seplim, &
                     reaccept)
    
end subroutine py_find_buddies

subroutine py_cleanup_buddy_check()
    use BuddyCheck_Mod
    implicit none
    call cleanup_buddy_check()
end subroutine py_cleanup_buddy_check    
