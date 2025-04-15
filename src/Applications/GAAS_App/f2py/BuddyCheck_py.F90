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
                          xobs, yobs, zobs, lev, omf, varF, varO, &
                          ls_h, ls_v, search_rad, single_level, nbuddy_max, &
                          iregbeg, ireglen, maxreg, seplim, &
                          reaccept)
    use BuddyCheck_Mod
    implicit none
    
    ! Input parameters
    integer, intent(in) :: nobs, n_susp
    integer, intent(in) :: ki_susp(n_susp), kr_susp(n_susp)
    real, intent(in) :: xobs(nobs), yobs(nobs), zobs(nobs), lev(nobs)
    real, intent(in) :: omf(nobs), varF(nobs), varO(nobs)
    real, intent(in) :: ls_h, ls_v, search_rad, seplim
    logical, intent(in) :: single_level
    integer, intent(in) :: nbuddy_max
    integer, intent(in) :: iregbeg(maxreg), ireglen(maxreg)
    integer, intent(in) :: maxreg
    ! Output parameters
    logical, intent(out) :: reaccept(nobs)
    
    ! Call the actual implementation
    call find_buddies(nobs, n_susp, ki_susp, kr_susp, &
                     xobs, yobs, zobs, lev, omf, varF, varO, &
                     ls_h, ls_v, search_rad, single_level, nbuddy_max, &
                     iregbeg, ireglen, maxreg, seplim, &
                     reaccept)
    
end subroutine py_find_buddies

subroutine py_cleanup_buddy_check()
    use BuddyCheck_Mod
    implicit none
    call cleanup_buddy_check()
end subroutine py_cleanup_buddy_check    
