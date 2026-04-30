module Icosahedron_Mod
    use m_Spherical_Partition_Simplified
    implicit none
    private 
    public icosahedron_regions
    public find_buddies 
    
    ! Keep ONLY the partitioner saved, as it takes real time to initialize
    type(Spherical_Partition), save :: partition

contains

    subroutine icosahedron_regions(nobs, xobs, yobs, zobs, kr, maxreg)
        implicit none
        integer, intent(in) :: nobs
        real*8, intent(in) :: xobs(nobs), yobs(nobs), zobs(nobs)
        integer, intent(out) :: kr(nobs)
        integer, intent(out) :: maxreg
        
        ! Initialize Partitioner
        call Initialize(n_levels=0, partition=partition, compress=.true.)
        maxreg = NumberOfRegions(partition)
        call XYZ2REG(nobs, xobs, yobs, zobs, kr, partition)  
    end subroutine icosahedron_regions

    subroutine find_buddies(nobs, n_susp, ki_susp, kr_susp, &
                           xobs, yobs, zobs, lats, lons, lev, omf, varF, varO, qcx, &
                           ls_h, ls_v, search_rad, tau_buddy, single_level, nbuddy_max, &
                           iregbeg, ireglen, maxreg, seplim, reaccept)
        implicit none

        ! Input parameters
        integer, intent(in) :: nobs, n_susp, nbuddy_max, maxreg
        integer, intent(in) :: ki_susp(n_susp), kr_susp(n_susp)
        real*8, intent(in) :: xobs(nobs), yobs(nobs), zobs(nobs)
        real*8, intent(in) :: lats(nobs), lons(nobs), lev(nobs)
        real*8, intent(in) :: omf(nobs), varF(nobs), varO(nobs)
        integer, intent(in) :: qcx(nobs)
        real*8, intent(in) :: ls_h, ls_v, search_rad, tau_buddy, seplim
        logical, intent(in) :: single_level
        integer, intent(in) :: iregbeg(maxreg), ireglen(maxreg)

        ! Output parameters
        logical, intent(out) :: reaccept(nobs)

        ! Local variables
        integer :: i, j, nbuddy, ireg, ibeg, iend, kis, krs, lvs, is, ib, ibb
        real*8 :: exponent, scgain, dist2, z_dist2, accum_de1, accum_de2, accum_wgt, accum_var
        real*8 :: tol_rel = 1.0e-5  
        real*8, parameter :: radius_earth = 6371000.0  

        ! ---------------------------------------------------------
        ! PURE LOCAL ALLOCATABLE ARRAYS (Clean, safe, no memory leaks)
        ! ---------------------------------------------------------
        integer, allocatable :: buddy_indices(:,:)
        real*8, allocatable :: buddy_weights(:,:)
        real*8, allocatable :: mysep(:,:)
        logical, allocatable :: suspect(:)
        real*8, allocatable :: weighted_pred, var_scaling_factor, tol2, obs_dev_squared
        
        ! Temporary sorting arrays
        integer, allocatable :: temp_buddy_indices(:), indx(:)
        real*8, allocatable :: temp_buddy_weights(:)

        ! Allocate everything cleanly for this specific run
        allocate(buddy_indices(nobs, nbuddy_max))
        allocate(buddy_weights(nobs, nbuddy_max))
        allocate(mysep(n_susp, maxreg))
        allocate(suspect(nobs))
        allocate(temp_buddy_indices(nobs), temp_buddy_weights(nobs), indx(nobs))

        ! Initialize standard arrays
        suspect = .false.
        reaccept = .false.

        do is = 1, n_susp
            suspect(ki_susp(is)) = .true.
        end do

        if (.not. InitializedQ(partition)) then
            print*, "Error: partition not initialized"
            return
        end if

        ! Pre-calculate separation between regions
        do is = 1, n_susp
            krs = kr_susp(is)
            do ireg = 1, maxreg
               mysep(is,ireg) = Separation(GetRegion(ireg,partition), GetRegion(krs,partition), SEPANG_MIN)
            end do
        end do

        ! BUDDY CHECK LOOP
        do is = 1, n_susp  
            kis = ki_susp(is)   
            krs = kr_susp(is)   
            lvs = lev(kis)      

            nbuddy = 0
            do ireg = 1, maxreg
                if ((mysep(is, ireg) <= seplim) .and. (ireglen(ireg) > 0)) then  
                    ibeg = iregbeg(ireg)
                    iend = ibeg + ireglen(ireg) - 1
                    
                    do i = ibeg, iend  
                        if (single_level .and. (abs(lvs-lev(i)) > tol_rel*lvs)) cycle

                        if (.not. suspect(i) .and. qcx(i) == 0) then
                            dist2 = (radius_earth**2) * ((xobs(kis)-xobs(i))**2 + &
                                                        (yobs(kis)-yobs(i))**2 + &
                                                        (zobs(kis)-zobs(i))**2)
                            exponent = dist2 / ls_h**2

                            if (ls_v > tol_rel) then  
                                z_dist2 = (lev(kis) - lev(i))**2
                                exponent = exponent + z_dist2 / ls_v**2
                            end if

                            if (exponent < search_rad**2) then 
                                nbuddy = nbuddy + 1
                                temp_buddy_indices(nbuddy) = i
                                scgain = varF(i) / (varF(i) + varO(i))
                                temp_buddy_weights(nbuddy) = scgain * exp(-0.5*exponent)
                            end if  
                        end if  
                    end do  
                end if  
            end do  

            ! Process buddies if found
            if (nbuddy > 0) then
                do i = 1, nbuddy
                    indx(i) = i  
                end do

                ! Sort buddies by weight (descending)
                do i = 1, nbuddy-1
                    do j = i+1, nbuddy
                        if (temp_buddy_weights(indx(j)) > temp_buddy_weights(indx(i))) then
                            ibb = indx(i)
                            indx(i) = indx(j)
                            indx(j) = ibb
                        end if
                    end do
                end do

                nbuddy = min(nbuddy, nbuddy_max)

                do ib = 1, nbuddy
                    ibb = indx(ib)
                    buddy_indices(kis, ib) = temp_buddy_indices(ibb)
                    buddy_weights(kis, ib) = temp_buddy_weights(ibb)
                end do

                accum_de1 = 0.0d0
                accum_de2 = 0.0d0
                accum_wgt = 0.0d0
                accum_var = 0.0d0

                do ib = 1, nbuddy
                    i = buddy_indices(kis, ib)
                    accum_wgt = accum_wgt + buddy_weights(kis, ib)
                    accum_de1 = accum_de1 + buddy_weights(kis, ib) * omf(i)
                    accum_var = accum_var + (varF(i) + varO(i))
                    accum_de2 = accum_de2 + omf(i)**2
                end do

                if (accum_wgt > 0.0d0 .and. accum_var > 0.0d0) then
                    weighted_pred = (accum_de1 / accum_wgt)
                    var_scaling_factor = (accum_de2 / accum_var)
                    
                    tol2 = tau_buddy**2 * var_scaling_factor * (varF(kis) + varO(kis))
                    obs_dev_squared = (omf(kis) - weighted_pred)**2

                    if (obs_dev_squared < tol2) then
                        reaccept(kis) = .true.
                    end if
                end if
            end if
        end do  

        ! Clean up memory before returning to Python
        deallocate(buddy_indices, buddy_weights, mysep, suspect, temp_buddy_indices, temp_buddy_weights, indx)

    end subroutine find_buddies

end module Icosahedron_Mod
