module BuddyCheck_Mod
    use m_Spherical_Partition
    use m_Spherical_Triangle
    use omp_lib
    
    implicit none
    private 
    public icosahedron_regions
    public find_buddies
    public cleanup_buddy_check 
    
    ! Keep internal variables private
    type(Spherical_Partition), save :: partition
    
    ! Variables for buddy statistics
    integer, allocatable, save :: n_buddies(:)
    real, allocatable, save :: weighted_avg(:), weighted_var(:)
    integer, allocatable, save :: buddy_indices(:,:)
    real, allocatable, save :: buddy_weights(:,:)
    logical, allocatable, save :: suspect(:)
    
contains

    subroutine icosahedron_regions(nobs, xobs, yobs, zobs, kr, maxreg)
        implicit none
        
        ! Arguments:
        ! ----------
        integer, intent(in) :: nobs
        real*8, intent(in) :: xobs(nobs)
        real*8, intent(in) :: yobs(nobs)
        real*8, intent(in) :: zobs(nobs)
        integer, intent(out) :: kr(nobs)
        integer, intent(out) :: maxreg
        
        integer:: n_levels
        
        n_levels = 0
        
        !Initialize Partitioner
        !----------------------
        call Initialize(n_levels=0, partition=partition, compress=.true.)
        maxreg = NumberOfRegions(partition)
        
        ! Compute index region from (x,y,z) after bin obs
        !-------------------------
        call XYZ2REG(nobs, xobs, yobs, zobs, kr, partition)  ! bin obs into regions
        
        ! Initialize arrays for buddy check
        if (.not. allocated(suspect)) allocate(suspect(nobs))
        suspect = .false.
        
        if (.not. allocated(n_buddies)) allocate(n_buddies(nobs))
        if (.not. allocated(weighted_avg)) allocate(weighted_avg(nobs))
        if (.not. allocated(weighted_var)) allocate(weighted_var(nobs))
        
    end subroutine icosahedron_regions
    
    subroutine find_buddies(nobs, n_susp, ki_susp, kr_susp, &
                           xobs, yobs, zobs, lev, omf, varF, varO, &
                           ls_h, ls_v, search_rad, single_level, nbuddy_max, &
                           iregbeg, ireglen, maxreg, seplim, &
                           reaccept)
        implicit none
        
        ! Input parameters
        integer, intent(in) :: nobs
        integer, intent(in) :: n_susp
        integer, intent(in) :: ki_susp(n_susp)
        integer, intent(in) :: kr_susp(n_susp)
        real, intent(in) :: xobs(nobs)
        real, intent(in) :: yobs(nobs)
        real, intent(in) :: zobs(nobs)
        real, intent(in) :: lev(nobs)
        real, intent(in) :: omf(nobs)
        real, intent(in) :: varF(nobs)
        real, intent(in) :: varO(nobs)
        real, intent(in) :: ls_h
        real, intent(in) :: ls_v
        real, intent(in) :: search_rad
        real, intent(in) :: seplim
        logical, intent(in) :: single_level
        integer, intent(in) :: nbuddy_max
        integer, intent(in) :: iregbeg(maxreg)
        integer, intent(in) :: ireglen(maxreg)
        integer, intent(in) :: maxreg
       
        ! Output parameters
        logical, intent(out) :: reaccept(nobs) 
       
        ! Local variables
        integer :: i, j, nbuddy, ireg, ibeg, iend, kis, krs, lvs
        real :: exponent, scgain, dist2, z_dist2
        real :: tol_rel = 1.0e-5  ! Small tolerance value
        real, parameter :: radius_earth = 6371000.0  ! Earth radius in meters
        integer :: is, ib, ibb
        real :: weight_sum, weighted_sum, diff, std_dev, accum_del, accum_de2, accum_wgt, accum_var
        
        ! Temporary arrays for all potential buddies (before sorting)
        integer :: temp_buddy_indices(nobs)
        real :: temp_buddy_weights(nobs)
        integer :: indx(nobs)
        
        ! Allocate buddy arrays if not already allocated
        if (.not. allocated(buddy_indices)) allocate(buddy_indices(nobs, nbuddy_max))
        if (.not. allocated(buddy_weights)) allocate(buddy_weights(nobs, nbuddy_max))
        
        print*, 'nobs, nsusp', shape(xobs), shape(ki_susp)


        ! Mark suspect observations
        do is = 1, n_susp
            suspect(ki_susp(is)) = .true.
        end do
        
        ! Initialize output arrays
        reaccept = .false.
        
        ! Now do buddy check for each suspect obs
        !$omp parallel do default(shared) &
        !$omp private(is, kis, krs, lvs, nbuddy, ireg, ibeg, iend, i, j, exponent, scgain, &
        !$omp         weight_sum, weighted_sum, diff, std_dev, temp_buddy_indices, temp_buddy_weights, &
        !$omp         indx, ib, ibb, accum_del, accum_de2, accum_wgt, accum_var, dist2, z_dist2)
        do is = 1, n_susp  ! for each suspect obs
            kis = ki_susp(is)   ! this suspect's index
            krs = kr_susp(is)   ! this suspect's region
            lvs = lev(kis)      ! this suspect's level
            
            ! Find all candidate buddies
            nbuddy = 0
            do ireg = 1, maxreg
                ! Check if regions are close enough
                if ((Separation(GetRegion(ireg, partition), GetRegion(krs, partition), SEPANG_MIN) <= seplim) &
                   .and. (ireglen(ireg) > 0)) then  ! nearby region with data
                    
                    ibeg = iregbeg(ireg)
                    iend = ibeg + ireglen(ireg) - 1
                    
                    do i = ibeg, iend  ! look within a region
                        ! Skip other levels if doing single level
                        if (single_level .and. (abs(lvs-lev(i)) > tol_rel*lvs)) cycle
                        
                        ! only if not suspect, not excluded:
                        if (.not. suspect(i)) then
                            ! Calculate horizontal distance squared
                            dist2 = (radius_earth**2) * ((xobs(kis)-xobs(i))**2 + &
                                                        (yobs(kis)-yobs(i))**2 + &
                                                        (zobs(kis)-zobs(i))**2)
                            
                            exponent = dist2 / ls_h**2
                            
                            ! Add vertical distance component if applicable
                            if (ls_v > tol_rel) then  ! upper-air data
             !                   print*, 'I am in the vertical'
                                z_dist2 = (lev(kis) - lev(i))**2
                                exponent = exponent + z_dist2 / ls_v**2
                            end if
                            
                            ! only if not too far (search_rad length scales):
                            if (exponent < search_rad**2) then ! found a candidate buddy:
                                nbuddy = nbuddy + 1
                                ! Store in temporary arrays
                                temp_buddy_indices(nbuddy) = i
                                    
                                ! associate weight with this candidate:
                                scgain = varF(i) / (varF(i) + varO(i))
                                temp_buddy_weights(nbuddy) = scgain * exp(-0.5*exponent)
                            end if  ! not too far
                        end if  ! not suspect
                    end do  ! within a region
                end if  ! nearby regions
            end do  ! over regions
            !$OMP MASTER
             if (mod(is, 10) == 0) then  ! Print every 100th iteration
                 print *, "Progress: iteration", is, "of", nbuddy
             end if
            !$OMP END MASTER 
             !   print*, 'region', krs, nbuddy, 'for supsect', is
            
            ! Process buddies if we found any
            if (nbuddy > 0) then
                ! Sort buddies by weight (descending order)
                do i = 1, nbuddy
                    indx(i) = i  ! Initialize index array
                end do
                
                ! Index sort
                do i = 1, nbuddy-1
                    do j = i+1, nbuddy
                        if (temp_buddy_weights(indx(j)) > temp_buddy_weights(indx(i))) then
                            ! Swap indices
                            ibb = indx(i)
                            indx(i) = indx(j)
                            indx(j) = ibb
                        end if
                    end do
                end do
                
                ! Limit number of buddies to nbuddy_max
                nbuddy = min(nbuddy, nbuddy_max)
                
                ! Store the best buddies in the output arrays
                do ib = 1, nbuddy
                    ibb = indx(ib)
                    buddy_indices(kis, ib) = temp_buddy_indices(ibb)
                    buddy_weights(kis, ib) = temp_buddy_weights(ibb)
                end do
                
                ! Calculate statistics using the best buddies
                accum_del = 0.0
                accum_de2 = 0.0
                accum_wgt = 0.0
                accum_var = 0.0
                
                do ib = 1, nbuddy
                    i = buddy_indices(kis, ib)
                    accum_wgt = accum_wgt + buddy_weights(kis, ib)
                    accum_del = accum_del + buddy_weights(kis, ib) * omf(i)
                    accum_var = accum_var + buddy_weights(kis, ib) * (varF(i) + varO(i))
                end do
                
                if (accum_wgt > tol_rel) then
                    ! Calculate weighted average
                    weighted_avg(kis) = accum_del / accum_wgt
                    
                    ! Calculate weighted variance
                    do ib = 1, nbuddy
                        i = buddy_indices(kis, ib)
                        diff = omf(i) - weighted_avg(kis)
                        accum_de2 = accum_de2 + buddy_weights(kis, ib) * diff * diff
                    end do
                    
                    weighted_var(kis) = accum_de2 / accum_wgt
                    
                    ! Add background error variance
                    weighted_var(kis) = weighted_var(kis) + accum_var / accum_wgt
                    
                    ! Decide whether to reaccept based on criteria
                    std_dev = sqrt(weighted_var(kis))
                    if (abs(omf(kis) - weighted_avg(kis)) < 3.0 * std_dev) then
                        reaccept(kis) = .true.
                    end if
                end if
            end if
            
            n_buddies(kis) = nbuddy
        end do  ! for each suspect obs
        !$omp end parallel do
        
    end subroutine find_buddies
    
    subroutine cleanup_buddy_check()
        ! Clean up allocated memory
        if (allocated(n_buddies)) deallocate(n_buddies)
        if (allocated(weighted_avg)) deallocate(weighted_avg)
        if (allocated(weighted_var)) deallocate(weighted_var)
        if (allocated(buddy_indices)) deallocate(buddy_indices)
        if (allocated(buddy_weights)) deallocate(buddy_weights)
        if (allocated(suspect)) deallocate(suspect)
        
        ! Clean up partition
        call Clean(partition)
    end subroutine cleanup_buddy_check
    
end module BuddyCheck_Mod
