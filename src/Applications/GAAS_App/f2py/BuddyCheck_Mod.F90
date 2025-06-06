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
    real*8, allocatable, save :: weighted_pred(:), var_scaling_factor(:)
    real*8, allocatable, save :: tol2(:), obs_dev_squared(:)
    integer, allocatable, save :: buddy_indices(:,:)
    real*8, allocatable, save :: buddy_weights(:,:)
    real*8, allocatable, save :: omf_buddy(:,:)
    real*8, allocatable, save :: buddy_dist2(:,:)
    real*8, allocatable, save :: buddy_distz2(:,:)
    real*8, allocatable, save :: buddy_lats(:,:)
    real*8, allocatable, save :: buddy_lons(:,:)

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
        
    end subroutine icosahedron_regions
    
    subroutine find_buddies(nobs, n_susp, ki_susp, kr_susp, &
                           xobs, yobs, zobs, lats, lons,lev, omf, varF, varO, qcx,&
                           ls_h, ls_v, search_rad,tau_buddy, single_level, nbuddy_max, &
                           iregbeg, ireglen, maxreg, seplim, &
                           reaccept)
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
       
        ! Local variables
        integer :: i, j, nbuddy, ireg, ibeg, iend, kis, krs, lvs
        real*8 :: exponent, scgain, dist2, z_dist2
        real*8 :: tol_rel = 1.0e-5  ! Small tolerance value
        real*8, parameter :: radius_earth = 6371000.0  ! Earth radius in meters
        integer :: is, ib, ibb, n_reacc
        real*8 ::  accum_de1, accum_de2, accum_wgt, accum_var
        
        ! Temporary arrays for all potential buddies (before sorting)
        integer:: temp_buddy_indices(nobs)
        real*8 :: temp_buddy_weights(nobs)
        real*8 :: temp_dist2(nobs)
        real*8 :: temp_distz2(nobs)
        integer :: indx(nobs)
      !  real*8 :: min_weight, max_weight, min_omf, max_omf
        
        ! Allocate buddy arrays if not already allocated
        if (.not. allocated(buddy_indices)) allocate(buddy_indices(nobs, nbuddy_max))
        if (.not. allocated(buddy_weights)) allocate(buddy_weights(nobs, nbuddy_max))
        if (.not. allocated(omf_buddy)) allocate(omf_buddy(nobs, nbuddy_max))
        if (.not. allocated(buddy_dist2)) allocate(buddy_dist2(nobs, nbuddy_max))
        if (.not. allocated(buddy_distz2)) allocate(buddy_distz2(nobs, nbuddy_max))
        if (.not. allocated(buddy_lats)) allocate(buddy_lats(nobs, nbuddy_max))
        if (.not. allocated(buddy_lons)) allocate(buddy_lons(nobs, nbuddy_max))
        
        
        print*, 'nobs, nsusp', shape(xobs), shape(ki_susp)

        ! Initialize arrays for buddy check
        if (.not. allocated(suspect)) allocate(suspect(nobs))
        suspect = .false.
        
        if (.not. allocated(weighted_pred)) allocate(weighted_pred(nobs))
        if (.not. allocated(var_scaling_factor)) allocate(var_scaling_factor(nobs))
        if (.not. allocated(tol2)) allocate(tol2(nobs))
        if (.not. allocated(obs_dev_squared)) allocate(obs_dev_squared(nobs))
        ! Mark suspect observations
        do is = 1, n_susp
            suspect(ki_susp(is)) = .true.
        end do
        
        ! Initialize output arrays
        reaccept = .false.
        n_reacc = 0 
        ! Now do buddy check for each suspect obs
        !$omp parallel do default(shared) &
        !$omp private(is, kis, krs, lvs, nbuddy, ireg, ibeg, iend, i, j, exponent, scgain, &
        !$omp         weighted_pred, var_scaling_factor, temp_buddy_indices, temp_buddy_weights, &
        !$omp         indx, ib, ibb, accum_de1, accum_de2, accum_wgt, accum_var, obs_dev_squared, tol2, dist2, z_dist2), &
        !$omp reduction(+:n_reacc)
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
                        if (.not. suspect(i) .and. qcx(i) ==0) then
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
                                temp_dist2(nbuddy) = dist2
                                temp_distz2(nbuddy) = z_dist2    
                                ! associate weight with this candidate:
                                scgain = varF(i) / (varF(i) + varO(i))
                                temp_buddy_weights(nbuddy) = scgain * exp(-0.5*exponent)
                            end if  ! not too far
                        end if  ! not suspect
                    end do  ! within a region
                end if  ! nearby regions
            end do  ! over regions
            
            ! Process buddies if we found any
            if (nbuddy > 0) then
                ! Sort buddies by weight (descending order)
                do i = 1, nbuddy
                indx(i) = i  ! Initialize index array for old array type omf, varF, varO to get the right one
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
                            
                nbuddy = min(nbuddy, nbuddy_max)
                
                ! Store the best buddies in the output arrays
                do ib = 1, nbuddy
                    ibb = indx(ib)
                    i = buddy_indices(kis, ib)
                    buddy_indices(kis, ib) = temp_buddy_indices(ibb)
                    buddy_weights(kis, ib) = temp_buddy_weights(ibb)
                !    buddy_dist2(kis, ib) = temp_dist2(ibb)
                !    buddy_distz2(kis, ib) = temp_distz2(ibb)
                    omf_buddy(kis,ib) = omf(i)
                    buddy_lats(kis, ib) = lats(i)
                    buddy_lons(kis, ib) = lons(i)
                enddo 

           !     !$OMP MASTER
            !        if (mod(is, 1) == 0) then  ! Print
               !        print*, 'buddy indices', is, kis, buddy_indices(kis,1:10)
               !        print*, 'buddy weights', is, kis, buddy_weights(kis,1:10)
               !        print*, 'buddy dist2', is, kis, buddy_dist2(kis,1:20)
               !        print*, 'buddy end dist2', is, kis, buddy_dist2(kis,80:100)
             !          print*, 'buddy omf', is, kis, omf_buddy(kis,1:10)
               !        print*, 'buddy omf test', is, kis, omf_buddy_(kis,1:10)
              !         print*, 'buddy lats', is, kis, buddy_lats(kis,1:20)
              !         print*, 'buddy lons', is, kis, buddy_lons(kis,1:20) 
               !        print*, 'buddy end lats', is, kis, buddy_lats(kis,80:100)
               !        print*, 'buddy end lons', is, kis, buddy_lons(kis,80:100)
              !      endif
              !   !$OMP END MASTER   


                ! Calculate statistics using the best buddies
                accum_de1 = 0.0
                accum_de2 = 0.0
                accum_wgt = 0.0
                accum_var = 0.0
                
                do ib = 1, nbuddy
                    i = buddy_indices(kis, ib)
                    accum_wgt = accum_wgt + buddy_weights(kis, ib)
                    accum_de1 = accum_de1 + buddy_weights(kis, ib) * omf(i)
                    accum_var = accum_var + (varF(i) + varO(i))
                    accum_de2 = accum_de2 + omf(i)**2
                end do
                
                    
                    ! Calculate weighted prediction (expected value for a suspect) calculate from its buddies based
                    ! on obs minus forecast value
                    
                    !$OMP MASTER
                    if (mod(is, 100) == 0) then  ! Print
                       print *, "LAT LON for ind", is, kis, "value of", lats(kis), lons(kis)
                    end if
                    !$OMP END MASTER
                    weighted_pred(kis) = (accum_de1 / accum_wgt)
                    !$OMP MASTER
                    if (mod(is, 100) == 0) then  ! Print
                       print *, "Progress2: weighted pred for ind", is, kis, "value of", weighted_pred(kis)
                    end if
                    !$OMP END MASTER
                    ! variance scaling factor used to adjust the tolerance for accepting/rejectig obs 
                    !(derived from the buddies var)
                    var_scaling_factor(kis) = (accum_de2 / accum_var)
                     !$OMP MASTER
                    if (mod(is, 100) == 0) then  ! Print
                       print *, "Progress3: var for scale factor for ind", is, kis, "value of", var_scaling_factor(kis)
                    end if
                    !$OMP END MASTER
                    ! Decide whether to reaccept based on criteria-> tolerance squared
                    tol2(kis) = tau_buddy**2 * var_scaling_factor(kis)  * (varF(kis) + varO(kis))

                    obs_dev_squared(kis) = (omf(kis) - weighted_pred(kis))**2 
                    !$OMP MASTER
                    if (mod(is, 100) == 0) then  ! Print
                       print *, "Progress5: obs dev squared", is, kis, "value of omf susp", omf(kis)
                   !    print *, "Progress5: obs dev squared", is, kis, "value of obs_dev_squared", obs_dev_squared(kis)
                   !    print *, "Progress5: obs dev squared compared to", is, kis, "value of tolerance", tol2(kis)
                    end if
                    !$OMP END MASTER
                    if (obs_dev_squared(kis) < tol2(kis)) then
                        reaccept(kis) = .true.
                        n_reacc = n_reacc +1
                    end if
            end if
        end do  ! for each suspect obs

        !$omp end parallel do
        
    end subroutine find_buddies
    
    subroutine cleanup_buddy_check()
        ! Clean up allocated memory
        if (allocated(weighted_pred)) deallocate(weighted_pred)
        if (allocated(var_scaling_factor)) deallocate(var_scaling_factor)
        if (allocated(tol2)) deallocate(tol2)
        if (allocated(obs_dev_squared)) deallocate(obs_dev_squared)
        if (allocated(buddy_indices)) deallocate(buddy_indices)
        if (allocated(buddy_weights)) deallocate(buddy_weights)
        if (allocated(omf_buddy)) deallocate(omf_buddy)
        if (allocated(buddy_dist2)) deallocate(buddy_dist2)
        if (allocated(buddy_distz2)) deallocate(buddy_distz2)
        if (allocated(buddy_lats)) deallocate(buddy_lats)
        if (allocated(buddy_lons)) deallocate(buddy_lons)
        if (allocated(suspect)) deallocate(suspect)
        
        ! Clean up partition
        call Clean(partition)
    end subroutine cleanup_buddy_check
    
end module BuddyCheck_Mod
