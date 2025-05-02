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
    real, allocatable, save :: weighted_pred(:), var_scaling_factor(:)
    real, allocatable, save :: tol2(:), obs_dev_squared(:)
    integer, allocatable, save :: buddy_indices(:,:)
    real, allocatable, save :: buddy_weights(:,:)
    real, allocatable, save :: omf_buddy(:,:)

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
        real :: tau_buddy = 1   ! buddy tol param
        real, parameter :: radius_earth = 6371000.0  ! Earth radius in meters
        integer :: is, ib, ibb, n_reacc
        real ::  accum_de1, accum_de2, accum_wgt, accum_var
        
        ! Temporary arrays for all potential buddies (before sorting)
        integer :: temp_buddy_indices(nobs)
        real :: temp_buddy_weights(nobs)
        integer :: indx(nobs)
        
        ! Allocate buddy arrays if not already allocated
        if (.not. allocated(buddy_indices)) allocate(buddy_indices(nobs, nbuddy_max))
        if (.not. allocated(buddy_weights)) allocate(buddy_weights(nobs, nbuddy_max))
        if (.not. allocated(omf_buddy)) allocate(omf_buddy(nobs, nbuddy_max))
        
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
                        if (.not. suspect(i)) then
                            ! Calculate horizontal distance squared
                            dist2 = (radius_earth**2) * ((xobs(kis)-xobs(i))**2 + &
                                                        (yobs(kis)-yobs(i))**2 + &
                                                        (zobs(kis)-zobs(i))**2)
                            
                            exponent = dist2 / ls_h**2
                             !$OMP MASTER
                                if (mod(is, 10000) == 0) then  ! Print
                                 print *, "exponent for suspect:", is, "and buddy index", i, "value=", exponent
                                 end if
                                 !$OMP END MASTER 

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

                                 !$OMP MASTER
                                if (mod(is, 10000) == 0) then  ! Print
                                 print *, "scgain  for suspect:", is, "of", nbuddy, "value=", scgain
                                 end if
                                 !$OMP END MASTER 

                                temp_buddy_weights(nbuddy) = scgain * exp(-0.5*exponent)
                                !$OMP MASTER
                                if (mod(is, 10000) == 0) then  ! Print
                                 print *, "Poids for suspect:", is, "of", nbuddy, "value=", temp_buddy_weights(nbuddy)
                                 end if
                                !$OMP END MASTER  


                            end if  ! not too far
                        end if  ! not suspect
                    end do  ! within a region
                end if  ! nearby regions
            end do  ! over regions
            !$OMP MASTER
             if (mod(is, 1000) == 0) then  ! Print
                 print *, "Progress: iteration", is, "of", nbuddy
             end if
            !$OMP END MASTER 
             !   print*, 'region', krs, nbuddy, 'for supsect', is
            
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
                
                ! Limit number of buddies to nbuddy_max
                nbuddy = min(nbuddy, nbuddy_max)
                
                ! Store the best buddies in the output arrays
                do ib = 1, nbuddy
                    ibb = indx(ib)
                    buddy_indices(kis, ib) = temp_buddy_indices(ibb)
                    buddy_weights(kis, ib) = temp_buddy_weights(ibb)
                    omf_buddy(kis, ib) = omf(ibb)
                end do
                 
                
                  !$OMP MASTER
                  if (mod(is, 100) == 0) then  ! Print
                  print *, "buddy weight 10:", is, "of", kis, "value=", buddy_weights(kis, 1:10), 'omf bud',omf_buddy(kis,1:10)
                  end if
                  !$OMP END MASTER 

                ! Calculate statistics using the best buddies
                accum_de1 = 0.0
                accum_de2 = 0.0
                accum_wgt = 0.0
                accum_var = 0.0
                
                do ib = 1, nbuddy
                    i = buddy_indices(kis, ib)
                    accum_wgt = accum_wgt + buddy_weights(kis, ib)
                    accum_de1 = accum_de1 + buddy_weights(kis, ib) * omf(i)
                    accum_var = accum_var + buddy_weights(kis, ib) * (varF(i) + varO(i))
                    accum_de2 = accum_de2 + (omf(i))**2
                end do
                
                    
                    ! Calculate weighted prediction (expected value for a suspect) calculate from its buddies based
                    ! on obs minus forecast value
                    
                    weighted_pred(kis) = accum_de1 / accum_wgt
                    !$OMP MASTER
                    if (mod(is, 100) == 0) then  ! Print
                       print *, "Progress2: weighted pred for ind", is, kis, "value of", weighted_pred(kis)
                    end if
                    !$OMP END MASTER 
                    ! variance scaling factor used to adjust the tolerance for accepting/rejectig obs 
                    !(derived from the buddies var)
                    var_scaling_factor(kis) = accum_de2 / accum_var
                    !$OMP MASTER
                    if (mod(is, 100) == 0) then  ! Print
                       print *, "Progress3: var for scale factor for ind", is, kis, "value of", var_scaling_factor(kis)
                    end if
                    !$OMP END MASTER
                    
                    ! Decide whether to reaccept based on criteria-> tolerance squared
                    tol2(kis) = tau_buddy**2 * var_scaling_factor(kis)  * (varF(kis) + varO(kis))
                     !$OMP MASTER
                    if (mod(is, 100) == 0) then  ! Print
                       print *, "Progress4: tol function of tau_buddy", is, kis, "value of", tol2(kis)
                    end if
                    !$OMP END MASTER

                    obs_dev_squared(kis) = (omf(kis) - weighted_pred(kis))**2 

                     !$OMP MASTER
                    if (mod(is, 100) == 0) then  ! Print
                       print *, "Progress45: obs dev squared", is, kis, "value of", omf(kis), weighted_pred(kis), obs_dev_squared(kis)
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
        if (allocated(suspect)) deallocate(suspect)
        
        ! Clean up partition
        call Clean(partition)
    end subroutine cleanup_buddy_check
    
end module BuddyCheck_Mod
