#include "ErrorHandle.h"

module disort_ancillary

	use modis_numerical_module, only: linearinterpolation
	use geos5_io_module
		
	use hdf5

	implicit none

	real, parameter :: GEOS5_Ptop = 0.01 ! top pressure in GEOS5 0.01 hPa
	integer, parameter :: num_levels = 27 ! 35
	integer, parameter :: num_layers = num_levels - 1 
	integer, parameter :: nchan = 25
	integer, parameter :: nPol = 1
!	integer, parameter :: nMom = 17 ! 16-stream calculations
!	integer, parameter :: nMom = 33 ! 32-stream calculations
	integer, parameter :: nMom = 65 ! 64-stream calculations
	
	real :: disort_heights(num_levels) 

	character(len=200) :: rcfile 


! there is no longer a concept of cloud phase because we are going to mix the phase 
! on an as-needed basis. 
	type disort_type
	
		integer*1 :: surface_level
		integer*1 :: have_cloud	
		real, dimension(nchan) :: albedo_factor
	
		real, dimension(num_levels) :: pressure, temperature, heights, ozone, mix_ratio
		real :: Ts
		integer*1, dimension(num_layers) :: cloud_layers
		real*8, dimension(num_layers) :: cloud_tau_ice, cloud_re_ice, cloud_tau_liq, cloud_re_liq
		real*8, dimension(num_layers, nchan) :: aero_tau, aero_ssa, aero_g
		real*4, dimension(:,:,:), allocatable :: pMom
	

	
	end type disort_type

	
	type(disort_type), allocatable, dimension(:,:) :: granule_info
	
	integer*2, dimension(:,:,:), allocatable :: surface_albedo

	logical :: AEROSOL_yes 


contains


	subroutine obtain_geos5(grid_xsize, grid_ysize, grid_xstart, grid_ystart, &
				geos5_name, geos5_aerosol_name, model_info)

		use general_array_io
		use ErrorHandleMod
		use AOP_Mod

! arguments	
		character(*), intent(in) :: geos5_name, geos5_aerosol_name
		integer, intent(in) :: 	grid_xsize, grid_ysize, grid_xstart, grid_ystart
! output goes here
		type(disort_type), dimension(:,:), intent(inout) :: model_info

! local variables
		integer :: i,j,k,l, ii
		real, dimension(:,:,:), allocatable :: met_temp1, met_temp2, met_temp3, met_temp4
		integer, dimension(:,:), allocatable :: met_temp1_2d
		real, dimension(:,:), allocatable :: met_temp2_2d, met_temp3_2d
		real, dimension(:,:), allocatable :: Pmsl, W2m
		real, dimension(:,:,:), allocatable :: delP, cf, tau_ice, tau_liq
		real, dimension(:,:,:), allocatable :: qlls, qlan, qils, qian, iwp, lwp	
			
		real, dimension(:,:,:), allocatable :: cf_out
		real, dimension(:,:,:,:), allocatable :: aero_temp1, aero_temp2, aero_temp3
		integer :: file_id(1), err_code, var_id
		integer :: gfile_id(2), start(3), edge(3), stride(3)
		real :: mxr, o3, ti, tl
		integer :: file_levels
		real :: total_tau

		integer, parameter :: cloud_bands(nchan) = &
			(/ 9,10,3,4,1,2,11,12,13,5,6,16, 7, 8,14,15,17,18,19,20,21,22,23,24,25 /) 
		integer, parameter :: aerosol_bands(nchan) = &
			(/ 1, 2,3,4,5,6, 7, 8, 8,9,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23 /)

! for the GEOS-5 aerosols
		integer, parameter :: nq   = 15  ! number of tracers
  		! tracer names
  		character(len=*), parameter :: tracer(nq) = (/ &
  		  'BCPHILIC', 'BCPHOBIC', 'DU001', 'DU002', 'DU003', &
  		  'DU004', 'DU005', 'OCPHILIC', 'OCPHOBIC', 'SO4', &
  		  'SS001', 'SS002', & 'SS003', 'SS004', 'SS005' /)
  		integer :: n, km, nsc, nfc, nch, verbose 
		real, allocatable :: local_channels(:)         ! (nch)                  ! wavelengths (nm)

	  real, allocatable :: local_tau  (:,:,:)        ! (km,nfc,nsc)           ! aerosol optical depth
	  real, allocatable :: local_ssa  (:,:,:)        ! (km,nfc,nsc)           ! single scattering albedo
	  real, allocatable :: local_g    (:,:,:)        ! (km,nfc,nsc)           ! asymmetry factor
	  real, allocatable :: local_pmom (:,:,:,:,:)    ! (km,nfc,nsc,nMom,nPol) ! scattering phase matrix
		real, allocatable :: rh(:,:,:)           ! (km,nfc,nsc)           ! relative humidity
	  	real, allocatable :: qm(:,:,:,:)         ! (km,nq,nfc,nsc)        ! (mixing ratio) * delp/g

		type (AOP) :: tables

! ******* for the SSA test ********

	real :: index_tau, ssa_use(5), taus550(5), ssa_test_table(5,5)
	real :: tt_data(25)


! ********* end SSA test ********



		__Iam__('obtain_geos5')


	tt_data = (/0.8804, 0.8836, 0.8866, 0.8880, 0.8887, &
		0.8644, 0.8693, 0.8748, 0.8798, 0.8807, &
		0.8450, 0.8515, 0.8597, 0.8688, 0.8698, &
		0.7954, 0.8045, 0.8182, 0.8371, 0.8385, &
		0.7183, 0.7021, 0.6850, 0.6801, 0.6761 /)

	ssa_test_table = reshape (tt_data, (/5,5/))
	
	print*, ssa_test_table (1, :)
	print*, ssa_test_table (2, :)


! read the 3D fields needed for DISORT work	
		allocate(met_temp1(num_levels, grid_xsize, grid_ysize))
		allocate(met_temp2(num_levels, grid_xsize, grid_ysize))
		allocate(met_temp3(num_levels, grid_xsize, grid_ysize))
		allocate(met_temp4(num_levels, grid_xsize, grid_ysize))

		call h5fopen_f(geos5_name, H5F_ACC_RDONLY_F, gfile_id(1), err_code)	

		call read_float_array_g5(gfile_id, 1, "O3", (/ 0, grid_xstart, grid_ystart, 0 /), (/1,1,1,1/), &
										(/ num_levels, grid_xsize, grid_ysize, 1 /), 3, met_temp1, err_code)
		
		call read_float_array_g5(gfile_id, 1, "T", (/ 0, grid_xstart, grid_ystart,0 /), (/1,1,1,1/), &
										(/ num_levels, grid_xsize, grid_ysize,1 /), 3, met_temp2, err_code)

		call read_float_array_g5(gfile_id, 1, "QV", (/ 0, grid_xstart, grid_ystart,0 /), (/1,1,1,1/), &
										(/ num_levels, grid_xsize, grid_ysize,1 /), 3, met_temp3, err_code)

		call read_float_array_g5(gfile_id, 1, "PE", (/ 0, grid_xstart, grid_ystart,0 /), (/1,1,1,1/), &
										(/ num_levels, grid_xsize, grid_ysize,1 /), 3, met_temp4, err_code)


		do k=1, num_levels
			do i=1, grid_xsize
				do j=1, grid_ysize
					o3 = met_temp1(k,i,j) 

					model_info(i,j)%ozone(k) = o3
					
					model_info(i,j)%temperature(k) = met_temp2(k,i,j)
					model_info(i,j)%mix_ratio(k) = met_temp3(k,i,j) / (1. - met_temp3(k,i,j)) 
					model_info(i,j)%pressure(k) = met_temp4(k,i,j) / 100. ! convert to mb
					model_info(i,j)%heights(k) = disort_heights(k)					

				end do
			end do
		end do

		deallocate(met_temp1, met_temp2, met_temp3, met_temp4)

! now read the k_surface SDS to get the surface level for each column 
! also read Ts SDS. 
		allocate (met_temp1_2d(grid_xsize, grid_ysize))
		allocate (met_temp2_2d(grid_xsize, grid_ysize))
		allocate (met_temp3_2d(grid_xsize, grid_ysize))
		
		call read_int_array_g5(gfile_id, 1, "K_SURFACE", (/ grid_xstart, grid_ystart,0 /), (/1,1,1/), &
										(/ grid_xsize, grid_ysize,1 /), 2, met_temp1_2d, err_code)
		call read_float_array_g5(gfile_id, 1, "TS", (/ grid_xstart, grid_ystart,0 /), (/1,1,1/), &
										(/ grid_xsize, grid_ysize, 1/), 2, met_temp2_2d, err_code)
		call read_float_array_g5(gfile_id, 1, "ZS", (/ grid_xstart, grid_ystart,0 /), (/1,1,1/), &
										(/ grid_xsize, grid_ysize, 1/), 2, met_temp3_2d, err_code)
				
		do j=1, grid_ysize
			do i=1, grid_xsize
				model_info(i,j)%Ts = met_temp2_2d(i,j)
				model_info(i,j)%surface_level = met_temp1_2d(i,j)
				model_info(i,j)%heights(model_info(i,j)%surface_level) = met_temp3_2d(i,j) / 1000. ! in km
			end do
		end do
		
		deallocate(met_temp1_2d, met_temp2_2d, met_temp3_2d)
! now it's time to ingest cloud information

		allocate(met_temp1(num_layers, grid_xsize, grid_ysize))
		allocate(met_temp2(num_layers, grid_xsize, grid_ysize))
		allocate(met_temp3(num_layers, grid_xsize, grid_ysize))
		allocate(met_temp4(num_layers, grid_xsize, grid_ysize))

		call read_float_array_g5(gfile_id, 1, "TAUI", (/ 0, grid_xstart, grid_ystart,0 /), (/1,1,1,1/), &
										(/ num_layers, grid_xsize, grid_ysize,1 /), 3, met_temp1, err_code)

		call read_float_array_g5(gfile_id, 1, "TAUL", (/ 0, grid_xstart, grid_ystart,0 /), (/1,1,1,1/), &
										(/ num_layers, grid_xsize, grid_ysize,1 /), 3, met_temp2, err_code)

		call read_float_array_g5(gfile_id, 1, "REI", (/ 0, grid_xstart, grid_ystart,0 /), (/1,1,1,1/), &
										(/ num_layers, grid_xsize, grid_ysize,1 /), 3, met_temp3, err_code)

		call read_float_array_g5(gfile_id, 1, "REL", (/ 0, grid_xstart, grid_ystart,0 /), (/1,1,1,1/), &
										(/ num_layers, grid_xsize, grid_ysize,1 /), 3, met_temp4, err_code)
		
		do j=1, grid_ysize
			do i=1, grid_xsize
				model_info(i,j)%have_cloud = 0
				model_info(i,j)%cloud_layers(1:num_layers) = 0
				model_info(i,j)%albedo_factor = 0
			end do
		end do
			
		do j=1, grid_ysize
			do i=1, grid_xsize

				total_tau = 0.
				model_info(i,j)%cloud_tau_ice = 0.
				model_info(i,j)%cloud_tau_liq = 0.

				do k=1, num_layers

					model_info(i,j)%cloud_tau_ice(k) = met_temp1(k,i,j)
					model_info(i,j)%cloud_tau_liq(k) = met_temp2(k,i,j)


					if (model_info(i,j)%cloud_tau_ice(k) > 0. .and. model_info(i,j)%cloud_tau_ice(k) < 1000. .or. &
						model_info(i,j)%cloud_tau_liq(k) > 0. .and. model_info(i,j)%cloud_tau_liq(k) < 1000. ) then

						model_info(i,j)%cloud_layers(k) = 1
						model_info(i,j)%have_cloud = 1

						total_tau = total_tau + model_info(i,j)%cloud_tau_ice(k) + model_info(i,j)%cloud_tau_liq(k)

						model_info(i,j)%cloud_re_ice(k) = 30. !met_temp3(k,i,j)	
						model_info(i,j)%cloud_re_liq(k) = 12. !met_temp4(k,i,j)	

					else

						model_info(i,j)%cloud_tau_ice(k) = 0.
						model_info(i,j)%cloud_tau_liq(k) = 0.

						model_info(i,j)%cloud_re_ice(k) = 0.	
						model_info(i,j)%cloud_re_liq(k) = 0.	

					endif
					
					
				end do

! **** CLOUD-FREE run only
!				model_info(i,j)%albedo_factor = 0.					
				model_info(i,j)%albedo_factor = total_tau 


			end do
		end do



!		print*, model_info(1,1)%cloud_tau_ice(:)
!		print*, model_info(1,1)%cloud_tau_liq(:)

!		print*, model_info(1,1)%cloud_re_ice(:)
!		print*, model_info(1,1)%cloud_re_liq(:)
	
!		print*, model_info(1,1)%cloud_layers


		call h5fclose_f(gfile_id(1), err_code)		
		deallocate(met_temp1, met_temp2, met_temp3, met_temp4)
	
! now read the aerosol information	
		if (AEROSOL_yes) then

	tables = AOP_create (rcfile,__RCM__)

! the creature calls for SDS indexing to start from 1, unlike regular HDF
	call AOP_PrepareFromMCS(tables, &
    		geos5_name, geos5_aerosol_name, nq, tracer, grid_ystart+1, grid_ystart+grid_ysize, &
			grid_xstart+1, grid_xstart + grid_xsize, &
    		nsc, nfc, nch, local_channels, km, rh, qm, __RCM__)

  ! reserve output space
  allocate( local_tau  (km,nfc,nsc)           ,__STATM__)
  allocate( local_ssa  (km,nfc,nsc)           ,__STATM__)
  allocate( local_g    (km,nfc,nsc)           ,__STATM__)
  allocate( local_pmom (km,nfc,nsc,nMom,nPol) ,__STATM__)

	do j=1, grid_ysize
		do i=1, grid_xsize
			allocate(model_info(i,j)%pMom(num_layers, nMom, nchan))
		end do
	end do

!	rh = 0.8 ! **** fix the relative humidity to 80% for these test


 ! get AOPs for different channels
	  verbose = 0
  	  do k = 1, nchan ! nch

!             if ( k < 9 .or. ich > 12 .and. ich < 17 ) then


    	print *, 'calculating for ', local_channels(aerosol_bands(k)), ' (nm) ...'
    	call AOP_calculate (tables, &
				      km, nsc*nfc, local_channels(aerosol_bands(k)), nMom, nPol, nq, tracer, verbose, &
     				qm, rh, local_tau, local_ssa, local_g, local_pmom, __RCM__)
      				
      ! now pull apart the information to where it's supposed to go				
		do j=1, grid_ysize
			do i=1, grid_xsize
				do l=1, num_layers
!				local_tau(l,i,j) = 0. ! AEROSOL-FREE RUN
					if (local_tau(l,i,j) > 0. .and. local_tau(l,i,j) < 50.) then 
						model_info(i,j)%aero_tau(l,cloud_bands(k)) = local_tau(l, i, j)
						model_info(i,j)%aero_ssa(l,cloud_bands(k)) = local_ssa(l, i, j)
						model_info(i,j)%aero_g(l,cloud_bands(k)) = local_g(l, i, j)
						model_info(i,j)%albedo_factor(k) = model_info(i,j)%albedo_factor(k) + &
																				model_info(i,j)%aero_tau(l,cloud_bands(k))

						do ii =1, nMom
							model_info(i,j)%pMom(l,ii,cloud_bands(k)) = local_pmom(l,i,j,ii,1) / (2*(ii-1) + 1)				
						end do
					else
						model_info(i,j)%aero_tau(l,cloud_bands(k)) = 0.
						model_info(i,j)%aero_ssa(l,cloud_bands(k)) = 0.
						model_info(i,j)%aero_g(l,cloud_bands(k)) = 0.
						model_info(i,j)%pMom(l,:,cloud_bands(k)) = 0.
					endif	

				end do
				
			end do
		end do

	  end do

!*********** TEST CONDITIONS ONLY ******************
#if 0
	 taus550 = (/ 0.25, 0.5, 1.0, 2.0, 3.0 /)

	 do j=1, grid_ysize
		do i=1, grid_xsize


			index_tau = model_info(i,j)%albedo_factor(4) ! total AOD550

	  		if (index_tau <= 0.25) then 
	  			ssa_use(:) = ssa_test_table(1, :)
	  		else if (index_tau >= 3.0) then 
	  			ssa_use(:) = ssa_test_table(5, :)
	  		else

		  		do l=2, 5
					if (index_tau < taus550(l)) exit				
				end do
		  
			   do k=1, 5
			  		ssa_use(k) = linearinterpolation ( (/taus550(l-1), taus550(l) /), &
		  												 (/ ssa_test_table(l-1, k), ssa_test_table(l,k) /), &
		  												 index_tau )
			  	end do
			endif
			
		  	do k=1, num_layers
		  		
		 		if (model_info(i,j)%aero_tau(k, 1) > 0.) then 
	  					model_info(i,j)%aero_ssa(k, 3) = ssa_use(1) !470
 				 		model_info(i,j)%aero_ssa(k, 4) = ssa_use(2) !550
	 			 		model_info(i,j)%aero_ssa(k, 1) = ssa_use(3) !660
	 				 	model_info(i,j)%aero_ssa(k, 2) = ssa_use(4) !860
	  				 	model_info(i,j)%aero_ssa(k, 8) = ssa_use(5) !2.1 (1.2 is listed twice)
	 			endif
		  		
			end do

		  
		end do
	end do
#endif


!************** END TEST ***************************



		
!	stop

	deallocate(local_tau, local_ssa, local_g, local_pmom, local_channels, rh, qm)

	call AOP_destroy (tables, __RCM__)

   endif ! end aerosol_yes


   do j=1, grid_ysize
     do i=1, grid_xsize
        do k=1, nchan           
			if (model_info(i,j)%albedo_factor(k) > 3.) then
 				model_info(i,j)%albedo_factor(k) = 1.
			else
				model_info(i,j)%albedo_factor(k) = model_info(i,j)%albedo_factor(k) / 3.
			endif
		  end do
		end do
	end do



	
	end subroutine obtain_geos5


	subroutine read_surface_albedo(file_id, phi, ocean_mask, start, edge)
	
		use general_array_io, only : read_int16_array, read_byte_array, read_float_array
	
		integer, intent(in) :: file_id(:), start(:), edge(:)
		real, dimension(:,:), intent(inout) :: phi
		integer*1, dimension(:,:), intent(inout) :: ocean_mask
	
		integer :: stride(2), err_code
		
		stride = 1
		

		call read_byte_array(file_id, "Ocean_Mask", start, stride, edge, ocean_mask, err_code)

		call read_int16_array(file_id, "Surface_Albedo_1", start, stride, edge, surface_albedo(1,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_2", start, stride, edge, surface_albedo(2,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_3", start, stride, edge, surface_albedo(3,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_4", start, stride, edge, surface_albedo(4,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_5", start, stride, edge, surface_albedo(5,:,:), err_code)
! this is because two halves of channel 5 use same albedo. Channel 5 is in two pieces because of Corr-K stuff
		surface_albedo(6,:,:) = surface_albedo(5,:,:)
		call read_int16_array(file_id, "Surface_Albedo_6", start, stride, edge, surface_albedo(7,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_7", start, stride, edge, surface_albedo(8,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_8", start, stride, edge, surface_albedo(9,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_9", start, stride, edge, surface_albedo(10,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_17", start, stride, edge, surface_albedo(11,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_18", start, stride, edge, surface_albedo(12,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_19", start, stride, edge, surface_albedo(13,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_20", start, stride, edge, surface_albedo(14,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_22", start, stride, edge, surface_albedo(15,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_26", start, stride, edge, surface_albedo(16,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_27", start, stride, edge, surface_albedo(17,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_28", start, stride, edge, surface_albedo(18,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_29", start, stride, edge, surface_albedo(19,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_31", start, stride, edge, surface_albedo(20,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_32", start, stride, edge, surface_albedo(21,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_33", start, stride, edge, surface_albedo(22,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_34", start, stride, edge, surface_albedo(23,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_35", start, stride, edge, surface_albedo(24,:,:), err_code)
		call read_int16_array(file_id, "Surface_Albedo_36", start, stride, edge, surface_albedo(25,:,:), err_code)
	
		call read_float_array(file_id, "Phi", start, stride, edge, phi, err_code)
	
	
	end subroutine read_surface_albedo


end module disort_ancillary
