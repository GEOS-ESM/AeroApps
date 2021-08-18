module modis_simulator

	use dis_global_pars
	use dis_rayleigh
	use DISORT_solver
	use disort_ancillary
	
!	use ifport
!	USE RDI1MACH_f90, ONLY : R1MACH, D1MACH, I1MACH
	
	implicit none
	
! this is the main driver that contains all the inputs
! this is the interface between model and DISORT
! this code will be replaced when integration to model occurs
! its only function is to fill in arrays as needed and pass them along to the
! actual simulator. 


! for now set everything up for one point. Test one point, does it make sense? 
! once the point calculation looks okay, then set up for entire granule. 


	integer, parameter :: numOfChannels =  25
	integer, parameter :: numberLayers = 26 !34
	integer, parameter :: numberLevels = 27 !35
	integer, parameter :: maxKs = 30
!	integer, parameter :: numberStreams = 16 !32 !64
	integer, parameter :: numberStreams = 64 !32 !64
	integer, parameter :: channelNum(numOfChannels) = &
			(/ 1,2,3,4,5,5,6,7,8,9,17,18,19,20,22,26,27,28,29,31,32,33,34,35,36 /) 
	real, parameter :: waveLengths(numOfChannels) = &
			(/ 0.65, 0.86, 0.47, 0.55, 1.24, 1.24, 1.63, 2.13, 0.41, 0.44, 0.91, 0.94, 0.94, &
				3.7, 3.9, 1.38, 6.2, 7.3, 8.5, 11.0, 12.0, 13.2, 13.4, 13.8, 14.2 /)
	real, parameter :: solarFlux(numOfChannels) = &
			(/ 80.09, 33.90, 41.52, 37.52, 11.26, 11.19, 5.56, 4.67, 26.17, 18.39, &
			26.81, 8.47, 41.15, 2.02, 0.57, 10.83, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /) 

	integer, parameter :: wnum_low(numOfChannels) = &
			(/14925, 11416, 20877, 17699, 8052, 7893, 6053, 4640, 23810, 22321, 10870, 10627, 10363, 2604, &
			2507, 7194, 1450, 1338, 1149, 887, 815, 742, 725, 710, 695 /)
	integer, parameter :: wnum_hi(numOfChannels) = &
			(/ 16129, 11891, 21786, 18349, 8210, 8052, 6143, 4751, 24691, 22831, 11236, 10741, 10929, 2732, &
			2545, 7353, 1530, 1394, 1190, 928, 850, 758, 742, 725, 710/)


	! reflectance = 1, radiance = 0
	integer*1, parameter :: refl_or_rad(numOfChannels) = &
			(/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
	! 1 = do rayleigh, 0 = don't
	integer*1, parameter :: rayleigh(numOfChannels) = &
			(/ 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
	! 1 = do emission, 0 = don't
	integer*1, parameter :: planck(numOfChannels) = &
			(/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
	
	logical, parameter :: lambertian = .true.
	
	real, parameter :: d2r = 0.017453292519943295
	
	real, parameter :: f_delta_ice = 0.2206
	real, parameter :: f_delta_liquid = 0.0
	
	real, parameter :: pi = 3.141592653589793238


	! this is the Legendre coefficient storage
	! this data persists forever, like a constant, but has to be read from file. 
	character (len=100) :: ice_coef_name, liquid_coef_name

	integer, parameter :: num_ice_radii = 23
	integer :: ice_radii(num_ice_radii)
	real :: legendre_ssa_ice(numOfChannels-1, num_ice_radii)
	real :: legendre_ce_ice(numOfChannels-1, num_ice_radii)
	real :: legendre_tf_ice(numOfChannels-1, num_ice_radii)
	real :: legendre_lib_ice(numOfChannels-1, num_ice_radii, numberStreams)
	
	integer, parameter :: num_water_radii = 24
	integer :: water_radii(num_water_radii)
	real :: legendre_ssa_water(numOfChannels-1, num_water_radii)
	real :: legendre_ce_water(numOfChannels-1, num_water_radii)
	real :: legendre_tf_water(numOfChannels-1, num_water_radii)
	real :: legendre_lib_water(numOfChannels-1, num_water_radii, numberStreams)
	



! run variables, one set for each granule / run
					   
	real :: surfaceAlbedo(numOfChannels)			
	real :: solar_zenith, sensor_zenith(1), relative_azimuth(1)
	integer :: use_layers, use_levels
		
	integer :: num_cloud_layers
	integer*1 :: cloud_layers(numberLayers)
	real :: cloud_tau_layers(numberLayers)
	real :: cloud_effective_radius(numberLayers)
	integer*1 :: cloud_phase(numberLayers)

	real :: rayleigh_tau(numOfChannels, numberLayers)

	! this is where the output will go to	
	real :: output_data(numOfChannels)
	
	! this will get filled in for each layer only
	real, dimension(:,:), allocatable :: single_scatter_albedo,  extinction_efficiency
	real, dimension(:,:), allocatable :: truncation_factor
	real, dimension(:,:,:), allocatable :: legendre_coefficients 
	
	! this is for the correlated-K's storage
	! These correlated-K's should come from the correlated-K module
	! that was previously integrated into GEOS-5 processing
	real :: correlated_Ks_weights(numOfChannels, maxKs)
	real :: correlated_Ks_taus(numOfChannels, numberLayers, maxKs)
	integer*1 :: number_Ks(numOfChannels)	
	
	! this is the atmospheric profile information used by both DISORT and rayleigh code
	real :: Tatm(0:numberLayers)
		
	! this is the aerosol variable stuff
	integer :: myexp(numberStreams+1)
	real :: omegaaero(numberLayers, numOfChannels), gaero(numberLayers, numOfChannels), aerotau(numberLayers, numOfChannels)
	real :: pMom_aero(numberLayers, 0:numberStreams, numOfChannels)	
		
	private
	
	public :: init_legendre, init_disort, run_disort, cleanup_disort, ice_coef_name, liquid_coef_name
	
	
	
	
contains


	subroutine init_legendre
		
		include "hdf.f90"
		include "dffunc.f90"
	
		integer ::  file_id, var_id, err_code, start(3), stride(3), edge(3)
		integer :: i
			
    	myexp = (/(i,i=0,numberStreams)/)

	
!		ice_coef_name = "IceLegendreCoeffs.hdf"
!		liquid_coef_name = "WaterLegendreCoeffs.hdf"

		file_id = sfstart(ice_coef_name, DFACC_READ)
	
		start = 0
		stride = 1
		edge(1) = num_ice_radii

		var_id = sfselect(file_id, sfn2index(file_id, "Effective_Radius"))
		err_code = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), ice_radii)
		err_code = sfendacc(var_id)

		edge(1) = numOfChannels-1
		edge(2) = num_ice_radii
	
		var_id = sfselect(file_id, sfn2index(file_id, "Extinction_Coefficient"))
		err_code = sfrdata(var_id, start(1:2), stride(1:2), edge(1:2), legendre_ce_ice)
		err_code = sfendacc(var_id)

		var_id = sfselect(file_id, sfn2index(file_id, "Truncation_Factor"))
		err_code = sfrdata(var_id, start(1:2), stride(1:2), edge(1:2), legendre_tf_ice)
		err_code = sfendacc(var_id)

	
		var_id = sfselect(file_id, sfn2index(file_id, "Single_Scatter_Albedo"))
		err_code = sfrdata(var_id, start(1:2), stride(1:2), edge(1:2), legendre_ssa_ice)
		err_code = sfendacc(var_id)

		edge(3) = numberStreams
	
		var_id = sfselect(file_id, sfn2index(file_id, "Legendre_Coefficients"))
		err_code = sfrdata(var_id, start, stride, edge, legendre_lib_ice)
		err_code = sfendacc(var_id)
	
	
		err_code = sfend(file_id)
	
	
		file_id = sfstart(liquid_coef_name, DFACC_READ)
	
		start = 0
		stride = 1
		edge(1) = num_water_radii

		var_id = sfselect(file_id, sfn2index(file_id, "Effective_Radius"))
		err_code = sfrdata(var_id, start(1:1), stride(1:1), edge(1:1), water_radii)

		err_code = sfendacc(var_id)

		edge(1) = numOfChannels-1
		edge(2) = num_water_radii
	
		var_id = sfselect(file_id, sfn2index(file_id, "Extinction_Coefficient"))
		err_code = sfrdata(var_id, start(1:2), stride(1:2), edge(1:2), legendre_ce_water)
		err_code = sfendacc(var_id)
	
		var_id = sfselect(file_id, sfn2index(file_id, "Truncation_Factor"))
		err_code = sfrdata(var_id, start(1:2), stride(1:2), edge(1:2), legendre_tf_water)
		err_code = sfendacc(var_id)

		var_id = sfselect(file_id, sfn2index(file_id, "Single_Scatter_Albedo"))
		err_code = sfrdata(var_id, start(1:2), stride(1:2), edge(1:2), legendre_ssa_water)
		err_code = sfendacc(var_id)

		edge(3) = numberStreams
	
		var_id = sfselect(file_id, sfn2index(file_id, "Legendre_Coefficients"))
		err_code = sfrdata(var_id, start, stride, edge, legendre_lib_water)
		err_code = sfendacc(var_id)
		
		err_code = sfend(file_id)
		
	end subroutine init_legendre

	subroutine init_disort(modis_point, As_vector, sza, vza, phi)

! in init_disort we will read in the legendre coefficient business
! and initialize all the variables to whatever values they will be. 
! Somehow init_disort needs to be made aware of the existence of correlated-K 
! module. For now we will simply read the standard ASCII output 
! of the corr-k code. Again, most of this code will be replaced

		use corrK_profile, only : modis_channel
		
		type(disort_type), intent(in) :: modis_point
		integer*2, intent(in) :: As_vector(:)
		real*8, intent(in) :: sza, vza, phi
		
		real :: Patm(0:numberLayers), Zatm(0:numberLayers)



		integer :: re_index_liq, re_index_ice, layer_index
		real :: tau_liq, tau_ice
		integer :: i, j, k, istr
		character(len=100) dummy
		real :: ang
		real, parameter :: approx_limit = 0.3

	
	! 1 = cloud, 0 = no cloud
		cloud_layers = 0
		cloud_tau_layers = 0.
		cloud_effective_radius = 0.

		do i=1, numberLayers
			cloud_tau_layers(i) = modis_point%cloud_tau_ice(i) + modis_point%cloud_tau_liq(i)
			if (cloud_tau_layers(i) > 0.) then 
				cloud_layers(i) = 1

			! effective radius weighed by optical thickness to obtain combined phase
				cloud_effective_radius(i) = ( modis_point%cloud_re_ice(i) * modis_point%cloud_tau_ice(i) + &
										modis_point%cloud_re_liq(i) * modis_point%cloud_tau_liq(i) ) / &
										cloud_tau_layers(i)

				if (cloud_effective_radius(i) > 60.) cloud_effective_radius(i) = 60. 
			endif
		end do
		
		ang = sza*d2r
		if (ang < approx_limit) then 
			solar_zenith = 1. - ang*ang/2.
		else
			solar_zenith = cos(ang)
		endif

		ang = vza*d2r
		if (ang < approx_limit) then 
			sensor_zenith(1) = 1. - ang*ang/2.
		else
			sensor_zenith(1) = cos(ang)
		endif
		
!		solar_zenith = cos(sza * d2r)
!		sensor_zenith(1) = cos(vza * d2r)

		relative_azimuth(1) = phi * d2r
	
	! IR surface albedo is 0.015 for ocean surface
	! for land it will be whatever the emissivity map prescribes
		surfaceAlbedo = As_vector*0.001d0

		do i=1, numOfChannels
		
			number_Ks(i) = modis_channel(i)%num_weights
			correlated_Ks_weights(i,1:number_Ks(i)) = modis_channel(i)%weights(:)
			
			do j=1, number_Ks(i)
				do k=1, numberLayers
				
					correlated_Ks_taus(i,k,j) = modis_channel(i)%taus(j,k)
				
				end do
			end do
		end do 


		use_levels = modis_point%surface_level
		use_layers = use_levels-1

	! cloud layer information where we quiery the libraries and assemble the information
	! for each layer as needed. Only cloudy layers are stored here in order to save space
	
		num_cloud_layers = 0
		
		do i=1, numberLayers
			if (cloud_layers(i) == 1) num_cloud_layers = num_cloud_layers + 1
		end do
		
		allocate(single_scatter_albedo(num_cloud_layers, numOfChannels))
		allocate(extinction_efficiency(num_cloud_layers, numOfChannels))
		allocate(truncation_factor (num_cloud_layers, numOfChannels))
		allocate(legendre_coefficients(0:numberStreams, num_cloud_layers, numOfChannels))

		legendre_coefficients = 0.

		layer_index = 1
		do i=1, numberLayers
	
			if (cloud_layers(i) == 1) then 

					tau_liq = modis_point%cloud_tau_liq(i)
					tau_ice = modis_point%cloud_tau_ice(i)

					re_index_liq = locate_re(modis_point%cloud_re_liq(i), water_radii)
					re_index_ice = locate_re(modis_point%cloud_re_ice(i)*2., ice_radii)
					
					! prevent illegal memory access. OK, because we would be multiplying whatever we get from arrays by 0.
					if (tau_liq <= 0.) then
						tau_liq = 0.
						re_index_liq = 1
					endif
					if (tau_ice <= 0.) then 
						tau_ice = 0.
						re_index_ice = 1
					endif

					
!					print*, water_radii
					
!					print*, tau_liq, tau_ice, cloud_layers(i), i
					
!					print*, re_index_liq, re_index_ice, modis_point%cloud_re_liq(i), modis_point%cloud_re_ice(i)*2.

					single_scatter_albedo(layer_index,1:5) = (legendre_ssa_water(1:5, re_index_liq)*tau_liq + &
															  	legendre_ssa_ice(1:5, re_index_ice)*tau_ice) / (tau_liq+tau_ice)
					single_scatter_albedo(layer_index,6) = (legendre_ssa_water(5, re_index_liq)*tau_liq + &
																legendre_ssa_ice(5, re_index_ice)*tau_ice) / (tau_liq+tau_ice)
					single_scatter_albedo(layer_index,7:numOfChannels) = (legendre_ssa_water(6:numOfChannels-1, re_index_liq)*tau_liq + &
																legendre_ssa_ice(6:numOfChannels-1, re_index_ice)*tau_ice) / (tau_liq+tau_ice)

					truncation_factor(layer_index,1:5) = (legendre_tf_water(1:5, re_index_liq)*tau_liq + &
																legendre_tf_ice(1:5, re_index_ice)*tau_ice) / (tau_liq+tau_ice)
					truncation_factor(layer_index,6) = (legendre_tf_water(5, re_index_liq)*tau_liq + &
																legendre_tf_ice(5, re_index_ice)*tau_ice) / (tau_liq+tau_ice)
					truncation_factor(layer_index,7:numOfChannels) = (legendre_tf_water(6:numOfChannels-1, re_index_liq)*tau_liq + &
																		legendre_tf_ice(6:numOfChannels-1, re_index_ice)*tau_ice) / (tau_liq+tau_ice)

					
					extinction_efficiency(layer_index,1:5) = (legendre_ce_water(1:5, re_index_liq)*tau_liq + &
																legendre_ce_ice(1:5, re_index_ice)*tau_ice) / (tau_liq+tau_ice)
					extinction_efficiency(layer_index,6) = (legendre_ce_water(5, re_index_liq)*tau_liq + &
																legendre_ce_ice(5, re_index_ice)*tau_ice) / (tau_liq+tau_ice)
					extinction_efficiency(layer_index,7:numOfChannels) = (legendre_ce_water(6:numOfChannels-1, re_index_liq)*tau_liq + &
																legendre_ce_ice(6:numOfChannels-1, re_index_ice)*tau_ice) / (tau_liq+tau_ice)
					
					do j=1, numOfChannels
						do istr = 1, numberStreams
							if (j <= 5) &
								legendre_coefficients(istr-1, layer_index, j) = (legendre_lib_water(j, re_index_liq, istr)*tau_liq + &
																				legendre_lib_ice(j, re_index_ice, istr)*tau_ice)/ (tau_liq+tau_ice)
							if (j==6) &
								legendre_coefficients(istr-1, layer_index, j) = (legendre_lib_water(5, re_index_liq, istr)*tau_liq + &
																				legendre_lib_ice(5, re_index_ice, istr)*tau_ice)/ (tau_liq+tau_ice)
							if (j > 6) &
								legendre_coefficients(istr-1, layer_index, j) = (legendre_lib_water(j-1, re_index_liq, istr)*tau_liq + &
																				legendre_lib_ice(j-1, re_index_ice, istr)*tau_ice)/ (tau_liq+tau_ice)
						end do
					end do
					legendre_coefficients(0, layer_index, :) = 1.					
		   
				layer_index = layer_index + 1
		   
		   
			end if
	
	
	
		end do

! ingest the atmospheric profile. 
		Zatm = modis_point%heights
		Patm = modis_point%pressure
		Tatm = modis_point%temperature

! now that we've ingested everything, it's time to set basic DISORT run variables up

!		rayleigh_tau = 0.
		
		do i=1, numOfChannels
			if (rayleigh(i) == 1) then 
!				call rayleigh_tau_calc (numberLayers,i,Patm,Zatm,wnum_low(i)*1.,wnum_hi(i)*1.,rayleigh_tau(i,:))
! zero the thing out so we don't have residual taus from other points hanging around as altitude changes
				rayleigh_tau(i,:) = 0.
				call rayleigh_tau_calc (use_layers,i,Patm,Zatm,wnum_low(i)*1.,wnum_hi(i)*1.,rayleigh_tau(i,:))
			endif		
		end do

! aerosol property arrays

		omegaaero = modis_point%aero_ssa
		aerotau = modis_point%aero_tau
		gaero = modis_point%aero_g
		if (allocated (modis_point%pMom) ) then 
			pMom_aero = modis_point%pMom
		else 
			pMom_aero = 0.
		endif

	end subroutine init_disort
	
	
	
	
	subroutine run_disort(Tsfc, Ttoa, take_out)
	
		real*8, dimension(:), intent(inout) :: take_out
		real, intent(in) :: Tsfc, Ttoa
	
		integer :: ich, ik, ilyr, tauint, icld
		real :: use_tau(numberLayers), use_omega(numberLayers)
		real :: tauin, otau, omega
		real :: utau(numberLevels)
		character(len=127) :: header
		integer, parameter :: maxumu = 5
		integer, parameter :: maxphi = 5
		REAL, DIMENSION(numberLevels)    :: dfdt, flup, rfldir, rfldn, uavg
		real, dimension(maxumu, numberLevels, maxphi) :: uu, uus	
		REAL, DIMENSION(maxumu)    :: umudeg
		INTEGER ::   btype = 2,szindex(4)=(/6,12,13,29/), quiet = 3
		REAL ::      brho0=0.0, bk=0.0, btheta=0.0, bsigma=0.0, bt1=0.0, bt2=0.0, bscale=0.0
		REAL ::     bu10=3.0, bpcl=0.15, bsal=34.3
		real :: radiance, reflectance
		real, dimension(0:maxcmu, numberLayers) :: pmom
		real, dimension(0:maxumu) :: umu
		real, dimension(maxphi) :: phi
		real, dimension(maxumu, numberLevels)    :: u0u
		logical :: local_planck
		logical, dimension(7)      :: prnt = .false.
		real, dimension(0:maxcmu)  :: hl
		integer, dimension(numberLevels) :: iflg  

		logical :: usrtau = .true.
		logical :: usrang = .true.
		logical :: deltam = .true.
		logical :: onlyfl = .false.
	
		real :: fisot, phi0, temis
		integer :: ibcnd, numu, nphi
		real :: accur = 1.0e-6
		integer :: NLYR, NLEV, NSTR, NTAU, MUMU, MPHI, MCMU
		real :: WNL, WNH, SFX, SFA
		
		
		fisot = 0.
		ibcnd = 0
		numu = 1
		nphi = 1
		phi0 = 0.
		temis = 0.0
		
		
		! set up the DISORT's compute array
		utau = 0.
		u0u = 0.
		
		header = " "
		
		do ich = 1, numOfChannels
		
			radiance = 0. 
			reflectance = 0.
		
			output_data(ich) = 0.0
					
			use_tau = 0.
			use_omega = 0.
			
			if (planck(ich) == 0) then
				local_planck = .false.
			else
				local_planck = .true.
			endif


                        ! have to skip channels 8, 9, 17, 18 and 22  because we don't really need them and also because 
                        ! we can't seem to afford them with a 64-stream run

                       if ( ich < 9 .or. ich > 12 .and. ich < 17 ) then 


			do ik=1, number_Ks(ich)
	
	
				! this is the Rayleigh scattering part, not the aerosol scattering, just the atmosphere
				do ilyr = 1, use_layers !numberLayers
		
			
						use_tau(ilyr) = correlated_Ks_taus( ich, ilyr, ik)
						use_omega(ilyr) = 0. 
						pmom(:, ilyr) = 0.
!						pmom = 0.
						
						if (rayleigh(ich) == 1) then 

							use_omega(ilyr) = 1. / (1. + (use_tau(ilyr) / rayleigh_tau(ich, ilyr)))
							use_tau(ilyr) = use_tau(ilyr) + rayleigh_tau(ich, ilyr)
							pmom(:,ilyr) = 0.1
!							pmom = 0.1
						
						endif

				end do
!     if (ich == 1) then
!		print*, ich, "CORR_K: ", correlated_Ks_taus(ich, :, ik) 
!     endif
							
!				pmom = 0.
				icld = 1
				
				do ilyr = 1, use_layers !numberLayers

			
				    if (AEROSOL_yes .and. aerotau(ilyr, ich) > 0.001 ) then 

					   pmom(:, ilyr) = pMom_aero(ilyr, :, ich) 

! This is the H-G aerosol phase function code
!						pmom(0:numberStreams, ilyr) = gaero(ilyr, ich) ** myexp
!						pmom(1:numberStreams, ilyr) = (pmom(1:numberStreams, ilyr) - pmom(numberStreams, ilyr)) / &
!													   (1. - pmom(numberStreams, ilyr))
				
						! we assume 0.1 for the rayleigh legendre coefficient thing. 
						if (rayleigh(ich) == 1) then 
							pmom(1:numberStreams, ilyr) = ( use_tau(ilyr)*use_omega(ilyr)*0.1 + &
													  aerotau(ilyr, ich) * omegaaero(ilyr, ich) * pmom(1:numberStreams, ilyr) ) / &
													  ( use_tau(ilyr)*use_omega(ilyr) + aerotau(ilyr, ich) * omegaaero(ilyr, ich) )
						endif
				
	
						use_omega(ilyr) = ( use_tau(ilyr)*use_omega(ilyr) + aerotau(ilyr, ich)*omegaaero(ilyr,ich) ) / &
												( use_tau(ilyr) + aerotau(ilyr, ich) )
						use_tau(ilyr) = use_tau(ilyr) + aerotau(ilyr, ich)

					endif

				!!!!!!!!! CLOUD-FREE RUN, DO NOT ADD THEM  
			!		   cloud_layers(ilyr) = 0
				

					if (cloud_layers(ilyr) == 1) then 
						
								
					! scaling optical thickness to reference wavelength of 0.65um
						tauin  = cloud_tau_layers(ilyr)
						otau = tauin * extinction_efficiency(icld, ich) / extinction_efficiency(icld, 1)
						omega = single_scatter_albedo(icld, ich)
						
					! scaling optical thickness for phase function truncation
						otau = otau * ( 1. - truncation_factor(icld, ich) * single_scatter_albedo(icld, ich))
						omega = omega * (1. - truncation_factor(icld, ich) ) / ( 1. - truncation_factor(icld, ich)  * omega)

						
					! combine cloud and atmospheric tau, original omega may or may not be 0, so can't just skip
					! uniformly like we used to do. 	

						pmom(1:numberStreams, ilyr) = ( use_tau(ilyr) * use_omega(ilyr) * pmom(1:numberStreams, ilyr) + &
														otau * omega * legendre_coefficients(1:numberStreams, icld, ich) ) / &
														( use_tau(ilyr)*use_omega(ilyr) + otau*omega )


						use_omega(ilyr) = ( use_tau(ilyr)*use_omega(ilyr) + otau*omega) / &
												( use_tau(ilyr) + otau)
!						use_omega(ilyr) = omega / ( 1. + use_tau(ilyr) / otau)
						use_tau(ilyr) = use_tau(ilyr) + otau 
!						pmom(0:numberStreams, ilyr) = legendre_coefficients(0:numberStreams, icld, ich)

						icld = icld + 1

					endif	
	
				end do ! end layer loop
					
				pmom(0, :) = 1.
					
					

			! now it's time to actually call disort
				uu = 0.
				uus = 0.


					umu(0) = sensor_zenith(1)
					phi(1) = relative_azimuth(1)			
																
					NLYR = use_layers ! numberLayers		
					NLEV = use_levels !numberLevels	
					NSTR = numberStreams
					WNL = wnum_low(ich)*1.
					WNH = wnum_hi(ich)*1.
					SFX = solarFlux(ich)
					SFA = surfaceAlbedo(ich)
					NTAU = 1
					MUMU = maxumu
					MPHI = maxphi
					MCMU = maxcmu
								
				call DISORT_DP_GC (NLYR, use_tau, use_omega, NSTR, &
								pmom, Tatm, WNL, &
								WNH, usrtau, NTAU, utau, NSTR, usrang, &
								numu, sensor_zenith, nphi, relative_azimuth,  ibcnd, SFX, &
								solar_zenith, phi0, fisot, lambertian, SFA, Tsfc, Ttoa, &
								temis, local_planck, onlyfl, accur, prnt, header, NLYR, &
								NLEV, MUMU, MPHI, MCMU, uu, uus)

	!			if (ich == 19) print*, ich, NLYR, use_tau, use_omega
								
				! the middle index of uu is tau, so if more cloud taus, we need to put in an extra loop here				
				radiance = radiance + uu(1,1,1)*correlated_Ks_weights(ich, ik)
					
				end do ! end k-loop


				if (refl_or_rad(ich) == 1) reflectance = pi * radiance / (solarFlux(ich)*solar_zenith)
				radiance = radiance / ( (1./wnum_low(ich)) - (1./wnum_hi(ich))) / 1.0e4

				if (refl_or_rad(ich) == 1) then 
					output_data(ich) = reflectance
				else
					output_data(ich) = radiance
				endif

                     endif ! end the skip channels for 64-streams loop. 


			
		end do ! end ch-loop


!		print*, "channels:", waveLengths
!		print*, "output:", output_data

		take_out = output_data
	
	
	end subroutine run_disort
	
	
	
	integer function locate_re(radius, radii_array)
	
		real, intent(in) :: radius
		integer, dimension(:), intent(in) :: radii_array
		
		integer :: i, nrad
		
		nrad = size(radii_array)
		
		do i=1, nrad
			if ( radius <= radii_array(i)) exit
		end do
	
		locate_re = i
	
	end function locate_re
	
	
	
	
	
	subroutine cleanup_disort
	
		deallocate(single_scatter_albedo)
		deallocate(truncation_factor)
		deallocate(extinction_efficiency)
		deallocate(legendre_coefficients)
	
	
	
	
	end subroutine cleanup_disort
	
	
	
	
	
	


end module modis_simulator
