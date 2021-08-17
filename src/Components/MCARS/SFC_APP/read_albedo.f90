! This program is an example of use of the albedo module. 


! To-do: 

! Separate executable that will now use this file
! 1. run the correlated-K module, no output to disk
! 2. take MOD06 cloud information to fill in the cloud positions
! 3. run DISORT (new one that supposedly does GFORTRAN) to simulate the radiances
!    over an entire granule's worth. 
! 4. Output that data to a file and visualize the resulting reflectance/radiance
! 5. send it to Arlindo and grin evilly




program read_albedo
	
	use general_array_io, only:read_float_array, write_float_array, write_int_array, read_int_array, GetWidHt, write_byte_array
	use global_model_grids
	use modis_albedo
	use emissivity
	use modis_numerical_module, only: linearinterpolation
	use geos5_io_module
	
	use hdf5
	
	implicit none
	
	
	
	character(len=200) :: MOD03_name, output_name, &
		emissivity_name, ocean_name, IGBP_name, &
		geos5_file
		
! there are 7 albedo files
	integer, parameter :: nbands = 7
	integer, parameter :: num_albs = 46
	integer, parameter :: albedo_tags(num_albs) = (/1, 9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89, &
		97, 105, 113, 121, 129, 137, 145, 153, 161, 169, 177, 185, 193, 201, 209, 217, 225, 233, &
		241, 249, 257, 265, 273, 281, 289, 297, 305, 313, 321, 329, 337, 345, 353, 361 /)
	real, parameter :: albedo_factor = 0.001

	integer*2, parameter :: season_break(4) = (/ 60, 152, 244, 335 /)
	
	integer, parameter :: emis_channels(11) = (/ 20, 22, 27, 28, 29, 31, 32, 33, 34, 35, 36 /)

	
	integer :: modis_wid, modis_ht, start(2), stride(2), edge(2)
	integer :: year, day_of_year, file_id(1), err_code, i, fid(1), edg(2), albedo_year

	integer*1, dimension(:,:), allocatable :: ecosystem_flag
	real, dimension(:,:), allocatable :: latitude, longitude, snow_ice_fraction, solar_zenith, sensor_zenith
	real, dimension(:,:), allocatable :: solar_azimuth, sensor_azimuth, relative_azimuth_angle
	integer*2, dimension(:,:,:), allocatable :: surface_albedo, emissive_albedo
	real, dimension(:,:,:), allocatable :: ocean_albedo
	real, dimension(:,:,:), allocatable :: temp_albedo
	
	character(len=200) :: albedo_files(nbands)
	character(len=200) :: albedo_path, albedo_dir, albedo_single
	character(len=3) :: doy_tag
	character(len=400) :: snow_albedo_filename, db_albedo_name
	
	integer :: num_Iter, lat_start, lat_end, j
	integer :: start_time, end_time, crate, cmax
	character(len=10) :: platform_name
	logical :: useoffset
	real*8, parameter :: d2r = 0.017453292519943295d0
	integer*1 :: leap_offset
	integer*1, dimension(:,:), allocatable :: land_mask
	

	integer, parameter :: nszen = 33, nvzen = 28, nraz = 37, nws = 3
	real :: lib_szen(nszen), lib_vzen(nvzen), lib_phi(nraz)
	real :: library_ocean(nszen, nvzen, nraz, 15, nws)
	real :: albedo_vector(15)
	
	character(len=300) :: parfile
	
	character(len=10) :: band_dir(nbands)

! GEOS-5 QUIERY

	integer :: grid_xsize, grid_xstart, iend
	integer :: grid_ysize, grid_ystart, jend
	
    real, allocatable, dimension(:,:)  :: wind_speed
	
	real, dimension(:,:), allocatable :: met_temp1, met_temp2, met_temp3, met_temp4
	real :: xx(2), yy1(2), yy2(2)
	real :: met_date(2), frac_time
   
   real, dimension(:,:), allocatable :: Tsfc, colO3
	integer :: MYTIME, MYMONTH
	integer, parameter :: model_layers = 101
	integer :: model_i, model_j
	real :: snow, ice 
	
	integer :: gfile_id(2)

! END GEOS-5 QUIERY	

	band_dir = (/ "MCD43GS066", "MCD43GS086", "MCD43GS047", "MCD43GS056",   &
		      "MCD43GS124", "MCD43GS164", "MCD43GS213"  /)

	
	call system_clock(start_time, crate, cmax)
	
	
	
	
	lib_szen(1:13)   = (/ (0.05 * i+0.15, i=0,12) /) !numOfSolarAngs = 33, numOfViewAngs = 28
	lib_vzen(1:8)     = (/(0.05 * i + 0.40,i=0,7) /)
	lib_szen(14:nszen)   = (/ (0.0125 * i+0.7625, i=0,19) /) !numOfSolarAngs = 33, numOfViewAngs = 28
	lib_vzen(9:nvzen)     = (/(0.0125 * i + 0.7625,i=0,19) /)
    lib_phi(1:nraz)= (/ (5.0 * i,i=0,nraz-1)/)  

	call getarg(1, parfile)
	call read_par_file(parfile)
	
	print*, year, day_of_year, MYTIME
	print*, "mod03: ", MOD03_name
	print*, "output: ", output_name
	print*, "GEOS: ", geos5_file
	print*, "alb_path: ", albedo_path
	
! select the albedo directory to read from based on year/date information
! this code needs to go wherever the initial decision on input files is being made
	albedo_year = year
	if (albedo_year > 2016) albedo_year = 2016

	do i=1, num_albs
		if (day_of_year < albedo_tags(i)) exit
	end do
	i=i-1

	if (day_of_year < 17) then 
		write(doy_tag, '(a2,i1)') "00", albedo_tags(i)
	else if (day_of_year >= 17 .and. day_of_year < 105) then 
		write(doy_tag, '(a1,i2)') "0", albedo_tags(i)
	else
		write(doy_tag, '(i3)') albedo_tags(i)
	endif

	write(albedo_dir, '(i4,a1,a3,a1)') albedo_year, "/", doy_tag , "/"
	print*, trim(albedo_dir)

! generate the filenames as they contain year/date information
	do i=1, nbands
		write(albedo_single, '(a16,i1,a1,a3,a1,i4,a9)') "MCD43GF_wsa_Band", i, "_", doy_tag, "_", albedo_year, "_V006.hdf"
		albedo_files(i) = trim(albedo_path) // trim(band_dir(i)) // "/" // trim(albedo_dir) // trim(albedo_single)
		print*, trim(albedo_files(i))
	end do

	

! now we know what to read, grab a granule and get its information
! in reality this information will come from somewhere else within the model
! so this is technically junk code.

	file_id(1) = sfstart(MOD03_name, DFACC_READ)
	call GetWidHt(file_id, "Latitude", modis_wid, modis_ht)

! we will go to 10 km
	modis_wid = modis_wid 
	modis_ht = modis_ht 

	allocate(latitude(modis_wid, modis_ht), &
			 longitude(modis_wid, modis_ht), &
			 ecosystem_flag(modis_wid, modis_ht), &
			 surface_albedo(modis_wid, modis_ht, 13), &
			 emissive_albedo(modis_wid, modis_ht, 2), &
			 ocean_albedo(modis_wid, modis_ht, 7), &
			 snow_ice_fraction(modis_wid, modis_ht))
			 
	allocate(solar_zenith(modis_wid, modis_ht), &
			 sensor_zenith(modis_wid, modis_ht), &
			 solar_azimuth(modis_wid, modis_ht), &
			 sensor_azimuth(modis_wid, modis_ht), &
			 relative_azimuth_angle(modis_wid, modis_ht))
			  
			 
	allocate(land_mask(modis_wid, modis_ht))		 
			 
	start = 0
	stride = 1
	edge = (/ modis_wid, modis_ht /)

	call read_float_array(file_id, "Latitude", start, stride, edge, latitude, err_code)
	call read_float_array(file_id, "Longitude", start, stride, edge, longitude, err_code)
	

	useoffset = .false.
	call read_int_array(file_id, "SensorZenith", start, stride, edge, useoffset, sensor_zenith, err_code)
	call read_int_array(file_id, "SensorAzimuth", start, stride, edge, useoffset, sensor_azimuth, err_code)
	call read_int_array(file_id, "SolarZenith", start, stride, edge, useoffset, solar_zenith, err_code)
	call read_int_array(file_id, "SolarAzimuth", start, stride, edge, useoffset, solar_azimuth, err_code) 

	err_code = sfend(file_id(1))

	solar_zenith = cos(solar_zenith*d2r)
	sensor_zenith = cos(sensor_zenith*d2r)


!  calculate the relative azimuth (From MOD06)

   relative_azimuth_angle =  solar_azimuth + 180. - sensor_azimuth

	deallocate(solar_azimuth)
	deallocate(sensor_azimuth)


	do j = 1,  modis_ht
		do i = 1, modis_wid
	 
       if (relative_azimuth_angle(i,j) <= 0.) relative_azimuth_angle(i,j) = -relative_azimuth_angle(i,j) 
       if (relative_azimuth_angle(i,j) > 180.) relative_azimuth_angle(i,j) = 360. - relative_azimuth_angle(i,j)  
     enddo
   enddo

   relative_azimuth_angle = abs(relative_azimuth_angle)
   where(relative_azimuth_angle > 180.) relative_azimuth_angle =360. - relative_azimuth_angle
	
! now we are read in the ocean surface albedo table from our snazzy lookup table that 
! we've calculated from a DISORT run. 

	file_id(1) = sfstart(ocean_name, DFACC_READ)
	call read_float_array(file_id, "Ocean_Reflectance", (/ 0,0,0,0,0 /), (/ 1,1,1,1,1 /), &
								(/ nszen, nvzen, nraz, 13, nws /), library_ocean, err_code)
	err_code = sfend(file_id(1))

	print*, library_ocean(1,1,1,1,:)

! GEOS-5 QUIERY
	file_id(1) = sfstart(output_name, DFACC_CREATE)


    allocate(wind_speed(modis_wid, modis_ht))
    allocate(met_temp1(modis_wid, modis_ht))
    allocate(met_temp2(modis_wid, modis_ht))
    
   
	call h5open_f(err_code)

	call h5fopen_f(geos5_file, H5F_ACC_RDONLY_F, gfile_id(1), err_code)	
   
	call read_float_array_g5(gfile_id, 1, "U10M", (/0,0,0/), (/1,1,1/), (/ modis_wid, modis_ht, 1/), 2, met_temp1, err_code)
	call read_float_array_g5(gfile_id, 1, "V10M", (/0,0,0/), (/1,1,1/), (/ modis_wid, modis_ht, 1/), 2, met_temp2, err_code)
   
	do j=1, modis_ht
		do i=1, modis_wid
		
			wind_speed(i,j) = sqrt( met_temp1(i,j) ** 2 + met_temp2(i,j) ** 2 )
		
		end do
	end do	

	call read_float_array_g5(gfile_id, 1, "FRSEAICE", (/ grid_xstart, grid_ystart, 0/), (/1,1,1/), (/  modis_wid, modis_ht, 1/), 2, met_temp1, err_code)
	call read_float_array_g5(gfile_id, 1, "FRSNO", (/ grid_xstart, grid_ystart, 0/), (/1,1,1/), (/  modis_wid, modis_ht, 1/), 2, met_temp2, err_code)

	snow_ice_fraction = 0.
	
	do j=1, modis_ht
		do i=1, modis_wid
				
			ice = met_temp1(i,j)
			snow = met_temp2(i,j)
			
			if (ice > 0. .and. ice <= 1.) snow_ice_fraction(i,j) = ice
			if (snow > 0.5 .and. snow <= 1.) snow_ice_fraction(i,j) = snow
			
		end do
	end do
		   
	call h5fclose_f(gfile_id(1), err_code)
	
	call write_float_array(file_id, "snow_ice_fraction", (/0,0/), (/1,1/), (/modis_wid, modis_ht/), snow_ice_fraction, err_code)

	deallocate(met_temp1, met_temp2)
		
	call h5close_f(err_code)

! END GEOS-5 QUIERY

! at this point, we are ready to read the albedo. The albedo has to be read in pieces in order to not overwhelm
! the system memory. The albedo files are very large. 

! we will take out the NISE stuff later once we know what exactly GEOS-5 can give us. 

	num_iter = modis_ht / 100 + 1

	do i=1, num_iter
  
		lat_start = (i-1)*100 + 1
		if (i==num_iter) then 
			lat_end = modis_ht
		else
			lat_end = lat_start + 99
		endif

		print*, lat_start, lat_end


		call getAlbedoEco(latitude(:, lat_start:lat_end), &
						longitude (:, lat_start:lat_end), &
						day_of_year, albedo_files, IGBP_name, &
						snow_albedo_filename, &
						snow_ice_fraction(:, lat_start:lat_end), &
						surface_albedo(:, lat_start:lat_end, 1:7), &
						land_mask(:, lat_start:lat_end) )

	end do

! write the output to disk. 
	stride = 1
	
	call write_float_array(file_id, "Phi", start, stride, edge, relative_azimuth_angle, err_code)
	

! now calculate surface albedo for channels 17, 18, 19 and 26 based on information we already have
! 18 and 19 share the albedo because they are basically the same channel. 

	do j=1, modis_ht
		do i=1, modis_wid
		
			surface_albedo(i,j,10) = int( linearinterpolation( (/ 0.858, 1.24 /), (/ surface_albedo(i,j,2)*1., surface_albedo(i,j,5)*1. /), &
										0.905 ) )										
			surface_albedo(i,j,11) = int( linearinterpolation( (/ 0.858, 1.24 /), (/ surface_albedo(i,j,2)*1., surface_albedo(i,j,5)*1. /), &
										0.94 ) )
			surface_albedo(i,j,12) = int( linearinterpolation( (/ 0.858, 1.24 /), (/ surface_albedo(i,j,2)*1., surface_albedo(i,j,5)*1. /), &
										0.94 ) )
			surface_albedo(i,j,13) = linearinterpolation( (/ 1.24, 1.64 /), (/ surface_albedo(i,j,5)*1., surface_albedo(i,j,6)*1. /), &
										1.375 )
		
		end do
	end do


! now we will get surface albedo for 412 and 430nm channels 8 and 9
	if ( mod(year,4) == 0 .and. mod(year,400) /= 0) then 
		leap_offset = 1
	else
		leap_offset = 0
	endif

	if ( day_of_year+leap_offset < season_break(1)+leap_offset .or. day_of_year+leap_offset >= season_break(4)+leap_offset) &
		db_albedo_name = "LIB/aqua_modis_surfdb_winter_20110913.hdf"
	if ( day_of_year+leap_offset >= season_break(1)+leap_offset .and. day_of_year+leap_offset < season_break(2)+leap_offset) &
		db_albedo_name = "LIB/aqua_modis_surfdb_spring_20110913.hdf"
	if ( day_of_year+leap_offset >= season_break(2)+leap_offset .and. day_of_year+leap_offset < season_break(3)+leap_offset) &
		db_albedo_name = "LIB/aqua_modis_surfdb_summer_20110913.hdf"
	if ( day_of_year+leap_offset >= season_break(3)+leap_offset .and. day_of_year+leap_offset < season_break(4)+leap_offset) &
		db_albedo_name = "LIB/aqua_modis_surfdb_fall_20110913.hdf"
		
		
		
	allocate(temp_albedo(modis_wid, modis_ht, 1))	
	
	call get_db_albedo(db_albedo_name, latitude, longitude, temp_albedo(:,:,1))

	surface_albedo(:,:,8) = int(temp_albedo(:,:,1)*10.)	
! we need to match resolutions on these suckers as much as possible
	where(land_mask(:,:) == 0) surface_albedo(:,:,8) = 32767 
		
	do j=1, modis_ht
		do i=1, modis_wid
			surface_albedo(i,j,9) = int( linearinterpolation( (/ 0.421, 0.469 /), (/ surface_albedo(i,j,1)*1., surface_albedo(i,j,3)*1. /), &
										0.443 ) )
		end do
	end do
		
	
! this is where we fill in the surface albedo. Wind speed has been read from GEOS-5
	print*, lib_szen
	print*, lib_vzen
	print*, lib_phi
		
	do j=1, modis_ht
		do i=1, modis_wid


			if (land_mask(i,j) == 0 .and. snow_ice_fraction(i,j) == 0. ) then ! this means ocean and no snow/ice
		


				call interpolate_albedo(library_ocean, &
							lib_szen, lib_vzen, lib_phi, &
							wind_speed(i,j), &
							solar_zenith(i,j), sensor_zenith(i,j), relative_azimuth_angle(i,j), &
							 albedo_vector)
		
				surface_albedo(i,j,:) = int( albedo_vector(1:13)*1000. )
				emissive_albedo(i,j,1:2) = int( albedo_vector(14:15)*1000. )
		
			endif
		
		end do
		
	end do
		
	call write_float_array(file_id, "wind_speed", start, stride, edge, wind_speed, err_code)	
		
	deallocate(wind_speed)
	
	where (snow_ice_fraction /= 0.) land_mask = 1
	call write_byte_array(file_id, "Ocean_Mask", start, stride, edge, land_mask, err_code)

	
	call write_int_array(file_id, "Surface_Albedo_1", start, stride, edge, surface_albedo(:,:,1), err_code)	
	call write_int_array(file_id, "Surface_Albedo_2", start, stride, edge, surface_albedo(:,:,2), err_code)	
	call write_int_array(file_id, "Surface_Albedo_3", start, stride, edge, surface_albedo(:,:,3), err_code)	
	call write_int_array(file_id, "Surface_Albedo_4", start, stride, edge, surface_albedo(:,:,4), err_code)	
	call write_int_array(file_id, "Surface_Albedo_5", start, stride, edge, surface_albedo(:,:,5), err_code)	
	call write_int_array(file_id, "Surface_Albedo_6", start, stride, edge, surface_albedo(:,:,6), err_code)	
	call write_int_array(file_id, "Surface_Albedo_7", start, stride, edge, surface_albedo(:,:,7), err_code)	
	call write_int_array(file_id, "Surface_Albedo_8", start, stride, edge, surface_albedo(:,:,8), err_code)	
	call write_int_array(file_id, "Surface_Albedo_9", start, stride, edge, surface_albedo(:,:,9), err_code)	
	call write_int_array(file_id, "Surface_Albedo_17", start, stride, edge, surface_albedo(:,:,10), err_code)	
	call write_int_array(file_id, "Surface_Albedo_18", start, stride, edge, surface_albedo(:,:,11), err_code)	
	call write_int_array(file_id, "Surface_Albedo_19", start, stride, edge, surface_albedo(:,:,12), err_code)	
	call write_int_array(file_id, "Surface_Albedo_26", start, stride, edge, surface_albedo(:,:,13), err_code)	

	print*, "1,2,3,4,5,6,7,8,9,17,18,19,26"
	print*, surface_albedo(630, 305, 1:13)



	err_code = sfend(file_id(1))


	deallocate(snow_ice_fraction)
	deallocate(surface_albedo)
	deallocate(temp_albedo)




! now that we've extracted the surface albedo from MOD43 files, we need to add the surface albedo 
! for all the emissive channels, of which there are 11: 20, 22, 27, 28, 29, 31, 32, 33, 34, 35 and 36

	allocate(surface_albedo(modis_wid, modis_ht, 11), temp_albedo(modis_wid, modis_ht, 11))

	platform_name = "Aqua"
	call get_emissivity(emissivity_name, emis_channels, platform_name, latitude, longitude, temp_albedo)
	
	surface_albedo = int( (1. - temp_albedo) * 1000. )
! match the resolutions with original surface albedo
	do j=1, modis_ht
		do i=1, modis_wid
		
			if (land_mask(i,j) == 0) then
				surface_albedo(i,j,1:2) = emissive_albedo(i,j,1:2)
				surface_albedo(i,j,3:11) = int( (1-0.985)*1000.)
			endif
		
		end do
	end do

	deallocate(temp_albedo)
	
	file_id(1) = sfstart(output_name, DFACC_WRITE)
	call write_int_array(file_id, "Surface_Albedo_20", start, stride, edge, surface_albedo(:,:,1), err_code)	
	call write_int_array(file_id, "Surface_Albedo_22", start, stride, edge, surface_albedo(:,:,2), err_code)	
	call write_int_array(file_id, "Surface_Albedo_27", start, stride, edge, surface_albedo(:,:,3), err_code)	
	call write_int_array(file_id, "Surface_Albedo_28", start, stride, edge, surface_albedo(:,:,4), err_code)	
	call write_int_array(file_id, "Surface_Albedo_29", start, stride, edge, surface_albedo(:,:,5), err_code)	
	call write_int_array(file_id, "Surface_Albedo_31", start, stride, edge, surface_albedo(:,:,6), err_code)	
	call write_int_array(file_id, "Surface_Albedo_32", start, stride, edge, surface_albedo(:,:,7), err_code)	
	call write_int_array(file_id, "Surface_Albedo_33", start, stride, edge, surface_albedo(:,:,8), err_code)	
	call write_int_array(file_id, "Surface_Albedo_34", start, stride, edge, surface_albedo(:,:,9), err_code)	
	call write_int_array(file_id, "Surface_Albedo_35", start, stride, edge, surface_albedo(:,:,10), err_code)	
	call write_int_array(file_id, "Surface_Albedo_36", start, stride, edge, surface_albedo(:,:,11), err_code)	
	
	print*, "20,22,27,28, 29, 31, 32, 33, 34, 35, 36"
	print*, surface_albedo(630, 305, 1:11)
	
	err_code = sfend(file_id(1))

	deallocate(latitude, longitude, ecosystem_flag, surface_albedo, ocean_albedo)
	deallocate(land_mask)
			
	call system_clock(end_time, crate, cmax)
	
	print*, "Time: ", (end_time - start_time) / (crate*1.0) / 60., "min"


contains

	subroutine read_par_file(parfile_name)
	
		character*(*), intent(in) :: parfile_name
		character(len=200) :: dummy, input_dir, temp
		character(len=4) :: ytemp, doytemp, ttemp
		
		open (unit=20, file=parfile_name)
		
		read(20,*) dummy
		read(20, *) ytemp, doytemp, ttemp

		read(ytemp, *) year
		read(doytemp, *) day_of_year
		read(ttemp, *) MYTIME

		input_dir = "inputs-" // trim(ytemp) // trim(doytemp) // "." // trim(ttemp) // "/"
			
		read(20,*) dummy
		read(20,'(a)') temp
		MOD03_name = trim(input_dir) // trim(temp)
		
		read(20,*) dummy
		read(20,'(a)') temp
		output_name = trim(input_dir) // trim(temp)

		read(20,*) dummy
		read(20,*) dummy

		read(20,*) dummy
		read(20,'(a)') temp
		geos5_file = trim(input_dir) // trim(temp)

		read(20,*) dummy
		read(20,'(a)') temp
		albedo_path = trim(temp)
		
! The ecosystem flag should come from IGBP

		snow_albedo_filename = "LIB/AlbSnwSts.ByNISE.W90.D90.WS.Hemi.2000-2004.YrAvg.hdf"
		emissivity_name = "LIB/global_emiss_intABI_2005001.hdf"
		ocean_name = "LIB/Ocean_Reflectance.hdf"
		IGBP_name = "LIB/IGBP.EcoMap.NtoS.2004.149.v004.hdf"

		close(20)
	
	end subroutine read_par_file	
	
	

end program read_albedo
