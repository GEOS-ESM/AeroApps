module emissivity

	implicit none
	
	
	
! bands list is a list of channels for each platform that has emissivity map associated with it 

contains

	subroutine get_emissivity (emissivity_name, channels, platform_name, lat, lon, sfc_emis_array)


	include "hdf.f90"
	include "dffunc.f90"

	
		character(*), intent(in) :: emissivity_name
		integer, dimension(:), intent(in) :: channels
		character(*), intent(in) :: platform_name
		real, dimension(:,:), intent(in) :: lat, lon
		real, dimension(:,:,:), intent(inout) :: sfc_emis_array
			
			
		integer :: file_id, var_id, err_code, start(2), stride(2), edge(2)
		integer :: channel_number, num_channels, i , j

		character(len=20) :: sds_name


	! emissivity map variables
		real, parameter :: emissivity_res = 0.05
		integer, parameter :: num_emiss_cols = 7200
		integer, parameter :: num_emiss_rows = 3600
		real, dimension(:,:,:), allocatable :: emissivity_map
		integer*2, dimension(:,:), allocatable :: emissivity_map_int
		integer :: emiss_x, emiss_y, NumberOfXPoints, NumberOfYPoints, NumMapCols, NumMapRows
		real :: emiss_west, emiss_north
		real, parameter :: emissivity_fill = 0.985 ! according to Wisconsin
		integer :: startx, starty, stopy, stopx
		real, parameter :: ScaleFactor = 0.001
	
		real :: minLat, minLong, maxLat, maxLong
		integer :: wid, ht
	
		num_channels = size(channels)

		wid = size(lat, 1)
		ht = size(lat, 2)

    !Determine min and max:
		if (maxval(Lat) == -999. .or. maxval(Lon) == -999.) then
		!there is no gelocation data in the scan.
			sfc_emis_array(:,:,:) = -999.
			return
		else
			minLat  = minval( Lat,  mask = ( (Lat  >= -90.0000 .and. Lat  <= 90.00000) ) )
			maxLat  = maxval( Lat,  mask = ( (Lat  >= -90.0000 .and. Lat  <= 90.00000) ) )
			minLong = minval( Lon, mask = ( (Lon >= -180.000 .and. Lon <= 180.0000) ) )
			maxLong = maxval( Lon, mask = ( (Lon >= -180.000 .and. Lon <= 180.0000) ) )
		endif


    !Compute the WesternMostLongitude and NothernMostLatitude of the maps:
		emiss_west = emissivity_res / 2.0 - 180.0
		emiss_north = 90.0 - emissivity_res / 2.0
    
    !Compute the start/stop x/y indices based on min/max lat/lon values:
	!use anint so that it rounds to nearest index.
		startx = anint( (minLong-emiss_west) / (emissivity_res) ) + 1
		stopx  = anint( (maxLong-emiss_west) / (emissivity_res) ) + 1
		starty = anint( (emiss_north-maxLat ) / (emissivity_res) ) + 1
		stopy  = anint( (emiss_north-minLat ) / (emissivity_res) ) + 1
    
    !Due to roundoff, geocoordinates of -90, 90, -180, 180 can result in 
    ! indices of one beyond the map bounds, so correct that here:
		if (startx == 0) startx = 1
		if (stopx  == num_emiss_cols+1 ) stopx = num_emiss_cols
		if (starty == 0) starty = 1
		if (stopy  == num_emiss_rows+1 ) stopy = num_emiss_rows
    
    !Compute the number of points:
		NumberOfXPoints = stopx - startx + 1
		NumberOfYPoints = stopy - starty + 1

    !Error Checking:
		if ( (startx < 1) .or. (stopx > num_emiss_cols) .or. &
			 (starty < 1) .or. (stopy > num_emiss_rows) .or. &
			 (startx > stopx) .or. (starty > stopy)      ) then

			if (stopx > NumMapCols) stopx = NumMapCols
			if (stopy > NumMapRows) stopy = NumMapRows
			if (startx > stopx) startx = stopx-1
			if (starty > stopy) starty = stopy-1
			if (startx < 1) startx = 1
			if (starty < 1) starty = 1

		end if
	
    !Define the starting point
		start (1) = startX - 1
		start (2) = startY - 1
    
    !Define the stride = 1 so that it doesn't skip any data.
		stride(:) = 1
    
    !Define the Number of Points to read:
		edge  (1) = NumberOfXPoints
		edge  (2) = NumberOfYPoints	
	
		allocate( emissivity_map    (startx:stopx, starty:stopy, num_channels) )
		allocate( emissivity_map_int (startx:stopx, starty:stopy      ) )

		file_id = sfstart(trim(emissivity_name), DFACC_READ)

		do i=1, num_channels
		
			channel_number = channels(i)
		
			if (platform_name(1:5) .eq. 'Terra' .or. &
				platform_name(1:5) .eq. 'terra' .or. &
				platform_name(1:5) .eq. 'TERRA' .or. &
				platform_name(1:4) .eq. 'Aqua' .or. & 
				platform_name(1:4) .eq. 'aqua' .or. &
				platform_name(1:4) .eq. 'AQUA'	 ) then 
	 
				if (channel_number >= 20 .and. channel_number <= 25 ) sds_name = "emiss7" ! 3.7 
				if (channel_number == 27 ) sds_name = "emiss9" ! 6.2
				if (channel_number == 28 ) sds_name = "emiss10" ! 7.3
				if (channel_number == 29 ) sds_name = "emiss11" ! 8.5
				if (channel_number == 31 ) sds_name = "emiss14" ! 11.0
				if (channel_number == 32 ) sds_name = "emiss15" ! 12.0
				if (channel_number > 32) sds_name = "emiss16"   ! any CO2 bands
			
			endif
	
			if (platform_name(1:3) .eq. 'Mas' .or. &
				platform_name(1:3) .eq. 'mas' .or. &
				platform_name(1:3) .eq. 'MAS') then
		
				if (channel_number == 30) sds_name = "emiss7" ! 3.7
				if (channel_number == 42) sds_name = "emiss11" ! 8.5
				if (channel_number == 45) sds_name = "emiss14" ! 11
				if (channel_number == 47) sds_name = "emiss15" ! 12
				if (channel_number > 47) sds_name = "emiss16"   ! any CO2 bands
	
			endif

			if (platform_name(1:6) .eq. 'Master' .or. &
				platform_name(1:6) .eq. 'master' .or. &
				platform_name(1:6) .eq. 'MASTER') then
		
				if (channel_number == 30) sds_name = "emiss7" ! 3.7
				if (channel_number == 43) sds_name = "emiss11" ! 8.5
				if (channel_number == 47) sds_name = "emiss14" ! 11
				if (channel_number == 49) sds_name = "emiss15" ! 12
				if (channel_number > 49) sds_name = "emiss16" ! SHIS bands may be present
	
			endif

			if (platform_name(1:6) .eq. 'Seviri' .or. &
				platform_name(1:6) .eq. 'seviri' .or. &
				platform_name(1:6) .eq. 'SEVIRI') then
			
				if (channel_number == 4 ) sds_name = "emiss7" ! 3.7 
				if (channel_number == 6 ) sds_name = "emiss10" ! 7.3
				if (channel_number == 7 ) sds_name = "emiss11" ! 8.5
				if (channel_number == 9 ) sds_name = "emiss14" ! 11.0
				if (channel_number == 10 ) sds_name = "emiss15" ! 12.0						
	
			endif


			var_id = sfselect(file_id, sfn2index(file_id, trim(sds_name)))
			err_code = sfrdata(var_id, start, stride, edge, emissivity_map_int)
			err_code = sfendacc(var_id)
			
			emissivity_map(:,:,i) = emissivity_map_int * ScaleFactor
			where (emissivity_map(:,:,i) < 0 ) emissivity_map(:,:,i) = emissivity_fill
			
		end do
		
		err_code = sfend(file_id)
		deallocate(emissivity_map_int)


		do j=1, ht
			do i=1, wid
		
				if ( (Lat (i,j) >= -90.0000) .and. &
					(Lat (i,j) <=  90.0000) .and. &
					(Lon (i,j) >= -180.000) .and. &
					(Lon (i,j) <=  180.000) ) then 
		
! put the emissivity values where they belong 
					emiss_x = anint( (Lon(i,j)-emiss_west) / (emissivity_res) ) + 1
					emiss_y = anint( (emiss_north-Lat (i,j)) / (emissivity_res) ) + 1

					if (emiss_x == 0) emiss_x = 1
					if (emiss_x  == num_emiss_cols+1 ) emiss_x = num_emiss_cols
					if (emiss_y == 0) emiss_y = 1
					if (emiss_y  == num_emiss_rows+1 ) emiss_y = num_emiss_rows


					sfc_emis_array(i,j, :) = emissivity_map(emiss_x, emiss_y, :)
		
				endif
		
			end do
		end do
		
	
	end subroutine get_emissivity


	subroutine get_seviri_emissivity(emissivity_name, channels, start, stride, edge, sfc_emis_array)


	include "hdf.f90"
	include "dffunc.f90"

	
		character(*), intent(in) :: emissivity_name
		integer, dimension(:), intent(in) :: channels
		real, dimension(:,:,:), intent(inout) :: sfc_emis_array
		integer, intent(in) :: start(:), stride(:), edge(:)
			
		integer :: file_id, var_id, err_code
		integer :: channel_number, num_channels, i , j

		character(len=20) :: sds_name


	! emissivity map variables
		
		real, parameter :: emissivity_fill = 0.985 ! according to Wisconsin
		real, parameter :: ScaleFactor = 0.001
		integer*2, dimension(:,:), allocatable :: emissivity_map_int
		
		num_channels = size(channels)

		allocate(emissivity_map_int(edge(1), edge(2)))

		file_id = sfstart(trim(emissivity_name), DFACC_READ)

		
		do i=1, num_channels
		
			channel_number = channels(i)
		
			
			if (channel_number == 4 ) sds_name = "emiss7" ! 3.7 
			if (channel_number == 6 ) sds_name = "emiss10" ! 7.3
			if (channel_number == 7 ) sds_name = "emiss11" ! 8.5
			if (channel_number == 9 ) sds_name = "emiss14" ! 11.0
			if (channel_number == 10 ) sds_name = "emiss15" ! 12.0						


			var_id = sfselect(file_id, sfn2index(file_id, trim(sds_name)))
			err_code = sfrdata(var_id, start, stride, edge, emissivity_map_int)
			err_code = sfendacc(var_id)
			
			sfc_emis_array(:,:,i) = emissivity_map_int * ScaleFactor
			where (sfc_emis_array(:,:,i) < 0. ) sfc_emis_array(:,:,i) = emissivity_fill
			
		end do
		
		err_code = sfend(file_id)
		deallocate(emissivity_map_int)

	
	
	
	
	end subroutine get_seviri_emissivity

	

end module emissivity