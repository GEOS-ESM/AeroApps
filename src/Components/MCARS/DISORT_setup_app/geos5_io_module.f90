module geos5_io_module

	use hdf5
	
	implicit none
	
contains


	subroutine find_geos5_bounds(istart, iend, jstart, jend, lat, lon)

		integer, intent(inout) :: istart, iend, jstart, jend
		real, dimension(:,:), intent(in) :: lat, lon

		real :: minLat, maxLat, minLon, maxLon
		real, parameter :: dlon = 5./16.
		real, parameter :: dlat = 0.25


		minLat  = minval( Lat,  mask = ( (Lat  >= -90.0000 .and. Lat  <= 90.00000) ) )
		maxLat  = maxval( Lat,  mask = ( (Lat  >= -90.0000 .and. Lat  <= 90.00000) ) )
		minLon = minval( Lon, mask = ( (Lon >= -180.000 .and. Lon <= 180.0000) ) )
		maxLon = maxval( Lon, mask = ( (Lon >= -180.000 .and. Lon <= 180.0000) ) )
		
		
		istart = int ((minlon + 180) / dlon + 1 )
		iend = int ((maxlon + 180) / dlon + 1 )
		
		jstart = int ((minlat + 90) / dlat + 1 )
		jend = int ((maxlat + 90) / dlat + 1 )
		
		print*, minlon, maxlon, minlat, maxlat
		print*, istart, iend, jstart, jend
		
		
	end subroutine find_geos5_bounds


	subroutine get_model_idx_geos5(grid_xstart, grid_ystart, lat, lon, model_i, model_j)

		integer, intent(in) :: grid_xstart, grid_ystart
		real, intent(in) :: lat, lon
		integer, intent(inout) :: model_i, model_j

		real :: minLat, maxLat, minLon, maxLon
		real, parameter :: dlon = 5./16.
		real, parameter :: dlat = 0.25

		
		model_i = int ((lon + 180) / dlon + 1 ) - grid_xstart + 1
		model_j = int ((lat + 90) / dlat + 1 ) - grid_ystart + 1
		
	end subroutine get_model_idx_geos5

 subroutine read_int_array_g5(filedata, file_idx, sdsname, start, stride, edge, rank, dataholder, status) 
 
	integer, dimension(:), intent(in) :: filedata
	integer, intent(in) :: rank, file_idx
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	integer, dimension(*), intent(inout) :: dataholder
 
 
	integer :: file_id
	
	integer :: dataspace_id, memspace_id
	integer(HSIZE_T) :: edge1d(1)
	integer(HID_T) :: var_id, err_code

	real :: temp_fac(2)
	integer(HSIZE_T) :: offset(4), count(4), localstride(4), offset_mem(4)
	
	character(len=2) :: band_tag

	
	file_id = filedata(file_idx)


	call h5dopen_f(file_id, trim(sdsname), var_id, err_code)
	call h5dget_space_f(var_id, dataspace_id, err_code)

	offset(1:rank+1) = start
	if (rank == 2) then 
		offset_mem(1:2) = (/0,0/)
	else if (rank == 3) then
		offset_mem(1:3) = (/0,0,0/)
	endif
	
	count(1:rank+1) = edge
	localstride(1:rank+1) = stride
	
!	print*, trim(sdsname), offset(1:rank)
!	print*, trim(sdsname), edge(1:rank)

	call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset(1:rank+1), count(1:rank+1), err_code, localstride(1:rank+1))
	call h5screate_simple_f(rank, count(1:rank+1), memspace_id, err_code)
	call h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, offset_mem(1:rank), count(1:rank+1), err_code, localstride(1:rank+1))

	call h5dread_f(var_id, H5T_NATIVE_INTEGER, dataholder, count(1:rank+1), err_code, memspace_id, dataspace_id)
		
	call h5sclose_f(dataspace_id, err_code)
	call h5sclose_f(memspace_id, err_code)
		
	call h5dclose_f(var_id, err_code)
		
	status = 0

 end subroutine read_int_array_g5 




 subroutine read_float_array_g5(filedata, file_idx, sdsname, start, stride, edge, rank, dataholder, status) 
 
	integer, dimension(:), intent(in) :: filedata
	integer, intent(in) :: rank, file_idx
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real, dimension(*), intent(inout) :: dataholder
 
 
	integer :: file_id
	
	integer :: dataspace_id, memspace_id
	integer(HSIZE_T) :: edge1d(1)
	integer(HID_T) :: var_id, err_code

	real :: temp_fac(2)
	integer(HSIZE_T) :: offset(5), count(5), localstride(5), offset_mem(5)
	
	character(len=2) :: band_tag

	
	file_id = filedata(file_idx)


	call h5dopen_f(file_id, trim(sdsname), var_id, err_code)
	call h5dget_space_f(var_id, dataspace_id, err_code)

	print*, rank+1, start, trim(sdsname)

	offset(1:rank+1) = start
	if (rank == 2) then 
		offset_mem(1:2) = (/0,0/)
	else if (rank == 3) then
		offset_mem(1:3) = (/0,0,0/)
	else if (rank == 4) then 
		offset_mem(1:4) = (/0,0,0,0/)
	endif
	
	
	count(1:rank+1) = edge
	localstride(1:rank+1) = stride
	
	call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset(1:rank+1), count(1:rank+1), err_code, localstride(1:rank+1))
	call h5screate_simple_f(rank, count(1:rank+1), memspace_id, err_code)
	call h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, offset_mem(1:rank), count(1:rank+1), err_code, localstride(1:rank+1))

	call h5dread_f(var_id, H5T_NATIVE_REAL, dataholder, count(1:rank+1), err_code, memspace_id, dataspace_id)
		
	call h5sclose_f(dataspace_id, err_code)
	call h5sclose_f(memspace_id, err_code)
		
	call h5dclose_f(var_id, err_code)
		
	status = 0

 end subroutine read_float_array_g5 
 
 
  subroutine check_datasize(filedata, start, stride, edge, status)


   integer, dimension(:), intent(in)                   :: filedata
   integer, dimension (:), intent(inout) :: start, edge, stride
   integer, intent(out)                  :: status


	integer :: file_id, err_code, var_id, dataspace_id
	integer*8 :: dims(10), rank(10)
		
	call h5dopen_f(filedata(1),"PS", var_id, err_code)
	call h5dget_space_f(var_id, dataspace_id, err_code)
		
	call h5sget_simple_extent_dims_f(dataspace_id, dims, rank, err_code)

	call h5sclose_f(dataspace_id, err_code)
		
	call h5dclose_f(var_id, err_code)
		
	edge(1) = dims(1)
	edge(2) = dims(2)

	print*, dims(1:3)

 end subroutine check_datasize

 
 end module geos5_io_module
