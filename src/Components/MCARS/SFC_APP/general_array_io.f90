 module general_array_io
 
 implicit none
 
 include "hdf.f90"
 include "dffunc.f90"
 
 contains
 
 subroutine GetWidHt(filedata, sdsname, wid, ht)
 
 	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, intent(inout) :: wid, ht
	
	integer :: file_id, var_id, err_code, dummy, dims(2)
	character(len=200) :: dummy_name

	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
	
	err_code = sfginfo(var_id, dummy_name, dummy, dims, dummy, dummy)
	err_code = sfendacc(var_id)
 
 
	wid = dims(1)
	ht = dims(2)
 
 
 end subroutine GetWidHt
 
 subroutine read_byte_array(filedata, sdsname, start, stride, edge, dataholder, status) 
  
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	integer*1, dimension(:,:), intent(inout) :: dataholder
 
 
	integer :: file_id, var_id, err_code, attr_id
	integer*1 :: fill_val
	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
			
	fill_val = 0
	
	err_code = sfrdata(var_id, start, stride, edge, dataholder)
	
	err_code = sfendacc(var_id)
		
	
	status = 0

 end subroutine read_byte_array	 


 subroutine read_int16_array(filedata, sdsname, start, stride, edge, dataholder, status) 
 
 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	integer*2, dimension(:,:), intent(inout) :: dataholder
 
 
	integer :: file_id, var_id, err_code, attr_id
	integer*1 :: fill_val
	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
			
	fill_val = 0
	
	err_code = sfrdata(var_id, start, stride, edge, dataholder)
	
	err_code = sfendacc(var_id)
		
	
	status = 0

 end subroutine read_int16_array	 

 
 
 
  subroutine read_int_array(filedata, sdsname, start, stride, edge, use_offset, dataholder, status) 
 
	integer, parameter :: fillvalue_real = -999
 
	integer, dimension(:), intent(in) :: filedata
	logical, intent(in) :: use_offset
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real, dimension(:,:), intent(inout) :: dataholder
 
 
	integer :: file_id, var_id, err_code, attr_id
	real*8 :: scale_factor, offset
	integer*2 :: fill_val
	integer*2, dimension(:,:), allocatable :: temp
	integer:: i, j
	
	allocate(temp(edge(1), edge(2)))
	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
		
	attr_id = sffattr(var_id, "scale_factor")
	err_code = sfrattr(var_id, attr_id, scale_factor)

	if (use_offset) then 
		attr_id = sffattr(var_id, "add_offset")
		err_code = sfrattr(var_id, attr_id, offset)
	else
		offset = 0.
	endif
		
	attr_id = sffattr(var_id, "_FillValue")
	err_code = sfrattr(var_id, attr_id, fill_val)

	err_code = sfrdata(var_id, start, stride, edge, temp)

	
	err_code = sfendacc(var_id)

	do i=1, edge(1)
	  do j=1, edge(2)

		if (temp(i,j) /= fill_val) then 
			dataholder(i,j) = (temp(i,j) - offset)*scale_factor
		else
			dataholder(i,j) = fillvalue_real
  		endif

	 end do
       end do
	   
	deallocate(temp)
	
	status = 0

 end subroutine read_int_array	 

 subroutine read_float_array(filedata, sdsname, start, stride, edge, dataholder, status) 
 
 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real, dimension(*), intent(inout) :: dataholder
 
 
	integer :: file_id, var_id, err_code, attr_id
	
	integer :: rank, dims(10), dtype, attrs
	character(len=100) :: name
	integer :: i, j
	
	real :: fill_val
	
	file_id = filedata(1)

	var_id = sfselect(file_id, sfn2index(file_id, trim(sdsname)))
		
	attr_id = sffattr(var_id, "_FillValue")
	err_code = sfrattr(var_id, attr_id, fill_val)
	
	err_code = sfrdata(var_id, start, stride, edge, dataholder)
	
	err_code = sfendacc(var_id)
	
!	where (dataholder == fill_val)
!		dataholder = fillvalue_real
!	end where
	
	status = 0

 end subroutine read_float_array	 

 subroutine read_double_array(filedata, sdsname, start, stride, edge, dataholder, status) 
 
 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real*8, dimension(*), intent(inout) :: dataholder
 
 
	integer :: file_id, var_id, err_code, attr_id
	
	integer :: rank, dims(10), dtype, attrs
	character(len=100) :: name
	integer :: i, j
	
	real :: fill_val
	
	file_id = filedata(1)

	var_id = sfselect(file_id, sfn2index(file_id, trim(sdsname)))
		
	attr_id = sffattr(var_id, "_FillValue")
	err_code = sfrattr(var_id, attr_id, fill_val)
	
	err_code = sfrdata(var_id, start, stride, edge, dataholder)
	
	err_code = sfendacc(var_id)
	
!	where (dataholder == fill_val)
!		dataholder = fillvalue_real
!	end where
	
	status = 0

 end subroutine read_double_array	 




subroutine write_scaled_array(filedata, sdsname, start, stride, edge, data_to_write, status) 

!   Subroutine added to covert OD floating point SDS's to scaled integers (to reduce output file size)
!       GTA 08-06-08

 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real, dimension(:,:), intent(in) :: data_to_write
	integer*2, dimension(:,:), allocatable :: scaled_data_to_write
	integer :: xsize, ysize
 
 
	integer :: file_id, var_id, err_code, attr_id, dim_id
	integer*2 :: fill_val
	
	real*8 :: scale_factor, offset
		
!   find dimesnions of input data (data to be scaled), then allocate array for output SDS

	xsize = size(data_to_write,1)
	ysize = size(data_to_write,2)

	allocate(scaled_data_to_write(xsize,ysize))
	  
!   Retrieve appropriate scale factor and offset (set previously in CR step)

	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
	if (var_id == -1) var_id = sfcreate(file_id, trim(sdsname), DFNT_INT16, 2, edge(1:2))
			
	attr_id = sffattr(var_id, "scale_factor")
	err_code = sfrattr(var_id, attr_id, scale_factor)

	attr_id = sffattr(var_id, "add_offset")
	err_code = sfrattr(var_id, attr_id, offset)
	
	attr_id = sffattr(var_id, "_FIllValue")
	err_code = sfrattr(var_id, attr_id, fill_val)
	
	if (err_code == -1) fill_val=-999  ! fill_val should have the proper value here but may so force it if necessary

!   apply scale factor and offset to input data, but only if data not fill value

	
	scaled_data_to_write = fill_val

	where (data_to_write /= fill_val*1.)
		scaled_data_to_write = NINT((data_to_write / scale_factor + offset)) 
	end where


	

	
!   Define dimension names for SDS (all scaled SDS's are two-dimensional)

	dim_id = sfdimid(var_id, 0)
	err_code=sfsdmname(dim_id, 'NumberOfPixels')

	dim_id = sfdimid(var_id, 1)
	err_code=sfsdmname(dim_id, 'NumberOfLines')

 
	err_code = sfwdata(var_id, start, stride, edge, scaled_data_to_write)
	
	err_code = sfendacc(var_id)
	
	status = 0

 end subroutine write_scaled_array	 



 subroutine write_int_array(filedata, sdsname, start, stride, edge, data_to_write, status) 
 
 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	integer*2, dimension(:,:), intent(inout) :: data_to_write
 
 
	integer :: file_id, var_id, err_code, attr_id, dim_id
	

	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
	if (var_id == -1) var_id = sfcreate(file_id, trim(sdsname), DFNT_INT16, 2, edge(1:2))
	
!   Define dimension names for two-dimensional SDS's and for Cloud Mask and Quality Assurance SDS
!     (all 2-D SDS's have same dimensions, CM and QA SDS's are 3-D) GAT 8-6-08 

	err_code = sfwdata(var_id, start, stride, edge, data_to_write)
	
	err_code = sfendacc(var_id)
	
	status = 0

 end subroutine write_int_array

 subroutine write_byte_array(filedata, sdsname, start, stride, edge, data_to_write, status) 
 
 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	integer*1, dimension(*), intent(inout) :: data_to_write
 
 
	integer :: file_id, var_id, err_code, attr_id, dim_id
	

	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
	if (var_id == -1) var_id = sfcreate(file_id, trim(sdsname), DFNT_INT8, 2, edge(1:2))
	

	err_code = sfwdata(var_id, start, stride, edge, data_to_write)
	
	err_code = sfendacc(var_id)
	
	status = 0

 end subroutine write_byte_array	 



 
 subroutine write_float_array(filedata, sdsname, start, stride, edge, data_to_write, status) 
 
 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real, dimension(:,:), intent(inout) :: data_to_write
 
 
	integer :: file_id, var_id, err_code, attr_id, dim_id
	
	real :: fill_val
	
	fill_val = -999.
	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
	if (var_id == -1) &
		var_id = sfcreate(file_id, trim(sdsname), DFNT_FLOAT, 2, edge(1:2))
	
!   Define dimension names for one or two-dimensional SDS's  GTA 12-18-08

	dim_id = sfdimid(var_id, 1)

	if (dim_id == -1) then             !   If one dimension name then assign name one such varible)

	   dim_id = sfdimid(var_id, 0)
	   err_code=sfsdmname(dim_id, 'NumberOfBands')
	   
	else        !   Define dimension names for 2-D SDS (most float SDS's are two-dimensional) 

	   dim_id = sfdimid(var_id, 0)
	   err_code=sfsdmname(dim_id, 'NumberOfPixels')

	   dim_id = sfdimid(var_id, 1)
	   err_code=sfsdmname(dim_id, 'NumberOfLines')
	 
	endif
		
		
	attr_id = sffattr(var_id, "_FillValue")
	if (attr_id == -1) then
		attr_id = sfsattr(var_id, "_FillValue", DFNT_FLOAT, 1, fill_val)
	endif
	
	err_code = sfwdata(var_id, start, stride, edge, data_to_write)
	
	err_code = sfendacc(var_id)
	
	status = 0

 end subroutine write_float_array	 


 
 end module general_array_io
