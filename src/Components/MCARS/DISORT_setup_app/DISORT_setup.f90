! This program is a setup for a DISORT granule test run that includes all the profile and cloud information. 
! This code will be thrown away once we have GEOS-5 interface working 


! To-do: 

! 1. Process GDAS so that it's at the correct model levels
! 2. take MOD06 cloud information to fill in the cloud positions and information
! 3. output everything to disk where a complete setup exists for every pixel that we would attempt


program DISORT_setup
	
	use general_array_io, only: GetWidHt, read_int_array, read_byte_array, read_float_array, write_float_array, read_int16_array
	use global_model_grids
	use disort_ancillary
	use corrK_profile
	use DISORT_solver
	use modis_simulator

	use geos5_io_module
	
	use hdf5
	
	implicit none
		
	include "hdf.f90"
	include "dffunc.f90"
	
	character(len=200) :: parfile_name, out_prefix
	character(len=200) :: MOD03_name, output_name, albedo_name
	character(len=200) :: geos5_name, geos5_aerosol_name
	character(len=10) :: arg_str

	
	integer :: modis_wid, modis_ht, start(2), stride(2), edge(2)
	integer :: year, day_of_year, file_id(1), err_code, i, fid(1), edg(2)

	real, dimension(:,:), allocatable :: latitude, longitude
	real, dimension(:,:), allocatable :: CTH, COT, CER
	integer*1, dimension(:,:), allocatable :: CPOP
	
	
	integer :: num_Iter, lat_start, lat_end, j, k
	integer :: start_time, end_time, crate, cmax
	character(len=10) :: platform_name
	logical :: useoffset
	real*8, parameter :: d2r = 0.017453292519943295d0
	integer*1 :: leap_offset
	integer*1, dimension(:,:), allocatable :: ocean_mask
	
! DISORT business
	real, dimension(:,:), allocatable :: sza, vza, phi
	integer, parameter :: nch = 25
	real*8 :: sza_pass, vza_pass, phi_pass
	real*8 :: Ts, Tt 
	real*8 :: rad(nch)
	

! GEOS-5 QUIERY

! GEOS-5 QUIERY

	integer :: grid_xsize, grid_xstart, iend
	integer :: grid_ysize, grid_ystart
		
	real, dimension(:,:), allocatable :: met_temp1, met_temp2, met_temp3, met_temp4
	real :: xx(2), yy1(2), yy2(2)
	real :: met_time(4), frac_time
   
   real, dimension(:,:), allocatable :: Tsfc, colO3, Ttoa
	integer :: MYTIME, MYMONTH
	integer :: model_i, model_j
	real :: snow, ice 
	
	integer :: gfile_id(2)

! END GEOS-5 QUIERY	

! MP VARIABLES

	integer :: mpistart, mpjstart, mpiend, mpjend, mpisize, mpjsize
	integer :: mpthread, mp_num_threads
	character(len=3) :: mptag

	integer :: jstart, jend

! END MP VARIABLES

! TEST variables
	integer*2 :: test_albedo(25)

! end TEST variables

! OUTPUT YAY!
	real, dimension(:,:,:), allocatable :: final_answer
	integer :: start3(3), stride3(3), edge3(3)
	
! TEST settings
! 1,2,3,4,5,6,7,8,9,17,18,19,26
! 20,22,27,28, 29, 31, 32, 33, 34, 35, 36
 
	test_albedo = (/27, 288, 17, 37, 252, 252, 146, 54, 14, 22, 283, 280, 280, 38, 38, 216, &
				25, 24, 29, 31, 25, 22, 22, 22, 22 /)

        
! end TEST settings	
	
	

	disort_heights(1:10) = (/80., 60., 50., 45., 40., 35., 30., 25., 20., 18. /)
	disort_heights(11:27) = (/ (17. - i, i=1, 17) /)
		
	
! get the arguments	: parameter file name, number of threads and thread number
	call getarg(1, parfile_name)

        call getarg(2, arg_str)
        read(arg_str, *) jstart
        call getarg(3, arg_str)
        read(arg_str, *) jend

        call getarg(4, arg_str)
        read(arg_str, *) mp_num_threads
	call getarg(5, arg_str)
	read(arg_str, *) mpthread

	
	if (mpthread > mp_num_threads .or. mpthread < 1) then
		write(*,'(a,i4)') "Invalid thread number. Valid thread numbers are: 1 through", mp_num_threads
		print*, "Stopping...."
		stop
	endif
	
! this is so the domain split works. Prund.pl doesn't like thread number 0 for some reason, 
! but the domain splitting relies on thread number running from 0 to nthreads-1
	mpthread = mpthread - 1

	if(mpthread < 10) then 
		write(mptag, '(a2,i1)') "00", mpthread
	else if (mpthread >= 10 .and. mpthread < 100) then 
		write(mptag, '(a1,i2)') "0", mpthread
	else
		write(mptag, '(i3)') mpthread
	endif	

! parse out the information in parameter file	
	call read_par_file(parfile_name)

	if (jstart < 10) then 
		write(arg_str, '(a3i1)') "000", jstart
	else if (jstart >= 10 .and. jstart < 100) then 
		write(arg_str, '(a2i2)') "00", jstart
	else if (jstart >= 100 .and. jstart < 1000) then 
		write(arg_str, '(a1,i3)') "0", jstart
	else
		write(arg_str, '(i4)') jstart
	endif

	
	output_name = trim(out_prefix) // ".DISORT." // mptag // "." // trim(arg_str) // ".hdf"
	print*, trim(output_name)

! we will process in thin strips of data, shouldn't be a problem I suppose. 
	file_id(1) = sfstart(MOD03_name, DFACC_READ)
	call GetWidHt(file_id, "Latitude", modis_wid, modis_ht)
	
	mpisize = nint( modis_wid*1.0 / mp_num_threads) 
	mpjsize = modis_ht !nint( modis_ht*1.0 / mp_num_threads) 

	mpistart = mpthread * mpisize + 1
	mpjstart = 1 !mpthread * mpjsize + 1
	
	mpiend = mpistart + mpisize - 1
	mpjend = modis_ht  !mpjstart + mpjsize - 1
			
	if (mpiend > modis_wid .or. mpthread == mp_num_threads-1 ) then
		mpiend = modis_wid
		mpisize = mpiend - mpistart + 1
	endif
	if (mpjend > modis_ht) mpjend = modis_ht

	if (mpistart > modis_wid) then 
		print*, "out of bounds. Stopping..."
		stop
	endif

	print*, mpthread, "|", mpistart, mpiend, mpisize, mpjsize

	allocate(sza(mpisize, mpjsize), vza(mpisize, mpjsize))
			 
	start = (/ mpistart-1, mpjstart-1 /)
	stride = 1
	edge = (/ mpisize, mpjsize /)

	call read_int_array(file_id, "SolarZenith", start, stride, edge, .false., sza, err_code)
	call read_int_array(file_id, "SensorZenith", start, stride, edge, .false., vza, err_code)
	
	err_code = sfend(file_id(1))


! now we've bounded the model, we can start reading data from it. Things we really need here are: 
	allocate(granule_info(mpisize, mpjsize))


	call h5open_f(err_code)

	call obtain_geos5(mpisize, mpjsize, mpistart-1, mpjstart-1, geos5_name, &
								geos5_aerosol_name, granule_info)			
				
	call h5close_f(err_code)
	
! now we need to read in the surface albedo from that file we made

	allocate(surface_albedo(nch, mpisize, mpjsize))
	allocate(phi(mpisize, mpjsize))
	allocate(ocean_mask(mpisize, mpjsize))
		
	file_id(1) = sfstart(albedo_name, DFACC_READ)
	call read_surface_albedo(file_id, phi, ocean_mask, (/mpistart-1, mpjstart-1/), (/mpisize, mpjsize/))
	err_code = sfend(file_id(1))

! GEOS-5 QUIERY 

	allocate(final_answer(nch, mpisize, mpjsize))

	call init_corrK

	call init_legendre

	call init_solver

	call system_clock(start_time, crate, cmax)


!	do j=1, mpjsize
	do j=jstart, jend
		if ( mod(j,10) == 0) &
		print*, j
		do i= 1, mpisize

! Now we will call the Corr-K code with the profile that's been made here. 
			call corrK_driver(i,j)
			
! Now we will set up the stuff for the DISORT driver
			sza_pass = sza(i,j)
			vza_pass = vza(i,j)
			phi_pass = phi(i,j)
	
			Ts = granule_info(i,j)%Ts
			Tt = granule_info(i,j)%temperature(1)
			

			if (ocean_mask(i,j) == 0 ) then 
				surface_albedo(1:13, i, j) = nint ( 50*granule_info(i,j)%albedo_factor(1:13) + &
								surface_albedo(1:13,i,j)*(1.-granule_info(i,j)%albedo_factor(1:13)))
								
                surface_albedo(14, i, j) = nint ( 25*granule_info(i,j)%albedo_factor(14) + &
                                                   surface_albedo(14,i,j)*(1.-granule_info(i,j)%albedo_factor(14)))
                                                   
                surface_albedo(16, i, j) = nint ( 50*granule_info(i,j)%albedo_factor(16) + &
                                                   surface_albedo(16,i,j)*(1.-granule_info(i,j)%albedo_factor(16)))

			endif
					
										
			call init_disort(granule_info(i,j), surface_albedo(:,i,j), &
					sza_pass, vza_pass, phi_pass)

!			call init_disort(granule_info(i,j), test_albedo, &
!					sza_pass, vza_pass, phi_pass)

			call run_disort(Ts, Tt, rad )
			final_answer(:,i,j) = rad(:)

!			write(*, '(i4,i3,i4,25f9.5)') mpthread, i, j, final_answer(:,i,j)
			call cleanup_disort
!			stop

		end do 
	end do


	call system_clock(end_time, crate, cmax)
	
	print*, "Time: ", (end_time - start_time) / (crate*1.0) , "sec"


	call cleanup_corrK


	file_id(1) = sfstart(output_name, DFACC_CREATE)
	
	start3 = 0
	stride3 = 1
	edge3 = (/ nch, mpisize, mpjsize /)
	call write_float_array(file_id, "Refl/Rad", start3, stride3, edge3, final_answer, err_code)
	
	err_code = sfend(file_id(1))

	do j=1, mpjsize
		do i=1, mpisize
				if (allocated( granule_info(i,j)%pMom)) deallocate(granule_info(i,j)%pMom)
		end do
	end do
	deallocate(granule_info)
	deallocate(surface_albedo, phi, ocean_mask)
	
	deallocate(final_answer)

	deallocate(sza, vza)
				
contains 


	subroutine read_par_file(fname)	

		character*(*) :: fname
		
		
		character(len=200) :: temp, dummy, input_dir, output_dir
		character(len=4) :: mydate, mt
		character(len=3) :: doy
		integer :: i
		
		
		open (unit=10, file = fname)
		
		read(10, '(a)') dummy
		read(10, *) year, doy, mt, mydate
		read(mt, '(i4)') MYTIME
		read(doy, '(i3)') day_of_year
		read(mydate, '(i2)') MYMONTH
		
		write(input_dir, '(a,i4,a3,a1,a4,a1)') "inputs-", year, doy, "." , mt, "/"
		write(output_dir, '(a,i4,a3,a1,a4,a1)') "outputs-", year, doy, "." , mt, "/"
				
		read(10, '(a)') dummy
		read(10, '(a)') temp

		MOD03_name = trim(input_dir) // trim(temp) 
		
		read(10, '(a)') dummy
		read(10, '(a)') temp

		albedo_name = trim(input_dir) // trim(temp)

		read(10, '(a)') dummy
		read(10, '(a)') temp
		
		out_prefix = trim(output_dir) // trim(temp)

		read(10, '(a)') dummy
		read(10, '(a)') temp
		geos5_name = trim(input_dir) // trim(temp)

		! skip albedo path
		read(10, '(a)') dummy
		read(10, '(a)') temp

		read(10, '(a)') dummy
		read(10, '(a)') temp
		if (temp == "none") then
			AEROSOL_yes = .false.
		else
			AEROSOL_yes = .true.
		endif
		geos5_aerosol_name = trim(input_dir) // trim(temp)


		read(10, '(a)') dummy
		read(10, '(a)') ice_coef_name
		read(10, '(a)') liquid_coef_name

		if (AEROSOL_yes) then 
			read(10, '(a)') dummy
			read(10, '(a)') rcfile
		endif


		close(10)	
	
	
	end subroutine read_par_file		
		
		
		

end program DISORT_setup
