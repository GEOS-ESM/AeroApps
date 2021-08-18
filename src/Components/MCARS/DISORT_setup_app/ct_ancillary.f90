module ct_ancillary

 implicit none
 private

 public :: read_ancillary_grids
 
 
contains

 subroutine get_bilinear_idx(lat, lon, i, j, idx_x, idx_y)
 
	real, intent(in) :: lat, lon
	integer, intent(in ) :: i,j
	integer, intent(inout) :: idx_x(2), idx_y(2)
	
	 real :: x,y, x0, dx, y0, dy
	
      x = min( max( lon,  -179.99 ), 179.99 )
      if( x > -999. .and. x < 0.0 ) x = lon+ 360.0

      y = min( max( lat, -89.99 ), 89.99 )
      y0 = 90.0
      dy = -1.0

		if ( x < (i*1.) ) then !+0.5) then 
			idx_x(1) = i-1
			idx_x(2) = i
		else
			idx_x(1) = i
			idx_x(2) = i+1
		endif

		if (idx_x(1) < 0 .or. idx_x(1) > 359) idx_x(1) = 0
		if (idx_x(2) < 0 .or. idx_x(2) > 359) idx_x(2) = 0
		
		if ( (y-y0)/dy < (j*1.) ) then !+0.5) then 
			idx_y(1) = j-1
			idx_y(2) = j
		else
			idx_y(1) = j
			idx_y(2) = j+1
		endif

		if (idx_y(1) < 0 ) idx_y(1) = 0
		if (idx_y(1) > 180) idx_y(1) = 180

		if (idx_y(2) < 0 ) idx_y(2) = 0
		if (idx_y(2) > 180) idx_y(2) = 180


		idx_x = idx_x + 1
		idx_y = idx_y + 1
 
 end subroutine get_bilinear_idx


real function get_W(RH, T, P)

! After Rogers and Yau. 

	real*8, intent(in) :: RH, T, P

	real, parameter :: c = 4187.
	real, parameter :: cpv = 1870.
	real, parameter :: L0 = 2.501e6
	real, parameter :: T0 = 273.
	
	real :: L
	real :: esat
	real :: ws
	real :: W
	
	if (T < 235.) then 
		L = 2.83e6   ! latent heat of ice, doesn't change very much
	else
		L = L0 - (c - cpv) * (T - T0) ! latent heat of water
	endif
	
	esat = 6.11657*exp((L/461.51)*(1/273. - 1/T) ) ! saturation vapor pressure in mb

	ws = esat * 0.622 / (P - esat) ! saturation mixing ratio in kg/kg

	W = RH / 100. * ws ! mixing ratio in kg/kg
	
	get_W = W * 1000. ! mixing ratio in g/kg
	


end function get_W

real function get_W_from_Td( Td, T, P)

	real*8, intent(in) :: Td, T, P

	real, parameter :: c = 4187.
	real, parameter :: cpv = 1870.
	real, parameter :: L0 = 2.501e6
	real, parameter :: T0 = 273.
	
	real :: L
	real :: esat, esatW, RH
	real :: ws
	real :: W
	
	if (T < 235.) then 
		L = 2.83e6   ! latent heat of ice, doesn't change very much
	else
		L = L0 - (c - cpv) * (T - T0) ! latent heat of water
	endif
	
	esat = 6.11657*exp((L/461.51)*(1/273. - 1/T) ) ! saturation vapor pressure in mb
	esatW = 6.11657*exp((L/461.51)*(1/273. - 1/Td) ) ! saturation vapor pressure in mb

	RH = esat / esatW * 100.

	ws = esat * 0.622 / (P - esat) ! saturation mixing ratio in kg/kg

	W = RH / 100. * ws ! mixing ratio in kg/kg
	
	get_W_from_Td = W * 1000. ! mixing ratio in g/kg

end function get_W_from_Td


integer function find_trop(temp_prof, p_prof, sfc_lev)

	real, dimension(:), intent(in) :: temp_prof, p_prof
	integer*1, intent(in) :: sfc_lev

	real :: xmin
	integer :: ilev, imin
	real, parameter :: PTOP = 100.
	
	xmin = 999999.
	imin = 1
	
	do ilev = 1, sfc_lev-4 !-5+1 A C adjustment
		if (temp_prof(ilev) < xmin .and. p_prof(ilev) >= PTOP) then
			xmin = temp_prof(ilev)
			imin = ilev
		endif
	end do


! don't allow trop height > 400 mb (level 71)
	find_trop = -99
	do ilev = imin, 71
	
		if (temp_prof(ilev-1) >= temp_prof(ilev) .and. &
			temp_prof(ilev+1) > temp_prof(ilev)) then 
				find_trop = ilev
				exit
		endif
	end do

	if (find_trop == -99) find_trop = imin

#ifdef MAS_INST	
	if (find_trop < 36) find_trop = 36
#endif

end function find_trop





subroutine read_ancillary_grids( model_info, lat, lon, ncepgdas_name,ncepgdas_name2,  &
                                 sfctmp , MYTIME, MYMONTH )
  use global_model_grids
  use modis_numerical_module, only: linearinterpolation, bilinear_interpolation
  implicit none

	include "hdf.f90"
	include "dffunc.f90"

!  include 'Atmos_AncData.f90.inc'

  !-----------------------------------------------------------------------
  ! !f90
  !
  ! !description:
  !      retrieve ancillary data items for a given set of 
  !      latitude and longitude pairs.
  !
  ! !input parameters:
  !      lat       latitude (degrees, -90s to +90.0n)
  !      lon       longitude (degrees, -180w to +180e, greenwich=0)
  !
  ! !output parameters:
  !      pres      array of pressure levels (hpa)
  !      temp      array of atmospheric temperatures (k) at pres(0:15)
  !      mixr      array of water vapor mixing ratios (g/kg) at pres(0:15)
  !      land      land mask (0=water, 1=land)
  !      sfctmp    surface temperature (k)
  !      prmsl     pressure (hpa) at mean sea level
  !      pwat      precipitable water (g/cm**2)
  !      ugrd      surface wind u component (m/s)
  !      vgrd      surface wind v component (m/s)
  !      Ozone     total column Ozone (dobson units)
  !      icec      ice concentration (fraction)
  !
  ! !notes
  !      Modified from the MOD_PR06OD Collection 4 ReadNCEP.f algorithm written
  !      by Liam Gumley and Jason Li.

  real, intent(in)     :: lat(:,:),lon(:,:)
  character(*), intent(in)     :: ncepgdas_name, ncepgdas_name2
  type(ancillary_type), dimension(:,:), intent(inout) :: model_info

  real, intent(inout)     ::   sfctmp(:,:)

  ! ... parameter definitions
  integer     npoints_x,        npoints_y  
  parameter ( npoints_x = 360,  npoints_y = 180 )

  real        missing
  parameter ( missing = -999.0 )

  ! ... declaration of local variables and  arrays
  character*8    esdt_name_ncep, esdt_name_ozone, esdt_name
  character*7    esdt_name_ice
  character*160  errmsg
  character*13   ncepgdastemp_name, ozonetemp_name
  character*12   icetemp_name
  character*400  temp_output_name, temp_input_name
	integer, intent(inout) :: MYTIME, MYMONTH



  integer header( 0:7 ), i, ios, j, k, level, lun, pcfnum, reclen, status, ii, jj

  integer data_xsize, data_ysize, grid_i, grid_j, kstart

  real  ::    satmix, x, x0, xlon, dx, y, y0, dy

        real :: yy2(2), yy(2), xx(2)
		integer :: idx_x(2), idx_y(2)

#ifdef USE_GDAS
	integer, parameter :: num_gdas_vars = 60
	real :: met_grid( 0:359, 0:180, 0:num_gdas_vars-1 ), met_grid2( 0:359, 0:180, 0:num_gdas_vars-1 )
#else
	integer, parameter :: num_gdas_vars = 52
	real :: met_grid( 0:359, 0:180, 0:num_gdas_vars-1 ), met_grid2( 0:359, 0:180, 0:num_gdas_vars-1 )
#endif

	integer met_year, met_month, met_day, met_hour, &
          ice_year, ice_month, ice_day, ice_hour, &
          ozn_year, ozn_month, ozn_day, ozn_hour

	integer met_date(4), met_date2(4), ice_date(4),  ozn_date(4)
	
	real*8 :: o3_temp(6)
	integer :: o3_start

	logical met_success, ice_success, ozn_success


	real :: Ts2m, W2m, Pmsl, Ts

  ! ... declaration of external functions
  integer modis_grib_driver 
  external modis_grib_driver 

	integer make_profile_101
	external make_profile_101
	
	integer profile_to_101
	external profile_to_101

	integer height_profile
	external height_profile
	
	external climoz_101
	
	integer adjo3
	external adjo3

	real*8 :: model_lat
	real*8 :: heights(model_levels)

#ifdef USE_GDAS
	integer, parameter :: model_coarse = 27
	integer, parameter :: wind_start = 47
	integer, parameter :: rh_start = 6
#else
	integer, parameter :: model_coarse = 21
	integer, parameter :: wind_start = 40	
	integer, parameter :: rh_start = 1
#endif

	real :: a_factor
	integer*1 :: sfc_level


    real*8 :: p(model_coarse), p_source(model_coarse)

	real*8 :: mixr(model_coarse), temp(model_coarse)
	
	real*8 :: pressures(model_levels), temp_hires(model_levels), mixr_hires(model_levels)
	real*8 :: pressures_source(model_levels)

	real :: Ts_forint(2,2), lat_set(2), lon_set(2)

	integer :: file_id, var_id, err_code, start(2), stride(2), edge(2)

#ifdef USE_GDAS
	p_source = (/ 10., 20., 30., 50., 70., 100., 150., 200., 250., 300., 350., 400., 450., 500., &
			550., 600., 650., 700., 750., 800., 850., 900., 925., 950., 975., 1000., 1100. /)
#else
	p_source = (/ 10., 20., 30., 50., 70., 100., 150., 200., 250., 300.,       400.,       500., &
			      600.,       700.,       800., 850., 900., 925., 950.,       1000., 1100. /)

#endif

  !-----------------------------------------------------------------------
  !      begin executable code 
  !-----------------------------------------------------------------------
  ! ... read input data files if this is the first call

	lun = 222


  data_xsize = size(lat,1)
  data_ysize = size(lat,2)

  if(.not. grids_are_read) then
	grids_are_read = .true. 
	! ..... set data ingest success/fail flags
	met_success  = .false.
	ice_success  = .false.
	ozn_success  = .false.

	do i = 1, 4
		met_date( i )  = int( missing )
		ice_date( i )  = int( missing )
		ozn_date( i )  = int( missing )
	end do


  !-----------------------------------------------------------------------
  !      get ncep meteorological data
  !-----------------------------------------------------------------------
  ! ...   unpack grib met file and write to binary file
	errmsg    = ' '
	ESDT_name = 'GDAS_0ZF'
  
	temp_input_name = trim(ncepgdas_name) // char(0)
	temp_output_name = trim(ncepgdas_name) // ".bin" // char(0)
  
	status    = modis_grib_driver( temp_input_name, &
                                 ESDT_name, errmsg, &
                                 met_year, met_month, &
                                 met_day, met_hour )



	if ( status /= success ) then
		status = failure
		print* , 'error reported from grib driver, ncep read, ', &
                                'read_ancillary_grids'
		stop
    ! ....  open unpacked met file
	else

		reclen = 360*181*num_gdas_vars*4

		open (unit = lun, file = temp_output_name, form = 'unformatted', &
	       access = 'direct', status = 'old', &
             recl = reclen, iostat = ios )

		if ( status /= 0 ) then
			status = 2
			print*, 'error reported from openr, ncep read ', &
                                'read_ancillary_grids'
			stop
		else 
		!  read the unpacked met file
			if ( status == 0 ) then
				ios= 0
				read( lun, rec = 1, iostat = ios ) met_grid


				if ( ios /= 0 ) then
					level = 2
				else
					met_success   = .true.
					met_date( 1 ) = met_year
					met_date( 2 ) = met_month
					met_date( 3 ) = met_day
					met_date( 4 ) = met_hour * 100
					
				endif

				close(lun)
			endif
			
		endif
		
	endif

	errmsg    = ' '
	ESDT_name = 'GDAS_0ZF'
  
	temp_input_name = trim(ncepgdas_name2) // char(0)
	temp_output_name = trim(ncepgdas_name2) // ".bin" // char(0)
  
	status    = modis_grib_driver( temp_input_name, &
                                 ESDT_name, errmsg, &
                                 met_year, met_month, &
                                 met_day, met_hour )
 
	if ( status /= success ) then
		status = failure
		print*, 'error reported from grib driver, ncep read, Check PCF file', &
                                'read_ancillary_grids'

		stop
    ! ....  open unpacked met file
	else

		reclen = 360*181*num_gdas_vars*4
		open (unit = lun, file = temp_output_name, form = 'unformatted', &
	       access = 'direct', status = 'old', &
             recl = reclen, iostat = ios )

		if ( status /= 0 ) then
			status = 2
			print*, 'error reported from openr, ncep read, Check PCF file', &
                                'read_ancillary_grids'
			stop
		else 
		!  read the unpacked met file
			if ( status == 0 ) then
				ios= 0
				read( lun, rec = 1, iostat = ios ) met_grid2

				if ( ios /= 0 ) then
					level = 2
				else
					met_success   = .true.
					met_date2( 1 ) = met_year
					met_date2( 2 ) = met_month
					met_date2( 3 ) = met_day
					met_date2( 4 ) = met_hour * 100
					if (met_date2(4) == 0) met_date2(4) = 2400
				endif

				close(lun)
			endif
			
		endif
		
	endif

! capture month of the granule
	MYMONTH = met_date(2)
	


        xx(1) = met_date(4) * 1.0
        xx(2) = met_date2(4) * 1.0


! calculate the pressure levels for the 101-level profile, courtesy UW-Madison
    status = make_profile_101(pressures_source)

	do j=1, grid_ysize
	
		jj = j-1

		model_lat = 89.5 - jj*1.0
		
		call climoz_101( model_lat, MYMONTH, model_info(1,j)%o3_profile)

		do i=1, grid_xsize
			ii = i-1

			if (i > 1) &
				model_info(i,j)%o3_profile = model_info(1,j)%o3_profile

			p = p_source
			pressures = pressures_source
    
			
			yy(1) = met_grid(ii,jj,wind_start)
			yy(2) = met_grid2(ii,jj, wind_start) 

			yy2(1) = met_grid(ii,jj,wind_start+1)
			yy2(2) = met_grid2(ii,jj,wind_start+1)       

			model_info(i,j)%wind_speed = sqrt( linearinterpolation( xx, yy, MYTIME*1.) **2 + &
									linearinterpolation( xx, yy2, MYTIME*1.) **2 )
									
									
			kstart = 0
			do k = 1, model_coarse-1
				temp(k) = linearinterpolation( xx, &
									(/ met_grid(ii,jj,kstart), met_grid2(ii,jj, kstart) /), MYTIME*1.) !met_grid( i, j, k )
				kstart = kstart + 1
			end do


			kstart = model_coarse-1
			do k = rh_start, model_coarse-1
				mixr(k) = linearinterpolation( xx, &
						(/ met_grid(ii,jj, kstart), met_grid2(ii,jj, kstart) /), MYTIME*1.)  !met_grid( i, j, k + 16 )
				if (mixr(k) < 0.) mixr(k) = 0.
				kstart = kstart + 1 
			end do

!     convert relative humidity profile (%) to mixing ratio (g/kg)

			do k = rh_start, model_coarse-1
				mixr(k) = get_W (mixr(k), temp(k), p(k))
			end do



!     extrapolate mixing ratio profile from 100 hPa to 10 hPa
#ifdef USE_GDAS
			do k = 1, 5
				mixr(k) = max( mixr(6), 0.003d0 ) * ( p( k ) / 100.0 )**3 ! was 21 
				mixr(k) = max( mixr(k), 0.003d0 )
			end do
#endif


            model_info(i,j)%Ps = linearinterpolation( xx, (/ met_grid(ii,jj,wind_start+2), &
													met_grid2(ii,jj, wind_start+2) /), MYTIME*1.) * 0.01 !met_grid( i, j, 76 )
			
! get the surface parameters and convert the RH at 2m to mixing ratio
			Ts2m = linearinterpolation( xx, (/ met_grid(ii,jj,wind_start+3), &
											met_grid2(ii,jj, wind_start+3) /), MYTIME*1.)
			
! if we have ECMWF this W2m is not an RH, but a dew point temperature.
			W2m = linearinterpolation( xx, (/ met_grid(ii,jj,wind_start+4),&
											met_grid2(ii,jj, wind_start+4) /), MYTIME*1.)

#ifdef USE_GDAS
			W2m = get_W(W2m*1.0d0, Ts2m*1.0d0, model_info(i,j)%Ps*1.0d0)
			Pmsl = linearinterpolation( xx, (/ met_grid(ii,jj,wind_start+5), &
											met_grid2(ii,jj, wind_start+5) /), MYTIME*1.)* 0.01
#else
			W2m = get_W_from_Td(W2m*1.0d0, Ts2m*1.0d0, model_info(i,j)%Ps*1.0d0)
			Pmsl = model_info(i,j)%Ps
#endif

			if (model_info(i,j)%Ps > 0. .and. p(model_coarse-1) <= model_info(i,j)%Ps) then 
				model_info(i,j)%surface_level = model_coarse
			else
				model_info(i,j)%surface_level = 0
			endif

			model_info(i,j)%Ts = Ts2m
#ifdef USE_GDAS
			model_info(i,j)%col_o3 = linearinterpolation( xx, (/ met_grid(ii,jj,wind_start+6), &
											met_grid2(ii,jj, wind_start+6) /), MYTIME*1.)
			o3_start = wind_start + 6
#else
			model_info(i,j)%col_o3 = 0.0
			o3_start = wind_start + 4
#endif

			
			do k=1, 6
				o3_temp(k) = linearinterpolation( xx, (/ met_grid(ii,jj,o3_start+k), &
											met_grid2(ii,jj, o3_start+k) /), MYTIME*1.)			
			end do
			
#ifdef USE_ECMWF
			model_info(i,j)%LSM = met_grid(ii,jj, num_gdas_vars-1)
#endif
			
						
			kstart = model_coarse / 2
			do k=kstart, model_coarse
			
				if (model_info(i,j)%Ps > 0. .and. p(k) > model_info(i,j)%Ps) then 
					if (model_info(i,j)%surface_level == 0) then
						
						model_info(i,j)%surface_level = k
						temp(k) = Ts2m
						if ( (model_info(i,j)%Ps - p(k-1)) < 5. .or. &
							 (p(k) - model_info(i,j)%Ps) < 5. ) then 
							p(k) = (p(k) + p(k-1)) / 2.
						else
							p(k) = model_info(i,j)%Ps
						endif
						mixr(k) = W2m
					
					else
					
						temp(k) = Ts2m
						mixr(k) = W2m
						p(k) = p_source(k-1)
					
					endif
				
				endif

			end do


! add the surface level into the coarse profile
		   if (model_info(i,j)%surface_level /= model_coarse) then 
				p(model_coarse) = p_source(model_coarse)
			else
				p(model_coarse) = model_info(i,j)%Ps
			endif
			temp(model_coarse) = Ts2m
			mixr(model_coarse) = W2m



! now we determine the lowest valid level of the new high-res profile
			kstart = model_levels / 2
			do k = kstart, model_levels
				if (pressures(k) >= model_info(i,j)%Ps) then
					model_info(i,j)%surface_level = k
					exit
				endif
			end do


	

			
! interpolate the profile to 101 levels, courtesy UW-Madison

			status = profile_to_101(p, temp, mixr, model_coarse, model_lat, pressures, temp_hires, mixr_hires)



			sfc_level = model_info(i,j)%surface_level
			a_factor = (model_info(i,j)%Ps - pressures(sfc_level-1)) / &
						( pressures(sfc_level) - pressures(sfc_level-1))
			temp_hires(sfc_level) = temp_hires(sfc_level-1) + a_factor * (temp_hires(sfc_level) - temp_hires(sfc_level-1))
			mixr_hires(sfc_level) = mixr_hires(sfc_level-1) + a_factor * (mixr_hires(sfc_level) - mixr_hires(sfc_level-1))
			
			pressures(sfc_level) = model_info(i,j)%Ps

			model_info(i,j)%temp_profile = temp_hires
			model_info(i,j)%mixr_profile = mixr_hires
							
! calculate the height profile, courtesy UW-Madison			

			status = height_profile(pressures, temp_hires, mixr_hires, &
								heights, model_levels, Pmsl*1.0d0)


			model_info(i,j)%height_profile = heights	
			model_info(i,j)%pressure_profile = pressures				

	
			model_info(i,j)%trop_level = find_trop (model_info(i,j)%temp_profile, &
											model_info(i,j)%pressure_profile, sfc_level)


			status = adjo3(pressures, o3_temp, model_info(i,j)%o3_profile, model_info(i,j)%o3_profile)

			
!			if (i==1 .and. j==91) then 
!				do k=1, model_levels
!					print*, k, model_info(i,j)%pressure_profile(k), model_info(i,j)%mixr_profile(k), model_info(i,j)%temp_profile(k)
!				end do
!			end if
	
			
		end do
	end do



 endif 
	


  
! get lat/long data from grids

	do grid_j = 1, data_ysize
		do grid_i = 1, data_xsize

! don't waste time processing and setting ancillary if there is no cloud. Why bother?
      if (lon(grid_i,grid_j) <= -999. .or. lat(grid_i,grid_j) <= -999.0 ) then 
		
		sfctmp(grid_i, grid_j) =  -999.
		
		cycle
		
      endif
	
	  call get_model_idx(lat(grid_i,grid_j), lon(grid_i,grid_j), i, j)

	  i = i-1
	  j = j-1

	  call get_bilinear_idx(lat(grid_i,grid_j), lon(grid_i,grid_j), i, j, idx_x, idx_y)
	
		Ts_forint(1,1) = model_info(idx_x(1), idx_y(1))%Ts
		Ts_forint(1,2) = model_info(idx_x(1), idx_y(2))%Ts
		Ts_forint(2,1) = model_info(idx_x(2), idx_y(1))%Ts
		Ts_forint(2,2) = model_info(idx_x(2), idx_y(2))%Ts


                if (idx_x(1)-1 < 180) then
                        lon_set(1) = (idx_x(1) - 1) * 1.0 !+ 0.5
                else
                        lon_set(1) = (idx_x(1) - 1) * 1.0 - 360. !+ 0.5
                endif

                if (idx_x(2)-1 < 180) then
                        lon_set(2) = (idx_x(2) - 1) * 1.0 !+ 0.5
                else
                        lon_set(2) = (idx_x(2) - 1) * 1.0 - 360. !+ 0.5
                endif


                if (lon(grid_i, grid_j) > 179 ) then !.5) then
                        lon_set(1) = 179 !.5
                        lon_set(2) = 180 !.5
                endif

                if (lon(grid_i, grid_j) < -179 ) then !.5) then
                        lon_set(1) = -180 !.5
                        lon_set(2) = -179 !.5
                endif


                lat_set = 90-(idx_y-1)*1.0 ! - 0.5
                if (lat(grid_i, grid_j) > 89 ) then !.5) then
                        lat_set(1) = 90.0
                        lat_set(2) = 89 !.5
                endif

                if (lat(grid_i, grid_j) < -89 ) then !.5) then
                        lat_set(1) = -89. !.5
                        lat_set(2) = -90.
                endif


      sfctmp(grid_i, grid_j) = bilinear_interpolation(  lon_set, &
					lat_set, lon(grid_i, grid_j), lat(grid_i, grid_j), Ts_forint, 1)


    enddo
  enddo
  
  
 end subroutine read_ancillary_grids



 end module ct_ancillary
