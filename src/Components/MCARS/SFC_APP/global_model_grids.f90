module global_model_grids
  
   implicit none
  
  real ::  ice_grid( 0:719, 0:359)
  real :: ozn_grid(0:359, 0:180)


  logical  :: grids_are_read
  
  type ancillary_type

    real, dimension(:), allocatable :: mixr_profile
	real, dimension(:), allocatable :: temp_profile 
	real, dimension(:), allocatable :: height_profile
	real, dimension(:), allocatable :: pressure_profile
	real, dimension(:), allocatable :: o3_profile
  
    
	real :: wind_speed
	real :: Ts
	real :: Ps
	real :: col_o3
	
	integer*1 :: surface_level
    integer*1 :: trop_level
	integer*1 :: LSM
	
  end type ancillary_type
  
contains

 subroutine get_model_idx(lat, lon, i, j) 


	real, intent(in) :: lat, lon
	integer, intent(inout ) :: i,j

	 real :: x,y, x0, dx, y0, dy
	
      x = min( max( lon,  -179.99 ), 179.99 )
      if( x > -999. .and. x < 0.0 ) x = lon+ 360.0
      x0 = 0.0
      dx = 1.0
      i = int( ( x - x0 + 0.5*dx ) / dx ) 

      if( i .eq. 360 ) i = 0

      y = min( max( lat, -89.99 ), 89.99 )
      y0 = 90.0
      dy = -1.0
      j = int( ( y - y0 + 0.5*dy ) / dy )

	  i = i+1
	  j = j+1


end subroutine get_model_idx
  
  
  

end module global_model_grids
