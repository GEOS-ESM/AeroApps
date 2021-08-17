module modis_albedo


	implicit none

	include "hdf.f90"
	include "dffunc.f90"

	real, parameter :: sza_th=0.025, vza_th=0.025, phi_th=3.0, ws_th = 1.0
	real :: prev_sza, prev_vza, prev_phi, prev_ws
	
	real :: mid_albedo(15, 3)

contains


! this subroutine returns the needed array bounds for a resized array
 subroutine find_bounds(array, min_val, max_val, min_bnd, max_bnd) 
 
	implicit none
	
	real, dimension(:), intent(in) :: array
	real, intent(in) :: min_val, max_val
	integer, intent(inout) :: min_bnd, max_bnd
	
	integer :: i, N
	
	N = size(array)
	
	do i=1, N
		if (array(i) >= min_val) then 
			min_bnd = i-1
			exit
		endif
	end do
 
	if (min_bnd < 1) min_bnd = 1
	if (min_bnd > N) min_bnd = N
 
	do i=1, N
		if (array(i) >= max_val) then 
			max_bnd = i
			exit
		endif
	end do

	if (max_bnd >= N) max_bnd = N
	if (max_bnd <= 1) max_bnd = 1

 end subroutine find_bounds


		subroutine interpolate_albedo(albedo_lib, &
						miu0_lib, miu_lib, phi_lib, &
						wind_speed, miu0, miu, phi, albedo_vec)
		
			use modis_numerical_module, only: BisectionSimple
		
			real, dimension(:,:,:,:,:), intent(in) :: albedo_lib
			real, intent(in) :: wind_speed, miu0, miu, phi
			real, dimension(:), intent(in) :: miu0_lib, miu_lib, phi_lib
			real, dimension(:), intent(inout) :: albedo_vec
		
			integer :: n1, n2, n3, ierror,i
			integer :: islow, ishi, ivlow, ivhi, iazlow, iazhi
	
			real :: dmu, dmu0, dazm, dmu_1, dmu_2, dmu0_1, dmu0_2, dazm_1, dazm_2	
			real :: volume, prod11, prod12, prod21, prod22						
			real :: vol1, vol2, vol3, vol4, vol5, vol6, vol7, vol8
											
		    real :: deltau															
			real, parameter :: lib_ws(3) = (/ 3.0, 7.0, 15.0 /)
			integer :: iuhi, iulow


			n1 = size(miu0_lib)
			n2 = size(miu_lib)
			n3 = size(phi_lib)

			CALL BisectionSimple(miu0,miu0_lib,n1,islow,ishi)
			CALL BisectionSimple(miu,miu_lib,n2,ivlow,ivhi)
			CALL BisectionSimple(phi,phi_lib,n3,iazlow,iazhi)




			dmu = miu0_lib(ishi) - miu0_lib(islow)
			dmu0  = miu_lib(ivhi)  - miu_lib(ivlow)
			dazm = phi_lib(iazhi) - phi_lib(iazlow)


	dmu_1 = miu0 - miu0_lib(islow) ; dmu_2 = miu0_lib(ishi) - miu0
	dmu0_1  = miu  - miu_lib(ivlow)  ; dmu0_2  = miu_lib(ivhi)  - miu
	dazm_1 = phi  -  phi_lib(iazlow); dazm_2 = phi_lib(iazhi) - phi

		
	volume = dmu0 * dmu * dazm
	IF(volume == 0.0)THEN
		print*, "zero_volume: invalid data"		!invalid data 
		RETURN
	ENDIF


	prod22 = dmu0_2 * dmu_2
	prod12 = dmu0_1 * dmu_2
	prod11 = dmu0_1 * dmu_1
	prod21 = dmu0_2 * dmu_1

	vol8 =  abs(prod22 * dazm_2)
	vol7 =  abs(prod22 * dazm_1)
	vol6 =  abs(prod12 * dazm_1)
	vol5 =  abs(prod12 * dazm_2)

	vol1 =  abs(prod11 * dazm_1)
	vol2 =  abs(prod21 * dazm_1)
	vol3 =  abs(prod21 * dazm_2)
	vol4 =  abs(prod11 * dazm_2)
	!volume = vol1 + vol2 + vol3 + vol4 + vol5 + vol6 + vol7 + vol8
	
	mid_albedo(:,:) = ( albedo_lib(islow,ivlow,iazlow,:,:) * vol8  + &
						albedo_lib(ishi,ivlow,iazlow,:,:)  * vol3  + &
						albedo_lib(ishi,ivlow,iazhi,:,:)   * vol2  + &
						albedo_lib(islow,ivlow,iazhi,:,:)  * vol7  + &
						albedo_lib(islow,ivhi,iazlow,:,:)  * vol5  + &
						albedo_lib(ishi,ivhi,iazlow,:,:)   * vol4  + &
						albedo_lib(ishi,ivhi,iazhi,:,:)    * vol1  + &
						albedo_lib(islow,ivhi,iazhi,:,:)   * vol6 )/volume






!			if (.not. wind_speed_only) then 
			
!				idx1 = minloc(abs(miu0_lib - miu0))
!				idx2 = minloc(abs(miu_lib - miu))
!				idx3 = minloc(abs(phi_lib - phi))

!				mid_albedo(:,:) = albedo_lib(idx1(1), idx2(1), idx3(1), :,:)


!			endif
		
		
		if (wind_speed - lib_ws(1) <= 0.2) then 
			albedo_vec = mid_albedo(:, 1)
		elseif (wind_speed - lib_ws(3) >= -0.2) then 
			albedo_vec = mid_albedo(:, 3)
		elseif ( abs(wind_speed - lib_ws(2)) <= 0.2) then 
			albedo_vec = mid_albedo(:,2)
		
		else
		
			if (wind_speed - lib_ws(2) < 0.) then 
			
				deltau = lib_ws(2) - lib_ws(1)
				iuhi = 2
				iulow = 1
				
				albedo_vec = ( mid_albedo(:,1)  * (lib_ws(iuhi)- wind_speed) + &
                    		mid_albedo(:,2)  * (wind_speed - lib_ws(iulow) ))/deltau
			else
			
				deltau = lib_ws(3) - lib_ws(2)
				iuhi = 3
				iulow = 2
				
				albedo_vec = ( mid_albedo(:,2)  * (lib_ws(iuhi)- wind_speed) + &
                    		mid_albedo(:,3)  * (wind_speed - lib_ws(iulow) ))/deltau
			
			endif 
			
		 endif
		 
		
		end subroutine interpolate_albedo
	
	

	subroutine get_db_albedo (albedo_name, lat, lon, albedo_array)


	include "hdf.f90"
	include "dffunc.f90"

	
		character(*), intent(in) :: albedo_name
		real, dimension(:,:), intent(in) :: lat, lon
		real, dimension(:,:), intent(inout) :: albedo_array
			
			
		integer :: file_id, var_id, err_code, start(2), stride(2), edge(2)
		integer ::  i , j

		character(len=20) :: sds_name


	! emissivity map variables
		real, parameter :: emissivity_res = 0.1
		integer, parameter :: num_emiss_cols = 3600
		integer, parameter :: num_emiss_rows = 1800
		real, dimension(:,:), allocatable :: emissivity_map
		integer :: emiss_x, emiss_y, NumberOfXPoints, NumberOfYPoints, NumMapCols, NumMapRows
		real :: emiss_west, emiss_north
		integer :: startx, starty, stopy, stopx
	
		real :: minLat, minLong, maxLat, maxLong
		integer :: wid, ht
	

		wid = size(lat, 1)
		ht = size(lat, 2)

		allocate(emissivity_map(num_emiss_cols, num_emiss_rows))
		start = 0
		stride = 1
		edge = (/num_emiss_cols, num_emiss_rows /)
		

		file_id = sfstart(trim(albedo_name), DFACC_READ)
		sds_name = "412_all"
		var_id = sfselect(file_id, sfn2index(file_id, trim(sds_name)))
		err_code = sfrdata(var_id, start, stride, edge, emissivity_map)
		err_code = sfendacc(var_id)
		err_code = sfend(file_id)


		do j=1, ht
			do i=1, wid
		
				if ( (Lat (i,j) >= -90.0000) .and. &
					(Lat (i,j) <=  90.0000) .and. &
					(Lon (i,j) >= -180.000) .and. &
					(Lon (i,j) <=  180.000) ) then 
		
! put the emissivity values where they belong 
					emiss_x = anint((180. + lon(i,j)) / emissivity_res) + 1
					emiss_y = anint((90. + lat(i,j)) / emissivity_res) + 1

!					print*, emiss_x, emiss_y, lon(i,j), lat(i,j)

					if (emiss_x < 1) emiss_x = 1
					if (emiss_x > num_emiss_cols) emiss_x = num_emiss_cols
					if (emiss_y < 1) emiss_y = 1
					if (emiss_y  > num_emiss_rows ) emiss_y = num_emiss_rows

					


					albedo_array(i,j) = emissivity_map(emiss_x, emiss_y)
					if (albedo_array(i,j) < 0.) albedo_array(i,j) = 2.5
		
				endif
		
			end do
		end do
		
	
	end subroutine get_db_albedo



  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------
  subroutine getAlbedoEco (Lat, Long, JulianDay, AlbedoWavelengthFN, EcosystemFN,           &
                           SnowAlbedoFN, snow_ice_fraction,   &
                           Albedo, land_mask                      )


    real, dimension(:,:),                    intent(in ) :: Lat, Long
    character (len = *),      &
              dimension(:),                      intent(in ) :: AlbedoWavelengthFN
    character (len = *),                         intent(in ) :: EcosystemFN
    character (len = *),                         intent(in ) :: SnowAlbedoFN
    integer,  intent(in ) :: JulianDay
	real, dimension(:,:), intent(in) :: snow_ice_fraction
	integer*1, dimension(:,:), intent(inout) :: land_mask

	integer*2, dimension(:,:,:), intent(inout) :: Albedo 


    !Local Variables

    !
    !Size variables:
    !
    !Number of Albedo Wavelengths to be processed.
    integer                       :: NumAlbWavesProcess
    !Number of Columns and Rows of the granule to be processed
    integer                       :: NumberOfCols, NumberOfRows
	!Number of Ecosystem classifications
    integer , parameter             :: NumEcosystems = 18

    !
    !Min/Max Lat/Lon:
    !
    real                                   :: minLat, maxLat, minLong, maxLong

    ! 
    !Albedo/Ecosystem map variables:
    !
    !Global map dims and maximum number of albedo wavelengths:
    integer                          :: NumMapCols, NumMapRows
    integer                         :: NumMapCols_alb, NumMapRows_alb
    integer  , parameter             ::  NumSnowTypes      = 2
                                                                 !dry = 1, wet = 2
    integer *2                        :: ii
                                                                 
    !Stores the resolution (in degrees) of the global map
    real                                      :: Map_Resolution
    real                                      :: Map_Resolution_alb
    !Stores the upper right corner global map coords, used to compute rest of map coords.
    real                                       :: westernMostLongitude,      &
                                                                 northernMostLatitude
    real                                      :: westernMostLongitude_alb,      &
                                                                 northernMostLatitude_alb
    !Indices for section of global map that is needed for supplying values for granule:
    integer                         :: startX, stopX,             &
                                                                 startY, stopY,             &
                                                                 NumberOfXPoints,           &
                                                                 NumberOfYPoints  
    integer                         :: startX_alb, stopX_alb,             &
                                                                 startY_alb, stopY_alb
	real                                     :: startXval, stopXval,       &
                                                                 startYval, stopYval 
    !Stores the real and integer values of the section of the global albedo map.
    integer*2,   &
             allocatable, dimension(:,:,:)                      :: AlbedoMap
    !Used to scale the albedo values read in as integers from the HDF files.
    real     , parameter                      :: ScaleFactor = 0.00100
    !Stores the snow albedo statistics:
    real    ,   &
             allocatable, dimension(:,:,:)                    :: SnowAlbedoStats,    &
                                                                 ProcessSnowAlbStats
    integer  , parameter             :: NumSnowAlbStatWaves = 10
    real  , dimension(NumSnowAlbStatWaves)    :: PossStatWavelengths
    logical                                                   :: FoundWave
    !Store the NISE snow type:
    INTEGER, allocatable, DIMENSION(:,:)    :: snowIceType
    !Store the section of the global ecosystem map needed for supplying values for granule:
    integer*1,   &
             allocatable, dimension(:,:)                      :: EcosystemMap 
    !Dummy index variables:
    integer                         :: xIndex, yIndex 
    integer                         :: xIndex_alb, yIndex_alb 
	!HDF variables:
    !File IDs
    integer                        :: EcoMapFID,      &
                                                                 AlbMapFID
    !SDS IDs
    integer                        :: EcoMapSDSID,    &
															AlbMapSDSID    
    !Stores the SDS name
    character(len = 64)                        :: SDSName
	
	integer :: Ecosystem 

    !Error handling:
    INTEGER   :: errorLevel

    !Wavelength descriptions used in snow statistics:
    
    !Variables used for sea ice albedos:
    integer  , parameter             :: NumSeaIceAlbWaves = 8
    real   , dimension(NumSeaIceAlbWaves)    :: PossSeaIceWavelengths
    real  , dimension(NumSeaIceAlbWaves)    :: meltSeaIceAlbedos,        &
                                                                 drySeaIceAlbedos
    real   , allocatable, dimension(:)       :: ProcessMeltSeaIceAlbedos, &
                                                                 ProcessDrySeaIceAlbedos,  &
                                                                 TransitionAlbedo,         &
                                                                 Slope, Intercept
    LOGICAL                                                   :: arcticMeltingSeason,      &
                                                                 WinterToMeltTransition,   &
                                                                 MeltToWinterTransition
    integer , parameter             :: TransitionDays = 10,      &
                                                                 NHMeltBeg = 152,          &
                                                                 NHMeltEnd = 244,          &
                                                                 SHMeltBeg = 334,          &
                                                                 SHMeltEnd = 61
    real   , parameter                       :: PercentageOfDry = 0.800d0
    
  
	character (len = 10),  dimension(7)       :: WavelengthText

    !Water albedo values:
    real   , parameter                        :: WaterAlbedo = 0.0250000
	real, dimension(:), allocatable :: albedo_real4
	
	integer :: file_id, var_id, err_code, start(2), stride(2), edge(2)
	
	
  !local variables:
  !counters:
  integer                          :: i,j,k,l,m,n

  !HDF variables:
  integer                          :: HDFstatus
  integer ,  &
             dimension(10)                         :: hdfStart, hdfStride, hdfEdge
  integer                        :: sds_id, sds_index
	
	
	
	
	WavelengthText(1) = "0.659"
	WavelengthText(2) = "0.858"
	WavelengthText(3) = "0.47"
	WavelengthText(4) = "0.555"
	WavelengthText(5) = "1.24"
	WavelengthText(6) = "1.64"
	WavelengthText(7) = "2.13"

	

    !*****************************************************************************
    ! Determine the size of the arrays:  
    !*****************************************************************************
    !Determine dimensions:
    NumberOfCols      = size( Albedo(:,:,:),   dim = 1 )
    NumberOfRows      = size( Albedo(:,:,:),   dim = 2 )
    NumAlbWavesProcess = size( Albedo(:,:,:),   dim = 3 )  
	
	allocate(albedo_real4(NumAlbWavesProcess))

    
    !*****************************************************************************
    ! Determine the min/max of lat/lon, which will be used to determine what portion
    !  of the albedo/eco maps need to be read in. 
    !*****************************************************************************
    !Determine min and max:
    if (maxval(Lat) == -999. .or. maxval(Long) == -999.) then
      !there is no gelocation data in the scan.
      Albedo(:,:,:) = -999.  
      return
    else
      minLat  = minval( Lat,  mask = ( (Lat  >= -90.0000 .and. Lat  <= 90.00000) ) )
      maxLat  = maxval( Lat,  mask = ( (Lat  >= -90.0000 .and. Lat  <= 90.00000) ) )
      minLong = minval( Long, mask = ( (Long >= -180.000 .and. Long <= 180.0000) ) )
      maxLong = maxval( Long, mask = ( (Long >= -180.000 .and. Long <= 180.0000) ) )
    endif
 
    
    do ii = 1,2
    
      if (ii == 1) then 
        NumMapCols = 43200
        NumMapRows = 21600
      else 
        NumMapCols = 21600
        NumMapRows = 10800
      end if 
      
      Map_Resolution = 360.000d0 / float(NumMapCols)
    
      !Compute the WesternMostLongitude and NothernMostLatitude of the maps:
      westernMostLongitude = Map_Resolution / 2.0 - 180.0
      northernMostLatitude = 90.0 - Map_Resolution / 2.0
      !print*,'mapres=',Map_Resolution,westernMostLongitude,northernMostLatitude
      !print*,'minLong=',minLong,maxLong,minLat,maxLat
   
      ! Because the albedo map res. is different than the ecosystem map, save the albedo map
      !  parameters for later determination surface albedo for each granule pixel  - GTA
      
      if (ii == 1) then 
        Map_Resolution_alb = Map_Resolution
        westernMostLongitude_alb = westernMostLongitude
        northernMostLatitude_alb = northernMostLatitude
        NumMapCols_alb = 43200
        NumMapRows_alb = 21600
      end if


      !Compute the start/stop x/y indices based on min/max lat/lon values:
      !use anint so that it rounds to nearest index.
      startx = anint( (minLong-westernMostLongitude) / (Map_Resolution) ) + 1
      stopx  = anint( (maxLong-westernMostLongitude) / (Map_Resolution) ) + 1
      starty = anint( (northernMostLatitude-maxLat ) / (Map_Resolution) ) + 1
      stopy  = anint( (northernMostLatitude-minLat ) / (Map_Resolution) ) + 1
    
      !Due to roundoff, geocoordinates of -90, 90, -180, 180 can result in 
      ! indices of one beyond the map bounds, so correct that here:
      if (startx == 0) startx = 1
      if (stopx  == NumMapCols+1 ) stopx = NumMapCols
      if (starty == 0) starty = 1
      if (stopy  == NumMapRows+1 ) stopy = NumMapRows
    
      !Compute the number of points:
      NumberOfXPoints = stopx - startx + 1
      NumberOfYPoints = stopy - starty + 1


      !Error Checking:
      if ( (startx < 1) .or. (stopx > NumMapCols) .or. &
           (starty < 1) .or. (stopy > NumMapRows) .or. &
           (startx > stopx) .or. (starty > stopy)      ) then
	    if (stopx > NumMapCols) stopx = NumMapCols
		if (stopy > NumMapRows) stopy = NumMapRows
		if (startx > stopx) startx = stopx-1
		if (starty > stopy) starty = stopy-1
		if (startx < 1) startx = 1
		if (starty < 1) starty = 1
		
      end if

      !*****************************************************************************
      ! Allocate the global albedo and ecosytem map storage arrays based upon these
      !  start and stop values.  Note that the arrays will be indexed using start/stop
     !  variables, instead of 1-NumberOfZPoints, for ease of selecting the nearest
      !  neighbor to the inputted lat/lon grid.
      !*****************************************************************************
      if (ii == 1) allocate( AlbedoMap    (startx:stopx, starty:stopy, 1:NumAlbWavesProcess) )
      if (ii == 2) allocate( EcosystemMap (startx:stopx, starty:stopy                      ) )
	
		
			    
      !*****************************************************************************
      ! Set up the HDF read variables, making note that HDF arrays start at 0
      !  instead of F90's 1:
      !*****************************************************************************
      !Define the starting point
      hdfStart (:) = 0
      hdfStart (1) = startX - 1
      hdfStart (2) = startY - 1
    
      !Define the stride = 1 so that it doesn't skip any data.
      hdfStride(:) = 1
    
      !Define the Number of Points to read:
      hdfEdge  (:) = 1
      hdfEdge  (1) = NumberOfXPoints
      hdfEdge  (2) = NumberOfYPoints

      
      if (ii == 1) then    ! read in the Surface Albedo Map data
        !*****************************************************************************
        ! Loop over each wavelength, except for 3.7, and read in the specified area
        !  of the albedo map for those wavelengths flagged for land processing.
        ! Also compute the 3.7 from the 2.1.
        !*****************************************************************************
        !Set the albedo map to fill value, so that wavelengths not needed for land 
        ! processing will be set to fill value (i.e. not read in)
        AlbedoMap(:,:,:) = -999.0000

        !Read in the portions of the albedo maps from their corresponding files:
        do i = 1, NumAlbWavesProcess 
          !Test to see if this wavelength should be read in:
            !Define the SDS Name:
            SDSName = "Albedo_Map_" // trim(WavelengthText(i))
            !Open the HDF file for reading:
            AlbMapFID = Sfstart(Trim(AlbedoWavelengthFN(i)), Dfacc_Read)
            !Obtain the data set ID's:
            AlbMapSDSID   = SFselect(AlbMapFID, SFn2index(AlbMapFID, trim(SDSName)))
            !Read in the data:
            HDFstatus = SFrdata(AlbMapSDSID, hdfStart, hdfStride, hdfEdge, AlbedoMap(:,:,i))
            !End Access to the datasets and files:
            HDFstatus = SFendacc( AlbMapSDSID )
            HDFstatus = SFend ( AlbMapFID   )   
        end do
	
      else
    
        !*****************************************************************************
        ! Read in the portion of the ecosystem map from the corresponding file:
        !*****************************************************************************
        !Define the SDS Name:
        SDSName = "IGBP_Land_Cover_Type"
    
        ! Open the HDF file for reading:
        EcoMapFID = Sfstart(Trim(EcosystemFN), Dfacc_Read)
        ! Obtain the data set ID's:
        EcoMapSDSID   = SFselect(EcoMapFID, SFn2index(EcoMapFID, trim(SDSName)))
        ! Read in the data:
        HDFstatus = SFrdata(EcoMapSDSID, hdfStart, hdfStride, hdfEdge, EcosystemMap)
        !End Access to the datasets and files:
        HDFstatus = SFendacc( EcoMapSDSID )
        HDFstatus = SFend ( EcoMapFID   )   
      end if
    end do

    !*****************************************************************************
    ! Allocate the Snow Albedo stats and types from 1 to NumberofRows and Cols:
    !*****************************************************************************
    allocate( SnowAlbedoStats    (1:NumSnowAlbStatWaves, 0:NumEcosystems-1, 1:NumSnowTypes) )
    allocate( ProcessSnowAlbStats(1:NumAlbWavesProcess,  0:NumEcosystems-1, 1:NumSnowTypes) )
	


    !*****************************************************************************
    ! Read in the snow albedo statistics from the corresponding file:
    !*****************************************************************************
    call ReadSnowAlbStats  ( SnowAlbedoFN, NumSnowTypes, NumSnowAlbStatWaves,  &
                             NumEcosystems, SnowAlbedoStats, errorLevel        )

    allocate( ProcessDrySeaIceAlbedos (1:NumAlbWavesProcess) )
    allocate( ProcessMeltSeaIceAlbedos(1:NumAlbWavesProcess) )
    allocate( TransitionAlbedo        (1:NumAlbWavesProcess) )
    allocate( Slope                   (1:NumAlbWavesProcess) )
    allocate( Intercept               (1:NumAlbWavesProcess) )

    !*****************************************************************************
    ! For each wavelength to be processed, store the albedo snow statistic in the
    !  cooresponding index in a temp array. Note that 8-10 are broad band, and 
    !  not included.
    !*****************************************************************************
    PossStatWavelengths(1)  = 0.470000
    PossStatWavelengths(2)  = 0.555000
    PossStatWavelengths(3)  = 0.659000
    PossStatWavelengths(4)  = 0.858000
    PossStatWavelengths(5)  = 1.240000
    PossStatWavelengths(6)  = 1.640000
    PossStatWavelengths(7)  = 2.130000
    

	! the snow albedo bands are ordered by wavelength rather than by channel, 
	! so we need to reshuffle a bit.
    ProcessSnowAlbStats= 0.

	ProcessSnowAlbStats(1,:,:) = SnowAlbedoStats(3,:,:)
	ProcessSnowAlbStats(2,:,:) = SnowAlbedoStats(4,:,:)
	ProcessSnowAlbStats(3,:,:) = SnowAlbedoStats(1,:,:)
	ProcessSnowAlbStats(4,:,:) = SnowAlbedoStats(2,:,:)
	ProcessSnowAlbStats(5,:,:) = SnowAlbedoStats(5,:,:)
	ProcessSnowAlbStats(6,:,:) = SnowAlbedoStats(6,:,:)
	ProcessSnowAlbStats(7,:,:) = SnowAlbedoStats(7,:,:)
	
                        
    !*****************************************************************************
    ! For each wavelength to be processed, store the albedo sea ica albedo in the
    !  cooresponding index in a temp array. 
    !*****************************************************************************
    !Load in the order of the statistical wavelengths from the file:
    PossSeaIceWavelengths(1) = 0.470000
    PossSeaIceWavelengths(2) = 0.555000
    PossSeaIceWavelengths(3) = 0.659000
    PossSeaIceWavelengths(4) = 0.858000
    PossSeaIceWavelengths(5) = 1.240000
    PossSeaIceWavelengths(6) = 1.640000
    PossSeaIceWavelengths(7) = 2.130000
    
    !Set the sea ice values to be equal to the dry perm snow/ice values or percentage of them:
    drySeaIceAlbedos(1:7)  = SnowAlbedoStats(1:7,15,1)
    meltSeaIceAlbedos(1:7) = SnowAlbedoStats(1:7,15,1) * PercentageOfDry
                         
	ProcessDrySeaIceAlbedos(1)  = drySeaIceAlbedos(3)
	ProcessMeltSeaIceAlbedos(1) = meltSeaIceAlbedos(3)
	ProcessDrySeaIceAlbedos(2)  = drySeaIceAlbedos(4)
	ProcessMeltSeaIceAlbedos(2) = meltSeaIceAlbedos(4)
	ProcessDrySeaIceAlbedos(3)  = drySeaIceAlbedos(1)
	ProcessMeltSeaIceAlbedos(3) = meltSeaIceAlbedos(1)
	ProcessDrySeaIceAlbedos(4)  = drySeaIceAlbedos(2)
	ProcessMeltSeaIceAlbedos(4) = meltSeaIceAlbedos(2)	
	ProcessDrySeaIceAlbedos(5:7)  = drySeaIceAlbedos(5:7)
	ProcessMeltSeaIceAlbedos(5:7) = meltSeaIceAlbedos(5:7)
    
    !*****************************************************************************
    ! Read in the NISE snow type from the corresponding file:
    !*****************************************************************************
                                      
  !print*,'tp latlon loop',NumberOfCols,NumberOfRows
    !*****************************************************************************
    ! Loop over each pixel in the lat/lon grid, compute the global coordinates and
    !  set the pixel's albedo and ecosystem values to the corresponding global
    !  albedo and ecosystem map values.  This is essentially a nearest neighbor 
    !  calculation.
    !*****************************************************************************
    do i = 1, NumberOfCols
       
      do j = 1, NumberOfRows
      
	  
        !Check to see if the geolocation is valid, also check we have a cloud. 
        if ( (Lat (i,j) >= -90.0000) .and. &
             (Lat (i,j) <=  90.0000) .and. &
             (Long(i,j) >= -180.000) .and. &
             (Long(i,j) <=  180.000) ) then 


   !print*,'rows=',j,i,Map_Resolution,Long(i,j),westernMostLongitude,northernMostLatitude,Lat (i,j)
         
          !Valid geolocation, so proceed.
          
          !Compute the ecosystem map x,y indices from this pixel's lat/long values.
          xIndex = anint( (Long(i,j)-westernMostLongitude) / (Map_Resolution) ) + 1
          yIndex = anint( (northernMostLatitude-Lat (i,j)) / (Map_Resolution) ) + 1

          !Due to roundoff, geocoordinates of -90, 90, -180, 180 can result in 
          ! indices of one beyond the map bounds, so correct that here:
          if (xIndex == 0) xIndex = 1
          if (xIndex  == NumMapCols+1 ) xIndex = NumMapCols
          if (yIndex == 0) yIndex = 1
          if (yIndex  == NumMapRows+1 ) yIndex = NumMapRows
          !print*,'xIndex, NumMapCols,yIndex, NumMapRows=',xIndex, NumMapCols,yIndex, NumMapRows

          !Compute the surface albedo map x,y indices from this pixel's lat/long values.
          xIndex_alb = anint( (Long(i,j)-westernMostLongitude_alb) / (Map_Resolution_alb) ) + 1
          yIndex_alb = anint( (northernMostLatitude_alb-Lat (i,j)) / (Map_Resolution_alb) ) + 1

          !Due to roundoff, geocoordinates of -90, 90, -180, 180 can result in 
          ! indices of one beyond the map bounds, so correct that here:
          if (xIndex_alb == 0) xIndex_alb = 1
          if (xIndex_alb  == NumMapCols_alb+1 ) xIndex_alb = NumMapCols_alb
          if (yIndex_alb == 0) yIndex_alb = 1
          if (yIndex_alb  == NumMapRows_alb+1 ) yIndex_alb = NumMapRows_alb
          !print*,'xIndex, NumMapCols,yIndex, NumMapRows=',xIndex, NumMapCols,yIndex, NumMapRows
          !Checks to make sure indices are within bounds of the map:

          !Store the ecosystem value:
          Ecosystem = EcosystemMap(xIndex,yIndex)

		  if (Ecosystem == 0) then 
			land_mask(i,j) = 0
		  else
			land_mask(i,j) = 1
		  endif


		  !print*,xIndex,yIndex,EcosystemMap(xIndex,yIndex)
		  albedo_real4(:) = AlbedoMap(xIndex_alb, yIndex_alb, :) * ScaleFactor
		  
		  
          !Set the albedo to either land values, water, or snow/ice, depending
          ! on the scenario:
          
          
          
          if (      (snow_ice_fraction(i,j) > 0.)   &
              .and. (snow_ice_fraction(i,j) <= 1.) &
              .and. (Ecosystem == 0)     ) then

              Ecosystem  = 15

              !This is sea ice, so set the either wet or dry scenario:
              ! Determine arctic melting season and set albedo values:
              if (lat(i,j) >= 0.0) then
                ! Arctic melting season: June 1st (day 152) to September 1st (day 244):
                SELECT CASE (JulianDay)
                  CASE ( (NHMeltBeg-TransitionDays):NHMeltBeg-1)
                    !transition period between winter and melt, so linear fit between them
                    ! to get albedo value:
                    Slope (:) = (ProcessMeltSeaIceAlbedos(:)-ProcessDrySeaIceAlbedos(:))   &
                              / real(TransitionDays) 
                    Intercept(:) = ProcessDrySeaIceAlbedos(:) - Slope(:) * real(NHMeltBeg-TransitionDays) 
                    TransitionAlbedo(:) = Slope(:) * real(JulianDay) + Intercept(:)    
                    albedo_real4(:) = TransitionAlbedo(:) * snow_ice_fraction(i,j) + WaterAlbedo * (1.- snow_ice_fraction(i,j))
                  CASE ( NHMeltBeg:NHMeltEnd )
                    !Melt Sea Ice Values:
                    albedo_real4(:) = ProcessMeltSeaIceAlbedos(:) * snow_ice_fraction(i,j) + WaterAlbedo * (1. - snow_ice_fraction(i,j))
                  CASE ( NHMeltEnd+1:(NHMeltEnd+TransitionDays) )
                    !transition period between melt and winter, so linear fit between them
                    ! to get albedo value:
                    Slope (:) = (ProcessDrySeaIceAlbedos(:)-ProcessMeltSeaIceAlbedos(:))   &
                              / real(TransitionDays) 
                    Intercept(:) = ProcessMeltSeaIceAlbedos(:) - Slope(:)  &
                                 * real(NHMeltEnd) 
                    TransitionAlbedo(:) = Slope(:) * real(JulianDay) + Intercept(:)    
                    albedo_real4(:) = TransitionAlbedo(:) * snow_ice_fraction(i,j) + (WaterAlbedo * (1. - snow_ice_fraction(i,j)))
                  CASE DEFAULT
                    !Dry Sea Ice values::
                    albedo_real4(:) = ProcessDrySeaIceAlbedos(:) * snow_ice_fraction(i,j) + (WaterAlbedo * (1. - snow_ice_fraction(i,j)))
                END SELECT
				
              else           
                ! Determine southern hemisphere melting season:
                ! Arctic melting season: December 1st (day 334) to March 1st (day 61):
                SELECT CASE (JulianDay)
                  CASE ( (SHMeltBeg-TransitionDays):SHMeltBeg-1)
                    !transition period between winter and melt, so linear fit between them
                    ! to get albedo value:
                    Slope (:) = (ProcessMeltSeaIceAlbedos(:)-ProcessDrySeaIceAlbedos(:))   &
                              / real(TransitionDays) 
                    Intercept(:) = ProcessDrySeaIceAlbedos(:) - Slope(:)  &
                                 * real(SHMeltBeg-TransitionDays) 
                    TransitionAlbedo(:) = Slope(:) * real(JulianDay) + Intercept(:)    
                    albedo_real4(:) = TransitionAlbedo(:) * snow_ice_fraction(i,j) + (WaterAlbedo * (1.-snow_ice_fraction(i,j)))
                    
                  CASE ( SHMeltBeg:366 )
                    !Melt Sea Ice Values:
                    albedo_real4(:) = ProcessMeltSeaIceAlbedos(:) * snow_ice_fraction(i,j) + (WaterAlbedo * (1.-snow_ice_fraction(i,j)))
                  CASE ( 1:SHMeltEnd )
                    !Melt Sea Ice Values:
                    albedo_real4(:) = ProcessMeltSeaIceAlbedos(:) * snow_ice_fraction(i,j)+ (WaterAlbedo * (1.-snow_ice_fraction(i,j)))
                  CASE (SHMeltEnd+1:(SHMeltEnd+TransitionDays) )
                    !transition period between melt and winter, so linear fit between them
                    ! to get albedo value:
                    Slope (:) = (ProcessDrySeaIceAlbedos(:)-ProcessMeltSeaIceAlbedos(:))   &
                              / real(TransitionDays) 
                    Intercept(:) = ProcessMeltSeaIceAlbedos(:) - Slope(:)  &
                                 * real(SHMeltEnd) 
                    TransitionAlbedo(:) = Slope(:) * real(JulianDay) + Intercept(:)    
                    albedo_real4(:) = TransitionAlbedo(:) * snow_ice_fraction(i,j) + (WaterAlbedo * (1.-snow_ice_fraction(i,j)))
	
                  CASE DEFAULT
                    !Dry Sea Ice values::
                    albedo_real4(:) = ProcessDrySeaIceAlbedos(:) * snow_ice_fraction(i,j) + (WaterAlbedo * (1.-snow_ice_fraction(i,j)))
                END SELECT


				
              end if

				Albedo(i,j,:) = nint(albedo_real4(:)*1000.)
			


          else if ( Ecosystem == 15 ) then
              !This is permanent snow/ice, so set to either map or stat value:
              if ((Ecosystem == 15) .and. &
                  (      AlbedoMap(xIndex_alb,yIndex_alb,1) > 0 &
                   .and. AlbedoMap(xIndex_alb,yIndex_alb,1) <= 1000) ) then
                !This pixel is perm snow class and has a valid map value,
                ! so set to the map's value:
                Albedo(i,j,:) = AlbedoMap(xIndex_alb,yIndex_alb,:)
              else
                !This pixel isn't flag as perm ice eco, so set it to the
                ! statistical value for perm ice:
                Ecosystem = 15
                Albedo(i,j,:) = nint (ProcessSnowAlbStats(:, Ecosystem, 1)*1000.)
              end if
			  
		  else if (snow_ice_fraction(i,j) > 0. .and. snow_ice_fraction(i,j) <= 1.0 .and. &
				Ecosystem /= 0 .and. Ecosystem /= 15.) then 
				
				Albedo(i,j,:) = nint( ProcessSnowAlbStats(:, Ecosystem, 1) * 1000.)
				Ecosystem = 15
			  
 
		else if (Ecosystem == 0) then
              !Water, so set all wavelengths to the water albedo value, 5%:
              Albedo(i,j,:) = nint(WaterAlbedo*1000.)
          else
              !This is a snow-free, ocean free, and not perm snow/ice pixel,
              ! so just store the value from the albedo map:
              Albedo(i,j,:) = AlbedoMap(xIndex_alb,yIndex_alb,:)

        end if

        else
          
          !Bad geolocation, so set albedos and emissivity to fill value
          Albedo(i,j,:) = -999
   
        end if
      
	if (Albedo(i,j,1) > 1000 .or. Albedo(i,j,1) < 0 ) land_mask(i,j) = 0

		
       end do
   
    end do
    





    !*****************************************************************************
    ! Deallocate arrays:
    !*****************************************************************************
    deallocate( AlbedoMap        )
    deallocate( EcosystemMap     ) 
    deallocate( ProcessSnowAlbStats )
    deallocate( SnowAlbedoStats  )
    deallocate( ProcessDrySeaIceAlbedos )
    deallocate( ProcessMeltSeaIceAlbedos ) 
    deallocate( TransitionAlbedo )
    deallocate( Slope            )
    deallocate( Intercept        )
	deallocate(albedo_real4)
	
		
  end subroutine getAlbedoEco
  
  
  ! ----------------------------------------------------------------------------------

  subroutine ReadSnowAlbStats  ( StatsFN, NumSnowTypes, NumAlbBnds, numEco, &
                                 AlbedoMean,  errorLevel                    )
                                 
    character (len = *),                 intent( in)  :: StatsFN
    real      ,  &
              dimension(:,0:,:),         intent(inout):: AlbedoMean
    integer   , intent( in)  :: NumAlbBnds,       &
                                                         numEco,           &
                                                         NumSnowTypes
    INTEGER,            INTENT(OUT) :: errorLevel
                                                         
    ! !Description:
    !   This routine will output only a single band of albedo, for the specified amount
    !    of data, to the specified position within the hdf file.
    !  
    ! !Input Parameters:
    !    StatsFN : The name, only, of the new HDF file.
    !    Albedo          : Contains the single band of albedo to be stored.
    !    NumSnowTypes    : 2, Dry or wet, via NISE classes.
    !    NumAlbBnds      : Number of Albedo wavelengths.
    !    numEco          : Number of Ecosystem Classes.
    !
    ! !Output Parameters:
    !    AlbedoMean      : Mean statistic
    !
    ! !Revision History:
    !   See Module revision history at the beginning of the file.
    !
    ! !Team-Unique Header:
    !   Cloud Retrieval Group, NASA Goddard Space Flight Center
    !
    ! !References and Credits:
    !   Written by
    !    Eric Moody
    !    Climate and Radiation Branch, Code 913
    !    NASA/GSFC
    !    Greenbelt MD 20771
    !    Eric.Moody@gsfc.nasa.gov
    !
    ! !Design Notes:
    !
    ! !END
    
    !local variables: 
    !HDF variables:
    integer                 :: status, i, j, k, HDFstatus
    integer ,  &
               dimension(5)  :: hdfStart, hdfStride, hdfEdge
    integer                 :: sds_id, newHDFID, sds_index
    character(len =  60)                :: SDSName
    real      ,  &
              dimension(1:2,1:NumAlbBnds,0:numEco,1:1,1:2)   :: DummyAlb
    !  status    : Used for error checking.
    !  hdfstart  : Follows HDF conventions. The array element from which to begin reading the
    !               HDF array. Zero based.
    !  hdfstride : Follows HDF conventions. The frequency with which data should be read in
    !               each dimension. Stride 1 means read every element.
    !  hdfedge   : Follows HDF conventions. The number of elements to read in each dimension. 
    !               If start(:) == 0 and stride(:) == 1, setting edge equal to the shape of 
    !               the underlying array reads all elements.
    !  sds_id,
    !   newHDFID : HDF SDS ID.
    !  sds_index : Index of the SDS in the HDF file.
    ! SDSName    : Stores the name of the SDS being procesed.


  !Set error level
  errorLevel = 0 

    !************************************************************************************
    ! Open the input file:
    !************************************************************************************
    !Open the albedo HDF file:
    newHDFID = SFstart( trim(StatsFN), DFACC_READ )


    !************************************************************************************
    ! Input the Albedo mean statistic:
    !************************************************************************************   
    !Set up the output hdf variables, in this case we are writing a Belt:
    ! hdfStart, note that we are storing the entire Longitude, Albedo Bands and Ecosystems
    ! so these start at 0, however we are storing a box of Latitude, so this starts at the
    ! starty counter.
    hdfStart ( : ) = 0
    hdfStart ( 5 ) = 0  !Latitude, 0 = NH, 1 = SH
    hdfStart ( 4 ) = 0  !Longitude, only 1, so 0
    hdfStart ( 3 ) = 0  !Eco
    hdfStart ( 2 ) = 0  !AlbBands
    hdfStart ( 1 ) = 0  !NumSnowTypes

    !Strides = 1, since we are not skiping points
    hdfStride( : ) = 1
    
    !Edge is the total number of points being read:
    ! It is the number of Albedo bands, the number of Ecosystems, the Number of Longitude
    ! boxes, and 1 lat box:
    hdfEdge  ( : ) = 1
    hdfEdge  ( 5 ) = 2
    hdfEdge  ( 4 ) = 1
    hdfEdge  ( 3 ) = numEco
    hdfEdge  ( 2 ) = NumAlbBnds
    hdfEdge  ( 1 ) = NumSnowTypes
    
    !Read the Mean data:
    SDSName = 'Snow_Albedo_Year_Mean'
    !determine the sds_index for the SDS:
    sds_index = SFn2index( newHDFID, trim(SDSName) )
    !get access to this SDS:
    sds_id = SFselect(newHDFID,sds_index)
        
    !Read the data:
    status = SFrdata(sds_id, hdfStart, hdfStride, hdfEdge, DummyAlb)

    !Store the dummy data in the final array:
    do i = 1, NumAlbBnds
      do j = 0, numEco-1
        do k = 1, NumSnowTypes
           AlbedoMean(i,j,k) = DummyAlb(k,i,j,1,1)
        end do
      end do
    end do

    
    !End Access to the SDS
    status = sfendacc(sds_id)
    !************************************************************************************
    ! Close the file:
    !************************************************************************************
    !Close the  HDF file:
    status = SFend(newHDFID)
 
  end subroutine ReadSnowAlbStats
  
  

end module modis_albedo
