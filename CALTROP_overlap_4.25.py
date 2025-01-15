""""
1/7/2025 - note ver _4.25.py is the version to use with the wrapper_CALTROP_4.1.py
Based on ver _4.1.py.

1/1/2025 - same as 4.1 but modified to deal with cases where there is a single IOWA file. Previos cases
had sevearl IOWA files in the directory and the code sorted it out. But with one, it gets confused.
So, this is code fixes this. Esentially I split the decision tree of matching orbits if 
there is only one Iowa file

Dec/30/2024 - ver 4.1 py. Running the rest of IOWA cases. 

Dec/30/2024 - ver4.0.py this code works in Discover with the same cases I used for AGU.
an important upgrade from Calculon is that all modules have been moved to 
an external code f_collocate_modules_1.py

12/10/204  Copy of ver 3.2 usef in CAlculon.
This copy is adapted to work in Discover
nonbackup/Python/AOCH. Essentially I had to adjust the path to Discover locations

12/3/2024 - ver 3.2 - Same as 3.1 with comestic changes such as removing print statement
and grouping information in different modules.


#-------------------------------
### #First, check the type of the object:
print(type(CAL)) 
# To see all attributes and methods of the object:
print(dir(CAL)) 
### If it's a dictionary-like object, you can try:
print(CAL.keys())  # if it's a dict-like object with .keys() method
# or
print(list(CAL))   # to see all keys/attributes
####For detailed information about the object's structure:
import inspect
print(inspect.__doc__)
## If it's a dictionary:
data = CAL['array1']
## If it's an object with attributes:
data = CAL.array1
###if you're still unsure, you can also try printing a small portion:
# For dictionary-like objects
print(list(CAL.items())[:5])  # First 5 items
# For general objects
print(vars(CAL))  # Shows object's attributes as dictionary
#### it it is a CLASS object use:
print(list(caldata.__dict__.keys())) # print the nonbuilt in attributes
print(dir(caldata)) # print all (ie builtin and non-built) in attributes

--------

Explore structures
TROP=TROPO()
vars(TROPO)
TROP.__dict__.items()
TROP.__dict__.keys() or print(CALI.keys())
help(TROP)
TROP.keys()

%who 

"""  
import os, sys,time,platform,re
import numpy as np
import f_collocate_modules_1 as f

from datetime import date, datetime, timedelta
from pyobs.calipso_l2 import CALIPSO_L2
from f_IOWA_loader import get_Iowa_heights
from netCDF4 import Dataset, date2num

###### -------- Main Code
if __name__ == "__main__":
   ### 
   # For CAL2.time[b[2] (2020-09-10 21:04:18), closest TROPO times are: 2020-09-10 20:34:03
   # 2021-08-28 two good orbits 
   #CAL_ 2021-08-28T14-33-47ZD.hdf : V. good dust sahara land
   #CAL_ 2021-08-28T21-07-48ZD.hdf : OK Smoke over Rockies

   #### set date  and hh,mm UTC time from CALIOP orbit to focus on
   # Comment out next line for operation with wrapper
   #    yyyy, mm, dd = '2020', '10', '01' ; cal_orb_time=[19,30] ;out_label='test'# Post AGU
   #   yyyy, mm, dd = '2020','08','21';cal_orb_time=[19,58];out_label='batch'
   sel_date_ymd     = datetime.strptime(yyyy + mm + dd, '%Y%m%d')
   cal_time_selected= datetime(year=int(yyyy),month=int(mm),day=int(dd), \
                               hour=cal_orb_time[0],minute=cal_orb_time[1])

   os_base_pth, tro_pthin, cal_full_path, pth_out = f.get_paths(yyyy, mm, dd)

   print('\n This code is operating from this folder ',os.getcwd())
   print('')
   ############
   ############ Get ALL CALIPSO orbits (and respective data) for this date
   print('Getting CALIPSO Level 2 data')
   ### Now sort filenames in ascending order in date and time
   cal_full_path.sort(key=lambda x: datetime.strptime(x.split('.')[-2].replace('ZD','Z'), '%Y-%m-%dT%H-%M-%SZ'))
   #### Get all orbits for selected day
   CAL2 = CALIPSO_L2(cal_full_path) #
#    print(list(CAL2.__dict__.keys()))
#    print(CAL2.filelist)
#    sys.exit()
   
   ########### Convert CAL2.time from list to numpy array. 
   CAL2.time = np.array(CAL2.time) 
   
   ### get some basic stuff 
   b=np.cumsum(CAL2.NOBS) # b contains the indexes in CAL2 where each CALIOP orbit starts
   b=np.insert(b,0,0)
   CAL_NORBS=len(b)-1
   for i in range(CAL_NORBS):
      print(f"   Orbit {i} starts at {CAL2.time[b[i]]} , ends at {CAL2.time[b[i+1]-1]}") 
   print('DONE with reading CALIOP orbits\n')

   ########### now get Iowa data 
   print('GET Tropomi data')
   tro_path2file = os_base_pth + tro_pthin +  yyyy + '/' + yyyy + '-' + mm + '-' + dd + '/'
   print("  Reading files from ",tro_path2file)
   var_names2read = ['lat', 'lon', 'height','aod']
   IWVAR, IWFILES= get_Iowa_heights(tro_path2file, var_names2read)
   TRO_NORBS=len(IWFILES)
   if TRO_NORBS == 0:sys.exit('Zero TROPO orbits read')    ## sanity check
   print("  Data loaded successfully.")
   print(f"  Number of orbits: {len(IWFILES)}")
   ### Now loop over the TROPO times (using filename since there is no Time SDS in the file)
   for i in range(TRO_NORBS):
      print(f"  {i} ,time hhmmss:{IWFILES[i][15:21]}  , array size lines x cols {np.shape(IWVAR['lat'][i])} ")
     # print(f"{var} shape: {np.shape(VAR[var])}")
   print('DONE Loading TROPO orbits\n')
   
   print("Now match each CALIOP orbit with closest TROPomi orbit in time") 
   ### get times from Tropomi files
   # Define a list to store the datetime objects
   TROP_time = []
   # Extract Date from each IOWA filename and convert to datetime element
   for filename in IWFILES:
      # Extract date and time components based on fixed positions in the string
      year = int('20' + filename[8:10])   # YY -> 20YY
      month = int(filename[10:12])        # MM
      day = int(filename[12:14])          # DD
      hour = int(filename[15:17])         # hh
      minute = int(filename[17:19])       # mm
      second = int(filename[19:21])       # ss
      # Create a datetime object and add to the list
      dt = datetime(year, month, day, hour, minute, second)
      TROP_time.append(dt)
   
   #### now match orbits
   #### Since ALL CALIOP orbits for this date have been read and stored in CAL2,
   #### then loop over each CAL orbit and compare the respective start time with
   #### with IOWA filenames (stored in TROP_time) and save closest (or two closest as in previous versions)
   
   base_date = CAL2.time[b[0]].date()  # Get the date from first TROPO time
   ## Next set CAL orbit time in datetime format
   cal_datetime = datetime.combine(base_date, datetime.min.time().replace(hour=cal_orb_time[0], minute=cal_orb_time[1]))
   ### initialize a couple of temporary vars
   min_time_diff = timedelta.max
   closest_index = 0
   
   # Now compute time difference between CALIOP orbits with each of the TROP orbits available
   # This next loop is necessarily if dealing with multiple CAL and IOWA orbits, but if using case
   # studies where usually it is one or two caliop orbits per day, this loop is not necessary
   # List to store the indexes of the two closest TROP orbits for each orbit in CAL2.time[b[i]]
   TROP_indx = []
   for i in range(CAL_NORBS):   
      CALstart=CAL2.time[b[i]]
      # Calculate absolute differences between `cal` and each element in `dt`
      differences = [(abs((t - CALstart).total_seconds()), idx) for idx, t in enumerate(TROP_time)]
      differences.sort()  # Sort by the absolute difference (first item in each tuple)
#       breakpoint()
      # Get the indexes of the two closest elements
#       IOW_indexes = [differences[0][1], differences[1][1]]
      IOW_indexes = [differences[0][1]]
      # IOW_indexes = differences.min()
      #len(TROP_index)= CAL_NORBS , TROP_INDEX contains the index in the list of IOWA
      #   files that has the closest time to each CALIOP orbit.
      #   if there is only 1 IOWA file, all values will be 0, meaning that only (and first) 
      #   IOWA file in the list matches the respective CAL orbit.
      TROP_indx.append(IOW_indexes)
      ### Find in the list of CAL orbits read, which one is closest to input cal_orb_time
      tmp0=CAL2.time[b[i]]
      time_diff = abs(tmp0 - cal_datetime)
      if time_diff < min_time_diff:
         min_time_diff = time_diff
         closest_index = i
   ###### 
   b_sel=closest_index # index in list of CAL orbits that corresponds to the user 
                       # selected CAL orbit as reported in cal_orb_time
   #TROP_indx=TROP_indx[b_sel]
   ##-----------------------------------------------------------
   # now test if the user provided a CALIOP time (ie if cal_orb_time exists). If yes,
   # it means that we need to use the closest orbit found in the previous statement. 
   # If not, it need to evaluate all time differences and decide whether to ignore because 
   # it is too much time between orbits. This is step is not considered next, as it is now 
   # (Jan/7/2025) it only work for user provided selected CALIOP orbits (ie cal_orb_time exists)
   # 
   if 'cal_orb_time' not in locals() :
      if TRO_NORBS != 1 :
         #### Now for selected caliop orbit (stored in cal_orb_time) ,find the closest TROP orbit 
         #### in time and also select the TROP orbit that will be used for comparison   
         # First convert cal_orb_time to a datetime object
         # Assuming that we will compare only hours and minutes for the same day (may miss cases over Dateline)
         base_date = TROP_time[TROP_indx[0][0]].date()  # Get the date from first TROPO time of TROP files avaialbe
         ## Next set CAL orbit time in datetime format
         cal_datetime = datetime.combine(base_date, datetime.min.time().replace(hour=cal_orb_time[0], minute=cal_orb_time[1]))
         ### initialize a couple of temporary vars
         min_time_diff = timedelta.max
         closest_index = 0
         #### Loop over TROP orbits
         for i, indexes in enumerate(TROP_indx):
             tropo_time = TROP_time[indexes[0]]
             time_diff = abs(tropo_time - cal_datetime)
             print('     ',min_time_diff,time_diff,closest_index,TROP_indx[closest_index],cal_datetime)
             if time_diff < min_time_diff:
                 min_time_diff = time_diff
                 closest_index = i
             print(f"  For CAL2.time[b[{i}]] ({CAL2.time[b[i]]}), closest TROPO time is: {TROP_time[indexes[0]]}")
         print("  Among these CALIOP and TROP Orbits , the closest match found in time")
         print(f"                         is for at CAL orbit with index b[{closest_index}]")
         
         print(f" Time difference: {min_time_diff}")
         print(f" The matching TROPO orbit time:ß {TROP_time[TROP_indx[closest_index][0]]}")   
      else :
         print('   *** NOTE there a single TROPOMI in IOWA directory ***')
         TROP_indx=TROP_indx[b_sel]
         print(f"  CALIOP orbit selected is # {b_sel} as ")
         print(f"  CAL2.time[b[{b_sel}]] ({CAL2.time[b[b_sel]]}), closest TROPO time is: {TROP_time[TROP_indx[0]]}")
   else : 
       print('User selected CALIOP orbit', cal_orb_time)
       print( "      The closest match CALIOP-TROPO found in time")
       print(f"              is for at CAL orbit with index b[{closest_index}]")
       print(f"      The time difference: {abs(TROP_time[TROP_indx[closest_index][0]]-cal_time_selected)}")
       print(f"      The matching TROPO orbit time: {TROP_time[TROP_indx[closest_index][0]]}")   

    
   ##########################----------------------------------------------------- 
   ########### GEO Collocation of tracks ---------------------- 
   ##########################----------------------------------------------------- 
   print("\nNow getting the TROPO pixels that contain the CALIOP orbit") 
   ## now loop over each CAL orbit. 
   ## Note an important consideration in the line of code in doing the matching CAL-TROP is that the 
   ## force-brute approache of directly comparing the CALIOP lat/longs with TROP lat/lon involves arrays too large
   ## for a simple distance computation. For comparing with CALIOP Level 2 array has ~4000x2 elements 
   ## whereas IOWA TROP is ~1000x900 elemetns. The latter can be larger is using TROP files that 
   ## span from pole to pole (maybe future versions of Iowa ?)
   ## 
   ## Also, this assumes using the existing Iowa dataset which only spans 0-50N, 120W to 30E , so there is not need to match the 
   ## orbits everywhere in the globe.
   ## So, the following code uses some shorcuts. The logic is as follows:
   ## 1) Take first and last scan line in TROP and find the first and last closest CALIOP pixel inside the TRO swath
   ##    This will return the first and last line of CALIOP data and the TROP column ranges that contain CALIOP pixels
   ##    So now the size of CALIOP and TROP arrays can be subsetted and the actuall distance CALIOP to TRO 
   ##    can be considered for each CAL pixel
   ## 2) Loop over each CALIOP pixel, find the closest TRO pixel and make a note of the lines and row 
   ## 
   ##Logic of next section
   ## Assume the temporal matching has been done then 
   ## 1) subset CALIOP orbit to lat range of TROPOMI data
   ## 2) then find CALIOP orbit inside the TROPOMI swath
   ## 3) quickly evaluate if there is TROPMI data along the CALIOP track 
   ##   (there should be some because I selected the case study have such data)  
   
   ### 
   if TRO_NORBS != 1 :
       ical = closest_index# 2  # For CAL2.time[b[2] (2020-09-10 21:04:18), closest TROPO times are: 2020-09-10 20:34:03
       ###### Get index for overlap TROP file  
       itro = TROP_indx[ical][0] # for selected CAL orbit, select the closest TROPMI orbit
   else :
       ical=b_sel
       ###### Get index for overlap TROP file  
       itro = TROP_indx[0][0]
   
   
   ### Set output file string with a label set by user
   output_file = 'IOWA-CALIOP_' + CAL2.filelist[ical][35:-4] + '_' + out_label + '.nc' 
#    if TRO_NORBS != 1 :
   print(f"\nWorking with CAL orbit ({CAL2.time[b[b_sel]]}), with closest TROPO time : {TROP_time[itro]}" )
#    else :
#       print(f"\nWorking with CAL orbit ({CAL2.time[b[b_sel]]}), with closest TROPO time : {TROP_time[itro[0]]}" )
   print("      CALIOP L2 filename :" , CAL2.filelist[ical])
   
#    sys.exit()

   ### To make life easier, rename TROPO orbit data
   tro_lat=IWVAR['lat'][itro]
   tro_lon=IWVAR['lon'][itro]
   tro_height=IWVAR['height'][itro]
   tro_aod=IWVAR['aod'][itro]
   
   ####
   #### 
   ###  Collocation Search
   ####
   ####
   
   print(' ')
   #### Indexes in CAL2 that corresponds to CALIOP orbit that is closest in time to the TROPO swath
   #### but still need to do the collocation of each pixel
   cal_orbit_indx=np.arange(b[ical],b[ical+1])
   cal_lat_data=CAL2.lat[cal_orbit_indx,1]
   cal_lon_data=CAL2.lon[cal_orbit_indx,1]
   tro_lines,tro_cols = tro_lat.shape
   
   ### Because IOWA data only spans 0 to 60N, subset the CALIOP data to these ranges.
   ### Take first scan line in TROP and find it inside the CALIOP track 
   ###(this may be done differently if doing full TRO orbits)
  
   ### Loop in every colum to find the colum in the first line contains the caliop track
   for icol in range(tro_cols):
      current_value=np.square(cal_lat_data[1000:2500]-tro_lat[0,icol])+np.square(cal_lon_data[1000:2500]-tro_lon[0,icol])
      min_index = np.argmin(current_value)
      min_val=current_value[min_index]
      if icol > 0 and min_val > prev_value:
         break  # Stop the loop if the value is lower
      prev_value = min_val  # Update for the next iteration
   icol=icol-1
#    sys.exit()
   #### now that the column was found, 
   print('First CALIOP pixels in first TROP line: TRO_Col, CAL_line ', min_index+1000,icol) # trop_colum=814,line=1910
   print('      CALIOP lat/lon ',[cal_lat_data[min_index+1000],cal_lon_data[min_index+1000]])
   print('      TROP   lat/lon ',[tro_lat[0,icol] , tro_lon[0,icol]])
   if abs(cal_lat_data[min_index+1000]-tro_lat[0,icol]) > 1:
      print('      ** NOTE Matching collocation is off by more than 1deg, likely CALIOP orbit is at edge of TROPOMI')
   index0=[min_index+1000, icol] #First CALIOP pixel that contain in first TROP line, column # in TROPO
 
   # sys.exit()
   #### do the same for the last TROP scam line
   for icol in range(tro_cols):
      current_value=np.square(cal_lat_data[2500:]-tro_lat[-1,icol])+np.square(cal_lon_data[2500:]-tro_lon[-1,icol])
      min_index = np.argmin(current_value)
      min_val=current_value[min_index]
      if icol > 0 and min_val > prev_value:
         break  # Stop the loop if the value is lower
      prev_value = min_val  # Update for the next iteration
   #### now that the column was found, 
   print('Last CALIOP pixel that contain in Last TROP line, column # in TROPO', min_index+2500,icol) # 
   print('      CALIOP ',[cal_lat_data[min_index+2500],cal_lon_data[min_index+2500]])
   print('      TROP   ',[tro_lat[-1,icol] , tro_lon[-1,icol]])
   index1=[min_index+2500, icol]
   print('')

   ####  Now gather information for next section
   cal_range=[index0[0],index1[0]] # indexes in CALIOP that spans TROP data
   tro_range=[[0,index0[1]],[tro_lines,index1[1]]] #indexes of rows&cols in TROP that includes CALIOP data
   
   ##### OK Now I have the first and last CAL pixel that are insider the TROP swath (stored in cal_range)
   ##### and I have the range of columns in TROP that include all CALIOP pixels. So next:
   ##### 1) subset TRO and CAL data 2) Loop over each CAL pixel and find it in the TROP subset 3)
   ##### save respective line and col number 
   tro_mincol=min([index0[1],index1[1]])
   tro_maxcol=max([index0[1],index1[1]])
   tro_cols =[tro_mincol,tro_maxcol]
   tro_scan=[0,tro_lines]
   #### Subset. if using the full tropo orbit , the first dimension should be constrained
   tro_sub_lat=tro_lat[:,tro_mincol:tro_maxcol+1]
   tro_sub_lon=tro_lon[:,tro_mincol:tro_maxcol+1]
   cal_Nlines=cal_range[1]-cal_range[0]+1
   # cal_ind_fine=np.arange(0, cal_Nlines )
   cal_ind_fine=np.arange(cal_range[0],cal_range[1]+1)
   
   ##  Loop over each CALIOP point , find shortest distance
   ## and save each tropomi pixel to which this CALIOP pixel is assigned
   ##
   start_time = time.perf_counter()
   ## Initialize a few things before loop. This array will save TROPOMI granule
   ## scan line that contains the cal_ind_fine[ical] profile
   saveit_cal=-999*np.ones((cal_Nlines, 2), dtype=int)
   tro_shift_lines=20
   ### optimize the loop by making the distance computation with a smaller TROPO
   ### array and use the fact that the index from the previous iteration
   ### can be used in current iteration to subset the TROPO lat/lon grid
   tro_start_line=tro_range[0][0] # for Iowa this should be 0, not for full orbit
   tro_Nlines_sub=tro_lines
   print("Now that first and last have been found,") 
   print("  loop over each CALIOP pixel and match it with each TROPO")
   for icount,ical2 in enumerate(cal_ind_fine):
      cal0=[cal_lat_data[ical2],cal_lon_data[ical2]]
      if tro_start_line < tro_Nlines_sub+1:
         #### Subset TROPOMI array to a smaller array
         #### Check if first , last and too close to the end by tro_shift_lines
         if tro_start_line==0:
            trlat_small=tro_sub_lat;
            trlon_small=tro_sub_lon
         elif tro_start_line < tro_Nlines_sub - tro_shift_lines:
            trlat_small=tro_sub_lat[tro_start_line:tro_start_line+tro_shift_lines,:]
            trlon_small=tro_sub_lon[tro_start_line:tro_start_line+tro_shift_lines,:]
         else:
            trlat_small=tro_sub_lat[tro_start_line:-1,:]
            trlon_small=tro_sub_lon[tro_start_line:-1,:]
      #### Compute actual distance
         dist=f.haversine(trlat_small,trlon_small,cal0[0],cal0[1])
         tmp=np.unravel_index(np.argmin(dist), dist.shape)
         if icount==0:
           saveit_cal[icount,:]=[tmp[0]+tro_scan[0],tmp[1]+tro_cols[0]]
         else:
            saveit_cal[icount,:]=[tmp[0]+saveit_cal[icount-1,0],tmp[1]+tro_cols[0]]
            tro_start_line=saveit_cal[icount,0]-tro_scan[0]
      else:
            print("Reached last line of TROPOMI daylight data. Ending matching process ",ical2)
      if (icount % 200) == 0 : print('   N = %8i of %8i'  %(icount,cal_Nlines))


    # End timer
   end_time = time.perf_counter()
    # Calculate elapsed time
   elapsed_time = end_time - start_time
   print("Secs. Elapsed time for mapping each Caliop pixe with TROPO: %10.4f \n" %elapsed_time)

   
   ##########################----------------------------------------------------- 
   ########### DONE with Collocations ---------------------- 
   ##########################----------------------------------------------------- 
   ## Now that we have the ranges spanned by both datasets, 
   ## we subset the CALIOP orbit that is exactly inside the TROPOMI
   
   #### CALIOP, use a dictionary to store the data 
   cal_sel = {} # initialize dictionary
   cal_sel['lat']=cal_lat_data[cal_range[0]:cal_range[1]+1] # 1-D array length = len(cal_range[0]:cal_range[1]+1)
   cal_sel['lon']=cal_lon_data[cal_range[0]:cal_range[1]+1] # 1-D array length = len(cal_range[0]:cal_range[1]+1)
   cal_sel['zlidarkm']=CAL2.zlidar[0:399] #1-D array length = 400


   #### Now within the selected CAL orbit, subset those CAL pixels that are collocated with
   #### TROPO data
   cal_collocated_indx=cal_orbit_indx[cal_range[0]:cal_range[1]+1]

   ## extract from the original dataset
   cal_sel['ext']  = CAL2.ext[cal_collocated_indx,:]  #2-D array shape = len(cal_collocated_indx) x 400
   cal_sel['ext2'] = CAL2.ext2[cal_collocated_indx,:] #2-D array shape = len(cal_collocated_indx) x 400
   cal_sel['tback']= CAL2.tback[cal_collocated_indx,:]#2-D array shape = len(cal_collocated_indx) x 400
   cal_sel['time'] = CAL2.time[cal_collocated_indx]   #1-D array shape = len(cal_collocated_indx) 
   # sys.exit()
   #### TROPOMI
   ### all next 2-D arrays are shape = len(cal_range[0]:cal_range[1]+1) x 2
   tro_sel={}
   tro_sel['lat'] = tro_lat[saveit_cal[:,0],saveit_cal[:,1]] 
   tro_sel['lon'] = tro_lon[saveit_cal[:,0],saveit_cal[:,1]]
   tro_sel['height'] = tro_height[saveit_cal[:,0],saveit_cal[:,1]]
   tro_sel['aod']    = tro_aod[saveit_cal[:,0],saveit_cal[:,1]]
   
   
   ##### Miscellanous information to save
   calipsol2_filename=CAL2.filelist[ical] # string with filename for CALIPSO
   iowa_trop_filename=IWFILES[itro]      #string with filename for IOWA TropOMI 
   calipso_indexes= [cal_range[0],cal_range[1]+1]  # start and end pixel number in selected CALIOP orbit that spans IOWA data 
   matching_indexes = saveit_cal # Array N lines x 2 columns where N = length(cal_sel('lat')) and
                                 # saveit_cal[:,0] = IOWA line number that contains respective CALIOP pixel
                                 # saveit_cal[:,0] = IOWA col number that contains respective CALIOP pixel    
   
   ####### Now save data
   # Create a new NetCDF file
   f.save_matched_data(pth_out, output_file, cal_sel, tro_sel,\
                     matching_indexes, calipso_indexes, 
                     calipsol2_filename, iowa_trop_filename)
   print(f"Data saved to {pth_out+output_file}")
   print(f"\n CAL_IOWA_File ='{output_file}' \n ") # cut and past as input in collocated_MERRA-CLIOP-IOWA_3.1.py

   print("\n DONE running script ", os.path.basename(__file__)) 
   ##### do some memory clean up 
   
   #### display arrays currently in memory
   # for name, obj in globals().items():
    # if not name.startswith('_') and isinstance(obj, np.ndarray):
        # print(f"{name}: {type(obj)}, Shape: {obj.shape}, Size: {obj.nbytes} bytes")
   # #### display lists currently in memory
   # for name, obj in globals().items():
       # if not name.startswith('_') and isinstance(obj, (np.ndarray, list)):
           # print(f"{name}: {type(obj)}, Length: {len(obj)}")
           
   ###### Erase everythign except noted variables
   # import numpy as np  # Example module import
   ##### Define some example arrays and variables
   # array1 = np.array([1, 2, 3])
   # array2 = np.array([4, 5, 6])
   # other_var = 42
   ##### List of variable names you want to keep
   # keep_vars = ['array1', 'array2'] #so variable other_var will be erased.
   ###### Delete everything else in globals() except loaded modules and specified arrays
   # for name in list(globals().keys()):
       # if name not in keep_vars and not name.startswith('_') and not isinstance(globals()[name], type(np)):
           # del globals()[name]

   # # Verify which variables remain
   # print("Remaining variables:", globals().keys())
   