# -*- coding: utf-8 -*-
"""
Plot attenuated backscatter or Bext curtains from CALIPOS Level 2 files.

-----------
Jan/13/2025 - original from Mac, created there and then moved here.

Jan/10/2024 - Code updated to work in Discover and generate a list of curtains for
an input list of cases as defined by the user.

Aug/23/2024 - Based on Leve1.5 code, this plots Level 2 data


Aug/19/2024: using as template plotcurtain_log_v2.py

it works. Just need to adjust the paths and file names.

    # f_get_track(path2file,filename)
    #caliop_plot2d(param, lat, lon, z, window_ave=0, colorange=(1e-6, 1e-3))
    # Good examples here
    #https://www.analyticsvidhya.com/blog/2020/09/colormaps-matplotlib/
    #https://github.com/NASA-DEVELOP/VOCAL/blob/master/calipso/dat/calipso-backscatter.cmap
    #http://meteothink.org/examples/meteoinfolab/satellite/calipso.html
    #https://pyhogs.github.io/colormap-examples.html good for understanding differences pcolor, pcountour...

@author: sgasso
"""

import os,sys,platform
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import ndimage
from pyhdf import HDF
from pyhdf.SD import SD, SDC
import pyhdf.VS

#### for plotting satellite track
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def load_hdf_arrays(file_path, array_names):
    """
    Loads specified arrays from an HDF file into memory using pyhdf.
    Args:
        file_path (str): Path to the HDF file.
        array_names (list of str): Names of variables to read.
    Returns:
        dict: A dictionary containing the loaded arrays.
    """
    try:
        # Attempt to open the HDF file
        try:
            hdf = SD(file_path, SDC.READ)
        except Exception as e:
            print(f"Error: Unable to open the HDF file '{file_path}'.")
            print(f"Details: {str(e)}")
            sys.exit(1)  # Terminate the program with an error code

        # Initialize an empty dictionary to store the arrays
        loaded_arrays = {}

        # Get the list of datasets in the file
        datasets = hdf.datasets()

        # Load each specified array
        for array_name in array_names:
            if array_name in datasets:
                sds = hdf.select(array_name)
                loaded_arrays[array_name] = sds.get()  # Load data
                print(f'Read {array_name}')
            elif array_name == 'Lidar_Data_Altitudes':
                # Altitude is not a dataset, it is stored in a VD structure
                hdf0 = HDF.HDF(file_path)
                vs = hdf0.vstart()
                xid = vs.find('metadata')
                altid = vs.attach(xid)
                altid.setfields('Lidar_Data_Altitudes')
                nrecs, _, _, _, _ = altid.inquire()
                altitude = altid.read(nRec=nrecs)
                altid.detach()
                alti = np.array(altitude[0][0])
                loaded_arrays[array_name] = np.flip(alti)
                print(f'Read {array_name}')
            else:
                print(f"Array '{array_name}' not found in the HDF file.")

        # Close the HDF file
        hdf.end()

        return loaded_arrays
    except Exception as e:
        print(f"Error loading HDF arrays: {e}")
        return None

#--------------------------------------------------------------------------
def f_caliop_plot2d(param, lat, lon, z,title_strings, colorange=(1e-6, 1e-3),\
                    scale='lin', window_ave=0,esY=10):
      """
      This function plots a 2D image from CALIOP data.
      Inspired by http://meteothink.org/examples/meteoinfolab/satellite/calipso.html
      and adapted/modified accordingly for this application by @sgassoumd

      Args:
          param  : n x m numpy array containing the data to be plotted. NOTE: param[0,:] = uppermost data point
          lat,lon: n element numpy array containing the latitude coordinates.
          z      : m element numpy array containing the altitude coordinates. z[0]=Bottom of profile
          colorange: tuple of two floats representing the color range for the plot (default: (1e-6, 1e-3))
          scale: 'lin' or 'log' use linear for data spanning within a x10 factor, use log for large magnitude range
          window_ave: integer specifying the window size for averaging (default: 0 - no averaging)


      Returns:
          himg_pos: numpy array containing the position of the image object in the figure.
      """
      title_str=title_strings[0]
      title_colbar=title_strings[1]



     # select spacing between points to display , useful for debuggin because speed
      esX = 1  # X-axis spacing for plotting , this is the height dimension
      # esY = 10  # Y-axis spacing for plotting , this is the lattitude dimension

      # Determine last index based on data orientation (assuming older format)
      # lastZ = len(param[:, 0])

      # Extract data with specified spacing
      var0 = param[::esX,::esY]
      # var= param[::esX,::-esY]
      varlat = lat[::esY]
      varlon = lon[::esY]
      # varz = np.flipud(z[::esX])
      varz = z[::-esX] # the minus sign reverses the other , so it matches the profile



      # Apply filter to reduce noise
      if window_ave == 0:
        var0[var0 < 1e-12] = np.nan
        var=var0
      else:
        # Apply averaging filter along the second dimension (X-axis)
        var = ndimage.uniform_filter(var0, size=window_ave, mode='constant')
        var[var <= 0] = np.nan

      plt_var=np.log10(var);vmin0=np.log10(colorange[0]);vmax0=np.log10(colorange[1]);
      space_col=30 # fixed number of colors to resemble official Quicklool

      #### Setup colors
      # Get the colormap
      base_cmap = plt.get_cmap('jet')
      # base_cmap = plt.get_cmap('RdYlBu_r')
      # Create color levels
      levels = np.linspace(vmin0, vmax0, space_col + 1)  # 30 spaces between 31 levels
      # Create a list of the first 20 colors from the base colormap
      calipso_cmap = [base_cmap(i) for i in np.linspace(0, 1, space_col-10)]
      # Add 10 shades of gray, from dark gray to white
      gray_colors = [colors.to_rgba(str(min(i/10, 1))) for i in range(1, 12)]
      calipso_cmap.extend(gray_colors) ## add the new color to exisint map
      # Create a ListedColormap
      cmap = colors.ListedColormap(calipso_cmap)
      # Create a BoundaryNorm
      norm = colors.BoundaryNorm(levels, len(calipso_cmap))


      ### Now plot
      # Create the figure
      # plt.figure(figsize=(13, 4))
      # Create the plot
      fig, ax = plt.subplots(figsize=(16, 6))

      plt.pcolor(plt_var,norm=norm,cmap=cmap)

      # Add colorbar with extended range
      cbar = plt.colorbar(extend='max')


      if scale=='log':
         # Place ticks in colorbars at specifcid locations
         # a=np.arange(1e-2,1.1e-1,0.01)
         # b=np.arange(1e-3,1e-2,0.001)
         # c=np.arange(1e-4,1e-3,0.0001)
         # bar_ticks=np.concatenate((c,b,a))
         bar_ticks= np.array(\
                   [1e-6, 2e-6, 3e-6, 5e-6,  7e-6,
                    1e-5, 2e-5, 3e-5, 5e-5,  7e-5,
                    1e-4, 2e-4, 3e-4, 5e-4,  7e-4,
                    1e-3, 2e-3, 3e-3, 5e-3,  7e-3,1e-2 ])
         # bar_ticks= np.array(\
         #              [1e-6, 2e-6, 3e-6, 5e-6,  7e-6,
         #               1e-5, 2e-5, 3e-5, 5e-5,  7e-5,
         #               1e-4, 2e-4, 3e-4, 5e-4,  7e-4,])

         log_bar = np.log10(bar_ticks)
         # Ensure first and last ticks are displayed
         # log_bar_with_ends = np.concatenate(([log_bar[0]], log_bar[::2], [log_bar[-1]]))
         log_bar_with_ends = np.concatenate(([log_bar[0]], log_bar, [log_bar[-1]]))
         cbar.ax.minorticks_off() # remove minor ticks
         cbar.set_ticks(log_bar_with_ends)
         tick_labels = [f"{10**label:.1e}" for label in log_bar_with_ends]
         cbar.set_ticklabels(tick_labels)
         # Adjust colorbar label positioning
         cbar.ax.tick_params(labelsize=8)
         cbar.set_label(title_colbar, rotation=270, labelpad=15)

      #### Now add a lat/Longs to x-axis
      # Ensure x-axis ticks include first and last data points
      xt = np.arange(len(varlat))
      if xt[0] != 0:               xt = np.insert(xt, 0, 0)
      if xt[-1] != len(varlat) - 1:xt = np.append(xt, len(varlat) - 1)
      #Same for y-axis
      yt = np.arange(len(varz))
      if yt[0] != 0:             yt = np.insert(yt, 0, 0)
      if yt[-1] != len(varz) - 1:yt = np.append(yt, len(varz) - 1)
      nXTicks=20;nYTicks=20
      # pre-select ~ 10-15 ticks in the horizontal make last and first are included.
      ixt=np.round(np.linspace(0, len(xt) - 1, nXTicks)).astype(int)
      iyt=np.round(np.linspace(0, len(yt) - 1, nYTicks)).astype(int)
      ## X axis ticks
      xvals = []
      xstrs = []
      # nx = atype.shape[1]
      for i in ixt:
            xvals.append(i)
            if i == 0:
                xstrs.append('Lat: %.2f\nLon: %.2f' % (varlat[i],varlon[i]))
            else:
                xstrs.append('%.2f\n%.2f' % (varlat[i],varlon[i]))
      ## Y -ticks
      # Set x and y axis labels and ticks
      # plt.xticks(xt[ixt], [f"{a:.2f}, {b:.2f}" for a, b in zip(varlat[xt[ixt]], varlon[xt[ixt]])])
      plt.xticks(xvals, xstrs,fontsize=8)
      plt.yticks(iyt[::-1], [f"{c:.1f}" for c in varz[iyt]],fontsize=8)


      # Miscellanous
      plt.ylabel('Height (km)')
      plt.xlabel('(Lat, Long)')
      if window_ave == 0 :
         plt.title(title_str )
      else :
         plt.title(title_str + ' , Moving average = ' + str(window_ave) )

      # Get the position of the image object in the figure
      # himg_pos = plt.gca().get_images()[0].get_extent()

      return

def plot_track(latitude, longitude, title):
    """
    Plots the satellite track on a global map with coastline contours.

    Parameters:
    latitude (list of float): List of latitude values.
    longitude (list of float): List of longitude values.
    title (str): Title of the plot.
    """
    # Create a plot with coastlines
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()

    # Add the satellite track to the plot
    ax.plot(longitude, latitude, color='red', linewidth=2, marker='o', transform=ccrs.Geodetic())

    # Add gridlines and labels for latitudes and longitudes
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

    # Add the title
    ax.set_title(title)

    plt.show()


# Example usage
if __name__ == "__main__":
   current_os = sys.platform
    # Get the operating system name and computer name
   current_os    = platform.system()
   computer_name = platform.node()
   print('Currently working from ', computer_name)
   # plot_track='yes' # If yes, add the CALIOP track
   print_fig='yes'

   if current_os == 'win32':  # PC office
       pthin = "D:\\Satellite\\Tropomi\\Level2\\226\\"
   elif current_os == 'darwin' or current_os == 'Darwin':  # Laptop
       pth_fig_out  = '/Users/sgasso/Pyfiles/Py3/Figures/'
       pth_cal_root = "/Volumes/ExtData1/SatData/Calipso/Level2/"
   elif current_os == 'linux' and "calculon" in computer_name:  # Calculon
       os_base_pth = '/nobackup/CALIPSO/Level1.5/Santiago'
       cal_pthin = '/Level2/'
       tro_pthin = '/IOWA/'
       pth_out = '/nobackup/CALIPSO/Level1.5/Santiago/Output/Collocations/'
   elif current_os == 'Linux' and "discover" in computer_name:  # Discover
       pth_cal_root = '/css/calipso-l2-aer/data/LID_L2_05kmAPro-Standard-V4-51/'
       pth_fig_out  = '/home/sgasso/Python/Figures/'
       # cal_folder_path = cal_pthin + yyyy + '/' + mm + '/'

       # # Find matching files
       # cal_full_path = []
       # for filename in os.listdir(cal_folder_path):
       #     if re.search(f"{yyyy}-{mm}-{dd}..........D\.hdf$", filename):
       #         cal_full_path.append(os.path.join(cal_folder_path, filename))
       # pth_out = os_base_pth + '/SliceTROCAL/'

   else:
       print('Current operating system not recognized.')
       print('Cannot set path to Level1 files. Terminate Run')
       sys.exit()

   filenames = [
      "AlongTrack_MERRA-IOWA-CALIOP_2020-08-21T19-58-50ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-08-26T13-16-35ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-07T20-49-39ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-12T20-41-20ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-14T18-39-46ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-14T20-18-16ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-15T17-38-59ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-15T19-17-29ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-16T16-38-17ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-16T18-16-47ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-16T19-55-17ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-17T20-33-00ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-18T19-32-13ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-20T19-09-10ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-09-27T15-20-46ZD_batch.nc",
      "AlongTrack_MERRA-IOWA-CALIOP_2020-10-01T19-30-09ZD_batch.nc"
   ]


#    filenames = [
#       "AlongTrack_MERRA-IOWA-CALIOP_2021-08-28T14-33-47ZD_AGU.nc",
#       "AlongTrack_MERRA-IOWA-CALIOP_2020-09-10T21-04-23ZD_AGU.nc"]


   pth_fig_out

   for file in filenames:
      ### Get date to get folder name
      yyyy=file[29:33];mm=file[34:36];dd=file[37:39]
      str_folder = yyyy+'-'+mm+'-'+dd
      if current_os == 'Linux' and "discover" in computer_name:  # Discover
         root_pth = pth_cal_root + yyyy+"/"+mm+"/"
      else: #Mac laptop
      # root_pth="/Volumes/ExtData1/SatData/Calipso/Level2/2019/08/2019-08-14/"
      # hdf_file_path = "CAL_LID_L2_05kmAPro-Standard-V4-20.2019-08-14T14-08-35ZD.hdf"
         root_pth = pth_cal_root + yyyy + "/" + str_folder +"/"
      
      # hdf_file_path = "CAL_LID_L2_05kmAPro-Standard-V4-20.2019-08-14T14-08-35ZD.hdf"
      hdf_file_path = "CAL_LID_L2_05kmAPro-Standard-V4-51." + file[29:48] + "ZD.hdf"
      # Specify the array names to load
      var_names = ['Latitude', 'Longitude','Lidar_Data_Altitudes','Profile_Time', 'Total_Backscatter_Coefficient_532','Extinction_Coefficient_532']
      in_data = load_hdf_arrays(root_pth+hdf_file_path,var_names)
      ## get original filename
      # print('Loaded data from original file:       \n',in_data['filename'].data.tobytes().decode('utf-8'))
      if in_data:
          print("\nLoaded arrays:")
          for array_name, array_data in in_data.items():
              print(f"    {array_name}: Shape {array_data.shape}")

   # sys.exit()
      ### reanme just compabitility with old code.
      # Define the new variable names
      new_var_names = ['lat', 'lon', 'alt','Time','tbc532','ext532']
      # Check if the number of elements in new_var_names matches var_names
      if len(new_var_names) != len(var_names):
          raise ValueError(f"Number of elements in new_var_names ({len(new_var_names)}) does not match var_names ({len(var_names)})")
      # Create a mapping between old and new names
      name_mapping = dict(zip(var_names, new_var_names))
      # Create a new dictionary with updated keys
      updated_data = {}
      for old_name, array_data in in_data.items():
          new_name = name_mapping.get(old_name, old_name)
          updated_data[new_name] = array_data
      # Replace the original dictionary with the updated one
      in_data = updated_data
      # Print the updated arrays with new names
      print("\nRenamed arrays:")
      for array_name, array_data in in_data.items():
          print(f"    {array_name}: Shape {array_data.shape}")


      # #### old
      lat = in_data['lat'][:,1]# Select Lat/Lon of center pulse
      lon = in_data['lon'][:,1]# Select Lat/Lon of center pulse
      zkm = in_data['alt'] # equally space Altitudeds
      tbc532 = in_data['tbc532'] # in 1/km.str units
      ext532= in_data['ext532']  # in 1/km units
      # ####

      # sys.exit()
      # title_str=in_data['filename'].data.tobytes().decode('utf-8') # bizzare command but it works
      title_str=hdf_file_path
      print('Data from ' , title_str)
      colbar_str=var_names[5]
      strings_plot= [title_str,colbar_str]


      #### now read the ouput txt file and get the rest of the data
      #### set ranges to display
      lat_min=lat.min();lat_max=lat.max()
      # lat_min=-50.5;lat_max=50.5;

      z_min=0;z_max=30.1;

      #subset indexes for subsetsin z and latitude dimensions
      idx_lat = np.where((lat >= lat_min) & (lat <= lat_max))
      idx_subset_lat = idx_lat[0]
      idx_z = np.where((zkm >= z_min) & (zkm<= z_max))
      idx_subset_z=len(zkm)-idx_z[0]

      ### get all subsets
      subset_lat = np.array([lat[i] for i in idx_subset_lat])
      subset_lon = np.array([lon[i] for i in idx_subset_lat])# otherwise the output is a 2D array
      ##### sub_back532=data.tab532[0][idx_subset_lat,idx_z]
      # sub_back532=tab532[idx_subset_lat[0]:idx_subset_lat[-1],idx_subset_z[0]:idx_subset_z[-1]]
      sub_ext532=1e-3*ext532[idx_subset_lat[:,None],idx_subset_z] # same as above but fancier, convert to 1/m units
      sub_z=zkm[idx_z]

      # f_get_track(path2file,filename)
      #caliop_plot2d(param, lat, lon, z, window_ave=0, colorange=(1e-6, 1e-3))
      # Good examples here
      #https://www.analyticsvidhya.com/blog/2020/09/colormaps-matplotlib/
      #https://github.com/NASA-DEVELOP/VOCAL/blob/master/calipso/dat/calipso-backscatter.cmap
      #http://meteothink.org/examples/meteoinfolab/satellite/calipso.html
      #https://pyhogs.github.io/colormap-examples.html good for understanding differences pcolor, pcountour...
      f_caliop_plot2d(sub_ext532.T,subset_lat,subset_lon,sub_z,strings_plot,\
                       colorange=(1e-6, 1e-2),scale='log',window_ave=0,esY=1)
      # use np.nanmin/nanamax to set values
      # plt.figure();plt.imshow(tab532[0::10,idx_subset_z].T,cmap='plasma',vmin=1e-5,vmax=1e-3)

      ##### save figure
      if print_fig=='yes':
         output_filename = 'BEXT_'+file[29:50]+'.png'
         print('   Saving to ',pth_fig_out+output_filename)
         plt.savefig(pth_fig_out+output_filename)

      ##### create map of satellite track
      # title = "Satellite Track, file " +  title_str
      # Call the function to plot the track
      # plot_track(lat, lon, title)
