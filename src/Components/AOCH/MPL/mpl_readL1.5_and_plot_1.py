#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Feb/24/2025 - Same as v2 but the plotting routine is now inside a module. it Works!

This script reads MPLnet data into memory:

A NetCDFData class to store the data and attributes

A read_netcdf function that:
    Takes a filename (required) and
    (optional) list of variables to read (default variables if none specified)
    Reads data & attributes for each variable
    Stores everything in the class structure
    Verb=0(default), no print,=1 variables names and sizes, =2 var names and attributes

Main Subroutines:
    read_file : actual read single file
    get_path  : set path where files are stored (needed if working in different clusters)


Version2 (Feb/19/2025) same as V1 with the addition of simple plot creationg

Created on Fri Jan 17 10:51:10 2025

@author: sgasso + chat.gsfc
"""
from netCDF4 import Dataset
import numpy as np
import platform
import sys
from datetime import datetime, timedelta # need to get date from filename

# import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#from matplotlib.colors import ListedColormap
from scipy import ndimage


class MPL_L15:
    def __init__(self):
        """Initialize the class with empty attributes"""
        self.time = None
        self.time_attributes = {}

        self.altitude = None
        self.altitude_attributes = {}

        self.extinction = None
        self.extinction_attributes = {}

        self.extinction_err = None
        self.extinction_err_attributes = {}

        self.qa_extinction = None
        self.qa_extinction_attributes = {}

    @classmethod
    def read_file(cls, filename, variables=None, verb=0):
        """
        Read data from NetCDF file

        Parameters:
        filename (str)  : Path to the NetCDF file
        variables (list): List of variables to read. If None, reads default variables
        verb (int): Verbosity level
                   0: No printing
                   1: Print variable name, type, and dimensions
                   2: Print level 1 info plus attributes

        Returns:
        NetCDFData: Class instance containing the requested data and attributes
        """
        # Default variables to read if none specified
        default_variables = ['time', 'altitude', 'extinction', 'extinction_err']

        # Use default variables if none provided
        variables_to_read = variables if variables is not None else default_variables

        # Initialize the class
        instance = cls()

        try:
            # Open the NetCDF file
            with Dataset(filename, 'r') as nc:
                # Read each requested variable
                for var_name in variables_to_read:
                    if var_name in nc.variables:
                        # Get the variable
                        var = nc.variables[var_name]

                        # Read the data
                        var_data = var[:]

                        # Read all attributes
                        var_attributes = {attr: var.getncattr(attr)
                                       for attr in var.ncattrs()}

                        # Print information based on verbosity level
                        if verb >= 1:
                            print("\nVariable Information:")
                            print(f"Name: {var_name}")
                            print(f"Type: {var.dtype}")
                            print(f"Dimensions: {var.dimensions}")
                            print(f"Shape: {var.shape}")

                            if verb >= 2:
                                print("\nAttributes:")
                                for attr, value in var_attributes.items():
                                    print(f"  {attr}: {value}")
                            print("=" * 50)

                        # Store data and attributes in the class
                        setattr(instance, var_name, var_data)
                        setattr(instance, f"{var_name}_attributes", var_attributes)
                    else:
                        print(f"Warning: Variable {var_name} not found in file")

        except Exception as e:
            print(f"Error reading NetCDF file: {str(e)}")
            return None

        return instance

    def get_variable_names(self):
        """Return list of available variables"""
        return [attr for attr in self.__dict__.keys() if not attr.endswith('_attributes')]

    def get_variable_info(self, variable_name):
        """Print information about a specific variable"""
        if hasattr(self, variable_name):
            data = getattr(self, variable_name)
            attrs = getattr(self, f"{variable_name}_attributes", {})
            print(f"\nVariable: {variable_name}")
            print(f"Shape: {data.shape}")
            print(f"Type: {data.dtype}")
            print("Attributes:", attrs)
        else:
            print(f"Variable {variable_name} not found")

def get_path(now_os, now_computer):
    # Your existing get_path function remains the same
    if now_os=='win32': # PC
        base_path    = 'C:/Users/sgasso/Downloads/'
        pth_fig_out='C:/Users/sgasso/OneDrive - NASA/ToShare/2025/GEOS/Pyfigures/'
    elif now_os == 'darwin': # My Laptop
#        base_path='/Users/sgasso/Downloads/'
        base_path='/Volumes/ExtData1/Surface/MPL/V3_L15/'
        pth_fig_out='/Users/sgasso/Library/CloudStorage/OneDrive-NASA/ToShare/2025/AOCH/PyFigures/'
    elif now_os == 'linux' and "calculon" in now_computer:
        base_path = '/nobackup/CALIPSO/Level1.5/Santiago/'
    elif now_os == 'linux' and "discover" in now_computer:
        base_path = '/nobackup/CALIPSO/Level1.5/Santiago/'
    else:
        print('Current operating system no recognized.')
        print('Cannot set path to MPL  files. Terminate Run')
        sys.exit()

    return base_path,pth_fig_out

#### Read list of MPL sites
def read_mpl_site_list(filename):
#Expected format: comma separated
#line 1 :Name_MPL_sites , name,lat, lon, altiude(km)
#line 2 : GSFC           ,38.9930,-76.8400,0.050
#.............
#
    sites = []
    with open(filename, 'r') as file:
        header = file.readline()
        for line in file:
            data = [item.strip() for item in line.split(',')]
            site = (
                data[0],
                float(data[1]),
                float(data[2]),
                float(data[3])
            )
            sites.append(site)
    return sites #this is a tuple


def f_mpl_plot2d(param, xtime, z, title_strings, colorange=(1e-6, 1e-3), scale='lin', window_ave=0, esY=10):
    """
    Create 2D plot of MPL data

    Parameters:
    param (array): Data to plot (extinction)
    xtime (array): Time array
    z (array): Height array
    title_strings (list): [title_str, colbar_str]
    colorange (tuple): Min and max values for colorbar
    scale (str): Scale type ('lin' or 'log')
    window_ave (int): Window size for averaging
    esY (int): Spacing for time axis
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    from scipy import ndimage

    plot_title = title_strings[0]
    colorbar_title = title_strings[1]
    esX = 1  # Fixed X-axis spacing for height dimension

    # Data preparation
    lastZ = len(param[0,:])
    var0 = param[::esX,::esY]
    plot_time = xtime[::esY]
    plot_z = z[::-esX]

    # Apply filtering
    if window_ave == 0:
        var0[var0 < 1e-12] = np.nan
        var = var0
    else:
        var = ndimage.uniform_filter(var0, size=window_ave, mode='constant')
        var[var <= 0] = np.nan

    # Prepare plot variables
    plot_var = np.log10(var)
    vmin = np.log10(colorange[0])
    vmax = np.log10(colorange[1])
    space_col = 30

    # Color setup
    base_cmap = plt.get_cmap('jet')
    levels = np.linspace(vmin, vmax, space_col + 1)
    calipso_colors = [base_cmap(i) for i in np.linspace(0, 1, space_col-10)]
    gray_colors = [colors.to_rgba(str(min(i/10, 1))) for i in range(1, 12)]
    calipso_colors.extend(gray_colors)
    cmap = colors.ListedColormap(calipso_colors)
    norm = colors.BoundaryNorm(levels, len(calipso_colors))

    # Create plot
    fig, ax = plt.subplots(figsize=(16, 6))
    plt.pcolor(plot_var, norm=norm, cmap=cmap)
    cbar = plt.colorbar(extend='max')

    # Colorbar setup for log scale
    if scale == 'log':
        bar_ticks = np.array([1e-6, 2e-6, 3e-6, 5e-6, 7e-6,
                             1e-5, 2e-5, 3e-5, 5e-5, 7e-5,
                             1e-4, 2e-4, 3e-4, 5e-4, 7e-4,
                             1e-3, 2e-3, 3e-3, 5e-3, 7e-3, 1e-2])
        log_bar = np.log10(bar_ticks)
        log_bar_with_ends = np.concatenate(([log_bar[0]], log_bar, [log_bar[-1]]))
        cbar.ax.minorticks_off()
        cbar.set_ticks(log_bar_with_ends)
        tick_labels = [f"{10**label:.1e}" for label in log_bar_with_ends]
        cbar.set_ticklabels(tick_labels)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label(colorbar_title, rotation=270, labelpad=15)

    # Axis setup
    xt = np.arange(0, 24) * 60
    yt = np.arange(len(plot_z))
    if yt[0] != 0:
        yt = np.insert(yt, 0, 0)
    if yt[-1] != len(plot_z) - 1:
        yt = np.append(yt, len(plot_z) - 1)
    nYTicks = 20

    iyt = np.round(np.linspace(0, len(yt) - 1, nYTicks)).astype(int)

    # X axis labels
    xvals = []
    xstrs = []
    for i in xt:
        xvals.append(i)
        xstrs.append('%s' % int(i/60))

    # Set axis labels and ticks
    plt.xticks(xvals, xstrs, fontsize=8)
    plt.yticks(iyt[::-1], [f"{c:.1f}" for c in plot_z[iyt]], fontsize=8)

    # Final labels and title
    plt.ylabel('Height (km)')
    plt.xlabel('Time (UTC)')
    if window_ave == 0:
        plt.title(plot_title)
    else:
        plt.title(f"{plot_title} , Moving average = {window_ave}")

    return fig, ax

######----------------------------------------------


if __name__ == "__main__":
    # Option 1: Direct call to the class method
    current_os    = platform.system()
    computer_name = platform.node()
    current_pth,pth_fig = get_path(current_os.lower(),computer_name)
    filename = "MPLNET_V3_L15_AER_20230820_MPL44258_GSFC.nc4"
    full_pathname = current_pth + filename

    # Define variables to read
    variables_to_read = ['time', 'altitude', 'extinction']

    # Read the data
    data = MPL_L15.read_file(full_pathname, variables=variables_to_read, verb=0)

    # ### now print some misc data directly
    # print(data.time.shape)
    # print(data.altitude.shape)
    # print(data.extinction.shape)
    # # Access attributes
    # print(data.time_attributes)
    # print(data.extinction_attributes.get('units'))

    #------------------------------------------------
    ### now do a quick plot along the track
    ###
    ### Get date from file
    # time is in the Julian Date format of number of days since January 1st, 4713 BC , noon UTC
    # where integer par is Ndays since this date and decimal time is hh:mm:ss since noon UTC
    # so for example : 2440588.000000 is A.D. 1970 January 1	12:00:00.0
    # this number can be substracted from data.time and use it as new reference point
    # altenrnativel the python module astropy can deal with it. NOTE the python module
    # datetime does not deal with this reference point .
    # I used a different approach based on the date in the filename.

    x     =data.time
    yyyy,mm,dd=[filename[18:22],filename[22:24],filename[24:26]]
    base_date = datetime(int(yyyy),int(mm),int(dd),0,0,0) -timedelta(days=0.5) # Your base date
    xtime = []
    for jd in x:
        # Get just the fractional part since we know the date
        day_fraction = jd % 1
        # Convert fractional day to minutes
        #seconds = int(day_fraction * 86400)
        minutes = int(day_fraction * 1440)
        # Create datetime by adding the time to base date
        dt = base_date + timedelta(minutes=minutes)
        xtime.append(dt)
    ### Because MPL data is reported every 60 min,
    ### one could create a time array of size 60*24 = 1440 and it should be the same.

    ### Assign inputs to module
    site_name = filename[36:-4]
    z_in    =data.altitude[0,:]
    varin   =1e-3*data.extinction[0,:].T # 1440 x 400 to 400 x 1440
    strings_plot= ['MPL Extinction 532nm, Site ' + site_name,'1/m']
    ### Plot
    f_mpl_plot2d(varin, xtime, z_in, strings_plot,
             colorange=(1e-6, 1e-2), scale='log', window_ave=0, esY=1)