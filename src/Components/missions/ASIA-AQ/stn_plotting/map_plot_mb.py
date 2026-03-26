import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# Define the custom colormap matching the image
# white -> yellow -> orange -> red -> blue -> dark blue
colors = [
    '#FFFFFF',  # white       (lowest density)
    '#FFFFCC',  # light yellow
    '#FFEDA0',  # yellow
    '#FED976',  # light orange
    '#FEB24C',  # orange
    '#FD8D3C',  # dark orange
    '#FC4E2A',  # orange-red
    '#E31A1C',  # red
    '#BD0026',  # dark red
    '#053061',  # dark blue   (highest density)
]

custom_cmap = mcolors.LinearSegmentedColormap.from_list(
    'density_cmap', colors, N=256
)


def plot_station_map_mb(df_list, spcname, lon_col='Longitude', lat_col='Latitude', 
                     val_col='PM25', titles=['','',''], units=r'$\mu g/m^3$',
                     cmap=plt.cm.gist_earth_r, vmin=None, vmax=None,
                     figsize=(9, 12), outdir=None):
    """
    Plot a map scatter plot where the color of each point 
    represents the absolute mean value at each station.

    Parameters
    ----------
    df       : pandas DataFrame containing station data
    lon_col  : column name for longitude
    lat_col  : column name for latitude
    val_col  : column name for the variable to plot
    title    : plot title
    units    : units string for colorbar label
    cmap     : colormap name
    vmin     : minimum value for colorbar (default: data min)
    vmax     : maximum value for colorbar (default: data max)
    figsize  : figure size tuple
    outfile  : if provided, save figure to this path
    """

    # --- Compute absolute mean at each station ---
    # --- Compute absolute mean at each station for each dataframe ---
    station_means = []
    for df in df_list:
        stn_mean = (df.groupby(['station', lat_col, lon_col])[val_col]
                      .apply(lambda x: np.nanmean(x))
                      .reset_index()
                      .rename(columns={val_col: 'abs_mean'}))
        station_means.append(stn_mean)

    # --- Set colorbar limits ---
    if vmin is None:
#        vmin = min(sm['abs_mean'].min() for sm in station_means)
        vmin = 0
    if vmax is None:
#        vmax = max(sm['abs_mean'].max() for sm in station_means)
        vmax = 10

    # --- Compute map extent across all 3 dataframes ---
    all_lons = np.concatenate([sm[lon_col].values for sm in station_means])
    all_lats = np.concatenate([sm[lat_col].values for sm in station_means])
    lon_pad  = 5
    lat_pad  = 5
    extent   = [all_lons.min() - lon_pad, 
                #all_lons.max() + lon_pad,
                -60,
                all_lats.min() - lat_pad, 
                all_lats.max() + lat_pad]


    # --- Create figure with 3 subplots ---
    fig, axes = plt.subplots(3, 1,
                             figsize=figsize,
                             subplot_kw={'projection': ccrs.PlateCarree()})

    # --- Plot each dataframe ---
    for idx, (ax, station_mean, title) in enumerate(zip(axes, station_means, titles)):

        # --- Add map features ---
        ax.add_feature(cfeature.LAND,       facecolor='lightgray', alpha=0.3)
        ax.add_feature(cfeature.OCEAN,      facecolor='lightblue', alpha=0.3)
        ax.add_feature(cfeature.COASTLINE,  linewidth=0.5)
        ax.add_feature(cfeature.BORDERS,    linewidth=0.5, linestyle='--')
        ax.add_feature(cfeature.STATES,     linewidth=0.3, linestyle=':')

        # --- Set map extent ---
        ax.set_extent(extent, crs=ccrs.PlateCarree())

        # --- Plot scatter points ---
        sc = ax.scatter(station_mean[lon_col],
                        station_mean[lat_col],
                        c=station_mean['abs_mean'],
                        cmap=cmap,
                        vmin=vmin,
                        vmax=vmax,
                        s=40,
                        edgecolors='black',
                        linewidths=0.5,
                        transform=ccrs.PlateCarree(),
                        zorder=5)

        # --- Add gridlines ---
        gl = ax.gridlines(draw_labels=True, linewidth=0.5,
                          color='gray', alpha=0.5, linestyle='--')
        gl.top_labels   = False
        gl.right_labels = False

        # --- Add title ---
        ax.set_title(title, fontsize=13, fontweight='bold', pad=10)


    # --- After the loop, add colorbar to the right of the last axis ---
    # Create a new axes on the right side of the figure
    fig.subplots_adjust(right=0.85)  # make room on the right

    # Add a new axes for the colorbar
    cbar_ax = fig.add_axes([0.88,   # left position
                             0.1,   # bottom position
                             0.02,  # width
                             0.8])  # height

    cbar = fig.colorbar(sc, cax=cbar_ax, orientation='vertical')
    cbar.set_label(f'Mean Bias {spcname} ({units})', fontsize=11)


#    plt.tight_layout()

    # --- Save or show ---
    if outdir is not None:
        outfile = f'{outdir}/{spcname}_map_mb.png'
        plt.savefig(outfile, dpi=150, bbox_inches='tight')
        print(f'  Saved: {outfile}')
    else:
        plt.show()

    plt.close()
