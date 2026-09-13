# GigaTraj Python Interface

This directory has modules implementating a Python interface to the C++ trajectory code ``GigaTraj`` developed by Leslie Lait . Although much of the code here is general purpose, what you find below is motivated by use cases encountered during NASA's [INSPYRE](https://www-air.larc.nasa.gov/missions/inspyre/index.html) airborne campaign. Some basic clean up is necessary before this becomes a truly general purpose Python interface to ``GigaTraj``.

## About  ``GigaTraj`` C++ Library and Application

``GigaTraj`` is a C++ software library, plus some sample programs, for calculating the trajectories of air parcels in the lower stratosphere. Its purpose is purely for scientific research. It allows for several different sources of meteorological data, and it can calculate either kinematic (in pressure coordinates using omega) or quasi-isentropic (in potential temperature coordinates using heating rates) trajectories.

Please note that ``GigaTraj`` is not a general Lagrangian transport model. It does not attempt to simulate the effects of turbulence or small-scale diffusion. This makes it unsuitable for use in regions where those effects are important--for example, at low altitudes.


``Gigatraj`` is available on [GitHub](https://github.com/nasa/gigatraj). 


## Main module and configuration


- ``gigatraj.py:`` Defines the main class GIGATRAJ with methods for generating parcel initial locations and actually running the ``gigatraj`` C++ code. By design, this class does not write any files and does not include any visualization functionality. Notice that the `GigaTraj` C++ applications invoked internally do write NetCDF files. Therefore, methods such as ``genParcel()`` and ``genTrajectories`` produce NetCDF file output.

- ``gigatraj.yaml:`` This is the main YAML configuration driving the executionn of ``GigaTraj``. The main sections are:

  * ``Fires:`` defines the sites where to launch trajectories from (refered to as **fires**, an INSPYRE bias). Probably this should be renamed ``Sites``.
  * ``Parcels:`` where the initial location of parcels are defined. Parcels are randomly generated around a 3D location given by the triplet ``(latitude,longitude,altitude)``. Notice that this is not the only way to initialize the parcels. See utility ``erthcare_particles.py`` below for another way to sample particles based on lidar backscatter.
  * ``Trajectories:`` controls the calculation of the trajectories: trajectory duration and other configuration parameters need by the C++ app. Examples includes the particular forecast cycle to be used, and which meteorological product to use (GEOS or IFS).
    * ``Releases:`` here is where the initial date and particular sites to start trajectories are defined. You can have multiple releases, each with its own starting date and fires.

  * ``Plots:`` specific configuration for plotting (not used in main class ``GIGATRAJ``). Notice that there is an additional YAML file called ``inspyre.yaml`` where you can specify campaign specific aircrafts, range radii, and airports; this is particulalry important for configuring your plots during suitcase flights.

- Meteorological data catalog needed by ``GigaTraj``. These are the configuration files specifying the file locations and variable names for driving the trajectory code. There are versions specific to GEOS and ECMWF IFS forecasts:
  * ``MetGEOS.cat:`` customization for GEOS at AWS, running on the parallel cluster ``(pcluster)`` associated with NASA's Airborne SMCE. 
  * ``MetIFS.cat:`` corresponding custorization for IFS forecasts
  * ``MetGEOSfp.cat:`` this is the generic name expected by the ``GigaTraj`` C++ application. Although you should be able to specifiy the name of the catalog on the command line, I could never get it to work. As a stop gap measure, this file is symlink'ed to ``MetGEOS.cat`` or ``MetIFS.cat`` when running with each of those meteorological fields.

## Other python supporting modules

- ``earthcare.py:`` handy utility for finding EarthCARE orbits and frames. Useful when making simulations for EarthCARE orbits/frames as seen in [JAXA's EarthCARE Quick Looks](https://www.eorc.jaxa.jp/EARTHCARE/Quicklook). See this [Notebook](https://github.com/Asticou-Earth/AeroApps/blob/develop/src/Components/missions/INSPYRE/Notebooks/Sat_plane.ipynb) for an example.
- ``earthcare_particles.py:`` utility for generating particle initial conditions based on Quick Look images. This code will eventually have an option to work with Level 2 HDF-5 as well. See this [Notebook](https://github.com/Asticou-Earth/AeroApps/blob/develop/src/Components/missions/INSPYRE/Notebooks/earthcare_particles_example.ipynb) for an example.
- ``sat_plane.py:`` utility to find the intersection of trajectories with specific EarthCARE orbit/frame. See this [Notebook](https://github.com/Asticou-Earth/AeroApps/blob/develop/src/Components/missions/INSPYRE/Notebooks/Sat_plane.ipynb) for an example.
- ``sat_tracks.py:`` this is a refactoring of module ``satellite_groundtracks``. For a set of satellites, built in (easily extended if you know the NORAD catalog number for the bird), it can download Two-Line Element (TLE) files and compute the satellite ground tracks given a time sequence.
- ``traj_plot.py:`` this is the workhorse for plotting the trajectories on a map based on ``cartopy.`` This is a refactoring of ``traj_make_plot.py`` and ``plot_parcel_forecast_density.py``  to remove file I/O and hardwired file name conventions. Now it only includes the basic plotting functionality.

## Driving scripts used during INSPYRE 


- ``do_alltrajs:`` taking configuration information from ``gigatraj.yaml``, this script computes the full trajectory workflow: initial conditions, trajectory calculation and plotting. This process is parallelized using Python's ``ThreadPool`` module.

- ``plot_allfires:`` while ``do_alltrajs`` can do some basic plotting for a given map region, this script allows a more comprehensive set of plots to be generated, in parallel. It works by modifying the YAML configuration in memory.

- ``plot_onefire:`` this is a helper script to allow parallelization with ``ThreadPool``. This is needed because ``matplotlib`` and ``cartopy`` are notoriously not thread safe. See this [Notebook](https://github.com/Asticou-Earth/AeroApps/blob/develop/src/Components/missions/INSPYRE/Notebooks/Merge_fires.ipynb) for an example of how to merge multiple fires in a single plot.

## Generic download scripts

These scripts are needed for downloading the necessary meteorological data when running on AWS. On Discover, the GEOS files are already present on disk and no data download is needed. However, IFS files are not routnely available on Discover and they will need to be downloaded.

- ``download_url:`` generic script for downloading files from an URL with options to skip some files and limit the number of files to download.

- ``download_fp:`` download GEOS-FP forecasts from NCCS

- ``download_ana:`` download GEOS assimilation fields

- ``download_ifs:`` download IFS forecasts from ECMWF; it also symlinks files in the analysis directory ``(diag/)`` to enable calculation of trajectories released before the forecast initial time. This script downloads GRIB2 files from ECMWF and convert them to be GEOS-like NetCDF files using CDO and NCO utilities. These utilities can be installed from ``conda-forge.``

- ``download_earthcare_quicklooks:`` download images from JAXA's site.

## Python Environment
Python dependencies can be foundn in ``environment.yaml``. To create a dedicated conda environment do this:
```
conda env create -f environment.yaml -n inspyre
conda activate inspyre
```

## Typical workflow during INSPYRE 2026 Deployment

1. Start the day by downloading the latest forecasts and analysis using the download scripts above.
2. Update ``gigatraj.yaml`` with the addition of new fires to the ``Fires`` dictionary, and enter the desired particle release dates/specific fires under ``Trajectories``. 
   - **Note**: Make sure ``maxConcurrent`` is consistent with the numbers of cores available.
3. During suitcase flights, also update the airports in ``inspyre.yaml``.
4. Customize ``plot_allfires``, updating regions and valid times to create multiple density plots.

## Credit

[Arlindo da Silva](mailto:arlindo.dasilva@asticou-earth.com), Asticou Earth Systems, designed/wrote/refactored this code during the INSPYRE 2026 deployment. Some of the refactoring aided by CODEX.



