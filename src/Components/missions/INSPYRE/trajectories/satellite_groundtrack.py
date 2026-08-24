#!/usr/bin/env python3
"""
Generate predicted satellite ground tracks from current TLEs.

Dependencies:
    pip install sgp4 astropy requests netCDF4 numpy

Inputs:
    days, timestep (minutes), satellite name
"""

import argparse
import requests
import numpy as np
from datetime import datetime, timedelta, timezone

from sgp4.api import Satrec, jday
from astropy.time import Time
from astropy.coordinates import TEME, ITRS, CartesianRepresentation, EarthLocation, AltAz, get_sun
import astropy.units as u
import os
#from netCDF4 import Dataset
import csv
from pathlib import Path

# ------------------------------------------------------------
# Satellites by NORAD ID
# ------------------------------------------------------------

SATELLITES = {
    "Suomi-NPP": 37849,
    "NOAA-20": 43013,
    "NOAA-21": 54234,
    "Terra": 25994,
    "Aqua": 27424,
    "Sentinel-5P": 42969,
    "EarthCARE": 59908,
    "PACE": 58928,
}


# ------------------------------------------------------------
# Parse Input Args
# ------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(
        description="Generate satellite ground tracks."
    )

    parser.add_argument(
        "--days",
        type=float,
        default=1.0,
        help="Number of days to propagate"
    )

    parser.add_argument(
        "--dt",
        type=float,
        default=5.0,
        help="Time step (minutes)"
    )

    parser.add_argument(
        "--satellite"
    )

    return parser.parse_args()

# ------------------------------------------------------------
# Download TLE
# ------------------------------------------------------------

def download_tle(
    catnr,
    cache_dir="./TLE/",
    force_download=False,
):
    """
    Load a satellite TLE from a local cache or download it from CelesTrak.

    Parameters
    ----------
    catnr : int or str
        NORAD catalog number.
    cache_dir : str or Path
        Directory used to store downloaded TLE files.
    force_download : bool
        If True, download a fresh TLE even when a cached file exists.

    Returns
    -------
    line1, line2 : tuple[str, str]
        The two TLE orbital-element lines.
    """

    catnr = str(catnr).strip()

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    tle_file = cache_dir / f"{catnr}.tle"

    # ---------------------------------------------------------------
    # Load the cached TLE
    # ---------------------------------------------------------------
    if tle_file.exists() and not force_download:
        try:
            lines = tle_file.read_text(encoding="utf-8").strip().splitlines()
            line1, line2 = _extract_tle_lines(lines)

            print(f"Loaded cached TLE for NORAD {catnr}: {tle_file}")
            return line1, line2

        except (OSError, RuntimeError) as exc:
            print(
                f"Cached TLE for NORAD {catnr} is invalid: {exc}. "
                "Downloading a fresh copy."
            )

    # ---------------------------------------------------------------
    # Download a fresh TLE
    # ---------------------------------------------------------------
    url = (
        "https://celestrak.org/NORAD/elements/gp.php"
        f"?CATNR={catnr}&FORMAT=TLE"
    )

    response = requests.get(url, timeout=100)
    response.raise_for_status()

    text = response.text.strip()

    if not text:
        raise RuntimeError(f"No TLE returned for NORAD {catnr}")

    lines = text.splitlines()
    line1, line2 = _extract_tle_lines(lines)

    # Save the complete response, including the satellite name when provided
    tle_file.write_text(text + "\n", encoding="utf-8")

    print(f" - Downloaded and cached TLE for NORAD {catnr}: {tle_file}")

    return line1, line2


def _extract_tle_lines(lines):
    """
    Extract and validate TLE line 1 and line 2.

    Handles either:
        satellite name
        line 1
        line 2

    or:
        line 1
        line 2
    """

    lines = [line.strip() for line in lines if line.strip()]

    line1 = next(
        (line for line in lines if line.startswith("1 ")),
        None,
    )
    line2 = next(
        (line for line in lines if line.startswith("2 ")),
        None,
    )

    if line1 is None or line2 is None:
        raise RuntimeError("The file does not contain valid TLE lines")

    # Confirm that both lines describe the same NORAD object
    if line1[2:7].strip() != line2[2:7].strip():
        raise RuntimeError(
            "TLE line 1 and line 2 have different NORAD catalog numbers"
        )

    return line1, line2

# ------------------------------------------------------------
# Time vector
# ------------------------------------------------------------

def create_times(start, stop, interval_minutes):

    times = []

    current = start

    while current < stop:
        times.append(current)
        current += timedelta(minutes=interval_minutes)

    return times

# ------------------------------------------------------------
# Propagation
# ------------------------------------------------------------

def propagate_satellite(line1, line2, times):
    """
    Propagate a satellite and calculate its ground position and
    local solar elevation.

    Returns
    -------
    lat : ndarray
        Satellite subpoint latitude, degrees.

    lon : ndarray
        Satellite subpoint longitude, degrees.

    solar_elevation : ndarray
        Solar elevation at the satellite subpoint, degrees.
        Values greater than zero indicate daytime.
    """
    sat = Satrec.twoline2rv(line1, line2)

    lat = np.full(len(times), np.nan, dtype=float)
    lon = np.full(len(times), np.nan, dtype=float)
    solar_elevation = np.full(len(times),np.nan,dtype=float)

    for i, current_time in enumerate(times):

        jd, fr = jday( current_time.year, 
                       current_time.month,
                       current_time.day,
                       current_time.hour,
                       current_time.minute,
                       current_time.second
                       + current_time.microsecond / 1.0e6)

        error, position, velocity = sat.sgp4(jd, fr)

        if error != 0:
            continue

        obstime = Time(current_time)

        teme = TEME(
            CartesianRepresentation(
                np.asarray(position) * u.km
            ),
            obstime=obstime,
        )

        itrs = teme.transform_to(ITRS(obstime=obstime))

        satellite_location = EarthLocation(x=itrs.x,y=itrs.y,z=itrs.z)

        lat[i] = satellite_location.lat.degree
        lon[i] = satellite_location.lon.degree

        # Use the satellite subpoint at ground level when calculating
        # the local solar elevation.
        ground_location = EarthLocation.from_geodetic(
            lon=lon[i] * u.deg,
            lat=lat[i] * u.deg,
            height=0.0 * u.m,
        )

        sun_altaz = get_sun(obstime).transform_to(
            AltAz(obstime=obstime,location=ground_location)
        )

        solar_elevation[i] = sun_altaz.alt.degree

    return lat, lon, solar_elevation
# --------------------------------------------------------------
# Split data by UTC day
# --------------------------------------------------------------

def split_by_day(times, lat, lon):

    daily_data = {}

    for i, t in enumerate(times):

        date = t.strftime("%Y%m%d")

        if date not in daily_data:
            daily_data[date] = {
                "times": [],
                "lat": [],
                "lon": []
            }

        daily_data[date]["times"].append(t)
        daily_data[date]["lat"].append(lat[i])
        daily_data[date]["lon"].append(lon[i])

    return daily_data

#--------------------------------------------------------------
#Write to csv
#--------------------------------------------------------------

def write_csv(filename, times, lat, lon):
    file_path = Path(filename) 
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with open(filename, 'w', newline='') as f:
        
        writer = csv.writer(f)
        
        writer.writerow([
            'Timestamp',
            'Latitude',
            'Longitude'
        ])

        for i, t in enumerate(times):
            writer.writerow([t,lat[i],lon[i]])

# ------------------------------------------------------------
def get_daytime_ground_track(
    satellite,
    start,
    stop,
    interval_minutes=5.0,
    min_solar_elevation=0.0,
    force_download=False,
    TLE_dir = './TLE/'
):
	"""
	Generate the daytime ground track for one satellite.
	
	Parameters
	----------
	satellite : str
		Satellite name listed in SATELLITES.
	
	start, stop : datetime
		Time range. UTC-aware datetime objects are preferred.
	
	interval_minutes : float
		Orbit propagation interval.
	
	min_solar_elevation : float
		Minimum solar elevation defining daytime. The normal
		daytime definition is zero degrees.
	
	Returns
	-------
	times : ndarray of datetime
	lat : ndarray
	lon : ndarray
	solar_elevation : ndarray
	daytime : ndarray of bool
	"""
	if satellite not in SATELLITES:
		raise ValueError(
			"Unknown satellite {!r}. Available satellites: {}".format(
				satellite,
				", ".join(SATELLITES.keys()),
			)
		)
	
	if start.tzinfo is None:
		start = start.replace(tzinfo=timezone.utc)
	
	if stop.tzinfo is None:
		stop = stop.replace(tzinfo=timezone.utc)
	
	times = create_times(
		start=start,
		stop=stop,
		interval_minutes=interval_minutes,
	)
	
	catnr = SATELLITES[satellite]
	
# 	line1, line2 = download_tle(catnr)
	line1, line2 = download_tle(catnr=catnr,
	                            cache_dir=TLE_dir,
	                            force_download=force_download)	
	                            
	lat, lon, solar_elevation = propagate_satellite(
		line1=line1,
		line2=line2,
		times=times,
	)
	
	times = np.asarray(times, dtype=object)
	
	daytime = (
		np.isfinite(lat)
		& np.isfinite(lon)
		& np.isfinite(solar_elevation)
		& (solar_elevation >= min_solar_elevation)
	)

	output = {}
	output['times'] = times
	output['lat'] = lat
	output['lon'] = lon
	output['solar_elevation'] = solar_elevation
	output['daytime'] = daytime

	return output
    
# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    args = parse_args()

    start = datetime.now(timezone.utc).replace(hour=0,minute=0,second=0,microsecond=0)

    stop = start + timedelta(days=args.days)

    interval = args.dt
    
    sat_name = args.satellite

    times = create_times(start, stop, interval)

    ntime = len(times)

    latitude = np.zeros(ntime)
    longitude = np.zeros(ntime)
    
    print(f"Processing {sat_name}")

    catnr = SATELLITES[sat_name]

    line1, line2 = download_tle(catnr)

    lat, lon = propagate_satellite(
            line1,
            line2,
            times,
    )

    latitude = lat
    longitude = lon

    #filename = './Orbits/' + start.strftime('%Y%m%d') + '/'+sat_name + '_' + start.strftime('%Y%m%d') + '.csv'
    #write_csv(filename, times, latitude, longitude)
    daily_data = split_by_day(
        times,
        latitude,
        longitude
    )
	
    os.mkdir('./Orbits/', exist_ok=True)
    for date, data in daily_data.items():

        filename = (
            './Orbits/'
            + 'GeneratedOn' + start.strftime('%Y%m%d')
            + '/'
            + sat_name
            + '_'
            + date
            + '.csv'
        )

        write_csv(
            filename,
            data["times"],
            data["lat"],
            data["lon"]
        )

        print(f"Wrote {filename}")
        print("Finished.")


if __name__ == "__main__":
    main()
