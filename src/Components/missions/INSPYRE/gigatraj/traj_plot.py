"""
Plot parcel trajectory forecasts from in-memory xarray datasets. This is a refactoring
of these packages by Pete Colarco and others:

- plot_parcel_forecast_density.py
- traj_make_plot.py

This module deliberately contains no trajectory-file discovery, NetCDF reading,
filename parsing, or satellite-orbit logic.  The caller owns trajectory data
loading and plot output; campaign geometry is loaded from a required YAML file.
This will hopefully make this code more PEP 8 compliant and reusable to other campaigns.

Arlindo da Silva, August 2026. Refactoring aided by CODEX.

"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from itertools import cycle
from pathlib import Path
from typing import Sequence

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.geodesic import Geodesic
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from scipy.stats import gaussian_kde
from shapely.geometry import Polygon
import xarray as xr
import yaml

from sat_tracks import SATELLITES, classify_daylight, propagate_ground_track


_COLOR_SCHEMES = (
    ("YlOrRd", "maroon"),
    ("Blues", "cornflowerblue"),
    ("YlGn", "lime"),
)
_AIRCRAFT_RING_STYLES = (
    ("black", 2),
    ("white", 3),
    ("tab:purple", 2),
    ("tab:cyan", 2),
    ("tab:orange", 2),
)
_SATELLITE_COLORS = {
    "EarthCARE": "crimson",
    "NOAA-20": "magenta",
    "PACE": "green",
}
_FALLBACK_SATELLITE_COLORS = (
    "tab:blue",
    "tab:orange",
    "tab:purple",
    "tab:brown",
    "tab:pink",
)
DEFAULT_GEOGRAPHIC_BOUNDS = (-120.0, -70.0, 22.5, 60.0)
_SATELLITE_FILTERS = {"day_only", "night_only", None}
_KDE_GRID_SIZE = 140
_KDE_PADDING_BANDWIDTHS = 3.0
_KDE_MINIMUM_RELATIVE_DENSITY = 0.01


def _as_utc_naive(value: datetime, name: str) -> datetime:
    if not isinstance(value, datetime):
        raise TypeError(f"{name} must be a datetime object")
    if value.tzinfo is not None:
        return value.astimezone(timezone.utc).replace(tzinfo=None)
    return value


def _altitude_km(dataset: xr.Dataset) -> float:
    if "altitude_km" not in dataset.attrs:
        raise ValueError("Every trajectory dataset must define attrs['altitude_km']")
    try:
        altitude = float(dataset.attrs["altitude_km"])
    except (TypeError, ValueError) as exc:
        raise ValueError("attrs['altitude_km'] must be numeric") from exc
    if not np.isfinite(altitude):
        raise ValueError("attrs['altitude_km'] must be finite")
    return altitude


def _validated_geographic_bounds(geographic_bounds):
    if not isinstance(geographic_bounds, Sequence) or isinstance(
        geographic_bounds, (str, bytes)
    ):
        raise TypeError(
            "geographic_bounds must be a sequence of "
            "[west, east, south, north]"
        )
    if len(geographic_bounds) != 4:
        raise ValueError(
            "geographic_bounds must contain exactly four values: "
            "[west, east, south, north]"
        )
    try:
        west, east, south, north = map(float, geographic_bounds)
    except (TypeError, ValueError) as exc:
        raise TypeError("geographic_bounds values must be numeric") from exc
    if not all(np.isfinite(value) for value in (west, east, south, north)):
        raise ValueError("geographic_bounds values must be finite")
    if not (-180.0 <= west <= 180.0 and -180.0 <= east <= 180.0):
        raise ValueError(
            "west and east must each be between -180 and 180 degrees"
        )
    if west == east:
        raise ValueError("west and east must not be equal")
    if not (-90.0 <= south < north <= 90.0):
        raise ValueError(
            "geographic_bounds must satisfy -90 <= south < north <= 90"
        )
    return west, east, south, north


def _trajectory_arrays(dataset: xr.Dataset):
    required = ("time", "lat", "lon", "PAlt")
    missing = [name for name in required if name not in dataset]
    if missing:
        raise ValueError(f"Trajectory dataset is missing variables: {', '.join(missing)}")

    time = dataset["time"]
    if time.ndim != 1:
        raise ValueError("The 'time' coordinate must be one-dimensional")
    time_dim = time.dims[0]

    variables = {name: dataset[name] for name in ("lat", "lon", "PAlt")}
    for name, variable in variables.items():
        if variable.ndim != 2 or time_dim not in variable.dims:
            raise ValueError(f"'{name}' must be two-dimensional and include the time dimension")

    particle_dims = [dim for dim in variables["lat"].dims if dim != time_dim]
    if len(particle_dims) != 1:
        raise ValueError("Unable to identify the parcel dimension")
    particle_dim = particle_dims[0]
    expected_dims = {time_dim, particle_dim}
    for name, variable in variables.items():
        if set(variable.dims) != expected_dims:
            raise ValueError("'lat', 'lon', and 'PAlt' must use the same dimensions")

    times = np.asarray(time.values).astype("datetime64[us]")
    lat = np.asarray(variables["lat"].transpose(time_dim, particle_dim).values, dtype=float)
    lon = np.asarray(variables["lon"].transpose(time_dim, particle_dim).values, dtype=float)
    altitude = np.asarray(
        variables["PAlt"].transpose(time_dim, particle_dim).values, dtype=float
    )
    if not (lat.shape == lon.shape == altitude.shape):
        raise ValueError("'lat', 'lon', and 'PAlt' must have identical shapes")

    # Cartopy accepts either convention, but the density calculation must not
    # straddle the -180/180 discontinuity.
    lon = np.mod(lon, 360.0)
    return times, lat, lon, altitude


def _fire_name(datasets: Sequence[xr.Dataset]) -> str:
    missing = [
        index
        for index, dataset in enumerate(datasets)
        if not str(dataset.attrs.get("fire", "")).strip()
    ]
    if missing:
        positions = ", ".join(str(index) for index in missing)
        raise ValueError(
            "Every trajectory dataset must define a non-empty attrs['fire']; "
            f"missing at Traj index(es): {positions}"
        )
    names = {str(dataset.attrs["fire"]).strip() for dataset in datasets}
    if len(names) != 1:
        raise ValueError(
            "All trajectory datasets must describe the same fire; found: "
            + ", ".join(sorted(names))
        )
    return names.pop()


def _release_time(datasets: Sequence[xr.Dataset]) -> datetime:
    release_times = []
    for index, dataset in enumerate(datasets):
        value = dataset.attrs.get("Trajectory_start")
        if value is None or not str(value).strip():
            raise ValueError(
                "Every trajectory dataset must define attrs['Trajectory_start']; "
                f"missing at Traj index {index}"
            )
        if isinstance(value, datetime):
            parsed = value
        else:
            try:
                parsed = datetime.fromisoformat(str(value).strip().replace("Z", "+00:00"))
            except ValueError as exc:
                raise ValueError(
                    f"Invalid attrs['Trajectory_start'] at Traj index {index}: {value!r}"
                ) from exc
        release_times.append(_as_utc_naive(parsed, "Trajectory_start"))

    unique_times = set(release_times)
    if len(unique_times) != 1:
        values = ", ".join(sorted(time.isoformat(timespec="minutes") for time in unique_times))
        raise ValueError(
            "All trajectory datasets must have the same Trajectory_start; found: "
            + values
        )
    return unique_times.pop()


def _sorted_release_datasets(datasets: Sequence[xr.Dataset]):
    altitude_dataset_pairs = [(_altitude_km(dataset), dataset) for dataset in datasets]
    altitudes = [altitude for altitude, _ in altitude_dataset_pairs]
    duplicate_altitudes = sorted(
        {altitude for altitude in altitudes if altitudes.count(altitude) > 1}
    )
    if duplicate_altitudes:
        values = ", ".join(f"{altitude:g}" for altitude in duplicate_altitudes)
        raise ValueError(
            "Each trajectory dataset must have a different release altitude; "
            f"duplicate altitude_km value(s): {values}"
        )
    return [dataset for _, dataset in sorted(altitude_dataset_pairs, key=lambda item: item[0])]


def _load_campaign(campaign_file):
    try:
        path = Path(campaign_file)
    except TypeError as exc:
        raise TypeError("CampaignFile must be a path to a YAML file") from exc
    if not path.is_file():
        raise FileNotFoundError(f"Campaign YAML file not found: {path}")

    try:
        with path.open("r", encoding="utf-8") as stream:
            campaign = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid campaign YAML in {path}: {exc}") from exc
    if not isinstance(campaign, dict):
        raise ValueError("Campaign YAML must contain a top-level mapping")

    markers_config = campaign.get("Markers")
    aircraft_config = campaign.get("Aircrafts")
    if not isinstance(markers_config, dict) or not markers_config:
        raise ValueError("Campaign YAML requires a non-empty 'Markers' mapping")
    if not isinstance(aircraft_config, dict) or not aircraft_config:
        raise ValueError("Campaign YAML requires a non-empty 'Aircrafts' mapping")

    markers = {}
    for name, coordinates in markers_config.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Every marker must have a non-empty string name")
        if not isinstance(coordinates, (list, tuple)) or len(coordinates) != 2:
            raise ValueError(f"Marker '{name}' must be [longitude, latitude]")
        try:
            longitude, latitude = (float(value) for value in coordinates)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Marker '{name}' coordinates must be numeric") from exc
        if not (-180 <= longitude <= 180 and -90 <= latitude <= 90):
            raise ValueError(f"Marker '{name}' coordinates are outside valid bounds")
        markers[name] = (longitude, latitude)

    aircraft = []
    for name, settings in aircraft_config.items():
        if not isinstance(name, str) or not name.strip() or not isinstance(settings, dict):
            raise ValueError("Every aircraft must have a name and a settings mapping")
        airport = settings.get("airport")
        if airport not in markers:
            raise ValueError(
                f"Aircraft '{name}' references unknown marker '{airport}' as its airport"
            )
        radii = settings.get("range_radius")
        if not isinstance(radii, (list, tuple)) or not radii:
            raise ValueError(f"Aircraft '{name}' requires a non-empty range_radius list")
        try:
            radii = tuple(float(radius) for radius in radii)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Aircraft '{name}' range radii must be numeric") from exc
        if any(not np.isfinite(radius) or radius <= 0 for radius in radii):
            raise ValueError(f"Aircraft '{name}' range radii must be positive and finite")
        aircraft.append((name, airport, radii))

    return campaign, markers, aircraft


def _add_density(ax, time_mask, lon, lat, cmap, zorder):
    x = lon[time_mask].ravel()
    y = lat[time_mask].ravel()
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if x.size < 3:
        return
    try:
        density = gaussian_kde(np.vstack((x, y)))
    except np.linalg.LinAlgError:
        return

    # gaussian_kde has non-zero support outside the parcel extrema. Evaluate
    # several kernel bandwidths beyond those extrema so that contourf does not
    # introduce artificial straight edges at x.min/x.max or y.min/y.max.
    longitude_bandwidth, latitude_bandwidth = np.sqrt(
        np.diag(density.covariance)
    )
    longitude_padding = _KDE_PADDING_BANDWIDTHS * longitude_bandwidth
    latitude_padding = _KDE_PADDING_BANDWIDTHS * latitude_bandwidth
    longitude_start = x.min() - longitude_padding
    longitude_stop = x.max() + longitude_padding
    latitude_start = y.min() - latitude_padding
    latitude_stop = y.max() + latitude_padding
    grid_step = complex(0, _KDE_GRID_SIZE)
    xi, yi = np.mgrid[
        longitude_start:longitude_stop:grid_step,
        latitude_start:latitude_stop:grid_step,
    ]
    zi = density(np.vstack((xi.ravel(), yi.ravel())))
    maximum_density = float(np.nanmax(zi))
    if not np.isfinite(maximum_density) or maximum_density <= 0.0:
        return
    minimum_density = _KDE_MINIMUM_RELATIVE_DENSITY * maximum_density
    zi = np.ma.masked_less(zi, minimum_density)
    ax.contourf(
        xi,
        yi,
        zi.reshape(xi.shape),
        levels=np.linspace(minimum_density, maximum_density, 10),
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        zorder=zorder,
    )


def _validated_satellites(satellites):
    if satellites is None:
        return []
    if isinstance(satellites, (str, bytes)) or not isinstance(satellites, Sequence):
        raise TypeError("satellites must be a sequence of satellite names or None")
    names = list(satellites)
    if any(not isinstance(name, str) or not name.strip() for name in names):
        raise TypeError("Every satellites item must be a non-empty string")
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError("Duplicate satellites are not allowed: " + ", ".join(duplicates))
    unknown = sorted(set(names) - set(SATELLITES))
    if unknown:
        raise ValueError(
            "Unknown satellite(s): "
            + ", ".join(unknown)
            + ". Available satellites: "
            + ", ".join(sorted(SATELLITES))
        )
    return names


def _validated_satellite_filter(satellite_filter):
    if satellite_filter not in _SATELLITE_FILTERS:
        raise ValueError(
            "satellite_filter must be 'day_only', 'night_only', or None"
        )
    return satellite_filter


def _inside_geographic_bounds(longitude, latitude, geographic_bounds):
    west, east, south, north = geographic_bounds
    longitude = (longitude + 180.0) % 360.0 - 180.0
    if west < east:
        inside_longitude = (longitude >= west) & (longitude <= east)
    else:
        # A decreasing west/east pair denotes an antimeridian-crossing box.
        inside_longitude = (longitude >= west) | (longitude <= east)
    return (
        np.isfinite(longitude)
        & np.isfinite(latitude)
        & inside_longitude
        & (latitude >= south)
        & (latitude <= north)
    )


def _add_satellite_tracks(
    ax,
    satellites,
    valid_time,
    geographic_bounds,
    satellite_filter,
):
    day_start = valid_time.replace(hour=0, minute=0, second=0, microsecond=0)
    day_stop = day_start + timedelta(days=1)
    fallback_colors = cycle(_FALLBACK_SATELLITE_COLORS)
    line_handles = []

    for satellite in satellites:
        color = _SATELLITE_COLORS.get(satellite)
        if color is None:
            color = next(fallback_colors)
        track = propagate_ground_track(
            satellite,
            day_start,
            day_stop,
            interval=timedelta(minutes=1),
        )
        daylight = classify_daylight(track, minimum_solar_elevation=0.0)

        valid_track = (
            track.successful
            & np.isfinite(track.latitude)
            & np.isfinite(track.longitude)
        )
        if satellite_filter == "day_only":
            selected = valid_track & daylight.daytime
        elif satellite_filter == "night_only":
            selected = (
                valid_track
                & np.isfinite(daylight.solar_elevation)
                & ~daylight.daytime
            )
        else:
            selected = valid_track

        plotted_latitude = np.where(selected, track.latitude, np.nan)
        plotted_longitude = np.where(selected, track.longitude, np.nan)
        longitude_jumps = np.abs(np.diff(plotted_longitude)) > 180.0
        plotted_latitude[1:][longitude_jumps] = np.nan
        plotted_longitude[1:][longitude_jumps] = np.nan

        line, = ax.plot(
            plotted_longitude,
            plotted_latitude,
            color=color,
            linewidth=3,
            transform=ccrs.PlateCarree(),
            zorder=80,
            clip_on=True,
            label=satellite,
        )
        line_handles.append(line)

        marker_times = np.asarray(
            [
                current_time.minute % 5 == 0 and current_time.second == 0
                for current_time in track.times
            ],
            dtype=bool,
        )
        marker_mask = selected & marker_times
        marker_longitude = track.longitude[marker_mask]
        marker_latitude = track.latitude[marker_mask]
        times = track.times[marker_mask]

        inside_map = _inside_geographic_bounds(
            marker_longitude,
            marker_latitude,
            geographic_bounds,
        )
        marker_longitude = marker_longitude[inside_map]
        marker_latitude = marker_latitude[inside_map]
        times = times[inside_map]

        ax.scatter(
            marker_longitude,
            marker_latitude,
            s=100,
            marker="o",
            facecolor=color,
            edgecolor="black",
            linewidth=1,
            transform=ccrs.PlateCarree(),
            zorder=300,
            clip_on=True,
        )
        for index, (longitude, latitude, current_time) in enumerate(
            zip(marker_longitude, marker_latitude, times)
        ):
            offset = (7, 7) if index % 2 == 0 else (7, -7)
            vertical_alignment = "bottom" if index % 2 == 0 else "top"
            annotation = ax.annotate(
                current_time.strftime("%H:%M"),
                xy=(longitude, latitude),
                xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
                xytext=offset,
                textcoords="offset points",
                fontsize=20,
                color=color,
                horizontalalignment="left",
                verticalalignment=vertical_alignment,
                zorder=310,
                annotation_clip=True,
                clip_on=True,
            )
            annotation.set_clip_path(ax.patch)

    if line_handles:
        legend = ax.legend(
            handles=line_handles,
            labels=satellites,
            loc="upper right",
            fontsize=20,
            framealpha=1,
            facecolor="white",
        )
        legend.set_zorder(999)


def _make_axes(
    release_time: datetime,
    valid_time: datetime,
    fire_name: str,
    markers,
    aircraft,
    geographic_bounds,
):
    west, east, south, north = geographic_bounds
    crosses_antimeridian = west > east
    if geographic_bounds == DEFAULT_GEOGRAPHIC_BOUNDS:
        central_longitude = -100.0
        central_latitude = 40.0
    else:
        unwrapped_east = east + 360.0 if crosses_antimeridian else east
        central_longitude = (west + unwrapped_east) / 2.0
        if central_longitude > 180.0:
            central_longitude -= 360.0
        central_latitude = (south + north) / 2.0

    projection = ccrs.LambertConformal(
        central_longitude=central_longitude,
        central_latitude=central_latitude,
    )
    figure = plt.figure(figsize=(16, 20))
    grid = GridSpec(2, 1, height_ratios=[3.5, 1], figure=figure)
    figure.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=0.02)

    map_axis = figure.add_subplot(grid[0], projection=projection)
    if crosses_antimeridian:
        extent_crs = ccrs.PlateCarree(central_longitude=180.0)
        shifted_west = west - 180.0
        shifted_east = east + 180.0
        map_axis.set_extent(
            (shifted_west, shifted_east, south, north),
            crs=extent_crs,
        )
    else:
        map_axis.set_extent(geographic_bounds, crs=ccrs.PlateCarree())
    map_axis.coastlines(resolution="50m", zorder=100)
    map_axis.gridlines(
        draw_labels=True, dms=True, x_inline=False, y_inline=False,
        linewidth=2, color="brown"
    )
    map_axis.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=2, zorder=100)
    map_axis.add_feature(
        cfeature.STATES, linestyle="--", edgecolor="black", linewidth=1, zorder=100
    )
    provinces = cfeature.NaturalEarthFeature(
        category="cultural", name="admin_1_states_provinces_lines",
        scale="50m", facecolor="none", edgecolor="black"
    )
    map_axis.add_feature(provinces, linestyle="--", linewidth=1, zorder=10)
    try:
        map_axis.background_img(name="NE", resolution="high")
    except (ValueError, OSError):
        # The optional Natural Earth raster is installation-specific.
        map_axis.add_feature(cfeature.LAND, facecolor="0.9", zorder=0)
        map_axis.add_feature(cfeature.OCEAN, facecolor="lightblue", zorder=0)

    figure.suptitle(
        f"Fire: {fire_name}  Release: {release_time.isoformat(timespec='hours')}  "
        f"Valid: {valid_time.isoformat(timespec='hours')}",
        fontsize=28,
        fontweight="bold",
        y=0.965,
    )

    # Campaign markers are numbered in YAML insertion order, matching the
    # compact marker style used by the original plot.
    locations = [(*coordinates, str(index)) for index, coordinates in enumerate(markers.values(), 1)]
    xs, ys, labels = zip(*locations)
    map_axis.plot(xs, ys, markersize=28, marker="o", color="black",
                  transform=ccrs.PlateCarree(), linestyle="")
    map_axis.plot(xs, ys, markersize=24, marker="o", color="red",
                  transform=ccrs.PlateCarree(), linestyle="")
    for x, y, label in locations:
        map_axis.text(x, y, label, color="black", transform=ccrs.PlateCarree(),
                      ha="center", va="center", size=16, fontweight="bold")

    ring_styles = cycle(_AIRCRAFT_RING_STYLES)
    for (_, airport, radii), (color, width) in zip(aircraft, ring_styles):
        center = markers[airport]
        for radius_km in radii:
            ring = Polygon(Geodesic().circle(center[0], center[1], radius_km * 1000, 80))
            map_axis.add_feature(cfeature.ShapelyFeature(
                [ring], ccrs.PlateCarree(), facecolor="none", edgecolor=color,
                linewidth=width, linestyle="-", zorder=101
            ))

    altitude_axis = figure.add_subplot(grid[1])
    altitude_axis.set_ylim(6, 16)
    altitude_axis.set_ylabel("Altitude [km]", fontsize=24)
    altitude_axis.set_xlabel("Time", fontsize=24)
    altitude_axis.set_title("Parcel Altitude Over Time", fontsize=26)
    altitude_axis.tick_params(axis="both", labelsize=20)
    altitude_axis.tick_params(axis="x", labelrotation=0)
    return figure, map_axis, altitude_axis


def plot_traj(
    Traj,
    ValidTime,
    CampaignFile="inspyre.yaml",
    satellites=None,
    geographic_bounds=DEFAULT_GEOGRAPHIC_BOUNDS,
    satellite_filter="day_only",
):
    """Create one parcel-trajectory density plot.

    Parameters
    ----------
    Traj : sequence of xarray.Dataset
        One dataset per release altitude.  Each dataset must contain 1-D
        ``time`` and 2-D ``lat``, ``lon``, and ``PAlt`` variables.  The 2-D
        variables must share time and parcel dimensions.  Release altitude is
        read from the numeric global attribute ``altitude_km``.
    ValidTime : datetime.datetime
        UTC plot-valid time. Naive values are interpreted as UTC. The release
        time is read from every dataset's required ``Trajectory_start`` global
        attribute, and the values must agree. The density layer uses ValidTime
        +/- one hour, matching the original implementation's sampling window.
    CampaignFile : path-like
        Required YAML campaign resource.  Its ``Markers`` mapping supplies
        named longitude/latitude locations.  Each entry in ``Aircrafts`` must
        reference one of those marker names as ``airport`` and provide a
        ``range_radius`` list in kilometers.
    satellites : sequence of str or None, optional
        Named satellites whose daytime ground tracks should be added for the
        UTC day containing ``ValidTime``. If None, no orbit data are acquired
        or plotted.
    geographic_bounds : sequence of four numbers, optional
        Map extent in ``[west, east, south, north]`` degrees. The default is
        ``(-120, -70, 22.5, 60)``, matching the original plot. A west value
        greater than east denotes a box crossing the antimeridian; for example,
        ``(22, -155, 40, 80)`` runs eastward from 22 E to 155 W.
    satellite_filter : {"day_only", "night_only", None}, optional
        Select which portions of satellite ground tracks to plot. The default
        is ``"day_only"``. Use ``"night_only"`` for nighttime subpoints or
        ``None`` to apply no solar-elevation filter.

    Returns
    -------
    figure, (map_axis, altitude_axis)
        Matplotlib objects.  The caller is responsible for saving or closing
        the figure.
    """
    if not isinstance(Traj, (list, tuple)) or not Traj:
        raise ValueError("Traj must be a non-empty list or tuple of xarray datasets")
    if not all(isinstance(dataset, xr.Dataset) for dataset in Traj):
        raise TypeError("Every item in Traj must be an xarray.Dataset")

    valid_time = _as_utc_naive(ValidTime, "ValidTime")
    satellite_names = _validated_satellites(satellites)
    geographic_bounds = _validated_geographic_bounds(geographic_bounds)
    satellite_filter = _validated_satellite_filter(satellite_filter)
    datasets = _sorted_release_datasets(Traj)
    release_time = _release_time(datasets)
    fire_name = _fire_name(datasets)
    _, markers, aircraft = _load_campaign(CampaignFile)
    figure, map_axis, altitude_axis = _make_axes(
        release_time,
        valid_time,
        fire_name,
        markers,
        aircraft,
        geographic_bounds,
    )

    window_start = np.datetime64(valid_time - timedelta(hours=1), "us")
    window_stop = np.datetime64(valid_time + timedelta(hours=1), "us")
    schemes = cycle(_COLOR_SCHEMES)
    start_plotted = False

    for layer, (dataset, (cmap, line_color)) in enumerate(zip(datasets, schemes)):
        times, lat, lon, altitude = _trajectory_arrays(dataset)
        time_mask = (times >= window_start) & (times <= window_stop)
        if np.any(time_mask):
            _add_density(map_axis, time_mask, lon, lat, cmap, zorder=97 + layer)

        if not start_plotted and lat.shape[0] and lat.shape[1]:
            map_axis.plot(lon[0, 0], lat[0, 0], marker="p", color="black",
                          transform=ccrs.PlateCarree(), zorder=110)
            start_plotted = True

        for parcel_index in range(0, lat.shape[1], 25):
            map_axis.plot(lon[:, parcel_index], lat[:, parcel_index], color=line_color,
                          transform=ccrs.PlateCarree(), zorder=50, alpha=0.5)
            altitude_axis.plot(times, altitude[:, parcel_index], color=line_color)

    if satellite_names:
        _add_satellite_tracks(
            map_axis,
            satellite_names,
            valid_time,
            geographic_bounds,
            satellite_filter,
        )

    return figure, (map_axis, altitude_axis)
