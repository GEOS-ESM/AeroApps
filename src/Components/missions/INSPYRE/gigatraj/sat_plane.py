"""Plot parcel intersections with satellite ground-track curtains.

A satellite curtain is the vertical surface whose footprint is a propagated
ground track. Parcel and satellite times do not need to coincide: parcel
history through ``t2`` is intersected geometrically with the satellite curtain
between ``t1`` and ``t2``.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from itertools import cycle
from pathlib import Path
from typing import Sequence

import cartopy.crs as ccrs
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import gaussian_kde
import xarray as xr

from sat_tracks import SATELLITES, classify_daylight, propagate_ground_track


# Keep these synchronized with the release-height colors in traj_plot.py.
_COLOR_SCHEMES = (
    ("YlOrRd", "maroon"),
    ("Blues", "cornflowerblue"),
    ("YlGn", "lime"),
)
_EARTH_RADIUS_KM = 6371.0088
_ARC_TOLERANCE_RAD = 2.0e-6
_MAP_CENTERS = {
    "atlantic": 0.0,
    "pacific": 180.0,
}


@dataclass(frozen=True)
class CurtainIntersections:
    """Geometric intersections for one release height and satellite curtain."""

    release_altitude_km: float
    satellite_times: np.ndarray
    parcel_times: np.ndarray
    latitude: np.ndarray
    longitude: np.ndarray
    altitude_km: np.ndarray
    parcel_indices: np.ndarray

    def __len__(self) -> int:
        return self.altitude_km.size


def _as_utc(value: datetime, name: str) -> datetime:
    if not isinstance(value, datetime):
        raise TypeError(f"{name} must be a datetime")
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _validate_inputs(trajs, t1, t2, satellites):
    if not isinstance(trajs, (list, tuple)) or not trajs:
        raise ValueError("Trajs must be a non-empty list or tuple")
    if not all(isinstance(dataset, xr.Dataset) for dataset in trajs):
        raise TypeError("Every item in Trajs must be an xarray.Dataset")

    t1 = _as_utc(t1, "t1")
    t2 = _as_utc(t2, "t2")
    if t2 <= t1:
        raise ValueError("t2 must be later than t1")

    if isinstance(satellites, (str, bytes)) or not isinstance(satellites, Sequence):
        raise TypeError("satellites must be a non-empty sequence of names")
    satellites = list(satellites)
    if not satellites:
        raise ValueError("satellites must not be empty")
    if any(not isinstance(name, str) or not name.strip() for name in satellites):
        raise TypeError("Every satellite name must be a non-empty string")
    duplicates = sorted({name for name in satellites if satellites.count(name) > 1})
    if duplicates:
        raise ValueError("Duplicate satellites are not allowed: " + ", ".join(duplicates))
    unknown = sorted(set(satellites) - set(SATELLITES))
    if unknown:
        raise ValueError(
            "Unknown satellite(s): "
            + ", ".join(unknown)
            + ". Available satellites: "
            + ", ".join(sorted(SATELLITES))
        )
    return t1, t2, satellites


def _release_altitude(dataset: xr.Dataset) -> float:
    try:
        altitude = float(dataset.attrs["altitude_km"])
    except KeyError as exc:
        raise ValueError("Every dataset must define attrs['altitude_km']") from exc
    except (TypeError, ValueError) as exc:
        raise ValueError("attrs['altitude_km'] must be numeric") from exc
    if not np.isfinite(altitude):
        raise ValueError("attrs['altitude_km'] must be finite")
    return altitude


def _sorted_datasets(trajs):
    pairs = [(_release_altitude(dataset), dataset) for dataset in trajs]
    altitudes = [altitude for altitude, _ in pairs]
    if len(set(altitudes)) != len(altitudes):
        raise ValueError("Each dataset must have a unique altitude_km")
    return sorted(pairs, key=lambda pair: pair[0])


def _shared_fire_name(trajs) -> str:
    names = [str(dataset.attrs.get("fire", "")).strip() for dataset in trajs]
    if any(not name for name in names):
        raise ValueError("Every dataset must define a non-empty attrs['fire']")
    if len(set(names)) != 1:
        raise ValueError("All trajectory datasets must describe the same fire")
    return names[0]


def _shared_release_time(trajs) -> datetime:
    release_times = []
    for index, dataset in enumerate(trajs):
        value = dataset.attrs.get("Trajectory_start")
        if value is None or not str(value).strip():
            raise ValueError(
                "Every dataset must define attrs['Trajectory_start']; "
                f"missing at Trajs index {index}"
            )
        if isinstance(value, datetime):
            release_time = _as_utc(value, "Trajectory_start")
        else:
            try:
                parsed = datetime.fromisoformat(str(value).strip().replace("Z", "+00:00"))
            except ValueError as exc:
                raise ValueError(
                    f"Invalid Trajectory_start at Trajs index {index}: {value!r}"
                ) from exc
            release_time = _as_utc(parsed, "Trajectory_start")
        release_times.append(release_time)
    if len(set(release_times)) != 1:
        values = ", ".join(sorted(time.isoformat() for time in set(release_times)))
        raise ValueError("All datasets must have the same Trajectory_start; found: " + values)
    return release_times[0]


def _trajectory_arrays(dataset: xr.Dataset):
    missing = [name for name in ("time", "lat", "lon", "PAlt") if name not in dataset]
    if missing:
        raise ValueError("Dataset is missing variables: " + ", ".join(missing))
    time = dataset["time"]
    if time.ndim != 1:
        raise ValueError("'time' must be one-dimensional")
    time_dimension = time.dims[0]
    latitude = dataset["lat"]
    parcel_dimensions = [dim for dim in latitude.dims if dim != time_dimension]
    if latitude.ndim != 2 or len(parcel_dimensions) != 1:
        raise ValueError("'lat' must have time and parcel dimensions")
    parcel_dimension = parcel_dimensions[0]
    required_dimensions = {time_dimension, parcel_dimension}

    arrays = []
    for name in ("lat", "lon", "PAlt"):
        variable = dataset[name]
        if variable.ndim != 2 or set(variable.dims) != required_dimensions:
            raise ValueError("'lat', 'lon', and 'PAlt' must use the same two dimensions")
        arrays.append(
            np.asarray(variable.transpose(time_dimension, parcel_dimension).values, dtype=float)
        )
    times = np.asarray(time.values).astype("datetime64[us]")
    if times.size < 2:
        raise ValueError("Each trajectory dataset requires at least two time samples")
    if np.any(np.diff(times).astype("timedelta64[us]") <= np.timedelta64(0, "us")):
        raise ValueError("Trajectory times must be strictly increasing")
    return times, arrays[0], arrays[1], arrays[2]


def _unit_vectors(longitude, latitude):
    longitude = np.deg2rad(longitude)
    latitude = np.deg2rad(latitude)
    cosine_latitude = np.cos(latitude)
    return np.stack(
        (
            cosine_latitude * np.cos(longitude),
            cosine_latitude * np.sin(longitude),
            np.sin(latitude),
        ),
        axis=-1,
    )


def _angular_distance(first, second):
    cross = np.linalg.norm(np.cross(first, second), axis=-1)
    dot = np.sum(first * second, axis=-1)
    return np.arctan2(cross, dot)


def _normalized(vector):
    magnitude = np.linalg.norm(vector, axis=-1, keepdims=True)
    return np.divide(
        vector,
        magnitude,
        out=np.full_like(vector, np.nan, dtype=float),
        where=magnitude > 1.0e-14,
    )


def _on_minor_arc(point, start, stop, arc_length):
    residual = (
        _angular_distance(start, point)
        + _angular_distance(point, stop)
        - arc_length
    )
    return abs(float(residual)) <= _ARC_TOLERANCE_RAD


def _great_circle_intersection(parcel_start, parcel_stop, sat_start, sat_stop):
    parcel_length = float(_angular_distance(parcel_start, parcel_stop))
    satellite_length = float(_angular_distance(sat_start, sat_stop))
    if parcel_length < 1.0e-12 or satellite_length < 1.0e-12:
        return None

    parcel_normal = np.cross(parcel_start, parcel_stop)
    satellite_normal = np.cross(sat_start, sat_stop)
    crossing = np.cross(parcel_normal, satellite_normal)
    magnitude = np.linalg.norm(crossing)
    if magnitude < 1.0e-12:
        # Coincident great circles do not define one unique intersection.
        return None

    crossing /= magnitude
    for point in (crossing, -crossing):
        if _on_minor_arc(point, parcel_start, parcel_stop, parcel_length) and _on_minor_arc(
            point, sat_start, sat_stop, satellite_length
        ):
            parcel_fraction = float(_angular_distance(parcel_start, point) / parcel_length)
            satellite_fraction = float(
                _angular_distance(sat_start, point) / satellite_length
            )
            return point, parcel_fraction, satellite_fraction
    return None


def _satellite_segment_index(track):
    valid = (
        track.successful[:-1]
        & track.successful[1:]
        & np.isfinite(track.latitude[:-1])
        & np.isfinite(track.latitude[1:])
        & np.isfinite(track.longitude[:-1])
        & np.isfinite(track.longitude[1:])
    )
    indices = np.flatnonzero(valid)
    if indices.size == 0:
        raise ValueError(f"Satellite {track.satellite} has no valid ground-track segments")

    points = _unit_vectors(track.longitude, track.latitude)
    starts = points[indices]
    stops = points[indices + 1]
    lengths = _angular_distance(starts, stops)
    usable = np.isfinite(lengths) & (lengths > 1.0e-12) & (lengths < np.pi)
    indices = indices[usable]
    starts = starts[usable]
    stops = stops[usable]
    lengths = lengths[usable]
    midpoints = _normalized(starts + stops)
    return {
        "indices": indices,
        "starts": starts,
        "stops": stops,
        "lengths": lengths,
        "tree": cKDTree(midpoints),
        "maximum_half_length": float(np.max(lengths) / 2.0),
    }


def _interpolate_datetime(start, stop, fraction):
    return start + (stop - start) * float(fraction)


def _vector_to_lon_lat(vector):
    longitude = np.rad2deg(np.arctan2(vector[1], vector[0]))
    latitude = np.rad2deg(np.arctan2(vector[2], np.hypot(vector[0], vector[1])))
    return float(longitude), float(latitude)


def find_curtain_intersections(dataset, track, t2) -> CurtainIntersections:
    """Find intersections between one trajectory dataset and one curtain."""
    release_altitude = _release_altitude(dataset)
    times, latitude, longitude, altitude = _trajectory_arrays(dataset)
    t2_naive = np.datetime64(_as_utc(t2, "t2").replace(tzinfo=None), "us")
    satellite_index = _satellite_segment_index(track)
    satellite_times = track.times

    output_satellite_times = []
    output_parcel_times = []
    output_latitude = []
    output_longitude = []
    output_altitude = []
    output_parcel_indices = []
    seen = set()

    for time_index in range(times.size - 1):
        if times[time_index] >= t2_naive:
            break
        interval_fraction = 1.0
        if times[time_index + 1] > t2_naive:
            interval_fraction = float(
                (t2_naive - times[time_index]) / (times[time_index + 1] - times[time_index])
            )

        for parcel_index in range(latitude.shape[1]):
            coordinates = (
                latitude[time_index, parcel_index],
                longitude[time_index, parcel_index],
                latitude[time_index + 1, parcel_index],
                longitude[time_index + 1, parcel_index],
            )
            if not np.all(np.isfinite(coordinates)):
                continue
            parcel_start = _unit_vectors(coordinates[1], coordinates[0])
            original_stop = _unit_vectors(coordinates[3], coordinates[2])
            parcel_length = float(_angular_distance(parcel_start, original_stop))
            if parcel_length < 1.0e-12 or parcel_length >= np.pi:
                continue

            if interval_fraction < 1.0:
                normal = _normalized(
                    np.sin((1.0 - interval_fraction) * parcel_length) * parcel_start
                    + np.sin(interval_fraction * parcel_length) * original_stop
                )
                parcel_stop = normal
                segment_altitude_stop = altitude[time_index, parcel_index] + interval_fraction * (
                    altitude[time_index + 1, parcel_index] - altitude[time_index, parcel_index]
                )
                segment_time_stop = t2_naive
            else:
                parcel_stop = original_stop
                segment_altitude_stop = altitude[time_index + 1, parcel_index]
                segment_time_stop = times[time_index + 1]

            parcel_length = float(_angular_distance(parcel_start, parcel_stop))
            midpoint = _normalized(parcel_start + parcel_stop)
            search_angle = parcel_length / 2.0 + satellite_index["maximum_half_length"]
            search_radius = 2.0 * np.sin(min(search_angle, np.pi) / 2.0) + 1.0e-10
            candidates = satellite_index["tree"].query_ball_point(midpoint, search_radius)

            for candidate in candidates:
                crossing = _great_circle_intersection(
                    parcel_start,
                    parcel_stop,
                    satellite_index["starts"][candidate],
                    satellite_index["stops"][candidate],
                )
                if crossing is None:
                    continue
                point, parcel_fraction, satellite_fraction = crossing
                satellite_segment = satellite_index["indices"][candidate]
                satellite_time = _interpolate_datetime(
                    satellite_times[satellite_segment],
                    satellite_times[satellite_segment + 1],
                    satellite_fraction,
                )
                parcel_time_us = times[time_index].astype(np.int64) + parcel_fraction * (
                    segment_time_stop.astype(np.int64) - times[time_index].astype(np.int64)
                )
                parcel_time = np.datetime64(int(round(parcel_time_us)), "us")
                crossing_altitude = altitude[time_index, parcel_index] + parcel_fraction * (
                    segment_altitude_stop - altitude[time_index, parcel_index]
                )
                crossing_longitude, crossing_latitude = _vector_to_lon_lat(point)

                key = (
                    parcel_index,
                    int(round(satellite_time.timestamp() * 1000.0)),
                    int(round(crossing_latitude * 1.0e5)),
                    int(round(crossing_longitude * 1.0e5)),
                )
                if key in seen:
                    continue
                seen.add(key)
                output_satellite_times.append(satellite_time)
                output_parcel_times.append(parcel_time)
                output_latitude.append(crossing_latitude)
                output_longitude.append(crossing_longitude)
                output_altitude.append(crossing_altitude)
                output_parcel_indices.append(parcel_index)

    order = np.argsort(np.asarray(output_satellite_times, dtype=object))
    return CurtainIntersections(
        release_altitude_km=release_altitude,
        satellite_times=np.asarray(output_satellite_times, dtype=object)[order],
        parcel_times=np.asarray(output_parcel_times, dtype="datetime64[us]")[order],
        latitude=np.asarray(output_latitude, dtype=float)[order],
        longitude=np.asarray(output_longitude, dtype=float)[order],
        altitude_km=np.asarray(output_altitude, dtype=float)[order],
        parcel_indices=np.asarray(output_parcel_indices, dtype=int)[order],
    )


def _night_intervals(times, nighttime):
    if len(times) < 2:
        return []
    edges = [times[0]]
    for first, second in zip(times[:-1], times[1:]):
        edges.append(first + (second - first) / 2)
    edges.append(times[-1])

    intervals = []
    start = None
    for index, is_night in enumerate(nighttime):
        if is_night and start is None:
            start = edges[index]
        if start is not None and (not is_night or index == len(nighttime) - 1):
            stop_index = index + 1 if is_night else index
            intervals.append((start, edges[stop_index]))
            start = None
    return intervals


def _plot_density(axis, intersections, cmap, point_color, t1, t2):
    finite = np.isfinite(intersections.altitude_km)
    if not np.any(finite):
        return
    times = mdates.date2num(list(intersections.satellite_times[finite]))
    altitudes = intersections.altitude_km[finite]
    if times.size < 4 or np.ptp(times) == 0 or np.ptp(altitudes) == 0:
        axis.scatter(times, altitudes, s=12, color=point_color, alpha=0.65, zorder=20)
        return

    try:
        density = gaussian_kde(np.vstack((times, altitudes)))
    except np.linalg.LinAlgError:
        axis.scatter(times, altitudes, s=12, color=point_color, alpha=0.65, zorder=20)
        return

    x_grid = np.linspace(mdates.date2num(t1), mdates.date2num(t2), 240)
    padding = max(0.25, np.ptp(altitudes) * 0.15)
    y_grid = np.linspace(max(0.0, altitudes.min() - padding), altitudes.max() + padding, 140)
    xx, yy = np.meshgrid(x_grid, y_grid)
    zz = density(np.vstack((xx.ravel(), yy.ravel()))).reshape(xx.shape)
    positive = zz[zz > 0]
    if positive.size == 0:
        return
    levels = np.linspace(positive.max() * 0.08, positive.max(), 9)
    axis.contourf(xx, yy, zz, levels=levels, cmap=cmap, alpha=0.82, zorder=10)


def _format_latitude(latitude):
    if not np.isfinite(latitude):
        return ""
    hemisphere = "N" if latitude >= 0 else "S"
    return f"{abs(latitude):.1f}°{hemisphere}"


def _add_latitude_axis(axis, track):
    top_axis = axis.twiny()
    top_axis.set_xlim(axis.get_xlim())
    ticks = axis.get_xticks()
    lower, upper = axis.get_xlim()
    ticks = ticks[(ticks >= lower) & (ticks <= upper)]
    satellite_numbers = mdates.date2num(list(track.times))
    valid = np.isfinite(track.latitude) & np.isfinite(satellite_numbers)
    latitudes = np.interp(ticks, satellite_numbers[valid], track.latitude[valid])
    top_axis.set_xticks(ticks)
    top_axis.set_xticklabels([_format_latitude(value) for value in latitudes])
    top_axis.set_xlabel(f"{track.satellite} ground-track latitude")
    top_axis.tick_params(axis="x", labelsize=9, pad=2)
    return top_axis


def _normalize_map_center(map_center: str) -> tuple[str, float]:
    if not isinstance(map_center, str):
        raise TypeError("map_center must be 'Atlantic' or 'Pacific'")
    name = map_center.strip().lower()
    try:
        central_longitude = _MAP_CENTERS[name]
    except KeyError as exc:
        raise ValueError("map_center must be 'Atlantic' or 'Pacific'") from exc
    return name, central_longitude


def _wrap_longitudes(longitude, central_longitude):
    """Wrap longitudes about the selected map center and break at its seam."""
    western_edge = central_longitude - 180.0
    wrapped = (np.asarray(longitude, dtype=float) - western_edge) % 360.0 + western_edge
    return wrapped


def _add_orbit_inset(axis, track, map_center):
    _, central_longitude = _normalize_map_center(map_center)
    projection = ccrs.PlateCarree(central_longitude=central_longitude)
    inset = axis.inset_axes(
        [0.82, 0.69, 0.15, 0.24],
        projection=projection,
        zorder=50,
    )
    inset.set_global()
    inset.stock_img()
    inset.coastlines(linewidth=0.45, color="0.2")
    inset.gridlines(
        xlocs=np.arange(-180.0, 181.0, 60.0),
        ylocs=np.arange(-90.0, 91.0, 30.0),
        color="0.25",
        linewidth=0.35,
        alpha=0.55,
    )

    longitude = _wrap_longitudes(track.longitude, central_longitude)
    latitude = track.latitude.copy()
    jumps = np.abs(np.diff(longitude)) > 180.0
    longitude[1:][jumps] = np.nan
    latitude[1:][jumps] = np.nan
    inset.plot(
        longitude,
        latitude,
        color="blue",
        linewidth=1.2,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )
    valid = np.flatnonzero(np.isfinite(track.longitude) & np.isfinite(track.latitude))
    if valid.size:
        inset.scatter(
            track.longitude[valid[0]],
            track.latitude[valid[0]],
            s=18,
            color="red",
            transform=ccrs.PlateCarree(),
            zorder=6,
        )
        arrow_indices = valid[
            np.linspace(0, max(0, valid.size - 2), min(6, valid.size), dtype=int)
        ]
        for index in arrow_indices:
            stop_index = min(index + 3, track.times.size - 1)
            if stop_index == index:
                continue
            start_longitude = longitude[index]
            stop_longitude = longitude[stop_index]
            longitude_change = stop_longitude - start_longitude
            if not np.isfinite(longitude_change) or abs(longitude_change) > 180.0:
                continue
            start = inset.projection.transform_point(
                start_longitude,
                track.latitude[index],
                ccrs.PlateCarree(),
            )
            stop = inset.projection.transform_point(
                stop_longitude,
                track.latitude[stop_index],
                ccrs.PlateCarree(),
            )
            if np.all(np.isfinite((*start, *stop))):
                inset.annotate(
                    "",
                    xy=stop,
                    xytext=start,
                    arrowprops={"arrowstyle": "->", "color": "blue", "lw": 1.0},
                    zorder=7,
                )
    inset.set_title(track.satellite, fontsize=8, pad=1)
    return inset


def plot_satplane(
    Trajs,
    t1,
    t2,
    satellites,
    *,
    tle_directory="TLE",
    propagation_interval=timedelta(minutes=1),
    minimum_solar_elevation=0.0,
    map_center="Atlantic",
):
    """Plot trajectory intersection densities for satellite curtains.

    One panel is produced per satellite. Bottom labels show satellite time;
    synchronized top labels show the possibly non-monotonic ground-track
    latitude. Nighttime is shaded but intersections are not filtered by light.

    The orbit inset uses a global Plate Carree projection. ``map_center`` may
    be ``"Atlantic"`` (the backward-compatible default, with a 0-degree
    central meridian) or ``"Pacific"`` (a 180-degree central meridian and an
    Atlantic seam, corresponding to a 0--360-degree longitude orientation).
    """
    t1, t2, satellites = _validate_inputs(Trajs, t1, t2, satellites)
    _normalize_map_center(map_center)
    if not isinstance(propagation_interval, timedelta) or propagation_interval <= timedelta(0):
        raise ValueError("propagation_interval must be a positive timedelta")
    datasets = _sorted_datasets(Trajs)
    fire_name = _shared_fire_name(Trajs)
    _shared_release_time(Trajs)

    figure, axes = plt.subplots(
        len(satellites),
        1,
        figsize=(16, max(5.0, 5.0 * len(satellites))),
        squeeze=False,
        sharex=True,
    )
    axes = list(axes[:, 0])
    color_schemes = list(_COLOR_SCHEMES)

    for axis, satellite in zip(axes, satellites):
        # sat_tracks uses an exclusive stop; extend one interval so t2 itself
        # is available as the final curtain endpoint.
        track = propagate_ground_track(
            satellite,
            t1,
            t2 + propagation_interval,
            interval=propagation_interval,
            tle_directory=Path(tle_directory),
        )
        within_window = np.asarray([time <= t2 for time in track.times])
        track = type(track)(
            satellite=track.satellite,
            norad_id=track.norad_id,
            tle_epoch=track.tle_epoch,
            times=track.times[within_window],
            latitude=track.latitude[within_window],
            longitude=track.longitude[within_window],
            sgp4_errors=track.sgp4_errors[within_window],
        )
        daylight = classify_daylight(track, minimum_solar_elevation)
        valid_solar = np.isfinite(daylight.solar_elevation)
        nighttime = valid_solar & ~daylight.daytime
        for start, stop in _night_intervals(track.times, nighttime):
            axis.axvspan(
                mdates.date2num(start),
                mdates.date2num(stop),
                color="0.85",
                zorder=0,
            )

        for layer, (release_altitude, dataset) in enumerate(datasets):
            cmap, point_color = color_schemes[layer % len(color_schemes)]
            intersections = find_curtain_intersections(dataset, track, t2)
            _plot_density(axis, intersections, cmap, point_color, t1, t2)
            axis.plot([], [], color=point_color, linewidth=3, label=f"{release_altitude:g} km")

        axis.set_xlim(mdates.date2num(t1), mdates.date2num(t2))
        axis.set_ylabel("Parcel altitude [km]")
        axis.grid(True, color="0.3", alpha=0.35, linewidth=0.7)
        legend = axis.legend(loc="upper left", title="Release altitude")
        legend.set_zorder(100)
        locator = mdates.AutoDateLocator(minticks=4, maxticks=10)
        axis.xaxis.set_major_locator(locator)
        axis.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        _add_latitude_axis(axis, track)
        _add_orbit_inset(axis, track, map_center)

    axes[-1].set_xlabel("Satellite ground-track time [UTC]")
    figure.suptitle(
        f"Fire: {fire_name} — trajectory intersections with satellite curtains",
        fontsize=18,
        fontweight="bold",
        y=0.995,
    )
    figure.subplots_adjust(hspace=0.55, top=0.86, bottom=0.09, left=0.08, right=0.97)
    return figure, axes


__all__ = ["CurtainIntersections", "find_curtain_intersections", "plot_satplane"]
