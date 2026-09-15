"""Generate GigaTraj launch particles from EarthCARE BS and LR curtains."""

from __future__ import annotations

from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path
import re

import numpy as np
from PIL import Image
from scipy.spatial import cKDTree
import xarray as xr


_ORBIT = re.compile(r"^(\d+)([A-Ha-h])$")
_EARTHCARE_NORAD_ID = 59908


def _parse_orbit(orbit: str) -> tuple[int, str]:
    match = _ORBIT.fullmatch(str(orbit).strip())
    if match is None:
        raise ValueError("orbit must have the form '12938B'")
    return int(match.group(1)), match.group(2).upper()


def _plot_rectangle(rgb: np.ndarray) -> tuple[int, int, int, int]:
    dark = rgb.max(axis=2) < 80
    rows = np.flatnonzero(dark.sum(axis=1) > 0.65 * rgb.shape[1])
    if rows.size < 2:
        raise ValueError("could not locate the horizontal plot axes")
    groups = np.split(rows, np.flatnonzero(np.diff(rows) > 1) + 1)
    top, bottom = (int(round(groups[i].mean())) for i in (0, -1))
    columns = np.flatnonzero(
        dark[top : bottom + 1].sum(axis=0) > 0.65 * (bottom - top + 1)
    )
    if columns.size < 2:
        raise ValueError("could not locate the vertical plot axes")
    groups = np.split(columns, np.flatnonzero(np.diff(columns) > 1) + 1)
    left, right = (int(round(groups[i].mean())) for i in (0, -1))
    if bottom - top < 20 or right - left < 100:
        raise ValueError("detected plotting rectangle is implausible")
    return left, right, top, bottom


def _longest_run(mask: np.ndarray) -> tuple[int, int] | None:
    changes = np.diff(np.pad(mask.astype(np.int8), 1))
    starts, stops = np.flatnonzero(changes == 1), np.flatnonzero(changes == -1)
    if not starts.size:
        return None
    index = int(np.argmax(stops - starts))
    return int(starts[index]), int(stops[index])


def _colorbar(rgb: np.ndarray, bottom: int) -> np.ndarray:
    best = None
    for row in range(bottom + 15, rgb.shape[0]):
        pixels = rgb[row]
        chroma = pixels.max(axis=1).astype(int) - pixels.min(axis=1).astype(int)
        run = _longest_run(chroma > 20)
        # The color bar occupies the lower-left portion of the quick look.
        # Reject long, lightly tinted JPEG-background runs spanning the page.
        plausible = (
            run is not None
            and run[0] < 0.25 * rgb.shape[1]
            and run[1] < 0.55 * rgb.shape[1]
            and 50 < run[1] - run[0] < 0.5 * rgb.shape[1]
        )
        if plausible and (best is None or run[1] - run[0] > best[2] - best[1]):
            best = row, *run
    if best is None or best[2] - best[1] < 50:
        raise ValueError("could not locate the color bar")
    row, start, stop = best
    return rgb[row, start:stop].astype(float)


def _decode(
    pixels: np.ndarray,
    palette: np.ndarray,
    limits: tuple[float, float],
    logarithmic: bool,
) -> np.ndarray:
    flat = pixels.reshape(-1, 3).astype(float)
    distances, indices = cKDTree(palette).query(flat, workers=-1)
    chroma = flat.max(axis=1) - flat.min(axis=1)
    fraction = indices / max(len(palette) - 1, 1)
    low, high = limits
    if logarithmic:
        values = 10 ** (np.log10(low) + fraction * np.log10(high / low))
    else:
        values = low + fraction * (high - low)
    values[(chroma <= 18) | (distances >= 70)] = np.nan
    return values.reshape(pixels.shape[:2])


def _line_groups(indices: np.ndarray) -> list[np.ndarray]:
    if indices.size == 0:
        return []
    return list(np.split(indices, np.flatnonzero(np.diff(indices) > 1) + 1))


def _interpolate_band(values: np.ndarray, axis: int, first: int, last: int) -> None:
    """Bridge one inclusive grid-line band using its nearest clean neighbors."""
    before, after = first - 1, last + 1
    limit = values.shape[axis]
    if before < 0 or after >= limit:
        return
    if axis == 0:
        a, b = values[before].copy(), values[after].copy()
        usable = np.isfinite(a) & np.isfinite(b)
        for index in range(first, last + 1):
            fraction = (index - before) / (after - before)
            row = np.full(values.shape[1], np.nan)
            row[usable] = a[usable] + fraction * (b[usable] - a[usable])
            values[index] = row
    else:
        a, b = values[:, before].copy(), values[:, after].copy()
        usable = np.isfinite(a) & np.isfinite(b)
        for index in range(first, last + 1):
            fraction = (index - before) / (after - before)
            column = np.full(values.shape[0], np.nan)
            column[usable] = a[usable] + fraction * (b[usable] - a[usable])
            values[:, index] = column


def _remove_gridlines(
    values: np.ndarray,
    rgb: np.ndarray,
    box: tuple[int, int, int, int],
) -> np.ndarray:
    """Remove plot-grid artifacts and bridge valid data across their bands."""
    left, right, top, bottom = box
    dark = rgb.max(axis=2) < 80
    row_indices = np.flatnonzero(
        dark[:, left : right + 1].sum(axis=1) > 0.65 * (right - left + 1)
    )
    column_indices = np.flatnonzero(
        dark[top : bottom + 1].sum(axis=0) > 0.65 * (bottom - top + 1)
    )
    repaired = values.copy()
    # Expanding by two pixels removes the dark core and JPEG antialiasing halo.
    for group in _line_groups(row_indices):
        if group[-1] <= top or group[0] >= bottom:
            continue
        first = max(0, int(group[0]) - top - 1 - 2)
        last = min(repaired.shape[0] - 1, int(group[-1]) - top - 1 + 2)
        _interpolate_band(repaired, 0, first, last)
    for group in _line_groups(column_indices):
        if group[-1] <= left or group[0] >= right:
            continue
        first = max(0, int(group[0]) - left - 1 - 2)
        last = min(repaired.shape[1] - 1, int(group[-1]) - left - 1 + 2)
        _interpolate_band(repaired, 1, first, last)
    return repaired


def _arrow_points_left(
    rgb: np.ndarray, left: int, right: int, bottom: int
) -> bool:
    red = (
        (rgb[:, :, 0] > 180)
        & (rgb[:, :, 0] > 1.7 * rgb[:, :, 1])
        & (rgb[:, :, 0] > 1.7 * rgb[:, :, 2])
    )
    band = red[max(0, bottom - 8) : min(len(rgb), bottom + 18)]
    nl = band[:, max(0, left - 25) : left].sum()
    nr = band[:, right + 1 : min(rgb.shape[1], right + 26)].sum()
    if nl == nr == 0:
        raise ValueError("could not determine time direction from the red arrow")
    return bool(nl > nr)


def _frame_track(
    orbit: int, frame: str, count: int, tle_directory: str | Path
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from earthcare import get_frametimes
    from sat_tracks import acquire_tle, propagate_tle

    start, stop = get_frametimes(orbit, frame, tle_directory)
    first = np.datetime64(start.replace(tzinfo=None), "ns").astype(np.int64)
    last = np.datetime64(stop.replace(tzinfo=None), "ns").astype(np.int64)
    nanoseconds = np.linspace(first, last, count).astype(np.int64)
    datetimes = [
        datetime.fromtimestamp(value / 1e9, tz=timezone.utc)
        for value in nanoseconds
    ]
    tle = acquire_tle(
        _EARTHCARE_NORAD_ID,
        cache_directory=tle_directory,
        maximum_cache_age=None,
    )
    track = propagate_tle(
        tle,
        datetimes,
        satellite_name="EarthCARE",
        maximum_epoch_distance=None,
    )
    if not np.all(track.successful):
        raise RuntimeError("EarthCARE ground-track propagation failed")
    return nanoseconds.astype("datetime64[ns]"), track.longitude, track.latitude


def _product_specification(
    filename: str | Path,
) -> tuple[str, str, tuple[float, float], bool]:
    """Return name, description, limits, and scale inferred from a filename."""
    name = Path(filename).name.lower()
    if "backscatter" in name:
        return "bs", "backscatter", (1e-7, 1e-4), True
    if "lidar_ratio" in name:
        return "lr", "lidar ratio", (0.0, 100.0), False
    raise ValueError(
        f"cannot infer the quick-look quantity and color scale from {filename!s}"
    )


def parse_images(
    filenames: Sequence[str | Path],
    orbit: str,
    *,
    tle_directory: str | Path = "TLE",
) -> tuple[xr.DataArray, ...]:
    """Decode compatible quick looks, returning arrays in the supplied order."""
    if isinstance(filenames, (str, Path)) or not isinstance(filenames, Sequence):
        raise TypeError("filenames must be a sequence of image paths")
    if not filenames:
        raise ValueError("filenames cannot be empty")
    orbit_number, frame = _parse_orbit(orbit)
    images = [np.asarray(Image.open(path).convert("RGB")) for path in filenames]
    shapes = {image.shape for image in images}
    if len(shapes) != 1:
        raise ValueError("all images must have identical pixel dimensions")
    boxes = [_plot_rectangle(image) for image in images]
    if any(box != boxes[0] for box in boxes[1:]):
        raise ValueError("image plotting rectangles do not match")
    left, right, top, bottom = boxes[0]
    directions = [
        _arrow_points_left(image, left, right, bottom) for image in images
    ]
    if any(direction != directions[0] for direction in directions[1:]):
        raise ValueError("image time directions do not match")

    decoded = []
    specifications = []
    for filename, image, box in zip(filenames, images, boxes, strict=True):
        specification = _product_specification(filename)
        name, _, limits, logarithmic = specification
        values = _decode(
            image[top + 1 : bottom, left + 1 : right],
            _colorbar(image, bottom),
            limits,
            logarithmic,
        )
        values = _remove_gridlines(values, image, box)[::-1]
        if directions[0]:
            values = values[:, ::-1]
        decoded.append(values)
        specifications.append(specification)

    nalt, ntime = decoded[0].shape
    times, longitude, latitude = _frame_track(
        orbit_number, frame, ntime, tle_directory
    )
    coords = {
        "time": times,
        "altitude": (
            "altitude",
            np.linspace(0.0, 20.0, nalt),
            {"units": "km", "positive": "up"},
        ),
        "longitude": ("time", longitude, {"units": "degrees_east"}),
        "latitude": ("time", latitude, {"units": "degrees_north"}),
    }
    common = {
        "orbit": orbit_number,
        "frame": frame,
        "orbit_frame": f"{orbit_number}{frame}",
        "source_kind": "quick-look image",
    }
    arrays = []
    for filename, values, specification in zip(
        filenames, decoded, specifications, strict=True
    ):
        name, long_name, _, logarithmic = specification
        arrays.append(
            xr.DataArray(
                values,
                dims=("altitude", "time"),
                coords=coords,
                name=name,
                attrs={
                    **common,
                    "long_name": long_name,
                    "source": str(filename),
                    "color_scale": "logarithmic" if logarithmic else "linear",
                },
            )
        )
    return tuple(arrays)


def parse_datafiles(
    filenames: Sequence[str | Path], orbit: str
) -> tuple[xr.DataArray, ...]:
    """Native-product entry point reserved until HDF5 paths are specified."""
    raise NotImplementedError(
        "EarthCARE HDF5 dataset paths and quality flags have not been specified"
    )


def time_bounds(
    lat1: float,
    lat2: float,
    orbit: str,
    *,
    tle_directory: str | Path = "TLE",
) -> tuple[datetime, datetime]:
    """Convert two latitudes on a monotonic orbit frame to UTC time bounds.

    Parameters
    ----------
    lat1, lat2
        Latitude endpoints in degrees north. Their order does not matter.
    orbit
        Compact EarthCARE orbit/frame identifier, for example ``"12938B"``.
    tle_directory
        Directory containing the cached EarthCARE ``59908.tle`` file.

    Returns
    -------
    tuple of datetime
        Timezone-aware UTC datetimes ``(t1, t2)`` with ``t2 > t1``.

    Raises
    ------
    ValueError
        If the frame latitude is not monotonic, either latitude lies outside
        the frame, or the two latitudes map to the same time.
    """
    latitudes_requested = np.asarray([lat1, lat2], dtype=float)
    if not np.isfinite(latitudes_requested).all():
        raise ValueError("lat1 and lat2 must be finite")
    if np.any(np.abs(latitudes_requested) > 90.0):
        raise ValueError("lat1 and lat2 must be between -90 and 90 degrees")

    orbit_number, frame = _parse_orbit(orbit)
    times, _, latitudes = _frame_track(
        orbit_number, frame, 4097, tle_directory
    )
    latitudes = np.asarray(latitudes, dtype=float)
    if not np.isfinite(latitudes).all():
        raise ValueError("ground track contains invalid latitudes")

    differences = np.diff(latitudes)
    increasing = np.all(differences >= 0.0) and np.any(differences > 0.0)
    decreasing = np.all(differences <= 0.0) and np.any(differences < 0.0)
    if not (increasing or decreasing):
        raise ValueError(
            f"latitude is not monotonic within orbit frame {orbit_number}{frame}"
        )

    interpolation_latitudes = latitudes if increasing else latitudes[::-1]
    time_ns = times.astype("datetime64[ns]").astype(np.int64)
    time_origin = int(time_ns[0])
    relative_times = time_ns - time_origin
    interpolation_times = relative_times if increasing else relative_times[::-1]
    minimum, maximum = interpolation_latitudes[[0, -1]]
    if np.any(latitudes_requested < minimum) or np.any(
        latitudes_requested > maximum
    ):
        raise ValueError(
            f"requested latitude is outside frame range [{minimum:g}, {maximum:g}]"
        )

    requested_ns = time_origin + np.rint(np.interp(
        latitudes_requested, interpolation_latitudes, interpolation_times
    )).astype(np.int64)
    requested_ns.sort()
    if requested_ns[0] == requested_ns[1]:
        raise ValueError("lat1 and lat2 must define a nonzero time interval")
    result = tuple(
        datetime.fromtimestamp(value / 1e9, tz=timezone.utc)
        for value in requested_ns
    )
    return result


def _utc(value: datetime, name: str) -> np.datetime64:
    if not isinstance(value, datetime):
        raise TypeError(f"{name} must be a UTC datetime")
    if value.tzinfo is None or value.utcoffset() is None:
        raise ValueError(f"{name} must be timezone-aware UTC")
    if value.utcoffset().total_seconds() != 0:
        raise ValueError(f"{name} must be expressed in UTC")
    return np.datetime64(value.replace(tzinfo=None), "ns")


def _edges(centers: np.ndarray) -> np.ndarray:
    centers = np.asarray(centers)
    if centers.size < 2:
        raise ValueError("at least two grid centers are required")
    if np.issubdtype(centers.dtype, np.datetime64):
        numeric = centers.astype("datetime64[ns]").astype(np.int64)
    else:
        numeric = centers.astype(float)
    edges = np.empty(numeric.size + 1)
    edges[1:-1] = (numeric[:-1] + numeric[1:]) / 2
    edges[0] = numeric[0] - np.diff(numeric[:2])[0] / 2
    edges[-1] = numeric[-1] + np.diff(numeric[-2:])[0] / 2
    return edges


def generate_particles(
    P: xr.DataArray,
    mask: xr.DataArray,
    N: int,
    bounds: tuple[datetime, datetime, float, float],
    *,
    seed: int | None = None,
) -> xr.Dataset:
    """Sample particles using a positive field as an unnormalized PDF proxy.

    Bounds are ``(UTC start, UTC stop, minimum altitude km, maximum altitude
    km)``. ``mask`` must be Boolean and exactly aligned with ``P``. Sampled
    orbital times locate particles on the ground track but are not retained in
    the fixed-location output.
    """
    if isinstance(N, bool) or not isinstance(N, (int, np.integer)) or N <= 0:
        raise ValueError("N must be a positive integer")
    if len(bounds) != 4:
        raise ValueError("bounds must be (UTC start, UTC stop, min km, max km)")
    start, stop = _utc(bounds[0], "bounds[0]"), _utc(bounds[1], "bounds[1]")
    amin, amax = float(bounds[2]), float(bounds[3])
    if stop <= start:
        raise ValueError("the ending UTC bound must follow the starting bound")
    if not np.isfinite([amin, amax]).all() or amax <= amin:
        raise ValueError("altitude bounds must be finite and increasing")

    required = {"time", "altitude", "longitude", "latitude"}
    for name, array in (("P", P), ("mask", mask)):
        missing = required - set(array.coords)
        if missing:
            raise ValueError(f"{name} is missing coordinates: {sorted(missing)}")
        if set(array.dims) != {"altitude", "time"}:
            raise ValueError(f"{name} must have altitude and time dimensions")
    P, mask = xr.align(
        P.transpose("altitude", "time"),
        mask.transpose("altitude", "time"),
        join="exact",
    )
    if not np.issubdtype(mask.dtype, np.bool_):
        raise TypeError("mask must be a Boolean xarray.DataArray")
    times = P.time.values.astype("datetime64[ns]")
    altitudes = P.altitude.values.astype(float)
    selected = ((altitudes >= amin) & (altitudes <= amax))[:, None] & (
        (times >= start) & (times <= stop)
    )[None, :]
    values = np.asarray(P, float)
    valid = selected & np.isfinite(values) & (values > 0)
    valid &= np.asarray(mask, dtype=bool)
    candidates = np.flatnonzero(valid)
    if not candidates.size:
        raise ValueError("no positive, finite P cells satisfy the bounds and mask")

    retained = values.ravel()[candidates]
    mean_P = float(retained.mean())
    relative = retained / mean_P
    probabilities = relative / relative.sum()
    rng = np.random.default_rng(seed)
    chosen = rng.choice(candidates, int(N), replace=True, p=probabilities)
    iz, it = np.unravel_index(chosen, values.shape)
    tedges, zedges = _edges(times), _edges(altitudes)
    sampled_time = rng.uniform(
        np.maximum(tedges[it], start.astype(np.int64)),
        np.minimum(tedges[it + 1], stop.astype(np.int64)),
    )
    sampled_altitude = rng.uniform(
        np.maximum(zedges[iz], amin), np.minimum(zedges[iz + 1], amax)
    )
    grid_time = times.astype(np.int64)
    longitude = np.interp(sampled_time, grid_time, P.longitude.values)
    latitude = np.interp(sampled_time, grid_time, P.latitude.values)

    particles = xr.Dataset(
        {
            "lon": (("time", "id"), longitude[None].astype(np.float32)),
            "lat": (("time", "id"), latitude[None].astype(np.float32)),
            "PAlt": (("time", "id"), sampled_altitude[None].astype(np.float32)),
        },
        coords={"time": [0.0], "id": np.arange(N, dtype=float)},
        attrs={
            "Contents": "gigatraj trajectories",
            "Trajectory_start": "0",
            "source_orbit_frame": P.attrs.get("orbit_frame", "unknown"),
            "probability_proxy": P.name or "unnamed",
            "P_mean": mean_P,
            "random_seed": "random" if seed is None else int(seed),
        },
    )
    particles.time.attrs.update(
        long_name="time", standard_name="time", units="days since 0:00"
    )
    particles.id.attrs.update(long_name="parcel id", units="1")
    particles.lon.attrs.update(long_name="longitude", units="degrees_east")
    particles.lat.attrs.update(long_name="latitude", units="degrees_north")
    particles.PAlt.attrs.update(
        units="km", positive="up", vertical_coordinate="yes"
    )
    for variable in (particles.lon, particles.lat, particles.PAlt):
        variable.encoding["_FillValue"] = np.float32(np.nan)
    particles.encoding["unlimited_dims"] = {"time"}
    return particles


__all__ = [
    "parse_images",
    "parse_datafiles",
    "time_bounds",
    "generate_particles",
]
