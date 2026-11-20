"""
Acquire TLEs and propagate satellite ground tracks with SGP4.

This code is a refactoring of satellite_groundtracks.py (origin unknown) for clarity
and to disentangle functionality. Orbit propagation and daylight classification are
intentionally separate.

``get_daytime_ground_track`` remains as a compatibility adapter for existing
plotting code.

Arlindo da Silva, August 2026. Refactoring aided by CODEX.

"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import logging
from pathlib import Path
from typing import Mapping, Sequence

import astropy.units as u
from astropy.coordinates import AltAz, CartesianRepresentation, EarthLocation, ITRS, TEME, get_sun
from astropy.time import Time
import numpy as np
import requests
from sgp4.api import SGP4_ERRORS, Satrec
from sgp4.conveniences import sat_epoch_datetime

from pathlib import Path

LOGGER = logging.getLogger(__name__)
CELESTRAK_GP_URL = "https://celestrak.org/NORAD/elements/gp.php"
DEFAULT_TLE_MAX_AGE = timedelta(hours=24)
DEFAULT_MAX_EPOCH_DISTANCE = timedelta(days=14)

SATELLITES: Mapping[str, int] = {
    "Suomi-NPP": 37849,
    "NOAA-20": 43013,
    "NOAA-21": 54234,
    "Terra": 25994,
    "Aqua": 27424,
    "Sentinel-5P": 42969,
    "EarthCARE": 59908,
    "PACE": 58928,
}


@dataclass(frozen=True)
class TwoLineElements:
    """Validated orbital elements for one NORAD catalog object."""

    norad_id: int
    line1: str
    line2: str
    epoch: datetime
    source: str


@dataclass(frozen=True)
class GroundTrack:
    """Satellite subpoints at UTC timestamps.

    Failed SGP4 samples contain NaN coordinates. ``sgp4_errors`` contains zero
    for successful samples and the documented SGP4 error code otherwise.
    """

    satellite: str
    norad_id: int
    tle_epoch: datetime
    times: np.ndarray
    latitude: np.ndarray
    longitude: np.ndarray
    sgp4_errors: np.ndarray

    @property
    def successful(self) -> np.ndarray:
        return self.sgp4_errors == 0


@dataclass(frozen=True)
class DaylightClassification:
    """Solar elevation and daylight status for a ground track."""

    solar_elevation: np.ndarray
    daytime: np.ndarray
    minimum_solar_elevation: float


def _as_utc(value: datetime, parameter_name: str) -> datetime:
    if not isinstance(value, datetime):
        raise TypeError(f"{parameter_name} must be a datetime")
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _extract_tle_lines(lines: Sequence[str], expected_norad_id: int) -> tuple[str, str]:
    cleaned = [line.strip() for line in lines if line.strip()]
    line1 = next((line for line in cleaned if line.startswith("1 ")), None)
    line2 = next((line for line in cleaned if line.startswith("2 ")), None)
    if line1 is None or line2 is None:
        raise ValueError("TLE content does not contain valid line 1 and line 2 records")

    try:
        line1_id = int(line1[2:7])
        line2_id = int(line2[2:7])
    except ValueError as exc:
        raise ValueError("TLE records contain an invalid NORAD catalog number") from exc
    if line1_id != line2_id:
        raise ValueError("TLE line 1 and line 2 refer to different NORAD objects")
    if line1_id != expected_norad_id:
        raise ValueError(
            f"Requested NORAD {expected_norad_id}, but TLE describes NORAD {line1_id}"
        )
    return line1, line2


def _parse_tle(text: str, norad_id: int, source: str) -> TwoLineElements:
    line1, line2 = _extract_tle_lines(text.splitlines(), norad_id)
    satellite = Satrec.twoline2rv(line1, line2)
    epoch = sat_epoch_datetime(satellite).astimezone(timezone.utc)
    return TwoLineElements(norad_id, line1, line2, epoch, source)


def acquire_tle(
    norad_id: int,
    cache_directory: str | Path = "TLE",
    *,
    force_download: bool = False,
    maximum_cache_age: timedelta | None = DEFAULT_TLE_MAX_AGE,
    timeout: float = 30.0,
    session: requests.Session | None = None,
) -> TwoLineElements:
    """Return a validated TLE, refreshing a missing or stale cache entry.

    Set ``maximum_cache_age=None`` to accept a valid cached TLE regardless of
    file age. The cache is updated atomically after a successful download.
    """
    try:
        norad_id = int(norad_id)
    except (TypeError, ValueError) as exc:
        raise TypeError("norad_id must be an integer") from exc
    if norad_id <= 0:
        raise ValueError("norad_id must be positive")
    if maximum_cache_age is not None and maximum_cache_age < timedelta(0):
        raise ValueError("maximum_cache_age cannot be negative")

    cache_directory = Path(cache_directory)
    cache_file = cache_directory / f"{norad_id}.tle"
    if cache_file.is_file() and not force_download:
        age = datetime.now(timezone.utc) - datetime.fromtimestamp(
            cache_file.stat().st_mtime, tz=timezone.utc
        )
        cache_is_fresh = maximum_cache_age is None or age <= maximum_cache_age
        if cache_is_fresh:
            try:
                return _parse_tle(
                    cache_file.read_text(encoding="utf-8"),
                    norad_id,
                    str(cache_file),
                )
            except (OSError, ValueError) as exc:
                LOGGER.warning("Ignoring invalid cached TLE %s: %s", cache_file, exc)

    client = session or requests
    response = client.get(
        CELESTRAK_GP_URL,
        params={"CATNR": str(norad_id), "FORMAT": "TLE"},
        timeout=timeout,
    )
    response.raise_for_status()
    text = response.text.strip()
    if not text:
        raise RuntimeError(f"CelesTrak returned no TLE for NORAD {norad_id}")
    tle = _parse_tle(text, norad_id, CELESTRAK_GP_URL)

    cache_directory.mkdir(parents=True, exist_ok=True)
    temporary_file = cache_file.with_suffix(".tle.tmp")
    temporary_file.write_text(text + "\n", encoding="utf-8")
    temporary_file.replace(cache_file)
    return tle


def create_time_grid(
    start: datetime,
    stop: datetime,
    interval: timedelta = timedelta(minutes=5),
) -> np.ndarray:
    """Create an inclusive-start, exclusive-stop array of UTC datetimes."""
    start = _as_utc(start, "start")
    stop = _as_utc(stop, "stop")
    if stop <= start:
        raise ValueError("stop must be later than start")
    if interval <= timedelta(0):
        raise ValueError("interval must be positive")

    count = (stop - start + interval - timedelta(microseconds=1)) // interval
    return np.asarray([start + index * interval for index in range(count)], dtype=object)


def _validate_tle_epoch(
    tle: TwoLineElements,
    times: np.ndarray,
    maximum_epoch_distance: timedelta | None,
) -> None:
    if maximum_epoch_distance is None:
        return
    if maximum_epoch_distance < timedelta(0):
        raise ValueError("maximum_epoch_distance cannot be negative")
    farthest = max(abs(current_time - tle.epoch) for current_time in times)
    if farthest > maximum_epoch_distance:
        raise ValueError(
            f"TLE epoch {tle.epoch.isoformat()} is {farthest} from the farthest "
            f"requested time, exceeding the allowed {maximum_epoch_distance}. "
            "Acquire an epoch-appropriate TLE or increase maximum_epoch_distance."
        )


def propagate_tle(
    tle: TwoLineElements,
    times: Sequence[datetime],
    *,
    satellite_name: str = "unknown",
    maximum_epoch_distance: timedelta | None = DEFAULT_MAX_EPOCH_DISTANCE,
) -> GroundTrack:
    """Propagate a validated TLE and return Earth-fixed geodetic subpoints."""
    utc_times = np.asarray([_as_utc(value, "times item") for value in times], dtype=object)
    if utc_times.size == 0:
        raise ValueError("times must not be empty")
    _validate_tle_epoch(tle, utc_times, maximum_epoch_distance)

    astropy_times = Time(list(utc_times), scale="utc")
    satellite = Satrec.twoline2rv(tle.line1, tle.line2)
    errors, positions, _ = satellite.sgp4_array(astropy_times.jd1, astropy_times.jd2)
    errors = np.asarray(errors, dtype=np.int16)

    latitude = np.full(utc_times.size, np.nan, dtype=float)
    longitude = np.full(utc_times.size, np.nan, dtype=float)
    valid = errors == 0
    if np.any(valid):
        teme = TEME(
            CartesianRepresentation(positions[valid].T * u.km),
            obstime=astropy_times[valid],
        )
        earth_fixed = teme.transform_to(ITRS(obstime=astropy_times[valid]))
        location = EarthLocation.from_geocentric(
            earth_fixed.x, earth_fixed.y, earth_fixed.z
        )
        latitude[valid] = location.lat.to_value(u.deg)
        longitude[valid] = location.lon.to_value(u.deg)

    if np.any(~valid):
        descriptions = {
            int(code): SGP4_ERRORS.get(int(code), "unknown error")
            for code in np.unique(errors[~valid])
        }
        LOGGER.warning("SGP4 failures for %s: %s", satellite_name, descriptions)

    return GroundTrack(
        satellite=satellite_name,
        norad_id=tle.norad_id,
        tle_epoch=tle.epoch,
        times=utc_times,
        latitude=latitude,
        longitude=longitude,
        sgp4_errors=errors,
    )


def propagate_ground_track(
    satellite: str,
    start: datetime,
    stop: datetime,
    *,
    interval: timedelta = timedelta(minutes=5),
    tle_directory: str | Path = "TLE",
    force_tle_download: bool = False,
    maximum_cache_age: timedelta | None = DEFAULT_TLE_MAX_AGE,
    maximum_epoch_distance: timedelta | None = DEFAULT_MAX_EPOCH_DISTANCE,
    satellite_catalog: Mapping[str, int] = SATELLITES,
) -> GroundTrack:
    """Acquire orbital elements and propagate one named satellite."""
    if satellite not in satellite_catalog:
        available = ", ".join(sorted(satellite_catalog))
        raise ValueError(f"Unknown satellite {satellite!r}. Available satellites: {available}")
    times = create_time_grid(start, stop, interval)
    tle = acquire_tle(
        satellite_catalog[satellite],
        tle_directory,
        force_download=force_tle_download,
        maximum_cache_age=maximum_cache_age,
    )
    return propagate_tle(
        tle,
        times,
        satellite_name=satellite,
        maximum_epoch_distance=maximum_epoch_distance,
    )


def classify_daylight(
    track: GroundTrack,
    minimum_solar_elevation: float = 0.0,
) -> DaylightClassification:
    """Classify daylight at each satellite subpoint.

    This evaluates illumination at the ground subpoint, not spacecraft eclipse.
    """
    minimum_solar_elevation = float(minimum_solar_elevation)
    solar_elevation = np.full(track.times.size, np.nan, dtype=float)
    valid = track.successful & np.isfinite(track.latitude) & np.isfinite(track.longitude)
    if np.any(valid):
        times = Time(list(track.times[valid]), scale="utc")
        locations = EarthLocation.from_geodetic(
            lon=track.longitude[valid] * u.deg,
            lat=track.latitude[valid] * u.deg,
            height=np.zeros(np.count_nonzero(valid)) * u.m,
        )
        sun = get_sun(times).transform_to(AltAz(obstime=times, location=locations))
        solar_elevation[valid] = sun.alt.to_value(u.deg)
    daytime = np.isfinite(solar_elevation) & (solar_elevation >= minimum_solar_elevation)
    return DaylightClassification(solar_elevation, daytime, minimum_solar_elevation)


def get_daytime_ground_track(
    satellite,
    start,
    stop,
    interval_minutes=5.0,
    min_solar_elevation=0.0,
    force_download=False,
    TLE_dir="./TLE/",
):
    """Compatibility wrapper returning the legacy dictionary structure."""
    track = propagate_ground_track(
        satellite,
        start,
        stop,
        interval=timedelta(minutes=float(interval_minutes)),
        tle_directory=TLE_dir,
        force_tle_download=force_download,
    )
    daylight = classify_daylight(track, min_solar_elevation)
    return {
        "times": track.times,
        "lat": track.latitude,
        "lon": track.longitude,
        "solar_elevation": daylight.solar_elevation,
        "daytime": daylight.daytime,
    }


def write_ground_track_csv(path: str | Path, track: GroundTrack) -> None:
    """Write a propagated ground track to CSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("timestamp_utc", "latitude", "longitude", "sgp4_error"))
        writer.writerows(
            zip(track.times, track.latitude, track.longitude, track.sgp4_errors)
        )


def _parse_datetime(value: str) -> datetime:
    try:
        return _as_utc(datetime.fromisoformat(value.replace("Z", "+00:00")), "datetime")
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid ISO-8601 datetime: {value}") from exc


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Propagate a satellite ground track.")
    parser.add_argument("satellite", choices=sorted(SATELLITES))
    parser.add_argument("--start", required=True, type=_parse_datetime)
    parser.add_argument("--stop", required=True, type=_parse_datetime)
    parser.add_argument("--interval-minutes", type=float, default=5.0)
    parser.add_argument("--tle-directory", type=Path, default=Path("TLE"))
    parser.add_argument("--force-tle-download", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    track = propagate_ground_track(
        args.satellite,
        args.start,
        args.stop,
        interval=timedelta(minutes=args.interval_minutes),
        tle_directory=args.tle_directory,
        force_tle_download=args.force_tle_download,
    )
    write_ground_track_csv(args.output, track)
    LOGGER.info("Wrote %d samples to %s", track.times.size, args.output)
    return 0

def refresh_all_tles(tle_directory: str | Path = "TLE") -> None:
    tle_directory = Path(tle_directory)
    failures = []

    for satellite_name, norad_id in SATELLITES.items():
        try:
            tle = acquire_tle(
                norad_id,
                cache_directory=tle_directory,
                force_download=True,
            )

            print(
                f"Refreshed {satellite_name}: "
                f"NORAD {norad_id}, epoch {tle.epoch.isoformat()}"
            )

        except Exception as error:
            failures.append((satellite_name, error))
            print(f"Failed to refresh {satellite_name}: {error}")

    if failures:
        names = ", ".join(name for name, _ in failures)
        raise RuntimeError(f"Failed to refresh TLEs for: {names}")


if __name__ == "__main__":
    raise SystemExit(main())
