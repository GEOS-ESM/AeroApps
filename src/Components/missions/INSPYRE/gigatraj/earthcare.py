"""Convert between EarthCARE orbit/frame identifiers and UTC times.

EarthCARE divides each orbit into eight geographically defined frames.  Frame
boundaries occur where the sub-satellite point crosses 22.5 or 67.5 degrees
latitude in a specified direction.  This module finds those crossings by
propagating an EarthCARE TLE with SGP4.

JAXA EarthCARE orbit numbers are two greater than the revolution numbers in
the public NORAD TLE for catalog object 59908.  For example, NORAD revolution
12877 corresponds to JAXA orbit 12879.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from sat_tracks import SATELLITES, TwoLineElements, acquire_tle, propagate_tle


EARTHCARE_NORAD_ID = SATELLITES["EarthCARE"]
JAXA_ORBIT_OFFSET = 2

# Each entry is (latitude in degrees, direction), where direction is +1 for
# ascending and -1 for descending.  The ninth boundary closes frame H and is
# also the start of frame A in the following orbit.
_BOUNDARIES = (
    (-22.5, +1),
    (+22.5, +1),
    (+67.5, +1),
    (+67.5, -1),
    (+22.5, -1),
    (-22.5, -1),
    (-67.5, -1),
    (-67.5, +1),
    (-22.5, +1),
)
_FRAMES = "ABCDEFGH"
_SEARCH_STEP = timedelta(seconds=30)
_ROOT_TOLERANCE = timedelta(milliseconds=50)
_MAX_TLE_DISTANCE = timedelta(days=14)


def _as_utc(value: datetime, name: str = "time") -> datetime:
    if not isinstance(value, datetime):
        raise TypeError(f"{name} must be a datetime")
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _load_tle(tle_directory: str | Path) -> TwoLineElements:
    """Load the cached EarthCARE TLE, downloading it if it is absent."""
    return acquire_tle(
        EARTHCARE_NORAD_ID,
        cache_directory=tle_directory,
        maximum_cache_age=None,
    )


def _satellite(tle: TwoLineElements):
    # Kept local so callers need not know about sgp4's Satrec representation.
    from sgp4.api import Satrec

    return Satrec.twoline2rv(tle.line1, tle.line2)


def _orbital_period(tle: TwoLineElements) -> timedelta:
    satellite = _satellite(tle)
    return timedelta(minutes=2.0 * np.pi / satellite.no_kozai)


def _latitude(tle: TwoLineElements, when: datetime) -> float:
    track = propagate_tle(
        tle,
        [when],
        satellite_name="EarthCARE",
        maximum_epoch_distance=None,
    )
    if track.sgp4_errors[0] != 0 or not np.isfinite(track.latitude[0]):
        raise RuntimeError(f"SGP4 failed at {when.isoformat()}")
    return float(track.latitude[0])


def _crosses(y0: float, y1: float, direction: int) -> bool:
    if direction > 0:
        return y0 <= 0.0 < y1
    return y0 >= 0.0 > y1


def _refine_crossing(
    tle: TwoLineElements,
    left: datetime,
    right: datetime,
    target_latitude: float,
    direction: int,
) -> datetime:
    """Bisect a bracketed latitude crossing to sub-second precision."""
    while right - left > _ROOT_TOLERANCE:
        middle = left + (right - left) / 2
        residual = _latitude(tle, middle) - target_latitude
        if (direction > 0 and residual > 0.0) or (
            direction < 0 and residual < 0.0
        ):
            right = middle
        else:
            left = middle
    return left + (right - left) / 2


def _find_next_crossing(
    tle: TwoLineElements,
    start: datetime,
    target_latitude: float,
    direction: int,
    *,
    search_duration: timedelta,
) -> datetime:
    """Return the first requested crossing at or after *start*."""
    left = start
    y_left = _latitude(tle, left) - target_latitude
    stop = start + search_duration

    while left < stop:
        right = min(left + _SEARCH_STEP, stop)
        y_right = _latitude(tle, right) - target_latitude
        if _crosses(y_left, y_right, direction):
            return _refine_crossing(
                tle, left, right, target_latitude, direction
            )
        left, y_left = right, y_right

    motion = "ascending" if direction > 0 else "descending"
    raise RuntimeError(
        f"Could not find the {target_latitude:g}-degree {motion} crossing "
        f"within {search_duration} of {start.isoformat()}"
    )


def _orbit_boundaries(
    tle: TwoLineElements,
    orbit: int,
) -> tuple[datetime, ...]:
    """Return the nine ordered boundaries enclosing a JAXA orbit."""
    try:
        orbit = int(orbit)
    except (TypeError, ValueError) as exc:
        raise TypeError("orbit must be an integer") from exc
    if orbit <= 0:
        raise ValueError("orbit must be positive")

    satellite = _satellite(tle)
    tle_orbit = orbit - JAXA_ORBIT_OFFSET
    delta_orbits = tle_orbit - satellite.revnum
    estimate = tle.epoch + delta_orbits * _orbital_period(tle)
    if abs(estimate - tle.epoch) > _MAX_TLE_DISTANCE:
        raise ValueError(
            f"Orbit {orbit} is more than {_MAX_TLE_DISTANCE.days} days from "
            f"the TLE epoch {tle.epoch.isoformat()}; use a closer TLE"
        )

    # The TLE revolution epoch is close to the ascending equator crossing.
    # Frame A begins several minutes earlier at 22.5 S ascending.
    first = _find_next_crossing(
        tle,
        estimate - timedelta(minutes=20),
        *_BOUNDARIES[0],
        search_duration=timedelta(minutes=40),
    )
    boundaries = [first]
    for target_latitude, direction in _BOUNDARIES[1:]:
        boundaries.append(
            _find_next_crossing(
                tle,
                boundaries[-1] + _ROOT_TOLERANCE,
                target_latitude,
                direction,
                search_duration=timedelta(minutes=30),
            )
        )
    return tuple(boundaries)


def get_frametimes(
    orbit: int,
    frame: str,
    tle_directory: str | Path,
) -> tuple[datetime, datetime]:
    """Return the UTC beginning and ending time of an EarthCARE frame.

    Parameters
    ----------
    orbit
        Five-digit JAXA EarthCARE orbit number (leading zeroes are optional).
    frame
        Frame letter A through H; matching is case-insensitive.
    tle_directory
        Directory containing ``59908.tle``.  The current TLE is downloaded
        with :func:`sat_tracks.acquire_tle` if that file is absent.

    Returns
    -------
    tuple of datetime
        Timezone-aware UTC datetimes calculated from the propagated orbit.
    """
    if not isinstance(frame, str):
        raise TypeError("frame must be a string")
    frame = frame.strip().upper()
    if frame not in _FRAMES:
        raise ValueError("frame must be one of A, B, C, D, E, F, G, or H")

    boundaries = _orbit_boundaries(_load_tle(tle_directory), orbit)
    index = _FRAMES.index(frame)
    return boundaries[index], boundaries[index + 1]


def get_frame(
    time: datetime,
    tle_directory: str | Path,
) -> tuple[int, str]:
    """Return the JAXA EarthCARE orbit number and frame containing *time*.

    Naive datetimes are interpreted as UTC.  Frame intervals are half-open:
    a timestamp exactly on a boundary belongs to the frame that begins there.
    """
    when = _as_utc(time)
    tle = _load_tle(tle_directory)
    if abs(when - tle.epoch) > _MAX_TLE_DISTANCE:
        raise ValueError(
            f"time is more than {_MAX_TLE_DISTANCE.days} days from the TLE "
            f"epoch {tle.epoch.isoformat()}; use a closer TLE"
        )

    satellite = _satellite(tle)
    elapsed_orbits = (when - tle.epoch) / _orbital_period(tle)
    estimated_orbit = (
        satellite.revnum + JAXA_ORBIT_OFFSET + int(np.floor(elapsed_orbits))
    )

    # The estimate is referenced to the ascending node whereas JAXA frames
    # start at 22.5 S, so inspect adjacent cycles at the boundary.
    for orbit in range(estimated_orbit - 1, estimated_orbit + 2):
        boundaries = _orbit_boundaries(tle, orbit)
        if boundaries[0] <= when < boundaries[-1]:
            index = next(
                i
                for i in range(len(_FRAMES))
                if boundaries[i] <= when < boundaries[i + 1]
            )
            return orbit, _FRAMES[index]

    raise RuntimeError(f"Could not associate {when.isoformat()} with a frame")


__all__ = ["get_frame", "get_frametimes"]
