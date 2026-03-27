#!/usr/bin/env python3
"""Compute daily possible sunlight exposure hours from sleep/wake times and location.

For each day with sleep data:
  1. Determine location (home=Knoxville unless on a trip from travel/trips.md)
  2. Compute sunrise/sunset for that location and date
  3. Awake sunlight = time between wake-up and next sleep, clipped to [sunrise, sunset]

Output: steps-sleep/sunlight.csv
"""

import re
import numpy as np
import pandas as pd
from datetime import datetime, date, timedelta, time as dt_time
from pathlib import Path
from astral import LocationInfo
from astral.sun import sun

ROOT = Path(__file__).resolve().parent.parent

HOME_LAT = 35.94
HOME_LON = -83.96

# Approximate coordinates for trip locations (from trips.md)
TRIP_COORDS = {
    "Nashville": (36.16, -86.78),
    "Houston": (29.76, -95.37),
    "West University Place": (29.72, -95.43),
    "Fruit Cove": (30.11, -81.64),
    "Seabrook Island": (32.58, -80.17),
    "Savannah": (32.08, -81.09),
    "Charleston": (32.78, -79.93),
    "Greenville": (34.85, -82.40),
    "North Atlanta": (33.86, -84.34),
    "Atlanta": (33.75, -84.39),
    "Nassau": (25.05, -77.35),
    "Jordan": (31.95, 35.93),
    "Greensboro": (36.07, -79.79),
    "Johnson City": (36.31, -82.35),
    "Bristol": (36.60, -82.19),
    "Jacksonville": (30.33, -81.66),
    "Jackson": (35.61, -88.81),
    "St. Marys": (30.73, -81.55),
    "Memphis": (35.15, -90.05),
    "Las Vegas": (36.17, -115.14),
    "Birmingham": (33.52, -86.81),
    "Rochester": (44.02, -92.47),
    "New Orleans": (29.95, -90.07),
    "Miami": (25.76, -80.19),
    "Denver": (39.74, -104.98),
    "Outer Banks": (35.56, -75.47),
    "Chattanooga": (35.05, -85.31),
    "East Lansing": (42.74, -84.48),
    "Minneapolis": (44.98, -93.27),
    "Cincinnati": (39.10, -84.51),
    "Berkeley": (37.87, -122.27),
    "Locust Grove": (33.35, -84.11),
}


def parse_trips(path):
    """Parse trips.md table into a list of (start_date, end_date, lat, lon)."""
    trips = []
    text = path.read_text()
    # Match table rows: | date | date | days | location | ...
    for line in text.split("\n"):
        m = re.match(r"\|\s*(\d{4}-\d{2}-\d{2})\s*\|\s*(\d{4}-\d{2}-\d{2})\s*\|\s*\d+\s*\|\s*(.+?)\s*\|", line)
        if not m:
            continue
        start = datetime.strptime(m.group(1), "%Y-%m-%d").date()
        end = datetime.strptime(m.group(2), "%Y-%m-%d").date()
        location = m.group(3).strip()

        # Match location to coordinates
        lat, lon = HOME_LAT, HOME_LON
        for key, (klat, klon) in TRIP_COORDS.items():
            if key.lower() in location.lower():
                lat, lon = klat, klon
                break

        trips.append((start, end, lat, lon, location))

    return trips


def get_location(d, trips):
    """Return (lat, lon) for a given date, checking trips first."""
    for start, end, lat, lon, _ in trips:
        if start <= d <= end:
            return lat, lon
    return HOME_LAT, HOME_LON


def compute_sunrise_sunset(d, lat, lon):
    """Compute sunrise and sunset times (as hours since midnight, local solar)."""
    loc = LocationInfo(latitude=lat, longitude=lon)
    try:
        s = sun(loc.observer, date=d)
        # These are UTC datetime objects
        sunrise_utc = s["sunrise"]
        sunset_utc = s["sunset"]
        return sunrise_utc, sunset_utc
    except ValueError:
        # Polar regions, midnight sun, etc.
        return None, None


def utc_offset_hours(offset_str):
    """Parse 'UTC-0500' to -5.0."""
    sign = 1 if offset_str[3] == "+" else -1
    hours = int(offset_str[4:6])
    minutes = int(offset_str[6:8]) if len(offset_str) >= 8 else 0
    return sign * (hours + minutes / 60)


def main():
    sleep = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])
    trips = parse_trips(ROOT / "travel" / "trips.md")

    rows = []
    for _, row in sleep.iterrows():
        d = row["date"].date() if hasattr(row["date"], "date") else row["date"]

        # Skip garbage rows (naps, very short)
        if row["sleep_hours"] < 2:
            continue

        lat, lon = get_location(d, trips)
        result = compute_sunrise_sunset(d, lat, lon)
        if result[0] is None:
            continue

        sunrise_utc, sunset_utc = result

        # Convert sleep_end (local time) to a datetime
        # sleep.csv date is the wake-up date, sleep_end is the wake time in local
        offset_h = utc_offset_hours(row["time_offset"])

        # Parse wake time
        wake_parts = row["sleep_end"].split(":")
        wake_hour, wake_min = int(wake_parts[0]), int(wake_parts[1])
        wake_local = datetime(d.year, d.month, d.day, wake_hour, wake_min)

        # Parse sleep start for next sleep
        # For sunlight window, we need: wake_time to next_sleep_start, clipped to sunrise-sunset
        # Since we don't know next sleep start, use sunset as the upper bound
        # (conservative: assumes awake until sunset)
        sleep_parts = row["sleep_start"].split(":")
        sleep_hour, sleep_min = int(sleep_parts[0]), int(sleep_parts[1])

        # sleep_start is typically after midnight (e.g. 03:00) meaning the previous evening
        # or could be before midnight. Since date = wake date, sleep_start is the night before.
        # For the NEXT sleep, we'd need the next row. Instead, use a fixed assumption:
        # awake from wake_time until ~3am next day (typical). But for sunlight, sunset caps it.

        # Convert sunrise/sunset to local time
        sunrise_local = sunrise_utc + timedelta(hours=offset_h)
        sunset_local = sunset_utc + timedelta(hours=offset_h)

        sunrise_hour = sunrise_local.hour + sunrise_local.minute / 60
        sunset_hour = sunset_local.hour + sunset_local.minute / 60
        wake_hour_dec = wake_hour + wake_min / 60
        daylight_total = sunset_hour - sunrise_hour

        # Possible sunlight = overlap of [wake, sunset] and [sunrise, sunset]
        window_start = max(wake_hour_dec, sunrise_hour)
        window_end = sunset_hour  # assume awake until sunset (3am sleep >> sunset)
        sunlight_hours = max(0, window_end - window_start)

        # Fraction of available daylight captured
        sunlight_fraction = sunlight_hours / daylight_total if daylight_total > 0 else 0

        rows.append({
            "date": d,
            "wake_time": f"{wake_hour:02d}:{wake_min:02d}",
            "sunrise": f"{sunrise_local.hour:02d}:{sunrise_local.minute:02d}",
            "sunset": f"{sunset_local.hour:02d}:{sunset_local.minute:02d}",
            "daylight_total_hours": round(daylight_total, 2),
            "sunlight_hours": round(sunlight_hours, 2),
            "sunlight_fraction": round(sunlight_fraction, 3),
            "lat": lat,
            "lon": lon,
        })

    out = pd.DataFrame(rows)
    out = out.sort_values("date").drop_duplicates(subset=["date"], keep="last")

    out_path = ROOT / "steps-sleep" / "sunlight.csv"
    out.to_csv(out_path, index=False)

    print(f"Wrote {len(out)} days to {out_path}")
    print(f"Date range: {out['date'].min()} to {out['date'].max()}")
    print(f"\nSunlight hours summary:")
    print(f"  Mean: {out['sunlight_hours'].mean():.1f}")
    print(f"  Std:  {out['sunlight_hours'].std():.1f}")
    print(f"  Min:  {out['sunlight_hours'].min():.1f}")
    print(f"  Max:  {out['sunlight_hours'].max():.1f}")
    print(f"\nBy season (month):")
    out["month"] = pd.to_datetime(out["date"]).dt.month
    for m in range(1, 13):
        sub = out[out["month"] == m]
        if len(sub) > 0:
            print(f"  {m:2d}: mean={sub['sunlight_hours'].mean():.1f}h, "
                  f"daylight={sub['daylight_total_hours'].mean():.1f}h, "
                  f"wake={sub['wake_time'].mode().iloc[0] if len(sub) > 0 else '?'}, "
                  f"n={len(sub)}")

    # Show travel days
    travel = out[((out["lat"] - HOME_LAT).abs() > 0.5) | ((out["lon"] - HOME_LON).abs() > 0.5)]
    if len(travel) > 0:
        print(f"\nTravel days: {len(travel)}")
        print(travel[["date", "wake_time", "sunrise", "sunset", "sunlight_hours", "lat", "lon"]].to_string(index=False))


if __name__ == "__main__":
    main()
