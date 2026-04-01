#!/usr/bin/env python3
"""Merge body_temperature.csv + body_temperature2.csv → temperature.csv

Samsung Health exports temperature in two files with different schemas.
Both use the phone's home timezone (Eastern US) for timestamps, NOT local
time during travel. The canonical output preserves the original timestamps
and adds a local_hour column corrected for travel using the sleep CSV's
time_offset (which reflects the phone's actual timezone setting each day).

Travel periods during the temperature era (Dec 2023 - Mar 2026):
  - Minneapolis (Oct 2024, 6d): Central, -1h
  - New Orleans (Mar 2025, 6d): Central, -1h
  - Berkeley, CA (Oct 31 - Dec 1 2025, 32d): Pacific, -3h
  - All others: Eastern (same as home)

Source of truth for timezone: steps-sleep/sleep.csv time_offset column,
which records the phone's actual UTC offset each day.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def parse_utc_offset_hours(offset_str):
    """Convert 'UTC-0500' to -5.0, 'UTC+0300' to +3.0."""
    if pd.isna(offset_str):
        return np.nan
    s = str(offset_str).replace("UTC", "")
    sign = -1 if s[0] == "-" else 1
    s = s.lstrip("+-")
    return sign * (int(s[:2]) + int(s[2:]) / 60)


def home_offset_for_date(date):
    """US Eastern: EDT (UTC-4) Mar-Nov, EST (UTC-5) Nov-Mar.
    Approximate — exact DST transitions vary by a week or two."""
    m = date.month
    return -4.0 if 4 <= m <= 10 else -5.0


def main():
    temp_dir = Path(__file__).resolve().parent

    # Load both source files
    bt1 = pd.read_csv(temp_dir / "body_temperature.csv")
    bt1.columns = ["timestamp", "temp_f"]

    bt2 = pd.read_csv(temp_dir / "body_temperature2.csv")
    bt2.columns = ["timestamp", "temp_f"]

    # Combine and deduplicate
    combined = pd.concat([bt1, bt2], ignore_index=True)
    combined["timestamp"] = pd.to_datetime(combined["timestamp"])
    combined = combined.drop_duplicates(subset=["timestamp"]).sort_values("timestamp").reset_index(drop=True)
    # Day boundary at 5am: readings before 5am belong to the previous calendar day's
    # waking period. Subject always sleeps before 5am and never wakes before 5am.
    combined["date"] = (combined["timestamp"] - pd.Timedelta(hours=5)).dt.floor("D")

    print(f"Source 1: {len(bt1)} readings")
    print(f"Source 2: {len(bt2)} readings")
    print(f"Combined (deduplicated): {len(combined)} readings")
    print(f"Date range: {combined['timestamp'].min()} to {combined['timestamp'].max()}")

    # Load sleep timezone offsets
    sleep = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])
    sleep["local_offset_h"] = sleep["time_offset"].apply(parse_utc_offset_hours)

    # Merge timezone info onto each reading
    combined = combined.merge(sleep[["date", "local_offset_h"]], on="date", how="left")

    # Compute timezone correction
    # Timestamps are in home Eastern time. Local hour = eastern_hour + (local_offset - home_offset)
    combined["home_offset_h"] = combined["timestamp"].apply(
        lambda ts: home_offset_for_date(ts))
    combined["tz_correction_h"] = (
        combined["local_offset_h"].fillna(combined["home_offset_h"]) - combined["home_offset_h"])

    raw_hour = combined["timestamp"].dt.hour + combined["timestamp"].dt.minute / 60
    combined["local_hour"] = (raw_hour + combined["tz_correction_h"]) % 24

    n_corrected = (combined["tz_correction_h"].abs() > 0.1).sum()
    print(f"Timezone-corrected readings: {n_corrected}")

    # Write canonical output
    # Keep original timestamp (Eastern) for reproducibility, add local_hour
    out = combined[["timestamp", "temp_f", "local_hour"]].copy()
    out.columns = ["date", "temp_f", "local_hour"]
    out.to_csv(temp_dir / "temperature.csv", index=False)
    print(f"Wrote {len(out)} readings to temperature.csv")


if __name__ == "__main__":
    main()
