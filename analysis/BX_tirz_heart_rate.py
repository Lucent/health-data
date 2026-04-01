#!/usr/bin/env python3
"""Check whether tirzepatide altered heart rate in Samsung data.

FINDINGS
- Exercise-side live-data is the correct source for per-workout HR traces.
- Tirzepatide heart-rate checks should separate:
  1. pace-adjusted walking HR during exercise
  2. passive overnight heart-rate proxy from tracker tag 21313
- Both comparisons are observational and can still be confounded by time.

Artifacts
- analysis/BX_tirz_heart_rate_summary.csv
"""

from __future__ import annotations

import csv
import json
import os
import statistics as stats
from collections import defaultdict
from datetime import date, datetime, timedelta
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SAMSUNG_ARCHIVE = Path(
    "/mnt/c/Users/Lucent/OneDrive/Documents/Backup/Samsung/samsunghealth_michael_20260324183778.7z"
)
TMP_ROOT = Path("/tmp")


def parse_ts(value: str) -> datetime:
    return datetime.strptime(value, "%Y-%m-%d %H:%M:%S.%f")


def extract_needed_files() -> Path:
    import subprocess
    import tempfile

    tmpdir = Path(tempfile.mkdtemp(prefix="bx_tirz_hr_", dir=str(TMP_ROOT)))
    subprocess.run(
        [
            "7z",
            "x",
            "-y",
            f"-o{tmpdir}",
            str(SAMSUNG_ARCHIVE),
            "samsunghealth_michael_20260324183778/com.samsung.shealth.exercise.20260324183778.csv",
            "samsunghealth_michael_20260324183778/jsons/com.samsung.shealth.exercise/*/*.com.samsung.health.exercise.live_data.json",
            "samsunghealth_michael_20260324183778/com.samsung.shealth.tracker.heart_rate.20260324183778.csv",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return tmpdir / "samsunghealth_michael_20260324183778"


def load_tirz_map() -> dict[date, float]:
    out: dict[date, float] = {}
    with (ROOT / "drugs" / "tirzepatide.csv").open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out[date.fromisoformat(row["date"])] = float(row.get("effective_level") or 0.0)
    return out


def fit_line(xs: list[float], ys: list[float]) -> tuple[float, float]:
    mx, my = stats.mean(xs), stats.mean(ys)
    varx = sum((x - mx) ** 2 for x in xs)
    if varx == 0:
        return my, 0.0
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / varx
    intercept = my - slope * mx
    return intercept, slope


def load_walks(base: Path, tirz: dict[date, float]) -> list[dict[str, float | int | date | datetime]]:
    live_map: dict[str, Path] = {}
    for root, _, files in os.walk(base / "jsons" / "com.samsung.shealth.exercise"):
        for fn in files:
            suffix = ".com.samsung.health.exercise.live_data.json"
            if fn.endswith(suffix):
                live_map[fn[: -len(suffix)]] = Path(root) / fn

    rows: list[dict[str, float | int | date | datetime]] = []
    with (base / "com.samsung.shealth.exercise.20260324183778.csv").open(
        newline="", encoding="utf-8-sig"
    ) as fh:
        next(fh)
        reader = csv.DictReader(fh)
        for row in reader:
            if row["com.samsung.health.exercise.exercise_type"] != "1001":
                continue
            uuid = row["com.samsung.health.exercise.datauuid"]
            live = live_map.get(uuid)
            if live is None:
                continue
            arr = json.loads(live.read_text(encoding="utf-8"))
            hrs = [
                float(obj["heart_rate"])
                for obj in arr
                if "heart_rate" in obj and 35 <= float(obj["heart_rate"]) <= 220
            ]
            if len(hrs) < 10:
                continue
            dist_m = float(row.get("total_distance") or 0)
            if dist_m <= 0:
                dist_m = sum(float(obj["distance"]) for obj in arr if "distance" in obj)
            dur_min = float(row["com.samsung.health.exercise.duration"] or 0) / 60000.0
            if dur_min < 15 or dist_m <= 0:
                continue
            kph = (dist_m / 1000.0) / (dur_min / 60.0)
            if not (2.5 <= kph <= 7.5):
                continue
            d = parse_ts(row["com.samsung.health.exercise.start_time"]).date()
            rows.append(
                {
                    "date": d,
                    "start_time": parse_ts(row["com.samsung.health.exercise.start_time"]),
                    "kph": kph,
                    "hr_mean": stats.fmean(hrs),
                    "hr_samples": len(hrs),
                    "on_tirz": 1 if tirz.get(d, 0.0) > 0 else 0,
                    "effective_level": tirz.get(d, 0.0),
                }
            )
    rows.sort(key=lambda r: r["start_time"])  # type: ignore[index]
    return rows


def load_nights(base: Path, tirz: dict[date, float]) -> list[dict[str, float | int | date]]:
    rows = []
    with (base / "com.samsung.shealth.tracker.heart_rate.20260324183778.csv").open(
        newline="", encoding="utf-8-sig"
    ) as fh:
        next(fh)
        reader = csv.DictReader(fh)
        tracker = []
        for row in reader:
            if row["tag_id"] != "21313":
                continue
            st = parse_ts(row["com.samsung.health.heart_rate.start_time"])
            hr = float(row["com.samsung.health.heart_rate.heart_rate"] or 0)
            if st.year < 2021 or not (35 <= hr <= 220):
                continue
            tracker.append((st, hr))

    bydate: defaultdict[date, list[float]] = defaultdict(list)
    for st, hr in tracker:
        if st.hour <= 6:
            d = st.date()
        elif st.hour >= 22:
            d = st.date() + timedelta(days=1)
        else:
            continue
        bydate[d].append(hr)

    for d, vals in sorted(bydate.items()):
        if len(vals) < 3:
            continue
        rows.append(
            {
                "date": d,
                "n_points": len(vals),
                "min_hr": min(vals),
                "median_hr": stats.median(vals),
                "on_tirz": 1 if tirz.get(d, 0.0) > 0 else 0,
                "effective_level": tirz.get(d, 0.0),
            }
        )
    return rows


def summarize_group(name: str, values: list[float]) -> dict[str, str]:
    return {
        "metric": name,
        "n": str(len(values)),
        "mean": f"{stats.mean(values):.2f}",
        "median": f"{stats.median(values):.2f}",
    }


def main() -> None:
    base = extract_needed_files()
    tirz = load_tirz_map()
    walks = load_walks(base, tirz)
    nights = load_nights(base, tirz)

    modern_walks = [r for r in walks if int(r["date"].year) >= 2021]  # type: ignore[attr-defined]
    pre_walks = [r for r in modern_walks if int(r["on_tirz"]) == 0]
    inter, slope = fit_line(
        [float(r["kph"]) for r in pre_walks],
        [float(r["hr_mean"]) for r in pre_walks],
    )
    for row in modern_walks:
        row["resid_pre_fit"] = float(row["hr_mean"]) - (inter + slope * float(row["kph"]))

    walk_pre = [float(r["resid_pre_fit"]) for r in modern_walks if int(r["on_tirz"]) == 0]
    walk_on = [float(r["resid_pre_fit"]) for r in modern_walks if int(r["on_tirz"]) == 1]
    night_pre = [float(r["median_hr"]) for r in nights if int(r["on_tirz"]) == 0]
    night_on = [float(r["median_hr"]) for r in nights if int(r["on_tirz"]) == 1]

    summary_rows = [
        summarize_group("walk_resid_pre_tirz", walk_pre),
        summarize_group("walk_resid_on_tirz", walk_on),
        summarize_group("night_median_hr_pre_tirz", night_pre),
        summarize_group("night_median_hr_on_tirz", night_on),
    ]
    summary_rows.append(
        {
            "metric": "walk_on_minus_pre_median_bpm",
            "n": str(min(len(walk_pre), len(walk_on))),
            "mean": "",
            "median": f"{stats.median(walk_on) - stats.median(walk_pre):.2f}",
        }
    )
    summary_rows.append(
        {
            "metric": "night_on_minus_pre_median_bpm",
            "n": str(min(len(night_pre), len(night_on))),
            "mean": "",
            "median": f"{stats.median(night_on) - stats.median(night_pre):.2f}",
        }
    )

    out_path = ROOT / "analysis" / "BX_tirz_heart_rate_summary.csv"
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["metric", "n", "mean", "median"])
        writer.writeheader()
        writer.writerows(summary_rows)

    print("=== BX. Tirzepatide and Heart Rate ===")
    print(f"Pre-tirz walk fit: HR = {inter:.2f} + {slope:.2f} * kph (2021+ walking only)")
    print(
        f"Walk residual median: pre={stats.median(walk_pre):+.2f}, "
        f"on={stats.median(walk_on):+.2f}, delta={stats.median(walk_on) - stats.median(walk_pre):+.2f} bpm"
    )
    print(
        f"Night median HR: pre={stats.median(night_pre):.2f}, "
        f"on={stats.median(night_on):.2f}, delta={stats.median(night_on) - stats.median(night_pre):+.2f} bpm"
    )


if __name__ == "__main__":
    main()
