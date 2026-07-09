#!/usr/bin/env python3
"""Detect bath events from KUB water_house hourly meter data.

Reads from the KUB data store (which collects raw utility meter data for its
own purpose) and writes the derived bath-event series alongside this script.
The output is consumed by analyses correlating bath frequency with thermo-
regulation hypotheses (e.g. Lamming Lab's claim that protein intake drives
cold intolerance, which would predict more winter bathing).

Idempotent: re-running overwrites baths.csv with the same content given the
same input.

Bath signature (user-supplied):
  - Fill: any time in the afternoon/evening (14:00-23:00), bathtub-volume
    spike of 7-12 CF. Usually 20-23; occasionally as early as 14:00.
  - Refill: ~3h later, small bump.
  - Short bedtime shower: 2-7h after fill, but ALWAYS before 05:00.
    A bath without a follow-up shower before 5am is rejected.

Calibration:
  - Nov 2025: user was away in California, the resident does not take baths.
    A correct detector returns 0 for Nov 2025.
  - Dec months: user reports taking baths mostly in winter.

Tuning notes:
  - Fill window 7-12 CF was identified by eyeballing the clean Dec-2022/Dec-2024
    detections; 5-6 CF events overlap with dishwasher/laundry-only days.
  - 'Quiet hour before' (≤2 CF) excludes washing-machine/sprinkler runs.
  - Leak/sprinkler filter: no more than 2 hours ≥4 CF in [H..H+7].
  - 2 CF shower floor: midnight-5am water use is common (94% of nights),
    mostly toilet flushes; 2 CF separates a brief real shower (~15 gal) from
    that background. Validated by Nov 2025 false-positive count.
"""
import csv
from datetime import datetime, timedelta
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).parent
KUB_DIR = Path("/mnt/c/Users/Lucent/OneDrive/Documents/Finance/Statements/KUB")
SRC = KUB_DIR / "water_house_hourly.csv"
OUT = HERE / "baths.csv"

FILL_LO, FILL_HI = 7, 12     # CF
PRIOR_MAX = 2                # quiet hour before fill (CF)
REFILL_WINDOW = (2, 3, 4)    # hours after fill (small bump)
SHOWER_MIN_CF = 2            # bedtime shower must be ≥ this (rules out a stray flush)
SHOWER_MAX_CF = 10           # …and not so big it's another fill
SHOWER_CUTOFF_HOUR = 5       # shower must be before 05:00
LEAK_HIGH_HOURS_MAX = 2      # max hours ≥4 CF in [H..H+8]

def load():
    rows = {}
    for r in csv.DictReader(open(SRC)):
        t = datetime.fromisoformat(r["timestamp"])
        rows[t] = float(r["read_value"]) if r["read_value"] not in ("", "None") else 0.0
    return rows

def detect_night(rows, date):
    """Return best match dict for this date's night, or None.
    Identifies the single bedtime shower hour (largest spike in the post-fill
    window, before 05:00), and computes total bath timespan.

    Fills may be single-hour OR split across two consecutive hours (when the
    actual tub-fill happens to straddle a meter-clock minute). About 19% of
    fills split this way.
    """
    best = None
    for H in range(14, 24):
        t_fill = datetime.combine(date, datetime.min.time()).replace(hour=H)
        v0 = rows.get(t_fill, 0)
        v1 = rows.get(t_fill + timedelta(hours=1), 0)

        # Single-hour fill, or two adjacent hours each ≥3 CF summing to fill range.
        if FILL_LO <= v0 <= FILL_HI:
            fill = v0
            fill_span = 1   # hours occupied by the fill
        elif v0 >= 3 and v1 >= 3 and FILL_LO <= v0 + v1 <= FILL_HI:
            fill = v0 + v1
            fill_span = 2
        else:
            continue

        if rows.get(t_fill - timedelta(hours=1), 0) > PRIOR_MAX: continue

        refill = sum(rows.get(t_fill + timedelta(hours=h), 0) for h in REFILL_WINDOW)
        if refill < 1: continue

        # Find the bedtime shower: largest hour from (fill_end + 2) up to 05:00
        # next day. SHOWER_MIN_CF separates a real shower from a stray flush.
        cutoff = (datetime.combine(date, datetime.min.time())
                  .replace(hour=SHOWER_CUTOFF_HOUR) + timedelta(days=1))
        t = t_fill + timedelta(hours=fill_span + 1)
        shower_hr = None
        shower_cf = 0
        while t < cutoff:
            v = rows.get(t, 0)
            if SHOWER_MIN_CF <= v <= SHOWER_MAX_CF and v > shower_cf:
                shower_cf = v
                shower_hr = t
            t += timedelta(hours=1)
        if shower_hr is None: continue

        # Leak filter — at most LEAK_HIGH_HOURS_MAX hours ≥4 CF after the fill ends.
        window = [rows.get(t_fill + timedelta(hours=h), 0) for h in range(fill_span, 9)]
        if sum(1 for v in window if v >= 4) > LEAK_HIGH_HOURS_MAX: continue

        span_hours = (shower_hr - t_fill).total_seconds() / 3600 + 1
        score = fill + refill + shower_cf
        full_pattern = [rows.get(t_fill + timedelta(hours=h), 0) for h in range(0, 9)]
        if best is None or score > best["score"]:
            best = dict(date=date, hour=H, fill=fill, fill_span=fill_span,
                        refill=refill, shower=shower_cf, shower_hr=shower_hr,
                        span_hours=span_hours, score=score, pattern=full_pattern)
    return best

def main():
    rows = load()
    dates = sorted({t.date() for t in rows})
    events = [e for e in (detect_night(rows, d) for d in dates) if e]

    # Write CSV
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date", "fill_hour", "fill_CF", "fill_span_h", "refill_CF",
                    "shower_hour", "shower_CF", "span_hours",
                    "score", "pattern_H..H+8"])
        for e in events:
            w.writerow([e["date"], e["hour"], e["fill"], e["fill_span"], e["refill"],
                        e["shower_hr"].strftime("%Y-%m-%d %H:%M"),
                        e["shower"], e["span_hours"], e["score"],
                        " ".join(str(int(v)) for v in e["pattern"])])

    # Per-month summary
    by_month = defaultdict(int)
    for e in events:
        by_month[(e["date"].year, e["date"].month)] += 1

    print(f"Wrote {OUT.name}  ({len(events)} bath events)")
    print()
    months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    years = sorted({y for y,_ in by_month})
    print(f"{'year':5}", *[f'{m:>3}' for m in months], " total")
    for y in years:
        row = [by_month.get((y, m+1), 0) for m in range(12)]
        print(f"{y:5}", *[f'{v:>3}' if v else '  .' for v in row], f" {sum(row):>5}")
    print()
    nov25 = by_month.get((2025, 11), 0)
    print(f"Nov 2025 (user away in CA, expected 0): {nov25}  "
          f"{'✓' if nov25 == 0 else '✗ FALSE POSITIVE'}")
    print(f"Winter (Dec/Jan/Feb) total: "
          f"{sum(n for (y,m), n in by_month.items() if m in (12,1,2))}")
    print(f"Summer (Jun/Jul/Aug) total: "
          f"{sum(n for (y,m), n in by_month.items() if m in (6,7,8))}")

    # Timespan and shower-volume distributions
    import statistics as stat
    spans = [e["span_hours"] for e in events]
    showers = [e["shower"] for e in events]
    print()
    print(f"Timespan (fill→end-of-shower): min={min(spans):.0f}h  "
          f"median={stat.median(spans):.0f}h  max={max(spans):.0f}h  "
          f"mean={stat.mean(spans):.1f}h")
    print(f"Bedtime shower volume (CF):   min={min(showers):.0f}  "
          f"median={stat.median(showers):.0f}  max={max(showers):.0f}  "
          f"mean={stat.mean(showers):.1f}  "
          f"(1 CF ≈ 7.5 gal; 2-min shower at 2 gpm ≈ 0.5 CF)")
    # Histogram of shower volumes
    from collections import Counter
    h = Counter(int(s) for s in showers)
    print("Bedtime shower volume distribution:")
    for v in sorted(h):
        bar = "█" * (h[v] * 50 // max(h.values()))
        print(f"  {v:>2} CF: {h[v]:>3}  {bar}")

if __name__ == "__main__":
    main()
