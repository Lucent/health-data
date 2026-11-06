#!/usr/bin/env python3
"""Extract rides from GPS Cycle Computer .gcc tracks (Windows Mobile app, Jul–Aug 2010).

The raw .gcc files live here (moved from OneDrive Backup/GPS Logs, now Backup/Google/GPS Logs, 2026-10-16).
Format (little-endian, reverse-engineered; distance matches the app's own
KML summaries exactly: 7.43 and 5.82 mi):
  "GCC" 0x01 | 6 bytes YY MM DD hh mm ss (local start) | double lat0 | double lon0
  then int16 groups of 5:
    (a, b, type, -1, -1)  block: type 2 header (1,1,2), type 1 battery % in a,
                          type 3 note of a UTF-16 chars following the block
    (x, y, alt, ?, t)     point: metres east/north of lat0/lon0, altitude m,
                          unknown, seconds since start (uint16)

Output: rides.csv — date, start_time, elapsed_min, distance_mi, max_alt_m, min_alt_m, note, file
"""

import csv
import math
import struct
from datetime import datetime
from pathlib import Path

DIR = Path(__file__).parent


def parse(path):
    b = path.read_bytes()
    assert b[:4] == b"GCC\x01", path
    start = datetime(2000 + b[4], *b[5:10])
    s = struct.unpack(f"<{(len(b) - 26) // 2}h", b[26:26 + (len(b) - 26) // 2 * 2])
    assert s[:5] == (1, 1, 2, -1, -1), path
    points, notes, i = [], [], 5
    while i + 5 <= len(s):
        a, _, typ, m1, m2 = s[i:i + 5]
        if m1 == -1 and m2 == -1:
            if typ == 1:
                i += 5
            elif typ == 3:
                notes.append("".join(chr(c) for c in s[i + 5:i + 5 + a]))
                i += 5 + a
            else:
                raise ValueError(f"{path.name}: unknown block type {typ} at short {i}")
            continue
        x, y, alt, _, t = s[i:i + 5]
        points.append((t & 0xFFFF, x, y, alt))
        i += 5
    return start, points, notes


def main():
    rows = []
    for path in sorted(DIR.glob("*.gcc")):
        start, pts, notes = parse(path)
        dist = sum(math.hypot(q[1] - p[1], q[2] - p[2]) for p, q in zip(pts, pts[1:]))
        rows.append({
            "date": start.strftime("%Y-%m-%d"),
            "start_time": start.strftime("%H:%M:%S"),
            "elapsed_min": round((pts[-1][0] - pts[0][0]) / 60, 1),
            "distance_mi": round(dist / 1609.344, 2),
            "max_alt_m": max(p[3] for p in pts),
            "min_alt_m": min(p[3] for p in pts),
            "note": "; ".join(notes),
            "file": path.name,
        })
    out_path = DIR / "rides.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} rides → {out_path}")


if __name__ == "__main__":
    main()
