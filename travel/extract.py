#!/usr/bin/env python3

import json
import math
import sys
from datetime import datetime
from pathlib import Path

import rgeocoder


HOME_LAT = 35.94
HOME_LON = -83.96
THRESHOLD_MILES = 100


def parse_latlon(value):
    value = value.replace("\u00b0", "")
    lat, lon = [x.strip() for x in value.split(",", 1)]
    return float(lat), float(lon)


def parse_date(value):
    return datetime.fromisoformat(value.replace("Z", "+00:00")).date()


def haversine_miles(lat1, lon1, lat2, lon2):
    r = 3958.7613
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    return r * (2 * math.atan2(math.sqrt(a), math.sqrt(1 - a)))


def label_for(rg, lat, lon):
    place = rg.nearest(lat, lon)
    if place.name and place.admin1:
        return f"{place.name}, {place.admin1}"
    if place.name and place.cc:
        return f"{place.name}, {place.cc}"
    if place.name:
        return place.name
    return f"{lat:.2f}, {lon:.2f}"


def main():
    path = sys.argv[1]
    with open(path) as f:
        data = json.load(f)

    rg = rgeocoder.ReverseGeocoder()
    visits = []

    for seg in data.get("semanticSegments", []):
        visit = seg.get("visit")
        if not visit:
            continue
        candidate = visit.get("topCandidate", {})
        latlng = candidate.get("placeLocation", {}).get("latLng")
        if not latlng:
            continue

        lat, lon = parse_latlon(latlng)
        dist = haversine_miles(HOME_LAT, HOME_LON, lat, lon)
        if dist < THRESHOLD_MILES:
            continue

        start = parse_date(seg["startTime"])
        end = parse_date(seg["endTime"])

        visits.append(
            {
                "start": start,
                "end": end,
                "label": label_for(rg, lat, lon),
                "distance": round(dist),
            }
        )

    visits.sort(key=lambda x: (x["start"], x["end"], x["label"]))

    trips = []
    current = None
    for visit in visits:
        if current is None:
            current = visit.copy()
            continue
        gap = (visit["start"] - current["end"]).days
        if gap <= 1 and visit["label"] == current["label"]:
            if visit["end"] > current["end"]:
                current["end"] = visit["end"]
            if visit["distance"] > current["distance"]:
                current["distance"] = visit["distance"]
        elif gap <= 1 and visit["distance"] <= current["distance"] + 25:
            if visit["end"] > current["end"]:
                current["end"] = visit["end"]
            if visit["distance"] > current["distance"]:
                current["distance"] = visit["distance"]
        else:
            trips.append(current)
            current = visit.copy()
    if current is not None:
        trips.append(current)

    size_mb = Path(path).stat().st_size / (1024 * 1024)
    segment_count = len(data.get("semanticSegments", []))

    print("# Travel >100 miles from Knoxville")
    print()
    print(
        f"Extracted from Google Timeline ({size_mb:.1f}MB, {segment_count:,} segments). "
        f"Trips identified by visit segments >{THRESHOLD_MILES} miles from home coordinates "
        f"({HOME_LAT:.2f}\u00b0N, {abs(HOME_LON):.2f}\u00b0W) grouped into consecutive days."
    )
    print()
    print("| Start | End | Days | Location | Distance |")
    print("|-------|-----|------|----------|----------|")
    for trip in trips:
        days = (trip["end"] - trip["start"]).days + 1
        print(
            f"| {trip['start']} | {trip['end']} | {days} | {trip['label']} | ~{trip['distance']} mi |"
        )


if __name__ == "__main__":
    main()
