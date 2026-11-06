"""Tiny results registry so one script's fitted numbers feed the next.

    from results import save, load
    save("DG_gap_model", s=-80.3, P=42.8, ...)      -> analysis/results/DG_gap_model.json
    p = load("DG_gap_model")                          -> dict

Every number that another script or a doc paragraph depends on should be
written here by the script that estimates it and read back by the consumer,
never retyped. Files are small JSON, sorted keys, committed to git so the
dependency chain is visible in diffs.
"""
from pathlib import Path
import json
import datetime as _dt
import numpy as np

RESULTS_DIR = Path(__file__).resolve().parent / "results"


def _clean(v):
    if isinstance(v, dict):
        return {str(k): _clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, np.ndarray)):
        return [_clean(x) for x in v]
    if isinstance(v, (np.floating, float)):
        f = float(v)
        return None if np.isnan(f) else round(f, 6)
    if isinstance(v, (np.integer, int, bool)):
        return int(v) if not isinstance(v, bool) else v
    if isinstance(v, (np.bool_,)):
        return bool(v)
    return v


def save(name, **values):
    RESULTS_DIR.mkdir(exist_ok=True)
    path = RESULTS_DIR / f"{name}.json"
    payload = _clean(values)
    payload["_written"] = _dt.date.today().isoformat()
    payload["_source"] = f"analysis/{name.split('__')[0]}.py" if "_source" not in payload else payload["_source"]
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return path


def load(name, required=True):
    path = RESULTS_DIR / f"{name}.json"
    if not path.exists():
        if required:
            raise FileNotFoundError(
                f"{path} missing — run analysis/{name}.py first (results chain)")
        return None
    with open(path) as fh:
        return json.load(fh)
