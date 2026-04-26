"""CD. Sanity-check secondary findings against the current Kalman.

Re-extracts load-bearing numbers from secondary findings (X, AT, AE, AF, U, AB, AK)
and prints them alongside the values cited in FINDINGS. Flags material drift (sign flip
or > 0.05 in correlation) so the docs can be updated when filter changes ripple downstream.
"""

import subprocess
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Each entry: (name, script, list of (label, regex, cited_value)).
CHECKS = [
    ("X temperature", "X_temperature_phase.py", [
        ("rising vs stable Δ°F", r"rising[^=]*=\s*([+\-\d.]+)", "+0.30"),
        ("intake → temp partial|FM", r"intake[^|]*\|FM[^=]*=\s*([+\-\d.]+)", "+0.10"),
    ]),
    ("AT temp retest", "AT_temperature_retest.py", [
        ("TDEE → baseline|FM", r"TDEE[^|]*\|FM[^=]*=\s*([+\-\d.]+)", "+0.03"),
    ]),
    ("AE sleep", "AE_sleep_null.py", [
        ("sleep → next-day cal", r"next[- ]day cal[^=]*=\s*([+\-\d.]+)", "-0.01"),
        ("30d bedtime → TDEE/RMR", r"30[- ]day bedtime[^=]*=\s*([+\-\d.]+)", "-0.10"),
    ]),
    ("AF intake variance", "AF_intake_variance.py", [
        ("std → fat-gain partial r", r"partial r\s*=\s*([+\-\d.]+)", "-0.20"),
    ]),
    ("U steps", "U_steps_compensation.py", [
        ("7d steps → 7d intake β", r"β\s*=\s*([+\-\d.]+)", "+0.103"),
    ]),
    ("AB running vs walking", "AB_running_vs_walking.py", [
        ("2014-2016 same-day cal Δ", r"same[- ]day cal[^=]*=\s*(-?\d+)", "-50"),
    ]),
    ("AK sunlight", "AK_sunlight_exposure.py", [
        ("sunlight → TDEE|FM", r"sunlight[^|]*\|FM[^=]*r\s*=\s*([+\-\d.]+)", "+0.11"),
    ]),
]


def run(script):
    p = subprocess.run(
        ["python", str(ROOT / "analysis" / script)],
        capture_output=True, text=True, timeout=180,
    )
    return p.stdout + p.stderr


def main():
    print(f"{'Finding':<30} {'metric':<35} {'cited':>10} {'now':>10} {'flag':>5}")
    print("-" * 95)
    for name, script, items in CHECKS:
        try:
            out = run(script)
        except Exception as e:
            print(f"  {name:<30} ERROR: {e}")
            continue
        for label, pat, cited in items:
            m = re.search(pat, out, re.IGNORECASE | re.MULTILINE)
            current = m.group(1) if m else "—"
            try:
                cv = float(cited)
                vv = float(current) if current != "—" else None
                if vv is None:
                    flag = "?"
                elif (cv > 0) != (vv > 0) and abs(vv) > 0.02:
                    flag = "FLIP"
                elif abs(vv - cv) > 0.05 and abs(cv) > 0.02:
                    flag = "DRIFT"
                else:
                    flag = "ok"
            except ValueError:
                flag = "?"
            print(f"{name:<30} {label:<35} {cited:>10} {current:>10} {flag:>5}")


if __name__ == "__main__":
    main()
