#!/usr/bin/env python3
"""DI. Can a heterogeneous cohort of LINEAR responders produce SURMOUNT-1's concave
dose-response?

DG: a single linear drug term cannot bend the cohort mean the way the trial does
(5 mg → 71% of the 15 mg loss at one-third the exposure), while an Emax term that
does so predicts a flat within-week sawtooth that this subject's data contradict.
Hypothesis to test: every individual is linear, but the cohort MEAN saturates
because (a) s and P vary across people and (b) intake cannot fall below a floor
(you cannot eat less than ~nothing; GI intolerance forces dose reductions and
discontinuations in the treatment-regimen estimand).

Vectorised cohort simulation (N individuals per arm, frozen SP — DG section 6 shows
the latch barely matters on-drug):

    intake_i(t) = max(I₀_i + s_i·B(t) + P_i·(FM₀_i − FM_i(t)),  floor_i)
    TDEE_i(t)   = MSJ(BW_i)·1.4
    partition   75% fat

Scenarios: H0 homogeneous; H1 s, P log-normal (CV 0.5); H2 H1 + intake floor at
50% of TDEE₀; H3 H1 + floor at 35% (near-starvation tolerated); H4 H2 with s and P
negatively correlated (strong responders defend less). For each scenario the mean
drug coefficient is scaled by a cohort multiplier m chosen so the 15 mg arm matches
−22.5%; the test is whether the 5 and 10 mg arms then land at −16.0 / −21.4 (ratio
5/15 = 0.71) without any Emax.
"""

from pathlib import Path
import sys
import numpy as np
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))
import DG_gap_model as dg  # noqa: E402

N = 400
WEEKS = 72
S_SUBJ, P_SUBJ = -80.1, 43.0        # true-cal subject values (DG)
PUB = {5: -16.0, 10: -21.4, 15: -22.5}
S1 = dict(start_kg=104.8, height_cm=167, age=45, pct_female=0.67, fm_frac=0.40, fat_frac_loss=0.75)


def cohort_arm(dose, m, cv_s, cv_p, floor_frac, rho=0.0, seed=0, n=N):
    rng = np.random.default_rng(seed)
    z1, z2 = rng.standard_normal(n), rng.standard_normal(n)
    z2 = rho * z1 + np.sqrt(1 - rho ** 2) * z2
    sig_s, sig_p = np.sqrt(np.log(1 + cv_s ** 2)), np.sqrt(np.log(1 + cv_p ** 2))
    s = m * S_SUBJ * np.exp(sig_s * z1 - sig_s ** 2 / 2)
    P = P_SUBJ * np.exp(sig_p * z2 - sig_p ** 2 / 2)
    B = dg.trial_exposure(dose, WEEKS, 0)
    bw0 = np.full(n, S1["start_kg"] * dg.KG_TO_LBS) * np.exp(0.15 * rng.standard_normal(n) - 0.15 ** 2 / 2)
    fm = bw0 * S1["fm_frac"]
    lm = bw0 - fm
    fm0 = fm.copy()
    tdee0 = dg.mifflin_tdee(bw0 * dg.LBS_TO_KG, S1["height_cm"], S1["age"], S1["pct_female"], 1.4)
    i0 = tdee0.copy()
    floor = floor_frac * tdee0
    ff = S1["fat_frac_loss"]
    e_per_lb = ff * dg.E_FAT + (1 - ff) * dg.E_LEAN
    for day in range(WEEKS * 7):
        intake = np.maximum(i0 + s * B[day] + P * (fm0 - fm), floor)
        bw = fm + lm
        tdee = dg.mifflin_tdee(bw * dg.LBS_TO_KG, S1["height_cm"], S1["age"], S1["pct_female"], 1.4)
        dW = (intake - tdee) / e_per_lb
        fm, lm = fm + ff * dW, lm + (1 - ff) * dW
    pct = ((fm + lm) - bw0) / bw0 * 100
    return pct.mean(), np.percentile(pct, [25, 75]), (intake <= floor + 1e-9).mean()


def fit_m(cv_s, cv_p, floor_frac, rho):
    return minimize_scalar(lambda m: (cohort_arm(15, m, cv_s, cv_p, floor_frac, rho)[0] - PUB[15]) ** 2,
                           bounds=(0.3, 6.0), method="bounded").x


def main():
    print("=" * 96)
    print("DI. Heterogeneous linear responders vs SURMOUNT-1's concave dose-response")
    print("=" * 96)
    print(f"  N = {N} per arm, 72 wk, frozen SP, subject s = {S_SUBJ} cal/unit, P = {P_SUBJ} cal/lb scaled by m")
    print(f"  Published: 5 mg {PUB[5]}, 10 mg {PUB[10]}, 15 mg {PUB[15]}  → 5/15 ratio {PUB[5] / PUB[15]:.2f}, 10/15 {PUB[10] / PUB[15]:.2f}")
    scenarios = [
        ("H0 homogeneous, no floor", 0.0, 0.0, 0.0, 0.0),
        ("H1 s,P CV 0.5, no floor", 0.5, 0.5, 0.0, 0.0),
        ("H2 H1 + floor 50% TDEE₀", 0.5, 0.5, 0.50, 0.0),
        ("H3 H1 + floor 35% TDEE₀", 0.5, 0.5, 0.35, 0.0),
        ("H4 H2 + corr(s,P) = −0.7", 0.5, 0.5, 0.50, -0.7),
        ("H5 s CV 1.0, P CV 0.5, floor 50%", 1.0, 0.5, 0.50, 0.0),
        ("H6 H5 + floor 65% TDEE₀", 1.0, 0.5, 0.65, 0.0),
    ]
    print(f"\n  {'scenario':<34} {'m':>5} {'5 mg':>7} {'10 mg':>7} {'15 mg':>7} {'5/15':>6} {'10/15':>6} {'floored@15':>11} {'IQR@15':>14}")
    for name, cv_s, cv_p, fl, rho in scenarios:
        m = fit_m(cv_s, cv_p, fl, rho)
        res = {d: cohort_arm(d, m, cv_s, cv_p, fl, rho) for d in (5, 10, 15)}
        a5, a10, a15 = (res[d][0] for d in (5, 10, 15))
        print(f"  {name:<34} {m:>5.2f} {a5:>+7.1f} {a10:>+7.1f} {a15:>+7.1f} {a5 / a15:>6.2f} {a10 / a15:>6.2f} "
              f"{res[15][2] * 100:>10.0f}% {res[15][1][0]:>+6.1f}/{res[15][1][1]:>+6.1f}")
    print("\n  Reading: heterogeneity alone does not bend the mean (each individual is linear in exposure, so the")
    print("  mean is too). A floor does bend it — the strongest responders are clipped at high dose — and the")
    print("  floor needed to reach the published 5/15 ratio tells you how much of the 15 mg arm has to be")
    print("  intake-limited (or dose-reduced / discontinued) for a linear-per-person world to look concave.")
    print("  SURMOUNT-1 published IQR of % loss at 15 mg is wide (roughly −10 to −35), consistent with large s spread.")


if __name__ == "__main__":
    main()
