"""CR. Identify the intake multiplier by profile innovation likelihood.

The premise of CI/CJ — that the multiplier is pinned by driving the Kalman's
cumulative eb_residual_lb to zero — is close to vacuous: TDEE is a free random
walk, so scaling intake by m raises TDEE by (m-1)*intake almost day for day.
Measured leakage is dEB/dm = 24.9 lb per unit m against the 3330.7 lb/unit that
would obtain if TDEE were fixed, i.e. the filter absorbs 99.25% of any level
change. CJ's protein sweep contains two configs differing by 8.4% of lifetime
intake that both close EB to within 0.002 lb.

But the level is NOT unidentifiable. The information lives in the weight
innovation sequence, which nothing in the pipeline was scoring. Getting the
scale wrong mis-predicts the weight RESPONSE to intake variation: too small an m
under-predicts the swings, too large over-predicts them. Both inflate the
innovations. That produces a genuine interior maximum in

    LL(m) = -0.5 * sum_t [ log(2*pi*S_t) + innov_t^2 / S_t ]

evaluated over weigh-in days, with S_t the filter's own innovation variance.

This is the same signal CQ_intake_scale_response.py estimates by direct
regression outside the filter. The two are therefore CONSISTENCY checks on each
other, not independent confirmations - both exploit the weight response to
intake variation.

Tuning note: tightening Q_BUF and sigma_TDEE narrows the interval substantially,
but degrades the fit to the composition scans (FM_MAE 1.08 -> 1.49 lb). That is
buying precision by forcing real water/gut variation into the FM state. The
canonical tuning is retained here; the tighter variants are reported so the
tradeoff is visible.
"""
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))
import P4_kalman_filter as P4
from results import save

CHI2_95 = 1.92          # half-width of the 95% profile-likelihood interval, 1 df
GRID = np.round(np.arange(0.94, 1.46, 0.02), 2)

CONFIGS = [
    ("canonical",      P4.Q_BUF, 50, 100),
    ("tight BUF",      0.02,     50, 100),
    ("tight TDEE",     P4.Q_BUF, 10,  20),
    ("tight both",     0.02,     10,  20),
]


def run_one(data, mult, q_buf, sigma_tdee, max_dtdee):
    smooth, intake, strength, comp_scans, rmr_df, p2, temp_df = data
    P4.Q_BUF = q_buf
    dates, cal, wobs, fm_scan, lm_scan, erm, wpw = \
        P4.build_inputs(smooth, intake, strength, comp_scans, p2,
                        calories_true=intake["calories"].values * mult)
    fs = int(np.where(~np.isnan(fm_scan))[0][0])
    coll = {}
    xf, Pf, xp, Pp, innov, _ = P4.kalman_v4_forward(
        cal, wobs, fm_scan, lm_scan, erm, wpw, fm_scan[fs], lm_scan[fs],
        sigma_tdee, max_dtdee, collect=coll)
    S = coll["innov_var"]
    ok = ~np.isnan(innov) & ~np.isnan(S)
    ll = float(-0.5 * np.sum(np.log(2 * np.pi * S[ok]) + innov[ok] ** 2 / S[ok]))
    xs, _ = P4.rts_smoother_v4(xf, Pf, xp, Pp, fs, wpw)
    d = P4.diagnostics(dates, xs, innov, fm_scan, lm_scan, cal, rmr_df, wobs)
    return ll, d["fm_scan_mae"], d["tdee_rmr_mean"], d["tdee_mean"]


def main():
    data = P4.load_data()
    q_buf_0 = P4.Q_BUF

    print("=== CR. Intake multiplier by profile innovation likelihood ===")
    from intake_model import describe
    print(f"Canonical intake model (analysis/intake_model.py): {describe()}")
    print("Cumulative-EB closure cannot select this: free-drift TDEE absorbs 99.25% of any level change.\n")

    summary = []
    for name, qb, st, mdt in CONFIGS:
        lls, fmaes, ratios, tdees = [], [], [], []
        for m in GRID:
            ll, fmae, ratio, tdee = run_one(data, m, qb, st, mdt)
            lls.append(ll); fmaes.append(fmae); ratios.append(ratio); tdees.append(tdee)
        lls = np.array(lls)
        rel = lls - lls.max()
        b = int(np.argmax(lls))
        inside = GRID[rel >= -CHI2_95]
        print(f"{name:<12} Q_BUF={qb:<5} sigma_TDEE={st:<3}")
        print(f"    m_hat = {GRID[b]:.2f}   95% profile CI [{inside.min():.2f}, {inside.max():.2f}]"
              f"   {'INTERIOR max' if 0 < b < len(GRID) - 1 else 'BOUNDARY (not identified)'}")
        print(f"    at m_hat:  FM_scan_MAE = {fmaes[b]:.3f} lb   "
              f"TDEE = {tdees[b]:.0f} cal/day   TDEE/RMR = {ratios[b]:.2f}")
        print()
        summary.append(dict(tuning=name, q_buf=qb, sigma_tdee=st, max_dtdee=mdt, m_hat=float(GRID[b]),
                            ci_lo=float(inside.min()), ci_hi=float(inside.max()),
                            interior=bool(0 < b < len(GRID) - 1), fm_scan_mae=fmaes[b],
                            tdee_rmr=ratios[b], tdee=tdees[b]))

    P4.Q_BUF = q_buf_0
    save("CR_multiplier_likelihood", by_tuning=summary, canonical_m=summary[0]["m_hat"],
         canonical_ci=[summary[0]["ci_lo"], summary[0]["ci_hi"]])
    print("Reading: tightening the noise states narrows the interval at a real cost in")
    print("composition-scan fit; the canonical tuning is the one that fits the scans, so its")
    print("wider interval is the honest one. The form of the correction (uniform vs era,")
    print("level, macro) is chosen in DM_intake_calibration.py; this script is the")
    print("tuning-sensitivity check for the uniform multiplier only.")
    print("Compare: python analysis/CQ_intake_scale_response.py (regression outside the filter)")
    print("Wrote analysis/results/CR_multiplier_likelihood.json")


if __name__ == "__main__":
    main()
