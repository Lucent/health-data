"""CJ. Fine Atwater-factor sweep with EB-closing multiplier.

Two improvements over CI:

  1. 0.05 cal/g resolution (vs CI's 0.25). Granularity to detect a real per-gram
     bias of a tenth of a calorie.

  2. EB-closing multiplier per config. CI held the canonical 1.09 intake multiplier
     fixed and let cumulative EB drift, which confounded "shape" (macro distribution
     bias) with "level" (overall undercount). Here the multiplier is solved per
     config so 15-year cumulative EB closes to ~0, and FM_MAE reports the clean
     daily-trajectory cost of that macro shape.

Output: CSV plus a 3-panel line plot (FM_MAE and EB-closing multiplier vs each
macro value). The macro whose FM_MAE curve has a clear minimum away from
textbook indicates a real per-gram efficiency offset; the multiplier curve
shows how much the per-gram shift trades off against the global level.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))

import P4_kalman_filter as P4
from P4_kalman_filter import (
    kalman_v4_forward, rts_smoother_v4, diagnostics, build_inputs, load_data,
)

# Cache loaded data so each Kalman run isn't re-reading CSVs.
_DATA = None
def _load_once():
    global _DATA
    if _DATA is None:
        _DATA = load_data()
    return _DATA


def make_intake(a_p, b_c, c_f, multiplier):
    """Adjust intake calories under (a, b, c) Atwater factors and a global multiplier."""
    _, intake_daily, _, _, _, _, _ = _load_once()
    intake = intake_daily.copy()
    adj = ((a_p - 4.0) * intake.protein_g.fillna(0) +
           (b_c - 4.0) * intake.carbs_g.fillna(0) +
           (c_f - 9.0) * intake.fat_g.fillna(0))
    intake["calories"] = ((intake["calories"] + adj) * multiplier).clip(lower=0)
    return intake


def run_kalman(intake_alt):
    smooth, _, strength, comp_scans, rmr_df, p2, temp_df = _load_once()
    # Disable build_inputs' own multiplier — we apply it manually in make_intake.
    P4.INTAKE_MULTIPLIER = 1.0
    dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk = \
        build_inputs(smooth, intake_alt, strength, comp_scans, p2)
    first_scan = np.where(~np.isnan(fm_scan))[0][0]
    fm_start = fm_scan[first_scan]
    lm_start = lm_scan[first_scan]
    x_filt, P_filt, x_pred, P_pred, innov, _ = kalman_v4_forward(
        calories, weight_obs, fm_scan, lm_scan, expected_rmr,
        workouts_per_wk, fm_start, lm_start, 50, 100,
    )
    x_smooth, _ = rts_smoother_v4(
        x_filt, P_filt, x_pred, P_pred, first_scan, workouts_per_wk)
    return diagnostics(dates, x_smooth, innov, fm_scan, lm_scan,
                       calories, rmr_df, weight_obs, temp_df)


def find_eb_closing(a_p, b_c, c_f, m_init=1.09, tol=0.02):
    """Return (multiplier, diagnostics) at which 15-yr EB residual ≈ 0."""
    intake1 = make_intake(a_p, b_c, c_f, m_init)
    d1 = run_kalman(intake1)
    eb1 = d1["eb_residual_lb"]
    if abs(eb1) < tol:
        return m_init, d1
    # Linear interpolation: probe a second multiplier and solve EB(m) = 0
    m2 = m_init - 0.02 * np.sign(eb1)
    intake2 = make_intake(a_p, b_c, c_f, m2)
    d2 = run_kalman(intake2)
    eb2 = d2["eb_residual_lb"]
    if abs(eb2 - eb1) < 1e-6:
        return (m_init if abs(eb1) < abs(eb2) else m2), (d1 if abs(eb1) < abs(eb2) else d2)
    m_zero = m_init + (m2 - m_init) * (-eb1) / (eb2 - eb1)
    intake3 = make_intake(a_p, b_c, c_f, m_zero)
    d3 = run_kalman(intake3)
    # If still off, one more bisection step
    eb3 = d3["eb_residual_lb"]
    if abs(eb3) > tol:
        # Pick the closest pair and re-interpolate
        pairs = sorted([(eb1, m_init), (eb2, m2), (eb3, m_zero)], key=lambda p: abs(p[0]))
        (eba, ma), (ebb, mb) = pairs[0], pairs[1]
        if abs(ebb - eba) > 1e-6:
            m_zero = ma + (mb - ma) * (-eba) / (ebb - eba)
            intake4 = make_intake(a_p, b_c, c_f, m_zero)
            d3 = run_kalman(intake4)
    return m_zero, d3


def sweep_one(label, fixed, var_values, var_pos):
    """1D sweep of one macro at 0.05 res, others at textbook (4, 4, 9).

    var_pos: 0=protein, 1=carb, 2=fat. fixed: tuple of the other two textbook values.
    """
    rows = []
    for v in var_values:
        if var_pos == 0:
            a, b, c = v, 4.0, 9.0
        elif var_pos == 1:
            a, b, c = 4.0, v, 9.0
        else:
            a, b, c = 4.0, 4.0, v
        m_zero, d = find_eb_closing(a, b, c)
        rows.append({
            "macro": label, "value": round(v, 3),
            "multiplier": round(m_zero, 5),
            "fm_mae": d["fm_scan_mae"],
            "lm_mae": d["lm_scan_mae"],
            "eb_residual_lb": d["eb_residual_lb"],
            "innov_acf1": d["innov_acf1"],
            "innov_std": d["innov_std"],
        })
        print(f"  {label} {v:>5.2f}  m={m_zero:.5f}  FM_MAE={d['fm_scan_mae']:.4f}  EB={d['eb_residual_lb']:+.3f}")
    return pd.DataFrame(rows)


def joint_optimum(start=(4.0, 4.0, 9.0)):
    """Nelder-Mead optimization over (a, b, c) with EB-closing multiplier inside."""
    def objective(x):
        a, b, c = x
        # Plausibility bounds: keep within ±2 cal/g of textbook
        if not (3.0 <= a <= 5.5 and 3.0 <= b <= 5.5 and 7.0 <= c <= 11.0):
            return 1e6
        m, d = find_eb_closing(a, b, c)
        return d["fm_scan_mae"]
    print(f"\nNelder-Mead joint optimization starting from {start}...")
    res = minimize(objective, start, method="Nelder-Mead",
                   options={"xatol": 0.02, "fatol": 1e-5, "maxiter": 80, "disp": True})
    a, b, c = res.x
    m, d = find_eb_closing(a, b, c)
    return a, b, c, m, d


def main():
    print("=" * 80)
    print("CJ. Fine Atwater sweep, 0.05 cal/g, EB-closing multiplier per config")
    print("=" * 80)

    # Reference
    print("\nReference: textbook (4, 4, 9) at canonical m=1.09")
    m_ref, d_ref = find_eb_closing(4.0, 4.0, 9.0)
    print(f"  EB-closing m: {m_ref:.5f}  FM_MAE: {d_ref['fm_scan_mae']:.4f}  EB: {d_ref['eb_residual_lb']:+.3f}")

    print("\n--- Protein 3.00 → 5.50 step 0.05 ---")
    df_p = sweep_one("protein", (4.0, 9.0), np.arange(3.00, 5.51, 0.05), 0)

    print("\n--- Carbs 3.00 → 5.50 step 0.05 ---")
    df_c = sweep_one("carb", (4.0, 9.0), np.arange(3.00, 5.51, 0.05), 1)

    print("\n--- Fat 7.00 → 10.50 step 0.05 ---")
    df_f = sweep_one("fat", (4.0, 4.0), np.arange(7.00, 10.51, 0.05), 2)

    # Joint Nelder-Mead optimum
    a_opt, b_opt, c_opt, m_opt, d_opt = joint_optimum()
    print(f"\nJoint optimum:")
    print(f"  Atwater (P, C, F):  ({a_opt:.3f}, {b_opt:.3f}, {c_opt:.3f})")
    print(f"  Multiplier:          {m_opt:.5f}")
    print(f"  FM_MAE:              {d_opt['fm_scan_mae']:.4f}")
    print(f"  EB residual:         {d_opt['eb_residual_lb']:+.3f} lb")
    print(f"  Δ vs textbook:       {d_opt['fm_scan_mae'] - d_ref['fm_scan_mae']:+.4f} lb")

    out = pd.concat([df_p, df_c, df_f], ignore_index=True)
    out_path = ROOT / "analysis" / "CJ_atwater_fine.csv"
    out.to_csv(out_path, index=False)
    print(f"\nSaved {out_path}")

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=False)
    panels = [("protein", df_p, 4.0), ("carb", df_c, 4.0), ("fat", df_f, 9.0)]
    for col, (label, df, textbook) in enumerate(panels):
        ax_top = axes[0, col]
        ax_bot = axes[1, col]
        ax_top.plot(df["value"], df["fm_mae"], "o-", color="tab:blue")
        i_min = df["fm_mae"].idxmin()
        ax_top.axvline(df.loc[i_min, "value"], color="tab:red", ls="--",
                       label=f"min @ {df.loc[i_min, 'value']:.2f}")
        ax_top.axvline(textbook, color="black", ls=":", label=f"textbook {textbook}")
        ax_top.set_title(f"{label}: FM_MAE")
        ax_top.set_xlabel("cal/g")
        ax_top.set_ylabel("FM_MAE (lb)")
        ax_top.legend(fontsize=8)
        ax_top.grid(alpha=0.3)

        ax_bot.plot(df["value"], df["multiplier"], "s-", color="tab:green")
        ax_bot.axvline(textbook, color="black", ls=":")
        ax_bot.set_title(f"{label}: EB-closing multiplier")
        ax_bot.set_xlabel("cal/g")
        ax_bot.set_ylabel("multiplier")
        ax_bot.grid(alpha=0.3)

    fig.suptitle("CJ — fine Atwater sweep, EB-closed multiplier per config", y=1.00)
    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CJ_atwater_fine.png"
    fig.savefig(plot_path, dpi=140, bbox_inches="tight")
    print(f"Saved {plot_path}")

    # Headline
    best_p = df_p.loc[df_p["fm_mae"].idxmin()]
    best_c = df_c.loc[df_c["fm_mae"].idxmin()]
    best_f = df_f.loc[df_f["fm_mae"].idxmin()]
    print("\n--- Best per-macro (FM_MAE) ---")
    print(f"  Protein: {best_p['value']:.2f} cal/g  (Δ={best_p['value']-4.0:+.2f}), "
          f"FM_MAE={best_p['fm_mae']:.4f}, m={best_p['multiplier']:.4f}")
    print(f"  Carb:    {best_c['value']:.2f} cal/g  (Δ={best_c['value']-4.0:+.2f}), "
          f"FM_MAE={best_c['fm_mae']:.4f}, m={best_c['multiplier']:.4f}")
    print(f"  Fat:     {best_f['value']:.2f} cal/g  (Δ={best_f['value']-9.0:+.2f}), "
          f"FM_MAE={best_f['fm_mae']:.4f}, m={best_f['multiplier']:.4f}")
    print(f"\n  Reference textbook FM_MAE: {d_ref['fm_scan_mae']:.4f}")


if __name__ == "__main__":
    main()
