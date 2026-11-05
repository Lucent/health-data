"""DL: Kevin Hall's NIH body-weight dynamic model, run FORWARD on 15 years of logged intake.

A physiologically constrained alternative to P4's free-drift Kalman TDEE. The model is
driven by logged intake only (x intake multiplier) from the first composition scan
(2011-04-28, FM 99.9 lb, lean 159.3 lb); nothing after day 0 is assimilated.

Model (Hall 2010 AJP-Endo 298:E449 "Predicting metabolic adaptation, body weight change,
and energy intake in humans"; Hall et al. 2011 Lancet 378:826 web appendix; NIH Body
Weight Planner). Constants used (source values, kcal/kg unless noted):
  rho_F = 9440, rho_L = 1807 (energy density of fat / lean tissue change)
  RMR = K + 3.2*F + 22*L  (gamma_F, gamma_L; K set so RMR at baseline = Mifflin-St Jeor)
  physical activity = delta*BW, delta = ((1-0.1)*PAL - 1)*RMR_b/BW_b  (Lancet appendix)
  TEF = 0.10 * EI
  adaptive thermogenesis: tau*dAT/dt = beta_AT*(EI - EI_b) - AT, beta_AT = 0.14, tau = 14 d
  turnover costs eta_F = 180, eta_L = 230 (EE depends on dF/dt, dL/dt; solved algebraically)
  Forbes/Hall partition: p = C/(C+F), C = 10.4 kg * rho_L/rho_F (fraction of imbalance to lean)
  glycogen: dG/dt = (CI - k_G*G^2)/rho_G, k_G = CI_b/G_b^2, G_b = 0.5 kg, rho_G = 4180, 2.7 g water/g
  ECW: dECW/dt = (dNa - 3000*ECW - 4000*(1 - CI/CI_b)) / 3220   (mg/L/day, mg/day, mg/L)
Approximations stated here: EI_b (baseline maintenance intake) = PAL*RMR_b as in the Body
Weight Planner (the pre-diet intake is not logged); CI_b = record-mean carb fraction * EI_b;
Na_b = record-mean logged sodium (keeps ECW centred); glycogen energy (rho_G*dG) is debited
from the F/L imbalance for conservation; Euler daily steps, 4 sub-steps for G/ECW; no ageing.

Command: python analysis/DL_hall_model.py
Artifact: analysis/DL_hall_daily.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent.parent
KG = 0.45359237
RHO_F, RHO_L, RHO_G = 9440.0, 1807.0, 4180.0
GAM_F, GAM_L = 3.2, 22.0
ETA_F, ETA_L = 180.0, 230.0
B_TEF, TAU_AT = 0.10, 14.0
FORBES_C = 10.4 * RHO_L / RHO_F
G_B, GLY_WATER = 0.5, 2.7
NA_CONC, XI_NA, XI_CI = 3220.0, 3000.0, 4000.0
SUBSTEPS = 4
PAL_DEFAULT, BETA_DEFAULT, M_CANON = 1.5, 0.14, 1.09
HEIGHT_CM, BIRTH = 182.25, pd.Timestamp("1982-10-15")
START = pd.Timestamp("2011-04-28")
F0_LB, L0_LB = 99.9, 159.3
ERAS = {"loss 2011-14": ("2011-04-28", "2014-12-31"),
        "regain 2015-24": ("2015-01-01", "2024-09-16"),
        "tirzepatide 2024-26": ("2024-09-17", "2026-02-28"),
        "  on-drug plateau 2025-10..": ("2025-10-01", "2026-02-28")}
PLATEAU = ("2025-10-01", "2026-02-28")


# ── data ──────────────────────────────────────────────────────────────
def load():
    d = {}
    it = pd.read_csv(ROOT / "intake/intake_daily.csv", parse_dates=["date"])
    d["dates"] = it["date"].values
    d["idx"] = {pd.Timestamp(x): i for i, x in enumerate(d["dates"])}
    d["cal"], d["carb"], d["na"] = (it["calories"].values.astype(float),
                                    it["carbs_g"].values * 4.0, it["sodium_mg"].values.astype(float))
    d["i0"] = d["idx"][START]
    w = pd.read_csv(ROOT / "weight/weight.csv", parse_dates=["date"])
    w = w[(w["date"] >= START) & (w["date"] <= it["date"].iloc[-1])]
    d["w_idx"] = np.array([d["idx"][t] for t in w["date"]]); d["w_obs"] = w["weight_lbs"].values
    p1 = pd.read_csv(ROOT / "analysis/P1_smoothed_weight.csv", parse_dates=["date"])
    p1 = p1[(p1["date"] >= START) & p1["smoothed_weight_lbs"].notna()]
    d["p1_idx"] = np.array([d["idx"][t] for t in p1["date"]]); d["p1_obs"] = p1["smoothed_weight_lbs"].values
    sc = pd.read_csv(ROOT / "composition/composition.csv", parse_dates=["date"])
    d["s_idx"] = np.array([d["idx"][t] for t in sc["date"]])
    d["s_fm"] = sc["fat_mass_lbs"].values; d["s_lm"] = (sc["weight_lbs"] - sc["fat_mass_lbs"]).values
    r = pd.read_csv(ROOT / "RMR/rmr.csv", parse_dates=["date"])
    d["r_idx"] = np.array([d["idx"][t] for t in r["date"]]); d["r_obs"] = r["rmr_kcal"].values.astype(float)
    d["r_dev"] = r["device"].values
    k = pd.read_csv(ROOT / "analysis/P4_kalman_daily.csv", parse_dates=["date"]).set_index("date")
    k = k.reindex(d["dates"])
    d["k_fm"], d["k_lm"], d["k_buf"], d["k_tdee"] = (k["fat_mass_lbs"].values, k["lm_lbs"].values,
                                                     k["buf_lbs"].values, k["tdee"].values)
    p2 = pd.read_csv(ROOT / "analysis/P2_daily_composition.csv", parse_dates=["date"]).set_index("date")
    p2 = p2.reindex(d["dates"])
    d["p2_fm"], d["p2_ffm"], d["p2_rmr"] = p2["fm_lbs"].values, p2["ffm_lbs"].values, p2["expected_rmr"].values
    d["carb_frac"] = d["carb"].sum() / d["cal"].sum()
    d["na_mean"] = d["na"].mean()
    return d


def baseline(pal, m, d):
    F0, L0 = F0_LB * KG, L0_LB * KG
    age = (START - BIRTH).days / 365.25
    rmr_b = 10 * (F0 + L0) + 6.25 * HEIGHT_CM - 5 * age + 5          # Mifflin-St Jeor, male
    b = dict(F0=F0, L0=L0, rmr_b=rmr_b, K=rmr_b - GAM_F * F0 - GAM_L * L0, ei_b=pal * rmr_b,
             delta=((1 - B_TEF) * pal - 1) * rmr_b / (F0 + L0))
    b["ci_b"] = d["carb_frac"] * b["ei_b"]; b["kG"] = b["ci_b"] / G_B ** 2; b["na_b"] = m * d["na_mean"]
    return b


def water_step(G, W, ci, dna, b):
    """Hall glycogen + extracellular-water dynamics over one day (sub-stepped)."""
    h = 1.0 / SUBSTEPS
    for _ in range(SUBSTEPS):
        G += h * (ci - b["kG"] * G * G) / RHO_G
        W += h * (dna - XI_NA * W - XI_CI * (1 - ci / b["ci_b"])) / NA_CONC
    return G, W


# ── forward model ────────────────────────────────────────────────────
def simulate(pal, beta, m, d, i0=None, F_init=None, L_init=None, at_init=0.0):
    """Forward Euler run. Arrays indexed by date; state at index t = morning of date t,
    intake of day t moves the state to t+1 (project convention)."""
    b = baseline(pal, m, d)
    n = len(d["cal"]); i0 = d["i0"] if i0 is None else i0
    F = b["F0"] if F_init is None else F_init; L = b["L0"] if L_init is None else L_init
    G, W, AT = G_B, 0.0, at_init
    cal, carb, na = d["cal"], d["carb"], d["na"]
    out = {k: np.full(n, np.nan) for k in ("bw", "fm", "ffm", "ee", "rmr", "at", "pa", "gly", "ecw")}
    K, delta, ei_b, na_b = b["K"], b["delta"], b["ei_b"], b["na_b"]
    for t in range(i0, n):
        bw = F + L + (1 + GLY_WATER) * (G - G_B) + W
        rmr = K + GAM_F * F + GAM_L * L
        out["bw"][t], out["fm"][t], out["ffm"][t] = bw / KG, F / KG, (bw - F) / KG
        out["rmr"][t], out["at"][t], out["gly"][t], out["ecw"][t] = rmr, AT, G, W
        if t == n - 1:
            break
        ei, ci, dna = m * cal[t], m * carb[t], m * na[t] - na_b
        G1, W1 = water_step(G, W, ci, dna, b)
        dG = G1 - G
        p = FORBES_C / (FORBES_C + F)
        q = (1 - p) * ETA_F / RHO_F + p * ETA_L / RHO_L
        pa = delta * bw
        e0 = rmr + pa + B_TEF * ei + AT
        ee = (e0 + q * (ei - RHO_G * dG)) / (1 + q)
        avail = ei - ee - RHO_G * dG
        out["ee"][t], out["pa"][t] = ee, pa
        F = max(F + (1 - p) * avail / RHO_F, 2.0)
        L = max(L + p * avail / RHO_L, 30.0)
        AT += (beta * (ei - ei_b) - AT) / TAU_AT
        G, W = G1, W1
    out["b"] = b
    return out


def invert(pal, beta, d, m_hyd=M_CANON, smooth_days=15):
    """Polidori/Hall 2016 inversion: intake implied by the observed weight path.
    Observed raw weigh-ins -> daily interpolation -> minus Hall's own glycogen/ECW
    deviations (driven by logged carbs/sodium) -> smoothed -> Forbes-partitioned into
    F/L -> EI = (RMR + PA + AT + eta*dM + rho*dM + rho_G*dG) / (1 - beta_TEF)."""
    b = baseline(pal, m_hyd, d)
    n = len(d["cal"]); i0 = d["i0"]
    gly, ecw = np.full(n, np.nan), np.full(n, np.nan)
    G, W = G_B, 0.0
    for t in range(i0, n):
        gly[t], ecw[t] = G, W
        if t < n - 1:
            G, W = water_step(G, W, m_hyd * d["carb"][t], m_hyd * d["na"][t] - b["na_b"], b)
    days = np.arange(n)
    bw_raw = np.interp(days, d["w_idx"], d["w_obs"]) * KG
    water = (1 + GLY_WATER) * (gly - G_B) + ecw
    M = pd.Series(bw_raw - water).rolling(smooth_days, center=True, min_periods=1).mean().values
    F, L, AT = b["F0"], b["L0"], 0.0
    ei_impl, fm_path = np.full(n, np.nan), np.full(n, np.nan)
    for t in range(i0, n - 1):
        fm_path[t] = F / KG
        p = FORBES_C / (FORBES_C + F)
        dM = M[t + 1] - M[t]
        avail = dM / (p / RHO_L + (1 - p) / RHO_F)
        dF, dL = (1 - p) * avail / RHO_F, p * avail / RHO_L
        dG = gly[t + 1] - gly[t]
        e0 = b["K"] + GAM_F * F + GAM_L * L + b["delta"] * bw_raw[t] + AT
        ei = (e0 + ETA_F * dF + ETA_L * dL + avail + RHO_G * dG) / (1 - B_TEF)
        ei_impl[t] = ei
        AT += (beta * (ei - b["ei_b"]) - AT) / TAU_AT
        F, L = F + dF, L + dL
    return ei_impl, fm_path, b


# ── scoring ──────────────────────────────────────────────────────────
def metrics(s, d):
    e_w = s["bw"][d["w_idx"]] - d["w_obs"]
    e_p1 = s["bw"][d["p1_idx"]] - d["p1_obs"]
    rmr_at = s["rmr"] + s["at"]
    return dict(w_rmse=np.sqrt(np.mean(e_w ** 2)), w_mean=e_w.mean(), w_max=np.abs(e_w).max(),
                w_end=e_w[-1], p1_rmse=np.sqrt(np.mean(e_p1 ** 2)),
                fm_mae=np.mean(np.abs(s["fm"][d["s_idx"]] - d["s_fm"])),
                ffm_mae=np.mean(np.abs(s["ffm"][d["s_idx"]] - d["s_lm"])),
                rmr_mae=np.mean(np.abs(s["rmr"][d["r_idx"]] - d["r_obs"])),
                rmr_at_mae=np.mean(np.abs(rmr_at[d["r_idx"]] - d["r_obs"])),
                rmr_bias=np.mean(s["rmr"][d["r_idx"]] - d["r_obs"]),
                tdee_rmr=np.mean(s["ee"][d["r_idx"]] / d["r_obs"]))


PARAMS = {"PAL": (PAL_DEFAULT, 1.0, 2.2), "beta_AT": (BETA_DEFAULT, 0.0, 0.5), "m": (M_CANON, 0.9, 1.4)}


def fit(d, free, fixed=None):
    """Least squares on raw weigh-ins only. `free` = names of parameters to fit; the rest
    stay at Hall defaults (or `fixed`). SEs: naive OLS and inflated by n/n_eff from the
    residual lag-1 autocorrelation (forward-run residuals are nearly a random walk)."""
    fixed = {**{k: v[0] for k, v in PARAMS.items()}, **(fixed or {})}

    def unpack(th):
        p = dict(fixed); p.update(zip(free, th))
        return p["PAL"], p["beta_AT"], p["m"]

    def resid(th):
        return simulate(*unpack(th), d)["bw"][d["w_idx"]] - d["w_obs"]
    r = least_squares(resid, [fixed[k] for k in free], bounds=([PARAMS[k][1] for k in free],
                      [PARAMS[k][2] for k in free]), diff_step=1e-3)
    n, k = len(r.fun), len(r.x)
    cov = np.linalg.inv(r.jac.T @ r.jac) * (r.fun ** 2).sum() / (n - k)
    se = np.sqrt(np.diag(cov))
    rr = np.corrcoef(r.fun[:-1], r.fun[1:])[0, 1]
    n_eff = n * (1 - rr) / (1 + rr)
    return unpack(r.x), dict(zip(free, se)), dict(zip(free, se * np.sqrt(n / max(n_eff, k + 1)))), rr, n_eff, \
        np.sqrt(np.mean(r.fun ** 2))


def era_mask(d, lo, hi):
    return (d["dates"] >= np.datetime64(lo)) & (d["dates"] <= np.datetime64(hi))


def slope_per_30d(idx, y):
    ok = ~np.isnan(y)
    A = np.column_stack([idx[ok], np.ones(ok.sum())])
    return np.linalg.lstsq(A, y[ok], rcond=None)[0][0] * 30


# ── main ─────────────────────────────────────────────────────────────
def main():
    d = load()
    n = len(d["cal"]); dates = d["dates"]
    print(f"Hall model forward run: {str(dates[d['i0']])[:10]} .. {str(dates[-1])[:10]} "
          f"({n - d['i0']} days), {len(d['w_obs'])} weigh-ins, {len(d['s_fm'])} scans, "
          f"{len(d['r_obs'])} calorimetry RMRs")
    b0 = baseline(PAL_DEFAULT, M_CANON, d)
    print(f"Baseline (Hall): Mifflin RMR_b={b0['rmr_b']:.0f}, K={b0['K']:.0f}, measured RMR 2011-04-30 = "
          f"{d['r_obs'][0]:.0f}; EI_b=PAL*RMR_b; carb frac={d['carb_frac']:.2f}; Na_b={b0['na_b']:.0f} mg")

    # Kalman reference numbers (P4) on the same targets
    kw = d["k_fm"] + d["k_lm"] + d["k_buf"]
    kal = dict(p1_rmse=np.sqrt(np.nanmean((kw[d["p1_idx"]] - d["p1_obs"]) ** 2)),
               fm_mae=np.nanmean(np.abs(d["k_fm"][d["s_idx"]] - d["s_fm"])),
               ffm_mae=np.nanmean(np.abs(d["k_lm"][d["s_idx"]] - d["s_lm"])),
               tdee_rmr=np.nanmean(d["k_tdee"][d["r_idx"]] / d["r_obs"]))

    # ── 0. Hall defaults, 1. fits ──
    runs = {"Hall defaults (PAL 1.5, beta .14, m 1.09)": (PAL_DEFAULT, BETA_DEFAULT, M_CANON)}
    print("\n=== Fit on weigh-ins only (scans + calorimetry held out) ===")
    fits = {}
    for label, free in (("fit PAL | beta=.14, m=1.09", ["PAL"]), ("fit PAL,beta | m=1.09", ["PAL", "beta_AT"]),
                        ("fit PAL,beta,m", ["PAL", "beta_AT", "m"])):
        th, se, se_eff, rr, n_eff, rm = fit(d, free)
        fits[label] = th; runs[label] = th
        vals = dict(zip(["PAL", "beta_AT", "m"], th))
        txt = "  ".join(f"{k}={vals[k]:.3f}±{se[k]:.3f} (acf-adj ±{se_eff[k]:.3f})" for k in free)
        bnd = [k for k in free if min(abs(vals[k] - PARAMS[k][1]), abs(vals[k] - PARAMS[k][2])) < 1e-3]
        print(f"  {label}: {txt}   [W rmse {rm:.2f}; resid lag-1 ACF {rr:.2f}, n_eff≈{n_eff:.0f}"
              f"{'; AT BOUND: ' + ','.join(bnd) if bnd else ''}]")
    print("  beta_AT profile (PAL refit at each beta, m=1.09) — the PAL/beta ridge:")
    for beta in (0.0, 0.14, 0.3, 0.5):
        th, _, _, _, _, rm = fit(d, ["PAL"], {"beta_AT": beta})
        s = simulate(*th, d)
        print(f"    beta={beta:.2f}: PAL={th[0]:.3f}, W rmse {rm:.2f}, mean AT {np.nanmean(s['at']):+.0f} kcal/d, "
              f"EI_b={s['b']['ei_b']:.0f}, delta={s['b']['delta']:.1f} kcal/kg/d")
    print("\n=== Forward-run errors (nothing assimilated after 2011-04-28) ===")
    hdr = f"{'run':<40}{'W rmse':>8}{'W mean':>8}{'W max':>7}{'W end':>7}{'P1 rmse':>8}{'FM mae':>7}" \
          f"{'FFM mae':>8}{'RMR mae':>8}{'+AT':>6}{'RMRbias':>8}{'EE/RMR':>7}"
    print(hdr)
    sims = {}
    for label, th in runs.items():
        s = simulate(*th, d); sims[label] = s; mt = metrics(s, d)
        print(f"{label:<40}{mt['w_rmse']:>8.2f}{mt['w_mean']:>+8.2f}{mt['w_max']:>7.1f}{mt['w_end']:>+7.1f}"
              f"{mt['p1_rmse']:>8.2f}{mt['fm_mae']:>7.2f}{mt['ffm_mae']:>8.2f}{mt['rmr_mae']:>8.0f}"
              f"{mt['rmr_at_mae']:>6.0f}{mt['rmr_bias']:>+8.0f}{mt['tdee_rmr']:>7.2f}")
    print(f"{'Kalman P4 (assimilates weights+scans)':<40}{'-':>8}{'-':>8}{'-':>7}{'-':>7}{kal['p1_rmse']:>8.2f}"
          f"{kal['fm_mae']:>7.2f}{kal['ffm_mae']:>8.2f}{'-':>8}{'-':>6}{'-':>8}{kal['tdee_rmr']:>7.2f}")
    print("  (W = raw weigh-ins, lb; P1 = glycogen/sodium-corrected weigh-ins; FM/FFM at 70 scans, lb;"
          "\n   RMR = Hall K+3.2F+22L vs 25 calorimetry values, kcal; +AT adds adaptive thermogenesis;"
          "\n   EE/RMR = model TDEE / measured RMR, cf. Kalman TDEE/RMR 1.28)")

    prim_label = "fit PAL | beta=.14, m=1.09"
    prim, alt = sims[prim_label], sims["fit PAL,beta | m=1.09"]
    pal_p, beta_p, _ = fits[prim_label]
    yrs = pd.DatetimeIndex(dates).year
    print("\nYear-by-year forward error (sim − observed weigh-ins, lb), primary fit:")
    ew = prim["bw"][d["w_idx"]] - d["w_obs"]; wy = yrs[d["w_idx"]]
    print("  " + "  ".join(f"{y}:{ew[wy == y].mean():+.1f}" for y in sorted(set(wy))))

    # ── 2. inversion: implied intake vs logged ──
    print("\n=== Implied intake (Hall inverted on the observed weight path) vs logged ===")
    ei_def, fm_def, _ = invert(PAL_DEFAULT, BETA_DEFAULT, d)
    ei_fit, fm_fit, _ = invert(pal_p, beta_p, d)
    print(f"  inverted-path FM error at scans: PAL 1.5 MAE {np.nanmean(np.abs(fm_def[d['s_idx']] - d['s_fm'])):.1f} lb,"
          f" fitted PAL MAE {np.nanmean(np.abs(fm_fit[d['s_idx']] - d['s_fm'])):.1f} lb (Forbes drift, no scan anchor)")
    print(f"  {'year':<6}{'logged':>8}{'x1.09':>8}{'impl@1.5':>10}{'ratio':>7}{'impl@fit':>10}{'ratio':>7}{'days':>6}")
    rows = []
    for y in sorted(set(yrs[d["i0"]:n - 1])):
        mk = (yrs == y) & ~np.isnan(ei_def)
        lg, i1, i2 = d["cal"][mk].mean(), ei_def[mk].mean(), ei_fit[mk].mean()
        rows.append((y, lg, i1, i2, mk.sum()))
        print(f"  {y:<6}{lg:>8.0f}{lg * M_CANON:>8.0f}{i1:>10.0f}{i1 / lg:>7.3f}{i2:>10.0f}{i2 / lg:>7.3f}{mk.sum():>6}")
    mk = ~np.isnan(ei_def)
    lg = d["cal"][mk].mean()
    print(f"  {'all':<6}{lg:>8.0f}{lg * M_CANON:>8.0f}{ei_def[mk].mean():>10.0f}{ei_def[mk].mean() / lg:>7.3f}"
          f"{ei_fit[mk].mean():>10.0f}{ei_fit[mk].mean() / lg:>7.3f}{mk.sum():>6}")
    r_y = np.array([r[2] / r[1] for r in rows]); r_f = np.array([r[3] / r[1] for r in rows])
    print(f"  implied/logged ratio across years: PAL 1.5 mean {r_y.mean():.3f} sd {r_y.std():.3f} "
          f"range [{r_y.min():.3f},{r_y.max():.3f}]; fitted PAL mean {r_f.mean():.3f} sd {r_f.std():.3f}")
    # Polidori/Hall 2016 appetite feedback: dEI/dBW relative to baseline (they report ~ -95 kcal/d per kg lost)
    mo = pd.DataFrame({"ei": ei_fit, "lg": d["cal"] * M_CANON, "bw": prim["bw"] * KG},
                      index=pd.DatetimeIndex(dates)).dropna().resample("MS").mean()
    bw_obs_m = pd.Series(np.interp(np.arange(n), d["w_idx"], d["w_obs"]) * KG,
                         index=pd.DatetimeIndex(dates)).iloc[d["i0"]:].resample("MS").mean()
    x = (bw_obs_m - (F0_LB + L0_LB) * KG).reindex(mo.index).values
    bI = np.polyfit(x, mo["ei"].values - ei_fit[mk].mean(), 1)[0]; bL = np.polyfit(x, mo["lg"].values - mo["lg"].mean(), 1)[0]
    print(f"  appetite-feedback slope dEI/dBW (monthly, vs baseline): implied {bI:+.0f}, logged x1.09 {bL:+.0f} kcal/d per kg"
          f"  (Polidori 2016 free-living appetite: ≈ −95 kcal/d per kg; positive here = intake rises with weight)")

    # ── 3. adaptive thermogenesis by era vs Kalman ──
    print("\n=== Adaptive thermogenesis by era (kcal/day) ===")
    bb = prim["b"]
    kal_struct = (bb["K"] + GAM_F * d["k_fm"] * KG + GAM_L * d["k_lm"] * KG
                  + bb["delta"] * (d["k_fm"] + d["k_lm"]) * KG + B_TEF * M_CANON * d["cal"])
    kal_at = d["k_tdee"] - kal_struct
    print(f"  {'era':<30}{'Hall AT':>9}{'Hall AT(β fit)':>16}{'Hall EE':>9}{'Kalman TDEE':>12}{'Kal−struct':>11}"
          f"{'EE−TDEE':>9}{'EI x1.09':>9}")
    for name, (lo, hi) in ERAS.items():
        mk = era_mask(d, lo, hi) & ~np.isnan(d["k_tdee"]) & ~np.isnan(prim["ee"])
        print(f"  {name:<30}{prim['at'][mk].mean():>9.0f}{alt['at'][mk].mean():>16.0f}{prim['ee'][mk].mean():>9.0f}"
              f"{d['k_tdee'][mk].mean():>12.0f}{kal_at[mk].mean():>11.0f}"
              f"{(prim['ee'][mk] - d['k_tdee'][mk]).mean():>+9.0f}{(M_CANON * d['cal'][mk]).mean():>9.0f}")
    print("  (Hall AT = beta*(EI − EI_b) lagged 14 d, EI_b = PAL*Mifflin at 2011 baseline; Kal−struct = Kalman TDEE"
          "\n   minus Hall's structural EE [K+3.2F+22L+delta*BW+0.1EI] evaluated at the Kalman's own FM/LM)")

    print("\n=== Calorimetry residuals: measured RMR − RMR predicted from scan-interpolated composition ===")
    rmr_hall_p2 = bb["K"] + GAM_F * d["p2_fm"] * KG + GAM_L * d["p2_ffm"] * KG
    print(f"  {'date':<12}{'device':<16}{'meas':>6}{'Hall(F,L)':>10}{'resid':>7}{'P2 fit':>8}{'resid':>7}"
          f"{'Hall AT':>9}{'Kal−struct':>11}")
    res_h, res_p = [], []
    for i, ri in enumerate(d["r_idx"]):
        rh, rp = d["r_obs"][i] - rmr_hall_p2[ri], d["r_obs"][i] - d["p2_rmr"][ri]
        res_h.append(rh); res_p.append(rp)
        print(f"  {str(dates[ri])[:10]:<12}{d['r_dev'][i]:<16}{d['r_obs'][i]:>6.0f}{rmr_hall_p2[ri]:>10.0f}{rh:>+7.0f}"
              f"{d['p2_rmr'][ri]:>8.0f}{rp:>+7.0f}{prim['at'][ri]:>9.0f}{kal_at[ri]:>11.0f}")
    res_h, res_p = np.array(res_h), np.array(res_p)
    rd = dates[d["r_idx"]]
    for name, mk in (("lab 2011-16", rd < np.datetime64("2020-01-01")),
                     ("fitmate pre-drug", (rd >= np.datetime64("2020-01-01")) & (rd < np.datetime64("2024-09-17"))),
                     ("on-drug (2024-09-17..)", rd >= np.datetime64("2024-09-17"))):
        print(f"  mean residual {name:<24} n={mk.sum()}: vs Hall(F,L) {res_h[mk].mean():+.0f}, vs P2 fit "
              f"{res_p[mk].mean():+.0f}; Hall AT {prim['at'][d['r_idx']][mk].mean():+.0f}, "
              f"Kal−struct {np.nanmean(kal_at[d['r_idx']][mk]):+.0f}")

    # ── 4. on-drug plateau: intake fact or expenditure fact? ──
    print(f"\n=== On-drug plateau {PLATEAU[0]}..{PLATEAU[1]} (observed FM ≈ 60-65 lb, weight ≈ 203-206) ===")
    mk_w = era_mask(d, *PLATEAU)[d["w_idx"]]
    wi, wo = d["w_idx"][mk_w], d["w_obs"][mk_w]
    i_s = d["idx"][pd.Timestamp("2025-10-17")]; i_e = n - 1
    s_fm, s_lm = d["s_fm"][d["s_idx"] == i_s][0], d["s_lm"][d["s_idx"] == i_s][0]
    restart = simulate(pal_p, beta_p, M_CANON, d, i0=i_s, F_init=s_fm * KG, L_init=s_lm * KG,
                       at_init=prim["at"][i_s])
    restart_noat = simulate(pal_p, 0.0, M_CANON, d, i0=i_s, F_init=s_fm * KG, L_init=s_lm * KG)
    mk_p = era_mask(d, *PLATEAU)
    ee_state = (bb["K"] + GAM_F * d["p2_fm"] * KG + GAM_L * d["p2_ffm"] * KG + bb["delta"] * kw * KG
                + B_TEF * M_CANON * d["cal"] + prim["at"])          # Hall EE at the *observed* state
    print(f"  observed weigh-ins: n={len(wo)}, mean {wo.mean():.1f}, slope {slope_per_30d(wi, wo):+.2f} lb/30d; "
          f"scans 2025-10-17 FM {s_fm:.1f} -> 2026-02-28 FM {d['s_fm'][-1]:.1f} (weight {d['s_fm'][-1] + d['s_lm'][-1]:.1f})")
    print(f"  logged intake x1.09 in window: {M_CANON * d['cal'][mk_p].mean():.0f} kcal/d; "
          f"Hall-implied intake (inversion, fitted PAL): {np.nanmean(ei_fit[mk_p]):.0f}; "
          f"gap {np.nanmean(ei_fit[mk_p]) - M_CANON * d['cal'][mk_p].mean():+.0f} kcal/d")
    print(f"  Hall EE evaluated at the observed state (scan composition, Kalman weight, logged EI): "
          f"{np.nanmean(ee_state[mk_p]):.0f} kcal/d -> model-predicted balance {M_CANON * d['cal'][mk_p].mean() - np.nanmean(ee_state[mk_p]):+.0f} "
          f"kcal/d ≈ {(M_CANON * d['cal'][mk_p].mean() - np.nanmean(ee_state[mk_p])) * 30 / RHO_F / KG:+.1f} lb fat/30d")
    for name, s in (("full 15-y forward run", prim), ("restart at 2025-10-17 scan", restart),
                    ("restart, AT off", restart_noat)):
        print(f"  {name:<28}: weight slope {slope_per_30d(wi, s['bw'][wi]):+.2f} lb/30d, "
              f"mean weight {np.nanmean(s['bw'][wi]):.1f}, end weight {s['bw'][i_e]:.1f}, end FM {s['fm'][i_e]:.1f}, "
              f"EE {np.nanmean(s['ee'][mk_p]):.0f}, AT {np.nanmean(s['at'][mk_p]):+.0f}; "
              f"Kalman TDEE {np.nanmean(d['k_tdee'][mk_p]):.0f}, Kalman end FM {d['k_fm'][i_e]:.1f}")

    # ── artifact ──
    out = pd.DataFrame({"date": pd.DatetimeIndex(dates).strftime("%Y-%m-%d"),
                        "logged_cal": d["cal"], "ei_x109": M_CANON * d["cal"],
                        "hall_bw_lb": prim["bw"], "hall_fm_lb": prim["fm"], "hall_ffm_lb": prim["ffm"],
                        "hall_ee": prim["ee"], "hall_rmr": prim["rmr"], "hall_at": prim["at"],
                        "hall_pa": prim["pa"], "hall_glycogen_kg": prim["gly"], "hall_ecw_kg": prim["ecw"],
                        "implied_ei_pal15": ei_def, "implied_ei_palfit": ei_fit,
                        "kalman_tdee": d["k_tdee"], "kalman_minus_struct": kal_at}).round(2)
    path = ROOT / "analysis/DL_hall_daily.csv"
    out.to_csv(path, index=False)
    print(f"\nWrote {len(out)} rows to {path}")


if __name__ == "__main__":
    main()
