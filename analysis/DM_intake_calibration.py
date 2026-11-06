"""DM. Calibrate the intake correction (form + parameters) and the filter's noise
hyperparameters under one stated criterion each, and write both for everyone.

Why the criterion matters. Cumulative energy-balance closure cannot see the
intake scale (free-drift TDEE absorbs any level change). The one-step weight
innovation likelihood can, but it rewards matching DAY-TO-DAY weight swings,
which are water and gut contents as much as energy: scored that way the
multiplier runs away upward against Q_BUF (m ↑, Q_BUF ↓, likelihood ↑ without
an interior optimum). The k-day-ahead prediction error of the filter — TDEE
held at its filtered value, FM/LM driven by the next k days of intake, no
intermediate updates — identifies the scale where the transient has washed out
but TDEE has not had time to drift. Scanning k shows exactly that: the optimum
is ~1.4 at one week and settles to ~1.1-1.2 from two weeks on. The criterion
used here is therefore the summed k-day-ahead Gaussian log-likelihood over
k ∈ {14, 21, 28}. Windows overlap, so intervals come from a block bootstrap of
the per-start-day likelihood terms (90-day blocks), not from the raw curvature.

Two stages, iterated once:
  A. noise hyperparameters (σ_TDEE, Q_BUF, R_WEIGHT, BUF_HL, BETA_BUF) by the
     one-step prediction-error likelihood over ALL assimilated observations
     (weigh-ins + FM + LM scans) at the current intake model;
  B. the intake model by the multi-week criterion at those hyperparameters.

Candidate forms for consumed = f(logged):
  multiplicative  m · logged                                      (1 param)
  affine          a + m · logged                                  (2)  a absorbed by TDEE → expected flat
  level           logged · (m + m1·(logged − 2000)/1000)          (2)
  macro           m·logged + a_fiber·fiber_g + a_protein·protein_g (3)
  era             m per calendar era (loss ≤2016, regain 2017-21, plateau 2022-24, drug)  (4)
  direction       m by sign of the 90-d scan-interpolated FM slope (losing/stable/gaining) (3)
  + per-year m profile (diagnostic)

Selection: a richer form replaces the uniform multiplier only if its raw ΔAIC
is ≥ 4 AND the block bootstrap puts P(ΔLL > 0) ≥ 0.95. Output is written to
analysis/results/DM_intake_calibration.json and served by analysis/intake_model.py.

Caveat on meaning: what the weight response identifies is the NET scale — the
energy that actually moves fat and lean per logged calorie, after thermic effect
and any within-window expenditure response — which is what the filter's energy
equation needs. It is a lower bound on the gross diary undercount.

Command: python analysis/DM_intake_calibration.py
Artifact: analysis/results/DM_intake_calibration.json, analysis/DM_profiles.csv
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))
import P4_kalman_filter as P4                      # noqa: E402
from intake_model import apply_model, MODEL_PARAM_NAMES, ERA_BOUNDS, era_of, _fm_direction  # noqa: E402
from results import save                           # noqa: E402

CHI2_95 = 1.92
MAX_DTDEE = 100                      # hard physiological cap on daily TDEE change (kept)
HORIZONS = (14, 21, 28)              # days; the multi-week criterion
DIAG_HORIZONS = (7, 14, 21, 28, 42)  # for the horizon scan printed as a diagnostic
DELTA_AIC_TO_SWITCH, P_BOOT_TO_SWITCH = 4.0, 0.95
BOOT_BLOCK, N_BOOT, SEED = 90, 1000, 0
M_GRID = np.round(np.arange(0.90, 1.601, 0.05), 3)
HP_NAMES = ["sigma_tdee", "q_buf", "r_weight", "buf_hl", "beta_buf"]
HP_LOG = np.array([True, True, True, True, False])
HP_LO = np.array([np.log(5.0), np.log(0.005), np.log(0.2), np.log(1.0), -0.5])
HP_HI = np.array([np.log(150.0), np.log(2.0), np.log(3.0), np.log(30.0), 3.0])
HP_CANON = dict(sigma_tdee=50.0, q_buf=P4.Q_BUF, r_weight=P4.R_WEIGHT, buf_hl=P4.BUF_HL, beta_buf=0.0)  # hand-set, for comparison


def gauss_terms(innov, var):
    ok = ~np.isnan(innov) & ~np.isnan(var)
    out = np.zeros(len(innov))
    out[ok] = -0.5 * (np.log(2 * np.pi * var[ok]) + innov[ok] ** 2 / var[ok])
    return out, ok


class Scorer:
    def __init__(self):
        smooth, intake, strength, comp, rmr, p2, temp = P4.load_data()
        self.intake, self.rmr = intake, rmr
        (self.dates, _, self.wobs, self.fm_scan, self.lm_scan,
         self.erm, self.wpw) = P4.build_inputs(smooth, intake, strength, comp, p2,
                                               calories_true=intake["calories"].values)
        self.fs = int(np.where(~np.isnan(self.fm_scan))[0][0])
        self.n_w = int(np.sum(~np.isnan(self.wobs)))
        self.n_fm = int(np.sum(~np.isnan(self.fm_scan)))
        self.n_lm = int(np.sum(~np.isnan(self.lm_scan)))
        self.years = pd.DatetimeIndex(self.dates).year
        self.block = (np.arange(len(self.dates)) // BOOT_BLOCK)
        self.n_eval = 0

    def forward(self, cal, hp):
        coll = {}
        out = P4.kalman_v4_forward(cal, self.wobs, self.fm_scan, self.lm_scan, self.erm, self.wpw,
                                   self.fm_scan[self.fs], self.lm_scan[self.fs],
                                   hp["sigma_tdee"], MAX_DTDEE, collect=coll,
                                   q_buf=hp["q_buf"], r_weight=hp["r_weight"],
                                   beta_buf=hp["beta_buf"], buf_hl=hp["buf_hl"])
        self.n_eval += 1
        return out, coll

    # one-step, all observations (stage A)
    def ll1(self, cal, hp, parts=False):
        (xf, Pf, xp, Pp, innov, _), c = self.forward(cal, hp)
        tw, _ = gauss_terms(innov, c["innov_var"])
        tf, _ = gauss_terms(c["fm_innov"], c["fm_var"])
        tl, _ = gauss_terms(c["lm_innov"], c["lm_var"])
        if parts:
            return float(tw.sum() + tf.sum() + tl.sum()), dict(weight=float(tw.sum()), fm=float(tf.sum()), lm=float(tl.sum()))
        return float(tw.sum() + tf.sum() + tl.sum())

    # multi-week prediction error (stage B); returns total and per-start-day terms
    def llk(self, cal, hp, horizons=HORIZONS, terms=False):
        (xf, Pf, xp, Pp, innov, _), c = self.forward(cal, hp)
        tot = np.zeros(len(cal)); n_by_k = {}
        for k in horizons:
            ik, vk = P4.kstep_prediction_errors(xf, Pf, cal, self.wobs, self.wpw, k, hp["sigma_tdee"],
                                                q_buf=hp["q_buf"], r_weight=hp["r_weight"],
                                                beta_buf=hp["beta_buf"], buf_hl=hp["buf_hl"])
            t, ok = gauss_terms(ik, vk); tot += t; n_by_k[k] = int(ok.sum())
        return (float(tot.sum()), tot, n_by_k) if terms else float(tot.sum())

    def fit_quality(self, cal, hp):
        (xf, Pf, xp, Pp, innov, _), c = self.forward(cal, hp)
        xs, _ = P4.rts_smoother_v4(xf, Pf, xp, Pp, self.fs, self.wpw)
        d = P4.diagnostics(self.dates, xs, innov, self.fm_scan, self.lm_scan, cal, self.rmr, self.wobs)
        return dict(fm_scan_mae=d["fm_scan_mae"], lm_scan_mae=d["lm_scan_mae"], innov_acf1=d["innov_acf1"],
                    innov_std=d["innov_std"], tdee_rmr_mean=d["tdee_rmr_mean"], tdee_mean=d["tdee_mean"],
                    eb_residual_lb=d["eb_residual_lb"])

    def bootstrap_argmax(self, grid, term_rows, rng):
        """term_rows: list (per grid value) of per-day term arrays. Block-bootstrap
        the argmax over the grid; returns (m_hat samples)."""
        blocks = np.unique(self.block)
        sums = np.array([[t[self.block == b].sum() for b in blocks] for t in term_rows])   # grid × blocks
        out = []
        for _ in range(N_BOOT):
            pick = rng.integers(0, len(blocks), len(blocks))
            out.append(grid[int(np.argmax(sums[:, pick].sum(axis=1)))])
        return np.array(out)

    def bootstrap_delta(self, t_a, t_b, rng):
        """P(Σ(t_a − t_b) > 0) under the block bootstrap."""
        blocks = np.unique(self.block)
        d = np.array([(t_a - t_b)[self.block == b].sum() for b in blocks])
        return float(np.mean([d[rng.integers(0, len(d), len(d))].sum() > 0 for _ in range(N_BOOT)]))


def hp_from_z(z):
    z = np.clip(np.asarray(z, dtype=float), HP_LO, HP_HI)
    return dict(zip(HP_NAMES, np.where(HP_LOG, np.exp(z), z)))


def hp_to_z(hp):
    v = np.array([hp[k] for k in HP_NAMES], dtype=float)
    return np.where(HP_LOG, np.log(np.maximum(v, 1e-9)), v)


def fit_hp(sc, cal, hp0):
    r = minimize(lambda z: -sc.ll1(cal, hp_from_z(z)), hp_to_z(hp0), method="Nelder-Mead",
                 options=dict(xatol=3e-3, fatol=0.05, maxiter=600, initial_simplex=None))
    return hp_from_z(r.x), -float(r.fun)


def fit_form(sc, model, names, x0, bounds, hp):
    def neg(x):
        for (lo, hi), v in zip(bounds, x):
            if v < lo or v > hi:
                return 1e9
        return -sc.llk(apply_model(model, dict(zip(names, x)), sc.intake), hp)
    r = minimize(neg, np.array(x0, dtype=float), method="Nelder-Mead",
                 options=dict(xatol=2e-3, fatol=0.05, maxiter=120 + 80 * len(x0)))
    return dict(zip(names, map(float, r.x))), -float(r.fun)


def profile_cond(sc, model, params, hp, name, grid):
    rows = [sc.llk(apply_model(model, {**params, name: float(v)}, sc.intake), hp, terms=True)[1] for v in grid]
    lls = np.array([r.sum() for r in rows])
    inside = grid[lls >= lls.max() - CHI2_95]
    return lls, rows, float(inside.min()), float(inside.max())


def main():
    sc = Scorer(); rng = np.random.default_rng(SEED)
    profiles = []
    logged = sc.intake["calories"].values
    print("=== DM. Intake correction + filter hyperparameters ===")
    print(f"observations: {sc.n_w} weigh-ins, {sc.n_fm} FM scans, {sc.n_lm} LM scans; max_dTDEE = {MAX_DTDEE} (fixed)")
    print(f"criterion for the intake model: summed {HORIZONS}-day-ahead weight prediction likelihood; "
          f"for the noise hyperparameters: one-step likelihood over all observations\n")

    # ── 0. horizon scan at the hand-set tuning: why the criterion is multi-week ──
    print("Horizon scan, uniform m, hand-set tuning (optimum m per horizon; one-step shown for contrast):")
    scan = {}
    for k in DIAG_HORIZONS:
        lls = np.array([sc.llk(m * logged, HP_CANON, horizons=(k,)) for m in M_GRID])
        scan[k] = float(M_GRID[int(np.argmax(lls))])
        for m, l in zip(M_GRID, lls):
            profiles.append(dict(model=f"horizon_{k}", param="m", value=m, ll=l))
    l1 = np.array([sc.ll1(m * logged, HP_CANON) for m in M_GRID])
    scan["one_step"] = float(M_GRID[int(np.argmax(l1))])
    for m, l in zip(M_GRID, l1):
        profiles.append(dict(model="one_step_fixed_tuning", param="m", value=m, ll=l))
    print("  " + "  ".join(f"k={k}: m={v:.2f}" for k, v in scan.items()))

    # ── A/B iteration ────────────────────────────────────────────────────
    hp = dict(HP_CANON); m_hat = scan[HORIZONS[0]]
    for it in range(2):
        hp, ll1_hp = fit_hp(sc, m_hat * logged, hp)
        lls = np.array([sc.llk(m * logged, hp) for m in M_GRID]); m_hat = float(M_GRID[int(np.argmax(lls))])
        print(f"  iteration {it + 1}: hyperparameters σ_TDEE {hp['sigma_tdee']:.1f}, Q_BUF {hp['q_buf']:.3f}, "
              f"R_WEIGHT {hp['r_weight']:.2f}, BUF_HL {hp['buf_hl']:.1f}, BETA_BUF {hp['beta_buf']:.3f} "
              f"(one-step LL {ll1_hp:.1f}); uniform m → {m_hat:.2f}")
    ll1_canon = sc.ll1(m_hat * logged, HP_CANON)
    ll1_nob = sc.ll1(m_hat * logged, {**hp, "beta_buf": 0.0})
    print(f"  one-step LL at the hand-set tuning {ll1_canon:.1f} vs fitted {ll1_hp:.1f}; fitted with BETA_BUF=0: {ll1_nob:.1f}")

    # ── 1. uniform multiplier: profile + bootstrap ───────────────────────
    print("\nUniform multiplier (multi-week criterion, fitted hyperparameters):")
    rows = [sc.llk(m * logged, hp, terms=True) for m in M_GRID]
    lls = np.array([r[0] for r in rows]); n_by_k = rows[0][2]
    for m, l in zip(M_GRID, lls):
        profiles.append(dict(model="multiplicative", param="m", value=m, ll=l))
    b = int(np.argmax(lls)); inside = M_GRID[lls >= lls.max() - CHI2_95]
    p_u, ll_u = fit_form(sc, "multiplicative", ["m"], [M_GRID[b]], [(0.8, 1.8)], hp)
    boot = sc.bootstrap_argmax(M_GRID, [r[1] for r in rows], rng)
    ci_boot = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]
    print(f"  {'m':>6}{'LL':>10}{'ΔLL':>8}")
    for m, l in zip(M_GRID, lls):
        print(f"  {m:>6.2f}{l:>10.1f}{l - lls.max():>+8.2f}")
    print(f"  → m = {p_u['m']:.3f};  raw-curvature 95% [{inside.min():.2f}, {inside.max():.2f}] (overstated: windows overlap);"
          f"  block-bootstrap 95% [{ci_boot[0]:.2f}, {ci_boot[1]:.2f}]  (windows per horizon {n_by_k})")
    base = {"model": "multiplicative", "params": p_u, "ll": ll_u, "k": 1,
            "ci": {"m": ci_boot}, "ci_raw": {"m": [float(inside.min()), float(inside.max())]}}
    base_terms = sc.llk(apply_model("multiplicative", p_u, sc.intake), hp, terms=True)[1]
    fits = [base]

    def report(model, names, p, ll, ci, k, extra=""):
        terms = sc.llk(apply_model(model, p, sc.intake), hp, terms=True)[1]
        pboot = sc.bootstrap_delta(terms, base_terms, rng)
        fits.append({"model": model, "params": p, "ll": ll, "k": k, "ci": ci, "p_boot_beats_uniform": pboot})
        txt = ", ".join(f"{n}={p[n]:.3f}" + (f" [{ci[n][0]:.2f}, {ci[n][1]:.2f}]" if n in ci else "") for n in names)
        print(f"{model:<15} {txt}  LL {ll:.1f} (ΔLL vs uniform {ll - ll_u:+.1f}; bootstrap P(better) {pboot:.2f}){extra}")

    print("\nOther forms (same criterion and hyperparameters; CIs are block-bootstrap, conditional on the other parameters):")

    def boot_ci(model, p, name, grid):
        lls_, rows_, lo, hi = profile_cond(sc, model, p, hp, name, grid)
        for v, l in zip(grid, lls_):
            profiles.append(dict(model=model, param=name, value=v, ll=l))
        bs = sc.bootstrap_argmax(grid, rows_, rng)
        return [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], lls_

    # affine
    p, ll = fit_form(sc, "affine", ["a", "m"], [0.0, p_u["m"]], [(-600, 600), (0.8, 1.8)], hp)
    g = np.arange(-400, 401, 100.0); ci_a, pl = boot_ci("affine", p, "a", g)
    report("affine", ["a", "m"], p, ll, {"a": ci_a}, 2,
           f"\n                LL range over a∈[-400,400]: {pl.max() - pl.min():.2f} → "
           f"{'flat: a constant offset is absorbed by TDEE' if pl.max() - pl.min() < CHI2_95 else 'identified'}")
    # level
    p, ll = fit_form(sc, "level", ["m", "m1"], [p_u["m"], 0.0], [(0.8, 1.8), (-0.6, 0.6)], hp)
    g = np.round(np.arange(-0.4, 0.41, 0.1), 2); ci_m1, _ = boot_ci("level", p, "m1", g)
    report("level", ["m", "m1"], p, ll, {"m1": ci_m1}, 2,
           f"\n                m at 1200 cal {p['m'] + p['m1'] * -0.8:.2f}, at 2000 {p['m']:.2f}, at 3500 {p['m'] + p['m1'] * 1.5:.2f}")
    # macro
    p, ll = fit_form(sc, "macro", ["m", "a_fiber", "a_protein"], [p_u["m"], 0.0, 0.0], [(0.8, 1.8), (-6, 6), (-3, 3)], hp)
    ci_f, _ = boot_ci("macro", p, "a_fiber", np.arange(-4.0, 4.01, 1.0))
    ci_p, _ = boot_ci("macro", p, "a_protein", np.arange(-2.0, 2.01, 0.5))
    report("macro", ["m", "a_fiber", "a_protein"], p, ll, {"a_fiber": ci_f, "a_protein": ci_p}, 3,
           "\n                (fiber −2 = labels count 4 cal/g, ~2 metabolizable; protein −1 = 25% TEF charged to intake)")
    # era
    names = MODEL_PARAM_NAMES["era"]
    p, ll = fit_form(sc, "era", names, [p_u["m"]] * 4, [(0.8, 1.8)] * 4, hp)
    ci = {nm: boot_ci("era", p, nm, M_GRID)[0] for nm in names}
    report("era", names, p, ll, ci, 4)
    eras = era_of(sc.intake["date"].values)
    for nm, (e, lo_d, hi_d) in zip(names, ERA_BOUNDS):
        print(f"                {e:<8} {lo_d[:4]}-{hi_d[:4]}  weigh-ins {int(np.sum(~np.isnan(sc.wobs) & (eras == e)))}")
    print(f"                regain/loss ratio {p['m_regain'] / p['m_loss']:.3f}")
    era_fit = dict(p)
    # direction
    names = MODEL_PARAM_NAMES["direction"]
    p, ll = fit_form(sc, "direction", names, [p_u["m"]] * 3, [(0.8, 1.8)] * 3, hp)
    ci = {nm: boot_ci("direction", p, nm, M_GRID)[0] for nm in names}
    report("direction", names, p, ll, ci, 3, "  (90-d slope of scan-interpolated FM, |slope| > 0.01 lb/d)")
    dirn = _fm_direction(sc.intake["date"].values)
    for nm in names:
        print(f"                {nm[2:]:<8} weigh-ins {int(np.sum(~np.isnan(sc.wobs) & (dirn == nm[2:])))}")
    print(f"                gaining/losing ratio {p['m_gaining'] / p['m_losing']:.3f}")

    # ── per-year profile ─────────────────────────────────────────────────
    print("\nPer-year multiplier (that year's m varied, other years at the uniform optimum; conditional, raw-curvature CI):")
    print(f"  {'year':<6}{'m':>7}{'95% CI':>16}{'weigh-ins':>11}")
    year_rows = []
    cal_base = p_u["m"] * logged
    for y in sorted(set(sc.years)):
        n_w = int(np.sum(~np.isnan(sc.wobs) & (sc.years == y)))
        if n_w < 20:
            continue
        lls_y = np.array([sc.llk(np.where(sc.years == y, m * logged, cal_base), hp) for m in M_GRID])
        inside = M_GRID[lls_y >= lls_y.max() - CHI2_95]
        my = float(M_GRID[int(np.argmax(lls_y))])
        edge = "" if 0 < int(np.argmax(lls_y)) < len(M_GRID) - 1 else "  (boundary)"
        year_rows.append(dict(year=int(y), m=my, lo=float(inside.min()), hi=float(inside.max()), weigh_ins=n_w))
        print(f"  {y:<6}{my:>7.2f}   [{inside.min():.2f}, {inside.max():.2f}]{n_w:>11}{edge}")

    # ── table, selection, write ──────────────────────────────────────────
    for f in fits:
        f["aic"] = 2 * f["k"] - 2 * f["ll"]
    for f in fits:
        f["d_aic_vs_uniform"] = f["aic"] - base["aic"]
        f.update(sc.fit_quality(apply_model(f["model"], f["params"], sc.intake), hp))
    print(f"\n{'model':<15}{'k':>3}{'LL':>10}{'ΔAIC':>8}{'P(boot)':>8}{'FM MAE':>8}{'LM MAE':>8}{'acf1':>7}{'TDEE/RMR':>10}")
    for f in fits:
        print(f"{f['model']:<15}{f['k']:>3}{f['ll']:>10.1f}{f['d_aic_vs_uniform']:>+8.1f}{f.get('p_boot_beats_uniform', float('nan')):>8.2f}"
              f"{f['fm_scan_mae']:>8.3f}{f['lm_scan_mae']:>8.3f}{f['innov_acf1']:>7.3f}{f['tdee_rmr_mean']:>10.3f}")
    cands = [f for f in fits if f is not base and f["d_aic_vs_uniform"] <= -DELTA_AIC_TO_SWITCH
             and f.get("p_boot_beats_uniform", 0) >= P_BOOT_TO_SWITCH]
    winner = min(cands, key=lambda f: f["aic"]) if cands else base
    cal_w = apply_model(winner["model"], winner["params"], sc.intake)
    eff = float(cal_w.sum() / logged.sum())
    print(f"\nSelected: {winner['model']} (rule: richer form needs ΔAIC ≤ −{DELTA_AIC_TO_SWITCH:.0f} and bootstrap P ≥ {P_BOOT_TO_SWITCH})")
    print(f"  intake params: " + ", ".join(f"{k}={v:.4g}" for k, v in winner["params"].items()))
    print(f"  filter hyperparameters: σ_TDEE {hp['sigma_tdee']:.1f} cal/d, Q_BUF {hp['q_buf']:.3f} lb², R_WEIGHT {hp['r_weight']:.2f} lb², "
          f"BUF_HL {hp['buf_hl']:.1f} d, BETA_BUF {hp['beta_buf']:.3f} lb/1000 cal, max_dTDEE {MAX_DTDEE}")
    print(f"  effective record-wide multiplier Σtrue/Σlogged = {eff:.3f}")
    print(f"  filter evaluations: {sc.n_eval}")

    save("DM_intake_calibration",
         model=winner["model"], params=winner["params"], ci=winner["ci"], ll=winner["ll"],
         criterion=dict(intake_model=f"summed {HORIZONS}-day-ahead weight prediction log-likelihood, block bootstrap {BOOT_BLOCK} d",
                        hyperparameters="one-step prediction-error log-likelihood over weigh-ins + FM + LM scans"),
         filter=dict(sigma_tdee=hp["sigma_tdee"], q_buf=hp["q_buf"], r_weight=hp["r_weight"],
                     buf_hl=hp["buf_hl"], beta_buf=hp["beta_buf"], max_dtdee=MAX_DTDEE),
         one_step_ll=dict(fitted=ll1_hp, hand_set=ll1_canon, fitted_without_beta_buf=ll1_nob),
         horizon_scan_m=scan, effective_multiplier=eff,
         selection_rule=f"richer form only if dAIC <= -{DELTA_AIC_TO_SWITCH} and bootstrap P(better) >= {P_BOOT_TO_SWITCH}",
         candidates=fits, uniform_m=p_u["m"], uniform_ci=ci_boot, uniform_ci_raw=base["ci_raw"]["m"],
         era_fit=era_fit, per_year=year_rows,
         n_obs=dict(weigh_ins=sc.n_w, fm_scans=sc.n_fm, lm_scans=sc.n_lm, windows_per_horizon=n_by_k))
    pd.DataFrame(profiles).to_csv(ROOT / "analysis" / "DM_profiles.csv", index=False)
    print("Wrote analysis/results/DM_intake_calibration.json and analysis/DM_profiles.csv")


if __name__ == "__main__":
    main()
