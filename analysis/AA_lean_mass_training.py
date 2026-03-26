"""Test: does strength training add lean mass that decays without training?

Model: each workout adds Δ grams of lean mass. Without training, the
accumulated training-effect decays exponentially with half-life τ.

    training_effect(t) = Δ × Σ exp(-ln2/τ × days_since_workout_i)

Lean mass = baseline(fat_mass, weight) + training_effect.

The baseline accounts for the fact that lean mass varies with total weight
(more weight = more supporting tissue). The training effect is what's left.

Fit Δ and τ from 70 composition measurements × 393 workout sessions.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent.parent


def load():
    comp = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    workouts = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    return comp, workouts


def training_effect(comp_dates, workout_dates, delta_lbs, half_life_days):
    """Compute accumulated training effect at each composition date."""
    decay = np.log(2) / half_life_days
    effects = []
    workout_days = workout_dates.values.astype("datetime64[D]")
    for comp_date in comp_dates:
        comp_day = np.datetime64(comp_date, "D")
        days_since = (comp_day - workout_days).astype(float)
        # Only count past workouts
        past = days_since[days_since > 0]
        effect = delta_lbs * np.sum(np.exp(-decay * past))
        effects.append(effect)
    return np.array(effects)


def fit_model(comp, workouts):
    """Fit delta and half-life to minimize lean mass residuals."""
    # Baseline: lean mass depends on weight (heavier = more lean tissue)
    # We fit: lean = a * weight + b + training_effect(delta, half_life)

    lean = comp["lean_mass_lbs"].values
    weight = comp["weight_lbs"].values
    comp_dates = comp["date"]
    workout_dates = workouts["date"]

    def objective(params):
        delta, half_life, a, b = params
        if half_life < 7 or half_life > 365 or delta < 0 or delta > 2:
            return 1e6
        te = training_effect(comp_dates, workout_dates, delta, half_life)
        predicted = a * weight + b + te
        return np.sum((lean - predicted) ** 2)

    # Grid search for starting point
    best_cost = np.inf
    best_params = None
    for delta in [0.05, 0.1, 0.2, 0.5]:
        for half_life in [30, 60, 90, 120, 180]:
            for a in [0.3, 0.4, 0.5]:
                b = np.mean(lean) - a * np.mean(weight)
                cost = objective([delta, half_life, a, b])
                if cost < best_cost:
                    best_cost = cost
                    best_params = [delta, half_life, a, b]

    result = minimize(objective, best_params, method="Nelder-Mead",
                      options={"maxiter": 10000, "xatol": 0.001})
    return result.x, result


def fit_baseline_only(weight, lean):
    """Best linear baseline without any training term."""
    X = np.column_stack([weight, np.ones(len(weight))])
    a, b = np.linalg.lstsq(X, lean, rcond=None)[0]
    return a, b


def main():
    comp, workouts = load()
    lean = comp["lean_mass_lbs"].values
    weight = comp["weight_lbs"].values

    # Simple correlations first
    print("=== Simple correlations: trailing workouts vs lean mass ===")
    for window in [30, 60, 90]:
        counts = []
        for _, row in comp.iterrows():
            n = len(workouts[(workouts["date"] >= row["date"] - pd.Timedelta(days=window)) &
                             (workouts["date"] < row["date"])])
            counts.append(n)
        r = np.corrcoef(counts, comp["lean_mass_lbs"])[0, 1]
        print(f"  {window:3d}-day trailing workouts vs lean mass: r={r:.4f}")

    # Partial: controlling for weight
    print("\n=== Partial (controlling for weight) ===")
    from numpy.linalg import lstsq
    for window in [30, 60, 90]:
        counts = []
        for _, row in comp.iterrows():
            n = len(workouts[(workouts["date"] >= row["date"] - pd.Timedelta(days=window)) &
                             (workouts["date"] < row["date"])])
            counts.append(n)
        counts = np.array(counts, dtype=float)
        X = np.column_stack([comp["weight_lbs"].values, np.ones(len(comp))])
        res_c = counts - X @ lstsq(X, counts, rcond=None)[0]
        res_l = comp["lean_mass_lbs"].values - X @ lstsq(X, comp["lean_mass_lbs"].values, rcond=None)[0]
        r = np.corrcoef(res_c, res_l)[0, 1]
        print(f"  {window:3d}-day trailing workouts vs lean mass (partial): r={r:.4f}")

    # Fit the decay model
    print("\n=== Decay model: lean = a×weight + b + Δ×Σexp(-ln2/τ × days) ===")
    params, result = fit_model(comp, workouts)
    delta, half_life, a, b = params

    print(f"  Δ (lean mass per workout): {delta:.3f} lbs ({delta * 453.6:.0f} g)")
    print(f"  Half-life: {half_life:.0f} days")
    print(f"  Baseline: lean = {a:.3f} × weight + {b:.1f}")
    print(f"  Optimizer converged: {result.success}")

    # Validation: predicted vs actual
    te = training_effect(comp["date"], workouts["date"], delta, half_life)
    predicted = a * weight + b + te
    residuals = lean - predicted
    rmse = np.sqrt(np.mean(residuals ** 2))

    # Compare to baseline-only model (no training effect)
    a0, b0 = fit_baseline_only(weight, lean)
    pred_baseline = a0 * weight + b0
    rmse_baseline = np.sqrt(np.mean((lean - pred_baseline) ** 2))

    print(f"\n  RMSE with training effect: {rmse:.2f} lbs")
    print(f"  RMSE baseline only (weight): {rmse_baseline:.2f} lbs")
    print(f"  Improvement: {(1 - rmse / rmse_baseline) * 100:.1f}%")
    print(f"  Baseline-only fit: lean = {a0:.3f} × weight + {b0:.1f}")

    # Show the effect at different training densities
    print(f"\n=== Steady-state training effect ===")
    for sessions_per_week in [0, 1, 2, 3]:
        # At steady state with weekly sessions, effect = delta * sum(exp(-decay * 7*k) for k=0..inf)
        # = delta / (1 - exp(-decay * 7))
        if sessions_per_week == 0:
            ss = 0
        else:
            decay = np.log(2) / half_life
            interval = 7 / sessions_per_week
            ss = delta / (1 - np.exp(-decay * interval))
        print(f"  {sessions_per_week}×/week: +{ss:.1f} lbs ({ss * 453.6:.0f} g) lean mass")

    # Build and save per-measurement fit
    rows = []
    for i, (_, row) in enumerate(comp.iterrows()):
        w30 = len(workouts[(workouts["date"] >= row["date"] - pd.Timedelta(days=30)) &
                           (workouts["date"] < row["date"])])
        rows.append({
            "date": row["date"].strftime("%Y-%m-%d"),
            "lean_mass_lbs": row["lean_mass_lbs"],
            "predicted_lbs": round(predicted[i], 1),
            "training_effect_lbs": round(te[i], 1),
            "workouts_30d": w30,
            "residual_lbs": round(residuals[i], 1),
        })

    out = pd.DataFrame(rows)
    out_path = ROOT / "analysis" / "AA_lean_mass_training.csv"
    out.to_csv(out_path, index=False)

    print(f"\n=== Per-measurement: actual vs predicted lean mass ===")
    print(f"{'Date':>12} {'Actual':>7} {'Pred':>7} {'Train':>6} {'W30':>4} {'Err':>6}")
    for _, r in out.iterrows():
        print(f"{r['date']:>12} {r['lean_mass_lbs']:7.1f} "
              f"{r['predicted_lbs']:7.1f} {r['training_effect_lbs']:6.1f} "
              f"{r['workouts_30d']:4d} {r['residual_lbs']:+6.1f}")

    print(f"\nWrote {len(out)} rows to {out_path}")


if __name__ == "__main__":
    main()
