# New findings

This folder is the living, audit-first analysis of the weight-loss project. It separates direct observations from derived-model claims and external evidence. The legacy `README.md`, `FINDINGS.md`, and `RECOMMENDATIONS.md` remain useful records of hypotheses, but their central set-point and drug-schedule conclusions are not currently decision-grade.

## Reports

- [`00_initial_integrity_and_safety_audit.md`](00_initial_integrity_and_safety_audit.md) — first-pass data audit, model audit, immediate safety boundary, and practical working strategy.
- [`01_claim_registry.md`](01_claim_registry.md) — claim-by-claim evidence status.
- [`02_obesity_hypothesis_evidence_map.md`](02_obesity_hypothesis_evidence_map.md) — mainstream and speculative obesity theories, current external evidence, and what this repository can actually test.
- [`03_direct_behavioral_signals.md`](03_direct_behavioral_signals.md) — held-out menu novelty/repetition results and a 41-trip paired event study.
- [`04_tirzepatide_observed_outcomes_and_support.md`](04_tirzepatide_observed_outcomes_and_support.md) — raw treatment-era intake/weight, the weekday-confounded cycle, same-device composition, protein, fiber, and strength.
- [`05_genetics_and_fatty_acids_audit.md`](05_genetics_and_fatty_acids_audit.md) — correction of common-SNP overclaims, seed-oil/FADS interpretation, supplements, and useful follow-up.
- [`06_current_menu_leverage.md`](06_current_menu_leverage.md) — the final-90-day foods, calorie-neutral protein substitutions, and a concrete two-stage menu plan.
- [`07_next_data_and_preregistered_experiments.md`](07_next_data_and_preregistered_experiments.md) — refresh boundary, low-burden tracking templates, ABAB menu test, travel playbook, and training outcomes.
- [`08_cardiometabolic_and_clinical_monitoring_audit.md`](08_cardiometabolic_and_clinical_monitoring_audit.md) — clinical-data gaps, saturated fat, lipids/ApoB, glycemia, BP, liver/kidney/thyroid, micronutrients, and medication/supplement reconciliation.
- [`09_sleep_activity_reaudit.md`](09_sleep_activity_reaudit.md) — exact-date sleep/steps/exercise re-analysis, legacy future-window errors, source-regime validity, and prospective metrics.
- [`10_current_action_plan.md`](10_current_action_plan.md) — integrated near-term action priorities and decision rules.
- [`11_P4b_forward_model_audit.md`](11_P4b_forward_model_audit.md) — leakage-resistant opt-in state-space model, deterministic tests, first temporal results, calibration failure, and naive baselines.
- [`12_diet_epoch_direct_audit.md`](12_diet_epoch_direct_audit.md) — raw-observation pre/during/post audit of keto, potatoes, weekend fasting, low protein, Enlightened products, travel, and annotation overlap.
- [`13_temporal_semantics_audit.md`](13_temporal_semantics_audit.md) — repository-wide forward-window audit, invalidation of the AF protective-variance claim, and corrected AB/next-record semantics.
- [`14_hyperpalatable_food_audit.md`](14_hyperpalatable_food_audit.md) — published fat+sugar criterion, enrichment and food-mass limits, replicated same-day marker, and null held-out future-risk result.
- [`15_long_run_trajectory_audit.md`](15_long_run_trajectory_audit.md) — direct observed loss/regain episodes, informative weigh-in gaps, water-sensitive endpoints, and a negative internal locked warning test.
- [`16_calendar_and_extreme_clustering_audit.md`](16_calendar_and_extreme_clustering_audit.md) — calendar effects, extreme-day runs, next-day continuation, and a locked but calibration-sensitive calorie-momentum warning model.
- [`17_energy_density_and_nutrient_leverage_audit.md`](17_energy_density_and_nutrient_leverage_audit.md) — exact-mass and nutrient-field coverage, same-day coupling, locked future-prediction nulls, beverage calories, and current meal leverage.

## Reproducer

Run the descriptive repository audit with:

```bash
uv run --with pandas python analysis/EA_repository_audit.py
uv run --with pandas python analysis/EH_cardiometabolic_snapshot.py
uv run --with pandas --with numpy python analysis/EI_sleep_activity_reaudit.py
uv run --with numpy --with pandas python analysis/test_P4b_forward_state_space.py
uv run --with pandas python analysis/EJ_diet_epoch_direct_audit.py
uv run --with numpy --with pandas python analysis/EK_temporal_semantics_audit.py
uv run --with numpy --with pandas python analysis/EJ_hyperpalatable_food_audit.py
uv run --with numpy --with pandas python analysis/EL_long_run_trajectory_audit.py
uv run --with numpy --with pandas python analysis/EM_calendar_extreme_clustering.py
uv run --with numpy --with pandas python analysis/EN_energy_density_nutrient_leverage.py
```

The repository audit modifies no source or derived data. The cardiometabolic
snapshot modifies no source data and rebuilds the three `analysis/EH_*.csv`
artifacts documented in report 08.

## Current priority

P4b now provides genuinely forward-filtered estimates, separate retrospective
smooths, exported uncertainty, input isolation, and temporal validation. It is
an opt-in engineering scaffold, not a replacement decision engine: later-period
uncertainty is too broad, and its conditional fat-mass forecast does not beat a
simple weight-adjusted persistence baseline. Legacy analyses that depend on
`analysis/P4_kalman_daily.csv` remain exploratory. Any P4b v2 tuning must call
2025–2026 calibration and reserve newly refreshed observations as a new holdout.

The strongest new raw-data intervention target is travel/context disruption. A second, smaller signal suggests that a stable menu reduces extreme-intake risk, while food novelty marks it. These are descriptive and should be tested prospectively.

The named diet experiments reinforce the menu-structure interpretation. All
four potato attempts reduced intake by 524–748 kcal/day versus their immediate
pre-windows, but rapid scale loss was substantially reversible; weekend fasting
also produced a large immediate deficit followed by a sharp intake/scale
rebound. Keto and the Enlightened era lack sufficient observed-scale coverage
for weight claims. The practical target is therefore a repeatable,
high-satiety, protein-adequate menu—not a potato monodiet or recurrent 36-hour
fast.

The fixed Fazzino fat+sugar criterion also marks higher-calorie current days,
but it does not reliably predict next-day calories or extreme-intake risk after
current-intake and calendar/context controls. Full hyper-palatable-food
exposure cannot be reconstructed because exact food weight covers only about
17% of solid-food calories. This does not justify elevating a blanket
hyper-palatable-food rule above the travel and stable-menu interventions.

The long trajectory is consistent with sustained intake regimes rather than a
single universal maintenance threshold. Three clean loss-to-regain reversals
coincided with large, sustained logged-intake increases, but no one-month
intake/context score reliably ranked next-quarter weight change in the internal
locked test. All five weigh-in gaps of at least 90 days ended higher, making
measurement staleness an informative process marker rather than evidence that
weighing itself causes control.

Before treatment, extreme-intake days frequently continued into the next day,
including at home. Recent calorie momentum improved one-day-ahead ranking in
both untouched pre-treatment periods, but calibration shifted across eras.
Use this as a cue to reinstate an ordinary default meal, not as a fasting,
exercise, or medication trigger.

Whole-diet energy density cannot be reconstructed retrospectively: exact mass
covers only 17.0% of solid calories overall and 10.8% in the final 90 days.
Protein/fiber density, beverage share, snack share, item count, variety, and
meal count did not add stable future-intake prediction across both holdouts.
Breakfast and snacks remain good calorie-neutral nutrition targets because
they supply 31% of current calories while having the lowest protein and fiber
density—not because a regression proved a satiety mechanism.

The treatment itself is producing a large, durable observed benefit. The current support gap is nutritional: late-treatment protein and fiber are low even though same-device bioimpedance is reassuring about lean mass so far.

The largest non-weight information gap is the absence of auditable conventional clinical measurements. A current lipid panel/ApoB/Lp(a), HbA1c/fasting glucose, standardized blood pressure, CBC/CMP, and symptom-directed thyroid/iron/B12/vitamin D review are more decision-relevant than further latent-state tuning. Recent matched foods also contribute about 15% of energy from known saturated fat, making lipid measurement and calorie-neutral fat-quality substitutions high priority.
