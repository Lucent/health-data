# Project context for AI assistants

Read these files in order:
1. [new_findings/README.md](new_findings/README.md) — current audit status and report index
2. [new_findings/01_claim_registry.md](new_findings/01_claim_registry.md) — current evidence status for major claims
3. [README.md](README.md) — legacy dataset narrative and historical findings
4. [BACKGROUND.md](BACKGROUND.md) — health history, genetics, substances, locations, weighing protocol
5. [FINDINGS.md](FINDINGS.md) — historical hypotheses organized by data requirements
6. [ROADMAP.md](ROADMAP.md) — completed work, next steps, known data quality notes
7. [prospective/README.md](prospective/README.md) — low-burden tracking templates and current refresh boundary

## Key principles
- **The intake date series is complete.** Every day is represented. A 10-calorie day is a real fasting day, not automatically an incomplete day.
- **Logging discipline is unusually high, but absolute error remains uncertain.** Food-scale use and full date coverage are direct facts. The claimed uniform ~10% undercount is model-dependent and must not be treated as independently validated.
- **Energy balance is an accounting identity, not a complete causal explanation.** Appetite, treatment, food environment, activity, sleep, and physiology can all affect its inputs. The repository's specific SmoothLatch/set-point mechanism is unvalidated.
- **Raw observations are not latent states.** P2/P3/P4 fat mass, lean mass, TDEE, blood level, effective level, and set point are derived/model quantities.
- **Current practical plan.** Read [`new_findings/10_current_action_plan.md`](new_findings/10_current_action_plan.md) before extending historical intervention recommendations.

## Documentation standard
- Every new or revised result paragraph added to `README.md`, `FINDINGS.md`, or `ROADMAP.md` must include the exact command that emits the key numbers immediately below the prose block.
- When a script writes a stable file artifact, add an `Artifact:` line below the `Command:` line naming the output CSV or plot.
- Prefer standalone analysis scripts in `analysis/` that print their headline numbers directly, so prose claims remain mechanically reproducible.

## Environment
Use `uv` for Python dependencies and script execution. Example:

`uv run --with pandas python analysis/EA_repository_audit.py`

## Extraction pipeline
All extractors are idempotent. Re-run to regenerate CSVs from raw data.

Raw source locations used by the extractors:
- Samsung Health export: newest `*.7z` in `/mnt/c/Users/Lucent/OneDrive/Documents/Backup/Samsung`
- `steps-sleep/extract.py` lists the archive with `7z l`, extracts only the needed Samsung CSVs to a temp dir, then writes `steps_samsung.csv`, `sleep.csv`, and `exercises_samsung.csv`

```
# Data extraction (order doesn't matter)
python intake/merge.py              → intake/intake_foods.csv + intake/intake_daily.csv + steps-sleep/mfp_exercises.csv
python intake/verify_checksums.py   — checksums, continuity, duplicates (run after every merge)
python weight/extract.py            → weight/weight.csv
python steps-sleep/extract.py       → steps-sleep/steps_samsung.csv + steps-sleep/sleep.csv + steps-sleep/exercises_samsung.csv
python composition/extract.py       → composition/composition.csv
python RMR/extract.py               → RMR/rmr.csv
python drugs/extract.py             → drugs/medicine.csv + drugs/tirzepatide.csv
python intake/atwater_check.py      — Atwater factor validation (per-item and per-day)
python intake/compare_extractors.py — cross-validates OXPS vs HTML extractors
python temperature/merge.py        → temperature/temperature.csv (body_temperature + body_temperature2, TZ-corrected)

# Merge pipeline (runs after extraction; each merges source-specific CSVs into canonical files)
python steps-sleep/merge_steps.py      → steps-sleep/steps.csv (Samsung + MFP/calendar hospital+walk+run backfill)
python steps-sleep/merge_exercises.py  → steps-sleep/exercises.csv (Samsung + pre-Samsung MFP runs)
python steps-sleep/sunlight.py         → steps-sleep/sunlight.csv (possible sunlight hours from sleep/wake + solar position)
python workout/merge.py                → workout/strength.csv (PDFs + Chloe xlsx + MFP circuit training)

# Analysis pipeline (order matters; the current version has dependency and validation leakage)
#
#   P1: weight + intake(carbs, sodium) → water corrections
#   P2: composition scans + RMR calorimetry → interpolated daily FM/FFM + fitted RMR
#   AA: composition scans + strength sessions → training delta + half-life
#   P4: 4-state Kalman (FM, LM, BUF, TDEE) anchored to weight + composition scans only.
#       Lean mass is now a state, not an interpolated input. TDEE drifts freely with rate
#       caps and no anchor. BUF is a bounded transient buffer for non-glycogen
#       water/gut variation.
#
# AUDIT WARNING: P4 initializes TDEE from P2's RMR-fitted expected_rmr, uses RMR and
# temperature in hyperparameter scoring, and writes retrospective RTS-smoothed states
# into the *_filtered columns. RMR/temperature are therefore not held out, and the
# aliases are not forward-filtered predictions. Do not use P4 for predictive or causal
# claims until the pipeline is refactored and temporally validated.

python analysis/P0_tune_glycogen.py    — parameter tuning (run once to derive P1 constants)
python analysis/P1_glycogen_smooth.py  → analysis/P1_smoothed_weight.csv
python analysis/P2_rmr_model.py        → analysis/P2_daily_composition.csv
python analysis/P3_interpolate_weight.py → analysis/P3_daily_weight.csv
python analysis/P4_kalman_filter.py    → analysis/P4_kalman_daily.csv + analysis/P4_kalman_sweep.csv
python analysis/P5_plot_models.py      → P5_plot_*.png diagnostic plots

# Standalone claim reproducers (each prints the numbers cited in README/THEORIES)
python analysis/A_energy_balance_quality.py  — ±5 lbs cumulative, undercount %
python analysis/Z_tdee_by_year.py            — TDEE by year table
python analysis/F_tirzepatide_pk.py          — r=-0.50, sawtooth, tachyphylaxis, lean mass
python analysis/B_weekend_fasting.py         — fasting microcosm (deficit recovery)
python analysis/D_food_noise_variance.py     — CV reduction, distance→intake, rebound
python analysis/E_weekly_invariance.py       — ratio 1.64, autocorrelation
python analysis/N_dietary_predictors.py      — protein leverage, meal timing, gravitostat
python analysis/Y_set_point_intake_tests.py  — 5 negative intake-side set point tests
python analysis/C_binge_analysis.py          — binge prediction AUC comparison
python analysis/AD_tdee_formula_sweep.py     — TDEE formula sweep (null: no formula beats Fitmate noise)
python analysis/AE_sleep_null.py             — sleep vs energy balance (null: r ≈ 0 for everything)
python analysis/AF_intake_variance.py        — legacy non-overlapping variance/future-FM test; protective interpretation invalidated
python analysis/AG_binge_set_point.py        — binge frequency reveals set point distance (r=-0.64)
python analysis/AH_set_point_properties.py   — inverted ratchet, dual defense, restriction prediction
python analysis/AI_expenditure_arm_timescale.py — expenditure arm adapts faster (HL≤10d) than appetite (50d)
python analysis/AJ_tirz_expenditure_defense.py  — tirzepatide suppresses the falling-phase TDEE bonus
python analysis/AK_sunlight_exposure.py          — sunlight exposure vs sleep; AD's walk effect survives sunlight control
python analysis/AL_walk_rescue_expenditure.py    — walking partially rescues tirz-suppressed expenditure defense (+95 cal/day)
python analysis/AM_lipostat_sensitivity.py       — parameter sweeps + bootstrap CIs for AG/AH/AI set point model
python analysis/AN_ratchet_profile.py            — profile likelihood for ratchet ratio, expenditure asymmetry
python analysis/AO_tirz_asymmetric_defense.py    — tirz suppresses below-SP expenditure defense selectively
python analysis/AP_overshoot_shape.py            — 90d mean surplus (r=-0.94) beats all binary thresholds
python analysis/AQ_tirz_set_point_coverage.py    — drug coverage in lbs/cal per arm, tachyphylaxis erosion
python analysis/AU_noncircular_set_point.py      — non-circular SP test from observed weigh-ins only (HL 35-50d)
python analysis/AX_drug_model.py                 — clean drug model: -74.5 cal/unit within-week, 35wk cumulative tachy
python analysis/BZ_drug_twocomponent.py          — two-component drug model: acute (-126/unit) + tonic (-63/unit) with AUC-based tachyphylaxis (K_AUC ≈ 0.0034/unit-wk)
python analysis/AV_surmount_simulation.py        — preliminary SURMOUNT simulation (15mg within 2.2% of published)
python analysis/AY_sp_from_regain.py             — derive SP adaptation from SURMOUNT-4 regain
python analysis/AZ_sp_model_search.py            — exhaustive search: SmoothLatch fits subject + SURMOUNT + STEP-1
python analysis/BD_smooth_latch.py               — SmoothLatch model: all key numbers, trial predictions, scenarios
python analysis/BG_trial_regain_search.py        — trial-only SP fit (no subject data): same latch family wins
python analysis/BR_sp_surplus_identifiability.py — 2D ridge: HL × surplus window are partially interchangeable
python analysis/BK_control_stock_model.py        — willpower-era stock model; 26% impossible days under one stock
python analysis/BA_q3d_prediction.py             — preregistered q3d intake prediction from two-component drug model
python analysis/AS_temperature_baseline.py       — empirical circadian curve (0.90°F, wake-anchored), per-reading baseline
python analysis/AT_temperature_retest.py         — re-tests all X claims with baseline-corrected temperatures

# Audit-first direct analyses (do not consume P2/P3/P4)
python analysis/EA_repository_audit.py            — coverage, continuity, staleness, enrichment
python analysis/EB_direct_future_intake.py        — locked future-intake tests for menu/context features
python analysis/EC_travel_event_study.py          — paired 14-day pre/post raw-intake travel study
python analysis/ED_tirzepatide_direct_audit.py    — raw intake/weight treatment windows + cycle confounding
python analysis/EE_tirzepatide_composition_support.py — same-device H30 + protein/training support
python analysis/EF_current_protein_menu.py         — final-90 foods, meal protein, calorie-neutral swaps
python analysis/EG_observation_provenance.py       — observed/source/device flags + nearest-measurement distance
python analysis/EI_sleep_activity_reaudit.py       — exact-date 2022+ sleep/activity tests + temporal holdout
python analysis/EJ_hyperpalatable_food_audit.py    — Fazzino FS exposure, coverage, and ordered same/next-day tests
python analysis/EK_temporal_semantics_audit.py      — static temporal inventory + corrected AF/AB/N checks
python analysis/EL_long_run_trajectory_audit.py      — direct episodes, weight-gap missingness, locked quarterly warning nulls
python analysis/EM_calendar_extreme_clustering.py    — weekday/weekend, extreme-day runs, locked calorie-momentum ranking
python analysis/EN_energy_density_nutrient_leverage.py — exact-mass coverage, meal leverage, locked future-prediction nulls
python analysis/P4b_forward_state_space.py          — opt-in raw-input forward filter + separate RTS smoother
python analysis/test_P4b_forward_state_space.py     — semantic, leakage, covariance, and determinism tests
python analysis/P4b_posthoc_baselines.py            — frozen-v1 persistence baselines; explicitly post-hoc
```

### P4b opt-in outputs

`P4b_state_daily.csv`, `P4b_covariance_daily.csv`,
`P4b_temporal_validation.csv`, and `P4b_config_selection.csv` are separate from
legacy P4 and are not consumed downstream. P4b reads only canonical intake,
observed scale weight, and composition scans. Its predicted, filtered, and
smoothed states have distinct semantics and complete uncertainty.

The fixed v1 passes all semantic tests but is not decision-grade:
later-period intervals are over-wide and the locked fat-mass forecast does not
beat a simple weight-adjusted persistence baseline. The 2025–2026 interval has
been opened; any tuning in response turns it into calibration, not a reusable
locked test. See `new_findings/11_P4b_forward_model_audit.md`.

## CSV column reference

### intake/intake_foods.csv
`date, meal, food, calories, carbs_g, fat_g, protein_g, cholest_mg, sodium_mg, sugars_g, fiber_g, source`
- `meal`: Breakfast, Lunch, Dinner, Snacks, Supper, or TOTAL
- `source`: `html`, `mhtml`, or `oxps`
- Empty nutrient = unknown (`--` in MFP), not zero
- Priority waterfall: HTML > MHTML > OXPS

### intake/intake_daily.csv
`date, calories, carbs_g, fat_g, protein_g, cholest_mg, sodium_mg, sugars_g, fiber_g`

### intake/diet_epochs.csv
`start, end, label, detail` — 27 manually annotated, overlapping diet phases. They do not fully partition 2011-2026. Labels include `initial_restriction`, `keto_phase_1`, `weekend_fasting`, `enlightened_era`, and `tirzepatide_era`.

### weight/weight.csv
`date, weight_lbs, time` — first reading of day (fasted, post-sleep)

### steps-sleep/mfp_exercises.csv
`date, name, calories, minutes, source` — all MFP exercise entries (extracted alongside food by intake/merge.py). Consumed by merge_steps.py, merge_exercises.py, and workout/merge.py.

### steps-sleep/steps_samsung.csv
`date, steps, distance, speed` — Samsung Health deduplicated phone+watch (starts 2014-04-24)

### steps-sleep/steps_calendar.csv
`date, steps` — hospital shifts from calendar not in MFP (16 dates, 5500 steps each)

### steps-sleep/steps.csv
`date, steps, distance, speed` — canonical merged steps (Samsung + MFP/calendar backfill). Use this for analysis.

### steps-sleep/exercises_samsung.csv
`date, start_time, end_time, duration_min, exercise_type, exercise_label, title, source_type, count, distance, calorie, time_offset, pkg_name, datauuid`

### steps-sleep/exercises.csv
`date, type, duration_min, distance, calorie, source` — canonical merged cardio exercises (Samsung + pre-Samsung MFP runs). Types: walking, running, bike, hiking, indoor_bike, pilates, yoga, other.

### workout/strength.csv
`date, duration_min, source` — canonical strength training dates. All 30 min. Sources: pdf (ActivTrax), chloe (personal trainer), mfp (circuit training).

### steps-sleep/sleep.csv
`date, sleep_start, sleep_end, sleep_hours, time_offset` — assigned to wake-up date

### composition/composition.csv
`date, measured_at, source, device, era, weight_lbs, fat_mass_lbs, fat_pct, lean_mass_lbs, smm_lbs, bmi, [extended columns]`
- `era`: `bodpod`, `inbody_partial`, `inbody_full`, `inbody_summary`
- `lean_mass_lbs` = `weight_lbs - fat_mass_lbs` (always derived, not raw InBody "Lean Mass" column which is dry lean mass)
- Extended columns include water compartments, visceral fat, segmental muscle/fat, ECW ratios, waist circumference, phase angle
- Extracted from InBody CSV exports + BOD POD rows in XLSX. InBody preferred when dates overlap.

### RMR/rmr.csv
`date, rmr_kcal, device, fasted`

### drugs/tirzepatide.csv
`date, dose_mg, days_since_injection, injection_date, blood_level, effective_level`
- `blood_level`: PK model (FDA: t½=5d, Tmax=24h, ka=3.31/day), sums all prior injections
- `effective_level`: `blood_level × exp(-0.0217 × weeks_on_current_dose)` — tachyphylaxis-adjusted

### analysis/P1_smoothed_weight.csv
`date, weight_lbs, glycogen_g, glycogen_correction_lbs, sodium_correction_lbs, smoothed_weight_lbs`
- Corrections centered on median (typical day ≈ 0 correction)
- 1-day lag: morning weight reflects previous day's intake

### analysis/P2_daily_composition.csv
`date, weight_lbs, fm_lbs, ffm_lbs, fat_pct, expected_rmr`
- FM/FFM linearly interpolated between composition anchors
- `expected_rmr`: individually fitted model (RMR = a×FFM_kg + b×FM_kg + c)

### analysis/P3_daily_weight.csv
`date, calories, carbs_g, sodium_mg, observed_weight_lbs, glycogen_g, glycogen_correction_lbs, interpolated_weight_lbs, smoothed_weight_lbs, tdee, window_id`
- Window-based TDEE derivation (constant per 7+ day window)

### analysis/P4_kalman_daily.csv

The Kalman filter underlies many TDEE, FM, and LM estimates cited in FINDINGS. It is currently an exploratory retrospective model, not an independently validated predictive series. See `new_findings/00_initial_integrity_and_safety_audit.md`.

**4 latent states.** FM (lbs), LM (lbs), BUF (lbs — bounded transient buffer for non-glycogen, non-sodium water/gut variation, ±3 lb cap, HL=5d, decays toward 0), TDEE (cal/day, free random walk with rate cap, no anchor).

**Forbes-coupled energy partition.** `deficit = dFM × 3500 + dLM × 1100`, `dFM/dW = ff(FM, workouts/wk)` where `ff_base = 1 − C²/(C+FM_kg)²`, C = 10.4 kg, plus a workout modulation per CB (α = 0.084 per workout/wk). Hard daily-rate caps from physiology: FM gain ≤ +0.5 lb/day; FM loss ≤ FM × 31/3500 lb/day (Alpert ceiling); |dTDEE/dt| ≤ 100 cal/day.

**Calorimetry is not held out in the current search.** P2 fits `expected_rmr` from RMR, P4 uses it for initialization, and the P4 composite score includes RMR and temperature statistics.

`date, fat_mass_lbs_filtered, fat_mass_std_filtered, tdee_filtered, tdee_std_filtered, fat_mass_lbs, fat_mass_std, tdee, tdee_std, innovation, lm_lbs, buf_lbs`
- `fat_mass_lbs`, `tdee`, and their current `_filtered` aliases all contain retrospective RTS-smoothed states.
- The four exported FM/TDEE standard-deviation columns are currently blank placeholders.
- `lm_lbs` and `buf_lbs` expose the lean-mass and transient-buffer states for downstream analyses (CB, CC, etc.).
- Innovation is the weight prediction residual (approximately white noise; current lag-1 ACF ≈ 0.08).

## Known pitfalls
- **Dry lean mass vs total lean mass**: InBody XLSX column 5 ("Lean Mass") contains dry lean mass (protein + minerals, ~38 lbs), not total lean mass (~140 lbs including water). Always derive lean mass as `weight - fat_mass`. The InBody CSV exports have the correct column names.
- **Off-by-one in TDEE derivation**: Morning weight reflects food eaten *yesterday* (you weigh before eating). TDEE windows must use intake from days s..e-1, not s..e. Including day e's intake caused 200 lbs of cumulative drift before the fix.

## Raw data notes
- OXPS `--` (unknown) nutrient values render as `0`. In the merged CSV, these are empty strings from HTML/MHTML but `0` from OXPS-only years (2013-2017). Biases cholesterol, sugars, fiber low in those years.
- 4.9% of food items fail Atwater validation. Sugar alcohol products (expected) and alcohol (7 cal/g not tracked as any macro). A handful are bad MFP database entries.
- Scale body fat % is unreliable. Use composition/ for body fat data.
- The reliable sleep regime is 2022 onward (1,504/1,520 dates through 2026-02-28). Thirty-two pre-2022 rows under 2 hours are nap-like; legacy AE/BS include them and use next-record rather than exact next-date joins.
- `steps-sleep/sunlight.csv` is possible wake-to-sunset daylight, not measured outdoor exposure or illuminance.
- On a forward-sorted series, `shift(-1).rolling(h)` is mostly trailing, not a future window. Use explicit shifts for `t+1..t+h`, and do not call a row shift “next day” after filtering or anchoring on a sparse table.
- Normal sleep cycle is approximately 2am to 10am in the reliable 2022+ regime.
- Temperature mean ~97.7°F, nearly 1°F below textbook. Aug 30-31 2024 readings >101°F are COVID vaccine reaction.
