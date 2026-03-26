# Project context for AI assistants

Read these files in order:
1. [README.md](README.md) — what this dataset is, quality case, key findings
2. [BACKGROUND.md](BACKGROUND.md) — health history, genetics, substances, locations, weighing protocol
3. [FINDINGS.md](FINDINGS.md) — testable hypotheses organized by data requirements
4. [ROADMAP.md](ROADMAP.md) — completed work, next steps, known data quality notes

## Key principles
- **The data is complete.** Every day is fully logged. A 10-calorie day is a real fasting day, not incomplete logging. Never assume missing entries.
- **Accuracy is unusually high.** Food scale always used. ~10% consistent undercount from uncounted snacking (quantified by energy balance model, uniform across phases). Incentive structure rewards precision over undercounting.
- **CICO is effect, not cause.** The working theory is a meandering set point controlled by an unknown homeostat. Simple calorie correlation has been tried and drifts.

## Documentation standard
- Every new or revised result paragraph added to `README.md`, `FINDINGS.md`, or `ROADMAP.md` must include the exact command that emits the key numbers immediately below the prose block.
- When a script writes a stable file artifact, add an `Artifact:` line below the `Command:` line naming the output CSV or plot.
- Prefer standalone analysis scripts in `analysis/` that print their headline numbers directly, so prose claims remain mechanically reproducible.

## Environment
All Python scripts run in the project virtualenv: `source .venv/bin/activate`

## Extraction pipeline
All extractors are idempotent. Re-run to regenerate CSVs from raw data.
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

# Merge pipeline (runs after extraction; each merges source-specific CSVs into canonical files)
python steps-sleep/merge_steps.py      → steps-sleep/steps.csv (Samsung + MFP/calendar hospital+walk+run backfill)
python steps-sleep/merge_exercises.py  → steps-sleep/exercises.csv (Samsung + pre-Samsung MFP runs)
python workout/merge.py                → workout/strength.csv (PDFs + Chloe xlsx + MFP circuit training)

# Analysis pipeline (order matters, no circular dependencies)
#
# Each component is estimated from independent raw data:
#   P1: weight + intake(carbs, sodium) → water corrections
#   P2: composition scans + RMR calorimetry → daily FM/FFM + fitted RMR
#   AA: composition scans + strength sessions → training delta + half-life
#   P4: combines P1 observations + P2 lean mass + AA training + intake → Kalman FM + TDEE
#
# P2 interpolates FM/FFM linearly between scans — it ignores energy balance.
# This is the weakest link: between scans, it misses weight changes that
# intake + weight observations capture. P4's Kalman FM is better between
# scans, but feeding P4 output back into P2 would create circularity.
# The 70 composition scans (median 39 days apart) limit the interpolation error.

python analysis/P0_tune_glycogen.py    — parameter tuning (run once to derive P1 constants)
python analysis/P1_glycogen_smooth.py  → analysis/P1_smoothed_weight.csv
python analysis/P2_rmr_model.py        → analysis/P2_daily_composition.csv
python analysis/P3_interpolate_weight.py → analysis/P3_daily_weight.csv
python analysis/P4_kalman_filter.py    → analysis/P4_kalman_daily.csv (v1: composition-interpolated lean)
python analysis/P4_kalman_filter_v2.py → analysis/P4_kalman_daily.csv (v2: training-adjusted lean, current default)
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
```

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
`start, end, label, detail` — manually annotated diet phases from food-level inspection. 22 epochs covering the full 2011-2026 period. Labels: `initial_restriction`, `keto_phase_1`, `weekend_fasting`, `enlightened_era`, `tirzepatide_era`, etc.

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

The Kalman filter is the basis of all TDEE and fat mass estimates used in FINDINGS.md. Every finding that references TDEE, TDEE/RMR ratio, fat mass trajectory, or metabolic adaptation depends on this model.

**Model.** Two latent states: fat mass (lbs) and TDEE (cal/day). Lean mass is a known input, not a state variable — it comes from linear interpolation between 70 InBody/BOD POD scans, adjusted for strength training effects (finding AA: 41g per session, 275-day half-life). The filter observes `smoothed_weight - known_lean = fat_mass + noise` on the 1,544 days with weigh-ins, and predicts forward through gaps using the energy balance.

**Process model (daily):**
- `fat_mass(t+1) = fat_mass(t) + (logged_intake(t) - tdee(t)) × forbes_fat_fraction / 3500`
- `tdee(t+1) = tdee(t) + 0.005 × (expected_rmr(t) × 1.15 - tdee(t)) + noise`
- Forbes fat fraction partitions surplus/deficit between fat and lean based on current fat mass. At high FM (80 lbs), ~95% of weight change is fat. At low FM (20 lbs), ~80%.
- TDEE mean-reverts toward the composition-aware expected RMR × 1.15 (sedentary activity factor) at 0.5% per day. This prevents unbounded drift during long gaps without weight observations but allows the filter to discover TDEE from the data when observations are dense.

**Constants.** These control the tradeoff between TDEE responsiveness and smoothness:
- `Q_FAT = 0.005 lbs²/day` — fat mass process noise. Small but nonzero. Absorbs the ~7% intake undercount (logged calories are lower than actual; the filter lets observations correct the resulting drift).
- `Q_TDEE = 500 cal²/day` — TDEE process noise (~22 cal/day std). Allows TDEE to respond to restriction runs and metabolic adaptation within their timescales (days to weeks). Chosen over Q=200 (smoother but too stiff to detect restriction effects) and Q=1000 (responsive but TDEE wanders too much on quiet days).
- `R = 0.97 lbs²` — observation noise. Measured from consecutive-day weight variance after glycogen+sodium correction. Represents gut contents, non-dietary hydration, and scale positioning.

**Why Q_TDEE=500.** The restriction archetype and hysteresis findings (J, K, L, M) examine TDEE changes around 3-30 day restriction runs. At Q=200, TDEE can shift at most ~17 cal over a 7-day run — too little to detect the effects being tested. At Q=500, TDEE can shift ~27 cal over 7 days, allowing the filter to track within-run metabolic changes. The cost is 0.2 lbs more FM error (3.0 vs 2.8) and slightly noisier TDEE on stable days. The innovation lag-1 autocorrelation is 0.35 at Q=500 vs 0.39 at Q=200 — closer to white noise, indicating better model fit.

**Intake is treated as logged, not corrected.** The filter does not adjust for the ~7% undercount. TDEE is therefore "logged-calorie TDEE" — the expenditure that balances logged intake against observed weight. This is systematically lower than true TDEE by the undercount. The 25 Cosmed calorimetry measurements validate the result: TDEE/RMR ratios of 1.01-1.31 (expect ~1.15-1.20 for sedentary). The gap below 1.20 at most dates is the undercount, visible but not fed back into the model. A v3 filter that estimates intake bias as a third state was tested but the bias and TDEE are nearly degenerate — they can't be separated from weight observations alone without more frequent calorimetry.

**Known limitation.** During the 2017-2018 keto era, BIA composition scans overreport fat loss due to keto-induced dehydration (reduced insulin → sodium/water excretion beyond what the glycogen/sodium model captures). The filter's FM during keto is 6-12 lbs higher than BIA readings. Outside keto, FM MAE at scans is 2.2 lbs.

**Two versions available:**
- `P4_kalman_filter.py` (v1): original, composition-interpolated lean
- `P4_kalman_filter_v2.py` (v2): training-adjusted lean, current default
- `P4_kalman_filter_v3.py` (v3): experimental, estimates intake bias as third state

`date, fat_mass_lbs_filtered, fat_mass_std_filtered, tdee_filtered, tdee_std_filtered, fat_mass_lbs, fat_mass_std, tdee, tdee_std, innovation`
- `*_filtered` columns are causal estimates available on that date (forward filter only)
- Unsuffixed `fat_mass_lbs` / `tdee` are retrospective RTS-smoothed estimates using future weigh-ins
- Innovation is the prediction error at each weight observation (should be approximately white noise)

## Known pitfalls
- **Dry lean mass vs total lean mass**: InBody XLSX column 5 ("Lean Mass") contains dry lean mass (protein + minerals, ~38 lbs), not total lean mass (~140 lbs including water). Always derive lean mass as `weight - fat_mass`. The InBody CSV exports have the correct column names.
- **Off-by-one in TDEE derivation**: Morning weight reflects food eaten *yesterday* (you weigh before eating). TDEE windows must use intake from days s..e-1, not s..e. Including day e's intake caused 200 lbs of cumulative drift before the fix.

## Raw data notes
- OXPS `--` (unknown) nutrient values render as `0`. In the merged CSV, these are empty strings from HTML/MHTML but `0` from OXPS-only years (2013-2017). Biases cholesterol, sugars, fiber low in those years.
- 4.9% of food items fail Atwater validation. Sugar alcohol products (expected) and alcohol (7 cal/g not tracked as any macro). A handful are bad MFP database entries.
- Scale body fat % is unreliable. Use composition/ for body fat data.
- Normal sleep cycle is 3am to 11am. Very little sunlight exposure.
- Temperature mean ~97.7°F, nearly 1°F below textbook. Aug 30-31 2024 readings >101°F are COVID vaccine reaction.
