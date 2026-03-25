I accept any question, no matter how personal, in service of our goal. You will not offend me.

# [Theories to test](theories.md)

# [Arguments for quality](quality.md)

# [Health background](background.md)

# [Roadmap](ROADMAP.md)

# Extracted CSVs (analysis-ready)

## intake_foods.csv
Every food item eaten from April 20, 2011 to present. 54,925 rows.
Columns: `date, meal, food, calories, carbs_g, fat_g, protein_g, cholest_mg, sodium_mg, sugars_g, fiber_g, source`

- `meal`: Breakfast, Lunch, Dinner, Snacks, Supper, or TOTAL (daily sum)
- `source`: `html`, `mhtml`, or `oxps` — provenance of extraction
- Empty nutrient values mean unknown (`--` in MFP), not zero
- Priority waterfall for overlapping dates: HTML > MHTML > OXPS
- Regenerate with `python intake/merge.py`

## intake_daily.csv
Daily nutrient totals. One row per day, 5,074 days.
Columns: `date, calories, carbs_g, fat_g, protein_g, cholest_mg, sodium_mg, sugars_g, fiber_g`

## weight/weight.csv
Daily weight. First reading of the day (fasted, post-sleep). 1,638 days from 1999 to 2025.
Columns: `date, weight_lbs, time`
Regenerate with `python weight/extract.py`

## steps-sleep/steps.csv
Daily step counts from Samsung Health. Deduplicated phone+watch (`source_type=-2`). 3,533 days from 2014 to 2024.
Columns: `date, steps, distance, speed`
Regenerate with `python steps-sleep/extract.py`

## steps-sleep/sleep.csv
Sleep periods from Samsung Health. Assigned to wake-up date. 1,279 days from 2016 to 2024.
Columns: `date, sleep_start, sleep_end, sleep_hours, time_offset`
Times are local. Regenerate with `python steps-sleep/extract.py`

## analysis/smoothed_weight.csv
Daily weight with glycogen-water correction. Removes ~4 lb swings from glycogen depletion/refill. 1,544 days.
Columns: `date, weight_lbs, glycogen_g, glycogen_correction_lbs, smoothed_weight_lbs`

- `glycogen_g`: estimated glycogen stores at time of weigh-in (lagged 1 day from intake)
- `glycogen_correction_lbs`: water weight added back (0 on normal days, up to +4 lbs after fasting)
- `smoothed_weight_lbs`: `weight_lbs + glycogen_correction_lbs` — use this for analysis
- Regenerate with `python analysis/glycogen_smooth.py`

## analysis/daily_weight.csv
Complete daily series with interpolated weight and derived TDEE. 5,429 days (every day with intake data).
Columns: `date, calories, carbs_g, observed_weight_lbs, glycogen_g, glycogen_correction_lbs, interpolated_weight_lbs, smoothed_weight_lbs, tdee, window_id`

- `observed_weight_lbs`: raw scale reading (NaN most days)
- `interpolated_weight_lbs`: simulated scale weight (matches observed where available)
- `smoothed_weight_lbs`: underlying fat mass trajectory (glycogen-corrected)
- `tdee`: derived total daily energy expenditure for the enclosing window
- TDEE derived by inverting energy balance between weigh-ins (not from a formula)
- Validated against 3 indirect calorimetry measurements: TDEE/RMR ratios of 1.14-1.27
- Regenerate with `python analysis/interpolate_weight.py`

## drugs/tirzepatide.csv
Daily tirzepatide pharmacokinetic state. 560 days from Sep 2024 to Mar 2026.
Columns: `date, dose_mg, days_since_injection, injection_date, blood_level, effective_level`

- `blood_level`: modeled serum concentration (arbitrary units) from one-compartment SC PK model. FDA parameters: t½=5.0d, Tmax=24h. Sums contributions from all prior injections (accumulation + weekly sawtooth).
- `effective_level`: blood_level adjusted for tachyphylaxis. `= blood_level × exp(-0.0217 × weeks_on_current_dose)`. Half-life of effectiveness: 32 weeks.
- Correlates with daily intake at r=-0.50 (partial). Injection day: 1652 cal. Trough day: 2220 cal.
- Regenerate with `python drugs/extract.py`

## composition/composition.csv
Body composition measurements. 49 rows across 3 eras (BOD POD, InBody partial, InBody full).
Columns: `date, weight_lbs, fat_mass_lbs, fat_pct, lean_mass_lbs, smm_lbs, bmi, era, [extended InBody columns]`
Regenerate with `python composition/extract.py`

## RMR/rmr.csv
Resting metabolic rate from indirect calorimetry. 21 measurements (3 lab Cosmed, 18 home Cosmed Fitmate).
Columns: `date, rmr_kcal, device, fasted`
Regenerate with `python RMR/extract.py`

# Raw Data Directories

## /intake
MyFitnessPal "Printable Diary" exports in OXPS, HTML, and MHTML formats spanning 2011-2026. Extractors in this directory parse all three formats.
- Column order is identical across all formats and years: Calories, Carbs, Fat, Protein, Cholest, Sodium, Sugars/Sugar, Fiber
- "Sugars" in old layout, "Sugar" in new MUI layout (2022-10+) — same column, renamed
- Limited diet variety. Same foods repeat constantly. A food name database can be built.
- MFP API can enrich food names with micronutrients (iron, potassium, vitamin D, saturated fat, added sugars). See `/myfitnesspal`.

## /composition
Body composition measurements in XLSX. Scans are backups.
- 2 BOD POD measurements (2011, 2016) — air displacement, gold standard
- 47 InBody measurements (2017-2024) — bioelectric impedance, includes segmental muscle/fat, visceral fat, water compartments

## /drugs
medicine.xlsx — injection/medication log. Extracted to `drugs/medicine.csv` + `drugs/tirzepatide.csv`.
- 80 Zepbound (tirzepatide) weekly injections, 2024-09-17 to present. Dose escalation: 2.5→5→7.5→10→12.5mg. Includes subjective strength ratings and injection sites.
- Pharmacokinetic blood level modeled from FDA parameters (t½=5d, Tmax=24h). Blood level correlates with daily intake at r=-0.50. Weekly appetite swing: 568 cal from injection day to trough.
- Tachyphylaxis modeled: effectiveness half-life 32 weeks. After 20 weeks on same dose, 65% effective.
- 3 days metformin XR 500, 2025-03-31 (prophylactic COVID exposure per COVID-OUT trial)

## /RMR
Resting metabolic rate from indirect calorimetry. XLSX + PDF scans.
- Any measure before noon is fasted
- 3 lab measurements (2011: 2415, 2012: 1956, 2016: 1700) on Cosmed. 2017 measurements marked INVALID.
- 18 home measurements (2022-2023) on personally owned Cosmed Fitmate, in clusters allowing averaging. Range 1700-2300. More data exportable from device.

## /steps-sleep
Samsung Health export. Raw CSVs with millisecond epoch timestamps.
- Steps: use `step_daily_trend` with `source_type=-2` (deduplicated phone+watch)
- Sleep: start/end times in UTC, convert using `time_offset` field. Represents phone put-down/pick-up, well-correlated with actual sleep.
- Normal sleep cycle is 3am to 11am
- Steps heuristics for outdoor detection documented in raw data notes
- Very little sunlight exposure

## /weight
bathroom-scale.xlsx — all weigh-ins from 1999 to present.
- Any measure before noon is fasted
- Weight data is reliable
- All % composition data from this scale is **UNRELIABLE** — use /composition

## /temperature
Withings thermometer. Two raw account exports (`body_temperature.csv`, `body_temperature2.csv`) merged into `temperature.csv`: 1,315 readings across 364 days, Dec 2023 – Mar 2026. Nearly daily for 15 months including the entire tirzepatide period. Multiple readings per day. Mean body temp ~97.7°F — nearly 1°F below textbook. Aug 30-31 2024 readings over 101°F are COVID vaccine reaction (Aug 28 Pfizer).

## /fatty-acids
3 blood fatty acid panels (Cleveland HeartLab 2019, OmegaQuant 2024, OmegaQuant 2025). `fatty_acids.csv` has all values. Key trend: n6:n3 ratio improving (14.8→7.5→6.1, desirable 3–5) but AA:EPA still severely elevated (36→14→20, desirable <11). Omega-3 Index rising (4.93→6.06%) but still below desirable (8–12%). Genetic FADS1 CC means efficient omega-6→AA conversion — the worst genotype for this ratio. Palmitoleic acid (carb marker) halved between 2024–2025, likely from tirzepatide-driven calorie reduction.

## /workout
Gym session PDFs (sets/reps) from 2018-2025. Chloe Workout.xlsx has 30-minute session dates.

## /myfitnesspal
Notes and script for enriching food names via MFP API to get micronutrients beyond the 8 macros in the printable diary. The food names in intake_foods.csv are the lookup keys.
