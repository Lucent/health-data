# Roadmap

## Completed

### Data extraction
- `intake/merge.py` → `intake_foods.csv` (54,925 items), `intake_daily.csv` (5,429 days)
- `weight/extract.py` → `weight.csv` (1,693 days)
- `steps-sleep/extract.py` → `steps_samsung.csv`, `sleep.csv`, `exercises_samsung.csv`
- `steps-sleep/merge_steps.py` → `steps.csv` (canonical merged, 4,540 days including MFP/calendar backfill)
- `composition/extract.py` → `composition.csv` (70 measurements, 4 eras)
- `RMR/extract.py` → `rmr.csv` (25 Cosmed Fitmate measurements: 3 at sports medicine clinic, 22 at home)
- `drugs/extract.py` → `medicine.csv`, `tirzepatide.csv`
- `workout/merge.py` → `workout/strength.csv` (session dates)

### Analysis pipeline
- P0: parameter tuning → P1: glycogen + sodium smoothing → P2: RMR model → P3: weight interpolation → P4: Kalman filter → P5: diagnostic plots

### Findings
[FINDINGS.md](FINDINGS.md) — 26 lettered analyses (A-Z), each with a standalone script.

## Data gaps
- **2018-02-20 to 2018-02-28**: Still OXPS-only. Re-exportable from MFP.

## Untested hypotheses

### Ready now (narrative-driven)

These follow directly from the restructured FINDINGS.md narrative — each sits at a seam between two established findings.

**Expenditure arm timescale.** AG gives the appetite arm a 50-day half-life. K shows the expenditure arm exists. Sweep EMA half-lives for TDEE residual (Kalman TDEE minus composition-predicted RMR) against fat mass distance, as AG did for binge rate. If the expenditure arm has a longer half-life, that explains why K's hysteresis appears persistent. Connects AG↔K.

**Archetype × set point interaction.** J shows high-step and low-protein cuts fail. AG shows SP distance drives binges. Regress post-run outcome on both archetype and SP distance at run end. If archetype effects vanish after controlling for SP distance, they're proxies. If both survive, they're independent failure modes. Connects J↔AG.

**Protein leverage × SP distance.** N shows within-day protein leverage (r=-0.34). AG shows SP distance drives appetite. Bin days by SP distance and compute protein leverage slope within each bin. If the set point overrides protein satiety at large deficits, the slope flattens. Connects N↔AG.

**Walk sessions around the 2013 inflection.** AC says the trigger was behavioral. AD says walks raise RMR. Did walk frequency change around the inflection? MFP exercise entries predate Samsung (Apr 2014). If walks dropped, the RMR collapse (TDEE/RMR 0.977) could partly be walk-frequency, not just set point defense. Connects AC↔AD.

**Tirzepatide and the expenditure defense direction.** F shows 206 cal metabolic clawback. K says falling phases have elevated TDEE. But tirz is falling — shouldn't TDEE be elevated? Compare Kalman TDEE/RMR during tirz-era falling vs pre-tirz falling at matched FM. If tirz falling shows lower TDEE/RMR, the drug suppresses the expenditure arm. Connects F↔K.

**Set point movement during restriction runs.** AH#4 shows runs ending above SP stick, below SP rebound. At 50-day HL, a 7-day run moves SP ~0.45 lbs. Do some archetypes move the SP faster? Connects AH↔J.

**Discontinuation forecast.** AG gives current SP (~63 lbs FM), sigmoid curve, two candidate rates (50d off-drug, 165d on-drug). Simulate discontinuation binge trajectories under both. Becomes testable if/when drug stops. Connects AG→future.

**Weekend fasting × set point.** B showed 7 weekend fasts where deficit vanished through expenditure. What was SP distance during those weeks (Oct-Nov 2019)? If FM was above SP, expenditure defense should have been minimal — implying acute expenditure defense (B) operates on a faster timescale than the SP. Connects B↔AG.

### Ready now (standalone)

**Set point shift conditions.** When does the set point move? Find periods of sustained overeating (>2800 cal/day for >5 consecutive days) and test whether the subsequent weight floor permanently increases. Is the ratchet one-way?

**Plateau dynamics on GLP-1s.** Weight drops, plateaus, then either (1) resumes from the plateau or (2) jumps down as if the plateau never happened. Glycogen smoothing should distinguish — if plateaus are water masking fat loss, model 2 is correct.

**Calorie misestimation by food.** Scale per-food calorie entries up or down to find which foods produce better fit with weight. Frequently underestimated restaurant meals would show up.

**Sunlight-intake correlation.** Sleep/wake times × local sunrise/sunset → hours of daylight overlap. Very little sunlight exposure makes any effect easier to detect.

**Circadian misalignment.** 3am-11am sleep, ~4 hours shifted from solar cortisol rhythm, near-zero zeitgeber correction. This is constant in the data — invisible as a variable but extremely unusual vs study populations.

**Cold intolerance and protein.** Temperature data + daily protein → does protein intake predict next-day body temperature? Reported: "couple days of low protein and already less cold intolerance."

**Low-protein / BCAA restriction.** The 47-day low-protein period (Jul-Aug 2023, 45g/day) is visible in the data. Correlate with temperature and weight trajectory.

### Needs MFP API enrichment

**Omega-6:3 ratio and RMR.** Does trailing 30-day omega-6:3 > 10:1 correspond to lower derived TDEE? Requires fat subtype data from MFP API for the top ~200 foods.

**NOVA classification and binges.** Ultra-processed food fraction → binge prediction. Feasible from food names.

**Flavorless oil undercounting.** Are oil-heavy cooking days more fattening per calorie than flavorful foods? Classify food items as oil-heavy.

**AGEs as food badness metric.** A reductionist candidate beyond NOVA/UPF. Classify foods by AGE content.

### Needs external data

**Gut microbiome.** Heavy childhood erythromycin + 35 years no animal protein = textbook persistent dysbiosis setup. Not testable from current data. A 16S or shotgun assay would reveal whether composition is atypical.

**Antihistamine contribution.** [H1 antihistamines associated with obesity: 10 kg (NHANES)](https://doi.org/10.1038/oby.2010.176). Cetirizine 10mg→5mg→levocetirizine 2.5mg transition occurred Dec 2024 – Apr 2025, confounded by tirzepatide dose escalation. See BACKGROUND.md substance timeline for exact dates.

## Next analysis steps

1. **Food monotony and repetition structure** — distinct-food count, recurrence interval, menu compression, transition entropy from `intake_foods.csv`. Test whether these predict binge persistence, rebound, and metabolic failure.
2. **Regime-switch transition models** — treat potato, keto, travel, baseline, and tirzepatide as distinct transition systems. Compare escalation rates across named regimes.
3. **Branch memory with stronger controls** — does the rising/falling branch effect survive controls for travel, prior cut sequencing, monotony, and meal-pattern structure?
4. **Food-template effects at matched macros** — do recurring meal templates produce different fat/TDEE outcomes than calorie-equivalent alternatives?
5. **Restriction run sequencing** — prior 30-day deficit, days since last cut, refeed intensity, recent high-intake streaks.
6. **Session-level walk regimes** — deliberate daylight walks vs incidental steps, validate Samsung labels.
7. **Temperature by phase** — needs more pre-drug falling-phase data to be identified.
8. **MFP API enrichment** — micronutrients for top 200 foods.
9. **Join table** — single daily table with all derived variables.
