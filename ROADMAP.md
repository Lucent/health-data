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
[FINDINGS.md](FINDINGS.md) — reorganized 2026-04-25. The Set Point section is now identified trial-first (SmoothLatch on SURMOUNT-4 / STEP-1, BG trial-only fit) with three subject-side independent confirmations (AG binge rate, P3 intake-free FM, AU non-circular); the autocorrelated mean-surplus regression (AM/AP) is demoted to consistent supporting evidence with the BR identifiability ridge made explicit. 2013 inflection and BK control-stock attempt reframed under a leptin-floor / receptor-recalibration reading. New asymptote section under Tirzepatide explains why GLP-1 weight loss flattens. 21 dead-end search scripts and ~18 MB of intermediate artifacts removed; nothing referenced in docs was deleted.

### Narrative-driven hypotheses (tested)
- **Expenditure arm timescale** → finding AI. HL ≤10 days, faster than appetite (50d).
- **Tirzepatide and expenditure defense** → finding AJ. Drug suppresses 63% of falling-phase TDEE bonus.
- **Session-level walk regimes** → finding AD (V folded in). Walk sessions predict RMR; paired daylight walks characterized.
- **Temperature by phase** → finding X. Confound corrected; drug does not override thermostat.
- **Sunlight-intake correlation** → finding AK. Sunlight weak for intake; AD walk effect independent of sunlight.

## Data gaps
- **2018-02-20 to 2018-02-28**: Still OXPS-only. Re-exportable from MFP.

## Untested hypotheses

### Receptor-dynamics frame (highest leverage)

The trial-first reorganization makes the SmoothLatch + leptin-floor reading the working model. These tests would confirm, refute, or sharpen it.

**Latch firing during the drug era.** The on-drug set-point HL is currently bracketed at 40-120 d (AQ surface is flat). The SmoothLatch says it should be much faster *during stable holds* and zero *during active loss*. Mark every drug-era day by whether FM had been within ±3 lbs of any reference for 14+ consecutive days. Estimate SP movement separately on latched vs unlatched days. The model predicts an order-of-magnitude difference. Distinguishes stability-gated latch from a generic slow EMA on-drug.

**Conservation-mode signature in 2012-13 composition.** The leptin-floor reading predicts BOD POD scans during deep restriction (FM ≲ 18 lb) should show water-compartment / FFM-preservation patterns characteristic of starvation-mode regulation, not uniform fat partitioning. Extract ECW/ICW and segmental data from `composition.csv`'s extended columns for 2012-14. If the floor is real, the composition signatures should switch on near FM ≈ 18 lb.

**Receptor-occupancy curvature in pressure.** AM finds 55 cal/lb is constant from -10 to +20 lb gap. The leptin-floor reading predicts the curve should *saturate at very low FM* because absolute leptin levels approach the receptor's lower detection limit. Re-bin the pressure curve including FM as a covariate; look for a flattening at FM < 25 lb that the linear model misses.

**Asymptote heterogeneity in SURMOUNT-1.** Receptor-dynamics says asymptote depth scales with `D₀ / 55 cal/lb`, so between-subject variation should sort by per-lb pressure (which depends on baseline SP stability). SURMOUNT-1 published responder distribution. If asymptote depth correlates with pre-trial FM stability (proxy for how recently the SP recalibrated), the latch interpretation is supported.

**Trial transfer to STEP-4 / SURMOUNT-3 maintenance.** BD/AZ used SURMOUNT-4 and STEP-1 only. Add STEP-4 and SURMOUNT-3 maintenance arms. The SmoothLatch predicts maintenance arms regain less than withdrawal arms because the drug keeps FM stable long enough for the SP to latch. Tightens parameters without using subject data.

### Cross-finding seams

**Archetype × set point interaction.** J shows high-step and low-protein cuts fail. AG shows SP distance drives binges. Regress post-run outcome on both archetype and SP distance at run end. If archetype effects vanish after controlling for SP distance, they're proxies. If both survive, they're independent failure modes.

**Protein leverage × SP distance.** N shows within-day protein leverage (r=-0.34). AG shows SP distance drives appetite. Bin days by SP distance and compute protein leverage slope within each bin. If the set point overrides protein satiety at large deficits, the slope flattens.

**Walk sessions around the 2013 inflection.** AC says the trigger was behavioral. AD says walks raise RMR. Did walk frequency change around the inflection? MFP exercise entries predate Samsung (Apr 2014). If walks dropped, the RMR collapse (TDEE/RMR 0.977) could partly be walk-frequency, not just set point defense.

**Set point movement during restriction runs.** AH shows runs ending above SP stick, below SP rebound. At 50-day HL, a 7-day run moves SP ~0.45 lbs. Do some archetypes move the SP faster? Under the SmoothLatch the answer is bounded: only stable-tail runs should move the SP at all.

**Discontinuation forecast.** Apply the SmoothLatch with current SP/FM. The on-drug SP has had ~6 mo of FM-stable holds at 60 lb to latch onto, so the gap at hypothetical discontinuation is small and predicted regain is single-digit %. Becomes testable if/when drug stops.

**Weekend fasting × set point.** B showed 7 weekend fasts (Oct-Nov 2019) where deficit vanished through expenditure. SP distance during those weeks would say whether expenditure defense (B, 9d HL) and the slow eating-pressure SP (50d HL) are decoupled at that timescale.

### Standalone (ready now)

**Conditional intake-discard via within-week instrument.** Earlier attempts to estimate an unconditional intake-bias state were degenerate with TDEE. The clean instrument is *within-week injection-cycle variation*: does a discard fraction differ between drug-peak days and trough days? If yes, the "calories above set point are dissipated" hypothesis becomes identifiable using the same within-week structure AX uses for the appetite coefficient. Adds a conditional bias term to the Kalman: `effective_intake = intake × (1 − φ × indicator)` where the indicator is exogenous (drug-peak day, walk day, surplus magnitude). Test on the q3d data once enough has accumulated, where the indicator structure changes meaningfully.

**Sustained-overeating asymmetry.** Does sustained overeating (>2800 cal/day, >5 days) move the SP faster than symmetric latching predicts? AH retracted the 2011-13 ratchet but that test was confounded by the willpower era; this is the cleaner restatement on 2014+ data only.

**Plateau dynamics on GLP-1s.** Weight drops, plateaus, then either (1) resumes from the plateau or (2) jumps down as if the plateau never happened. Glycogen smoothing should distinguish — if plateaus are water masking fat loss, model 2 is correct.

**Calorie misestimation by food.** Scale per-food calorie entries up or down to find which foods produce better fit with weight. Frequently underestimated restaurant meals would show up.

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
3. **Branch memory with stronger controls** — K's core branch effect is well-established (S, T robustness folded in). Remaining: does it survive controls for travel, prior cut sequencing, monotony, and meal-pattern structure?
4. **Food-template effects at matched macros** — do recurring meal templates produce different fat/TDEE outcomes than calorie-equivalent alternatives?
5. **MFP API enrichment** — micronutrients for top 200 foods.
6. **Join table** — single daily table with all derived variables.
