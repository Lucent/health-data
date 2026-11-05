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

### Audit 2026-11 (DD, DE, DF, BD, P4 fixes)
- **DD** — AG's binge-rate set point signal is trailing-window reverse causation (r −0.64 → −0.12 forward, anti-lipostat residual); forward CCF inside an intake-only AR null; +12 cal/day/lb residual at 90 d only.
- **DE** — `drugs/extract.py` now uses injection timestamps: daily-mean exposure over the eating window, calendar-day cycle labels (day 1 = peak). AX coefficient −74.5 → −89.9 cal/unit; BZ acute:tonic 2:1 → 1:1.
- **DF** — drug-era intake rise is frozen-set-point pressure (30-40 cal/lb, ΔAIC −38 vs AUC-tachyphylaxis; only form that transfers to the 12.5 mg era). Subject-fit latch rejected by the drug era; trial-fit latch supported.
- **BD** — regain denominator fixed (was `FM_end/0.40` ≈ 120 lb). "+14.1% vs +14.0%" withdrawn; subject-fit latch gives +5.5%.
- **P4** — energy-balance closure now includes lean-mass energy (+4.5 lb, not ~0) and is documented as non-identifying for the intake multiplier.
- **DG** — canonical gap-form model (supersedes BZ/CE/AV): fit once on the drug era, subject FM path to 1.2 lb RMSE, SURMOUNT-4 regain shape at m = 1, withdrawal regain ~2/3 in all arms with a frozen SP; cohort dose-response concavity unresolved; q3d forecast written to `analysis/DG_q3d_forecast.csv`.
- Known stale: `tests/test_steps_sleep_extract.py` (2 tests) references functions that no longer exist in `steps-sleep/extract.py`.

### Narrative-driven hypotheses (tested)
- **Expenditure arm timescale** → finding AI. HL ≤10 days, faster than appetite (50d).
- **Tirzepatide and expenditure defense** → finding AJ. Drug suppresses 63% of falling-phase TDEE bonus.
- **Session-level walk regimes** → finding AD (V folded in). Walk sessions predict RMR; paired daylight walks characterized.
- **Temperature by phase** → finding X. Confound corrected; drug does not override thermostat.
- **Sunlight-intake correlation** → finding AK. Sunlight weak for intake; AD walk effect independent of sunlight.

## Data gaps
- **2018-02-20 to 2018-02-28**: Still OXPS-only. Re-exportable from MFP.

## Untested hypotheses

### Frozen-set-point frame (highest leverage after DF)

DF makes the working model: drug suppression linear in exposure (~−85 cal/unit, no erosion) against eating pressure P·(SP₀ − FM) with P ≈ 30-40 cal/day/lb and SP₀ frozen at the pre-drug fat mass (equivalently the trial-fit SmoothLatch). These tests would confirm, refute, or sharpen it.

**q3d era (from 2026-04-01) as the decisive test.** Switching to 5.2 mg every 3 days roughly halves the sawtooth while keeping mean exposure near 12 units. Tachyphylaxis (AUC-driven) and the frozen-SP gap make different predictions for the *mean*: the gap model says mean intake tracks FM only (flat if FM stays at 61), tachyphylaxis says it keeps drifting up with AUC. BA's preregistered q3d prediction uses the tachyphylaxis form; add a gap-form prediction before the q3d intake data are merged, then score both.

**Pressure constant across fat mass.** P is estimated at FM 61-85 lb. The leptin-floor reading (AC+BK) says it should fall toward zero below ~18 lb FM. Not testable on-drug; the 2012-14 willpower era is the only low-FM window and has no identifiable gap (DD). Leave as a prediction.

**Discontinuation forecast, gap form.** At FM 61 and SP₀ 85, pressure ≈ 35 × 24 ≈ 850 cal/day at drug stop, decaying as FM returns. Forward-integrate with Forbes partition and this subject's TDEE: predicts the time to recover half the loss. Becomes testable if the drug stops; write the forecast down now.

**Rebuild BZ/CE/AV on the gap form — done (DG).** Subject FM path reproduced to 1.2 lb RMSE; SURMOUNT-4 regain shape reproduced at m = 1; cohort m = 1.56 (SURMOUNT-1) / 0.81 (SURMOUNT-2). Two things DG leaves open:

**Cohort dose-response concavity — partly resolved.** Heterogeneity is out (DI). The subject's drug era mildly prefers EC50 ≈ 12-20 units to linear and excludes EC50 ≈ 5; at EC50 = 12 the SURMOUNT-1 arms are within ~2 points. Remaining: no single (m, P_cohort) fits SURMOUNT-1, the SURMOUNT-4 lead-in and the SURMOUNT-4 regain together (DG §5b). Pin EC50 from Lilly's published exposure-response (ng/mL → dose units), and treat SURMOUNT-4's faster-than-SURMOUNT-1 lead-in as a between-trial fact before asking the model to fit it. The q3d era (6.25 mg every 3 days, mean exposure ≈ 12 with a small swing) is the within-subject linearity test: at equal mean exposure to the 12.5 mg weekly era, Emax at EC50 = 12 predicts a slightly *stronger* mean suppression than linear would, because the smaller swing wastes less of the curve's top.

**Lean mass as the regulated variable (Dulloo).** DH: lean-mass deficit at given FM is the only forward appetite signal pre-drug (partial +0.4-0.6). Discriminate protein-stat from cumulative-restriction history with training-induced LM (CB/CC: +30-40 g per session): if LM gained by training lowers subsequent intake at matched FM and past intake, the LM level is the signal; if not, it is the history. Also: does the LM deficit predict RMR residuals (Dulloo's thermogenesis side) and protein-seeking? DJ runs these.

**The 2014 switch.** No set point gap peaks at the 2014Q3-Q4 eruption (DH). Date the restriction → free-eating switch objectively (HMM/changepoint on daily intake) and fit a termination hazard on FM gap, LM deficit, time-in-regime, variance and season. DK runs this.

**Seasonality.** Pre-drug 2014+ binge rate is 1% in February and 8% in October (DH appendix). Every pooled analysis in FINDINGS carries this unmodelled; add month-of-year to the forward tests and check whether the autumn peak is light (AK), walks (AD) or something else.

**Cohort P from regain timing.** The SURMOUNT-4 trajectory wants P ≈ 15 with a frozen SP; the subject's drug era says 43. Re-fit P per withdrawal trial from the published regain curves (AY's digitised SURMOUNT-4 curve; STEP-1 extension and STEP-4 curves need digitising) and test whether P correlates with baseline BMI or trial drug. The set point does not need to move on-drug to fit any of them; drop the latch from the trial model unless something requires it.

**Latch firing during the drug era.** Partly answered by DF: the subject-fit latch (tol 3, hold 14) adds nothing, the trial-fit latch (tol 1, hold 30) matches the frozen SP. What remains is whether tol/hold can be pinned more tightly: FM held within ±1.5 lb for 106 days at the plateau, which under the trial-fit latch moves SP ~5 lb toward FM; whether the plateau's intake drift (−185 cal/day across the window, the opposite of tachyphylaxis) is consistent with that much SP descent is a direct test.

**Conservation-mode signature in 2012-13 composition.** The leptin-floor reading predicts BOD POD scans during deep restriction (FM ≲ 18 lb) should show water-compartment / FFM-preservation patterns characteristic of starvation-mode regulation, not uniform fat partitioning. Extract ECW/ICW and segmental data from `composition.csv`'s extended columns for 2012-14. If the floor is real, the composition signatures should switch on near FM ≈ 18 lb.

**Receptor-occupancy curvature in pressure.** AM finds 55 cal/lb is constant from -10 to +20 lb gap. The leptin-floor reading predicts the curve should *saturate at very low FM* because absolute leptin levels approach the receptor's lower detection limit. Re-bin the pressure curve including FM as a covariate; look for a flattening at FM < 25 lb that the linear model misses.

**Asymptote heterogeneity in SURMOUNT-1.** Receptor-dynamics says asymptote depth scales with `D₀ / 55 cal/lb`, so between-subject variation should sort by per-lb pressure (which depends on baseline SP stability). SURMOUNT-1 published responder distribution. If asymptote depth correlates with pre-trial FM stability (proxy for how recently the SP recalibrated), the latch interpretation is supported.

**Trial transfer to STEP-4 / SURMOUNT-3 maintenance.** BD/AZ used SURMOUNT-4 and STEP-1 only. Add STEP-4 and SURMOUNT-3 maintenance arms. The SmoothLatch predicts maintenance arms regain less than withdrawal arms because the drug keeps FM stable long enough for the SP to latch. Tightens parameters without using subject data.

### Cross-finding seams

**Archetype × set point interaction.** J shows high-step and low-protein cuts fail. AG shows SP distance drives binges. Regress post-run outcome on both archetype and SP distance at run end. If archetype effects vanish after controlling for SP distance, they're proxies. If both survive, they're independent failure modes.

**Protein leverage × SP distance.** N shows within-day protein leverage (r=-0.34). AG shows SP distance drives appetite. Bin days by SP distance and compute protein leverage slope within each bin. If the set point overrides protein satiety at large deficits, the slope flattens.

**Walk sessions around the 2013 inflection.** AC says the trigger was behavioral. AD says walks raise RMR. Did walk frequency change around the inflection? MFP exercise entries predate Samsung (Apr 2014). If walks dropped, the RMR collapse (TDEE/RMR 0.977) could partly be walk-frequency, not just set point defense.

**Set point movement during restriction runs.** AH shows runs ending above SP stick, below SP rebound. At 50-day HL, a 7-day run moves SP ~0.45 lbs. Do some archetypes move the SP faster? Under the SmoothLatch the answer is bounded: only stable-tail runs should move the SP at all.

**Discontinuation forecast.** Superseded by the gap-form forecast above: DF shows the on-drug SP did *not* latch onto the 61-lb plateau (the subject-fit latch that would have is rejected), so the gap at discontinuation is ~20-25 lb and predicted regain is large, not single-digit.

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
