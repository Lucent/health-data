# Findings and testable theories

Tested hypotheses include results with `Command:` lines pointing to the standalone script that reproduces the numbers. Untested hypotheses are organized by data requirements.

## Ready now (intake_daily.csv + weight.csv)

**Glycogen-water smoothing.** Derive a formula to smooth daily weight perturbations from glycogen-water binding. Each gram of carbs retains ~3g of water, depleted over ~2 days of restriction, recaptured within hours of refeeding. Validate against the Oct-Nov 2019 weekend fasts (weight drop should far exceed caloric deficit, with immediate bounce-back on refeeding). This is prerequisite infrastructure for everything below — without it, daily weight noise drowns real signal.

**Hidden set point estimation.** Find long runs (months) of stable weight (±3 lbs) where intake does not support no change. This is direct evidence for a set point that can only drift under certain circumstances. Extract the set point as a slow-drifting latent variable — Kalman filter or similar state-space model rather than the breakpoint approach in analysis-bmr-fit (which didn't work well).

**Variable RMR trendline.** Fit intake vs. smoothed weight with a slow-changing RMR variable. Find what high/low RMR correlates with — a macronutrient, a ratio, a pattern. Graph derived RMR over time. The 3 indirect calorimetry measurements (2011, 2012, 2016) are anchor points.

**Set point shift conditions.** Look for the set point abruptly increasing or decreasing and examine the days/weeks when it appeared to happen. Can only multi-day binges shift set point up (e.g. 3 consecutive days of 3000 calories)? Do the same calories spread out cause no shift? What conditions push it down — are there magic months where 3 lbs are lost and not regained for a long time?

**Weekend fasting as set point microcosm — TESTED.** Seven consecutive weekend 36-hour fasts (Oct-Nov 2019) with daily weigh-ins provide a controlled experiment in miniature.

**Results:** Each weekend fast creates a ~3,300 cal deficit (TDEE ~1,750 × 2 minus ~200 consumed). The Kalman filter shows 0.76 lbs of real fat loss by Monday — 80% of the expected 0.94 lbs. But by the following Friday, fat mass is back to starting level (+0.15 lbs). Seven fasting weekends × 3,300 cal deficit = 23,100 cal total = 6.6 lbs expected. Actual net fat change: approximately zero.

There is no compensatory overeating. Post-fast weekday intake (2,369 cal) is actually 78 cal/day BELOW the pre-fast baseline (2,447). The deficit disappears through **reduced expenditure**, not increased intake. The body recaptures the weekend deficit by running more efficiently during the week. Over 7 cycles, ~23,000 cal of deficit produced zero net fat loss — a metabolic efficiency adjustment of ~470 cal/day.
Command: `python analysis/B_weekend_fasting.py`

This is the set point defense in microcosm. Acute restriction produces real short-term fat loss. The body recovers it not by making you eat more, but by burning less. The TDEE/RMR ratio dropping to 1.02 during sustained restriction is the same mechanism operating at a longer timescale.

**Set point and protein.** Does set point only move up in the presence of protein, to capture high-effort hunted meals? Studies show 20% protein diets result in same intake as 10% protein + 500 calories. Test whether high-protein days are more likely to precede set point increases.

**Binge prediction and food noise — TESTED, CENTRAL PREDICTION FALSIFIED.** What best predicts whether a day will become a binge (>2800 cal)? The [food noise hypothesis](https://lucent.substack.com/p/craving-food-noise) predicts that binge probability correlates with cumulative distance below set point — not yesterday's intake, not a single macronutrient, but the gap between current weight and the weight the body is defending.

**Result (617 binges across 5,428 days):** Distance from set point does NOT predict binges. AUC=0.49 — literally a coin flip, across all set point window sizes tested (90d, 180d, 365d). Binge rate shows no monotonic pattern: 14.9% at 5 lbs below set point, 10.1% at 2 lbs above, 15.7% at 5 lbs above.

What actually predicts binges: yesterday's calories (AUC=0.74), trailing 3-day calories (0.74), protein % (0.71), cumulative 30-day deficit (0.67). All dietary variables combined: AUC=0.81. Adding set point distance improves this to 0.82 — marginal. **Binges are predicted by recent eating patterns and low protein, not by any set point distance.**

**The tirzepatide experiment is definitive.** At the same distance from set point (-2 to +2 lbs), binge rate drops from 11.7% to 3.4%. Logistic regression: distance β=+0.009 (zero), tirzepatide β=-0.497 (strong). The drug doesn't move the set point. It silences food noise at the same metabolic position.
Command: `python analysis/C_binge_analysis.py`

**Refined tests (inspired by the [food noise essay](https://lucent.substack.com/p/craving-food-noise)):**

The essay distinguishes hunger (daily deficit punishment) from food noise (persistent background awareness proportional to distance below defended weight). Binges are a crude proxy for food noise. Better tests:

1. **Distance → continuous intake: -30 cal/day per kg below set point.** Direction matches the essay's claim (~100 cal/kg) but magnitude is 3× weaker. Pre-tirzepatide, intake is nearly flat across ±5 lbs of set point (2156-2177 cal).

2. **Restriction duration → rebound: r=-0.065.** The breathing analogy (longer restriction → bigger gasp) is not supported. 3-5 day and 8-14 day restriction runs produce identical rebounds (~2000 cal).

3. **Intake VARIANCE is where food noise lives.** Below set point: std=615, CV=0.278. At set point: std=528, CV=0.244. **16% more erratic eating when below defended weight.** Tirzepatide reduces CV from 0.24-0.28 to 0.19-0.20 — a **25-30% reduction in intake variability.** The drug doesn't just lower mean intake; it makes eating dramatically more consistent. Food noise manifests as chaos, not as a constant upward pressure.

4. **Post-restriction binges halved by tirzepatide.** 26% binge rate in the 7 days after ending a restriction run → 13% on drug.
Command: `python analysis/D_food_noise_variance.py`

**Revised model:** The essay's core insight — that food noise is distinct from hunger — is confirmed. But food noise is not proportional to set point distance (the original prediction). It manifests as **intake variance** (erratic eating, behavioral momentum, clustering binges) rather than a constant caloric drive. The set point defends weight through expenditure adaptation (TDEE/RMR ratio); food noise is a separate system operating through behavioral turbulence. Tirzepatide silences the turbulence (CV -25%), dampens the variance, and halves post-restriction binges. If food noise is a generalized resource-acquisition drive (same circuit as alcohol, gambling — all silenced by GLP-1 agonists), then the trigger is opportunity and behavioral momentum, not the metabolic gap itself.

**Week-scale intake invariance — REJECTED.** [Preregistered claim](https://bsky.app/profile/lucent.substack.com): week-scale food intake is uncorrelated with emotional state or individual meal interventions. A week where you "cut back" ends up within a few hundred calories of a week where you feasted then compensated without noticing. Test: compare weekly calorie variance to daily variance. If weekly variance is dramatically lower than you'd expect from independent daily draws, there's a homeostatic weekly regulator.

**Result:** Weekly total calorie std is 2383 cal/week. If days were independent, expected std = 550 × √7 = 1456. Ratio = **1.64** — weekly variance is 64% *higher* than independent draws. Binges cluster rather than compensate. This is anti-homeostatic: high days pull neighboring days up. The intake autocorrelation structure confirms this: r=0.40 at 1-day lag, 0.22 at 7 days, 0.17 at 14 days, dropping to 0.04 by 30 days. Eating patterns persist for ~2 weeks then reset — behavioral cycles, not a homeostatic regulator.
Command: `python analysis/E_weekly_invariance.py`

**Tirzepatide dose-response quantification — PARTIALLY ANSWERED.** "Very curious to learn how many mg correspond to kcal subtracted from daily satiety point." The [study 128-OR](https://doi.org/10.2337/db23-128-OR) showed 15mg lowered intake by 900 kcal. At 2.5mg, observed reduction was 200-400 kcal. Firsthand: "each mg shaves 100 kcal/day from satiety." With dose escalation data (2.5→5→7.5→10→12.5mg) and daily intake logged, fit the dose-response curve. Key observation: at 12.5mg, weight stabilized at 202 then bounced to 205 over Christmas and the drug "started working again" at the higher weight — as if each dose corresponds to a specific set point or loss delta, not a fixed kcal reduction.

**Results from pharmacokinetic modeling (2026-03-24):** Modeled blood concentration using FDA PK parameters (t½=5d, Tmax=24h) with superposition of all prior injection curves. Key findings:
- **Blood level → daily intake: r = -0.50** (partial, controlling time trend). The strongest single predictor of daily calorie intake in the entire dataset.
- **Weekly sawtooth**: injection day 1652 cal → trough (day 5) 2220 cal. A 568 cal/day appetite swing directly tracking the PK curve, confirming that appetite suppression is pharmacokinetically mediated, not psychological.
- **Intake model**: `calories = 2345 - 49 × effective_blood_level`. At zero drug: 2345 cal/day. Fresh 12.5mg peak: 1504 cal/day. This gives ~49 cal reduction per arbitrary blood level unit, or roughly **35 cal/day per mg at steady state** (lower than the self-estimated 100 cal/mg, likely because the self-estimate compared peak effect vs pre-drug baseline rather than steady-state average).
- **Tachyphylaxis**: dose effectiveness decays with half-life of 32 weeks. After 20 weeks on the same dose, 65% effective. This explains the firsthand observation of the drug "starting to work again" after a weight plateau — the plateau is tolerance building, and dose escalation resets it partially.
- **Overall reduction**: 456 cal/day (18.6%) from the pre-tirzepatide year. Study 128-OR at 15mg showed 900 kcal reduction; this data at 12.5mg steady-state shows ~534 cal reduction, scaling linearly.
- **Christmas 2025 spike visible**: weeks 16-18 at 12.5mg show 2672, 2793 cal — environmental override of pharmacological suppression, consistent with food noise being a drive that can be overcome by sufficiently strong external cues.
Command: `python analysis/F_tirzepatide_pk.py`
Artifact: `drugs/tirzepatide.csv`

**New transition-dynamics result (2026-03-25):** The drug's strongest effect is on **state transitions**, not just mean intake. Daily intake states were classified as `restriction` (<1800), `typical` (1800-2399), `high` (2400-2799), `binge` (2800+). Pre-tirzepatide, a `high` day escalates to a `binge` the next day 19.6% of the time; on tirzepatide, 4.8%. A `binge` day is followed by another `binge` 31.6% of the time pre-drug and 0% on drug. This fits the "food noise as behavioral turbulence" model better than a simple mean-shift model: tirzepatide mainly breaks escalation and persistence.
Command: `python analysis/G_tirzepatide_dynamics.py`
Artifact: `analysis/G_tirzepatide_transition_summary.csv`, `analysis/G_tirzepatide_rebound_summary.csv`

**Trajectory phase changes the binge landscape, but adds little predictive power once recent intake is known.** In the pre-tirzepatide era, binge rate is 8.9% in `falling` phase, 12.2% in `stable`, 15.1% in `rising`. At similar set-point distances, `rising` remains the most binge-prone branch. But in walk-forward prediction, adding phase to distance does not beat distance alone (AUC 0.525 vs 0.537), and adding phase to yesterday's calories slightly hurts (0.721 vs 0.724). Phase seems to organize the background landscape more than it improves day-ahead prediction.
Command: `python analysis/H_phase_binge_landscape.py`
Artifact: `analysis/H_phase_binge_summary.csv`, `analysis/H_phase_binge_distance_bins.csv`, `analysis/H_phase_binge_model_auc.csv`

Remaining questions: Does each dose correspond to a specific set point, or does the tachyphylaxis model fully explain the plateaus? The 32-week half-life of effectiveness means the drug never fully stops working — it just asymptotically approaches a reduced effect. At 12.5mg after 1 year, effective level would be ~36% of initial. Is this consistent with the observed weight trajectory?

**Plateau dynamics: resume vs. jump-to-catch-up.** Two models of weight loss on GLP-1s: (1) weight drops, plateaus, then resumes linear loss from the plateau; (2) weight drops, plateaus, then jumps down as if the plateau never happened, resuming the original trajectory. GLP-1 Discord chose model 1, but firsthand experience suggests model 2. The glycogen-water smoothing should distinguish these — if plateaus are water/glycogen masking ongoing fat loss, model 2 is correct and the "jump" is just water finally releasing.

**Cold intolerance as protein-mediated.** "Couple days of low protein and already feeling hypoglycemic and much less cold intolerance." Cold resolved by first bench presses of a workout, corroborating that protein should be eaten before workouts (mTOR activation). The temperature data + daily protein intake can test whether protein intake predicts next-day body temperature.

**Set point ratchet asymmetry.** "A month of consistent overeating ratchets your set point up a couple pounds that never, ever come off no matter how hard you fight." Test: find periods of sustained overeating (>2800 cal/day for >5 consecutive days) and measure whether the subsequent weight floor permanently increases. Is the ratchet truly one-way, or can sustained undereating ratchet it down? "Weight gain not being reversible by simply removing the element that caused it would be a devastating violation of an expected symmetry."

**Dead zone — REFRAMED AS PHASE-DEPENDENT.** "I have a dead zone between 2000-2500 kcal where my metabolism adapts to maintain." Pre-tirzepatide data do not support a single universal dead zone. Instead, the intake→TDEE relationship depends on trajectory phase. In `stable` phase, higher intake weakly raises TDEE (coef `+0.018` cal TDEE per intake cal). In `rising`, the slope is near zero (`+0.005`). In `falling`, it is slightly negative overall (`-0.018`), with TDEE nearly flat from 1600 through 3000 kcal. The more accurate model is not one dead zone for all states, but **branch-specific expenditure compression**.
Command: `python analysis/I_deadzone_phase.py`
Artifact: `analysis/I_deadzone_phase_bin_summary.csv`, `analysis/I_deadzone_phase_regression.csv`

**Restriction archetypes — PARTIALLY ANSWERED.** Pre-tirzepatide restriction runs (`<1800 kcal` for `>=3` days) leave an average post-run TDEE/RMR penalty of `-0.0063` versus the 30 days before the run (`-19 cal/day`). But the penalty varies sharply by numeric archetype. `Long` runs (`>=6` days) recover best: post-pre ratio `+0.0038`. `Low-carb` runs (`<170g/day`) also recover relatively well: `+0.0016`. `Low-protein` runs (`<58g/day`) do worst: `-0.0119`. `High-step` cuts (`>=4200 steps/day`) show the largest penalty: `-0.0177` (`-41 cal/day`). So the old data does not support “restriction is restriction.” The hysteresis cost depends on how the restriction is structured.
Command: `python analysis/J_restriction_archetypes.py`, `python analysis/K_tdee_hysteresis.py`
Artifact: `analysis/J_restriction_archetype_summary.csv`, `analysis/K_tdee_hysteresis_phase_summary.csv`
Artifact: `analysis/J_restriction_runs.csv`, `analysis/J_restriction_archetype_summary.csv`, `analysis/J_restriction_archetype_regression.csv`

**The restriction-archetype penalties are not mainly rebound-eating effects.** The worst post-run expenditure penalties do **not** line up with the biggest 7-day rebound. `Low-protein` runs have a larger penalty (`-0.0119`) despite *lower* next-7-day calories (1934) and fewer binges (20.6%) than the higher-recovery non-low-protein runs. `High-step` cuts are worst of all (`-0.0177`) with only middling rebound. In a regression of post-run penalty on archetype flags plus next-7-day mean calories and any binge, the rebound terms are near zero while `long_run` stays strongly positive (`+0.0153`), `low_carb` stays positive (`+0.0068`), and `high_steps` stays strongly negative (`-0.0187`). That points to a real metabolic branch effect, not just post-diet overeating.
Command: `python analysis/L_restriction_rebound.py`
Artifact: `analysis/restriction_rebound_summary.csv`, `analysis/restriction_rebound_regression.csv`

**The bad archetypes fail in different ways.** `High-step` cuts are bad **during** the run and after it: they shift TDEE/RMR downward during the run itself (`run-pre -0.0065`) and keep falling afterward (`post-run -0.0112`). `Low-protein` cuts look different: the during-run shift is near zero (`-0.0015`), but recovery after the run is poor (`post-run -0.0104`). By contrast, `long` and `low-carb` runs start on a higher branch (`run-pre +0.0076`, `+0.0042`) and then lose relatively little afterward. So the mechanism is not one generic “diet penalty”: some cuts trigger an immediate downward branch shift, others mainly impair recovery.
Command: `python analysis/M_restriction_branch_shift.py`
Artifact: `analysis/restriction_branch_shift_summary.csv`, `analysis/restriction_branch_shift_regression.csv`

**Antihistamine contribution.** [H1 antihistamines associated with obesity: 10 kg difference (NHANES)](https://doi.org/10.1038/oby.2010.176).

Timeline:
- **2011–2024-12**: Daily 10mg cetirizine (Zyrtec) mornings. Entire dataset under this influence.
- **~2024-12-17**: Started weaning. Switched to 5mg cetirizine at night.
- **2024-12-31**: "Been on 5mg of Zyrtec at night instead of 10mg in the morning for 2 weeks now and my quality of life is down about 25% from general malaise, but the 10 kg difference in this study is incredible."
- **2025-01-25**: "After weening down, if I had to estimate the difference between 10mg cetirizine in the morning and 5mg at night, it is almost an additional 2mg of serum tirzepatide."
- **2025-04-06**: Amazon purchase of levocetirizine (generic Xyzal). Switch from cetirizine to levocetirizine complete.
- **2025-06-18**: Costco Xyzal 110ct restock ($37.99).
- **2025-07-07**: Confirmed nightly routine: "2.5mg (halved) levocetirizine."

Net change: 5mg active enantiomer mornings → 2.5mg active enantiomer nightly.

Confounded by tirzepatide dose escalation (7.5→10mg ~2025-05-27) in the same period. Open question from Daniel Quinn/Ishmael framing: does Zyrtec "artificially" increase appetite, or does allergy misery decrease it? Which is baseline?

**Flavorless oil calorie undercounting.** "Maybe we have a running counter that uses flavor (to estimate macros) × chew/crunch/swallow to tally calories and ancestrally absent flavorless oil undershoots the calculation." If cooking oil is systematically more fattening per calorie than flavorful foods, days heavy in oil-cooked foods should show more weight gain than their calorie count predicts. Testable by classifying food items as oil-heavy vs not.

**Advanced glycation end-products (AGEs).** A reductionist-friendly candidate for quantifying food "badness" beyond NOVA/UPF classification. "Seems way more likely than hyperpalatability" since GLP-1 users aren't especially tempted by hyperpalatable foods. Would require classifying food items by AGE content.

**Breakfast/lunch front-loading — TESTED, MIXED.** Does eating over 1000 calories for breakfast and lunch combined predict a lower total for the day?

**Result:** Morning *percentage* of daily calories correlates with lower total (r=-0.19, 437 cal range from lowest to highest bin). But morning *absolute* calories correlate with HIGHER total (r=+0.48). Eating 1000+ cal before dinner → 2412 cal/day, vs <1000 → 2024. The percentage is circular — it inversely reflects the total. Front-loading in absolute terms does not reduce daily intake.

**Fiber satiety — WEAK SIGNAL.** Fiber consumed in breakfast+lunch inversely correlated with daily intake? Morning fiber controlling for morning calories: r=-0.094 partial. Directionally correct — at the same morning calorie level, more fiber predicts slightly lower daily total (~100 cal at the extreme). Weak but potentially actionable.
Command: `python analysis/N_dietary_predictors.py` (meal timing and fiber sections)

**Protein leverage — WITHIN-DAY ONLY.** Same-day protein % vs total intake: r=-0.34. The 4-8% protein bin eats 2373 cal/day vs 21-30% at 1783 cal — a 590 cal range. But next-day prediction controlling for today's calories: r=-0.04. Protein leverage regulates within each day (low protein meals → more eating same day) but does not carry over to the next day. Tirzepatide weakens the effect (r drops from -0.19 to -0.09 on drug), suggesting protein leverage acts through the same appetite circuits the drug suppresses.
Command: `python analysis/N_dietary_predictors.py`

**Diet experiment detection — CURATED EPOCHS DO SEPARATE REAL REGIMES.** The manually annotated epochs in `intake/diet_epochs.csv` are not decorative; they map onto clear numeric regimes. `Travel_binge` is the cleanest high-intake family (`2716` cal/day, binge rate `32.1%`, `+10.4` lbs fat across the two epochs). `Keto_phase` is distinct mainly by high protein and low binge frequency (`2163` cal/day, `110.5g` protein, binge rate `2.1%`), not by unusually high TDEE/RMR (`1.042`). `Potato_diet` is the strongest monotony-cut family: across 4 epochs / 69 days it averages `1930` cal/day, `48.6g` protein, `0%` binge days, and `-7.85` lbs fat. So the epoch annotations are good enough to use as regime labels without rereading food names.
Command: `python analysis/O_diet_epoch_analysis.py`
Artifact: [analysis/diet_epoch_summary.csv](/home/lucent/health-data/analysis/diet_epoch_summary.csv), [analysis/diet_epoch_family_summary.csv](/home/lucent/health-data/analysis/diet_epoch_family_summary.csv)

**Potato diets look like low-noise restriction, not a magic metabolic branch.** Relative to the 14 days before each potato attempt, the pooled potato windows cut calories from `2519` to `1877` and protein from `74.9g` to `38.8g`, while binge rate falls from `14.3%` to `0%` and average fat mass drops `1.96` lbs over the epoch windows. But TDEE/RMR is basically unchanged during the potato windows (`1.168` before vs `1.159` during), then lower afterward (`1.098` in the 14-day post windows) as calories and binges rebound (`2629` cal/day, binge rate `19.6%`). The useful interpretation is: potatoes behaved like a highly monotone binge-suppressing cut, not a uniquely expenditure-preserving intervention.
Command: `python analysis/O_diet_epoch_analysis.py`
Artifact: [analysis/potato_epoch_window_summary.csv](/home/lucent/health-data/analysis/potato_epoch_window_summary.csv), [analysis/potato_epoch_contrast.csv](/home/lucent/health-data/analysis/potato_epoch_contrast.csv)

**Jordan trips spike calories without matching TDEE.** The two month-long Jordan stays are annotated in `travel/trips.md` and show the same pattern twice: intake jumps to `3262`/`2529` cal/day while the Kalman TDEE estimate barely budges (`2071`/`2021` kcal, which is `-35`/`-21` below the 28 days before each trip) and the TDEE/RMR ratio drifts downward. No compensatory TDEE surge appears, so the trips end up with big fat gains (`+10.3` lbs in 2015, `+4.4` lbs in 2019) even though the high post-trip intake quickly eases back to regular levels. They are textbook jet-lagged binge windows without a temporary "energy out" bonus.
Command: `python analysis/P_jordan_trip_analysis.py`
Artifact: [analysis/jordan_trip_summary.csv](/home/lucent/health-data/analysis/jordan_trip_summary.csv), [analysis/jordan_trip_delta.csv](/home/lucent/health-data/analysis/jordan_trip_delta.csv)

**Latent control capacity is detectable, but only modestly.** A simple hidden-state proxy for "willpower" can be built from short-memory restriction pressure, monotony, low-calorie days, and high steps, with partial recovery from easier refeed days. In the pre-tirzepatide era from 2017-01-01 onward, that latent control-debt score reaches AUC `0.696` for next-day binge prediction, nearly matching yesterday's calories alone at `0.704`. Combining the two improves to `0.7174`, so the latent model adds some real information, but not enough to call it a dominant hidden driver. The strongest version is fast-moving rather than slow: best half-life is `2` days, which suggests short-term control depletion is more defensible than a long-burn "ego depletion" story.
Command: `python analysis/Q_latent_control_capacity.py`
Artifact: [analysis/latent_control_model_summary.csv](/home/lucent/health-data/analysis/latent_control_model_summary.csv), [analysis/latent_control_feature_search.csv](/home/lucent/health-data/analysis/latent_control_feature_search.csv), [analysis/latent_control_daily.csv](/home/lucent/health-data/analysis/latent_control_daily.csv)

**Metabolic failure is better predicted by branch position than by a hidden willpower story.** Using 133 pre-tirzepatide restriction runs and taking `TDEE/RMR` collapse as the endpoint, the strongest broad predictor of both immediate and net metabolic failure is the phase you enter the cut from. For `run_minus_pre_ratio`, the best single predictor is `phase_code` with leave-one-out `R^2=0.167`; for `post_minus_pre_ratio`, it is again `phase_code` at `R^2=0.197`. Adding `mean_steps` gives the best 2-feature models (`R^2=0.200` and `0.239`). The concrete pattern is stark: falling-phase high-step cuts are the worst (`post-pre -0.0598`), while rising-phase cuts are metabolically resilient or even favorable (`+0.0305` low-step, `+0.0122` high-step). So for the metabolic endpoint, branch state and activity load matter more than a generic latent-control narrative.
Command: `python analysis/R_metabolic_failure_predictors.py`
Artifact: [analysis/metabolic_failure_feature_search.csv](/home/lucent/health-data/analysis/metabolic_failure_feature_search.csv), [analysis/metabolic_failure_phase_summary.csv](/home/lucent/health-data/analysis/metabolic_failure_phase_summary.csv), [analysis/metabolic_failure_runs.csv](/home/lucent/health-data/analysis/metabolic_failure_runs.csv)

**The restriction-phase effect survives matching on starting fat mass.** Within matched start-fat bands, falling-phase cuts still underperform rising-phase cuts: at `35-50` lbs FM, post-pre `TDEE/RMR` is `-0.0279` in falling vs `+0.0006` in rising; at `50-65`, `-0.0480` vs `+0.0340`; at `65-80`, `-0.0320` vs `+0.0448`. In a simple regression controlling for start fat mass directly, `falling_vs_rising` remains strongly negative: `-0.0295` for during-run branch shift, `-0.0534` for net post-run failure, and `-0.0240` for recovery. So the branch effect is not just a disguised fat-mass effect; it looks like real path memory.
Command: `python analysis/S_metabolic_failure_matched.py`
Artifact: [analysis/metabolic_failure_matched_bands.csv](/home/lucent/health-data/analysis/metabolic_failure_matched_bands.csv), [analysis/metabolic_failure_matched_regression.csv](/home/lucent/health-data/analysis/metabolic_failure_matched_regression.csv)

**The branch-memory effect also survives matching on pre-run metabolic state.** Matching each falling-phase restriction run to the nearest rising-phase run in the 2D space of `fat_mass_start` and pre-run `TDEE/RMR` (`pre_ratio`) still leaves a large penalty for falling runs: across all 43 matched pairs, falling minus rising is `-0.0327` for during-run shift and `-0.0442` for net post-run failure; in the tighter matches (`distance <= 0.50`, `n=34`), the net penalty is still `-0.0462`. A regression controlling for both start fat mass and `pre_ratio` also keeps the phase term negative: `falling_vs_rising = -0.0281` for `run_minus_pre_ratio`, `-0.0398` for `post_minus_pre_ratio`, and `-0.0116` for `post_minus_run_ratio`. So the falling-branch disadvantage is not explained away by starting the run fatter or already more metabolically adapted.
Command: `python analysis/T_metabolic_failure_state_matched.py`
Artifact: [analysis/metabolic_failure_state_match_summary.csv](/home/lucent/health-data/analysis/metabolic_failure_state_match_summary.csv), [analysis/metabolic_failure_state_match_pairs.csv](/home/lucent/health-data/analysis/metabolic_failure_state_match_pairs.csv), [analysis/metabolic_failure_state_match_regression.csv](/home/lucent/health-data/analysis/metabolic_failure_state_match_regression.csv)

**Calorie misestimation detection.** Create a database of distinct foods eaten (there aren't many) and see if scaling their calories up or down produces a better fit with weight. A frequent restaurant meal that's consistently underestimated would show up. More interestingly, foods that appear super-fattening or near-zero-calorie beyond their label — like pistachios were once predicted to be.

## Ready now (+ steps.csv and sleep.csv)

**Gravitostat — WEAK, WRONG DIRECTION.** Compute daily foot-pounds exerted (current weight × steps) and test if this predicts next-day intake or weight change. Firsthand report: "+16 lb weighted vest 3 mile walk kills hunger 12-36 hours after so I can easily cut calories" — but like everything, "stop and it all comes back." Specific hypothesis: the top end of the set point control system discards all excess calories if >x ft-lbs (~8000 steps?) are exerted on legs, and this is the mechanism behind "no weight gain while traveling Europe." The [gravitostat study](https://doi.org/10.1016/j.eclinm.2020.100338) supports this. Impressive if the model can detect the weighted vest period. Concern: seasonal confounding (more steps in summer, different eating too).

**Result:** Foot-pounds → next-day intake: r=0.11 raw, **r=0.055 partial** (controlling today's intake). Positive — more steps predict *more* eating the next day, the opposite of the gravitostat prediction. The signal is weak and likely confounded by activity patterns (active days = higher overall eating). Steps → next-day weight change: r=0.014 partial. No detectable signal at any timescale (daily, weekly, monthly). Steps were also tested as a TDEE covariate: literature-rate activity correction (0.5 kcal/kg/km) produces only 0.6% variance reduction in derived TDEE. At median 3,400 steps/day, the activity signal (~150 cal) is buried in ~1 lb/day of non-caloric weight noise. The gravitostat may operate on longer timescales or require a threshold (the "8000 steps" hypothesis) that can't be tested with this step count distribution (median 3,400, 75th percentile 6,372).
Command: `python analysis/N_dietary_predictors.py` (gravitostat section)

**Refined step result: acute walking is not protective, but sustained walking load looks compensated and branch-dependent.** Using 7-day average steps and foot-pounds instead of a one-day lag, higher recent walking predicts **more** future intake (`steps_7d -> future7 calories beta = +0.103`, `future14 calories = +0.075`; foot-pounds slightly stronger at `+0.116` and `+0.093`) and a tiny reduction in future `TDEE/RMR` overall (`-0.0039` for `steps_7d -> future14 ratio`, `-0.0051` for foot-pounds). The interesting part is the phase split: in the `falling` branch, higher recent steps predict a stronger future ratio drop (`beta = -0.0478`), while in `stable` the effect is nearly zero and in `rising` it is slightly positive. Thresholded summaries show the same thing: in falling phase, the top 20% of 7-day steps (`>=6311/day`) is followed by **lower** 14-day intake (`1974` vs `2088`) and **lower** 14-day `TDEE/RMR` (`1.0918` vs `1.1569`). That is much closer to a constrained-energy / compensated-activity model than to a simple gravitostat benefit.
Command: `python analysis/U_steps_compensation.py`
Artifact: [analysis/steps_compensation_regression.csv](/home/lucent/health-data/analysis/steps_compensation_regression.csv), [analysis/steps_compensation_phase_thresholds.csv](/home/lucent/health-data/analysis/steps_compensation_phase_thresholds.csv)

**Exercise sessions add a better movement variable than blunt daily steps.** The Samsung exercise export is usable enough to isolate deliberate walks. Type `1001` is almost certainly walking from the session geometry: `1,636` sessions, mean `26.63` min, mean `1,798` m, mean `2,507` steps, with `324` recent sessions in 2024+. Type `1002` is run-like (`31.49` min, `2,877` m, `4,089` steps, `10.15` cal/min). More importantly, there are `99` post-2020 spring/summer `1-5pm` double-walk days matching the pattern you described: two `20-45` min walking sessions separated by a `5-40` min break. Those days have median pair totals of `67.32` min, `5.32` km, and `6,688` session steps, with median full-day steps `8,204`.

These deliberate walk-pair days do **not** behave like generic high-step days. Compared with other `>=5800` step days, they show lower same-day intake (`2066` vs `2308` cal), lower 14-day future intake (`2174` vs `2235`), and slightly higher future `TDEE/RMR` (`1.0994` vs `1.0855`). Even after nearest-neighbor matching on total daily steps (mean step difference `4.86`), paired-walk days still show `-27.89` cal lower future 14-day intake and `+0.0236` higher future `TDEE/RMR`. So the first pass suggests these deliberate daylight walks are a meaningfully different regime than incidental high-step load.
Command: `python analysis/V_exercise_walk_analysis.py`
Artifact: [analysis/exercise_type_summary.csv](/home/lucent/health-data/analysis/exercise_type_summary.csv), [analysis/daylight_walk_pair_days.csv](/home/lucent/health-data/analysis/daylight_walk_pair_days.csv), [analysis/daylight_walk_pair_summary.csv](/home/lucent/health-data/analysis/daylight_walk_pair_summary.csv), [analysis/daylight_walk_pair_contrast.csv](/home/lucent/health-data/analysis/daylight_walk_pair_contrast.csv)

**Steps as set point shifter.** Are there conditions where restriction without steps produced no loss, but restriction with >5000 steps did?

**Sunlight-intake correlation.** Use sleep/wake times to count hours of overlap with local sunrise/sunset (Knoxville coordinates, historical weather data) and see if that corresponds to intake. Very little sunlight exposure overall, making any effect easier to detect.

## Needs MFP API enrichment (iron, potassium, vitamin D, saturated fat, added sugars, PUFA/MUFA)

**Omega-6:3 ratio and RMR.** The contamination-adjacent theory. If omega-6 to 3 ratio exceeds some threshold, does apparent RMR drop? Requires fat subtype breakdown, not available from the printable diary. The top ~200 foods cover ~80% of intake and can be enriched via MFP API with serving-size calorie checksums to verify matches.

**Nutrient-specific set point control.** Can a control system predict RMR or next day's intake from micronutrient ratios? Example (pure whimsy): if trailing 30-day omega-6:3 > 10:1 and yesterday's protein > 50g, lower RMR by 500 to store all fat consumed, else raise RMR to maintain set point.

**NOVA classification and binges.** Correlate ultra-processed food consumption with subsequent binges. Requires classifying foods by processing level, feasible from the food names.

## Needs additional data export

**Body temperature during restriction — PARTIALLY ANSWERED.** Background notes body temp below 97° when consuming below RMR. Withings data: 1,315 readings, 364 days, Dec 2023 – Mar 2026.

**Results (2026-03-25):** Pre-tirzepatide (9 months), 14-day trailing intake → temperature: r=+0.22. Eating more → warmer, confirming the documented cold intolerance during restriction. The 2600-3200 cal/day bin is 0.10°F warmer than 2000-2300 cal. But temperature only weakly tracks Kalman TDEE (r=0.04) — it's not a substitute for calorimetry.

**Tirzepatide does NOT override the thermostat defense — CORRECTED by calorimetry.** An earlier analysis found temperature correlating with tirzepatide blood level at r=+0.45 (5mg → 97.40°F, 12.5mg → 97.75°F), suggesting the drug raised body temperature despite restriction. This was confounded by dose escalation timing (low doses in fall/winter, high doses in spring/summer).

**Direct calorimetry disproves the override hypothesis.** Three new Cosmed Fitmate measurements from the tirzepatide era:
- 2024-09-17 (tirz day 0, FM=83 lbs): RMR=1930 cal/day, temp=97.7°F
- 2024-09-19 (tirz day 2, FM=83 lbs): RMR=1836 cal/day, temp=97.7°F
- 2025-09-03 (tirz 1 year, FM=68 lbs): RMR=1750 cal/day, temp=97.5°F

RMR dropped 180 cal/day over one year. Temperature dropped 0.20°F. The composition-aware model (RMR = 20.5×FFM_kg + 4.6×FM_kg + 500, fitted to 25 measurements) predicts 1956 cal at the 2025 body composition; actual is 1750 — a **206 cal/day metabolic adaptation** beyond tissue loss (~10% suppression). The body IS conserving on tirzepatide. The drug reduces intake by ~450 cal/day; the body claws back ~200 through reduced RMR, leaving a net deficit of ~250 cal/day driving the observed fat loss rate.
Command: `python analysis/P2_rmr_model.py`
Artifact: `RMR/rmr.csv`, `analysis/P2_daily_composition.csv`

**Temperature by trajectory phase — PARTIALLY IDENTIFIABLE ONLY.** The overlap is structurally limited: temperature starts Dec 2023, so pre-tirzepatide coverage is mostly `stable` plus a small `rising` tail at 75-85 lbs fat mass, while the on-drug era is dominated by `falling` phase during weight loss. The temperature→blood_level correlation (r=+0.45) is now known to be seasonally confounded. Temperature tracks restriction (r=+0.22 with trailing 14-day intake) but does not independently track the drug's effect once RMR is measured directly.
Command: `python analysis/X_temperature_phase.py`
Artifact: `analysis/X_temperature_daily.csv`, `analysis/X_temperature_phase_overlap.csv`, `analysis/X_temperature_phase_band_summary.csv`, `analysis/X_temperature_phase_regression.csv`

**Expenditure-side set point — confirmed by calorimetry.** The set point defense operates through TDEE, not intake. During restriction, TDEE/RMR ratio drops to 1.02 and temperature drops below 97°F — the body conserves. Tirzepatide does NOT override this: RMR drops 206 cal below composition-predicted on the drug, confirming the metabolic defense is active. The drug's advantage is purely on the intake side (silencing food noise, breaking binge persistence), not the expenditure side.

**RMR now has 25 measurements** (3 lab Cosmed + 22 home Cosmed Fitmate, 2011-2025) extracted from the device's DBF database. The fitted model (RMR = 20.5×FFM_kg + 4.6×FM_kg + 500) is close to Cunningham (22, 0, 500) and physiologically stable. RMSE=170 cal/day, at the instrument's noise floor.

**Tirzepatide as set point intervention — PARTIALLY ANSWERED.** The medication log (Sep 2024 onward) combined with the 13-year pre-intervention baseline.

Pre-tirzepatide binge rate (>2800 cal): 12.2% of days. On tirzepatide: 3.4%. A 3.6× reduction. Blood level predicts daily intake at r=-0.50 with a 568 cal/day weekly sawtooth (injection day 1652 cal, trough day 2220 cal). Dose tolerance builds with 32-week half-life.

**The set point operates through expenditure, not intake.** Analysis of 5,429 days shows:
- Intake autocorrelation dies at 30 days (no persistent homeostatic signal)
- Weekly intake anti-compensates (binges cluster rather than correct — ratio 1.64 vs expected <1.0)
- Recent weight gain predicts higher future intake (r=+0.26), not compensatory restriction
- Binge rate peaks at FM=50 lbs (~195 lbs), not at the lowest weight
- But TDEE/RMR ratio clearly shows metabolic adaptation: 1.02 during sustained restriction, 1.14 when restriction eases

The body doesn't defend its weight by making you eat more. It defends by burning less. Food noise (the thing tirzepatide silences) is a separate, non-homeostatic drive — a resource-acquisition impulse that compounds rather than compensates. The drug attacks both: silencing food noise (intake r=-0.50) AND overriding the metabolic defense (temperature up despite restriction).

**Circadian misalignment.** Consistent 3am-11am sleep with very little sunlight exposure means meals are shifted ~4 hours relative to the solar cortisol rhythm, with essentially no zeitgeber correction. Meal timing relative to circadian phase affects metabolic rate and fat storage independently of calories. This is so constant in the data it's invisible as a variable but extremely unusual compared to study populations. There was an early bedtime experiment — check the sleep data for periods of deviation from the 3am-11am pattern and correlate with weight trajectory. The [housekeepers study](https://doi.org/10.1111/j.1467-9280.2007.01867.x) (tell people their activity is exercise and they lose weight) suggests belief/awareness matters too.

**Hunger taxonomy.** "We have Sapir-Whorf'd ourselves out of understanding the varied nature of hunger, giving it one word while over-granulating 'cravings' by food." At minimum three distinct signals: (1) hunger-the-stick (caloric deficit punishment, same-day), (2) food noise (resource acquisition drive, set-point-distance-proportional), (3) protein-ravenous ("the immediate craving after a protein meal that goes away quickly if you wait"). GLP-1s prove hunger is not purely mechanical stomach-stretch signaling. These may be separately quantifiable from the intake data — days where protein was high but calories low should show different next-day behavior than days where both were low.

**Gut microbiome as missing variable.** Heavy childhood erythromycin (macrolide that reshapes gut microbiome) combined with 35 years of no animal protein recolonization source is a textbook setup for persistent dysbiosis. Not testable from current data, but a microbiome assay (16S or shotgun sequencing) could reveal whether composition is atypical and whether it changes in response to dietary interventions already in the data.

**Low-protein / BCAA restriction.** Tested a low-protein period targeting foods with <2g isoleucine per 2000 kcal. Reported less cold intolerance within days. The intake data should show this period clearly (protein dropping to unusually low levels). Correlate with body temperature and weight trajectory. Related to the protein leverage hypothesis (20% protein diets = 10% protein + 500 cal).

**GLP-1 weekly cycle as mood/intake predictor — CONFIRMED.** "Interesting it's not obvious what my weekly shot day is given the buzzing high at peak and depressive low on day 6." With weekly injection dates from the medicine log and daily intake, test whether intake varies by day-of-injection-cycle. Day 1-2 post-injection should show lowest intake, day 6-7 highest. The subjective strength ratings already capture this partially.

**Result**: The pharmacokinetic blood level model confirms this precisely. By day of injection cycle: day 0 = 1652 cal (blood level 11.7), day 1 = 1741 (11.4), day 2 = 2015 (9.9), day 3 = 2148 (8.6), day 4 = 2139 (7.5), day 5 = 2220 (6.5). The 568 cal swing from peak to trough tracks the PK curve almost exactly — appetite is pharmacokinetically mediated at the daily timescale.
Command: `python analysis/F_tirzepatide_pk.py` (weekly cycle section)

## Modeling approach

Slime Mold recommends [modeling as a control system](https://slimemoldtimemold.com/2022/03/15/control-and-correlation/) where control input has no correlation with modeled variables, but is instead a homeostatic mechanism. A day's macros or total consumed may have no correlation with previous or next days or any other available variable, but instead be correlated with a hidden variable, set point, that meanders slowly up or down. Much like RMR, a secret set point could be determined from the data and binges occur when below the set point.

The [food noise essay](https://lucent.substack.com/p/craving-food-noise) refines this: the set point is defended by two distinct mechanisms. Hunger (the stick) punishes today's shortfall — it's specific in timescale and proportional to deficit. Food noise (the carrot) rewards resource acquisition — it's proportional to distance below set point and compounds over time like an unclaimed package. The control system has two feedback loops, not one. GLP-1 agonists silence food noise but may not affect hunger, which is why the tirzepatide data can distinguish them.

## Cross-cutting findings (2026-03-25)

### The set point operates through expenditure, not intake

Five tests of set point behavior on the intake side all come back negative:

1. **Recent weight change vs future intake: r=+0.26** (positive, not negative). Weight gain predicts *more* eating, not compensatory restriction.
2. **Weekly intake anti-compensates**: ratio 1.64. Binges cluster rather than correct.
3. **Intake autocorrelation dies at 30 days**: r=0.40→0.22→0.17→0.04 over 1→7→14→30 day lags. No persistent homeostatic signal.
4. **Binge rate peaks at FM 50 (17.1%), not at minimum weight** (FM 20: 10.4%). Not monotonic with distance from any putative set point.
5. **Binge rate drops 3.6× on tirzepatide** (12.2% → 3.4%), confirming binges are pharmacologically mediated, not homeostatic.
Command: `python analysis/Y_set_point_intake_tests.py`

The expenditure side tells a different story: TDEE/RMR ratio clearly shows metabolic adaptation (1.02 during sustained restriction, 1.14 when restriction eases). The body defends weight by burning less, not by making you eat more. Food noise is a separate, non-homeostatic drive that compounds rather than compensates.

### Tirzepatide attacks only the intake side

The drug silences food noise (intake r=-0.50, binge persistence 31.6% → 0%). But direct calorimetry confirms it does NOT override the metabolic defense: RMR dropped from 1930 to 1750 cal/day over one year on the drug, 206 cal below composition-predicted — the body is conserving at the same rate it would during any restriction. The drug's advantage is purely behavioral: reduced intake and eliminated binge cascades, not increased expenditure.

**Lean mass preserved.** 70 composition measurements show FM 83.6→60.0 lbs (-23.6) with lean mass 141.1→143.0 lbs (+1.9). Only ~7% of weight loss was lean tissue, vs typical 20-25% during caloric restriction. The drug produces nearly pure fat loss.
Command: `python analysis/F_tirzepatide_pk.py` (lean mass section)

### What predicts daily weight change (after glycogen + sodium correction)

Tested every available variable as a predictor of next-day smoothed weight change, controlling for calories, carbs, and sodium:

| Variable | Partial r | Signal? |
|---|---|---|
| Calories | -0.004 | Absorbed by TDEE derivation |
| Carbs | -0.002 | Absorbed by glycogen model |
| Sodium | +0.001 | Absorbed by sodium model |
| Sleep hours | +0.094 | Tested as correction — failed (negative variance reduction) |
| Steps | +0.014 | No signal at any timescale |
| Estimated food mass | +0.031 | Confounded with calories |
| Fat | -0.015 | No independent signal |
| Protein | -0.019 | No independent signal |
| Fiber | -0.022 | No independent signal |

Only carbs (glycogen) and sodium have actionable signal for weight smoothing. The remaining ~85% of daily weight variance is unmeasured (gut contents, non-dietary hydration, scale positioning) and requires a state-space model (Kalman filter) rather than additional correction layers.
Command: `python analysis/N_dietary_predictors.py` (all predictors section)

### TDEE by year

| Year | TDEE | Intake | Gap | FM | ~Weight |
|---|---|---|---|---|---|
| 2011 | 2604 | 1772 | +832 | 70 | 215 |
| 2012 | 2200 | 1917 | +283 | 35 | 180 |
| 2013 | 2038 | 1847 | +191 | 20 | 165 |
| 2014 | 2111 | 2105 | +6 | 22 | 165 |
| 2015 | 2100 | 2130 | -30 | 37 | 182 |
| 2016 | 2120 | 1954 | +166 | 36 | 181 |
| 2017 | 2084 | 2140 | -56 | 53 | 198 |
| 2018 | 2056 | 2005 | +51 | 61 | 205 |
| 2019 | 2002 | 2018 | -16 | 62 | 208 |
| 2020 | 2160 | 2120 | +40 | 66 | 211 |
| 2021 | 2154 | 2136 | +18 | 71 | 216 |
| 2022 | 2242 | 2224 | +18 | 72 | 217 |
| 2023 | 2316 | 2267 | +49 | 79 | 224 |
| 2024 | 2378 | 2274 | +103 | 83 | 227 |
| 2025 | 2150 | 1943 | +207 | 71 | 214 |
| 2026 | 2047 | 1844 | +203 | 61 | 204 |

The body defends ~2100-2200 cal/day TDEE across a 60 lb fat mass range. The 2011 outlier (2604) reflects the high metabolic cost of carrying 100 lbs of fat. By 2013, TDEE had dropped to 2038 despite continued restriction — metabolic adaptation in action.
Command: `python analysis/Z_tdee_by_year.py`
Artifact: `analysis/P4_kalman_daily.csv`

## Previous attempt

An earlier breakpoint model (`analysis-bmr-fit/`, now removed) tried to split the timeline into segments with different constant BMR values. Did not work well — the breakpoint approach is too rigid for a continuously drifting set point. The Kalman filter in `analysis/P4_kalman_filter.py` is its successor.
