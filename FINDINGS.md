# Findings

Each section has a `Command:` line pointing to the standalone script that reproduces the numbers. Untested hypotheses are in [ROADMAP.md](ROADMAP.md).

# Expenditure

## A. Energy balance quality

Does the intake data close against weight and calorimetry?

Cumulative energy balance residual: ±5 lbs over 15 years. Undercount by trajectory phase (composition-aware model, 25 calorimetry × 70 composition measurements): gaining 14.7%, losing 8.6%, stable 12.4%. Approximately uniform.

Command: `python analysis/A_energy_balance_quality.py`

## Z. TDEE by year

Does derived TDEE vary with weight, or does the body defend a band?

| Year | TDEE | Intake | Gap | Fat (lbs) | Lean (lbs) | TDEE/RMR |
|---|---|---|---|---|---|---|
| 2011 | 2661 | 1772 | +890 | 106→53 | 159→157 | 1.223 |
| 2012 | 2230 | 1917 | +313 | 53→36 | 157→154 | 1.049 |
| 2013 | 2044 | 1847 | +197 | 36→25 | 154→152 | 0.991 |
| 2014 | 2050 | 2105 | -55 | 25→37 | 152→149 | 1.015 |
| 2015 | 2084 | 2130 | -46 | 37→50 | 149→146 | 1.057 |
| 2016 | 2137 | 1954 | +183 | 50→48 | 146→141 | 1.107 |
| 2017 | 2138 | 2140 | -2 | 49→56 | 141→144 | 1.084 |
| 2018 | 2025 | 2005 | +20 | 57→63 | 144→150 | 1.029 |
| 2019 | 2015 | 2018 | -3 | 63→67 | 150→146 | 0.982 |
| 2020 | 2125 | 2120 | +6 | 68→71 | 146→146 | 1.069 |
| 2021 | 2142 | 2136 | +6 | 71→75 | 146→143 | 1.076 |
| 2022 | 2248 | 2224 | +24 | 75→78 | 143→143 | 1.138 |
| 2023 | 2264 | 2267 | -3 | 79→86 | 143→144 | 1.142 |
| 2024 | 2388 | 2274 | +114 | 86→80 | 144→140 | 1.189 |
| 2025 | 2148 | 1943 | +205 | 80→62 | 140→142 | 1.095 |
| 2026 | 2041 | 1844 | +197 | 62→60 | 142→143 | 1.045 |

Start→end of year, Kalman v2. Fat and lean in lbs. Lost 81 lbs fat and 18 lbs lean (106→25 fat, 159→141 lean) over 2011-2016, regained 61 lbs fat but only 3 lbs lean (25→86 fat, 141→144 lean) over 2014-2024. On tirzepatide: fat 86→60 while lean stable at 140-143.

Command: `python analysis/Z_tdee_by_year.py`
Artifact: `analysis/P4_kalman_daily.csv`

## B. Weekend fasting as expenditure defense microcosm

Do acute caloric deficits produce lasting fat loss, or does the body recover through reduced expenditure?

Seven consecutive Sat-Sun 36-hour fasts (Oct-Nov 2019). Mean deficit per fast: 3,300 cal. Kalman FM change Fri→Mon: -0.76 lbs (80% of expected). Kalman FM change Fri→Fri+7: +0.08 lbs. Post-fast weekday intake: 2,369 cal (78 cal/day below pre-fast 2,447). Zero compensatory overeating. The deficit disappears through reduced expenditure within a week.

Command: `python analysis/B_weekend_fasting.py`

## K. TDEE hysteresis

At the same fat mass, does TDEE differ depending on whether weight was reached from above (falling) or below (rising)?

Matched fat-mass bands (pre-tirzepatide, retrospective Kalman states):
- FM 25-45 lbs: rising 2095 vs falling 2207 (+112)
- FM 45-65 lbs: rising 2082 vs falling 2347 (+266)
- FM 65-85 lbs: rising 2166 vs falling 2498 (+332)

Regression: TDEE = 1892 + 4.28×FM + 202×falling - 19×rising.

Command: `python analysis/K_tdee_hysteresis.py`
Artifact: `analysis/K_tdee_hysteresis_phase_summary.csv`, `analysis/K_tdee_hysteresis_band_summary.csv`, `analysis/K_tdee_hysteresis_regression.csv`

## I. Dead zone is phase-dependent

Is there a single intake band (2000-2500 cal) where metabolism adapts to maintain weight?

No universal dead zone. In falling phase, TDEE is flat across all intake levels (1600-3000 cal). In stable phase, higher intake weakly raises TDEE (coef +0.018). In rising, near zero (+0.005). The pattern is branch-specific expenditure compression, not one dead zone.

Command: `python analysis/I_deadzone_phase.py`
Artifact: `analysis/I_deadzone_phase_bin_summary.csv`, `analysis/I_deadzone_phase_regression.csv`

## J. Restriction archetypes

Do all restriction runs (<1800 cal, ≥3 days) produce the same metabolic penalty?

Average post-run TDEE/RMR penalty: -0.0063 (-19 cal/day). By archetype:
- Long runs (≥6 days): +0.0005 (recover best)
- Low-carb (<170g/day): -0.0005 (recover well)
- Low-protein (<58g/day): -0.0080 (worst recovery)
- High-steps (≥4200/day): -0.0112 (-29 cal/day, worst overall)

Command: `python analysis/J_restriction_archetypes.py`
Artifact: `analysis/J_restriction_runs.csv`, `analysis/J_restriction_archetype_summary.csv`, `analysis/J_restriction_archetype_regression.csv`

## L. Restriction penalties are not rebound eating

Do the worst archetype penalties line up with the biggest post-restriction rebound?

Low-protein runs have larger penalty despite lower next-7-day calories and fewer binges. High-step cuts are worst with only middling rebound. In regression, rebound terms are near zero while long_run stays positive and high_steps stays negative.

Command: `python analysis/L_restriction_rebound.py`
Artifact: `analysis/L_restriction_rebound_summary.csv`, `analysis/L_restriction_rebound_regression.csv`

## M. Bad archetypes fail in different ways

Do high-step and low-protein cuts fail through the same mechanism?

High-step cuts shift TDEE/RMR downward during the run (run-pre -0.0065) and keep falling afterward (post-run -0.0112). Low-protein cuts show near-zero during-run shift (-0.0015) but poor recovery (post-run -0.0104). Long and low-carb runs start on a higher branch and lose little afterward.

Command: `python analysis/M_restriction_branch_shift.py`
Artifact: `analysis/M_restriction_branch_shift_summary.csv`, `analysis/M_restriction_branch_shift_regression.csv`

## R. Metabolic failure predictors

What predicts TDEE/RMR collapse after restriction — branch position or willpower depletion?

Best single predictor of post-run TDEE/RMR failure: phase_code (R²=0.197). Adding mean_steps gives R²=0.239. Falling-phase high-step cuts: post-pre -0.0598. Rising-phase cuts: +0.0305 (low-step), +0.0122 (high-step). Branch state and activity load matter more than control depletion.

Command: `python analysis/R_metabolic_failure_predictors.py`
Artifact: `analysis/R_metabolic_failure_runs.csv`, `analysis/R_metabolic_failure_feature_search.csv`, `analysis/R_metabolic_failure_phase_summary.csv`

## S. Branch effect survives matching on fat mass

Is the falling-phase penalty just a disguised fat-mass effect?

Within matched start-fat bands: FM 50-65, falling post-pre -0.0480 vs rising +0.0340. Regression controlling for start fat mass: falling_vs_rising = -0.0534 for net post-run failure.

Command: `python analysis/S_metabolic_failure_matched.py`
Artifact: `analysis/S_metabolic_failure_matched_bands.csv`, `analysis/S_metabolic_failure_matched_regression.csv`

## T. Branch effect survives matching on metabolic state

Is the falling-phase penalty explained by starting the run already adapted?

Matching on both fat_mass_start and pre-run TDEE/RMR: 43 pairs, falling minus rising = -0.0442 for net post-run failure. Tighter matches (distance ≤0.50, n=34): -0.0462. Regression controlling for both: falling_vs_rising = -0.0398.

Command: `python analysis/T_metabolic_failure_state_matched.py`
Artifact: `analysis/T_metabolic_failure_state_match_summary.csv`, `analysis/T_metabolic_failure_state_match_pairs.csv`, `analysis/T_metabolic_failure_state_match_regression.csv`

# Body composition

## AA. Lean mass response to strength training

Does strength training add lean mass that decays without continued training?

Model: each workout adds Δ lbs of lean mass. Accumulated effect decays exponentially. Lean mass = baseline(weight) + training_effect. Fitted to 70 composition measurements × 393 workout sessions.

The objective surface is flat: 14 parameter combinations within 1% of the best RMSE. Δ ranges from 0.06-0.12 lbs (27-54g) per workout, half-life 240-365 days. Point estimate: 41g, 275 days. RMSE improves from 4.01 (weight-only baseline) to 3.80 (F≈3.9, p≈0.025). Partial correlation (trailing 60-day workouts vs lean mass, controlling weight): r=0.18.

At the point estimate, steady-state lean mass above baseline: 1×/week +5.2 lbs, 2×/week +10.3 lbs. The 8-12 month half-life means training effects persist nearly a year after stopping, consistent with myonuclear persistence. The effect is real but the exact parameters are weakly identified — more composition measurements during detraining/retraining transitions would tighten the fit.

Command: `python analysis/AA_lean_mass_training.py`
Artifact: `analysis/AA_lean_mass_training.csv`

# Intake

## Y. Five negative intake-side set point tests

Does the set point defend through intake, or only through expenditure?

1. Recent weight change vs future intake: r=+0.26. Weight gain predicts more eating, not compensatory restriction.
2. Weekly intake anti-compensates: ratio 1.64. Binges cluster rather than correct.
3. Intake autocorrelation dies at 30 days: r=0.40→0.22→0.17→0.04 over 1→7→14→30 day lags.
4. Binge rate peaks at FM 50 (17.1%), not at minimum weight (FM 20: 10.4%).
5. Binge rate drops 3.6× on tirzepatide (12.2%→3.4%), confirming binges are pharmacologically mediated.

Command: `python analysis/Y_set_point_intake_tests.py`

## C. Binge prediction from set point distance

Does cumulative distance below set point predict binge probability (>2800 cal/day)?

572 binges across 4,807 pre-tirzepatide days (walk-forward, causal predictors only). Distance from set point (180-day trailing average): AUC=0.55. Yesterday's calories: AUC=0.73. Yesterday's protein %: AUC=0.59. Distance adds weak signal but yesterday's intake dominates.

At the same distance from set point (-2 to +2 lbs), tirzepatide reduces binge rate from 11.7% to 3.4%. Logistic regression: distance β=+0.20, tirzepatide β=-0.44.

Command: `python analysis/C_binge_analysis.py`

## D. Food noise as intake variance

The [food noise essay](https://lucent.substack.com/p/craving-food-noise) distinguishes hunger (daily deficit) from food noise (persistent awareness proportional to distance below defended weight). Four tests:

1. Distance → continuous intake: -30 cal/day per kg below set point. Direction correct, 3× weaker than claimed 100 cal/kg.
2. Restriction duration → rebound: r=-0.065. No compounding effect.
3. Intake variance below set point: std=615, CV=0.278 vs at set point std=528, CV=0.244. Tirzepatide reduces CV from 0.24-0.28 to 0.19-0.20 (25-30% reduction).
4. Post-restriction binge rate: 26% pre-drug → 13% on drug.

Food noise manifests as intake variance (erratic eating, clustering binges), not as a constant upward pressure proportional to set point distance.

Command: `python analysis/D_food_noise_variance.py`

## E. Week-scale intake invariance

[Preregistered claim](https://bsky.app/profile/lucent.substack.com): weekly calorie variance is lower than independent daily draws would predict, implying a homeostatic weekly regulator.

Weekly total std: 2383 cal. Expected if independent: 1456 cal. Ratio: 1.64. Rejected — variance is 64% higher, not lower. Binges cluster rather than compensate. Autocorrelation: r=0.40 at 1 day, 0.04 at 30 days.

Command: `python analysis/E_weekly_invariance.py`

## N. Dietary predictors

Do protein leverage, meal timing, fiber, or the [gravitostat](https://doi.org/10.1073/pnas.1800033115) predict intake or weight change?

Protein leverage (same-day): r=-0.34 (protein % vs total intake). 590 cal range across bins. Next-day partial (controlling today's calories): r=-0.04. Within-day only, no carryover. Tirzepatide weakens the effect (r -0.19 → -0.09).

Front-loading: morning % vs daily total r=-0.19, but morning absolute calories vs total r=+0.48. The percentage is circular. Absolute front-loading does not reduce intake.

Fiber: morning fiber controlling for morning calories: r=-0.094 partial. Weak.

Gravitostat: foot-pounds → next-day intake r=0.050 partial. Positive (wrong direction). Steps → weight change: r=0.010. No signal at any timescale.

Command: `python analysis/N_dietary_predictors.py`

## Q. Latent control capacity

Can a hidden "willpower battery" predict next-day binges?

Latent control-debt AUC: 0.696. Yesterday's calories alone: 0.704. Combined: 0.717. Best half-life: 2 days. The latent model adds marginal information. Short-term control depletion is more defensible than long-burn ego depletion.

Command: `python analysis/Q_latent_control_capacity.py`
Artifact: `analysis/Q_latent_control_model_summary.csv`, `analysis/Q_latent_control_feature_search.csv`, `analysis/Q_latent_control_daily.csv`

# Tirzepatide

## F. Pharmacokinetics

One-compartment SC model (FDA: t½=5d, Tmax=24h, ka=3.31/day). 80 weekly injections, dose escalation 2.5→12.5mg.

Blood level → daily intake (partial, controlling time): r=-0.50. Weekly sawtooth: injection day 1652 cal, trough (day 5) 2220 cal. 568 cal/day swing. Intake model: `cal = 2345 - 49 × effective_level`.

Tachyphylaxis: effectiveness decays with 32-week half-life. After 20 weeks on the same dose, 65% effective.

Overall reduction: 456 cal/day (18.6%) from pre-tirzepatide year. Weight loss is 93% fat: FM 83→60 lbs, lean 141→143 lbs.

Direct calorimetry (25 measurements, 3 during tirzepatide era): RMR dropped from 1930 to 1750 over one year. Composition-aware model predicts 1956; actual 1750. 206 cal/day metabolic adaptation beyond tissue loss. The body IS conserving on the drug.

Energy budget: -450 intake, +200 metabolic clawback, net -250 cal/day deficit.

Command: `python analysis/F_tirzepatide_pk.py`, `python analysis/P2_rmr_model.py` (RMR/calorimetry numbers)
Artifact: `drugs/tirzepatide.csv`, `RMR/rmr.csv`

## G. Transition dynamics

Does the drug reduce mean intake or break escalation cascades?

State machine: restriction (<1800), typical (1800-2399), high (2400-2799), binge (2800+).
- High → binge: 19.6% pre-drug → 4.8% on drug
- Binge → binge: 31.6% → 0%
- Post-restriction binge (7 days): 25.9% → 13.0%

The drug breaks escalation and persistence, not just mean level.

Command: `python analysis/G_tirzepatide_dynamics.py`
Artifact: `analysis/G_tirzepatide_transition_summary.csv`, `analysis/G_tirzepatide_rebound_summary.csv`

## H. Phase × binge landscape

Does trajectory phase (falling/stable/rising) add predictive power for binges?

Pre-tirzepatide binge rate: falling 8.9%, stable 12.2%, rising 15.1%. But in walk-forward prediction, adding phase to yesterday's calories slightly hurts (AUC 0.721 vs 0.724). Phase organizes the background landscape but does not improve day-ahead prediction.

Command: `python analysis/H_phase_binge_landscape.py`
Artifact: `analysis/H_phase_binge_summary.csv`, `analysis/H_phase_binge_distance_bins.csv`, `analysis/H_phase_binge_model_auc.csv`

## X. Temperature

Does body temperature track metabolic adaptation? Does tirzepatide override the thermostat?

Pre-tirzepatide: 14-day trailing intake → temperature r=+0.22. The 2600-3200 cal bin is 0.10°F warmer than 2000-2300. Temperature weakly tracks Kalman TDEE (r=0.04).

An earlier analysis found temperature correlating with tirzepatide blood level at r=+0.45. This was confounded by dose escalation timing (low doses in fall/winter, high doses in spring/summer). Direct calorimetry disproves the override hypothesis: RMR 1930→1750, temperature 97.7→97.5°F. Both dropped. The drug does not override the thermostat.

Command: `python analysis/X_temperature_phase.py`, `python analysis/P2_rmr_model.py`
Artifact: `analysis/X_temperature_daily.csv`, `analysis/X_temperature_phase_overlap.csv`, `analysis/X_temperature_phase_band_summary.csv`, `analysis/X_temperature_phase_regression.csv`

# Diet experiments

## O. Epoch analysis and potato diets

Do the 22 manually annotated diet epochs in `intake/diet_epochs.csv` separate meaningful regimes? Do potato diets have special metabolic properties?

Travel binges: 2716 cal/day, 32.1% binge rate, +9.5 lbs fat. Keto: 2163 cal, 110.5g protein, 2.1% binges. Potato diet (4 epochs, 69 days): 1930 cal, 48.6g protein, 0% binges, -7.3 lbs fat.

Potato diets: TDEE/RMR 1.168 before, 1.159 during (unchanged), 1.098 after. Post-potato rebound: 2629 cal/day, 19.6% binge rate. Potatoes are a monotone binge-suppressing cut, not a uniquely expenditure-preserving intervention.

Command: `python analysis/O_diet_epoch_analysis.py`
Artifact: `analysis/O_diet_epoch_summary.csv`, `analysis/O_diet_epoch_family_summary.csv`, `analysis/O_potato_epoch_window_summary.csv`, `analysis/O_potato_epoch_contrast.csv`

## P. Jordan trips

Do month-long travel periods produce compensatory TDEE increases?

2015: intake 3262 cal/day, TDEE 2071 (-35 vs pre-trip), +10.3 lbs fat. 2019: intake 2529, TDEE 2021 (-21 vs pre-trip), +4.4 lbs fat. No compensatory TDEE surge. The trips are binge windows without an energy-out bonus.

Command: `python analysis/P_jordan_trip_analysis.py`
Artifact: `analysis/P_jordan_trip_summary.csv`, `analysis/P_jordan_trip_delta.csv`

# Steps and activity

## U. Steps compensation

Does higher step load reduce future intake or raise TDEE?

7-day average steps → future 7-day intake: β=+0.103 (more steps → more eating). In falling phase, top 20% of 7-day steps (≥6311/day) is followed by lower 14-day intake (1974 vs 2088) and lower TDEE/RMR (1.0918 vs 1.1569). Closer to constrained-energy compensation than gravitostat benefit.

Command: `python analysis/U_steps_compensation.py`
Artifact: `analysis/U_steps_compensation_regression.csv`, `analysis/U_steps_compensation_phase_thresholds.csv`

## V. Deliberate walk sessions

Do paired daylight walks differ from generic high-step days?

1,636 walking sessions (type 1001). Post-2020 noon-7pm year-round: 149 paired-daylight-walk days, 153 single-daylight-walk days. Step-matched paired walks vs other days: -57 cal lower future 14-day intake, +0.019 higher future TDEE/RMR. Single walks: +17 cal (no intake benefit), +0.027 TDEE/RMR. The paired-walk structure matters more than raw step count.

Command: `python analysis/V_exercise_walk_analysis.py`
Artifact: `analysis/V_exercise_type_summary.csv`, `analysis/V_daylight_walk_regime_days.csv`, `analysis/V_daylight_walk_regime_summary.csv`, `analysis/V_daylight_walk_regime_contrast.csv`, `analysis/V_daylight_walk_pair_days.csv`, `analysis/V_daylight_walk_pair_summary.csv`, `analysis/V_daylight_walk_pair_contrast.csv`

## W. Walk compactness

Does concentrated walking (compact sessions) differ from smeared-out steps?

Single-daylight-walk days vs step-matched smeared days separated by compactness metrics (steps-per-span-minute, top-60-minute share). Details in artifact.

Command: `python analysis/W_daylight_walk_compactness.py`
Artifact: `analysis/W_daylight_walk_compactness_daily.csv`, `analysis/W_daylight_walk_compactness_matches.csv`, `analysis/W_daylight_walk_compactness_summary.csv`

## AB. Running vs walking

Does running suppress appetite more than walking at the same step count? (Exercise-induced anorexia hypothesis.)

121 running sessions (2013-2023), 1,640 walking sessions. All-eras step-matched comparison showed -338 cal same-day for running — but this was entirely era confounding (running clustered in the 2014-2016 restriction era).

Era-matched within 2014-2016 only, step-matched (90 pairs, mean 9650 steps): same-day cal difference -50, next-day +30, next-3d -46. All within noise. TDEE/RMR: run 1.072 vs walk 1.098 (-0.026). No exercise-induced anorexia detected. Running and walking produce equivalent intake effects at matched steps within the same era.

Command: `python analysis/AB_running_vs_walking.py`

## AC. The 2013 inflection

Fat mass bottomed at 17 lbs (Oct 2013) and rose every year for a decade. What triggered the regain — metabolic adaptation forcing intake up, or behavioral momentum?

Bottom (2013 Jul-Sep): intake 1864, TDEE 2021, ratio 0.977, binge rate 1%. Regain (2014 Jul-Dec): intake 2359, TDEE 2058, ratio 1.026, binge rate 20%.

TDEE recovered +38 cal from bottom to regain. Intake increased +495 cal. Intake outran TDEE recovery by 457 cal/day. The metabolic adaptation was easing — ratio went from 0.98 to 1.03. But binge frequency exploded: 0% (Jul-Aug 2013) → 13% (Oct) → 20-23% (Nov-Dec) → sustained 12-26% through all of 2014. By 2014-H2, binges were running at 20%.

The trigger was behavioral, not metabolic. The body's expenditure defense was releasing. Food noise — the variance, the binge clustering — took hold at FM=17 and never let go for a decade.

Command: `python analysis/AC_inflection_2013.py`

## AD. What predicts measured RMR beyond body composition?

Exhaustive sweep of trailing dietary, activity, and sleep features at 7 windows (7-60 days) against 23 Cosmed Fitmate calorimetry measurements with real Samsung Health data (2016+, dropping 2 pre-Samsung measurements with backfilled steps). Leave-one-cluster-out cross-validation with ridge regression. 12 independent clusters.

**Null results.** Baseline (expected_rmr only): CV RMSE = 183. Dietary features (calories, protein %, carbs, fat, sodium) at every window: best CV RMSE = 169, no improvement over the Fitmate noise floor (~170 cal). Step counts at every window: best CV RMSE = 151 (steps_14d). Sleep at every window: negative R² throughout. Strength training count: no signal.

**Walk sessions predict RMR.** The count of deliberate walks Samsung Health logged as exercise in the prior 30 days predicts RMR at CV RMSE = 116 (R² = 0.49), well below the Fitmate noise floor. The coefficient: +14 cal RMR per walk session. Going from 3 walks/month (Sep 2025, RMR 1750) to 33 walks/month (May 2022, RMR 2292) corresponds to a 420 cal/day difference — at nearly identical body composition (FM 66-68 in 2022 vs 68 in 2025).

Walk sessions beat walk minutes (CV RMSE 135) and total steps (CV RMSE 179) at the same 30-day window. The count of distinct outings matters more than total duration or total movement.

**Not season.** Walk sessions and summer are correlated at r = 0.83 (more walks in warm months). But: (1) Season alone barely helps — `is_summer` reduces CV RMSE from 183 to only 175, while walk sessions reduce it to 116. (2) Adding season on top of walk sessions adds nothing (CV RMSE 116 → 116). (3) Controlling for season, the partial correlation between walk sessions and RMR is still r = 0.47 (steps drop to r = 0.13). (4) Within winter only (n=12), walk sessions still predict RMR at r = 0.63. (5) Within 2022-2023 where FM was nearly constant (66-75 lbs, n=18), walk sessions vs RMR: r = 0.68, partial controlling for expected_rmr: r = 0.62. (6) Within the May-Jun 2022 cluster alone (n=8, same FM=66-68, same season), walk sessions vs RMR: r = 0.69.

**Not tirzepatide, not body composition.** Only 2 of 23 measurements are on the drug. Adding tirz_level to the model doesn't change the walk session coefficient (CV RMSE 117 vs 116). Walk sessions are uncorrelated with fat mass (r = -0.02) and negatively correlated with expected_rmr (r = -0.28) — the effect is not mediated by composition changes.

**Interpretation.** Deliberate sustained walks (typically 20+ min continuous, enough for Samsung Health to log as an exercise session) raise resting metabolic rate in a way that total step count — which includes all incidental movement — does not. The mechanism is consistent with NEAT upregulation: intentional exercise sessions may activate a metabolic afterburn that persists at rest, while shuffling around the house does not. The effect size (+14 cal/session, ~420 cal/day at 30 vs 3 sessions/month) is large but consistent across subgroups. The remaining confound is that with 23 measurements clustered by month, we cannot fully exclude an unmeasured seasonal factor that drives both walking and RMR.

Command: `python analysis/AD_tdee_formula_sweep.py`
