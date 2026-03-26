# Findings

Each section has a `Command:` line pointing to the standalone script that reproduces the numbers. Untested hypotheses are in [ROADMAP.md](ROADMAP.md).

# Expenditure

## A. Energy balance quality

Does the intake data close against weight and calorimetry?

Cumulative energy balance residual: ±5 lbs over 15 years. Undercount by trajectory phase (composition-aware model, 25 calorimetry × 70 composition measurements): gaining 13.2%, losing 9.7%, stable 13.9%. Approximately uniform.

Command: `python analysis/A_energy_balance_quality.py`

## Z. TDEE by year

Does derived TDEE vary with weight, or does the body defend a band?

| Year | TDEE | Intake | Gap | Fat (lbs) | Lean (lbs) | TDEE/RMR |
|---|---|---|---|---|---|---|
| 2011 | 2604 | 1772 | +832 | 102→49 | 159→157 | 1.198 |
| 2012 | 2200 | 1917 | +283 | 49→32 | 157→154 | 1.037 |
| 2013 | 2038 | 1847 | +191 | 32→20 | 154→152 | 0.983 |
| 2014 | 2111 | 2105 | +6 | 20→33 | 152→149 | 1.048 |
| 2015 | 2100 | 2130 | -30 | 33→46 | 149→146 | 1.061 |
| 2016 | 2120 | 1954 | +166 | 46→47 | 146→141 | 1.097 |
| 2017 | 2084 | 2140 | -56 | 47→58 | 141→144 | 1.065 |
| 2018 | 2056 | 2005 | +51 | 58→64 | 144→150 | 1.040 |
| 2019 | 2002 | 2018 | -16 | 64→69 | 150→146 | 0.977 |
| 2020 | 2160 | 2120 | +40 | 70→69 | 146→146 | 1.087 |
| 2021 | 2154 | 2136 | +18 | 69→74 | 146→143 | 1.079 |
| 2022 | 2242 | 2224 | +18 | 74→76 | 143→143 | 1.137 |
| 2023 | 2316 | 2267 | +49 | 77→83 | 143→144 | 1.166 |
| 2024 | 2378 | 2274 | +103 | 83→80 | 144→140 | 1.183 |
| 2025 | 2150 | 1943 | +207 | 80→62 | 140→142 | 1.099 |
| 2026 | 2047 | 1844 | +203 | 62→60 | 142→143 | 1.048 |

Start→end of year, Kalman-filtered. Fat and lean are lbs (fat-free mass includes water, bone, organs). Lost 80 lbs fat and 18 lbs lean (102→20 fat, 159→141 lean) over 2011-2016, regained 63 lbs fat but only 3 lbs lean (20→83 fat, 141→144 lean) over 2014-2024. On tirzepatide: fat 83→60 while lean stable at 140-143.

Command: `python analysis/Z_tdee_by_year.py`
Artifact: `analysis/P4_kalman_daily.csv`

## B. Weekend fasting as expenditure defense microcosm

Do acute caloric deficits produce lasting fat loss, or does the body recover through reduced expenditure?

Seven consecutive Sat-Sun 36-hour fasts (Oct-Nov 2019). Mean deficit per fast: 3,300 cal. Kalman FM change Fri→Mon: -0.76 lbs (80% of expected). Kalman FM change Fri→Fri+7: +0.15 lbs. Post-fast weekday intake: 2,369 cal (78 cal/day below pre-fast 2,447). Zero compensatory overeating. The deficit disappears through reduced expenditure within a week.

Command: `python analysis/B_weekend_fasting.py`

## K. TDEE hysteresis

At the same fat mass, does TDEE differ depending on whether weight was reached from above (falling) or below (rising)?

Matched fat-mass bands (pre-tirzepatide, retrospective Kalman states):
- FM 25-45 lbs: rising 2091 vs falling 2196 (+106)
- FM 45-65 lbs: rising 2073 vs falling 2358 (+285)
- FM 65-85 lbs: rising 2208 vs falling 2490 (+282)

Regression: TDEE = 1943 + 3.56×FM + 169×falling - 7×rising.

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
- Long runs (≥6 days): +0.0038 (recover best)
- Low-carb (<170g/day): +0.0016 (recover well)
- Low-protein (<58g/day): -0.0119 (worst recovery)
- High-steps (≥4200/day): -0.0167 (-39 cal/day, worst overall)

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

Command: `python analysis/F_tirzepatide_pk.py`
Artifact: `drugs/tirzepatide.csv`

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

Travel binges: 2716 cal/day, 32.1% binge rate, +10.4 lbs fat. Keto: 2163 cal, 110.5g protein, 2.1% binges. Potato diet (4 epochs, 69 days): 1930 cal, 48.6g protein, 0% binges, -7.85 lbs fat.

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
