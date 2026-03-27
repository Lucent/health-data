# Findings

Each section has a `Command:` line pointing to the standalone script that reproduces the numbers. Untested hypotheses are in [ROADMAP.md](ROADMAP.md).

# Data quality

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

# Expenditure defense

The year-by-year table shows TDEE compressing from 2661 to 2044 as fat mass fell, then barely recovering during a decade of regain. This raises the question: is the compression symmetric, or does direction matter?

## K. TDEE hysteresis

At the same fat mass, does TDEE differ depending on whether weight was reached from above (falling) or below (rising)?

Matched fat-mass bands (pre-tirzepatide, retrospective Kalman states):
- FM 25-45 lbs: rising 2095 vs falling 2207 (+112)
- FM 45-65 lbs: rising 2082 vs falling 2347 (+266)
- FM 65-85 lbs: rising 2166 vs falling 2498 (+332)

Regression: TDEE = 1892 + 4.28×FM + 202×falling - 19×rising.

Robustness: 1,482 FM-matched pairs (within 2 lbs), mean difference +179 cal (falling - rising). Pair lag-1 autocorrelation 0.83, effective n = 142. Bootstrap 95% CI: [+104, +267]. p < 10^-22. This is the most statistically robust finding in the dataset.

**No universal dead zone.** A corollary: there is no single intake band (2000-2500 cal) where metabolism adapts to maintain weight. In falling phase, TDEE is flat across all intake levels (1600-3000 cal). In stable phase, higher intake weakly raises TDEE (coef +0.018). In rising, near zero (+0.005). The pattern is branch-specific expenditure compression, not one dead zone.

**Survives matching on fat mass.** Within matched start-fat bands: FM 50-65, falling post-pre -0.0480 vs rising +0.0340. Regression controlling for start fat mass: falling_vs_rising = -0.0534 for net post-run failure.

**Survives matching on metabolic state.** Matching on both fat_mass_start and pre-run TDEE/RMR: 43 pairs, falling minus rising = -0.0442 for net post-run failure. Tighter matches (distance ≤0.50, n=34): -0.0462. Regression controlling for both: falling_vs_rising = -0.0398.

Command: `python analysis/K_tdee_hysteresis.py`, `python analysis/I_deadzone_phase.py`, `python analysis/S_metabolic_failure_matched.py`, `python analysis/T_metabolic_failure_state_matched.py`
Artifact: `analysis/K_tdee_hysteresis_phase_summary.csv`, `analysis/K_tdee_hysteresis_band_summary.csv`, `analysis/K_tdee_hysteresis_regression.csv`, `analysis/I_deadzone_phase_bin_summary.csv`, `analysis/I_deadzone_phase_regression.csv`, `analysis/S_metabolic_failure_matched_bands.csv`, `analysis/S_metabolic_failure_matched_regression.csv`, `analysis/T_metabolic_failure_state_match_summary.csv`, `analysis/T_metabolic_failure_state_match_pairs.csv`, `analysis/T_metabolic_failure_state_match_regression.csv`

## B. Weekend fasting as expenditure defense microcosm

If TDEE is direction-dependent, acute deficits should vanish through expenditure compression rather than compensatory eating. A natural experiment tests this.

Do acute caloric deficits produce lasting fat loss, or does the body recover through reduced expenditure?

Seven consecutive Sat-Sun 36-hour fasts (Oct-Nov 2019). Mean deficit per fast: 3,300 cal. Kalman FM change Fri→Mon: -0.76 lbs (80% of expected). Kalman FM change Fri→Fri+7: +0.08 lbs. Post-fast weekday intake: 2,369 cal (78 cal/day below pre-fast 2,447). Zero compensatory overeating. The deficit disappears through reduced expenditure within a week.

Command: `python analysis/B_weekend_fasting.py`

## J. Restriction archetypes

The hysteresis in K operates at the macro scale. Within individual restriction runs, does the type of cut matter?

Do all restriction runs (<1800 cal, ≥3 days) produce the same metabolic penalty?

Average post-run TDEE/RMR penalty: -0.0063 (-19 cal/day). By archetype:
- Long runs (≥6 days): +0.0005 (recover best)
- Low-carb (<170g/day): -0.0005 (recover well)
- Low-protein (<58g/day): -0.0080 (worst recovery)
- High-steps (≥4200/day): -0.0112 (-29 cal/day, worst overall)

**Not explained by rebound eating.** Low-protein runs have larger penalty despite lower next-7-day calories and fewer binges. High-step cuts are worst with only middling rebound. In regression, rebound terms are near zero while long_run stays positive and high_steps stays negative.

**Different failure mechanisms.** High-step cuts shift TDEE/RMR downward during the run (run-pre -0.0065) and keep falling afterward (post-run -0.0112). Low-protein cuts show near-zero during-run shift (-0.0015) but poor recovery (post-run -0.0104). Long and low-carb runs start on a higher branch and lose little afterward.

**Phase dominates.** Best single predictor of post-run TDEE/RMR failure: phase_code (R²=0.197). Adding mean_steps gives R²=0.239. Falling-phase high-step cuts: post-pre -0.0598. Rising-phase cuts: +0.0305 (low-step), +0.0122 (high-step). Branch state and activity load matter more than control depletion.

Command: `python analysis/J_restriction_archetypes.py`, `python analysis/L_restriction_rebound.py`, `python analysis/M_restriction_branch_shift.py`, `python analysis/R_metabolic_failure_predictors.py`
Artifact: `analysis/J_restriction_runs.csv`, `analysis/J_restriction_archetype_summary.csv`, `analysis/J_restriction_archetype_regression.csv`, `analysis/L_restriction_rebound_summary.csv`, `analysis/L_restriction_rebound_regression.csv`, `analysis/M_restriction_branch_shift_summary.csv`, `analysis/M_restriction_branch_shift_regression.csv`, `analysis/R_metabolic_failure_runs.csv`, `analysis/R_metabolic_failure_feature_search.csv`, `analysis/R_metabolic_failure_phase_summary.csv`

# Set point

The expenditure defense in K and B raises a deeper question: is there a single mechanism — a defended fat mass — driving both the expenditure compression and the appetite pressure that ultimately overwhelms it?

## AG. Binge frequency reveals a hidden, moving fat mass set point

Binge rate (days > TDEE + 1000 cal) can be used to reverse-engineer the body's defended weight. The set point that best predicts binge frequency is an exponential moving average of **fat mass** with a **50-day half-life** (~7 weeks). Distance below this set point predicts 90-day binge rate at r = -0.62, and at partial r = -0.66 after controlling for absolute fat mass.

**It's a lipostat, not a gravitostat.** Fat mass predicts binge rate better than total weight (r = -0.62 vs -0.54) or scale weight (r = -0.53) at every half-life tested. The body tracks its own fat stores, not total body mass.

**It moves, not fixed.** A moving EMA (r = -0.62) crushes any fixed set point (best fixed: r = +0.25 at FM = 30 lbs, barely above zero). The defended weight adapts to sustained changes. After ~150 days (3 half-lives, 5 months) at a new weight, the set point has 87% adapted.

**Half-life: 50 days.** The plateau from 35-65 days is broad (r = -0.611 to -0.618), so the half-life is not precisely identified, but the timescale is clear: weeks, not months or years. This explains why the first few months of regain have the most binges (2014-2015: 10-12%) while by 2018-2021 (FM stable at 60-71 lbs) binges had calmed to 3-6%.

**Sigmoid response curve.** Binge rate by distance below the reconstructed set point:

| Distance from SP | Binge rate |
|---|---|
| 5-7.5 lbs below | 14.8% |
| 2.5-5 lbs below | 11.4% |
| 0-2.5 lbs below | 6.0% |
| 0-2.5 lbs above | 3.3% |
| 2.5-5 lbs above | 2.8% |
| 5-10 lbs above | 1.3% |
| 10+ lbs above | 0.8% |

Baseline ~3% at/above set point, rising to ~15% at 5+ lbs below. The gradient is continuous — every pound below the set point adds roughly 0.5% to the binge rate. No hard threshold.

**Reconstructed set point trajectory:**

| Year | FM | Set Point | Distance | Binge % | State |
|---|---|---|---|---|---|
| 2011 | 73 | 79 | +10 | 0.8% | Losing (SP chasing down) |
| 2012 | 36 | 42 | +5 | 4.6% | Losing |
| 2013 | 21 | 23 | +3 | 4.4% | Near bottom |
| 2014 | 24 | 22 | -2 | 9.6% | FM overshot SP → binges |
| 2015 | 38 | 35 | -3 | 11.8% | Peak binges, gaining |
| 2016 | 37 | 38 | +1 | 10.1% | SP catching up |
| 2017 | 53 | 50 | -3 | 8.8% | Still gaining |
| 2018 | 61 | 60 | -1 | 3.8% | Converging |
| 2019-21 | 63-71 | 62-70 | -1 | 3-7% | At set point |
| 2022-24 | 71-82 | 71-82 | 0 | 2-7% | At set point |
| 2025 | 70 | 74 | +3 | 0.0% | Tirzepatide: below SP, zero binges |

The 2025 data point is the cleanest test: FM = 70 (same as 2020) but 0% binges because 70 is *above* the trailing set point of 74 (tirzepatide drove FM down from 84). In 2020, FM 66 was *at* the set point and binge rate was 3.8%. Same absolute weight, opposite relationship to the set point, completely different binge behavior.

**Intake anti-compensation.** Earlier tests (Y, C, E) asked whether the set point defends through intake and found the answer was "yes, but not through homeostatic correction." Weekly intake variance is 64% *higher* than independent daily draws would predict (weekly std 2383 cal vs expected 1456 if independent; ratio 1.64; [preregistered](https://bsky.app/profile/lucent.substack.com) as the opposite) — binges cluster rather than compensate. Intake autocorrelation dies at 30 days (r=0.40→0.22→0.17→0.04 over 1→7→14→30 day lags). Weight gain predicts more eating (r=+0.26), not compensatory restriction. Binge rate peaks at FM 50 (17.1%), not at minimum weight (FM 20: 10.4%) — explained by the moving set point, since FM 50 was the peak distance below. Tirzepatide reduces overall binge rate 3.6× (12.2%→3.4%), confirming binges are pharmacologically mediated.

In walk-forward binge prediction (572 binges across 4,807 pre-tirzepatide days, causal predictors only): distance from set point AUC=0.55, yesterday's calories AUC=0.73, yesterday's protein % AUC=0.59. Distance adds weak signal but yesterday's intake dominates. Pre-tirzepatide binge rates by trajectory phase — falling 8.9%, stable 12.2%, rising 15.1% — organize the background landscape but do not improve day-ahead prediction (AUC 0.721 vs 0.724 for yesterday's calories alone). FM velocity correlates with binge rate at r = +0.62, confirming the mechanism: binges are how the weight comes back.

**Tirzepatide interaction.** The drug suppresses binges independently of the set point: at matched SP distance (0 to +2.5 lbs above SP), off-tirz binge rate is 3.9%, on-tirz is 0.5% (partial r = -0.23 controlling for SP distance). Each unit of effective drug level reduces binge rate by 0.4%, equivalent to being ~2.5 lbs further above the set point. A drug level of 10 offsets ~5 lbs of set point deficit.

The drug slows set point adaptation. Off tirz, the optimal SP half-life is 50 days. On tirz it is 165 days — the set point chases FM down more slowly because the drug is overriding the appetite signal that normally drives convergence. During 17 months on tirz, FM dropped 84 → 62 while the set point dropped only 83 → 63, maintaining a persistent 2-4 lb gap. Without the drug, this gap would produce 4-6% binge rate. With it, 0%.

The injection cycle is visible in daily calories: 1640-1740 cal on days 0-1 post-injection, rising to 2140-2220 cal by days 3-5. A ~500 cal/day sawtooth within each weekly cycle.

The implication for discontinuation: if the drug's appetite suppression is removed, the set point will be wherever it has adapted to (currently ~63 lbs FM). If the 165-day on-drug half-life reflects genuine slower adaptation (not just suppressed appetite masking the signal), the set point may not have fully adapted and binge risk could spike. If the half-life difference is entirely due to suppressed binges (the mechanism the set point normally uses to drive convergence), then the set point has been adapting at the normal 50-day rate and discontinuation risk is lower.

Command: `python analysis/AG_binge_set_point.py`, `python analysis/Y_set_point_intake_tests.py`, `python analysis/C_binge_analysis.py`, `python analysis/E_weekly_invariance.py`, `python analysis/H_phase_binge_landscape.py`
Artifact: `analysis/H_phase_binge_summary.csv`, `analysis/H_phase_binge_distance_bins.csv`, `analysis/H_phase_binge_model_auc.csv`

## AH. Set point mechanics: ratchet, dual defense, and restriction prediction

Six tests characterizing the hidden fat mass set point identified in finding AG.

**1. Inverted ratchet (suggestive, not confirmed).** An asymmetric model (HL_down=20d, HL_up=100d) improves in-sample r from -0.62 to -0.71, suggesting the set point adapts faster going down than up. However, 90-day block bootstrap shows the improvement is not robust: asymmetric wins only 53% of resamples, Δ|r| 95% CI [-0.054, +0.054] includes zero. Leave-one-year-out CV: wins 9/16 years. With ~59 independent 90-day windows, the data cannot reliably distinguish one extra parameter. The direction is interesting — if real, it means the body quickly accepts lower weight but is slow to raise its defended weight during regain, the opposite of the feared "ratchet" — but this needs more data to confirm.

**2. Dual defense: appetite + expenditure (appetite confirmed, expenditure suggestive).** TDEE residual (actual minus composition-predicted) correlates with SP distance at partial r = +0.38 (controlling for FM). However, the Kalman TDEE has residual lag-1 autocorrelation of 0.9995, giving effective n = 4 by Bartlett correction — the TDEE defense arm is not statistically robust at this resolution (p = 0.69 at n_eff=4, p = 0.002 at n_windows=60). Finding K separately confirms TDEE hysteresis (+179 cal falling vs rising, p < 10^-22) through FM-matched pairs, but the direct link between SP distance and TDEE defense is not formally demonstrated. The appetite arm (binge rate, p < 0.04 at n_eff=11) is the only statistically confirmed arm of the set point defense.

**Reconciling AG and K.** AG's 50-day half-life is appetite-specific — it was optimized to predict binge frequency and says nothing about expenditure adaptation timescales. K's hysteresis (+179 cal) measures expenditure during active weight change: a day is classified "falling" if FM dropped 3+ lbs in the prior 90 days. This captures the *transient* defense response while the set point lags behind FM, not a permanent memory at matched fat mass. The two findings describe the same underlying asymmetry from different measurement angles: during a falling phase, FM is below the set point, which elevates TDEE (K) and suppresses binges relative to stable phases at the same FM (AG). Neither contradicts the other — the appetite set point adapts in ~50 days, the expenditure defense fires during the lag period, and both quiet down once the set point catches up.

**3. Frequency not magnitude.** Once a binge fires, its size is constant: ~1421 cal surplus regardless of SP distance (r = -0.04). But non-binge days show continuous upward pressure: r = -0.35 between SP distance and daily surplus. At 5+ lbs below SP, even non-binge days average +113 cal surplus. At 5+ lbs above: -533 cal deficit. The set point doesn't modulate binge intensity — it tilts the background drift on every day, plus the probability of a discrete binge event.

**4. Restriction run prediction.** 237 restriction runs (3+ days, calories < TDEE - 200). Runs ending with FM above the set point continue losing: -1.23 lbs over the next 30 days. Runs ending below the set point rebound: +0.35 lbs over 30 days. SP distance at end of run vs 30-day rebound: r = -0.48. The set point model predicts which calorie cuts will stick.

**5. Exercise independence.** Walking raises RMR (finding AD) but does not change set point dynamics. High-walk and low-walk periods both show optimal SP half-life of 40 days. The two mechanisms are independent: exercise changes metabolic rate, the set point changes appetite. Neither modulates the other.

**6. Floor effect.** The set point reached a minimum of 18.6 lbs FM (Nov 2013), near essential body fat. Adaptation rate at SP 20-30 lbs (1.1 lbs/month absolute) was comparable to higher levels but the *net* rate was near zero — the set point stalled. In the 6 months after the floor: FM rose 19→23, binge rate 8.3%. The floor is where the 2014 inflection began (finding AC).

**Robustness of the core set point model (finding AG).** The 90-day smoothed binge rate has residual lag-1 autocorrelation of 0.996, giving an effective n of only 11 by the Bartlett formula. Even at n=11: r = -0.62, p = 0.036. Using raw daily binge events (0/1) instead of the smoothed rate: r = -0.12, residual ACF = 0.19, effective n = 3639, p < 10^-13. The set point effect is real under any reasonable correction. Block bootstrap (90-day blocks, 5000 resamples) 95% CI for r_sym: excludes zero. The core finding survives; the asymmetric refinement does not yet.

Command: `python analysis/AH_set_point_properties.py`

## AI. Expenditure arm timescale

AG shows the appetite arm adapts with a 50-day half-life. Does the expenditure arm have a measurable adaptation timescale, and is it different?

Sweep EMA half-lives for fat mass, correlate (EMA - FM) with TDEE residual (Kalman TDEE minus composition-predicted RMR), controlling for absolute FM.

**The expenditure arm is faster, not slower.** Optimal half-life: 10 days (partial r = +0.52). The correlation is monotonically stronger at shorter half-lives down to ~9 days, then plateaus. For comparison, appetite peaks at 50 days (r = -0.60). At HL=50: appetite r = -0.60, expenditure partial r = +0.41.

| HL (days) | r (appetite) | r (expenditure, partial) |
|---|---|---|
| 10 | — | +0.52 |
| 30 | -0.59 | +0.47 |
| 50 | -0.60 | +0.41 |
| 100 | -0.56 | +0.33 |
| 200 | -0.48 | +0.29 |

TDEE residual by SP distance bin (HL=10d): 2.5-5 lbs below SP, TDEE-RMR = +49 cal. 0-2.5 above: +198 cal. 2.5-5 above: +502 cal. The gradient is steep — every lb above the set point adds ~100 cal/day to the TDEE residual.

**Not tautological.** Short-HL correlations risk measuring Kalman filter mechanics (TDEE inferred from recent weight change). Lagging SP distance by 7-60 days: partial r drops from 0.52 (lag 0) to 0.47 (lag 7) to 0.30 (lag 30) to 0.19 (lag 60). The signal decays but does not collapse — a purely mechanical artifact would vanish by lag 7-14. SP distance 30 days ago still predicts today's TDEE residual, which is physiological.

**Statistical robustness.** Lag-1 ACF of residual product: 0.993. Bartlett effective n = 17. At n_eff=17: p = 0.03. At n_blocks(90d) = 54: p < 0.001.

**Interpretation.** The expenditure arm fires fast (HL ≤10 days) and dissipates quickly. The appetite arm (HL ~50 days) is slower to engage but persists longer. This explains the 2013 inflection (AC): expenditure defense had already adapted (TDEE/RMR recovered from 0.98 to 1.03) while appetite defense was still active (binge rate 0%→20%). The two arms operate on different timescales — expenditure is a sprint, appetite is a marathon.

Command: `python analysis/AI_expenditure_arm_timescale.py`
Artifact: `analysis/AI_expenditure_arm_sweep.csv`

# Intake characterization

The set point in AG drives appetite through binge frequency. But what shapes day-to-day intake beyond the set point? Several dietary hypotheses predict intake modulation that should be visible in this dataset.

## D. Food noise as intake variance

The [food noise essay](https://lucent.substack.com/p/craving-food-noise) distinguishes hunger (daily deficit) from food noise (persistent awareness proportional to distance below defended weight). Food noise manifests as intake variance (erratic eating, clustering binges), not as a constant upward pressure proportional to set point distance:

- Restriction duration → rebound: r=-0.065. No compounding effect.
- Intake variance below set point: std=615, CV=0.278 vs at set point std=528, CV=0.244. Tirzepatide reduces CV from 0.24-0.28 to 0.19-0.20 (25-30% reduction).

The distance-to-intake gradient (-30 cal/day per kg below set point) is 3× weaker than the essay's claimed 100 cal/kg, consistent with the set point acting through discrete binge events (AG) rather than continuous pressure.

Command: `python analysis/D_food_noise_variance.py`

## N. Dietary predictors

Do protein leverage, meal timing, fiber, or the [gravitostat](https://doi.org/10.1073/pnas.1800033115) predict intake or weight change?

Protein leverage (same-day): r=-0.34 (protein % vs total intake). 590 cal range across bins. Next-day partial (controlling today's calories): r=-0.04. Within-day only, no carryover. Tirzepatide weakens the effect (r -0.19 → -0.09).

Front-loading: morning % vs daily total r=-0.19, but morning absolute calories vs total r=+0.48. The percentage is circular. Absolute front-loading does not reduce intake.

Fiber: morning fiber controlling for morning calories: r=-0.094 partial. Weak.

Gravitostat: foot-pounds → next-day intake r=0.050 partial. Positive (wrong direction). Steps → weight change: r=0.010. No signal at any timescale.

Command: `python analysis/N_dietary_predictors.py`

## AF. Intake variance is mildly protective, not fattening

If food noise manifests as variance (finding D), does that variance cause metabolic damage, or is it neutral?

At the same caloric surplus, does higher day-to-day calorie variance cause more fat gain ("metabolic damage from yo-yo dieting") or less?

Controlling for surplus (mean intake - TDEE), trailing 14-day calorie standard deviation predicts slightly *less* fat gain: partial r = -0.20, coefficient -0.0011 lbs/day per cal of std (bootstrap 95% CI: -0.0013 to -0.0010, significant). At 30 days: partial r = -0.10, also significant.

At matched calorie levels (30-day tertiles): low-calorie + low-variance periods lose 1.3 lbs/month; low-calorie + high-variance periods lose only 0.5 lbs/month. But high-calorie + high-variance periods gain 0.5 lbs/month vs 0.7 for low-variance. The protective effect is strongest at the high end — variable eating while in surplus is less fattening than consistent surplus.

The mechanism is likely finding B: high-variance periods include fasting days that transiently raise TDEE through expenditure defense. The variance itself is not protective — it proxies for intermittent acute deficits. Diet epochs confirm: weekend fasting (CV = 0.68) lost weight despite moderate mean; COVID lockdown (CV = 0.18, very consistent eating) gained despite similar mean.

The effect is small — adding variance to a surplus-only model reduces RMSE from 1.019 to 0.998 (2%). But the sign is the opposite of "yo-yo dieting is metabolically damaging." In this dataset, consistent overeating is more fattening than variable overeating at the same mean.

Command: `python analysis/AF_intake_variance.py`

# Body composition

## AA. Lean mass response to strength training

Does strength training add lean mass that decays without continued training?

Model: each workout adds Δ lbs of lean mass. Accumulated effect decays exponentially. Lean mass = baseline(weight) + training_effect. Fitted to 70 composition measurements × 393 workout sessions.

The objective surface is flat: 14 parameter combinations within 1% of the best RMSE. Δ ranges from 0.06-0.12 lbs (27-54g) per workout, half-life 240-365 days. Point estimate: 41g, 275 days. RMSE improves from 4.01 (weight-only baseline) to 3.80 (F≈3.9, p≈0.025). Partial correlation (trailing 60-day workouts vs lean mass, controlling weight): r=0.18.

At the point estimate, steady-state lean mass above baseline: 1×/week +5.2 lbs, 2×/week +10.3 lbs. The 8-12 month half-life means training effects persist nearly a year after stopping, consistent with myonuclear persistence. The effect is real but the exact parameters are weakly identified — more composition measurements during detraining/retraining transitions would tighten the fit.

Command: `python analysis/AA_lean_mass_training.py`
Artifact: `analysis/AA_lean_mass_training.csv`

# Steps and activity

The set point controls appetite (AG) and the expenditure branch depends on trajectory direction (K). But finding AD shows a third lever — deliberate walking — that modulates RMR independently of both. This section characterizes how activity interacts with the energy balance system.

## AD. Walk sessions predict RMR

Exhaustive sweep of trailing dietary, activity, and sleep features at 7 windows (7-60 days) against 23 Cosmed Fitmate calorimetry measurements with real Samsung Health data (2016+, dropping 2 pre-Samsung measurements with backfilled steps). Leave-one-cluster-out cross-validation with ridge regression. 12 independent clusters.

**Null results.** Baseline (expected_rmr only): CV RMSE = 183. Dietary features (calories, protein %, carbs, fat, sodium) at every window: best CV RMSE = 169, no improvement over the Fitmate noise floor (~170 cal). Step counts at every window: best CV RMSE = 151 (steps_14d). Sleep at every window: negative R² throughout. Strength training count: no signal.

**Walk sessions predict RMR.** The count of deliberate walks Samsung Health logged as exercise in the prior 30 days predicts RMR at CV RMSE = 116 (R² = 0.49), well below the Fitmate noise floor. The coefficient: +14 cal RMR per walk session. Going from 3 walks/month (Sep 2025, RMR 1750) to 33 walks/month (May 2022, RMR 2292) corresponds to a 420 cal/day difference — at nearly identical body composition (FM 66-68 in 2022 vs 68 in 2025).

Walk sessions beat walk minutes (CV RMSE 135) and total steps (CV RMSE 179) at the same 30-day window. The count of distinct outings matters more than total duration or total movement.

**Paired daylight walks.** 1,636 walking sessions (type 1001). Post-2020 noon-7pm year-round: 149 paired-daylight-walk days, 153 single-daylight-walk days. Step-matched paired walks vs other days: -57 cal lower future 14-day intake, +0.019 higher future TDEE/RMR. Single walks: +17 cal (no intake benefit), +0.027 TDEE/RMR. The paired-walk structure matters more than raw step count.

**Not season.** Walk sessions and summer are correlated at r = 0.83 (more walks in warm months). But: (1) Season alone barely helps — `is_summer` reduces CV RMSE from 183 to only 175, while walk sessions reduce it to 116. (2) Adding season on top of walk sessions adds nothing (CV RMSE 116 → 116). (3) Controlling for season, the partial correlation between walk sessions and RMR is still r = 0.47 (steps drop to r = 0.13). (4) Within winter only (n=12), walk sessions still predict RMR at r = 0.63. (5) Within 2022-2023 where FM was nearly constant (66-75 lbs, n=18), walk sessions vs RMR: r = 0.68, partial controlling for expected_rmr: r = 0.62. (6) Within the May-Jun 2022 cluster alone (n=8, same FM=66-68, same season), walk sessions vs RMR: r = 0.69.

**Not tirzepatide, not body composition.** Only 2 of 23 measurements are on the drug. Adding tirz_level to the model doesn't change the walk session coefficient (CV RMSE 117 vs 116). Walk sessions are uncorrelated with fat mass (r = -0.02) and negatively correlated with expected_rmr (r = -0.28) — the effect is not mediated by composition changes.

**Interpretation.** Deliberate sustained walks (typically 20+ min continuous, enough for Samsung Health to log as an exercise session) raise resting metabolic rate in a way that total step count — which includes all incidental movement — does not. The mechanism is consistent with NEAT upregulation: intentional exercise sessions may activate a metabolic afterburn that persists at rest, while shuffling around the house does not. The effect size (+14 cal/session, ~420 cal/day at 30 vs 3 sessions/month) is large but consistent across subgroups. The remaining confound is that with 23 measurements clustered by month, we cannot fully exclude an unmeasured seasonal factor that drives both walking and RMR.

Command: `python analysis/AD_tdee_formula_sweep.py`, `python analysis/V_exercise_walk_analysis.py`
Artifact: `analysis/V_exercise_type_summary.csv`, `analysis/V_daylight_walk_regime_days.csv`, `analysis/V_daylight_walk_regime_summary.csv`, `analysis/V_daylight_walk_regime_contrast.csv`, `analysis/V_daylight_walk_pair_days.csv`, `analysis/V_daylight_walk_pair_summary.csv`, `analysis/V_daylight_walk_pair_contrast.csv`

## U. Steps compensation

Does higher step load reduce future intake or raise TDEE?

7-day average steps → future 7-day intake: β=+0.103 (more steps → more eating). In falling phase, top 20% of 7-day steps (≥6311/day) is followed by lower 14-day intake (1974 vs 2088) and lower TDEE/RMR (1.0918 vs 1.1569). Closer to constrained-energy compensation than gravitostat benefit. The falling-phase high-step penalty is consistent with restriction archetype J: high-step cuts produce the worst TDEE/RMR outcomes.

Command: `python analysis/U_steps_compensation.py`
Artifact: `analysis/U_steps_compensation_regression.csv`, `analysis/U_steps_compensation_phase_thresholds.csv`

## AB. Running vs walking

If walk *sessions* predict RMR (finding AD), does exercise intensity matter, or just the discrete outing?

Does running suppress appetite more than walking at the same step count? (Exercise-induced anorexia hypothesis.)

121 running sessions (2013-2023), 1,640 walking sessions. All-eras step-matched comparison showed -338 cal same-day for running — but this was entirely era confounding (running clustered in the 2014-2016 restriction era).

Era-matched within 2014-2016 only, step-matched (90 pairs, mean 9650 steps): same-day cal difference -50, next-day +30, next-3d -46. All within noise. TDEE/RMR: run 1.072 vs walk 1.098 (-0.026). No exercise-induced anorexia detected. Running and walking produce equivalent intake effects at matched steps within the same era.

Command: `python analysis/AB_running_vs_walking.py`

# Tirzepatide

The set point (AG) and expenditure defense (K) describe the body's endogenous weight regulation. Tirzepatide intervenes pharmacologically. How does it interact with the system described above?

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

The PK model (F) shows tirzepatide reduces mean intake. But the set point (AG) predicts that appetite defense operates through binge *frequency*, not mean level. Does the drug reduce the mean or break the escalation pattern?

State machine: restriction (<1800), typical (1800-2399), high (2400-2799), binge (2800+).
- High → binge: 19.6% pre-drug → 4.8% on drug
- Binge → binge: 31.6% → 0%
- Post-restriction binge (7 days): 25.9% → 13.0%

The drug breaks escalation and persistence, not just mean level.

Command: `python analysis/G_tirzepatide_dynamics.py`
Artifact: `analysis/G_tirzepatide_transition_summary.csv`, `analysis/G_tirzepatide_rebound_summary.csv`

## AJ. Tirzepatide suppresses the expenditure defense

Finding K shows falling phases have elevated TDEE (+0.075 on TDEE/RMR ratio). Finding F shows 206 cal metabolic clawback on tirzepatide. But tirzepatide is a falling phase — shouldn't TDEE be *elevated*, not suppressed? Does the drug override the expenditure defense?

**Yes.** At matched fat mass (FM 61-71 lbs), tirz falling shows TDEE/RMR = 1.098 vs pre-tirz falling = 1.182 (Δ = -0.084). The drug eliminates most of the falling-phase TDEE bonus.

Regression (all days, TDEE/RMR ~ FM + falling + rising + on_tirz + falling×tirz):
- Pre-tirz falling effect: +0.075
- Tirz falling total effect: +0.012
- The falling×tirz interaction is -0.047, wiping out 63% of the falling bonus.

Within falling days only, controlling for FM: tirz coefficient = -0.095.

**The drug loses on both sides of the energy balance.** F shows -450 cal/day intake reduction. But the expenditure defense that normally accompanies weight loss (+75 on the ratio, roughly +150 cal/day) is suppressed to +12. The net deficit is smaller than the intake reduction alone would predict: -450 intake + 150 expenditure defense lost = net ~-300, consistent with the observed -250 cal/day net deficit and -1.64 lbs/month fat loss (vs -3.06 pre-tirz).

**Direct calorimetry confirms.** Pre-tirz falling calorimetry (n=8): measured RMR 90 cal above composition-predicted (elevated, K's defense). On-tirz calorimetry (n=2): measured RMR 193 cal below predicted (suppressed). The sign flip is visible in direct measurement.

**Mechanism.** The drug does not merely reduce appetite while leaving expenditure defense intact. It suppresses both arms: appetite (G: breaks binge escalation) and expenditure (this finding: eliminates the falling-phase TDEE bonus). This is consistent with GLP-1 agonists acting on hypothalamic energy homeostasis rather than purely on satiety circuits. The body's normal response to weight loss — burning more to defend the set point — is pharmacologically blunted.

Command: `python analysis/AJ_tirz_expenditure_defense.py`

## AL. Walking partially rescues the suppressed expenditure defense

AJ shows tirzepatide eliminates 63% of the falling-phase TDEE bonus (~125 cal/day lost). AD shows walk sessions raise RMR independently of composition. AH#5 says exercise and the set point are independent mechanisms. If all three are correct, walking should be additive — partially offsetting the pharmacologically suppressed defense.

**Walking predicts TDEE within the tirz era.** Walks (30d) vs TDEE/RMR: r=+0.21, controlling for FM and drug level (n=529 tirz days). During tirz falling phases (n=357), above-median walks (≥3 sessions/30d, n=211) show TDEE/RMR 1.130 vs below-median (n=146) at 1.082 — a +95 cal/day difference at the same fat mass (70.5 vs 70.8 lbs).

**The drug barely dampens the walk effect.** Regression (TDEE/RMR ~ FM + on_tirz + walks_30d + walks×tirz): walk coefficient pre-tirz = +0.00177 per session, on tirz = +0.00133 (interaction -0.00044, small). The drug preserves ~75% of the walk effect. At mean RMR: +3.5 cal/day per walk-session pre-tirz, +2.6 cal/day on tirz.

**Partial rescue, not full recovery.** The lost defense is 125 cal/day. At a realistic 15 walks/month: 15 × 2.6 = 39 cal/day recovered (~31% of the loss). Full recovery would require ~47 sessions/month (11/week) — unrealistic. But within the achievable range, walking is the largest single behavioral lever for expenditure during GLP-1 treatment.

**Pre-tirz comparison.** Walk effect on Kalman TDEE is actually *weaker* pre-tirz after controlling for FM (partial r=0.04 vs 0.19 on tirz). This is because pre-tirz, the endogenous expenditure defense dominates — walks are a rounding error on top of the body's own response. On the drug, the endogenous defense is suppressed, making walks the primary remaining lever.

Command: `python analysis/AL_walk_rescue_expenditure.py`

# Nulls and minor effects

## X. Temperature

Does body temperature track metabolic adaptation? Does tirzepatide override the thermostat?

Pre-tirzepatide: 14-day trailing intake → temperature r=+0.22. The 2600-3200 cal bin is 0.10°F warmer than 2000-2300. Temperature weakly tracks Kalman TDEE (r=0.04).

An earlier analysis found temperature correlating with tirzepatide blood level at r=+0.45. This was confounded by dose escalation timing (low doses in fall/winter, high doses in spring/summer). Direct calorimetry disproves the override hypothesis: RMR 1930→1750, temperature 97.7→97.5°F. Both dropped. The drug does not override the thermostat.

Command: `python analysis/X_temperature_phase.py`, `python analysis/P2_rmr_model.py`
Artifact: `analysis/X_temperature_daily.csv`, `analysis/X_temperature_phase_overlap.csv`, `analysis/X_temperature_phase_band_summary.csv`, `analysis/X_temperature_phase_regression.csv`

## AE. Sleep and energy balance (null)

Does sleep duration predict next-day intake, steps, weight change, or TDEE?

2,057 Samsung Health sleep measurements (2016-2026). Sleep → next-day calories: r = -0.01. Sleep → next-day steps: r = -0.02. Trailing 7-day sleep vs TDEE: r = +0.22 (confounded with era and fat mass). No trailing window of sleep hours predicts any energy balance variable.

Extreme short sleep (<5h, n=57): -80 cal and +757 steps next day — opposite the direction the sleep-obesity literature predicts. These are likely unusual days (travel, appointments), not chronic poor sleep.

The subject sleeps 3am-11am consistently (median 7.8h, std 1.5h). Weekend effect: +0.9h on Sat/Sun. Autocorrelation lag-1: r = 0.21. The schedule is too consistent to detect effects — there is not enough variation to separate signal from noise. Most sleep-obesity research uses self-reported cross-sectional data with far more variation.

Command: `python analysis/AE_sleep_null.py`

## AK. Sunlight exposure vs sleep duration

AE found sleep duration predicts nothing (r≈0 for all targets). But sleep barely varies (std 1.3h, CV 0.164). Possible sunlight exposure — hours awake between sunrise and sunset, computed from wake time + solar position at the subject's location — varies more (std 1.8h, CV 0.196) and is only weakly correlated with sleep (r=-0.13). Is sunlight a better predictor, and does it explain AD's walk-session → RMR effect?

**Sunlight beats sleep for daily energy balance, but both are weak.** Sunlight → same-day TDEE: r=+0.14 (sleep: +0.05). Sunlight → same-day steps: r=+0.14 (sleep: -0.02). After controlling for FM: sunlight→TDEE r=+0.11, sleep→TDEE r=-0.02. At trailing 30 days vs TDEE/RMR: sunlight r=+0.14, sleep r=+0.03 (controlling for FM). No sunlight window predicts intake (r=0.03-0.04 after FM control). Sunlight is a marginally better predictor than sleep, but neither explains meaningful variance in daily energy balance.

**AD's walk-session finding is not a sunlight confound.** Among 22 calorimetry measurements: walk sessions (30d) vs measured RMR: r=+0.73. Sunlight (30d) vs RMR: r=+0.53. Controlling for sunlight, walks retain r=+0.66. Controlling for walks, sunlight drops to r=+0.36. In LOO-CV: walks alone RMSE=116, sunlight alone RMSE=144, walks+sunlight RMSE=121. Adding sunlight to walks adds nothing — walks are the signal, sunlight is along for the ride.

**Sleep shows a surprising raw correlation with RMR** (r=-0.71: less sleep = higher RMR), but this mediates the walk effect — active periods have both more walks and less sleep. After controlling for walks + sunlight + composition: sleep partial r=+0.21. The full model (walks + sunlight + sleep + expected_rmr) achieves CV RMSE=115, only marginally better than walks alone (116).

Command: `python analysis/AK_sunlight_exposure.py`
Artifact: `steps-sleep/sunlight.csv`

# Diet experiments

## O. Epoch analysis, potato diets, and travel

Do the 22 manually annotated diet epochs in `intake/diet_epochs.csv` separate meaningful regimes? Do potato diets have special metabolic properties?

Travel binges: 2716 cal/day, 32.1% binge rate, +9.5 lbs fat. Keto: 2163 cal, 110.5g protein, 2.1% binges. Potato diet (4 epochs, 69 days): 1930 cal, 48.6g protein, 0% binges, -7.3 lbs fat.

Potato diets: TDEE/RMR 1.168 before, 1.159 during (unchanged), 1.098 after. Post-potato rebound: 2629 cal/day, 19.6% binge rate. Potatoes are a monotone binge-suppressing cut, not a uniquely expenditure-preserving intervention.

**Jordan trips.** Two month-long travel periods tested whether high intake produces compensatory TDEE increases. 2015: intake 3262 cal/day, TDEE 2071 (-35 vs pre-trip), +10.3 lbs fat. 2019: intake 2529, TDEE 2021 (-21 vs pre-trip), +4.4 lbs fat. No compensatory TDEE surge. The trips are binge windows without an energy-out bonus.

Command: `python analysis/O_diet_epoch_analysis.py`, `python analysis/P_jordan_trip_analysis.py`
Artifact: `analysis/O_diet_epoch_summary.csv`, `analysis/O_diet_epoch_family_summary.csv`, `analysis/O_potato_epoch_window_summary.csv`, `analysis/O_potato_epoch_contrast.csv`, `analysis/P_jordan_trip_summary.csv`, `analysis/P_jordan_trip_delta.csv`

# Narrative case

The set point and expenditure defense together explain the 2013 inflection — the moment the decade of regain began.

## AC. The 2013 inflection

Fat mass bottomed at 17 lbs (Oct 2013) and rose every year for a decade. What triggered the regain — metabolic adaptation forcing intake up, or behavioral momentum?

Bottom (2013 Jul-Sep): intake 1864, TDEE 2021, ratio 0.977, binge rate 1%. Regain (2014 Jul-Dec): intake 2359, TDEE 2058, ratio 1.026, binge rate 20%.

TDEE recovered +38 cal from bottom to regain. Intake increased +495 cal. Intake outran TDEE recovery by 457 cal/day. The metabolic adaptation was easing — ratio went from 0.98 to 1.03. But binge frequency exploded: 0% (Jul-Aug 2013) → 13% (Oct) → 20-23% (Nov-Dec) → sustained 12-26% through all of 2014. By 2014-H2, binges were running at 20%.

The trigger was behavioral, not metabolic. The body's expenditure defense was releasing. Food noise — the variance, the binge clustering — took hold at FM=17 and never let go for a decade. The set point model (AG) explains the mechanism: at FM=17, the set point had stalled near 18.6 lbs (AH test #6, floor effect). Any weight above the floor put FM below the set point, triggering the binge frequency gradient. The expenditure arm (K) was simultaneously releasing — TDEE/RMR rose from 0.98 to 1.03 — but this was too little, too late. The appetite arm had already taken over.

Command: `python analysis/AC_inflection_2013.py`
