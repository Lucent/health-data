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

# Metabolic rate adjustment

The year-by-year table shows TDEE falling from 2661 to 2044 as fat mass dropped from 106 to 25 lbs — partly from losing metabolically active tissue, partly from metabolic adaptation. TDEE barely recovered during a decade of regain. This raises the question: at the same body composition, does TDEE differ depending on direction?

## K. TDEE hysteresis

At the same fat mass, does TDEE differ depending on whether weight was reached from above (falling) or below (rising)? [Leibel et al. (1995)](https://doi.org/10.1056/NEJM199503093321001) demonstrated this in a metabolic ward (n=18, ±10% body weight). This dataset replicates the finding free-living with 1,482 FM-matched pairs.

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

## B. Weekend fasting as metabolic adjustment microcosm

Do acute caloric deficits produce lasting fat loss? Seven consecutive Sat-Sun 36-hour fasts (Oct-Nov 2019). Mean deficit per fast: 3,300 cal. Kalman FM change Fri→Mon: -0.76 lbs (80% of expected). Kalman FM change Fri→Fri+7: +0.08 lbs. Post-fast weekday intake: 2,369 cal (78 cal/day below pre-fast 2,447). Zero compensatory overeating — yet the deficit vanishes within a week. The Kalman TDEE drops on fasting days and stays low through the following week, absorbing the deficit through reduced metabolic rate rather than increased intake.

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

K and B show the body adjusts metabolic rate by direction. Is there a single defended fat mass driving both the metabolic and appetite responses? [Speakman (2011)](https://doi.org/10.1242/dmm.008698) proposes dual intervention points with a "zone of indifference." The data finds no such zone — the response is graded and continuous.

## AG. A hidden, moving fat mass set point

Binge rate (days > TDEE + 1000 cal) reverse-engineers the body's defended weight. **The set point is an exponential moving average of fat mass with a ~50-day half-life.** Distance below it predicts 90-day binge rate at r = -0.60, partial r = -0.63 controlling for absolute FM.

**It tracks fat mass, not weight** (r = -0.60 fat vs -0.53 total weight vs -0.52 scale at every HL). Consistent with [Kennedy (1953)](https://doi.org/10.1098/rspb.1953.0009) and the leptin system. Gravitostat: null (N: r=+0.05, wrong direction).

**It moves, not fixed** (r = -0.60 moving vs +0.20 best fixed). After ~150 days at a new weight, the set point has 87% adapted. The literature has no comparable estimate — clinical sources cite vague "1-6 years" with no measurement. [Lowe et al. (2007)](https://doi.org/10.1002/eat.20405) showed weight suppression predicts binge frequency in bulimia, but used a fixed proxy.

**Sigmoid binge response.** Baseline ~3% at/above set point, rising to ~15% at 5+ lbs below. Continuous gradient, no hard threshold.

| Distance from SP | Binge rate |
|---|---|
| 5-7.5 lbs below | 14.8% |
| 2.5-5 lbs below | 11.4% |
| 0-2.5 lbs below | 6.0% |
| 0-2.5 lbs above | 3.3% |
| 2.5-5 lbs above | 2.8% |
| 5-10 lbs above | 1.3% |

**Tirzepatide overrides the set point's eating pressure.** At matched SP distance, off-tirz binge rate is 3.9%, on-tirz is 0.5% (partial r = -0.23 controlling for distance). The 2025 row of the trajectory table (FM 70, SP 74, 4 lbs below SP — should produce ~8-10% binges, actual 0%) demonstrates the drug's independent suppression. Injection cycle visible: 1640-1740 cal on days 0-1, rising to 2140-2220 by days 3-5 (~500 cal/day sawtooth).

**On-drug set point adaptation: not frozen, but uncertain.** The 165-day on-drug HL from binge rate was an artifact of sparse data. With mean surplus (AP), the on-drug HL surface is flat: r = -0.92 at 40d, r = -0.91 at 120d. The SP adapts on-drug (HL 40-120d vs 40-55d pre-tirz), but the exact rate cannot be resolved.

Command: `python analysis/AG_binge_set_point.py`, `python analysis/Y_set_point_intake_tests.py`, `python analysis/C_binge_analysis.py`, `python analysis/E_weekly_invariance.py`, `python analysis/H_phase_binge_landscape.py`
Artifact: `analysis/H_phase_binge_summary.csv`, `analysis/H_phase_binge_distance_bins.csv`, `analysis/H_phase_binge_model_auc.csv`

## AH. Set point mechanics

**Apparent ratchet: retracted.** An asymmetric model (HL_down=20d, HL_up=100d) improved binge-rate r from -0.60 to -0.69 on the full dataset. However, **this asymmetry does not replicate on 2014+ natural-dynamics data** (best asymmetric Δr = 0.0004 vs symmetric, AM). The asymmetry was driven by the 2011-2013 willpower period: the EMA mechanically chased a rapid externally-forced descent, creating the appearance of fast downward adaptation. On data where intake responds naturally to the set point, the set point is symmetric at ~45 days in both directions.

**Binge size is constant; the set point tilts the daily mean.** Binge surplus ~1427 cal regardless of SP distance (r = -0.02). Non-binge days: r = -0.32 between SP distance and surplus (+80 cal at 5+ below SP, -481 cal at 5+ above). The set point modulates background drift on every day, not just discrete events — confirmed at a deeper level by AP.

**Restriction runs ending above SP stick; below SP, they rebound.** 237 runs (3+ days, cal < TDEE-200). Above SP at end: -1.2 lbs/30d. Below: +0.4 lbs/30d. SP distance vs rebound: r = -0.48.

**Exercise independence.** Walking raises RMR (AD) but does not change the SP half-life (40d for both high-walk and low-walk periods).

**Floor effect.** SP minimum: 18.6 lbs FM (Nov 2013), near essential body fat. Post-floor: FM 19→23 over 6 months, binge rate 8.3%. This is where the 2014 inflection began (AC).

Command: `python analysis/AH_set_point_properties.py`

## BK. A first control-stock model of willpower

The `2011-2013` crash from ~260 lbs to ~180 lbs does not look like natural set-point-following intake. A first attempt to model this as a **depletable control stock** says the data do contain a latent "willpower" signal, but one smooth reservoir is not enough.

Model:
- `pressure_t = 55 * (SP_t - FM_t)`
- observed surplus = baseline + pressure - control exertion
- control comes from a latent stock that depletes when exerted and partially recovers over time

**The 45-day SP timescale reappears again.** In the coarse search over SP half-life, baseline offset, control capacity, recovery, and depletion, the best run again used **SP HL = 45d**. This same ~45-50d timescale has now appeared in binge-rate set-point fits (AG), mean-surplus fits on 2014+ natural dynamics (AM/AP), intake-free FM variants, and this first latent-control attempt. The most defensible interpretation is that `~45d` is the timescale of the **expressed appetite-pressure layer** rather than necessarily the whole structural defended state.

Best coarse run:
- SP HL = 45d
- baseline = -400 cal/day
- control capacity = 1200 cal/day-equivalent
- recovery = 2%/day
- depletion cost = 0.4 × exerted control
- extra low-demand recovery = +20/day

What it found by phase:
- `2011-2012`: huge required control (`528 cal/day` mean) and large mean shortfall (`467 cal/day`) while FM fell `99→32 lbs`
- `2013`: required control drops to `132 cal/day`, shortfall `70`, FM `32→20`
- `2014-2015`: required control near zero (`23-28 cal/day`), stock mostly refilled, FM regains from `20→47`
- `2016`: positive pressure returns, required control rises to `123 cal/day`, shortfall `89`

**Interpretation:** the model does "spot" real latent control pressure in the history. It supports the idea that early restriction involved sustained suppression of a strong biological drive, and that regain bouts can be understood partly as periods where required control exceeds available stock.

**But one stock is too simple.** Even the best coarse run still had `25.6%` "impossible" days where required control exceeded available stock by more than `100 cal/day`. So the `2011-2012` period is too extreme to be explained by one smooth depletable reservoir alone. A better next model would likely separate:
- `diet intent / rule commitment`
- `control fatigue / depletion`

Command: `python analysis/BK_control_stock_model.py`
Artifact: `analysis/BK_control_stock_daily.csv`

## AI. The metabolic channel operates on a 9-day timescale

The eating channel adapts with a 50-day half-life (AG). **The metabolic channel is faster: optimal HL = 9 days** (partial r = +0.39, p = 0.03 at n_eff=17). The correlation is monotonically stronger at shorter half-lives, then plateaus. K's hysteresis (+179 cal falling vs rising, p < 10^-22) confirms the same signal through FM-matched pairs.

| HL (days) | r (eating) | r (metabolic, partial) |
|---|---|---|
| 9 | — | +0.39 |
| 30 | -0.59 | +0.32 |
| 50 | -0.60 | +0.27 |

**Not a Kalman artifact.** Lagging SP distance by 7-60 days: partial r decays from 0.39 → 0.35 → 0.23 → 0.15. A purely mechanical artifact would vanish by lag 7-14; this persists to 30+ days.

**The two channels operate on different timescales** — metabolic rate is a sprint (~9 days), eating pressure is a marathon (~50 days). This explains the 2013 inflection (AC): the metabolic boost had faded while eating pressure was still fully engaged. It reframes ["persistent metabolic adaptation" in Biggest Loser contestants](https://doi.org/10.1002/oby.21538): if the metabolic channel adapts in ≤10 days, persistent RMR suppression after 6 years means the subjects remained chronically below a set point that itself had moved.

Command: `python analysis/AI_expenditure_arm_timescale.py`
Artifact: `analysis/AI_expenditure_arm_sweep.csv`

## AM/AP. The set point controls mean daily energy balance (r = -0.92)

**Methodological note: 2014+ data only.** All parameters in this section are derived from January 1, 2014 onward — the date when 90-day binge rate first exceeded 5%, indicating sustained set-point-driven eating pressure. FM had bottomed at 17 lbs (Oct 2013) and the first binge occurred the next day. The 2011-2013 period of aggressive willpower-driven restriction (1200 cal/day, FM 106→17 lbs) is excluded: intake was externally controlled and did not respond to the set point. Including 2011-2013 inflated the per-lb pressure from -27 to -45 cal/lb and created a spurious 3-5x asymmetric "ratchet" (the EMA mechanically chasing a rapid externally-forced descent). On 2014+ natural-dynamics data, the set point is symmetric and the pressure is weaker.

AG discovered the set point using binge rate (r = -0.60). AP tested 24 functional forms and found **the true signal is 90-day mean surplus** (daily calories minus TDEE, averaged): **r = -0.92** on 2014+ data. This beats every binary threshold and every transform.

**Identifiability caveat: half-life and outcome window partially trade off.** A 2D sweep over symmetric SP half-life and surplus lookback window produces a broad diagonal ridge rather than a sharp point optimum. Examples from the 2014+ pre-tirz data: HL=20d with 45d mean surplus gives `r=-0.939`; HL=40d with 75d gives `r=-0.929`; HL=45d with 90d gives `r=-0.926`; HL=75d with 120d gives `r=-0.912`. That means the mean-surplus fit alone does **not** uniquely identify a structural 45d half-life. What it identifies is a slow appetite-pressure timescale. The argument for `~45-50d` specifically comes from convergence with AG's binge-rate fit and the intake-free P3/P4 variants, not from the `r=-0.92` surplus regression in isolation.

**The 45-day half-life survives an intake-free fat mass estimate.** The Kalman filter uses logged intake in its process model, raising a circularity concern. **P3 derives fat mass from weight interpolation alone — no intake data enters the FM estimate.** Using P3's intake-free FM on 2014+ data, the set point still predicts mean surplus at **HL = 50d (r = -0.73)**. The half-life is the same; only the correlation strength drops because P3's FM is noisier between weigh-ins. **The 45-50 day timescale is a property of fat mass dynamics, not a Kalman filter artifact.** It appears in intake-free weight interpolation, Kalman-filtered FM, binary binge rate, and continuous mean surplus — four independent measurement combinations converging on the same number. In other words: the mean-surplus surface has a ridge, but the ridge intersects the intake-free and binge-rate estimates in the same 45-50d band.

The set point tilts the entire daily distribution. It controls **how often** you eat above maintenance (r = -0.91) and **how deeply** you restrict on deficit days (r = -0.82), but **not how large** individual overshoots are (conditional magnitude r = -0.40; adding magnitude to rate improves R² by 0.01).

**Per-lb pressure: ~55 cal/day per lb below SP** (from binned data across all gap sizes, -10 to +20 lbs, remarkably constant). Five pounds below = ~275 cal/day of upward pressure. The regression-derived estimate of -27 cal/lb is diluted by autocorrelation; the binned estimate is more reliable.

**Symmetric half-life: 45 days.** The plateau spans 40-50d. **Asymmetric models show zero improvement on 2014+ data** (best asymmetric Δr = 0.0004 vs symmetric). The 3-5x ratchet reported in AH was driven by the excluded 2011-2013 period. On natural-dynamics data, the set point adapts at the same rate in both directions.

**Metabolic channel asymmetry (survives the exclusion).** The metabolic adjustment (AI: partial r = +0.39 overall) is **4.2x stronger when FM is above SP** (above-SP subset partial r = +0.52) than below (+0.12). The body actively helps fat loss by burning more when above the set point, but barely adjusts when below. This asymmetry is measured on 2014+ data using the metabolic channel's own HL (9d) and is not affected by the 2011-2013 exclusion.

**Restriction run prediction: r = -0.48 (AH) to -0.72 (AM, stricter definition).** Robust across definitions. Runs ending above SP continue losing (-1.2 lbs/30d); below SP, they rebound (+0.3 lbs/30d).

Command: `python analysis/AM_lipostat_sensitivity.py`, `python analysis/AP_overshoot_shape.py`, `python analysis/AN_ratchet_profile.py`

## AQ. Tirzepatide set point coverage — what the drug buys and how it erodes

Using the mean surplus metric (AP), we can decompose the drug's effect into set-point-dependent and set-point-independent components, and quantify coverage in lbs and cal/day for each arm.

**The set point is not frozen on the drug, but the on-drug HL is uncertain.** AG estimated the on-drug SP half-life at 165 days using the 1000-cal binge threshold — but binge events are too rare on-drug to track SP adaptation reliably. With mean surplus, the optimal on-drug HL is 40 days (r = -0.92), but the surface is flat: r = -0.91 at HL=120d. The data cannot distinguish 40d from 120d on-drug. What is clear: the 165d estimate was an artifact of sparse binge data, and the SP is not frozen — it is adapting, with HL somewhere in the 40-120d range (vs 40-55d pre-tirz).

**Appetite arm: the drug operates independently of SP distance.** Regression (90d mean surplus ~ SP distance + FM + effective level): SP distance = -77 cal/lb [CI excludes zero], effective level = +5.4 cal/unit [CI includes zero]. Adding a distance × drug interaction improves R² by < 0.001. The drug does not attenuate the per-lb pressure from the set point — it adds a **constant offset** of ~49 cal/day (Model A, on-tirz intercept). Drug equivalence: 1 unit of effective level offsets only 0.1 lbs of SP distance [CI: -0.0, 0.1]. At current levels (eff = 6.9): the drug suppresses 37 cal/day against a gap pressure of -102 cal/day from the 1.3 lb gap, leaving a net -65 cal/day (continued slow loss).

**Metabolic channel: the drug suppresses the body's calorie-burning boost.** Regression (TDEE residual ~ SP distance + FM + effective level): SP distance = +19 cal/lb (body burns more when above SP, accelerating loss), effective level = -10.2 cal/unit (drug reduces this burning). Each unit of drug level suppresses 10 cal/day of the metabolic boost that normally helps fat loss. At eff = 7: -71 cal/day of lost calorie burning. The metabolic drug equivalence is larger: 1 unit = 0.5 lbs of SP offset. With the distance × drug interaction (R² +0.009): pre-tirz gets +18.4 cal/lb of metabolic boost; each unit of drug reduces this by 8.7 cal/lb. At eff = 7, the per-lb boost drops from +18 to -43 cal/lb — the drug *reverses* the metabolic cooperation into metabolic *resistance*.

**Coverage table — drug effect on each channel at various drug levels:**

| Eff level | Appetite (less eating) | Metabolic (less burning) | Combined |
|---|---|---|---|
| 2 | +11 cal/day | -20 cal/day | -9 cal/day |
| 4 | +21 | -41 | -20 |
| 6 | +32 | -61 | -29 |
| 8 | +43 | -81 | -38 |
| 10 | +54 | -102 | -48 |
| 12 | +64 | -122 | -58 |

The appetite column is the drug's direct effect on eating (positive = less restriction, helping avoid regain). The metabolic column is the drug suppressing the body's natural calorie-burning boost during weight loss (negative = less burning, slowing loss). The drug reduces appetite (pro-loss) while also reducing the metabolic rate increase that normally accelerates loss (anti-loss). The net is still pro-loss because the appetite effect on total energy balance is larger.

**Tachyphylaxis erosion.** At 32-week half-life, both effects degrade: metabolic suppression from -102 cal/day (week 0) to -51 cal/day (week 32) to -25 cal/day (week 64); appetite reduction degrades proportionally. The SP adapts simultaneously (HL = 45d symmetric), shrinking the gap. At the current ~1 lb gap, the SP reaches FM within weeks regardless of drug level.

**Injection cycle.** The weekly sawtooth is visible in daily surplus: injection day (day 0) -525 cal, rising to +51 cal by day 5. A 576 cal/day swing. Effective drug level drops from 9.6 (day 0) to 5.2 (day 5) through the sawtooth.

Command: `python analysis/AQ_tirz_set_point_coverage.py`

## AW. Reconciling set point and tachyphylaxis

AQ decomposed the drug into set-point-dependent and independent components but couldn't separate tachyphylaxis from set point adaptation. AW resolves this by holding the set point parameters fixed at their pre-tirz values (HL=45d, -27 cal/lb, from AM 2014+ data) and fitting only the drug parameters (per-unit appetite effect and tachyphylaxis half-life) to the 529 on-drug days.

**Reconciled formula:** `daily_intake = TDEE + (-27 × SP_gap) + (-25 × effective_level)` where effective_level = blood_level × exp(-0.069 × weeks_on_current_dose).

**Drug appetite effect: -25 cal per unit effective level** (vs F's original -49). F's estimate was confounded — it attributed some of the set point's eating pressure reduction (as SP converged toward FM) to the drug. With the SP pressure held fixed, the drug's independent contribution is half as large.

**Tachyphylaxis half-life: 10 weeks** (vs F's original 32 weeks). The drug loses half its effectiveness every 10 weeks, not 32. This is consistent with the trough-day intake trend at stable 12.5mg: +12 cal/week observed, +7.4 cal/week predicted by the model (62% explained). The remaining +4.5 cal/week residual trend may reflect behavioral drift or nonlinear tachyphylaxis.

**The injection sawtooth is captured but understated.** The model predicts a ~75 cal amplitude between day 0 and day 5; the actual amplitude is ~580 cal. The per-unit effect may be nonlinear — stronger at peak blood levels (day 0-1) than the linear model assumes.

**Sensitivity to SP half-life.** The drug parameters are stable across SP HLs of 30-80d: appetite effect ranges -24 to -26 cal/unit, tachyphylaxis HL 6-13 weeks. RMSE is flat at 376-386 cal/day (against daily intake noise of ~500 cal). The drug parameters are not sensitive to the exact SP HL.

**Distinguishing tachyphylaxis from set point adaptation.** At stable 12.5mg dose (28 weeks), trough-day (day 4-6) intake rises at +12 cal/week — the tachyphylaxis signature. Set point adaptation alone would predict *falling* trough intake as the gap closes. The rising trend confirms **tachyphylaxis dominates set point adaptation during stable dosing**, with the SP pressure (-27 × ~2 lb gap = ~54 cal) too weak to offset the tachyphylaxis-driven erosion.

## AX. Clean drug effect identification

AW's reconciliation was structurally confounded: the set point was recomputed on the full series (including drug-era FM decline), absorbing drug effect into SP adaptation and producing artificially weak drug coefficients (-25 cal/unit, 10-week tachy HL). AX fixes this:

**Set point frozen from pre-drug data.** The EMA is computed on pre-drug FM only, then propagated forward using the same HL=45d rule. During the drug era, the SP chases observed FM (retrospective decomposition) but was not fitted to drug-era data. At drug start: SP = 83.5, FM = 84.2, gap = -0.7 lbs (essentially at equilibrium). At end: SP = 62.5, FM = 59.9, gap = +2.6 lbs (FM below SP — the drug pushed past the defended weight).

**Drug appetite effect: -74.5 cal per unit blood level.** Identified from within-week injection cycle variation (week fixed effects absorb SP pressure, tachyphylaxis, and all slow trends). This is 1.5x F's -49 and 3x AW's confounded -25. The naive cross-sectional estimate is -37, confirming that cross-week confounds dilute the signal by half.

**Tachyphylaxis: 35-week half-life, cumulative exposure.** Cumulative weeks on drug (r = -0.57) fits much better than per-dose-reset (r = -0.22). The per-dose reset model is not portable — this subject's 20 weeks at 5mg creates deep tachyphylaxis that a 4-week trial escalation never develops. Cumulative exposure avoids this problem.

**Metabolic arm: zero within-week.** TDEE residual does not respond to blood level variation within the injection cycle (+0.1 cal/unit, r = +0.04). The metabolic effect operates on longer timescales than the 7-day cycle can capture — consistent with AI's 9-day HL and K's 90-day phase classification.

**Forward validation on this subject: RMSE = 364 cal/day, r = +0.41.** The model decomposes each on-drug day into TDEE + SP pressure (27 × gap) + drug effect (-74.5 × blood × tachyphylaxis). Residuals by injection day range from -215 (day 0, model undershoots suppression) to +201 (day 5, model undershoots rebound) — the actual sawtooth is ~100 cal wider than predicted on both ends.

**Preliminary SURMOUNT comparison (AV).** Using AX's parameters in a forward simulator (not fully validated — uses retrospective SP, appetite-only energy balance):

| Dose | Simulated | SURMOUNT-1 | Δ | SURMOUNT-2 (diabetic) |
|---|---|---|---|---|
| 5mg | -8.3% | -16.0% | +7.7% | — |
| 10mg | -15.0% | -21.4% | +6.4% | -12.8% |
| 15mg | -20.3% | -22.5% | +2.2% | -14.7% |

The 15mg non-diabetic prediction is within 2.2 percentage points of the published result with zero trial-fitted parameters. The model overpredicts diabetic weight loss (simulated -21% vs published -15% at 15mg), consistent with the drug's appetite effect being partially diverted to glucose control in T2D — a known mechanism that the non-diabetic-derived -74.5 cal/unit does not account for. The 5mg arm undershoots because 35-week tachyphylaxis erodes the low dose too quickly relative to trial subjects who escalate past 5mg within 4 weeks.

## AY/AZ. The set point only adapts when weight is stable (SmoothLatch model)

The 45d EMA predicts ~97% SP adaptation over 36 weeks on drug, leaving minimal regain (+2.8% simulated vs +14% published SURMOUNT-4). No single EMA half-life, two-component mixture, rolling mean, floor constraint, or pressure function fits both this subject and trial regain data (AZ: exhaustive search, 268 configurations).

**The SmoothLatch model resolves this.** The set point only adapts toward fat mass when FM has been stable — within ±3 lbs of a reference level for at least 14 consecutive days. When FM is changing rapidly, the SP freezes. When FM holds steady, the SP approaches it at rate 0.01/day. Pressure: 55 cal/lb (from binned data across all gap sizes, remarkably constant from -10 to +20 lbs).

**The SmoothLatch and EMA are observationally equivalent on this subject.** Both fit at r ≈ -0.94 on 180-day surplus. An EMA at HL=80d on 180-day surplus gives r = -0.94; the SmoothLatch(tol=3, hold=14, rate=0.01) on 180-day surplus gives r = -0.94. On 2014+ pre-tirz data: SmoothLatch r = -0.90, EMA 45d r = -0.92. **The models cannot be distinguished from this subject's slow-gain natural dynamics alone** — FM changed at ~6 lbs/year, always within the latch tolerance, so the SP adapted continuously in both models.

**They diverge on rapid-loss scenarios where only trial data can break the tie:**

| Scenario | EMA (45-80d) | SmoothLatch | Published |
|---|---|---|---|
| **SURMOUNT-4** regain (52wk post-drug) | +2.8% | **+14.1%** | **+14.0%** |
| **STEP-1** regain (52wk post-semaglutide) | ~+2% | **+9.3%** | **~+10%** |
| This subject 2014+ (r, 180d surplus) | -0.94 | -0.94 | — |

The SmoothLatch matches SURMOUNT-4 within 0.1 percentage points and STEP-1 within 1 percentage point. **The EMA undershoots regain by 5-10x.** Trial data decisively favors the SmoothLatch.

**Why it works:** during drug-driven loss (~1.2 lbs/week), FM is never within ±3 lbs of any reference for 14 consecutive days. The SP freezes — it retains the pre-treatment defended weight. At discontinuation, the gap is ~24 lbs, generating +1320 cal/day initial surplus. During regain, FM eventually stabilizes enough for the SP to begin latching onto the higher weight, gradually closing the gap over months.

**Testable prediction: regain depends on the speed of loss, not just the amount.** At the same total fat loss, faster loss produces more regain because the SP has less time to latch:

| Loss rate | Weeks | Total loss | Gap at stop | Predicted regain |
|---|---|---|---|---|
| 0.5 lb/week | 100 | ~12% | 2 lbs | +0.7% |
| 1.0 lb/week | 53 | ~13% | 5 lbs | +2.3% |
| 1.5 lb/week | 35 | ~13% | 11 lbs | +4.7% |
| 2.0 lb/week | 26 | ~13% | 17 lbs | +7.0% |

This differentiates the SmoothLatch from any EMA: **the EMA predicts regain depends only on time, while the SmoothLatch predicts it depends on speed.** Consistent with clinical observations that surgical and very-low-calorie patients regain more than gradual losers.

**Biological plausibility.** Leptin receptor adaptation, adipocyte number remodeling, and hypothalamic circuitry all require sustained exposure to a new adiposity level. Rapid FM changes may not trigger the cellular mechanisms that reset the defended weight. The 14-day stability requirement and 3-lb tolerance window are consistent with the timescale of leptin signaling equilibration.

**This is the preferred model.** It costs nothing on subject fit (r = -0.94 on 180d surplus, tied with best EMA), correctly predicts two independent trial regain datasets, and generates a novel testable prediction. The EMA remains valid as a computationally simpler approximation for slow-changing weight scenarios.

Command: `python analysis/AY_sp_from_regain.py`, `python analysis/AZ_sp_model_search.py`

Command: `python analysis/AX_drug_model.py`, `python analysis/AV_surmount_simulation.py`

# Intake characterization

The set point (AG/AM) controls mean daily surplus at r = -0.92. But what shapes day-to-day intake beyond the set point? Several dietary hypotheses predict intake modulation that should be visible in this dataset.

## D. Food noise as intake variance

The [food noise essay](https://lucent.substack.com/p/craving-food-noise) distinguishes hunger (daily deficit) from food noise (persistent awareness proportional to distance below defended weight). Food noise manifests as intake variance (erratic eating, clustering binges), not as a constant upward pressure proportional to set point distance:

- Restriction duration → rebound: r=-0.065. No compounding effect.
- Intake variance below set point: std=615, CV=0.278 vs at set point std=528, CV=0.244. Tirzepatide reduces CV from 0.24-0.28 to 0.19-0.20 (25-30% reduction).

The distance-to-intake gradient (-30 cal/day per kg below set point) was measured here using a symmetric SP and 180-day rolling mean, which diluted the signal. Binned data across all gap sizes (AM) shows a consistent **55 cal/lb (121 cal/kg)** — close to the [100 cal/kg claim](https://doi.org/10.1002/oby.21653) and consistent with the food noise essay's estimate of ~100 cal/kg. The essay's gravitostat mechanism is wrong (N: null), but the magnitude of eating pressure it described is approximately correct.

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

The mechanism is likely finding B: high-variance periods include fasting days that transiently raise TDEE through metabolic adjustment. The variance itself is not protective — it proxies for intermittent acute deficits. Diet epochs confirm: weekend fasting (CV = 0.68) lost weight despite moderate mean; COVID lockdown (CV = 0.18, very consistent eating) gained despite similar mean.

The effect is small — adding variance to a surplus-only model reduces RMSE from 1.019 to 0.998 (2%). But the sign is the opposite of "yo-yo dieting is metabolically damaging." In this dataset, consistent overeating is more fattening than variable overeating at the same mean. (Cf. the [Biggest Loser reanalysis](https://doi.org/10.1002/oby.23308) which found 500 cal/day metabolic adaptation from rapid cycling — but that study's subjects were in continuous severe deficit, not the intermittent variance measured here.)

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

The set point controls appetite (AG) and the metabolic rate adjustment depends on trajectory direction (K). But finding AD shows a third lever — deliberate walking — that modulates RMR independently of both. This section characterizes how activity interacts with the energy balance system.

## AD. Walk sessions predict RMR

Exhaustive sweep of trailing dietary, activity, and sleep features at 7 windows (7-60 days) against 23 Cosmed Fitmate calorimetry measurements with real Samsung Health data (2016+, dropping 2 pre-Samsung measurements with backfilled steps). Leave-one-cluster-out cross-validation with ridge regression. 12 independent clusters.

**Null results.** Baseline (expected_rmr only): CV RMSE = 183. Dietary features (calories, protein %, carbs, fat, sodium) at every window: best CV RMSE = 169, no improvement over the Fitmate noise floor (~170 cal). Step counts at every window: best CV RMSE = 151 (steps_14d). Sleep at every window: negative R² throughout. Strength training count: no signal.

**Walk sessions predict RMR.** The count of deliberate walks Samsung Health logged as exercise in the prior 30 days predicts RMR at CV RMSE = 116 (R² = 0.49), well below the Fitmate noise floor. The coefficient: +14 cal RMR per walk session. Going from 3 walks/month (Sep 2025, RMR 1750) to 33 walks/month (May 2022, RMR 2292) corresponds to a 420 cal/day difference — at nearly identical body composition (FM 66-68 in 2022 vs 68 in 2025).

Walk sessions beat walk minutes (CV RMSE 135) and total steps (CV RMSE 179) at the same 30-day window. The count of distinct outings matters more than total duration or total movement.

**Paired daylight walks.** 1,636 walking sessions (type 1001). Post-2020 noon-7pm year-round: 149 paired-daylight-walk days, 153 single-daylight-walk days. Step-matched paired walks vs other days: -57 cal lower future 14-day intake, +0.019 higher future TDEE/RMR. Single walks: +17 cal (no intake benefit), +0.027 TDEE/RMR. The paired-walk structure matters more than raw step count.

**Not season.** Walk sessions and summer are correlated at r = 0.83 (more walks in warm months). But: (1) Season alone barely helps — `is_summer` reduces CV RMSE from 183 to only 175, while walk sessions reduce it to 116. (2) Adding season on top of walk sessions adds nothing (CV RMSE 116 → 116). (3) Controlling for season, the partial correlation between walk sessions and RMR is still r = 0.47 (steps drop to r = 0.13). (4) Within winter only (n=12), walk sessions still predict RMR at r = 0.63. (5) Within 2022-2023 where FM was nearly constant (66-75 lbs, n=18), walk sessions vs RMR: r = 0.68, partial controlling for expected_rmr: r = 0.62. (6) Within the May-Jun 2022 cluster alone (n=8, same FM=66-68, same season), walk sessions vs RMR: r = 0.69.

**Not tirzepatide, not body composition.** Only 2 of 23 measurements are on the drug. Adding tirz_level to the model doesn't change the walk session coefficient (CV RMSE 117 vs 116). Walk sessions are uncorrelated with fat mass (r = -0.02) and negatively correlated with expected_rmr (r = -0.28) — the effect is not mediated by composition changes.

**Interpretation.** Deliberate sustained walks (typically 20+ min continuous, enough for Samsung Health to log as an exercise session) raise resting metabolic rate in a way that total step count — which includes all incidental movement — does not. This contradicts Pontzer's [constrained total energy expenditure model](https://doi.org/10.1016/j.cub.2015.12.046), which predicts TEE plateaus at higher activity levels. The effect here is on *resting* metabolic rate measured by calorimetry, not TEE redistribution. The mechanism is consistent with [NEAT upregulation](https://doi.org/10.1126/science.283.5399.212): intentional exercise sessions may activate a metabolic afterburn that persists at rest, while shuffling around the house does not. The effect size (+14 cal/session, ~420 cal/day at 30 vs 3 sessions/month) is large but consistent across subgroups. The remaining confound is that with 23 measurements clustered by month, we cannot fully exclude an unmeasured seasonal factor that drives both walking and RMR.

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

## BW. Walking fitness from exercise live-data heart rate

The earlier hourly `tracker.heart_rate` fallback was the wrong source for per-workout heart-rate. Samsung stores the chart shown in the exercise UI in `jsons/com.samsung.shealth.exercise/*/*.com.samsung.health.exercise.live_data.json`, attached to exercise sessions in `com.samsung.shealth.exercise.*.csv`.

Using the correct source gives 448 usable walking sessions with live heart-rate, pace, and distance from 2016-04-18 to 2026-03-10. The early 2016-2017 regime is on a different scale (walking HR median 153 in 2016, 134 in 2017) and should not be pooled with the later watch-era data.

Restricting to the modern 2021+ era, the pace relation is:

`walk HR = 81.48 + 5.89 * kph`

Pace-adjusted yearly medians (2021+ fit residuals):
- 2021: `+3.76 bpm`
- 2022: `+1.09`
- 2023: `-1.83`
- 2024: `+0.75`
- 2025: `-1.29`
- 2026: `-1.17` (sparse)

**Conclusion.** The corrected exercise-side HR data does **not** show a clear deterioration in walking efficiency after 2021. If anything, 2023-2025 are slightly better than 2021-2022 at similar walking pace. The earlier tracker-based worsening signal was an artifact of using the wrong table.

Command: `python analysis/BW_exercise_live_hr_walks.py`
Artifact: `analysis/BW_walk_live_hr_sessions.csv`, `analysis/BW_walk_live_hr_year_summary.csv`

# Tirzepatide

The set point (AG) and metabolic rate adjustment (K) describe the body's endogenous weight regulation. Tirzepatide intervenes pharmacologically. How does it interact with the system described above?

## F. Pharmacokinetics

One-compartment SC model (FDA: t½=5d, Tmax=24h, ka=3.31/day). 80 weekly injections, dose escalation 2.5→12.5mg.

**Appetite suppression: -49 cal per unit of effective level (F), revised to -74.5 cal per unit blood level (AX).** F's cross-sectional regression confounded drug effect with set point pressure. AX's within-week identification (week fixed effects) isolates the pure drug effect: **-74.5 cal per unit blood level**, 1.5x larger. Tachyphylaxis revised from 32-week per-dose HL (F) to **35-week cumulative HL** (AX). At 12.5mg: blood level cycles 8-17 through the weekly sawtooth, giving -600 to -1270 cal/day raw drug effect at peak, attenuated by cumulative tachyphylaxis. Observed intake: day 0: 1643 cal, day 5: 2222 cal (579 cal/day swing).

**Metabolic suppression: -9 cal per unit of effective level (AQ).** The body normally burns more during weight loss (+90 cal/day above composition-predicted, K). On the drug, direct calorimetry (3 on-drug measurements) shows RMR 206 cal *below* composition-predicted — **the metabolic boost is not just absent but reversed**. At effective level 7: -63 cal/day of metabolic cooperation removed. This is the drug working against itself on the calorie-burning side.

**Net energy budget:** -450 cal intake reduction, +200 cal metabolic cooperation lost, net ~-250 cal/day deficit. Overall reduction: 456 cal/day (18.6%) from pre-tirzepatide year. Weight loss is 93% fat: FM 83→60 lbs, lean 141→143 lbs. The [15mg trial](https://doi.org/10.2337/db23-128-OR) reported -900 cal/day; this dataset's -456 at 12.5mg with tachyphylaxis is consistent.

Command: `python analysis/F_tirzepatide_pk.py`, `python analysis/P2_rmr_model.py` (RMR/calorimetry numbers)
Artifact: `drugs/tirzepatide.csv`, `RMR/rmr.csv`

## G. Transition dynamics

The PK model (F) shows tirzepatide reduces mean intake. The set point (AG/AM) operates through mean daily surplus, which manifests as binge frequency at higher thresholds. Does the drug reduce the mean or break the escalation pattern?

State machine: restriction (<1800), typical (1800-2399), high (2400-2799), binge (2800+).
- High → binge: 19.6% pre-drug → 4.8% on drug
- Binge → binge: 31.6% → 0%
- Post-restriction binge (7 days): 25.9% → 13.0%

The drug breaks escalation and persistence, not just mean level.

Command: `python analysis/G_tirzepatide_dynamics.py`
Artifact: `analysis/G_tirzepatide_transition_summary.csv`, `analysis/G_tirzepatide_rebound_summary.csv`

## BX. Tirzepatide and heart rate

Two separate heart-rate checks:

1. **Walking exercise HR.** Using the exercise live-data charts (not the passive tracker), fit pace-adjusted walking HR on pre-tirz 2021+ walks:

`walk HR = 81.73 + 5.91 * kph`

Applying that pre-tirz fit to all 2021+ walks:
- pre-tirz walk residual median: `+0.32 bpm`
- on-tirz walk residual median: `-1.51 bpm`
- delta: `-1.83 bpm`

So at matched walking pace, on-drug walks are slightly **lower** heart-rate, not higher.

2. **Passive overnight tracker HR.** Using `tracker.heart_rate` tag `21313`, overnight windows with at least 3 hourly points:
- pre-tirz nightly median HR: `83.5`
- on-tirz nightly median HR: `84.0`
- delta: `+0.5 bpm`

**Conclusion.** Tirzepatide does not show a large heart-rate shift here. Exercise-side walking HR is slightly lower on-drug at matched pace (`-1.8 bpm`), while passive overnight HR is essentially unchanged to slightly higher (`+0.5 bpm`). Any strong claim of a clinically meaningful HR increase from this subject's Samsung data would be overstated.

Command: `python analysis/BX_tirz_heart_rate.py`
Artifact: `analysis/BX_tirz_heart_rate_summary.csv`

## AJ. Tirzepatide suppresses the metabolic cooperation with weight loss

Finding K shows falling phases have elevated TDEE (+0.075 on TDEE/RMR ratio) — the body burns more when losing weight, which normally accelerates fat loss by widening the gap between intake and expenditure. Finding F shows 206 cal lower metabolic rate on tirzepatide. But tirzepatide is a falling phase — shouldn't TDEE be *elevated*? Does the drug suppress the body's own assistance?

**Yes.** At matched fat mass (FM 61-71 lbs), tirz falling shows TDEE/RMR = 1.098 vs pre-tirz falling = 1.182 (Δ = -0.084). The drug eliminates most of the falling-phase metabolic boost.

Regression (all days, TDEE/RMR ~ FM + falling + rising + on_tirz + falling×tirz):
- Pre-tirz falling effect: +0.075
- Tirz falling total effect: +0.012
- The falling×tirz interaction is -0.047, wiping out 63% of the falling-phase boost.

Within falling days only, controlling for FM: tirz coefficient = -0.095.

**The drug gives on appetite and takes on metabolism.** F shows -450 cal/day intake reduction. But the metabolic boost that normally accompanies weight loss (+75 on the ratio, roughly +150 cal/day of extra burning that accelerates loss) is suppressed to +12. This means less of the dietary restriction translates into actual fat loss: -450 intake reduction + 150 metabolic cooperation lost = net ~-300 cal/day effective deficit, consistent with the observed -250 cal/day and -1.64 lbs/month fat loss (vs -3.06 pre-tirz when the metabolic boost was intact).

**Direct calorimetry confirms.** Pre-tirz falling calorimetry (n=8): measured RMR 90 cal above composition-predicted (body burning more during loss). On-tirz calorimetry (n=2): measured RMR 193 cal below predicted (body burning less). The sign flip is visible in direct measurement.

**Mechanism.** The drug does not merely reduce appetite while leaving the metabolic response intact. It affects both channels: appetite (G: breaks binge escalation) and metabolism (this finding: eliminates the falling-phase calorie-burning boost). This is consistent with GLP-1 agonists acting on hypothalamic energy homeostasis rather than purely on satiety circuits. The body's normal response to weight loss — burning more, which accelerates the loss while FM is above the set point — is pharmacologically suppressed.

**FM-matched confirmation (AO).** Finding AN showed the metabolic rate adjustment is 4.2x stronger when FM is above SP (where it assists loss by burning more) than below (where it barely adjusts). Tirz suppresses this adjustment. At matched fat mass and matched below-SP status: FM 60-70, pre-tirz TDEE-RMR = +191, on-tirz = +155, Δ = -37 cal. FM 70-84: pre-tirz = +344, on-tirz = +246, Δ = -98 cal. The drug effect is larger at higher fat mass.

Command: `python analysis/AJ_tirz_expenditure_defense.py`, `python analysis/AO_tirz_asymmetric_defense.py`

## AL. Walking partially rescues the suppressed metabolic boost

AJ shows tirzepatide eliminates 63% of the falling-phase metabolic boost (~125 cal/day of extra burning lost). AD shows walk sessions raise RMR independently of composition. AH#5 says exercise and the set point are independent mechanisms. If all three are correct, walking should be additive — partially restoring the metabolic cooperation that the drug suppresses.

**Walking predicts TDEE within the tirz era.** Walks (30d) vs TDEE/RMR: r=+0.21, controlling for FM and drug level (n=529 tirz days). During tirz falling phases (n=357), above-median walks (≥3 sessions/30d, n=211) show TDEE/RMR 1.130 vs below-median (n=146) at 1.082 — a +95 cal/day difference at the same fat mass (70.5 vs 70.8 lbs).

**The drug barely dampens the walk effect.** Regression (TDEE/RMR ~ FM + on_tirz + walks_30d + walks×tirz): walk coefficient pre-tirz = +0.00177 per session, on tirz = +0.00133 (interaction -0.00044, small). The drug preserves ~75% of the walk effect. At mean RMR: +3.5 cal/day per walk-session pre-tirz, +2.6 cal/day on tirz.

**Partial rescue, not full recovery.** The lost metabolic boost is 125 cal/day. At a realistic 15 walks/month: 15 × 2.6 = 39 cal/day recovered (~31% of the loss). Full recovery would require ~47 sessions/month (11/week) — unrealistic. But within the achievable range, walking is the largest single behavioral lever for calorie burning during GLP-1 treatment.

**Pre-tirz comparison.** Walk effect on Kalman TDEE is actually *weaker* pre-tirz after controlling for FM (partial r=0.04 vs 0.19 on tirz). This is because pre-tirz, the body's own metabolic adjustment dominates — walks are a rounding error on top of the body's natural response. On the drug, that natural response is suppressed, making walks the primary remaining lever for calorie burning.

Command: `python analysis/AL_walk_rescue_expenditure.py`

# Nulls and minor effects

## X. Temperature and metabolic state

1,305 temperature readings (Dec 2023 - Mar 2026), timezone-corrected for travel, with a 5am day boundary so late-night readings map to the correct waking period. Each reading is converted to a daily baseline estimate by subtracting the empirical circadian offset for its hours-since-wake.

**Circadian curve.** 0.90°F peak-to-trough, wake-anchored (R² = 0.072, beats sunrise 0.056 and clock 0.045 on downstream signal strength). The curve shows a morning rise, afternoon dip at 5-6h, main peak at 11-12h after wake, and pre-sleep decline:

| Hours after wake | Mean °F | Offset | n |
|---|---|---|---|
| 0-1 (wake) | 97.42 | -0.25 | 98 |
| 2-3 (morning) | 97.72 | +0.06 | 83 |
| 5-6 (afternoon dip) | 97.53 | -0.14 | 72 |
| 11-12 (peak) | 97.99 | +0.32 | 92 |
| 16-17 (pre-sleep) | 97.35 | -0.31 | 28 |

**Injection-day sawtooth.** The strongest temperature signal. Baseline-corrected temperature by day post-injection: day 0 = 97.41°F, day 5 = 97.74°F (0.33°F swing, r = +0.22, partial|calories r = +0.22). The drug modulates body temperature through a non-caloric pathway — the signal strengthens after controlling for same-day intake.

**Rising phase is hot.** Weight gain phases are +0.30°F warmer than stable phases (regression controlling for FM, drug, phase). Falling phase: +0.06°F. Consistent with diet-induced thermogenesis during sustained surplus.

**Intake → temperature: weak.** 14-day trailing intake, baseline-corrected, partial|FM: r = +0.10. Real but small.

**Nulls.** Drug blood level → baseline temperature|FM: r = -0.03 (zero; raw r = +0.45 was season-confounded). Pre-tirz vs on-tirz baseline at matched FM: Δ = -0.07°F (negligible). TDEE residual → baseline|FM: r = +0.03 (temperature does not track the metabolic rate adjustment to the set point at daily resolution).

Command: `python analysis/AS_temperature_baseline.py`, `python analysis/AT_temperature_retest.py`
Artifact: `analysis/AS_temperature_baseline.csv`, `analysis/AS_circadian_curve.csv`

## AE. Sleep and energy balance (null)

2,057 sleep measurements (median 7.8h, std 1.5h, schedule 3am-11am). Sleep duration → next-day calories: r = -0.01. Sleep → next-day steps: r = -0.02. No trailing window of sleep duration predicts intake, TDEE, or steps. The schedule is too consistent (CV 0.19) to detect effects that require cross-sectional variation.

**An extra hour of sleep still does not show an independent effect.** Recasting the question as accumulated sleep debt does not rescue the signal. Raw sleep-debt correlations with metabolic proxies are negative, but they vanish after controlling for bedtime timing. Example: 30-day sleep debt vs TDEE residual is raw r = -0.24, partial r = +0.03 controlling for bedtime, FM, weekend, and year. The apparent "sleep debt lowers RMR" pattern was mostly shared structure with later bedtimes.

**Late bedtime matters a little; sleep duration does not.** Bedtime timing carries a small but persistent signal on the daily metabolic proxies even after controlling for sleep duration. Later bedtimes predict slightly lower TDEE residual / TDEE-to-RMR ratio, with the strongest effect on trailing windows: same-day bedtime partial r = -0.04, 7-day bedtime -0.06, 14-day -0.08, 30-day -0.10. Baseline temperature does not track the same signal (baseline temp vs TDEE residual partial r = +0.01).

**Interaction null.** Combining the two ideas does not help. The late-bed × short-sleep interaction is still null for TDEE residual (partial r = +0.01), TDEE/RMR ratio (+0.01), and same-day calories (-0.02). In this dataset, the independent signal is weak circadian timing, not quantity or debt.

Command: `python analysis/AE_sleep_null.py`
Artifact: `analysis/BS_sleep_timing_cutoff_sweep.csv`

## BY. Samsung exercise weather coverage

Samsung's export includes a per-exercise weather table with `exercise_id`, start-time weather snapshot, latitude, longitude, temperature, humidity, phrase, wind, UV, and provider. This is useful local metadata, but it is **not** hourly historical weather through the workout. It is one start-time snapshot row per exercise.

Coverage in the newest archive is uneven: 589 matched weather rows out of 1,809 exercises overall (32.6%). Coverage is strong in 2018 (87%), 2020 (62%), and 2021 (70%), weak in 2022 (17%) and 2023 (5%), and absent in 2024-2026. By exercise type, walking has the bulk of the rows (558 weather-tagged walks; 34% coverage), with smaller coverage for running (24 rows; 21%) and biking (7 rows; 22%).

Providers are mostly TWC (461 rows) with a smaller AccuWeather block (128 rows counting one capitalization variant). The table appears usable for start-time weather adjustment, but it needs light cleaning: at least one humidity entry is invalid (`-1`), so climate analyses should clamp or drop out-of-range values. For true hourly weather by route and time, you would still need an external historical weather source keyed by exercise timestamp and coordinates.

Command: `python analysis/BY_exercise_weather.py`
Artifacts: `analysis/BY_exercise_weather_joined.csv`, `analysis/BY_exercise_weather_summary.csv`

## AK. Sunlight exposure (null, but walk-session signal survives it)

Sunlight hours (wake-to-sunset, from solar position at subject's location) vary more than sleep (CV 0.20 vs 0.16). Sunlight → TDEE|FM: r = +0.11. No sunlight window predicts intake. AD's walk-session → RMR finding is not a sunlight confound: controlling for sunlight, walks retain r = +0.66 (vs r = +0.73 raw); controlling for walks, sunlight drops to r = +0.36. In LOO-CV: walks alone RMSE = 116, sunlight alone 144, walks + sunlight 121. Sunlight adds nothing beyond walks.

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

The set point and metabolic rate adjustment together explain the 2013 inflection — the moment the decade of regain began.

## AC. The 2013 inflection

Fat mass bottomed at 17 lbs (Oct 2013) and rose every year for a decade. What triggered the regain — metabolic adaptation forcing intake up, or behavioral momentum?

Bottom (2013 Jul-Sep): intake 1864, TDEE 2021, ratio 0.977, binge rate 1%. Regain (2014 Jul-Dec): intake 2359, TDEE 2058, ratio 1.026, binge rate 20%.

TDEE recovered +38 cal from bottom to regain. Intake increased +495 cal. Intake outran TDEE recovery by 457 cal/day. The metabolic adaptation was easing — ratio went from 0.98 to 1.03. But binge frequency exploded: 0% (Jul-Aug 2013) → 13% (Oct) → 20-23% (Nov-Dec) → sustained 12-26% through all of 2014. By 2014-H2, binges were running at 20%.

The trigger was behavioral, not metabolic. The body's metabolic rate was recovering — TDEE/RMR rose from 0.98 to 1.03 — meaning the metabolic channel was returning to normal and no longer helping weight loss. But food noise — the variance, the binge clustering — took hold at FM=17 and never let go for a decade. This is one case study of the [1-in-210 statistic](https://doi.org/10.2105/AJPH.2015.302773): of men who pass BMI 30, only 0.5% reach normal weight within 9 years. The set point model (AG) explains the mechanism: at FM=17, the set point had stalled near 18.6 lbs (AH test #6, floor effect). Any weight above the floor put FM below the set point, triggering the eating-pressure gradient. The metabolic channel (K) had already normalized — it was no longer burning extra — but the eating channel was fully engaged and overwhelmed any remaining metabolic contribution.

Command: `python analysis/AC_inflection_2013.py`
