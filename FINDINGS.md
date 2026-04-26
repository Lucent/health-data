# Findings

Each section has a `Command:` line pointing to the standalone script that reproduces the numbers. Untested hypotheses are in [ROADMAP.md](ROADMAP.md).

# Data quality

## A. Energy balance quality

Does the intake data close against weight and calorimetry?

Cumulative energy balance residual: ±5 lbs over 15 years. Undercount by trajectory phase (composition-aware model, 25 calorimetry × 70 composition measurements): gaining 14.7%, losing 8.6%, stable 12.4%. Approximately uniform.

Command: `python analysis/A_energy_balance_quality.py`

## Z. TDEE by year

Does derived TDEE vary with weight, or does the body defend a band?

| Year | TDEE | Intake | Gap | Fat (lbs) | TDEE/RMR |
|---|---|---|---|---|---|
| 2011 | 2575 | 1772 | +804 | 100→54 | 1.182 |
| 2012 | 2204 | 1917 | +288 | 54→38 | 1.039 |
| 2013 | 2148 | 1847 | +301 | 38→22 | 1.030 |
| 2014 | 2087 | 2105 | -18 | 22→35 | 1.033 |
| 2015 | 2116 | 2130 | -14 | 35→48 | 1.080 |
| 2016 | 2131 | 1954 | +177 | 48→47 | 1.104 |
| 2017 | 2239 | 2140 | +99 | 47→50 | 1.147 |
| 2018 | 1913 | 2005 | -92 | 50→64 | 0.966 |
| 2019 | 1955 | 2018 | -63 | 64→68 | 0.961 |
| 2020 | 2140 | 2120 | +21 | 68→70 | 1.076 |
| 2021 | 2155 | 2136 | +19 | 70→74 | 1.083 |
| 2022 | 2220 | 2224 | -4 | 74→77 | 1.123 |
| 2023 | 2293 | 2267 | +26 | 78→86 | 1.154 |
| 2024 | 2394 | 2274 | +120 | 86→79 | 1.195 |
| 2025 | 2141 | 1943 | +198 | 79→62 | 1.089 |
| 2026 | 1977 | 1844 | +133 | 62→61 | 1.012 |

Start→end of year. Notable: **2018-19 keto era TDEE/RMR drops to 0.96** — keto is metabolically suppressed, not neutral. 2017 sits at 1.15 (higher than the surrounding years). 2024 drug-era at 1.20 is the all-time high.

Command: `python analysis/Z_tdee_by_year.py`
Artifact: `analysis/P4_kalman_daily.csv`

# Metabolic rate adjustment

The year-by-year table shows TDEE falling from 2661 to 2044 as fat mass dropped from 106 to 25 lbs — partly from losing metabolically active tissue, partly from metabolic adaptation. TDEE barely recovered during a decade of regain. This raises the question: at the same body composition, does TDEE differ depending on direction?

## K. TDEE hysteresis

At the same fat mass, does TDEE differ depending on whether weight was reached from above (falling) or below (rising)? [Leibel et al. (1995)](https://doi.org/10.1056/NEJM199503093321001) demonstrated this in a metabolic ward (n=18, ±10% body weight). This dataset replicates the finding free-living with 1,482 FM-matched pairs.

Matched fat-mass bands (pre-tirzepatide):
- FM 25-45 lbs: rising 2027 vs falling 2152 (+124)
- FM 45-65 lbs: rising 2028 vs falling 2289 (+261)
- FM 65-85 lbs: rising 2157 vs falling 2465 (+308)

Regression: TDEE = 1916 + 4.25×FM + 169×falling − 85×rising. Falling phase adds ~170 cal; rising phase compresses by ~85 cal — the metabolic channel is asymmetric across direction.

Robustness: 1,482 FM-matched pairs (within 2 lbs), mean difference +179 cal (falling − rising). Pair lag-1 autocorrelation 0.83, effective n = 142. Bootstrap 95% CI: [+104, +267]. p < 10⁻²². This is the most statistically robust finding in the dataset.

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

K and B show the body adjusts metabolic rate by direction. Is there a single defended fat mass driving both responses? [Speakman (2011)](https://doi.org/10.1242/dmm.008698) proposes dual intervention points with a "zone of indifference." The data finds no such zone — the response is graded.

**Methodological frame.** Any model that defines the set point as a smoothing of recent fat mass and is then evaluated against a smoothing of recent surplus shares input bandwidth with what it tries to predict (BR makes this explicit, below). To avoid that ridge, the set-point form is identified primarily from **published GLP-1 withdrawal trials** — where this subject's intake autocorrelation cannot bias the answer — and only afterwards confirmed against this subject's binge rate, intake-free FM, and weigh-in-only outcomes. The continuous-surplus regression (AM/AP) is treated as consistent supporting evidence, not a primary identification.

**Receptor-dynamics reading.** The model that survives this protocol has the shape of a two-receptor adaptation system: the set point is the equilibrium of slow leptin/melanocortin receptor populations that up-regulate when chronic exposure is sparse and down-regulate when it is dense. The 45-50 d half-life is the receptor-turnover and AgRP/POMC neuropeptide-equilibration timescale. Adaptation requires the receptor to be repeatedly stimulated near a new equilibrium — explicitly *not* during rapid drug-driven loss, where the receptor is bypassed by exogenous GLP-1/GIP agonism. This frames every section below.

## What this is, by elimination

The SmoothLatch + receptor-dynamics reading is what's left after testing every other candidate set-point form against this dataset. The remainder of this section establishes each elimination — but the chain is short enough to summarize first.

| Candidate | Source | Falsifying finding | Verdict |
|---|---|---|---|
| Pure CICO, no homeostat | implicit in textbook | Energy balance closes (A: ±5 lb/15 yr) — but TDEE depends on direction (K: +179 cal falling vs rising at matched FM, p < 10⁻²²); ratio drops to 0.97 in restriction, recovers to 1.14 when restriction eases (Z). The body actively defends. | Insufficient as a model |
| Gravitostat (foot-pounds drive intake) | [Jansson et al. 2018](https://doi.org/10.1073/pnas.1800033115) | Foot-pounds → next-day intake: r = +0.050 partial. Wrong direction (more weight-bearing → *more* eating). Steps → weight change: r = +0.010 at every timescale (N, U). | Falsified |
| Dual intervention points with "zone of indifference" | [Speakman 2011](https://doi.org/10.1242/dmm.008698) | Sigmoid binge response across the gap — continuous gradient, no flat zone. Surplus is graded across all gap sizes (AM). | Falsified |
| Set point tracks scale weight | gravitostat-adjacent | At every half-life, fat mass beats total weight beats scale: r = -0.60 / -0.53 / -0.52. The body weighs *fat*, not total mass. Consistent with [Kennedy 1953](https://doi.org/10.1098/rspb.1953.0009) and the leptin system. | Falsified |
| Fixed set point at any specific weight | classical "set point theory" | Best fixed defended weight: r = +0.20 (and wrong sign). Best moving 50d EMA: r = -0.60. Difference 4×, sign-reversed. | Falsified |
| Fixed weight-suppression proxy | [Lowe et al. 2007](https://doi.org/10.1002/eat.20405) (bulimia binges) | Lowe's directional finding (suppression below historical max → binges) replicates here, but their fixed-historical-max proxy is dominated by a moving 50d EMA at every half-life tested (AG). Right intuition, wrong reference. | Form falsified, signal real |
| Symmetric EMA (any half-life) of FM | extension of classical lipostat | Fits this subject's slow-gain decade at r = -0.94, but undershoots SURMOUNT-4 regain at +2.8% (published +14.0%) and STEP-1 at ~+2% (~+10%). Off by 5-10×. | Falsified by trials |
| Asymmetric ratchet (downward HL faster than upward) | data-driven hypothesis (AH original) | Improved subject fit r −0.60 → −0.69 on full dataset. Does not replicate on 2014+ natural-dynamics data (best Δr = 0.0004). The asymmetry was the EMA mechanically chasing the willpower-era descent. | Retracted |
| One-stock depletable willpower reservoir | data-driven hypothesis (BK) | Best coarse fit produces 25.6% "impossible" days where required willpower exceeds any plausible reservoir, concentrated in 2011-12. | Falsified |
| Yo-yo / variance as metabolic damage | folk wisdom; cf. [Biggest Loser reanalysis](https://doi.org/10.1002/oby.23308) | At matched mean surplus, higher 14d intake std predicts *less* fat gain (partial r = -0.20), driven by transient acute-deficit days raising TDEE. Sign-reversed. | Falsified |
| Constrained total energy expenditure | [Pontzer et al. 2016](https://doi.org/10.1016/j.cub.2015.12.046) | Walk *sessions* (not steps) raise *resting* metabolic rate at +14 cal/session (AD: CV RMSE 116 vs noise floor 170, n=23 calorimetry). RMR rises with intentional exercise; not just TEE redistribution. | Violated on RMR |
| "Metabolic adaptation" as the GLP-1 plateau | clinical/lay consensus | Names a phenomenon, not a mechanism. The plateau is fully accounted for by the three-timescale race between drug suppression D(t), gap pressure 55 × (SP − FM), and SmoothLatch recalibration (see Tirzepatide § Asymptote). | Reframed |

What's left after these eliminations: a fat-mass-tracking, slow-moving, stability-gated defended state with a 45-50 d adaptation half-life that adapts only when FM holds within a narrow band — the SmoothLatch — interacting with a fast (~13 d) metabolic channel that's asymmetric above vs below SP, and with exogenous GLP-1/GIP agonism that bypasses the endogenous receptor and therefore prevents the latch from firing during drug-driven loss. The next subsections establish each piece.

Each row of the falsification table cites the finding (and script) that produces the falsifying number — gravitostat (N), zone of indifference (AM/AP), scale-weight tracking (AG), fixed SP (AG), Lowe proxy (AG), symmetric EMA (AY/AZ/BD), asymmetric ratchet (AH), one-stock willpower (BK), variance-as-damage (AF), constrained TEE (AD), GLP-1 plateau ("Asymptote of GLP-1 weight loss" under Tirzepatide). Reproducer scripts are listed under each named finding's section.

Command: `python analysis/N_dietary_predictors.py`, `python analysis/AG_binge_set_point.py`, `python analysis/AH_set_point_properties.py`, `python analysis/AM_lipostat_sensitivity.py`, `python analysis/AP_overshoot_shape.py`, `python analysis/AF_intake_variance.py`, `python analysis/AD_tdee_formula_sweep.py`, `python analysis/BK_control_stock_model.py`, `python analysis/BD_smooth_latch.py`, `python analysis/AY_sp_from_regain.py`, `python analysis/AZ_sp_model_search.py`

## SmoothLatch: the set point only adapts when fat mass is stable

Identified by inverting published GLP-1 trial regain (AY) and confirmed by exhaustive search across 268 configurations (AZ). **The set point adapts toward fat mass only when FM has been within ±3 lbs of a reference for 14+ consecutive days, at rate 0.01/day.** When FM changes rapidly, the latch is open and the SP freezes. Pressure on intake: ~55 cal/day per lb of gap, symmetric (binned data across -10 to +20 lbs).

**Subject fit decisive vs EMA.** No EMA half-life beats the SmoothLatch on this subject's 2014+ pre-tirz data (both r ≈ -0.94 on 180d surplus). On the subject's slow-gain decade the two are observationally equivalent — FM changed at ~6 lbs/yr, always within latch tolerance, so SP adapted continuously under either rule.

**Trial-only refit reaches the same family (BG).** Searching across SURMOUNT-4, STEP-1 extension, and STEP-4 with no subject data, the lowest worst-case error in the simple-model family lands on a HoldMean latch at tol≈2 lbs, hold≈42-56 d, rate≈0.0025/day — same family, somewhat tighter parameters. The latch family wins in trial space *without ever touching this subject's intake*.

**Speed-of-loss prediction (preserved across parameterisations).** At the same total fat loss, faster loss produces more regain because the SP had less time to latch:

| Loss rate | Weeks | Total loss | Gap at stop | Predicted regain |
|---|---|---|---|---|
| 0.5 lb/week | 100 | ~12% | 2 lbs | +0.7% |
| 1.0 lb/week | 53 | ~13% | 5 lbs | +2.3% |
| 1.5 lb/week | 35 | ~13% | 11 lbs | +4.7% |
| 2.0 lb/week | 26 | ~13% | 17 lbs | +7.0% |

This survives every (tol, hold, rate) combination tested in CF (r = +0.6 to +0.95 between trial loss rate and predicted regain asymptote). Differentiates SmoothLatch from any EMA: EMA predicts regain depends only on time; SmoothLatch predicts it depends on speed. Consistent with the clinical observation that surgical and very-low-calorie patients regain more than gradual losers. **The published meta-regression of 6 GLP-1RA withdrawal trials explicitly flags speed-of-loss as an untested moderator (CF, [Wilding 2025](https://www.medrxiv.org/content/10.1101/2025.06.09.25328726v1.full)).**

**Magnitude limitation: a single latch parameter set does not fit both subject and trial-cohort regain magnitudes.** The published meta-regression has R(52) = 60% mean across 6 RCTs; carrying our subject-fit latch + pressure forward predicts ~14% mean (CF). The most plausible structural fix is that **pressure constant scales with central-pathway responsiveness `m`** (CE): subject `m = 0.52` ↔ pressure ≈ 55 cal/lb; non-diabetic trial cohort `m ≈ 1.0` ↔ pressure ≈ 106 cal/lb. Under this scaling, the same loose latch fits both populations because the pressure × gap term is population-specific. See CF for details.

**Receptor-dynamics reading.** The 14-day stability + 3-lb tolerance is what receptor recalibration requires: leptin/melanocortin receptors adjust density and downstream signalling only after the input has held within a narrow band long enough for receptor turnover to track it. During drug-driven loss the satiety circuit is being driven by exogenous GLP-1/GIP, so the endogenous receptor isn't getting the new chronic-leptin baseline to recalibrate against. The same `m` factor that scales drug response and lipostat pressure is the central-pathway "gain" of this entire feedback loop.

Command: `python analysis/AY_sp_from_regain.py`, `python analysis/AZ_sp_model_search.py`, `python analysis/BD_smooth_latch.py`, `python analysis/BG_trial_regain_search.py`, `python analysis/CF_regain_vs_meta.py`

## Independent confirmations from this subject's data

Three measurement combinations find the same 45-50 d timescale without depending on rolling-surplus autocorrelation. They converge on the same answer the trial fit gives.

**Binge rate (AG; binary outcome, no surplus involved).** Distance below the moving SP predicts 90-day binge rate at r = -0.62, partial r = -0.72 controlling for absolute FM. Half-life sweep peaks at HL ≈ 50 d.

*Tracks fat mass, not weight.* At every HL: r = -0.60 (fat mass) > -0.53 (total weight) > -0.52 (scale weight). The defended quantity is fat, consistent with [Kennedy (1953)](https://doi.org/10.1098/rspb.1953.0009) and the leptin system; gravitostat-style scale-weight tracking is dominated.

*Moves, not fixed.* Best moving 50d-EMA SP: r = -0.60. Best *fixed* defended weight (any single value): r = +0.20 — wrong sign and 4× weaker. The literature has no comparable measurement: clinical sources cite vague "1-6 years" for set-point convergence with no quantification. [Lowe et al. (2007)](https://doi.org/10.1002/eat.20405) showed weight suppression predicts binge frequency in bulimia using a fixed historical-max proxy; the directional finding replicates here, but the fixed proxy is dominated by the moving EMA at every HL.

*Sigmoid response, no threshold or zone.* Continuous gradient across the gap:

| Distance from SP | Binge rate |
|---|---|
| 5-7.5 lbs below | 14.8% |
| 2.5-5 lbs below | 11.4% |
| 0-2.5 lbs below | 6.0% |
| 0-2.5 lbs above | 3.3% |
| 2.5-5 lbs above | 2.8% |
| 5-10 lbs above | 1.3% |

**Intake-free fat mass (P3; Kalman uses logged intake, P3 does not).** P3 derives daily FM by interpolating between weigh-ins only — no intake feeds the FM estimate. Using P3's intake-free FM on 2014+ data, the set point still predicts mean surplus at HL = 50 d (r = -0.73). The half-life is the same; only the correlation strength drops because P3's FM is noisier between weigh-ins. So the 45-50 d timescale is a property of fat mass dynamics, not a Kalman filter artifact.

**Non-circular weigh-in-only test (AU).** Building a trailing EMA "defended weight" from corrected morning weigh-ins only (no Kalman, no RMR, no intake-fed interpolation), the half-life sweep peaks at HL ≈ 35-40 d on next-30-day intake (r = -0.41, partial r = -0.40 controlling for absolute weight). The SP-distance signal beats both fixed defended weight (r = -0.13) and recent weight change. Same direction, same band of timescales, with no intake-side circularity at all.

Three independent legs — binary outcome, intake-free FM, weigh-in-only — converge on 35-50 d. Combined with trial regain (BD/AZ/BG) above, this is what supports the timescale, not the autocorrelated surplus regression.

Command: `python analysis/AG_binge_set_point.py`, `python analysis/AU_noncircular_set_point.py`, `python analysis/H_phase_binge_landscape.py`
Artifact: `analysis/H_phase_binge_summary.csv`, `analysis/H_phase_binge_distance_bins.csv`, `analysis/H_phase_binge_model_auc.csv`

## BR. Why the surplus regression is consistent but doesn't pin the half-life

A 2-D sweep over symmetric SP half-life and surplus lookback window produces a broad diagonal ridge, not a sharp optimum. On 2014+ pre-tirz data: HL=20d × 45d window gives r=-0.939; HL=40d × 75d gives r=-0.929; HL=45d × 90d gives r=-0.926; HL=75d × 120d gives r=-0.912. **The mean-surplus fit alone does not identify a structural 45 d half-life.** It identifies *some* slow appetite-pressure timescale.

This is mechanical, not a defect of the data. SP is a low-pass of FM, and FM is the cumulative integral of surplus, so (SP − FM) is mathematically a high-pass of cumulative surplus while the rolling-surplus window is a band-pass of surplus. Their correlation reflects shared input bandwidth. The r = -0.92 says the slow-timescale signals line up; it does not say the timescale is 45 d.

What pins the timescale is the *intersection* of this ridge with the trial regain anchor (BD/AZ/BG, no autocorrelation), the binge-rate fit (AG, binary), the intake-free FM fit (P3, no intake), and the weigh-in-only test (AU, no Kalman). All four hit the same 35-50 d band. The surplus regression is consistent with that intersection but should not be cited as the primary evidence for it.

**AM/AP findings that survive this demotion** (the regression isn't load-bearing for any of them):
- Per-lb pressure ~55 cal/day per lb (binned, robust across -10 to +20 lb gap; the regression-derived -27 cal/lb is autocorrelation-diluted).
- The set point tilts the entire daily distribution: it controls how often you eat above maintenance (r = -0.91) and how deeply you restrict on deficit days (r = -0.82), but not individual binge size (r = -0.40 conditional magnitude, +0.01 R² added when included).
- Restriction-run prediction r = -0.48 to -0.72 (above-SP runs continue losing -1.2 lb/30d; below-SP runs rebound +0.3 lb/30d).

Command: `python analysis/BR_sp_surplus_identifiability.py`, `python analysis/AM_lipostat_sensitivity.py`, `python analysis/AP_overshoot_shape.py`, `python analysis/AN_ratchet_profile.py`
Artifact: `analysis/BR_sp_surplus_identifiability.csv`

## AH. Mechanics: pressure constant, ratchet retraction

**Pressure: ~55 cal/day per lb of gap, symmetric.** From binned data across -10 to +20 lbs, remarkably constant; five pounds below = ~275 cal/day of upward pressure. Asymmetric models show zero improvement on 2014+ data (best Δr = 0.0004 vs symmetric).

**Apparent ratchet: retracted.** An asymmetric model (HL_down=20d, HL_up=100d) improved binge-rate r from -0.60 to -0.69 on the full dataset. It does not replicate on 2014+ natural-dynamics data. The asymmetry was driven by 2011-13: the EMA mechanically chased a willpower-driven descent that was not actually responding to the SP. On data where intake responds naturally, the SP adapts at the same rate in both directions.

**Binge size is constant; the SP tilts the daily mean.** Binge surplus ~1427 cal regardless of distance (r = -0.02). Non-binge days: r = -0.32 between distance and surplus (+80 cal at 5+ below SP, -481 cal at 5+ above). The SP modulates background drift on every day, not just discrete events.

**Exercise independence.** Walking raises RMR (AD) but does not change the SP half-life (40 d for both high-walk and low-walk periods).

**Floor.** SP minimum: 14.8 lb FM on 2014-08-01. Post-floor: FM 15 → 38 over six months, binge rate 17.8% — the eruption is the receptor system coming back online above the leptin floor (see AC + BK).

Command: `python analysis/AH_set_point_properties.py`

## AI. The metabolic channel runs on a ~13-day timescale

The eating channel adapts with a 50 d half-life (above). **The metabolic channel is faster: optimal HL = 13 days** (partial r = +0.47 controlling for absolute FM, p < 0.001 at n_eff=54). Correlation is monotonically stronger at shorter half-lives, then plateaus. K's hysteresis (+179 cal falling vs rising, p < 10⁻²²) confirms the same signal through FM-matched pairs.

| HL (days) | r (metabolic, partial) |
|---|---|
| 10 | +0.47 |
| 15 | +0.47 |
| 30 | +0.44 |
| 50 | +0.39 |
| 100 | +0.31 |

**Not a Kalman artifact.** Lagging SP distance by 7-60 days at HL=13: partial r decays 0.47 → 0.43 → 0.39 → 0.31 → 0.22. A purely mechanical artifact would vanish by lag 7-14; this persists to 60 days.

**Metabolic asymmetry.** The metabolic adjustment is **4.2× stronger when FM is above SP** (above-SP partial r = +0.52) than below (+0.12). The body actively helps fat loss by burning more when above the set point, but barely adjusts when below.

**Different timescales** — metabolic rate is a sprint (~13 days), eating pressure is a marathon (~50 days). This is the engine of the 2013 inflection: the eating channel ramped up faster than the metabolic channel could rebound, locking in the regain. It reframes ["persistent metabolic adaptation" in Biggest Loser contestants](https://doi.org/10.1002/oby.21538): if the metabolic channel adapts on a two-week timescale, persistent RMR suppression after 6 years means the subjects remained chronically below a set point that itself had moved.

Command: `python analysis/AI_expenditure_arm_timescale.py`
Artifact: `analysis/AI_expenditure_arm_sweep.csv`

## AC + BK. The 2013 inflection and the leptin floor

The 2011-13 willpower era is the cleanest test of the set-point pressure form. Under any EMA-style model, FM falling 100 → 22 lb on sustained 1200-cal restriction should produce an enormous SP gap and explosive eating pressure. It does not — binges stay near 0% through the deepest restriction (FM 13-22 lb in mid-2014), then erupt as FM begins to rise. Under a one-stock control-stock attempt to model that period (BK), 25.6% of days have "impossible" required-control values exceeding any plausible reservoir. Neither smooth-EMA nor one-stock willpower fits the willpower era.

The receptor-dynamics reading resolves both:

**Below the leptin floor (FM ≲ 15-18 lb), the eating-pressure signal is gated, not maximal.** Leptin secretion is roughly proportional to FM and saturates low at very low fat mass. Below the floor, hypothalamic input from the leptin axis is near-zero regardless of receptor density — the system enters conservation mode (TDEE/RMR collapses to 0.96-0.98) but also stops generating large eating-pressure gradients, because the gradient operates on receptor-occupancy *differences*, not absolute hormone level. The willpower era didn't have hundreds of cal/day of biological pressure for willpower to fight; it had a floored signal under a starvation-mode regulatory regime.

**The 2013-14 inflection is leptin coming back online.** FM trajectory bottoms at 13 lb in June 2014. Binge rate progression: 0% through Jul-Aug 2013 → 13% (Oct) → 20-23% (Nov-Dec) → sustained 12-26% through 2014. AH#6 (floor effect) lands at SP minimum 14.8 lb FM on 2014-08-01; FM rises 15 → 38 over the next six months, binge rate 17.8%. **The bottom and the eruption co-occur: as soon as FM crossed back above the leptin floor, the receptor — calibrated by years of high-FM exposure — read a large deficit and binges began.**

The TDEE story is what one would expect once the metabolic and eating channels are read together: **TDEE *fell* by 136 cal between bottom (1864 intake / 2019 TDEE / ratio 0.977) and regain phase (2359 intake / 1883 TDEE / ratio 0.938).** Intake outran TDEE by 631 cal/day. The metabolic channel (~13 d HL) did not get a chance to rebound because the eating channel ramped up immediately as the floor cleared, keeping fat mass climbing and the body in the "below SP" regime that never engages the metabolic boost.

This is one case study of the [1-in-210 statistic](https://doi.org/10.2105/AJPH.2015.302773): of men who pass BMI 30, only 0.5% reach normal weight within 9 years. The mechanism is now concrete: bottom-out at the leptin floor, conservation mode masks the brewing defence; on any FM increase above the floor the un-recalibrated SP at high-FM-baseline reads a large gap and the eating channel fires. Binges aren't a willpower failure — they're the receptor system finally generating the signal it had been holding back during starvation mode.

**BK reframed.** The control-stock model still fits 2014+ data well — required control collapses to 23-28 cal/day in 2014-15 once FM is above the floor. The "26% impossible days" failure of the one-stock model is concentrated in 2011-12, exactly the period the floor explanation removes from the smooth-EMA pressure picture. BK survives as the cleanest formal statement of the willpower-era anomaly; it does not survive as evidence that biological pressure was 528 cal/day in 2011-12. The cleaner reading is that the floor capped the biological signal during the early restriction, and willpower carried a real but more moderate load.

This is the prediction the receptor-dynamics frame was always going to make: the system's pressure is a function of *receptor occupancy gradients* against a slow-adapting baseline, not a linear function of FM gap regardless of absolute level.

**Residual support from the energy-balance closure.** Year-by-year EB residuals (cumulative implied weight change from `(intake − TDEE)/3500` minus observed FM change) are tightly closed in modern years (typically ±0.5 lb), but **2013 alone holds a -3.45 lb residual** — meaning the filter over-predicts weight loss in 2013 by 3.5 lb beyond what the observed weight trajectory shows. Possible attributions: TDEE was even lower than the filter's free-drift estimate during deep restriction, intake was under-logged more than usual, or the Forbes partition mis-attributed the trajectory. The conservation-mode reading of the floor predicts the first: TDEE drops further than any energy-balance-anchored filter can naturally find, because the model has no signal pulling it that low. The residual sits exactly where the floor mechanism predicts the body is doing something the filter can't see.

Command: `python analysis/AC_inflection_2013.py`, `python analysis/BK_control_stock_model.py`
Artifact: `analysis/BK_control_stock_daily.csv`

# Intake characterization

The set point sets the slow drift in mean daily surplus. But what shapes day-to-day intake on top of that? Several dietary hypotheses predict intake modulation that should be visible in this dataset.

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

At the point estimate, steady-state lean mass above baseline: 1×/week +5.2 lbs, 2×/week +10.3 lbs. **The amplitude per session (~41 g) replicates at daily resolution (CB: lag-1 +39 g). The steady-state extrapolation does not** — composition scans through 2017-2026 hold LM at 138-150 lb without a clean +5 lb shift tied to workout frequency. The 275-day half-life inferred from scan correlations was likely too long; the per-session response is real, the long chronic stack appears to be a fitting artifact. See CB.

Command: `python analysis/AA_lean_mass_training.py`
Artifact: `analysis/AA_lean_mass_training.csv`

## CB. Workouts: per-session lean response and lean preservation during cuts

Using the Kalman's daily lean-mass state (which AA could not — AA had only 70 scans to fit), two questions become well-posed: how much lean does each session add, and how much do workouts protect lean during deficit runs.

**Per-session lean response peaks at lag-1 day, integrating to ~46 g/session.** Regression of daily ΔLM on workout indicator with deficit and FM controls, by lag:

| Lag (days after session) | Coefficient (g/session) |
|---|---|
| 0 (same day) | **−4** |
| **1** | **+39** |
| 3 | +11 |
| 7 | −2 |
| 14 | −1 |

Same-day workouts produce a small LM *decrease* (sweat / glycogen depletion that the P1 correction doesn't fully absorb). The +39 g rebound at lag-1 is the muscle protein synthesis response. Integrated over the week: ~46 g per session. AA's +41 g/session point estimate (averaged over 60-day windows) was right in magnitude; daily resolution localizes the response to the day after.

**Workouts substantially preserve lean during cuts.** 37 cuts identified (≥14 days, rolling 7d intake < TDEE − 200 cal/day). For each cut, observed lean fraction of weight loss vs base-Forbes-predicted lean fraction at that FM:

| Workout density during cut | n runs | Lean % of dW | Forbes-predicted | Δ from Forbes |
|---|---|---|---|---|
| **None** (0 sessions) | 23 | +17% | +13% | **+4%** (more lean lost than Forbes) |
| Low (<1/wk) | 4 | 0% | +6% | -6% |
| **High (≥1/wk)** | 10 | -5% | +8% | **-13%** (less lean lost than Forbes) |

Regression `lean_frac_minus_forbes ~ workouts_per_wk + dW + FM`: **workouts/wk coefficient = -13.7 percentage points** [bootstrap 95% CI: -17.9, -6.6]. Each workout per week reduces the lean fraction of weight loss by ~14 percentage points beyond what plain Forbes predicts.

**Circularity note.** The Kalman uses workout density to modulate Forbes (α = 0.084 per workout/wk; see Kalman appendix). That modulation accounts for ~8 of the ~14 pp; the remaining ~6 pp is data-driven residual lean protection beyond the modulation. The lag-1 +39 g per session result above is unaffected — it's measured directly on dLM and does not depend on the Forbes partition.

Concrete: at a typical post-2014 cut (FM 60-80 lb, Forbes-lean ~6%), zero workouts → lean is 10-15% of weight lost; ≥2 workouts/week → lean is roughly zero or slightly negative (LM grows during the cut). Most striking individual cuts:
- 2018-01-07 (19 d, 7 workouts during keto): -1.3 lb total, **+0.3 lb LM**.
- 2024-10-16 (37 d, 9 workouts, tirz era): -4.0 lb total, **+0.1 lb LM**.
- 2025-08-24 (54 d, 11 workouts, tirz era): -4.9 lb total, only -0.3 lb LM.

The pre-2014 willpower-era no-workout cuts *exceed* Forbes by 1.6-6.2 percentage points — consistent with starvation-mode lean catabolism that Forbes (free-living-calibrated) underestimates. Workouts cancel that excess and then some.

**The per-session spike does not accumulate the way AA's steady-state extrapolation suggested.** AA fit a 41 g/session response with a 275-day half-life from scan-only data, implying steady-state +5-10 lb LM above baseline at 1-2×/week training. The daily LM trajectory does not show this accumulation — composition scans (138-150 lb LM across 2017-2026) do not contain a clean +5 lb baseline shift correlated with workout frequency, so the Kalman cannot manufacture one between scans. The +39 g lag-1 spike is real and localized; the slow chronic stack predicted by AA is likely overstated. The mechanism through which workouts protect lean during cuts is therefore not chronic accumulation but **trajectory interruption**: each session's lag-1 spike interrupts the deficit-driven lean slide; ≥1 spike per week keeps the running trajectory near flat.

**Caveat.** Workout-rich cuts cluster post-2014; the regression controls for FM and dW but the era confound can't be fully broken. The lag-1 +39 g per-session result is more directly identified than the run-level result.

Command: `python analysis/CB_workouts_lean.py`
Artifact: `analysis/CB_workouts_lean_runs.csv`

## CC. Trajectory interruption confirmed: workout protection lasts ~1 day

CB suggested workouts protect lean during cuts via discrete lag-1 spikes rather than chronic stock building. CC tests this directly: bin daily ΔLM during cuts (rolling-7d intake < TDEE − 200) by `days_since_last_workout`. If workouts build a stock that decays slowly, ΔLM should rise gradually with days since the last session. If workouts work via discrete lag-1 spikes, only the day-1 bin should be positive, and all later bins should sit at the baseline cut rate.

| days since last workout | n | mean ΔLM (g/day) | mean deficit (cal) |
|---|---|---|---|
| 0 (workout day) | 128 | -10 | -323 |
| **1** | **131** | **+31** | **-353** |
| 2 | 85 | -11 | -328 |
| 3-4 | 114 | -10 | -259 |
| 5-6 | 41 | -10 | -381 |
| 7-13 | 29 | -12 | -253 |
| 14-29 | 50 | -9 | -456 |
| 30-89 | 80 | -10 | -364 |
| 90-365 | 91 | -7 | -318 |

The lag-1 bin is the only positive entry. Every other bin sits at -10 g/day (the deficit-driven baseline cut rate), regardless of whether the last workout was 2 days ago or 200. **Workout protection is a discrete +30-40 g pulse on the day after each session, then ΔLM returns to baseline.** No gradual decay; no chronic stock.

Recent (≤5 d since last workout) vs distant (≥14 d) at matched deficit:

| Deficit bin | Recent ΔLM | Distant ΔLM | Δ (g/day) |
|---|---|---|---|
| -1500 to -700 cal | +9 g | -10 g | **+19** |
| -700 to -400 | +0 g | -9 g | +9 |
| -400 to -200 | +1 g | -7 g | +8 |

The per-session protection scales with deficit depth — deeper deficits benefit more from a recent session. Consistent with a fixed-amplitude spike interrupting deeper catabolism.

**Implication for protocol design.** Frequency dominates volume. A workout's protection vanishes in ~24 hours; weeks-long gaps re-expose lean to deficit catabolism immediately. Pre-cut "loading" with extra training accomplishes nothing because no stock builds. The relevant quantity is sessions-per-week, not total-sessions-this-cycle.

Command: `python analysis/CC_workout_trajectory.py`
Artifact: `analysis/CC_workout_trajectory_cuts.csv`

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

## AL. Walks do not rescue tirz-suppressed expenditure (retracted)

Tirz-era walk coefficient on TDEE/RMR is **−0.00044, indistinguishable from zero**. An earlier version of this analysis reported a +95 cal/day rescue at high walk frequency; that signal does not replicate and is retracted.

**Mechanism preserved.** AD (calorimetry, RMR-side, n=23 measurements) is unaffected — walks raise resting metabolic rate at +14 cal/session/30d. AJ's "tirzepatide suppresses the falling-phase metabolic boost" finding is unaffected — that comparison uses direct calorimetry, not Kalman state. **What's retracted is the specific claim that walks partially compensate during tirz.**

Command: `python analysis/AL_walk_rescue_expenditure.py`

## AQ. Drug coverage on each channel

Decomposing the drug's effect into appetite and metabolic components, both as functions of effective blood level.

**Appetite arm: constant offset.** Regression (90 d mean surplus ~ SP distance + FM + effective level): SP distance −77 cal/lb [CI excludes zero], effective level +5.4 cal/unit [CI includes zero]. Adding a distance × drug interaction improves R² by < 0.001. The drug does not attenuate the per-lb pressure from the set point — it adds a constant offset of ~49 cal/day. Drug equivalence: 1 unit ≈ 0.1 lb of SP offset.

**Metabolic arm: the drug suppresses the body's calorie-burning boost.** Regression (TDEE residual ~ SP distance + FM + effective level): SP distance +19 cal/lb (body burns more above SP), effective level −10.2 cal/unit (drug reduces this burning). Pre-tirz gets +18.4 cal/lb of metabolic boost; each unit of drug reduces this by 8.7 cal/lb. At eff = 7, the per-lb boost drops from +18 to −43 cal/lb — the drug *reverses* metabolic cooperation into metabolic *resistance*.

| Eff level | Appetite (less eating) | Metabolic (less burning) | Combined |
|---|---|---|---|
| 2 | +11 cal/day | −20 cal/day | −9 cal/day |
| 4 | +21 | −41 | −20 |
| 6 | +32 | −61 | −29 |
| 8 | +43 | −81 | −38 |
| 10 | +54 | −102 | −48 |
| 12 | +64 | −122 | −58 |

Net is still pro-loss because the appetite effect on total energy balance is larger than the lost metabolic boost.

**Tachyphylaxis erosion.** Both effects degrade with cumulative exposure (timescale identified cleanly in AX): metabolic suppression from −102 cal/day (week 0) to −51 (week 32) to −25 (week 64). Appetite reduction degrades proportionally.

**Injection cycle.** The weekly sawtooth is visible in daily surplus: injection day (day 0) −525 cal, rising to +51 cal by day 5. A 576 cal/day swing. Effective drug level drops from 9.6 (day 0) to 5.2 (day 5) through the sawtooth.

Command: `python analysis/AQ_tirz_set_point_coverage.py`

## AX. Clean within-week drug identification

The cross-sectional regressions in F and AQ are confounded with the moving set point: as SP converges toward FM, eating pressure falls, and the cross-sectional fit attributes some of that to the drug. AX isolates the drug effect by exploiting **within-week injection-cycle variation** — week fixed effects absorb the SP, the tachyphylaxis trend, and every slow-moving confound. What's left is the pure drug coefficient.

**Drug appetite effect: −74.5 cal per unit blood level.** Within-week identification. This is 1.5× F's −49 (cross-sectional). The naive cross-sectional estimate on the same data is −37, confirming that cross-week confounds dilute the signal by half.

**Tachyphylaxis: cumulative exposure beats per-dose reset.** Cumulative weeks on drug (r = −0.57) fits much better than per-dose-reset (r = −0.22). Per-dose-reset is not portable — this subject's 20 weeks at 5mg creates deep tachyphylaxis that a 4-week trial escalation never develops. The calendar-cumulative form used here gives a 35-week half-life summary; BZ supersedes this with an AUC-cumulative form (occupancy-driven rather than time-driven), which is the form used in trial simulations.

**Metabolic arm zero within-week.** TDEE residual does not respond to blood level variation within the injection cycle (+0.1 cal/unit, r = +0.04). The metabolic effect operates on longer timescales than the 7-day cycle can capture — consistent with AI's 9 d HL.

**Set-point pressure is held fixed at the SmoothLatch values during the drug-era fit**, so no part of the drug coefficient absorbs SP adaptation. At drug start: SP = 83.5, FM = 84.2, gap = −0.7 lbs (essentially at equilibrium). At end: SP = 62.5, FM = 59.9, gap = +2.6 lbs (FM below SP — the drug pushed past the defended weight, but only by ~2 lbs because the latch chases stable holds).

**Forward validation on this subject: RMSE = 364 cal/day, r = +0.41.** Residuals by injection day range from −215 (day 0, model undershoots peak suppression) to +201 (day 5, model undershoots rebound) — the actual sawtooth is ~100 cal wider than the linear model predicts on both ends, motivating BZ below.

Command: `python analysis/AX_drug_model.py`

## BZ. Two-component drug effect model with AUC-based tachyphylaxis

The drug effect on intake decomposes into an acute response to blood-level change plus a tonic suppression that erodes with cumulative receptor exposure (AUC), not with calendar weeks.

**Model.**

    drug_effect = 2×tonic × (blood − EMA₁ᵈ(blood)) + tonic × blood × exp(−K_AUC × ∫blood dτ)

The acute term responds to blood-level *changes* relative to a 1-day adapted baseline. It captures the injection-day sawtooth and does not decay — consistent with GLP-1 receptor recycling (β-arrestin desensitisation recovers within ~6 hours; Jones et al., doi:10.1038/s41467-018-03941-2). The tonic term captures sustained mean suppression; **its erosion is proportional to the integral of blood level over time**, not to elapsed calendar time. Receptor desensitisation scales with how much agonist the receptor has been exposed to, not with how many weeks have passed.

**Subject fit (1 free param + 1 decay param).** tonic = −63 cal/unit, K_AUC = 0.00337 / unit-week of cumulative exposure (AUC half-life ≈ 206 unit-weeks). RMSE ≈ 330. Sawtooth prediction by era is preserved within 20% across the 80-week subject history.

**Why AUC, not calendar weeks.** The earlier calendar-time form (`exp(−0.0126 × weeks)`) gave every dose arm the same tachyphylaxis at any given week, regardless of how much drug was actually present. That predicted a slope across SURMOUNT-1 dose arms (sim − published of +7.3, +5.4, +0.5 at 5/10/15mg) that has the exact pattern of an AUC mis-specification: low-dose arms over-eroded, high-dose arms on the nose only because their AUC happened to match the calibration AUC. Switching to AUC fixes the slope.

**Trial validation (no set point, pure drug + MSJ-scaled TDEE; subject K_AUC carried over).**

| Trial | Dose | Simulated | Published | Δ |
|---|---|---|---|---|
| SURMOUNT-1 | 5mg | −13.9% | −16.0% | +2.1% |
| SURMOUNT-1 | 10mg | −18.1% | −21.4% | +3.3% |
| SURMOUNT-1 | 15mg | −19.5% | −22.5% | +3.0% |
| SURMOUNT-2 | 10mg | −18.8% | −12.8% | −6.0% |
| SURMOUNT-2 | 15mg | −20.3% | −14.7% | −5.6% |
| SURMOUNT-4 | 15mg | −14.4% | −20.0% | +5.6% |

Best trial tonic: −113 cal/unit. SSE = 123. The dose-response shape across SURMOUNT-1 is now flat at ~3% under-prediction across all three arms — a uniform shift rather than the steep ramp the calendar-time form produced. SURMOUNT-2 over-predictions persist (T2D diverts GLP-1 to glucose control; not in the model). SURMOUNT-4 under-predicts by 5.6% over a 36-week window — short relative to the AUC erosion timescale.

**Honest acknowledgement of the 15mg shift.** Under calendar-time tachyphylaxis the 15mg arm fit within 0.5% of published. Under AUC-based tachyphylaxis it sits at 3.0% under. The earlier near-perfect fit at 15mg was partly accidental — the calibration AUC happened to coincide with the 15mg arm's AUC. AUC-form is more biologically defensible (occupancy-driven receptor desensitisation) and produces a more honest distribution of error across doses.

**Subject vs trial tonic.** Subject −63 (logged-calorie basis); corrected for 7% undercount: −68. Trial best −113. CE quantifies this gap as a central-GLP-1-pathway dampening factor (below).

Command: `python analysis/BZ_drug_twocomponent.py`, `python analysis/AV_surmount_simulation.py`

## CE. Drug-response dampening: a one-number summary of individual variation

A single multiplier `m` applied to the central-pathway tonic (with acute = 2 × tonic following) characterises individual variation in apparent GLP-1 weight-loss response on a continuous scale.

**Method.** Hold the BZ AUC-form drug model fixed (K_AUC = 0.0034 / unit-week, acute/tonic ratio = 2:1, all trial-population body composition and demographics). Find the scalar `m` such that `tonic = m × tonic_NonDiabetic` best fits each cohort's published weight loss across all dose arms.

**Subject result.** The subject's implied `m = 0.52` (intake-fit, undercount-corrected). At their `m`, the model predicts ≈ −11% body-weight loss at 72 weeks under SURMOUNT-1 conditions and 10 mg dose; the subject's actual tirz-era loss was −13% over ~80 weeks. Cross-modality check passes within 2 pp — the dampening factor is real, not a fitting artifact of intake or weight alone.

**SURMOUNT-2 fits at `m = 0.605`** but this is mechanistically misleading and should not be over-interpreted. SURMOUNT-2's HbA1c-stratified post hoc analysis ([Patoulias et al. 2025](https://diabetesjournals.org/care/article/48/11/e136/163253/Weight-Reduction-With-Tirzepatide-Varied)) shows the T2D weight-loss attenuation is **dominated by glucosuria reduction**, not central-pathway scaling: at HbA1c <7% (low glucosuria) T2D weight loss matches non-diabetic SURMOUNT-1 (−17.7% vs −17.8%); the apparent "T2D dampening" is the disappearance of urinary calorie loss as treatment normalises HbA1c. **A scalar on `tonic` is the wrong functional form for the T2D effect** — it numerically fits the cohort means by absorbing the glucosuria correction into a single parameter, but the underlying mechanism is peripheral, not central. The subject is non-diabetic and has no glucosuria; their `m = 0.52` is genuinely a central-pathway summary, but the previous framing of "more dampened than T2D" was a category error.

**What `m` is and isn't.** A useful collapse of an individual's observed cohort-relative weight-loss response into one number. Candidates for what generates the subject's particular `m`: MC4R heterozygosity and FTO homozygous risk (both documented in BACKGROUND); hypothalamic GLP-1R density variation; baseline ghrelin/leptin profile (not available for this subject); genetic variants identified in [Lasky-Su et al. 2026 GWAS, n=27,885](https://www.nature.com/articles/s41586-026-10330-z). Not testable from this dataset alone; flagged for future bloodwork.

**Generalisation.** A patient's `m` could in principle be estimated from (a) intake response over the first few weeks of dosing, (b) early weight-loss trajectory vs trial cohort mean, (c) genotype + phenotype priors. Low-`m` patients (the subject's category) may need higher doses or alternative agents. Whether `m` is a stable individual property or shifts with metabolic state is open.

Command: `python analysis/CE_drug_dampening.py`

## CF. Trial regain vs SmoothLatch: the model needs `m`-scaled pressure

The eClinicalMedicine 2026 meta-regression ([Wilding et al.](https://www.medrxiv.org/content/10.1101/2025.06.09.25328726v1.full), n=3,236 across 6 GLP-1RA withdrawal RCTs) fits a universal exponential recovery `R(t) = A(1 − exp(−kt))` with **A = 75.6%** and **k = 0.0302/wk (HL = 23 wk)**. They report no per-trial variation, no moderator analysis, and explicitly flag speed-of-loss-during-treatment as an untested predictor.

Carrying our SmoothLatch (`tol=3, hold=14, rate=0.01, pressure=55 cal/lb`) onto these 6 trials produces a **mean R(52) ≈ 14% — about 3× under the published 60%**. Diagnosis: at typical trial loss rates (~0.15 lb/day FM), FM moves only ~2 lb in 14 days — well within the tol=3 tolerance. The latch fires every day during loss; SP follows FM down; the SP-FM gap at cessation stays small (~10-20 lb); pressure × that gap can't drive the published regain.

**No single SmoothLatch parameter set fits both the subject and the trials.**

| latch (tol, hold, rate) | subject r (180d surplus) | mean trial R(52) |
|---|---|---|
| 3, 14, 0.01 (subject-fit) | **−0.94** | 14% (under) |
| 1, 30, 0.0025 (trial-fit) | −0.60 (broken) | 60% (matches meta) |
| 1.5, 21, 0.005 (compromise) | −0.84 | 28% |

**Pressure-only scaling can't bridge the gap** (raising pressure 55 → 220 cal/lb under the loose latch only moves trial R(52) from 14% to 20% — the latch closes the gap before pressure can act).

**Velocity-gated latch is partial** (SP fires only when `|dFM/dt| < ~0.05 lb/day`): subject r ≈ −0.84 (worse than original), trial mean R(52) ≈ 46-58%. Better than the rigid alternative but still not clean.

**Most plausible structural reconciliation: pressure scales with central-pathway responsiveness** — the same dampening factor `m` from CE governs both drug tonic and lipostat pressure. Subject `m = 0.52` ↔ subject pressure ≈ 55 cal/lb; non-diabetic trial cohort `m ≈ 1.0` ↔ trial pressure ≈ 106 cal/lb. Under population-specific pressure with the loose latch retained, both the subject's intake fit and the trial cohorts' regain magnitudes become explainable by a single mechanism: **central GLP-1 / melanocortin pathway responsiveness scales both directions of the lipostat signal — the body's pressure when below SP and its sensitivity to drug-induced satiety from above.**

This is a structural enhancement to the set point theory, not a free parameter: the dampening multiplier was already in the model from the drug side (CE), and it predicts that the lipostat side should be scaled identically. The eClinicalMedicine meta-regression's 60% population mean is what `m = 1.0` cohort behaviour produces; this subject's predicted regain on hypothetical discontinuation is much smaller because their lipostat is also dampened.

**Speed-of-loss → asymptote correlation survives every parameterisation tested** (r = +0.6 to +0.95 across all latch/pressure choices). Faster loss → larger SP-FM gap at cessation → larger eventual regain. The published meta-regression treats this as an untested gap; our model says it's a real moderator that should explain per-trial variation in future analyses.

**Per-trial predicted regain at 52 wk** (default subject params, illustrative — population-scaled prediction would be ~3× higher):

| Trial | lead-in (wk) | loss rate (lb/wk) | predicted R(52) |
|---|---|---|---|
| SCALE liraglutide 3mg | 56 | 0.32 | 11% |
| STEP-1 sema 2.4mg | 68 | 0.51 | 9% |
| SURMOUNT-1 15mg | 72 | 0.72 | 10% |
| STEP-4 sema 2.4mg | 20 | 1.23 | 23% |
| SURMOUNT-4 15mg | 36 | 1.38 | 18% |

The fastest-loss trials (STEP-4, SURMOUNT-4) predict the largest regain even at this subject's dampened parameters — consistent with the published trend that SURMOUNT-4 sees the largest absolute regain in the meta-set.

Command: `python analysis/CF_regain_vs_meta.py`

## Asymptote of GLP-1 weight loss

Why does GLP-1 / GIP weight loss flatten? Trials almost always plateau by week 60-72 — the SURMOUNT-1 15mg curve runs ~−15% by week 24, ~−21% by week 48, ~−22.5% by week 72; STEP-1 sema 2.4mg runs ~−10% by week 20, ~−14.9% by week 68. The asymptote is dose-dependent and largely durable while the drug is continued. Mainstream framings ("metabolic adaptation", "the body fights back") name a phenomenon without naming a mechanism. This dataset's parameters give a quantitative one.

**Three timescales race.**

1. **Drug suppression** D(t) = tonic × blood(t) × exp(−K_AUC × ∫blood dτ). The drug provides exogenous satiety on a millisecond timescale, eroded over months by GLP-1R β-arrestin desensitisation and downstream hypothalamic adaptation as cumulative receptor exposure (AUC) grows. Subject K_AUC ≈ 0.0034 / unit-week; AUC half-life ≈ 206 unit-weeks (BZ).
2. **Endogenous eating pressure** P(t) = 55 × (SP(t) − FM(t)), symmetric. Instantaneous response of leptin / melanocortin signalling to current FM against the un-recalibrated receptor baseline. Positive below SP (drives intake up); negative above (suppresses intake).
3. **Set-point recalibration** dSP/dt = 0.01 × (FM − SP), gated on FM having held within ±3 lbs for 14 days (SmoothLatch, above). Receptor turnover only resets to a new baseline when chronically stimulated near that baseline — and chronic GLP-1/GIP agonism *bypasses* the endogenous receptor, so during active loss the latch is open and the SP doesn't move.

**The asymptote in one equation.** At equilibrium FM is stable, daily intake = TDEE, and pressure equals drug:

    55 × (SP − FM_eq) = D₀ × tachy(t_eq)

So at the asymptote, the gap that the drug holds open is

    SP − FM_eq = D₀ × tachy(t_eq) / 55 cal/lb

The equation gives the *residual gap held by the drug at equilibrium*, not the total weight loss. For SURMOUNT-1 15mg (tonic = −113 cal/unit trial fit, mean blood ~13 units at steady state): D₀ ≈ 1470 cal/day at week 0. Tachyphylaxis is now AUC-based: at week 48, accumulated AUC ≈ 600 unit-weeks, tachy ≈ exp(−0.0034 × 600) ≈ 0.13, so D_eff ≈ 190 cal/day at the time the curve flattens. Residual gap = 190 / 55 ≈ 3.5 lb FM. Total FM loss is much larger — it's the cumulative deficit along the trajectory plus whatever the SmoothLatch has managed to drag the SP down during the stable phase near the asymptote. BZ's full forward integration produces SURMOUNT-1 15mg at −19.5% body weight (published −22.5%, Δ +3.0%) with the AUC-form tachyphylaxis carried over from the subject fit. The toy equation captures *why the curve flattens*; BZ captures *how the trajectory got there*.

**Why dose maps to asymptote depth.** D₀ scales with mean steady-state blood level, which scales with dose. Each doubling of dose roughly doubles D₀, which roughly doubles the gap held at equilibrium, which deepens FM_eq by the same multiple. SURMOUNT-1: −9% at 5mg, −16% at 10mg, −22.5% at 15mg. Sublinear because tachyphylaxis erodes the higher-dose advantage faster, and because the SmoothLatch latches fastest near the asymptote (when FM is most stable).

**Why duration on drug doesn't deepen the asymptote much past 60-72 weeks.** Tachyphylaxis erodes D(t) at the same rate the latch closes the SP-to-FM gap. The two cancel: the drug loses effective potency, but the SP also stops fighting back as it slowly recalibrates. The system hovers near equilibrium. SURMOUNT-1 gains <2% additional loss between week 48 and week 72.

**Why discontinuation produces large regain.** The latch only fired during the brief stable equilibrium near the asymptote; for the bulk of the loss, FM was changing fast enough to keep the latch open. So at discontinuation D drops to zero immediately, but SP is anchored not far below the pre-treatment level. Pressure = 55 × (SP − FM) reactivates at full strength against an un-recalibrated baseline. SURMOUNT-4 +14% in 52 weeks. STEP-1 ~+10%. The SmoothLatch's regain numbers (+14.1%, +9.3%) hit both within a percentage point.

**The asymptote is the receptor-equilibrium gap that the drug can hold open, against a set point that's barely moved.** The patient is not "stuck"; the patient is at the dynamic equilibrium between exogenous satiety and endogenous defence of an un-recalibrated baseline. This makes three concrete, falsifiable predictions that distinguish this account from "metabolic adaptation":

- *Asymptote depth scales with D₀ / (55 cal/lb).* Lower-pressure phenotypes (set point that adapts faster) reach deeper asymptotes at the same dose. Tested in part by between-subject heterogeneity in SURMOUNT-1.
- *Asymptote depth is independent of baseline TDEE adaptation.* The drug's asymptote does not depend on the body "fighting back" via metabolic rate suppression. AJ in this dataset shows the metabolic channel is reversed on-drug — yet the asymptote is reached.
- *Maintaining the asymptote for 12+ months should produce a deeper effective set point at discontinuation.* The latch needs sustained FM stability to fire. SURMOUNT-4 stopped at peak loss; SURMOUNT-3 maintenance arms with longer stable holds should regain less. (Not yet published with the timescale needed to test; testable now in clinical practice.)

The receptor-dynamics frame says the asymptote isn't a failure of the drug or a triumph of homeostasis — it's the equilibrium *both* sides settle into when one is acting through an exogenous bypass and the other is waiting for chronic receptor exposure that never arrives.

Command: `python analysis/BZ_drug_twocomponent.py`, `python analysis/AV_surmount_simulation.py`, `python analysis/AZ_sp_model_search.py`

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

**Nulls.** Drug blood level → baseline temperature|FM: r = -0.03 (zero; raw r = +0.45 was season-confounded). Pre-tirz vs on-tirz baseline at matched FM: Δ = -0.07°F (negligible).

**TDEE → baseline temperature: real but small.** Partial correlation controlling for FM is **r = +0.11-0.13** for TDEE residual, TDEE/RMR ratio, and absolute TDEE at daily resolution. The signal isn't strong but it's real — temperature does track the metabolic rate state to a modest degree.

Command: `python analysis/CD_secondary_recheck.py` (re-runs secondary analyses against the current Kalman)

Command: `python analysis/AS_temperature_baseline.py`, `python analysis/AT_temperature_retest.py`, `python analysis/X_temperature_phase.py`
Artifact: `analysis/AS_temperature_baseline.csv`, `analysis/AS_circadian_curve.csv`, `analysis/X_temperature_daily.csv`, `analysis/X_temperature_phase_overlap.csv`

## AE. Sleep and energy balance (null)

2,057 sleep measurements (median 7.8h, std 1.5h, schedule 3am-11am). Sleep duration → next-day calories: r = -0.01. Sleep → next-day steps: r = -0.02. No trailing window of sleep duration predicts intake, TDEE, or steps. The schedule is too consistent (CV 0.19) to detect effects that require cross-sectional variation.

**An extra hour of sleep still does not show an independent effect.** Recasting the question as accumulated sleep debt does not rescue the signal. Raw sleep-debt correlations with metabolic proxies are negative, but they vanish after controlling for bedtime timing. Example: 30-day sleep debt vs TDEE residual is raw r = -0.24, partial r = +0.03 controlling for bedtime, FM, weekend, and year. The apparent "sleep debt lowers RMR" pattern was mostly shared structure with later bedtimes.

**Late bedtime matters a little; sleep duration does not.** Bedtime timing carries a small but persistent signal on the daily metabolic proxies even after controlling for sleep duration. Later bedtimes predict slightly lower TDEE residual / TDEE-to-RMR ratio, with the strongest effect on trailing windows: same-day bedtime partial r = -0.04, 7-day bedtime -0.06, 14-day -0.08, 30-day -0.10. Baseline temperature does not track the same signal (baseline temp vs TDEE residual partial r = +0.01).

**Interaction null.** Combining the two ideas does not help. The late-bed × short-sleep interaction is still null for TDEE residual (partial r = +0.01), TDEE/RMR ratio (+0.01), and same-day calories (-0.02). In this dataset, the independent signal is weak circadian timing, not quantity or debt.

Command: `python analysis/AE_sleep_null.py`, `python analysis/BS_sleep_timing_split.py`
Artifact: `analysis/BS_sleep_timing_cutoff_sweep.csv`
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

## BQ. Creatine water-retention signature (underpowered)

Creatine monohydrate causes intracellular water retention (~2-5 lb at 5 g/day), which the Kalman should attribute to lean mass (the LM state observes composition scans during loading). Documented creatine periods: 2 g from 2018-10-29 to 2021-04-28, then 5 g from 2025-05-09. Tests:

- **2018-10-29 onset (2 g):** not testable — only 1 weight reading within ±60 days.
- **2021-04-28 cessation (2 g):** marginal — regression discontinuity −1.75 lb at offset −10d (p=0.057). Direction correct, near detection threshold for the small dose.
- **2025-05-09 onset (5 g):** confounded — full-window RD +2.92 lb (p<0.001), but the tirzepatide dose escalated 7.5→10 mg on 2025-05-27 (+18 d), dominating the signal. Creatine loading kinetics (~2-4 weeks) overlap the dose change.

Verdict: cannot confirm or reject the expected ~2-5 lb signature. A clean detection would require a creatine washout during stable weight and stable medication, with composition scans bracketing the loading window so the LM state can absorb the water rather than letting it leak into FM.

Command: `python analysis/BQ_creatine_signature.py`

## CA. Samsung stress data (sparse, suggestive only)

Samsung Health logs a daily stress score, but coverage is sparse: 962 session rows collapse to only 208 daily-stress days, almost all in 2021-2022. Within that window, same-day binge risk shows a small positive association with stress, but it rests on only 16 binge days. Stress correlates more with later bedtime than with sleep duration — consistent with AE's "late bedtime is the live circadian variable" reading. The series is not dense enough to anchor major conclusions; flagged here so future Samsung pulls with denser stress coverage can be re-evaluated.

Command: `python analysis/CA_stress_daily.py`

# Diet experiments

## O. Epoch analysis, potato diets, and travel

Do the 22 manually annotated diet epochs in `intake/diet_epochs.csv` separate meaningful regimes? Do potato diets have special metabolic properties?

Travel binges: 2716 cal/day, 32.1% binge rate, +9.5 lbs fat. Keto: 2163 cal, 110.5g protein, 2.1% binges. Potato diet (4 epochs, 69 days): 1930 cal, 48.6g protein, 0% binges, -7.3 lbs fat.

Potato diets: TDEE/RMR 1.168 before, 1.159 during (unchanged), 1.098 after. Post-potato rebound: 2629 cal/day, 19.6% binge rate. Potatoes are a monotone binge-suppressing cut, not a uniquely expenditure-preserving intervention.

**Jordan trips.** Two month-long travel periods tested whether high intake produces compensatory TDEE increases. 2015: intake 3262 cal/day, TDEE 2071 (-35 vs pre-trip), +10.3 lbs fat. 2019: intake 2529, TDEE 2021 (-21 vs pre-trip), +4.4 lbs fat. No compensatory TDEE surge. The trips are binge windows without an energy-out bonus.

Command: `python analysis/O_diet_epoch_analysis.py`, `python analysis/P_jordan_trip_analysis.py`
Artifact: `analysis/O_diet_epoch_summary.csv`, `analysis/O_diet_epoch_family_summary.csv`, `analysis/O_potato_epoch_window_summary.csv`, `analysis/O_potato_epoch_contrast.csv`, `analysis/P_jordan_trip_summary.csv`, `analysis/P_jordan_trip_delta.csv`

# Kalman filter appendix

The Kalman filter underlies every TDEE, FM, and LM estimate cited above.

**Four states.** FM (lbs), LM (lbs), BUF (lbs — bounded transient buffer for non-glycogen, non-sodium water/gut variation), TDEE (cal/day, free random walk).

**Key constants.**
- Forbes-coupled energy partition: `deficit = dFM × 3500 + dLM × 1100`, `dFM/dW = ff(FM, workouts/wk)` where `ff_base = 1 − C²/(C+FM_kg)²`, C = 10.4 kg, and `ff = clip(ff_base + 0.084 × workouts_per_week_30d, 0, 1)` per CB. The workout modulation shifts the partition toward fat preservation of lean during cuts (and toward lean during surplus), with α = 0.084 derived from CB's run-level regression. Set α = 0 to recover plain Forbes.
- FM rate caps: gain ≤ +0.5 lb/day (Levine NEAT, [doi:10.1126/science.283.5399.212](https://doi.org/10.1126/science.283.5399.212)); loss ≤ FM × 31/3500 lb/day (Alpert mobilization ceiling — auto-tightens at low FM).
- LM rate caps: gain ≤ +0.15 lb/day (AA training term + headroom); loss ≤ -0.5 lb/day (extreme protein turnover, rarely hit).
- BUF: HL = 5 days, hard cap ±3 lb, daily noise std ≈ 0.4 lb. Mean ≈ 0, 99% within ±1 lb in practice; max excursion in 15 years was 1.65 lb.
- TDEE: σ ≈ 50 cal/day, hard daily-step cap ±100 cal/day. Free random walk (no mean reversion, no anchor).
- Composition scans observe FM and LM directly with σ ≈ 1.5 lb each.
- Calorimetry RMR observations are **held out** — never used by the filter; they're a posterior validation channel.

**Diagnostics on current fit.** Innovation lag-1 ACF = 0.08 (close to white). FM scan MAE = 1.23 lb (max error 5.1 lb). TDEE/RMR fraction in [0.95, 1.30] band = 92% across 25 calorimetry measurements. r(TDEE, baseline temperature) partial⏐FM = 0.14 daily, 0.26 on 7-day rolling. Cumulative energy-balance closure = -3.0 lb over 15 years.

Command: `python analysis/P4_kalman_filter.py`
Artifact: `analysis/P4_kalman_daily.csv`, `analysis/P4_kalman_sweep.csv`
