# Findings

Each section has a `Command:` line pointing to the standalone script that reproduces the numbers. Untested hypotheses are in [ROADMAP.md](ROADMAP.md).

# Data quality

## A. Energy balance quality

Does the intake data close against weight and calorimetry?

Cumulative energy balance residual: ±5.9 lbs over 15 years (per-year drift ≈ 0). Residual undercount after the canonical 1.09 intake multiplier (composition-aware model, 25 calorimetry × 70 composition measurements): gaining 3.0%, losing 1.6%, stable 1.3% — approximately uniform across phases.

Command: `python analysis/A_energy_balance_quality.py`

## Z. TDEE by year

Does derived TDEE vary with weight, or does the body defend a band?

| Year | TDEE | Intake | Gap | Fat (lbs) | TDEE/RMR |
|---|---|---|---|---|---|
| 2011 | 2732 | 1772 | +961 | 100→55 | 1.255 |
| 2012 | 2390 | 1917 | +473 | 55→39 | 1.124 |
| 2013 | 2306 | 1847 | +459 | 39→22 | 1.110 |
| 2014 | 2288 | 2105 | +183 | 22→35 | 1.133 |
| 2015 | 2318 | 2130 | +188 | 35→47 | 1.184 |
| 2016 | 2308 | 1954 | +354 | 47→47 | 1.194 |
| 2017 | 2399 | 2140 | +259 | 47→51 | 1.229 |
| 2018 | 2099 | 2005 | +94 | 51→64 | 1.061 |
| 2019 | 2096 | 2018 | +78 | 64→68 | 1.031 |
| 2020 | 2357 | 2120 | +238 | 68→71 | 1.186 |
| 2021 | 2334 | 2136 | +198 | 71→74 | 1.176 |
| 2022 | 2460 | 2224 | +236 | 74→77 | 1.241 |
| 2023 | 2505 | 2267 | +238 | 78→86 | 1.260 |
| 2024 | 2592 | 2274 | +318 | 83→79 | 1.289 |
| 2025 | 2320 | 1943 | +377 | 79→62 | 1.181 |
| 2026 | 2163 | 1844 | +319 | 62→61 | 1.108 |

Start→end of year. **TDEE numbers are scaled to actual energy expenditure** (canonical filter applies a 1.09 intake multiplier per CI to absorb the ~9% logged-intake undercount documented in finding A; cumulative EB closes within ±1 lb over 15 years). Notable: **2018-19 keto era TDEE/RMR drops to 1.03-1.06** — markedly suppressed relative to surrounding years (2017 at 1.23, 2020 at 1.19), though no longer below sedentary maintenance. 2024 drug-era at 1.29 is the all-time high (subject more active in walks/strength + drug-era TDEE artefacts). The "Gap" column is intake − TDEE; negative gaps imply weight loss, positive imply gain.

Command: `python analysis/Z_tdee_by_year.py`
Artifact: `analysis/P4_kalman_daily.csv`

# Metabolic rate adjustment

The year-by-year table shows TDEE falling from 2661 to 2044 as fat mass dropped from 106 to 25 lbs — partly from losing metabolically active tissue, partly from metabolic adaptation. TDEE barely recovered during a decade of regain. This raises the question: at the same body composition, does TDEE differ depending on direction?

> **Kalman attribution caveat (CY/CZ3, applies to K, B, AI, AJ, AL, and other direction-conditional TDEE analyses below).** P4 has TDEE as a free random-walk state with no anchor; over any window it equals `intake − ΔFM·3500 − ΔLM·1100 − ΔBUF` mechanically. When ΔFM is negative, TDEE is mechanically attributed up; when ΔFM is positive, attributed down. Given the subject's asymmetric rise/fall rates, this builds in ~500 cal/day of "K-like" asymmetry by construction. Under a constrained filter (`TDEE = c·RMR + d·steps + small drift`), the K direction effect collapses to chance (CY: 10/19 cycles vs 14/19 baseline; median Δ TDEE/RMR +0.001 vs +0.032). CZ3 confirmed no step-anchored filter beats baseline FM tracking, so the canonical filter at σ_TDEE=50 stays — but direction-conditional analyses on its TDEE inherit this bias. Findings below that depend on filter TDEE direction (K, AI most directly) need separate validation against held-out calorimetry, which we have not run.

## K. TDEE hysteresis

At the same fat mass, does TDEE differ depending on whether weight was reached from above (falling) or below (rising)? [Leibel et al. (1995)](https://doi.org/10.1056/NEJM199503093321001) demonstrated this in a metabolic ward (n=18, ±10% body weight). This dataset replicates the finding free-living with 1,482 FM-matched pairs — **but the magnitude reported below is upper-bounded; CY shows much of it is Kalman attribution rather than real metabolic asymmetry.**

Matched fat-mass bands (pre-tirzepatide):
- FM 25-45 lbs: rising 2235 vs falling 2336 (+101)
- FM 45-65 lbs: rising 2212 vs falling 2474 (+262)
- FM 65-85 lbs: rising 2374 vs falling 2659 (+286)

Regression: TDEE = 2145 + 3.90×FM + 124×falling − 90×rising. Falling phase adds ~125 cal; rising phase compresses by ~90 cal — the metabolic channel is asymmetric across direction.

Robustness: 1,482 FM-matched pairs (within 2 lbs), mean difference +179 cal (falling − rising). Pair lag-1 autocorrelation 0.83, effective n = 142. Bootstrap 95% CI: [+104, +267]. p < 10⁻²². This is the most statistically robust finding in the dataset.

**No universal dead zone.** A corollary: there is no single intake band (2000-2500 cal) where metabolism adapts to maintain weight. In falling phase, TDEE is flat across all intake levels (1600-3000 cal). In stable phase, higher intake weakly raises TDEE (coef +0.018). In rising, near zero (+0.005). The pattern is branch-specific expenditure compression, not one dead zone.

**Survives matching on fat mass.** Within matched start-fat bands: FM 50-65, falling post-pre -0.0480 vs rising +0.0340. Regression controlling for start fat mass: falling_vs_rising = -0.0534 for net post-run failure.

**Survives matching on metabolic state.** Matching on both fat_mass_start and pre-run TDEE/RMR: 43 pairs, falling minus rising = -0.0442 for net post-run failure. Tighter matches (distance ≤0.50, n=34): -0.0462. Regression controlling for both: falling_vs_rising = -0.0398.

Command: `python analysis/K_tdee_hysteresis.py`, `python analysis/I_deadzone_phase.py`, `python analysis/S_metabolic_failure_matched.py`, `python analysis/T_metabolic_failure_state_matched.py`
Artifact: `analysis/K_tdee_hysteresis_phase_summary.csv`, `analysis/K_tdee_hysteresis_band_summary.csv`, `analysis/K_tdee_hysteresis_regression.csv`, `analysis/I_deadzone_phase_bin_summary.csv`, `analysis/I_deadzone_phase_regression.csv`, `analysis/S_metabolic_failure_matched_bands.csv`, `analysis/S_metabolic_failure_matched_regression.csv`, `analysis/T_metabolic_failure_state_match_summary.csv`, `analysis/T_metabolic_failure_state_match_pairs.csv`, `analysis/T_metabolic_failure_state_match_regression.csv`

## B. Weekend fasting as metabolic adjustment microcosm

Do acute caloric deficits produce lasting fat loss? Seven consecutive Sat-Sun 36-hour fasts (Oct-Nov 2019). Mean deficit per fast: 3,300 cal. Kalman FM change Fri→Mon: -0.76 lbs (80% of expected). Kalman FM change Fri→Fri+7: +0.08 lbs. Post-fast weekday intake: 2,369 cal (78 cal/day below pre-fast 2,447). Zero compensatory overeating — yet the deficit vanishes within a week. **Reading**: deficit is absorbed somewhere within ~7 days; original prose said "TDEE drops to absorb the deficit", but per the Kalman attribution caveat that "TDEE drop" is partly mechanical. The robust observation is that intake stays low (no compensatory overeating) and FM returns to baseline by Fri+7. *Where* the energy goes (TDEE drop vs BUF re-equilibration vs intake mis-counted on fast days) cannot be cleanly disentangled from filter TDEE alone.

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

> **Status: framework under revision. The receptor-dynamics SmoothLatch reading below is preserved for reference but several load-bearing pieces are now known to be artifactual or misinterpreted.** Specifically:
>
> - **AG sign convention is reversed in the prose.** AG's binned data show binges at NEGATIVE sp_distance (FM ABOVE its trailing EMA = regain phase) at 8.7%, falling monotonically to <1% at POSITIVE sp_distance (below EMA = restriction phase). The narrative throughout this section says "below SP drives binges"; the data say "above the trailing EMA drives binges." This is consistent with AC+BK's post-floor regain story, but it inverts the textbook lipostat reading. The sigmoid AG fits (with positive k and inflection at +19.5 below SP) is mathematically backwards from the bin data.
> - **CN: no fixed defended weight.** Maintenance FM (90-day windows with mean surplus < 50 cal/day) drifts monotonically from 30 lb (2013-14) → 55 lb (2015-19) → 72 lb (2020-24). The body settles wherever the lifestyle of the decade puts it.
> - **CU/CV: SP is mathematically unidentifiable from cycle mass-balance.** Subject's BW max in this dataset is 257 lb; if SP > 257, the threshold `max(0, SP − BW)` never activates and SP becomes degenerate in the linear drain model. The user's hypothesis SP ≈ 285-300 (subject's pre-data peak weight) is *consistent* but *not preferred* over arbitrary higher values.
> - **CW: no constant willpower budget.** Across 19 falls, no simple drain function (deficit, FM, FM×deficit, deficit/FM, BW×deficit) produces a constant total drain across cycles. CV of cumulative drain at trough ranges 0.71-2.05 — not even close to constant. SP-and-W free-per-cycle fits succeed only by force-fitting and SP_i tracks mean BW per cycle.
> - **AG/AM autocorrelation problem.** EMA-distance and rolling-surplus share input bandwidth (BR ridge); much of AG's r=−0.62 is inherited from this. AR(intake)+FM controls reduce SP-distance partial r to −0.18 (CN Test 1) — meaningful but ~⅓ the marginal headline.
> - **K hysteresis is partly Kalman bookkeeping** (CY/CZ3): under constrained TDEE, the K direction asymmetry collapses to chance.
>
> What survives: a slow-moving operational anchor (mean BW per cycle drifts upward over 13 years), trial-regain incompleteness pointing at *some* defended-state behavior in published cohorts, and the AC+BK post-floor binge eruption as a real phenomenon. What's open: whether any of the SmoothLatch parameters or the EMA-SP framework correspond to anything biologically real for this subject. The text below is the pre-revision framing; treat with the caveats above.

K and B show the body adjusts metabolic rate by direction. Is there a single defended fat mass driving both responses? [Speakman (2011)](https://doi.org/10.1242/dmm.008698) proposes dual intervention points with a "zone of indifference." The data finds no such zone — the response is graded.

**Methodological frame.** Any model that defines the set point as a smoothing of recent fat mass and is then evaluated against a smoothing of recent surplus shares input bandwidth with what it tries to predict (BR makes this explicit, below). To avoid that ridge, the set-point form is identified primarily from **published GLP-1 withdrawal trials** — where this subject's intake autocorrelation cannot bias the answer — and only afterwards confirmed against this subject's binge rate, intake-free FM, and weigh-in-only outcomes. The continuous-surplus regression (AM/AP) is treated as consistent supporting evidence, not a primary identification.

**Receptor-dynamics reading.** The model that survives this protocol has the shape of a two-receptor adaptation system: the set point is the equilibrium of slow leptin/melanocortin receptor populations that up-regulate when chronic exposure is sparse and down-regulate when it is dense. The 45-50 d half-life is the receptor-turnover and AgRP/POMC neuropeptide-equilibration timescale. Adaptation requires the receptor to be repeatedly stimulated near a new equilibrium — explicitly *not* during rapid drug-driven loss, where the receptor is bypassed by exogenous GLP-1/GIP agonism. This frames every section below.

## What this is, by elimination

The SmoothLatch + receptor-dynamics reading is what's left after testing every other candidate set-point form against this dataset. The remainder of this section establishes each elimination — but the chain is short enough to summarize first.

| Candidate | Source | Falsifying finding | Verdict |
|---|---|---|---|
| Pure CICO, no homeostat | implicit in textbook | Energy balance closes (A: ±5.9 lb/15 yr) — but TDEE depends on direction (K: regression coefficients +124 cal during falling phases, −90 cal during rising at matched FM, spread = 214 cal/day, p < 10⁻²²); ratio drops to 0.97 in restriction, recovers to 1.19 when restriction eases (Z). The body actively defends. | Insufficient as a model |
| Gravitostat (foot-pounds drive intake) | [Jansson et al. 2018](https://doi.org/10.1073/pnas.1800033115) | Foot-pounds → next-day intake: r = +0.050 partial. Wrong direction (more weight-bearing → *more* eating). Steps → weight change: r = +0.010 at every timescale (N, U). | Falsified |
| Dual intervention points with "zone of indifference" | [Speakman 2011](https://doi.org/10.1242/dmm.008698) | Sigmoid binge response across the gap — continuous gradient, no flat zone. Surplus is graded across all gap sizes (AM). | Falsified |
| Set point tracks scale weight | gravitostat-adjacent | At every half-life, fat mass beats total weight beats scale: r = -0.62 / -0.57 / -0.46. The body weighs *fat*, not total mass. Consistent with [Kennedy 1953](https://doi.org/10.1098/rspb.1953.0009) and the leptin system. | Falsified |
| Fixed set point at any specific weight | classical "set point theory" | Best fixed defended weight: r = +0.23 (and wrong sign). Best moving 60d EMA: r = -0.62. Difference 3×, sign-reversed. | Falsified |
| Fixed weight-suppression proxy | [Lowe et al. 2007](https://doi.org/10.1002/eat.20405) (bulimia binges) | Lowe's directional finding (suppression below historical max → binges) replicates here, but their fixed-historical-max proxy is dominated by a moving 50d EMA at every half-life tested (AG). Right intuition, wrong reference. | Form falsified, signal real |
| Symmetric EMA (any half-life) of FM | extension of classical lipostat | Fits this subject's slow-gain decade at r = -0.94, but the subject-best EMA (HL=45d) undershoots SURMOUNT-4 regain at +4.1% (published +14.0%) and STEP-1 at +1.4% (~+10%). Stretching HL out to 120 d closes some of the trial gap (+12.9% / +5.1%) but degrades subject fit. No single half-life satisfies both. | Falsified by trials |
| Asymmetric ratchet (downward HL faster than upward) | data-driven hypothesis (AH original) | Improved subject fit r −0.60 → −0.69 on full dataset. Does not replicate on 2014+ natural-dynamics data (best Δr = 0.0004). The asymmetry was the EMA mechanically chasing the willpower-era descent. | Retracted |
| One-stock depletable willpower reservoir | data-driven hypothesis (BK) | Best coarse fit produces 25.6% "impossible" days where required willpower exceeds any plausible reservoir, concentrated in 2011-12. | Falsified |
| Yo-yo / variance as metabolic damage | folk wisdom; cf. [Biggest Loser reanalysis](https://doi.org/10.1002/oby.23308) | At matched mean surplus, higher 14d intake std predicts *less* fat gain (partial r = -0.20), driven by transient acute-deficit days raising TDEE. Sign-reversed. | Falsified |
| Constrained total energy expenditure | [Pontzer et al. 2016](https://doi.org/10.1016/j.cub.2015.12.046) | Walk *sessions* (not steps) raise *resting* metabolic rate at +14 cal/session (AD: CV RMSE 116 vs noise floor 170, n=23 calorimetry). RMR rises with intentional exercise; not just TEE redistribution. | Violated on RMR |
| "Metabolic adaptation" as the GLP-1 plateau | clinical/lay consensus | Names a phenomenon, not a mechanism. The plateau is fully accounted for by the three-timescale race between drug suppression D(t), gap pressure 55 × (SP − FM), and SmoothLatch recalibration (see Tirzepatide § Asymptote). | Reframed |

What's left after these eliminations: a fat-mass-tracking, slow-moving, stability-gated defended state with a 45-60 d adaptation half-life that adapts only when FM holds within a narrow band — the SmoothLatch — interacting with a fast (~15 d) metabolic channel that's asymmetric above vs below SP, and with exogenous GLP-1/GIP agonism that bypasses the endogenous receptor and therefore prevents the latch from firing during drug-driven loss. The next subsections establish each piece.

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

**Magnitude limitation: a single latch parameter set does not fit both subject and trial-cohort regain magnitudes.** The published meta-regression has R(52) = 60% mean across 6 RCTs; carrying our subject-fit latch + 55 cal/lb pressure forward predicts ~14% mean across the same trials (CF). The structural cause of the 3-4× gap is unresolved. CE rules out the simplest candidate: subject `m = 0.95` is essentially trial-matching, so an `m`-scaled pressure constant does not separate the two populations. Loosening the latch parameters (CF sweep: tol≈1.0, hold≈30, rate≈0.0025) reaches mean R(52) = 60%, but then the subject correlation degrades. The remaining open question is whether the trial cohorts and this subject differ on a parameter the current model doesn't expose — most likely the latch tolerance or hold itself, since the subject's slow-gain decade kept FM continuously stable while trial loss was rapid.

**Receptor-dynamics reading.** The 14-day stability + 3-lb tolerance is what receptor recalibration requires: leptin/melanocortin receptors adjust density and downstream signalling only after the input has held within a narrow band long enough for receptor turnover to track it. During drug-driven loss the satiety circuit is being driven by exogenous GLP-1/GIP, so the endogenous receptor isn't getting the new chronic-leptin baseline to recalibrate against. The same `m` factor that scales drug response and lipostat pressure is the central-pathway "gain" of this entire feedback loop.

Command: `python analysis/AY_sp_from_regain.py`, `python analysis/AZ_sp_model_search.py`, `python analysis/BD_smooth_latch.py`, `python analysis/BG_trial_regain_search.py`, `python analysis/CF_regain_vs_meta.py`

## Independent confirmations from this subject's data

Three measurement combinations find the same 35-60 d band without depending on rolling-surplus autocorrelation. They converge on the same answer the trial fit gives.

**Binge rate (AG; binary outcome, no surplus involved).** Distance below the moving SP predicts 90-day binge rate at r = -0.62, partial r = -0.67 controlling for absolute FM. Half-life sweep peaks at HL ≈ 60 d.

*Tracks fat mass, not weight.* At every HL: r = -0.62 (fat mass) > -0.57 (total weight) > -0.46 (scale weight). The defended quantity is fat, consistent with [Kennedy (1953)](https://doi.org/10.1098/rspb.1953.0009) and the leptin system; gravitostat-style scale-weight tracking is dominated.

*Moves, not fixed.* Best moving 60d-EMA SP: r = -0.62. Best *fixed* defended weight (any single value): r = +0.23 — wrong sign and 3× weaker. The literature has no comparable measurement: clinical sources cite vague "1-6 years" for set-point convergence with no quantification. [Lowe et al. (2007)](https://doi.org/10.1002/eat.20405) showed weight suppression predicts binge frequency in bulimia using a fixed historical-max proxy; the directional finding replicates here, but the fixed proxy is dominated by the moving EMA at every HL.

*Sigmoid response, no threshold or zone.* Continuous gradient across the gap:

| Distance from SP | Binge rate |
|---|---|
| 7.5-10 lbs below | 8.7% |
| 5-7.5 lbs below | 7.9% |
| 2.5-5 lbs below | 5.9% |
| 0-2.5 lbs below | 3.5% |
| 0-2.5 lbs above | 2.1% |
| 2.5-5 lbs above | 1.7% |
| 5-7.5 lbs above | 1.4% |
| 7.5-10 lbs above | 1.0% |

**Intake-free fat mass (P3; Kalman uses logged intake, P3 does not).** P3 derives daily FM by interpolating between weigh-ins only — no intake feeds the FM estimate. Using P3's intake-free FM on 2014+ data, the set point still predicts mean surplus at HL = 50 d (r = -0.73). The half-life is the same; only the correlation strength drops because P3's FM is noisier between weigh-ins. So the 45-60 d timescale is a property of fat mass dynamics, not a Kalman filter artifact.

**Non-circular weigh-in-only test (AU).** Building a trailing EMA "defended weight" from corrected morning weigh-ins only (no Kalman, no RMR, no intake-fed interpolation), the half-life sweep peaks at HL ≈ 35-40 d on next-30-day intake (r = -0.41, partial r = -0.40 controlling for absolute weight). The SP-distance signal beats both fixed defended weight (r = -0.13) and recent weight change. Same direction, same band of timescales, with no intake-side circularity at all.

Three independent legs — binary outcome, intake-free FM, weigh-in-only — converge on 35-60 d. Combined with trial regain (BD/AZ/BG) above, this is what supports the timescale, not the autocorrelated surplus regression.

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

## AH. Mechanics: pressure constant, symmetric adaptation

**Pressure: ~55 cal/day per lb of gap, symmetric.** From binned data across -10 to +20 lbs, remarkably constant; five pounds below = ~275 cal/day of upward pressure. Asymmetric models show zero improvement on 2014+ data (best Δr = 0.0004 vs symmetric).

**Symmetric in both directions.** On 2014+ natural-dynamics data, an asymmetric SP model (HL_down ≠ HL_up) shows zero improvement over a symmetric ~45-day half-life. (An apparent ratchet on the full dataset including 2011-13 was the EMA mechanically chasing willpower-driven descent rather than a real biological asymmetry.)

**Binge size is constant; the SP tilts the daily mean.** Binge surplus ~1427 cal regardless of distance (r = -0.02). Non-binge days: r = -0.32 between distance and surplus (+80 cal at 5+ below SP, -481 cal at 5+ above). The SP modulates background drift on every day, not just discrete events.

**Exercise independence.** Walking raises RMR (AD) but does not change the SP half-life (40-50 d for both high-walk and low-walk periods).

**Floor.** SP minimum: 16.1 lb FM on 2014-07-20. Post-floor: FM rises and binge rate ~9% over the next six months — the eruption is the receptor system coming back online above the leptin floor (see AC + BK).

Command: `python analysis/AH_set_point_properties.py`

## AI. The metabolic channel runs on a ~15-day timescale

> **Kalman attribution caveat.** This analysis uses Kalman TDEE residual against expected RMR. The TDEE residual is mechanically related to ΔFM (above), so direction-conditional TDEE residuals inherit the K bookkeeping bias documented in CY. The 15-day HL and lag-decay shape may be real, but the *magnitude* of the asymmetry is contaminated. Treat this finding as direction-suggesting, not magnitude-quantifying, until validated against held-out calorimetry stratified by direction.

The eating channel adapts with a 50-60 d half-life (above). **The metabolic channel is faster: optimal HL = 15 days** (partial r = +0.385 controlling for absolute FM, p < 0.005 at n_eff=54). Correlation is monotonically stronger at shorter half-lives, then plateaus. K's hysteresis (regression spread +214 cal/day between falling and rising at matched FM) confirms the same signal through FM-matched pairs — though K itself inherits the same Kalman attribution caveat.

| HL (days) | r (metabolic, partial) |
|---|---|
| 10 | +0.379 |
| 15 | +0.385 |
| 30 | +0.365 |
| 50 | +0.322 |
| 100 | +0.243 |

**Lag survives.** Lagging SP distance by 7-60 days at HL=15: partial r decays 0.385 → 0.357 → 0.326 → 0.268 → 0.182. A *direct* mechanical attribution from ΔFM at the same time would vanish by lag 7-14; this persists to 60 days. (Earlier prose claimed this rules out Kalman artifact; CY/CZ3 show the asymmetry magnitude is still contaminated by attribution. The persistence across lag suggests *some* signal is real, but the absolute numbers above the row are upper bounds.)

**Metabolic asymmetry.** Mean TDEE-RMR residual rises monotonically with above-SP gap: 207 cal at 2.5-5 below → 329 at 0-2.5 below → 385 at 0-2.5 above → 476 at 2.5-5 above → 612 at 5-10 above. The body actively burns more when FM is above the set point.

**Different timescales** — metabolic rate is a sprint (~15 days), eating pressure is a marathon (~50-60 days). This is the engine of the 2013 inflection: the eating channel ramped up faster than the metabolic channel could rebound, locking in the regain. It reframes ["persistent metabolic adaptation" in Biggest Loser contestants](https://doi.org/10.1002/oby.21538): if the metabolic channel adapts on a two-week timescale, persistent RMR suppression after 6 years means the subjects remained chronically below a set point that itself had moved.

Command: `python analysis/AI_expenditure_arm_timescale.py`
Artifact: `analysis/AI_expenditure_arm_sweep.csv`

## AC + BK. The 2013 inflection and the leptin floor

The 2011-13 willpower era is the cleanest test of the set-point pressure form. Under any EMA-style model, FM falling 100 → 22 lb on sustained 1200-cal restriction should produce an enormous SP gap and explosive eating pressure. It does not — binges stay near 0% through the deepest restriction (FM 13-22 lb in mid-2014), then erupt as FM begins to rise. Under a one-stock control-stock attempt to model that period (BK), 25.6% of days have "impossible" required-control values exceeding any plausible reservoir. Neither smooth-EMA nor one-stock willpower fits the willpower era.

The receptor-dynamics reading resolves both:

**Below the leptin floor (FM ≲ 15-18 lb), the eating-pressure signal is gated, not maximal.** Leptin secretion is roughly proportional to FM and saturates low at very low fat mass. Below the floor, hypothalamic input from the leptin axis is near-zero regardless of receptor density — the system enters conservation mode (TDEE/RMR collapses to 0.96-0.98) but also stops generating large eating-pressure gradients, because the gradient operates on receptor-occupancy *differences*, not absolute hormone level. The willpower era didn't have hundreds of cal/day of biological pressure for willpower to fight; it had a floored signal under a starvation-mode regulatory regime.

**The 2013-14 inflection is leptin coming back online.** FM trajectory bottoms at 14 lb in June 2014. Binge rate progression: 0% through Jul-Aug 2013 → 13% (Oct) → 20-23% (Nov-Dec) → sustained 12-26% through 2014. AH#6 (floor effect) lands at SP minimum 16.1 lb FM on 2014-07-20; FM rises 16 → 35 over the next six months, binge rate ~9%. **The bottom and the eruption co-occur: as soon as FM crossed back above the leptin floor, the receptor — calibrated by years of high-FM exposure — read a large deficit and binges began.**

The TDEE story is what one would expect once the metabolic and eating channels are read together: **TDEE *fell* by 56 cal between bottom (1864 intake / 2194 TDEE / ratio 1.061) and regain phase (2359 intake / 2138 TDEE / ratio 1.065).** Intake outran TDEE by 551 cal/day. The metabolic channel (~13 d HL) did not get a chance to rebound because the eating channel ramped up immediately as the floor cleared, keeping fat mass climbing and the body in the "below SP" regime that never engages the metabolic boost.

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

## AA. Per-session lean addition channel: searched and not found (footnote)

A scan-only fit of "each workout adds Δ lb of lean, decaying exponentially" produced a flat objective surface with point estimate ~41 g/session and HL ~275 days, suggesting steady-state lean above baseline of +5-10 lb at 1-2× weekly training. CH tested this channel directly in the daily-resolution Kalman: at every (Δ, HL) combination, the channel either fails to improve any fit metric or reduces FM scan MAE marginally while worsening EB closure and innovation ACF. The candidate is below detection — composition scans through 2017-2026 hold LM at 138-150 lb without a clean +5 lb shift tied to workout frequency, so the chronic-accumulation mechanism implied by AA is not in the data. **The canonical filter does not include a per-session LM addition.** The lean-preservation-during-cuts effect is real but operates through Forbes-partition modulation (CB), not chronic stock-building.

Command: `python analysis/AA_lean_mass_training.py`, `python analysis/CH_strength_training_test.py`
Artifact: `analysis/AA_lean_mass_training.csv`

## CB. Workouts protect lean during cuts via Forbes-partition modulation

Using the Kalman's daily lean-mass state, two questions about strength training become well-posed: do workouts add lean per session, and do they preserve lean during deficit runs?

**Per-session lean response (regression of daily ΔLM on workout indicator with deficit + FM controls, by lag):**

| Lag (days after session) | Coefficient (g/session) |
|---|---|
| 0 (same day) | −4 |
| **1** | **+39** |
| 3 | +11 |
| 7 | −2 |
| 14 | −1 |

Lag-1 has a clean +39 g spike; every other lag is at the baseline cut rate or near zero. Same-day −4 g is sweat / glycogen depletion not fully absorbed by P1; the +39 g rebound is the day-after muscle protein synthesis response. The integrated weekly amplitude is ~46 g, but **the channel does not improve the canonical Kalman as a per-session LM addition** (CH searched (Δ, HL) jointly and found no improvement). The lag-1 spike is observed in dLM but is captured by the Forbes-partition modulation below — it does not need a separate per-session injection term to be modelled.

**Workouts preserve lean during cuts.** 37 cuts identified (≥14 days, rolling 7d intake < TDEE − 200 cal/day). For each cut, observed lean fraction of weight loss vs base-Forbes-predicted lean fraction at that FM:

| Workout density during cut | n runs | Lean % of dW | Forbes-predicted | Δ from Forbes |
|---|---|---|---|---|
| **None** (0 sessions) | 23 | +17% | +13% | **+4%** (more lean lost than Forbes) |
| Low (<1/wk) | 4 | 0% | +6% | −6% |
| **High (≥1/wk)** | 10 | −5% | +8% | **−13%** (less lean lost than Forbes) |

Regression `lean_frac_minus_forbes ~ workouts_per_wk + dW + FM`: **workouts/wk coefficient = −13.7 percentage points** [bootstrap 95% CI: −17.9, −6.6]. Each workout per week reduces the lean fraction of weight loss by ~14 percentage points beyond what plain Forbes predicts. The canonical Kalman incorporates this as `ff_modified = ff_base + α × workouts_per_wk` with α = 0.084 (see Kalman appendix); CH's parameter sweep favours α ≈ 0.12-0.18 with marginal additional improvement.

Concrete: at a typical post-2014 cut (FM 60-80 lb, Forbes-lean ~6%), zero workouts → lean is 10-15% of weight lost; ≥2 workouts/week → lean is roughly zero or slightly negative (LM grows during the cut). Most striking individual cuts:
- 2018-01-07 (19 d, 7 workouts during keto): −1.3 lb total, **+0.3 lb LM**.
- 2024-10-16 (37 d, 9 workouts, tirz era): −4.0 lb total, **+0.1 lb LM**.
- 2025-08-24 (54 d, 11 workouts, tirz era): −4.9 lb total, only −0.3 lb LM.

Pre-2014 willpower-era no-workout cuts *exceed* Forbes by 1.6-6.2 percentage points — consistent with starvation-mode lean catabolism that Forbes (free-living-calibrated) underestimates. Workouts cancel that excess and then some.

**Mechanism: trajectory interruption, not chronic stock-building.** Composition scans 2017-2026 hold LM at 138-150 lb without a clean +5 lb shift tied to workout frequency. There is no chronic LM stock to "build"; the protection during cuts is per-session partition modulation that interrupts the deficit's pull on lean. CC confirms this directly (workout protection lasts ~1 day, not weeks).

**Caveat.** Workout-rich cuts cluster post-2014; the regression controls for FM and dW but the era confound can't be fully broken.

Command: `python analysis/CB_workouts_lean.py`, `python analysis/CH_strength_training_test.py`
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

Per-workout heart-rate comes from the chart Samsung Health attaches to exercise sessions: `jsons/com.samsung.shealth.exercise/*/*.com.samsung.health.exercise.live_data.json` keyed by exercise_id. 448 usable walking sessions with live heart-rate, pace, and distance from 2016-04-18 to 2026-03-10. The 2016-2017 sessions sit on a different HR scale (walking HR median 153 in 2016 vs 134 in 2017) and aren't pooled with later watch-era data.

Restricting to the 2021+ watch era, pace-adjusted walking heart-rate:

`walk HR = 81.48 + 5.89 × kph`

Pace-adjusted yearly residuals: 2021 +3.76 bpm, 2022 +1.09, 2023 −1.83, 2024 +0.75, 2025 −1.29, 2026 −1.17 (sparse).

**Conclusion.** No deterioration in walking efficiency 2021+. If anything, 2023-2025 are slightly better than 2021-2022 at matched walking pace.

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

## AL. Walks do not rescue tirz-suppressed expenditure (null)

Tirz-era walk coefficient on TDEE/RMR is **indistinguishable from zero**. AD's calorimetry-anchored RMR boost (+14 cal/session/30d) is unaffected — walks still raise RMR per direct calorimetry — and AJ's tirz-suppression finding (calorimetry-based) holds. The specific claim that walks compensate for the drug-era metabolic suppression does not survive in Kalman state.

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

**Subject fit (1 free param + 1 decay param).** tonic = −69 cal/unit (intake-multiplier-scaled), K_AUC = 0.00187 / unit-week of cumulative exposure (AUC half-life ≈ 371 unit-weeks). Sawtooth prediction by era is preserved within 20% across the 80-week subject history.

**Why AUC, not calendar weeks.** Receptor desensitisation scales with cumulative agonist exposure, not elapsed calendar time. A calendar-time form gives every dose arm the same tachyphylaxis at any given week regardless of how much drug was present — predicting a wrong dose-response slope across SURMOUNT-1 arms (low-dose arms over-eroded). The AUC form correctly differentiates trajectories by cumulative occupancy and fits the dose-response slope.

**Trial validation (subject K_AUC carried over to all trials).** Best trial tonic re-derived under the canonical filter (INTAKE_MULTIPLIER=1.09, K_AUC = 0.00187): **−78 cal/unit**, SSE = 190.

**Subject vs trial tonic.** Subject −69 cal/unit (intake-multiplier-scaled), trial best −78. CE quantifies the dampening factor as `m ≈ 0.95` — the subject's central GLP-1 response is within 5% of the trial cohort baseline.

Command: `python analysis/BZ_drug_twocomponent.py`, `python analysis/AV_surmount_simulation.py`

## CE. Drug-response dampening: subject is nearly normal centrally

A single multiplier `m` on the central-pathway tonic characterises individual variation in apparent GLP-1 weight-loss response. Hold the BZ AUC-form drug model fixed and find the scalar `m` such that `tonic = m × tonic_NonDiabetic` best fits each cohort's published weight loss.

| Cohort | Dampening m |
|---|---|
| Non-diabetic (SURMOUNT-1) | 1.000 (calibration target) |
| T2D (SURMOUNT-2) | 0.69 (mechanistic note below) |
| **This subject (intake-fit, undercount-corrected)** | **0.95** |

**The subject's central GLP-1 response is within 5% of the non-diabetic trial cohort.** Renormalised against the SURMOUNT-1-perfect-fit reference, subject `m = 0.82` — still in the normal range. No structural dampening relative to the trial population.

**SURMOUNT-2 m = 0.69 is mechanistically misleading.** [Patoulias et al. 2025](https://diabetesjournals.org/care/article/48/11/e136/163253/Weight-Reduction-With-Tirzepatide-Varied) shows T2D's weight-loss attenuation is dominated by **glucosuria reduction**, not central-pathway scaling: at HbA1c <7% (low glucosuria), T2D weight loss matches non-diabetic SURMOUNT-1. A scalar on `tonic` numerically absorbs the glucosuria correction at the cohort level, but the underlying mechanism is peripheral. The subject is non-diabetic and has no glucosuria (CG); their `m = 0.95` is genuinely a central-pathway summary.

**Cross-modality consistency check.** Subject's `m = 0.95` predicts SURMOUNT-1-equivalent loss of −16.9% at 10 mg (proxy for ~12.5 mg average) and −19.7% at 15 mg under the trial simulator. Actual tirz-era body-weight loss: −14% over ~80 weeks. Simulator over-predicts by 3-6 pp, consistent with cohort vs subject body-composition differences (subject 27% body fat vs cohort 40%) and tirz-era activity factor.

**Generalisation.** Patient `m` could in principle be estimated from (a) intake response over the first few weeks of dosing, (b) early weight-loss trajectory vs trial cohort mean, (c) genotype + phenotype priors. The subject's near-normal `m` aligns with their typical clinical response to tirzepatide despite genetic risk variants (FTO, MC4R) documented in BACKGROUND.

Command: `python analysis/CE_drug_dampening.py`

## CF. Trial regain vs SmoothLatch: 3× magnitude gap, mechanism unresolved

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

**The m-scaling reconciliation is ruled out by CE.** A natural candidate fix would be to scale lipostat pressure by the same dampening factor `m` that the drug model uses (CE), so that population-level pressure is `m × 55 cal/lb`. Under the canonical filter, however, the subject's `m ≈ 0.95` is essentially trial-matching, so this fix can't separate the two populations. The structural cause of the 3× gap between the subject's predicted regain (~14%) and the trial cohort's published 60% remains unresolved; the candidate parameter most likely to differ is the latch tolerance or hold itself, since the subject's slow-gain decade kept FM continuously stable while trial loss is rapid by comparison.

**Speed-of-loss → asymptote correlation survives every parameterisation tested** (r = +0.6 to +0.95 across all latch/pressure choices). Faster loss → larger SP-FM gap at cessation → larger eventual regain. The published meta-regression treats this as an untested gap; our model says it's a real moderator that should explain per-trial variation in future analyses.

**Per-trial predicted regain at 52 wk** (subject latch + 55 cal/lb pressure carried over to each trial's loss trajectory):

| Trial | lead-in (wk) | loss rate (lb/wk) | predicted R(52) |
|---|---|---|---|
| SCALE liraglutide 3mg | 56 | 0.32 | 11% |
| STEP-1 sema 2.4mg | 68 | 0.51 | 9% |
| SURMOUNT-1 15mg | 72 | 0.72 | 10% |
| STEP-4 sema 2.4mg | 20 | 1.23 | 23% |
| SURMOUNT-4 15mg | 36 | 1.38 | 18% |

The fastest-loss trials (STEP-4, SURMOUNT-4) predict the largest regain at this subject's latch parameters — consistent with the published trend that SURMOUNT-4 sees the largest absolute regain in the meta-set.

Command: `python analysis/CF_regain_vs_meta.py`

## Asymptote of GLP-1 weight loss

Why does GLP-1 / GIP weight loss flatten? Trials almost always plateau by week 60-72 — the SURMOUNT-1 15mg curve runs ~−15% by week 24, ~−21% by week 48, ~−22.5% by week 72; STEP-1 sema 2.4mg runs ~−10% by week 20, ~−14.9% by week 68. The asymptote is dose-dependent and largely durable while the drug is continued. Mainstream framings ("metabolic adaptation", "the body fights back") name a phenomenon without naming a mechanism. This dataset's parameters give a quantitative one.

**Three timescales race.**

1. **Drug suppression** D(t) = tonic × blood(t) × exp(−K_AUC × ∫blood dτ). The drug provides exogenous satiety on a millisecond timescale, eroded over months by GLP-1R β-arrestin desensitisation and downstream hypothalamic adaptation as cumulative receptor exposure (AUC) grows. Subject K_AUC ≈ 0.00187 / unit-week; AUC half-life ≈ 371 unit-weeks (BZ).
2. **Endogenous eating pressure** P(t) = 55 × (SP(t) − FM(t)), symmetric. Instantaneous response of leptin / melanocortin signalling to current FM against the un-recalibrated receptor baseline. Positive below SP (drives intake up); negative above (suppresses intake).
3. **Set-point recalibration** dSP/dt = 0.01 × (FM − SP), gated on FM having held within ±3 lbs for 14 days (SmoothLatch, above). Receptor turnover only resets to a new baseline when chronically stimulated near that baseline — and chronic GLP-1/GIP agonism *bypasses* the endogenous receptor, so during active loss the latch is open and the SP doesn't move.

**The asymptote in one equation.** At equilibrium FM is stable, daily intake = TDEE, and pressure equals drug:

    55 × (SP − FM_eq) = D₀ × tachy(t_eq)

So at the asymptote, the gap that the drug holds open is

    SP − FM_eq = D₀ × tachy(t_eq) / 55 cal/lb

The equation gives the *residual gap held by the drug at equilibrium*, not the total weight loss. For SURMOUNT-1 15mg (tonic = −78 cal/unit trial fit, mean blood ~13 units at steady state): D₀ ≈ 1014 cal/day at week 0. Tachyphylaxis is AUC-based: at week 48, accumulated AUC ≈ 600 unit-weeks, tachy ≈ exp(−0.00187 × 600) ≈ 0.33, so D_eff ≈ 330 cal/day at the time the curve flattens. Residual gap = 330 / 55 ≈ 6 lb FM. Total FM loss is much larger — it's the cumulative deficit along the trajectory plus whatever the SmoothLatch has managed to drag the SP down during the stable phase near the asymptote. BZ's full forward integration produces SURMOUNT-1 15mg at −18.3% body weight (published −22.5%, Δ +4.2%) with the AUC-form tachyphylaxis carried over from the subject fit. The toy equation captures *why the curve flattens*; BZ captures *how the trajectory got there*.

**Why dose maps to asymptote depth.** D₀ scales with mean steady-state blood level, which scales with dose. Each doubling of dose roughly doubles D₀, which roughly doubles the gap held at equilibrium, which deepens FM_eq by the same multiple. SURMOUNT-1: −9% at 5mg, −16% at 10mg, −22.5% at 15mg. Sublinear because tachyphylaxis erodes the higher-dose advantage faster, and because the SmoothLatch latches fastest near the asymptote (when FM is most stable).

**Why duration on drug doesn't deepen the asymptote much past 60-72 weeks.** Tachyphylaxis erodes D(t) at the same rate the latch closes the SP-to-FM gap. The two cancel: the drug loses effective potency, but the SP also stops fighting back as it slowly recalibrates. The system hovers near equilibrium. SURMOUNT-1 gains <2% additional loss between week 48 and week 72.

**Why discontinuation produces large regain.** The latch only fired during the brief stable equilibrium near the asymptote; for the bulk of the loss, FM was changing fast enough to keep the latch open. So at discontinuation D drops to zero immediately, but SP is anchored not far below the pre-treatment level. Pressure = 55 × (SP − FM) reactivates at full strength against an un-recalibrated baseline. SURMOUNT-4 +14% in 52 weeks. STEP-1 ~+10%. The SmoothLatch's per-trial forward integration (CF) lands at +18.2% on SURMOUNT-4 and +8.9% on STEP-1 — directionally correct, with SURMOUNT-4 over-shot by ~4 points and STEP-1 within 1. The mean across the 6-trial meta-set is +14.2%, ~3× under the published 60% (the magnitude limitation discussed in the SmoothLatch section); the *qualitative* prediction that faster-loss trials regain more is preserved (rate→A: r = +0.86).

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

## CG. No detectable glucosuria at HbA1c ~5%

The SURMOUNT-2 HbA1c-stratified analysis ([Diabetes Care 2025](https://diabetesjournals.org/care/article/48/11/e136/163253/Weight-Reduction-With-Tirzepatide-Varied)) attributes T2D's attenuated tirzepatide weight loss largely to glucosuria — high HbA1c → urinary glucose → passive caloric drain that disappears as treatment lowers HbA1c. This subject is non-diabetic (HbA1c ~5.0%); the renal threshold is ~180 mg/dL. Could very large rapid-sugar meals still produce transient glucosuria?

**Initial test (biased).** For each threshold T (per-day or per-meal sugar in grams), subtract `max(0, sugar − T) × 4 cal/g` from logged daily intake and re-run the Kalman filter with the canonical `INTAKE_MULTIPLIER = 1.09`. The 1.09 multiplier was calibrated against the *unclipped* intake to close cumulative EB residual to ~0 (finding A); applying it to a clipped series mechanically under-counts intake by the discarded amount and drives EB drift that has nothing to do with glucose. The first version of this test (CG) drew its glucosuria-rejection conclusion in part from that artificial drift.

**Corrected test (DC).** Compensate the magnitude loss by rescaling the multiplier per clipping scheme: `m_clip = 1.09 × Σ(intake_original) / Σ(intake_clipped)`, so total consumed calories integrate to the same value as the canonical baseline. The Kalman comparison then reflects only the *timing* benefit of clipping — does the shape of the redistribution improve fit?

**Per-day sweep, compensated multiplier:**

| T_day | discarded (cal, 15 yr) | m_eff | FM MAE (lb) | innov ACF | EB closure (lb) |
|---|---|---|---|---|---|
| baseline (no clip) | 0 | 1.0900 | **1.081** | **+0.052** | **−0.04** |
| 200 | 11,676 | 1.0911 | 1.082 | +0.055 | −0.58 |
| 150 | 51,176 | 1.0948 | 1.080 | +0.055 | −1.40 |
| 125 | 104,208 | 1.0998 | 1.078 | +0.055 | −1.74 |
| 100 | 206,960 | 1.1097 | 1.074 | +0.057 | −1.82 |
| 75 | 395,376 | 1.1283 | 1.071 | +0.061 | −1.87 |
| 50 | 702,164 | 1.1599 | **1.068** | +0.066 | −1.76 |

Per-meal clipping (T_meal = 25-100 g) gives the same shape: monotonic small FM-MAE improvement, monotonic worsening of innovation ACF, monotonic worsening of cumulative EB closure.

**Verdict: still no glucosuria signal, now properly defended.**
- FM MAE drops at most 0.013 lb (1.2%), well within the identification-noise floor characterised by CJ (0.4% MAE shifts across plausible Atwater regimes).
- Innovation ACF *worsens* monotonically (+0.052 → +0.066) — residuals become more autocorrelated, the opposite of what a real physiological clip would produce.
- Cumulative EB closure drifts from a near-perfect −0.04 lb (baseline at 1.09) to −1.76 lb at the most aggressive clip, even with magnitude compensated. The drift is smaller than under uncompensated clipping (which reached −4.75 lb) but still wrong-direction.
- The "improvement" in FM-MAE scales monotonically with clipping aggressiveness; there is no bowl-shaped optimum at a physiological threshold. The model is just smoothing out a few high-variance days, not detecting a renal threshold.
- The d50 case discards 702,164 cal over 15 years (≈129 cal/day on average; 6% of consumed energy). That magnitude of urinary glucose loss is implausible for a non-diabetic subject.

**Conclusion: glucosuria is effectively ruled out at all tested thresholds.** Confidence is high (~95%) that meaningful (>50 cal/day average) glucose-in-urine caloric loss is not present in this subject. Sugar logging is not a meaningful uncertainty source — sugar is a subfield of carbohydrate, and the per-day macro Atwater closure (CI: median residual −19 cal, IQR ±60 cal) verifies that macro totals (and thus the database-derived sugar field that propagates from them) are accurate. The residual 5% uncertainty allows for: sub-hourly acute spikes invisible to daily/meal aggregates, and small (<50 cal/day average) steady effects below the test's detection floor.

**Implications.**
- The subject's CE dampening factor `m ≈ 0.95` (essentially trial-matching) is not a glucosuria artifact; the central-pathway reading from CE is what it is, and it isn't being inflated upward by hidden urinary calorie loss.
- Whatever residual gap there is between the subject and the trial cohort cannot be explained by the SURMOUNT-2 peripheral (glucosuria) mechanism. The subject is non-diabetic, urinary calorie loss is ruled out, and any difference operates centrally.
- The 2013 willpower era's −3.45 lb EB residual is NOT explained by glucosuria either (low-cal era, sugar threshold rarely exceeded). The conservation-mode reading of AC + BK (TDEE drops further than the free-drift Kalman naturally finds at the leptin floor) remains the live explanation for that residual.

**Note on the residual-regression false positive.** A naive regression of `(observed dW − predicted dW)` on excess sugar produced a strongly negative coefficient (−0.005 to −0.010 lb/g of excess sugar) even after controlling for total calories, day-of-week, and FM. The magnitude was 4-8× the glucosuria physical upper bound (`-4 cal/g ÷ 3140 cal/lb = -0.00127 lb/g`), indicating the signal was gut-content / non-electrolyte water inflation that the P1 glycogen+sodium correction does not capture but the Kalman BUF state does. **The clip-and-rerun-with-compensation protocol is the proper test**; the residual regression confounds energy loss with measurement noise.

Command: `python analysis/DC_glucosuria_compensated.py` (compensated retest; supersedes `analysis/CG_sugar_glucosuria.py`)
Artifact: `analysis/DC_glucosuria_compensated.csv`

## CI. Textbook Atwater factors are right; the EB drift is the intake undercount

Atwater factors (4 cal/g protein, 4 cal/g carb, 9 cal/g fat) are the textbook conversion from macros to calories. Could this subject's tissue absorb at different per-gram rates — protein closer to 3 cal/g (as the Atwater 4 number includes incomplete amino-acid utilisation), or fat closer to 10 cal/g for high saturated-fat days?

**Method.** Per-day macro Atwater closure is tight (median residual −19 cal, IQR ±60 cal) — clean enough to test. Replace logged calories with `alt = logged + (a−4)P + (b−4)C + (c−9)F`, preserving the residual filler (alcohol, sugar alcohols) untouched. Re-run the Kalman, sweep (a, b, c).

**1D sweeps under canonical filter:**

| Variant | Σ adjustment (cal) | FM MAE | innov ACF | EB (lb) |
|---|---|---|---|---|
| Textbook (4, 4, 9) | 0 | **1.084** | **+0.049** | −2.22 |
| (3.0, 4, 9) — protein at 3 | −403k | 1.087 | +0.049 | −2.37 |
| (5.0, 4, 9) — protein at 5 | +403k | 1.081 | +0.051 | −2.12 |
| (4, 3.0, 9) — carbs at 3 | −1,303k | 1.085 | +0.057 | −4.23 |
| (4, 5.0, 9) — carbs at 5 | +1,303k | 1.078 | +0.058 | +0.46 |
| (4, 4, 7.5) — fat at 7.5 | −796k | 1.086 | +0.044 | −3.40 |
| (4, 4, 11) — fat at 11 | +1,062k | 1.082 | +0.059 | −0.45 |

The pattern is monotonic: shifting any factor up tightens cumulative EB closure (because it adds calories) and shifting any down loosens it. **The shifts that close EB best add ~1M cal cumulatively, which is almost exactly the documented intake undercount.**

**Confirming with a uniform multiplier.** A flat `intake × 1.10` (textbook Atwater throughout) adds +1,166k cal cumulatively and produces FM_MAE 1.080, EB +0.23 — virtually identical to Atwater (4.4, 4.4, 9.9) which adds +1,160k cal and produces FM_MAE 1.080, EB −0.27. The two formulations are mathematically equivalent (4.4 = 4 × 1.10, 9.9 = 9 × 1.10). **Any apparent improvement from tweaking Atwater factors is the Kalman discovering the undercount through the macros, not a metabolisable-energy difference.**

**Verdict.** Textbook (4, 4, 9) is the correct per-gram conversion. The 8-14% intake undercount documented in finding A explains the cumulative EB drift, and applying it as a uniform multiplier on intake is the right place to absorb it.

**Fine sweep on uniform intake multiplier** (textbook Atwater, intake × `m`):

| `m` | FM MAE | innov ACF | EB (lb) | TDEE/RMR mean |
|---|---|---|---|---|
| 1.00 | 1.084 | +0.049 | **−2.22** | 1.185 |
| 1.04 | 1.083 | +0.051 | −1.34 | 1.214 |
| 1.07 | 1.082 | +0.051 | −0.60 | 1.241 |
| **1.09** | **1.081** | +0.052 | **−0.04** | **1.265** |
| 1.12 | 1.080 | +0.053 | +0.77 | 1.281 |

Composite minimum lands at `m = 1.09`, matching finding A's mean undercount estimate. **`INTAKE_MULTIPLIER = 1.09` is the canonical-filter default.** Cumulative EB closes within ±1 lb over 15 years; FM scan MAE marginally tighter; innovation ACF essentially unchanged.

**TDEE/RMR mean is 1.28** under the multiplier — reasonable for a subject who does up to 33 walks/month in walk-rich periods, 1-3 strength sessions/week, and historical hospital-shift cycling. The high-end calorimetry-day TDEE/RMR ratios (1.39-1.47) cluster in walk-rich periods (May 2022, on-tirz with walks 2024-09) where elevated activity factors are physiologically expected.

**The undercount is real, the multiplier is the right knob.** Building it into the canonical filter at 1.09 makes the model internally consistent (EB closes near zero) and produces physiologically defensible TDEE estimates. Atwater per se is unchanged: innovation ACF minimises at textbook (+0.049), confirming the *relative* per-gram values are right; only the absolute intake level is biased low, and that bias goes into the multiplier.

A phase-conditional version (1.125 gain / 1.105 loss / 1.102 stable per finding A) was also tested: no improvement over the uniform multiplier on composite score.

Command: `python analysis/CI_atwater_test.py`

## CJ. Fine Atwater sweep — textbook is consistent within identification noise

Refines CI with 0.05 cal/g resolution and per-config EB-closing multiplier. Joint Nelder-Mead optimum lands at biologically implausible bounds (P=5.50, C=4.43, F=7.00 with multiplier 1.117), giving FM_MAE = 1.0766 — only 0.0041 lb (0.4%) better than textbook. Two of three macros at the box edges and the optimum sliding along a degeneracy.

The four parameters (a_protein, b_carb, c_fat, m) span only ~2 effective directions. Textbook + 1.09 multiplier sits in a flat basin of nearly-equivalent solutions. **Textbook Atwater is consistent with the data within identification noise; the 9% gap is all multiplier (uniform per-calorie).** No carb-specific signal that would suggest glucosuria, no protein-specific signal that would suggest TEF asymmetry. The undercount is uniform.

Command: `python analysis/CJ_atwater_fine.py`
Artifact: `analysis/CJ_atwater_fine.csv`, `analysis/CJ_atwater_fine.png`

## CK. Set-point-mediated clipping: registered negative

Tested whether storage efficiency η = predicted_ΔFM − actual_ΔFM has a regime change at some FM threshold (revealed-SP), as the metabolic-clipping-above-SP hypothesis would predict. Per-90d window fits with free SP* per window: 4 of 5 preregistered checks fail. β > 0 in only 45% of windows (predicted >70%); SP* hits the window's FM range bound 54% of the time (predicted few); per-cycle SP* SD = 67 lb (predicted < 10 lb); corr(SP*, EMA-SP-50d) = +0.976 (predicted 0.4-0.9 — at +0.97 we just rediscovered EMA, not a real threshold).

**No metabolic regime change at any FM threshold.** Either clipping is continuous-not-thresholded (consistent with AM-style continuous pressure), or this signal is too noisy to detect a threshold. Storage-efficiency analysis cannot be the path to a hidden SP for this dataset.

Command: `python analysis/CK_revealed_set_point.py`
Artifact: `analysis/CK_revealed_sp.csv`, `analysis/CK_revealed_sp.png`

## CN. No fixed defended weight; SP-distance partial r ≈ −0.18 after AR

Two diagnostics on whether a defended SP is necessary:

1. **AR(intake) + FM + sp_distance regression**: weekly intake autoregression (lag 7/14/28/56) gets R²=0.192. Adding FM lifts to 0.195. Adding sp_distance_60d lifts to 0.221 (ΔR² = +0.025). Marginal r(intake, sp_dist) = −0.35 collapses to partial r = −0.18 after AR + FM controls. **About half of AG's apparent SP signal was intake autocorrelation; about half is genuine residual.** The headline AG marginal r = −0.62 is misleading.

2. **Maintenance-FM attractor test**: 38 windows (90-day mean surplus < 50 cal/day). Per decade: 2011-14 maintenance FM 15-36 lb (median 30.8), 2015-19 maintenance 36-66 lb (median 55.7), 2020-24 maintenance 63-85 lb (median 72.0). Maintenance IQR = 19 lb across a lifetime range of 87 lb. **The body settles at whatever FM the lifestyle of the decade puts it; there is no single defended weight.** Slow drift over years is real but observationally identical to lifestyle drift in this dataset.

Command: `python analysis/CN_setpoint_necessity.py`
Artifact: `analysis/CN_setpoint_necessity.png`

## CW. No constant willpower budget under simple drain functions

Across 19 falls (peak_a → trough), tested whether cumulative willpower drain at the trough is constant under various simple drain rate functions. CV of drain at trough across cycles:

| Drain function | CV |
|---|---|
| (SP_i − BW), SP free per cycle | 0.081 (force-fit) |
| FM × deficit | 0.712 |
| FM × time | 0.765 |
| BW × deficit | 1.069 |
| deficit alone | 1.209 |
| time alone | 1.228 |

No simple drain function gives consistent total drain across cycles. The "constant willpower budget that exhausts at troughs" hypothesis fails. The SP-free fits get CV=0.08 only by force-fitting per-cycle SP_i, which then track mean BW per cycle (i.e. just rediscover the body's running average).

Command: `python analysis/CW_willpower_exhaustion.py`
Artifact: `analysis/CW_willpower_exhaustion.png`

## CY/CZ. Kalman attribution bias and tuning verdict

CY tested whether tightening TDEE drift (or anchoring it deterministically to RMR + steps + workouts) changes K hysteresis. Under deterministic TDEE = c·RMR + d·steps, K direction asymmetry collapses (10/19 cycles fall>rise vs 14/19 baseline; median Δ TDEE/RMR +0.001 vs +0.032). **Most of K's directional signal is Kalman bookkeeping** given the subject's asymmetric rise/fall rates.

CZ2/CZ3 swept σ_TDEE with co-fit intake multiplier and tested step-anchored variants. Tightening σ_TDEE forces intake multiplier above 1.20 (implausible) and pushes calorimetry TDEE/RMR mean above 1.40 (implausibly high activity factor). Step-anchored filters track FM ~0.8 lb worse than baseline at every σ. **The canonical filter (σ_TDEE=50, intake_multiplier=1.09) sits at the calibration sweet spot.** Tighter is worse, not better. The directional bias is a downstream-analysis problem, not a filter-tuning problem.

Command: `python analysis/CY_constrained_tdee.py`, `python analysis/CZ2_kalman_tuning.py`, `python analysis/CZ3_kalman_steps.py`

# Thermoregulation

## DB. Protein → cold behavior, behavioral test from utility water meter

Lamming-Lab framing: incomplete dietary-isoleucine processing depresses body temperature, so higher protein intake should be followed by colder feelings and warming behavior. A new behavioral endpoint comes from the household water meter (KUB AMI, hourly, sourced from `~/Finance/Statements/KUB/water_house_hourly.csv`) processed by `baths/baths.py` into 447 detected bath events from 2020-09 onward. Detector calibration: Nov 2025 (subject was in California; co-resident does not bathe) returns 0 events; the winter-to-summer bath ratio is 2.5:1.

In winter (Dec/Jan/Feb, where baths concentrate at 152 of 440 events), prior-days protein intake significantly predicts a bath that night:

| Predictor (winter) | r (point-biserial) | p | Δ protein, bath vs no-bath day |
|---|---|---|---|
| protein same day | +0.036 | 0.40 | +1.9 g |
| protein lag-1 day | +0.109 | **0.011** | +5.7 g (95 % CI +1.2 / +10.1) |
| protein lag-2 days | +0.167 | **<0.001** | +8.7 g (95 % CI +4.2 / +13.4) |
| protein 3-day trailing mean | +0.176 | **<0.001** | +7.5 g (95 % CI +3.7 / +11.2) |

Direction matches the hypothesis. Summer flips sign on small numbers (n=60 baths) and probably reflects digestive thermogenesis from a same-day high-protein meal rather than the lagged metabolic process the hypothesis names. Spring and fall are flat.

The direct body-temperature endpoint (n=355 circadian-corrected daily baselines since Dec 2023, from `AS_temperature_baseline.csv`) goes the predicted direction at lag-1, lag-2, and 3-day mean (all negative slopes: −0.005 / −0.013 / −0.020 °F per +10 g protein) but no individual lag reaches significance — directionally consistent, underpowered. Same-day temperature is mildly positive, again consistent with digestive thermogenesis dominating short-horizon thermal effects.

Behavioral data is the more powerful test for now; body temp will tighten as more daily readings accumulate. Both endpoints point the same way.

Command: `python analysis/DB_bath_protein.py`
Artifact: `analysis/DB_bath_protein_summary.csv`, `analysis/DB_temp_protein_summary.csv`

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
- LM rate caps: gain ≤ +0.15 lb/day; loss ≤ -0.5 lb/day (extreme protein turnover, rarely hit).
- BUF: HL = 5 days, hard cap ±3 lb, daily noise std ≈ 0.57 lb. Mean ≈ 0, max excursion in 15 years ≈ 2.2 lb (well under cap).
- TDEE: σ ≈ 50 cal/day, hard daily-step cap ±100 cal/day. Free random walk (no mean reversion, no anchor).
- Composition scans observe FM and LM directly with σ ≈ 1.5 lb each.
- Calorimetry RMR observations are **held out** — never used by the filter; they're a posterior validation channel.

**Intake undercount correction (CI, finding A).** Logged intake is biased low by ~9% on this subject (composition-aware EB check: 8-14% by phase). The canonical filter multiplies logged intake by 1.09 before consumption. This closes cumulative 15-year EB residual to ~0 lb (vs −2.2 lb without). TDEE estimates are scaled to actual expenditure rather than "logged-cal expenditure."

**Diagnostics on current fit.** Innovation lag-1 ACF = 0.05 (close to white). FM scan MAE = 1.08 lb. LM scan MAE = 1.78 lb. TDEE/RMR mean across 25 calorimetry measurements = 1.28 (range 1.10-1.47), reflecting the subject's actual activity level (walks, strength, cycling). r(TDEE, baseline temperature) partial⏐FM = 0.09. Cumulative energy-balance closure = ~0 lb over 15 years.

**Tuning verdict (CZ2, CZ3).** Tighter σ_TDEE forces intake multiplier above 1.20 (implausible vs CI's 1.09 snacking-undercount estimate) and pushes calorimetry TDEE/RMR mean to >1.40 (implausibly high activity factor; in_band drops 96% → 20%). Step-anchored variants (TDEE = c·RMR + d·steps + small drift) degrade FM scan MAE by ~0.8 lb at every σ tested. The canonical σ_TDEE=50 / m=1.09 sits at the calibration sweet spot.

**Attribution-bias caveat (CY).** Free-drift TDEE is constrained only by mass balance over windows: `TDEE ≈ intake − ΔFM·3500 − ΔLM·1100 − ΔBUF`. When ΔFM is negative, the filter mechanically attributes upward TDEE; when ΔFM is positive, downward. This produces a built-in directional asymmetry of ~500 cal/day given the subject's asymmetric rise/fall rates. The asymmetry **averages out across calorimetry days** (which mix directions), so calorimetry validation passes — but **direction-conditional analyses (K, AI, B's TDEE drop, AJ, AL) inherit it**. Findings that depend on `TDEE during falling` vs `TDEE during rising` at matched FM should be considered upper-bounded; the *true* direction-dependent metabolic effect for this subject is at most a fraction of the filter's reported magnitude. Cleanly resolving this requires direct calorimetry stratified by direction at matched FM, which we do not have at sufficient density.

Command: `python analysis/P4_kalman_filter.py`
Artifact: `analysis/P4_kalman_daily.csv`, `analysis/P4_kalman_sweep.csv`
