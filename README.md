# 15 years of daily food intake, weight, body composition, and metabolic data

Fifteen years ago, I lost 100 lbs [in a fairly boring manner](https://www.reddit.com/r/loseit/comments/498afc/) by eating 1200 calories a day. It's creeping back. There are no low-hanging fruits — no sugarwater, no fast food, near-zero restaurant meals. Since then I've recorded every bite eaten every day, whether 10 calories of broth on a weekend fast or a 4000 calorie binge, using a food scale for everything. 54,925 food items across 5,429 days. Zero gaps.

This is not self-reported-from-memory dietary recall. The cumulative energy balance — intake minus expenditure, accumulated daily — closes to ±5 lbs over the full 15 years. A composition-aware model fitted to 25 indirect calorimetry measurements and 70 body composition scans quantifies the undercount at ~10-14%, uniform across weight loss, gain, and stable phases. **The best metabolic ward studies achieve ~5% error over 2-4 weeks. Self-reported dietary recalls average 30-50% error with missing days. This dataset sustains ~10-14% error across 5,429 consecutive days.**

**The dataset contains hundreds of overlapping natural experiments**: months of keto, weekend 36-hour fasts, low protein, high fiber, potatoes-only, waves of monotonous meals, daily ice cream, different cooking oils. Combined with simultaneous weight, body composition (BOD POD + InBody), resting metabolic rate (Cosmed indirect calorimetry), steps, sleep, body temperature, blood fatty acid panels, and 80 weekly tirzepatide injections logged from the first dose. **If the answer to obesity requires a complex overlap or sequence of conditions, it may be hidden within and first discovered through data mining rather than invented by a brilliant hypothesizer.** Very likely, the answer is already in the data — a few weeks spread across years where something worked, masked by the noise in daily weigh-ins.

I am [mostly](https://ptable.com) [retired](https://flightaware.com), 43, no prescriptions, no stressors, no emotional eating. No meat ([35 years](https://www.reddit.com/r/vegetarian/comments/35d7iq/)). No workplace contaminant exposure. Reverse-osmosis water. The worst-case common-variant obesity genotype: FTO homozygous risk at all five loci, MC4R heterozygous, UCP2 reduced thermogenesis — not [easily-solved obesity](https://doi.org/10.2105/AJPH.2015.302773). See [BACKGROUND.md](BACKGROUND.md) for the full genetic profile and health history.

## What the analysis has found

**The body defends a narrow expenditure band regardless of weight.** At 165 lbs (20 lbs fat, 2013) and 225 lbs (82 lbs fat, 2024), derived TDEE is ~2100-2200 cal/day. The TDEE/RMR ratio drops to 1.02 during sustained restriction and recovers to 1.14 when restriction eases. **My set point moves up or down only when some conditions are met**, and weight loss above it is trivial; below it, [fiendishly difficult](https://doi.org/10.2105/AJPH.2015.302773).

**A defended fat mass set point is identified in two places: published GLP-1/GIP withdrawal trials, and this subject's own tirzepatide era.** In the trials, 60-70% of the lost weight returns within a year of stopping. In this subject's 18 months on tirzepatide, the slow rise of weekday-adjusted intake at constant dose is fit best by eating pressure proportional to fat mass lost below the *pre-drug* level — a frozen set point at **30-40 cal/day per lb** — and that form, not receptor tachyphylaxis, is the one that predicts the 12.5 mg era out of sample ([DF](FINDINGS.md#df-what-erodes-the-drug-effect-frozen-set-point-pressure-not-tachyphylaxis)). The SmoothLatch rule (the set point adapts only while FM holds within a narrow band) survives, but with the *trial-fit* parameters (tol ≈ 1 lb, hold ≈ 30 d, rate ≈ 0.0025/day): under those the set point barely moves during drug-driven loss, which is what both the trials and this subject's drug era show. The earlier subject-fit latch (±3 lb, 14 d, 1%/day) lets the set point track fat mass down during the drug era and is rejected by it. **The earlier claim that the subject-fit rule reproduces SURMOUNT-4 regain at +14.1% rested on a mis-scaled body-weight denominator in BD; correctly scaled, that lax latch predicts +5.5% of cessation weight against the published +14.0%.**

**The 45-50 day timescale does not survive forward-looking tests on this subject ([DD](FINDINGS.md#dd-forward-window-test-the-subject-side-set-point-signal-is-mostly-reverse-causation)).** Three of the four "independent" identifications were not independent of intake persistence. The binge-rate fit (AG) correlated a *trailing* 90-day binge rate with today's EMA distance, so past binges mechanically put fat mass above its EMA: r = −0.64 trailing becomes −0.12 with a forward window, and the residual sign is anti-lipostat. The weigh-in-only test (AU) predicts *lower* future intake below the EMA, the persistence direction. The intake-free FM fit (P3) sits on the same surplus-regression ridge as AM (BR). The lead-lag correlation of EMA distance with future intake sits inside an intake-only AR(60) null at every lead from 1 to 365 days. What remains is small: adding EMA distance to a past-intake forecast of the next 90 days' intake improves RMSE by 3%, with a coefficient of about +12 cal/day per lb below the EMA (8 of 9 yearly folds positive) — a fifth of the 55 cal/lb claimed, and absent at 30 days. The pre-drug decade mostly cannot identify a set point because fat mass never left its own trailing average; the drug era can, and does (above).

**Two channels, different speeds.** The body closes the gap on two fronts. An eating channel (~50-day timescale, ~55 cal/lb) shifts daily intake toward the set point. A metabolic rate channel (~9-day timescale, confirmed by 25 calorimetry measurements) burns ~90 cal/day more during weight loss than weight gain at the same body composition — actively assisting loss when fat mass is above the set point, passively permitting regain when below. The eating channel is stronger but slower. Binge size is constant regardless of distance from the set point; what shifts is the *frequency* and the depth of deficit days. The set point tilts the entire daily distribution rather than triggering discrete events.

**Receptor-dynamics reading.** The latch's stability requirement and the 45-50 day half-life are what hypothalamic leptin/melanocortin receptor recalibration would look like at a behavioral level: receptor populations adjust to chronic exposure, and the chronic exposure has to actually persist for the system to adapt. During drug-driven loss the satiety circuit is being stimulated by exogenous GLP-1/GIP — the *endogenous* receptor isn't being challenged at the new low FM, so it doesn't recalibrate. The set point stays anchored to pre-treatment, which is exactly what produces SURMOUNT-4-style regain on discontinuation. The same frame explains why GLP-1 trials asymptote where they do (see [FINDINGS.md § Asymptote](FINDINGS.md#asymptote-of-glp-1-weight-loss)) and why this subject's 2011-13 willpower era didn't produce binges until FM crossed back above ~18 lbs (the leptin-floor reading of [AC + BK](FINDINGS.md#ac--bk-the-2013-inflection-and-the-leptin-floor)).

**Tirzepatide suppresses appetite at about −90 cal per unit of daily-mean blood level** (AX, identified from within-week injection cycle variation with the timestamp-aware exposure series; −80 after removing the pre-drug weekday pattern, since every injection is a Tuesday). The weekly injection cycle: the day after injection (peak exposure, Tmax ≈ 24 h) averages 1595 cal, day 6 (trough) averages 2219 cal, a ~620 cal/day swing driven by pharmacokinetics. The slow rise of intake across the drug era, previously read as cumulative-exposure tachyphylaxis, is fit better (ΔAIC −38) by pressure proportional to fat mass lost below the pre-drug level, and only that form predicts the 12.5 mg era out of sample (DF). Binge-to-binge escalation drops from 31.6% to 0%.

**The drug also suppresses the body's metabolic cooperation with weight loss.** Direct calorimetry shows RMR 206 cal below composition-predicted on the drug — the metabolic boost that normally accelerates fat loss is pharmacologically eliminated. This effect operates on longer timescales than the injection cycle and cannot be separated from the appetite effect within-week. Weight loss is 93% fat; the net energy budget is approximately -450 cal intake reduction, +200 cal metabolic cooperation lost, net ~-250 cal/day deficit.

**One gap-form model, fit once on the drug era, run against everything ([DG](FINDINGS.md#dg-one-gap-form-model-fit-once-run-against-everything)).** `intake = I₀ + s·B + P·(SP − FM)` with s = −80 cal/unit, P = +43 cal/day per lb (true calories), SP frozen at the pre-drug fat mass by the trial-fit latch. Fed only the exposure series, it reproduces the subject's own 18-month fat-mass path to 1.2 lb RMSE (every dose era within ~1 lb) and, with no cohort parameters, the *shape* of SURMOUNT-4's regain curve (+8/+11/+12/+12% vs published +7/+11/+13/+14% at 12/24/36/52 weeks off). What it does not do cleanly: the published dose-response is concave (5 mg gets 71% of the 15 mg loss at a third of the exposure). Cohort heterogeneity cannot produce that from linear responders (DI: 5/15 ratio ≤ 0.45); pharmacological saturation with EC50 ≈ 12-20 exposure units can, and the subject's own drug era mildly prefers that EC50 to linear while excluding the cohort-exact EC50 ≈ 5. Published withdrawal regain, meanwhile, is 80-90% explained by intake simply returning to habit with TDEE lower at the lower weight (P = 0 regains 55% of the loss in a year; published 65-80%); the trials need almost no active defense, and none of them depends on how fast the weight came off.

**Lean mass, not fat mass, carries the forward appetite signal in the pre-drug decade ([DH](FINDINGS.md#dh-set-point-dynamics-zoo-does-the-set-point-move-down-and-what-moved-in-2014)).** Across 19 set point rules scored on the next 90 days, no fat-mass gap predicts future intake beyond past intake (partial r ≤ +0.1), but lean mass below its running maximum at a given fat mass predicts higher intake over the next 90 days at partial r +0.4 to +0.6, with scan-only lean mass, year fixed effects, workouts and TDEE controlled. Lean mass fell 157 → 142 lb in 2011-14 and never recovered. This is what Dulloo's protein-stat predicts and what a cumulative-restriction history would also produce; the two are being separated with training-induced lean mass (DJ). The 2014 binge eruption is not predicted by any model's pressure: binges switched on in 2014Q3-Q4 after two years in which a frozen set point would have exerted 2,850 cal/day with no leak.

**Why GLP-1 weight loss asymptotes** is a question the literature usually answers with "metabolic adaptation" without specifying a mechanism. This dataset gives a simpler one than the earlier tachyphylaxis story: the plateau is the point where drug suppression `s × B` balances eating pressure `P × (SP₀ − FM)` against a set point that has not moved off its pre-drug value. In this subject, P ≈ 30-40 cal/day per lb and the plateau at FM ≈ 61 lb sits 25 lb below the pre-drug level, i.e. ~800 cal/day of pressure against ~1000 cal/day of suppression at 12.5 mg. No receptor erosion is needed; during the 15-week FM plateau intake *fell* 185 cal/day while cumulative exposure kept rising. See [FINDINGS.md § Asymptote of GLP-1 weight loss](FINDINGS.md#asymptote-of-glp-1-weight-loss).

**Not all restriction is equal.** Long runs (≥6 days) and low-carb restriction recover best. Low-protein and high-step restriction produce the largest metabolic penalties. Potato diets show zero binges across 69 days and high TDEE/RMR ratio, but severe rebound after stopping.

**What didn't survive.** The model above is what's left after testing every other candidate. Each elimination has a number and a reproducer in FINDINGS:
- [Gravitostat](https://doi.org/10.1073/pnas.1800033115) (foot-pounds drives intake): r = +0.05, wrong direction.
- [Speakman's dual-intervention zone of indifference](https://doi.org/10.1242/dmm.008698): no zone, continuous sigmoid.
- Set point tracking scale weight: dominated by FM tracking at every HL (-0.60 fat vs -0.52 scale).
- Fixed set point: r = +0.20 (and wrong sign) vs r = -0.60 for the moving 50d EMA.
- [Lowe et al. (2007)](https://doi.org/10.1002/eat.20405) fixed weight-suppression proxy: directional finding replicates, but their fixed reference is dominated by a moving 50d EMA.
- Symmetric EMA at any HL: undershoots SURMOUNT-4 regain by 5-10×.
- Subject-fit SmoothLatch (±3 lb, 14 d, 1%/day): tracks FM down during the drug era (DF), predicts +5.5% SURMOUNT-4 regain vs +14.0% (BD, corrected denominator).
- Subject-side "confirmations" of the 50-day timescale (AG binge rate, AU weigh-in test, CN partial r): trailing-window reverse causation and intake persistence (DD); forward residual ≈ +12 cal/day per lb at 90 d only.
- AUC-tachyphylaxis as the drug-era erosion mechanism: dominated by a frozen-set-point gap (ΔAIC −38) and over-predicts 12.5 mg-era intake by 340 cal/day out of sample (DF).
- SmoothLatch speed-of-loss prediction (faster loss → more regain): all four published withdrawal arms regain about two-thirds of the loss in a year regardless of lead-in speed (SURMOUNT-4 67%, STEP-1 65%, STEP-4 66%, SCALE ~80%). Intake returning to habit explains 55% with no set point; P ≈ 0-5 cal/day/lb and a set point that never moves on-drug explain the rest (DG).
- Cohort heterogeneity as the source of the concave dose-response: a mean of linear responders is linear; even an intake floor pinning 11% of the 15 mg arm reaches only a 5/15 ratio of 0.45 vs 0.71 (DI).
- Flatt's glycogenostat (depleted glycogen → eat more): partial r +0.09, wrong sign (DH).
- Asymmetric ratchet: EMA artifact of the willpower era; doesn't replicate on natural-dynamics data.
- One-stock willpower reservoir (BK): 25.6% impossible days during 2011-12.
- Yo-yo / variance-as-metabolic-damage: sign-reversed; cf. [Biggest Loser reanalysis](https://doi.org/10.1002/oby.23308).
- [Pontzer's constrained TEE](https://doi.org/10.1016/j.cub.2015.12.046): walk *sessions* raise *resting* metabolic rate (AD); not just TEE redistribution.
- "Metabolic adaptation" as the GLP-1 plateau: names a phenomenon, not a mechanism. The plateau is the three-timescale equilibrium between drug suppression, gap pressure, and the SmoothLatch — see [FINDINGS.md § Asymptote](FINDINGS.md#asymptote-of-glp-1-weight-loss).

Steps don't predict weight change at any timescale. Sleep hours, protein, fat, and fiber have no independent signal after controlling for calories, carbs, and sodium. The week-scale intake invariance claim was falsified.

Every finding with numbers is reproduced by a standalone script in [FINDINGS.md](FINDINGS.md). Omega-6:3 ratios, circadian misalignment, and other untested theories are also documented there.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Data

| File | Rows | Description |
|---|---|---|
| `intake/intake_foods.csv` | 54,925 | Every food item, 8 nutrients per item |
| `intake/intake_daily.csv` | 5,429 | Daily nutrient totals |
| `weight/weight.csv` | 1,693 | Daily fasted weight |
| `steps-sleep/steps.csv` | 4,275 | Daily step counts |
| `steps-sleep/sleep.csv` | 2,057 | Sleep periods |
| `composition/composition.csv` | 70 | Body composition (FM, FFM, segmental) |
| `RMR/rmr.csv` | 25 | Indirect calorimetry RMR |
| `drugs/tirzepatide.csv` | 666 | Daily-mean PK exposure (timestamp-aware) + tachyphylaxis |

Raw data in `intake/`, `weight/`, `composition/`, `RMR/`, `steps-sleep/`, `drugs/`, `temperature/`, `fatty-acids/`, `workout/`, `travel/`. All extractors idempotent. Pipeline in [CLAUDE.md](CLAUDE.md).

I accept any question, no matter how personal, in service of our goal. You will not offend me.
