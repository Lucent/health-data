# 15 years of daily food intake, weight, body composition, and metabolic data

Fifteen years ago, I lost 100 lbs [in a fairly boring manner](https://www.reddit.com/r/loseit/comments/498afc/) by eating 1200 calories a day. It's creeping back. There are no low-hanging fruits — no sugarwater, no fast food, near-zero restaurant meals. Since then I've recorded every bite eaten every day, whether 10 calories of broth on a weekend fast or a 4000 calorie binge, using a food scale for everything. 54,925 food items across 5,429 days. Zero gaps.

This is not self-reported-from-memory dietary recall. The cumulative energy balance — intake minus expenditure, accumulated daily — closes to ±5 lbs over the full 15 years. A composition-aware model fitted to 25 indirect calorimetry measurements and 70 body composition scans quantifies the undercount at ~10-14%, uniform across weight loss, gain, and stable phases. **The best metabolic ward studies achieve ~5% error over 2-4 weeks. Self-reported dietary recalls average 30-50% error with missing days. This dataset sustains ~10-14% error across 5,429 consecutive days.**

**The dataset contains hundreds of overlapping natural experiments**: months of keto, weekend 36-hour fasts, low protein, high fiber, potatoes-only, waves of monotonous meals, daily ice cream, different cooking oils. Combined with simultaneous weight, body composition (BOD POD + InBody), resting metabolic rate (Cosmed indirect calorimetry), steps, sleep, body temperature, blood fatty acid panels, and 80 weekly tirzepatide injections logged from the first dose. **If the answer to obesity requires a complex overlap or sequence of conditions, it may be hidden within and first discovered through data mining rather than invented by a brilliant hypothesizer.** Very likely, the answer is already in the data — a few weeks spread across years where something worked, masked by the noise in daily weigh-ins.

I am [mostly](https://ptable.com) [retired](https://flightaware.com), 43, no prescriptions, no stressors, no emotional eating. No meat ([35 years](https://www.reddit.com/r/vegetarian/comments/35d7iq/)). No workplace contaminant exposure. Reverse-osmosis water. The worst-case common-variant obesity genotype: FTO homozygous risk at all five loci, MC4R heterozygous, UCP2 reduced thermogenesis — not [easily-solved obesity](https://doi.org/10.2105/AJPH.2015.302773). See [BACKGROUND.md](BACKGROUND.md) for the full genetic profile and health history.

## What the analysis has found

**The body defends a narrow expenditure band regardless of weight.** At 165 lbs (20 lbs fat, 2013) and 225 lbs (82 lbs fat, 2024), derived TDEE is ~2100-2200 cal/day. The TDEE/RMR ratio drops to 1.02 during sustained restriction and recovers to 1.14 when restriction eases. **My set point moves up or down only when some conditions are met**, and weight loss above it is trivial; below it, [fiendishly difficult](https://doi.org/10.2105/AJPH.2015.302773).

**The body defends a moving fat mass set point that only updates when weight is stable.** Identified primarily from published GLP-1/GIP withdrawal trials — where this subject's intake autocorrelation cannot bias the answer — and confirmed by three independent subject-side measurements. The rule: the set point adapts toward fat mass at 1% per day **only when FM has held within ±3 lbs of a reference for at least 14 consecutive days**. When FM is changing rapidly, the set point freezes. Pressure on intake while there's a gap: ~55 cal/day per lb, symmetric in both directions. This single rule predicts SURMOUNT-4 post-tirzepatide regain at +14.1% (published +14.0%) and STEP-1 post-semaglutide regain at +9.3% (published ~+10%). A simpler trailing-EMA model fits this subject's slow-gain decade equally well — but undershoots trial regain by 5-10×. **Trial data is what pins the model; this subject's data is the special case.**

**The 45-50 day timescale converges from four independent measurements.** Binge rate from set-point distance (binary outcome, no surplus involved): HL ≈ 50d. Intake-free fat mass from weight interpolation alone: HL = 50d. Weigh-in-only EMA test of distance against next-30-day intake: HL ≈ 35-45d. Trial regain inversion: same family, slightly slower latching. The continuous mean-surplus correlation also supports a slow timescale (best r = -0.92), but on its own it cannot pin a specific half-life: a 2-D sweep over SP half-life and surplus lookback window produces a broad ridge (HL=20d × 45d window, HL=45d × 90d window, and HL=75d × 120d window all give r ≈ -0.91 to -0.94). The mean-surplus regression and the rolling-FM SP share input bandwidth. What the timescale rests on is the *intersection* of that ridge with the four non-autocorrelated identifications above.

**Two channels, different speeds.** The body closes the gap on two fronts. An eating channel (~50-day timescale, ~55 cal/lb) shifts daily intake toward the set point. A metabolic rate channel (~9-day timescale, confirmed by 25 calorimetry measurements) burns ~90 cal/day more during weight loss than weight gain at the same body composition — actively assisting loss when fat mass is above the set point, passively permitting regain when below. The eating channel is stronger but slower. Binge size is constant regardless of distance from the set point; what shifts is the *frequency* and the depth of deficit days. The set point tilts the entire daily distribution rather than triggering discrete events.

**Receptor-dynamics reading.** The latch's stability requirement and the 45-50 day half-life are what hypothalamic leptin/melanocortin receptor recalibration would look like at a behavioral level: receptor populations adjust to chronic exposure, and the chronic exposure has to actually persist for the system to adapt. During drug-driven loss the satiety circuit is being stimulated by exogenous GLP-1/GIP — the *endogenous* receptor isn't being challenged at the new low FM, so it doesn't recalibrate. The set point stays anchored to pre-treatment, which is exactly what produces SURMOUNT-4-style regain on discontinuation. The same frame explains why GLP-1 trials asymptote where they do (see [FINDINGS.md § Asymptote](FINDINGS.md#asymptote-of-glp-1-weight-loss)) and why this subject's 2011-13 willpower era didn't produce binges until FM crossed back above ~18 lbs (the leptin-floor reading of [AC + BK](FINDINGS.md#ac--bk-the-2013-inflection-and-the-leptin-floor)).

**Tirzepatide suppresses appetite at -74.5 cal per unit of blood level** (AX, identified from within-week injection cycle variation — the cleanest estimate, free of set point and tachyphylaxis confounds). The set point's per-lb eating pressure continues at the same rate on and off drug; the drug subtracts from the total independently. This is visible in the weekly injection cycle: day 0 intake averages 1643 cal (drug peak), rising to 2222 cal by day 5 (drug trough), a 579 cal/day swing driven purely by drug pharmacokinetics. Tachyphylaxis erodes effectiveness as cumulative receptor exposure (AUC of blood level over time) accumulates — `tachy = exp(−K_AUC × ∫blood dτ)`, K_AUC ≈ 0.0034 / unit-week. Binge-to-binge escalation drops from 31.6% to 0%.

**The drug also suppresses the body's metabolic cooperation with weight loss.** Direct calorimetry shows RMR 206 cal below composition-predicted on the drug — the metabolic boost that normally accelerates fat loss is pharmacologically eliminated. This effect operates on longer timescales than the injection cycle and cannot be separated from the appetite effect within-week. Weight loss is 93% fat; the net energy budget is approximately -450 cal intake reduction, +200 cal metabolic cooperation lost, net ~-250 cal/day deficit.

**Using these parameters to simulate the SURMOUNT-1 trial** (Jastreboff et al. NEJM 2022, n=2539, zero parameters fitted to trial data): the model predicts 15mg weight loss of −19.5% vs published −22.5%, with a roughly uniform 2-3 percentage point under-prediction across all three dose arms (5mg, 10mg, 15mg). The under-prediction reflects an honest cross-subject difference in central drug sensitivity; what's notable is that the *shape* of the dose-response is now flat-error rather than ramping. SURMOUNT-2 over-predicts (T2D diverts GLP-1 to glucose control, not in the model). Post-discontinuation regain is captured by the SmoothLatch (+14.1% vs +14.0% in SURMOUNT-4).

**Why GLP-1 weight loss asymptotes** is a question the literature usually answers with "metabolic adaptation" without specifying a mechanism. This dataset's parameters give a quantitative one: at the asymptote, the eroded drug suppression `D₀ × tachy(t)` equals the per-lb pressure `55 × (SP − FM)` against a set point that has barely moved. Dose sets `D₀`; tachyphylaxis sets the timescale of asymptote drift; the SmoothLatch determines how much the set point has actually descended (very little during rapid loss, more during a long stable hold). See [FINDINGS.md § Asymptote of GLP-1 weight loss](FINDINGS.md#asymptote-of-glp-1-weight-loss).

**Not all restriction is equal.** Long runs (≥6 days) and low-carb restriction recover best. Low-protein and high-step restriction produce the largest metabolic penalties. Potato diets show zero binges across 69 days and high TDEE/RMR ratio, but severe rebound after stopping.

**What didn't survive.** The model above is what's left after testing every other candidate. Each elimination has a number and a reproducer in FINDINGS:
- [Gravitostat](https://doi.org/10.1073/pnas.1800033115) (foot-pounds drives intake): r = +0.05, wrong direction.
- [Speakman's dual-intervention zone of indifference](https://doi.org/10.1242/dmm.008698): no zone, continuous sigmoid.
- Set point tracking scale weight: dominated by FM tracking at every HL (-0.60 fat vs -0.52 scale).
- Fixed set point: r = +0.20 (and wrong sign) vs r = -0.60 for the moving 50d EMA.
- [Lowe et al. (2007)](https://doi.org/10.1002/eat.20405) fixed weight-suppression proxy: directional finding replicates, but their fixed reference is dominated by a moving 50d EMA.
- Symmetric EMA at any HL: undershoots SURMOUNT-4 regain by 5-10×.
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
| `drugs/tirzepatide.csv` | 560 | Daily PK blood level + tachyphylaxis |

Raw data in `intake/`, `weight/`, `composition/`, `RMR/`, `steps-sleep/`, `drugs/`, `temperature/`, `fatty-acids/`, `workout/`, `travel/`. All extractors idempotent. Pipeline in [CLAUDE.md](CLAUDE.md).

I accept any question, no matter how personal, in service of our goal. You will not offend me.
