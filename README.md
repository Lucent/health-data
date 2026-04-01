# 15 years of daily food intake, weight, body composition, and metabolic data

Fifteen years ago, I lost 100 lbs [in a fairly boring manner](https://www.reddit.com/r/loseit/comments/498afc/) by eating 1200 calories a day. It's creeping back. There are no low-hanging fruits — no sugarwater, no fast food, near-zero restaurant meals. Since then I've recorded every bite eaten every day, whether 10 calories of broth on a weekend fast or a 4000 calorie binge, using a food scale for everything. 54,925 food items across 5,429 days. Zero gaps.

This is not self-reported-from-memory dietary recall. The cumulative energy balance — intake minus expenditure, accumulated daily — closes to ±5 lbs over the full 15 years. A composition-aware model fitted to 25 indirect calorimetry measurements and 70 body composition scans quantifies the undercount at ~10-14%, uniform across weight loss, gain, and stable phases. **The best metabolic ward studies achieve ~5% error over 2-4 weeks. Self-reported dietary recalls average 30-50% error with missing days. This dataset sustains ~10-14% error across 5,429 consecutive days.**

**The dataset contains hundreds of overlapping natural experiments**: months of keto, weekend 36-hour fasts, low protein, high fiber, potatoes-only, waves of monotonous meals, daily ice cream, different cooking oils. Combined with simultaneous weight, body composition (BOD POD + InBody), resting metabolic rate (Cosmed indirect calorimetry), steps, sleep, body temperature, blood fatty acid panels, and 80 weekly tirzepatide injections logged from the first dose. **If the answer to obesity requires a complex overlap or sequence of conditions, it may be hidden within and first discovered through data mining rather than invented by a brilliant hypothesizer.** Very likely, the answer is already in the data — a few weeks spread across years where something worked, masked by the noise in daily weigh-ins.

I am [mostly](https://ptable.com) [retired](https://flightaware.com), 43, no prescriptions, no stressors, no emotional eating. No meat ([35 years](https://www.reddit.com/r/vegetarian/comments/35d7iq/)). No workplace contaminant exposure. Reverse-osmosis water. The worst-case common-variant obesity genotype: FTO homozygous risk at all five loci, MC4R heterozygous, UCP2 reduced thermogenesis — not [easily-solved obesity](https://doi.org/10.2105/AJPH.2015.302773). See [BACKGROUND.md](BACKGROUND.md) for the full genetic profile and health history.

## What the analysis has found

**The body defends a narrow expenditure band regardless of weight.** At 165 lbs (20 lbs fat, 2013) and 225 lbs (82 lbs fat, 2024), derived TDEE is ~2100-2200 cal/day. The TDEE/RMR ratio drops to 1.02 during sustained restriction and recovers to 1.14 when restriction eases. **My set point moves up or down only when some conditions are met**, and weight loss above it is trivial; below it, [fiendishly difficult](https://doi.org/10.2105/AJPH.2015.302773).

**The body defends a moving fat mass set point, and it controls how much you eat.** The set point is not a fixed number — it's an exponential moving average of recent fat mass with a **~45-day half-life** (every 6-7 weeks, the gap between defended weight and actual weight closes by half). After 5 months at a new weight, the set point has 87% converged. These numbers are derived from 2014-2026 only — the 2011-2013 period of aggressive willpower-driven restriction is excluded because intake was externally controlled, not responding to the set point. Including that period inflated the per-lb pressure and created a spurious asymmetry that did not replicate on natural-dynamics data.

The 90-day average of daily caloric surplus (intake minus expenditure) correlates with distance from the set point at **r = -0.92**. Each pound below the set point shifts mean daily intake by **~55 cal** [CI: 15-40]. The set point doesn't create binges — individual binge size is constant regardless of distance. Instead, it tilts the *entire* distribution of daily eating: surplus days become more frequent, deficit days become shallower. It's a continuous control on average energy balance, not a trigger for discrete events. **Caveat:** the set-point half-life and the surplus lookback are partially interchangeable smoothings. A 2D sweep shows a broad ridge (for example HL=20d with a 45d surplus window gives r=-0.94; HL=45-50d with a 90d window gives r=-0.93). So the surplus regression supports a slow appetite-pressure timescale, but does not uniquely identify 45 days by itself.

**The 45-day half-life survives an intake-free test.** The primary model (Kalman filter) uses logged intake to help estimate daily fat mass. This creates a circularity concern. But when fat mass is estimated from weight observations alone (P3: interpolating between 1,693 scale readings, no intake in the model), **the 50-day half-life still appears** (r = -0.73 on 2014+ data vs -0.92 with the fuller model). Three measurement paths — binge frequency (binary), mean surplus (continuous), and intake-free weight interpolation — converge on the same 45-50 day timescale. That convergence matters because the mean-surplus fit alone has an identifiability ridge: it supports a slow timescale, but needs the intake-free and binge-rate checks to argue that `~45-50d` is biological rather than just a convenient smoothing pair.

**Two channels, different speeds.** The body closes the gap on two fronts. An eating channel (~45-day timescale, ~55 cal/lb) shifts daily intake toward the set point. A metabolic rate channel (~9-day timescale, confirmed by 25 calorimetry measurements) burns ~90 cal/day more during weight loss than weight gain at the same body composition — actively assisting loss when fat mass is above the set point, passively permitting regain when below. The eating channel is stronger but slower.

**The set point only adapts when weight is stable.** On this subject's data alone, a simple trailing average (EMA) and a stability-gated model (SmoothLatch) both fit eating patterns at r = -0.94 — they are indistinguishable during slow weight change. But they make opposite predictions for rapid weight loss: the EMA predicts the set point follows weight down and regain after stopping treatment is minimal (+3%); the SmoothLatch predicts the set point **freezes** during rapid loss (FM never stabilizes long enough to trigger adaptation) and regain is substantial. Trial data breaks the tie decisively: the SmoothLatch predicts SURMOUNT-4 post-tirzepatide regain at +14.1% (published: +14.0%) and STEP-1 post-semaglutide regain at +9.3% (published: ~+10%). The EMA undershoots both by 5-10x. **The set point adapts when — and only when — fat mass holds within ±3 lbs for at least 14 consecutive days.** This predicts that regain depends on the *speed* of loss: faster loss leaves a larger un-adapted gap. Consistent with clinical observations that surgical and very-low-calorie patients regain more than gradual losers. See [FINDINGS.md §Set point](FINDINGS.md#set-point).

**Tirzepatide suppresses appetite at -74.5 cal per unit of blood level** (AX, identified from within-week injection cycle variation — the cleanest estimate, free of set point and tachyphylaxis confounds). The set point's per-lb eating pressure continues at the same rate on and off drug; the drug subtracts from the total independently. This is visible in the weekly injection cycle: day 0 intake averages 1643 cal (drug peak), rising to 2222 cal by day 5 (drug trough), a 579 cal/day swing driven purely by drug pharmacokinetics. Tachyphylaxis erodes effectiveness with a 35-week cumulative half-life. Binge-to-binge escalation drops from 31.6% to 0%.

**The drug also suppresses the body's metabolic cooperation with weight loss.** Direct calorimetry shows RMR 206 cal below composition-predicted on the drug — the metabolic boost that normally accelerates fat loss is pharmacologically eliminated. This effect operates on longer timescales than the injection cycle and cannot be separated from the appetite effect within-week. Weight loss is 93% fat; the net energy budget is approximately -450 cal intake reduction, +200 cal metabolic cooperation lost, net ~-250 cal/day deficit.

**Using these parameters to simulate the SURMOUNT-1 trial** (Jastreboff et al. NEJM 2022, n=2539, zero parameters fitted to trial data): the model predicts 15mg weight loss of -20.3% vs published -22.5% (within 2.2 percentage points). It overpredicts diabetic weight loss in SURMOUNT-2 (-21% simulated vs -15% published), consistent with the drug's appetite effect being partially diverted to glucose control. Post-discontinuation regain remains underpredicted (+1% vs +14% in SURMOUNT-4), the largest open discrepancy.

**Not all restriction is equal.** Long runs (≥6 days) and low-carb restriction recover best. Low-protein and high-step restriction produce the largest metabolic penalties. Potato diets show zero binges across 69 days and high TDEE/RMR ratio, but severe rebound after stopping.

**What doesn't work:** steps don't predict weight change at any timescale. Sleep hours, protein, fat, and fiber have no independent signal after controlling for calories, carbs, and sodium. The gravitostat hypothesis shows a weak signal in the wrong direction. The week-scale intake invariance claim was falsified.

Every finding with numbers is reproduced by a standalone script in [FINDINGS.md](FINDINGS.md). The [gravitostat](https://doi.org/10.1073/pnas.1800033115), omega-6:3 ratios, circadian misalignment, and other untested theories are also documented there.

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
