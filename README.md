# 15 years of daily food intake, weight, body composition, and metabolic data

Fifteen years ago, I lost 100 lbs [in a fairly boring manner](https://www.reddit.com/r/loseit/comments/498afc/) by eating 1200 calories a day. It's creeping back. There are no low-hanging fruits — no sugarwater, no fast food, near-zero restaurant meals. Since then I've recorded every bite eaten every day, whether 10 calories of broth on a weekend fast or a 4000 calorie binge, using a food scale for everything. 54,925 food items across 5,429 days. Zero gaps.

This is not self-reported-from-memory dietary recall. The cumulative energy balance — intake minus expenditure, accumulated daily — closes to ±5 lbs over the full 15 years. A composition-aware model fitted to 25 indirect calorimetry measurements and 70 body composition scans quantifies the undercount at ~10-13%, uniform across weight loss, gain, and stable phases. **The best metabolic ward studies achieve ~5% error over 2-4 weeks. Self-reported dietary recalls average 30-50% error with missing days. This dataset sustains ~10-13% error across 5,429 consecutive days.**

**The dataset contains hundreds of overlapping natural experiments**: months of keto, weekend 36-hour fasts, low protein, high fiber, potatoes-only, waves of monotonous meals, daily ice cream, different cooking oils. Combined with simultaneous weight, body composition (BOD POD + InBody), resting metabolic rate (Cosmed indirect calorimetry), steps, sleep, body temperature, blood fatty acid panels, and 80 weekly tirzepatide injections logged from the first dose. **If the answer to obesity requires a complex overlap or sequence of conditions, it may be hidden within and first discovered through data mining rather than invented by a brilliant hypothesizer.** Very likely, the answer is already in the data — a few weeks spread across years where something worked, masked by the noise in daily weigh-ins.

I am [mostly](https://ptable.com) [retired](https://flightaware.com), 43, no prescriptions, no stressors, no emotional eating. No meat ([35 years](https://www.reddit.com/r/vegetarian/comments/35d7iq/)). No workplace contaminant exposure. Reverse-osmosis water. The worst-case common-variant obesity genotype: FTO homozygous risk at all five loci, MC4R heterozygous, UCP2 reduced thermogenesis — not [easily-solved obesity](https://doi.org/10.2105/AJPH.2015.302773). See [BACKGROUND.md](BACKGROUND.md) for the full genetic profile and health history.

## What the analysis has found

**The body defends a narrow expenditure band regardless of weight.** At 165 lbs (20 lbs fat, 2013) and 225 lbs (82 lbs fat, 2024), derived TDEE is ~2100-2200 cal/day. The TDEE/RMR ratio drops to 1.02 during sustained restriction and recovers to 1.14 when restriction eases. **My set point moves up or down only when some conditions are met**, and weight loss above it is trivial; below it, [fiendishly difficult](https://doi.org/10.2105/AJPH.2015.302773).

**The set point defends through expenditure, not intake.** Five tests of intake-side homeostasis all come back negative. Binges cluster rather than compensate. Distance from set point does not predict binge probability (AUC=0.49). After weekend fasts, there is zero compensatory overeating — but the deficit disappears through reduced expenditure within a week.

**Tirzepatide works by silencing behavioral turbulence.** Blood level predicts daily intake at r=-0.50 with a 568 cal/day weekly appetite swing. But the drug's strongest effect is on state transitions: binge-to-binge persistence drops from 31.6% to 0%. Direct calorimetry confirms the body still claws back ~200 cal/day through reduced RMR on the drug — the metabolic defense is NOT overridden. Weight loss is 93% fat.

**Not all restriction is equal.** Long runs (≥6 days) and low-carb restriction recover best. Low-protein and high-step restriction produce the largest metabolic penalties. Potato diets show zero binges across 69 days and high TDEE/RMR ratio, but severe rebound after stopping.

**What doesn't work:** steps don't predict weight change at any timescale. Sleep hours, protein, fat, and fiber have no independent signal after controlling for calories, carbs, and sodium. The gravitostat hypothesis shows a weak signal in the wrong direction. The week-scale intake invariance claim was falsified.

Every finding with numbers is reproduced by a standalone script in [FINDINGS.md](FINDINGS.md). The [gravitostat](https://doi.org/10.1073/pnas.1800033115), omega-6:3 ratios, circadian misalignment, and other untested theories are also documented there.

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
