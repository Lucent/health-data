# Direct audit of named diet epochs

## Bottom line

The named diet experiments reveal a practical pattern, but they do **not**
identify a special fat-loss physiology:

- All four potato attempts lowered logged intake by **524–748 kcal/day versus
  their immediately preceding calendar windows**. Three short attempts showed
  rapid scale declines, but **3.3–3.9 lb returned in each measured immediate
  post-window**. The most defensible interpretation is that potatoes were an
  effective short-term food-environment/satiety constraint, while much of the
  rapid scale movement was water, glycogen, sodium, and gut content.
- Seven weekend fasts lowered the 44-day mean to **1,624 kcal/day**, but ordinary
  eating days still averaged **2,335 kcal/day** and the subsequent 28-day period
  averaged **2,364 kcal/day**. Observed weight fell 4.2 lb during the measured
  portion and then rose 7.3 lb in the post-window. This is a strong warning
  about rebound and short-term scale illusions, not proof that fasting damaged
  metabolism.
- The two keto annotations have **zero observed scale weights**. Keto phase 1
  also overlaps the high-protein experiment for 14 days. Intake was 440
  kcal/day above its pre-window in phase 1 and essentially unchanged in phase
  2, so this repository cannot support a keto-specific weight-loss claim.
- The high-protein experiment raised logged protein, but it has only one scale
  measurement and overlaps keto. It is evidence that a protein target was
  behaviorally achievable, not that protein caused weight loss in this person.
- Low protein did not suppress intake or weight: 34 of 37 days met the
  prespecified low-protein rule, intake averaged 2,340 kcal on those days, and
  observed weight rose 1.3 lb across the measured experiment.
- Enlightened-brand products appeared on 338 of 426 days in the “ice cream”
  era. Those days had lower intake than non-product days, but there was only one
  scale observation in the entire era. This contradicts a simplistic “any
  ultra-processed dessert inevitably causes overeating” story, but selection
  and the surrounding menu make it non-causal.
- Both labeled travel binges raised intake sharply—by **795** and **632
  kcal/day** versus their pre-windows—while providing too little scale coverage
  to measure weight effects. Travel remains a high-confidence intake-risk
  context.

The useful lesson is therefore to borrow the **repeatable, low-choice,
high-satiety structure** of the potato experiments without copying their very
low protein or treating a rapid water-weight drop as the objective.

## What this audit uses

The analysis reads only:

- `intake/diet_epochs.csv` for the annotations;
- `intake/intake_daily.csv` for observed daily totals;
- `intake/intake_foods.csv` for logged product/food confirmation; and
- `weight/weight.csv` for exact observed scale measurements.

It does not read P1, P2, P3, P4, modeled TDEE, modeled fat mass, modeled
medication exposure, or interpolated weights.

For every annotation, the script compares:

1. the full annotated epoch;
2. an immediately preceding window of equal length, capped at 28 days; and
3. an immediately following window of equal length, capped at 28 days.

Long epochs therefore get local 28-day flanks, while short interventions get
same-length flanks. All intake dates are present. Scale summaries retain the
first and last **actually observed** values inside each window; a blank means
there is no usable anchor.

The food-log confirmation rules were fixed before reading the contrasts:

- potato: at least 500 logged kcal from a food name containing potato, Yukon,
  or russet;
- weekend fasting: Saturday or Sunday with less than 200 kcal;
- low protein: protein at or below 10% of calories;
- high protein: protein at or above 20% of calories;
- keto: carbohydrate at or below 25% and fat at or above 50% of calories;
- Enlightened era: any logged Enlightened-brand calories;
- travel binge: at least 3,000 kcal; and
- post-travel restriction: at most 1,800 kcal.

These are audit definitions, not clinical thresholds. MyFitnessPal food
entries—including unusually high historical fiber values—inherit the accuracy
limits of the user-entered food database.

## Annotation overlap and contamination

The annotation table is not a partition. There are two literal overlaps:

| Annotation 1 | Annotation 2 | Shared dates |
|---|---|---:|
| high-protein experiment | keto phase 1 | 14 |
| Slate chocolate-milk era | potato diet 1 | 13 |

The high-protein epoch has only 17 calendar-exclusive days. Its full-period
average was 1,866 kcal and 119.8 g protein; its 17 exclusive days averaged
1,745 kcal and 120.1 g protein. Keto phase 1 has 72 exclusive days, averaging
2,231 kcal. The overlap is large enough that neither label should be credited
with the October 2017 result by itself.

Potato attempt 1 has zero annotation-exclusive days because it is entirely
nested in the two-year Slate era. That broad era is context rather than a
competing acute intervention, but the overlap still prevents a clean
single-label interpretation.

## Potato attempts

| Attempt | Logged exposure | Mean kcal before | Mean kcal during | During minus before | Mean kcal after | Exact weights during | First → last observed | Exact weights after | First → last observed after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2022-11-13 to 11-25 | 13/13 days | 2,378 | 1,716 | -663 | 2,769 | 13 | 219.4 → 211.9 (-7.5) | 8 | 212.4 → 215.7 (+3.3) |
| 2023-01-02 to 02-10 | 15/40 days | 2,754 | 2,006 | -748 | 2,326 | 40 | 222.8 → 215.7 (-7.1) | 22 | 214.5 → 217.9 (+3.4) |
| 2023-08-20 to 08-25 | 6/6 days | 2,379 | 1,822 | -557 | 2,936 | 5 | 219.1 → 219.4 (+0.3) | 3 | 218.8 → 222.7 (+3.9) |
| 2024-02-26 to 03-06 | 10/10 days | 2,491 | 1,967 | -525 | 2,472 | 10 | 226.1 → 222.4 (-3.7) | 9 | 221.6 → 224.9 (+3.3) |

Important qualifications:

- Attempt 2 was not a continuous 40-day potato diet. Food logs confirm 15
  potato days interspersed among 25 ordinary days. Confirmed days averaged
  **1,505 kcal**, compared with **2,307 kcal** on the other labeled days.
- Exact next-calendar-day scale changes in attempt 2 averaged **-1.31 lb after
  a confirmed potato day** and **+0.45 lb after another labeled day**. A
  biological fat change of that size overnight is implausible; this is direct
  evidence that water/gut-content effects dominate the short-run scale signal.
- Attempts 1, 3, and 4 were continuous food-log-confirmed attempts. Logged
  potato foods supplied about 56%, 77%, and 68% of total calories,
  respectively. Attempt 2's confirmed potato days got about 53% of their
  calories from logged potato foods.
- Mean protein was only 26 g/day in attempt 1, 31 g/day in attempt 3, and 36
  g/day in attempt 4. Repeating these as a long-term diet would conflict with
  the present goal of preserving lean tissue.
- The attempts changed several things simultaneously: calories, menu variety,
  carbohydrate, protein, sodium, food mass, and sometimes added fat type. They
  cannot distinguish potato-specific effects, olive oil versus butter, “seed
  oils,” carbohydrate-insulin mechanisms, or metabolic advantage.

What survives the caveats is unusually consistent: a highly repetitive
potato-centered menu made lower intake feasible on all four attempts. A better
prospective translation is a calorie-controlled, high-protein, high-fiber
repetitive menu—not a prolonged potato monodiet.

## Weekend fasting

From 2019-10-12 through 2019-11-24:

- all 14 Saturdays and Sundays were food-log-confirmed fast days;
- fast days averaged **99 kcal** and 3.1 g protein;
- the other 30 days averaged **2,335 kcal** and 85.8 g protein;
- the full epoch averaged **1,624 kcal/day**, versus **2,084** in the preceding
  28 days and **2,364** in the following 28 days;
- there were 31 observed weights during the epoch, from 213.2 lb on October 23
  to 209.0 lb on November 24; and
- the post-window ran from 208.3 lb on November 25 to 215.6 lb on December 21,
  a 7.3 lb rise.

Among exact next-day scale pairs, a fast day was followed by a mean **-1.93
lb** change, while an ordinary epoch day was followed by **+0.79 lb**. Those
changes are far too rapid to interpret as body fat. They show just how strongly
this schedule moved water, glycogen, sodium, and gut contents.

The fasting schedule created a real weekly calorie deficit, so it may have
reduced tissue mass during the period. But the large post-period intake and
scale rebound make it a poor default strategy for this individual unless a
future preregistered trial shows that a gentler schedule prevents compensation.

## Keto and high protein

| Epoch | Rule-confirmed days | Mean kcal before | Mean kcal during | Difference | Observed scale days |
|---|---:|---:|---:|---:|---:|
| high protein | 25/31 | 2,032 | 1,866 | -166 | 1 |
| keto phase 1 | 63/86 | 1,756 | 2,195 | +440 | 0 |
| keto phase 2 | 42/60 | 2,114 | 2,118 | +4 | 0 |

On confirmed high-protein days, intake averaged 1,823 kcal and 126.6 g
protein, versus 2,048 kcal and 91.5 g on the other six days. This is compatible
with protein leverage or improved satiety, but the result is not separable from
date, menu, and the keto transition. There is no weight trajectory.

Keto phase 1 did not coincide with lower intake relative to its short
pre-window. Phase 2 was calorie-neutral. With no observed scale weights in
either phase, the data cannot test whether carbohydrate restriction altered
weight, appetite, energy expenditure, or water independently of calories.
This is missing evidence, not evidence that keto failed.

## Low protein

The 37-day low-protein annotation averaged 2,324 kcal and 45.1 g protein,
compared with 2,376 kcal and 50.6 g protein in the preceding 28 days. The
pre-window was already low in protein, making the contrast only 5.5 g/day.

Thirty-four of 37 days met the at-most-10%-of-calories rule; those days averaged
2,340 kcal and 43.7 g protein. Across 29 observed scale dates, weight went from
221.5 to 222.8 lb (+1.3 lb). The following window is contaminated by the
six-day potato attempt and a subsequent intake rebound, so it cannot isolate a
lagged low-protein effect.

This does not support low protein as a weight-loss tool here. It also does not
provide a clean test of the protein-leverage hypothesis because the pre-period
was already protein-poor and there was no controlled calorie target.

## Enlightened / “ice cream” era

During the 426-day epoch:

- Enlightened products were logged on 338 days (79.3%);
- those products supplied 146,658 kcal, about 16.7% of all epoch calories;
- product days averaged 2,030 kcal, compared with 2,184 kcal on the other 88
  days; and
- at-least-3,000-kcal days occurred on 3.3% of product days versus 11.4% of
  non-product days.

Only one scale weight exists in the entire epoch. Product choice was not
random: Enlightened may have been deliberately used on controlled-menu days,
and non-product days may disproportionately contain travel, celebrations, or
logging differences. The result therefore says only that this portioned
dessert was compatible with lower-intake days. It cannot establish that the
product caused satiety or weight loss.

This is still a useful counterexample to theories that classify foods solely
by “processed” status. In this history, menu context, portioning, repetition,
and total intake appear more predictive than a universal processed/unprocessed
binary.

## Travel and post-travel restriction

| Epoch | Mean kcal before | Mean kcal during | Difference | Days at least 3,000 kcal | Observed-scale limitation |
|---|---:|---:|---:|---:|---|
| September 2015 travel | 2,152 | 2,947 | +795 | 13/23 | 5 weights, all within the first 6 days |
| June 2019 travel | 1,907 | 2,539 | +632 | 4/30 | no weights |
| July–September 2019 restriction | 2,546 | 1,845 | -701 | 2/92 | no weights |

The 2015 within-epoch scale values declined 2.2 lb, but because all five occur
at the beginning of the trip, they do not measure the full exposure and may
reflect travel hydration. The direct result worth acting on is the replicated
intake rise, not a weight estimate.

## How the findings bear on obesity theories

- **Energy intake / food environment:** strongly supported as the immediate
  behavioral pathway. Repetitive potato menus and fasting schedules lowered
  recorded intake; travel raised it.
- **Satiety and protein leverage:** plausible but not cleanly identified.
  High-protein days were lower-calorie, while low-protein days were not, but
  calendar/menu confounding is large. The potato periods show that food volume
  and monotony can overwhelm low protein in the short run.
- **Carbohydrate-insulin or keto advantage:** untested because there are no
  weights during keto and calories did not fall consistently.
- **Seed oils / added-fat type:** untested. Potato attempts differed in oil,
  butter, cream, and alcohol while all lowered calories. No comparison isolates
  linoleic acid or a specific fat.
- **Ultra-processed-food inevitability:** too simple for these observations.
  A portion-controlled processed dessert co-occurred with lower intake, though
  causality is unknown.
- **“Metabolic damage” from fasting:** not established. The rebound is directly
  observed, but behavioral compensation, holidays, water, glycogen, and sodium
  are sufficient explanations.

## Actionable translation

The most evidence-aligned experiment is not another extreme diet. It is a
stable, repetitive menu that retains the potato trials' low-choice/high-volume
properties while correcting their weaknesses:

1. keep calories in a modest, prespecified range rather than alternating near
   zero and rebound days;
2. include a meaningful protein source at each meal;
3. pre-plan travel defaults, because travel produced the most replicated
   upward intake shock;
4. judge progress by multiweek observed-weight averages plus waist/body
   composition under consistent conditions, not next-day scale drops; and
5. change one major variable at a time and preregister the analysis window.

Any medication change, altered injection interval, or dosing experiment remains
clinician/pharmacist-directed and is outside this diet-epoch audit.

## Reproduce

```bash
uv run --with pandas python analysis/EJ_diet_epoch_direct_audit.py
```

Outputs:

- `analysis/EJ_epoch_window_summary.csv`
- `analysis/EJ_epoch_contrasts.csv`
- `analysis/EJ_epoch_overlaps.csv`
- `analysis/EJ_confirmed_exposure_summary.csv`
