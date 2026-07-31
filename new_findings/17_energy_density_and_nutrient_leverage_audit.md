# Energy density and nutrient leverage: coverage first, prediction mostly null

Date: 2026-07-30
Status: direct food-ledger audit with ordered holdouts

## Bottom line

This repository does **not** measure whole-diet energy density well enough to
test the usual “eat lower-energy-density food” hypothesis retrospectively.
Exact consumed mass is recoverable for only **17.0% of solid-food calories**
across the full ledger and **10.8% in the final 90 days**. Restaurant foods,
recipes, pizza, and generic “serving” units are disproportionately missing, so
the weighed subset is not representative of the whole diet.

The other major findings are more useful:

- Food-row calories, daily calories, protein, and fiber reconcile exactly on
  all 5,429 days. Every meal label is one of the five expected labels.
- Protein and fiber field availability is high in the modern MHTML source, but
  old OXPS zeroes are ambiguous because that extractor rendered unknown values
  as zero.
- Same-day nutrient-density, item-count, snack-share, and meal-count
  associations are mostly mathematical coupling or descriptions of how a
  calorie total was assembled. They are not evidence of satiety.
- None of the prespecified features added a stable, practically useful
  next-day calorie forecast across the 2022 validation and 2023–2024
  confirmation periods.
- The directly weighed subset's energy density improved prediction in one
  sparse holdout and failed in the next. Its direction also changed across
  same-day splits. That is an unstable selected-subset result, not a usable
  intervention signal.
- The final-90-day meal table reinforces the existing nutrition plan:
  breakfast and snacks have the lowest protein and fiber density. Beverage-like
  calories average only **72 kcal/day (3.6%)**, so they are not the main current
  leverage point.

The practical implication is subtle but important: improve protein and fiber
because the current absolute intakes and meal distribution are weak, not
because a historical nutrient-density regression predicts that doing so will
automatically lower tomorrow's calories.

## Reproducer and artifacts

Command:

```bash
uv run --with numpy --with pandas python analysis/EN_energy_density_nutrient_leverage.py
```

Artifacts:

- `analysis/EN_data_coverage_by_year_source.csv`
- `analysis/EN_daily_features.csv`
- `analysis/EN_feature_split_summary.csv`
- `analysis/EN_same_day_markers.csv`
- `analysis/EN_locked_prediction_results.csv`
- `analysis/EN_current_meal_summary.csv`
- `analysis/EN_liquid_food_summary.csv`

The script reads only the canonical food/day ledgers, the checked-in MFP
enrichment table for parsed serving units, and the checked-in travel table. It
does not read P1, P2, P3, P4, modeled TDEE, modeled body composition, or drug
PK quantities.

## Three possible food-mass approaches

1. **Use only explicit mass units.** This is the conservative approach used
   here. It avoids fabricated weights but leaves poor, nonrandom coverage.
2. **Infer grams from generic servings, recipes, or current database matches.**
   This would increase apparent coverage, but historical product sizes,
   restaurant portions, recipe yields, and serving definitions are not stable.
   It would replace missingness with undocumented imputation.
3. **Measure prospectively.** Record served edible grams, leftovers, and
   beverage volume for a short preregistered experiment. This is the only
   defensible route to a subject-specific whole-diet energy-density test.

The first approach is best for the retrospective audit; the third is best if
the hypothesis is important enough to test.

## What is actually measured

There is no food-device field. The relevant provenance is extraction source:
HTML, OXPS, or MHTML.

| Source | Days | Food rows | Protein field present | Fiber field present | Exact-mass share of solid kcal |
|---|---:|---:|---:|---:|---:|
| HTML | 622 | 6,040 | 95.4% | 80.1% | 13.2% |
| OXPS | 1,885 | 17,595 | 97.6% | 96.3% | 17.5% |
| MHTML | 2,922 | 30,818 | 96.7% | 88.5% | 17.4% |
| All | 5,429 | 54,453 | 96.8% | 90.1% | 17.0% |

The OXPS percentages overstate semantic completeness. In that source, **1,756
protein zeroes** and **5,222 fiber zeroes** could mean either a true zero or an
unknown field. If all ambiguous zeroes are excluded, the lower-bound OXPS row
coverage is 87.6% for protein and 66.7% for fiber. This is a data-quality bound,
not an estimate that all excluded values were actually missing.

Additional checks:

- Calories are numeric on 100% of food rows.
- Food-row calories reconcile to the daily ledger with maximum absolute error
  **0 kcal**.
- Food-row protein and fiber reconcile to the daily ledger with maximum
  absolute error **0 g**.
- Meal-label validity is 100%: Breakfast, Lunch, Dinner, Snacks, or Supper.
- The label does not contain a clock time. “Supper” and “Snacks” are user
  categories, not validated circadian exposures.
- Current MFP enrichment resolves the food strings for 89.5% of food rows, but
  the nutrient-density analysis uses the canonical diary values,
  not substituted enrichment values.

### Why the fiber series needs interpretation

Historical fiber is not a uniform whole-food measure. Several modern high-fiber
days are dominated by low-carbohydrate tortillas, ChocZero syrups, and similar
products reporting very large amounts of functional fiber. For example,
several syrup entries report 44–56 g fiber for about 105–128 kcal. Those may
reflect the product labels, but they make historical “fiber density” a mixture
of ordinary foods, added isolated fibers, formulation changes, and database
entries. A fiber regression therefore cannot identify the effect of vegetables,
legumes, whole grains, or a particular fiber type.

## Whole-diet energy density is unavailable

Direct-mass solid-calorie coverage is:

| Ordered period | Exact-mass share of solid kcal | Days with at least 25% coverage |
|---|---:|---:|
| Discovery, 2018-03-01–2021-12-31 | 17.1% | 370 |
| Validation, 2022 | 18.2% | 99 |
| Confirmation, 2023-01-01–2024-09-16 | 18.6% | 187 |
| Treatment train, 2024-09-17–2025-08-18 | 17.5% | 94 |
| Treatment holdout, 2025-08-19–2026-02-28 | 13.0% | 39 |
| Final 90 days | 10.8% | — |

The audit calculates kcal/g only for the explicitly weighed solid-food subset
and always exports its coverage. Calling that value “daily diet energy density”
would be wrong. A day with weighed candy and an unweighed restaurant meal can
look high-density; a day with weighed vegetables and an unweighed pizza can
look low-density. The missingness mechanism is part of the exposure.

## Same-day markers are not intervention effects

The table below gives adjusted kcal/day per one discovery-period standard
deviation of each feature. Controls are prior-28-day intake, weekday, season,
linear time, and travel. Intervals use deterministic Newey-West standard
errors.

| Feature | Discovery | Validation | Confirmation | Main interpretation problem |
|---|---:|---:|---:|---|
| Protein g/1,000 kcal | -208 | -138 | -161 | Calories are the denominator |
| Fiber g/1,000 kcal | -203 | -248 | -242 | Calories are the denominator; fiber type changes |
| Beverage-like calorie share | +9 | +47 | +124 | Beverage calories are part of the same total; first two intervals cross zero |
| Snack calorie share | +47 | +95 | +81 | Snack calories are part of the same total |
| Food-row count | +142 | +213 | +171 | More entries are being summed into the outcome |
| Food-cluster count | +135 | +183 | +151 | More food sources are being summed into the outcome |
| Number of meal labels used | +70 | +159 | +156 | More eating occasions contribute calories |
| Weighed-subset kcal/g | -12 | +247 | +75 | Sparse selected subset; direction is unstable |

Protein and fiber density look “protective” on the current day, while item
count looks “harmful,” almost by construction. These results can describe the
menu but cannot establish that adding protein, adding fiber, or deleting an
item would reproduce the coefficient.

## Locked future prediction is null or unstable

The primary chronological design was fixed as:

- train: 2018-03-01–2021-12-31;
- validation: 2022;
- confirmation: 2023-01-01–2024-09-16.

Each extended model added one prespecified feature to the same baseline:
current calories, prior-28-day mean, weekday, season, time, and present/future
travel flags. All scaling and coefficients were fit on the training period
only. The score below is extended-model MAE minus baseline MAE; negative is
better. Baseline next-day MAE was 371 kcal in validation and 317 kcal in
confirmation.

| Feature | Validation ΔMAE | Confirmation ΔMAE | Treatment holdout ΔMAE |
|---|---:|---:|---:|
| Protein g/1,000 kcal | +3.4 | -1.1 | +0.3 |
| Fiber g/1,000 kcal | +0.1 | +1.7 | +5.2 |
| Beverage-like calorie share | +1.4 | +2.1 | -2.2 |
| Snack calorie share | -0.3 | +3.3 | +5.5 |
| Food-row count | -0.6 | +1.0 | +0.0 |
| Food-cluster count | -1.2 | +0.6 | -0.0 |
| Number of meal labels used | -1.2 | +2.1 | -0.1 |
| Exact-mass coverage | +1.6 | -0.3 | +0.1 |
| Weighed-subset kcal/g | -26.8 | +3.5 | unavailable |

No general feature improves next-day MAE in both pre-treatment holdouts. The
treatment transport check trains through 2025-08-18 and scores the later 12.5
mg-assigned period; it is a regime-transport sensitivity, not a new randomized
validation. The largest favorable treatment result is only about 2 kcal MAE
against a 287-kcal baseline.

The three-day-mean outcome does not rescue a practical claim. Snack share and
food-cluster count improve MAE in both pre-treatment holdouts, but only by
0.2–1.8 kcal against 232–308-kcal baselines. Weighed-subset kcal/g improves
three-day MAE by 16.8 kcal in the 98-day validation sample and 2.3 kcal in the
187-day confirmation sample, yet fails next-day confirmation, has opposite
same-day directions, and exists only on a nonrepresentative subset. This is an
exploratory lead at most.

The home-only sensitivity produces the same overall conclusion. Travel
stratification does not turn any of these food-composition variables into a
replicated predictor, and the independent paired travel analysis remains much
stronger.

## Current meal leverage

The final 90 days, 2025-12-01 through 2026-02-28, average 2,004 kcal, 54.6 g
protein, and 14.1 g fiber. Meal labels reveal where calorie-neutral
substitutions are easiest:

| Meal label | Share of calories | kcal on logged days | Protein g | Fiber g | Protein g/1,000 kcal | Fiber g/1,000 kcal |
|---|---:|---:|---:|---:|---:|---:|
| Dinner | 26% | 518 | 19.1 | 4.6 | 36.9 | 8.9 |
| Lunch | 25% | 515 | 15.7 | 4.0 | 30.6 | 7.9 |
| Supper | 18% | 399 | 11.5 | 2.8 | 29.0 | 7.1 |
| Breakfast | 18% | 364 | 6.4 | 1.9 | 17.5 | 5.2 |
| Snacks | 13% | 316 | 4.1 | 1.3 | 12.8 | 4.0 |

Breakfast and snacks supply about 31% of calories but have by far the lowest
protein and fiber density. This independently supports the existing plan to
start calorie-neutral substitutions there rather than redesigning every meal.

The conservative beverage heuristic identifies 72 kcal/day, or 3.6% of final
90-day energy. Its largest historical sources are milk, Slate protein/chocolate
milk, cocktails, coffee drinks, beer, soda, juice, and protein shakes; the
full list is exported for manual inspection. Liquid calories can matter on an
individual restaurant, alcohol, or travel day, but they are not the dominant
average current source and did not add a stable locked predictor.

## Practical conclusion

1. Keep the current protein/fiber substitution priority, especially at
   breakfast and snacks.
2. Do not claim that nutrient density has a demonstrated next-day satiety
   effect in this subject. The rationale is adequacy, lean-mass support,
   cardiometabolic quality, and menu leverage.
3. Do not use retrospective kcal/g estimates to rank foods or diets. Historical
   mass coverage is too low and too selected.
4. Treat beverage calories as a context-specific secondary lever, not the main
   explanation for the current plateau.
5. If energy density is worth testing, preregister a short prospective protocol
   that captures total served edible grams, leftovers, beverages, travel,
   medication stability, hunger, unplanned additions, and total calories.

This audit narrows the hypothesis space: the repo supports changing the current
menu's nutrient allocation, but not a broad claim that energy density, liquid
calories, meal count, or fiber density is the hidden driver of future intake.
