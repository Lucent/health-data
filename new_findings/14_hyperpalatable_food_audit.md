# Hyper-palatable-food audit: a same-day marker, not a replicated future-risk signal

Date: 2026-07-30

## Bottom line

The repository can measure **one** published hyper-palatable-food cluster reasonably well in the modern pre-tirzepatide period: foods with both more than 20% of calories from fat and more than 20% from total sugar. This is the Fazzino fat + simple-sugar (FS) criterion.

It cannot validly measure total Fazzino-style hyper-palatable-food exposure. The two sodium clusters require the weight in grams of each consumed food, and exact consumed weight is recoverable for only **17.0% of solid-food calories**. Missingness is plainly non-random: pizza, recipes, restaurant foods, and serving-based foods are disproportionately hard to weigh. Those entries are left unclassified rather than silently called non-hyper-palatable.

Within the assessable FS cluster:

- Higher exposure consistently marked higher **same-day** calories.
- It did **not** produce a decision-grade next-day calorie or extreme-intake signal after controlling for current intake, prior intake, weekday, linear time, season, and travel.
- The result does not show that fat+sugar combinations, ultra-processing, seed oils, sugar, or any individual ingredient caused the higher same-day intake. It also does not show that reducing them would cause weight loss.

This is useful mostly as a negative prioritization result: travel and menu disruption remain better prospective targets than a generic “hyper-palatable food” rule.

## Reproducer and evidence boundary

Command:

```bash
uv run --with numpy --with pandas python analysis/EJ_hyperpalatable_food_audit.py
```

Artifacts:

- `analysis/EJ_hpf_daily.csv`
- `analysis/EJ_hpf_results.csv`
- `analysis/EJ_hpf_split_summary.csv`
- `analysis/EJ_hpf_food_summary.csv`

The script reads only:

- `intake/intake_foods.csv`
- `intake/intake_daily.csv`
- `intake/mfp_food_dict.csv`
- `travel/trips.md`

It does not read P1, P2, P3, P4, modeled TDEE, modeled body composition, or modeled drug levels.

## Prespecified definition

[Fazzino et al. (2019)](https://doi.org/10.1002/oby.22639) proposed three solid-food clusters:

| Cluster | Published criterion |
|---|---|
| Fat + sodium (FSOD) | >25% kcal from fat and ≥0.30% sodium by food weight |
| Fat + simple sugar (FS) | >20% kcal from fat and >20% kcal from total sugar |
| Carbohydrate + sodium (CSOD) | >40% kcal from nonsugar, nonfiber carbohydrate and ≥0.20% sodium by food weight |

Consistent with the paper, carbohydrate calories for CSOD subtract both sugar and fiber, and liquids are excluded. The paper itself describes these cut points as a starting point that needs predictive validation and warns against assuming they are fixed or final.

The diary provides calories, fat, and sugar directly for most modern entries. Enrichment fills some missing nutrients. Diary values take priority; there is no food-name classification. For the sodium criteria, this audit uses only entries recorded directly in grams, ounces, or kilograms. It does not guess a weight for “one slice,” “one serving,” a restaurant item, or a recipe.

## Coverage

The analysis period ends before tirzepatide to avoid conflating exposure with the treatment transition.

| Split | Calendar days | Primary eligible days | Calorie-weighted FS coverage | Exact-mass sodium-cluster coverage | Mean FS share on eligible days |
|---|---:|---:|---:|---:|---:|
| Discovery, 2018-03-01–2021-12-31 | 1,402 | 1,274 | 93.8% | 16.3% | 14.3% |
| Validation, 2022 | 365 | 356 | 96.9% | 18.5% | 12.9% |
| Confirmation, 2023-01-01–2024-09-16 | 625 | 585 | 95.0% | 17.4% | 20.1% |

Primary days require at least 500 logged kcal and at least 80% of solid-food calories assessable for FS. The discovery-period 75th percentile, **22.89% FS calories**, was locked as a descriptive high-exposure threshold before applying it to validation and confirmation.

Across the whole analysis period, FS classification covers 94.6% of solid-food calories. A stricter sensitivity requiring the recorded macros to agree with stated calories within 20% covers 85.7%.

## Adjusted results

Effects below are per **10 percentage points** more FS calories among assessable solid foods. Calorie outcomes are kcal/day. Extreme-intake outcomes are absolute probability changes for a ≥3,000-kcal day, so `0.003` means 0.3 percentage points.

### Same day

| Outcome | Discovery | Validation | Confirmation |
|---|---:|---:|---:|
| Total calories | +50.0 `[+31.8, +68.3]` | +114.1 `[+60.4, +167.8]` | +65.8 `[+37.3, +94.4]` |
| ≥3,000-kcal probability | +0.014 `[+0.004, +0.023]` | +0.019 `[-0.007, +0.044]` | +0.006 `[-0.010, +0.022]` |

The calorie association replicates in direction and interval exclusion. It is still not a clean test of causality: the exposure is calculated from foods contributing to the same calorie total, and a high-FS day may also be a dessert, restaurant, social, novelty, or permissive-eating day.

### Next day

The next-day models additionally control for current-day calories.

| Outcome | Discovery | Validation | Confirmation |
|---|---:|---:|---:|
| Next-day calories | +13.9 `[-2.7, +30.5]` | +34.0 `[-12.1, +80.2]` | +16.6 `[-11.2, +44.4]` |
| Next-day ≥3,000-kcal probability | +0.003 `[-0.005, +0.011]` | +0.004 `[-0.024, +0.031]` | +0.002 `[-0.016, +0.020]` |

The calorie point estimate is positive in all three splits but small, and every interval includes zero. Partial correlations are only `0.045`, `0.094`, and `0.054`. There is no credible replication for next-day extreme-intake risk.

Removing all travel days and pre-travel days leaves the same conclusion:

| Outcome | Discovery | Validation | Confirmation |
|---|---:|---:|---:|
| Next-day calories | +12.7 `[-4.1, +29.6]` | +34.0 `[-12.1, +80.2]` | +13.8 `[-13.1, +40.8]` |
| Next-day ≥3,000-kcal probability | +0.002 `[-0.005, +0.010]` | +0.004 `[-0.024, +0.031]` | +0.004 `[-0.013, +0.020]` |

The ≥90%-coverage sensitivity and Atwater-consistent sensitivity also fail to produce a held-out future-risk result. The Atwater sensitivity gives next-day calorie effects of +15.5, +23.4, and +15.4 kcal per 10 percentage points; only the discovery interval narrowly excludes zero, while validation and confirmation do not.

## What the FS label captures here

The largest classified sources include vegan cupcakes, several ice creams and frozen yogurts, cookies and cookie dough, Catalina dressing, orange tofu, sweet potato casserole, kunefe, and burfi. That face-valid list is reassuring, but it also demonstrates that the criterion is a nutrient-combination label, not an ultra-processing or food-category label.

Several high-calorie pizza entries have missing sugar and no exact gram weight. They cannot be classified by FS, FSOD, or CSOD without guessing. Recipe nutrient panels can also contain aggregation mistakes, sugar-alcohol accounting differences, reformulations, or user-entered database errors. Approximate MFP matches remain in the broad exposure, which is why the Atwater-consistent sensitivity is important.

## Interpretation

### Supported

- In this record, a higher share of assessable calories from foods meeting the fixed FS criterion is a reproducible marker of a higher-calorie **current day**.
- The future effect, if any, is much smaller than the same-day association.
- A broad full-HPF exposure cannot be estimated from the current files because exact food mass is too sparse.

### Not supported

- That FS exposure reliably increases the next day's extreme-intake risk.
- That the same-day association is causal.
- That ultra-processing, seed oils, sodium, sugar, fat, or a particular ingredient is the mechanism.
- That eliminating every food meeting the FS criterion is a validated weight-loss strategy for this person.

### Practical use

Do not elevate a blanket “avoid hyper-palatable foods” rule above the already stronger direct findings:

1. plan explicitly for travel and context changes;
2. use a stable, satisfying menu to reduce novelty and extreme-intake risk;
3. improve protein and fiber through calorie-neutral substitutions;
4. keep treatment decisions clinician-directed.

If the FS hypothesis is still interesting, test it prospectively as a calorie-neutral, protein/fiber-matched ABAB substitution. Replace a fixed set of common FS foods with familiar alternatives while keeping planned calories, meal timing, travel, and medication stable. The informative outcomes would be hunger, unplanned additions, adherence, and ≥3,000-kcal days—not the planned meal calories themselves.

That experiment would distinguish “the nutrient combination changes appetite” from the simpler explanation already consistent with this analysis: FS foods occur on days when more calories were chosen.
