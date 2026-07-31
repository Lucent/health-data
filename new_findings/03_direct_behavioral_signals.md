# Direct behavioral signals: travel, novelty, and repeatable menus

Date: 2026-07-30
Status: raw-data result with calendar-ordered confirmation; descriptive, not yet randomized

## Bottom line

The clearest actionable pattern in the repo is not a nutrient, genotype, temperature, or latent set point. It is **environmental disruption**:

- Travel increased intake by about 600 kcal/day relative to each trip's own preceding two weeks, across decades and after excluding the two Jordan trips.
- A highly repeatable recent menu consistently marked lower-intake days and lower next-day extreme-intake risk.
- High food novelty consistently marked higher-intake days and higher next-day extreme-intake risk.
- Domino's was a reliable same-day high-intake event, but not a stable next-day trigger.
- Potato-majority days were low-intake templates in 2022–2024, but the dataset cannot tell whether potatoes themselves, monotony, precommitment, energy density, or deliberate restriction caused that.

This points toward a practical strategy: preserve structure when the environment changes and use a small menu of satisfying default meals. It does **not** justify a potato-only diet or avoiding every novel food.

## Analysis 1: calendar-ordered food-pattern confirmation

The analysis uses only `intake_daily.csv`, `intake_foods.csv`, and checked-in travel dates. It does not read P2, P3, P4, modeled TDEE, fat mass, set point, or modeled drug levels.

The periods were fixed before the script was run:

- Discovery: 2018-03-01 through 2021-12-31
- Validation: 2022-01-01 through 2022-12-31
- Confirmation: 2023-01-01 through 2024-09-16, ending before tirzepatide

Generic food-pattern features were calculated from normalized diary strings. Repetition is the calorie share from foods seen in the preceding 14 days; novelty is the calorie share not seen in the preceding 30 days. Partial correlations control for current-day intake, the prior 28-day intake mean, weekday, season, and linear time. This is a deliberately small, pre-specified feature set.

Command:

```bash
uv run --with numpy --with pandas python analysis/EB_direct_future_intake.py
```

Artifacts: `analysis/EB_direct_daily_features.csv`, `analysis/EB_direct_future_intake_results.csv`, `analysis/EB_direct_event_summaries.csv`, `analysis/EB_direct_event_year_summaries.csv`

### Repetition is a stable marker of lower intake

Days on which at least 90% of calories came from foods seen in the preceding 14 days were lower than other days in every split:

| Period | Event days | Same-day difference | Next-day difference | Next-day ≥3,000 kcal |
|---|---:|---:|---:|---:|
| Discovery | 162 | -211 kcal | -112 kcal | 1.9% vs 5.2% |
| Validation | 44 | -431 kcal | -298 kcal | 4.5% vs 11.2% |
| Confirmation | 114 | -338 kcal | -124 kcal | 2.6% vs 11.2% |

The same-day difference was negative in all seven calendar years from 2018–2024. The next-day difference was also negative in all seven.

After removing travel days and days immediately before travel, the association with **mean** next-day calories mostly disappears. The association with next-day ≥3,000-kcal intake remains in the same direction in all three splits: partial `r = -0.070, -0.051, -0.055`.

Interpretation: repetition is not a powerful continuous appetite drug. It is a modest, replicating marker of protection against extreme days.

### Novelty is the mirror image

Days on which at least 50% of calories came from foods not seen in the preceding 30 days were higher-intake days:

| Period | Event days | Same-day difference | Next-day difference | Next-day ≥3,000 kcal |
|---|---:|---:|---:|---:|
| Discovery | 380 | +168 kcal | +77 kcal | 8.2% vs 3.6% |
| Validation | 75 | +245 kcal | +113 kcal | 14.7% vs 9.3% |
| Confirmation | 103 | +461 kcal | +162 kcal | 19.4% vs 7.7% |

The same-day difference was positive in all seven years; the next-day difference was positive in six of seven and essentially zero in 2020.

At home, novelty no longer predicts the next day's mean intake, but it continues to mark next-day extreme-intake risk in every split: partial `r = +0.056, +0.027, +0.064`.

Interpretation: much of the calorie magnitude belongs to travel and context change, but a smaller extreme-intake signal remains at home. Novelty probably bundles availability, ordering, reward, social setting, and loss of precommitment. The data does not isolate sensory-specific satiety.

### Specific-food results are weaker

- **Domino's:** +404, +492, and +346 kcal on the day across the three splits. Next-day differences were +108, +479, and -127 kcal. It is a high-calorie event, not a replicating future trigger.
- **Potato-majority days:** only one occurred in discovery, so this cannot pass the locked design. In validation and confirmation, 48 combined days averaged about 700 kcal less than other days and had almost no next-day ≥3,000-kcal events. The template is promising; the potato-specific mechanism is unidentified.
- **Dessert ≥25% of calories:** changed sign across periods. No stable protective or harmful conclusion survives.
- **Menu entropy, number of foods, and largest-food share:** small and inconsistent future-intake associations. The simple repetition/novelty history is more informative.

## Analysis 2: paired travel event study

Forty-one trips lasting at least two days had complete 14-day pre/post intake windows. Each trip was compared with its own preceding 14 days.

Command:

```bash
uv run --with numpy --with pandas python analysis/EC_travel_event_study.py
```

Artifacts: `analysis/EC_travel_trip_results.csv`, `analysis/EC_travel_summary.csv`

| Sensitivity | Trips | Travel days | Mean change vs pre-trip | Median change | Trips higher / lower | Exact sign-test p |
|---|---:|---:|---:|---:|---:|---:|
| All eligible | 41 | 271 | +605 kcal/day | +632 | 36 / 5 | 0.00000078 |
| Pre-tirzepatide | 37 | 223 | +639 | +659 | 33 / 4 | 0.0000011 |
| Pre-tirzepatide, excluding Jordan | 35 | 163 | +615 | +659 | 31 / 4 | 0.0000035 |
| Modern MHTML pre-tirzepatide | 14 | 86 | +622 | +708 | 13 / 1 | 0.0018 |
| Short trips, 2–6 days | 29 | 105 | +664 | +692 | 26 / 3 | 0.000015 |
| Long trips, 7+ days | 12 | 166 | +465 | +360 | 10 / 2 | 0.039 |
| Tirzepatide era | 4 | 48 | +299 | +293 | 3 / 1 | 0.625 |

The mean post-trip difference versus the pre-trip baseline was -48 kcal/day overall and -62 kcal/day pre-tirzepatide. That is too small to compensate for the during-trip excess, although overlapping windows and changing baselines prevent exact energy accounting.

This is not just Jordan. It is not just the old OXPS source. It appears in short trips, recent pre-drug trips, and most locations. Tirzepatide may attenuate it, but four trips are far too few for a treatment interaction claim.

Observed scale-weight changes are reported in the artifact only as anchors. They should not be interpreted as fat gain because weigh-ins can be up to 14 days from the trip and are heavily affected by water, sodium, gut contents, and missing-on-trip measurement.

## What this changes about the weight-loss plan

### 1. Treat travel as a planned intervention, not an exception

A plausible target is avoiding part of the historical +600 kcal/day excursion, not “burning it off.” Before each trip:

- choose a default breakfast and one or two default high-satiety meals that are easy to obtain at the destination;
- identify groceries or delivery options before arrival;
- pre-log or pre-decide the first 24 hours, when the normal food environment disappears;
- allow planned novelty rather than making every meal a novelty event;
- keep the prescribed medication schedule unchanged unless the prescriber says otherwise;
- continue logging hunger, GI effects, and intake without treating a high day as failure.

The exact template should be vegetarian, protein-adequate, tolerable on the medication, and enjoyable enough to repeat.

### 2. Build a small “boring but good” home menu

Use a rotating set of satisfying defaults, not a single-food diet. The result suggests that menu stability helps prevent extremes, but it does not say less variety is nutritionally superior.

A clean prospective test would alternate predeclared two-week blocks:

- **Stable-menu block:** a small set of repeatable meals, with calorie target and protein kept ordinary.
- **Usual-choice block:** normal variety.

Do not deliberately starve during the stable block. Measure hunger, daily intake, ≥3,000-kcal events, eating friction, GI tolerance, and adherence. Keep medication and travel constant. An ABAB sequence is deterministic and avoids pretending that post hoc potato periods were randomized.

### 3. Use potatoes as one template component, not the theory

Potato periods demonstrate that a highly structured vegetarian menu can be filling and low-intake. They do not prove a unique metabolic potato effect, and the historical versions were sometimes low in protein or followed by rebound. A safer test is a nutritionally complete repeatable menu that may include potatoes.

### 4. Stop blaming isolated ingredients for context bundles

Domino's days combine ordering, palatability, energy density, sodium, refined carbohydrate, fat, portion size, and occasion. The current data cannot separate those. The high-leverage move is changing the ordering template or portion/environment, not asserting that cheese, seed oil, carbohydrate, or sodium is the sole cause.

## Remaining uncertainty

- The travel table is derived from an unavailable Google Timeline export, so dates can be audited only against the checked-in Markdown.
- Food-name novelty can reflect renamed diary entries as well as genuinely new foods. Year-level stability and the home-only sensitivity reduce, but do not eliminate, this concern.
- Repetition and novelty are partly chosen in response to appetite or planned restriction. The next-day and calendar-ordered results reduce same-day reverse causality but remain observational.
- The cutoffs of 90% repetition and 50% novelty are interpretable event definitions, not biological thresholds.
- The decisive next step is prospective replication, not more cut-point searching.
