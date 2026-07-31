# Calendar and extreme-intake clustering audit

Date: 2026-07-30
Status: direct-observation, calendar-ordered audit; descriptive and predictive, not causal

## Bottom line

Days at or above 3,000 kcal were not independent one-off events before
tirzepatide. They often continued into the next day, including at home:

- across the full pre-treatment home record, the next day was also extreme
  after 22.2% of extreme days versus 5.4% of other days;
- after weekday, season, time, prior-28-day intake, and travel adjustment, an
  extreme day marked another **+359 kcal the next day**
  `[+266, +452]` and **+16.2 percentage points** of next-day extreme risk
  `[+10.9, +21.5]`;
- 176 of 391 extreme days occurred in multi-day runs, and 47 of 68 multi-day
  episodes contained no checked-in travel day.

The continuation signal was weaker in the middle validation period, so this is
not a deterministic “binge begets binge” rule. There is also no stable
multi-day compensation rule: days 2–7 after an extreme day were above the
recent baseline in discovery and confirmation but below it in validation.

Simple recent calorie history did provide useful one-day-ahead **ranking**
information in both held-out pre-treatment periods. The useful signal was not
the binary extreme label alone; it was whether yesterday and the prior week
were running above the preceding 28-day level. Calibration moved across eras,
so the model is a warning-state prototype, not a calorie prescription or a
deployment-ready risk score.

Weekend intake was somewhat higher in the long pre-treatment average, but the
effect was not stable across ordered periods. The very large treatment-era
weekday pattern cannot be interpreted as a weekend mechanism because injection
day was almost perfectly tied to weekday in that era.

## Reproducer and evidence boundary

Command:

```bash
uv run --with numpy --with pandas python analysis/EM_calendar_extreme_clustering.py
```

Artifacts:

- `analysis/EM_calendar_dayofweek.csv`
- `analysis/EM_epoch_extreme_summary.csv`
- `analysis/EM_extreme_transitions.csv`
- `analysis/EM_extreme_runs.csv`
- `analysis/EM_forward_prediction.csv`
- `analysis/EM_daily_features.csv`

The script reads only:

- `intake/intake_daily.csv`;
- the overlapping annotations in `intake/diet_epochs.csv`;
- the checked-in date ranges in `travel/trips.md`.

It does not read food features, body weight, P1/P2/P3/P4/P4b, modeled fat mass,
modeled TDEE, set point, RMR, sleep, steps, or modeled drug levels. “Extreme”
means only the prespecified diary threshold of at least 3,000 kcal; it is not a
clinical binge diagnosis.

The ledger contains 5,429 consecutive days and 391 extreme days (7.2%). There
were 384 extreme days among 4,899 pre-treatment days (7.8%) and seven among 530
treatment-era days (1.3%).

## Calendar-ordered design

The complete ledger permits long, non-overlapping calendar periods:

| Period | Dates | Role |
|---|---|---|
| Discovery | 2011-04-20–2017-12-31 | Fit preprocessing and prediction coefficients |
| Validation | 2018-01-01–2021-12-31 | First untouched pre-treatment check |
| Confirmation | 2022-01-01–2024-09-16 | Second untouched pre-treatment check |
| Treatment shift | 2024-09-17–2026-02-28 | Descriptive regime-shift audit only |

One-day-ahead features for target day `t` use only days through `t-1`.
Transition analyses use exposure on day `t` and outcomes beginning at `t+1`.
The six-day later outcome uses only `t+2..t+7`. No exposure includes its own
outcome day.

## Day of week: a weak and unstable pre-treatment cue

Across all pre-treatment home days, weekends averaged 2,168 versus 2,118 kcal
on weekdays, a difference of +51 kcal/day. Extreme-day rates were 8.5% versus
5.6%. Within-calendar-week centering gives the same direction.

The ordered periods show why the long-run average should not be promoted to a
strong rule:

| Home-only period | Weekend minus weekday kcal | Weekend extreme rate | Weekday extreme rate |
|---|---:|---:|---:|
| Discovery | +62 | 10.4% | 5.8% |
| Validation | +4 | 4.3% | 3.8% |
| Confirmation | +92 | 9.9% | 7.9% |
| Treatment shift | +208 | 1.5% | 0.0% |

Saturday was the clearest high day in discovery but not in validation.
Confirmation again had higher weekends, but the difference was modest relative
to daily variation.

The treatment-era average is dominated by a different pattern: home intake
averaged about 1,626 kcal on Wednesday and 1,778 on Thursday, versus 2,110 on
Saturday and 2,167 on Sunday. The direct medication audit already established
that injection cycle day mapped to weekday with 99.1% purity. Weekday,
medication cycle, calendar time, and behavior therefore cannot be separated
here. The treatment-era pattern is real as a description but is not evidence
for a biological weekend effect or a dosing change.

Travel remains a much larger context than weekend. Checked-in travel days
averaged 2,843 kcal on weekends versus 2,599 on weekdays, with extreme rates of
40.4% and 27.9%. At home, the pre-treatment weekend difference was only 51
kcal/day.

## Extreme days cluster, but the strength varies

There were 283 extreme-intake episodes:

| Consecutive run length | Episodes | Extreme days represented |
|---|---:|---:|
| 1 day | 215 | 215 |
| 2 days | 47 | 94 |
| 3 days | 13 | 39 |
| 4–7 days | 8 | 43 |

The longest run was seven days. All seven treatment-era extreme days were
isolated; all multi-day runs occurred before treatment.

The overlapping named epochs help locate, but not explain, the clusters.
`travel_binge_1` had 13 extreme days in 23 days (56.5%), all part of multi-day
runs. But clustering was not confined to a named travel period:
`pre_tirzepatide_baseline_1` had 24 extreme days in 184 days (13.0%), 15 of
them clustered, while only six epoch days were checked-in travel days.
`sweet_earth_era` had 24 extreme days in 153 days (15.7%), 16 clustered. In
contrast, the 530-day tirzepatide annotation had seven extreme days (1.3%) and
none clustered.

These are overlapping, manually declared context annotations, not a complete
causal partition. In particular, the epoch names cannot establish that a food
brand, macronutrient pattern, or medication caused the rates.

### Exact next-day continuation at home

| Period | Extreme anchors | Extreme tomorrow after extreme | Extreme tomorrow after other day | Adjusted next-day kcal effect |
|---|---:|---:|---:|---:|
| Discovery | 156 | 25.6% | 5.7% | +425 `[+288, +562]` |
| Validation | 55 | 9.1% | 3.8% | +152 `[-31, +335]` |
| Confirmation | 82 | 24.4% | 7.1% | +361 `[+208, +514]` |
| Treatment shift | 2 | 0.0% | 0.4% | -70 `[-176, +36]` |

The adjusted estimates use split-specific ordinary least squares with
Newey–West standard errors at seven lags. Controls are the prior-28-day mean,
target weekday, season, linear time, travel today, and travel tomorrow.

The direction replicated in all three pre-treatment periods, but the
validation interval includes zero and the magnitude is much smaller. This
supports a continuation-prone state, not a fixed fourfold causal multiplier.

### Compensation is not reliable

Relative to the trailing 28-day baseline, the next day after a home extreme was
still higher by about +303 kcal in discovery, +99 in validation, and +303 in
confirmation. That is inconsistent with immediate full compensation.

For days 2–7, however, the corresponding differences were +71, -68, and +132
kcal/day. The middle period shows partial compensation while the other two show
continued elevation. The repository therefore does not support either of
these blanket claims:

- “one high day is automatically offset over the next week”;
- “every high day causes a week-long cascade.”

## What can be predicted one day ahead?

Two discovery-fitted models were compared:

1. **Context only:** prior-28-day mean, target weekday, season, time, and known
   checked-in travel status.
2. **Context plus calorie history:** context plus yesterday's deviation from
   the prior-28-day mean and the prior-seven-day mean's deviation from that
   baseline.

All preprocessing and coefficients were fit in discovery only.

### All days

| Split | Context AUC | + calorie-history AUC | Context Brier | + history Brier |
|---|---:|---:|---:|---:|
| Discovery | 0.721 | 0.806 | 0.0727 | 0.0696 |
| Validation | 0.645 | 0.700 | 0.0469 | 0.0467 |
| Confirmation | 0.634 | 0.775 | 0.0830 | 0.0787 |

### Home only

| Split | Context AUC | + calorie-history AUC | Context Brier | + history Brier |
|---|---:|---:|---:|---:|
| Discovery | 0.659 | 0.771 | 0.0644 | 0.0614 |
| Validation | 0.588 | 0.665 | 0.0402 | 0.0409 |
| Confirmation | 0.573 | 0.746 | 0.0794 | 0.0747 |

The ranking improvement replicated. Probability calibration did not: the
home-only history model predicted 7.8% in validation when 3.9% was observed,
and its Brier score was slightly worse than context alone. In confirmation it
predicted 11.7% versus 8.5% observed and improved Brier score.

Using the discovery-locked 90th-percentile risk cutoff, home extreme rates were:

| Split | Flagged days | Flagged extreme rate | Other-days extreme rate |
|---|---:|---:|---:|
| Discovery | 229 | 24.5% | 5.1% |
| Validation | 132 | 7.6% | 3.6% |
| Confirmation | 228 | 20.2% | 4.9% |

The flagged fraction expanded to 23.6% in confirmation because the calorie
distribution shifted. This is another calibration warning, not a reason to
retune the cutoff after seeing confirmation.

Adding only binary extreme-history variables did not generalize as well. Its
home AUCs were 0.720, 0.584, and 0.678; in validation it was slightly worse
than context alone. The more portable signal is **recent calorie momentum**,
not yesterday's threshold label by itself.

## Interpretation and practical use

### Directly supported, descriptive

- Pre-treatment extreme days frequently clustered, including outside travel.
- The next day after an extreme day was usually above the recent baseline, not
  fully compensated.
- A week running above its recent calorie baseline carried useful
  one-day-ahead ranking information in both held-out periods.
- Treatment coincided with a major collapse in both event frequency and
  persistence.

### Not established

- That one extreme day causes the next; permissive context, food availability,
  stress, intentional restriction cycles, and latent appetite state can persist
  across both days.
- That weekend itself is a strong causal driver.
- That a prediction flag should trigger fasting, medication changes, or
  punitive exercise.
- That the discovery probabilities remain calibrated after treatment or a
  major lifestyle shift.

### Low-risk translation

The data favor a **same-day interruption plan** over waiting for automatic
compensation. When intake or unplanned additions are already running above the
recent baseline:

1. keep the next meal ordinary rather than “making the day count”;
2. reinstate the familiar protein-forward default at the next eating
   opportunity;
3. remove leftover delivery/party food from immediate availability when
   practical;
4. preserve normal medication and exercise plans;
5. avoid a compensatory fast, since the weekend-fasting audit found a separate
   rebound pattern.

This is best treated as an implementation cue for the stable-menu and travel
plans, not as a new diet theory. A prospective log could record one binary
field—“above-plan momentum noticed before the final meal”—and whether the
ordinary-default reset was used. That would test whether an early cue can
interrupt continuation without turning a retrospective risk score into a
self-fulfilling restriction cycle.
