# Sleep, activity, and circadian re-audit

Date: 2026-07-30

## Bottom line

Sleep and activity matter for health, but neither is currently the highest-leverage explanation for this subject's historical weight pattern.

- The reliable sleep regime is 2022 onward, not the full 2016–2026 file. In that regime, 1,504 of 1,520 calendar days have a plausible main-sleep record, mean sleep is 8.03 hours, and median sleep is 8.00 hours. Only six nights are under five hours.
- Sleeping longer than the subject's own prior-28-day mean has a small association with lower subsequent intake, not a true null: about −38 kcal the next day and −36 kcal/day over the following seven days per additional hour. Partial correlations are only −0.071 and −0.104, the effect does not persist in the tirzepatide era, and sleep features do not materially improve temporal holdout prediction.
- Later bedtime does not predict higher intake here. The small association actually points toward lower intake, which is likely schedule/context confounding and is not a reason to delay sleep.
- Once current intake is controlled, step deviations do not independently predict meaningful future intake across 2022–2026. The pre-tirzepatide coefficients are essentially zero. A small positive association appears during tirzepatide treatment, but it is not causal and does not survive as a useful general predictor.
- Exact observed-weight pairs show no credible beneficial step effect over 1, 7, or 14 days. The 28-day coefficient points the wrong direction. This is a noisy observational test, not evidence that exercise is useless.
- Walking and resistance exercise should remain in the plan for fitness, health, function, and lean-mass support. The repo does not justify treating walking as a precisely quantified RMR drug, an appetite suppressant, or a license to eat back estimated activity calories.

The action change is therefore narrow: preserve the already-adequate sleep routine; do not spend scarce adherence effort forcing an earlier clock time; use exercise for its established benefits; and prospectively watch for extra hunger/intake after unusually active days rather than assuming either perfect compensation or zero compensation.

## Why the legacy results are not reliable

### AE: “sleep predicts nothing”

[`AE_sleep_null.py`](../analysis/AE_sleep_null.py) anchors the table on recorded sleep rows and then uses `shift(-1)` for “next-day” calories, steps, and weight. Before 2022, sleep dates are sparse. The next row can be many calendar days later, so this is not an exact next-day join. Its rolling windows also count sleep records rather than calendar days.

AE includes all 2,057 exported rows. Thirty-two are under two hours, all before 2022, and many are plainly daytime naps. Those rows account for more than half of its 57 “under five hour” observations. In the near-complete 2022–2026 regime, only six nights are under five hours. The published extreme-sleep comparison is therefore mostly a source-regime artifact.

The corrected statement is:

> In a subject who usually gets about eight hours of recorded sleep, ordinary night-to-night duration variation has little predictive value for intake, steps, or observed weight. The repository cannot estimate the effect of chronic sleep restriction because that exposure is rare in the reliable regime.

### BS: bedtime timing and sleep debt

[`BS_sleep_timing_split.py`](../analysis/BS_sleep_timing_split.py) inherits the nap problem, uses record-count rolling windows, and uses the same non-calendar `shift(-1)` outcomes. It also searches 25 bedtime cutoffs across multiple outcomes and then reports the largest absolute partial correlation without a held-out selection period.

Its TDEE endpoints come from the quarantined P4 latent-state model, not a direct metabolic measurement. The reported 30-day bedtime-to-TDEE partial correlation of about −0.10 is therefore exploratory twice over: cutoff/window search plus model-derived outcome.

The reliable post-2022 schedule is:

| Metric | Mean | Median | SD |
|---|---:|---:|---:|
| Sleep duration | 8.03 h | 8.00 h | 1.18 h |
| Bedtime, night clock | 02:08 | 02:10 | 1.29 h |
| Wake time | 10:10 | 10:15 | 1.04 h |

Weekend sleep is about 0.84 hours longer, composed of an approximately 0.65-hour earlier bedtime and 0.19-hour later wake time. This is a real schedule pattern. It is not evidence that the weekend itself changes metabolism.

### U and V: future activity outcomes

[`U_steps_compensation.py`](../analysis/U_steps_compensation.py) and [`V_exercise_walk_analysis.py`](../analysis/V_exercise_walk_analysis.py) compute forward outcomes as:

```python
series.shift(-1).rolling(horizon).mean()
```

On a forward-sorted daily table, a nominal 14-day value at date *t* averages roughly *t−12 through t+1*. It is mostly historical, not future. The V step-matched comparison also matches against self-selected “other” days without time, season, treatment-era, or reuse safeguards.

Consequences:

- The “high steps predict future intake/TDEE” result in U is not a valid future-outcome test.
- The paired-daylight-walk future-14-day contrasts in V are not future-14-day contrasts.
- These artifacts should not support constrained-expenditure, gravitostat, sunlight, or appetite-suppression conclusions.

### AK: “sunlight exposure”

[`sunlight.py`](../steps-sleep/sunlight.py) does not measure exposure. It calculates the maximum daylight available between recorded wake time and sunset, assuming the subject remains awake and potentially exposed until sunset. It has no indoor/outdoor, illuminance, weather, shade, or screen data.

[`AK_sunlight_exposure.py`](../analysis/AK_sunlight_exposure.py) also uses next-record `shift(-1)` outcomes. Its RMR section has only a small, clustered calorimetry sample and selects among related walk/sunlight models. The exact `+14 kcal/day per walk session` claim remains exploratory and is not a sound exercise prescription.

## Corrected source and regime audit

| Source | Reliable analysis boundary | Result |
|---|---|---|
| Sleep | 2022-01-01 to 2026-02-28 | 1,504/1,520 dates, 98.9% |
| Main-sleep plausibility | 2–14 hours | All 1,504 post-2022 rows retained |
| Early nap-like episodes | Full sleep file | 32/2,057 under 2 h; all pre-2022 |
| Steps | 2022-01-01 to 2026-02-28 | 1,520/1,520 dates, all from Samsung |
| Very low steps | Post-2022 | 14 days under 100 steps, likely partial/non-wear |
| Exercise reconciliation | Post-2022 exercise days with session step count | Session steps exceed the daily total on 4/404 days |
| Observed weight | Post-2022 | 1,059/1,520 dates, 69.7%; no interpolation used |

“Samsung source” is not identical to “complete wear.” Step analyses should retain a wear-quality flag. The sharp decline in median steps during 2025–2026 may contain real behavior change, watch/phone carry change, or both.

## Corrected direct results

The primary analysis uses a complete calendar index, exact future dates, direct calories/steps/observed weight, and deviations from each variable's prior-28-day mean. Calendar controls are weekday, smooth annual seasonality, time trend, and tirzepatide era. No P2/P3/P4 outcomes are used.

### Sleep and timing

| Exposure | Outcome | n | Coefficient | Partial r |
|---|---|---:|---:|---:|
| +1 h sleep vs local baseline | Next-day intake deviation | 1,489 | −38 kcal | −0.071 |
| +1 h sleep vs local baseline | Following 7-day intake deviation | 1,489 | −36 kcal/day | −0.104 |
| +1 h sleep vs local baseline | Following 14-day intake deviation | 1,489 | −33 kcal/day | −0.111 |
| +1 h later bedtime vs local baseline | Next-day intake deviation | 1,489 | −26 kcal | −0.051 |
| +1 h later bedtime vs local baseline | Following 7-day intake deviation | 1,489 | −31 kcal/day | −0.096 |

The sleep-intake direction is compatible with the general sleep-extension literature, but this observational effect is small, may reflect illness/work/travel/restriction contexts, and is regime-dependent:

- Before tirzepatide, +1 hour of sleep predicts about −57 kcal/day over the following week, partial r −0.147.
- During tirzepatide, the corresponding coefficient is +8 kcal/day, partial r +0.033.

Sleep does not predict subsequent steps in a stable direction. Exact observed-weight coefficients are small and inconsistent: −0.08 lb over 7 days, −0.16 lb over 14 days, and −0.08 lb over 28 days per extra hour relative to baseline, with absolute partial correlations at or below 0.068. These are not decision-grade effects.

### Steps and exercise

After controlling for current intake deviation, sleep, bedtime, calendar, time, and treatment era:

| Exposure | Outcome | n | Coefficient | Partial r |
|---|---|---:|---:|---:|
| +1,000 current-day steps vs local baseline | Next-day intake | 1,489 | +6 kcal | +0.041 |
| +1,000 current-day steps vs local baseline | Following 7-day intake | 1,489 | +0 kcal/day | +0.003 |
| +1,000/day in recent 7-day steps vs preceding baseline | Next-day intake | 1,476 | +7 kcal | +0.028 |
| +1,000/day in recent 7-day steps vs preceding baseline | Following 7-day intake | 1,476 | +3 kcal/day | +0.020 |

The pre-tirzepatide coefficients are essentially zero. During tirzepatide, current-day step deviation is associated with about +15 kcal/day per 1,000 steps over the following three days and +10 kcal/day over the following week; recent seven-day step deviation gives +22 and +13 kcal/day, respectively. Even if causal, these are partial rather than complete compensation estimates. More likely explanations include active social/travel days, planned higher-food days, non-wear structure, and residual time-varying behavior.

Recorded exercise sessions do not show exercise-induced overeating. After step adjustment, a walk session is associated with −36 next-day kcal across the full period but essentially 0 kcal during tirzepatide. Partial correlations are at or below 0.054. This also is not evidence that walking pharmacologically suppresses appetite.

For exact observed-weight pairs, recent steps have coefficients of −0.02, −0.00, −0.00, and +0.09 lb per extra 1,000 steps/day over 1, 7, 14, and 28 days. The signs and magnitudes do not form a credible weight-loss dose-response.

## Temporal holdout

Models were fit through 2024-12-31 and scored without refitting on 416 dates from 2025-01-01 through 2026-02-28. The baseline includes calendar terms, time, tirzepatide era, and current intake deviation.

| Outcome | Model | Test MAE | Test RMSE |
|---|---|---:|---:|
| Next-day calorie deviation | Baseline | 283.3 | 365.1 |
| Next-day calorie deviation | + sleep duration and bedtime | 282.3 | 363.7 |
| Next-day calorie deviation | + current steps | 283.0 | 364.5 |
| Next-day calorie deviation | + recent 7-day steps | 283.2 | 364.8 |
| Next-day calorie deviation | Full sleep + activity | 281.9 | 362.9 |
| Following 7-day calorie deviation | Baseline | 141.3 | 180.9 |
| Following 7-day calorie deviation | + sleep duration and bedtime | 143.5 | 183.9 |
| Following 7-day calorie deviation | + current steps | 141.5 | 181.2 |
| Following 7-day calorie deviation | + recent 7-day steps | 140.8 | 180.9 |
| Following 7-day calorie deviation | Full sleep + activity | 144.0 | 184.8 |

The best next-day improvement is 1.4 kcal MAE, less than 0.5%. At seven days, the full feature set is worse. Sleep and activity do not add useful predictive power for intake in the current feature formulation.

## What should change now

1. **Maintain, do not optimize, the sleep duration.** Around eight hours is already adequate for most nights. Preventing erosion is sensible; making sleep the central weight-loss project is not supported.
2. **Do not infer that a 02:00 bedtime is metabolically harmful from this repo.** If the schedule is socially workable and consistent, forced phase-advancement is unlikely to beat the established food-context and treatment levers.
3. **Do not infer that the sleep record rules out sleep pathology.** Samsung duration does not diagnose sleep apnea, fragmentation, hypoxemia, or sleep quality. Loud snoring, witnessed pauses, morning headaches, excessive sleepiness, or hypertension would justify clinical screening regardless of these correlations.
4. **Keep walking and resistance training.** Their strongest rationale here is health, conditioning, function, and lean-mass preservation—not the repo's exact RMR coefficient.
5. **Do not automatically eat back activity calories.** There is no evidence of full compensation. During unusually active treatment-era periods, pre-plan a satisfying protein- and fiber-forward meal and observe hunger rather than improvising.
6. **Treat very low step totals as missing/partial wear until reconciled.** A zero is not sedentary truth if an exercise session exists on the same date.

## Prospective measurement plan

Use fixed definitions before looking at outcomes:

| Domain | Metric | Prespecified definition |
|---|---|---|
| Weight | Weekly observed weight | Median of at least 4 standardized morning weights; retain individual readings |
| Intake | Weekly intake | Seven exact calendar days; report mean, median, and ≥3,000-kcal days |
| Sleep quantity | Main sleep | Longest plausible 2–14 h episode assigned to wake date; flag missing |
| Sleep regularity | Bedtime variability | Night-clock bedtime SD over seven days; do not use coefficient of variation |
| Short sleep | Event | Under 6 h, chosen prospectively because under-5 h events are too rare |
| Steps | Daily movement | Samsung total plus wear-quality flag; flag under 100 steps or exercise/daily mismatch |
| Exercise | Deliberate activity | Sessions and minutes by type, separate from daily steps |
| Hunger | Optional direct outcome | 0–10 rating before dinner and next morning, recorded before seeing calorie totals |
| Context | Confound flags | Travel, illness, restaurant/order-in, major schedule disruption, and medication change |

For an activity-compensation check:

- Define a high-activity day before analysis as at least 3,000 steps above the preceding 28-day mean.
- Primary outcome: exact mean calorie deviation on days +1 through +3.
- Secondary outcomes: hunger next morning and days +1 through +7 calories.
- Do not repeatedly change the threshold or horizon.
- Accumulate at least 20 high-activity events and 20 eligible ordinary days in each treatment regime before interpreting.
- Report the event dates and wear-quality exclusions, not only a regression coefficient.

For sleep, do not deliberately restrict sleep. Treat naturally occurring under-6-hour nights as events and use the same exact-date outcome logic. A useful prospective effect would need to be large enough to change behavior, not merely achieve a small correlation.

## Reproducer and artifacts

Run:

```bash
uv run --with pandas --with numpy python analysis/EI_sleep_activity_reaudit.py
```

Artifacts:

- [`EI_sleep_activity_quality.csv`](../analysis/EI_sleep_activity_quality.csv)
- [`EI_sleep_activity_results.csv`](../analysis/EI_sleep_activity_results.csv)
- [`EI_sleep_activity_holdout.csv`](../analysis/EI_sleep_activity_holdout.csv)

The script is deterministic, takes no arguments, and does not use latent-state outcomes.
