# Long-run observed-weight trajectory audit

Date: 2026-07-30
Status: direct-observation descriptive audit plus a small internal time-ordered test

## Bottom line

The 15-year record supports a plain but useful trajectory story:

1. **The clearest loss-to-regain reversals coincide with sustained increases in logged intake.** The three clean paired reversals were approximately 1,962→2,233, 1,797→2,326, and 2,103→2,366 kcal/day. Their frequencies of days at or above 3,000 kcal rose by 7.6, 8.8, and 6.5 percentage points.
2. **Tirzepatide produced the opposite transition.** The 2023-to-pre-treatment drift averaged 2,400 kcal/day with 9.7% of days at or above 3,000 kcal. The treatment-era loss averaged 1,999 kcal/day with 1.3% high-intake days, while observed weight changed by -25.6 lb.
3. **There is no single timeless “maintenance calorie” in the direct record.** The 2019–2021 broad-maintenance interval averaged 2,175 kcal/day around 213–217 lb. Late-treatment maintenance averaged 2,005 kcal/day around 203 lb. Body size, activity, treatment, logging error, and physiology differ; the observed values are state-specific, not portable TDEE estimates.
4. **No proposed 28-day leading indicator reliably ranked the next quarter's weight change in the internal locked test.** Intake, high-intake frequency, prior travel, repetition, strength frequency, and weighing gaps all failed to preserve the expected development-period direction in 14 later pre-treatment quarters. The actionable warning is therefore a *sustained contemporaneous rise* in intake and weight trend, not a claimed one-month crystal ball.
5. **Long weigh-in gaps are a useful process alarm.** All five gaps of at least 90 days ended at a higher observed weight, with a median +11.6 lb. This does not mean not weighing causes gain; the record says the missingness was informative, plausibly because weighing was avoided during regain.

This analysis never reads P1, P2, P3, P4, P4b, modeled fat mass, or inferred TDEE.

## Reproducer and artifacts

Run:

```bash
uv run --with numpy --with pandas python analysis/EL_long_run_trajectory_audit.py
```

Artifacts:

- `analysis/EL_long_run_episodes.csv`
- `analysis/EL_quarterly_warning_anchors.csv`
- `analysis/EL_locked_warning_results.csv`
- `analysis/EL_weight_gaps.csv`
- `analysis/EL_weekly_weight_swings.csv`

Inputs are canonical daily intake, observed scale weight, checked-in travel dates, canonical strength dates, canonical steps, and a deterministic food-name repetition feature.

## The observed trajectory

The catalogue below uses predeclared observed anchors that summarize conspicuous long movements in the raw scale record. It is retrospective description, not model-estimated change-point detection.

| Episode | Observed change | Local-median change | Mean kcal/day | ≥3,000 kcal | Weight coverage |
|---|---:|---:|---:|---:|---|
| 2011-04-21→2013-09-23 long loss | -93.2 lb | -91.6 lb | 1,962 | 4.8% | 128 observations; max gap 89 d |
| 2013-09-23→2015-10-17 rebound | +24.0 | +22.8 | 2,233 | 12.5% | 16 observations; max gap 517 d |
| 2015-10-17→2016-04-08 renewed loss | -17.8 | -17.4 | 1,797 | 4.6% | 29 observations; max gap 43 d |
| 2016-04-08→2017-09-11 regain | +32.0 | +31.0 | 2,326 | 13.4% | 33 observations; max gap 99 d |
| 2017-10-06→2019-10-23 sparse regain | +14.8 | +14.8 | 2,062 | 4.3% | 3 observations; max gap 399 d |
| 2019-10-23→2021-12-21 broad maintenance | +3.8 | +4.3 | 2,175 | 4.8% | 279 observations; max gap 61 d |
| 2022-01-04→2022-06-06 loss | -13.3 | -7.3 | 2,103 | 5.8% | 99 observations; max gap 17 d |
| 2022-06-06→2023-04-01 regain | +19.7 | +14.2 | 2,366 | 12.3% | 219 observations; max gap 19 d |
| 2023-04-01→2024-09-16 pre-treatment drift | +5.3 | +6.1 | 2,400 | 9.7% | 370 observations; max gap 16 d |
| 2024-09-16→2026-02-28 treatment-era loss | -25.6 | -25.2 | 1,999 | 1.3% | 374 observations; max gap 32 d |
| 2025-12-02→2026-02-28 late maintenance | -0.7 | +0.2 | 2,005 | 1.1% | 49 observations; max gap 7 d |

“Local-median change” subtracts the median observed weight within ±7 days of each anchor. It is included to show how much a single endpoint can be moved by water, glycogen, gut contents, sodium, or measurement noise.

### What preceded, sustained, and reversed the major movements

The most repeatable direct pattern is not a special food, macro ratio, or inferred metabolic state. It is the duration of the intake regime:

- The 2015 travel spike preceded the 2015–2016 loss. The 56 days before the loss anchor averaged 2,590 kcal/day and had 41.1% high-intake days; the loss interval averaged 1,797 and 4.6%. In the 56 days after the loss endpoint, intake rose to 2,279 and the longer regain followed.
- Early-2022 loss averaged 2,103 kcal/day. The next 56 days rose to 2,395, and the subsequent regain interval averaged 2,366.
- The pre-treatment drift averaged 2,400 kcal/day. The treatment-era loss averaged 1,999, a roughly 401-kcal/day transition, while high-intake frequency fell by 8.4 percentage points.
- Late-treatment weight was stable by local medians even though intake rose from 1,898 in the preceding 56 days to 2,005 during the final interval. That is consistent with approaching current maintenance, not with continued rapid loss.

These comparisons do not hold expenditure constant and do not prove that the numerical intake differences caused every pound. They do establish that the successful regimes were sustained lower-intake regimes, not short deficits whose calories immediately returned.

The 2017–2019 sparse regain is an important counterexample to simplistic calorie thresholding: it occurred at 2,062 logged kcal/day, lower than the later 2019–2021 maintenance mean. Only three weights span the interval, body size and activity differ, and absolute logging error is unknown. A universal cutoff such as “2,100 always loses” is not supported.

## Context, repetition, and strength

### Travel remains an acute hazard, not a long-horizon oracle

Travel occupied 11.1% of the 2016–2017 regain and 8.6% of the sparse 2017–2019 regain, versus 0% of the intervening 2015–2016 loss. But it also occupied 9.0% of the treatment-era loss. Broad episode overlap therefore cannot isolate travel from treatment or era.

The stronger travel result remains the paired event study: trips averaged about +605 kcal/day versus their own prior 14 days. The present quarterly test asks a different question—whether travel in the *preceding* 28 days predicts weight over the *following* quarter—and finds no reliable signal. That null does not erase travel's contemporaneous intake effect. It says to intervene during the trip rather than claim a three-month aftereffect.

### Repetition is not a validated quarterly weight predictor

Repetition moved in inconsistent directions across episodes:

- 2015–2016 loss versus the following regain: 61.6% versus 48.1% repeated calories, consistent with the stable-menu hypothesis;
- early-2022 loss versus the following regain: 60.1% versus 62.7%, the opposite ordering;
- pre-treatment drift versus treatment loss: 67.4% versus 62.1%, also the opposite ordering.

The internal later-period test likewise gives higher repetition a positive correlation with next-quarter weight change (`r=+0.24`, AUC `0.725`), opposite the proposed risk direction. This does not overturn the replicated daily result that very high repetition marks fewer extreme-intake days. It means that daily risk marker has not translated into a reliable 91-day weight forecast.

### Strength frequency does not identify energy-balance direction

Logged strength frequency was about 1.13 sessions/week during 2016–2017 regain, 1.23 during sparse 2017–2019 regain, and 1.24 during treatment-era loss. The later-period quarterly correlation is favorable (`r=-0.36`), but the development coefficient has the opposite sign and the fitted model is worse than baseline. Strength remains worthwhile for function and lean-mass support; this dataset does not show that session count alone determines whether total weight rises or falls.

## Internal time-ordered warning test

The test uses quarter starts with observed weights within ±21 days and an observed endpoint about 91 days later. Every predictor is calculated from the 28 calendar days before the starting weight. Scaling and univariate coefficients are fitted on 17 development anchors through 2020 and then frozen for 14 anchors from 2021 through the pre-treatment boundary.

The test period is mechanically locked inside the script, but it is not epistemically pristine: this repository and several of these features have already been analyzed. Endpoint noise, only 14 later quarters, overlapping physiological timescales, and treatment exclusion also limit power.

| Candidate | Development slope, lb/SD | Later `r` with weight change | Later AUC for ≥2-lb gain | MAE | Expected direction preserved? |
|---|---:|---:|---:|---:|---|
| Development-mean baseline | 0.00 | — | 0.500 | 5.27 lb | — |
| Prior mean calories | +3.01 | -0.07 | 0.550 | 3.73 | No |
| Prior high-intake rate | +1.30 | -0.36 | 0.300 | 5.02 | No |
| Prior travel days | -0.38 | -0.16 | 0.550 | 5.16 | No |
| Prior repetition | +0.27 | +0.24 | 0.725 | 5.12 | No; sign is opposite |
| Prior strength frequency | +2.21 | -0.36 | 0.300 | 5.78 | Development sign is opposite |
| Days since previous weight | -1.40 | -0.17 | 0.600 | 4.57 | No; sign is opposite |
| Prior-90-day weight count | -1.94 | +0.25 | 0.400 | 10.27 | No |
| Equal-weight context/tracking composite | +0.02 | -0.11 | 0.400 | 5.26 | No |

The prior-calorie model's lower MAE should not be mistaken for a validated warning. Development quarters average -4.08 lb because they include the rapid 2011–2013 loss, while later quarters average +0.94 lb. All later calorie values are above the development mean, so the model mostly corrects a large era-level intercept error; within the later period its ranking correlation is `-0.07`.

The honest conclusion is negative: **none of these one-month antecedents currently earns decision-rule status for next-quarter weight change.**

## Informative missingness

There are 23 observed-weight gaps of at least 30 days:

- 14/23 end at a higher weight;
- median observed change is +2.4 lb;
- all 5 gaps of at least 90 days end higher;
- the ≥90-day median change is +11.6 lb.

The five very long gaps are 2013–2015 (+11.8 lb), spring 2015 (+6.0), winter 2016–2017 (+12.6), 2017–2018 (+11.6), and 2018–2019 (+3.2).

This is not evidence that weighing causes weight control. It is evidence against treating missing weights as random. Operationally, 30 days without a measurement—especially when daily intake remains available—should trigger a neutral “resume measurement” prompt. A model should carry an informative-missingness flag rather than silently interpolate through the gap.

## Water weight and sparse endpoints

The early-2022 episode illustrates why short-window scale loss cannot be called fat loss:

- exact-anchor change: -13.3 lb;
- ±7-day local-median change: -7.3 lb.

The following regain is +19.7 lb by exact anchors but +14.2 lb by local medians. Roughly 5–6 lb of the apparent trough/rebound magnitude is endpoint sensitivity, even though the longer directional reversal is real.

In the dense 2019–2026 record, 237 usable week-to-week changes in weekly median weight have:

- median absolute change 0.80 lb;
- 95th percentile absolute change 2.66 lb;
- only one weekly-median change at or above 5 lb.

For long episodes, exact and local-median changes generally agree closely. For sparse 2013–2015 and 2017–2019 intervals, the endpoint change is observed but the path is unknown. No interpolation should be presented as direct evidence of when the gain occurred.

## Practical translation

1. **Track regime shifts, not isolated days.** Use a 28-day intake mean, high-intake frequency, and a four-week observed-weight trend as a dashboard. Treat it as contemporaneous monitoring, not a forecast.
2. **Make a missing weight visible.** Flag 14 days as stale and 30 days as a process alert. Do not punish or moralize; the historical record suggests avoidance is itself information.
3. **Keep the travel playbook.** Travel's strongest evidence is immediate. Stabilize defaults before and during the trip; do not wait for a predicted quarterly gain.
4. **Keep strength for body composition and function.** Do not “eat back” a presumed weight-loss effect from session count.
5. **Judge loss over at least four weeks and preferably with local medians.** Short potato/fasting-like scale drops can substantially reverse through water, glycogen, and gut contents.

## What this changes

- It strengthens the claim that sustained intake environment—not a unique named diet—best distinguishes the major direct loss/regain reversals.
- It adds informative weigh-in missingness as a supported process marker.
- It does **not** validate a new predictive feature, calorie threshold, maintenance TDEE, metabolic-adaptation parameter, or causal effect of repetition, strength, or travel on long-horizon weight.
