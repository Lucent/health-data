# Temporal semantics audit of legacy analyses

Date: 2026-07-30

## Bottom line

A repository-wide static scan and targeted semantic recheck found one newly
materially invalidated promoted claim:

> The repository does **not** show that higher day-to-day intake variance is
> protective, or less fattening, at the same contemporaneous calorie surplus.

The legacy AF analysis compares intake mean, intake variance, and modeled TDEE
from the *preceding* 14 or 30 days with modeled fat-mass change over the
*following* 14 or 30 days. Those intervals do not overlap. When fat-mass change
is aligned to the same interval as the intake measurements, the association
vanishes:

| Horizon | Legacy non-overlapping partial r | Same-window partial r |
|---|---:|---:|
| 14 days | -0.204 | +0.003 |
| 30 days | -0.123 | -0.006 |

Both versions still use legacy P4 modeled fat mass and TDEE, so the aligned null
is not a definitive biological estimate. It is sufficient to show that AF's
code does not support the stated causal interpretation.

A second high-severity defect affects an unpromoted diagnostic:
`CL_revealed_sp_from_binges.py` constructs its purported current predictor from
a centered 30- or 90-day binge-rate window and then “predicts” the next seven
days of binges. Every one of those seven outcome days is already inside the
predictor. That diagnostic is direct target leakage and must not be used.

The audit also confirms the already-documented sleep and U/V future-window
errors, corrects the AB next-three-day horizon, and inventories several
retrospective centered transforms that are acceptable for description only.

## Scope and method

The scan covered all 114 Python analysis files that existed before this audit.
It searched for negative shifts, centered rolling windows, filtered-table row
shifts, forward-window constructions, and claims framed as prediction in
`README.md`, `FINDINGS.md`, `ROADMAP.md`, and `RECOMMENDATIONS.md`. Twenty-four
candidate scripts contained a negative shift or centered window. Plot-only
transforms, correctly aligned windows, and the already-documented P3/P4 defects
were separated from material claim paths.

The deterministic reproducer locates every audited expression in its current
source file, exports the live line number and source text, and fails if a
declared expression disappears:

```bash
uv run --with numpy --with pandas python analysis/EK_temporal_semantics_audit.py
```

Artifacts:

- `analysis/EK_temporal_semantics_inventory.csv`
- `analysis/EK_temporal_semantics_checks.csv`

The quantitative comparisons use the current repository state. They are
semantic checks, not attempts to rehabilitate P4 physiology.

## High-severity paths

| Script | Implemented time semantics | Impact |
|---|---|---|
| `AF_intake_variance.py:37-42` | Trailing intake/TDEE at `t-w+1..t`; FM outcome at `t..t+w` | Invalidates “variance is protective at the same surplus” |
| `CL_revealed_sp_from_binges.py:118,142` | Centered predictor includes future binge outcomes | Invalidates CL's next-seven-day diagnostic |
| `U_steps_compensation.py:47-48` | `shift(-1).rolling(h)` is mostly trailing | Already invalidated in report 09 |
| `V_exercise_walk_analysis.py:238-243` | Same reversed window, on available step rows | Already invalidated in report 09 |
| `AE_sleep_null.py:42-44` | Next sleep record, not necessarily next date | Already corrected in report 09 |
| `BS_sleep_timing_split.py:122-123` | Next sleep record, not necessarily next date | Already corrected in report 09 |

The rolling-window algebra is exact:

- `shift(-1).rolling(7)` at date *t* uses *t-5* through *t+1*: five prior
  days, today, and one future day.
- `shift(-1).rolling(14)` uses *t-12* through *t+1*.
- `rolling(3).shift(-1)` uses *t-1*, *t*, and *t+1*.
- In contrast, `rolling(7).shift(-7)` uses *t+1* through *t+7* on a complete
  daily index and is temporally correct. The scanner therefore does not treat
  every negative shift as an error.

## AF: the protective-variance result disappears when intervals align

The legacy code computes:

```python
cal_std[t] = std(calories[t-w+1:t])
surplus[t] = mean(calories[t-w+1:t] - tdee[t-w+1:t])
fm_change[t] = fm[t+w] - fm[t]
```

This can ask whether earlier intake variability predicts a later period, but it
cannot ask whether variability during a calorie-surplus window changes fat gain
*during that same window*. The prose repeatedly makes the latter claim.

Using the same P4 values solely to isolate this indexing difference:

| Horizon/outcome | Partial r, controlling modeled surplus | Variance coefficient |
|---|---:|---:|
| 14d future, legacy | -0.2042 | -0.001272 lb per kcal SD |
| 14d same window | +0.0029 | +0.000004 |
| 30d future, legacy | -0.1227 | -0.001607 |
| 30d same window | -0.0063 | -0.000011 |

The same-window coefficients are practically zero. The old bootstrap confidence
interval is also not decision-grade because it resamples thousands of heavily
overlapping daily windows as if they were independent and sets a fixed random
seed. No new bootstrap was needed to establish the temporal mismatch.

Correct interpretation:

- past intake variability may mark a later phase of the subject's weight cycle;
- AF does not identify a protective physiological effect of variable eating;
- the data do not establish either that dietary variability causes “metabolic
  damage” or that variability prevents it;
- stable-menu advice should rest on the direct ordered menu-repetition results,
  not on AF.

## AB: the exact next-three-day result is null, but the old number is wrong

The current legacy script prints a run-minus-walk difference of `-62 kcal/day`
for `next_3d_cal`, but that variable averages the prior day, exercise day, and
next day. Recomputing the exact mean of days +1, +2, and +3 for the same 90
step-matched pairs gives `+15 kcal/day`.

The conclusion remains a null: neither number is a credible appetite effect,
and the current same-day and next-day differences are also small (`-25` and
`+41 kcal`). The audit changes the number and the permitted interpretation, not
the practical exercise recommendation.

## Next-record is not next-day

Row shifts are only calendar-day shifts on a complete daily index:

| Source table | Unique dates | Consecutive next-row transitions | Maximum gap |
|---|---:|---:|---:|
| Sleep | 2,057 | 91.4% | 80 days |
| Possible-sunlight | 2,012 | 91.5% | 117 days |
| Steps | 4,592 | 95.1% | 14 days |

Consequences:

- AE/BS sleep results required the exact-calendar repair already in report 09.
- AK's next-day sunlight outcomes are next-available-record outcomes. Its
  same-day and clustered RMR sections are separate, though still exploratory.
- CA's 208 stress days are even sparser. Its promoted same-day stress/binge
  association is unaffected; its exported `cal_next` and `binge_next` are not
  next-calendar-day variables.
- N's gravitostat coefficient is technically mislabeled after filtering to
  available step rows, but the current exact-date repair barely changes it:
  partial `r=+0.0552` becomes `+0.0574`. The qualitative null survives.

## Centered transforms: retrospective, not automatically wrong

Several scripts deliberately use future data to describe completed history.
That is permissible if the result is not called a prediction or a
contemporaneous biological state.

| Script | Centered transform | Permitted interpretation |
|---|---|---|
| `D_food_noise_variance.py` | 180-day centered modeled-FM “set point” | Retrospective trend normalization only |
| `DA5_intake_pressure_foods.py` | 30-day centered intake residual | Associative food ranking, not food pressure |
| `DA6_food_macros.py` | Same centered residual in trigger/protective labels | Exploratory retrospective clustering |
| `CR_jump_events.py` | FM slope from 15 prior and 15 future records | Retrospective event segmentation |
| `CS_cycle_peaks_troughs.py` | 14-day centered peak/trough detector | Retrospective cycle description |
| `CT_cycle_pairing.py` | Same centered extrema | Post-hoc cycle pairing |
| `A_energy_balance_quality.py` | 90-day centered phase label | Retrospective phase annotation |

The D “set point” includes up to about 90 future days of modeled fat mass.
Therefore its distance-to-intake gradient is not a forward appetite-pressure
test. DA5/DA6 similarly cannot establish that a named food prospectively
triggers or protects against later intake. Report 03's time-ordered direct
analysis is the appropriate evidence boundary for food-context claims.

## Post-assimilation scoring boundary

This pass intentionally does not duplicate the known P3/P4 audit:

- P3 is intake-derived and therefore not an intake-free confirmation.
- Legacy P4 exports retrospective smooths under filtered names and uses
  outcome channels during selection.

The opt-in P4b model is separately tested for predicted-versus-filtered
semantics and scores observations before the relevant composition measurement
is assimilated. Report 11 documents its remaining calibration limitations.
No additional promoted post-assimilation scoring path was found outside those
already documented model families.

## What changes for the weight-loss plan

Very little operationally:

1. Do not deliberately add calorie variability in hopes that it protects
   metabolism or reduces fat gain.
2. Stable, familiar menu defaults remain reasonable because their lower
   extreme-intake association replicated in ordered raw-data splits—not because
   AF proved that consistency is metabolically harmful or helpful.
3. Running versus walking can be chosen for fitness, preference, injury risk,
   and adherence; this repository does not show a meaningful three-day appetite
   difference.
4. Keep using exact calendar joins and explicit day +1 through day +*h*
   constructions for all future analyses.
