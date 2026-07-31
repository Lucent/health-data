# P4b: leakage-resistant forward model, first audit

Date: 2026-07-30
Status: structurally sound opt-in v1; not calibrated well enough for causal or medical use

## Bottom line

P4b fixes the most serious semantic defects in the legacy P4 pipeline:

- it reads only canonical intake, observed scale weight, and composition scans;
- it does not read P1/P2/P3/P4, RMR, temperature, strength, or downstream results;
- future observations do not alter earlier forward-filtered states or covariance;
- predicted, filtered, and retrospectively smoothed estimates are separate;
- all state uncertainties and covariances are exported;
- scale weight is not assimilated twice on composition days;
- lean mass derived from the same scan is not treated as an independent observation;
- no configuration sweep used the locked period.

This is an important repair, but it does **not** make the latent states established physiology. The v1 uncertainty is too broad in later periods, and its fat-mass prediction does not beat a simple baseline that adds observed total-weight change to the previous scan's fat mass. The expenditure state remains an intake-conditional latent accounting state, not independently measured TDEE.

P4b therefore stays opt-in. It does not overwrite P4 and no downstream analysis has been switched to it.

## Reproducer and artifacts

Run:

```bash
uv run --with numpy --with pandas python analysis/P4b_forward_state_space.py
uv run --with numpy --with pandas python analysis/test_P4b_forward_state_space.py
uv run --with numpy --with pandas python analysis/P4b_posthoc_baselines.py
```

Artifacts:

- `analysis/P4b_state_daily.csv`
- `analysis/P4b_covariance_daily.csv`
- `analysis/P4b_temporal_validation.csv`
- `analysis/P4b_config_selection.csv`
- `analysis/P4b_posthoc_baselines.csv`

The 5,429-row daily output begins at the intake boundary and initializes its state at the first composition scan. The covariance artifact contains the ten unique covariance elements for each of predicted, filtered, and smoothed states on 5,421 initialized days: 16,287 rows.

## Fixed v1 model

The state is:

1. fat mass;
2. fat-free mass;
3. transient buffer mass;
4. expenditure conditional on logged intake.

The primary intake multiplier is fixed at 1.0. The expenditure state is therefore named `tdee_logged_kcal`; it must not be presented as true physiological expenditure without an independent calibration.

Energy partition uses the Forbes relationship without the legacy workout coefficient. The initial expenditure prior is 2,200 ± 600 kcal/day. Device-error and process-error assumptions are fixed and recorded in `P4b_config_selection.csv`; they were not chosen by looking at 2025–2026 performance.

Each daily artifact row records:

- raw observations and device;
- which channels were assimilated;
- days since the prior scale observation;
- pre-update innovations and standardized innovations;
- predicted, filtered, and smoothed states and standard deviations;
- covariance-repair and device-fallback flags.

No covariance repair was required.

## Semantic verification

All 11 deterministic tests pass.

| Requirement | Result |
|---|---|
| Raw-input isolation | Pass |
| Intake multiplier fixed at 1.0 | Pass |
| Morning-state energy indexing | Pass |
| No future leakage in filtered states/covariance | Pass to `1e-12` |
| Deterministic repeatability | Pass |
| Covariance finite, symmetric, positive semidefinite | Pass |
| Filtered and smoothed semantics differ before terminal day | Pass |
| No double counting of total weight | Pass |
| No derived lean-mass assimilation | Pass |
| Complete exported state SDs | Pass |
| Temporal scores use pre-channel-update innovations | Pass |

Repeated serialized runs produced identical outputs.

## First temporal results

The periods were fixed as:

- development: through 2023;
- calibration: 2024;
- first locked test: 2025–2026.

### Point prediction

| Period | Channel | n | MAE | RMSE |
|---|---|---:|---:|---:|
| Development | Total weight | 1,010 | 0.973 lb | 1.449 lb |
| Development | Composition fat mass | 48 | 3.144 lb | 4.152 lb |
| 2024 calibration | Total weight | 256 | 0.632 lb | 0.778 lb |
| 2024 calibration | Composition fat mass | 8 | 1.909 lb | 2.124 lb |
| 2025–2026 locked test | Total weight | 294 | 0.648 lb | 0.825 lb |
| 2025–2026 locked test | Composition fat mass | 13 | 1.833 lb | 2.223 lb |

These are predictions immediately before the relevant observation channel is assimilated. On a day with scale weight and a composition scan, the fat-mass prediction may condition on the same-day scale-weight update but not on the scan fat mass.

### Uncertainty calibration fails conservatively

In the locked period:

- total-weight standardized residual RMS was 0.45, with 99.0% inside the nominal 80% interval;
- composition-fat-mass standardized residual RMS was 0.52, with 100% inside the nominal 80% interval.

A calibrated model would have standardized RMS closer to 1 and coverage closer to the advertised level. Here, the intervals are substantially too wide. This is safer than false precision, but it fails the proposed calibration gate.

### Simple post-hoc baselines

The locked interval had already been opened before these baselines were run, so this is explicitly post-hoc:

| Locked-test channel | P4b MAE | Last-observation MAE | Last FM + observed total-weight change MAE |
|---|---:|---:|---:|
| Total weight | 0.648 lb | 0.716 lb | — |
| Composition fat mass | 1.833 lb | 2.085 lb | 1.785 lb |

P4b modestly improves over carrying the last observed weight forward. For fat mass, it improves over carrying the last scan FM forward but is slightly worse than the simple weight-adjusted scan baseline. With only 13 locked fat-mass scans, the difference is not stable evidence of superiority.

This is exactly why a naive comparator belongs in model validation: a complicated state-space model should have to beat a simple forecast, not merely fit observations plausibly.

## What P4b can and cannot be used for

### Reasonable current uses

- provenance-aware visualization of observed versus predicted/filtered/smoothed states;
- explicitly exploratory sensitivity analyses;
- checking whether a conclusion survives use of only forward-filtered information;
- serving as an engineering scaffold for a future preregistered model.

### Not justified

- claiming an independently measured daily TDEE;
- validating a 9–14% intake multiplier by energy closure;
- estimating a precise biological set point or adaptation rate;
- attributing weight change to metabolic suppression versus logging error;
- selecting a medication schedule;
- presenting retrospective smooths as predictions;
- switching the legacy downstream findings wholesale from P4 to P4b.

## Next model step

The 2025–2026 interval has now been observed and should remain the frozen v1 test in the historical record. Any configuration changed in response to these results must call that interval calibration, not locked validation.

Before loading newly refreshed post-February outcomes:

1. predeclare v2 configuration and acceptance gates;
2. add naive comparators to the primary scorecard;
3. predeclare intake-multiplier sensitivities such as 1.00, 1.05, 1.09, and 1.15 without selecting by closure;
4. calibrate device/process errors using development data only;
5. add sensitivity to initial expenditure and state uncertainty;
6. test informative missingness, because long weigh-in gaps are not random;
7. specify what level of point-error improvement and interval calibration is required.

Then use genuinely new post-February weight/composition data once, as the next locked test. Until that happens, P4b is a much more honest model representation, but still not a decision engine.
