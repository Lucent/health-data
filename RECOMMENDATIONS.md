# Recommendations

Data-driven only. Every claim below cites a specific finding with reproducible numbers. Nothing here comes from external nutrition literature.

**Current state (Feb 2026):** FM 60 lbs, set point 63 lbs (FM is 3 lbs *below* the set point — the body wants to be heavier). Tirzepatide 12.5mg at 50% effectiveness after 32 weeks (tachyphylaxis HL=32 weeks). Effective drug level ~5-10 through the weekly sawtooth, which offsets 12-25 lbs of set point deficit (AG: each unit of effective level = 2.5 lbs equivalent). Without the drug, 3 lbs below SP would produce ~6% daily binge risk (AG sigmoid). The drug is what's keeping binge risk at zero.

## 1. Walk sessions — the only validated expenditure lever

+14 cal RMR per walk session (30d trailing count), validated by 23 calorimetry measurements (AD: CV RMSE 116, R²=0.49). On tirzepatide specifically, walks are the primary expenditure lever because the drug has suppressed the endogenous defense (AL: partial r=0.19 vs 0.04 pre-tirz). At 15 walks/month: ~39 cal/day recovered of the ~125 cal/day the drug suppresses (AJ). At 20 walks/month: ~52 cal/day.

The count of discrete outings matters more than total minutes (AD: sessions CV RMSE 116 beats minutes at 135) or total steps (179). Running adds nothing over walking at matched steps (AB: null). A 20-minute walk that Samsung logs as an exercise session counts the same as a 60-minute walk. Not a season confound (AK: walks r=0.66 controlling for sunlight; sunlight r=0.36 controlling for walks).

**Priority: the single most important behavior.**

Command: `python analysis/AD_tdee_formula_sweep.py`, `python analysis/AL_walk_rescue_expenditure.py`

## 2. Hold FM stable while the set point converges — the race condition

The set point adapts to current FM with a 50-day half-life (AG: r=-0.62). FM=60 is 3 lbs *below* SP=63 — the gap that would produce ~6% binge risk without the drug. Every month FM stays at 60, the SP drops ~1 lb toward it. After 3 months (May 2026), SP will be ~61 — a 1 lb gap. After 5 months (July 2026), SP ~60.5 — near convergence. In 2018-2021, FM stabilized at 60-71 and binge rate dropped from 10-12% to 3-6% as the SP caught up (AG trajectory table).

**The drug is buying time for this convergence.** The race is between tachyphylaxis (drug losing ~5% remaining effectiveness per month at 32-week HL) and set point adaptation (closing ~1 lb/month). Current drug coverage: effective level 5-10 through the weekly sawtooth = 12-25 lbs of equivalent offset, covering a 3 lb gap with large margin. By August 2026, even at 40% effectiveness (another 16 weeks of decay), the trough still provides ~8 lbs of offset on what will be a <1 lb gap. The margin is adequate but not unlimited.

**The vulnerable scenario is discontinuation before convergence.** AG's on-drug SP half-life is 165 days (vs 50 off-drug), meaning the SP may be adapting slower than estimated — the true SP could still be above 63. If the drug stops while SP > FM, the full binge-rate gradient activates immediately. AH#4 confirms: restriction runs ending below SP rebound (r=-0.48).

**Priority: don't chase further loss and don't discontinue early. Stability is the investment. Time closes the gap.**

Command: `python analysis/AG_binge_set_point.py`, `python analysis/AH_set_point_properties.py`

## 3. Protect protein during any restriction

Low-protein cuts (<58g/day) have the worst TDEE/RMR recovery (J: -0.008 penalty, -24 cal/day). The failure mechanism is distinct from high-step cuts: normal TDEE during the run, poor recovery afterward (M subsection of J). Protein leverage within a day is real (N: r=-0.34, 590 cal range across protein % bins) but has zero next-day carryover (r=-0.04). Tirzepatide weakens the within-day effect (r -0.19 to -0.09).

**Priority: maintain >58g/day protein, especially during calorie restriction. Don't expect protein to control multi-day appetite.**

Command: `python analysis/J_restriction_archetypes.py`, `python analysis/N_dietary_predictors.py`

## 4. Avoid combining restriction with high step counts

High-step calorie restriction is the worst archetype (J: -0.011 TDEE/RMR, -29 cal/day penalty). It suppresses TDEE during the run AND keeps falling afterward (M subsection of J). This conflicts with recommendation #1: walk sessions raise RMR, but walking *during severe restriction* produces the worst metabolic penalty. Falling-phase high-step cuts specifically: post-pre -0.060 (R subsection of J).

**Priority: walk during maintenance or moderate deficit, not during hard cuts. Separate the levers.**

Command: `python analysis/J_restriction_archetypes.py`, `python analysis/R_metabolic_failure_predictors.py`

## 5. Manage the injection-day sawtooth

Intake swings ~500 cal/day within each weekly cycle (F): injection day 1650 cal, trough (day 4-5) 2220 cal. At 50% tachyphylaxis the swings are smaller but still present. The drug's appetite suppression on the trough day is roughly half the peak-day level. With FM 3 lbs below SP, the trough is where the set point's binge pressure is least opposed — the drug breaks binge-to-binge escalation (G: 31.6% → 0%), but this protection is weakest at the trough.

**Priority: awareness, not action. Day 4-5 post-injection is the vulnerable window.**

Command: `python analysis/F_tirzepatide_pk.py`, `python analysis/G_tirzepatide_dynamics.py`

## What the data says doesn't matter

| Intervention | Finding | Result |
|---|---|---|
| Meal timing / front-loading | N | Absolute front-loading r=+0.48 (wrong direction) |
| Fiber | N | Partial r=-0.094 |
| Sleep optimization | AE, AK | r=-0.01 for everything; sunlight also weak |
| Running instead of walking | AB | Null at matched steps and era |
| Gravitostat / weight-bearing | N | r=+0.050 (wrong direction) |
| Willpower / ego depletion | Removed (Q) | AUC 0.696 vs yesterday's calories 0.704 |
| Yo-yo / variance harm | AF | Wrong sign — variance mildly protective |

## Summary

Walk frequently, don't cut hard, keep protein up, and let time close the 3 lb gap between FM and the set point. The drug's most valuable function now is suppressing the ~6% binge risk from being below the set point while the 50-day EMA converges to FM=60. The margin is adequate through mid-2026 even with ongoing tachyphylaxis — the gap is closing faster than the drug is fading. Walking is the only behavioral lever that directly raises measured RMR, and it matters more on the drug than off it because the endogenous expenditure defense has been pharmacologically cleared (AL: partial r 0.04 → 0.19).
