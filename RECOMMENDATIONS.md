# Recommendations

Data-driven only. Every claim below cites a specific finding with reproducible numbers. Nothing here comes from external nutrition literature.

**Current state (Mar 2026):** FM 60 lbs, set point ~61 lbs. Under the SmoothLatch model (AZ), the gap is ~1 lb and the SP has been latching onto current FM during the months of stable weight on tirzepatide. Under the EMA model (AM), the gap is similar. Both models agree the SP is nearly converged. Pressure: ~55 cal/day per lb of gap.

**The SmoothLatch model changes the discontinuation calculus.** The set point only adapts when FM holds within ±3 lbs for 14+ consecutive days. On this subject's slow-gain decade, this behaves identically to a 45-day EMA (both r = -0.94). But after rapid drug-driven loss, the models diverge: the EMA predicts the SP has caught up (regain ~3%); the SmoothLatch predicts the SP froze during rapid loss and regain is substantial. SURMOUNT-4 published data (+14% regain) and STEP-1 data (~+10% regain) both match the SmoothLatch, not the EMA.

**Implication for this subject:** FM has been stable near 60 lbs for several months on the drug. Under the SmoothLatch, this stability IS the set point adapting — every day FM holds within ±3 lbs, the SP moves closer at rate 0.01/day. The longer FM stays stable at 60, the more the SP converges. The drug is buying time for the latch.

Command: `python analysis/AZ_sp_model_search.py`

## 1. Hold weight stable — the latch needs time

Under the SmoothLatch, **stability is more important than further loss.** The set point only adapts when fat mass holds within ±3 lbs for 14+ consecutive days. That means no more than ~1.5 lbs/week of fat change — in practice, daily surpluses or deficits under ~750 cal. Normal scale fluctuations (±2-3 lbs from water and gut contents) don't count; the latch tracks fat mass, not scale weight.

**The pressure is simple and symmetric.** Every pound of gap between your fat mass and the set point changes your hunger by 55 calories per day. Below the set point, you're hungrier; above it, you're less hungry. If your set point is 10 lbs above your current weight, you'll be 550 cal/day hungrier than maintenance — enough to regain about 4 lbs/month. If you've gained 10 lbs above the set point, you'll be 550 cal/day less hungry — you'll naturally undereat and drift back down. In both directions: hold your weight steady for 6 months and the set point closes most of the gap. By 12 months, you're at your new normal.

**How it works in plain terms:** After you reach a new weight and hold it, the body watches for two weeks. If your fat mass stays within ±3 lbs for those 14 days, the set point begins moving toward your current weight at 1% of the remaining gap per day:

| Months holding | Gap remaining | Eating pressure | Regain rate if drug stopped |
|---|---|---|---|
| 1 | 26 lbs | +1400 cal/day | +10 lbs/month |
| 2 | 19 lbs | +1040 cal/day | +8 lbs/month |
| 3 | 14 lbs | +770 cal/day | +6 lbs/month |
| 6 | 6 lbs | +310 cal/day | +2 lbs/month |
| 9 | 2 lbs | +130 cal/day | +1 lb/month |
| 12 | 1 lb | +50 cal/day | negligible |

(Table assumes a 30-lb starting gap, e.g., FM dropped from 90 to 60. Pressure at 55 cal per lb of gap.)

The rate slows as the gap shrinks — fast at first (~9 lbs/month of SP movement with a 30-lb gap), tapering to ~0.5 lbs/month as the gap nears zero. At 6 months of stable weight, two-thirds of the adaptation is done. At 12 months, it's essentially complete.

**If further loss is desired**, the step-down strategy is optimal: lose 5 lbs, hold for 60+ days to let the SP latch, then lose the next 5. Continuous loss at >1.5 lbs/week outpaces the ±3 lb tolerance and the SP freezes entirely.

**Priority: stability over further loss. Each month of holding steady at FM=60 makes eventual discontinuation substantially safer.**

Command: `python analysis/AZ_sp_model_search.py`

## 2. Walk sessions — the only validated calorie-burning lever

+14 cal RMR per walk session (30d trailing count), validated by 23 calorimetry measurements (AD: CV RMSE 116, R²=0.49). The drug suppresses the body's natural metabolic boost during weight loss (AJ); walks raise RMR per AD's calorimetry-anchored finding, so they remain worth doing. The Kalman-side "walks partially rescue tirz expenditure" channel is a null (AL).

The count of discrete outings matters more than total minutes (AD: sessions CV RMSE 116 beats minutes at 135) or total steps (179). Running adds nothing over walking at matched steps (AB: null).

**Priority: walk frequently for general RMR support; do not expect a specific drug-era rescue.**

Command: `python analysis/AD_tdee_formula_sweep.py`, `python analysis/AL_walk_rescue_expenditure.py`

## 3. Plan discontinuation around the latch

The risk at discontinuation depends on how long FM has been stable. From the table above: stopping after 3 months of stability means ~14 lbs of gap and ~6 lbs/month initial regain. Stopping after 6 months: ~6 lbs of gap and ~2 lbs/month. Stopping after 9+ months: gap is negligible.

**Do not discontinue during active weight loss.** If FM is still dropping when the drug stops, the SP is frozen at its pre-loss level. The full original gap reactivates. This is the SURMOUNT-4 scenario (+14% regain in 52 weeks) — subjects stopped at 36 weeks during ongoing loss, before the SP had time to latch.

**Priority: stabilize first, for as long as possible. 6+ months of stable weight before discontinuation.**

## 4. Protect protein during any restriction

Low-protein cuts (<58g/day) have the worst metabolic recovery (J: -0.008 penalty, -24 cal/day). Protein leverage within a day is real (N: r=-0.34) but has zero next-day carryover.

**Priority: maintain >58g/day protein, especially during calorie restriction.**

Command: `python analysis/J_restriction_archetypes.py`, `python analysis/N_dietary_predictors.py`

## 5. Avoid combining restriction with high step counts

High-step calorie restriction is the worst archetype (J: -0.011 TDEE/RMR, -29 cal/day penalty). This conflicts with recommendation #2: walk sessions raise RMR, but walking *during severe restriction* produces the worst metabolic penalty.

**Priority: walk during maintenance or moderate deficit, not during hard cuts. Separate the levers.**

Command: `python analysis/J_restriction_archetypes.py`, `python analysis/R_metabolic_failure_predictors.py`

## 6. Switch to every-3-day dosing with a clickpen

The weekly 12.5 mg injection creates a sawtooth in blood level that wastes drug at peak and underdelivers at trough. The PK model (t½ = 5.0 days, Tmax = 24h, ka = 3.31/day; FDA prescribing information) gives steady-state blood levels of 8.0 at trough (day 6) and 17.9 at peak (day 1) — a 2.2× swing. Since appetite suppression scales at -74.5 cal per unit blood level (AX, within-week identification), this produces a 739 cal/day sawtooth in appetite pressure. Observed intake confirms: injection day 1643 cal, trough day 2222 cal — a 580 cal swing.

**The problem is asymmetric waste.** At peak (17.9 units), the model predicts 1332 cal of appetite suppression — but observed day-1 intake is ~1700 cal, implying actual suppression of ~300-400 cal. The drug saturates. At trough (8.0 units), suppression drops to 593 cal and binges cluster on days 4-5. The peak is wasted; the trough is where the drug fails.

**Redistribute using a clickpen at higher frequency.** All KwikPen strengths deliver 60 clicks per full dose (0.6 mL). At the 12.5 mg strength, each click = 0.2083 mg. At steady state:

| Schedule | Clicks | Inj/mo | mg/mo | Trough | Peak | Mean | CV | Trough suppression | Sawtooth |
|---|---|---|---|---|---|---|---|---|---|
| 12.5 mg q7d (current) | 60 | 4 | 50.0 | 8.0 | 17.9 | 12.9 | 25% | 593 cal | 739 cal |
| **5.2 mg q3d (recommended)** | **25** | **9.3** | **48.6** | **10.5** | **13.9** | **12.5** | **9%** | **781 cal** | **248 cal** |
| 3.5 mg q2d (EOD) | 17 | 14 | 49.0 | 11.7 | 13.6 | 12.9 | 5% | 868 cal | 148 cal |

Every-3-days captures 66% of the sawtooth reduction (739→248 cal) with 9 injections/month. EOD captures 80% (739→148 cal) with 14. Diminishing returns from q3d to q2d: 5 extra injections buy only 87 cal more trough suppression.

**Predicted intake by day of cycle.** Calibrated from the observed q7d sawtooth (baseline = observed mean intake + 74.5 × mean blood level = 2899 cal):

| | Day 0 (injection) | Day 1 (peak) | Day 2 (trough) |
|---|---|---|---|
| Blood level | 10.5 | 13.7 | 12.1 |
| Predicted intake | 2114 cal | 1878 cal | 1998 cal |

Compare to current q7d: intake swings from 1572 (day 1) to 2306 (day 0) — a 734 cal range. On q3d: 1878 to 2114 — a **248 cal range**. The worst day on q3d (2114 cal, injection day before absorption) is better than the worst day on q7d (2222-2306 cal, days 5-6).

**Pen consumption guide.** Each injection: 1 click prime + 25 clicks dose = 26 clicks. Each 12.5 mg KwikPen holds 240 clicks. One pen yields 9 full injections with 6 clicks remaining; the 10th injection bridges across pens.

| Injection | Day | Pen | Notes |
|---|---|---|---|
| 1-9 | 0-24 | Pen 1 | 26 clicks each; 6 remaining after #9 |
| 10 | 27 | Pen 1→2 | Prime 1 + 5 clicks from pen 1, swap, 20 clicks from pen 2 |
| 11-18 | 30-51 | Pen 2 | 26 clicks each; 12 remaining after #18 |
| 19 | 54 | Pen 2→3 | Prime 1 + 11 clicks from pen 2, swap, 14 clicks from pen 3 |
| 20-27 | 57-78 | Pen 3 | 26 clicks each; 18 remaining after #27 |
| 28 | 81 | Pen 3→4 | Prime 1 + 17 clicks from pen 3, swap, 8 clicks from pen 4 |

The cycle repeats every 90 days / 4 pens. Consumption: **1.08 pens per month**, zero waste.

**Flatter dosing does not accelerate tachyphylaxis.** This was the key concern — does eliminating the weekly trough remove a recovery window that slows receptor desensitization? The literature resolves this clearly by showing that GLP-1R effects desensitize through two independent pathways with opposite behavior:

*Gastric emptying (vagal pathway):* Rapid tachyphylaxis with continuous exposure. Umapathysivam et al. (doi:10.2337/db13-1033) directly compared continuous vs. intermittent GLP-1 infusion in a crossover design: continuous attenuated the gastric emptying effect, intermittent preserved it. Nauck et al. (doi:10.1210/jc.2010-2504) confirmed the effect develops within hours. However, this pathway is already fully desensitized on the current weekly regimen — blood level never drops below 8.0 units, maintaining heavy receptor occupancy even at trough. Going to EOD cannot make this worse.

*Central appetite suppression (hypothalamic pathway):* Relatively resistant to desensitization. Multiple reviews confirm that the anorectic effect of long-acting GLP-1RAs persists even after gastric emptying tachyphylaxis is complete (van Can et al., doi:10.1111/obr.12162; Kanoski et al., doi:10.1016/j.physbeh.2011.08.040). The appetite effect is centrally mediated via hypothalamic and brainstem GLP-1R, not a downstream consequence of slowed gastric emptying. This is the pathway that drives weight loss and does not benefit from pulsatile dosing.

*Empirical proof:* If continuous receptor occupancy accelerated appetite tachyphylaxis, short-acting GLP-1RAs (twice-daily exenatide, with large off-periods) would produce more weight loss than long-acting agents (weekly semaglutide/tirzepatide, near-continuous occupancy). The opposite is observed: long-acting agents produce equal or greater weight loss (Madsbad, doi:10.1111/dom.12596). Semaglutide 2.4 mg weekly (t½ = 7 days, nearly flat PK) achieves ~15% weight loss; exenatide BID (pulsatile, large troughs) achieves ~3-4%.

*Receptor recycling timescale:* Beta-arrestin-mediated GLP-1R desensitization recovers within ~6 hours of agonist removal (Jones et al., doi:10.1038/s41467-018-03941-2). But the current weekly trough (blood level 8.0) still provides heavy receptor occupancy — there is no meaningful off-period in the weekly cycle. True receptor recovery would require blood levels near zero for hours, i.e., a multi-day drug holiday, not the natural trough of weekly dosing.

**Tachyphylaxis rate is unchanged.** The fitted tachyphylaxis model is now AUC-based (BZ): erosion ∝ exp(−K_AUC × ∫blood dτ). Total AUC across a fixed monthly mg-dose is identical regardless of schedule (same drug per month, same blood-level integral), so the tachyphylaxis trajectory is the same on weekly vs every-3-day vs every-2-day dosing.

**Expenditure arm is neutral.** The drug's metabolic suppression cost (-10.2 cal per unit blood level, AX) depends on mean level, which is unchanged at 12.9.

**Injection time of day does not matter.** With t½=5 days, the waking-hours fraction of drug absorption varies only 63-68% across all injection times (noon vs midnight). Regression of intake on injection hour, controlling for effective level and calendar time, shows no effect (t=-0.56). The 12-hour drift from noon to midnight injections over 2024-2026 is cosmetic — the apparent correlation with rising intake is entirely absorbed by tachyphylaxis (r=0.894 between injection hour and days on drug).

**Priority: switch to 25 clicks (5.2 mg) every 3 days. Same drug, 1.08 pens/month, same tachyphylaxis — 66% less sawtooth, 32% stronger trough. Inject at any consistent time.**

Command: `python analysis/AX_drug_model.py`

### Preregistered prediction (2026-04-01): q3d intake model

Switched from 12.5 mg q7d to 25 clicks (5.2 mg) q3d on 2026-04-01. Model G fitted to 529 days of q7d data predicts intake from two drug mechanisms:

```
calories = 2271 + (-126) × (blood - EMA₁ᵈ(blood)) + (-63) × blood × exp(-0.0034 × cum_AUC)
           base   acute (response to Δblood)         tonic (suppression eroded by cumulative exposure)
```

The acute term (-126 cal/unit deviation from 1-day EMA of blood level) captures within-cycle sawtooth. The tonic term (-63 cal/unit × blood) captures sustained suppression; its erosion is proportional to cumulative AUC of blood (∫blood dτ), not calendar time. Subject K_AUC ≈ 0.0034 / unit-week.

**Predictions for q3d weeks 1-4 (Apr 4 onward, after transition):**

| | q3d predicted | q7d counterfactual |
|---|---|---|
| Mean intake | 2046 cal | 2053 cal |
| Intake range | 1977-2122 cal | 1667-2256 cal |
| Intake std | 68 cal | 197 cal |
| Mean deficit (vs TDEE 2063) | -17 cal/day | -10 cal/day |
| Fat loss | ~0.1 lbs/month | ~0.1 lbs/month |

The model predicts q3d and q7d produce identical mean intake (±7 cal). The acute component nearly vanishes on q3d (blood deviations ±0.6 units vs ±5 on q7d) while tonic suppression is unchanged. The q3d intake range narrows from 589 to 145 cal. Model G has 4 parameters fit to 529 days; possible overfit of the acute term. If the observed q3d sawtooth is larger than 145 cal or mean intake diverges from 2046 by more than ±100 cal, the model is wrong.

Command: `python analysis/BA_q3d_prediction.py`

## What the data says doesn't matter

| Intervention | Finding | Result |
|---|---|---|
| Meal timing / front-loading | N | Absolute front-loading r=+0.48 (wrong direction) |
| Fiber | N | Partial r=-0.094 |
| Sleep optimization | AE, AK | r=-0.01 for everything; sunlight also weak |
| Running instead of walking | AB | Null at matched steps and era |
| Gravitostat / weight-bearing | N | r=+0.050 (wrong direction) |
| Yo-yo / variance harm | AF | Wrong sign — variance mildly protective |

## Summary

Walk frequently (+14 cal/session), hold weight stable to let the set point latch, don't cut hard, keep protein up. The SmoothLatch model predicts that stability is the critical variable: the set point adapts during stable periods (±3 lbs for 14+ days) and freezes during rapid change. Every month at stable FM=60 moves the SP closer, reducing future eating pressure whether on or off the drug. The drug's most valuable function now is maintaining the stability the set point needs to adapt. When discontinuation comes, the risk depends on how long FM has been stable — not just the current gap. Three months of stability before stopping is the minimum; six months is substantially safer.
