# Tirzepatide: observed outcomes and support priorities

Date: 2026-07-30
Status: raw within-person observations plus current external guidance; not a dose recommendation

## Bottom line

Tirzepatide is the strongest successful intervention in the repository:

- Raw intake fell from 2,454 kcal/day in the pre-drug year to 2,063 in the first 90 treatment days and 2,004 in the last 90 logged days.
- Days at or above 3,000 kcal fell from 11.0% pre-drug to 2.2% in the first 90 days and 1.1% in the last 90.
- Observed scale weight fell from 228.5 lb immediately pre-drug to 202.9 lb at the intake endpoint: -25.6 lb, or -11.2%.
- The final 90 days were close to weight maintenance: 203.6 to 202.9 lb.

The immediate improvement opportunity is **supporting the successful treatment**, not inventing a more aggressive dosing schedule:

- Protein averaged 54.6 g/day in the last 90 days; 60/90 days were below 60 g and only 4/90 reached 80 g.
- Fiber averaged 14.1 g/day, down from 20.7 g/day in the pre-drug 90 days.
- Strength training was logged about 1.24 times/week, but the repo records no exercises, sets, reps, load, or effort.
- Same-device bioimpedance is reassuring about lean mass so far, but it is not precise enough to prove that no muscle was lost.

## Direct intake and weight audit

The script reads raw daily intake, the injection log, and observed scale weights. It does not use modeled blood levels, effective levels, P2/P3/P4, set point, or a calorie multiplier.

Command:

```bash
uv run --with numpy --with pandas python analysis/ED_tirzepatide_direct_audit.py
```

Artifacts: `analysis/ED_tirzepatide_window_summary.csv`, `analysis/ED_tirzepatide_phase_summary.csv`, `analysis/ED_tirzepatide_cycle_summary.csv`

| Window | Mean kcal | Median kcal | ≥3,000 kcal | Protein | Fiber | Observed weight change |
|---|---:|---:|---:|---:|---:|---:|
| Pre-drug 365 days | 2,454 | 2,362 | 11.0% | 64.0 g | 22.9 g | +5.3 lb |
| Pre-drug 90 days | 2,435 | 2,350 | 7.8% | 70.2 g | 20.7 g | +0.1 lb |
| First 90 on-drug days | 2,063 | 2,021 | 2.2% | 58.7 g | 16.8 g | -10.1 lb |
| Last 90 logged on-drug days | 2,004 | 1,974 | 1.1% | 54.6 g | 14.1 g | -0.7 lb |

This is a strong within-person intervention signal and is directionally consistent with randomized tirzepatide trials. It still cannot assign every change to the drug: dose escalation, calendar time, travel, food choices, metformin entries, and behavior change overlap. The magnitude and abruptness make “no useful treatment effect” implausible.

### Dose phases are not a dose-response experiment

Mean intake by most recently injected dose was:

| Assigned phase | Days | Mean kcal | ≥3,000 kcal | Observed weight change |
|---|---:|---:|---:|---:|
| 2.5 mg | 28 | 2,270 | 7.1% | -3.1 lb |
| 5 mg | 140 | 2,003 | 0.7% | -8.8 lb |
| 7.5 mg | 84 | 2,092 | 1.2% | -3.7 lb |
| 10 mg | 84 | 1,984 | 1.2% | -1.1 lb |
| 12.5 mg | 194 | 1,920 | 1.0% | -9.8 lb |

Dose only rose with calendar time, while body weight, travel, food environment, and treatment duration changed. The non-monotonic intake means are a warning against estimating an exact causal dose slope from these phases.

## The weekly sawtooth exists, but its cause is not identified

Raw mean intake rises from 1,595 kcal on cycle day 1 to 2,219 on cycle day 6. That is a real descriptive pattern. However:

- 79 of 80 weekly injections were Tuesday; one was Thursday.
- Cycle day maps to weekday with 99.1% purity.
- Wednesday is therefore almost always “cycle day 1,” Monday almost always “cycle day 6,” and so on.
- Week fixed effects do not solve this because cycle day still equals weekday within almost every week.
- Dose, body weight, cumulative treatment time, and seasonal food patterns also change over the series.

The data can say that intake is lower early and higher late in the usual Tuesday-to-Tuesday week. It cannot identify `-74.5 kcal per modeled blood-level unit`, prove receptor tachyphylaxis, or predict that redistributing the same monthly dose will preserve the mean while removing only variance.

The q3d click entries start on 2026-04-01. Canonical intake ends 2026-02-28 and observed weight ends 2026-03-22. There is therefore **no outcome data in the repo for the q3d experiment**. Any claim that it worked, failed, or preserved mean intake is currently unevaluable.

The current [FDA Zepbound label](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/217806s042lbl.pdf) specifies once-weekly dosing. The [KwikPen patient leaflet](https://www.medicines.org.uk/emc/product/15483/pil) describes fixed doses and says not to use leftover medicine. Schedule and pen-use decisions belong with the prescriber and pharmacist.

## Same-device composition and training support

Command:

```bash
uv run --with numpy --with pandas python analysis/EE_tirzepatide_composition_support.py
```

Artifacts: `analysis/EE_tirzepatide_composition_summary.csv`, `analysis/EE_tirzepatide_support_periods.csv`

The first H30 scan was 33 days after treatment began, so it misses the first month. Across 20 H30 scans from 2024-10-20 to 2026-02-28:

| H30 measure | First | Last | Change |
|---|---:|---:|---:|
| Weight | 224.7 lb | 203.0 lb | -21.7 lb |
| Fat mass | 83.6 lb | 60.0 lb | -23.6 lb |
| Lean mass | 141.1 lb | 143.0 lb | +1.9 lb |
| Skeletal muscle estimate | 80.2 lb | 82.0 lb | +1.8 lb |
| Visceral-fat level | 18 | 11 | -7 |
| Waist-to-hip ratio | 1.07 | 0.98 | -0.09 |

The declines in the same device's visceral-fat level and waist-to-hip ratio are directionally favorable, though they remain device-derived rather than tape-measured clinical endpoints.

The lean result is reassuring: the same device does not show progressive lean loss over the measured interval. It is not proof of a 109% fat-loss fraction. Bioimpedance lean mass includes water and is affected by hydration, sodium, glycogen, timing, and device algorithms; its 137.6–144.2 lb range is much larger than the endpoint difference.

Creatine 5 g/day restarted on 2025-05-09, midway through this interval. The last H30 scan before restart was 139.4 lb lean mass on 2025-01-31; the first afterward was 141.6 lb on 2025-05-15, six days after restart. With no intervening scan, this does not quantify a creatine effect, but it is an important confound because creatine can increase body water and BIA/DXA lean-mass estimates without equivalent new muscle protein: [BIA study](https://pubmed.ncbi.nlm.nih.gov/37675500/), [water/manipulation study](https://pubmed.ncbi.nlm.nih.gov/28410328/). This is a measurement caveat, not a reason to stop a useful supplement.

The result should be phrased as “no obvious measured same-device lean decline,” not “all muscle definitely preserved.”

There were 94 logged strength dates during the 530-day intake-covered treatment era: 1.24 sessions/week. Every session is normalized to 30 minutes. The current log cannot tell whether the program provided progressive overload or whether strength/function improved.

## Protein is the clearest support gap

During the full logged treatment era:

- Mean protein: 56.1 g/day
- Days below 60 g: 316/530
- Days at or above 80 g: 47/530
- Mean relative to the final measured lean mass: 0.87 g/kg lean mass/day

In the final 90 days, mean protein was 54.6 g/day, or 0.84 g/kg of final measured lean mass.

The 2025 joint advisory from four obesity/nutrition organizations emphasizes nutrient adequacy, sufficient protein, and structured strength training during GLP-1 therapy. It explicitly says there is no consensus on the best body-size denominator. It discusses proposed targets of 1.2–1.6 g/kg during active weight reduction, 1.5 g/kg lean mass, or an easier absolute target of 80–120 g/day, all requiring individualization. At 143 lb measured lean mass, 1.5 g/kg lean mass is about 97 g/day. This is a **discussion range**, not an automatic prescription. [Full joint advisory](https://pmc.ncbi.nlm.nih.gov/articles/PMC12264624/).

The gap between 55 g and roughly 80–100 g is large enough to address with the clinician or a registered dietitian, especially because this is a long-term vegetarian diet. Kidney function, GI tolerance, food preferences, total calories, and any history of disordered eating matter.

## Practical support plan

### 1. Preserve the successful calorie level; substitute for protein

Do not simply add several hundred calories of protein food. Replace lower-satiety calories with an additional 25–40 g/day of tolerable protein, then reassess. Vegetarian options include Greek yogurt or cottage cheese if dairy is used, eggs if used, tofu, tempeh, seitan, edamame, textured vegetable protein, lentils/beans paired with more concentrated protein, and a tolerable protein drink.

A reasonable clinician/dietitian discussion target is 80–100 g/day because it is inside the advisory's pragmatic range and near the lean-mass calculation. The right target may be lower or higher.

### 2. Record and progress resistance training

Aim for a clinician-cleared, structured program rather than merely adding session dates. Log exercise, sets, reps, load, and effort. Two or more whole-body sessions per week is a reasonable discussion point, but the program should fit recovery, joint health, experience, and preferences. Track strength or function as an outcome; body-composition numbers alone are insufficient.

### 3. Restore nutrient density and fiber gradually

Fiber fell to 14 g/day in the last 90 days. Increase legumes, vegetables, fruit, whole grains, nuts, and seeds gradually and with adequate fluid, adjusting for nausea, constipation, diarrhea, or early satiety. The goal is nutrient adequacy and tolerability, not forcing a generic gram target through GI symptoms.

Because long-term energy intake is lower and the diet is vegetarian, a clinician/dietitian can review B12, iron/ferritin, vitamin D, calcium, and other intake/labs based on history rather than adding supplements blindly.

### 4. Monitor what BIA cannot

Useful longitudinal outcomes are:

- gym strength and training volume;
- grip strength or a simple chair-stand test;
- waist circumference;
- same-device, similarly hydrated composition readings;
- clinically indicated labs and symptoms;
- function, fatigue, falls, and recovery.

A DXA is optional, not automatically necessary. If there is declining strength, unexplained fatigue, a large further loss, or clinical concern about bone/muscle, a clinician can decide whether DXA or other evaluation is warranted.

### 5. Choose the next weight goal clinically

The last 90 days look like maintenance, not failure. Whether to lose more should depend on waist, blood pressure, glycemia/lipids, sleep apnea, mobility, side effects, preferences, and clinician assessment—not the repo's estimated fat-mass set point.

## What is known versus open

| Question | Current answer |
|---|---|
| Did intake and extreme-intake frequency fall? | Yes, strongly and directly observed |
| Did weight fall? | Yes, about 11% from immediate pre-drug weight |
| Has the effect completely disappeared? | No; late intake remains about 450 kcal/day below the pre-drug year |
| Is late weight stability evidence of a latched set point? | No |
| Is the exact weekly PK appetite coefficient identified? | No; cycle day is weekday |
| Did same-device H30 lean mass decline after day 33? | No obvious measured decline; creatine and hydration confound the estimate |
| Does that prove no muscle loss? | No |
| Is protein/nutrient support worth improving? | Yes; it is the clearest low-risk support opportunity |
| Did q3d work? | Unknown; outcome data ends before q3d begins |
