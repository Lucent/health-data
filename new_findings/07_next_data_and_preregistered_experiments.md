# Next data and preregistered experiments

Date: 2026-07-30
Status: prospective plan; medication remains clinician-directed

## First: close the current data gap

The repository cannot evaluate its most recent intervention period:

- Intake ends 2026-02-28.
- Observed weight ends 2026-03-22.
- q3d click entries begin 2026-04-01.
- Medication entries end 2026-04-10.
- Strength dates continue through 2026-07-07.

The first task is not another model. It is exporting MFP from March onward, refreshing scale history, and updating the actual medication log. If the q3d schedule was continued, changed, or stopped, record what happened and whether it was clinician-directed. Do not reconstruct it from memory when an app/pharmacy record exists.

No conclusion about q3d efficacy, mean intake, weight, adverse effects, or “tachyphylaxis” is possible until this boundary is repaired.

## Minimum viable prospective data

Templates are in `prospective/`.

### Daily, about one minute beyond existing logging

- morning weight;
- morning/evening hunger, 0–10;
- nausea, constipation, diarrhea, 0–10;
- travel and order-in flags;
- menu block;
- injection occurrence/time and prescribed dose;
- illness/unusual event.

MFP remains canonical for food, calories, protein, and fiber.

### Each strength session

- exercise;
- sets and reps;
- load;
- effort/RPE;
- duration/distance where relevant.

### Weekly

- waist under one fixed tape protocol;
- blood pressure/resting heart rate if already available and appropriate;
- a consistent strength/function measure.

Blank data is acceptable. Invented or backfilled ratings are not.

## Experiment 1: calorie-neutral protein and menu stability

Goal: determine whether the current diet can reach adequate protein and reduce extreme-intake risk without increasing calories or relying on a single-food diet.

### Design

Use a deterministic ABAB sequence if medication, travel, and health are stable:

- A1: 14 days usual menu
- B1: 14 days stable protein-forward menu
- A2: 14 days usual menu
- B2: 14 days stable protein-forward menu

Do not start a block during travel, illness, dose change, or major GI symptoms. Delay the entire sequence rather than silently changing inclusion rules.

### B-menu rule

- Keep planned calories near the current maintenance level.
- Replace 250–350 kcal/day of recurring low-protein foods with familiar protein-dense foods.
- Use a small rotating menu, not potatoes alone.
- Choose the protein goal with a clinician/dietitian; 80–100 g/day is a reasonable discussion range from the current 55 g baseline.

### Outcomes fixed in advance

Primary:

1. Share of days at or above the agreed protein target.
2. Share of days at or above 3,000 kcal.
3. Completion/adherence to the assigned menu.

Secondary:

- mean calories;
- hunger ratings;
- GI symptoms;
- strength-session performance;
- seven-day mean weight;
- fiber.

Success is not “the B block lost more weight in two weeks.” Success is that it is tolerable, nutritionally better, and does not provoke more extreme-intake days. Weight is secondary at this short horizon.

## Experiment 2: travel playbook

Goal: reduce the historically observed +600 kcal/day travel excursion.

For each future trip, predeclare:

- default breakfast;
- first-day grocery/delivery plan;
- one or two default protein-forward meals;
- which meal(s) are intentionally novel;
- the historical comparison window: 14 days immediately before the trip.

Primary outcome: during-trip mean intake minus the trip's own pre-14-day mean.

Secondary:

- ≥3,000-kcal rate;
- protein target attainment;
- hunger/GI symptoms;
- post-trip 14-day intake;
- adherence to the plan.

Do not call one trip proof. After four planned trips, compare the distribution with the 41 historical trips and report every trip, including failures.

## Experiment 3: progressive resistance training

Goal: preserve or improve function while treatment continues.

This is a training program, not a body-composition model. Choose a clinician-cleared program and record actual load/reps/effort.

Primary outcomes:

- strength progression in a small set of repeatable exercises;
- session completion.

Secondary:

- grip/chair stand;
- same-device composition under standardized conditions;
- protein adherence;
- fatigue/recovery.

Do not infer muscle gain from BIA alone, especially while using creatine.

## Experiment 4: medication observation

Do not manipulate dose or frequency for identification.

If the prescriber changes the plan, predefine an observation window and record:

- actual prescribed schedule;
- injection date/time;
- hunger and GI symptoms;
- intake;
- weight;
- travel/illness;
- other medication changes.

The existing Tuesday schedule makes cycle day and weekday nearly identical. A future clinician-directed schedule that naturally breaks that alignment may be more informative, but research convenience is not a reason to alter treatment.

## Secondary experiments only after the priorities above

- **Sleep:** first clean primary sleeps and focus on genuinely short-sleep deviations. Current median duration is already adequate and historical records contain nap-like artifacts.
- **Fatty acids:** repeat the same lab after at least 16 weeks on a stable verified EPA/DHA dose; this tests uptake, not weight loss.
- **RMR/thermogenesis:** only with standardized indirect calorimetry under thermoneutral conditions. Do not use cold exposure or peripheral temperature as a substitute.
- **Food processing:** compare repeatable meal templates prospectively; do not infer seed-oil causality from incomplete historical PUFA enrichment.

## Analysis rules

1. Write block dates and primary outcomes before looking at results.
2. Keep raw observations separate from derived features.
3. Use calendar-ordered comparisons.
4. Report all blocks/trips, not just successful ones.
5. Do not search many thresholds and report the best.
6. Treat medication, illness, travel, and source/device changes as explicit regime changes.
7. Stop or seek clinical help for concerning symptoms rather than completing an experiment.
