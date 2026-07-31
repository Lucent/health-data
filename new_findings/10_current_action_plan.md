# Current weight-loss and health action plan

Date: 2026-07-30
Status: evidence-ranked working plan, not individualized medical care

## Bottom line

The repository does not support a secret metabolic defect, a precise fat-mass set point, a special genotype diet, or an improvised tirzepatide schedule. It does support a simpler and more useful story:

1. **The current prescribed treatment has worked.** Logged intake fell by roughly 400–450 kcal/day, extreme-intake days nearly disappeared, and observed weight fell 25.6 lb (11.2%) through the last intake-covered weight.
2. **The late-treatment plateau is not failure.** The final 90 logged days averaged about 2,004 kcal/day and weight was nearly stable. That is a plausible current maintenance level under the current environment, not proof that a biological set point has latched.
3. **The biggest support gaps are protein, fiber, resistance-training quality, and ordinary clinical risk data.** These matter more than cold exposure, supplement timing, or genetic anecdotes.
4. **Travel and food-environment disruption are the strongest behavioral hazards.** Historical trips averaged about 600 kcal/day above their own preceding two weeks. A stable familiar menu marked fewer extreme-intake days.
5. **If further loss is desirable, first improve diet composition at the current calorie level.** Then consider a modest additional calorie reduction with the treating clinician, while protecting protein, strength, and tolerability.
6. **High-intake states can persist, but they are interruptible.** Before treatment, an extreme day often continued into the next day, including at home. The low-risk response is an ordinary default meal at the next opportunity—not fasting, punitive exercise, or a medication change.

## The next four weeks

### 1. Repair the present-tense record

Before drawing any conclusion about the latest medication pattern:

- export MyFitnessPal intake from 2026-03-01 onward;
- refresh observed scale weight after 2026-03-22;
- reconcile the medication log after 2026-04-10 with pharmacy/app records;
- record the actual prescribed dose and schedule, not modeled blood levels or reconstructed “clicks.”
- verify whether “500 mg magnesium” means 500 mg elemental magnesium, and reconcile the metformin indication/current use, creatine, and all other supplements with the clinician or pharmacist.

The q3d entries begin after canonical intake and weight end. The repo therefore cannot yet say whether that period worked, failed, or changed adverse effects. Medication frequency and pen use belong with the prescriber and pharmacist.

### 2. Hold calories roughly steady while repairing nutrition

For two to four weeks, keep intake near the established late-treatment level and make calorie-neutral substitutions:

- replace about 250–350 kcal/day of recurring low-protein foods with familiar protein-dense foods;
- discuss a practical protein target around 80–100 g/day with a clinician or dietitian, starting from the observed 55 g/day;
- restore fiber gradually from the observed 14 g/day, adjusting for early satiety, nausea, constipation, or diarrhea;
- preferentially replace some butter, cream, cheese, pastry, and chocolate calories with protein-rich foods and unsaturated-fat sources rather than adding calories.

The last 90 days contained about 914 kcal/day from recurring foods below 3 g protein per 100 kcal. A 300-kcal substitution from roughly 1.5 to 10 g protein per 100 kcal would add about 25 g protein without changing total calories.

Familiar options already present in the log include seitan, tofu, plant-based chick'n products, and vegetarian chili. Dairy protein can be useful if desired and tolerated. The point is a repeatable menu the person likes, not nutritional purity.

Start with breakfast and snacks. They supply about 31% of final-90 calories but
have the lowest protein and fiber density. Beverage-like calories average only
about 72 kcal/day, so liquid calories are a secondary, context-specific lever.
Do not rank the historical diet by kcal/g: exact food mass covers only 10.8% of
final-90 solid calories.

### 3. Make resistance training measurable

The treatment-era log records only about 1.24 strength sessions per week and assigns every session 30 minutes. It cannot establish progressive overload or preserved function.

Use a clinician-cleared whole-body program and record:

- exercise;
- sets and repetitions;
- load;
- effort or reps in reserve;
- any pain or limiting symptom.

Two or more sessions per week is a reasonable discussion point, but progression and adherence matter more than the label “strength workout.” Follow a few repeatable performance measures. Bioimpedance is secondary.

### 4. Obtain an actual cardiometabolic baseline

The repository does not contain an auditable current blood pressure series or conventional lab table. The prose-only LDL and HbA1c statements have no date, source, or complete panel attached.

Discuss ordinary clinical measurements rather than inferring them from food or a Kalman model:

- seated blood pressure and resting pulse;
- fasting lipid panel, with ApoB when clinically appropriate;
- HbA1c or fasting glucose;
- comprehensive metabolic panel, including kidney and liver markers;
- CBC;
- TSH when clinically indicated, especially given severe cold intolerance;
- vegetarian/low-intake adequacy labs guided by history, commonly including B12, ferritin/iron studies, and vitamin D.

This is a missing-baseline list, not a claim that any result is abnormal. The exact set should be chosen by the treating clinician.

See the source-by-source audit and medication/supplement reconciliation in
[`08_cardiometabolic_and_clinical_monitoring_audit.md`](08_cardiometabolic_and_clinical_monitoring_audit.md).

## The next eight to twelve weeks

### If further weight loss is clinically worthwhile

Once the protein-forward menu is easy, strength is being tracked, and medication is stable:

- remove a modest additional 150–250 kcal/day from remaining low-protein discretionary sources;
- keep protein and resistance training stable;
- judge the result from four-week weight trends, waist, strength/function, hunger, and side effects;
- do not respond to a noisy week by repeatedly cutting calories or changing medication.

The repo does not identify a single correct target weight. The decision should use waist, blood pressure, metabolic risk, sleep apnea if relevant, function, preferences, and treatment tolerability. At 5 ft 11.75 in and 202.9 lb, BMI is about 27.7, but BMI alone does not determine the goal.

### Run the menu experiment

Use the preregistered ABAB plan in [`07_next_data_and_preregistered_experiments.md`](07_next_data_and_preregistered_experiments.md):

- 14 days usual menu;
- 14 days stable protein-forward menu;
- repeat once;
- delay blocks during travel, illness, or a clinician-directed dose change.

Primary outcomes are adherence, protein attainment, and the frequency of days at or above 3,000 kcal. Short-block weight loss is secondary.

This translation is supported by the named-epoch audit: all four potato
attempts cut intake by 524–748 kcal/day versus their immediate pre-windows, but
rapid scale changes substantially reversed and protein was very low in several
attempts. Weekend fasting likewise produced a sharp deficit followed by a
2,364-kcal/day post-window and a 7.3 lb observed scale rebound. Borrow the
structure; do not repeat the monodiet or extreme fast.

### Treat travel as an intervention

Before each trip, decide:

- the default breakfast;
- a first-day grocery/delivery option;
- one or two protein-forward default meals;
- which meal will be intentionally novel.

The realistic first target is eliminating part of the historical +600 kcal/day travel excursion. Do not plan to “burn it off” afterward; post-trip intake historically did not compensate.

### Break a high-intake run without punishment

Before treatment, 22.2% of home extreme-intake days were followed by another
extreme day, versus 5.4% after other days. Recent calorie momentum also ranked
next-day risk in both held-out pre-treatment periods, although its probabilities
were not stable enough for an automated score.

If intake is already running above the recent baseline:

- make the next eating opportunity the familiar protein-forward default;
- remove leftover delivery or party food from immediate availability when practical;
- continue the normal medication and exercise plan;
- do not turn the reset into a compensatory fast.

### Monitor sustained regimes, not isolated scale noise

Use the 28-day intake mean, high-intake frequency, and a four-week local-median
weight trend as a simple contemporaneous dashboard. No one-month feature score
reliably predicted the next quarter in the internal locked test, so do not
pretend this is a forecast.

Mark weight data stale at 14 days and use 30 days without an observed weight as
a neutral “resume measurement” alert. All five historical gaps of at least 90
days ended higher, but that is informative missingness—not proof that weighing
itself prevents gain.

## Hypotheses ranked by present usefulness

| Hypothesis | Subject-specific status | Action |
|---|---|---|
| Appetite/reward biology responds to anti-obesity medication | Strongly consistent with observed response and trials | Preserve effective supervised treatment; monitor nutrition and side effects |
| Food environment, novelty, and travel disrupt intake control | Strongest replicated behavioral signal in the repo | Stable defaults, planned novelty, travel playbook |
| High-intake momentum can continue into the next day | Replicated descriptive pre-treatment signal; risk-score calibration shifts | Reinstate the ordinary default meal early; do not fast or punish |
| Protein/fiber dilution matters during low-appetite treatment | Direct current gap plus external support | Calorie-neutral protein and nutrient-density substitution |
| Highly palatable/ultra-processed food can raise intake | Externally plausible; the measurable fat+sugar cluster marks higher same-day intake but has no replicated next-day extreme-risk signal | Change defaults and portions; a blanket ban is not supported |
| Saturated fat is relevant to conventional lipid risk | Plausible current dietary concern; current labs missing | Replace some sources and measure an actual lipid/ApoB baseline |
| Sleep loss impairs appetite and metabolic health | Externally supported; reliable 2022+ data average 8.03 h with only six under-5 h nights, and add negligible holdout prediction | Maintain the routine; screen clinically for apnea symptoms rather than forcing an earlier bedtime |
| Energy-expenditure adaptation occurs during weight loss | Real generally, not quantified by the old models here | Use observed trends; do not chase cold exposure or inferred TDEE |
| Seed oils/linoleic acid are the main cause | Not demonstrated; PUFA coverage is incomplete and human inflammation evidence is unsupportive | No seed-oil-specific prescription |
| Common obesity SNPs dictate treatment | Not supported | Ignore for day-to-day treatment decisions |
| A precise latched set point explains the plateau | Not established; old pipeline leaks future/outcome information | Do not use it for calorie or medication decisions |
| Deliberately varying intake protects metabolism | Contradicted as established; the legacy effect disappears when exposure/outcome windows align | Prefer consistency for adherence, not a metabolic-variance trick |
| Cold exposure, taurine, or magnesium timing unlocks loss | Weak or contradicted, and cold is a poor safety fit | Deprioritize |

## What would change the plan

Escalate clinical review rather than the calorie deficit if there is:

- declining strength or function;
- persistent vomiting, severe constipation, dehydration, or inability to meet basic nutrition needs;
- marked fatigue, orthostasis, syncope, or worsening cold intolerance;
- concerning blood pressure, glycemia, kidney/liver markers, or lipids;
- symptoms of gallbladder disease, pancreatitis, or another urgent adverse event.

Conversely, if nutrition, strength, cardiometabolic markers, waist, and quality of life are all good, maintaining the current loss may be a valid success even if the scale is not falling.

## Decision hierarchy

1. **Do not improvise medication.**
2. **Refresh missing current data and get conventional health markers.**
3. **Protect protein, fiber, and progressive strength at roughly current calories.**
4. **Engineer travel and menu defaults, and interrupt high-intake momentum early.**
5. **Only then test a modest calorie reduction if more loss is wanted and clinically appropriate.**
6. **Spend research effort on prospective outcomes and a valid forward model, not additional stories built on the quarantined pipeline.**

The opt-in P4b model now has correct forward/smoothed semantics and explicit
uncertainty, but its intervals are over-wide and its fat-mass prediction does
not beat a simple weight-adjusted scan baseline. See
[`11_P4b_forward_model_audit.md`](11_P4b_forward_model_audit.md). It remains a
research scaffold rather than a decision engine.
