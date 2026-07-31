# Cardiometabolic and clinical monitoring audit

Date: 2026-07-30

This is an observation-first audit, not a diagnosis or a prescription. The
purpose is to identify substantial non-weight health information in the
repository, distinguish current evidence from stale or modeled claims, and
define the smallest useful clinical-data refresh.

## Bottom line

The repository is unusually strong for longitudinal intake, weight, medication
dates, and repeated consumer body-composition measurements. It is unusually
weak for the ordinary clinical measurements needed to interpret those data
safely.

No auditable current conventional lipid panel, ApoB, Lp(a), HbA1c, fasting
glucose, blood pressure, CBC, metabolic panel, creatinine/eGFR, liver enzymes,
TSH/free T4, ferritin, B12, or vitamin D result was found. “Not found” is not
“normal.” The historical prose claims of LDL 130 and HbA1c around 5.0 have no
date, units, source report, or trend table in the repository.

The clearest new diet signal is high enough to act on without invoking a
special obesity theory: in the final 90 logged days, matched foods contributed
**33.5 g/day of known saturated fat, or 15.0% of logged calories**. Saturated-fat
values covered 90.3% of logged calories, so 15.0% is a known-food estimate, not
an upper bound. The 2026 American Heart Association statement says
heart-healthy patterns are unlikely to exceed 10% of energy from saturated fat
and emphasizes replacing saturated with unsaturated fat, not merely subtracting
food: [AHA dietary guidance](https://professional.heart.org/en/science-news/2026-dietary-guidance-to-improve-cardiovascular-health/top-things-to-know).

This does **not** show that saturated fat caused the subject's weight history.
It does make a current lipid/ApoB measurement and calorie-neutral replacement
of some butter, cream, cheese, pastry, chocolate, and queso calories unusually
high-value.

## Evidence inventory and freshness

| Domain | Best repository evidence | Last date | What it can establish | Important limitation |
|---|---|---:|---|---|
| Daily intake | `intake/intake_daily.csv` and item log | 2026-02-28 | Complete diary totals and foods | Absolute intake error remains uncertain |
| Enriched nutrients | `intake/mfp_food_dict.csv` | fetched in 2026 | Saturated-fat and potassium estimates for matched foods | Database enrichment, incomplete coverage, not biomarkers |
| Body composition | `composition/composition.csv`, H30 | 2026-02-28 | Same-device directional change | BIA is hydration-sensitive; no tape waist |
| RMR | `RMR/rmr.csv` | 2025-09-03 | Direct calorimetry readings | Sparse; substantial test-to-test variation |
| Fatty acids | `fatty-acids/fatty_acids.csv` | 2025-09-05 | Manually transcribed fatty-acid values | Source lab PDFs absent; 2019 used a different lab/method |
| Temperature | `temperature/temperature.csv` | 2026-03-23 | 1,315 logged readings | Site, device, and standardized protocol are not recorded |
| Heart-rate proxy | `analysis/BX_tirz_heart_rate_summary.csv` | era summary | Derived overnight median comparison | Raw Samsung archive is external and unavailable here |
| Medication log | `drugs/medicine.csv` | 2026-04-10 | Recorded administrations | It does not prove the current regimen |
| Conventional labs | no auditable source identified | unknown | Nothing current | Prose snippets cannot substitute for reports |
| Blood pressure | no source identified | unknown | Nothing current | Cannot calculate a defensible cardiovascular-risk estimate |

“Current” above means current relative to the repository, not necessarily
current for clinical decisions in July 2026. The intake and composition series
end five months earlier; the latest RMR and fatty-acid tests are about eleven
months old.

## 1. Lipids and saturated fat

`BACKGROUND.md` says “Elevated LDL (130) which lowered in response to 6 months
keto.” There is no date, unit, HDL, triglyceride, total cholesterol, ApoB,
non-HDL cholesterol, Lp(a), treatment context, or source report. The direction
of a historical response cannot be checked.

The 2019 Cleveland HeartLab source apparently included more than fatty acids,
but only fatty-acid values were transcribed and the source PDF is absent.
Omega-3 Index, AA:EPA, and omega-6:omega-3 ratios are **not** substitutes for a
standard lipid panel or ApoB.

### Direct current-window diet result

For 2025-12-01 through 2026-02-28:

| Measure | Result | Evidence class |
|---|---:|---|
| Calories | 2,004 kcal/day | direct diary total |
| Known saturated fat | 33.5 g/day | database-enriched food estimate |
| Known saturated-fat energy | 15.0% of calories | deterministic calculation |
| Saturated-fat calorie coverage | 90.3% | deterministic coverage measure |
| Sodium | 2,598 mg/day | direct diary total |
| Fiber | 14.1 g/day | direct diary total |
| Protein | 54.6 g/day | direct diary total |

Known saturated-fat energy was 13.2% to 17.5% in every complete calendar year
from 2018 through 2025. This is descriptive and subject to food-database error,
but its persistence makes it unlikely that the final-90 result is a one-week
artifact.

The largest final-90 matched sources were distributed rather than dominated by
one food:

| Food entry | Saturated fat across 90 days |
|---|---:|
| Saag paneer | 167.9 g |
| Heavy cream | 157.5 g |
| White queso | 91.0 g |
| Liege waffles | 89.5 g |
| Bolt Bites chocolates | 70.0 g |
| Peanut noodles | 63.0 g |
| Butter cookies | 54.0 g |
| Shredded cheddar | 51.8 g |
| Peanut-butter almond cups | 51.8 g |
| Cinnamon rolls | 51.0 g |

The practical implication is a portfolio of substitutions, not banning a
single food. Preserve calories and menu stability initially while replacing
some heavy cream/butter/cheese/confection calories with unsaturated-fat and
higher-protein familiar foods. Recheck lipids after an interval chosen with the
clinician.

### Measurements to discuss now

- Fasting or nonfasting lipid panel as the clinician prefers: total
  cholesterol, LDL-C, HDL-C, and triglycerides; calculate non-HDL cholesterol.
- ApoB. The 2026 dyslipidemia guideline describes ApoB as useful when the
  standard profile may underestimate residual particle risk, especially with
  elevated triglycerides, diabetes, or low achieved LDL-C.
- Lp(a) once in adulthood. The same guideline recommends at least one
  measurement.
- Blood pressure, smoking status, family history, kidney function, and glycemia
  so a PREVENT risk estimate can be calculated rather than guessed.
- Coronary artery calcium only if the ordinary risk discussion remains
  uncertain and the result would change treatment. The 2026 guideline supports
  selective use beginning at age 40 in men and 45 in women, not automatic
  scanning: [2026 dyslipidemia guideline](https://professional.heart.org/en/science-news/2026-guideline-on-the-management-of-dyslipidemia/top-things-to-know).

## 2. Glycemia

`FINDINGS.md` calls the subject non-diabetic with “HbA1c ~5.0%,” but no
laboratory report or date is present. The repository's glucosuria exercise asks
whether modifying logged caloric intake improves an energy-balance model. It
cannot measure glucose, urine glucose, or establish current normoglycemia.
Failure of an outcome-fitted calorie perturbation to improve a model is not a
diabetes screen.

The subject is 43 with a latest BMI of 27.8. The USPSTF recommends screening
asymptomatic adults aged 35–70 who have overweight or obesity. Accepted screens
include fasting plasma glucose, HbA1c, or an oral glucose-tolerance test:
[USPSTF recommendation](https://www.uspreventiveservicestaskforce.org/uspstf/recommendation/screening-for-prediabetes-and-type-2-diabetes).

Discuss obtaining **HbA1c plus fasting plasma glucose**. An oral glucose
tolerance test is a second-line possibility if results are discordant or the
clinician has a specific concern. Also clarify why metformin appears in the
medication record, whether it is still taken, and what baseline abnormality or
goal prompted it. The log alone cannot answer those questions.

## 3. Blood pressure, resting heart rate, and sodium

No blood-pressure series was found. That is a larger information gap than most
of the elaborate modeled variables in the repository.

The final-90 diary sodium mean was 2,598 mg/day. This is a direct diary total,
not urinary sodium. Known potassium was only 977 mg/day, but potassium values
covered 79.6% of calories, so the repository cannot quantify total potassium
adequacy.

Use a validated upper-arm cuff and a standardized home protocol agreed with the
clinician—commonly repeated morning and evening readings across seven days—plus
resting pulse. Do not substitute a cuffless smartwatch. The 2025 AHA/ACC
guideline emphasizes home monitoring with standardized protocols and advises
against reliance on cuffless devices until their accuracy improves:
[blood-pressure guideline](https://professional.heart.org/en/science-news/2025-high-blood-pressure-guideline/top-things-to-know).

The derived overnight heart-rate medians were 83.5 bpm before tirzepatide and
84.0 bpm on tirzepatide, a difference of only 0.5 bpm. This is not evidence of a
large drug-associated increase. However, the raw archive is not in the
repository, and the “overnight” window may include wake time. The current FDA
label reports a mean 1–3 bpm increase in pooled trials:
[FDA Zepbound label](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/217806s042lbl.pdf).

## 4. Central adiposity

The same H30 device moved in a favorable direction:

| H30 measure | 2024-10-20 | 2026-02-28 | Change |
|---|---:|---:|---:|
| Weight | 224.7 lb | 203.0 lb | -21.7 lb |
| Fat mass | 83.6 lb | 60.0 lb | -23.6 lb |
| Fat percentage | 37.2% | 29.5% | -7.7 points |
| Lean mass | 141.1 lb | 143.0 lb | +1.9 lb |
| Visceral-fat level | 18 | 11 | -7 levels |
| Waist-to-hip ratio | 1.07 | 0.98 | -0.09 |

This is encouraging direct device evidence, not proof of exact tissue change.
The file has no tape-measured waist and the `visceral_fat_area_cm2` field is
blank. Consumer BIA “visceral-fat level” and waist-to-hip estimates are vendor
outputs, not actual waist circumference or imaging.

Add a standardized tape waist measurement, ideally repeated at the same
anatomic landmark, time, and respiratory state. It is cheap, interpretable, and
less model-dependent than trying to infer visceral risk from the latent-state
pipeline.

## 5. Liver and kidney

No auditable ALT, AST, alkaline phosphatase, bilirubin, albumin, platelet,
creatinine, eGFR, electrolyte, or urine albumin result was found. Consequently,
the repository cannot evaluate fatty liver, fibrosis risk, kidney function, or
electrolyte status.

Discuss:

- a comprehensive metabolic panel and CBC;
- urine albumin-to-creatinine ratio if blood pressure, glycemia, or other kidney
  risk warrants it;
- FIB-4 only if there is a clinical reason to assess metabolic liver risk. It
  requires age, AST, ALT, and platelets; none are available here. AASLD describes
  it as a primary risk-stratification tool, not a diagnosis:
  [AASLD clinical guidance](https://www.aasld.org/liver-fellow-network/core-series/clinical-pearls/spare-me-jab-noninvasive-assessment-patients-masld).

Tirzepatide does not require indiscriminate kidney testing after every dose.
The FDA label specifically says to monitor renal function when adverse effects
could cause volume depletion, especially during initiation or escalation. It
also flags gallbladder disease and pancreatitis symptoms. Persistent severe
abdominal pain, severe or persistent GI symptoms, inability to maintain fluids,
or gallbladder-type symptoms belong with the prescriber promptly; routine
amylase/lipase screening in an asymptomatic person is not proposed here.

## 6. Cold intolerance, thyroid questions, and RMR

The history of severe cold intolerance deserves an ordinary clinical workup
before more thermogenesis modeling or cold-exposure experiments. No TSH, free
T4, CBC, ferritin, iron studies, or B12 result was found.

The latest measured RMR was 1,750 kcal/day on 2025-09-03. The mean of two
FitMate tests on 2024-09-17 and 2024-09-19 was 1,883 kcal/day, an absolute
decline of 7.1% while body weight also fell substantially. This does not by
itself establish abnormal metabolic suppression. The same instrument produced
2,241 and 1,884 kcal/day on consecutive days in June 2023, demonstrating enough
short-term variability to make one-versus-one comparisons unsafe.

The temperature log cannot validate a thyroid diagnosis or caloric threshold.
Measurement site, device, and protocol are absent, readings are opportunistic,
and peripheral/oral temperature is not whole-body energy expenditure.

Discuss:

- TSH and free T4;
- CBC, ferritin, iron/TIBC or transferrin saturation, and B12 as clinically
  appropriate;
- a repeat RMR only if the answer would change management, ideally with two or
  three standardized tests under matched fasting, activity, caffeine,
  temperature, and rest conditions.

The Zepbound boxed warning concerns medullary thyroid carcinoma and MEN2; it is
not an explanation of routine hypothyroidism or cold intolerance. Confirm
personal/family MTC or MEN2 history with the prescriber. The FDA label says
routine calcitonin or thyroid-ultrasound monitoring is of uncertain value in
asymptomatic users.

## 7. Micronutrients and a data-dictionary error

The long-term vegetarian pattern, dairy as the only regular animal food, recent
low protein/fiber intake, metformin entries, and weight-loss treatment make
nutrient review sensible. The repository cannot currently quantify it.

Vitamin B12 is absent from the enrichment schema. Folate, zinc, iodine, and
selenium are also not captured consistently. More importantly, fields named
`calcium_mg` and `iron_mg` appear to contain MyFitnessPal percent-Daily-Value
numbers in many source records, not milligrams. Examples:

- “Beechers Smoked Cheddar, 1 ounce” has `calcium_mg=20`, a plausible 20% DV
  label value but an implausible 20 mg total for an ounce of cheddar.
- “Kroger Tender Spinach, 1.38 cup” has calcium 4, iron 7.5, vitamin A 80, and
  vitamin C 20—the characteristic shape of percent-DV label fields.

`intake/mfp_enrich.py` scales these values unconditionally without a unit
conversion. Therefore calcium, iron, vitamin A, and vitamin C totals should not
be used for adequacy claims until the schema is repaired and records are
unit-normalized. This does not invalidate the script's saturated-fat result,
which uses a gram field, or direct diary sodium/fiber/protein totals.

Discuss CBC, ferritin plus iron/TIBC or transferrin saturation, B12 (with
methylmalonic acid if the clinician finds a borderline result), and 25-OH
vitamin D. Calcium is included in a metabolic panel, but serum calcium is not a
dietary-adequacy test; a dietitian can assess intake after the source fields are
fixed. Do not reflexively supplement every unmeasured nutrient.

Long-term metformin use is associated with B12 deficiency; ADA standards
support periodic B12 assessment, particularly with anemia or neuropathy:
[ADA evidence-based recommendation](https://diabetesjournals.org/care/article/47/Supplement_1/S43/153945/3-Prevention-or-Delay-of-Diabetes-and-Associated).

## 8. Fatty acids

The two same-lab OmegaQuant panels show:

| Analyte | 2024-09-05 | 2025-09-05 |
|---|---:|---:|
| Omega-3 Index | 4.93% | 6.06% |
| EPA | 0.54% | 0.45% |
| DHA | 2.63% | 3.71% |
| Arachidonic acid | 7.52% | 8.93% |
| AA:EPA | 14.0 | 19.8 |
| Omega-6:omega-3 | 7.0 | 6.1 |

The daily 500 mg EPA plus 500 mg DHA record began only 20 days before the 2025
test. That is too soon for a stable long-term interpretation, but not too soon
for any incorporation. A repeat from the same lab after a stable verified
exposure can answer an uptake question. It is lower priority than obtaining
conventional lipids, ApoB, blood pressure, and glycemia, and it will not test a
weight-loss mechanism.

## 9. Medication and supplement reconciliation

### Tirzepatide

The log contains 12.5 mg weekly entries through 2026-03-24, followed by four
“25 click” entries recorded as 5.2083 mg on 2026-04-01, 04-04, 04-07, and 04-10.
The intake and observed-weight outcome series end before this three-day pattern,
so the repository has no outcome validation for it.

The FDA regimen is once weekly. The label's 72-hour language allows changing
the weekly administration day while keeping doses separated; it does not
describe chronic every-three-day administration. Do not operationalize
click-count dosing from this repository. Ask the prescriber and dispensing
pharmacist to verify the exact product, actual delivered dose, interval, and
rationale against the current label:
[FDA prescribing information](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/217806s042lbl.pdf).

### Metformin

There are duplicate/ambiguous 500 mg XR entries on 2025-03-31 and 2025-04-01,
but no reliable stop date, indication, or current regimen. Reconcile it. If it
is chronic, B12 monitoring is particularly relevant given the vegetarian diet.

### Magnesium

The background lists 500 mg/day but does not say whether that is 500 mg of a
compound or 500 mg **elemental magnesium**. Supplement Facts labels declare
elemental magnesium. If it is 500 mg elemental from a supplement, it exceeds the
U.S. adult supplemental upper level of 350 mg/day. High supplemental doses can
cause diarrhea and cramping, and toxicity risk rises with impaired kidney
function: [NIH Office of Dietary Supplements](https://ods.od.nih.gov/factsheets/Magnesium-HealthProfessional/).
Verify label, form, elemental dose, indication, GI effects, and kidney function
with the clinician or pharmacist.

### Creatine

Creatine 5 g/day can modestly increase serum creatinine without a corresponding
fall in measured GFR, so disclose use when kidney labs are interpreted. It is
not evidence that creatine is damaging the kidney. If creatinine-based eGFR is
unexpected or discordant with the clinical picture, the clinician may consider
cystatin C or another confirmatory assessment:
[2026 randomized-trial meta-analysis](https://pubmed.ncbi.nlm.nih.gov/42035842/).

## Prioritized action list

### Priority 1: close safety and risk-assessment gaps

1. Bring the actual medication bottles/pens and Supplement Facts panels to a
   prescriber or pharmacist reconciliation. Resolve the tirzepatide interval
   and click-count entries, metformin indication/current use, elemental
   magnesium dose, and personal/family MTC or MEN2 history.
2. Obtain standardized upper-arm home blood pressure and resting pulse, plus an
   actual tape waist.
3. Discuss a lipid panel, ApoB, Lp(a), HbA1c, and fasting glucose.
4. Discuss CBC and comprehensive metabolic panel, including creatinine/eGFR,
   electrolytes, calcium, and liver enzymes.
5. Because of cold intolerance and vegetarian/metformin context, discuss TSH,
   free T4, ferritin/iron studies, B12, and 25-OH vitamin D.

### Priority 2: use results to choose targeted follow-up

- Calculate PREVENT risk only after lipids and blood pressure exist.
- Consider urine albumin/creatinine if blood pressure, glycemia, or kidney risk
  warrants it.
- Consider cystatin C if creatinine-based eGFR is unexpectedly low or difficult
  to interpret while taking creatine.
- Calculate FIB-4 only when AST, ALT, and platelets exist and a liver-risk
  assessment is clinically relevant.
- Assess sleep apnea clinically if snoring, witnessed apneas, morning headache,
  nonrestorative sleep, or daytime sleepiness is present; the repository cannot
  infer it from weight alone.
- Consider CAC only if the clinician-patient lipid decision remains uncertain
  after ordinary risk factors are measured.

### Priority 3: interesting but lower-yield refinements

- Repeat OmegaQuant after a stable verified supplement exposure.
- Repeat standardized RMR only if it would change the plan.
- Consider DXA only if precise body composition would materially change
  training or nutrition management; otherwise strength outcomes, tape waist,
  and same-device trends are sufficient.

## Concrete clinician questions

1. “The repo has no source lab reports. Can we establish a baseline with
   lipids, ApoB, one-time Lp(a), HbA1c/fasting glucose, CBC, CMP, and the
   symptom-directed thyroid/iron/B12/vitamin D tests?”
2. “Does the recorded every-three-day, click-count tirzepatide pattern match
   the prescribed product and intended regimen?”
3. “What was metformin prescribed for, and am I still supposed to take it?”
4. “Is the magnesium label 500 mg elemental magnesium, and is that dose
   appropriate given GI effects and kidney function?”
5. “Does creatine make creatinine-based eGFR hard to interpret in my case, and
   would cystatin C be useful if the result is unexpected?”
6. “Given the high known saturated-fat pattern, what lipid response would make
   us intensify food substitutions or discuss medication?”
7. “What standardized BP and waist protocol should I use at home?”
8. “Could cold intolerance merit thyroid, anemia/iron, or B12 evaluation before
   further metabolic testing?”

## Reproducer and artifacts

Run:

```bash
uv run --with pandas python analysis/EH_cardiometabolic_snapshot.py
```

It writes:

- `analysis/EH_clinical_snapshot.csv` — direct and derived snapshot rows with
  source and evidence notes;
- `analysis/EH_nutrition_year_summary.csv` — annual direct diary measures and
  known saturated-fat estimates with coverage;
- `analysis/EH_top_saturated_fat_sources.csv` — final-90 food sources.

The script has no random component, takes no command-line arguments, and does
not modify source data.

## Limits

- This audit can identify missing tests and risk-relevant patterns; it cannot
  diagnose dyslipidemia, diabetes, hypertension, liver disease, kidney disease,
  thyroid disease, anemia, or a nutrient deficiency.
- Food-database enrichment has match and unit error. Known saturated fat is
  more usable than the mislabeled calcium/iron/vitamin fields, but it remains an
  estimate.
- BIA directional trends are more credible than exact tissue masses.
- Medication logs show recorded events, not necessarily the prescribed or
  currently followed regimen.
- Clinical test selection and interpretation belong with a clinician who knows
  symptoms, family history, examination, and the full medication list.
