# Claim registry

Date: 2026-07-30

This is the decision boundary for the project. “Supported” does not mean universal or causal unless explicitly stated.

## Status definitions

| Status | Meaning |
|---|---|
| **Directly supported** | Follows from raw or minimally transformed observations; causal interpretation may still be limited |
| **Externally supported** | Supported by relevant trials, meta-analyses, guidelines, or labels; subject-specific effect still needs testing |
| **Exploratory** | Plausible signal, but confounding, leakage, multiplicity, sparse data, or lack of holdout prevents reliance |
| **Contradicted** | The stated interpretation conflicts with its own result, independent evidence, or the underlying code/data |
| **Quarantined** | Should not guide action because the downside of error is material |

## Data and model claims

| Claim | Status | Current assessment |
|---|---|---|
| Every intake day from 2011-04-20 through 2026-02-28 is represented | **Directly supported** | 5,429 consecutive dates; zero gaps |
| Daily calorie totals reconcile with food items | **Directly supported** | All dates reconcile within 1 kcal |
| The ledger contains 54,925 food rows | **Contradicted** | Current canonical file contains 54,453 rows |
| Logging error is uniformly about 9–14% | **Exploratory** | Depends on a model with sparse observed weights, intake-derived TDEE, and uncertain composition interpolation |
| Cumulative energy balance closes to ±5 lb and proves intake accuracy | **Exploratory** | Closure is model-dependent; it is not independent validation when TDEE/latent states are tuned using the same channels |
| P3 is an intake-free fat-mass confirmation | **Contradicted** | P3 derives TDEE from intake and propagates missing weight using energy balance |
| P4 RMR and temperature are held out | **Contradicted** | RMR enters the initial TDEE prior and both channels enter hyperparameter scoring |
| P4 `*_filtered` columns are usable forward estimates | **Contradicted** | They contain retrospective RTS-smoothed states |
| P4 exports usable uncertainty | **Contradicted** | All FM/TDEE standard-deviation columns are blank |
| Daily P4 FM/TDEE/LM can validate downstream causal claims | **Exploratory** | Sixty-seven scripts consume a leaky, retrospectively smoothed model |
| P4b has genuine forward-filtered and separate smoothed semantics | **Directly supported, implementation-level** | Eleven deterministic tests verify input isolation, no future leakage to `1e-12`, distinct filter/smoother outputs, complete covariance, and no double-counted scan weight |
| P4b v1 is ready to replace P4 downstream | **Contradicted** | Later-period intervals are substantially over-wide; locked FM MAE was 1.83 lb versus 1.78 lb for a simple weight-adjusted persistence baseline |
| P4b `tdee_logged_kcal` is independently measured physiological TDEE | **Contradicted** | It is an intake-conditional latent accounting state under a fixed 1.0 logged-intake multiplier and sparse observation channels |

## Set point and behavior

| Claim | Status | Current assessment |
|---|---|---|
| Weight loss invokes homeostatic appetite and expenditure responses | **Externally supported** | Strong general evidence; magnitude and architecture vary |
| A set point updates only after 14 stable days within ±3 lb | **Contradicted as established fact** | Selected in-sample across many configurations; no locked trial holdout |
| Adaptation occurs at exactly 1%/day | **Exploratory** | Not independently identified |
| Appetite pressure is exactly 55 kcal/day per lb of fat-mass gap | **Contradicted as established fact** | Derived from shared modeled inputs; sign/mechanism problems remain |
| The 45–50 day timescale converges from four independent measurements | **Contradicted** | Purported confirmations share intake-derived, smoothed, or outcome-linked inputs |
| Holding weight for 6–12 months makes medication withdrawal safe | **Quarantined** | Not supported by the repo or withdrawal trials |
| Tirzepatide-era intake is lower and binge persistence is reduced | **Directly supported, descriptive** | Strong within-person era difference; treatment causality is plausible but not isolated from calendar and behavior changes |
| Tirzepatide-era intake benefit has fully worn off | **Contradicted** | Last-90 mean intake remained about 450 kcal/day below the pre-drug-year mean; ≥3,000-kcal days were 1.1% vs 11.0% |
| Tirzepatide suppresses exactly 74.5 kcal per modeled blood-level unit | **Exploratory** | Injection day is nearly identical to weekday and multiple time trends overlap |
| Cumulative receptor exposure causes the fitted tachyphylaxis constant | **Exploratory** | Outcome-fitted construct; not a measured receptor state |
| The observed weekly intake sawtooth proves a PK dose-response | **Exploratory** | Cycle day maps to weekday with 99.1% purity; only one weekly injection was not Tuesday |
| The q3d experiment has validated the model | **Contradicted** | Canonical intake and observed weight end before q3d begins |
| Potato-template days had lower intake/no binges in the annotated periods | **Directly supported, descriptive** | All four attempts averaged 524–748 kcal/day below their immediate pre-windows; each of the four measurable immediate post-windows regained 3.3–3.9 lb, so rapid scale loss is not equivalent to fat loss |
| The potato experiments prove a potato-specific metabolic advantage | **Contradicted as established fact** | Calories, protein, carbohydrate, food variety, sodium, food mass, added fats, and water/gut contents changed together; the raw data identify a repeatable low-intake menu, not its mechanism |
| Weekend 36-hour fasting produced durable loss without compensation | **Contradicted** | Fourteen fast days lowered the 44-day mean to 1,624 kcal/day and measured weight fell 4.2 lb, but the next 28 days averaged 2,364 kcal and observed weight rose 7.3 lb |
| Keto produced demonstrated weight loss in the annotated phases | **Exploratory/open** | Neither keto phase contains an observed scale weight; phase 1 intake was 440 kcal/day above its local pre-window and phase 2 was unchanged |
| The low-protein experiment reduced intake or weight | **Contradicted** | Thirty-four of 37 days met the direct low-protein rule; intake averaged 2,340 kcal on those days and observed weight rose 1.3 lb across the measured epoch |
| Enlightened-brand dessert days prove processed foods reduce appetite | **Exploratory** | Product days were 154 kcal/day lower than other days in the era, but product choice was not randomized and only one scale weight exists in the 426-day epoch |
| Domino's-template days are high-intake/high-binge days | **Directly supported, descriptive** | Likely useful context marker; does not identify cheese, oil, refined carbohydrate, sodium, or ordering as the cause |
| Travel periods are higher-intake contexts | **Directly supported, descriptive** | In a 41-trip paired study, 36/41 trips exceeded their own prior 14-day intake; mean +605 kcal/day and median +632 |
| Travel reliably causes fat gain of a known magnitude | **Exploratory** | Intake effect is strong, but sparse off-trip weigh-ins and water/gut shifts prevent direct fat-gain attribution |
| A highly repeated recent menu marks lower extreme-intake risk | **Directly supported, descriptive** | Direction replicated in discovery, validation, confirmation, all seven calendar years, and a home-only sensitivity |
| Food novelty marks higher extreme-intake risk | **Directly supported, descriptive** | Direction replicated across ordered splits; much of the mean-calorie effect is travel/context |
| Menu monotony itself causally suppresses appetite | **Exploratory** | Reverse causality and intentional restriction remain; needs prospective ABAB replication |
| Major observed loss/regain reversals coincide with sustained logged-intake regime shifts | **Directly supported, descriptive** | Three clean loss-to-regain reversals moved from 1,962→2,233, 1,797→2,326, and 2,103→2,366 kcal/day, while the pre-treatment-to-treatment transition moved from 2,400→1,999 |
| A 28-day intake/context score reliably predicts next-quarter weight gain | **Not supported by the internal locked test** | None of intake, high-intake frequency, travel, repetition, strength frequency, or weighing history preserved its expected development-period direction in 14 later pre-treatment quarters |
| Very long weigh-in gaps are uninformative or missing at random | **Contradicted** | All five gaps of at least 90 days ended at a higher observed weight, median +11.6 lb; this is an informative process marker, not evidence that weighing prevents gain |
| An extreme-intake day often continues into the next day before treatment | **Directly supported, descriptive** | At home, 22.2% of pre-treatment extreme days were followed by another versus 5.4% after other days; adjusted next-day difference was +359 kcal and +16.2 percentage points |
| Recent calorie momentum is a deployment-ready extreme-intake risk score | **Exploratory** | One-day-ahead ranking improved in two held-out pre-treatment periods, but probability calibration and flagged-day prevalence shifted across eras |
| Weekends are a stable causal driver of high intake | **Not supported** | The home weekend difference was +62 kcal/day in discovery, +4 in validation, and +92 in confirmation; treatment-era weekday patterns are inseparable from injection weekday |
| Higher day-to-day intake variance is mildly protective at the same calorie surplus | **Contradicted as established** | AF pairs trailing intake variance/surplus with a non-overlapping future modeled-FM interval; on the aligned same window, partial r is +0.003 at 14 days and -0.006 at 30 days |
| The current menu has no obvious calorie-neutral nutrition leverage | **Contradicted** | In the last 90 days, 45.6% of calories came from recurring foods below 3 g protein/100 kcal; breakfast and snacks were especially protein-light |
| Adding protein necessarily requires adding calories | **Contradicted** | About 914 kcal/day currently comes from recurring low-protein sources, while familiar 7–13 g/100-kcal alternatives already exist in the log |
| Whole-diet energy density can be reconstructed from the historical ledger | **Contradicted** | Exact mass covers 17.0% of solid-food calories overall and 10.8% in the final 90 days, with nonrandom missingness |
| Protein/fiber density, beverage share, or item count reliably predicts future intake | **Not supported by locked tests** | No prespecified feature added stable, practically meaningful next-day prediction across both ordered pre-treatment holdouts |
| Breakfast and snacks are the clearest current nutrient-allocation targets | **Directly supported, descriptive** | They supply about 31% of final-90 calories and have the lowest protein and fiber density; this supports substitution but does not prove a satiety mechanism |

## Medication and intervention claims

| Claim | Status | Current assessment |
|---|---|---|
| Long-term anti-obesity pharmacotherapy can produce substantial weight loss | **Externally supported** | Strong randomized-trial and guideline evidence |
| Stopping GLP-1/GIP therapy often leads to regain | **Externally supported** | Seen in SURMOUNT-4 and STEP withdrawal/extension data |
| Switch to 25 clicks every three days and bridge pens | **Quarantined** | Off-label, unvalidated, contrary to intended fixed-dose pen use; requires prescriber/pharmacist decision |
| `>58 g/day` is the subject's validated protein threshold | **Exploratory** | A dichotomized in-sample association; no causal or held-out support |
| Adequate protein plus resistance training can help preserve lean mass | **Externally supported** | Broad evidence; exact target must be individualized |
| Current treatment-era protein support is clearly adequate | **Contradicted** | Last-90 mean was 54.6 g/day; 60/90 days were below 60 g and only 4/90 reached 80 g |
| Tirzepatide caused major observed lean loss after the H30 device began | **Not supported by available same-device data** | Twenty H30 scans show no obvious decline from day 33 onward, but creatine restarted mid-series and BIA is hydration-sensitive |
| Central-adiposity indicators improved during treatment | **Directly supported, device-derived** | Same-device H30 visceral-fat level fell 18→11 and waist-to-hip ratio 1.07→0.98; tape-measured waist is unavailable |
| Weight loss was exactly 93% fat | **Exploratory** | Same-device H30 gives an implausible 108.8% measured fat share, illustrating BIA error; first treatment month is missing |
| Each walk session raises RMR by exactly 14 kcal/day | **Exploratory** | Only 25 clustered RMR observations; model selection was not nested |
| High steps during restriction cause a metabolic penalty | **Contradicted as established** | U and V miscompute nominal future windows as mostly trailing windows; the corrected direct analysis finds no overall future-intake or observed-weight penalty, with only a small treatment-era intake association left exploratory |
| Sleep does not matter for this subject | **Contradicted as a general conclusion** | The reliable 2022+ regime averages 8.03 h and has only six under-5 h nights; longer-than-usual sleep has a small lower-intake association but negligible holdout value, so maintenance is reasonable and aggressive optimization is low priority |
| A later bedtime lowers metabolism in this subject | **Exploratory** | The legacy endpoint is P4-derived and window-selected; direct intake results do not support earlier bedtime as a weight-loss lever |
| Calculated wake-to-sunset hours measure sunlight exposure | **Contradicted** | `sunlight.py` estimates possible available daylight, not outdoor exposure, illuminance, or dose |

## Genetic, fatty-acid, and speculative claims

| Claim | Status | Current assessment |
|---|---|---|
| The subject has a “worst-case common-variant obesity genotype” | **Contradicted** | Linked SNPs are double-counted; no validated ancestry-aware polygenic score |
| FTO common variants modestly raise average BMI/appetite risk | **Externally supported** | Small, non-deterministic effect; does not prescribe a special diet |
| Common MC4R-region tags imply receptor dysfunction | **Contradicted** | These are common nearby GWAS markers, not pathogenic MC4R deficiency |
| UCP2/FTO genotype proves low RMR or cold intolerance | **Contradicted** | Mechanistic extrapolation unsupported |
| GLP1R rs6923761 predicts tirzepatide response | **Contradicted as clinically actionable** | No consistent pharmacogenomic effect or genotype-based dosing guidance |
| IRS1 rs2943641 T is the insulin-resistance risk allele | **Contradicted** | The risk direction is reversed; C is the reported risk allele |
| APOA2 heterozygosity creates saturated-fat sensitivity | **Contradicted** | Foundational interaction is recessive in CC homozygotes |
| FADS variants affect PUFA profiles | **Externally supported** | Biomarker association does not prove inflammation or obesity causation |
| Seed oils/linoleic acid are a major demonstrated cause here | **Contradicted as established** | Human trials do not show the claimed inflammatory effect; repo PUFA coverage is inadequate |
| High linoleic acid generally raises inflammatory markers | **Contradicted** | Randomized-trial reviews do not support the claim |
| Ultra-processed food can increase ad-libitum intake | **Externally supported** | Controlled inpatient evidence; exact mechanism and subject effect remain open |
| Fazzino fat+sugar exposure marks higher same-day intake in this subject | **Directly supported, descriptive** | Across ordered discovery, validation, and confirmation splits, +10 percentage points of assessable FS calories marked about +50, +114, and +66 same-day kcal after calendar/context controls; the exposure is part of the same day's calorie total, so causality is not identified |
| Fazzino fat+sugar exposure reliably predicts next-day extreme intake | **Not supported by the current direct test** | Adjusted ≥3,000-kcal probability effects were +0.3, +0.4, and +0.2 percentage points per +10 percentage points FS exposure, with every 95% interval crossing zero; home-only and data-quality sensitivities remain null |
| The repository measures total Fazzino-style hyper-palatable-food exposure | **Contradicted as a complete measurement** | FS has 94.6% calorie-weighted coverage in the modern pre-drug analysis, but sodium clusters require food weight and have only 17.0% exact-mass coverage; missing foods cannot be called non-HPF |
| Cold exposure is a meaningful genotype-directed weight-loss intervention | **Quarantined** | Weak weight-loss evidence, compensatory appetite possible, and poor fit for severe cold intolerance |
| Taurine 1–3 g is a genetically indicated weight-loss treatment | **Contradicted** | Preclinical rationale without demonstrated human weight-loss benefit |
| Magnesium should be timed to the highest-carbohydrate meal because of genotype | **Contradicted** | No validated genotype-based rationale; elemental dose may itself need safety review |
| A hand-selected combination of common SNPs gives a meaningful combined prevalence | **Contradicted** | Linked loci were double-counted and no ancestry-calibrated polygenic score was used |
| The 2025 Omega-3 Index rise was entirely dietary because supplementation had only run 20 days | **Contradicted** | Too early for full RBC steady state is not the same as too early for any incorporation |
| Microbiome, circadian timing, or food-order effects are impossible | **Exploratory/open** | Not ruled out; require better-defined prospective or time-ordered tests |

## Clinical and cardiometabolic claims

| Claim | Status | Current assessment |
|---|---|---|
| The repository establishes a current normal conventional clinical baseline | **Contradicted** | No auditable current lipid panel, ApoB, glucose/HbA1c, blood pressure, CBC/CMP, kidney, liver, or thyroid result was found |
| The subject's current HbA1c is about 5.0% and the repository establishes that he is non-diabetic | **Contradicted as established fact** | The value appears only in prose without a date or source report; the energy-balance glucosuria perturbation is not a glucose measurement or diabetes screen |
| Meaningful glucosuria was ruled out with 95% confidence | **Contradicted as a clinical claim** | Re-fitting an intake/energy-balance model after hypothetical calorie removal cannot observe glucose or urine glucose and inherits the latent model's identification problems |
| Recent saturated-fat intake is within a heart-healthy range | **Contradicted** | In the final 90 days, matched foods supplied 33.5 g/day, or 15.0% of logged energy, with 90.3% calorie coverage; current AHA guidance describes heart-healthy patterns as unlikely to exceed 10% |
| `calcium_mg` and `iron_mg` can be summed as milligrams to assess adequacy | **Contradicted** | Many underlying MyFitnessPal records contain percent-Daily-Value numbers in these fields; the enrichment pipeline scales them without unit normalization |
| One RMR reading of 1,750 kcal/day proves abnormal metabolic suppression | **Not supported** | Sparse calorimetry has large short-term variability, including 2,241 and 1,884 kcal/day on consecutive 2023 tests; matched repeat testing and ordinary clinical evaluation are needed |
| Current H30 central-adiposity indicators improved | **Directly supported, device-derived** | H30 visceral-fat level fell 18→11 and estimated waist-to-hip ratio 1.07→0.98, but no tape waist or imaging-derived visceral-fat area is available |
| The recorded every-three-day tirzepatide pattern is an FDA weekly regimen | **Contradicted** | The label specifies once-weekly use; its 72-hour language permits changing the weekly administration day, not chronic every-three-day dosing |

## Sources for the highest-stakes corrections

- [FDA Zepbound prescribing information, 2026](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/217806s042lbl.pdf)
- [2026 AHA dietary guidance](https://professional.heart.org/en/science-news/2026-dietary-guidance-to-improve-cardiovascular-health/top-things-to-know)
- [2026 dyslipidemia guideline](https://professional.heart.org/en/science-news/2026-guideline-on-the-management-of-dyslipidemia/top-things-to-know)
- [USPSTF diabetes screening recommendation](https://www.uspreventiveservicestaskforce.org/uspstf/recommendation/screening-for-prediabetes-and-type-2-diabetes)
- [Mounjaro KwikPen patient leaflet](https://www.medicines.org.uk/emc/product/15483/pil)
- [SURMOUNT-4 randomized withdrawal trial](https://pubmed.ncbi.nlm.nih.gov/38078870/)
- [STEP 1 extension after semaglutide withdrawal](https://pubmed.ncbi.nlm.nih.gov/35441470/)
- [FTO variant mechanism, Claussnitzer et al.](https://pubmed.ncbi.nlm.nih.gov/26287746/)
- [FTO and intervention response meta-analysis](https://pubmed.ncbi.nlm.nih.gov/27650503/)
- [Linoleic acid and inflammatory markers review](https://pubmed.ncbi.nlm.nih.gov/28752873/)
- [Ultra-processed diet inpatient crossover trial](https://pubmed.ncbi.nlm.nih.gov/31105044/)
- [Joint nutrition advisory for GLP-1 therapy](https://pmc.ncbi.nlm.nih.gov/articles/PMC12264624/)
