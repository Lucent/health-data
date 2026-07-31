# Initial integrity and safety audit

Date: 2026-07-30
Status: first checkpoint; living report, not medical diagnosis or individualized prescribing

## Bottom line

The repository contains an unusually valuable 5,429-day uninterrupted food ledger. That makes it capable of answering questions that ordinary diet recalls cannot. It does **not** validate the repo's flagship claims about a stability-gated fat-mass set point, exact hunger pressure, tirzepatide tachyphylaxis, or every-three-day click dosing.

The most reliable current weight-loss inference is simpler:

1. Tirzepatide-era intake and binge persistence are descriptively lower, but the exact `-74.5 cal/unit` coefficient and modeled receptor erosion are not identified.
2. Long-term anti-obesity pharmacotherapy, when clinically appropriate and supervised, has much stronger external evidence than the repo's proposed discontinuation timing or off-label dosing schedule.
3. A repeatable, minimally processed, protein-adequate food environment and resistance training are sensible supports for appetite control and lean-mass retention. The repo-specific `>58 g/day` threshold is not validated.
4. Travel, order-in meals, and the transition out of unusually restrictive food templates look more actionable than steps, sleep duration, a genetic “thermogenesis defect,” or seed oils. These are priority targets for prospective testing, not yet causal conclusions.

## What is genuinely strong

| Source | Directly observed coverage | Audit result |
|---|---:|---|
| Daily intake | 5,429 consecutive days, 2011-04-20 to 2026-02-28 | Zero calendar gaps |
| Food items | 54,453 rows | Item calories reconcile with daily totals within 1 kcal on every date |
| Scale weight | 1,693 days, 1999-05-17 to 2026-03-22 | Useful but sparse; 71.6% of intake-period days lack an observed weigh-in |
| Body composition | 70 scans | Six device regimes and gaps up to 1,730 days |
| Indirect calorimetry | 25 tests | Valuable, but too few and clustered for aggressive feature selection |
| Steps | 4,592 days | 4,327 Samsung totals plus 265 non-comparable exercise/shift estimates |
| Sleep | 2,057 days | Only 60.5% coverage across its span; several retained “primary sleeps” are likely naps |
| Strength | 437 dates through 2026-07-07 | Presence only; all sessions normalized to 30 minutes |
| Fatty-acid labs | 70 analyte rows | Only three blood-draw dates |

The audit script is `analysis/EA_repository_audit.py`:

```bash
uv run --with pandas python analysis/EA_repository_audit.py
```

## Why the flagship model is not validated

### P3 is not intake-free

`analysis/P3_interpolate_weight.py` derives TDEE from intake and uses energy-balance propagation to fill long periods without an observed weight. It therefore cannot independently validate an intake/set-point model. Of 5,429 intake days, 3,885 have no observed weight. Calling the resulting daily series “weight interpolation alone” materially overstates the evidence.

### P4 leaks validation channels and future information

`analysis/P4_kalman_filter.py`:

- initializes TDEE from `expected_rmr`, which is fitted to the 25 RMR measurements;
- uses RMR and temperature statistics in the hyperparameter composite score even while calling them held-out validation;
- writes retrospective RTS-smoothed states into both ordinary and `*_filtered` columns, so the supposed filtered series sees future observations;
- exports four uncertainty columns that are entirely blank;
- partitions modeled tissue change using workout density, after which downstream scripts test workouts against the modeled lean-mass result.

At least 67 analysis scripts consume P4. Their results are not all false, but predictive, causal, or “independently confirmed” claims need to be rerun using a repaired pipeline.

### The set-point result has a sign and identification problem

The repo's own `FINDINGS.md` records that the binge/set-point association has the opposite sign from the proposed mechanism. The analysis then forces a positive-slope sigmoid. Other purportedly independent confirmations share intake-derived or smoothed inputs, while the fixed-set-point comparison is not identified by a correlation that is invariant to an additive constant.

The exact SmoothLatch rule — 14 stable days, ±3 lb, 1% adaptation/day, and 55 kcal/day/lb pressure — was selected across many configurations without a locked held-out trial set. It predicts only about 14% recovery of the original loss for this subject, while withdrawal trials commonly show much larger recovery. It is an interesting hypothesis, not a basis for discontinuation planning.

### The drug model is heavily confounded

Eighty of 84 weekly injections occurred on Tuesdays, so injection-cycle day is also weekday. Week fixed effects cannot fully recover a treatment-cycle effect when almost every injection has the same weekday. Dose escalation, calendar time, cumulative exposure, antihistamine changes, metformin, weight loss, and food-pattern changes also overlap.

The q3d period begins after the intake ledger ends, so the preregistered prediction cannot yet be evaluated from the canonical daily intake data.

## Immediate medication safety boundary

The every-three-day click-count protocol in the legacy `RECOMMENDATIONS.md` is quarantined.

- The current FDA Zepbound label specifies once-weekly dosing and describes a 72-hour interval only for changing the weekly dosing day, not for converting to q3d dosing.
- The Mounjaro KwikPen patient leaflet describes four fixed doses, says not to use leftover medicine, and says not to change dose unless instructed by a doctor.
- The repo's click arithmetic and pen-bridging instructions are therefore both unvalidated and contrary to the manufacturer's intended use.

Any schedule, dose, discontinuation, or pen-handling change should be decided with the prescriber and dispensing pharmacist. The useful question for them is whether appetite recurrence, adverse effects, access, or treatment goals warrant an **approved** dose adjustment or a different clinically supported maintenance strategy.

Primary sources:

- [FDA Zepbound prescribing information, 2026](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/217806s042lbl.pdf)
- [Mounjaro KwikPen patient leaflet](https://www.medicines.org.uk/emc/product/15483/pil)
- [WHO guideline announcement on GLP-1 medicines and comprehensive long-term obesity care](https://www.who.int/news/item/01-12-2025-who-issues-global-guideline-on-the-use-of-GLP-1-medicines-in-treating-obesity)

## Data limitations that will otherwise create false discoveries

- The food ledger changes source regimes: HTML, OXPS, and MHTML. OXPS renders unknown nutrients as zero, creating false historical nutrient shifts.
- Exact duplicate food rows include plausible repeated servings and must not be mechanically deleted.
- About 10.4% of item rows remain unresolved in the MFP enrichment. Polyunsaturated-fat coverage is only 46.9% of rows and about 42% of logged calories; seed-oil or omega-6 conclusions are premature.
- Body-composition devices change six times. The H30 regime begins only a month after tirzepatide, confounding treatment with device.
- Pre-Samsung “steps” are logged exercise or estimated hospital shifts, not total daily movement.
- `diet_epochs.csv` has 27 overlapping annotations and uncovered intervals; it is not a full partition.
- Several important raw sources are external or absent from the repo, including Samsung exports, bath-meter data, Fitmate DBFs, and Google Timeline.

## Practical working strategy while the models are rebuilt

This is the provisional strategy most consistent with both the direct data and current evidence:

1. **Keep effective pharmacotherapy clinically supervised.** Do not infer that weight stability has made withdrawal safe. Withdrawal studies show substantial regain after semaglutide or tirzepatide cessation: [SURMOUNT-4](https://pubmed.ncbi.nlm.nih.gov/38078870/) and the [STEP 1 extension](https://pubmed.ncbi.nlm.nih.gov/35441470/).
2. **Protect lean mass without worshipping a repo-specific threshold.** Track actual resistance-training volume and progression, not a 30-minute yes/no placeholder. Set protein intake with a clinician or dietitian using dietary pattern, target/adjusted weight, kidney status, tolerance, and medication side effects. See the [joint obesity-society nutrition advisory for GLP-1 therapy](https://www.obesity.org/wp-content/uploads/2025/06/Joint_Nutritional-Priorities-to-Support-GLP-1-Therapy-for-Obesity.pdf).
3. **Design the food environment around low-friction satiety.** Ultra-processed diets caused higher ad-libitum intake and weight gain in a controlled inpatient crossover trial: [Hall et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31105044/). The repository can test whether its own high-risk meal templates predict later intake under matched calories, medication, weekday, and travel.
4. **Precommit around travel and ordering.** Legacy direct summaries repeatedly identify travel and Domino's-style order-in templates as high-intake contexts. Test substitutions prospectively rather than attributing the effect to a single ingredient.
5. **Use exercise for health and lean-mass retention, not as a precise calorie-offset formula.** The repo's `+14 kcal per walk session` estimate is too fragile. Resistance exercise has broader evidence for preserving lean mass during weight loss: [systematic review](https://pubmed.ncbi.nlm.nih.gov/35191588/).
6. **Do not overread the sleep null.** Median sleep is already relatively long, while historical coverage is sparse and contaminated by nap-like records. Sleep extension reduced intake in short sleepers in a randomized trial, but that result may not apply here: [Tasali et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35129580/).
7. **Deprioritize genotype-directed cold, taurine, magnesium timing, and seed-oil avoidance as weight-loss levers.** The genetic claims are mostly small associations or unsupported extrapolations, and the nutrient data is not ready to isolate linoleic acid.

## Next analysis gates

1. **Completed:** build a provenance-aware daily table with observed/source/device flags, time since nearest measurement, and uncertainty tier. See `analysis/EG_observation_provenance.py`.
2. Replace P4 with separately exported forward-filtered and retrospectively smoothed states; export covariance.
3. Lock a temporal holdout before selecting hyperparameters or food hypotheses.
4. Re-estimate direct, pre-specified outcomes: next-day and next-7-day intake, binge onset/persistence, weight trend, and lean-mass change where same-device scans exist.
5. Test temporal food antecedents, menu diversity/entropy, order-in/travel contexts, and medication changes before pursuing molecular ingredient stories.
