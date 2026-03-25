# Theories to test

Organized by what data is needed. "Ready now" means the extracted CSVs are sufficient.

## Ready now (intake_daily.csv + weight.csv)

**Glycogen-water smoothing.** Derive a formula to smooth daily weight perturbations from glycogen-water binding. Each gram of carbs retains ~3g of water, depleted over ~2 days of restriction, recaptured within hours of refeeding. Validate against the Oct-Nov 2019 weekend fasts (weight drop should far exceed caloric deficit, with immediate bounce-back on refeeding). This is prerequisite infrastructure for everything below — without it, daily weight noise drowns real signal.

**Hidden set point estimation.** Find long runs (months) of stable weight (±3 lbs) where intake does not support no change. This is direct evidence for a set point that can only drift under certain circumstances. Extract the set point as a slow-drifting latent variable — Kalman filter or similar state-space model rather than the breakpoint approach in analysis-bmr-fit (which didn't work well).

**Variable RMR trendline.** Fit intake vs. smoothed weight with a slow-changing RMR variable. Find what high/low RMR correlates with — a macronutrient, a ratio, a pattern. Graph derived RMR over time. The 3 indirect calorimetry measurements (2011, 2012, 2016) are anchor points.

**Set point shift conditions.** Look for the set point abruptly increasing or decreasing and examine the days/weeks when it appeared to happen. Can only multi-day binges shift set point up (e.g. 3 consecutive days of 3000 calories)? Do the same calories spread out cause no shift? What conditions push it down — are there magic months where 3 lbs are lost and not regained for a long time?

**Set point and protein.** Does set point only move up in the presence of protein, to capture high-effort hunted meals? Studies show 20% protein diets result in same intake as 10% protein + 500 calories. Test whether high-protein days are more likely to precede set point increases.

**Binge prediction and food noise.** What best predicts whether a day will become a binge (>2800 cal)? The [food noise hypothesis](https://lucent.substack.com/p/craving-food-noise) predicts that binge probability correlates with cumulative distance below set point — not yesterday's intake, not a single macronutrient, but the gap between current weight and the weight the body is defending. Test: estimate set point from weight plateau data, compute distance below it on each day, and see if that predicts binges better than any single-day dietary variable. If food noise is a generalized resource-acquisition drive (same circuit as alcohol, gambling — all silenced by GLP-1 agonists), then the trigger is the gap itself, regardless of what was eaten to get there.

**Week-scale intake invariance.** [Preregistered claim](https://bsky.app/profile/lucent.substack.com): week-scale food intake is uncorrelated with emotional state or individual meal interventions. A week where you "cut back" ends up within a few hundred calories of a week where you feasted then compensated without noticing. Test: compare weekly calorie variance to daily variance. If weekly variance is dramatically lower than you'd expect from independent daily draws, there's a homeostatic weekly regulator.

**Tirzepatide dose-response quantification — PARTIALLY ANSWERED.** "Very curious to learn how many mg correspond to kcal subtracted from daily satiety point." The [study 128-OR](https://doi.org/10.2337/db23-128-OR) showed 15mg lowered intake by 900 kcal. At 2.5mg, observed reduction was 200-400 kcal. Firsthand: "each mg shaves 100 kcal/day from satiety." With dose escalation data (2.5→5→7.5→10→12.5mg) and daily intake logged, fit the dose-response curve. Key observation: at 12.5mg, weight stabilized at 202 then bounced to 205 over Christmas and the drug "started working again" at the higher weight — as if each dose corresponds to a specific set point or loss delta, not a fixed kcal reduction.

**Results from pharmacokinetic modeling (2026-03-24):** Modeled blood concentration using FDA PK parameters (t½=5d, Tmax=24h) with superposition of all prior injection curves. Key findings:
- **Blood level → daily intake: r = -0.50** (partial, controlling time trend). The strongest single predictor of daily calorie intake in the entire dataset.
- **Weekly sawtooth**: injection day 1652 cal → trough (day 5) 2220 cal. A 568 cal/day appetite swing directly tracking the PK curve, confirming that appetite suppression is pharmacokinetically mediated, not psychological.
- **Intake model**: `calories = 2345 - 49 × effective_blood_level`. At zero drug: 2345 cal/day. Fresh 12.5mg peak: 1504 cal/day. This gives ~49 cal reduction per arbitrary blood level unit, or roughly **35 cal/day per mg at steady state** (lower than the self-estimated 100 cal/mg, likely because the self-estimate compared peak effect vs pre-drug baseline rather than steady-state average).
- **Tachyphylaxis**: dose effectiveness decays with half-life of 32 weeks. After 20 weeks on the same dose, 65% effective. This explains the firsthand observation of the drug "starting to work again" after a weight plateau — the plateau is tolerance building, and dose escalation resets it partially.
- **Overall reduction**: 456 cal/day (18.6%) from the pre-tirzepatide year. Study 128-OR at 15mg showed 900 kcal reduction; this data at 12.5mg steady-state shows ~534 cal reduction, scaling linearly.
- **Christmas 2025 spike visible**: weeks 16-18 at 12.5mg show 2672, 2793 cal — environmental override of pharmacological suppression, consistent with food noise being a drive that can be overcome by sufficiently strong external cues.

Remaining questions: Does each dose correspond to a specific set point, or does the tachyphylaxis model fully explain the plateaus? The 32-week half-life of effectiveness means the drug never fully stops working — it just asymptotically approaches a reduced effect. At 12.5mg after 1 year, effective level would be ~36% of initial. Is this consistent with the observed weight trajectory?

**Plateau dynamics: resume vs. jump-to-catch-up.** Two models of weight loss on GLP-1s: (1) weight drops, plateaus, then resumes linear loss from the plateau; (2) weight drops, plateaus, then jumps down as if the plateau never happened, resuming the original trajectory. GLP-1 Discord chose model 1, but firsthand experience suggests model 2. The glycogen-water smoothing should distinguish these — if plateaus are water/glycogen masking ongoing fat loss, model 2 is correct and the "jump" is just water finally releasing.

**Cold intolerance as protein-mediated.** "Couple days of low protein and already feeling hypoglycemic and much less cold intolerance." Cold resolved by first bench presses of a workout, corroborating that protein should be eaten before workouts (mTOR activation). The temperature data + daily protein intake can test whether protein intake predicts next-day body temperature.

**Set point ratchet asymmetry.** "A month of consistent overeating ratchets your set point up a couple pounds that never, ever come off no matter how hard you fight." Test: find periods of sustained overeating (>2800 cal/day for >5 consecutive days) and measure whether the subsequent weight floor permanently increases. Is the ratchet truly one-way, or can sustained undereating ratchet it down? "Weight gain not being reversible by simply removing the element that caused it would be a devastating violation of an expected symmetry."

**Dead zone.** "I have a dead zone between 2000-2500 kcal where my metabolism adapts to maintain." Test: find periods of sustained intake in this range and compute apparent RMR. Does it compress toward intake, as if the body adjusts expenditure to match? Compare with periods above or below this range.

**Antihistamine contribution.** [H1 antihistamines associated with obesity: 10 kg difference (NHANES)](https://doi.org/10.1038/oby.2010.176).

Timeline:
- **2011–2024-12**: Daily 10mg cetirizine (Zyrtec) mornings. Entire dataset under this influence.
- **~2024-12-17**: Started weaning. Switched to 5mg cetirizine at night.
- **2024-12-31**: "Been on 5mg of Zyrtec at night instead of 10mg in the morning for 2 weeks now and my quality of life is down about 25% from general malaise, but the 10 kg difference in this study is incredible."
- **2025-01-25**: "After weening down, if I had to estimate the difference between 10mg cetirizine in the morning and 5mg at night, it is almost an additional 2mg of serum tirzepatide."
- **2025-04-06**: Amazon purchase of levocetirizine (generic Xyzal). Switch from cetirizine to levocetirizine complete.
- **2025-06-18**: Costco Xyzal 110ct restock ($37.99).
- **2025-07-07**: Confirmed nightly routine: "2.5mg (halved) levocetirizine."

Net change: 5mg active enantiomer mornings → 2.5mg active enantiomer nightly.

Confounded by tirzepatide dose escalation (7.5→10mg ~2025-05-27) in the same period. Open question from Daniel Quinn/Ishmael framing: does Zyrtec "artificially" increase appetite, or does allergy misery decrease it? Which is baseline?

**Flavorless oil calorie undercounting.** "Maybe we have a running counter that uses flavor (to estimate macros) × chew/crunch/swallow to tally calories and ancestrally absent flavorless oil undershoots the calculation." If cooking oil is systematically more fattening per calorie than flavorful foods, days heavy in oil-cooked foods should show more weight gain than their calorie count predicts. Testable by classifying food items as oil-heavy vs not.

**Advanced glycation end-products (AGEs).** A reductionist-friendly candidate for quantifying food "badness" beyond NOVA/UPF classification. "Seems way more likely than hyperpalatability" since GLP-1 users aren't especially tempted by hyperpalatable foods. Would require classifying food items by AGE content.

**Breakfast/lunch front-loading.** Does eating over 1000 calories for breakfast and lunch combined predict a lower total for the day?

**Fiber satiety.** Fiber consumed in breakfast+lunch inversely correlated with daily intake?

**Diet experiment detection.** There are runs of 6 months of keto, 1 month of low protein, 2 or 3 tries at potato only, weekend 60-hour fasts. Algorithmically detect these regime changes from the food data and test whether any stand out in weight trajectory.

**Calorie misestimation detection.** Create a database of distinct foods eaten (there aren't many) and see if scaling their calories up or down produces a better fit with weight. A frequent restaurant meal that's consistently underestimated would show up. More interestingly, foods that appear super-fattening or near-zero-calorie beyond their label — like pistachios were once predicted to be.

## Ready now (+ steps.csv and sleep.csv)

**Gravitostat.** Compute daily foot-pounds exerted (current weight × steps) and test if this predicts next-day intake or weight change. Firsthand report: "+16 lb weighted vest 3 mile walk kills hunger 12-36 hours after so I can easily cut calories" — but like everything, "stop and it all comes back." Specific hypothesis: the top end of the set point control system discards all excess calories if >x ft-lbs (~8000 steps?) are exerted on legs, and this is the mechanism behind "no weight gain while traveling Europe." The [gravitostat study](https://doi.org/10.1016/j.eclinm.2020.100338) supports this. Impressive if the model can detect the weighted vest period. Concern: seasonal confounding (more steps in summer, different eating too).

**Steps as set point shifter.** Are there conditions where restriction without steps produced no loss, but restriction with >5000 steps did?

**Sunlight-intake correlation.** Use sleep/wake times to count hours of overlap with local sunrise/sunset (Knoxville coordinates, historical weather data) and see if that corresponds to intake. Very little sunlight exposure overall, making any effect easier to detect.

## Needs MFP API enrichment (iron, potassium, vitamin D, saturated fat, added sugars, PUFA/MUFA)

**Omega-6:3 ratio and RMR.** The contamination-adjacent theory. If omega-6 to 3 ratio exceeds some threshold, does apparent RMR drop? Requires fat subtype breakdown, not available from the printable diary. The top ~200 foods cover ~80% of intake and can be enriched via MFP API with serving-size calorie checksums to verify matches.

**Nutrient-specific set point control.** Can a control system predict RMR or next day's intake from micronutrient ratios? Example (pure whimsy): if trailing 30-day omega-6:3 > 10:1 and yesterday's protein > 50g, lower RMR by 500 to store all fat consumed, else raise RMR to maintain set point.

**NOVA classification and binges.** Correlate ultra-processed food consumption with subsequent binges. Requires classifying foods by processing level, feasible from the food names.

## Needs additional data export

**Body temperature during restriction.** Background notes body temp below 97° when consuming below RMR. Withings data in `temperature/temperature.csv` (1,315 readings, 364 days, Dec 2023 – Mar 2026). Directly quantifies thermic adaptation — the body turning down the furnace to defend set point. May be a more useful daily signal than derived RMR.

**More RMR measurements.** Additional Cosmed Fitmate readings exist beyond the 3 currently in `/RMR`. Need to export from device.

**Tirzepatide as set point intervention.** The medication log (Sep 2024 onward) combined with the 13-year pre-intervention baseline is extraordinarily valuable. Tirzepatide silences food noise — the resource-acquisition drive that defends the set point. If hunger (the stick, punishing today's shortfall) and food noise (the carrot, rewarding movement toward set point) are distinct mechanisms, the data should show them separating: pre-tirzepatide restriction had low intake + misery (both firing). On tirzepatide, the same low intake should come without the compounding urge to binge. Test: compare binge frequency at equivalent caloric deficits pre- vs. post-tirzepatide. If binges disappear at the same deficit that previously guaranteed them, food noise was the cause, not hunger. The subjective strength ratings in the medicine log are already tracking this. This is arguably the most important analysis available — the 13 years before are the control, the months on it are the experiment.

**Circadian misalignment.** Consistent 3am-11am sleep with very little sunlight exposure means meals are shifted ~4 hours relative to the solar cortisol rhythm, with essentially no zeitgeber correction. Meal timing relative to circadian phase affects metabolic rate and fat storage independently of calories. This is so constant in the data it's invisible as a variable but extremely unusual compared to study populations. There was an early bedtime experiment — check the sleep data for periods of deviation from the 3am-11am pattern and correlate with weight trajectory. The [housekeepers study](https://doi.org/10.1111/j.1467-9280.2007.01867.x) (tell people their activity is exercise and they lose weight) suggests belief/awareness matters too.

**Hunger taxonomy.** "We have Sapir-Whorf'd ourselves out of understanding the varied nature of hunger, giving it one word while over-granulating 'cravings' by food." At minimum three distinct signals: (1) hunger-the-stick (caloric deficit punishment, same-day), (2) food noise (resource acquisition drive, set-point-distance-proportional), (3) protein-ravenous ("the immediate craving after a protein meal that goes away quickly if you wait"). GLP-1s prove hunger is not purely mechanical stomach-stretch signaling. These may be separately quantifiable from the intake data — days where protein was high but calories low should show different next-day behavior than days where both were low.

**Gut microbiome as missing variable.** Heavy childhood erythromycin (macrolide that reshapes gut microbiome) combined with 35 years of no animal protein recolonization source is a textbook setup for persistent dysbiosis. Not testable from current data, but a microbiome assay (16S or shotgun sequencing) could reveal whether composition is atypical and whether it changes in response to dietary interventions already in the data.

**Low-protein / BCAA restriction.** Tested a low-protein period targeting foods with <2g isoleucine per 2000 kcal. Reported less cold intolerance within days. The intake data should show this period clearly (protein dropping to unusually low levels). Correlate with body temperature and weight trajectory. Related to the protein leverage hypothesis (20% protein diets = 10% protein + 500 cal).

**GLP-1 weekly cycle as mood/intake predictor — CONFIRMED.** "Interesting it's not obvious what my weekly shot day is given the buzzing high at peak and depressive low on day 6." With weekly injection dates from the medicine log and daily intake, test whether intake varies by day-of-injection-cycle. Day 1-2 post-injection should show lowest intake, day 6-7 highest. The subjective strength ratings already capture this partially.

**Result**: The pharmacokinetic blood level model confirms this precisely. By day of injection cycle: day 0 = 1652 cal (blood level 11.7), day 1 = 1741 (11.4), day 2 = 2015 (9.9), day 3 = 2148 (8.6), day 4 = 2139 (7.5), day 5 = 2220 (6.5). The 568 cal swing from peak to trough tracks the PK curve almost exactly — appetite is pharmacokinetically mediated at the daily timescale.

## Modeling approach

Slime Mold recommends [modeling as a control system](https://slimemoldtimemold.com/2022/03/15/control-and-correlation/) where control input has no correlation with modeled variables, but is instead a homeostatic mechanism. A day's macros or total consumed may have no correlation with previous or next days or any other available variable, but instead be correlated with a hidden variable, set point, that meanders slowly up or down. Much like RMR, a secret set point could be determined from the data and binges occur when below the set point.

The [food noise essay](https://lucent.substack.com/p/craving-food-noise) refines this: the set point is defended by two distinct mechanisms. Hunger (the stick) punishes today's shortfall — it's specific in timescale and proportional to deficit. Food noise (the carrot) rewards resource acquisition — it's proportional to distance below set point and compounds over time like an unclaimed package. The control system has two feedback loops, not one. GLP-1 agonists silence food noise but may not affect hunger, which is why the tirzepatide data can distinguish them.

Very open to [use of LLMs for estimating nutrient data](https://chat.openai.com/share/b77dc121-0580-4d66-9e67-131fb3b18a8a). "The following food named x has nutrient profile y. Estimate ingredients and ratio of PUFA to MUFA."

## Previous attempt

The `analysis-bmr-fit/fit.py` script tried to find breakpoints of changing BMR by splitting a time period until the left and right segments' best fit lines had the highest R² value. Did not work well — the breakpoint approach is too rigid for a continuously drifting set point. A state-space model (Kalman filter) is the natural successor.
