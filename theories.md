# Theories to test

* Derive a formula to smooth daily perturbations in weight from glycogen-water binding (water weight). A few days below RMR expell it and a few days above recapture it. Useful for more precise analysis of other points below.
* Make trendline fit for intake vs. weight that has slow-changing RMR variable to force a fit. Find what high/low RMR correlates with, like a macronutrient. Graph RMR over time.
* Make trendline fit for intake vs. weight but hold RMR constant and change the calories per gram of carbs, protein, fat (4, 4, 9) to force a fit (5, 3, 7?).
* Look for set point abruptly increasing or decreasing and examine the few days or weeks when it appeared to happen. Is it possible only multi-day binges can shift set-point up, like 3 days of 3000 calories? Do days of 3000 calories spread out cause no shift in set point?
* Possible that set point only moves up in the presence of protein to capture high-effort hunted meals? Studies showing 20% protein = 10% protein + 500 calories.
* Look for long runs (months) of stable weight (±3 lbs) and see if intake does not support no change, which would be evidence for a set point that can only drift under certain circumstances.
* Are there conditions that must be met to shift set point down? Days of restriction without steps resulted in no loss but with > 5000 steps did?
* Test gravitostat theory by computing foot-pounds exerted per day (current weight * steps) and seeing if that leads to less intake the following day. Impressive if this model can detect anomalies like days I wore a weighted vest.
* Use sleep/wake times to count hours of overlap with local sunrise/sunset and see if that corresponds to intake.
* Create a database of distinct foods eaten (there aren't many) and see if scaling their calories up or down produces a better fit, suggesting a misestimate of a frequent restaurant meal, or more interestingly, super fattening or near-zero calorie type foods like pistachios were predicted to be.
* Correlate a nutrient, food, or food classification (NOVA) with subsequent binges.
* Track an independent "willpower" variable that when exhausted predicts a binge.
* Does eating over 1000 calories for breakfast and lunch combined predict a lower total for the day?
* There are runs of 6 months of keto, 1 month of low protein, 2 or 3 tries at potato only, weekend 60-hour fasts. Do these stand out in any way?
* Fiber consumed in breakfast+lunch inversely correlated with daily intake?
* What best predicts whether a day will become a binge (> 2800)? Sugar intake the previous day? Being below current set point?
* What conditions push the set point down? Are there magic months where I lose 3 lbs and don't regain it for a long time?
* Can a control system predict RMR/next day's intake? If omega-6 to 3 ratio exceeds 10:1 and protein exceeds 50g, lower RMR to store all fat consumed, else raise RMR to maintain current set point.

Very open to [use of LLMs for the purpose of estimating other nutrient data](https://chat.openai.com/share/b77dc121-0580-4d66-9e67-131fb3b18a8a). "The following food named x has nutrient profile y. Estimate ingredients and ratio of PUFA to MUFA." OpenAI API key available for these requests.

Slime Mold recommends [modeling as a control system](https://slimemoldtimemold.com/2022/03/15/control-and-correlation/) where control input has no correlation with modeled variables, but is instead a homeostatic mechanism. A day's macros or total consumed may have no correlation with previous or next days or any other available variable, but instead be correlated with a hidden variable, set point, that meanders slowly up or down. Much like RMR, a secret set point could be determined from the data and binges occur when I am below the set point. 

## [/analysis-bmr-fit](/analysis-bmr-fit/fit.py)
* Tries to find breakpoints of changing BMR by splitting a time period until the left and right segments' best fit lines have the highest R^2 value.
* Does not appear to be a great approach.
* Next step: revert script to NUM_LINES=1 and instead add modifier columns and target an increase of a single trendline's R^2 by:
  * smoothing weight (control system that sheds water weight or increases it from previous days' rolling consumption)
  * giving calorie values to daily steps
