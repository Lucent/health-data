# Obesity hypothesis evidence map

Date: 2026-07-30

The useful question is not “which camp is right?” It is which mechanisms have strong causal evidence, which could explain this person's history, and which predictions the repository can distinguish without fooling itself.

## High-confidence foundations

### Energy intake and expenditure still determine weight change

This is an accounting identity, not a full behavioral explanation. Long-term weight change reflects energy imbalance, but appetite, food environment, medication, adaptive physiology, activity, and measurement error determine both sides. Self-reported intake is often biased, making this repository's uninterrupted daily food ledger unusually valuable even though portion accuracy and absolute error are not independently validated.

Useful sources:

- [Lichtman et al., discrepancy between self-reported intake/exercise and obesity](https://pubmed.ncbi.nlm.nih.gov/1454084/)
- [Schoeller, doubly labeled water and self-report validity](https://pubmed.ncbi.nlm.nih.gov/2200586/)
- [Hall and Chow, quantitative physiology of weight change](https://pubmed.ncbi.nlm.nih.gov/21562087/)

Repo prediction: locked, forward models using logged intake should predict future weight better than calorie-blind baselines, with error varying across measurement and food-source regimes.

### Weight loss triggers homeostatic opposition

After weight loss, appetite often increases and expenditure can decrease beyond simple body-size changes. This supports a defended-weight or settling-point family, but not one exact SmoothLatch architecture.

- [Polidori et al., feedback control of energy intake](https://pubmed.ncbi.nlm.nih.gov/27804272/)
- [Sumithran et al., persistent hormonal adaptation](https://pubmed.ncbi.nlm.nih.gov/22029981/)
- [Rosenbaum and Leibel, adaptive thermogenesis](https://pubmed.ncbi.nlm.nih.gov/22237063/)
- [Systematic review: adaptive thermogenesis is often small and attenuates after stabilization](https://pubmed.ncbi.nlm.nih.gov/32844188/)
- [Review of multiple body-weight control architectures](https://pubmed.ncbi.nlm.nih.gov/37661735/)

Repo prediction: hunger/intake should rise following sustained loss even when weekday, medication, food availability, and recent intake are controlled. Several model families should be compared on a future time block.

### Food processing can alter spontaneous intake

In a controlled inpatient crossover trial, an ultra-processed diet increased ad-libitum intake by about 500 kcal/day and produced short-term weight gain despite diets being presented as matched on several nutrients.

- [Hall et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31105044/)

Repo prediction: prospectively classified meal templates should predict later same-day and next-day intake after matching on calories already consumed, protein, energy density, medication, weekday, travel, and availability. “Domino's day” is a context bundle, not an ingredient diagnosis.

### Modern anti-obesity pharmacotherapy is highly effective but usually chronic

Tirzepatide and other incretin therapies have the strongest treatment signal relevant to this dataset. Withdrawal commonly produces regain, so the burden of proof is on any personalized claim that a holding period makes cessation safe.

- [2026 BMJ pharmacotherapy network meta-analysis](https://pubmed.ncbi.nlm.nih.gov/42419792/)
- [SURMOUNT-4 withdrawal trial](https://pubmed.ncbi.nlm.nih.gov/38078870/)
- [STEP 1 withdrawal extension](https://pubmed.ncbi.nlm.nih.gov/35441470/)
- [WHO guideline announcement](https://www.who.int/news/item/01-12-2025-who-issues-global-guideline-on-the-use-of-GLP-1-medicines-in-treating-obesity)

Repo prediction: treatment-era reductions should be visible in meal size, daily intake distribution, binge onset, and binge persistence. Exact dose-response and cycle effects require variation not locked to weekday.

## Moderate evidence and likely support levers

### Protein and resistance training

Higher-protein diets can modestly improve satiety and lean-mass retention for some people, while resistance training helps preserve lean mass during weight loss. These are support strategies, not magic independent fat-loss mechanisms.

- [Higher-protein diet meta-analysis](https://pubmed.ncbi.nlm.nih.gov/23221572/)
- [Protein and appetite/weight-control review](https://pubmed.ncbi.nlm.nih.gov/34120735/)
- [Resistance exercise during dietary weight loss](https://pubmed.ncbi.nlm.nih.gov/35191588/)
- [Joint nutrition advisory for GLP-1 therapy](https://pubmed.ncbi.nlm.nih.gov/40452753/)

Repo prediction: a prospectively chosen protein target and measured lifting volume should improve hunger, adherence, strength, or same-device lean-mass outcomes relative to matched prior periods. Do not derive a threshold by repeatedly searching cut points.

### Exercise

Exercise has major health benefits and can aid maintenance, but average weight loss from exercise alone is often modest because of compensation. The repository cannot reliably distinguish a precise RMR bonus from 25 clustered calorimetry tests.

- [ACSM consensus statement](https://pubmed.ncbi.nlm.nih.gov/39277776/)

Repo prediction: use measured training dose and cardiovascular fitness/strength outcomes. Test weight effects prospectively and avoid treating device-estimated calorie burn as ground truth.

### Sleep

Sleep restriction can increase intake, and extending sleep reduced intake in a randomized trial of habitual short sleepers. This subject's median duration is relatively long, so duration may not be a large lever even though the current repo “null” is methodologically weak.

- [Sleep extension randomized trial](https://pubmed.ncbi.nlm.nih.gov/35129580/)

Repo prediction: first clean primary sleeps and restrict historical inference to high-coverage years. Test schedule regularity and short-sleep deviations within person, not a broad claim that sleep never matters.

## Plausible but unresolved here

### Low-carbohydrate diets and the carbohydrate-insulin model

Low-carbohydrate diets may help some people through satiety, food preference, glycemic control, or adherence. Strong claims that carbohydrate uniquely drives fat gain independent of calories have not held up well in controlled feeding studies.

- [Plant-based low-fat versus animal-based ketogenic ad-libitum trial](https://pubmed.ncbi.nlm.nih.gov/33479499/)
- [Isocaloric ketogenic diet metabolic-ward study](https://pubmed.ncbi.nlm.nih.gov/27385608/)
- [DIETFITS randomized trial](https://pubmed.ncbi.nlm.nih.gov/29466592/)

Repo prediction: compare pre-specified low-carb templates on later intake and persistence, not short-term water-weight loss. Match medication era and food monotony.

### Menu monotony, sensory-specific satiety, and decision friction

The potato and repetitive-meal periods make this especially plausible. Monotony may reduce choice, reward, and decision load without potatoes having a unique metabolic property.

Repo prediction: compute pre-specified menu entropy and repetition metrics from food strings, then test future intake on a locked period. Separate within-day reverse causality from next-day effects.

### Travel and environmental disruption

Travel bundles availability, meal sourcing, social context, sleep timing, movement, and measurement changes. It is likely actionable even if no single mechanism is identified.

Repo prediction: event-study each trip with comparable pre/post windows; model missing weigh-ins and restaurant/order-in templates explicitly.

### Circadian timing

Meal timing may matter through behavior or circadian physiology, but historical meal timestamps are incomplete and legacy “front-loading” tests are vulnerable to reverse causality.

Repo prediction: collect timestamps prospectively and compare time-consistent periods without changing total calories or medication.

### Microbiome

Biologically plausible and difficult to infer from this repository. Without stool sampling or a randomized dietary perturbation, microbiome explanations are post hoc stories.

Repo prediction: none from existing data that can identify a microbiome mechanism.

## Weak, overstated, or poor-priority theories

### Seed oils as a special obesity toxin

Seed-oil-rich ultra-processed food can promote overconsumption, but that does not establish linoleic acid as the causal agent. Randomized-trial reviews do not show that increasing linoleic acid raises inflammatory markers, and one controlled overfeeding trial found polyunsaturated fat produced less liver and visceral fat than saturated fat.

- [15-trial review of linoleic acid and arachidonic acid](https://pubmed.ncbi.nlm.nih.gov/22889633/)
- [30-trial review of linoleic acid and inflammation](https://pubmed.ncbi.nlm.nih.gov/28752873/)
- [LIPOGAIN trial](https://pubmed.ncbi.nlm.nih.gov/24550191/)

Repo limitation: polyunsaturated-fat enrichment covers only about 42% of logged calories and does not identify cooking-oil ingredients reliably. Testing “seed oils” now would mainly rediscover source regime and restaurant/UPF patterns.

### Common obesity SNPs as personalized treatment instructions

FTO associations are real but small; linked SNPs do not add independent hits. Common MC4R-region tags are not MC4R deficiency, and UCP2/ADRB2/APOA2/GLP1R variants in the repo do not validate personalized cold, taurine, magnesium, saturated-fat, or tirzepatide prescriptions.

- [FTO mechanism](https://pubmed.ncbi.nlm.nih.gov/26287746/)
- [FTO does not materially alter weight-loss intervention response](https://pubmed.ncbi.nlm.nih.gov/27650503/)
- [MC4R-region common variant](https://pubmed.ncbi.nlm.nih.gov/18454148/)
- [GLP1R pharmacogenetics systematic review](https://pubmed.ncbi.nlm.nih.gov/36999540/)

Repo prediction: genotype is fixed in one person and cannot explain time variation. Analyze the phenotype directly.

### Cold exposure and “brown-fat rescue”

Cold can acutely activate thermogenesis, but sustained clinically meaningful weight loss is unproven and appetite compensation is plausible. Severe cold intolerance makes this a particularly poor experiment without clinical evaluation.

- [Cold-induced thermogenesis review](https://pubmed.ncbi.nlm.nih.gov/36568070/)

### Peripheral temperature and baths as TDEE meters

Wearable/peripheral temperature, bath events, and core thermogenesis are not interchangeable. Measurement clustering, sparse date coverage, seasonal confounding, and lack of raw bath-meter coverage make current claims weak.

## Best discriminating experiments

1. **Prospective high-risk context substitution:** predefine travel/order-in occasions and a satisfying replacement template; compare total intake, hunger, and binge persistence.
2. **Menu-entropy natural experiment:** keep calories/protein/medication stable while alternating a small repeatable menu with a varied menu in blocks.
3. **Food-processing crossover:** use two calorie-unrestricted but palatability-matched patterns differing in processing and energy density, with clinician/dietitian review.
4. **Lean-mass protection block:** record lifting exercises, sets, reps, load, proximity to failure, and protein; use strength plus same-device composition as outcomes.
5. **Medication observation, not improvisation:** log prescribed dose, injection time, hunger, GI effects, and intake prospectively; do not alter frequency to identify a model.
6. **Same-lab fatty-acid follow-up:** after at least 16 weeks of a stable verified EPA/DHA exposure, repeat the same lab. This tests uptake, not weight loss.
