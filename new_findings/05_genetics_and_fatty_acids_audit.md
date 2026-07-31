# Genetics and fatty-acid audit

Date: 2026-07-30
Status: evidence correction; genotype calls should be verified from raw build/strand before further analysis

## Bottom line

The genetic data does not reveal a special missing key to weight loss. It contains several real common-variant associations, especially an FTO risk haplotype, but the repo previously multiplied linked signals and converted population associations into individual mechanistic diagnoses.

What survives:

- FTO risk variants modestly raise average BMI/appetite risk.
- FADS variants are associated with fatty-acid profile/desaturase measures.
- The actual tirzepatide response is strong.
- EPA/DHA status can be monitored directly.

What does not survive:

- “Worst-case common-variant profile”
- MC4R deficiency
- genetically proven low RMR/cold intolerance
- GLP1R-predicted tirzepatide response
- genotype-directed cold exposure, taurine escalation, or magnesium timing
- APOA2 heterozygote saturated-fat sensitivity
- seed-oil/linoleic-acid inflammation as an established cause

The clearest factual error was IRS1 rs2943641: the repo said T was the insulin-resistance risk allele. Published evidence reports **C**, which the subject carries. The effect is common and small; it still does not prescribe magnesium timing.

## FTO: one modest locus, not five independent strikes

The five listed FTO SNPs are in strong linkage disequilibrium. They should be interpreted mainly as one inherited haplotype.

- rs1421085 C has a supported molecular effect on ARID5B/IRX3/IRX5 adipocyte programming, but that in vitro mechanism does not establish clinical low core temperature or cold intolerance: [Claussnitzer et al.](https://pubmed.ncbi.nlm.nih.gov/26287746/).
- Large analyses put the common per-allele BMI effect in the range of a few tenths of a BMI unit, explaining very little individual variance—not a deterministic obesity syndrome.
- In a 9,563-person intervention meta-analysis, FTO genotype did not materially change weight-loss response: [Livingstone et al.](https://pubmed.ncbi.nlm.nih.gov/27650503/).

Practical conclusion: ordinary effective obesity treatment still works. The genotype does not select low-carb, cold exposure, or a drug dose.

## MC4R-region tags are not MC4R deficiency

rs17782313 is roughly 188 kb downstream of MC4R and has a small population association with BMI. The other listed markers are also common GWAS tags. They are not equivalent to pathogenic MC4R coding variants or monogenic hyperphagia: [primary GWAS](https://pubmed.ncbi.nlm.nih.gov/18454148/).

This matters because “MC4R risk” can sound like a receptor disorder. Nothing in these calls justifies setmelanotide or a monogenic-obesity treatment frame.

## FADS and omega-6: biomarker effect, not inflammatory verdict

rs174546 C and rs174547 T are associated with desaturase measures and fatty-acid profiles. That does not establish that this person converts dietary linoleic acid into clinically harmful inflammation.

Human evidence cuts against the strong repo claim:

- Increasing linoleic acid generally does not substantially raise tissue arachidonic acid: [systematic review](https://pubmed.ncbi.nlm.nih.gov/21663641/).
- Randomized trials do not show higher inflammatory markers from increased linoleic acid: [30-trial review](https://pubmed.ncbi.nlm.nih.gov/28752873/).
- Higher circulating linoleic acid has been associated with lower diabetes risk, while arachidonic acid was broadly neutral: [pooled analysis](https://pubmed.ncbi.nlm.nih.gov/29032079/).

The repository is also not ready for a seed-oil analysis: polyunsaturated-fat enrichment covers only 46.9% of food rows and about 42% of logged calories, with no reliable ingredient-level oil decomposition.

## GLP1R, UCP2, ADRB2, and APOA2

- **GLP1R rs6923761:** small studies are inconsistent and a systematic review found no consistent treatment-response association. No FDA/CPIC/DPWG genotype dosing exists: [systematic review](https://pubmed.ncbi.nlm.nih.gov/36999540/).
- **UCP2:** UCP2 is not UCP1 and is not established as a simple whole-body thermogenic proton-leak switch. Common heterozygous calls do not diagnose low RMR: [2024 review](https://pubmed.ncbi.nlm.nih.gov/38927514/).
- **ADRB2:** association meta-analyses are inconsistent or small. “Impaired fat breakdown” is too strong for an individual.
- **ADRB3:** a normal receptor does not mean beta-3 agonists have “no target”; the receptor is the target. Experimental BAT activation does not make these drugs established obesity treatment.
- **APOA2 rs5082:** the replicated saturated-fat interaction was recessive in CC homozygotes. The subject's reverse-strand AG corresponds to CT and was grouped with TT in the foundational model: [primary evidence](https://pubmed.ncbi.nlm.nih.gov/19901143/).

## Intervention implications

### Keep what is working for clinical reasons

Tirzepatide response is directly observed and consistent with randomized evidence. rs6923761 adds no useful decision information.

### EPA/DHA is a nutrient-status question, not a weight-loss theory

For a long-term vegetarian with low EPA/DHA status, algal EPA/DHA may be reasonable after reviewing dose, diet, bleeding/medication context, and clinical goals. Omega-3 supplements have not reliably produced meaningful weight loss.

The best test is simple: after at least 16 weeks on a stable verified dose, repeat the same OmegaQuant assay under similar conditions. A higher EPA/DHA/Omega-3 Index and lower AA:EPA would show uptake. It would not prove a weight-loss mechanism.

### Do not use cold exposure as a genotype experiment

Cold can acutely raise expenditure, but sustained weight-loss evidence is weak and appetite compensation can occur. Severe cold intolerance and 48-hour vasomotor rhinitis make the risk/benefit especially poor. Evaluate ordinary causes of cold intolerance clinically if warranted rather than treating a common SNP story.

### Do not escalate taurine or retime magnesium based on these SNPs

The small human taurine trial did not demonstrate a between-group weight advantage: [trial](https://pubmed.ncbi.nlm.nih.gov/24065043/). There is no validated PPARG/IRS1 reason to take magnesium with the highest-carbohydrate meal.

If the listed “500 mg magnesium” means 500 mg **elemental** magnesium from supplements, it exceeds the U.S. adult supplemental upper limit of 350 mg/day. The label may instead list 500 mg of a compound with much less elemental magnesium. Check the Supplement Facts panel, indication, kidney function, and GI effects with a clinician: [NIH Office of Dietary Supplements](https://ods.od.nih.gov/factsheets/Magnesium-HealthProfessional/).

## Fatty-acid panels: what the three draws say

The repo has 70 analyte rows but only three draw dates:

- 2019 Cleveland HeartLab
- 2024 OmegaQuant
- 2025 OmegaQuant

The two same-lab panels show Omega-3 Index 4.93% to 6.06%, EPA 0.54% to 0.45%, DHA 2.63% to 3.71%, LA 25.03% to 22.93%, and AA 7.52% to 8.93%. The 2019 panel is less comparable because lab, matrix, and method may differ.

The commonly cited Omega-3 Index, AA:EPA, and n6:n3 ranges are proposed risk or lab ranges, not established obesity-treatment targets. See the [NIH omega-3 review](https://ods.od.nih.gov/factsheets/Omega3FattyAcids-HealthProfessional/).

The 2025 draw occurred only 20 days after supplementation began. That is too early for full red-cell steady state, but not too early for any incorporation. Therefore the year-over-year rise cannot be assigned to diet alone.

## Useful follow-up, if desired

1. Verify strand, build, and imputation quality for every publicly discussed SNP.
2. Do not compute a hand-selected “combined prevalence.” If genetic risk quantification matters, use a validated ancestry-aware genome-wide score with appropriate counseling; it still will not dictate routine treatment.
3. If low thermogenesis is a genuine clinical question, use repeated standardized indirect calorimetry adjusted for body composition. Peripheral temperature is not a BAT measurement.
4. If insulin resistance is a question, test the phenotype with conventional clinical markers or a clinician-directed mixed-meal/CGM evaluation. Do not infer it from one SNP.
5. Repeat fatty acids only to answer a pre-specified nutrient-status question, with the same lab and stable exposure.
