# Health background

> **Genetics audit, 2026-07-30:** The genotype calls below may be useful raw annotations, but the former “worst-case profile” and genotype-directed intervention story was not clinically valid. Linked SNPs were double-counted, several mechanisms were overextended, and the IRS1 risk-allele direction was reversed. The interpretation and intervention sections have been corrected; see [`new_findings/05_genetics_and_fatty_acids_audit.md`](new_findings/05_genetics_and_fatty_acids_audit.md).

## Health conditions
* Born October 1982
* Height 5'11.75" without shoes
* Palmoplantar hyperhidrosis. Generally overactive sympathetic nervous system (sinus drainage). 48-hour vasomotor rhinitis from cold exposure.
* Extreme cold intolerance proportional to calories consumed. Body temp gets below 97° when consuming below RMR.
* Recurrent childhood tonsilitis/sinusitis. Heavy use of erythromycin as a child. No non-topical antibiotics used for 6+ years.
* No dental interventions from 1995 to 2015. Enamel/near-dentine cavities filled 2015. Pulp untouched and all 32 teeth intact.
* Elevated LDL (130) which lowered in response to 6 months keto
* COVID 2021-01-05. No illnesses since.
* No other interesting health history, deformities, or abnormalities.
* Overall body shape is slim arms and legs with central obesity, more "hanging off" than the visceral fat look of those barrel-shaped older guys.
* Vegetarian since childhood/early teens. No meat for decades; dairy is the only regular animal product.
* Biracial, half Jordanian (paternal haplogroup E-M78), half European (maternal haplogroup U3a1). Neither parent obese.
* 23andMe and GenesForGood data available for subject and both parents. Imputation/Promethease report available. Raw genetic data is private; analysis below is the public-facing summary.

## Genetic profile (obesity/metabolism SNPs)

Raw 23andMe data for subject and both parents. Genes for Good imputed genotypes fill gaps.
The previous analysis called this a “worst-case common-variant obesity profile.” That is not defensible: linked SNPs within a locus are not independent hits, the list is not an ancestry-calibrated polygenic score, and common variants do not determine treatment.

### FTO — Fat Mass and Obesity-Associated (homozygous risk at all loci)

Both parents are heterozygous carriers at every FTO locus. Subject drew the risk allele from both parents at all five. These five SNPs are in strong linkage disequilibrium (r²>0.85 in Europeans) — they coinherit as a haplotype. ~16% of Europeans carry the full homozygous-risk FTO haplotype.

| SNP | Subject | Father | Mother | EUR hom. risk | Function |
|-----|---------|--------|--------|---------------|----------|
| rs9939609 | **AA** | AT | AT | ~17% | Common obesity association; effect is modest and non-deterministic (about +0.34 BMI units per allele in a large analysis). |
| rs1558902 | **AA** | AT | AT | ~17% | Strongest FTO GWAS signal. +0.39 kg/m² BMI per allele. |
| rs1421085 | **CC** | CT | CT | ~17% | Alters ARID5B/IRX3/IRX5 adipocyte programming in vitro; does not establish clinical low body temperature or cold intolerance. |
| rs17817449 | **GG** | GT | GT | ~17% | Increased ghrelin, reduced satiety signaling. |
| rs8050136 | **AA** | AC | AC | ~17% | Obesity and type 2 diabetes risk. |

### MC4R — Melanocortin 4 Receptor (heterozygous risk, all from father)

These are common nearby GWAS tags in moderate LD, not pathogenic MC4R variants or evidence of MC4R deficiency. Risk allele frequency ~23-28% in Europeans.

| SNP | Subject | Father | Mother | EUR carrier freq | Function |
|-----|---------|--------|--------|------------------|----------|
| rs17782313 | **TC** | CT | TT | ~36% | +0.22 kg/m² BMI per allele. Increased hunger drive. |
| rs571312 | **AC** | AC | CC | ~36% | Strong obesity association, central appetite regulation. |
| rs12970134 | **AG** | AG | GG | ~40% | Waist circumference and insulin resistance. |

### FADS1/FADS2 — Fatty Acid Desaturase (high activity — efficient omega-6 converter)

These variants are associated with fatty-acid profile and desaturase measures. They do not establish that dietary linoleic acid creates systemic inflammation or obesity; human trials generally do not show higher inflammatory markers from increased linoleic acid.

| SNP | Subject | Father | Mother | EUR hom. freq | Function |
|-----|---------|--------|--------|---------------|----------|
| rs174546 | **CC** | CC | CT | ~32% | High delta-5 desaturase activity. Both C alleles from father. |
| rs174547 | **TT** | TT | CT | ~32% | Confirms high FADS1 activity. (Genes for Good imputed.) |
| rs1535 | **AA** | AG | AG | | |
| rs174575 | **CC** | CG | CC | | |

### GLP1R — GLP-1 Receptor

| SNP | Subject | Father | Mother | EUR freq | Function |
|-----|---------|--------|--------|----------|----------|
| rs6923761 | **GG** | GG | AG | ~58% | No consistent, clinically actionable prediction of GLP-1 drug response; actual treatment response is more informative. |

### UCP — Uncoupling Proteins / Thermogenesis

The heterozygous UCP2 calls do not provide a clear genetic diagnosis of low RMR, reduced whole-body thermogenesis, or cold intolerance. UCP2 is not interchangeable with UCP1, and association studies are inconsistent.

| SNP | Gene | Subject | Father | Mother | EUR carrier freq | Function |
|-----|------|---------|--------|--------|------------------|----------|
| rs1800592 | UCP1 | not typed | TT | CT | | Father has no thermogenesis variant. |
| rs659366 | UCP2 | **CT** | CC | CT | ~44% | Common heterozygous call; no established individual-level RMR or obesity prescription. |
| rs660339 | UCP2 | **AG** | GG | AG | ~44% | Common heterozygous call; does not establish reduced energy expenditure or cold intolerance. |
| rs1800849 | UCP3 | GG | GG | GG | | Normal. |

### Adrenergic Receptors — Lipolysis / Fat Breakdown

~12-15% of Europeans are homozygous risk at both ADRB2 loci (partially linked).

| SNP | Gene | Subject | Father | Mother | EUR hom. freq | Function |
|-----|------|---------|--------|--------|---------------|----------|
| rs1042713 | ADRB2 | **GG** | AG | GG | ~37% | Gly/Gly. Altered lipolysis response to catecholamines. |
| rs1042714 | ADRB2 | **GG** | CG | GG | ~36% | Glu/Glu. Associated with increased fat mass. Both risk alleles from mother. |
| rs4994 | ADRB3 | AA | AA | AA | | Normal. No additional RMR reduction from this locus. |

### Other Relevant SNPs

| SNP | Gene | Subject | EUR freq | Risk? | Function |
|-----|------|---------|----------|-------|----------|
| rs1801282 | PPARG | **CC** | ~82% | Pro/Pro — lacks protective Ala allele | Common risk genotype for insulin resistance and fat storage. |
| rs5082 | APOA2 | AG | ~39% | Heterozygous carrier | Foundational saturated-fat interaction was recessive in CC homozygotes; no validated heterozygote sensitivity. |
| rs6265 | BDNF | CC | ~80% | Normal | No additional appetite dysregulation. |
| rs1137101 | LEPR | AA | ~35% | Normal | Leptin receptor signaling intact. |
| rs4988235 | LCT | AG | ~41% | Lactase persistent | Can digest dairy — relevant since dairy is the only animal product consumed. |
| rs1801260 | CLOCK | AA | ~65% | Normal | No circadian gene penalty (though behavior is shifted 4 hours regardless). |
| rs5443 | GNB3 | CT | ~37% | Heterozygous carrier | One copy of enhanced signal transduction variant, associated with obesity and hypertension. |
| rs2943641 | IRS1 | CC | ~56% | Association-level | C, not T, is the reported insulin-resistance/T2D-risk allele; the effect is small and does not prescribe magnesium timing. |

### Genetic summary

The FTO haplotype is associated with a modestly higher average BMI/appetite risk, but all five listed FTO calls mostly represent one linked signal. The common MC4R-region markers are not receptor deficiency. FADS variants can affect fatty-acid profiles without proving inflammation. UCP2, ADRB2, APOA2, and GLP1R results do not justify personalized treatment prescriptions.

Normal single SNPs also do not prove that an entire pathway is intact or has “no target.” The clinical phenotype and response to treatment matter more than these annotations.

### Intervention implications after audit

- **Tirzepatide:** supported by the observed clinical response and external trial evidence, not by rs6923761. Dosing remains clinical, not genotype-based.
- **EPA/DHA:** reasonable to discuss for nutrient adequacy in a long-term vegetarian with low EPA/DHA status. It is not a FADS-specific weight-loss treatment.
- **Cold exposure:** not genetically indicated and a poor fit for severe cold intolerance and prolonged vasomotor rhinitis. Ordinary clinical causes of cold intolerance deserve consideration before a cold experiment.
- **Taurine:** human weight-loss evidence is insufficient to justify escalation for this genotype.
- **Magnesium:** there is no genotype-based reason to time it to the highest-carbohydrate meal. If “500 mg” is elemental magnesium, verify the label and indication with a clinician because it exceeds the U.S. supplemental upper limit for adults.
- **Saturated fat:** investigate for conventional lipid outcomes if relevant, not because APOA2 AG establishes special sensitivity.
- **Further genetic work:** first verify strand, genome build, and imputation quality. A validated ancestry-aware genome-wide score would be required to quantify polygenic risk, but would still not choose routine obesity treatment.

## Weighing protocol
* Always upon waking
* Gaps of months without weigh-ins (especially 2013-09 to 2015-02, 2017-10 to 2019-10) are avoidance — weight was going up. Readings on either side of gaps are valid. Intake is fully logged throughout.

## Substance use
* No exposure to deodorant or lotions (Unusually low BO. Socks or underwear worn for days are indistinguishable by smell from freshly laundered) 
* Head & Shoulders shampoo (zinc pyrithione)
* No prescription drug use until tirzepatide late 2024
* Prescription acne topical cream history (basic stuff, no Accutane)
* Daily 10mg cetirizine (Zyrtec) mornings through ~April 2025, then switched to 2.5mg levocetirizine (Xyzal) nightly (half a 5mg pill). Net: 5mg active enantiomer → 2.5mg, morning → night.
* Nightly Rogaine 5% foam front temples
* 500mg magnesium, 500mg taurine mornings 2022-01-01 to present
* Creatine monohydrate history:
  2018-10-29 to 2021-04-28: 2g/day
  2025-05-09 to present: 5g/day
* Invisalign multilayer aromatic thermoplastic polyurethane retainer use full time 2015 then nightly 2016 to present

## Locations lived
Extensive records of location history for the duration of the intake data. Individual vacations/restaurants and menu items can be reverse engineered from this.

### UT Knoxville adjacent apartment building
* 2002 to 2015-03-01
* Electric range
* Window unit AC/heat
* 3rd floor apartment
* Brita filtered water

### UT Knoxville adjacent duplex
* 2015-03-01 to 2020-03-30
* Electric range
* Central AC/heat
* Top half of duplex
* Brita filtered water

### Knoxville, Sequoyah Hills home
* 2020-02-01 to present
* Gas range
* Whole house UV and MERV11+ filtration
* Whole house carbon filter + RO drinking water
