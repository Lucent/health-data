# References

How the findings in this dataset relate to published literature. Organized by finding, not by paper.

The broad architecture of these findings is not alien to the literature. A defended-adiposity framework, coordinated appetite and expenditure compensation, and regain pressure after weight loss are all well established. What this dataset adds is a more specific operationalization: a moving defended fat-mass signal, a fitted adaptation timescale, a binge-frequency readout of distance below that signal, and a possible separation between appetite and expenditure clocks. The strongest overall claim is not that the existing literature was wrong, but that it had the right qualitative structure and that this dataset may sharpen, parameterize, and in a few places extend it.

## Set point, defended adiposity, and lipostat vs gravitostat

**AG is more consistent with a lipostatic than a gravitostatic account.** The classical lipostatic view proposes that adipose tissue generates a signal that the brain uses to regulate intake and expenditure around a defended level of adiposity. Leptin provided the strongest molecular support for that framework. In this dataset, fat mass predicts binge rate better than total weight or scale weight across all tested half-lives, which is the pattern one would expect if the regulated variable is adiposity rather than gross body mass. By contrast, finding N does not support a gravitostat-like signal in these data (foot-pounds → next-day intake: r=+0.05, wrong direction).

This does not prove the full lipostat theory in the strong sense, nor does it settle the broader debate experimentally. But it places AG much closer to the defended-adiposity / leptin tradition than to a body-load or weight-per-se mechanism.

- Kennedy GC. The role of depot fat in the hypothalamic control of food intake in the rat. *Proc R Soc Lond B* 140:578-592, 1953.
- Zhang Y, Proenca R, Maffei M, et al. Positional cloning of the mouse obese gene and its human homologue. *Nature* 372:425-432, 1994. [doi:10.1038/372425a0](https://doi.org/10.1038/372425a0)
- Jansson JO, Palsdottir V, Hägg DA, et al. Body weight homeostat that regulates fat mass independently of leptin in rats and mice. *PNAS* 115:427-432, 2018. [doi:10.1073/pnas.1800033115](https://doi.org/10.1073/pnas.1800033115) (Gravitostat hypothesis. The weighted-vest [RCT](https://doi.org/10.1016/j.eclinm.2020.100338) supports it in humans; finding N does not support it in these data.)

**AG supports a moving defended level rather than a fixed set point or dual intervention points.** The literature has long entertained the possibility that the defended level can shift with chronic state, but usually in conceptual terms rather than as an explicitly fitted dynamic variable. AG goes further by showing that a moving exponential average of fat mass predicts binge frequency much better than any fixed defended level. Speakman's framework proposes weak regulation within a zone of indifference and stronger responses near upper and lower intervention points. In this dataset, binge probability changes continuously with distance below the reconstructed set point, with no obvious flat region — arguing against a large free-drift zone in this particular outcome, though not ruling out all intervention-point formulations.

- Speakman JR, Levitsky DA, Allison DB, et al. Set points, settling points and some alternative models: theoretical options to understand how genes and environments combine to regulate body adiposity. *Dis Model Mech* 4:733-745, 2011. [doi:10.1242/dmm.008698](https://doi.org/10.1242/dmm.008698)

## Set point adaptation timescale

**AG provides a fitted timescale for a moving defended-weight signal inferred from binge-frequency data.** The broader literature includes dynamic and first-order energy-balance models with explicit time constants, so it would be too strong to say no half-life-like quantities have ever been estimated. What appears more specific here is the use of a moving defended fat-mass signal, inferred from binge behavior, with an empirically fitted half-life.

The estimate is broad rather than razor-sharp: the plateau across nearby half-lives suggests the important result is the order of magnitude (weeks) more than the exact point estimate. AG supports an appetite-related defended-weight timescale on the order of several weeks, centered near 50 days in this dataset (CI [20, 195]).

The Biggest Loser literature is relevant mainly as a contrast. That work documents persistent metabolic adaptation in reduced-weight individuals, but it does not directly measure set-point location. It therefore cannot by itself distinguish a slowly adapting set point from a subject who remains chronically below a defended level. AI's expenditure arm HL ≤10 days suggests the "persistent" adaptation may be the latter: the subjects' set points may have moved, but they kept chasing them down.

- Fothergill E, Guo J, Howard L, et al. Persistent metabolic adaptation 6 years after "The Biggest Loser" competition. *Obesity* 24:1612-1619, 2016. [doi:10.1002/oby.21538](https://doi.org/10.1002/oby.21538)
- Hall KD, Sacks G, Chandramohan D, et al. Quantification of the effect of energy imbalance on bodyweight. *Lancet* 378:826-837, 2011. [doi:10.1016/S0140-6736(11)60812-X](https://doi.org/10.1016/S0140-6736(11)60812-X)

## Binge frequency and defended weight

**AG extends the weight-suppression literature by replacing a crude historical proxy with a moving model-derived signal.** Lowe and colleagues showed that weight suppression (the gap between current and highest previous weight) predicts binge eating in bulimia nervosa. That finding is conceptually close to the idea that being below a defended level increases pressure toward overeating. AG sharpens that framework by replacing "highest-ever weight" with a moving defended fat-mass estimate and by showing a graded relation (r=-0.62, continuous sigmoid) between distance below that estimate and binge frequency. The underlying idea is similar; the level of mechanistic precision is not.

More broadly, the reduced-weight-state literature is strongly consistent with the direction of AG. Maintaining reduced weight is associated with increased hunger, altered neural responses to food cues, and reduced expenditure; several of these changes are at least partly leptin-reversible.

- Lowe MR, Thomas JG, Safer DL, Butryn ML. The relationship of weight suppression and dietary restraint to binge eating in bulimia nervosa. *Int J Eat Disord* 40:640-644, 2007. [doi:10.1002/eat.20405](https://doi.org/10.1002/eat.20405)
- Rosenbaum M, Goldsmith R, Bloomfield D, et al. Low-dose leptin reverses skeletal muscle, autonomic, and neuroendocrine adaptations to maintenance of reduced weight. *J Clin Invest* 115:3579-3586, 2005. [doi:10.1172/JCI25977](https://doi.org/10.1172/JCI25977)
- Rosenbaum M, Sy M, Pavlovich K, Leibel RL, Hirsch J. Leptin reverses weight loss-induced changes in regional neural activity responses to visual food stimuli. *J Clin Invest* 118:2583-2591, 2008. [doi:10.1172/JCI35055](https://doi.org/10.1172/JCI35055)

## Expenditure hysteresis

**K replicates the direction and asymmetry of the controlled literature under noisier but ecologically richer conditions.** Leibel, Rosenbaum, and Hirsch showed in a controlled ward setting (n=18) that maintaining body weight below usual levels is associated with lower expenditure, and maintaining it above usual levels is associated with higher expenditure. K finds +179 cal/day (falling vs rising at matched FM, 1,482 pairs, p<10^-22) across 15 years of free-living data. The asymmetry is the same: falling-phase TDEE is elevated, rising-phase is not. K additionally shows the effect grows with fat mass (+112 at FM 25-45, +332 at FM 65-85). The strongest claim is not that K proves a new theory of expenditure, but that the same branch-dependent compensation can still be detected outside a ward.

- Leibel RL, Rosenbaum M, Hirsch J. Changes in energy expenditure resulting from altered body weight. *N Engl J Med* 332:621-628, 1995. [doi:10.1056/NEJM199503093321001](https://doi.org/10.1056/NEJM199503093321001)

(The Rosenbaum/Leibel leptin reversal papers cited above under binge frequency are from the same group extending this same ward work into mechanistic territory.)

## Dual-arm timescale separation

**AI is best presented as a novel inference from this dataset.** The body-weight literature already supports the existence of both appetite-side and expenditure-side compensation during reduced-weight maintenance. What was not identified in the literature is prior work that explicitly fits separate adaptation clocks for the two arms and concludes that the expenditure arm adapts materially faster than the appetite arm (expenditure HL ≤10 days, appetite HL ~50 days). The existence of two arms is established. The precise timescale separation appears to be the new part.

## GLP-1 agonists and the set point

**AG and AJ are more consistent with suppression than durable resetting.** When GLP-1 receptor agonists are stopped, a substantial fraction of lost weight is regained, often fairly quickly. That pattern fits better with continued underlying biological pressure that was being pharmacologically suppressed than with a full and durable reset of the defended level. AG's longer fitted on-drug half-life (165 vs 50 days) should be interpreted cautiously — it may reflect genuinely slower adaptation of the defended signal under pharmacologic suppression, or it may partly reflect the fact that the behavioral output used to infer the signal (binge frequency) is itself directly suppressed by the drug.

AJ's finding that the drug suppresses 63% of the falling-phase expenditure bonus is a within-dataset result without clear precedent in the GLP-1 literature, which emphasizes appetite and weight outcomes more than branch-specific expenditure defense.

- Tzang CC, Chen CC, Tsai CF, et al. Metabolic rebound after GLP-1 receptor agonist discontinuation: a systematic review and meta-analysis. *EClinicalMedicine* 80:103680, 2025. [doi:10.1016/j.eclinm.2025.103680](https://doi.org/10.1016/j.eclinm.2025.103680)
- Rubino D, Abrahamsson N, Davies M, et al. Effect of continued weekly subcutaneous semaglutide vs placebo on weight loss maintenance in adults with overweight or obesity: the STEP 4 randomized clinical trial. *JAMA* 325:1414-1425, 2021. [doi:10.1001/jama.2021.3224](https://doi.org/10.1001/jama.2021.3224)
- [Many patients maintain weight loss a year after stopping semaglutide and liraglutide](https://epicresearch-prd.azurewebsites.net/articles/many-patients-maintain-weight-loss-a-year-after-stopping-semaglutide-and-liraglutide). Epic Health Research Network, 2024. (Real-world EHR data showing better maintenance than trials predict — consistent with AG's set point adaptation during treatment. If the set point converges toward current FM while on the drug, discontinuation after sufficient time should produce less regain than trial protocols that stop at peak weight loss.)

## Exercise, NEAT, and RMR

**AD is not well described as a direct contradiction of constrained-energy theory.** Pontzer's constrained-energy model concerns total energy expenditure across broad activity ranges, with compensation at higher activity levels. AD concerns resting metabolic rate and a signal linked to discrete walking sessions rather than total step volume. Those are related but not identical claims. The safer interpretation is that AD identifies an RMR-associated effect of discrete walking sessions that is not captured by a simple "more steps = higher burn" view, and that is not obviously explained by the specific TEE plateau result in Pontzer.

The Levine NEAT paper remains a useful supporting reference, especially if the interpretation is that deliberate walk sessions may trigger broader activity-linked thermogenic responses rather than merely adding their direct exercise cost.

- Pontzer H, Durazo-Arvizu R, Dugas LR, et al. Constrained total energy expenditure and metabolic adaptation to physical activity in adult humans. *Curr Biol* 26:410-417, 2016. [doi:10.1016/j.cub.2015.12.046](https://doi.org/10.1016/j.cub.2015.12.046)
- Levine JA, Eberhardt NL, Jensen MD. Role of nonexercise activity thermogenesis in resistance to fat gain in humans. *Science* 283:212-214, 1999. [doi:10.1126/science.283.5399.212](https://doi.org/10.1126/science.283.5399.212)

## Intake variance and "yo-yo dieting"

**AF does not support the popular claim that day-to-day intake variability is intrinsically metabolically damaging.** Greater short-run intake variance is associated with slightly less fat gain at matched mean surplus (partial r=-0.20), not more. The most plausible interpretation is not that variability itself is beneficial, but that higher variance captures intermittent acute deficits that transiently engage expenditure-side compensation (finding B).

The literature on weight cycling is mixed and much less unanimous than popular discussion implies. Modern reviews generally do not support a simple consensus that weight cycling reliably causes major long-term metabolic damage independent of net weight change. AF is better framed as consistent with the weaker modern literature than as a dramatic overthrow of a settled scientific position.

- Montani JP, Schutz Y, Dulloo AG. Dieting and weight cycling as risk factors for cardiometabolic diseases: who is really at risk? *Obes Rev* 16 Suppl 1:7-18, 2015. [doi:10.1111/obr.12251](https://doi.org/10.1111/obr.12251)
- Schreiner AD, Zhang M, Hu EA, et al. Weight cycling and its cardiometabolic impact: a systematic review. *Nutr Rev* 82:1510-1526, 2024. [doi:10.1093/nutrit/nuad159](https://doi.org/10.1093/nutrit/nuad159)

## Appetite compensation per unit weight lost

**D finds a weaker appetite gradient than the most-cited estimate.** The often-quoted figure of ~100 kcal/day increased appetite per kg of weight lost comes from a model-derived estimate. In this dataset, the distance-to-intake gradient is -30 cal/day per kg below set point — 3x weaker. The difference is consistent with the set point acting through discrete binge events (AG) rather than continuous upward pressure, so part of the signal is missed by looking at mean intake alone.

- Polidori D, Sanghvi A, Seeley RJ, Hall KD. How strongly does appetite counter weight loss? Quantification of the feedback control of human energy intake. *Obesity* 24:2289-2295, 2016. [doi:10.1002/oby.21653](https://doi.org/10.1002/oby.21653)

## Probability of sustained weight loss

**AC documents one case of the base-rate failure.** The 2013 inflection — FM bottomed at 17 lbs and rose every year for a decade — is one case study of the population-level statistic that fewer than 1 in 200 men who pass BMI 30 reach normal weight within 9 years.

- Fildes A, Charlton J, Rudisill C, et al. Probability of an obese person attaining normal body weight. *Am J Public Health* 105:e54-e59, 2015. [doi:10.2105/AJPH.2015.302773](https://doi.org/10.2105/AJPH.2015.302773)

## Summary: what this dataset adds

| Finding | Literature status | This dataset |
|---|---|---|
| Fat mass is the key regulated signal | Long proposed (Kennedy 1953); supported by leptin biology | AG is more consistent with adiposity tracking than scale-weight or gravitostat tracking |
| Defended level can move | Conceptually proposed, usually not operationalized | AG fits a moving defended fat-mass signal better than any fixed set point |
| Timescale of defended-weight adaptation | Dynamic timescales exist in energy-balance models, but not usually as a fitted binge-derived defended signal | AG supports an appetite-related timescale on the order of weeks, centered near 50 days |
| Distance below defended level predicts overeating | Supported indirectly by weight-suppression literature (Lowe 2007) | AG links a model-derived defended signal to binge frequency with a continuous graded curve |
| Expenditure compensation opposes reduced-weight maintenance | Established in ward studies (Leibel 1995, n=18) | K replicates the asymmetry free-living (+179 cal, 1,482 pairs, p<10^-22) |
| Appetite and expenditure may adapt on different clocks | Not explicitly parameterized in prior literature, as far as identified | AI suggests faster expenditure adaptation (≤10d) than appetite adaptation (~50d) |
| GLP-1 discontinuation leads to regain | Supported by discontinuation meta-analyses | AG/AJ are more consistent with suppression than durable reset |
| Discrete walking sessions may influence resting metabolism | Not a standard claim in exercise-energy models | AD suggests a session-linked RMR signal not reducible to step totals (+14 cal/session, p<0.0001) |
| Intake variance is intrinsically harmful | Popular claim; scientific support is mixed | AF shows a small association in the opposite direction at matched mean surplus |
| Appetite compensation is ~100 kcal/day per kg lost | Model-derived estimate (Polidori 2016) | D finds 30 cal/day per kg — 3x weaker, consistent with binge-mediated rather than continuous pressure |
