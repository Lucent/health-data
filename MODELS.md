# Published quantitative models of body-weight / appetite regulation — digest for the health-data project

Reference point for every entry: the project's DG model `intake = I0 + s·B + P·(SP − FM)`, s = −80 kcal/d per exposure unit, **P = +43 kcal/d per lb of FM below a set point frozen at pre-drug FM (= 95 kcal/d per kg)**; pre-drug, no forward fat-mass feedback on intake, but a *lean-mass* deficit at given FM predicts higher intake over the next 90 d (partial r +0.4 to +0.6).

Unit crib: 1 kcal/d per kg = 0.4536 kcal/d per lb. Hall's ρ_F = 39.5 MJ/kg = 9,440 kcal/kg = 4,280 kcal/lb (the project uses 3,500 kcal/lb FM); ρ_L = 7.6 MJ/kg = 1,816 kcal/kg = 824 kcal/lb (project uses 1,100 kcal/lb LM). Any P re-derived under Hall's densities will be ~20 % larger per lb than under the project's.

Sources I could read in full are marked [full]; abstract/secondary-only are marked [abs].

---

## 1. Polidori, Sanghvi, Seeley, Hall 2016 — proportional appetite feedback from the canagliflozin trial [full]

**What it is.** Energy intake is back-calculated from group-mean weight curves (the Sanghvi/Hall method, entry 2) in a 52-wk placebo-controlled canagliflozin 300 mg trial (n = 153 drug, 89 placebo, all T2D). Canagliflozin removes ~90 g glucose/d in urine (UGE ≈ 360 kcal/d) without the patient knowing, so the intake rise that follows is pure physiology. Weight fell "several kg" and plateaued; calculated intake rose to a steady state **~+350 kcal/d** above baseline, exactly offsetting UGE. Placebo intake dipped ~100 kcal/d and was back to baseline by week 15. A proportional controller on weight fits; integral control (zero steady-state error) does not.

**Equations.**
```
ΔEI(t) = ρ·dBW/dt + ε·(BW − BW0) + Δδ·(1−β)·BW0 + UGE        (intake from weight)
ρ = (η_FM + ρ_FM + α·η_FFM + α·ρ_FFM) / ((1−β)(1+α)),  α = C/FM, C = 10.4 kg  (Forbes)
ε = [γ_FM + α·γ_FFM/(1+α) + δ0 + Δδ] / (1−β)
Proportional control:   ΔEI(t) = −k_P · ΔBW(t)
Integral alternative:   ΔEI(t) = −k_P·ΔBW − k_I·∫ΔBW dt   (k_I = 1 kcal/kg/d², rejected)
```
**Parameters (Table 1).** γ_FM 3.2, γ_FFM 22, δ0 10 kcal/kg/d; η_FM 180, η_FFM 230 kcal/kg; ρ_FM 9,300, ρ_FFM 1,100 kcal/kg; β 0.24; UGE 360 kcal/d. **Fitted k_P = 95 kcal/d per kg lost (abstract: "~100")**; no CI is printed in the text (figures show ±95 % CI bands only). Expenditure adaptation for comparison: ~30 kcal/kg/d at 10–20 % loss, so appetite feedback is ">3-fold larger".

**Test here.** k_P = 95 kcal/d/kg is **43.1 kcal/d/lb — numerically identical to P**. Direct replication: regress weekday-adjusted intake on (BW0 − BW) in the drug era with the drug term held at s·B; also fit on *weight* rather than FM to check whether the per-kg-BW and per-lb-FM framings diverge when LM moves.
DOI 10.1002/oby.21653 (PMC5108589; bioRxiv 10.1101/051045).

## 2. Hall's NIH body-weight model (Hall et al. 2011 Lancet; Hall 2010 AJP; Chow & Hall 2008) and the Sanghvi/Hall intake-from-weight method [full appendix]

**Structure.** Two-compartment (F, L) energy partition, plus glycogen and ECF for the first weeks. No appetite term at all: EI is an input; the "intake-estimation" method is the model run backwards.
```
ρ_F dF/dt = (1−p)(EI − EE − glycogen flux),  ρ_L dL/dt = p(…),  p = C/(C+F), C = 10.4 kg·ρ_L/ρ_F
EE = K + γ_F·F + γ_L·L + δ·BW + TEF + AT + η_F dF/dt + η_L dL/dt
TEF = β_TEF·ΔEI (β_TEF = 0.1);   τ_AT dAT/dt = β_AT·ΔEI − AT,  β_AT = 0.14, τ_AT = 14 d
Glycogen: dG/dt = CI − k_G·G², ρ_G = 17.6 MJ/kg, 2.7 g water/g G, baseline ~500 g
ECF: d(ECF)/dt = (ΔNa_diet − ξ_Na·ΔECF − ξ_CI·(1 − CI/CI_b))/[Na];  [Na] 3.22 mg/ml, ξ_Na 3000 mg/L/d, ξ_CI 4000 mg/d
Linearised long run: ρ·dΔBW/dt = (1−β)ΔEI − ε·ΔBW,  β = β_AT + β_TEF,  τ = ρ/ε (half-time ≈ 1 yr)
```
**Parameters.** ρ_F 39.5, ρ_L 7.6 MJ/kg; γ_F 13, γ_L 92 kJ/kg/d (3.1 / 22 kcal/kg/d); η_F 750, η_L 960 kJ/kg; δ ≈ 30 kJ/kg/d (PAL 1.5). Adaptive thermogenesis is **0.14 × ΔEI, 14-d time constant** — i.e. a 1,000 kcal/d cut costs 140 kcal/d of EE beyond composition.

**Sanghvi et al. 2015 (CALERIE, n = 140, 2 y):** model ΔEI vs DLW/DXA ΔEI, per-subject RMSD 215 kcal/d, most within 132 kcal/d, group mean within 40 kcal/d.

**Test here.** (a) The project's P1/P4 glycogen+sodium water corrections are the same construction; compare ξ_Na, ξ_CI and the k_G·G² form against P0's tuned constants. (b) AT: regress calorimetry RMR residual (vs composition) on a 14-d EMA of ΔEI; Hall predicts slope 0.14 (on RMR alone, lower). (c) Run Hall's inversion on the subject's weight alone and compare the inferred ΔEI with the *logged* ΔEI — this dataset is the validation set that Sanghvi lacked at the individual level.
DOIs 10.1016/S0140-6736(11)60812-X (PMC3880593; appendix at niddk.nih.gov/-/media/Files/BWP/Hall_Lancet_Web_Appendix.pdf); 10.1152/ajpendo.00559.2009; 10.1371/journal.pcbi.1000045; Sanghvi 10.3945/ajcn.115.111070.

## 3. Dulloo's protein-stat / dual-feedback model and "collateral fattening" (Dulloo, Jacquet & Girardier 1997; Dulloo 1997 PNS; Dulloo 2017; Jacquet et al. 2020) [1997 PNS full; 2020 full]

**Evidence.** Minnesota Experiment re-analysis, 12 men (24 wk semistarvation: −25 % BW, −70 % fat, −27 % FFM; 12 wk restricted then 8 wk ad-lib refeeding). Hyperphagia index = ∫(EI − control EI) over the 8 ad-lib weeks. Simple and partial correlations: hyperphagia is inversely related to **both** fat recovery and FFM recovery at R12, "and the latter are two independent determinants"; the correlation with each "persists after adjusting (by partial correlation) for … the other". Stepwise regression: (1) fat recovery (strongest), (2) prior intake deficit, (3) FFM recovery — together ~80 % of inter-individual variance. Hyperphagia "persisted until FFM were completely restored", after fat and weight had fully recovered — that is the "collateral fattening" / fat-overshoot mechanism. No kcal/d-per-kg-FFM coefficient is printed anywhere; the 1997 AJCN paper's regression slopes are in its figures only.

**Formal model (Jacquet et al. 2020, weight cycling).**
```
P_SS = 1/(1 + ΔFAT_SS/ΔFFM_SS)  (lean fraction of loss);  P_RF likewise for regain
γ = P_SS/P_RF;   FAT_overshoot = (γ − 1)·ΔW_SS          (exact, by construction)
P_SS(%FAT0) = [(100 − %FAT0)/100]·exp(−0.015·%FAT0)
γ = 1 + 0.92·exp(−0.11·%FAT0)   →  overshoot ≈ 0.92·e^(−0.11·%FAT0)·ΔW_SS
```
At 10/20/30 % body fat a 5 kg cycle leaves ≈1.5 / 0.5 / 0.2 kg of extra fat. Minnesota overshoot averaged ~4 kg (0–9). Blundell/Hopkins add the "passive" FFM pathway: FFM → RMR → daily EI (positive, mediated ~96 % by RMR), FM weakly negative or nil (Hopkins 2016 IJO; Hopkins & Blundell 2023).

**Test here.** This is the only published model with the sign of the project's lean-mass finding. Fit `intake = I0 + k_F·(FM_ref − FM) + k_L·(LM_ref − LM)` on the pre-drug record, where LM_ref is the LM at the same FM on the subject's own Forbes curve; report k_L in kcal/d per lb and whether k_F ≈ 0 pre-drug while k_L > 0 (project's finding) or k_F > k_L (Dulloo's ranking). Also check the post-fast (weekend-fast era) and post-restriction overshoot magnitudes against γ(%FAT0) with the 70 scans.
DOIs 10.1093/ajcn/65.3.717 (1997 AJCN); 10.1079/PNS19970007 (1997 PNS); 10.1038/ejcn.2016.256 (2017 EJCN); 10.1002/oby.21734 (2017 Obesity "collateral fattening"); 10.1038/s41366-020-0547-1 (Jacquet 2020, PMC7260129); Hopkins 10.1038/ijo.2015.155.

## 4. Speakman: settling point, dual-intervention-point (DIP) and the drifty gene (Speakman et al. 2011; Speakman 2018; Speakman & Hall 2023) [2011 full; 2018 full; 2023 partial]

**Settling point:** no reference; one flow independent of stores, the other scaling with them (EE ≈ +20–30 kcal/d per kg); a 300 kcal/d cut produces exponential loss taking years; Minnesota's rapid regain is cited as evidence *against* a pure settling point. **Control-theory framing (2023):**
```
e(t) = SP − PV;   u(t) = K1·e(t) + (K2/Ti)·∫e dt;   human K1 ≈ 100 kcal/kg/d, K2 = 0 (from the SGLT2 study)
```
gain → 0 gives dynamic equilibrium, strong gain gives set point; weak proportional-only gain yields "settling points" far from SP. **DIP:** two independently set thresholds, LIP (starvation risk) and UIP (predation risk); "no set-point at all within the zone of indifference"; weight drifts passively with environment between them; UIP alleles drifted after predation release ~2 Ma ("drifty gene"). Quantitative DIP (2018 Mol Metab, mouse-scaled):
```
M_s = a·e^(−b·x),  M_p = c·e^(g·x),  M_tot = M_s + M_p,  optimum  x* = ln(ab/(cg))/(b+g)
```
with 1 % CV on (a,b,c,g) the optimum's CV is 20–90 %, so a single set point is "highly unlikely to evolve"; the DIP zone = fatness range where M_tot is within 1 % of its minimum (e.g. 4.0–5.4 g in scenario B). Polygenic effects: 0.17 kg/m² per allele (range 0.06–0.39), "well below 150 g/allele".

**Distinguishing in longitudinal data (2023):** set point predicts reactive hyperphagia after any loss; DIP predicts it only once FM has been pushed below LIP; within the zone, no forward feedback. Individual-level series, not cohort means, are needed (SD of yearly change 2.6 kg in 192 people).

**Test here.** The subject's pattern — no FM feedback across 2011–2024 while FM wandered 20→85 lb, then P = 43 kcal/d/lb once the drug pushed FM below pre-drug — is exactly DIP's signature **if** the pressure switches on at a fixed threshold rather than scaling from SP0. Fit `intake = I0 + s·B + P·max(0, LIP − FM)` vs the DG linear form and compare AIC; the LIP would be ~pre-drug FM only if the subject sat on his own LIP for a decade, which is DIP's weak point. (README already rejects a *zone* in the EMA framing; this re-tests it in the frozen-SP framing.)
DOIs 10.1242/dmm.008698 (2011 DMM); 10.1016/j.molmet.2017.10.007 (2018, PMC5784320); 10.1098/rstb.2022.0231 (2023, PMC10475878); 10.1242/jeb.167254 (2018 JEB).

## 5. Leibel / Rosenbaum: metabolic adaptation and the leptin threshold [1995 abs; 2005 full; 2010 full]

**Numbers.** Leibel 1995 NEJM (in-patient, weight clamped by formula titration): at −10 % BW, TEE falls **6 ± 3 kcal/kg FFM/d (never-obese) and 8 ± 5 (obese)**; at +10 %, rises 9 ± 7 / 8 ± 4. REE and NREE each fall ~3–4 kcal/kg FFM/d. For a 65-kg FFM that is ≈ −400 to −520 kcal/d TEE, of which REE ≈ −200 to −260. Rosenbaum 2008 AJCN: the deficit is the same in people who have held −10 % for >1 y. Rosenbaum & Leibel 2010 table (−10 % maintained): TEE −15 %, NREE −30 %, SNS tone −40 %, T3 −7 %, muscle efficiency +20 %; "~300–400 fewer calories per day" than a never-reduced person of the same composition. **Leptin repletion (Rosenbaum 2005 JCI, n = 10):** 0.068 ± 0.01 mg/kg/d s.c. (restoring 8 am leptin to pre-loss levels) **reversed** the TEE/NREE deficit ("300–500 kcal/d greater than predicted"), the +23 % muscle efficiency at 10 W, SNS and T3 changes; REE/TEF were not significantly beyond composition.

**Model (verbal).** Leptin acts as a threshold: below the individual's leptin "floor" the full adaptation package (and hunger) engages; above it, extra leptin does little (weight loss needs >10× normal leptin). No equation published; Tam et al. 2009 formalise it (entry 9).

**Test here.** 25 calorimetry RMRs + Kalman FM/LM: (a) regress RMR residual on fractional weight change from the local max — Leibel predicts ≈ −3 to −4 kcal/kg FFM/d per 10 % loss, i.e. ≈ −200 to −250 kcal/d, which is the size of the project's "−206 on drug" and the ~90 kcal/d loss/gain asymmetry (AI/AO); (b) does the deficit persist in the stable 2016–2019 years (Rosenbaum 2008 says yes)? (c) the 2013 inflection at FM ≈ 18 lb (AC/BK) is a leptin-floor prediction.
DOIs 10.1056/NEJM199503093321001; 10.1093/ajcn/88.4.906; 10.1172/JCI25977 (PMC1297250); 10.1038/ijo.2010.184 (PMC3673773).

## 6. The weight-loss plateau: appetite feedback (Hall 2024) vs intake creep / adherence (Thomas 2014; Polidori lifestyle arm) [Hall 2024 full preprint; Thomas full]

**Hall 2024 Obesity** ("Physiology of the weight-loss plateau…"):
```
ΔI(t) = −P(t) − k·ΔW(t);   P(t) ramps P0→Pmax over T (dose escalation), then constant
τ = ρ/(k + ε),  ε ≈ 25 kcal/d/kg;  plateau time depends on k (and ε), "nothing to do with P"
```
Fitted k (kcal/d per kg): **SGLT2/baseline 95; CALERIE diet 82–83 (P = 830); semaglutide 2.4 mg 49 (P 610→1300 over 16 wk); tirzepatide 10 mg 48 (P 830→1560 over 12 wk); RYGB 58 (P = 3600)**. Plateau ≈ 12 mo for diet, ≈ 24 mo for drugs/surgery, explained by 40–50 % weaker feedback, not by expenditure. The model does not predict regain without extra assumptions.

**Thomas et al. 2014 AJCN:** plateau at ~6 mo reproduced by intermittent non-adherence (daily EI = prescription when U(0,1) > adherence, else deviates); adherence matched to CALERIE: women 80→60→50→40 %, men 80 % (mo 1–5) → 70 %; a 5–10 % EE reduction shifts plateau *weight* but moves plateau *time* to 365–730 d. **Polidori's lifestyle arm:** EI −700 kcal/d at start, exponential relaxation; within 100 kcal/d of baseline at the 6–8 mo plateau; at 1 y back to baseline with weight still −5 kg.

**Test here.** The trials cannot separate "appetite creep" from "adherence creep" because intake is inferred from weight; this dataset logs intake, so the two are separable: (a) fit Hall's form with measured ΔI on the drug era and read k directly (project gets 95/kg-equivalent, not 48); (b) during the willpower era (2011–13) test Thomas's generative model — is the day-level intake distribution bimodal (adherent days at prescription + a nonadherent tail) with adherence decaying monthly, or unimodal with a drifting mean? BK's "impossible days" analysis is the right place.
DOIs 10.1002/oby.24027 (bioRxiv 10.1101/2023.11.05.565699, PMC10705578); 10.3945/ajcn.113.079822 (PMC4135489); Hall & Guo 2017 Gastroenterology 10.1053/j.gastro.2017.01.052 (EE slope 20–30 kcal/d/kg; intake slope ≈ −100 kcal/d/kg; energy gap 200–250 kcal/d).

## 7. Lowe's weight suppression (WS) and binge pressure (Lowe 2007; Lowe et al. 2011; Keel group 2025) [2011 full; 2007 abs]

**Model.** WS = highest past adult weight − current weight. Biobehavioural chain (Keel/Lowe): WS → lower leptin (beyond what BMI predicts) → blunted post-prandial GLP-1 → reduced reward satiation → larger, out-of-control binges; a 2025 SEM in bulimic syndromes reports excellent fit to that chain (binge *size*, not frequency). Quantitative results are regressions, not dynamics: Lowe 2007 (n = 182 BN): WS directly, and dietary restraint inversely, related to objective binge frequency. Lowe 2011 (n = 110 BN, 8-y follow-up): Cox HR for first full remission **0.89 per unit WS (95 % CI 0.82–0.97, p = .01)**, robust to current BMI, highest weight and weight change (HR 0.88 adjusted). No equation for intake as a function of WS and no gain in kcal/d per kg is published.

**Test here.** The project has already shown a fixed-maximum reference is dominated by a moving 50-d EMA (README), and that binge *size* is constant while *frequency* moves. Remaining Lowe-specific tests: (a) WS = 2011 peak weight (~100 lb above) as the reference for the willpower era — Lowe predicts binge onset scales with WS, the project found onset only when FM crossed ~18 lb (a floor, not a suppression distance); (b) Keel's chain predicts binge *size* should track WS — the data say constant size, which is a disconfirmation worth stating.
DOIs 10.1002/eat.20405 (2007 IJED); 10.1037/a0025714 (2011 JCCP); Keel 2025 PMC12234017.

## 8. GLP-1/GIP plateau and regain models with an appetite term [Hall 2024 full; Horn 2025 abs; 2026 modelling secondary]

**Hall 2024** is the only published model with an explicit appetite-feedback term on incretin drugs (entry 6): fixed drug effect P(t), proportional gain k ≈ 48–49 kcal/d/kg, set point implicitly frozen at baseline weight (ΔW measured from pre-treatment), no erosion/tachyphylaxis term; discontinuation = P → 0 and return toward baseline on τ ≈ ρ/(k+ε). **Speakman & Hall 2023** note that the same model gives ≈ 4 y for a semaglutide-withdrawal placebo arm to return within 1 % of baseline. **Horn et al. 2025 (SURMOUNT-1/-4):** empirical time-to-plateau — longer for class II/III than overweight, most participants plateaued by week 72, later with 10/15 mg, younger age and female sex; no appetite parameters. **2026 modelling report (Medscape, primary not located):** Hall model applied to SURMOUNT-1 weight data with 1,000 Monte-Carlo runs: EI −32.1 % in 4 wk (2,761 → 1,875 kcal/d), back to 2,500 kcal/d by week 98 where it equals EE (−9.2 %, 2,509 kcal/d, flat to week 176); predicted regain 5.7 kg in 17 wk after stopping at week 176. The project's own README decomposition (frozen SP, P ≈ 30–43 kcal/d/lb, plateau where s·B = P·(SP0 − FM)) is structurally the same equation as Hall 2024 with P(t) = −s·B(t).

**Test here.** (a) Hall's cohort k = 48 kcal/d/kg (22/lb) vs the subject's 95/kg (43/lb): the project's DG already finds cohort pressure 15–25 kcal/d/lb reproduces the trial regain — so the *subject-vs-cohort gap replicates Hall's own result*; ask whether Hall's lower on-drug k is an artefact of forcing a constant P when the real s·B is linear in exposure and the trial dose–response is concave (DG's open problem). (b) Hall's intake trajectory for SURMOUNT-1 (−32 % at 4 wk, recovering to EE by ~2 y) is a direct template for the subject's weekly-mean logged intake at matched exposure. (c) Horn's covariates (dose, age, sex) → the subject's plateau time per dose era.
DOIs 10.1002/oby.24027; 10.1111/cob.12734 (PMC12096058); 10.1098/rstb.2022.0231.

## 9. Other quantitative frames

**Tam, Fukumura & Jain 2009 (Cell Metab) — leptin settling point vs set point vs threshold [full].** Mouse-parameterised but transferable in form:
```
Settling: FI = k4·(1 − L_b/(k5 + L_b)),  E_out = k6·BM·(1 + k7·L_b/(k8 + L_b)),  dE/dt = ρ_food·FI − E_out
          leptin: d(L_p·V)/dt = FM·R_syn − GFR·Cl·L_p;  brain L_b = k1·L_p/(k2+L_p) + k3·L_p
Set point: FI = a1·e + a2·∫e + c1,  E_out = BM·(a3·e + a4·∫e + c2),  e = L_b − SetPt   (PI control)
Unified:   integral term × H(L_thresh − L_b)  — integral (zero-error) defence only BELOW a leptin threshold
```
k4 5.6 g/d, k5 0.55 ng/g, k6 244 cal/g/d, k8 0.22 ng/g, R_syn 51.8 ng/g fat/d; a1 −0.24, a2 −288, a3 7.2 cal/ng/d, a4 86.4 kcal/ng/d². Pure set point cannot produce diet-induced obesity; the threshold model (defend below, drift above) reproduces it and gives multiple steady states under leptin resistance (k2 rises above L_p = 10 ng/ml). **This is the formal version of the project's "dual defence / leptin floor" reading**: test whether intake pressure below SP is better fit by proportional (P·gap) or proportional-plus-integral (∫gap) control — the integral form predicts intake keeps rising at constant gap during the 15-wk FM plateau, which the project says did *not* happen (intake fell 185 kcal/d). DOI 10.1016/j.cmet.2008.11.005 (PMC2643084).

**Flatt glycogenostat (Flatt 1987 AJCN; Flatt 2004 Obes Res two-compartment model) [abs].** Intake is regulated to hold hepatic/muscle glycogen in a habitual band; fat mass settles where fat oxidation = fat intake, so steady-state FM rises with dietary fat fraction and with the glycogen band. 2004 simulation: a 10 % shift in dietary-fat energy shifts steady-state body fat by the amount offset by ~26 g of mean glycogen. Human evidence: previous-day carbohydrate balance explains ~5.5 % of next-day energy balance. **Test here:** P1 already carries a daily glycogen state — regress next-day intake (and next-day carbohydrate intake) on today's glycogen deviation, and test whether keto/potato epochs (near-zero vs saturating carb balance) show the predicted intake shift. Flatt 1987 doi 10.1093/ajcn/45.1.296; 2004 Obes Res 12:2013–2022 (PMID 15687403).

**Brobeck 1948 thermostat [abs].** Intake falls with heat load and rises with cold ("animals eat to keep warm and stop eating to prevent hyperthermia"); the specific dynamic action of food is the satiety signal. No parameters. **Test:** regress intake on the baseline-corrected body-temperature series (AS/AT) and on ambient temperature; the project's ~97.7 °F mean is a standing deviation this model would read as a chronic hunger signal. Yale J Biol Med 1948;20:545–552 (PMC2602371).

**Pontzer 2016 constrained TEE [full].** DLW TEE vs accelerometry in n = 332: below 230 CPM/d (95 % CI 44–428) TEE = 1.12·CPM + 2,336 kcal/d; above it slope 0.21 ± 0.35 (p = 0.54), plateau ≈ 2,600 kcal/d; RMR uncorrelated with activity. Cross-sectional only. The project's AD/AK (walk sessions raise *resting* rate) is already the within-person counter-test; the remaining check is whether the subject's TDEE–steps relation has a breakpoint. DOI 10.1016/j.cub.2015.12.046.

**Energy flux (Hume, Yokum & Stice 2016 AJCN) [abs].** In 154 adolescents and 75 young women with DLW, *low* energy flux (low intake + low expenditure), not surplus, predicted future fat gain; higher flux protected partly via higher RMR. No model equation. **Test:** classify months by flux (intake + TDEE) at matched balance and look at subsequent 90-d FM change; the project's "not all restriction is equal" (high-step restriction penalised) is a related but opposite-sign observation worth reconciling. DOI 10.3945/ajcn.115.127753.

**MC4R / FTO and a genetic set point [abs].** No published quantitative set-point model exists for carriers; the literature offers effect sizes only: ~0.17 kg/m² per common risk allele (FTO ~0.39), MC4R heterozygous loss-of-function → hyperphagia at test meals and higher fat mass (Farooqi 2003 NEJM 10.1056/NEJMoa022050); Speakman's DIP reads these as UIP-shifting alleles. The subject's genotype (FTO homozygous ×5, MC4R het, UCP2) makes him a single-case test of whether the pre-drug "no feedback" decade is a high/absent UIP.

**Hall–Guo "sliding set point" (in Speakman & Hall 2023).** EE line rising ≈ proportionally with weight; intake response line negative in (W − W*); W* "not fixed but able to move up and down the expenditure line" with the food environment. No numbers; it is the verbal form of the project's SmoothLatch.

---

## Ranked shortlist — the five models most worth fitting

1. **Polidori/Hall proportional-control model (entries 1, 6, 8).** Same equation as DG; k_P = 95 kcal/d/kg is P = 43 kcal/d/lb to the first decimal. Fit it with *measured* intake (which no published fit has had), report τ = ρ/(k+ε), and resolve whether on-drug k really falls to ~48 (Hall) or stays at 95 with a linear drug term (DG). Adds the one missing piece in DG: Hall's AT term (0.14·ΔEI, τ 14 d) for the expenditure side.
2. **Dulloo/Jacquet dual-feedback (adipostat + proteinstat) (entry 3).** The only model with the lean-mass sign the project found; Minnesota gives a cohort anchor (fat strongest, FFM independent, ~80 % R²) and the γ(%FAT0) overshoot equation gives a second, composition-based prediction testable with 70 scans. Worth fitting k_L explicitly in kcal/d per lb and checking against the Hopkins/Blundell *passive* FFM→RMR→intake pathway, which has the opposite sign.
3. **Tam 2009 threshold (PI-below-threshold) model (entry 9) = Leibel leptin floor (entry 5).** Formalises "defend below, drift above"; discriminates proportional vs integral pressure using the 15-wk constant-gap plateau and the 2013 inflection at FM ≈ 18 lb. Calorimetry gives the expenditure half: Leibel predicts ≈ −3–4 kcal/kg FFM/d of REE per 10 % loss, which is the size of the observed −206 kcal/d.
4. **Speakman DIP (entry 4).** The one framework in which "no fat-mass feedback for a decade, then P = 43 once pushed below pre-drug FM" is the *expected* pattern rather than a puzzle; but it requires the pressure to switch on at a fixed LIP, not scale from SP0. A piecewise `max(0, LIP − FM)` vs linear `(SP0 − FM)` comparison on the drug era is cheap and decisive.
5. **Thomas 2014 adherence model (entry 6).** The null for the willpower era: random non-adherence with monthly-decaying adherence reproduces a 6-month plateau with no physiology. BK already found a one-stock reservoir fails; Thomas's two-state generative model is the natural next null before invoking set-point pressure in 2011–13.

(Runner-up: Flatt's glycogenostat, because P1's glycogen state is already computed and the test is one regression.)

## Where the literature agrees / disagrees with P = 43 kcal/d/lb and the lean-mass finding

- **Agrees:** Polidori 2016 k_P = 95 kcal/d/kg (43.1/lb), Hall & Guo 2017 "≈ 100 kcal/d per kg", Speakman & Hall 2023 K1 ≈ 100 kcal/kg/d — all per kg *body weight* in people losing ~75 % fat; the subject's number is per lb *fat mass* in 93 %-fat loss, so the framings coincide on-drug and diverge only where LM moves, which is exactly where the lean finding lives.
- **Disagrees (factor ~2):** Hall 2024 on-drug k = 48–49 kcal/d/kg (22/lb) for tirzepatide/semaglutide, 58 for RYGB. The project's cohort fit (15–25 kcal/d/lb) lands on Hall's values, so the disagreement is subject-vs-cohort, not model-vs-model. Hall also places the set point at baseline weight with no erosion term — identical to DG's frozen SP — and his model does not predict regain without additional assumptions, whereas DG's frozen SP reproduces SURMOUNT-4's shape.
- **Expenditure side:** Hall ε ≈ 25 kcal/d/kg (11/lb) and AT = 0.14·ΔEI; Leibel ≈ −6 to −8 kcal/kg FFM/d TEE at −10 % (REE share ≈ −200–260 kcal/d), persistent >1 y, leptin-reversible. The project's −206 kcal/d on drug and ~90 kcal/d loss/gain asymmetry sit inside these ranges; nobody has published the drug *suppressing* the falling-phase bonus (AJ/AO) — that is new.
- **Lean mass:** Dulloo 1997 is the only quantitative support (FFM deficit an independent determinant of hyperphagia after partialling fat; fat the stronger of the two; no kcal/kg coefficient published). Hall, Polidori, Speakman and Tam have no lean term at all. Hopkins/Blundell find the *opposite* sign cross-sectionally (more FFM → more intake, via RMR), so the project's within-person deficit effect (partial r +0.4–0.6) and their between-person level effect are not in conflict but must be stated separately. Dulloo's ranking (fat > FFM) is contradicted by the subject's pre-drug decade (FFM only); DIP would explain that as "inside the zone for fat, but no zone for lean".
- **Timescales:** Hall's τ = ρ/(k+ε) ≈ 65–105 d for k = 48–95 and ε = 25 is the same order as the project's 45–50 d eating-channel estimate (which DD then shows is mostly reverse causation pre-drug). Hall's AT τ = 14 d matches the project's ~9-d expenditure channel.

### Primary sources (DOIs)
Polidori 2016 10.1002/oby.21653 · Sanghvi 2015 10.3945/ajcn.115.111070 · Hall 2011 10.1016/S0140-6736(11)60812-X · Hall 2010 10.1152/ajpendo.00559.2009 · Chow & Hall 2008 10.1371/journal.pcbi.1000045 · Hall & Guo 2017 10.1053/j.gastro.2017.01.052 · Hall 2024 10.1002/oby.24027 · Thomas 2014 10.3945/ajcn.113.079822 · Dulloo 1997 AJCN 10.1093/ajcn/65.3.717 · Dulloo 1997 PNS 10.1079/PNS19970007 · Dulloo 2017 10.1038/ejcn.2016.256 · Dulloo 2017 Obesity 10.1002/oby.21734 · Jacquet 2020 10.1038/s41366-020-0547-1 · Hopkins 2016 10.1038/ijo.2015.155 · Speakman 2011 10.1242/dmm.008698 · Speakman 2018 10.1016/j.molmet.2017.10.007 · Speakman & Hall 2023 10.1098/rstb.2022.0231 · Leibel 1995 10.1056/NEJM199503093321001 · Rosenbaum 2005 10.1172/JCI25977 · Rosenbaum 2008 10.1093/ajcn/88.4.906 · Rosenbaum & Leibel 2010 10.1038/ijo.2010.184 · Lowe 2007 10.1002/eat.20405 · Lowe 2011 10.1037/a0025714 · Horn 2025 10.1111/cob.12734 · Tam 2009 10.1016/j.cmet.2008.11.005 · Flatt 1987 10.1093/ajcn/45.1.296 · Pontzer 2016 10.1016/j.cub.2015.12.046 · Hume 2016 10.3945/ajcn.115.127753 · Farooqi 2003 10.1056/NEJMoa022050.

### Access notes
Read in full: Polidori (PMC), Hall 2024 preprint (PMC), Hall 2011 appendix (NIDDK PDF), Jacquet 2020 (PMC), Dulloo 1997 PNS (Cambridge PDF), Speakman 2011 (Biologists), Speakman 2018 (PMC), Speakman & Hall 2023 (PMC, first half), Thomas 2014 (PMC), Hall & Guo 2017 (PMC), Rosenbaum 2005/2010 (PMC), Tam 2009 (PMC), Pontzer 2016 (PMC), Chow & Hall 2008 (PMC), Lowe 2011 (Drexel PDF). Abstract/secondary only: Leibel 1995, Flatt 2004, Lowe 2007, Horn 2025, Dulloo 1997 AJCN, Hume 2016, and the 2026 SURMOUNT-1 Hall-model report (Medscape paywalled; primary paper not found). Polidori prints no CI for k_P in the text; Dulloo prints no kcal-per-kg coefficients anywhere.
