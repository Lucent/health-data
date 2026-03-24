# Fatty Acid Composition

Three blood fatty acid panels spanning 6 years. All values are % of total fatty acids in whole blood unless noted.

## Trend summary

| Metric | 2019-03-26 (CHL) | 2024-09-05 (OQ) | 2025-09-05 (OQ) | Desirable |
|--------|-------------------|------------------|------------------|-----------|
| Omega-3 Index (RBC) | — | 4.93% | 6.06% | 8–12% |
| Omega-3 total | 2.8% | 4.78% | 5.68% | |
| EPA | 0.3% | 0.54% | 0.45% | |
| DHA | 1.1% | 2.63% | 3.71% | |
| Omega-6 total | 41.3% | 35.80% | 34.83% | |
| Linoleic Acid (LA) | 26.9% | 25.03% | 22.93% | |
| Arachidonic Acid (AA) | 10.9% | 7.52% | 8.93% | |
| AA:EPA ratio | 36.3 | 14.0 | 19.8 | 2.5–11 |
| n6:n3 ratio | 14.8 | 7.5 | 6.1 | 3–5 |
| Trans Fat Index | — | 0.67% | 0.44% | <1% |
| Palmitoleic (carb marker) | — | 0.92% | 0.53% | |

## Key observations

- **Omega-3 improving over time** but still below desirable (8%). Vegetarian with no fish — omega-3 comes from ALA conversion and possibly algal supplements. Genetic profile (FADS1 CC) means efficient desaturase, but conversion of ALA to EPA/DHA is still poor (as expected from the literature).
- **AA:EPA ratio severely elevated** at all 3 time points (36, 14, 20 — desirable <11). The 2019 value of 36:1 is extreme. Improved by 2024 but worsened from 2024→2025 (EPA dropped from 0.54 to 0.45, AA rose from 7.52 to 8.93).
- **n6:n3 ratio improving** (14.8 → 7.5 → 6.1) but still above desirable (3–5). Progress is coming from omega-3 increasing more than omega-6 decreasing.
- **Linoleic acid trending down** (26.9 → 25.0 → 22.9) — possible dietary oil changes.
- **Palmitoleic acid halved** (0.92 → 0.53) between 2024–2025 — a marker of excess carbohydrate intake. The tirzepatide-driven calorie reduction may explain this.
- **Trans fat index low and dropping** (0.67 → 0.44) — consistent with no meat and limited processed food.
- **AA rose** from 2024→2025 (7.52→8.93) despite lower overall omega-6. This is consistent with the FADS1 CC genotype efficiently converting LA to AA. As total LA drops, the conversion machinery has less substrate but the same efficiency, so the AA:LA ratio increases.

## Context for omega-6:3 theories

The genetic profile (FADS1 rs174546 CC — homozygous high-activity desaturase) means this subject converts dietary omega-6 linoleic acid to arachidonic acid more efficiently than most people. The blood AA level is in the normal range but the AA:EPA ratio is severely elevated because EPA is so low (vegetarian, no fish). This is the worst combination for the omega-6:3 ratio theory: maximum conversion of dietary omega-6 to inflammatory AA, minimum intake of counterbalancing EPA/DHA.

## Supplementation timeline

- **2025-08-17**: Started daily 500mg EPA + 500mg DHA (algal oil). The 2025-09-05 OmegaQuant test was only 20 days after starting — too early to reflect in RBC membranes (3-4 month turnover). The 4.93→6.06% Omega-3 Index improvement from 2024→2025 was from dietary changes alone. The next test (~Sep 2026) should capture the full supplementation effect.

## Source files

Data extracted from PDFs. See `fatty_acids.csv` for machine-readable format.
- 2019-03-26: Cleveland HeartLab cardiometabolic panel (via Quest Diagnostics). Also includes lipids, inflammation markers, LDL fractionation.
- 2024-09-05: OmegaQuant Omega-3 Index Complete (7-page report)
- 2025-09-05: OmegaQuant Omega-3 Index Complete (7-page report)
