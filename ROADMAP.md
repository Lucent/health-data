# Roadmap

## Completed

### Data extraction
- `intake/merge.py` — OXPS + HTML + MHTML → `intake_foods.csv` (54,925 items), `intake_daily.csv` (5,429 days). Cross-validated: 100% match on 1,109 overlapping dates.
- `weight/extract.py` → `weight.csv` (1,693 days, 1999-2026)
- `steps-sleep/extract.py` → `steps.csv` (4,275 days), `sleep.csv` (2,057 days)
- `composition/extract.py` → `composition.csv` (49 measurements, 3 eras)
- `RMR/extract.py` → `rmr.csv` (21 calorimetry measurements)
- `drugs/extract.py` → `medicine.csv`, `tirzepatide.csv` (PK blood level + tachyphylaxis)
- Removed superseded parsers (Node.js, Go, Python scraper)

### Analysis pipeline
- `analysis/P1_glycogen_smooth.py` — Glycogen + sodium water-weight corrections. 14.3% variance reduction.
- `analysis/P2_rmr_model.py` — Individual RMR coefficients from 21 calorimetry × 49 composition measurements. Linear FM/FFM interpolation to all 5,429 days.
- `analysis/P3_interpolate_weight.py` — Window-based TDEE derivation, complete daily weight series.
- `analysis/P4_kalman_filter.py` — Kalman fat mass + TDEE with uncertainty bands, now emitting both causal filtered and retrospective smoothed estimates. 152× smoother TDEE than window method.
- `analysis/plot_models.py` — Diagnostic plots.

### Key findings
Documented in [README.md](README.md) and [FINDINGS.md](FINDINGS.md). Each finding has a `Command:` line in FINDINGS.md pointing to the standalone script that reproduces it.

## Data gaps to fill
- **2025-11**: Missing MHTML file. Meals at Lighthaven need entering.
- **2018-02-20 to 2018-02-28**: Re-exportable from MFP. Last OXPS-only gap in HTML era.
- **Cosmed Fitmate RMR**: Additional measurements exportable from device.
- **Temperature**: 1,315 readings ready to join. Not yet extracted to analysis pipeline.
- **Workout PDFs**: ~200 gym sessions. Need OCR/PDF extraction.

## Next analysis steps
1. **Exploit item-level purity, not generic proxies** — this is a full 15-year food-sequence dataset, not ordinary diary epidemiology. Prioritize analyses that use actual food identity, repetition, regime boundaries, and transition structure rather than weak defaults like weekday/weekend.
2. **Food monotony / novelty / repetition structure** — quantify actual monotony from `intake_foods.csv`: distinct-food count, recurrence interval, menu compression, transition entropy, and meal-template persistence. Test whether these predict binge persistence, post-restriction rebound, and metabolic failure better than calorie variance alone.
3. **Regime-switch transition models** — treat `potato_diet`, `keto_phase`, `travel_binge`, baseline periods, and tirzepatide as distinct transition systems rather than mean-intake eras. Compare `typical -> high`, `high -> binge`, and `binge -> binge` across named regimes.
4. **Branch memory with stronger controls** — continue the main thread by testing whether the rising/falling branch effect survives additional controls for travel exposure, sequencing of prior cuts/refeeds, item-level monotony, and meal-pattern structure.
5. **Food-template effects at matched macros/calories** — identify recurring meal or day templates and ask whether some systematically produce better or worse fat/TDEE outcomes than calorie-equivalent alternatives.
6. **Restriction run sequencing** — quantify prior 30-day deficit, days since last cut, prior refeed intensity, and recent high-intake streaks. Test whether these explain restriction success or failure beyond branch phase.
7. **Session-level walk regimes** — use `steps-sleep/exercises.csv` instead of blunt daily steps. Separate deliberate daylight walk pairs from incidental high-step days, validate ambiguous Samsung labels against the app, and test whether walk structure or timing is what matters.
8. **Temperature by phase** — still worth doing, but only after stronger historical controls are exhausted. Current overlap limits a clean pre-drug falling-phase thermal test.
9. **MFP API enrichment** — micronutrients for top 200 foods (iron, potassium, vitamin D, sat fat, PUFA/MUFA) once the sequence/regime work above is mined.
10. **Join table** — single daily and run-level analysis tables with all derived variables, especially item-sequence metrics and exercise-session metrics, so future theory tests stop rebuilding ad hoc joins.
