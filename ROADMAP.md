# Roadmap

## Completed (2026-03-24 session)

### Intake extraction pipeline
- Built `extract_oxps.py`: parses OXPS files using Glyphs X/Y coordinates and horizontal `<Path>` separator lines to reconstruct table cells. Handles orphaned nutrients across page breaks.
- Built `extract_html.py`: parses HTML and MHTML across old layout (2011–2022-09) and new MUI layout (2022-10+). Decodes quoted-printable MHTML, unescapes HTML entities.
- Built `compare_extractors.py`: cross-validates OXPS vs HTML on 1,109 overlapping dates. 100% daily total match (±1 tolerance). 100% item count match on overlapping dates.
- Built `merge.py`: combines all sources with priority waterfall (HTML > MHTML > OXPS) into `intake_foods.csv` (54,925 items) and `intake_daily.csv` (5,074 days).
- Built `sanity_check.py`: Atwater factor validation at item and daily level. 4.9% of items are outliers — split between sugar alcohol products (expected) and bad MFP database entries.

### Weight extraction
- Built `weight/extract.py`: extracts first daily reading from bathroom-scale.xlsx. 1,638 days, 1999–2025.

### Steps and sleep extraction
- Built `steps-sleep/extract.py`: extracts from Samsung Health export using `source_type=-2` deduplication. Verified against 5 monthly totals and 4 individual calibration points. 3,533 step days, 1,279 sleep days.

### Cleanup
- Removed superseded parsers: intake-parse-html (Node.js, calories-only), intake-parsed (Go text dumps + incomplete CSV), intake-scrape (Python API scraper, calories-only), parser-xps-go (Go XPS text extractor), models/ and utils/ (Go dependencies).
- Removed superseded outputs: lm_nutrient_weight_updated_modified_v2.csv, intake/lucent_nutrient_data.csv.
- Moved intake/quarter-depends up a directory and symlinked from 2011, 2012, 2019.

## Next up

### Data gaps to fill
- **2025-11**: Missing MHTML file. Export pending (meals at Lighthaven need to be entered first).
- **2018-02-20 to 2018-02-28**: Can be re-exported from MFP (data available from Feb 20, 2018 onward). Would fill the last OXPS-only gap in the HTML era.
- **Samsung Health**: New export needed to extend steps/sleep past 2024-01-28.
- **Weight**: Recent weigh-ins may not yet be in the spreadsheet. Potential split of bathroom scale export from hand-recorded weights.
- **Cosmed Fitmate RMR**: Additional indirect calorimetry measurements exist beyond the 3 in RMR/. Need to export from device.
- **Withings thermometer**: Two account exports merged into `temperature/temperature.csv`. 1,315 readings, 364 days, Dec 2023 – Mar 2026. Ready to join with intake data.

### Extraction not yet done
- **Body composition** (`composition/Body composition.xlsx`): 49 measurements across 3 eras (BOD POD, InBody partial, InBody full). Needs era-aware column mapping to CSV.
- **Medicine** (`drugs/medicine.xlsx`): 83 entries. Simple extraction to CSV. Critical for modeling — tirzepatide starting 2024-09 is a massive intervention.
- **Workout PDFs** (`workout/`): ~200 gym session PDFs. Need OCR or PDF text extraction to structured sets/reps data.
- **RMR** (`RMR/RMR.xlsx`): Already has 21 measurements (3 lab + 18 home). Small enough to use directly from the XLSX. More readings exportable from the Cosmed Fitmate device.

### Join into analysis table
Build a single daily table joining: intake, weight, steps, sleep, medicine (as binary/dose columns), composition (interpolated). This is the input to all theories in theories.md.

### MFP API micronutrient enrichment
Use food names from intake_foods.csv to query MFP API for iron, potassium, vitamin D, saturated fat, added sugars. Match by name + calorie checksum per the approach in `/myfitnesspal/README.md`. The 1,092 macro-mismatch foods from sanity_check.py are highest priority (their existing data is unreliable).

### Analysis (theories.md) — NEXT UP
**Glycogen-water smoothing** is the immediate next task. Prerequisite for all other analysis.

Inputs:
- `intake/intake_daily.csv` — daily calories and carbs
- `weight/weight.csv` — daily weight (with gaps to interpolate)

Model: each gram of glycogen retains ~3g water. Glycogen depleted over ~2 days of restriction, recaptured within hours of refeeding. Daily carb intake drives glycogen stores; the model predicts water weight fluctuations, which when subtracted from scale weight reveal true fat mass trajectory.

Validation: the Oct-Nov 2019 weekend fasts (36-hour, visible in intake data as <150 cal Sat/Sun) should show weight drops far exceeding caloric deficit, with immediate bounce-back on refeeding days. If the model correctly predicts the bounce magnitude, it's working.

Then:
2. **Set point estimation** — Kalman filter / state-space model extracting the hidden variable
3. **Hunger vs. food noise separation** — compare binge frequency at equivalent deficits pre- vs. post-tirzepatide
4. **Binge prediction from set point distance** — does cumulative gap below set point predict binges better than any dietary variable?
5. **Body temperature as real-time metabolic adaptation** — correlate daily temp with intake deficit
6. **Gravitostat** — daily foot-pounds (weight × steps) predicting next-day intake

## Known data quality notes
- The intake data is genuinely complete. Low-calorie days (10 cal sugar cube, 90 cal coconut water) are real fasting/restriction days, not incomplete logging.
- OXPS `--` (unknown) nutrient values render as `0`. In the merged CSV, these are empty strings from HTML/MHTML sources but `0` from OXPS-only years (2013-2017). This biases cholesterol, sugars, and fiber low in those years.
- 4.9% of food items fail Atwater validation (fat×9 + carbs×4 + protein×4 ≠ calories). These are sugar alcohol products (expected — erythritol subtracted from calories but not carbs) and alcohol (7 cal/g not tracked as any macro). A handful are bad MFP database entries.
- ±1 calorie rounding differences between OXPS and HTML exports of the same data are expected (MFP internal rounding varies between export times).
- Scale body fat % is unreliable. Use /composition for body fat data.
