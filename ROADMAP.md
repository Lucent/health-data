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
- **Samsung Health**: Export current through 2026-03-24 (4,275 step days, 2,057 sleep days).
- **Weight**: Recent weigh-ins may not yet be in the spreadsheet. Potential split of bathroom scale export from hand-recorded weights.
- **Cosmed Fitmate RMR**: Additional indirect calorimetry measurements exist beyond the 3 in RMR/. Need to export from device.
- **Withings thermometer**: Two account exports merged into `temperature/temperature.csv`. 1,315 readings, 364 days, Dec 2023 – Mar 2026. Ready to join with intake data.

### Body composition extraction — DONE
Built `composition/extract.py` → `composition/composition.csv`. 49 measurements across 3 eras (bodpod, inbody_partial, inbody_full). Era-aware column mapping: core fields always populated, extended fields (segmental, visceral, water compartments) only for InBody eras.

### RMR extraction — DONE
Built `RMR/extract.py` → `RMR/rmr.csv`. 21 measurements (3 lab Cosmed, 18 home Cosmed Fitmate). More readings exportable from device.

### Composition-aware RMR model — DONE
Built `analysis/rmr_model.py`. Fits individual-specific RMR coefficients (RMR = 32.5×FFM_kg + 13.0×FM_kg - 524) from 21 calorimetry measurements against composition-interpolated FM/FFM. Forbes curve partitions daily weight changes into fat/lean between composition anchors. RMSE=151 kcal/day (at instrument noise floor). Outputs `analysis/daily_composition.csv` (5,429 days: FM, FFM, fat%, expected RMR). Used by interpolate_weight.py for long-gap TDEE derivation.

### Medicine extraction + tirzepatide pharmacokinetics — DONE
Built `drugs/extract.py` → `drugs/medicine.csv` (89 entries) + `drugs/tirzepatide.csv` (560 daily rows). Dose escalation: 2.5→5→7.5→10→12.5mg over 80 weekly injections, Sep 2024 to Mar 2026.

**Pharmacokinetic blood level model** using FDA label parameters (t½=5.0 days, Tmax=24h, ka=3.31/day). Sums one-compartment SC absorption curves from all prior injections, capturing weekly sawtooth and multi-week accumulation (steady state ~week 5, accumulation ratio 1.6×).

**Blood level → intake: r = -0.50** (partial, controlling for time trend). The weekly sawtooth directly tracks appetite: injection day = 1652 cal, day 5 (trough) = 2220 cal — a 568 cal/day swing within each week. Dose-response: each mg reduces intake ~35 cal/day.

**Tachyphylaxis (dose tolerance)**: effectiveness decays exponentially with time on current dose. Decay rate = 0.0217/week → half-life of 32 weeks. After 20 weeks at the same dose, effectiveness falls to 65%. Correlation with intake improves from r=-0.40 (raw blood level) to r=-0.43 (with tachyphylaxis). Resets partially on dose escalation.

**Intake prediction**: `calories = 2345 - 49 × effective_level`
- Zero drug: 2345 cal/day
- Fresh 12.5mg (effective ~17): 1504 cal/day
- 12.5mg after 20 weeks (effective ~11): 1800 cal/day

### Extraction not yet done
- **Workout PDFs** (`workout/`): ~200 gym session PDFs. Need OCR or PDF text extraction to structured sets/reps data.
- **Additional RMR**: More readings exportable from Cosmed Fitmate device.

### Join into analysis table
Build a single daily table joining: intake, weight, steps, sleep, medicine (as binary/dose columns), composition (interpolated). This is the input to all theories in theories.md.

### MFP API micronutrient enrichment
Use food names from intake_foods.csv to query MFP API for iron, potassium, vitamin D, saturated fat, added sugars. Match by name + calorie checksum per the approach in `/myfitnesspal/README.md`. The 1,092 macro-mismatch foods from sanity_check.py are highest priority (their existing data is unreliable).

### Water weight smoothing (glycogen + sodium) — DONE
Built `analysis/glycogen_smooth.py` with two independent corrections, both validated by the same criterion (nutrient→weight partial correlation drops to ~zero after correction):
- **Glycogen**: target-seeking model, tuned via grid search (G_max=600g, carb_ref=350g, rate_up=0.60, rate_down=0.45). Carb→weight partial r: 0.27 → ~0.
- **Sodium**: linear model, k=0.00030 lbs/mg (=136ml/g Na, literature 130-150). Sodium→weight partial r: 0.175 → -0.015.
- Both use 1-day lag (morning weight reflects previous day's intake), centered on median.
- **Combined: 14.3% global variance reduction**. Per-window: 40% (fasts), 62% (high-sodium-variance periods), 9-23% (other dense periods).
- Steps tested: r=0.014 partial correlation with weight change. No detectable signal at any timescale. Not included.
- Output: `analysis/smoothed_weight.csv`

### Weight interpolation + TDEE derivation — DONE
Built `analysis/interpolate_weight.py`. Derives TDEE by inverting the energy balance between each pair of weigh-ins, then simulates daily weight through gaps using that TDEE.
- 360 TDEE windows (min 7 days each): short windows get constant TDEE, long windows (>60d) get weight-dependent TDEE via Mifflin-St Jeor scale factor.
- Snaps to every observed weight (0.00 lbs max interpolation error at observations).
- Complete 5,429-day series: `analysis/daily_weight.csv` with interpolated weight, glycogen-corrected smoothed weight, and per-window TDEE.
- Validated against 3 indirect calorimetry measurements: derived TDEE/measured RMR ratios of 1.21, 1.14, 1.27 (expected 1.1-1.3 for TDEE/RMR).
- TDEE distribution: median 2158, P5=1792, P95=3207.

### Analysis — NEXT UP
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
