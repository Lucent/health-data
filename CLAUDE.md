# Project context for AI assistants

Read these files in order:
1. [README.md](README.md) — data directories, extracted CSVs, how to regenerate
2. [background.md](background.md) — health history, genetics, substances, locations, weighing protocol
3. [quality.md](quality.md) — arguments for and against data quality
4. [theories.md](theories.md) — testable hypotheses organized by data requirements
5. [ROADMAP.md](ROADMAP.md) — completed work, next steps, known data quality notes

## Key principles
- **The data is complete.** Every day is fully logged. A 10-calorie day is a real fasting day, not incomplete logging. Never assume missing entries.
- **Accuracy is unusually high.** Food scale always used. ~5% consistent undercount from uncounted snacking acknowledged. Incentive structure rewards precision over undercounting.
- **CICO is effect, not cause.** The working theory is a meandering set point controlled by an unknown homeostat. Simple calorie correlation has been tried and drifts.

## Extraction pipeline
All extractors are idempotent. Re-run to regenerate CSVs from raw data.
```
python intake/merge.py              → intake/intake_foods.csv + intake/intake_daily.csv
python intake/verify.py             — full integrity check (run after every merge)
python weight/extract.py            → weight/weight.csv
python steps-sleep/extract.py       → steps-sleep/steps.csv + steps-sleep/sleep.csv
python intake/sanity_check.py       — Atwater factor validation
python intake/compare_extractors.py — cross-validates OXPS vs HTML extractors
```
