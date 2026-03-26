# Exercise and step data pipeline

```
intake/merge.py → steps-sleep/mfp_exercises.csv (all MFP exercise entries)
                                ↓
            ┌───────────────────┼───────────────────┐
            ↓                   ↓                   ↓
   steps-sleep/            workout/            steps-sleep/
   merge_exercises.py      merge.py            merge_steps.py
            ↓                   ↓                   ↓
     exercises.csv          strength.csv         steps.csv
  (samsung + mfp runs)   (pdf + chloe + mfp)  (samsung + mfp/cal hospital)
```

## Sources

### Steps
- **steps_samsung.csv** — Samsung Health daily step counts (2014-04+)
- **mfp_exercises.csv** — Hospital (5500 steps), walking (mph x time), running (mph x time), treadmill
- **steps_calendar.csv** — 16 hospital shifts from calendar that MFP missed

### Cardio exercises
- **exercises_samsung.csv** — auto-detected walk/run/bike/hike sessions (2014-04+)
- **mfp_exercises.csv** — 7 pre-Samsung running entries (Apr 2013)

### Strength training
- **workout/*.pdf** — ActivTrax YMCA sessions (2018-01 to 2025-10, 284 dates)
- **workout/Chloe Workout.xlsx** — personal trainer sessions (2016-09 to 2017-01, 37 dates)
- **mfp_exercises.csv** — "Circuit training, general" fills gaps (71 unique dates, mostly 2017)

## Deduplication

Steps: Samsung takes priority. Backfill only covers dates with no Samsung data.

Exercises: Samsung is canonical for all dates it covers. MFP adds only pre-Samsung runs.

Strength: PDF > Chloe > MFP Circuit on date overlap. MFP dates where the next day has a PDF are off-by-one duplicates (after-midnight logging) and are excluded.
