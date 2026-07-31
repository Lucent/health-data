# Prospective tracking

These templates fill the gaps that the retrospective repository cannot solve. Keep burden low: MyFitnessPal already supplies calories, protein, fiber, and food items.

## Core daily fields

Use `daily.csv` for:

- morning weight under the existing protocol;
- morning and evening hunger, 0–10;
- GI symptoms, 0–10;
- travel, order-in, and planned menu block;
- whether above-plan calorie momentum was noticed before the final meal and whether the ordinary-default reset was used;
- whether an injection occurred, using the prescribed schedule;
- a short note for illness or unusual events.

Do not copy calories/protein/fiber manually if the MFP export will remain canonical.
The reset fields test a low-risk implementation cue; they must not trigger
fasting, punitive exercise, or medication changes.

Treat the weight series as stale after 14 days and show a neutral
resume-measurement prompt after 30 days. Do this in reporting code rather than
adding another daily manual field.

## Device-derived sleep and steps

Continue importing Samsung sleep and steps rather than copying them into the
daily template. For prospective analyses:

- define main sleep as the longest plausible 2–14 hour episode assigned to wake date;
- flag sleep under 6 hours rather than deliberately restricting sleep;
- flag daily steps under 100 or a session-step/daily-total mismatch as possible non-wear;
- use exact calendar dates for future outcomes.

The reliable historical sleep regime already averages about eight hours, so
maintenance and apnea-symptom screening are higher value than forcing an
earlier bedtime.

## Strength sessions

Use `strength_sessions.csv` for exercises, sets, reps, load, and effort. The existing `workout/strength.csv` stores only a date and an assumed 30 minutes; it cannot test progressive overload or function.

## Weekly health checks

Use `weekly.csv` for waist, home blood pressure if already available/clinically appropriate, resting heart rate, grip strength if a consistent dynamometer is available, and a standardized chair-stand time.

Use the same equipment, conditions, and measurement protocol. Blank is better than guessed.

## Medication boundary

Medication changes belong to the prescriber and pharmacist. These files observe a prescribed plan; they are not a protocol for changing dose, frequency, or pen handling.

## Data refresh

Before analyzing anything after February 2026, refresh:

1. MFP export from 2026-03-01 onward.
2. Scale history after 2026-03-22.
3. Medication/injection log after 2026-04-10.
4. Same-device composition after 2026-02-28, if any.
