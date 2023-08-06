I accept any question, no matter how personal, in service of our goal. You will not offend me.

# [Theories to test](theories.md)

# [Arguments for quality](quality.md)

# [Health background](background.md)

# Data Directories
## /composition
* All data in XLSX. Scans are backups.
* BOD POD measurements (2011, 2016) are air displacement
* InBody is bioelectric impedance

## /intake
* XPS non-vector and HTML saves of MyFitnessPal data. Last 2 years can be queried directly for full nutrition data.
* Identical data in multiple formats to assist in import or verify accuracy
* MyFitnessPal API holds full data on entries >2 years old, but exact matches must be ensured due to large amounts of garbage entries in their database.
* Pretty limited diet variety. You'll see the same foods over and over and can create a database of them.

## /RMR
* All data in XLSX. Scans are backups.
* Any measure before noon is fasted
* Data from indirect calorimetry, Cosmed Fitmate

## /steps-sleep (also contains heart rate, detected exercise of > 10 min walking)
* Samsung Health export is ugly but some CSV contain JavaScript milliseconds since epoch timestamps
* Very little sunlight exposure
* Normal sleep cycle is 3am to 11am. Samsung Health export will show deviations with wake/sleep timestamps (incomplete).
* Steps over 2000/hr under 2 hrs total before sunset most likely taken outdoors, sunlight determinable from historical hourly weather
* Steps under 2000/hr over 2 hrs total before sunset most likely taken indoors
* Steps over 2000/hr under 2 hrs total after sunset most likely taken indoors

## /weight
* Any measure before noon is fasted
* XLSX weight data is reliable
* All % composition data is **UNRELIABLE**. Use /composition

## /workout
* All gym visits with sets/reps
* List of 30 minute sessions by date in Chloe Workout.xlsx
