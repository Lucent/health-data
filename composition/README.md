# Composition Data Notes

`extract.py` now builds `composition.csv` from two sources:

- `InBody-*.csv`: imported exports from the InBody source system
- `Body composition.xlsx`: static BOD POD rows only

The workbook should retain only the pre-InBody BOD POD rows. Historical PDFs remain the raw archive for older scans.

## Fields Missing From The New InBody Export

The imported InBody CSV does not preserve all historical fields for the early `2017-03-17` through `2017-10-18` InBody 370 measurements.

These fields are blank in the new export for that period even though they exist in the older local records:

- `arm_muscle_r_lbs`
- `arm_muscle_l_lbs`
- `trunk_muscle_lbs`
- `leg_muscle_r_lbs`
- `leg_muscle_l_lbs`

The new export also lacks these early-era fields for those rows:

- intracellular water
- extracellular water
- ECW ratio
- visceral fat level
- segmental fat mass fields
- InBody score

So:

- Keep `2017 InBody composition.pdf` as the raw source archive for the `2017-03-17` to `2017-10-18` scans.
- Keep the two BOD POD PDFs for `2011-04-28` and `2016-01-23`.
- The deleted dated PNGs are redundant with the imported InBody CSV and do not need to be restored.

## Output Notes

`composition.csv` keeps the existing core columns used by downstream analysis and adds:

- `measured_at`: ISO timestamp from the InBody export, or midnight for BOD POD rows
- `source`: input filename used for the kept row
- `device`: exported device code / `BOD POD`

If multiple InBody measurements exist on the same day, `extract.py` keeps one row per date and prefers the richer row, then the later timestamp.
