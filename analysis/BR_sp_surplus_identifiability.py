#!/usr/bin/env python3
"""BR. Identifiability of SP half-life vs surplus lookback window.

Question:
  Does the ~45-50 day set-point half-life remain preferred when the surplus
  lookback window is also allowed to vary, or are these two smoothing choices
  partially interchangeable?

Method:
  - Use 2014-01-01 to 2023-12-31 pre-tirzepatide data only.
  - Compute symmetric EMA set point over fat mass for a grid of half-lives.
  - Compute rolling mean surplus for a grid of lookback windows.
  - Correlate SP distance (SP - FM) with rolling surplus.

Output:
  - analysis/BR_sp_surplus_identifiability.csv
  - analysis/BR_sp_surplus_identifiability.png
"""

from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent.parent
OUT_CSV = ROOT / "analysis" / "BR_sp_surplus_identifiability.csv"
OUT_PNG = ROOT / "analysis" / "BR_sp_surplus_identifiability.png"

DATE_START = "2014-01-01"
DATE_END = "2024-01-01"
HALF_LIVES = list(range(15, 121, 5))
WINDOWS = list(range(30, 181, 5))


def ema(series: pd.Series, half_life: int) -> pd.Series:
    alpha = 1 - np.exp(-np.log(2) / half_life)
    return series.ewm(alpha=alpha, min_periods=30).mean()


def load_pre_drug_daily() -> pd.DataFrame:
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]],
        on="date",
        how="inner",
    )
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df["surplus"] = df["calories"] - df["tdee"]

    mask = (
        (df["effective_level"] == 0)
        & (df["date"] >= DATE_START)
        & (df["date"] < DATE_END)
    )
    return df.loc[mask].copy().reset_index(drop=True)


def evaluate_grid(df: pd.DataFrame) -> pd.DataFrame:
    fm = df["fat_mass_lbs"]
    rows = []

    for hl in HALF_LIVES:
        sp = ema(fm, hl)
        dist = sp - fm
        for win in WINDOWS:
            rolled = df["surplus"].rolling(win, min_periods=win).mean()
            valid = dist.notna() & rolled.notna()
            if valid.sum() < 300:
                continue
            r = np.corrcoef(dist[valid], rolled[valid])[0, 1]
            rows.append(
                {
                    "half_life_days": hl,
                    "surplus_window_days": win,
                    "r": r,
                    "abs_r": abs(r),
                    "n": int(valid.sum()),
                }
            )

    return pd.DataFrame(rows).sort_values("abs_r", ascending=False).reset_index(drop=True)


def make_plot(grid: pd.DataFrame) -> None:
    pivot = grid.pivot(
        index="half_life_days",
        columns="surplus_window_days",
        values="r",
    ).sort_index().sort_index(axis=1)

    x_vals = list(pivot.columns)
    y_vals = list(pivot.index)
    z = pivot.values
    best = grid.iloc[0]

    ridge = (
        grid.sort_values(["surplus_window_days", "abs_r"], ascending=[True, False])
        .drop_duplicates("surplus_window_days")
        .sort_values("surplus_window_days")
    )

    cell = 18
    left = 100
    top = 55
    right = 170
    bottom = 120
    width = left + len(x_vals) * cell + right
    height = top + len(y_vals) * cell + bottom

    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()

    z_min = float(np.nanmin(z))
    z_max = float(np.nanmax(z))

    def blend(c1, c2, t):
        return tuple(int(round(a + (b - a) * t)) for a, b in zip(c1, c2))

    def color_for(value):
        t = 0.0 if z_max == z_min else (value - z_min) / (z_max - z_min)
        return blend((35, 87, 164), (255, 245, 157), t)

    # Title and labels
    draw.text((left, 15), "SP Half-Life vs Surplus Lookback: Identifiability Ridge", fill="black", font=font)
    draw.text((left, height - 95), "Surplus lookback window (days)", fill="black", font=font)
    draw.text((15, top - 25), "SP half-life", fill="black", font=font)

    # Heatmap cells
    for yi, hl in enumerate(y_vals):
        for xi, win in enumerate(x_vals):
            value = pivot.loc[hl, win]
            x0 = left + xi * cell
            y0 = top + yi * cell
            draw.rectangle([x0, y0, x0 + cell, y0 + cell], fill=color_for(value))

    # Grid and axes
    draw.rectangle(
        [left, top, left + len(x_vals) * cell, top + len(y_vals) * cell],
        outline="black",
        width=1,
    )
    for xi, win in enumerate(x_vals):
        if xi % 3 == 0:
            x = left + xi * cell
            draw.line([x, top, x, top + len(y_vals) * cell], fill=(220, 220, 220))
            draw.text((x - 6, top + len(y_vals) * cell + 6), str(win), fill="black", font=font)
    for yi, hl in enumerate(y_vals):
        if yi % 2 == 0:
            y = top + yi * cell
            draw.line([left, y, left + len(x_vals) * cell, y], fill=(220, 220, 220))
            draw.text((left - 30, y - 4), str(hl), fill="black", font=font)

    # Ridge line
    ridge_points = []
    for _, row in ridge.iterrows():
        xi = x_vals.index(int(row["surplus_window_days"]))
        yi = y_vals.index(int(row["half_life_days"]))
        ridge_points.append((left + xi * cell + cell // 2, top + yi * cell + cell // 2))
    if len(ridge_points) > 1:
        draw.line(ridge_points, fill="white", width=2)

    # Best point
    best_xi = x_vals.index(int(best["surplus_window_days"]))
    best_yi = y_vals.index(int(best["half_life_days"]))
    bx = left + best_xi * cell + cell // 2
    by = top + best_yi * cell + cell // 2
    draw.ellipse([bx - 4, by - 4, bx + 4, by + 4], fill="red", outline="black")

    # Legend / note
    legend_x = left + len(x_vals) * cell + 20
    legend_y = top
    draw.text((legend_x, legend_y), "Best ridge", fill="black", font=font)
    draw.line([legend_x, legend_y + 16, legend_x + 30, legend_y + 16], fill="white", width=2)
    draw.text((legend_x, legend_y + 32), "Best overall", fill="black", font=font)
    draw.ellipse([legend_x, legend_y + 48, legend_x + 8, legend_y + 56], fill="red", outline="black")
    draw.text(
        (legend_x, legend_y + 70),
        f"HL={int(best['half_life_days'])}d,\nwin={int(best['surplus_window_days'])}d,\nr={best['r']:+.3f}",
        fill="black",
        font=font,
    )

    bar_x0 = legend_x
    bar_y0 = legend_y + 125
    bar_h = 120
    bar_w = 18
    for i in range(bar_h):
        t = 1 - i / max(1, bar_h - 1)
        value = z_min + t * (z_max - z_min)
        draw.rectangle(
            [bar_x0, bar_y0 + i, bar_x0 + bar_w, bar_y0 + i + 1],
            fill=color_for(value),
        )
    draw.rectangle([bar_x0, bar_y0, bar_x0 + bar_w, bar_y0 + bar_h], outline="black")
    draw.text((bar_x0 + 28, bar_y0 - 4), f"{z_max:+.3f}", fill="black", font=font)
    draw.text((bar_x0 + 28, bar_y0 + bar_h - 8), f"{z_min:+.3f}", fill="black", font=font)
    draw.text((legend_x, bar_y0 + bar_h + 10), "Correlation r", fill="black", font=font)

    note = (
        "Broad diagonal ridge: shorter SP half-lives pair with shorter surplus windows; "
        "longer half-lives pair with longer windows. The surplus regression supports a slow "
        "appetite-pressure timescale, but does not uniquely identify 45d by itself."
    )
    draw.multiline_text((left, height - 65), note, fill="black", font=font, spacing=3)

    img.save(OUT_PNG)


def main() -> None:
    df = load_pre_drug_daily()
    grid = evaluate_grid(df)
    grid.to_csv(OUT_CSV, index=False)
    make_plot(grid)

    print("=" * 88)
    print("SP HALF-LIFE vs SURPLUS LOOKBACK IDENTIFIABILITY")
    print("=" * 88)
    print(f"Period: {DATE_START} to {pd.Timestamp(DATE_END) - pd.Timedelta(days=1):%Y-%m-%d}, pre-tirzepatide")
    print(f"Rows: {len(df)} days")
    print("\nTop 10 combinations by |r|:")
    print(
        grid.head(10).to_string(
            index=False,
            formatters={
                "r": lambda x: f"{x:+.4f}",
                "abs_r": lambda x: f"{x:.4f}",
            },
        )
    )

    by_window = (
        grid.sort_values(["surplus_window_days", "abs_r"], ascending=[True, False])
        .drop_duplicates("surplus_window_days")
        .sort_values("surplus_window_days")
    )
    print("\nBest half-life at each surplus window:")
    print(
        by_window.to_string(
            index=False,
            formatters={
                "r": lambda x: f"{x:+.4f}",
                "abs_r": lambda x: f"{x:.4f}",
            },
        )
    )
    print(f"\nArtifacts: {OUT_CSV}, {OUT_PNG}")


if __name__ == "__main__":
    main()
