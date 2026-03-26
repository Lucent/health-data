"""Evaluate manually annotated diet epochs against daily intake and latent states.

Questions:
Do the curated labels in intake/diet_epochs.csv separate meaningful numeric
regimes? Do the potato diets behave like a distinct intervention when pooled
across attempts?

Outputs:
    analysis/diet_epoch_summary.csv
    analysis/diet_epoch_family_summary.csv
    analysis/potato_epoch_window_summary.csv
    analysis/potato_epoch_contrast.csv
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
POTATO_WINDOW_DAYS = 14
MATCH_LOW_CAL = 1500
MATCH_HIGH_CAL = 2200
BINGE_THRESHOLD = 3000


def normalize_family(label: str) -> str:
    """Collapse numbered variants into a common family label."""
    return re.sub(r"_[0-9]+$", "", label)


def load_daily() -> pd.DataFrame:
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    daily = intake.merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    daily = daily.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily["binge"] = daily["calories"] >= BINGE_THRESHOLD
    daily["tdee_rmr_ratio"] = daily["tdee"] / daily["expected_rmr"]
    return daily.sort_values("date").reset_index(drop=True)


def load_epochs() -> pd.DataFrame:
    epochs = pd.read_csv(ROOT / "intake" / "diet_epochs.csv", parse_dates=["start", "end"])
    epochs["family"] = epochs["label"].map(normalize_family)
    epochs["days"] = (epochs["end"] - epochs["start"]).dt.days + 1
    return epochs


def summarize_slice(df: pd.DataFrame) -> dict[str, float]:
    return {
        "days": int(len(df)),
        "mean_calories": round(df["calories"].mean(), 1),
        "mean_carbs_g": round(df["carbs_g"].mean(), 1),
        "mean_fat_g": round(df["fat_g"].mean(), 1),
        "mean_protein_g": round(df["protein_g"].mean(), 1),
        "binge_rate": round(df["binge"].mean(), 4),
        "mean_tdee_rmr_ratio": round(df["tdee_rmr_ratio"].mean(), 4),
        "fat_start_lbs": round(df["fat_mass_lbs"].iloc[0], 2),
        "fat_end_lbs": round(df["fat_mass_lbs"].iloc[-1], 2),
        "fat_delta_lbs": round(df["fat_mass_lbs"].iloc[-1] - df["fat_mass_lbs"].iloc[0], 2),
        "fat_delta_per_30d_lbs": round(
            30.0 * (df["fat_mass_lbs"].iloc[-1] - df["fat_mass_lbs"].iloc[0]) / max(len(df) - 1, 1),
            2,
        ),
    }


def epoch_summary(daily: pd.DataFrame, epochs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, epoch in epochs.iterrows():
        mask = (daily["date"] >= epoch["start"]) & (daily["date"] <= epoch["end"])
        window = daily.loc[mask].copy()
        if window.empty:
            continue
        row = {
            "label": epoch["label"],
            "family": epoch["family"],
            "start": epoch["start"].strftime("%Y-%m-%d"),
            "end": epoch["end"].strftime("%Y-%m-%d"),
            "detail": epoch["detail"],
        }
        row.update(summarize_slice(window))
        rows.append(row)
    return pd.DataFrame(rows)


def family_summary(epoch_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    grouped = epoch_df.groupby("family", sort=True)
    for family, grp in grouped:
        total_days = grp["days"].sum()
        rows.append(
            {
                "family": family,
                "epochs": int(len(grp)),
                "days": int(total_days),
                "mean_calories": round((grp["mean_calories"] * grp["days"]).sum() / total_days, 1),
                "mean_carbs_g": round((grp["mean_carbs_g"] * grp["days"]).sum() / total_days, 1),
                "mean_fat_g": round((grp["mean_fat_g"] * grp["days"]).sum() / total_days, 1),
                "mean_protein_g": round((grp["mean_protein_g"] * grp["days"]).sum() / total_days, 1),
                "binge_rate": round((grp["binge_rate"] * grp["days"]).sum() / total_days, 4),
                "mean_tdee_rmr_ratio": round((grp["mean_tdee_rmr_ratio"] * grp["days"]).sum() / total_days, 4),
                "fat_delta_lbs": round(grp["fat_delta_lbs"].sum(), 2),
                "fat_delta_per_30d_lbs": round((grp["fat_delta_per_30d_lbs"] * grp["days"]).sum() / total_days, 2),
            }
        )
    return pd.DataFrame(rows)


def potato_window_summary(daily: pd.DataFrame, epochs: pd.DataFrame) -> pd.DataFrame:
    potato_epochs = epochs[epochs["family"] == "potato_diet"].copy()
    rows = []
    for _, epoch in potato_epochs.iterrows():
        windows = [
            ("pre", epoch["start"] - pd.Timedelta(days=POTATO_WINDOW_DAYS), epoch["start"] - pd.Timedelta(days=1)),
            ("epoch", epoch["start"], epoch["end"]),
            ("post", epoch["end"] + pd.Timedelta(days=1), epoch["end"] + pd.Timedelta(days=POTATO_WINDOW_DAYS)),
        ]
        for phase, start, end in windows:
            mask = (daily["date"] >= start) & (daily["date"] <= end) & (daily["effective_level"] == 0)
            window = daily.loc[mask].copy()
            if window.empty:
                continue
            row = {
                "label": epoch["label"],
                "phase": phase,
                "window_start": start.strftime("%Y-%m-%d"),
                "window_end": end.strftime("%Y-%m-%d"),
            }
            row.update(summarize_slice(window))
            rows.append(row)
    return pd.DataFrame(rows)


def potato_contrast(daily: pd.DataFrame, potato_windows: pd.DataFrame) -> pd.DataFrame:
    pre_tirz = daily[daily["effective_level"] == 0].copy()
    potato_dates = set()
    epoch_rows = potato_windows[potato_windows["phase"] == "epoch"]
    for _, row in epoch_rows.iterrows():
        window_dates = pd.date_range(row["window_start"], row["window_end"], freq="D")
        potato_dates.update(window_dates)

    matched = pre_tirz[
        (pre_tirz["calories"] >= MATCH_LOW_CAL)
        & (pre_tirz["calories"] <= MATCH_HIGH_CAL)
        & (~pre_tirz["date"].isin(potato_dates))
    ].copy()
    potato = pre_tirz[pre_tirz["date"].isin(potato_dates)].copy()

    rows = []
    for label, frame in [("potato_epoch_days", potato), ("matched_non_potato_days", matched)]:
        row = {
            "group": label,
            "days": int(len(frame)),
            "mean_calories": round(frame["calories"].mean(), 1),
            "mean_carbs_g": round(frame["carbs_g"].mean(), 1),
            "mean_fat_g": round(frame["fat_g"].mean(), 1),
            "mean_protein_g": round(frame["protein_g"].mean(), 1),
            "binge_rate": round(frame["binge"].mean(), 4),
            "mean_tdee_rmr_ratio": round(frame["tdee_rmr_ratio"].mean(), 4),
            "fat_start_lbs": None,
            "fat_end_lbs": None,
            "fat_delta_lbs": None,
            "fat_delta_per_30d_lbs": None,
        }
        rows.append(row)

    pooled = potato_windows.groupby("phase")[["mean_calories", "mean_protein_g", "mean_carbs_g", "binge_rate", "mean_tdee_rmr_ratio", "fat_delta_lbs"]].mean().reset_index()
    for _, row in pooled.iterrows():
        rows.append(
            {
                "group": f"potato_{row['phase']}_window_mean",
                "days": None,
                "mean_calories": round(row["mean_calories"], 1),
                "mean_carbs_g": round(row["mean_carbs_g"], 1),
                "mean_fat_g": None,
                "mean_protein_g": round(row["mean_protein_g"], 1),
                "binge_rate": round(row["binge_rate"], 4),
                "mean_tdee_rmr_ratio": round(row["mean_tdee_rmr_ratio"], 4),
                "fat_start_lbs": None,
                "fat_end_lbs": None,
                "fat_delta_lbs": round(row["fat_delta_lbs"], 2),
                "fat_delta_per_30d_lbs": None,
            }
        )
    return pd.DataFrame(rows)


def save_outputs(epoch_df: pd.DataFrame, family_df: pd.DataFrame, potato_df: pd.DataFrame, contrast_df: pd.DataFrame) -> None:
    epoch_df.to_csv(ROOT / "analysis" / "O_diet_epoch_summary.csv", index=False)
    family_df.to_csv(ROOT / "analysis" / "O_diet_epoch_family_summary.csv", index=False)
    potato_df.to_csv(ROOT / "analysis" / "O_potato_epoch_window_summary.csv", index=False)
    contrast_df.to_csv(ROOT / "analysis" / "O_potato_epoch_contrast.csv", index=False)


def print_report(epoch_df: pd.DataFrame, family_df: pd.DataFrame, potato_df: pd.DataFrame, contrast_df: pd.DataFrame) -> None:
    print("\n=== Diet Epoch Analysis ===")
    print(f"Epochs analyzed: {len(epoch_df)}")

    print("\nLowest-calorie named epochs:")
    cols = ["label", "days", "mean_calories", "mean_protein_g", "binge_rate", "mean_tdee_rmr_ratio", "fat_delta_lbs"]
    print(epoch_df.sort_values("mean_calories").head(6)[cols].to_string(index=False))

    print("\nPooled families with multiple epochs:")
    fam_cols = ["family", "epochs", "days", "mean_calories", "mean_protein_g", "binge_rate", "mean_tdee_rmr_ratio", "fat_delta_lbs"]
    multi = family_df[family_df["epochs"] > 1].sort_values("days", ascending=False)
    if multi.empty:
        print("  none")
    else:
        print(multi[fam_cols].to_string(index=False))

    print("\nPotato before/during/after means:")
    phases = potato_df.groupby("phase")[["days", "mean_calories", "mean_protein_g", "binge_rate", "mean_tdee_rmr_ratio", "fat_delta_lbs"]].mean()
    print(phases.round(4).to_string())

    print("\nPotato contrast:")
    contrast_cols = ["group", "mean_calories", "mean_protein_g", "binge_rate", "mean_tdee_rmr_ratio", "fat_delta_lbs"]
    print(contrast_df[contrast_cols].to_string(index=False))


def main() -> None:
    daily = load_daily()
    epochs = load_epochs()
    epoch_df = epoch_summary(daily, epochs)
    family_df = family_summary(epoch_df)
    potato_df = potato_window_summary(daily, epochs)
    contrast_df = potato_contrast(daily, potato_df)
    save_outputs(epoch_df, family_df, potato_df, contrast_df)
    print_report(epoch_df, family_df, potato_df, contrast_df)


if __name__ == "__main__":
    main()
