import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Config
# -----------------------------
INT_LUMI_CSV = "csvs/grouped_integrated_lumi/combined_lumi_Run3_integrated.csv"
DATA_CSV     = "csvs/event_data/combined_data_fit_Run3.csv"
OUTDIR       = "plots/ML_fits_summary"
os.makedirs(OUTDIR, exist_ok=True)

LS_SECONDS   = 23.4
BIN_LS       = 50
LUMI_REL_ERR = 0.013

APPLY_PU_EFFICIENCY = True
PU_COEFF = [-0.00322, 1]   # coefficients for pileup efficiency correction

CUSTOM_SVIS_BY_YEAR = {
    2022: 7.2,
    2023: 7.2,
    2024: 7.2,
    2025: 7.2
}

OUTPUT_CSV   = "csvs/final_fit_Run3.csv"

# -----------------------------
# Functions
# -----------------------------
def load_inputs(int_lumi_csv: str, data_csv: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load luminosity and event data with initial filtering."""
    int_lumi = pd.read_csv(int_lumi_csv)
    int_lumi = int_lumi[int_lumi["delivered"] > 0][["bin_key", "integrated_lumi"]]

    df = pd.read_csv(data_csv)
    # Apply cuts
    df = df[
        (df["fit_status"] == "success")
        & (df["lumi"] > 0)
        & (df["nlumis"] > 40)
        & (df["pileup"].between(0, 80))
    ]
    return int_lumi, df


def plot_histograms(df: pd.DataFrame, outdir: str) -> None:
    """Save histograms for all columns in df."""
    for col in df.columns:
        plt.figure(figsize=(8, 6))
        plt.hist(df[col], bins=50, alpha=0.7, color="blue", edgecolor="black")
        plt.title(col, fontsize=16)
        plt.xlabel(col, fontsize=14)
        plt.ylabel("No of 50 LS bins", fontsize=14)
        plt.grid(axis="y", alpha=0.75)
        plt.tight_layout()
        plt.savefig(f"{outdir}/{col}_histogram.png")
        plt.close()


def apply_conversions(df: pd.DataFrame) -> pd.DataFrame:
    """Convert events and lumi to per-second and Hz/nb, apply PU efficiency correction if requested."""
    df = df[["bin_key", "run", "lumiblock", "events", "events_error", "lumi", "pileup", "nPV", "nPV_err", "year", "era"]].copy()

    # Unit conversions
    df["events"]       = df["events"] / (BIN_LS * LS_SECONDS)
    df["events_error"] = df["events_error"] / (BIN_LS * LS_SECONDS)
    df["lumi"]         = df["lumi"] / LS_SECONDS

    # PU efficiency correction
    if APPLY_PU_EFFICIENCY:
        print("Applying pileup efficiency correction")
        denom = np.polyval(PU_COEFF, df["pileup"].values)
        # denom = np.clip(denom, 0.10, None)  # avoid tiny/negative
        df["events"]       = df["events"] / denom
        df["events_error"] = df["events_error"] / denom

    return df


def compute_scaling(df, custom_svis) -> tuple[float, float]:
    """Compute scale factor and sigvis used for a given year's data."""
    sum_events = df["events"].sum()
    sum_lumi   = df["lumi"].sum()

    if custom_svis is not None:
        return 1.0 / custom_svis, custom_svis
    scale_factor = sum_lumi / sum_events
    return scale_factor, 1.0 / scale_factor


def apply_scaling(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply scaling by year and return updated df + summary."""
    df = df.copy()
    df["events_scaled"] = np.nan
    df["final_ratio"]   = np.nan

    summary = []
    for year in sorted(df["year"].dropna().unique().astype(int)):
        mask = df["year"] == year
        df_year = df.loc[mask]

        scale_factor, svis_used = compute_scaling(df_year, CUSTOM_SVIS_BY_YEAR.get(year))
        df.loc[mask, "events_scaled"] = df_year["events"] * scale_factor
        df.loc[mask, "final_ratio"]   = df.loc[mask, "events_scaled"] / df_year["lumi"]

        print(f"[{year}] rows={len(df_year)} | svis_used={svis_used:.4f} | scale_factor={scale_factor:.4f}")

    return df


def main():
    # Load inputs
    int_lumi, event_df = load_inputs(INT_LUMI_CSV, DATA_CSV)
    print(event_df.columns)

    # Plot histograms for diagnostics
    # plot_histograms(event_df, OUTDIR)

    # Apply conversions + PU efficiency correction
    event_df = apply_conversions(event_df)

    # Merge lumi
    event_df = event_df.merge(int_lumi, on="bin_key", how="left")

    # Apply scaling by year
    event_df = apply_scaling(event_df)

    # Save outputs
    event_df.to_csv(OUTPUT_CSV, index=False)

    print(f"Saved prepared data → {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
