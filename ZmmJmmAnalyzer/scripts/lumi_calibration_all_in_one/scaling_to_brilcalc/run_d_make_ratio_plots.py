# run3_rate_vs_lumi_by_year.py
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep
plt.style.use(hep.style.CMS)

# -----------------------------
# Config you may tweak
# -----------------------------
INT_LUMI_CSV = "csvs/grouped_integrated_lumi/combined_lumi_Run3_integrated.csv"
DATA_CSV     = "csvs/combined_data_fit_Run3.csv"

OUTDIR = "plots/run3"
os.makedirs(OUTDIR, exist_ok=True)

LS_SECONDS = 23.4           # length of a luminosity section
BIN_LS     = 50             # your binning in LS
LUMI_REL_ERR = 0.013        # 1.3% lumi scale uncertainty

APPLY_PU_EFFICIENCY = True
PU_COEFF = -0.00322         # poly = 1 + PU_COEFF * pileup  (your coeffs [-0.00322, 1])

# Set a custom svis per year if desired; use None to auto-compute from data
# Example keeps your earlier choice 7.2 for 2023
CUSTOM_SVIS_BY_YEAR = {
    2022: 7.2,
    2023: 7.2,   # ← your custom svis
    2024: 7.2,
    2025: 7.2
}

# Quality cuts
REQUIRE_FIT_STATUS = "success"
MIN_NLUMIS = 40

# y-limit for the rate/lumi overlay plot
YLIM_RATE_PLOT = (0, 25)


# -----------------------------
# Utilities
# -----------------------------
def load_inputs(int_lumi_csv: str, data_csv: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    int_lumi = pd.read_csv(int_lumi_csv)[["bin_key", "integrated_lumi"]]
    df = pd.read_csv(data_csv)

    # Basic filtering
    df = df[df["fit_status"] == REQUIRE_FIT_STATUS]
    df = df[df["nlumis"] > MIN_NLUMIS]

    # Compute “per-second” and “Hz/nb” conversions
    df = df.copy()
    df["events"]       = df["events"] / (BIN_LS * LS_SECONDS)   # counts per second
    df["events_error"] = df["events_error"] / (BIN_LS * LS_SECONDS)
    df["lumi"]         = df["lumi"] / LS_SECONDS                # Hz/nb

    # Optional pileup efficiency correction: divide by (1 + a*PU)
    if APPLY_PU_EFFICIENCY:
        denom = 1.0 + PU_COEFF * df["pileup"].values
        # Safety: avoid tiny/negative denominators
        denom = np.clip(denom, 0.10, None)
        df["events"]       = df["events"] / denom
        df["events_error"] = df["events_error"] / denom

    # Useful diagnostic columns (unchanged behavior from your snippet)
    df["Ratio"]     = df["events"] / df["lumi"]
    df["Ratio_err"] = df["Ratio"] * np.sqrt((df["events_error"]/df["events"])**2 + LUMI_REL_ERR**2)

    # Merge integrated lumi
    df = df.merge(int_lumi, on="bin_key", how="left")
    return int_lumi, df


def compute_scaling(df_year: pd.DataFrame, custom_svis) -> tuple[float, float]:
    """
    Returns (scale_factor, svis_used). If custom_svis is not None, uses 1/custom_svis.
    Else computes scale factor so that sum(events_scaled) == sum(lumi), and reports svis_required.
    """
    sum_events = df_year["events"].sum()
    sum_lumi   = df_year["lumi"].sum()
    if sum_events <= 0 or sum_lumi <= 0:
        # Degenerate case: fall back safely
        scale_factor = 1.0
        svis_used = 1.0
        return scale_factor, svis_used

    if custom_svis is not None:
        scale_factor = 1.0 / float(custom_svis)
        svis_used = float(custom_svis)
    else:
        # scale so that events_scaled ~ lumi → svis_required = 1 / scale_factor
        scale_factor = sum_lumi / sum_events
        svis_used = 1.0 / scale_factor

    return scale_factor, svis_used


def sort_bins_like_input(df_year: pd.DataFrame) -> pd.DataFrame:
    """
    Keep a stable, readable order. If bin_key exists as string "run_lsbin",
    simple string sort is usually OK; otherwise preserve the original order.
    """
    if "bin_key" in df_year.columns:
        return df_year.sort_values("bin_key", kind="stable").reset_index(drop=True)
    return df_year.reset_index(drop=True)


def plot_rate_vs_lumi(df_year: pd.DataFrame, year: int, outdir: str) -> str:
    fig, ax1 = plt.subplots(figsize=(20, 4))

    x = np.arange(len(df_year))
    ax1.plot(x, df_year["events_scaled"], label=r"$J/\psi$ rate", marker=".", linestyle="none")
    ax1.plot(x, df_year["lumi"],          label="Ref. Luminosity", marker=".", linestyle="none")

    ax1.set_xlabel("Run_LS")
    ax1.set_ylabel("[Hz/nb]")
    ax1.grid(True)
    ax1.set_ylim(*YLIM_RATE_PLOT)

    ax1.xaxis.set_major_locator(ticker.MaxNLocator(50))
    visible_ticks = [int(t) for t in ax1.get_xticks() if 0 <= t < len(df_year)]
    ax1.set_xticks(visible_ticks)
    ax1.set_xticklabels(df_year.iloc[visible_ticks]["bin_key"], rotation=90, fontsize=7)

    # fig.legend(loc="best")
    plt.title(f"{year} $J/\\psi$ Rate and Ref Luminosity")
    fig.tight_layout()

    outpath = os.path.join(outdir, f"{year}_Jpsi_rate_vs_lumi.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return outpath


def plot_ratio_vs_intlumi(df_year: pd.DataFrame, year: int, outdir: str) -> str:
    # Only keep rows with valid integrated luminosity
    mask = df_year["integrated_lumi"].notna() & df_year["lumi"].ne(0)
    data = df_year.loc[mask].copy()
    if "final_ratio" not in data:
        data["final_ratio"] = data["events_scaled"] / data["lumi"]

    plt.figure(figsize=(20, 6))
    plt.plot(data["integrated_lumi"], data["final_ratio"], label="Ratio", marker=".", linestyle="none", color="purple")
    plt.ylim(0.7, 1.2)
    plt.xlabel(r"Integrated Luminosity [fb$^{-1}$]", fontsize=14)
    plt.ylabel(r"$J/\psi$ Rate / Ref Luminosity", fontsize=14)
    plt.grid(True)
    plt.title(f"Run 3 — {year}", fontsize=14)
    plt.tight_layout()

    outpath = os.path.join(outdir, f"{year}_ratio_vs_intlumi.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    return outpath


def process_year(df_all: pd.DataFrame, year: int) -> dict:
    df_year = df_all[df_all["year"] == year].copy()
    if df_year.empty:
        return {"year": year, "n": 0, "rate_plot": None, "ratio_plot": None, "svis": None}

    # Scale to match lumi (or use custom svis)
    scale_factor, svis_used = compute_scaling(df_year, CUSTOM_SVIS_BY_YEAR.get(year))
    df_year["events_scaled"] = df_year["events"] * scale_factor
    df_year["final_ratio"]   = df_year["events_scaled"] / df_year["lumi"]

    # Sort for nicer x-axis ordering
    df_year = sort_bins_like_input(df_year)

    # Make plots
    rate_png  = plot_rate_vs_lumi(df_year, year, OUTDIR)
    ratio_png = plot_ratio_vs_intlumi(df_year, year, OUTDIR)

    # Console summary
    print(f"[{year}] rows={len(df_year)} | svis_used={svis_used:.4f} | scale_factor={scale_factor:.4f}")

    return {"year": year, "n": len(df_year), "rate_plot": rate_png, "ratio_plot": ratio_png, "svis": svis_used}


def main():
    _, df_all = load_inputs(INT_LUMI_CSV, DATA_CSV)

    years = sorted(df_all["year"].dropna().unique().astype(int))
    summary = []
    for y in years:
        summary.append(process_year(df_all, y))

    # Optional: write a little summary CSV of svis per year
    pd.DataFrame(summary).to_csv(os.path.join(OUTDIR, "per_year_summary.csv"), index=False)
    print("Saved:", os.path.join(OUTDIR, "per_year_summary.csv"))


if __name__ == "__main__":
    main()
