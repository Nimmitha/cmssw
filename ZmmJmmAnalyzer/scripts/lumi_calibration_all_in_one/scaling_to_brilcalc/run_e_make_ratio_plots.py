import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep
plt.style.use(hep.style.CMS)

# -----------------------------
# Config
# -----------------------------
PREPARED_CSV = "csvs/final_fit_Run3.csv"   # from Script 1
OUTDIR = "plots/run3"
os.makedirs(OUTDIR, exist_ok=True)

YLIM_RATE_PLOT = (0, 25)


def plot_rate_vs_lumi(df_year: pd.DataFrame, year: int):
    """Plot scaled event rate vs lumi for one year."""
    fig, ax1 = plt.subplots(figsize=(20, 4))

    x = np.arange(len(df_year))
    ax1.plot(x, df_year["events_scaled"], label=r"$J/\psi$ rate", marker=".", linestyle="none")
    ax1.plot(x, df_year["lumi"],          label="Ref. Luminosity", marker=".", linestyle="none")

    ax1.set_xlabel("Run_LS", fontsize=14)
    ax1.set_ylabel("[Hz/nb]", fontsize=14)
    ax1.grid(True)
    ax1.set_ylim(*YLIM_RATE_PLOT)

    ax1.xaxis.set_major_locator(ticker.MaxNLocator(50))
    visible_ticks = [int(t) for t in ax1.get_xticks() if 0 <= t < len(df_year)]
    ax1.set_xticks(visible_ticks)
    ax1.set_xticklabels(df_year.iloc[visible_ticks]["bin_key"], rotation=90, fontsize=7)

    plt.title(f"{year} $J/\\psi$ Rate and Ref Luminosity", fontsize=14)
    fig.tight_layout()

    outpath = os.path.join(OUTDIR, f"{year}_Jpsi_rate_vs_lumi.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_ratio_vs_intlumi(df_year: pd.DataFrame, year: int):
    """Plot final ratio vs integrated lumi for one year."""
    mask = df_year["integrated_lumi"].notna() & df_year["lumi"].ne(0)
    data = df_year.loc[mask].copy()

    if "final_ratio" not in data:
        data["final_ratio"] = data["events_scaled"] / data["lumi"]

    plt.figure(figsize=(20, 6))
    plt.plot(
        data["integrated_lumi"], data["final_ratio"],
        label="Ratio", marker=".", linestyle="none", color="purple"
    )
    plt.ylim(0.7, 1.2)
    plt.xlabel(r"Integrated Luminosity [fb$^{-1}$]", fontsize=14)
    plt.ylabel(r"$J/\psi$ Rate / Ref Luminosity", fontsize=14)
    plt.grid(True)
    plt.title(f"Run 3 — {year}", fontsize=14)
    plt.tight_layout()

    outpath = os.path.join(OUTDIR, f"{year}_ratio_vs_intlumi.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    # Load prepared file
    df_all = pd.read_csv(PREPARED_CSV)

    # Loop over years
    years = sorted(df_all["year"].dropna().unique().astype(int))
    df_all.sort_values(by=["run", "lumiblock"], inplace=True)

    for year in years:
        df_year = df_all[df_all["year"] == year].copy()
        if df_year.empty:
            continue
        plot_rate_vs_lumi(df_year, year)
        plot_ratio_vs_intlumi(df_year, year)

        print(f"[{year}] Processed {len(df_year)} rows.")


if __name__ == "__main__":
    main()
