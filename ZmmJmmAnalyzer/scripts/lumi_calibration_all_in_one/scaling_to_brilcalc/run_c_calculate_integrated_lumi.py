import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

"""Use the brilcalc CSV files to calculate integrated luminosity in 50 lumiblock bins.
Generates plots of source distribution and combined CSV files.
"""

def read_lumi_data(lumi_file: Path) -> pd.DataFrame:
    """Read luminosity data from brilcalc CSV."""
    df = pd.read_csv(lumi_file, skiprows=1, engine="python", skipfooter=5)
    df["time"] = pd.to_datetime(df["time"], format="%m/%d/%y %H:%M:%S")
    df[["run", "fill"]] = df["#run:fill"].str.split(":", expand=True)
    df = df.drop(columns=["#run:fill"])
    df = df[[ "run", "fill", "time", "ls", "beamstatus", "recorded(/ub)", "delivered(/ub)", "avgpu", "source"]]
    df["run"] = df["run"].astype(int)
    df["ls"] = df["ls"].str.split(":").str[0].astype(int)
    return df


def process_lumi_data(df: pd.DataFrame) -> pd.DataFrame:
    """Filter and group lumi data by 50 lumiblock bins."""
    df = df[df["beamstatus"] == "STABLE BEAMS"].copy()
    df["lumiblock"] = (df["ls"] // 50) * 50
    df["delivered"] = df["delivered(/ub)"] * 1e-9  # Convert to /fb
    df["bin_key"] = df["run"].astype(str) + "_" + df["lumiblock"].astype(str)

    grouped = (
        df.groupby("bin_key", as_index=False)
        .agg(run=("run", "first"),
             lumiblock=("lumiblock", "first"),
             delivered=("delivered", "sum"))
    )
    return grouped.sort_values(by=["run", "lumiblock"])


def plot_source_distribution(df: pd.DataFrame, year: int, outdir: Path):
    """Plot source distribution for a given year."""
    outdir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(5, 4))
    df["source"].value_counts().plot(kind="bar")
    plt.yscale("log")
    plt.xlabel("Source")
    plt.ylabel("Number of Lumisections")
    plt.title(f"Lumi Source Distribution - {year}")
    plt.grid(axis="y")
    plt.tight_layout()
    plt.savefig(outdir / f"source_distribution_{year}.png")
    plt.close()


def combine_csv(files, outfile: Path):
    """Combine multiple CSV files into one sorted output."""
    combined = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    combined.sort_values(by=["run", "lumiblock"], inplace=True)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(outfile, index=False)
    print(f"Combined CSV saved to: {outfile}")


def main():
    base_path = Path("/home/nimmitha/LPCfiles/run3/jmm_counting/cmssw/ZmmJmmAnalyzer/preselection/parkingDoubleMuonLowMass")
    output_dir = Path("csvs/grouped_integrated_lumi")
    plot_dir = Path("plots/lumi_source")

    year_files = {
        2022: "lumi_normtagPHYSICS_2022.csv",
        2023: "lumi_normtagPHYSICS_2023.csv",
        2024: "lumi_normtagBRIL_2024.csv",
        2025: "lumi_normtagBRIL_2025.csv",
    }

    # Process each year
    for year, fname in year_files.items():
        df = read_lumi_data(base_path / fname)
        grouped = process_lumi_data(df)

        out_csv = output_dir / f"grouped_{fname}"
        grouped.to_csv(out_csv, index=False)

        plot_source_distribution(df, year, plot_dir)

    # Combine all grouped CSVs
    grouped_files = sorted(output_dir.glob("grouped_lumi_*.csv"))
    combined_file = output_dir / "combined_lumi_Run3.csv"
    combine_csv(grouped_files, combined_file)

    # Add cumulative lumi
    df = pd.read_csv(combined_file)
    df.sort_values(by=["run", "lumiblock"], inplace=True)
    df["integrated_lumi"] = df["delivered"].cumsum()
    df.to_csv(output_dir / "combined_lumi_Run3_integrated.csv", index=False)


if __name__ == "__main__":
    main()
