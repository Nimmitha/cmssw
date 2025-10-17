import pandas as pd
import numpy as np
import uproot as up
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path

# ---------- CONFIGURATION ----------
# ROOT_FILE_PATH = '/home/nimmitha/LPCfiles/run3/jmm_counting/cmssw/ZmmJmmAnalyzer/preselection/parkingDoubleMuonLowMass/2023/emit/PDMLM_mm_2023D_Run370293_370579_emit_v1.root'
# ROOT_FILE_PATH = '/home/nimmitha/LPCfiles/run3/jmm_counting/cmssw/ZmmJmmAnalyzer/preselection/parkingDoubleMuonLowMass/2023/emit/PDMLM_mm_2023D1_emit_v3.root'
ROOT_FILE_PATH = '/home/nimmitha/LPCfiles/run3/jmm_counting/cmssw/ZmmJmmAnalyzer/preselection/parkingDoubleMuonLowMass/2024/emit/PDMLM_mm_2024F1_emit_v3.root'
BRANCHES_TO_READ = ["eventTime", "B_J1_mass", "B_Mu1_pt", "B_Mu2_pt", "B_J1_pt"]
SCAN_DIR = Path('output/scan_info')
OUTPUT_PLOT_DIR = Path('output/plots')
OUTPUT_CSV_DIR = Path('output/csvs')
TIME_BIN = '1s'
SAVE_PLOTS = True
MUON_PT_CUT = 5000  # example threshold, change as needed
MASS_CUT = (2.7, 3.5)  # example mass cut range

OUTPUT_PLOT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_CSV_DIR.mkdir(parents=True, exist_ok=True)

# ---------- FUNCTIONS ----------
def load_event_data(root_file_path, branches):
    """Load and preprocess ROOT data."""
    tree = up.open(root_file_path)["ntuple"]
    print(f"[INFO] available branches: {tree.keys()}")
    print(f"[INFO] loading branches: {branches}")
    arrays = tree.arrays(branches, library="np")
    df = pd.DataFrame(arrays)
    df['Time'] = pd.to_datetime(df['eventTime'], unit='ms')
    return df

def apply_cuts(df, pt_threshold, mass_range):
    """Apply muon pt selection."""
    df = df[(df['B_Mu1_pt'] > pt_threshold) & (df['B_Mu2_pt'] > pt_threshold)]
    # df = df[(df['B_J1_pt'] > pt_threshold*2)]
    df = df[(df['B_J1_mass'] > mass_range[0]) & (df['B_J1_mass'] < mass_range[1])]
    return df

def load_scan_file(csv_file):
    """Load scan step info from CSV."""
    try:
        df = pd.read_csv(csv_file)
        df['tStart'] = pd.to_datetime(df['tStart'], unit='s')
        df['tStop'] = pd.to_datetime(df['tStop'], unit='s')
        return df
    except Exception as e:
        print(f"[ERROR] Failed to load {csv_file.name}: {e}")
        return None

def count_events_per_step(df, scan_df):
    """Count passing events and collect B_J1_mass values per scan step."""
    event_counts = []
    bj1_mass_lists = []

    for _, row in scan_df.iterrows():
        t0 = row['tStart']
        t1 = row['tStop']
        mask = (df['Time'] >= t0) & (df['Time'] <= t1)
        selected = df.loc[mask]
        event_counts.append(selected.shape[0])
        bj1_mass_lists.append(selected['B_J1_mass'].tolist())

    scan_df = scan_df.copy()
    scan_df["eventCount"] = event_counts
    scan_df["B_J1_mass_list"] = bj1_mass_lists

    return scan_df


def plot_counts_with_scan_steps(counts_df, scan_df, output_path):
    """Plot event counts with scan steps highlighted."""
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(counts_df['Time'], counts_df['Counts'], '.-', label='Event Count')
    for _, row in scan_df.iterrows():
        ax.axvspan(row['tStart'], row['tStop'], alpha=0.3, color='orange')

    ax.xaxis.set_major_locator(mdates.SecondLocator(bysecond=range(0, 60, 10)))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax.set_xlabel('Time')
    ax.set_ylabel('Event Count')
    ax.legend()
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.tight_layout()
    plt.title(f"Event Count with Scan Steps: {output_path.stem}")

    if SAVE_PLOTS:
        plt.savefig(output_path)
        plt.close()
    else:
        plt.show()

def process_all_scan_files():
    """Main processing loop."""
    df = load_event_data(ROOT_FILE_PATH, BRANCHES_TO_READ)
    df = apply_cuts(df, MUON_PT_CUT, MASS_CUT)

    # Build a time-binned count DataFrame for plotting
    df = df.set_index('Time')
    counts_df = df.resample(TIME_BIN).count()[["eventTime"]]
    counts_df.rename(columns={"eventTime": "Counts"}, inplace=True)
    counts_df = counts_df.reset_index()
    counts_df.to_csv(OUTPUT_CSV_DIR / "aggregated_counts.csv", index=False)

    for csv_file in SCAN_DIR.glob("scan_*_Run*.csv"):
        print(f"[INFO] Processing {csv_file.name}")
        scan_df = load_scan_file(csv_file)
        if scan_df is None:
            print(f"[ERROR] Failed to load scan file {csv_file.name}. Skipping...")
            continue

        # Filter the counts_df for plotting
        time_padding = pd.Timedelta(seconds=30)

        time_window_mask = (counts_df['Time'] >= (scan_df['tStart'].min() - time_padding)) & \
                           (counts_df['Time'] <= (scan_df['tStop'].max() + time_padding))
        filtered_counts_df = counts_df[time_window_mask]

        if filtered_counts_df['Counts'].sum() == 0:
            print(f"[WARN] No data in time range for {csv_file.name}. Skipping plot...")
            continue
        
        # Plot the counts with scan steps
        plot_counts_with_scan_steps(filtered_counts_df, scan_df, OUTPUT_PLOT_DIR / f"{csv_file.stem}_plot.png")

        # Count events per step after cut
        step_summary_df = count_events_per_step(df.reset_index(), scan_df)
        step_summary_df.to_csv(OUTPUT_CSV_DIR / f"{csv_file.stem}_summary.csv", index=False)

# ---------- MAIN ----------
if __name__ == "__main__":
    process_all_scan_files()
