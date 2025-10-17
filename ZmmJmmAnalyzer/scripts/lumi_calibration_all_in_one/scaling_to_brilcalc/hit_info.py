import os
import glob
import uproot
import numpy as np
import pandas as pd


# -----------------------------
# Config
# -----------------------------
ROOT_PATH = "/home/nimmitha/LPCfiles/run3/jmm_counting/cmssw/ZmmJmmAnalyzer/preselection"
TREE_NAME = "ntuple"
cuts = {
    "B_J1_mass": (2.7, 3.5),
    "B_Mu1_pt": (5000, None),
    "B_Mu2_pt": (5000, None)
}
branches_to_read = ["Run", "LumiBlock", "B_J1_mass", "B_Mu1_pt", "B_Mu2_pt"]
hit_branches = [
    "mu1_B1", "mu1_B2", "mu1_B3", "mu1_B4", "mu1_F",
    "mu2_B1", "mu2_B2", "mu2_B3", "mu2_B4", "mu2_F"
]
branches = branches_to_read + hit_branches

# -----------------------------
# Helper function to analyze one file
# -----------------------------
import uproot
import numpy as np
import pandas as pd

def analyze_file_streaming(filepath, tree_name="ntuple", step_size=500_000):
    """
    Process ROOT file in chunks to reduce memory usage.
    Returns hit fractions as a dictionary.
    """
    cuts = {
        "B_J1_mass": (2.7, 3.5),
        "B_Mu1_pt": (5000, None),
        "B_Mu2_pt": (5000, None),
    }
    branches_to_read = ["Run", "LumiBlock", "B_J1_mass", "B_Mu1_pt", "B_Mu2_pt"]
    hit_branches = [
        "mu1_B1", "mu1_B2", "mu1_B3", "mu1_B4", "mu1_F",
        "mu2_B1", "mu2_B2", "mu2_B3", "mu2_B4", "mu2_F"
    ]
    branches = branches_to_read + hit_branches

    # Counters
    total_muons = 0
    counts = {
        "First_L1": 0,
        "First_L2": 0,
        "First_L3": 0,
        "First_L4": 0,
        "Had_B1": 0,
        "Had_B2": 0,
        "Had_B3": 0,
        "Had_B4": 0,
    }

    # Loop over file in chunks
    for arrays in uproot.iterate(
        f"{filepath}:{tree_name}",
        branches,
        step_size=step_size,
        library="np"
    ):
        # Apply cuts
        mask = np.ones(len(arrays["Run"]), dtype=bool)
        for branch, (low, high) in cuts.items():
            if low is not None:
                mask &= arrays[branch] > low
            if high is not None:
                mask &= arrays[branch] < high

        if mask.sum() == 0:
            continue

        # Build muon arrays
        mu1 = np.stack([arrays["mu1_B1"][mask], arrays["mu1_B2"][mask],
                        arrays["mu1_B3"][mask], arrays["mu1_B4"][mask],
                        arrays["mu1_F"][mask]], axis=1)
        mu2 = np.stack([arrays["mu2_B1"][mask], arrays["mu2_B2"][mask],
                        arrays["mu2_B3"][mask], arrays["mu2_B4"][mask],
                        arrays["mu2_F"][mask]], axis=1)
        muons = np.vstack([mu1, mu2])  # shape = (n_muons, 5)

        # Only muons without FPIX hits
        # no_fpix = muons[muons[:, 4] == False]
        # temporary change to include FPIX hits as well
        no_fpix = muons  # include FPIX hits as well
        n = len(no_fpix)
        if n == 0:
            continue

        total_muons += n

        B1, B2, B3, B4 = no_fpix[:, 0], no_fpix[:, 1], no_fpix[:, 2], no_fpix[:, 3]

        # Increment counters
        counts["First_L1"] += np.sum(B1)
        counts["First_L2"] += np.sum((~B1) & B2)
        counts["First_L3"] += np.sum((~B1) & (~B2) & B3)
        counts["First_L4"] += np.sum((~B1) & (~B2) & (~B3) & B4)
        counts["Had_B1"]  += np.sum(B1)
        counts["Had_B2"]  += np.sum(B2)
        counts["Had_B3"]  += np.sum(B3)
        counts["Had_B4"]  += np.sum(B4)

    # Normalize to fractions
    if total_muons == 0:
        return {k: 0 for k in counts}
    else:
        return {k: v / total_muons for k, v in counts.items()}


# -----------------------------
# Main loop over ROOT files
# -----------------------------
summary = {}
for filepath in glob.glob(os.path.join(ROOT_PATH, "PDMLM_*.root")):
    filename = os.path.basename(filepath)
    print(f"Processing {filename} ...")
    column_name = filename.split("_")[2]
    summary[column_name] = analyze_file_streaming(filepath)

# -----------------------------
# Build summary DataFrame
# -----------------------------
summary_df = pd.DataFrame(summary)
print(summary_df.round(4))

# Optionally save
summary_df.to_csv("hit_summary.csv")
