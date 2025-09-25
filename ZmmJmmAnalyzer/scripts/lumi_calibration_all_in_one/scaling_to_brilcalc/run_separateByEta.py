import uproot
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path


def load_tree(root_file: Path, tree_name: str = "ntuple"):
    """Open ROOT file and return the tree."""
    return 


def process_chunk(arrays: dict,
                  mass_range=(2.7, 3.5), pt_cut=5000,
                  barrel_max_eta=120, forward_max_eta=240,
                  bin_width=50) -> pd.DataFrame:
    """Process a single chunk of arrays and return category counts per bin."""
    run = arrays["Run"]
    ls = arrays["LumiBlock"]
    mass = arrays["B_J1_mass"]
    mu1_pt, mu2_pt = arrays["B_Mu1_pt"], arrays["B_Mu2_pt"]
    mu1_eta, mu2_eta = np.abs(arrays["B_Mu1_eta"]), np.abs(arrays["B_Mu2_eta"])

    # Apply cuts
    mask = (
        (mass > mass_range[0]) & (mass < mass_range[1]) &
        (mu1_pt > pt_cut) & (mu2_pt > pt_cut)
    )
    run, ls, mu1_eta, mu2_eta = run[mask], ls[mask], mu1_eta[mask], mu2_eta[mask]

    if len(run) == 0:
        return pd.DataFrame()

    # Lumiblock binning
    lumiblock_bin = (ls // bin_width) * bin_width
    keys = np.char.add(np.char.add(run.astype(str), "_"), lumiblock_bin.astype(str))

    # Classify by eta regions
    mu1_region = np.where(mu1_eta < barrel_max_eta, "barrel",
                   np.where(mu1_eta < forward_max_eta, "forward", "out"))
    mu2_region = np.where(mu2_eta < barrel_max_eta, "barrel",
                   np.where(mu2_eta < forward_max_eta, "forward", "out"))

    categories = np.full(len(mu1_region), "other", dtype=object)
    categories[(mu1_region == "barrel") & (mu2_region == "barrel")] = "barrel-barrel"
    categories[(mu1_region == "forward") & (mu2_region == "forward")] = "forward-forward"
    categories[((mu1_region == "barrel") & (mu2_region == "forward")) |
               ((mu1_region == "forward") & (mu2_region == "barrel"))] = "mixed"

    # Build dataframe
    chunk_df = pd.DataFrame({"key": keys, "category": categories})
    chunk_df = chunk_df[chunk_df["category"].isin(
        ["barrel-barrel", "forward-forward", "mixed"]
    )]
    if chunk_df.empty:
        return pd.DataFrame()

    counts = chunk_df.value_counts(["key", "category"]).unstack(fill_value=0)
    counts["total"] = counts.sum(axis=1)
    return counts


def classify_candidates(root_file: Path,
                        mass_range=(2.7, 3.5), pt_cut=5000,
                        barrel_max_eta=120, forward_max_eta=240,
                        chunk_size=1_000_000, bin_width=50,
                        max_chunks=None) -> pd.DataFrame:
    """
    Classify J/psi candidates by muon eta regions.
    Returns aggregated DataFrame per (run, lumiblock_bin).
    """
    tree = uproot.open(root_file)['ntuple']
    branches = ["Run", "LumiBlock", "B_J1_mass",
                "B_Mu1_pt", "B_Mu2_pt", "B_Mu1_eta", "B_Mu2_eta"]

    results = pd.DataFrame()
    nentries = tree.num_entries

    with tqdm(total=nentries, desc="Processing", unit="entries") as pbar:
        for i, arrays in enumerate(tree.iterate(branches, step_size=chunk_size, library="np")):
            counts = process_chunk(
                arrays, mass_range, pt_cut, barrel_max_eta, forward_max_eta, bin_width
            )
            if not counts.empty:
                results = results.add(counts, fill_value=0) if not results.empty else counts

            pbar.update(len(arrays["Run"]))

            if max_chunks and (i + 1) >= max_chunks:  # useful for testing
                break

    # Finalize results
    if results.empty:
        return pd.DataFrame()

    results.index.name = "bin_key"
    for col in ["barrel-barrel", "forward-forward", "mixed"]:
        if col not in results.columns:
            results[col] = 0

    for col in ["barrel-barrel", "forward-forward", "mixed"]:
        results[col + "_frac"] = results[col] / results["total"]

    return results.reset_index()


def main():
    base_path = Path(
        "/home/nimmitha/LPCfiles/run3/jmm_counting/cmssw/"
        "ZmmJmmAnalyzer/preselection/parkingDoubleMuonLowMass/2024/v6/"
    )
    file_name = "PDMLM_mm_2024G1_v6.root"
    file_path = base_path / file_name

    df = classify_candidates(
        file_path,
        mass_range=(2.7, 3.5),
        pt_cut=5000,
        barrel_max_eta=120,   # |eta| < 1.2
        forward_max_eta=240,  # |eta| < 2.4
        max_chunks=10       # set small int for testing
    )

    print(df.head())
    out_csv = f"csvs/fractionsPerEta/jpsi_regions_{file_name.replace('.root', '')}.csv"
    df.to_csv(out_csv, index=False)
    print(f"Results saved to {out_csv}")


if __name__ == "__main__":
    main()
