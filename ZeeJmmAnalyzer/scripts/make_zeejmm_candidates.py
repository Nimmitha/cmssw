#!/usr/bin/env python3
"""
Build flat ZeeJmm candidate trees from miniAODeemm preselection ntuples using PyROOT.

This is intentionally simple and reviewable:
  * configure input/output paths at the top;
  * loop over events and candidate index i = 0..nB-1;
  * apply targeted ZeeJmm cuts;
  * write one scalar row per selected candidate.

Expected input tree: ntuple
Expected branch convention: miniAODeemm analyzer output.

Nominal recommended use:
  python3 make_zeejmm_candidates.py

Before large production, run on a small file and check printed cutflow/yields.
"""

from __future__ import annotations

from array import array
from pathlib import Path
from typing import Dict, List

import ROOT

ROOT.gROOT.SetBatch(True)

# =============================================================================
# User configuration
# =============================================================================

TREE_NAME = "ntuple"
OUTDIR = "selection"

MAKE_SIGNAL = True
MAKE_BACKGROUND = True
MAKE_FINAL_BLINDED = True
MAKE_FINAL_UNBLINDED = True  # keep False until ready to inspect/unblind

SIGNAL_PRESELECTION = "preselection/zeejmm_mc_2018_v2.root"
DATA_PRESELECTION = "preselection/TTree_13TeV_eemm_UL_Run2.root"

OUTPUTS = {
    "signal": "signal_candidates.root",
    "background": "background_candidates.root",
    "final_blinded": "final_blinded_candidates.root",
    "final_unblinded": "final_unblinded_candidates.root",
}

# Analysis and blinded Higgs windows
ANALYSIS_LOW = 112.0
ANALYSIS_HIGH = 142.0
MASK_LOW = 120.0
MASK_HIGH = 130.0

# Core physics cuts from ZeeJmm studies / AN-style preselection
REQUIRE_ELE_TRIGGER_MC = True
REQUIRE_TRIGGER_MATCH = False  # set True only if you explicitly want offline electron-HLT matching
REQUIRE_SOFT_MUONS = True
REQUIRE_ELECTRON_ID = True
ELECTRON_ID = "Loose"  # options: "Loose", "WP90", "WP80"

MUON_PT_MIN = 3.0
MUON_ABS_ETA_MAX = 2.4
ELE1_PT_MIN = 32.0
ELE2_PT_MIN = 5.0
ELE_ABS_ETA_MAX = 2.5

PAIR_VTXPROB_MIN = 0.01
FOURL_VTXPROB_MIN = 0.01
FOURL_PT_MIN = 5.0

SIGNAL_JPSI_MASS = (3.0, 3.2)
DATA_JPSI_MASS = (3.0, 3.2)
Z_MASS = (80.0, 100.0)
SIGNAL_FOURL_MASS = (112.0, 142.0)

# Candidate variants to write.
#   onecand: keep one best candidate per input TTree entry
#   allcand: keep all selected candidates
# Best-candidate rule: highest fourL_vtxProb only.
MAKE_ONECAND = True
MAKE_ALLCAND = True

CANDIDATE_VARIANTS = []
if MAKE_ONECAND:
    CANDIDATE_VARIANTS.append(("onecand", True))
if MAKE_ALLCAND:
    CANDIDATE_VARIANTS.append(("allcand", False))

# =============================================================================
# Output content
# =============================================================================

FLOAT_BRANCHES = [
    "fourL_mass", "fourL_pt", "fourL_eta", "fourL_phi", "fourL_rapidity", "fourL_vtxProb",
    "Z_mass", "Z_pt", "Z_eta", "Z_phi", "Z_rapidity", "Z_vtxProb", "Z_dR_ee",
    "Jpsi_mass", "Jpsi_pt", "Jpsi_eta", "Jpsi_phi", "Jpsi_rapidity", "Jpsi_vtxProb", "Jpsi_dR_mumu",
    "Z_Jpsi_dR", "Z_Jpsi_dPhi", "Z_Jpsi_dEta", "Z_Jpsi_dY", "pt_balance",
    "cosTheta_Z_ePlus", "cosTheta_Jpsi_muPlus", "phi_decayPlane_Z_Jpsi",
    "e1_pt", "e1_eta", "e1_phi", "e1_dxy", "e1_dz", "e1_mvaRaw", "e1_triggerDR",
    "e2_pt", "e2_eta", "e2_phi", "e2_dxy", "e2_dz", "e2_mvaRaw", "e2_triggerDR",
    "mu1_pt", "mu1_eta", "mu1_phi", "mu1_pfRelIso03", "mu1_dxy", "mu1_dz", "mu1_dB3D",
    "mu2_pt", "mu2_eta", "mu2_phi", "mu2_pfRelIso03", "mu2_dxy", "mu2_dz", "mu2_dB3D",
    "Jpsi_trackIso03", "Jpsi_relIso03", "Z_trackIso03", "Z_relIso03",
]

INT_BRANCHES = [
    "Run", "LumiBlock", "Event", "nPV", "label",
    "passEleTrigger", "passEleTriggerMatch",
    "e1_charge", "e2_charge", "e1_passLooseID", "e2_passLooseID", "e1_passWP90", "e2_passWP90", "e1_passWP80", "e2_passWP80",
    "e1_triggerMatched", "e2_triggerMatched",
    "mu1_charge", "mu2_charge", "mu1_soft", "mu2_soft", "mu1_loose", "mu2_loose", "mu1_tight", "mu2_tight",
    "nExtraLooseElectrons", "nExtraLooseMuons",
]

# =============================================================================
# Helpers
# =============================================================================


def in_window(x: float, lo: float, hi: float) -> bool:
    return lo < x < hi


def outside_mask_window(mass: float) -> bool:
    return not in_window(mass, MASK_LOW, MASK_HIGH)


def analysis_window(mass: float) -> bool:
    return in_window(mass, ANALYSIS_LOW, ANALYSIS_HIGH)


def analysis_sideband(mass: float) -> bool:
    return analysis_window(mass) and outside_mask_window(mass)


def get_vec_value(tree: ROOT.TTree, branch: str, idx: int):
    return getattr(tree, branch).at(idx)


def require_branch(tree: ROOT.TTree, branch: str) -> None:
    if not tree.GetBranch(branch):
        raise RuntimeError(f"Missing required branch: {branch}")


def electron_id_pass(tree: ROOT.TTree, idx: int) -> bool:
    if ELECTRON_ID == "Loose":
        return bool(get_vec_value(tree, "e1_passLooseID", idx)) and bool(get_vec_value(tree, "e2_passLooseID", idx))
    if ELECTRON_ID == "WP90":
        return bool(get_vec_value(tree, "e1_passWP90", idx)) or bool(get_vec_value(tree, "e2_passWP90", idx))
    if ELECTRON_ID == "WP80":
        return bool(get_vec_value(tree, "e1_passWP80", idx)) or bool(get_vec_value(tree, "e2_passWP80", idx))
    raise ValueError(f"Unknown ELECTRON_ID = {ELECTRON_ID}")


def common_selection(tree: ROOT.TTree, idx: int, sample: str) -> bool:
    """Cuts common to signal, background, and final data candidates."""
    require_ele_trigger = sample != "signal" or REQUIRE_ELE_TRIGGER_MC
    if require_ele_trigger and not bool(get_vec_value(tree, "passEleTrigger", idx)):
        return False

    if REQUIRE_TRIGGER_MATCH and not bool(get_vec_value(tree, "passEleTriggerMatch", idx)):
        return False

    if REQUIRE_SOFT_MUONS:
        if not bool(get_vec_value(tree, "mu1_soft", idx)):
            return False
        if not bool(get_vec_value(tree, "mu2_soft", idx)):
            return False

    if REQUIRE_ELECTRON_ID and not electron_id_pass(tree, idx):
        return False

    if get_vec_value(tree, "mu1_pt", idx) <= MUON_PT_MIN:
        return False
    if get_vec_value(tree, "mu2_pt", idx) <= MUON_PT_MIN:
        return False
    if abs(get_vec_value(tree, "mu1_eta", idx)) >= MUON_ABS_ETA_MAX:
        return False
    if abs(get_vec_value(tree, "mu2_eta", idx)) >= MUON_ABS_ETA_MAX:
        return False

    # e1/e2 are pT-ordered by the analyzer. With the EGamma-corrected analyzer,
    # these pT values are based on ecalTrkEnergyPostCorr.
    if get_vec_value(tree, "e1_pt", idx) <= ELE1_PT_MIN:
        return False
    if get_vec_value(tree, "e2_pt", idx) <= ELE2_PT_MIN:
        return False
    if abs(get_vec_value(tree, "e1_eta", idx)) >= ELE_ABS_ETA_MAX:
        return False
    if abs(get_vec_value(tree, "e2_eta", idx)) >= ELE_ABS_ETA_MAX:
        return False

    if get_vec_value(tree, "Z_vtxProb", idx) <= PAIR_VTXPROB_MIN:
        return False
    if get_vec_value(tree, "Jpsi_vtxProb", idx) <= PAIR_VTXPROB_MIN:
        return False
    if get_vec_value(tree, "fourL_vtxProb", idx) <= FOURL_VTXPROB_MIN:
        return False
    if get_vec_value(tree, "fourL_pt", idx) <= FOURL_PT_MIN:
        return False

    if not in_window(get_vec_value(tree, "Z_mass", idx), *Z_MASS):
        return False

    return True


def sample_selection(tree: ROOT.TTree, idx: int, sample: str) -> bool:
    if not common_selection(tree, idx, sample):
        return False

    jmass = get_vec_value(tree, "Jpsi_mass", idx)
    four_mass = get_vec_value(tree, "fourL_mass", idx)

    if sample == "signal":
        return (
            in_window(jmass, *SIGNAL_JPSI_MASS)
            and in_window(four_mass, *SIGNAL_FOURL_MASS)
        )

    if sample == "background":
        return (
            in_window(jmass, *DATA_JPSI_MASS)
            and outside_mask_window(four_mass)
        )

    if sample == "final_blinded":
        return (
            in_window(jmass, *DATA_JPSI_MASS)
            and analysis_sideband(four_mass)
        )

    if sample == "final_unblinded":
        return (
            in_window(jmass, *DATA_JPSI_MASS)
            and analysis_window(four_mass)
        )

    raise ValueError(f"Unknown sample: {sample}")


def candidate_rank(tree: ROOT.TTree, idx: int) -> tuple:
    """Lower tuple is better for optional per-entry deduplication.

    Per your current choice, the only ranking variable is fourL_vtxProb.
    We use the negative value because Python's min(...) returns the smallest tuple.
    Do not include closeness to 125 GeV here, because that would sculpt the
    four-lepton mass distribution.
    """
    four_vtx = get_vec_value(tree, "fourL_vtxProb", idx)
    return (-four_vtx,)


def make_output_tree() -> tuple[ROOT.TFile, ROOT.TTree, Dict[str, array]]:
    arrays: Dict[str, array] = {}
    fout = None  # created by caller after path is known; kept for typing clarity
    tout = ROOT.TTree("ntuple", "selected ZeeJmm candidates")

    for name in FLOAT_BRANCHES:
        arrays[name] = array("f", [0.0])
        tout.Branch(name, arrays[name], f"{name}/F")

    for name in INT_BRANCHES:
        # Event can exceed signed 32-bit; store event identifiers as unsigned long long for safety.
        if name in {"Run", "LumiBlock", "Event"}:
            arrays[name] = array("L", [0])
            tout.Branch(name, arrays[name], f"{name}/l")
        else:
            arrays[name] = array("i", [0])
            tout.Branch(name, arrays[name], f"{name}/I")

    return fout, tout, arrays


def fill_output(tree: ROOT.TTree, idx: int, out_arrays: Dict[str, array], label: int) -> None:
    out_arrays["label"][0] = label

    for name in FLOAT_BRANCHES:
        out_arrays[name][0] = float(get_vec_value(tree, name, idx))

    for name in INT_BRANCHES:
        if name == "label":
            continue
        out_arrays[name][0] = int(get_vec_value(tree, name, idx))


def selected_indices_for_entry(tree: ROOT.TTree, sample: str, deduplicate: bool) -> List[int]:
    """Return selected candidate indices for one input TTree entry.

    If deduplicate=True, keep only the candidate with the highest fourL_vtxProb
    within this entry. This avoids using Run/LumiBlock/Event as a global key, which is
    important for your private MC where those numbers can repeat after merging.
    """
    n_cands = int(tree.nB)
    selected = [i for i in range(n_cands) if sample_selection(tree, i, sample)]

    if not deduplicate or len(selected) <= 1:
        return selected

    best = min(selected, key=lambda i: candidate_rank(tree, i))
    return [best]


def process_sample(sample: str, input_path: str, output_path: str, label: int, deduplicate: bool) -> None:
    print(f"\nProcessing {sample}")
    print(f"  input : {input_path}")
    print(f"  output: {output_path}")

    fin = ROOT.TFile.Open(input_path)
    if not fin or fin.IsZombie():
        raise RuntimeError(f"Could not open input file: {input_path}")

    tree = fin.Get(TREE_NAME)
    if not tree:
        raise RuntimeError(f"Could not find tree '{TREE_NAME}' in {input_path}")

    required = ["nB"] + FLOAT_BRANCHES + [b for b in INT_BRANCHES if b != "label"]
    for branch in sorted(set(required)):
        require_branch(tree, branch)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fout = ROOT.TFile.Open(output_path, "RECREATE")
    _, tout, out_arrays = make_output_tree()

    n_events = int(tree.GetEntries())
    n_candidates = 0
    n_selected_entries_before_dedup = 0
    n_selected_candidates_before_dedup = 0
    n_saved_entries = 0
    n_saved_candidates = 0

    for iev in range(n_events):
        tree.GetEntry(iev)
        n_candidates += int(tree.nB)

        selected_before = selected_indices_for_entry(tree, sample, deduplicate=False)
        if selected_before:
            n_selected_entries_before_dedup += 1
            n_selected_candidates_before_dedup += len(selected_before)

        indices_to_save = selected_indices_for_entry(
            tree, sample, deduplicate=deduplicate
        )
        if indices_to_save:
            n_saved_entries += 1

        for idx in indices_to_save:
            fill_output(tree, idx, out_arrays, label)
            tout.Fill()
            n_saved_candidates += 1

        if iev > 0 and iev % 100000 == 0:
            print(
                f"    entries {iev}/{n_events}, "
                f"selected candidates before dedup {n_selected_candidates_before_dedup}, "
                f"saved {n_saved_candidates}",
                flush=True,
            )

    fout.cd()
    tout.Write()
    output_entries = int(tout.GetEntries())
    fout.Close()
    fin.Close()

    print(f"  input TTree entries              : {n_events}")
    print(f"  input candidates                 : {n_candidates}")
    print(f"  selected entries before dedup    : {n_selected_entries_before_dedup}")
    print(f"  selected candidates before dedup : {n_selected_candidates_before_dedup}")
    print(f"  saved entries                    : {n_saved_entries}")
    print(f"  saved candidates                 : {n_saved_candidates}")
    print(f"  output tree entries              : {output_entries}")
    print(f"  saved file                       : {output_path}")


def main() -> None:
    if not CANDIDATE_VARIANTS:
        raise RuntimeError("CANDIDATE_VARIANTS is empty.")

    if not (MAKE_SIGNAL or MAKE_BACKGROUND or MAKE_FINAL_BLINDED or MAKE_FINAL_UNBLINDED):
        raise RuntimeError("No outputs requested. Enable at least one MAKE_* switch.")

    print("ZeeJmm targeted selection")
    print(f"  Electron ID: {ELECTRON_ID}")
    print("  Ele trigger required for data: True")
    print(f"  Ele trigger required for MC: {REQUIRE_ELE_TRIGGER_MC}")
    print(f"  Trigger match required: {REQUIRE_TRIGGER_MATCH}")
    print("  Best-candidate rule: highest fourL_vtxProb")

    base_outdir = Path(OUTDIR)

    for variant_name, deduplicate in CANDIDATE_VARIANTS:
        outdir = base_outdir / variant_name
        print("\n" + "=" * 80)
        print(f"Candidate variant: {variant_name}")
        print(f"  output directory: {outdir}")
        print(f"  one candidate per input entry: {deduplicate}")
        print("=" * 80)

        if MAKE_SIGNAL:
            process_sample("signal", SIGNAL_PRESELECTION, str(outdir / OUTPUTS["signal"]), label=1, deduplicate=deduplicate)

        # Read data once per output. This is slightly less efficient than the Zmm script,
        # but much easier to review and debug.
        if MAKE_BACKGROUND:
            process_sample("background", DATA_PRESELECTION, str(outdir / OUTPUTS["background"]), label=0, deduplicate=deduplicate)

        if MAKE_FINAL_BLINDED:
            process_sample("final_blinded", DATA_PRESELECTION, str(outdir / OUTPUTS["final_blinded"]), label=-1, deduplicate=deduplicate)

        if MAKE_FINAL_UNBLINDED:
            process_sample("final_unblinded", DATA_PRESELECTION, str(outdir / OUTPUTS["final_unblinded"]), label=-1, deduplicate=deduplicate)


if __name__ == "__main__":
    main()
