#!/usr/bin/env python3
"""
Build flat Z(mm)J/psi(mm) candidate trees from the cleaned miniAODmmmm preselection ntuples.

Design choices:
  * One script with hardcoded paths/switches at the top.
  * Reads two preselection files: merged Run-2 data and signal MC.
  * Produces signal once, and processes the data file once to make background/final outputs.
  * Does not write string branches, for compatibility with older uproot writers.
  * Uses the old-compatible pairing choice by default.
  * Applies the combined old preselection+selection cuts in Task 2, while Task 1 remains broad.
  * Can keep all candidates or reduce to one candidate per event.

Required packages:
  pip install uproot awkward numpy

Expected input tree:
  ntuple

Expected preselection branch convention:
  fourMu_*, pair12_*, pair34_*, pair23_*, pair14_*, muP1_*, muM1_*, muP2_*, muM2_*

Pair groups from preselection:
  Group A: pair12 = muP1 + muM1, pair34 = muP2 + muM2
  Group B: pair14 = muP1 + muM2, pair23 = muM1 + muP2

The selected J/psi/Z assignment may be pair12/pair34, pair34/pair12,
pair14/pair23, or pair23/pair14, depending on masses and vertex quality.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from time import time
from typing import Dict, Iterable, List, Optional, Tuple

import awkward as ak
import numpy as np
import uproot

# =============================================================================
# User configuration
# =============================================================================

TREE_NAME = "ntuple"
OUTDIR = "selection"

MAKE_SIGNAL = True
MAKE_BACKGROUND = True
MAKE_FINAL_BLINDED = True
MAKE_FINAL_UNBLINDED = True  # keep False until you are ready to inspect/unblind

SIGNAL_PRESELECTION = "preselection/zmmjmm_mc_2018_v2.root"
DATA_PRESELECTION = "prep/skimmed_TTree_13TeV_mmmm_UL_Run2.root"

OUTPUTS = {
    "signal": "signal_candidates.root",
    "background": "background_candidates.root",
    "final_blinded": "final_blinded_candidates.root",
    "final_unblinded": "final_unblinded_candidates.root",
}

# Candidate variants to write.
#   onecand: reduce to one candidate per (Run, LumiBlock, Event) after final cuts
#   allcand: keep all selected candidates
MAKE_ONECAND = True
MAKE_ALLCAND = True

# Analysis and blinded Higgs windows
ANALYSIS_LOW = 112.0
ANALYSIS_HIGH = 162.0
MASK_LOW = 120.0
MASK_HIGH = 130.0

# Optional broad four-muon range for background training files.
# Currently the default background sample uses analysis_sideband, so this is only
# used if you switch a sample fourmu_region to background_masked later.
BACKGROUND_FOURMU_RANGE = None  # example: (80.0, 200.0)

# Pairing score constants. These are practical relative weights for old-compatible
# J/psi/Z assignment, not detector-resolution measurements.
M_JPSI = 3.0969
M_Z = 91.1876
SIGMA_JPSI_SCORE = 0.10
SIGMA_Z_SCORE = 10.0
MASS_SCORE_CLOSE = 0.5

# Broad topology cuts. These mirror the miniAOD preselector.
BROAD_LOW_MASS = (2.0, 4.0)
BROAD_Z_MASS = (60.0, 120.0)
BROAD_PAIR_VTXPROB_MIN = 0.005
BROAD_FOURMU_VTXPROB_MIN = 0.005

# Final physics-selection cuts.
SIGNAL_JPSI_MASS = (3.0, 3.2)
DATA_JPSI_MASS = (2.8, 3.4)
SIGNAL_Z_MASS = (80.0, 100.0)
DATA_Z_MASS = (70.0, 110.0)

MUON_PT_MIN = 3.0
MUON_ABS_ETA_MAX = 2.4
REQUIRE_SOFT_MUONS = True

PAIR_VTXPROB_MIN = 0.01
FOURMU_VTXPROB_MIN = 0.01
DIMUON_PT_MIN = 5.0
FOURMU_PT_MIN = 5.0

BACKGROUND_FOURMU_REGION = "analysis_sideband"
FINAL_BLINDED_FOURMU_REGION = "analysis_sideband"
FINAL_UNBLINDED_FOURMU_REGION = "analysis_window"
SIGNAL_FOURMU_REGION = "analysis_window"

# Chunk size for uproot.iterate. Increase if memory is fine.
STEP_SIZE = "100 MB"

# Progress printing for large files.
PRINT_PROGRESS = True
PROGRESS_EVERY_CHUNKS = 1

# Derived list of output variants. Usually you should only edit MAKE_ONECAND and MAKE_ALLCAND above.
CANDIDATE_VARIANTS = []
if MAKE_ONECAND:
    CANDIDATE_VARIANTS.append(("onecand", True))
if MAKE_ALLCAND:
    CANDIDATE_VARIANTS.append(("allcand", False))

# =============================================================================
# Internal helpers
# =============================================================================


@dataclass(frozen=True)
class SampleConfig:
    name: str
    input_path: str
    output_path: str
    jpsi_mass_window: Tuple[float, float]
    z_mass_window: Tuple[float, float]
    fourmu_region: str
    label: int


def delta_phi(phi1: np.ndarray, phi2: np.ndarray) -> np.ndarray:
    """Return signed delta phi in [-pi, pi)."""
    return (phi1 - phi2 + np.pi) % (2.0 * np.pi) - np.pi


def in_window(x: np.ndarray, window: Tuple[float, float]) -> np.ndarray:
    lo, hi = window
    return (x > lo) & (x < hi)


def mass_score(jpsi_mass: np.ndarray, z_mass: np.ndarray) -> np.ndarray:
    return np.abs(jpsi_mass - M_JPSI) / SIGMA_JPSI_SCORE + np.abs(z_mass - M_Z) / SIGMA_Z_SCORE


def safe_vertex_product(jpsi_vtx: np.ndarray, z_vtx: np.ndarray) -> np.ndarray:
    j = np.where(np.isfinite(jpsi_vtx) & (jpsi_vtx > 0.0), jpsi_vtx, 0.0)
    z = np.where(np.isfinite(z_vtx) & (z_vtx > 0.0), z_vtx, 0.0)
    return j * z


def fourmu_region_mask(fourmu_mass: np.ndarray, region: str) -> np.ndarray:
    """Four-muon mass-region masks.

    Regions:
      all               : no four-muon mass cut
      analysis_window   : 112--162 GeV
      mask_window       : 120--130 GeV
      analysis_sideband : 112--162 GeV excluding 120--130 GeV
      outside_mask           : all masses excluding 120--130 GeV
      outside_analysis       : outside 112--162 GeV
      background_masked      : optional background range, excluding 120--130
    """
    analysis = (fourmu_mass > ANALYSIS_LOW) & (fourmu_mass < ANALYSIS_HIGH)
    mask = (fourmu_mass > MASK_LOW) & (fourmu_mass < MASK_HIGH)

    if BACKGROUND_FOURMU_RANGE is None:
        background_range = np.ones_like(fourmu_mass, dtype=bool)
    else:
        blo, bhi = BACKGROUND_FOURMU_RANGE
        background_range = (fourmu_mass > blo) & (fourmu_mass < bhi)

    if region == "all":
        return np.ones_like(fourmu_mass, dtype=bool)
    if region == "analysis_window":
        return analysis
    if region == "mask_window":
        return mask
    if region == "analysis_sideband":
        return analysis & (~mask)
    if region == "outside_mask":
        return ~mask
    if region == "outside_analysis":
        return ~analysis
    if region == "background_masked":
        return background_range & (~mask)
    raise ValueError(f"Unknown fourmu_region: {region}")


def has_branch(arrays: ak.Array, name: str) -> bool:
    """Return True if an uproot/awkward record array has this branch/field."""
    return name in arrays.fields


def flatten_branch(arrays: ak.Array, name: str, dtype=None) -> np.ndarray:
    if not has_branch(arrays, name):
        available = ", ".join(arrays.fields[:20])
        raise KeyError(f"Missing required input branch: {name}. First available fields: {available}")
    values = ak.to_numpy(ak.flatten(arrays[name], axis=None))
    if dtype is not None:
        values = values.astype(dtype, copy=False)
    return values


def event_level_flatten(arrays: ak.Array, branch: str, dtype=None) -> np.ndarray:
    """Repeat event-level vector values to match candidate-level flattened arrays.

    In the preselection, Run/Lumi/Event/nPV are stored as vectors filled once per candidate,
    so a normal flatten is enough. This helper exists in case that changes later.
    """
    return flatten_branch(arrays, branch, dtype=dtype)


def required_branches() -> List[str]:
    branches = [
        "Run", "LumiBlock", "Event", "nPV", "passMuonTrigger",
        "fourMu_mass", "fourMu_pt", "fourMu_eta", "fourMu_phi", "fourMu_rapidity", "fourMu_vtxProb",
        "pair12_mass", "pair12_pt", "pair12_eta", "pair12_phi", "pair12_rapidity", "pair12_vtxProb",
        "pair34_mass", "pair34_pt", "pair34_eta", "pair34_phi", "pair34_rapidity", "pair34_vtxProb",
        "pair23_mass", "pair23_pt", "pair23_eta", "pair23_phi", "pair23_rapidity", "pair23_vtxProb",
        "pair14_mass", "pair14_pt", "pair14_eta", "pair14_phi", "pair14_rapidity", "pair14_vtxProb",
        "pair12_dR", "pair34_dR", "pair23_dR", "pair14_dR",
        "pair12_trackIso03", "pair12_trackIso04", "pair12_relIso03", "pair12_relIso04",
        "pair34_trackIso03", "pair34_trackIso04", "pair34_relIso03", "pair34_relIso04",
        "pair23_trackIso03", "pair23_trackIso04", "pair23_relIso03", "pair23_relIso04",
        "pair14_trackIso03", "pair14_trackIso04", "pair14_relIso03", "pair14_relIso04",
        "pair12_cosThetaMu", "pair34_cosThetaMu", "pair23_cosThetaMu", "pair14_cosThetaMu",
        "pair12_34_phiDecayPlane", "pair23_14_phiDecayPlane",
    ]

    for mu in ("muP1", "muM1", "muP2", "muM2"):
        branches += [
            f"{mu}_pt", f"{mu}_eta", f"{mu}_phi", f"{mu}_charge",
            f"{mu}_soft", f"{mu}_tight", f"{mu}_loose",
            f"{mu}_pfRelIso03", f"{mu}_pfRelIso04", f"{mu}_trackAbsIso03", f"{mu}_trackRelIso03", f"{mu}_dB3D",
            f"{mu}_dxy", f"{mu}_dz", f"{mu}_normChi2", f"{mu}_nValidHits", f"{mu}_nValidPixelHits",
        ]
    return branches


def pairing_options(flat: Dict[str, np.ndarray]) -> List[Dict[str, object]]:
    """Return the four possible J/psi/Z assignments.

    Codes:
      0: Jpsi=pair12, Z=pair34
      1: Jpsi=pair34, Z=pair12
      2: Jpsi=pair14, Z=pair23
      3: Jpsi=pair23, Z=pair14
    """
    return [
        {"code": 0, "jpsi": "pair12", "z": "pair34", "plane": "pair12_34_phiDecayPlane",
         "jpsiP": "muP1", "jpsiM": "muM1", "zP": "muP2", "zM": "muM2"},
        {"code": 1, "jpsi": "pair34", "z": "pair12", "plane": "pair12_34_phiDecayPlane",
         "jpsiP": "muP2", "jpsiM": "muM2", "zP": "muP1", "zM": "muM1"},
        {"code": 2, "jpsi": "pair14", "z": "pair23", "plane": "pair23_14_phiDecayPlane",
         "jpsiP": "muP1", "jpsiM": "muM2", "zP": "muP2", "zM": "muM1"},
        {"code": 3, "jpsi": "pair23", "z": "pair14", "plane": "pair23_14_phiDecayPlane",
         "jpsiP": "muP2", "jpsiM": "muM1", "zP": "muP1", "zM": "muM2"},
    ]


def broad_topology_mask(flat: Dict[str, np.ndarray]) -> np.ndarray:
    """Preselection topology: one low-mass dimuon and one Z-like dimuon.

    This is applied on either disjoint pair group, with either pair allowed to be
    the low-mass object.  It is deliberately broader than the final J/psi mass
    selection.
    """
    p12_low = in_window(flat["pair12_mass"], BROAD_LOW_MASS)
    p34_low = in_window(flat["pair34_mass"], BROAD_LOW_MASS)
    p14_low = in_window(flat["pair14_mass"], BROAD_LOW_MASS)
    p23_low = in_window(flat["pair23_mass"], BROAD_LOW_MASS)

    p12_z = in_window(flat["pair12_mass"], BROAD_Z_MASS)
    p34_z = in_window(flat["pair34_mass"], BROAD_Z_MASS)
    p14_z = in_window(flat["pair14_mass"], BROAD_Z_MASS)
    p23_z = in_window(flat["pair23_mass"], BROAD_Z_MASS)

    group_a_mass = (p12_low & p34_z) | (p34_low & p12_z)
    group_b_mass = (p14_low & p23_z) | (p23_low & p14_z)

    group_a_vtx = (flat["pair12_vtxProb"] > BROAD_PAIR_VTXPROB_MIN) & (flat["pair34_vtxProb"] > BROAD_PAIR_VTXPROB_MIN)
    group_b_vtx = (flat["pair14_vtxProb"] > BROAD_PAIR_VTXPROB_MIN) & (flat["pair23_vtxProb"] > BROAD_PAIR_VTXPROB_MIN)
    four_vtx = flat["fourMu_vtxProb"] > BROAD_FOURMU_VTXPROB_MIN

    return four_vtx & ((group_a_mass & group_a_vtx) | (group_b_mass & group_b_vtx))


def compute_pairing_arrays(flat: Dict[str, np.ndarray], cfg: SampleConfig, opts: List[Dict[str, object]]):
    scores = []
    valid = []
    vtx_products = []
    for opt in opts:
        jp = str(opt["jpsi"])
        zz = str(opt["z"])
        scores.append(mass_score(flat[f"{jp}_mass"], flat[f"{zz}_mass"]))
        valid.append(in_window(flat[f"{jp}_mass"], cfg.jpsi_mass_window) & in_window(flat[f"{zz}_mass"], cfg.z_mass_window))
        vtx_products.append(safe_vertex_product(flat[f"{jp}_vtxProb"], flat[f"{zz}_vtxProb"]))

    return np.vstack(scores), np.vstack(valid), np.vstack(vtx_products)


def choose_best_pairing_indices(score_arr: np.ndarray, valid_arr: np.ndarray, vtx_arr: np.ndarray) -> np.ndarray:
    n_options, n_candidates = score_arr.shape

    # Score-mode reference, kept here for later studies but not active:
    #   penalized = score_arr + np.where(valid_arr, 0.0, 1e6)
    #   no_valid = ~np.any(valid_arr, axis=0)
    #   penalized[:, no_valid] = score_arr[:, no_valid]
    #   return np.argmin(penalized, axis=0)

    # Old-compatible default: prefer assignments passing the final mass windows.
    # Among valid options, lower mass score wins unless close, then higher
    # vertex-product wins.
    best_idx = np.zeros(n_candidates, dtype=np.int32)
    for i in range(n_candidates):
        candidates = np.nonzero(valid_arr[:, i])[0]
        if len(candidates) == 0:
            candidates = np.arange(n_options)

        best = int(candidates[0])
        for cand in candidates[1:]:
            cand = int(cand)
            if score_arr[cand, i] < score_arr[best, i] - MASS_SCORE_CLOSE:
                best = cand
            elif abs(score_arr[cand, i] - score_arr[best, i]) <= MASS_SCORE_CLOSE:
                if vtx_arr[cand, i] > vtx_arr[best, i]:
                    best = cand
        best_idx[i] = best
    return best_idx


def choose_pair_branch(flat: Dict[str, np.ndarray], opts: List[Dict[str, object]], best_idx: np.ndarray, suffix: str, kind: str) -> np.ndarray:
    vals = np.empty(len(best_idx), dtype=np.float32)
    for idx, opt in enumerate(opts):
        pair = str(opt[kind])
        mask = best_idx == idx
        vals[mask] = flat[f"{pair}_{suffix}"][mask]
    return vals


def choose_muon_branch(
    flat: Dict[str, np.ndarray],
    opts: List[Dict[str, object]],
    best_idx: np.ndarray,
    var: str,
    role: str,
    dtype,
) -> np.ndarray:
    vals = np.empty(len(best_idx), dtype=dtype)
    for idx, opt in enumerate(opts):
        mu = str(opt[role])
        mask = best_idx == idx
        vals[mask] = flat[f"{mu}_{var}"][mask]
    return vals


def build_selected_pair_branches(flat: Dict[str, np.ndarray], opts: List[Dict[str, object]], best_idx: np.ndarray) -> Dict[str, np.ndarray]:
    selected: Dict[str, np.ndarray] = {}

    for suffix in ("mass", "pt", "eta", "phi", "rapidity", "vtxProb"):
        selected[f"jpsi_{suffix}"] = choose_pair_branch(flat, opts, best_idx, suffix, "jpsi")
        selected[f"z_{suffix}"] = choose_pair_branch(flat, opts, best_idx, suffix, "z")

    selected["jpsi_dR_mumu"] = choose_pair_branch(flat, opts, best_idx, "dR", "jpsi")
    selected["z_dR_mumu"] = choose_pair_branch(flat, opts, best_idx, "dR", "z")
    for suffix in ("trackIso03", "trackIso04", "relIso03", "relIso04", "cosThetaMu"):
        selected[f"jpsi_{suffix}"] = choose_pair_branch(flat, opts, best_idx, suffix, "jpsi")
        selected[f"z_{suffix}"] = choose_pair_branch(flat, opts, best_idx, suffix, "z")

    selected["cosTheta_jpsiMuP"] = selected.pop("jpsi_cosThetaMu")
    selected["cosTheta_zMuP"] = selected.pop("z_cosThetaMu")

    phi_plane = np.empty(len(best_idx), dtype=np.float32)
    for idx, opt in enumerate(opts):
        mask = best_idx == idx
        phi_plane[mask] = flat[str(opt["plane"])][mask]
    selected["phi_decayPlane_z_jpsi"] = phi_plane

    selected["z_jpsi_dR"] = np.sqrt((selected["z_eta"] - selected["jpsi_eta"]) ** 2 + delta_phi(selected["z_phi"], selected["jpsi_phi"]) ** 2)
    selected["z_jpsi_dPhi"] = delta_phi(selected["z_phi"], selected["jpsi_phi"])
    selected["z_jpsi_dEta"] = selected["z_eta"] - selected["jpsi_eta"]
    selected["z_jpsi_dY"] = selected["z_rapidity"] - selected["jpsi_rapidity"]
    selected["pt_balance"] = np.abs(selected["z_pt"] - selected["jpsi_pt"]) / (selected["z_pt"] + selected["jpsi_pt"])
    return selected


def build_selected_daughter_branches(flat: Dict[str, np.ndarray], opts: List[Dict[str, object]], best_idx: np.ndarray) -> Dict[str, np.ndarray]:
    selected: Dict[str, np.ndarray] = {}

    charge_ordered_prefixes = ("jpsi_muP", "jpsi_muM", "z_muP", "z_muM")
    role_for_prefix = {
        "jpsi_muP": "jpsiP",
        "jpsi_muM": "jpsiM",
        "z_muP": "zP",
        "z_muM": "zM",
    }

    float_vars = (
        "pt", "eta", "phi",
        "pfRelIso03", "pfRelIso04",
        "trackAbsIso03", "trackRelIso03",
        "dB3D", "dxy", "dz", "normChi2",
    )
    int_vars = ("charge", "nValidHits", "nValidPixelHits")
    bool_vars = ("soft", "tight", "loose")

    # Keep the charge-labelled daughter branches for compatibility/debugging.
    for var in float_vars:
        for out_prefix in charge_ordered_prefixes:
            selected[f"{out_prefix}_{var}"] = choose_muon_branch(flat, opts, best_idx, var, role_for_prefix[out_prefix], np.float32)

    for var in int_vars:
        for out_prefix in charge_ordered_prefixes:
            selected[f"{out_prefix}_{var}"] = choose_muon_branch(flat, opts, best_idx, var, role_for_prefix[out_prefix], np.int32)

    for var in bool_vars:
        for out_prefix in charge_ordered_prefixes:
            selected[f"{out_prefix}_{var}"] = choose_muon_branch(flat, opts, best_idx, var, role_for_prefix[out_prefix], np.bool_)

    # Add pT-ordered daughter branches after the selected J/psi and Z assignment.
    # mu1 is the leading-pT daughter within that resonance; mu2 is subleading.
    all_vars = float_vars + int_vars + bool_vars
    resonance_pairs = (
        ("jpsi", "jpsi_muP", "jpsi_muM"),
        ("z", "z_muP", "z_muM"),
    )
    for resonance, first_charge_prefix, second_charge_prefix in resonance_pairs:
        first_is_leading = selected[f"{first_charge_prefix}_pt"] >= selected[f"{second_charge_prefix}_pt"]
        for var in all_vars:
            first_values = selected[f"{first_charge_prefix}_{var}"]
            second_values = selected[f"{second_charge_prefix}_{var}"]
            selected[f"{resonance}_mu1_{var}"] = np.where(first_is_leading, first_values, second_values)
            selected[f"{resonance}_mu2_{var}"] = np.where(first_is_leading, second_values, first_values)

    return selected

def choose_pairing(flat: Dict[str, np.ndarray], cfg: SampleConfig) -> Dict[str, np.ndarray]:
    """Choose and materialize the selected J/psi/Z assignment."""
    opts = pairing_options(flat)
    score_arr, valid_arr, vtx_arr = compute_pairing_arrays(flat, cfg, opts)
    best_idx = choose_best_pairing_indices(score_arr, valid_arr, vtx_arr)

    selected: Dict[str, np.ndarray] = {
        "pairing": np.asarray([int(opts[i]["code"]) for i in best_idx], dtype=np.int32),
        "pairing_massScore": score_arr[best_idx, np.arange(len(best_idx))],
        "pairing_vertexProduct": vtx_arr[best_idx, np.arange(len(best_idx))],
    }
    selected.update(build_selected_pair_branches(flat, opts, best_idx))
    selected.update(build_selected_daughter_branches(flat, opts, best_idx))
    return selected

def flatten_input_chunk(arrays: ak.Array) -> Dict[str, np.ndarray]:
    flat: Dict[str, np.ndarray] = {}
    for branch in required_branches():
        if branch in ("Run", "LumiBlock", "nPV"):
            flat[branch] = event_level_flatten(arrays, branch, dtype=np.uint64 if branch != "nPV" else np.int32)
        elif branch == "Event":
            flat[branch] = event_level_flatten(arrays, branch, dtype=np.uint64)
        elif branch.startswith("mu") and branch.endswith(("charge", "nValidHits", "nValidPixelHits")):
            flat[branch] = flatten_branch(arrays, branch, dtype=np.int32)
        elif branch.endswith(("soft", "tight", "loose", "passMuonTrigger")) or branch == "passMuonTrigger":
            flat[branch] = flatten_branch(arrays, branch, dtype=bool)
        else:
            flat[branch] = flatten_branch(arrays, branch, dtype=np.float32)
    return flat


def build_selection_mask(flat: Dict[str, np.ndarray], chosen: Dict[str, np.ndarray], cfg: SampleConfig) -> np.ndarray:
    topology_mask = broad_topology_mask(flat)

    # Trigger is required for data-derived samples only.
    trigger_mask = np.ones_like(flat["fourMu_mass"], dtype=bool)
    if cfg.name != "signal":
        trigger_mask &= flat["passMuonTrigger"]

    # Four basic muon cuts. Use the original charge-ordered muons, independent of pairing choice.
    muon_mask = np.ones_like(flat["fourMu_mass"], dtype=bool)
    for mu in ("muP1", "muM1", "muP2", "muM2"):
        muon_mask &= flat[f"{mu}_pt"] > MUON_PT_MIN
        muon_mask &= np.abs(flat[f"{mu}_eta"]) < MUON_ABS_ETA_MAX
        if REQUIRE_SOFT_MUONS:
            muon_mask &= flat[f"{mu}_soft"]

    pair_mass_mask = in_window(chosen["jpsi_mass"], cfg.jpsi_mass_window) & in_window(chosen["z_mass"], cfg.z_mass_window)

    vtx_mask = flat["fourMu_vtxProb"] > FOURMU_VTXPROB_MIN
    vtx_mask &= chosen["jpsi_vtxProb"] > PAIR_VTXPROB_MIN
    vtx_mask &= chosen["z_vtxProb"] > PAIR_VTXPROB_MIN

    pt_mask = chosen["jpsi_pt"] > DIMUON_PT_MIN
    pt_mask &= chosen["z_pt"] > DIMUON_PT_MIN
    pt_mask &= flat["fourMu_pt"] > FOURMU_PT_MIN

    fourmu_mask = fourmu_region_mask(flat["fourMu_mass"], cfg.fourmu_region)
    finite_mask = np.isfinite(chosen["pairing_massScore"]) & np.isfinite(flat["fourMu_mass"])

    return topology_mask & trigger_mask & muon_mask & pair_mass_mask & vtx_mask & pt_mask & fourmu_mask & finite_mask

def output_template() -> Dict[str, List]:
    keys = [
        "Run", "LumiBlock", "Event", "nPV", "label", "pairing",
        "fourMu_mass", "fourMu_pt", "fourMu_eta", "fourMu_phi", "fourMu_rapidity", "fourMu_vtxProb",
        "passMuonTrigger",
        "jpsi_mass", "jpsi_pt", "jpsi_eta", "jpsi_phi", "jpsi_rapidity", "jpsi_vtxProb", "jpsi_dR_mumu",
        "z_mass", "z_pt", "z_eta", "z_phi", "z_rapidity", "z_vtxProb", "z_dR_mumu",
        "pairing_massScore", "pairing_vertexProduct",
        "z_jpsi_dR", "z_jpsi_dPhi", "z_jpsi_dEta", "z_jpsi_dY", "pt_balance",
        "jpsi_trackIso03", "jpsi_trackIso04", "jpsi_relIso03", "jpsi_relIso04",
        "z_trackIso03", "z_trackIso04", "z_relIso03", "z_relIso04",
        "cosTheta_jpsiMuP", "cosTheta_zMuP", "phi_decayPlane_z_jpsi",
    ]
    for prefix in ("jpsi_muP", "jpsi_muM", "z_muP", "z_muM", "jpsi_mu1", "jpsi_mu2", "z_mu1", "z_mu2"):
        keys += [
            f"{prefix}_pt", f"{prefix}_eta", f"{prefix}_phi", f"{prefix}_charge",
            f"{prefix}_soft", f"{prefix}_loose", f"{prefix}_tight",
            f"{prefix}_pfRelIso03", f"{prefix}_pfRelIso04",
            f"{prefix}_trackAbsIso03", f"{prefix}_trackRelIso03",
            f"{prefix}_dB3D", f"{prefix}_dxy", f"{prefix}_dz", f"{prefix}_normChi2",
            f"{prefix}_nValidHits", f"{prefix}_nValidPixelHits",
        ]
    return {k: [] for k in keys}


def append_selected_rows(out: Dict[str, List], flat: Dict[str, np.ndarray], chosen: Dict[str, np.ndarray], indices: np.ndarray, cfg: SampleConfig) -> None:
    # Event identifiers and four-muon quantities
    for key in (
        "Run", "LumiBlock", "Event", "nPV", "passMuonTrigger",
        "fourMu_mass", "fourMu_pt", "fourMu_eta", "fourMu_phi", "fourMu_rapidity", "fourMu_vtxProb",
    ):
        out[key].extend(flat[key][indices].tolist())

    out["label"].extend([cfg.label] * len(indices))

    selected_keys = [k for k in out.keys() if k in chosen]
    for key in selected_keys:
        out[key].extend(chosen[key][indices].tolist())


def keep_best_candidate_per_event(out: Dict[str, List]) -> Dict[str, List]:
    """Keep one candidate per event within the already accumulated output.

    Ranking:
      1. lowest pairing_massScore
      2. highest fourMu_vtxProb when mass scores are close
      3. highest pairing_vertexProduct
    """
    n = len(out["Event"])
    if n == 0:
        return out

    events = np.asarray(out["Event"], dtype=np.uint64)
    runs = np.asarray(out["Run"], dtype=np.uint64)
    lumis = np.asarray(out["LumiBlock"], dtype=np.uint64)
    mass_score_arr = np.asarray(out["pairing_massScore"], dtype=float)
    four_vtx_arr = np.asarray(out["fourMu_vtxProb"], dtype=float)
    pair_vtx_arr = np.asarray(out["pairing_vertexProduct"], dtype=float)

    best_by_event: Dict[Tuple[int, int, int], int] = {}

    def is_better(i: int, j: int) -> bool:
        """Return True if i is better than current j."""
        if mass_score_arr[i] < mass_score_arr[j] - MASS_SCORE_CLOSE:
            return True
        if mass_score_arr[i] > mass_score_arr[j] + MASS_SCORE_CLOSE:
            return False
        if four_vtx_arr[i] > four_vtx_arr[j] + 1e-12:
            return True
        if four_vtx_arr[i] < four_vtx_arr[j] - 1e-12:
            return False
        return pair_vtx_arr[i] > pair_vtx_arr[j]

    for i in range(n):
        key = (int(runs[i]), int(lumis[i]), int(events[i]))
        if key not in best_by_event or is_better(i, best_by_event[key]):
            best_by_event[key] = i

    keep = sorted(best_by_event.values())
    return {key: [values[i] for i in keep] for key, values in out.items()}


def convert_for_uproot(out: Dict[str, List]) -> Dict[str, np.ndarray]:
    arrays: Dict[str, np.ndarray] = {}
    int_keys = {"Run", "LumiBlock", "Event", "nPV", "label", "pairing"}
    for key, values in out.items():
        is_bool = key == "passMuonTrigger" or key.endswith(("_soft", "_tight", "_loose"))
        is_int = key in int_keys or key.endswith(("_charge", "_nValidHits", "_nValidPixelHits"))
        if key in {"Run", "LumiBlock", "Event"}:
            arrays[key] = np.asarray(values, dtype=np.uint64)
        elif is_int:
            arrays[key] = np.asarray(values, dtype=np.int32)
        elif is_bool:
            arrays[key] = np.asarray(values, dtype=np.bool_)
        else:
            arrays[key] = np.asarray(values, dtype=np.float32)
    return arrays


def is_remote_input(path: str) -> bool:
    return path.startswith(("root://", "file:"))


def ensure_input_exists(path: str) -> None:
    input_path = Path(path)
    if not input_path.exists() and not is_remote_input(path):
        raise FileNotFoundError(f"Input file not found: {path}")


def input_tree_path(path: str) -> str:
    return f"{path}:{TREE_NAME}"


def count_tree_entries(path: str) -> Optional[int]:
    try:
        with uproot.open(input_tree_path(path)) as tree_for_count:
            return int(tree_for_count.num_entries)
    except Exception:
        return None


def iter_input_chunks(path: str):
    return uproot.iterate(input_tree_path(path), expressions=required_branches(), step_size=STEP_SIZE, library="ak")


def reduce_candidates_and_write(cfg: SampleConfig, out: Dict[str, List], keep_one_candidate_per_event: bool) -> Tuple[int, int]:
    before = len(out["Event"])
    if keep_one_candidate_per_event:
        out = keep_best_candidate_per_event(out)
    after = len(out["Event"])

    output_path = Path(cfg.output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(output_path) as fout:
        fout[TREE_NAME] = convert_for_uproot(out)

    return before, after


def process_sample(cfg: SampleConfig, keep_one_candidate_per_event: bool) -> None:
    ensure_input_exists(cfg.input_path)

    print(f"\nProcessing {cfg.name}")
    print(f"  input : {cfg.input_path}")
    print(f"  output: {cfg.output_path}")
    print(f"  J/psi window: {cfg.jpsi_mass_window}")
    print(f"  Z window    : {cfg.z_mass_window}")
    print(f"  fourMu mode : {cfg.fourmu_region}")
    print(f"  analysis window: ({ANALYSIS_LOW}, {ANALYSIS_HIGH})")
    print(f"  mask window    : ({MASK_LOW}, {MASK_HIGH})")
    print(f"  step size      : {STEP_SIZE}")
    print("  pairing mode: old_compatible")
    print(f"  one candidate per event: {keep_one_candidate_per_event}")

    out = output_template()
    total_flat = 0
    total_selected_before_reduction = 0
    total_events_seen = 0
    chunk_index = 0
    t0 = time()

    total_entries = count_tree_entries(cfg.input_path)

    for arrays in iter_input_chunks(cfg.input_path):
        chunk_index += 1
        n_events_chunk = len(arrays["fourMu_mass"])
        total_events_seen += n_events_chunk

        flat = flatten_input_chunk(arrays)
        n_flat_chunk = len(flat["fourMu_mass"])
        total_flat += n_flat_chunk

        if n_flat_chunk == 0:
            if PRINT_PROGRESS and (chunk_index % PROGRESS_EVERY_CHUNKS == 0):
                denom = f"/{total_entries}" if total_entries is not None else ""
                print(f"    chunk {chunk_index:5d}: events {total_events_seen}{denom}, flat +0, selected +0, kept so far {len(out['Event'])}", flush=True)
            continue

        chosen = choose_pairing(flat, cfg)
        mask = build_selection_mask(flat, chosen, cfg)
        indices = np.nonzero(mask)[0]

        total_selected_before_reduction += len(indices)

        if len(indices) > 0:
            append_selected_rows(out, flat, chosen, indices, cfg)

        if PRINT_PROGRESS and (chunk_index % PROGRESS_EVERY_CHUNKS == 0):
            elapsed = max(time() - t0, 1e-9)
            rate = total_events_seen / elapsed
            if total_entries is not None and total_entries > 0:
                frac = 100.0 * total_events_seen / total_entries
                denom = f"/{total_entries} ({frac:5.1f}%)"
            else:
                denom = ""
            print(
                f"    chunk {chunk_index:5d}: events {total_events_seen}{denom}, "
                f"flat +{n_flat_chunk}, selected +{len(indices)}, "
                f"selected total {total_selected_before_reduction}, kept so far {len(out['Event'])}, "
                f"rate {rate:.1f} events/s",
                flush=True,
            )

    n_before, n_after = reduce_candidates_and_write(cfg, out, keep_one_candidate_per_event)

    elapsed = time() - t0
    print(f"  elapsed seconds              : {elapsed:.1f}")
    print(f"  input events processed       : {total_events_seen}")
    print(f"  flattened candidates         : {total_flat}")
    print(f"  selected before reduction: {total_selected_before_reduction}")
    print(f"  written before reduction : {n_before}")
    print(f"  written after reduction  : {n_after}")


def make_sample_config(
    name: str,
    input_path: str,
    output_name: str,
    jpsi_mass_window: Tuple[float, float],
    z_mass_window: Tuple[float, float],
    fourmu_region: str,
    label: int,
    outdir: Path,
) -> SampleConfig:
    return SampleConfig(
        name=name,
        input_path=input_path,
        output_path=str(outdir / output_name),
        jpsi_mass_window=jpsi_mass_window,
        z_mass_window=z_mass_window,
        fourmu_region=fourmu_region,
        label=label,
    )


def build_signal_config(outdir: Path) -> Optional[SampleConfig]:
    if not MAKE_SIGNAL:
        return None
    return make_sample_config(
        "signal",
        SIGNAL_PRESELECTION,
        OUTPUTS["signal"],
        SIGNAL_JPSI_MASS,
        SIGNAL_Z_MASS,
        SIGNAL_FOURMU_REGION,
        1,
        outdir,
    )


def build_data_configs(outdir: Path) -> List[SampleConfig]:
    configs: List[SampleConfig] = []
    if MAKE_BACKGROUND:
        configs.append(make_sample_config(
            "background",
            DATA_PRESELECTION,
            OUTPUTS["background"],
            DATA_JPSI_MASS,
            DATA_Z_MASS,
            BACKGROUND_FOURMU_REGION,
            0,
            outdir,
        ))
    if MAKE_FINAL_BLINDED:
        configs.append(make_sample_config(
            "final_blinded",
            DATA_PRESELECTION,
            OUTPUTS["final_blinded"],
            DATA_JPSI_MASS,
            DATA_Z_MASS,
            FINAL_BLINDED_FOURMU_REGION,
            -1,
            outdir,
        ))
    if MAKE_FINAL_UNBLINDED:
        configs.append(make_sample_config(
            "final_unblinded",
            DATA_PRESELECTION,
            OUTPUTS["final_unblinded"],
            DATA_JPSI_MASS,
            DATA_Z_MASS,
            FINAL_UNBLINDED_FOURMU_REGION,
            -1,
            outdir,
        ))
    return configs


def validate_data_configs(configs: List[SampleConfig]) -> None:
    """All data outputs are processed in one pass.

    This assumes the data-derived outputs use the same J/psi/Z mass windows,
    which is true for the current background/final definitions. The four-muon
    region may differ output-by-output.
    """
    if not configs:
        return

    jpsi0 = configs[0].jpsi_mass_window
    z0 = configs[0].z_mass_window
    for cfg in configs[1:]:
        if cfg.jpsi_mass_window != jpsi0 or cfg.z_mass_window != z0:
            raise ValueError(
                "One-pass data processing requires all data configs to use the "
                "same J/psi/Z mass windows. Split the pass if you change this."
            )


def process_data_once(configs: List[SampleConfig], keep_one_candidate_per_event: bool) -> None:
    """Read the data preselection file once and fill all requested data outputs.

    This avoids reading the same large Run-2 data file separately for background,
    final_blinded, and final_unblinded samples.
    """
    if not configs:
        return

    validate_data_configs(configs)

    ensure_input_exists(DATA_PRESELECTION)

    print("\nProcessing data once for:")
    for cfg in configs:
        print(f"  {cfg.name:15s} -> {cfg.output_path}")
    print(f"  input : {DATA_PRESELECTION}")
    print(f"  J/psi window: {configs[0].jpsi_mass_window}")
    print(f"  Z window    : {configs[0].z_mass_window}")
    print(f"  analysis window: ({ANALYSIS_LOW}, {ANALYSIS_HIGH})")
    print(f"  mask window    : ({MASK_LOW}, {MASK_HIGH})")
    print(f"  background region: {BACKGROUND_FOURMU_REGION if MAKE_BACKGROUND else 'disabled'}")
    print(f"  optional background mass range: {BACKGROUND_FOURMU_RANGE}")
    print(f"  step size      : {STEP_SIZE}")
    print("  pairing mode   : old_compatible")
    print(f"  one candidate per event: {keep_one_candidate_per_event}")

    outputs = {cfg.name: output_template() for cfg in configs}
    total_flat = 0
    totals_selected = {cfg.name: 0 for cfg in configs}
    total_events_seen = 0
    chunk_index = 0
    t0 = time()

    total_entries = count_tree_entries(DATA_PRESELECTION)

    # Pairing choice is common across the data outputs because J/psi/Z windows are common.
    pair_cfg = configs[0]

    for arrays in iter_input_chunks(DATA_PRESELECTION):
        chunk_index += 1
        n_events_chunk = len(arrays["fourMu_mass"])
        total_events_seen += n_events_chunk

        flat = flatten_input_chunk(arrays)
        n_flat_chunk = len(flat["fourMu_mass"])
        total_flat += n_flat_chunk

        if n_flat_chunk == 0:
            if PRINT_PROGRESS and (chunk_index % PROGRESS_EVERY_CHUNKS == 0):
                denom = f"/{total_entries}" if total_entries is not None else ""
                print(f"    chunk {chunk_index:5d}: events {total_events_seen}{denom}, flat +0", flush=True)
            continue

        chosen = choose_pairing(flat, pair_cfg)

        selected_this_chunk = {}
        for cfg in configs:
            mask = build_selection_mask(flat, chosen, cfg)
            indices = np.nonzero(mask)[0]
            selected_this_chunk[cfg.name] = len(indices)
            totals_selected[cfg.name] += len(indices)

            if len(indices) > 0:
                append_selected_rows(outputs[cfg.name], flat, chosen, indices, cfg)

        if PRINT_PROGRESS and (chunk_index % PROGRESS_EVERY_CHUNKS == 0):
            elapsed = max(time() - t0, 1e-9)
            rate = total_events_seen / elapsed
            if total_entries is not None and total_entries > 0:
                frac = 100.0 * total_events_seen / total_entries
                denom = f"/{total_entries} ({frac:5.1f}%)"
            else:
                denom = ""

            per_chunk = ", ".join(
                f"{name} +{selected_this_chunk[name]} total {totals_selected[name]}"
                for name in selected_this_chunk
            )
            kept_so_far = ", ".join(
                f"{name} {len(outputs[name]['Event'])}"
                for name in outputs
            )
            print(
                f"    chunk {chunk_index:5d}: events {total_events_seen}{denom}, "
                f"flat +{n_flat_chunk}, {per_chunk}, kept [{kept_so_far}], "
                f"rate {rate:.1f} events/s",
                flush=True,
            )

    for cfg in configs:
        before, after = reduce_candidates_and_write(cfg, outputs[cfg.name], keep_one_candidate_per_event)

        print(f"\nFinished {cfg.name}")
        print(f"  output: {cfg.output_path}")
        print(f"  selected before reduction: {totals_selected[cfg.name]}")
        print(f"  written before reduction : {before}")
        print(f"  written after reduction  : {after}")

    elapsed = max(time() - t0, 1e-9)
    print("\nFinished one-pass data processing")
    print(f"  flattened candidates: {total_flat}")
    print(f"  runtime             : {elapsed:.1f} s")
    print(f"  rate                : {total_events_seen / elapsed:.1f} events/s")


def main() -> None:
    if not CANDIDATE_VARIANTS:
        raise RuntimeError("CANDIDATE_VARIANTS is empty.")

    if not (MAKE_SIGNAL or MAKE_BACKGROUND or MAKE_FINAL_BLINDED or MAKE_FINAL_UNBLINDED):
        raise RuntimeError("No outputs requested. Enable at least one MAKE_* switch.")

    base_outdir = Path(OUTDIR)

    for variant_name, keep_one_candidate_per_event in CANDIDATE_VARIANTS:
        variant_outdir = base_outdir / variant_name
        print("\n" + "=" * 80)
        print(f"Candidate variant: {variant_name}")
        print(f"  output directory: {variant_outdir}")
        print(f"  one candidate per event: {keep_one_candidate_per_event}")
        print("=" * 80)

        signal_cfg = build_signal_config(variant_outdir)
        data_cfgs = build_data_configs(variant_outdir)

        if signal_cfg is not None:
            process_sample(signal_cfg, keep_one_candidate_per_event=keep_one_candidate_per_event)

        # The data file is read once per candidate variant for background/final outputs.
        process_data_once(data_cfgs, keep_one_candidate_per_event=keep_one_candidate_per_event)


if __name__ == "__main__":
    main()
