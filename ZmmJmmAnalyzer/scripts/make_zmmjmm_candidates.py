#!/usr/bin/env python3
"""
Build flat Z(mm)J/psi(mm) candidate trees from the cleaned miniAODmmmm preselection ntuples.

Design choices:
  * One script with hardcoded paths/switches at the top.
  * Reads two preselection files: merged Run-2 data and signal MC.
  * Produces signal once, and processes the data file once to make background/final outputs.
  * Does not write string branches, for compatibility with older uproot writers.
  * Can run in old-compatible reproduction mode or optimized score mode.
  * Applies the combined old preselection+selection cuts in Task 2, while Task 1 remains broad.
  * Can keep all candidates or deduplicate to one candidate per event.

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

MAKE_SIGNAL = True
MAKE_BACKGROUND = True
MAKE_FINAL_BLINDED = True
MAKE_FINAL_UNBLINDED = True  # keep False until you are ready to inspect/unblind

DATA_PRESELECTION = "preselection/skimmed_with_ROOT_v2/skimmed_TTree_13TeV_mmmm_UL_Run2.root"
SIGNAL_PRESELECTION = "preselection/run2_mc/zmmjmm/zmmjmm_mc_v3_2018_ss.root"
OUTDIR = "selection"
TREE_NAME = "ntuple"

# Data analysis window and blinded/masked Higgs window
# The data analysis region is the broad four-muon mass range used for candidate production.
# The mask window is the narrow Higgs window to exclude from blinded data/background training.
ANALYSIS_LOW = 112.0
ANALYSIS_HIGH = 162.0
MASK_LOW = 120.0
MASK_HIGH = 130.0

# -----------------------------------------------------------------------------
# Selection mode switches
# -----------------------------------------------------------------------------
# Start with this old-compatible mode to reproduce the previous AnalysisWithIso_*
# candidate yields. Once that is understood, switch PAIRING_MODE to "score" and
# DEDUPLICATE to True for the cleaner ML production.
PAIRING_MODE = "old_compatible"  # "old_compatible" or "score"
DEDUPLICATE = False              # old scripts wrote candidates, not one/event

# Broad topology layer: reproduces the old miniAOD-level requirement that one
# dimuon is low-mass and the other is Z-like.  This is intentionally loose.
APPLY_OLD_BROAD_TOPOLOGY = True
LOW_MASS_WINDOW = (0.0, 12.0)
BROAD_Z_WINDOW = (70.0, 110.0)
BROAD_PAIR_VTXPROB_MIN = 0.001
BROAD_FOURMU_VTXPROB_MIN = 0.001

# Final physics-selection layer. These are the old post-selection style cuts.
APPLY_TRIGGER_FOR_DATA = True     # applies to background/final samples only
APPLY_SOFT_MUON = True
APPLY_MUON_KINEMATICS = True
APPLY_VERTEX_CUTS = True
APPLY_DIMUON_PT_CUTS = True
APPLY_FOURMU_PT_CUT = True
APPLY_ISOLATION = False           # keep off until the old isolation definition is confirmed

MUON_PT_MIN = 3.0
MUON_ABS_ETA_MAX = 2.4
FOURMU_VTXPROB_MIN = 0.01
PAIR_VTXPROB_MIN = 0.01
DIMUON_PT_MIN = 5.0
FOURMU_PT_MIN = 5.0

# Placeholder if APPLY_ISOLATION=True later.  This uses the maximum per-muon
# pfRelIso03 among the selected four muons. Do not enable until validated.
MAX_MUON_RELISO03 = 0.35

# Pairing score constants. These are not detector resolutions; they are practical
# relative weights for choosing the Z/Jpsi assignment.
M_JPSI = 3.0969
M_Z = 91.1876
SIGMA_JPSI_SCORE = 0.10
SIGMA_Z_SCORE = 10.0
MASS_SCORE_CLOSE = 0.5

# Sample-specific mass windows. Pairing choice and selection use these windows.
CUTS = {
    "signal": {
        "jpsi_mass": (3.0, 3.2),
        "z_mass": (80.0, 100.0),
        "fourmu_region": "analysis_window",  # old signal selection used 112--162
        "label": 1,
    },
    "background": {
        "jpsi_mass": (2.8, 3.4),
        "z_mass": (70.0, 110.0),
        "fourmu_region": "outside_mask",  # no fourMu mass limit except the 120--130 mask
        "label": 0,
    },
    "final_blinded": {
        "jpsi_mass": (2.8, 3.4),
        "z_mass": (70.0, 110.0),
        "fourmu_region": "analysis_sideband",
        "label": -1,
    },
    "final_unblinded": {
        "jpsi_mass": (2.8, 3.4),
        "z_mass": (70.0, 110.0),
        "fourmu_region": "analysis_window",
        "label": -1,
    },
}

# Chunk size for uproot.iterate. Increase if memory is fine.
STEP_SIZE = "100 MB"

# Progress printing for large files.
PRINT_PROGRESS = True
PROGRESS_EVERY_CHUNKS = 1

# Output file names
OUTPUTS = {
    "signal": "signal_candidates.root",
    "background": "background_candidates.root",
    "final_blinded": "final_blinded_candidates.root",
    "final_unblinded": "final_unblinded_candidates.root",
}

# =============================================================================
# Internal helpers
# =============================================================================

SCALAR_DEFAULT_FLOAT = -999.0
SCALAR_DEFAULT_INT = -999


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
      outside_mask      : all masses excluding 120--130 GeV
      outside_analysis  : outside 112--162 GeV
    """
    analysis = (fourmu_mass > ANALYSIS_LOW) & (fourmu_mass < ANALYSIS_HIGH)
    mask = (fourmu_mass > MASK_LOW) & (fourmu_mass < MASK_HIGH)

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


def event_level_flatten(arrays: ak.Array, branch: str, n_cands_per_event: ak.Array, dtype=None) -> np.ndarray:
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
            f"{mu}_pfRelIso03", f"{mu}_pfRelIso04", f"{mu}_dB3D",
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
    """Old loose topology: one low-mass dimuon and one Z-like dimuon.

    This is applied on either disjoint pair group, with either pair allowed to be
    the low-mass object.  It is deliberately broader than the final J/psi mass
    selection.
    """
    if not APPLY_OLD_BROAD_TOPOLOGY:
        return np.ones_like(flat["fourMu_mass"], dtype=bool)

    p12_low = in_window(flat["pair12_mass"], LOW_MASS_WINDOW)
    p34_low = in_window(flat["pair34_mass"], LOW_MASS_WINDOW)
    p14_low = in_window(flat["pair14_mass"], LOW_MASS_WINDOW)
    p23_low = in_window(flat["pair23_mass"], LOW_MASS_WINDOW)

    p12_z = in_window(flat["pair12_mass"], BROAD_Z_WINDOW)
    p34_z = in_window(flat["pair34_mass"], BROAD_Z_WINDOW)
    p14_z = in_window(flat["pair14_mass"], BROAD_Z_WINDOW)
    p23_z = in_window(flat["pair23_mass"], BROAD_Z_WINDOW)

    group_a_mass = (p12_low & p34_z) | (p34_low & p12_z)
    group_b_mass = (p14_low & p23_z) | (p23_low & p14_z)

    group_a_vtx = (flat["pair12_vtxProb"] > BROAD_PAIR_VTXPROB_MIN) & (flat["pair34_vtxProb"] > BROAD_PAIR_VTXPROB_MIN)
    group_b_vtx = (flat["pair14_vtxProb"] > BROAD_PAIR_VTXPROB_MIN) & (flat["pair23_vtxProb"] > BROAD_PAIR_VTXPROB_MIN)
    four_vtx = flat["fourMu_vtxProb"] > BROAD_FOURMU_VTXPROB_MIN

    return four_vtx & ((group_a_mass & group_a_vtx) | (group_b_mass & group_b_vtx))


def choose_pairing(flat: Dict[str, np.ndarray], cfg: SampleConfig) -> Dict[str, np.ndarray]:
    """Choose the selected J/psi/Z assignment for each flattened 4mu candidate.

    In old_compatible mode, valid final mass-window assignments are preferred,
    and vertex-product breaks near mass-score ties. This preserves the old idea
    of choosing between disjoint pair groups using vertex quality, while also
    allowing either pair in a group to be the J/psi.

    In score mode, the lowest J/psi/Z mass score dominates, with vertex-product
    only as the tie-breaker.
    """
    opts = pairing_options(flat)
    n = len(flat["fourMu_mass"])

    scores = []
    valid = []
    vtx_products = []
    for opt in opts:
        jp = str(opt["jpsi"])
        zz = str(opt["z"])
        scores.append(mass_score(flat[f"{jp}_mass"], flat[f"{zz}_mass"]))
        valid.append(in_window(flat[f"{jp}_mass"], cfg.jpsi_mass_window) & in_window(flat[f"{zz}_mass"], cfg.z_mass_window))
        vtx_products.append(safe_vertex_product(flat[f"{jp}_vtxProb"], flat[f"{zz}_vtxProb"]))

    score_arr = np.vstack(scores)          # shape: (4, n)
    valid_arr = np.vstack(valid)           # shape: (4, n)
    vtx_arr = np.vstack(vtx_products)      # shape: (4, n)

    if PAIRING_MODE == "score":
        # Invalid assignments get a large penalty, but are still selectable if no
        # assignment is valid; this keeps diagnostics possible.
        penalized = score_arr + np.where(valid_arr, 0.0, 1e6)
        no_valid = ~np.any(valid_arr, axis=0)
        penalized[:, no_valid] = score_arr[:, no_valid]
        best_idx = np.argmin(penalized, axis=0)
    elif PAIRING_MODE == "old_compatible":
        # Prefer assignments passing the final mass windows. Among valid options,
        # lower mass score wins unless close, then higher vertex-product wins.
        best_idx = np.zeros(n, dtype=np.int32)
        for i in range(n):
            candidates = np.nonzero(valid_arr[:, i])[0]
            if len(candidates) == 0:
                candidates = np.arange(len(opts))

            best = int(candidates[0])
            for cand in candidates[1:]:
                cand = int(cand)
                if score_arr[cand, i] < score_arr[best, i] - MASS_SCORE_CLOSE:
                    best = cand
                elif abs(score_arr[cand, i] - score_arr[best, i]) <= MASS_SCORE_CLOSE:
                    if vtx_arr[cand, i] > vtx_arr[best, i]:
                        best = cand
            best_idx[i] = best
    else:
        raise ValueError(f"Unknown PAIRING_MODE: {PAIRING_MODE}")

    selected: Dict[str, np.ndarray] = {
        "pairing": np.asarray([int(opts[i]["code"]) for i in best_idx], dtype=np.int32),
    }

    def choose_from_pairs(suffix: str, kind: str) -> np.ndarray:
        vals = np.empty(n, dtype=np.float32)
        for idx, opt in enumerate(opts):
            pair = str(opt[kind])
            mask = best_idx == idx
            vals[mask] = flat[f"{pair}_{suffix}"][mask]
        return vals

    for suffix in ("mass", "pt", "eta", "phi", "rapidity", "vtxProb"):
        selected[f"jpsi_{suffix}"] = choose_from_pairs(suffix, "jpsi")
        selected[f"z_{suffix}"] = choose_from_pairs(suffix, "z")

    selected["jpsi_dR_mumu"] = choose_from_pairs("dR", "jpsi")
    selected["z_dR_mumu"] = choose_from_pairs("dR", "z")
    for suffix in ("trackIso03", "trackIso04", "relIso03", "relIso04", "cosThetaMu"):
        selected[f"jpsi_{suffix}"] = choose_from_pairs(suffix, "jpsi")
        selected[f"z_{suffix}"] = choose_from_pairs(suffix, "z")

    selected["cosTheta_jpsiMu"] = selected.pop("jpsi_cosThetaMu")
    selected["cosTheta_zMu"] = selected.pop("z_cosThetaMu")

    phi_plane = np.empty(n, dtype=np.float32)
    for idx, opt in enumerate(opts):
        mask = best_idx == idx
        phi_plane[mask] = flat[str(opt["plane"])][mask]
    selected["phi_decayPlane"] = phi_plane

    selected["pairing_massScore"] = score_arr[best_idx, np.arange(n)]
    selected["pairing_vertexProduct"] = vtx_arr[best_idx, np.arange(n)]

    selected["dR_jpsi_z"] = np.sqrt((selected["jpsi_eta"] - selected["z_eta"]) ** 2 + delta_phi(selected["jpsi_phi"], selected["z_phi"]) ** 2)
    selected["dPhi_jpsi_z"] = delta_phi(selected["jpsi_phi"], selected["z_phi"])
    selected["dEta_jpsi_z"] = selected["jpsi_eta"] - selected["z_eta"]
    selected["dY_jpsi_z"] = selected["jpsi_rapidity"] - selected["z_rapidity"]

    # Daughter assignment for the selected Jpsi/Z option.
    for var in ("pt", "eta", "phi", "pfRelIso03", "pfRelIso04", "dB3D", "dxy", "dz", "normChi2"):
        for out_prefix, role in (("jpsi_muP", "jpsiP"), ("jpsi_muM", "jpsiM"), ("z_muP", "zP"), ("z_muM", "zM")):
            vals = np.empty(n, dtype=np.float32)
            for idx, opt in enumerate(opts):
                mu = str(opt[role])
                mask = best_idx == idx
                vals[mask] = flat[f"{mu}_{var}"][mask]
            selected[f"{out_prefix}_{var}"] = vals

    for var in ("nValidHits", "nValidPixelHits"):
        for out_prefix, role in (("jpsi_muP", "jpsiP"), ("jpsi_muM", "jpsiM"), ("z_muP", "zP"), ("z_muM", "zM")):
            vals = np.empty(n, dtype=np.int32)
            for idx, opt in enumerate(opts):
                mu = str(opt[role])
                mask = best_idx == idx
                vals[mask] = flat[f"{mu}_{var}"][mask]
            selected[f"{out_prefix}_{var}"] = vals

    return selected

def flatten_input_chunk(arrays: ak.Array) -> Dict[str, np.ndarray]:
    n_cands = ak.num(arrays["fourMu_mass"], axis=1)
    flat: Dict[str, np.ndarray] = {}
    for branch in required_branches():
        if branch in ("Run", "LumiBlock", "nPV"):
            flat[branch] = event_level_flatten(arrays, branch, n_cands, dtype=np.uint64 if branch != "nPV" else np.int32)
        elif branch == "Event":
            flat[branch] = event_level_flatten(arrays, branch, n_cands, dtype=np.uint64)
        elif branch.startswith("mu") and branch.endswith(("charge", "nValidHits", "nValidPixelHits")):
            flat[branch] = flatten_branch(arrays, branch, dtype=np.int32)
        elif branch.endswith(("soft", "tight", "loose", "passMuonTrigger")) or branch == "passMuonTrigger":
            flat[branch] = flatten_branch(arrays, branch, dtype=bool)
        else:
            flat[branch] = flatten_branch(arrays, branch, dtype=np.float32)
    return flat


def build_selection_mask(flat: Dict[str, np.ndarray], chosen: Dict[str, np.ndarray], cfg: SampleConfig) -> np.ndarray:
    # Broad old-preselection-equivalent topology, if enabled.
    topology_mask = broad_topology_mask(flat)

    # Trigger: apply only to data-derived samples, not signal MC, unless you choose otherwise.
    trigger_mask = np.ones_like(flat["fourMu_mass"], dtype=bool)
    if APPLY_TRIGGER_FOR_DATA and cfg.name != "signal":
        trigger_mask &= flat["passMuonTrigger"]

    # Four basic muon cuts. Use the original charge-ordered muons, independent of pairing choice.
    muon_mask = np.ones_like(flat["fourMu_mass"], dtype=bool)
    if APPLY_MUON_KINEMATICS or APPLY_SOFT_MUON:
        for mu in ("muP1", "muM1", "muP2", "muM2"):
            if APPLY_MUON_KINEMATICS:
                muon_mask &= flat[f"{mu}_pt"] > MUON_PT_MIN
                muon_mask &= np.abs(flat[f"{mu}_eta"]) < MUON_ABS_ETA_MAX
            if APPLY_SOFT_MUON:
                muon_mask &= flat[f"{mu}_soft"]

    pair_mass_mask = in_window(chosen["jpsi_mass"], cfg.jpsi_mass_window) & in_window(chosen["z_mass"], cfg.z_mass_window)

    vtx_mask = np.ones_like(flat["fourMu_mass"], dtype=bool)
    if APPLY_VERTEX_CUTS:
        vtx_mask &= flat["fourMu_vtxProb"] > FOURMU_VTXPROB_MIN
        vtx_mask &= chosen["jpsi_vtxProb"] > PAIR_VTXPROB_MIN
        vtx_mask &= chosen["z_vtxProb"] > PAIR_VTXPROB_MIN

    pt_mask = np.ones_like(flat["fourMu_mass"], dtype=bool)
    if APPLY_DIMUON_PT_CUTS:
        pt_mask &= chosen["jpsi_pt"] > DIMUON_PT_MIN
        pt_mask &= chosen["z_pt"] > DIMUON_PT_MIN
    if APPLY_FOURMU_PT_CUT:
        pt_mask &= flat["fourMu_pt"] > FOURMU_PT_MIN

    iso_mask = np.ones_like(flat["fourMu_mass"], dtype=bool)
    if APPLY_ISOLATION:
        max_iso = np.maximum.reduce([
            flat["muP1_pfRelIso03"], flat["muM1_pfRelIso03"],
            flat["muP2_pfRelIso03"], flat["muM2_pfRelIso03"],
        ])
        iso_mask &= max_iso < MAX_MUON_RELISO03

    fourmu_mask = fourmu_region_mask(flat["fourMu_mass"], cfg.fourmu_region)
    finite_mask = np.isfinite(chosen["pairing_massScore"]) & np.isfinite(flat["fourMu_mass"])

    return topology_mask & trigger_mask & muon_mask & pair_mass_mask & vtx_mask & pt_mask & iso_mask & fourmu_mask & finite_mask

def output_template() -> Dict[str, List]:
    keys = [
        "run", "lumi", "event", "nPV", "label", "pairing",
        "fourMu_mass", "fourMu_pt", "fourMu_eta", "fourMu_phi", "fourMu_rapidity", "fourMu_vtxProb",
        "passMuonTrigger",
        "jpsi_mass", "jpsi_pt", "jpsi_eta", "jpsi_phi", "jpsi_rapidity", "jpsi_vtxProb",
        "z_mass", "z_pt", "z_eta", "z_phi", "z_rapidity", "z_vtxProb",
        "pairing_massScore", "pairing_vertexProduct",
        "dR_mumu_jpsi", "dR_mumu_z", "dR_jpsi_z", "dPhi_jpsi_z", "dEta_jpsi_z", "dY_jpsi_z",
        "jpsi_trackIso03", "jpsi_trackIso04", "jpsi_relIso03", "jpsi_relIso04",
        "z_trackIso03", "z_trackIso04", "z_relIso03", "z_relIso04",
        "cosTheta_jpsiMu", "cosTheta_zMu", "phi_decayPlane",
    ]
    for prefix in ("jpsi_muP", "jpsi_muM", "z_muP", "z_muM"):
        keys += [
            f"{prefix}_pt", f"{prefix}_eta", f"{prefix}_phi",
            f"{prefix}_pfRelIso03", f"{prefix}_pfRelIso04", f"{prefix}_dB3D",
            f"{prefix}_dxy", f"{prefix}_dz", f"{prefix}_normChi2",
            f"{prefix}_nValidHits", f"{prefix}_nValidPixelHits",
        ]
    return {k: [] for k in keys}


def append_selected_rows(out: Dict[str, List], flat: Dict[str, np.ndarray], chosen: Dict[str, np.ndarray], indices: np.ndarray, cfg: SampleConfig) -> None:
    # Event identifiers and four-muon quantities
    mapping = {
        "run": flat["Run"],
        "lumi": flat["LumiBlock"],
        "event": flat["Event"],
        "nPV": flat["nPV"],
        "passMuonTrigger": flat["passMuonTrigger"],
        "fourMu_mass": flat["fourMu_mass"],
        "fourMu_pt": flat["fourMu_pt"],
        "fourMu_eta": flat["fourMu_eta"],
        "fourMu_phi": flat["fourMu_phi"],
        "fourMu_rapidity": flat["fourMu_rapidity"],
        "fourMu_vtxProb": flat["fourMu_vtxProb"],
    }
    for key, values in mapping.items():
        out[key].extend(values[indices].tolist())

    out["label"].extend([cfg.label] * len(indices))

    chosen_mapping_keys = [k for k in out.keys() if k in chosen]
    for key in chosen_mapping_keys:
        out[key].extend(chosen[key][indices].tolist())

    # Rename chosen keys to stable output names where needed
    rename = {
        "jpsi_dR_mumu": "dR_mumu_jpsi",
        "z_dR_mumu": "dR_mumu_z",
    }
    for src, dst in rename.items():
        out[dst].extend(chosen[src][indices].tolist())


def deduplicate_one_chunk(out: Dict[str, List]) -> Dict[str, List]:
    """Keep one candidate per event within the already accumulated output.

    Ranking:
      1. lowest pairing_massScore
      2. highest fourMu_vtxProb when mass scores are close
      3. highest pairing_vertexProduct
    """
    n = len(out["event"])
    if n == 0:
        return out

    events = np.asarray(out["event"], dtype=np.uint64)
    runs = np.asarray(out["run"], dtype=np.uint64)
    lumis = np.asarray(out["lumi"], dtype=np.uint64)
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
    int_keys = {"run", "lumi", "event", "nPV", "label", "pairing"}
    bool_keys = {"passMuonTrigger"}
    for key, values in out.items():
        if key in int_keys:
            # ROOT handles uint64 fine for event ids; use signed int for smaller labels/pairing.
            if key in {"run", "lumi", "event"}:
                arrays[key] = np.asarray(values, dtype=np.uint64)
            else:
                arrays[key] = np.asarray(values, dtype=np.int32)
        elif key in bool_keys:
            arrays[key] = np.asarray(values, dtype=np.bool_)
        else:
            arrays[key] = np.asarray(values, dtype=np.float32)
    return arrays


def process_sample(cfg: SampleConfig) -> None:
    input_path = Path(cfg.input_path)
    if not input_path.exists() and not str(input_path).startswith(("root://", "file:")):
        raise FileNotFoundError(f"Input file not found: {cfg.input_path}")

    print(f"\nProcessing {cfg.name}")
    print(f"  input : {cfg.input_path}")
    print(f"  output: {cfg.output_path}")
    print(f"  J/psi window: {cfg.jpsi_mass_window}")
    print(f"  Z window    : {cfg.z_mass_window}")
    print(f"  fourMu mode : {cfg.fourmu_region}")
    print(f"  analysis window: ({ANALYSIS_LOW}, {ANALYSIS_HIGH})")
    print(f"  mask window    : ({MASK_LOW}, {MASK_HIGH})")
    print(f"  step size      : {STEP_SIZE}")
    print(f"  pairing mode: {PAIRING_MODE}")
    print(f"  deduplicate : {DEDUPLICATE}")

    out = output_template()
    total_flat = 0
    total_selected_before_dedup = 0
    total_events_seen = 0
    chunk_index = 0
    t0 = time()

    tree_path = f"{cfg.input_path}:{TREE_NAME}"
    try:
        with uproot.open(tree_path) as tree_for_count:
            total_entries = int(tree_for_count.num_entries)
    except Exception:
        total_entries = None

    for arrays in uproot.iterate(tree_path, expressions=required_branches(), step_size=STEP_SIZE, library="ak"):
        chunk_index += 1
        n_events_chunk = len(arrays["fourMu_mass"])
        total_events_seen += n_events_chunk

        flat = flatten_input_chunk(arrays)
        n_flat_chunk = len(flat["fourMu_mass"])
        total_flat += n_flat_chunk

        if n_flat_chunk == 0:
            if PRINT_PROGRESS and (chunk_index % PROGRESS_EVERY_CHUNKS == 0):
                denom = f"/{total_entries}" if total_entries is not None else ""
                print(f"    chunk {chunk_index:5d}: events {total_events_seen}{denom}, flat +0, selected +0, kept so far {len(out['event'])}", flush=True)
            continue

        chosen = choose_pairing(flat, cfg)
        mask = build_selection_mask(flat, chosen, cfg)
        indices = np.nonzero(mask)[0]

        total_selected_before_dedup += len(indices)

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
                f"selected total {total_selected_before_dedup}, kept so far {len(out['event'])}, "
                f"rate {rate:.1f} events/s",
                flush=True,
            )

    n_before = len(out["event"])
    if DEDUPLICATE:
        out = deduplicate_one_chunk(out)
    n_after = len(out["event"])

    output_path = Path(cfg.output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(output_path) as fout:
        fout[TREE_NAME] = convert_for_uproot(out)

    elapsed = time() - t0
    print(f"  elapsed seconds              : {elapsed:.1f}")
    print(f"  input events processed       : {total_events_seen}")
    print(f"  flattened candidates         : {total_flat}")
    print(f"  selected before deduplication: {total_selected_before_dedup}")
    print(f"  written before deduplication : {n_before}")
    print(f"  written after deduplication  : {n_after}")


def build_signal_config() -> Optional[SampleConfig]:
    outdir = Path(OUTDIR)
    if not MAKE_SIGNAL:
        return None
    c = CUTS["signal"]
    return SampleConfig(
        name="signal",
        input_path=SIGNAL_PRESELECTION,
        output_path=str(outdir / OUTPUTS["signal"]),
        jpsi_mass_window=c["jpsi_mass"],
        z_mass_window=c["z_mass"],
        fourmu_region=c["fourmu_region"],
        label=c["label"],
    )


def build_data_configs() -> List[SampleConfig]:
    outdir = Path(OUTDIR)
    configs: List[SampleConfig] = []

    if MAKE_BACKGROUND:
        c = CUTS["background"]
        configs.append(SampleConfig(
            name="background",
            input_path=DATA_PRESELECTION,
            output_path=str(outdir / OUTPUTS["background"]),
            jpsi_mass_window=c["jpsi_mass"],
            z_mass_window=c["z_mass"],
            fourmu_region=c["fourmu_region"],
            label=c["label"],
        ))

    if MAKE_FINAL_BLINDED:
        c = CUTS["final_blinded"]
        configs.append(SampleConfig(
            name="final_blinded",
            input_path=DATA_PRESELECTION,
            output_path=str(outdir / OUTPUTS["final_blinded"]),
            jpsi_mass_window=c["jpsi_mass"],
            z_mass_window=c["z_mass"],
            fourmu_region=c["fourmu_region"],
            label=c["label"],
        ))

    if MAKE_FINAL_UNBLINDED:
        c = CUTS["final_unblinded"]
        configs.append(SampleConfig(
            name="final_unblinded",
            input_path=DATA_PRESELECTION,
            output_path=str(outdir / OUTPUTS["final_unblinded"]),
            jpsi_mass_window=c["jpsi_mass"],
            z_mass_window=c["z_mass"],
            fourmu_region=c["fourmu_region"],
            label=c["label"],
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


def process_data_once(configs: List[SampleConfig]) -> None:
    """Read the data preselection file once and fill all requested data outputs.

    This avoids reading the same large Run-2 data file separately for background,
    final_blinded, and final_unblinded samples.
    """
    if not configs:
        return

    validate_data_configs(configs)

    input_path = Path(DATA_PRESELECTION)
    if not input_path.exists() and not str(input_path).startswith(("root://", "file:")):
        raise FileNotFoundError(f"Input file not found: {DATA_PRESELECTION}")

    print("\nProcessing data once for:")
    for cfg in configs:
        print(f"  {cfg.name:15s} -> {cfg.output_path}")
    print(f"  input : {DATA_PRESELECTION}")
    print(f"  J/psi window: {configs[0].jpsi_mass_window}")
    print(f"  Z window    : {configs[0].z_mass_window}")
    print(f"  analysis window: ({ANALYSIS_LOW}, {ANALYSIS_HIGH})")
    print(f"  mask window    : ({MASK_LOW}, {MASK_HIGH})")
    print(f"  background region: {CUTS['background']['fourmu_region']}")
    print(f"  step size      : {STEP_SIZE}")
    print(f"  pairing mode   : {PAIRING_MODE}")
    print(f"  deduplicate    : {DEDUPLICATE}")

    outputs = {cfg.name: output_template() for cfg in configs}
    total_flat = 0
    totals_selected = {cfg.name: 0 for cfg in configs}
    total_events_seen = 0
    chunk_index = 0
    t0 = time()

    tree_path = f"{DATA_PRESELECTION}:{TREE_NAME}"
    try:
        with uproot.open(tree_path) as tree_for_count:
            total_entries = int(tree_for_count.num_entries)
    except Exception:
        total_entries = None

    # Pairing choice is common across the data outputs because J/psi/Z windows are common.
    pair_cfg = configs[0]

    for arrays in uproot.iterate(tree_path, expressions=required_branches(), step_size=STEP_SIZE, library="ak"):
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
                f"{name} {len(outputs[name]['event'])}"
                for name in outputs
            )
            print(
                f"    chunk {chunk_index:5d}: events {total_events_seen}{denom}, "
                f"flat +{n_flat_chunk}, {per_chunk}, kept [{kept_so_far}], "
                f"rate {rate:.1f} events/s",
                flush=True,
            )

    for cfg in configs:
        out = outputs[cfg.name]
        before = len(out["event"])
        if DEDUPLICATE:
            out = deduplicate_one_chunk(out)
        after = len(out["event"])

        Path(cfg.output_path).parent.mkdir(parents=True, exist_ok=True)
        with uproot.recreate(cfg.output_path) as fout:
            fout[TREE_NAME] = convert_for_uproot(out)

        print(f"\nFinished {cfg.name}")
        print(f"  output: {cfg.output_path}")
        print(f"  selected before deduplication: {totals_selected[cfg.name]}")
        print(f"  written before deduplication : {before}")
        print(f"  written after deduplication  : {after}")

    elapsed = max(time() - t0, 1e-9)
    print("\nFinished one-pass data processing")
    print(f"  flattened candidates: {total_flat}")
    print(f"  runtime             : {elapsed:.1f} s")
    print(f"  rate                : {total_events_seen / elapsed:.1f} events/s")


def main() -> None:
    signal_cfg = build_signal_config()
    data_cfgs = build_data_configs()

    if signal_cfg is None and not data_cfgs:
        raise RuntimeError("No outputs requested. Enable at least one MAKE_* switch.")

    if signal_cfg is not None:
        process_sample(signal_cfg)

    # The data file is read once for background/final outputs.
    process_data_once(data_cfgs)


if __name__ == "__main__":
    main()
