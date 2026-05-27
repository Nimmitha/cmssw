#!/usr/bin/env python3
"""
Memory-efficient CMSSW-style skimmer for zmmjmm preselection ntuples.

Why this version
----------------
This script avoids the two common memory problems:

  1. It does NOT use uproot to write the skimmed tree, so it does not create
     automatic n<branch> counter branches.

  2. It does NOT collect selected candidates in Python lists. It streams
     event-by-event from the input file to the output file.

It writes a ROOT TTree with:
  - scalar nB
  - std::vector branches with the same branch names as the input
  - only the passing candidates inside each event
  - events with zero passing candidates dropped

Skim logic
----------
No trigger requirement. passMuonTrigger is kept as a branch only.

Keep candidate i if:

  fourMu_vtxProb[i] > 0.001

and either pair group A or B passes:

  group A: pair12 + pair34
  group B: pair14 + pair23

For a group to pass:
  both dimuon vtxProb > 0.001
  one dimuon mass in 1--5 GeV
  the other dimuon mass in 70--110 GeV

Recommended run:
  python scripts/skim_zmmjmm_preselection_root_stream.py

If running on LPC/condor and the job gets killed, reduce PRINT_EVERY only affects
printing, not memory. The script itself is already streaming.
"""

import os
import time
from array import array

import ROOT
ROOT.gInterpreter.GenerateDictionary("std::vector<unsigned int>", "vector")
ROOT.gInterpreter.GenerateDictionary("std::vector<unsigned long long>", "vector")
ROOT.gInterpreter.GenerateDictionary("std::vector<bool>", "vector")

ROOT.gROOT.SetBatch(True)

# =============================================================================
# User config
# =============================================================================

INPUT_FILE = "/uscms/home/wkarunar/nobackup/analysis/run2/reboot/CMSSW_10_6_30/src/ZmmJmmAnalyzer/preselection/run2_data/v2/TTree_13TeV_mmmm_UL_2018.root"
OUTPUT_FILE = "preselection/skimmed_with_ROOT_v2/skimmed_TTree_13TeV_mmmm_UL_2018.root"

TREE_NAME = "ntuple"

LOW_MASS_WINDOW = (1.0, 5.0)
Z_MASS_WINDOW = (70.0, 110.0)

PAIR_VTX_MIN = 0.001
FOURMU_VTX_MIN = 0.001

PRINT_EVERY = 100000

# Optional compression. Higher compression saves space but may be slower.
# ROOT default is usually fine. 404 means ZLIB level 4.
COMPRESSION_SETTINGS = 404

# If True, prints branch type summary at start.
PRINT_BRANCH_SUMMARY = False


# =============================================================================
# Helpers
# =============================================================================

def in_window(x, window):
    return window[0] < float(x) < window[1]


def group_passes(m_a, m_b, vtx_a, vtx_b):
    if float(vtx_a) <= PAIR_VTX_MIN or float(vtx_b) <= PAIR_VTX_MIN:
        return False
    return (
        (in_window(m_a, LOW_MASS_WINDOW) and in_window(m_b, Z_MASS_WINDOW))
        or
        (in_window(m_b, LOW_MASS_WINDOW) and in_window(m_a, Z_MASS_WINDOW))
    )


def get_input_tree(fin):
    tree = fin.Get(TREE_NAME)
    if tree:
        return tree

    tree = fin.Get(f"rootuple/{TREE_NAME}")
    if tree:
        return tree

    raise RuntimeError(f"Could not find '{TREE_NAME}' or 'rootuple/{TREE_NAME}' in input file")


def required_branches_exist(tree):
    required = [
        "fourMu_mass",
        "fourMu_vtxProb",
        "pair12_mass", "pair12_vtxProb",
        "pair34_mass", "pair34_vtxProb",
        "pair14_mass", "pair14_vtxProb",
        "pair23_mass", "pair23_vtxProb",
    ]
    missing = [b for b in required if not tree.GetBranch(b)]
    if missing:
        raise RuntimeError("Missing required branches: " + ", ".join(missing))


def n_candidates(tree):
    # Prefer vector size. It is safer than trusting nB in files produced by
    # temporary skimmers.
    return int(tree.fourMu_mass.size())


def candidate_passes(tree, i):
    if float(tree.fourMu_vtxProb[i]) <= FOURMU_VTX_MIN:
        return False

    group_a = group_passes(
        tree.pair12_mass[i],
        tree.pair34_mass[i],
        tree.pair12_vtxProb[i],
        tree.pair34_vtxProb[i],
    )

    group_b = group_passes(
        tree.pair14_mass[i],
        tree.pair23_mass[i],
        tree.pair14_vtxProb[i],
        tree.pair23_vtxProb[i],
    )

    return group_a or group_b


def is_vector_branch(branch):
    cls = branch.GetClassName().replace("std::", "")
    return cls.startswith("vector<")


def vector_element_type(branch):
    cls = branch.GetClassName().replace("std::", "")
    elem = cls[len("vector<"):]
    if elem.endswith(">"):
        elem = elem[:-1]
    return elem.strip()


def normalize_vector_type(elem_type):
    aliases = {
        "Bool_t": "bool",
        "Int_t": "int",
        "UInt_t": "unsigned int",
        "Float_t": "float",
        "Double_t": "double",
        "ULong64_t": "unsigned long long",
        "Long64_t": "long long",
    }
    return aliases.get(elem_type, elem_type)


def make_vector(elem_type):
    return ROOT.std.vector(normalize_vector_type(elem_type))()


def scalar_leaf_code(leaf):
    typename = leaf.GetTypeName()
    mapping = {
        "UInt_t": ("I", "i"),   # Python array code, ROOT leaflist code
        "unsigned int": ("I", "i"),
        "Int_t": ("i", "I"),
        "int": ("i", "I"),
        "Float_t": ("f", "F"),
        "float": ("f", "F"),
        "Double_t": ("d", "D"),
        "double": ("d", "D"),
        "ULong64_t": ("Q", "l"),
        "unsigned long long": ("Q", "l"),
        "Long64_t": ("q", "L"),
        "long long": ("q", "L"),
        "Bool_t": ("b", "O"),
        "bool": ("b", "O"),
    }
    return mapping.get(typename)


def build_output_tree(input_tree):
    """
    Create output tree and output buffers.

    For vector branches, output buffers are std::vector<T> objects.
    For scalar branches, output buffers are array.array objects.

    nB is always written as UInt_t and set to the number of kept candidates.
    """
    out_tree = ROOT.TTree(TREE_NAME, TREE_NAME)

    vector_buffers = {}
    scalar_buffers = {}

    # Always create nB once, even if input has nB.
    nB = array("I", [0])
    out_tree.Branch("nB", nB, "nB/i")
    scalar_buffers["nB"] = (nB, None)

    for branch in input_tree.GetListOfBranches():
        name = branch.GetName()

        if name == "nB":
            continue

        if is_vector_branch(branch):
            elem = vector_element_type(branch)
            vec = make_vector(elem)
            vector_buffers[name] = vec
            out_tree.Branch(name, vec)
            continue

        # Rare scalar branches other than nB. Copy if simple.
        leaves = branch.GetListOfLeaves()
        if leaves.GetEntries() != 1:
            continue

        leaf = leaves.At(0)
        codes = scalar_leaf_code(leaf)
        if codes is None:
            continue

        py_code, root_code = codes
        if py_code in ("b",):
            holder = array(py_code, [0])
        elif py_code in ("i", "I", "q", "Q"):
            holder = array(py_code, [0])
        else:
            holder = array(py_code, [0.0])

        scalar_buffers[name] = (holder, root_code)
        out_tree.Branch(name, holder, f"{name}/{root_code}")

    if PRINT_BRANCH_SUMMARY:
        print(f"Vector branches: {len(vector_buffers)}")
        print(f"Scalar branches: {len(scalar_buffers)}")
        for k in list(vector_buffers)[:20]:
            print("  vector", k)

    return out_tree, vector_buffers, scalar_buffers


def copy_selected_vectors(input_tree, vector_buffers, keep_indices):
    """
    Clear output vectors and copy only passing candidate indices.
    """
    for name, out_vec in vector_buffers.items():
        out_vec.clear()

        in_vec = getattr(input_tree, name)
        n = int(in_vec.size())

        # Candidate-level vector branches should have same length as fourMu_mass.
        # If a vector has a different length for some reason, copy only if index exists.
        for idx in keep_indices:
            if idx < n:
                out_vec.push_back(in_vec[idx])


def copy_scalars(input_tree, scalar_buffers, n_kept):
    """
    Fill scalar buffers. nB is overwritten with number of kept candidates.
    """
    for name, (holder, root_code) in scalar_buffers.items():
        if name == "nB":
            holder[0] = int(n_kept)
            continue

        if not hasattr(input_tree, name):
            continue

        value = getattr(input_tree, name)
        try:
            holder[0] = value
        except TypeError:
            # Bool_t branch support
            holder[0] = int(bool(value))


def main():
    if not os.path.exists(INPUT_FILE):
        raise FileNotFoundError(INPUT_FILE)

    outdir = os.path.dirname(OUTPUT_FILE)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    fin = ROOT.TFile.Open(INPUT_FILE, "READ")
    if not fin or fin.IsZombie():
        raise RuntimeError(f"Could not open input file: {INPUT_FILE}")

    tin = get_input_tree(fin)
    required_branches_exist(tin)

    fout = ROOT.TFile.Open(OUTPUT_FILE, "RECREATE")
    if not fout or fout.IsZombie():
        raise RuntimeError(f"Could not create output file: {OUTPUT_FILE}")
    fout.SetCompressionSettings(COMPRESSION_SETTINGS)

    fout.cd()
    tout, vector_buffers, scalar_buffers = build_output_tree(tin)

    n_entries = int(tin.GetEntries())

    events_written = 0
    events_dropped = 0
    candidates_in = 0
    candidates_out = 0

    t0 = time.time()

    print(f"Input : {INPUT_FILE}")
    print(f"Output: {OUTPUT_FILE}")
    print(f"Tree  : {TREE_NAME}")
    print(f"Entries: {n_entries}")
    print(f"Vector branches copied: {len(vector_buffers)}")
    print("Starting streaming skim...")

    for entry in range(n_entries):
        tin.GetEntry(entry)

        ncand = n_candidates(tin)
        candidates_in += ncand

        keep = []
        for i in range(ncand):
            if candidate_passes(tin, i):
                keep.append(i)

        if keep:
            copy_selected_vectors(tin, vector_buffers, keep)
            copy_scalars(tin, scalar_buffers, len(keep))
            tout.Fill()

            events_written += 1
            candidates_out += len(keep)
        else:
            events_dropped += 1

        if PRINT_EVERY and (entry + 1) % PRINT_EVERY == 0:
            dt = max(time.time() - t0, 1e-9)
            frac = 100.0 * (entry + 1) / n_entries
            print(
                f"  processed {entry + 1}/{n_entries} ({frac:5.1f}%), "
                f"events kept {events_written}, "
                f"cand in {candidates_in}, cand kept {candidates_out}, "
                f"rate {(entry + 1) / dt:.1f} ev/s",
                flush=True,
            )

    fout.cd()
    tout.Write("", ROOT.TObject.kOverwrite)
    fout.Close()
    fin.Close()

    dt = max(time.time() - t0, 1e-9)

    print("\nDone.")
    print(f"Events written : {events_written}")
    print(f"Events dropped : {events_dropped}")
    print(f"Candidates in  : {candidates_in}")
    print(f"Candidates out : {candidates_out}")
    print(f"Runtime        : {dt:.1f} s")
    print(f"Rate           : {n_entries / dt:.1f} events/s")


if __name__ == "__main__":
    main()
