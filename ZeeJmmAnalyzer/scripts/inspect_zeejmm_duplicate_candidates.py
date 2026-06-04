#!/usr/bin/env python3
"""
Inspect duplicate selected ZeeJmm candidates in the output of make_zeejmm_candidates_*.py.

Purpose
-------
This script is for checking cases where the selected flat candidate tree contains
more candidates than unique data events. It groups selected candidates by
(run, lumi, event), lists events with multiplicity > 1, and writes simple plots
for the basic kinematic/vertex variables that allowed the candidates to pass.

Important MC caveat
-------------------
For private MC made by merging files that reuse run/lumi/event numbers, the
(run,lumi,event) key is not globally unique. In that case, this script may report
fake duplicates across different generated files. Use this mainly for data, or
for MC only when event ids are known to be unique.

Nominal use
-----------
  python3 inspect_zeejmm_duplicate_candidates.py

Outputs
-------
  duplicate_event_summary.csv     one row per duplicate event key
  duplicate_candidate_rows.csv    one row per candidate in duplicate events
  duplicate_candidates.root       histograms + optional duplicate-only tree
  plots/*.png                     quick-look plots
"""

from __future__ import annotations

import csv
from array import array
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import ROOT

ROOT.gROOT.SetBatch(True)

# =============================================================================
# User configuration
# =============================================================================

INPUT_FILE = "selection/signal_candidates.root"
TREE_NAME = "ntuple"
OUTDIR = "duplicate_check"

# Use True for data. For merged private MC with repeated run/lumi/event, these
# duplicates are not physically meaningful unless you made the event ids unique.
GROUP_BY_RUN_LUMI_EVENT = True

# Print only the first N duplicate groups in the terminal. CSV files contain all.
MAX_PRINT_GROUPS = 20
MAX_PRINT_CANDIDATES_PER_GROUP = 10

# Write a duplicate-only ROOT tree in addition to histograms.
WRITE_DUPLICATE_ONLY_TREE = True

# Variables to print/write for understanding why multiple candidates passed.
# Keep this focused: masses, pT, eta, and vertex probabilities.
BASIC_FLOAT_VARS = [
    "fourL_mass", "fourL_pt", "fourL_vtxProb",
    "Z_mass", "Z_pt", "Z_vtxProb",
    "Jpsi_mass", "Jpsi_pt", "Jpsi_vtxProb",
    "e1_pt", "e2_pt", "e1_eta", "e2_eta",
    "mu1_pt", "mu2_pt", "mu1_eta", "mu2_eta",
]

BASIC_INT_VARS = [
    "run", "lumi", "event", "label",
    "passEleTrigger", "passEleTriggerMatch",
    "e1_passWP90", "e2_passWP90",
    "mu1_soft", "mu2_soft",
]

# Variables to histogram for all candidates in duplicate events.
HIST_SPECS = {
    "candidate_multiplicity": (20, 0.5, 20.5),
    "fourL_mass": (50, 112.0, 162.0),
    "fourL_vtxProb": (50, 0.0, 1.0),
    "Z_mass": (50, 70.0, 110.0),
    "Z_vtxProb": (50, 0.0, 1.0),
    "Jpsi_mass": (60, 2.8, 4.0),
    "Jpsi_vtxProb": (50, 0.0, 1.0),
    "e1_pt": (60, 0.0, 120.0),
    "e2_pt": (60, 0.0, 80.0),
    "mu1_pt": (60, 0.0, 60.0),
    "mu2_pt": (60, 0.0, 60.0),
}

# =============================================================================
# Helpers
# =============================================================================


def require_branch(tree: ROOT.TTree, branch: str) -> None:
    if not tree.GetBranch(branch):
        raise RuntimeError(f"Missing required branch: {branch}")


def has_branch(tree: ROOT.TTree, branch: str) -> bool:
    return bool(tree.GetBranch(branch))


def get_value(tree: ROOT.TTree, branch: str):
    return getattr(tree, branch)


def event_key(tree: ROOT.TTree) -> Tuple[int, int, int]:
    return (int(get_value(tree, "run")), int(get_value(tree, "lumi")), int(get_value(tree, "event")))


def candidate_record(tree: ROOT.TTree, entry_index: int, available_floats: List[str], available_ints: List[str]) -> Dict[str, float | int]:
    rec: Dict[str, float | int] = {"entry_index": int(entry_index)}
    for var in available_ints:
        rec[var] = int(get_value(tree, var))
    for var in available_floats:
        rec[var] = float(get_value(tree, var))
    return rec


def format_candidate(rec: Dict[str, float | int]) -> str:
    return (
        f"entry={rec['entry_index']} "
        f"fourL_mass={rec.get('fourL_mass', -999):.3f} fourL_vtx={rec.get('fourL_vtxProb', -999):.4f} "
        f"Z_mass={rec.get('Z_mass', -999):.3f} Z_vtx={rec.get('Z_vtxProb', -999):.4f} "
        f"J_mass={rec.get('Jpsi_mass', -999):.3f} J_vtx={rec.get('Jpsi_vtxProb', -999):.4f} "
        f"ept=({rec.get('e1_pt', -999):.2f},{rec.get('e2_pt', -999):.2f}) "
        f"mupt=({rec.get('mu1_pt', -999):.2f},{rec.get('mu2_pt', -999):.2f}) "
        f"WP90=({rec.get('e1_passWP90', -999)},{rec.get('e2_passWP90', -999)})"
    )


def write_csvs(outdir: Path, duplicates: Dict[Tuple[int, int, int], List[Dict[str, float | int]]], available_floats: List[str], available_ints: List[str]) -> None:
    summary_path = outdir / "duplicate_event_summary.csv"
    rows_path = outdir / "duplicate_candidate_rows.csv"

    with summary_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "group", "run", "lumi", "event", "multiplicity",
            "best_entry_by_fourL_vtxProb", "best_fourL_vtxProb",
            "fourL_mass_min", "fourL_mass_max",
            "Jpsi_mass_min", "Jpsi_mass_max",
            "Z_mass_min", "Z_mass_max",
        ])
        for ig, (key, cands) in enumerate(sorted(duplicates.items()), start=1):
            run, lumi, event = key
            best = max(cands, key=lambda r: float(r.get("fourL_vtxProb", -999.0)))
            def vals(name: str) -> List[float]:
                return [float(c[name]) for c in cands if name in c]
            writer.writerow([
                ig, run, lumi, event, len(cands),
                best["entry_index"], best.get("fourL_vtxProb", -999.0),
                min(vals("fourL_mass")) if vals("fourL_mass") else -999.0,
                max(vals("fourL_mass")) if vals("fourL_mass") else -999.0,
                min(vals("Jpsi_mass")) if vals("Jpsi_mass") else -999.0,
                max(vals("Jpsi_mass")) if vals("Jpsi_mass") else -999.0,
                min(vals("Z_mass")) if vals("Z_mass") else -999.0,
                max(vals("Z_mass")) if vals("Z_mass") else -999.0,
            ])

    row_fields = ["group", "candidate_in_group", "multiplicity", "entry_index"] + available_ints + available_floats
    # Avoid duplicated entry_index if someone adds it to variable lists later.
    seen = set()
    row_fields = [x for x in row_fields if not (x in seen or seen.add(x))]

    with rows_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=row_fields, extrasaction="ignore")
        writer.writeheader()
        for ig, (_key, cands) in enumerate(sorted(duplicates.items()), start=1):
            for ic, rec in enumerate(cands, start=1):
                row = dict(rec)
                row["group"] = ig
                row["candidate_in_group"] = ic
                row["multiplicity"] = len(cands)
                writer.writerow(row)

    print(f"  wrote duplicate event summary : {summary_path}")
    print(f"  wrote duplicate candidate rows: {rows_path}")


def make_histograms(outdir: Path, duplicates: Dict[Tuple[int, int, int], List[Dict[str, float | int]]], available_floats: List[str]) -> None:
    root_path = outdir / "duplicate_candidates.root"
    plot_dir = outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    fout = ROOT.TFile.Open(str(root_path), "RECREATE")

    hists: Dict[str, ROOT.TH1F] = {}
    for name, (nbins, lo, hi) in HIST_SPECS.items():
        h = ROOT.TH1F(f"h_{name}", f"{name};{name};Candidates", nbins, lo, hi)
        h.SetLineWidth(2)
        hists[name] = h

    for _key, cands in duplicates.items():
        hists["candidate_multiplicity"].Fill(len(cands))
        for rec in cands:
            for name in HIST_SPECS:
                if name == "candidate_multiplicity":
                    continue
                if name in rec:
                    hists[name].Fill(float(rec[name]))

    canvas = ROOT.TCanvas("c", "c", 900, 700)
    for name, hist in hists.items():
        fout.cd()
        hist.Write()
        canvas.Clear()
        hist.Draw("hist")
        canvas.SaveAs(str(plot_dir / f"{name}.png"))

    fout.Close()
    print(f"  wrote histogram ROOT file     : {root_path}")
    print(f"  wrote plots                  : {plot_dir}/*.png")


def write_duplicate_only_tree(input_file: str, tree_name: str, outdir: Path, duplicate_entries: Iterable[int]) -> None:
    if not WRITE_DUPLICATE_ONLY_TREE:
        return

    entry_set = set(int(x) for x in duplicate_entries)
    out_path = outdir / "duplicate_only_tree.root"

    fin = ROOT.TFile.Open(input_file)
    tree = fin.Get(tree_name)
    fout = ROOT.TFile.Open(str(out_path), "RECREATE")
    clone = tree.CloneTree(0)

    for i in range(tree.GetEntries()):
        if i in entry_set:
            tree.GetEntry(i)
            clone.Fill()

    fout.cd()
    clone.Write()
    fout.Close()
    fin.Close()
    print(f"  wrote duplicate-only tree     : {out_path}")


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    outdir = Path(OUTDIR)
    outdir.mkdir(parents=True, exist_ok=True)

    print("ZeeJmm duplicate-candidate inspection")
    print(f"  input : {INPUT_FILE}")
    print(f"  tree  : {TREE_NAME}")
    print(f"  outdir: {OUTDIR}")
    print("  grouping key: run/lumi/event")

    fin = ROOT.TFile.Open(INPUT_FILE)
    if not fin or fin.IsZombie():
        raise RuntimeError(f"Could not open input file: {INPUT_FILE}")

    tree = fin.Get(TREE_NAME)
    if not tree:
        raise RuntimeError(f"Could not find tree '{TREE_NAME}' in {INPUT_FILE}")

    for branch in ["run", "lumi", "event"]:
        require_branch(tree, branch)

    available_floats = [x for x in BASIC_FLOAT_VARS if has_branch(tree, x)]
    available_ints = [x for x in BASIC_INT_VARS if has_branch(tree, x)]

    missing_floats = [x for x in BASIC_FLOAT_VARS if x not in available_floats]
    missing_ints = [x for x in BASIC_INT_VARS if x not in available_ints]
    if missing_floats or missing_ints:
        print("  warning: some requested branches are missing and will be skipped")
        if missing_floats:
            print("    missing floats:", ", ".join(missing_floats))
        if missing_ints:
            print("    missing ints  :", ", ".join(missing_ints))

    groups: Dict[Tuple[int, int, int], List[Dict[str, float | int]]] = defaultdict(list)
    n_entries = int(tree.GetEntries())

    for i in range(n_entries):
        tree.GetEntry(i)
        key = event_key(tree)
        rec = candidate_record(tree, i, available_floats, available_ints)
        groups[key].append(rec)

    duplicates = {key: cands for key, cands in groups.items() if len(cands) > 1}
    duplicate_candidate_count = sum(len(cands) for cands in duplicates.values())
    extra_candidate_count = sum(len(cands) - 1 for cands in duplicates.values())

    print("\nSummary")
    print(f"  selected candidate rows       : {n_entries}")
    print(f"  unique run/lumi/event keys    : {len(groups)}")
    print(f"  duplicate event keys          : {len(duplicates)}")
    print(f"  candidates in duplicate keys  : {duplicate_candidate_count}")
    print(f"  extra candidates beyond 1/key : {extra_candidate_count}")

    if not duplicates:
        print("\nNo duplicate selected events found.")
        fin.Close()
        return

    print(f"\nFirst {min(MAX_PRINT_GROUPS, len(duplicates))} duplicate groups")
    for ig, (key, cands) in enumerate(sorted(duplicates.items()), start=1):
        if ig > MAX_PRINT_GROUPS:
            break
        run, lumi, event = key
        best = max(cands, key=lambda r: float(r.get("fourL_vtxProb", -999.0)))
        print(f"\n  group {ig}: run={run} lumi={lumi} event={event} multiplicity={len(cands)}")
        print(f"    best by highest fourL_vtxProb: entry={best['entry_index']} fourL_vtxProb={best.get('fourL_vtxProb', -999):.4f}")
        for rec in cands[:MAX_PRINT_CANDIDATES_PER_GROUP]:
            print("    " + format_candidate(rec))
        if len(cands) > MAX_PRINT_CANDIDATES_PER_GROUP:
            print(f"    ... {len(cands) - MAX_PRINT_CANDIDATES_PER_GROUP} more candidates in this group")

    write_csvs(outdir, duplicates, available_floats, available_ints)
    make_histograms(outdir, duplicates, available_floats)

    duplicate_entries = [int(rec["entry_index"]) for cands in duplicates.values() for rec in cands]
    fin.Close()
    write_duplicate_only_tree(INPUT_FILE, TREE_NAME, outdir, duplicate_entries)

    print("\nDone.")


if __name__ == "__main__":
    main()
