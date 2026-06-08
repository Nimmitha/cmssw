#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import uproot


def load_tree(path, tree_name):
    with uproot.open(path) as f:
        tree = f[tree_name]
        df = tree.arrays(library="pd")
    return df


def asimov_significance(s, b):
    if s <= 0 or b <= 0:
        return 0.0
    return np.sqrt(2.0 * ((s + b) * np.log(1.0 + s / b) - s))


def clean_numeric_columns(df_mc, df_data, protected_columns):
    common = sorted(set(df_mc.columns).intersection(df_data.columns))

    usable = []
    for col in common:
        if col in protected_columns:
            continue

        if not np.issubdtype(df_mc[col].dtype, np.number):
            continue

        if not np.issubdtype(df_data[col].dtype, np.number):
            continue

        mc_values = df_mc[col].to_numpy()
        data_values = df_data[col].to_numpy()

        if np.all(~np.isfinite(mc_values)) or np.all(~np.isfinite(data_values)):
            continue

        # Skip constant or almost-constant variables
        finite_mc = mc_values[np.isfinite(mc_values)]
        finite_data = data_values[np.isfinite(data_values)]

        if len(finite_mc) < 10 or len(finite_data) < 10:
            continue

        if np.nanstd(finite_mc) == 0 and np.nanstd(finite_data) == 0:
            continue

        usable.append(col)

    return usable


def evaluate_cut(sig, bkg, variable, threshold, direction):
    if direction == ">":
        sig_pass = sig[variable] > threshold
        bkg_pass = bkg[variable] > threshold
    elif direction == "<":
        sig_pass = sig[variable] < threshold
        bkg_pass = bkg[variable] < threshold
    else:
        raise ValueError("direction must be '>' or '<'")

    s0 = len(sig)
    b0 = len(bkg)

    s = int(sig_pass.sum())
    b = int(bkg_pass.sum())

    sig_eff = s / s0 if s0 > 0 else 0.0
    bkg_eff = b / b0 if b0 > 0 else 0.0
    bkg_rej = 1.0 - bkg_eff

    s_over_sqrt_b = s / np.sqrt(b) if b > 0 else np.inf
    za = asimov_significance(s, b)

    return {
        "variable": variable,
        "direction": direction,
        "threshold": threshold,
        "s_before": s0,
        "b_before": b0,
        "s_after": s,
        "b_after": b,
        "signal_eff": sig_eff,
        "background_eff": bkg_eff,
        "background_rej": bkg_rej,
        "s_over_sqrt_b": s_over_sqrt_b,
        "asimov_Z": za,
    }


def scan_variable(sig, bkg, variable, n_thresholds, min_signal_eff):
    values = pd.concat([sig[variable], bkg[variable]], ignore_index=True)
    values = values.to_numpy()
    values = values[np.isfinite(values)]

    if len(values) < 20:
        return []

    # Avoid extreme unstable thresholds
    qs = np.linspace(0.02, 0.98, n_thresholds)
    thresholds = np.unique(np.quantile(values, qs))

    results = []

    for threshold in thresholds:
        for direction in [">", "<"]:
            result = evaluate_cut(sig, bkg, variable, threshold, direction)

            if result["signal_eff"] >= min_signal_eff:
                results.append(result)

    return results


def apply_cut(df, variable, threshold, direction):
    if direction == ">":
        return df[df[variable] > threshold].copy()
    elif direction == "<":
        return df[df[variable] < threshold].copy()
    else:
        raise ValueError("direction must be '>' or '<'")


def greedy_scan(sig, bkg, variables, n_steps, n_thresholds, min_signal_eff_per_step):
    selected_cuts = []

    current_sig = sig.copy()
    current_bkg = bkg.copy()

    for step in range(n_steps):
        all_results = []

        for variable in variables:
            results = scan_variable(
                current_sig,
                current_bkg,
                variable,
                n_thresholds=n_thresholds,
                min_signal_eff=min_signal_eff_per_step,
            )
            all_results.extend(results)

        if len(all_results) == 0:
            break

        ranked = pd.DataFrame(all_results)
        ranked = ranked.replace([np.inf, -np.inf], np.nan).dropna()

        if len(ranked) == 0:
            break

        # Main ranking metric.
        # You can change this to "s_over_sqrt_b".
        ranked = ranked.sort_values("asimov_Z", ascending=False)

        best = ranked.iloc[0].to_dict()
        selected_cuts.append(best)

        current_sig = apply_cut(
            current_sig,
            best["variable"],
            best["threshold"],
            best["direction"],
        )

        current_bkg = apply_cut(
            current_bkg,
            best["variable"],
            best["threshold"],
            best["direction"],
        )

        print()
        print(f"Step {step + 1}")
        print(
            f"  Cut: {best['variable']} {best['direction']} {best['threshold']:.6g}"
        )
        print(
            f"  Signal: {len(current_sig)} / {len(sig)} "
            f"= {len(current_sig) / len(sig):.4f}"
        )
        print(
            f"  Bkg:    {len(current_bkg)} / {len(bkg)} "
            f"= {len(current_bkg) / len(bkg):.4f}"
        )

    return pd.DataFrame(selected_cuts), current_sig, current_bkg


def main():
    parser = argparse.ArgumentParser(
        description="Scan selection variables to find signal-efficient cuts."
    )

    parser.add_argument("--mc", required=True, help="Signal MC ROOT file")
    parser.add_argument("--data", required=True, help="Data ROOT file")
    parser.add_argument("--tree", default="Events", help="Tree name")

    parser.add_argument("--mass-branch", default="fourL_mass")
    parser.add_argument("--blind-low", type=float, default=120.0)
    parser.add_argument("--blind-high", type=float, default=130.0)

    parser.add_argument("--min-signal-eff", type=float, default=0.85)
    parser.add_argument("--n-thresholds", type=int, default=80)
    parser.add_argument("--greedy-steps", type=int, default=0)

    parser.add_argument("--output", default="cut_scan_results.csv")
    parser.add_argument("--greedy-output", default="greedy_cuts.csv")

    parser.add_argument(
    "--variables",
    nargs="+",
    default=None,
    help="Specific variables to scan. Example: --variables mu1_pfRelIso03 mu2_pfRelIso03 Jpsi_relIso03",
    )

    args = parser.parse_args()

    mc = load_tree(args.mc, args.tree)
    data = load_tree(args.data, args.tree)

    if args.mass_branch not in mc.columns:
        raise RuntimeError(f"Missing mass branch in MC: {args.mass_branch}")

    if args.mass_branch not in data.columns:
        raise RuntimeError(f"Missing mass branch in data: {args.mass_branch}")

    # Signal sample: use all MC candidates by default.
    # You can tighten this if the MC contains multiple components.
    sig = mc.copy()

    # Background proxy: data sidebands only.
    bkg = data[
        (data[args.mass_branch] < args.blind_low)
        | (data[args.mass_branch] > args.blind_high)
    ].copy()

    print(f"Loaded MC candidates:   {len(mc)}")
    print(f"Loaded data candidates: {len(data)}")
    print(f"Signal candidates used: {len(sig)}")
    print(f"Data sideband used:     {len(bkg)}")

    protected_columns = {
        args.mass_branch,
        "run",
        "luminosityBlock",
        "event",
        "event_number",
        "candidate_index",
        "weight",
        "genWeight",
    }

    if args.variables is not None:
        missing = [v for v in args.variables if v not in sig.columns or v not in bkg.columns]
        if missing:
            raise RuntimeError(f"Requested variables are missing from MC or data: {missing}")

        variables = args.variables
    else:
        variables = clean_numeric_columns(sig, bkg, protected_columns)

    print(f"Scanning {len(variables)} variables")

    all_results = []

    for variable in variables:
        results = scan_variable(
            sig,
            bkg,
            variable,
            n_thresholds=args.n_thresholds,
            min_signal_eff=args.min_signal_eff,
        )
        all_results.extend(results)

    results_df = pd.DataFrame(all_results)

    if len(results_df) == 0:
        print("No cuts passed the requested signal-efficiency threshold.")
        return

    results_df = results_df.replace([np.inf, -np.inf], np.nan)
    results_df = results_df.sort_values(
        ["asimov_Z", "background_rej", "signal_eff"],
        ascending=[False, False, False],
    )

    results_df.to_csv(args.output, index=False)

    print()
    print("Top single-variable cuts:")
    print(
        results_df[
            [
                "variable",
                "direction",
                "threshold",
                "signal_eff",
                "background_eff",
                "background_rej",
                "s_after",
                "b_after",
                "asimov_Z",
                "s_over_sqrt_b",
            ]
        ]
        .head(20)
        .to_string(index=False)
    )

    print()
    print(f"Saved full scan to: {args.output}")

    if args.greedy_steps > 0:
        greedy_df, final_sig, final_bkg = greedy_scan(
            sig,
            bkg,
            variables,
            n_steps=args.greedy_steps,
            n_thresholds=args.n_thresholds,
            min_signal_eff_per_step=args.min_signal_eff,
        )

        greedy_df.to_csv(args.greedy_output, index=False)

        print()
        print(f"Saved greedy cuts to: {args.greedy_output}")
        print()
        print("Final greedy efficiency:")
        print(f"  Signal efficiency: {len(final_sig) / len(sig):.4f}")
        print(f"  Bkg efficiency:    {len(final_bkg) / len(bkg):.4f}")


if __name__ == "__main__":
    main()