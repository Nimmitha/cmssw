#!/usr/bin/env python3
"""
Full pipeline script to process all eras:
1. Convert ROOT files to HDF5
2. Run ML analysis on HDF5 files
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path

# Era configuration
era_config = {
    "2022": [
        "C1", "D1", "E1", "F1", "G1"
    ],
    "2023": [
        "B1", "C1", "C2", "C3", "C4", "D1", "D2",
    ],
    "2024": [
        "C1", "D1", "E1", "E2", "F1", "G1", "H1", "I1", "I2"
    ],
    "2025": [
        "B1", "C1", "C2", "D1"
    ]
}

# Base path pattern
# base_path_pattern = "/uscms/home/wkarunar/nobackup/datasets/data/run3/parkingDoubleMuonLowMass/{year}/v6"
base_path_pattern = "/home/nimmitha/LPCfiles/run3/jmm_counting/cmssw/ZmmJmmAnalyzer/preselection/parkingDoubleMuonLowMass/{year}/v6"

# Cuts definition
cuts_config = {
    "B_J1_mass": (2.7, 3.5),
    "B_Mu1_pt": (5000, None),
    "B_Mu2_pt": (5000, None)
}


def run_command(cmd, description=""):
    """Run a shell command and handle errors."""
    if description:
        print(f"\n{description}")
    
    print(f"Running: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, text=True)  # stderr still goes to terminal
    
    if result.returncode != 0:
        print(f"Error running command: {' '.join(cmd)}")
        # print(f"STDOUT: {result.stdout}")
        # print(f"STDERR: {result.stderr}")
        sys.exit(1)
    
    return result


def process_era(year, era, base_path, args):
    """Process a single era through both stages."""
    
    print(f"\n{'='*60}")
    print(f"Processing {year} Era {era}")
    print(f"{'='*60}")
    
    # File paths
    root_file = f"{base_path}/PDMLM_mm_{year}{era}_v6.root"
    hdf5_file = f"{args.hdf5_dir}/{year}_{era}.h5"

    if year == "2022" or year == "2023":
        normtag = "PHYSICS"
    elif year == "2024" or year == "2025":
        normtag = "BRIL"
    else:
        raise ValueError(f"Unknown year {year}. Please define the lumi file name.")

    lumi_file = f"{base_path}/lumi_normtag{normtag}_PDMLM_{year}{era}.csv"
    output_csv = f"{args.csv_dir}/event_data_fit_{year}_{era}_{normtag}.csv"

    if args.no_convert:
        print(f"Skipping conversion for {year} {era}")
        if not os.path.exists(hdf5_file):
            print(f"Error: HDF5 file does not exist, cannot skip conversion: {hdf5_file}")
            return False
    else:
        # Check if ROOT file exists
        if not os.path.exists(root_file) and not args.no_convert:
            print(f"Warning: ROOT file not found: {root_file}")
            return False
        
        # Stage 1: ROOT to HDF5 conversion (if needed)
        if not os.path.exists(hdf5_file) or args.force_convert:
            cmd = [
                "python", "root_to_hdf5.py",
                root_file,
                hdf5_file,
                "--bin-width", str(args.bin_width),
                "--chunk-size", args.chunk_size,
                "--compression", args.compression,
                "--verify"  # Add verification after conversion
            ]
            
            run_command(cmd, f"Stage 1: Converting ROOT to HDF5 for {year} {era}")
        else:
            print(f"HDF5 file already exists: {hdf5_file}")
        
        # Check if HDF5 was created successfully
        if not os.path.exists(hdf5_file):
            print(f"Error: HDF5 file was not created: {hdf5_file}")
            return False
    
    # Stage 2: ML Analysis
    if args.no_analyze:
        print(f"Skipping ML analysis for {year} {era}")
    else:    
        if not os.path.exists(output_csv) or args.force_analyze:
            cmd = [
                "python", "hdf5_ml_analysis.py",
                hdf5_file,
                output_csv,
                "--mass-min", str(args.mass_min),
                "--mass-max", str(args.mass_max),
                "--n-cores", str(args.n_cores),
                "--year", year,
                "--era", era
            ]
            
            # Add luminosity file if it exists
            if os.path.exists(lumi_file):
                cmd.extend(["--lumi-file", lumi_file])
            else:
                print(f"Warning: Luminosity file not found: {lumi_file}")
            
            # Add plot saving options
            if args.save_plots:
                cmd.append("--save-plots")
                cmd.extend(["--plot-dir", args.plot_dir])
            
            run_command(cmd, f"Stage 2: ML Analysis for {year} {era}")
        else:
            print(f"Output CSV already exists: {output_csv}")
    
    return True


def main():
    parser = argparse.ArgumentParser(description='Process all eras through the full pipeline')
    parser.add_argument('--years', nargs='+', help='Years to process (default: 2024)')
    parser.add_argument('--eras', nargs='+', 
                       help='Specific eras to process (overrides config)')
    parser.add_argument('--hdf5-dir', default='hdf5_files', 
                       help='Directory for HDF5 files')
    parser.add_argument('--csv-dir', default='csvs/event_data', 
                       help='Directory for output CSV files')
    parser.add_argument('--plot-dir', default='plots/ML_fits', 
                       help='Directory for fit plots')
    parser.add_argument('--bin-width', type=int, default=50, 
                       help='Lumiblock bin width')
    parser.add_argument('--chunk-size', default='200 MB', 
                       help='Chunk size for ROOT reading')
    parser.add_argument('--compression', default='gzip', 
                       choices=['gzip', 'lzf'],
                       help='HDF5 compression algorithm')
    parser.add_argument('--mass-min', type=float, default=2.7, 
                       help='Minimum mass for fitting')
    parser.add_argument('--mass-max', type=float, default=3.5, 
                       help='Maximum mass for fitting')
    parser.add_argument('--n-cores', type=int, default=1, 
                       help='Number of cores for parallel processing')
    parser.add_argument('--save-plots', default=True, action='store_true', 
                       help='Save fit plots')
    parser.add_argument('--force-convert', action='store_true', 
                       help='Force re-conversion of ROOT to HDF5')
    parser.add_argument('--force-analyze', action='store_true', 
                       help='Force re-analysis of HDF5 files')
    parser.add_argument('--no-analyze', default=False, action='store_true',
                       help='Skip analysis of HDF5 files')
    parser.add_argument('--no-convert', default=False, action='store_true',
                       help='Skip conversion of ROOT to HDF5')

    args = parser.parse_args()
    
    # Create output directories
    os.makedirs(args.hdf5_dir, exist_ok=True)
    os.makedirs(args.csv_dir, exist_ok=True)
    if args.save_plots:
        os.makedirs(args.plot_dir, exist_ok=True)
    
    # Check that required scripts exist
    required_scripts = ['root_to_hdf5.py', 'hdf5_ml_analysis.py']
    for script in required_scripts:
        if not os.path.exists(script):
            print(f"Error: Required script not found: {script}")
            print("Please ensure all pipeline scripts are in the current directory")
            sys.exit(1)
    
    # Process each year and era
    processed_count = 0
    failed_count = 0

    if not args.years:
        # Get all years from the configuration
        years = era_config.keys()
    else:
        years = args.years

    for year in years:
        # Get eras for this year
        if args.eras:
            eras = args.eras
        else:
            eras = era_config.get(year, [])
        
        if not eras:
            print(f"No eras configured for year {year}")
            continue
        
        base_path = base_path_pattern.format(year=year)
        
        for era in eras:
            success = process_era(year, era, base_path, args)
            if success:
                processed_count += 1
            else:
                failed_count += 1        
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Processing Complete")
    print(f"  Successfully processed: {processed_count} era(s)")
    print(f"  Failed: {failed_count} era(s)")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()