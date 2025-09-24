#!/usr/bin/env python3
"""
Stage 2: HDF5 ML Analysis
Reads HDF5 files and performs ML fits with multiprocessing
"""

import h5py
import pandas as pd
import numpy as np
import ROOT
import os
import json
import argparse
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm
import gc

# Suppress ROOT messages
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


def perform_ml_fit(mass_values, mass_range, bin_key=None, output_dir=None):
    """
    Perform ML fit on mass values for a single bin.
    
    Parameters
    ----------
    mass_values : array-like
        Mass values for fitting
    mass_range : tuple
        (min, max) mass range
    bin_key : str
        Identifier for this bin (for filename)
    output_dir : str
        Directory to save plots
    
    Returns
    -------
    dict
        Fit results including parameters and errors
    """
    
    if len(mass_values) < 10000:
        return {"fit_status": "insufficient_data"}

    try:
        mass_min, mass_max = mass_range
        mass_var = ROOT.RooRealVar("J_mass", "J mass (GeV)", mass_min, mass_max)
        
        # Create dataset
        dataset = ROOT.RooDataSet("data", "Unbinned dataset", ROOT.RooArgSet(mass_var))
        for value in mass_values:
            if mass_min < value < mass_max:
                mass_var.setVal(value)
                dataset.add(ROOT.RooArgSet(mass_var))
        
        # Signal: Crystal Ball + Gaussian
        mean = ROOT.RooRealVar("mean", "Mean", 3.09, 2.7, 3.5)
        sigma = ROOT.RooRealVar("sigma", "Sigma", 0.02, 0.001, 0.1)
        alpha = ROOT.RooRealVar("alpha", "Alpha", 2.0, 0.5, 5.0)
        n = ROOT.RooRealVar("n", "n", 3, 1, 10)
        sigma2 = ROOT.RooRealVar("sigma2", "Sigma2", 0.04, 0.01, 0.1)
        frac = ROOT.RooRealVar("frac", "Fraction", 0.7, 0.0, 1.0)
        
        cb = ROOT.RooCrystalBall("cb", "Crystal Ball", mass_var, mean, sigma, alpha, n)
        gauss2 = ROOT.RooGaussian("gauss2", "Gaussian 2", mass_var, mean, sigma2)
        signal = ROOT.RooAddPdf("signal", "CB + Gaussian", 
                               ROOT.RooArgList(cb, gauss2), 
                               ROOT.RooArgList(frac))
        
        # Background: Chebychev
        c0 = ROOT.RooRealVar("c0", "c0", -0.1, -1, 0)
        background = ROOT.RooChebychev("background", "Background", 
                                       mass_var, ROOT.RooArgList(c0))
        
        # Yields
        n_total = len(mass_values)
        nsig = ROOT.RooRealVar("nsig", "Signal yield", n_total/2, 0, n_total*2)
        nbkg = ROOT.RooRealVar("nbkg", "Background yield", n_total/2, 0, n_total*2)
        
        # Full model
        model = ROOT.RooAddPdf("model", "Signal + Background",
                              ROOT.RooArgList(signal, background),
                              ROOT.RooArgList(nsig, nbkg))
        
        # Perform fit
        fit_result = model.fitTo(
            dataset,
            ROOT.RooFit.Extended(),
            ROOT.RooFit.Save(),
            ROOT.RooFit.NumCPU(1),  # Single CPU per process
            ROOT.RooFit.PrintLevel(-1),
            ROOT.RooFit.Verbose(False),
            ROOT.RooFit.Warnings(False)
        )
        
        # Create plot if requested
        chi2 = -1
        if output_dir and bin_key:
            chi2 = create_fit_plot(mass_var, dataset, model, nsig.getVal(), 
                                 nbkg.getVal(), bin_key, output_dir)
        
        # Extract results
        result = {
            "fit_status": "success",
            "chi2": chi2,
        }
        
        for param in model.getParameters(ROOT.RooArgSet(mass_var)):
            result[param.GetName()] = param.getVal()
            result[param.GetName() + "_err"] = param.getError()
        
        return result
        
    except Exception as e:
        print(f"Error in ML fit for bin {bin_key}: {e}")
        return {"fit_status": "error"}


def create_fit_plot(mass_var, dataset, model, nsig, nbkg, bin_key, output_dir):
    """Create and save fit plot."""
    
    c1 = ROOT.TCanvas("c1", "Fit Result", 800, 600)
    frame = mass_var.frame(ROOT.RooFit.Title(f"Mass Fit - Bin {bin_key}"))
    
    dataset.plotOn(frame, ROOT.RooFit.Name('dataset'), ROOT.RooFit.Binning(50))
    model.plotOn(frame, ROOT.RooFit.Name('model'))
    model.plotOn(frame, ROOT.RooFit.Components("background"),
                ROOT.RooFit.LineStyle(ROOT.kDashed),
                ROOT.RooFit.LineColor(ROOT.kRed))
    model.plotOn(frame, ROOT.RooFit.Components("signal"),
                ROOT.RooFit.LineStyle(ROOT.kDashed),
                ROOT.RooFit.LineColor(ROOT.kGreen))
    
    chi2 = frame.chiSquare("model", "dataset", 7)
    
    model.paramOn(frame, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
    frame.getAttText().SetTextSize(0.03)
    
    pt = frame.findObject("model_paramBox")
    if pt:
        pt.AddText(ROOT.Form(f"Chi2/ndof = {chi2:.2f}"))
        pt.AddText(ROOT.Form(f"sigFrac = {nsig/(nsig + nbkg):.3f}"))
        pt.AddText(ROOT.Form(f"Entries = {dataset.numEntries()}"))
    
    frame.Draw()
    
    # Save plot
    run = bin_key.split('_')[0]
    plot_subdir = os.path.join(output_dir, run)
    os.makedirs(plot_subdir, exist_ok=True)
    
    filename = os.path.join(plot_subdir, f"ML_fit_bin_{bin_key}.png")
    c1.SaveAs(filename)
    
    del c1, frame
    
    return chi2


def process_single_bin(args):
    """Worker function for multiprocessing."""
    bin_key, hdf5_path, mass_range, save_plots, plot_dir = args
    
    with h5py.File(hdf5_path, 'r') as hf:
        bin_group = hf[f'bins/{bin_key}']
        masses = bin_group['events']['mass']
        npvs = bin_group['events']['npv']

        # Get bin metadata
        run = bin_group.attrs['run']
        lumiblock = bin_group.attrs['lumiblock']
    
    # Add bin information
    fit_result = {}
    fit_result['run'] = run
    fit_result['lumiblock'] = lumiblock
    fit_result['bin_key'] = bin_key
    fit_result['nPV'] = npvs.mean() if len(npvs) > 0 else -1
    fit_result['nPV_err'] = npvs.std() if len(npvs) > 1 else 0
    fit_result['totalEvents'] = len(masses)

    # Perform ML fit
    ml_fit_results = perform_ml_fit(
        masses,
        mass_range,
        bin_key=bin_key if save_plots else None,
        output_dir=plot_dir if save_plots else None
    )

    fit_result.update(ml_fit_results)
    
    return fit_result


def analyze_hdf5_with_ml(
    hdf5_path,
    lumi_df=None,
    mass_range=(2.7, 3.5),
    n_cores=1,
    save_plots=True,
    plot_dir="plots/ML_fits"
):
    """
    Analyze HDF5 file with ML fitting using multiprocessing.
    
    Parameters
    ----------
    hdf5_path : str
        Path to HDF5 file
    lumi_df : DataFrame, optional
        Luminosity data
    mass_range : tuple
        Mass range for fitting
    n_cores : int, optional
        Number of cores to use (None = all available)
    save_plots : bool
        Whether to save fit plots
    plot_dir : str
        Directory for plots
    
    Returns
    -------
    DataFrame
        Analysis results
    """
    
    print(f"Analyzing HDF5 file: {hdf5_path}")
    
    # Read metadata and get bin list
    with h5py.File(hdf5_path, 'r') as hf:
        metadata = dict(hf['metadata'].attrs)
        bin_width = metadata['bin_width']
        bin_keys = list(hf['bins'].keys())
        
        print(f"Found {len(bin_keys)} bins")
        print(f"Bin width: {bin_width}")
    
    # Setup multiprocessing
    n_cores = min(n_cores, len(bin_keys), cpu_count())
    print(f"Using {n_cores} cores for parallel processing")
    
    # Prepare arguments for workers
    worker_args = [
        (bin_key, hdf5_path, mass_range, save_plots, plot_dir)
        for bin_key in bin_keys
    ]
    
    # Process bins in parallel
    results = []
    with Pool(n_cores) as pool:
        with tqdm(total=len(bin_keys), desc="Processing bins") as pbar:
            for result in pool.imap_unordered(process_single_bin, worker_args):
                results.append(result)
                pbar.update(1)
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Add luminosity data if provided
    if lumi_df is not None:
        lumi_df['lumiblock'] = (lumi_df['ls'] // bin_width) * bin_width

        lumi_agg = lumi_df.groupby(['run', 'lumiblock']).agg(
            lumi=('recorded(/ub)', lambda x: x.mean() / 1000),  # ub -> nb
            pileup=('avgpu', 'mean'),
            nlumis=('ls', 'count')
        ).reset_index()

        df = df.merge(
            lumi_agg,
            on=['run', 'lumiblock'],
            how='inner'
        )
    
    # Sort by run and lumiblock
    df = df.sort_values(['run', 'lumiblock'])
    
    # Add event counts and errors
    df['events'] = df['nsig']
    df['events_error'] = df['nsig_err']

    # move some columns to front
    front_cols = ['run', 'lumiblock', 'bin_key',
                  'fit_status', 
                  'totalEvents', 'events', 'events_error', 
                  'nPV', 'nPV_err', 'lumi', 'pileup', 'nlumis', 'chi2']
    other_cols = [col for col in df.columns if col not in front_cols]
    df = df[front_cols + other_cols]

    print(f"Analysis complete. Processed {len(df)} bins.")
    
    # Print summary of fit results
    successful_fits = df[df['fit_status'] == 'success']
    print(f"Successful fits: {len(successful_fits)}/{len(df)}")
    
    return df


def read_lumi_data(lumi_file):
    """Read luminosity data from brilcalc CSV."""
    df = pd.read_csv(lumi_file, skiprows=1, engine='python', skipfooter=5)
    df['time'] = pd.to_datetime(df['time'], format='%m/%d/%y %H:%M:%S')
    df[['run', 'fill']] = df['#run:fill'].str.split(':', expand=True)
    df = df.drop(columns=['#run:fill'])
    df = df[['run', 'fill', 'time', 'ls', 'beamstatus', 'recorded(/ub)', 'delivered(/ub)', 'avgpu']]
    df['run'] = df['run'].astype(int)
    df['ls'] = df['ls'].apply(lambda x: x.split(':')[0]).astype(int)
    return df


def main():
    parser = argparse.ArgumentParser(description='Analyze HDF5 file with ML fitting')
    parser.add_argument('hdf5_file', help='Input HDF5 file path')
    parser.add_argument('output_csv', help='Output CSV file path')
    parser.add_argument('--lumi-file', help='Luminosity CSV file path')
    parser.add_argument('--mass-min', type=float, default=2.7, help='Minimum mass')
    parser.add_argument('--mass-max', type=float, default=3.5, help='Maximum mass')
    parser.add_argument('--n-cores', type=int, help='Number of cores to use')
    parser.add_argument('--save-plots', default=True, action='store_true', help='Save fit plots')
    parser.add_argument('--plot-dir', default='plots/ML_fits', help='Directory for plots')
    parser.add_argument('--year', help='Year label')
    parser.add_argument('--era', help='Era label')
    
    args = parser.parse_args()
    
    # Read luminosity data if provided
    lumi_df = None
    if args.lumi_file:
        print(f"Reading luminosity data from {args.lumi_file}")
        lumi_df = read_lumi_data(args.lumi_file)
    
    # Perform analysis
    df = analyze_hdf5_with_ml(
        args.hdf5_file,
        lumi_df=lumi_df,
        mass_range=(args.mass_min, args.mass_max),
        n_cores=args.n_cores,
        save_plots=args.save_plots,
        plot_dir=args.plot_dir
    )
    
    # Add year and era if provided
    if args.year:
        df['year'] = args.year
    if args.era:
        df['era'] = args.era
    
    # Save results
    os.makedirs(os.path.dirname(args.output_csv) if os.path.dirname(args.output_csv) else '.', exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    print(f"Results saved to {args.output_csv}")


if __name__ == "__main__":
    main()