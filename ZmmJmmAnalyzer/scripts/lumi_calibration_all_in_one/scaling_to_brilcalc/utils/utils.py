import uproot as up
import pandas as pd
import numpy as np
import ROOT
import os
import glob
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
import gc

def read_lumi_data(lumi_file):
    """Read luminosity data from brilcalc CSV."""
    df = pd.read_csv(lumi_file, skiprows=1, engine='python', skipfooter=5)#, parse_dates=['time'])
    df['time'] = pd.to_datetime(df['time'], format='%m/%d/%y %H:%M:%S')
    df[['run', 'fill']] = df['#run:fill'].str.split(':', expand=True)
    df = df.drop(columns=['#run:fill'])
    df = df[['run', 'fill', 'time', 'ls', 'beamstatus', 'recorded(/ub)', 'delivered(/ub)', 'avgpu']]
    df['run'] = df['run'].astype(int)
    df['ls'] = df['ls'].apply(lambda x: x.split(':')[0]).astype(int)
    return df


def weighted_average(values, errors):
    """
    Calculate the weighted average and its error.
    
    Parameters:
    values (array-like): The values to average (e.g., fit_slope)
    errors (array-like): The errors associated with each value (e.g., fit_slope_err)
    
    Returns:
    tuple: (weighted_average, error_on_weighted_average)
    """
    # Calculate weights as 1/error^2 (inverse variance weighting)
    weights = 1.0 / (errors**2)
    
    # Calculate weighted average
    weighted_avg = np.sum(values * weights) / np.sum(weights)
    
    # Calculate error on the weighted average
    weighted_avg_error = np.sqrt(1.0 / np.sum(weights))
    
    return weighted_avg, weighted_avg_error
