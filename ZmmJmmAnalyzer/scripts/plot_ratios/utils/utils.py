import uproot as up
import pandas as pd
import numpy as np
import ROOT
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

def read_root_file(file_path, cuts):
    """Read ROOT file and apply cuts."""
    tree = up.open(file_path)["ntuple"]
    
    print(f"Available branches: {tree.keys()}")

    branches_to_read = ["Run", "LumiBlock", "B_J1_mass", "B_Mu1_pt", "B_Mu2_pt"]
    arrays = tree.arrays(branches_to_read, library="np")

    # Apply cuts (e.g., for mass and other properties)
    mask = np.ones(len(arrays["Run"]), dtype=bool)
    for var, (low, high) in cuts.items():
        if low is not None:
            mask &= (arrays[var] > low)
        if high is not None:
            mask &= (arrays[var] < high)

    # Combine trigger and cuts
    final_mask = mask

    # Apply the mask
    run = arrays["Run"][final_mask].astype(int)
    lumiblock = arrays["LumiBlock"][final_mask].astype(int)
    mass = arrays["B_J1_mass"][final_mask].astype(float)
    return run, lumiblock, mass

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

def prepare_data(run, lumiblock, mass_range, mass=None, lumi_df=None, prescale_points=None, bin_width=50, 
                 count_method="direct"):
    print("Using mass range", mass_range)
    prescale_points = prescale_points or []

    # Prepare structured data array
    names = "run,lumiblock" + (",mass" if mass is not None else "")
    arrays = [run, lumiblock] + ([mass] if mass is not None else [])
    data = np.rec.fromarrays(arrays, names=names)
    data.sort(order=["run", "lumiblock"])

    # Bin data by (run, lumiblock // bin_width)
    bin_data = {}
    for i in range(len(data)):
        key = (data.run[i], data.lumiblock[i] // bin_width)
        if key not in bin_data:
            bin_data[key] = []
        if mass is not None:
            bin_data[key].append(data.mass[i])

    bin_label_map = {key: f"{key[0]}_{key[1] * bin_width}" for key in bin_data}

    # Accumulate luminosity and pileup
    lumi_sums, pileup_sums, lumi_counts = {}, {}, {}
    if lumi_df is not None:
        for _, row in lumi_df.iterrows():
            key = (row['run'], row['ls'] // bin_width)
            lumi_sums[key] = lumi_sums.get(key, 0) + row['recorded(/ub)'] / 1000  # nb
            pileup_sums[key] = pileup_sums.get(key, 0) + row['avgpu']
            lumi_counts[key] = lumi_counts.get(key, 0) + 1

    event_counts, event_counts_error, fit_rows = {}, {}, []

    if count_method == "direct":
        print("Counting events directly")
        for key, values in bin_data.items():
            filtered = [v for v in values if mass_range[0] < v < mass_range[1]]
            count = len(filtered)
            event_counts[key] = count
            event_counts_error[key] = np.sqrt(count)

    elif count_method == "fit" and mass is not None:
        print("Performing ML fit")
        for key, mass_values in bin_data.items():
            if len(mass_values) <= 10000:
                print(f"bin {key} — Not enough statistics for the ML fit")
                event_counts[key] = -1
                event_counts_error[key] = -1
                continue

            print(f"bin {key}")
            plot_filename = f"ML_fit_bin_{key[0]}_{key[1] * bin_width}.png"
            fit_row = perform_ml_fit(mass_values, mass_range=mass_range, filename=plot_filename)

            if fit_row:
                fit_row["run"] = key[0]
                fit_row["lumiblock"] = key[1] * bin_width
                fit_rows.append(fit_row)
                event_counts[key] = fit_row.get("nsig", -1)
                event_counts_error[key] = fit_row.get("nsig_err", -1)
            else:
                event_counts[key] = -1
                event_counts_error[key] = -1

    else:
        raise ValueError("Invalid count method or missing mass for fit")

    # Assemble all rows
    all_keys = sorted(set(event_counts.keys()) | (set(lumi_sums) if lumi_df is not None else set()))
    rows = []
    for key in all_keys:
        label = bin_label_map.get(key, f"{key[0]}_{key[1] * bin_width}")
        lc = lumi_counts.get(key, 1)
        rows.append({
            'run': key[0],
            'lumiblock': key[1] * bin_width,
            'x_label': label,
            'event_y': event_counts.get(key, 0),
            'event_y_error': event_counts_error.get(key, 0),
            'lumi_y': lumi_sums.get(key, 0) / lc if lumi_df is not None else 0,
            'pileup_y': pileup_sums.get(key, 0) / lc if lumi_df is not None else 0,
            'nlumis': lc if lumi_df is not None else 0
        })

    df = pd.DataFrame(rows)

    if fit_rows:
        df_fit = pd.DataFrame(fit_rows)
        df = pd.merge(df, df_fit, on=["run", "lumiblock"], how="left")

    prescale_lines = {
        bin_label_map.get((r, l // bin_width))
        for r, l in prescale_points if bin_label_map.get((r, l // bin_width))
    }

    df.to_csv("event_data.csv", index=False)
    return df, sorted(prescale_lines)


def perform_ml_fit(mass_values, mass_range, filename=None):
    if len(mass_values) < 10000:
        raise ValueError("Not enough events for a reasonable fit")

    mass_min, mass_max = mass_range
    mass_var = ROOT.RooRealVar("J mass", "J mass (GeV)", mass_min, mass_max)

    # Fill dataset directly
    dataset = ROOT.RooDataSet("data", "Unbinned dataset", ROOT.RooArgSet(mass_var))
    for value in mass_values:
        if mass_min < value < mass_max:
            mass_var.setVal(value)
            dataset.add(ROOT.RooArgSet(mass_var))

    # Signal: Crystal Ball + Gaussian
    mean   = ROOT.RooRealVar("mean", "Mean", 3.09, 2.7, 3.5)
    sigma  = ROOT.RooRealVar("sigma", "Sigma", 0.02, 0.001, 0.1)
    alpha  = ROOT.RooRealVar("alpha", "Alpha (Tail Slope)", 2.0, 0.5, 5.0)
    n      = ROOT.RooRealVar("n", "n (Tail Exponent)", 3, 1, 10)
    sigma2 = ROOT.RooRealVar("sigma2", "Width of Gaussian 2", 0.04, 0.01, 0.1)
    frac   = ROOT.RooRealVar("frac", "Fraction of Gauss2", 0.7, 0.0, 1.0)

    cb     = ROOT.RooCrystalBall("cb", "Crystal Ball", mass_var, mean, sigma, alpha, n)
    gauss2 = ROOT.RooGaussian("gauss2", "Gaussian 2", mass_var, mean, sigma2)
    signal = ROOT.RooAddPdf("signal", "CB + Gaussian", ROOT.RooArgList(cb, gauss2), ROOT.RooArgList(frac))

    # Background: 1st order Chebychev
    c0 = ROOT.RooRealVar("c0", "c0", -0.1, -1, 0)
    background = ROOT.RooChebychev("background", "Background", mass_var, ROOT.RooArgList(c0))

    # Yields
    n_total = len(mass_values)
    nsig = ROOT.RooRealVar("nsig", "Signal yield", n_total/2, 0, n_total*2)
    nbkg = ROOT.RooRealVar("nbkg", "Background yield", n_total/2, 0, n_total*2)

    # Full model
    model = ROOT.RooAddPdf("model", "Signal + Background",
                           ROOT.RooArgList(signal, background),
                           ROOT.RooArgList(nsig, nbkg))

    # Fit
    fit_result = model.fitTo(
        dataset,
        ROOT.RooFit.Extended(),
        ROOT.RooFit.Save(),
        ROOT.RooFit.NumCPU(10),
        ROOT.RooFit.PrintLevel(-1),
        ROOT.RooFit.Verbose(False),
        ROOT.RooFit.Warnings(False)
    )

    # Optional: plot + chi2
    chi2 = None
    if filename:
        _, chi2 = create_fit_plot(mass_var, dataset, model, nsig.getVal(), nbkg.getVal(), filename=filename)

    # Prepare output dictionary
    row = {
        "chi2": chi2 if chi2 is not None else -1
    }

    for param in model.getParameters(ROOT.RooArgSet(mass_var)):
        row[param.GetName()] = param.getVal()
        row[param.GetName() + "_err"] = param.getError()

    return row


def create_fit_plot(mass_var, dataset, model, nsig, nbkg, filename=None):
    # Create canvas
    c1 = ROOT.TCanvas("c1", "Fit Result", 800, 600)
    
    # Create frame
    frame = mass_var.frame(ROOT.RooFit.Title("Mass Fit"))
    
    # Plot data and fit
    dataset.plotOn(frame, ROOT.RooFit.Name('dataset'), ROOT.RooFit.Binning(50))
    model.plotOn(frame, ROOT.RooFit.Name('model'))
    model.plotOn(frame, ROOT.RooFit.Components("background"), 
                ROOT.RooFit.LineStyle(ROOT.kDashed), 
                ROOT.RooFit.LineColor(ROOT.kRed))
    model.plotOn(frame, ROOT.RooFit.Components("signal"), 
                ROOT.RooFit.LineStyle(ROOT.kDashed), 
                ROOT.RooFit.LineColor(ROOT.kGreen))
    
    # Calculate chi2
    chi2 = frame.chiSquare("model", "dataset", 7)
    
    # Add stats box
    model.paramOn(frame, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
    frame.getAttText().SetTextSize(0.03)
    frame.GetYaxis().SetTitleSize(0.03)
    frame.GetYaxis().SetTitleOffset(1.7)
    
    # Find parameter box and add additional text
    pt = frame.findObject("model_paramBox")
    if pt:
        pt.AddText(ROOT.Form(f"Chi2/ndof = {chi2:.2f}"))
        pt.AddText(ROOT.Form(f"sigFrac = {nsig/(nsig + nbkg):.2f}"))
        pt.AddText(ROOT.Form(f"Entries = {dataset.numEntries()}"))
    
    # Draw frame
    frame.Draw()
    
    # Save plot if filename is provided
    if filename:
        c1.SaveAs(f"plots/ML_fits/{filename}")
    
    return c1, chi2


def analyze_data(run, lumiblock, mass=None, lumi_df=None, prescale_points=None, bin_width=50, 
                count_method="direct", mass_range=(2.6, 3.6)):
    """
    Wrapper function to analyze data with different counting methods
    
    Parameters:
    -----------
    run : array-like
        Array of run numbers
    lumiblock : array-like
        Array of lumiblock numbers
    mass : array-like, optional
        Array of mass values for ML fitting
    lumi_df : DataFrame
        DataFrame containing luminosity information
    prescale_points : list, optional
        List of (run, lumiblock) tuples where prescale changes occurred
    bin_width : int, default=50
        Width of lumiblock bins
    count_method : str, default="direct"
        Method for counting events: "direct" or "fit"
    mass_range : tuple, default=(2.6, 3.6)
        Min and max values for the mass range
    
    Returns:
    --------
    tuple
        Containing event counts, luminosity values, pileup values, bin labels,
        number of lumisections per bin, and prescale line positions
    """
    return prepare_data(
        run, lumiblock, mass, lumi_df, prescale_points, bin_width, 
        count_method, mass_range
    )