import uproot as up
import pandas as pd
import numpy as np
import ROOT

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
    # event = arrays["Event"][final_mask].astype(int)
    return run, lumiblock, mass

def read_lumi_data(lumi_file):
    """Read luminosity data from brilcalc CSV."""
    df = pd.read_csv(lumi_file, skiprows=1, engine='python', skipfooter=5, parse_dates=['time'])
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
    """
    Prepare data for analysis with option to use direct counting or ML fit.
    
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
        Min and max values for the mass range in GeV
    
    Returns:
    --------
    tuple
        Containing event counts, luminosity values, pileup values, bin labels,
        number of lumisections per bin, and prescale line positions
    """
    print("Using mass range", mass_range)
    prescale_points = prescale_points or []

    data = np.rec.fromarrays([run, lumiblock], names="run,lumiblock")
    data.sort(order=["run", "lumiblock"])

    # Create a recarray that includes mass if provided
    if mass is not None:
        data = np.rec.fromarrays([run, lumiblock, mass], names="run,lumiblock,mass")
        data.sort(order=["run", "lumiblock"])

    # Group data by binned lumiblocks
    bin_data = {}
    for i in range(len(data)):
        r, l = data.run[i], data.lumiblock[i]
        key = (r, l // bin_width)
        if key not in bin_data:
            bin_data[key] = []
        
        if mass is not None:
            bin_data[key].append(data.mass[i])

    # Create bin label mapping
    bin_label_map = {key: f"{key[0]}_{key[1]*bin_width}" for key in bin_data.keys()}

    # Process luminosity data
    lumi_sums = {}
    lumi_counts = {}
    pileup_sums = {}
    
    if lumi_df is not None:
        for _, row in lumi_df.iterrows():
            r, l = row['run'], row['ls']
            recorded = row['recorded(/ub)'] / 1000  # nb
            avgpu = row['avgpu']
            key = (r, l // bin_width)
            if key not in lumi_sums:
                lumi_sums[key] = 0
                pileup_sums[key] = 0
                lumi_counts[key] = 0
            lumi_sums[key] += recorded
            pileup_sums[key] += avgpu
            lumi_counts[key] += 1

    # Count events based on selected method
    event_counts = {}
    event_counts_error = {}
    
    if count_method == "direct":
        # Direct counting method
        for key, values in bin_data.items():
            values = [v for v in values if v > mass_range[0] and v < mass_range[1]]
            event_counts[key] = len(values)
            event_counts_error[key] = np.sqrt(len(values))  # Poisson error
    elif count_method == "fit" and mass is not None:
        # ML fit method
        for key, mass_values in bin_data.items():
            if len(mass_values) > 10000:  # Minimum events needed for a reasonable fit
                print(f"Performing ML fit for bin {key}")
                # plot_filename = None
                plot_filename = f"ML_fit_bin_{key[0]}_{key[1]*bin_width}.png"
                
                # Perform the fit with the specified mass range
                fit_count, fit_count_error = perform_ml_fit(
                    mass_values, 
                    mass_range=mass_range,
                    filename = plot_filename
                )
                event_counts[key] = fit_count
                event_counts_error[key] = fit_count_error
            else:
                event_counts[key] = -1
                event_counts_error[key] = -1
                print("Not enough statistcs for the ML fit")
    else:
        raise ValueError("Invalid count method or missing mass data for fitting")

    # Finalize bins
    all_keys = sorted(set(event_counts.keys()) | set(lumi_sums.keys() if lumi_df is not None else set()))
    event_y, event_y_error, lumi_y, pileup_y, x_labels, nlumis = [], [], [], [], [], []

    for key in all_keys:
        label = bin_label_map.get(key, f"{key[0]}_{key[1]*bin_width}")
        x_labels.append(label)
        event_y.append(event_counts.get(key, 0))
        event_y_error.append(event_counts_error.get(key, 0))
        
        if lumi_df is not None:
            lc = lumi_counts.get(key, 1)
            lumi_y.append(lumi_sums.get(key, 0) / lc)
            pileup_y.append(pileup_sums.get(key, 0) / lc)
            nlumis.append(lc)
        else:
            lumi_y.append(0)
            pileup_y.append(0)
            nlumis.append(0)

    # Handle prescale lines
    prescale_lines = {
        bin_label_map.get((r, l // bin_width))
        for r, l in prescale_points if bin_label_map.get((r, l // bin_width))
    }

    return event_y, event_y_error, lumi_y, pileup_y, x_labels, nlumis, sorted(prescale_lines)


def perform_ml_fit(mass_values, mass_range, filename=None):
    """
    Perform unbinned Maximum Likelihood fit on mass values using RooFit
    with double Gaussian signal and polynomial background
    
    Parameters:
    -----------
    mass_values : array-like
        Array of mass values to fit
    mass_range : tuple, default=(2.6, 3.6)
        Min and max values for the mass range
    
    Returns:
    --------
    float
        Number of signal events from the fit
    """
    if len(mass_values) < 10000:
        # Not enough statistics for a reliable fit
        return -1, -1
    
    # Create the mass variable
    mass_min, mass_max = mass_range
    mass_var = ROOT.RooRealVar("J mass", "J mass (GeV)", mass_min, mass_max)
    
    # Create unbinned dataset
    data_array = np.array(mass_values, dtype=float)
    data_list = ROOT.RooArgList()
    data_list.add(mass_var)
    dataset = ROOT.RooDataSet("data", "Unbinned dataset", data_list)
    
    # Fill the dataset
    for value in data_array:
        if mass_min < value < mass_max:  # Ensure value is within range
            mass_var.setVal(value)
            dataset.add(data_list)

    # ########## Model 1: Double Gaussian + Polynomial Background ##########
    # # Define double Gaussian signal model
    # mean = ROOT.RooRealVar("mean", "Mean of Gaussians", 3.1, 3.0, 3.2)
    # sigma1 = ROOT.RooRealVar("sigma1", "Width of Gaussian 1", 0.02, 0.001, 0.1)
    # sigma2 = ROOT.RooRealVar("sigma2", "Width of Gaussian 2", 0.05, 0.001, 0.1)
    # gauss1 = ROOT.RooGaussian("gauss1", "Gaussian 1", mass_var, mean, sigma1)
    # gauss2 = ROOT.RooGaussian("gauss2", "Gaussian 2", mass_var, mean, sigma2)
    # frac = ROOT.RooRealVar("frac", "Fraction of Gauss1", 0.5, 0.0, 1.0)
    # signal = ROOT.RooAddPdf("signal", "Double Gaussian", ROOT.RooArgList(gauss1, gauss2), ROOT.RooArgList(frac))
    
    # # Define background model (2nd degree polynomial)
    # a0 = ROOT.RooRealVar("a0", "a0", 0.0, -1.0, 1.0)
    # a1 = ROOT.RooRealVar("a1", "a1", 0.0, -1.0, 1.0)
    # background = ROOT.RooPolynomial("background", "polynomial", mass_var, ROOT.RooArgList(a0, a1))
    
    # # Combine signal and background
    # nsig = ROOT.RooRealVar("nsig", "Number of signal events", len(data_array)/2, 0, len(data_array)*2)
    # nbkg = ROOT.RooRealVar("nbkg", "Number of background events", len(data_array)/2, 0, len(data_array)*2)
    
    # model = ROOT.RooAddPdf("model", "Signal + Background", 
    #                      ROOT.RooArgList(signal, background),
    #                      ROOT.RooArgList(nsig, nbkg))
    # ############ END Model 1 ##########

    ########### Model 2: CB + Gaussian + Polynomial Background ##########
    mean = ROOT.RooRealVar("mean", "Mean", 3.09, 2.7, 3.5)
    sigma = ROOT.RooRealVar("sigma", "Sigma", 0.02, 0.001, 0.1)
    alpha = ROOT.RooRealVar("alpha", "Alpha (Tail Slope)", 2.0, 0.5, 5.0)
    n = ROOT.RooRealVar("n", "n (Tail Exponent)", 3, 1, 10)
    cb = ROOT.RooCrystalBall("cb", "Crystal Ball", mass_var, mean, sigma, alpha, n)
    
    sigma2 = ROOT.RooRealVar("sigma2", "Width of Gaussian 2", 0.04, 0.01, 0.1)
    gauss2 = ROOT.RooGaussian("gauss2", "Gaussian 2", mass_var, mean, sigma2)
    frac = ROOT.RooRealVar("frac", "Fraction of Gauss2", 0.7, 0.0, 1.0)
    signal = ROOT.RooAddPdf("signal", "CB + Gaussian", ROOT.RooArgList(cb, gauss2), ROOT.RooArgList(frac))

    # Define background model (2nd degree polynomial)
    c0 = ROOT.RooRealVar("c0", "c0", -0.1, -1, 0)
    # c1 = ROOT.RooRealVar("c1", "c1", 0.1, -5, 5)
    # c2 = ROOT.RooRealVar("c2", "c2", 0.05, -5, 5)
    background = ROOT.RooChebychev("background", "Polynomial Background", mass_var, ROOT.RooArgList(c0))

    # Combine signal and background
    nsig = ROOT.RooRealVar("nsig", "Number of signal events", len(data_array)/2, 0, len(data_array)*2)
    nbkg = ROOT.RooRealVar("nbkg", "Number of background events", len(data_array)/2, 0, len(data_array)*2)

    model = ROOT.RooAddPdf("model", "Signal + Background",
                         ROOT.RooArgList(signal, background),
                         ROOT.RooArgList(nsig, nbkg))
    ############# END Model 2 ##########
   
    
    # Perform the extended ML fit
    model.fitTo(dataset, ROOT.RooFit.Extended(), ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(-1))
    
    # Optional: Create and save plot if needed
    if filename is not None:
        create_fit_plot(mass_var, dataset, model, nsig.getVal(), nbkg.getVal(), filename=filename)
    
    # Return the number of signal events
    return nsig.getVal(), nsig.getError()


def create_fit_plot(mass_var, dataset, model, nsig, nbkg, filename=None):
    """
    Create and save a plot of the fit result
    
    Parameters:
    -----------
    mass_var : RooRealVar
        Mass variable
    dataset : RooDataSet
        Dataset with mass values
    model : RooAddPdf
        Combined signal+background model
    signal : RooAddPdf
        Signal component
    background : RooPolynomial
        Background component
    nsig : float
        Number of signal events
    nbkg : float
        Number of background events
    filename : str, optional
        Filename to save the plot to
    """
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
    
    return c1


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