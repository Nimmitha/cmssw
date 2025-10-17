import ROOT
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def fit_gaussian(data, bins=50, ax=None, plot=True):
    """
    Fit a Gaussian to a 1D NumPy array of data.
    Can optionally plot into a provided matplotlib Axes.

    Returns:
    amplitude, mean, sigma, (x_centers, y_hist, popt)
    """
    def gaussian(x, amplitude, mean, sigma):
        return amplitude * np.exp(-(x - mean)**2 / (2 * sigma**2))
    
    # Histogram data
    y_hist, bin_edges = np.histogram(data, bins=bins)
    x_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    x_fit_centers = np.linspace(min(x_centers), max(x_centers), 200)

    # Initial guess
    initial_guess = [max(y_hist), np.mean(data), np.std(data)]
    
    # Fit
    popt, _ = curve_fit(gaussian, x_centers, y_hist, p0=initial_guess)
    amplitude, mean, sigma = popt
    
    # Select where to plot
    if plot:
        if ax is None:
            fig, ax = plt.subplots()
        ax.hist(data, bins=bins, alpha=0.6, label='Data')
        ax.plot(x_fit_centers, gaussian(x_fit_centers, *popt), 'r-', label='Gaussian Fit')
        # ax.set_xlabel('Ratio')
        # ax.set_ylabel('Counts')
        # ax.set_title(f'Gaussian Fit: μ={mean:.3f}, σ={sigma:.3f}')
        ax.text(0.2, 0.95, f'μ={mean:.3f}\nσ={sigma:.3f}', 
                horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        # draw mean line and sigma shaded area
        ax.axvline(mean, color='k', linestyle='--', label='Mean')
        ax.fill_betweenx([0, max(y_hist)], mean - sigma, mean + sigma, color='gray', alpha=0.2, label='1σ range')
    
    return amplitude, mean, sigma#, (x_centers, y_hist, popt)


def root_polyfit(x, y, degree=1):
    """
    Perform a polynomial fit using ROOT (TGraph.Fit).

    Parameters
    ----------
    x : array-like
        Independent variable
    y : array-like
        Dependent variable
    degree : int
        Degree of polynomial (1, 2, ...)

    Returns
    -------
    coeffs : list
        Polynomial coefficients [highest degree first]
    errors : list
        Errors of the coefficients
    fit_result : ROOT.TFitResultPtr
        Full fit result (contains covariance matrix, chi2, etc.)
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    # Build TGraph
    npoints = len(x)
    graph = ROOT.TGraph(npoints, x, y)

    # Define polynomial function
    func = ROOT.TF1("poly", f"pol{degree}", float(min(x)), float(max(x)))

    # Perform fit
    fit_result = graph.Fit(func, "S")  # "S" = return full fit result

    # Extract coefficients and errors
    coeffs = []
    errors = []
    for i in range(degree + 1):
        coeffs.append(func.GetParameter(i))
        errors.append(func.GetParError(i))

    # Match np.polyval convention (highest degree first)
    coeffs = coeffs[::-1]
    errors = errors[::-1]

    # Optional: plot
    if True:
        c = ROOT.TCanvas("c", "c", 800, 600)
        graph.Draw("AP")
        func.Draw("same")
        c.Draw()
        c.SaveAs("root_polyfit.png")
        # input("Press Enter to continue...")  # Pause to view the plot
        c.Close()

    return coeffs, errors, fit_result



def roofit_polyfit(x, y, degree=1):
    """
    Perform a maximum likelihood polynomial fit using RooFit.

    Parameters
    ----------
    x : array-like
        Independent variable
    y : array-like
        Dependent variable
    degree : int
        Degree of polynomial (1 or 2)

    Returns
    -------
    coeffs : list
        Best-fit polynomial coefficients (highest degree first, like np.polyval)
    errors : list
        Errors of the coefficients
    """
    # Define RooRealVars
    xvar = ROOT.RooRealVar("x", "x", float(min(x)), float(max(x)))
    yvar = ROOT.RooRealVar("y", "y", float(min(y)), float(max(y)))

    # Put data into a RooDataSet
    data = ROOT.RooDataSet("data", "data",
                           ROOT.RooArgSet(xvar, yvar))
    for xi, yi in zip(x, y):
        xvar.setVal(float(xi))
        yvar.setVal(float(yi))
        data.add(ROOT.RooArgSet(xvar, yvar))

    # Build polynomial model
    if degree == 1:
        a0 = ROOT.RooRealVar("a0", "intercept", 1, -1, 2)
        a1 = ROOT.RooRealVar("a1", "slope", 0, -1e6, 1e6)
        model = ROOT.RooPolynomial("model", "linear", xvar, ROOT.RooArgList(a1, a0), 1)
    elif degree == 2:
        a0 = ROOT.RooRealVar("a0", "const", 0, -1e6, 1e6)
        a1 = ROOT.RooRealVar("a1", "linear", 0, -1e6, 1e6)
        a2 = ROOT.RooRealVar("a2", "quadratic", 0, -1e6, 1e6)
        model = ROOT.RooPolynomial("model", "quadratic", xvar,
                                   ROOT.RooArgList(a2, a1, a0), 2)
    else:
        raise ValueError("Only degree 1 or 2 supported")

    # Fit to data
    fit_result = model.fitTo(data, ROOT.RooFit.Save(True))

    # plot data and fit
    if True:
        frame = xvar.frame()
        data.plotOn(frame)
        model.plotOn(frame)
        c = ROOT.TCanvas("c", "c", 800, 600)
        frame.Draw()
        c.Draw()
        c.SaveAs("roofit_polyfit.png")
        # input("Press Enter to continue...")  # Pause to view the plot
        c.Close()

    # Extract coefficients and errors
    coeffs = []
    errors = []
    for par in [a2, a1, a0] if degree == 2 else [a1, a0]:
        coeffs.append(par.getVal())
        errors.append(par.getError())

    return coeffs, errors, fit_result



def mle_polyfit(x, y, degree=1):
    """
    Perform maximum likelihood estimation for polynomial fit.
    
    Parameters:
        x : array-like, independent variable
        y : array-like, dependent variable
        degree : int, degree of polynomial (1 or 2)
    
    Returns:
        coeffs : fitted polynomial coefficients
    """
    x = np.asarray(x)
    y = np.asarray(y)
    
    # Define negative log-likelihood
    def nll(params):
        coeffs = params[:-1]
        sigma = params[-1]
        model = np.polyval(coeffs, x)
        residuals = y - model
        # Gaussian log-likelihood
        return 0.5 * np.sum(np.log(2 * np.pi * sigma**2) + (residuals**2) / sigma**2)
    
    # Initial guess: polyfit + std of residuals
    init_coeffs = np.polyfit(x, y, degree)
    init_sigma = np.std(y - np.polyval(init_coeffs, x))
    init_params = np.append(init_coeffs, init_sigma)
    
    # Optimize
    result = minimize(nll, init_params, bounds=[(None, None)]*len(init_coeffs) + [(1e-6, None)])
    
    coeffs = result.x[:-1]
    sigma = result.x[-1]
    return coeffs, sigma