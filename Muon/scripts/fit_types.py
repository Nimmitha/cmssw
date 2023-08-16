import ROOT

def fit_unbinned_double_gauss(data, Z_mass):
    # Define the mean, sigma, and background parameters for the first Gaussian
    mean1 = ROOT.RooRealVar("mean1", "Mean1", 91.2, 85, 95)
    sigma1 = ROOT.RooRealVar("sigma1", "Sigma1", 2, 0.1, 10)

    # Define the mean, sigma, and background parameters for the second Gaussian
    mean2 = ROOT.RooRealVar("mean2", "Mean2", 91.2, 85, 95)
    sigma2 = ROOT.RooRealVar("sigma2", "Sigma2", 2, 0.1, 10)

    # Define the fractions of the first and second Gaussian
    f_sig1 = ROOT.RooRealVar("f_sig1", "Signal Fraction 1", 0.5, 0, 1)
    f_sig2 = ROOT.RooRealVar("f_sig2", "Signal Fraction 2", 0.5, 0, 1)

    # Create the first Gaussian component
    gaussian1 = ROOT.RooGaussian("gaussian1", "Gaussian1", Z_mass, mean1, sigma1)

    # Create the second Gaussian component
    gaussian2 = ROOT.RooGaussian("gaussian2", "Gaussian2", Z_mass, mean2, sigma2)

    # Combine the Gaussian components with their corresponding fractions
    # signal_model = ROOT.RooAddPdf("signal_model", "Signal Model", ROOT.RooArgList(gaussian1, gaussian2), ROOT.RooArgList(f_sig1, f_sig2))
    model = ROOT.RooAddPdf("model", "Signal Model", ROOT.RooArgList(gaussian1, gaussian2), ROOT.RooArgList(f_sig1))

    # # Define the background parameters
    # background = ROOT.RooRealVar("background", "Background", -0.1, -10, 10)

    # # Create the background model
    # background_model = ROOT.RooPolynomial("background_model", "Background Model", Z_mass, ROOT.RooArgList(background))

    # # Define the overall signal and background fractions
    # f_sig = ROOT.RooRealVar("f_sig", "Total Signal Fraction", 0.5, 0, 1)
    # f_bkg = ROOT.RooRealVar("f_bkg", "Total Background Fraction", 0.5, 0, 1)

    # # Combine the signal and background models with their corresponding fractions
    # model = ROOT.RooAddPdf("model", "Total Model", ROOT.RooArgList(signal_model, background_model), ROOT.RooArgList(f_sig, f_bkg))

    # Perform the fit to the data
    model.fitTo(data, ROOT.RooFit.Save())

    # Plot data
    frame = Z_mass.frame()
    data.plotOn(frame, ROOT.RooFit.Name('dataset'), ROOT.RooFit.Binning(70))

    # Plot the total fit
    model.plotOn(frame, ROOT.RooFit.Name('model'))

    # Plot the background fit
    model.plotOn(frame, ROOT.RooFit.Components("background_model"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen))

    # Plot the first Gaussian fit
    model.plotOn(frame, ROOT.RooFit.Components("gaussian1"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

    # Plot the second Gaussian fit
    model.plotOn(frame, ROOT.RooFit.Components("gaussian2"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue))

    # Calculate chi2
    chi2 = frame.chiSquare("model", "dataset", 5)

    # Stat box
    model.paramOn(frame, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
    frame.getAttText().SetTextSize(0.03)
    pt = frame.findObject("model_paramBox")
    pt.AddText(ROOT.Form(f"Chi2/ndof =  {chi2:.2f}"))

    return frame

    
def fit_unbinned_gauss(data, Z_mass):
    # Define the mean, sigma, and background parameters
    mean = ROOT.RooRealVar("mean", "Mean", 91.2, 85, 95)
    sigma = ROOT.RooRealVar("sigma", "Sigma", 2, 0.1, 10)
    background = ROOT.RooRealVar("background", "Background", -0.1, -10, 10)

    # Create the Gaussian and background models
    gaussian = ROOT.RooGaussian("gaussian", "Gaussian", Z_mass, mean, sigma)
    background_model = ROOT.RooPolynomial("background_model", "Background Model", Z_mass, ROOT.RooArgList(background))

    f_sig = ROOT.RooRealVar("f_sig", "Signal Fraction", 0.5, 0, 1)
    f_bkg = ROOT.RooRealVar("f_bkg", "Background Fraction", 0.5, 0, 1)

    # Combine the models
    model = ROOT.RooAddPdf("model", "Total Model", ROOT.RooArgList(gaussian, background_model), ROOT.RooArgList(f_sig, f_bkg)) 

    # Perform the fit to the data
    model.fitTo(data, ROOT.RooFit.Save())

    # Plot data
    frame = Z_mass.frame()
    data.plotOn(frame, ROOT.RooFit.Name('dataset'), ROOT.RooFit.Binning(70))

    # Plot the total fit
    model.plotOn(frame, ROOT.RooFit.Name('model'))

    # Plot the background fit
    model.plotOn(frame, ROOT.RooFit.Components("background_model"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen))

    # Plot the signal fit
    model.plotOn(frame, ROOT.RooFit.Components("gaussian"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

    # Calculate chi2
    chi2 = frame.chiSquare("model", "dataset", 5)

    # Stat box
    # frame.SetStats(ROOT.kTRUE)
    model.paramOn(frame, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
    frame.getAttText().SetTextSize(0.03)
    pt = frame.findObject("model_paramBox")
    pt.AddText(ROOT.Form(f"Chi2/ndof =  {chi2:.2f}"))

    return frame