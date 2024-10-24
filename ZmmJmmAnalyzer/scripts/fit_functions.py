import ROOT


def fit_unbinned_double_gauss_Z(data, Zcand_mass, nbins):
    # Define the mean, sigma, and background parameters for the first Gaussian
    mean1 = ROOT.RooRealVar("mean1", "Mean1", 91.2, 80, 100)
    sigma1 = ROOT.RooRealVar("sigma1", "Sigma1", 2, 0.1, 50)

    # Define the mean, sigma, and background parameters for the second Gaussian
    # mean2 = ROOT.RooRealVar("mean2", "Mean2", 91.2, 80, 100)
    sigma2 = ROOT.RooRealVar("sigma2", "Sigma2", 2, 0.1, 50)

    # Define the fractions of the first and second Gaussian
    f_sig1 = ROOT.RooRealVar("f_sig1", "Signal Fraction 1", 0.5, 0, 1)
    f_sig2 = ROOT.RooRealVar("f_sig2", "Signal Fraction 2", 0.5, 0, 1)

    # Create the first Gaussian component
    gaussian1 = ROOT.RooGaussian("gaussian1", "Gaussian1", Zcand_mass, mean1, sigma1)

    # Create the second Gaussian component
    gaussian2 = ROOT.RooGaussian("gaussian2", "Gaussian2", Zcand_mass, mean1, sigma2)

    # Combine the Gaussian components with their corresponding fractions
    # signal_model = ROOT.RooAddPdf("signal_model", "Signal Model", ROOT.RooArgList(gaussian1, gaussian2), ROOT.RooArgList(f_sig1, f_sig2))
    model = ROOT.RooAddPdf("model", "Signal Model", ROOT.RooArgList(gaussian1, gaussian2), ROOT.RooArgList(f_sig1))

    # # Define the background parameters
    # background = ROOT.RooRealVar("background", "Background", -0.1, -10, 10)

    # # Create the background model
    # background_model = ROOT.RooPolynomial("background_model", "Background Model", Zcand_mass, ROOT.RooArgList(background))

    # # Define the overall signal and background fractions
    # f_sig = ROOT.RooRealVar("f_sig", "Total Signal Fraction", 0.5, 0, 1)
    # f_bkg = ROOT.RooRealVar("f_bkg", "Total Background Fraction", 0.5, 0, 1)

    # # Combine the signal and background models with their corresponding fractions
    # model = ROOT.RooAddPdf("model", "Total Model", ROOT.RooArgList(signal_model, background_model), ROOT.RooArgList(f_sig, f_bkg))

    # Perform the fit to the data
    model.fitTo(data, ROOT.RooFit.Save())

    # Plot data
    frame = Zcand_mass.frame()
    data.plotOn(frame, ROOT.RooFit.Name('dataset'), ROOT.RooFit.Binning(nbins))

    # Plot the total fit
    model.plotOn(frame, ROOT.RooFit.Name('model'))

    # Plot the background fit
    model.plotOn(frame, ROOT.RooFit.Components("background_model"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen))

    # Plot the first Gaussian fit
    model.plotOn(frame, ROOT.RooFit.Components("gaussian1"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

    # Plot the second Gaussian fit
    model.plotOn(frame, ROOT.RooFit.Components("gaussian2"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue))

    # Calculate chi2
    chi2 = frame.chiSquare("model", "dataset", 4)

    # Stat box
    model.paramOn(frame, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
    frame.getAttText().SetTextSize(0.03)
    pt = frame.findObject("model_paramBox")
    pt.AddText(ROOT.Form(f"Chi2/ndof =  {chi2:.2f}"))
    pt.AddText(ROOT.Form("Entries =  {0:.0f}".format(data.numEntries())))

    return frame


def fit_unbinned_gauss_Z_with_background(data, Zcand_mass, nbins):
    # Define the mean, sigma, and background parameters
    mean = ROOT.RooRealVar("mean", "Mean", 91, 80, 100)
    sigma = ROOT.RooRealVar("sigma", "Sigma", 2, 0.1, 10)

    # create a polynomial of degree 2 for the background
    a0 = ROOT.RooRealVar("a0", "a0", 0.1, -10, 10)
    a1 = ROOT.RooRealVar("a1", "a1", 0.1, -10, 10)
    a2 = ROOT.RooRealVar("a2", "a2", 0.1, -10, 10)

    # Create the Gaussian and background models
    gaussian = ROOT.RooGaussian("gaussian", "Gaussian", Zcand_mass, mean, sigma)
    # background = ROOT.RooChebychev("background", "Background", Zcand_mass, ROOT.RooArgList(a0, a1, a2))
    background = ROOT.RooPolynomial("background", "Background", Zcand_mass, ROOT.RooArgList(a0, a1, a2))

    # Sum the composite signal and background models
    nbkg = ROOT.RooRealVar("nbkg", "number of background events", 10000, 0, 25000)
    nsig = ROOT.RooRealVar("nsig", "number of signal events", 15000, 0, 25000)

    # Combine the models
    model = ROOT.RooAddPdf("model", "g1 + chebychev", ROOT.RooArgList(background, gaussian), ROOT.RooArgList(nbkg, nsig))

    # Perform the fit to the data
    model.fitTo(data, ROOT.RooFit.Save())

    # Plot data
    frame = Zcand_mass.frame(ROOT.RooFit.Title(""))
    data.plotOn(frame, ROOT.RooFit.Name('dataset'), ROOT.RooFit.Binning(nbins))

    # Plot the total fit
    model.plotOn(frame, ROOT.RooFit.Name('model'))

    # Plot the background fit
    model.plotOn(frame, ROOT.RooFit.Components("background"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen))

    # Plot the signal fit
    model.plotOn(frame, ROOT.RooFit.Components("gaussian"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

    # Calculate chi2
    chi2 = frame.chiSquare("model", "dataset", 7)

    # Stat box
    model.paramOn(frame, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
    frame.getAttText().SetTextSize(0.03)
    pt = frame.findObject("model_paramBox")
    pt.AddText(ROOT.Form(f"Chi2/ndof =  {chi2:.2f}"))
    pt.AddText(ROOT.Form("Entries =  {0:.0f}".format(data.numEntries())))

    return frame


def fit_unbinned_double_gauss_J(data, Zcand_mass, nbins):
    # Define the mean, sigma, and background parameters for the first Gaussian
    mean1 = ROOT.RooRealVar("mean1", "Mean1", 3.1, 2.5, 3.6)
    sigma1 = ROOT.RooRealVar("sigma1", "Sigma1", 0.01, 0.001, 0.1)

    # Define the mean, sigma, and background parameters for the second Gaussian
    # mean2 = ROOT.RooRealVar("mean2", "Mean2", 91.2, 80, 100)
    sigma2 = ROOT.RooRealVar("sigma2", "Sigma2", 0.01, 0.001, 0.1)

    # Define the fractions of the first and second Gaussian
    f_sig1 = ROOT.RooRealVar("f_sig1", "Signal Fraction 1", 0.5, 0, 1)
    f_sig2 = ROOT.RooRealVar("f_sig2", "Signal Fraction 2", 0.5, 0, 1)

    # Create the first Gaussian component
    gaussian1 = ROOT.RooGaussian("gaussian1", "Gaussian1", Zcand_mass, mean1, sigma1)

    # Create the second Gaussian component
    gaussian2 = ROOT.RooGaussian("gaussian2", "Gaussian2", Zcand_mass, mean1, sigma2)

    # Combine the Gaussian components with their corresponding fractions
    # signal_model = ROOT.RooAddPdf("signal_model", "Signal Model", ROOT.RooArgList(gaussian1, gaussian2), ROOT.RooArgList(f_sig1, f_sig2))
    model = ROOT.RooAddPdf("model", "Signal Model", ROOT.RooArgList(gaussian1, gaussian2), ROOT.RooArgList(f_sig1))

    # # Define the background parameters
    # background = ROOT.RooRealVar("background", "Background", -0.1, -10, 10)

    # # Create the background model
    # background_model = ROOT.RooPolynomial("background_model", "Background Model", Zcand_mass, ROOT.RooArgList(background))

    # # Define the overall signal and background fractions
    # f_sig = ROOT.RooRealVar("f_sig", "Total Signal Fraction", 0.5, 0, 1)
    # f_bkg = ROOT.RooRealVar("f_bkg", "Total Background Fraction", 0.5, 0, 1)

    # # Combine the signal and background models with their corresponding fractions
    # model = ROOT.RooAddPdf("model", "Total Model", ROOT.RooArgList(signal_model, background_model), ROOT.RooArgList(f_sig, f_bkg))

    # Perform the fit to the data
    model.fitTo(data, ROOT.RooFit.Save())

    # Plot data
    frame = Zcand_mass.frame()
    data.plotOn(frame, ROOT.RooFit.Name('dataset'), ROOT.RooFit.Binning(nbins))

    # Plot the total fit
    model.plotOn(frame, ROOT.RooFit.Name('model'))

    # Plot the background fit
    model.plotOn(frame, ROOT.RooFit.Components("background_model"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen))

    # Plot the first Gaussian fit
    model.plotOn(frame, ROOT.RooFit.Components("gaussian1"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

    # Plot the second Gaussian fit
    model.plotOn(frame, ROOT.RooFit.Components("gaussian2"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue))

    # Calculate chi2
    chi2 = frame.chiSquare("model", "dataset", 4)

    # Stat box
    model.paramOn(frame, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
    frame.getAttText().SetTextSize(0.03)
    pt = frame.findObject("model_paramBox")
    pt.AddText(ROOT.Form(f"Chi2/ndof =  {chi2:.2f}"))
    pt.AddText(ROOT.Form("Entries =  {0:.0f}".format(data.numEntries())))

    return frame


def fit_unbinned_gauss_J_with_background(data, Zcand_mass, nbins):
    # Define the mean, sigma, and background parameters
    mean = ROOT.RooRealVar("mean", "Mean", 3.1, 3.0, 3.2)
    sigma = ROOT.RooRealVar("sigma", "Sigma", 0.03, 0.01, 0.1)

    # create a polynomial of degree 2 for the background
    a0 = ROOT.RooRealVar("a0", "a0", 0.01, -10, 10)
    a1 = ROOT.RooRealVar("a1", "a1", 0.01, -10, 10)
    a2 = ROOT.RooRealVar("a2", "a2", 0.01, -10, 10)

    # Create the Gaussian and background models
    gaussian = ROOT.RooGaussian("gaussian", "Gaussian", Zcand_mass, mean, sigma)
    # background = ROOT.RooChebychev("background", "Background", Zcand_mass, ROOT.RooArgList(a0, a1, a2))
    background = ROOT.RooPolynomial("background", "Background", Zcand_mass, ROOT.RooArgList(a0, a1, a2))

    # Sum the composite signal and background models
    nbkg = ROOT.RooRealVar("nbkg", "number of background events", 10000, 0, 25000)
    nsig = ROOT.RooRealVar("nsig", "number of signal events", 15000, 0, 25000)

    # Combine the models
    model = ROOT.RooAddPdf("model", "g1 + chebychev", ROOT.RooArgList(background, gaussian), ROOT.RooArgList(nbkg, nsig))

    # Perform the fit to the data
    model.fitTo(data, ROOT.RooFit.Save())

    # Plot data
    frame = Zcand_mass.frame(ROOT.RooFit.Title(""))
    data.plotOn(frame, ROOT.RooFit.Name('dataset'), ROOT.RooFit.Binning(nbins))

    # Plot the total fit
    model.plotOn(frame, ROOT.RooFit.Name('model'))

    # Plot the background fit
    model.plotOn(frame, ROOT.RooFit.Components("background"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen))

    # Plot the signal fit
    model.plotOn(frame, ROOT.RooFit.Components("gaussian"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

    # Calculate chi2
    chi2 = frame.chiSquare("model", "dataset", 7)

    # Stat box
    model.paramOn(frame, ROOT.RooFit.Layout(0.6, 0.9, 0.9))
    frame.getAttText().SetTextSize(0.03)
    pt = frame.findObject("model_paramBox")
    pt.AddText(ROOT.Form(f"Chi2/ndof =  {chi2:.2f}"))
    pt.AddText(ROOT.Form("Entries =  {0:.0f}".format(data.numEntries())))

    return frame
