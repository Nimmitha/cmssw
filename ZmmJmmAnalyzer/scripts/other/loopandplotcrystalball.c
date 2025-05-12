#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <TH1F.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <cmath>
#include <TLatex.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooFit.h>
#include <fstream>

void loopandplotcrystalball() {
  // Open the ROOT file
  TFile *file = TFile::Open("output_data.root");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Unable to open file!" << std::endl;
    return;
  }

  // Select desired Tree in ROOT file
  TTree *fTree = (TTree *)file->Get("ntuple;1");
  if (!fTree) {
    std::cerr << "Error: TTree 'ntuple' not found in file!" << std::endl;
    return;
  }

  // Set branch addresses
  float Run;
  float LumiBlock;
  float B_J1_mass;
  float B_Mu1_pt;
  float B_Mu2_pt;
  float B_J1_VtxProb;

  fTree->SetBranchAddress("Run", &Run);
  fTree->SetBranchAddress("B_J1_mass", &B_J1_mass);
  fTree->SetBranchAddress("LumiBlock", &LumiBlock);
  fTree->SetBranchAddress("B_Mu1_pt", &B_Mu1_pt);
  fTree->SetBranchAddress("B_Mu2_pt", &B_Mu2_pt);

  // Define LumiBlock intervals
  std::vector<std::pair<int, int>> lumiIntervals = {
      {150, 200}  //choice
  };

  std::vector<double> storedmasses;
  // Loop over LumiBlock intervals
  for (size_t k = 0; k < lumiIntervals.size(); k++) {
    int lumiMin = lumiIntervals[k].first;
    int lumiMax = lumiIntervals[k].second;

    // Create a histogram for this LumiBlock range
    TH1D *hist = new TH1D(Form("hist_%zu", k), Form("Fill 9043 J/#psi Mass (LumiSections %d-%d) Crystal Ball + Gaussian + Poly 2", lumiMin, lumiMax), 100, 2.7, 3.5);
    hist->GetYaxis()->SetTitle("Events / [0.008 GeV/c^{2}]");
    hist->GetXaxis()->SetTitle("M(#mu^{+}#mu^{-}) [GeV/c^{2}]");
    hist->SetFillColor(kOrange + 6);

    // Loop through entries and fill histogram
    Long64_t ntotEntries = fTree->GetEntries();
    int globalLumiBlock;
    for (Long64_t i = 0; i < ntotEntries; i++) {
      if (i % 1000000 == 0)
        std::cout << "Processing entry: " << i << std::endl;
      fTree->GetEntry(i);

      if (B_Mu1_pt > 7 && B_Mu2_pt > 7) {
        storedmasses.push_back(B_J1_mass);
        hist->Fill(B_J1_mass);
      }
    }
    // Define mass variable
    RooRealVar mass("mass", "M(#mu^{+}#mu^{-}) [GeV/c^{2}]", 2.7, 3.5);
    RooDataHist data("data", "J/#psi mass dataset", RooArgSet(mass), hist);

    // Crystal Ball parameters
    RooRealVar mean("mean", "Mean", 3.09, 2.7, 3.5);
    RooRealVar sigma("sigma", "Sigma", 0.02, 0.001, 0.1);
    RooRealVar alpha("alpha", "Alpha (Tail Slope)", 2.0, 0.5, 5.0);
    RooRealVar n("n", "n (Tail Exponent)", 3, 1, 10);
    RooCrystalBall cb("cb", "Crystal Ball", mass, mean, sigma, alpha, n);

    // Second Gaussian to refine peak shape
    RooRealVar sigma2("sigma2", "Second Gaussian width", 0.04, 0.01, 0.1);
    RooGaussian gauss2("gauss2", "Second Gaussian", mass, mean, sigma2);

    // Combine Crystal Ball + Gaussian
    RooRealVar fracSig("fracSig", "Fraction of CB", 0.7, 0, 1);
    RooAddPdf signal("signal", "CB + Gaussian", RooArgList(cb, gauss2), RooArgList(fracSig));

    // chebychev polynomial for background
    RooRealVar c0("c0", "Background constant", -0.1, -1., 0.);
    RooChebychev background("background", "Polynomial Background", mass, RooArgList(c0));

    // === FINAL MODEL: (Signal + Background) ===
    RooRealVar frac("frac", "Fraction of Signal", 0.7, 0, 1.0);
    RooAddPdf model("model", "Signal + Background", RooArgList(signal, background), RooArgList(frac));

    // === FIT THE MODEL ===
    // model.fitTo(data, RooFit::Minimizer("Minuit2", "MIGRAD"));
    model.fitTo(data, RooFit::Save(), RooFit::Minimizer("Minuit", "migradimproved"));

    // === PLOTTING ===
    RooPlot *frame = mass.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    // model.plotOn(frame, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));       // Signal (CB + Gaussian)
    // model.plotOn(frame, RooFit::Components("background"), RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));  // Background

    // Set custom titles for the X and Y axes
    //frame->GetXaxis()->SetTitle("M(#mu^{+}#mu^{-})[GeV/c^{2}]");
    //frame->GetYaxis()->SetTitle("Entries / 0.008 GeV/c^{2}] ");

    // Set a custom title for the RooPlot
    //frame->SetTitle("Fill 9043 J/#psi Invariant Mass Distrbution: LS 90-140");

    //Remove display of stat box (set to 1 for checking)
    gStyle->SetOptStat(0);

    // Error bar color and marker style (20 is filled circles)
    hist->SetMarkerColor(kBlack);
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(20);

    //**********Create canvas, add visual elements, and draw***************/

    // Create canvas to draw histogram on
    TCanvas *canvas1 = new TCanvas("canvas", "Run 304144", 800, 700);
    hist->Draw();
    frame->Draw("same");

    // Create legend
    TLegend *legend = new TLegend(0.13, 0.77, 0.35, 0.87);
    legend->AddEntry(hist, "J/#psi #rightarrow #mu^{+}#mu^{-}", "F");
    legend->AddEntry(frame->getObject(0), "Data", "LP");
    legend->SetBorderSize(0);
    legend->Draw();

    // Use LateX to draw on additional elements
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextFont(42);
    latex.DrawLatex(0.14, 0.735, "p_{T}^{#mu_{1}} > 7 GeV/c");
    latex.DrawLatex(0.14, 0.675, "p_{T}^{#mu_{2}} > 7 GeV/c");
    latex.DrawLatex(0.14, 0.615, "VtxProb > 10%");
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.80, 0.91, "(#sqrt{13} TeV)");

    // Get fit results
    //NOTE: hist->GetVal(); gets raw data value, while value.GetVal(); gets fit value!!
    double mean_val = mean.getVal();
    double mean_error = mean.getError();
    double sigma_val = sigma.getVal();
    double sigma_err = sigma.getError();
    double alpha_val = alpha.getVal();
    double alpha_err = alpha.getError();
    double n_val = n.getVal();
    double n_err = n.getError();
    double frac_val = frac.getVal();
    double frac_err = frac.getError();
    double chi2_val = frame->chiSquare();
    int nEntries = hist->GetEntries();

    //Error Propagation
    double N_signal = frac_val * nEntries;
    double sigma_N_signal = nEntries * sqrt(pow(frac_err, 2) + pow((sqrt(nEntries) / nEntries) * frac_val, 2));  //wrong

    // Create statistics box
    TPaveText *statBox = new TPaveText(0.67, 0.67, 0.87, 0.87, "NDC");
    statBox->SetFillColor(0);
    statBox->SetTextAlign(12);
    statBox->SetBorderSize(1);
    statBox->SetTextSize(0.02);

    // Add text with values
    TText *text;
    text = statBox->AddText(Form("Entries = %d", nEntries));
    text->SetTextSize(0.025);
    //text = statBox->AddText(Form("Mean = %.3f #pm %.3f", mean_val, mean_error));
    text = statBox->AddText(Form("Mean = %.3f #pm %.3f", mean_val, mean_error));
    text->SetTextSize(0.025);
    text = statBox->AddText(Form("#sigma = %.3f #pm %.3f", sigma_val, sigma_err));
    text->SetTextSize(0.025);
    text = statBox->AddText(Form("#alpha = %.3f #pm %.3f", alpha_val, alpha_err));
    text->SetTextSize(0.025);
    text = statBox->AddText(Form("n = %.3f #pm %.3f", n_val, n_err));
    text->SetTextSize(0.025);
    text = statBox->AddText(Form("Frac(s) = %.3f #pm %.3f ", frac_val, frac_err));
    text->SetTextSize(0.025);
    text = statBox->AddText(Form("#chi^{2} = %.3f", chi2_val));
    text->SetTextSize(0.025);

    // Draw on plot
    statBox->Draw();

    // Set ticks on all sides of histogram
    gPad->Update();
    gPad->SetTicks(1, 1);

    // Update canvas with new visuals
    canvas1->Update();

    // Save plot
    canvas1->SaveAs(Form("Fill9043LS%d-%dpt3pt3vtx10pctcrystalball.png", lumiMin, lumiMax));

    //delete canvas1;
    //delete frame;
    //delete statBox;
  }

  // Close file
  file->Close();
  delete file;
}
