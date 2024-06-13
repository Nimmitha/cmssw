#define AnalysisWithIso_cxx
#include "AnalysisWithIso.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnalysisWithIso::Loop() {
  //   In a ROOT session, you can do:
  //      root> .L AnalysisWithIso.C
  //      root> AnalysisWithIso t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  // by  b_branchname->GetEntry(ientry); //read only this branch

  // Ups2 -> z & Ups1 -> Ups
  if (fChain == 0)
    return;

  ofstream myfile;
  ofstream myfile1;
  myfile.precision(17);
  myfile.open("ZmmSelected.dat");
  myfile1.open("ZmmSelected_candidatesEventNumberOfMultiples.dat");
  TFile *fFile = new TFile("ZmmSelected.root", "recreate");
  TTree *fTree = new TTree("ntuple", "ntuple");

  // Now create the branches on tree
  // Define Tree name

  vector<float> *B_TriggerDelta_pt_b;
  vector<float> *Events_b;
  vector<float> *B_Ups1_mass_b;
  vector<float> *B_Ups2_mass_b;
  vector<float> *B_Ups1_VtxProb_b;
  vector<float> *B_Ups2_VtxProb_b;
  vector<float> *B_Ups1To2dY_b;
  vector<float> *B_Ups1_Pt_b;
  vector<float> *B_Ups2_Pt_b;
  vector<float> *B_Ups1_Phi_b;
  vector<float> *B_Ups2_Phi_b;
  vector<float> *B_Mu1_pt_b;
  vector<float> *B_Mu2_pt_b;
  vector<float> *B_Mu3_pt_b;
  vector<float> *B_Mu4_pt_b;
  vector<float> *B_Mu1_eta_b;
  vector<float> *B_Mu2_eta_b;
  vector<float> *B_Mu3_eta_b;
  vector<float> *B_Mu4_eta_b;
  vector<float> *B_Mu1_phi_b;
  vector<float> *B_Mu2_phi_b;
  vector<float> *B_Mu3_phi_b;
  vector<float> *B_Mu4_phi_b;

  vector<float> *mu1mC2_b;
  vector<float> *mu1pC2_b;
  vector<float> *mu2mC2_b;
  vector<float> *mu2pC2_b;
  vector<int> *mu1mNHits_b;
  vector<int> *mu1mNPHits_b;
  vector<int> *mu1pNHits_b;
  vector<int> *mu1pNPHits_b;
  vector<int> *mu2mNHits_b;
  vector<int> *mu2mNPHits_b;
  vector<int> *mu2pNHits_b;
  vector<int> *mu2pNPHits_b;

  vector<bool> *B_Mu1_soft_b;
  vector<bool> *B_Mu1_loose_b;
  vector<bool> *B_Mu1_tight_b;
  vector<bool> *B_Mu2_soft_b;
  vector<bool> *B_Mu2_loose_b;
  vector<bool> *B_Mu2_tight_b;
  vector<bool> *B_Mu3_soft_b;
  vector<bool> *B_Mu3_loose_b;
  vector<bool> *B_Mu3_tight_b;
  vector<bool> *B_Mu4_soft_b;
  vector<bool> *B_Mu4_loose_b;
  vector<bool> *B_Mu4_tight_b;

  vector<float> *Z_xyP1_b;
  vector<float> *Z_xyM1_b;
  vector<float> *Z_zP1_b;
  vector<float> *Z_zM1_b;
  vector<float> *Z_xyP2_b;
  vector<float> *Z_xyM2_b;
  vector<float> *Z_zP2_b;
  vector<float> *Z_zM2_b;

  vector<float> *Y_mass_b;
  vector<float> *Z_mass_b;
  vector<float> *B_J1_mass_b;
  vector<float> *B_J2_mass_b;
  vector<float> *B_J3_mass_b;
  vector<float> *B_J4_mass_b;
  vector<float> *FourL_mass_b;
  vector<float> *FourL_VtxProb_b;
  vector<float> *FourL_pt_b;
  vector<float> *FourL_eta_b;
  vector<float> *FourL_phi_b;
  vector<float> *FourL_rapidity_b;
  vector<float> *Y_pt1_b;
  vector<float> *Y_eta1_b;
  vector<float> *Y_phi1_b;
  vector<float> *Y_pt2_b;
  vector<float> *Y_eta2_b;
  vector<float> *Y_phi2_b;
  vector<float> *Z_pt1_b;
  vector<float> *Z_eta1_b;
  vector<float> *Z_phi1_b;
  vector<float> *Z_pt2_b;
  vector<float> *Z_eta2_b;
  vector<float> *Z_phi2_b;

  vector<bool> *Y_mvaIsoWP90_1_b;
  vector<bool> *Y_mvaIsoWP90_2_b;

  B_TriggerDelta_pt_b = 0;
  Events_b = 0;
  B_Ups1_mass_b = 0;
  B_Ups2_mass_b = 0;
  B_Ups1_Pt_b = 0;
  B_Ups2_Pt_b = 0;
  B_Ups1_Phi_b = 0;
  B_Ups2_Phi_b = 0;
  B_Ups1_VtxProb_b = 0;
  B_Ups2_VtxProb_b = 0;
  B_Ups1To2dY_b = 0;
  Y_mass_b = 0;
  Z_mass_b = 0;
  FourL_mass_b = 0;
  FourL_VtxProb_b = 0;
  FourL_pt_b = 0;
  FourL_eta_b = 0;
  FourL_phi_b = 0;
  FourL_rapidity_b = 0;
  B_Mu1_pt_b = 0;
  B_Mu2_pt_b = 0;
  B_Mu3_pt_b = 0;
  B_Mu4_pt_b = 0;
  B_Mu1_eta_b = 0;
  B_Mu2_eta_b = 0;
  B_Mu3_eta_b = 0;
  B_Mu4_eta_b = 0;
  B_Mu1_phi_b = 0;
  B_Mu2_phi_b = 0;
  B_Mu3_phi_b = 0;
  B_Mu4_phi_b = 0;
  mu1mC2_b = 0;
  mu1pC2_b = 0;
  mu1mNHits_b = 0;
  mu1pNHits_b = 0;
  mu1mNPHits_b = 0;
  mu1pNPHits_b = 0;
  mu2mC2_b = 0;
  mu2pC2_b = 0;
  mu2mNHits_b = 0;
  mu2pNHits_b = 0;
  mu2mNPHits_b = 0;
  mu2pNPHits_b = 0;
  B_Mu1_soft_b = 0;
  B_Mu1_loose_b = 0;
  B_Mu1_tight_b = 0;
  B_Mu2_soft_b = 0;
  B_Mu2_loose_b = 0;
  B_Mu2_tight_b = 0;
  B_Mu3_soft_b = 0;
  B_Mu3_loose_b = 0;
  B_Mu3_tight_b = 0;
  B_Mu4_soft_b = 0;
  B_Mu4_loose_b = 0;
  B_Mu4_tight_b = 0;
  Z_xyP1_b = 0;
  Z_xyM1_b = 0;
  Z_zP1_b = 0;
  Z_zM1_b = 0;
  Z_xyP2_b = 0;
  Z_xyM2_b = 0;
  Z_zP2_b = 0;
  Z_zM2_b = 0;
  B_J1_mass_b = 0;
  B_J2_mass_b = 0;
  B_J3_mass_b = 0;
  B_J4_mass_b = 0;
  Y_pt1_b = 0;
  Y_eta1_b = 0;
  Y_phi1_b = 0;
  Y_pt2_b = 0;
  Y_eta2_b = 0;
  Y_phi2_b = 0;
  Z_pt1_b = 0;
  Z_eta1_b = 0;
  Z_phi1_b = 0;
  Z_pt2_b = 0;
  Z_eta2_b = 0;
  Z_phi2_b = 0;

  b_Y_mvaIsoWP90_1 = 0;
  b_Y_mvaIsoWP90_2 = 0;

  // initialize
  fTree->Branch("B_TriggerDelta_pt", &B_TriggerDelta_pt_b);

  fTree->Branch("Event", &Events_b);
  fTree->Branch("B_Ups1_mass", &B_Ups1_mass_b);
  fTree->Branch("B_Ups2_mass", &B_Ups2_mass_b);
  fTree->Branch("B_Ups1_VtxProb", &B_Ups1_VtxProb_b);
  fTree->Branch("B_Ups2_VtxProb", &B_Ups2_VtxProb_b);
  fTree->Branch("B_Ups1To2dY", &B_Ups1To2dY_b);
  fTree->Branch("B_Ups1_Pt", &B_Ups1_Pt_b);
  fTree->Branch("B_Ups2_Pt", &B_Ups2_Pt_b);
  fTree->Branch("B_Ups1_Phi", &B_Ups1_Phi_b);
  fTree->Branch("B_Ups2_Phi", &B_Ups2_Phi_b);
  fTree->Branch("Y_mass", &Y_mass_b);
  fTree->Branch("Z_mass", &Z_mass_b);
  fTree->Branch("B_J3_mass", &B_J3_mass_b);
  fTree->Branch("B_J4_mass", &B_J4_mass_b);
  fTree->Branch("Y_pt1", &Y_pt1_b);
  fTree->Branch("Y_pt2", &Y_pt2_b);
  fTree->Branch("Y_eta1", &Y_eta1_b);
  fTree->Branch("Y_eta2", &Y_eta2_b);
  fTree->Branch("Y_phi1", &Y_phi1_b);
  fTree->Branch("Y_phi2", &Y_phi2_b);

  fTree->Branch("Z_pt1", &Z_pt1_b);
  fTree->Branch("Z_pt2", &Z_pt2_b);
  fTree->Branch("Z_eta1", &Z_eta1_b);
  fTree->Branch("Z_eta2", &Z_eta2_b);
  fTree->Branch("Z_phi1", &Z_phi1_b);
  fTree->Branch("Z_phi2", &Z_phi2_b);
  fTree->Branch("FourL_VtxProb", &FourL_VtxProb_b);
  fTree->Branch("FourL_mass", &FourL_mass_b);
  fTree->Branch("FourL_pt", &FourL_pt_b);
  fTree->Branch("FourL_eta", &FourL_eta_b);
  fTree->Branch("FourL_phi", &FourL_phi_b);
  fTree->Branch("FourL_rapidity", &FourL_rapidity_b);
  fTree->Branch("B_Mu1_pt", &B_Mu1_pt_b);
  fTree->Branch("B_Mu2_pt", &B_Mu2_pt_b);
  fTree->Branch("B_Mu3_pt", &B_Mu3_pt_b);
  fTree->Branch("B_Mu4_pt", &B_Mu4_pt_b);
  fTree->Branch("B_Mu1_eta", &B_Mu1_eta_b);
  fTree->Branch("B_Mu2_eta", &B_Mu2_eta_b);
  fTree->Branch("B_Mu3_eta", &B_Mu3_eta_b);
  fTree->Branch("B_Mu4_eta", &B_Mu4_eta_b);
  fTree->Branch("B_Mu1_phi", &B_Mu1_phi_b);
  fTree->Branch("B_Mu2_phi", &B_Mu2_phi_b);
  fTree->Branch("B_Mu3_phi", &B_Mu3_phi_b);
  fTree->Branch("B_Mu4_phi", &B_Mu4_phi_b);

  fTree->Branch("mu1mC2", &mu1mC2_b);
  fTree->Branch("mu1pC2", &mu1pC2_b);
  fTree->Branch("mu1mNHits", &mu1mNHits_b);
  fTree->Branch("mu1pNHits", &mu1pNHits_b);
  fTree->Branch("mu1mNPHits", &mu1mNPHits_b);
  fTree->Branch("mu1pNPHits", &mu1pNPHits_b);
  fTree->Branch("mu2mC2", &mu2mC2_b);
  fTree->Branch("mu2pC2", &mu2pC2_b);
  fTree->Branch("mu2mNHits", &mu2mNHits_b);
  fTree->Branch("mu2pNHits", &mu2pNHits_b);
  fTree->Branch("mu2mNPHits", &mu2mNPHits_b);
  fTree->Branch("mu2pNPHits", &mu2pNPHits_b);
  fTree->Branch("B_Mu1_soft", &B_Mu1_soft_b);
  fTree->Branch("B_Mu1_loose", &B_Mu1_loose_b);
  fTree->Branch("B_Mu1_tight", &B_Mu1_tight_b);
  fTree->Branch("B_Mu2_soft", &B_Mu2_soft_b);
  fTree->Branch("B_Mu2_loose", &B_Mu2_loose_b);
  fTree->Branch("B_Mu2_tight", &B_Mu2_tight_b);
  fTree->Branch("B_Mu3_soft", &B_Mu3_soft_b);
  fTree->Branch("B_Mu3_loose", &B_Mu3_loose_b);
  fTree->Branch("B_Mu3_tight", &B_Mu3_tight_b);
  fTree->Branch("B_Mu4_soft", &B_Mu4_soft_b);
  fTree->Branch("B_Mu4_loose", &B_Mu4_loose_b);
  fTree->Branch("B_Mu4_tight", &B_Mu4_tight_b);

  fTree->Branch("Z_xyP1", &Z_xyP1_b);
  fTree->Branch("Z_xyM1", &Z_xyM1_b);
  fTree->Branch("Z_zP1", &Z_zP1_b);
  fTree->Branch("Z_zM1", &Z_zM1_b);
  fTree->Branch("Z_xyP2", &Z_xyP2_b);
  fTree->Branch("Z_xyM2", &Z_xyM2_b);
  fTree->Branch("Z_zP2", &Z_zP2_b);
  fTree->Branch("Z_zM2", &Z_zM2_b);

  fTree->Branch("Z_zP2", &Z_zP2_b);
  fTree->Branch("Z_zM2", &Z_zM2_b);

  fTree->Branch("Y_mvaIsoWP90_1", &b_Y_mvaIsoWP90_1);
  fTree->Branch("Y_mvaIsoWP90_2", &b_Y_mvaIsoWP90_2);
  fTree->Branch("Z_soft1", &b_Z_soft1);
  fTree->Branch("Z_soft2", &b_Z_soft2);

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t ncandiPreSelection(0), nevPreSelection(0), ncandiTrigger(0), nevTrigger(0), ncandiSoft(0), nevSoft(0), ncandiDetector(0), nevDetector(0), ncandiLite4Vertex(0), nevLite4Vertex(0),
      ncandiCombine(0), nevCombine(0), ncandiIsolation1(0), nevIsolation1(0), ncandiIsolation2(0), nevIsolation2(0), ncandiIsolation3(0), nevIsolation3(0), ncandiIsolation4(0), nevIsolation4(0),
      ncandiMuonPtGreaterThan(0), nevMuonPtGreaterThan(0), ncandiYZVertexing(0), nevYZVertexing(0), ncandiYZPt(0), nevYZPt(0), ncandiYZMass(0), nevYZMass(0), ncandi4Vertex(0), nev4Vertex(0),
      ncandi4Pt(0), nev4Pt(0), ncandi4Mass(0), nev4Mass(0), ncandiYmvaiso(0), nevYmvaiso(0), ncandiYmass(0), nevYmass(0), ncandiZmass(0), nevZmass(0);
  Long64_t temp_eventPreSelection(0), temp_eventTrigger(0), temp_eventSoft(0), temp_eventDetector(0), temp_eventLite4Vertex(0), temp_eventCombine(0), temp_eventIsolation1(0), temp_eventIsolation2(0),
      temp_eventIsolation3(0), temp_eventIsolation4(0), temp_eventMuonPtGreaterThan(0), temp_eventYZVertexing(0), temp_eventYZPt(0), temp_eventYZMass(0), temp_event4Vertex(0), temp_event4Pt(0),
      temp_event4Mass(0), temp_eventYmvaiso(0), temp_eventYmass(0), temp_eventZmass(0);
  Long64_t nbytes = 0, nb = 0;
  Double_t Events = 0;
  int doublecount = 0;
  TLorentzVector FourMu, M1, M2, M3, M4;
  float FourL_rapidity;

  cout << "Total Entries: " << nentries << endl;

  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    // Multiple Candidate
    // event no. of those having multiple candidate

    for (int i = 0; i < Event->size(); i++) {
      // if (Cut(ientry) < 0) continue;

      // preselection cuts
      ncandiPreSelection++;
      Events = Event->at(i);

      if (temp_eventPreSelection != Event->at(i)) {
        temp_eventPreSelection = Event->at(i);
        nevPreSelection++;
      }

      // cut 2 Trigger
      if (Y_TriggerPath->at(i) == 0)
        continue;

      ncandiTrigger++;
      if (temp_eventTrigger != Event->at(i)) {
        temp_eventTrigger = Event->at(i);
        nevTrigger++;
      }

      if (Z_soft1->at(i) < 1 || Z_soft2->at(i) < 1)
        continue;

      ncandiSoft++;
      if (temp_eventSoft != Event->at(i)) {
        temp_eventSoft = Event->at(i);
        nevSoft++;
      }

      if (Y_mvaIsoWP90_1->at(i) < 1 || Y_mvaIsoWP90_2->at(i) < 1)
        continue;

      ncandiYmvaiso++;
      if (temp_eventYmvaiso != Event->at(i)) {
        temp_eventYmvaiso = Event->at(i);
        nevYmvaiso++;
      }

      if (Z_pt1->at(i) < 3.0 || Z_pt2->at(i) < 3.0 || Y_pt1->at(i) < 5.0 || Y_pt2->at(i) < 5.0)
        continue;

      if (abs(Z_eta1->at(i)) > 2.4 || abs(Z_eta2->at(i)) > 2.4 || abs(Y_eta1->at(i)) > 2.5 || abs(Y_eta2->at(i)) > 2.5)
        continue;

      ncandiDetector++;
      if (temp_eventDetector != Event->at(i)) {
        temp_eventDetector = Event->at(i);
        nevDetector++;
      }

      if (Y_VtxProb->at(i) < 0.01 || Z_VtxProb->at(i) < 0.01)
        continue;

      ncandiYZVertexing++;
      if (temp_eventYZVertexing != Event->at(i)) {
        temp_eventYZVertexing = Event->at(i);
        nevYZVertexing++;
      }

      if (FourL_VtxProb->at(i) < 0.01)
        continue;

      ncandi4Vertex++;
      if (temp_event4Vertex != Event->at(i)) {
        temp_event4Vertex = Event->at(i);
        nev4Vertex++;
      }

      // find the Upsilon//
      float Ups_VtxProb1 = 0;
      float Ups_VtxProb2 = 0;
      float Ups_Pt1 = 0;
      float Ups_Pt2 = 0;
      float Ups_Phi1 = 0;
      float Ups_Phi2 = 0;
      float Ups1_mass = 0;
      float Ups2_mass = 0;
      float Ups_Rapidity1 = 0;
      float Ups_Rapidity2 = 0;
      float Ups1To2_dY = 0;

      // Now Find

      Ups1_mass = Y_mass->at(i);
      Ups_VtxProb1 = Y_VtxProb->at(i);
      Ups_Pt1 = Y_pt->at(i);
      Ups_Rapidity1 = Y_rapidity->at(i);
      Ups_Phi1 = Y_phi->at(i);
      Ups2_mass = Z_mass->at(i);
      Ups_VtxProb2 = Z_VtxProb->at(i);
      Ups_Pt2 = Z_pt->at(i);
      Ups_Rapidity2 = Z_rapidity->at(i);
      Ups_Phi2 = Z_phi->at(i);

      // Dilepton mass cut 5
      if (Ups1_mass < 8 || Ups1_mass > 10)
        continue;

      ncandiYmass++;
      if (temp_eventYmass != Event->at(i)) {
        temp_eventYmass = Event->at(i);
        nevYmass++;
      }

      if (Ups2_mass < 80.0 || Ups2_mass > 110)
        continue;

      ncandiZmass++;
      if (temp_eventZmass != Event->at(i)) {
        temp_eventZmass = Event->at(i);
        nevZmass++;
      }

      // FourL_mass Cut
      if (FourL_mass->at(i) < 112 || FourL_mass->at(i) > 162)
        continue;

      ncandi4Mass++;
      // Events=Event->at(i);
      if (temp_event4Mass != Event->at(i)) {
        temp_event4Mass = Event->at(i);
        nev4Mass++;

        // continue;
      } else {
        continue;
        // doublecount++;
        // myfile<<doublecount<<" double Event number detected # "<<Event->at(i)<<endl;
        // myfile<<" Event number of previous"<<Event->at(i-1)<<endl;
        // myfile<<Event->at(i)<<endl;

        cout << "Event no for multiple Candidate " << Event->at(i) << endl;
        myfile1 << Event->at(i) << endl;
      }

      myfile << Event->at(i) << " " << Ups1_mass << " " << Ups2_mass << " " << FourL_mass->at(i) << " " << FourL_VtxProb->at(i) << endl;

      // myfile<<Ups1_mass<<endl;

      // if (ncandi>1) continue;
      Events_b->push_back(Events);
      B_Ups1_mass_b->push_back(Ups1_mass);
      B_Ups2_mass_b->push_back(Ups2_mass);
      B_Ups1_VtxProb_b->push_back(Ups_VtxProb1);
      B_Ups2_VtxProb_b->push_back(Ups_VtxProb2);
      B_Ups1To2dY_b->push_back(Ups1To2_dY);
      B_Ups1_Pt_b->push_back(Ups_Pt1);
      B_Ups1_Phi_b->push_back(Ups_Phi1);
      B_Ups2_Pt_b->push_back(Ups_Pt2);
      B_Ups2_Phi_b->push_back(Ups_Phi2);
      FourL_VtxProb_b->push_back(FourL_VtxProb->at(i));
      Y_mass_b->push_back(Y_mass->at(i));
      Z_mass_b->push_back(Z_mass->at(i));
      B_TriggerDelta_pt_b->push_back(Y_TriggerPt1->at(i) - Y_pt1->at(i));

      FourL_mass_b->push_back(FourL_mass->at(i));
      FourL_pt_b->push_back(FourL_pt->at(i));
      FourL_eta_b->push_back(FourL_eta->at(i));
      FourL_phi_b->push_back(FourL_phi->at(i));
      FourL_rapidity_b->push_back(FourL_rapidity);

      Y_pt1_b->push_back(Y_pt1->at(i));
      Y_pt2_b->push_back(Y_pt2->at(i));
      Y_eta1_b->push_back(Y_eta1->at(i));
      Y_eta2_b->push_back(Y_eta2->at(i));

      Z_pt1_b->push_back(Z_pt1->at(i));
      Z_pt2_b->push_back(Z_pt2->at(i));
      Z_eta1_b->push_back(Z_eta1->at(i));
      Z_eta2_b->push_back(Z_eta2->at(i));
    }
  }
  cout << "Number of event Preselection=" << nevPreSelection << endl;
  cout << "Number of Candidate Preselection=" << ncandiPreSelection << endl;

  cout << "Number of event Trigger=" << nevTrigger << endl;
  cout << "Number of Candidate Trigger=" << ncandiTrigger << endl;

  cout << "Number of event Soft=" << nevSoft << endl;
  cout << "Number of Candidate Soft=" << ncandiSoft << endl;

  cout << "Number of event electron ID=" << nevYmvaiso << endl;
  cout << "Number of Candidate electron ID=" << ncandiYmvaiso << endl;

  cout << "Number of event Detector=" << nevDetector << endl;
  cout << "Number of Candidate Detector=" << ncandiDetector << endl;

  cout << "Number of event Dilepton Vertexing=" << nevYZVertexing << endl;
  cout << "Number of Candidate Dilepton Vertexing=" << ncandiYZVertexing << endl;

  cout << "Number of events 4Vertex=" << nev4Vertex << endl;
  cout << "Number of Candidate 4Vertex=" << ncandi4Vertex << endl;

  cout << "Number of event Y mass=" << nevYmass << endl;
  cout << "Number of Candidate Y mass=" << ncandiYmass << endl;

  cout << "Number of event Z mass=" << nevZmass << endl;
  cout << "Number of Candidate Z mass=" << ncandiZmass << endl;

  cout << "Number of event 4Mass=" << nev4Mass << endl;
  cout << "Number of Candidate 4Mass=" << ncandi4Mass << endl;

  if (ncandiPreSelection > 0 || ncandiTrigger > 0 || ncandiSoft > 0 || ncandiDetector > 0 || ncandiLite4Vertex > 0 || ncandiCombine > 0 || ncandiMuonPtGreaterThan > 0 || ncandiYZVertexing > 0 ||
      ncandiYZPt > 0 || ncandiYZMass > 0 || ncandi4Vertex > 0 || ncandi4Mass > 0) {
    fTree->Fill();
  }

  fFile->Write();
  myfile.close();
  myfile1.close();
  fFile->Close();

  delete fFile;
}