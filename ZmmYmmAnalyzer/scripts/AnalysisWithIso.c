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
  //by  b_branchname->GetEntry(ientry); //read only this branch

  // Ups2 -> z & Ups1 -> Ups
  if (fChain == 0)
    return;
  ofstream myfile;
  ofstream myfile1;
  myfile.precision(17);
  //myfile.open("SingleMuon_Run_2016_2017take2_2018_most_HimalsCuts_candidatesE.dat");
  //myfile1.open("SingleMuon_Run_2016_2017take2_2018_most_HimalsCuts_candidatesEventNumberOfMultiples.dat");
  //TFile *fFile = new TFile("SingleMuon_Run_2016_2017take2_2018_most_HimalsCuts_candidatesE.root","recreate");
  //myfile.open("ZuuYuu_singleMuonRun2_HimalsCutsBlindedRegion_candidatesE.dat");
  //myfile1.open("ZuuYuu_singleMuonRun2_HimalsCutsBlindedRegion_candidatesEventNumberOfMultiples.dat");
  //TFile *fFile = new TFile("ZuuYuu_singleMuonRun2_BlindedRegion_candidatesE.root","recreate");
  myfile.open("temp.dat");
  myfile1.open("tempNumberOfMultiples.dat");
  TFile *fFile = new TFile("temp.root", "recreate");
  TTree *fTree = new TTree("ntuple", "ntuple");
  //Now create the branches on tree
  //Define Tree name
  vector<float> *Events_b;
  vector<float> *B_Ups1_mass_b;
  vector<float> *B_Ups2_mass_b;
  vector<float> *B_Ups1_VtxProb_b;
  vector<float> *B_Ups2_VtxProb_b;
  vector<float> *B_Ups1To2dY_b;
  vector<float> *B_Ups1_Pt_b;
  vector<float> *B_Ups2_Pt_b;
  vector<float> *B_Ups1_Eta_b;
  vector<float> *B_Ups2_Eta_b;
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

  vector<float> *B_J_xyP1_b;
  vector<float> *B_J_xyM1_b;
  vector<float> *B_J_zP1_b;
  vector<float> *B_J_zM1_b;
  vector<float> *B_J_xyP2_b;
  vector<float> *B_J_xyM2_b;
  vector<float> *B_J_zP2_b;
  vector<float> *B_J_zM2_b;

  vector<float> *B_Z_mass_b;
  vector<float> *B_J_mass_b;
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
  vector<float> *B_Z_pt1_b;
  vector<float> *B_Z_eta1_b;
  vector<float> *B_Z_phi1_b;
  vector<float> *B_Z_pt2_b;
  vector<float> *B_Z_eta2_b;
  vector<float> *B_Z_phi2_b;
  vector<float> *B_J_pt1_b;
  vector<float> *B_J_eta1_b;
  vector<float> *B_J_phi1_b;
  vector<float> *B_J_pt2_b;
  vector<float> *B_J_eta2_b;
  vector<float> *B_J_phi2_b;

  Events_b = 0;
  B_Ups1_mass_b = 0;
  B_Ups2_mass_b = 0;
  B_Ups1_Pt_b = 0;
  B_Ups2_Pt_b = 0;
  B_Ups1_Eta_b = 0;
  B_Ups2_Eta_b = 0;
  B_Ups1_Phi_b = 0;
  B_Ups2_Phi_b = 0;
  B_Ups1_VtxProb_b = 0;
  B_Ups2_VtxProb_b = 0;
  B_Ups1To2dY_b = 0;
  B_Z_mass_b = 0;
  B_J_mass_b = 0;
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
  B_J_xyP1_b = 0;
  B_J_xyM1_b = 0;
  B_J_zP1_b = 0;
  B_J_zM1_b = 0;
  B_J_xyP2_b = 0;
  B_J_xyM2_b = 0;
  B_J_zP2_b = 0;
  B_J_zM2_b = 0;
  B_J1_mass_b = 0;
  B_J2_mass_b = 0;
  B_J3_mass_b = 0;
  B_J4_mass_b = 0;
  B_Z_pt1_b = 0;
  B_Z_eta1_b = 0;
  B_Z_phi1_b = 0;
  B_Z_pt2_b = 0;
  B_Z_eta2_b = 0;
  B_Z_phi2_b = 0;
  B_J_pt1_b = 0;
  B_J_eta1_b = 0;
  B_J_phi1_b = 0;
  B_J_pt2_b = 0;
  B_J_eta2_b = 0;
  B_J_phi2_b = 0;

  //initialize

  fTree->Branch("Event", &Events_b);
  fTree->Branch("B_Ups1_mass", &B_Ups1_mass_b);
  fTree->Branch("B_Ups2_mass", &B_Ups2_mass_b);
  fTree->Branch("B_Ups1_VtxProb", &B_Ups1_VtxProb_b);
  fTree->Branch("B_Ups2_VtxProb", &B_Ups2_VtxProb_b);
  fTree->Branch("B_Ups1To2dY", &B_Ups1To2dY_b);
  fTree->Branch("B_Ups1_Pt", &B_Ups1_Pt_b);
  fTree->Branch("B_Ups2_Pt", &B_Ups2_Pt_b);
  fTree->Branch("B_Ups1_Eta", &B_Ups1_Eta_b);
  fTree->Branch("B_Ups2_Eta", &B_Ups2_Eta_b);
  fTree->Branch("B_Ups1_Phi", &B_Ups1_Phi_b);
  fTree->Branch("B_Ups2_Phi", &B_Ups2_Phi_b);
  fTree->Branch("B_J1_mass", &B_J1_mass_b);
  fTree->Branch("B_J2_mass", &B_J2_mass_b);
  fTree->Branch("B_J3_mass", &B_J3_mass_b);
  fTree->Branch("B_J4_mass", &B_J4_mass_b);
  //fTree->Branch("B_Z_pt1",             &B_Z_pt1_b);
  //fTree->Branch("B_Z_pt2",             &B_Z_pt2_b);
  //fTree->Branch("B_Z_eta1",             &B_Z_eta1_b);
  //fTree->Branch("B_Z_eta2",             &B_Z_eta2_b);
  //fTree->Branch("B_Z_phi1",             &B_Z_phi1_b);
  //fTree->Branch("B_Z_phi2",             &B_Z_phi2_b);

  //fTree->Branch("B_J_pt1",             &B_J_pt1_b);
  //fTree->Branch("B_J_pt2",             &B_J_pt2_b);
  //fTree->Branch("B_J_eta1",             &B_J_eta1_b);
  //fTree->Branch("B_J_eta2",             &B_J_eta2_b);
  //fTree->Branch("B_J_phi1",             &B_J_phi1_b);
  //fTree->Branch("B_J_phi2",             &B_J_phi2_b);
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

  fTree->Branch("B_J_xyP1", &B_J_xyP1_b);
  fTree->Branch("B_J_xyM1", &B_J_xyM1_b);
  fTree->Branch("B_J_zP1", &B_J_zP1_b);
  fTree->Branch("B_J_zM1", &B_J_zM1_b);
  fTree->Branch("B_J_xyP2", &B_J_xyP2_b);
  fTree->Branch("B_J_xyM2", &B_J_xyM2_b);
  fTree->Branch("B_J_zP2", &B_J_zP2_b);
  fTree->Branch("B_J_zM2", &B_J_zM2_b);

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t ncandiPreSelection(0), nevPreSelection(0), ncandiTrigger(0), nevTrigger(0), ncandiSoft(0), nevSoft(0), ncandiDetector(0), nevDetector(0),
      ncandiLite4Vertex(0), nevLite4Vertex(0), ncandiCombine(0), nevCombine(0), ncandiIsolation1(0), nevIsolation1(0), ncandiIsolation2(0),
      nevIsolation2(0), ncandiIsolation3(0), nevIsolation3(0), ncandiIsolation4(0), nevIsolation4(0), ncandiMuonPtGreaterThan(0),
      nevMuonPtGreaterThan(0), ncandiYZVertexing(0), nevYZVertexing(0), ncandiYZPt(0), nevYZPt(0), ncandiYZMass(0), nevYZMass(0), ncandi4Vertex(0),
      nev4Vertex(0), ncandi4Pt(0), nev4Pt(0), ncandi4Mass(0), nev4Mass(0);
  Long64_t temp_eventPreSelection(0), temp_eventTrigger(0), temp_eventSoft(0), temp_eventDetector(0), temp_eventLite4Vertex(0), temp_eventCombine(0),
      temp_eventIsolation1(0), temp_eventIsolation2(0), temp_eventIsolation3(0), temp_eventIsolation4(0), temp_eventMuonPtGreaterThan(0),
      temp_eventYZVertexing(0), temp_eventYZPt(0), temp_eventYZMass(0), temp_event4Vertex(0), temp_event4Pt(0), temp_event4Mass(0);
  Long64_t nbytes = 0, nb = 0;
  Double_t Events = 0;
  int doublecount = 0;
  TLorentzVector FourMu, M1, M2, M3, M4;
  float FourL_rapidity;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    //Multiple Candidate
    //event no. of those havinf multiple candidate
    //first 39 from 2018, then 19 from 2017, then 22 from 2016

    for (int i = 0; i < Event->size(); i++) {
      // if (Cut(ientry) < 0) continue;

      //preselection cuts

      ncandiPreSelection++;
      Events = Event->at(i);
      if (temp_eventPreSelection != Event->at(i)) {
        temp_eventPreSelection = Event->at(i);
        nevPreSelection++;
      }

      // cut 2 Trigger
      if (B_U_TriggerPath->at(i) == 0)
        continue;

      ncandiTrigger++;
      //Events=Event->at(i);
      if (temp_eventTrigger != Event->at(i)) {
        temp_eventTrigger = Event->at(i);
        nevTrigger++;
      }

      if (B_Mu1_soft->at(i) < 1 || B_Mu2_soft->at(i) < 1 || B_Mu3_soft->at(i) < 1 || B_Mu4_soft->at(i) < 1)
        continue;

      ncandiSoft++;
      //Events=Event->at(i);
      if (temp_eventSoft != Event->at(i)) {
        temp_eventSoft = Event->at(i);
        nevSoft++;
      }
      //if ( ( B_Mu1_tight->at(i) + B_Mu2_tight->at(i) + B_Mu3_tight->at(i) + B_Mu4_tight->at(i) ) < 3) continue;
      if (B_Mu1_pt->at(i) < 3.0 || B_Mu2_pt->at(i) < 3.0 || B_Mu3_pt->at(i) < 3.0 || B_Mu4_pt->at(i) < 3.0)
        continue;
      if (abs(B_Mu1_eta->at(i)) > 2.4 || abs(B_Mu2_eta->at(i)) > 2.4 || abs(B_Mu3_eta->at(i)) > 2.4 || abs(B_Mu4_eta->at(i)) > 2.4)
        continue;

      ncandiDetector++;
      //Events=Event->at(i);
      if (temp_eventDetector != Event->at(i)) {
        temp_eventDetector = Event->at(i);
        nevDetector++;
      }

      int UpsVtx = 0;
      if (B_J1_VtxProb->at(i) > 0.01) {
        if (B_J2_VtxProb->at(i) > 0.01) {
          UpsVtx++;
        }
      }
      if (B_J3_VtxProb->at(i) > 0.01) {
        if (B_J4_VtxProb->at(i) > 0.01) {
          UpsVtx++;
        }
      }
      if (UpsVtx < 1)
        continue;
      if (FourL_VtxProb->at(i) < 0.01)
        continue;

      ncandiLite4Vertex++;
      //Events=Event->at(i);
      if (temp_eventLite4Vertex != Event->at(i)) {
        temp_eventLite4Vertex = Event->at(i);
        nevLite4Vertex++;
      }

      int check1, check2, UpsMass;
      UpsMass = check1 = check2 = 0;
      if ((B_J1_mass->at(i) > 9.0 && B_J1_mass->at(i) < 9.7) || (B_J1_mass->at(i) > 85.0 && B_J1_mass->at(i) < 110.0)) {
        if ((B_J2_mass->at(i) > 9.0 && B_J2_mass->at(i) < 9.7) || (B_J2_mass->at(i) > 85.0 && B_J2_mass->at(i) < 110.0)) {
          UpsMass++;
          check1++;
        }
      }
      if ((B_J3_mass->at(i) > 9.0 && B_J3_mass->at(i) < 9.7) || (B_J3_mass->at(i) > 85.0 && B_J3_mass->at(i) < 110.0)) {
        if ((B_J4_mass->at(i) > 9.0 && B_J4_mass->at(i) < 9.7) || (B_J4_mass->at(i) > 85.0 && B_J4_mass->at(i) < 110.0)) {
          UpsMass++;
          check2++;
        }
      }
      //if (UpsMass<1) continue;

      //find the Upsilon//
      float Ups_VtxProb1 = 0;
      float Ups_VtxProb2 = 0;
      float Ups_Pt1 = 0;
      float Ups_Pt2 = 0;
      float Ups_Phi1 = 0;
      float Ups_Phi2 = 0;
      float Ups1_mass = 0;
      float Ups2_mass = 0;
      float Ups_Eta1 = 0;
      float Ups_Eta2 = 0;
      float Ups_Rapidity1 = 0;
      float Ups_Rapidity2 = 0;
      float Ups1To2_dY = 0;

      //Now Find
      //if we find 4 Upsilon i.e
      
      if (UpsMass == 2) {
        //choose the upsilon having best vertex probability
        //cout<<"Are more than 2 upsilon in a event"<<endl;
        if ((B_J1_VtxProb->at(i) + B_J2_VtxProb->at(i)) > (B_J3_VtxProb->at(i) + B_J4_VtxProb->at(i))) {
          if ((B_J1_mass->at(i) > 9.0 && B_J1_mass->at(i) < 9.7)) {
            Ups1_mass = B_J1_mass->at(i);
            Ups_VtxProb1 = B_J1_VtxProb->at(i);
            Ups_Pt1 = B_J1_pt->at(i);
            Ups_Eta1 = B_J1_eta->at(i);
            Ups_Rapidity1 = B_J1_rapidity->at(i);
            Ups_Phi1 = B_J1_phi->at(i);
            Ups2_mass = B_J2_mass->at(i);
            Ups_VtxProb2 = B_J2_VtxProb->at(i);
            Ups_Pt2 = B_J2_pt->at(i);
            Ups_Eta2 = B_J2_eta->at(i);
            Ups_Rapidity2 = B_J2_rapidity->at(i);
            Ups_Phi2 = B_J2_phi->at(i);
          } else {
            Ups1_mass = B_J2_mass->at(i);
            Ups_VtxProb1 = B_J2_VtxProb->at(i);
            Ups_Pt1 = B_J2_pt->at(i);
            Ups_Eta1 = B_J2_eta->at(i);
            Ups_Rapidity1 = B_J2_rapidity->at(i);
            Ups_Phi1 = B_J2_phi->at(i);
            Ups2_mass = B_J1_mass->at(i);
            Ups_VtxProb2 = B_J1_VtxProb->at(i);
            Ups_Pt2 = B_J1_pt->at(i);
            Ups_Eta2 = B_J1_eta->at(i);
            Ups_Rapidity2 = B_J1_rapidity->at(i);
            Ups_Phi2 = B_J1_phi->at(i);
          }
        } else {
          if ((B_J3_mass->at(i) > 9.0 && B_J3_mass->at(i) < 9.7)) {
            Ups1_mass = B_J3_mass->at(i);
            Ups_VtxProb1 = B_J3_VtxProb->at(i);
            Ups_Pt1 = B_J3_pt->at(i);
            Ups_Eta1 = B_J3_eta->at(i);
            Ups_Rapidity1 = B_J3_rapidity->at(i);
            Ups_Phi1 = B_J3_phi->at(i);
            Ups2_mass = B_J4_mass->at(i);
            Ups_VtxProb2 = B_J4_VtxProb->at(i);
            Ups_Pt2 = B_J4_pt->at(i);
            Ups_Eta2 = B_J4_eta->at(i);
            Ups_Rapidity2 = B_J4_rapidity->at(i);
            Ups_Phi2 = B_J4_phi->at(i);
          } else {
            Ups1_mass = B_J4_mass->at(i);
            Ups_VtxProb1 = B_J4_VtxProb->at(i);
            Ups_Pt1 = B_J4_pt->at(i);
            Ups_Eta1 = B_J4_eta->at(i);
            Ups_Rapidity1 = B_J4_rapidity->at(i);
            Ups_Phi1 = B_J4_phi->at(i);
            Ups2_mass = B_J3_mass->at(i);
            Ups_VtxProb2 = B_J3_VtxProb->at(i);
            Ups_Pt2 = B_J3_pt->at(i);
            Ups_Eta2 = B_J3_eta->at(i);
            Ups_Rapidity2 = B_J3_rapidity->at(i);
            Ups_Phi2 = B_J3_phi->at(i);
          }
        }
      }
      //if have only one pairs of upsislon
      else if (UpsMass == 1) {
        if (check1 > 0) {
          if (B_J1_mass->at(i) > 9.0 && B_J1_mass->at(i) < 9.7) {
            Ups1_mass = B_J1_mass->at(i);
            Ups_VtxProb1 = B_J1_VtxProb->at(i);
            Ups_Pt1 = B_J1_pt->at(i);
            Ups_Eta1 = B_J1_eta->at(i);
            Ups_Rapidity1 = B_J1_rapidity->at(i);
            Ups_Phi1 = B_J1_phi->at(i);
            Ups2_mass = B_J2_mass->at(i);
            Ups_VtxProb2 = B_J2_VtxProb->at(i);
            Ups_Pt2 = B_J2_pt->at(i);
            Ups_Eta2 = B_J2_eta->at(i);
            Ups_Rapidity2 = B_J2_rapidity->at(i);
            Ups_Phi2 = B_J2_phi->at(i);
          } else if (B_J2_mass->at(i) > 9.0 && B_J2_mass->at(i) < 9.7) {
            Ups1_mass = B_J2_mass->at(i);
            Ups_VtxProb1 = B_J2_VtxProb->at(i);
            Ups_Pt1 = B_J2_pt->at(i);
            Ups_Eta1 = B_J2_eta->at(i);
            Ups_Rapidity1 = B_J2_rapidity->at(i);
            Ups_Phi1 = B_J2_phi->at(i);
            Ups2_mass = B_J1_mass->at(i);
            Ups_VtxProb2 = B_J1_VtxProb->at(i);
            Ups_Pt2 = B_J1_pt->at(i);
            Ups_Eta2 = B_J1_eta->at(i);
            Ups_Rapidity2 = B_J1_rapidity->at(i);
            Ups_Phi2 = B_J1_phi->at(i);
          }
        } else if (check2 > 0) {
          if (B_J3_mass->at(i) > 9.0 && B_J3_mass->at(i) < 9.7) {
            Ups1_mass = B_J3_mass->at(i);
            Ups_VtxProb1 = B_J3_VtxProb->at(i);
            Ups_Pt1 = B_J3_pt->at(i);
            Ups_Rapidity1 = B_J3_rapidity->at(i);
            Ups_Phi1 = B_J3_phi->at(i);
            Ups2_mass = B_J4_mass->at(i);
            Ups_VtxProb2 = B_J4_VtxProb->at(i);
            Ups_Pt2 = B_J4_pt->at(i);
            Ups_Rapidity2 = B_J4_rapidity->at(i);
            Ups_Phi2 = B_J4_phi->at(i);
          } else if (B_J4_mass->at(i) > 9.0 && B_J4_mass->at(i) < 9.7) {
            Ups1_mass = B_J4_mass->at(i);
            Ups_VtxProb1 = B_J4_VtxProb->at(i);
            Ups_Pt1 = B_J4_pt->at(i);
            Ups_Rapidity1 = B_J4_rapidity->at(i);
            Ups_Phi1 = B_J4_phi->at(i);
            Ups2_mass = B_J3_mass->at(i);
            Ups_VtxProb2 = B_J3_VtxProb->at(i);
            Ups_Pt2 = B_J3_pt->at(i);
            Ups_Rapidity2 = B_J3_rapidity->at(i);
            Ups_Phi2 = B_J3_phi->at(i);
          }
          //else continue;
        }
      }
      //if have no upsilon
      else {
        continue;
      }

      float Mu_mass = 0.1056583745;
      M1.SetPtEtaPhiM(B_Mu1_pt->at(i), B_Mu1_eta->at(i), B_Mu1_phi->at(i), Mu_mass);
      M2.SetPtEtaPhiM(B_Mu2_pt->at(i), B_Mu2_eta->at(i), B_Mu2_phi->at(i), Mu_mass);
      M3.SetPtEtaPhiM(B_Mu3_pt->at(i), B_Mu3_eta->at(i), B_Mu3_phi->at(i), Mu_mass);
      M4.SetPtEtaPhiM(B_Mu4_pt->at(i), B_Mu4_eta->at(i), B_Mu4_phi->at(i), Mu_mass);
      FourMu = M1 + M2 + M3 + M4;

      ncandiCombine++;
      //Events=Event->at(i);
      if (temp_eventCombine != Event->at(i)) {
        temp_eventCombine = Event->at(i);
        nevCombine++;
      }

      //Calculation of Isolation 1
      //Calculate Delta R between the muons
      float DeltaR12 = M1.DeltaR(M2);
      float DeltaR13 = M1.DeltaR(M3);
      float DeltaR14 = M1.DeltaR(M4);
      float B_Mu_IsoTrackCorr1 = B_Mu1_IsoTrack->at(i);
      if (DeltaR12 < 0.3) {
        //      cout<<" Himal deltaR12 "<< DeltaR12<< "  " << Mu_TrakIso[0] << " - " << P1.Pt() << endl;
        // really subtract?
        B_Mu_IsoTrackCorr1 = B_Mu_IsoTrackCorr1 - M2.Pt();
        //if(P1.Pt()/Mu_TrakIso[0] >10.) Mu_TIcorr = Mu_TrakIso[0] ;
      }
      if (DeltaR13 < 0.3) {
        B_Mu_IsoTrackCorr1 = B_Mu_IsoTrackCorr1 - M3.Pt();
      }
      if (DeltaR14 < 0.3) {
        B_Mu_IsoTrackCorr1 = B_Mu_IsoTrackCorr1 - M4.Pt();
      }

      // FourL_rapidity=FourMu.Rapidity();
      // Ups1To2_dY = abs(Ups_Rapidity2-Ups_Rapidity1); // ?

      //Calculation of Isolation 2
      //Calculate Delta R between the muons
      float DeltaR21 = M2.DeltaR(M1);
      float DeltaR23 = M2.DeltaR(M3);
      float DeltaR24 = M2.DeltaR(M4);
      float B_Mu_IsoTrackCorr2 = B_Mu2_IsoTrack->at(i);
      if (DeltaR21 < 0.3) {
        //      cout<<" Himal deltaR12 "<< DeltaR12<< "  " << Mu_TrakIso[0] << " - " << P1.Pt() << endl;
        // really subtract?
        B_Mu_IsoTrackCorr2 = B_Mu_IsoTrackCorr2 - M1.Pt();
        //if(P1.Pt()/Mu_TrakIso[0] >10.) Mu_TIcorr = Mu_TrakIso[0] ;
      }
      if (DeltaR23 < 0.3) {
        B_Mu_IsoTrackCorr2 = B_Mu_IsoTrackCorr2 - M3.Pt();
      }
      if (DeltaR24 < 0.3) {
        B_Mu_IsoTrackCorr2 = B_Mu_IsoTrackCorr2 - M4.Pt();
      }

      // FourL_rapidity=FourMu.Rapidity();
      // Ups1To2_dY = abs(Ups_Rapidity2-Ups_Rapidity1); // ?

      /*
      //Calculation of Isolation 3
      //Calculate Delta R between the muons
      float DeltaR32 = M3.DeltaR(M2);
      float DeltaR31 = M3.DeltaR(M1);
      float DeltaR34 = M3.DeltaR(M4);
      float B_Mu_IsoTrackCorr3 = B_Mu3_IsoTrack->at(i);
      if (DeltaR32<0.3) {
	//      cout<<" Himal deltaR12 "<< DeltaR12<< "  " << Mu_TrakIso[0] << " - " << P1.Pt() << endl;
	// really subtract?                                                                                                                                                                             
	B_Mu_IsoTrackCorr3 = B_Mu_IsoTrackCorr3 - M2.Pt();	
	//if(P1.Pt()/Mu_TrakIso[0] >10.) Mu_TIcorr = Mu_TrakIso[0] ;
      }
      if (DeltaR31<0.3) {
	B_Mu_IsoTrackCorr3 = B_Mu_IsoTrackCorr3 - M1.Pt();	
      }
      if (DeltaR34<0.3) {
	B_Mu_IsoTrackCorr3 = B_Mu_IsoTrackCorr3 - M4.Pt();	
      }
      
      // FourL_rapidity=FourMu.Rapidity();
      // Ups1To2_dY = abs(Ups_Rapidity2-Ups_Rapidity1); // ?


      //Calculation of Isolation 4
      //Calculate Delta R between the muons
      float DeltaR42 = M4.DeltaR(M2);
      float DeltaR43 = M4.DeltaR(M3);
      float DeltaR41 = M4.DeltaR(M1);
      float B_Mu_IsoTrackCorr4 = B_Mu4_IsoTrack->at(i);
      
      if (DeltaR42<0.3) {
	//      cout<<" Himal deltaR12 "<< DeltaR12<< "  " << Mu_TrakIso[0] << " - " << P1.Pt() << endl;
	// really subtract?                                                                                                                                                                             
	B_Mu_IsoTrackCorr4 = B_Mu_IsoTrackCorr4 - M2.Pt();	
	//if(P1.Pt()/Mu_TrakIso[0] >10.) Mu_TIcorr = Mu_TrakIso[0] ;
      }
      if (DeltaR43<0.3) {
	B_Mu_IsoTrackCorr4 = B_Mu_IsoTrackCorr4 - M3.Pt();	
      }
      if (DeltaR41<0.3) {
	B_Mu_IsoTrackCorr4 = B_Mu_IsoTrackCorr4 - M1.Pt();	
      }
      */
      FourL_rapidity = FourMu.Rapidity();
      Ups1To2_dY = abs(Ups_Rapidity2 - Ups_Rapidity1);

      //isolation cuts

      if ((B_Mu_IsoTrackCorr1 / B_Mu1_pt->at(i)) > 0.50)
        continue;
      ncandiIsolation1++;
      //Events=Event->at(i);
      if (temp_eventIsolation1 != Event->at(i)) {
        temp_eventIsolation1 = Event->at(i);
        nevIsolation1++;
      }

      // if ( (B_Mu_IsoTrackCorr2/B_Mu2_pt->at(i) ) >0.50) continue;

      ncandiIsolation2++;
      //Events=Event->at(i);
      if (temp_eventIsolation2 != Event->at(i)) {
        temp_eventIsolation2 = Event->at(i);
        nevIsolation2++;
      }
      /*
      if ( (B_Mu_IsoTrackCorr3/B_Mu3_pt->at(i) ) >0.35) continue;

      ncandiIsolation3++;
      //Events=Event->at(i);
      if (temp_eventIsolation3 != Event->at(i)  ){
	      temp_eventIsolation3= Event->at(i);
	      nevIsolation3++;
      }



      if ( (B_Mu_IsoTrackCorr4/B_Mu4_pt->at(i) ) >0.35) continue;
      ncandiIsolation4++;
      //Events=Event->at(i);
      if (temp_eventIsolation4 != Event->at(i)  ){
	      temp_eventIsolation4= Event->at(i);
	      nevIsolation4++;
      }

*/

      /*
      //cut 3 Muon Acceptance+ Blinding
      if (B_Mu1_pt->at(i) < 3.2 || B_Mu2_pt->at(i) < 3.2 || B_Mu3_pt->at(i) < 3.2 || B_Mu4_pt->at(i) < 3.2) continue;


      ncandiMuonPtGreaterThan++;
      //Events=Event->at(i);
      if (temp_eventMuonPtGreaterThan != Event->at(i)  ){
	      temp_eventMuonPtGreaterThan= Event->at(i);
	      nevMuonPtGreaterThan++;
      }
   */
      //blinding
      //if (FourL_mass->at(i) > 120. && FourL_mass->at(i)  < 130.) continue;
      //if (FourL_mass->at(i) > 70. && FourL_mass->at(i)  < 100.) continue;
      //cut 4 Z, J Vtx Prob
      if (Ups_VtxProb1 < 0.01 || Ups_VtxProb2 < 0.01)
        continue;

      ncandiYZVertexing++;
      //Events=Event->at(i);
      if (temp_eventYZVertexing != Event->at(i)) {
        temp_eventYZVertexing = Event->at(i);
        nevYZVertexing++;
      }
      //cut 4 J, Z Pt

      if (Ups_Pt1 < 5. || Ups_Pt2 < 5.)
        continue;

      ncandiYZPt++;
      //Events=Event->at(i);
      if (temp_eventYZPt != Event->at(i)) {
        temp_eventYZPt = Event->at(i);
        nevYZPt++;
      }

      //Dilepton mass cut 5
      if (Ups1_mass < 9.0 || Ups1_mass > 9.7)
        continue;
      if (Ups2_mass < 70.0 || Ups2_mass > 110)
        continue;
      // if (Ups2_mass < 85.0 || Ups2_mass>110) continue;

      ncandiYZMass++;
      //Events=Event->at(i);
      if (temp_eventYZMass != Event->at(i)) {
        temp_eventYZMass = Event->at(i);
        nevYZMass++;
      }

      //OnlyUps(1S)
      //if (Ups1_mass > 9.7) continue;
      //if (Ups2_mass > 9.7) continue;
      //rapidity cut cut 6
      // if (Ups1To2_dY > 3.0 ) continue;
      //delta phi cut cut 7
      // if (abs(Ups_Phi1-Ups_Phi2) < 1) continue;
      //4 lepton Vertex Pt cut 8
      if (FourL_pt->at(i) < 5.)
        continue;

      ncandi4Pt++;
      //Events=Event->at(i);
      if (temp_event4Pt != Event->at(i)) {
        temp_event4Pt = Event->at(i);
        nev4Pt++;
      }

      //cut 7 4 lepton Pt
      if (FourL_VtxProb->at(i) < 0.01)
        continue;
      //Four mu rapidity cut
      // if (abs(FourL_rapidity) > 1.7) continue;

      ncandi4Vertex++;
      //Events=Event->at(i);
      if (temp_event4Vertex != Event->at(i)) {
        temp_event4Vertex = Event->at(i);
        nev4Vertex++;
      }

      //FourL_mass Cut
      if (FourL_mass->at(i) < 112 || FourL_mass->at(i) > 162)
        continue;

      // blinding cut

      // if (FourL_mass->at(i) < 120 || FourL_mass->at(i) > 130) continue;

      ncandi4Mass++;
      //Events=Event->at(i);
      if (temp_event4Mass != Event->at(i)) {
        temp_event4Mass = Event->at(i);
        nev4Mass++;

        //continue;
      } else {
        // continue;
        doublecount++;
        myfile << doublecount << " double Event number detected # " << Event->at(i) << endl;
        myfile << " Event number of previous" << Event->at(i - 1) << endl;
        myfile << Event->at(i) << endl;

        cout << "Event no for multiple Candidate " << Event->at(i) << endl;
        myfile1 << Event->at(i) << endl;

        myfile << "Four muon vertex prob " << FourL_VtxProb->at(i) << endl;
        myfile << "Previous Four muon vertex prob " << FourL_VtxProb->at(i - 1) << endl << endl;

        if (FourL_VtxProb->at(i) > FourL_VtxProb->at(i - 1)) {
          cout << "thrown away the better" << endl;
        }
      }

      myfile << Event->at(i) << " " << Ups1_mass << " " << Ups2_mass << " " << FourL_mass->at(i) << " " << FourL_VtxProb->at(i) << endl;

      //myfile<<Ups1_mass<<endl;

      //if (ncandi>1) continue;
      Events_b->push_back(Events);
      B_Ups1_mass_b->push_back(Ups1_mass);
      B_Ups2_mass_b->push_back(Ups2_mass);
      B_Ups1_VtxProb_b->push_back(Ups_VtxProb1);
      B_Ups2_VtxProb_b->push_back(Ups_VtxProb2);
      B_Ups1To2dY_b->push_back(Ups1To2_dY);
      B_Ups1_Pt_b->push_back(Ups_Pt1);
      B_Ups1_Eta_b->push_back(Ups_Eta1);
      B_Ups1_Phi_b->push_back(Ups_Phi1);
      B_Ups2_Pt_b->push_back(Ups_Pt2);
      B_Ups2_Eta_b->push_back(Ups_Eta2);
      B_Ups2_Phi_b->push_back(Ups_Phi2);
      FourL_VtxProb_b->push_back(FourL_VtxProb->at(i));
      B_J1_mass_b->push_back(B_J1_mass->at(i));
      B_J2_mass_b->push_back(B_J2_mass->at(i));
      B_J3_mass_b->push_back(B_J3_mass->at(i));
      B_J4_mass_b->push_back(B_J4_mass->at(i));
      FourL_mass_b->push_back(FourL_mass->at(i));
      FourL_pt_b->push_back(FourL_pt->at(i));
      FourL_eta_b->push_back(FourL_eta->at(i));
      FourL_phi_b->push_back(FourL_phi->at(i));
      FourL_rapidity_b->push_back(FourL_rapidity);
      B_Mu1_pt_b->push_back(B_Mu1_pt->at(i));
      B_Mu2_pt_b->push_back(B_Mu2_pt->at(i));
      B_Mu3_pt_b->push_back(B_Mu3_pt->at(i));
      B_Mu4_pt_b->push_back(B_Mu4_pt->at(i));
      B_Mu1_eta_b->push_back(B_Mu1_eta->at(i));
      B_Mu2_eta_b->push_back(B_Mu2_eta->at(i));
      B_Mu3_eta_b->push_back(B_Mu3_eta->at(i));
      B_Mu4_eta_b->push_back(B_Mu4_eta->at(i));
      B_Mu1_phi_b->push_back(B_Mu1_phi->at(i));
      B_Mu2_phi_b->push_back(B_Mu2_phi->at(i));
      B_Mu3_phi_b->push_back(B_Mu3_phi->at(i));
      B_Mu4_phi_b->push_back(B_Mu4_phi->at(i));
      mu1mC2_b->push_back(mu1mC2->at(i));
      mu1pC2_b->push_back(mu1pC2->at(i));
      mu1mNHits_b->push_back(mu1mNHits->at(i));
      mu1pNHits_b->push_back(mu1pNHits->at(i));
      mu1mNPHits_b->push_back(mu1mNPHits->at(i));
      mu1pNPHits_b->push_back(mu1pNPHits->at(i));
      mu2mC2_b->push_back(mu2mC2->at(i));
      mu2pC2_b->push_back(mu2pC2->at(i));
      mu2mNHits_b->push_back(mu2mNHits->at(i));
      mu2pNHits_b->push_back(mu2pNHits->at(i));
      mu2mNPHits_b->push_back(mu2mNPHits->at(i));
      mu2pNPHits_b->push_back(mu2pNPHits->at(i));
      B_Mu1_soft_b->push_back(B_Mu1_soft->at(i));
      B_Mu1_loose_b->push_back(B_Mu1_loose->at(i));
      B_Mu1_tight_b->push_back(B_Mu1_tight->at(i));
      B_Mu2_soft_b->push_back(B_Mu2_soft->at(i));
      B_Mu2_loose_b->push_back(B_Mu2_loose->at(i));
      B_Mu2_tight_b->push_back(B_Mu2_tight->at(i));
      B_Mu3_soft_b->push_back(B_Mu3_soft->at(i));
      B_Mu3_loose_b->push_back(B_Mu3_loose->at(i));
      B_Mu3_tight_b->push_back(B_Mu3_tight->at(i));
      B_Mu4_soft_b->push_back(B_Mu4_soft->at(i));
      B_Mu4_loose_b->push_back(B_Mu4_loose->at(i));
      B_Mu4_tight_b->push_back(B_Mu4_tight->at(i));

      B_J_xyP1_b->push_back(B_J_xyP1->at(i));
      B_J_xyM1_b->push_back(B_J_xyM1->at(i));
      B_J_zP1_b->push_back(B_J_zP1->at(i));
      B_J_zM1_b->push_back(B_J_zM1->at(i));

      B_J_xyP2_b->push_back(B_J_xyP2->at(i));
      B_J_xyM2_b->push_back(B_J_xyM2->at(i));
      B_J_zP2_b->push_back(B_J_zP2->at(i));
      B_J_zM2_b->push_back(B_J_zM2->at(i));

      //B_Z_pt1_b->push_back(B_Z_pt1->at(i));
      //B_Z_pt2_b->push_back(B_Z_pt2->at(i));
      //B_Z_eta1_b->push_back(B_Z_eta1->at(i));
      //B_Z_eta2_b->push_back(B_Z_eta2->at(i));
      //B_Z_phi1_b->push_back(B_Z_phi1->at(i));
      //B_Z_phi2_b->push_back(B_Z_phi2->at(i));

      //B_J_pt1_b->push_back(B_J_pt1->at(i));
      //B_J_pt2_b->push_back(B_J_pt2->at(i));
      //B_J_eta1_b->push_back(B_J_eta1->at(i));
      //B_J_eta2_b->push_back(B_J_eta2->at(i));
      //B_J_phi1_b->push_back(B_J_phi1->at(i));
      //B_J_phi2_b->push_back(B_J_phi2->at(i));
    }
  }
  cout << "Number of event Preselection=" << nevPreSelection << endl;
  cout << "Number of Candidate Preselection=" << ncandiPreSelection << endl;

  cout << "Number of event Trigger=" << nevTrigger << endl;
  cout << "Number of Candidate Trigger=" << ncandiTrigger << endl;

  cout << "Number of event Soft=" << nevSoft << endl;
  cout << "Number of Candidate Soft=" << ncandiSoft << endl;
  cout << "Number of event Detector=" << nevDetector << endl;
  cout << "Number of Candidate Detector=" << ncandiDetector << endl;
  cout << "Number of event Lite4Vertex=" << nevLite4Vertex << endl;
  cout << "Number of Candidate Lite4Vertex=" << ncandiLite4Vertex << endl;
  cout << "Number of event Combine=" << nevCombine << endl;
  cout << "Number of Candidate Combine=" << ncandiCombine << endl;

  cout << "Number of event Isolation1=" << nevIsolation1 << endl;
  cout << "Number of Candidate Isolation1=" << ncandiIsolation1 << endl;

  cout << "Number of event Isolation2=" << nevIsolation2 << endl;
  cout << "Number of Candidate Isolation2=" << ncandiIsolation2 << endl;

  cout << "Number of event Isolation3=" << nevIsolation3 << endl;
  cout << "Number of Candidate Isolation3=" << ncandiIsolation3 << endl;

  cout << "Number of event Isolation4=" << nevIsolation4 << endl;
  cout << "Number of Candidate Isolation4=" << ncandiIsolation4 << endl;

  cout << "Number of event MuonPtGreaterThan=" << nevMuonPtGreaterThan << endl;
  cout << "Number of Candidate MuonPtGreaterThan=" << ncandiMuonPtGreaterThan << endl;
  cout << "Number of event YZVertexing=" << nevYZVertexing << endl;
  cout << "Number of Candidate YZVertexing=" << ncandiYZVertexing << endl;
  cout << "Number of event YZPt=" << nevYZPt << endl;
  cout << "Number of Candidate YZPt=" << ncandiYZPt << endl;
  cout << "Number of event YZMass=" << nevYZMass << endl;
  cout << "Number of Candidate YZMass=" << ncandiYZMass << endl;
  cout << "Number of event 4Pt=" << nev4Pt << endl;
  cout << "Number of Candidate 4Pt=" << ncandi4Pt << endl;
  cout << "Number of event 4Vertex=" << nev4Vertex << endl;
  cout << "Number of Candidate 4Vertex=" << ncandi4Vertex << endl;
  cout << "Number of event 4Mass=" << nev4Mass << endl;
  cout << "Number of Candidate 4Mass=" << ncandi4Mass << endl;

  if (ncandiPreSelection > 0 || ncandiTrigger > 0 || ncandiSoft > 0 || ncandiDetector > 0 || ncandiLite4Vertex > 0 || ncandiCombine > 0 ||
      ncandiMuonPtGreaterThan > 0 || ncandiYZVertexing > 0 || ncandiYZPt > 0 || ncandiYZMass > 0 || ncandi4Vertex > 0 || ncandi4Mass > 0) {
    fTree->Fill();
  }

  fFile->Write();
  myfile.close();
  myfile1.close();
  fFile->Close();

  delete fFile;
}
