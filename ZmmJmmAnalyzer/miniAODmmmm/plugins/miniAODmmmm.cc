// -*- C++ -*-
//
// Package:    miniAODmmmm
// Class:      miniAODmmmm
//
/**\class miniAODmmmm miniAODmmmm.cc ZmmJmmAnalyzer/miniAODmmmm/plugins/miniAODmmmm.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nimmitha Karunarathna
//         Created:  Tue, 21 May 2024 21:55:18 GMT
//
//

// system include files
#include <memory>

// user include files
#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/miniAODmmmm.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"
#include <TFile.h>
#include "TLorentzVector.h"

// Kinematic vertex fitter
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"

// // trigger
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

// // packedCandidates
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// // From Kalman example
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoVertex/KalmanVertexFit/interface/SimpleVertexTree.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include <iostream>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// #include <cmath>
// constants, enums and typedefs
//

typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//

miniAODmmmm::miniAODmmmm(const edm::ParameterSet &iConfig)
    : muonsToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      TriggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      prunedGenToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned"))),
      triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
      triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
      primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      // estoken_MF(esConsumes()),
      // estoken_TTB(esConsumes<TransientTrackBuilder, TransientTrackRecord>()),
      estoken_TTB(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      MuonTriggerString(iConfig.getParameter<std::string>("MuonTrigger")),

      isMC_(iConfig.getParameter<bool>("isMC")),

      tree_(0),

      Run(0),
      LumiBlock(0),
      Event(0),
      eventTime(0),
      TriggerFired(false),

      nPV(0),
      B_J1_mass(0),
      B_J1_pt(0),
      B_J1_rapidity(0),

      B_J1_VtxPt(0),
      B_J1_VtxMass(0),
      B_J1_VtxProb(0),

      mu1_B1(false), mu1_B2(false), mu1_B3(false), mu1_B4(false), mu1_F(false),
      mu2_B1(false), mu2_B2(false), mu2_B3(false), mu2_B4(false), mu2_F(false),

      B_Mu1_pt(0),
      B_Mu2_pt(0),
      B_Mu1_eta(0),
      B_Mu2_eta(0),
      isBestCandidate(false) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

miniAODmmmm::~miniAODmmmm() {}

//
// member functions
//

// ------------ method called for each event  ------------
void miniAODmmmm::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  const auto &theB = iSetup.getData(estoken_TTB);

  edm::Handle<pat::MuonCollection> thePATMuonHandle;
  iEvent.getByToken(muonsToken_, thePATMuonHandle);

  edm::Handle<edm::TriggerResults> TriggerResults;
  iEvent.getByToken(TriggerResultsToken_, TriggerResults);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  if (!thePATMuonHandle.isValid()) {
    edm::LogWarning("miniAODmmmm") << "No pat::Muon found on Event!";
    return;
  }
  if (!TriggerResults.isValid()) {
    edm::LogWarning("miniAODmmmm") << "No TriggerResults found on Event!";
    return;
  }
  if (!triggerPrescales.isValid()) {
    edm::LogWarning("miniAODmmmm") << "no Trigger prescale in event!";
    return;
  }

  const edm::TriggerNames &names = iEvent.triggerNames(*TriggerResults);

  // Check if the trigger fired
  bool triggerFlag = false;
  // Loop over all triggers to see if the desired one fired
  for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
    if (names.triggerName(i).find(MuonTriggerString.c_str()) != string::npos) {
      if (TriggerResults->accept(i)) {
        triggerFlag = true;  // Trigger fired!
      }
      break;  // No need to check further if we already found the trigger
    }
  }
  if (!triggerFlag) {
    return;
  }

  // Struct to hold dimuon candidates
  struct DimuonCandidate {
    const pat::Muon *mu1;
    const pat::Muon *mu2;
    float vtxProb;
    TLorentzVector p4;
    // XYZTLorentzVectorD J_vtx;
    float mu1_pt, mu2_pt, mu1_eta, mu2_eta;

    bool operator<(const DimuonCandidate &other) const {
      return vtxProb > other.vtxProb;  // sort descending by vtxProb
    }
  };

  std::vector<DimuonCandidate> validCandidates;

  //*********************************
  //Now we get the primary vertex
  //*********************************

  reco::Vertex bestVtx;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  bestVtx = *(primaryVertices_handle->begin());

  UShort_t goodPVCount = 0;

  if (primaryVertices_handle.isValid()) {
    for (const auto &vtx : *primaryVertices_handle) {
      if (!vtx.isFake() && vtx.ndof() > 4 && fabs(vtx.z()) <= 24.0 && fabs(vtx.position().Rho()) <= 2.0) {
        ++goodPVCount;
      }
    }
    if (!primaryVertices_handle->empty()) {
      bestVtx = *(primaryVertices_handle->begin());
    }
  } else {
    edm::LogWarning("miniAODmmmm") << "Primary vertex collection is not valid!";
  }

  nPV = goodPVCount;

  //****************************************************
  //*********Now we get the muons***********************
  //****************************************************
  for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
    for (pat::MuonCollection::const_iterator iMuon2 = iMuon1 + 1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
      //make sure all muons are diferent
      if (iMuon1 == iMuon2)
        continue;

      if (!(abs(iMuon1->charge()) == 1))
        continue;
      if (!(abs(iMuon2->charge()) == 1))
        continue;

      //opposite charge
      if (!((iMuon1->charge()) + (iMuon2->charge())) == 0)
        continue;

      if (iMuon1->pt() < 3.0)
        continue;
      if (iMuon2->pt() < 3.0)
        continue;

      if (iMuon1->eta() < -2.4 || iMuon1->eta() > 2.4)
        continue;
      if (iMuon2->eta() < -2.4 || iMuon2->eta() > 2.4)
        continue;

      TrackRef glbTrackP1;
      TrackRef glbTrackM1;

      if (iMuon1->charge() == 1 && iMuon2->charge() == -1) {
        glbTrackP1 = iMuon1->track();
        glbTrackM1 = iMuon2->track();
      } else if (iMuon1->charge() == -1 && iMuon2->charge() == 1) {
        glbTrackP1 = iMuon2->track();
        glbTrackM1 = iMuon1->track();
      } else {
        cout << "Something is wrong while making glb track ref" << endl;
      }

      if (glbTrackP1.isNull() || glbTrackM1.isNull()) {
        //std::cout << "continue due to no track ref" << endl;
        continue;
      }

      TLorentzVector M1, M2, MM1;
      float mu_mass = 0.1056583745;  //[PDG mass]

      //make muon 4 vectors
      if (iMuon1->charge() == 1 && iMuon2->charge() == -1) {
        M1.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
        M2.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
      } else if (iMuon1->charge() == -1 && iMuon2->charge() == 1) {
        M1.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
        M2.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
      } else {
        cout << "Something is wrong while making charge lorentz vector" << endl;
      }

      //make netral dimuon combination
      MM1 = M1 + M2;

      //cout<<"Start looking muon track quality"<<endl;
      if (!(glbTrackM1->quality(reco::TrackBase::highPurity)))
        continue;
      if (!(glbTrackP1->quality(reco::TrackBase::highPurity)))
        continue;

      reco::TransientTrack muon1TT;
      reco::TransientTrack muon2TT;

      muon1TT = theB.build(glbTrackP1);
      muon2TT = theB.build(glbTrackM1);

      //Kalman Vtx----------------------//
      vector<TransientTrack> mu_tks;

      KalmanVertexFitter kvfM(true);
      mu_tks.clear();
      mu_tks.push_back(muon1TT);
      mu_tks.push_back(muon2TT);
      TransientVertex J_candi1 = kvfM.vertex(mu_tks);

      if (!J_candi1.isValid()) {
        // cout << "continue because no vertexed dimuon" << endl;
        continue;
      }

      reco::Vertex JPsi_Vtx1 = J_candi1;

      float B_Prob_tmp1 = TMath::Prob(J_candi1.totalChiSquared(), J_candi1.degreesOfFreedom());
      // const math::XYZTLorentzVectorD JPsi_mom1 = JPsi_Vtx1.p4(mu_mass, 0.0);

      if (B_Prob_tmp1 < 0.1) {
        continue;
      }

      // Remove events with mass outside J/Psi or Z mass window
      if ((MM1.M() < 2.6 || MM1.M() > 3.6)) {
        continue;
      }

      validCandidates.push_back({&(*iMuon1),
                                 &(*iMuon2),
                                 B_Prob_tmp1,
                                 MM1,
                                 //  JPsi_mom1,
                                 static_cast<float>(iMuon1->pt()),
                                 static_cast<float>(iMuon2->pt()),
                                 static_cast<float>(iMuon1->eta()),
                                 static_cast<float>(iMuon2->eta())});
    }
  }

  //Event Information
  Run = iEvent.id().run();
  LumiBlock = iEvent.luminosityBlock();
  Event = iEvent.id().event();

  edm::Timestamp timestamp = iEvent.eventAuxiliary().time();
  unsigned int seconds = timestamp.unixTime();
  unsigned int microseconds = timestamp.microsecondOffset();

  unsigned long long milliseconds = static_cast<unsigned long long>(seconds) * 1000 + static_cast<unsigned long long>(microseconds) / 1000;

  eventTime = milliseconds;

  // cout << "milliseconds: " << milliseconds << endl;
  // cout << "Run: " << Run << " LumiBlock: " << LumiBlock << " Event: " << Event << " eventTime: " << eventTime << endl;
  TriggerFired = triggerFlag;

  std::sort(validCandidates.begin(), validCandidates.end());

  // Track used muons
  std::set<const pat::Muon *> usedMuons;
  usedMuons.clear();

  UShort_t nRecoDistinctJWVtx = 0;
  for (const auto &cand : validCandidates) {
    if (usedMuons.count(cand.mu1) || usedMuons.count(cand.mu2)) {
      continue;
    }
    nRecoDistinctJWVtx++;
    usedMuons.insert(cand.mu1);
    usedMuons.insert(cand.mu2);

    mu1_B1 = mu1_B2 = mu1_B3 = mu1_B4 = mu1_F = false;
    mu2_B1 = mu2_B2 = mu2_B3 = mu2_B4 = mu2_F = false;

    // Analyze hit patterns for both muons
    const reco::Track* tk1 = cand.mu1->innerTrack().get();
    const reco::Track* tk2 = cand.mu2->innerTrack().get();
    if (!tk1 || !tk2) continue;

    const reco::HitPattern& hp1 = tk1->hitPattern();
    const reco::HitPattern& hp2 = tk2->hitPattern();

    for (int i = 0; i < hp1.numberOfAllHits(reco::HitPattern::TRACK_HITS); i++) {
      uint32_t hit = hp1.getHitPattern(reco::HitPattern::TRACK_HITS, i);
      if (!hp1.validHitFilter(hit)) continue;

      int subdet = hp1.getSubStructure(hit);
      int layer  = hp1.getLayer(hit);

      if (subdet == 1) { // Pixel barrel
          if (layer == 1) mu1_B1 = true;
          else if (layer == 2) mu1_B2 = true;
          else if (layer == 3) mu1_B3 = true;
          else if (layer == 4) mu1_B4 = true;
      } else if (subdet == 2) { // Pixel forward
          mu1_F = true;
      }
    }

    for (int i = 0; i < hp2.numberOfAllHits(reco::HitPattern::TRACK_HITS); i++) {
      uint32_t hit = hp2.getHitPattern(reco::HitPattern::TRACK_HITS, i);
      if (!hp2.validHitFilter(hit)) continue;

      int subdet = hp2.getSubStructure(hit);
      int layer  = hp2.getLayer(hit);

      if (subdet == 1) {
          if (layer == 1) mu2_B1 = true;
          else if (layer == 2) mu2_B2 = true;
          else if (layer == 3) mu2_B3 = true;
          else if (layer == 4) mu2_B4 = true;
      } else if (subdet == 2) {
          mu2_F = true;
      }
    }

    // Fill the tree with all the candidates
    B_J1_mass = cand.p4.M();
    B_J1_pt = std::lround(cand.p4.Pt() * 1000);
    B_J1_rapidity = std::lround(cand.p4.Rapidity() * 1000);
    // B_J1_VtxPt = cand.J_vtx.Pt();
    // B_J1_VtxMass = cand.J_vtx.mass();
    B_J1_VtxProb = cand.vtxProb;
    B_Mu1_pt = std::lround(cand.mu1->pt() * 1000);
    B_Mu2_pt = std::lround(cand.mu2->pt() * 1000);
    B_Mu1_eta = std::lround(cand.mu1->eta() * 100);
    B_Mu2_eta = std::lround(cand.mu2->eta() * 100);

    if (nRecoDistinctJWVtx == 1) {
      isBestCandidate = true;  // Mark the first candidate as the best
    } else {
      isBestCandidate = false;  // Subsequent candidates are not the best
    }

    tree_->Fill();

    B_J1_mass = -999;
    B_J1_pt = 1000;
    B_J1_rapidity = -999;

    B_J1_VtxPt = -999;
    B_J1_VtxMass = -999;
    B_J1_VtxProb = -999;

    B_Mu1_pt = 1000;
    B_Mu2_pt = 1000;
    B_Mu1_eta = 10;
    B_Mu2_eta = 10;
  }

  // Reset variables for the next event
  Run = 0;
  LumiBlock = 0;
  Event = 0;
  eventTime = 0;
  TriggerFired = false;
  nPV = 0;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------

void miniAODmmmm::beginJob() {
  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  tree_ = new TTree("ntuple", "ntuple");

  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  // tree_->Branch("Event", &Event);
  tree_->Branch("eventTime", &eventTime);
  // tree_->Branch("TriggerFired", &TriggerFired);

  tree_->Branch("nPV", &nPV);
  tree_->Branch("B_J1_mass", &B_J1_mass);
  tree_->Branch("B_J1_pt", &B_J1_pt);
  // tree_->Branch("B_J1_rapidity", &B_J1_rapidity);

  // tree_->Branch("B_J1_VtxPt", &B_J1_VtxPt);
  // tree_->Branch("B_J1_VtxMass", &B_J1_VtxMass);
  // tree_->Branch("B_J1_VtxProb", &B_J1_VtxProb);

  // Muon 1
  tree_->Branch("mu1_B1", &mu1_B1, "mu1_B1/O");
  tree_->Branch("mu1_B2", &mu1_B2, "mu1_B2/O");
  tree_->Branch("mu1_B3", &mu1_B3, "mu1_B3/O");
  tree_->Branch("mu1_B4", &mu1_B4, "mu1_B4/O");
  tree_->Branch("mu1_F",  &mu1_F,  "mu1_F/O");

  // Muon 2
  tree_->Branch("mu2_B1", &mu2_B1, "mu2_B1/O");
  tree_->Branch("mu2_B2", &mu2_B2, "mu2_B2/O");
  tree_->Branch("mu2_B3", &mu2_B3, "mu2_B3/O");
  tree_->Branch("mu2_B4", &mu2_B4, "mu2_B4/O");
  tree_->Branch("mu2_F",  &mu2_F,  "mu2_F/O");

  tree_->Branch("B_Mu1_pt", &B_Mu1_pt);
  tree_->Branch("B_Mu2_pt", &B_Mu2_pt);
  tree_->Branch("B_Mu1_eta", &B_Mu1_eta);
  tree_->Branch("B_Mu2_eta", &B_Mu2_eta);
  tree_->Branch("isBestCandidate", &isBestCandidate);
}

// ------------ method called once each job just after ending the event loop  ------------
void miniAODmmmm::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void miniAODmmmm::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  //  Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(miniAODmmmm);