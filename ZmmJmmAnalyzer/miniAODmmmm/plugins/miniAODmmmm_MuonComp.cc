// -*- C++ -*-
//
// Package:    miniAODmmmm_MuonComp
// Class:      miniAODmmmm_MuonComp
//
/**\class miniAODmmmm_MuonComp miniAODmmmm_MuonComp.cc ZmmJmmAnalyzer/miniAODmmmm_MuonComp/plugins/miniAODmmmm_MuonComp.cc

 Description: v6 J/psi selection (miniAODmmmm, production_reduced_size_v6) with per-muon type, ID and
              inner-track hit information, for the muon composition and pixel-layer study (track_hits_v2).

 Implementation:
     Selection is identical to v6: no muon ID cut. Soft, loose, medium and tight ID are stored as flags
     (mu*_type), so "v6 as run" and "v6 + soft ID" both come from the same ntuple.
     Per muon: type/ID bits, valid and missing hits per pixel layer/disk, first valid hit,
     pixel/tracker layers with measurement and missing inner hits (see miniAODmmmm_MuonComp.h).
*/
//
// Original Author:  Nimmitha Karunarathna
//         Created:  Tue, 21 May 2024 21:55:18 GMT
//
//

// system include files
#include <memory>

// user include files
#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/miniAODmmmm_MuonComp.h"

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
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include <algorithm>

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

miniAODmmmm_MuonComp::miniAODmmmm_MuonComp(const edm::ParameterSet &iConfig)
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

      mu1_type(0), mu2_type(0),
      mu1_pixValid(0), mu2_pixValid(0),
      mu1_pixLost(0), mu2_pixLost(0),
      mu1_firstHit(0), mu2_firstHit(0),
      mu1_nPixLay(0), mu2_nPixLay(0),
      mu1_nTrkLay(0), mu2_nTrkLay(0),
      mu1_nLostInner(0), mu2_nLostInner(0),
      mu1_nLostInnerPix(0), mu2_nLostInnerPix(0),

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

miniAODmmmm_MuonComp::~miniAODmmmm_MuonComp() {}

// Type and ID bits: 0 global, 1 tracker, 2 standalone, 3 PF, 4 soft, 5 loose, 6 medium, 7 tight.
// Soft and tight use the same vertex as the v6 code (first offline primary vertex).
UChar_t miniAODmmmm_MuonComp::muonType(const pat::Muon &mu, const reco::Vertex &vtx) const {
  UChar_t t = 0;
  if (mu.isGlobalMuon())
    t |= 1 << 0;
  if (mu.isTrackerMuon())
    t |= 1 << 1;
  if (mu.isStandAloneMuon())
    t |= 1 << 2;
  if (mu.isPFMuon())
    t |= 1 << 3;
  if (mu.isSoftMuon(vtx))
    t |= 1 << 4;
  if (mu.isLooseMuon())
    t |= 1 << 5;
  if (mu.isMediumMuon())
    t |= 1 << 6;
  if (mu.isTightMuon(vtx))
    t |= 1 << 7;
  return t;
}

// Hit bits: 0-3 BPix layers 1-4, 4-6 FPix disks 1-3.
static int pixelBit(uint32_t subdet, uint32_t layer) {
  if (subdet == PixelSubdetector::PixelBarrel && layer >= 1 && layer <= 4)
    return layer - 1;
  if (subdet == PixelSubdetector::PixelEndcap && layer >= 1 && layer <= 3)
    return 3 + layer;
  return -1;
}

void miniAODmmmm_MuonComp::fillHits(const reco::Track &tk, UChar_t &valid, UChar_t &lost, UChar_t &first, UChar_t &nPixLay,
                                    UChar_t &nTrkLay, UChar_t &nLostInner, UChar_t &nLostInnerPix) const {
  const reco::HitPattern &hp = tk.hitPattern();
  valid = lost = first = 0;
  for (int i = 0; i < hp.numberOfAllHits(reco::HitPattern::TRACK_HITS); i++) {
    uint16_t hit = hp.getHitPattern(reco::HitPattern::TRACK_HITS, i);
    if (!reco::HitPattern::trackerHitFilter(hit))
      continue;
    uint32_t subdet = reco::HitPattern::getSubStructure(hit);
    uint32_t layer = reco::HitPattern::getLayer(hit);
    if (reco::HitPattern::validHitFilter(hit)) {
      if (first == 0)
        first = 10 * subdet + layer;
      int b = pixelBit(subdet, layer);
      if (b >= 0)
        valid |= 1 << b;
    } else if (reco::HitPattern::missingHitFilter(hit)) {
      int b = pixelBit(subdet, layer);
      if (b >= 0)
        lost |= 1 << b;
    }
  }
  nPixLay = hp.pixelLayersWithMeasurement();
  nTrkLay = hp.trackerLayersWithMeasurement();
  nLostInner = std::min(hp.numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS), 255);
  nLostInnerPix = std::min(hp.numberOfLostPixelHits(reco::HitPattern::MISSING_INNER_HITS), 255);
}

//
// member functions
//

// ------------ method called for each event  ------------
void miniAODmmmm_MuonComp::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
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
    edm::LogWarning("miniAODmmmm_MuonComp") << "No pat::Muon found on Event!";
    return;
  }
  if (!TriggerResults.isValid()) {
    edm::LogWarning("miniAODmmmm_MuonComp") << "No TriggerResults found on Event!";
    return;
  }
  if (!triggerPrescales.isValid()) {
    edm::LogWarning("miniAODmmmm_MuonComp") << "no Trigger prescale in event!";
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

  if (primaryVertices_handle.isValid() && !primaryVertices_handle->empty()) {
    for (const auto &vtx : *primaryVertices_handle) {
      if (!vtx.isFake() && vtx.ndof() > 4 && fabs(vtx.z()) <= 24.0 && fabs(vtx.position().Rho()) <= 2.0) {
        ++goodPVCount;
      }
    }
    bestVtx = *(primaryVertices_handle->begin());
  } else {
    edm::LogWarning("miniAODmmmm_MuonComp") << "Primary vertex collection is not valid!";
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
      if ( (iMuon1->charge()) + (iMuon2->charge()) != 0 )
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

    // Type/ID and inner-track hits for both muons (track() is non-null by the selection above)
    mu1_type = muonType(*cand.mu1, bestVtx);
    mu2_type = muonType(*cand.mu2, bestVtx);
    fillHits(*cand.mu1->track(), mu1_pixValid, mu1_pixLost, mu1_firstHit, mu1_nPixLay, mu1_nTrkLay, mu1_nLostInner,
             mu1_nLostInnerPix);
    fillHits(*cand.mu2->track(), mu2_pixValid, mu2_pixLost, mu2_firstHit, mu2_nPixLay, mu2_nTrkLay, mu2_nLostInner,
             mu2_nLostInnerPix);

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

void miniAODmmmm_MuonComp::beginJob() {
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

  // Per muon, bit layouts in miniAODmmmm_MuonComp.h
  tree_->Branch("mu1_type", &mu1_type, "mu1_type/b");
  tree_->Branch("mu2_type", &mu2_type, "mu2_type/b");
  tree_->Branch("mu1_pixValid", &mu1_pixValid, "mu1_pixValid/b");
  tree_->Branch("mu2_pixValid", &mu2_pixValid, "mu2_pixValid/b");
  tree_->Branch("mu1_pixLost", &mu1_pixLost, "mu1_pixLost/b");
  tree_->Branch("mu2_pixLost", &mu2_pixLost, "mu2_pixLost/b");
  tree_->Branch("mu1_firstHit", &mu1_firstHit, "mu1_firstHit/b");
  tree_->Branch("mu2_firstHit", &mu2_firstHit, "mu2_firstHit/b");
  tree_->Branch("mu1_nPixLay", &mu1_nPixLay, "mu1_nPixLay/b");
  tree_->Branch("mu2_nPixLay", &mu2_nPixLay, "mu2_nPixLay/b");
  tree_->Branch("mu1_nTrkLay", &mu1_nTrkLay, "mu1_nTrkLay/b");
  tree_->Branch("mu2_nTrkLay", &mu2_nTrkLay, "mu2_nTrkLay/b");
  tree_->Branch("mu1_nLostInner", &mu1_nLostInner, "mu1_nLostInner/b");
  tree_->Branch("mu2_nLostInner", &mu2_nLostInner, "mu2_nLostInner/b");
  tree_->Branch("mu1_nLostInnerPix", &mu1_nLostInnerPix, "mu1_nLostInnerPix/b");
  tree_->Branch("mu2_nLostInnerPix", &mu2_nLostInnerPix, "mu2_nLostInnerPix/b");

  tree_->Branch("B_Mu1_pt", &B_Mu1_pt);
  tree_->Branch("B_Mu2_pt", &B_Mu2_pt);
  tree_->Branch("B_Mu1_eta", &B_Mu1_eta);
  tree_->Branch("B_Mu2_eta", &B_Mu2_eta);
  tree_->Branch("isBestCandidate", &isBestCandidate);
}

// ------------ method called once each job just after ending the event loop  ------------
void miniAODmmmm_MuonComp::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void miniAODmmmm_MuonComp::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  //  Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(miniAODmmmm_MuonComp);