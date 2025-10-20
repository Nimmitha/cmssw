// system include files
#include <memory>

// user include files
#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/TriggerEffAnalyzer.h"

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

TriggerEffAnalyzer::TriggerEffAnalyzer(const edm::ParameterSet &iConfig)
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

      B_Mu1_pt(0),
      B_Mu1_eta(0),
      B_Mu2_pt(0),
      B_Mu2_eta(0),
      firedHLT_L1DoubleMu(false),
      hasDimuonInRange(false),
      dm_triggerMatched_m1(false),
      dm_triggerMatched_m2(false) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

TriggerEffAnalyzer::~TriggerEffAnalyzer() {}

//
// member functions
//

// ------------ method called for each event  ------------
void TriggerEffAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
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
    edm::LogWarning("TriggerEffAnalyzer") << "No pat::Muon found on Event!";
    return;
  }
  if (!TriggerResults.isValid()) {
    edm::LogWarning("TriggerEffAnalyzer") << "No TriggerResults found on Event!";
    return;
  }
  if (!triggerPrescales.isValid()) {
    edm::LogWarning("TriggerEffAnalyzer") << "no Trigger prescale in event!";
    return;
  }

  Run = iEvent.id().run();
  LumiBlock = iEvent.id().luminosityBlock();
  Event = iEvent.id().event();

  const edm::TriggerNames &names = iEvent.triggerNames(*TriggerResults);

  const pat::Muon *sm_muon_ptr = nullptr;
  float B_Mu1_phi = -99;

  B_Mu2_pt = -1;
  B_Mu2_eta = -99;

  // -------------------------------
  // Selection 1:
  // -------------------------------
  bool firedHLT_Barrel = false;
  for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
    std::string name = names.triggerName(i);
    if (name.find("HLT_Mu0_Barrel_v") != std::string::npos && TriggerResults->accept(i)) {
      firedHLT_Barrel = true;
      break;
    }
  }

  if (!firedHLT_Barrel)
    return;

  bool sm_triggerMatched = false;
  for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
    if (iMuon1->pt() < 5.0)
      continue;

    if (abs(iMuon1->eta()) > 2.4)
      continue;

    // Find the trigger object that matches to this muon
    for (const auto &obj : *triggerObjects) {
      pat::TriggerObjectStandAlone unpackedObj = obj;
      unpackedObj.unpackPathNames(names);

      if (unpackedObj.hasPathName("HLT_Mu0_Barrel_v*", true, true)) {
        float dR_temp = reco::deltaR(iMuon1->eta(), iMuon1->phi(), unpackedObj.eta(), unpackedObj.phi());
        if (dR_temp < 0.01) {
          sm_triggerMatched = true;

          sm_muon_ptr = &(*iMuon1);
          B_Mu1_pt = iMuon1->pt();
          B_Mu1_eta = iMuon1->eta();
          B_Mu1_phi = iMuon1->phi();

          float_t sm_obj_pt = unpackedObj.pt();
          float_t sm_obj_eta = unpackedObj.eta();
          float_t sm_obj_phi = unpackedObj.phi();
        }
      }
    }
  }

  if (!sm_triggerMatched)
    return;

  // --------------------------------------
  // Selection 2:
  // --------------------------------------

  // --------------------------------------
  // 1. Check if event fired HLT_Mu0_L1DoubleMu_v*
  // --------------------------------------
  firedHLT_L1DoubleMu = false;
  hasDimuonInRange = false;
  dm_triggerMatched_m1 = false;
  dm_triggerMatched_m2 = false;
  for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
    std::string name = names.triggerName(i);
    if (name.find("HLT_Mu0_L1DoubleMu_v") != std::string::npos && TriggerResults->accept(i)) {
      firedHLT_L1DoubleMu = true;
      break;
    }
  }

  // --------------------------------------
  // 2. Build dimuon with another muon
  // --------------------------------------
  if (firedHLT_L1DoubleMu) {
    for (pat::MuonCollection::const_iterator iMuon2 = thePATMuonHandle->begin(); iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
      // Skip same muon (pointer comparison)
      if (&(*iMuon2) == sm_muon_ptr)
        continue;

      if (iMuon2->pt() < 5.0)
        continue;

      if (abs(iMuon2->eta()) > 2.4)
        continue;

      TLorentzVector M1, M2, MM;
      float mu_mass = 0.1056583745;  //[PDG mass]

      M1.SetPtEtaPhiM(B_Mu1_pt, B_Mu1_eta, B_Mu1_phi, mu_mass);
      M2.SetPtEtaPhiM(iMuon2->pt(), iMuon2->eta(), iMuon2->phi(), mu_mass);

      MM = M1 + M2;

      if (MM.M() > 2.0 && MM.M() < 4.0) {
        hasDimuonInRange = true;
        B_Mu2_pt = iMuon2->pt();
        B_Mu2_eta = iMuon2->eta();

        // check trigger object matching
        for (const auto &obj : *triggerObjects) {
          pat::TriggerObjectStandAlone unpackedObj = obj;
          unpackedObj.unpackPathNames(names);

          if (unpackedObj.hasPathName("HLT_Mu0_L1DoubleMu_v*", true, true)) {
            if (reco::deltaR(sm_muon_ptr->eta(), sm_muon_ptr->phi(), unpackedObj.eta(), unpackedObj.phi()) < 0.01) {
              dm_triggerMatched_m1 = true;
            }
            if (reco::deltaR(iMuon2->eta(), iMuon2->phi(), unpackedObj.eta(), unpackedObj.phi()) < 0.01) {
              dm_triggerMatched_m2 = true;
            }
            break;
          }
        }

        break;  // break after first dimuon candidate found
      }
    }
  }

  tree_->Fill();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------

void TriggerEffAnalyzer::beginJob() {
  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  tree_ = new TTree("ntuple", "ntuple");

  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  // tree_->Branch("Event", &Event);

  tree_->Branch("B_Mu1_pt", &B_Mu1_pt);
  tree_->Branch("B_Mu1_eta", &B_Mu1_eta);
  tree_->Branch("B_Mu2_pt", &B_Mu2_pt);
  tree_->Branch("B_Mu2_eta", &B_Mu2_eta);
  tree_->Branch("firedHLT_L1DoubleMu", &firedHLT_L1DoubleMu);
  tree_->Branch("hasDimuonInRange", &hasDimuonInRange);
  tree_->Branch("dm_triggerMatched_m1", &dm_triggerMatched_m1);
  tree_->Branch("dm_triggerMatched_m2", &dm_triggerMatched_m2);
}

// ------------ method called once each job just after ending the event loop  ------------
void TriggerEffAnalyzer::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TriggerEffAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  //  Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(TriggerEffAnalyzer);