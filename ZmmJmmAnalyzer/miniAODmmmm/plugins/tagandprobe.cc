// system include files
#include <memory>

// user include files
#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/tagandprobe.h"

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

tagandprobe::tagandprobe(const edm::ParameterSet& iConfig)
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

      // B_Mu1_pt(0),
      // B_Mu1_eta(0),
      // B_Mu2_pt(0),
      // B_Mu2_eta(0),
      firedHLT_IsoMu24(false),
      sm_triggerMatched(false),
      nIsoMu24Objs(0),
      nIsoMu24MatchedObjs(0),
      matchMask(0),
      firedHLT_L1DoubleMu(false),
      dm_matched_offline(false),
      dm_matched_sm(false)
// dm_triggerMatched_m1(false),
// dm_triggerMatched_m2(false)
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

tagandprobe::~tagandprobe() {}

//
// member functions
//

// ------------ method called for each event  ------------
void tagandprobe::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  const auto& theB = iSetup.getData(estoken_TTB);

  edm::Handle<pat::MuonCollection> thePATMuonHandle;
  iEvent.getByToken(muonsToken_, thePATMuonHandle);

  edm::Handle<edm::TriggerResults> TriggerResults;
  iEvent.getByToken(TriggerResultsToken_, TriggerResults);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  if (!thePATMuonHandle.isValid()) {
    edm::LogWarning("tagandprobe") << "No pat::Muon found on Event!";
    return;
  }
  if (!TriggerResults.isValid()) {
    edm::LogWarning("tagandprobe") << "No TriggerResults found on Event!";
    return;
  }
  if (!triggerPrescales.isValid()) {
    edm::LogWarning("tagandprobe") << "no Trigger prescale in event!";
    return;
  }

  Run = iEvent.id().run();
  LumiBlock = iEvent.id().luminosityBlock();
  Event = iEvent.id().event();

  const edm::TriggerNames& names = iEvent.triggerNames(*TriggerResults);

  firedHLT_IsoMu24 = false;
  for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
    std::string name = names.triggerName(i);
    if (name.find("HLT_IsoMu24_v") != std::string::npos && TriggerResults->accept(i)) {
      firedHLT_IsoMu24 = true;
      break;
    }
  }

  if (!firedHLT_IsoMu24)
    return;

  // collect objects to a vector
  std::vector<pat::TriggerObjectStandAlone> isoMu24_objects;

  for (const auto& obj : *triggerObjects) {
    pat::TriggerObjectStandAlone unpackedObj = obj;
    unpackedObj.unpackPathNames(names);
    unpackedObj.unpackFilterLabels(iEvent, *TriggerResults);

    if (unpackedObj.hasPathName("HLT_IsoMu24_v*", true, true)) {
      isoMu24_objects.push_back(unpackedObj);
    }
  }

  sm_triggerMatched = false;
  nIsoMu24Objs = isoMu24_objects.size();
  nIsoMu24MatchedObjs = 0;
  matchMask = 0;
  vector<float> mu1_pt;
  vector<float> mu1_eta;
  float dR_temp = 0.0;

  // check if each  isoMu24 object matches to any slimmed muon individually
  for (size_t i = 0; i < isoMu24_objects.size(); ++i) {
    for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
      dR_temp = reco::deltaR(iMuon1->eta(), iMuon1->phi(), isoMu24_objects[i].eta(), isoMu24_objects[i].phi());
      if (dR_temp < 0.1) {
        nIsoMu24MatchedObjs++;
        matchMask |= (1 << i);

        sm_triggerMatched = true;

        mu1_pt.push_back(iMuon1->pt());
        mu1_eta.push_back(iMuon1->eta());

        // B_Mu1_pt = iMuon1->pt();
        // B_Mu1_eta = iMuon1->eta();

        break;
      }
    }
  }
  // // cout << "nIsoMu24Objs: " << nIsoMu24Objs << ", nIsoMu24MatchedObjs: " << nIsoMu24MatchedObjs << ", matchMask: " << matchMask << " - " << std::bitset<8>(matchMask) << std::endl;

  // // --------------------------------------
  // // Selection 2:
  // // --------------------------------------
  firedHLT_L1DoubleMu = false;

  // // 1. Check if event fired HLT_Mu0_L1DoubleMu_v*
  for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
    std::string name = names.triggerName(i);
    if (name.find("HLT_Mu0_L1DoubleMu_v") != std::string::npos && TriggerResults->accept(i)) {
      firedHLT_L1DoubleMu = true;
      break;
    }
  }

  dm_matched_offline = false;
  dm_matched_sm = false;
  if (firedHLT_L1DoubleMu == true) {
    std::vector<pat::TriggerObjectStandAlone> DoubleMu_objects;
    for (const auto& obj : *triggerObjects) {
      pat::TriggerObjectStandAlone unpackedObj = obj;
      unpackedObj.unpackPathNames(names);
      unpackedObj.unpackFilterLabels(iEvent, *TriggerResults);

      if (unpackedObj.hasFilterLabel("hltL3fL1sDoubleMu0SQL1f0L2PreFilteres0L3Filtered0")) {
        DoubleMu_objects.push_back(unpackedObj);
      }
    }

    if (DoubleMu_objects.size() < 2) {
      cout << "Recovering the two muons" << std::endl;
      DoubleMu_objects.clear();
      for (const auto& obj : *triggerObjects) {
        pat::TriggerObjectStandAlone unpackedObj = obj;
        unpackedObj.unpackPathNames(names);
        unpackedObj.unpackFilterLabels(iEvent, *TriggerResults);

        if (unpackedObj.hasFilterLabel("hltL2fL1sL1DoubleMuL1f0L2PreFiltered0ForLowMassInclusive")) {
          DoubleMu_objects.push_back(unpackedObj);
        }
      }
      if (DoubleMu_objects.size() >= 2) {
        cout << "Successfully recovered the two muons" << std::endl;
      } else {
        cout << "Failed to recover the two muons" << std::endl;
      }
    }

    // Check if any of the two muons match to any of the slimmed muons and to the IsoMu24 trigger objects
    for (size_t i = 0; i < DoubleMu_objects.size(); ++i) {
      for (pat::MuonCollection::const_iterator iMuon2 = thePATMuonHandle->begin(); iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
        dR_temp = reco::deltaR(iMuon2->eta(), iMuon2->phi(), DoubleMu_objects[i].eta(), DoubleMu_objects[i].phi());
        if (dR_temp < 0.1) {
          dm_matched_offline = true;

          // Check if this muon also matched to the IsoMu24 trigger object
          for (size_t j = 0; j < isoMu24_objects.size(); ++j) {
            dR_temp = reco::deltaR(iMuon2->eta(), iMuon2->phi(), isoMu24_objects[j].eta(), isoMu24_objects[j].phi());

            if (dR_temp < 0.1) {
              dm_matched_sm = true;
              break;
            }
          }
        }
      }
    }
  }

  // if (!firedHLT_L1DoubleMu)
  //   return;

  // 2. Find the triggar objects that has the last filter of HLT_Mu0_L1DoubleMu_v*

  // for (const auto &obj : *triggerObjects) {
  //   pat::TriggerObjectStandAlone unpackedObj = obj;
  //   unpackedObj.unpackPathNames(names);
  //   unpackedObj.unpackFilterLabels(iEvent, *TriggerResults);
  //   if (unpackedObj.hasFilterLabel("hltL3fL1sDoubleMu0SQL1f0L2PreFilteres0L3Filtered0")) {

  //   }
  // }

  // n_dm_matched = 0;
  // if (firedHLT_L1DoubleMu) {
  //   for (pat::MuonCollection::const_iterator iMuon2 = thePATMuonHandle->begin(); iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
  //     for (const auto &obj : *triggerObjects) {
  //       pat::TriggerObjectStandAlone unpackedObj = obj;
  //       unpackedObj.unpackPathNames(names);
  //       unpackedObj.unpackFilterLabels(iEvent, *TriggerResults);
  //       if (unpackedObj.hasFilterLabel("hltL3fL1sDoubleMu0SQL1f0L2PreFilteres0L3Filtered0")) {
  //         float dR_temp = reco::deltaR(iMuon2->eta(), iMuon2->phi(), unpackedObj.eta(), unpackedObj.phi());
  //         // if (dR_temp < 0.01) {
  //           n_dm_matched += 1;
  //           // break;
  //         // }
  //       }
  //     }
  //   }
  // }
  // cout << "Number of muons matched to dimuon trigger objects: " << n_dm_matched << std::endl;

  // // Find the trigger object that matches to this muon

  // // cout << "### Trigger objects: " << triggerObjects->size() << std::endl;
  // int fcount = 0;
  // int pcount = 0;
  // for (const auto &obj : *triggerObjects) {
  //   // cout << "Examining trigger object with pt=" << obj.pt() << ", eta=" << obj.eta() << ", phi=" << obj.phi() << std::endl;
  //   pat::TriggerObjectStandAlone unpackedObj = obj;
  //   unpackedObj.unpackPathNames(names);
  //   unpackedObj.unpackFilterLabels(iEvent, *TriggerResults);

  //   if (unpackedObj.hasFilterLabel("hltL3fL1sDoubleMu0SQL1f0L2PreFilteres0L3Filtered0")) {
  //     // std::cout << "Muon matched to last filter: pt=" << unpackedObj.pt() << ", eta=" << unpackedObj.eta() << ", phi=" << unpackedObj.phi() << std::endl;
  //     fcount++;
  //   }

  //   // true, true = last filter in path
  //   if (unpackedObj.hasPathName("HLT_IsoMu24_v*", true, true)) {
  //     // std::cout << "Found trigger object with pt=" << obj.pt() << ", eta=" << obj.eta() << ", phi=" << obj.phi() << std::endl;
  //     pcount++;
  //     // for (auto &label : unpackedObj.filterLabels()) {
  //     //   std::cout << "  filter label: " << label << std::endl;
  //     // }
  //   }
  // }
  // cout << pcount << " : " << fcount << std::endl;

  tree_->Fill();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------

void tagandprobe::beginJob() {
  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  tree_ = new TTree("ntuple", "ntuple");

  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  // tree_->Branch("Event", &Event);

  // tree_->Branch("B_Mu1_pt", &B_Mu1_pt);
  // tree_->Branch("B_Mu1_eta", &B_Mu1_eta);
  // tree_->Branch("B_Mu2_pt", &B_Mu2_pt);
  // tree_->Branch("B_Mu2_eta", &B_Mu2_eta);
  // tree_->Branch("dR_temp", &dR_temp);

  tree_->Branch("firedHLT_IsoMu24", &firedHLT_IsoMu24);
  tree_->Branch("sm_triggerMatched", &sm_triggerMatched);
  tree_->Branch("nIsoMu24Objs", &nIsoMu24Objs);
  tree_->Branch("nIsoMu24MatchedObjs", &nIsoMu24MatchedObjs);
  tree_->Branch("matchMask", &matchMask);
  tree_->Branch("firedHLT_L1DoubleMu", &firedHLT_L1DoubleMu);
  tree_->Branch("dm_matched_offline", &dm_matched_offline);
  tree_->Branch("dm_matched_sm", &dm_matched_sm);
  // tree_->Branch("hasDimuonInRange", &hasDimuonInRange);
  // tree_->Branch("n_dm_matched", &n_dm_matched);

  // tree_->Branch("dm_triggerMatched_m1", &dm_triggerMatched_m1);
  // tree_->Branch("dm_triggerMatched_m2", &dm_triggerMatched_m2);
}

// ------------ method called once each job just after ending the event loop  ------------
void tagandprobe::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void tagandprobe::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  //  Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(tagandprobe);