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

      firedHLT_IsoMu24(false),
      oneSingleMuObject(false),
      dR_SMO_offM(-1),
      sm_triggerMatched(false),
      firedHLT_L1DoubleMu(false),
      atLeastTwoOfflineAtDMT(false),
      nDMT_objs(0),
      DMobjMatchedToSMobj(false),
      dR_dmO_smO(-1.0),
      dm_matched_offline(false),
      dR_dmO_OffM(-1.0),
      dmMatchedIsSM(false)
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

  oneSingleMuObject = false;
  dR_SMO_offM = -1;
  sm_triggerMatched = false;
  firedHLT_L1DoubleMu = false;
  atLeastTwoOfflineAtDMT = false;
  nDMT_objs = 0;
  DMobjMatchedToSMobj = false;
  dR_dmO_smO = -1.0;
  dm_matched_offline = false;
  dR_dmO_OffM = -1.0;
  dmMatchedIsSM = false;
  const pat::Muon* matchedToSM = nullptr;

  // Ignore if there are more than 2 objects
  if (isoMu24_objects.size() == 1) {
    oneSingleMuObject = true;

    // And that object matches to a slimmed muon
    for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
      dR_SMO_offM = reco::deltaR(iMuon1->eta(), iMuon1->phi(), isoMu24_objects[0].eta(), isoMu24_objects[0].phi());
      if (dR_SMO_offM < 0.1) {
        sm_triggerMatched = true;
        // save iMuon1 for later use
        matchedToSM = &(*iMuon1);
        break;
      }
    }

    // 1. Check if event fired HLT_Mu0_L1DoubleMu_v*
    for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
      std::string name = names.triggerName(i);
      if (name.find("HLT_Mu0_L1DoubleMu_v") != std::string::npos && TriggerResults->accept(i)) {
        firedHLT_L1DoubleMu = true;
        break;
      }
    }

    // Find dimuon trigger objects
    if (firedHLT_L1DoubleMu == true) {
      // 2. Are there at least two offline muons
      if (thePATMuonHandle->size() >= 2) {
        atLeastTwoOfflineAtDMT = true;
      }

      std::vector<pat::TriggerObjectStandAlone> DoubleMu_objects;
      for (const auto& obj : *triggerObjects) {
        pat::TriggerObjectStandAlone unpackedObj = obj;
        unpackedObj.unpackPathNames(names);
        unpackedObj.unpackFilterLabels(iEvent, *TriggerResults);

        if (unpackedObj.hasFilterLabel("hltL3fL1sDoubleMu0SQL1f0L2PreFilteres0L3Filtered0")) {
          DoubleMu_objects.push_back(unpackedObj);
        }
      }

      // if failed to find two objects, try recovering them from a different filter
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

      nDMT_objs = DoubleMu_objects.size();
      // cout << "Number of dimuon trigger objects: " << nDMT_objs << std::endl;

      if (nDMT_objs > 0) {
        // Check if one of these DMT objects matches to the IsoMu24 objects
        for (size_t i = 0; i < DoubleMu_objects.size(); ++i) {
          dR_dmO_smO = reco::deltaR(DoubleMu_objects[i].eta(), DoubleMu_objects[i].phi(), isoMu24_objects[0].eta(), isoMu24_objects[0].phi());
          if (dR_dmO_smO < 0.15) {
            DMobjMatchedToSMobj = true;
            break;
          }
        }

        // Check if one of the DMT objects matches to an offline muon
        for (size_t i = 0; i < DoubleMu_objects.size(); ++i) {
          for (pat::MuonCollection::const_iterator iMuon2 = thePATMuonHandle->begin(); iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
            float temp_dR_dmO_OffM = reco::deltaR(iMuon2->eta(), iMuon2->phi(), DoubleMu_objects[i].eta(), DoubleMu_objects[i].phi());
            if (temp_dR_dmO_OffM < 0.1) {
              dm_matched_offline = true;
              dR_dmO_OffM = temp_dR_dmO_OffM;

              // Check if this muon is the same as the one matched to the IsoMu24 trigger object
              if (&(*iMuon2) == matchedToSM) {
                dmMatchedIsSM = true;
                break;
              }
            }
          }
        }
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

void tagandprobe::beginJob() {
  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  tree_ = new TTree("ntuple", "ntuple");

  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  // tree_->Branch("Event", &Event);

  tree_->Branch("firedHLT_IsoMu24", &firedHLT_IsoMu24);
  tree_->Branch("oneSingleMuObject", &oneSingleMuObject);
  tree_->Branch("dR_SMO_offM", &dR_SMO_offM);
  tree_->Branch("sm_triggerMatched", &sm_triggerMatched);
  tree_->Branch("firedHLT_L1DoubleMu", &firedHLT_L1DoubleMu);
  tree_->Branch("atLeastTwoOfflineAtDMT", &atLeastTwoOfflineAtDMT);
  tree_->Branch("nDMT_objs", &nDMT_objs);
  tree_->Branch("DMobjMatchedToSMobj", &DMobjMatchedToSMobj);
  tree_->Branch("dR_dmO_smO", &dR_dmO_smO);
  tree_->Branch("dm_matched_offline", &dm_matched_offline);
  tree_->Branch("dR_dmO_OffM", &dR_dmO_OffM);
  tree_->Branch("dmMatchedIsSM", &dmMatchedIsSM);
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