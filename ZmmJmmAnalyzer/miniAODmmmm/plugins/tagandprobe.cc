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
      firedHLT_L1DoubleMu(false),
      nIsoMu24_objects(0),
      oneSingleMuObject(false),
      dR_SMO_offM(-1),
      sm_triggerMatched(false),
      SMinAcceptance(false),
      dPt(10),
      nOffline_Muons(0),
      OtherMuinAcceptanceI(false),
      OtherMuinAcceptanceII(false),
      foundJPsi(false),
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

  // Single Muon Trigger
  firedHLT_IsoMu24 = false;
  for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
    std::string name = names.triggerName(i);
    if (name.find("HLT_IsoMu24_v") != std::string::npos && TriggerResults->accept(i)) {
      firedHLT_IsoMu24 = true;
      break;
    }
  }

  // Double Muon Trigger
  firedHLT_L1DoubleMu = false;
  for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
    const std::string name = names.triggerName(i);
    if (name.find("HLT_Mu0_L1DoubleMu_v") != std::string::npos && TriggerResults->accept(i)) {
      firedHLT_L1DoubleMu = true;
      break;
    }
  }

  std::vector<pat::TriggerObjectStandAlone> isoMu24_objects;
  nIsoMu24_objects = 0;

  oneSingleMuObject = false;
  dR_SMO_offM = -1;
  sm_triggerMatched = false;

  SMinAcceptance = false;
  dPt = -2;
  OtherMuinAcceptanceI = false;
  OtherMuinAcceptanceII = false;
  foundJPsi = false;

  nOffline_Muons = 0;
  nDMT_objs = 0;
  DMobjMatchedToSMobj = false;
  dR_dmO_smO = -1.0;
  dm_matched_offline = false;
  dR_dmO_OffM = -1.0;
  dmMatchedIsSM = false;
  const pat::Muon* matchedToSM = nullptr;

  int charge = 0;

  if (firedHLT_IsoMu24) {
    // Find all IsoMu24 trigger objects
    for (const auto& obj : *triggerObjects) {
      pat::TriggerObjectStandAlone unpackedObj = obj;
      unpackedObj.unpackPathNames(names);
      unpackedObj.unpackFilterLabels(iEvent, *TriggerResults);

      if (unpackedObj.hasPathName("HLT_IsoMu24_v*", true, true)) {
        isoMu24_objects.push_back(unpackedObj);
      }
    }

    nIsoMu24_objects = isoMu24_objects.size();

    // Ignore if there are exactly 1 object
    if (nIsoMu24_objects == 1) {
      // cout << "oneSingleMuObject" << endl;
      oneSingleMuObject = true;

      // And the singleMuon object matches to a slimmed muon
      float min_dR_SMO_offM = 10;

      for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
        dR_SMO_offM = reco::deltaR(iMuon1->eta(), iMuon1->phi(), isoMu24_objects[0].eta(), isoMu24_objects[0].phi());

        float temp_dPt = (iMuon1->pt() - isoMu24_objects[0].pt()) / isoMu24_objects[0].pt();
        float gap_0 = abs(temp_dPt);

        if (dR_SMO_offM < 0.1 && (temp_dPt < 1 && temp_dPt > -1) && (gap_0 < abs(dPt))) {
          dPt = temp_dPt;

          sm_triggerMatched = true;

          // save iMuon1 for later use
          matchedToSM = &(*iMuon1);

          if (iMuon1->charge() == 1) {
            charge = 1;
          } else {
            charge = -1;
          }
        }
      }

      if (sm_triggerMatched) {  // the triggered muons is in the acceptance
        if ((matchedToSM->pt() > 3) && (abs(matchedToSM->eta()) < 2.4)) {
          // cout << "SMinAcceptance" << endl;
          SMinAcceptance = true;

          // Check the number of offline muons in the acceptance
          nOffline_Muons = 0;
          for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
            if (&(*iMuon1) == matchedToSM) {
              // cout << "Skipping the same muon" << endl;
              continue;
            }

            if (iMuon1->charge() == charge)
              continue;

            if (iMuon1->pt() < 3.0)
              continue;

            if (abs(iMuon1->eta()) > 2.4)
              continue;

            nOffline_Muons++;
            OtherMuinAcceptanceI = true;
            // cout << "OtherMuinAcceptanceI = True" << endl;

            if (iMuon1->pt() < 25) {
              // cout << "OtherMuinAcceptanceII = True" << endl;
              OtherMuinAcceptanceII = true;

              TLorentzVector M1, M2, MM1;
              float mu_mass = 0.1056583745;  //[PDG mass]
              M1.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M2.SetXYZM(matchedToSM->px(), matchedToSM->py(), matchedToSM->pz(), mu_mass);

              MM1 = M1 + M2;

              if ((MM1.M() > 2) && (MM1.M() < 4))
                foundJPsi = true;
            }
          }
        }
      }
      // cout << "Completed Trigger matched muon in acceptance" << endl;

      // Find dimuon trigger objects
      if (firedHLT_L1DoubleMu == true) {
        // 2. Are there at least two offline muons
        // if (thePATMuonHandle->size() >= 2) { atLeastTwoOfflineAtDMT = true; }

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
          // cout << "Recovering the two muons" << std::endl;
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
            // cout << "Successfully recovered the two muons" << std::endl;
          } else {
            // cout << "Failed to recover the two muons" << std::endl;
          }
        }

        nDMT_objs = DoubleMu_objects.size();

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
      // cout << "Completed DoubleMu Loop" << endl;
    }
    // cout << "Completed objects == 1" << endl;
  }
  // cout << "Filling Tree" << endl;

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
  tree_->Branch("firedHLT_L1DoubleMu", &firedHLT_L1DoubleMu);
  // tree_->Branch("nIsoMu24_objects", &nIsoMu24_objects);
  tree_->Branch("oneSingleMuObject", &oneSingleMuObject);
  // tree_->Branch("dR_SMO_offM", &dR_SMO_offM);
  tree_->Branch("sm_triggerMatched", &sm_triggerMatched);
  tree_->Branch("SMinAcceptance", &SMinAcceptance);
  tree_->Branch("dPt", &dPt);
  tree_->Branch("nOffline_Muons", &nOffline_Muons);
  tree_->Branch("OtherMuinAcceptanceI", &OtherMuinAcceptanceI);
  tree_->Branch("OtherMuinAcceptanceII", &OtherMuinAcceptanceII);
  tree_->Branch("foundJPsi", &foundJPsi);
  tree_->Branch("nDMT_objs", &nDMT_objs);
  tree_->Branch("DMobjMatchedToSMobj", &DMobjMatchedToSMobj);
  // tree_->Branch("dR_dmO_smO", &dR_dmO_smO);
  tree_->Branch("dm_matched_offline", &dm_matched_offline);
  // tree_->Branch("dR_dmO_OffM", &dR_dmO_OffM);
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