// -*- C++ -*-
//
// Package:    Muon/MyAnalyzer
// Class:      MyAnalyzer
//
/**\class MyAnalyzer MyAnalyzer.cc Muon/MyAnalyzer/plugins/MyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nimmitha Karunarathna
//         Created:  Thu, 15 Jun 2023 23:43:35 GMT
//
//

#include "CommonTools/CandAlgos/interface/CandMatcher.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"
#include <TFile.h>
#include "TLorentzVector.h"

// Kinematic vertex fitter
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"

//trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

// packedCandidates
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// From Kalman example
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

#include "DataFormats/PatCandidates/interface/Muon.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;
// edm::EDGetTokenT<std::vector<pat::Muon> > muonCollLabel;

class miniAODhzpkk : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit miniAODhzpkk(const edm::ParameterSet&);
  ~miniAODhzpkk() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<pat::MuonCollection> muonCollLabel;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitLabel;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> estoken_MF;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> estoken_TTB;
  edm::ParameterSet theConfig;

  TTree* tree_;
  float Run, LumiBlock, Event;

  float Z_mass, Z_px, Z_py, Z_pz, Z_pt;
  float Z_eta, Z_phi, Z_vx, Z_vy, Z_vz, Z_rapidity;
  float Z_Vtx_mass, Z_Vtx_Px, Z_Vtx_Py, Z_Vtx_Pz, Z_Vtx_Pt;
  float Z_Vtx_Eta, Z_Vtx_Phi, Z_Vtx_rapidity;
  float Z_Vtx_x, Z_Vtx_y, Z_Vtx_z, Z_Vtx_xError, Z_Vtx_yError, Z_Vtx_zError;
  float Z_Vtx_chi2, Z_Vtx_ndof, Z_Vtx_Prob;
  float muP_pt, muP_eta, muP_phi, muP_charge, muP_trackIso;
  float muM_pt, muM_eta, muM_phi, muM_charge, muM_trackIso;
  float muP_fit_pt, muP_fit_ptError, muP_fit_eta, muP_fit_phi, muP_fit_charge;
  float muM_fit_pt, muM_fit_ptError, muM_fit_eta, muM_fit_phi, muM_fit_charge;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
using namespace reco;
using namespace edm;
using namespace std;

miniAODhzpkk::miniAODhzpkk(const edm::ParameterSet& iConfig)
    // : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
    : muonCollLabel(consumes<pat::MuonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("muons"))),
      triggerBitLabel(consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("trigbits"))),
      estoken_MF(esConsumes()),
      estoken_TTB(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      theConfig(iConfig),

      // initialize the variables for the tree
      tree_(0),
      Run(0),
      LumiBlock(0),
      Event(0),
      Z_mass(0),
      Z_px(0),
      Z_py(0),
      Z_pz(0),
      Z_pt(0),
      Z_eta(0),

      Z_phi(0),
      Z_rapidity(0),

      muP_pt(0),
      muP_eta(0),
      muP_phi(0),
      muP_charge(0),

      muM_pt(0),
      muM_eta(0),
      muM_phi(0),
      muM_charge(0)
 {
  //now do what ever initialization is needed
}

miniAODhzpkk::~miniAODhzpkk() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void miniAODhzpkk::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto& theB = &iSetup.getData(estoken_TTB);

  edm::Handle<pat::MuonCollection> thePATMuonHandle;
  iEvent.getByToken(muonCollLabel, thePATMuonHandle);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitLabel, triggerBits);

  if (!theB) {
    cout << "No TransientTrackBuilder in event!" << endl;
    return;
  }

  if (!thePATMuonHandle.isValid()) {
    cout << "No PAT Muons in event!" << endl;
    return;
  }

  if (!triggerBits.isValid()) {
    cout << "No triggerBits in event!" << endl;
    return;
  }

  // work on the trigger information
  const edm::TriggerNames& names = iEvent.triggerNames(*triggerBits);
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    // cout << i << " " << names.triggerName(i) << endl;
    bool accept = triggerBits->accept(i);
    // cout << "Trigger " << names.triggerName(i) << ", accept = " << accept << endl;
    // if (names.triggerName(i).find("HLT_IsoMu27_v") != std::string::npos) {
    //   cout << "Trigger fired: " << names.triggerName(i) << endl;
    // }
  }

  // length of the muon collection
  cout << "thePATMuonHandle->size() = " << thePATMuonHandle->size() << endl;

  const pat::Muon* iMuonP = nullptr;
  const pat::Muon* iMuonM = nullptr;
  // construct Z candidates
  for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
    for (pat::MuonCollection::const_iterator iMuon2 = iMuon1 + 1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
      if (iMuon1->pt() < 3 || iMuon2->pt() < 3)
        continue;

      if (fabs(iMuon1->eta()) > 2.4 || fabs(iMuon2->eta()) > 2.4)
        continue;

      if (iMuon1->charge() * iMuon2->charge() > 0)
        continue;
    
      if (iMuon1->charge() == 1) {
        iMuonP = &(*iMuon1);
        iMuonM = &(*iMuon2);
      } else {
        iMuonP = &(*iMuon2);
        iMuonM = &(*iMuon1);
      }

      TLorentzVector M_p, M_m;
      TLorentzVector MM;
      float mu_mass = 0.1056583745;  //GeV [PDG]

      // Set M_p as positive muon 4-vector and M_m as negative muon 4-vector
      M_p.SetXYZM(iMuonP->px(), iMuonP->py(), iMuonP->pz(), mu_mass);
      M_m.SetXYZM(iMuonM->px(), iMuonM->py(), iMuonM->pz(), mu_mass);

      // Set MM as the dimuon 4-vector
      MM = M_p + M_m;

      // Fill the tree
      Run = iEvent.id().run();
      LumiBlock = iEvent.luminosityBlock();
      Event = iEvent.id().event();

      Z_mass = MM.M();
      Z_px = MM.Px();
      Z_py = MM.Py();
      Z_pz = MM.Pz();
      Z_pt = MM.Pt();
      Z_eta = MM.Eta();
      Z_phi = MM.Phi();
      Z_rapidity = MM.Rapidity();

      muP_pt = iMuonP->pt();
      muP_eta = iMuonP->eta();
      muP_phi = iMuonP->phi();
      muP_charge = iMuonP->charge();

      muM_pt = iMuonM->pt();
      muM_eta = iMuonM->eta();
      muM_phi = iMuonM->phi();
      muM_charge = iMuonM->charge();

      tree_->Fill();

      // clear the variables
      Run = 0;
      LumiBlock = 0;
      Event = 0;

      Z_mass = 0;
      Z_px = 0;
      Z_py = 0;
      Z_pz = 0;
      Z_pt = 0;
      Z_eta = 0;
      Z_phi = 0;
      Z_rapidity = 0;

      muP_pt = 0;
      muP_eta = 0;
      muP_phi = 0;
      muP_charge = 0;

      muM_pt = 0;
      muM_eta = 0;
      muM_phi = 0;
      muM_charge = 0;

    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void miniAODhzpkk::beginJob() {
  // data tree
  tree_ = new TTree("ntuple", "ntuple with info");
  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  tree_->Branch("Event", &Event);
  tree_->Branch("Z_mass", &Z_mass);
  tree_->Branch("Z_px", &Z_px);
  tree_->Branch("Z_py", &Z_py);
  tree_->Branch("Z_pz", &Z_pz);
  tree_->Branch("Z_pt", &Z_pt);
  tree_->Branch("Z_eta", &Z_eta);
  tree_->Branch("Z_phi", &Z_phi);
  tree_->Branch("Z_rapidity", &Z_rapidity);

  tree_->Branch("muP_pt", &muP_pt);
  tree_->Branch("muP_eta", &muP_eta);
  tree_->Branch("muP_phi", &muP_phi);
  tree_->Branch("muP_charge", &muP_charge);

  tree_->Branch("muM_pt", &muM_pt);
  tree_->Branch("muM_eta", &muM_eta);
  tree_->Branch("muM_phi", &muM_phi);
  tree_->Branch("muM_charge", &muM_charge);
}

// ------------ method called once each job just after ending the event loop  ------------
void miniAODhzpkk::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void miniAODhzpkk::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(miniAODhzpkk);
