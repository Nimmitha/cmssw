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

class MyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MyAnalyzer(const edm::ParameterSet&);
  ~MyAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<pat::MuonCollection> muonCollLabel;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitLabel;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleLabel;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandsLabel;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> estoken_MF;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> estoken_TTB;
  edm::ParameterSet theConfig;

  TTree* genTree_;
  float genZ_mass, genZ_px, genZ_py, genZ_pz, genZ_pt;
  float genZ_eta, genZ_phi;
  float genPhi_mass, genPhi_px, genPhi_py, genPhi_pz, genPhi_pt;
  float genPhi_eta, genPhi_phi;
  float genMuP_pt, genMuP_eta, genMuP_phi, genMuP_charge;
  float genMuM_pt, genMuM_eta, genMuM_phi, genMuM_charge;
  float genKp_pt, genKp_eta, genKp_phi, genKp_charge;
  float genKm_pt, genKm_eta, genKm_phi, genKm_charge;

  TTree* trackTree_;
  float Kp_Px, Kp_Py, Kp_Pz, Kp_Pt, Kp_Eta, Kp_Phi, Kp_charge;
  float Km_Px, Km_Py, Km_Pz, Km_Pt, Km_Eta, Km_Phi, Km_charge;
  float Phi_mass, Phi_Px, Phi_Py, Phi_Pz, Phi_Pt, Phi_Eta, Phi_Phi;

  TTree* tree_;
  float Run, LumiBlock, Event;

  float Z_mass, Z_px, Z_py, Z_pz, Z_pt;
  float Z_eta, Z_phi, Z_vx, Z_vy, Z_vz, Z_rapidity;
  float Z_Vtx_mass, Z_Vtx_Px, Z_Vtx_Py, Z_Vtx_Pz, Z_Vtx_Pt;
  float Z_Vtx_Eta, Z_Vtx_Phi, Z_Vtx_rapidity;
  float Z_Vtx_x, Z_Vtx_y, Z_Vtx_z, Z_Vtx_xError, Z_Vtx_yError, Z_Vtx_zError;
  float Z_Vtx_chi2, Z_Vtx_ndof, Z_Vtx_Prob;
  float muP_pt, muP_eta, muP_phi, muP_charge;
  float muM_pt, muM_eta, muM_phi, muM_charge;
  float muP_fit_pt, muP_fit_ptError, muP_fit_eta, muP_fit_phi, muP_fit_charge;
  float muM_fit_pt, muM_fit_ptError, muM_fit_eta, muM_fit_phi, muM_fit_charge;

  TH1D* mumuMass_;

  // edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  // double trackPtMin_;
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

MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig)
    // : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
    : muonCollLabel(consumes<pat::MuonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("muons"))),
      triggerBitLabel(consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("trigbits"))),
      genParticleLabel(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles"))),
      pfCandsLabel(consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfCands"))),
      // genParticleLabel = consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));
      estoken_MF(esConsumes()),
      estoken_TTB(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      theConfig(iConfig),

      // initialize the genTree
      genTree_(0),
      genZ_mass(0),
      genZ_px(0),
      genZ_py(0),
      genZ_pz(0),
      genZ_pt(0),
      genZ_eta(0),
      genZ_phi(0),
      genPhi_mass(0),
      genPhi_px(0),
      genPhi_py(0),
      genPhi_pz(0),
      genPhi_pt(0),
      genPhi_eta(0),
      genPhi_phi(0),
      genMuP_pt(0),
      genMuP_eta(0),
      genMuP_phi(0),
      genMuP_charge(0),
      genMuM_pt(0),
      genMuM_eta(0),
      genMuM_phi(0),
      genMuM_charge(0),
      genKp_pt(0),
      genKp_eta(0),
      genKp_phi(0),
      genKp_charge(0),
      genKm_pt(0),
      genKm_eta(0),
      genKm_phi(0),
      genKm_charge(0),

      // initialize the pfCandsTree
      trackTree_(0),
      Kp_Px(0),
      Kp_Py(0),
      Kp_Pz(0),
      Kp_Pt(0),
      Kp_Eta(0),
      Kp_Phi(0),
      Kp_charge(0),
      Km_Px(0),
      Km_Py(0),
      Km_Pz(0),
      Km_Pt(0),
      Km_Eta(0),
      Km_Phi(0),
      Km_charge(0),
      Phi_mass(0),
      Phi_Px(0),
      Phi_Py(0),
      Phi_Pz(0),
      Phi_Pt(0),
      Phi_Eta(0),
      Phi_Phi(0),

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
      Z_vx(0),
      Z_vy(0),
      Z_vz(0),
      Z_rapidity(0),
      Z_Vtx_mass(0),
      Z_Vtx_Px(0),
      Z_Vtx_Py(0),
      Z_Vtx_Pz(0),
      Z_Vtx_Pt(0),
      Z_Vtx_Eta(0),
      Z_Vtx_Phi(0),

      Z_Vtx_rapidity(0),

      Z_Vtx_x(0),
      Z_Vtx_y(0),
      Z_Vtx_z(0),

      Z_Vtx_xError(0),
      Z_Vtx_yError(0),
      Z_Vtx_zError(0),

      Z_Vtx_chi2(0),
      Z_Vtx_ndof(0),
      Z_Vtx_Prob(0),

      muP_pt(0),
      muP_eta(0),
      muP_phi(0),
      muP_charge(0),

      muM_pt(0),
      muM_eta(0),
      muM_phi(0),
      muM_charge(0),

      muP_fit_pt(0),
      muP_fit_ptError(0),
      muP_fit_eta(0),
      muP_fit_phi(0),
      muP_fit_charge(0),

      muM_fit_pt(0),
      muM_fit_ptError(0),
      muM_fit_eta(0),
      muM_fit_phi(0),
      muM_fit_charge(0) {
  //now do what ever initialization is needed
  // edm::InputTag muonTag("slimmedMuons");
  // muonCollLabel = consumes<pat::MuonCollection>(muonTag);
  // edm::InputTag tracksTag("generalTracks");
  // token_tracks = consumes<TrackCollection>(tracksTag);

  edm::Service<TFileService> fs;
  mumuMass_ = fs->make<TH1D>("mumuMass", "Z Candidates in Di-Muon Channel", 100, 0, 120);
}

MyAnalyzer::~MyAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void MyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto& theB = &iSetup.getData(estoken_TTB);

  edm::Handle<pat::MuonCollection> thePATMuonHandle;
  iEvent.getByToken(muonCollLabel, thePATMuonHandle);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitLabel, triggerBits);

  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genParticleLabel, pruned);

  edm::Handle<pat::PackedCandidateCollection> pfCands;
  iEvent.getByToken(pfCandsLabel, pfCands);

  // variables for the gen particles
  const reco::Candidate* genHiggs = 0;
  const reco::Candidate* genZ = 0;
  const reco::Candidate* genMuP = 0;
  const reco::Candidate* genMuM = 0;
  const reco::Candidate* genPhi = 0;
  const reco::Candidate* genKp = 0;
  const reco::Candidate* genKm = 0;

  float genZ_candidates = 0;

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

  bool GenInfo = true;
  if (GenInfo == true) {
    // loop over the pruned collection
    for (size_t i = 0; i < pruned->size(); i++) {
      // if not a muon continue
      if ((*pruned)[i].pdgId() != 13)
        continue;

      // if the mother is not a Z continue
      if ((*pruned)[i].mother()->pdgId() != 23)
        continue;

      // if the mother of the mother is not a Higgs continue
      if ((*pruned)[i].mother()->mother()->pdgId() != 35)
        continue;

      // Now we have a muon from a Z from a Higgs
      // Assign Higgs to a variable
      genHiggs = (*pruned)[i].mother()->mother();

      // if the Higgs has more than 2 daughters continue
      if (genHiggs->numberOfDaughters() != 2)
        continue;

      if (genHiggs->daughter(0)->pdgId() == 23 && genHiggs->daughter(1)->pdgId() == 333) {
        genZ = genHiggs->daughter(0);
        genPhi = genHiggs->daughter(1);
      } else if (genHiggs->daughter(0)->pdgId() == 333 && genHiggs->daughter(1)->pdgId() == 23) {
        genZ = genHiggs->daughter(1);
        genPhi = genHiggs->daughter(0);
      } else {
        continue;
      }

      // check if the genZ has muons as daughters
      if (genZ->numberOfDaughters() != 2)
        continue;

      // variables for the daughters of the genZ
      if (genZ->daughter(0)->pdgId() == 13 && genZ->daughter(1)->pdgId() == -13) {
        genMuP = genZ->daughter(1);
        genMuM = genZ->daughter(0);
      } else if (genZ->daughter(0)->pdgId() == -13 && genZ->daughter(1)->pdgId() == 13) {
        genMuP = genZ->daughter(0);
        genMuM = genZ->daughter(1);
      } else {
        continue;
      }

      // check if the phi has kaons as daughters
      if (genPhi->numberOfDaughters() != 2)
        continue;

      if (genPhi->daughter(0)->pdgId() == 321 && genPhi->daughter(1)->pdgId() == -321) {
        genKp = genPhi->daughter(0);
        genKm = genPhi->daughter(1);
      } else if (genPhi->daughter(0)->pdgId() == -321 && genPhi->daughter(1)->pdgId() == 321) {
        genKp = genPhi->daughter(1);
        genKm = genPhi->daughter(0);
      } else {
        continue;
      }

      // print the pdgId of the particles
      // cout << "genHiggs->pdgId() = " << genHiggs->pdgId() << endl;
      // cout << "genZ->pdgId() = " << genZ->pdgId() << endl;
      // cout << "genPhi->pdgId() = " << genPhi->pdgId() << endl;
      // cout << "genMuP->pdgId() = " << genMuP->pdgId() << endl;
      // cout << "genMuM->pdgId() = " << genMuM->pdgId() << endl;
      // cout << "genKp->pdgId() = " << genKp->pdgId() << endl;
      // cout << "genKm->pdgId() = " << genKm->pdgId() << endl;

      genZ_candidates++;

      // Fill the gen tree
      Event = iEvent.id().event();
      genZ_mass = genZ->mass();
      genZ_px = genZ->px();
      genZ_py = genZ->py();
      genZ_pz = genZ->pz();
      genZ_pt = genZ->pt();
      genZ_eta = genZ->eta();
      genZ_phi = genZ->phi();

      genMuP_pt = genMuP->pt();
      genMuP_eta = genMuP->eta();
      genMuP_phi = genMuP->phi();
      genMuP_charge = genMuP->charge();

      genMuM_pt = genMuM->pt();
      genMuM_eta = genMuM->eta();
      genMuM_phi = genMuM->phi();
      genMuM_charge = genMuM->charge();

      genPhi_mass = genPhi->mass();
      genPhi_px = genPhi->px();
      genPhi_py = genPhi->py();
      genPhi_pz = genPhi->pz();
      genPhi_pt = genPhi->pt();
      genPhi_eta = genPhi->eta();
      genPhi_phi = genPhi->phi();

      genKp_pt = genKp->pt();
      genKp_eta = genKp->eta();
      genKp_phi = genKp->phi();
      genKp_charge = genKp->charge();

      genKm_pt = genKm->pt();
      genKm_eta = genKm->eta();
      genKm_phi = genKm->phi();
      genKm_charge = genKm->charge();
      genTree_->Fill();
      // clear gen variables
      Event = 0;
      genZ_mass = 0;
      genZ_px = 0;
      genZ_py = 0;
      genZ_pz = 0;
      genZ_pt = 0;
      genZ_eta = 0;
      genZ_phi = 0;

      genPhi_mass = 0;
      genPhi_px = 0;
      genPhi_py = 0;
      genPhi_pz = 0;
      genPhi_pt = 0;
      genPhi_eta = 0;
      genPhi_phi = 0;

      genMuP_pt = 0;
      genMuP_eta = 0;
      genMuP_phi = 0;
      genMuP_charge = 0;

      genMuM_pt = 0;
      genMuM_eta = 0;
      genMuM_phi = 0;
      genMuM_charge = 0;

      genKp_pt = 0;
      genKp_eta = 0;
      genKp_phi = 0;
      genKp_charge = 0;

      genKm_pt = 0;
      genKm_eta = 0;
      genKm_phi = 0;
      genKm_charge = 0;
    }

    if (genZ_candidates > 1) {
      cout << "More than one Z candidate" << endl;
    }
  }

  bool pfCandsInfo = true;
  if (pfCandsInfo == true){
    // loop over the pfCands collection
    for (pat::PackedCandidateCollection::const_iterator cTrack1 = pfCands->begin(); cTrack1 != pfCands->end(); cTrack1++) {
      if (cTrack1->pt() < 5 || fabs(cTrack1->eta()) > 2.4 || cTrack1->pdgId() != 211)
        continue;

      for (pat::PackedCandidateCollection::const_iterator cTrack2 = cTrack1 + 1; cTrack2 != pfCands->end(); cTrack2++) {
        if (cTrack1->charge() * cTrack2->charge() > 0) {
          continue;
        }
        if (cTrack2->pt() < 5 || fabs(cTrack2->eta()) > 2.4 || cTrack1->pdgId() != 211)
          continue;

        // print track isolation information
        // cout << "cTrack2->trackIso() = " << cTrack2->trackIso() << endl;

        // Set P as the positive Kaon and M as the negative Kaon
        TLorentzVector Kp_Lorentz, Km_Lorentz;
        TLorentzVector Phi;
        float K_mass = 0.493677;  //K+/- GeV [PDG]

        if (cTrack1->charge() == 1){
          Kp_Lorentz.SetXYZM(cTrack1->px(), cTrack1->py(), cTrack1->pz(), K_mass);
          Km_Lorentz.SetXYZM(cTrack2->px(), cTrack2->py(), cTrack2->pz(), K_mass);
        } else {
          Km_Lorentz.SetXYZM(cTrack1->px(), cTrack1->py(), cTrack1->pz(), K_mass);
          Kp_Lorentz.SetXYZM(cTrack2->px(), cTrack2->py(), cTrack2->pz(), K_mass);
        }

        // Set KK as the diKaon 4-vector
        Phi = Kp_Lorentz + Km_Lorentz;

        Event = iEvent.id().event();
        Kp_Px = Kp_Lorentz.Px();
        Kp_Py = Kp_Lorentz.Py();
        Kp_Pz = Kp_Lorentz.Pz();
        Kp_Pt = Kp_Lorentz.Pt();
        Kp_Eta = Kp_Lorentz.Eta();
        Kp_Phi = Kp_Lorentz.Phi();
        Kp_charge = cTrack1->charge();

        Km_Px = Km_Lorentz.Px();
        Km_Py = Km_Lorentz.Py();
        Km_Pz = Km_Lorentz.Pz();
        Km_Pt = Km_Lorentz.Pt();
        Km_Eta = Km_Lorentz.Eta();
        Km_Phi = Km_Lorentz.Phi();
        Km_charge = cTrack2->charge();

        Phi_mass = Phi.M();
        Phi_Px = Phi.Px();
        Phi_Py = Phi.Py();
        Phi_Pz = Phi.Pz();
        Phi_Pt = Phi.Pt();
        Phi_Eta = Phi.Eta();
        Phi_Phi = Phi.Phi();

        trackTree_->Fill();

        // clear the variables
        Event = 0;
        Kp_Px = 0;
        Kp_Py = 0;
        Kp_Pz = 0;
        Kp_Pt = 0;
        Kp_Eta = 0;
        Kp_Phi = 0;
        Kp_charge = 0;

        Km_Px = 0;
        Km_Py = 0;
        Km_Pz = 0;
        Km_Pt = 0;
        Km_Eta = 0;
        Km_Phi = 0;
        Km_charge = 0;

        Phi_mass = 0;
        Phi_Px = 0;
        Phi_Py = 0;
        Phi_Pz = 0;
        Phi_Pt = 0;
        Phi_Eta = 0;
        Phi_Phi = 0;
      }
    }

  //   for (size_t i = 0; i < pfCands->size(); i++) {
  //     // if charge is 0 continue
  //     // if ((*pfCands)[i].charge() == 0)
  //       // continue;

  //     if the abs(pdgId) is not 211
  //     // if (abs((*pfCands)[i].pdgId()) != 211)
  //       // continue;

  //     if pt < 0.5 continue
  //     // if ((*pfCands)[i].pt() < 0.5)
  //       // continue;

  //     if eta > 2.4 continue
  //     // if (abs((*pfCands)[i].eta()) > 2.4)
  //       // continue;

  //     // TLorentzVector TpfCand;
  //     // float K_mass = 0.493677;  //K+/- GeV [PDG]
  //     // TpfCand.SetXYZM((*pfCands)[i].px(), (*pfCands)[i].py(), (*pfCands)[i].pz(), K_mass);

  //     cout << "pfCand pt = " << (*pfCands)[i].pt() << endl;
  //     cout << "TpfCand pt = " << TpfCand.Pt() << endl;

  //     cout << "pfCand eta = " << (*pfCands)[i].eta() << endl;
  //     cout << "TpfCand eta = " << TpfCand.Eta() << endl;

  //     cout << "pfCand phi = " << (*pfCands)[i].phi() << endl;
  //     cout << "TpfCand phi = " << TpfCand.Phi() << endl;

  //     // 
  //     Fill the pfCands tree
  //     // Event = iEvent.id().event();
  //     // pfCands_px = (*pfCands)[i].px();
  //     // pfCands_py = (*pfCands)[i].py();
  //     // pfCands_pz = (*pfCands)[i].pz();
  //     // pfCands_pt = (*pfCands)[i].pt();
  //     // pfCands_eta = (*pfCands)[i].eta();
  //     // pfCands_phi = (*pfCands)[i].phi();
  //     // pfCands_charge = (*pfCands)[i].charge();
  //     // pfCands_vx = (*pfCands)[i].vx();
  //     // pfCands_vy = (*pfCands)[i].vy();
  //     // pfCands_vz = (*pfCands)[i].vz();
  //     // trackTree_->Fill();

  //     clear pfCands variables
  //     // Event = 0;
  //     // pfCands_px = 0;
  //     // pfCands_py = 0;
  //     // pfCands_pz = 0;
  //     // pfCands_pt = 0;
  //     // pfCands_eta = 0;
  //     // pfCands_phi = 0;
  //     // pfCands_charge = 0;
  //     // pfCands_vx = 0;
  //     // pfCands_vy = 0;
  //     // pfCands_vz = 0;
  //   }
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
  // cout << "thePATMuonHandle->size() = " << thePATMuonHandle->size() << endl;

  // construct Z candidates
  for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
    // print track parameters
    reco::TrackRef glbTrackP;
    reco::TrackRef glbTrackM;

    // if (trackRef->quality(reco::TrackBase::highPurity)) {
    //   cout << "Track is high purity" << endl;
    // } else {
    //   cout << "Track is not high purity" << endl;
    // }

    // cout << tr

    // print track quality
    // cout << "Track quality: " << iMuon1->track()->quality(reco::TrackBase::highPurity) << endl;

    // eta cut and pt cut
    if (iMuon1->pt() < 3 || fabs(iMuon1->eta()) > 2.4)
      continue;
    for (pat::MuonCollection::const_iterator iMuon2 = iMuon1 + 1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
      if (iMuon1->charge() * iMuon2->charge() > 0) {
        continue;
      }
      // eta cut
      if (iMuon2->pt() < 3 || fabs(iMuon2->eta()) > 2.4)
        continue;
      mumuMass_->Fill((iMuon1->p4() + iMuon2->p4()).mass());

      // Set P as positive muon and M as negative muon
      if (iMuon1->charge() == 1) {
        glbTrackP = iMuon1->track();
        glbTrackM = iMuon2->track();
      } else {
        glbTrackP = iMuon2->track();
        glbTrackM = iMuon1->track();
      }

      // check if track is available
      if (glbTrackP.isNull() || glbTrackM.isNull()) {
        // cout << "Continue due ot no track ref" << endl;
        continue;
      }

      TLorentzVector M_p, M_m;
      TLorentzVector MM;
      float mu_mass = 0.1056583745;  //GeV [PDG]

      // Set M_p as positive muon 4-vector and M_m as negative muon 4-vector
      if (iMuon1->charge() == 1) {
        M_p.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
        M_m.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
      } else {
        M_m.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
        M_p.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
      }

      // Set MM as the dimuon 4-vector
      MM = M_p + M_m;

      reco::TransientTrack glbTrackP_TT = theB->build(glbTrackP);
      reco::TransientTrack glbTrackM_TT = theB->build(glbTrackM);

      std::vector<reco::TransientTrack> mu_tks;
      KalmanVertexFitter kvf(true);

      // Transform Track to TransientTrack
      mu_tks.clear();
      mu_tks.push_back(glbTrackP_TT);
      mu_tks.push_back(glbTrackM_TT);

      TransientVertex Z_candi = kvf.vertex(mu_tks);
      // fit parameters
      // cout << "TESTING REFITTED TRACKS" << Z_candi.refittedTrack(mu_tks[0]) << endl;

      if (!Z_candi.isValid()) {
        // cout << "Continue due to invalid vertex" << endl;
        continue;
      }

      reco::Vertex Z_vertex = Z_candi;

      float Z_Prob = TMath::Prob(Z_candi.totalChiSquared(), Z_candi.degreesOfFreedom());
      const math::XYZTLorentzVector Z_mom = Z_vertex.p4(mu_mass, 0.0);

      // cout << "Z total chi squared = " << Z_candi.totalChiSquared() << endl;
      // cout << "Z degrees of freedom = " << Z_candi.degreesOfFreedom() << endl;
      // cout << "Z probability = " << Z_Prob << endl;

      // cout << "Z mass before fit = " << MM.M() << endl;
      // cout << "Z mass after fit = " << Z_mom.mass() << endl;

      // cout << "Z pt before fit = " << MM.Pt() << endl;
      // cout << "Z pt after fit = " << Z_mom.Pt() << endl;

      // Get the refitted tracks
      reco::TransientTrack fit_glbTrackP_TT = Z_candi.refittedTrack(glbTrackP_TT);
      reco::TransientTrack fit_glbTrackM_TT = Z_candi.refittedTrack(glbTrackM_TT);

      // Check if refitted tracks are valid
      if (!fit_glbTrackP_TT.isValid() || !fit_glbTrackM_TT.isValid()) {
        cout << "Continue due to invalid refitted tracks" << endl;
        continue;
      }

      // cout << "Mu1 pt before = " << glbTrackP_TT.track().pt() << endl;
      // cout << "Mu1 pt after = " << fit_glbTrackP_TT.track().pt() << endl;

      // // errors
      // cout << "Mu1 pt Error before = " << glbTrackP_TT.track().ptError() << endl;
      // cout << "Mu1 pt Error after = " << fit_glbTrackP_TT.track().ptError() << endl;

      // cout << "Mu2 pt before = " << glbTrackM_TT.track().pt() << endl;
      // cout << "Mu2 pt after = " << fit_glbTrackM_TT.track().pt() << endl;

      // //errors
      // cout << "Mu2 pt Error before = " << glbTrackM_TT.track().ptError() << endl;
      // cout << "Mu2 pt Error after = " << fit_glbTrackM_TT.track().ptError() << endl;

      // cout << "Mu1 eta before = " << glbTrackP_TT.track().eta() << endl;
      // cout << "Mu1 eta after = " << fit_glbTrackP_TT.track().eta() << endl;

      // cout << "Mu2 eta before = " << glbTrackM_TT.track().eta() << endl;
      // cout << "Mu2 eta after = " << fit_glbTrackM_TT.track().eta() << endl;

      // cout << "Mu1 phi before = " << glbTrackP_TT.track().phi() << endl;
      // cout << "Mu1 phi after = " << fit_glbTrackP_TT.track().phi() << endl;

      // cout << "Mu2 phi before = " << glbTrackM_TT.track().phi() << endl;
      // cout << "Mu2 phi after = " << fit_glbTrackM_TT.track().phi() << endl;

      // cout << "Mu1 charge before = " << glbTrackP_TT.track().charge() << endl;
      // cout << "Mu1 charge after = " << fit_glbTrackP_TT.track().charge() << endl;

      // cout << "Mu2 charge before = " << glbTrackM_TT.track().charge() << endl;
      // cout << "Mu2 charge after = " << fit_glbTrackM_TT.track().charge() << endl;

      // cout << "Mu1 chi after = " << fit_glbTrackP_TT.chi2() << endl;
      // cout << "Mu2 chi after = " << fit_glbTrackM_TT.chi2() << endl;

      bool KinFit = false;

      if (KinFit) {
        ParticleMass muon_mass = 0.10565837;  // pdg mass
        float muon_sigma = muon_mass * 1.e-6;
        vector<RefCountedKinematicParticle> muonParticles;

        // Creating a KinematicParticleFactory
        KinematicParticleFactoryFromTransientTrack pFactory;

        // initial chi2 and ndf before kinematic fits.
        float chi = 0.;
        float ndf = 0.;
        muonParticles.push_back(pFactory.particle(glbTrackP_TT, muon_mass, chi, ndf, muon_sigma));
        muonParticles.push_back(pFactory.particle(glbTrackM_TT, muon_mass, chi, ndf, muon_sigma));

        KinematicParticleVertexFitter fitter;
        RefCountedKinematicTree ZTree;

        try {
          ZTree = fitter.fit(muonParticles);
        } catch (...) {
          cout << "ZTree failed" << endl;
          continue;
        }

        if (!ZTree->isValid()) {
          cout << "ZTree is not valid" << endl;
          continue;
        }

        ZTree->movePointerToTheTop();

        RefCountedKinematicParticle Z_vfit = ZTree->currentParticle();
        RefCountedKinematicVertex Z_vfit_vertex = ZTree->currentDecayVertex();

        // refitted dimuon pt
        // float Z_vfit_pt = Z_vfit->currentState().kinematicParameters().momentum().perp();

        // cout << "Z_vfit_pt = " << Z_vfit_pt << endl;
        // cout << "Z pt before (using previous state) = " << Z_vfit->initialState().kinematicParameters().momentum().perp() << endl;

        // // // refitted muon eta
        // // float Z_vfit_eta = Z_vfit->currentState().kinematicParameters().momentum().eta();

        // // cout << "Z_vfit_eta = " << Z_vfit_eta << endl;

        // cout << "Z_vfit_vertex->chiSquared() = " << Z_vfit_vertex->chiSquared() << endl;
        // cout << "Z_vfit_vertex->degreesOfFreedom() = " << Z_vfit_vertex->degreesOfFreedom() << endl;

        // // Get Z mass from the tree
        // float Z_vfit_mass = Z_vfit->currentState().mass();
        // // float Z_vfit_massErr = sqrt(Z_vfit->currentState().kinematicParametersError().matrix()(6, 6));

        // cout << "Z_vfit_mass = " << Z_vfit_mass << endl;
        // // cout << "Z_vfit_massErr = " << Z_vfit_massErr << endl;

        // // Get Z momentum from the tree
        // GlobalVector Z_vfit_momentum2 = Z_vfit->currentState().globalMomentum();
        // float Z_vfit_pt2 = Z_vfit_momentum2.perp();

        // cout << "Z_vfit_pt = " << Z_vfit_pt2 << endl;

        // // probability of the fit
        // float Z_vfit_prob = TMath::Prob(Z_vfit_vertex->chiSquared(), Z_vfit_vertex->degreesOfFreedom());

        // cout << "Z_vfit_prob = " << Z_vfit_prob << endl;

        // Get Z vertex from the tree
        // GlobalPoint Z_vfit_vertex_pos = Z_vfit_vertex->position();
        // float Z_vfit_vertex_x = Z_vfit_vertex_pos.x();
        // float Z_vfit_vertex_y = Z_vfit_vertex_pos.y();
        // float Z_vfit_vertex_z = Z_vfit_vertex_pos.z();

        // cout << "Z_vfit_vertex_x = " << Z_vfit_vertex_x << endl;
        // cout << "Z_vfit_vertex_y = " << Z_vfit_vertex_y << endl;
        // cout << "Z_vfit_vertex_z = " << Z_vfit_vertex_z << endl;

        // Get Z vertex chi2 from the tree
        // float Z_vfit_vertex_chi2 = Z_vfit_vertex->chiSquared();
        // float Z_vfit_vertex_ndof = Z_vfit_vertex->degreesOfFreedom();

        // cout << "Z_vfit_vertex_chi2 = " << Z_vfit_vertex_chi2 << endl;

        // // Get Z vertex normalized chi2 from the tree
        // float Z_vfit_vertex_normChi2 = Z_vfit_vertex_chi2 / Z_vfit_vertex_ndof;

        // cout << "Z_vfit_vertex_normChi2 = " << Z_vfit_vertex_normChi2 << endl;

        // move to a child node
        // bool child = ZTree->movePointerToTheFirstChild();

        // // Get first muon from the tree
        // if (!child){
        //   cout << "ZTree->movePointerToTheFirstChild() failed" << endl;
        //   continue;
        // }

        // RefCountedKinematicParticle fitMu1 = ZTree->currentParticle();
        // float fitMu1_pt = fitMu1->currentState().kinematicParameters().momentum().perp();
        // float fitMu1_eta = fitMu1->currentState().kinematicParameters().momentum().eta();
        // float fitMu1_phi = fitMu1->currentState().kinematicParameters().momentum().phi();
        // float fitMu1_mass = fitMu1->currentState().mass();

        // cout << "before fit Mu1" << endl;
        // cout << "Mu1_charge = " << iMuon1->charge() << endl;
        // cout << "Mu1_pt = " << iMuon1->pt() << endl;
        // cout << "Mu1_eta = " << iMuon1->eta() << endl;
        // cout << "Mu1_phi = " << iMuon1->phi() << endl;
        // cout << "Mu1_mass = " << mu_mass << endl;
        // cout << "Mu1 chi2 = " << iMuon1->track()->chi2() << endl;

        // cout << "before fit Mu2" << endl;
        // cout << "Mu2_charge = " << iMuon2->charge() << endl;
        // cout << "Mu2_pt = " << iMuon2->pt() << endl;
        // cout << "Mu2_eta = " << iMuon2->eta() << endl;
        // cout << "Mu2_phi = " << iMuon2->phi() << endl;
        // cout << "Mu2_mass = " << mu_mass << endl;
        // cout << "Mu2 chi2 = " << iMuon2->track()->chi2() << endl;

        // cout << "after fit Mu1" << endl;
        // cout << "fitMu1_charge = " << fitMu1->currentState().particleCharge() << endl;
        // cout << "fitMu1_pt = " << fitMu1_pt << endl;
        // cout << "fitMu1_eta = " << fitMu1_eta << endl;
        // cout << "fitMu1_phi = " << fitMu1_phi << endl;
        // cout << "fitMu1_mass = " << fitMu1_mass << endl;
        // cout << "fitMu1 chi2 = " << fitMu1->currentState().kinematicParametersError().matrix()(6, 6) << endl;
        // fitMu1->currentState()..kinematicParameters().chiSquared();

        // // move to a child node
        // child = ZTree->movePointerToTheNextChild();

        // if (!child){
        //   cout << "ZTree->movePointerToTheNextChild() failed" << endl;
        //   continue;
        // }

        // cout << "after fit Mu2" << endl;
        // RefCountedKinematicParticle fitMu2 = ZTree->currentParticle();
        // float fitMu2_pt = fitMu2->currentState().kinematicParameters().momentum().perp();
        // float fitMu2_eta = fitMu2->currentState().kinematicParameters().momentum().eta();
        // float fitMu2_phi = fitMu2->currentState().kinematicParameters().momentum().phi();
        // float fitMu2_mass = fitMu2->currentState().mass();

        // cout << "fitMu2_charge = " << fitMu2->currentState().particleCharge() << endl;
        // cout << "fitMu2_pt = " << fitMu2_pt << endl;
        // cout << "fitMu2_eta = " << fitMu2_eta << endl;
        // cout << "fitMu2_phi = " << fitMu2_phi << endl;
        // cout << "fitMu2_mass = " << fitMu2_mass << endl;
        // cout << "fitMu2 chi2 = " << fitMu2->currentState().kinematicParametersError().matrix()(6, 6) << endl;
      }

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

      Z_Vtx_mass = Z_mom.mass();
      Z_Vtx_Px = Z_mom.Px();
      Z_Vtx_Py = Z_mom.Py();
      Z_Vtx_Pz = Z_mom.Pz();
      Z_Vtx_Pt = Z_mom.Pt();
      Z_Vtx_Eta = Z_mom.Eta();
      Z_Vtx_Phi = Z_mom.Phi();
      Z_Vtx_rapidity = Z_mom.Rapidity();

      Z_Vtx_x = Z_vertex.x();
      Z_Vtx_y = Z_vertex.y();
      Z_Vtx_z = Z_vertex.z();
      Z_Vtx_xError = Z_vertex.xError();
      Z_Vtx_yError = Z_vertex.yError();
      Z_Vtx_zError = Z_vertex.zError();

      Z_Vtx_chi2 = Z_candi.totalChiSquared();
      Z_Vtx_ndof = Z_candi.degreesOfFreedom();
      Z_Vtx_Prob = Z_Prob;

      muP_pt = glbTrackP->pt();
      muP_eta = glbTrackP->eta();
      muP_phi = glbTrackP->phi();
      muP_charge = glbTrackP->charge();

      muM_pt = glbTrackM->pt();
      muM_eta = glbTrackM->eta();
      muM_phi = glbTrackM->phi();
      muM_charge = glbTrackM->charge();

      muP_fit_pt = fit_glbTrackP_TT.track().pt();
      muP_fit_ptError = fit_glbTrackP_TT.track().ptError();
      muP_fit_eta = fit_glbTrackP_TT.track().eta();
      muP_fit_phi = fit_glbTrackP_TT.track().phi();
      muP_fit_charge = fit_glbTrackP_TT.track().charge();

      muM_fit_pt = fit_glbTrackM_TT.track().pt();
      muM_fit_ptError = fit_glbTrackM_TT.track().ptError();
      muM_fit_eta = fit_glbTrackM_TT.track().eta();
      muM_fit_phi = fit_glbTrackM_TT.track().phi();
      muM_fit_charge = fit_glbTrackM_TT.track().charge();

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

      Z_Vtx_mass = 0;
      Z_Vtx_Px = 0;
      Z_Vtx_Py = 0;
      Z_Vtx_Pz = 0;
      Z_Vtx_Pt = 0;
      Z_Vtx_Eta = 0;
      Z_Vtx_Phi = 0;
      Z_Vtx_rapidity = 0;

      Z_Vtx_x = 0;
      Z_Vtx_y = 0;
      Z_Vtx_z = 0;
      Z_Vtx_xError = 0;
      Z_Vtx_yError = 0;
      Z_Vtx_zError = 0;

      Z_Vtx_chi2 = 0;
      Z_Vtx_ndof = 0;
      Z_Vtx_Prob = 0;

      muP_pt = 0;
      muP_eta = 0;
      muP_phi = 0;
      muP_charge = 0;

      muM_pt = 0;
      muM_eta = 0;
      muM_phi = 0;
      muM_charge = 0;

      muP_fit_pt = 0;
      muP_fit_ptError = 0;
      muP_fit_eta = 0;
      muP_fit_phi = 0;
      muP_fit_charge = 0;

      muM_fit_pt = 0;
      muM_fit_ptError = 0;
      muM_fit_eta = 0;
      muM_fit_phi = 0;
      muM_fit_charge = 0;
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void MyAnalyzer::beginJob() {
  // gen tree
  genTree_ = new TTree("genTree", "genTree with info");
  genTree_->Branch("Event", &Event);
  genTree_->Branch("genZ_mass", &genZ_mass);
  genTree_->Branch("genZ_px", &genZ_px);
  genTree_->Branch("genZ_py", &genZ_py);
  genTree_->Branch("genZ_pz", &genZ_pz);
  genTree_->Branch("genZ_pt", &genZ_pt);
  genTree_->Branch("genZ_eta", &genZ_eta);
  genTree_->Branch("genZ_phi", &genZ_phi);

  genTree_->Branch("genPhi_mass", &genPhi_mass);
  genTree_->Branch("genPhi_px", &genPhi_px);
  genTree_->Branch("genPhi_py", &genPhi_py);
  genTree_->Branch("genPhi_pz", &genPhi_pz);
  genTree_->Branch("genPhi_pt", &genPhi_pt);
  genTree_->Branch("genPhi_eta", &genPhi_eta);
  genTree_->Branch("genPhi_phi", &genPhi_phi);

  genTree_->Branch("genMuP_pt", &genMuP_pt);
  genTree_->Branch("genMuP_eta", &genMuP_eta);
  genTree_->Branch("genMuP_phi", &genMuP_phi);
  genTree_->Branch("genMuP_charge", &genMuP_charge);

  genTree_->Branch("genMuM_pt", &genMuM_pt);
  genTree_->Branch("genMuM_eta", &genMuM_eta);
  genTree_->Branch("genMuM_phi", &genMuM_phi);
  genTree_->Branch("genMuM_charge", &genMuM_charge);

  genTree_->Branch("genKp_pt", &genKp_pt);
  genTree_->Branch("genKp_eta", &genKp_eta);
  genTree_->Branch("genKp_phi", &genKp_phi);
  genTree_->Branch("genKp_charge", &genKp_charge);

  genTree_->Branch("genKm_pt", &genKm_pt);
  genTree_->Branch("genKm_eta", &genKm_eta);
  genTree_->Branch("genKm_phi", &genKm_phi);
  genTree_->Branch("genKm_charge", &genKm_charge);

  // pfCands tree
  trackTree_ = new TTree("trackTree", "trackTree with info");
  trackTree_->Branch("Event", &Event);
  trackTree_->Branch("Kp_Px", &Kp_Px);
  trackTree_->Branch("Kp_Py", &Kp_Py);
  trackTree_->Branch("Kp_Pz", &Kp_Pz);
  trackTree_->Branch("Kp_Pt", &Kp_Pt);
  trackTree_->Branch("Kp_Eta", &Kp_Eta);
  trackTree_->Branch("Kp_Phi", &Kp_Phi);
  trackTree_->Branch("Kp_charge", &Kp_charge);

  trackTree_->Branch("Km_Px", &Km_Px);
  trackTree_->Branch("Km_Py", &Km_Py);
  trackTree_->Branch("Km_Pz", &Km_Pz);
  trackTree_->Branch("Km_Pt", &Km_Pt);
  trackTree_->Branch("Km_Eta", &Km_Eta);
  trackTree_->Branch("Km_Phi", &Km_Phi);
  trackTree_->Branch("Km_charge", &Km_charge);

  trackTree_->Branch("Phi_mass", &Phi_mass);
  trackTree_->Branch("Phi_Px", &Phi_Px);
  trackTree_->Branch("Phi_Py", &Phi_Py);
  trackTree_->Branch("Phi_Pz", &Phi_Pz);
  trackTree_->Branch("Phi_Pt", &Phi_Pt);
  trackTree_->Branch("Phi_Eta", &Phi_Eta);
  trackTree_->Branch("Phi_Phi", &Phi_Phi);

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

  tree_->Branch("Z_Vtx_mass", &Z_Vtx_mass);
  tree_->Branch("Z_Vtx_Px", &Z_Vtx_Px);
  tree_->Branch("Z_Vtx_Py", &Z_Vtx_Py);
  tree_->Branch("Z_Vtx_Pz", &Z_Vtx_Pz);
  tree_->Branch("Z_Vtx_Pt", &Z_Vtx_Pt);
  tree_->Branch("Z_Vtx_Eta", &Z_Vtx_Eta);
  tree_->Branch("Z_Vtx_Phi", &Z_Vtx_Phi);
  tree_->Branch("Z_Vtx_rapidity", &Z_Vtx_rapidity);

  tree_->Branch("Z_Vtx_x", &Z_Vtx_x);
  tree_->Branch("Z_Vtx_y", &Z_Vtx_y);
  tree_->Branch("Z_Vtx_z", &Z_Vtx_z);
  tree_->Branch("Z_Vtx_xError", &Z_Vtx_xError);
  tree_->Branch("Z_Vtx_yError", &Z_Vtx_yError);
  tree_->Branch("Z_Vtx_zError", &Z_Vtx_zError);

  tree_->Branch("Z_Vtx_chi2", &Z_Vtx_chi2);
  tree_->Branch("Z_Vtx_ndof", &Z_Vtx_ndof);
  tree_->Branch("Z_Vtx_Prob", &Z_Vtx_Prob);

  tree_->Branch("muP_pt", &muP_pt);
  tree_->Branch("muP_eta", &muP_eta);
  tree_->Branch("muP_phi", &muP_phi);
  tree_->Branch("muP_charge", &muP_charge);

  tree_->Branch("muM_pt", &muM_pt);
  tree_->Branch("muM_eta", &muM_eta);
  tree_->Branch("muM_phi", &muM_phi);
  tree_->Branch("muM_charge", &muM_charge);

  tree_->Branch("muP_fit_pt", &muP_fit_pt);
  tree_->Branch("muP_fit_ptError", &muP_fit_ptError);
  tree_->Branch("muP_fit_eta", &muP_fit_eta);
  tree_->Branch("muP_fit_phi", &muP_fit_phi);
  tree_->Branch("muP_fit_charge", &muP_fit_charge);

  tree_->Branch("muM_fit_pt", &muM_fit_pt);
  tree_->Branch("muM_fit_ptError", &muM_fit_ptError);
  tree_->Branch("muM_fit_eta", &muM_fit_eta);
  tree_->Branch("muM_fit_phi", &muM_fit_phi);
  tree_->Branch("muM_fit_charge", &muM_fit_charge);
}

// ------------ method called once each job just after ending the event loop  ------------
void MyAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(MyAnalyzer);
