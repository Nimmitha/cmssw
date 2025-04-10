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
      savedtriggerNames(0),
      savedtriggerBits(0),
      savedtriggerPrescales(0),

      Run(0),
      LumiBlock(0),
      Event(0),

      B_J1_mass(0),
      B_J1_px(0),
      B_J1_py(0),
      B_J1_pz(0),
      B_J1_pt(0),
      B_J1_eta(0),
      B_J1_phi(0),
      B_J1_rapidity(0),

      B_J1_VtxPx(0),
      B_J1_VtxPy(0),
      B_J1_VtxPz(0),
      B_J1_VtxPt(0),
      B_J1_VtxEta(0),
      B_J1_VtxPhi(0),
      B_J1_VtxRapidity(0),
      B_J1_VtxMass(0),
      B_J1_PVx(0),
      B_J1_PVy(0),
      B_J1_PVz(0),
      B_J1_PVxError(0),
      B_J1_PVyError(0),
      B_J1_PVzError(0),

      B_J1_VtxProb(0),

      B_Mu1_px(0),
      B_Mu1_py(0),
      B_Mu1_pz(0),
      B_Mu1_pt(0),
      B_Mu1_eta(0),
      B_Mu1_phi(0),
      B_Mu1_soft(0),
      B_Mu1_tight(0),
      B_Mu1_loose(0),
      B_Mu1_charge(0),

      B_Mu2_px(0),
      B_Mu2_py(0),
      B_Mu2_pz(0),
      B_Mu2_pt(0),
      B_Mu2_eta(0),
      B_Mu2_phi(0),
      B_Mu2_soft(0),
      B_Mu2_tight(0),
      B_Mu2_loose(0),
      B_Mu2_charge(0),
      
      B_M1_pt(0),
      B_M1_eta(0),
      B_M1_phi(0),
      B_M1_px(0),
      B_M1_py(0),
      B_M1_pz(0),
  
      B_M2_pt(0),
      B_M2_eta(0),
      B_M2_phi(0),
      B_M2_px(0),
      B_M2_py(0),
      B_M2_pz(0),

      nB(0)
{
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

  // Only save info if the trigger fired.
  if (MuonTriggerString != "pass") {
    bool triggerFired = false;
    // Loop over all triggers to see if the desired one fired
    for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
      if (names.triggerName(i).find(MuonTriggerString.c_str()) != string::npos) {
        if (TriggerResults->accept(i)) {
          triggerFired = true;  // Trigger fired!
        }
        break;  // No need to check further if we already found the trigger
      }
    }
    if (!triggerFired) {
      return;
    }
  }

  //*********************************
  //Now we get the primary vertex
  //*********************************

  reco::Vertex bestVtx;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  bestVtx = *(primaryVertices_handle->begin());

  //****************************************************
  //*********Now we get the muons***********************
  //****************************************************

  for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
    for (pat::MuonCollection::const_iterator iMuon2 = iMuon1 + 1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
      //make sure all muons are diferent
      if (iMuon1 == iMuon2)
        continue;

      //opposite charge
      //only look for neutral charge combination

      if (!(abs(iMuon1->charge()) == 1))
        continue;
      if (!(abs(iMuon2->charge()) == 1))
        continue;

      if (!((iMuon1->charge()) + (iMuon2->charge())) == 0)
        continue;
      if (iMuon1->pt() < 3.0)
        continue;
      if (iMuon2->pt() < 3.0)
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
      //initialize 4 lepton mass
      float mu_mass = 0.1056583745;  //[PDG mass]
      //float ele_mass =  0.000510998928;//PDG mass

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
      //now M1 is first positive muon 4 vector
      //then M2 is first negative muon 4 vector

      //****************************************************************************************
      //make netral dimuon combination
      MM1 = M1 + M2;
      //****************************************************************************************

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
      //First neutral combination
      mu_tks.clear();
      mu_tks.push_back(muon1TT);
      mu_tks.push_back(muon2TT);
      TransientVertex J_candi1 = kvfM.vertex(mu_tks);

      if (!J_candi1.isValid()) {
        cout << "continue because no vertexed dimuon" << endl;
        continue;
      }

      reco::Vertex JPsi_Vtx1 = J_candi1;

      float B_Prob_tmp1 = TMath::Prob(J_candi1.totalChiSquared(), J_candi1.degreesOfFreedom());
      const math::XYZTLorentzVectorD JPsi_mom1 = JPsi_Vtx1.p4(mu_mass, 0.0);

      if (B_Prob_tmp1 < 0.001) {
        continue;
      }

      // Remove events with mass outside J/Psi or Z mass window
      if ((MM1.M() < 2.6 || MM1.M() > 3.6)) {
        continue;
      }

      //****************************************************************************************
      //Event Information
      Run = iEvent.id().run();
      LumiBlock = iEvent.luminosityBlock();
      Event = iEvent.id().event();

      for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
        savedtriggerNames.push_back(names.triggerName(i));
        savedtriggerBits->push_back(TriggerResults->accept(i));
        savedtriggerPrescales->push_back(triggerPrescales->getPrescaleForIndex<double>(i));
      }

      B_J1_mass = MM1.M();
      B_J1_px = MM1.Px();
      B_J1_py = MM1.Py();
      B_J1_pz = MM1.Pz();
      B_J1_pt = MM1.Pt();
      B_J1_eta = MM1.Eta();
      B_J1_phi = MM1.Phi();
      B_J1_rapidity = MM1.Rapidity();

      B_J1_VtxPx = JPsi_mom1.Px();
      B_J1_VtxPy = JPsi_mom1.Py();
      B_J1_VtxPz = JPsi_mom1.Pz();
      B_J1_VtxPt = JPsi_mom1.Pt();
      B_J1_VtxEta = JPsi_mom1.Eta();
      B_J1_VtxPhi = JPsi_mom1.Phi();
      B_J1_VtxRapidity = JPsi_mom1.Rapidity();
      B_J1_VtxMass = JPsi_mom1.mass();
      B_J1_PVx = JPsi_Vtx1.x();
      B_J1_PVy = JPsi_Vtx1.y();
      B_J1_PVz = JPsi_Vtx1.z();
      B_J1_PVxError = JPsi_Vtx1.xError();
      B_J1_PVyError = JPsi_Vtx1.yError();
      B_J1_PVzError = JPsi_Vtx1.zError();

      //dimuon vtx prob
      B_J1_VtxProb = B_Prob_tmp1;

      //new branch defn for muons
      B_Mu1_px = iMuon1->px();
      B_Mu1_py = iMuon1->py();
      B_Mu1_pz = iMuon1->pz();
      B_Mu1_pt = iMuon1->pt();
      B_Mu1_eta = iMuon1->eta();
      B_Mu1_phi = iMuon1->phi();
      B_Mu1_soft = iMuon1->isSoftMuon(bestVtx);
      B_Mu1_tight = iMuon1->isTightMuon(bestVtx);
      B_Mu1_loose = muon::isLooseMuon(*iMuon1);
      B_Mu1_charge = iMuon1->charge();

      B_Mu2_px = iMuon2->px();
      B_Mu2_py = iMuon2->py();
      B_Mu2_pz = iMuon2->pz();
      B_Mu2_pt = iMuon2->pt();
      B_Mu2_eta = iMuon2->eta();
      B_Mu2_phi = iMuon2->phi();
      B_Mu2_soft = iMuon2->isSoftMuon(bestVtx);
      B_Mu2_tight = iMuon2->isTightMuon(bestVtx);
      B_Mu2_loose = muon::isLooseMuon(*iMuon2);
      B_Mu2_charge = iMuon2->charge();

      B_M1_pt = M1.Pt();
      B_M1_eta = M1.Eta();
      B_M1_phi = M1.Phi();
      B_M1_px = M1.Px();
      B_M1_py = M1.Py();
      B_M1_pz = M1.Pz();

      B_M2_pt = M2.Pt();
      B_M2_eta = M2.Eta();
      B_M2_phi = M2.Phi();
      B_M2_px = M2.Px();
      B_M2_py = M2.Py();
      B_M2_pz = M2.Pz();

      tree_->Fill();

      nB = 0;
      savedtriggerNames.clear();
      savedtriggerBits->clear();
      savedtriggerPrescales->clear();

      Run = -999;
      LumiBlock = -999;
      Event = -999;

      B_J1_mass = -999;
      B_J1_px = -999;
      B_J1_py = -999;
      B_J1_pz = -999;
      B_J1_pt = -999;
      B_J1_eta = -999;
      B_J1_phi = -999;
      B_J1_rapidity = -999;

      B_J1_VtxPx = -999;
      B_J1_VtxPy = -999;
      B_J1_VtxPz = -999;
      B_J1_VtxPt = -999;
      B_J1_VtxEta = -999;
      B_J1_VtxPhi = -999;
      B_J1_VtxRapidity = -999;
      B_J1_VtxMass = -999;
      B_J1_PVx = -999;
      B_J1_PVy = -999;
      B_J1_PVz = -999;
      B_J1_PVxError = -999;
      B_J1_PVyError = -999;
      B_J1_PVzError = -999;
      
      B_J1_VtxProb = -999;

      B_Mu1_px = -999;
      B_Mu1_py = -999;
      B_Mu1_pz = -999;
      B_Mu1_pt = -999;
      B_Mu1_eta = -999;
      B_Mu1_phi = -999;
      B_Mu1_soft = -999;
      B_Mu1_tight = -999;
      B_Mu1_loose = -999;
      B_Mu1_charge = -999;
      
      B_Mu2_px = -999;
      B_Mu2_py = -999;
      B_Mu2_pz = -999;
      B_Mu2_pt = -999;
      B_Mu2_eta = -999;
      B_Mu2_phi = -999;
      B_Mu2_soft = -999;
      B_Mu2_tight = -999;
      B_Mu2_loose = -999;
      B_Mu2_charge = -999;

      B_M1_pt = -999;
      B_M1_eta = -999;
      B_M1_phi = -999;
      B_M1_px = -999;
      B_M1_px = -999;
      B_M1_pz = -999;

      B_M2_pt = -999;
      B_M2_eta = -999;
      B_M2_phi = -999;
      B_M2_px = -999;
      B_M2_px = -999;
      B_M2_pz = -999;
    }
  }
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

  //edm::Service<TFileService> fs;
  //tree_ = fs->make<TTree>("ntuple"," J/psi ntuple");

  tree_ = new TTree("ntuple", "ntuple");

  tree_->Branch("nB", &nB, "nB/i");
  tree_->Branch("savedtriggerNames", &savedtriggerNames);
  tree_->Branch("savedtriggerBits", &savedtriggerBits);
  tree_->Branch("savedtriggerPrescales", &savedtriggerPrescales);

  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  tree_->Branch("Event", &Event);

  tree_->Branch("B_J1_mass", &B_J1_mass);
  tree_->Branch("B_J1_px", &B_J1_px);
  tree_->Branch("B_J1_py", &B_J1_py);
  tree_->Branch("B_J1_pz", &B_J1_pz);
  tree_->Branch("B_J1_pt", &B_J1_pt);
  tree_->Branch("B_J1_eta", &B_J1_eta);
  tree_->Branch("B_J1_phi", &B_J1_phi);
  tree_->Branch("B_J1_rapidity", &B_J1_rapidity);

  tree_->Branch("B_J1_VtxPx", &B_J1_VtxPx);
  tree_->Branch("B_J1_VtxPy", &B_J1_VtxPy);
  tree_->Branch("B_J1_VtxPz", &B_J1_VtxPz);
  tree_->Branch("B_J1_VtxPt", &B_J1_VtxPt);
  tree_->Branch("B_J1_VtxEta", &B_J1_VtxEta);
  tree_->Branch("B_J1_VtxPhi", &B_J1_VtxPhi);
  tree_->Branch("B_J1_VtxRapidity", &B_J1_VtxRapidity);
  tree_->Branch("B_J1_VtxMass", &B_J1_VtxMass);
  tree_->Branch("B_J1_PVx", &B_J1_PVx);
  tree_->Branch("B_J1_PVy", &B_J1_PVy);
  tree_->Branch("B_J1_PVz", &B_J1_PVz);
  tree_->Branch("B_J1_PVxError", &B_J1_PVxError);
  tree_->Branch("B_J1_PVyError", &B_J1_PVyError);
  tree_->Branch("B_J1_PVzError", &B_J1_PVzError);

  tree_->Branch("B_J1_VtxProb", &B_J1_VtxProb);

  tree_->Branch("B_Mu1_px", &B_Mu1_px);
  tree_->Branch("B_Mu1_py", &B_Mu1_py);
  tree_->Branch("B_Mu1_pz", &B_Mu1_pz);
  tree_->Branch("B_Mu1_pt", &B_Mu1_pt);
  tree_->Branch("B_Mu1_eta", &B_Mu1_eta);
  tree_->Branch("B_Mu1_phi", &B_Mu1_phi);
  tree_->Branch("B_Mu1_soft", &B_Mu1_soft);
  tree_->Branch("B_Mu1_tight", &B_Mu1_tight);
  tree_->Branch("B_Mu1_loose", &B_Mu1_loose);
  tree_->Branch("B_Mu1_charge", &B_Mu1_charge);

  tree_->Branch("B_Mu2_px", &B_Mu2_px);
  tree_->Branch("B_Mu2_py", &B_Mu2_py);
  tree_->Branch("B_Mu2_pz", &B_Mu2_pz);
  tree_->Branch("B_Mu2_pt", &B_Mu2_pt);
  tree_->Branch("B_Mu2_eta", &B_Mu2_eta);
  tree_->Branch("B_Mu2_phi", &B_Mu2_phi);
  tree_->Branch("B_Mu2_soft", &B_Mu2_soft);
  tree_->Branch("B_Mu2_tight", &B_Mu2_tight);
  tree_->Branch("B_Mu2_loose", &B_Mu2_loose);
  tree_->Branch("B_Mu2_charge", &B_Mu2_charge);

  tree_->Branch("B_M1_pt", &B_M1_pt);
  tree_->Branch("B_M1_eta", &B_M1_eta);
  tree_->Branch("B_M1_phi", &B_M1_phi);
  tree_->Branch("B_M1_px", &B_M1_px);
  tree_->Branch("B_M1_py", &B_M1_py);
  tree_->Branch("B_M1_pz", &B_M1_pz);

  tree_->Branch("B_M2_pt", &B_M2_pt);
  tree_->Branch("B_M2_eta", &B_M2_eta);
  tree_->Branch("B_M2_phi", &B_M2_phi);
  tree_->Branch("B_M2_px", &B_M2_px);
  tree_->Branch("B_M2_py", &B_M2_py);
  tree_->Branch("B_M2_pz", &B_M2_pz);
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

  // Specify that only 'tracks' is allowed
  // To use, remove the default given above and uncomment below
  // ParameterSetDescription desc;
  // desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  // descriptions.addWithDefaultLabel(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(miniAODmmmm);
