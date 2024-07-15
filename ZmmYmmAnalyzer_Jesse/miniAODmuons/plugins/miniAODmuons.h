#ifndef _miniAODmuons_h
#define _miniAODmuons_h

// system include files
#include <memory>
#include <map>
// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"  // muy importante para MiniAOD

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"

//
// class decleration
//

class miniAODmuons : public edm::EDAnalyzer {
public:
  explicit miniAODmuons(const edm::ParameterSet &);
  ~miniAODmuons();

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event &, const edm::EventSetup &);
  virtual void endJob();
  void printout(const RefCountedKinematicVertex &myVertex) const;
  void printout(const RefCountedKinematicParticle &myParticle) const;
  void printout(const RefCountedKinematicTree &myTree) const;

  // ----------member data ---------------------------

  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  edm::EDGetTokenT<edm::View<pat::Electron>> dielectron_Label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;

  //trigger------------
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  //edm::EDGetTokenT<reco::GenParticleCollection> GenCollection_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> prunedGenToken_;

  bool isMC_;

  TTree *tree_;

  std::vector<float> *Run, *LumiBlock, *Event;
  std::vector<float> *FourL_mass;
  std::vector<float> *FourL_px, *FourL_py, *FourL_pz;
  std::vector<float> *FourL_pt, *FourL_eta, *FourL_phi;
  std::vector<float> *FourL_VtxProb;
  std::vector<float> *FourL_PVx, *FourL_PVy, *FourL_PVz;
  std::vector<float> *FourL_PVxError, *FourL_PVyError, *FourL_PVzError;

  std::vector<bool> *B_U_TriggerPath;
  std::vector<float> *B_U_TriggerPt1, *B_U_TriggerEta1, *B_U_TriggerPhi1;
  std::vector<float> *B_U_TriggerPt2, *B_U_TriggerEta2, *B_U_TriggerPhi2;
  std::vector<float> *B_U_TriggerPt3, *B_U_TriggerEta3, *B_U_TriggerPhi3;
  std::vector<float> *B_U_TriggerPt4, *B_U_TriggerEta4, *B_U_TriggerPhi4;
  std::vector<float> *B_U_TriggerPt5, *B_U_TriggerEta5, *B_U_TriggerPhi5;

  std::vector<float> *B_J1_mass, *B_J1_px, *B_J1_py, *B_J1_pz;
  std::vector<float> *B_J1_pt, *B_J1_eta, *B_J1_phi, *B_J1_rapidity;
  std::vector<float> *B_J2_mass, *B_J2_px, *B_J2_py, *B_J2_pz;
  std::vector<float> *B_J2_pt, *B_J2_eta, *B_J2_phi, *B_J2_rapidity;
  std::vector<float> *B_J3_mass, *B_J3_px, *B_J3_py, *B_J3_pz;
  std::vector<float> *B_J3_pt, *B_J3_eta, *B_J3_phi, *B_J3_rapidity;
  std::vector<float> *B_J4_mass, *B_J4_px, *B_J4_py, *B_J4_pz;
  std::vector<float> *B_J4_pt, *B_J4_eta, *B_J4_phi, *B_J4_rapidity;

  std::vector<float> *B_J1_VtxPx, *B_J1_VtxPy, *B_J1_VtxPz;
  std::vector<float> *B_J1_VtxPt, *B_J1_VtxEta, *B_J1_VtxPhi, *B_J1_VtxRapidity, *B_J1_VtxMass;

  std::vector<float> *B_J1_PVx, *B_J1_PVy, *B_J1_PVz;
  std::vector<float> *B_J1_PVxError, *B_J1_PVyError, *B_J1_PVzError;
  std::vector<float> *B_J2_VtxPx, *B_J2_VtxPy, *B_J2_VtxPz;
  std::vector<float> *B_J2_VtxPt, *B_J2_VtxEta, *B_J2_VtxPhi, *B_J2_VtxRapidity, *B_J2_VtxMass;
  std::vector<float> *B_J2_PVx, *B_J2_PVy, *B_J2_PVz;
  std::vector<float> *B_J2_PVxError, *B_J2_PVyError, *B_J2_PVzError;

  std::vector<float> *B_J3_VtxPx, *B_J3_VtxPy, *B_J3_VtxPz;
  std::vector<float> *B_J3_VtxPt, *B_J3_VtxEta, *B_J3_VtxPhi, *B_J3_VtxRapidity, *B_J3_VtxMass;
  std::vector<float> *B_J3_PVx, *B_J3_PVy, *B_J3_PVz;
  std::vector<float> *B_J3_PVxError, *B_J3_PVyError, *B_J3_PVzError;
  std::vector<float> *B_J4_VtxPx, *B_J4_VtxPy, *B_J4_VtxPz;
  std::vector<float> *B_J4_VtxPt, *B_J4_VtxEta, *B_J4_VtxPhi, *B_J4_VtxRapidity, *B_J4_VtxMass;
  std::vector<float> *B_J4_PVx, *B_J4_PVy, *B_J4_PVz;
  std::vector<float> *B_J4_PVxError, *B_J4_PVyError, *B_J4_PVzError;

  std::vector<float> *B_Mu1_px, *B_Mu1_py, *B_Mu1_pz;
  std::vector<float> *B_Mu1_pt, *B_Mu1_eta, *B_Mu1_phi;
  std::vector<bool> *B_Mu1_soft, *B_Mu1_tight, *B_Mu1_loose;
  std::vector<float> *B_Mu1_IsoTrack, *B_Mu1_IsoHcal, *B_Mu1_IsoEcal, *B_Mu1_IsoCalo;

  std::vector<float> *B_Mu1_PaperIsoTrackRF04, *B_Mu1_PaperIsoTrackRF03;
  std::vector<float> *B_Mu1_Paper3DIP;

  std::vector<float> *B_Mu2_px, *B_Mu2_py, *B_Mu2_pz;
  std::vector<float> *B_Mu2_pt, *B_Mu2_eta, *B_Mu2_phi;
  std::vector<int> *B_Mu1_charge, *B_Mu2_charge;
  std::vector<bool> *B_Mu2_soft, *B_Mu2_tight, *B_Mu2_loose;
  std::vector<float> *B_Mu2_IsoTrack, *B_Mu2_IsoHcal, *B_Mu2_IsoEcal, *B_Mu2_IsoCalo;

  std::vector<float> *B_Mu2_PaperIsoTrackRF04, *B_Mu2_PaperIsoTrackRF03;
  std::vector<float> *B_Mu2_Paper3DIP;

  std::vector<float> *B_Mu3_px, *B_Mu3_py, *B_Mu3_pz;
  std::vector<float> *B_Mu3_pt, *B_Mu3_eta, *B_Mu3_phi;
  std::vector<bool> *B_Mu3_soft, *B_Mu3_tight, *B_Mu3_loose;
  std::vector<int> *B_Mu3_charge;
  std::vector<float> *B_Mu3_IsoTrack, *B_Mu3_IsoHcal, *B_Mu3_IsoEcal, *B_Mu3_IsoCalo;

  std::vector<float> *B_Mu3_PaperIsoTrackRF04, *B_Mu3_PaperIsoTrackRF03;
  std::vector<float> *B_Mu3_Paper3DIP;

  std::vector<float> *B_Mu4_px, *B_Mu4_py, *B_Mu4_pz;
  std::vector<float> *B_Mu4_pt, *B_Mu4_eta, *B_Mu4_phi;
  std::vector<bool> *B_Mu4_soft, *B_Mu4_tight, *B_Mu4_loose;
  std::vector<int> *B_Mu4_charge;
  std::vector<float> *B_Mu4_IsoTrack, *B_Mu4_IsoHcal, *B_Mu4_IsoEcal, *B_Mu4_IsoCalo;

  std::vector<float> *B_Mu4_PaperIsoTrackRF04, *B_Mu4_PaperIsoTrackRF03;
  std::vector<float> *B_Mu4_Paper3DIP;

  std::vector<float> *B_J1_VtxProb, *B_J2_VtxProb, *B_J3_VtxProb, *B_J4_VtxProb;
  std::vector<float> *B_J_xyP1, *B_J_xyM1, *B_J_zP1, *B_J_zM1;
  std::vector<float> *B_J_xyP2, *B_J_xyM2, *B_J_zP2, *B_J_zM2;

  std::vector<float> *mu1mC2;
  std::vector<int> *mu1mNHits, *mu1mNPHits;
  std::vector<float> *mu1pC2;
  std::vector<int> *mu1pNHits, *mu1pNPHits;
  std::vector<float> *mu2mC2;
  std::vector<int> *mu2mNHits, *mu2mNPHits;
  std::vector<float> *mu2pC2;
  std::vector<int> *mu2pNHits, *mu2pNPHits;
  std::vector<float> *B_M1_pt, *B_M1_eta, *B_M1_phi;
  std::vector<float> *B_M1_px, *B_M1_py, *B_M1_pz;
  std::vector<float> *B_M2_pt, *B_M2_eta, *B_M2_phi;
  std::vector<float> *B_M2_px, *B_M2_py, *B_M2_pz;
  std::vector<float> *B_M3_pt, *B_M3_eta, *B_M3_phi;
  std::vector<float> *B_M3_px, *B_M3_py, *B_M3_pz;
  std::vector<float> *B_M4_pt, *B_M4_eta, *B_M4_phi;
  std::vector<float> *B_M4_px, *B_M4_py, *B_M4_pz;
  std::vector<float> *B_J_GenMuonPt, *B_J_GenMuonEta, *B_J_GenMuonPhi;
  std::vector<float> *B_Z_GenMuonPt, *B_Z_GenMuonEta, *B_Z_GenMuonPhi;
  unsigned int nB;
};
#endif
