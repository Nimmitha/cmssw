#ifndef _miniAODmmmm_h
#define _miniAODmmmm_h

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

class miniAODmmmm : public edm::EDAnalyzer {
public:
  explicit miniAODmmmm(const edm::ParameterSet &);
  ~miniAODmmmm();

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

  std::string MuonTriggerString;
  bool isMC_;

  TTree *tree_;

  std::vector<float> *Run, *LumiBlock, *Event;

  std::vector<bool> *Mu_TriggerPath;
  std::vector<float> *B_U_TriggerPt1, *B_U_TriggerEta1, *B_U_TriggerPhi1;
  std::vector<float> *B_U_TriggerPt2, *B_U_TriggerEta2, *B_U_TriggerPhi2;
  std::vector<float> *B_U_TriggerPt3, *B_U_TriggerEta3, *B_U_TriggerPhi3;
  std::vector<float> *B_U_TriggerPt4, *B_U_TriggerEta4, *B_U_TriggerPhi4;
  std::vector<float> *B_U_TriggerPt5, *B_U_TriggerEta5, *B_U_TriggerPhi5;

  std::vector<float> *B_J1_mass, *B_J1_px, *B_J1_py, *B_J1_pz;
  std::vector<float> *B_J1_pt, *B_J1_eta, *B_J1_phi, *B_J1_rapidity;

  std::vector<float> *B_J1_VtxPx, *B_J1_VtxPy, *B_J1_VtxPz;
  std::vector<float> *B_J1_VtxPt, *B_J1_VtxEta, *B_J1_VtxPhi, *B_J1_VtxRapidity, *B_J1_VtxMass;

  std::vector<float> *B_J1_PVx, *B_J1_PVy, *B_J1_PVz;
  std::vector<float> *B_J1_PVxError, *B_J1_PVyError, *B_J1_PVzError;

  std::vector<float> *B_Mu1_px, *B_Mu1_py, *B_Mu1_pz;
  std::vector<float> *B_Mu1_pt, *B_Mu1_eta, *B_Mu1_phi;
  std::vector<bool> *B_Mu1_soft, *B_Mu1_tight, *B_Mu1_loose;
  std::vector<float> *B_Mu1_IsoTrack, *B_Mu1_IsoHcal, *B_Mu1_IsoEcal, *B_Mu1_IsoCalo;

  std::vector<float> *B_Mu1_PaperIsoTrackRF04, *B_Mu1_PaperIsoTrackRF03;
  std::vector<float> *B_Mu1_Paper3DIP;

  std::vector<int> *B_Mu1_charge;

  std::vector<float> *B_J1_VtxProb, *B_J2_VtxProb, *B_J3_VtxProb, *B_J4_VtxProb;
  std::vector<float> *B_J_xyP1, *B_J_xyM1, *B_J_zP1, *B_J_zM1;
  std::vector<float> *B_J_xyP2, *B_J_xyM2, *B_J_zP2, *B_J_zM2;

  std::vector<float> *mu1mC2;
  std::vector<int> *mu1mNHits, *mu1mNPHits;
  std::vector<float> *mu1pC2;
  std::vector<int> *mu1pNHits, *mu1pNPHits;
  std::vector<float> *B_M1_pt, *B_M1_eta, *B_M1_phi;
  std::vector<float> *B_M1_px, *B_M1_py, *B_M1_pz;
  std::vector<float> *B_J_GenMuonPt, *B_J_GenMuonEta, *B_J_GenMuonPhi;
  std::vector<float> *B_Z_GenMuonPt, *B_Z_GenMuonEta, *B_Z_GenMuonPhi;
  unsigned int nB;
};
#endif
