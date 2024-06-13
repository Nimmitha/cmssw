#ifndef _miniAODzmmyee_h
#define _miniAODzmmyee_h

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

class miniAODzmmyee : public edm::EDAnalyzer {
public:
  explicit miniAODzmmyee(const edm::ParameterSet &);
  ~miniAODzmmyee();

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
  std::string ElectronTriggerString;
  std::string DataTypeString;
  bool isMC_;

  TTree *tree_;

  std::vector<float> *Y_GenMuonPt, *Y_GenMuonEta, *Y_GenMuonPhi;
  std::vector<float> *Z_GenMuonPt, *Z_GenMuonEta, *Z_GenMuonPhi;

  std::vector<float> *Run, *LumiBlock, *Event;
  std::vector<float> *FourL_mass;
  std::vector<float> *FourL_px, *FourL_py, *FourL_pz;
  std::vector<float> *FourL_pt, *FourL_eta, *FourL_phi, *FourL_rapidity;
  std::vector<float> *FourL_Vtx_Prob;
  std::vector<float> *FourL_Vtx_x, *FourL_Vtx_y, *FourL_Vtx_z;
  std::vector<float> *FourL_Vtx_xError, *FourL_Vtx_yError, *FourL_Vtx_zError;

  std::vector<float> *FourL_Vtx_Px, *FourL_Vtx_Py, *FourL_Vtx_Pz;
  std::vector<float> *FourL_Vtx_Pt, *FourL_Vtx_Eta, *FourL_Vtx_Phi;
  std::vector<float> *FourL_Vtx_Rapidity, *FourL_Vtx_Mass;

  std::vector<float> *Y_dca;
  std::vector<bool> *Mu_TriggerPath, *Ele_TriggerPath;
  std::vector<float> *Mu_TriggerPt1, *Mu_TriggerEta1, *Mu_TriggerPhi1;
  std::vector<float> *Mu_TriggerPt2, *Mu_TriggerEta2, *Mu_TriggerPhi2;
  std::vector<float> *Mu_TriggerPt3, *Mu_TriggerEta3, *Mu_TriggerPhi3;
  std::vector<float> *Mu_TriggerPt4, *Mu_TriggerEta4, *Mu_TriggerPhi4;
  std::vector<float> *Mu_TriggerPt5, *Mu_TriggerEta5, *Mu_TriggerPhi5;
  std::vector<float> *Ele_TriggerPt1, *Ele_TriggerEta1, *Ele_TriggerPhi1;
  std::vector<float> *Ele_TriggerPt2, *Ele_TriggerEta2, *Ele_TriggerPhi2;
  std::vector<float> *Ele_TriggerPt3, *Ele_TriggerEta3, *Ele_TriggerPhi3;
  std::vector<float> *Ele_TriggerPt4, *Ele_TriggerEta4, *Ele_TriggerPhi4;
  std::vector<float> *Ele_TriggerPt5, *Ele_TriggerEta5, *Ele_TriggerPhi5;
  std::vector<float> *Y_mass, *Y_Vtx_Prob, *Y_px, *Y_py, *Y_pz;
  std::vector<float> *Y_pt, *Y_eta, *Y_phi, *Y_rapidity;
  std::vector<float> *Y_Vtx_Px, *Y_Vtx_Py, *Y_Vtx_Pz;
  std::vector<float> *Y_Vtx_Pt, *Y_Vtx_Eta, *Y_Vtx_Phi, *Y_Vtx_Rapidity, *Y_Vtx_Mass;
  std::vector<float> *Y_Vtx_x, *Y_Vtx_y, *Y_Vtx_z;
  std::vector<float> *Y_Vtx_xError, *Y_Vtx_yError, *Y_Vtx_zError;

  std::vector<float> *Y_px1, *Y_py1, *Y_pz1;
  std::vector<float> *Y_pt1, *Y_eta1, *Y_SCeta1, *Y_phi1, *Y_energy1, *Y_energyCorr1;

  std::vector<float> *Y_ecalIso1, *Y_hcalIso1, *Y_trackIso1, *Z_trackIso1;
  std::vector<bool> *Y_CutBaseLoose1, *Y_CutBaseVeto1, *Y_mvaIsoWP90_1, *Y_mvaIsoWP80_1;
  std::vector<float> *Y_px2, *Y_py2, *Y_pz2;
  std::vector<float> *Y_pt2, *Y_eta2, *Y_SCeta2, *Y_phi2, *Y_energy2, *Y_energyCorr2;
  std::vector<float> *Y_ecalIso2, *Y_hcalIso2, *Y_trackIso2, *Z_trackIso2;
  std::vector<bool> *Y_CutBaseLoose2, *Y_CutBaseVeto2, *Y_mvaIsoWP90_2, *Y_mvaIsoWP80_2;
  std::vector<int> *Y_charge1, *Y_charge2;
  std::vector<float> *Y_fit_pt1, *Y_fit_ptError1, *Y_fit_eta1, *Y_fit_phi1;
  std::vector<float> *Y_fit_pt2, *Y_fit_ptError2, *Y_fit_eta2, *Y_fit_phi2;
  std::vector<float> *Y_dxy1, *Y_dxy2, *Y_dz1, *Y_dz2;

  std::vector<float> *Y_lowPt;
  std::vector<float> *Y_highPt;
  std::vector<float> *Z_lowPt;
  std::vector<float> *Z_highPt;
  std::vector<float> *Mu_dCA;
  std::vector<float> *Ele_dCA;

  std::vector<float> *Z_mass, *Z_px, *Z_py, *Z_pz;
  std::vector<float> *Z_pt, *Z_eta, *Z_phi, *Z_rapidity;
  std::vector<float> *Z_Vtx_Px, *Z_Vtx_Py, *Z_Vtx_Pz;
  std::vector<float> *Z_Vtx_Pt, *Z_Vtx_Eta, *Z_Vtx_Phi, *Z_Vtx_Rapidity, *Z_Vtx_Mass;
  std::vector<float> *Z_Vtx_x, *Z_Vtx_y, *Z_Vtx_z;
  std::vector<float> *Z_Vtx_xError, *Z_Vtx_yError, *Z_Vtx_zError;

  std::vector<float> *Z_px1, *Z_py1, *Z_pz1;
  std::vector<float> *Z_pt1, *Z_eta1, *Z_phi1;
  std::vector<bool> *Z_soft1, *Z_tight1, *Z_loose1;
  std::vector<float> *Z_fit_pt1, *Z_fit_ptError1, *Z_fit_eta1, *Z_fit_phi1;

  std::vector<float> *Z_px2, *Z_py2, *Z_pz2;
  std::vector<float> *Z_pt2, *Z_eta2, *Z_phi2;
  std::vector<int> *Z_charge1, *Z_charge2;
  std::vector<bool> *Z_soft2, *Z_tight2, *Z_loose2;
  std::vector<float> *Z_fit_pt2, *Z_fit_ptError2, *Z_fit_eta2, *Z_fit_phi2;
  std::vector<float> *Z_Vtx_Prob;
  std::vector<float> *Z_xy1, *Z_xy2, *Z_z2, *Z_z1;

  std::vector<float> *Mu1_C2;
  std::vector<int> *Mu1_NHits, *Mu1_NPHits;
  std::vector<float> *Mu2_C2;
  std::vector<int> *Mu2_NHits, *Mu2_NPHits;

  unsigned int nCandi;
  unsigned int nEleKalman;
  unsigned int nEleRefitTrack;
};
#endif
