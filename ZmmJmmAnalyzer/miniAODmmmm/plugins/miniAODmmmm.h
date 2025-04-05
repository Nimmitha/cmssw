#ifndef _miniAODmmmm_h
#define _miniAODmmmm_h

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

// trigger
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
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

// class declaration
//
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class miniAODmmmm : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit miniAODmmmm(const edm::ParameterSet &);
  ~miniAODmmmm() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event &, const edm::EventSetup &) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<edm::TriggerResults> TriggerResultsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;

  // edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> estoken_TTB;
  // edm::ESHandle<TransientTrackBuilder> theB;

  // const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> estoken_MF;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> estoken_TTB;
  std::string MuonTriggerString;
  bool isMC_;

  TTree *tree_;

  std::vector<std::string> savedtriggerNames;
  std::vector<bool> *savedtriggerBits;
  std::vector<double> *savedtriggerPrescales;

  std::vector<float> *Run, *LumiBlock, *Event;

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

  std::vector<float> *B_Mu2_px, *B_Mu2_py, *B_Mu2_pz;
  std::vector<float> *B_Mu2_pt, *B_Mu2_eta, *B_Mu2_phi;
  std::vector<int> *B_Mu1_charge, *B_Mu2_charge;
  std::vector<bool> *B_Mu2_soft, *B_Mu2_tight, *B_Mu2_loose;
  std::vector<float> *B_Mu2_IsoTrack, *B_Mu2_IsoHcal, *B_Mu2_IsoEcal, *B_Mu2_IsoCalo;

  std::vector<float> *B_Mu2_PaperIsoTrackRF04, *B_Mu2_PaperIsoTrackRF03;
  std::vector<float> *B_Mu2_Paper3DIP;

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
  std::vector<float> *B_J_GenMuonPt, *B_J_GenMuonEta, *B_J_GenMuonPhi;
  std::vector<float> *B_Z_GenMuonPt, *B_Z_GenMuonEta, *B_Z_GenMuonPhi;
  unsigned int nB;

  // #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  //   edm::ESGetToken<SetupData, SetupRecord> setupToken_;
  // #endif
};
#endif