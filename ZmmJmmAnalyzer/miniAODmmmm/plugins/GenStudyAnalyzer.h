#ifndef _GenStudyAnalyzer_h
#define _GenStudyAnalyzer_h

// system include files
#include <memory>

// user include files
#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/GenStudyAnalyzer.h"

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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// for pileup information
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// class declaration
//
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class GenStudyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit GenStudyAnalyzer(const edm::ParameterSet &);
  ~GenStudyAnalyzer() override;

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
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupToken_;

  // edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> estoken_TTB;
  // edm::ESHandle<TransientTrackBuilder> theB;

  // const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> estoken_MF;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> estoken_TTB;
  std::string MuonTriggerString;
  bool isMC_;

  TTree *tree_;

  UInt_t Run;
  UShort_t LumiBlock;
  ULong64_t Event;
  float pileup;

  bool isFiducialGen;
  UShort_t nGenFidJpsi;
  UInt_t B_J1_Gen_pt;
  Short_t B_J1_Gen_eta;

  bool TriggerFired;

  UShort_t nRecoJpsiNoVtx;
  UShort_t nRecoDistinctJWVtx;

  float B_J1_mass;
  UInt_t B_J1_pt;
  Short_t B_J1_rapidity;

  float B_J1_vtxProb;

  UInt_t B_Mu1_pt;
  UInt_t B_Mu2_pt;
  Short_t B_Mu1_eta;
  Short_t B_Mu2_eta;
  bool B_Mu1_soft;
  bool B_Mu2_soft;

  UShort_t nMatchedMuons;
};
#endif