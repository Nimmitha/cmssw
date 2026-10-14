#ifndef _miniAODmmmm_MuonComp_h
#define _miniAODmmmm_MuonComp_h

// system include files
#include <memory>

// user include files
#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/miniAODmmmm_MuonComp.h"

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

class miniAODmmmm_MuonComp : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit miniAODmmmm_MuonComp(const edm::ParameterSet &);
  ~miniAODmmmm_MuonComp() override;

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

  UInt_t Run;
  UShort_t LumiBlock;
  ULong64_t Event;
  ULong64_t eventTime;
  bool TriggerFired;

  UShort_t nPV;
  float B_J1_mass;
  UInt_t B_J1_pt;
  Short_t B_J1_rapidity;

  float B_J1_VtxPt, B_J1_VtxMass;
  float B_J1_VtxProb;

  // Per muon (mu1, mu2 of the candidate), see muonType() and fillHits() in the .cc
  UChar_t mu1_type, mu2_type;              // bits: 0 global, 1 tracker, 2 standalone, 3 PF, 4 soft, 5 loose, 6 medium, 7 tight
  UChar_t mu1_pixValid, mu2_pixValid;      // valid hit, bits 0-3: BPix L1-L4, bits 4-6: FPix disks 1-3
  UChar_t mu1_pixLost, mu2_pixLost;        // hit expected on the track but missing, same bits (TRACK_HITS)
  UChar_t mu1_firstHit, mu2_firstHit;      // first valid hit of the inner track: 10*subdet + layer (1 BPix, 2 FPix, 3 TIB, 4 TID, 5 TOB, 6 TEC)
  UChar_t mu1_nPixLay, mu2_nPixLay;        // pixelLayersWithMeasurement
  UChar_t mu1_nTrkLay, mu2_nTrkLay;        // trackerLayersWithMeasurement
  UChar_t mu1_nLostInner, mu2_nLostInner;  // missing hits before the first valid hit (MISSING_INNER_HITS), all tracker
  UChar_t mu1_nLostInnerPix, mu2_nLostInnerPix;  // the same, pixel only

  UInt_t B_Mu1_pt;
  UInt_t B_Mu2_pt;
  Short_t B_Mu1_eta;
  Short_t B_Mu2_eta;
  bool isBestCandidate;

  void fillHits(const reco::Track &tk, UChar_t &valid, UChar_t &lost, UChar_t &first, UChar_t &nPixLay,
                UChar_t &nTrkLay, UChar_t &nLostInner, UChar_t &nLostInnerPix) const;
  UChar_t muonType(const pat::Muon &mu, const reco::Vertex &vtx) const;
};
#endif