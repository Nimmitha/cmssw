#ifndef _mytagAndProbeV1_h
#define _mytagAndProbeV1_h

#include <memory>
#include <vector>
#include <iostream>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include "TLorentzVector.h"

// Trigger includes
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// Vertex includes
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class mytagAndProbeV1 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit mytagAndProbeV1(const edm::ParameterSet&);
  ~mytagAndProbeV1() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // Tokens
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<edm::TriggerResults> TriggerResultsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;

  bool isMC_;

  // TTree and Variables
  TTree* tree_;

  // Event Info
  UInt_t Run;
  UShort_t LumiBlock;
  ULong64_t Event;

  // Tag Info (The high pT muon matching IsoMu24)
  float tag_pt;
  float tag_eta;
  float tag_phi;

  // Probe Info (The other muon in the J/Psi pair)
  float probe_pt;
  float probe_eta;
  float probe_phi;

  // Pair Info
  float mass;

  // The Efficiency Bit
  // true if the event passed HLT_Mu0_L1DoubleMu_v*
  bool probe_passAlgo;
};

#endif