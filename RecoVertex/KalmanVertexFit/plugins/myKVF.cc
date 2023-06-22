// system include files
#include <memory>

// user include files
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/SimpleVertexTree.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include <TFile.h>


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include <iostream>

/**
   * This is a very simple test analyzer mean to test the KalmanVertexFitter
   */

class myKVF : public edm::one::EDAnalyzer<> {
public:
  explicit myKVF(const edm::ParameterSet&);
  ~myKVF() override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  void beginJob() override;
  void endJob() override;

private:
  TrackingVertex getSimVertex(const edm::Event& iEvent) const;

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> estoken_MF;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> estoken_TTB;

  edm::ParameterSet theConfig;
  edm::ParameterSet kvfPSet;
  std::unique_ptr<SimpleVertexTree> tree;
  TFile* rootFile_;

  std::string outputFile_;  // output file
  edm::EDGetTokenT<reco::TrackCollection> token_tracks;
};

using namespace reco;
using namespace edm;
using namespace std;

myKVF::myKVF(const edm::ParameterSet& iConfig)
    : estoken_MF(esConsumes()),
      estoken_TTB(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      theConfig(iConfig) {
  token_tracks = consumes<TrackCollection>(iConfig.getParameter<string>("TrackLabel"));
  outputFile_ = iConfig.getUntrackedParameter<std::string>("outputFile");
  kvfPSet = iConfig.getParameter<edm::ParameterSet>("KVFParameters");
  rootFile_ = TFile::Open(outputFile_.c_str(), "RECREATE");
  edm::LogInfo("RecoVertex/myKVF") << "Initializing KVF TEST analyser  - Output file: " << outputFile_ << "\n";
}

myKVF::~myKVF() { delete rootFile_; }

void myKVF::beginJob() {}

void myKVF::endJob() {}

//
// member functions
//

void myKVF::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::LogInfo("RecoVertex/myKVF") << "Reconstructing event number: " << iEvent.id() << "\n";

  // get RECO tracks from the event
  // `tks` can be used as a ptr to a reco::TrackCollection
//   edm::Handle<edm::View<reco::Track> > tks;
  edm::Handle<reco::TrackCollection> tks;
//   edm::Handle<vector<reco::Track>> tks;
  iEvent.getByToken(token_tracks, tks);

  if (!tks.isValid()) {
    edm::LogPrint("RecoVertex/myKVF") << "Exception during event number: " << iEvent.id() << "\n";
  } else {
    edm::LogPrint("RecoVertex/myKVF") << "got " << (*tks).size() << " tracks " << std::endl;

    // Transform Track to TransientTrack

    //get the builder:
    const auto& theB = &iSetup.getData(estoken_TTB);
    //do the conversion:
    std::vector<TransientTrack> t_tks = theB->build(tks);

    edm::LogPrint("RecoVertex/myKVF") << "Found: " << t_tks.size() << " reconstructed tracks"
                                       << "\n";

    // Call the KalmanVertexFitter if more than 1 track
    if (t_tks.size() > 1) {
      //      KalmanVertexFitter kvf(kvfPSet);
      KalmanVertexFitter kvf(true);
      TransientVertex tv = kvf.vertex(t_tks);

      edm::LogPrint("RecoVertex/myKVF") << "Position: " << Vertex::Point(tv.position()) << "\n";

    }
  }
}
DEFINE_FWK_MODULE(myKVF);