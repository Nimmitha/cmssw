// -*- C++ -*-
//
// Package:    Muon/MyProducer
// Class:      MyProducer
//
/**\class MyProducer MyProducer.cc Muon/MyProducer/plugins/MyProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nimmitha Karunarathna
//         Created:  Tue, 25 Jul 2023 23:26:33 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"


#include <vector>
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"

//
// class declaration
//

class MyProducer : public edm::stream::EDProducer<> {
public:
  explicit MyProducer(const edm::ParameterSet&);
  ~MyProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::Candidate>> prunedGenParticlesToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> slimmedMuonsToken_;
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
MyProducer::MyProducer(const edm::ParameterSet& config) {
    // Input collection token
    prunedGenParticlesToken_ = consumes<edm::View<reco::Candidate>>(config.getParameter<edm::InputTag>("prunedGenParticles"));
    slimmedMuonsToken_ = consumes<edm::View<reco::Candidate>>(config.getParameter<edm::InputTag>("slimmedMuons"));

    // Register the output
    // produces<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate>>>();
    produces<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate>>>("prunedGenParticles");
    produces<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate>>>("slimmedMuons");
}

MyProducer::~MyProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void MyProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
    // Get the prunedGenParticles collection from the event
    edm::Handle<edm::View<reco::Candidate>> prunedGenParticlesHandle;
    event.getByToken(prunedGenParticlesToken_, prunedGenParticlesHandle);

    // Get the slimmedMuons collection from the event
    edm::Handle<edm::View<reco::Candidate>> slimmedMuonsHandle;
    event.getByToken(slimmedMuonsToken_, slimmedMuonsHandle);

    // Create the output collection
    std::unique_ptr<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate>>> prunedGenParticlesOutput(new edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate>>());
    std::unique_ptr<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate>>> slimmedMuonsOutput(new edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate>>());


    // Copy prunedGenParticles to the output collection
    for (const auto& particle : *prunedGenParticlesHandle) {
        prunedGenParticlesOutput->push_back(particle.clone());
    }

    // Copy slimmedMuons to the output collection
    for (const auto& muon : *slimmedMuonsHandle) {
        slimmedMuonsOutput->push_back(muon.clone());
    }

    // Put the output collections into the event
    event.put(std::move(prunedGenParticlesOutput), "prunedGenParticles");
    event.put(std::move(slimmedMuonsOutput), "slimmedMuons");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void MyProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void MyProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
MyProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
MyProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MyProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MyProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MyProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyProducer);

