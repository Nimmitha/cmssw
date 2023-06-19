// -*- C++ -*-
//
// Package:    Muon/MyAnalyzer
// Class:      MyAnalyzer
//
/**\class MyAnalyzer MyAnalyzer.cc Muon/MyAnalyzer/plugins/MyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nimmitha Karunarathna
//         Created:  Thu, 15 Jun 2023 23:43:35 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include <TFile.h> //??

//kalman vertexing
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;
edm::EDGetTokenT<std::vector<pat::Muon> > muonCollToken;

class MyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MyAnalyzer(const edm::ParameterSet&);
  ~MyAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  TH1D* mumuMass_;
  TH1D* pt_;
  TH1D* eta_;
  TH1D* phi_;

  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  // double trackPtMin_;
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
MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))) {
  //now do what ever initialization is needed
  edm::InputTag muonTag("slimmedMuons");
  muonCollToken = consumes<pat::MuonCollection>(muonTag);

  edm::Service<TFileService> fs;
  mumuMass_ = fs->make<TH1D>("mumuMass", "Z Candidates in Di-Muon Channel", 100, 0, 120);
  pt_ = fs->make<TH1D>("pt", "pt", 100, 0, 120);
  eta_ = fs->make<TH1D>("eta", "eta", 100, -4, 4);
  phi_ = fs->make<TH1D>("phi", "phi", 100, -4, 4);
}

MyAnalyzer::~MyAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void MyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;

  //Kalman Vtx----------------------//
	// vector <TransientTrack> mu_tks;
	// vector <TransientTrack> mmmm_tks;
	KalmanVertexFitter kvf(true);

  edm::Handle<pat::MuonCollection> thePATMuonHandle;
  iEvent.getByToken(muonCollToken, thePATMuonHandle);

  // get the collection of muons from the event
  // edm::Handle<pat::MuonCollection> thePATMuonHandle;
  // edm::Handle<std::vector<Muon> > thePATMuonHandle;
  // iEvent.getByLabel("slimmedMuons", thePATMuonHandle);

  if (!thePATMuonHandle.isValid()) {
    cout << "No PAT Muons in event!" << endl;
    return;
  }

  // length of the muon collection
  cout << "thePATMuonHandle->size() = " << thePATMuonHandle->size() << endl;

  // construct Z candidates
  for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
    // eta cut and pt cut
        pt_->Fill(iMuon1->pt());
        eta_->Fill(iMuon1->eta());
        phi_->Fill(iMuon1->phi());
    if (iMuon1->pt()<3 || fabs(iMuon1->eta()) > 2.4)
      continue;
    for (pat::MuonCollection::const_iterator iMuon2 = iMuon1 + 1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
      if (iMuon1->charge() * iMuon2->charge() < 0) {
        // eta cut
        if (iMuon2->pt()<3 || fabs(iMuon2->eta()) > 2.4)
          continue;
        mumuMass_->Fill((iMuon1->p4() + iMuon2->p4()).mass());
        // cout << "muon1 pt = " << iMuon1->pt() << ", muon2 pt = " << iMuon2->pt() << endl;
      }
    }
  }

  // int nTrack = 0;
  // for (const auto& track : iEvent.get(tracksToken_)) {
  //   // do something with track parameters, e.g, plot the charge.
  //   // int charge = track.charge();
  //   if (track.pt() < trackPtMin_)
  //     continue;
  //   nTrack++;
  // }
  // cout << "nTrack = " << nTrack << endl;
}

// ------------ method called once each job just before starting event loop  ------------
void MyAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void MyAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyAnalyzer);
