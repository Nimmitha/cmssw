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
#include "TTree.h"
#include <TFile.h>
#include "TLorentzVector.h"

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
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> estoken_MF;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> estoken_TTB;
  edm::ParameterSet theConfig;

  TTree* tree_;
  std::vector<float>*Run, *LumiBlock, *Event;

  std::vector<float>*Z_mass, *Z_px, *Z_py, *Z_pz, *Z_pt;
  std::vector<float>*Z_eta, *Z_phi, *Z_vx, *Z_vy, *Z_vz, *Z_rapidity;
  std::vector<float>*Z_Vtx_mass, *Z_Vtx_Px, *Z_Vtx_Py, *Z_Vtx_Pz, *Z_Vtx_Pt;
  std::vector<float>*Z_Vtx_Eta, *Z_Vtx_Phi, *Z_Vtx_rapidity;
  std::vector<float>*Z_Vtx_x, *Z_Vtx_y, *Z_Vtx_z, *Z_Vtx_xError, *Z_Vtx_yError, *Z_Vtx_zError;
  std::vector<float>*Z_Vtx_chi2, *Z_Vtx_ndof, *Z_Vtx_Prob, *Z_Vtx_pull;

  TH1D* mumuMass_;
  TH1D* pt_;
  TH1D* eta_;
  TH1D* phi_;

  // edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
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
using namespace reco;
using namespace edm;
using namespace std;

MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig)
    // : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
    : estoken_MF(esConsumes()),
      estoken_TTB(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      theConfig(iConfig),

      // initialize the variables for the tree
      tree_(0),
      Run(0),
      LumiBlock(0),
      Event(0),
      Z_mass(0),
      Z_px(0),
      Z_py(0),
      Z_pz(0),
      Z_pt(0),
      Z_eta(0),

      Z_phi(0),
      Z_vx(0),
      Z_vy(0),
      Z_vz(0),
      Z_rapidity(0),
      Z_Vtx_mass(0),
      Z_Vtx_Px(0),
      Z_Vtx_Py(0),
      Z_Vtx_Pz(0),
      Z_Vtx_Pt(0),
      Z_Vtx_Eta(0),
      Z_Vtx_Phi(0),

      Z_Vtx_rapidity(0),

      Z_Vtx_x(0),
      Z_Vtx_y(0),
      Z_Vtx_z(0),

      Z_Vtx_xError(0),
      Z_Vtx_yError(0),
      Z_Vtx_zError(0),

      Z_Vtx_chi2(0),
      Z_Vtx_ndof(0),
      Z_Vtx_Prob(0),
      Z_Vtx_pull(0)

{
  //now do what ever initialization is needed
  edm::InputTag muonTag("slimmedMuons");
  muonCollToken = consumes<pat::MuonCollection>(muonTag);
  // edm::InputTag tracksTag("generalTracks");
  // token_tracks = consumes<TrackCollection>(tracksTag);

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
  edm::Handle<pat::MuonCollection> thePATMuonHandle;
  iEvent.getByToken(muonCollToken, thePATMuonHandle);

  const auto& theB = &iSetup.getData(estoken_TTB);

  if (!theB) {
    cout << "No TransientTrackBuilder in event!" << endl;
    return;
  }

  if (!thePATMuonHandle.isValid()) {
    cout << "No PAT Muons in event!" << endl;
    return;
  }

  // length of the muon collection
  cout << "thePATMuonHandle->size() = " << thePATMuonHandle->size() << endl;

  // construct Z candidates
  for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end();
       ++iMuon1) {
    // print track parameters
    reco::TrackRef glbTrackP;
    reco::TrackRef glbTrackM;

    // if (trackRef->quality(reco::TrackBase::highPurity)) {
    //   cout << "Track is high purity" << endl;
    // } else {
    //   cout << "Track is not high purity" << endl;
    // }

    // cout << tr

    // print track quality
    // cout << "Track quality: " << iMuon1->track()->quality(reco::TrackBase::highPurity) << endl;

    // eta cut and pt cut
    pt_->Fill(iMuon1->pt());
    eta_->Fill(iMuon1->eta());
    phi_->Fill(iMuon1->phi());
    if (iMuon1->pt() < 3 || fabs(iMuon1->eta()) > 2.4)
      continue;
    for (pat::MuonCollection::const_iterator iMuon2 = iMuon1 + 1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
      if (iMuon1->charge() * iMuon2->charge() < 0) {
        // eta cut
        if (iMuon2->pt() < 3 || fabs(iMuon2->eta()) > 2.4)
          continue;
        mumuMass_->Fill((iMuon1->p4() + iMuon2->p4()).mass());

        // Set P as positive muon and M as negative muon
        if (iMuon1->charge() == 1) {
          glbTrackP = iMuon1->track();
          glbTrackM = iMuon2->track();
        } else {
          glbTrackP = iMuon2->track();
          glbTrackM = iMuon1->track();
        }

        // check if track is available
        if (glbTrackP.isNull() || glbTrackM.isNull()) {
          cout << "Continue due ot no track ref" << endl;
          continue;
        }

        TLorentzVector M1, M2;
        TLorentzVector MM;
        float mu_mass = 0.1056583745;  //GeV [PDG]

        // Set M1 as positive muon 4-vector and M2 as negative muon 4-vector
        if (iMuon1->charge() == 1) {
          M1.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
          M2.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
        } else {
          M2.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
          M1.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
        }

        // Set MM as the dimuon 4-vector
        MM = M1 + M2;

        reco::TransientTrack glbTrackP_TT = theB->build(glbTrackP);
        reco::TransientTrack glbTrackM_TT = theB->build(glbTrackM);

        std::vector<reco::TransientTrack> mu_tks;
        KalmanVertexFitter kvf(true);

        // Transform Track to TransientTrack
        mu_tks.clear();
        mu_tks.push_back(glbTrackP_TT);
        mu_tks.push_back(glbTrackM_TT);

        TransientVertex Z_candi = kvf.vertex(mu_tks);

        if (!Z_candi.isValid()) {
          cout << "Continue due to invalid vertex" << endl;
          continue;
        }

        reco::Vertex Z_vertex = Z_candi;

        float Z_Prob = TMath::Prob(Z_candi.totalChiSquared(), Z_candi.degreesOfFreedom());
        const math::XYZTLorentzVector Z_mom = Z_vertex.p4(mu_mass, 0.0);

        cout << "Z total chi squared = " << Z_candi.totalChiSquared() << endl;
        cout << "Z degrees of freedom = " << Z_candi.degreesOfFreedom() << endl;
        cout << "Z probability = " << Z_Prob << endl;

        cout << "Z mass before fit = " << MM.M() << endl;
        cout << "Z mass after fit = " << Z_mom.mass() << endl;

        // Fill the tree
        Run->push_back(iEvent.id().run());
        LumiBlock->push_back(iEvent.luminosityBlock());
        Event->push_back(iEvent.id().event());

        Z_mass->push_back(MM.M());
        Z_px->push_back(MM.Px());
        Z_py->push_back(MM.Py());
        Z_pz->push_back(MM.Pz());
        Z_pt->push_back(MM.Pt());
        Z_eta->push_back(MM.Eta());
        Z_phi->push_back(MM.Phi());
        Z_rapidity->push_back(MM.Rapidity());

        Z_Vtx_mass->push_back(Z_mom.mass());
        Z_Vtx_Px->push_back(Z_mom.Px());
        Z_Vtx_Py->push_back(Z_mom.Py());
        Z_Vtx_Pz->push_back(Z_mom.Pz());
        Z_Vtx_Pt->push_back(Z_mom.Pt());
        Z_Vtx_Eta->push_back(Z_mom.Eta());
        Z_Vtx_Phi->push_back(Z_mom.Phi());
        Z_Vtx_rapidity->push_back(Z_mom.Rapidity());

        Z_Vtx_x->push_back(Z_vertex.x());
        Z_Vtx_y->push_back(Z_vertex.y());
        Z_Vtx_z->push_back(Z_vertex.z());
        Z_Vtx_xError->push_back(Z_vertex.xError());
        Z_Vtx_yError->push_back(Z_vertex.yError());
        Z_Vtx_zError->push_back(Z_vertex.zError());

        Z_Vtx_chi2->push_back(Z_candi.totalChiSquared());
        Z_Vtx_ndof->push_back(Z_candi.degreesOfFreedom());
        Z_Vtx_Prob->push_back(Z_Prob);
        Z_Vtx_pull->push_back((Z_mom.Pt() - MM.Pt())/ Z_candi.totalChiSquared());


        tree_->Fill();

        // clear the vectors
        Run->clear();
        LumiBlock->clear();
        Event->clear();

        Z_mass->clear();
        Z_px->clear();
        Z_py->clear();
        Z_pz->clear();
        Z_pt->clear();
        Z_eta->clear();
        Z_phi->clear();
        Z_rapidity->clear();

        Z_Vtx_mass->clear();
        Z_Vtx_Px->clear();
        Z_Vtx_Py->clear();
        Z_Vtx_Pz->clear();
        Z_Vtx_Pt->clear();
        Z_Vtx_Eta->clear();
        Z_Vtx_Phi->clear();
        Z_Vtx_rapidity->clear();

        Z_Vtx_x->clear();
        Z_Vtx_y->clear();
        Z_Vtx_z->clear();
        Z_Vtx_xError->clear();
        Z_Vtx_yError->clear();
        Z_Vtx_zError->clear();

        Z_Vtx_chi2->clear();
        Z_Vtx_ndof->clear();
        Z_Vtx_Prob->clear();
        Z_Vtx_pull->clear();
      }
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void MyAnalyzer::beginJob() {
  // please remove this method if not needed
  tree_ = new TTree("ntuple", "ntuple with info");
  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  tree_->Branch("Event", &Event);
  tree_->Branch("Z_mass", &Z_mass);
  tree_->Branch("Z_px", &Z_px);
  tree_->Branch("Z_py", &Z_py);
  tree_->Branch("Z_pz", &Z_pz);
  tree_->Branch("Z_pt", &Z_pt);
  tree_->Branch("Z_eta", &Z_eta);
  tree_->Branch("Z_phi", &Z_phi);
  tree_->Branch("Z_rapidity", &Z_rapidity);

  tree_->Branch("Z_Vtx_mass", &Z_Vtx_mass);
  tree_->Branch("Z_Vtx_Px", &Z_Vtx_Px);
  tree_->Branch("Z_Vtx_Py", &Z_Vtx_Py);
  tree_->Branch("Z_Vtx_Pz", &Z_Vtx_Pz);
  tree_->Branch("Z_Vtx_Pt", &Z_Vtx_Pt);
  tree_->Branch("Z_Vtx_Eta", &Z_Vtx_Eta);
  tree_->Branch("Z_Vtx_Phi", &Z_Vtx_Phi);
  tree_->Branch("Z_Vtx_rapidity", &Z_Vtx_rapidity);

  tree_->Branch("Z_Vtx_x", &Z_Vtx_x);
  tree_->Branch("Z_Vtx_y", &Z_Vtx_y);
  tree_->Branch("Z_Vtx_z", &Z_Vtx_z);
  tree_->Branch("Z_Vtx_xError", &Z_Vtx_xError);
  tree_->Branch("Z_Vtx_yError", &Z_Vtx_yError);
  tree_->Branch("Z_Vtx_zError", &Z_Vtx_zError);

  tree_->Branch("Z_Vtx_chi2", &Z_Vtx_chi2);
  tree_->Branch("Z_Vtx_ndof", &Z_Vtx_ndof);
  tree_->Branch("Z_Vtx_Prob", &Z_Vtx_Prob);
  tree_->Branch("Z_Vtx_pull", &Z_Vtx_pull);
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
