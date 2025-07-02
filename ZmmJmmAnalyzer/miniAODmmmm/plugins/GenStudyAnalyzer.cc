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

// #include <cmath>
// constants, enums and typedefs
//

typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//

GenStudyAnalyzer::GenStudyAnalyzer(const edm::ParameterSet &iConfig)
    : muonsToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      TriggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      prunedGenToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned"))),
      triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
      triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
      primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      // estoken_MF(esConsumes()),
      // estoken_TTB(esConsumes<TransientTrackBuilder, TransientTrackRecord>()),
      estoken_TTB(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      MuonTriggerString(iConfig.getParameter<std::string>("MuonTrigger")),

      isMC_(iConfig.getParameter<bool>("isMC")) {
  //now do what ever initialization is needed
  pileupToken_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));
}

GenStudyAnalyzer::~GenStudyAnalyzer() {}

//
// member functions
//

// ------------ method called for each event  ------------
void GenStudyAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  const auto &theB = iSetup.getData(estoken_TTB);

  edm::Handle<pat::MuonCollection> thePATMuonHandle;
  iEvent.getByToken(muonsToken_, thePATMuonHandle);

  edm::Handle<edm::TriggerResults> TriggerResults;
  iEvent.getByToken(TriggerResultsToken_, TriggerResults);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
  iEvent.getByToken(pileupToken_, puInfo);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(prunedGenToken_, genParticles);
  if (!genParticles.isValid()) {
    edm::LogWarning("GenStudyAnalyzer") << "No pruned GenParticles found on Event!";
    return;
  }

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  if (!thePATMuonHandle.isValid()) {
    edm::LogWarning("GenStudyAnalyzer") << "No pat::Muon found on Event!";
    return;
  }
  if (!TriggerResults.isValid()) {
    edm::LogWarning("GenStudyAnalyzer") << "No TriggerResults found on Event!";
    return;
  }
  if (!triggerPrescales.isValid()) {
    edm::LogWarning("GenStudyAnalyzer") << "no Trigger prescale in event!";
    return;
  }

  // edm::LogInfo("GenStudyAnalyzer") << "Processing event: " << iEvent.id().event() << " in run: " << iEvent.id().run();
  // cout << "LumiBlock: " << iEvent.id().luminosityBlock() << std::endl;

  // Pileup Info
  pileup = -1;
  if (puInfo.isValid()) {
    for (const auto &pu : *puInfo) {
      if (pu.getBunchCrossing() == 0) {  // Only take in-time PU
        pileup = pu.getTrueNumInteractions();
        break;
      }
    }
  }

  // check if we have a valid J/Psi candidate in the genParticles
  isFiducialGen = false;
  nGenFidJpsi = 0;
  B_J1_Gen_pt = 0;

  for (const auto &jpsi : *genParticles) {
    if (jpsi.pdgId() != 443)
      continue;

    if (nGenFidJpsi > 1) {  // Sanity check
      nGenFidJpsi++;
      continue;
    }

    // Adding some sanity check logging
    if (jpsi.status() != 2) {  // J/Psi should have status 2
      edm::LogInfo("GenStudyAnalyzer") << "Skipped J/Psi with status " << jpsi.status();
      continue;
    }

    if (jpsi.numberOfDaughters() < 2 || jpsi.numberOfDaughters() > 2) {
      edm::LogInfo("GenStudyAnalyzer") << "Skipped J/Psi with " << jpsi.numberOfDaughters() << " daughters";
      continue;
    }

    if (jpsi.isPromptDecayed() == false) {
      edm::LogInfo("GenStudyAnalyzer") << "Skipped J/Psi that is not prompt decayed";
      continue;
    }

    // assgin the two daughters to muons
    if (fabs(jpsi.daughter(0)->pdgId()) != 13 || fabs(jpsi.daughter(1)->pdgId()) != 13) {
      edm::LogInfo("GenStudyAnalyzer") << "Skipped J/Psi with daughters not muons";
      continue;
    }
    if (jpsi.daughter(0)->status() != 1 || jpsi.daughter(1)->status() != 1) {
      edm::LogInfo("GenStudyAnalyzer") << "Skipped J/Psi with daughters not stable";
      continue;
    }

    // if we reach here, we have a valid J/Psi candidate

    // Gen Muon fiducial cuts
    TLorentzVector genMu1, genMu2;
    const reco::GenParticle *muon1ptr = nullptr;
    const reco::GenParticle *muon2ptr = nullptr;

    if (jpsi.daughter(0)->pt() > jpsi.daughter(1)->pt()) {
      genMu1.SetPxPyPzE(jpsi.daughter(0)->px(), jpsi.daughter(0)->py(), jpsi.daughter(0)->pz(), jpsi.daughter(0)->energy());
      genMu2.SetPxPyPzE(jpsi.daughter(1)->px(), jpsi.daughter(1)->py(), jpsi.daughter(1)->pz(), jpsi.daughter(1)->energy());
      muon1ptr = dynamic_cast<const reco::GenParticle *>(jpsi.daughter(0));
      muon2ptr = dynamic_cast<const reco::GenParticle *>(jpsi.daughter(1));
    } else {
      genMu2.SetPxPyPzE(jpsi.daughter(0)->px(), jpsi.daughter(0)->py(), jpsi.daughter(0)->pz(), jpsi.daughter(0)->energy());
      genMu1.SetPxPyPzE(jpsi.daughter(1)->px(), jpsi.daughter(1)->py(), jpsi.daughter(1)->pz(), jpsi.daughter(1)->energy());
      muon1ptr = dynamic_cast<const reco::GenParticle *>(jpsi.daughter(1));
      muon2ptr = dynamic_cast<const reco::GenParticle *>(jpsi.daughter(0));
    }

    if (muon1ptr->charge() + muon2ptr->charge() != 0)
      continue;
    if (genMu1.Pt() < 5.0 || genMu2.Pt() < 5.0)
      continue;
    if (fabs(genMu1.Eta()) > 2.4 || fabs(genMu2.Eta()) > 2.4)
      continue;

    B_J1_Gen_pt = std::lround(jpsi.pt() * 1000);

    isFiducialGen = true;
    nGenFidJpsi++;
  }

  // Sanity check
  if (nGenFidJpsi > 1) {
    edm::LogInfo("GenStudyAnalyzer") << "Found " << nGenFidJpsi << " J/Psi Gen candidates";
  }

  const edm::TriggerNames &names = iEvent.triggerNames(*TriggerResults);

  // Check if the trigger fired
  TriggerFired = false;
  for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
    if (names.triggerName(i).find(MuonTriggerString.c_str()) != string::npos) {
      if (TriggerResults->accept(i)) {
        TriggerFired = true;  // Trigger fired!
      }
      break;  // No need to check further if we already found the trigger
    }
  }

  //Now we get the primary vertex
  reco::Vertex bestVtx;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  bestVtx = *(primaryVertices_handle->begin());

  //***************Now we check reco***********************
  B_J1_mass = -1;
  B_J1_pt = 0;
  B_J1_rapidity = -2999;

  B_Mu1_pt = 0;
  B_Mu2_pt = 0;
  B_Mu1_eta = -399;
  B_Mu2_eta = -399;

  Run = iEvent.id().run();
  LumiBlock = iEvent.luminosityBlock();
  Event = iEvent.id().event();

  B_J1_vtxProb = -1;
  nRecoJpsi = 0;

  for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
    for (pat::MuonCollection::const_iterator iMuon2 = iMuon1 + 1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
      //make sure all muons are diferent
      if (iMuon1 == iMuon2)
        continue;

      if (!(abs(iMuon1->charge()) == 1))
        continue;
      if (!(abs(iMuon2->charge()) == 1))
        continue;

      //opposite charge
      if (!((iMuon1->charge()) + (iMuon2->charge())) == 0)
        continue;

      if (iMuon1->pt() < 5.0)
        continue;
      if (iMuon2->pt() < 5.0)
        continue;

      if (iMuon1->eta() < -2.4 || iMuon1->eta() > 2.4)
        continue;
      if (iMuon2->eta() < -2.4 || iMuon2->eta() > 2.4)
        continue;

      TrackRef glbTrackP1;
      TrackRef glbTrackM1;

      if (iMuon1->charge() == 1 && iMuon2->charge() == -1) {
        glbTrackP1 = iMuon1->track();
        glbTrackM1 = iMuon2->track();
      } else if (iMuon1->charge() == -1 && iMuon2->charge() == 1) {
        glbTrackP1 = iMuon2->track();
        glbTrackM1 = iMuon1->track();
      } else {
        edm::LogInfo("GenStudyAnalyzer") << "Something is wrong while making charge track reference";
        continue;
      }

      if (glbTrackP1.isNull() || glbTrackM1.isNull()) {
        //std::cout << "continue due to no track ref" << endl;
        continue;
      }

      TLorentzVector m1, m2, MM1;
      float mu_mass = 0.1056583745;  //[PDG mass]

      //make muon 4 vectors
      if (iMuon1->charge() == 1 && iMuon2->charge() == -1) {
        m1.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
        m2.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
      } else if (iMuon1->charge() == -1 && iMuon2->charge() == 1) {
        m1.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
        m2.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
      } else {
        edm::LogInfo("GenStudyAnalyzer") << "Something is wrong while making muon 4 vectors";
        continue;
      }

      //make netral dimuon combination
      MM1 = m1 + m2;

      //cout<<"Start looking muon track quality"<<endl;
      if (!(glbTrackM1->quality(reco::TrackBase::highPurity)))
        continue;
      if (!(glbTrackP1->quality(reco::TrackBase::highPurity)))
        continue;

      reco::TransientTrack muon1TT;
      reco::TransientTrack muon2TT;

      muon1TT = theB.build(glbTrackP1);
      muon2TT = theB.build(glbTrackM1);

      //Kalman Vtx----------------------//
      vector<TransientTrack> mu_tks;

      KalmanVertexFitter kvfM(true);
      mu_tks.clear();
      mu_tks.push_back(muon1TT);
      mu_tks.push_back(muon2TT);
      TransientVertex J_candi1 = kvfM.vertex(mu_tks);

      if (!J_candi1.isValid())
        continue;

      reco::Vertex JPsi_Vtx1 = J_candi1;

      float Jpsi_vtxProb = TMath::Prob(J_candi1.totalChiSquared(), J_candi1.degreesOfFreedom());
      // const math::XYZTLorentzVectorD JPsi_mom1 = JPsi_Vtx1.p4(mu_mass, 0.0);

      // Remove events with mass outside J/Psi mass window
      if ((MM1.M() < 2.6 || MM1.M() > 3.6))
        continue;

      nRecoJpsi++;

      if (Jpsi_vtxProb > B_J1_vtxProb) {
        B_J1_vtxProb = Jpsi_vtxProb;

        B_J1_mass = MM1.M();
        B_J1_pt = std::lround(MM1.Pt() * 1000);
        B_J1_rapidity = std::lround(MM1.Rapidity() * 1000);

        B_Mu1_pt = std::lround(iMuon1->pt() * 1000);
        B_Mu2_pt = std::lround(iMuon2->pt() * 1000);
        B_Mu1_eta = std::lround(iMuon1->eta() * 100);
        B_Mu2_eta = std::lround(iMuon2->eta() * 100);
      }
    }
  }

  // Fill the tree with the event information
  tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------

void GenStudyAnalyzer::beginJob() {
  tree_ = new TTree("ntuple", "ntuple");

  tree_->Branch("Run", &Run, "Run/i");
  tree_->Branch("LumiBlock", &LumiBlock, "LumiBlock/s");
  tree_->Branch("Event", &Event, "Event/l");
  tree_->Branch("pileup", &pileup, "pileup/F");

  tree_->Branch("isFiducialGen", &isFiducialGen, "isFiducialGen/O");
  tree_->Branch("nGenFidJpsi", &nGenFidJpsi, "nGenFidJpsi/s");
  tree_->Branch("nRecoJpsi", &nRecoJpsi, "nRecoJpsi/s");

  tree_->Branch("B_J1_Gen_pt", &B_J1_Gen_pt, "B_J1_Gen_pt/i");

  tree_->Branch("TriggerFired", &TriggerFired, "TriggerFired/O");

  tree_->Branch("B_J1_mass", &B_J1_mass, "B_J1_mass/F");
  tree_->Branch("B_J1_pt", &B_J1_pt, "B_J1_pt/i");
  tree_->Branch("B_J1_rapidity", &B_J1_rapidity, "B_J1_rapidity/S");

  tree_->Branch("B_J1_vtxProb", &B_J1_vtxProb, "B_J1_vtxProb/F");

  tree_->Branch("B_Mu1_pt", &B_Mu1_pt, "B_Mu1_pt/i");
  tree_->Branch("B_Mu2_pt", &B_Mu2_pt, "B_Mu2_pt/i");
  tree_->Branch("B_Mu1_eta", &B_Mu1_eta, "B_Mu1_eta/S");
  tree_->Branch("B_Mu2_eta", &B_Mu2_eta, "B_Mu2_eta/S");
}

// ------------ method called once each job just after ending the event loop  ------------
void GenStudyAnalyzer::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenStudyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  //  Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(GenStudyAnalyzer);