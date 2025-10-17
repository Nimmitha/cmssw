// system include files
#include <memory>

// user include files
#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/TriggerEffAnalyzer.h"

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

// // trigger
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

// // packedCandidates
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// // From Kalman example
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

TriggerEffAnalyzer::TriggerEffAnalyzer(const edm::ParameterSet &iConfig)
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

      isMC_(iConfig.getParameter<bool>("isMC")),

      tree_(0),

      Run(0),
      LumiBlock(0),
      Event(0),
      passSel1(false),
      passSel2(false),
      sameMuonTrigMatch(false),
      nRecoDistinctJWVtx(0),
      // prescaleBarrel(0), // optional
      // prescaleL1DM(0), // optional
      B_J1_mass(0),
      B_J1_pt(0),
      B_J1_rapidity(0),

      B_Mu1_pt(0),
      B_Mu2_pt(0),
      B_Mu1_eta(0),
      B_Mu2_eta(0) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

TriggerEffAnalyzer::~TriggerEffAnalyzer() {}

//
// member functions
//

// ------------ method called for each event  ------------
void TriggerEffAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
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

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  if (!thePATMuonHandle.isValid()) {
    edm::LogWarning("TriggerEffAnalyzer") << "No pat::Muon found on Event!";
    return;
  }
  if (!TriggerResults.isValid()) {
    edm::LogWarning("TriggerEffAnalyzer") << "No TriggerResults found on Event!";
    return;
  }
  if (!triggerPrescales.isValid()) {
    edm::LogWarning("TriggerEffAnalyzer") << "no Trigger prescale in event!";
    return;
  }

  // Struct to hold dimuon candidates
  struct DimuonCandidate {
    const pat::Muon *mu1;
    const pat::Muon *mu2;
    float vtxProb;
    TLorentzVector p4;
    // XYZTLorentzVectorD J_vtx;
    float mu1_pt, mu2_pt, mu1_eta, mu2_eta;

    bool operator<(const DimuonCandidate &other) const {
      return vtxProb > other.vtxProb;  // sort descending by vtxProb
    }
  };

  std::vector<DimuonCandidate> validCandidates;

  // Reset variables for the next event
  Run = iEvent.id().run();
  LumiBlock = iEvent.id().luminosityBlock();
  Event = iEvent.id().event();

  const edm::TriggerNames &names = iEvent.triggerNames(*TriggerResults);

  // -------------------------------
  // Selection 1: HLT_Mu0_Barrel_v*, muon pT>5, |eta|<1.5
  // -------------------------------
  passSel1 = false;
  passSel2 = false;
  bool firedHLT_Barrel = false;
  for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
    std::string name = names.triggerName(i);
    if (name.find("HLT_Mu0_Barrel_v") != std::string::npos && TriggerResults->accept(i)) {
      firedHLT_Barrel = true;
      break;
    }
  }

  if (firedHLT_Barrel) {
    for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
      if (iMuon1->pt() > 5.0 && std::abs(iMuon1->eta()) < 1.5) {
        passSel1 = true;
        break;
      }
    }
  }

  if (passSel1) {
    // cout << "Event passed selection 1" << endl;
  } else {
    return;
  }

  // -------------------------------
  // Selection 2: HLT_Mu0_L1DoubleMu_v*, dimuon selection
  // -------------------------------
  bool firedHLT_L1DM = false;
  for (unsigned int i = 0; i < TriggerResults->size(); ++i) {
    std::string name = names.triggerName(i);
    if (name.find("HLT_Mu0_L1DoubleMu_v") != std::string::npos && TriggerResults->accept(i)) {
      firedHLT_L1DM = true;
      break;
    }
  }

  if (firedHLT_L1DM) {
    // cout << "Event fired HLT_Mu0_L1DoubleMu" << endl;
  } else {
    // cout << "Event did not fire HLT_Mu0_L1DoubleMu" << endl;
    // return;
  }

  if (firedHLT_L1DM) {
    cout << "Start making dimuon candidates" << endl;
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

        if (iMuon1->eta() < -1.5 || iMuon1->eta() > 1.5)
          continue;
        if (iMuon2->eta() < -1.5 || iMuon2->eta() > 1.5)
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
          cout << "Something is wrong while making glb track ref" << endl;
        }

        if (glbTrackP1.isNull() || glbTrackM1.isNull()) {
          //std::cout << "continue due to no track ref" << endl;
          continue;
        }

        TLorentzVector M1, M2, MM1;
        float mu_mass = 0.1056583745;  //[PDG mass]

        //make muon 4 vectors
        if (iMuon1->charge() == 1 && iMuon2->charge() == -1) {
          M1.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
          M2.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
        } else if (iMuon1->charge() == -1 && iMuon2->charge() == 1) {
          M1.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
          M2.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
        } else {
          cout << "Something is wrong while making charge lorentz vector" << endl;
        }

        //make netral dimuon combination
        MM1 = M1 + M2;

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

        if (!J_candi1.isValid()) {
          // cout << "continue because no vertexed dimuon" << endl;
          continue;
        }

        reco::Vertex JPsi_Vtx1 = J_candi1;

        float B_Prob_tmp1 = TMath::Prob(J_candi1.totalChiSquared(), J_candi1.degreesOfFreedom());
        // const math::XYZTLorentzVectorD JPsi_mom1 = JPsi_Vtx1.p4(mu_mass, 0.0);

        // if (B_Prob_tmp1 < 0.1) {
        //   continue;
        // }

        // Remove events with mass outside J/Psi or Z mass window
        if ((MM1.M() < 2.6 || MM1.M() > 3.6)) {
          continue;
        }

        validCandidates.push_back({&(*iMuon1),
                                   &(*iMuon2),
                                   B_Prob_tmp1,
                                   MM1,
                                   //  JPsi_mom1,
                                   static_cast<float>(iMuon1->pt()),
                                   static_cast<float>(iMuon2->pt()),
                                   static_cast<float>(iMuon1->eta()),
                                   static_cast<float>(iMuon2->eta())});
      }
    }
    cout << "Number of valid dimuon candidates: " << validCandidates.size() << endl;

    if (validCandidates.size() > 0) {
      // Process valid candidates
      std::sort(validCandidates.begin(), validCandidates.end());

      // Track used muons
      std::set<const pat::Muon *> usedMuons;
      usedMuons.clear();
      nRecoDistinctJWVtx = 0;
      for (const auto &cand : validCandidates) {
        if (usedMuons.count(cand.mu1) || usedMuons.count(cand.mu2)) {
          continue;
        }
        nRecoDistinctJWVtx++;
        cout << "Working on candidate: " << nRecoDistinctJWVtx << endl;
        usedMuons.insert(cand.mu1);
        usedMuons.insert(cand.mu2);

        if (nRecoDistinctJWVtx == 1) {
          cout << "Saving candidate variables" << endl;
          // Fill the tree with all the candidates
          B_J1_mass = cand.p4.M();
          B_J1_pt = std::lround(cand.p4.Pt() * 1000);
          B_J1_rapidity = std::lround(cand.p4.Rapidity() * 1000);
          B_J1_VtxProb = cand.vtxProb;
          B_Mu1_pt = std::lround(cand.mu1->pt() * 1000);
          B_Mu2_pt = std::lround(cand.mu2->pt() * 1000);
          B_Mu1_eta = std::lround(cand.mu1->eta() * 100);
          B_Mu2_eta = std::lround(cand.mu2->eta() * 100);

        }
      }
      passSel2 = true;
    }
  }

  if (!passSel2) {
    cout << "Event did not pass selection 2" << endl;
    B_J1_mass = -999;
    B_J1_pt = 1000;
    B_J1_rapidity = -999;
  
    B_J1_VtxProb = -999;
  
    B_Mu1_pt = 1000;
    B_Mu2_pt = 1000;
    B_Mu1_eta = 10;
    B_Mu2_eta = 10;
  }


  // -------------------------------
  // ΔR matching between triggers (0.15)
  // -------------------------------
  sameMuonTrigMatch = false;
  if (passSel1 && passSel2) {
    std::cout << "**Start trigger object matching**" << std::endl;

    const pat::Muon *mu1 = validCandidates.front().mu1;
    const pat::Muon *mu2 = validCandidates.front().mu2;

    bool mu1_Barrel = false, mu1_DoubleMu = false;
    bool mu2_Barrel = false, mu2_DoubleMu = false;

    for (auto trigObj : *triggerObjects) {
      trigObj.unpackPathNames(names);

      // --- mu1 matching ---
      if (reco::deltaR(*mu1, trigObj) < 0.15) {
        if (trigObj.hasPathName("HLT_Mu0_Barrel_v*", true, true))
          mu1_Barrel = true;
        if (trigObj.hasPathName("HLT_Mu0_L1DoubleMu_v*", true, true))
          mu1_DoubleMu = true;
      }

      // --- mu2 matching ---
      if (reco::deltaR(*mu2, trigObj) < 0.15) {
        if (trigObj.hasPathName("HLT_Mu0_Barrel_v*", true, true))
          mu2_Barrel = true;
        if (trigObj.hasPathName("HLT_Mu0_L1DoubleMu_v*", true, true))
          mu2_DoubleMu = true;
      }
    }

    // Check if either muon fired both triggers
    if ((mu1_Barrel && mu1_DoubleMu) || (mu2_Barrel && mu2_DoubleMu)) {
      sameMuonTrigMatch = true;
      std::cout << "→ One of the two muons matched both triggers!" << std::endl;
    } else {
      std::cout << "→ No muon matched both triggers." << std::endl;
    }
  }

  tree_->Fill();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------

void TriggerEffAnalyzer::beginJob() {
  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  tree_ = new TTree("ntuple", "ntuple");

  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  // tree_->Branch("Event", &Event);

  // --- New tag-and-probe helper branches ---
  tree_->Branch("passSel1", &passSel1, "passSel1/O");
  tree_->Branch("passSel2", &passSel2, "passSel2/O");
  tree_->Branch("sameMuonTrigMatch", &sameMuonTrigMatch, "sameMuonTrigMatch/O");
  tree_->Branch("nRecoDistinctJWVtx", &nRecoDistinctJWVtx, "nRecoDistinctJWVtx/I");
  // tree_->Branch("HLT_Mu0_Barrel_prescale", &prescaleBarrel, "HLT_Mu0_Barrel_prescale/I"); // optional
  // tree_->Branch("HLT_Mu0_L1DoubleMu_prescale", &prescaleL1DM, "HLT_Mu0_L1DoubleMu_prescale/I");

  tree_->Branch("B_J1_mass", &B_J1_mass);
  tree_->Branch("B_J1_pt", &B_J1_pt);
  tree_->Branch("B_J1_VtxProb", &B_J1_VtxProb);

  tree_->Branch("B_Mu1_pt", &B_Mu1_pt);
  tree_->Branch("B_Mu2_pt", &B_Mu2_pt);
  tree_->Branch("B_Mu1_eta", &B_Mu1_eta);
  tree_->Branch("B_Mu2_eta", &B_Mu2_eta);
}

// ------------ method called once each job just after ending the event loop  ------------
void TriggerEffAnalyzer::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TriggerEffAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  //  Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(TriggerEffAnalyzer);