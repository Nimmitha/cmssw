// -*- C++ -*-
//
// Package:    trigFireRateAnalyzer/miniAODtrigRate
// Class:      miniAODtrigRate
//
/**\class miniAODtrigRate miniAODtrigRate.cc trigFireRateAnalyzer/miniAODtrigRate/plugins/miniAODtrigRate.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nimmitha Karunarathna
//         Created:  Tue, 21 May 2024 21:55:18 GMT
//
//

// system include files
#include <memory>

// user include files
#include "trigFireRateAnalyzer/miniAODtrigRate/plugins/miniAODtrigRate.h"

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
// #include "DataFormats/TrackReco/interface/Track.h"
// #include "DataFormats/TrackReco/interface/TrackFwd.h"

// class declaration
//
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class miniAODtrigRate : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit miniAODtrigRate(const edm::ParameterSet &);
  ~miniAODtrigRate() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event &, const edm::EventSetup &) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> estoken_TTB;
  std::string MuonTriggerString;
 
  bool isMC_;

  TTree *tree_;

  unsigned int Run, LumiBlock, Event;

  unsigned int ntrigBits;
  bool TrigPathAvailable;
  bool TriggerPath;
  unsigned int nMuons;
  unsigned int nB;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

typedef math::Error<3>::type CovarianceMatrix;

using namespace reco;
using namespace edm;
using namespace std;

miniAODtrigRate::miniAODtrigRate(const edm::ParameterSet &iConfig)
    : muonsToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
      // estoken_MF(esConsumes()),
      // estoken_TTB(esConsumes<TransientTrackBuilder, TransientTrackRecord>()),
      estoken_TTB(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      MuonTriggerString(iConfig.getParameter<std::string>("MuonTrigger")),

      isMC_(iConfig.getParameter<bool>("isMC")),

      tree_(0),
      Run(0), LumiBlock(0), Event(0),
      ntrigBits(0),
      TrigPathAvailable(false), TriggerPath(false),
      nMuons(0),
      nB(0)
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  // now do what ever initialization is needed
}

miniAODtrigRate::~miniAODtrigRate()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void miniAODtrigRate::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  using namespace edm;

  edm::Handle<pat::MuonCollection> thePATMuonHandle;
  iEvent.getByToken(muonsToken_, thePATMuonHandle);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  if (!thePATMuonHandle.isValid())
  {
    edm::LogWarning("miniAODtrigRate") << "No pat::Muon found on Event!";
    return;
  }

  if (!triggerBits.isValid())
  {
    edm::LogWarning("miniAODtrigRate") << "No TriggerResults found on Event!";
    return;
  }

  if (!triggerObjects.isValid())
  {
    edm::LogWarning("miniAODtrigRate") << "No TriggerObjectStandAlone found on Event!";
    return;
  }


  // First Work on Trigger Information//
  //**********************************
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  // std::cout << "\n === TRIGGER PATHS === " << std::endl;

  ntrigBits = triggerBits->size();

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
  {
    if (names.triggerName(i).find("HLT_IsoMu24_v") != string::npos) // || names.triggerName(i).find("HLT_IsoTkMu24_v") != string::npos)
    {
      TrigPathAvailable = true;
      if (triggerBits->accept(i))
      {
        TriggerPath = true;
      }
    }
  }

  //****************************************************
  //*********Now we get the muons***********************
  //****************************************************

  // length of the muon collection

  Run = iEvent.id().run();
  LumiBlock = iEvent.luminosityBlock();
  Event = iEvent.id().event();

  nMuons = thePATMuonHandle->size();

  for (pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1)
  {
    for (pat::MuonCollection::const_iterator iMuon2 = iMuon1 + 1; iMuon2 != thePATMuonHandle->end(); ++iMuon2)
    {
      for (pat::MuonCollection::const_iterator iMuon3 = iMuon2 + 1; iMuon3 != thePATMuonHandle->end(); ++iMuon3)
      {
        for (pat::MuonCollection::const_iterator iMuon4 = iMuon3 + 1; iMuon4 != thePATMuonHandle->end(); ++iMuon4)
        {
          // cout << "Start of 4 muon loop" << endl;

          // make sure all muons are diferent
          if (iMuon1 == iMuon2)
            continue;
          if (iMuon1 == iMuon3)
            continue;
          if (iMuon1 == iMuon4)
            continue;
          if (iMuon2 == iMuon3)
            continue;
          if (iMuon2 == iMuon4)
            continue;
          if (iMuon3 == iMuon4)
            continue;
          // opposite charge

          // only look for neutral charge combination
          if (!(abs(iMuon1->charge()) == 1))
            continue;
          if (!(abs(iMuon2->charge()) == 1))
            continue;
          if (!(abs(iMuon3->charge()) == 1))
            continue;
          if (!(abs(iMuon4->charge()) == 1))
            continue;

          if (!((iMuon1->charge()) + (iMuon2->charge()) + (iMuon3->charge()) + (iMuon4->charge()) == 0))
            continue;
          // cout<<"Addition of 4 muon charge "<<(iMuon1->charge()) + (iMuon2->charge()) + (iMuon3->charge()) + (iMuon4->charge())<<endl;
          // preselection ptcut
          if (iMuon1->pt() < 2.0)
            continue;
          if (iMuon2->pt() < 2.0)
            continue;
          if (iMuon3->pt() < 2.0)
            continue;
          if (iMuon4->pt() < 2.0)
            continue;
          TrackRef glbTrackP1;
          TrackRef glbTrackP2;
          TrackRef glbTrackM1;
          TrackRef glbTrackM2;

          if (iMuon1->charge() == 1)
          {
            glbTrackP1 = iMuon1->track();
            if (iMuon2->charge() == 1)
            {
              glbTrackP2 = iMuon2->track();
              glbTrackM1 = iMuon3->track();
              glbTrackM2 = iMuon4->track();
            }
            else if (iMuon3->charge() == 1)
            {
              glbTrackP2 = iMuon3->track();
              glbTrackM1 = iMuon2->track();
              glbTrackM2 = iMuon4->track();
            }
            else if (iMuon4->charge() == 1)
            {
              glbTrackP2 = iMuon4->track();
              glbTrackM1 = iMuon2->track();
              glbTrackM2 = iMuon3->track();
            }
            else
            {
              cout << "Something is wrong while making +1 glb track ref" << endl;
            }
          }

          if (iMuon1->charge() == -1)
          {
            glbTrackM1 = iMuon1->track();
            if (iMuon2->charge() == -1)
            {
              glbTrackM2 = iMuon2->track();
              glbTrackP1 = iMuon3->track();
              glbTrackP2 = iMuon4->track();
            }
            else if (iMuon3->charge() == -1)
            {
              glbTrackM2 = iMuon3->track();
              glbTrackP1 = iMuon2->track();
              glbTrackP2 = iMuon4->track();
            }
            else if (iMuon4->charge() == -1)
            {
              glbTrackM2 = iMuon4->track();
              glbTrackP1 = iMuon2->track();
              glbTrackP2 = iMuon3->track();
            }
            else
            {
              cout << "Something is wrong while making -1 glb track ref" << endl;
            }
          }

          if (glbTrackP1.isNull() || glbTrackM1.isNull() || glbTrackP2.isNull() || glbTrackM2.isNull())
          {
            // std::cout << "continue due to no track ref" << endl;
            continue;
          }

          TLorentzVector M1, M2, M3, M4, MM1, MM2, MM3, MM4, MMMM;
          // initialize 4 lepton mass
          float mu_mass = 0.1056583745; //[PDG mass]
          // float ele_mass =  0.000510998928;//PDG mass

          // make muon 4 vectors

          if (iMuon1->charge() == 1)
          {
            if (iMuon2->charge() == 1)
            {
              // P1->Mu1;P2->Mu2;M1->Mu3;M2->Mu4
              M1.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M3.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M2.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M4.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            }
            else if (iMuon3->charge() == 1)
            {
              // P1->Mu1;P2->Mu3;M1->Mu2;M2->Mu4
              M1.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M2.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M3.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M4.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            }
            else if (iMuon4->charge() == 1)
            {
              // P1->Mu1;P2->Mu4;M1->Mu2;M2->Mu3
              M1.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M2.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M4.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M3.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            }
            else
            {
              cout << "Something is wrong while making +1 charge lorentz vector" << endl;
            }
          }
          else if (iMuon1->charge() == -1)
          {
            if (iMuon2->charge() == -1)
            {
              // P1->Mu3;P2->Mu4;M1->Mu1;M2->Mu2
              M3.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M1.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M4.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M2.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            }
            else if (iMuon3->charge() == -1)
            {
              // P1->Mu2;P2->Mu4;M1->Mu1;M2->Mu3
              M2.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M1.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M4.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M3.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            }
            else if (iMuon4->charge() == -1)
            {
              // P1->Mu2;P2->Mu3;M1->Mu1;M2->Mu4
              M2.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M1.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M3.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M4.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            }
            else
            {
              cout << "Something is wrong while making -1 charge lorentz vector" << endl;
            }
          }
          // now M1 is first positive muon 4 vector
          // then M2 is first negative muon 4 vector
          // then M3 is second positive muon 4 vector
          // finally, M4 is second negative muon
          //****************************************************************************************
          // make netral dimuon combination
          MM1 = M1 + M2;
          MM2 = M3 + M4;
          MM3 = M2 + M3;
          MM4 = M1 + M4;

          MMMM = M1 + M2 + M3 + M4;

          // cout<<"Start looking muon track quality"<<endl;
          if (!(glbTrackM1->quality(reco::TrackBase::highPurity)))
            continue;
          if (!(glbTrackP1->quality(reco::TrackBase::highPurity)))
            continue;
          if (!(glbTrackM2->quality(reco::TrackBase::highPurity)))
            continue;
          if (!(glbTrackP2->quality(reco::TrackBase::highPurity)))
            continue;
          // cout<<"Start Building Track"<<endl;
          nB++;
        }
      }
    }
  }


  tree_->Fill();

  Run = 0;
  LumiBlock = 0;
  Event = 0;

  ntrigBits = 0;
  TrigPathAvailable = false;
  TriggerPath = false;
  nMuons = 0;
  nB = 0;


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void miniAODtrigRate::beginJob()
{
  // please remove this method if not needed

  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  tree_ = new TTree("ntuple", "ntuple");

  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  tree_->Branch("Event", &Event);

  tree_->Branch("ntrigBits", &ntrigBits);
  tree_->Branch("TrigPathAvailable", &TrigPathAvailable);
  tree_->Branch("TriggerPath", &TriggerPath);
  tree_->Branch("nMuons", &nMuons);
  tree_->Branch("nB", &nB, "nB/i");
}

// ------------ method called once each job just after ending the event loop  ------------
void miniAODtrigRate::endJob()
{
  // please remove this method if not needed
  tree_->GetDirectory()->cd();
  tree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void miniAODtrigRate::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
  // The following says we do not know what parameters are allowed so do no validation
  //  Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  // Specify that only 'tracks' is allowed
  // To use, remove the default given above and uncomment below
  // ParameterSetDescription desc;
  // desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  // descriptions.addWithDefaultLabel(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(miniAODtrigRate);
