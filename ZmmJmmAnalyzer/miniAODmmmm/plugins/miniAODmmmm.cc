// -*- C++ -*-
//
// Package:    miniAODmmmm
// Class:      miniAODmmmm
//

//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Monday Aug 28 (2017)         |
//         <jhovanny.andres.mejia.guisao@cern.ch> |
//=================================================

// system include files
#include <memory>

// #include "myAnalyzers/JPsiKsPAT/src/miniAODmmmm.h"
#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/miniAODmmmm.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
//kalman vertexing
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//GenInfo
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH2F.h"

//
// constants, enums and typedefs
//

typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//

miniAODmmmm::miniAODmmmm(const edm::ParameterSet& iConfig)
    : dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
      dielectron_Label(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("dielectron"))),
      trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
      primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      //trigger
      triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
      //GenLevel Info
      prunedGenToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("pruned"))),
      MuonTriggerString(iConfig.getParameter<std::string>("MuonTrigger")),

      isMC_(iConfig.getParameter<bool>("isMC")),

      tree_(0),
      Run(0),
      LumiBlock(0),
      Event(0),

      //Z
      Mu_TriggerPath(0),
      B_U_TriggerPt1(0),
      B_U_TriggerEta1(0),
      B_U_TriggerPhi1(0),
      B_U_TriggerPt2(0),
      B_U_TriggerEta2(0),
      B_U_TriggerPhi2(0),
      B_U_TriggerPt3(0),
      B_U_TriggerEta3(0),
      B_U_TriggerPhi3(0),
      B_U_TriggerPt4(0),
      B_U_TriggerEta4(0),
      B_U_TriggerPhi4(0),
      B_U_TriggerPt5(0),
      B_U_TriggerEta5(0),
      B_U_TriggerPhi5(0),

      B_J1_mass(0),
      B_J1_px(0),
      B_J1_py(0),
      B_J1_pz(0),
      B_J1_pt(0),
      B_J1_eta(0),
      B_J1_phi(0),
      B_J1_rapidity(0),

      B_J1_VtxPx(0),
      B_J1_VtxPy(0),
      B_J1_VtxPz(0),
      B_J1_VtxPt(0),
      B_J1_VtxEta(0),
      B_J1_VtxPhi(0),
      B_J1_VtxRapidity(0),
      B_J1_VtxMass(0),
      B_J1_PVx(0),
      B_J1_PVy(0),
      B_J1_PVz(0),
      B_J1_PVxError(0),
      B_J1_PVyError(0),
      B_J1_PVzError(0),

      B_Mu1_px(0),
      B_Mu1_py(0),
      B_Mu1_pz(0),
      B_Mu1_pt(0),
      B_Mu1_eta(0),
      B_Mu1_phi(0),
      B_Mu1_soft(0),
      B_Mu1_tight(0),
      B_Mu1_loose(0),
      B_Mu1_IsoTrack(0),
      B_Mu1_IsoHcal(0),
      B_Mu1_IsoEcal(0),
      B_Mu1_IsoCalo(0),

      B_Mu1_PaperIsoTrackRF04(0),
      B_Mu1_PaperIsoTrackRF03(0),
      B_Mu1_Paper3DIP(0),

      B_Mu1_charge(0),

      B_J1_VtxProb(0),
      B_J_xyP1(0),
      B_J_xyM1(0),
      B_J_zP1(0),
      B_J_zM1(0),
      B_J_xyP2(0),
      B_J_xyM2(0),
      B_J_zP2(0),
      B_J_zM2(0),

      mu1mC2(0),
      mu1mNHits(0),
      mu1mNPHits(0),
      mu1pC2(0),
      mu1pNHits(0),
      mu1pNPHits(0),

      B_M1_pt(0),
      B_M1_eta(0),
      B_M1_phi(0),
      B_M1_px(0),
      B_M1_py(0),
      B_M1_pz(0),
      
      B_J_GenMuonPt(0),
      B_J_GenMuonEta(0),
      B_J_GenMuonPhi(0),
      B_Z_GenMuonPt(0),
      B_Z_GenMuonEta(0),
      B_Z_GenMuonPhi(0),

      nB(0)

{
  //now do what ever initialization is needed
}

miniAODmmmm::~miniAODmmmm() {}

//
// member functions
//

// ------------ method called to for each event  ------------
void miniAODmmmm::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  //*********************************
  // Get event content information
  //*********************************

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);

  edm::Handle<View<pat::PackedCandidate>> thePATTrackHandle;
  iEvent.getByToken(trakCollection_label, thePATTrackHandle);

  edm::Handle<View<pat::Muon>> thePATMuonHandle;
  iEvent.getByToken(dimuon_Label, thePATMuonHandle);

  edm::Handle<View<pat::Electron>> thePATElectronHandle;
  iEvent.getByToken(dielectron_Label, thePATElectronHandle);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  //gen Particles
  // Pruned particles are the one containing "important" stuff
  edm::Handle<edm::View<reco::GenParticle>> pruned;
  iEvent.getByToken(prunedGenToken_, pruned);
  //some  cross checks
  if (!theB.isValid()) {
    edm::LogWarning("miniAODmmmm") << "no Transient Track in event";
    return;
  }

  if (!thePATElectronHandle.isValid()) {
    edm::LogWarning("miniAODmmmm") << "no pat::Electrons in event";
    return;
  }
  if (!thePATMuonHandle.isValid()) {
    edm::LogWarning("miniAODmmmm") << "no pat::Muons in event";
    return;
  }
  if (!triggerBits.isValid()) {
    edm::LogWarning("miniAODmmmm") << "no Trigger path in event";
    return;
  }
  if (!triggerObjects.isValid()) {
    edm::LogWarning("miniAODmmmm") << "no Trigger Object in event";
    return;
  }

  //Before Begining lets get Gen level Information
  int h = 0;
  float EET_pt[5] = {-999, -999, -999, -999, -999};
  float EET_eta[5] = {-999, -999, -999, -999, -999};
  float EET_phi[5] = {-999, -999, -999, -999, -999};

  bool firedMuTrigger = false;

  //First Work on Trigger Information//
  //**********************************

  const edm::TriggerNames& names = iEvent.triggerNames(*triggerBits);
  //std::cout << "\n === TRIGGER PATHS === " << std::endl;

  // Loop though the trigger bits
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    // until we find a fired Muon
    if (!firedMuTrigger) {
      if (names.triggerName(i).find(MuonTriggerString.c_str()) != string::npos) {
        if (triggerBits->accept(i))
          firedMuTrigger = true;  // Muon trigger was fired
      }
    }
  }

    //****************************************************
  //***************Trigger Object***********************
  //****************************************************

  // if any of the triggers fired, check the trigger objects
  if (firedMuTrigger) {
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {  // note: not "const &" since we want to call unpackPathNames
      obj.unpackPathNames(names);

      if (firedMuTrigger > 0) {
        if (obj.hasPathName(MuonTriggerString.c_str(), true, true) > 0) {
          EET_pt[h] = obj.pt();
          EET_eta[h] = obj.eta();
          EET_phi[h] = obj.phi();
          h++;
        }
      }
    }
  }


  //*********************************
  //Now we get the primary vertex
  //*********************************

  reco::Vertex bestVtx;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  bestVtx = *(primaryVertices_handle->begin());

  //nVtx = primaryVertices_handle->size();

  //*****************************************
  //Let's begin by looking for J/psi

  //unsigned int nMu_tmp = thePATMuonHandle->size();

  for (View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
    for (View<pat::Muon>::const_iterator iMuon2 = iMuon1 + 1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
          //make sure all muons are diferent
          if (iMuon1 == iMuon2)
            continue;
          
          //opposite charge
          //only look for neutral charge combination

          if (!(abs(iMuon1->charge()) == 1))
            continue;
          if (!(abs(iMuon2->charge()) == 1))
            continue;

          if (!((iMuon1->charge()) + (iMuon2->charge())) == 0)
            continue;
          if (iMuon1->pt() < 2.0)
            continue;
          if (iMuon2->pt() < 2.0)
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
          //initialize 4 lepton mass
          float mu_mass = 0.1056583745;  //[PDG mass]
          //float ele_mass =  0.000510998928;//PDG mass

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
          //now M1 is first positive muon 4 vector
          //then M2 is first negative muon 4 vector

          //****************************************************************************************
          //make netral dimuon combination
          MM1 = M1 + M2;
          //****************************************************************************************

          //cout<<"Start looking muon track quality"<<endl;
          if (!(glbTrackM1->quality(reco::TrackBase::highPurity)))
            continue;
          if (!(glbTrackP1->quality(reco::TrackBase::highPurity)))
            continue;

          reco::TransientTrack muon1TT((*theB).build(glbTrackP1));
          reco::TransientTrack muon2TT((*theB).build(glbTrackM1));

          //Kalman Vtx----------------------//
          vector<TransientTrack> mu_tks;

          KalmanVertexFitter kvfM(true);
          //First neutral combination
          mu_tks.clear();
          mu_tks.push_back(muon1TT);
          mu_tks.push_back(muon2TT);
          TransientVertex J_candi1 = kvfM.vertex(mu_tks);

          if (!J_candi1.isValid()) {
            cout << "continue because no vertexed dimuon" << endl;
            continue;
          }

          reco::Vertex JPsi_Vtx1 = J_candi1;

          float B_Prob_tmp1 = TMath::Prob(J_candi1.totalChiSquared(), J_candi1.degreesOfFreedom());
          const math::XYZTLorentzVectorD JPsi_mom1 = JPsi_Vtx1.p4(mu_mass, 0.0);

          if (B_Prob_tmp1 < 0.001) {
            continue;
          }

          // Remove events with mass outside J/Psi or Z mass window
          if ((MM1.M() < 2.6 || MM1.M() > 3.6)){
            continue;
          }

          //****************************************************************************************
          //Event Information
          Run->push_back(iEvent.id().run());
          LumiBlock->push_back(iEvent.luminosityBlock());
          Event->push_back(iEvent.id().event());

          Mu_TriggerPath->push_back(firedMuTrigger);

          B_U_TriggerPt1->push_back(EET_pt[0]);
          B_U_TriggerEta1->push_back(EET_eta[0]);
          B_U_TriggerPhi1->push_back(EET_phi[0]);
          B_U_TriggerPt2->push_back(EET_pt[1]);
          B_U_TriggerEta2->push_back(EET_eta[1]);
          B_U_TriggerPhi2->push_back(EET_phi[1]);
          B_U_TriggerPt3->push_back(EET_pt[2]);
          B_U_TriggerEta3->push_back(EET_eta[2]);
          B_U_TriggerPhi3->push_back(EET_phi[2]);
          B_U_TriggerPt4->push_back(EET_pt[3]);
          B_U_TriggerEta4->push_back(EET_eta[3]);
          B_U_TriggerPhi4->push_back(EET_phi[3]);
          B_U_TriggerPt5->push_back(EET_pt[4]);
          B_U_TriggerEta5->push_back(EET_eta[4]);
          B_U_TriggerPhi5->push_back(EET_phi[4]);

          B_J1_mass->push_back(MM1.M());
          B_J1_px->push_back(MM1.Px());
          B_J1_py->push_back(MM1.Py());
          B_J1_pz->push_back(MM1.Pz());
          B_J1_pt->push_back(MM1.Pt());
          B_J1_eta->push_back(MM1.Eta());
          B_J1_phi->push_back(MM1.Phi());
          B_J1_rapidity->push_back(MM1.Rapidity());

          B_J1_VtxPx->push_back(JPsi_mom1.Px());
          B_J1_VtxPy->push_back(JPsi_mom1.Py());
          B_J1_VtxPz->push_back(JPsi_mom1.Pz());
          B_J1_VtxPt->push_back(JPsi_mom1.Pt());
          B_J1_VtxEta->push_back(JPsi_mom1.Eta());
          B_J1_VtxPhi->push_back(JPsi_mom1.Phi());
          B_J1_VtxRapidity->push_back(JPsi_mom1.Rapidity());
          B_J1_VtxMass->push_back(JPsi_mom1.mass());
          B_J1_PVx->push_back(JPsi_Vtx1.x());
          B_J1_PVy->push_back(JPsi_Vtx1.y());
          B_J1_PVz->push_back(JPsi_Vtx1.z());
          B_J1_PVxError->push_back(JPsi_Vtx1.xError());
          B_J1_PVyError->push_back(JPsi_Vtx1.yError());
          B_J1_PVzError->push_back(JPsi_Vtx1.zError());

          //dimuon vtx prob
          B_J1_VtxProb->push_back(B_Prob_tmp1);

          //new branch defn for muons
          B_Mu1_px->push_back(iMuon1->px());
          B_Mu1_py->push_back(iMuon1->py());
          B_Mu1_pz->push_back(iMuon1->pz());
          B_Mu1_pt->push_back(iMuon1->pt());
          B_Mu1_eta->push_back(iMuon1->eta());
          B_Mu1_phi->push_back(iMuon1->phi());
          B_Mu1_charge->push_back(iMuon1->charge());
          B_Mu1_soft->push_back(iMuon1->isSoftMuon(bestVtx));
          B_Mu1_tight->push_back(iMuon1->isTightMuon(bestVtx));
          B_Mu1_loose->push_back(muon::isLooseMuon(*iMuon1));
          B_Mu1_IsoTrack->push_back(iMuon1->trackIso());
          B_Mu1_IsoEcal->push_back(iMuon1->ecalIso());
          B_Mu1_IsoHcal->push_back(iMuon1->hcalIso());
          B_Mu1_IsoCalo->push_back(iMuon1->caloIso());

          B_Mu1_PaperIsoTrackRF04->push_back(
              (iMuon1->pfIsolationR04().sumChargedHadronPt +
               std::max(
                   0., iMuon1->pfIsolationR04().sumNeutralHadronEt + iMuon1->pfIsolationR04().sumPhotonEt - iMuon1->pfIsolationR04().sumPUPt * 0.5)) /
              iMuon1->pt());
          B_Mu1_PaperIsoTrackRF03->push_back(
              (iMuon1->pfIsolationR03().sumChargedHadronPt +
               std::max(
                   0., iMuon1->pfIsolationR03().sumNeutralHadronEt + iMuon1->pfIsolationR03().sumPhotonEt - iMuon1->pfIsolationR03().sumPUPt * 0.5)) /
              iMuon1->pt());

          B_Mu1_Paper3DIP->push_back(iMuon1->dB(pat::Muon::PV3D) / iMuon1->edB(pat::Muon::PV3D));

          B_J_xyP1->push_back(glbTrackP1->dxy(bestVtx.position()));
          B_J_xyM1->push_back(glbTrackM1->dxy(bestVtx.position()));
          B_J_zP1->push_back(glbTrackM1->dz(bestVtx.position()));
          B_J_zM1->push_back(glbTrackP1->dz(bestVtx.position())); 

          //cout<<"End of all loop"<<endl;

          mu1mC2->push_back(glbTrackM1->normalizedChi2());
          //mumAngT->push_back( muon::isGoodMuon(*iMuon1,muon::TMLastStationAngTight) ); //
          mu1mNHits->push_back(glbTrackM1->numberOfValidHits());
          mu1mNPHits->push_back(glbTrackM1->hitPattern().numberOfValidPixelHits());
          mu1pC2->push_back(glbTrackP1->normalizedChi2());
          //mupAngT->push_back( muon::isGoodMuon(*iMuon2,muon::TMLastStationAngTight) );  //
          mu1pNHits->push_back(glbTrackP1->numberOfValidHits());
          mu1pNPHits->push_back(glbTrackP1->hitPattern().numberOfValidPixelHits());

          B_M1_pt->push_back(M1.Pt());
          B_M1_eta->push_back(M1.Eta());
          B_M1_phi->push_back(M1.Phi());
          B_M1_px->push_back(M1.Px());
          B_M1_py->push_back(M1.Py());
          B_M1_pz->push_back(M1.Pz());
          nB++;
    }
  }

  if (nB > 0) {
    //std::cout << "filling tree" << endl;
    tree_->Fill();
  }

  nB = 0;
  Run->clear();
  LumiBlock->clear();
  Event->clear();

  //B_Z_dca->clear();
  Mu_TriggerPath->clear();
  B_U_TriggerPt1->clear();
  B_U_TriggerEta1->clear();
  B_U_TriggerPhi1->clear();
  B_U_TriggerPt2->clear();
  B_U_TriggerEta2->clear();
  B_U_TriggerPhi2->clear();
  B_U_TriggerPt3->clear();
  B_U_TriggerEta3->clear();
  B_U_TriggerPhi3->clear();
  B_U_TriggerPt4->clear();
  B_U_TriggerEta4->clear();
  B_U_TriggerPhi4->clear();
  B_U_TriggerPt5->clear();
  B_U_TriggerEta5->clear();
  B_U_TriggerPhi5->clear();

  B_J1_mass->clear();
  B_J1_px->clear();
  B_J1_py->clear();
  B_J1_pz->clear();
  B_J1_pt->clear();
  B_J1_eta->clear();
  B_J1_phi->clear();
  B_J1_rapidity->clear();

  B_J1_VtxPx->clear();
  B_J1_VtxPy->clear();
  B_J1_VtxPz->clear();
  B_J1_VtxPt->clear();
  B_J1_VtxEta->clear();
  B_J1_VtxPhi->clear();
  B_J1_VtxRapidity->clear();
  B_J1_VtxMass->clear();
  B_J1_PVx->clear();
  B_J1_PVy->clear();
  B_J1_PVz->clear();
  B_J1_PVxError->clear();
  B_J1_PVyError->clear();
  B_J1_PVzError->clear();

  B_Mu1_px->clear();
  B_Mu1_py->clear();
  B_Mu1_pz->clear();
  B_Mu1_charge->clear();
  B_Mu1_pt->clear();
  B_Mu1_eta->clear();
  B_Mu1_phi->clear();
  B_Mu1_soft->clear();
  B_Mu1_tight->clear();
  B_Mu1_loose->clear();
  B_Mu1_IsoTrack->clear();
  B_Mu1_IsoHcal->clear();
  B_Mu1_IsoEcal->clear();
  B_Mu1_IsoCalo->clear();
  B_J1_VtxProb->clear();
  B_J_xyP1->clear();
  B_J_xyM1->clear();
  B_J_zP1->clear();
  B_J_zM1->clear();
  B_Mu1_PaperIsoTrackRF04->clear();
  B_Mu1_PaperIsoTrackRF03->clear();

  B_Mu1_Paper3DIP->clear();

  mu1mC2->clear();
  mu1mNHits->clear();
  mu1mNPHits->clear();
  mu1pC2->clear();
  mu1pNHits->clear();
  mu1pNPHits->clear();

  B_M1_pt->clear();
  B_M1_eta->clear();
  B_M1_phi->clear();
  B_M1_px->clear();
  B_M1_px->clear();
  B_M1_pz->clear();
}

// ------------ method called once each job just before starting event loop  ------------

void miniAODmmmm::beginJob() {
  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  //edm::Service<TFileService> fs;
  //tree_ = fs->make<TTree>("ntuple"," J/psi ntuple");

  //tree_->Branch("nB",&nB,"nB/i");
  tree_ = new TTree("ntuple", "ntuple");

  tree_->Branch("nB", &nB, "nB/i");
  //electron channels
  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  tree_->Branch("Event", &Event);

  tree_->Branch("Mu_TriggerPath", &Mu_TriggerPath);
  tree_->Branch("B_U_TriggerPt1", &B_U_TriggerPt1);
  tree_->Branch("B_U_TriggerEta1", &B_U_TriggerEta1);
  tree_->Branch("B_U_TriggerPhi1", &B_U_TriggerPhi1);
  tree_->Branch("B_U_TriggerPt2", &B_U_TriggerPt2);
  tree_->Branch("B_U_TriggerEta2", &B_U_TriggerEta2);
  tree_->Branch("B_U_TriggerPhi2", &B_U_TriggerPhi2);
  tree_->Branch("B_U_TriggerPt3", &B_U_TriggerPt3);
  tree_->Branch("B_U_TriggerEta3", &B_U_TriggerEta3);
  tree_->Branch("B_U_TriggerPhi3", &B_U_TriggerPhi3);
  tree_->Branch("B_U_TriggerPt4", &B_U_TriggerPt4);
  tree_->Branch("B_U_TriggerEta4", &B_U_TriggerEta4);
  tree_->Branch("B_U_TriggerPhi4", &B_U_TriggerPhi4);
  tree_->Branch("B_U_TriggerPt5", &B_U_TriggerPt5);
  tree_->Branch("B_U_TriggerEta5", &B_U_TriggerEta5);
  tree_->Branch("B_U_TriggerPhi5", &B_U_TriggerPhi5);

  tree_->Branch("B_J1_mass", &B_J1_mass);
  tree_->Branch("B_J1_px", &B_J1_px);
  tree_->Branch("B_J1_py", &B_J1_py);
  tree_->Branch("B_J1_pz", &B_J1_pz);
  tree_->Branch("B_J1_pt", &B_J1_pt);
  tree_->Branch("B_J1_eta", &B_J1_eta);
  tree_->Branch("B_J1_phi", &B_J1_phi);
  tree_->Branch("B_J1_rapidity", &B_J1_rapidity);

  tree_->Branch("B_J1_VtxPx", &B_J1_VtxPx);
  tree_->Branch("B_J1_VtxPy", &B_J1_VtxPy);
  tree_->Branch("B_J1_VtxPz", &B_J1_VtxPz);
  tree_->Branch("B_J1_VtxPt", &B_J1_VtxPt);
  tree_->Branch("B_J1_VtxEta", &B_J1_VtxEta);
  tree_->Branch("B_J1_VtxPhi", &B_J1_VtxPhi);
  tree_->Branch("B_J1_VtxRapidity", &B_J1_VtxRapidity);
  tree_->Branch("B_J1_VtxMass", &B_J1_VtxMass);
  tree_->Branch("B_J1_PVx", &B_J1_PVx);
  tree_->Branch("B_J1_PVy", &B_J1_PVy);
  tree_->Branch("B_J1_PVz", &B_J1_PVz);
  tree_->Branch("B_J1_PVxError", &B_J1_PVxError);
  tree_->Branch("B_J1_PVyError", &B_J1_PVyError);
  tree_->Branch("B_J1_PVzError", &B_J1_PVzError);

  tree_->Branch("B_Mu1_px", &B_Mu1_px);
  tree_->Branch("B_Mu1_py", &B_Mu1_py);
  tree_->Branch("B_Mu1_pz", &B_Mu1_pz);
  tree_->Branch("B_Mu1_pt", &B_Mu1_pt);
  tree_->Branch("B_Mu1_eta", &B_Mu1_eta);
  tree_->Branch("B_Mu1_phi", &B_Mu1_phi);
  tree_->Branch("B_Mu1_charge", &B_Mu1_charge);
  tree_->Branch("B_Mu1_soft", &B_Mu1_soft);
  tree_->Branch("B_Mu1_tight", &B_Mu1_tight);
  tree_->Branch("B_Mu1_loose", &B_Mu1_loose);
  tree_->Branch("B_Mu1_IsoTrack", &B_Mu1_IsoTrack);
  tree_->Branch("B_Mu1_IsoHcal", &B_Mu1_IsoHcal);
  tree_->Branch("B_Mu1_IsoEcal", &B_Mu1_IsoEcal);
  tree_->Branch("B_Mu1_IsoCalo", &B_Mu1_IsoCalo);

  tree_->Branch("B_Mu1_PaperIsoTrackRF03", &B_Mu1_PaperIsoTrackRF03);
  tree_->Branch("B_Mu1_PaperIsoTrackRF04", &B_Mu1_PaperIsoTrackRF04);

  tree_->Branch("B_Mu1_Paper3DIP", &B_Mu1_Paper3DIP);

  tree_->Branch("B_J1_VtxProb", &B_J1_VtxProb);

  tree_->Branch("B_J_xyP1", &B_J_xyP1);
  tree_->Branch("B_J_xyM1", &B_J_xyM1);
  tree_->Branch("B_J_zP1", &B_J_zP1);
  tree_->Branch("B_J_zM1", &B_J_zM1);

  tree_->Branch("B_J_xyP2", &B_J_xyP2);
  tree_->Branch("B_J_xyM2", &B_J_xyM2);
  tree_->Branch("B_J_zP2", &B_J_zP2);
  tree_->Branch("B_J_zM2", &B_J_zM2);

  tree_->Branch("mu1mC2", &mu1mC2);
  tree_->Branch("mu1mNHits", &mu1mNHits);
  tree_->Branch("mu1mNPHits", &mu1mNPHits);
  tree_->Branch("mu1pC2", &mu1pC2);
  tree_->Branch("mu1pNHits", &mu1pNHits);
  tree_->Branch("mu1pNPHits", &mu1pNPHits);

  tree_->Branch("B_M1_pt", &B_M1_pt);
  tree_->Branch("B_M1_eta", &B_M1_eta);
  tree_->Branch("B_M1_phi", &B_M1_phi);
  tree_->Branch("B_M1_px", &B_M1_px);
  tree_->Branch("B_M1_py", &B_M1_py);
  tree_->Branch("B_M1_pz", &B_M1_pz);


  tree_->Branch("B_J_GenMuonPt", &B_J_GenMuonPt);
  tree_->Branch("B_J_GenMuonEta", &B_J_GenMuonEta);
  tree_->Branch("B_J_GenMuonPhi", &B_J_GenMuonPhi);
  tree_->Branch("B_Z_GenMuonPt", &B_Z_GenMuonPt);
  tree_->Branch("B_Z_GenMuonEta", &B_Z_GenMuonEta);
  tree_->Branch("B_Z_GenMuonPhi", &B_Z_GenMuonPhi);
}

// ------------ method called once each job just after ending the event loop  ------------
void miniAODmmmm::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(miniAODmmmm);
