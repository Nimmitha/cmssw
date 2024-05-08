// -*- C++ -*-
//
// Package:    miniAODmuons
// Class:      miniAODmuons
// 

//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Monday Aug 28 (2017)         |
//         <jhovanny.andres.mejia.guisao@cern.ch> |
// updated by: Nimmitha Karunarathna              |
//         on:  04/18/2024                        |
//=================================================

// system include files
#include <memory>

#include "myAnalyzers/JPsiKsPAT/src/miniAODmuons.h"

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

miniAODmuons::miniAODmuons(const edm::ParameterSet& iConfig)
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
      ElectronTriggerString(iConfig.getParameter<std::string>("ElectronTrigger")),
      DataTypeString(iConfig.getParameter<std::string>("DataType")),
      isMC_(iConfig.getParameter<bool>("isMC")),

      tree_(0),
      //Gen add
      Y_GenMuonPt(0),
      Y_GenMuonEta(0),
      Y_GenMuonPhi(0),
      Z_GenMuonPt(0),
      Z_GenMuonEta(0),
      Z_GenMuonPhi(0),

      Run(0),
      LumiBlock(0),
      Event(0),
      FourL_mass(0),
      FourL_px(0),
      FourL_py(0),
      FourL_pz(0),
      FourL_pt(0),
      FourL_eta(0),
      FourL_phi(0),
      FourL_rapidity(0),

      FourL_Vtx_Prob(0),
      FourL_Vtx_x(0),
      FourL_Vtx_y(0),
      FourL_Vtx_z(0),
      FourL_Vtx_xError(0),
      FourL_Vtx_yError(0),
      FourL_Vtx_zError(0),

      FourL_Vtx_Px(0),
      FourL_Vtx_Py(0),
      FourL_Vtx_Pz(0),
      FourL_Vtx_Pt(0),
      FourL_Vtx_Eta(0),
      FourL_Vtx_Phi(0),
      FourL_Vtx_Rapidity(0),
      FourL_Vtx_Mass(0),

      //Z
      Mu_TriggerPath(0),
      Ele_TriggerPath(0),
      Mu_TriggerPt1(0),
      Mu_TriggerEta1(0),
      Mu_TriggerPhi1(0),
      Mu_TriggerPt2(0),
      Mu_TriggerEta2(0),
      Mu_TriggerPhi2(0),
      Mu_TriggerPt3(0),
      Mu_TriggerEta3(0),
      Mu_TriggerPhi3(0),
      Mu_TriggerPt4(0),
      Mu_TriggerEta4(0),
      Mu_TriggerPhi4(0),
      Mu_TriggerPt5(0),
      Mu_TriggerEta5(0),
      Mu_TriggerPhi5(0),

      Ele_TriggerPt1(0),
      Ele_TriggerEta1(0),
      Ele_TriggerPhi1(0),
      Ele_TriggerPt2(0),
      Ele_TriggerEta2(0),
      Ele_TriggerPhi2(0),
      Ele_TriggerPt3(0),
      Ele_TriggerEta3(0),
      Ele_TriggerPhi3(0),
      Ele_TriggerPt4(0),
      Ele_TriggerEta4(0),
      Ele_TriggerPhi4(0),
      Ele_TriggerPt5(0),
      Ele_TriggerEta5(0),
      Ele_TriggerPhi5(0),

      Y_mass(0),
      Y_Vtx_Prob(0),
      Y_px(0),
      Y_py(0),
      Y_pz(0),
      Y_pt(0),
      Y_eta(0),
      Y_phi(0),
      Y_rapidity(0),
      Y_Vtx_Px(0),
      Y_Vtx_Py(0),
      Y_Vtx_Pz(0),
      Y_Vtx_Pt(0),
      Y_Vtx_Eta(0),
      Y_Vtx_Phi(0),
      Y_Vtx_Rapidity(0),
      Y_Vtx_Mass(0),

      Y_Vtx_x(0),
      Y_Vtx_y(0),
      Y_Vtx_z(0),
      Y_Vtx_xError(0),
      Y_Vtx_yError(0),
      Y_Vtx_zError(0),
      Y_px1(0),
      Y_py1(0),
      Y_pz1(0),
      Y_pt1(0),
      Y_eta1(0),
      Y_SCeta1(0),
      Y_phi1(0),
      Y_energy1(0),
      Y_energyCorr1(0),
      Y_ecalIso1(0),
      Y_hcalIso1(0),
      Y_trackIso1(0),
      Z_trackIso1(0),
      Y_CutBaseLoose1(0),
      Y_CutBaseVeto1(0),
      Y_mvaIsoWP90_1(0),
      Y_mvaIsoWP80_1(0),
      Y_px2(0),
      Y_py2(0),
      Y_pz2(0),
      Y_pt2(0),
      Y_eta2(0),
      Y_SCeta2(0),
      Y_phi2(0),
      Y_energy2(0),
      Y_energyCorr2(0),
      Y_ecalIso2(0),
      Y_hcalIso2(0),
      Y_trackIso2(0),
      Z_trackIso2(0),
      Y_CutBaseLoose2(0),
      Y_CutBaseVeto2(0),
      Y_mvaIsoWP90_2(0),
      Y_mvaIsoWP80_2(0),
      Y_charge1(0),
      Y_charge2(0),
      Y_fit_pt1(0),
      Y_fit_ptError1(0),
      Y_fit_eta1(0),
      Y_fit_phi1(0),
      Y_fit_pt2(0),
      Y_fit_ptError2(0),
      Y_fit_eta2(0),
      Y_fit_phi2(0),
      Y_dxy1(0),
      Y_dxy2(0),
      Y_dz1(0),
      Y_dz2(0),

      Y_lowPt(0),
      Y_highPt(0),
      Z_lowPt(0),
      Z_highPt(0),
      Mu_dCA(0),
      Ele_dCA(0),
      Z_mass(0),
      Z_px(0),
      Z_py(0),
      Z_pz(0),
      Z_pt(0),
      Z_eta(0),
      Z_phi(0),
      Z_rapidity(0),
      Z_Vtx_Px(0),
      Z_Vtx_Py(0),
      Z_Vtx_Pz(0),
      Z_Vtx_Pt(0),
      Z_Vtx_Eta(0),
      Z_Vtx_Phi(0),
      Z_Vtx_Rapidity(0),
      Z_Vtx_Mass(0),
      Z_Vtx_x(0),
      Z_Vtx_y(0),
      Z_Vtx_z(0),
      Z_Vtx_xError(0),
      Z_Vtx_yError(0),
      Z_Vtx_zError(0),
      Z_px1(0),
      Z_py1(0),
      Z_pz1(0),
      Z_pt1(0),
      Z_eta1(0),
      Z_phi1(0),
      Z_soft1(0),
      Z_tight1(0),
      Z_loose1(0),
      Z_fit_pt1(0),
      Z_fit_ptError1(0),
      Z_fit_eta1(0),
      Z_fit_phi1(0),
      Z_px2(0),
      Z_py2(0),
      Z_pz2(0),
      Z_pt2(0),
      Z_eta2(0),
      Z_phi2(0),
      Z_charge1(0),
      Z_charge2(0),
      Z_soft2(0),
      Z_tight2(0),
      Z_loose2(0),
      Z_fit_pt2(0),
      Z_fit_ptError2(0),
      Z_fit_eta2(0),
      Z_fit_phi2(0),
      Z_Vtx_Prob(0),
      Z_xy1(0),
      Z_xy2(0),
      Z_z2(0),
      Z_z1(0),

      Mu1_C2(0),
      Mu1_NHits(0),
      Mu1_NPHits(0),
      Mu2_C2(0),
      Mu2_NHits(0),
      Mu2_NPHits(0),

      nCandi(0)

{
  //now do what ever initialization is needed
}

miniAODmuons::~miniAODmuons() {}

//
// member functions
//

// ------------ method called to for each event  ------------
void miniAODmuons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  //*********************************
  // Get event content informationF
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
    edm::LogWarning("miniAODmuons") << "no Transient Track in event";
    return;
  }

  if (!thePATElectronHandle.isValid()) {
    edm::LogWarning("miniAODmuons") << "no pat::Electrons in event";
    return;
  }
  if (!thePATMuonHandle.isValid()) {
    edm::LogWarning("miniAODmuons") << "no pat::Muons in event";
    return;
  }
  if (!triggerBits.isValid()) {
    edm::LogWarning("miniAODmuons") << "no Trigger path in event";
    return;
  }
  if (!triggerObjects.isValid()) {
    edm::LogWarning("miniAODmuons") << "no Trigger Object in event";
    return;
  }

  //Before Begining lets get Gen level Information
  int h = 0;
  int h1 = 0;

  float EET_pt[5] = {-999, -999, -999, -999, -999};
  float EET_eta[5] = {-999, -999, -999, -999, -999};
  float EET_phi[5] = {-999, -999, -999, -999, -999};
  float EET_pt1[5] = {-999, -999, -999, -999, -999};
  float EET_eta1[5] = {-999, -999, -999, -999, -999};
  float EET_phi1[5] = {-999, -999, -999, -999, -999};
  int nG1 = 0;
  int nG2 = 0;

  bool firedMuTrigger = false;
  bool firedEleTrig = false;

  bool GenInfo = false;
  if (GenInfo) {
    //Gen Level Info

    float Z_GenMuon_pt = -999;
    float Z_GenMuon_eta = -999;
    float Z_GenMuon_phi = -999;
    float Y_GenMuon_pt = -999;
    float Y_GenMuon_eta = -999;
    float Y_GenMuon_phi = -999;

    for (size_t i = 0; i < pruned->size(); i++) {
      //if( (*pruned)[i].isPromptFinalState()  && abs((*pruned)[i].pdgId() ==13) ){

      if (abs((*pruned)[i].pdgId()) == 13) {
        //cout<<"Found gen level muon"<<endl;
        if ((*pruned)[i].mother()->pdgId() == 553) {
          //cout<<"Found gen level muon from Jpsi "<<endl;
          //if ( (*pruned)[i].->pdgId()==10443 ) {
          //cout<<"Found gen level muon with number of mother"<<(*pruned)[i].numberOfMothers()<<endl;
          Z_GenMuon_pt = (*pruned)[i].pt();
          Z_GenMuon_eta = (*pruned)[i].eta();
          Z_GenMuon_phi = (*pruned)[i].phi();
          Z_GenMuonPt->push_back(Z_GenMuon_pt);
          Z_GenMuonEta->push_back(Z_GenMuon_eta);
          Z_GenMuonPhi->push_back(Z_GenMuon_phi);
          nG1++;
          //}
        }
      }
      if (abs((*pruned)[i].pdgId()) == 11) {
        //cout<<"Found gen level muon"<<endl;
        if ((*pruned)[i].mother()->pdgId() == 23) {
          //cout<<"Found gen level muon from Z"<<endl;
          Y_GenMuon_pt = (*pruned)[i].pt();
          Y_GenMuon_eta = (*pruned)[i].eta();
          Y_GenMuon_phi = (*pruned)[i].phi();
          Y_GenMuonPt->push_back(Y_GenMuon_pt);
          Y_GenMuonEta->push_back(Y_GenMuon_eta);
          Y_GenMuonPhi->push_back(Y_GenMuon_phi);
          nG2++;
        }
      }
    }
    // cout<<"End Gen Work for this event"<<endl;
  }

  //First Work on Trigger Information//

  //****************************************************
  //***************Trigger Path************************
  //****************************************************
  const edm::TriggerNames& names = iEvent.triggerNames(*triggerBits);

  // Loop though the trigger bits
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    // until we find a fired Muon and Electron triggers
    if (!firedMuTrigger || !firedEleTrig) {
      if (names.triggerName(i).find(MuonTriggerString.c_str()) != string::npos) {
        if (triggerBits->accept(i))
          firedMuTrigger = true;  // Muon trigger was fired
      }

      if (names.triggerName(i).find(ElectronTriggerString.c_str()) != string::npos) {
        if (triggerBits->accept(i))
          firedEleTrig = true;  // Electron trigger was fired
      }
    }
  }

  //****************************************************
  //***************Trigger Object***********************
  //****************************************************

  // if any of the triggers fired, check the trigger objects
  if (firedMuTrigger || firedEleTrig) {
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

      if (firedEleTrig > 0) {
        if (obj.hasPathName(ElectronTriggerString.c_str(), true, true) > 0) {
          EET_pt1[h1] = obj.pt();
          EET_eta1[h1] = obj.eta();
          EET_phi1[h1] = obj.phi();
          h1++;
        }
      }
    }
  }

  //****************************************************
  //*********Now we get the primary vertex**************
  //****************************************************
  reco::Vertex bestVtx;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  bestVtx = *(primaryVertices_handle->begin());

  //****************************************************
  //*******Identify the muon and electron pairs*********
  //****************************************************
  for (View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
    for (View<pat::Muon>::const_iterator iMuon2 = iMuon1 + 1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
      if ((iMuon1->charge()) * (iMuon2->charge()) == 1)
        continue;

      if (!(abs(iMuon1->charge()) == 1) || !(abs(iMuon2->charge()) == 1))
        continue;

      TrackRef glbTrack1;
      TrackRef glbTrack2;

      glbTrack1 = iMuon1->track();
      glbTrack2 = iMuon2->track();

      if (glbTrack1.isNull() || glbTrack2.isNull())
        continue;

      if (iMuon1->track()->pt() < 2.0 || iMuon2->track()->pt() < 2.0)
        continue;

      if (!(glbTrack1->quality(reco::TrackBase::highPurity)))
        continue;

      if (!(glbTrack2->quality(reco::TrackBase::highPurity)))
        continue;

      for (View<pat::Electron>::const_iterator iEle1 = thePATElectronHandle->begin(); iEle1 != thePATElectronHandle->end(); ++iEle1) {
        for (View<pat::Electron>::const_iterator iEle2 = iEle1 + 1; iEle2 != thePATElectronHandle->end(); ++iEle2) {
          if ((iEle1->charge()) * (iEle2->charge()) == 1)
            continue;

          // Electron charge check
          if (!(abs(iEle1->charge()) == 1) || !(abs(iEle2->charge()) == 1))
            continue;

          if (!iEle1->gsfTrack().isAvailable() || !iEle1->gsfTrack().isNonnull())
            continue;

          if (!iEle2->gsfTrack().isAvailable() || !iEle2->gsfTrack().isNonnull())
            continue;

          TLorentzVector M1, M2, E1, E2, MM, EE, EEMM;
          float mu_mass = 0.1056583745;     //PDG mass
          float ele_mass = 0.000510998928;  //PDG mass

          M1.SetXYZM(iMuon1->track()->px(), iMuon1->track()->py(), iMuon1->track()->pz(), mu_mass);
          M2.SetXYZM(iMuon2->track()->px(), iMuon2->track()->py(), iMuon2->track()->pz(), mu_mass);

          E1.SetPtEtaPhiE(iEle1->pt(), iEle1->eta(), iEle1->phi(), iEle1->energy());
          E2.SetPtEtaPhiE(iEle2->pt(), iEle2->eta(), iEle2->phi(), iEle2->energy());
          //energy scale smeared
          // E1.SetPtEtaPhiE(iEle1->pt(),iEle1->eta(),iEle1->phi(),iEle1->userFloat("ecalTrkEnergyPostCorr"))
          // E2.SetPtEtaPhiE(iEle2->pt(),iEle2->eta(),iEle2->phi(),iEle2->userFloat("ecalTrkEnergyPostCorr"));

          MM = M1 + M2;
          EE = E1 + E2;
          EEMM = MM + EE;

          if (MM.M() < 50 || MM.M() > 130)
            continue;
          if (EE.M() < 0 || EE.M() > 12)
            continue;

          // Start building the tracks
          reco::TransientTrack muon1_TT((*theB).build(glbTrack1));
          reco::TransientTrack muon2_TT((*theB).build(glbTrack2));
          //reco::TransientTrack electron1_TT((*theB).build(kfTrackRefP));//working
          //reco::TransientTrack electron2_TT((*theB).build(kfTrackRefM));//working
          reco::TransientTrack electron1_TT((*theB).build(iEle1->gsfTrack()));
          reco::TransientTrack electron2_TT((*theB).build(iEle2->gsfTrack()));

          // Validating the impact points
          if (!muon1_TT.impactPointTSCP().isValid() || !muon2_TT.impactPointTSCP().isValid())
            continue;
          if (!electron1_TT.impactPointTSCP().isValid() || !electron2_TT.impactPointTSCP().isValid())
            continue;

          // Trajectory states to calculate DCA for the 2 muons
          // FreeTrajectoryState mu1_State = muon1_TT.impactPointTSCP().theState();
          // FreeTrajectoryState mu2_State = muon2_TT.impactPointTSCP().theState();
          // FreeTrajectoryState electron1_State = muon1_TT.impactPointTSCP().theState();
          // FreeTrajectoryState electron2_State = muon1_TT.impactPointTSCP().theState();

          // Measure distance between tracks at their closest approach
          // ClosestApproachInRPhi cAppMu;
          // cAppMu.calculate(mu1_State, mu2_State);
          // if (!cAppMu.status())
          //   continue;

          float MuCA = -1;
          // float MuCA = fabs(cAppMu.distance());
          //if (MuCA < 0. || MuCA > 0.5) continue;
          //cout<<"Closest approach  Mu "<<MuCA<<endl;

          // ClosestApproachInRPhi cAppEle;
          // cAppEle.calculate(electron1_State, electron2_State);
          // if (!cAppEle.status()){
          //   continue;
          // }

          float EleCA = -1;
          // float EleCA = fabs(cAppEle.distance());
          // if (EleCA < 0. || EleCA > 0.5) continue;
          // cout << "Closest approach  Ele " << EleCA << endl;

          //************************************************************
          //******************Kalman Vertex Fitting*********************
          //************************************************************

          vector<TransientTrack> ele_tks;
          vector<TransientTrack> mu_tks;
          vector<TransientTrack> mmee_tks;

          // Dielectron Vertex Fit
          KalmanVertexFitter kvf(true);
          ele_tks.clear();
          ele_tks.push_back(electron1_TT);
          ele_tks.push_back(electron2_TT);

          TransientVertex Y_candi = kvf.vertex(ele_tks);
          if (!Y_candi.isValid()) {
            //cout<<"Z candidate non validated by kalman fitter"<<endl;
            continue;
          }
          float tmpY_Prob = TMath::Prob(Y_candi.totalChiSquared(), Y_candi.degreesOfFreedom());

          reco::Vertex Y_Vtx = Y_candi;
          const math::XYZTLorentzVectorD Y_mom = Y_Vtx.p4(ele_mass, 0.0);

          // Get the refitted tracks
          reco::TransientTrack fit_electron1_TT = Y_candi.refittedTrack(electron1_TT);
          reco::TransientTrack fit_electron2_TT = Y_candi.refittedTrack(electron2_TT);

          if (!fit_electron1_TT.isValid() || !fit_electron2_TT.isValid())
            continue;

          // Dimuon Vertex Fit
          KalmanVertexFitter kvfM(true);
          mu_tks.clear();
          mu_tks.push_back(muon1_TT);
          mu_tks.push_back(muon2_TT);

          TransientVertex Z_candi = kvfM.vertex(mu_tks);
          if (!Z_candi.isValid()) {
            //cout<<"Z candidate non validated by kalman fitter"<<endl;
            continue;
          }
          float tmpZ_Prob = TMath::Prob(Z_candi.totalChiSquared(), Z_candi.degreesOfFreedom());

          reco::Vertex Z_Vtx = Z_candi;
          const math::XYZTLorentzVectorD Z_mom = Z_Vtx.p4(mu_mass, 0.0);

          // Get the refitted tracks
          reco::TransientTrack fit_muon1_TT = Z_candi.refittedTrack(muon1_TT);
          reco::TransientTrack fit_muon2_TT = Z_candi.refittedTrack(muon2_TT);

          if (!fit_muon1_TT.isValid() || !fit_muon2_TT.isValid())
            continue;

          // Four Lepton Vertex Fit
          KalmanVertexFitter kvfEM(true);
          mmee_tks.clear();
          mmee_tks.push_back(electron1_TT);
          mmee_tks.push_back(electron2_TT);
          mmee_tks.push_back(muon1_TT);
          mmee_tks.push_back(muon2_TT);

          TransientVertex FourL_candi = kvfEM.vertex(mmee_tks);
          if (!FourL_candi.isValid()) {
            //cout<<"FourL candidate non validated by kalman fitter"<<endl;
            continue;
          }
          float tmpFourL_Prob = TMath::Prob(FourL_candi.totalChiSquared(), FourL_candi.degreesOfFreedom());
          if (tmpFourL_Prob < 0.001)
            continue;

          reco::Vertex FourL_Vtx = FourL_candi;
          // const math::XYZTLorentzVectorD FourL_mom = FourL_Vtx.p4();

          // Get the refitted tracks
          // reco::TransientTrack fit_FourL_electron1_TT = FourL_candi.refittedTrack(electron1_TT);
          // reco::TransientTrack fit_FourL_electron2_TT = FourL_candi.refittedTrack(electron2_TT);
          // reco::TransientTrack fit_FourL_muon1_TT = FourL_candi.refittedTrack(muon1_TT);
          // reco::TransientTrack fit_FourL_muon2_TT = FourL_candi.refittedTrack(muon2_TT);

          // if (!fit_FourL_electron1_TT.isValid() || !fit_FourL_electron2_TT.isValid() || !fit_FourL_muon1_TT.isValid() || !fit_FourL_muon2_TT.isValid())
          //   continue;

          //******************************************************
          //******************Kinematic Fit***********************
          //******************************************************
          bool KinFit = false;
          if (KinFit == true) {
            ParticleMass muon_mass = 0.10565837;  //pdg mass
            float muon_sigma = muon_mass * 1.e-6;
            vector<RefCountedKinematicParticle> muonParticles;
            muonParticles.clear();

            //Creating a KinematicParticleFactory
            KinematicParticleFactoryFromTransientTrack pFactory;

            //initial chi2 and ndf before kinematic fits.
            float chi = 0.;
            float ndf = 0.;
            try {
              muonParticles.push_back(pFactory.particle(muon1_TT, muon_mass, chi, ndf, muon_sigma));
              muonParticles.push_back(pFactory.particle(muon2_TT, muon_mass, chi, ndf, muon_sigma));
            } catch (...) {
              std::cout << " Exception caught ... continuing 1 " << std::endl;
              continue;
            }

            KinematicParticleVertexFitter fitter;

            RefCountedKinematicTree psiVertexFitTree;
            try {
              psiVertexFitTree = fitter.fit(muonParticles);
            } catch (...) {
              std::cout << " Exception caught ... continuing 2 " << std::endl;
              continue;
            }

            if (!psiVertexFitTree->isValid()) {
              //std::cout << "caught an exception in the psi vertex fit" << std::endl;
              continue;
            }

            psiVertexFitTree->movePointerToTheTop();

            RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
            RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();

            if (psi_vFit_vertex_noMC->chiSquared() < 0) {
              std::cout << "negative chisq from psi fit" << endl;
              continue;
            }
          }

          //Event Information
          Run->push_back(iEvent.id().run());
          LumiBlock->push_back(iEvent.luminosityBlock());
          Event->push_back(iEvent.id().event());

          //FourL Information
          FourL_mass->push_back(EEMM.M());
          FourL_px->push_back(EEMM.Px());
          FourL_py->push_back(EEMM.Py());
          FourL_pz->push_back(EEMM.Pz());
          FourL_pt->push_back(EEMM.Pt());
          FourL_eta->push_back(EEMM.Eta());
          FourL_phi->push_back(EEMM.Phi());
          FourL_rapidity->push_back(EEMM.Rapidity());

          FourL_Vtx_Prob->push_back(tmpFourL_Prob);
          FourL_Vtx_x->push_back(FourL_Vtx.x());
          FourL_Vtx_y->push_back(FourL_Vtx.y());
          FourL_Vtx_z->push_back(FourL_Vtx.z());
          FourL_Vtx_xError->push_back(FourL_Vtx.xError());
          FourL_Vtx_yError->push_back(FourL_Vtx.yError());
          FourL_Vtx_zError->push_back(FourL_Vtx.zError());

          // FourL_Vtx_Px->push_back(FourL_mom.Px());
          // FourL_Vtx_Py->push_back(FourL_mom.Py());
          // FourL_Vtx_Pz->push_back(FourL_mom.Pz());
          // FourL_Vtx_Pt->push_back(FourL_mom.Pt());
          // FourL_Vtx_Eta->push_back(FourL_mom.Eta());
          // FourL_Vtx_Phi->push_back(FourL_mom.Phi());
          // FourL_Vtx_Rapidity->push_back(FourL_mom.Rapidity());
          // FourL_Vtx_Mass->push_back(FourL_mom.mass());

          // Trigger Path Information
          Mu_TriggerPath->push_back(firedMuTrigger);
          Ele_TriggerPath->push_back(firedEleTrig);

          // Trigger Object Information
          Mu_TriggerPt1->push_back(EET_pt[0]);
          Mu_TriggerEta1->push_back(EET_eta[0]);
          Mu_TriggerPhi1->push_back(EET_phi[0]);
          Mu_TriggerPt2->push_back(EET_pt[1]);
          Mu_TriggerEta2->push_back(EET_eta[1]);
          Mu_TriggerPhi2->push_back(EET_phi[1]);
          Mu_TriggerPt3->push_back(EET_pt[2]);
          Mu_TriggerEta3->push_back(EET_eta[2]);
          Mu_TriggerPhi3->push_back(EET_phi[2]);
          Mu_TriggerPt4->push_back(EET_pt[3]);
          Mu_TriggerEta4->push_back(EET_eta[3]);
          Mu_TriggerPhi4->push_back(EET_phi[3]);
          Mu_TriggerPt5->push_back(EET_pt[4]);
          Mu_TriggerEta5->push_back(EET_eta[4]);
          Mu_TriggerPhi5->push_back(EET_phi[4]);

          Ele_TriggerPt1->push_back(EET_pt1[0]);
          Ele_TriggerEta1->push_back(EET_eta1[0]);
          Ele_TriggerPhi1->push_back(EET_phi1[0]);
          Ele_TriggerPt2->push_back(EET_pt1[1]);
          Ele_TriggerEta2->push_back(EET_eta1[1]);
          Ele_TriggerPhi2->push_back(EET_phi1[1]);
          Ele_TriggerPt3->push_back(EET_pt1[2]);
          Ele_TriggerEta3->push_back(EET_eta1[2]);
          Ele_TriggerPhi3->push_back(EET_phi1[2]);
          Ele_TriggerPt4->push_back(EET_pt1[3]);
          Ele_TriggerEta4->push_back(EET_eta1[3]);
          Ele_TriggerPhi4->push_back(EET_phi1[3]);
          Ele_TriggerPt5->push_back(EET_pt1[4]);
          Ele_TriggerEta5->push_back(EET_eta1[4]);
          Ele_TriggerPhi5->push_back(EET_phi1[4]);

          // Y Information
          Y_mass->push_back(EE.M());
          Y_px->push_back(EE.Px());
          Y_py->push_back(EE.Py());
          Y_pz->push_back(EE.Pz());
          Y_pt->push_back(EE.Pt());
          Y_eta->push_back(EE.Eta());
          Y_phi->push_back(EE.Phi());
          Y_rapidity->push_back(EE.Rapidity());

          Y_Vtx_Prob->push_back(tmpY_Prob);
          Y_Vtx_x->push_back(Y_Vtx.x());
          Y_Vtx_y->push_back(Y_Vtx.y());
          Y_Vtx_z->push_back(Y_Vtx.z());
          Y_Vtx_xError->push_back(Y_Vtx.xError());
          Y_Vtx_yError->push_back(Y_Vtx.yError());
          Y_Vtx_zError->push_back(Y_Vtx.zError());

          Y_Vtx_Px->push_back(Y_mom.Px());
          Y_Vtx_Py->push_back(Y_mom.Py());
          Y_Vtx_Pz->push_back(Y_mom.Pz());
          Y_Vtx_Pt->push_back(Y_mom.Pt());
          Y_Vtx_Eta->push_back(Y_mom.Eta());
          Y_Vtx_Phi->push_back(Y_mom.Phi());
          Y_Vtx_Rapidity->push_back(Y_mom.Rapidity());
          Y_Vtx_Mass->push_back(Y_mom.mass());

          // separate the two electrons by low and high pT
          if (iEle1->pt() > iEle2->pt()) {
            Y_lowPt->push_back(iEle2->pt());
            Y_highPt->push_back(iEle1->pt());
          } else {
            Y_lowPt->push_back(iEle1->pt());
            Y_highPt->push_back(iEle2->pt());
          }

          Y_px1->push_back(iEle1->px());
          Y_py1->push_back(iEle1->py());
          Y_pz1->push_back(iEle1->pz());
          Y_pt1->push_back(iEle1->pt());
          Y_eta1->push_back(iEle1->eta());
          Y_SCeta1->push_back(iEle1->superCluster()->eta());
          Y_phi1->push_back(iEle1->phi());
          Y_energy1->push_back(iEle1->energy());
          Y_energyCorr1->push_back(iEle1->userFloat("ecalTrkEnergyPostCorr"));
          Y_ecalIso1->push_back(iEle1->ecalIso());
          Y_hcalIso1->push_back(iEle1->hcalIso());
          Y_trackIso1->push_back(iEle1->trackIso());
          Z_trackIso1->push_back(iMuon1->trackIso());
          Y_charge1->push_back(iEle1->charge());

          Y_fit_pt1->push_back(fit_electron1_TT.track().pt());
          Y_fit_ptError1->push_back(fit_electron1_TT.track().ptError());
          Y_fit_eta1->push_back(fit_electron1_TT.track().eta());
          Y_fit_phi1->push_back(fit_electron1_TT.track().phi());

          // Electron IDs
          Y_CutBaseLoose1->push_back(iEle1->electronID("cutBasedElectronID-Fall17-94X-V2-loose"));
          Y_CutBaseVeto1->push_back(iEle1->electronID("cutBasedElectronID-Fall17-94X-V2-veto"));
          Y_mvaIsoWP90_1->push_back(iEle1->electronID("mvaEleID-Fall17-iso-V2-wp90"));
          Y_mvaIsoWP80_1->push_back(iEle1->electronID("mvaEleID-Fall17-iso-V2-wp80"));

          Y_px2->push_back(iEle2->px());
          Y_py2->push_back(iEle2->py());
          Y_pz2->push_back(iEle2->pz());
          Y_pt2->push_back(iEle2->pt());
          Y_eta2->push_back(iEle2->eta());
          Y_SCeta2->push_back(iEle2->superCluster()->eta());
          Y_phi2->push_back(iEle2->phi());
          Y_energy2->push_back(iEle2->energy());
          Y_energyCorr2->push_back(iEle2->userFloat("ecalTrkEnergyPostCorr"));
          Y_ecalIso2->push_back(iEle2->ecalIso());
          Y_hcalIso2->push_back(iEle2->hcalIso());
          Y_trackIso2->push_back(iEle2->trackIso());
          Z_trackIso2->push_back(iMuon2->trackIso());
          Y_charge2->push_back(iEle2->charge());

          Y_fit_pt2->push_back(fit_electron2_TT.track().pt());
          Y_fit_ptError2->push_back(fit_electron2_TT.track().ptError());
          Y_fit_eta2->push_back(fit_electron2_TT.track().eta());
          Y_fit_phi2->push_back(fit_electron2_TT.track().phi());

          // Electron IDs
          Y_CutBaseLoose2->push_back(iEle2->electronID("cutBasedElectronID-Fall17-94X-V2-loose"));
          Y_CutBaseVeto2->push_back(iEle2->electronID("cutBasedElectronID-Fall17-94X-V2-veto"));
          Y_mvaIsoWP90_2->push_back(iEle2->electronID("mvaEleID-Fall17-iso-V2-wp90"));
          Y_mvaIsoWP80_2->push_back(iEle2->electronID("mvaEleID-Fall17-iso-V2-wp80"));

          Y_dxy1->push_back(iEle1->gsfTrack()->dxy(bestVtx.position()));
          Y_dz1->push_back(iEle1->gsfTrack()->dz(bestVtx.position()));
          Y_dxy2->push_back(iEle2->gsfTrack()->dxy(bestVtx.position()));
          Y_dz2->push_back(iEle2->gsfTrack()->dz(bestVtx.position()));

          Mu_dCA->push_back(MuCA);
          Ele_dCA->push_back(EleCA);

          // Z Information
          Z_mass->push_back(MM.M());
          Z_px->push_back(MM.Px());
          Z_py->push_back(MM.Py());
          Z_pz->push_back(MM.Pz());
          Z_pt->push_back(MM.Pt());
          Z_eta->push_back(MM.Eta());
          Z_phi->push_back(MM.Phi());
          Z_rapidity->push_back(MM.Rapidity());

          Z_Vtx_Prob->push_back(tmpZ_Prob);
          Z_Vtx_x->push_back(Z_Vtx.x());
          Z_Vtx_y->push_back(Z_Vtx.y());
          Z_Vtx_z->push_back(Z_Vtx.z());
          Z_Vtx_xError->push_back(Z_Vtx.xError());
          Z_Vtx_yError->push_back(Z_Vtx.yError());
          Z_Vtx_zError->push_back(Z_Vtx.zError());

          Z_Vtx_Px->push_back(Z_mom.Px());
          Z_Vtx_Py->push_back(Z_mom.Py());
          Z_Vtx_Pz->push_back(Z_mom.Pz());
          Z_Vtx_Pt->push_back(Z_mom.Pt());
          Z_Vtx_Eta->push_back(Z_mom.Eta());
          Z_Vtx_Phi->push_back(Z_mom.Phi());
          Z_Vtx_Rapidity->push_back(Z_mom.Rapidity());
          Z_Vtx_Mass->push_back(Z_mom.mass());

          // separate the two muons by low and high pT
          if (iMuon1->track()->pt() > iMuon2->track()->pt()) {
            Z_lowPt->push_back(iMuon2->track()->pt());
            Z_highPt->push_back(iMuon1->track()->pt());
          } else {
            Z_lowPt->push_back(iMuon1->track()->pt());
            Z_highPt->push_back(iMuon2->track()->pt());
          }

          Z_px1->push_back(iMuon1->track()->px());
          Z_py1->push_back(iMuon1->track()->py());
          Z_pz1->push_back(iMuon1->track()->pz());
          Z_pt1->push_back(iMuon1->track()->pt());
          Z_eta1->push_back(iMuon1->track()->eta());
          Z_phi1->push_back(iMuon1->track()->phi());
          Z_charge1->push_back(iMuon1->charge());
          Z_soft1->push_back(iMuon1->isSoftMuon(bestVtx));
          Z_tight1->push_back(iMuon1->isTightMuon(bestVtx));
          Z_loose1->push_back(muon::isLooseMuon(*iMuon1));

          Z_fit_pt1->push_back(fit_muon1_TT.track().pt());
          Z_fit_ptError1->push_back(fit_muon1_TT.track().ptError());
          Z_fit_eta1->push_back(fit_muon1_TT.track().eta());
          Z_fit_phi1->push_back(fit_muon1_TT.track().phi());

          Z_px2->push_back(iMuon2->track()->px());
          Z_py2->push_back(iMuon2->track()->py());
          Z_pz2->push_back(iMuon2->track()->pz());
          Z_pt2->push_back(iMuon2->track()->pt());
          Z_eta2->push_back(iMuon2->track()->eta());
          Z_phi2->push_back(iMuon2->track()->phi());
          Z_charge2->push_back(iMuon2->charge());
          Z_soft2->push_back(iMuon2->isSoftMuon(bestVtx));
          Z_tight2->push_back(iMuon2->isTightMuon(bestVtx));
          Z_loose2->push_back(muon::isLooseMuon(*iMuon2));

          Z_fit_pt2->push_back(fit_muon2_TT.track().pt());
          Z_fit_ptError2->push_back(fit_muon2_TT.track().ptError());
          Z_fit_eta2->push_back(fit_muon2_TT.track().eta());
          Z_fit_phi2->push_back(fit_muon2_TT.track().phi());

          Z_xy1->push_back(glbTrack1->dxy(bestVtx.position()));
          Z_xy2->push_back(glbTrack2->dxy(bestVtx.position()));
          Z_z1->push_back(glbTrack1->dz(bestVtx.position()));
          Z_z2->push_back(glbTrack2->dz(bestVtx.position()));

          Mu1_C2->push_back(glbTrack1->normalizedChi2());
          //mumAngT->push_back( muon::isGoodMuon(*iMuon1,muon::TMLastStationAngTight) ); //
          Mu1_NHits->push_back(glbTrack1->numberOfValidHits());
          Mu1_NPHits->push_back(glbTrack1->hitPattern().numberOfValidPixelHits());
          Mu2_C2->push_back(glbTrack2->normalizedChi2());
          //mupAngT->push_back( muon::isGoodMuon(*iMuon2,muon::TMLastStationAngTight) );  //
          Mu2_NHits->push_back(glbTrack2->numberOfValidHits());
          Mu2_NPHits->push_back(glbTrack2->hitPattern().numberOfValidPixelHits());

          nCandi++;
        }
      }
    }
  }

  if ((nCandi + nG1 + nG2) > 0) {
    tree_->Fill();
  }

  nCandi = 0;
  //Gen Info Clear
  Y_GenMuonPt->clear();
  Y_GenMuonEta->clear();
  Y_GenMuonPhi->clear();
  Z_GenMuonPt->clear();
  Z_GenMuonEta->clear();
  Z_GenMuonPhi->clear();

  Run->clear();
  LumiBlock->clear();
  Event->clear();
  FourL_mass->clear();
  FourL_px->clear();
  FourL_py->clear();
  FourL_pz->clear();
  FourL_pt->clear();
  FourL_eta->clear();
  FourL_phi->clear();
  FourL_rapidity->clear();
  FourL_Vtx_Prob->clear();
  FourL_Vtx_x->clear();
  FourL_Vtx_y->clear();
  FourL_Vtx_z->clear();
  FourL_Vtx_xError->clear();
  FourL_Vtx_yError->clear();
  FourL_Vtx_zError->clear();

  FourL_Vtx_Px->clear();
  FourL_Vtx_Py->clear();
  FourL_Vtx_Pz->clear();
  FourL_Vtx_Pt->clear();
  FourL_Vtx_Eta->clear();
  FourL_Vtx_Phi->clear();
  FourL_Vtx_Rapidity->clear();
  FourL_Vtx_Mass->clear();

  Mu_TriggerPath->clear();
  Ele_TriggerPath->clear();
  Mu_TriggerPt1->clear();
  Mu_TriggerEta1->clear();
  Mu_TriggerPhi1->clear();
  Mu_TriggerPt2->clear();
  Mu_TriggerEta2->clear();
  Mu_TriggerPhi2->clear();
  Mu_TriggerPt3->clear();
  Mu_TriggerEta3->clear();
  Mu_TriggerPhi3->clear();
  Mu_TriggerPt4->clear();
  Mu_TriggerEta4->clear();
  Mu_TriggerPhi4->clear();
  Mu_TriggerPt5->clear();
  Mu_TriggerEta5->clear();
  Mu_TriggerPhi5->clear();

  Ele_TriggerPt1->clear();
  Ele_TriggerEta1->clear();
  Ele_TriggerPhi1->clear();
  Ele_TriggerPt2->clear();
  Ele_TriggerEta2->clear();
  Ele_TriggerPhi2->clear();
  Ele_TriggerPt3->clear();
  Ele_TriggerEta3->clear();
  Ele_TriggerPhi3->clear();
  Ele_TriggerPt4->clear();
  Ele_TriggerEta4->clear();
  Ele_TriggerPhi4->clear();
  Ele_TriggerPt5->clear();
  Ele_TriggerEta5->clear();
  Ele_TriggerPhi5->clear();

  Y_mass->clear();
  Y_Vtx_Prob->clear();
  Y_px->clear();
  Y_py->clear();
  Y_pz->clear();
  Y_pt->clear();
  Y_eta->clear();
  Y_phi->clear();
  Y_rapidity->clear();
  Y_Vtx_Px->clear();
  Y_Vtx_Py->clear();
  Y_Vtx_Pz->clear();
  Y_Vtx_Pt->clear();
  Y_Vtx_Eta->clear();
  Y_Vtx_Phi->clear();
  Y_Vtx_Mass->clear();
  Y_Vtx_x->clear();
  Y_Vtx_y->clear();
  Y_Vtx_z->clear();
  Y_Vtx_xError->clear();
  Y_Vtx_yError->clear();
  Y_Vtx_zError->clear();
  Y_Vtx_Rapidity->clear();
  Y_px1->clear();
  Y_py1->clear();
  Y_pz1->clear();
  Y_charge1->clear();
  Y_fit_pt1->clear();
  Y_fit_ptError1->clear();
  Y_fit_eta1->clear();
  Y_fit_phi1->clear();
  Y_pt1->clear();
  Y_eta1->clear();
  Y_SCeta1->clear();
  Y_phi1->clear();
  Y_energy1->clear();
  Y_energyCorr1->clear();
  Y_ecalIso1->clear();
  Y_hcalIso1->clear();
  Y_trackIso1->clear();
  Z_trackIso1->clear();
  Y_CutBaseLoose1->clear();
  Y_CutBaseVeto1->clear();
  Y_mvaIsoWP90_1->clear();
  Y_mvaIsoWP80_1->clear();
  Y_px2->clear();
  Y_py2->clear();
  Y_pz2->clear();
  Y_charge2->clear();
  Y_fit_pt2->clear();
  Y_fit_ptError2->clear();
  Y_fit_eta2->clear();
  Y_fit_phi2->clear();
  Y_pt2->clear();
  Y_eta2->clear();
  Y_SCeta2->clear();
  Y_phi2->clear();
  Y_energy2->clear();
  Y_energyCorr2->clear();
  Y_ecalIso2->clear();
  Y_hcalIso2->clear();
  Y_trackIso2->clear();
  Z_trackIso2->clear();
  Y_CutBaseLoose2->clear();
  Y_CutBaseVeto2->clear();
  Y_mvaIsoWP90_2->clear();
  Y_mvaIsoWP80_2->clear();
  Y_dxy1->clear();
  Y_dxy2->clear();
  Y_dz1->clear();
  Y_dz2->clear();

  Y_lowPt->clear();
  Y_highPt->clear();
  Z_lowPt->clear();
  Z_highPt->clear();
  Mu_dCA->clear();
  Ele_dCA->clear();
  Z_mass->clear();
  Z_px->clear();
  Z_py->clear();
  Z_pz->clear();
  Z_pt->clear();
  Z_eta->clear();
  Z_phi->clear();
  Z_rapidity->clear();
  Z_Vtx_Px->clear();
  Z_Vtx_Py->clear();
  Z_Vtx_Pz->clear();
  Z_Vtx_Pt->clear();
  Z_Vtx_Eta->clear();
  Z_Vtx_Phi->clear();
  Z_Vtx_Rapidity->clear();
  Z_Vtx_Mass->clear();
  Z_Vtx_x->clear();
  Z_Vtx_y->clear();
  Z_Vtx_z->clear();
  Z_Vtx_xError->clear();
  Z_Vtx_yError->clear();
  Z_Vtx_zError->clear();
  Z_px1->clear();
  Z_py1->clear();
  Z_pz1->clear();
  Z_charge1->clear();
  Z_pt1->clear();
  Z_eta1->clear();
  Z_phi1->clear();
  Z_soft1->clear();
  Z_tight1->clear();
  Z_loose1->clear();

  Z_fit_pt1->clear();
  Z_fit_ptError1->clear();
  Z_fit_eta1->clear();
  Z_fit_phi1->clear();

  Z_px2->clear();
  Z_py2->clear();
  Z_pz2->clear();
  Z_charge2->clear();
  Z_pt2->clear();
  Z_eta2->clear();
  Z_phi2->clear();
  Z_Vtx_Prob->clear();
  Z_soft2->clear();
  Z_tight2->clear();
  Z_loose2->clear();

  Z_fit_pt2->clear();
  Z_fit_ptError2->clear();
  Z_fit_eta2->clear();
  Z_fit_phi2->clear();

  Z_xy1->clear();
  Z_xy2->clear();
  Z_z2->clear();
  Z_z1->clear();

  Mu1_C2->clear();
  Mu1_NHits->clear();
  Mu1_NPHits->clear();
  Mu2_C2->clear();
  Mu2_NHits->clear();
  Mu2_NPHits->clear();
}

// ------------ method called once each job just before starting event loop  ------------

void miniAODmuons::beginJob() {
  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  std::string title = DataTypeString + " with " + MuonTriggerString + " " + ElectronTriggerString;

  // Add informetion to the title of the TTree
  tree_ = new TTree("ntuple", title.c_str());

  tree_->Branch("nCandi", &nCandi, "nCandi/i");
  //gen branches
  tree_->Branch("Y_GenMuonPt", &Y_GenMuonPt);
  tree_->Branch("Y_GenMuonEta", &Y_GenMuonEta);
  tree_->Branch("Y_GenMuonPhi", &Y_GenMuonPhi);
  tree_->Branch("Z_GenMuonPt", &Z_GenMuonPt);
  tree_->Branch("Z_GenMuonEta", &Z_GenMuonEta);
  tree_->Branch("Z_GenMuonPhi", &Z_GenMuonPhi);
  //electron channels
  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  tree_->Branch("Event", &Event);
  tree_->Branch("FourL_mass", &FourL_mass);
  tree_->Branch("FourL_px", &FourL_px);
  tree_->Branch("FourL_py", &FourL_py);
  tree_->Branch("FourL_pz", &FourL_pz);
  tree_->Branch("FourL_pt", &FourL_pt);
  tree_->Branch("FourL_eta", &FourL_eta);
  tree_->Branch("FourL_phi", &FourL_phi);
  tree_->Branch("FourL_rapidity", &FourL_rapidity);
  tree_->Branch("FourL_Vtx_Prob", &FourL_Vtx_Prob);
  tree_->Branch("FourL_Vtx_x", &FourL_Vtx_x);
  tree_->Branch("FourL_Vtx_y", &FourL_Vtx_y);
  tree_->Branch("FourL_Vtx_z", &FourL_Vtx_z);
  tree_->Branch("FourL_Vtx_xError", &FourL_Vtx_xError);
  tree_->Branch("FourL_Vtx_yError", &FourL_Vtx_yError);
  tree_->Branch("FourL_Vtx_zError", &FourL_Vtx_zError);

  tree_->Branch("FourL_Vtx_Px", &FourL_Vtx_Px);
  tree_->Branch("FourL_Vtx_Py", &FourL_Vtx_Py);
  tree_->Branch("FourL_Vtx_Pz", &FourL_Vtx_Pz);
  tree_->Branch("FourL_Vtx_Pt", &FourL_Vtx_Pt);
  tree_->Branch("FourL_Vtx_Eta", &FourL_Vtx_Eta);
  tree_->Branch("FourL_Vtx_Phi", &FourL_Vtx_Phi);
  tree_->Branch("FourL_Vtx_Rapidity", &FourL_Vtx_Rapidity);
  tree_->Branch("FourL_Vtx_Mass", &FourL_Vtx_Mass);

  tree_->Branch("Mu_TriggerPath", &Mu_TriggerPath);
  tree_->Branch("Ele_TriggerPath", &Ele_TriggerPath);
  tree_->Branch("Mu_TriggerPt1", &Mu_TriggerPt1);
  tree_->Branch("Mu_TriggerEta1", &Mu_TriggerEta1);
  tree_->Branch("Mu_TriggerPhi1", &Mu_TriggerPhi1);
  tree_->Branch("Mu_TriggerPt2", &Mu_TriggerPt2);
  tree_->Branch("Mu_TriggerEta2", &Mu_TriggerEta2);
  tree_->Branch("Mu_TriggerPhi2", &Mu_TriggerPhi2);
  tree_->Branch("Mu_TriggerPt3", &Mu_TriggerPt3);
  tree_->Branch("Mu_TriggerEta3", &Mu_TriggerEta3);
  tree_->Branch("Mu_TriggerPhi3", &Mu_TriggerPhi3);
  tree_->Branch("Mu_TriggerPt4", &Mu_TriggerPt4);
  tree_->Branch("Mu_TriggerEta4", &Mu_TriggerEta4);
  tree_->Branch("Mu_TriggerPhi4", &Mu_TriggerPhi4);
  tree_->Branch("Mu_TriggerPt5", &Mu_TriggerPt5);
  tree_->Branch("Mu_TriggerEta5", &Mu_TriggerEta5);
  tree_->Branch("Mu_TriggerPhi5", &Mu_TriggerPhi5);

  tree_->Branch("Ele_TriggerPt1", &Ele_TriggerPt1);
  tree_->Branch("Ele_TriggerEta1", &Ele_TriggerEta1);
  tree_->Branch("Ele_TriggerPhi1", &Ele_TriggerPhi1);
  tree_->Branch("Ele_TriggerPt2", &Ele_TriggerPt2);
  tree_->Branch("Ele_TriggerEta2", &Ele_TriggerEta2);
  tree_->Branch("Ele_TriggerPhi2", &Ele_TriggerPhi2);
  tree_->Branch("Ele_TriggerPt3", &Ele_TriggerPt3);
  tree_->Branch("Ele_TriggerEta3", &Ele_TriggerEta3);
  tree_->Branch("Ele_TriggerPhi3", &Ele_TriggerPhi3);
  tree_->Branch("Ele_TriggerPt4", &Ele_TriggerPt4);
  tree_->Branch("Ele_TriggerEta4", &Ele_TriggerEta4);
  tree_->Branch("Ele_TriggerPhi4", &Ele_TriggerPhi4);
  tree_->Branch("Ele_TriggerPt5", &Ele_TriggerPt5);
  tree_->Branch("Ele_TriggerEta5", &Ele_TriggerEta5);
  tree_->Branch("Ele_TriggerPhi5", &Ele_TriggerPhi5);

  tree_->Branch("Y_mass", &Y_mass);
  tree_->Branch("Y_Vtx_Prob", &Y_Vtx_Prob);

  tree_->Branch("Y_px", &Y_px);
  tree_->Branch("Y_py", &Y_py);
  tree_->Branch("Y_pz", &Y_pz);
  tree_->Branch("Y_pt", &Y_pt);
  tree_->Branch("Y_eta", &Y_eta);
  tree_->Branch("Y_phi", &Y_phi);
  tree_->Branch("Y_rapidity", &Y_rapidity);
  tree_->Branch("Y_Vtx_Px", &Y_Vtx_Px);
  tree_->Branch("Y_Vtx_Py", &Y_Vtx_Py);
  tree_->Branch("Y_Vtx_Pz", &Y_Vtx_Pz);
  tree_->Branch("Y_Vtx_Pt", &Y_Vtx_Pt);
  tree_->Branch("Y_Vtx_Eta", &Y_Vtx_Eta);
  tree_->Branch("Y_Vtx_Phi", &Y_Vtx_Phi);
  tree_->Branch("Y_Vtx_Rapidity", &Y_Vtx_Rapidity);
  tree_->Branch("Y_Vtx_Mass", &Y_Vtx_Mass);

  tree_->Branch("Y_Vtx_x", &Y_Vtx_x);
  tree_->Branch("Y_Vtx_y", &Y_Vtx_y);
  tree_->Branch("Y_Vtx_z", &Y_Vtx_z);
  tree_->Branch("Y_Vtx_xError", &Y_Vtx_xError);
  tree_->Branch("Y_Vtx_yError", &Y_Vtx_yError);
  tree_->Branch("Y_Vtx_zError", &Y_Vtx_zError);

  tree_->Branch("Y_px1", &Y_px1);
  tree_->Branch("Y_py1", &Y_py1);
  tree_->Branch("Y_pz1", &Y_pz1);
  tree_->Branch("Y_pt1", &Y_pt1);
  tree_->Branch("Y_eta1", &Y_eta1);
  tree_->Branch("Y_SCeta1", &Y_SCeta1);
  tree_->Branch("Y_phi1", &Y_phi1);
  tree_->Branch("Y_energy1", &Y_energy1);
  tree_->Branch("Y_energyCorr1", &Y_energyCorr1);
  tree_->Branch("Y_ecalIso1", &Y_ecalIso1);
  tree_->Branch("Y_hcalIso1", &Y_hcalIso1);
  tree_->Branch("Y_trackIso1", &Y_trackIso1);
  tree_->Branch("Z_trackIso1", &Z_trackIso1);
  tree_->Branch("Y_CutBaseLoose1", &Y_CutBaseLoose1);
  tree_->Branch("Y_CutBaseVeto1", &Y_CutBaseVeto1);
  tree_->Branch("Y_mvaIsoWP90_1", &Y_mvaIsoWP90_1);
  tree_->Branch("Y_mvaIsoWP80_1", &Y_mvaIsoWP80_1);
  tree_->Branch("Y_charge1", &Y_charge1);
  tree_->Branch("Y_fit_pt1", &Y_fit_pt1);
  tree_->Branch("Y_fit_ptError1", &Y_fit_ptError1);
  tree_->Branch("Y_fit_eta1", &Y_fit_eta1);
  tree_->Branch("Y_fit_phi1", &Y_fit_phi1);
  tree_->Branch("Y_px2", &Y_px2);
  tree_->Branch("Y_py2", &Y_py2);
  tree_->Branch("Y_pz2", &Y_pz2);
  tree_->Branch("Y_pt2", &Y_pt2);
  tree_->Branch("Y_eta2", &Y_eta2);
  tree_->Branch("Y_SCeta2", &Y_SCeta2);
  tree_->Branch("Y_phi2", &Y_phi2);
  tree_->Branch("Y_energy2", &Y_energy2);
  tree_->Branch("Y_energyCorr2", &Y_energyCorr2);
  tree_->Branch("Y_ecalIso2", &Y_ecalIso2);
  tree_->Branch("Y_hcalIso2", &Y_hcalIso2);
  tree_->Branch("Y_trackIso2", &Y_trackIso2);
  tree_->Branch("Z_trackIso2", &Z_trackIso2);
  tree_->Branch("Y_CutBaseLoose2", &Y_CutBaseLoose2);
  tree_->Branch("Y_CutBaseVeto2", &Y_CutBaseVeto2);
  tree_->Branch("Y_mvaIsoWP90_2", &Y_mvaIsoWP90_2);
  tree_->Branch("Y_mvaIsoWP80_2", &Y_mvaIsoWP80_2);
  tree_->Branch("Y_charge2", &Y_charge2);
  tree_->Branch("Y_fit_pt2", &Y_fit_pt2);
  tree_->Branch("Y_fit_ptError2", &Y_fit_ptError2);
  tree_->Branch("Y_fit_eta2", &Y_fit_eta2);
  tree_->Branch("Y_fit_phi2", &Y_fit_phi2);
  tree_->Branch("Y_dxy1", &Y_dxy1);
  tree_->Branch("Y_dxy2", &Y_dxy2);
  tree_->Branch("Y_dz1", &Y_dz1);
  tree_->Branch("Y_dz2", &Y_dz2);

  tree_->Branch("Y_lowPt", &Y_lowPt);
  tree_->Branch("Y_highPt", &Y_highPt);
  tree_->Branch("Z_lowPt", &Z_lowPt);
  tree_->Branch("Z_highPt", &Z_highPt);

  tree_->Branch("Mu_dCA", &Mu_dCA);
  tree_->Branch("Ele_dCA", &Ele_dCA);
  tree_->Branch("Z_mass", &Z_mass);
  tree_->Branch("Z_px", &Z_px);
  tree_->Branch("Z_py", &Z_py);
  tree_->Branch("Z_pz", &Z_pz);
  tree_->Branch("Z_pt", &Z_pt);
  tree_->Branch("Z_eta", &Z_eta);
  tree_->Branch("Z_phi", &Z_phi);
  tree_->Branch("Z_rapidity", &Z_rapidity);
  tree_->Branch("Z_Vtx_Px", &Z_Vtx_Px);
  tree_->Branch("Z_Vtx_Py", &Z_Vtx_Py);
  tree_->Branch("Z_Vtx_Pz", &Z_Vtx_Pz);
  tree_->Branch("Z_Vtx_Pt", &Z_Vtx_Pt);
  tree_->Branch("Z_Vtx_Eta", &Z_Vtx_Eta);
  tree_->Branch("Z_Vtx_Phi", &Z_Vtx_Phi);
  tree_->Branch("Z_Vtx_Rapidity", &Z_Vtx_Rapidity);
  tree_->Branch("Z_Vtx_Mass", &Z_Vtx_Mass);
  tree_->Branch("Z_Vtx_x", &Z_Vtx_x);
  tree_->Branch("Z_Vtx_y", &Z_Vtx_y);
  tree_->Branch("Z_Vtx_z", &Z_Vtx_z);
  tree_->Branch("Z_Vtx_xError", &Z_Vtx_xError);
  tree_->Branch("Z_Vtx_yError", &Z_Vtx_yError);
  tree_->Branch("Z_Vtx_zError", &Z_Vtx_zError);

  tree_->Branch("Z_px1", &Z_px1);
  tree_->Branch("Z_py1", &Z_py1);
  tree_->Branch("Z_pz1", &Z_pz1);
  tree_->Branch("Z_pt1", &Z_pt1);
  tree_->Branch("Z_eta1", &Z_eta1);
  tree_->Branch("Z_phi1", &Z_phi1);
  tree_->Branch("Z_charge1", &Z_charge1);
  tree_->Branch("Z_soft1", &Z_soft1);
  tree_->Branch("Z_tight1", &Z_tight1);
  tree_->Branch("Z_loose1", &Z_loose1);

  tree_->Branch("Z_fit_pt1", &Z_fit_pt1);
  tree_->Branch("Z_fit_ptError1", &Z_fit_ptError1);
  tree_->Branch("Z_fit_eta1", &Z_fit_eta1);
  tree_->Branch("Z_fit_phi1", &Z_fit_phi1);

  tree_->Branch("Z_px2", &Z_px2);
  tree_->Branch("Z_py2", &Z_py2);
  tree_->Branch("Z_pz2", &Z_pz2);
  tree_->Branch("Z_pt2", &Z_pt2);
  tree_->Branch("Z_eta2", &Z_eta2);
  tree_->Branch("Z_phi2", &Z_phi2);
  tree_->Branch("Z_charge2", &Z_charge2);
  tree_->Branch("Z_soft2", &Z_soft2);
  tree_->Branch("Z_tight2", &Z_tight2);
  tree_->Branch("Z_loose2", &Z_loose2);

  tree_->Branch("Z_fit_pt2", &Z_fit_pt2);
  tree_->Branch("Z_fit_ptError2", &Z_fit_ptError2);
  tree_->Branch("Z_fit_eta2", &Z_fit_eta2);
  tree_->Branch("Z_fit_phi2", &Z_fit_phi2);

  tree_->Branch("Z_Vtx_Prob", &Z_Vtx_Prob);
  tree_->Branch("Z_xy1", &Z_xy1);
  tree_->Branch("Z_xy2", &Z_xy2);
  tree_->Branch("Z_z2", &Z_z2);
  tree_->Branch("Z_z1", &Z_z1);

  tree_->Branch("Mu1_C2", &Mu1_C2);
  tree_->Branch("Mu1_NHits", &Mu1_NHits);
  tree_->Branch("Mu1_NPHits", &Mu1_NPHits);
  tree_->Branch("Mu2_C2", &Mu2_C2);
  tree_->Branch("Mu2_NHits", &Mu2_NHits);
  tree_->Branch("Mu2_NPHits", &Mu2_NPHits);
}

// ------------ method called once each job just after ending the event loop  ------------
void miniAODmuons::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(miniAODmuons);
