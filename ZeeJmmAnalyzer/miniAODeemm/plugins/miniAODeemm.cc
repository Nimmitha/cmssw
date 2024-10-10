// -*- C++ -*-
//
// Package:    miniAODeemm
// Class:      miniAODeemm
//

//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Monday Aug 28 (2017)         |
//         <jhovanny.andres.mejia.guisao@cern.ch> |
//=================================================

// system include files
#include <memory>

#include "ZeeJmmAnalyzer/miniAODeemm/plugins/miniAODeemm.h"

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

//the functions which actually match the trigger objects and see if it passes
namespace {
  std::vector<const pat::TriggerObjectStandAlone*> getMatchedObjs(const float eta,
                                                                  const float phi,
                                                                  const std::vector<pat::TriggerObjectStandAlone>& trigObjs,
                                                                  const float maxDeltaR = 0.1) {
    std::vector<const pat::TriggerObjectStandAlone*> matchedObjs;
    const float maxDR2 = maxDeltaR * maxDeltaR;
    for (auto& trigObj : trigObjs) {
      const float dR2 = reco::deltaR2(eta, phi, trigObj.eta(), trigObj.phi());
      if (dR2 < maxDR2)
        matchedObjs.push_back(&trigObj);
    }
    return matchedObjs;
  }
}  // namespace

//
// constants, enums and typedefs
//
Double_t Events = 0;

typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//

miniAODeemm::miniAODeemm(const edm::ParameterSet& iConfig)
    : dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
      dielectron_Label(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("dielectron"))),
      trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
      primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      //trigger
      triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
      //GenLevel Info
      prunedGenToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("pruned"))),
      ElectronTriggerString(iConfig.getParameter<std::string>("ElectronTrigger")),
      DataTypeString(iConfig.getParameter<std::string>("DataType")),
      isMC_(iConfig.getParameter<bool>("isMC")),

      tree_(0),
      //Gen add
      B_Z_GenMuonPt(0),
      B_Z_GenMuonEta(0),
      B_Z_GenMuonPhi(0),
      B_J_GenMuonPt(0),
      B_J_GenMuonEta(0),
      B_J_GenMuonPhi(0),

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

      FourL_VtxProb(0),
      FourL_PVx(0),
      FourL_PVy(0),
      FourL_PVz(0),
      FourL_PVxError(0),
      FourL_PVyError(0),
      FourL_PVzError(0),

      //Z
      Ele_TriggerPath(0),
      B_Z_TriggerPt1(0),
      B_Z_TriggerEta1(0),
      B_Z_TriggerPhi1(0),
      B_Z_TriggerPt2(0),
      B_Z_TriggerEta2(0),
      B_Z_TriggerPhi2(0),
      B_Z_TriggerPt3(0),
      B_Z_TriggerEta3(0),
      B_Z_TriggerPhi3(0),
      B_Z_TriggerPt4(0),
      B_Z_TriggerEta4(0),
      B_Z_TriggerPhi4(0),
      B_Z_TriggerPt5(0),
      B_Z_TriggerEta5(0),
      B_Z_TriggerPhi5(0),

      B_Z_Trigger32Pt1(0),
      B_Z_Trigger32Eta1(0),
      B_Z_Trigger32Phi1(0),
      B_Z_Trigger32Pt2(0),
      B_Z_Trigger32Eta2(0),
      B_Z_Trigger32Phi2(0),
      B_Z_Trigger32Pt3(0),
      B_Z_Trigger32Eta3(0),
      B_Z_Trigger32Phi3(0),
      B_Z_Trigger32Pt4(0),
      B_Z_Trigger32Eta4(0),
      B_Z_Trigger32Phi4(0),
      B_Z_Trigger32Pt5(0),
      B_Z_Trigger32Eta5(0),
      B_Z_Trigger32Phi5(0),

      B_Z_mass(0),
      B_Z_VtxProb(0),
      B_Z_px(0),
      B_Z_py(0),
      B_Z_pz(0),
      B_Z_pt(0),
      B_Z_eta(0),
      B_Z_phi(0),
      B_Z_rapidity(0),
      B_Z_VtxPx(0),
      B_Z_VtxPy(0),
      B_Z_VtxPz(0),
      B_Z_VtxPt(0),
      B_Z_VtxEta(0),
      B_Z_VtxPhi(0),
      B_Z_VtxRapidity(0),
      B_Z_VtxMass(0),

      B_Z_PVx(0),
      B_Z_PVy(0),
      B_Z_PVz(0),
      B_Z_PVxError(0),
      B_Z_PVyError(0),
      B_Z_PVzError(0),
      B_Z_px1(0),
      B_Z_py1(0),
      B_Z_pz1(0),
      B_Z_pt1(0),
      B_Z_eta1(0),
      B_Z_SCeta1(0),
      B_Z_phi1(0),
      B_Z_energy1(0),
      B_Z_energyCorr1(0),
      B_Z_ecalIso1(0),
      B_Z_hcalIso1(0),
      B_Z_trackIso1(0),
      B_Z_looseCutBase1(0),
      B_Z_mediumCutBase1(0),
      B_Z_tightCutBase1(0),
      B_Z_mvaIsoWP80_1(0),
      B_Z_mvaIsoWP90_1(0),
      B_Z_px2(0),
      B_Z_py2(0),
      B_Z_pz2(0),
      B_Z_pt2(0),
      B_Z_eta2(0),
      B_Z_SCeta2(0),
      B_Z_phi2(0),
      B_Z_energy2(0),
      B_Z_energyCorr2(0),
      B_Z_ecalIso2(0),
      B_Z_hcalIso2(0),
      B_Z_trackIso2(0),
      B_Z_looseCutBase2(0),
      B_Z_mediumCutBase2(0),
      B_Z_tightCutBase2(0),
      B_Z_mvaIsoWP80_2(0),
      B_Z_mvaIsoWP90_2(0),
      B_Z_charge1(0),
      B_Z_charge2(0),
      B_Z_dxy1(0),
      B_Z_dxy2(0),
      B_Z_dz1(0),
      B_Z_dz2(0),
      B_J_lowPt(0),
      B_J_highPt(0),
      B_Z_lowPt(0),
      B_Z_highPt(0),
      B_J_dca(0),
      B_J_mass(0),
      B_J_px(0),
      B_J_py(0),
      B_J_pz(0),
      B_J_pt(0),
      B_J_eta(0),
      B_J_phi(0),
      B_J_rapidity(0),
      B_J_VtxPx(0),
      B_J_VtxPy(0),
      B_J_VtxPz(0),
      B_J_VtxPt(0),
      B_J_VtxEta(0),
      B_J_VtxPhi(0),
      B_J_VtxRapidity(0),
      B_J_VtxMass(0),
      B_J_PVx(0),
      B_J_PVy(0),
      B_J_PVz(0),
      B_J_PVxError(0),
      B_J_PVyError(0),
      B_J_PVzError(0),
      B_J_px1(0),
      B_J_py1(0),
      B_J_pz1(0),
      B_J_pt1(0),
      B_J_eta1(0),
      B_J_phi1(0),
      B_J_soft1(0),
      B_J_tight1(0),
      B_J_loose1(0),
      B_J_px2(0),
      B_J_py2(0),
      B_J_pz2(0),
      B_J_pt2(0),
      B_J_eta2(0),
      B_J_phi2(0),
      B_J_charge1(0),
      B_J_charge2(0),
      B_J_soft2(0),
      B_J_tight2(0),
      B_J_loose2(0),
      B_J_VtxProb(0),
      B_J_xyP(0),
      B_J_xyM(0),
      B_J_zP(0),
      B_J_zM(0),

      mumC2(0),
      mumNHits(0),
      mumNPHits(0),
      mupC2(0),
      mupNHits(0),
      mupNPHits(0),

      nB(0)

{
  //now do what ever initialization is needed
}

miniAODeemm::~miniAODeemm() {}

//
// member functions
//

// ------------ method called to for each event  ------------
void miniAODeemm::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
    edm::LogWarning("miniAODeemm") << "no Transient Track in event";
    return;
  }

  if (!thePATElectronHandle.isValid()) {
    edm::LogWarning("miniAODeemm") << "no pat::Electrons in event";
    return;
  }
  if (!thePATMuonHandle.isValid()) {
    edm::LogWarning("miniAODeemm") << "no pat::Muons in event";
    return;
  }
  if (!triggerBits.isValid()) {
    edm::LogWarning("miniAODeemm") << "no Trigger path in event";
    return;
  }
  if (!triggerObjects.isValid()) {
    edm::LogWarning("miniAODeemm") << "no Trigger Object in event";
    return;
  }

  //Before Begining lets get Gen level Information
  int h = 0;

  float EET_pt[5] = {-999, -999, -999, -999, -999};
  float EET_eta[5] = {-999, -999, -999, -999, -999};
  float EET_phi[5] = {-999, -999, -999, -999, -999};
  float EET_pt1[5] = {-999, -999, -999, -999, -999};
  float EET_eta1[5] = {-999, -999, -999, -999, -999};
  float EET_phi1[5] = {-999, -999, -999, -999, -999};
  int nG1 = 0;
  int nG2 = 0;

  bool firedEleTrig = false;
  bool passedL1seed = false;

  bool GenInfo = false;
  if (GenInfo) {
    //Gen Level Info

    float B_J_GenMuon_pt = -999;
    float B_J_GenMuon_eta = -999;
    float B_J_GenMuon_phi = -999;
    float B_Z_GenMuon_pt = -999;
    float B_Z_GenMuon_eta = -999;
    float B_Z_GenMuon_phi = -999;

    for (size_t i = 0; i < pruned->size(); i++) {
      //if( (*pruned)[i].isPromptFinalState()  && abs((*pruned)[i].pdgId() ==13) ){

      if (abs((*pruned)[i].pdgId()) == 13) {
        //cout<<"Found gen level muon"<<endl;
        if ((*pruned)[i].mother()->pdgId() == 553) {
          //cout<<"Found gen level muon from Jpsi "<<endl;
          //if ( (*pruned)[i].->pdgId()==10443 ) {
          //cout<<"Found gen level muon with number of mother"<<(*pruned)[i].numberOfMothers()<<endl;
          B_J_GenMuon_pt = (*pruned)[i].pt();
          B_J_GenMuon_eta = (*pruned)[i].eta();
          B_J_GenMuon_phi = (*pruned)[i].phi();
          B_J_GenMuonPt->push_back(B_J_GenMuon_pt);
          B_J_GenMuonEta->push_back(B_J_GenMuon_eta);
          B_J_GenMuonPhi->push_back(B_J_GenMuon_phi);
          nG1++;
          //}
        }
      }
      if (abs((*pruned)[i].pdgId()) == 11) {
        //cout<<"Found gen level muon"<<endl;
        if ((*pruned)[i].mother()->pdgId() == 23) {
          //cout<<"Found gen level muon from Z"<<endl;
          B_Z_GenMuon_pt = (*pruned)[i].pt();
          B_Z_GenMuon_eta = (*pruned)[i].eta();
          B_Z_GenMuon_phi = (*pruned)[i].phi();
          B_Z_GenMuonPt->push_back(B_Z_GenMuon_pt);
          B_Z_GenMuonEta->push_back(B_Z_GenMuon_eta);
          B_Z_GenMuonPhi->push_back(B_Z_GenMuon_phi);
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
    // until we find a fired Electron triggers
    if (!firedEleTrig) {
      if (names.triggerName(i).find(ElectronTriggerString.c_str()) != string::npos) {
        if (triggerBits->accept(i))
          firedEleTrig = true;  // Electron trigger was fired
      }
    }
  }

  // Special treatment for 2017 when using HLT_Ele32_WPTight_Gsf_L1DoubleEG_v
  // It needs to pass hltEGL1SingleEGOrFilter filter as well
  // Together this will emulate HLT_Ele32_WPTight_Gsf which was not in the menu in 2017
  // https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
  if (firedEleTrig && ElectronTriggerString.find("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v") != string::npos) {
    // std::cout << "Special treatment for 2017 trigger" << std::endl;
    //so the filter names are all packed in miniAOD so we need to create a new collection of them which are unpacked
    std::vector<pat::TriggerObjectStandAlone> unpackedTrigObjs;
    for (auto& trigObj : *triggerObjects) {
      unpackedTrigObjs.push_back(trigObj);
      unpackedTrigObjs.back().unpackFilterLabels(iEvent, *triggerBits);
    }

    for (auto& ele : *thePATElectronHandle) {
      const float eta = ele.superCluster()->eta();
      const float phi = ele.superCluster()->phi();

      //now match ALL objects in a cone of DR<0.1
      //it is important to match all objects as there are different ways to reconstruct the same electron
      //eg, L1 seeded, unseeded, as a jet etc
      //and so you want to be sure you get all possible objects
      std::vector<const pat::TriggerObjectStandAlone*> matchedTrigObjs = getMatchedObjs(eta, phi, unpackedTrigObjs, 0.1);

      for (const auto& trigObj : unpackedTrigObjs) {
        // if (trigObj->hasFilterLabel("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter") && trigObj->hasFilterLabel("hltEGL1SingleEGOrFilter"))
        // I think we already checked the first filter by using the trigger name
        if (trigObj.hasFilterLabel("hltEGL1SingleEGOrFilter")) {
          passedL1seed = true;
          break;
        }
      }
      if (passedL1seed)
        break;
    }
    firedEleTrig = firedEleTrig && passedL1seed;
  }

  //****************************************************
  //***************Trigger Object***********************
  //****************************************************

  // if any of the triggers fired, check the trigger objects
  if (firedEleTrig) {
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {  // note: not "const &" since we want to call unpackPathNames
      obj.unpackPathNames(names);
      if (obj.hasPathName(ElectronTriggerString.c_str(), true, true) > 0) {
        // This needs to be fixed for 2017
        EET_pt1[h] = obj.pt();
        EET_eta1[h] = obj.eta();
        EET_phi1[h] = obj.phi();
        h++;
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

      if (abs(iMuon1->charge()) != 1 || abs(iMuon2->charge()) != 1)
        continue;

      for (View<pat::Electron>::const_iterator iEle1 = thePATElectronHandle->begin(); iEle1 != thePATElectronHandle->end(); ++iEle1) {
        for (View<pat::Electron>::const_iterator iEle2 = iEle1 + 1; iEle2 != thePATElectronHandle->end(); ++iEle2) {
          // std::cout << "Begining of cut on Ele loop" << std::endl;
          Events = iEvent.id().event();

          //  cout<<"Begining of cut on Ele 2"<<endl;
          if ((iEle1->charge()) * (iEle2->charge()) == 1)
            continue;

          if (abs(iEle1->charge()) != 1 || abs(iEle2->charge()) != 1)
            continue;

          TrackRef glbTrackP;
          TrackRef glbTrackM;
          // cout<<"Electron Track Pt"<<iEle1->gsfTrack()->pt()<<endl;

          if (iMuon1->charge() == 1) {
            glbTrackP = iMuon1->track();
            glbTrackM = iMuon2->track();
          } else if (iMuon1->charge() == -1) {
            glbTrackM = iMuon1->track();
            glbTrackP = iMuon2->track();
          }

          if (glbTrackP.isNull() || glbTrackM.isNull()) {
            //  std::cout << "continue due to no track ref" << endl;
            continue;
          }

          reco::TrackRef kfTrackRefP;
          reco::TrackRef kfTrackRefM;

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

          if (MM.M() < 0 || MM.M() > 12)
            continue;
          if (EE.M() < 70 || EE.M() > 110)
            continue;

          // cout << "Start looking for track quality" << endl;
          int tkquality = 0;

          if (iEle1->gsfTrack().isAvailable() && iEle1->gsfTrack().isNonnull()) {
            if (iEle2->gsfTrack().isAvailable() && iEle2->gsfTrack().isNonnull()) {
              tkquality++;
            }
          }
          if (tkquality == 0)
            continue;

          if (iMuon1->track()->pt() < 2.0)
            continue;
          if (iMuon2->track()->pt() < 2.0)
            continue;

          //cout<<"Start looking muon track quality"<<endl;
          if (!(glbTrackM->quality(reco::TrackBase::highPurity)))
            continue;
          if (!(glbTrackP->quality(reco::TrackBase::highPurity)))
            continue;
          // cout<<"Start Building Track"<<endl;

          reco::TransientTrack muon1TT((*theB).build(glbTrackP));
          reco::TransientTrack muon2TT((*theB).build(glbTrackM));
          reco::TransientTrack electron1TT((*theB).build(iEle1->gsfTrack()));
          reco::TransientTrack electron2TT((*theB).build(iEle2->gsfTrack()));

          // *****  Trajectory states to calculate DCA for the 2 muons *********************
          FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
          FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();
          // FreeTrajectoryState electron1State = muon1TT.impactPointTSCP().theState();
          // FreeTrajectoryState electron2State = muon2TT.impactPointTSCP().theState();
          // cout<<"Start validating impact point"<<endl;
          if (!muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid())
            continue;
          // cout<<"Validated  impact point for muon"<<endl;
          if (!electron1TT.impactPointTSCP().isValid() || !electron2TT.impactPointTSCP().isValid())
            continue;
          // cout<<"Start validating impact point for electron"<<endl;
          // Measure distance between tracks at their closest approach
          ClosestApproachInRPhi cApp;
          cApp.calculate(mu1State, mu2State);
          if (!cApp.status()) {
            // std::cout << "cApp failed" << std::endl;
            continue;
          }
          float dca = fabs(cApp.distance());
          // if (dca < 0. || dca > 0.5) continue;
          // cout<<"dca"<<dca<<endl;
          // cout<<" closest approach  "<<dca<<endl;
          // ClosestApproachInRPhi cAppE;
          // cApp.calculate(electron1State, electron2State);
          // if( !cAppE.status() ) continue;
          // float dcaE = fabs( cAppE.distance() );
          // if (dcaE < 0. || dcaE > 0.5) continue;
          // cout<<"dca E "<<dcaE<<endl;

          //************************************************************
          //******************Kalman Vertex Fitting*********************
          //************************************************************

          vector<TransientTrack> ele_tks;
          vector<TransientTrack> mu_tks;
          vector<TransientTrack> mmee_tks;

          // Dielectron Vertex Fit
          KalmanVertexFitter kvf(true);
          ele_tks.clear();
          ele_tks.push_back(electron1TT);
          ele_tks.push_back(electron2TT);

          TransientVertex Z_candi = kvf.vertex(ele_tks);
          reco::Vertex Z_Vtx = Z_candi;
          const math::XYZTLorentzVectorD Z_mom = Z_Vtx.p4(ele_mass, 0.0);
          if (!Z_candi.isValid()) {
            // cout<<"Z candidate non validated by kalman fitter"<<endl;
            continue;
          }
          float B_Prob_tmp1 = TMath::Prob(Z_candi.totalChiSquared(), Z_candi.degreesOfFreedom());
          // cout << "vertex" << endl;
          //if (B_Prob_tmp1 < 0.001) continue;

          KalmanVertexFitter kvfM(true);
          mu_tks.clear();
          mu_tks.push_back(muon1TT);
          mu_tks.push_back(muon2TT);
          TransientVertex J_candi = kvfM.vertex(mu_tks);
          if (!J_candi.isValid()) {
            // cout<<"J candidate non validated by kalman fitter"<<endl;
            continue;
          }
          reco::Vertex JPsi_Vtx = J_candi;
          float B_Prob_tmp = TMath::Prob(J_candi.totalChiSquared(), J_candi.degreesOfFreedom());
          const math::XYZTLorentzVectorD JPsi_mom = JPsi_Vtx.p4(mu_mass, 0.0);

          // const math::XYZTLorentzVectorD jpsi_mom = jpsi_vtx.p4(0.1056583,0.0);

          // if (B_Prob_tmp < 0.001) continue;
          KalmanVertexFitter kvfEM(true);
          mmee_tks.clear();
          mmee_tks.push_back(electron1TT);
          mmee_tks.push_back(electron2TT);
          mmee_tks.push_back(muon1TT);
          mmee_tks.push_back(muon2TT);
          TransientVertex FourL_candi = kvfEM.vertex(mmee_tks);
          if (!FourL_candi.isValid()) {
            // cout<<"H candidate non validated by kalman fitter"<<endl;
            continue;
          }
          float B_Prob_tmp4L = TMath::Prob(FourL_candi.totalChiSquared(), FourL_candi.degreesOfFreedom());
          reco::Vertex FourL_Vtx = FourL_candi;
          // cout<<"vertex 4l"<<endl;
          if (B_Prob_tmp4L < 0.001) {
            // cout<<"B_Prob_tmp4L "<<B_Prob_tmp4L<<endl;
            continue;
          }
          std::cout << "FINAL VERTEX 4L" << std::endl;

          // ******   Let's check the vertex and mass ****

          //The mass of a muon and the insignificant mass sigma
          //to avoid singularities in the covariance matrix.

          bool KinFit = false;

          //cout<<"Start Kin Loop"<<endl;

          if (KinFit == true) {
            ParticleMass muon_mass = 0.10565837;  //pdg mass
            //ParticleMass psi_mass = 3.096916;
            float muon_sigma = muon_mass * 1.e-6;
            //float psi_sigma = psi_mass*1.e-6;
            vector<RefCountedKinematicParticle> muonParticles;
            //Creating a KinematicParticleFactory
            KinematicParticleFactoryFromTransientTrack pFactory;

            //initial chi2 and ndf before kinematic fits.
            float chi = 0.;
            float ndf = 0.;
            //vector<RefCountedKinematicParticle> muonParticles;
            try {
              muonParticles.push_back(pFactory.particle(muon1TT, muon_mass, chi, ndf, muon_sigma));
              muonParticles.push_back(pFactory.particle(muon2TT, muon_mass, chi, ndf, muon_sigma));
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

            //some loose cuts go here

            // if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
            // if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>12) continue;

            //fill variables?iMuon1->track()->pt()

            //B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
            //B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
            //B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
            //B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );
          }

          //Event Information
          Run->push_back(iEvent.id().run());
          LumiBlock->push_back(iEvent.luminosityBlock());
          Event->push_back(iEvent.id().event());

          //FourL Information
          // std::cout << "FourL Information" << std::endl;
          FourL_mass->push_back(EEMM.M());
          FourL_px->push_back(EEMM.Px());
          FourL_py->push_back(EEMM.Py());
          FourL_pz->push_back(EEMM.Pz());
          FourL_pt->push_back(EEMM.Pt());
          FourL_eta->push_back(EEMM.Eta());
          FourL_phi->push_back(EEMM.Phi());

          FourL_VtxProb->push_back(B_Prob_tmp4L);

          // Trigger Path Information
          // std::cout << "Trigger Path Information" << std::endl;
          Ele_TriggerPath->push_back(firedEleTrig);

          // Trigger Object Information
          // std::cout << "Trigger Object Information" << std::endl;
          B_Z_TriggerPt1->push_back(EET_pt[0]);
          B_Z_TriggerEta1->push_back(EET_eta[0]);
          B_Z_TriggerPhi1->push_back(EET_phi[0]);
          B_Z_TriggerPt2->push_back(EET_pt[1]);
          B_Z_TriggerEta2->push_back(EET_eta[1]);
          B_Z_TriggerPhi2->push_back(EET_phi[1]);
          B_Z_TriggerPt3->push_back(EET_pt[2]);
          B_Z_TriggerEta3->push_back(EET_eta[2]);
          B_Z_TriggerPhi3->push_back(EET_phi[2]);
          B_Z_TriggerPt4->push_back(EET_pt[3]);
          B_Z_TriggerEta4->push_back(EET_eta[3]);
          B_Z_TriggerPhi4->push_back(EET_phi[3]);
          B_Z_TriggerPt5->push_back(EET_pt[4]);
          B_Z_TriggerEta5->push_back(EET_eta[4]);
          B_Z_TriggerPhi5->push_back(EET_phi[4]);

          B_Z_Trigger32Pt1->push_back(EET_pt1[0]);
          B_Z_Trigger32Eta1->push_back(EET_eta1[0]);
          B_Z_Trigger32Phi1->push_back(EET_phi1[0]);
          B_Z_Trigger32Pt2->push_back(EET_pt1[1]);
          B_Z_Trigger32Eta2->push_back(EET_eta1[1]);
          B_Z_Trigger32Phi2->push_back(EET_phi1[1]);
          B_Z_Trigger32Pt3->push_back(EET_pt1[2]);
          B_Z_Trigger32Eta3->push_back(EET_eta1[2]);
          B_Z_Trigger32Phi3->push_back(EET_phi1[2]);
          B_Z_Trigger32Pt4->push_back(EET_pt1[3]);
          B_Z_Trigger32Eta4->push_back(EET_eta1[3]);
          B_Z_Trigger32Phi4->push_back(EET_phi1[3]);
          B_Z_Trigger32Pt5->push_back(EET_pt1[4]);
          B_Z_Trigger32Eta5->push_back(EET_eta1[4]);
          B_Z_Trigger32Phi5->push_back(EET_phi1[4]);

          B_Z_mass->push_back(EE.M());
          B_Z_VtxProb->push_back(B_Prob_tmp1);

          B_Z_px->push_back(EE.Px());
          B_Z_py->push_back(EE.Py());
          B_Z_pz->push_back(EE.Pz());
          B_Z_pt->push_back(EE.Pt());
          B_Z_eta->push_back(EE.Eta());
          B_Z_phi->push_back(EE.Phi());
          B_Z_rapidity->push_back(EE.Rapidity());
          B_Z_VtxPx->push_back(Z_mom.Px());
          B_Z_VtxPy->push_back(Z_mom.Py());
          B_Z_VtxPz->push_back(Z_mom.Pz());
          B_Z_VtxPt->push_back(Z_mom.Pt());
          B_Z_VtxEta->push_back(Z_mom.Eta());
          B_Z_VtxPhi->push_back(Z_mom.Phi());
          B_Z_VtxRapidity->push_back(Z_mom.Rapidity());
          B_Z_VtxMass->push_back(Z_mom.mass());

          B_Z_PVx->push_back(Z_Vtx.x());
          B_Z_PVy->push_back(Z_Vtx.y());
          B_Z_PVz->push_back(Z_Vtx.z());
          B_Z_PVxError->push_back(Z_Vtx.xError());
          B_Z_PVyError->push_back(Z_Vtx.yError());
          B_Z_PVzError->push_back(Z_Vtx.zError());

          // separate the two electrons by low and high pT
          // std::cout << "Separate the two electrons by low and high pT" << std::endl;
          if (iEle1->pt() > iEle2->pt()) {
            B_Z_lowPt->push_back(iEle2->pt());
            B_Z_highPt->push_back(iEle1->pt());
          } else {
            B_Z_lowPt->push_back(iEle1->pt());
            B_Z_highPt->push_back(iEle2->pt());
          }

          // std::cout << "electron1 information" << std::endl;
          B_Z_px1->push_back(iEle1->px());
          B_Z_py1->push_back(iEle1->py());
          B_Z_pz1->push_back(iEle1->pz());
          B_Z_pt1->push_back(iEle1->pt());
          B_Z_eta1->push_back(iEle1->eta());
          B_Z_SCeta1->push_back(iEle1->superCluster()->eta());
          B_Z_phi1->push_back(iEle1->phi());
          B_Z_energy1->push_back(iEle1->energy());
          // B_Z_energyCorr1->push_back(iEle1->userFloat("ecalTrkEnergyPostCorr"));

          B_Z_ecalIso1->push_back(iEle1->ecalIso());
          B_Z_hcalIso1->push_back(iEle1->hcalIso());
          B_Z_trackIso1->push_back(iEle1->trackIso());

          // std::cout << "electron1 id" << std::endl;

          B_Z_looseCutBase1->push_back(iEle1->electronID("cutBasedElectronID-Fall17-94X-V2-loose"));
          B_Z_mediumCutBase1->push_back(iEle1->electronID("cutBasedElectronID-Fall17-94X-V2-medium"));
          B_Z_tightCutBase1->push_back(iEle1->electronID("cutBasedElectronID-Fall17-94X-V2-tight"));
          // B_Z_mvaIsoWP80_1->push_back(iEle1->userFloat("mvaEleID-Fall17-iso-V2-wp80"));
          B_Z_mvaIsoWP90_1->push_back(iEle1->electronID("mvaEleID-Fall17-iso-V2-wp90"));

          B_Z_charge1->push_back(iEle1->charge());

          B_Z_px2->push_back(iEle2->px());
          B_Z_py2->push_back(iEle2->py());
          B_Z_pz2->push_back(iEle2->pz());
          B_Z_pt2->push_back(iEle2->pt());
          B_Z_eta2->push_back(iEle2->eta());
          B_Z_SCeta2->push_back(iEle2->superCluster()->eta());
          B_Z_phi2->push_back(iEle2->phi());
          B_Z_energy2->push_back(iEle2->energy());
          B_Z_energyCorr2->push_back(iEle2->userFloat("ecalTrkEnergyPostCorr"));
          B_Z_ecalIso2->push_back(iEle2->ecalIso());
          B_Z_hcalIso2->push_back(iEle2->hcalIso());
          B_Z_trackIso2->push_back(iEle2->trackIso());

          // std::cout << "electron2 id" << std::endl;
          B_Z_looseCutBase2->push_back(iEle2->electronID("cutBasedElectronID-Fall17-94X-V2-loose"));
          B_Z_mediumCutBase2->push_back(iEle2->electronID("cutBasedElectronID-Fall17-94X-V2-medium"));
          B_Z_tightCutBase2->push_back(iEle2->electronID("cutBasedElectronID-Fall17-94X-V2-tight"));
          // B_Z_mvaIsoWP80_2->push_back(iEle2->userFloat("mvaEleID-Fall17-iso-V2-wp80"));
          B_Z_mvaIsoWP90_2->push_back(iEle2->electronID("mvaEleID-Fall17-iso-V2-wp90"));

          B_Z_charge2->push_back(iEle2->charge());
          B_Z_dxy1->push_back(iEle1->gsfTrack()->dxy(bestVtx.position()));
          B_Z_dz1->push_back(iEle1->gsfTrack()->dz(bestVtx.position()));
          B_Z_dxy2->push_back(iEle2->gsfTrack()->dxy(bestVtx.position()));
          B_Z_dz2->push_back(iEle2->gsfTrack()->dz(bestVtx.position()));

          B_J_dca->push_back(dca);
          B_J_mass->push_back(MM.M());
          B_J_px->push_back(MM.Px());
          B_J_py->push_back(MM.Py());
          B_J_pz->push_back(MM.Pz());

          B_J_pt->push_back(MM.Pt());
          B_J_eta->push_back(MM.Eta());
          B_J_phi->push_back(MM.Phi());
          B_J_rapidity->push_back(MM.Rapidity());

          B_J_VtxPx->push_back(JPsi_mom.Px());
          B_J_VtxPy->push_back(JPsi_mom.Py());
          B_J_VtxPz->push_back(JPsi_mom.Pz());
          B_J_VtxPt->push_back(JPsi_mom.Pt());
          B_J_VtxEta->push_back(JPsi_mom.Eta());
          B_J_VtxPhi->push_back(JPsi_mom.Phi());
          B_J_VtxRapidity->push_back(JPsi_mom.Rapidity());
          B_J_VtxMass->push_back(JPsi_mom.mass());
          B_J_PVx->push_back(JPsi_Vtx.x());
          B_J_PVy->push_back(JPsi_Vtx.y());
          B_J_PVz->push_back(JPsi_Vtx.z());
          B_J_PVxError->push_back(JPsi_Vtx.xError());
          B_J_PVyError->push_back(JPsi_Vtx.yError());
          B_J_PVzError->push_back(JPsi_Vtx.zError());

          // separate the two muons by low and high pT
          // std::cout << "separate the two muons by low and high pT" << std::endl;
          if (iMuon1->track()->pt() > iMuon2->track()->pt()) {
            B_J_lowPt->push_back(iMuon2->track()->pt());
            B_J_highPt->push_back(iMuon1->track()->pt());
          } else {
            B_J_lowPt->push_back(iMuon1->track()->pt());
            B_J_highPt->push_back(iMuon2->track()->pt());
          }

          B_J_px1->push_back(iMuon1->track()->px());
          B_J_py1->push_back(iMuon1->track()->py());
          B_J_pz1->push_back(iMuon1->track()->pz());
          B_J_pt1->push_back(iMuon1->track()->pt());
          B_J_eta1->push_back(iMuon1->track()->eta());
          B_J_phi1->push_back(iMuon1->track()->phi());
          B_J_charge1->push_back(iMuon1->charge());
          B_J_soft1->push_back(iMuon1->isSoftMuon(bestVtx));
          B_J_tight1->push_back(iMuon1->isTightMuon(bestVtx));
          B_J_loose1->push_back(muon::isLooseMuon(*iMuon1));

          B_J_px2->push_back(iMuon2->track()->px());
          B_J_py2->push_back(iMuon2->track()->py());
          B_J_pz2->push_back(iMuon2->track()->pz());
          B_J_pt2->push_back(iMuon2->track()->pt());
          B_J_eta2->push_back(iMuon2->track()->eta());
          B_J_phi2->push_back(iMuon2->track()->phi());
          B_J_charge2->push_back(iMuon2->charge());
          B_J_VtxProb->push_back(B_Prob_tmp);
          B_J_soft2->push_back(iMuon2->isSoftMuon(bestVtx));
          B_J_tight2->push_back(iMuon2->isTightMuon(bestVtx));
          B_J_loose2->push_back(muon::isLooseMuon(*iMuon2));
          B_J_xyP->push_back(glbTrackP->dxy(bestVtx.position()));
          B_J_xyM->push_back(glbTrackM->dxy(bestVtx.position()));
          B_J_zP->push_back(glbTrackM->dz(bestVtx.position()));
          B_J_zM->push_back(glbTrackP->dz(bestVtx.position()));

          // std::cout << "mumC" << std::endl;
          mumC2->push_back(glbTrackP->normalizedChi2());
          //mumAngT->push_back( muon::isGoodMuon(*iMuon1,muon::TMLastStationAngTight) ); //
          mumNHits->push_back(glbTrackP->numberOfValidHits());
          mumNPHits->push_back(glbTrackP->hitPattern().numberOfValidPixelHits());
          mupC2->push_back(glbTrackM->normalizedChi2());
          //mupAngT->push_back( muon::isGoodMuon(*iMuon2,muon::TMLastStationAngTight) );  //
          mupNHits->push_back(glbTrackM->numberOfValidHits());
          mupNPHits->push_back(glbTrackM->hitPattern().numberOfValidPixelHits());

          nB++;
          // std::cout << "End of all loop" << std::endl;
          if (KinFit == true) {
            //muonParticles.clear();
          }
        }
      }
    }
  }

  if ((nB + nG1 + nG2) > 0) {
    std::cout << "filling tree" << endl;
    tree_->Fill();
  }

  nB = 0;
  //Gen Info Clear
  B_Z_GenMuonPt->clear();
  B_Z_GenMuonEta->clear();
  B_Z_GenMuonPhi->clear();
  B_J_GenMuonPt->clear();
  B_J_GenMuonEta->clear();
  B_J_GenMuonPhi->clear();

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
  FourL_VtxProb->clear();
  FourL_PVx->clear();
  FourL_PVy->clear();
  FourL_PVz->clear();
  FourL_PVxError->clear();
  FourL_PVyError->clear();
  FourL_PVzError->clear();

  //B_Z_dca->clear();
  Ele_TriggerPath->clear();
  B_Z_TriggerPt1->clear();
  B_Z_TriggerEta1->clear();
  B_Z_TriggerPhi1->clear();
  B_Z_TriggerPt2->clear();
  B_Z_TriggerEta2->clear();
  B_Z_TriggerPhi2->clear();
  B_Z_TriggerPt3->clear();
  B_Z_TriggerEta3->clear();
  B_Z_TriggerPhi3->clear();
  B_Z_TriggerPt4->clear();
  B_Z_TriggerEta4->clear();
  B_Z_TriggerPhi4->clear();
  B_Z_TriggerPt5->clear();
  B_Z_TriggerEta5->clear();
  B_Z_TriggerPhi5->clear();

  B_Z_Trigger32Pt1->clear();
  B_Z_Trigger32Eta1->clear();
  B_Z_Trigger32Phi1->clear();
  B_Z_Trigger32Pt2->clear();
  B_Z_Trigger32Eta2->clear();
  B_Z_Trigger32Phi2->clear();
  B_Z_Trigger32Pt3->clear();
  B_Z_Trigger32Eta3->clear();
  B_Z_Trigger32Phi3->clear();
  B_Z_Trigger32Pt4->clear();
  B_Z_Trigger32Eta4->clear();
  B_Z_Trigger32Phi4->clear();
  B_Z_Trigger32Pt5->clear();
  B_Z_Trigger32Eta5->clear();
  B_Z_Trigger32Phi5->clear();

  B_Z_mass->clear();
  B_Z_VtxProb->clear();
  B_Z_px->clear();
  B_Z_py->clear();
  B_Z_pz->clear();
  B_Z_pt->clear();
  B_Z_eta->clear();
  B_Z_phi->clear();
  B_Z_rapidity->clear();
  B_Z_VtxPx->clear();
  B_Z_VtxPy->clear();
  B_Z_VtxPz->clear();
  B_Z_VtxPt->clear();
  B_Z_VtxEta->clear();
  B_Z_VtxPhi->clear();
  B_Z_VtxMass->clear();
  B_Z_PVx->clear();
  B_Z_PVy->clear();
  B_Z_PVz->clear();
  B_Z_PVxError->clear();
  B_Z_PVyError->clear();
  B_Z_PVzError->clear();
  B_Z_VtxRapidity->clear();
  B_Z_px1->clear();
  B_Z_py1->clear();
  B_Z_pz1->clear();
  B_Z_charge1->clear();
  B_Z_pt1->clear();
  B_Z_eta1->clear();
  B_Z_SCeta1->clear();
  B_Z_phi1->clear();
  B_Z_energy1->clear();
  B_Z_energyCorr1->clear();
  B_Z_ecalIso1->clear();
  B_Z_hcalIso1->clear();
  B_Z_trackIso1->clear();
  B_Z_looseCutBase1->clear();
  B_Z_mediumCutBase1->clear();
  B_Z_tightCutBase1->clear();
  B_Z_mvaIsoWP80_1->clear();
  B_Z_mvaIsoWP90_1->clear();
  B_Z_px2->clear();
  B_Z_py2->clear();
  B_Z_pz2->clear();
  B_Z_charge2->clear();
  B_Z_pt2->clear();
  B_Z_eta2->clear();
  B_Z_SCeta2->clear();
  B_Z_phi2->clear();
  B_Z_energy2->clear();
  B_Z_energyCorr2->clear();
  B_Z_ecalIso2->clear();
  B_Z_hcalIso2->clear();
  B_Z_trackIso2->clear();
  B_Z_looseCutBase2->clear();
  B_Z_mediumCutBase2->clear();
  B_Z_tightCutBase2->clear();
  B_Z_mvaIsoWP80_2->clear();
  B_Z_mvaIsoWP90_2->clear();
  B_Z_dxy1->clear();
  B_Z_dxy2->clear();
  B_Z_dz1->clear();
  B_Z_dz2->clear();

  B_J_lowPt->clear();
  B_J_highPt->clear();
  B_Z_lowPt->clear();
  B_Z_highPt->clear();
  B_J_dca->clear();
  B_J_mass->clear();
  B_J_px->clear();
  B_J_py->clear();
  B_J_pz->clear();
  B_J_pt->clear();
  B_J_eta->clear();
  B_J_phi->clear();
  B_J_rapidity->clear();
  B_J_VtxPx->clear();
  B_J_VtxPy->clear();
  B_J_VtxPz->clear();
  B_J_VtxPt->clear();
  B_J_VtxEta->clear();
  B_J_VtxPhi->clear();
  B_J_VtxMass->clear();
  B_J_PVx->clear();
  B_J_PVy->clear();
  B_J_PVz->clear();
  B_J_PVxError->clear();
  B_J_PVyError->clear();
  B_J_PVzError->clear();
  B_J_px1->clear();
  B_J_py1->clear();
  B_J_pz1->clear();
  B_J_charge1->clear();
  B_J_pt1->clear();
  B_J_eta1->clear();
  B_J_phi1->clear();
  B_J_soft1->clear();
  B_J_tight1->clear();
  B_J_loose1->clear();
  B_J_px2->clear();
  B_J_py2->clear();
  B_J_pz2->clear();
  B_J_charge2->clear();
  B_J_pt2->clear();
  B_J_eta2->clear();
  B_J_phi2->clear();
  B_J_VtxProb->clear();
  B_J_soft2->clear();
  B_J_tight2->clear();
  B_J_loose2->clear();
  B_J_xyP->clear();
  B_J_xyM->clear();
  B_J_zP->clear();
  B_J_zM->clear();

  mumC2->clear();
  mumNHits->clear();
  mumNPHits->clear();
  mupC2->clear();
  mupNHits->clear();
  mupNPHits->clear();
}

// ------------ method called once each job just before starting event loop  ------------

void miniAODeemm::beginJob() {
  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  std::string title = DataTypeString + " with " + ElectronTriggerString;

  //edm::Service<TFileService> fs;
  //tree_ = fs->make<TTree>("ntuple"," J/psi ntuple");

  //tree_->Branch("nB",&nB,"nB/i");
  // tree_ = new TTree("ntuple", "ntuple");
  // Add informetion to the title of the TTree
  tree_ = new TTree("ntuple", title.c_str());

  tree_->Branch("nB", &nB, "nB/i");
  //gen branches
  tree_->Branch("B_Z_GenMuonPt", &B_Z_GenMuonPt);
  tree_->Branch("B_Z_GenMuonEta", &B_Z_GenMuonEta);
  tree_->Branch("B_Z_GenMuonPhi", &B_Z_GenMuonPhi);
  tree_->Branch("B_J_GenMuonPt", &B_J_GenMuonPt);
  tree_->Branch("B_J_GenMuonEta", &B_J_GenMuonEta);
  tree_->Branch("B_J_GenMuonPhi", &B_J_GenMuonPhi);
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
  tree_->Branch("FourL_VtxProb", &FourL_VtxProb);
  tree_->Branch("FourL_PVx", &FourL_PVx);
  tree_->Branch("FourL_PVy", &FourL_PVy);
  tree_->Branch("FourL_PVz", &FourL_PVz);
  tree_->Branch("FourL_PVxError", &FourL_PVxError);
  tree_->Branch("FourL_PVyError", &FourL_PVyError);
  tree_->Branch("FourL_PVzError", &FourL_PVzError);

  //tree_->Branch("B_Z_dca", &B_Z_dca);
  tree_->Branch("Ele_TriggerPath", &Ele_TriggerPath);
  tree_->Branch("B_Z_TriggerPt1", &B_Z_TriggerPt1);
  tree_->Branch("B_Z_TriggerEta1", &B_Z_TriggerEta1);
  tree_->Branch("B_Z_TriggerPhi1", &B_Z_TriggerPhi1);
  tree_->Branch("B_Z_TriggerPt2", &B_Z_TriggerPt2);
  tree_->Branch("B_Z_TriggerEta2", &B_Z_TriggerEta2);
  tree_->Branch("B_Z_TriggerPhi2", &B_Z_TriggerPhi2);
  tree_->Branch("B_Z_TriggerPt3", &B_Z_TriggerPt3);
  tree_->Branch("B_Z_TriggerEta3", &B_Z_TriggerEta3);
  tree_->Branch("B_Z_TriggerPhi3", &B_Z_TriggerPhi3);
  tree_->Branch("B_Z_TriggerPt4", &B_Z_TriggerPt4);
  tree_->Branch("B_Z_TriggerEta4", &B_Z_TriggerEta4);
  tree_->Branch("B_Z_TriggerPhi4", &B_Z_TriggerPhi4);
  tree_->Branch("B_Z_TriggerPt5", &B_Z_TriggerPt5);
  tree_->Branch("B_Z_TriggerEta5", &B_Z_TriggerEta5);
  tree_->Branch("B_Z_TriggerPhi5", &B_Z_TriggerPhi5);

  tree_->Branch("B_Z_Trigger32Pt1", &B_Z_Trigger32Pt1);
  tree_->Branch("B_Z_Trigger32Eta1", &B_Z_Trigger32Eta1);
  tree_->Branch("B_Z_Trigger32Phi1", &B_Z_Trigger32Phi1);
  tree_->Branch("B_Z_Trigger32Pt2", &B_Z_Trigger32Pt2);
  tree_->Branch("B_Z_Trigger32Eta2", &B_Z_Trigger32Eta2);
  tree_->Branch("B_Z_Trigger32Phi2", &B_Z_Trigger32Phi2);
  tree_->Branch("B_Z_Trigger32Pt3", &B_Z_Trigger32Pt3);
  tree_->Branch("B_Z_Trigger32Eta3", &B_Z_Trigger32Eta3);
  tree_->Branch("B_Z_Trigger32Phi3", &B_Z_Trigger32Phi3);
  tree_->Branch("B_Z_Trigger32Pt4", &B_Z_Trigger32Pt4);
  tree_->Branch("B_Z_Trigger32Eta4", &B_Z_Trigger32Eta4);
  tree_->Branch("B_Z_Trigger32Phi4", &B_Z_Trigger32Phi4);
  tree_->Branch("B_Z_Trigger32Pt5", &B_Z_Trigger32Pt5);
  tree_->Branch("B_Z_Trigger32Eta5", &B_Z_Trigger32Eta5);
  tree_->Branch("B_Z_Trigger32Phi5", &B_Z_Trigger32Phi5);

  tree_->Branch("B_Z_mass", &B_Z_mass);
  tree_->Branch("B_Z_VtxProb", &B_Z_VtxProb);

  tree_->Branch("B_Z_px", &B_Z_px);
  tree_->Branch("B_Z_py", &B_Z_py);
  tree_->Branch("B_Z_pz", &B_Z_pz);
  tree_->Branch("B_Z_pt", &B_Z_pt);
  tree_->Branch("B_Z_eta", &B_Z_eta);
  tree_->Branch("B_Z_phi", &B_Z_phi);
  tree_->Branch("B_Z_rapidity", &B_Z_rapidity);
  tree_->Branch("B_Z_VtxPx", &B_Z_VtxPx);
  tree_->Branch("B_Z_VtxPy", &B_Z_VtxPy);
  tree_->Branch("B_Z_VtxPz", &B_Z_VtxPz);
  tree_->Branch("B_Z_VtxPt", &B_Z_VtxPt);
  tree_->Branch("B_Z_VtxEta", &B_Z_VtxEta);
  tree_->Branch("B_Z_VtxPhi", &B_Z_VtxPhi);
  tree_->Branch("B_Z_VtxRapidity", &B_Z_VtxRapidity);
  tree_->Branch("B_Z_VtxMass", &B_Z_VtxMass);

  tree_->Branch("B_Z_PVx", &B_Z_PVx);
  tree_->Branch("B_Z_PVy", &B_Z_PVy);
  tree_->Branch("B_Z_PVz", &B_Z_PVz);
  tree_->Branch("B_Z_PVxError", &B_Z_PVxError);
  tree_->Branch("B_Z_PVyError", &B_Z_PVyError);
  tree_->Branch("B_Z_PVzError", &B_Z_PVzError);

  tree_->Branch("B_Z_px1", &B_Z_px1);
  tree_->Branch("B_Z_py1", &B_Z_py1);
  tree_->Branch("B_Z_pz1", &B_Z_pz1);
  tree_->Branch("B_Z_pt1", &B_Z_pt1);
  tree_->Branch("B_Z_eta1", &B_Z_eta1);
  tree_->Branch("B_Z_SCeta1", &B_Z_SCeta1);
  tree_->Branch("B_Z_phi1", &B_Z_phi1);
  tree_->Branch("B_Z_energy1", &B_Z_energy1);
  tree_->Branch("B_Z_energyCorr1", &B_Z_energyCorr1);
  tree_->Branch("B_Z_ecalIso1", &B_Z_ecalIso1);
  tree_->Branch("B_Z_hcalIso1", &B_Z_hcalIso1);
  tree_->Branch("B_Z_trackIso1", &B_Z_trackIso1);
  tree_->Branch("B_Z_looseCutBase1", &B_Z_looseCutBase1);
  tree_->Branch("B_Z_mediumCutBase1", &B_Z_mediumCutBase1);
  tree_->Branch("B_Z_tightCutBase1", &B_Z_tightCutBase1);
  tree_->Branch("B_Z_mvaIsoWP80_1", &B_Z_mvaIsoWP80_1);
  tree_->Branch("B_Z_mvaIsoWP90_1", &B_Z_mvaIsoWP90_1);
  tree_->Branch("B_Z_charge1", &B_Z_charge1);
  tree_->Branch("B_Z_px2", &B_Z_px2);
  tree_->Branch("B_Z_py2", &B_Z_py2);
  tree_->Branch("B_Z_pz2", &B_Z_pz2);
  tree_->Branch("B_Z_pt2", &B_Z_pt2);
  tree_->Branch("B_Z_eta2", &B_Z_eta2);
  tree_->Branch("B_Z_SCeta2", &B_Z_SCeta2);
  tree_->Branch("B_Z_phi2", &B_Z_phi2);
  tree_->Branch("B_Z_energy2", &B_Z_energy2);
  tree_->Branch("B_Z_energyCorr2", &B_Z_energyCorr2);
  tree_->Branch("B_Z_ecalIso2", &B_Z_ecalIso2);
  tree_->Branch("B_Z_hcalIso2", &B_Z_hcalIso2);
  tree_->Branch("B_Z_trackIso2", &B_Z_trackIso2);
  tree_->Branch("B_Z_looseCutBase2", &B_Z_looseCutBase2);
  tree_->Branch("B_Z_mediumCutBase2", &B_Z_mediumCutBase2);
  tree_->Branch("B_Z_tightCutBase2", &B_Z_tightCutBase2);
  tree_->Branch("B_Z_mvaIsoWP80_2", &B_Z_mvaIsoWP80_2);
  tree_->Branch("B_Z_mvaIsoWP90_2", &B_Z_mvaIsoWP90_2);
  tree_->Branch("B_Z_charge2", &B_Z_charge2);
  tree_->Branch("B_Z_dxy1", &B_Z_dxy1);
  tree_->Branch("B_Z_dxy2", &B_Z_dxy2);
  tree_->Branch("B_Z_dz1", &B_Z_dz1);
  tree_->Branch("B_Z_dz2", &B_Z_dz2);

  tree_->Branch("B_J_lowPt", &B_J_lowPt);
  tree_->Branch("B_J_highPt", &B_J_highPt);
  tree_->Branch("B_Z_lowPt", &B_Z_lowPt);
  tree_->Branch("B_Z_highPt", &B_Z_highPt);

  tree_->Branch("B_J_dca", &B_J_dca);
  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);
  tree_->Branch("B_J_pt", &B_J_pt);
  tree_->Branch("B_J_eta", &B_J_eta);
  tree_->Branch("B_J_phi", &B_J_phi);
  tree_->Branch("B_J_rapidity", &B_J_rapidity);
  tree_->Branch("B_J_VtxPx", &B_J_VtxPx);
  tree_->Branch("B_J_VtxPy", &B_J_VtxPy);
  tree_->Branch("B_J_VtxPz", &B_J_VtxPz);
  tree_->Branch("B_J_VtxPt", &B_J_VtxPt);
  tree_->Branch("B_J_VtxEta", &B_J_VtxEta);
  tree_->Branch("B_J_VtxPhi", &B_J_VtxPhi);
  tree_->Branch("B_J_VtxRapidity", &B_J_VtxRapidity);
  tree_->Branch("B_J_VtxMass", &B_J_VtxMass);
  tree_->Branch("B_J_PVx", &B_J_PVx);
  tree_->Branch("B_J_PVy", &B_J_PVy);
  tree_->Branch("B_J_PVz", &B_J_PVz);
  tree_->Branch("B_J_PVxError", &B_J_PVxError);
  tree_->Branch("B_J_PVyError", &B_J_PVyError);
  tree_->Branch("B_J_PVzError", &B_J_PVzError);

  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_pt1", &B_J_pt1);
  tree_->Branch("B_J_eta1", &B_J_eta1);
  tree_->Branch("B_J_phi1", &B_J_phi1);
  tree_->Branch("B_J_charge1", &B_J_charge1);
  tree_->Branch("B_J_soft1", &B_J_soft1);
  tree_->Branch("B_J_tight1", &B_J_tight1);
  tree_->Branch("B_J_loose1", &B_J_loose1);

  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_pt2", &B_J_pt2);
  tree_->Branch("B_J_eta2", &B_J_eta2);
  tree_->Branch("B_J_phi2", &B_J_phi2);
  tree_->Branch("B_J_charge2", &B_J_charge2);
  tree_->Branch("B_J_soft2", &B_J_soft2);
  tree_->Branch("B_J_tight2", &B_J_tight2);
  tree_->Branch("B_J_loose2", &B_J_loose2);
  tree_->Branch("B_J_VtxProb", &B_J_VtxProb);
  tree_->Branch("B_J_xyP", &B_J_xyP);
  tree_->Branch("B_J_xyM", &B_J_xyM);
  tree_->Branch("B_J_zP", &B_J_zP);
  tree_->Branch("B_J_zM", &B_J_zM);

  tree_->Branch("mumC2", &mumC2);
  tree_->Branch("mumNHmvaisoWP90 on both electron hits", &mumNHits);
  tree_->Branch("mumNPHits", &mumNPHits);
  tree_->Branch("mupC2", &mupC2);
  tree_->Branch("mupNHits", &mupNHits);
  tree_->Branch("mupNPHits", &mupNPHits);
}

// ------------ method called once each job just after ending the event loop  ------------
void miniAODeemm::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(miniAODeemm);
