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

int nuuuu = 0;
float isZuuZuu = 0;
int counter = 0;

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

      isMC_(iConfig.getParameter<bool>("isMC")),

      tree_(0),
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
      B_U_TriggerPath(0),
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
      B_J2_mass(0),
      B_J2_px(0),
      B_J2_py(0),
      B_J2_pz(0),
      B_J2_pt(0),
      B_J2_eta(0),
      B_J2_phi(0),
      B_J2_rapidity(0),
      B_J3_mass(0),
      B_J3_px(0),
      B_J3_py(0),
      B_J3_pz(0),
      B_J3_pt(0),
      B_J3_eta(0),
      B_J3_phi(0),
      B_J3_rapidity(0),
      B_J4_mass(0),
      B_J4_px(0),
      B_J4_py(0),
      B_J4_pz(0),
      B_J4_pt(0),
      B_J4_eta(0),
      B_J4_phi(0),
      B_J4_rapidity(0),

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
      B_J2_VtxPx(0),
      B_J2_VtxPy(0),
      B_J2_VtxPz(0),
      B_J2_VtxPt(0),
      B_J2_VtxEta(0),
      B_J2_VtxPhi(0),
      B_J2_VtxRapidity(0),
      B_J2_VtxMass(0),
      B_J2_PVx(0),
      B_J2_PVy(0),
      B_J2_PVz(0),
      B_J2_PVxError(0),
      B_J2_PVyError(0),
      B_J2_PVzError(0),
      B_J3_VtxPx(0),
      B_J3_VtxPy(0),
      B_J3_VtxPz(0),
      B_J3_VtxPt(0),
      B_J3_VtxEta(0),
      B_J3_VtxPhi(0),
      B_J3_VtxRapidity(0),
      B_J3_VtxMass(0),
      B_J3_PVx(0),
      B_J3_PVy(0),
      B_J3_PVz(0),
      B_J3_PVxError(0),
      B_J3_PVyError(0),
      B_J3_PVzError(0),
      B_J4_VtxPx(0),
      B_J4_VtxPy(0),
      B_J4_VtxPz(0),
      B_J4_VtxPt(0),
      B_J4_VtxEta(0),
      B_J4_VtxPhi(0),
      B_J4_VtxRapidity(0),
      B_J4_VtxMass(0),
      B_J4_PVx(0),
      B_J4_PVy(0),
      B_J4_PVz(0),
      B_J4_PVxError(0),
      B_J4_PVyError(0),
      B_J4_PVzError(0),

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

      B_Mu2_px(0),
      B_Mu2_py(0),
      B_Mu2_pz(0),
      B_Mu2_pt(0),
      B_Mu2_eta(0),
      B_Mu2_phi(0),
      B_Mu1_charge(0),
      B_Mu2_charge(0),
      B_Mu2_soft(0),
      B_Mu2_tight(0),
      B_Mu2_loose(0),
      B_Mu2_IsoTrack(0),
      B_Mu2_IsoHcal(0),
      B_Mu2_IsoEcal(0),
      B_Mu2_IsoCalo(0),

      B_Mu2_PaperIsoTrackRF04(0),
      B_Mu2_PaperIsoTrackRF03(0),
      B_Mu2_Paper3DIP(0),

      B_Mu3_px(0),
      B_Mu3_py(0),
      B_Mu3_pz(0),
      B_Mu3_pt(0),
      B_Mu3_eta(0),
      B_Mu3_phi(0),
      B_Mu3_soft(0),
      B_Mu3_tight(0),
      B_Mu3_loose(0),
      B_Mu3_charge(0),
      B_Mu3_IsoTrack(0),
      B_Mu3_IsoHcal(0),
      B_Mu3_IsoEcal(0),
      B_Mu3_IsoCalo(0),

      B_Mu3_PaperIsoTrackRF04(0),
      B_Mu3_PaperIsoTrackRF03(0),
      B_Mu3_Paper3DIP(0),

      B_Mu4_px(0),
      B_Mu4_py(0),
      B_Mu4_pz(0),
      B_Mu4_pt(0),
      B_Mu4_eta(0),
      B_Mu4_phi(0),
      B_Mu4_soft(0),
      B_Mu4_tight(0),
      B_Mu4_loose(0),
      B_Mu4_charge(0),
      B_Mu4_IsoTrack(0),
      B_Mu4_IsoHcal(0),
      B_Mu4_IsoEcal(0),
      B_Mu4_IsoCalo(0),

      B_Mu4_PaperIsoTrackRF04(0),
      B_Mu4_PaperIsoTrackRF03(0),
      B_Mu4_Paper3DIP(0),

      B_J1_VtxProb(0),
      B_J2_VtxProb(0),
      B_J3_VtxProb(0),
      B_J4_VtxProb(0),
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

      mu2mC2(0),
      mu2mNHits(0),
      mu2mNPHits(0),
      mu2pC2(0),
      mu2pNHits(0),
      mu2pNPHits(0),

      B_M1_pt(0),
      B_M1_eta(0),
      B_M1_phi(0),
      B_M1_px(0),
      B_M1_py(0),
      B_M1_pz(0),
      B_M2_pt(0),
      B_M2_eta(0),
      B_M2_phi(0),
      B_M2_px(0),
      B_M2_py(0),
      B_M2_pz(0),
      B_M3_pt(0),
      B_M3_eta(0),
      B_M3_phi(0),
      B_M3_px(0),
      B_M3_py(0),
      B_M3_pz(0),
      B_M4_pt(0),
      B_M4_eta(0),
      B_M4_phi(0),
      B_M4_px(0),
      B_M4_py(0),
      B_M4_pz(0),
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
  int pass = 0;
  int passTrig = 0;
  bool EET = false;
  float EET_pt[5] = {-999, -999, -999, -999, -999};
  float EET_eta[5] = {-999, -999, -999, -999, -999};
  float EET_phi[5] = {-999, -999, -999, -999, -999};
  int nG1 = 0;
  int nG2 = 0;
  bool GenInfo = false;
  if (GenInfo == true) {
    //Gen Level Info

    float B_J_GenMuon_pt = -999;
    float B_J_GenMuon_eta = -999;
    float B_J_GenMuon_phi = -999;
    float B_Z_GenMuon_pt = -999;
    float B_Z_GenMuon_eta = -999;
    //float B_Z_GenMuon_phi=-999;

    //cout<<"Start Gen Work"<<endl;
    for (size_t i = 0; i < pruned->size(); i++) {
      if ((*pruned)[i].isPromptFinalState() && abs((*pruned)[i].pdgId()) == 13) {
        //if( abs((*pruned)[i].pdgId()) ==13 ){
        //cout<<"Found gen level muon"<<endl;
        //if ( (*pruned)[i].mother()->pdgId()==553 ) {
        //cout<<"Have Ups1 in Gen"<<endl;
        B_J_GenMuon_pt = (*pruned)[i].pt();
        B_J_GenMuon_eta = (*pruned)[i].eta();
        B_J_GenMuon_phi = (*pruned)[i].phi();
        B_J_GenMuonPt->push_back(B_J_GenMuon_pt);
        B_J_GenMuonEta->push_back(B_J_GenMuon_eta);
        B_J_GenMuonPhi->push_back(B_J_GenMuon_phi);
        //cout<<"Pt of Ups1 in Gen"<<(*pruned)[i].pt()<<endl;
        nG1++;
        //}
      }
      if (abs((*pruned)[i].pdgId()) == 13) {
        //cout<<"Found gen level muon"<<endl;
        if ((*pruned)[i].mother()->pdgId() == 23) {
          B_Z_GenMuon_pt = (*pruned)[i].pt();
          B_Z_GenMuon_eta = (*pruned)[i].eta();
          //B_Z_GenMuon_phi = (*pruned)[i].phi();
          B_Z_GenMuonPt->push_back(B_Z_GenMuon_pt);
          B_Z_GenMuonEta->push_back(B_Z_GenMuon_eta);
          B_Z_GenMuonPhi->push_back(nuuuu);
          nG2++;
        }
      }
    }
    //cout<<"End Gen Work"<<endl;
    isZuuZuu = 0;
    counter++;
    if (nG2 == 4) {
      isZuuZuu = 1;
      nuuuu++;
      cout << "Number of ZZ " << nuuuu << endl;
      cout << "Number of events " << counter << endl;
    }
  }

  //First Work on Trigger Information//
  //**********************************

  const edm::TriggerNames& names = iEvent.triggerNames(*triggerBits);
  //std::cout << "\n === TRIGGER PATHS === " << std::endl;

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    bool acceptE = triggerBits->accept(i);
    if (names.triggerName(i).find("HLT_IsoMu24_v") != string::npos) {
      if (triggerBits->accept(i)) {
        //std::cout << "Trigger " << names.triggerName(i) <<
        // ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
        //        << std::endl;
        pass = acceptE;
        if (pass == 1) {
          passTrig++;
        }
      }
    }
    /*
    if (names.triggerName(i).find("HLT_IsoMu20_v")!=string::npos)  {                                                                                                                                           
      if (triggerBits->accept(i)) { 
        pass=acceptE;                                                                                                                                                                                       
        if (pass==1){                                                                                                                                                                                           
         passTrig++; 
        }
      } 
    }
    
    if (names.triggerName(i).find("HLT_IsoMu22_v")!=string::npos)  {                                                                                                                                                     if (triggerBits->accept(i)) {                                                                                                                                                                                        pass=acceptE;                                                                                                                                                                                                      if (pass==1){                                                                                                                                                                                                       passTrig++;                                                                                                                                                                                                       }                                                                                                                                                                                                                }                                                                                                                                                                                                                }
    if (names.triggerName(i).find("HLT_IsoMu24_v")!=string::npos)  {                                                                                                                                                     if (triggerBits->accept(i)) {                                                                                                                                                                                        pass=acceptE;                                                                                                                                                                                                      if (pass==1){                                                                                                                                                                                                       passTrig++;                                                                                                                                                                                                       }                                                                                                                                                                                                                }                                                                                                                                                                                                                }
    if (names.triggerName(i).find("HLT_Mu50_v")!=string::npos)  {                                                                                                                                                     if (triggerBits->accept(i)) {                                                                                                                                                                                        pass=acceptE;                                                                                                                                                                                                      if (pass==1){                                                                                                                                                                                                       passTrig++;                                                                                                                                                                                                       }                                                                                                                                                                                                                }                                                                                                                                                                                                                }
    if (names.triggerName(i).find("HLT_Mu55_v")!=string::npos)  {                                                                                                                                                     if (triggerBits->accept(i)) {                                                                                                                                                                                        pass=acceptE;                                                                                                                                                                                                      if (pass==1){                                                                                                                                                                                                       passTrig++;                                                                                                                                                                                                       }                                                                                                                                                                                                                }                                                                                                                                                                                                                }
*/
  }
  //****************************************************
  //***************Trigger Object***********************
  //****************************************************
  if (passTrig == 0) {
    //    cout<<"trigger didnt work "<<endl;
  }
  if (passTrig > 0) {
    EET = true;
    //float ElectronTriggerPt;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {  // note: not "const &" since we want to call unpackPathNames
      obj.unpackPathNames(names);
      //for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
      //bool isBoth = obj.hasPathName("HLT_Ele35_WPTight_Gsf_v*", true, true );
      //bool isL3 = obj.hasFilterLabel("HLTEle32WPTightGsfSequence");
      //cout<<"The L3 Filter is : "<<isL3<<endl;
      //bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
      //bool isL3   = obj.hasPathName( "HLT_Ele35_WPTight_Gsf_v*", false, true );
      //bool isLF   = obj.hasPathName( "HLT_Ele35_WPTight_Gsf_v*", true, false );
      //bool isNone = obj.hasPathName( "HLT_Ele35_WPTight_Gsf_v*", false, false );
      //**************************************************************************************************/
      //****************Definition of hasPathName*********************************************************//
      //**************************************************************************************************//
      //bool hasPathName(const std::string &pathName,
      //               bool pathLastFilterAccepted = false,
      //               bool pathL3FilterAccepted = true) const {
      //return hasPathOrAlgorithm(pathName, pathLastFilterAccepted, pathL3FilterAccepted);
      //};
      //**************************************************************************************************//
      //cout<< "obj.hasPathL3FilterAccepted() "<<obj.hasPathL3FilterAccepted()<<endl;
      //std::cout << "   " << pathNamesAll[h];

      //if (obj.hasPathName("HLT_IsoMu27_vfake", true, true )>0||obj.hasPathName("HLT_IsoMu24_vfake", true, true )>0||obj.hasPathName("HLT_IsoMu22_vfake", true, true )>0||obj.hasPathName("HLT_IsoMu20_vfake", true, true )>0||obj.hasPathName("HLT_Mu55_vfake", true, true )>0||obj.hasPathName("HLT_Mu50_vfake", true, true )>0) {
      if (obj.hasPathName("HLT_IsoMu24_v", true, true) > 0) {
        //std::cout << "\tTrigger objectisBoth:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
        //open from here
        EET_pt[h] = obj.pt();
        EET_eta[h] = obj.eta();
        EET_phi[h] = obj.phi();
        h++;
      }
      //if (isL3>0) {
      //std::cout << "\tTrigger objectisL3:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      //}
      //if (isLF>0) {
      //      std::cout << "\tTrigger objectisLF:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      //}
      //if (isNone>0) {
      //std::cout << "\tTrigger objectNone :  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      //}
      //h++;
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
      for (View<pat::Muon>::const_iterator iMuon3 = iMuon2 + 1; iMuon3 != thePATMuonHandle->end(); ++iMuon3) {
        for (View<pat::Muon>::const_iterator iMuon4 = iMuon3 + 1; iMuon4 != thePATMuonHandle->end(); ++iMuon4) {
          //make sure all muons are diferent
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
          //opposite charge
          //only look for neutral charge combination

          if (!(abs(iMuon1->charge()) == 1))
            continue;
          if (!(abs(iMuon2->charge()) == 1))
            continue;
          if (!(abs(iMuon3->charge()) == 1))
            continue;
          if (!(abs(iMuon4->charge()) == 1))
            continue;
          //	  cout<<"paperIsolationR04: "<<(iMuon1->pfIsolationR03().sumChargedHadronPt+ std::max(0., iMuon1->pfIsolationR03().sumNeutralHadronEt+ iMuon1->pfIsolationR03().sumPhotonEt - iMuon1->pfIsolationR03().sumPUPt*0.5))/iMuon1->pt()<<endl;
          //cout<<"charge of 2nd muon"<<iMuon2->charge()<<endl;
          //cout<<"charge of 3rd muon"<<iMuon3->charge()<<endl;
          //cout<<"charge of 4th muon"<<iMuon4->charge()<<endl;

          if (!((iMuon1->charge()) + (iMuon2->charge()) + (iMuon3->charge()) + (iMuon4->charge()) == 0))
            continue;
          //cout<<"Addition of 4 muon charge "<<(iMuon1->charge()) + (iMuon2->charge()) + (iMuon3->charge()) + (iMuon4->charge())<<endl;
          //preselection ptcut
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

          if (iMuon1->charge() == 1) {
            glbTrackP1 = iMuon1->track();
            if (iMuon2->charge() == 1) {
              glbTrackP2 = iMuon2->track();
              glbTrackM1 = iMuon3->track();
              glbTrackM2 = iMuon4->track();
            } else if (iMuon3->charge() == 1) {
              glbTrackP2 = iMuon3->track();
              glbTrackM1 = iMuon2->track();
              glbTrackM2 = iMuon4->track();
            } else if (iMuon4->charge() == 1) {
              glbTrackP2 = iMuon4->track();
              glbTrackM1 = iMuon2->track();
              glbTrackM2 = iMuon3->track();
            } else {
              cout << "Something is wrong while making +1 glb track ref" << endl;
            }
          }

          if (iMuon1->charge() == -1) {
            glbTrackM1 = iMuon1->track();
            if (iMuon2->charge() == -1) {
              glbTrackM2 = iMuon2->track();
              glbTrackP1 = iMuon3->track();
              glbTrackP2 = iMuon4->track();
            } else if (iMuon3->charge() == -1) {
              glbTrackM2 = iMuon3->track();
              glbTrackP1 = iMuon2->track();
              glbTrackP2 = iMuon4->track();

            } else if (iMuon4->charge() == -1) {
              glbTrackM2 = iMuon4->track();
              glbTrackP1 = iMuon2->track();
              glbTrackP2 = iMuon3->track();

            } else {
              cout << "Something is wrong while making -1 glb track ref" << endl;
            }
          }

          if (glbTrackP1.isNull() || glbTrackM1.isNull() || glbTrackP2.isNull() || glbTrackM2.isNull()) {
            //std::cout << "continue due to no track ref" << endl;
            continue;
          }

          TLorentzVector M1, M2, M3, M4, MM1, MM2, MM3, MM4, MMMM;
          //initialize 4 lepton mass
          float mu_mass = 0.1056583745;  //[PDG mass]
          //float ele_mass =  0.000510998928;//PDG mass

          //make muon 4 vectors

          if (iMuon1->charge() == 1) {
            if (iMuon2->charge() == 1) {
              //P1->Mu1;P2->Mu2;M1->Mu3;M2->Mu4
              M1.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M3.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M2.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M4.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            } else if (iMuon3->charge() == 1) {
              //P1->Mu1;P2->Mu3;M1->Mu2;M2->Mu4
              M1.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M2.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M3.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M4.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            } else if (iMuon4->charge() == 1) {
              //P1->Mu1;P2->Mu4;M1->Mu2;M2->Mu3
              M1.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M2.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M4.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M3.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            } else {
              cout << "Something is wrong while making +1 charge lorentz vector" << endl;
            }
          } else if (iMuon1->charge() == -1) {
            if (iMuon2->charge() == -1) {
              //P1->Mu3;P2->Mu4;M1->Mu1;M2->Mu2
              M3.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M1.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M4.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M2.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            } else if (iMuon3->charge() == -1) {
              //P1->Mu2;P2->Mu4;M1->Mu1;M2->Mu3
              M2.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M1.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M4.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M3.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            } else if (iMuon4->charge() == -1) {
              //P1->Mu2;P2->Mu3;M1->Mu1;M2->Mu4
              M2.SetXYZM(iMuon1->px(), iMuon1->py(), iMuon1->pz(), mu_mass);
              M1.SetXYZM(iMuon2->px(), iMuon2->py(), iMuon2->pz(), mu_mass);
              M3.SetXYZM(iMuon3->px(), iMuon3->py(), iMuon3->pz(), mu_mass);
              M4.SetXYZM(iMuon4->px(), iMuon4->py(), iMuon4->pz(), mu_mass);
            } else {
              cout << "Something is wrong while making -1 charge lorentz vector" << endl;
            }
          }
          //now M1 is first positive muon 4 vector
          //then M2 is first negative muon 4 vector
          //then M3 is second positive muon 4 vector
          //finally, M4 is second negative muon
          //****************************************************************************************
          //make netral dimuon combination
          MM1 = M1 + M2;
          MM2 = M3 + M4;
          MM3 = M2 + M3;
          MM4 = M1 + M4;
          //****************************************************************************************
          //check all the combinetion again
          //if (MM.M()<2.8) continue;
          //if (MM.M()>3.4) continue;
          //if (EE.M()<70) continue;
          //if (EE.M()>110) continue;

          MMMM = M1 + M2 + M3 + M4;

          //cout<<"Start looking muon track quality"<<endl;
          if (!(glbTrackM1->quality(reco::TrackBase::highPurity)))
            continue;
          if (!(glbTrackP1->quality(reco::TrackBase::highPurity)))
            continue;
          if (!(glbTrackM2->quality(reco::TrackBase::highPurity)))
            continue;
          if (!(glbTrackP2->quality(reco::TrackBase::highPurity)))
            continue;
          //cout<<"Start Building Track"<<endl;

          reco::TransientTrack muon1TT((*theB).build(glbTrackP1));
          reco::TransientTrack muon2TT((*theB).build(glbTrackM1));
          reco::TransientTrack muon3TT((*theB).build(glbTrackP2));
          reco::TransientTrack muon4TT((*theB).build(glbTrackM2));

          // *****  Trajectory states to calculate DCA for the 2 muons *********************
          //FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
          //FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();
          //FreeTrajectoryState mu3State = muon3TT.impactPointTSCP().theState();
          //FreeTrajectoryState mu4State = muon4TT.impactPointTSCP().theState();

          //cout<<"Start validating impact point"<<endl;
          //if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() || !muon3TT.impactPointTSCP().isValid() || !muon4TT.impactPointTSCP().isValid() ) continue;
          //cout<<"Validated  impact point for muon"<<endl;

          // Measure distance between tracks at their closest approach
          //ClosestApproachInRPhi cApp,cApp1,cApp2,cApp3;
          //cApp.calculate(mu1State, mu2State);
          //cApp1.calculate(mu3State, mu4State);
          //cApp2.calculate(mu2State, mu3State);
          //cApp3.calculate(mu1State, mu4State);
          //if( !cApp1.status() ) continue;
          //if( !cApp2.status() ) continue;
          //if( !cApp3.status() ) continue;
          //if( !cApp4.status() ) continue;
          //float dca = fabs( cApp1.distance() );
          //if (dca < 0. || dca > 0.5) continue;
          //cout<<" closest approach  "<<dca<<endl;

          //Kalman Vtx----------------------//
          vector<TransientTrack> mu_tks;
          vector<TransientTrack> mmmm_tks;
          KalmanVertexFitter kvf(true);

          KalmanVertexFitter kvfM(true);
          //First neutral combination
          mu_tks.clear();
          mu_tks.push_back(muon1TT);
          mu_tks.push_back(muon2TT);
          TransientVertex J_candi1 = kvfM.vertex(mu_tks);
          //second combination
          mu_tks.clear();
          mu_tks.push_back(muon3TT);
          mu_tks.push_back(muon4TT);
          TransientVertex J_candi2 = kvfM.vertex(mu_tks);
          //3rd combination
          mu_tks.clear();
          mu_tks.push_back(muon2TT);
          mu_tks.push_back(muon3TT);
          TransientVertex J_candi3 = kvfM.vertex(mu_tks);
          //4th combination
          mu_tks.clear();
          mu_tks.push_back(muon1TT);
          mu_tks.push_back(muon4TT);
          TransientVertex J_candi4 = kvfM.vertex(mu_tks);

          int psi_candi = 0;

          if (J_candi1.isValid()) {
            //cout<<"J candidate non validated by kalman fitter"<<endl;
            if (J_candi2.isValid()) {
              //cout<<"J candidate non validated by kalman fitter"<<endl;
              psi_candi++;
            }
          }
          if (J_candi3.isValid()) {
            //cout<<"J candidate non validated by kalman fitter"<<endl;
            if (J_candi4.isValid()) {
              //cout<<"J candidate non validated by kalman fitter"<<endl;
              psi_candi++;
            }
          }
          if (psi_candi < 1) {
            cout << "continue because no more than 2 vertexed dimuon" << endl;
            continue;
          }

          reco::Vertex JPsi_Vtx1 = J_candi1;
          reco::Vertex JPsi_Vtx2 = J_candi2;
          reco::Vertex JPsi_Vtx3 = J_candi3;
          reco::Vertex JPsi_Vtx4 = J_candi4;

          float B_Prob_tmp1 = TMath::Prob(J_candi1.totalChiSquared(), J_candi1.degreesOfFreedom());
          const math::XYZTLorentzVectorD JPsi_mom1 = JPsi_Vtx1.p4(mu_mass, 0.0);
          float B_Prob_tmp2 = TMath::Prob(J_candi2.totalChiSquared(), J_candi2.degreesOfFreedom());
          const math::XYZTLorentzVectorD JPsi_mom2 = JPsi_Vtx2.p4(mu_mass, 0.0);
          float B_Prob_tmp3 = TMath::Prob(J_candi3.totalChiSquared(), J_candi3.degreesOfFreedom());
          const math::XYZTLorentzVectorD JPsi_mom3 = JPsi_Vtx3.p4(mu_mass, 0.0);
          float B_Prob_tmp4 = TMath::Prob(J_candi4.totalChiSquared(), J_candi4.degreesOfFreedom());
          const math::XYZTLorentzVectorD JPsi_mom4 = JPsi_Vtx4.p4(mu_mass, 0.0);

          int DimuonVtx = 0;
          if (B_Prob_tmp1 > 0.001) {
            if (B_Prob_tmp2 > 0.001) {
              DimuonVtx++;
            }
          }
          if (B_Prob_tmp3 > 0.001) {
            if (B_Prob_tmp4 < 0.001) {
              DimuonVtx++;
            }
          }
          //cout << "Dimuon vetexed number" << DimuonVtx << endl;
          if (DimuonVtx < 1)
            continue;
          KalmanVertexFitter kvf4M(true);
          mmmm_tks.clear();
          mmmm_tks.push_back(muon1TT);
          mmmm_tks.push_back(muon2TT);
          mmmm_tks.push_back(muon3TT);
          mmmm_tks.push_back(muon4TT);

          TransientVertex FourL_candi = kvf4M.vertex(mmmm_tks);
          if (!FourL_candi.isValid()) {
            cout << "4 mu candidate non validated by kalman fitter" << endl;
            continue;
          }
          float B_Prob_tmp4L = TMath::Prob(FourL_candi.totalChiSquared(), FourL_candi.degreesOfFreedom());
          reco::Vertex FourL_Vtx = FourL_candi;

          if (B_Prob_tmp4L < 0.001) {
            //cout<<"4l vertexing failed"<<endl;
            continue;
          }

          //Final check for # of ups*****
          // edit ** oldmuon cut 8.0-12.0 ** Jesse Harris
          // edit ** widemuon cut 2.8-110 ** Jesse Harris
          int UpsNumber = 0;
          if ((MM1.M() > 2.8 && MM1.M() < 12) || (MM1.M() > 70 && MM1.M() < 110)) {
            if ((MM2.M() > 2.8 && MM2.M() < 12) || (MM2.M() > 70 && MM1.M() < 110)) {
              if ((MM1.M() > 2.8 && MM1.M() < 12) && (MM2.M() > 1 && MM2.M() < 12)) {
              } else {
                UpsNumber++;
              }
              if ((MM1.M() > 70 && MM1.M() < 110) && (MM2.M() > 70 && MM2.M() < 110)) {
              } else {
                UpsNumber++;
              }
            }
          }
          if ((MM3.M() > 2.8 && MM3.M() < 12) || (MM3.M() > 70 && MM3.M() < 110)) {
            if ((MM4.M() > 2.8 && MM4.M() < 12) || (MM4.M() > 70 && MM4.M() < 110)) {
              if ((MM3.M() > 2.8 && MM3.M() < 12) && (MM4.M() > 1 && MM4.M() < 12)) {
              } else {
                UpsNumber++;
              }
              if ((MM3.M() > 70 && MM3.M() < 110) && (MM4.M() > 70 && MM4.M() < 110)) {
              } else {
                UpsNumber++;
              }
            }
          }
          if (UpsNumber < 1)
            continue;
          //cout<<"Upsilon Mumber is : "<<UpsNumber<<endl;

          //*********************
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
              //std::cout << "negative chisq from psi fit" << endl;
              continue;
            }

            //some loose cuts go here

            // if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
            // if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;

            //fill variables?iMuon1->track()->pt()

            //B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
            //B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
            //B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
            //B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );
          }
          //cout<<JPsi_mom.mass()<<" : Is JPsi Vtx Mass "<<MM.M()<<"Is JPsi mass"<<endl;
          cout << iMuon1->dB(pat::Muon::PV3D) / iMuon1->edB(pat::Muon::PV3D) << endl;
          //Event Information
          Run->push_back(iEvent.id().run());
          LumiBlock->push_back(iEvent.luminosityBlock());
          Event->push_back(iEvent.id().event());

          FourL_mass->push_back(MMMM.M());
          FourL_px->push_back(MMMM.Px());
          FourL_py->push_back(MMMM.Py());
          FourL_pz->push_back(MMMM.Pz());
          FourL_pt->push_back(MMMM.Pt());
          FourL_eta->push_back(MMMM.Eta());
          FourL_phi->push_back(MMMM.Phi());
          FourL_VtxProb->push_back(B_Prob_tmp4L);
          FourL_PVx->push_back(FourL_Vtx.x());
          FourL_PVy->push_back(FourL_Vtx.y());
          FourL_PVz->push_back(FourL_Vtx.z());
          FourL_PVxError->push_back(FourL_Vtx.xError());
          FourL_PVyError->push_back(FourL_Vtx.yError());
          FourL_PVzError->push_back(FourL_Vtx.zError());

          B_U_TriggerPath->push_back(EET);

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
          B_J1_px->push_back(isZuuZuu);
          B_J1_py->push_back(MM1.Py());
          B_J1_pz->push_back(MM1.Pz());
          B_J1_pt->push_back(MM1.Pt());
          B_J1_eta->push_back(MM1.Eta());
          B_J1_phi->push_back(MM1.Phi());
          B_J1_rapidity->push_back(MM1.Rapidity());
          B_J2_mass->push_back(MM2.M());
          B_J2_px->push_back(MM2.Px());
          B_J2_py->push_back(MM2.Py());
          B_J2_pz->push_back(MM2.Pz());
          B_J2_pt->push_back(MM2.Pt());
          B_J2_eta->push_back(MM2.Eta());
          B_J2_phi->push_back(MM2.Phi());
          B_J2_rapidity->push_back(MM2.Rapidity());
          B_J3_mass->push_back(MM3.M());
          B_J3_px->push_back(MM3.Px());
          B_J3_py->push_back(MM3.Py());
          B_J3_pz->push_back(MM3.Pz());
          B_J3_pt->push_back(MM3.Pt());
          B_J3_eta->push_back(MM3.Eta());
          B_J3_phi->push_back(MM3.Phi());
          B_J3_rapidity->push_back(MM3.Rapidity());
          B_J4_mass->push_back(MM4.M());
          B_J4_px->push_back(MM4.Px());
          B_J4_py->push_back(MM4.Py());
          B_J4_pz->push_back(MM4.Pz());
          B_J4_pt->push_back(MM4.Pt());
          B_J4_eta->push_back(MM4.Eta());
          B_J4_phi->push_back(MM4.Phi());
          B_J4_rapidity->push_back(MM4.Rapidity());

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

          B_J2_VtxPx->push_back(JPsi_mom2.Px());
          B_J2_VtxPy->push_back(JPsi_mom2.Py());
          B_J2_VtxPz->push_back(JPsi_mom2.Pz());
          B_J2_VtxPt->push_back(JPsi_mom2.Pt());
          B_J2_VtxEta->push_back(JPsi_mom2.Eta());
          B_J2_VtxPhi->push_back(JPsi_mom2.Phi());
          B_J2_VtxRapidity->push_back(JPsi_mom2.Rapidity());
          B_J2_VtxMass->push_back(JPsi_mom2.mass());
          B_J2_PVx->push_back(JPsi_Vtx2.x());
          B_J2_PVy->push_back(JPsi_Vtx2.y());
          B_J2_PVz->push_back(JPsi_Vtx2.z());
          B_J2_PVxError->push_back(JPsi_Vtx2.xError());
          B_J2_PVyError->push_back(JPsi_Vtx2.yError());
          B_J2_PVzError->push_back(JPsi_Vtx2.zError());

          B_J3_VtxPx->push_back(JPsi_mom3.Px());
          B_J3_VtxPy->push_back(JPsi_mom3.Py());
          B_J3_VtxPz->push_back(JPsi_mom3.Pz());
          B_J3_VtxPt->push_back(JPsi_mom3.Pt());
          B_J3_VtxEta->push_back(JPsi_mom3.Eta());
          B_J3_VtxPhi->push_back(JPsi_mom3.Phi());
          B_J3_VtxRapidity->push_back(JPsi_mom3.Rapidity());
          B_J3_VtxMass->push_back(JPsi_mom3.mass());
          B_J3_PVx->push_back(JPsi_Vtx3.x());
          B_J3_PVy->push_back(JPsi_Vtx3.y());
          B_J3_PVz->push_back(JPsi_Vtx3.z());
          B_J3_PVxError->push_back(JPsi_Vtx3.xError());
          B_J3_PVyError->push_back(JPsi_Vtx3.yError());
          B_J3_PVzError->push_back(JPsi_Vtx3.zError());

          B_J4_VtxPx->push_back(JPsi_mom4.Px());
          B_J4_VtxPy->push_back(JPsi_mom4.Py());
          B_J4_VtxPz->push_back(JPsi_mom4.Pz());
          B_J4_VtxPt->push_back(JPsi_mom4.Pt());
          B_J4_VtxEta->push_back(JPsi_mom4.Eta());
          B_J4_VtxPhi->push_back(JPsi_mom4.Phi());
          B_J4_VtxRapidity->push_back(JPsi_mom4.Rapidity());
          B_J4_VtxMass->push_back(JPsi_mom4.mass());
          B_J4_PVx->push_back(JPsi_Vtx4.x());
          B_J4_PVy->push_back(JPsi_Vtx4.y());
          B_J4_PVz->push_back(JPsi_Vtx4.z());
          B_J4_PVxError->push_back(JPsi_Vtx4.xError());
          B_J4_PVyError->push_back(JPsi_Vtx4.yError());
          B_J4_PVzError->push_back(JPsi_Vtx4.zError());
          //dimuon vtx prob
          B_J1_VtxProb->push_back(B_Prob_tmp1);
          B_J2_VtxProb->push_back(B_Prob_tmp2);
          B_J3_VtxProb->push_back(B_Prob_tmp3);
          B_J4_VtxProb->push_back(B_Prob_tmp4);
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
          //cout<<"Tracker Isolation work"<<iMuon1->trackIso()<<endl;
          //cout<<"Ecal Isolation work"<<iMuon1->ecalIso()<<endl;
          //cout<<"Hcal Isolation work"<<iMuon1->hcalIso()<<endl;
          //cout<<"Calo Isolation work"<<iMuon1->caloIso()<<endl;
          B_Mu2_px->push_back(iMuon2->px());
          B_Mu2_py->push_back(iMuon2->py());
          B_Mu2_pz->push_back(iMuon2->pz());
          B_Mu2_pt->push_back(iMuon2->pt());
          B_Mu2_eta->push_back(iMuon2->eta());
          B_Mu2_phi->push_back(iMuon2->phi());
          B_Mu2_charge->push_back(iMuon2->charge());
          B_Mu2_soft->push_back(iMuon2->isSoftMuon(bestVtx));
          B_Mu2_tight->push_back(iMuon2->isTightMuon(bestVtx));
          B_Mu2_loose->push_back(muon::isLooseMuon(*iMuon2));
          B_Mu2_IsoTrack->push_back(iMuon2->trackIso());
          B_Mu2_IsoEcal->push_back(iMuon2->ecalIso());
          B_Mu2_IsoHcal->push_back(iMuon2->hcalIso());
          B_Mu2_IsoCalo->push_back(iMuon2->caloIso());

          B_Mu2_PaperIsoTrackRF04->push_back(
              (iMuon2->pfIsolationR04().sumChargedHadronPt +
               std::max(
                   0., iMuon2->pfIsolationR04().sumNeutralHadronEt + iMuon2->pfIsolationR04().sumPhotonEt - iMuon2->pfIsolationR04().sumPUPt * 0.5)) /
              iMuon2->pt());
          B_Mu2_PaperIsoTrackRF03->push_back(
              (iMuon2->pfIsolationR03().sumChargedHadronPt +
               std::max(
                   0., iMuon2->pfIsolationR03().sumNeutralHadronEt + iMuon2->pfIsolationR03().sumPhotonEt - iMuon2->pfIsolationR03().sumPUPt * 0.5)) /
              iMuon2->pt());

          B_Mu2_Paper3DIP->push_back(iMuon2->dB(pat::Muon::PV3D) / iMuon2->edB(pat::Muon::PV3D));

          B_Mu3_px->push_back(iMuon3->px());
          B_Mu3_py->push_back(iMuon3->py());
          B_Mu3_pz->push_back(iMuon3->pz());
          B_Mu3_pt->push_back(iMuon3->pt());
          B_Mu3_eta->push_back(iMuon3->eta());
          B_Mu3_phi->push_back(iMuon3->phi());
          B_Mu3_soft->push_back(iMuon3->isSoftMuon(bestVtx));
          B_Mu3_tight->push_back(iMuon3->isTightMuon(bestVtx));
          B_Mu3_loose->push_back(muon::isLooseMuon(*iMuon3));
          B_Mu3_charge->push_back(iMuon3->charge());
          B_Mu3_IsoTrack->push_back(iMuon3->trackIso());
          B_Mu3_IsoEcal->push_back(iMuon3->ecalIso());
          B_Mu3_IsoHcal->push_back(iMuon3->hcalIso());
          B_Mu3_IsoCalo->push_back(iMuon3->caloIso());

          B_Mu3_PaperIsoTrackRF04->push_back(
              (iMuon3->pfIsolationR04().sumChargedHadronPt +
               std::max(
                   0., iMuon3->pfIsolationR04().sumNeutralHadronEt + iMuon3->pfIsolationR04().sumPhotonEt - iMuon3->pfIsolationR04().sumPUPt * 0.5)) /
              iMuon3->pt());
          B_Mu3_PaperIsoTrackRF03->push_back(
              (iMuon3->pfIsolationR03().sumChargedHadronPt +
               std::max(
                   0., iMuon3->pfIsolationR03().sumNeutralHadronEt + iMuon3->pfIsolationR03().sumPhotonEt - iMuon3->pfIsolationR03().sumPUPt * 0.5)) /
              iMuon3->pt());

          B_Mu3_Paper3DIP->push_back(iMuon3->dB(pat::Muon::PV3D) / iMuon3->edB(pat::Muon::PV3D));

          B_Mu4_px->push_back(iMuon4->px());
          B_Mu4_py->push_back(iMuon4->py());
          B_Mu4_pz->push_back(iMuon4->pz());
          B_Mu4_pt->push_back(iMuon4->pt());
          B_Mu4_eta->push_back(iMuon4->eta());
          B_Mu4_phi->push_back(iMuon4->phi());
          B_Mu4_soft->push_back(iMuon4->isSoftMuon(bestVtx));
          B_Mu4_tight->push_back(iMuon4->isTightMuon(bestVtx));
          B_Mu4_loose->push_back(muon::isLooseMuon(*iMuon4));
          B_Mu4_charge->push_back(iMuon4->charge());
          B_Mu4_IsoTrack->push_back(iMuon4->trackIso());
          B_Mu4_IsoEcal->push_back(iMuon4->ecalIso());
          B_Mu4_IsoHcal->push_back(iMuon4->hcalIso());
          B_Mu4_IsoCalo->push_back(iMuon4->caloIso());

          B_Mu4_PaperIsoTrackRF04->push_back(
              (iMuon4->pfIsolationR04().sumChargedHadronPt +
               std::max(
                   0., iMuon4->pfIsolationR04().sumNeutralHadronEt + iMuon4->pfIsolationR04().sumPhotonEt - iMuon4->pfIsolationR04().sumPUPt * 0.5)) /
              iMuon4->pt());
          B_Mu4_PaperIsoTrackRF03->push_back(
              (iMuon4->pfIsolationR03().sumChargedHadronPt +
               std::max(
                   0., iMuon4->pfIsolationR03().sumNeutralHadronEt + iMuon4->pfIsolationR03().sumPhotonEt - iMuon4->pfIsolationR03().sumPUPt * 0.5)) /
              iMuon4->pt());

          B_Mu4_Paper3DIP->push_back(iMuon4->dB(pat::Muon::PV3D) / iMuon4->edB(pat::Muon::PV3D));

          B_J_xyP1->push_back(glbTrackP1->dxy(bestVtx.position()));
          B_J_xyM1->push_back(glbTrackM1->dxy(bestVtx.position()));
          B_J_zP1->push_back(glbTrackM1->dz(bestVtx.position()));
          B_J_zM1->push_back(glbTrackP1->dz(bestVtx.position()));

          B_J_xyP2->push_back(glbTrackP2->dxy(bestVtx.position()));
          B_J_xyM2->push_back(glbTrackM2->dxy(bestVtx.position()));
          B_J_zP2->push_back(glbTrackM2->dz(bestVtx.position()));
          B_J_zM2->push_back(glbTrackP2->dz(bestVtx.position()));

          //cout<<"End of all loop"<<endl;

          mu1mC2->push_back(glbTrackM1->normalizedChi2());
          //mumAngT->push_back( muon::isGoodMuon(*iMuon1,muon::TMLastStationAngTight) ); //
          mu1mNHits->push_back(glbTrackM1->numberOfValidHits());
          mu1mNPHits->push_back(glbTrackM1->hitPattern().numberOfValidPixelHits());
          mu1pC2->push_back(glbTrackP1->normalizedChi2());
          //mupAngT->push_back( muon::isGoodMuon(*iMuon2,muon::TMLastStationAngTight) );  //
          mu1pNHits->push_back(glbTrackP1->numberOfValidHits());
          mu1pNPHits->push_back(glbTrackP1->hitPattern().numberOfValidPixelHits());

          mu2mC2->push_back(glbTrackM2->normalizedChi2());
          //mumAngT->push_back( muon::isGoodMuon(*iMuon1,muon::TMLastStationAngTight) ); //
          mu2mNHits->push_back(glbTrackM2->numberOfValidHits());
          mu2mNPHits->push_back(glbTrackM2->hitPattern().numberOfValidPixelHits());
          mu2pC2->push_back(glbTrackP2->normalizedChi2());
          //mupAngT->push_back( muon::isGoodMuon(*iMuon2,muon::TMLastStationAngTight) );  //
          mu2pNHits->push_back(glbTrackP2->numberOfValidHits());
          mu2pNPHits->push_back(glbTrackP2->hitPattern().numberOfValidPixelHits());
          //muon pt eta phi arranged acc to P and M as M1,M2,M3,M4
          //M1_pt->push_back(M1.Pt());
          B_M1_pt->push_back(M1.Pt());
          B_M1_eta->push_back(M1.Eta());
          B_M1_phi->push_back(M1.Phi());
          B_M1_px->push_back(M1.Px());
          B_M1_py->push_back(M1.Py());
          B_M1_pz->push_back(M1.Pz());
          B_M2_pt->push_back(M2.Pt());
          B_M2_eta->push_back(M2.Eta());
          B_M2_phi->push_back(M2.Phi());
          B_M2_px->push_back(M2.Px());
          B_M2_py->push_back(M2.Py());
          B_M2_pz->push_back(M2.Pz());
          B_M3_pt->push_back(M3.Pt());
          B_M3_eta->push_back(M3.Eta());
          B_M3_phi->push_back(M3.Phi());
          B_M3_px->push_back(M3.Px());
          B_M3_py->push_back(M3.Py());
          B_M3_pz->push_back(M3.Pz());
          B_M4_pt->push_back(M4.Pt());
          B_M4_eta->push_back(M4.Eta());
          B_M4_phi->push_back(M4.Phi());
          B_M4_px->push_back(M4.Px());
          B_M4_py->push_back(M4.Py());
          B_M4_pz->push_back(M4.Pz());
          nB++;
          if (KinFit == true) {
            //muonParticles.clear();
          }
        }
      }
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
  B_U_TriggerPath->clear();
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
  B_J2_mass->clear();
  B_J2_px->clear();
  B_J2_py->clear();
  B_J2_pz->clear();
  B_J2_pt->clear();
  B_J2_eta->clear();
  B_J2_phi->clear();
  B_J2_rapidity->clear();
  B_J3_mass->clear();
  B_J3_px->clear();
  B_J3_py->clear();
  B_J3_pz->clear();
  B_J3_pt->clear();
  B_J3_eta->clear();
  B_J3_phi->clear();
  B_J3_rapidity->clear();
  B_J4_mass->clear();
  B_J4_px->clear();
  B_J4_py->clear();
  B_J4_pz->clear();
  B_J4_pt->clear();
  B_J4_eta->clear();
  B_J4_phi->clear();
  B_J4_rapidity->clear();

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
  B_J2_VtxPx->clear();
  B_J2_VtxPy->clear();
  B_J2_VtxPz->clear();
  B_J2_VtxPt->clear();
  B_J2_VtxEta->clear();
  B_J2_VtxPhi->clear();
  B_J2_VtxRapidity->clear();
  B_J2_VtxMass->clear();
  B_J2_PVx->clear();
  B_J2_PVy->clear();
  B_J2_PVz->clear();
  B_J2_PVxError->clear();
  B_J2_PVyError->clear();
  B_J2_PVzError->clear();
  B_J3_VtxPx->clear();
  B_J3_VtxPy->clear();
  B_J3_VtxPz->clear();
  B_J3_VtxPt->clear();
  B_J3_VtxEta->clear();
  B_J3_VtxPhi->clear();
  B_J3_VtxRapidity->clear();
  B_J3_VtxMass->clear();
  B_J3_PVx->clear();
  B_J3_PVy->clear();
  B_J3_PVz->clear();
  B_J3_PVxError->clear();
  B_J3_PVyError->clear();
  B_J3_PVzError->clear();
  B_J4_VtxPx->clear();
  B_J4_VtxPy->clear();
  B_J4_VtxPz->clear();
  B_J4_VtxPt->clear();
  B_J4_VtxEta->clear();
  B_J4_VtxPhi->clear();
  B_J4_VtxRapidity->clear();
  B_J4_VtxMass->clear();
  B_J4_PVx->clear();
  B_J4_PVy->clear();
  B_J4_PVz->clear();
  B_J4_PVxError->clear();
  B_J4_PVyError->clear();
  B_J4_PVzError->clear();

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
  B_Mu2_px->clear();
  B_Mu2_py->clear();
  B_Mu2_pz->clear();
  B_Mu2_charge->clear();
  B_Mu2_pt->clear();
  B_Mu2_eta->clear();
  B_Mu2_phi->clear();
  B_J1_VtxProb->clear();
  B_J2_VtxProb->clear();
  B_J3_VtxProb->clear();
  B_J4_VtxProb->clear();
  B_Mu2_soft->clear();
  B_Mu2_tight->clear();
  B_Mu2_loose->clear();
  B_Mu2_IsoTrack->clear();
  B_Mu2_IsoHcal->clear();
  B_Mu2_IsoEcal->clear();
  B_Mu2_IsoCalo->clear();
  B_Mu3_px->clear();
  B_Mu3_py->clear();
  B_Mu3_pz->clear();
  B_Mu3_pt->clear();
  B_Mu3_eta->clear();
  B_Mu3_phi->clear();
  B_Mu3_soft->clear();
  B_Mu3_tight->clear();
  B_Mu3_loose->clear();
  B_Mu3_charge->clear();
  B_Mu3_IsoTrack->clear();
  B_Mu3_IsoHcal->clear();
  B_Mu3_IsoEcal->clear();
  B_Mu3_IsoCalo->clear();
  B_Mu4_px->clear();
  B_Mu4_py->clear();
  B_Mu4_pz->clear();
  B_Mu4_pt->clear();
  B_Mu4_eta->clear();
  B_Mu4_phi->clear();
  B_Mu4_soft->clear();
  B_Mu4_tight->clear();
  B_Mu4_loose->clear();
  B_Mu4_charge->clear();
  B_Mu4_IsoTrack->clear();
  B_Mu4_IsoHcal->clear();
  B_Mu4_IsoEcal->clear();
  B_Mu4_IsoCalo->clear();
  B_J_xyP1->clear();
  B_J_xyM1->clear();
  B_J_zP1->clear();
  B_J_zM1->clear();
  B_J_xyP2->clear();
  B_J_xyM2->clear();
  B_J_zP2->clear();
  B_J_zM2->clear();
  B_Mu1_PaperIsoTrackRF04->clear();
  B_Mu1_PaperIsoTrackRF03->clear();
  B_Mu2_PaperIsoTrackRF04->clear();
  B_Mu2_PaperIsoTrackRF03->clear();
  B_Mu3_PaperIsoTrackRF04->clear();
  B_Mu3_PaperIsoTrackRF03->clear();
  B_Mu4_PaperIsoTrackRF04->clear();
  B_Mu4_PaperIsoTrackRF03->clear();

  B_Mu1_Paper3DIP->clear();
  B_Mu2_Paper3DIP->clear();
  B_Mu3_Paper3DIP->clear();
  B_Mu4_Paper3DIP->clear();

  mu1mC2->clear();
  mu1mNHits->clear();
  mu1mNPHits->clear();
  mu1pC2->clear();
  mu1pNHits->clear();
  mu1pNPHits->clear();
  mu2mC2->clear();
  mu2mNHits->clear();
  mu2mNPHits->clear();
  mu2pC2->clear();
  mu2pNHits->clear();
  mu2pNPHits->clear();

  B_M1_pt->clear();
  B_M1_eta->clear();
  B_M1_phi->clear();
  B_M1_px->clear();
  B_M1_px->clear();
  B_M1_pz->clear();
  B_M2_pt->clear();
  B_M2_eta->clear();
  B_M2_phi->clear();
  B_M2_px->clear();
  B_M2_px->clear();
  B_M2_pz->clear();
  B_M3_pt->clear();
  B_M3_eta->clear();
  B_M3_phi->clear();
  B_M3_px->clear();
  B_M3_px->clear();
  B_M3_pz->clear();
  B_M4_pt->clear();
  B_M4_eta->clear();
  B_M4_phi->clear();
  B_M4_px->clear();
  B_M4_px->clear();
  B_M4_pz->clear();

  //Gen muon
  B_J_GenMuonPt->clear();
  B_J_GenMuonEta->clear();
  B_J_GenMuonPhi->clear();
  B_Z_GenMuonPt->clear();
  B_Z_GenMuonEta->clear();
  B_Z_GenMuonPhi->clear();
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

  tree_->Branch("B_U_TriggerPath", &B_U_TriggerPath);
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

  tree_->Branch("B_J2_mass", &B_J2_mass);
  tree_->Branch("B_J2_px", &B_J2_px);
  tree_->Branch("B_J2_py", &B_J2_py);
  tree_->Branch("B_J2_pz", &B_J2_pz);
  tree_->Branch("B_J2_pt", &B_J2_pt);
  tree_->Branch("B_J2_eta", &B_J2_eta);
  tree_->Branch("B_J2_phi", &B_J2_phi);
  tree_->Branch("B_J2_rapidity", &B_J2_rapidity);

  tree_->Branch("B_J3_mass", &B_J3_mass);
  tree_->Branch("B_J3_px", &B_J3_px);
  tree_->Branch("B_J3_py", &B_J3_py);
  tree_->Branch("B_J3_pz", &B_J3_pz);
  tree_->Branch("B_J3_pt", &B_J3_pt);
  tree_->Branch("B_J3_eta", &B_J3_eta);
  tree_->Branch("B_J3_phi", &B_J3_phi);
  tree_->Branch("B_J3_rapidity", &B_J3_rapidity);

  tree_->Branch("B_J4_mass", &B_J4_mass);
  tree_->Branch("B_J4_px", &B_J4_px);
  tree_->Branch("B_J4_py", &B_J4_py);
  tree_->Branch("B_J4_pz", &B_J4_pz);
  tree_->Branch("B_J4_pt", &B_J4_pt);
  tree_->Branch("B_J4_eta", &B_J4_eta);
  tree_->Branch("B_J4_phi", &B_J4_phi);
  tree_->Branch("B_J4_rapidity", &B_J4_rapidity);

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
  tree_->Branch("B_J2_VtxPx", &B_J2_VtxPx);
  tree_->Branch("B_J2_VtxPy", &B_J2_VtxPy);
  tree_->Branch("B_J2_VtxPz", &B_J2_VtxPz);
  tree_->Branch("B_J2_VtxPt", &B_J2_VtxPt);
  tree_->Branch("B_J2_VtxEta", &B_J2_VtxEta);
  tree_->Branch("B_J2_VtxPhi", &B_J2_VtxPhi);
  tree_->Branch("B_J2_VtxRapidity", &B_J2_VtxRapidity);
  tree_->Branch("B_J2_VtxMass", &B_J2_VtxMass);
  tree_->Branch("B_J2_PVx", &B_J2_PVx);
  tree_->Branch("B_J2_PVy", &B_J2_PVy);
  tree_->Branch("B_J2_PVz", &B_J2_PVz);
  tree_->Branch("B_J2_PVxError", &B_J2_PVxError);
  tree_->Branch("B_J2_PVyError", &B_J2_PVyError);
  tree_->Branch("B_J2_PVzError", &B_J2_PVzError);
  tree_->Branch("B_J3_VtxPx", &B_J3_VtxPx);
  tree_->Branch("B_J3_VtxPy", &B_J3_VtxPy);
  tree_->Branch("B_J3_VtxPz", &B_J3_VtxPz);
  tree_->Branch("B_J3_VtxPt", &B_J3_VtxPt);
  tree_->Branch("B_J3_VtxEta", &B_J3_VtxEta);
  tree_->Branch("B_J3_VtxPhi", &B_J3_VtxPhi);
  tree_->Branch("B_J3_VtxRapidity", &B_J3_VtxRapidity);
  tree_->Branch("B_J3_VtxMass", &B_J3_VtxMass);
  tree_->Branch("B_J3_PVx", &B_J3_PVx);
  tree_->Branch("B_J3_PVy", &B_J3_PVy);
  tree_->Branch("B_J3_PVz", &B_J3_PVz);
  tree_->Branch("B_J3_PVxError", &B_J3_PVxError);
  tree_->Branch("B_J3_PVyError", &B_J3_PVyError);
  tree_->Branch("B_J3_PVzError", &B_J3_PVzError);

  tree_->Branch("B_J4_VtxPx", &B_J4_VtxPx);
  tree_->Branch("B_J4_VtxPy", &B_J4_VtxPy);
  tree_->Branch("B_J4_VtxPz", &B_J4_VtxPz);
  tree_->Branch("B_J4_VtxPt", &B_J4_VtxPt);
  tree_->Branch("B_J4_VtxEta", &B_J4_VtxEta);
  tree_->Branch("B_J4_VtxPhi", &B_J4_VtxPhi);
  tree_->Branch("B_J4_VtxRapidity", &B_J4_VtxRapidity);
  tree_->Branch("B_J4_VtxMass", &B_J4_VtxMass);
  tree_->Branch("B_J4_PVx", &B_J4_PVx);
  tree_->Branch("B_J4_PVy", &B_J4_PVy);
  tree_->Branch("B_J4_PVz", &B_J4_PVz);
  tree_->Branch("B_J4_PVxError", &B_J4_PVxError);
  tree_->Branch("B_J4_PVyError", &B_J4_PVyError);
  tree_->Branch("B_J4_PVzError", &B_J4_PVzError);

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

  tree_->Branch("B_Mu2_px", &B_Mu2_px);
  tree_->Branch("B_Mu2_py", &B_Mu2_py);
  tree_->Branch("B_Mu2_pz", &B_Mu2_pz);
  tree_->Branch("B_Mu2_pt", &B_Mu2_pt);
  tree_->Branch("B_Mu2_eta", &B_Mu2_eta);
  tree_->Branch("B_Mu2_phi", &B_Mu2_phi);
  tree_->Branch("B_Mu2_charge", &B_Mu2_charge);
  tree_->Branch("B_Mu2_soft", &B_Mu2_soft);
  tree_->Branch("B_Mu2_tight", &B_Mu2_tight);
  tree_->Branch("B_Mu2_loose", &B_Mu2_loose);
  tree_->Branch("B_Mu2_IsoTrack", &B_Mu2_IsoTrack);
  tree_->Branch("B_Mu2_IsoHcal", &B_Mu2_IsoHcal);
  tree_->Branch("B_Mu2_IsoEcal", &B_Mu2_IsoEcal);
  tree_->Branch("B_Mu2_IsoCalo", &B_Mu2_IsoCalo);

  tree_->Branch("B_Mu2_PaperIsoTrackRF03", &B_Mu2_PaperIsoTrackRF03);
  tree_->Branch("B_Mu2_PaperIsoTrackRF04", &B_Mu2_PaperIsoTrackRF04);

  tree_->Branch("B_Mu2_Paper3DIP", &B_Mu2_Paper3DIP);

  tree_->Branch("B_Mu3_px", &B_Mu3_px);
  tree_->Branch("B_Mu3_py", &B_Mu3_py);
  tree_->Branch("B_Mu3_pz", &B_Mu3_pz);
  tree_->Branch("B_Mu3_pt", &B_Mu3_pt);
  tree_->Branch("B_Mu3_eta", &B_Mu3_eta);
  tree_->Branch("B_Mu3_phi", &B_Mu3_phi);
  tree_->Branch("B_Mu3_soft", &B_Mu3_soft);
  tree_->Branch("B_Mu3_tight", &B_Mu3_tight);
  tree_->Branch("B_Mu3_loose", &B_Mu3_loose);
  tree_->Branch("B_Mu3_charge", &B_Mu3_charge);
  tree_->Branch("B_Mu3_IsoTrack", &B_Mu3_IsoTrack);
  tree_->Branch("B_Mu3_IsoHcal", &B_Mu3_IsoHcal);
  tree_->Branch("B_Mu3_IsoEcal", &B_Mu3_IsoEcal);
  tree_->Branch("B_Mu3_IsoCalo", &B_Mu3_IsoCalo);

  tree_->Branch("B_Mu3_PaperIsoTrackRF03", &B_Mu3_PaperIsoTrackRF03);
  tree_->Branch("B_Mu3_PaperIsoTrackRF04", &B_Mu3_PaperIsoTrackRF04);

  tree_->Branch("B_Mu3_Paper3DIP", &B_Mu3_Paper3DIP);

  tree_->Branch("B_Mu4_px", &B_Mu4_px);
  tree_->Branch("B_Mu4_py", &B_Mu4_py);
  tree_->Branch("B_Mu4_pz", &B_Mu4_pz);
  tree_->Branch("B_Mu4_pt", &B_Mu4_pt);
  tree_->Branch("B_Mu4_eta", &B_Mu4_eta);
  tree_->Branch("B_Mu4_phi", &B_Mu4_phi);
  tree_->Branch("B_Mu4_soft", &B_Mu4_soft);
  tree_->Branch("B_Mu4_tight", &B_Mu4_tight);
  tree_->Branch("B_Mu4_loose", &B_Mu4_loose);
  tree_->Branch("B_Mu4_charge", &B_Mu4_charge);
  tree_->Branch("B_Mu4_IsoTrack", &B_Mu4_IsoTrack);
  tree_->Branch("B_Mu4_IsoHcal", &B_Mu4_IsoHcal);
  tree_->Branch("B_Mu4_IsoEcal", &B_Mu4_IsoEcal);
  tree_->Branch("B_Mu4_IsoCalo", &B_Mu4_IsoCalo);

  tree_->Branch("B_Mu4_PaperIsoTrackRF03", &B_Mu4_PaperIsoTrackRF03);
  tree_->Branch("B_Mu4_PaperIsoTrackRF04", &B_Mu4_PaperIsoTrackRF04);

  tree_->Branch("B_Mu4_Paper3DIP", &B_Mu4_Paper3DIP);

  tree_->Branch("B_J1_VtxProb", &B_J1_VtxProb);
  tree_->Branch("B_J2_VtxProb", &B_J2_VtxProb);
  tree_->Branch("B_J3_VtxProb", &B_J3_VtxProb);
  tree_->Branch("B_J4_VtxProb", &B_J4_VtxProb);

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

  tree_->Branch("mu2mC2", &mu2mC2);
  tree_->Branch("mu2mNHits", &mu2mNHits);
  tree_->Branch("mu2mNPHits", &mu2mNPHits);
  tree_->Branch("mu2pC2", &mu2pC2);
  tree_->Branch("mu2pNHits", &mu2pNHits);
  tree_->Branch("mu2pNPHits", &mu2pNPHits);

  tree_->Branch("B_M1_pt", &B_M1_pt);
  tree_->Branch("B_M1_eta", &B_M1_eta);
  tree_->Branch("B_M1_phi", &B_M1_phi);
  tree_->Branch("B_M1_px", &B_M1_px);
  tree_->Branch("B_M1_py", &B_M1_py);
  tree_->Branch("B_M1_pz", &B_M1_pz);

  tree_->Branch("B_M2_pt", &B_M2_pt);
  tree_->Branch("B_M2_eta", &B_M2_eta);
  tree_->Branch("B_M2_phi", &B_M2_phi);
  tree_->Branch("B_M2_px", &B_M2_px);
  tree_->Branch("B_M2_py", &B_M2_py);
  tree_->Branch("B_M2_pz", &B_M2_pz);

  tree_->Branch("B_M3_pt", &B_M3_pt);
  tree_->Branch("B_M3_eta", &B_M3_eta);
  tree_->Branch("B_M3_phi", &B_M3_phi);
  tree_->Branch("B_M3_px", &B_M3_px);
  tree_->Branch("B_M3_py", &B_M3_py);
  tree_->Branch("B_M3_pz", &B_M3_pz);

  tree_->Branch("B_M4_pt", &B_M4_pt);
  tree_->Branch("B_M4_eta", &B_M4_eta);
  tree_->Branch("B_M4_phi", &B_M4_phi);
  tree_->Branch("B_M4_px", &B_M4_px);
  tree_->Branch("B_M4_py", &B_M4_py);
  tree_->Branch("B_M4_pz", &B_M4_pz);

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
