// -*- C++ -*-
// miniAODeemm: Run-2 H -> Z(e e) J/psi(mu mu) MiniAOD preselection ntuplizer.
// Design target: open-ended ML ntuple. No m_ee_mumu cut and no candidate-score pruning.

#include "ZeeJmmAnalyzer/miniAODeemm/plugins/miniAODeemm.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaR2.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"

namespace {
constexpr float kElectronMass = 0.00051099895f;
constexpr float kMuonMass = 0.1056583745f;
constexpr float kBad = -999.0f;

// Hardcoded Run-2 electron IDs embedded in MiniAOD.
const std::string kLooseElectronID = "mvaEleID-Fall17-iso-V2-wpLoose";
const std::string kElectronWP90ID = "mvaEleID-Fall17-iso-V2-wp90";
const std::string kElectronWP80ID = "mvaEleID-Fall17-iso-V2-wp80";
const std::vector<std::string> kElectronMVARawNames = {
    "ElectronMVAEstimatorRun2Fall17IsoV2Values",
    "ElectronMVAEstimatorRun2Fall17NoIsoV2Values",
    "ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"};


float safeRatio(float numerator, float denominator) {
  if (std::abs(denominator) < 1e-6f) return kBad;
  return numerator / denominator;
}

bool pathMatches(const std::string& pathName, const std::string& pattern) {
  return pathName.find(pattern) != std::string::npos;
}

bool firedElectronTrigger(const edm::TriggerResults& triggerBits,
                          const edm::TriggerNames& triggerNames,
                          const std::string& triggerPattern) {
  if (triggerPattern.empty()) return false;
  for (unsigned int i = 0; i < triggerBits.size(); ++i) {
    if (!triggerBits.accept(i)) continue;
    const std::string& name = triggerNames.triggerName(i);
    if (pathMatches(name, triggerPattern)) return true;
  }
  return false;
}

float vertexProbability(const TransientVertex& vtx) {
  if (!vtx.isValid()) return kBad;
  return ChiSquaredProbability(vtx.totalChiSquared(), static_cast<int>(std::round(vtx.degreesOfFreedom())));
}

float pfAbsIso03(const pat::Muon& mu) {
  const auto& iso = mu.pfIsolationR03();
  return iso.sumChargedHadronPt + std::max(0.0f, static_cast<float>(iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5 * iso.sumPUPt));
}

float pfRelIso03(const pat::Muon& mu) {
  return safeRatio(pfAbsIso03(mu), mu.pt());
}

float trackChi2(const reco::TrackRef& trk) {
  if (trk.isNull()) return kBad;
  return trk->normalizedChi2();
}

int nValidHits(const reco::TrackRef& trk) {
  if (trk.isNull()) return -999;
  return trk->numberOfValidHits();
}

int nValidPixelHits(const reco::TrackRef& trk) {
  if (trk.isNull()) return -999;
  return trk->hitPattern().numberOfValidPixelHits();
}

float dxy(const reco::TrackRef& trk, const reco::Vertex& pv) {
  if (trk.isNull()) return kBad;
  return trk->dxy(pv.position());
}

float dz(const reco::TrackRef& trk, const reco::Vertex& pv) {
  if (trk.isNull()) return kBad;
  return trk->dz(pv.position());
}

float dxy(const reco::GsfTrackRef& trk, const reco::Vertex& pv) {
  if (trk.isNull()) return kBad;
  return trk->dxy(pv.position());
}

float dz(const reco::GsfTrackRef& trk, const reco::Vertex& pv) {
  if (trk.isNull()) return kBad;
  return trk->dz(pv.position());
}

float electronID(const pat::Electron& ele, const std::string& idName) {
  if (idName.empty()) return 0.0f;
  try {
    return ele.electronID(idName);
  } catch (...) {
    return 0.0f;
  }
}

float electronMVARaw(const pat::Electron& ele, const std::vector<std::string>& rawNames) {
  for (const auto& name : rawNames) {
    if (name.empty()) continue;
    if (ele.hasUserFloat(name)) return ele.userFloat(name);
  }
  return kBad;
}

bool hasElectronTriggerPath(const pat::TriggerObjectStandAlone& obj,
                            const std::string& triggerPattern) {
  if (triggerPattern.empty()) return false;
  for (const auto& pathName : obj.pathNames()) {
    if (pathMatches(pathName, triggerPattern)) return true;
  }
  return false;
}

std::vector<const pat::TriggerObjectStandAlone*> getMatchedObjs(
    float eta,
    float phi,
    const std::vector<pat::TriggerObjectStandAlone>& trigObjs,
    float drMax) {
  std::vector<const pat::TriggerObjectStandAlone*> matched;
  for (const auto& obj : trigObjs) {
    if (reco::deltaR(eta, phi, obj.eta(), obj.phi()) < drMax) {
      matched.push_back(&obj);
    }
  }
  return matched;
}

float minTriggerDR(const pat::Electron& ele,
                   const std::vector<pat::TriggerObjectStandAlone>& trigObjs,
                   const std::string& triggerPattern) {
  float best = 999.0f;
  for (const auto& obj : trigObjs) {
    if (!hasElectronTriggerPath(obj, triggerPattern)) continue;
    const float dr = reco::deltaR(ele.eta(), ele.phi(), obj.eta(), obj.phi());
    if (dr < best) best = dr;
  }
  return best;
}

float pairTrackIso(const edm::View<pat::PackedCandidate>& packed,
                   const TLorentzVector& pair,
                   const std::vector<TLorentzVector>& selectedLeptons,
                   float cone) {
  float sumPt = 0.0f;
  for (const auto& cand : packed) {
    if (cand.charge() == 0) continue;
    if (cand.pt() < 0.5f) continue;
    if (cand.fromPV() <= 1) continue;
    if (reco::deltaR(cand.eta(), cand.phi(), pair.Eta(), pair.Phi()) > cone) continue;

    bool overlapsSelectedLepton = false;
    for (const auto& lep : selectedLeptons) {
      if (reco::deltaR(cand.eta(), cand.phi(), lep.Eta(), lep.Phi()) < 1.0e-3) {
        overlapsSelectedLepton = true;
        break;
      }
    }
    if (!overlapsSelectedLepton) sumPt += cand.pt();
  }
  return sumPt;
}

float cosThetaLeptonInPairRest(const TLorentzVector& leptonPlus,
                               const TLorentzVector& pair,
                               const TLorentzVector& fourL) {
  if (pair.P() <= 0.0 || fourL.M() <= 0.0) return kBad;

  TLorentzVector lepRest = leptonPlus;
  lepRest.Boost(-pair.BoostVector());

  TLorentzVector pairInFourRest = pair;
  pairInFourRest.Boost(-fourL.BoostVector());

  const TVector3 axis = pairInFourRest.Vect().Unit();
  const TVector3 lepDir = lepRest.Vect().Unit();
  if (axis.Mag() == 0.0 || lepDir.Mag() == 0.0) return kBad;
  return lepDir.Dot(axis);
}

float decayPlaneAngle(const TLorentzVector& a1,
                      const TLorentzVector& a2,
                      const TLorentzVector& b1,
                      const TLorentzVector& b2,
                      const TLorentzVector& fourL) {
  TLorentzVector aa1 = a1;
  TLorentzVector aa2 = a2;
  TLorentzVector bb1 = b1;
  TLorentzVector bb2 = b2;
  const TVector3 boost = -fourL.BoostVector();
  aa1.Boost(boost);
  aa2.Boost(boost);
  bb1.Boost(boost);
  bb2.Boost(boost);

  TVector3 n1 = aa1.Vect().Cross(aa2.Vect());
  TVector3 n2 = bb1.Vect().Cross(bb2.Vect());
  if (n1.Mag() == 0.0 || n2.Mag() == 0.0) return kBad;
  n1 = n1.Unit();
  n2 = n2.Unit();
  const double dot = std::max(-1.0, std::min(1.0, n1.Dot(n2)));
  return std::acos(dot);
}

void pushP4(std::vector<float>& mass,
            std::vector<float>& px,
            std::vector<float>& py,
            std::vector<float>& pz,
            std::vector<float>& pt,
            std::vector<float>& eta,
            std::vector<float>& phi,
            std::vector<float>& rapidity,
            const TLorentzVector& p4) {
  mass.push_back(p4.M());
  px.push_back(p4.Px());
  py.push_back(p4.Py());
  pz.push_back(p4.Pz());
  pt.push_back(p4.Pt());
  eta.push_back(p4.Eta());
  phi.push_back(p4.Phi());
  rapidity.push_back(p4.Rapidity());
}

void pushObjP4(std::vector<float>& px,
               std::vector<float>& py,
               std::vector<float>& pz,
               std::vector<float>& pt,
               std::vector<float>& eta,
               std::vector<float>& phi,
               const TLorentzVector& p4) {
  px.push_back(p4.Px());
  py.push_back(p4.Py());
  pz.push_back(p4.Pz());
  pt.push_back(p4.Pt());
  eta.push_back(p4.Eta());
  phi.push_back(p4.Phi());
}

}  // namespace

miniAODeemm::miniAODeemm(const edm::ParameterSet& iConfig)
    : electronToken_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
      muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      packedCandToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedCandidates"))),
      primaryVerticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
      ElectronTriggerString_(iConfig.getParameter<std::string>("ElectronTrigger")),
      requireTrigger_(iConfig.getParameter<bool>("requireTrigger")),
      requireTriggerMatch_(iConfig.getParameter<bool>("requireTriggerMatch")),
      keepEmptyEvents_(iConfig.getParameter<bool>("keepEmptyEvents")),
      requireLooseElectronID_(iConfig.getParameter<bool>("requireLooseElectronID")),
      isMC_(iConfig.getParameter<bool>("isMC")) {}

void miniAODeemm::beginJob() {
  edm::Service<TFileService> fs;
  if (!fs) throw cms::Exception("TFileService") << "TFileService is required.";
  fs->file().cd();
  tree_ = new TTree("ntuple", "ntuple");
  tree_->SetDirectory(&fs->file());

  tree_->Branch("nB", &nB, "nB/i");

  branch("Run", Run); branch("LumiBlock", LumiBlock); branch("Event", Event); branch("nPV", nPV);
  branch("passEleTrigger", passEleTrigger); branch("passEleTriggerMatch", passEleTriggerMatch);

  branch("fourL_mass", fourL_mass); branch("fourL_px", fourL_px); branch("fourL_py", fourL_py); branch("fourL_pz", fourL_pz); branch("fourL_pt", fourL_pt); branch("fourL_eta", fourL_eta); branch("fourL_phi", fourL_phi); branch("fourL_rapidity", fourL_rapidity);
  branch("fourL_vtxProb", fourL_vtxProb); branch("fourL_pvx", fourL_pvx); branch("fourL_pvy", fourL_pvy); branch("fourL_pvz", fourL_pvz); branch("fourL_pvxError", fourL_pvxError); branch("fourL_pvyError", fourL_pvyError); branch("fourL_pvzError", fourL_pvzError);

  branch("Z_mass", Z_mass); branch("Z_px", Z_px); branch("Z_py", Z_py); branch("Z_pz", Z_pz); branch("Z_pt", Z_pt); branch("Z_eta", Z_eta); branch("Z_phi", Z_phi); branch("Z_rapidity", Z_rapidity); branch("Z_vtxProb", Z_vtxProb); branch("Z_dR_ee", Z_dR_ee);
  branch("Z_trackIso03", Z_trackIso03); branch("Z_trackIso04", Z_trackIso04); branch("Z_relIso03", Z_relIso03); branch("Z_relIso04", Z_relIso04);

  branch("Jpsi_mass", Jpsi_mass); branch("Jpsi_px", Jpsi_px); branch("Jpsi_py", Jpsi_py); branch("Jpsi_pz", Jpsi_pz); branch("Jpsi_pt", Jpsi_pt); branch("Jpsi_eta", Jpsi_eta); branch("Jpsi_phi", Jpsi_phi); branch("Jpsi_rapidity", Jpsi_rapidity); branch("Jpsi_vtxProb", Jpsi_vtxProb); branch("Jpsi_dR_mumu", Jpsi_dR_mumu);
  branch("Jpsi_trackIso03", Jpsi_trackIso03); branch("Jpsi_trackIso04", Jpsi_trackIso04); branch("Jpsi_relIso03", Jpsi_relIso03); branch("Jpsi_relIso04", Jpsi_relIso04);

  branch("Z_Jpsi_dR", Z_Jpsi_dR); branch("Z_Jpsi_dPhi", Z_Jpsi_dPhi); branch("Z_Jpsi_dEta", Z_Jpsi_dEta); branch("Z_Jpsi_dY", Z_Jpsi_dY); branch("pt_balance", pt_balance);
  branch("cosTheta_Z_ePlus", cosTheta_Z_ePlus); branch("cosTheta_Jpsi_muPlus", cosTheta_Jpsi_muPlus); branch("phi_decayPlane_Z_Jpsi", phi_decayPlane_Z_Jpsi);

  branch("e1_px", e1_px); branch("e1_py", e1_py); branch("e1_pz", e1_pz); branch("e1_pt", e1_pt); branch("e1_eta", e1_eta); branch("e1_phi", e1_phi); branch("e1_charge", e1_charge);
  branch("e2_px", e2_px); branch("e2_py", e2_py); branch("e2_pz", e2_pz); branch("e2_pt", e2_pt); branch("e2_eta", e2_eta); branch("e2_phi", e2_phi); branch("e2_charge", e2_charge);
  branch("e1_dxy", e1_dxy); branch("e1_dz", e1_dz); branch("e2_dxy", e2_dxy); branch("e2_dz", e2_dz);
  branch("e1_mvaRaw", e1_mvaRaw); branch("e2_mvaRaw", e2_mvaRaw);
  branch("e1_passLooseID", e1_passLooseID); branch("e2_passLooseID", e2_passLooseID); branch("e1_passWP90", e1_passWP90); branch("e2_passWP90", e2_passWP90); branch("e1_passWP80", e1_passWP80); branch("e2_passWP80", e2_passWP80);
  branch("e1_triggerMatched", e1_triggerMatched); branch("e2_triggerMatched", e2_triggerMatched); branch("e1_triggerDR", e1_triggerDR); branch("e2_triggerDR", e2_triggerDR);

  branch("mu1_px", mu1_px); branch("mu1_py", mu1_py); branch("mu1_pz", mu1_pz); branch("mu1_pt", mu1_pt); branch("mu1_eta", mu1_eta); branch("mu1_phi", mu1_phi); branch("mu1_charge", mu1_charge);
  branch("mu2_px", mu2_px); branch("mu2_py", mu2_py); branch("mu2_pz", mu2_pz); branch("mu2_pt", mu2_pt); branch("mu2_eta", mu2_eta); branch("mu2_phi", mu2_phi); branch("mu2_charge", mu2_charge);
  branch("mu1_soft", mu1_soft); branch("mu2_soft", mu2_soft); branch("mu1_loose", mu1_loose); branch("mu2_loose", mu2_loose); branch("mu1_tight", mu1_tight); branch("mu2_tight", mu2_tight);
  branch("mu1_pfRelIso03", mu1_pfRelIso03); branch("mu2_pfRelIso03", mu2_pfRelIso03); branch("mu1_pfAbsIso03", mu1_pfAbsIso03); branch("mu2_pfAbsIso03", mu2_pfAbsIso03);
  branch("mu1_dxy", mu1_dxy); branch("mu2_dxy", mu2_dxy); branch("mu1_dz", mu1_dz); branch("mu2_dz", mu2_dz); branch("mu1_dB3D", mu1_dB3D); branch("mu2_dB3D", mu2_dB3D);
  branch("mu1_normChi2", mu1_normChi2); branch("mu2_normChi2", mu2_normChi2); branch("mu1_nValidHits", mu1_nValidHits); branch("mu2_nValidHits", mu2_nValidHits); branch("mu1_nValidPixelHits", mu1_nValidPixelHits); branch("mu2_nValidPixelHits", mu2_nValidPixelHits);

  branch("nExtraLooseElectrons", nExtraLooseElectrons); branch("nExtraLooseMuons", nExtraLooseMuons);
}

void miniAODeemm::clearVectors() {
#define CLR(x) x.clear()
  nB = 0;
  CLR(Run); CLR(LumiBlock); CLR(Event); CLR(nPV); CLR(passEleTrigger); CLR(passEleTriggerMatch);
  CLR(fourL_mass); CLR(fourL_px); CLR(fourL_py); CLR(fourL_pz); CLR(fourL_pt); CLR(fourL_eta); CLR(fourL_phi); CLR(fourL_rapidity); CLR(fourL_vtxProb); CLR(fourL_pvx); CLR(fourL_pvy); CLR(fourL_pvz); CLR(fourL_pvxError); CLR(fourL_pvyError); CLR(fourL_pvzError);
  CLR(Z_mass); CLR(Z_px); CLR(Z_py); CLR(Z_pz); CLR(Z_pt); CLR(Z_eta); CLR(Z_phi); CLR(Z_rapidity); CLR(Z_vtxProb); CLR(Z_dR_ee); CLR(Z_trackIso03); CLR(Z_trackIso04); CLR(Z_relIso03); CLR(Z_relIso04);
  CLR(Jpsi_mass); CLR(Jpsi_px); CLR(Jpsi_py); CLR(Jpsi_pz); CLR(Jpsi_pt); CLR(Jpsi_eta); CLR(Jpsi_phi); CLR(Jpsi_rapidity); CLR(Jpsi_vtxProb); CLR(Jpsi_dR_mumu); CLR(Jpsi_trackIso03); CLR(Jpsi_trackIso04); CLR(Jpsi_relIso03); CLR(Jpsi_relIso04);
  CLR(Z_Jpsi_dR); CLR(Z_Jpsi_dPhi); CLR(Z_Jpsi_dEta); CLR(Z_Jpsi_dY); CLR(pt_balance); CLR(cosTheta_Z_ePlus); CLR(cosTheta_Jpsi_muPlus); CLR(phi_decayPlane_Z_Jpsi);
  CLR(e1_px); CLR(e1_py); CLR(e1_pz); CLR(e1_pt); CLR(e1_eta); CLR(e1_phi); CLR(e2_px); CLR(e2_py); CLR(e2_pz); CLR(e2_pt); CLR(e2_eta); CLR(e2_phi); CLR(e1_charge); CLR(e2_charge); CLR(e1_dxy); CLR(e1_dz); CLR(e2_dxy); CLR(e2_dz); CLR(e1_mvaRaw); CLR(e2_mvaRaw); CLR(e1_passLooseID); CLR(e2_passLooseID); CLR(e1_passWP90); CLR(e2_passWP90); CLR(e1_passWP80); CLR(e2_passWP80); CLR(e1_triggerMatched); CLR(e2_triggerMatched); CLR(e1_triggerDR); CLR(e2_triggerDR);
  CLR(mu1_px); CLR(mu1_py); CLR(mu1_pz); CLR(mu1_pt); CLR(mu1_eta); CLR(mu1_phi); CLR(mu2_px); CLR(mu2_py); CLR(mu2_pz); CLR(mu2_pt); CLR(mu2_eta); CLR(mu2_phi); CLR(mu1_charge); CLR(mu2_charge); CLR(mu1_soft); CLR(mu2_soft); CLR(mu1_loose); CLR(mu2_loose); CLR(mu1_tight); CLR(mu2_tight); CLR(mu1_pfRelIso03); CLR(mu2_pfRelIso03); CLR(mu1_pfAbsIso03); CLR(mu2_pfAbsIso03); CLR(mu1_dxy); CLR(mu2_dxy); CLR(mu1_dz); CLR(mu2_dz); CLR(mu1_dB3D); CLR(mu2_dB3D); CLR(mu1_normChi2); CLR(mu2_normChi2); CLR(mu1_nValidHits); CLR(mu2_nValidHits); CLR(mu1_nValidPixelHits); CLR(mu2_nValidPixelHits);
  CLR(nExtraLooseElectrons); CLR(nExtraLooseMuons);
#undef CLR
}

void miniAODeemm::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  clearVectors();

  edm::Handle<edm::View<pat::Electron>> electronHandle;
  iEvent.getByToken(electronToken_, electronHandle);
  if (!electronHandle.isValid()) return;

  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);
  if (!muonHandle.isValid()) return;

  edm::Handle<edm::View<pat::PackedCandidate>> packedCandHandle;
  iEvent.getByToken(packedCandToken_, packedCandHandle);
  if (!packedCandHandle.isValid()) return;

  edm::Handle<reco::VertexCollection> primaryVerticesHandle;
  iEvent.getByToken(primaryVerticesToken_, primaryVerticesHandle);
  if (!primaryVerticesHandle.isValid() || primaryVerticesHandle->empty()) return;
  const reco::Vertex& bestVtx = primaryVerticesHandle->front();

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  if (!triggerBits.isValid()) return;
  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerBits);
  bool firedEleTrigger = firedElectronTrigger(*triggerBits, triggerNames, ElectronTriggerString_);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjectsToken_, triggerObjects);
  std::vector<pat::TriggerObjectStandAlone> unpackedTriggerObjects;
  if (triggerObjects.isValid()) {
    unpackedTriggerObjects.reserve(triggerObjects->size());
    for (auto obj : *triggerObjects) {
      obj.unpackPathNames(triggerNames);
      obj.unpackFilterLabels(iEvent, *triggerBits);
      unpackedTriggerObjects.push_back(obj);
    }
  }

  // Special treatment for 2017 when using HLT_Ele32_WPTight_Gsf_L1DoubleEG_v.
  // This emulates HLT_Ele32_WPTight_Gsf, which was not in the 2017 menu, by
  // additionally requiring the matched trigger object to pass hltEGL1SingleEGOrFilter.
  if (firedEleTrigger &&
      ElectronTriggerString_.find("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v") != std::string::npos) {
    bool passedL1seed = false;
    for (const auto& ele : *electronHandle) {
      if (!ele.superCluster().isNonnull()) continue;
      const float eta = ele.superCluster()->eta();
      const float phi = ele.superCluster()->phi();
      const auto matchedTrigObjs = getMatchedObjs(eta, phi, unpackedTriggerObjects, 0.1f);
      for (const auto* trigObj : matchedTrigObjs) {
        if (trigObj->hasFilterLabel("hltEGL1SingleEGOrFilter")) {
          passedL1seed = true;
          break;
        }
      }
      if (passedL1seed) break;
    }
    firedEleTrigger = firedEleTrigger && passedL1seed;
  }

  if (requireTrigger_ && !firedEleTrigger) return;

  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);
  KalmanVertexFitter kvf(true);

  std::vector<size_t> looseElectronIdx;
  looseElectronIdx.reserve(electronHandle->size());
  for (size_t i = 0; i < electronHandle->size(); ++i) {
    const auto& e = electronHandle->at(i);
    if (e.pt() < 5.0) continue;
    if (std::abs(e.eta()) > 2.5) continue;
    if (e.gsfTrack().isNull()) continue;
    if (requireLooseElectronID_ && electronID(e, kLooseElectronID) <= 0.5f) continue;
    looseElectronIdx.push_back(i);
  }

  std::vector<size_t> looseMuonIdx;
  looseMuonIdx.reserve(muonHandle->size());
  for (size_t i = 0; i < muonHandle->size(); ++i) {
    const auto& mu = muonHandle->at(i);
    if (mu.pt() < 2.0) continue;
    if (std::abs(mu.eta()) > 2.4) continue;
    if (mu.innerTrack().isNull()) continue;
    if (!mu.innerTrack()->quality(reco::TrackBase::highPurity)) continue;
    looseMuonIdx.push_back(i);
  }

  for (size_t ie1 = 0; ie1 < electronHandle->size(); ++ie1) {
    for (size_t ie2 = ie1 + 1; ie2 < electronHandle->size(); ++ie2) {
      const auto& e1 = electronHandle->at(ie1);
      const auto& e2 = electronHandle->at(ie2);

      // Keep the two-electron quality and preselection together for readability.
      if (e1.gsfTrack().isNull() || e2.gsfTrack().isNull()) continue;
      if (e1.pt() < 5.0 || e2.pt() < 5.0) continue;
      if (std::abs(e1.eta()) > 2.5 || std::abs(e2.eta()) > 2.5) continue;
      if (requireLooseElectronID_ &&
          (electronID(e1, kLooseElectronID) <= 0.5f || electronID(e2, kLooseElectronID) <= 0.5f)) continue;
      if (e1.charge() * e2.charge() >= 0) continue;
      if (std::max(e1.pt(), e2.pt()) < 20.0) continue;

      const pat::Electron* ePlusPtr = (e1.charge() > 0) ? &e1 : &e2;
      const pat::Electron* eMinusPtr = (e1.charge() > 0) ? &e2 : &e1;
      const pat::Electron* eLeadPtr = (e1.pt() >= e2.pt()) ? &e1 : &e2;
      const pat::Electron* eSubPtr = (e1.pt() >= e2.pt()) ? &e2 : &e1;

      TLorentzVector EPlus, EMinus, ELead, ESub;
      EPlus.SetPtEtaPhiM(ePlusPtr->pt(), ePlusPtr->eta(), ePlusPtr->phi(), kElectronMass);
      EMinus.SetPtEtaPhiM(eMinusPtr->pt(), eMinusPtr->eta(), eMinusPtr->phi(), kElectronMass);
      ELead.SetPtEtaPhiM(eLeadPtr->pt(), eLeadPtr->eta(), eLeadPtr->phi(), kElectronMass);
      ESub.SetPtEtaPhiM(eSubPtr->pt(), eSubPtr->eta(), eSubPtr->phi(), kElectronMass);
      const TLorentzVector Zee = EPlus + EMinus;
      if (!(Zee.M() > 60.0 && Zee.M() < 120.0)) continue;

      std::vector<reco::TransientTrack> zeeTracks{ttBuilder->build(ePlusPtr->gsfTrack()), ttBuilder->build(eMinusPtr->gsfTrack())};
      const TransientVertex zeeVtx = kvf.vertex(zeeTracks);
      const float zeeProb = vertexProbability(zeeVtx);
      if (!(zeeProb > 0.001)) continue;

      const float eLeadTrigDR = minTriggerDR(*eLeadPtr, unpackedTriggerObjects, ElectronTriggerString_);
      const float eSubTrigDR = minTriggerDR(*eSubPtr, unpackedTriggerObjects, ElectronTriggerString_);
      const bool eLeadTrigMatch = eLeadTrigDR < 0.3;
      const bool eSubTrigMatch = eSubTrigDR < 0.3;
      const bool hasEleTriggerMatch = eLeadTrigMatch || eSubTrigMatch;
      if (requireTriggerMatch_ && !hasEleTriggerMatch) continue;

      for (size_t im1 = 0; im1 < muonHandle->size(); ++im1) {
        for (size_t im2 = im1 + 1; im2 < muonHandle->size(); ++im2) {
          const auto& m1 = muonHandle->at(im1);
          const auto& m2 = muonHandle->at(im2);

          // Keep the two-muon quality and preselection together for readability.
          if (m1.innerTrack().isNull() || m2.innerTrack().isNull()) continue;
          if (!m1.innerTrack()->quality(reco::TrackBase::highPurity) ||
              !m2.innerTrack()->quality(reco::TrackBase::highPurity)) continue;
          if (m1.pt() < 2.0 || m2.pt() < 2.0) continue;
          if (std::abs(m1.eta()) > 2.4 || std::abs(m2.eta()) > 2.4) continue;
          if (m1.charge() * m2.charge() >= 0) continue;

          const pat::Muon* muPlusPtr = (m1.charge() > 0) ? &m1 : &m2;
          const pat::Muon* muMinusPtr = (m1.charge() > 0) ? &m2 : &m1;
          const pat::Muon* muLeadPtr = (m1.pt() >= m2.pt()) ? &m1 : &m2;
          const pat::Muon* muSubPtr = (m1.pt() >= m2.pt()) ? &m2 : &m1;

          TLorentzVector MUP, MUM, MuLead, MuSub;
          MUP.SetPtEtaPhiM(muPlusPtr->pt(), muPlusPtr->eta(), muPlusPtr->phi(), kMuonMass);
          MUM.SetPtEtaPhiM(muMinusPtr->pt(), muMinusPtr->eta(), muMinusPtr->phi(), kMuonMass);
          MuLead.SetPtEtaPhiM(muLeadPtr->pt(), muLeadPtr->eta(), muLeadPtr->phi(), kMuonMass);
          MuSub.SetPtEtaPhiM(muSubPtr->pt(), muSubPtr->eta(), muSubPtr->phi(), kMuonMass);
          const TLorentzVector Jpsi = MUP + MUM;
          if (!(Jpsi.M() > 2.8 && Jpsi.M() < 3.4)) continue;

          std::vector<reco::TransientTrack> jpsiTracks{ttBuilder->build(muPlusPtr->innerTrack()), ttBuilder->build(muMinusPtr->innerTrack())};
          const TransientVertex jpsiVtx = kvf.vertex(jpsiTracks);
          const float jpsiProb = vertexProbability(jpsiVtx);
          if (!(jpsiProb > 0.001)) continue;

          std::vector<reco::TransientTrack> fourTracks{ttBuilder->build(ePlusPtr->gsfTrack()), ttBuilder->build(eMinusPtr->gsfTrack()), ttBuilder->build(muPlusPtr->innerTrack()), ttBuilder->build(muMinusPtr->innerTrack())};
          const TransientVertex fourVtx = kvf.vertex(fourTracks);
          const float fourProb = vertexProbability(fourVtx);
          if (!(fourProb > 0.001)) continue;

          const TLorentzVector fourL = Zee + Jpsi;
          const std::vector<TLorentzVector> selectedLeptons{EPlus, EMinus, MUP, MUM};
          const float zIso03 = pairTrackIso(*packedCandHandle, Zee, selectedLeptons, 0.3f);
          const float zIso04 = pairTrackIso(*packedCandHandle, Zee, selectedLeptons, 0.4f);
          const float jIso03 = pairTrackIso(*packedCandHandle, Jpsi, selectedLeptons, 0.3f);
          const float jIso04 = pairTrackIso(*packedCandHandle, Jpsi, selectedLeptons, 0.4f);

          int nExtraE = 0;
          for (const auto idx : looseElectronIdx) if (idx != ie1 && idx != ie2) ++nExtraE;
          int nExtraMu = 0;
          for (const auto idx : looseMuonIdx) if (idx != im1 && idx != im2) ++nExtraMu;

          Run.push_back(iEvent.id().run());
          LumiBlock.push_back(iEvent.luminosityBlock());
          Event.push_back(static_cast<unsigned long long>(iEvent.id().event()));
          nPV.push_back(primaryVerticesHandle->size());
          passEleTrigger.push_back(firedEleTrigger ? 1 : 0);
          passEleTriggerMatch.push_back(hasEleTriggerMatch ? 1 : 0);

          pushP4(fourL_mass, fourL_px, fourL_py, fourL_pz, fourL_pt, fourL_eta, fourL_phi, fourL_rapidity, fourL);
          fourL_vtxProb.push_back(fourProb);
          fourL_pvx.push_back(bestVtx.x()); fourL_pvy.push_back(bestVtx.y()); fourL_pvz.push_back(bestVtx.z());
          fourL_pvxError.push_back(bestVtx.xError()); fourL_pvyError.push_back(bestVtx.yError()); fourL_pvzError.push_back(bestVtx.zError());

          pushP4(Z_mass, Z_px, Z_py, Z_pz, Z_pt, Z_eta, Z_phi, Z_rapidity, Zee);
          Z_vtxProb.push_back(zeeProb);
          Z_dR_ee.push_back(reco::deltaR(EPlus.Eta(), EPlus.Phi(), EMinus.Eta(), EMinus.Phi()));
          Z_trackIso03.push_back(zIso03); Z_trackIso04.push_back(zIso04); Z_relIso03.push_back(safeRatio(zIso03, Zee.Pt())); Z_relIso04.push_back(safeRatio(zIso04, Zee.Pt()));

          pushP4(Jpsi_mass, Jpsi_px, Jpsi_py, Jpsi_pz, Jpsi_pt, Jpsi_eta, Jpsi_phi, Jpsi_rapidity, Jpsi);
          Jpsi_vtxProb.push_back(jpsiProb);
          Jpsi_dR_mumu.push_back(reco::deltaR(MUP.Eta(), MUP.Phi(), MUM.Eta(), MUM.Phi()));
          Jpsi_trackIso03.push_back(jIso03); Jpsi_trackIso04.push_back(jIso04); Jpsi_relIso03.push_back(safeRatio(jIso03, Jpsi.Pt())); Jpsi_relIso04.push_back(safeRatio(jIso04, Jpsi.Pt()));

          Z_Jpsi_dR.push_back(reco::deltaR(Zee.Eta(), Zee.Phi(), Jpsi.Eta(), Jpsi.Phi()));
          Z_Jpsi_dPhi.push_back(reco::deltaPhi(Zee.Phi(), Jpsi.Phi()));
          Z_Jpsi_dEta.push_back(Zee.Eta() - Jpsi.Eta());
          Z_Jpsi_dY.push_back(Zee.Rapidity() - Jpsi.Rapidity());
          pt_balance.push_back(safeRatio(std::abs(static_cast<float>(Zee.Pt() - Jpsi.Pt())), static_cast<float>(Zee.Pt() + Jpsi.Pt())));
          cosTheta_Z_ePlus.push_back(cosThetaLeptonInPairRest(EPlus, Zee, fourL));
          cosTheta_Jpsi_muPlus.push_back(cosThetaLeptonInPairRest(MUP, Jpsi, fourL));
          phi_decayPlane_Z_Jpsi.push_back(decayPlaneAngle(EPlus, EMinus, MUP, MUM, fourL));

          pushObjP4(e1_px, e1_py, e1_pz, e1_pt, e1_eta, e1_phi, ELead);
          pushObjP4(e2_px, e2_py, e2_pz, e2_pt, e2_eta, e2_phi, ESub);
          e1_charge.push_back(eLeadPtr->charge()); e2_charge.push_back(eSubPtr->charge());
          e1_dxy.push_back(dxy(eLeadPtr->gsfTrack(), bestVtx)); e1_dz.push_back(dz(eLeadPtr->gsfTrack(), bestVtx));
          e2_dxy.push_back(dxy(eSubPtr->gsfTrack(), bestVtx)); e2_dz.push_back(dz(eSubPtr->gsfTrack(), bestVtx));
          e1_mvaRaw.push_back(electronMVARaw(*eLeadPtr, kElectronMVARawNames));
          e2_mvaRaw.push_back(electronMVARaw(*eSubPtr, kElectronMVARawNames));
          e1_passLooseID.push_back(electronID(*eLeadPtr, kLooseElectronID) > 0.5f ? 1 : 0);
          e2_passLooseID.push_back(electronID(*eSubPtr, kLooseElectronID) > 0.5f ? 1 : 0);
          e1_passWP90.push_back(electronID(*eLeadPtr, kElectronWP90ID) > 0.5f ? 1 : 0);
          e2_passWP90.push_back(electronID(*eSubPtr, kElectronWP90ID) > 0.5f ? 1 : 0);
          e1_passWP80.push_back(electronID(*eLeadPtr, kElectronWP80ID) > 0.5f ? 1 : 0);
          e2_passWP80.push_back(electronID(*eSubPtr, kElectronWP80ID) > 0.5f ? 1 : 0);
          e1_triggerMatched.push_back(eLeadTrigMatch ? 1 : 0); e2_triggerMatched.push_back(eSubTrigMatch ? 1 : 0);
          e1_triggerDR.push_back(eLeadTrigDR); e2_triggerDR.push_back(eSubTrigDR);

          pushObjP4(mu1_px, mu1_py, mu1_pz, mu1_pt, mu1_eta, mu1_phi, MuLead);
          pushObjP4(mu2_px, mu2_py, mu2_pz, mu2_pt, mu2_eta, mu2_phi, MuSub);
          mu1_charge.push_back(muLeadPtr->charge()); mu2_charge.push_back(muSubPtr->charge());
          mu1_soft.push_back(muLeadPtr->isSoftMuon(bestVtx) ? 1 : 0); mu2_soft.push_back(muSubPtr->isSoftMuon(bestVtx) ? 1 : 0);
          mu1_loose.push_back(muLeadPtr->isLooseMuon() ? 1 : 0); mu2_loose.push_back(muSubPtr->isLooseMuon() ? 1 : 0);
          mu1_tight.push_back(muLeadPtr->isTightMuon(bestVtx) ? 1 : 0); mu2_tight.push_back(muSubPtr->isTightMuon(bestVtx) ? 1 : 0);
          mu1_pfAbsIso03.push_back(pfAbsIso03(*muLeadPtr)); mu2_pfAbsIso03.push_back(pfAbsIso03(*muSubPtr));
          mu1_pfRelIso03.push_back(pfRelIso03(*muLeadPtr)); mu2_pfRelIso03.push_back(pfRelIso03(*muSubPtr));
          mu1_dxy.push_back(dxy(muLeadPtr->innerTrack(), bestVtx)); mu2_dxy.push_back(dxy(muSubPtr->innerTrack(), bestVtx));
          mu1_dz.push_back(dz(muLeadPtr->innerTrack(), bestVtx)); mu2_dz.push_back(dz(muSubPtr->innerTrack(), bestVtx));
          mu1_dB3D.push_back(muLeadPtr->dB(pat::Muon::PV3D)); mu2_dB3D.push_back(muSubPtr->dB(pat::Muon::PV3D));
          mu1_normChi2.push_back(trackChi2(muLeadPtr->innerTrack())); mu2_normChi2.push_back(trackChi2(muSubPtr->innerTrack()));
          mu1_nValidHits.push_back(nValidHits(muLeadPtr->innerTrack())); mu2_nValidHits.push_back(nValidHits(muSubPtr->innerTrack()));
          mu1_nValidPixelHits.push_back(nValidPixelHits(muLeadPtr->innerTrack())); mu2_nValidPixelHits.push_back(nValidPixelHits(muSubPtr->innerTrack()));

          nExtraLooseElectrons.push_back(nExtraE);
          nExtraLooseMuons.push_back(nExtraMu);

          ++nB;
        }
      }
    }
  }

  if (nB > 0 || keepEmptyEvents_) tree_->Fill();
}

void miniAODeemm::endJob() {}

DEFINE_FWK_MODULE(miniAODeemm);
