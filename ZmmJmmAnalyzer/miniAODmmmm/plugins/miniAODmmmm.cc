// -*- C++ -*-
// Clean zmmjmm MiniAOD preselection ntuplizer.
// Electron, trigger-object, kinematic-fit, and gen-particle dependencies removed.

#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/miniAODmmmm.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
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
constexpr float kMuonMass = 0.1056583745f;
constexpr float kBad = -999.0f;

float safeRatio(float numerator, float denominator) {
  if (std::abs(denominator) < 1e-6f) return kBad;
  return numerator / denominator;
}

float pfRelIso03(const pat::Muon& mu) {
  const auto& iso = mu.pfIsolationR03();
  const float absIso = iso.sumChargedHadronPt + std::max(0.0f, static_cast<float>(iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5 * iso.sumPUPt));
  return safeRatio(absIso, mu.pt());
}

float pfRelIso04(const pat::Muon& mu) {
  const auto& iso = mu.pfIsolationR04();
  const float absIso = iso.sumChargedHadronPt + std::max(0.0f, static_cast<float>(iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5 * iso.sumPUPt));
  return safeRatio(absIso, mu.pt());
}

float pfAbsIso03(const pat::Muon& mu) {
  const auto& iso = mu.pfIsolationR03();
  return iso.sumChargedHadronPt + std::max(0.0f, static_cast<float>(iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5 * iso.sumPUPt));
}

float pfAbsIso04(const pat::Muon& mu) {
  const auto& iso = mu.pfIsolationR04();
  return iso.sumChargedHadronPt + std::max(0.0f, static_cast<float>(iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5 * iso.sumPUPt));
}

float vertexProbability(const TransientVertex& vtx) {
  if (!vtx.isValid()) return kBad;
  return ChiSquaredProbability(vtx.totalChiSquared(), static_cast<int>(std::round(vtx.degreesOfFreedom())));
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

float pairTrackIso(const edm::View<pat::PackedCandidate>& packed,
                   const TLorentzVector& pair,
                   const std::vector<TLorentzVector>& selectedMuons,
                   float cone) {
  float sumPt = 0.0f;
  for (const auto& cand : packed) {
    if (cand.charge() == 0) continue;
    if (cand.pt() < 0.5f) continue;
    if (cand.fromPV() <= 1) continue;
    if (reco::deltaR(cand.eta(), cand.phi(), pair.Eta(), pair.Phi()) > cone) continue;

    bool overlapsSelectedMuon = false;
    for (const auto& mu : selectedMuons) {
      if (reco::deltaR(cand.eta(), cand.phi(), mu.Eta(), mu.Phi()) < 1.0e-3) {
        overlapsSelectedMuon = true;
        break;
      }
    }
    if (overlapsSelectedMuon) continue;
    sumPt += cand.pt();
  }
  return sumPt;
}

float cosThetaMuonInPairRest(const TLorentzVector& muPlus, const TLorentzVector& pair, const TLorentzVector& fourMu) {
  if (pair.P() <= 0.0 || fourMu.M() <= 0.0) return kBad;

  TLorentzVector muRest = muPlus;
  muRest.Boost(-pair.BoostVector());

  TLorentzVector pairInFourRest = pair;
  pairInFourRest.Boost(-fourMu.BoostVector());

  const TVector3 axis = pairInFourRest.Vect().Unit();
  const TVector3 muDir = muRest.Vect().Unit();
  if (axis.Mag() == 0.0 || muDir.Mag() == 0.0) return kBad;
  return muDir.Dot(axis);
}

float decayPlaneAngle(const TLorentzVector& a1, const TLorentzVector& a2, const TLorentzVector& b1, const TLorentzVector& b2, const TLorentzVector& fourMu) {
  TLorentzVector aa1 = a1;
  TLorentzVector aa2 = a2;
  TLorentzVector bb1 = b1;
  TLorentzVector bb2 = b2;
  const TVector3 boost = -fourMu.BoostVector();
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

template <typename T>
void clearVector(T*& ptr) {
  ptr->clear();
}

}  // namespace

miniAODmmmm::miniAODmmmm(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
      packedCandToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
      primaryVerticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      muonTriggerString_(iConfig.getParameter<std::string>("MuonTrigger")),
      requireTrigger_(iConfig.existsAs<bool>("requireTrigger") ? iConfig.getParameter<bool>("requireTrigger") : false),
      keepEmptyEvents_(iConfig.existsAs<bool>("keepEmptyEvents") ? iConfig.getParameter<bool>("keepEmptyEvents") : false) {}

void miniAODmmmm::beginJob() {
  edm::Service<TFileService> fs;
  if (!fs) throw cms::Exception("TFileService") << "TFileService is required.";
  // Put the tree at the ROOT file top level instead of inside the module-label TDirectory.
  fs->file().cd();
  tree_ = new TTree("ntuple", "ntuple");
  tree_->SetDirectory(&fs->file());

#define MAKE_VEC(NAME, TYPE) \
  NAME = new std::vector<TYPE>(); \
  tree_->Branch(#NAME, &NAME)

  tree_->Branch("nB", &nB, "nB/i");

  MAKE_VEC(Run, unsigned int);
  MAKE_VEC(LumiBlock, unsigned int);
  MAKE_VEC(Event, unsigned long long);
  MAKE_VEC(nPV, unsigned int);

  MAKE_VEC(fourMu_mass, float);
  MAKE_VEC(fourMu_px, float);
  MAKE_VEC(fourMu_py, float);
  MAKE_VEC(fourMu_pz, float);
  MAKE_VEC(fourMu_pt, float);
  MAKE_VEC(fourMu_eta, float);
  MAKE_VEC(fourMu_phi, float);
  MAKE_VEC(fourMu_rapidity, float);
  MAKE_VEC(fourMu_vtxProb, float);
  MAKE_VEC(fourMu_pvx, float);
  MAKE_VEC(fourMu_pvy, float);
  MAKE_VEC(fourMu_pvz, float);
  MAKE_VEC(fourMu_pvxError, float);
  MAKE_VEC(fourMu_pvyError, float);
  MAKE_VEC(fourMu_pvzError, float);

  MAKE_VEC(passMuonTrigger, bool);

  MAKE_VEC(pair12_mass, float); MAKE_VEC(pair12_px, float); MAKE_VEC(pair12_py, float); MAKE_VEC(pair12_pz, float); MAKE_VEC(pair12_pt, float); MAKE_VEC(pair12_eta, float); MAKE_VEC(pair12_phi, float); MAKE_VEC(pair12_rapidity, float);
  MAKE_VEC(pair34_mass, float); MAKE_VEC(pair34_px, float); MAKE_VEC(pair34_py, float); MAKE_VEC(pair34_pz, float); MAKE_VEC(pair34_pt, float); MAKE_VEC(pair34_eta, float); MAKE_VEC(pair34_phi, float); MAKE_VEC(pair34_rapidity, float);
  MAKE_VEC(pair23_mass, float); MAKE_VEC(pair23_px, float); MAKE_VEC(pair23_py, float); MAKE_VEC(pair23_pz, float); MAKE_VEC(pair23_pt, float); MAKE_VEC(pair23_eta, float); MAKE_VEC(pair23_phi, float); MAKE_VEC(pair23_rapidity, float);
  MAKE_VEC(pair14_mass, float); MAKE_VEC(pair14_px, float); MAKE_VEC(pair14_py, float); MAKE_VEC(pair14_pz, float); MAKE_VEC(pair14_pt, float); MAKE_VEC(pair14_eta, float); MAKE_VEC(pair14_phi, float); MAKE_VEC(pair14_rapidity, float);

  MAKE_VEC(pair12_vtxProb, float); MAKE_VEC(pair34_vtxProb, float); MAKE_VEC(pair23_vtxProb, float); MAKE_VEC(pair14_vtxProb, float);
  MAKE_VEC(pair12_dR, float); MAKE_VEC(pair34_dR, float); MAKE_VEC(pair23_dR, float); MAKE_VEC(pair14_dR, float);
  MAKE_VEC(pair12_34_dR, float); MAKE_VEC(pair23_14_dR, float); MAKE_VEC(pair12_34_dPhi, float); MAKE_VEC(pair23_14_dPhi, float); MAKE_VEC(pair12_34_dEta, float); MAKE_VEC(pair23_14_dEta, float); MAKE_VEC(pair12_34_dY, float); MAKE_VEC(pair23_14_dY, float);

  MAKE_VEC(pair12_trackIso03, float); MAKE_VEC(pair12_trackIso04, float); MAKE_VEC(pair12_relIso03, float); MAKE_VEC(pair12_relIso04, float);
  MAKE_VEC(pair34_trackIso03, float); MAKE_VEC(pair34_trackIso04, float); MAKE_VEC(pair34_relIso03, float); MAKE_VEC(pair34_relIso04, float);
  MAKE_VEC(pair23_trackIso03, float); MAKE_VEC(pair23_trackIso04, float); MAKE_VEC(pair23_relIso03, float); MAKE_VEC(pair23_relIso04, float);
  MAKE_VEC(pair14_trackIso03, float); MAKE_VEC(pair14_trackIso04, float); MAKE_VEC(pair14_relIso03, float); MAKE_VEC(pair14_relIso04, float);

  MAKE_VEC(pair12_cosThetaMu, float); MAKE_VEC(pair34_cosThetaMu, float); MAKE_VEC(pair23_cosThetaMu, float); MAKE_VEC(pair14_cosThetaMu, float);
  MAKE_VEC(pair12_34_phiDecayPlane, float); MAKE_VEC(pair23_14_phiDecayPlane, float);

  MAKE_VEC(muP1_px, float); MAKE_VEC(muP1_py, float); MAKE_VEC(muP1_pz, float); MAKE_VEC(muP1_pt, float); MAKE_VEC(muP1_eta, float); MAKE_VEC(muP1_phi, float); MAKE_VEC(muP1_charge, int); MAKE_VEC(muP1_soft, bool); MAKE_VEC(muP1_tight, bool); MAKE_VEC(muP1_loose, bool); MAKE_VEC(muP1_isoTrack, float); MAKE_VEC(muP1_isoHcal, float); MAKE_VEC(muP1_isoEcal, float); MAKE_VEC(muP1_isoCalo, float); MAKE_VEC(muP1_pfAbsIso03, float); MAKE_VEC(muP1_pfAbsIso04, float); MAKE_VEC(muP1_pfRelIso03, float); MAKE_VEC(muP1_pfRelIso04, float); MAKE_VEC(muP1_dB3D, float);
  MAKE_VEC(muM1_px, float); MAKE_VEC(muM1_py, float); MAKE_VEC(muM1_pz, float); MAKE_VEC(muM1_pt, float); MAKE_VEC(muM1_eta, float); MAKE_VEC(muM1_phi, float); MAKE_VEC(muM1_charge, int); MAKE_VEC(muM1_soft, bool); MAKE_VEC(muM1_tight, bool); MAKE_VEC(muM1_loose, bool); MAKE_VEC(muM1_isoTrack, float); MAKE_VEC(muM1_isoHcal, float); MAKE_VEC(muM1_isoEcal, float); MAKE_VEC(muM1_isoCalo, float); MAKE_VEC(muM1_pfAbsIso03, float); MAKE_VEC(muM1_pfAbsIso04, float); MAKE_VEC(muM1_pfRelIso03, float); MAKE_VEC(muM1_pfRelIso04, float); MAKE_VEC(muM1_dB3D, float);
  MAKE_VEC(muP2_px, float); MAKE_VEC(muP2_py, float); MAKE_VEC(muP2_pz, float); MAKE_VEC(muP2_pt, float); MAKE_VEC(muP2_eta, float); MAKE_VEC(muP2_phi, float); MAKE_VEC(muP2_charge, int); MAKE_VEC(muP2_soft, bool); MAKE_VEC(muP2_tight, bool); MAKE_VEC(muP2_loose, bool); MAKE_VEC(muP2_isoTrack, float); MAKE_VEC(muP2_isoHcal, float); MAKE_VEC(muP2_isoEcal, float); MAKE_VEC(muP2_isoCalo, float); MAKE_VEC(muP2_pfAbsIso03, float); MAKE_VEC(muP2_pfAbsIso04, float); MAKE_VEC(muP2_pfRelIso03, float); MAKE_VEC(muP2_pfRelIso04, float); MAKE_VEC(muP2_dB3D, float);
  MAKE_VEC(muM2_px, float); MAKE_VEC(muM2_py, float); MAKE_VEC(muM2_pz, float); MAKE_VEC(muM2_pt, float); MAKE_VEC(muM2_eta, float); MAKE_VEC(muM2_phi, float); MAKE_VEC(muM2_charge, int); MAKE_VEC(muM2_soft, bool); MAKE_VEC(muM2_tight, bool); MAKE_VEC(muM2_loose, bool); MAKE_VEC(muM2_isoTrack, float); MAKE_VEC(muM2_isoHcal, float); MAKE_VEC(muM2_isoEcal, float); MAKE_VEC(muM2_isoCalo, float); MAKE_VEC(muM2_pfAbsIso03, float); MAKE_VEC(muM2_pfAbsIso04, float); MAKE_VEC(muM2_pfRelIso03, float); MAKE_VEC(muM2_pfRelIso04, float); MAKE_VEC(muM2_dB3D, float);

  MAKE_VEC(muP1_dxy, float); MAKE_VEC(muM1_dxy, float); MAKE_VEC(muP1_dz, float); MAKE_VEC(muM1_dz, float); MAKE_VEC(muP2_dxy, float); MAKE_VEC(muM2_dxy, float); MAKE_VEC(muP2_dz, float); MAKE_VEC(muM2_dz, float);

  MAKE_VEC(muM1_normChi2, float); MAKE_VEC(muM1_nValidHits, int); MAKE_VEC(muM1_nValidPixelHits, int); MAKE_VEC(muP1_normChi2, float); MAKE_VEC(muP1_nValidHits, int); MAKE_VEC(muP1_nValidPixelHits, int);
  MAKE_VEC(muM2_normChi2, float); MAKE_VEC(muM2_nValidHits, int); MAKE_VEC(muM2_nValidPixelHits, int); MAKE_VEC(muP2_normChi2, float); MAKE_VEC(muP2_nValidHits, int); MAKE_VEC(muP2_nValidPixelHits, int);

  MAKE_VEC(muP1_p4_pt, float); MAKE_VEC(muP1_p4_eta, float); MAKE_VEC(muP1_p4_phi, float); MAKE_VEC(muP1_p4_px, float); MAKE_VEC(muP1_p4_py, float); MAKE_VEC(muP1_p4_pz, float);
  MAKE_VEC(muM1_p4_pt, float); MAKE_VEC(muM1_p4_eta, float); MAKE_VEC(muM1_p4_phi, float); MAKE_VEC(muM1_p4_px, float); MAKE_VEC(muM1_p4_py, float); MAKE_VEC(muM1_p4_pz, float);
  MAKE_VEC(muP2_p4_pt, float); MAKE_VEC(muP2_p4_eta, float); MAKE_VEC(muP2_p4_phi, float); MAKE_VEC(muP2_p4_px, float); MAKE_VEC(muP2_p4_py, float); MAKE_VEC(muP2_p4_pz, float);
  MAKE_VEC(muM2_p4_pt, float); MAKE_VEC(muM2_p4_eta, float); MAKE_VEC(muM2_p4_phi, float); MAKE_VEC(muM2_p4_px, float); MAKE_VEC(muM2_p4_py, float); MAKE_VEC(muM2_p4_pz, float);

#undef MAKE_VEC
}

void miniAODmmmm::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
#define CLEAR(NAME) NAME->clear()
  nB = 0;
  CLEAR(Run); CLEAR(LumiBlock); CLEAR(Event); CLEAR(nPV);
  CLEAR(fourMu_mass); CLEAR(fourMu_px); CLEAR(fourMu_py); CLEAR(fourMu_pz); CLEAR(fourMu_pt); CLEAR(fourMu_eta); CLEAR(fourMu_phi); CLEAR(fourMu_rapidity); CLEAR(fourMu_vtxProb); CLEAR(fourMu_pvx); CLEAR(fourMu_pvy); CLEAR(fourMu_pvz); CLEAR(fourMu_pvxError); CLEAR(fourMu_pvyError); CLEAR(fourMu_pvzError);
  CLEAR(passMuonTrigger);
  CLEAR(pair12_mass); CLEAR(pair12_px); CLEAR(pair12_py); CLEAR(pair12_pz); CLEAR(pair12_pt); CLEAR(pair12_eta); CLEAR(pair12_phi); CLEAR(pair12_rapidity); CLEAR(pair34_mass); CLEAR(pair34_px); CLEAR(pair34_py); CLEAR(pair34_pz); CLEAR(pair34_pt); CLEAR(pair34_eta); CLEAR(pair34_phi); CLEAR(pair34_rapidity); CLEAR(pair23_mass); CLEAR(pair23_px); CLEAR(pair23_py); CLEAR(pair23_pz); CLEAR(pair23_pt); CLEAR(pair23_eta); CLEAR(pair23_phi); CLEAR(pair23_rapidity); CLEAR(pair14_mass); CLEAR(pair14_px); CLEAR(pair14_py); CLEAR(pair14_pz); CLEAR(pair14_pt); CLEAR(pair14_eta); CLEAR(pair14_phi); CLEAR(pair14_rapidity);
  CLEAR(pair12_vtxProb); CLEAR(pair34_vtxProb); CLEAR(pair23_vtxProb); CLEAR(pair14_vtxProb);
  CLEAR(pair12_dR); CLEAR(pair34_dR); CLEAR(pair23_dR); CLEAR(pair14_dR); CLEAR(pair12_34_dR); CLEAR(pair23_14_dR); CLEAR(pair12_34_dPhi); CLEAR(pair23_14_dPhi); CLEAR(pair12_34_dEta); CLEAR(pair23_14_dEta); CLEAR(pair12_34_dY); CLEAR(pair23_14_dY);
  CLEAR(pair12_trackIso03); CLEAR(pair12_trackIso04); CLEAR(pair12_relIso03); CLEAR(pair12_relIso04); CLEAR(pair34_trackIso03); CLEAR(pair34_trackIso04); CLEAR(pair34_relIso03); CLEAR(pair34_relIso04); CLEAR(pair23_trackIso03); CLEAR(pair23_trackIso04); CLEAR(pair23_relIso03); CLEAR(pair23_relIso04); CLEAR(pair14_trackIso03); CLEAR(pair14_trackIso04); CLEAR(pair14_relIso03); CLEAR(pair14_relIso04);
  CLEAR(pair12_cosThetaMu); CLEAR(pair34_cosThetaMu); CLEAR(pair23_cosThetaMu); CLEAR(pair14_cosThetaMu); CLEAR(pair12_34_phiDecayPlane); CLEAR(pair23_14_phiDecayPlane);
  CLEAR(muP1_px); CLEAR(muP1_py); CLEAR(muP1_pz); CLEAR(muP1_pt); CLEAR(muP1_eta); CLEAR(muP1_phi); CLEAR(muP1_charge); CLEAR(muP1_soft); CLEAR(muP1_tight); CLEAR(muP1_loose); CLEAR(muP1_isoTrack); CLEAR(muP1_isoHcal); CLEAR(muP1_isoEcal); CLEAR(muP1_isoCalo); CLEAR(muP1_pfAbsIso03); CLEAR(muP1_pfAbsIso04); CLEAR(muP1_pfRelIso03); CLEAR(muP1_pfRelIso04); CLEAR(muP1_dB3D);
  CLEAR(muM1_px); CLEAR(muM1_py); CLEAR(muM1_pz); CLEAR(muM1_pt); CLEAR(muM1_eta); CLEAR(muM1_phi); CLEAR(muM1_charge); CLEAR(muM1_soft); CLEAR(muM1_tight); CLEAR(muM1_loose); CLEAR(muM1_isoTrack); CLEAR(muM1_isoHcal); CLEAR(muM1_isoEcal); CLEAR(muM1_isoCalo); CLEAR(muM1_pfAbsIso03); CLEAR(muM1_pfAbsIso04); CLEAR(muM1_pfRelIso03); CLEAR(muM1_pfRelIso04); CLEAR(muM1_dB3D);
  CLEAR(muP2_px); CLEAR(muP2_py); CLEAR(muP2_pz); CLEAR(muP2_pt); CLEAR(muP2_eta); CLEAR(muP2_phi); CLEAR(muP2_charge); CLEAR(muP2_soft); CLEAR(muP2_tight); CLEAR(muP2_loose); CLEAR(muP2_isoTrack); CLEAR(muP2_isoHcal); CLEAR(muP2_isoEcal); CLEAR(muP2_isoCalo); CLEAR(muP2_pfAbsIso03); CLEAR(muP2_pfAbsIso04); CLEAR(muP2_pfRelIso03); CLEAR(muP2_pfRelIso04); CLEAR(muP2_dB3D);
  CLEAR(muM2_px); CLEAR(muM2_py); CLEAR(muM2_pz); CLEAR(muM2_pt); CLEAR(muM2_eta); CLEAR(muM2_phi); CLEAR(muM2_charge); CLEAR(muM2_soft); CLEAR(muM2_tight); CLEAR(muM2_loose); CLEAR(muM2_isoTrack); CLEAR(muM2_isoHcal); CLEAR(muM2_isoEcal); CLEAR(muM2_isoCalo); CLEAR(muM2_pfAbsIso03); CLEAR(muM2_pfAbsIso04); CLEAR(muM2_pfRelIso03); CLEAR(muM2_pfRelIso04); CLEAR(muM2_dB3D);
  CLEAR(muP1_dxy); CLEAR(muM1_dxy); CLEAR(muP1_dz); CLEAR(muM1_dz); CLEAR(muP2_dxy); CLEAR(muM2_dxy); CLEAR(muP2_dz); CLEAR(muM2_dz);
  CLEAR(muM1_normChi2); CLEAR(muM1_nValidHits); CLEAR(muM1_nValidPixelHits); CLEAR(muP1_normChi2); CLEAR(muP1_nValidHits); CLEAR(muP1_nValidPixelHits); CLEAR(muM2_normChi2); CLEAR(muM2_nValidHits); CLEAR(muM2_nValidPixelHits); CLEAR(muP2_normChi2); CLEAR(muP2_nValidHits); CLEAR(muP2_nValidPixelHits);
  CLEAR(muP1_p4_pt); CLEAR(muP1_p4_eta); CLEAR(muP1_p4_phi); CLEAR(muP1_p4_px); CLEAR(muP1_p4_py); CLEAR(muP1_p4_pz); CLEAR(muM1_p4_pt); CLEAR(muM1_p4_eta); CLEAR(muM1_p4_phi); CLEAR(muM1_p4_px); CLEAR(muM1_p4_py); CLEAR(muM1_p4_pz); CLEAR(muP2_p4_pt); CLEAR(muP2_p4_eta); CLEAR(muP2_p4_phi); CLEAR(muP2_p4_px); CLEAR(muP2_p4_py); CLEAR(muP2_p4_pz); CLEAR(muM2_p4_pt); CLEAR(muM2_p4_eta); CLEAR(muM2_p4_phi); CLEAR(muM2_p4_px); CLEAR(muM2_p4_py); CLEAR(muM2_p4_pz);
#undef CLEAR

  edm::Handle<edm::View<pat::PackedCandidate>> packedCandHandle;
  iEvent.getByToken(packedCandToken_, packedCandHandle);
  if (!packedCandHandle.isValid()) {
    edm::LogWarning("miniAODmmmm") << "No packedPFCandidates in event";
    return;
  }

  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);
  if (!muonHandle.isValid()) {
    edm::LogWarning("miniAODmmmm") << "No slimmedMuons in event";
    return;
  }

  edm::Handle<reco::VertexCollection> primaryVerticesHandle;
  iEvent.getByToken(primaryVerticesToken_, primaryVerticesHandle);
  if (!primaryVerticesHandle.isValid() || primaryVerticesHandle->empty()) {
    edm::LogWarning("miniAODmmmm") << "No primary vertices in event";
    return;
  }
  const reco::Vertex& bestVtx = primaryVerticesHandle->front();

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  if (!triggerBits.isValid()) {
    edm::LogWarning("miniAODmmmm") << "No HLT TriggerResults in event";
    return;
  }

  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);
  if (!ttBuilder.isValid()) {
    edm::LogWarning("miniAODmmmm") << "No TransientTrackBuilder in event";
    return;
  }

  bool firedMuTrigger = false;
  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerBits);
  for (unsigned int i = 0; i < triggerBits->size(); ++i) {
    if (triggerNames.triggerName(i).find(muonTriggerString_) != std::string::npos && triggerBits->accept(i)) {
      firedMuTrigger = true;
      break;
    }
  }
  if (requireTrigger_ && !firedMuTrigger) return;

  if (muonHandle->size() < 4) {
    if (keepEmptyEvents_) tree_->Fill();
    return;
  }

  KalmanVertexFitter kvf(true);

  for (auto iMuon1 = muonHandle->begin(); iMuon1 != muonHandle->end(); ++iMuon1) {
    for (auto iMuon2 = iMuon1 + 1; iMuon2 != muonHandle->end(); ++iMuon2) {
      for (auto iMuon3 = iMuon2 + 1; iMuon3 != muonHandle->end(); ++iMuon3) {
        for (auto iMuon4 = iMuon3 + 1; iMuon4 != muonHandle->end(); ++iMuon4) {
          if (std::abs(iMuon1->charge()) != 1 || std::abs(iMuon2->charge()) != 1 || std::abs(iMuon3->charge()) != 1 || std::abs(iMuon4->charge()) != 1) continue;
          if (iMuon1->charge() + iMuon2->charge() + iMuon3->charge() + iMuon4->charge() != 0) continue;
          if (iMuon1->pt() < 2.0 || iMuon2->pt() < 2.0 || iMuon3->pt() < 2.0 || iMuon4->pt() < 2.0) continue;

          reco::TrackRef trkP1, trkP2, trkM1, trkM2;
          const pat::Muon *muP1 = nullptr, *muP2 = nullptr, *muM1 = nullptr, *muM2 = nullptr;

          std::vector<const pat::Muon*> plusMuons;
          std::vector<const pat::Muon*> minusMuons;
          for (const pat::Muon* mu : {&(*iMuon1), &(*iMuon2), &(*iMuon3), &(*iMuon4)}) {
            if (mu->charge() > 0) plusMuons.push_back(mu);
            else minusMuons.push_back(mu);
          }
          if (plusMuons.size() != 2 || minusMuons.size() != 2) continue;
          muP1 = plusMuons[0];
          muP2 = plusMuons[1];
          muM1 = minusMuons[0];
          muM2 = minusMuons[1];
          trkP1 = muP1->track();
          trkP2 = muP2->track();
          trkM1 = muM1->track();
          trkM2 = muM2->track();
          if (trkP1.isNull() || trkP2.isNull() || trkM1.isNull() || trkM2.isNull()) continue;
          if (!trkP1->quality(reco::TrackBase::highPurity) || !trkP2->quality(reco::TrackBase::highPurity) || !trkM1->quality(reco::TrackBase::highPurity) || !trkM2->quality(reco::TrackBase::highPurity)) continue;

          TLorentzVector M1, M2, M3, M4;
          M1.SetXYZM(muP1->px(), muP1->py(), muP1->pz(), kMuonMass);
          M2.SetXYZM(muM1->px(), muM1->py(), muM1->pz(), kMuonMass);
          M3.SetXYZM(muP2->px(), muP2->py(), muP2->pz(), kMuonMass);
          M4.SetXYZM(muM2->px(), muM2->py(), muM2->pz(), kMuonMass);
          const TLorentzVector MM1 = M1 + M2;
          const TLorentzVector MM2 = M3 + M4;
          const TLorentzVector MM3 = M2 + M3;
          const TLorentzVector MM4 = M1 + M4;
          const TLorentzVector MMMM = M1 + M2 + M3 + M4;

          reco::TransientTrack ttP1 = ttBuilder->build(trkP1);
          reco::TransientTrack ttM1 = ttBuilder->build(trkM1);
          reco::TransientTrack ttP2 = ttBuilder->build(trkP2);
          reco::TransientTrack ttM2 = ttBuilder->build(trkM2);

          std::vector<reco::TransientTrack> tracks12{ttP1, ttM1};
          std::vector<reco::TransientTrack> tracks34{ttP2, ttM2};
          std::vector<reco::TransientTrack> tracks23{ttM1, ttP2};
          std::vector<reco::TransientTrack> tracks14{ttP1, ttM2};
          std::vector<reco::TransientTrack> tracks4{ttP1, ttM1, ttP2, ttM2};
          TransientVertex vtx1 = kvf.vertex(tracks12);
          TransientVertex vtx2 = kvf.vertex(tracks34);
          TransientVertex vtx3 = kvf.vertex(tracks23);
          TransientVertex vtx4 = kvf.vertex(tracks14);
          TransientVertex vtx4mu = kvf.vertex(tracks4);

          const std::vector<TLorentzVector> selectedMuonP4{M1, M2, M3, M4};
          const float iso1_03 = pairTrackIso(*packedCandHandle, MM1, selectedMuonP4, 0.3f);
          const float iso1_04 = pairTrackIso(*packedCandHandle, MM1, selectedMuonP4, 0.4f);
          const float iso2_03 = pairTrackIso(*packedCandHandle, MM2, selectedMuonP4, 0.3f);
          const float iso2_04 = pairTrackIso(*packedCandHandle, MM2, selectedMuonP4, 0.4f);
          const float iso3_03 = pairTrackIso(*packedCandHandle, MM3, selectedMuonP4, 0.3f);
          const float iso3_04 = pairTrackIso(*packedCandHandle, MM3, selectedMuonP4, 0.4f);
          const float iso4_03 = pairTrackIso(*packedCandHandle, MM4, selectedMuonP4, 0.3f);
          const float iso4_04 = pairTrackIso(*packedCandHandle, MM4, selectedMuonP4, 0.4f);

          Run->push_back(iEvent.id().run());
          LumiBlock->push_back(iEvent.luminosityBlock());
          Event->push_back(iEvent.id().event());
          nPV->push_back(primaryVerticesHandle->size());

          fourMu_mass->push_back(MMMM.M());
          fourMu_px->push_back(MMMM.Px());
          fourMu_py->push_back(MMMM.Py());
          fourMu_pz->push_back(MMMM.Pz());
          fourMu_pt->push_back(MMMM.Pt());
          fourMu_eta->push_back(MMMM.Eta());
          fourMu_phi->push_back(MMMM.Phi());
          fourMu_rapidity->push_back(MMMM.Rapidity());
          fourMu_vtxProb->push_back(vertexProbability(vtx4mu));
          fourMu_pvx->push_back(bestVtx.x());
          fourMu_pvy->push_back(bestVtx.y());
          fourMu_pvz->push_back(bestVtx.z());
          fourMu_pvxError->push_back(bestVtx.xError());
          fourMu_pvyError->push_back(bestVtx.yError());
          fourMu_pvzError->push_back(bestVtx.zError());
          passMuonTrigger->push_back(firedMuTrigger);

#define PUSH_P4(PREFIX, P4) \
          PREFIX##_mass->push_back((P4).M()); \
          PREFIX##_px->push_back((P4).Px()); \
          PREFIX##_py->push_back((P4).Py()); \
          PREFIX##_pz->push_back((P4).Pz()); \
          PREFIX##_pt->push_back((P4).Pt()); \
          PREFIX##_eta->push_back((P4).Eta()); \
          PREFIX##_phi->push_back((P4).Phi()); \
          PREFIX##_rapidity->push_back((P4).Rapidity())
          PUSH_P4(pair12, MM1);
          PUSH_P4(pair34, MM2);
          PUSH_P4(pair23, MM3);
          PUSH_P4(pair14, MM4);
#undef PUSH_P4

          pair12_vtxProb->push_back(vertexProbability(vtx1));
          pair34_vtxProb->push_back(vertexProbability(vtx2));
          pair23_vtxProb->push_back(vertexProbability(vtx3));
          pair14_vtxProb->push_back(vertexProbability(vtx4));

          pair12_dR->push_back(reco::deltaR(M1.Eta(), M1.Phi(), M2.Eta(), M2.Phi()));
          pair34_dR->push_back(reco::deltaR(M3.Eta(), M3.Phi(), M4.Eta(), M4.Phi()));
          pair23_dR->push_back(reco::deltaR(M2.Eta(), M2.Phi(), M3.Eta(), M3.Phi()));
          pair14_dR->push_back(reco::deltaR(M1.Eta(), M1.Phi(), M4.Eta(), M4.Phi()));
          pair12_34_dR->push_back(reco::deltaR(MM1.Eta(), MM1.Phi(), MM2.Eta(), MM2.Phi()));
          pair23_14_dR->push_back(reco::deltaR(MM3.Eta(), MM3.Phi(), MM4.Eta(), MM4.Phi()));
          pair12_34_dPhi->push_back(reco::deltaPhi(MM1.Phi(), MM2.Phi()));
          pair23_14_dPhi->push_back(reco::deltaPhi(MM3.Phi(), MM4.Phi()));
          pair12_34_dEta->push_back(MM1.Eta() - MM2.Eta());
          pair23_14_dEta->push_back(MM3.Eta() - MM4.Eta());
          pair12_34_dY->push_back(MM1.Rapidity() - MM2.Rapidity());
          pair23_14_dY->push_back(MM3.Rapidity() - MM4.Rapidity());

          pair12_trackIso03->push_back(iso1_03); pair12_trackIso04->push_back(iso1_04); pair12_relIso03->push_back(safeRatio(iso1_03, MM1.Pt())); pair12_relIso04->push_back(safeRatio(iso1_04, MM1.Pt()));
          pair34_trackIso03->push_back(iso2_03); pair34_trackIso04->push_back(iso2_04); pair34_relIso03->push_back(safeRatio(iso2_03, MM2.Pt())); pair34_relIso04->push_back(safeRatio(iso2_04, MM2.Pt()));
          pair23_trackIso03->push_back(iso3_03); pair23_trackIso04->push_back(iso3_04); pair23_relIso03->push_back(safeRatio(iso3_03, MM3.Pt())); pair23_relIso04->push_back(safeRatio(iso3_04, MM3.Pt()));
          pair14_trackIso03->push_back(iso4_03); pair14_trackIso04->push_back(iso4_04); pair14_relIso03->push_back(safeRatio(iso4_03, MM4.Pt())); pair14_relIso04->push_back(safeRatio(iso4_04, MM4.Pt()));

          pair12_cosThetaMu->push_back(cosThetaMuonInPairRest(M1, MM1, MMMM));
          pair34_cosThetaMu->push_back(cosThetaMuonInPairRest(M3, MM2, MMMM));
          pair23_cosThetaMu->push_back(cosThetaMuonInPairRest(M3, MM3, MMMM));
          pair14_cosThetaMu->push_back(cosThetaMuonInPairRest(M1, MM4, MMMM));
          pair12_34_phiDecayPlane->push_back(decayPlaneAngle(M1, M2, M3, M4, MMMM));
          pair23_14_phiDecayPlane->push_back(decayPlaneAngle(M2, M3, M1, M4, MMMM));

#define PUSH_MUON(PREFIX, MU) \
          PREFIX##_px->push_back((MU)->px()); \
          PREFIX##_py->push_back((MU)->py()); \
          PREFIX##_pz->push_back((MU)->pz()); \
          PREFIX##_pt->push_back((MU)->pt()); \
          PREFIX##_eta->push_back((MU)->eta()); \
          PREFIX##_phi->push_back((MU)->phi()); \
          PREFIX##_charge->push_back((MU)->charge()); \
          PREFIX##_soft->push_back((MU)->isSoftMuon(bestVtx)); \
          PREFIX##_tight->push_back((MU)->isTightMuon(bestVtx)); \
          PREFIX##_loose->push_back((MU)->isLooseMuon()); \
          PREFIX##_isoTrack->push_back((MU)->isolationR03().sumPt); \
          PREFIX##_isoHcal->push_back((MU)->isolationR03().hadEt); \
          PREFIX##_isoEcal->push_back((MU)->isolationR03().emEt); \
          PREFIX##_isoCalo->push_back((MU)->isolationR03().hadEt + (MU)->isolationR03().emEt); \
          PREFIX##_pfAbsIso03->push_back(pfAbsIso03(*(MU))); \
          PREFIX##_pfAbsIso04->push_back(pfAbsIso04(*(MU))); \
          PREFIX##_pfRelIso03->push_back(pfRelIso03(*(MU))); \
          PREFIX##_pfRelIso04->push_back(pfRelIso04(*(MU))); \
          PREFIX##_dB3D->push_back((MU)->dB(pat::Muon::PV3D))
          PUSH_MUON(muP1, muP1);
          PUSH_MUON(muM1, muM1);
          PUSH_MUON(muP2, muP2);
          PUSH_MUON(muM2, muM2);
#undef PUSH_MUON

          muP1_dxy->push_back(dxy(trkP1, bestVtx)); muM1_dxy->push_back(dxy(trkM1, bestVtx)); muP1_dz->push_back(dz(trkP1, bestVtx)); muM1_dz->push_back(dz(trkM1, bestVtx));
          muP2_dxy->push_back(dxy(trkP2, bestVtx)); muM2_dxy->push_back(dxy(trkM2, bestVtx)); muP2_dz->push_back(dz(trkP2, bestVtx)); muM2_dz->push_back(dz(trkM2, bestVtx));

          muM1_normChi2->push_back(trackChi2(trkM1)); muM1_nValidHits->push_back(nValidHits(trkM1)); muM1_nValidPixelHits->push_back(nValidPixelHits(trkM1));
          muP1_normChi2->push_back(trackChi2(trkP1)); muP1_nValidHits->push_back(nValidHits(trkP1)); muP1_nValidPixelHits->push_back(nValidPixelHits(trkP1));
          muM2_normChi2->push_back(trackChi2(trkM2)); muM2_nValidHits->push_back(nValidHits(trkM2)); muM2_nValidPixelHits->push_back(nValidPixelHits(trkM2));
          muP2_normChi2->push_back(trackChi2(trkP2)); muP2_nValidHits->push_back(nValidHits(trkP2)); muP2_nValidPixelHits->push_back(nValidPixelHits(trkP2));

#define PUSH_MU_P4(PREFIX, P4) \
          PREFIX##_pt->push_back((P4).Pt()); \
          PREFIX##_eta->push_back((P4).Eta()); \
          PREFIX##_phi->push_back((P4).Phi()); \
          PREFIX##_px->push_back((P4).Px()); \
          PREFIX##_py->push_back((P4).Py()); \
          PREFIX##_pz->push_back((P4).Pz())
          PUSH_MU_P4(muP1_p4, M1);
          PUSH_MU_P4(muM1_p4, M2);
          PUSH_MU_P4(muP2_p4, M3);
          PUSH_MU_P4(muM2_p4, M4);
#undef PUSH_MU_P4

          ++nB;
        }
      }
    }
  }
  if (nB > 0 || keepEmptyEvents_) tree_->Fill();
}

void miniAODmmmm::endJob() {}

DEFINE_FWK_MODULE(miniAODmmmm);
