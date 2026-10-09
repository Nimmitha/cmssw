// mytagAndProbeV9: one MiniAOD pass for the J/psi long-term efficiency study (tracking, muon reconstruction, trigger,
// vertex), run on the Muon PDs (IsoMu24 tag: unbiased for every step) and/or on ParkingDoubleMuonLowMass (tag matched
// to HLT_Mu0_L1DoubleMu: analysis kinematics; unbiased for the tracking tree only).
//
// Tag: tight muon, pt > tagMinPt, matched (dR < 0.1) to a trigger object of one of tagPaths (bit i of tag_trig).
// Trees (one entry per tag-probe pair, opposite charge, pair mass in the window):
//   sa   probe = muon with a standalone (outer) track, kinematics from the outer track; pass = the analysis muon:
//        the same muon has a high-purity inner track (+ pixel and layer counts) -> tracking efficiency. Also the
//        nearest other muon with an inner track and the nearest high-purity packed/lost track (dR, pt ratio, pixel
//        hits), so the pass definition can be chosen offline (unmerged standalone duplicates)
//   trk  probe = high-purity packed track (packedPFCandidates + lostTracks) with pt > 3; pass = a pat::Muon with
//        that inner track -> muon reconstruction given a track; also soft/loose ID, L1/HLT matching, Kalman
//        vertex probability of tag + probe (as the analysis), event decision of the analysis path -> trigger, vertex
// Event: run, lumi, bx, good-PV count (analysis definition), analysis-path decision and HLT prescale.
// endJob prints events seen / with a tag / entries, events with a missing product, fired HLT_Mu*/HLT_IsoMu* paths.
// prescale: keep events with event % prescale == 0 (for the PDMLM pass; 1 = all). fillSA / fillTrk switch the trees
// (on PDMLM only sa is unbiased: the L1 DoubleMu seed needs the probe in the muon system).
#include <algorithm>
#include <iostream>
#include <map>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
#include "TTree.h"

namespace {
  constexpr double kMuMass = 0.1056583745;
}

class mytagAndProbeV9 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit mytagAndProbeV9(const edm::ParameterSet&);
  ~mytagAndProbeV9() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  struct Match {
    float dR = 99.f;
    int qual = -1;
  };
  Match matchL1(float eta, float phi, const l1t::MuonBxCollection&, float etaVeto, float phiVeto) const;
  float matchHLT(float eta, float phi, const std::vector<pat::TriggerObjectStandAlone>&) const;
  float vtxProb(const reco::Track&, const reco::Track&, const TransientTrackBuilder&) const;

  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostToken_;
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  edm::EDGetTokenT<edm::TriggerResults> bitsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> objectsToken_;
  edm::EDGetTokenT<l1t::MuonBxCollection> l1Token_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> prescalesToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;

  std::vector<std::string> tagPaths_;
  std::string analysisPath_;
  double tagMinPt_, probeMinPt_, saMinPt_, maxEta_;
  double massMin_, massMax_, saMassMin_, saMassMax_;
  unsigned prescale_;
  bool fillSA_, fillTrk_;

  TTree* sa_;
  TTree* trk_;
  // job summary (endJob): events seen, events with a missing product (per product), fired HLT_Mu*/HLT_IsoMu* paths
  unsigned long nSeen_ = 0, nTagEvents_ = 0, nSA_ = 0, nTrk_ = 0;
  std::map<std::string, unsigned long> missing_, fired_;

  // event (shared by both trees)
  UInt_t run_, lumi_;
  ULong64_t event_;
  Int_t bx_, nPV_;
  Bool_t passAnalysis_;
  Int_t psAnalysis_;
  // tag
  Float_t tag_pt_, tag_eta_, tag_phi_;
  Int_t tag_trig_;
  // pair
  Float_t mass_, pair_pt_, pair_dR_;
  // sa probe
  Float_t sa_pt_, sa_eta_, sa_phi_;
  Int_t sa_nStations_, sa_validMuonHits_;
  Bool_t sa_pass_, sa_hasInner_, sa_isTracker_, sa_isGlobal_;
  Float_t sa_inner_pt_, sa_massInner_;
  Int_t sa_pixelHits_, sa_pixelLayers_, sa_trackerLayers_, sa_bpixLayer1_;
  // alternatives to the same-object match: nearest other muon with an inner track, nearest high-purity track
  Float_t sa_otherMu_dR_, sa_otherMu_ptRatio_, sa_track_dR_, sa_track_ptRatio_;
  Int_t sa_otherMu_pixelHits_, sa_track_pixelHits_;
  Bool_t sa_otherMu_highPurity_;
  // trk probe
  Float_t p_pt_, p_eta_, p_phi_, p_dxy_, p_dz_;
  Int_t p_pixelHits_, p_trackerLayers_;
  Bool_t p_isLost_, p_isMuon_, p_isTracker_, p_isGlobal_, p_isSoft_, p_isLoose_;
  Float_t p_vtxProb_, p_l1dR_, p_hltdR_;
  Int_t p_l1qual_;
};

mytagAndProbeV9::mytagAndProbeV9(const edm::ParameterSet& cfg)
    : muonsToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
      pfToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("pfCands"))),
      lostToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
      verticesToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
      bitsToken_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("bits"))),
      objectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(cfg.getParameter<edm::InputTag>("objects"))),
      l1Token_(consumes<l1t::MuonBxCollection>(cfg.getParameter<edm::InputTag>("l1Muons"))),
      prescalesToken_(consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("prescales"))),
      ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      tagPaths_(cfg.getParameter<std::vector<std::string>>("tagPaths")),
      analysisPath_(cfg.getParameter<std::string>("analysisPath")),
      tagMinPt_(cfg.getParameter<double>("tagMinPt")),
      probeMinPt_(cfg.getParameter<double>("probeMinPt")),
      saMinPt_(cfg.getParameter<double>("saMinPt")),
      maxEta_(cfg.getParameter<double>("maxEta")),
      massMin_(cfg.getParameter<double>("massMin")),
      massMax_(cfg.getParameter<double>("massMax")),
      saMassMin_(cfg.getParameter<double>("saMassMin")),
      saMassMax_(cfg.getParameter<double>("saMassMax")),
      prescale_(cfg.getParameter<unsigned>("prescale")),
      fillSA_(cfg.getParameter<bool>("fillSA")),
      fillTrk_(cfg.getParameter<bool>("fillTrk")),
      sa_(nullptr),
      trk_(nullptr) {
  usesResource("TFileService");
}

mytagAndProbeV9::Match mytagAndProbeV9::matchL1(
    float eta, float phi, const l1t::MuonBxCollection& l1, float etaVeto, float phiVeto) const {
  // nearest BX-0 L1 muon (at the vertex), skipping the one nearest the tag (dR < 0.3)
  Match m;
  if (l1.isEmpty(0))
    return m;
  for (auto it = l1.begin(0); it != l1.end(0); ++it) {
    if (reco::deltaR(etaVeto, phiVeto, it->etaAtVtx(), it->phiAtVtx()) < 0.3)
      continue;
    const float d = reco::deltaR(eta, phi, it->etaAtVtx(), it->phiAtVtx());
    if (d < m.dR) {
      m.dR = d;
      m.qual = it->hwQual();
    }
  }
  return m;
}

float mytagAndProbeV9::matchHLT(float eta, float phi, const std::vector<pat::TriggerObjectStandAlone>& objs) const {
  float best = 99.f;
  for (const auto& o : objs)
    if (o.hasPathName(analysisPath_ + "*", false, true))
      best = std::min(best, static_cast<float>(reco::deltaR(eta, phi, o.eta(), o.phi())));
  return best;
}

float mytagAndProbeV9::vtxProb(const reco::Track& a, const reco::Track& b, const TransientTrackBuilder& ttb) const {
  std::vector<reco::TransientTrack> tt{ttb.build(a), ttb.build(b)};
  KalmanVertexFitter kvf(true);
  TransientVertex v = kvf.vertex(tt);
  if (!v.isValid())
    return -1.f;
  return ChiSquaredProbability(v.totalChiSquared(), v.degreesOfFreedom());
}

void mytagAndProbeV9::beginJob() {
  edm::Service<TFileService> fs;
  sa_ = fs->make<TTree>("sa", "standalone-muon probes: tracking efficiency");
  trk_ = fs->make<TTree>("trk", "track probes: muon reconstruction, ID, trigger, vertex");
  for (TTree* t : {sa_, trk_}) {
    t->Branch("run", &run_, "run/i");
    t->Branch("lumi", &lumi_, "lumi/i");
    t->Branch("event", &event_, "event/l");
    t->Branch("bx", &bx_, "bx/I");
    t->Branch("nPV", &nPV_, "nPV/I");
    t->Branch("passAnalysis", &passAnalysis_, "passAnalysis/O");
    t->Branch("psAnalysis", &psAnalysis_, "psAnalysis/I");
    t->Branch("tag_pt", &tag_pt_, "tag_pt/F");
    t->Branch("tag_eta", &tag_eta_, "tag_eta/F");
    t->Branch("tag_phi", &tag_phi_, "tag_phi/F");
    t->Branch("tag_trig", &tag_trig_, "tag_trig/I");
    t->Branch("mass", &mass_, "mass/F");
    t->Branch("pair_pt", &pair_pt_, "pair_pt/F");
    t->Branch("pair_dR", &pair_dR_, "pair_dR/F");
  }
  sa_->Branch("pt", &sa_pt_, "pt/F");
  sa_->Branch("eta", &sa_eta_, "eta/F");
  sa_->Branch("phi", &sa_phi_, "phi/F");
  sa_->Branch("nStations", &sa_nStations_, "nStations/I");
  sa_->Branch("validMuonHits", &sa_validMuonHits_, "validMuonHits/I");
  sa_->Branch("pass", &sa_pass_, "pass/O");
  sa_->Branch("hasInner", &sa_hasInner_, "hasInner/O");
  sa_->Branch("isTracker", &sa_isTracker_, "isTracker/O");
  sa_->Branch("isGlobal", &sa_isGlobal_, "isGlobal/O");
  sa_->Branch("inner_pt", &sa_inner_pt_, "inner_pt/F");
  sa_->Branch("massInner", &sa_massInner_, "massInner/F");
  sa_->Branch("pixelHits", &sa_pixelHits_, "pixelHits/I");
  sa_->Branch("pixelLayers", &sa_pixelLayers_, "pixelLayers/I");
  sa_->Branch("trackerLayers", &sa_trackerLayers_, "trackerLayers/I");
  sa_->Branch("bpixLayer1", &sa_bpixLayer1_, "bpixLayer1/I");
  sa_->Branch("otherMu_dR", &sa_otherMu_dR_, "otherMu_dR/F");
  sa_->Branch("otherMu_ptRatio", &sa_otherMu_ptRatio_, "otherMu_ptRatio/F");
  sa_->Branch("otherMu_pixelHits", &sa_otherMu_pixelHits_, "otherMu_pixelHits/I");
  sa_->Branch("otherMu_highPurity", &sa_otherMu_highPurity_, "otherMu_highPurity/O");
  sa_->Branch("track_dR", &sa_track_dR_, "track_dR/F");
  sa_->Branch("track_ptRatio", &sa_track_ptRatio_, "track_ptRatio/F");
  sa_->Branch("track_pixelHits", &sa_track_pixelHits_, "track_pixelHits/I");

  trk_->Branch("pt", &p_pt_, "pt/F");
  trk_->Branch("eta", &p_eta_, "eta/F");
  trk_->Branch("phi", &p_phi_, "phi/F");
  trk_->Branch("dxy", &p_dxy_, "dxy/F");
  trk_->Branch("dz", &p_dz_, "dz/F");
  trk_->Branch("pixelHits", &p_pixelHits_, "pixelHits/I");
  trk_->Branch("trackerLayers", &p_trackerLayers_, "trackerLayers/I");
  trk_->Branch("isLost", &p_isLost_, "isLost/O");
  trk_->Branch("isMuon", &p_isMuon_, "isMuon/O");
  trk_->Branch("isTracker", &p_isTracker_, "isTracker/O");
  trk_->Branch("isGlobal", &p_isGlobal_, "isGlobal/O");
  trk_->Branch("isSoft", &p_isSoft_, "isSoft/O");
  trk_->Branch("isLoose", &p_isLoose_, "isLoose/O");
  trk_->Branch("vtxProb", &p_vtxProb_, "vtxProb/F");
  trk_->Branch("l1dR", &p_l1dR_, "l1dR/F");
  trk_->Branch("l1qual", &p_l1qual_, "l1qual/I");
  trk_->Branch("hltdR", &p_hltdR_, "hltdR/F");
}

void mytagAndProbeV9::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (prescale_ > 1 && iEvent.id().event() % prescale_ != 0)
    return;
  edm::Handle<pat::MuonCollection> muons;
  edm::Handle<pat::PackedCandidateCollection> pf, lost;
  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<edm::TriggerResults> bits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> objects;
  edm::Handle<l1t::MuonBxCollection> l1;
  edm::Handle<pat::PackedTriggerPrescales> prescales;
  iEvent.getByToken(muonsToken_, muons);
  iEvent.getByToken(pfToken_, pf);
  iEvent.getByToken(lostToken_, lost);
  iEvent.getByToken(verticesToken_, vertices);
  iEvent.getByToken(bitsToken_, bits);
  iEvent.getByToken(objectsToken_, objects);
  iEvent.getByToken(l1Token_, l1);
  iEvent.getByToken(prescalesToken_, prescales);
  ++nSeen_;
  bool bad = false;
  for (const auto& [n, ok] : std::vector<std::pair<std::string, bool>>{
           {"muons", muons.isValid()}, {"pfCands", pf.isValid()}, {"lostTracks", lost.isValid()},
           {"vertices", vertices.isValid()}, {"bits", bits.isValid()}, {"objects", objects.isValid()},
           {"l1Muons", l1.isValid()}, {"prescales", prescales.isValid()}})
    if (!ok) {
      ++missing_[n];
      bad = true;
    }
  if (bad || vertices->empty())
    return;
  const reco::Vertex& pv = vertices->front();
  const auto& ttb = iSetup.getData(ttbToken_);

  run_ = iEvent.id().run();
  lumi_ = iEvent.id().luminosityBlock();
  event_ = iEvent.id().event();
  bx_ = iEvent.bunchCrossing();
  nPV_ = 0;
  for (const auto& v : *vertices)  // analysis definition (miniAODmmmm.cc)
    if (!v.isFake() && v.ndof() > 4 && std::fabs(v.z()) <= 24.0 && std::fabs(v.position().Rho()) <= 2.0)
      ++nPV_;

  const edm::TriggerNames& names = iEvent.triggerNames(*bits);
  passAnalysis_ = false;
  psAnalysis_ = -1;
  bool anyTag = false;
  for (unsigned i = 0; i < bits->size(); ++i) {
    const std::string& n = names.triggerName(i);
    if (n.rfind(analysisPath_, 0) == 0)
      psAnalysis_ = prescales->getPrescaleForIndex(i);
    if (!bits->accept(i))
      continue;
    if (n.rfind("HLT_Mu", 0) == 0 || n.rfind("HLT_IsoMu", 0) == 0)
      ++fired_[n.substr(0, n.rfind("_v"))];
    if (n.rfind(analysisPath_, 0) == 0)
      passAnalysis_ = true;
    for (const auto& p : tagPaths_)
      if (n.rfind(p, 0) == 0)
        anyTag = true;
  }
  if (!anyTag)
    return;
  std::vector<pat::TriggerObjectStandAlone> objs;
  objs.reserve(objects->size());
  for (auto o : *objects) {
    o.unpackPathNames(names);
    objs.push_back(o);
  }

  // tags
  struct Tag {
    const pat::Muon* mu;
    int trig;
  };
  std::vector<Tag> tags;
  for (const auto& mu : *muons) {
    if (mu.pt() < tagMinPt_ || std::fabs(mu.eta()) > maxEta_ || !muon::isTightMuon(mu, pv) || mu.innerTrack().isNull())
      continue;
    int trig = 0;
    for (unsigned k = 0; k < tagPaths_.size(); ++k)
      for (const auto& o : objs)
        if (o.hasPathName(tagPaths_[k] + "*", true, true) && reco::deltaR(mu, o) < 0.1) {
          trig |= (1 << k);
          break;
        }
    if (trig)
      tags.push_back({&mu, trig});
  }
  if (tags.empty())
    return;
  ++nTagEvents_;

  for (const auto& t : tags) {
    const pat::Muon& tag = *t.mu;
    tag_pt_ = tag.pt();
    tag_eta_ = tag.eta();
    tag_phi_ = tag.phi();
    tag_trig_ = t.trig;
    TLorentzVector vt;
    vt.SetPtEtaPhiM(tag.pt(), tag.eta(), tag.phi(), kMuMass);

    // --- sa tree: standalone probes
    for (const auto& mu : *muons) {
      if (!fillSA_)
        break;
      if (&mu == &tag || mu.outerTrack().isNull())
        continue;
      const reco::Track& o = *mu.outerTrack();
      if (o.pt() < saMinPt_ || std::fabs(o.eta()) > maxEta_ || o.charge() * tag.charge() >= 0)
        continue;
      TLorentzVector vp;
      vp.SetPtEtaPhiM(o.pt(), o.eta(), o.phi(), kMuMass);
      mass_ = (vt + vp).M();
      if (mass_ < saMassMin_ || mass_ > saMassMax_)
        continue;
      pair_pt_ = (vt + vp).Pt();
      pair_dR_ = reco::deltaR(tag.eta(), tag.phi(), o.eta(), o.phi());
      sa_pt_ = o.pt();
      sa_eta_ = o.eta();
      sa_phi_ = o.phi();
      sa_nStations_ = mu.numberOfMatchedStations();
      sa_validMuonHits_ = o.hitPattern().numberOfValidMuonHits();
      sa_isTracker_ = mu.isTrackerMuon();
      sa_isGlobal_ = mu.isGlobalMuon();
      sa_hasInner_ = mu.innerTrack().isNonnull();
      sa_pass_ = false;
      sa_inner_pt_ = sa_massInner_ = -1.f;
      sa_pixelHits_ = sa_pixelLayers_ = sa_trackerLayers_ = sa_bpixLayer1_ = -1;
      if (sa_hasInner_) {
        const reco::Track& in = *mu.innerTrack();
        const auto& hp = in.hitPattern();
        sa_pass_ = in.quality(reco::TrackBase::highPurity);
        sa_inner_pt_ = in.pt();
        TLorentzVector vi;
        vi.SetPtEtaPhiM(in.pt(), in.eta(), in.phi(), kMuMass);
        sa_massInner_ = (vt + vi).M();
        sa_pixelHits_ = hp.numberOfValidPixelHits();
        sa_pixelLayers_ = hp.pixelLayersWithMeasurement();
        sa_trackerLayers_ = hp.trackerLayersWithMeasurement();
        sa_bpixLayer1_ = hp.hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel, 1) ? 1 : 0;
      }
      sa_otherMu_dR_ = sa_track_dR_ = 99.f;
      sa_otherMu_ptRatio_ = sa_track_ptRatio_ = -1.f;
      sa_otherMu_pixelHits_ = sa_track_pixelHits_ = -1;
      sa_otherMu_highPurity_ = false;
      for (const auto& m2 : *muons) {
        if (&m2 == &mu || &m2 == &tag || m2.innerTrack().isNull())
          continue;
        const float d = reco::deltaR(o.eta(), o.phi(), m2.innerTrack()->eta(), m2.innerTrack()->phi());
        if (d < sa_otherMu_dR_) {
          sa_otherMu_dR_ = d;
          sa_otherMu_ptRatio_ = m2.innerTrack()->pt() / o.pt();
          sa_otherMu_pixelHits_ = m2.innerTrack()->hitPattern().numberOfValidPixelHits();
          sa_otherMu_highPurity_ = m2.innerTrack()->quality(reco::TrackBase::highPurity);
        }
      }
      for (const auto* coll : {pf.product(), lost.product()})
        for (const auto& c : *coll) {
          if (c.charge() == 0 || !c.hasTrackDetails() || c.pt() < 1.0)
            continue;
          const reco::Track& tk = c.pseudoTrack();
          if (!tk.quality(reco::TrackBase::highPurity))
            continue;
          const float d = reco::deltaR(o.eta(), o.phi(), tk.eta(), tk.phi());
          if (d < sa_track_dR_) {
            sa_track_dR_ = d;
            sa_track_ptRatio_ = tk.pt() / o.pt();
            sa_track_pixelHits_ = tk.hitPattern().numberOfValidPixelHits();
          }
        }
      sa_->Fill();
      ++nSA_;
    }

    // --- trk tree: track probes (packed PF candidates and lost tracks with track details)
    for (const auto* coll : {pf.product(), lost.product()}) {
      if (!fillTrk_)
        break;
      const bool isLost = (coll == lost.product());
      for (const auto& c : *coll) {
        if (c.charge() == 0 || !c.hasTrackDetails() || c.pt() < probeMinPt_ || std::fabs(c.eta()) > maxEta_ ||
            c.charge() * tag.charge() >= 0)
          continue;
        const reco::Track& tk = c.pseudoTrack();
        if (!tk.quality(reco::TrackBase::highPurity))
          continue;
        if (reco::deltaR(tk.eta(), tk.phi(), tag.innerTrack()->eta(), tag.innerTrack()->phi()) < 1e-3)
          continue;
        TLorentzVector vp;
        vp.SetPtEtaPhiM(tk.pt(), tk.eta(), tk.phi(), kMuMass);
        mass_ = (vt + vp).M();
        if (mass_ < massMin_ || mass_ > massMax_)
          continue;
        pair_pt_ = (vt + vp).Pt();
        pair_dR_ = reco::deltaR(tag.eta(), tag.phi(), tk.eta(), tk.phi());
        p_pt_ = tk.pt();
        p_eta_ = tk.eta();
        p_phi_ = tk.phi();
        p_dxy_ = tk.dxy(pv.position());
        p_dz_ = tk.dz(pv.position());
        p_pixelHits_ = tk.hitPattern().numberOfValidPixelHits();
        p_trackerLayers_ = tk.hitPattern().trackerLayersWithMeasurement();
        p_isLost_ = isLost;
        const pat::Muon* m = nullptr;  // the muon carrying this track
        for (const auto& mu : *muons)
          if (mu.innerTrack().isNonnull() && &mu != &tag &&
              reco::deltaR(tk.eta(), tk.phi(), mu.innerTrack()->eta(), mu.innerTrack()->phi()) < 0.005 &&
              std::fabs(mu.innerTrack()->pt() / tk.pt() - 1) < 0.02) {
            m = &mu;
            break;
          }
        p_isMuon_ = (m != nullptr);
        p_isTracker_ = m && m->isTrackerMuon();
        p_isGlobal_ = m && m->isGlobalMuon();
        p_isSoft_ = m && muon::isSoftMuon(*m, pv);
        p_isLoose_ = m && muon::isLooseMuon(*m);
        p_vtxProb_ = vtxProb(*tag.innerTrack(), tk, ttb);
        const Match l1m = matchL1(tk.eta(), tk.phi(), *l1, tag.eta(), tag.phi());
        p_l1dR_ = l1m.dR;
        p_l1qual_ = l1m.qual;
        p_hltdR_ = matchHLT(tk.eta(), tk.phi(), objs);
        trk_->Fill();
        ++nTrk_;
      }
    }
  }
}

void mytagAndProbeV9::endJob() {
  std::cout << "mytagAndProbeV9 summary: events " << nSeen_ << ", with a tag " << nTagEvents_ << ", sa entries " << nSA_
            << ", trk entries " << nTrk_ << std::endl;
  for (const auto& [n, c] : missing_)
    std::cout << "  MISSING product " << n << " in " << c << " events" << std::endl;
  for (const auto& [n, c] : fired_)
    std::cout << "  fired " << n << ": " << c << std::endl;
}

void mytagAndProbeV9::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("pfCands", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("bits", edm::InputTag("TriggerResults", "", "HLT"));
  desc.add<edm::InputTag>("objects", edm::InputTag("slimmedPatTrigger"));
  desc.add<edm::InputTag>("l1Muons", edm::InputTag("gmtStage2Digis", "Muon"));
  desc.add<edm::InputTag>("prescales", edm::InputTag("patTrigger"));
  desc.add<std::vector<std::string>>("tagPaths", {"HLT_IsoMu24_v"});
  desc.add<std::string>("analysisPath", "HLT_Mu0_L1DoubleMu_v");
  desc.add<double>("tagMinPt", 26.0);
  desc.add<double>("probeMinPt", 3.0);
  desc.add<double>("saMinPt", 2.0);
  desc.add<double>("maxEta", 2.4);
  desc.add<double>("massMin", 2.6);
  desc.add<double>("massMax", 3.6);
  desc.add<double>("saMassMin", 2.0);
  desc.add<double>("saMassMax", 4.5);
  desc.add<unsigned>("prescale", 1);
  desc.add<bool>("fillSA", true);
  desc.add<bool>("fillTrk", true);
  descriptions.add("mytagAndProbeV9", desc);
}

DEFINE_FWK_MODULE(mytagAndProbeV9);
