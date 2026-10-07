#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
#include "TTree.h"

class mytagAndProbeV8 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit mytagAndProbeV8(const edm::ParameterSet&);
  ~mytagAndProbeV8() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  bool eventPassesHLT(const edm::TriggerResults&, const edm::TriggerNames&) const;
  const reco::Muon* bestMatchedMuon(const reco::Track&, const reco::MuonCollection&) const;
  bool goodProbeTrack(const reco::Track&, const reco::Vertex&) const;
  bool goodTagMuon(const reco::Muon&, const reco::Vertex&) const;

  // tokens
  edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;

  // config
  std::vector<std::string> hltPaths_;
  double tagMinPt_;
  double probeMinPt_;
  double maxEta_;
  double massMin_;
  double massMax_;
  double maxTrackMuonDR_;
  double maxProbeDxy_;
  double maxProbeDz_;
  int minTrackerLayers_;
  int minValidPixelHits_;

  // tree
  TTree* tree_;

  // event branches
  UInt_t run_;
  UInt_t lumi_;
  ULong64_t event_;
  Int_t nvtx_;
  Bool_t event_pass_hlt_;

  // tag branches
  Float_t tag_pt_;
  Float_t tag_eta_;
  Float_t tag_phi_;
  Int_t tag_charge_;
  Bool_t tag_isSoft_;
  Bool_t tag_isTracker_;
  Bool_t tag_isGlobal_;
  Int_t tag_nStations_;
  Float_t tag_segmentComp_;
  Float_t tag_dxy_;
  Float_t tag_dz_;
  Int_t tag_validPixelHits_;
  Int_t tag_trackerLayers_;

  // probe-track branches
  Float_t probe_pt_;
  Float_t probe_eta_;
  Float_t probe_phi_;
  Int_t probe_charge_;
  Bool_t probe_highPurity_;
  Float_t probe_normChi2_;
  Float_t probe_dxy_;
  Float_t probe_dz_;
  Int_t probe_validPixelHits_;
  Int_t probe_validTrackerHits_;
  Int_t probe_trackerLayers_;
  Int_t probe_pixelLayers_;
  Int_t probe_missingInnerHits_;

  // matched-muon branches
  Bool_t probe_hasMatchedMuon_;
  Float_t probe_trkMu_dR_;
  Bool_t probe_passSoft_;
  Bool_t probe_matched_isTracker_;
  Bool_t probe_matched_isGlobal_;
  Int_t probe_matched_nStations_;
  Float_t probe_matched_segmentComp_;
  Int_t probe_matched_validMuonHits_;
  Float_t probe_matched_pt_;
  Float_t probe_matched_eta_;
  Float_t probe_matched_phi_;

  // pair branches
  Float_t pair_mass_;
  Float_t pair_dr_;
  Float_t pair_pt_;
  Int_t pair_os_;
};

mytagAndProbeV8::mytagAndProbeV8(const edm::ParameterSet& cfg)
    : muonsToken_(consumes<reco::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
      tracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
      verticesToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
      triggerBitsToken_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerBits"))),
      hltPaths_(cfg.getParameter<std::vector<std::string>>("hltPaths")),
      tagMinPt_(cfg.getParameter<double>("tagMinPt")),
      probeMinPt_(cfg.getParameter<double>("probeMinPt")),
      maxEta_(cfg.getParameter<double>("maxEta")),
      massMin_(cfg.getParameter<double>("massMin")),
      massMax_(cfg.getParameter<double>("massMax")),
      maxTrackMuonDR_(cfg.getParameter<double>("maxTrackMuonDR")),
      maxProbeDxy_(cfg.getParameter<double>("maxProbeDxy")),
      maxProbeDz_(cfg.getParameter<double>("maxProbeDz")),
      minTrackerLayers_(cfg.getParameter<int>("minTrackerLayers")),
      minValidPixelHits_(cfg.getParameter<int>("minValidPixelHits")),
      tree_(nullptr) {}

bool mytagAndProbeV8::eventPassesHLT(const edm::TriggerResults& triggerBits, const edm::TriggerNames& names) const {
  for (unsigned int i = 0; i < triggerBits.size(); ++i) {
    if (!triggerBits.accept(i))
      continue;
    const std::string& name = names.triggerName(i);

    for (const auto& pat : hltPaths_) {
      // treat entries like "HLT_Mu0_L1DoubleMu_v"
      if (name.rfind(pat, 0) == 0)
        return true;
    }
  }
  return false;
}

bool mytagAndProbeV8::goodTagMuon(const reco::Muon& mu, const reco::Vertex& pv) const {
  if (mu.pt() < tagMinPt_)
    return false;
  if (std::abs(mu.eta()) > maxEta_)
    return false;
  if (!mu.innerTrack().isNonnull())
    return false;
  if (!muon::isSoftMuon(mu, pv))
    return false;
  if (std::abs(mu.innerTrack()->dxy(pv.position())) > 0.3)
    return false;
  if (std::abs(mu.innerTrack()->dz(pv.position())) > 20.0)
    return false;
  if (mu.numberOfMatchedStations() < 1)
    return false;
  return true;
}

bool mytagAndProbeV8::goodProbeTrack(const reco::Track& trk, const reco::Vertex& pv) const {
  if (trk.pt() < probeMinPt_)
    return false;
  if (std::abs(trk.eta()) > maxEta_)
    return false;
  if (!trk.quality(reco::TrackBase::highPurity))
    return false;
  if (std::abs(trk.dxy(pv.position())) > maxProbeDxy_)
    return false;
  if (std::abs(trk.dz(pv.position())) > maxProbeDz_)
    return false;
  if (trk.hitPattern().trackerLayersWithMeasurement() < minTrackerLayers_)
    return false;
  if (trk.hitPattern().numberOfValidPixelHits() < minValidPixelHits_)
    return false;
  return true;
}

const reco::Muon* mytagAndProbeV8::bestMatchedMuon(const reco::Track& trk, const reco::MuonCollection& muons) const {
  const reco::Muon* best = nullptr;
  double bestDR = 1e9;

  for (const auto& mu : muons) {
    if (!mu.innerTrack().isNonnull())
      continue;

    double dR = reco::deltaR(trk.eta(), trk.phi(), mu.innerTrack()->eta(), mu.innerTrack()->phi());
    if (dR < maxTrackMuonDR_ && dR < bestDR) {
      best = &mu;
      bestDR = dR;
    }
  }
  return best;
}

void mytagAndProbeV8::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tpTree", "AOD Jpsi tag-and-probe tree");

  tree_->Branch("run", &run_, "run/i");
  tree_->Branch("lumi", &lumi_, "lumi/i");
  tree_->Branch("event", &event_, "event/l");
  tree_->Branch("nvtx", &nvtx_, "nvtx/I");
  tree_->Branch("event_pass_hlt", &event_pass_hlt_, "event_pass_hlt/O");

  tree_->Branch("tag_pt", &tag_pt_, "tag_pt/F");
  tree_->Branch("tag_eta", &tag_eta_, "tag_eta/F");
  tree_->Branch("tag_phi", &tag_phi_, "tag_phi/F");
  tree_->Branch("tag_charge", &tag_charge_, "tag_charge/I");
  tree_->Branch("tag_isSoft", &tag_isSoft_, "tag_isSoft/O");
  tree_->Branch("tag_isTracker", &tag_isTracker_, "tag_isTracker/O");
  tree_->Branch("tag_isGlobal", &tag_isGlobal_, "tag_isGlobal/O");
  tree_->Branch("tag_nStations", &tag_nStations_, "tag_nStations/I");
  tree_->Branch("tag_segmentComp", &tag_segmentComp_, "tag_segmentComp/F");
  tree_->Branch("tag_dxy", &tag_dxy_, "tag_dxy/F");
  tree_->Branch("tag_dz", &tag_dz_, "tag_dz/F");
  tree_->Branch("tag_validPixelHits", &tag_validPixelHits_, "tag_validPixelHits/I");
  tree_->Branch("tag_trackerLayers", &tag_trackerLayers_, "tag_trackerLayers/I");

  tree_->Branch("probe_pt", &probe_pt_, "probe_pt/F");
  tree_->Branch("probe_eta", &probe_eta_, "probe_eta/F");
  tree_->Branch("probe_phi", &probe_phi_, "probe_phi/F");
  tree_->Branch("probe_charge", &probe_charge_, "probe_charge/I");
  tree_->Branch("probe_highPurity", &probe_highPurity_, "probe_highPurity/O");
  tree_->Branch("probe_normChi2", &probe_normChi2_, "probe_normChi2/F");
  tree_->Branch("probe_dxy", &probe_dxy_, "probe_dxy/F");
  tree_->Branch("probe_dz", &probe_dz_, "probe_dz/F");
  tree_->Branch("probe_validPixelHits", &probe_validPixelHits_, "probe_validPixelHits/I");
  tree_->Branch("probe_validTrackerHits", &probe_validTrackerHits_, "probe_validTrackerHits/I");
  tree_->Branch("probe_trackerLayers", &probe_trackerLayers_, "probe_trackerLayers/I");
  tree_->Branch("probe_pixelLayers", &probe_pixelLayers_, "probe_pixelLayers/I");
  tree_->Branch("probe_missingInnerHits", &probe_missingInnerHits_, "probe_missingInnerHits/I");

  tree_->Branch("probe_hasMatchedMuon", &probe_hasMatchedMuon_, "probe_hasMatchedMuon/O");
  tree_->Branch("probe_trkMu_dR", &probe_trkMu_dR_, "probe_trkMu_dR/F");
  tree_->Branch("probe_passSoft", &probe_passSoft_, "probe_passSoft/O");
  tree_->Branch("probe_matched_isTracker", &probe_matched_isTracker_, "probe_matched_isTracker/O");
  tree_->Branch("probe_matched_isGlobal", &probe_matched_isGlobal_, "probe_matched_isGlobal/O");
  tree_->Branch("probe_matched_nStations", &probe_matched_nStations_, "probe_matched_nStations/I");
  tree_->Branch("probe_matched_segmentComp", &probe_matched_segmentComp_, "probe_matched_segmentComp/F");
  tree_->Branch("probe_matched_validMuonHits", &probe_matched_validMuonHits_, "probe_matched_validMuonHits/I");
  tree_->Branch("probe_matched_pt", &probe_matched_pt_, "probe_matched_pt/F");
  tree_->Branch("probe_matched_eta", &probe_matched_eta_, "probe_matched_eta/F");
  tree_->Branch("probe_matched_phi", &probe_matched_phi_, "probe_matched_phi/F");

  tree_->Branch("pair_mass", &pair_mass_, "pair_mass/F");
  tree_->Branch("pair_dr", &pair_dr_, "pair_dr/F");
  tree_->Branch("pair_pt", &pair_pt_, "pair_pt/F");
  tree_->Branch("pair_os", &pair_os_, "pair_os/I");
}

void mytagAndProbeV8::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<reco::MuonCollection> muons;
  edm::Handle<reco::TrackCollection> tracks;
  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<edm::TriggerResults> triggerBits;

  iEvent.getByToken(muonsToken_, muons);
  iEvent.getByToken(tracksToken_, tracks);
  iEvent.getByToken(verticesToken_, vertices);
  iEvent.getByToken(triggerBitsToken_, triggerBits);

  if (!muons.isValid() || !tracks.isValid() || !vertices.isValid() || !triggerBits.isValid())
    return;
  if (vertices->empty())
    return;

  const reco::Vertex& pv = vertices->front();
  if (pv.isFake())
    return;

  run_ = iEvent.id().run();
  lumi_ = iEvent.id().luminosityBlock();
  event_ = iEvent.id().event();
  nvtx_ = static_cast<int>(vertices->size());

  const edm::TriggerNames& names = iEvent.triggerNames(*triggerBits);
  event_pass_hlt_ = eventPassesHLT(*triggerBits, names);
  if (!event_pass_hlt_)
    return;

  std::vector<const reco::Muon*> tags;
  tags.reserve(muons->size());
  for (const auto& mu : *muons) {
    if (goodTagMuon(mu, pv))
      tags.push_back(&mu);
  }
  if (tags.empty())
    return;

  constexpr double muMass = 0.105658;

  for (const auto* tag : tags) {
    TLorentzVector p4tag;
    p4tag.SetPtEtaPhiM(tag->pt(), tag->eta(), tag->phi(), muMass);

    for (const auto& trk : *tracks) {
      if (!goodProbeTrack(trk, pv))
        continue;

      // avoid using the tag track itself as the probe
      if (tag->innerTrack().isNonnull()) {
        if (reco::deltaR(trk.eta(), trk.phi(), tag->innerTrack()->eta(), tag->innerTrack()->phi()) < 1e-4) {
          continue;
        }
      }

      int qprod = tag->charge() * trk.charge();
      if (qprod >= 0)
        continue;

      TLorentzVector p4probe;
      p4probe.SetPtEtaPhiM(trk.pt(), trk.eta(), trk.phi(), muMass);
      TLorentzVector p4pair = p4tag + p4probe;

      const float mass = p4pair.M();
      if (mass < massMin_ || mass > massMax_)
        continue;

      const reco::Muon* matchedMu = bestMatchedMuon(trk, *muons);

      // tag
      tag_pt_ = tag->pt();
      tag_eta_ = tag->eta();
      tag_phi_ = tag->phi();
      tag_charge_ = tag->charge();
      tag_isSoft_ = muon::isSoftMuon(*tag, pv);
      tag_isTracker_ = tag->isTrackerMuon();
      tag_isGlobal_ = tag->isGlobalMuon();
      tag_nStations_ = tag->numberOfMatchedStations();
      tag_segmentComp_ = muon::segmentCompatibility(*tag);
      tag_dxy_ = tag->innerTrack().isNonnull() ? tag->innerTrack()->dxy(pv.position()) : 999.f;
      tag_dz_ = tag->innerTrack().isNonnull() ? tag->innerTrack()->dz(pv.position()) : 999.f;
      tag_validPixelHits_ = tag->innerTrack().isNonnull() ? tag->innerTrack()->hitPattern().numberOfValidPixelHits() : -1;
      tag_trackerLayers_ = tag->innerTrack().isNonnull() ? tag->innerTrack()->hitPattern().trackerLayersWithMeasurement() : -1;

      // probe track
      probe_pt_ = trk.pt();
      probe_eta_ = trk.eta();
      probe_phi_ = trk.phi();
      probe_charge_ = trk.charge();
      probe_highPurity_ = trk.quality(reco::TrackBase::highPurity);
      probe_normChi2_ = trk.normalizedChi2();
      probe_dxy_ = trk.dxy(pv.position());
      probe_dz_ = trk.dz(pv.position());
      probe_validPixelHits_ = trk.hitPattern().numberOfValidPixelHits();
      probe_validTrackerHits_ = trk.hitPattern().numberOfValidTrackerHits();
      probe_trackerLayers_ = trk.hitPattern().trackerLayersWithMeasurement();
      probe_pixelLayers_ = trk.hitPattern().pixelLayersWithMeasurement();
      probe_missingInnerHits_ = trk.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);

      // matched reco muon
      probe_hasMatchedMuon_ = (matchedMu != nullptr);
      probe_trkMu_dR_ = 999.f;
      probe_passSoft_ = false;
      probe_matched_isTracker_ = false;
      probe_matched_isGlobal_ = false;
      probe_matched_nStations_ = -1;
      probe_matched_segmentComp_ = -99.f;
      probe_matched_validMuonHits_ = -1;
      probe_matched_pt_ = -99.f;
      probe_matched_eta_ = -99.f;
      probe_matched_phi_ = -99.f;

      if (matchedMu != nullptr && matchedMu->innerTrack().isNonnull()) {
        probe_trkMu_dR_ = reco::deltaR(trk.eta(), trk.phi(), matchedMu->innerTrack()->eta(), matchedMu->innerTrack()->phi());
        probe_passSoft_ = muon::isSoftMuon(*matchedMu, pv);
        probe_matched_isTracker_ = matchedMu->isTrackerMuon();
        probe_matched_isGlobal_ = matchedMu->isGlobalMuon();
        probe_matched_nStations_ = matchedMu->numberOfMatchedStations();
        probe_matched_segmentComp_ = muon::segmentCompatibility(*matchedMu);
        probe_matched_pt_ = matchedMu->pt();
        probe_matched_eta_ = matchedMu->eta();
        probe_matched_phi_ = matchedMu->phi();

        if (matchedMu->globalTrack().isNonnull()) {
          probe_matched_validMuonHits_ = matchedMu->globalTrack()->hitPattern().numberOfValidMuonHits();
        }
      }

      // pair
      pair_mass_ = mass;
      pair_dr_ = reco::deltaR(tag->eta(), tag->phi(), trk.eta(), trk.phi());
      pair_pt_ = p4pair.Pt();
      pair_os_ = (qprod < 0) ? 1 : 0;

      tree_->Fill();
    }
  }
}

void mytagAndProbeV8::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("triggerBits", edm::InputTag("TriggerResults", "", "HLT"));

  desc.add<std::vector<std::string>>("hltPaths", {"HLT_Mu0_L1DoubleMu_v"});
  desc.add<double>("tagMinPt", 7.0);
  desc.add<double>("probeMinPt", 3.0);
  desc.add<double>("maxEta", 2.4);
  desc.add<double>("massMin", 2.6);
  desc.add<double>("massMax", 3.5);
  desc.add<double>("maxTrackMuonDR", 0.03);
  desc.add<double>("maxProbeDxy", 0.3);
  desc.add<double>("maxProbeDz", 20.0);
  desc.add<int>("minTrackerLayers", 6);
  desc.add<int>("minValidPixelHits", 1);

  descriptions.add("mytagAndProbeV8", desc);
}

DEFINE_FWK_MODULE(mytagAndProbeV8);