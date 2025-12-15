#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/mytagAndProbeV7.h"

mytagAndProbeV7::mytagAndProbeV7(const edm::ParameterSet& iConfig)
    : muonsToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfcands"))),
      TriggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
      primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      tree_(0) {
  usesResource("TFileService");
}

mytagAndProbeV7::~mytagAndProbeV7() {}

void mytagAndProbeV7::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;
  using namespace std;

  Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonsToken_, muons);

  Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfToken_, pfcands);

  Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(TriggerResultsToken_, triggerResults);

  Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(primaryVertices_Label, vertices);

  if (!muons.isValid() || !pfcands.isValid() || !triggerResults.isValid() || !triggerObjects.isValid() || !vertices.isValid())
    return;
  if (vertices->empty())
    return;
  const reco::Vertex& PV = vertices->front();

  Run = iEvent.id().run();
  LumiBlock = iEvent.id().luminosityBlock();
  Event = iEvent.id().event();

  const edm::TriggerNames& names = iEvent.triggerNames(*triggerResults);

  // 1. Tag Selection
  std::vector<const pat::Muon*> tagMuons;
  for (const auto& mu : *muons) {
    if (mu.pt() < 25.0 || std::abs(mu.eta()) > 2.4 || !mu.isTightMuon(PV))
      continue;

    bool isMatchedToTrigger = false;
    for (const auto& obj : *triggerObjects) {
      pat::TriggerObjectStandAlone unpackedObj = obj;
      unpackedObj.unpackPathNames(names);
      if (unpackedObj.hasPathName("HLT_IsoMu24_v*", true, true)) {
        if (reco::deltaR(mu.eta(), mu.phi(), unpackedObj.eta(), unpackedObj.phi()) < 0.1) {
          isMatchedToTrigger = true;
          break;
        }
      }
    }
    if (isMatchedToTrigger)
      tagMuons.push_back(&mu);
  }

  if (tagMuons.empty())
    return;

  // 2. Probe Selection: LOOP OVER PACKED PF CANDIDATES
  for (const auto& tag : tagMuons) {
    for (const auto& cand : *pfcands) {
      // Only Charged Candidates
      if (cand.charge() == 0)
        continue;

      // OPTIMIZATION: Only look for Probes near the Tag
      // High Pt J/Psi are boosted, so the muons are close.
      // This drastically reduces combinatorial background.
      if (reco::deltaR(tag->eta(), tag->phi(), cand.eta(), cand.phi()) > 1.2)
        continue;

      // Avoid self-matching
      if (reco::deltaR(tag->eta(), tag->phi(), cand.eta(), cand.phi()) < 0.01)
        continue;

      // Basic Kinematics for J/Psi Probe
      if (cand.pt() < 3.0 || std::abs(cand.eta()) > 2.4)
        continue;

      // Track Purity (Standard requirement)
      if (!cand.trackHighPurity())
        continue;

      // Opposite Charge
      if (tag->charge() * cand.charge() >= 0)
        continue;

      // Form J/Psi (Assume Muon Mass)
      TLorentzVector vTag, vProbe;
      vTag.SetPtEtaPhiM(tag->pt(), tag->eta(), tag->phi(), 0.105658);
      vProbe.SetPtEtaPhiM(cand.pt(), cand.eta(), cand.phi(), 0.105658);
      TLorentzVector vJPsi = vTag + vProbe;

      if (vJPsi.M() < 2.6 || vJPsi.M() > 3.5)
        continue;

      // --- CHECKING FOR MATCH ---
      // Does this Generic Track match a Soft Muon?
      bool isSoft = false;

      // Loop over muons to see if this track IS a Soft Muon
      for (const auto& mu : *muons) {
        // Match PF Candidate to Muon
        if (reco::deltaR(mu.eta(), mu.phi(), cand.eta(), cand.phi()) < 0.01) {
          // Check Soft ID
          if (mu.isSoftMuon(PV)) {
            isSoft = true;
          }
          break;
        }
      }

      // --- FILL TREE ---
      mass = vJPsi.M();
      probe_isSoftMuon = isSoft;

      tree_->Fill();
    }
  }
}

void mytagAndProbeV7::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple", "Total Efficiency Tree");

  tree_->Branch("Run", &Run, "Run/i");
  tree_->Branch("LumiBlock", &LumiBlock, "LumiBlock/s");
  tree_->Branch("mass", &mass, "mass/F");

  // probe_isSoftMuon = (Tracking * Matching * SoftID) Efficiency
  tree_->Branch("probe_passAlgo", &probe_isSoftMuon, "probe_passAlgo/O");
}

void mytagAndProbeV7::endJob() {}

void mytagAndProbeV7::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(mytagAndProbeV7);