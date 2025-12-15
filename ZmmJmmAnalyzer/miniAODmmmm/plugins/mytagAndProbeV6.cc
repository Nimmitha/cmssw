#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/mytagAndProbeV6.h"

mytagAndProbeV6::mytagAndProbeV6(const edm::ParameterSet& iConfig)
    : muonsToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      TriggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
      primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      tree_(0) {
  usesResource("TFileService");
}

mytagAndProbeV6::~mytagAndProbeV6() {}

void mytagAndProbeV6::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;
  using namespace std;

  Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonsToken_, muons);

  Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(TriggerResultsToken_, triggerResults);

  Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(primaryVertices_Label, vertices);

  if (!muons.isValid() || !triggerResults.isValid() || !triggerObjects.isValid() || !vertices.isValid())
    return;
  if (vertices->empty())
    return;
  const reco::Vertex& PV = vertices->front();

  Run = iEvent.id().run();
  LumiBlock = iEvent.id().luminosityBlock();
  Event = iEvent.id().event();

  const edm::TriggerNames& names = iEvent.triggerNames(*triggerResults);

  // 1. Tag Selection (High pT, Tight, Triggered)
  std::vector<const pat::Muon*> tagMuons;
  for (const auto& mu : *muons) {
    if (mu.pt() < 25.0 || std::abs(mu.eta()) > 2.4 || !mu.isTightMuon(PV))
      continue;

    // Match to HLT_IsoMu24
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

  // 2. Probe Selection (RELAXED to TrackerMuon)
  for (const auto& tag : tagMuons) {
    for (const auto& probe : *muons) {
      if (&probe == tag)
        continue;

      // --- THE CRITICAL CHANGE ---
      // Denominator: Must be a TrackerMuon (has inner track)
      // We do NOT require SoftID here.
      if (!probe.isTrackerMuon())
        continue;

      // Standard J/Psi Kinematics
      if (probe.pt() < 3.0 || std::abs(probe.eta()) > 2.4)
        continue;

      // Opposite Charge
      if (tag->charge() * probe.charge() >= 0)
        continue;

      // Mass Window
      TLorentzVector vTag, vProbe;
      vTag.SetPtEtaPhiM(tag->pt(), tag->eta(), tag->phi(), 0.105658);
      vProbe.SetPtEtaPhiM(probe.pt(), probe.eta(), probe.phi(), 0.105658);
      TLorentzVector vJPsi = vTag + vProbe;

      if (vJPsi.M() < 2.6 || vJPsi.M() > 3.5)
        continue;

      // --- FILL TREE ---
      mass = vJPsi.M();

      // NUMERATOR: Does this TrackerMuon PASS the SoftID?
      // If aging is happening, this boolean will become false more often over time.
      probe_passAlgo = probe.isSoftMuon(PV);

      tree_->Fill();
    }
  }
}

void mytagAndProbeV6::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple", "Tracking Efficiency Tree");

  tree_->Branch("Run", &Run, "Run/i");
  tree_->Branch("LumiBlock", &LumiBlock, "LumiBlock/s");
  tree_->Branch("mass", &mass, "mass/F");
  tree_->Branch("probe_passAlgo", &probe_passAlgo, "probe_passAlgo/O");
}

void mytagAndProbeV6::endJob() {}

void mytagAndProbeV6::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(mytagAndProbeV6);