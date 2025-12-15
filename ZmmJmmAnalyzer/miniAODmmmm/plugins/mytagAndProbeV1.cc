#include "ZmmJmmAnalyzer/miniAODmmmm/plugins/mytagAndProbeV1.h"

mytagAndProbeV1::mytagAndProbeV1(const edm::ParameterSet& iConfig)
    : muonsToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      TriggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
      triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
      primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      isMC_(iConfig.getParameter<bool>("isMC")),
      tree_(0) {
  usesResource("TFileService");
}

mytagAndProbeV1::~mytagAndProbeV1() {}

void mytagAndProbeV1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

  // Basic validity checks
  if (!muons.isValid() || !triggerResults.isValid() || !triggerObjects.isValid() || !vertices.isValid())
    return;

  // FIX 1: Require at least one good primary vertex and store it
  if (vertices->empty())
    return;
  const reco::Vertex& PV = vertices->front();

  Run = iEvent.id().run();
  LumiBlock = iEvent.id().luminosityBlock();
  Event = iEvent.id().event();

  const edm::TriggerNames& names = iEvent.triggerNames(*triggerResults);

  // 1. Check if Event passed the Single Muon Trigger (Reference Trigger)
  bool eventPassesIsoMu24 = false;
  for (unsigned int i = 0; i < triggerResults->size(); ++i) {
    std::string name = names.triggerName(i);
    if (name.find("HLT_IsoMu24_v") != std::string::npos && triggerResults->accept(i)) {
      eventPassesIsoMu24 = true;
      break;
    }
  }

  if (!eventPassesIsoMu24)
    return;  // Skip event if reference trigger didn't fire

  // 2. Check if Event passed the Double Muon Trigger (The Signal we want to measure)
  bool eventPassesDoubleMu = false;
  for (unsigned int i = 0; i < triggerResults->size(); ++i) {
    std::string name = names.triggerName(i);
    if (name.find("HLT_Mu0_L1DoubleMu_v") != std::string::npos && triggerResults->accept(i)) {
      eventPassesDoubleMu = true;
      break;
    }
  }

  // 3. Identify TAG muons (Matched to IsoMu24)
  std::vector<const pat::Muon*> tagMuons;

  for (const auto& mu : *muons) {
    // Tag Selection:
    // FIX 2: Passed PV to isTightMuon()
    if (mu.pt() < 25.0 || std::abs(mu.eta()) > 2.4 || !mu.isTightMuon(PV))
      continue;

    // Match to Trigger Object
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

    if (isMatchedToTrigger) {
      tagMuons.push_back(&mu);
    }
  }

  if (tagMuons.empty())
    return;

  // 4. Loop over Tags and find Probes
  for (const auto& tag : tagMuons) {
    for (const auto& probe : *muons) {
      // Don't match to itself
      if (&probe == tag)
        continue;

      // Probe Selection: J/Psi specific (Low pT allowed)
      if (probe.pt() < 3.0 || std::abs(probe.eta()) > 2.4)
        continue;

      // Opposite Charge
      if (tag->charge() * probe.charge() >= 0)
        continue;

      // Soft Muon ID (Standard for J/Psi analysis)
      // Using PV for consistency
      if (!probe.isSoftMuon(PV))
        continue;

      // Calculate Invariant Mass
      TLorentzVector vTag, vProbe;
      vTag.SetPtEtaPhiM(tag->pt(), tag->eta(), tag->phi(), 0.105658);
      vProbe.SetPtEtaPhiM(probe.pt(), probe.eta(), probe.phi(), 0.105658);
      TLorentzVector vJPsi = vTag + vProbe;

      // Wide J/Psi Window (e.g., 2.6 to 3.6)
      if (vJPsi.M() < 2.6 || vJPsi.M() > 3.6)
        continue;

      // --- FILL TREE ---

      tag_pt = tag->pt();
      tag_eta = tag->eta();
      tag_phi = tag->phi();

      probe_pt = probe.pt();
      probe_eta = probe.eta();
      probe_phi = probe.phi();

      mass = vJPsi.M();

      // This is the CRITICAL boolean.
      // If the event passed HLT_Mu0_L1DoubleMu, then this J/Psi (Tag+Probe)
      // contributed to the rate you are measuring in your main analysis.
      probe_passAlgo = eventPassesDoubleMu;

      tree_->Fill();
    }
  }
}

void mytagAndProbeV1::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple", "TagAndProbe Tree");

  tree_->Branch("Run", &Run, "Run/i");
  tree_->Branch("LumiBlock", &LumiBlock, "LumiBlock/s");
  tree_->Branch("Event", &Event, "Event/l");

  tree_->Branch("tag_pt", &tag_pt, "tag_pt/F");
  tree_->Branch("tag_eta", &tag_eta, "tag_eta/F");
  tree_->Branch("tag_phi", &tag_phi, "tag_phi/F");

  tree_->Branch("probe_pt", &probe_pt, "probe_pt/F");
  tree_->Branch("probe_eta", &probe_eta, "probe_eta/F");
  tree_->Branch("probe_phi", &probe_phi, "probe_phi/F");

  tree_->Branch("mass", &mass, "mass/F");

  // The boolean flag to calculate efficiency: N(pass) / N(total)
  tree_->Branch("probe_passAlgo", &probe_passAlgo, "probe_passAlgo/O");
}

void mytagAndProbeV1::endJob() {}

void mytagAndProbeV1::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// FIX 3: Ensure module name matches class name
DEFINE_FWK_MODULE(mytagAndProbeV1);