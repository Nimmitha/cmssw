#ifndef ZeeJmmAnalyzer_MiniAODeemm_miniAODeemm_h
#define ZeeJmmAnalyzer_MiniAODeemm_miniAODeemm_h

#include <string>
#include <vector>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TTree.h"

class miniAODeemm : public edm::EDAnalyzer {
public:
  explicit miniAODeemm(const edm::ParameterSet&);
  ~miniAODeemm() override = default;

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  template <typename T>
  void branch(const char* name, std::vector<T>& v) { tree_->Branch(name, &v); }

  void clearVectors();

  edm::EDGetTokenT<edm::View<pat::Electron>> electronToken_;
  edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> packedCandToken_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken_;

  std::string ElectronTriggerString_;

  bool requireTrigger_;
  bool requireTriggerMatch_;
  bool keepEmptyEvents_;
  bool requireLooseElectronID_;
  bool isMC_;

  TTree* tree_ = nullptr;
  unsigned int nB = 0;

  // Event-level branches.
  std::vector<unsigned int> Run, LumiBlock, nPV;
  std::vector<unsigned long long> Event;
  std::vector<int> passEleTrigger, passEleTriggerMatch;

  // Four-lepton system.
  std::vector<float> fourL_mass, fourL_px, fourL_py, fourL_pz, fourL_pt, fourL_eta, fourL_phi, fourL_rapidity;
  std::vector<float> fourL_vtxProb, fourL_pvx, fourL_pvy, fourL_pvz, fourL_pvxError, fourL_pvyError, fourL_pvzError;

  // Zee candidate.
  std::vector<float> Z_mass, Z_px, Z_py, Z_pz, Z_pt, Z_eta, Z_phi, Z_rapidity, Z_vtxProb, Z_dR_ee;
  std::vector<float> Z_trackIso03, Z_trackIso04, Z_relIso03, Z_relIso04;

  // J/psi candidate.
  std::vector<float> Jpsi_mass, Jpsi_px, Jpsi_py, Jpsi_pz, Jpsi_pt, Jpsi_eta, Jpsi_phi, Jpsi_rapidity, Jpsi_vtxProb, Jpsi_dR_mumu;
  std::vector<float> Jpsi_trackIso03, Jpsi_trackIso04, Jpsi_relIso03, Jpsi_relIso04;

  // Z-J/psi topology and angular variables.
  std::vector<float> Z_Jpsi_dR, Z_Jpsi_dPhi, Z_Jpsi_dEta, Z_Jpsi_dY, pt_balance;
  std::vector<float> cosTheta_Z_ePlus, cosTheta_Jpsi_muPlus, phi_decayPlane_Z_Jpsi;

  // Electron branches: e1/e2 are pT-ordered, with e1 leading in pT.
  std::vector<float> e1_px, e1_py, e1_pz, e1_pt, e1_eta, e1_phi;
  std::vector<float> e2_px, e2_py, e2_pz, e2_pt, e2_eta, e2_phi;
  std::vector<int> e1_charge, e2_charge;
  std::vector<float> e1_dxy, e1_dz, e2_dxy, e2_dz;
  std::vector<float> e1_mvaRaw, e2_mvaRaw;
  std::vector<int> e1_passLooseID, e2_passLooseID, e1_passWP90, e2_passWP90, e1_passWP80, e2_passWP80;
  std::vector<int> e1_triggerMatched, e2_triggerMatched;
  std::vector<float> e1_triggerDR, e2_triggerDR;

  // Muon branches: mu1/mu2 are pT-ordered, with mu1 leading in pT.
  std::vector<float> mu1_px, mu1_py, mu1_pz, mu1_pt, mu1_eta, mu1_phi;
  std::vector<float> mu2_px, mu2_py, mu2_pz, mu2_pt, mu2_eta, mu2_phi;
  std::vector<int> mu1_charge, mu2_charge;
  std::vector<int> mu1_soft, mu2_soft, mu1_loose, mu2_loose, mu1_tight, mu2_tight;
  std::vector<float> mu1_trackAbsIso03, mu2_trackAbsIso03, mu1_trackRelIso03, mu2_trackRelIso03;
  std::vector<float> mu1_pfRelIso03, mu2_pfRelIso03, mu1_pfAbsIso03, mu2_pfAbsIso03;
  std::vector<float> mu1_dxy, mu2_dxy, mu1_dz, mu2_dz, mu1_dB3D, mu2_dB3D;
  std::vector<float> mu1_normChi2, mu2_normChi2;
  std::vector<int> mu1_nValidHits, mu2_nValidHits, mu1_nValidPixelHits, mu2_nValidPixelHits;

  // Cleanliness counters.
  std::vector<int> nExtraLooseElectrons, nExtraLooseMuons;
};

#endif
