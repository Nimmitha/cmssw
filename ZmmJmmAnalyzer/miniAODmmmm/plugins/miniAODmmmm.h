#ifndef ZMMJMM_MINIAODMMMM_H
#define ZMMJMM_MINIAODMMMM_H

#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TTree.h"

class miniAODmmmm : public edm::EDAnalyzer {
public:
  explicit miniAODmmmm(const edm::ParameterSet&);
  ~miniAODmmmm() override = default;

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> packedCandToken_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;

  std::string muonTriggerString_;
  bool requireTrigger_;
  bool keepEmptyEvents_;
  bool applyBroadTopologySkim_;
  double lowMassMin_;
  double lowMassMax_;
  double zMassMin_;
  double zMassMax_;
  double broadDimuonVtxProbMin_;
  double broadFourMuVtxProbMin_;

  TTree* tree_ = nullptr;

  UInt_t nB = 0;

  std::vector<int>* Run = nullptr;
  std::vector<int>* LumiBlock = nullptr;
  std::vector<unsigned long long>* Event = nullptr;
  std::vector<unsigned int>* nPV = nullptr;

  std::vector<float>* fourMu_mass = nullptr;
  std::vector<float>* fourMu_px = nullptr;
  std::vector<float>* fourMu_py = nullptr;
  std::vector<float>* fourMu_pz = nullptr;
  std::vector<float>* fourMu_pt = nullptr;
  std::vector<float>* fourMu_eta = nullptr;
  std::vector<float>* fourMu_phi = nullptr;
  std::vector<float>* fourMu_rapidity = nullptr;
  std::vector<float>* fourMu_vtxProb = nullptr;
  std::vector<float>* fourMu_pvx = nullptr;
  std::vector<float>* fourMu_pvy = nullptr;
  std::vector<float>* fourMu_pvz = nullptr;
  std::vector<float>* fourMu_pvxError = nullptr;
  std::vector<float>* fourMu_pvyError = nullptr;
  std::vector<float>* fourMu_pvzError = nullptr;

  std::vector<bool>* passMuonTrigger = nullptr;

  std::vector<float>* pair12_mass = nullptr;
  std::vector<float>* pair12_px = nullptr;
  std::vector<float>* pair12_py = nullptr;
  std::vector<float>* pair12_pz = nullptr;
  std::vector<float>* pair12_pt = nullptr;
  std::vector<float>* pair12_eta = nullptr;
  std::vector<float>* pair12_phi = nullptr;
  std::vector<float>* pair12_rapidity = nullptr;
  std::vector<float>* pair34_mass = nullptr;
  std::vector<float>* pair34_px = nullptr;
  std::vector<float>* pair34_py = nullptr;
  std::vector<float>* pair34_pz = nullptr;
  std::vector<float>* pair34_pt = nullptr;
  std::vector<float>* pair34_eta = nullptr;
  std::vector<float>* pair34_phi = nullptr;
  std::vector<float>* pair34_rapidity = nullptr;
  std::vector<float>* pair23_mass = nullptr;
  std::vector<float>* pair23_px = nullptr;
  std::vector<float>* pair23_py = nullptr;
  std::vector<float>* pair23_pz = nullptr;
  std::vector<float>* pair23_pt = nullptr;
  std::vector<float>* pair23_eta = nullptr;
  std::vector<float>* pair23_phi = nullptr;
  std::vector<float>* pair23_rapidity = nullptr;
  std::vector<float>* pair14_mass = nullptr;
  std::vector<float>* pair14_px = nullptr;
  std::vector<float>* pair14_py = nullptr;
  std::vector<float>* pair14_pz = nullptr;
  std::vector<float>* pair14_pt = nullptr;
  std::vector<float>* pair14_eta = nullptr;
  std::vector<float>* pair14_phi = nullptr;
  std::vector<float>* pair14_rapidity = nullptr;

  std::vector<float>* pair12_vtxProb = nullptr;
  std::vector<float>* pair34_vtxProb = nullptr;
  std::vector<float>* pair23_vtxProb = nullptr;
  std::vector<float>* pair14_vtxProb = nullptr;

  std::vector<float>* pair12_dR = nullptr;
  std::vector<float>* pair34_dR = nullptr;
  std::vector<float>* pair23_dR = nullptr;
  std::vector<float>* pair14_dR = nullptr;
  std::vector<float>* pair12_34_dR = nullptr;
  std::vector<float>* pair23_14_dR = nullptr;
  std::vector<float>* pair12_34_dPhi = nullptr;
  std::vector<float>* pair23_14_dPhi = nullptr;
  std::vector<float>* pair12_34_dEta = nullptr;
  std::vector<float>* pair23_14_dEta = nullptr;
  std::vector<float>* pair12_34_dY = nullptr;
  std::vector<float>* pair23_14_dY = nullptr;

  std::vector<float>* pair12_trackIso03 = nullptr;
  std::vector<float>* pair12_trackIso04 = nullptr;
  std::vector<float>* pair12_relIso03 = nullptr;
  std::vector<float>* pair12_relIso04 = nullptr;
  std::vector<float>* pair34_trackIso03 = nullptr;
  std::vector<float>* pair34_trackIso04 = nullptr;
  std::vector<float>* pair34_relIso03 = nullptr;
  std::vector<float>* pair34_relIso04 = nullptr;
  std::vector<float>* pair23_trackIso03 = nullptr;
  std::vector<float>* pair23_trackIso04 = nullptr;
  std::vector<float>* pair23_relIso03 = nullptr;
  std::vector<float>* pair23_relIso04 = nullptr;
  std::vector<float>* pair14_trackIso03 = nullptr;
  std::vector<float>* pair14_trackIso04 = nullptr;
  std::vector<float>* pair14_relIso03 = nullptr;
  std::vector<float>* pair14_relIso04 = nullptr;

  std::vector<float>* pair12_cosThetaMu = nullptr;
  std::vector<float>* pair34_cosThetaMu = nullptr;
  std::vector<float>* pair23_cosThetaMu = nullptr;
  std::vector<float>* pair14_cosThetaMu = nullptr;
  std::vector<float>* pair12_34_phiDecayPlane = nullptr;
  std::vector<float>* pair23_14_phiDecayPlane = nullptr;

  std::vector<float>* muP1_px = nullptr;
  std::vector<float>* muP1_py = nullptr;
  std::vector<float>* muP1_pz = nullptr;
  std::vector<float>* muP1_pt = nullptr;
  std::vector<float>* muP1_eta = nullptr;
  std::vector<float>* muP1_phi = nullptr;
  std::vector<int>* muP1_charge = nullptr;
  std::vector<bool>* muP1_soft = nullptr;
  std::vector<bool>* muP1_tight = nullptr;
  std::vector<bool>* muP1_loose = nullptr;
  std::vector<float>* muP1_isoTrack = nullptr;
  std::vector<float>* muP1_isoHcal = nullptr;
  std::vector<float>* muP1_isoEcal = nullptr;
  std::vector<float>* muP1_isoCalo = nullptr;
  std::vector<float>* muP1_pfAbsIso03 = nullptr;
  std::vector<float>* muP1_pfAbsIso04 = nullptr;
  std::vector<float>* muP1_pfRelIso03 = nullptr;
  std::vector<float>* muP1_pfRelIso04 = nullptr;
  std::vector<float>* muP1_dB3D = nullptr;

  std::vector<float>* muM1_px = nullptr;
  std::vector<float>* muM1_py = nullptr;
  std::vector<float>* muM1_pz = nullptr;
  std::vector<float>* muM1_pt = nullptr;
  std::vector<float>* muM1_eta = nullptr;
  std::vector<float>* muM1_phi = nullptr;
  std::vector<int>* muM1_charge = nullptr;
  std::vector<bool>* muM1_soft = nullptr;
  std::vector<bool>* muM1_tight = nullptr;
  std::vector<bool>* muM1_loose = nullptr;
  std::vector<float>* muM1_isoTrack = nullptr;
  std::vector<float>* muM1_isoHcal = nullptr;
  std::vector<float>* muM1_isoEcal = nullptr;
  std::vector<float>* muM1_isoCalo = nullptr;
  std::vector<float>* muM1_pfAbsIso03 = nullptr;
  std::vector<float>* muM1_pfAbsIso04 = nullptr;
  std::vector<float>* muM1_pfRelIso03 = nullptr;
  std::vector<float>* muM1_pfRelIso04 = nullptr;
  std::vector<float>* muM1_dB3D = nullptr;

  std::vector<float>* muP2_px = nullptr;
  std::vector<float>* muP2_py = nullptr;
  std::vector<float>* muP2_pz = nullptr;
  std::vector<float>* muP2_pt = nullptr;
  std::vector<float>* muP2_eta = nullptr;
  std::vector<float>* muP2_phi = nullptr;
  std::vector<int>* muP2_charge = nullptr;
  std::vector<bool>* muP2_soft = nullptr;
  std::vector<bool>* muP2_tight = nullptr;
  std::vector<bool>* muP2_loose = nullptr;
  std::vector<float>* muP2_isoTrack = nullptr;
  std::vector<float>* muP2_isoHcal = nullptr;
  std::vector<float>* muP2_isoEcal = nullptr;
  std::vector<float>* muP2_isoCalo = nullptr;
  std::vector<float>* muP2_pfAbsIso03 = nullptr;
  std::vector<float>* muP2_pfAbsIso04 = nullptr;
  std::vector<float>* muP2_pfRelIso03 = nullptr;
  std::vector<float>* muP2_pfRelIso04 = nullptr;
  std::vector<float>* muP2_dB3D = nullptr;

  std::vector<float>* muM2_px = nullptr;
  std::vector<float>* muM2_py = nullptr;
  std::vector<float>* muM2_pz = nullptr;
  std::vector<float>* muM2_pt = nullptr;
  std::vector<float>* muM2_eta = nullptr;
  std::vector<float>* muM2_phi = nullptr;
  std::vector<int>* muM2_charge = nullptr;
  std::vector<bool>* muM2_soft = nullptr;
  std::vector<bool>* muM2_tight = nullptr;
  std::vector<bool>* muM2_loose = nullptr;
  std::vector<float>* muM2_isoTrack = nullptr;
  std::vector<float>* muM2_isoHcal = nullptr;
  std::vector<float>* muM2_isoEcal = nullptr;
  std::vector<float>* muM2_isoCalo = nullptr;
  std::vector<float>* muM2_pfAbsIso03 = nullptr;
  std::vector<float>* muM2_pfAbsIso04 = nullptr;
  std::vector<float>* muM2_pfRelIso03 = nullptr;
  std::vector<float>* muM2_pfRelIso04 = nullptr;
  std::vector<float>* muM2_dB3D = nullptr;

  std::vector<float>* muP1_dxy = nullptr;
  std::vector<float>* muM1_dxy = nullptr;
  std::vector<float>* muP1_dz = nullptr;
  std::vector<float>* muM1_dz = nullptr;
  std::vector<float>* muP2_dxy = nullptr;
  std::vector<float>* muM2_dxy = nullptr;
  std::vector<float>* muP2_dz = nullptr;
  std::vector<float>* muM2_dz = nullptr;

  std::vector<float>* muM1_normChi2 = nullptr;
  std::vector<int>* muM1_nValidHits = nullptr;
  std::vector<int>* muM1_nValidPixelHits = nullptr;
  std::vector<float>* muP1_normChi2 = nullptr;
  std::vector<int>* muP1_nValidHits = nullptr;
  std::vector<int>* muP1_nValidPixelHits = nullptr;
  std::vector<float>* muM2_normChi2 = nullptr;
  std::vector<int>* muM2_nValidHits = nullptr;
  std::vector<int>* muM2_nValidPixelHits = nullptr;
  std::vector<float>* muP2_normChi2 = nullptr;
  std::vector<int>* muP2_nValidHits = nullptr;
  std::vector<int>* muP2_nValidPixelHits = nullptr;

  std::vector<float>* muP1_p4_pt = nullptr;
  std::vector<float>* muP1_p4_eta = nullptr;
  std::vector<float>* muP1_p4_phi = nullptr;
  std::vector<float>* muP1_p4_px = nullptr;
  std::vector<float>* muP1_p4_py = nullptr;
  std::vector<float>* muP1_p4_pz = nullptr;
  std::vector<float>* muM1_p4_pt = nullptr;
  std::vector<float>* muM1_p4_eta = nullptr;
  std::vector<float>* muM1_p4_phi = nullptr;
  std::vector<float>* muM1_p4_px = nullptr;
  std::vector<float>* muM1_p4_py = nullptr;
  std::vector<float>* muM1_p4_pz = nullptr;
  std::vector<float>* muP2_p4_pt = nullptr;
  std::vector<float>* muP2_p4_eta = nullptr;
  std::vector<float>* muP2_p4_phi = nullptr;
  std::vector<float>* muP2_p4_px = nullptr;
  std::vector<float>* muP2_p4_py = nullptr;
  std::vector<float>* muP2_p4_pz = nullptr;
  std::vector<float>* muM2_p4_pt = nullptr;
  std::vector<float>* muM2_p4_eta = nullptr;
  std::vector<float>* muM2_p4_phi = nullptr;
  std::vector<float>* muM2_p4_px = nullptr;
  std::vector<float>* muM2_p4_py = nullptr;
  std::vector<float>* muM2_p4_pz = nullptr;
};

#endif
