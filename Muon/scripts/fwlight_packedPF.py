from DataFormats.FWLite import Handle, Events
import mplhep as hep
import ROOT
import matplotlib.pyplot as plt
import pandas as pd
# choose backend for matplotlib
plt.switch_backend('agg')
plt.style.use(hep.style.CMS)

ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.FWLiteEnabler.enable()


events = Events("/eos/home-n/nkarunar/2023Analysis/datasets/2018_MC/Zuu_2018_MINIAOD.root")

genParticleLabel = "prunedGenParticles"
genParticles = Handle("std::vector<reco::GenParticle>")

tracksLabel = "packedPFCandidates"
tracks = Handle("std::vector<pat::PackedCandidate>")

# recMuonLabel = "slimmedMuons"
# recMuons = Handle("std::vector<pat::Muon>")
# genMuons = Handle("std::vector<pat::PackedGenParticle>")
# l1MuonLabel = "gmtStage2Digis:Muon"
# l1Muons = Handle("BXVector<l1t::Muon>")
# triggerObjects, triggerObjectLabel = Handle("std::vector<pat::TriggerObjectStandAlone>"), "slimmedPatTrigger"
# triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")

Event = []
track_pt = []
track_eta = []


events.toBegin()
for iev,event in enumerate(events):

    event.getByLabel(tracksLabel, tracks)

    for j, track in enumerate(tracks.product()):
        if track.charge() == 0:
            continue

        if track.pt() < 0.5:
            continue

        if abs(track.pdgId()) != 211:
            continue

        if abs(track.eta()) > 2.4:
            continue

        if track.trackHighPurity() == False:
            continue


        # print track information
        # if a phi track is found, print it
        print("Track %d: pt %f, eta %f, phi %f, pdgId %d, charge %d" % (
        j, track.pt(), track.eta(), track.phi(), track.pdgId(), track.charge()))
    
        # print track vertex information
        vertex = track.vertex()
        print("Track %d: vx %f, vy %f, vz %f" % (
            j, vertex.x(), vertex.y(), vertex.z()))
        
        # print track hit pattern
        print(track.trackHighPurity())
