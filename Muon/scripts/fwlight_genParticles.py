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

# recMuonLabel = "slimmedMuons"
# recMuons = Handle("std::vector<pat::Muon>")
# genMuons = Handle("std::vector<pat::PackedGenParticle>")
# l1MuonLabel = "gmtStage2Digis:Muon"
# l1Muons = Handle("BXVector<l1t::Muon>")
# triggerObjects, triggerObjectLabel = Handle("std::vector<pat::TriggerObjectStandAlone>"), "slimmedPatTrigger"
# triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")

# diMu = []
# canvas = ROOT.TCanvas("canvas", "Histogram", 800, 600)
# histogram = ROOT.TH1F("histogram", "Histogram of Z Candidates in Di-Muon Channel", 90, 30, 120)

event_no = []

# Momentum
ar_H_Pt = []
ar_H_Pz = []
ar_Phi_Pt = []
ar_Phi_Pz = []
ar_Z_Pt = []
ar_Z_Pz = []

# lab frame angels
ar_ang_lab_Kp_Km = []
ar_ang_lab_Kp_Phi = []
ar_ang_lab_Km_Phi = []

ar_ang_lab_muP_muM = []
ar_ang_lab_muP_Z = []
ar_ang_lab_muM_Z = []

ar_ang_lab_Z_Phi = []

# Higgs frame angles
ar_ang_Hf_Kp_Km = []
ar_ang_Hf_Kp_Phi = []
ar_ang_Hf_Km_Phi = []

ar_ang_Hf_muP_muM = []
ar_ang_Hf_muP_Z = []
ar_ang_Hf_muM_Z = []

ar_ang_Hf_Z_Phi = []


# loop over events
events.toBegin()
for i, event in enumerate(events):
    # print "#%d: run %6d, lumi %4d, event %12d" % (i, event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(), event.eventAuxiliary().event())

    event.getByLabel(genParticleLabel, genParticles)
    for j, genParticle in enumerate(genParticles.product()):
        # print(j)
        # print(genParticle.pdgId(), genParticle.pt(), genParticle.eta(), genParticle.phi(), genParticle.mass())

        # if it is a muon with pt > 3 and |eta| < 2.4
        if (genParticle.pdgId() == 13) and (genParticle.pt() > 3) and (abs(genParticle.eta()) < 2.4):
            # if the mother is a Z boson
            if (abs(genParticle.mother().pdgId()) == 23):
                # continue if the mother of the mother is not a Higgs boson
                if (not abs(genParticle.mother().mother().pdgId()) == 35):
                    continue

                Higgs = genParticle.mother().mother()
                # number of daughters of the higgs
                print(Higgs.numberOfDaughters())

                # if the daughters of the higgs are Z bosons and phi
                if (Higgs.numberOfDaughters() == 2):
                    if (Higgs.daughterRefVector()[0].pdgId() == 23 and Higgs.daughterRefVector()[1].pdgId() == 333):
                        # print("Z and phi")
                        pass
                    elif (Higgs.daughterRefVector()[0].pdgId() == 333 and Higgs.daughterRefVector()[1].pdgId() == 23):
                        # print("phi and Z")
                        pass
                    else:
                        continue

                    # print the daughers of the higgs
                    for k, daughter in enumerate(Higgs.daughterRefVector()):
                        if (abs(daughter.pdgId()) == 333):
                            # print("phi")
                            Phi = daughter
                            # print the daughters of the phi
                            for l, daughter in enumerate(Phi.daughterRefVector()):
                                # print(daughter.pdgId(), daughter.mass())

                                # if the daughters of phi are K+ and K-
                                if (daughter.pdgId() == 321):
                                    # print("K+")
                                    Kp = daughter
                                elif (daughter.pdgId() == -321):
                                    # print("K-")
                                    Km = daughter

                        if (abs(daughter.pdgId()) == 23):
                            # print("Z")
                            Zboson = daughter

                            # print the daughters of the Z boson
                            for l, daughter in enumerate(Zboson.daughterRefVector()):
                                # print(daughter.pdgId(), daughter.mass())
                                if (daughter.pdgId() == 13):
                                    # print("mu-")
                                    muM = daughter
                                elif (daughter.pdgId() == -13):
                                    # print("mu+")
                                    muP = daughter
                        else:
                            # print("Not phi")
                            pass

                    # candidate_list.append([i, Phi, Kp, Km])
                    print("Event: ", i)
                    print("Z: ", Zboson.pdgId(), Zboson.pt(), Zboson.eta(), Zboson.phi(), Zboson.mass())
                    print("Phi: ", Phi.pdgId(), Phi.pt(), Phi.eta(), Phi.phi(), Phi.mass())
                    print("K+: ", Kp.pdgId(), Kp.pt(), Kp.eta(), Kp.phi(), Kp.mass())
                    print("K-: ", Km.pdgId(), Km.pt(), Km.eta(), Km.phi(), Km.mass())
                    print("mu+: ", muP.pdgId(), muP.pt(), muP.eta(), muP.phi(), muP.mass())
                    print("mu-: ", muM.pdgId(), muM.pt(), muM.eta(), muM.phi(), muM.mass())

                    # save the event number
                    event_no.append(i)

                    # Save the momentums
                    ar_H_Pt.append(Higgs.pt())
                    ar_H_Pz.append(Higgs.pz())
                    ar_Phi_Pt.append(Phi.pt())
                    ar_Phi_Pz.append(Phi.pz())
                    ar_Z_Pt.append(Zboson.pt())
                    ar_Z_Pz.append(Zboson.pz())

                    # get the Lorentz vectors of the particles
                    Kp_lorentz = ROOT.TLorentzVector(Kp.p4().px(), Kp.p4().py(), Kp.p4().pz(), Kp.p4().energy())
                    Km_lorentz = ROOT.TLorentzVector(Km.p4().px(), Km.p4().py(), Km.p4().pz(), Km.p4().energy())
                    Phi_lorentz = ROOT.TLorentzVector(Phi.p4().px(), Phi.p4().py(), Phi.p4().pz(), Phi.p4().energy())
                    Z_lorentz = ROOT.TLorentzVector(Zboson.p4().px(), Zboson.p4().py(), Zboson.p4().pz(), Zboson.p4().energy())
                    Higgs_lorentz = ROOT.TLorentzVector(Higgs.p4().px(), Higgs.p4().py(), Higgs.p4().pz(), Higgs.p4().energy())
                    muP_lorentz = ROOT.TLorentzVector(muP.p4().px(), muP.p4().py(), muP.p4().pz(), muP.p4().energy())
                    muM_lorentz = ROOT.TLorentzVector(muM.p4().px(), muM.p4().py(), muM.p4().pz(), muM.p4().energy())

                    # compute the opening angles in the lab frame
                    ang_lab_Kp_Km = ROOT.TMath.Cos(Kp_lorentz.Angle(Km_lorentz.Vect()))
                    ang_lab_Kp_Phi = ROOT.TMath.Cos(Kp_lorentz.Angle(Phi_lorentz.Vect()))
                    ang_lab_Km_Phi = ROOT.TMath.Cos(Km_lorentz.Angle(Phi_lorentz.Vect()))
                    ang_lab_muP_muM = ROOT.TMath.Cos(muP_lorentz.Angle(muM_lorentz.Vect()))
                    ang_lab_muP_Z = ROOT.TMath.Cos(muP_lorentz.Angle(Z_lorentz.Vect()))
                    ang_lab_muM_Z = ROOT.TMath.Cos(muM_lorentz.Angle(Z_lorentz.Vect()))
                    ang_lab_Z_Phi = ROOT.TMath.Cos(Z_lorentz.Angle(Phi_lorentz.Vect()))

                    # save the lab frame angles
                    ar_ang_lab_Km_Phi.append(ang_lab_Km_Phi)
                    ar_ang_lab_Kp_Km.append(ang_lab_Kp_Km)
                    ar_ang_lab_Kp_Phi.append(ang_lab_Kp_Phi)
                    ar_ang_lab_muP_muM.append(ang_lab_muP_muM)
                    ar_ang_lab_muP_Z.append(ang_lab_muP_Z)
                    ar_ang_lab_muM_Z.append(ang_lab_muM_Z)
                    ar_ang_lab_Z_Phi.append(ang_lab_Z_Phi)

                    # boost the particles to the rest frame of the Higgs boson
                    Kp_lorentz.Boost(-Higgs_lorentz.BoostVector())
                    Km_lorentz.Boost(-Higgs_lorentz.BoostVector())
                    Phi_lorentz.Boost(-Higgs_lorentz.BoostVector())
                    muP_lorentz.Boost(-Higgs_lorentz.BoostVector())
                    muM_lorentz.Boost(-Higgs_lorentz.BoostVector())
                    Z_lorentz.Boost(-Higgs_lorentz.BoostVector())

                    # compute the opening angles in the Higgs frame
                    ang_Hf_Kp_Km = ROOT.TMath.Cos(Kp_lorentz.Angle(Km_lorentz.Vect()))
                    ang_Hf_Kp_Phi = ROOT.TMath.Cos(Kp_lorentz.Angle(Phi_lorentz.Vect()))
                    ang_Hf_Km_Phi = ROOT.TMath.Cos(Km_lorentz.Angle(Phi_lorentz.Vect()))
                    ang_Hf_muP_muM = ROOT.TMath.Cos(muP_lorentz.Angle(muM_lorentz.Vect()))
                    ang_Hf_muP_Z = ROOT.TMath.Cos(muP_lorentz.Angle(Z_lorentz.Vect()))
                    ang_Hf_muM_Z = ROOT.TMath.Cos(muM_lorentz.Angle(Z_lorentz.Vect()))
                    ang_Hf_Z_Phi = ROOT.TMath.Cos(Z_lorentz.Angle(Phi_lorentz.Vect()))

                    # save the Higgs frame angles
                    ar_ang_Hf_Km_Phi.append(ang_Hf_Km_Phi)
                    ar_ang_Hf_Kp_Km.append(ang_Hf_Kp_Km)
                    ar_ang_Hf_Kp_Phi.append(ang_Hf_Kp_Phi)
                    ar_ang_Hf_muP_muM.append(ang_Hf_muP_muM)
                    ar_ang_Hf_muP_Z.append(ang_Hf_muP_Z)
                    ar_ang_Hf_muM_Z.append(ang_Hf_muM_Z)
                    ar_ang_Hf_Z_Phi.append(ang_Hf_Z_Phi)

                    # opening angle between K+ and K-
                    # opening_angle = Kp_momentum.Angle(Km_momentum)
                    # print("Opening angle: ", opening_angle)
                    # print("Opening angle in degrees: ", opening_angle * 180 / ROOT.TMath.Pi())

# make a dataframe from the arrays
df = pd.DataFrame({'event': event_no, 'H_Pt': ar_H_Pt, 'H_Pz': ar_H_Pz, 'Phi_Pt': ar_Phi_Pt, 'Phi_Pz': ar_Phi_Pz, 'Z_Pt': ar_Z_Pt, 'Z_Pz': ar_Z_Pz, 'ang_lab_Kp_Km': ar_ang_lab_Kp_Km, 'ang_lab_Kp_Phi': ar_ang_lab_Kp_Phi, 'ang_lab_Km_Phi': ar_ang_lab_Km_Phi, 'ang_lab_muP_muM': ar_ang_lab_muP_muM, 'ang_lab_muP_Z': ar_ang_lab_muP_Z, 'ang_lab_muM_Z': ar_ang_lab_muM_Z, 'ang_lab_Z_Phi': ar_ang_lab_Z_Phi, 'ang_Hf_Kp_Km': ar_ang_Hf_Kp_Km, 'ang_Hf_Kp_Phi': ar_ang_Hf_Kp_Phi, 'ang_Hf_Km_Phi': ar_ang_Hf_Km_Phi, 'ang_Hf_muP_muM': ar_ang_Hf_muP_muM, 'ang_Hf_muP_Z': ar_ang_Hf_muP_Z, 'ang_Hf_muM_Z': ar_ang_Hf_muM_Z, 'ang_Hf_Z_Phi': ar_ang_Hf_Z_Phi})

df.to_csv("genParticles.csv", index=False)

# for i, candidate in enumerate(candidate_list):
#   print("Event: ", candidate[0])
#   print(candidate[1].pdgId())
#   # print("Phi: ", candidate[1].pdgId(), candidate[1].pt(), candidate[1].eta(), candidate[1].phi(), candidate[1].mass())
#   # print("K+: ", candidate[2].pdgId(), candidate[2].pt(), candidate[2].eta(), candidate[2].phi(), candidate[2].mass())
#   # print("K-: ", candidate[3].pdgId(), candidate[3].pt(), candidate[3].eta(), candidate[3].phi(), candidate[3].mass())
#   break


# plt.hist(array_angle_Kp_Km, bins=100)
# plt.xlabel("Opening angle K+ and K- (Cosine)")
# # disable x axis offset
# plt.ticklabel_format(useOffset=False)

# plt.ylabel("Number of events")
# # plt.xlim(-1, -0.999999996)
# plt.savefig("opening_angle_Kp_Km.png")

# plt.clf()

# plt.hist(array_angle_Kp_Phi, bins=100, label="K+ and Phi")
# plt.hist(array_angle_Km_Phi, bins=100, label="K- and Phi")
# # plt.xlabel("Opening angle K- and Phi (Cosine)")
# plt.xlabel("Opening angles (Cosine)")
# plt.ylabel("Number of events")
# plt.legend()
# plt.savefig("opening_angle_K_Phi.png")

# plt.clf()

# plt.hist(array_angle_Z_Phi, bins=100)
# plt.xlabel("Opening angle Z and Phi (Cosine)")
# plt.ylabel("Number of events")
# plt.savefig("opening_angle_Z_Phi.png")
