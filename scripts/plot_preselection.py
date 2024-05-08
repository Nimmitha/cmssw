import ROOT
import uproot as up
import matplotlib.pyplot as plt
import numpy as np
import awkward as ak
import pandas as pd
import mplhep as hep
hep.style.use(hep.style.CMS)

FIGURE_PATH = 'plots/Data/pre/'

# open the file
file = up.open('../output/Data/TTree_13TeV_fourmuon_2018.root')
events = file['ntuple;1']
columns = ['Y_TriggerPath',
           'Z_soft1', 'Z_soft2',
           'Z_mass',
           'Z_pt1', 'Z_pt2', 'Z_eta1', 'Z_eta2',
           'Z_lowPt', 'Z_highPt',
           'Z_trackIso1', 'Z_trackIso2',
           'Y_mass',
           'Y_pt1', 'Y_pt2', 'Y_eta1', 'Y_eta2',
           'Y_mvaIsoWP90_1', 'Y_mvaIsoWP90_2',
           'Z_VtxProb', 'Y_VtxProb',
           'FourL_mass', 'FourL_VtxProb']
# branches = events.arrays(columns, entry_start=430, entry_stop=450)
branches = events.arrays(columns)

# convert awkward array to numpy array
Y_candidates = ak.flatten(branches['Y_mass']).to_numpy()
Z_candidates = ak.flatten(branches['Z_mass']).to_numpy()
lowPt_muon_candidates = ak.flatten(branches['Z_lowPt']).to_numpy()
highPt_muon_candidates = ak.flatten(branches['Z_highPt']).to_numpy()
Y_VtxProb = ak.flatten(branches['Y_VtxProb']).to_numpy()
Z_VtxProb = ak.flatten(branches['Z_VtxProb']).to_numpy()
FourL_mass = ak.flatten(branches['FourL_mass']).to_numpy()
FourL_VtxProb = ak.flatten(branches['FourL_VtxProb']).to_numpy()

NEvents = len(branches)
NCandidates = len(Z_candidates)

# plot the histograms
# Plot the Y mass
plt.figure()
plt.hist(Y_candidates, bins=80, range=(6, 12))

# Annotate mean and standard deviation on the plot
plt.text(0.05, 0.9, f'Mean: {np.mean(Y_candidates):.2f}', ha='left', va='center', transform=plt.gca().transAxes)
plt.text(0.05, 0.85, f'Std: {np.std(Y_candidates):.2f}', ha='left', va='center', transform=plt.gca().transAxes)
plt.text(0.05, 0.8, f'Candi/Eve: {NCandidates}/{NEvents}', ha='left', va='center', transform=plt.gca().transAxes)

plt.xlabel("Dielectron mass [GeV]")
plt.tight_layout()
plt.savefig(f"{FIGURE_PATH}/Y_mass.png")


# Plot the Z mass
plt.figure()
plt.hist(Z_candidates, bins=80)

# Annotate mean and standard deviation on the plot on the left top corner
plt.text(0.05, 0.9, f'Mean: {np.mean(Z_candidates):.2f}', ha='left', va='center', transform=plt.gca().transAxes)
plt.text(0.05, 0.85, f'Std: {np.std(Z_candidates):.2f}', ha='left', va='center', transform=plt.gca().transAxes)
plt.savefig(FIGURE_PATH + "Z_mass.png")
plt.text(0.05, 0.8, f'Candi/Eve: {NCandidates}/{NEvents}', ha='left', va='center', transform=plt.gca().transAxes)

plt.xlabel("Dimuon mass [GeV]")
plt.tight_layout()
plt.savefig(f"{FIGURE_PATH}/Z_mass.png")


# plot the lowPt and highPt muon candidates on the same plot
plt.figure()
plt.hist(lowPt_muon_candidates, bins=40, alpha=0.5, range=(0,75), label='LowPt muon')
plt.hist(highPt_muon_candidates, bins=40, alpha=0.5, range=(0,75), label='HighPt muon')

# Annotate mean and standard deviation on the plot on the left top corner
plt.text(0.05, 0.9, f'MeanLow: {np.mean(lowPt_muon_candidates):.2f}', ha='left', va='center', transform=plt.gca().transAxes)
plt.text(0.05, 0.85, f'StdLow: {np.std(lowPt_muon_candidates):.2f}', ha='left', va='center', transform=plt.gca().transAxes)

plt.text(0.05, 0.8, f'MeanHigh: {np.mean(highPt_muon_candidates):.2f}', ha='left', va='center', transform=plt.gca().transAxes)
plt.text(0.05, 0.75, f'StdHigh: {np.std(highPt_muon_candidates):.2f}', ha='left', va='center', transform=plt.gca().transAxes)

plt.xlabel("Muon pT [GeV]")
plt.legend()
plt.tight_layout()
plt.savefig(f"{FIGURE_PATH}/muon_pt.png")


# plot the Y and Z vertex probability
plt.figure()
plt.hist(Y_VtxProb, bins=40, alpha=0.5, range=(0,1), label='Y vertex probability')
plt.hist(Z_VtxProb, bins=40, alpha=0.5, range=(0,1), label='Z vertex probability')

plt.xlabel("Vertex probability")
plt.legend()
plt.tight_layout()
plt.savefig(f"{FIGURE_PATH}/vertex_prob.png")


# plot the four lepton mass
plt.figure()
plt.hist(FourL_mass, bins=50, alpha=0.5, range=(90,150), label='Four lepton mass')

# Annotate mean and standard deviation on the plot on the left top corner
plt.text(0.05, 0.9, f'Mean: {np.mean(FourL_mass):.2f}', ha='left', va='center', transform=plt.gca().transAxes)
plt.text(0.05, 0.85, f'Std: {np.std(FourL_mass):.2f}', ha='left', va='center', transform=plt.gca().transAxes)
plt.text(0.05, 0.8, f'Candi/Eve: {NCandidates}/{NEvents}', ha='left', va='center', transform=plt.gca().transAxes)

plt.xlabel("Four lepton mass [GeV]")
plt.legend()
plt.tight_layout()
plt.savefig(f"{FIGURE_PATH}/H_mass.png")


# plot the four lepton vertex probability
plt.figure()
plt.hist(FourL_VtxProb, bins=40, alpha=0.5, range=(0,1), label='Four lepton vertex probability')

# Annotate mean and standard deviation on the plot on the left top corner
plt.text(0.05, 0.9, f'Mean: {np.mean(FourL_VtxProb):.2f}', ha='left', va='center', transform=plt.gca().transAxes)
plt.text(0.05, 0.85, f'Std: {np.std(FourL_VtxProb):.2f}', ha='left', va='center', transform=plt.gca().transAxes)

plt.xlabel("Four lepton vertex probability")
plt.legend()
plt.tight_layout()
plt.savefig(f"{FIGURE_PATH}/H_vertex_prob.png")

