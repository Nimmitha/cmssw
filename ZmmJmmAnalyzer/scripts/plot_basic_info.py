import matplotlib.pyplot as plt
import uproot as up
import numpy as np
import pandas as pd
import awkward as ak
import mplhep as hep
hep.style.use(hep.style.CMS)

fname = 'zmmjmm_2024_final_unblinded'
fpath = '../selection/' + fname + '.root'
file = up.open(fpath)

events = file['ntuple;1']
columns = ['Event', 'B_JPsi_mass', 'B_Z_mass', 'FourL_mass', 'FourL_VtxProb']
branches = events.arrays(columns)

# convert to pandas dataframe
data_dict = {key: ak.to_list(branches[key][0]) for key in branches.fields}
df_candi = pd.DataFrame(data_dict)

# for duplicated events, keep the one with the highest FourL_VtxProb
df_events = df_candi.sort_values('FourL_VtxProb', ascending=False).drop_duplicates('Event').sort_index()

ncandi = len(df_candi)
nevents = len(df_events)

# show extra candidates
print(df_candi[df_candi.duplicated(subset='Event', keep=False)].sort_values('Event', ascending=False))

plt.figure(figsize=(8, 8))
nbins, xlow, xhigh = 10, 112, 162
# plt.hist(df_candi['FourL_mass'], bins=nbins, range=(xlow, xhigh), color='r', alpha=0.2, label='Candidates')
plt.hist(df_events['FourL_mass'], bins=nbins, range=(xlow, xhigh), color='b', alpha=0.5, label=f'Events ({nevents})')
plt.ylabel(f"Counts / {(xhigh-xlow)/nbins:.2f} GeV")
plt.xlabel("FourL_mass [GeV]")
plt.legend()
plt.tight_layout()
plt.savefig(f"{fname}_fourL.png")

plt.figure(figsize=(8, 8))
# nbins, xlow, xhigh = 7, 9, 9.7
nbins, xlow, xhigh = 8, 3.0, 3.2
# plt.hist(df_candi['B_JPsi_mass'], bins=nbins, range=(xlow, xhigh), color='r', alpha=0.2, label='Candidates')
plt.hist(df_events['B_JPsi_mass'], bins=nbins, range=(xlow, xhigh), color='b', alpha=0.5, label=f'Events ({nevents})')
plt.ylabel(f"Counts / {(xhigh-xlow)/nbins:.3f} GeV")
plt.xlabel("Dimuon inv. mass [GeV]")
plt.legend()
plt.tight_layout()
plt.savefig(f"{fname}_JPsi.png")

plt.figure(figsize=(8, 8))
nbins, xlow, xhigh = 10, 80, 100
# plt.hist(df_candi['B_Z_mass'], bins=nbins, range=(xlow, xhigh), color='r', alpha=0.2, label='Candidates')
plt.hist(df_events['B_Z_mass'], bins=nbins, range=(xlow, xhigh), color='b', alpha=0.5, label=f'Events ({nevents})')
plt.ylabel(f"Counts / {(xhigh-xlow)/nbins:.2f} GeV")
plt.xlabel("Dimuon inv. mass [GeV]")
plt.legend()
plt.tight_layout()
plt.savefig(f"{fname}_Z.png")
