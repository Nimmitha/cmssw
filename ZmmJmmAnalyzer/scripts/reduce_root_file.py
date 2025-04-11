import uproot
import awkward as ak

base_path = '/uscms/home/wkarunar/nobackup/datasets/data/run3/parkingDoubleMuonLowMass/2023/allTrig/'
# base_path = '/uscms/home/wkarunar/nobackup/datasets/data/run3/Muon/jmm/'

input_files = [
    # 'crab_TTree_136TeV_Muon0_mm_2023D_v1.root',
    # 'crab_TTree_136TeV_Muon1_mm_2023D_v1.root',
    'crab_TTree_136TeV_PDMLM0_mm_2023D_v1.root',
    'crab_TTree_136TeV_PDMLM4_mm_2023D_v1.root',
    'crab_TTree_136TeV_PDMLM5_mm_2023D_v1.root',
    'crab_TTree_136TeV_PDMLM7_mm_2023D_v1.root',
]

# Trigger to check
trigger_tag = "HLT_Mu0_L1DoubleMu_v5"
# trigger_tag = "HLT_IsoMu24_v17"
batch_size = 10000  # Process 10k events at a time

# List of variables you want to save
variables_to_save = [
    "Event", "Run", "LumiBlock",
    "B_J1_mass", "B_J1_pt", "B_J1_VtxProb", "B_J1_VtxMass",
    "B_Mu1_pt", "B_Mu2_pt", "B_M1_pt", "B_M2_pt",
    "B_Mu1_eta", "B_Mu2_eta", "B_M1_eta", "B_M2_eta",
    "B_Mu1_soft", "B_Mu2_soft", "B_Mu1_loose", "B_Mu2_loose",
    "B_Mu1_tight", "B_Mu2_tight"
]

# Data types (you can edit this easily later too)
output_branches_types = {
    "Event": "uint64",
    "Run": "uint32",
    "LumiBlock": "uint32",
    "B_J1_mass": "float32",
    "B_J1_pt": "float32",
    "B_J1_VtxProb": "float32",
    "B_J1_VtxMass": "float32",
    "B_Mu1_pt": "float32",
    "B_Mu2_pt": "float32",
    "B_M1_pt": "float32",
    "B_M2_pt": "float32",
    "B_Mu1_eta": "float32",
    "B_Mu2_eta": "float32",
    "B_M1_eta": "float32",
    "B_M2_eta": "float32",
    "B_Mu1_soft": "bool",
    "B_Mu2_soft": "bool",
    "B_Mu1_loose": "bool",
    "B_Mu2_loose": "bool",
    "B_Mu1_tight": "bool",
    "B_Mu2_tight": "bool",
    "TriggerFired": "bool",
}

# Branches needed from input file
branches_to_read = variables_to_save + ["savedtriggerNames", "savedtriggerBits"]

def process_file(base_path, file_path, trigger_tag):
    """
    Process a single ROOT file and filter events based on the trigger.
    """
    # Open the ROOT file
    file_in = uproot.open(base_path + '/' + file_path)
    tree_in = file_in["ntuple"]
    nEntries = tree_in.num_entries

    # Create output file
    output_file = f"{file_path[:-5]}_{trigger_tag}.root"

    with uproot.recreate(output_file) as fout:
        fout.mktree("tree", output_branches_types)

        for i, arrays in enumerate(tree_in.iterate(branches_to_read, step_size=batch_size)):
            print(f"Processing {i*batch_size}/{nEntries}")

            # Find if trigger fired (bool array)
            trigger_fired = ak.any(
                (arrays["savedtriggerNames"] == trigger_tag) & (arrays["savedtriggerBits"]),
                axis=1
            )

            output = {}
            for var in variables_to_save:
                output[var] = arrays[var]
            # Add the trigger status (always True because we select on it, but you could generalize)
            output["TriggerFired"] = trigger_fired

            fout["tree"].extend(output)

    file_in.close()
    print(f"Finished writing {output_file}")

for file_idx, file_path in enumerate(input_files):
    print(f"Opening file {file_idx+1}/{len(input_files)}: {file_path}")
    process_file(base_path, file_path, trigger_tag)

print("All files done!")
