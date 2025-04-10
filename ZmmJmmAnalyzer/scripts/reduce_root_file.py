import uproot
import awkward as ak

base_path = '/uscms/home/wkarunar/nobackup/datasets/data/run3/parkingDoubleMuonLowMass/2023/allTrig/'

# List of input files
input_files = [
    'crab_TTree_136TeV_PDMLM0_mm_2023D_v1.root',
    'crab_TTree_136TeV_PDMLM4_mm_2023D_v1.root',
    'crab_TTree_136TeV_PDMLM5_mm_2023D_v1.root',
    'crab_TTree_136TeV_PDMLM7_mm_2023D_v1.root'
]

# Variables to keep
branches = [
    "Event", "Run", "LumiBlock",
    "B_J1_mass", "B_J1_pt",
    "B_J1_VtxPt", "B_J1_VtxProb", "B_J1_VtxMass",
    "B_Mu1_pt", "B_Mu2_pt", "B_M1_pt", "B_M2_pt",
    "B_Mu1_soft", "B_Mu2_soft", "B_Mu1_loose", "B_Mu2_loose", 
    "savedtriggerNames", "savedtriggerBits"
]

# Trigger to check
trigger_tag = "HLT_Mu0_L1DoubleMu_v5"
batch_size = 10000  # Process 10k events at a time

def process_file(base_path, file_path, trigger_tag):
    """
    Process a single ROOT file and filter events based on the trigger.
    """
    # Open the ROOT file
    file_in = uproot.open(base_path + '/' + file_path)
    tree_in = file_in["ntuple"]
    nEntries = tree_in.num_entries

    # Create a new ROOT file for this input
    output_file = f"{file_path[:-5]}_{trigger_tag}.root"

    with uproot.recreate(output_file) as fout:
        fout.mktree("tree", {
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
            "B_Mu1_soft": "bool",
            "B_Mu2_soft": "bool",
            "B_Mu1_loose": "bool",
            "B_Mu2_loose": "bool",
        })

        for i, arrays in enumerate(tree_in.iterate(branches, step_size=batch_size)):
            print(f"Processing {i*batch_size}/{nEntries}")

            # Find events where the trigger fired
            trigger_fired = ak.any(
                (arrays["savedtriggerNames"] == trigger_tag) & (arrays["savedtriggerBits"]),
                axis=1
            )

            selected = arrays[trigger_fired]

            if len(selected) > 0:
                output = {
                    "Event": selected["Event"],
                    "Run": selected["Run"],
                    "LumiBlock": selected["LumiBlock"],
                    "B_J1_mass": selected["B_J1_mass"],
                    "B_J1_pt": selected["B_J1_pt"],
                    "B_J1_VtxProb": selected["B_J1_VtxProb"],
                    "B_J1_VtxMass": selected["B_J1_VtxMass"],
                    "B_Mu1_pt": selected["B_Mu1_pt"],
                    "B_Mu2_pt": selected["B_Mu2_pt"],
                    "B_M1_pt": selected["B_M1_pt"],
                    "B_M2_pt": selected["B_M2_pt"],
                    "B_Mu1_soft": selected["B_Mu1_soft"],
                    "B_Mu2_soft": selected["B_Mu2_soft"],
                    "B_Mu1_loose": selected["B_Mu1_loose"],
                    "B_Mu2_loose": selected["B_Mu2_loose"],
                }
                fout["tree"].extend(output)

    file_in.close()  # Close the input file
    print(f"Finished writing {output_file}")
    

for file_idx, file_path in enumerate(input_files):
    print(f"Opening file {file_idx+1}/{len(input_files)}: {file_path}")
    process_file(base_path, file_path, trigger_tag)

print("All files done!")
