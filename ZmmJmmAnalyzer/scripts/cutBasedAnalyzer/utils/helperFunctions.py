import sys
import yaml
import uproot
import awkward as ak
from coffea.nanoevents import TreeMakerSchema, BaseSchema, NanoEventsFactory


def load_config(config_path):
    """
    Load configuration from a YAML file.
    """
    try:
        with open(config_path, 'r') as config_file:
            return yaml.safe_load(config_file)
    except FileNotFoundError:
        print(f"Error: Configuration file '{config_path}' not found.")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing YAML configuration file: {e}")
        sys.exit(1)


def load_events(file_path, schema_class, entry_stop):
    """
    Load events from a ROOT file with error handling.
    """
    try:
        events = NanoEventsFactory.from_root(
            {file_path: "ntuple"},
            schemaclass=schema_class,
            entry_stop=entry_stop
        ).events()
        return events
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return None
    except Exception as e:
        print(f"Error loading events from '{file_path}': {e}")
        return None


def filter_and_flatten(events):
    """
    Filter non-empty fields and flatten them.
    """
    # Keep only records where at least one field is non-empty
    filtered = events[ak.any([ak.num(events[field], axis=1) > 0 for field in events.fields], axis=0)]

    # Flatten fields that are lists
    flattened = {key: ak.flatten(filtered[key], axis=1) for key in filtered.fields}

    return flattened


def save_events(events, isMC):
    """
    Save events to a CSV file and a ROOT file.
    """
    print("\nSaving final tree to CSV and ROOT files...")
    suffix = "mc" if isMC else "data"
    ak.to_dataframe(events).reset_index().to_csv(f'selection_{suffix}.csv', index=False)
    flattened_tree = filter_and_flatten(events)

    with uproot.recreate(f"selection_{suffix}.root") as file:
        file["ntuple"] = flattened_tree


def combine_dimuon_candidates(data, J_lower, J_upper, z_lower, z_upper):
    """
    Combine dimuon candidates into Upsilon candidates. It looks for possible combinations of the
    J1-J2 and J3-J4 pairs that could form Upsilon candidates. The selection is based on the vertex
    probability of the J1-J2 and J3-J4 pairs. The Upsilon candidate with the better vertex probability
    is selected. If there is only one Upsilon candidate, it is selected.
    The function adds columns to the events data and returns the updated array.
    """
    B_J1_mass = data.B_J1_mass
    B_J2_mass = data.B_J2_mass
    B_J3_mass = data.B_J3_mass
    B_J4_mass = data.B_J4_mass
    B_J1_VtxProb = data.B_J1_VtxProb
    B_J2_VtxProb = data.B_J2_VtxProb
    B_J3_VtxProb = data.B_J3_VtxProb
    B_J4_VtxProb = data.B_J4_VtxProb
    B_Mu1_pt = data.B_Mu1_pt
    B_Mu2_pt = data.B_Mu2_pt
    B_Mu3_pt = data.B_Mu3_pt
    B_Mu4_pt = data.B_Mu4_pt

    # Create masks for Upsilon mass ranges
    x = ((B_J1_mass > J_lower) & (B_J1_mass < J_upper)) | ((B_J1_mass > z_lower) & (B_J1_mass < z_upper))
    y = ((B_J2_mass > J_lower) & (B_J2_mass < J_upper)) | ((B_J2_mass > z_lower) & (B_J2_mass < z_upper))
    ups_mask1 = x & y

    x = ((B_J3_mass > J_lower) & (B_J3_mass < J_upper)) | ((B_J3_mass > z_lower) & (B_J3_mass < z_upper))
    y = ((B_J4_mass > J_lower) & (B_J4_mass < J_upper)) | ((B_J4_mass > z_lower) & (B_J4_mass < z_upper))

    ups_mask2 = x & y

    # Count Upsilon candidates
    JPsiMass = ak.values_astype(ups_mask1, "int") + ak.values_astype(ups_mask2, "int")

    # Initialize output arrays
    JPsi_mass = ak.zeros_like(B_J1_mass)
    Z_mass = ak.zeros_like(B_J1_mass)
    JPsi_VtxProb = ak.zeros_like(B_J1_mass)
    Z_VtxProb = ak.zeros_like(B_J1_mass)
    JPsi_Pt = ak.zeros_like(B_J1_mass)
    Z_Pt = ak.zeros_like(B_J1_mass)
    # JPsi_Eta = ak.zeros_like(B_J1_mass)
    # Z_Eta = ak.zeros_like(B_J1_mass)
    # JPsi_Rapidity = ak.zeros_like(B_J1_mass)
    # Z_Rapidity = ak.zeros_like(B_J1_mass)
    # JPsi_Phi = ak.zeros_like(B_J1_mass)
    # Z_Phi = ak.zeros_like(B_J1_mass)

    # Process events with 2 Upsilon candidates
    two_ups_mask = (JPsiMass == 2)
    vtx_prob_sum1 = B_J1_VtxProb + B_J2_VtxProb
    vtx_prob_sum2 = B_J3_VtxProb + B_J4_VtxProb
    better_vtx_mask = vtx_prob_sum1 > vtx_prob_sum2

    # For events with 2 Upsilon candidates and better vertex probability in J1-J2 pair
    mask_2ups_better12 = two_ups_mask & better_vtx_mask
    j1_is_ups = (B_J1_mass > J_lower) & (B_J1_mass < J_upper)
    JPsi_mass = ak.where(mask_2ups_better12 & j1_is_ups, B_J1_mass, JPsi_mass)
    Z_mass = ak.where(mask_2ups_better12 & j1_is_ups, B_J2_mass, Z_mass)
    JPsi_mass = ak.where(mask_2ups_better12 & ~j1_is_ups, B_J2_mass, JPsi_mass)
    Z_mass = ak.where(mask_2ups_better12 & ~j1_is_ups, B_J1_mass, Z_mass)

    # Similar assignments for other variables (VtxProb, Pt, Eta, Rapidity, Phi)
    # ... (omitted for brevity, but follow the same pattern as mass assignments)
    JPsi_VtxProb = ak.where(mask_2ups_better12 & j1_is_ups, B_J1_VtxProb, JPsi_VtxProb)
    Z_VtxProb = ak.where(mask_2ups_better12 & j1_is_ups, B_J2_VtxProb, Z_VtxProb)
    JPsi_VtxProb = ak.where(mask_2ups_better12 & ~j1_is_ups, B_J2_VtxProb, JPsi_VtxProb)
    Z_VtxProb = ak.where(mask_2ups_better12 & ~j1_is_ups, B_J1_VtxProb, Z_VtxProb)

    JPsi_Pt = ak.where(mask_2ups_better12 & j1_is_ups, B_Mu1_pt, JPsi_Pt)
    Z_Pt = ak.where(mask_2ups_better12 & j1_is_ups, B_Mu2_pt, Z_Pt)
    JPsi_Pt = ak.where(mask_2ups_better12 & ~j1_is_ups, B_Mu2_pt, JPsi_Pt)
    Z_Pt = ak.where(mask_2ups_better12 & ~j1_is_ups, B_Mu1_pt, Z_Pt)

    # For events with 2 Upsilon candidates and better vertex probability in J3-J4 pair
    mask_2ups_better34 = two_ups_mask & ~better_vtx_mask
    j3_is_ups = (B_J3_mass > J_lower) & (B_J3_mass < J_upper)
    JPsi_mass = ak.where(mask_2ups_better34 & j3_is_ups, B_J3_mass, JPsi_mass)
    Z_mass = ak.where(mask_2ups_better34 & j3_is_ups, B_J4_mass, Z_mass)
    JPsi_mass = ak.where(mask_2ups_better34 & ~j3_is_ups, B_J4_mass, JPsi_mass)
    Z_mass = ak.where(mask_2ups_better34 & ~j3_is_ups, B_J3_mass, Z_mass)

    # Similar assignments for other variables (VtxProb, Pt, Eta, Rapidity, Phi)
    # ... (omitted for brevity, but follow the same pattern as mass assignments)
    JPsi_VtxProb = ak.where(mask_2ups_better34 & j3_is_ups, B_J3_VtxProb, JPsi_VtxProb)
    Z_VtxProb = ak.where(mask_2ups_better34 & j3_is_ups, B_J4_VtxProb, Z_VtxProb)
    JPsi_VtxProb = ak.where(mask_2ups_better34 & ~j3_is_ups, B_J4_VtxProb, JPsi_VtxProb)
    Z_VtxProb = ak.where(mask_2ups_better34 & ~j3_is_ups, B_J3_VtxProb, Z_VtxProb)

    JPsi_Pt = ak.where(mask_2ups_better34 & j3_is_ups, B_Mu3_pt, JPsi_Pt)
    Z_Pt = ak.where(mask_2ups_better34 & j3_is_ups, B_Mu4_pt, Z_Pt)
    JPsi_Pt = ak.where(mask_2ups_better34 & ~j3_is_ups, B_Mu4_pt, JPsi_Pt)
    Z_Pt = ak.where(mask_2ups_better34 & ~j3_is_ups, B_Mu3_pt, Z_Pt)

    # Process events with 1 Upsilon candidate
    one_ups_mask = (JPsiMass == 1)
    j1j2_ups_mask = one_ups_mask & ups_mask1
    j3j4_ups_mask = one_ups_mask & ups_mask2

    # For J1-J2 pair
    j1_is_ups = (B_J1_mass > J_lower) & (B_J1_mass < J_upper)
    JPsi_mass = ak.where(j1j2_ups_mask & j1_is_ups, B_J1_mass, JPsi_mass)
    Z_mass = ak.where(j1j2_ups_mask & j1_is_ups, B_J2_mass, Z_mass)
    JPsi_mass = ak.where(j1j2_ups_mask & ~j1_is_ups, B_J2_mass, JPsi_mass)
    Z_mass = ak.where(j1j2_ups_mask & ~j1_is_ups, B_J1_mass, Z_mass)

    # Similar assignments for other variables (VtxProb, Pt, Eta, Rapidity, Phi)
    # ... (omitted for brevity, but follow the same pattern as mass assignments)
    JPsi_VtxProb = ak.where(j1j2_ups_mask & j1_is_ups, B_J1_VtxProb, JPsi_VtxProb)
    Z_VtxProb = ak.where(j1j2_ups_mask & j1_is_ups, B_J2_VtxProb, Z_VtxProb)
    JPsi_VtxProb = ak.where(j1j2_ups_mask & ~j1_is_ups, B_J2_VtxProb, JPsi_VtxProb)
    Z_VtxProb = ak.where(j1j2_ups_mask & ~j1_is_ups, B_J1_VtxProb, Z_VtxProb)

    JPsi_Pt = ak.where(j1j2_ups_mask & j1_is_ups, B_Mu1_pt, JPsi_Pt)
    Z_Pt = ak.where(j1j2_ups_mask & j1_is_ups, B_Mu2_pt, Z_Pt)
    JPsi_Pt = ak.where(j1j2_ups_mask & ~j1_is_ups, B_Mu2_pt, JPsi_Pt)
    Z_Pt = ak.where(j1j2_ups_mask & ~j1_is_ups, B_Mu1_pt, Z_Pt)

    # For J3-J4 pair
    j3_is_ups = (B_J3_mass > J_lower) & (B_J3_mass < J_upper)
    JPsi_mass = ak.where(j3j4_ups_mask & j3_is_ups, B_J3_mass, JPsi_mass)
    Z_mass = ak.where(j3j4_ups_mask & j3_is_ups, B_J4_mass, Z_mass)
    JPsi_mass = ak.where(j3j4_ups_mask & ~j3_is_ups, B_J4_mass, JPsi_mass)
    Z_mass = ak.where(j3j4_ups_mask & ~j3_is_ups, B_J3_mass, Z_mass)

    # Similar assignments for other variables (VtxProb, Pt, Eta, Rapidity, Phi)
    # ... (omitted for brevity, but follow the same pattern as mass assignments)
    JPsi_VtxProb = ak.where(j3j4_ups_mask & j3_is_ups, B_J3_VtxProb, JPsi_VtxProb)
    Z_VtxProb = ak.where(j3j4_ups_mask & j3_is_ups, B_J4_VtxProb, Z_VtxProb)
    JPsi_VtxProb = ak.where(j3j4_ups_mask & ~j3_is_ups, B_J4_VtxProb, JPsi_VtxProb)
    Z_VtxProb = ak.where(j3j4_ups_mask & ~j3_is_ups, B_J3_VtxProb, Z_VtxProb)

    JPsi_Pt = ak.where(j3j4_ups_mask & j3_is_ups, B_Mu3_pt, JPsi_Pt)
    Z_Pt = ak.where(j3j4_ups_mask & j3_is_ups, B_Mu4_pt, Z_Pt)
    JPsi_Pt = ak.where(j3j4_ups_mask & ~j3_is_ups, B_Mu4_pt, JPsi_Pt)
    Z_Pt = ak.where(j3j4_ups_mask & ~j3_is_ups, B_Mu3_pt, Z_Pt)

    # valid_mask = (JPsiMass > 0)
    data['JPsi_mass'] = JPsi_mass
    data['Z_mass'] = Z_mass
    data['JPsiMass'] = JPsiMass
    data['JPsi_VtxProb'] = JPsi_VtxProb
    data['Z_VtxProb'] = Z_VtxProb
    data['JPsi_Pt'] = JPsi_Pt
    data['Z_Pt'] = Z_Pt
    # data['JPsi_Eta'] = JPsi_Eta
    # data['Z_Eta'] = Z_Eta

    return data
