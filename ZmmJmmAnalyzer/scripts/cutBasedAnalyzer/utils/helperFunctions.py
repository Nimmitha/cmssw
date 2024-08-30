import sys
import yaml
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
    UpsMass = ak.values_astype(ups_mask1, "int") + ak.values_astype(ups_mask2, "int")

    # Initialize output arrays
    Ups1_mass = ak.zeros_like(B_J1_mass)
    Ups2_mass = ak.zeros_like(B_J1_mass)
    Ups_VtxProb1 = ak.zeros_like(B_J1_mass)
    Ups_VtxProb2 = ak.zeros_like(B_J1_mass)
    Ups_Pt1 = ak.zeros_like(B_J1_mass)
    Ups_Pt2 = ak.zeros_like(B_J1_mass)
    # Ups_Eta1 = ak.zeros_like(B_J1_mass)
    # Ups_Eta2 = ak.zeros_like(B_J1_mass)
    # Ups_Rapidity1 = ak.zeros_like(B_J1_mass)
    # Ups_Rapidity2 = ak.zeros_like(B_J1_mass)
    # Ups_Phi1 = ak.zeros_like(B_J1_mass)
    # Ups_Phi2 = ak.zeros_like(B_J1_mass)

    # Process events with 2 Upsilon candidates
    two_ups_mask = (UpsMass == 2)
    vtx_prob_sum1 = B_J1_VtxProb + B_J2_VtxProb
    vtx_prob_sum2 = B_J3_VtxProb + B_J4_VtxProb
    better_vtx_mask = vtx_prob_sum1 > vtx_prob_sum2

    # For events with 2 Upsilon candidates and better vertex probability in J1-J2 pair
    mask_2ups_better12 = two_ups_mask & better_vtx_mask
    j1_is_ups = (B_J1_mass > J_lower) & (B_J1_mass < J_upper)
    Ups1_mass = ak.where(mask_2ups_better12 & j1_is_ups, B_J1_mass, Ups1_mass)
    Ups2_mass = ak.where(mask_2ups_better12 & j1_is_ups, B_J2_mass, Ups2_mass)
    Ups1_mass = ak.where(mask_2ups_better12 & ~j1_is_ups, B_J2_mass, Ups1_mass)
    Ups2_mass = ak.where(mask_2ups_better12 & ~j1_is_ups, B_J1_mass, Ups2_mass)

    # Similar assignments for other variables (VtxProb, Pt, Eta, Rapidity, Phi)
    # ... (omitted for brevity, but follow the same pattern as mass assignments)
    Ups_VtxProb1 = ak.where(mask_2ups_better12 & j1_is_ups, B_J1_VtxProb, Ups_VtxProb1)
    Ups_VtxProb2 = ak.where(mask_2ups_better12 & j1_is_ups, B_J2_VtxProb, Ups_VtxProb2)
    Ups_VtxProb1 = ak.where(mask_2ups_better12 & ~j1_is_ups, B_J2_VtxProb, Ups_VtxProb1)
    Ups_VtxProb2 = ak.where(mask_2ups_better12 & ~j1_is_ups, B_J1_VtxProb, Ups_VtxProb2)

    Ups_Pt1 = ak.where(mask_2ups_better12 & j1_is_ups, B_Mu1_pt, Ups_Pt1)
    Ups_Pt2 = ak.where(mask_2ups_better12 & j1_is_ups, B_Mu2_pt, Ups_Pt2)
    Ups_Pt1 = ak.where(mask_2ups_better12 & ~j1_is_ups, B_Mu2_pt, Ups_Pt1)
    Ups_Pt2 = ak.where(mask_2ups_better12 & ~j1_is_ups, B_Mu1_pt, Ups_Pt2)

    # For events with 2 Upsilon candidates and better vertex probability in J3-J4 pair
    mask_2ups_better34 = two_ups_mask & ~better_vtx_mask
    j3_is_ups = (B_J3_mass > J_lower) & (B_J3_mass < J_upper)
    Ups1_mass = ak.where(mask_2ups_better34 & j3_is_ups, B_J3_mass, Ups1_mass)
    Ups2_mass = ak.where(mask_2ups_better34 & j3_is_ups, B_J4_mass, Ups2_mass)
    Ups1_mass = ak.where(mask_2ups_better34 & ~j3_is_ups, B_J4_mass, Ups1_mass)
    Ups2_mass = ak.where(mask_2ups_better34 & ~j3_is_ups, B_J3_mass, Ups2_mass)

    # Similar assignments for other variables (VtxProb, Pt, Eta, Rapidity, Phi)
    # ... (omitted for brevity, but follow the same pattern as mass assignments)
    Ups_VtxProb1 = ak.where(mask_2ups_better34 & j3_is_ups, B_J3_VtxProb, Ups_VtxProb1)
    Ups_VtxProb2 = ak.where(mask_2ups_better34 & j3_is_ups, B_J4_VtxProb, Ups_VtxProb2)
    Ups_VtxProb1 = ak.where(mask_2ups_better34 & ~j3_is_ups, B_J4_VtxProb, Ups_VtxProb1)
    Ups_VtxProb2 = ak.where(mask_2ups_better34 & ~j3_is_ups, B_J3_VtxProb, Ups_VtxProb2)

    Ups_Pt1 = ak.where(mask_2ups_better34 & j3_is_ups, B_Mu3_pt, Ups_Pt1)
    Ups_Pt2 = ak.where(mask_2ups_better34 & j3_is_ups, B_Mu4_pt, Ups_Pt2)
    Ups_Pt1 = ak.where(mask_2ups_better34 & ~j3_is_ups, B_Mu4_pt, Ups_Pt1)
    Ups_Pt2 = ak.where(mask_2ups_better34 & ~j3_is_ups, B_Mu3_pt, Ups_Pt2)

    # Process events with 1 Upsilon candidate
    one_ups_mask = (UpsMass == 1)
    j1j2_ups_mask = one_ups_mask & ups_mask1
    j3j4_ups_mask = one_ups_mask & ups_mask2

    # For J1-J2 pair
    j1_is_ups = (B_J1_mass > J_lower) & (B_J1_mass < J_upper)
    Ups1_mass = ak.where(j1j2_ups_mask & j1_is_ups, B_J1_mass, Ups1_mass)
    Ups2_mass = ak.where(j1j2_ups_mask & j1_is_ups, B_J2_mass, Ups2_mass)
    Ups1_mass = ak.where(j1j2_ups_mask & ~j1_is_ups, B_J2_mass, Ups1_mass)
    Ups2_mass = ak.where(j1j2_ups_mask & ~j1_is_ups, B_J1_mass, Ups2_mass)

    # Similar assignments for other variables (VtxProb, Pt, Eta, Rapidity, Phi)
    # ... (omitted for brevity, but follow the same pattern as mass assignments)
    Ups_VtxProb1 = ak.where(j1j2_ups_mask & j1_is_ups, B_J1_VtxProb, Ups_VtxProb1)
    Ups_VtxProb2 = ak.where(j1j2_ups_mask & j1_is_ups, B_J2_VtxProb, Ups_VtxProb2)
    Ups_VtxProb1 = ak.where(j1j2_ups_mask & ~j1_is_ups, B_J2_VtxProb, Ups_VtxProb1)
    Ups_VtxProb2 = ak.where(j1j2_ups_mask & ~j1_is_ups, B_J1_VtxProb, Ups_VtxProb2)

    Ups_Pt1 = ak.where(j1j2_ups_mask & j1_is_ups, B_Mu1_pt, Ups_Pt1)
    Ups_Pt2 = ak.where(j1j2_ups_mask & j1_is_ups, B_Mu2_pt, Ups_Pt2)
    Ups_Pt1 = ak.where(j1j2_ups_mask & ~j1_is_ups, B_Mu2_pt, Ups_Pt1)
    Ups_Pt2 = ak.where(j1j2_ups_mask & ~j1_is_ups, B_Mu1_pt, Ups_Pt2)

    # For J3-J4 pair
    j3_is_ups = (B_J3_mass > J_lower) & (B_J3_mass < J_upper)
    Ups1_mass = ak.where(j3j4_ups_mask & j3_is_ups, B_J3_mass, Ups1_mass)
    Ups2_mass = ak.where(j3j4_ups_mask & j3_is_ups, B_J4_mass, Ups2_mass)
    Ups1_mass = ak.where(j3j4_ups_mask & ~j3_is_ups, B_J4_mass, Ups1_mass)
    Ups2_mass = ak.where(j3j4_ups_mask & ~j3_is_ups, B_J3_mass, Ups2_mass)

    # Similar assignments for other variables (VtxProb, Pt, Eta, Rapidity, Phi)
    # ... (omitted for brevity, but follow the same pattern as mass assignments)
    Ups_VtxProb1 = ak.where(j3j4_ups_mask & j3_is_ups, B_J3_VtxProb, Ups_VtxProb1)
    Ups_VtxProb2 = ak.where(j3j4_ups_mask & j3_is_ups, B_J4_VtxProb, Ups_VtxProb2)
    Ups_VtxProb1 = ak.where(j3j4_ups_mask & ~j3_is_ups, B_J4_VtxProb, Ups_VtxProb1)
    Ups_VtxProb2 = ak.where(j3j4_ups_mask & ~j3_is_ups, B_J3_VtxProb, Ups_VtxProb2)

    Ups_Pt1 = ak.where(j3j4_ups_mask & j3_is_ups, B_Mu3_pt, Ups_Pt1)
    Ups_Pt2 = ak.where(j3j4_ups_mask & j3_is_ups, B_Mu4_pt, Ups_Pt2)
    Ups_Pt1 = ak.where(j3j4_ups_mask & ~j3_is_ups, B_Mu4_pt, Ups_Pt1)
    Ups_Pt2 = ak.where(j3j4_ups_mask & ~j3_is_ups, B_Mu3_pt, Ups_Pt2)

    # valid_mask = (UpsMass > 0)
    data['Ups1_mass'] = Ups1_mass
    data['Ups2_mass'] = Ups2_mass
    data['UpsMass'] = UpsMass
    data['Ups_VtxProb1'] = Ups_VtxProb1
    data['Ups_VtxProb2'] = Ups_VtxProb2
    data['Ups_Pt1'] = Ups_Pt1
    data['Ups_Pt2'] = Ups_Pt2
    # data['Ups_Eta1'] = Ups_Eta1
    # data['Ups_Eta2'] = Ups_Eta2

    return data
