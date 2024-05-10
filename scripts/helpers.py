import pandas as pd
import awkward as ak

def get_count(data):
    array = data.Y_mass
    
    nevents = len(array[ak.num(array, axis=1) > 0])
    ncandidates = ak.sum(ak.num(array, axis=1))

    return ncandidates, nevents

def get_summary_table(events, primary_cuts, mass_cuts, electronID_cuts):

    key_list = []
    ncandidates_list = []
    nevents_list = []
    aggMask_list = []

    agg_mask = primary_cuts['Preselection']

    for key, mask in primary_cuts.items():
        # print(key)
        
        agg_mask = agg_mask & mask

        ncandidates, nevents = get_count(events[agg_mask])
        # print(f"Number of candidates: {ncandidates}")
        # print(f"Number of events: {nevents}")

        key_list.append(key)
        ncandidates_list.append(ncandidates)
        nevents_list.append(nevents)
        aggMask_list.append(agg_mask)
        
    for key, mask in mass_cuts.items():
        # print(key)

        agg_mask = agg_mask & mask

        ncandidates, nevents = get_count(events[agg_mask])
        # print(f"Number of candidates: {ncandidates}")
        # print(f"Number of events: {nevents}")

        key_list.append(key)
        ncandidates_list.append(ncandidates)
        nevents_list.append(nevents)
        aggMask_list.append(agg_mask)

    temp_cuts = []
    temp_ncandidates = []
    temp_nevents = []
    temp_aggMask = []

    for key, mask in electronID_cuts.items():
        # print(key)

        temp_mask = agg_mask & mask

        ncandidates, nevents = get_count(events[temp_mask])
        # print(f"Number of candidates: {ncandidates}")
        # print(f"Number of events: {nevents}")

        temp_cuts.append(key)
        temp_ncandidates.append(ncandidates)
        temp_nevents.append(nevents)
        temp_aggMask.append(temp_mask)

    key_list.append(temp_cuts)
    ncandidates_list.append(temp_ncandidates)
    nevents_list.append(temp_nevents)
    aggMask_list.append(temp_aggMask)

    # Make a table to summarize cuts
    cut_summary_table = pd.DataFrame(
        {"Cuts": key_list,
        "Candidates": ncandidates_list,
        "Events": nevents_list}
    )

    return cut_summary_table, aggMask_list