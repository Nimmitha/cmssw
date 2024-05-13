import pandas as pd
import awkward as ak

def get_count(data):
    array = data.Y_mass
    
    nevents = len(array[ak.num(array, axis=1) > 0])
    ncandidates = ak.sum(ak.num(array, axis=1))

    return ncandidates, nevents


def get_summary_of_cuts(events, cut_order, cut_list):
    key_list = []
    ncandidates_list = []
    nevents_list = []
    aggMask_list = []
    var_list = []

    agg_mask = cut_list["0"]["mask"]

    for i in cut_order:
        # print(f"Cut {i}: {cut_list[str(i)]['name']}")
        mask = cut_list[str(i)]["mask"]

        agg_mask = agg_mask & mask
        ncandidates, nevents = get_count(events[agg_mask])

        key_list.append(cut_list[str(i)]["name"])
        ncandidates_list.append(ncandidates)
        nevents_list.append(nevents)   
        aggMask_list.append(agg_mask)
        var_list.append(cut_list[str(i)].get("var", None))

        # make a dictionary from the four lists
        summary_dict = {"Cut": key_list, "Candidates": ncandidates_list, "Events": nevents_list, "Aggregated mask": aggMask_list, "Var": var_list}
        # get selected items from the dictionary
        keys = ["Cut", "Candidates", "Events"]
        summary_table = pd.DataFrame(dict((key, value) for key, value in summary_dict.items() if key in keys))

    return summary_dict, summary_table



# class MaskClass:
#     def __init__(self, events, var_type, lower_lim=None, upper_lim=None):
#         self.var_type = var_type
#         self.lower_lim = lower_lim
#         self.upper_lim = upper_lim
#         self.mask = self.getMask(events)

#     def getMask(self, events):
#         if (self.lower_lim is not None) and (self.upper_lim is not None):
#             mask = (events[self.var_type] > self.lower_lim) & (events[self.var_type] < self.upper_lim)
#         elif self.upper_lim is not None:
#             mask = events[self.var_type] < self.upper_lim
#         elif self.lower_lim is not None:
#             mask = events[self.var_type] > self.lower_lim
#         else:
#             mask = events[self.var_type]
#         return mask

# my_mask = MaskClass(events, "Z_pt1", lower_lim=24)



# def get_summary_of_cuts_old(events, primary_cuts, mass_cuts, electronID_cuts):

#     key_list = []
#     ncandidates_list = []
#     nevents_list = []
#     aggMask_list = []

#     agg_mask = primary_cuts['Preselection']

#     for key, mask in primary_cuts.items():
#         # print(key)
        
#         agg_mask = agg_mask & mask

#         ncandidates, nevents = get_count(events[agg_mask])
#         # print(f"Number of candidates: {ncandidates}")
#         # print(f"Number of events: {nevents}")

#         key_list.append(key)
#         ncandidates_list.append(ncandidates)
#         nevents_list.append(nevents)
#         aggMask_list.append(agg_mask)
        
#     for key, mask in mass_cuts.items():
#         # print(key)

#         agg_mask = agg_mask & mask

#         ncandidates, nevents = get_count(events[agg_mask])
#         # print(f"Number of candidates: {ncandidates}")
#         # print(f"Number of events: {nevents}")

#         key_list.append(key)
#         ncandidates_list.append(ncandidates)
#         nevents_list.append(nevents)
#         aggMask_list.append(agg_mask)

#     temp_cuts = []
#     temp_ncandidates = []
#     temp_nevents = []
#     temp_aggMask = []

#     for key, mask in electronID_cuts.items():
#         # print(key)

#         temp_mask = agg_mask & mask

#         ncandidates, nevents = get_count(events[temp_mask])
#         # print(f"Number of candidates: {ncandidates}")
#         # print(f"Number of events: {nevents}")

#         temp_cuts.append(key)
#         temp_ncandidates.append(ncandidates)
#         temp_nevents.append(nevents)
#         temp_aggMask.append(temp_mask)

#     key_list.append(temp_cuts)
#     ncandidates_list.append(temp_ncandidates)
#     nevents_list.append(temp_nevents)
#     aggMask_list.append(temp_aggMask)

#     # Make a table to summarize cuts
#     cut_summary_table = pd.DataFrame(
#         {"Cuts": key_list,
#         "Candidates": ncandidates_list,
#         "Events": nevents_list}
#     )

#     return cut_summary_table, aggMask_list
