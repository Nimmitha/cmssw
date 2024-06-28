import pandas as pd
import awkward as ak

def get_count(data):
    array = data.B_J_mass
    
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