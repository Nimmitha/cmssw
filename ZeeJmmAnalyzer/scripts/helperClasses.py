import matplotlib.pyplot as plt
import awkward as ak
import mplhep as hep
hep.style.use(hep.style.CMS)

class CutDict:
    def __init__(self):
        self.cutdict = {}

    def add_cut(self, cut_id, cut):
        self.cutdict[cut_id] = cut

    def get_cut(self, cut_id):
        return self.cutdict.get(cut_id)
    
    def __str__(self):
        return "\n".join([f"{cut_id}: {cut}" for cut_id, cut in self.cutdict.items()])


class CutType:
    def __init__(self, name, long_name, mask, variables=None, plot=True, bins=50, x_range=None, labels=None, xlabel=None, cut_line=None):
        self.name = name
        self.long_name = long_name
        self.mask = mask
        self.variables = variables if variables else []
        self.plot = plot
        self.bins = bins
        self.x_range = x_range
        self.labels = labels if labels else []
        self.xlabel = xlabel
        self.cut_line = cut_line

    def __str__(self):
        return f"Cut: {self.name} ({self.long_name})"
    

class CutAnalysis:
    def __init__(self, events, cutdict):
        self.events = events
        self.cutdict = cutdict

    def get_stats(self, data):
        array = data.B_J_mass
        nevents = len(array[ak.num(array, axis=1) > 0])
        ncandidates = ak.sum(ak.num(array, axis=1))
        return nevents, ncandidates
    
    def prepare_masks(self, myorder):
        masks_list = []
        cutobj_list = []
        ncandidates_list = []
        nevents_list = []

        aggregate_mask = self.cutdict.cutdict[0].mask
        
        for cut_id in myorder:
            cut = self.cutdict.cutdict[cut_id]
            aggregate_mask = aggregate_mask & cut.mask
            masks_list.append(aggregate_mask)
            cutobj_list.append(cut)
            ncandidates, nevents = self.get_stats(self.events[aggregate_mask])
            ncandidates_list.append(ncandidates)
            nevents_list.append(nevents)

            print(f"Cut {cut_id}: {cut.long_name} - {nevents} events, {ncandidates} candidates")

        summary_dict = {
            "cut_id": myorder,
            "cut_name": [self.cutdict.cutdict[cut_id].long_name for cut_id in myorder],
            "nevents": nevents_list,
            "ncandidates": ncandidates_list,
            "mask": masks_list,
            "cutobj": cutobj_list
        }

        return summary_dict


class Plotter:
    def __init__(self, events, cutdict, save_path):
        self.events = events
        self.cutdict = cutdict
        self.save_path = save_path

    def plot_preselection(self):
        for cut_id in self.cutdict.cutdict:
            cut = self.cutdict.get_cut(cut_id)
            if not cut.plot:
                print(f"Cut {cut_id} does not have a plot.")
                continue
            
            print(f"Plotting cut {cut_id}: {cut.long_name}")
            if len(cut.variables) == 0:
                print(f"Cut {cut_id} does not have any variables to plot.")
                continue

            nbins = cut.bins
            xrange = cut.x_range
            labels = cut.labels
            xlabel = cut.xlabel
            variables = cut.variables
            unit = 'GeV' if 'eta' not in xlabel else ''
            
            plt.figure()
            for i, variable in enumerate(variables):
                plt.hist(ak.flatten(self.events[variable][cut.mask]), bins=nbins, range=xrange, label=labels[i], alpha=0.5)
            plt.xlabel(cut.labels[0])
            plt.ylabel("Events")
            plt.title(f"{cut.long_name}")
            plt.legend()
            plt.savefig(f"{self.save_path}/pre_{cut.long_name}.png")
            plt.close()

    def plot_summary_at(self, cut_id, cutSummary):
        
        idx_after_plot = cutSummary["cut_id"].index(cut_id)
        name_after_plot = cutSummary["cut_name"][idx_after_plot]
        n_eve_after = cutSummary["nevents"][idx_after_plot]
        print(f"Next cut: {idx_after_plot} {name_after_plot} ({n_eve_after})")

        cut_obj = cutSummary["cutobj"][idx_after_plot]

        if not cut_obj.plot:
            print(f"Cut {cut_id} does not have a plot.")
            return
        
        variables = cut_obj.variables
        nbins = cut_obj.bins
        xlow, xhigh = cut_obj.x_range
        labels = cut_obj.labels
        xlabel = cut_obj.xlabel
        unit = 'GeV' if 'eta' not in xlabel else ''
        fileName = f"{self.save_path}/cut_{name_after_plot}.png"

        idx_at_plot = idx_after_plot - 1
        name_at_plot = cutSummary["cut_name"][idx_at_plot]
        n_eve_at_plot = cutSummary["nevents"][idx_at_plot]
        mask = cutSummary["mask"][idx_at_plot]
        print(f"Drawn at: {idx_at_plot} {name_at_plot} ({n_eve_at_plot})")

        for i, variable in enumerate(variables):
            plt.hist(ak.flatten(self.events[variable][mask]), bins=nbins, range=(xlow, xhigh), label=labels[i])

        plt.text(0.5, 0.5, f"Drawn at: {name_at_plot} ({n_eve_at_plot})", fontsize=12, transform=plt.gca().transAxes)
        plt.text(0.5, 0.45, f"Next cut: {name_after_plot} ({n_eve_after})", fontsize=12, transform=plt.gca().transAxes)

        plt.xlabel(xlabel)
        plt.ylabel(f"Counts / {(xhigh-xlow)/nbins:.3f} {unit}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(fileName)
        plt.close()

    def plot_for_single_variable(self, variable, cutdict, cutSummary):
        # find the cut which has the variable
        cut_id = None
        for cut_id in self.cutdict.cutdict:
            if variable in self.cutdict.cutdict[cut_id].variables:
                break

        if cut_id is None:
            print(f"Variable {variable} not found in any cutdict.")
            return
        
        print(f"Plotting variable {variable} for cutdict: {cutdict}")
        
        nbins = self.cutdict.cutdict[cut_id].bins
        xlow, xhigh = self.cutdict.cutdict[cut_id].x_range
        xlabel = self.cutdict.cutdict[cut_id].xlabel
        unit = 'GeV' if 'eta' not in variable else ''

        cut_names = [self.cutdict.get_cut(cut_id).long_name for cut_id in cutdict]
        print(f"Plotting variable {variable} for cutdict: {cut_names}")

        plt.figure()
        for cut_id in cutdict:
            local_idx = cutSummary["cut_id"].index(cut_id)
            mask = cutSummary["mask"][local_idx]
            nevents = cutSummary["nevents"][local_idx]
            label = f"{cutSummary['cut_name'][local_idx]}({nevents})"
            plt.hist(ak.flatten(self.events[variable][mask]), bins=50, label=label, alpha=0.5)

        plt.xlabel(xlabel)
        plt.ylabel(f"Counts / {(xhigh-xlow)/nbins:.3f} {unit}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{self.save_path}/var_{variable}.png")
        plt.close()