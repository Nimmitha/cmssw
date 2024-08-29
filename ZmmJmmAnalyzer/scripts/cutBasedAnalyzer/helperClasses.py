import matplotlib.pyplot as plt
import awkward as ak
import mplhep as hep
hep.style.use(hep.style.CMS)


def process_ups_candidates(data):
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

    J_lower = 2.8
    J_upper = 3.4
    z_lower = 70.0
    z_upper = 110.0
    
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
    
    
class CutDict:
    def __init__(self):
        self.cutdict = {}

    def add_cut(self, cut_id, cut):
        self.cutdict[cut_id] = cut

    def get_cut(self, cut_id):
        return self.cutdict.get(cut_id)
    
    def __str__(self):
        return "\n".join([f"{cut_id}: {cut}" for cut_id, cut in self.cutdict.items()])


class CutAnalysis:
    def __init__(self, events, cutdict):
        self.events = events
        self.cutdict = cutdict

    def get_stats(self, data):
        array = data.B_J1_mass
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
            # print(f"Processing cut {cut_id}: {cut.long_name}")
            # print(f"Applying mask: {cut.mask}")
            # if type(cut.mask) == str and cut.mask == 'Special_Combine':
            #     print("Processing Upsilon candidates")
            #     self.events = process_ups_candidates(self.events)
            #     cut.mask = (self.events.UpsMass > 0)
            #     print(f"Upsilon candidates processed. {ak.sum(cut.mask)} events")
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
            xlow, xhigh = cut.x_range
            labels = cut.labels
            xlabel = cut.xlabel
            variables = cut.variables
            unit = 'GeV' if 'eta' not in xlabel else ''
            
            plt.figure()
            for i, variable in enumerate(variables):
                plt.hist(ak.flatten(self.events[variable][cut.mask]), bins=nbins, range=(xlow, xhigh), label=labels[i], alpha=0.5)
            plt.xlabel(cut.xlabel)
            plt.ylabel(f"Counts / {(xhigh-xlow)/nbins:.3f} {unit}")
            plt.title(f"{cut.long_name}")
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"{self.save_path}/pre_{cut.long_name}.png")
            plt.close()

    def plot_summary(self, cutSummary):
        for cut_id in cutSummary["cut_id"]:
            self.plot_summary_at(cut_id, cutSummary)
            

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


    def plot_single_variable(self, cutsToShow, cutSummary):
        for cut_id in self.cutdict.cutdict:
            if not self.cutdict.cutdict[cut_id].plot:
                print(f"Cut {cut_id} does not have a plot.")
                continue
            for variable in self.cutdict.cutdict[cut_id].variables:
                self.plot_single_variable_for(variable, cutsToShow, cutSummary)

    def plot_single_variable_for(self, variable, cutsToShow, cutSummary):
        # find the cut which has the variable
        cut_id = None
        for cut_id in self.cutdict.cutdict:
            if variable in self.cutdict.cutdict[cut_id].variables and self.cutdict.cutdict[cut_id].plot:
                break

        if cut_id is None:
            print(f"Info for variable {variable} not found in any cutdict.")
            return
        
        if not self.cutdict.cutdict[cut_id].plot:
            print(f"Again, Cut {cut_id} does not have a plot.")
            return
        
        print(f"Plotting variable {variable} at steps: {cutsToShow}")
        
        nbins = self.cutdict.cutdict[cut_id].bins
        xlow, xhigh = self.cutdict.cutdict[cut_id].x_range
        xlabel = self.cutdict.cutdict[cut_id].xlabel
        unit = 'GeV' if 'eta' not in variable else ''

        cut_names = [self.cutdict.get_cut(cut_id).long_name for cut_id in cutsToShow]
        # print(f"Plotting variable {variable} at: {cut_names}")

        plt.figure()
        for cut_id in cutsToShow:
            local_idx = cutSummary["cut_id"].index(cut_id)
            mask = cutSummary["mask"][local_idx]
            nevents = cutSummary["nevents"][local_idx]
            label = f"{cutSummary['cut_name'][local_idx]}({nevents})"
            plt.hist(ak.flatten(self.events[variable][mask]), bins=50, range=(xlow, xhigh), label=label, alpha=0.5)

        plt.xlabel(variable)
        plt.ylabel(f"Counts / {(xhigh-xlow)/nbins:.3f} {unit}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{self.save_path}/var_{variable}.png")
        plt.close()