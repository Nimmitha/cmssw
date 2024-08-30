import os
import sys
import pandas as pd
import argparse

from utils.helperClasses import *
from utils.helperFunctions import *


def create_cutdict(events):
    cutdict = CutDict()
    cutdict.add_cut(0,
                    CutType(name="0",
                            long_name="Preselection",
                            mask=events['B_J1_mass'] > -1,
                            variables=['B_J1_mass'],
                            plot=False))

    cutdict.add_cut(1,
                    CutType(name="1",
                            long_name="Muon Trigger",
                            mask=events['Mu_TriggerPath'] == 1,
                            variables=['Mu_TriggerPath'],
                            plot=False))

    cutdict.add_cut(2,
                    CutType(name="2",
                            long_name="Soft Muons",
                            mask=(events['B_Mu1_soft'] == 1) & (events['B_Mu2_soft'] == 1) & (events['B_Mu3_soft'] == 1) & (events['B_Mu4_soft'] == 1),
                            variables=['B_J1_soft', 'B_J2_soft', 'B_J3_soft', 'B_J4_soft'],
                            plot=False))

    cutdict.add_cut(3,
                    CutType(name="3",
                            long_name="Muon pT",
                            mask=(events['B_Mu1_pt'] > 3) & (events['B_Mu2_pt'] > 3) & (events['B_Mu3_pt'] > 3) & (events['B_Mu4_pt'] > 3),
                            variables=['B_Mu1_pt', 'B_Mu2_pt', 'B_Mu3_pt', 'B_Mu4_pt'],
                            plot=True,
                            bins=100,
                            x_range=(0, 100),
                            labels=["Mu1_pt", "Mu2_pt", "Mu3_pt", "Mu4_pt"],
                            xlabel="pT (GeV)"))

    cutdict.add_cut(4,
                    CutType(name="4",
                            long_name="Detector acceptance",
                            mask=(abs(events['B_Mu1_eta']) < 2.4) & (abs(events['B_Mu2_eta']) < 2.4) & (abs(events['B_Mu3_eta']) < 2.4) & (abs(events['B_Mu4_eta']) < 2.4),
                            variables=['B_Mu1_eta', 'B_Mu2_eta', 'B_Mu3_eta', 'B_Mu4_eta'],
                            plot=True,
                            bins=100,
                            x_range=(-3, 3),
                            labels=["Mu1_eta", "Mu2_eta", "Mu3_eta", "Mu4_eta"],
                            xlabel="Eta"))

    cutdict.add_cut(5,
                    CutType(name="5",
                            long_name="Any dimuon pair vertex",
                            mask=((events['B_J1_VtxProb'] > 0.01) & (events['B_J2_VtxProb'] > 0.01)) | ((events['B_J3_VtxProb'] > 0.01) & (events['B_J4_VtxProb'] > 0.01)),
                            variables=['B_J1_VtxProb', 'B_J2_VtxProb', 'B_J3_VtxProb', 'B_J4_VtxProb'],
                            plot=True,
                            bins=100,
                            x_range=(0, 1),
                            labels=["J1_VtxProb", "J2_VtxProb", "J3_VtxProb", "J4_VtxProb"],
                            xlabel="Vtx Prob"))

    cutdict.add_cut(6,
                    CutType(name="6",
                            long_name="Four muon vertex",
                            mask=(events['FourL_VtxProb'] > 0.01),
                            variables=['FourL_VtxProb'],
                            plot=True,
                            bins=100,
                            x_range=(0, 1),
                            labels=["FourL_VtxProb"],
                            xlabel="Vtx Prob"))

    cutdict.add_cut(7,
                    CutType(name="7",
                            long_name="Candidate Combine",
                            mask=events['UpsMass'] > 0,
                            variables=['UpsMass'],
                            plot=False))

    cutdict.add_cut(8,
                    CutType(name="8",
                            long_name="Dimuon Vertex Prob",
                            mask=(events['Ups_VtxProb1'] > 0.01) & (events['Ups_VtxProb2'] > 0.01),
                            variables=['Ups_VtxProb1', 'Ups_VtxProb2'],
                            plot=True,
                            bins=100,
                            x_range=(0, 1),
                            labels=["VtxProb1", "VtxProb2"],
                            xlabel="Vtx Prob"))

    cutdict.add_cut(9,
                    CutType(name="9",
                            long_name="Dimuon pT",
                            mask=(events['Ups_Pt1'] > 5) & (events['Ups_Pt2'] > 5),
                            variables=['Ups_Pt1', 'Ups_Pt2'],
                            plot=True,
                            bins=100,
                            x_range=(0, 20),
                            labels=["Ups_Pt1", "Ups_Pt2"],
                            xlabel="pT (GeV)"))

    cutdict.add_cut(10,
                    CutType(name="10",
                            long_name="Dimuon mass",
                            mask=(events['Ups1_mass'] > 2.8) & (events['Ups1_mass'] < 3.4) & (events['Ups2_mass'] > 70) & (events['Ups2_mass'] < 110),
                            variables=['Ups1_mass', 'Ups2_mass'],
                            plot=True,
                            bins=100,
                            x_range=(0, 20),
                            labels=["Ups1_mass", "Ups2_mass"],
                            xlabel="mass (GeV)"))

    cutdict.add_cut(11,
                    CutType(name="11",
                            long_name="Four muon pT",
                            mask=(events['FourL_pt'] > 5),
                            variables=['FourL_pt'],
                            plot=True,
                            bins=100,
                            x_range=(0, 100),
                            labels=["FourL_pt"],
                            xlabel="pT (GeV)"))

    cutdict.add_cut(12,
                    CutType(name="12",
                            long_name="Four muon mass signal",
                            mask=(events['FourL_mass'] > 112) & (events['FourL_mass'] < 162),
                            variables=['FourL_mass'],
                            plot=True,
                            bins=50,
                            x_range=(112, 162),
                            labels=["FourL_mass"],
                            xlabel="mass (GeV)"))

    cutdict.add_cut(13,
                    CutType(name="13",
                            long_name="Four muon mass background",
                            mask=(events['FourL_mass'] < 120) | (events['FourL_mass'] > 130),
                            variables=['FourL_mass'],
                            plot=True,
                            bins=50,
                            x_range=(0, 300),
                            labels=["FourL_mass"],
                            xlabel="mass (GeV)"))

    return cutdict


def run_analysis(config, isMC):
    if isMC:
        print("Running analysis on MC events...")
        save_path = os.path.join(config['savepath_base'], 'mc')
        file_path = config['mc_path']
    else:
        print("Running analysis on data events...")
        save_path = os.path.join(config['savepath_base'], 'data')
        file_path = config['data_path']

    os.makedirs(save_path, exist_ok=True)

    print(f"Loading events from {file_path}...")
    events = load_events(file_path, BaseSchema, config['entry_stop'])
    if events is None:
        sys.exit(1)

    print("Computing events with selected columns...")
    events = events[config['columns']].compute()

    print("Combining dimuon candidates...")
    events = combine_dimuon_candidates(events, J_lower=2.8, J_upper=3.4, z_lower=70, z_upper=110)

    print("\nThese are the available cuts:")
    cutdict = create_cutdict(events)
    for key, item in cutdict.cutdict.items():
        print(key, item)

    myorder = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    print(f"\nRunning cut analysis in {myorder}")

    cut_analysis = CutAnalysis(events, cutdict)
    cutSummary = cut_analysis.prepare_masks(myorder)

    cutSummarydf = {k: v for k, v in cutSummary.items() if k not in ['mask', 'cutobj']}
    print(pd.DataFrame(cutSummarydf))

    plotter = Plotter(events, cutdict, save_path)
    plotter.plot_preselection()
    plotter.plot_summary(cutSummary)

    cuts_to_show = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    plotter.plot_single_variable(cuts_to_show, cutSummary)

    print("Analysis complete!!\n")
    return events, cutdict, cutSummary


def main(config_path, run_data, run_mc):
    config = load_config(config_path)

    if run_data:
        events_data, cutdict, cutSummary = run_analysis(config, isMC=False)
    if run_mc:
        events_mc, cutdict, cutSummary = run_analysis(config, isMC=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run particle physics analysis on data and/or MC events.")
    parser.add_argument("--config", default="config.yaml", help="Path to the configuration file")
    parser.add_argument("--data", action="store_true", help="Run analysis on data events")
    parser.add_argument("--mc", action="store_true", help="Run analysis on MC events")
    args = parser.parse_args()

    if not args.data and not args.mc:
        print("Error: You must specify at least one of --data or --mc")
        sys.exit(1)

    main(config_path=args.config, run_data=args.data, run_mc=args.mc)
