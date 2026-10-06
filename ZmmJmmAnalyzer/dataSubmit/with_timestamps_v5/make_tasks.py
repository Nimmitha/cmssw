"""Create the with_timestamps_v5 CRAB task folders from the emittance-scan lumi masks.

Masks: dimuon_vdm_calibration/workflows/00_scan_inventory/lumimasks/lumi_mask_emit_<year>.json
(selected scans, scan LS +- 120 s). Copied to <year>/lumi_mask_emit_<year>.json here, and into
every task folder (crab needs it next to the config).

One task per dataset (ParkingDoubleMuonLowMass0-7 x era x processing version) that holds runs of
the mask; folder <era><PD><version>, e.g. 2023/C03 = PD0, Run2023C-22Sep2023_v3-v1. Datasets and
Global Tags as in production_reduced_size_v6; run ranges per dataset from the processedLumis of
that production. A mask run outside every range stops the script.

Usage: python3 make_tasks.py [--mask-dir DIR] [--years 2023 2024 2025]
Rerunning overwrites the configs and masks of existing folders, but not crab_* directories.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
MASK_DIR = Path("/home/nimmitha/research/jpsi_lumi/dimuon_vdm_calibration/workflows/00_scan_inventory/lumimasks")
VERSION = "v5"

GT_2023 = "130X_dataRun3_PromptAnalysis_v1"
GT_PROMPT = "141X_dataRun3_Prompt_v3"
# (era, version digit, dataset name after /ParkingDoubleMuonLowMass<n>/, first run, last run, GT)
DATASETS = {
    "2023": [
        ("C", 1, "Run2023C-22Sep2023_v1-v1", 367095, 367515, GT_2023),
        ("C", 2, "Run2023C-22Sep2023_v2-v1", 367516, 367619, GT_2023),
        ("C", 3, "Run2023C-22Sep2023_v3-v1", 367620, 367758, GT_2023),
        ("C", 4, "Run2023C-22Sep2023_v4-v1", 367770, 368823, GT_2023),
        ("D", 1, "Run2023D-22Sep2023_v1-v1", 369927, 370580, GT_2023),
        ("D", 2, "Run2023D-22Sep2023_v2-v1", 370667, 370790, GT_2023),
    ],
    "2024": [
        ("C", 1, "Run2024C-PromptReco-v1", 379416, 380238, GT_PROMPT),
        ("D", 1, "Run2024D-PromptReco-v1", 380306, 380947, GT_PROMPT),
        ("E", 1, "Run2024E-PromptReco-v1", 380963, 381380, GT_PROMPT),
        ("E", 2, "Run2024E-PromptReco-v2", 381384, 381544, GT_PROMPT),
        ("F", 1, "Run2024F-PromptReco-v1", 382229, 383779, GT_PROMPT),
        ("G", 1, "Run2024G-PromptReco-v1", 383811, 385801, GT_PROMPT),
        ("H", 1, "Run2024H-PromptReco-v1", 385836, 386319, GT_PROMPT),
        ("I", 1, "Run2024I-PromptReco-v1", 386478, 386693, GT_PROMPT),
        ("I", 2, "Run2024I-PromptReco-v2", 386694, 386951, GT_PROMPT),
    ],
    "2025": [
        ("C", 1, "Run2025C-PromptReco-v1", 392293, 393087, GT_PROMPT),
        ("C", 2, "Run2025C-PromptReco-v2", 393111, 393461, GT_PROMPT),
        ("D", 1, "Run2025D-PromptReco-v1", 394637, 395948, GT_PROMPT),
        ("E", 1, "Run2025E-PromptReco-v1", 395982, 396422, GT_PROMPT),
        ("F", 1, "Run2025F-PromptReco-v1", 396733, 397596, GT_PROMPT),
        ("F", 2, "Run2025F-PromptReco-v2", 397619, 397817, GT_PROMPT),
        ("G", 1, "Run2025G-PromptReco-v1", 398011, 398860, GT_PROMPT),
    ],
}

CRAB = """from CRABClient.UserUtilities import config
config = config()

# user specific generic parameters
config.General.requestName = '{request}'   # Used as the task/Project directory name
config.General.transferOutputs = True                               # Transfer output files to the storage site
config.General.transferLogs = False

# job type and related configurables
config.JobType.pluginName = 'Analysis'                              # Specify: analysis or MC generation
config.JobType.psetName = 'miniAODmuonsRootupler.py'           # parameter-set config file
config.JobType.allowUndistributedCMSSW = True                       # Allow CMSSW release possibly not available at sites
config.JobType.outputFiles = ['{output}']   # List of output files that needs to be collected
config.JobType.maxJobRuntimeMin = 60                                  # Maximum job runtime in minutes
# config.JobType.maxMemoryMB = 1000

# data to be analyzed
config.Data.inputDBS = 'global'
config.Data.inputDataset = '{dataset}'               # Name of the dataset
config.Data.lumiMask = '{mask}' # Lumi-section filter
config.Data.splitting = 'LumiBased'                                                         # Split the task based on
config.Data.unitsPerJob = 50                                                                # Number of splitted units per job
config.Data.totalUnits = -1                                                                 # Number of untis to analyze
config.Data.outLFNDirBase = '/store/user/nkarunar/emit/'
config.Data.publication = False                                                             # Whether to publish the EDM output files in DBS

# Grid site parameters
config.Site.storageSite = 'T3_US_FNALLPC'         # Place to copy the output files
"""

CMSSW = """import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '{gt}')

process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  allowUnscheduled = cms.untracked.bool(True),
  # SkipEvent = cms.untracked.vstring('ProductNotFound')
  )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:test.root')
)

process.rootuple = cms.EDAnalyzer('miniAODmmmm',
                          muons = cms.InputTag("slimmedMuons"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bits = cms.InputTag("TriggerResults::HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          prescales = cms.InputTag("patTrigger"),
                          pruned = cms.InputTag("prunedGenParticles"),
                          MuonTrigger = cms.string("HLT_Mu0_L1DoubleMu_v"),
                          isMC = cms.bool(False),
                          )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('{output}'),
)

process.p = cms.Path(process.rootuple)
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mask-dir", type=Path, default=MASK_DIR)
    ap.add_argument("--years", nargs="+", default=list(DATASETS))
    args = ap.parse_args()
    for year in args.years:
        mask_name = f"lumi_mask_emit_{year}.json"
        mask = json.loads((args.mask_dir / mask_name).read_text())
        runs = sorted(int(r) for r in mask)
        ydir = HERE / year
        ydir.mkdir(exist_ok=True)
        shutil.copy2(args.mask_dir / mask_name, ydir / mask_name)
        used = {}
        for run in runs:
            ds = [d for d in DATASETS[year] if d[3] <= run <= d[4]]
            if not ds:
                raise SystemExit(f"{year}: run {run} is outside every dataset range; extend DATASETS")
            used.setdefault(ds[0], []).append(run)
        for (era, ver, name, *_rest, gt), ds_runs in used.items():
            for pd in range(8):
                tdir = ydir / f"{era}{pd}{ver}"
                tdir.mkdir(exist_ok=True)
                output = f"PDMLM{pd}_Run{year}{era}{ver}_Data.root"
                (tdir / "crabConfig_TTree.py").write_text(CRAB.format(
                    request=f"PDMLM{pd}_mm_{year}{era}{ver}_emit_{VERSION}", output=output,
                    dataset=f"/ParkingDoubleMuonLowMass{pd}/{name}/MINIAOD", mask=mask_name))
                (tdir / "miniAODmuonsRootupler.py").write_text(CMSSW.format(gt=gt, output=output))
                shutil.copy2(ydir / mask_name, tdir / mask_name)
            print(f"{year} {era}{ver} {name}: {len(ds_runs)} runs -> {era}0{ver} ... {era}7{ver}")
        n_ls = sum(b - a + 1 for v in mask.values() for a, b in v)
        print(f"{year}: {len(runs)} runs, {n_ls} LS, {len(used)} datasets x 8 PDs = {8 * len(used)} tasks")


if __name__ == "__main__":
    main()
