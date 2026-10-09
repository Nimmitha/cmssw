"""Create the tnp_v9 CRAB task folders (mytagAndProbeV9): one production for the long-term efficiency study.

Two input streams, same analyzer, one submission round:
  muon   /Muon0, /Muon1 (IsoMu24 + non-isolated single-muon tags, tight, pt > 8; bitmask tag_trig): sa + trk trees, all events. Unbiased for tracking, muon
         reconstruction, ID, trigger (L1 DoubleMu + HLT) and vertex; boosted J/psi kinematics.
  pdmlm  /ParkingDoubleMuonLowMass0-7 (tag matched to HLT_Mu0_L1DoubleMu, tight, pt > 3): sa tree only (the L1 seed
         biases the track probes), every PRESCALE-th event. Tracking efficiency with the analysis kinematics.
Eras, processing versions, Global Tags and golden JSONs as the analysis production (production_reduced_size_v6):
the PDMLM dataset names are read from its configs, the Muon names follow the same processing strings and must be
checked on DAS (check_datasets.sh, written here) before submitting.
Folders: <stream>/<year>/<era><PD><version>, e.g. muon/2024/G01, pdmlm/2024/G51.
Usage: python3 make_tasks.py [--streams muon pdmlm] [--years 2023 2024 2025]
Rerunning overwrites configs and JSONs, never crab_* directories.
"""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROD = HERE.parent / "production_reduced_size_v6"
VERSION = "v9"
PRESCALE = 20
MUON_PDS = (0, 1)
# Muon stream: IsoMu24 plus non-isolated single-muon paths (bits 0..; a name the menu lacks just never fires): the
# isolated tag can veto J/psi whose probe track is close (HLT isolation cone ~ opening angle) -> compare by tag bit.
# bits 0..11 of tag_trig, in this order (names checked on a 2024G Muon0 file, 2026-10-09). Mu7p5_L2Mu2_Jpsi: standard
# tracking TnP path (second leg = L2 standalone muon: unbiased for tracking), prescaled. Double-muon / TkMu paths left
# out: their second leg biases the probe.
MUON_TAGS = ('"HLT_IsoMu24_v", "HLT_Mu50_v", "HLT_Mu8_v", "HLT_Mu17_v", "HLT_Mu3_PFJet40_v", "HLT_Mu12eta2p3_v", '
             '"HLT_Mu3_L1SingleMu5orSingleMu7_v", "HLT_Mu7p5_L2Mu2_Jpsi_v", "HLT_Mu15_v", "HLT_Mu19_v", "HLT_Mu20_v", '
             '"HLT_Mu27_v", "HLT_Mu55_v"')
TAG = {"muon": dict(tagPaths=MUON_TAGS, tagMinPt=8.0, prescale=1, fillSA=True, fillTrk=True),
       "pdmlm": dict(tagPaths='"HLT_Mu0_L1DoubleMu_v"', tagMinPt=3.0, prescale=PRESCALE, fillSA=True, fillTrk=False)}

CRAB = """from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '{request}'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniAODmuonsRootupler.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['{output}']
config.JobType.maxJobRuntimeMin = 240

config.Data.inputDBS = 'global'
config.Data.inputDataset = '{dataset}'
config.Data.lumiMask = '{mask}'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = {units}
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/nkarunar/tnp_v9/'
config.Data.publication = False

config.Site.storageSite = 'T3_US_FNALLPC'
"""

CMSSW = """import FWCore.ParameterSet.Config as cms
process = cms.Process("TnPv9")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '{gt}')

process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:test.root'))

process.rootuple = cms.EDAnalyzer('mytagAndProbeV9',
    muons = cms.InputTag("slimmedMuons"),
    pfCands = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    bits = cms.InputTag("TriggerResults", "", "HLT"),
    objects = cms.InputTag("slimmedPatTrigger"),
    l1Muons = cms.InputTag("gmtStage2Digis", "Muon"),
    prescales = cms.InputTag("patTrigger"),
    tagPaths = cms.vstring({tagPaths}),
    analysisPath = cms.string("HLT_Mu0_L1DoubleMu_v"),
    tagMinPt = cms.double({tagMinPt}),
    probeMinPt = cms.double(3.0),
    saMinPt = cms.double(2.0),
    maxEta = cms.double(2.4),
    massMin = cms.double(2.6),
    massMax = cms.double(3.6),
    saMassMin = cms.double(2.0),
    saMassMax = cms.double(4.5),
    prescale = cms.uint32({prescale}),
    fillSA = cms.bool({fillSA}),
    fillTrk = cms.bool({fillTrk}),
)
process.TFileService = cms.Service("TFileService", fileName = cms.string('{output}'))
process.p = cms.Path(process.rootuple)
"""


def analysis_tasks(year: str):
    """(folder, PD index, dataset, GT, mask path) of the analysis production, from its configs."""
    out = []
    for d in sorted((PROD / year).iterdir()):
        cfg = d / "crabConfig_TTree.py"
        if not cfg.exists():
            continue
        c = cfg.read_text()
        ds = re.search(r"inputDataset = '([^']+)'", c).group(1)
        mask = re.search(r"lumiMask = '([^']+)'", c).group(1)
        gt = re.search(r"GlobalTag, '([^']+)'", (d / "miniAODmuonsRootupler.py").read_text()).group(1)
        pd = int(re.match(r"/ParkingDoubleMuonLowMass(\d)/", ds).group(1))
        out.append((d.name, pd, ds, gt, d / mask))
    return out


def write(tdir: Path, stream: str, request: str, dataset: str, gt: str, mask: Path, units: int):
    tdir.mkdir(parents=True, exist_ok=True)
    output = f"{request}.root"
    (tdir / "crabConfig_TTree.py").write_text(CRAB.format(request=request, output=output, dataset=dataset,
                                                          mask=mask.name, units=units))
    t = {k: (str(v) if not isinstance(v, bool) else str(v)) for k, v in TAG[stream].items()}
    (tdir / "miniAODmuonsRootupler.py").write_text(CMSSW.format(gt=gt, output=output, **t))
    shutil.copy2(mask, tdir / mask.name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--streams", nargs="+", default=["muon", "pdmlm"])
    ap.add_argument("--years", nargs="+", default=["2023", "2024", "2025"])
    a = ap.parse_args()
    muon_datasets = []
    for year in a.years:
        tasks = analysis_tasks(year)
        if "pdmlm" in a.streams:
            for name, pd, ds, gt, mask in tasks:
                write(HERE / "pdmlm" / year / name, "pdmlm", f"PDMLM{pd}_{year}{name[0]}{name[2]}_tnp_{VERSION}", ds, gt,
                      mask, 200)
            print(f"pdmlm {year}: {len(tasks)} tasks")
        if "muon" in a.streams:
            n = 0
            for name, pd, ds, gt, mask in tasks:
                if pd != 0:
                    continue
                proc = ds.split("/")[2]            # e.g. Run2024G-PromptReco-v1
                for mpd in MUON_PDS:
                    mds = f"/Muon{mpd}/{proc}/MINIAOD"
                    muon_datasets.append(mds)
                    write(HERE / "muon" / year / f"{name[0]}{mpd}{name[2]}", "muon",
                          f"Muon{mpd}_{year}{name[0]}{name[2]}_tnp_{VERSION}", mds, gt, mask, 25)
                    n += 1
            print(f"muon {year}: {n} tasks")
    if muon_datasets:
        sh = HERE / "check_datasets.sh"
        sh.write_text("#!/bin/bash\n# Run on LPC/lxplus with a grid proxy: every Muon dataset must exist (else fix the\n"
                      "# processing string in the task's crabConfig_TTree.py; DAS: dataset=/Muon*/<era>*/MINIAOD).\n"
                      + "".join(f'n=$(dasgoclient -query="dataset={d}" | wc -l); [ "$n" = 1 ] && echo "ok      {d}" '
                                f'|| echo "MISSING {d}"\n' for d in muon_datasets))
        sh.chmod(0o755)
        print(f"{sh.name}: {len(muon_datasets)} Muon datasets to check")


if __name__ == "__main__":
    main()
