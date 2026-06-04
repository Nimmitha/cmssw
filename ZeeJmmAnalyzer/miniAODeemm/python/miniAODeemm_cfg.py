import FWCore.ParameterSet.Config as cms

process = cms.Process("MINIAODEEMM")

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

# Data UL GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "106X_dataRun2_v37", "")

# 2018 MC GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, "106X_upgrade2018_realistic_v16_L1v1", "")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True),
    allowUnscheduled=cms.untracked.bool(True),
)

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(-1)
)

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(
        # Replace with your Run-2 SingleElectron/EGamma MiniAOD files.
        "file:input.root"
    ),
)

process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string("miniAODeemm.root"),
)

process.miniAODeemm = cms.EDAnalyzer(
    "miniAODeemm",
    electrons=cms.InputTag("slimmedElectrons"),
    muons=cms.InputTag("slimmedMuons"),
    packedCandidates=cms.InputTag("packedPFCandidates"),
    primaryVertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    bits=cms.InputTag("TriggerResults", "", "HLT"),
    objects=cms.InputTag("slimmedPatTrigger"),

    # Pick the trigger for the data-taking year by commenting/uncommenting one line.
    # 2016
    # ElectronTrigger=cms.string("HLT_Ele27_WPTight_Gsf_v"),
    # 2017: this analyzer applies the hltEGL1SingleEGOrFilter treatment automatically.
    # ElectronTrigger=cms.string("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"),
    # 2018
    ElectronTrigger=cms.string("HLT_Ele32_WPTight_Gsf_v"),

    isMC=cms.bool(False),
    requireTrigger=cms.bool(True),
    requireTriggerMatch=cms.bool(False),
    keepEmptyEvents=cms.bool(False),
    requireLooseElectronID=cms.bool(False),
)

process.p = cms.Path(process.miniAODeemm)
