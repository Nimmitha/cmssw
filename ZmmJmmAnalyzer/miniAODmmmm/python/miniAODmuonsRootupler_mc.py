import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")

options.setDefault(
    "inputFiles",
    "file:/uscms/home/wkarunar/nobackup/datasets/mc/miniAOD/run2/zmmjmm/2018_ss/MiniAOD/MiniAOD_1.root"
)

options.setDefault(
    "outputFile",
    "selection/2018_ss/zmmjmm_mc_2018_v1_1.root"
)

options.parseArguments()


process = cms.Process("Rootuple")

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
# process.GlobalTag = GlobalTag(process.GlobalTag, "106X_dataRun2_v37", "")

# 2018 MC GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "106X_upgrade2018_realistic_v16_L1v1", "")

process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.rootuple = cms.EDAnalyzer(
    "miniAODmmmm",
    dimuons = cms.InputTag("slimmedMuons"),
    Trak = cms.InputTag("packedPFCandidates"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    bits = cms.InputTag("TriggerResults::HLT"),

    # 2016/2018 data
    MuonTrigger = cms.string("HLT_IsoMu24_v"),

    # For 2017 data use:
    # MuonTrigger = cms.string("HLT_IsoMu27_v"),

    requireTrigger = cms.bool(True),
    keepEmptyEvents = cms.bool(False),
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile),
)

process.p = cms.Path(process.rootuple)