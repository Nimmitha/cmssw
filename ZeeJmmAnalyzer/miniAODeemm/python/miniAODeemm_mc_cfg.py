import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")

options.setDefault(
    "inputFiles",
    "file:/uscms/home/wkarunar/nobackup/datasets/mc/miniAOD/run2/zeejmm/2018_ss/MiniAOD/MiniAOD_1.root"
)

options.setDefault(
    "outputFile",
    "selection/2018_ss/zeejmm_mc_2018_v1_1.root"
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

# Run 2 UL MC GlobalTags. Uncomment the one matching the MC campaign/year.
# 2016 MC
# process.GlobalTag = GlobalTag(process.GlobalTag, "106X_mcRun2_asymptotic_v17", "")

# 2017 MC
# process.GlobalTag = GlobalTag(process.GlobalTag, "106X_mc2017_realistic_v10", "")

# 2018 MC
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
    "miniAODeemm",
    electrons = cms.InputTag("slimmedElectrons"),
    muons = cms.InputTag("slimmedMuons"),
    packedCandidates = cms.InputTag("packedPFCandidates"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    bits = cms.InputTag("TriggerResults", "", "HLT"),
    objects = cms.InputTag("slimmedPatTrigger"),

    # Pick the trigger for the MC campaign/year by commenting/uncommenting one line.
    # 2016 MC
    # ElectronTrigger = cms.string("HLT_Ele27_WPTight_Gsf_v"),

    # 2017 MC: the analyzer applies the hltEGL1SingleEGOrFilter treatment automatically.
    # ElectronTrigger = cms.string("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"),

    # 2018 MC
    ElectronTrigger = cms.string("HLT_Ele32_WPTight_Gsf_v"),

    isMC = cms.bool(True),
    requireTrigger = cms.bool(True),
    requireTriggerMatch = cms.bool(False),
    keepEmptyEvents = cms.bool(False),
    requireLooseElectronID = cms.bool(False),
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile),
)

process.p = cms.Path(process.rootuple)
