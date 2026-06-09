import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")

options.setDefault(
    "inputFiles",
    "file:/uscms/home/wkarunar/nobackup/datasets/mc/miniAOD/run2/zeejmm_2018/zeejmm_ss/MiniAOD/MiniAOD_1.root"
)

options.setDefault(
    "outputFile",
    "preselection/2018_ss/zeejmm_mc_2018_v2_1.root"
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

# 2018 MC GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "106X_upgrade2018_realistic_v16_L1v1", "")

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True,
                       runVID=True,
                       era='2018-UL') # 2018-UL  2017-UL  2016postVFP-UL  2016preVFP-UL


process.MessageLogger.cerr.FwkReport.reportEvery = 1000

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

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile),
)

process.rootuple = cms.EDAnalyzer(
    "miniAODeemm",
    electrons = cms.InputTag("slimmedElectrons"),
    muons = cms.InputTag("slimmedMuons"),
    packedCandidates = cms.InputTag("packedPFCandidates"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    bits = cms.InputTag("TriggerResults", "", "HLT"),
    objects = cms.InputTag("slimmedPatTrigger"),
    ElectronTrigger = cms.string("HLT_Ele32_WPTight_Gsf_v"), # 2018 MC
    isMC = cms.bool(True),
    requireTrigger = cms.bool(True),
    requireTriggerMatch = cms.bool(False),
    keepEmptyEvents = cms.bool(False),
    requireLooseElectronID = cms.bool(False),
)

process.p = cms.Path(process.egammaPostRecoSeq + process.rootuple)
