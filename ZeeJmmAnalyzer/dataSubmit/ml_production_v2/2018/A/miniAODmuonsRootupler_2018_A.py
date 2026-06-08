import FWCore.ParameterSet.Config as cms

process = cms.Process("rootuple")

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, "106X_dataRun2_v37", "")

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True,
                       runVID=True,
                       era="2018-UL")

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
    fileNames = cms.untracked.vstring(
        "file:test.root"
    ),
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("EGamma_Run2018A_UL_v2_v1_Data.root"),
)

process.rootuple = cms.EDAnalyzer(
    "miniAODeemm",
    electrons = cms.InputTag("slimmedElectrons"),
    muons = cms.InputTag("slimmedMuons"),
    packedCandidates = cms.InputTag("packedPFCandidates"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    bits = cms.InputTag("TriggerResults", "", "HLT"),
    objects = cms.InputTag("slimmedPatTrigger"),
    ElectronTrigger = cms.string("HLT_Ele32_WPTight_Gsf_v"),
    isMC = cms.bool(False),
    requireTrigger = cms.bool(True),
    requireTriggerMatch = cms.bool(False),
    keepEmptyEvents = cms.bool(False),
    requireLooseElectronID = cms.bool(False),
)

process.p = cms.Path(process.egammaPostRecoSeq + process.rootuple)
