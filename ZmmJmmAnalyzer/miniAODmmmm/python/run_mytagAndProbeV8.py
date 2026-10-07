import FWCore.ParameterSet.Config as cms

process = cms.Process("JPsiTP")

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "141X_dataRun3_Prompt_v3", "")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        "file:/path/to/your_AOD.root"
    )
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("jpsi_aod_tagprobe.root")
)

process.tp = cms.EDAnalyzer(
    "mytagAndProbeV8",
    muons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    triggerBits = cms.InputTag("TriggerResults", "", "HLT"),

    # prefix match: "HLT_Mu0_L1DoubleMu_v"
    hltPaths = cms.vstring("HLT_Mu0_L1DoubleMu_v"),

    tagMinPt = cms.double(7.0),
    probeMinPt = cms.double(3.0),
    maxEta = cms.double(2.4),
    massMin = cms.double(2.6),
    massMax = cms.double(3.5),

    maxTrackMuonDR = cms.double(0.03),
    maxProbeDxy = cms.double(0.3),
    maxProbeDz = cms.double(20.0),
    minTrackerLayers = cms.int32(6),
    minValidPixelHits = cms.int32(1),
)

process.p = cms.Path(process.tp)