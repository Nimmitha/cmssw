import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackingEff")

# --- Standard Sequences ---
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# --- Global Tag ---
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_dataRun3_Prompt_v3')

# --- Message Logger ---
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# --- Source (IMPORTANT: Use SingleMuon dataset) ---
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # Example SingleMuon file
        'file:miniAODs/00b7fe9f-e5a2-42c5-92e6-5992e2238db2_MINI6_2024_Muon_C0.root'
    )
)

# --- Analyzer Configuration ---
# Name must match the DEFINE_FWK_MODULE in the .cc file
process.rootuple = cms.EDAnalyzer('mytagAndProbeV6',
    muons = cms.InputTag("slimmedMuons"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    bits = cms.InputTag("TriggerResults::HLT"),
    objects = cms.InputTag("slimmedPatTrigger"),
    # prescales = cms.InputTag("patTrigger"),
    # isMC = cms.bool(False)
)

# --- Output File ---
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('preselection/tracking_efficiency.root'),
)

process.p = cms.Path(process.rootuple)