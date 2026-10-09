import FWCore.ParameterSet.Config as cms
process = cms.Process("TnPv9")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_dataRun3_Prompt_v3')

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
    tagPaths = cms.vstring("HLT_IsoMu24_v", "HLT_Mu50_v", "HLT_Mu8_v", "HLT_Mu17_v", "HLT_Mu3_PFJet40_v", "HLT_Mu12eta2p3_v", "HLT_Mu3_L1SingleMu5orSingleMu7_v", "HLT_Mu7p5_L2Mu2_Jpsi_v", "HLT_Mu15_v", "HLT_Mu19_v", "HLT_Mu20_v", "HLT_Mu27_v", "HLT_Mu55_v"),
    analysisPath = cms.string("HLT_Mu0_L1DoubleMu_v"),
    tagMinPt = cms.double(8.0),
    probeMinPt = cms.double(3.0),
    saMinPt = cms.double(2.0),
    maxEta = cms.double(2.4),
    massMin = cms.double(2.6),
    massMax = cms.double(3.6),
    saMassMin = cms.double(2.0),
    saMassMax = cms.double(4.5),
    prescale = cms.uint32(1),
    fillSA = cms.bool(True),
    fillTrk = cms.bool(True),
)
process.TFileService = cms.Service("TFileService", fileName = cms.string('Muon1_mm_2025G1_tnp_v9.root'))
process.p = cms.Path(process.rootuple)
