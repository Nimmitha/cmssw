import FWCore.ParameterSet.Config as cms
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
process.GlobalTag = GlobalTag(process.GlobalTag, '133X_mcRun3_2024_realistic_v8') # JPsiToMuMu MC

# process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('messages'),
    categories = cms.untracked.vstring('GenStudyAnalyzer'),

    messages = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),  # or 'WARNING' if you want only warnings+
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)          # suppress all other categories
        ),
        GenStudyAnalyzer = cms.untracked.PSet(
            limit = cms.untracked.int32(100000)     # allow all messages from this category
        )
    )
)

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  allowUnscheduled = cms.untracked.bool(True),
  # SkipEvent = cms.untracked.vstring('ProductNotFound')
  )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:test.root')
)

process.rootuple = cms.EDAnalyzer('GenStudyAnalyzer',
                          muons = cms.InputTag("slimmedMuons"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bits = cms.InputTag("TriggerResults::HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          prescales = cms.InputTag("patTrigger"),
                          pruned = cms.InputTag("prunedGenParticles"),
                          MuonTrigger = cms.string("HLT_Mu0_L1DoubleMu_v"),
                          isMC = cms.bool(True),
                          )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('JPsiToMuMu_PT_0to100_Winter24_MC.root'),
)

process.p = cms.Path(process.rootuple)

