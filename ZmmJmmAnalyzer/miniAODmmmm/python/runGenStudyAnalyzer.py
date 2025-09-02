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
# process.GlobalTag = GlobalTag(process.GlobalTag, '133X_mcRun3_2024_realistic_v8') # JPsiToMuMu MC
process.GlobalTag = GlobalTag(process.GlobalTag, '142X_mcRun3_2025_realistic_v7') # JPsiToMuMu MC25


# process.MessageLogger.cerr.FwkReport.reportEvery = 500


process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout'),
    categories = cms.untracked.vstring('GenStudyAnalyzer'),

    cout = cms.untracked.PSet(
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
    fileNames = cms.untracked.vstring(       
        # 'file:miniAODs/1d0b3161-17d1-450d-bca4-1898c5268117_Run3WinterMC.root'
        '/store/mc/Run3Winter25MiniAOD/JPsiToMuMu_PT-0to100_pythia8-gun/MINIAODSIM/FlatPU0to120_MiniAODv6_150X_mcRun3_2025_realistic_v6-v2/120000/6b784ff3-0f02-478c-9ba6-e9e8f4c8ca7a.root'
 )
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
  fileName = cms.string('preselection/testout_2025WinterMC.root'),
)

process.p = cms.Path(process.rootuple)

