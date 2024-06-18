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

#process.load('RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16_ID_ISO_cff')
# process.load("Configuration.Geometry.GeometryIdeal_cff")

# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
# from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc') # Tried: 106X_dataRun2_v32

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v15')
# process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_v2')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data')

process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  allowUnscheduled = cms.untracked.bool(True),
  SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:test.root')
)

process.rootuple = cms.EDAnalyzer('miniAODtrigRate',
                          muons = cms.InputTag("slimmedMuons"),
                          bits = cms.InputTag("TriggerResults::HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          MuonTrigger = cms.string("HLT_IsoMu24_v"),
                          isMC = cms.bool(False),
                          )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('Muon_Run2022C_v1_Data.root'),
)

process.p = cms.Path(process.rootuple)

