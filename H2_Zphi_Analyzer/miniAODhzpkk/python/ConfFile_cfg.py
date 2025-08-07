import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff') # Jesse

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
            'file:miniAOD/MiniAOD1000.root'
                )
                            )

process.demo = cms.EDAnalyzer('miniAODhzpkk',
   muons      = cms.untracked.InputTag('slimmedMuons'),
   trigbits      = cms.untracked.InputTag('TriggerResults','','HLT')
                              )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('preselection/output.root'),
)

process.p = cms.Path(process.demo)
