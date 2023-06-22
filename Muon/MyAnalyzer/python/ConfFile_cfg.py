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
            # 'file:/eos/home-n/nkarunar/2023Analysis/datasets/2022C/slimMiniAOD_data_MuEle.root'
            'file:/eos/home-n/nkarunar/2023Analysis/datasets/2018_MC/Zuu_2018_MINIAOD.root'
                )
                            )

process.demo = cms.EDAnalyzer('MyAnalyzer',
   tracks    = cms.untracked.InputTag('generalTracks'),
   trackPtMin = cms.double(0.3)
                              )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('analyzer_output.root'),
)

# process.Tracer = cms.Service("Tracer")

process.p = cms.Path(process.demo)