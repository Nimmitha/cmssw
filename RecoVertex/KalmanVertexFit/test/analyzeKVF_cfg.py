import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.Services_cff")

# process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff") # Jesse

# process.load("Configuration.StandardSequences.FakeConditions_cff")

# process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff') # Jesse

process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByChi2_cfi")

process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #    '/store/relval/CMSSW_2_1_8/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v1/0002/04983078-9082-DD11-BB8C-0019DB2F3F9B.root'
    'file:/eos/home-n/nkarunar/2023Analysis/datasets/2022C/7afde7ee-e709-46e9-ae45-7725c28af2bb.root'   
       )
)

process.simpleVertexAnalysis = cms.EDAnalyzer("KVFTest",
    KVFParameters = cms.PSet(
        maxDistance = cms.double(0.01),
        maxNbrOfIterations = cms.int32(10)
    ),
    outputFile = cms.untracked.string('simpleVertexAnalyzer.root'),
    TrackLabel = cms.string('generalTracks')
)

process.p = cms.Path(process.simpleVertexAnalysis)


