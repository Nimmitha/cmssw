# import FWCore.ParameterSet.Config as cms

# process = cms.Process("MyProducerDemo")

# # Replace 'source' with your desired input source (e.g., input file)
# process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring('file:/eos/home-n/nkarunar/2023Analysis/datasets/2018_MC/Zuu_2018_MINIAOD.root')
# )

# # Define your producer module and provide input parameters (e.g., prunedGenParticles)
# process.myProducer = cms.EDProducer('MyProducer',
#     prunedGenParticles = cms.InputTag('prunedGenParticles'),  # Adjust the input tag if necessary
#     slimmedMuons = cms.InputTag('slimmedMuons')
# )

# # Output module to store the produced collection (if needed)
# process.out = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string('output.root'),
#     # outp
# )

# # Adjust the number of events to process as needed
# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# # Adjust the frequency of printouts as needed
# process.MessageLogger = cms.Service("MessageLogger",
#     destinations = cms.untracked.vstring('cout'),
#     cout = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
# )

# # Endpath to define the sequence of modules to run
# process.myPath = cms.EndPath(process.myProducer * process.out)

# # Uncomment this line if you want to write the output to a file
# # process.e = cms.EndPath(process.out)

# # Uncomment this line if you want to run the output to a file without producing any output
# # process.p = cms.EndPath()






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

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi") # ParticleListDrawer
process.load("PhysicsTools.HepMCCandAlgos.allMuonsGenParticlesMatch_cfi") # MCTruthDeltaRMatcherNew 

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
            # 'file:/eos/home-n/nkarunar/2023Analysis/datasets/2022C/slimMiniAOD_data_MuEle.root'
            # 'file:/eos/home-n/nkarunar/2023Analysis/datasets/2018_MC/Zuu_2018_MINIAOD.root'
            'file:output.root'
                )
                            )

process.demo = cms.EDAnalyzer('MyAnalyzer',
   muons      = cms.untracked.InputTag('slimmedMuons'),
   trigbits      = cms.untracked.InputTag('TriggerResults','','HLT'),
   tracks    = cms.untracked.InputTag('generalTracks'),
   trackPtMin = cms.double(0.3)
                              )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('analyzer_output.root'),
)

# process.printList = cms.EDAnalyzer("ParticleListDrawer",
#   maxEventsToPrint = cms.untracked.int32(1),
#   printVertex = cms.untracked.bool(False),
#   src = cms.InputTag("prunedGenParticles")
# )

# process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
#   src = cms.InputTag("prunedGenParticles"),
#   printP4 = cms.untracked.bool(False),
#   printPtEtaPhi = cms.untracked.bool(False),
#   printVertex = cms.untracked.bool(False),
#   printStatus = cms.untracked.bool(False),
#   printIndex = cms.untracked.bool(False),
#   status = cms.untracked.vint32( 3 )
# )

process.selectedMuonsGenParticlesMatch = cms.EDProducer("MCTruthDeltaRMatcher",
  src = cms.InputTag("myProducer", "slimmedMuons"),
  matched = cms.InputTag("myProducer", "prunedGenParticles"),
  distMin = cms.double(0.15),
  matchPDGId = cms.vint32(13)
)

# process.allMuonsGenParticlesMatch = cms.EDFilter("MCTruthDeltaRMatcher",
#      src = cms.InputTag("allMuons"),
#      distMin = cms.double(0.15),
#      matchPDGId = cms.vint32(13),
#      matched = cms.InputTag("genParticleCandidates")
#  )


# process.Tracer = cms.Service("Tracer")

# Output module to store the produced collection (if needed)
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('final.root'),
    # outp
)

# process.p = cms.Path(process.demo * process.selectedMuonsGenParticlesMatch)
process.p = cms.Path(process.selectedMuonsGenParticlesMatch)
process.e = cms.EndPath(process.out)