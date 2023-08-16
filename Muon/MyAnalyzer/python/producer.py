import FWCore.ParameterSet.Config as cms

process = cms.Process("MyProducerDemo")

# Replace 'source' with your desired input source (e.g., input file)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/home-n/nkarunar/2023Analysis/datasets/2018_MC/Zuu_2018_MINIAOD.root')
)

# Define your producer module and provide input parameters (e.g., prunedGenParticles)
process.myProducer = cms.EDProducer('MyProducer',
    prunedGenParticles = cms.InputTag('prunedGenParticles'),  # Adjust the input tag if necessary
    slimmedMuons = cms.InputTag('slimmedMuons')
)

# Output module to store the produced collection (if needed)
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('output.root'),
    # outp
)

# Adjust the number of events to process as needed
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Adjust the frequency of printouts as needed
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout'),
    cout = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
)

# Endpath to define the sequence of modules to run
process.myPath = cms.EndPath(process.myProducer * process.out)

# Uncomment this line if you want to write the output to a file
# process.e = cms.EndPath(process.out)

# Uncomment this line if you want to run the output to a file without producing any output
# process.p = cms.EndPath()