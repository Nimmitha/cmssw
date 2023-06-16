import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/cms/Tutorials/workbook_twiki2021/MinBias_pythia8_14TeV_100events.root'
                )
                            )

process.demo = cms.EDAnalyzer('MyAnalyzer',
   tracks    = cms.untracked.InputTag('generalTracks'),
   trackPtMin = cms.double(0.3)
                              )

process.p = cms.Path(process.demo)