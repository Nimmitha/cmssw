import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            # 'file:/afs/cern.ch/cms/Tutorials/workbook_twiki2021/MinBias_pythia8_14TeV_100events.root'
            # 'file:/eos/home-n/nkarunar/2023Analysis/datasets/2022C/slimMiniAOD_data_MuEle.root'
            'file:/eos/home-n/nkarunar/2023Analysis/datasets/2018_MC/Zuu_2018_MINIAOD.root'
                )
                            )

process.demo = cms.EDAnalyzer('MyAnalyzer',
   tracks    = cms.untracked.InputTag('generalTracks'),
   trackPtMin = cms.double(0.3)
                              )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('PatHistos.root'),
)

# process.Tracer = cms.Service("Tracer")

process.p = cms.Path(process.demo)