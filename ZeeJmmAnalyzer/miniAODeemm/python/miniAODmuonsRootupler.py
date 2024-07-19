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
# process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v15', '') # 2022 B-E
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1') # 2023

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  allowUnscheduled = cms.untracked.bool(True),
  SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(       
#  '/store/data/Run2022B/EGamma/MINIAOD/22Sep2023-v2/70000/01921504-9d3d-4993-a46b-e45dc979a2a4.root'
#  '/store/data/Run2022B/EGamma/MINIAOD/27Jun2023-v2/2520000/2cbda963-c777-4a15-9f7a-56b1be95099c.root'
 '/store/data/Run2023C/EGamma0/MINIAOD/22Sep2023_v1-v1/2530000/0180e051-34b5-47a0-a851-665aa47846fd.root'
 )
)

# from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# dataFormat = DataFormat.MiniAOD
# switchOnVIDElectronIdProducer(process, dataFormat)
# process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

# my_id_modules = [
#     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_RunIIIWinter22_V1_cff'
# ]


process.rootuple = cms.EDAnalyzer('miniAODeemm',
                          dimuons = cms.InputTag("slimmedMuons"),
                          dielectron = cms.InputTag("slimmedElectrons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bits = cms.InputTag("TriggerResults::HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          pruned = cms.InputTag("prunedGenParticles"),
                          isMC = cms.bool(False),
                          )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('Testout_2022B_27Jun_zeeymm.root'),
)

# process.p = cms.Path(process.electronMVAValueMapProducer + process.rootuple)
process.p = cms.Path(process.rootuple)

