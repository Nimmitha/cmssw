import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
#process.load('RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16_ID_ISO_cff')
# process.load("Configuration.Geometry.GeometryIdeal_cff")

# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
# from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc') # Tried: 106X_dataRun2_v32

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v15') # 2022 check?
# process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1') # 2023 check?
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_mcRun3_2022_realistic_v12') # used to generate 2022 MC
# process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_v2')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data')



process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  allowUnscheduled = cms.untracked.bool(True),
  SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(       
#  '/store/data/Run2022B/SingleMuon/MINIAOD/22Sep2023-v1/50000/2ad62630-b826-4686-a441-2099de25476a.root'
#  '/store/data/Run2023D/Muon0/MINIAOD/22Sep2023_v1-v1/2530000/0fa55ced-a7cc-4a8a-a0e5-0381c3ac8e37.root'
# 'file:/uscms/home/wkarunar/nobackup/datasets/mc/zmmymm_run3/mmmm_v1/MiniAOD/MiniAOD_10.root'
'file:/uscms/home/wkarunar/nobackup/datasets/mc/run3_zmmJpsimm/ZmmJpsimm/MiniAOD/MiniAOD_10.root'
#'file:../../../../datasets/ZmmYee/Y1S/MiniAOD/MiniAOD_1.root' 
 )
)

process.rootuple = cms.EDAnalyzer('miniAODmmmm',
                          muons = cms.InputTag("slimmedMuons"),
                          # dielectron = cms.InputTag("slimmedElectrons"),
                          # Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bits = cms.InputTag("TriggerResults::HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          pruned = cms.InputTag("prunedGenParticles"),
                          MuonTrigger = cms.string("HLT_IsoMu24_v"),
                          isMC = cms.bool(True),
                          )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('inputFiles/mc_zmmjpsimm_v1/mc_zmmjpsimm_v1_10.root'),
)

#process.p = cms.Path(process.egammaPostRecoSeq+process.rootuple)
process.p = cms.Path(process.rootuple)

