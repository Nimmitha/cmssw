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
# process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v15') # 2022 check?
# process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1') # 2023 check?
# process.GlobalTag = GlobalTag(process.GlobalTag, '124X_mcRun3_2022_realistic_v12') # used to generate 2022 MC
# process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_v2')
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_dataRun3_Prompt_v3')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data')



process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  allowUnscheduled = cms.untracked.bool(True),
  # SkipEvent = cms.untracked.vstring('ProductNotFound')
  )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(       
        'file:miniAODs/0f4cf3e8-65e8-470e-a587-7a676ffcfc89_2024H0.root'
# 'file:/uscms/home/wkarunar/nobackup/datasets/mc/zmmymm_run3/mmmm_v1/MiniAOD/MiniAOD_10.root'
# '/store/data/Run2022D/Muon/MINIAOD/22Sep2023-v1/2520000/00242360-6d04-4eb0-b75c-0d743850f2fc.root'
#'file:../../../../datasets/ZmmYee/Y1S/MiniAOD/MiniAOD_1.root' 
 )
)

process.rootuple = cms.EDAnalyzer('miniAODmmmm',
                          muons = cms.InputTag("slimmedMuons"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bits = cms.InputTag("TriggerResults::HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          prescales = cms.InputTag("patTrigger"),
                          pruned = cms.InputTag("prunedGenParticles"),
                          MuonTrigger = cms.string("HLT_Mu0_L1DoubleMu_v"),
                          isMC = cms.bool(False),
                          )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('preselection/testout_2024H0.root'),
)

process.p = cms.Path(process.rootuple)

