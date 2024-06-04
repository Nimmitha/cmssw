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

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v32', '')#2016,2017,2018
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v15_L1v1', '')#MC

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True,
                       runVID=True,
                       era='2018-UL')

process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))#,SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
        'file:/uscms/home/wkarunar/nobackup/datasets/ZmmYee/Y3S/MiniAOD/MiniAOD_5.root'
        # '/store/data/Run2018B/SingleMuon/MINIAOD/12Nov2019_UL2018-v3/00000/006332C4-470B-A44B-A281-8B0B10A7D591.root'
        # 'file:/uscms/home/wkarunar/nobackup/Analysis/ZeeYuu_Jesse/CMSSW_10_6_20/src/InputFiles/simWork/realFromTheBegining/RECOandMini/MiniAOD.root'
 )
)

process.rootuple = cms.EDAnalyzer('miniAODmuons',
                          dimuons = cms.InputTag("slimmedMuons"),
                          dielectron = cms.InputTag("slimmedElectrons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bits = cms.InputTag("TriggerResults::HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          pruned = cms.InputTag("prunedGenParticles"),
                          MuonTrigger = cms.string("HLT_IsoMu24_v"),
                          ElectronTrigger = cms.string("HLT_Ele27_WPTight_Gsf_v"),
                          # DataType = cms.string("2018B MC"),  # Title of the output ROOT file
                          # isMC = cms.bool(True),
                          DataType = cms.string("2018 MC Y3S"),  # Title of the output ROOT file
                          isMC = cms.bool(True),
                          )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('output/ZmmMCOut_Y3S_5.root'),
  # fileName = cms.string('output/.root'),
)

process.p = cms.Path(process.egammaPostRecoSeq+process.rootuple)


