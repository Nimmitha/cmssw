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
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v37', '') # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       #runVID=True, #if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       runEnergyCorrections=True,
                       runVID=True,
                       #eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                       era='2016postVFP-UL')                       
#                       )  #era is new to select between 2016 / 2017,  it defaults to 2017

process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True),
    # skipEvent = cms.untracked.vstring('ProductNotFound')
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:test.root')
)

process.rootuple = cms.EDAnalyzer('miniAODeemm',
                          dimuons = cms.InputTag("slimmedMuons"),
                          dielectron = cms.InputTag("slimmedElectrons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bits = cms.InputTag("TriggerResults::HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          pruned = cms.InputTag("prunedGenParticles"),
                          ElectronTrigger = cms.string("HLT_Ele27_WPTight_Gsf_v"), # 2016
                          DataType = cms.string("2016_UL"),  # Title of the output ROOT file
                          isMC = cms.bool(False),
                          )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('SingleElectron_Run2016Fpost_UL_v2_v2_Data.root'),
)

process.p = cms.Path(process.egammaPostRecoSeq+process.rootuple)

