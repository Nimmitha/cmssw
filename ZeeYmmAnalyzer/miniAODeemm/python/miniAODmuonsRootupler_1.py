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


# process.load('RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16_ID_ISO_cff')
# process.load('RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer17UL_ID_ISO_cff')
# process.load('RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer18UL_ID_ISO_cff')


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v37', '') # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis
# process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v32', '')#2016,2017,2018
# process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v15_L1v1', '')#MC

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True,
                       runVID=True,
                       era='2018-UL')

                     # eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer18UL_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                     #  eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Autumn18_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],

process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))#,SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(     
        # 'file:../../../prep/CMSSW_10_6_20/src/MiniAOD.root'
        # '/store/data/Run2018B/EGamma/MINIAOD/UL2018_MiniAODv2-v1/260000/00FAE4A8-E0E7-A941-AB7E-0DD0383C9CB5.root'
        '/store/data/Run2018B/EGamma/MINIAOD/UL2018_MiniAODv2-v1/260000/041F6B17-C831-5249-BD07-90756D24D544.root'
        # '/store/data/Run2018B/EGamma/MINIAOD/UL2018_MiniAODv2-v1/260000/041F6B17-C831-5249-BD07-90756D24D544.root'
#'file:simWork/retryFromTheBegining/RECOandMINI/MiniAOD.root'
#'file:simWork/realFromTheBegining/RECOandMini/MiniAOD.root'
#'file:ZeeY2S/a1/mini/MiniAOD.root'
# '/store/data/Run2017B/DoubleEG/MINIAOD/09Aug2019_UL2017-v1/00000/011FE2C9-C2E7-034B-B571-C30A5385F83F.root'
# 'file:ZeeY3S/a9/MiniAOD.root'
#'file:simWork/realFromTheBegining/ZeeJuu/s9/MiniAOD.root'
#'/DPS_ToYZ_YToMuMu_ZToMuMu_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1/MINIAODSIM'

# '/store/mc/RunIIFall17MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/60A88420-4E42-E811-B1E6-001E675813C4.root'

#'/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToZZTo4L_M125_CP5TuneDown_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/5EDD93D8-D09C-DD4D-B86E-4235A88648C0.root'# UL ZZ
#'file:ZeeYuu_2018_MC_real_MINIAOD.root'
#'/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/017E4F7F-26F4-2D4D-87C5-6E58954FDEC8.root'
#        'file:/uscms/home/jdharris/nobackup/YOURWORKINGAREA/testing/MCprocessing/again/CMSSW_10_2_18/src/forReal/ZeeJuu/s1/ZeeJuu_2018_MC_real_MINIAOD.root'
#'/store/data/Run2017B/SingleElectron/MINIAOD/09Aug2019_UL2017-v1/130000/005B3780-5B9B-3B4F-ABA9-DCC474BE15BE.root'
 )
)

#Prueba de fuego.
#process.load("myAnalyzers.JPsiKsPAT.miniAODeemmRootupler_cfi")
#process.rootuple.dimuons = cms.InputTag('miniaodPATMuonsWithTrigger') 

process.rootuple = cms.EDAnalyzer('miniAODeemm',
                          dimuons = cms.InputTag("slimmedMuons"),
                          dielectron = cms.InputTag("slimmedElectrons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bits = cms.InputTag("TriggerResults::HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          pruned = cms.InputTag("prunedGenParticles"),
                          MuonTrigger = cms.string("HLT_IsoMu24_v"),
                          ElectronTrigger = cms.string("HLT_Ele27_WPTight_Gsf_v"),
                          DataType = cms.string("2018B_UL"),  # Title of the output ROOT file
                          isMC = cms.bool(False),
                          )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('JesseTestMC_out.root'),
)

process.p = cms.Path(process.egammaPostRecoSeq+process.rootuple)


