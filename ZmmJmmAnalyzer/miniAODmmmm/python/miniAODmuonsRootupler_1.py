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
# process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v32', '') # data
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v37', '') # data UL
# process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v15_L1v1', '') # MC


from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       #runVID=True, #if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       runEnergyCorrections=True,
                       runVID=True,
                       #eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                       era='2018-UL')                 
                      #  era='2017-UL')                 
                      #  era='2016postVFP-UL')                 
#                       )  #era is new to select between 2016 / 2017,  it defaults to 2017


process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  allowUnscheduled = cms.untracked.bool(True),
  # SkipEvent = cms.untracked.vstring('ProductNotFound')
  )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:miniAODs/000464E1-1144-1641-BE88-4600BD58923C.root'
        # 'file:/eos/uscms/store/user/wkarunar/MC2018_UL/hzmmpsimm/MC2018_UL_hzmmpsimm_MINIAOD/241227_060243/0000/MINIAOD_2018_UL_13TeV_hzmmpsimm_mc_v1_10.root'
        # 'file:/eos/uscms/store/user/wkarunar/MC2018_UL/hzmmx2mm/MC2018_UL_hzmmx2mm_MINIAOD/241231_193736/0000/MINIAOD_2018_UL_13TeV_hzmmx2mm_mc_v1_3.root'
        # 'file:/uscms/home/wkarunar/nobackup/datasets/mc/miniAOD/run2/zmmjmm_2018/zmmjmm_ss/MiniAOD/MiniAOD_10.root'
 )
)

process.rootuple = cms.EDAnalyzer('miniAODmmmm',
                          dimuons = cms.InputTag("slimmedMuons"),
                          dielectron = cms.InputTag("slimmedElectrons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bits = cms.InputTag("TriggerResults::HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          pruned = cms.InputTag("prunedGenParticles"),
                          MuonTrigger = cms.string("HLT_IsoMu24_v"), # 2016 and 18
                          # MuonTrigger = cms.string("HLT_IsoMu27_v"), # 2017
                          isMC = cms.bool(True),
                          )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('preselection/jmm_output_23C.root'),
)

process.p = cms.Path(process.egammaPostRecoSeq+process.rootuple)

