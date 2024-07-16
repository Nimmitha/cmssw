import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                        #  filterEfficiency = cms.untracked.double(0.026),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                        #  crossSection = cms.untracked.double(174000000.0),
                         comEnergy = cms.double(13600.0),
                         maxEventsToPrint = cms.untracked.int32(1),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring(
            'Higgs:useBSM = on', 
            'HiggsBSM:gg2H2 = on', 
            'HiggsH2:coup2d = 10.0', 
            'HiggsH2:coup2u = 10.0', 
            'HiggsH2:coup2Z = 0.0', 
            'HiggsH2:coup2W = 0.0', 
            'HiggsA3:coup2H2Z = 0.0', 
            'HiggsH2:coup2A3A3 = 0.0', 
            'HiggsH2:coup2H1H1 = 0.0', 
            '443:onMode = off', 
            '443:onIfMatch 13 13', 
            '333:onMode = off', 
            '333:onIfMatch 13 13', 
            '553:onMode = off', 
            '553:onIfMatch 13 13',
            '23:onMode = off',
            '23:onIfMatch 13 13',   
            '35:mMin = 0', 
            '35:mMax = 200', 
            '35:m0   = 125.', 
            '35:mWidth = 0.0001', 
            '35:addChannel 1 1.00 100 23 553', 
            '35:spinType = 1', 
            '35:onMode = off', 
            '35:onIfMatch 23 553'
            ),
        parameterSets = cms.vstring(
                                    'pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'processParameters',
                                    )
        )
)

# Next two muon filter are derived from muon reconstruction
# mumugenfilter = cms.EDFilter("MCParticlePairFilter",
#     Status = cms.untracked.vint32(1, 1),
#     MinPt = cms.untracked.vdouble(0.5, 0.5),
#     MinP = cms.untracked.vdouble(0., 0.),
#     MaxEta = cms.untracked.vdouble(2.5, 2.5),
#     MinEta = cms.untracked.vdouble(-2.5, -2.5),
#     MinInvMass = cms.untracked.double(2.0),
#     MaxInvMass = cms.untracked.double(4.0),
#     ParticleCharge = cms.untracked.int32(-1),
#     ParticleID1 = cms.untracked.vint32(13),
#     ParticleID2 = cms.untracked.vint32(13)
# )

etafilter = cms.EDFilter("PythiaFilter",
    MaxEta = cms.untracked.double(9999.0),
    MinEta = cms.untracked.double(-9999.0),
    ParticleID = cms.untracked.int32(35)
)

# oniafilter = cms.EDFilter("PythiaFilter",
#     Status = cms.untracked.int32(2),
#     MaxEta = cms.untracked.double(1000.0),
#     MinEta = cms.untracked.double(-1000.0),
#     MinPt = cms.untracked.double(8.0),
#     ParticleID = cms.untracked.int32(443)
# )

# ProductionFilterSequence = cms.Sequence(generator*oniafilter*mumugenfilter)
ProductionFilterSequence = cms.Sequence(generator*etafilter)
