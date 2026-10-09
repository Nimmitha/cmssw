from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Muon1_mm_2024F1_tnp_v9'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniAODmuonsRootupler.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['Muon1_mm_2024F1_tnp_v9.root']
config.JobType.maxJobRuntimeMin = 240

config.Data.inputDBS = 'global'
config.Data.inputDataset = '/Muon1/Run2024F-PromptReco-v1/MINIAOD'
config.Data.lumiMask = 'Cert_Collisions2024_378981_386951_Golden.json'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 25
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/nkarunar/tag/'
config.Data.publication = False

config.Site.storageSite = 'T3_US_FNALLPC'
