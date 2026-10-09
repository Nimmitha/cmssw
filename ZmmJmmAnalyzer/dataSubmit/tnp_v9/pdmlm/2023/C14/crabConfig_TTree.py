from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'PDMLM1_2023C4_tnp_v9'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniAODmuonsRootupler.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['PDMLM1_2023C4_tnp_v9.root']
config.JobType.maxJobRuntimeMin = 240

config.Data.inputDBS = 'global'
config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2023C-22Sep2023_v4-v1/MINIAOD'
config.Data.lumiMask = 'Cert_Collisions2023_366442_370790_Golden.json'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 200
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/nkarunar/tnp_v9/'
config.Data.publication = False

config.Site.storageSite = 'T3_US_FNALLPC'
