from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'PDMLM6_mm_2023D1_tnp_v9'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniAODmuonsRootupler.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['PDMLM6_mm_2023D1_tnp_v9.root']
config.JobType.maxJobRuntimeMin = 240

config.Data.inputDBS = 'global'
config.Data.inputDataset = '/ParkingDoubleMuonLowMass6/Run2023D-22Sep2023_v1-v1/MINIAOD'
config.Data.lumiMask = 'Cert_Collisions2023_366442_370790_Golden.json'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 200
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/nkarunar/tag/'
config.Data.publication = False

config.Site.storageSite = 'T3_US_FNALLPC'
