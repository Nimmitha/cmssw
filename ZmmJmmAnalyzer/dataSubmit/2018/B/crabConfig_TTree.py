from CRABClient.UserUtilities import config
config = config()

# user specific generic parameters
config.General.requestName = 'TTree_13TeV_mmmm_UL_2018B'   # Used as the task/Project directory name
config.General.transferOutputs = True                               # Transfer output files to the storage site
config.General.transferLogs = False

# job type and related configurables
config.JobType.pluginName = 'Analysis'                              # Specify: analysis or MC generation
config.JobType.psetName = 'miniAODmmmmRootupler_2018_B.py'           # parameter-set config file
config.JobType.allowUndistributedCMSSW = True                       # Allow CMSSW release possibly not available at sites
config.JobType.outputFiles = ['SingleMuon_Run2018B_UL_v2_v2_Data.root']   # List of output files that needs to be collected
# config.JobType.maxMemoryMB = 1000

# data to be analyzed
config.Data.inputDBS = 'global'
config.Data.inputDataset = '/SingleMuon/Run2018B-UL2018_MiniAODv2-v2/MINIAOD'               # Name of the dataset
#config.Data.outputPrimaryDataset = 'DoubleJPsiToMuMu_RAWSIM_SPS_LO_may2016_largetest_FNAL' # Used when running private input files or MC generation
config.Data.lumiMask = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_MuonPhys.txt' # Lumi-section filter
config.Data.splitting = 'LumiBased'                                                         # Split the task based on 
config.Data.unitsPerJob = 50                                                                # Number of splitted units per job
#NJOBS = 500  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = -1                                                                 # Number of untis to analyze 
config.Data.outLFNDirBase = '/store/user/nkarunar/'
config.Data.publication = False                                                             # Whether to publish the EDM output files in DBS

# Grid site parameters
config.Site.storageSite = 'T3_US_FNALLPC'         # Place to copy the output files
