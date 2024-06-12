from CRABClient.UserUtilities import config
config = config()

# user specific generic parameters
config.General.requestName = 'TTree_13TeV_Muon1_mmmm_2023D_v2'   # Used as the task/Project directory name
config.General.transferOutputs = True                               # Transfer output files to the storage site
config.General.transferLogs = False

# job type and related configurables
config.JobType.pluginName = 'Analysis'                              # Specify: analysis or MC generation
config.JobType.psetName = 'miniAODmuonsRootupler_2023_M1_D_v2.py'           # parameter-set config file
config.JobType.allowUndistributedCMSSW = True                       # Allow CMSSW release possibly not available at sites
config.JobType.outputFiles = ['Muon1_Run2023D_v2_v1_Data.root']   # List of output files that needs to be collected
# config.JobType.maxMemoryMB = 1000

# data to be analyzed
config.Data.inputDBS = 'global'
config.Data.inputDataset = '/Muon1/Run2023D-22Sep2023_v2-v1/MINIAOD'               # Name of the dataset
#config.Data.outputPrimaryDataset = 'DoubleJPsiToMuMu_RAWSIM_SPS_LO_may2016_largetest_FNAL' # Used when running private input files or MC generation
config.Data.lumiMask = 'Cert_Collisions2023_366442_370790_Golden.json' # Lumi-section filter
config.Data.splitting = 'LumiBased'                                                         # Split the task based on 
config.Data.unitsPerJob = 50                                                                # Number of splitted units per job
#NJOBS = 500  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = -1                                                                 # Number of untis to analyze 
config.Data.outLFNDirBase = '/store/user/nkarunar/'
config.Data.publication = False                                                             # Whether to publish the EDM output files in DBS

# Grid site parameters
config.Site.storageSite = 'T3_US_FNALLPC'         # Place to copy the output files
