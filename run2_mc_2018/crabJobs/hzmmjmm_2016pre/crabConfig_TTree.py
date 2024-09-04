from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = '2016pre_13TeV_hzmmjpsimm_ss_mc'                   # Used as the task/Project directory name
config.General.transferOutputs = True                               # Transfer output files to the storage site
config.General.transferLogs = False
# config.General.workArea = 'crabsubmit'

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'                             # Specify: analysis or MC generation
config.JobType.psetName = 'hzj-2016-pre.py'              # parameter-set config file
# config.JobType.generator = 'lhe'
# config.JobType.eventsPerLumi = 100
config.JobType.outputFiles = ['SIM_2016pre_13TeV_hzmmjpsimm_mc_v1.root']   # List of output files that needs to be collected
config.JobType.allowUndistributedCMSSW = True                       # Allow CMSSW release possibly not available at sites

config.section_("Data")
config.Data.outputPrimaryDataset = 'MinBias'    
config.Data.splitting = 'EventBased'                                # Split the task based on 
config.Data.unitsPerJob = 1000                                        # Number of splitted units per job
NJOBS = 10  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS            # Number of untis to analyze 
config.Data.outLFNDirBase = '/store/user/nkarunar/'
config.Data.publication = False                                      # Whether to publish the EDM output files in DBS
config.Data.outputDatasetTag = 'MinBias_TuneCP5_13TeV-pythia8_2016pre_hzmmjpsimm_ss_mc'

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'                           # Place to copy the output files