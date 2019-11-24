from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.General.requestName = 'Ds_MC_Run2018A_WithoutCuts'
config.General.workArea = 'crab_projects' #folder - diferents projects may be put in the same work area
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'dstard0_cfg_MC.py' #config file used

config.Data.inputDataset = '/MinBias_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/MINIAODSIM'
#config.Data.inputDataset = '/DStarToD0Pi_D0K3Pi_DStarFilter_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/MINIAODSIM'
#config.Data.inputDataset = '/DStarToD0Pi_D0KPi_DStarFilter_TuneCP5_13TeV-pythia8-evtgen/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
#config.Data.inputDataset = '/DStarToD0Pi_D0KPi_DStarFilter_TuneCP5_13TeV-pythia8-evtgen/RunIIFall17MiniAODv2-NoPU_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10 #if it increases, the number of jobs decreases

#config.JobType.outputFiles = ['D0DstarMC.root']
config.Data.outLFNDirBase = '/store/user/ragomesd/crab' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag= 'Ds_MC_13TeV-pythia8'

#config.Site.storageSite     = 'T3_US_FNALLPC'
config.Site.storageSite = 'T2_CH_CERNBOX' #T2_CH_CERNBOX # T2_US_Nebraska
#config.Site.storageSite = 'T2_PL_Swierk'
#config.Site.blacklist = ['T2_US_Nebraska'] #Avoid that the job runs in the T2 that failed
