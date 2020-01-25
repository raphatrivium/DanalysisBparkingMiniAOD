#from CRABClient.UserUtiilities import config, getUsernameFromSiteDB

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).

from WMCore.Configuration import Configuration

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
#config.JobType.maxJobRuntimeMin = 60
config.JobType.psetName = 'null' #'test_cfg-el.py'
#config.JobType.maxMemoryMB = 760
config.JobType.inputFiles = ['MyDataPileupHistogram.root','PileupMC.root']
#config.JobType.outputFiles = ['output_HepMC.root', 'output_track_xy.root']

config.section_('Data')
#config.Data.useParent = True
#config.Data.ignoreLocality = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased' #'FileBased' 'Automatic' 'FileBased' 'LumiBased' 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
config.Data.allowNonValidInputDataset = True
#config.Data.totalUnits = 100
#config.Data.outLFNDirBase = 'gsiftp://eosuserftp.cern.ch/eos/user/m/malvesga' #%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'Ds_MC_Run2018A_Official_Tracks'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERNBOX' #'T2_CH_CERNBOX' 'T2_IT_Bari'
#config.Site.whitelist = ['T2_CH_CERN','T2_DE_*','T2_IT_Bari','T2_US_*']


def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

if __name__ == '__main__':

    config.General.workArea = 'crab_projects'
    config.General.requestName = 'Ds_MC_Run2018A_Official_Tracks'
    config.JobType.psetName = 'dstard0_cfg_MC.py'
    config.Data.inputDataset = '/DStarToD0Pi_D0KPi_NoMuFilter_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
    config.Data.outLFNDirBase = '/store/user/ragomesd/crab'
    submit(config)

    config.General.workArea = 'crab_projects'
    config.General.requestName = 'Ds_MC_Run2018A_MinBias_Tracks'
    config.JobType.psetName = 'dstard0_cfg_MC.py'
    config.Data.inputDataset = '/MinBias_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/MINIAODSIM'
    config.Data.outLFNDirBase = '/store/user/ragomesd/crab'
    submit(config)	
       
