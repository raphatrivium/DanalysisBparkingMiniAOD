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
#config.JobType.inputFiles = ['alignment.xml','direct_simu_reco_cff.py','customisation_cff.py','PUHistos_data.root','PUHistos_mc.root']
#config.JobType.outputFiles = ['output_HepMC.root', 'output_track_xy.root']

config.section_('Data')
#config.Data.useParent = True
#config.Data.ignoreLocality = True
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
config.Data.runRange = '314472-325175'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased' #'FileBased' 'Automatic' 'FileBased' 'LumiBased' 'EventAwareLumiBased'
config.Data.unitsPerJob = 5
#config.Data.totalUnits = 100
#config.Data.outLFNDirBase = 'gsiftp://eosuserftp.cern.ch/eos/user/m/malvesga' #%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'ParkingBPH_Run2018A_MINIAOD_triggers'

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
    config.General.requestName = 'Bparking1_Run2018A_triggers'
    config.JobType.psetName = 'dstard0_cfg1_triggers.py'
    config.Data.inputDataset = '/ParkingBPH1/Run2018A-14May2018-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/user/ragomesd/crab'
    submit(config)

    config.General.workArea = 'crab_projects'
    config.General.requestName = 'Bparking2_Run2018A_triggers'
    config.JobType.psetName = 'dstard0_cfg1_triggers.py'
    config.Data.inputDataset = '/ParkingBPH2/Run2018A-14May2018-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/user/ragomesd/crab'
    submit(config)

    config.General.workArea = 'crab_projects'
    config.General.requestName = 'Bparking3_Run2018A_triggers'
    config.JobType.psetName = 'dstard0_cfg1_triggers.py'
    config.Data.inputDataset = '/ParkingBPH3/Run2018A-14May2018-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/user/ragomesd/crab'
    submit(config)

    config.General.workArea = 'crab_projects'
    config.General.requestName = 'Bparking4_Run2018A_triggers'
    config.JobType.psetName = 'dstard0_cfg1_triggers.py'
    config.Data.inputDataset = '/ParkingBPH4/Run2018A-14May2018-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/user/ragomesd/crab'
    submit(config)

    config.General.workArea = 'crab_projects'
    config.General.requestName = 'Bparking5_Run2018A_triggers'
    config.JobType.psetName = 'dstard0_cfg1_triggers.py'
    config.Data.inputDataset = '/ParkingBPH5/Run2018A-14May2018-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/user/ragomesd/crab'
    submit(config)

    config.General.workArea = 'crab_projects'
    config.General.requestName = 'Bparking6_Run2018A_triggers'
    config.JobType.psetName = 'dstard0_cfg1_triggers.py'
    config.Data.inputDataset = '/ParkingBPH6/Run2018A-14May2018-v1/MINIAOD'
    config.Data.outLFNDirBase = '/store/user/ragomesd/crab'
    submit(config)

