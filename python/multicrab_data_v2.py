
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
#from CRABClient.UserUtilities #import config, getUsernameFromSiteDB
from httplib import HTTPException

from WMCore.Configuration import Configuration


# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
config = Configuration()
config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'ConfFile_data_cfg.py'
#config.JobType.inputFiles = ['MyDataPileupHistogram_2017.root','PileupMC_2017.root']
#config.JobType.outputFiles = ['out.root']
config.JobType.allowUndistributedCMSSW = True

   
config.section_("Data")
config.Data.useParent = True
config.Data.ignoreLocality = True
config.Data.splitting = 'LumiBased' #'FileBased' 'Automatic' 'FileBased' 'LumiBased' 'EventAwareLumiBased'
config.Data.unitsPerJob = 50
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'  #'Automatic'
#config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/ragomesd/crab' #%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag =  'data'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
config.Data.runRange = '314472-325175'

config.section_("Site")
config.Site.whitelist = ['T2_CH_*','T2_DE_*','T2_IT_*','T2_US_*']
config.Site.storageSite = "T2_CH_CERNBOX" #T2_CH_CERNBOX  #T2_BR_UERJ



def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

if __name__ == '__main__':


#example
#config.General.requestName = 'FSQJet1Bprompt'
#config.Data.inputDataset = '/FSQJet1/Run2018B-PromptReco-v2/MINIAOD'
#submit(config)

    config.JobType.psetName = 'dstard0_cfg2.py'
    config.General.requestName = 'Bparking2_Run2018A'
    config.Data.inputDataset = '/ParkingBPH2/Run2018A-14May2018-v1/MINIAOD'
    submit(config)

    config.JobType.psetName = 'dstard0_cfg3.py'
    config.General.requestName = 'Bparking3_Run2018A'
    config.Data.inputDataset = '/ParkingBPH3/Run2018A-14May2018-v1/MINIAOD'
    submit(config)

    config.JobType.psetName = 'dstard0_cfg4.py'
    config.General.requestName = 'Bparking4_Run2018A'
    config.Data.inputDataset = '/ParkingBPH4/Run2018A-14May2018-v1/MINIAOD'
    submit(config)

    config.JobType.psetName = 'dstard0_cfg5.py'
    config.General.requestName = 'Bparking5_Run2018A'
    config.Data.inputDataset = '/ParkingBPH5/Run2018A-14May2018-v1/MINIAOD'
    submit(config)

    config.JobType.psetName = 'dstard0_cfg6.py'
    config.General.requestName = 'Bparking6_Run2018A'
    config.Data.inputDataset = '/ParkingBPH6/Run2018A-14May2018-v1/MINIAOD'
    submit(config)



    

