import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")

#process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = "101X_dataRun2_Prompt_v9"

#number of events (-1 -> all)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
IgnoreCompletely = cms.untracked.vstring('ProductNotFound'),
SkipEvent = cms.untracked.vstring('Error: uninitialized ProxyBase used') )

#sourcefile = 'list259431.txt'

#USE THIS BLOCK FOR CONVENTIONAL FILE OPENING!===========
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
#'root://cmsxrootd.fnal.gov//store/data/Run2018A/ParkingBPH1/MINIAOD/14May2018-v1/30000/129EAF02-B85A-E811-AC4E-0CC47A1DF818.root'
#'file:129EAF02-B85A-E811-AC4E-0CC47A1DF818.root' #Data mniAOD
'file:EA334703-C759-E811-B440-008CFAE45030_BParking1.root'
#'file:F800AD57-C458-E811-AC54-0025904C67B4_BParking2.root'
#'file:E0B15542-285A-E811-A365-001E67E71D03_BParking3.root'
#'file:FCC64175-325A-E811-A1EF-FA163E17FD8B_BParking4.root'
#'file:FA6C34C5-8E58-E811-8048-008CFAFBE5CE_Bparking5.root'
#'file:FCEDA0BD-3F5B-E811-98D6-68B59972BF74_BParking6.root'
	   )
)
#===================================================================

process.analysis = cms.EDAnalyzer('DstarD0TTree',
	# Analysis
	doMC = cms.bool(False), #True to run on MC
	doRec = cms.bool(True),
	bits = cms.InputTag("TriggerResults","","HLT"),
	prescales = cms.InputTag("patTrigger"),#linked to pat::PackedTriggerPrescales
	objects = cms.InputTag("selectedPatTrigger"),
	#PathName = cms.untracked.string("HLT_Mu8p5_IP3p5_part0"), #ParkingBPH1
	#PathName = cms.untracked.string("HLT_Mu12_IP6_part0"),  #ParkingBPH1
	PathName = cms.untracked.string("HLT_Mu9_IP6_part2"),  #ParkingBPH1
	
	tracks = cms.InputTag('packedPFCandidates'),#linked to vector<pat::PackedCandidate
	recVtxs = cms.InputTag('offlineSlimmedPrimaryVertices'), #linked to vector<reco::Vertex> 
	gens = cms.InputTag("prunedGenParticles"), #linked to reco::GenParticleCollection
	gensD0 = cms.InputTag("prunedGenParticles"), #linked to reco::GenParticleCollection
	# Options
	comEnergy = cms.double(13000.),
	TTBIt = cms.int32(34),
	debug = cms.untracked.bool(False),
	SaveROOTTree = cms.untracked.bool(True)
)

process.TFileService = cms.Service("TFileService",
   fileName = cms.string('D0DstarDataBparking.root') #output files
)

process.p = cms.Path(process.analysis)
