import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")

#process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = "101X_dataRun2_Prompt_v9"

#number of events (-1 -> all)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
IgnoreCompletely = cms.untracked.vstring('ProductNotFound'),
SkipEvent = cms.untracked.vstring('Error: uninitialized ProxyBase used') )

#sourcefile = 'list259431.txt'

#USE THIS BLOCK FOR CONVENTIONAL FILE OPENING!===========
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
'file:FF14A12B-98F4-0F4E-AAF5-31FBA6F037EE_MC.root'
	   )
)
#===================================================================

process.analysis = cms.EDAnalyzer('DstarD0TTree',
	# Analysis
	doMC = cms.bool(True),
	doRec = cms.bool(True),
	bits = cms.InputTag("TriggerResults","","HLT"),
	prescales = cms.InputTag("patTrigger"),#linked to pat::PackedTriggerPrescales
	#PathName = cms.untracked.string("HLT_Mu8_IP3_part0_v"), #ParkingBPH1
	PathName = cms.untracked.string("HLT_Mu9_IP6_part"),  #ParkingBPH1
	
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
	fileName = cms.string('D0DstarMC.root') #output files
)

process.p = cms.Path(process.analysis)
