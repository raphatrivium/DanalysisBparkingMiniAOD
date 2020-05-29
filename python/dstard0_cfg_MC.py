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
'file:FE4C1BFF-3224-134B-94EB-49E603290280_OfficialMC.root'
#'file:F79E9F51-558B-3445-B2B6-7160DF5FE99E_MinBiasMC.root'	
	   )
)
#===================================================================
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
	vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
	minimumNDOF = cms.uint32(4) ,
	maxAbsZ = cms.double(24), 
	maxd0 = cms.double(2) 
)

#process.noscraping = cms.EDFilter("FilterOutScraping",
#applyfilter = cms.untracked.bool(True),
#debugOn = cms.untracked.bool(False),
#numtrack = cms.untracked.uint32(10),
#thresh = cms.untracked.double(0.25)
#)

process.analysis = cms.EDAnalyzer('DstarD0TTree',
	# Analysis
	doMC = cms.bool(True),  #True to run on MC
	doRec = cms.bool(True),
	bits = cms.InputTag("TriggerResults","","HLT"),
	prescales = cms.InputTag("patTrigger"),#linked to pat::PackedTriggerPrescales
	objects = cms.InputTag("selectedPatTrigger"),
	PathName = cms.untracked.string("HLT_Mu9_IP6"),  #triggerName	HLT_Mu8_IP3 / HLT_Mu8p5_IP3p5 / HLT_Mu9_IP6 / HLT_Mu10p5_IP3p5
	ReferencePathName = cms.untracked.string("HLT_Mu8_IP3"),
	tracks = cms.InputTag('packedPFCandidates'),#linked to vector<pat::PackedCandidate
	recVtxs = cms.InputTag('offlineSlimmedPrimaryVertices'), #linked to vector<reco::Vertex> 
	gens = cms.InputTag("prunedGenParticles"), #linked to reco::GenParticleCollection
	gensD0 = cms.InputTag("prunedGenParticles"), #linked to reco::GenParticleCollection
	#Weight
	PileupSumInfoInputTag = cms.InputTag("slimmedAddPileupInfo"), 
	mEventInfo = cms.InputTag("generator"),
	# Options
	comEnergy = cms.double(13000.),
	#TTBIt = cms.int32(34),
	debug = cms.untracked.bool(False),
	DstarSignificance3D = cms.double(3.0),
	D0Significance3D = cms.double(3.0),
	selectionCuts = cms.bool(False), #Apply Cuts Online
	triggerOn = cms.bool(False), #Apply Trigger
	triggerReferenceOn = cms.bool(False), #Apply Reference Trigger
	TracksOn = cms.bool(True)
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('D0DstarMC.root') #output files
)

#process.p = cms.Path(process.analysis)
process.p = cms.Path(process.primaryVertexFilter+process.analysis)
#process.p = cms.Path(process.primaryVertexFilter+process.noscraping+process.analysis)
