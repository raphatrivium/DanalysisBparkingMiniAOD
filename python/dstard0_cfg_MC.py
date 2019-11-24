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
#'file:A6DB8D2B-D159-E811-9EDA-0025905D1D02_MCminiAOD_noPU.root'
#'file:E8452CB8-2A97-E811-8826-0CC47A78A3B4_MCminiAOD_PU.root'
#'file:F4D6028F-1A27-CB49-BEA0-6771E00F049F_MCminiAOD_D0K3Pi.root'
'file:FD9CB40F-D40A-CE44-8FBC-82A8820D06F1_MinBias_TuneCP5_13TeV_pythia8.root'	

	   )
)
#===================================================================
process.analysis = cms.EDAnalyzer('DstarD0TTree',
	# Analysis
	doMC = cms.bool(True),
	doRec = cms.bool(True),
	bits = cms.InputTag("TriggerResults","","HLT"),
	prescales = cms.InputTag("patTrigger"),#linked to pat::PackedTriggerPrescales
	objects = cms.InputTag("selectedPatTrigger"),
	#PathName = cms.untracked.string("HLT_Mu8p5_IP3p5_part0"), #triggerName
	#PathName = cms.untracked.string("HLT_Mu12_IP6_part0"),  #triggerName
	PathName = cms.untracked.string("HLT_Mu10p5_IP3p5"),  #triggerName
	tracks = cms.InputTag('packedPFCandidates'),#linked to vector<pat::PackedCandidate
	recVtxs = cms.InputTag('offlineSlimmedPrimaryVertices'), #linked to vector<reco::Vertex> 
	gens = cms.InputTag("prunedGenParticles"), #linked to reco::GenParticleCollection
	gensD0 = cms.InputTag("prunedGenParticles"), #linked to reco::GenParticleCollection
	# Options
	comEnergy = cms.double(13000.),
	#TTBIt = cms.int32(34),
	debug = cms.untracked.bool(False),
	DstarSignificance3D = cms.double(1.0),
	D0Significance3D = cms.double(1.0),
	selectionCuts = cms.bool(False), #Apply Cuts Online
	triggerOn = cms.bool(False) #Apply Trigger Selection
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('D0DstarMC.root') #output files
)

process.p = cms.Path(process.analysis)
