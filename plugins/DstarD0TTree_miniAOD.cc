// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <map>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector.h"

#include <TString.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <iostream>
#include "TMath.h"
#include "TTree.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoVertex/KalmanVertexFit/interface/SingleTrackVertexConstraint.h"

//Weighting PileUp
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//GEN MC Matching
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

//JetCollection (CaloTower) - Gapside
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

//-------------------------------------------------------------
//#include "DataFormats/CaloTowers/interface/CaloTower.h"
//#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
//--------------------------------------------------------------

//Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DstarD0TTree_miniAOD.h"

//Triggers
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "HLTrigger/HLTanalyzers/plugins/HLTInfo.h"

// L1 related
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtTriggerMenuLite.h"

#include "TLorentzVector.h"

//MiniAOD
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


DstarD0TTree::DstarD0TTree(const edm::ParameterSet& iConfig):
	doMC(iConfig.getParameter<bool>("doMC")),
	doRec(iConfig.getParameter<bool>("doRec")),
	debug(iConfig.getUntrackedParameter<bool>("debug",false)),
	selectionCuts(iConfig.getParameter<bool>("selectionCuts")),
	triggerOn(iConfig.getParameter<bool>("triggerOn")),
	triggerReferenceOn(iConfig.getParameter<bool>("triggerReferenceOn")),
	TracksOn(iConfig.getParameter<bool>("TracksOn")),
	triggerName_(iConfig.getUntrackedParameter<std::string>("PathName","HLT_Mu9_IP6_part0_v2")),
	triggerReferenceName_(iConfig.getUntrackedParameter<std::string>("ReferencePathName","HLT_Mu9_IP6_part0_v2")),
	triggerReferenceName2_(iConfig.getUntrackedParameter<std::string>("ReferencePathName2","HLT_Mu9_IP6_part0_v2")),
	//Weighting PileUp and Gen 
	PileupSumInfoInputTag_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PileupSumInfoInputTag"))), 
	LumiWeights_(0), //
	mEventInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("mEventInfo"))),
 	//Triggers
	triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("objects"))),
	triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
	//Tracks and Vertex
	trkToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("tracks"))), //MiniAOD
	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("recVtxs"))),
	//Gen
	genParticlesTokenDstar_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gens"))),
	genParticlesTokenD0_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gensD0"))),
	//Values
	comEnergy_(iConfig.getParameter<double>("comEnergy")),
	DstarSignificance3D_(iConfig.getParameter<double>("DstarSignificance3D")),
	D0Significance3D_(iConfig.getParameter<double>("D0Significance3D"))
{     
	Total_Events_test = 0;
	Triggered_Event_test = 0;
	TriggeredReference_Event_test = 0;
	TriggeredReference_Event_test2 = 0;
	B0test = 0;
	Ebeam_ = comEnergy_/2.;
	edm::Service<TFileService> fs;

	//Creating trees
	data = fs->make<TTree>("data","data");
	
	// Pile-up re-weighting setup:
	//const edm::InputTag& PileupSumInfoInputTags = edm::InputTag("slimmedAddPileupInfo");
	//LumiWeights_ = edm::LumiReWeighting("PileupMC.root", "MyDataPileupHistogram.root", "input_Event/N_TrueInteractions", "pileup", PileupSumInfoInputTags);

}

DstarD0TTree::~DstarD0TTree()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//



//***********************************************************************************
// ------------ method called for each event  ------------
void DstarD0TTree::analyze(const edm::Event& iEvent, 
									const edm::EventSetup& iSetup)
{
	pi_mass=0.13957061; k_mass=0.493677;

	Total_Events_test++; 
	//To clear and initialize variables
	initialize();

	//run, event, lumi section
	runNumber= iEvent.id().run();
	eventNumber= iEvent.id().event();
	lumi = iEvent.luminosityBlock();

	Handle<double> lumiWeight;
	iEvent.getByLabel("lumiWeight",lumiWeight);
	lumiWeight_= lumi;

	if(doMC){
	//-------------------------------------------------------------
	//Weighting PileUp - Look in BeginJob for Files(Data, MC) input
	//-------------------------------------------------------------
	Handle<std::vector<PileupSummaryInfo>> PupInfo;
	iEvent.getByToken(PileupSumInfoInputTag_, PupInfo);
	const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
	if (LumiWeights_) PUWeight = LumiWeights_->weight( (*iEventB) );
	std::vector<PileupSummaryInfo>::const_iterator PVI;
	int npv = -1; // True number of primary vertices
	for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
	{	int BXs = PVI->getBunchCrossing();
		if(BXs == 0)
		{	npv = PVI->getTrueNumInteractions();	
			continue;
		}
	}
	nPU = npv;
	//cout << "PUWeight: " << PUWeight << endl;
	//-------------------------------------------------------------
	//WEIGHTING GEN
	//-------------------------------------------------------------
	Handle<GenEventInfoProduct> hEventInfo;
	iEvent.getByToken(mEventInfo_, hEventInfo);
	GenWeight = hEventInfo->weight();
	//cout << "hEventInfo->weight(): "<< GenWeight << endl;
	}
	//-------------------------------------------------------------
	//TRIGGERS
	//-------------------------------------------------------------
	Handle<edm::TriggerResults> triggerBits;
	iEvent.getByToken(triggerBits_, triggerBits);

	//Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
	//iEvent.getByToken(triggerObjects_, triggerObjects);

	Handle<pat::PackedTriggerPrescales> triggerPrescales;
	iEvent.getByToken(triggerPrescales_, triggerPrescales);
	//-------------------------------------------------------------
	//PRIMARY VERTEX COLLECTION	  
	//-------------------------------------------------------------
	Handle<VertexCollection> recVtxs;
	iEvent.getByToken(vtxToken_,recVtxs);

	Total_Events++;

	//Only events in which the path actually fired had stored the filter results and products:    
   bool triggerFired = TriggerInfo(iEvent,triggerBits,triggerPrescales,triggerName_);
   if(triggerFired) Triggered_Event_test++;
	bool triggerReferenceFired = TriggerInfo(iEvent,triggerBits,triggerPrescales,triggerReferenceName_);
   if(triggerReferenceFired) TriggeredReference_Event_test++;
	bool triggerReferenceFired2 = TriggerInfo(iEvent,triggerBits,triggerPrescales,triggerReferenceName2_);
	if(triggerReferenceFired2) TriggeredReference_Event_test2++;

	// Getting tracks from vert(ex)ices
	Handle< View < pat::PackedCandidate >> tracks; //access that PackedCandidate
	iEvent.getByToken(trkToken_,tracks);

	//miniAOD
	edm::ESHandle<TransientTrackBuilder> theB; 
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

	//process_pileup(iEvent, LumiWeights_, PileupSumInfoInputTag_);
	
	//To do the analisys with trigger or without trigger
	//if( (triggerComb == 0 and triggerOn and triggerFired  ) || //With trigger
   // (triggerComb == 1 and triggerOn )                 )  )  //Without Trigger
   // {canDo = true; }

	//Flag - trigger or without trigger
	bool canDo= false;
	for( int triggerComb = 0; triggerComb<3; triggerComb++) //For trigger selection on and off
	{
		canDo = false;//initiate the bool each loop.	
		//With trigger
		if (triggerComb == 0 and triggerOn and !triggerReferenceOn and triggerFired ){canDo = true;}
		//Without Trigger
		else if (triggerComb == 1 and !triggerOn and !triggerReferenceOn) {canDo = true; cout << "No Trigger Applied " << endl;	}
		// Trigger or Reference trigger
		else if (triggerComb == 2 and triggerOn and triggerReferenceOn and (triggerFired or triggerReferenceFired or triggerReferenceFired2) ) {canDo = true; /*cout << "Two Trigger Applied " << endl;*/}
		//security
		else{ canDo = false; }
	if (canDo){
			
	Triggered_Event++;
	NameOfFiredTriggers = TriggersFired(iEvent, triggerBits, triggerPrescales);
	
	//Loop in the trakcs of the PackedCandidate 
	for(View<pat::PackedCandidate>::const_iterator iTrack1 = tracks->begin(); iTrack1 != tracks->end(); ++iTrack1 ) 
	{   
		CounterTotalTracks++;
		if(!iTrack1->hasTrackDetails()) continue;
		CounterTracksHasTrackDetails++;
		if(!(iTrack1->trackHighPurity())) continue;
		CounterTracksHighPurity++;

		if(TracksOn){
		TracksCharge.push_back(iTrack1->charge());
		TracksEta.push_back(iTrack1->eta());
		TracksPt.push_back(iTrack1->pt());
		TracksPhi.push_back(iTrack1->phi());
		TracksChi2.push_back(iTrack1->pseudoTrack().normalizedChi2());
		TracksNumberOfHits.push_back(iTrack1->numberOfHits());
		TracksNumberOfPixelHits.push_back(iTrack1->numberOfPixelHits());
		Tracksdxy.push_back(iTrack1->dxy());
		Tracksdz.push_back(iTrack1->dz());
		TracksdxyError.push_back(iTrack1->dxyError());
		TracksdzError.push_back(iTrack1->dzError());
		}

		if(iTrack1->charge()==0) continue;
		CounterTracksCharged++;
		if(fabs(iTrack1->eta())>2.4) continue;
		CounterTracksEta++;
		
		//if(fabs(iTrack1->pdgId()) != 211 ) continue;
		//TracksPDG211++;
		// SELECTING TRACKS FOR D0 Prompt   
		if( iTrack1->pt()>0.8 /*&& (fabs(iTrack1->pdgId()) == 211)*/)
		{
			CounterTracksPt0p8++;
			if( fabs(iTrack1->eta()) < 2.4 ){
				CounterD0PromptTracksEta2p5++;	
				if( iTrack1->pseudoTrack().normalizedChi2() < 5.0 ) {
					CounterD0PromptTracksChi5++;			
					//if( iTrack1->pseudoTrack().hitPattern().numberOfValidHits() <= 5 ) continue;
					if( iTrack1->numberOfHits() > 5) {
						CounterD0PromptTracksNumberHits5++;
						//if( iTrack1->pseudoTrack().hitPattern().numberOfValidPixelHits() <= 2 ) continue;
						if( iTrack1->numberOfPixelHits() > 2) {
							CounterD0PromptTracksPixelHits2++;
							//if( iTrack1->p() <1.0 ) continue;
							//D0PromptTracksMomentum1++;
							if( fabs(iTrack1->dxy()) < 0.1 ) {
								CounterD0PromptTracksDxy0p1++;
								if( fabs(iTrack1->dz()) < 1. ) {
									CounterD0PromptTracksDz0p5++;
									reco::TransientTrack  D0TT = theB->build(iTrack1->pseudoTrack());
									goodTracksD0.push_back(D0TT);//Fill Transient Vector
								}
							}
						}
					}
				}
			}
		}//End SELECTING TRACKS FOR D0      

		//Selection for slowpion from D* (The note is 0.5 but the MiniAOD works with pt > 0.5
		if( iTrack1->pt()>0.5 )
		{
			CounterTracksPtZeroFive++;
			if(iTrack1->pseudoTrack().normalizedChi2() > 3.) continue;				  
			CounterTracksChi3++;
			//if(iTrack1->pseudoTrack().hitPattern().numberOfValidHits() < 2) continue;
			if(iTrack1->numberOfHits() < 2) continue; // > and = 2
			CounterTracksNumberOfHits2++;
			if( ( fabs(iTrack1->dxy()) / fabs(iTrack1->dxyError()) ) >3. ) continue;			
			CounterTracksDxyThree++;
			if( ( fabs(iTrack1->dz()) / fabs(iTrack1->dzError()) ) >3. ) continue;
			CounterTracksDzThree++;
			reco::TransientTrack slowpionTT = theB->build(iTrack1->pseudoTrack());
			slowPiTracks.push_back(slowpionTT);	//Fill Transient Vector
		}//End Selection for slowpion from D*
			
		//Kaon and Pion Candidates for D*
		if( iTrack1->pt()>0.5 /*&& (fabs(iTrack1->pdgId()) == 211)*/)
		{
			CounterTracksPtZeroSix++;
			if(iTrack1->pseudoTrack().normalizedChi2() > 2.5) continue;
			CounterTracksChiTwoFive++;
			//cout << "iTrack1->pseudoTrack().normalizedChi2(): "<< iTrack1->pseudoTrack().normalizedChi2() << endl ;
			//cout << "iTrack1->pseudoTrack().Chi2(): "<< iTrack1->pseudoTrack().chi2() << endl ;		
			//if(iTrack1->pseudoTrack().hitPattern().numberOfValidHits() < 5) continue;
			//if(iTrack1->pseudoTrack().hitPattern().numberOfValidPixelHits() < 2) continue;
			if( iTrack1->numberOfHits() < 5) continue;
			CounterTracksNumberOfHits5++;
			if( iTrack1->numberOfPixelHits() < 2) continue;
			CounterTracksNumberOfPixelHits2++;				
			if( fabs(iTrack1->dxy()) > 0.1) continue;
			CounterTracksDxyZeroOne++;	
			if( fabs(iTrack1->dz()) > 1.) continue;
			CounterTracksDzOne++;
			reco::TransientTrack  PionTT = theB->build(iTrack1->pseudoTrack());
			if (debug)cout << " PionTT "  << PionTT.track().momentum() << endl;
			goodTracks.push_back(PionTT); //Fill Transient Vector
		}//End Kaon and Pion Candidates for D* 
			
   }//loop packed candidates
	
	if (debug) cout << " goodTracks size " << goodTracks.size() << endl;
 	ntracksDstar = slowPiTracks.size();
 	ntracksD0Kpi = goodTracksD0.size();

	//Vertex Informations
	size_t vtx_trk_size = (*recVtxs)[0].tracksSize();
	int VtxIn=0;
	for(size_t i = 0; i < recVtxs->size(); ++ i)
	{
		const Vertex &vtx = (*recVtxs)[i];
		if(vtx.tracksSize()>vtx_trk_size)
		{ VtxIn = i; vtx_trk_size = vtx.tracksSize(); }
	}

	const Vertex& RecVtx = (*recVtxs)[VtxIn];
	n_pVertex = recVtxs->size();
	PVx = RecVtx.x();
	PVy = RecVtx.y();
	PVz = RecVtx.z();
	PVerrx=RecVtx.xError();
	PVerry=RecVtx.yError();
	PVerrz=RecVtx.zError();

	//Calling Secondary Functions
	RecDstar(iEvent,iSetup,RecVtx); //Reconstruction of D*
	RecDstarWrongCombination(iEvent,iSetup,RecVtx); //Reconstruction of D* with wrong charge combination (for BG model)
	GenDstarInfo(iEvent,iSetup); //Stores information from D0 and its products.
	GenDstarMesonBInfo(iEvent,iSetup);
	RecD0(iEvent,iSetup,RecVtx); //Reconstruction of prompt D0
	GenD0Info(iEvent,iSetup); //Stores information from D0 and its products.
	GenD0MesonBInfo(iEvent,iSetup);
	//FindAngle(RecVtx,v_D0,d0kpi_p4); //Calculates the opening angle

	//Fill the tree
	data->Fill();

	}//End canDo	
	}//end trigger flag
}//End analyze
//*********************************************************************************
bool DstarD0TTree::TriggerInfo(const edm::Event& iEvent, 
										 edm::Handle<edm::TriggerResults> itriggerBits, 
										 edm::Handle<pat::PackedTriggerPrescales> itriggerPrescales,
										 TString trigname)
{
	const edm::TriggerNames &names = iEvent.triggerNames(*itriggerBits);
	
	//List of Trigger Names in the Dataset
	//for (unsigned int k = 0, nt = itriggerBits->size(); k < nt; ++k) 
	//{std::cout << "Trigger name: "<< names.triggerName(k) << std::endl;}

	bool foundtrigger = false;
	//std::cout << "\n == TRIGGER PATHS FOUND =" ;
	for (unsigned int i = 0, n = itriggerBits->size(); i < n; ++i) 
	{		
		TString trigName = names.triggerName(i);
		if(trigName.Contains(trigname)) 
		{  foundtrigger = true;
			//std::cout << "TRIGGER PATHS FOUND: " << names.triggerName(i) <<
			//", prescale " << itriggerPrescales->getPrescaleForIndex(i) <<
			//": " << (itriggerBits->accept(i) ? "PASS" : "fail (or not run) --") << endl;
			
			if (itriggerBits->accept(i) )
			{	
				//std::cout << "itriggerBits->accept(i): " << itriggerBits->accept(i) << std::endl; std::cout << "---" << std::endl;
				return true;}
		} 		
	}
		    	
	if (!foundtrigger) {std::cout << "**TRIGGER PATHS NOT FOUND**" << std::endl; std::cout << "---" << std::endl;}
	return false;
}
//*********************************************************************************
std::vector<std::string> DstarD0TTree::TriggersFired(const edm::Event& iEvent, 
																	  edm::Handle<edm::TriggerResults> itriggerBits,
	 																  edm::Handle<pat::PackedTriggerPrescales> itriggerPrescales)
{	//std::cout << "--------------: " << std::endl;
	const edm::TriggerNames &names = iEvent.triggerNames(*itriggerBits);
	std::vector<std::string> Name_Trigger; Name_Trigger.clear();
	std::string trigName;
	for (unsigned int i = 0, n = itriggerBits->size(); i < n; ++i) 
	{	
		trigName = names.triggerName(i);
		if(itriggerBits->accept(i)) //Check if the trigger fired
		{ 	Name_Trigger.push_back(trigName); //Fill Vector
			//std::cout << "TriggersFired: " << trigName <<
			//", prescale " << itriggerPrescales->getPrescaleForIndex(i) << std::endl;
		} 				 		
	}		    	
	return Name_Trigger;
}
//***********************************************************************************
/*void DstarD0TTree::process_pileup(const edm::Event& iEvent, 
												edm::LumiReWeighting weights, 
												edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupInfo){
	float tnpv = -1;
	float wpu = 1;
	Handle<std::vector<PileupSummaryInfo>> info;
	iEvent.getByToken(pileupInfo, info);
	for (std::vector<PileupSummaryInfo>::const_iterator pvi = info->begin(); pvi != info->end(); ++pvi)
	{	int bx = pvi->getBunchCrossing();
		if (bx == 0) 
		{	tnpv = pvi->getTrueNumInteractions();
			continue;
		}
	}
	wpu = weights.weight(tnpv);	
	cout << "wpu  = " << wpu << endl;
	cout << "tnpv = " << tnpv << endl;
	//branches["event"]["wpu"].push_back(wpu);
	//branches["event"]["tnpv"].push_back(tnpv);
}*/
//***********************************************************************************
void DstarD0TTree::RecDstar(const edm::Event& iEvent, 
									 const edm::EventSetup& iSetup, 
									 const reco::Vertex& RecVtx){
   if( goodTracks.size() >=2 ){
	for(size_t i=0; i<(goodTracks.size()-1);i++)
	{  
		TransientTrack trk1 = goodTracks[i];
		for(size_t j=i+1;j<goodTracks.size();j++)
		{			    
			TransientTrack trk2 = goodTracks[j];
			if(trk1.charge() == trk2.charge()) continue;
			for(size_t k=0;k<slowPiTracks.size();k++)
			{
				TransientTrack trkS = slowPiTracks[k];
				if(trkS == trk1 || trkS == trk2) continue;			
				TransientTrack K;	TransientTrack pi;
                                
				//Right Combination Charges
				if(trk1.charge() == trkS.charge()){pi = trk1; K = trk2;}
				else{	K = trk1; pi = trk2;}

				//Wrong Combination Charges - background
				/*if(trk1.charge() == trkS.charge()){pi = trk2;	K = trk1;}
				else{K = trk2;pi = trk1;}*/
                               
				//D0 4-momentum Reconstruction (kaon and pion)
				math::XYZTLorentzVector ip4_K(K.track().px(),
														K.track().py(),
														K.track().pz(),
														sqrt(pow(K.track().p(),2)+pow(k_mass,2)) );
				math::XYZTLorentzVector ip4_pi(	pi.track().px(),
															pi.track().py(),
															pi.track().pz(),
															sqrt(pow(pi.track().p(),2)+pow(pi_mass,2)) );        
				math::XYZTLorentzVector ip4_D0 = ip4_K + ip4_pi;
				CounterD0AfterLorentzVector++;

				//SlowPion 4-momentum Reconstruction
				math::XYZTLorentzVector p4_S(trkS.track().px(),
														trkS.track().py(),
														trkS.track().pz(),
														sqrt(pow(trkS.track().p(),2)+pow(pi_mass,2)));

				//D* 4-momentum Reconstruction (D0 + SlowPion)
				math::XYZTLorentzVector ip4_DS = ip4_D0 + p4_S;
		
				//To descrase statistc to gain speed. we have more tight cut below
				if( fabs(ip4_D0.M()-1.86484)  > 0.2) continue;
				CounterD0MinusPDGzero2++;

				//To descrase statistc to gain speed. we have more tight cut below
				if((ip4_DS.M() - ip4_D0.M()) > 0.25) continue;
 				//cout << "ip4_DS.M() :" << ip4_DS.M() << " ip4_D0.M(): " << ip4_D0.M() << endl;
				CounterDsMinusD0Zerothree++;
				                      
				//KalmanVertexFitter
				vector<TransientTrack> tks;
				tks.push_back(K);
				tks.push_back(pi);
				KalmanVertexFitter kalman(true);
				TransientVertex v = kalman.vertex(tks);
				if(!v.isValid() || !v.hasRefittedTracks()) continue;

				//Vertex Probability
				double vtxProb =TMath::Prob( (Double_t) v.totalChiSquared(), (Int_t) v.degreesOfFreedom());

				TransientTrack K_f = v.refittedTrack(K);
				TransientTrack pi_f = v.refittedTrack(pi);
				CounterTransientTrackOfpiK++;
				
				//SV Confidence Level
				if(vtxProb < 0.01) continue;
				CounterSVConfidenceLevel++;
							
				//D* 4-momentum Reconstruction after KalmanVertexFitter
				math::XYZTLorentzVector p4_K(K_f.track().px(),
														K_f.track().py(),
														K_f.track().pz(),
														sqrt(pow(K_f.track().p(),2)+pow(k_mass,2)));
				math::XYZTLorentzVector p4_pi(pi_f.track().px(),
														pi_f.track().py(),
														pi_f.track().pz(),
														sqrt(pow(pi_f.track().p(),2)+pow(pi_mass,2)));
				math::XYZTLorentzVector d0_p4 = p4_K + p4_pi;
				double d0mass = d0_p4.M();
				CounterD0AfterLorentzVectorKarman++;
				
				//Angle between formation and decay of D0
  				double anglephi = FindAngle(RecVtx,v,d0_p4);
				double cosPhi = cos(anglephi);
				if( cosPhi < 0.99 ) continue;
				CounterPointingcosPhi++;

				//D0(from D*) Siginificance xy
				VertexDistanceXY vD0fromDSdXY ;
				double D0fromDSdXY = vD0fromDSdXY.distance(RecVtx,v).value() ;
				double D0fromDSeXY = vD0fromDSdXY.distance(RecVtx,v).error() ;
				double D0fromDSsXY =  D0fromDSdXY / D0fromDSeXY;

				//D0(from D*) Siginificance 3D
				VertexDistance3D vD0Kpid3D ;
				double D0fromDSd3D = vD0Kpid3D.distance(RecVtx,v).value() ;
				double D0fromDSe3D = vD0Kpid3D.distance(RecVtx,v).error() ;
				double D0fromDSs3D = D0fromDSd3D / D0fromDSe3D;

				if(selectionCuts) {if ( D0fromDSs3D < DstarSignificance3D_ ) continue;}
				if(selectionCuts) {CounterSignificance3D++;}

				//Difference between D0 and D0PDG
				if( fabs(d0mass - 1.86484) > 0.1 ) continue;
				CounterD0MinusPDG++;

				//Pt Cut of mesons D0
				if(selectionCuts) { if ( d0_p4.Pt() < 3. ) continue;}
				if(selectionCuts) { CounterD0pTThree++;}
				
				//Definition of the D* quadri-vector					
				math::XYZTLorentzVector dS_p4 = d0_p4 + p4_S;
				double dsmass = dS_p4.M();
				CounterDsAfterLorentzVector++;
								                      											
				//Difference between D* and D0 -> Must be close to pion
				if(selectionCuts) {if( (dsmass - d0mass) > 0.16) continue;}
				if(selectionCuts) {CounterDsMinusD0++;}
				CounterDsCandidates++; //Number of D* candidates
						                       		
				//Fill vectors
				Dsmass.push_back(dS_p4.M());
				Dspt.push_back(dS_p4.Pt());
				Dseta.push_back(dS_p4.eta());
				Dsphi.push_back(dS_p4.phi());

				D0mass.push_back(d0_p4.M());
				D0pt.push_back(d0_p4.Pt());
				D0eta.push_back(d0_p4.eta());
				D0phi.push_back(d0_p4.phi());
				D0_VtxProb.push_back(vtxProb);
				D0_VtxPosx.push_back(v.position().x());
				D0_VtxPosy.push_back(v.position().y());
				D0_VtxPosz.push_back(v.position().z());
				D0_Vtxerrx.push_back(v.positionError().cxx());
				D0_Vtxerry.push_back(v.positionError().cyy());
				D0_Vtxerrz.push_back(v.positionError().czz());
				D0fromDSd3D_vec.push_back(D0fromDSd3D);
				D0fromDSe3D_vec.push_back(D0fromDSe3D);	
				D0fromDSs3D_vec.push_back(D0fromDSs3D);
				D0fromDSsXY_vec.push_back(D0fromDSsXY);
				D0fromDSdXY_vec.push_back(D0fromDSdXY);
				D0fromDSeXY_vec.push_back(D0fromDSeXY);
				Anglephi_vec.push_back(cosPhi);   

				TrkKmass.push_back(p4_K.M());
				TrkKdxy.push_back(K_f.track().dxy(RecVtx.position()));
				TrkKdz.push_back(K_f.track().dz(RecVtx.position()));
				TrkKnhits.push_back(K.track().numberOfValidHits());
				TrkKchi2.push_back(K.track().normalizedChi2());
				TrkKpt.push_back(K_f.track().pt());
				TrkKeta.push_back(K_f.track().eta());
				TrkKphi.push_back(K_f.track().phi());

				Trkpimass.push_back(p4_pi.M());
				Trkpidxy.push_back(pi_f.track().dxy(RecVtx.position()));
				Trkpidz.push_back(pi_f.track().dz(RecVtx.position()));
				Trkpinhits.push_back(pi.track().numberOfValidHits());
				Trkpichi2.push_back(pi.track().normalizedChi2());
				Trkpipt.push_back(pi_f.track().pt());
				Trkpieta.push_back(pi_f.track().eta());
				Trkpiphi.push_back(pi_f.track().phi());

				TrkSmass.push_back(p4_S.M());
				TrkSdxy.push_back(trkS.track().dxy(RecVtx.position()));
				TrkSdz.push_back(trkS.track().dz(RecVtx.position()));
				TrkSnhits.push_back(trkS.track().numberOfValidHits());
				TrkSchi2.push_back(trkS.track().normalizedChi2());
				TrkSpt.push_back(trkS.track().pt());
				TrkSeta.push_back(trkS.track().eta());
				TrkSphi.push_back(trkS.track().phi());
				TrkScharge.push_back(trkS.charge());

				DSDeltaR.push_back(deltaR(d0_p4.eta(),d0_p4.phi(),trkS.track().eta(),trkS.track().phi()));
                       
				NKpiCand++;  
									
				if(NKpiCand>999) break;// this a control counter. We do not expect high rate of D*
			}
			if(NKpiCand>999) break;
	 } 
	if(NKpiCand>999) break;
	}
	}
}//End RecDstar

//***********************************************************************************
void DstarD0TTree::RecDstarWrongCombination(const edm::Event& iEvent,
														  const edm::EventSetup& iSetup, 
														  const reco::Vertex& RecVtx){	
   if( goodTracks.size() >=2 ){
	for(size_t i=0;i<(goodTracks.size()-1);i++)
	{  
		TransientTrack trk1 = goodTracks[i];

		for(size_t j=i+1;j<goodTracks.size();j++)
		{			    
			TransientTrack trk2 = goodTracks[j];
			//cout << goodTracks[i].track().momentum() << goodTracks[j].track().momentum() << endl;

			if(trk1.charge() == trk2.charge()) continue;
			//cout <<   trk1.track().momentum() <<    trk2.track().momentum() << endl;

			for(size_t k=0;k<slowPiTracks.size();k++)
			{
				TransientTrack trkS = slowPiTracks[k];
				//cout << "Slow Pions Tracks(px,py,pz): "  << trkS.track().momentum() << endl;                
				if(trkS == trk1 || trkS == trk2) continue;
				
				TransientTrack K;	TransientTrack pi;
                                
				//Wrong Combination Charges - background
				if(trk1.charge() == trkS.charge())	{ pi = trk2; K = trk1;}
				else{K = trk2;	pi = trk1;}
                               
				//D0 4-momentum Reconstruction (kaon and pion)
				math::XYZTLorentzVector ip4_K(K.track().px(),
														K.track().py(),
														K.track().pz(),
														sqrt(pow(K.track().p(),2)+pow(k_mass,2)) );
				math::XYZTLorentzVector ip4_pi(	pi.track().px(),
															pi.track().py(),
															pi.track().pz(),
															sqrt(pow(pi.track().p(),2)+pow(pi_mass,2)) );        
				math::XYZTLorentzVector ip4_D0 = ip4_K + ip4_pi;
				//CounterD0AfterLorentzVector++;

				//SlowPion 4-momentum Reconstruction
				math::XYZTLorentzVector p4_S(trkS.track().px(),
														trkS.track().py(),
														trkS.track().pz(),
														sqrt(pow(trkS.track().p(),2)+pow(pi_mass,2)));

				//D* 4-momentum Reconstruction (D0 + SlowPion)
				math::XYZTLorentzVector ip4_DS = ip4_D0 + p4_S;
		
				//To descrase statistc to gain speed. we have more tight cut below
				if( fabs(ip4_D0.M()-1.86484)  > 0.2) continue;
				//CounterD0MinusPDGzero2++;

				//To descrase statistc to gain speed. we have more tight cut below
				if((ip4_DS.M() - ip4_D0.M()) > 0.2) continue;
 				//cout << "ip4_DS.M() :" << ip4_DS.M() << " ip4_D0.M(): " << ip4_D0.M() << endl;
				//CounterDsMinusD0Zerothree++;
				                      
				//KalmanVertexFitter
				vector<TransientTrack> tks;
				tks.push_back(K);
				tks.push_back(pi);
				KalmanVertexFitter kalman(true);
				TransientVertex v = kalman.vertex(tks);
				if(!v.isValid() || !v.hasRefittedTracks()) continue;

				//Vertex Probability
				double vtxProb =TMath::Prob( (Double_t) v.totalChiSquared(), (Int_t) v.degreesOfFreedom());

				TransientTrack K_f = v.refittedTrack(K);
				TransientTrack pi_f = v.refittedTrack(pi);
				//CounterTransientTrackOfpiK++;
				
				//SV Confidence Level
				if(selectionCuts) {if(vtxProb < 0.01) continue;}
				//CounterSVConfidenceLevel++;
				
				//D* 4-momentum Reconstruction after KalmanVertexFitter
				math::XYZTLorentzVector p4_K(K_f.track().px(),
														K_f.track().py(),
														K_f.track().pz(),
														sqrt(pow(K_f.track().p(),2)+pow(k_mass,2)));
				math::XYZTLorentzVector p4_pi(pi_f.track().px(),
														pi_f.track().py(),
														pi_f.track().pz(),
														sqrt(pow(pi_f.track().p(),2)+pow(pi_mass,2)));
				math::XYZTLorentzVector d0_p4 = p4_K + p4_pi;
				double d0mass = d0_p4.M();
				//CounterD0AfterLorentzVectorKarman++;
				
				//Angle between formation and decay of D0(from D*)
  				double anglephi = FindAngle(RecVtx,v,d0_p4);
				double cosPhi = cos(anglephi);
				if( cosPhi < 0.99 ) continue;
				//CounterPointingcosPhi++;

				//D0 from D* Siginificance
				VertexDistanceXY vD0fromDSdXY ;
				double D0fromDSdXY = vD0fromDSdXY.distance(RecVtx,v).value() ;
				double D0fromDSeXY = vD0fromDSdXY.distance(RecVtx,v).error() ;
				double D0fromDSsXY =  D0fromDSdXY / D0fromDSeXY;

				//D0 significance 3D
				VertexDistance3D vD0Kpid3D ;
				double D0fromDSd3D = vD0Kpid3D.distance(RecVtx,v).value() ;
				double D0fromDSe3D = vD0Kpid3D.distance(RecVtx,v).error() ;
				double D0fromDSs3D = D0fromDSd3D / D0fromDSe3D;

				//cout << "significance 3D: "<< D0fromDSs3D << endl;
				if(selectionCuts) {if( D0fromDSs3D < DstarSignificance3D_ ) continue;}
				//CounterSignificance3D++;

				//Difference between D0 and D0PDG 
				if( fabs(d0mass - 1.86484) > 0.1 ) continue;
				//CounterD0MinusPDG++;

				//Pt Cut of mesons D0
				if(selectionCuts) {if( d0_p4.Pt() < 3. ) continue;}
				//CounterD0pTThree++;
								
				math::XYZTLorentzVector dS_p4 = d0_p4 + p4_S;
				double dsmass = dS_p4.M();
				CounterDsAfterLorentzVector++;
								                      											
				//Difference between D* and D0 -> Must be close to pion
				if(selectionCuts) {if( (dsmass - d0mass) > 0.16) continue;}
				//CounterDsMinusD0++;
				//CounterDsCandidates++; //Number of D* candidates
			                       		
				//Fill vectors
				D0_VtxProbWrong.push_back(vtxProb);
				D0massWrong.push_back(d0_p4.M());
				DsmassWrong.push_back(dS_p4.M());
				D0ptWrong.push_back(d0_p4.Pt());
				DsptWrong.push_back(dS_p4.Pt());
				D0etaWrong.push_back(d0_p4.eta());
				D0phiWrong.push_back(d0_p4.phi());
				DsetaWrong.push_back(dS_p4.eta());
				DsphiWrong.push_back(dS_p4.phi());

				D0_VtxPosxWrong.push_back(v.position().x());
				D0_VtxPosyWrong.push_back(v.position().y());
				D0_VtxPoszWrong.push_back(v.position().z());
				D0_VtxerrxWrong.push_back(v.positionError().cxx());
				D0_VtxerryWrong.push_back(v.positionError().cyy());
				D0_VtxerrzWrong.push_back(v.positionError().czz());

				TrkKmassWrong.push_back(p4_K.M());
				TrkpimassWrong.push_back(p4_pi.M());
				TrkSmassWrong.push_back(p4_S.M());

				TrkKdxyWrong.push_back(K_f.track().dxy(RecVtx.position()));
				TrkpidxyWrong.push_back(pi_f.track().dxy(RecVtx.position()));
				TrkSdxyWrong.push_back(trkS.track().dxy(RecVtx.position()));

				TrkKdzWrong.push_back(K_f.track().dz(RecVtx.position()));
				TrkpidzWrong.push_back(pi_f.track().dz(RecVtx.position()));
				TrkSdzWrong.push_back(trkS.track().dz(RecVtx.position()));

				TrkKnhitsWrong.push_back(K.track().numberOfValidHits());
				TrkpinhitsWrong.push_back(pi.track().numberOfValidHits());
				TrkSnhitsWrong.push_back(trkS.track().numberOfValidHits());

				TrkKchi2Wrong.push_back(K.track().normalizedChi2());
				Trkpichi2Wrong.push_back(pi.track().normalizedChi2());
				TrkSchi2Wrong.push_back(trkS.track().normalizedChi2());

				DSDeltaRWrong.push_back(deltaR(d0_p4.eta(),d0_p4.phi(),trkS.track().eta(),trkS.track().phi()));

				TrkKptWrong.push_back(K_f.track().pt());
				TrkpiptWrong.push_back(pi_f.track().pt());
				TrkSptWrong.push_back(trkS.track().pt());

				TrkKetaWrong.push_back(K_f.track().eta());
				TrkpietaWrong.push_back(pi_f.track().eta());
				TrkSetaWrong.push_back(trkS.track().eta());

				TrkKphiWrong.push_back(K_f.track().phi());
				TrkpiphiWrong.push_back(pi_f.track().phi());
				TrkSphiWrong.push_back(trkS.track().phi());

				TrkSchargeWrong.push_back(trkS.charge());

				D0fromDSs3D_vecWrong.push_back(D0fromDSs3D);
				D0fromDSsXY_vecWrong.push_back(D0fromDSsXY);
				Anglephi_vecWrong.push_back(cosPhi);   
                               
				//NKpiCandWrong++;  
									
				if(NKpiCand>999) break;
			}
			if(NKpiCand>999) break;
	 } 
	if(NKpiCand>999) break;
	}
	} //Protection
}//RecDstarWrongCombination

//***********************************************************************************
//Evaluate if a particle is stable (status 1). If not it repeat the process until find the a stable one.
void DstarD0TTree::assignStableDaughters(const reco::Candidate* p,
													  std::vector<int> & pids){
	for(size_t i=0; i< p->numberOfDaughters(); i++)
	{
		if(p->daughter(i)->status() == 1) pids.push_back(abs(p->daughter(i)->pdgId()));
		else	assignStableDaughters(p->daughter(i), pids);
		//cout << p->daughter(i)->pdgId() << endl;
	}
	return;
}

//***********************************************************************************
void DstarD0TTree::GenDstarInfo(const edm::Event& iEvent,
										  const edm::EventSetup& iSetup){
	if(doMC) //MC flag - see config file
	{ 
		edm::Handle<GenParticleCollection> gens; //linked to "prunedGenParticles" in config file
		iEvent.getByToken(genParticlesTokenDstar_, gens);
		int index0,index1;
		int IDdau1,IDdau2;
		int index00,index11;

		for(size_t i=0;i<gens->size();i++)
		{
     		const GenParticle & p = (*gens)[i];
			// D*
			if( fabs(p.pdgId()) != 413 ) continue;
			// 2 daughters
			if( p.numberOfDaughters()!=2 ) continue;		  
			//2 daughters and the first is D0         					
			if( fabs(p.daughter(0)->pdgId()) == 421 
			&& p.daughter(0)->numberOfDaughters()==2 
			&& fabs(p.daughter(1)->pdgId()) == 211 ){ index0 =0; index1 =1; }
			//2 daughters and the first is Pion
			else if( fabs(p.daughter(0)->pdgId()) == 211 
			&& p.daughter(1)->numberOfDaughters()==2 
			&& fabs(p.daughter(1)->pdgId()) == 421 ){ index0 =1; index1 =0; }
			else continue;

			const Candidate* dau = p.daughter(index0);
			const Candidate* slow = p.daughter(index1);
			IDdau1 = dau->daughter(0)->pdgId(); 
			IDdau2 = dau->daughter(1)->pdgId();
			//Fisrt Particle is a Kaon
			if(  (fabs(IDdau1) == 321) 
			&&	  (fabs(IDdau2) == 211) 
			&&	  (IDdau1*IDdau2 < 0  ) ){index00 =0; index11 =1;}
			//Fisrt Particle is a Pion
			else if( (fabs(IDdau1) == 321) 
			&&			(fabs(IDdau2) == 211) 
			&&			(IDdau1*IDdau2 < 0)   ){index00 =1; index11 =0;}
			else continue;
			
			//Dstar Quantities
			MCDseta.push_back(p.eta());
			MCDsphi.push_back(p.phi());
			MCDspt.push_back(p.pt());
			MCDsenergy.push_back(p.energy());
			MCDsp.push_back(p.p());
			MCDset.push_back(p.et());
			MCDsrapidity.push_back(p.rapidity());
			MCDsmass.push_back(p.mass());
			//D0 Quantities
			MCD0eta.push_back(dau->eta());
			MCD0phi.push_back(dau->phi());
			MCD0pt.push_back(dau->pt());
			MCD0energy.push_back(dau->energy());
			MCD0p.push_back(dau->p());
			MCD0et.push_back(dau->et());
			MCD0rapidity.push_back(dau->rapidity());   
			MCD0mass.push_back(dau->mass());
			//Find D0 displacement and Lifetime
			double D0Lxy = sqrt( pow(dau->vx()-dau->daughter(0)->vx(),2) + 
									   pow(dau->vy()-dau->daughter(0)->vy(),2) );
			double D0lifeTime = dau->mass() * D0Lxy*pow(10,-2) / (c * dau->pt() );
			//cout << "MC D0(fromD*) displacement-xy[cm]: " << D0Lxy << " # LifeTime [s]: " << D0lifeTime << endl;
			MCD0dispXY.push_back(D0Lxy);
			MCD0lifetime.push_back(D0lifeTime);
			//SlowPion Quantities
			MCDsSeta.push_back(slow->eta());
			MCDsSphi.push_back(slow->phi());
			MCDsSpt.push_back(slow->pt());
			MCDsSenergy.push_back(slow->energy());
			MCDsSp.push_back(slow->p());
			MCDsSet.push_back(slow->et());
			MCDsSrapidity.push_back(slow->rapidity());
			MCDsSmass.push_back(slow->mass());					
			//Kaon Quantities
			MCDsKeta.push_back(dau->daughter(index00)->eta());
			MCDsKphi.push_back(dau->daughter(index00)->phi());
			MCDsKpt.push_back(dau->daughter(index00)->pt());
			MCDsKenergy.push_back(dau->daughter(index00)->energy());
			MCDsKp.push_back(dau->daughter(index00)->p());
			MCDsKet.push_back(dau->daughter(index00)->et());
			MCDsKrapidity.push_back(dau->daughter(index00)->rapidity());
			MCDsKmass.push_back(dau->daughter(index00)->mass());
			//Pion Quantities
			MCDsPieta.push_back(dau->daughter(index11)->eta());
			MCDsPiphi.push_back(dau->daughter(index11)->phi());
			MCDsPipt.push_back(dau->daughter(index11)->pt());
			MCDsPienergy.push_back(dau->daughter(index11)->energy());
			MCDsPip.push_back(dau->daughter(index11)->p());
			MCDsPiet.push_back(dau->daughter(index11)->et());
			MCDsPirapidity.push_back(dau->daughter(index11)->rapidity());
			MCDsPimass.push_back(dau->daughter(index11)->mass());			 		    	
    	}// loop Gen Particle
   } //MC flag
}//End GenDstarInfo

//***********************************************************************************
void DstarD0TTree::GenDstarMesonBInfo(const edm::Event& iEvent,
										  const edm::EventSetup& iSetup){
	if(doMC) //MC flag - see config file
	{ 
		edm::Handle<GenParticleCollection> gens; //"prunedGenParticles" module in config file
		iEvent.getByToken(genParticlesTokenDstar_, gens);
		int index0, index1;
		int IDdau1, IDdau2;
		int index00, index11;
		for(size_t i=0;i<gens->size();i++)
		{		
     		const GenParticle & p = (*gens)[i];
			//B0 or B+-
      	if( (fabs(p.pdgId()) != 511) && (fabs(p.pdgId()) != 521) ) continue; 			 
     		for(size_t j=0; j < p.numberOfDaughters(); j++)
			{	
          	const Candidate* dau = p.daughter(j);
				// D* 
				if( fabs(dau->pdgId()) != 413) continue;
				// 2 daughtes
				if( dau->numberOfDaughters() !=2 ) continue;
    			// If the first particle is a D0
				B0test++;
				if( ( fabs(dau->daughter(0)->pdgId()) == 421  ) 
				&&  ( dau->daughter(0)->numberOfDaughters()==2) 
				&&  ( fabs(dau->daughter(1)->pdgId()) == 211  ) ){ index0 =0; index1 =1; }
				// If the first particle is a Pion
				else if( (fabs(dau->daughter(0)->pdgId()) == 211  )
				&& 		(dau->daughter(1)->numberOfDaughters()==2) 
				&& 		(fabs(dau->daughter(1)->pdgId()) == 421  ) ) { index0 =1; index1 =0; }
				else continue;
					
				const Candidate* granddau = dau->daughter(index0);
				const Candidate* slow = dau->daughter(index1);
				IDdau1 = granddau->daughter(0)->pdgId(); 
				IDdau2 = granddau->daughter(1)->pdgId();
				// If the first particle is a Kaon
				if( (fabs(IDdau1) == 321 ) 
				&&  (fabs(IDdau2) == 211 ) 
				&&  (IDdau1*IDdau2 < 0   ) )  { index00 =0; index11 =1;}
				// If the first particle is a Pion
				else if( ( fabs(IDdau1) == 211 )
				&& 		( fabs(IDdau2) == 321 ) 
				&& 		( (IDdau1*IDdau2) < 0 ) ){ index00 =1; index11 =0;}
				else continue;
				
				// Dstar Quantities
				MCDstarB_eta.push_back(dau->eta());
				MCDstarB_phi.push_back(dau->phi());
				MCDstarB_pt.push_back(dau->pt());
				MCDstarB_energy.push_back(dau->energy());
				MCDstarB_p.push_back(dau->p());
				MCDstarB_et.push_back(dau->et());
				MCDstarB_rapidity.push_back(dau->rapidity());
				MCDstarB_mass.push_back(dau->mass());
 				// D0 Quantities
				MCDstarB_D0eta.push_back(granddau->eta());
				MCDstarB_D0phi.push_back(granddau->phi());
				MCDstarB_D0pt.push_back(granddau->pt());
				MCDstarB_D0energy.push_back(granddau->energy());
				MCDstarB_D0p.push_back(granddau->p());
				MCDstarB_D0et.push_back(granddau->et());
				MCDstarB_D0rapidity.push_back(granddau->rapidity());   
				MCDstarB_D0mass.push_back(granddau->mass());
				// Find D0 displacement and Lifetime
				double D0Lxy = sqrt( pow(granddau->vx()-granddau->daughter(0)->vx(),2) + 
											pow(granddau->vy()-granddau->daughter(0)->vy(),2) );
				double D0lifeTime = granddau->mass() * D0Lxy*pow(10,-2) / (c * granddau->pt() );
				// cout << "MC D0(fromD*) displacement-xy[cm]: " << D0Lxy << " # LifeTime [s]: " << D0lifeTime << endl;
				MCDstarB_D0dispXY.push_back(D0Lxy);   
				MCDstarB_D0lifetime.push_back(D0lifeTime);
				//SlowPion Quantities
				MCDstarB_Seta.push_back(slow->eta());
				MCDstarB_Sphi.push_back(slow->phi());
				MCDstarB_Spt.push_back(slow->pt());
				MCDstarB_Senergy.push_back(slow->energy());
				MCDstarB_Sp.push_back(slow->p());
				MCDstarB_Set.push_back(slow->et());
				MCDstarB_Srapidity.push_back(slow->rapidity());
				MCDstarB_Smass.push_back(slow->mass());					
				//Kaon Quantities
				MCDstarB_Keta.push_back(granddau->daughter(index00)->eta());
				MCDstarB_Kphi.push_back(granddau->daughter(index00)->phi());
				MCDstarB_Kpt.push_back(granddau->daughter(index00)->pt());
				MCDstarB_Kenergy.push_back(granddau->daughter(index00)->energy());
				MCDstarB_Kp.push_back(granddau->daughter(index00)->p());
				MCDstarB_Ket.push_back(granddau->daughter(index00)->et());
				MCDstarB_Krapidity.push_back(granddau->daughter(index00)->rapidity());
				MCDstarB_Kmass.push_back(granddau->daughter(index00)->mass());
				//Pion Quantities
				MCDstarB_Pieta.push_back(granddau->daughter(index11)->eta());
				MCDstarB_Piphi.push_back(granddau->daughter(index11)->phi());
				MCDstarB_Pipt.push_back(granddau->daughter(index11)->pt());
				MCDstarB_Pienergy.push_back(granddau->daughter(index11)->energy());
				MCDstarB_Pip.push_back(granddau->daughter(index11)->p());
				MCDstarB_Piet.push_back(granddau->daughter(index11)->et());
				MCDstarB_Pirapidity.push_back(granddau->daughter(index11)->rapidity());
				MCDstarB_Pimass.push_back(granddau->daughter(index11)->mass());						
			}// B number Of Daughters	
    	}// Gen Particle Loop
   } // doMC Flag
}//End GenDstarMesonBInfo

//***********************************************************************************
void DstarD0TTree::RecD0(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::Vertex& RecVtx){

	if( goodTracksD0.size() >= 2 ) {
	for(size_t i = 0; i < (goodTracksD0.size()-1); i++)
	{
		TransientTrack trk1D0 = goodTracksD0[i];
		for(size_t j=i+1;j<goodTracksD0.size();j++)
		{
			CounterTracksD0combination++;
			TransientTrack trk2D0 = goodTracksD0[j];
			//Testing charge and if tracks are equal
			if(trk1D0 == trk2D0) continue;
			if(trk1D0.charge() == trk2D0.charge()) continue;

			//TransientTrack KfromD0 = 0, PifromD0 = 0;	
			TransientTrack KfromD0, PifromD0;	        
            			
			mass1 = 0.;
			//Reconstructio of 4-vector D0 Prompt for track1 = Kaon and track2 = pion
			math::XYZTLorentzVector vD0_1(trk1D0.track().px(),
											trk1D0.track().py(),
											trk1D0.track().pz(),
											sqrt(pow(trk1D0.track().p(),2)+pow(k_mass,2)) );
			math::XYZTLorentzVector vD0_2(trk2D0.track().px(),
									trk2D0.track().py(),
									trk2D0.track().pz(),
									sqrt(pow(trk2D0.track().p(),2)+pow(pi_mass,2)) );        

			math::XYZTLorentzVector vD0 = vD0_1 + vD0_2;
			mass1 = vD0.M();

			if( fabs(mass1-1.86484) > 0.2) continue;

			// tracks combinations	
			if( vD0_1.M() > vD0_2.M() ) {KfromD0 = trk1D0; PifromD0 = trk2D0;}
			else {KfromD0 = trk2D0; PifromD0 = trk1D0;} 

			CounterD0PromptminusPDG1p0++;

			//KalmanVertexFitter
			vector<TransientTrack> tksD0;
			tksD0.push_back(KfromD0);
			tksD0.push_back(PifromD0);
			KalmanVertexFitter kalman(true);
			TransientVertex v_D0 = kalman.vertex(tksD0);
			if(!v_D0.isValid() || !v_D0.hasRefittedTracks()) continue;

			//Vertex Probability
			double D0KpivtxProb =TMath::Prob( (Double_t) v_D0.totalChiSquared(), (Int_t) v_D0.degreesOfFreedom());

			TransientTrack KfromD0_f = v_D0.refittedTrack(KfromD0);
			TransientTrack pifromD0_f = v_D0.refittedTrack(PifromD0);
			CounterD0PromptKpiAfterTransientp0++;

			//SV Confidence Level
			if(D0KpivtxProb < 0.01) continue;
			CounterD0PromptSVConfidenceLevel++;

			//D* 4-momentum Reconstruction after KalmanVertexFitter		
			math::XYZTLorentzVector p4_KfromD0( KfromD0_f.track().px(),
															KfromD0_f.track().py(),
															KfromD0_f.track().pz(),
															sqrt(pow(KfromD0_f.track().p(),2)+pow(k_mass,2))  );
			math::XYZTLorentzVector p4_pifromD0(pifromD0_f.track().px(),
															pifromD0_f.track().py(),
															pifromD0_f.track().pz(),
															sqrt(pow(pifromD0_f.track().p(),2)+pow(pi_mass,2)) );  
			math::XYZTLorentzVector d0kpi_p4 = p4_KfromD0 + p4_pifromD0;
			double d0kpimass = d0kpi_p4.M();

			//Angle between formation and decay of D0
			double dispAngle = FindAngle(RecVtx,v_D0,d0kpi_p4);
			double D0cosPhi = cos(dispAngle);
			if( D0cosPhi < 0.99 ) continue;
			CounterD0PromptPointingcosPhi++;
			
			//D0 Significance
			VertexDistanceXY vD0KpidXY ;			
			double D0KpidXY = vD0KpidXY.distance(RecVtx,v_D0).value() ;
			double D0KpieXY = vD0KpidXY.distance(RecVtx,v_D0).error() ;
			double D0KpisXY = D0KpidXY / D0KpieXY;

			//D0 significance 3D
			VertexDistance3D vD0Kpid3D ;
			double D0Kpid3D = vD0Kpid3D.distance(RecVtx,v_D0).value() ;
			double D0Kpie3D = vD0Kpid3D.distance(RecVtx,v_D0).error() ;
			double D0Kpis3D = D0Kpid3D / D0Kpie3D;

			if(selectionCuts) {if( D0Kpis3D < D0Significance3D_ ) continue;}
			if(selectionCuts) {CounterD0PromptSignificance3D++;}
			
			// Kaon momentum
			math::XYZVector K_p = trk2D0.track().momentum();
			// D0 momentum
			math::XYZVector D0_p = trk1D0.track().momentum() + trk2D0.track().momentum();			
			//Vetorial product of 4-momentum Kaon and 4-momentum D0
			double D0_kT = sqrt( (K_p).Cross(D0_p).Mag2() / D0_p.Mag2() ) ;
			
			if(selectionCuts) {if(fabs(d0kpimass - 1.86484)>0.15) continue;}
			if(selectionCuts) {CounterD0PromptCandidates++;}
			
			D0Kpi_VtxProb.push_back(D0KpivtxProb);
			D0Kpimass.push_back(d0kpi_p4.M());
			D0Kpipt.push_back(d0kpi_p4.Pt());
			D0Kpieta.push_back(d0kpi_p4.eta());
			D0Kpiphi.push_back(d0kpi_p4.phi());

			D0Kpi_VtxPosx.push_back(v_D0.position().x());
			D0Kpi_VtxPosy.push_back(v_D0.position().y());
			D0Kpi_VtxPosz.push_back(v_D0.position().z());
			D0Kpi_Vtxerrx.push_back(v_D0.positionError().cxx());
			D0Kpi_Vtxerry.push_back(v_D0.positionError().cyy());
			D0Kpi_Vtxerrz.push_back(v_D0.positionError().czz());
			D0Kpi_DispAngle.push_back(D0cosPhi);
			TrkD0Kdxy.push_back(KfromD0_f.track().dxy(RecVtx.position()));
			TrkD0pidxy.push_back(pifromD0_f.track().dxy(RecVtx.position()));

			TrkD0Kdz.push_back(KfromD0_f.track().dz(RecVtx.position()));
			TrkD0Knhits.push_back(KfromD0.track().numberOfValidHits());
			TrkD0Kchi2.push_back(KfromD0.track().normalizedChi2());
			TrkD0Kpt.push_back(KfromD0_f.track().pt());
			TrkD0Keta.push_back(KfromD0_f.track().eta());
			TrkD0Kphi.push_back(KfromD0_f.track().phi());

			TrkD0pidz.push_back(pifromD0_f.track().dz(RecVtx.position()));
			TrkD0pinhits.push_back(PifromD0.track().numberOfValidHits());
			TrkD0pichi2.push_back(PifromD0.track().normalizedChi2());
			TrkD0pipt.push_back(pifromD0_f.track().pt());
			TrkD0pieta.push_back(pifromD0_f.track().eta());
			TrkD0piphi.push_back(pifromD0_f.track().phi());

			D0KpisXY_vec.push_back(D0KpisXY);
			D0KpidXY_vec.push_back(D0KpidXY);
			D0KpieXY_vec.push_back(D0KpieXY);  
			D0Kpis3D_vec.push_back(D0Kpis3D);
			D0Kpid3D_vec.push_back(D0Kpid3D);
			D0Kpie3D_vec.push_back(D0Kpie3D);
			D0_kT_vec.push_back(D0_kT);
					
		}//second track loop
	}//first track loop
} // protection
} //END OF RECD0 

//***********************************************************************************
void DstarD0TTree::GenD0Info(const edm::Event& iEvent,
									  const edm::EventSetup& iSetup){
	if(doMC)
	{
		edm::Handle<GenParticleCollection> gensD0;
		iEvent.getByToken(genParticlesTokenD0_, gensD0);
		int IDdau1, IDdau2;
		int index00, index11;

		for(size_t i=0; i< gensD0->size(); i++)
		{
			const GenParticle & p = (*gensD0)[i];		
			// D0
			if( fabs(p.pdgId()) != 421 ) continue;
			// 2 daughters
			if( p.numberOfDaughters() != 2) continue;
			 
			IDdau1 = p.daughter(0)->pdgId();
			IDdau2 = p.daughter(1)->pdgId();
			// If the first particle is a Kaon
			if( ( fabs(IDdau1) == 321) 
			&&  ( fabs(IDdau2) == 211) 
			&&  ( (IDdau1*IDdau2) < 0) ){ index00 =0; index11 =1;}
			// If the first particle is a Pion
			else if( ( fabs(IDdau1) == 211) 
			&& 		( fabs(IDdau2) == 321) 
			&& 		( (IDdau1*IDdau2) < 0) ){ index00 =1; index11 =0;}
			else continue;
			// D0 Quantites		
			MCpromptD0eta.push_back(p.eta());
			MCpromptD0phi.push_back(p.phi());
			MCpromptD0pt.push_back(p.pt());
			MCpromptD0energy.push_back(p.energy());
			MCpromptD0p.push_back(p.p());
			MCpromptD0et.push_back(p.et());
			MCpromptD0rapidity.push_back(p.rapidity());
			MCpromptD0mass.push_back(p.mass());
			//Find angle between Primary and Secondary Vertex
			CLHEP::Hep3Vector MCdisplacement(p.mother(0)->vx() - p.vx(), 
														p.mother(0)->vy() - p.vy(), 
														p.mother(0)->vz() - p.vz() );
			CLHEP::Hep3Vector momentum( p.px(), p.py(), p.pz() );
			double MCD0dispAngle = momentum.angle(MCdisplacement); 						
			MCpromptD0_DispAngle.push_back(MCD0dispAngle);
			//Find D0 displacement and Lifetime
			double D0Lxy = sqrt( pow(p.vx()-p.daughter(0)->vx(),2) + 
									 	pow(p.vy()-p.daughter(0)->vy(),2) );
			double D0lifeTime = p.mass() * D0Lxy*pow(10,-2) / (c * p.pt() );
			//cout << "MC D0 displacement[cm]: " << D0Lxy << " # LifeTime [s]: " << D0lifeTime << endl;
			MCpromptD0dispXY.push_back(D0Lxy);   
			MCpromptD0lifetime.push_back(D0lifeTime);				
			//Kaon Quantities
			MCpromptD0_Keta.push_back(p.daughter(index00)->eta());
			MCpromptD0_Kphi.push_back(p.daughter(index00)->phi());
			MCpromptD0_Kpt.push_back(p.daughter(index00)->pt());
			MCpromptD0_Kenergy.push_back(p.daughter(index00)->energy());
			MCpromptD0_Kp.push_back(p.daughter(index00)->p());
			MCpromptD0_Ket.push_back(p.daughter(index00)->et());
			MCpromptD0_Krapidity.push_back(p.daughter(index00)->rapidity());
			MCpromptD0_Kmass.push_back(p.daughter(index00)->mass());
			//Pion Quantities
			MCpromptD0_Pieta.push_back(p.daughter(index11)->eta());
			MCpromptD0_Piphi.push_back(p.daughter(index11)->phi());
			MCpromptD0_Pipt.push_back(p.daughter(index11)->pt());
			MCpromptD0_Pienergy.push_back(p.daughter(index11)->energy());
			MCpromptD0_Pip.push_back(p.daughter(index11)->p());
			MCpromptD0_Piet.push_back(p.daughter(index11)->et());
			MCpromptD0_Pirapidity.push_back(p.daughter(index11)->rapidity());
			MCpromptD0_Pimass.push_back(p.daughter(index11)->mass());									
			//If D0 has a boson W
			//if(fabs(IDdau1) == 24 or fabs(IDdau2) == 24 )
			//{ cout << "-------IDdau1 = W: "<< IDdau1 << " # IDdau2 = W: "<< IDdau2 << endl;
			//}			
		}// loop Gen Particle
	}// Flag MC
}//End GenD0Info
//***********************************************************************************
void DstarD0TTree::GenD0MesonBInfo(const edm::Event& iEvent,
										  const edm::EventSetup& iSetup){
	if(doMC) //MC flag - see config file
	{ 
		//Handle<GenParticleCollection> genParticles; 
		edm::Handle<GenParticleCollection> gens; //linked to "prunedGenParticles" in config file
		iEvent.getByToken(genParticlesTokenDstar_, gens);
		int IDdau1, IDdau2;
		int index00, index11;

		for(size_t i=0;i<gens->size();i++)
		{		
     		const GenParticle & p = (*gens)[i];
			// B0 or B+-
      	if( (fabs(p.pdgId()) != 511) && (fabs(p.pdgId()) != 521) ) continue; 			
     		for(size_t j=0; j < p.numberOfDaughters(); j++)
			{	
          	const Candidate* dau = p.daughter(j);
				// D0
				if( fabs(dau->pdgId()) != 421 ) continue;
				// 2 daughters
				if( dau->numberOfDaughters()!=2 ) continue;
    						
				IDdau1 = dau->daughter(0)->pdgId(); 
				IDdau2 = dau->daughter(1)->pdgId();
				// If the first particle is a Kaon
				if( ( fabs(IDdau1) == 321) 
				&&  ( fabs(IDdau2) == 211) 
				&&  ( (IDdau1*IDdau2) < 0) ){ index00 =0; index11 =1;}
				// If the first particle is a Pion
				else if( ( fabs(IDdau1) == 211) 
				&& 		( fabs(IDdau2) == 321) 
				&& 		( (IDdau1*IDdau2) < 0) ){ index00 =1; index11 =0;}
				else continue;									
 				// D0 Quantities
				MCD0B_D0eta.push_back(dau->eta());
				MCD0B_D0phi.push_back(dau->phi());
				MCD0B_D0pt.push_back(dau->pt());
				MCD0B_D0energy.push_back(dau->energy());
				MCD0B_D0p.push_back(dau->p());
				MCD0B_D0et.push_back(dau->et());
				MCD0B_D0rapidity.push_back(dau->rapidity());   
				MCD0B_D0mass.push_back(dau->mass());
				//Find angle between Primary and Secondary Vertex
				CLHEP::Hep3Vector MCdisplacement(p.mother(0)->vx() - p.vx(), 
															p.mother(0)->vy() - p.vy(), 
															p.mother(0)->vz() - p.vz() );
				CLHEP::Hep3Vector momentum( p.px(), p.py(), p.pz() );
				double MCD0dispAngle = momentum.angle(MCdisplacement); 						
				MCD0B_D0DispAngle.push_back(MCD0dispAngle);
				// Find D0 displacement and Lifetime
				double D0Lxy = sqrt( pow(dau->vx()-dau->daughter(0)->vx(),2) + 
											pow(dau->vy()-dau->daughter(0)->vy(),2) );
				double D0lifeTime = dau->mass() * D0Lxy*pow(10,-2) / (c * dau->pt() );
				// cout << "MC D0(fromD*) displacement-xy[cm]: " << D0Lxy << " # LifeTime [s]: " << D0lifeTime << endl;
				MCD0B_D0dispXY.push_back(D0Lxy);   
				MCD0B_D0lifetime.push_back(D0lifeTime);
				//Kaon Quantities
				MCD0B_Keta.push_back(dau->daughter(index00)->eta());
				MCD0B_Kphi.push_back(dau->daughter(index00)->phi());
				MCD0B_Kpt.push_back(dau->daughter(index00)->pt());
				MCD0B_Kenergy.push_back(dau->daughter(index00)->energy());
				MCD0B_Kp.push_back(dau->daughter(index00)->p());
				MCD0B_Ket.push_back(dau->daughter(index00)->et());
				MCD0B_Krapidity.push_back(dau->daughter(index00)->rapidity());
				MCD0B_Kmass.push_back(dau->daughter(index00)->mass());
				//Pion Quantities
				MCD0B_Pieta.push_back(dau->daughter(index11)->eta());
				MCD0B_Piphi.push_back(dau->daughter(index11)->phi());
				MCD0B_Pipt.push_back(dau->daughter(index11)->pt());
				MCD0B_Pienergy.push_back(dau->daughter(index11)->energy());
				MCD0B_Pip.push_back(dau->daughter(index11)->p());
				MCD0B_Piet.push_back(dau->daughter(index11)->et());
				MCD0B_Pirapidity.push_back(dau->daughter(index11)->rapidity());
				MCD0B_Pimass.push_back(dau->daughter(index11)->mass());														
			}// B number Of Daughters	
    	}// Gen Particle Loop
   } // doMC Flag
}//End GenD0MesonBInfo

//***********************************************************************************
//Calculate the angle between two vectors
double DstarD0TTree::FindAngle(const reco::Vertex& pv , 
										 const TransientVertex& sv, const 
										 math::XYZTLorentzVector& vD0angle)
{
	CLHEP::Hep3Vector displacement(sv.position().x() - pv.position().x(),
											 sv.position().y() - pv.position().y(),
											 sv.position().z() - pv.position().z() );
	CLHEP::Hep3Vector momentum( vD0angle.Px(), 
										 vD0angle.Py(), 
										 vD0angle.Pz() ); 

	return momentum.angle(displacement) ;
}
//***********************************************************************************
void DstarD0TTree::initialize(){

	goodTracks.clear(); 
	goodTracksD0.clear(); slowPiTracks.clear();
	NKpiCand=0; CounterDsCandidates=0; 
	PVx = PVy = PVz = PVerrx = PVerry = PVerrz = -999.;
	ntracksD0Kpi = -999; ntracksDstar = -999; n_pVertex = -999; 
	runNumber=0; eventNumber=0; lumi=0; lumiWeight_=0.;

	NameOfFiredTriggers.clear();

	//counters
	Total_Events = 0; CounterTotalTracks = 0; CounterTracksHasTrackDetails =0; CounterTracksCharged = 0; CounterTracksEta = 0; CounterTracksHighPurity = 0; CounterTracksPDG211 = 0; CounterTracksPtZeroFive = 0; CounterTracksChi3 = 0; CounterTracksNumberOfHits2 = 0; CounterTracksDxyThree = 0; CounterTracksDzThree = 0; CounterTracksPtZeroSix = 0; CounterTracksChiTwoFive = 0; CounterTracksNumberOfHits5 = 0; CounterTracksNumberOfPixelHits2 = 0; CounterTracksDxyZeroOne = 0; CounterTracksDzOne = 0;

	//D0 Conters
	CounterTracksPt0p8 = 0; CounterD0PromptTracksEta2p5 = 0; CounterD0PromptTracksChi5 =0; CounterD0PromptTracksNumberHits5 = 0;
	CounterD0PromptTracksPixelHits2 = 0; CounterD0PromptTracksDxy0p1 = 0; CounterD0PromptTracksDz0p5 =0; CounterTracksD0combination = 0; CounterD0PromptminusPDG1p0 = 0; CounterD0PromptKpiAfterTransientp0 =0; CounterD0PromptSVConfidenceLevel = 0; CounterD0PromptPointingcosPhi = 0; CounterD0PromptSignificance3D =0; CounterD0PromptCandidates = 0;
	
	CounterD0AfterLorentzVector = 0; CounterD0MinusPDGzero2 = 0; CounterDsMinusD0Zerothree = 0; CounterTransientTrackOfpiK = 0; CounterSVConfidenceLevel = 0; CounterD0AfterLorentzVectorKarman = 0; CounterPointingcosPhi = 0;	CounterSignificance3D = 0; CounterD0pTThree = 0; CounterDsAfterLorentzVector = 0; CounterD0MinusPDG = 0; CounterDsMinusD0 = 0;

	D0_VtxProb.clear(); D0mass.clear(); TrkKmass.clear(); Trkpimass.clear(); TrkSmass.clear(); Dsmass.clear(); D0pt.clear(); Dspt.clear(); D0eta.clear(); D0phi.clear(); Dseta.clear(); Dsphi.clear(); D0_VtxPosx.clear(); D0_VtxPosy.clear(); D0_VtxPosz.clear(); D0_Vtxerrx.clear(); D0_Vtxerry.clear(); D0_Vtxerrz.clear(); TrkKdxy.clear(); Trkpidxy.clear(); TrkSdxy.clear(); TrkKdz.clear(); Trkpidz.clear();TrkSdz.clear(); TrkKnhits.clear(); Trkpinhits.clear(); TrkSnhits.clear();	TrkKchi2.clear(); Trkpichi2.clear(); TrkSchi2.clear(); DSDeltaR.clear(); TrkKpt.clear(); Trkpipt.clear(); TrkSpt.clear(); TrkKeta.clear(); Trkpieta.clear(); TrkSeta.clear(); TrkKphi.clear(); Trkpiphi.clear(); TrkSphi.clear(); TrkScharge.clear(); D0fromDSsXY_vec.clear(); D0fromDSs3D_vec.clear(); Anglephi_vec.clear(); D0fromDSd3D_vec.clear(); D0fromDSe3D_vec.clear(); D0fromDSdXY_vec.clear(); D0fromDSeXY_vec.clear();
	//----------------------------
	//D* QUANTITIES WRONG COMBINATION
	//----------------------------
	D0_VtxProbWrong.clear(); D0massWrong.clear(); TrkKmassWrong.clear(); TrkpimassWrong.clear();  TrkSmassWrong.clear(); DsmassWrong.clear(); D0ptWrong.clear(); DsptWrong.clear();D0etaWrong.clear(); D0phiWrong.clear(); DsetaWrong.clear(); DsphiWrong.clear(); D0_VtxPosxWrong.clear(); D0_VtxPosyWrong.clear(); D0_VtxPoszWrong.clear(); D0_VtxerrxWrong.clear(); D0_VtxerryWrong.clear(); D0_VtxerrzWrong.clear(); TrkKdxyWrong.clear(); TrkpidxyWrong.clear(); TrkSdxyWrong.clear(); TrkKdzWrong.clear();TrkpidzWrong.clear();TrkSdzWrong.clear(); TrkKnhitsWrong.clear(); TrkpinhitsWrong.clear(); TrkSnhitsWrong.clear(); TrkKchi2Wrong.clear(); Trkpichi2Wrong.clear(); TrkSchi2Wrong.clear(); DSDeltaRWrong.clear(); TrkKptWrong.clear(); TrkpiptWrong.clear(); TrkSptWrong.clear(); TrkKetaWrong.clear(); TrkpietaWrong.clear(); TrkSetaWrong.clear(); TrkKphiWrong.clear();TrkpiphiWrong.clear(); TrkSphiWrong.clear();TrkSchargeWrong.clear(); D0fromDSsXY_vecWrong.clear(); D0fromDSs3D_vecWrong.clear(); Anglephi_vecWrong.clear();

	D0Kpi_VtxProb.clear();D0Kpimass.clear();D0Kpipt.clear();D0Kpieta.clear(); D0Kpiphi.clear(); D0Kpi_VtxPosx.clear(); D0Kpi_VtxPosy.clear(); D0Kpi_VtxPosz.clear();D0Kpi_Vtxerrx.clear(); D0Kpi_Vtxerry.clear(); D0Kpi_Vtxerrz.clear(); D0Kpi_DispAngle.clear(); TrkD0pidxy.clear(); TrkD0Kdxy.clear(); TrkD0Kdz.clear(); TrkD0pidz.clear(); TrkD0Knhits.clear(); TrkD0pinhits.clear(); TrkD0Kchi2.clear();TrkD0pichi2.clear(); TrkD0Kpt.clear(); TrkD0pipt.clear(); TrkD0Keta.clear(); TrkD0pieta.clear(); TrkD0Kphi.clear(); TrkD0piphi.clear(); D0KpisXY_vec.clear(); D0KpidXY_vec.clear(); D0KpieXY_vec.clear(); D0Kpis3D_vec.clear(); D0Kpid3D_vec.clear(); D0Kpie3D_vec.clear(); D0_kT_vec.clear();	
	//	D0KpisXY_.clear();
	//----------------------------
	//D* MC
	//----------------------------	
	MCDseta.clear(); MCDsphi.clear(); MCDspt.clear(); MCDsenergy.clear(); MCDsp.clear(); MCDset.clear(); MCDsrapidity.clear(); MCDsmass.clear(); MCD0eta.clear(); MCD0phi.clear(); MCD0pt.clear(); MCD0energy.clear(); MCD0p.clear(); MCD0et.clear(); MCD0rapidity.clear(); MCD0mass.clear(); MCD0dispXY.clear(); MCD0lifetime.clear(); MCDsKeta.clear(); MCDsKphi.clear(); MCDsKpt.clear(); MCDsKenergy.clear(); MCDsKp.clear(); MCDsKet.clear(); MCDsKrapidity.clear(); MCDsKmass.clear(); MCDsPieta.clear(); MCDsPiphi.clear(); MCDsPipt.clear(); MCDsPienergy.clear(); MCDsPip.clear(); MCDsPiet.clear(); MCDsPirapidity.clear(); MCDsPimass.clear(); MCDsSeta.clear(); MCDsSphi.clear(); MCDsSpt.clear(); MCDsSenergy.clear(); MCDsSp.clear(); MCDsSet.clear(); MCDsSrapidity.clear(); MCDsSmass.clear(); 
	//----------------------------
	//D* from B MC
	//----------------------------
	MCDstarB_eta.clear(); MCDstarB_phi.clear(); MCDstarB_pt.clear(); MCDstarB_energy.clear(); MCDstarB_p.clear(); MCDstarB_et.clear(); MCDstarB_rapidity.clear(); MCDstarB_mass.clear(); MCDstarB_D0eta.clear(); MCDstarB_D0phi.clear(); MCDstarB_D0pt.clear(); MCDstarB_D0energy.clear(); MCDstarB_D0p.clear(); MCDstarB_D0et.clear(); MCDstarB_D0rapidity.clear(); MCDstarB_D0mass.clear(); MCDstarB_D0dispXY.clear(); MCDstarB_D0lifetime.clear(); MCDstarB_Keta.clear(); MCDstarB_Kphi.clear(); MCDstarB_Kpt.clear(); MCDstarB_Kenergy.clear(); MCDstarB_Kp.clear(); MCDstarB_Ket.clear(); MCDstarB_Krapidity.clear(); MCDstarB_Kmass.clear(); MCDstarB_Pieta.clear(); MCDstarB_Piphi.clear(); MCDstarB_Pipt.clear(); MCDstarB_Pienergy.clear(); MCDstarB_Pip.clear(); MCDstarB_Piet.clear(); MCDstarB_Pirapidity.clear(); MCDstarB_Pimass.clear(); MCDstarB_Seta.clear(); MCDstarB_Sphi.clear(); MCDstarB_Spt.clear(); MCDstarB_Senergy.clear(); MCDstarB_Sp.clear(); MCDstarB_Set.clear(); MCDstarB_Srapidity.clear(); MCDstarB_Smass.clear(); 
	//----------------------------
	//D0 MC
	//----------------------------	
	MCpromptD0eta.clear(); MCpromptD0phi.clear(); MCpromptD0pt.clear(); MCpromptD0energy.clear(); MCpromptD0p.clear(); MCpromptD0et.clear(); MCpromptD0rapidity.clear(); MCpromptD0mass.clear(); MCpromptD0dispXY.clear(); MCpromptD0lifetime.clear(); MCpromptD0_Keta.clear(); MCpromptD0_Kphi.clear(); MCpromptD0_Kpt.clear(); MCpromptD0_Kenergy.clear(); MCpromptD0_Kp.clear(); MCpromptD0_Ket.clear(); MCpromptD0_Krapidity.clear(); MCpromptD0_Kmass.clear();MCpromptD0_Pieta.clear(); MCpromptD0_Piphi.clear(); MCpromptD0_Pipt.clear(); MCpromptD0_Pienergy.clear(); MCpromptD0_Pip.clear(); MCpromptD0_Piet.clear(); MCpromptD0_Pirapidity.clear(); MCpromptD0_Pimass.clear();MCpromptD0_DispAngle.clear(); MCpromptD0_Kt.clear();
	//----------------------------
	//D0 from B MC
	//----------------------------
	MCD0B_D0eta.clear();	MCD0B_D0phi.clear(); MCD0B_D0pt.clear(); MCD0B_D0energy.clear(); MCD0B_D0p.clear(); MCD0B_D0et.clear(); MCD0B_D0rapidity.clear(); MCD0B_D0mass.clear(); MCD0B_D0dispXY.clear(); MCD0B_D0lifetime.clear(); MCD0B_Keta.clear(); MCD0B_Kphi.clear(); MCD0B_Kpt.clear(); MCD0B_Kenergy.clear(); MCD0B_Kp.clear(); MCD0B_Ket.clear(); MCD0B_Krapidity.clear(); MCD0B_Kmass.clear(); MCD0B_Pieta.clear(); MCD0B_Piphi.clear(); MCD0B_Pipt.clear(); MCD0B_Pienergy.clear(); MCD0B_Pip.clear(); MCD0B_Piet.clear(); MCD0B_Pirapidity.clear(); MCD0B_Pimass.clear(); MCD0B_D0DispAngle.clear();;
	
	TracksCharge.clear(); TracksEta.clear(); TracksPt.clear(); TracksChi2.clear(); TracksNumberOfHits.clear(); TracksdxyError.clear(); TracksdzError.clear(); TracksNumberOfPixelHits.clear(); Tracksdxy.clear(); Tracksdz.clear(); TracksPhi.clear();

}
//++++++++++++++++++
void DstarD0TTree::endJob(){
	cout <<"######################################################################"<<endl;
	cout << "Number of Events: " << eventNumber << " Run Number: " << runNumber << endl;
	cout << "Total events: " << Total_Events_test << " #events triggered by " << triggerName_ 
		<< ": " << Triggered_Event_test << " #events triggered by " << triggerReferenceName_ 
		<< ": " << TriggeredReference_Event_test << " #events triggered by " << triggerReferenceName2_ << ": "
		<< TriggeredReference_Event_test2 << endl;
	cout << "Number of B0's: " << B0test << endl;   
}

//***********************************************************************************
void DstarD0TTree::beginJob(){
		
	//CREATING BRANCHES
	if(doMC)
	{
	//Weighting PileUp 
	const edm::InputTag& PileupSumInfoInputTags = edm::InputTag( "slimmedAddPileupInfo" ) ;
	LumiWeights_ = new edm::LumiReWeighting("PileupMC.root", "MyDataPileupHistogram.root", "input_Event/N_TrueInteractions", "pileup", PileupSumInfoInputTags);
	data->Branch("PUWeight",&PUWeight,"PUWeight/D");
	data->Branch("nPU",&nPU,"nPU/D");
	data->Branch("GenWeight",&GenWeight,"GenWeight/D");
	}
       	
	data->Branch("Total_Events",&Total_Events,"Total_Events/I"); 
	data->Branch("Triggered_Event",&Triggered_Event,"Triggered_Event/I");

	if(TracksOn)
	{
	data->Branch("TracksCharge",&TracksCharge);
	data->Branch("TracksEta",&TracksEta);
	data->Branch("TracksPt",&TracksPt);
	data->Branch("TracksPhi",&TracksPhi);
	data->Branch("TracksChi2",&TracksChi2);
	data->Branch("TracksNumberOfHits",&TracksNumberOfHits);
	data->Branch("TracksNumberOfPixelHits",&TracksNumberOfPixelHits);
	data->Branch("Tracksdxy",&Tracksdxy);
	data->Branch("Tracksdz",&Tracksdz);
	data->Branch("TracksdxyError",&TracksdxyError);
	data->Branch("TracksdzError",&TracksdzError);
	}

	data->Branch("runNumber",&runNumber,"runNumber/I");
	data->Branch("eventNumber",&eventNumber,"eventNumber/I");
	data->Branch("lumi",&lumi,"lumi/I");
	data->Branch("lumiWeight_",&lumiWeight_,"lumiWeight_/D");
	data->Branch("NameOfFiredTriggers",&NameOfFiredTriggers);
	data->Branch("n_pVertex",&n_pVertex,"n_pVertex/I");
	data->Branch("PVx",&PVx,"PVx/D");
	data->Branch("PVy",&PVy,"PVy/D");
	data->Branch("PVz",&PVz,"PVz/D");
	data->Branch("PVerrx",&PVerrx,"PVerrx/D");
	data->Branch("PVerry",&PVerry,"PVerry/D");
	data->Branch("PVerrz",&PVerrz,"PVerrz/D");

	//Counters
	data->Branch("CounterTotalTracks",&CounterTotalTracks,"CounterTotalTracks/L");
 	data->Branch("CounterTracksHasTrackDetails",&CounterTracksHasTrackDetails,"CounterTracksHasTrackDetails/L");
	data->Branch("CounterTracksHighPurity",&CounterTracksHighPurity,"CounterTracksHighPurity/L");
	data->Branch("CounterTracksCharged",&CounterTracksCharged,"CounterTracksCharged/L");
	data->Branch("CounterTracksEta",&CounterTracksEta,"CounterTracksEta/L");

	data->Branch("CounterTracksPtZeroFive",&CounterTracksPtZeroFive,"CounterTracksPtZeroFive/L");
 	data->Branch("CounterTracksChi3",&CounterTracksChi3,"CounterTracksChi3/L");
	data->Branch("CounterTracksNumberOfHits2",&CounterTracksNumberOfHits2,"CounterTracksNumberOfHits2/L");
	data->Branch("CounterTracksDxyThree",&CounterTracksDxyThree,"CounterTracksDxyThree/L");
 	data->Branch("CounterTracksDzThree",&CounterTracksDzThree,"CounterTracksDzThree/L");

	data->Branch("CounterTracksPtZeroSix",&CounterTracksPtZeroSix,"CounterTracksPtZeroSix/L");
	data->Branch("CounterTracksChiTwoFive",&CounterTracksChiTwoFive,"CounterTracksChiTwoFive/L");
	data->Branch("CounterTracksNumberOfHits5",&CounterTracksNumberOfHits5,"CounterTracksNumberOfHits5/L");
	data->Branch("CounterTracksNumberOfPixelHits2",&CounterTracksNumberOfPixelHits2,"CounterTracksNumberOfPixelHits2/L");
	data->Branch("CounterTracksDxyZeroOne",&CounterTracksDxyZeroOne,"CounterTracksDxyZeroOne/L");
	data->Branch("CounterTracksDzOne",&CounterTracksDzOne,"CounterTracksDzOne/L");
	
	//======================================================
	// D* -> D0 + pi_slow  Variables
	//======================================================
	//Counters
	data->Branch("CounterD0AfterLorentzVector",&CounterD0AfterLorentzVector,"CounterD0AfterLorentzVector/L");
	data->Branch("CounterD0MinusPDGzero2",&CounterD0MinusPDGzero2,"CounterD0MinusPDGzero2/L");
	data->Branch("CounterDsMinusD0Zerothree",&CounterDsMinusD0Zerothree,"CounterDsMinusD0Zerothree/L");
	data->Branch("CounterTransientTrackOfpiK",&CounterTransientTrackOfpiK,"CounterTransientTrackOfpiK/L");
	data->Branch("CounterSVConfidenceLevel",&CounterSVConfidenceLevel,"CounterSVConfidenceLevel/L");
	data->Branch("CounterD0AfterLorentzVectorKarman",&CounterD0AfterLorentzVectorKarman,"CounterD0AfterLorentzVectorKarman/L");
	data->Branch("CounterPointingcosPhi",&CounterPointingcosPhi,"CounterPointingcosPhi/L");
	data->Branch("CounterSignificance3D",&CounterSignificance3D,"CounterSignificance3D/L");
	data->Branch("CounterD0pTThree",&CounterD0pTThree,"CounterD0pTThree/L");
	data->Branch("CounterDsAfterLorentzVector",&CounterDsAfterLorentzVector,"CounterDsAfterLorentzVector/L");
	data->Branch("CounterD0MinusPDG",&CounterD0MinusPDG,"CounterD0MinusPDG/L");
	data->Branch("CounterDsMinusD0",&CounterDsMinusD0,"CounterDsMinusD0/L");
	data->Branch("CounterDsCandidates",&CounterDsCandidates,"CounterDsCandidates/L");

	data->Branch("Dsmass",&Dsmass);
	data->Branch("Dseta",&Dseta);
	data->Branch("Dspt",&Dspt);
	data->Branch("Dsphi",&Dsphi);
	data->Branch("DSDeltaR",&DSDeltaR);

	data->Branch("D0mass",&D0mass);
	data->Branch("D0eta",&D0eta);
	data->Branch("D0phi",&D0phi);
	data->Branch("D0pt",&D0pt);

	data->Branch("D0_VtxProb",&D0_VtxProb);
	data->Branch("D0_VtxPosx",&D0_VtxPosx);
	data->Branch("D0_VtxPosy",&D0_VtxPosy);
	data->Branch("D0_VtxPosz",&D0_VtxPosz);
	data->Branch("D0_Vtxerrx",&D0_Vtxerrx);
	data->Branch("D0_Vtxerry",&D0_Vtxerry);
	data->Branch("D0_Vtxerrz",&D0_Vtxerrz);

	data->Branch("Anglephi",&Anglephi_vec);
	data->Branch("D0fromDSsXY",&D0fromDSsXY_vec);
	data->Branch("D0fromDSdXY",&D0fromDSdXY_vec);
	data->Branch("D0fromDSeXY",&D0fromDSeXY_vec);
	data->Branch("D0fromDSs3D",&D0fromDSs3D_vec);
	data->Branch("D0fromDSd3D",&D0fromDSd3D_vec);
	data->Branch("D0fromDSe3D",&D0fromDSe3D_vec);

	data->Branch("TrkKmass",&TrkKmass);
	data->Branch("TrkKpt",&TrkKpt);
	data->Branch("TrkKnhits",&TrkKnhits);
	data->Branch("TrkKchi2",&TrkKchi2);
	data->Branch("TrkKdxy",&TrkKdxy);
	data->Branch("TrkKdz",&TrkKdz);
	data->Branch("TrkKeta",&TrkKeta);      
	data->Branch("TrkKphi",&TrkKphi);

	data->Branch("Trkpimass",&Trkpimass);
	data->Branch("Trkpipt",&Trkpipt);
	data->Branch("Trkpinhits",&Trkpinhits);
	data->Branch("Trkpichi2",&Trkpichi2);
	data->Branch("Trkpidxy",&Trkpidxy);
	data->Branch("Trkpidz",&Trkpidz);
	data->Branch("Trkpieta",&Trkpieta);
	data->Branch("Trkpiphi",&Trkpiphi);

	data->Branch("TrkSmass",&TrkSmass);
	data->Branch("TrkSpt",&TrkSpt);
	data->Branch("TrkSnhits",&TrkSnhits);
	data->Branch("TrkSchi2",&TrkSchi2);
	data->Branch("TrkSdxy",&TrkSdxy);
	data->Branch("TrkSdz",&TrkSdz); 
	data->Branch("TrkSeta",&TrkSeta);
	data->Branch("TrkSphi",&TrkSphi);
	data->Branch("TrkScharge",&TrkScharge);

	//======================================================
	// D* -> D0 + pi_slow  WRONG COMBINATION
	//======================================================
	data->Branch("DsmassWrong",&DsmassWrong);
	data->Branch("DsetaWrong",&DsetaWrong);
	data->Branch("DsphiWrong",&DsphiWrong);
	data->Branch("DsptWrong",&DsptWrong);

	data->Branch("D0massWrong",&D0massWrong);
	data->Branch("D0etaWrong",&D0etaWrong);
	data->Branch("D0phiWrong",&D0phiWrong);
	data->Branch("D0ptWrong",&D0ptWrong);

	data->Branch("AnglephiWrong",&Anglephi_vecWrong);
	data->Branch("D0fromDSsXYWrong",&D0fromDSsXY_vecWrong);
	data->Branch("D0fromDSs3DWrong",&D0fromDSs3D_vecWrong);
	data->Branch("D0_VtxProbWrong",&D0_VtxProbWrong);
	/*data->Branch("DSDeltaRWrong",&DSDeltaRWrong);

	data->Branch("D0_VtxPosxWrong",&D0_VtxPosxWrong);
	data->Branch("D0_VtxPosyWrong",&D0_VtxPosyWrong);
	data->Branch("D0_VtxPoszWrong",&D0_VtxPoszWrong);
	data->Branch("D0_VtxerrxWrong",&D0_VtxerrxWrong);
	data->Branch("D0_VtxerryWrong",&D0_VtxerryWrong);
	data->Branch("D0_VtxerrzWrong",&D0_VtxerrzWrong);

	data->Branch("TrkKptWrong",&TrkKptWrong);
	data->Branch("TrkKnhitsWrong",&TrkKnhitsWrong);
	data->Branch("TrkKchi2Wrong",&TrkKchi2Wrong);
	data->Branch("TrkKdxyWrong",&TrkKdxyWrong);
	data->Branch("TrkKdzWrong",&TrkKdzWrong);
	data->Branch("TrkKetaWrong",&TrkKetaWrong);      
	data->Branch("TrkKphiWrong",&TrkKphiWrong);
	data->Branch("TrkKmassWrong",&TrkKmassWrong);

	data->Branch("TrkpimassWrong",&TrkpimassWrong);
	data->Branch("TrkpiptWrong",&TrkpiptWrong);
	data->Branch("TrkpinhitsWrong",&TrkpinhitsWrong);
	data->Branch("Trkpichi2Wrong",&Trkpichi2Wrong);
	data->Branch("TrkpidxyWrong",&TrkpidxyWrong);
	data->Branch("TrkpidzWrong",&TrkpidzWrong);
	data->Branch("TrkpietaWrong",&TrkpietaWrong);
	data->Branch("TrkpiphiWrong",&TrkpiphiWrong);

	data->Branch("TrkSmassWrong",&TrkSmassWrong);
	data->Branch("TrkSptWrong",&TrkSptWrong);
	data->Branch("TrkSnhitsWrong",&TrkSnhitsWrong);
	data->Branch("TrkSchi2Wrong",&TrkSchi2Wrong);
	data->Branch("TrkSdxyWrong",&TrkSdxyWrong);
	data->Branch("TrkSdzWrong",&TrkSdzWrong); 
	data->Branch("TrkSetaWrong",&TrkSetaWrong);
	data->Branch("TrkSphiWrong",&TrkSphiWrong);
	data->Branch("TrkSchargeWrong",&TrkSchargeWrong);*/

    //======================================================
	// D0 Prompt Variables
	//======================================================
	//Counters
	data->Branch("CounterTracksPt0p8",&CounterTracksPt0p8,"CounterTracksPt0p8/L");
	data->Branch("CounterD0PromptTracksEta2p5",&CounterD0PromptTracksEta2p5,"CounterD0PromptTracksEta2p5/L");
	data->Branch("CounterD0PromptTracksChi5",&CounterD0PromptTracksChi5,"CounterD0PromptTracksChi5/L");
	data->Branch("CounterD0PromptTracksNumberHits5",&CounterD0PromptTracksNumberHits5,"CounterD0PromptTracksNumberHits5/L");
	data->Branch("CounterD0PromptTracksPixelHits2",&CounterD0PromptTracksPixelHits2,"CounterD0PromptTracksPixelHits2/L");
	data->Branch("CounterD0PromptTracksDxy0p1",&CounterD0PromptTracksDxy0p1,"CounterD0PromptTracksDxy0p1/L");
	data->Branch("CounterD0PromptTracksDz0p5",&CounterD0PromptTracksDz0p5,"CounterD0PromptTracksDz0p5/L");
	data->Branch("CounterTracksD0combination",&CounterTracksD0combination,"CounterTracksD0combination/L");
	data->Branch("CounterD0PromptminusPDG1p0",&CounterD0PromptminusPDG1p0,"CounterD0PromptminusPDG1p0/L");
	data->Branch("CounterD0PromptKpiAfterTransientp0",&CounterD0PromptKpiAfterTransientp0,"CounterD0PromptKpiAfterTransientp0/L");
	data->Branch("CounterD0PromptSVConfidenceLevel",&CounterD0PromptSVConfidenceLevel,"CounterD0PromptSVConfidenceLevel/L");
	data->Branch("CounterD0PromptPointingcosPhi",&CounterD0PromptPointingcosPhi,"CounterD0PromptPointingcosPhi/L");
	data->Branch("CounterD0PromptSignificance3D",&CounterD0PromptSignificance3D,"CounterD0PromptSignificance3D/L");
	data->Branch("CounterD0PromptCandidates",&CounterD0PromptCandidates,"CounterD0PromptCandidates/L");

	//Variables data
	data->Branch("D0Kpimass",&D0Kpimass);
	data->Branch("D0Kpipt",&D0Kpipt);
	data->Branch("D0Kpieta",&D0Kpieta);
	data->Branch("D0Kpiphi",&D0Kpiphi);

	data->Branch("D0Kpi_VtxProb",&D0Kpi_VtxProb);
	data->Branch("D0KpiDispAngle",&D0Kpi_DispAngle);
	data->Branch("D0KpisXY",&D0KpisXY_vec);
	data->Branch("D0KpidXY",&D0KpidXY_vec);
	data->Branch("D0KpieXY",&D0KpieXY_vec);
	data->Branch("D0Kpis3D",&D0Kpis3D_vec);
	data->Branch("D0Kpid3D",&D0Kpid3D_vec);
	data->Branch("D0Kpie3D",&D0Kpie3D_vec);
	data->Branch("D0KpikT",&D0_kT_vec);

	data->Branch("D0Kpi_VtxPosx",&D0Kpi_VtxPosx);
	data->Branch("D0Kpi_VtxPosy",&D0Kpi_VtxPosy);
	data->Branch("D0Kpi_VtxPosz",&D0Kpi_VtxPosz);
	data->Branch("D0Kpi_Vtxerrx",&D0Kpi_Vtxerrx);
	data->Branch("D0Kpi_Vtxerry",&D0Kpi_Vtxerry);
	data->Branch("D0Kpi_Vtxerrz",&D0Kpi_Vtxerrz);

	data->Branch("TrkD0Kdxy",&TrkD0Kdxy);
	data->Branch("TrkD0Kdz",&TrkD0Kdz);
	data->Branch("TrkD0Kchi2",&TrkD0Kchi2);
	data->Branch("TrkD0Kpt",&TrkD0Kpt);
	data->Branch("TrkD0Keta",&TrkD0Keta);
	data->Branch("TrkD0kphi",&TrkD0Kphi);
	data->Branch("TrkD0Knhits",&TrkD0Knhits);

	data->Branch("TrkD0pidxy",&TrkD0pidxy);
	data->Branch("TrkD0pidz",&TrkD0pidz);
	data->Branch("TrkD0pichi2",&TrkD0pichi2);
	data->Branch("TrkD0pipt",&TrkD0pipt);
	data->Branch("TrkD0pieta",&TrkD0pieta);
	data->Branch("TrkD0piphi",&TrkD0piphi);
	data->Branch("TrkD0pinhits",&TrkD0pinhits);

	
	if(doMC){	        
	//======================================================
	// MC Variables - D_star
	//======================================================
	data->Branch("MCDseta",&MCDseta);
	data->Branch("MCDsphi",&MCDsphi);
	data->Branch("MCDspt",&MCDspt);
	data->Branch("MCDsenergy",&MCDsenergy);
	data->Branch("MCDsp",&MCDsp);
	data->Branch("MCDset",&MCDset);
	data->Branch("MCDsrapidity",&MCDsrapidity);
	data->Branch("MCDsmass",&MCDsmass);

	data->Branch("MCD0eta",&MCD0eta);
	data->Branch("MCD0phi",&MCD0phi);
	data->Branch("MCD0pt",&MCD0pt);
	data->Branch("MCD0energy",&MCD0energy);
	data->Branch("MCD0p",&MCD0p);
	data->Branch("MCD0et",&MCD0et);
	data->Branch("MCD0rapidity",&MCD0rapidity);
	data->Branch("MCD0mass",&MCD0mass);
	data->Branch("MCD0dispXY",&MCD0dispXY);
	data->Branch("MCD0lifetime",&MCD0lifetime);

	data->Branch("MCDsSeta",&MCDsSeta);
	data->Branch("MCDsSphi",&MCDsSphi);
	data->Branch("MCDsSpt",&MCDsSpt);
	data->Branch("MCDsSenergy",&MCDsSenergy);
	data->Branch("MCDsSp",&MCDsSp);
	data->Branch("MCDsSet",&MCDsSet);
	data->Branch("MCDsSrapidity",&MCDsSrapidity);
	data->Branch("MCDsSmass",&MCDsSmass);
	
	data->Branch("MCDsKeta",&MCDsKeta);
	data->Branch("MCDsKphi",&MCDsKphi);
	data->Branch("MCDsKpt",&MCDsKpt);
	data->Branch("MCDsKenergy",&MCDsKenergy);
	data->Branch("MCDsKp",&MCDsKp);
	data->Branch("MCDsKet",&MCDsKet);
	data->Branch("MCDsKrapidity",&MCDsKrapidity);
	data->Branch("MCDsKmass",&MCDsKmass);

	data->Branch("MCDsPieta",&MCDsPieta);
	data->Branch("MCDsPiphi",&MCDsPiphi);
	data->Branch("MCDsPipt",&MCDsPipt);
	data->Branch("MCDsPienergy",&MCDsPienergy);
	data->Branch("MCDsPip",&MCDsPip);
	data->Branch("MCDsPiet",&MCDsPiet);
	data->Branch("MCDsPirapidity",&MCDsPirapidity);
	data->Branch("MCDsPimass",&MCDsPimass);
	//======================================================
	// MC Variables - D* from B
	//======================================================									
	data->Branch("MCDstarB_eta",&MCDstarB_eta);
	data->Branch("MCDstarB_phi",&MCDstarB_phi);
	data->Branch("MCDstarB_pt",&MCDstarB_pt);
	data->Branch("MCDstarB_energy",&MCDstarB_energy);
	data->Branch("MCDstarB_p",&MCDstarB_p);
	data->Branch("MCDstarB_et",&MCDstarB_et);
	data->Branch("MCDstarB_rapidity",&MCDstarB_rapidity);
	data->Branch("MCDstarB_mass",&MCDstarB_mass);

	data->Branch("MCDstarB_D0eta",&MCDstarB_D0eta);
	data->Branch("MCDstarB_D0phi",&MCDstarB_D0phi);
	data->Branch("MCDstarB_D0pt",&MCDstarB_D0pt);
	data->Branch("MCDstarB_D0energy",&MCDstarB_D0energy);
	data->Branch("MCDstarB_D0p",&MCDstarB_D0p);
	data->Branch("MCDstarB_D0et",&MCDstarB_D0et);
	data->Branch("MCDstarB_D0rapidity",&MCDstarB_D0rapidity);
	data->Branch("MCDstarB_D0mass",&MCDstarB_D0mass);
	data->Branch("MCDstarB_D0dispXY",&MCDstarB_D0dispXY);
	data->Branch("MCDstarB_D0lifetime",&MCDstarB_D0lifetime);

	data->Branch("MCDstarB_Seta",&MCDstarB_Seta);
	data->Branch("MCDstarB_Sphi",&MCDstarB_Sphi);
	data->Branch("MCDstarB_Spt",&MCDstarB_Spt);
	data->Branch("MCDstarB_Senergy",&MCDstarB_Senergy);
	data->Branch("MCDstarB_Sp",&MCDstarB_Sp);
	data->Branch("MCDstarB_Set",&MCDstarB_Set);
	data->Branch("MCDstarB_Srapidity",&MCDstarB_Srapidity);
	data->Branch("MCDstarB_Smass",&MCDstarB_Smass);
	
	data->Branch("MCDstarB_Keta",&MCDstarB_Keta);
	data->Branch("MCDstarB_Kphi",&MCDstarB_Kphi);
	data->Branch("MCDstarB_Kpt",&MCDstarB_Kpt);
	data->Branch("MCDstarB_Kenergy",&MCDstarB_Kenergy);
	data->Branch("MCDstarB_Kp",&MCDstarB_Kp);
	data->Branch("MCDstarB_Ket",&MCDstarB_Ket);
	data->Branch("MCDstarB_Krapidity",&MCDstarB_Krapidity);
	data->Branch("MCDstarB_Kmass",&MCDstarB_Kmass);

	data->Branch("MCDstarB_Pieta",&MCDstarB_Pieta);
	data->Branch("MCDstarB_Piphi",&MCDstarB_Piphi);
	data->Branch("MCDstarB_Pipt",&MCDstarB_Pipt);
	data->Branch("MCDstarB_Pienergy",&MCDstarB_Pienergy);
	data->Branch("MCDstarB_Pip",&MCDstarB_Pip);
	data->Branch("MCDstarB_Piet",&MCDstarB_Piet);
	data->Branch("MCDstarB_Pirapidity",&MCDstarB_Pirapidity);
	data->Branch("MCDstarB_Pimass",&MCDstarB_Pimass);															
	//======================================================
	// MC Variables - D0 
	//======================================================
	data->Branch("MCpromptD0eta",&MCpromptD0eta);
	data->Branch("MCpromptD0phi",&MCpromptD0phi);
	data->Branch("MCpromptD0pt",&MCpromptD0pt);
	data->Branch("MCpromptD0energy",&MCpromptD0energy);
	data->Branch("MCpromptD0p",&MCpromptD0p);
	data->Branch("MCpromptD0et",&MCpromptD0et);
	data->Branch("MCpromptD0rapidity",&MCpromptD0rapidity);
	data->Branch("MCpromptD0mass",&MCpromptD0mass);
	data->Branch("MCpromptD0dispXY",&MCpromptD0dispXY);
	data->Branch("MCpromptD0lifetime",&MCpromptD0lifetime);
	data->Branch("MCpromptD0_DispAngle",&MCpromptD0_DispAngle);

	data->Branch("MCpromptD0_Keta",&MCpromptD0_Keta);
	data->Branch("MCpromptD0_Kphi",&MCpromptD0_Kphi);
	data->Branch("MCpromptD0_Kpt",&MCpromptD0_Kpt);
	data->Branch("MCpromptD0_Kenergy",&MCpromptD0_Kenergy); 
	data->Branch("MCpromptD0_Kp",&MCpromptD0_Kp);
	data->Branch("MCpromptD0_Ket",&MCpromptD0_Ket);
	data->Branch("MCpromptD0_Krapidity",&MCpromptD0_Krapidity);
	data->Branch("MCpromptD0_Kmass",&MCpromptD0_Kmass);

	data->Branch("MCpromptD0_Pieta",&MCpromptD0_Pieta);
	data->Branch("MCpromptD0_Piphi",&MCpromptD0_Piphi);
	data->Branch("MCpromptD0_Pipt",&MCpromptD0_Pipt);
	data->Branch("MCpromptD0_Pienergy",&MCpromptD0_Pienergy);
	data->Branch("MCpromptD0_Pip",&MCpromptD0_Pip);
	data->Branch("MCpromptD0_Piet",&MCpromptD0_Piet);
	data->Branch("MCpromptD0_Pirapidity",&MCpromptD0_Pirapidity);
	data->Branch("MCpromptD0_Pimass",&MCpromptD0_Pimass);
	//======================================================
	// MC Variables - D0 From B
	//======================================================
	data->Branch("MCD0B_D0eta",&MCD0B_D0eta);
	data->Branch("MCD0B_D0phi",&MCD0B_D0phi);
	data->Branch("MCD0B_D0pt",&MCD0B_D0pt);
	data->Branch("MCD0B_D0energy",&MCD0B_D0energy);
	data->Branch("MCD0B_D0p",&MCD0B_D0p);
	data->Branch("MCD0B_D0et",&MCD0B_D0et);
	data->Branch("MCD0B_D0rapidity",&MCD0B_D0rapidity);
	data->Branch("MCD0B_D0mass",&MCD0B_D0mass);
	data->Branch("MCD0B_D0dispXY",&MCD0B_D0dispXY);	
	data->Branch("MCD0B_D0lifetime",&MCD0B_D0lifetime);
	data->Branch("MCD0B_D0DispAngle",&MCD0B_D0DispAngle);

	data->Branch("MCD0B_Keta",&MCD0B_Keta);
	data->Branch("MCD0B_Kphi",&MCD0B_Kphi);
	data->Branch("MCD0B_Kpt",&MCD0B_Kpt);
	data->Branch("MCD0B_Kenergy",&MCD0B_Kenergy); 
	data->Branch("MCD0B_Kp",&MCD0B_Kp);
	data->Branch("MCD0B_Ket",&MCD0B_Ket);
	data->Branch("MCD0B_Krapidity",&MCD0B_Krapidity);
	data->Branch("MCD0B_Kmass",&MCD0B_Kmass);

	data->Branch("MCD0B_Pieta",&MCD0B_Pieta);
	data->Branch("MCD0B_Piphi",&MCD0B_Piphi);
	data->Branch("MCD0B_Pipt",&MCD0B_Pipt);
	data->Branch("MCD0B_Pienergy",&MCD0B_Pienergy);
	data->Branch("MCD0B_Pip",&MCD0B_Pip);
	data->Branch("MCD0B_Piet",&MCD0B_Piet);
	data->Branch("MCD0B_Pirapidity",&MCD0B_Pirapidity);
	data->Branch("MCD0B_Pimass",&MCD0B_Pimass);
	
	}

}

//define this as a plug-in
DEFINE_FWK_MODULE(DstarD0TTree);
