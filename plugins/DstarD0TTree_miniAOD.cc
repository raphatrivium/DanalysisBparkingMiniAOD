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

//GEN MC Matching
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
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

using namespace std; 
using namespace reco;
using namespace edm;

DstarD0TTree::DstarD0TTree(const edm::ParameterSet& iConfig):
	doMC(iConfig.getParameter<bool>("doMC")),
	doRec(iConfig.getParameter<bool>("doRec")),
	debug(iConfig.getUntrackedParameter<bool>("debug",false)),
	selectionCuts(iConfig.getParameter<bool>("selectionCuts")),
	triggerOn(iConfig.getParameter<bool>("triggerOn")),
	triggerName_(iConfig.getUntrackedParameter<std::string>("PathName","HLT_Mu9_IP6_part0_v2")),        
 	//Triggers
	triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("objects"))),
	triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
	trkToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("tracks"))), //MiniAOD
	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("recVtxs"))),
	genParticlesTokenDstar_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gens"))),
	genParticlesTokenD0_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gensD0"))),
	comEnergy_(iConfig.getParameter<double>("comEnergy")),
	DstarSignificance3D_(iConfig.getParameter<double>("DstarSignificance3D")),
	D0Significance3D_(iConfig.getParameter<double>("D0Significance3D"))
{     
	Total_Events_test = 0;
	Triggered_Event_test = 0;
	Ebeam_ = comEnergy_/2.;
	edm::Service<TFileService> fs;

	//Creating trees
	data = fs->make<TTree>("data","data");

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
void DstarD0TTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace std;
	using namespace reco;
	pi_mass=0.13957018; k_mass=0.493677;
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
	
	//Triggers
	Handle<edm::TriggerResults> triggerBits;
	iEvent.getByToken(triggerBits_, triggerBits);

	//Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
	//iEvent.getByToken(triggerObjects_, triggerObjects);

	Handle<pat::PackedTriggerPrescales> triggerPrescales;
	iEvent.getByToken(triggerPrescales_, triggerPrescales);

	//Primary VertexCollection	  
	edm::Handle<VertexCollection> recVtxs;
	iEvent.getByToken(vtxToken_,recVtxs);

	Total_Events++;

	//Only events in which the path actually fired had stored the filter results and products:    
   bool triggerFired = TriggerInfo(iEvent,triggerBits,triggerPrescales,triggerName_);
   if(triggerFired) Triggered_Event_test++;

	NameTrigger.push_back(triggerName_);
	
	// Getting tracks from vert(ex)ices
	edm::Handle< View < pat::PackedCandidate >> tracks; //access that PackedCandidate
	iEvent.getByToken(trkToken_,tracks);

	//miniAOD
	edm::ESHandle<TransientTrackBuilder> theB; 
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
	
	bool canDo= false;
	for( int triggerComb = 0; triggerComb<2; triggerComb++) //For trigger selection on and off
	{
	canDo = false;//initiate the bool each loop.	
	//With trigger
	if (triggerComb == 0 and triggerOn==true and triggerFired ){canDo = true;}
	//Without Trigger
	else if (triggerComb == 1 and triggerOn==false ) {canDo = true; cout << "--------------No Trigger Applied " << endl; }
	//security
	else{ canDo = false; }
	if (canDo){
			
	//Loop in the trakcs of the PackedCandidate 
	for(View<pat::PackedCandidate>::const_iterator iTrack1 = tracks->begin(); iTrack1 != tracks->end(); ++iTrack1 ) 
	{   TotalTracks++;

		if(!iTrack1->hasTrackDetails()) continue;
		TracksHasTrackDetails++;
		if(iTrack1->charge()==0) continue;
		TracksChargeZero++;
		if(fabs(iTrack1->eta())>2.4) continue;
		TracksEta++;
		if(!(iTrack1->trackHighPurity())) continue;
		TracksHighPurity++;
		//if(fabs(iTrack1->pdgId()) != 211 ) continue;
		//TracksPDG211++;

		//Selection for slowpion from D* (The note is 0.5 but the MiniAOD works with pt > 0.5
		if( iTrack1->pt()>0.3 )
		{
			TracksPtZeroFive++;
			if(iTrack1->pseudoTrack().normalizedChi2() > 3.) continue;				  
			TracksChi3++;
			//if(iTrack1->pseudoTrack().hitPattern().numberOfValidHits() < 2) continue;
			if(iTrack1->numberOfHits() < 2) continue; // > and = 2
			TracksNumberOfHits2++;
			if( ( fabs(iTrack1->dxy()) / fabs(iTrack1->dxyError()) ) >3. ) continue;			
			TracksDxyThree++;
			if( ( fabs(iTrack1->dz()) / fabs(iTrack1->dzError()) ) >3. ) continue;
			TracksDzThree++;
			reco::TransientTrack slowpionTT = theB->build(iTrack1->pseudoTrack());
			TrackSlowPionCandidates++; 
			slowPiTracks.push_back(slowpionTT);	//Fill Transient Vector
		}//End Selection for slowpion from D*
			
		Observation++;
		//Kaon and Pion Candidates for D*
		if( iTrack1->pt()>0.5 /*&& (fabs(iTrack1->pdgId()) == 211)*/)
		{
			TracksPtZeroSix++;
			if(iTrack1->pseudoTrack().normalizedChi2() > 2.5) continue;
			TracksChiTwoFive++;
			//cout << "iTrack1->pseudoTrack().normalizedChi2(): "<< iTrack1->pseudoTrack().normalizedChi2() << endl ;
			//cout << "iTrack1->pseudoTrack().Chi2(): "<< iTrack1->pseudoTrack().chi2() << endl ;		
			//if(iTrack1->pseudoTrack().hitPattern().numberOfValidHits() < 5) continue;
			//if(iTrack1->pseudoTrack().hitPattern().numberOfValidPixelHits() < 2) continue;
			if( iTrack1->numberOfHits() < 5) continue;
			TracksNumberOfHits5++;
			if( iTrack1->numberOfPixelHits() < 2) continue;
			TracksNumberOfPixelHits2++;				
			if( fabs(iTrack1->dxy()) > 0.1) continue;
			TracksDxyZeroOne++;	
			if( fabs(iTrack1->dz()) > 1.) continue;
			TracksDzOne++;
			reco::TransientTrack  PionTT = theB->build(iTrack1->pseudoTrack());
			if (debug)cout << " PionTT "  << PionTT.track().momentum() << endl;
			TrackKaonPionCandidates++;
			goodTracks.push_back(PionTT); //Fill Transient Vector
		}//End Kaon and Pion Candidates for D8 

		// SELECTING TRACKS FOR D0 Prompt   
		if( iTrack1->pt()>0.8 /*&& (fabs(iTrack1->pdgId()) == 211)*/)
		{
			TracksPt0p8++;
			if( fabs(iTrack1->eta())>2.5 ) continue;
			D0PromptTracksEta2p5++;	
			if( iTrack1->pseudoTrack().normalizedChi2() > 5.0 ) continue;
			D0PromptTracksChi5++;			
			//if( iTrack1->pseudoTrack().hitPattern().numberOfValidHits() <= 5 ) continue;
			if( iTrack1->numberOfHits() < 5) continue;
			D0PromptTracksNumberHits5++;
			//if( iTrack1->pseudoTrack().hitPattern().numberOfValidPixelHits() <= 2 ) continue;
			if( iTrack1->numberOfPixelHits() < 2) continue;
			D0PromptTracksPixelHits2++;
			//if( iTrack1->p() <1.0 ) continue;
			//D0PromptTracksMomentum1++;
			if( fabs(iTrack1->dxy())>0.1 ) continue;
			D0PromptTracksDxy0p1++;
			if( fabs(iTrack1->dz())>1. ) continue;
			D0PromptTracksDz0p5++;
			reco::TransientTrack  D0TT = theB->build(iTrack1->pseudoTrack());
			goodTracksD0.push_back(D0TT);//Fill Transient Vector
		}//End SELECTING TRACKS FOR D0 Prompt            
				
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
         {
         	VtxIn = i;
         	vtx_trk_size = vtx.tracksSize();
         }
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
      RecD0(iEvent,iSetup,RecVtx); //Reconstruction of prompt D0
      GenDstarInfo(iEvent,iSetup); //Stores information from D0 and its products.
      GenD0Info(iEvent,iSetup); //Stores information from D0 and its products.

		//FindAngle(RecVtx,v_D0,d0kpi_p4); //Calculates the opening angle
      //FindAngleMCpromptD0(p);

		data->Fill();
		}//End canDo	
    }//end trigger flag
}//End analyze

//*********************************************************************************
bool DstarD0TTree::TriggerInfo(const edm::Event& iEvent, edm::Handle<edm::TriggerResults> itriggerBits, edm::Handle<pat::PackedTriggerPrescales> itriggerPrescales, TString trigname){

	const edm::TriggerNames &names = iEvent.triggerNames(*itriggerBits);
	
	//List of Trigger Names in the Dataset
	//for (unsigned int k = 0, nt = itriggerBits->size(); k < nt; ++k) 
	//{std::cout << "Trigger name: "<< names.triggerName(k) << std::endl;}

	//std::cout << "\n == TRIGGER PATHS FOUND =" ;
	for (unsigned int i = 0, n = itriggerBits->size(); i < n; ++i) 
	{	
		TString trigName = names.triggerName(i);
		if(trigName.Contains(trigname)) 
		{ 
			std::cout << "TRIGGER PATHS FOUND: " << names.triggerName(i) <<
			", prescale " << itriggerPrescales->getPrescaleForIndex(i) <<
			": " << (itriggerBits->accept(i) ? "PASS" : "fail (or not run) --");
			
			if (itriggerBits->accept(i) )
			{
				std::cout << "itriggerBits->accept(i): " << itriggerBits->accept(i) << std::endl; std::cout << "---" << std::endl;
				return true;
			}
			else
			{
				std::cout << "itriggerBits->accept(i): " << itriggerBits->accept(i) << std::endl; std::cout << "---" << std::endl;
				return false;
			} 		
      }     	
	}
	std::cout << "**TRIGGER PATHS NOT FOUND**" << std::endl; std::cout << "---" << std::endl;
	return false;
}

//***********************************************************************************
void DstarD0TTree::RecDstar(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::Vertex& RecVtx){

	using namespace std;
	using namespace reco;
	using namespace edm;

   //cout << " RecDstar using goodTracks is " << goodTracks.size() << endl;
	for(size_t i=0;i<goodTracks.size();i++)
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
                                
				//Right Combination Charges
				if(trk1.charge() == trkS.charge()){pi = trk1; K = trk2;}
				else{	K = trk1; pi = trk2;}

				//Wrong Combination Charges - background
				/*if(trk1.charge() == trkS.charge()){pi = trk2;	K = trk1;}
				else{K = trk2;pi = trk1;}	*/
                               
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
				D0AfterLorentzVector++;

				//SlowPion 4-momentum Reconstruction
				math::XYZTLorentzVector p4_S(trkS.track().px(),
														trkS.track().py(),
														trkS.track().pz(),
														sqrt(pow(trkS.track().p(),2)+pow(pi_mass,2)));

				//D* 4-momentum Reconstruction (D0 + SlowPion)
				math::XYZTLorentzVector ip4_DS = ip4_D0 + p4_S;
		
				//To descrase statistc to gain speed. we have more tight cut below
				if( fabs(ip4_D0.M()-1.86484)  > 0.2) continue;
				D0MinusPDGzero2++;

				//To descrase statistc to gain speed. we have more tight cut below
				if((ip4_DS.M() - ip4_D0.M()) > 0.3) continue;
 				//cout << "ip4_DS.M() :" << ip4_DS.M() << " ip4_D0.M(): " << ip4_D0.M() << endl;
				DsMinusD0Zerothree++;
				                      
				//KalmanVertexFitter
				vector<TransientTrack> tks;
				tks.push_back(K);
				tks.push_back(pi);
				KalmanVertexFitter kalman(true);
				TransientVertex v = kalman.vertex(tks);
				if(!v.isValid() || !v.hasRefittedTracks()) continue;
				double vtxProb =TMath::Prob( (Double_t) v.totalChiSquared(), (Int_t) v.degreesOfFreedom());
				TransientTrack K_f = v.refittedTrack(K);
				TransientTrack pi_f = v.refittedTrack(pi);
				TransientTrackOfpiK++;
				
				//SV Confidence Level
				if(selectionCuts) {if(vtxProb < 0.01) continue;}
				SVConfidenceLevel++;
				
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
				D0AfterLorentzVectorKarman++;
				
				//Angle between formation and decay of D0
  				double anglephi = FindAngle(RecVtx,v,d0_p4);
				double cosPhi = cos(anglephi);
				if(selectionCuts) {if( cosPhi < 0.99 ) continue;}
				//cout << "cos(anglephi)	: " << cos(anglephi) << endl;
				PointingcosPhi++;

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
				if(selectionCuts) {if ( D0fromDSs3D < DstarSignificance3D_ ) continue;}
				Significance++;

				//Difference between D0 and D0PDG 5sigmas
				if( fabs(d0mass - 1.86484) > 0.1 ) continue;
				D0MinusPDG++;

				//Pt Cut of mesnons D0
				if(selectionCuts) {if ( d0_p4.Pt() < 3. ) continue;}
				D0pTThree++;
								
				D0Candidates++;

				math::XYZTLorentzVector dS_p4 = d0_p4 + p4_S;
				double dsmass = dS_p4.M();
				DsAfterLorentzVector++;
								                      											
				//Difference between D* and D0 -> Must be close to pion
				if(selectionCuts) {if( (dsmass - d0mass) > 0.16) continue;}
				DsMinusD0++;
				DsCandidates++; //Number of D* candidates

				//Flag
				if(doMC)
				{
        			Handle<GenParticleCollection> gens; //linked to "prunedGenParticles" in config file
        			iEvent.getByToken(genParticlesTokenDstar_, gens);	 

		  			for(size_t i=0; i<dScandsKpi.size();i++)
					{
			    		const GenParticle & ds = gens ->at(dScandsKpi.at(i));
			    		double massDS = ds.mass();
			    		if (massDS>0) FlagMC=1;
			  		}
				}
						                       		
				//Fill vectors
				D0_VtxProb.push_back(vtxProb);
				D0mass.push_back(d0_p4.M());
				Dsmass.push_back(dS_p4.M());
				D0pt.push_back(d0_p4.Pt());
				Dspt.push_back(dS_p4.Pt());
				D0eta.push_back(d0_p4.eta());
				D0phi.push_back(d0_p4.phi());
				Dseta.push_back(dS_p4.eta());
				Dsphi.push_back(dS_p4.phi());

				D0_VtxPosx.push_back(v.position().x());
				D0_VtxPosy.push_back(v.position().y());
				D0_VtxPosz.push_back(v.position().z());
				D0_Vtxerrx.push_back(v.positionError().cxx());
				D0_Vtxerry.push_back(v.positionError().cyy());
				D0_Vtxerrz.push_back(v.positionError().czz());

				TrkKmass.push_back(p4_K.M());
				Trkpimass.push_back(p4_pi.M());
				TrkSmass.push_back(p4_S.M());

				TrkKdxy.push_back(K_f.track().dxy(RecVtx.position()));
				Trkpidxy.push_back(pi_f.track().dxy(RecVtx.position()));
				TrkSdxy.push_back(trkS.track().dxy(RecVtx.position()));

				TrkKdz.push_back(K_f.track().dz(RecVtx.position()));
				Trkpidz.push_back(pi_f.track().dz(RecVtx.position()));
				TrkSdz.push_back(trkS.track().dz(RecVtx.position()));

				TrkKnhits.push_back(K.track().numberOfValidHits());
				Trkpinhits.push_back(pi.track().numberOfValidHits());
				TrkSnhits.push_back(trkS.track().numberOfValidHits());

				TrkKchi2.push_back(K.track().normalizedChi2());
				Trkpichi2.push_back(pi.track().normalizedChi2());
				TrkSchi2.push_back(trkS.track().normalizedChi2());

				DSDeltaR.push_back(deltaR(d0_p4.eta(),d0_p4.phi(),trkS.track().eta(),trkS.track().phi()));

				TrkKpt.push_back(K_f.track().pt());
				Trkpipt.push_back(pi_f.track().pt());
				TrkSpt.push_back(trkS.track().pt());

				TrkKeta.push_back(K_f.track().eta());
				Trkpieta.push_back(pi_f.track().eta());
				TrkSeta.push_back(trkS.track().eta());

				TrkKphi.push_back(K_f.track().phi());
				Trkpiphi.push_back(pi_f.track().phi());
				TrkSphi.push_back(trkS.track().phi());

				TrkScharge.push_back(trkS.charge());

				D0fromDSs3D_vec.push_back(D0fromDSs3D);
				D0fromDSsXY_vec.push_back(D0fromDSsXY);
				Anglephi_vec.push_back(cosPhi);   
                               
				NKpiCand++;  
									
				if(NKpiCand>999) break;// this a control counter. We do not expect high rate of D*
			}
			if(NKpiCand>999) break;
	 } 
	if(NKpiCand>999) break;
	}
}//End RecDstar

//***********************************************************************************
void DstarD0TTree::RecDstarWrongCombination(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::Vertex& RecVtx){

	using namespace std;
	using namespace reco;
	using namespace edm;

   //cout << " RecDstar using goodTracks is " << goodTracks.size() << endl;
	for(size_t i=0;i<goodTracks.size();i++)
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
				//D0AfterLorentzVector++;

				//SlowPion 4-momentum Reconstruction
				math::XYZTLorentzVector p4_S(trkS.track().px(),
														trkS.track().py(),
														trkS.track().pz(),
														sqrt(pow(trkS.track().p(),2)+pow(pi_mass,2)));

				//D* 4-momentum Reconstruction (D0 + SlowPion)
				math::XYZTLorentzVector ip4_DS = ip4_D0 + p4_S;
		
				//To descrase statistc to gain speed. we have more tight cut below
				if( fabs(ip4_D0.M()-1.86484)  > 0.2) continue;
				//D0MinusPDGzero2++;

				//To descrase statistc to gain speed. we have more tight cut below
				if((ip4_DS.M() - ip4_D0.M()) > 0.3) continue;
 				//cout << "ip4_DS.M() :" << ip4_DS.M() << " ip4_D0.M(): " << ip4_D0.M() << endl;
				//DsMinusD0Zerothree++;
				                      
				//KalmanVertexFitter
				vector<TransientTrack> tks;
				tks.push_back(K);
				tks.push_back(pi);
				KalmanVertexFitter kalman(true);
				TransientVertex v = kalman.vertex(tks);
				if(!v.isValid() || !v.hasRefittedTracks()) continue;
				double vtxProb =TMath::Prob( (Double_t) v.totalChiSquared(), (Int_t) v.degreesOfFreedom());
				TransientTrack K_f = v.refittedTrack(K);
				TransientTrack pi_f = v.refittedTrack(pi);
				//TransientTrackOfpiK++;
				
				//SV Confidence Level
				if(selectionCuts) {if(vtxProb < 0.01) continue;}
				//SVConfidenceLevel++;
				
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
				//D0AfterLorentzVectorKarman++;
				
				//Angle between formation and decay of D0
  				double anglephi = FindAngle(RecVtx,v,d0_p4);
				double cosPhi = cos(anglephi);
				if(selectionCuts) {if( cosPhi < 0.99 ) continue;}
				//cout << "cos(anglephi)	: " << cos(anglephi) << endl;
				//PointingcosPhi++;

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
				//Significance++;

				//Difference between D0 and D0PDG 5sigmas
				if( fabs(d0mass - 1.86484) > 0.1 ) continue;
				//D0MinusPDG++;

				//Pt Cut of mesnons D0
				if(selectionCuts) {if( d0_p4.Pt() < 3. ) continue;}
				//D0pTThree++;
								
				//D0Candidates++;

				math::XYZTLorentzVector dS_p4 = d0_p4 + p4_S;
				double dsmass = dS_p4.M();
				DsAfterLorentzVector++;
								                      											
				//Difference between D* and D0 -> Must be close to pion
				if(selectionCuts) {if( (dsmass - d0mass) > 0.16) continue;}
				//DsMinusD0++;
				//DsCandidates++; //Number of D* candidates
			                       		
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

}//RecDstarWrongCombination

//***********************************************************************************
//Evaluate if a particle is stable (status 1). If not it repeat the process until find the a stable one.
void DstarD0TTree::assignStableDaughters(const reco::Candidate* p, std::vector<int> & pids){
	for(size_t i=0; i< p->numberOfDaughters(); i++)
	{
		if(p->daughter(i)->status() == 1)	pids.push_back(abs(p->daughter(i)->pdgId()));
		else	assignStableDaughters(p->daughter(i), pids);
		//cout << p->daughter(i)->pdgId() << endl;
	}
	return;
}

//***********************************************************************************
void DstarD0TTree::GenDstarInfo(const edm::Event& iEvent,const edm::EventSetup& iSetup){

	if(doMC)
	{ 
		using namespace std;
		using namespace reco;
  		using namespace edm;
		//Handle<GenParticleCollection> genParticles; 
		edm::Handle<GenParticleCollection> gens; //linked to "prunedGenParticles" in config file
		iEvent.getByToken(genParticlesTokenDstar_, gens);

		for(size_t i=0;i<gens->size();i++)
		{
     		const GenParticle & p = (*gens)[i];
      	if(fabs(p.pdgId()) == 413) //D*
			{ 
        		for(size_t j=0; j < p.numberOfDaughters(); j++)
				{	
	          	const Candidate* dau = p.daughter(j);
	    			if(fabs(dau->pdgId()) == 421) //D0
					{ 
						std::vector<int> d0dauspids;
			        	assignStableDaughters(dau,d0dauspids);
        	  			int K_num = 0, pi_num = 0, n_daughters = d0dauspids.size(); 
	  				   while (!d0dauspids.empty())
						{
			    	   	int pid = d0dauspids.back(); //Access last element 
            			if(pid==321) K_num++;
			           	if(pid==211) pi_num++;
            			d0dauspids.pop_back(); //Delete last element 
						}		
						if(K_num == 1 && pi_num == 1 && n_daughters == 2)
						{         
							dScandsKpi.push_back(i); 
							MCDseta.push_back(p.eta()); 
							MCDsphi.push_back(p.phi());
							MCDspt.push_back(p.pt());
							MCDsenergy.push_back(p.energy());
							MCDsp.push_back(p.p());
							MCDset.push_back(p.et());
							MCDsrapidity.push_back(p.rapidity());
							MCDsmass.push_back(p.mass());           

							MCD0eta.push_back(dau->eta());
							MCD0phi.push_back(dau->phi());
							MCD0pt.push_back(dau->pt());
							MCD0energy.push_back(dau->energy());
							MCD0p.push_back(dau->p());
							MCD0et.push_back(dau->et());
				        	MCD0rapidity.push_back(dau->rapidity());   
							MCD0mass.push_back(dau->mass());
            		}		     
          		}
				}
      	}
    	}

		NdsKpiMC = dScandsKpi.size();
		//cout << "NdsKpiMC: "<< NdsKpiMC << endl;
   }

}//GenDstarInfo

//***********************************************************************************
void DstarD0TTree::RecD0(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::Vertex& RecVtx){

	using namespace edm;
	using namespace std;
	using namespace reco;

	for(size_t i = 0; i < goodTracksD0.size(); i++)
	{
		TransientTrack trk1D0 = goodTracksD0[i];

		for(size_t j=i+1;j<goodTracksD0.size();j++)
		{
			TracksD0combination++;
			TransientTrack trk2D0 = goodTracksD0[j];
			//Testing charge and if tracks are equal
			if(trk1D0 == trk2D0) continue;
			if(trk1D0.charge() == trk2D0.charge()) continue;

			//Recosntruction of D0 momentum
			math::XYZVector D0_p = trk1D0.track().momentum() + trk2D0.track().momentum();

			//TransientTrack KfromD0 = 0, PifromD0 = 0;	
			TransientTrack KfromD0 , PifromD0 ;	        
            	
			mass1 = mass2 = 0. ;

			//Reconstructio of 4-vector D0 Prompt for track1 = Kaon and track2 = pion
			vD0kaon.SetPtEtaPhiM(trk1D0.track().pt(), 
										trk1D0.track().eta(), 
										trk1D0.track().phi(), 
										k_mass);
			vD0pion.SetPtEtaPhiM(trk2D0.track().pt(), 
										trk2D0.track().eta(), 
										trk2D0.track().phi(), 
										pi_mass);
			vD0_1 = vD0kaon + vD0pion;
			mass1 = vD0_1.M();

			//Reconstruction of 4-vector D0 Prompt track1 = pion and track2 = Kaon
			vD0pion.SetPtEtaPhiM(trk1D0.track().pt(),
										trk1D0.track().eta(),
							 			trk1D0.track().phi(), 
										pi_mass);
			vD0kaon.SetPtEtaPhiM(trk2D0.track().pt(), 
										trk2D0.track().eta(), 
										trk2D0.track().phi(), 
										k_mass);
			vD0_2 = vD0pion + vD0kaon;
			mass2 = vD0_2.M();

			comb1 = comb2 = false ;
			if( fabs(mass1-1.86484) < 0.2) {comb1 = true;}		
			if( fabs(mass2-1.86484) < 0.2) {comb2 = true;}

			//=============================
			//combinations						
			for( int icomb = 0; icomb<3; icomb++) 
			{	//Combination 1
				if (icomb == 0 and comb1 and !comb2)
				{  KfromD0 = trk1D0; PifromD0 = trk2D0;}
				//Combination 2
				else if (icomb == 1 and comb2 and !comb1)	
				{  KfromD0 = trk2D0; PifromD0 = trk1D0;}
				//Combination 1 and 2
				else if (icomb == 2 and comb1 and comb2)		
				{ 	if(fabs( mass1-1.864 ) < fabs( mass2-1.864 ))
					{KfromD0 = trk1D0; PifromD0 = trk2D0;}
					else
					{KfromD0 = trk2D0; PifromD0 = trk1D0;} 
				}		
				else continue;

				D0PromptminusPDG1p0++;
				math::XYZVector K_p = trk2D0.track().momentum();

				vector<TransientTrack> tksD0;
				tksD0.push_back(KfromD0);
				tksD0.push_back(PifromD0);
				KalmanVertexFitter kalman(true);
				v_D0 = kalman.vertex(tksD0);
				//TransientVertex v_D0 = kalman.vertex(tksD0);
				if(!v_D0.isValid() || !v_D0.hasRefittedTracks()) continue;
				double D0KpivtxProb =TMath::Prob( (Double_t) v_D0.totalChiSquared(), (Int_t) v_D0.degreesOfFreedom());
				TransientTrack KfromD0_f = v_D0.refittedTrack(KfromD0);
				TransientTrack pifromD0_f = v_D0.refittedTrack(PifromD0);
				D0PromptKpiAfterTransientp0++;

				//SV Confidence Level
				if(selectionCuts) {if(D0KpivtxProb < 0.01) continue;}
				D0PromptSVConfidenceLevel++;


				math::XYZTLorentzVector p4_KfromD0(KfromD0_f.track().px(),
																KfromD0_f.track().py(),
																KfromD0_f.track().pz(),
																sqrt(pow(KfromD0_f.track().p(),2)+pow(k_mass,2))  );
				math::XYZTLorentzVector p4_pifromD0(pifromD0_f.track().px(),
																pifromD0_f.track().py(),
																pifromD0_f.track().pz(),
																sqrt(pow(pifromD0_f.track().p(),2)+pow(pi_mass,2)) );  
				d0kpi_p4 = p4_KfromD0 + p4_pifromD0;

				//Angle between formation and decay of D0
				double dispAngle = FindAngle(RecVtx,v_D0,d0kpi_p4);
				double D0cosPhi = cos(dispAngle);
				if(selectionCuts) {if( D0cosPhi < 0.99 ) continue;}
				D0PromptPointingCosPhi++;
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

				//cout << "significance 3D: "<< D0Kpis3D << endl;
				if(selectionCuts) {if( D0Kpis3D < D0Significance3D_ ) continue;}
				D0PromptSignificance++;

				//Vetorial product of 4-momentum Kaon and 4-momentum D0
				double D0_kT = sqrt( (K_p).Cross(D0_p).Mag2() / D0_p.Mag2() ) ;
				double d0kpimass = d0kpi_p4.M();
				//cout << "d0kpimass comb1: " << d0kpimass << endl;
				if(selectionCuts) {if(fabs(d0kpimass - 1.86484)>0.15) continue;}
				D0PromptCandidates++;
				ND0KpiCand++;

				cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$parte9 "<< "mass1: "<< mass1 << "  ## mass2: "<< mass2 << endl;

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
				D0Kpi_DispAngle.push_back(dispAngle);
				TrkD0Kdxy.push_back(KfromD0_f.track().dxy(RecVtx.position()));
				TrkD0pidxy.push_back(pifromD0_f.track().dxy(RecVtx.position()));

				TrkD0Kdz.push_back(KfromD0_f.track().dz(RecVtx.position())); 
				TrkD0pidz.push_back(pifromD0_f.track().dz(RecVtx.position()));
				TrkD0Knhits.push_back(KfromD0.track().numberOfValidHits());
				TrkD0Knhits.push_back(KfromD0.track().numberOfValidHits());
				TrkD0pinhits.push_back(PifromD0.track().numberOfValidHits());

				TrkD0Kchi2.push_back(KfromD0.track().normalizedChi2());
				TrkD0pichi2.push_back(PifromD0.track().normalizedChi2());
				TrkD0Kpt.push_back(KfromD0_f.track().pt());
				TrkD0pipt.push_back(pifromD0_f.track().pt());

				TrkD0Keta.push_back(KfromD0_f.track().eta());
				TrkD0pieta.push_back(pifromD0_f.track().eta());
				TrkD0Kphi.push_back(KfromD0_f.track().phi());
				TrkD0piphi.push_back(pifromD0_f.track().phi());
				D0KpisXY_vec.push_back(D0KpisXY);      
				D0Kpis3D_vec.push_back(D0Kpis3D);
				D0_kT_vec.push_back(D0_kT);	
			
			}//end Combinations
		}//second track loop
	}//first track loop
} //END OF RECD0 

//***********************************************************************************
void DstarD0TTree::GenD0Info(const edm::Event& iEvent,const edm::EventSetup& iSetup){

	if(doMC)
	{
		using namespace std;
		using namespace reco;
		using namespace edm;

		edm::Handle<GenParticleCollection> gensD0;
		iEvent.getByToken(genParticlesTokenD0_, gensD0);

		for(size_t i=0; i< gensD0->size(); i++)
		{
			const GenParticle & p = (*gensD0)[i];
			//const Candidate * mom = p.mother();
			//cout << "gensD0 pdgID " << p.pdgId() << endl;
				
			if(fabs(p.pdgId()) == 421 && p.numberOfDaughters() == 2) //PDG D0
			{ 
				int IDdau1 = p.daughter(0)->pdgId() ; int IDdau2 = p.daughter(1)->pdgId(); 

				if(fabs(IDdau1) == 321 && fabs(IDdau2) == 211 && IDdau1*IDdau2 < 0)
				{
					//cout << "Passa if 1" << endl;
	
					MCpromptD0eta.push_back(p.eta());
					MCpromptD0phi.push_back(p.phi());
					MCpromptD0pt.push_back(p.pt());
					MCpromptD0energy.push_back(p.energy());
					MCpromptD0p.push_back(p.p());
					MCpromptD0et.push_back(p.et());
					MCpromptD0rapidity.push_back(p.rapidity());   
					MCpromptD0mass.push_back(p.mass());					

					MCpromptD0_Keta.push_back(p.daughter(0)->eta());
					MCpromptD0_Kphi.push_back(p.daughter(0)->phi());
					MCpromptD0_Kpt.push_back(p.daughter(0)->pt());
					MCpromptD0_Kenergy.push_back(p.daughter(0)->energy());
					MCpromptD0_Kp.push_back(p.daughter(0)->p());
					MCpromptD0_Ket.push_back(p.daughter(0)->et());
					MCpromptD0_Krapidity.push_back(p.daughter(0)->rapidity());
					MCpromptD0_Kmass.push_back(p.daughter(0)->mass());
 			
					MCpromptD0_Pieta.push_back(p.daughter(1)->eta());
					MCpromptD0_Piphi.push_back(p.daughter(1)->phi());
					MCpromptD0_Pipt.push_back(p.daughter(1)->pt());
					MCpromptD0_Pienergy.push_back(p.daughter(1)->energy());
					MCpromptD0_Pip.push_back(p.daughter(1)->p());
					MCpromptD0_Piet.push_back(p.daughter(1)->et());
					MCpromptD0_Pirapidity.push_back(p.daughter(1)->rapidity());
					MCpromptD0_Pimass.push_back(p.daughter(1)->mass());

					//cout << "p.mother(0) = " << p.mother(0)->pdgId() << "p.pdgId() = " << p.pdgId() << endl;					

					CLHEP::Hep3Vector MCdisplacement(p.mother(0)->vx() - p.vx(), 
																p.mother(0)->vy() - p.vy(), 
																p.mother(0)->vz() - p.vz() );
					CLHEP::Hep3Vector momentum( p.px(), p.py(), p.pz() );
					//cout << " MCdisplacement = " << MCdisplacement << " momentum = " << momentum << endl;
					double MCD0dispAngle = momentum.angle(MCdisplacement); 						

					MCpromptD0_DispAngle.push_back(MCD0dispAngle);
					
					//cout << "passa if 2" << endl;
				}

				if(fabs(IDdau1) == 211 && fabs(IDdau2) == 321 && IDdau1*IDdau2 < 0 )
				{
					//cout << "Passa if 3" << endl;

					MCpromptD0eta.push_back(p.eta());
					MCpromptD0phi.push_back(p.phi());
					MCpromptD0pt.push_back(p.pt());
					MCpromptD0energy.push_back(p.energy());
					MCpromptD0p.push_back(p.p());
					MCpromptD0et.push_back(p.et());
					MCpromptD0rapidity.push_back(p.rapidity());  
					MCpromptD0mass.push_back(p.mass());

					MCpromptD0_Keta.push_back(p.daughter(1)->eta());
					MCpromptD0_Kphi.push_back(p.daughter(1)->phi());
					MCpromptD0_Kpt.push_back(p.daughter(1)->pt());
					MCpromptD0_Kenergy.push_back(p.daughter(1)->energy());
					MCpromptD0_Kp.push_back(p.daughter(1)->p());
					MCpromptD0_Ket.push_back(p.daughter(1)->et());
					MCpromptD0_Krapidity.push_back(p.daughter(1)->rapidity());
					MCpromptD0_Kmass.push_back(p.daughter(1)->mass());
	
					MCpromptD0_Pieta.push_back(p.daughter(0)->eta());
					MCpromptD0_Piphi.push_back(p.daughter(0)->phi());
					MCpromptD0_Pipt.push_back(p.daughter(0)->pt());
					MCpromptD0_Pienergy.push_back(p.daughter(0)->energy());
					MCpromptD0_Pip.push_back(p.daughter(0)->p());
					MCpromptD0_Piet.push_back(p.daughter(0)->et());
					MCpromptD0_Pirapidity.push_back(p.daughter(0)->rapidity());
					MCpromptD0_Pimass.push_back(p.daughter(0)->mass());

					CLHEP::Hep3Vector MCdisplacement(p.mother(0)->vx() - p.vx(),
																p.mother(0)->vy() - p.vy(),
																p.mother(0)->vz() - p.vz() );
					CLHEP::Hep3Vector momentum( p.px() , p.py() , p.pz() );
					double MCD0dispAngle = momentum.angle(MCdisplacement) ;

					MCpromptD0_DispAngle.push_back(MCD0dispAngle);

					//cout << "Passa if 4" << endl;   
				}
			}
		}
	}
}//GenD0Info

//***********************************************************************************
double DstarD0TTree::FindAngle(const reco::Vertex& pv , const TransientVertex& sv, const math::XYZTLorentzVector& vD0angle )
{
	CLHEP::Hep3Vector displacement( 	sv.position().x() - pv.position().x(),
												sv.position().y() - pv.position().y(),
												sv.position().z() - pv.position().z() ); 
	CLHEP::Hep3Vector momentum( vD0angle.Px() , vD0angle.Py() , vD0angle.Pz() ) ; 
	return momentum.angle(displacement) ;
}

/*double DstarD0TTree::FindAngleMCpromptD0(const GenParticle& genD0)
  {
  CLHEP::Hep3Vector MCdisplacement( genD0.xProd() - genD0.xDec() ,
  genD0.yProd() - genD0.yProd() ,
  genD0.zProd() - genD0.zProd() ) ;
  CLHEP::Hep3Vector momentum( genD0.px() , genD0.py() , genD0.pz() ) ;
  return momentum.angle(MCdisplacement) ;	
  }*/

//***********************************************************************************
void DstarD0TTree::initialize(){

	dScandsKpi.clear(); goodTracks.clear(); 
	goodTracksD0.clear(); slowPiTracks.clear();
	NKpiCand=0; D0Candidates=0; DsCandidates=0; FlagMC=0; NdsKpiMC=0; 
	//PVx = PVy = PVz = PVerrx = PVerry = PVerrz = -999.;
	//ntracksD0Kpi = -999; ntracksDstar = -999; n_pVertex = -999; 
	runNumber=0; eventNumber=0; lumi=0; HLTPath_=0; 
	HFEnergyPlus=0.; HFEnergyMinus=0.;  lumiWeight_=0.;

	//counters
	Total_Events = 0; TotalTracks = 0; TracksHasTrackDetails =0; TracksChargeZero = 0; TracksEta = 0; 
	TracksHighPurity = 0; TracksPDG211 = 0; TracksPtZeroFive = 0; TracksChi3 = 0; TracksNumberOfHits2 = 0; 
	TracksDxyThree = 0; TracksDzThree = 0; TrackSlowPionCandidates = 0; TracksPtZeroSix = 0; TracksChiTwoFive = 0;
	TracksNumberOfHits5 = 0; TracksNumberOfPixelHits2 = 0; TracksDxyZeroOne = 0; TracksDzOne = 0; TrackKaonPionCandidates = 0;
	Observation = 0;

	//D0 Prompt Conters
	TracksPt0p8 = 0; D0PromptTracksEta2p5 = 0; D0PromptTracksChi5 =0; D0PromptTracksNumberHits5 = 0;
	D0PromptTracksPixelHits2 = 0; D0PromptTracksDxy0p1 = 0; D0PromptTracksDz0p5 =0;
	TracksD0combination = 0; D0PromptminusPDG1p0 = 0; D0PromptKpiAfterTransientp0 =0;
	D0PromptSVConfidenceLevel = 0; D0PromptPointingCosPhi = 0; D0PromptSignificance =0;
	D0PromptCandidates = 0; ND0KpiCand = 0;
	
	D0AfterLorentzVector = 0; D0MinusPDGzero2 = 0; DsMinusD0Zerothree = 0; TransientTrackOfpiK = 0; SVConfidenceLevel = 0; D0AfterLorentzVectorKarman = 0; 
	PointingcosPhi = 0;	Significance = 0; D0pTThree = 0; DsAfterLorentzVector = 0; D0MinusPDG = 0; DsMinusD0 = 0;

	NameTrigger.clear();
	D0_VtxProb.clear(); D0mass.clear(); TrkKmass.clear(); Trkpimass.clear(); TrkSmass.clear(); 
	Dsmass.clear(); D0pt.clear(); Dspt.clear(); D0eta.clear(); D0phi.clear(); Dseta.clear(); 
	Dsphi.clear(); D0_VtxPosx.clear(); D0_VtxPosy.clear(); D0_VtxPosz.clear(); D0_Vtxerrx.clear(); 
	D0_Vtxerry.clear(); D0_Vtxerrz.clear(); TrkKdxy.clear(); Trkpidxy.clear();
	TrkSdxy.clear(); TrkKdz.clear(); Trkpidz.clear();TrkSdz.clear(); TrkKnhits.clear(); 
	Trkpinhits.clear(); TrkSnhits.clear();	TrkKchi2.clear(); Trkpichi2.clear(); TrkSchi2.clear();
	DSDeltaR.clear(); TrkKpt.clear(); Trkpipt.clear(); TrkSpt.clear(); TrkKeta.clear(); 
	Trkpieta.clear(); TrkSeta.clear(); TrkKphi.clear(); Trkpiphi.clear(); TrkSphi.clear(); TrkScharge.clear();
	D0fromDSsXY_vec.clear(); D0fromDSs3D_vec.clear(); Anglephi_vec.clear();

	D0_VtxProbWrong.clear(); D0massWrong.clear(); TrkKmassWrong.clear(); TrkpimassWrong.clear(); 
	TrkSmassWrong.clear(); DsmassWrong.clear(); D0ptWrong.clear(); DsptWrong.clear();D0etaWrong.clear();
	D0phiWrong.clear(); DsetaWrong.clear(); DsphiWrong.clear(); D0_VtxPosxWrong.clear(); D0_VtxPosyWrong.clear();
	D0_VtxPoszWrong.clear(); D0_VtxerrxWrong.clear(); D0_VtxerryWrong.clear();
	D0_VtxerrzWrong.clear(); TrkKdxyWrong.clear(); TrkpidxyWrong.clear();
	TrkSdxyWrong.clear(); TrkKdzWrong.clear();TrkpidzWrong.clear();TrkSdzWrong.clear();
	TrkKnhitsWrong.clear(); TrkpinhitsWrong.clear(); TrkSnhitsWrong.clear();
	TrkKchi2Wrong.clear(); Trkpichi2Wrong.clear(); TrkSchi2Wrong.clear();
	DSDeltaRWrong.clear(); TrkKptWrong.clear(); TrkpiptWrong.clear(); TrkSptWrong.clear();
	TrkKetaWrong.clear();TrkpietaWrong.clear();TrkSetaWrong.clear();TrkKphiWrong.clear();TrkpiphiWrong.clear();TrkSphiWrong.clear();TrkSchargeWrong.clear();
	D0fromDSsXY_vecWrong.clear(); D0fromDSs3D_vecWrong.clear(); Anglephi_vecWrong.clear();


	D0Kpi_VtxProb.clear();D0Kpimass.clear();D0Kpipt.clear();D0Kpieta.clear();
	D0Kpiphi.clear();D0Kpi_VtxPosx.clear();D0Kpi_VtxPosy.clear();
	D0Kpi_VtxPosz.clear();D0Kpi_Vtxerrx.clear();
	D0Kpi_Vtxerry.clear();D0Kpi_Vtxerrz.clear();
	D0Kpi_DispAngle.clear();TrkD0Kdxy.clear();
	TrkD0pidxy.clear();TrkD0Kdz.clear();
	TrkD0pidz.clear();TrkD0Knhits.clear();
	TrkD0pinhits.clear();TrkD0Kchi2.clear();TrkD0pichi2.clear();
	TrkD0Kpt.clear();TrkD0pipt.clear();TrkD0Keta.clear();
	TrkD0pieta.clear();TrkD0Kphi.clear();TrkD0piphi.clear();
	D0KpisXY_vec.clear();
	D0Kpis3D_vec.clear();
	D0_kT_vec.clear();	
	//	D0KpisXY_.clear();

	MxFromPFCands_.clear(); EPlusPzFromPFCands_.clear(); EMinusPzFromPFCands_.clear();sumEHFPlusFromPFCands_.clear(); sumEHFMinusFromPFCands_.clear();
	xiPlusFromPFCands_.clear(); xiMinusFromPFCands_.clear(); etaMaxFromPFCands_.clear(); etaMinFromPFCands_.clear(); missingMassFromXiFromPFCands_.clear();
	etaMaxFromPFCandsNew_.clear(); etaMinFromPFCandsNew_.clear();  MCDseta.clear(); MCDsphi.clear(); MCDspt.clear(); MCDsenergy.clear(); MCDsp.clear(); MCDset.clear(); MCDsrapidity.clear(); MCDsmass.clear(); 
	MCDseta.clear(); MCDsphi.clear(); MCDspt.clear(); MCDsenergy.clear(); MCDsp.clear(); MCDset.clear(); MCDsrapidity.clear(); MCDsmass.clear(); 
	MCD0eta.clear(); MCD0phi.clear(); MCD0pt.clear(); MCD0energy.clear(); MCD0p.clear(); MCD0et.clear(); MCD0rapidity.clear(); MCD0mass.clear(); 

	MCpromptD0eta.clear(); MCpromptD0phi.clear(); MCpromptD0pt.clear(); MCpromptD0energy.clear(); MCpromptD0p.clear(); MCpromptD0et.clear(); MCpromptD0rapidity.clear(); MCpromptD0mass.clear();
	MCpromptD0_Keta.clear(); MCpromptD0_Kphi.clear(); MCpromptD0_Kpt.clear(); MCpromptD0_Kenergy.clear(); MCpromptD0_Kp.clear(); MCpromptD0_Ket.clear(); MCpromptD0_Krapidity.clear(); MCpromptD0_Kmass.clear();
	MCpromptD0_Pieta.clear(); MCpromptD0_Piphi.clear(); MCpromptD0_Pipt.clear(); MCpromptD0_Pienergy.clear(); MCpromptD0_Pip.clear(); MCpromptD0_Piet.clear(); MCpromptD0_Pirapidity.clear(); MCpromptD0_Pimass.clear();

	MCpromptD0_DispAngle.clear(); MCpromptD0_Kt.clear();
	/*  xiGenPlus_.clear(); xiGenMinus_.clear(); MxGen_.clear(); MxGenRange_.clear(); sumEnergyHEPlusGen_.clear(); 
	    sumEnergyHEMinusGen_.clear(); sumEnergyHFPlusGen_.clear(); MxGenMinus_.clear();
	    sumEnergyHFMinusGen_.clear(); etaMaxGen_.clear(); etaMinGen_.clear(); deltaEtaGen_.clear(); etaGapLow_.clear(); etaGapHigh_.clear(); MxGenPlus_.clear();*/

	pfsis1Eta_max.clear(); pfsis2Eta_max.clear(); pfsis1Eta_min.clear(); pfsis2Eta_min.clear(); deltaEtapf.clear();

	//triggers.clear();

}
//++++++++++++++++++
void DstarD0TTree::endJob(){

	cout <<"######################################################################"<<endl;
	cout << "Number of Events: " << eventNumber << " Run Number: " << runNumber << endl;
	cout << "Total events: " << Total_Events_test << " # Total events triggered by " << triggerName_ << ": " << Triggered_Event_test << endl;
	cout << "D* Candidates: " << DsCandidates << " # D0 Candidates: " << D0PromptCandidates << endl;
    
}

//***********************************************************************************
void DstarD0TTree::beginJob(){
        
    //cout << "DsCandidates: " << DsCandidates << " NKpiCand: " << NKpiCand << endl;
	
	data->Branch("Total_Events",&Total_Events,"Total_Events/I"); 
	data->Branch("Triggered_Event",&Triggered_Event,"Triggered_Event/I");
	data->Branch("TriggerName", &NameTrigger);
	data->Branch("NdsKpiMC",&NdsKpiMC,"NdsKpiMC/I");

	data->Branch("runNumber",&runNumber,"runNumber/I");
	data->Branch("eventNumber",&eventNumber,"eventNumber/I");
	data->Branch("lumi",&lumi,"lumi/I");
	data->Branch("lumiWeight_",&lumiWeight_,"lumiWeight_/D");

	data->Branch("n_pVertex",&n_pVertex,"n_pVertex/I");
	data->Branch("PVx",&PVx,"PVx/D");
	data->Branch("PVy",&PVy,"PVy/D");
	data->Branch("PVz",&PVz,"PVz/D");
	data->Branch("PVerrx",&PVerrx,"PVerrx/D");
	data->Branch("PVerry",&PVerry,"PVerry/D");
	data->Branch("PVerrz",&PVerrz,"PVerrz/D");
	data->Branch("ntracksD0Kpi",&ntracksD0Kpi,"ntracksD0Kpi/I");
	data->Branch("ntracksDstar",&ntracksDstar,"ntracksDstar/I");


	//Counters
	data->Branch("TotalTracks",&TotalTracks,"TotalTracks/L");
 	data->Branch("TracksHasTrackDetails",&TracksHasTrackDetails,"TracksHasTrackDetails/L");
	data->Branch("TracksChargeZero",&TracksChargeZero,"TracksChargeZero/L");
	data->Branch("TracksEta",&TracksEta,"TracksEta/L");
 	data->Branch("TracksHighPurity",&TracksHighPurity,"TracksHighPurity/L");
	data->Branch("TracksPDG211",&TracksPDG211,"TracksPDG211/L");
	data->Branch("TracksPtZeroFive",&TracksPtZeroFive,"TracksPtZeroFive/L");
 	data->Branch("TracksChi3",&TracksChi3,"TracksChi3/L");
	data->Branch("TracksNumberOfHits2",&TracksNumberOfHits2,"TracksNumberOfHits2/L");
	data->Branch("TracksDxyThree",&TracksDxyThree,"TracksDxyThree/L");
 	data->Branch("TracksDzThree",&TracksDzThree,"TracksDzThree/L");
	data->Branch("TrackSlowPionCandidates",&TrackSlowPionCandidates,"TrackSlowPionCandidates/L");
	data->Branch("Observation",&Observation,"Observation/L");
	data->Branch("TracksPtZeroSix",&TracksPtZeroSix,"TracksPtZeroSix/L");
	data->Branch("TracksChiTwoFive",&TracksChiTwoFive,"TracksChiTwoFive/L");
	data->Branch("TracksNumberOfHits5",&TracksNumberOfHits5,"TracksNumberOfHits5/L");
	data->Branch("TracksNumberOfPixelHits2",&TracksNumberOfPixelHits2,"TracksNumberOfPixelHits2/L");
	data->Branch("TracksDxyZeroOne",&TracksDxyZeroOne,"TracksDxyZeroOne/L");
	data->Branch("TracksDzOne",&TracksDzOne,"TracksDzOne/L");
	data->Branch("TrackKaonPionCandidates",&TrackKaonPionCandidates,"TrackKaonPionCandidates/L");
	
//	data->Branch("HLTPath_",&HLTPath_,"HLTPath_/I");

	//======================================================
	// D* -> D0 + pi_slow  Variables
	//======================================================

	//Counters
	data->Branch("D0AfterLorentzVector",&D0AfterLorentzVector,"D0AfterLorentzVector/L");
	data->Branch("D0MinusPDGzero2",&D0MinusPDGzero2,"D0MinusPDGzero2/L");
	data->Branch("DsMinusD0Zerothree",&DsMinusD0Zerothree,"DsMinusD0Zerothree/L");
	data->Branch("TransientTrackOfpiK",&TransientTrackOfpiK,"TransientTrackOfpiK/L");
	data->Branch("SVConfidenceLevel",&SVConfidenceLevel,"SVConfidenceLevel/L");
	data->Branch("D0AfterLorentzVectorKarman",&D0AfterLorentzVectorKarman,"D0AfterLorentzVectorKarman/L");
	data->Branch("PointingcosPhi",&PointingcosPhi,"PointingcosPhi/L");
	data->Branch("Significance",&Significance,"Significance/L");
	data->Branch("D0pTThree",&D0pTThree,"D0pTThree/L");
	data->Branch("D0Candidates",&D0Candidates,"D0Candidates/L");
	data->Branch("DsAfterLorentzVector",&DsAfterLorentzVector,"DsAfterLorentzVector/L");
	data->Branch("D0MinusPDG",&D0MinusPDG,"D0MinusPDG/L");
	data->Branch("DsMinusD0",&DsMinusD0,"DsMinusD0/L");
	data->Branch("DsCandidates",&DsCandidates,"DsCandidates/L");

	data->Branch("D0mass",&D0mass);
	data->Branch("Dsmass",&Dsmass);
	data->Branch("D0_VtxProb",&D0_VtxProb);
	data->Branch("D0_VtxPosx",&D0_VtxPosx);
	data->Branch("D0_VtxPosy",&D0_VtxPosy);
	data->Branch("D0_VtxPosz",&D0_VtxPosz);
	data->Branch("D0_Vtxerrx",&D0_Vtxerrx);
	data->Branch("D0_Vtxerry",&D0_Vtxerry);
	data->Branch("D0_Vtxerrz",&D0_Vtxerrz);
	data->Branch("D0eta",&D0eta);
	data->Branch("D0phi",&D0phi);
	data->Branch("Dseta",&Dseta);
	data->Branch("Dsphi",&Dsphi);
	data->Branch("TrkKmass",&TrkKmass);
	data->Branch("Trkpimass",&Trkpimass);
	data->Branch("TrkSmass",&TrkSmass);
	data->Branch("TrkKpt",&TrkKpt);
	data->Branch("Trkpipt",&Trkpipt);
	data->Branch("TrkSpt",&TrkSpt);
	data->Branch("D0pt",&D0pt);
	data->Branch("Dspt",&Dspt);
	data->Branch("DSDeltaR",&DSDeltaR);
	data->Branch("TrkKnhits",&TrkKnhits);
	data->Branch("Trkpinhits",&Trkpinhits);
	data->Branch("TrkSnhits",&TrkSnhits);
	data->Branch("TrkKchi2",&TrkKchi2);
	data->Branch("Trkpichi2",&Trkpichi2);
	data->Branch("TrkSchi2",&TrkSchi2);
	data->Branch("TrkKdxy",&TrkKdxy);
	data->Branch("Trkpidxy",&Trkpidxy);
	data->Branch("TrkSdxy",&TrkSdxy);
	data->Branch("TrkKdz",&TrkKdz);
	data->Branch("Trkpidz",&Trkpidz);
	data->Branch("TrkSdz",&TrkSdz); 
	data->Branch("TrkKeta",&TrkKeta);      
	data->Branch("Trkpieta",&Trkpieta);
	data->Branch("TrkSeta",&TrkSeta);
	data->Branch("TrkKphi",&TrkKphi);
	data->Branch("Trkpiphi",&Trkpiphi);
	data->Branch("TrkSphi",&TrkSphi);
	data->Branch("TrkScharge",&TrkScharge);
	data->Branch("Anglephi",&Anglephi_vec);
	data->Branch("D0fromDSsXY",&D0fromDSsXY_vec);
	data->Branch("D0fromDSs3D",&D0fromDSs3D_vec);

	//data->Branch("triggers",&triggers);

	//======================================================
	// D* -> D0 + pi_slow  WRONG COMBINATION
	//======================================================

	data->Branch("D0massWrong",&D0massWrong);
	data->Branch("DsmassWrong",&DsmassWrong);
	data->Branch("D0_VtxProbWrong",&D0_VtxProbWrong);
	data->Branch("D0_VtxPosxWrong",&D0_VtxPosxWrong);
	data->Branch("D0_VtxPosyWrong",&D0_VtxPosyWrong);
	data->Branch("D0_VtxPoszWrong",&D0_VtxPoszWrong);
	data->Branch("D0_VtxerrxWrong",&D0_VtxerrxWrong);
	data->Branch("D0_VtxerryWrong",&D0_VtxerryWrong);
	data->Branch("D0_VtxerrzWrong",&D0_VtxerrzWrong);
	data->Branch("D0etaWrong",&D0etaWrong);
	data->Branch("D0phiWrong",&D0phiWrong);
	data->Branch("DsetaWrong",&DsetaWrong);
	data->Branch("DsphiWrong",&DsphiWrong);
	data->Branch("TrkKmassWrong",&TrkKmassWrong);
	data->Branch("TrkpimassWrong",&TrkpimassWrong);
	data->Branch("TrkSmassWrong",&TrkSmassWrong);
	data->Branch("TrkKptWrong",&TrkKptWrong);
	data->Branch("TrkpiptWrong",&TrkpiptWrong);
	data->Branch("TrkSptWrong",&TrkSptWrong);
	data->Branch("D0ptWrong",&D0ptWrong);
	data->Branch("DsptWrong",&DsptWrong);
	data->Branch("DSDeltaRWrong",&DSDeltaRWrong);
	data->Branch("TrkKnhitsWrong",&TrkKnhitsWrong);
	data->Branch("TrkpinhitsWrong",&TrkpinhitsWrong);
	data->Branch("TrkSnhitsWrong",&TrkSnhitsWrong);
	data->Branch("TrkKchi2Wrong",&TrkKchi2Wrong);
	data->Branch("Trkpichi2Wrong",&Trkpichi2Wrong);
	data->Branch("TrkSchi2Wrong",&TrkSchi2Wrong);
	data->Branch("TrkKdxyWrong",&TrkKdxyWrong);
	data->Branch("TrkpidxyWrong",&TrkpidxyWrong);
	data->Branch("TrkSdxyWrong",&TrkSdxyWrong);
	data->Branch("TrkKdzWrong",&TrkKdzWrong);
	data->Branch("TrkpidzWrong",&TrkpidzWrong);
	data->Branch("TrkSdzWrong",&TrkSdzWrong); 
	data->Branch("TrkKetaWrong",&TrkKetaWrong);      
	data->Branch("TrkpietaWrong",&TrkpietaWrong);
	data->Branch("TrkSetaWrong",&TrkSetaWrong);
	data->Branch("TrkKphiWrong",&TrkKphiWrong);
	data->Branch("TrkpiphiWrong",&TrkpiphiWrong);
	data->Branch("TrkSphiWrong",&TrkSphiWrong);
	data->Branch("TrkSchargeWrong",&TrkSchargeWrong);
	data->Branch("AnglephiWrong",&Anglephi_vecWrong);
	data->Branch("D0fromDSsXYWrong",&D0fromDSsXY_vecWrong);
	data->Branch("D0fromDSs3DWrong",&D0fromDSs3D_vecWrong);

	//data->Branch("triggers",&triggers);


    //======================================================
	// D0 Prompt Variables
	//======================================================
	//Counters
	data->Branch("TracksPt0p8",&TracksPt0p8,"TracksPt0p8/L");
	data->Branch("D0PromptTracksEta2p5",&D0PromptTracksEta2p5,"D0PromptTracksEta2p5/L");
	data->Branch("D0PromptTracksChi5",&D0PromptTracksChi5,"D0PromptTracksChi5/L");
	data->Branch("D0PromptTracksNumberHits5",&D0PromptTracksNumberHits5,"D0PromptTracksNumberHits5/L");
	data->Branch("D0PromptTracksPixelHits2",&D0PromptTracksPixelHits2,"D0PromptTracksPixelHits2/L");
	data->Branch("D0PromptTracksDxy0p1",&D0PromptTracksDxy0p1,"D0PromptTracksDxy0p1/L");
	data->Branch("D0PromptTracksDz0p5",&D0PromptTracksDz0p5,"D0PromptTracksDz0p5/L");
	data->Branch("TracksD0combination",&TracksD0combination,"TracksD0combination/L");
	data->Branch("D0PromptminusPDG1p0",&D0PromptminusPDG1p0,"D0PromptminusPDG1p0/L");
	data->Branch("D0PromptKpiAfterTransientp0",&D0PromptKpiAfterTransientp0,"D0PromptKpiAfterTransientp0/L");
	data->Branch("D0PromptSVConfidenceLevel",&D0PromptSVConfidenceLevel,"D0PromptSVConfidenceLevel/L");
	data->Branch("D0PromptPointingCosPhi",&D0PromptPointingCosPhi,"D0PromptPointingCosPhi/L");
	data->Branch("D0PromptSignificance",&D0PromptSignificance,"D0PromptSignificance/L");
	data->Branch("D0PromptCandidates",&D0PromptCandidates,"D0PromptCandidates/L");
	data->Branch("ND0KpiCand",&ND0KpiCand,"ND0KpiCand/L");

	//Variables data
	data->Branch("D0Kpi_VtxProb",&D0Kpi_VtxProb);
    data->Branch("D0Kpimass",&D0Kpimass);
    data->Branch("D0Kpipt",&D0Kpipt);
    data->Branch("D0Kpieta",&D0Kpieta);
    data->Branch("D0Kpiphi",&D0Kpiphi);
    data->Branch("D0Kpi_VtxPosx",&D0Kpi_VtxPosx);
    data->Branch("D0Kpi_VtxPosy",&D0Kpi_VtxPosy);
    data->Branch("D0Kpi_VtxPosz",&D0Kpi_VtxPosz);
    data->Branch("D0Kpi_Vtxerrx",&D0Kpi_Vtxerrx);
    data->Branch("D0Kpi_Vtxerry",&D0Kpi_Vtxerry);
    data->Branch("D0Kpi_Vtxerrz",&D0Kpi_Vtxerrz);
    data->Branch("TrkD0Kdxy",&TrkD0Kdxy);
    data->Branch("TrkD0pidxy",&TrkD0pidxy);
    data->Branch("TrkD0Kdz",&TrkD0Kdz);
    data->Branch("TrkD0pidz",&TrkD0pidz);
	data->Branch("TrkD0Kchi2",&TrkD0Kchi2);
	data->Branch("TrkD0pichi2",&TrkD0pichi2);
	data->Branch("TrkD0Kpt",&TrkD0Kpt);
	data->Branch("TrkD0pipt",&TrkD0pipt);
	data->Branch("TrkD0Keta",&TrkD0Keta);
	data->Branch("TrkD0pieta",&TrkD0pieta);
	data->Branch("TrkD0kphi",&TrkD0Kphi);
	data->Branch("TrkD0piphi",&TrkD0piphi);
	data->Branch("TrkD0Knhits",&TrkD0Knhits);
	data->Branch("TrkD0pinhits",&TrkD0pinhits);
	data->Branch("D0KpisXY",&D0KpisXY_vec);
	data->Branch("D0Kpis3D",&D0Kpis3D_vec);
	data->Branch("D0KpikT",&D0_kT_vec);
	data->Branch("D0KpiDispAngle",&D0Kpi_DispAngle);
	
//======================================================
// MC Variables - D_star
//======================================================
	        
	data->Branch("dScandsKpi",&dScandsKpi);
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

//======================================================
// MC Variables - D_0 Prompt
//======================================================

	data->Branch("MCpromptD0eta",&MCpromptD0eta);
	data->Branch("MCpromptD0phi",&MCpromptD0phi);
	data->Branch("MCpromptD0pt",&MCpromptD0pt);
	data->Branch("MCpromptD0energy",&MCpromptD0energy);
	data->Branch("MCpromptD0p",&MCpromptD0p);
	data->Branch("MCpromptD0et",&MCpromptD0et);
	data->Branch("MCpromptD0rapidity",&MCpromptD0rapidity);
	data->Branch("MCpromptD0mass",&MCpromptD0mass);
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

	//data->Branch("MCpromptD0_DispAngle",&MCpromptD0_DispAngle);
	//data->Branch("MCpromptD0_Kt",&MCpromptD0_Kt); 

}

//define this as a plug-in
DEFINE_FWK_MODULE(DstarD0TTree);
